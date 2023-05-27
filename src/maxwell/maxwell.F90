!! Copyright (C) 2019-2020 Franco Bonafe, Heiko Appel, Rene Jestaedt
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module maxwell_oct_m
  use algorithm_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use debug_oct_m
  use distributed_oct_m
  use external_densities_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_mxll_oct_m
  use index_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use iso_c_binding
  use lattice_vectors_oct_m
  use linear_medium_to_em_field_oct_m
  use current_to_mxll_field_oct_m
  use loct_oct_m
  use lorentz_force_oct_m
  use maxwell_boundary_op_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mesh_interpolation_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use mxll_field_to_medium_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use propagator_mxll_oct_m
  use quantity_oct_m
  use regridding_oct_m
  use restart_oct_m
  use sort_oct_m
  use space_oct_m
  use system_oct_m
  use states_mxll_oct_m
  use states_mxll_restart_oct_m
  use electrons_oct_m
  use td_write_oct_m
  use unit_oct_m
  use unit_system_oct_m


  implicit none

  private
  public ::               &
    maxwell_t,        &
    maxwell_init

  integer, parameter, public ::           &
    MULTIGRID_MX_TO_MA_EQUAL   = 1,       &
    MULTIGRID_MX_TO_MA_LARGE   = 2

  type, extends(system_t) :: maxwell_t
    type(states_mxll_t)          :: st    !< the states
    type(hamiltonian_mxll_t)     :: hm
    type(grid_t)                 :: gr    !< the mesh
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators

    type(mesh_interpolation_t)   :: mesh_interpolate

    type(propagator_mxll_t)      :: tr_mxll   !< contains the details of the Maxwell time-evolution
    type(td_write_t)             :: write_handler
    type(c_ptr)                  :: output_handle

    CMPLX, allocatable           :: ff_rs_inhom_t1(:,:), ff_rs_inhom_t2(:,:)
    CMPLX, allocatable           :: rs_state_init(:,:)
    CMPLX, allocatable           :: current_density_medium(:,:,:) ! (np, 3, 2), last index is for two times (t1 and t2)
    FLOAT                        :: bc_bounds(2,MAX_DIM), dt_bounds(2,MAX_DIM)
    integer                      :: energy_update_iter
    type(restart_t)              :: restart_dump
    type(restart_t)              :: restart

  contains
    procedure :: init_interaction => maxwell_init_interaction
    procedure :: init_parallelization => maxwell_init_parallelization
    procedure :: initial_conditions => maxwell_initial_conditions
    procedure :: do_td_operation => maxwell_do_td
    procedure :: is_tolerance_reached => maxwell_is_tolerance_reached
    procedure :: update_quantity => maxwell_update_quantity
    procedure :: update_exposed_quantity => maxwell_update_exposed_quantity
    procedure :: init_interaction_as_partner => maxwell_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => maxwell_copy_quantities_to_interaction
    procedure :: output_start => maxwell_output_start
    procedure :: output_write => maxwell_output_write
    procedure :: output_finish => maxwell_output_finish
    procedure :: update_interactions_start => maxwell_update_interactions_start
    procedure :: update_interactions_finish => maxwell_update_interactions_finish
    procedure :: restart_write_data => maxwell_restart_write_data
    procedure :: restart_read_data => maxwell_restart_read_data
    procedure :: update_kinetic_energy => maxwell_update_kinetic_energy
    final :: maxwell_finalize
  end type maxwell_t

  interface maxwell_t
    procedure maxwell_constructor
  end interface maxwell_t

contains

  ! ---------------------------------------------------------
  function maxwell_constructor(namespace) result(sys)
    class(maxwell_t),   pointer    :: sys
    type(namespace_t),  intent(in) :: namespace

    PUSH_SUB(maxwell_constructor)

    SAFE_ALLOCATE(sys)

    call maxwell_init(sys, namespace)

    POP_SUB(maxwell_constructor)
  end function maxwell_constructor


  ! ---------------------------------------------------------
  subroutine maxwell_init(this, namespace)
    class(maxwell_t),     intent(inout) :: this
    type(namespace_t),    intent(in)    :: namespace

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_init)

    call profiling_in(prof,"MAXWELL_INIT")

    this%namespace = namespace

    call messages_obsolete_variable(this%namespace, 'SystemName')
    call space_init(this%space, this%namespace)
    if (this%space%is_periodic()) then
      call messages_not_implemented('Maxwell for periodic systems', namespace=namespace)
    end if

    call grid_init_stage_1(this%gr, this%namespace, this%space)
    call states_mxll_init(this%st, this%namespace, this%space)

    this%quantities(E_FIELD)%required = .true.
    this%quantities(E_FIELD)%protected = .true.
    this%quantities(B_FIELD)%required = .true.
    this%quantities(B_FIELD)%protected = .true.

    call mesh_interpolation_init(this%mesh_interpolate, this%gr%mesh)

    call this%supported_interactions_as_partner%add(LORENTZ_FORCE)
    call this%supported_interactions_as_partner%add(MXLL_FIELD_TO_MEDIUM)
    call this%supported_interactions%add(LINEAR_MEDIUM_TO_EM_FIELD)
    call this%supported_interactions%add(CURRENT_TO_MXLL_FIELD)

    call profiling_out(prof)

    POP_SUB(maxwell_init)
  end subroutine maxwell_init

  ! ---------------------------------------------------------
  subroutine maxwell_init_interaction(this, interaction)
    class(maxwell_t),     target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(maxwell_init_interaction)

    select type (interaction)
    type is (linear_medium_to_em_field_t)
      call interaction%init(this%gr)
    type is (current_to_mxll_field_t)
      call interaction%init(this%gr)
      this%hm%current_density_from_medium = .true.
      SAFE_ALLOCATE(this%current_density_medium(1:this%gr%mesh%np, 1:this%gr%box%dim, 1:2))
    class default
      message(1) = "Trying to initialize an unsupported interaction by Maxwell."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(maxwell_init_interaction)
  end subroutine maxwell_init_interaction

  ! ---------------------------------------------------------
  subroutine maxwell_init_parallelization(this, grp)
    class(maxwell_t),     intent(inout) :: this
    type(mpi_grp_t),      intent(in)    :: grp

    integer(i8) :: index_range(4)
    integer :: ierr

    PUSH_SUB(maxwell_init_parallelization)

    call system_init_parallelization(this, grp)

    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%mesh%np_global  ! Number of points in mesh
    index_range(2) = this%st%nst             ! Number of states
    index_range(3) = 1                      ! Number of k-points
    index_range(4) = 100000                 ! Some large number

    ! create index and domain communicators
    call multicomm_init(this%mc, this%namespace, mpi_world, calc_mode_par_parallel_mask(), &
    &calc_mode_par_default_parallel_mask(),mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

    call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc)
    call output_mxll_init(this%outp, this%namespace, this%space)
    call hamiltonian_mxll_init(this%hm, this%namespace, this%gr, this%st)

    this%st%energy_rate = M_ZERO
    this%st%delta_energy = M_ZERO
    this%st%energy_via_flux_calc = M_ZERO
    this%st%trans_energy_rate = M_ZERO
    this%st%trans_delta_energy = M_ZERO
    this%st%trans_energy_via_flux_calc = M_ZERO
    this%st%plane_waves_energy_rate = M_ZERO
    this%st%plane_waves_delta_energy = M_ZERO
    this%st%plane_waves_energy_via_flux_calc = M_ZERO

    SAFE_ALLOCATE(this%rs_state_init(1:this%gr%mesh%np_part, 1:this%st%dim))
    this%rs_state_init(:,:) = M_z0

    this%energy_update_iter = 1

    call propagator_mxll_init(this%gr, this%namespace, this%st, this%hm, this%tr_mxll)
    call states_mxll_allocate(this%st, this%gr%mesh)
    call external_current_init(this%st, this%namespace, this%gr%mesh)
    this%hm%propagation_apply = .true.

    if (parse_is_defined(this%namespace, 'MaxwellIncidentWaves') .and. (this%tr_mxll%bc_plane_waves)) then
      this%st%rs_state_plane_waves(:,:) = M_z0
    end if

    this%hm%plane_waves_apply = .true.
    this%hm%spatial_constant_apply = .true.

    ! TODO: Calculate rs_trans using Helmholtz decomposition
    this%st%rs_state_trans(:,:) = this%st%rs_state

    call bc_mxll_init(this%hm%bc, this%namespace, this%space, this%gr, this%st)
    this%bc_bounds(:,1:3) = this%hm%bc%bc_bounds(:,1:3)
    call inner_and_outer_points_mapping(this%gr%mesh, this%st, this%bc_bounds)
    this%dt_bounds(2, 1:3) = this%bc_bounds(1, 1:3)
    this%dt_bounds(1, 1:3) = this%bc_bounds(1, 1:3) - this%gr%der%order * this%gr%mesh%spacing(1:3)
    call surface_grid_points_mapping(this%gr%mesh, this%st, this%dt_bounds)

    call restart_init(this%restart, this%namespace, RESTART_TD, RESTART_TYPE_LOAD, this%mc, ierr, mesh=this%gr%mesh)
    call restart_init(this%restart_dump, this%namespace, RESTART_TD, RESTART_TYPE_DUMP, this%mc, ierr, mesh=this%gr%mesh)

    POP_SUB(maxwell_init_parallelization)
  end subroutine maxwell_init_parallelization

  ! ---------------------------------------------------------
  subroutine maxwell_initial_conditions(this)
    class(maxwell_t), intent(inout) :: this

    FLOAT   :: courant

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_initial_conditions)

    call profiling_in(prof,"MAXWELL_INITIAL_CONDITIONS")

    courant = M_ONE/(P_c * sqrt(M_ONE/this%gr%mesh%spacing(1)**2 + M_ONE/this%gr%mesh%spacing(2)**2 + &
      M_ONE/this%gr%mesh%spacing(3)**2))

    if (this%prop%dt > M_TWO * courant) then
      write(message(1),'(a,es9.2)') 'Time step seems too large, check this value'
      call messages_warning(1, namespace=this%namespace)
    end if

    if (parse_is_defined(this%namespace, 'UserDefinedInitialMaxwellStates')) then
      call states_mxll_read_user_def(this%namespace, this%space, this%gr%mesh, this%st, this%hm%bc, this%rs_state_init)
      call messages_print_stress(msg="Setting initial EM field inside box", namespace=this%namespace)
      ! TODO: add consistency check that initial state fulfills Gauss laws
      this%st%rs_state(:,:) = this%st%rs_state + this%rs_state_init
      if (this%tr_mxll%bc_plane_waves) then
        this%st%rs_state_plane_waves(:,:) = this%rs_state_init
      end if
    end if

    ! initialize the spatial constant field according to the conditions set in the
    ! UserDefinedConstantSpatialMaxwellField block
    if (this%tr_mxll%bc_constant) then
      call spatial_constant_calculation(this%tr_mxll%bc_constant, this%st, this%gr, this%hm, M_ZERO, &
        this%prop%dt, M_ZERO, this%st%rs_state, set_initial_state = .true.)
      ! for mesh parallelization, this needs communication!
      this%st%rs_state_const(:) = this%st%rs_state(mesh_global_index_from_coords(this%gr%mesh, [0,0,0]),:)
    end if

    if (parse_is_defined(this%namespace, 'UserDefinedInitialMaxwellStates')) then
      SAFE_DEALLOCATE_A(this%rs_state_init)
    end if

    call hamiltonian_mxll_update(this%hm, time = M_ZERO)

    ! calculate Maxwell energy
    call energy_mxll_calc(this%gr, this%st, this%hm, this%hm%energy, this%st%rs_state, &
      this%st%rs_state_plane_waves)

    ! TODO: Calculate rs_trans using Helmholtz decomposition
    this%st%rs_state_trans(:,:) = this%st%rs_state

    call get_rs_state_at_point(this%st%selected_points_rs_state(:,:), this%st%rs_state, &
      this%st%selected_points_coordinate(:,:), this%st, this%gr%mesh)

    call profiling_out(prof)

    POP_SUB(maxwell_initial_conditions)
  end subroutine maxwell_initial_conditions

  ! ---------------------------------------------------------
  subroutine maxwell_do_td(this, operation)
    class(maxwell_t),               intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    CMPLX, allocatable :: charge_density_ext(:)
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_do_td)

    select case (operation%id)
    case (SKIP)
      ! Do nothing

    case (STORE_CURRENT_STATUS)
      ! For the moment we do nothing

    case (EXPMID_START)
      SAFE_ALLOCATE(this%ff_rs_inhom_t1(1:this%gr%mesh%np_part, 1:this%hm%dim))
      SAFE_ALLOCATE(this%ff_rs_inhom_t2(1:this%gr%mesh%np_part, 1:this%hm%dim))

    case (EXPMID_FINISH)

      SAFE_DEALLOCATE_A(this%ff_rs_inhom_t1)
      SAFE_DEALLOCATE_A(this%ff_rs_inhom_t2)

    case (EXPMID_PREDICT_DT_2)  ! predict: psi(t+dt/2) = 0.5*(U_H(dt) psi(t) + psi(t)) or via extrapolation

      ! calculation of external RS density at time (time-dt)
      if (this%hm%current_density_ext_flag) then
        call get_rs_density_ext(this%st, this%gr%mesh, this%clock%time(), this%st%rs_current_density_t1)
        call get_rs_density_ext(this%st, this%gr%mesh, this%clock%time()+this%prop%dt, this%st%rs_current_density_t2)
      else
        this%st%rs_current_density_t1(:,:) = M_z0
        this%st%rs_current_density_t2(:,:) = M_z0
      end if
      if (this%hm%current_density_from_medium) then
        this%st%rs_current_density_t1(:,:) = this%st%rs_current_density_t1 + this%current_density_medium(:,:,1)
        this%st%rs_current_density_t2(:,:) = this%st%rs_current_density_t2 + this%current_density_medium(:,:,2)
      end if

      this%quantities(E_FIELD)%clock = this%quantities(E_FIELD)%clock + CLOCK_TICK
      this%quantities(B_FIELD)%clock = this%quantities(B_FIELD)%clock + CLOCK_TICK

    case (UPDATE_HAMILTONIAN)   ! update: H(t+dt/2) from psi(t+dt/2)
      ! Empty for the moment
    case (EXPMID_PREDICT_DT)    ! predict: psi(t+dt) = U_H(t+dt/2) psi(t)

      call profiling_in(prof, "SYSTEM_MXLL_DO_TD")

      ! Propagation

      !We first compute thre external charge and current densities and we convert them as RS vectors
      SAFE_ALLOCATE(charge_density_ext(1:this%gr%mesh%np))

      !No charge density at the moment
      charge_density_ext = M_z0

      call transform_rs_densities(this%hm, this%gr%mesh, charge_density_ext, &
        this%st%rs_current_density_t1, this%ff_rs_inhom_t1, RS_TRANS_FORWARD)
      call transform_rs_densities(this%hm, this%gr%mesh, charge_density_ext, &
        this%st%rs_current_density_t2, this%ff_rs_inhom_t2, RS_TRANS_FORWARD)

      SAFE_DEALLOCATE_A(charge_density_ext)

      ! Propagation dt with H_maxwell
      call mxll_propagation_step(this%hm, this%namespace, this%gr, this%space, this%st, this%tr_mxll,&
        this%st%rs_state, this%ff_rs_inhom_t1, this%ff_rs_inhom_t2, this%clock%time(), this%prop%dt)

      ! TODO: Calculate rs_trans using Helmholtz decomposition
      this%st%rs_state_trans(:,:) = this%st%rs_state

      ! calculate Maxwell energy
      call energy_mxll_calc(this%gr, this%st, this%hm, this%hm%energy, this%st%rs_state, &
        this%st%rs_state_plane_waves)

      ! get RS state values for selected points
      call get_rs_state_at_point(this%st%selected_points_rs_state(:,:), this%st%rs_state, this%st%selected_points_coordinate(:,:), &
        this%st, this%gr%mesh)

      this%quantities(E_FIELD)%clock = this%quantities(E_FIELD)%clock + CLOCK_TICK
      this%quantities(B_FIELD)%clock = this%quantities(B_FIELD)%clock + CLOCK_TICK

      call profiling_out(prof)

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(maxwell_do_td)
  end subroutine maxwell_do_td

  ! ---------------------------------------------------------
  logical function maxwell_is_tolerance_reached(this, tol) result(converged)
    class(maxwell_t),       intent(in)    :: this
    FLOAT,                  intent(in)    :: tol

    PUSH_SUB(maxwell_is_tolerance_reached)

    converged = .false.

    POP_SUB(maxwell_is_tolerance_reached)
  end function maxwell_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine maxwell_update_quantity(this, iq)
    class(maxwell_t),          intent(inout) :: this
    integer,                   intent(in)    :: iq

    PUSH_SUB(maxwell_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(maxwell_update_quantity)
  end subroutine maxwell_update_quantity

  ! ---------------------------------------------------------
  subroutine maxwell_update_exposed_quantity(partner, iq)
    class(maxwell_t),      intent(inout) :: partner
    integer,                   intent(in)    :: iq

    PUSH_SUB(maxwell_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(maxwell_update_exposed_quantity)
  end subroutine maxwell_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine maxwell_init_interaction_as_partner(partner, interaction)
    class(maxwell_t),           intent(in)    :: partner
    class(interaction_t),       intent(inout) :: interaction

    PUSH_SUB(maxwell_init_interaction_as_partner)

    select type (interaction)
    type is (lorentz_force_t)
      ! Nothing to be initialized for the Lorentz force.
    type is (mxll_field_to_medium_t)
      SAFE_ALLOCATE(interaction%partner_E_field(1:partner%gr%mesh%np, 1:3))
      interaction%partner_np = partner%gr%mesh%np
      call grid_transfer_mapping(interaction%system_gr, partner%gr, interaction%partner_to_system_map, &
           partner%namespace)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(maxwell_init_interaction_as_partner)
  end subroutine maxwell_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine maxwell_copy_quantities_to_interaction(partner, interaction)
    class(maxwell_t),           intent(inout) :: partner
    class(interaction_t),       intent(inout) :: interaction

    integer :: ip
    CMPLX :: interpolated_value(3)
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_copy_quantities_to_interaction)

    call profiling_in(prof, "MAXWELL_COPY_QUANTITIES_TO_INTERACTION")

    select type (interaction)
    type is (lorentz_force_t)
      do ip = 1, interaction%system_np
        call mesh_interpolation_evaluate(partner%mesh_interpolate, partner%st%rs_state(:,1), &
          interaction%system_pos(:, ip), interpolated_value(1))
        call mesh_interpolation_evaluate(partner%mesh_interpolate, partner%st%rs_state(:,2), &
          interaction%system_pos(:, ip), interpolated_value(2))
        call mesh_interpolation_evaluate(partner%mesh_interpolate, partner%st%rs_state(:,3), &
          interaction%system_pos(:, ip), interpolated_value(3))
        call get_electric_field_vector(interpolated_value, interaction%partner_E_field(:, ip))
        call get_magnetic_field_vector(interpolated_value, 1, interaction%partner_B_field(:, ip))
      end do

    type is (mxll_field_to_medium_t)
       call get_electric_field_state(partner%st%rs_state, partner%gr%mesh, interaction%partner_E_field)

    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    call profiling_out(prof)

    POP_SUB(maxwell_copy_quantities_to_interaction)
  end subroutine maxwell_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine maxwell_output_start(this)
    class(maxwell_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_output_start)

    call profiling_in(prof, "MAXWELL_OUTPUT_START")

    call td_write_mxll_init(this%write_handler, this%namespace, this%clock%get_tick(), this%prop%dt)
    if (this%st%fromScratch) then
      call td_write_mxll_iter(this%write_handler, this%space, this%gr, this%st, this%hm, this%prop%dt, this%clock%get_tick())
      call td_write_mxll_free_data(this%write_handler, this%namespace, this%space, this%gr, this%st, this%hm, &
        this%outp, this%clock)
    end if

    call profiling_out(prof)

    POP_SUB(maxwell_output_start)
  end subroutine maxwell_output_start

  ! ---------------------------------------------------------
  subroutine maxwell_output_write(this)
    class(maxwell_t), intent(inout) :: this

    logical :: stopping, reached_output_interval
    type(profile_t), save :: prof

    integer :: iout

    PUSH_SUB(maxwell_output_write)

    call profiling_in(prof, "MAXWELL_OUTPUT_WRITE")

    stopping = clean_stop(this%mc%master_comm)

    call td_write_mxll_iter(this%write_handler, this%space, this%gr, this%st, this%hm, this%prop%dt, this%clock%get_tick())

    reached_output_interval = .false.
    do iout = 1, OUT_MAXWELL_MAX
      if (this%outp%output_interval(iout) > 0) then
        if (mod(this%clock%get_tick(), this%outp%output_interval(iout)) == 0) then
          reached_output_interval = .true.
          exit
        end if
      end if
    end do

    if (reached_output_interval .or. stopping) then
      call td_write_mxll_free_data(this%write_handler, this%namespace, this%space, this%gr, this%st, this%hm, &
        this%outp, this%clock)
    end if

    call profiling_out(prof)

    POP_SUB(maxwell_output_write)
  end subroutine maxwell_output_write

  ! ---------------------------------------------------------
  subroutine maxwell_output_finish(this)
    class(maxwell_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_output_finish)

    call profiling_in(prof, "MAXWELL_OUTPUT_FINISH")

    call td_write_mxll_end(this%write_handler)

    call profiling_out(prof)

    POP_SUB(maxwell_output_finish)
  end subroutine maxwell_output_finish

  ! ---------------------------------------------------------
  subroutine maxwell_update_interactions_start(this)
    class(maxwell_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    integer :: int_counter

    PUSH_SUB(maxwell_update_interactions_start)

    int_counter = 0
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (linear_medium_to_em_field_t)
        int_counter = int_counter + 1
      class is (current_to_mxll_field_t)
        this%current_density_medium(:,:,1) = this%current_density_medium(:,:,2)
      end select
    end do

    if (int_counter /= 0 .and. .not. allocated(this%hm%medium_boxes)) then
      SAFE_ALLOCATE(this%hm%medium_boxes(1:int_counter))
      this%hm%calc_medium_box = .true.
    end if

    POP_SUB(maxwell_update_interactions_start)
  end subroutine maxwell_update_interactions_start

  ! ---------------------------------------------------------
  subroutine maxwell_update_interactions_finish(this)
    class(maxwell_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    integer :: iint

    PUSH_SUB(maxwell_update_interactions_finish)

    iint = 0
    if (this%hm%current_density_from_medium) this%current_density_medium(:,:,2) = M_z0

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (linear_medium_to_em_field_t)
        if (allocated(this%hm%medium_boxes) .and. .not. this%hm%medium_boxes_initialized) then
          iint = iint + 1
          this%hm%medium_boxes(iint) = interaction%medium_box
        end if
      class is (current_to_mxll_field_t)
        this%current_density_medium(:,:,2) = this%current_density_medium(:,:,2) + interaction%rs_current_p(:,:)
      end select
    end do

    if (allocated(this%hm%medium_boxes) .and. .not. this%hm%medium_boxes_initialized) then
      call set_medium_rs_state(this%st, this%gr, this%hm)
      this%hm%medium_boxes_initialized = .true.
    end if

    if (this%hm%medium_boxes_initialized .and. this%hm%operator /= FARADAY_AMPERE_MEDIUM) then
      message(1) = "A linear medium has been defined in the input file but the Hamiltonian"
      message(2) = "type you specified is not capable of dealing with the medium."
      message(3) = "Please use MaxwellHamiltonianOperator = faraday_ampere_medium to enable"
      message(4) = "the medium propagation."
      call messages_fatal(4, namespace=this%namespace)
    end if

    if (.not. this%hm%medium_boxes_initialized .and. this%hm%operator == FARADAY_AMPERE_MEDIUM) then
      message(1) = "The variable MaxwellHamiltonianOperator has been defined as faraday_ampere_medium"
      message(2) = "in the input file but no linear medium has been defined in the system block."
      message(3) = "Please either use a different option for MaxwellHamiltonianOperator or add"
      message(4) = "a linear medium to the system block."
      call messages_fatal(4, namespace=this%namespace)
    end if

    POP_SUB(maxwell_update_interactions_finish)
  end subroutine maxwell_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine maxwell_restart_write_data(this)
    class(maxwell_t), intent(inout) :: this

    integer :: ierr, err, zff_dim, id, id1, id2, ip_in
    logical :: pml_check
    CMPLX, allocatable :: zff(:,:)


    PUSH_SUB(maxwell_restart_write_data)
    ierr = 0

    pml_check = any(this%hm%bc%bc_ab_type(1:3) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML)

    if (debug%info) then
      message(1) = "Debug: Writing td_maxwell restart."
      call messages_info(1, namespace=this%namespace)
    end if

    if (this%tr_mxll%bc_plane_waves) then
      zff_dim = 2 * this%st%dim
    else
      zff_dim = 1 * this%st%dim
    end if
    if (pml_check) then
      zff_dim = zff_dim + 18
    end if

    SAFE_ALLOCATE(zff(1:this%gr%mesh%np,1:zff_dim))
    zff = M_z0

    if (this%tr_mxll%bc_plane_waves) then
      zff(1:this%gr%mesh%np, 1:this%st%dim) = this%st%rs_state(1:this%gr%mesh%np, 1:this%st%dim)
      zff(1:this%gr%mesh%np, this%st%dim+1:this%st%dim+this%st%dim) = &
        this%st%rs_state_plane_waves(1:this%gr%mesh%np, 1:this%st%dim)
      if (pml_check) then
        id = 0
        do id1 = 1, 3
          do id2 = 1, 3
            id = id + 1
            do ip_in = 1, this%hm%bc%pml%points_number
              zff(ip_in, 2*this%st%dim+id) = this%hm%bc%pml%conv_plus(ip_in, id1, id2)
              zff(ip_in, 2*this%st%dim+9+id) = this%hm%bc%pml%conv_minus(ip_in, id1, id2)
            end do
          end do
        end do
      end if
    else
      zff(1:this%gr%mesh%np, 1:this%st%dim) = this%st%rs_state(1:this%gr%mesh%np, 1:this%st%dim)
      if (pml_check) then
        id = 0
        do id1 = 1, 3
          do id2 = 1, 3
            id = id + 1
            do ip_in = 1, this%hm%bc%pml%points_number
              zff(ip_in, this%st%dim+id) = this%hm%bc%pml%conv_plus(ip_in, id1, id2)
              zff(ip_in, this%st%dim+9+id) = this%hm%bc%pml%conv_minus(ip_in, id1, id2)
            end do
          end do
        end do
      end if
    end if

    call states_mxll_dump(this%restart_dump, this%st, this%space, this%gr%mesh, zff, zff_dim, err, this%clock%get_tick())
    if (err /= 0) ierr = ierr + 1

    if (debug%info) then
      message(1) = "Debug: Writing td_maxwell restart done."
      call messages_info(1, namespace=this%namespace)
    end if

    SAFE_DEALLOCATE_A(zff)

    if (ierr /=0) then
      message(1) = "Unable to write time-dependent Maxwell restart information."
      call messages_warning(1, namespace=this%namespace)
    end if

    POP_SUB(maxwell_restart_write_data)
  end subroutine maxwell_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function maxwell_restart_read_data(this)
    class(maxwell_t), intent(inout) :: this

    integer :: ierr, err, zff_dim, id, id1, id2, ip_in
    logical :: pml_check
    CMPLX, allocatable :: zff(:,:)

    PUSH_SUB(maxwell_restart_read_data)

    if (.not. restart_skip(this%restart)) then
      ierr = 0
      pml_check = any(this%hm%bc%bc_ab_type(1:3) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML)

      if (restart_skip(this%restart)) then
        ierr = -1
        POP_SUB(td_load_mxll)
        return
      end if

      if (debug%info) then
        message(1) = "Debug: Reading td_maxwell restart."
        call messages_info(1, namespace=this%namespace)
      end if

      if (this%tr_mxll%bc_plane_waves) then
        zff_dim = 2 * this%st%dim
      else
        zff_dim = 1 * this%st%dim
      end if
      if (pml_check) then
        zff_dim = zff_dim + 18
      end if

      SAFE_ALLOCATE(zff(1:this%gr%mesh%np,1:zff_dim))

      call states_mxll_load(this%restart, this%st, this%gr%mesh, this%namespace, this%space, zff, &
        zff_dim, err, label = ": td_maxwell")
      this%st%rs_current_density_restart = .true.

      if (this%tr_mxll%bc_plane_waves) then
        this%st%rs_state(1:this%gr%mesh%np,1:this%st%dim) = zff(1:this%gr%mesh%np, 1:this%st%dim)
        this%st%rs_state_plane_waves(1:this%gr%mesh%np,1:this%st%dim) = &
          zff(1:this%gr%mesh%np,this%st%dim+1:this%st%dim+3)
        if (pml_check) then
          id = 0
          do id1 = 1, 3
            do id2 = 1, 3
              id = id+1
              do ip_in=1, this%hm%bc%pml%points_number
                this%hm%bc%pml%conv_plus(ip_in,id1,id2)  = zff(ip_in, 2*this%st%dim+  id)
                this%hm%bc%pml%conv_minus(ip_in,id1,id2) = zff(ip_in, 2*this%st%dim+9+id)
              end do
            end do
          end do
        end if
      else
        this%st%rs_state(1:this%gr%mesh%np,1:this%st%dim) = zff(1:this%gr%mesh%np, 1:this%st%dim)
        if (pml_check) then
          id = 0
          do id1 = 1, 3
            do id2 = 1, 3
              id = id+1
              do ip_in = 1, this%hm%bc%pml%points_number
                this%hm%bc%pml%conv_plus(ip_in,id1,id2)  = zff(ip_in, this%st%dim+  id)
                this%hm%bc%pml%conv_minus(ip_in,id1,id2) = zff(ip_in, this%st%dim+9+id)
              end do
            end do
          end do
        end if
      end if

      if (err /= 0) then
        ierr = ierr + 1
      end if

      if (debug%info) then
        message(1) = "Debug: Reading td restart done."
        call messages_info(1, namespace=this%namespace)
      end if
      SAFE_DEALLOCATE_A(zff)

      this%st%fromScratch = .false.
      maxwell_restart_read_data = .true.
    else
      message(1) = "Unable to read time-dependent Maxwell restart information: Starting from scratch"
      call messages_warning(1, namespace=this%namespace)
      maxwell_restart_read_data = .false.
    end if

    POP_SUB(maxwell_restart_read_data)
  end function maxwell_restart_read_data

  subroutine maxwell_update_kinetic_energy(this)
    class(maxwell_t), intent(inout) :: this

    PUSH_SUB(maxwell_update_kinetic_energy)

    ! calculate Maxwell energy
    call energy_mxll_calc(this%gr, this%st, this%hm, this%hm%energy, this%st%rs_state, &
      this%st%rs_state_plane_waves)

    ! the energy of the EM wave is computed and stored in energy_mxll%energy;
    ! energy_mxll%energy = energy_mxll%e_energy + energy_mxll%b_energy
    ! here this%hm%energy is 'energy_mxll'
    this%kinetic_energy = this%hm%energy%energy

    POP_SUB(maxwell_update_kinetic_energy)
  end subroutine maxwell_update_kinetic_energy

  ! ---------------------------------------------------------
  subroutine maxwell_finalize(this)
    type(maxwell_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_finalize)

    call profiling_in(prof, "MAXWELL_FINALIZE")

    call system_end(this)

    ! free memory
    SAFE_DEALLOCATE_A(this%rs_state_init)
    SAFE_DEALLOCATE_A(this%current_density_medium)
    call hamiltonian_mxll_end(this%hm)

    call multicomm_end(this%mc)

    call states_mxll_end(this%st)

    call grid_end(this%gr)

    call restart_end(this%restart)
    call restart_end(this%restart_dump)

    call profiling_out(prof)

    POP_SUB(maxwell_finalize)
  end subroutine maxwell_finalize

end module maxwell_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
