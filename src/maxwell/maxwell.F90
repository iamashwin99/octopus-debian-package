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
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use algorithm_oct_m
  use algorithm_factory_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use distributed_oct_m
  use external_densities_oct_m
  use field_transfer_oct_m
  use time_interpolation_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_mxll_oct_m
  use helmholtz_decomposition_m
  use index_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use iso_c_binding
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use linear_medium_to_em_field_oct_m
  use current_to_mxll_field_oct_m
  use loct_oct_m
  use lorentz_force_oct_m
  use math_oct_m
  use maxwell_boundary_op_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mesh_interpolation_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use mxll_field_to_medium_oct_m
  use namespace_oct_m
  use output_oct_m
  use output_mxll_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use propagator_mxll_oct_m
  use propagator_leapfrog_oct_m
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
  use types_oct_m
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

  !> @brief Class describing Maxwell systems
  !!
  type, extends(system_t) :: maxwell_t
    type(states_mxll_t)          :: st    !< the states
    type(hamiltonian_mxll_t)     :: hm    !< The Maxwell Hamiltonian (in Riemann-Silberstein formulation)
    type(grid_t)                 :: gr    !< the mesh
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators

    type(mesh_interpolation_t)      :: mesh_interpolate

    type(propagator_mxll_t)         :: tr_mxll   !< contains the details of the Maxwell time-evolution
    type(td_write_t)                :: write_handler
    type(c_ptr)                     :: output_handle

    CMPLX, allocatable           :: ff_rs_inhom_t1(:,:), ff_rs_inhom_t2(:,:)
    CMPLX, allocatable           :: rs_state_init(:,:)
    type(time_interpolation_t), pointer :: current_interpolation
    FLOAT                        :: bc_bounds(2,MAX_DIM), dt_bounds(2,MAX_DIM)
    integer                      :: energy_update_iter
    type(restart_t)              :: restart_dump
    type(restart_t)              :: restart

    type(lattice_vectors_t)         :: latt !< Maxwells Lattice is independent of any other system.

    type(helmholtz_decomposition_t) :: helmholtz !< Helmholtz decomposition

    logical                         :: write_previous_state = .false.

  contains
    procedure :: init_interaction => maxwell_init_interaction
    procedure :: init_parallelization => maxwell_init_parallelization
    procedure :: init_algorithm => maxwell_init_algorithm
    procedure :: initial_conditions => maxwell_initial_conditions
    procedure :: do_algorithmic_operation => maxwell_do_algorithmic_operation
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
    procedure :: get_current => maxwell_get_current
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
    this%latt = lattice_vectors_t(this%namespace, this%space)

    this%quantities(E_FIELD)%required = .true.
    this%quantities(E_FIELD)%updated_on_demand = .false.
    this%quantities(B_FIELD)%required = .true.
    this%quantities(B_FIELD)%updated_on_demand = .false.

    call mesh_interpolation_init(this%mesh_interpolate, this%gr)

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
      call interaction%init(this%gr, this%st%dim)
      call interaction%init_space_latt(this%space, this%latt)
      this%hm%current_density_from_medium = .true.
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
    integer :: ierr, ip, pos_index, rankmin
    FLOAT :: dmin

    PUSH_SUB(maxwell_init_parallelization)

    call system_init_parallelization(this, grp)

    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%np_global  ! Number of points in mesh
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

    SAFE_ALLOCATE(this%rs_state_init(1:this%gr%np_part, 1:this%st%dim))
    this%rs_state_init(:,:) = M_z0

    this%energy_update_iter = 1

    call poisson_init(this%st%poisson, this%namespace, this%space, this%gr%der, this%mc)

    call propagator_mxll_init(this%gr, this%namespace, this%st, this%hm, this%tr_mxll)
    call states_mxll_allocate(this%st, this%gr)
    call external_current_init(this%st, this%namespace, this%gr)
    this%hm%propagation_apply = .true.

    if (parse_is_defined(this%namespace, 'MaxwellIncidentWaves') .and. (this%tr_mxll%bc_plane_waves)) then
      this%st%rs_state_plane_waves(:,:) = M_z0
    end if

    ! set map for selected points
    do ip = 1, this%st%selected_points_number
      pos_index = mesh_nearest_point(this%gr, this%st%selected_points_coordinate(:,ip), dmin, rankmin)
      if (this%gr%mpi_grp%rank == rankmin) then
        this%st%selected_points_map(ip) = pos_index
      else
        this%st%selected_points_map(ip) = -1
      end if
    end do
    if (accel_is_enabled()) then
      call accel_create_buffer(this%st%buff_selected_points_map, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, &
        this%st%selected_points_number)
      call accel_write_buffer(this%st%buff_selected_points_map, this%st%selected_points_number, &
        this%st%selected_points_map)
    end if

    this%hm%plane_waves_apply = .true.
    this%hm%spatial_constant_apply = .true.

    call bc_mxll_init(this%hm%bc, this%namespace, this%space, this%gr, this%st)
    this%bc_bounds(:,1:3) = this%hm%bc%bc_bounds(:,1:3)
    call inner_and_outer_points_mapping(this%gr, this%st, this%bc_bounds)
    this%dt_bounds(2, 1:3) = this%bc_bounds(1, 1:3)
    this%dt_bounds(1, 1:3) = this%bc_bounds(1, 1:3) - this%gr%der%order * this%gr%spacing(1:3)
    call surface_grid_points_mapping(this%gr, this%st, this%dt_bounds)

    call restart_init(this%restart, this%namespace, RESTART_TD, RESTART_TYPE_LOAD, this%mc, ierr, mesh=this%gr)
    call restart_init(this%restart_dump, this%namespace, RESTART_TD, RESTART_TYPE_DUMP, this%mc, ierr, mesh=this%gr)

    ! initialize batches
    call zbatch_init(this%st%rs_stateb, 1, 1, this%st%dim, this%gr%np_part)
    if (this%st%pack_states) call this%st%rs_stateb%do_pack()
    call this%st%rs_stateb%copy_to(this%st%rs_state_prevb)
    call this%st%rs_stateb%copy_to(this%st%inhomogeneousb)
    if (this%tr_mxll%bc_plane_waves) then
      call this%st%rs_stateb%copy_to(this%st%rs_state_plane_wavesb)
    end if

    ! the order should depend on the propagator
    !this%current_interpolation => time_interpolation_t(this%gr%np, this%st%dim, 2, .true., "current")
    ! Initialize Helmholtz decomposition
    call this%helmholtz%init(this%namespace, this%gr, this%mc, this%space)

    POP_SUB(maxwell_init_parallelization)
  end subroutine maxwell_init_parallelization

  ! ---------------------------------------------------------
  subroutine maxwell_init_algorithm(this, factory)
    class(maxwell_t), intent(inout) :: this
    class(algorithm_factory_t), intent(in)    :: factory

    integer :: depth, i_interaction
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction
    character(len=256) :: label

    PUSH_SUB(maxwell_init_algorithm)

    call system_init_algorithm(this, factory)

    ! interpolation depth depends on the propagator
    select type (prop => this%algo)
    type is (propagator_exp_mid_t)
      depth = 2
    type is (propagator_leapfrog_t)
      depth = 1
    class default
      message(1) = "Propagator does not yet support initializing of interaction interpolation"
      call messages_fatal(1)
    end select
    ! set interpolation depth for interactions
    i_interaction = 0
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      select type (interaction)
      class is (current_to_mxll_field_t)
        write(label, "(A, I5.5)") "current_to_mxll_", i_interaction
        call interaction%init_interpolation(depth, trim(label), cmplx=.true.)
      end select
      i_interaction = i_interaction + 1
    end do

    POP_SUB(maxwell_init_algorithm)
  end subroutine maxwell_init_algorithm

  ! ---------------------------------------------------------
  subroutine maxwell_initial_conditions(this)
    class(maxwell_t), intent(inout) :: this

    FLOAT   :: courant

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_initial_conditions)

    call profiling_in(prof,"MAXWELL_INIT_CONDITIONS")

    courant = M_ONE/(P_c * sqrt(M_ONE/this%gr%spacing(1)**2 + M_ONE/this%gr%spacing(2)**2 + &
      M_ONE/this%gr%spacing(3)**2))

    if (this%algo%dt > M_TWO * courant) then
      write(message(1),'(a,es9.2)') 'Time step seems too large, check this value'
      call messages_warning(1, namespace=this%namespace)
    end if

    if (parse_is_defined(this%namespace, 'UserDefinedInitialMaxwellStates')) then
      call states_mxll_read_user_def(this%namespace, this%space, this%gr, this%st, this%hm%bc, this%rs_state_init)
      call messages_print_with_emphasis(msg="Setting initial EM field inside box", namespace=this%namespace)
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
        this%algo%dt, M_ZERO, this%st%rs_state, set_initial_state = .true.)
      ! for mesh parallelization, this needs communication!
      this%st%rs_state_const(:) = this%st%rs_state(mesh_global_index_from_coords(this%gr, [0,0,0]),:)
    end if

    if (parse_is_defined(this%namespace, 'UserDefinedInitialMaxwellStates')) then
      SAFE_DEALLOCATE_A(this%rs_state_init)
    end if

    call hamiltonian_mxll_update(this%hm, time = M_ZERO)

    ! calculate Maxwell energy
    call energy_mxll_calc(this%gr, this%st, this%hm, this%hm%energy, this%st%rs_state, &
      this%st%rs_state_plane_waves)

    ! Get RS states values for selected points
    call get_rs_state_at_point(this%st%selected_points_rs_state(:,:), this%st%rs_state, &
      this%st%selected_points_coordinate(:,:), this%st, this%gr)

    call mxll_set_batch(this%st%rs_stateb, this%st%rs_state, this%gr%np, this%st%dim)
    call batch_set_zero(this%st%inhomogeneousb)
    if (this%tr_mxll%bc_plane_waves) then
      call mxll_set_batch(this%st%rs_state_plane_wavesb, this%st%rs_state_plane_waves, this%gr%np, this%st%dim)
    end if

    call profiling_out(prof)

    POP_SUB(maxwell_initial_conditions)
  end subroutine maxwell_initial_conditions

  ! ---------------------------------------------------------
  logical function maxwell_do_algorithmic_operation(this, operation) result(done)
    class(maxwell_t),               intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    CMPLX, allocatable :: charge_density_ext(:)
    type(profile_t), save :: prof
    type(batch_t) :: rs_state_tmpb
    type(profile_t), save :: prof_full, prof_step, prof_inh

    PUSH_SUB(maxwell_do_algorithmic_operation)

    done = .true.
    select case (operation%id)
    case (STORE_CURRENT_STATUS)
      ! For the moment we do nothing

    case (EXPMID_START)
      SAFE_ALLOCATE(this%ff_rs_inhom_t1(1:this%gr%np_part, 1:this%hm%dim))
      SAFE_ALLOCATE(this%ff_rs_inhom_t2(1:this%gr%np_part, 1:this%hm%dim))

    case (EXPMID_FINISH)

      SAFE_DEALLOCATE_A(this%ff_rs_inhom_t1)
      SAFE_DEALLOCATE_A(this%ff_rs_inhom_t2)

    case (EXPMID_EXTRAPOLATE)
      if (this%hm%current_density_ext_flag .or. this%hm%current_density_from_medium) then
        call this%get_current(this%clock%time(), this%st%rs_current_density_t1)
        call this%get_current(this%clock%time()+this%algo%dt, this%st%rs_current_density_t2)
      end if

    case (EXPMID_PROPAGATE)

      call profiling_in(prof, "SYSTEM_MXLL_DO_TD")

      ! Propagation

      !We first compute three external charge and current densities and we convert them as RS vectors
      SAFE_ALLOCATE(charge_density_ext(1:this%gr%np))

      !No charge density at the moment
      charge_density_ext = M_z0

      call transform_rs_densities(this%hm, this%gr, charge_density_ext, &
        this%st%rs_current_density_t1, this%ff_rs_inhom_t1, RS_TRANS_FORWARD)
      call transform_rs_densities(this%hm, this%gr, charge_density_ext, &
        this%st%rs_current_density_t2, this%ff_rs_inhom_t2, RS_TRANS_FORWARD)

      SAFE_DEALLOCATE_A(charge_density_ext)

      ! Propagation dt with H_maxwell
      call mxll_propagation_step(this%hm, this%namespace, this%gr, this%space, this%st, this%tr_mxll,&
        this%st%rs_stateb, this%ff_rs_inhom_t1, this%ff_rs_inhom_t2, this%clock%time(), this%algo%dt)

      this%quantities(E_FIELD)%clock = this%quantities(E_FIELD)%clock + CLOCK_TICK
      this%quantities(B_FIELD)%clock = this%quantities(B_FIELD)%clock + CLOCK_TICK

      call profiling_out(prof)

    case (LEAPFROG_START)
      if (any(this%hm%bc%bc_ab_type == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML)) then
        call bc_mxll_initialize_pml_simple(this%hm%bc, this%space, this%gr, this%algo%dt)
      end if

    case (LEAPFROG_FINISH)

    case (LEAPFROG_PROPAGATE)
      call profiling_in(prof_full, "LEAPFROG_PROPAGATE")

      call profiling_in(prof_inh, "LEAPFROG_INHOMOGENEOUS")
      if (this%hm%current_density_ext_flag .or. this%hm%current_density_from_medium) then
        call this%get_current(this%clock%time(), this%st%rs_current_density_t1)
        call mxll_set_batch(this%st%inhomogeneousb, this%st%rs_current_density_t1, this%gr%np, this%st%dim)
        call batch_scal(this%gr%np, -M_ONE, this%st%inhomogeneousb)
      end if
      call profiling_out(prof_inh)

      call this%st%rs_stateb%copy_to(rs_state_tmpb)

      call hamiltonian_mxll_update(this%hm, this%algo%clock%time())

      ! do boundaries at the beginning
      call mxll_apply_boundaries(this%tr_mxll, this%st, this%hm, this%gr, this%namespace, this%algo%clock%time(), &
        this%algo%dt, this%st%rs_stateb)

      ! apply hamiltonian
      call hamiltonian_mxll_apply_simple(this%hm, this%namespace, this%gr, this%st%rs_stateb, rs_state_tmpb)

      ! add inhomogeneous terms
      call batch_xpay(this%gr%np, this%st%inhomogeneousb, M_ONE, rs_state_tmpb)

      call profiling_in(prof_step, "LEAPFROG_STEP")

      if (this%clock%get_tick() == 0) then
        ! for the first step, we do one forward Euler step
        call batch_xpay(this%gr%np, this%st%rs_stateb, this%algo%dt, rs_state_tmpb)
      else
        ! the leapfrog step depends on the previous state
        call batch_xpay(this%gr%np, this%st%rs_state_prevb, M_TWO*this%algo%dt, rs_state_tmpb)
      end if

      ! save the current rs state
      call this%st%rs_stateb%copy_data_to(this%gr%np, this%st%rs_state_prevb)
      ! update to new timestep
      call rs_state_tmpb%copy_data_to(this%gr%np, this%st%rs_stateb)

      call profiling_out(prof_step)

      ! update PML convolution values
      if (any(this%hm%bc%bc_ab_type == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML)) then
        call mxll_update_pml_simple(this%hm, this%st%rs_stateb)
      end if

      call rs_state_tmpb%end()

      call profiling_out(prof_full)

      this%quantities(E_FIELD)%clock = this%quantities(E_FIELD)%clock + CLOCK_TICK
      this%quantities(B_FIELD)%clock = this%quantities(B_FIELD)%clock + CLOCK_TICK

    case (STEP_DONE)
      call maxwell_exec_end_of_timestep_tasks(this)
      done = .false.

    case default
      done = .false.
    end select

    POP_SUB(maxwell_do_algorithmic_operation)
  end function maxwell_do_algorithmic_operation

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

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(this%quantities(iq)%updated_on_demand)

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

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(partner%quantities(iq)%updated_on_demand)

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
      call interaction%init_from_partner(partner%gr, partner%space, partner%namespace)
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
    FLOAT, allocatable :: b_field(:,:), vec_pot(:,:)

    PUSH_SUB(maxwell_copy_quantities_to_interaction)

    call profiling_in(prof, "MXLL_CPY_QUANTITIES_INT")

    select type (interaction)
    type is (lorentz_force_t)
      call mxll_get_batch(partner%st%rs_stateb, partner%st%rs_state, partner%gr%np, partner%st%dim)

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
      call mxll_get_batch(partner%st%rs_stateb, partner%st%rs_state, partner%gr%np, partner%st%dim)
      select case (interaction%type)

      case (MXLL_FIELD_TOTAL)
        call get_electric_field_state(partner%st%rs_state, partner%gr, interaction%partner_field)

      case (MXLL_FIELD_TRANS)
        call get_transverse_rs_state(partner%helmholtz, partner%st, partner%namespace)
        call get_electric_field_state(partner%st%rs_state_trans, partner%gr, interaction%partner_field)

      case (MXLL_FIELD_LONG)
        call partner%helmholtz%get_long_field(partner%namespace, partner%st%rs_state_long, total_field=partner%st%rs_state)
        call get_electric_field_state(partner%st%rs_state_long, partner%gr, interaction%partner_field)

      case (MXLL_VEC_POT_TRANS)
        SAFE_ALLOCATE(b_field(1:partner%gr%np_part, 1:partner%gr%box%dim))
        SAFE_ALLOCATE(vec_pot(1:partner%gr%np_part, 1:partner%gr%box%dim))
        ! magnetic field is always transverse
        call get_magnetic_field_state(partner%st%rs_state, partner%gr, partner%st%rs_sign, b_field, &
          partner%st%mu(1:partner%gr%np), partner%gr%np)
        ! vector potential stored in partner_field
        call partner%helmholtz%get_vector_potential(partner%namespace, vec_pot, b_field)
        ! in the convention used by the electronic system, the -1/c factor is included in the vector potential
        interaction%partner_field(1:partner%gr%np,1:partner%gr%box%dim) = - M_ONE / P_c * &
          vec_pot(1:partner%gr%np,1:partner%gr%box%dim)
        SAFE_DEALLOCATE_A(b_field)
        SAFE_DEALLOCATE_A(vec_pot)

      case default
        message(1) = "Unknown type of field requested by interaction."
        call messages_fatal(1, namespace=partner%namespace)
      end select

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

    call mxll_get_batch(this%st%rs_stateb, this%st%rs_state, this%gr%np, this%st%dim)

    call td_write_mxll_init(this%write_handler, this%namespace, this%clock%get_tick(), this%algo%dt)
    if (this%st%fromScratch) then
      call td_write_mxll_iter(this%write_handler, this%space, this%gr, this%st, this%hm, this%helmholtz, this%algo%dt, &
        this%clock%get_tick(), this%namespace)
      call td_write_mxll_free_data(this%write_handler, this%namespace, this%space, this%gr, this%st, this%hm, this%helmholtz, &
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

    call td_write_mxll_iter(this%write_handler, this%space, this%gr, this%st, this%hm, this%helmholtz, this%algo%dt, &
      this%clock%get_tick(), this%namespace)

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
      call mxll_get_batch(this%st%rs_stateb, this%st%rs_state, this%gr%np, this%st%dim)

      call td_write_mxll_free_data(this%write_handler, this%namespace, this%space, this%gr, this%st, this%hm, this%helmholtz, &
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

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (linear_medium_to_em_field_t)
        if (allocated(this%hm%medium_boxes) .and. .not. this%hm%medium_boxes_initialized) then
          iint = iint + 1
          this%hm%medium_boxes(iint) = interaction%medium_box
        end if
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

    integer :: ierr, err, zff_dim, id, id1, id2, ip_in, offset, iout
    logical :: pml_check, write_previous_state
    CMPLX, allocatable :: zff(:,:)
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction


    PUSH_SUB(maxwell_restart_write_data)
    ierr = 0

    if (mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_MAXWELL_MAX
        if (this%write_handler%out(iout)%write) then
          call write_iter_flush(this%write_handler%out(iout)%handle)
        end if
      end do
    end if

    if (.not. restart_skip(this%restart_dump)) then
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
      select type (prop => this%algo)
      type is (propagator_leapfrog_t)
        write_previous_state = .true.
        zff_dim = zff_dim + this%st%dim
      class  default
        write_previous_state = .false.
      end select
      if (pml_check .and. accel_is_enabled()) then
        call accel_read_buffer(this%hm%bc%pml%buff_conv_plus, &
          int(this%hm%bc%pml%points_number, i8)*3*3, this%hm%bc%pml%conv_plus)
        call accel_read_buffer(this%hm%bc%pml%buff_conv_minus, &
          int(this%hm%bc%pml%points_number, i8)*3*3, this%hm%bc%pml%conv_minus)
      end if

      call mxll_get_batch(this%st%rs_stateb, this%st%rs_state, this%gr%np, this%st%dim)
      if (write_previous_state) then
        call mxll_get_batch(this%st%rs_state_prevb, this%st%rs_state_prev, this%gr%np, this%st%dim)
      end if

      SAFE_ALLOCATE(zff(1:this%gr%np,1:zff_dim))
      zff = M_z0

      zff(1:this%gr%np, 1:this%st%dim) = this%st%rs_state(1:this%gr%np, 1:this%st%dim)
      if (this%tr_mxll%bc_plane_waves) then
        call mxll_get_batch(this%st%rs_state_plane_wavesb, this%st%rs_state_plane_waves, &
          this%gr%np, this%st%dim)
        zff(1:this%gr%np, this%st%dim+1:this%st%dim+this%st%dim) = &
          this%st%rs_state_plane_waves(1:this%gr%np, 1:this%st%dim)
        offset = 2*this%st%dim
      else
        offset = this%st%dim
      end if
      if (pml_check) then
        id = 0
        do id1 = 1, 3
          do id2 = 1, 3
            id = id + 1
            do ip_in = 1, this%hm%bc%pml%points_number
              zff(ip_in, offset+id) = this%hm%bc%pml%conv_plus(ip_in, id1, id2)
              zff(ip_in, offset+9+id) = this%hm%bc%pml%conv_minus(ip_in, id1, id2)
            end do
          end do
        end do
        offset = offset + 18
      end if
      if (write_previous_state) then
        zff(1:this%gr%np, offset+1:offset+this%st%dim) = &
          this%st%rs_state_prev(1:this%gr%np, 1:this%st%dim)
      end if

      call states_mxll_dump(this%restart_dump, this%st, this%space, this%gr, zff, zff_dim, err, this%clock%get_tick())
      if (err /= 0) ierr = ierr + 1

      if (this%hm%current_density_from_medium) then
        !call this%current_interpolation%write_restart(this%gr, this%space, this%restart_dump, err)
        call iter%start(this%interactions)
        do while (iter%has_next())
          interaction => iter%get_next()
          select type (interaction)
          class is (current_to_mxll_field_t)
            call interaction%write_restart(this%gr, this%space, this%restart_dump, err)
          end select
        end do
        if (err /= 0) ierr = ierr + 1
      end if

      if (debug%info) then
        message(1) = "Debug: Writing td_maxwell restart done."
        call messages_info(1, namespace=this%namespace)
      end if

      SAFE_DEALLOCATE_A(zff)

      if (ierr /=0) then
        message(1) = "Unable to write time-dependent Maxwell restart information."
        call messages_warning(1, namespace=this%namespace)
      end if
    end if

    POP_SUB(maxwell_restart_write_data)
  end subroutine maxwell_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function maxwell_restart_read_data(this)
    class(maxwell_t), intent(inout) :: this

    integer :: ierr, err, zff_dim, id, id1, id2, ip_in, offset
    logical :: pml_check, read_previous_state
    CMPLX, allocatable :: zff(:,:)
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

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
      select type (prop => this%algo)
      type is (propagator_leapfrog_t)
        read_previous_state = .true.
        zff_dim = zff_dim + this%st%dim
      class  default
        read_previous_state = .false.
      end select

      SAFE_ALLOCATE(zff(1:this%gr%np,1:zff_dim))

      call states_mxll_load(this%restart, this%st, this%gr, this%namespace, this%space, zff, &
        zff_dim, err, label = ": td_maxwell")
      this%st%rs_current_density_restart = .true.

      this%st%rs_state(1:this%gr%np,1:this%st%dim) = zff(1:this%gr%np, 1:this%st%dim)
      if (this%tr_mxll%bc_plane_waves) then
        this%st%rs_state_plane_waves(1:this%gr%np,1:this%st%dim) = &
          zff(1:this%gr%np,this%st%dim+1:this%st%dim+3)
        offset = 2*this%st%dim
      else
        offset = this%st%dim
      end if
      if (pml_check) then
        id = 0
        do id1 = 1, 3
          do id2 = 1, 3
            id = id+1
            do ip_in = 1, this%hm%bc%pml%points_number
              this%hm%bc%pml%conv_plus(ip_in,id1,id2)  = zff(ip_in, offset+  id)
              this%hm%bc%pml%conv_minus(ip_in,id1,id2) = zff(ip_in, offset+9+id)
            end do
          end do
        end do
        this%hm%bc%pml%conv_plus_old = this%hm%bc%pml%conv_plus
        this%hm%bc%pml%conv_minus_old = this%hm%bc%pml%conv_minus
        offset = offset + 18
      end if
      if (read_previous_state) then
        this%st%rs_state_prev(1:this%gr%np, 1:this%st%dim) = &
          zff(1:this%gr%np, offset+1:offset+this%st%dim)
      end if

      if (err /= 0) then
        ierr = ierr + 1
      end if

      if (debug%info) then
        message(1) = "Debug: Reading td restart done."
        call messages_info(1, namespace=this%namespace)
      end if
      SAFE_DEALLOCATE_A(zff)

      if (pml_check .and. accel_is_enabled()) then
        call accel_write_buffer(this%hm%bc%pml%buff_conv_plus, &
          int(this%hm%bc%pml%points_number, i8)*3*3, this%hm%bc%pml%conv_plus)
        call accel_write_buffer(this%hm%bc%pml%buff_conv_minus, &
          int(this%hm%bc%pml%points_number, i8)*3*3, this%hm%bc%pml%conv_minus)
        call accel_write_buffer(this%hm%bc%pml%buff_conv_plus_old, &
          int(this%hm%bc%pml%points_number, i8)*3*3, this%hm%bc%pml%conv_plus_old)
      end if

      if (this%hm%current_density_from_medium) then
        !call this%current_interpolation%read_restart(this%gr, this%space, this%restart, err)
        call iter%start(this%interactions)
        do while (iter%has_next())
          interaction => iter%get_next()
          select type (interaction)
          class is (current_to_mxll_field_t)
            call interaction%read_restart(this%gr, this%space, this%restart, err)
          end select
        end do
        if (err /= 0) then
          ierr = ierr + 1
        end if
      end if

      ! set batches
      call mxll_set_batch(this%st%rs_stateb, this%st%rs_state, this%gr%np, this%st%dim)
      if (read_previous_state) then
        call mxll_set_batch(this%st%rs_state_prevb, this%st%rs_state_prev, this%gr%np, this%st%dim)
      end if
      call batch_set_zero(this%st%inhomogeneousb)
      if (this%tr_mxll%bc_plane_waves) then
        call mxll_set_batch(this%st%rs_state_plane_wavesb, this%st%rs_state_plane_waves, this%gr%np, this%st%dim)
      end if

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

    ! the energy has already been computed at the end of the timestep

    ! the energy of the EM wave is computed and stored in energy_mxll%energy;
    ! energy_mxll%energy = energy_mxll%e_energy + energy_mxll%b_energy
    ! here this%hm%energy is 'energy_mxll'
    this%kinetic_energy = this%hm%energy%energy

    POP_SUB(maxwell_update_kinetic_energy)
  end subroutine maxwell_update_kinetic_energy

  ! ---------------------------------------------------------
  subroutine maxwell_exec_end_of_timestep_tasks(this)
    class(maxwell_t), intent(inout) :: this

    PUSH_SUB(maxwell_exec_end_of_timestep_tasks)

    ! calculate Maxwell energy
    call energy_mxll_calc_batch(this%gr, this%st, this%hm, this%hm%energy, this%st%rs_stateb, this%st%rs_state_plane_wavesb)

    ! get RS state values for selected points
    call get_rs_state_batch_selected_points(this%st%selected_points_rs_state(:,:), this%st%rs_stateb, &
      this%st)

    POP_SUB(maxwell_exec_end_of_timestep_tasks)
  end subroutine maxwell_exec_end_of_timestep_tasks

  ! ---------------------------------------------------------
  subroutine maxwell_get_current(this, time, current)
    class(maxwell_t), intent(inout) :: this
    FLOAT,            intent(in)    :: time
    CMPLX,            intent(inout) :: current(:, :)

    type(interaction_iterator_t) :: iter
    CMPLX, allocatable :: current_density_ext(:, :)

    PUSH_SUB(maxwell_get_current)

    SAFE_ALLOCATE(current_density_ext(1:this%gr%np, 1:this%st%dim))
    current = M_z0
    if (this%hm%current_density_from_medium) then
      ! interpolate current from interaction
      call iter%start(this%interactions)
      do while (iter%has_next())
        select type (interaction => iter%get_next())
        class is (current_to_mxll_field_t)
          call interaction%interpolate(time, current_density_ext)
          call lalg_axpy(this%gr%np, 3, M_ONE, current_density_ext, current)
        end select
      end do
    end if
    ! calculation of external RS density
    if (this%hm%current_density_ext_flag) then
      call get_rs_density_ext(this%st, this%gr, time, current_density_ext)
      call lalg_axpy(this%gr%np, 3, M_ONE, current_density_ext, current)
    end if
    SAFE_DEALLOCATE_A(current_density_ext)

    POP_SUB(maxwell_get_current)
  end subroutine maxwell_get_current

  ! ---------------------------------------------------------
  subroutine maxwell_finalize(this)
    type(maxwell_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_finalize)

    call profiling_in(prof, "MAXWELL_FINALIZE")

    call system_end(this)

    ! free memory
    SAFE_DEALLOCATE_A(this%rs_state_init)
    call hamiltonian_mxll_end(this%hm)

    call multicomm_end(this%mc)

    call this%st%rs_stateb%end()
    call this%st%rs_state_prevb%end()
    call this%st%inhomogeneousb%end()
    if (this%tr_mxll%bc_plane_waves) then
      call this%st%rs_state_plane_wavesb%end()
    end if

    call states_mxll_end(this%st)

    call grid_end(this%gr)

    call restart_end(this%restart)
    call restart_end(this%restart_dump)

    call poisson_end(this%st%poisson)

    call profiling_out(prof)

    POP_SUB(maxwell_finalize)
  end subroutine maxwell_finalize

end module maxwell_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
