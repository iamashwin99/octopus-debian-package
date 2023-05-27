!! Copyright (C) 2021 F. Bonafé
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

module dispersive_medium_oct_m
  use algorithm_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use current_to_mxll_field_oct_m
  use debug_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use grid_oct_m
  use linear_medium_to_em_field_oct_m
  use linear_medium_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use mxll_field_to_medium_oct_m
  use namespace_oct_m
  use parser_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use regridding_oct_m
  use space_oct_m
  use system_oct_m
  use unit_system_oct_m
  use unit_oct_m
  use varinfo_oct_m
  use write_iter_oct_m

  implicit none

  private
 public ::           &
    dispersive_medium_t,    &
    dispersive_medium_init

  type, extends(system_t) :: dispersive_medium_t
    FLOAT                 :: omega_p  !< pole frequency
    FLOAT                 :: gamma_p  !< inverse relaxation time
    FLOAT                 :: strength_p  !< pole strength
    FLOAT, allocatable    :: current_p(:,:) !< polarization current
    FLOAT, allocatable    :: e_field(:,:)
    FLOAT, allocatable    :: e_field_prev(:,:)
    integer, allocatable  :: points_map(:)
    integer               :: points_number
    FLOAT, allocatable    :: current_at_point(:,:)
    FLOAT, allocatable    :: selected_points_coordinate(:,:)
    integer               :: n_output_points
    integer               :: medium_type
    type(grid_t)          :: gr    !< the mesh
    type(multicomm_t)     :: mc    !< index and domain communicators
    type(single_medium_box_t) :: medium_box
    type(c_ptr)           :: write_handle

  contains
    procedure :: init_interaction => dispersive_medium_init_interaction
    procedure :: init_interaction_as_partner => dispersive_medium_init_interaction_as_partner
    procedure :: initial_conditions => dispersive_medium_initial_conditions
    procedure :: do_td_operation => dispersive_medium_do_td
    procedure :: is_tolerance_reached => dispersive_medium_is_tolerance_reached
    procedure :: update_quantity => dispersive_medium_update_quantity
    procedure :: update_exposed_quantity => dispersive_medium_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => dispersive_medium_copy_quantities_to_interaction
    procedure :: update_interactions_start => dispersive_medium_update_interactions_start
    procedure :: update_interactions_finish => dispersive_medium_update_interactions_finish
    procedure :: restart_write_data => dispersive_medium_restart_write_data
    procedure :: restart_read_data => dispersive_medium_restart_read_data
    procedure :: update_kinetic_energy => dispersive_medium_update_kinetic_energy
    procedure :: output_start => dispersive_medium_output_start
    procedure :: output_write => dispersive_medium_output_write
    procedure :: output_finish => dispersive_medium_output_finish
    procedure :: init_parallelization => dispersive_medium_init_parallelization
    final :: dispersive_medium_finalize
  end type dispersive_medium_t

  interface dispersive_medium_t
    procedure dispersive_medium_constructor
  end interface dispersive_medium_t

  integer, public, parameter ::      &
    DRUDE_MEDIUM                = 0

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function dispersive_medium_constructor(namespace) result(sys)
    class(dispersive_medium_t), pointer    :: sys
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(dispersive_medium_constructor)

    SAFE_ALLOCATE(sys)
    
    call dispersive_medium_init(sys, namespace)

    POP_SUB(dispersive_medium_constructor)
  end function dispersive_medium_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine dispersive_medium_init(this, namespace)
    class(dispersive_medium_t), target, intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace

    integer :: nlines, ncols, idim, il
    FLOAT :: pos(3)
    type(block_t) :: blk
    type(profile_t), save :: prof

    PUSH_SUB(dispersive_medium_init)

    this%namespace = namespace

    call profiling_in(prof, 'DISPERSIVE_MEDIUM_INIT')

    call space_init(this%space, this%namespace)
    if (this%space%is_periodic()) then
      call messages_not_implemented('Linear medium for periodic systems', namespace=namespace)
    end if
    call grid_init_stage_1(this%gr, this%namespace, this%space)

    call medium_box_init(this%medium_box, namespace)

    ! Parse electromagnetic properties of dispersive media...
    
    !%Variable MediumDispersionType
    !%Type integer
    !%Default drude_medium
    !%Section Maxwell
    !%Description
    !% Dispersion model used for the medium (only Drude model available for the moment).
    !%Option drude_medium 0
    !% Drude type of dispersion.
    !%End
    call parse_variable(namespace, 'MediumDispersionType', DRUDE_MEDIUM, this%medium_type)

    ! Parse electromagnetic properties of Dispersive media...
    !%Variable MediumPoleEnergy
    !%Type float
    !%Default 0
    !%Section Maxwell
    !%Description
    !% Energy of the pole.
    !%End
    call parse_variable(namespace, 'MediumPoleEnergy', M_ZERO, this%omega_p, unit_one/units_inp%time)

    !%Variable MediumPoleDamping
    !%Type float
    !%Default 0
    !%Section Maxwell
    !%Description
    !% Damping factor (inverse relaxation time) of the medium.
    !%End
    call parse_variable(namespace, 'MediumPoleDamping', M_ZERO, this%gamma_p, unit_one/units_inp%time)

    !%Variable MediumPoleStrength
    !%Type float
    !%Default 1.0
    !%Section Maxwell
    !%Description
    !% Strength of the pole (unitless).
    !%End
    call parse_variable(namespace, 'MediumPoleStrength', M_ONE, this%strength_p, unit_one)

    !%Variable MediumCurrentCoordinates
    !%Type block
    !%Section Maxwell
    !%Description
    !%  This allows to output phasor current vectors at particular points in space.
    !%
    !% <tt>%MediumCurrentCoordinates
    !% <br>&nbsp;&nbsp;    -1.0 | 2.0 |  4.0
    !% <br>&nbsp;&nbsp;     0.0 | 1.0 | -2.0
    !% <br>%</tt>
    !%
    !%End

    if (parse_block(namespace, 'MediumCurrentCoordinates', blk) == 0) then
      nlines = parse_block_n(blk)
      this%n_output_points = nlines
      SAFE_ALLOCATE(this%selected_points_coordinate(1:nlines,1:3))
      SAFE_ALLOCATE(this%current_at_point(1:nlines,1:3))
      do il = 1, nlines
        ncols = parse_block_cols(blk,0)
        do idim = 1, 3
          call parse_block_float(blk, il-1, idim-1, pos(idim), units_inp%length)
        end do
        this%selected_points_coordinate(il,:) = pos(:)
        this%current_at_point(il,:)  = M_ZERO
      end do
      call parse_block_end(blk)
    else
      this%n_output_points = 1
      SAFE_ALLOCATE(this%selected_points_coordinate(1,3))
      SAFE_ALLOCATE(this%current_at_point(1,3))
      this%selected_points_coordinate = M_ZERO
      this%current_at_point = M_ZERO
    end if

    call this%supported_interactions_as_partner%add(CURRENT_TO_MXLL_FIELD)
    call this%supported_interactions%add(MXLL_FIELD_TO_MEDIUM)
    this%quantities(CURRENT)%required = .true.
    this%quantities(CURRENT)%protected = .true.

    call profiling_out(prof)

    POP_SUB(dispersive_medium_init)
  end subroutine dispersive_medium_init

  ! ---------------------------------------------------------
  subroutine dispersive_medium_init_parallelization(this, grp)
    class(dispersive_medium_t),     intent(inout) :: this
    type(mpi_grp_t),      intent(in)    :: grp

    integer(i8) :: index_range(4)

    PUSH_SUB(dispersive_medium_init_parallelization)

    call system_init_parallelization(this, grp)
    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%mesh%np_global  ! Number of points in mesh
    index_range(2) = 1                      ! Number of states
    index_range(3) = 1                      ! Number of k-points
    index_range(4) = 100000                 ! Some large number

    ! create index and domain communicators
    call multicomm_init(this%mc, this%namespace, mpi_world, calc_mode_par_parallel_mask(), &
         &calc_mode_par_default_parallel_mask(), mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))
    call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc)

    call get_medium_box_points_map(this%medium_box, this%gr)
    this%points_number = this%medium_box%points_number
    SAFE_ALLOCATE(this%points_map(1:this%points_number))
    this%points_map = this%medium_box%points_map

    SAFE_ALLOCATE(this%current_p(1:this%gr%mesh%np, 1:3))
    SAFE_ALLOCATE(this%e_field(1:this%gr%mesh%np, 1:3))
    SAFE_ALLOCATE(this%e_field_prev(1:this%gr%mesh%np, 1:3))
    this%current_p(:,:) = M_ZERO

    POP_SUB(dispersive_medium_init_parallelization)
  end subroutine dispersive_medium_init_parallelization

  ! ---------------------------------------------------------
  subroutine dispersive_medium_init_interaction(this, interaction)
    class(dispersive_medium_t), target, intent(inout) :: this
    class(interaction_t),          intent(inout) :: interaction

    PUSH_SUB(dispersive_medium_init_interaction)

    select type (interaction)
    type is (mxll_field_to_medium_t)
      call interaction%init(this%gr)
    class default
      message(1) = "Trying to initialize an unsupported interaction by a Dispersive medium."
      call messages_fatal(1)
    end select

    POP_SUB(dispersive_medium_init_interaction)
  end subroutine dispersive_medium_init_interaction

  ! ---------------------------------------------------------
  subroutine dispersive_medium_init_interaction_as_partner(partner, interaction)
    class(dispersive_medium_t),       intent(in)    :: partner
    class(interaction_t),        intent(inout) :: interaction

    PUSH_SUB(dispersive_medium_init_interaction_as_partner)

    select type (interaction)
    type is (current_to_mxll_field_t)
      SAFE_ALLOCATE(interaction%partner_points_map(partner%points_number))
      SAFE_ALLOCATE(interaction%partner_current_p(partner%gr%mesh%np, partner%gr%box%dim)) 
      interaction%partner_points_map(:) = partner%points_map(1:partner%points_number)
      interaction%partner_points_number = partner%points_number
      call grid_transfer_mapping(interaction%system_gr, partner%gr, interaction%partner_to_system_map, &
           partner%namespace)
    class default
      message(1) = "Trying to initialize an unsupported interaction by a linear medium."
      call messages_fatal(1)
    end select

    POP_SUB(dispersive_medium_init_interaction_as_partner)
  end subroutine dispersive_medium_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine dispersive_medium_initial_conditions(this)
    class(dispersive_medium_t), intent(inout) :: this

    PUSH_SUB(dispersive_medium_initial_conditions)
    POP_SUB(dispersive_medium_initial_conditions)
  end subroutine dispersive_medium_initial_conditions

  ! ---------------------------------------------------------
  subroutine dispersive_medium_do_td(this, operation)
    class(dispersive_medium_t),    intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    integer :: ip, ip_in
    FLOAT :: k1(1:3), k2(1:3), k3(1:3), k4(1:3), e_field_midpoint(1:3)

    PUSH_SUB(dispersive_medium_do_td)

    ! calculation of the current using ADE
    ! \partial_t J_P(t) = - \gamma_p * J_P(t) + \epsilon_0 \omega_p^2 E(t)
    ! (Computational Electrodynamics, Taflov and Hagness, 3rd Ed., section 9.4.3, eq. 9.56c)
    ! Analysis of Units:
    ! [e/time^2*a0^2] = (1/time) * (e/time*a0^2) + (e^2/hbar*c) * (1/time)^2 * Ha / (e * a0)
    ! [e/time^2*a0^2] = (e/time^2*a0^2) + (e/hbar*c) * 1/time^2 * Ha/a0, and [c]=a0/time, so [hbar*c]=Ha*a0
    ! [e/time^2*a0^2] = (e/time^2*a0^2) + (e/time^2*a0^2)

    select case (operation%id)
    case (SKIP)
      ! Do nothing

    case (STORE_CURRENT_STATUS)
      ! For the moment we do nothing

    case (EXPMID_START)
      ! Do nothing

    case (EXPMID_FINISH)

    case (EXPMID_PREDICT_DT_2)
      this%e_field_prev(:,:) = this%e_field

      this%quantities(CURRENT)%clock = this%quantities(CURRENT)%clock + CLOCK_TICK
    case (UPDATE_HAMILTONIAN)
      ! Empty for the moment
    case (EXPMID_PREDICT_DT)

      if (this%medium_type == DRUDE_MEDIUM) then

        do ip_in = 1, this%points_number
           ip = this%points_map(ip_in)
           e_field_midpoint = (this%e_field_prev(ip,1:3) + this%e_field(ip,1:3))/M_TWO

           k1(1:3) = current_derivative(this%current_p(ip,1:3), this%e_field_prev(ip,1:3), &
                this%gamma_p, this%omega_p, this%strength_p)

           k2(1:3) = current_derivative(this%current_p(ip,1:3) + this%prop%dt * k1(1:3) / M_TWO, e_field_midpoint(1:3), &
                this%gamma_p, this%omega_p, this%strength_p)

           k3(1:3) = current_derivative(this%current_p(ip,1:3) + this%prop%dt * k2(1:3) / M_TWO, e_field_midpoint(1:3), &
                this%gamma_p, this%omega_p, this%strength_p)

           k4(1:3) = current_derivative(this%current_p(ip,1:3) + this%prop%dt * k3(1:3), this%e_field(ip,1:3), &
                this%gamma_p, this%omega_p, this%strength_p)

           this%current_p(ip,1:3) = this%current_p(ip,1:3) + this%prop%dt / CNST(6.0) * &
                (k1(1:3) + M_TWO * k2(1:3) + M_TWO * k3(1:3) + k4(1:3))
        end do
        this%quantities(CURRENT)%clock = this%quantities(CURRENT)%clock + CLOCK_TICK

      end if

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(dispersive_medium_do_td)

  contains

    function current_derivative(current_p, e_field, gamma_p, omega_p, strength_p) result(current_dot)
      FLOAT, intent(in)      :: current_p(1:3)
      FLOAT, intent(in)      :: e_field(1:3)
      FLOAT, intent(in)      :: gamma_p
      FLOAT, intent(in)      :: omega_p
      FLOAT, intent(in)      :: strength_p
      FLOAT                  :: current_dot(1:3)

      current_dot(1:3) = - gamma_p * current_p(1:3) + strength_p * P_ep * omega_p**2 * e_field(1:3)
    end function current_derivative

  end subroutine dispersive_medium_do_td

  ! ---------------------------------------------------------
  logical function dispersive_medium_is_tolerance_reached(this, tol) result(converged)
    class(dispersive_medium_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(dispersive_medium_is_tolerance_reached)

    ! this routine is never called at present, no reason to be here
    ASSERT(.false.)
    converged = .false.

    POP_SUB(dispersive_medium_is_tolerance_reached)
  end function dispersive_medium_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine dispersive_medium_update_quantity(this, iq)
    class(dispersive_medium_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(dispersive_medium_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(dispersive_medium_update_quantity)
  end subroutine dispersive_medium_update_quantity

  ! ---------------------------------------------------------
  subroutine dispersive_medium_update_exposed_quantity(partner, iq)
    class(dispersive_medium_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(dispersive_medium_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(dispersive_medium_update_exposed_quantity)
  end subroutine dispersive_medium_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine dispersive_medium_copy_quantities_to_interaction(partner, interaction)
    class(dispersive_medium_t),          intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(dispersive_medium_copy_quantities_to_interaction)

    select type (interaction)
    type is (current_to_mxll_field_t)
      interaction%partner_current_p(:,:) = partner%current_p
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(dispersive_medium_copy_quantities_to_interaction)
  end subroutine dispersive_medium_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine dispersive_medium_update_interactions_start(this)
    class(dispersive_medium_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(dispersive_medium_update_interactions_start)

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (mxll_field_to_medium_t)
        this%e_field(:,:) = interaction%system_E_field(:,:)
     end select
    end do

    POP_SUB(dispersive_medium_update_interactions_start)
  end subroutine dispersive_medium_update_interactions_start

  ! ---------------------------------------------------------
  subroutine dispersive_medium_update_interactions_finish(this)
    class(dispersive_medium_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(dispersive_medium_update_interactions_finish)

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (mxll_field_to_medium_t)
         ! Nothing to do
      end select
    end do

    POP_SUB(dispersive_medium_update_interactions_finish)
  end subroutine dispersive_medium_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine dispersive_medium_restart_write_data(this)
    class(dispersive_medium_t), intent(inout) :: this

    integer :: restart_file_unit

    PUSH_SUB(dispersive_medium_restart_write_data)

    call io_mkdir('restart/'//TD_DIR, this%namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_dispersive_medium', this%namespace, action='write')
    write(restart_file_unit, *) this%current_p(:,:)
    write(restart_file_unit, *) this%e_field(:,:)
    call io_close(restart_file_unit)

    message(1) = "Successfully wrote restart data for system "//trim(this%namespace%get())
    call messages_info(1, namespace=this%namespace)

    POP_SUB(dispersive_medium_restart_write_data)
  end subroutine dispersive_medium_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function dispersive_medium_restart_read_data(this)
    class(dispersive_medium_t), intent(inout) :: this

    integer :: restart_file_unit

    PUSH_SUB(dispersive_medium_restart_read_data)

    call io_mkdir('restart/'//TD_DIR, this%namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_dispersive_medium', this%namespace, action='read', die=.false.)
    if (restart_file_unit > 0) then
      read(restart_file_unit, *) this%current_p(:,:)
      read(restart_file_unit, *) this%e_field(:,:)
      call io_close(restart_file_unit)
      dispersive_medium_restart_read_data = .true.
    else
      ! could not open file
      dispersive_medium_restart_read_data = .false.
    end if

    if (dispersive_medium_restart_read_data) then
      message(1) = "Successfully read restart data for system "//trim(this%namespace%get())
      call messages_info(1, namespace=this%namespace)
    end if

    POP_SUB(dispersive_medium_restart_read_data)
  end function dispersive_medium_restart_read_data

  ! ---------------------------------------------------------
  subroutine dispersive_medium_update_kinetic_energy(this)
    class(dispersive_medium_t), intent(inout) :: this

    PUSH_SUB(dispersive_medium_update_kinetic_energy)

    ! TODO: evaluate proper energy associated with the current distribution
    ! For Drude model: check Giuliani/Vignale book, section 4.6.1
    this%kinetic_energy = M_ZERO

    POP_SUB(dispersive_medium_update_kinetic_energy)

  end subroutine dispersive_medium_update_kinetic_energy

  ! ---------------------------------------------------------
  subroutine dispersive_medium_output_start(this)
    class(dispersive_medium_t), intent(inout) :: this

    type(profile_t), save :: prof
    integer :: first, id, idir
    character(len=130) :: aux

    PUSH_SUB(dispersive_medium_output_start)

    call profiling_in(prof, "DISPERSIVE_MEDIUM_OUTPUT_START")

    if (this%clock%get_tick() == 0) then
      first = 0
    else
      first = this%clock%get_tick() + 1
    end if

    call io_mkdir('td.general', this%namespace)
    call write_iter_init(this%write_handle, first, units_from_atomic(units_out%time, this%prop%dt), &
         trim(io_workpath("td.general/current_at_points.dat", this%namespace)))

    if (mpi_grp_is_root(mpi_world)) then
      if (this%clock%get_tick() == 0) then
        call write_iter_clear(this%write_handle)
        call write_iter_string(this%write_handle,&
             '################################################################################')
        call write_iter_nl(this%write_handle)
        call write_iter_string(this%write_handle,'# HEADER')
        call write_iter_nl(this%write_handle)

        ! first line
        write(aux, '(a7,e20.12,3a)') '# dt = ', units_from_atomic(units_out%time, this%prop%dt), &
             " [", trim(units_abbrev(units_out%time)), "]"
        call write_iter_string(this%write_handle, aux)
        call write_iter_nl(this%write_handle)

        call write_iter_header_start(this%write_handle)

        do id = 1, this%n_output_points
          do idir = 1, 3
            write(aux, '(a,i1,a,i1,a)') 'j(', id, ',', idir, ')'
            call write_iter_header(this%write_handle, aux)
           end do
         end do

        call write_iter_nl(this%write_handle)
        call write_iter_header(this%write_handle, '#          [' // trim(units_abbrev(units_out%time)) // ']')

        !FIXME: this is not printing the proper unit to output yet, for some reason
        aux = '          [' // trim(units_abbrev(unit_one/(units_out%time*units_out%length**2))) // ']'
        do id = 1, this%n_output_points
          do idir = 1, 3
            call write_iter_header(this%write_handle, aux)
          end do
        end do
        call write_iter_nl(this%write_handle)
        call write_iter_string(this%write_handle,&
             '################################################################################')
        call write_iter_nl(this%write_handle)
      end if
    end if

    if (first == 0) call this%output_write()

    call profiling_out(prof)

    POP_SUB(dispersive_medium_output_start)
  end subroutine dispersive_medium_output_start

  ! ---------------------------------------------------------
  subroutine dispersive_medium_output_write(this)
    class(dispersive_medium_t), intent(inout) :: this

    FLOAT   :: dmin, dtmp(3)
    integer :: ip, pos_index, rankmin
    type(profile_t), save :: prof

    PUSH_SUB(dispersive_medium_output_write)

    call profiling_in(prof, "DISPERSIVE_MEDIUM_OUTPUT_WRITE")

    do ip = 1, this%n_output_points
      pos_index = mesh_nearest_point(this%gr%mesh, this%selected_points_coordinate(ip,:), dmin, rankmin)
      if (this%gr%mesh%mpi_grp%rank == rankmin) then
        dtmp(:) = this%current_p(pos_index,:)
      end if
      if (this%gr%mesh%parallel_in_domains) then
        call this%gr%mesh%mpi_grp%bcast(dtmp(:), 3, MPI_FLOAT, rankmin)
      end if
      this%current_at_point(ip,:) = units_from_atomic((unit_one/units_out%time)/(units_out%length**2), dtmp(:))
    end do

    if (.not. mpi_grp_is_root(mpi_world)) then
      POP_SUB(dispersive_medium_output_write)
      return ! only first node outputs
    end if

    call write_iter_start(this%write_handle)
    do ip = 1, this%n_output_points
      call write_iter_double(this%write_handle, this%current_at_point(ip,1:3), 3)
    end do

    call write_iter_nl(this%write_handle)
    call write_iter_flush(this%write_handle)

    call profiling_out(prof)

    POP_SUB(dispersive_medium_output_write)
  end subroutine dispersive_medium_output_write

  ! ---------------------------------------------------------
  subroutine dispersive_medium_output_finish(this)
    class(dispersive_medium_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(dispersive_medium_output_finish)

    call profiling_in(prof, "DISPERSIVE_MEDIUM_OUTPUT_FINISH")

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%write_handle)
    end if

    call profiling_out(prof)

    POP_SUB(dispersive_medium_output_finish)
  end subroutine dispersive_medium_output_finish

  ! ---------------------------------------------------------
  subroutine dispersive_medium_finalize(this)
    type(dispersive_medium_t), intent(inout) :: this

    PUSH_SUB(dispersive_medium_finalize)
    call system_end(this)
    call single_medium_box_end(this%medium_box)
    SAFE_DEALLOCATE_A(this%current_p)
    SAFE_DEALLOCATE_A(this%e_field)
    SAFE_DEALLOCATE_A(this%e_field_prev)
    SAFE_DEALLOCATE_A(this%selected_points_coordinate)
    SAFE_DEALLOCATE_A(this%current_at_point)
    SAFE_DEALLOCATE_A(this%points_map)
    call multicomm_end(this%mc)
    call grid_end(this%gr)
    POP_SUB(dispersive_medium_finalize)
  end subroutine dispersive_medium_finalize

end module dispersive_medium_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
