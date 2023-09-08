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
  use algorithm_factory_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use current_to_mxll_field_oct_m
  use debug_oct_m
  use field_transfer_oct_m
  use external_source_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use interaction_partner_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use grid_oct_m
  use lalg_basic_oct_m
  use linear_medium_to_em_field_oct_m
  use linear_medium_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use mxll_field_to_medium_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use propagator_rk4_oct_m
  use propagator_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use regridding_oct_m
  use restart_oct_m
  use space_oct_m
  use system_oct_m
  use unit_system_oct_m
  use unit_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::           &
    dispersive_medium_t,    &
    dispersive_medium_init

  !> @brief dispersive medium for classical electrodynamics calculations
  type, extends(system_t) :: dispersive_medium_t
    private
    FLOAT                 :: omega_p  !< pole frequency
    FLOAT                 :: gamma_p  !< inverse relaxation time
    FLOAT                 :: strength_p  !< pole strength
    FLOAT, allocatable    :: current_p(:,:) !< polarization current
    FLOAT, allocatable    :: e_field(:,:)
    FLOAT, allocatable    :: e_field_dt_half(:,:)
    FLOAT, allocatable    :: e_field_dt_full(:,:)
    FLOAT, allocatable    :: current_at_point(:,:)
    FLOAT, allocatable    :: selected_points_coordinate(:,:)
    integer               :: n_output_points
    integer               :: medium_type
    type(grid_t)          :: gr    !< the mesh
    type(multicomm_t)     :: mc    !< index and domain communicators
    type(c_ptr)           :: write_handle
    type(output_t)        :: outp  !< the output
    logical               :: from_scratch = .true.
    type(restart_t)       :: restart_load
    type(restart_t)       :: restart_dump

  contains
    procedure :: init_interaction => dispersive_medium_init_interaction
    procedure :: init_interaction_as_partner => dispersive_medium_init_interaction_as_partner
    procedure :: init_algorithm => dispersive_medium_init_algorithm
    procedure :: initial_conditions => dispersive_medium_initial_conditions
    procedure :: do_algorithmic_operation => dispersive_medium_do_algorithmic_operation
    procedure :: is_tolerance_reached => dispersive_medium_is_tolerance_reached
    procedure :: update_quantity => dispersive_medium_update_quantity
    procedure :: update_exposed_quantity => dispersive_medium_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => dispersive_medium_copy_quantities_to_interaction
    procedure :: restart_write_data => dispersive_medium_restart_write_data
    procedure :: restart_read_data => dispersive_medium_restart_read_data
    procedure :: update_kinetic_energy => dispersive_medium_update_kinetic_energy
    procedure :: output_start => dispersive_medium_output_start
    procedure :: output_write => dispersive_medium_output_write
    procedure :: output_finish => dispersive_medium_output_finish
    procedure :: init_parallelization => dispersive_medium_init_parallelization
    procedure :: get_efield => dispersive_medium_get_efield
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
    type(namespace_t),          intent(in) :: namespace

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

    call profiling_in(prof, 'DISP_MEDIUM_INIT')

    call space_init(this%space, this%namespace)
    if (this%space%is_periodic()) then
      call messages_not_implemented('Linear medium for periodic systems', namespace=namespace)
    end if
    if (this%space%dim /= 3) then
      call messages_not_implemented('Linear medium for dimensions other than 3', namespace=namespace)
    end if
    call grid_init_stage_1(this%gr, this%namespace, this%space)

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
    if (.not. varinfo_valid_option('MediumDispersionType', this%medium_type)) then
      call messages_input_error(namespace, 'MediumDispersionType')
    end if

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
    this%quantities(CURRENT)%updated_on_demand = .false.

    call profiling_out(prof)

    POP_SUB(dispersive_medium_init)
  end subroutine dispersive_medium_init

  ! ---------------------------------------------------------
  subroutine dispersive_medium_init_parallelization(this, grp)
    class(dispersive_medium_t),     intent(inout) :: this
    type(mpi_grp_t),      intent(in)    :: grp

    integer(i8) :: index_range(4)
    integer :: ierr

    PUSH_SUB(dispersive_medium_init_parallelization)

    call system_init_parallelization(this, grp)
    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%np_global  ! Number of points in mesh
    index_range(2) = 1                      ! Number of states
    index_range(3) = 1                      ! Number of k-points
    index_range(4) = 100000                 ! Some large number

    ! create index and domain communicators
    call multicomm_init(this%mc, this%namespace, mpi_world, calc_mode_par_parallel_mask(), &
    &calc_mode_par_default_parallel_mask(), mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))
    call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc)

    call restart_init(this%restart_dump, this%namespace, RESTART_TD, RESTART_TYPE_DUMP, &
      this%mc, ierr, mesh=this%gr)
    call restart_init(this%restart_load, this%namespace, RESTART_TD, RESTART_TYPE_LOAD, &
      this%mc, ierr, mesh=this%gr)

    SAFE_ALLOCATE(this%current_p(1:this%gr%np, 1:3))
    SAFE_ALLOCATE(this%e_field(1:this%gr%np, 1:3))
    this%e_field(:,:) = M_ZERO

    POP_SUB(dispersive_medium_init_parallelization)
  end subroutine dispersive_medium_init_parallelization

  ! ---------------------------------------------------------
  subroutine dispersive_medium_init_algorithm(this, factory)
    class(dispersive_medium_t), intent(inout) :: this
    class(algorithm_factory_t), intent(in)    :: factory

    integer :: depth, i_interaction
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction
    character(len=256) :: label

    PUSH_SUB(dispersive_medium_init_algorithm)

    call system_init_algorithm(this, factory)

    ! interpolation depth depends on the propagator
    select type (prop => this%algo)
    type is (propagator_rk4_t)
      ! TODO: should be five, two is only for backwards compatibility
      depth = 2
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
      class is (field_transfer_t)
        write(label, "(A, I5.5)") "field_transfer_medium_", i_interaction
        call interaction%init_interpolation(depth, trim(label))
      end select
      i_interaction = i_interaction + 1
    end do

    POP_SUB(dispersive_medium_init_algorithm)
  end subroutine dispersive_medium_init_algorithm

  ! ---------------------------------------------------------
  subroutine dispersive_medium_init_interaction(this, interaction)
    class(dispersive_medium_t), target, intent(inout) :: this
    class(interaction_t),          intent(inout) :: interaction

    PUSH_SUB(dispersive_medium_init_interaction)

    select type (interaction)
    type is (mxll_field_to_medium_t)
      call interaction%init(this%gr, 3)
      interaction%type = MXLL_FIELD_TOTAL
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
      call interaction%init_from_partner(partner%gr, partner%space, partner%namespace)
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

    this%from_scratch = .true.
    this%current_p(:,:) = M_ZERO

    POP_SUB(dispersive_medium_initial_conditions)
  end subroutine dispersive_medium_initial_conditions

  ! ---------------------------------------------------------
  logical function dispersive_medium_do_algorithmic_operation(this, operation) result(done)
    class(dispersive_medium_t),     intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    integer :: ip
    FLOAT :: k1(1:3), k2(1:3), k3(1:3), k4(1:3)

    PUSH_SUB(dispersive_medium_do_algorithmic_operation)

    ! calculation of the current using ADE
    ! \partial_t J_P(t) = - \gamma_p * J_P(t) + \epsilon_0 \omega_p^2 E(t)
    ! (Computational Electrodynamics, Taflov and Hagness, 3rd Ed., section 9.4.3, eq. 9.56c)
    ! Analysis of Units:
    ! [e/time^2*a0^2] = (1/time) * (e/time*a0^2) + (e^2/hbar*c) * (1/time)^2 * Ha / (e * a0)
    ! [e/time^2*a0^2] = (e/time^2*a0^2) + (e/hbar*c) * 1/time^2 * Ha/a0, and [c]=a0/time, so [hbar*c]=Ha*a0
    ! [e/time^2*a0^2] = (e/time^2*a0^2) + (e/time^2*a0^2)

    done = .true.
    select case (operation%id)
    case (STORE_CURRENT_STATUS)
      ! For the moment we do nothing

    case (RK4_START)
      SAFE_ALLOCATE(this%e_field_dt_half(1:this%gr%np, 1:3))
      SAFE_ALLOCATE(this%e_field_dt_full(1:this%gr%np, 1:3))

    case (RK4_FINISH)
      SAFE_DEALLOCATE_A(this%e_field_dt_half)
      SAFE_DEALLOCATE_A(this%e_field_dt_full)

    case (RK4_EXTRAPOLATE)
      call this%get_efield(this%clock%time(), this%e_field)
      call this%get_efield(this%clock%time()+this%algo%dt/M_TWO, this%e_field_dt_half)
      call this%get_efield(this%clock%time()+this%algo%dt, this%e_field_dt_full)

    case (RK4_PROPAGATE)
      !$omp parallel do private(k1, k2, k3, k4)
      do ip = 1, this%gr%np
        k1(1:3) = current_derivative(this%current_p(ip,1:3), &
          this%e_field(ip,1:3), this%gamma_p, this%omega_p, this%strength_p)

        k2(1:3) = current_derivative(this%current_p(ip,1:3) + this%algo%dt * k1(1:3) / M_TWO, &
          this%e_field_dt_half(ip,1:3), this%gamma_p, this%omega_p, this%strength_p)

        k3(1:3) = current_derivative(this%current_p(ip,1:3) + this%algo%dt * k2(1:3) / M_TWO, &
          this%e_field_dt_half(ip,1:3), this%gamma_p, this%omega_p, this%strength_p)

        k4(1:3) = current_derivative(this%current_p(ip,1:3) + this%algo%dt * k3(1:3), &
          this%e_field_dt_full(ip,1:3), this%gamma_p, this%omega_p, this%strength_p)

        this%current_p(ip,1:3) = this%current_p(ip,1:3) + this%algo%dt / CNST(6.0) * &
          (k1(1:3) + M_TWO * k2(1:3) + M_TWO * k3(1:3) + k4(1:3))
      end do
      this%quantities(CURRENT)%clock = this%quantities(CURRENT)%clock + CLOCK_TICK

    case default
      done = .false.
    end select

    POP_SUB(dispersive_medium_do_algorithmic_operation)
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

  end function dispersive_medium_do_algorithmic_operation

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

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(this%quantities(iq)%updated_on_demand)

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

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(partner%quantities(iq)%updated_on_demand)

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
      interaction%partner_field(:,:) = partner%current_p
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(dispersive_medium_copy_quantities_to_interaction)
  end subroutine dispersive_medium_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine dispersive_medium_restart_write_data(this)
    class(dispersive_medium_t), intent(inout) :: this

    character(len=256) :: filename
    integer :: idir, err
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(dispersive_medium_restart_write_data)

    if (.not. restart_skip(this%restart_dump)) then
      do idir = 1, this%space%dim
        write(filename, "(A,A)") "current_p-", index2axis(idir)
        call drestart_write_mesh_function(this%restart_dump, this%space, trim(filename), this%gr, &
          this%current_p(:, idir), err)
      end do

      call iter%start(this%interactions)
      do while (iter%has_next())
        interaction => iter%get_next()
        select type (interaction)
        class is (mxll_field_to_medium_t)
          call interaction%write_restart(this%gr, this%space, this%restart_dump, err)
        end select
      end do

      if (err == 0) then
        message(1) = "Successfully wrote restart data for system "//trim(this%namespace%get())
        call messages_info(1, namespace=this%namespace)
      end if
    end if

    POP_SUB(dispersive_medium_restart_write_data)
  end subroutine dispersive_medium_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function dispersive_medium_restart_read_data(this)
    class(dispersive_medium_t), intent(inout) :: this

    character(len=256) :: filename
    integer :: idir, err
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(dispersive_medium_restart_read_data)

    dispersive_medium_restart_read_data = .false.
    if (.not. restart_skip(this%restart_load)) then
      do idir = 1, 3
        write(filename, "(A,A)") "current_p-", index2axis(idir)
        call drestart_read_mesh_function(this%restart_load, this%space, trim(filename), this%gr, &
          this%current_p(:, idir), err)
      end do

      call iter%start(this%interactions)
      do while (iter%has_next())
        interaction => iter%get_next()
        select type (interaction)
        class is (mxll_field_to_medium_t)
          call interaction%read_restart(this%gr, this%space, this%restart_load, err)
        end select
      end do

      if (err == 0) then
        dispersive_medium_restart_read_data = .true.
        this%from_scratch = .false.
      else
        ! set to 0 again in case this was read incompletely
        this%current_p(:, :) = M_ZERO
      end if
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

    call profiling_in(prof, "DISP_MEDIUM_OUTPUT_START")

    if (this%clock%get_tick() == 0) then
      first = 0
    else
      first = this%clock%get_tick() + 1
    end if

    call io_mkdir('td.general', this%namespace)
    call write_iter_init(this%write_handle, first, units_from_atomic(units_out%time, this%algo%dt), &
      trim(io_workpath("td.general/current_at_points", this%namespace)))

    if (mpi_grp_is_root(mpi_world)) then
      if (this%clock%get_tick() == 0) then
        call write_iter_clear(this%write_handle)
        call write_iter_string(this%write_handle,&
          '################################################################################')
        call write_iter_nl(this%write_handle)
        call write_iter_string(this%write_handle,'# HEADER')
        call write_iter_nl(this%write_handle)

        ! first line
        write(aux, '(a7,e20.12,3a)') '# dt = ', units_from_atomic(units_out%time, this%algo%dt), &
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

    call profiling_in(prof, "DISP_MEDIUM_OUTPUT_WRITE")

    do ip = 1, this%n_output_points
      pos_index = mesh_nearest_point(this%gr, this%selected_points_coordinate(ip,:), dmin, rankmin)
      if (this%gr%mpi_grp%rank == rankmin) then
        dtmp(:) = this%current_p(pos_index,:)
      end if
      if (this%gr%parallel_in_domains) then
        call this%gr%mpi_grp%bcast(dtmp(:), 3, MPI_FLOAT, rankmin)
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

    call profiling_in(prof, "DISP_MEDIUM_OUTPUT_FINISH")

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%write_handle)
    end if

    call profiling_out(prof)

    POP_SUB(dispersive_medium_output_finish)
  end subroutine dispersive_medium_output_finish

  ! ---------------------------------------------------------
  subroutine dispersive_medium_get_efield(this, time, efield)
    class(dispersive_medium_t), intent(inout) :: this
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(inout) :: efield(:, :)

    type(interaction_iterator_t) :: iter
    FLOAT, allocatable :: efield_tmp(:, :)

    PUSH_SUB(dispersive_medium_get_efield)

    SAFE_ALLOCATE(efield_tmp(1:this%gr%np, 1:3))
    efield = M_z0
    ! interpolate efield from interaction
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (mxll_field_to_medium_t)
        call interaction%interpolate(time, efield_tmp)
        call lalg_axpy(this%gr%np, 3, M_ONE, efield_tmp, efield)
      end select
    end do
    SAFE_DEALLOCATE_A(efield_tmp)

    POP_SUB(dispersive_medium_get_efield)
  end subroutine dispersive_medium_get_efield

  ! ---------------------------------------------------------
  subroutine dispersive_medium_finalize(this)
    type(dispersive_medium_t), intent(inout) :: this

    PUSH_SUB(dispersive_medium_finalize)
    call system_end(this)
    SAFE_DEALLOCATE_A(this%current_p)
    SAFE_DEALLOCATE_A(this%e_field)
    SAFE_DEALLOCATE_A(this%selected_points_coordinate)
    SAFE_DEALLOCATE_A(this%current_at_point)
    call multicomm_end(this%mc)
    call grid_end(this%gr)
    POP_SUB(dispersive_medium_finalize)
  end subroutine dispersive_medium_finalize

end module dispersive_medium_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
