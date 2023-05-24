!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2019 N. Tancogne-Dejean
!! Copyright (C) 2020 M. Oliveira
!! Copyright (C) 2021 S. Ohlmann, I. Albar, A. Obzhirov, A. Geondzhian, M. Lawan
!! Copyright (C) 2022 M. Oliveira
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

module classical_particles_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
  use force_interaction_oct_m
  use global_oct_m
  use interaction_oct_m
  use io_oct_m
  use lalg_adv_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_data_classical_particles_oct_m
  use propagator_beeman_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m

  implicit none

  private
  public ::                                      &
    classical_particles_t,                       &
    classical_particles_init,                    &
    classical_particles_copy,                    &
    classical_particles_end,                     &
    classical_particles_init_interaction,        &
    classical_particles_update_quantity,         &
    classical_particles_update_exposed_quantity, &
    classical_particles_init_interaction_as_partner

  type, extends(system_t), abstract :: classical_particles_t
    private
    integer, public :: np                        !< Number of particles in the system
    FLOAT, allocatable, public :: mass(:)        !< Mass of the particles
    FLOAT, allocatable, public :: pos(:,:)       !< Position of the particles
    FLOAT, allocatable, public :: vel(:,:)       !< Velocity of the particles
    FLOAT, allocatable, public :: tot_force(:,:) !< Total force acting on each particle
    logical, allocatable, public :: fixed(:)     !< True if a giving particle is to be kept fixed during a propagation. The default is to let the particles move.
    type(propagator_data_t),public :: prop_data
  contains
    procedure :: do_td_operation => classical_particles_do_td
    procedure :: is_tolerance_reached => classical_particles_is_tolerance_reached
    procedure :: copy_quantities_to_interaction => classical_particles_copy_quantities_to_interaction
    procedure :: update_interactions_start => classical_particles_update_interactions_start
    procedure :: update_interactions_finish => classical_particles_update_interactions_finish
    procedure :: restart_write_data => classical_particles_restart_write_data
    procedure :: restart_read_data => classical_particles_restart_read_data
    procedure :: update_kinetic_energy => classical_particles_update_kinetic_energy
    procedure :: center_of_mass => classical_particles_center_of_mass
    procedure :: center_of_mass_vel => classical_particles_center_of_mass_vel
    procedure :: center => classical_particles_center
    procedure :: axis_large => classical_particles_axis_large
    procedure :: axis_inertia => classical_particles_axis_inertia
  end type classical_particles_t


contains

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine classical_particles_init(this, np)
    class(classical_particles_t), intent(inout) :: this
    integer,                      intent(in)    :: np !< Number of particles

    PUSH_SUB(classical_particles_init)

    this%np = np
    SAFE_ALLOCATE(this%mass(1:np))
    SAFE_ALLOCATE(this%pos(1:this%space%dim, 1:np))
    SAFE_ALLOCATE(this%vel(1:this%space%dim, 1:np))
    SAFE_ALLOCATE(this%tot_force(1:this%space%dim, 1:np))
    SAFE_ALLOCATE(this%fixed(1:np))
    this%fixed = .false. ! By default we let the particles move.

    this%quantities(POSITION)%required = .true.
    this%quantities(VELOCITY)%required = .true.
    this%quantities(POSITION)%protected = .true.
    this%quantities(VELOCITY)%protected = .true.

    POP_SUB(classical_particles_init)
  end subroutine classical_particles_init

  ! ---------------------------------------------------------
  subroutine classical_particles_copy(this, cp_in)
    class(classical_particles_t), intent(out) :: this
    class(classical_particles_t), intent(in)  :: cp_in

    PUSH_SUB(classical_particles_copy)

    this%np = cp_in%np
    SAFE_ALLOCATE_SOURCE_A(this%mass,      cp_in%mass)
    SAFE_ALLOCATE_SOURCE_A(this%pos,       cp_in%pos)
    SAFE_ALLOCATE_SOURCE_A(this%vel,       cp_in%vel)
    SAFE_ALLOCATE_SOURCE_A(this%tot_force, cp_in%tot_force)
    SAFE_ALLOCATE_SOURCE_A(this%fixed,     cp_in%fixed)

    this%kinetic_energy    = cp_in%kinetic_energy

    this%quantities(POSITION)%required = .true.
    this%quantities(VELOCITY)%required = .true.
    this%quantities(POSITION)%protected = .true.
    this%quantities(VELOCITY)%protected = .true.

    this%prop_data = cp_in%prop_data

    POP_SUB(classical_particles_copy)
  end subroutine classical_particles_copy

  ! ---------------------------------------------------------
  subroutine classical_particles_init_interaction(this, interaction)
    class(classical_particles_t), target, intent(inout) :: this
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(classical_particles_init_interaction)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by classical particles."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(classical_particles_init_interaction)
  end subroutine classical_particles_init_interaction

  ! ---------------------------------------------------------
  subroutine classical_particles_do_td(this, operation)
    class(classical_particles_t),    intent(inout) :: this
    class(algorithmic_operation_t),  intent(in)    :: operation

    integer :: ii, ip
    FLOAT, allocatable :: tmp_pos(:,:,:), tmp_vel(:,:,:)
    FLOAT :: factor

    PUSH_SUB(classical_particles_do_td)

    select case (operation%id)
    case (SKIP)
      ! Do nothing
    case (STORE_CURRENT_STATUS)
      this%prop_data%save_pos(:, 1:this%np) = this%pos(:, 1:this%np)
      this%prop_data%save_vel(:, 1:this%np) = this%vel(:, 1:this%np)

    case (VERLET_START)
      if (.not. this%prop_data%initialized) then
        call this%prop_data%initialize(this%prop, this%space%dim, this%np)
        do ip = 1, this%np
          if (this%fixed(ip)) then
            this%prop_data%acc(:, ip) = M_ZERO
          else
            this%prop_data%acc(:, ip) = this%tot_force(:, ip) / this%mass(ip)
          end if
        end do
      end if

    case (VERLET_FINISH)
      call this%prop_data%end()

    case (BEEMAN_FINISH)
      call this%prop_data%end()

    case (VERLET_UPDATE_POS)
      this%pos(:, 1:this%np) = this%pos(:, 1:this%np) + this%prop%dt * this%vel(:, 1:this%np) &
        + M_HALF * this%prop%dt**2 * this%prop_data%acc(:, 1:this%np)

      this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(POSITION), &
        this%quantities(POSITION)%clock, "tick"))


    case (VERLET_COMPUTE_ACC, BEEMAN_COMPUTE_ACC)
      do ii = size(this%prop_data%prev_acc, dim=3) - 1, 1, -1
        this%prop_data%prev_acc(:, 1:this%np, ii + 1) = this%prop_data%prev_acc(:, 1:this%np, ii)
      end do
      do ip = 1, this%np
        this%prop_data%prev_acc(:, ip, 1) = this%prop_data%acc(:, ip)
        if (this%fixed(ip)) then
          this%prop_data%acc(:, ip) = M_ZERO
        else
          this%prop_data%acc(:, ip) = this%tot_force(:, ip) / this%mass(ip)
        end if
      end do

    case (VERLET_COMPUTE_VEL)
      this%vel(:, 1:this%np) = this%vel(:, 1:this%np) &
        + M_HALF * this%prop%dt * (this%prop_data%prev_acc(:, 1:this%np, 1) + this%prop_data%acc(:, 1:this%np))

      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", &
        QUANTITY_LABEL(VELOCITY), this%quantities(VELOCITY)%clock, "tick"))

    case (BEEMAN_START)
      if (.not. this%prop_data%initialized) then
        call this%prop_data%initialize(this%prop, this%space%dim, this%np)
        do ip = 1, this%np
          if (this%fixed(ip)) then
            this%prop_data%acc(:, ip) = M_ZERO
          else
            this%prop_data%acc(:, ip) = this%tot_force(:, ip) / this%mass(ip)
          end if
          this%prop_data%prev_acc(:, ip, 1) = this%prop_data%acc(:, ip)
        end do
      end if

    case (BEEMAN_PREDICT_POS)
      this%pos(:, 1:this%np) = this%pos(:, 1:this%np) + this%prop%dt * this%vel(:, 1:this%np) + &
        M_ONE/CNST(6.0) * this%prop%dt**2 * (M_FOUR*this%prop_data%acc(:, 1:this%np) - this%prop_data%prev_acc(:, 1:this%np, 1))

      if (.not. this%prop%predictor_corrector) then
        this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(POSITION), &
          this%quantities(POSITION)%clock, "tick"))
      end if

    case (BEEMAN_PREDICT_VEL)
      this%vel(:, 1:this%np) = this%vel(:, 1:this%np) + M_ONE/CNST(6.0) * this%prop%dt * ( M_TWO * this%prop_data%acc(:, 1:this%np) + &
        CNST(5.0) * this%prop_data%prev_acc(:, 1:this%np, 1) - this%prop_data%prev_acc(:, 1:this%np, 2))

      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", &
        QUANTITY_LABEL(VELOCITY), this%quantities(VELOCITY)%clock, "tick"))

    case (BEEMAN_CORRECT_POS)
      this%pos(:, 1:this%np) = this%prop_data%save_pos(:, 1:this%np) + this%prop%dt * this%prop_data%save_vel(:, 1:this%np) &
        + M_ONE/CNST(6.0) * this%prop%dt**2 * (this%prop_data%acc(:, 1:this%np) + M_TWO * this%prop_data%prev_acc(:, 1:this%np, 1))

      ! We set it to the propagation time to avoid double increment
      call this%quantities(POSITION)%clock%set_time(this%prop%clock)

    case (BEEMAN_CORRECT_VEL)
      this%vel(:, 1:this%np) = this%prop_data%save_vel(:, 1:this%np) &
        + M_HALF * this%prop%dt * (this%prop_data%acc(:, 1:this%np) + this%prop_data%prev_acc(:, 1:this%np, 1))

      ! We set it to the propagation time to avoid double increment
      call this%quantities(VELOCITY)%clock%set_time(this%prop%clock)
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", &
        QUANTITY_LABEL(VELOCITY), this%quantities(VELOCITY)%clock, "set"))

    case (EXPMID_START)
      if (.not. this%prop_data%initialized) then
        call this%prop_data%initialize(this%prop, this%space%dim, this%np)
        this%prop_data%prev_pos(:, 1:this%np, 1) = this%pos(:, 1:this%np)
        this%prop_data%prev_vel(:, 1:this%np, 1) = this%vel(:, 1:this%np)
      end if

    case (EXPMID_FINISH)
      call this%prop_data%end()

    case (EXPMID_PREDICT_DT_2)
      this%pos(:, 1:this%np) = CNST(1.5)*this%prop_data%save_pos(:, 1:this%np) - M_HALF*this%prop_data%prev_pos(:, 1:this%np, 1)
      this%vel(:, 1:this%np) = CNST(1.5)*this%prop_data%save_vel(:, 1:this%np) - M_HALF*this%prop_data%prev_vel(:, 1:this%np, 1)
      this%prop_data%prev_pos(:, 1:this%np, 1) = this%prop_data%save_pos(:, 1:this%np)
      this%prop_data%prev_vel(:, 1:this%np, 1) = this%prop_data%save_vel(:, 1:this%np)
      this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(POSITION), &
        this%quantities(POSITION)%clock, "tick"))
      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(VELOCITY), &
        this%quantities(VELOCITY)%clock, "tick"))

    case (UPDATE_HAMILTONIAN)
      do ip = 1, this%np
        if (this%fixed(ip)) then
          this%prop_data%hamiltonian_elements(:, ip) = M_ZERO
        else
          this%prop_data%hamiltonian_elements(:, ip) = this%tot_force(:, ip) / (this%mass(ip) * this%pos(:, ip))
        end if
      end do

    case (EXPMID_PREDICT_DT)
      SAFE_ALLOCATE(tmp_pos(1:this%space%dim, 1:this%np, 1:2))
      SAFE_ALLOCATE(tmp_vel(1:this%space%dim, 1:this%np, 1:2))
      ! apply exponential - at some point this could use the machinery of
      !   exponential_apply (but this would require a lot of boilerplate code
      !   like a Hamiltonian class etc)
      ! prop_data%save_pos/vel contain the state at t - this is the state we want to
      !   apply the Hamiltonian to
      tmp_pos(:, 1:this%np, 1) = this%prop_data%save_pos(:, 1:this%np)
      tmp_vel(:, 1:this%np, 1) = this%prop_data%save_vel(:, 1:this%np)
      this%pos(:, 1:this%np) = this%prop_data%save_pos(:, 1:this%np)
      this%vel(:, 1:this%np) = this%prop_data%save_vel(:, 1:this%np)
      ! compute exponential with Taylor expansion
      factor = M_ONE
      do ii = 1, 4
        factor = factor * this%prop%dt / ii
        do ip = 1, this%np
          ! apply hamiltonian
          tmp_pos(:, ip, 2) = tmp_vel(:, ip, 1)
          tmp_vel(:, ip, 2) = this%prop_data%hamiltonian_elements(:, ip) * tmp_pos(:, ip, 1)
          ! swap temporary variables
          tmp_pos(:, ip, 1) = tmp_pos(:, ip, 2)
          tmp_vel(:, ip, 1) = tmp_vel(:, ip, 2)
          ! accumulate components of Taylor expansion
          this%pos(:, ip) = this%pos(:, ip) + factor * tmp_pos(:, ip, 1)
          this%vel(:, ip) = this%vel(:, ip) + factor * tmp_vel(:, ip, 1)
        end do
      end do
      SAFE_DEALLOCATE_A(tmp_pos)
      SAFE_DEALLOCATE_A(tmp_vel)
      this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(POSITION), &
        this%quantities(POSITION)%clock, "tick"))
      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(VELOCITY), &
        this%quantities(VELOCITY)%clock, "tick"))

    case (EXPMID_CORRECT_DT_2)
      ! only correct for dt/2 if not converged yet
      if (.not. this%is_tolerance_reached(this%prop%scf_tol)) then
        this%pos(:, 1:this%np) = M_HALF*(this%pos(:, 1:this%np) + this%prop_data%save_pos(:, 1:this%np))
        this%vel(:, 1:this%np) = M_HALF*(this%vel(:, 1:this%np) + this%prop_data%save_vel(:, 1:this%np))
        this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", &
          QUANTITY_LABEL(POSITION), this%quantities(POSITION)%clock, action="tick"))
        this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(VELOCITY), &
          this%quantities(VELOCITY)%clock, "tick"))
      end if

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(classical_particles_do_td)
  end subroutine classical_particles_do_td

  ! ---------------------------------------------------------
  logical function classical_particles_is_tolerance_reached(this, tol) result(converged)
    class(classical_particles_t), intent(in)    :: this
    FLOAT,                        intent(in)    :: tol

    integer :: ip
    FLOAT :: change, max_change

    PUSH_SUB(classical_particles_is_tolerance_reached)

    ASSERT(this%prop%predictor_corrector)

    ! Here we put the criterion that maximum acceleration change is below the tolerance
    max_change = M_ZERO
    do ip = 1, this%np
      change = sum((this%prop_data%prev_tot_force(1:this%space%dim, ip) - this%tot_force(1:this%space%dim, ip))**2)/this%mass(ip)
      if (change > max_change) then
        max_change = change
      end if
    end do
    converged = max_change < tol**2

    if (debug%info) then
      write(message(1), '(a, e13.6, a, e13.6)') "Debug: -- Maximum change in acceleration  ", &
        sqrt(max_change), " and tolerance ", tol
      call messages_info(1, namespace=this%namespace)
    end if

    POP_SUB(classical_particles_is_tolerance_reached)
  end function classical_particles_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine classical_particles_update_quantity(this, iq)
    class(classical_particles_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particles_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The classical particles have a mass, but it is not necessary to update them, as they do not change with time.
      ! We still need to set their clock, so we set it to be in sync with the particles position.
      call this%quantities(iq)%clock%set_time(this%quantities(POSITION)%clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(classical_particles_update_quantity)
  end subroutine classical_particles_update_quantity

  ! ---------------------------------------------------------
  subroutine classical_particles_update_exposed_quantity(partner, iq)
    class(classical_particles_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particles_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The classical particles have a mass, but they do not require any update, as they do not change with time.
      ! We still need to set their clock, so we set it to be in sync with the particles position.
      call partner%quantities(iq)%clock%set_time(partner%quantities(POSITION)%clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(classical_particles_update_exposed_quantity)
  end subroutine classical_particles_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine classical_particles_init_interaction_as_partner(partner, interaction)
    class(classical_particles_t), intent(in)    :: partner
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(classical_particles_init_interaction_as_partner)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(classical_particles_init_interaction_as_partner)
  end subroutine classical_particles_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine classical_particles_copy_quantities_to_interaction(partner, interaction)
    class(classical_particles_t),         intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(classical_particles_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(classical_particles_copy_quantities_to_interaction)
  end subroutine classical_particles_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine classical_particles_update_interactions_start(this)
    class(classical_particles_t), intent(inout) :: this

    PUSH_SUB(classical_particles_update_interactions_start)

    ! Store previous force, as it is used as SCF criterium
    if (this%prop%predictor_corrector) then
      this%prop_data%prev_tot_force(1:this%space%dim, 1:this%np) = this%tot_force(1:this%space%dim, 1:this%np)
    end if

    POP_SUB(classical_particles_update_interactions_start)
  end subroutine classical_particles_update_interactions_start

  ! ---------------------------------------------------------
  subroutine classical_particles_update_interactions_finish(this)
    class(classical_particles_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(classical_particles_update_interactions_finish)

    ! Compute the total force acting on the classical particles
    this%tot_force(1:this%space%dim, 1:this%np) = M_ZERO
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (force_interaction_t)
        this%tot_force(1:this%space%dim, 1:this%np) = this%tot_force(1:this%space%dim, 1:this%np) + &
          interaction%force(1:this%space%dim, 1:this%np)
      end select
    end do

    POP_SUB(classical_particles_update_interactions_finish)
  end subroutine classical_particles_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine classical_particles_restart_write_data(this)
    class(classical_particles_t), intent(inout) :: this

    integer :: restart_file_unit

    PUSH_SUB(classical_particles_restart_write_data)

    call io_mkdir('restart/'//TD_DIR, this%namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_classical_particles', this%namespace, action='write')
    write(restart_file_unit, *) this%np
    write(restart_file_unit, *) this%mass(:)
    write(restart_file_unit, *) this%pos(:,:)
    write(restart_file_unit, *) this%vel(:,:)
    write(restart_file_unit, *) this%tot_force(:,:)
    call io_close(restart_file_unit)

    if (this%clock%get_tick() > 0) then
      ! only initialized after the first time step
      call this%prop_data%restart_write(this%namespace, this%prop)
    end if

    message(1) = "Successfully wrote restart data for system "//trim(this%namespace%get())
    call messages_info(1, namespace=this%namespace)

    POP_SUB(classical_particles_restart_write_data)
  end subroutine classical_particles_restart_write_data

  ! ---------------------------------------------------------
  logical function classical_particles_restart_read_data(this)
    class(classical_particles_t), intent(inout) :: this

    integer :: restart_file_unit

    PUSH_SUB(classical_particles_restart_read_data)

    call io_mkdir('restart/'//TD_DIR, this%namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_classical_particles', this%namespace, action='read', die=.false.)
    if (restart_file_unit > 0) then
      read(restart_file_unit, *) this%np
      read(restart_file_unit, *) this%mass(:)
      read(restart_file_unit, *) this%pos(:,:)
      read(restart_file_unit, *) this%vel(:,:)
      read(restart_file_unit, *) this%tot_force(:,:)
      call io_close(restart_file_unit)
      call this%prop_data%initialize(this%prop, this%space%dim, this%np)
      if (this%clock%get_tick() > 0) then
        ! only initialized after the first time step
        classical_particles_restart_read_data = this%prop_data%restart_read(this%namespace, this%prop)
      else
        classical_particles_restart_read_data = .true.
      end if
    else
      ! could not open file
      classical_particles_restart_read_data = .false.
    end if

    if (classical_particles_restart_read_data) then
      message(1) = "Successfully read restart data for system "//trim(this%namespace%get())
      ! namespace should be added here at some point
      call messages_info(1)
    end if

    POP_SUB(classical_particles_restart_read_data)
  end function classical_particles_restart_read_data

  ! ---------------------------------------------------------
  subroutine classical_particles_update_kinetic_energy(this)
    class(classical_particles_t),     intent(inout) :: this

    integer :: ip

    PUSH_SUB(classical_particles_update_kinetic_energy)

    this%kinetic_energy = M_ZERO
    do ip = 1, this%np
      this%kinetic_energy = this%kinetic_energy + M_HALF * this%mass(ip) * sum(this%vel(:, ip)**2)
    end do

    POP_SUB(classical_particles_update_kinetic_energy)
  end subroutine classical_particles_update_kinetic_energy

  ! ---------------------------------------------------------
  function classical_particles_center_of_mass(this, mask, pseudo) result(pos)
    class(classical_particles_t), intent(in) :: this
    logical,            optional, intent(in) :: mask(:)
    logical,            optional, intent(in) :: pseudo !< calculate center considering all particles to have equal mass.
    FLOAT :: pos(this%space%dim)

    FLOAT :: mass, total_mass
    integer :: ip

    PUSH_SUB(classical_particles_center_of_mass)

    ASSERT(.not. this%space%is_periodic())

    pos = M_ZERO
    total_mass = M_ZERO
    mass = M_ONE
    do ip = 1, this%np
      if (present(mask)) then
        if (.not. mask(ip)) cycle
      end if
      if (.not. optional_default(pseudo, .false.)) then
        mass = this%mass(ip)
      end if
      pos = pos + mass*this%pos(:, ip)
      total_mass = total_mass + mass
    end do
    pos = pos/total_mass

    POP_SUB(classical_particles_center_of_mass)
  end function classical_particles_center_of_mass

  ! ---------------------------------------------------------
  function classical_particles_center_of_mass_vel(this) result(vel)
    class(classical_particles_t), intent(in) :: this
    FLOAT :: vel(this%space%dim)

    FLOAT :: mass, total_mass
    integer :: ip

    PUSH_SUB(classical_particles_center_of_mass_vel)

    vel = M_ZERO
    total_mass = M_ZERO
    do ip = 1, this%np
      mass = this%mass(ip)
      total_mass = total_mass + mass
      vel = vel + mass*this%vel(:, ip)
    end do
    vel = vel/total_mass

    POP_SUB(classical_particles_center_of_mass_vel)
  end function classical_particles_center_of_mass_vel

  ! ---------------------------------------------------------
  function classical_particles_center(this) result(pos)
    class(classical_particles_t), intent(in) :: this
    FLOAT :: pos(this%space%dim)

    FLOAT :: xmin(this%space%dim), xmax(this%space%dim)
    integer  :: ip, idir

    PUSH_SUB(classical_particles_center)

    xmin =  M_HUGE
    xmax = -M_HUGE
    do ip = 1, this%np
      do idir = 1, this%space%dim
        if (this%pos(idir, ip) > xmax(idir)) xmax(idir) = this%pos(idir, ip)
        if (this%pos(idir, ip) < xmin(idir)) xmin(idir) = this%pos(idir, ip)
      end do
    end do

    pos = (xmax + xmin)/M_TWO

    POP_SUB(classical_particles_center)
  end function classical_particles_center

  ! ---------------------------------------------------------
  subroutine classical_particles_axis_large(this, x, x2)
    class(classical_particles_t), intent(in)  :: this
    FLOAT,                        intent(out) :: x(this%space%dim)
    FLOAT,                        intent(out) :: x2(this%space%dim)

    integer  :: ip, jp
    FLOAT :: rmax, r, r2

    PUSH_SUB(classical_particles_axis_large)

    ASSERT(.not. this%space%is_periodic())

    ! first get the further apart atoms
    rmax = -M_HUGE
    do ip = 1, this%np
      do jp = 1, this%np/2 + 1
        r = norm2(this%pos(:, ip) - this%pos(:, jp))
        if (r > rmax) then
          rmax = r
          x = this%pos(:, ip) - this%pos(:, jp)
        end if
      end do
    end do
    x  = x /norm2(x)

    ! now let us find out what is the second most important axis
    rmax = -M_HUGE
    do ip = 1, this%np
      r2 = sum(x * this%pos(:, ip))
      r = norm2(this%pos(:, ip) - r2*x)
      if (r > rmax) then
        rmax = r
        x2 = this%pos(:, ip) - r2*x
      end if
    end do

    POP_SUB(classical_particles_axis_large)
  end subroutine classical_particles_axis_large

  ! ---------------------------------------------------------
  !> This subroutine assumes that the origin of the coordinates is the
  !! center of mass of the system
  subroutine classical_particles_axis_inertia(this, x, x2, pseudo)
    class(classical_particles_t), intent(in)  :: this
    FLOAT,                        intent(out) :: x(this%space%dim)
    FLOAT,                        intent(out) :: x2(this%space%dim)
    logical,                      intent(in)  :: pseudo !< calculate axis considering all particles to have equal mass.

    FLOAT :: mass, tinertia(this%space%dim, this%space%dim), eigenvalues(this%space%dim)
    integer :: ii, jj, ip
    type(unit_t) :: unit

    PUSH_SUB(classical_particles_axis_inertia)

    ASSERT(.not. this%space%is_periodic())

    ! first calculate the inertia tensor
    tinertia = M_ZERO
    mass = M_ONE
    do ip = 1, this%np
      if (.not. pseudo) mass = this%mass(ip)
      do ii = 1, this%space%dim
        tinertia(ii, :) = tinertia(ii, :) - mass*this%pos(ii, ip)*this%pos(:, ip)
        tinertia(ii, ii) = tinertia(ii, ii) + mass*sum(this%pos(:, ip)**2)
      end do
    end do

    unit = units_out%length**2
    ! note: we always use amu for atomic masses, so no unit conversion to/from atomic is needed.
    if (pseudo) then
      write(message(1),'(a)') 'Moment of pseudo-inertia tensor [' // trim(units_abbrev(unit)) // ']'
    else
      write(message(1),'(a)') 'Moment of inertia tensor [amu*' // trim(units_abbrev(unit)) // ']'
    end if
    call messages_info(1, namespace=this%namespace)
    call output_tensor(tinertia, this%space%dim, unit, write_average = .true., namespace=this%namespace)

    call lalg_eigensolve(this%space%dim, tinertia, eigenvalues)

    write(message(1),'(a,6f25.6)') 'Eigenvalues: ', units_from_atomic(unit, eigenvalues)
    call messages_info(1, namespace=this%namespace)

    ! make a choice to fix the sign of the axis.
    do ii = 1, 2
      jj = maxloc(abs(tinertia(:,ii)), dim = 1)
      if (tinertia(jj,ii) < M_ZERO) tinertia(:,ii) = -tinertia(:,ii)
    end do
    x  = tinertia(:,1)
    x2 = tinertia(:,2)

    POP_SUB(classical_particles_axis_inertia)
  end subroutine classical_particles_axis_inertia

  ! ---------------------------------------------------------
  subroutine classical_particles_end(this)
    class(classical_particles_t), intent(inout) :: this

    PUSH_SUB(classical_particles_end)

    SAFE_DEALLOCATE_A(this%mass)
    SAFE_DEALLOCATE_A(this%pos)
    SAFE_DEALLOCATE_A(this%vel)
    SAFE_DEALLOCATE_A(this%tot_force)
    SAFE_DEALLOCATE_A(this%fixed)

    call system_end(this)

    POP_SUB(classical_particles_end)
  end subroutine classical_particles_end

end module classical_particles_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
