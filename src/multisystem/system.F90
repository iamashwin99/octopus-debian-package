!! Copyright (C) 2019 N. Tancogne-Dejean
!! Copyright (C) 2020 M. Oliveira, Heiko Appel
!! Copyright (C) 2021 S. Ohlmann
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

!> This module implements the abstract system type.
!!
module system_oct_m
  use algorithm_oct_m
  use algorithm_factory_oct_m
  use clock_oct_m
  use debug_oct_m
  use ghost_interaction_oct_m
  use global_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use interaction_with_partner_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use multisystem_debug_oct_m
  use linked_list_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  implicit none

  private
  public ::                                    &
    system_t,                                  &
    system_execute_algorithm,                  &
    system_init_parallelization,               &
    system_init_algorithm,                     &
    system_init_clocks,                        &
    system_reset_clocks,                       &
    system_propagation_start,                  &
    system_propagation_finish,                 &
    system_restart_read,                       &
    system_restart_write,                      &
    system_update_potential_energy,            &
    system_update_total_energy,                &
    system_end,                                &
    system_list_t,                             &
    system_iterator_t

  type :: barrier_t
    logical :: active
    FLOAT :: target_time
  end type barrier_t

  integer, parameter, public :: &
    NUMBER_BARRIERS = 1,        &
    BARRIER_RESTART = 1

  !> @brief Abstract class for systems
  !!
  !! All explicit systems are derived from this class.
  type, extends(interaction_partner_t), abstract :: system_t
    private
    class(algorithm_t),  pointer, public :: algo => null()

    integer, public :: accumulated_loop_ticks

    integer, public :: interaction_timing  !< parameter to determine if interactions
    !!                                        should use the quantities at the exact time or if retardation is allowed

    type(integer_list_t), public :: supported_interactions
    type(interaction_list_t), public :: interactions !< List with all the interactions of this system

    type(mpi_grp_t), public :: grp  !< mpi group for this system

    type(barrier_t) :: barrier(NUMBER_BARRIERS)
    FLOAT, public :: kinetic_energy     !< Energy not from interactions, like the kinetic energy
    FLOAT, public :: potential_energy   !< Energy from the interactions with external systems
    FLOAT, public :: internal_energy    !< Energy from the interactions with itself and for containers the kinetic energy of its constituents
    FLOAT, public :: total_energy       !< Sum of internal, external, and self energy

  contains
    procedure :: execute_algorithm =>  system_execute_algorithm !< @copydoc system_oct_m::system_execute_algorithm
    procedure :: reset_clocks => system_reset_clocks !< @copydoc system_oct_m::system_reset_clocks
    procedure :: update_exposed_quantities => system_update_exposed_quantities !< @copydoc system_oct_m::system_update_exposed_quantities
    procedure :: init_algorithm => system_init_algorithm !< @copydoc system_oct_m::system_init_algorithm
    procedure :: algorithm_finished => system_algorithm_finished !< @copydoc system_oct_m::system_algorithm_finished
    procedure :: init_clocks => system_init_clocks !< @copydoc system_oct_m::system_init_clocks
    procedure :: init_all_interactions => system_init_all_interactions !< @copydoc system_oct_m::system_init_all_interactions
    procedure :: init_parallelization => system_init_parallelization !< @copydoc system_oct_m::system_init_parallelization
    procedure :: update_interactions => system_update_interactions !< @copydoc system_oct_m::system_update_interactions
    procedure :: update_interactions_start => system_update_interactions_start !< @copydoc system_oct_m::system_update_interactions_start
    procedure :: update_interactions_finish => system_update_interactions_finish !< @copydoc system_oct_m::system_update_interactions_finish
    procedure :: propagation_start => system_propagation_start !< @copydoc system_oct_m::system_propagation_start
    procedure :: propagation_finish => system_propagation_finish !< @copydoc system_oct_m::system_propagation_finish
    procedure :: iteration_info => system_iteration_info !< @copydoc system_oct_m::system_iteration_info
    procedure :: restart_write => system_restart_write !< @copydoc system_oct_m::system_restart_write
    procedure :: restart_read => system_restart_read !< @copydoc system_oct_m::system_restart_read
    procedure :: output_start => system_output_start !< @copydoc system_oct_m::system_output_start
    procedure :: output_write => system_output_write !< @copydoc system_oct_m::system_output_write
    procedure :: output_finish => system_output_finish !< @copydoc system_oct_m::system_output_finish
    procedure :: process_is_slave => system_process_is_slave !< @copydoc system_oct_m::system_process_is_slave
    procedure :: start_barrier => system_start_barrier !< @copydoc system_oct_m::system_start_barrier
    procedure :: end_barrier => system_end_barrier !< @copydoc system_oct_m::system_end_barrier
    procedure :: arrived_at_barrier => system_arrived_at_barrier !< @copydoc system_oct_m::system_arrived_at_barrier
    procedure :: arrived_at_any_barrier => system_arrived_at_any_barrier !< @copydoc system_oct_m::system_arrived_at_any_barrier
    procedure :: update_potential_energy => system_update_potential_energy !< @copydoc system_oct_m::system_update_potential_energy
    procedure :: update_internal_energy => system_update_internal_energy !< @copydoc system_oct_m::system_update_internal_energy
    procedure :: update_total_energy => system_update_total_energy !< @copydoc system_oct_m::system_update_total_energy
    procedure(system_init_interaction),          deferred :: init_interaction !< @copydoc system_oct_m::system_init_interaction
    procedure(system_initial_conditions),        deferred :: initial_conditions !< @copydoc system_oct_m::system_initial_conditions
    procedure(system_do_algorithmic_operation),  deferred :: do_algorithmic_operation !< @copydoc system_oct_m::system_do_algorithmic_operation
    procedure(system_is_tolerance_reached),      deferred :: is_tolerance_reached !< @copydoc system_oct_m::system_is_tolerance_reached
    procedure(system_update_quantity),           deferred :: update_quantity !< @copydoc system_oct_m::system_update_quantity
    procedure(system_update_exposed_quantity),   deferred :: update_exposed_quantity !< @copydoc system_oct_m::system_update_exposed_quantity
    procedure(system_restart_write_data),        deferred :: restart_write_data !< @copydoc system_oct_m::system_restart_write_data
    procedure(system_restart_read_data),         deferred :: restart_read_data !< @copydoc system_oct_m::system_restart_read_data
    procedure(system_update_kinetic_energy),     deferred :: update_kinetic_energy !< @copydoc system_oct_m::system_update_kinetic_energy
  end type system_t

  abstract interface

    ! ---------------------------------------------------------
    !> @brief initialize a given interaction of the system
    subroutine system_init_interaction(this, interaction)
      import system_t
      import interaction_t
      class(system_t), target, intent(inout) :: this
      class(interaction_t),    intent(inout) :: interaction
    end subroutine system_init_interaction

    ! ---------------------------------------------------------
    !> set initial conditions for a system
    subroutine system_initial_conditions(this)
      import system_t
      class(system_t), intent(inout) :: this
    end subroutine system_initial_conditions

    ! ---------------------------------------------------------
    !> Execute one operation that is part of a larger algorithm.  Returns true
    !! if the operation was successfully executed, false otherwise.
    logical function system_do_algorithmic_operation(this, operation) result(done)
      import system_t
      import algorithmic_operation_t
      class(system_t),                intent(inout) :: this
      class(algorithmic_operation_t), intent(in)    :: operation
    end function system_do_algorithmic_operation

    ! ---------------------------------------------------------
    !> @brief check whether a system has reached a given tolerance
    logical function system_is_tolerance_reached(this, tol)
      import system_t
      class(system_t), intent(in) :: this
      FLOAT,           intent(in) :: tol
    end function system_is_tolerance_reached

    ! ---------------------------------------------------------
    subroutine system_store_current_status(this)
      import system_t
      class(system_t), intent(inout) :: this
    end subroutine system_store_current_status

    ! ---------------------------------------------------------
    subroutine system_update_quantity(this, iq)
      import system_t
      import clock_t
      class(system_t),      intent(inout) :: this
      integer,              intent(in)    :: iq
    end subroutine system_update_quantity

    ! ---------------------------------------------------------
    subroutine system_update_exposed_quantity(partner, iq)
      import system_t
      import clock_t
      class(system_t), intent(inout) :: partner
      integer,         intent(in)    :: iq
    end subroutine system_update_exposed_quantity

    ! ---------------------------------------------------------
    subroutine system_restart_write_data(this)
      import system_t
      class(system_t), intent(inout) :: this
    end subroutine system_restart_write_data

    ! ---------------------------------------------------------
    ! this function returns true if restart data could be read
    logical function system_restart_read_data(this)
      import system_t
      class(system_t), intent(inout) :: this
    end function system_restart_read_data
    subroutine system_update_kinetic_energy(this)
      import system_t
      class(system_t),      intent(inout) :: this
    end subroutine system_update_kinetic_energy

  end interface

  !> @brief These classes extends the list and list iterator to create a system list.
  !!
  !! Since a list of systems is also a list of interaction partners, the system
  !! list is an extension of the partner list.
  type, extends(partner_list_t) :: system_list_t
    private
  contains
    procedure :: add => system_list_add_node !< @copydoc system_oct_m::system_list_add_node
    procedure :: contains => system_list_contains !< @copydoc system_oct_m::system_list_contains
  end type system_list_t

  type, extends(linked_list_iterator_t) :: system_iterator_t
    private
  contains
    procedure :: get_next => system_iterator_get_next !< @copydoc system_oct_m::system_iterator_get_next
  end type system_iterator_t

contains

  ! ---------------------------------------------------------
  !> @brief perform one or more algorithmic operations
  !!
  !! The following subroutine takes a system and performs as many algorithmic
  !! operations as possible on the system until a barrier is reached. There are
  !! two types of barriers:
  !!
  !!  - explicit barriers, implemented using the barrier_t type
  !!  - the interaction update
  !!
  !! The interaction update is always considered a barrier, even if the update
  !! was successful. This is to allow other system to also update their
  !! interactions before this system moves on to the next operations.
  subroutine system_execute_algorithm(this)
    class(system_t),     intent(inout) :: this

    type(algorithmic_operation_t) :: operation
    logical :: all_updated, at_barrier, operation_done
    type(event_handle_t) :: debug_handle

    PUSH_SUB(system_execute_algorithm)

    at_barrier = .false.

    do while (.not. at_barrier)

      operation = this%algo%get_current_operation()

      if (debug%info) then
        write(message(1), '(a,a,1X,a)') "Debug: Start  ", trim(operation%label), " for '" + trim(this%namespace%get()) + "'"
        call messages_info(1, namespace=this%namespace)
      end if

      debug_handle = multisystem_debug_write_event_in(this%namespace, event_function_call_t("dt_operation", operation),    &
        system_clock=this%clock, prop_clock=this%algo%clock)

      ! First try to execute the operation as a system specific operation
      operation_done = this%do_algorithmic_operation(operation)

      ! If not done, we try to execute it as an algorithm-specific operation.
      if (.not. operation_done) then
        operation_done = this%algo%do_operation(operation)
      else
        call this%algo%next()
      end if

      ! If still not done, the operation must be a generic operation
      if (.not. operation_done) then

        select case (operation%id)
        case (SKIP)
          ! Do nothing
          call this%algo%next()

        case (STEP_DONE)
          ! Increment the system clock by one time-step
          this%clock = this%clock + CLOCK_TICK
          call multisystem_debug_write_marker(this%namespace, event_clock_update_t("system",  "", this%clock, "tick"))

          ! Recompute the total energy
          call this%update_total_energy()

          ! Write output
          call this%output_write()

          ! Update elapsed time
          call this%algo%update_elapsed_time()

          ! Print information about the current iteration
          ! (NB: needs to be done after marking the propagation step as finished,
          ! so that the timings are correct)
          call this%iteration_info()

          call this%algo%next()

        case (REWIND_ALGORITHM)
          if (.not. this%arrived_at_any_barrier() .and. .not. this%algorithm_finished()) then
            ! Reset propagator for next step if not waiting at barrier and if
            ! the algorithm is not finished
            call this%algo%rewind()
          else
            at_barrier = .true.
          end if

        case (UPDATE_INTERACTIONS)
          ! We increment by one algorithmic step
          this%algo%clock = this%algo%clock + CLOCK_TICK
          call multisystem_debug_write_marker(this%namespace, event_clock_update_t("propagator", "", this%algo%clock, "tick"))

          ! Try to update all the interactions
          all_updated = this%update_interactions()

          ! Move to next algorithm step if all interactions have been
          ! updated. Otherwise try again later.
          if (all_updated) then
            this%accumulated_loop_ticks = this%accumulated_loop_ticks + 1
            call this%algo%next()
          else
            this%algo%clock = this%algo%clock - CLOCK_TICK
            call multisystem_debug_write_marker(this%namespace, event_clock_update_t("propagator", "", this%algo%clock, "reverse"))
          end if

          ! Interactions are implicit barriers
          at_barrier = .true.

        case default
          message(1) = "Unsupported algorithmic operation."
          write(message(2), '(A,A,A)') trim(operation%id), ": ", trim(operation%label)
          call messages_fatal(2, namespace=this%namespace)
        end select
      end if

      if (debug%info) then
        write(message(1), '(a,a,1X,a, l)') "Debug: Finish ", trim(operation%label), " for '" + trim(this%namespace%get()) + "' "
        call messages_info(1, namespace=this%namespace)
      end if

      call multisystem_debug_write_event_out(debug_handle, system_clock=this%clock, prop_clock=this%algo%clock)
    end do

    POP_SUB(system_execute_algorithm)
  end subroutine system_execute_algorithm

  ! ---------------------------------------------------------
  subroutine system_reset_clocks(this, accumulated_ticks)
    class(system_t),      intent(inout) :: this
    integer,              intent(in)    :: accumulated_ticks

    integer :: iq
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    character(len=MAX_INFO_LEN) :: extended_label

    PUSH_SUB(system_reset_clocks)

    ! Propagator clock
    this%algo%clock = this%algo%clock - accumulated_ticks*CLOCK_TICK
    call multisystem_debug_write_marker(this%namespace, event_clock_update_t("propagator", "", this%algo%clock, "reset"))

    ! Interaction clocks
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      interaction%clock = interaction%clock - accumulated_ticks*CLOCK_TICK

      select type (interaction)
      class is (interaction_with_partner_t)
        extended_label = trim(interaction%label)//"-"//trim(interaction%partner%namespace%get())
      class default
        extended_label = trim(interaction%label)
      end select
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t( "interaction", extended_label, &
        interaction%clock, "reset"))
    end do

    ! Internal quantities clocks
    do iq = 1, MAX_QUANTITIES
      if (this%quantities(iq)%required) then
        this%quantities(iq)%clock = this%quantities(iq)%clock - accumulated_ticks*CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", QUANTITY_LABEL(iq), &
          this%quantities(iq)%clock, "reset"))
      end if
    end do

    POP_SUB(system_reset_clocks)
  end subroutine system_reset_clocks

  ! ---------------------------------------------------------
  !> @brief update all exposed quantities of the system.
  !!
  !! this function is called as partner from the interaction
  logical function system_update_exposed_quantities(partner, requested_time, interaction) result(allowed_to_update)
    class(system_t),      intent(inout) :: partner
    type(clock_t),        intent(in)    :: requested_time
    class(interaction_t), intent(inout) :: interaction

    logical :: ahead_in_time, right_on_time, need_to_copy
    integer :: iq, q_id

    type(event_handle_t) :: debug_handle

    PUSH_SUB(system_update_exposed_quantities)

    if (debug%info) then
      write(message(1), '(a,a,a)') "Debug: -- Updating exposed quantities for partner '", trim(partner%namespace%get()), "'"
      call messages_info(1, namespace=partner%namespace)
    end if

    debug_handle = multisystem_debug_write_event_in(system_namespace = partner%namespace, &
      event = event_function_call_t("system_update_exposed_quantities"), &
      partner_clock = partner%clock, &
      requested_clock = requested_time, &
      interaction_clock = interaction%clock)

    select type (interaction)
    class is (interaction_with_partner_t)

      if (partner%algo%inside_scf .or. partner%algo%clock + CLOCK_TICK < requested_time) then
        ! we are inside an SCF cycle and therefore are not allowed to expose any quantities.
        ! or we are too much behind the requested time
        allowed_to_update = .false.
      else
        allowed_to_update = .true.
        need_to_copy = .true.
        do iq = 1, interaction%n_partner_quantities
          ! Get the requested quantity ID
          q_id = interaction%partner_quantities(iq)

          ! All needed quantities must have been marked as required. If not, then fix your code!
          ASSERT(partner%quantities(q_id)%required)

          ! First update the exposed quantities that are updated on demand
          if (partner%quantities(q_id)%updated_on_demand) then
            if (partner%quantities(q_id)%clock /= requested_time .and. &
              partner%quantities(q_id)%clock + CLOCK_TICK <= requested_time) then
              ! We can update because the partner will reach this time in the next sub-timestep
              call partner%update_exposed_quantity(q_id)

              partner%quantities(q_id)%clock = partner%quantities(q_id)%clock + CLOCK_TICK

              call updated_quantity_debug()
            else
              call not_updated_quantity_debug()
            end if
          else
            call auto_update_quantity_debug()
          end if

          if (partner%quantities(q_id)%available_at_any_time) then
            ! We ignore the clock times when a quantity is available at any time
            ahead_in_time = .false.
            right_on_time = .true.
          else
            ! Compare the clock times
            ahead_in_time = partner%quantities(q_id)%clock > requested_time
            right_on_time = partner%quantities(q_id)%clock == requested_time
          end if

          select case (partner%interaction_timing)
          case (OPTION__INTERACTIONTIMING__TIMING_EXACT)
            ! only allow interaction at exactly the same time
            allowed_to_update = allowed_to_update .and. right_on_time
            need_to_copy = allowed_to_update
          case (OPTION__INTERACTIONTIMING__TIMING_RETARDED)
            ! allow retarded interaction
            allowed_to_update = allowed_to_update .and. &
              (right_on_time .or. ahead_in_time)
            need_to_copy = need_to_copy .and. .not. ahead_in_time
          case default
            call messages_not_implemented("Method for interaction quantity timing", namespace=partner%namespace)
          end select

        end do

        ! If the quantities have been updated, we copy them to the interaction
        if (allowed_to_update .and. need_to_copy) then
          select type (interaction)
          type is (ghost_interaction_t)
            ! Nothing to copy. We still need to check that we are at the right
            ! time for the update though!
          class default
            call partner%copy_quantities_to_interaction(interaction)
          end select
        end if
      end if

    class default
      message(1) = "A system can only expose quantities to an interaction as a partner."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    if (debug%info) then
      write(message(1), '(a,a,a)') "Debug: -- Finished updating exposed quantities for partner '", &
        trim(partner%namespace%get()), "'"
      call messages_info(1, namespace=partner%namespace)
    end if

    call multisystem_debug_write_event_out(debug_handle, update=allowed_to_update, &
      partner_clock = partner%clock, &
      requested_clock = requested_time, &
      interaction_clock = interaction%clock)

    POP_SUB(system_update_exposed_quantities)
  contains

    subroutine updated_quantity_debug()

      if (debug%info) then
        write(message(1), '(a,a,a)') "Debug: ---- Updated exposed quantity ", trim(QUANTITY_LABEL(q_id)), "'"
        write(message(2), '(a,f16.6,a,f16.6)') "Debug: ------ Requested time is ", requested_time%time(), &
          ", quantity time is ", partner%quantities(q_id)%clock%time()
        call messages_info(2, namespace=partner%namespace)
      end if

    end subroutine updated_quantity_debug

    subroutine not_updated_quantity_debug()

      if (debug%info) then
        write(message(1), '(a,a,a)') "Debug: ---- Did not update exposed quantity '", trim(QUANTITY_LABEL(q_id)), "'"
        write(message(2), '(a,f16.6,a,f16.6,a,f16.6)') "Debug: ------ Requested time is ", requested_time%time(), &
          ", quantity time is ", partner%quantities(q_id)%clock%time(), &
          " and partner propagator time is ", partner%algo%clock%time()
        call messages_info(2, namespace=partner%namespace)
      end if

    end subroutine not_updated_quantity_debug

    subroutine auto_update_quantity_debug()

      if (debug%info) then
        write(message(1), '(a,a,a)') "Debug: ---- Skip update of quantity '", trim(QUANTITY_LABEL(q_id)), &
          "' as it is updated automatically"
        write(message(2), '(a,f16.6,a,f16.6)') "Debug: ------ Requested time is ", requested_time%time(), &
          ", quantity time is ", partner%quantities(q_id)%clock%time()
        call messages_info(2, namespace=partner%namespace)
      end if

    end subroutine auto_update_quantity_debug

  end function system_update_exposed_quantities

  ! ---------------------------------------------------------
  subroutine system_init_all_interactions(this)
    class(system_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(system_init_all_interactions)

    if (debug%info) then
      write(message(1), '(a,a,1X,a)') "Debug: Start  init_all_interactions for "+ trim(this%namespace%get()) + "'"
      call messages_info(1, namespace=this%namespace)
    end if

    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      select type (interaction)
      type is (ghost_interaction_t)
        ! Skip the ghost interactions
      class is (interaction_with_partner_t)
        call this%init_interaction(interaction)
        call interaction%partner%init_interaction_as_partner(interaction)
      class default
        call this%init_interaction(interaction)
      end select
    end do

    if (debug%info) then
      write(message(1), '(a,a,1X,a)') "Debug: Finish init_all_interactions for "+ trim(this%namespace%get()) + "'"
      call messages_info(1, namespace=this%namespace)
    end if

    POP_SUB(system_init_all_interactions)
  end subroutine system_init_all_interactions

  ! ---------------------------------------------------------
  logical function system_update_interactions(this) result(all_updated)
    class(system_t),      intent(inout) :: this

    logical :: none_updated
    integer :: iq, q_id
    class(interaction_t), pointer :: interaction
    type(interaction_iterator_t) :: iter

    PUSH_SUB(system_update_interactions)

    ! Some systems might need to perform some specific operations before the
    ! update. This should only be done if no interaction has been updated yet,
    ! so that it is only done once.
    none_updated = .true.
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      if (interaction%clock == this%algo%clock) then
        none_updated = .false.
        exit
      end if
    end do
    if (none_updated) then
      call this%update_interactions_start()
    end if

    !Loop over all interactions
    all_updated = .true.
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()

      if (.not. interaction%clock == this%algo%clock) then
        ! Update the system quantities that will be needed for computing the interaction
        do iq = 1, interaction%n_system_quantities
          ! Get requested quantity ID
          q_id = interaction%system_quantities(iq)

          ! All needed quantities must have been marked as required. If not, then fix your code!
          ASSERT(this%quantities(q_id)%required)

          ! We do not need to update quantities that are not updated on demand,
          ! as the algorithm takes care of doing that
          if (.not. this%quantities(q_id)%updated_on_demand) cycle

          if (.not. this%quantities(q_id)%clock == this%algo%clock) then
            ! The requested quantity is not at the requested time, so we try to update it

            ! Sanity check: it should never happen that the quantity is in advance
            ! with respect to the requested time.
            if (this%quantities(q_id)%clock > this%algo%clock) then
              message(1) = "The quantity clock is in advance compared to the requested time."
              call messages_fatal(1, namespace=this%namespace)
            end if

            call this%update_quantity(q_id)
          end if

        end do

        ! We can now try to update the interaction
        all_updated = interaction%update(this%algo%clock) .and. all_updated
      end if
    end do

    ! Some systems might need to perform some specific operations after all the
    ! interactions have been updated
    if (all_updated) then
      call this%update_interactions_finish()
    end if

    POP_SUB(system_update_interactions)
  end function system_update_interactions

  ! ---------------------------------------------------------
  subroutine system_update_interactions_start(this)
    class(system_t), intent(inout) :: this

    PUSH_SUB(system_update_interactions_start)

    ! By default nothing is done just before updating the interactions. Child
    ! classes that wish to change this behaviour should override this method.

    POP_SUB(system_update_interactions_start)
  end subroutine system_update_interactions_start

  ! ---------------------------------------------------------
  subroutine system_update_interactions_finish(this)
    class(system_t), intent(inout) :: this

    PUSH_SUB(system_update_interactions_finish)

    ! By default nothing is done just after updating the interactions. Child
    ! classes that wish to change this behaviour should override this method.

    POP_SUB(system_update_interactions_finish)
  end subroutine system_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine system_restart_write(this)
    class(system_t), intent(inout) :: this

    logical :: restart_write
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction
    integer :: ii

    PUSH_SUB(system_restart_write)

    call parse_variable(this%namespace, 'RestartWrite', .true., restart_write)

    if (restart_write) then
      ! do some generic restart steps here
      ! write clock restart data
      call this%clock%restart_write('restart_clock_system', this%namespace)
      call this%algo%clock%restart_write('restart_clock_propagator', this%namespace)
      call iter%start(this%interactions)
      do while (iter%has_next())
        interaction => iter%get_next()
        call interaction%restart_write(this%namespace)
      end do
      do ii = 1, MAX_QUANTITIES
        if (this%quantities(ii)%required) then
          call this%quantities(ii)%clock%restart_write('restart_clock_quantity_'//trim(QUANTITY_LABEL(ii)), &
            this%namespace)
        end if
      end do
      ! the following call is delegated to the corresponding system
      call this%restart_write_data()
      message(1) = "Wrote restart data for system "//trim(this%namespace%get())
      call messages_info(1, namespace=this%namespace)
    end if

    POP_SUB(system_restart_write)
  end subroutine system_restart_write

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function system_restart_read(this)
    class(system_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction
    integer :: ii

    PUSH_SUB(system_restart_read)

    ! do some generic restart steps here
    ! read clock data
    system_restart_read = this%clock%restart_read('restart_clock_system', this%namespace)
    system_restart_read = system_restart_read .and. &
      this%algo%clock%restart_read('restart_clock_propagator', this%namespace)
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      system_restart_read = system_restart_read .and. interaction%restart_read(this%namespace)
      ! reduce by one because of the first UPDATE_INTERACTIONS
      interaction%clock = interaction%clock - CLOCK_TICK
    end do
    do ii = 1, MAX_QUANTITIES
      if (this%quantities(ii)%required) then
        system_restart_read = system_restart_read .and. &
          this%quantities(ii)%clock%restart_read('restart_clock_quantity_'//trim(QUANTITY_LABEL(ii)), &
          this%namespace)
      end if
    end do
    ! the following call is delegated to the corresponding system
    system_restart_read = system_restart_read .and. this%restart_read_data()

    if (system_restart_read) then
      message(1) = "Successfully read restart data for system "//trim(this%namespace%get())
      call messages_info(1, namespace=this%namespace)
    end if

    POP_SUB(system_restart_read)
  end function system_restart_read

  ! ---------------------------------------------------------
  subroutine system_output_start(this)
    class(system_t), intent(inout) :: this

    PUSH_SUB(system_output_start)

    ! By default nothing is done to regarding output. Child classes that wish
    ! to change this behaviour should override this method.

    POP_SUB(system_output_start)
  end subroutine system_output_start

  ! ---------------------------------------------------------
  subroutine system_output_write(this)
    class(system_t), intent(inout) :: this

    PUSH_SUB(system_output_write)

    ! By default nothing is done to regarding output. Child classes that wish
    ! to change this behaviour should override this method.

    POP_SUB(system_output_write)
  end subroutine system_output_write

  ! ---------------------------------------------------------
  subroutine system_output_finish(this)
    class(system_t), intent(inout) :: this

    PUSH_SUB(system_output_finish)

    ! By default nothing is done to regarding output. Child classes that wish
    ! to change this behaviour should override this method.

    POP_SUB(system_output_finish)
  end subroutine system_output_finish

  ! ---------------------------------------------------------
  subroutine system_init_algorithm(this, factory)
    class(system_t),            intent(inout) :: this
    class(algorithm_factory_t), intent(in)    :: factory

    integer :: ii

    PUSH_SUB(system_init_algorithm)

    call messages_experimental('Multi-system framework')

    this%algo => factory%create(this)

    call this%init_clocks()

    !%Variable InteractionTiming
    !%Type integer
    !%Default timing_exact
    !%Section Time-Dependent::Propagation
    !%Description
    !% A parameter to determine if interactions should use the quantities
    !% at the exact time or if retardation is allowed.
    !%Option timing_exact 1
    !% Only allow interactions at exactly the same times
    !%Option timing_retarded 2
    !% Allow retarded interactions
    !%End
    call parse_variable(this%namespace, 'InteractionTiming', &
      OPTION__INTERACTIONTIMING__TIMING_EXACT, &
      this%interaction_timing)
    if (.not. varinfo_valid_option('InteractionTiming', this%interaction_timing)) then
      call messages_input_error(this%namespace, 'InteractionTiming')
    end if
    call messages_print_var_option('InteractionTiming', this%interaction_timing, namespace=this%namespace)

    do ii = 1, NUMBER_BARRIERS
      this%barrier(ii)%active = .false.
      this%barrier(ii)%target_time = M_ZERO
    end do

    POP_SUB(system_init_algorithm)
  end subroutine system_init_algorithm

  ! ---------------------------------------------------------------------------------------
  recursive function system_algorithm_finished(this) result(finished)
    class(system_t),       intent(in) :: this
    logical :: finished

    finished = this%algo%finished()

  end function system_algorithm_finished

  ! ---------------------------------------------------------
  subroutine system_init_clocks(this)
    class(system_t),            intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(system_init_clocks)

    ! Initialize propagator clock
    this%algo%clock = clock_t(time_step=this%algo%dt/this%algo%algo_steps)

    ! Initialize system clock
    this%clock = clock_t(time_step=this%algo%dt)

    ! Interactions clocks
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      interaction%clock = this%algo%clock - CLOCK_TICK
    end do

    ! Required quantities clocks
    where (this%quantities%required)
      this%quantities%clock = this%algo%clock
    end where

    POP_SUB(system_init_clocks)
  end subroutine system_init_clocks

  ! ---------------------------------------------------------
  subroutine system_propagation_start(this)
    class(system_t),      intent(inout) :: this

    logical :: all_updated
    type(event_handle_t) :: debug_handle

    PUSH_SUB(system_propagation_start)

    debug_handle = multisystem_debug_write_event_in(this%namespace, event_function_call_t("system_propagation_start"), &
      system_clock = this%clock, prop_clock = this%algo%clock)

    if (debug%info) then
      write(message(1), '(a,a,1X,a)') "Debug: Start  propagation_start for '" + trim(this%namespace%get()) + "'"
      call messages_info(1, namespace=this%namespace)
    end if

    ! Update interactions at initial time
    all_updated = this%update_interactions()
    if (.not. all_updated) then
      message(1) = "Unable to update interactions when initializing the propagation."
      call messages_fatal(1, namespace=this%namespace)
    end if

    ! System-specific and propagator-specific initialization step
    if (this%algo%start_step%id /= SKIP) then
      if (.not. this%do_algorithmic_operation(this%algo%start_step)) then
        message(1) = "Unsupported algorithmic operation."
        write(message(2), '(A,A,A)') trim(this%algo%start_step%id), ": ", trim(this%algo%start_step%label)
        call messages_fatal(2, namespace=this%namespace)
      end if
    end if

    ! Compute the total energy at the beginning of the simulation
    call this%update_total_energy()

    ! Start output
    call this%output_start()

    ! Write header for propagation log
    call messages_print_with_emphasis(msg="Multi-system propagation", namespace=this%namespace)
    write(message(1), '(a6,1x,a14,1x,a13,1x,a10,1x,a15)') 'Iter', 'Time', 'Energy', 'SC Steps', 'Elapsed Time'
    call messages_info(1, namespace=this%namespace)
    call messages_print_with_emphasis(namespace=this%namespace)

    ! Rewind propagator (will also set the initial time to compute the elapsed time)
    call this%algo%rewind()

    if (debug%info) then
      write(message(1), '(a,a,1X,a)') "Debug: Finish propagation_start for '" + trim(this%namespace%get()) + "'"
      call messages_info(1, namespace=this%namespace)
    end if

    call multisystem_debug_write_event_out(debug_handle, system_clock = this%clock, prop_clock = this%algo%clock)

    POP_SUB(system_propagation_start)
  end subroutine system_propagation_start

  ! ---------------------------------------------------------
  subroutine system_propagation_finish(this)
    class(system_t),      intent(inout) :: this
    type(event_handle_t) :: debug_handle

    PUSH_SUB(system_propagation_finish)

    debug_handle = multisystem_debug_write_event_in(this%namespace, event_function_call_t("system_propagation_finish"), &
      system_clock = this%clock, prop_clock = this%algo%clock)

    ! Finish output
    call this%output_finish()

    ! System-specific and propagator-specific finalization step
    if (this%algo%final_step%id /= SKIP) then
      if (.not.  this%do_algorithmic_operation(this%algo%final_step)) then
        message(1) = "Unsupported algorithmic operation."
        write(message(2), '(A,A,A)') trim(this%algo%final_step%id), ": ", trim(this%algo%final_step%label)
        call messages_fatal(2, namespace=this%namespace)
      end if
    end if

    call multisystem_debug_write_event_out(debug_handle, system_clock = this%clock, prop_clock = this%algo%clock)

    POP_SUB(system_propagation_finish)
  end subroutine system_propagation_finish

  ! ---------------------------------------------------------
  subroutine system_iteration_info(this)
    class(system_t), intent(in) :: this

    FLOAT :: energy
    character(len=40) :: fmt

    PUSH_SUB(system_iteration_info)

    energy = units_from_atomic(units_out%energy, this%total_energy)
    if (abs(energy) >= 1e5) then
      fmt = '(i7,1x,f14.6,1X,es13.6,1X,i9,1X,'
    else
      fmt = '(i7,1x,f14.6,1X,f13.6,1X,i9,1X,'
    end if
    if (this%algo%elapsed_time < 1e-3) then
      fmt = trim(fmt)//'es13.3)'
    else
      fmt = trim(fmt)//'f13.3)'
    end if

    write(message(1), fmt) this%clock%get_tick(), &
      units_from_atomic(units_out%time, this%clock%time()), energy, &
      0, this%algo%elapsed_time
    call messages_info(1, namespace=this%namespace)

    POP_SUB(system_iteration_info)
  end subroutine system_iteration_info

  ! ---------------------------------------------------------
  logical function system_process_is_slave(this)
    class(system_t), intent(in) :: this

    PUSH_SUB(system_process_is_slave)

    ! By default an MPI process is not a slave
    system_process_is_slave = .false.

    POP_SUB(system_process_is_slave)
  end function system_process_is_slave

  ! ---------------------------------------------------------
  subroutine system_end(this)
    class(system_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(system_end)

    ! No call to safe_deallocate macro here, as it gives an ICE with gfortran
    if (associated(this%algo)) then
      deallocate(this%algo)
    end if

    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      SAFE_DEALLOCATE_P(interaction)
    end do

    POP_SUB(system_end)
  end subroutine system_end

  ! ---------------------------------------------------------
  subroutine system_list_add_node(this, partner)
    class(system_list_t)         :: this
    class(interaction_partner_t), target :: partner

    PUSH_SUB(system_list_add_node)

    select type (partner)
    class is (system_t)
      call this%add_ptr(partner)
    class default
      ASSERT(.false.)
    end select

    POP_SUB(system_list_add_node)
  end subroutine system_list_add_node

  ! ---------------------------------------------------------
  recursive logical function system_list_contains(this, partner) result(contains)
    class(system_list_t)         :: this
    class(interaction_partner_t), target :: partner

    type(partner_iterator_t)  :: iterator
    class(interaction_partner_t),  pointer :: system

    PUSH_SUB(system_list_contains)

    contains = .false.

    select type (partner)
    class is (system_t)

      call iterator%start(this)
      do while (iterator%has_next() .and. .not. contains)
        system => iterator%get_next()
        contains = associated(system, partner)
      end do

    class default
      contains = .false.
    end select

    POP_SUB(system_list_contains)
  end function system_list_contains

  ! ---------------------------------------------------------
  function system_iterator_get_next(this) result(system)
    class(system_iterator_t), intent(inout) :: this
    class(system_t),          pointer       :: system

    PUSH_SUB(system_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (system_t)
      system => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(system_iterator_get_next)
  end function system_iterator_get_next

  ! ---------------------------------------------------------
  !> Basic functionality: copy the MPI group.
  !! This function needs to be implemented by extended types
  !! that need more initialization for their parallelization.
  subroutine system_init_parallelization(this, grp)
    class(system_t), intent(inout) :: this
    type(mpi_grp_t), intent(in)    :: grp

    PUSH_SUB(system_init_parallelization)

    call mpi_grp_copy(this%grp, grp)
    call messages_update_mpi_grp(this%namespace, grp)

    POP_SUB(system_init_parallelization)
  end subroutine system_init_parallelization



  ! ---------------------------------------------------------
  subroutine system_start_barrier(this, target_time, barrier_index)
    class(system_t), intent(inout) :: this
    FLOAT,           intent(in)    :: target_time
    integer,         intent(in)    :: barrier_index

    PUSH_SUB(system_start_barrier)

    this%barrier(barrier_index)%active = .true.
    this%barrier(barrier_index)%target_time = target_time

    POP_SUB(system_start_barrier)
  end subroutine system_start_barrier

  ! ---------------------------------------------------------
  subroutine system_end_barrier(this, barrier_index)
    class(system_t), intent(inout) :: this
    integer,         intent(in)    :: barrier_index

    PUSH_SUB(system_end_barrier)

    this%barrier(barrier_index)%active = .false.
    this%barrier(barrier_index)%target_time = M_ZERO

    POP_SUB(system_end_barrier)
  end subroutine system_end_barrier

  ! ---------------------------------------------------------
  logical function system_arrived_at_barrier(this, barrier_index)
    class(system_t), intent(inout) :: this
    integer,         intent(in)    :: barrier_index

    type(clock_t) :: clock

    PUSH_SUB(system_arrived_at_barrier)

    system_arrived_at_barrier = .false.
    if (this%barrier(barrier_index)%active) then
      clock = this%clock + CLOCK_TICK
      if (clock%time() > this%barrier(barrier_index)%target_time) then
        system_arrived_at_barrier = .true.
      end if
    end if

    POP_SUB(system_arrived_at_barrier)
  end function system_arrived_at_barrier

  ! ---------------------------------------------------------
  logical function system_arrived_at_any_barrier(this)
    class(system_t), intent(inout) :: this

    integer :: ii

    PUSH_SUB(system_arrived_at_any_barrier)

    system_arrived_at_any_barrier = .false.
    do ii = 1, NUMBER_BARRIERS
      system_arrived_at_any_barrier = system_arrived_at_any_barrier &
        .or. this%arrived_at_barrier(ii)
    end do

    POP_SUB(system_arrived_at_any_barrier)
  end function system_arrived_at_any_barrier


  ! ---------------------------------------------------------
  !> Calculate the potential energy of the system.
  !! The potential energy is defined as the sum of all energies
  !! arising from interactions with external systems.
  !! (Note that multisystems override this function)
  subroutine system_update_potential_energy(this)
    class(system_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(system_update_potential_energy)

    this%potential_energy = M_ZERO

    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      if(.not. interaction%intra_interaction) then
        call interaction%calculate_energy()
        this%potential_energy = this%potential_energy + interaction%energy
      end if
    end do

    POP_SUB(system_update_potential_energy)
  end subroutine system_update_potential_energy

  ! ---------------------------------------------------------
  !> Calculate the internal energy of the system.
  !! The internal energy is defined as the sum of all energies
  !! arising from intra-interactions and the entropy terms (if available).
  !! (Note that multisystems override this function)
  subroutine system_update_internal_energy(this)
    class(system_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(system_update_internal_energy)

    this%internal_energy = M_ZERO
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      if(interaction%intra_interaction) then
        call interaction%calculate_energy()
        this%internal_energy = this%internal_energy + interaction%energy
      end if
    end do

    POP_SUB(system_update_internal_energy)
  end subroutine system_update_internal_energy

  ! ---------------------------------------------------------
  !> Calculate the total energy of the system.
  !! The total energy is defined as the sum of
  !! the kinetic, the internal and the potential energies.
  subroutine system_update_total_energy(this)
    class(system_t), intent(inout) :: this

    PUSH_SUB(system_update_total_energy)

    call this%update_kinetic_energy()
    this%total_energy = this%kinetic_energy

    !> External energy as the sum over interaction energies
    call this%update_potential_energy()
    this%total_energy = this%total_energy + this%potential_energy

    !> Self energy arises from the interaction with itself
    call this%update_internal_energy()
    this%total_energy = this%total_energy + this%internal_energy

    POP_SUB(system_update_total_energy)
  end subroutine system_update_total_energy

end module system_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
