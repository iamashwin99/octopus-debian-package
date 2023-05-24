!! Copyright (C) 2019-2020 M. Oliveira, Heiko Appel
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

module multisystem_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use ghost_interaction_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use interaction_with_partner_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_oct_m
  use propagator_static_oct_m
  use system_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    multisystem_t,        &
    multisystem_init,     &
    multisystem_end

  type, extends(system_t), abstract :: multisystem_t
    type(system_list_t) :: list
  contains
    procedure :: dt_operation =>  multisystem_dt_operation
    procedure :: init_parallelization => multisystem_init_parallelization
    procedure :: smallest_algo_dt => multisystem_smallest_algo_dt
    procedure :: next_time_on_largest_dt => multisystem_next_time_on_largest_dt
    procedure :: reset_clocks => multisystem_reset_clocks
    procedure :: init_propagator => multisystem_init_propagator
    procedure :: propagation_start => multisystem_propagation_start
    procedure :: propagation_finish => multisystem_propagation_finish
    procedure :: has_reached_final_propagation_time => multisystem_has_reached_final_propagation_time
    procedure :: init_all_interactions => multisystem_init_all_interactions
    procedure :: init_interaction => multisystem_init_interaction
    procedure :: write_interaction_graph => multisystem_write_interaction_graph
    procedure :: initial_conditions => multisystem_initial_conditions
    procedure :: do_td_operation => multisystem_do_td_operation
    procedure :: is_tolerance_reached => multisystem_is_tolerance_reached
    procedure :: update_quantity => multisystem_update_quantity
    procedure :: update_exposed_quantity => multisystem_update_exposed_quantity
    procedure :: init_interaction_as_partner => multisystem_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => multisystem_copy_quantities_to_interaction
    procedure :: process_is_slave => multisystem_process_is_slave
    procedure :: start_barrier => multisystem_start_barrier
    procedure :: end_barrier => multisystem_end_barrier
    procedure :: arrived_at_barrier => multisystem_arrived_at_barrier
    procedure :: restart_write => multisystem_restart_write
    procedure :: restart_read => multisystem_restart_read
    procedure :: restart_write_data => multisystem_restart_write_data
    procedure :: restart_read_data => multisystem_restart_read_data
    procedure :: update_kinetic_energy => multisystem_update_kinetic_energy
    procedure :: update_potential_energy => multisystem_update_potential_energy
    procedure :: update_internal_energy => multisystem_update_internal_energy
    procedure :: get_flat_list => multisystem_get_flat_list
  end type multisystem_t

contains

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_init(this, namespace, factory)
    class(multisystem_t),      intent(inout) :: this
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory

    integer :: isys, system_type, ic
    character(len=128) :: system_name
    type(block_t) :: blk

    PUSH_SUB(multisystem_init)

    this%namespace = namespace

    if (parse_block(this%namespace, factory%block_name(), blk) == 0) then

      do isys = 1, parse_block_n(blk)
        ! Parse system name and type
        call parse_block_string(blk, isys - 1, 0, system_name)
        if (len_trim(system_name) == 0) then
          call messages_input_error(this%namespace, factory%block_name(), 'All systems must have a name')
        end if
        do ic = 1, len(parser_varname_excluded_characters)
          if (index(trim(system_name), parser_varname_excluded_characters(ic:ic)) /= 0) then
            call messages_input_error(this%namespace, factory%block_name(), &
              'Illegal character "' // parser_varname_excluded_characters(ic:ic) // '" in system name', row=isys-1, column=0)
          end if
        end do
        call parse_block_integer(blk, isys - 1, 1, system_type)

        call multisystem_create_system(this, system_name, system_type, isys, factory)
      end do
      call parse_block_end(blk)
    else
      message(1) = "Input error while reading block "//trim(this%namespace%get())//"."//trim(factory%block_name())
      call messages_fatal(1, namespace=this%namespace)
    end if

    POP_SUB(multisystem_init)
  end subroutine multisystem_init

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_create_system(this, system_name, system_type, isys, factory)
    class(multisystem_t),      intent(inout) :: this
    character(len=128),           intent(in) :: system_name
    integer,                      intent(in) :: system_type
    integer,                      intent(in) :: isys
    class(system_factory_abst_t), intent(in) :: factory

    type(system_iterator_t) :: iter
    class(system_t), pointer :: sys, other

    PUSH_SUB(multisystem_create_system)

    ! Create folder to store system files.
    ! Needs to be done before creating the system as this in turn might create subfolders.
    call io_mkdir(system_name, namespace=this%namespace)

    ! Create system
    sys => factory%create(this%namespace, system_name, system_type)
    if (.not. associated(sys)) then
      call messages_input_error(this%namespace, factory%block_name(), 'Unknown system type.')
    end if

    ! Check that the system is unique
    call iter%start(this%list)
    do while (iter%has_next())
      other => iter%get_next()
      if (sys%namespace == other%namespace) then
        call messages_input_error(this%namespace, factory%block_name(), 'Duplicated system in multisystem', &
          row=isys-1, column=0)
      end if
    end do

    ! Add system to list of systems
    call this%list%add(sys)

    POP_SUB(multisystem_create_system)
  end subroutine multisystem_create_system

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_init_parallelization(this, grp)
    class(multisystem_t), intent(inout) :: this
    type(mpi_grp_t),      intent(in)    :: grp

    type(system_iterator_t) :: iter
    class(system_t), pointer :: sys
    type(mpi_grp_t) :: sys_grp

    PUSH_SUB(multisystem_init_parallelization)

    call system_init_parallelization(this, grp)

    ! Now parallelize over systems in this multisystem
    call iter%start(this%list)
    do while (iter%has_next())
      sys => iter%get_next()
      ! for now, duplicate communicator - more complicated parallelization schemes can be implemented here
      call mpi_grp_duplicate(sys_grp, grp)
      call sys%init_parallelization(sys_grp)
    end do

    POP_SUB(multisystem_init_parallelization)
  end subroutine multisystem_init_parallelization

  ! ---------------------------------------------------------------------------------------
  recursive function multisystem_smallest_algo_dt(this) result(smallest_algo_dt)
    class(multisystem_t), intent(inout) :: this
    FLOAT :: smallest_algo_dt

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_smallest_algo_dt)

    smallest_algo_dt = M_HUGE
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      select type (system)
      class is (multisystem_t)
        smallest_algo_dt = min(smallest_algo_dt, system%smallest_algo_dt())
      class default
        smallest_algo_dt = min(smallest_algo_dt, system%prop%dt/system%prop%algo_steps)
      end select
    end do

    POP_SUB(multisystem_smallest_algo_dt)
  end function multisystem_smallest_algo_dt

  ! ---------------------------------------------------------------------------------------
  recursive function multisystem_next_time_on_largest_dt(this) result(next_time_on_largest_dt)
    class(multisystem_t), intent(inout) :: this
    FLOAT :: next_time_on_largest_dt

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system
    type(clock_t) :: clock

    PUSH_SUB(multisystem_next_time_on_largest_dt)

    next_time_on_largest_dt = M_ZERO
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      select type (system)
      class is (multisystem_t)
        next_time_on_largest_dt = max(next_time_on_largest_dt, system%next_time_on_largest_dt())
      class default
        clock = system%clock + CLOCK_TICK
        next_time_on_largest_dt = max(next_time_on_largest_dt, clock%time())
      end select
    end do

    POP_SUB(multisystem_next_time_on_largest_dt)
  end function multisystem_next_time_on_largest_dt

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_dt_operation(this)
    class(multisystem_t),     intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    type(event_handle_t) :: debug_handle

    PUSH_SUB(multisystem_dt_operation)

    if (debug%info) then
      write(message(1), '(a,a,1X,a)') "Debug: Start multisystem_dt_operation for '" + trim(this%namespace%get()) + "'"
      call messages_info(1, namespace=this%namespace)
    end if
    debug_handle = multisystem_debug_write_event_in(this%namespace, event_function_call_t("multisystem_dt_operation"), &
      system_clock = this%clock, prop_clock = this%prop%clock)

    ! Multisystem
    call system_dt_operation(this)

    ! Subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%dt_operation()
    end do

    if (debug%info) then
      write(message(1), '(a,a,1X,a)') "Debug: Finish multisystem_dt_operation for '" + trim(this%namespace%get()) + "'"
      call messages_info(1, namespace=this%namespace)
    end if

    call multisystem_debug_write_event_out(debug_handle, system_clock = this%clock, prop_clock = this%prop%clock)

    POP_SUB(multisystem_dt_operation)
  end subroutine multisystem_dt_operation

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_reset_clocks(this, accumulated_ticks)
    class(multisystem_t),      intent(inout) :: this
    integer,                   intent(in)    :: accumulated_ticks

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_reset_clocks)

    ! Multisystem clocks
    call system_reset_clocks(this, accumulated_ticks)

    ! Subsystems clocks
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%reset_clocks(accumulated_ticks)
    end do

    POP_SUB(multisystem_reset_clocks)
  end subroutine multisystem_reset_clocks

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_init_propagator(this)
    class(multisystem_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system
    type(interaction_iterator_t) :: inter_iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(multisystem_init_propagator)

    ! Now initialized the propagators of the subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%init_propagator()
    end do

    ! Initialize the propagator of the multisystem. By default the
    ! multisystem itself and its own quantities are kept unchaged
    ! by using the static propagator. However, the subsystems are allowed to have
    ! their own propagators and those do not have to be static.
    ! Needs to be done after initializing the subsystems propagators,
    ! as we use the smallest dt of the subsystems.
    this%prop => propagator_static_t(this%smallest_algo_dt())
    this%interaction_timing = OPTION__INTERACTIONTIMING__TIMING_EXACT
    call this%prop%rewind()

    ! Initialize propagator clock
    this%prop%clock = clock_t(time_step=this%prop%dt/this%prop%algo_steps)

    ! Initialize system clock
    this%clock = clock_t(time_step=this%prop%dt)

    ! Interaction clocks
    call inter_iter%start(this%interactions)
    do while (inter_iter%has_next())
      interaction => inter_iter%get_next()
      interaction%clock = this%prop%clock - CLOCK_TICK
    end do

    ! Required quantities clocks
    where (this%quantities%required)
      this%quantities%clock = this%prop%clock
    end where

    POP_SUB(multisystem_init_propagator)
  end subroutine multisystem_init_propagator

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_propagation_start(this)
    class(multisystem_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    type(event_handle_t) :: debug_handle

    PUSH_SUB(multisystem_propagation_start)

    debug_handle = multisystem_debug_write_event_in(this%namespace, event_function_call_t("multisystem_propagation_start"), &
      system_clock = this%clock, prop_clock = this%prop%clock)

    ! Start the propagation of the multisystem
    call system_propagation_start(this)

    ! Now start the propagation of the subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%propagation_start()
    end do

    call multisystem_debug_write_event_out(debug_handle, system_clock = this%clock, prop_clock = this%prop%clock)

    POP_SUB(multisystem_propagation_start)
  end subroutine multisystem_propagation_start

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_propagation_finish(this)
    class(multisystem_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    type(event_handle_t) :: debug_handle

    PUSH_SUB(multisystem_propagation_finish)

    debug_handle = multisystem_debug_write_event_in(this%namespace, event_function_call_t("multisystem_propagation_finish"), &
      system_clock = this%clock, prop_clock = this%prop%clock)

    ! Finish the propagation of the multisystem
    call system_propagation_finish(this)

    ! Now finish the propagation of the subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%propagation_finish()
    end do

    call multisystem_debug_write_event_out(debug_handle, system_clock = this%clock, prop_clock = this%prop%clock)

    POP_SUB(multisystem_propagation_finish)
  end subroutine multisystem_propagation_finish

  ! ---------------------------------------------------------------------------------------
  recursive logical function multisystem_has_reached_final_propagation_time(this, final_time)
    class(multisystem_t),      intent(inout) :: this
    FLOAT,                     intent(in)    :: final_time

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_has_reached_final_propagation_time)

    multisystem_has_reached_final_propagation_time = system_has_reached_final_propagation_time(this, final_time)
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      multisystem_has_reached_final_propagation_time = multisystem_has_reached_final_propagation_time .and. &
        system%has_reached_final_propagation_time(final_time)
    end do

    POP_SUB(multisystem_has_reached_final_propagation_time)
  end function multisystem_has_reached_final_propagation_time

  ! ---------------------------------------------------------
  recursive subroutine multisystem_init_all_interactions(this)
    class(multisystem_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter_i
    class(interaction_t), pointer :: interaction
    type(system_iterator_t) :: iter_s
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_init_all_interactions)

    ! Initialize interactions directly owned by the multisystem
    call iter_i%start(this%interactions)
    do while (iter_i%has_next())
      interaction => iter_i%get_next()
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

    ! Initialize interactions owned by the subsystems
    call iter_s%start(this%list)
    do while (iter_s%has_next())
      system => iter_s%get_next()
      call system%init_all_interactions()
    end do

    POP_SUB(multisystem_init_all_interactions)
  end subroutine multisystem_init_all_interactions

  ! ---------------------------------------------------------
  subroutine multisystem_init_interaction(this, interaction)
    class(multisystem_t), target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(multisystem_init_interaction)

    ! The multitystem class should never know about any specific interaction.
    ! Only classes that extend it can know about specific interactions.
    ! Such classes should override this method to add new supported interactions.
    message(1) = "Trying to initialize an interaction in the multisystem class"
    call messages_fatal(1, namespace=this%namespace)

    POP_SUB(multisystem_init_interaction)
  end subroutine multisystem_init_interaction

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_write_interaction_graph(this, iunit, include_ghosts)
    class(multisystem_t), intent(in) :: this
    integer,              intent(in) :: iunit
    logical,              intent(in) :: include_ghosts

    class(system_t), pointer :: system
    class(interaction_t), pointer :: interaction
    type(system_iterator_t) :: sys_iter
    type(interaction_iterator_t) :: inter_iter

    PUSH_SUB(multisystem_write_interaction_graph)

    ! Loop over all the subsystems
    call sys_iter%start(this%list)
    do while (sys_iter%has_next())
      system => sys_iter%get_next()

      ! Loop over the interactions that this subsystem has
      call inter_iter%start(system%interactions)
      do while (inter_iter%has_next())
        interaction => inter_iter%get_next()

        ! Write interaction to DOT graph if this interaction has a partner
        select type (interaction)
        type is (ghost_interaction_t)
          if (include_ghosts) then
            write(iunit, '(2x,a)') '"' + trim(system%namespace%get()) + '" -> "' + trim(interaction%partner%namespace%get()) + &
              '" [label="'+ interaction%label + '"];'
          end if
          ! Do not include systems connected by ghost interactions
        class is (interaction_with_partner_t)
          write(iunit, '(2x,a)') '"' + trim(system%namespace%get()) + '" -> "' + trim(interaction%partner%namespace%get()) + &
            '" [label="'+ interaction%label + '"];'
        end select
      end do

      ! If this subsystem is also a multisystem, then we also need to traverse it
      select type (system)
      class is (multisystem_t)
        call system%write_interaction_graph(iunit, include_ghosts)
      end select
    end do

    POP_SUB(multisystem_write_interaction_graph)
  end subroutine multisystem_write_interaction_graph

  ! ---------------------------------------------------------
  recursive subroutine multisystem_initial_conditions(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_initial_conditions)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%initial_conditions()
    end do

    POP_SUB(multisystem_initial_conditions)
  end subroutine multisystem_initial_conditions

  ! ---------------------------------------------------------
  subroutine multisystem_do_td_operation(this, operation)
    class(multisystem_t),           intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    PUSH_SUB(multisystem_do_td_operation)

    select case (operation%id)
    case (SKIP)
      ! Nothing to do
    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(multisystem_do_td_operation)
  end subroutine multisystem_do_td_operation

  ! ---------------------------------------------------------
  recursive logical function multisystem_is_tolerance_reached(this, tol) result(converged)
    class(multisystem_t), intent(in)    :: this
    FLOAT,                intent(in)    :: tol

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_is_tolerance_reached)

    converged = .true.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      if (.not. system%is_tolerance_reached(tol)) converged = .false.
    end do

    POP_SUB(multisystem_is_tolerance_reached)
  end function multisystem_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine multisystem_update_quantity(this, iq)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: iq

    PUSH_SUB(multisystem_update_quantity)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to update a quantity in the multisystem class"
    call messages_fatal(1, namespace=this%namespace)

    POP_SUB(multisystem_update_quantity)
  end subroutine multisystem_update_quantity

  ! ---------------------------------------------------------
  subroutine multisystem_update_exposed_quantity(partner, iq)
    class(multisystem_t), intent(inout) :: partner
    integer,              intent(in)    :: iq

    PUSH_SUB(multisystem_update_exposed_quantity)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to update an exposed quantity in the multisystem class"
    call messages_fatal(1, namespace=partner%namespace)

    POP_SUB(multisystem_update_exposed_quantity)
  end subroutine multisystem_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine multisystem_init_interaction_as_partner(partner, interaction)
    class(multisystem_t),         intent(in)    :: partner
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(multisystem_init_interaction_as_partner)

    ! The multitystem class should never know about any specific interaction.
    ! Only classes that extend it can know about specific interactions.
    ! Such classes should override this method to add new supported interactions.
    message(1) = "Trying to initialize an interaction as partner in the multisystem class"
    call messages_fatal(1, namespace=partner%namespace)

    POP_SUB(multisystem_init_interaction_as_partner)
  end subroutine multisystem_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine multisystem_copy_quantities_to_interaction(partner, interaction)
    class(multisystem_t),         intent(inout) :: partner
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(multisystem_copy_quantities_to_interaction)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to copy quantities to interaction in the multisystem class"
    call messages_fatal(1, namespace=partner%namespace)

    POP_SUB(multisystem_copy_quantities_to_interaction)
  end subroutine multisystem_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  recursive logical function multisystem_process_is_slave(this) result(is_slave)
    class(multisystem_t), intent(in) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_process_is_slave)

    is_slave = .false.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      if (system%process_is_slave()) is_slave = .true.
    end do

    POP_SUB(multisystem_process_is_slave)
  end function multisystem_process_is_slave

  !--------------------------------------------------------------------
  !> Calculate the kinetic energy:
  !! The kinetic energy of a container (multisystem) is defined by the
  !! kinetic energy with respect to the centre or mass motion.
  recursive subroutine multisystem_update_kinetic_energy(this)
    class(multisystem_t), intent(inout) :: this

    PUSH_SUB(multisystem_update_kinetic_energy)

    ! We currently do not have the centre of mass coordinates implemented for multisystems,
    ! hence we set the kinetic energy to zero.
    ! The kinetic energies of the constituents are contributing to the internal energy.

    this%kinetic_energy = M_ZERO

    POP_SUB(multisystem_update_kinetic_energy)
  end subroutine multisystem_update_kinetic_energy

  !---------------------------------------------------------
  recursive subroutine multisystem_update_internal_energy(this)
    class(multisystem_t), intent(inout) :: this

    class(system_t),              pointer :: system
    class(system_t),              pointer :: system_2
    type(system_iterator_t)      :: system_iter
    type(system_iterator_t)      :: system_iter_2

    PUSH_SUB(multisystem_update_internal_energy)

    ! The internal energy of the multisystem containes the kinetic and internal energies of the consistuents
    !TODO: the kinetic energy wrt the centre of mass motion should be subtracted.

    this%internal_energy = M_ZERO

    call system_iter%start(this%list)
    do while (system_iter%has_next())

      system => system_iter%get_next()

      !> First add the kinetic energies of the subsystems
      call system%update_kinetic_energy()
      this%internal_energy = this%internal_energy + system%kinetic_energy

      !> First add the internal energies of the subsystems
      call system%update_internal_energy()
      this%internal_energy = this%internal_energy + system%internal_energy

      !> Now add the (inter-) interactions between the systems in the container.
      call system_iter_2%start(this%list)
      do while(system_iter_2%has_next())

        system_2 => system_iter_2%get_next()

        !> exclude self-interactions (intra-interactions) as they are included in the internal energy
        !! of the subsystem, which was already added above.
        if(.not. associated(system, system_2)) then
          this%internal_energy = this%internal_energy + multisystem_pair_energy(system, system_2)
        endif
      end do ! system_iter_2

    end do ! system_iter

    POP_SUB(multisystem_update_internal_energy)
  end subroutine multisystem_update_internal_energy

  ! ---------------------------------------------------------
  !> Calculate the potential energy for a container:
  !! The potential energy accounts for the energy of the content of the container
  !! in the field of all systems, which are not part of the container.
  !!
  !! It is important to note, that the potential energy of a container is not the sum
  !! of the potential energies of the constituent system.
  !!
  !! Another complication arises when the container contains another container, as then
  !! we need to account for the interactions of its constituents with partners outside the
  !! current container.
  !!
  subroutine multisystem_update_potential_energy(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t)  :: system_iter
    class(system_t), pointer :: system
    type(interaction_iterator_t)  :: interaction_iter
    class(interaction_t), pointer :: interaction
    type(system_list_t) :: flat_list

    PUSH_SUB(multisystem_update_potential_energy)

    this%potential_energy = M_ZERO

    !> We need to handle interactions of the container itself:
    call system_update_potential_energy(this)

    !> generate a list of all systems inside the container and its subcontainers:
    call this%get_flat_list(flat_list)

    !> loop over all systems inside the container
    call system_iter%start(flat_list)
    do while (system_iter%has_next())

      system => system_iter%get_next()

      !> Even though we are not using the potential energy of the subsystems here, we need to trigger their calculation
      call system%update_potential_energy()

      !> loop over all interactions and discard those with partners inside the container
      call interaction_iter%start(system%interactions)
      do while (interaction_iter%has_next())
        interaction => interaction_iter%get_next()
        select type(interaction)
        class is (interaction_with_partner_t)
          if(.not. flat_list%contains(interaction%partner) .and. .not. interaction%intra_interaction) then
            call interaction%calculate_energy()
            this%potential_energy = this%potential_energy + interaction%energy
          end if
        class is (interaction_t)
          call interaction%calculate_energy()
          this%potential_energy = this%potential_energy + interaction%energy
        end select
      end do

    end do

    POP_SUB(multisystem_update_potential_energy)
  end subroutine multisystem_update_potential_energy

  ! ---------------------------------------------------------
  !> This function calculates the complete interaction energy between partner_A and partner_B, which
  !! means that for any container its constituents will be accounted for. This continues recursively
  !! up to the level of non-container systems.
  !!
  !! If partner_A and partner_B are not interacting, the routine returns zero.
  !!
  recursive FLOAT function multisystem_pair_energy(partner_A, partner_B) result(pair_energy)
    class(interaction_partner_t), intent(in) :: partner_A
    class(interaction_partner_t), intent(in) :: partner_B

    class(system_t), pointer :: system_A
    class(system_t), pointer :: system_B
    type(system_iterator_t) :: system_iterator_A
    type(system_iterator_t) :: system_iterator_B

    PUSH_SUB(multisystem_pair_energy)

    pair_energy = M_ZERO

    select type(partner_A)
    class is (multisystem_t) ! partner_A is container

      call system_iterator_A%start(partner_A%list)
      do while( system_iterator_A%has_next() )

        system_A => system_iterator_A%get_next()

        select type(partner_B)
        class is (multisystem_t)

          call system_iterator_B%start(partner_B%list)
          do while( system_iterator_B%has_next() )
            system_B => system_iterator_B%get_next()
            pair_energy = pair_energy + multisystem_pair_energy(system_A, system_B)
          end do

        class is (system_t)
          pair_energy = pair_energy + interaction_energy(partner_A, partner_B)
        class default
          ASSERT(.false.) ! partner_A must be a system_t
        end select
      end do

    class is (system_t) ! partner_A is non-container system

      select type(partner_B)
      class is (multisystem_t) ! partner_B is container

        call system_iterator_B%start(partner_B%list)
        do while( system_iterator_B%has_next() )
          system_B => system_iterator_B%get_next()
          pair_energy = pair_energy + multisystem_pair_energy(partner_A, system_B)
        end do

      class default ! both partner_A and partner_B are explicit: we need to calculate
        pair_energy = pair_energy + interaction_energy(partner_A, partner_B)
      end select

    class default
      ASSERT(.false.)
    end select

    POP_SUB(multisystem_pair_energy)

  contains

    FLOAT function interaction_energy(system, partner) result (energy)
      class(system_t), target, intent(in) :: system
      class(interaction_partner_t), target, intent(in) :: partner

      type(interaction_iterator_t) :: interaction_iterator
      class(interaction_t), pointer :: interaction

      energy = M_ZERO

      call interaction_iterator%start(system%interactions)
      do while(interaction_iterator%has_next())
        interaction => interaction_iterator%get_next()
        select type(interaction)
        class is (interaction_with_partner_t)
          if( associated(interaction%partner, partner)) then
            call interaction%calculate_energy()
            energy = energy + interaction%energy
          end if
        end select
      end do
    end function interaction_energy

  end function multisystem_pair_energy


  ! ---------------------------------------------------------
  !> Generate a list of all systems contained in a multisystem,
  !! including those inside child containers.
  recursive subroutine multisystem_get_flat_list(this, flat_list)
    class(multisystem_t), intent(in) :: this
    type(system_list_t), intent(out) :: flat_list

    class(interaction_partner_t), pointer :: partner
    type(partner_iterator_t) :: iterator

    PUSH_SUB(multisystem_get_flat_list)

    call iterator%start(this%list)
    do while (iterator%has_next())
      partner => iterator%get_next()

      call flat_list%add(partner)

      select type (partner)
      class is (multisystem_t)
        ! Also include the subsystems of a multisystem
        call partner%get_flat_list(flat_list)
      end select

    end do

    POP_SUB(multisystem_get_flat_list)

  end subroutine multisystem_get_flat_list

  ! ---------------------------------------------------------
  recursive subroutine multisystem_end(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_end)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      SAFE_DEALLOCATE_P(system)
    end do

    call system_end(this)

    POP_SUB(multisystem_end)
  end subroutine multisystem_end

  ! ---------------------------------------------------------
  recursive subroutine multisystem_start_barrier(this, target_time, barrier_index)
    class(multisystem_t), intent(inout) :: this
    FLOAT,                intent(in)    :: target_time
    integer,              intent(in)    :: barrier_index

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_start_barrier)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%start_barrier(target_time, barrier_index)
    end do

    POP_SUB(multisystem_start_barrier)
  end subroutine multisystem_start_barrier

  ! ---------------------------------------------------------
  recursive subroutine multisystem_end_barrier(this, barrier_index)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: barrier_index

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_end_barrier)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%end_barrier(barrier_index)
    end do

    POP_SUB(multisystem_end_barrier)
  end subroutine multisystem_end_barrier

  ! ---------------------------------------------------------
  recursive logical function multisystem_arrived_at_barrier(this, barrier_index)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: barrier_index

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_arrived_at_barrier)

    multisystem_arrived_at_barrier = .true.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      multisystem_arrived_at_barrier = multisystem_arrived_at_barrier .and. &
        system%arrived_at_barrier(barrier_index)
    end do

    POP_SUB(multisystem_arrived_at_barrier)
  end function multisystem_arrived_at_barrier

  ! ---------------------------------------------------------
  recursive subroutine multisystem_restart_write(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_restart_write)

    ! do generic restart steps
    call system_restart_write(this)

    ! loop over all subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%restart_write()
    end do

    POP_SUB(multisystem_restart_write)
  end subroutine multisystem_restart_write

  ! ---------------------------------------------------------
  recursive logical function multisystem_restart_read(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_restart_read)

    ! read generic restart data
    multisystem_restart_read = system_restart_read(this)
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      ! TODO: adapt logics here for consistent restarting
      multisystem_restart_read = multisystem_restart_read .and. &
        system%restart_read()
    end do

    POP_SUB(multisystem_restart_read)
  end function multisystem_restart_read

  ! ---------------------------------------------------------
  subroutine multisystem_restart_write_data(this)
    class(multisystem_t), intent(inout) :: this

    PUSH_SUB(multisystem_restart_write_data)

    ! do not write restart data for multisystem_t

    POP_SUB(multisystem_restart_write_data)
  end subroutine multisystem_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function multisystem_restart_read_data(this)
    class(multisystem_t), intent(inout) :: this

    PUSH_SUB(multisystem_restart_read_data)

    multisystem_restart_read_data = .true.

    POP_SUB(multisystem_restart_read_data)
  end function multisystem_restart_read_data

end module multisystem_oct_m
