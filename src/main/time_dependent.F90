!! Copyright (C) 2019-2020 M. Oliveira, H. Appel, N. Tancogne-Dejean
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

module time_dependent_oct_m
  use debug_oct_m
  use electrons_oct_m
  use global_oct_m
  use messages_oct_m
  use multisystem_basic_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use system_oct_m
  use td_oct_m
  use unit_system_oct_m
  use walltimer_oct_m

  implicit none

  private
  public :: time_dependent_run

contains

  ! ---------------------------------------------------------
  subroutine time_dependent_run(system, from_scratch)
    class(*), intent(inout) :: system
    logical,  intent(inout) :: from_scratch

    PUSH_SUB(time_dependent_run)

    select type (system)
    class is (multisystem_basic_t)
      call time_dependent_run_multisystem(system, from_scratch)
    type is (electrons_t)
      call time_dependent_run_legacy(system, from_scratch)
    end select

    POP_SUB(time_dependent_run)
  end subroutine time_dependent_run

  ! ---------------------------------------------------------
  subroutine time_dependent_run_multisystem(systems, from_scratch)
    class(multisystem_basic_t), intent(inout) :: systems
    logical,                    intent(in)    :: from_scratch

    FLOAT :: final_time
    logical :: trigger_restart, stop_code, stop_loop
    logical :: restart_read

    PUSH_SUB(time_dependent_run_multisystem)

    call multisystem_debug_init("debug/multisystem_propagation.log", global_namespace, systems%grp)

    call messages_write('Info: Running Multi-System time evolution')
    call messages_info(namespace=systems%namespace)

    ! Get final propagation time from input
    ! This variable is also defined (and properly documented) in td/td.F90.
    ! This is temporary, until all the propagators are moved to the new framework.
    call parse_variable(systems%namespace, 'TDPropagationTime', CNST(-1.0), final_time, unit = units_inp%time)
    if (final_time <= M_ZERO) then
      call messages_input_error(systems%namespace, 'TDPropagationTime', 'must be greater than zero')
    end if
    call messages_print_var_value('TDPropagationTime', final_time)

    ! Initialize all propagators
    call systems%init_propagator()

    ! Read restart files or set initial conditions
    if (.not. from_scratch) then
      restart_read = systems%restart_read()
    else
      restart_read = .false.
    end if
    if (.not. restart_read) then
      call systems%initial_conditions()
    end if

    call systems%propagation_start()

    call multisystem_debug_write_marker(systems%namespace, event_marker_t("propagation_start"))

    ! The full TD loop
    stop_loop = .false.
    call systems%start_barrier(final_time, BARRIER_TIME)
    do while (.not. systems%arrived_at_barrier(BARRIER_TIME))
      ! Do one algorithmic operation
      call systems%dt_operation()

      ! determine cases in which to trigger writing restart files
      stop_code = clean_stop(systems%grp%comm) .or. walltimer_alarm(systems%grp%comm)
      trigger_restart = stop_code .or. restart_walltime_period_alarm(systems%grp%comm)

      if (trigger_restart .and. .not. stop_loop) then
        call systems%start_barrier(systems%next_time_on_largest_dt(), BARRIER_RESTART)
        stop_loop = stop_code
        trigger_restart = .false.
      end if
      if (systems%arrived_at_barrier(BARRIER_RESTART)) then
        call systems%restart_write()
        call systems%end_barrier(BARRIER_RESTART)
        if (stop_loop) exit
      end if
    end do

    if (systems%arrived_at_barrier(BARRIER_TIME)) then
      call systems%restart_write()
    end if

    call multisystem_debug_write_marker(systems%namespace, event_marker_t("propagation_finish"))

    call systems%propagation_finish()

    call multisystem_debug_end()

    POP_SUB(time_dependent_run_multisystem)
  end subroutine time_dependent_run_multisystem

  ! ---------------------------------------------------------
  subroutine time_dependent_run_legacy(electrons, from_scratch)
    class(electrons_t), intent(inout) :: electrons
    logical,            intent(inout) :: from_scratch

    PUSH_SUB(time_dependent_run_legacy)

    call td_init(electrons%td, electrons%namespace, electrons%space, electrons%gr, electrons%ions, electrons%st, electrons%ks, &
      electrons%hm, electrons%ext_partners, electrons%outp)
    call td_init_run(electrons%td, electrons%namespace, electrons%mc, electrons%gr, electrons%ions, electrons%st, electrons%ks, &
      electrons%hm, electrons%ext_partners, electrons%outp, electrons%space, from_scratch)
    call td_run(electrons%td, electrons%namespace, electrons%mc, electrons%gr, electrons%ions, electrons%st, electrons%ks, &
      electrons%hm, electrons%ext_partners, electrons%outp, electrons%space, from_scratch)
    call td_end_run(electrons%td, electrons%st, electrons%hm)
    call td_end(electrons%td)

    POP_SUB(time_dependent_run_legacy)
  end subroutine time_dependent_run_legacy

end module time_dependent_oct_m
