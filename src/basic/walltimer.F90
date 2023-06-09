!! Copyright (C) 2019 M. Lueders, Heiko Appel
!! Copyright (C) 2021 I-Te Lu
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

!> This module provices a simple timer class which can be used to trigger
!! the writing of a restart file before the requested CPU time is up.
!!
!! It allows to take into account the time required for one iteration and
!! optionally a time margin for completing the restart dump process.

#include "global.h"

module walltimer_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m

  implicit none

  private

  FLOAT :: start_time       !< time the timer was started
  FLOAT :: last_tap         !< time of the last call to tap()
  FLOAT :: iteration_time   !< time difference of two calls to tap()
  FLOAT :: margin           !< additional time margin for writing the restart file
  FLOAT :: duration         !< time when the alarm should trigger
  FLOAT :: restart_last_tap !< time of the last call for restart
  FLOAT :: restart_duration !< time period to restart (trigger writing restart file)

  logical :: active         !< if .false. the timer will not issue an alarm.
  logical :: auto_tap       !< if .true., tap() is automatically called in every wakeUp() call.
  logical :: restart_active !< if .false. the timer will not issue an alarm for the restart period timer

  public ::         &
    walltimer_init, &
    walltimer_end,  &
    walltimer_tap,  &
    walltimer_alarm, &
    restart_walltime_period_alarm

contains

  !> initialize the timer
  subroutine walltimer_init(auto)
    logical, optional,   intent(in)   :: auto   !< automatically call walltimer_tap in walltimer_alarm() if .true.

    FLOAT  :: alarm_time, write_time, restart_alarm_time

    PUSH_SUB(walltimer_init)

    start_time = M_ZERO
    last_tap = M_ZERO
    iteration_time = M_ZERO
    margin = M_ZERO

    active = .false.
    auto_tap = optional_default(auto, .true.)

    ! The following have to be moved to the right place, after the names for the variables have been confirmed:

    !%Variable Walltime
    !%Type float
    !%Default 0
    !%Section Execution::IO
    !%Description
    !% Time in minutes before which the restart file will be written. This is to make sure that at least one restart
    !% file can be written before the code might be killed to to exceeding the given CPU time.
    !% If a finite time (in minutes) is specified, the code will write the restart file when the next
    !% iteration (plus the RestartWriteTime) would exceed the given time.
    !% A value less than 1 second (1/60 minutes) will disable the timer.
    !%End
    call parse_variable(global_namespace, 'Walltime', M_ZERO, alarm_time)
    call set_alarm(alarm_time*CNST(60.0))

    !%Variable RestartWriteTime
    !%Type float
    !%Default 5
    !%Section Execution::IO
    !%Description
    !% The RestartWriteTime (in minutes) will be subtracted from the WallTime to allow time for writing the restart file.
    !% In huge calculations, this value should be increased.
    !%End
    call parse_variable(global_namespace, 'RestartWriteTime', CNST(5.0), write_time)
    if (write_time > alarm_time/M_FOUR) write_time = alarm_time/M_FOUR
    call set_margin(write_time*CNST(60.0))

    !%Variable RestartWallTimePeriod
    !%Type float
    !%Default 120
    !%Section Execution::IO
    !%Description
    !% Period Time (in minutes) at which the restart file will be written.
    !% If a finite time (in minutes) is specified, the code will write the restart file every period.
    !%End
    call parse_variable(global_namespace, 'RestartWallTimePeriod', CNST(120.0), restart_alarm_time)
    call set_restart_alarm(restart_alarm_time*CNST(60.0))
    restart_active = .false.

    call start()

    POP_SUB(walltimer_init)
  end subroutine walltimer_init

  !> destructor
  subroutine walltimer_end()

    PUSH_SUB(walltimer_end)

    active = .false.

    POP_SUB(walltimer_end)
  end subroutine walltimer_end

  !> set alarm interval in seconds
  subroutine set_alarm(time)
    FLOAT :: time

    PUSH_SUB(set_alarm)

    duration = time

    POP_SUB(set_alarm)
  end subroutine set_alarm

  !> set restart alarm interval in seconds
  subroutine set_restart_alarm(time)
    FLOAT :: time

    PUSH_SUB(set_restart_alarm)

    restart_duration = time

    POP_SUB(set_restart_alarm)
  end subroutine set_restart_alarm

  !> set safty margin in seconds
  subroutine set_margin(time)

    FLOAT :: time

    PUSH_SUB(set_margin)

    margin = time

    POP_SUB(set_margin)
  end subroutine set_margin

  !> start the timer (save starting time)
  subroutine start()

    PUSH_SUB(start)

    start_time = loct_clock()
    last_tap = start_time
    restart_last_tap = start_time

    !> disable timer if duration is less than one second.
    if (duration > M_ONE) active = .true.
    !> disable restart period timer if period is less than one minute.
    if (restart_duration > CNST(60.0)) restart_active = .true.

    POP_SUB(start)
  end subroutine start

  !> measure time of one iteration
  subroutine walltimer_tap(print)
    logical, optional, intent(in) :: print

    FLOAT :: now

    PUSH_SUB(walltimer_tap)

    now = loct_clock()

    iteration_time = now - last_tap
    last_tap = now

    if (optional_default(print, .false.)) then
      write(message(1), '("Walltimer_tap:   elapsed time = ",F6.2," (", 3F10.5, "), active = ",L1 )')  &
        now - start_time, duration, iteration_time, margin, active
      call messages_info(1, all_nodes=.true.)
    end if

    POP_SUB(walltimer_tap)
  end subroutine walltimer_tap

  !> indicate whether time is up
  logical function walltimer_alarm(comm, print)
    integer,           intent(in) :: comm
    logical, optional, intent(in) :: print

#ifdef HAVE_MPI
    logical :: alarm
#endif
    FLOAT :: now

    PUSH_SUB(walltimer_alarm)

    now = loct_clock()

    if (optional_default(print, .false.)) then
      write(message(1), '("Walltimer_alarm: elapsed time = ",F6.2," (", 3F10.5, "), active = ",L1 )')  &
        now - start_time, duration, iteration_time, margin, active
      call messages_info(1, all_nodes=.true.)
    end if

    if (auto_tap) call walltimer_tap()

    walltimer_alarm = active .and. (now > start_time + duration - iteration_time - margin)

    ! All processes need to know that the walltimer was triggered.
#ifdef HAVE_MPI
    call MPI_Allreduce(walltimer_alarm, alarm, 1, MPI_LOGICAL, MPI_LOR, comm, mpi_err)
    walltimer_alarm = alarm
#endif

    if (walltimer_alarm) then
      write(message(1), '("Walltimer stopping execution after = ",F6.2," minutes.")') (now - start_time)/CNST(60.0)
      call messages_info(1)
      ! Switch status to indicate that the walltimer aborted the calculation.
      ! This can be used to communicate the status of octopus to a queuing system and to daisy-chain jobs.
      call messages_switch_status('walltimer-aborted')
    end if

    POP_SUB(walltimer_alarm)
  end function walltimer_alarm

  logical function restart_walltime_period_alarm(comm)
    integer,           intent(in) :: comm

#ifdef HAVE_MPI
    logical :: alarm
#endif
    FLOAT :: now

    PUSH_SUB(restart_walltime_period_alarm)

    now = loct_clock()

    restart_walltime_period_alarm = restart_active .and. (now > restart_last_tap + restart_duration)

    ! All processes need to know that the walltimer was triggered.
#ifdef HAVE_MPI
    call MPI_Allreduce(restart_walltime_period_alarm, alarm, 1, MPI_LOGICAL, MPI_LOR, comm, mpi_err)
    restart_walltime_period_alarm = alarm
#endif

    if (restart_walltime_period_alarm) then
      restart_last_tap = now
      write(message(1), '(a)') "Restart walltime period is reached: writing restart data"
      call messages_info(1)
    end if

    POP_SUB(restart_walltime_period_alarm)
  end function restart_walltime_period_alarm

end module walltimer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
