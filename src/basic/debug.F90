!! Copyright (C) 2016 X. Andrade
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

module debug_oct_m
  use global_oct_m
  use namespace_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use loct_oct_m
  use parser_oct_m

  implicit none

  private
  public ::             &
    debug_t,            &
    debug_init,         &
    debug_enable,       &
    debug_disable,      &
    debug_delete_trace, &
    debug_open_trace,   &
    debug,              &
    epoch_time_diff

  type debug_t
    private
    logical, public :: info
    logical, public :: trace
    logical, public :: trace_term
    logical, public :: trace_file
    logical :: extra_checks
    logical, public :: interaction_graph
    logical, public :: interaction_graph_full
    logical, public :: propagation_graph
    integer :: bits
  end type debug_t

  type(debug_t), save :: debug

  !> max_lun is currently 99, i.e. we can hardwire unit_offset above 1000
  integer, parameter :: unit_offset = 1000

contains

  subroutine debug_init(this, namespace)
    type(debug_t),     intent(out)   :: this
    type(namespace_t), intent(in)    :: namespace

    character(len=256) :: node_hook
    logical :: file_exists, mpi_debug_hook
    integer :: sec, usec

    !%Variable Debug
    !%Type flag
    !%Default no
    !%Section Execution::Debug
    !%Description
    !% This variable controls the amount of debugging information
    !% generated by Octopus. You can use include more than one option
    !% with the + operator.
    !%Option no 0
    !% (default) <tt>Octopus</tt> does not enter debug mode.
    !%Option info 1
    !% Octopus prints additional information to the terminal.
    !%Option trace 2
    !% Octopus generates a stack trace as it enters end exits
    !% subroutines. This information is reported if Octopus stops with
    !% an error.
    !%Option trace_term 4
    !% The trace is printed to the terminal as Octopus enters or exits subroutines. This slows down execution considerably.
    !%Option trace_file 8
    !% The trace is written to files in the <tt>debug</tt>
    !% directory. For each node (when running in parallel) there is a file called
    !% <tt>debug_trace.&lt;rank&gt;</tt>. Writing these files slows down the code by a huge factor and
    !% it is usually only necessary for parallel runs.
    !%Option extra_checks 16
    !% This enables Octopus to perform some extra checks, to ensure
    !% code correctness, that might be too costly for regular runs.
    !%Option interaction_graph 32
    !% Octopus generates a dot file containing the graph for a multisystem run.
    !%Option interaction_graph_full 64
    !% Octopus generates a dot file containing the graph for a multisystem run including ghost interactions.
    !%Option propagation_graph 128
    !% Octopus generates a file with information for the propagation diagram.
    !%End
    call parse_variable(namespace, 'Debug', OPTION__DEBUG__NO, this%bits)

    call from_bits(this)

    call mpi_debug_init(mpi_world%rank, this%info)

    if (this%info) then
      !%Variable MPIDebugHook
      !%Type logical
      !%Default no
      !%Section Execution::Debug
      !%Description
      !% When debugging the code in parallel it is usually difficult to find the origin
      !% of race conditions that appear in MPI communications. This variable introduces
      !% a facility to control separate MPI processes. If set to yes, all nodes will
      !% start up, but will get trapped in an endless loop. In every cycle of the loop
      !% each node is sleeping for one second and is then checking if a file with the
      !% name <tt>node_hook.xxx</tt> (where <tt>xxx</tt> denotes the node number) exists. A given node can
      !% only be released from the loop if the corresponding file is created. This allows
      !% to selectively run, <i>e.g.</i>, a compute node first followed by the master node. Or, by
      !% reversing the file creation of the node hooks, to run the master first followed
      !% by a compute node.
      !%End
      call parse_variable(global_namespace, 'MPIDebugHook', .false., mpi_debug_hook)
      if (mpi_debug_hook) then
        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(stdout,'(a,i6,a,i6.6,20x,a)') '* I ',sec,'.',usec,' | MPI debug hook'

        write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, ' In debug hook'
        write(node_hook,'(i3.3)') mpi_world%rank
        file_exists = .false.

        do while (.not. file_exists)
          inquire(file='node_hook.'//node_hook, exist=file_exists)
          call loct_nanosleep(1,0)
          write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, &
            ' - still sleeping. To release me touch: node_hook.'//trim(node_hook)
        end do

        write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, ' Leaving debug hook'
        ! remove possible debug hooks
        call loct_rm('node_hook.'//trim(node_hook))

        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(stdout,'(a,i6,a,i6.6,20x,a)') '* O ', sec, '.', usec,' | MPI debug hook'
      end if
    end if

  end subroutine debug_init

  !--------------------------------------------------

  subroutine debug_enable(this)
    type(debug_t), intent(inout) :: this

    this%info       = .true.
    this%trace      = .true.
    this%trace_term = .true.
    this%trace_file = .true.
    this%interaction_graph = .true.
    this%interaction_graph_full = .true.
    this%propagation_graph = .true.

  end subroutine debug_enable

  !--------------------------------------------------

  subroutine debug_disable(this)
    type(debug_t), intent(inout) :: this

    call from_bits(this)

  end subroutine debug_disable

  !--------------------------------------------------

  subroutine debug_delete_trace()

    integer :: iunit
    character(len=6) :: filenum

    iunit = mpi_world%rank + unit_offset
    write(filenum, '(i6.6)') iunit - unit_offset
    call loct_mkdir('debug')
    call loct_rm('debug/debug_trace.node.'//filenum)

  end subroutine debug_delete_trace

  ! ---------------------------------------------------------

  subroutine debug_open_trace(iunit)
    integer, intent(out) :: iunit

    character(len=6) :: filenum

    iunit = mpi_world%rank + unit_offset
    write(filenum, '(i6.6)') iunit - unit_offset
    call loct_mkdir('debug')
    open(iunit, file = 'debug/debug_trace.node.'//filenum, &
      action='write', status='unknown', position='append')

  end subroutine debug_open_trace

  ! ---------------------------------------------------------

  subroutine from_bits(this)
    type(debug_t), intent(inout) :: this

    this%info         = (bitand(this%bits, OPTION__DEBUG__INFO)         /= 0)
    this%trace_term   = (bitand(this%bits, OPTION__DEBUG__TRACE_TERM)   /= 0)
    this%trace_file   = (bitand(this%bits, OPTION__DEBUG__TRACE_FILE)   /= 0)
    this%trace        = (bitand(this%bits, OPTION__DEBUG__TRACE)        /= 0) .or. this%trace_term .or. this%trace_file
    this%extra_checks = (bitand(this%bits, OPTION__DEBUG__EXTRA_CHECKS) /= 0) .or. this%trace_term .or. this%trace_file
    this%interaction_graph      = (bitand(this%bits, OPTION__DEBUG__INTERACTION_GRAPH)      /= 0)
    this%interaction_graph_full = (bitand(this%bits, OPTION__DEBUG__INTERACTION_GRAPH_FULL) /= 0)
    this%propagation_graph      = (bitand(this%bits, OPTION__DEBUG__PROPAGATION_GRAPH)      /= 0)

  end subroutine from_bits


  ! ---------------------------------------------------------
  subroutine epoch_time_diff(sec, usec)
    integer, intent(inout) :: sec
    integer, intent(inout) :: usec

    ! this is called by push/pop so there cannot be a push/pop in this routine

    call time_diff(s_epoch_sec, s_epoch_usec, sec, usec)
  end subroutine epoch_time_diff


  ! ---------------------------------------------------------
  !> Computes t2 <- t2-t1. sec1,2 and usec1,2 are
  !! seconds,microseconds of t1,2
  subroutine time_diff(sec1, usec1, sec2, usec2)
    integer, intent(in)    :: sec1
    integer, intent(in)    :: usec1
    integer, intent(inout) :: sec2
    integer, intent(inout) :: usec2

    ! this is called by push/pop so there cannot be a push/pop in this routine

    ! Correct overflow.
    if (usec2 - usec1  <  0) then
      usec2 = 1000000 + usec2
      if (sec2 >= sec1) then
        sec2 = sec2 - 1
      end if
    end if

    ! Replace values.
    if (sec2 >= sec1) then
      sec2 = sec2 - sec1
    end if
    usec2 = usec2 - usec1

  end subroutine time_diff

end module debug_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
