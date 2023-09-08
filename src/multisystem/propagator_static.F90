!! Copyright (C)  2021 M. Oliveira
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

module propagator_static_oct_m
  use algorithm_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use propagator_oct_m

  implicit none

  private
  public ::                       &
    propagator_static_t

  !> Implements a propagator that keeps the state of the system constant.
  !! Note that a time-step is still required to specify at which times
  !! the system quantities are allowed to be exposed.
  type, extends(propagator_t) :: propagator_static_t
    private
  end type propagator_static_t

  interface propagator_static_t
    procedure propagator_static_constructor
  end interface propagator_static_t

contains

  ! ---------------------------------------------------------
  function propagator_static_constructor(dt, nsteps) result(this)
    FLOAT,                     intent(in) :: dt
    integer,                   intent(in) :: nsteps
    type(propagator_static_t), pointer    :: this

    integer :: istep

    PUSH_SUB(propagator_static_constructor)

    SAFE_ALLOCATE(this)

    this%is_static = .true.

    this%start_step = OP_SKIP
    this%final_step = OP_SKIP

    do istep = 1, nsteps
      call this%add_operation(OP_UPDATE_INTERACTIONS)
    end do
    call this%add_operation(OP_STEP_DONE)
    call this%add_operation(OP_REWIND_ALGORITHM)

    this%algo_steps = nsteps

    this%dt = dt

    POP_SUB(propagator_static_constructor)
  end function propagator_static_constructor

end module propagator_static_oct_m

!!o, Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
