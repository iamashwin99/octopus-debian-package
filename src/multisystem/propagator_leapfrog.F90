!! Copyright (C) 2023 Sebastian Ohlmann
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

module propagator_leapfrog_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_oct_m

  implicit none

  private
  public ::                            &
    propagator_leapfrog_t

  type, extends(propagator_t) :: propagator_leapfrog_t
    private
  end type propagator_leapfrog_t

  interface propagator_leapfrog_t
    procedure propagator_leapfrog_constructor
  end interface propagator_leapfrog_t

  !# doc_start leapfrog_propagation_operations
  character(len=ALGO_LABEL_LEN), public, parameter :: &
    LEAPFROG_START        = 'LEAPFROG_START',             &
    LEAPFROG_FINISH       = 'LEAPFROG_FINISH',            &
    LEAPFROG_PROPAGATE    = 'LEAPFROG_PROPAGATE'

  type(algorithmic_operation_t), public, parameter :: &
    OP_LEAPFROG_START        = algorithmic_operation_t(LEAPFROG_START,        'Starting leap frog'),  &
    OP_LEAPFROG_FINISH       = algorithmic_operation_t(LEAPFROG_FINISH,       'Finishing leap frog'), &
    OP_LEAPFROG_PROPAGATE    = algorithmic_operation_t(LEAPFROG_PROPAGATE,    'Propagation step for leap frog')
  !# doc_end

contains

  ! ---------------------------------------------------------
  function propagator_leapfrog_constructor(dt) result(this)
    FLOAT,                  intent(in) :: dt
    type(propagator_leapfrog_t), pointer   :: this

    PUSH_SUB(propagator_leapfrog_constructor)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = .false.
    this%start_step = OP_LEAPFROG_START
    this%final_step = OP_LEAPFROG_FINISH

    call this%add_operation(OP_LEAPFROG_PROPAGATE)
    call this%add_operation(OP_UPDATE_INTERACTIONS)
    call this%add_operation(OP_STEP_DONE)
    call this%add_operation(OP_REWIND_ALGORITHM)

    this%algo_steps = 1
    this%dt = dt

    POP_SUB(propagator_leapfrog_constructor)
  end function propagator_leapfrog_constructor

end module propagator_leapfrog_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
