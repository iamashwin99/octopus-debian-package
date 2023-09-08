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

module propagator_aetrs_oct_m
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
    propagator_aetrs_t

  type, extends(propagator_t) :: propagator_aetrs_t
    private
  end type propagator_aetrs_t

  interface propagator_aetrs_t
    procedure propagator_aetrs_constructor
  end interface propagator_aetrs_t

  !# doc_start aetrs_propagation_operations
  ! Specific exponential mid-point propagation operations identifiers
  character(len=ALGO_LABEL_LEN), public, parameter :: &
    AETRS_START = 'AETRS_START', &
    AETRS_FINISH = 'AETRS_FINISH', &
    AETRS_FIRST_HALF = 'AETRS_FIRST_HALF', &
    AETRS_SECOND_HALF = 'AETRS_SECOND_HALF', &
    AETRS_EXTRAPOLATE = 'AETRS_EXTRAPOLATE'

  ! Specific exponential mid-point propagation operations
  type(algorithmic_operation_t), public, parameter :: &
    OP_AETRS_START = algorithmic_operation_t(AETRS_START, 'AETRS: start'), &
    OP_AETRS_FINISH = algorithmic_operation_t(AETRS_FINISH, 'AETRS: finish'), &
    OP_AETRS_FIRST_HALF = algorithmic_operation_t(AETRS_FIRST_HALF, 'AETRS: first half'), &
    OP_AETRS_SECOND_HALF = algorithmic_operation_t(AETRS_SECOND_HALF, 'AETRS: second half'), &
    OP_AETRS_EXTRAPOLATE = algorithmic_operation_t(AETRS_EXTRAPOLATE, 'AETRS: extrapolate to dt')
  !# doc_end

contains

  ! ---------------------------------------------------------
  function propagator_aetrs_constructor(dt) result(this)
    FLOAT,                     intent(in) :: dt
    type(propagator_aetrs_t), pointer   :: this

    PUSH_SUB(propagator_aetrs_constructor)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = .false.
    this%start_step = OP_AETRS_START
    this%final_step = OP_AETRS_FINISH

    call this%add_operation(OP_AETRS_FIRST_HALF)
    call this%add_operation(OP_AETRS_EXTRAPOLATE)
    call this%add_operation(OP_AETRS_SECOND_HALF)
    call this%add_operation(OP_UPDATE_INTERACTIONS)
    call this%add_operation(OP_STEP_DONE)
    call this%add_operation(OP_REWIND_ALGORITHM)

    this%algo_steps = 1

    this%dt = dt

    POP_SUB(propagator_aetrs_constructor)
  end function propagator_aetrs_constructor

end module propagator_aetrs_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
