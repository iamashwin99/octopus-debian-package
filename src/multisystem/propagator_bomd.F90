!! Copyright (C) 2023 N. Tancogne-Dejean, M. Lueders, A. Buccheri
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

module propagator_bomd_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m

  implicit none

  private
  public ::                            &
    propagator_bomd_t

  type, extends(propagator_t) :: propagator_bomd_t
    private
  end type propagator_bomd_t

  interface propagator_bomd_t
    procedure propagator_bomd_constructor
  end interface propagator_bomd_t

  !# doc_start bomd_propagation_operations
  ! Specific exponential mid-point propagation operations identifiers
  character(len=ALGO_LABEL_LEN), public, parameter :: &
    BOMD_START        = 'BOMD_START',                 &
    BOMD_FINISH       = 'BOMD_FINISH',                &
    BOMD_ELEC_SCF     = 'BOMD_ELEC_SCF'

  ! Specific exponential mid-point propagation operations
  type(algorithmic_operation_t), public, parameter :: &
    OP_BOMD_START        = algorithmic_operation_t(BOMD_START,        'Starting Born-Oppenheimer MD'),  &
    OP_BOMD_FINISH       = algorithmic_operation_t(BOMD_FINISH,       'Finishing Born-Oppenheimer MD'), &
    OP_BOMD_ELEC_SCF     = algorithmic_operation_t(BOMD_ELEC_SCF,     'SCF for the electrons')
  !# doc_end

contains

  ! ---------------------------------------------------------
  function propagator_bomd_constructor(dt) result(this)
    FLOAT,                  intent(in) :: dt
    type(propagator_bomd_t), pointer   :: this

    PUSH_SUB(propagator_bomd_constructor)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = .false.
    this%start_step = OP_BOMD_START
    this%final_step = OP_BOMD_FINISH

    call this%add_operation(OP_VERLET_UPDATE_POS)
    call this%add_operation(OP_UPDATE_INTERACTIONS)
    call this%add_operation(OP_BOMD_ELEC_SCF)
    call this%add_operation(OP_VERLET_COMPUTE_ACC)
    call this%add_operation(OP_VERLET_COMPUTE_VEL)
    call this%add_operation(OP_STEP_DONE)
    call this%add_operation(OP_REWIND_ALGORITHM)

    this%algo_steps = 1
    this%dt = dt

    POP_SUB(propagator_bomd_constructor)
  end function propagator_bomd_constructor

end module propagator_bomd_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
