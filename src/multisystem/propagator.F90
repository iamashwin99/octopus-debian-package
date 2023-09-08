!! Copyright (C)  2019 N. Tancogne-Dejean
!! Copyright (C)  2020 M. Oliveira
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

!> @brief This module implements the basic propagator framework.
!!
module propagator_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use system_oct_m

  implicit none

  private
  public ::                       &
    propagator_t

  !> @brief Abstract class implementing propagators
  !!
  !! Propagators are implemented as a state machine. This abstract class defines the steps which
  !! are independent of the actual propatation algorithm, and also independent of the system, to
  !! which is is applied.
  type, extends(algorithm_t), abstract :: propagator_t
    private
    class(system_t), pointer, public :: system

    type(algorithm_iterator_t) :: scf_start

    !< Options related to predictor-corrector propagators
    logical, public :: predictor_corrector = .false.
    integer, public :: scf_count
    integer, public :: max_scf_count
    FLOAT, public   :: scf_tol

    FLOAT, public :: final_time = M_ZERO
  contains
    ! Below are the list of operations that needs to be implemented
    procedure :: do_operation => propagator_do_operation       !< @copydoc propagator_oct_m::propagator_do_operation
    procedure :: finished => propagator_finished               !< @copydoc propagator_oct_m::propagator_finished
    procedure :: save_scf_start => propagator_save_scf_start   !< @copydoc propagator_oct_m::propagator_save_scf_start
    procedure :: rewind_scf_loop => propagator_rewind_scf_loop !< @copydoc propagator_oct_m::propagator_rewind_scf_loop
  end type propagator_t

  !# doc_start general_propagation_operations
  ! Known propagation operations
  character(len=ALGO_LABEL_LEN), public, parameter ::   &
    START_SCF_LOOP       = 'START_SCF_LOOP',            &
    END_SCF_LOOP         = 'END_SCF_LOOP',              &
    STORE_CURRENT_STATUS = 'STORE_CURRENT_STATUS'

  type(algorithmic_operation_t), public, parameter :: &
    OP_START_SCF_LOOP       = algorithmic_operation_t(START_SCF_LOOP,       'Starting SCF loop'),         &
    OP_END_SCF_LOOP         = algorithmic_operation_t(END_SCF_LOOP,         'End of SCF iteration'),      &
    OP_STORE_CURRENT_STATUS = algorithmic_operation_t(STORE_CURRENT_STATUS, 'Store current status')
  !# doc_end


contains

  !> @brief perform one operation of the state machine
  !!
  !! this routine performs operations which are general (not system specific).
  !!
  logical function propagator_do_operation(this, operation) result(done)
    class(propagator_t),           intent(inout) :: this
    type(algorithmic_operation_t), intent(in)    :: operation

    done = .true.

    select case (operation%id)
    case (START_SCF_LOOP)
      ASSERT(this%predictor_corrector)

      call this%save_scf_start()
      this%inside_scf = .true.
      this%system%accumulated_loop_ticks = 0

      if (debug%info) then
        write(message(1), '(a,i3,a)') "Debug: -- SCF iter ", this%scf_count, " for '" + trim(this%system%namespace%get()) + "'"
        call messages_info(1, namespace=this%system%namespace)
      end if

    case (END_SCF_LOOP)
      ! Here we first check if we did the maximum number of steps.
      ! Otherwise, we need check the tolerance
      if (this%scf_count == this%max_scf_count) then
        if (debug%info) then
          message(1) = "Debug: -- Max SCF Iter reached for '" + trim(this%system%namespace%get()) + "'"
          call messages_info(1, namespace=this%system%namespace)
        end if
        this%inside_scf = .false.
        call this%next()
      else
        ! We reset the pointer to the beginning of the scf loop
        if (this%system%is_tolerance_reached(this%scf_tol)) then
          if (debug%info) then
            message(1) = "Debug: -- SCF tolerance reached for '" + trim(this%system%namespace%get()) + "'"
            call messages_info(1, namespace=this%system%namespace)
          end if
          this%inside_scf = .false.
          call this%next()
        else
          ! We rewind the instruction stack
          call this%rewind_scf_loop()

          ! We reset the clocks
          call this%system%reset_clocks(this%system%accumulated_loop_ticks)
          this%system%accumulated_loop_ticks = 0
          if (debug%info) then
            write(message(1), '(a,i3,a,a)') "Debug: -- SCF iter ", this%scf_count, " for '" + trim(this%system%namespace%get()), "'"
            call messages_info(1, namespace=this%system%namespace)
          end if
        end if
      end if

    case default
      done = .false.
    end select

  end function propagator_do_operation

  ! ---------------------------------------------------------
  !> @brief indicate whether a propagation has reached the final time
  !!
  logical function propagator_finished(this)
    class(propagator_t), intent(in) :: this

    type(clock_t) :: clock_

    clock_ = this%system%clock + CLOCK_TICK
    propagator_finished = clock_%time() > this%final_time

  end function propagator_finished

  !> @brief Save the current iteration state (START_SCF_LOOP) and move to next step
  !!
  subroutine propagator_save_scf_start(this)
    class(propagator_t), intent(inout) :: this

    PUSH_SUB(propagator_save_scf_start)

    this%scf_start = this%iter
    call this%next()
    this%scf_count = 0

    POP_SUB(propagator_save_scf_start)
  end subroutine propagator_save_scf_start

  ! ---------------------------------------------------------
  !> @brief Reset the iteration state to the beginning of the loop (START_SCF_LOOP) and move to next step
  !!
  subroutine propagator_rewind_scf_loop(this)
    class(propagator_t), intent(inout) :: this

    PUSH_SUB(propagator_rewind_scf_loop)

    this%iter = this%scf_start
    call this%next()
    this%scf_count = this%scf_count + 1

    POP_SUB(propagator_rewind_scf_loop)
  end subroutine propagator_rewind_scf_loop

end module propagator_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
