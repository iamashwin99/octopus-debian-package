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

!> The basic elements defining algorithms
module algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use linked_list_oct_m
  use loct_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                  &
    algorithmic_operation_t, &
    algorithm_t,             &
    algorithm_iterator_t

  integer, parameter, public :: ALGO_LABEL_LEN = 50

  !> @brief Descriptor of one algorithmic operation.
  !!
  !! Algorithms are a sequence of operations. Those operations are identified by a
  !! string identifier that can be used in select case statements.
  type :: algorithmic_operation_t
    character(len=ALGO_LABEL_LEN) :: id    !< Operation identifier. We use a string instead of
    !!                                        an integer to minimize the chance of having duplicated identifiers.
    character(len=ALGO_LABEL_LEN) :: label !< Label describing what the code is doing when performing this operation.
  end type algorithmic_operation_t

  !# doc_start generic_algorithmic_operations
  ! Operations that can be used by any algorithm and, therefore, should be
  ! implemented by all systems.
  character(len=ALGO_LABEL_LEN), public, parameter ::   &
    SKIP                 = 'SKIP',                      &
    UPDATE_INTERACTIONS  = 'UPDATE_INTERACTIONS',       &
    STEP_DONE            = 'STEP_DONE',                 &
    REWIND_ALGORITHM     = 'REWIND_ALGORITHM'

  type(algorithmic_operation_t), public, parameter :: &
    OP_SKIP                 = algorithmic_operation_t(SKIP,                 'Skipping algorithmic step'), &
    OP_UPDATE_INTERACTIONS  = algorithmic_operation_t(UPDATE_INTERACTIONS,  'Algorithmic step - Updating interactions'), &
    OP_STEP_DONE            = algorithmic_operation_t(STEP_DONE,            'Propagation step finished'), &
    OP_REWIND_ALGORITHM     = algorithmic_operation_t(REWIND_ALGORITHM,     'Rewind algorithm')
  !# doc_end


  !> Iterator to loop over the algorithmic operations of an algorithm
  type, extends(linked_list_iterator_t) :: algorithm_iterator_t
    private
  contains
    procedure :: get_next => algorithm_iterator_get_next
  end type algorithm_iterator_t


  !> An algorithm is a list of algorithmic operations executed
  !! sequentially. This is implemented as a linked list of algorithmic
  !! operations.
  type, extends(linked_list_t), abstract :: algorithm_t
    private
    type(algorithm_iterator_t), public :: iter
    type(algorithmic_operation_t) :: current_ops

    type(algorithmic_operation_t), public       :: start_step
    type(algorithmic_operation_t), public       :: final_step

    integer, public :: algo_steps
    FLOAT, public   :: dt

    logical :: step_done
    logical, public :: inside_scf = .false.

    type(clock_t), public :: clock
    FLOAT :: start_time = M_ZERO
    FLOAT, public :: elapsed_time = M_ZERO

    logical, public :: is_static = .false.
  contains
    procedure :: add_operation => algorithm_add_operation
    procedure :: do_operation => algorithm_do_operation
    procedure :: update_elapsed_time => algorithm_update_elapsed_time
    procedure :: rewind => algorithm_rewind
    procedure :: next => algorithm_next
    procedure :: get_current_operation => algorithm_get_current_operation
    procedure(algorithm_finished), deferred :: finished
  end type algorithm_t

  abstract interface
    logical function algorithm_finished(this)
      import :: algorithm_t
      class(algorithm_t), intent(in) :: this
    end function algorithm_finished
  end interface

contains

  ! ---------------------------------------------------------
  subroutine algorithm_add_operation(this, operation)
    class(algorithm_t),            intent(inout) :: this
    type(algorithmic_operation_t), intent(in)    :: operation

    PUSH_SUB(algorithm_add_operation)

    call this%add_copy(operation)

    POP_SUB(algorithm_add_operation)
  end subroutine algorithm_add_operation

  ! ---------------------------------------------------------
  logical function algorithm_do_operation(this, operation) result(done)
    class(algorithm_t),            intent(inout) :: this
    type(algorithmic_operation_t), intent(in)    :: operation

    ! By default no algorithm specific operation is implemented in the algorithm
    ! class. Child classes that wish to change this behaviour should override
    ! this method.
    done = .false.

  end function algorithm_do_operation

  ! ---------------------------------------------------------
  subroutine algorithm_update_elapsed_time(this)
    class(algorithm_t), intent(inout) :: this

    PUSH_SUB(algorithm_update_elapsed_time)

    this%elapsed_time = loct_clock() - this%start_time

    POP_SUB(algorithm_update_elapsed_time)
  end subroutine algorithm_update_elapsed_time

  ! ---------------------------------------------------------
  subroutine algorithm_rewind(this)
    class(algorithm_t), intent(inout) :: this

    PUSH_SUB(algorithm_rewind)

    call this%iter%start(this)
    call this%next()
    this%start_time = loct_clock()

    POP_SUB(algorithm_rewind)
  end subroutine algorithm_rewind

  ! ---------------------------------------------------------
  subroutine algorithm_next(this)
    class(algorithm_t), intent(inout) :: this

    PUSH_SUB(algorithm_next)

    this%current_ops = this%iter%get_next()

    POP_SUB(algorithm_next)
  end subroutine algorithm_next

  ! ---------------------------------------------------------
  type(algorithmic_operation_t) function algorithm_get_current_operation(this) result(operation)
    class(algorithm_t), intent(in) :: this

    PUSH_SUB(algorithm_get_current_operation)

    operation = this%current_ops

    POP_SUB(algorithm_get_current_operation)
  end function algorithm_get_current_operation

  ! ---------------------------------------------------------
  function algorithm_iterator_get_next(this) result(operation)
    class(algorithm_iterator_t),  intent(inout) :: this
    type(algorithmic_operation_t)               :: operation

    PUSH_SUB(algorithm_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (algorithmic_operation_t)
      operation = ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(algorithm_iterator_get_next)
  end function algorithm_iterator_get_next

end module algorithm_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
