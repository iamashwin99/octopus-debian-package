!! Copyright (C) 2008 X. Andrade
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

module stencil_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                        &
    stencil_t,                     &
    stencil_allocate,              &
    stencil_copy,                  &
    stencil_end,                   &
    stencil_init_center,           &
    stencil_union

  type stargeneral_arms_t
    ! Components are public by default
    integer          :: arms(1:3,1:3) = reshape([0,0,0,0,0,0,0,0,0], (/3, 3/))
    integer          :: narms = 0
  end type stargeneral_arms_t


  type stencil_t
    ! Components are public by default
    integer              :: center
    integer              :: size
    integer              :: npoly
    integer, allocatable :: points(:, :)

    ! The stargeneral arms
    type(stargeneral_arms_t) :: stargeneral
  end type stencil_t

contains

  !-------------------------------------------------------
  subroutine stencil_allocate(this, size)
    type(stencil_t), intent(inout) :: this
    integer,         intent(in)  :: size

    PUSH_SUB(stencil_allocate)

    this%size = size
    this%npoly = size

    SAFE_ALLOCATE(this%points(1:MAX_DIM, 1:size))

    this%points = 0

    POP_SUB(stencil_allocate)
  end subroutine stencil_allocate

  !-------------------------------------------------------
  subroutine stencil_copy(input, output)
    type(stencil_t), intent(in)  :: input
    type(stencil_t), intent(out) :: output

    PUSH_SUB(stencil_copy)

    call stencil_allocate(output, input%size)
    output%points(1:MAX_DIM, 1:output%size) = input%points(1:MAX_DIM, 1:output%size)
    output%center = input%center
    output%npoly = input%npoly

    output%stargeneral%narms = input%stargeneral%narms
    output%stargeneral%arms = input%stargeneral%arms

    POP_SUB(stencil_copy)
  end subroutine stencil_copy


  !-------------------------------------------------------
  subroutine stencil_end(this)
    type(stencil_t), intent(inout) :: this

    PUSH_SUB(stencil_end)

    SAFE_DEALLOCATE_A(this%points)

    POP_SUB(stencil_end)
  end subroutine stencil_end


  !-------------------------------------------------------
  subroutine stencil_union(dim, st1, st2, stu)
    integer,         intent(in)    :: dim
    type(stencil_t), intent(inout) :: st1
    type(stencil_t), intent(inout) :: st2
    type(stencil_t), intent(inout) :: stu

    integer :: idir, ii, jj, nstu
    logical :: not_in_st1

    PUSH_SUB(stencil_union)

    call stencil_allocate(stu, st1%size + st2%size)

    ! copy the first stencil
    do ii = 1, st1%size
      do idir = 1, dim
        stu%points(idir, ii) = st1%points(idir, ii)
      end do
    end do

    nstu = st1%size

    do ii = 1, st2%size

      not_in_st1 = .true.

      ! check whether that point was already in the stencil
      do jj = 1, st1%size
        if (all(st1%points(1:dim, jj) == st2%points(1:dim, ii))) then
          not_in_st1 = .false.
          exit
        end if
      end do

      if (not_in_st1) then !add it
        nstu = nstu + 1
        stu%points(1:dim, nstu) = st2%points(1:dim, ii)
      end if

    end do

    stu%points(dim + 1:MAX_DIM, 1:nstu) = 0

    stu%size = nstu

    ! this is not defined for a union, which could be any combination. The
    ! weights should already be known here
    stu%npoly = -1

    call stencil_init_center(stu)

    POP_SUB(stencil_union)
  end subroutine stencil_union


  !-------------------------------------------------------
  subroutine stencil_init_center(this)
    type(stencil_t), intent(inout) :: this

    integer :: ii

    PUSH_SUB(stencil_init_center)

    this%center = -1

    do ii = 1, this%size
      if (all(this%points(1:MAX_DIM, ii) == 0)) this%center = ii
    end do

    POP_SUB(stencil_init_center)
  end subroutine stencil_init_center

end module stencil_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
