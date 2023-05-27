!! Copyright (C) 2021 M. Oliveira, K. Lively, A. Obzhirov, I. Albar
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

module multibox_oct_m
  use basis_vectors_oct_m
  use box_oct_m
  use debug_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::      &
    multibox_t,  &
    multibox_end

  !> Abstract class for boxes that are made up of a list of boxes.
  type, abstract, extends(box_t) :: multibox_t
    private
    type(box_list_t), public :: list !< list containing the boxes that make up this multibox
  contains
    procedure :: bounds => multibox_bounds
    procedure :: add_box => multibox_add_box
  end type multibox_t

contains

  !--------------------------------------------------------------
  subroutine multibox_end(this)
    class(multibox_t), intent(inout) :: this

    type(box_iterator_t)  :: iter
    class(box_t), pointer :: box

    PUSH_SUB(multibox_end)

    SAFE_DEALLOCATE_A(this%bounding_box_l)

    call iter%start(this%list)
    do while (iter%has_next())
      box => iter%get_next()
      SAFE_DEALLOCATE_P(box)
    end do

    POP_SUB(multibox_end)
  end subroutine multibox_end

  !--------------------------------------------------------------
  ! Although formally this will return the bounds of a union of boxes, we will
  ! also use it for an intersection of boxes. This should be a reasonable (and
  ! much simpler to implement) approximation to the bounds of an intersection of
  ! boxes, as the intersection should be contained in the union (unless the box
  ! is unbound, but in that case it does not matter how we handle this).
  function multibox_bounds(this, axes) result(bounds)
    class(multibox_t),                intent(in)  :: this
    class(basis_vectors_t), optional, intent(in)  :: axes
    FLOAT :: bounds(2, this%dim)

    FLOAT :: box_bounds(2, this%dim)
    type(box_iterator_t) :: iter
    class(box_t), pointer :: box

    PUSH_SUB(multibox_bounds)

    bounds(1, :) =  M_HUGE
    bounds(2, :) = -M_HUGE
    call iter%start(this%list)
    do while (iter%has_next())
      box => iter%get_next()
      if (this%is_inside_out()) then
        call box%turn_inside_out()
      end if
      box_bounds = box%bounds(axes)

      where (box_bounds(1,:) < bounds(1,:))
        bounds(1, :) = box_bounds(1,:)
      elsewhere (box_bounds(2,:) > bounds(2,:))
        bounds(2, :) = box_bounds(2,:)
      end where
      if (this%is_inside_out()) then
        call box%turn_inside_out()
      end if
    end do

    POP_SUB(multibox_bounds)
  end function multibox_bounds

  !--------------------------------------------------------------
  subroutine multibox_add_box(this, new_box)
    class(multibox_t), intent(inout) :: this
    class(box_t),      intent(in)    :: new_box

    integer :: idir

    PUSH_SUB(multibox_add_box)

    ASSERT(this%dim == new_box%dim)

    do idir = 1, this%dim
      if (new_box%bounding_box_l(idir) > this%bounding_box_l(idir)) then
        this%bounding_box_l(idir) = new_box%bounding_box_l(idir)
      end if
    end do

    call this%list%add(new_box)

    POP_SUB(multibox_add_box)
  end subroutine multibox_add_box

end module multibox_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
