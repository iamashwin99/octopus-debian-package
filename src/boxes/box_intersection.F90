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

module box_intersection_oct_m
  use box_oct_m
  use debug_oct_m
  use multibox_oct_m
  use namespace_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m

  implicit none

  private
  public ::           &
    box_intersection_t

  !> Class implementing a box that is an intersection other boxes.
  !! Note that the we do not override the bounds method of multibox, although
  !! formally it returns the bounds of an union of boxes. The reason is that
  !! the proper calculation of the bounds of an intersection of boxes would be
  !! complicated and would require detailed knowleged of each box shape. Since
  !! the main purpose of the function is to create bounding boxes, using the
  !! bounds of the union of boxes is a reasonable approximation.
  type, extends(multibox_t) :: box_intersection_t
    private
  contains
    procedure :: contains_points => box_intersection_contains_points
    procedure :: write_info => box_intersection_write_info
    procedure :: short_info => box_intersection_short_info
    final     :: box_intersection_finalize
  end type box_intersection_t

  interface box_intersection_t
    procedure box_intersection_constructor
  end interface box_intersection_t

contains

  !--------------------------------------------------------------
  function box_intersection_constructor(dim) result(box)
    integer, intent(in) :: dim
    class(box_intersection_t), pointer :: box

    PUSH_SUB(box_intersection_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)
    SAFE_ALLOCATE(box%bounding_box_l(1:dim))

    ! Initialize box
    box%dim = dim
    box%bounding_box_l = M_ZERO

    POP_SUB(box_intersection_constructor)
  end function box_intersection_constructor

  !--------------------------------------------------------------
  subroutine box_intersection_finalize(this)
    type(box_intersection_t), intent(inout) :: this

    PUSH_SUB(box_intersection_finalize)

    call multibox_end(this)

    POP_SUB(box_intersection_finalize)
  end subroutine box_intersection_finalize

  !--------------------------------------------------------------
  recursive function box_intersection_contains_points(this, nn, xx) result(contained)
    class(box_intersection_t), intent(in) :: this
    integer,           intent(in) :: nn
    FLOAT,             intent(in) :: xx(:,:)
    logical :: contained(nn)

    integer :: ip
    FLOAT :: point(1:this%dim)
    type(box_iterator_t) :: iter
    class(box_t), pointer :: box

    ! A point must be inside all boxes to be considered inside an intersection of boxes
    do ip = 1, nn
      point(1:this%dim) = xx(ip, 1:this%dim)
      contained(ip) = .true.

      call iter%start(this%list)
      do while (iter%has_next())
        box => iter%get_next()
        contained(ip) = box%contains_point(point)
        if (.not. contained(ip)) exit
      end do

      contained(ip) = contained(ip) .neqv. this%is_inside_out()
    end do

  end function box_intersection_contains_points

  !--------------------------------------------------------------
  subroutine box_intersection_write_info(this, iunit, namespace)
    class(box_intersection_t),   intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(box_intersection_write_info)

    ! Todo: need to decide how best to display the information of the boxes that make the intersection

    POP_SUB(box_intersection_write_info)
  end subroutine box_intersection_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_intersection_short_info(this, unit_length) result(info)
    class(box_intersection_t), intent(in) :: this
    type(unit_t),              intent(in) :: unit_length

    PUSH_SUB(box_intersection_short_info)

    ! Todo: need to decide how best to display the information of the boxes that make the intersection
    info = ''

    POP_SUB(box_intersection_short_info)
  end function box_intersection_short_info

end module box_intersection_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
