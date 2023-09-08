!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 M. Oliveira, K. Lively, A. Obzhirov, I. Albar
!! Copyright (C) 2022 S. Ohlmann
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

module box_cgal_oct_m
  use box_oct_m
  use box_shape_oct_m
  use cgal_polyhedra_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use string_oct_m
  use unit_oct_m

  implicit none

  private
  public :: box_cgal_t

  !> Class implementing a box defined by a file read using the cgal library
  type, extends(box_shape_t) :: box_cgal_t
    private
    type(cgal_polyhedra_t) :: cgal_poly
    character(len=1024) :: filename
  contains
    procedure :: shape_contains_points => box_cgal_shape_contains_points
    procedure :: write_info => box_cgal_write_info
    procedure :: short_info => box_cgal_short_info
    final     :: box_cgal_finalize
  end type box_cgal_t

  interface box_cgal_t
    procedure box_cgal_constructor
  end interface box_cgal_t

contains

  !--------------------------------------------------------------
  function box_cgal_constructor(dim, center, filename, length, namespace) result(box)
    integer,             intent(in) :: dim
    FLOAT,               intent(in) :: center(dim)
    character(len=*),    intent(in) :: filename
    FLOAT,               intent(in) :: length(dim) !< length of the parallelepiped along each corresponding basis vector
    type(namespace_t),   intent(in) :: namespace
    class(box_cgal_t), pointer :: box

    PUSH_SUB(box_cgal_constructor)

    if (dim /= 3) then
      message(1) = "CGAL boxes can only be used in 3 dimensions."
      call messages_fatal(1)
    end if

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, namespace, dim, center, bounding_box_min=-M_HALF*length, &
      bounding_box_max=M_HALF*length)
    ! initialize cgal part
    box%filename = trim(filename)
    call cgal_polyhedron_init(box%cgal_poly, trim(box%filename), verbose = .false.)

    box%bounding_box_l = M_HALF*length + abs(center)

    POP_SUB(box_cgal_constructor)
  end function box_cgal_constructor

  !--------------------------------------------------------------
  subroutine box_cgal_finalize(this)
    type(box_cgal_t), intent(inout) :: this

    PUSH_SUB(box_cgal_finalize)

    call cgal_polyhedron_end(this%cgal_poly)
    call box_shape_end(this)

    POP_SUB(box_cgal_finalize)
  end subroutine box_cgal_finalize

  !--------------------------------------------------------------
  function box_cgal_shape_contains_points(this, nn, xx) result(contained)
    class(box_cgal_t), intent(in)  :: this
    integer,           intent(in)  :: nn
    FLOAT,             intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip

    ! no push_sub/pop_sub, called too often

    ip = 1
    contained(:) = .false.

    do ip = 1, nn
      contained(ip) = cgal_polyhedron_point_inside(this%cgal_poly, xx(ip, 1), xx(ip, 2), xx(ip, 3)) &
        .neqv. this%is_inside_out()
    end do

  end function box_cgal_shape_contains_points

  !--------------------------------------------------------------
  subroutine box_cgal_write_info(this, iunit, namespace)
    class(box_cgal_t),           intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(box_cgal_write_info)

    write(message(1),'(2x,a)') 'Type = cgal'
    write(message(2),'(2x,a,x,a)') 'Loaded shape data from file', trim(this%filename)
    call messages_info(2, iunit=iunit, namespace=namespace)

    POP_SUB(box_cgal_write_info)
  end subroutine box_cgal_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_cgal_short_info(this, unit_length) result(info)
    class(box_cgal_t), intent(in) :: this
    type(unit_t),      intent(in) :: unit_length

    PUSH_SUB(box_cgal_short_info)

    write(info,'(3a)') 'BoxShape = cgal; BoxCgalFile = ', trim(this%filename)

    POP_SUB(box_cgal_short_info)
  end function box_cgal_short_info

end module box_cgal_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
