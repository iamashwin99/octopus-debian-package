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

module box_oct_m
  use basis_vectors_oct_m
  use debug_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m

  implicit none

  private
  public ::        &
    box_t,         &
    box_list_t,    &
    box_iterator_t

  integer, parameter, public :: BOX_INFO_LEN=200
  FLOAT, parameter, public :: BOX_BOUNDARY_DELTA = CNST(1e-12)

  !> @brief class to tell whether a point is inside or outside
  !!
  !! The purpose of a box is to tell if something is inside or outside of it.
  !! To do that it provides a function that tells if a given list of points are
  !! inside or outside the box. Furthermore, a box might be turned inside out,
  !! i.e., in that case what is usually considered inside becomes outside and
  !! vice-versa.
  type, abstract :: box_t
    private
    integer, public :: dim                          !< dimensions of the space the box lives in
    logical :: inside_out = .false.                 !< if the box is inside out or not
    FLOAT, allocatable, public :: bounding_box_l(:) !< Half lengths of the bounding box that contains the
    !!                                                 box. Note that this box always contains the origin
    !!                                                 and is symmetrical with respect to it.
  contains
    procedure(box_contains_points),        deferred :: contains_points
    procedure(box_bounds),                 deferred :: bounds
    procedure(box_write_info),             deferred :: write_info
    procedure(box_short_info),             deferred :: short_info
    procedure :: get_surface_points => box_get_surface_points
    procedure :: get_surface_point_info => box_get_surface_point_info
    procedure, non_overridable :: contains_point => box_contains_point
    procedure, non_overridable :: is_inside_out => box_is_inside_out
    procedure, non_overridable :: turn_inside_out => box_turn_inside_out
  end type box_t

  abstract interface
    !> Given a list of points, this function should return an array indicating
    !! for each point if it is inside the box or not.
    recursive function box_contains_points(this, nn, xx) result(contained)
      import :: box_t
      class(box_t), intent(in) :: this
      integer,      intent(in) :: nn      !< number of points to check
      FLOAT,        intent(in) :: xx(:,:) !< points to check. The sizes are
      !!                                     (1:,1:this%dim), so that it is
      !!                                     possible to pass an array with more
      !!                                     points than the ones we are
      !!                                     checking.
      logical :: contained(1:nn)
    end function box_contains_points

    !> Box bounds along some axes
    function box_bounds(this, axes) result(bounds)
      import :: box_t
      import :: basis_vectors_t
      class(box_t),                     intent(in)  :: this
      class(basis_vectors_t), optional, intent(in)  :: axes
      FLOAT                                         :: bounds(2, this%dim) !< minimum and maximum coordinates along each axis.
    end function box_bounds

    !> Write the complete information about the box to a file.
    subroutine box_write_info(this, iunit, namespace)
      import :: box_t
      import :: namespace_t
      class(box_t),                intent(in) :: this
      integer,           optional, intent(in) :: iunit
      type(namespace_t), optional, intent(in) :: namespace
    end subroutine box_write_info

    !> Return a string containing a short description of the box.
    function box_short_info(this, unit_length)
      use unit_oct_m
      import :: box_t
      import :: BOX_INFO_LEN
      class(box_t), intent(in) :: this
      type(unit_t), intent(in) :: unit_length
      character(len=BOX_INFO_LEN) :: box_short_info
    end function box_short_info
  end interface

  !> These classes extends the list and list iterator to create a box list.
  type, extends(linked_list_t) :: box_list_t
    private
  contains
    procedure :: add => box_list_add_node
  end type box_list_t

  type, extends(linked_list_iterator_t) :: box_iterator_t
    private
  contains
    procedure :: get_next => box_iterator_get_next
  end type box_iterator_t

contains

  !!--------------------------------------------------------------
  !> Turn a box inside out.
  subroutine box_turn_inside_out(this)
    class(box_t), intent(inout) :: this

    this%inside_out = .not. this%inside_out

  end subroutine box_turn_inside_out

  !!--------------------------------------------------------------
  !> Is the box inside out?
  logical function box_is_inside_out(this)
    class(box_t), intent(in) :: this

    box_is_inside_out = this%inside_out

  end function box_is_inside_out

  !!---------------------------------------------------------------
  !> Convenience function to check if a single point is inside the box when that
  !! point is passed as a rank-one array.
  recursive logical function box_contains_point(this, xx) result(contained)
    class(box_t),         intent(in) :: this
    FLOAT,        target, intent(in) :: xx(1:this%dim)

    FLOAT, pointer :: xx_ptr(:,:)
    logical :: points_contained(1)

    xx_ptr(1:1, 1:this%dim) => xx(1:this%dim)
    points_contained = this%contains_points(1, xx_ptr)
    contained = points_contained(1)

  end function box_contains_point

  !--------------------------------------------------------------
  function box_get_surface_points(this, namespace, mesh_spacing, nn, xx, number_of_layers) result(surface_points)
    class(box_t),         intent(in)  :: this
    type(namespace_t),    intent(in)  :: namespace
    FLOAT,                intent(in)  :: mesh_spacing(:)
    integer,              intent(in)  :: nn
    FLOAT,                intent(in)  :: xx(:,:)
    integer, optional,    intent(in)  :: number_of_layers
    logical :: surface_points(1:nn)

    surface_points = .false.
    call messages_not_implemented("get_surface_points for box shape")

  end function box_get_surface_points

  !--------------------------------------------------------------
  subroutine box_get_surface_point_info(this, point_coordinates, mesh_spacing, normal_vector, surface_element)
    class(box_t), intent(in)  :: this
    FLOAT,        intent(in)  :: point_coordinates(:) !< (x,y,z) coordinates of the point
    FLOAT,        intent(in)  :: mesh_spacing(:)      !< spacing of the mesh
    FLOAT,        intent(out) :: normal_vector(:)     !< normal vector to the surface point
    FLOAT,        intent(out) :: surface_element      !< surface element (needed to compute the surface integral)

    PUSH_SUB(box_get_surface_point_info)

    call messages_not_implemented("get_surface_point_info for box shape")

    POP_SUB(box_get_surface_point_info)
  end subroutine box_get_surface_point_info

  ! ---------------------------------------------------------
  subroutine box_list_add_node(this, box)
    class(box_list_t)    :: this
    class(box_t), target :: box

    select type (box)
    class is (box_t)
      call this%add_ptr(box)
    class default
      ASSERT(.false.)
    end select

  end subroutine box_list_add_node

  ! ---------------------------------------------------------
  function box_iterator_get_next(this) result(box)
    class(box_iterator_t), intent(inout) :: this
    class(box_t),          pointer       :: box

    select type (ptr => this%get_next_ptr())
    class is (box_t)
      box => ptr
    class default
      ASSERT(.false.)
    end select

  end function box_iterator_get_next

end module box_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
