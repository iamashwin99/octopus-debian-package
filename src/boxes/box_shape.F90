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

module box_shape_oct_m
  use box_oct_m
  use basis_vectors_oct_m
  use debug_oct_m
  use global_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::         &
    box_shape_t,    &
    box_shape_init, &
    box_shape_end

  !> Base class for more specialized boxes that are defined by a shape and have
  !! a center and basis vectors.
  type, abstract, extends(box_t) :: box_shape_t
    private
    FLOAT, allocatable,    public :: center(:)         !< where is the box centered
    type(basis_vectors_t), public :: axes              !< Unit vectors used to define the orientation of the shape in space
    FLOAT, allocatable            :: bounding_box(:,:) !< Start and end points of the axes-aligned bounding box.
    !!                                                    The axes are the one that define the orientation of the
    !!                                                    shape in space (see above).
  contains
    procedure :: contains_points => box_shape_contains_points
    procedure :: bounds => box_shape_bounds
    procedure(box_shape_shape_contains_points), deferred :: shape_contains_points
  end type box_shape_t

  abstract interface
    function box_shape_shape_contains_points(this, nn, xx) result(contained)
      import :: box_shape_t
      class(box_shape_t), intent(in) :: this
      integer,            intent(in) :: nn
      FLOAT,              intent(in) :: xx(:,:)
      logical :: contained(1:nn)
    end function box_shape_shape_contains_points
  end interface

contains

  !--------------------------------------------------------------
  recursive function box_shape_contains_points(this, nn, xx) result(contained)
    class(box_shape_t), intent(in)  :: this
    integer,            intent(in)  :: nn
    FLOAT,              intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    ! no push_sub because this function is called very frequently

    contained = this%shape_contains_points(nn, this%axes%from_cartesian(nn, xx))

  end function box_shape_contains_points

  !--------------------------------------------------------------
  subroutine box_shape_init(this, namespace, dim, center, bounding_box_min, bounding_box_max, axes)
    class(box_shape_t), intent(inout) :: this
    type(namespace_t),  intent(in)    :: namespace
    integer,            intent(in)    :: dim
    FLOAT,              intent(in)    :: center(dim)
    FLOAT,              intent(in)    :: bounding_box_min(dim)
    FLOAT,              intent(in)    :: bounding_box_max(dim)
    FLOAT, optional,    intent(in)    :: axes(dim, dim)

    integer :: idir
    FLOAT :: normalized_axes(dim, dim)

    this%dim = dim

    SAFE_ALLOCATE(this%center(1:dim))
    SAFE_ALLOCATE(this%bounding_box_l(1:dim))
    this%bounding_box_l = M_ZERO

    this%center(1:dim) = center(1:dim)

    if (present(axes)) then
      do idir = 1, dim
        normalized_axes(:, idir) = axes(:, idir)/norm2(axes(:, idir))
      end do
    else
      normalized_axes = diagonal_matrix(dim, M_ONE)
    end if
    this%axes = basis_vectors_t(namespace, dim, normalized_axes)

    SAFE_ALLOCATE(this%bounding_box(1:2, 1:dim))
    this%bounding_box(1,:) = bounding_box_min + this%center
    this%bounding_box(2,:) = bounding_box_max + this%center

  end subroutine box_shape_init

  !--------------------------------------------------------------
  !> Returns the bounding box of the shape. This is a bounding box aligned with
  !! some given axes (if present) or with the Cartesian axes.
  !! Note that in both cases we are returning the bounding box of the bounding
  !! box aligned with the axes that define the orientation of the shape in
  !! space. This is not optimal, but should be acceptable in most cases and in
  !! any case it is always possible to override this method in the children
  !! class.
  function box_shape_bounds(this, axes) result(bounds)
    class(box_shape_t),               intent(in)  :: this
    class(basis_vectors_t), optional, intent(in)  :: axes
    FLOAT :: bounds(2, this%dim)

    integer :: idir, ib
    FLOAT :: vertex_red(this%dim), vertex_cart(this%dim)

    PUSH_SUB(box_shape_bounds)

    if (this%is_inside_out()) then
      ! The box is unbound.
      bounds(1, :) = -M_HUGE
      bounds(2, :) =  M_HUGE
    else
      ! The bounds must be located at the axis-aligned bounding box vertices, so
      ! we loop over the vertices, convert them to reduced coordinates, and find
      ! the maximum and minimum coordinates. This algorithm might not be
      ! optimal, but it will do for now.
      bounds(1, :) =  M_HUGE
      bounds(2, :) = -M_HUGE
      do idir = 1, this%dim
        do ib = 1, 2
          ! Vertex in box reduced coordinates
          vertex_red = M_ZERO
          vertex_red(idir) = this%bounding_box(ib, idir)

          ! Transform vertex to Cartesian coordinates
          vertex_cart = this%axes%to_cartesian(vertex_red)

          ! Convert to reduced coordinates of the input axes
          if (present(axes)) then
            vertex_red = axes%from_cartesian(vertex_cart)
          else
            vertex_red = vertex_cart
          end if

          where (vertex_red < bounds(1,:))
            bounds(1, :) = vertex_red
          elsewhere (vertex_red > bounds(2,:))
            bounds(2, :) = vertex_red
          end where
        end do
      end do
    end if

    POP_SUB(box_shape_bounds)
  end function box_shape_bounds

  !--------------------------------------------------------------
  subroutine box_shape_end(this)
    class(box_shape_t), intent(inout) :: this

    SAFE_DEALLOCATE_A(this%center)
    SAFE_DEALLOCATE_A(this%bounding_box)
    SAFE_DEALLOCATE_A(this%bounding_box_l)

  end subroutine box_shape_end

end module box_shape_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
