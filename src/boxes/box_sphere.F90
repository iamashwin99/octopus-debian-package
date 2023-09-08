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

module box_sphere_oct_m
  use box_oct_m
  use box_shape_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public :: box_sphere_t

  !> Class implementing a spherical box.
  type, extends(box_shape_t) :: box_sphere_t
    FLOAT   :: radius !< the radius of the sphere
  contains
    procedure :: shape_contains_points => box_sphere_shape_contains_points
    procedure :: get_surface_points => box_sphere_shape_get_surface_points
    procedure :: get_surface_point_info => box_sphere_shape_get_surface_point_info
    procedure :: write_info => box_sphere_write_info
    procedure :: short_info => box_sphere_short_info
    final     :: box_sphere_finalize
  end type box_sphere_t

  interface box_sphere_t
    procedure box_sphere_constructor
  end interface box_sphere_t

contains

  !--------------------------------------------------------------
  function box_sphere_constructor(dim, center, radius, namespace) result(box)
    integer,            intent(in) :: dim
    FLOAT,              intent(in) :: center(dim)
    FLOAT,              intent(in) :: radius
    type(namespace_t),  intent(in) :: namespace
    class(box_sphere_t), pointer :: box

    integer :: idir

    PUSH_SUB(box_sphere_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, namespace, dim, center, bounding_box_min=[(-radius, idir=1, dim)], &
      bounding_box_max=[(radius, idir=1, dim)])
    box%radius = radius
    box%bounding_box_l = abs(center) + radius

    POP_SUB(box_sphere_constructor)
  end function box_sphere_constructor

  !--------------------------------------------------------------
  subroutine box_sphere_finalize(this)
    type(box_sphere_t), intent(inout) :: this

    PUSH_SUB(box_sphere_finalize)

    call box_shape_end(this)

    POP_SUB(box_sphere_finalize)
  end subroutine box_sphere_finalize

  !--------------------------------------------------------------
  function box_sphere_shape_contains_points(this, nn, xx) result(contained)
    class(box_sphere_t), intent(in)  :: this
    integer,             intent(in)  :: nn
    FLOAT,               intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip

    do ip = 1, nn
      contained(ip) = sum((xx(ip, 1:this%dim) - this%center(1:this%dim))**2) <= (this%radius + BOX_BOUNDARY_DELTA)**2 &
        .neqv. this%is_inside_out()
    end do

  end function box_sphere_shape_contains_points

  !--------------------------------------------------------------
  !> Get a mask for the grid points telling which of them are surface points
  !! 1. Create a box that is fractionally smaller than the original.
  !! 2. Look for points that belong to the original box but not to the smaller one - These define surface points
  function box_sphere_shape_get_surface_points(this, namespace, mesh_spacing, nn, xx, number_of_layers) result(surface_points)
    class(box_sphere_t), intent(in)  :: this
    type(namespace_t),   intent(in)  :: namespace
    FLOAT,               intent(in)  :: mesh_spacing(:)
    integer,             intent(in)  :: nn
    FLOAT,               intent(in)  :: xx(:,:)
    integer, optional,   intent(in)  :: number_of_layers

    logical          :: surface_points(1:nn)
    FLOAT            :: shrink
    FLOAT, parameter :: compression = CNST(0.94)
    class(box_sphere_t), pointer :: shrinked_box
    ! Check that the spacing is the same in each direction
    ASSERT(abs(mesh_spacing(1) - mesh_spacing(2)) < M_EPSILON .and. abs(mesh_spacing(1) - mesh_spacing(3)) < M_EPSILON)

    ! First of all determine how much the shrink should be based on the spacing
    shrink = 1 - compression * (mesh_spacing(1) / this%radius)

    ! Then, we have to initialize the new sphere to be slightly smaller that the original one
    shrinked_box => box_sphere_t(this%dim, this%center, this%radius * shrink, namespace)

    ! Then we look for points of the old sphere that are not containted in the smaller sphere
    ! box_sphere_shape_contains_points will return the point contained in the shrinked box, but we are interested in the ones not
    ! contained, as they will the surface points of the original sphere
    surface_points = .not. box_sphere_shape_contains_points(shrinked_box, nn, xx)

    call box_shape_end(shrinked_box)
    SAFE_DEALLOCATE_P(shrinked_box)

  end function box_sphere_shape_get_surface_points

  !--------------------------------------------------------------
  subroutine box_sphere_shape_get_surface_point_info(this, point_coordinates, mesh_spacing, normal_vector, surface_element)
    class(box_sphere_t), intent(in)  :: this
    FLOAT,               intent(in)  :: point_coordinates(:) !< (x,y,z) coordinates of the point
    FLOAT,               intent(in)  :: mesh_spacing(:)      !< spacing of the mesh
    FLOAT,               intent(out) :: normal_vector(:)     !< normal vector to the surface point
    FLOAT,               intent(out) :: surface_element      !< surface element (needed to compute the surface integral)

    FLOAT :: dtheta, dphi, rr

    PUSH_SUB(box_sphere_shape_get_surface_point_info)
    dtheta = M_ZERO
    dphi = M_ZERO

    ! Compute the normal vector to the surface point

    ! For a sphere, the normal vector to a point is simply the vector representing its radius, but normalized.
    ! Therefore it is enough to normalize the vector containing the coordinates of the vector

    ! So first of all we have to determine in which face we are. To do that let us divide the vector containing the coordinates
    ! of the point by the maximum length of the parallelepiped
    normal_vector(:) = point_coordinates(:) - this%center(:)

    rr = sqrt(normal_vector(1)**2 + normal_vector(2)**2)
    if (rr > M_ZERO) then
      dtheta = normal_vector(1) * normal_vector(3) / (rr * this%radius**2) * mesh_spacing(1) + &
        normal_vector(2) * normal_vector(3) / (rr * this%radius**2) * mesh_spacing(2) - rr / this%radius**2 * mesh_spacing(3)

      dphi = (- normal_vector(2) * mesh_spacing(1) + normal_vector(1) * mesh_spacing(2)) / rr**2
    end if

    normal_vector(:) = normal_vector(:) / norm2(normal_vector)

    ! Compute the surface element, which for a sphere is r^2 * sin(theta) dxdy, where theta = acos(z / r)
    surface_element = abs(dtheta * dphi) * this%radius**2 * sin(acos(normal_vector(3)))

    POP_SUB(box_sphere_shape_get_surface_point_info)
  end subroutine box_sphere_shape_get_surface_point_info

  !--------------------------------------------------------------
  subroutine box_sphere_write_info(this, iunit, namespace)
    class(box_sphere_t),         intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(box_sphere_write_info)

    write(message(1),'(2x,a)') 'Type = sphere'
    write(message(2),'(2x,3a,f7.3)') 'Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%radius)
    call messages_info(2, iunit=iunit, namespace=namespace)

    POP_SUB(box_sphere_write_info)
  end subroutine box_sphere_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_sphere_short_info(this, unit_length) result(info)
    class(box_sphere_t), intent(in) :: this
    type(unit_t),        intent(in) :: unit_length

    PUSH_SUB(box_sphere_short_info)

    write(info,'(a,f11.6,a,a)') 'BoxShape = sphere; Radius =', units_from_atomic(unit_length, this%radius), ' ', &
      trim(units_abbrev(unit_length))

    POP_SUB(box_sphere_short_info)
  end function box_sphere_short_info

end module box_sphere_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
