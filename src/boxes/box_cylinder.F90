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

module box_cylinder_oct_m
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
  public :: box_cylinder_t

  !> Class implementing a cylinder box. The cylinder axis is always along the
  !> first direction defined by the box_shape_t basis vectors.
  type, extends(box_shape_t) :: box_cylinder_t
    private
    FLOAT, public :: radius      !< the radius of the cylinder
    FLOAT, public :: half_length !< half the length of the cylinder

    logical :: periodic_boundaries = .false. !< are the bases of the cylinder to be treated as periodic?
  contains
    procedure :: shape_contains_points => box_cylinder_shape_contains_points
    procedure :: get_surface_points => box_cylinder_shape_get_surface_points
    procedure :: get_surface_point_info => box_cylinder_shape_get_surface_point_info
    procedure :: write_info => box_cylinder_write_info
    procedure :: short_info => box_cylinder_short_info
    final     :: box_cylinder_finalize
  end type box_cylinder_t

  interface box_cylinder_t
    procedure box_cylinder_constructor
  end interface box_cylinder_t

contains

  !--------------------------------------------------------------
  function box_cylinder_constructor(dim, center, axes, radius, length, namespace, periodic_boundaries) result(box)
    integer,            intent(in) :: dim
    FLOAT,              intent(in) :: center(dim)
    FLOAT,              intent(in) :: axes(dim, dim)
    FLOAT,              intent(in) :: radius                  !< cylinder radius
    FLOAT,              intent(in) :: length                  !< lenght of the cylinder along the basis vectors
    type(namespace_t),  intent(in) :: namespace
    logical, optional,  intent(in) :: periodic_boundaries     !< are the bases of the cylinder to be treated as periodic?
    class(box_cylinder_t), pointer :: box

    integer :: idir

    PUSH_SUB(box_cylinder_constructor)

    ! Sanity checks
    if (dim <= 2) then
      message(1) = "Cannot create a cylinder in 1D or 2D. Use sphere if you want a circle."
      call messages_fatal(1, namespace=namespace)
    end if

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, namespace, dim, center, bounding_box_min=[-M_HALF*length, (-radius, idir=2,dim)], &
      bounding_box_max=[M_HALF*length, (radius, idir=2,dim)], axes=axes)
    box%radius = radius
    box%half_length = M_HALF*length

    if (present(periodic_boundaries)) then
      box%periodic_boundaries = periodic_boundaries
    end if

    box%bounding_box_l(1) = M_HALF*length + abs(center(1))  !< xsize
    box%bounding_box_l(2:dim) = radius + abs(center(2:dim)) !< rsize

    POP_SUB(box_cylinder_constructor)
  end function box_cylinder_constructor

  !--------------------------------------------------------------
  subroutine box_cylinder_finalize(this)
    type(box_cylinder_t), intent(inout) :: this

    PUSH_SUB(box_cylinder_finalize)

    call box_shape_end(this)

    POP_SUB(box_cylinder_finalize)
  end subroutine box_cylinder_finalize

  !--------------------------------------------------------------
  recursive function box_cylinder_shape_contains_points(this, nn, xx) result(contained)
    class(box_cylinder_t), intent(in)  :: this
    integer,               intent(in)  :: nn
    FLOAT,                 intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip
    FLOAT :: rr2, vv(this%dim)

    do ip = 1, nn
      vv = xx(ip, 1:this%dim) - this%center(1:this%dim)

      ! First check if we are "inside" along the axis direction. If not, do not bother checking the other directions.
      if (.not. this%periodic_boundaries) then
        contained(ip) = abs(vv(1)) <= this%half_length + BOX_BOUNDARY_DELTA .neqv. this%is_inside_out()
      else
        ! When periodic, we exclude one of the faces from the box. Also, we
        ! never consider the box to be inside out along the periodic dimension.
        contained(ip) = abs(vv(1) + BOX_BOUNDARY_DELTA)  <= this%half_length
      end if
      if (.not. contained(ip)) cycle

      ! Check if we are inside along the directions perpendicular to the axis
      rr2 = sum(vv(2:this%dim)**2)
      contained(ip) = rr2 <= (this%radius + BOX_BOUNDARY_DELTA)**2 .neqv. this%is_inside_out()
    end do

  end function box_cylinder_shape_contains_points

  !--------------------------------------------------------------
  !> Get a mask for the grid points telling which of them are surface points
  !! 1. Create a box that is fractionally smaller than the original.
  !! 2. Look for points that belong to the original box but not to the smaller one - These define surface points
  function box_cylinder_shape_get_surface_points(this, namespace, mesh_spacing, nn, xx, number_of_layers) result(surface_points)
    class(box_cylinder_t), intent(in)  :: this
    type(namespace_t),     intent(in)  :: namespace
    FLOAT,                 intent(in)  :: mesh_spacing(:)
    integer,               intent(in)  :: nn
    FLOAT,                 intent(in)  :: xx(:,:)
    integer, optional,     intent(in)  :: number_of_layers

    logical          :: surface_points(1:nn)
    FLOAT            :: shrink(this%dim)
    FLOAT, parameter :: compression = CNST(0.94)
    class(box_cylinder_t), pointer :: shrinked_box
    ! TODO: See Issue 705 (ftroisi) - Helmholtz decomposition refinement
    ! The axis of the cylinder is always in the first direction of the box_shape. So the following assert checks that the spacing
    ! in the other two direction (which define the circle) is the same (to avoid problems due to deformation effects)
    ASSERT(abs(mesh_spacing(2) - mesh_spacing(3)) < M_EPSILON)

    ! First of all determine how much the shrink should be based on the spacing in each direction
    shrink(1) = 1 - (1 - BOX_BOUNDARY_DELTA) * (mesh_spacing(1) / this%half_length)
    shrink(2) = 1 - compression * (mesh_spacing(2) / this%radius)

    ! First of all we create a shrunk cylindric box
    shrinked_box => box_cylinder_t(this%dim, this%center, this%axes%vectors, this%radius * shrink(2), &
      M_TWO * this%half_length * shrink(1), namespace, periodic_boundaries=this%periodic_boundaries)

    ! Then we look for points of the old box that are not containted in the smaller box
    ! box_parallelepiped_shape_contains_points will return the point contained in the shrinked box,
    ! but we are interested in the ones not contained, as they will the surface points of the original parallelepiped
    surface_points = .not. box_cylinder_shape_contains_points(shrinked_box, nn, xx)
    SAFE_DEALLOCATE_P(shrinked_box)

  end function box_cylinder_shape_get_surface_points

  !--------------------------------------------------------------
  subroutine box_cylinder_shape_get_surface_point_info(this, point_coordinates, mesh_spacing, normal_vector, surface_element)
    class(box_cylinder_t), intent(in)  :: this
    FLOAT,                 intent(in)  :: point_coordinates(:) !< (x,y,z) coordinates of the point
    FLOAT,                 intent(in)  :: mesh_spacing(:)      !< spacing of the mesh
    FLOAT,                 intent(out) :: normal_vector(:)     !< normal vector to the surface point
    FLOAT,                 intent(out) :: surface_element      !< surface element (needed to compute the surface integral)

    PUSH_SUB(box_cylinder_shape_get_surface_point_info)

    ! Compute the normal vector to the surface point

    ! For a cylinder, the normal vector to a point is the normal to a circle if we are on its side while if we are on one of
    ! its flat surfaces it is simply a vector of type (1,0,0)

    ! So first of all we have to determine in which face we are. We know that the cylinder is along the first direction
    ! defined by the box_shape_t basis vectors (which is x). To do that let us divide the vector containing the coordinates
    ! of the point by the maximum length of the parallelepiped
    normal_vector = (point_coordinates - this%center) / [this%half_length, this%radius, this%radius]

    if (abs(normal_vector(1)) > norm2(normal_vector(2:3))) then
      normal_vector(2:3) = M_ZERO
      surface_element = mesh_spacing(2) * mesh_spacing(3)
    else
      ! in this case we are on the circle
      normal_vector(1) = M_ZERO
      surface_element = mesh_spacing(1)**2
    end if

    POP_SUB(box_cylinder_shape_get_surface_point_info)
  end subroutine box_cylinder_shape_get_surface_point_info

  !--------------------------------------------------------------
  subroutine box_cylinder_write_info(this, iunit, namespace)
    class(box_cylinder_t),       intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(box_cylinder_write_info)

    write(message(1),'(2x,a)') 'Type = cylinder'
    write(message(2),'(2x,3a,f7.3)') 'Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%radius)
    write(message(3),'(2x,3a,f7.3)') 'Xlength [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%half_length)
    call messages_info(3, iunit=iunit, namespace=namespace)

    POP_SUB(box_cylinder_write_info)
  end subroutine box_cylinder_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_cylinder_short_info(this, unit_length) result(info)
    class(box_cylinder_t), intent(in) :: this
    type(unit_t),          intent(in) :: unit_length

    PUSH_SUB(box_cylinder_short_info)

    write(info, '(a,f11.6,a,a,a,f11.6,a,a)') 'BoxShape = cylinder, Radius =', units_from_atomic(unit_length, this%radius), ' ', &
      trim(units_abbrev(unit_length)), '; Xlength =', units_from_atomic(unit_length, this%half_length), ' ', &
      trim(units_abbrev(unit_length))

    POP_SUB(box_cylinder_short_info)
  end function box_cylinder_short_info

end module box_cylinder_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
