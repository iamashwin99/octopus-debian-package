!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module box_parallelepiped_oct_m
  use basis_vectors_oct_m
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
  public :: box_parallelepiped_t

  !> Class implementing a parallelepiped box. Currently this is restricted to a
  !! rectangular cuboid (all the faces must be rectangular) and the vectors
  !! generating the parallelepiped must be along the Cartesian axes.
  type, extends(box_shape_t) :: box_parallelepiped_t
    private
    FLOAT, allocatable, public :: half_length(:) !< half the length of the parallelepiped in each direction.

    integer, public :: n_periodic_boundaries = 0 !< in how many directions the parallelepiped boundaries are periodic
  contains
    procedure :: shape_contains_points => box_parallelepiped_shape_contains_points
    procedure :: get_surface_points => box_parallelepiped_shape_get_surface_points
    procedure :: get_surface_point_info => box_parallelepiped_shape_get_surface_point_info
    procedure :: write_info => box_parallelepiped_write_info
    procedure :: short_info => box_parallelepiped_short_info
    final     :: box_parallelepiped_finalize
  end type box_parallelepiped_t

  interface box_parallelepiped_t
    procedure box_parallelepiped_constructor
  end interface box_parallelepiped_t

contains

  !--------------------------------------------------------------
  function box_parallelepiped_constructor(dim, center, axes, length, namespace, n_periodic_boundaries) result(box)
    integer,            intent(in) :: dim
    FLOAT,              intent(in) :: center(dim)
    FLOAT,              intent(in) :: axes(dim, dim) !< the A matrix from Chelikowski PRB 78 075109 (2008)
    FLOAT,              intent(in) :: length(dim)             !< length of the parallelepiped along each basis vector
    type(namespace_t),  intent(in) :: namespace
    integer, optional,  intent(in) :: n_periodic_boundaries   !< in how many directions the parallelepiped boundaries are periodic
    class(box_parallelepiped_t), pointer :: box

    PUSH_SUB(box_parallelepiped_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    SAFE_ALLOCATE(box%half_length(1:dim))
    box%half_length = M_HALF * length
    if (present(n_periodic_boundaries)) then
      box%n_periodic_boundaries = n_periodic_boundaries
    end if
    call box_shape_init(box, namespace, dim, center, bounding_box_min=-box%half_length, bounding_box_max=box%half_length, axes=axes)

    box%bounding_box_l = box%half_length + abs(center)

    POP_SUB(box_parallelepiped_constructor)
  end function box_parallelepiped_constructor

  !--------------------------------------------------------------
  subroutine box_parallelepiped_finalize(this)
    type(box_parallelepiped_t), intent(inout) :: this

    PUSH_SUB(box_parallelepiped_finalize)

    call box_shape_end(this)
    SAFE_DEALLOCATE_A(this%half_length)

    POP_SUB(box_parallelepiped_finalize)
  end subroutine box_parallelepiped_finalize

  !--------------------------------------------------------------
  function box_parallelepiped_shape_contains_points(this, nn, xx) result(contained)
    class(box_parallelepiped_t), intent(in)  :: this
    integer,                     intent(in)  :: nn
    FLOAT,                       intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip, idir
    FLOAT :: ulimit(this%dim), llimit(this%dim)

    llimit = this%center - this%half_length - BOX_BOUNDARY_DELTA
    do idir = 1, this%dim
      if (idir <= this%n_periodic_boundaries) then
        ! When periodic, we exclude one of the faces from the box.
        ulimit(idir) = this%center(idir) + this%half_length(idir) - BOX_BOUNDARY_DELTA
      else
        ulimit(idir) = this%center(idir) + this%half_length(idir) + BOX_BOUNDARY_DELTA
      end if
    end do

    do ip = 1, nn
      contained(ip) = .true.
      do idir = 1, this%dim
        contained(ip) = contained(ip) .and. xx(ip, idir) >= llimit(idir) .and. xx(ip, idir) <= ulimit(idir)
        if (idir > this%n_periodic_boundaries) then
          ! We only consider the box to be inside out along the non-periodic directions.
          contained(ip) = contained(ip) .neqv. this%is_inside_out()
        end if
      end do
    end do

  end function box_parallelepiped_shape_contains_points

  !--------------------------------------------------------------
  !> Get a mask for the grid points telling which of them are surface points
  !! 1. Create a box that is fractionally smaller than the original.
  !! 2. Look for points that belong to the original box but not to the smaller one - These define surface points
  function box_parallelepiped_shape_get_surface_points(this,namespace,mesh_spacing,nn,xx, number_of_layers) result(surface_points)
    class(box_parallelepiped_t), intent(in)  :: this
    type(namespace_t),  intent(in)  :: namespace
    FLOAT,              intent(in)  :: mesh_spacing(:)
    integer,            intent(in)  :: nn
    FLOAT,              intent(in)  :: xx(:,:)
    integer, optional,  intent(in)  :: number_of_layers

    logical :: surface_points(1:nn)
    integer :: idir, number_of_layers_
    FLOAT   :: shrink(this%dim), test_axis(this%dim)
    class(box_parallelepiped_t), pointer :: shrinked_box

    ! Check that each axis is along the Cartersian axis
    do idir = 1, this%dim
      test_axis = M_ZERO
      test_axis(idir) = M_ONE
      ASSERT(all(abs(this%axes%vectors(:, idir) - test_axis(:)) < M_EPSILON))
    end do

    number_of_layers_ = 1
    if (present(number_of_layers)) number_of_layers_ = number_of_layers

    ! Determine how much the shrink should be based on the spacing in each direction
    shrink(:) = 1 - number_of_layers_ * (1 - BOX_BOUNDARY_DELTA) * (mesh_spacing(:) / this%half_length(:))

    ! Then we create a shrunk parallelepiped box
    shrinked_box => box_parallelepiped_t(this%dim, this%center, this%axes%vectors, &
      M_TWO*this%half_length(:)*shrink(:), namespace)

    ! Then we look for points of the old box that are not containted in the smaller box
    ! box_parallelepiped_shape_contains_points will return the point contained in the shrinked box,
    ! but we are interested in the ones not contained, as they will the surface points of the original parallelepiped
    if (SIZE(xx, 1) == 3) then
      surface_points = .not. box_parallelepiped_shape_contains_points(shrinked_box, nn, transpose(xx))
    else
      surface_points = .not. box_parallelepiped_shape_contains_points(shrinked_box, nn, xx)
    end if

    SAFE_DEALLOCATE_P(shrinked_box)
  end function box_parallelepiped_shape_get_surface_points

  !--------------------------------------------------------------
  subroutine box_parallelepiped_shape_get_surface_point_info(this, point_coordinates,mesh_spacing, normal_vector, surface_element)
    class(box_parallelepiped_t), intent(in)  :: this
    FLOAT,                       intent(in)  :: point_coordinates(:) !< (x,y,z) coordinates of the point
    FLOAT,                       intent(in)  :: mesh_spacing(:)      !< spacing of the mesh
    FLOAT,                       intent(out) :: normal_vector(:)     !< normal vector to the surface point
    FLOAT,                       intent(out) :: surface_element      !< surface element (needed to compute the surface integral)

    FLOAT   :: vector_norm
    integer :: idir, first_index, second_index, third_index, face_counter

    PUSH_SUB(box_parallelepiped_shape_get_surface_point_info)

    ! Compute the normal vector to the surface point

    ! For a parallelepiped, the normal vector to a point is simply a vector that has plus or minus one in the direction at
    ! which the face is pointing, and zero in the others. For instance, if the top face is facing the z direction, then
    ! the normal vector for all points in that face should be (0,0,1)

    ! So first of all we have to determine in which face we are. To do that let us divide the vector containing the coordinates
    ! of the point by the maximum length of the parallelepiped
    normal_vector(:) = (point_coordinates(:) - this%center(:)) / this%half_length(:)
    face_counter = M_ZERO
    surface_element = M_ZERO

    ! Loop over the three directions
    do idir = 1, this%dim
      ! Determine the indeces of the faces
      first_index = this%dim - mod(idir, this%dim)
      second_index = this%dim - mod(idir + 1, this%dim)
      third_index = this%dim - mod(idir + 2, this%dim)
      if (abs(normal_vector(first_index)) >= abs(normal_vector(second_index)) .and. &
        abs(normal_vector(first_index)) >= abs(normal_vector(third_index))) then
        ! Define the surface element for the first_index
        surface_element = mesh_spacing(second_index) * mesh_spacing(third_index)
        face_counter = face_counter + 1
        ! Now we check whether we are on the face or on a vertex or corner
        if (abs(normal_vector(second_index)) / abs(normal_vector(first_index)) < &
          M_ONE - M_HALF * mesh_spacing(second_index) / this%half_length(first_index)) then
          normal_vector(second_index) = M_ZERO
        end if
        if (abs(normal_vector(third_index)) / abs(normal_vector(first_index)) < &
          M_ONE - M_HALF * mesh_spacing(third_index) / this%half_length(first_index)) then
          normal_vector(third_index) = M_ZERO
        end if
      end if
    end do

    if (face_counter == 3) surface_element = surface_element * (M_THREE / M_FOUR)

    ! Now if the point lies in a face, then the normalized_coordinates has one component which is one, and all the others are
    ! zero. If the point is on a vertex, more than one coordinate is non zero. the final step is to normalize the vector
    vector_norm = norm2(normal_vector(:))
    ! Normal vector should never be zero
    ASSERT(vector_norm > M_EPSILON)
    normal_vector(:) = normal_vector(:) / vector_norm

    POP_SUB(box_parallelepiped_shape_get_surface_point_info)
  end subroutine box_parallelepiped_shape_get_surface_point_info

  !--------------------------------------------------------------
  subroutine box_parallelepiped_write_info(this, iunit, namespace)
    class(box_parallelepiped_t), intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer :: idir

    PUSH_SUB(box_parallelepiped_write_info)

    write(message(1),'(2x,a)') 'Type = parallelepiped'
    write(message(2),'(2x,3a, 99(f8.3,a))') 'Lengths [', trim(units_abbrev(units_out%length)), '] = (', &
      (units_from_atomic(units_out%length, M_TWO*this%half_length(idir)), ',', idir = 1, this%dim - 1), &
      units_from_atomic(units_out%length, M_TWO*this%half_length(this%dim)), ')'
    call messages_info(2, iunit=iunit, namespace=namespace)

    POP_SUB(box_parallelepiped_write_info)
  end subroutine box_parallelepiped_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_parallelepiped_short_info(this, unit_length) result(info)
    class(box_parallelepiped_t), intent(in) :: this
    type(unit_t),                intent(in) :: unit_length

    integer :: idir

    PUSH_SUB(box_parallelepiped_short_info)

    write(info, '(a,a,a,99(f11.6,a))') 'BoxShape = parallelepiped; Lengths [', trim(units_abbrev(unit_length)),'] = [', &
      (units_from_atomic(unit_length, M_TWO*this%half_length(idir)), ',', idir = 1, this%dim - 1), &
      units_from_atomic(unit_length, M_TWO*this%half_length(this%dim)), ']'

    POP_SUB(box_parallelepiped_short_info)
  end function box_parallelepiped_short_info

end module box_parallelepiped_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
