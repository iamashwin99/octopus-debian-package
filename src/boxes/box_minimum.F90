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

module box_minimum_oct_m
  use box_oct_m
  use box_shape_oct_m
  use debug_oct_m
  use global_oct_m
  use lookup_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public :: box_minimum_t

  !> Class implementing a box that is a union of spheres. We do this in a specific class
  !! instead of using the box_union class for performance reasons (although this
  !! should be benchmarked at some point).
  type, extends(box_shape_t) :: box_minimum_t
    private
    FLOAT, public :: radius = M_ZERO !< Radius of the spheres
    integer :: n_sites                       !< How many sites there are.
    FLOAT, allocatable :: site_position(:,:) !< Site coordinates.
    type(lookup_t) :: site_lookup
  contains
    procedure :: shape_contains_points => box_minimum_shape_contains_points
    procedure :: write_info => box_minimum_write_info
    procedure :: short_info => box_minimum_short_info
    final     :: box_minimum_finalize
  end type box_minimum_t

  interface box_minimum_t
    procedure box_minimum_constructor
  end interface box_minimum_t

contains

  !--------------------------------------------------------------
  function box_minimum_constructor(dim, radius, n_sites, site_position, namespace) result(box)
    integer,                  intent(in) :: dim
    FLOAT,                    intent(in) :: radius
    integer,                  intent(in) :: n_sites
    FLOAT,                    intent(in) :: site_position(1:dim,1:n_sites)
    type(namespace_t),        intent(in) :: namespace
    class(box_minimum_t), pointer :: box

    FLOAT :: center(dim)

    PUSH_SUB(box_minimum_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    center = M_ZERO
    call box_shape_init(box, namespace, dim, center, bounding_box_min=minval(abs(site_position), dim=2) - radius, &
      bounding_box_max=maxval(abs(site_position), dim=2) + radius)
    box%radius = radius
    box%n_sites = n_sites
    SAFE_ALLOCATE_SOURCE(box%site_position, site_position)

    box%bounding_box_l = maxval(abs(site_position), dim=2) + box%radius

    call lookup_init(box%site_lookup, box%dim, box%n_sites, box%site_position)

    POP_SUB(box_minimum_constructor)
  end function box_minimum_constructor

  !--------------------------------------------------------------
  subroutine box_minimum_finalize(this)
    type(box_minimum_t), intent(inout) :: this

    PUSH_SUB(box_minimum_finalize)

    call lookup_end(this%site_lookup)

    SAFE_DEALLOCATE_A(this%site_position)

    call box_shape_end(this)

    POP_SUB(box_minimum_finalize)
  end subroutine box_minimum_finalize

  !--------------------------------------------------------------
  function box_minimum_shape_contains_points(this, nn, xx) result(contained)
    class(box_minimum_t), intent(in)  :: this
    integer,              intent(in)  :: nn
    FLOAT,                intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip
    integer, allocatable :: nlist(:)

    SAFE_ALLOCATE(nlist(1:nn))

    call lookup_get_list(this%site_lookup, nn, xx, this%radius + BOX_BOUNDARY_DELTA, nlist)

    do ip = 1, nn
      contained(ip) = nlist(ip) /= 0 .neqv. this%is_inside_out()
    end do

    SAFE_DEALLOCATE_A(nlist)

  end function box_minimum_shape_contains_points

  !--------------------------------------------------------------
  subroutine box_minimum_write_info(this, iunit, namespace)
    class(box_minimum_t),        intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(box_minimum_write_info)

    write(message(1), '(2x,a)') 'Type = minimum'
    call messages_info(1, iunit, namespace=namespace)
    write(message(1),'(2x,3a,f7.3)') 'Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%radius)
    call messages_info(1, iunit, namespace=namespace)

    POP_SUB(box_minimum_write_info)
  end subroutine box_minimum_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_minimum_short_info(this, unit_length) result(info)
    class(box_minimum_t), intent(in) :: this
    type(unit_t),         intent(in) :: unit_length

    PUSH_SUB(box_minimum_short_info)

    write(info,'(a,f11.6,a,a)') 'BoxShape = minimum; Radius =', units_from_atomic(unit_length, this%radius), ' ', &
      trim(units_abbrev(unit_length))

    POP_SUB(box_minimum_short_info)
  end function box_minimum_short_info

end module box_minimum_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
