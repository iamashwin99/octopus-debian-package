!! Copyright (C) 2021 M. Oliveira
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

module coordinate_system_oct_m
  use namespace_oct_m
  implicit none

  private
  public :: coordinate_system_t

  type, abstract :: coordinate_system_t
    logical :: local_basis  !< Do the basis vectors depend on the position, i.e., is the basis local? (false for Cartesian and affine, true for curvilinear coordinates in general)
    logical :: orthogonal   !< Are the basis vectors orthogonal?
    integer :: dim          !< Dimension of the space
    FLOAT :: min_mesh_scaling_product !< product of the smallest scaling :: min(distance between the grid points / spacing)
  contains
    generic   :: vector_from_cartesian => dvector_from_cartesian, zvector_from_cartesian
    procedure :: dvector_from_cartesian => dcoordinate_system_vector_from_cartesian
    procedure :: zvector_from_cartesian => zcoordinate_system_vector_from_cartesian
    generic   :: covector_to_cartesian => dcovector_to_cartesian, zcovector_to_cartesian
    procedure :: dcovector_to_cartesian => dcoordinate_system_covector_to_cartesian
    procedure :: zcovector_to_cartesian => zcoordinate_system_covector_to_cartesian
    procedure(coordinate_system_to_cartesian),           deferred :: to_cartesian
    procedure(coordinate_system_from_cartesian),         deferred :: from_cartesian
    procedure(coordinate_system_det_jac),                deferred :: det_jac
    procedure(coordinate_system_write_info),             deferred :: write_info
    procedure(coordinates_surface_element),              deferred :: surface_element
  end type coordinate_system_t

  abstract interface
    ! ---------------------------------------------------------
    ! Convert coordinates given in this coordinate system to Cartesian
    ! coordinates
    function coordinate_system_to_cartesian(this, chi) result(xx)
      import coordinate_system_t
      class(coordinate_system_t), target, intent(in)  :: this
      FLOAT,                              intent(in)  :: chi(:)
      FLOAT :: xx(1:this%dim)
    end function coordinate_system_to_cartesian

    ! ---------------------------------------------------------
    ! Convert Cartesian coordinates to coordinates in this coordinate system
    function coordinate_system_from_cartesian(this, xx) result(chi)
      import coordinate_system_t
      class(coordinate_system_t), target, intent(in)  :: this
      FLOAT,                              intent(in)  :: xx(:)
      FLOAT :: chi(1:this%dim)
    end function coordinate_system_from_cartesian

    ! ---------------------------------------------------------
    FLOAT function coordinate_system_det_jac(this, xx, chi) result(jdet)
      import coordinate_system_t
      class(coordinate_system_t), intent(in)  :: this
      FLOAT,                      intent(in)  :: xx(:)
      FLOAT,                      intent(in)  :: chi(:)
    end function coordinate_system_det_jac

    ! ---------------------------------------------------------
    subroutine coordinate_system_write_info(this, iunit, namespace)
      import coordinate_system_t
      import namespace_t
      class(coordinate_system_t),           intent(in) :: this
      integer,                    optional, intent(in) :: iunit
      type(namespace_t),          optional, intent(in) :: namespace
    end subroutine coordinate_system_write_info

    ! ---------------------------------------------------------
    FLOAT function coordinates_surface_element(this, idir) result(ds)
      import coordinate_system_t
      class(coordinate_system_t), intent(in) :: this
      integer,                    intent(in) :: idir
    end function coordinates_surface_element
  end interface

contains

  ! ---------------------------------------------------------
  subroutine dcoordinate_system_vector_from_cartesian(this, xx, vv, src)
    class(coordinate_system_t), intent(in)    :: this
    FLOAT,                      intent(in)    :: xx(:)
    FLOAT,                      intent(inout) :: vv(:)
    FLOAT, optional,            intent(in)    :: src(:)

    ! Classes that want to provide this method must override it.
    ASSERT(.false.)

  end subroutine dcoordinate_system_vector_from_cartesian

  ! ---------------------------------------------------------
  subroutine zcoordinate_system_vector_from_cartesian(this, xx, vv, src)
    class(coordinate_system_t), intent(in)    :: this
    FLOAT,                      intent(in)    :: xx(:)
    CMPLX,                      intent(inout) :: vv(:)
    CMPLX, optional,            intent(in)    :: src(:)

    ! Classes that want to provide this method must override it.
    ! Note that when the basis is local, this transformantion depends on the
    ! position
    ASSERT(.false.)

  end subroutine zcoordinate_system_vector_from_cartesian

  ! ---------------------------------------------------------
  subroutine dcoordinate_system_covector_to_cartesian(this, xx, cv, src)
    class(coordinate_system_t), intent(in)    :: this
    FLOAT,                      intent(in)    :: xx(:)
    FLOAT,                      intent(inout) :: cv(:)
    FLOAT, optional,            intent(in)    :: src(:)

    ! Classes that want to provide this method must override it.
    ASSERT(.false.)

  end subroutine dcoordinate_system_covector_to_cartesian

  ! ---------------------------------------------------------
  subroutine zcoordinate_system_covector_to_cartesian(this, xx, cv, src)
    class(coordinate_system_t), intent(in)    :: this
    FLOAT,                      intent(in)    :: xx(:)
    CMPLX,                      intent(inout) :: cv(:)
    CMPLX, optional,            intent(in)    :: src(:)

    ! Classes that want to provide this method must override it.
    ASSERT(.false.)

  end subroutine zcoordinate_system_covector_to_cartesian

end module coordinate_system_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
