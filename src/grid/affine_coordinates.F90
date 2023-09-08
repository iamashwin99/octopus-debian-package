!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module affine_coordinates_oct_m
  use basis_vectors_oct_m
  use coordinate_system_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_adv_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::   &
    affine_coordinates_t,    &
    affine_coordinates_copy

  type, extends(coordinate_system_t) :: affine_coordinates_t
    private
    type(basis_vectors_t), public :: basis
  contains
    procedure :: to_cartesian => affine_coordinates_to_cartesian
    procedure :: from_cartesian => affine_coordinates_from_cartesian
    procedure :: dvector_from_cartesian => daffine_coordinates_vector_from_cartesian
    procedure :: zvector_from_cartesian => zaffine_coordinates_vector_from_cartesian
    procedure :: dcovector_to_cartesian => daffine_coordinates_covector_to_cartesian
    procedure :: zcovector_to_cartesian => zaffine_coordinates_covector_to_cartesian
    procedure :: det_jac => affine_coordinates_det_jac
    procedure :: write_info => affine_coordinates_write_info
    procedure :: surface_element => affine_coordinates_surface_element
  end type affine_coordinates_t

  interface affine_coordinates_t
    procedure affine_coordinates_constructor
  end interface affine_coordinates_t

contains

  function affine_coordinates_constructor(namespace, dim, basis_vectors) result(affine)
    type(namespace_t),          intent(in)  :: namespace
    integer,                    intent(in)  :: dim
    FLOAT,                      intent(in)  :: basis_vectors(1:dim, 1:dim)
    class(affine_coordinates_t), pointer :: affine

    PUSH_SUB(affine_coordinates_constructor)

    SAFE_ALLOCATE(affine)

    affine%dim = dim
    affine%local_basis = .false.
    affine%basis = basis_vectors_t(namespace, dim, basis_vectors)
    affine%orthogonal = affine%basis%orthogonal

    POP_SUB(affine_coordinates_constructor)
  end function affine_coordinates_constructor

  ! --------------------------------------------------------------
  subroutine affine_coordinates_copy(this_out, this_in)
    type(affine_coordinates_t), intent(inout) :: this_out
    type(affine_coordinates_t), intent(in)    :: this_in

    PUSH_SUB(affine_coordinates_copy)

    this_out%dim = this_in%dim
    this_out%local_basis = this_in%local_basis
    this_out%basis = this_in%basis

    POP_SUB(affine_coordinates_copy)
  end subroutine affine_coordinates_copy

  ! ---------------------------------------------------------
  pure function affine_coordinates_to_cartesian(this, chi) result(xx)
    class(affine_coordinates_t), target, intent(in)  :: this
    FLOAT,                       intent(in)  :: chi(:)
    FLOAT :: xx(1:this%dim)

    ! no PUSH_SUB, called too often

    xx(:) = this%basis%to_cartesian(chi(:))

  end function affine_coordinates_to_cartesian

  ! ---------------------------------------------------------
  pure function affine_coordinates_from_cartesian(this, xx) result(chi)
    class(affine_coordinates_t), target, intent(in)  :: this
    FLOAT,                               intent(in)  :: xx(:)
    FLOAT :: chi(1:this%dim)

    ! no PUSH_SUB, called too often

    chi(:) = this%basis%from_cartesian(xx(:))

  end function affine_coordinates_from_cartesian

  ! ---------------------------------------------------------
  FLOAT function affine_coordinates_det_jac(this, xx, chi) result(jdet)
    class(affine_coordinates_t), intent(in)  :: this
    FLOAT,                       intent(in)  :: xx(:)
    FLOAT,                       intent(in)  :: chi(:)

    FLOAT :: jac(1:this%dim, 1:this%dim)

    ! No PUSH_SUB, called too often

    jac(1:this%dim, 1:this%dim) = this%basis%vectors(1:this%dim,1:this%dim)
    jdet = lalg_determinant(this%dim, jac, preserve_mat = .false.)

  end function affine_coordinates_det_jac

  ! ---------------------------------------------------------
  subroutine affine_coordinates_write_info(this, iunit, namespace)
    class(affine_coordinates_t),           intent(in) :: this
    integer,                     optional, intent(in) :: iunit
    type(namespace_t),           optional, intent(in) :: namespace

    integer :: idir1, idir2

    PUSH_SUB(affine_coordinates_write_info)

    write(message(1), '(a)')  '  Using affine coordinates'
    write(message(2), '(a)')  '  Basis vectors [',  trim(units_abbrev(units_out%length)), ']:'
    call messages_info(2, iunit=iunit, namespace=namespace)
    do idir1 = 1, this%dim
      write(message(2), '(4x,a,99(f8.3,a))') '(', &
        (units_from_atomic(units_out%length, this%basis%vectors(idir2, idir1)), idir2 = 1, this%dim - 1), &
        units_from_atomic(units_out%length, this%basis%vectors(this%dim, idir1)), ')'
      call messages_info(1, iunit=iunit, namespace=namespace)
    end do

    POP_SUB(affine_coordinates_write_info)
  end subroutine affine_coordinates_write_info

  ! ---------------------------------------------------------
  FLOAT function affine_coordinates_surface_element(this, idir) result(ds)
    class(affine_coordinates_t), intent(in)  :: this
    integer,                     intent(in)  :: idir

    PUSH_SUB(affine_coordinates_surface_element)

    select case (this%dim)
    case (3)
      select case (idir)
      case (1)
        ds = norm2(dcross_product(this%basis%vectors(1:3, 2), this%basis%vectors(1:3, 3)))
      case (2)
        ds = norm2(dcross_product(this%basis%vectors(1:3, 3), this%basis%vectors(1:3, 1)))
      case (3)
        ds = norm2(dcross_product(this%basis%vectors(1:3, 1), this%basis%vectors(1:3, 2)))
      end select

    case (2)
      select case (idir)
      case (3)
        ds = this%basis%vectors(1, 1)*this%basis%vectors(2, 2) - this%basis%vectors(2, 1)*this%basis%vectors(1, 2)
      case default
        ! We can only get the surface element along z, as the only existing plane is the xy plane
        ASSERT(.false.)
      end select

    case default
      ! We only know how to do this for 2D and 3D
      ASSERT(.false.)
    end select

    POP_SUB(affine_coordinates_surface_element)
  end function affine_coordinates_surface_element

#include "undef.F90"
#include "real.F90"
#include "affine_coordinates_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "affine_coordinates_inc.F90"

end module affine_coordinates_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
