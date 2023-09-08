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

module cartesian_oct_m
  use affine_coordinates_oct_m
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
    cartesian_t,    &
    cartesian_copy

  type, extends(affine_coordinates_t) :: cartesian_t
    private
  contains
    procedure :: to_cartesian => cartesian_to_cartesian
    procedure :: from_cartesian => cartesian_from_cartesian
    procedure :: dvector_from_cartesian => dcartesian_vector_from_cartesian
    procedure :: zvector_from_cartesian => zcartesian_vector_from_cartesian
    procedure :: dcovector_to_cartesian => dcartesian_covector_to_cartesian
    procedure :: zcovector_to_cartesian => zcartesian_covector_to_cartesian
    procedure :: det_jac => cartesian_det_jac
    procedure :: write_info => cartesian_write_info
  end type cartesian_t

  interface cartesian_t
    procedure cartesian_constructor
  end interface cartesian_t

contains

  ! ---------------------------------------------------------
  function cartesian_constructor(namespace, dim) result(cart)
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(in)  :: dim
    class(cartesian_t), pointer :: cart

    PUSH_SUB(cartesian_constructor)

    SAFE_ALLOCATE(cart)

    cart%dim = dim
    cart%local_basis = .false.
    cart%orthogonal = .true.
    cart%basis = basis_vectors_t(namespace, dim, diagonal_matrix(dim, M_ONE))

    POP_SUB(cartesian_constructor)
  end function cartesian_constructor

  ! --------------------------------------------------------------
  subroutine cartesian_copy(this_out, this_in)
    type(cartesian_t), intent(inout) :: this_out
    type(cartesian_t), intent(in)    :: this_in

    PUSH_SUB(cartesian_copy)

    this_out%dim = this_in%dim
    this_out%local_basis = this_in%local_basis

    POP_SUB(cartesian_copy)
  end subroutine cartesian_copy

  ! ---------------------------------------------------------
  pure function cartesian_to_cartesian(this, chi) result(xx)
    class(cartesian_t), target, intent(in)  :: this
    FLOAT,                      intent(in)  :: chi(:)
    FLOAT :: xx(1:this%dim)

    ! no PUSH_SUB, called too often

    xx = chi

  end function cartesian_to_cartesian

  ! ---------------------------------------------------------
  pure function cartesian_from_cartesian(this, xx) result(chi)
    class(cartesian_t), target, intent(in)  :: this
    FLOAT,                      intent(in)  :: xx(:)
    FLOAT :: chi(1:this%dim)

    ! no PUSH_SUB, called too often

    chi = xx

  end function cartesian_from_cartesian

  ! ---------------------------------------------------------
  pure subroutine dcartesian_vector_from_cartesian(this, xx, vv, src)
    class(cartesian_t), intent(in)    :: this
    FLOAT,              intent(in)    :: xx(:)
    FLOAT,              intent(inout) :: vv(:)
    FLOAT, optional,    intent(in)    :: src(:)

    ! no PUSH_SUB, called too often

    ! We are already in Cartesian coordinates, so nothing to do
    if (present(src)) then
      vv = src
    end if

  end subroutine dcartesian_vector_from_cartesian

  ! ---------------------------------------------------------
  pure subroutine zcartesian_vector_from_cartesian(this, xx, vv, src)
    class(cartesian_t), intent(in)    :: this
    FLOAT,              intent(in)    :: xx(:)
    CMPLX,              intent(inout) :: vv(:)
    CMPLX, optional,    intent(in)    :: src(:)

    ! no PUSH_SUB, called too often

    ! We are already in Cartesian coordinates, so nothing to do
    if (present(src)) then
      vv = src
    end if

  end subroutine zcartesian_vector_from_cartesian

  ! ---------------------------------------------------------
  pure subroutine dcartesian_covector_to_cartesian(this, xx, cv, src)
    class(cartesian_t), intent(in)    :: this
    FLOAT,              intent(in)    :: xx(:)
    FLOAT,              intent(inout) :: cv(:)
    FLOAT, optional,    intent(in)    :: src(:)

    ! no PUSH_SUB, called too often

    ! We are already in Cartesian coordinates, so nothing to do
    if (present(src)) then
      cv = src
    end if

  end subroutine dcartesian_covector_to_cartesian

  ! ---------------------------------------------------------
  pure subroutine zcartesian_covector_to_cartesian(this, xx, cv, src)
    class(cartesian_t), intent(in)    :: this
    FLOAT,              intent(in)    :: xx(:)
    CMPLX,              intent(inout) :: cv(:)
    CMPLX, optional,    intent(in)    :: src(:)

    ! no PUSH_SUB, called too often

    ! We are already in Cartesian coordinates, so nothing to do
    if (present(src)) then
      cv = src
    end if

  end subroutine zcartesian_covector_to_cartesian

  ! ---------------------------------------------------------
  pure FLOAT function cartesian_det_jac(this, xx, chi) result(jdet)
    class(cartesian_t),    intent(in)  :: this
    FLOAT,                 intent(in)  :: xx(:)
    FLOAT,                 intent(in)  :: chi(:)

    ! No PUSH_SUB, called too often

    jdet = M_ONE

  end function cartesian_det_jac

  ! ---------------------------------------------------------
  subroutine cartesian_write_info(this, iunit, namespace)
    class(cartesian_t),           intent(in) :: this
    integer,            optional, intent(in) :: iunit
    type(namespace_t),  optional, intent(in) :: namespace

    PUSH_SUB(cartesian_write_info)

    write(message(1), '(a)')  '  Using Cartesian coordinates'
    call messages_info(1, iunit=iunit, namespace=namespace)

    POP_SUB(cartesian_write_info)
  end subroutine cartesian_write_info

end module cartesian_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
