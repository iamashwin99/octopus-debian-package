!! Copyright (C) 2021 N. Tancogne-Dejean, M. Oliveira
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

module basis_vectors_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                 &
    basis_vectors_t

  !> Vectors defining a basis in a vector space. This class provides methods to
  !! convert vector coordinates in this basis to and from Cartesian coordinates.
  !! Given a vector v, its coordinates in the basis are related to the Cartesian
  !! coordinates y_i in the following way:
  !! \f[
  !!  x_i = \sum_j a_{ij} y_j
  !! \f]
  !! where a_{ij} are the elements of the change of basis matrix A.  The
  !! previous expression can be expressed in matrix form:
  !! \f[
  !! X = A Y
  !! \f]
  !! Obtaining the Cartesian coordinates from the basis coordinates is then
  !! straighforward using the inverse of A:
  !! \f[
  !! Y = A^{-1} X
  !!\f]
  !! Note that the columns of \f$A^{-1}\f$ are simply the basis vectors.
  type basis_vectors_t
    ! Components are public by default
    integer, private :: dim
    FLOAT, allocatable :: vectors(:,:)  !< the vectors of the basis
    FLOAT, allocatable :: change_of_basis_matrix(:,:)  !< the change-of-basis matrix to convert from a Cartesian basis to this basis
    logical :: orthogonal
  contains
    procedure :: copy => basis_vectors_copy
    generic   :: assignment(=) => copy
    procedure :: scale => basis_vectors_scale
    generic   :: from_cartesian => from_cartesian1, from_cartesian2
    generic   :: to_cartesian => to_cartesian1, to_cartesian2
    procedure :: from_cartesian1 => basis_vectors_from_cartesian1
    procedure :: to_cartesian1 => basis_vectors_to_cartesian1
    procedure :: from_cartesian2 => basis_vectors_from_cartesian2
    procedure :: to_cartesian2 => basis_vectors_to_cartesian2
    final :: basis_vectors_finalize
  end type basis_vectors_t

  interface basis_vectors_t
    module procedure basis_vectors_constructor
  end interface basis_vectors_t

contains

  !--------------------------------------------------------------
  type(basis_vectors_t) function basis_vectors_constructor(namespace, dim, vectors) result(basis)
    type(namespace_t), intent(in) :: namespace
    integer,           intent(in) :: dim
    FLOAT,             intent(in) :: vectors(dim, dim)

    integer :: idir1, idir2

    PUSH_SUB(basis_vectors_constructor)

    basis%dim = dim

    SAFE_ALLOCATE(basis%vectors(1:dim, 1:dim))
    SAFE_ALLOCATE(basis%change_of_basis_matrix(1:dim, 1:dim))

    basis%vectors = vectors

    basis%orthogonal = .true.
    do idir1 = 1, dim
      do idir2 = idir1 + 1, dim
        if (abs(dot_product(basis%vectors(:, idir1), basis%vectors(:, idir2))) > M_EPSILON) then
          basis%orthogonal = .false.
          exit
        end if
      end do
    end do

    if (abs(lalg_determinant(dim, basis%vectors, preserve_mat = .true.)) < M_EPSILON) then
      message(1) = "Basis vectors are not linearly independent and therefore the change-of-basis matrix cannot be calculated."
      call messages_fatal(1, namespace=namespace)
    end if
    call calculate_change_of_basis_matrix(dim, basis%vectors, basis%change_of_basis_matrix)

    POP_SUB(basis_vectors_constructor)
  end function basis_vectors_constructor

  !--------------------------------------------------------------
  subroutine basis_vectors_copy(this, source)
    class(basis_vectors_t), intent(out) :: this
    class(basis_vectors_t), intent(in)  :: source

    PUSH_SUB(basis_vectors_copy)

    this%dim = source%dim
    SAFE_ALLOCATE_SOURCE_A(this%vectors, source%vectors)
    SAFE_ALLOCATE_SOURCE_A(this%change_of_basis_matrix, source%change_of_basis_matrix)
    this%orthogonal = source%orthogonal

    POP_SUB(basis_vectors_copy)
  end subroutine basis_vectors_copy

  !--------------------------------------------------------------
  subroutine basis_vectors_finalize(this)
    type(basis_vectors_t), intent(inout) :: this

    PUSH_SUB(basis_vectors_finalize)

    SAFE_DEALLOCATE_A(this%vectors)
    SAFE_DEALLOCATE_A(this%change_of_basis_matrix)

    POP_SUB(basis_vectors_finalize)
  end subroutine basis_vectors_finalize

  !--------------------------------------------------------------
  subroutine basis_vectors_scale(this, factor)
    class(basis_vectors_t), intent(inout) :: this
    FLOAT,                  intent(in)    :: factor(1:this%dim)

    integer :: idir

    PUSH_SUB(basis_vectors_scale)

    ! Scale the basis in real space
    do idir = 1, this%dim
      this%vectors(:, idir) = this%vectors(:, idir)*factor(idir)
    end do

    ! Calculate the new change-of-basis matrix
    call calculate_change_of_basis_matrix(this%dim, this%vectors, this%change_of_basis_matrix)

    POP_SUB(basis_vectors_scale)
  end subroutine basis_vectors_scale

  !--------------------------------------------------------------
  pure function basis_vectors_from_cartesian1(this, xx_cart) result(xx)
    class(basis_vectors_t), intent(in) :: this
    FLOAT,                  intent(in) :: xx_cart(1:this%dim)
    FLOAT :: xx(1:this%dim)

    ! no PUSH_SUB, called too often

    xx = matmul(this%change_of_basis_matrix, xx_cart)

  end function basis_vectors_from_cartesian1

  !--------------------------------------------------------------
  pure function basis_vectors_to_cartesian1(this, xx) result(xx_cart)
    class(basis_vectors_t), intent(in) :: this
    FLOAT,                  intent(in) :: xx(1:this%dim)
    FLOAT :: xx_cart(1:this%dim)

    ! no PUSH_SUB, called too often

    xx_cart = matmul(this%vectors, xx)

  end function basis_vectors_to_cartesian1

  !--------------------------------------------------------------
  function basis_vectors_from_cartesian2(this, np, xx_cart) result(xx)
    class(basis_vectors_t), intent(in)    :: this
    integer,                intent(in)    :: np
    FLOAT,                  intent(in)    :: xx_cart(:,:)
    FLOAT :: xx(np, this%dim)

    ! no PUSH_SUB, called too often

    call lalg_gemm_nc(np, this%dim, this%dim, M_ONE, xx_cart, this%change_of_basis_matrix, M_ZERO, xx)

  end function basis_vectors_from_cartesian2

  !--------------------------------------------------------------
  function basis_vectors_to_cartesian2(this, np, xx) result(xx_cart)
    class(basis_vectors_t), intent(in)    :: this
    integer,                intent(in)    :: np
    FLOAT,                  intent(in)    :: xx(:,:)
    FLOAT :: xx_cart(np, this%dim)

    ! no PUSH_SUB, called too often
    call lalg_gemm_nc(np, this%dim, this%dim, M_ONE, xx, this%vectors, M_ZERO, xx_cart)

  end function basis_vectors_to_cartesian2

  !--------------------------------------------------------------
  subroutine calculate_change_of_basis_matrix(dim, vectors, matrix)
    integer,           intent(in)  :: dim
    FLOAT,             intent(in)  :: vectors(1:dim, 1:dim)
    FLOAT,             intent(out) :: matrix(1:dim, 1:dim)

    FLOAT :: volume, cross(1:3)

    PUSH_SUB(calculate_change_of_basis_matrix)

    select case (dim)
    case (3)
      cross(1:3) = dcross_product(vectors(1:3, 2), vectors(1:3, 3))
      volume = dot_product(vectors(1:3, 1), cross(1:3))

      matrix(1, 1:3) = dcross_product(vectors(:, 2), vectors(:, 3))/volume
      matrix(2, 1:3) = dcross_product(vectors(:, 3), vectors(:, 1))/volume
      matrix(3, 1:3) = dcross_product(vectors(:, 1), vectors(:, 2))/volume
    case (2)
      volume = vectors(1, 1)*vectors(2, 2) - vectors(2, 1)*vectors(1, 2)
      matrix(1, 1) =  vectors(2, 2)/volume
      matrix(1, 2) = -vectors(1, 2)/volume
      matrix(2, 1) = -vectors(2, 1)/volume
      matrix(2, 2) =  vectors(1, 1)/volume
    case (1)
      volume = vectors(1, 1)
      matrix(1, 1) = M_ONE / vectors(1, 1)
    case default ! dim > 3
      matrix = vectors
      call lalg_inverter(dim, matrix)
    end select

    POP_SUB(calculate_change_of_basis_matrix)
  end subroutine calculate_change_of_basis_matrix

end module basis_vectors_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
