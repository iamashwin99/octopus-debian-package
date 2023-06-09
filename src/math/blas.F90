!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

! -----------------------------------------------------------------------
!> This module contains interfaces for BLAS routines
!! You should not use these routines directly. Please use the lalg_XXXX
! -----------------------------------------------------------------------
module blas_oct_m

  implicit none

  public ! only interfaces in this module

  ! ---------------------------------------------------------------------
  ! BLAS level I
  ! ---------------------------------------------------------------------

  !> ----------------- swap ------------------
  !! Interchanges two vectors.
  interface blas_swap
    subroutine dswap(n, dx, incx, dy, incy)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n, incx, incy
      real(r8),    intent(inout) :: dx, dy !< dx(n), dy(n)
    end subroutine dswap

    subroutine zswap(n, dx, incx, dy, incy)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n, incx, incy
      complex(r8), intent(inout) :: dx, dy !< dx(n), dy(n)
    end subroutine zswap
  end interface blas_swap

  !> ----------------- scal ------------------
  !! Scales a vector by a constant.
  interface blas_scal
    subroutine dscal(n, da, dx, incx)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n, incx
      real(r8),    intent(in)    :: da
      real(r8),    intent(inout) :: dx !< dx(n)
    end subroutine dscal

    subroutine zscal(n, da, dx, incx)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n, incx
      complex(r8), intent(in)    :: da
      complex(r8), intent(inout) :: dx !< dx(n)
    end subroutine zscal

    subroutine dazscal(n, da, dx)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n
      real(r8),    intent(in)    :: da
      complex(r8), intent(inout) :: dx !< dx(n)
    end subroutine dazscal
  end interface blas_scal

  !> ----------------- axpy ------------------
  !! Constant times a vector plus a vector.
  interface blas_axpy
    subroutine daxpy (n, da, dx, incx, dy, incy)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n, incx, incy
      real(r8),    intent(in)    :: da, dx !< dx(n)
      real(r8),    intent(inout) :: dy     !< dy(n)
    end subroutine daxpy

    subroutine zaxpy (n, da, dx, incx, dy, incy)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n, incx, incy
      complex(r8), intent(in)    :: da, dx !< dx(n)
      complex(r8), intent(inout) :: dy     !< dy(n)
    end subroutine zaxpy

    subroutine dazaxpy (n, da, dx, dy)
      use kind_oct_m
      implicit none
      integer,    intent(in)    :: n
      real(r8),    intent(in)    :: da
      complex(r8), intent(in)    :: dx     !< dx(n)
      complex(r8), intent(inout) :: dy     !< dy(n)
    end subroutine dazaxpy
  end interface blas_axpy

  !> ----------------- copy ------------------
  !! Copies a vector, x, to a vector, y.
  interface blas_copy
    subroutine dcopy(n, dx, incx, dy, incy)
      use kind_oct_m
      implicit none
      integer,    intent(in)  :: n, incx, incy
      real(r8),    intent(in)  :: dx !< dx(n)
      real(r8),    intent(out) :: dy !< dy(n)
    end subroutine dcopy

    subroutine zcopy(n, dx, incx, dy, incy)
      use kind_oct_m
      implicit none
      integer,    intent(in)  :: n, incx, incy
      complex(r8), intent(in)  :: dx !< dx(n)
      complex(r8), intent(out) :: dy !< dy(n)
    end subroutine zcopy
  end interface blas_copy

  !> ----------------- dot  ------------------
  !! Forms the dot product of two vectors.
  interface blas_dot
    real(r8) function ddot(n, dx, incx, dy, incy)
      use kind_oct_m
      integer,    intent(in) :: n, incx, incy
      real(r8),    intent(in) :: dx, dy !< dx(n), dy(n)
    end function ddot

    complex(r8) function zdotc(n, dx, incx, dy, incy)
      use kind_oct_m
      integer,    intent(in) :: n, incx, incy
      complex(r8), intent(in) :: dx, dy !< dx(n), dy(n)
    end function zdotc
  end interface blas_dot

  interface blas_dotu
    complex(r8) function zdotu(n, dx, incx, dy, incy)
      use kind_oct_m
      integer,    intent(in) :: n, incx, incy
      complex(r8), intent(in) :: dx, dy !< dx(n), dy(n)
    end function zdotu
  end interface blas_dotu

  !> ----------------- nrm2 ------------------
  !! Returns the euclidean norm of a vector via the function
  !! name, so that
  !! \f[
  !! SNRM2 := sqrt( x`*x )
  !! \f]
  interface blas_nrm2
    real(r8) function dnrm2(n, dx, incx)
      use kind_oct_m
      integer,    intent(in) :: n, incx
      real(r8),    intent(in) :: dx !< dx(n)
    end function dnrm2

    real(r8) function dznrm2(n, dx, incx)
      use kind_oct_m
      integer,    intent(in) :: n, incx
      complex(r8), intent(in) :: dx !< dx(n)
    end function dznrm2
  end interface blas_nrm2


  ! ------------------------------------------------------------------
  ! BLAS level II
  ! ------------------------------------------------------------------

  !> ----------------- symv ------------------
  !! performs the matrix-vector  operation
  !!
  !! \f[
  !!     y := \alpha A x + \beta y
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, x and y are n
  !!  element vectors and A is an \f$n\times n\f$ symmetric matrix.
  interface blas_symv
    subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, incx, incy
      real(r8),      intent(in)    :: alpha, beta
      real(r8),      intent(in)    :: a !< a(lda,n)
      real(r8),      intent(in)    :: x !< x(:)
      real(r8),      intent(inout) :: y !< y(:)
    end subroutine dsymv

    subroutine zsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, incx, incy
      complex(r8),   intent(in)    :: alpha, beta
      complex(r8),   intent(in)    :: a !< a(lda,n)
      complex(r8),   intent(in)    :: x !< x(:)
      complex(r8),   intent(inout) :: y !< y(:)
    end subroutine zsymv
  end interface blas_symv

  !> ----------------- gemv ------------------
  !! SGEMV  performs one of the matrix-vector operations
  !!
  !! \f[
  !!     y := \alpha A x + \beta y,
  !! \f]
  !! or
  !! \f[
  !!     y := \alpha A^Tx + \beta y
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, x and y are
  !!  vectors and A is an \f$m\times n\f$ matrix.
  interface blas_gemv
    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      real(r8),      intent(in)    :: alpha, beta
      real(r8),      intent(in)    :: a !< a(lda,n)
      real(r8),      intent(in)    :: x !< x(:)
      real(r8),      intent(inout) :: y !< y(:)
    end subroutine dgemv

    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      complex(r8),   intent(in)    :: alpha, beta
      complex(r8),   intent(in)    :: a !< a(lda,n)
      complex(r8),   intent(in)    :: x !< x(:)
      complex(r8),   intent(inout) :: y !< y(:)
    end subroutine zgemv
  end interface blas_gemv


  ! ------------------------------------------------------------------
  ! BLAS level III
  ! ------------------------------------------------------------------

  !> ----------------- gemm ------------------
  !! performs one of the matrix-matrix operations
  !!
  !! \f[
  !!     C := \alpha op( A ) op( B ) + \beta C,
  !! \f]
  !!
  !!  where  op(X) is one of
  !!
  !! \f[
  !!     op( X ) = X   \mbox{ or }   op( X ) = X^T,
  !! \f]
  !!
  !!  \f$ \alpha \f$ and \f$ \beta \f$ are scalars, and A, B and C are
  !!  matrices, with op(A) an \f$ m \times k \f$ matrix, op(B) a \f$ k
  !!  \times n \f$ matrix and C an \f$ m \times n \f$m matrix.
  interface blas_gemm
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      real(r8),      intent(in)    :: alpha, beta
      real(r8),      intent(in)    :: a !< a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      real(r8),      intent(in)    :: b !< b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      real(r8),      intent(inout) :: c !< c(ldc,n)
    end subroutine dgemm

    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      complex(r8),   intent(in)    :: alpha, beta
      complex(r8),   intent(in)    :: a !< a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      complex(r8),   intent(in)    :: b !< b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      complex(r8),   intent(inout) :: c !< c(ldc,n)
    end subroutine zgemm

    subroutine zdgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      real(r8),      intent(in)    :: alpha, beta
      complex(r8),   intent(in)    :: a !< a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      real(r8),      intent(in)    :: b !< b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      complex(r8),   intent(inout) :: c !< c(ldc,n)
    end subroutine zdgemm
  end interface

  !> ----------------- trmm ------------------
  !! Performs one of the matrix-matrix operations
  !!
  !! \f[
  !!     B := \alpha op( A )B,  \mbox{ or }  B := \alpha B op( A ),
  !! \f]
  !!
  !!  where \f$\alpha\f$ is a scalar, B is an \f$m\times n\f$matrix, A
  !!  is a unit, or non-unit, upper or lower triangular matrix and op(
  !!  A ) is one of
  !!
  !! \f[
  !!     op( A ) = A   \mbox{ or }   op( A ) = A^T.
  !! \f]
  interface blas_trmm
    subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: side, uplo, transa, diag
      integer,      intent(in)    :: m, n, lda, ldb
      real(r8),      intent(in)    :: a, alpha
      real(r8),      intent(inout) :: b
    end subroutine dtrmm

    subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: side, uplo, transa, diag
      integer,      intent(in)    :: m, n, lda, ldb
      complex(r8),   intent(in)    :: a, alpha
      complex(r8),   intent(inout) :: b
    end subroutine ztrmm
  end interface blas_trmm

  !> ----------------- symm, hemm ------------------
  !! performs one of the matrix-matrix operations
  !!
  !! \f[
  !!     C := \alpha A B + \beta C,
  !! \f]
  !!
  !!  or
  !!
  !! \f[
  !!     C := \alpha  B A + \beta C
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, A is a symmetric
  !!  matrix and B and C are \f$m\times n\f$ matrices.
  interface blas_symm
    subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: side, uplo
      integer,      intent(in)    :: m, n, lda, ldb, ldc
      real(r8),      intent(in)    :: alpha, beta, a, b
      real(r8),      intent(inout) :: c
    end subroutine dsymm

    subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: side, uplo
      integer,      intent(in)    :: m, n, lda, ldb, ldc
      complex(r8),   intent(in)    :: alpha, beta, a, b
      complex(r8),   intent(inout) :: c
    end subroutine zsymm
  end interface blas_symm

  !> ----------------- syrk, herk ------------------
  !! performs one of the symmetric rank k operations
  !!
  !! \f[
  !!     C := \alpha A A^T + \beta C,
  !!
  !! \f]
  !! or
  !! \f[
  !!     C := \alpha A^T A + \beta*C
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, C is an
  !!  \f$n\times n\f$ symmetric matrix and A is an \f$n\times k\f$
  !!  matrix in the first case and a \f$k\times n\f$ matrix in the
  !!  second case.
  interface blas_herk
    subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo, trans
      integer,      intent(in)    :: n, k, lda, ldc
      real(r8),      intent(in)    :: alpha, beta, a
      real(r8),      intent(inout) :: c
    end subroutine dsyrk

    subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo, trans
      integer,      intent(in)    :: n, k, lda, ldc
      real(r8),      intent(in)    :: alpha, beta
      complex(r8),   intent(in)    :: a
      complex(r8),   intent(inout) :: c
    end subroutine zherk
  end interface blas_herk

  !> -----------------------trsm-------------------------
  !! Solves one of the matrix equations
  !!
  !! \f[
  !!     op( A )X = \alpha B,
  !! \f]
  !! or
  !! \f[
  !! X op( A ) = \alpha B,
  !! \f]
  !!
  !!  where \f$\alpha\f$ is a scalar, X and B are \f$m\times n\f$
  !!  matrices, A is a unit, or non-unit, upper or lower triangular
  !!  matrix and op(A) is one of
  !!
  !! \f[
  !!     op( A ) = A   \mbox{ or }   op( A ) = A^T
  !! \f]
  !!
  !!  The matrix X is overwritten on B.
  interface blas_trsm
    subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      real(r8),      intent(in)    :: alpha
      real(r8),      intent(in)    :: a
      integer,      intent(in)    :: lda
      real(r8),      intent(inout) :: b
      integer,      intent(in)    :: ldb
    end subroutine dtrsm

    subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      complex(r8),   intent(in)    :: alpha
      complex(r8),   intent(in)    :: a
      integer,      intent(in)    :: lda
      complex(r8),   intent(inout) :: b
      integer,      intent(in)    :: ldb
    end subroutine ztrsm
  end interface blas_trsm

end module blas_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
