!! Copyright (C) 2011 X. Andrade
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

module pblas_oct_m
  implicit none

  private

  public ::            &
    pblas_gemm,        &
    pblas_herk,        &
    pblas_trsm,        &
    pblas_tran


  !> Performs any one of the following combined matrix computations:
  !!
  !! \f[
  !!   C := \alpha A B + \beta C
  !! \f]
  !! \f[
  !!  C := \alpha A B^T + \beta C
  !! \f]
  !! \f[
  !!  C := \alpha A^T B + \beta C
  !! \f]
  !! \f[
  !!  C := \alpha A^T B^T + \beta C
  !! \f]
  interface pblas_gemm
    subroutine pdgemm(transa, transb, m, n, k, alpha, &
      a, ia, ja, desca, b, ib, jb, descb, &
      beta, c, ic, jc, descc)
      use kind_oct_m
      implicit none

      character(len=1) :: transa, transb
      integer          :: m, n, k, ia, ja, ib, jb, ic, jc
      real(r8)          :: alpha, beta
      real(r8)          :: a, b, c
      integer          :: desca, descb, descc
    end subroutine pdgemm

    subroutine pzgemm(transa, transb, m, n, k, alpha, &
      a, ia, ja, desca, b, ib, jb, descb, &
      beta, c, ic, jc, descc)
      use kind_oct_m
      implicit none

      character(len=1) :: transa, transb
      integer          :: m, n, k, ia, ja, ib, jb, ic, jc
      complex(r8)       :: alpha, beta
      complex(r8)       :: a, b, c
      integer          :: desca, descb, descc
    end subroutine pzgemm
  end interface

  !>  Performs one of the symmetric rank k operations
  !!
  !! \f[
  !!     C := \alpha A A^T + \beta C,
  !! \f]
  !!  or
  !! \f[
  !!     C := \alpha A^T A + \beta C
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, C is an
  !!  \f$n\times n\f$ symmetric matrix and A is an \f$n\times k\f$
  !!  matrix in the first case and A \f$k \times n\f$ matrix in the
  !!  second case.
  interface pblas_herk
    subroutine pdsyrk(uplo, trans, n, k, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: n
      integer,      intent(in)    :: k
      real(r8),      intent(in)    :: alpha
      real(r8),      intent(in)    :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      real(r8),      intent(in)    :: beta
      real(r8),      intent(inout) :: c
      integer,      intent(in)    :: ic
      integer,      intent(in)    :: jc
      integer,      intent(in)    :: descc
    end subroutine pdsyrk

    subroutine pzherk(uplo, trans, n, k, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: n
      integer,      intent(in)    :: k
      complex(r8),   intent(in)    :: alpha
      complex(r8),   intent(in)    :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      complex(r8),   intent(in)    :: beta
      complex(r8),   intent(inout) :: c
      integer,      intent(in)    :: ic
      integer,      intent(in)    :: jc
      integer,      intent(in)    :: descc
    end subroutine pzherk

  end interface pblas_herk

  !> Solves one of the matrix equations
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
  interface pblas_trsm
    subroutine pdtrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
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
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      real(r8),      intent(inout) :: b
      integer,      intent(in)    :: ib
      integer,      intent(in)    :: jb
      integer,      intent(in)    :: descb
    end subroutine pdtrsm

    subroutine pztrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
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
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      complex(r8),   intent(inout) :: b
      integer,      intent(in)    :: ib
      integer,      intent(in)    :: jb
      integer,      intent(in)    :: descb
    end subroutine pztrsm

  end interface pblas_trsm

  !  PDTRAN  transposes a matrix
  !
  !     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
  !
  !  where
  !
  !     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
  !
  !     sub( A ) denotes A(IA:IA+N-1,JA:JA+M-1), and, op( X ) = X**T.
  !
  !  Thus, op( sub( A ) ) denotes A(IA:IA+N-1,JA:JA+M-1)**T.
  !
  !  Beta is a scalar, sub( C ) is an m by n submatrix, and sub( A ) is an
  !  n by m submatrix.
  interface pblas_tran
    subroutine pdtran(M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC)
      use kind_oct_m
      implicit none
      integer :: M, N, IA, JA, DESCA, IC, JC, DESCC
      real(r8) :: ALPHA, A, BETA, C
    end subroutine

    ! conjugated transpose
    subroutine pztranc(M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC)
      use kind_oct_m
      implicit none
      integer :: M, N, IA, JA, DESCA, IC, JC, DESCC
      complex(r8) :: ALPHA, A, BETA, C
    end subroutine
  end interface pblas_tran

end module pblas_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
