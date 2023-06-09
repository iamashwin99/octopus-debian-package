!! Copyright (C) 2011 D. Strubbe
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
!> This module contains interfaces for ScaLAPACK routines
!! Interfaces are from http://www.netlib.org/scalapack/tools, double, complex16
!!
!! \note
!!  Each global data object is described by an associated description
!!  vector.  This vector stores the information required to establish
!!  the mapping between an object element and its corresponding process
!!  and memory location.
!!
!! \par
!!  Let A be a generic term for any 2D block cyclicly distributed array.
!!  Such a global array has an associated description vector DESCA.
!!  In the following comments, the character _ should be read as
!!  "of the global array".
! -----------------------------------------------------------------------

module scalapack_oct_m
  implicit none
#ifdef HAVE_SCALAPACK

  interface
    integer function iceil(inum, idenom)
      implicit none
      integer, intent(in) :: inum
      integer, intent(in) :: idenom
    end function iceil
  end interface

  interface
    subroutine descinit(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld, info)
      implicit none
      integer, intent(in)  :: desc
      integer, intent(in)  :: m
      integer, intent(in)  :: n
      integer, intent(in)  :: mb
      integer, intent(in)  :: nb
      integer, intent(in)  :: irsrc
      integer, intent(in)  :: icsrc
      integer, intent(in)  :: ictxt
      integer, intent(in)  :: lld
      integer, intent(out) :: info
    end subroutine descinit
  end interface

  interface
    subroutine infog2l(grindx, gcindx, desc, nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc)
      implicit none
      integer, intent(in)  :: grindx
      integer, intent(in)  :: gcindx
      integer, intent(in)  :: desc
      integer, intent(in)  :: nprow
      integer, intent(in)  :: npcol
      integer, intent(in)  :: myrow
      integer, intent(in)  :: mycol
      integer, intent(out) :: lrindx
      integer, intent(out) :: lcindx
      integer, intent(out) :: rsrc
      integer, intent(out) :: csrc
    end subroutine infog2l
  end interface

  !>  Computes a QR factorization of a real distributed \f$ m \times n\f$
  !! \f[
  !!  matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = Q R
  !! \f]
  interface scalapack_geqrf
    subroutine pdgeqrf(m, n, a, ia, ja, desca, tau, work, lwork, info)
      implicit none
      integer            ia, info, ja, lwork, m, n
      integer            desca
      double precision   a, tau, work
    end subroutine pdgeqrf

    subroutine pzgeqrf(m, n, a, ia, ja, desca, tau, work, lwork, info)
      use kind_oct_m
      implicit none
      integer            ia, info, ja, lwork, m, n
      integer            desca
      complex(r8)         a, tau, work
    end subroutine pzgeqrf
  end interface scalapack_geqrf

  !>  Generates an \f$ m \times n\f$ real distributed matrix Q denoting
  !!  A(IA:IA+M-1,JA:JA+N-1) with orthonormal columns, which is defined as
  !!  the first N columns of a product of K elementary reflectors of order
  !!  M
  !! \f[
  !!        Q  =  H(1) H(2) . . . H(k)
  !! \f]
  !!
  !!  as returned by PDGEQRF.
  interface scalapack_orgqr
    subroutine pdorgqr(m, n, k, a, ia, ja, desca, tau, work, lwork, info)
      implicit none
      integer            ia, info, ja, k, lwork, m, n
      integer            desca
      double precision   a, tau, work
    end subroutine pdorgqr

    subroutine pzungqr(m, n, k, a, ia, ja, desca, tau, work, lwork, info)
      use kind_oct_m
      implicit none
      integer            ia, info, ja, k, lwork, m, n
      integer            desca
      complex(r8)         a, tau, work
    end subroutine pzungqr
  end interface scalapack_orgqr

  interface
    subroutine pdgesv(n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info)
      implicit none
      integer            ia, ib, info, ja, jb, n, nrhs
      integer            desca, descb, ipiv
      double precision   a, b
    end subroutine pdgesv
  end interface

  interface
    subroutine pzgesv(n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info)
      use kind_oct_m
      implicit none
      integer            ia, ib, info, ja, jb, n, nrhs
      integer            desca, descb, ipiv
      complex(r8)         a, b
    end subroutine pzgesv
  end interface

  !>  Computes all eigenvalues and, optionally, eigenvectors
  !!  of a real symmetric matrix A by calling the recommended sequence
  !!  of ScaLAPACK routines.
  !!
  !!  In its present form, PDSYEV assumes a homogeneous system and makes
  !!  no checks for consistency of the eigenvalues or eigenvectors across
  !!  the different processes.  Because of this, it is possible that a
  !!  heterogeneous system may return incorrect results without any error
  !!  messages.
  !!
  !! \note
  !!  A description vector is associated with each 2D block-cyclicly dis-
  !!  tributed matrix.  This vector stores the information required to
  !!  establish the mapping between a matrix entry and its corresponding
  !!  process and memory location.
  !!
  !! \par
  !!  In the following comments, the character _ should be read as
  !!  "of the distributed matrix".  Let A be a generic term for any 2D
  !!  block cyclicly distributed matrix.  Its description vector is DESCA:
  interface scalapack_syev
    subroutine pdsyev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info)
      use kind_oct_m
      implicit none
      character,        intent(in)    :: jobz
      character,        intent(in)    :: uplo
      integer,          intent(in)    :: n
      real(r8),          intent(inout) :: a
      integer,          intent(in)    :: ia
      integer,          intent(in)    :: ja
      integer,          intent(in)    :: desca
      real(r8),          intent(out)   :: w
      real(r8),          intent(out)   :: z
      integer,          intent(in)    :: iz
      integer,          intent(in)    :: jz
      integer,          intent(in)    :: descz
      real(r8),          intent(out)   :: work
      integer,          intent(in)    :: lwork
      integer,          intent(out)   :: info
    end subroutine pdsyev

    subroutine pzheev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, rwork, lrwork, info)
      use kind_oct_m
      implicit none
      character,        intent(in)    :: jobz
      character,        intent(in)    :: uplo
      integer,          intent(in)    :: n
      complex(r8),       intent(inout) :: a
      integer,          intent(in)    :: ia
      integer,          intent(in)    :: ja
      integer,          intent(in)    :: desca
      real(r8),          intent(out)   :: w
      complex(r8),       intent(out)   :: z
      integer,          intent(in)    :: iz
      integer,          intent(in)    :: jz
      integer,          intent(in)    :: descz
      complex(r8),       intent(out)   :: work
      integer,          intent(in)    :: lwork
      complex(r8),       intent(out)   :: rwork
      integer,          intent(in)    :: lrwork
      integer,          intent(out)   :: info
    end subroutine pzheev
  end interface scalapack_syev

  !>  Computes selected eigenvalues and, optionally, eigenvectors
  !!  of a real symmetric matrix A by calling the recommended sequence
  !!  of ScaLAPACK routines.  Eigenvalues/vectors can be selected by
  !!  specifying a range of values or a range of indices for the desired
  !!  eigenvalues.
  interface scalapack_syevx
    subroutine pdsyevx( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, &
      m, nz, w, orfac, z, iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info )
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: jobz
      character(1), intent(in)    :: range
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      real(r8),      intent(inout) :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      real(r8),      intent(in)    :: vl
      real(r8),      intent(in)    :: vu
      integer,      intent(in)    :: il
      integer,      intent(in)    :: iu
      real(r8) ,     intent(in)    :: abstol
      integer,      intent(out)   :: m
      integer,      intent(out)   :: nz
      real(r8),      intent(out)   :: w
      real(r8),      intent(in)    :: orfac
      real(r8),      intent(out)   :: z
      integer,      intent(in)    :: iz
      integer,      intent(in)    :: jz
      integer,      intent(in)    :: descz
      real(r8),      intent(out)   :: work
      integer,      intent(in)    :: lwork
      integer,      intent(inout) :: iwork
      integer,      intent(in)    :: liwork
      integer,      intent(out)   :: ifail
      integer,      intent(out)   :: iclustr
      real(r8),      intent(out)   :: gap
      integer,      intent(out)   :: info

    end subroutine pdsyevx

    subroutine pzheevx(jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, &
      jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: jobz
      character(1), intent(in)    :: range
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      complex(r8),   intent(inout) :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      real(r8),      intent(in)    :: vl
      real(r8),      intent(in)    :: vu
      integer,      intent(in)    :: il
      integer,      intent(in)    :: iu
      real(r8),      intent(in)    :: abstol
      integer,      intent(out)   :: m
      integer,      intent(out)   :: nz
      real(r8),      intent(out)   :: w
      real(r8),      intent(in)    :: orfac
      complex(r8),   intent(out)   :: z
      integer,      intent(in)    :: iz
      integer,      intent(in)    :: jz
      integer,      intent(in)    :: descz
      complex(r8),   intent(out)   :: work
      integer,      intent(in)    :: lwork
      real(r8),      intent(out)   :: rwork
      integer,      intent(in)    :: lrwork
      integer,      intent(inout) :: iwork
      integer,      intent(in)    :: liwork
      integer,      intent(out)   :: ifail
      integer,      intent(out)   :: iclustr
      real(r8),      intent(out)   :: gap
      integer,      intent(out)   :: info
    end subroutine pzheevx
  end interface scalapack_syevx

  !>  Computes all the eigenvalues, and optionally,
  !!  the eigenvectors
  !!  of a real generalized SY-definite eigenproblem, of the form
  !!  \f$ sub( A ) x=(\lambda) sub( B ) x,  sub( A ) sub( B ) x=(\lambda) x, \mbox{ or }
  !!  sub( B ) sub( A ) x=(\lambda) x \f$.
  !!  Here sub(A) denoting A(IA:IA+N-1, JA:JA+N-1) is assumed to be
  !!  SY, and sub(B) denoting B(IB:IB+N-1, JB:JB+N-1) is assumed
  !!  to be symmetric positive definite.
  interface scalapack_sygvx
    subroutine pdsygvx(ibtype, jobz, range, uplo, n, a, ia, ja,       &
      desca, b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz,  &
      work, lwork, iwork, liwork, ifail, iclustr, gap, info)
      use kind_oct_m
      implicit none
      integer,             intent(in)    :: ibtype
      character,           intent(in)    :: jobz
      character,           intent(in)    :: range
      character,           intent(in)    :: uplo
      integer,             intent(in)    :: n
      real(r8),             intent(inout) :: a
      integer,             intent(in)    :: ia
      integer,             intent(in)    :: ja
      integer,             intent(in)    :: desca
      real(r8),             intent(inout) :: b
      integer,             intent(in)    :: ib
      integer,             intent(in)    :: jb
      integer,             intent(in)    :: descb
      real(r8),             intent(in)    :: vl
      real(r8),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(r8),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      integer,             intent(out)   :: nz
      real(r8),             intent(in)    :: w
      real(r8),             intent(in)    :: orfac
      real(r8),             intent(out)   :: z
      integer,             intent(in)    :: iz
      integer,             intent(in)    :: jz
      integer,             intent(in)    :: descz
      real(r8),             intent(out)   :: work
      integer,             intent(in)    :: lwork
      integer,             intent(out)   :: iwork
      integer,             intent(in)    :: liwork
      integer,             intent(out)   :: ifail
      integer,             intent(out)   :: iclustr
      real(r8),             intent(out)   :: gap
      integer,             intent(out)   :: info
    end subroutine pdsygvx
  end interface scalapack_sygvx

  ! -------------------------------------------------------------
  !>  Computes all the eigenvalues, and optionally,
  !!  the eigenvectors
  !!  of a complex generalized Hermitian-definite eigenproblem, of the form
  !!  \f$ sub( A ) x=(\lambda) sub( B ) x,  sub( A ) sub( B ) x=(\lambda) x, \mbox{ or }
  !!  sub( B ) sub( A ) x=(\lambda) x \f$.
  !!  Here sub(A) denoting A(IA:IA+N-1, JA:JA+N-1) is assumed to be
  !!  Hermitian, and sub(B) denoting B(IB:IB+N-1, JB:JB+N-1) is assumed
  !!  to be Hermitian positive definite.
  interface scalapack_hegvx
    subroutine pzhegvx(ibtype, jobz, range, uplo, n, a, ia, ja,       &
      desca, b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz,  &
      work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
      use kind_oct_m
      implicit none
      integer,             intent(in)    :: ibtype
      character,           intent(in)    :: jobz
      character,           intent(in)    :: range
      character,           intent(in)    :: uplo
      integer,             intent(in)    :: n
      complex(r8),          intent(inout) :: a
      integer,             intent(in)    :: ia
      integer,             intent(in)    :: ja
      integer,             intent(in)    :: desca
      complex(r8),          intent(inout) :: b
      integer,             intent(in)    :: ib
      integer,             intent(in)    :: jb
      integer,             intent(in)    :: descb
      real(r8),             intent(in)    :: vl
      real(r8),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(r8),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      integer,             intent(out)   :: nz
      real(r8),             intent(in)    :: w
      real(r8),             intent(in)    :: orfac
      complex(r8),          intent(out)   :: z
      integer,             intent(in)    :: iz
      integer,             intent(in)    :: jz
      integer,             intent(in)    :: descz
      complex(r8),          intent(out)   :: work
      integer,             intent(in)    :: lwork
      real(r8),             intent(out)   :: rwork
      integer,             intent(in)    :: lrwork
      integer,             intent(out)   :: iwork
      integer,             intent(in)    :: liwork
      integer,             intent(out)   :: ifail
      integer,             intent(out)   :: iclustr
      real(r8),             intent(out)   :: gap
      integer,             intent(out)   :: info
    end subroutine pzhegvx
  end interface scalapack_hegvx

  ! -------------------------------------------------------------
  !>  Computes the Cholesky factorization of an \f$ n \times n \f$ real
  !!  symmetric positive definite distributed matrix sub(A) denoting
  !!  A(IA:IA+N-1, JA:JA+N-1).
  !!
  !!  The factorization has the form
  !!
  !! \f[
  !!            sub( A ) = U` U , \mbox{ if UPLO }= 'U'
  !! \f]
  !! or
  !! \f[
  !!            sub( A ) = L L` , \mbox{ if UPLO }= 'L',
  !! \f]
  !!
  !!  where U is an upper triangular matrix and L is lower triangular.
  interface scalapack_potrf
    subroutine pdpotrf(uplo, n, a, ia, ja, desca, info)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      real(r8),      intent(inout) :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      integer,      intent(out)   :: info
    end subroutine pdpotrf

    subroutine pzpotrf(uplo, n, a, ia, ja, desca, info)
      use kind_oct_m
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      complex(r8),   intent(inout) :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      integer,      intent(out)   :: info
    end subroutine pzpotrf
  end interface scalapack_potrf

  interface
    subroutine pzlacp3(m, i, a, desca, b, ldb, ii, jj, rev)
      use kind_oct_m
      implicit none
      integer            i, ii, jj, ldb, m, rev
      integer            desca
      complex(r8)         a, b
    end subroutine pzlacp3
  end interface

  interface
    subroutine pdlacp3(m, i, a, desca, b, ldb, ii, jj, rev)
      use kind_oct_m
      implicit none
      integer            i, ii, jj, ldb, m, rev
      integer            desca
      real(r8)            a, b
    end subroutine pdlacp3
  end interface

  interface
    integer function indxl2g(indxloc, nb, iproc, isrcproc, nprocs)
      implicit none
      integer            indxloc, iproc, isrcproc, nb, nprocs
    end function
  end interface

  interface
    integer function indxg2l(indxglob, nb, iproc, isrcproc, nprocs)
      implicit none
      integer            indxglob, iproc, isrcproc, nb, nprocs
    end function
  end interface

  interface
    integer function indxg2p(indxglob, nb, iproc, isrcproc, nprocs)
      implicit none
      integer            indxglob, iproc, isrcproc, nb, nprocs
    end function
  end interface

#endif
end module scalapack_oct_m

!! local Variables:
!! mode: f90
!! coding: utf-8
!! End:
