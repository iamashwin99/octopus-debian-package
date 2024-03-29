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

module lalg_adv_oct_m
  use blacs_proc_grid_oct_m
  use blacs_oct_m
  use debug_oct_m
  use global_oct_m
  use lapack_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use scalapack_oct_m
  use sort_oct_m
  use utils_oct_m
#ifdef HAVE_ELPA
  use elpa
#endif

  implicit none

  private
  public ::                       &
    lalg_cholesky,                &
    lalg_geneigensolve,           &
    lalg_eigensolve,              &
    lalg_eigensolve_nonh,         &
    lalg_eigensolve_parallel,     &
    lalg_determinant,             &
    lalg_inverter,                &
    lalg_sym_inverter,            &
    lalg_linsyssolve,             &
    lalg_singular_value_decomp,   &
    lalg_svd_inverse,             &
    lalg_invert_upper_triangular, &
    lalg_lowest_geneigensolve,    &
    lalg_lowest_eigensolve,       &
    zlalg_exp,                    &
    zlalg_phi,                    &
    lalg_zpseudoinverse,          &
    lalg_zeigenderivatives,       &
    lalg_check_zeigenderivatives, &
    lalg_zdni,                    &
    lalg_zduialpha,               &
    lalg_zd2ni,                   &
    lalg_least_squares,           &
    lalg_matrix_norm2


  type(profile_t), save :: cholesky_prof, eigensolver_prof

  interface lalg_cholesky
    module procedure dcholesky, zcholesky
  end interface lalg_cholesky

  interface lalg_geneigensolve
    module procedure dgeneigensolve, zgeneigensolve
  end interface lalg_geneigensolve

  interface lalg_eigensolve_nonh
    module procedure zeigensolve_nonh, deigensolve_nonh
  end interface lalg_eigensolve_nonh

  interface lalg_eigensolve
    module procedure deigensolve, zeigensolve
  end interface lalg_eigensolve

  interface lalg_eigensolve_parallel
    module procedure deigensolve_parallel, zeigensolve_parallel
  end interface lalg_eigensolve_parallel

  !> Note that lalg_determinant and lalg_inverter are just wrappers
  !! over the same routine.

  interface lalg_determinant
    module procedure ddeterminant, zdeterminant
  end interface lalg_determinant

  interface lalg_inverter
    module procedure dinverter, zinverter
  end interface lalg_inverter

  interface lalg_sym_inverter
    module procedure dsym_inverter, zsym_inverter
  end interface lalg_sym_inverter

  interface lalg_linsyssolve
    module procedure dlinsyssolve, zlinsyssolve
  end interface lalg_linsyssolve

  interface lalg_singular_value_decomp
    module procedure dsingular_value_decomp, zsingular_value_decomp
  end interface lalg_singular_value_decomp

  interface lalg_svd_inverse
    module procedure dsvd_inverse, zsvd_inverse
  end interface lalg_svd_inverse

  interface lalg_invert_upper_triangular
    module procedure dinvert_upper_triangular, zinvert_upper_triangular
  end interface lalg_invert_upper_triangular

  interface lalg_lowest_geneigensolve
    module procedure dlowest_geneigensolve, zlowest_geneigensolve
  end interface lalg_lowest_geneigensolve

  interface lalg_lowest_eigensolve
    module procedure dlowest_eigensolve, zlowest_eigensolve
  end interface lalg_lowest_eigensolve

  interface lapack_geev
    module procedure lalg_dgeev, lalg_zgeev
  end interface lapack_geev

  interface lalg_least_squares
    module procedure dleast_squares_vec, zleast_squares_vec
  end interface lalg_least_squares

  interface lalg_matrix_norm2
    module procedure dmatrix_norm2, zmatrix_norm2
  end interface lalg_matrix_norm2

contains

  ! ---------------------------------------------------------
  !> Auxiliary function
  FLOAT function sfmin()
    interface
      FLOAT function dlamch(cmach)
        implicit none
        character(1), intent(in) :: cmach
      end function dlamch
    end interface

    sfmin = dlamch('S')
  end function sfmin

  subroutine lalg_dgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    character(1), intent(in)    :: jobvl, jobvr
    integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
    real(r8),      intent(inout) :: a(:,:) !< a(lda,n)
    complex(r8),   intent(out)   :: w(:) !< w(n)
    real(r8),      intent(out)   :: vl(:,:), vr(:,:) !< vl(ldvl,n), vl(ldvr,n)
    real(r8),      intent(out)   :: rwork(:) !< rwork(max(1,2n))
    real(r8),      intent(out)   :: work(:)  !< work(lwork)
    integer,      intent(out)   :: info

    FLOAT, allocatable :: wr(:), wi(:)

    PUSH_SUB(lalg_dgeev)

    SAFE_ALLOCATE(wr(1:n))
    SAFE_ALLOCATE(wi(1:n))

    call dgeev(jobvl, jobvr, n, a(1, 1), lda, wr(1), wi(1), vl(1, 1), ldvl, vr(1, 1), ldvr, work(1), lwork, rwork(1), info)

    w(1:n) = TOCMPLX(wr(1:n), wi(1:n))

    SAFE_DEALLOCATE_A(wr)
    SAFE_DEALLOCATE_A(wi)

    POP_SUB(lalg_dgeev)
  end subroutine lalg_dgeev

  subroutine lalg_zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    character(1), intent(in)    :: jobvl, jobvr
    integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
    complex(r8),   intent(inout) :: a(:,:) !< a(lda,n)
    complex(r8),   intent(out)   :: w(:) !< w(n)
    complex(r8),   intent(out)   :: vl(:,:), vr(:,:) !< vl(ldvl,n), vl(ldvr,n)
    real(r8),      intent(out)   :: rwork(:) !< rwork(max(1,2n))
    complex(r8),   intent(out)   :: work(:)  !< work(lwork)
    integer,      intent(out)   :: info

    PUSH_SUB(lalg_zgeev)

    call zgeev(jobvl, jobvr, n, a(1, 1), lda, w(1), vl(1, 1), ldvl, vr(1, 1), ldvr, work(1), lwork, rwork(1), info)

    POP_SUB(lalg_zgeev)
  end subroutine lalg_zgeev

!>-------------------------------------------------
  !!
  !! This routine calculates the exponential of a matrix by using an
  !! eigenvalue decomposition.
  !!
  !! For the hermitian case:
  !!
  !!   A = V D V^T => exp(A) = V exp(D) V^T
  !!
  !! and in general
  !!
  !!   A = V D V^-1 => exp(A) = V exp(D) V^-1
  !!
  !! This is slow but it is simple to implement, and for the moment it
  !! does not affect performance.
  !!
  !!---------------------------------------------
  subroutine zlalg_exp(nn, pp, aa, ex, hermitian)
    integer,           intent(in)      :: nn
    CMPLX,             intent(in)      :: pp
    CMPLX,             intent(in)      :: aa(:, :)
    CMPLX,             intent(inout)   :: ex(:, :)
    logical,           intent(in)      :: hermitian

    CMPLX, allocatable :: evectors(:, :), zevalues(:)
    FLOAT, allocatable :: evalues(:)

    integer :: ii

    PUSH_SUB(zlalg_exp)

    SAFE_ALLOCATE(evectors(1:nn, 1:nn))

    if (hermitian) then
      SAFE_ALLOCATE(evalues(1:nn))
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors(1:nn, 1:nn) = aa(1:nn, 1:nn)

      call lalg_eigensolve(nn, evectors, evalues)

      zevalues(1:nn) = exp(pp*evalues(1:nn))

      do ii = 1, nn
        ex(1:nn, ii) = zevalues(1:nn)*conjg(evectors(ii, 1:nn))
      end do

      ex(1:nn, 1:nn) = matmul(evectors(1:nn, 1:nn), ex(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(evalues)
      SAFE_DEALLOCATE_A(zevalues)
    else
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors(1:nn, 1:nn) = aa(1:nn, 1:nn)

      call lalg_eigensolve_nonh(nn, evectors, zevalues)

      zevalues(1:nn) = exp(pp*zevalues(1:nn))

      ex(1:nn, 1:nn) = evectors(1:nn, 1:nn)

      call lalg_inverter(nn, evectors)

      do ii = 1, nn
        evectors(1:nn, ii) = zevalues(1:nn)*evectors(1:nn, ii)
      end do

      ex(1:nn, 1:nn) = matmul(ex(1:nn, 1:nn), evectors(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(zevalues)
    end if

    SAFE_DEALLOCATE_A(evectors)

    POP_SUB(zlalg_exp)
  end subroutine zlalg_exp


  !>-------------------------------------------------
  !!
  !! This routine calculates phi(pp*A), where A is a matrix,
  !! pp is any complex number, and phi is the function:
  !!
  !! phi(x) = (e^x - 1)/x
  !!
  !! For the Hermitian case, for any function f:
  !!
  !!   A = V D V^T => f(A) = V f(D) V^T
  !!
  !! and in general
  !!
  !!   A = V D V^-1 => f(A) = V f(D) V^-1
  !!
  !!---------------------------------------------
  subroutine zlalg_phi(nn, pp, aa, ex, hermitian)
    integer,           intent(in)      :: nn
    CMPLX,             intent(in)      :: pp
    CMPLX,             intent(in)      :: aa(:, :)
    CMPLX,             intent(inout)   :: ex(:, :)
    logical,           intent(in)      :: hermitian

    CMPLX, allocatable :: evectors(:, :), zevalues(:)
    FLOAT, allocatable :: evalues(:)

    integer :: ii

    PUSH_SUB(zlalg_phi)

    SAFE_ALLOCATE(evectors(1:nn, 1:nn))

    if (hermitian) then
      SAFE_ALLOCATE(evalues(1:nn))
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors = aa

      call lalg_eigensolve(nn, evectors, evalues)

      do ii = 1, nn
        zevalues(ii) = (exp(pp*evalues(ii)) - M_z1) / (pp*evalues(ii))
      end do

      do ii = 1, nn
        ex(1:nn, ii) = zevalues(1:nn)*conjg(evectors(ii, 1:nn))
      end do

      ex(1:nn, 1:nn) = matmul(evectors(1:nn, 1:nn), ex(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(evalues)
      SAFE_DEALLOCATE_A(zevalues)
    else
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors(1:nn, 1:nn) = aa(1:nn, 1:nn)

      call lalg_eigensolve_nonh(nn, evectors, zevalues)

      do ii = 1, nn
        zevalues(ii) = (exp(pp*zevalues(ii)) - M_z1) / (pp*zevalues(ii))
      end do

      ex(1:nn, 1:nn) = evectors(1:nn, 1:nn)

      call lalg_inverter(nn, evectors)

      do ii = 1, nn
        evectors(1:nn, ii) = zevalues(1:nn)*evectors(1:nn, ii)
      end do

      ex(1:nn, 1:nn) = matmul(ex(1:nn, 1:nn), evectors(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(zevalues)
    end if

    POP_SUB(zlalg_phi)
  end subroutine zlalg_phi


  !>-------------------------------------------------
  !! Computes the necessary ingredients to obtain,
  !! later, the first and second derivatives of the
  !! eigenvalues of a Hermitean complex matrix zmat,
  !! and the first derivatives of the eigenvectors.
  !!
  !! This follows the scheme of J. R. Magnus,
  !! Econometric Theory 1, 179 (1985), restricted to
  !! Hermitean matrices, although probably this can be
  !! found in other sources.
  !!---------------------------------------------
  subroutine lalg_zeigenderivatives(n, mat, zeigenvec, zeigenval, zmat)
    integer, intent(in) :: n
    CMPLX, intent(in) :: mat(:, :)
    CMPLX, intent(out) :: zeigenvec(:, :)
    CMPLX, intent(out) :: zeigenval(:)
    CMPLX, intent(out) :: zmat(:, :, :)

    integer :: i, alpha, beta
    CMPLX, allocatable :: knaught(:, :)
    CMPLX, allocatable :: lambdaminusdm(:, :)
    CMPLX, allocatable :: ilambdaminusdm(:, :)
    CMPLX, allocatable :: unit(:, :)

    PUSH_SUB(lalg_zeigenderivatives)

    SAFE_ALLOCATE(unit(1:n, 1:n))
    SAFE_ALLOCATE(knaught(1:n, 1:n))
    SAFE_ALLOCATE(lambdaminusdm(1:n, 1:n))
    SAFE_ALLOCATE(ilambdaminusdm(1:n, 1:n))

    zeigenvec = mat
    call lalg_eigensolve_nonh(n, zeigenvec, zeigenval)

    unit = M_z0
    do i = 1, n
      unit(i, i) = M_z1
    end do

    do i = 1, n
      do alpha = 1, n
        do beta = 1, n
          knaught(alpha, beta) = - zeigenvec(alpha, i)*conjg(zeigenvec(beta, i))
        end do
      end do
      knaught = knaught + unit
      lambdaminusdm = zeigenval(i)*unit - mat
      call lalg_zpseudoinverse(n, lambdaminusdm, ilambdaminusdm)
      zmat(:, :, i) = matmul(ilambdaminusdm, knaught)
    end do

    SAFE_DEALLOCATE_A(unit)
    SAFE_DEALLOCATE_A(knaught)
    SAFE_DEALLOCATE_A(lambdaminusdm)
    SAFE_DEALLOCATE_A(ilambdaminusdm)
    POP_SUB(lalg_zeigenderivatives)
  end subroutine lalg_zeigenderivatives


  !>-------------------------------------------------
  !! Computes the Moore-Penrose pseudoinverse of a
  !! complex matrix.
  !!---------------------------------------------
  subroutine lalg_zpseudoinverse(n, mat, imat)
    integer, intent(in) :: n
    CMPLX, intent(in) :: mat(:, :)
    CMPLX, intent(out) :: imat(:, :)

    integer :: i
    CMPLX, allocatable :: u(:, :), vt(:, :), sigma(:, :)
    FLOAT, allocatable :: sg_values(:)

    PUSH_SUB(lalg_zpseudoinverse)

    SAFE_ALLOCATE(u(1:n, 1:n))
    SAFE_ALLOCATE(vt(1:n, 1:n))
    SAFE_ALLOCATE(sigma(1:n, 1:n))
    SAFE_ALLOCATE(sg_values(1:n))

    imat = mat
    call lalg_singular_value_decomp(n, n, imat, u, vt, sg_values)

    sigma = M_z0
    do i = 1, n
      if (abs(sg_values(i)) <= CNST(1.0e-12) * maxval(abs(sg_values))) then
        sigma(i, i) = M_z0
      else
        sigma(i, i) = M_z1 / sg_values(i)
      end if
    end do

    vt = conjg(transpose(vt))
    u = conjg(transpose(u))
    imat = matmul(vt, matmul(sigma, u))

    ! Check if we truly have a pseudoinverse
    vt = matmul(mat, matmul(imat, mat)) - mat
    if (maxval(abs(vt)) > CNST(1.0e-10) * maxval(abs(mat))) then
      write(*, *) maxval(abs(vt))
      write(*, *) vt
      write(*, *)
      write(*, *) CNST(1.0e-10) * maxval(abs(mat))
      write(*, *) maxval(abs(vt)) > CNST(1.0e-10) * maxval(abs(mat))
      write(*, *) mat
      write(message(1), '(a)') 'Pseudoinverse failed.'
      call messages_fatal(1)
    end if

    SAFE_DEALLOCATE_A(u)
    SAFE_DEALLOCATE_A(vt)
    SAFE_DEALLOCATE_A(sigma)
    SAFE_DEALLOCATE_A(sg_values)
    POP_SUB(lalg_zpseudoinverse)
  end subroutine lalg_zpseudoinverse


  !>-------------------------------------------------
  !! The purpose of this routine is to check that "lalg_zeigenderivatives"
  !! is working properly, and therefore, it is not really called anywhere
  !! in the code. It is here only for debugging purposes (perhaps it will
  !! disappear in the future...)
  !!-------------------------------------------------
  subroutine lalg_check_zeigenderivatives(n, mat)
    integer, intent(in) :: n
    CMPLX, intent(in) :: mat(:, :)

    integer :: alpha, beta, gamma, delta
    CMPLX :: deltah, zuder_direct, zder_direct, zuder_directplus, zuder_directminus
    CMPLX, allocatable :: zeigenvec(:, :), dm(:, :), zeigref_(:, :), zeigenval(:), mmatrix(:, :, :)
    CMPLX, allocatable :: zeigplus(:), zeigminus(:), zeig0(:), zeigplusminus(:), zeigminusplus(:)

    PUSH_SUB(lalg_check_zeigenderivatives)

    SAFE_ALLOCATE(zeigenvec(1:n, 1:n))
    SAFE_ALLOCATE(dm(1:n, 1:n))
    SAFE_ALLOCATE(zeigref_(1:n, 1:n))
    SAFE_ALLOCATE(zeigenval(1:n))
    SAFE_ALLOCATE(mmatrix(1:n, 1:n, 1:n))
    SAFE_ALLOCATE(zeigplus(1:n))
    SAFE_ALLOCATE(zeigminus(1:n))
    SAFE_ALLOCATE(zeig0(1:n))
    SAFE_ALLOCATE(zeigplusminus(1:n))
    SAFE_ALLOCATE(zeigminusplus(1:n))

    ASSERT(n == 2)

    dm = mat
    call lalg_zeigenderivatives(2, dm, zeigref_, zeigenval, mmatrix)


    deltah = (CNST(0.000001), CNST(0.000))
    !deltah = M_z1 * maxval(abs(dm)) * CNST(0.001)
    do alpha = 1, n
      do beta = 1, n
        zder_direct = lalg_zdni(zeigref_(:, 2), alpha, beta)
        zuder_direct = lalg_zduialpha(zeigref_(:, 2), mmatrix(:, :, 2), 2, alpha, beta)

        zeigenvec = dm
        zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
        call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplus)
        zeigenvec(:, 1) = zeigenvec(:, 1)/sum(conjg(zeigref_(1:2, 1))*zeigenvec(1:2, 1))
        zeigenvec(:, 2) = zeigenvec(:, 2)/sum(conjg(zeigref_(1:2, 2))*zeigenvec(1:2, 2))
        zuder_directplus = zeigenvec(2, 2)

        zeigenvec = dm
        zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
        call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminus)
        zeigenvec(:, 1) = zeigenvec(:, 1)/sum(conjg(zeigref_(1:2, 1))*zeigenvec(1:2, 1))
        zeigenvec(:, 2) = zeigenvec(:, 2)/sum(conjg(zeigref_(1:2, 2))*zeigenvec(1:2, 2))
        zuder_directminus = zeigenvec(2, 2)

        write(*, '(2i1,4f24.12)') alpha, beta, zder_direct, (zeigplus(2) - zeigminus(2))/(M_TWO * deltah)
        write(*, '(2i1,4f24.12)') alpha, beta, &
          zuder_direct, (zuder_directplus - zuder_directminus) / (M_TWO * deltah)

        do gamma = 1, n
          do delta = 1, n
            if (alpha == gamma .and. beta == delta) then
              zder_direct = lalg_zd2ni(zeigref_(:, 1), mmatrix(:, :, 1), alpha, beta, gamma, delta)

              zeigenvec = dm
              call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeig0)

              zeigenvec = dm
              zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
              call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplus)

              zeigenvec = dm
              zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
              call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminus)

              write(*, '(4i1,4f24.12)') alpha, beta, gamma, delta, &
                zder_direct, (zeigplus(1) + zeigminus(1) - M_TWO*zeig0(1))/(deltah**2)
            else
              zder_direct = lalg_zd2ni(zeigref_(:, 1), mmatrix(:, :, 1), alpha, beta, gamma, delta)

              zeigenvec = dm
              zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
              zeigenvec(gamma, delta) = zeigenvec(gamma, delta) + deltah
              call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplus)

              zeigenvec = dm
              zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
              zeigenvec(gamma, delta) = zeigenvec(gamma, delta) - deltah
              call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminus)

              zeigenvec = dm
              zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
              zeigenvec(gamma, delta) = zeigenvec(gamma, delta) - deltah
              call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplusminus)

              zeigenvec = dm
              zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
              zeigenvec(gamma, delta) = zeigenvec(gamma, delta) + deltah
              call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminusplus)

              write(*, '(4i1,4f24.12)') alpha, beta, gamma, delta, &
                zder_direct, (zeigplus(1) + zeigminus(1) - zeigplusminus(1) - zeigminusplus(1))/(M_FOUR*deltah**2)


            end if
          end do
        end do

      end do
    end do

    SAFE_DEALLOCATE_A(zeigenval)
    SAFE_DEALLOCATE_A(mmatrix)
    SAFE_DEALLOCATE_A(zeigenvec)
    SAFE_DEALLOCATE_A(dm)
    SAFE_DEALLOCATE_A(zeigref_)
    SAFE_DEALLOCATE_A(zeigplus)
    SAFE_DEALLOCATE_A(zeigminus)
    SAFE_DEALLOCATE_A(zeig0)
    SAFE_DEALLOCATE_A(zeigplusminus)
    SAFE_DEALLOCATE_A(zeigminusplus)
    POP_SUB(lalg_check_zeigenderivatives)
  end subroutine lalg_check_zeigenderivatives

  CMPLX function lalg_zdni(eigenvec, alpha, beta)
    integer, intent(in) :: alpha, beta
    CMPLX :: eigenvec(2)
    lalg_zdni = conjg(eigenvec(alpha)) * eigenvec(beta)
  end function lalg_zdni

  CMPLX function lalg_zduialpha(eigenvec, mmatrix, alpha, gamma, delta)
    integer, intent(in) :: alpha, gamma, delta
    CMPLX :: eigenvec(2), mmatrix(2, 2)
    lalg_zduialpha = mmatrix(alpha, gamma) * eigenvec(delta)
  end function lalg_zduialpha

  CMPLX function lalg_zd2ni(eigenvec, mmatrix, alpha, beta, gamma, delta)
    integer, intent(in) :: alpha, beta, gamma, delta
    CMPLX :: eigenvec(2), mmatrix(2, 2)
    lalg_zd2ni  = conjg(mmatrix(alpha, delta) * eigenvec(gamma)) * eigenvec(beta) + &
      conjg(eigenvec(alpha)) * mmatrix(beta, gamma) * eigenvec(delta)
  end function lalg_zd2ni

#include "undef.F90"
#include "complex.F90"
#include "lalg_adv_lapack_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "lalg_adv_lapack_inc.F90"

end module lalg_adv_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
