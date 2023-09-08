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

!> @brief This module is intended to contain "only mathematical" functions and procedures.
module math_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use loct_math_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                     &
    is_close,                   &
    diagonal_matrix,            &
    ylmr_cmplx,                 &
    ylmr_real,                  &
    weights,                    &
    factorial,                  &
    hermite,                    &
    dcross_product,             &
    zcross_product,             &
    zdcross_product,            &
    dzcross_product,            &
    ddelta,                     &
    member,                     &
    make_idx_set,               &
    interpolation_coefficients, &
    interpolate,                &
    even,                       &
    odd,                        &
    cartesian2hyperspherical,   &
    hyperspherical2cartesian,   &
    hypersphere_grad_matrix,    &
    pad,                        &
    pad_pow2,                   &
    log2,                       &
    is_prime,                   &
    generate_rotation_matrix,   &
    numder_ridders,             &
    generalized_laguerre_polynomial


  interface is_close
    module procedure dis_close_scalar_0
  end interface

  interface diagonal_matrix
    module procedure idiagonal_matrix, ddiagonal_matrix, zdiagonal_matrix
  end interface diagonal_matrix

  ! ------------------------------------------------------------------------------
  !> This is the common interface to a simple-minded polynomical interpolation
  !! procedure (simple use of the classical formula of Lagrange).
  interface interpolate
    module procedure dinterpolate_0, dinterpolate_1, dinterpolate_2
    module procedure zinterpolate_0, zinterpolate_1, zinterpolate_2
  end interface interpolate
  ! ------------------------------------------------------------------------------

  interface log2
    module procedure dlog2, ilog2, llog2
  end interface log2

  interface pad
    module procedure pad4, pad8, pad48, pad88
  end interface pad

contains

  !> @brief Are \f$x\f$ and \f$y\f$ equal within a tolerance.
  !!
  !! The function evaluates the expression:
  !! \f[
  !!     |x - y| \leq (atol * rtol) * |y|
  !! \f]
  !!
  !! The tolerance values are positive, typically very small numbers.
  !! The relative difference \f$(rtol * |y|)\f$ and the absolute difference
  !! \f$ atol \f$ are added together to compare against the absolute difference
  !! between \f$x\f$ and \f$y\f$. Default tolerances are based on numpy''s
  !! [implementation](https://numpy.org/doc/stable/reference/generated/numpy.isclose.html).
  !!
  !! @param[in] x           Scalar.
  !! @param[in] y           Scalar.
  !! @param[in] rtol        Optional, relative tolerance.
  !! @param[in] atol        Optional, absolute tolerance.
  !! @return    is_close   .true. if \f$x\f$ and \f$y\f$ are close.
  pure logical function dis_close_scalar_0(x, y, rtol, atol)
    FLOAT, intent(in) :: x, y
    FLOAT, optional, intent(in) :: rtol
    FLOAT, optional, intent(in) :: atol
    FLOAT :: atol_, rtol_

    atol_ = optional_default(atol, CNST(1.e-8))
    rtol_ = optional_default(rtol, CNST(1.e-5))

    dis_close_scalar_0 = abs(x - y) <= (atol_ + rtol_ * abs(y))
  end function

  ! ---------------------------------------------------------
  !> Currently only returns a matrix whose diagonal elements are all the
  !! same. Note that the real and complex versions are in math_inc.F90.
  pure function idiagonal_matrix(dim, diag) result(matrix)
    integer, intent(in) :: dim
    integer, intent(in) :: diag
    integer :: matrix(dim, dim)

    integer :: ii

    matrix = 0
    do ii = 1, dim
      matrix(ii, ii) = diag
    end do

  end function idiagonal_matrix

  ! ---------------------------------------------------------
  recursive function hermite(n, x) result (h)
    integer, intent(in) :: n
    FLOAT,   intent(in) :: x

    FLOAT :: h

    ! no push_sub, called too frequently

    if (n <= 0) then
      h = M_ONE
    elseif (n == 1) then
      h = M_TWO*x
    else
      h = M_TWO*x*hermite(n-1,x) - M_TWO*(n-1)*hermite(n-2,x)
    end if

  end function hermite


  ! ---------------------------------------------------------
  recursive function factorial (n) RESULT (fac)
    integer, intent(in) :: n

    integer :: fac

    ! no push_sub for recursive function

    if (n <= 1) then
      fac = 1
    else
      fac = n * factorial(n-1)
    end if
  end function factorial

  ! ---------------------------------------------------------
  !> Computes spherical harmonics ylm at position (x, y, z)
  subroutine ylmr_cmplx(xx, li, mi, ylm)
    FLOAT,   intent(in) :: xx(3)
    integer, intent(in) :: li, mi
    CMPLX,   intent(out) :: ylm

    integer :: i
    FLOAT :: dx, dy, dz, r, plm, cosm, sinm, cosmm1, sinmm1, cosphi, sinphi
    FLOAT,   parameter :: tiny = CNST(1.e-15)

    !   no push_sub: this routine is called too many times

    ! if l=0, no calculations are required
    if (li == 0) then
      ylm = TOCMPLX(sqrt(M_ONE/(M_FOUR*M_PI)), M_ZERO)
      return
    end if

    r = sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)

    ! if r=0, direction is undefined => make ylm=0 except for l=0
    if (r <= tiny) then
      ylm = M_z0
      return
    end if
    dx = xx(1)/r
    dy = xx(2)/r
    dz = xx(3)/r

    ! get the associated Legendre polynomial (including the normalization factor)
    plm = loct_legendre_sphplm(li, abs(mi), dz)

    ! compute sin(|m|*phi) and cos(|m|*phi)
    r = sqrt(dx*dx + dy*dy)
    if(r > tiny) then
      cosphi = dx/r
      sinphi = dy/r
    else ! In this case the values are ill defined so we choose \phi=0
      cosphi = M_ZERO
      sinphi = M_ZERO
    end if

    cosm = M_ONE
    sinm = M_ZERO
    do i = 1, abs(mi)
      cosmm1 = cosm
      sinmm1 = sinm
      cosm = cosmm1*cosphi - sinmm1*sinphi
      sinm = cosmm1*sinphi + sinmm1*cosphi
    end do

    !And now ylm
    ylm = plm*TOCMPLX(cosm, sinm)

    if (mi < 0) then
      ylm = conjg(ylm)
      do i = 1, abs(mi)
        ylm = -ylm
      end do
    end if

  end subroutine ylmr_cmplx


  ! ---------------------------------------------------------
  !> This is a Numerical Recipes-based subroutine
  !! computes real spherical harmonics ylm at position (x, y, z):
  !!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
  !!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
  !! with (theta,phi) the polar angles of r, c a positive normalization
  subroutine ylmr_real(xx, li, mi, ylm)
    FLOAT,           intent(in)  :: xx(3)
    integer,         intent(in)  :: li, mi
    FLOAT,           intent(out) :: ylm

    integer, parameter :: lmaxd = 20
    FLOAT,   parameter :: tiny = CNST(1.e-15)
    integer :: i, ilm0, l, m, mabs
    integer, save :: lmax = -1

    FLOAT :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac, &
      fourpi, plgndr, phase, pll, pmm, pmmp1, sinm, &
      sinmm1, sinphi, rsize, Rx, Ry, Rz, xysize
    FLOAT, save :: c(0:(lmaxd+1)*(lmaxd+1))

    ! no push_sub: called too frequently

    ! evaluate normalization constants once and for all
    if (li > lmax) then
      fourpi = M_FOUR*M_PI
      do l = 0, li
        ilm0 = l*l + l
        do m = 0, l
          fac = (2*l+1)/fourpi
          do i = l - m + 1, l + m
            fac = fac/i
          end do
          c(ilm0 + m) = sqrt(fac)
          ! next line because ylm are real combinations of m and -m
          if (m /= 0) c(ilm0 + m) = c(ilm0 + m)*sqrt(M_TWO)
          c(ilm0 - m) = c(ilm0 + m)
        end do
      end do
      lmax = li
    end if

    ! if l=0, no calculations are required
    if (li == 0) then
      ylm = c(0)
      return
    end if

    ! if r=0, direction is undefined => make ylm=0 except for l=0
    rsize = sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)
    if (rsize < tiny) then
      ylm = M_ZERO
      return
    end if

    Rx = xx(1)/rsize
    Ry = xx(2)/rsize
    Rz = xx(3)/rsize

    ! explicit formulas for l=1 and l=2
    if (li == 1) then
      select case (mi)
      case (-1)
        ylm = (-c(1))*Ry
      case (0)
        ylm = c(2)*Rz
      case (1)
        ylm = (-c(3))*Rx
      end select
      return
    end if

    if (li == 2) then
      select case (mi)
      case (-2)
        ylm = c(4)*CNST(6.0)*Rx*Ry
      case (-1)
        ylm = (-c(5))*M_THREE*Ry*Rz
      case (0)
        ylm = c(6)*M_HALF*(M_THREE*Rz*Rz - M_ONE)
      case (1)
        ylm = (-c(7))*M_THREE*Rx*Rz
      case (2)
        ylm = c(8)*M_THREE*(Rx*Rx - Ry*Ry)
      end select
      return
    end if

    ! general algorithm based on routine plgndr of numerical recipes
    mabs = abs(mi)
    xysize = sqrt(max(Rx*Rx + Ry*Ry, tiny))
    cosphi = Rx/xysize
    sinphi = Ry/xysize
    cosm = M_ONE
    sinm = M_ZERO
    do m = 1, mabs
      cosmm1 = cosm
      sinmm1 = sinm
      cosm = cosmm1*cosphi - sinmm1*sinphi
      sinm = cosmm1*sinphi + sinmm1*cosphi
    end do

    if (mi < 0) then
      phase = sinm
      dphase = mabs*cosm
    else
      phase = cosm
      dphase = (-mabs)*sinm
    end if

    pmm = M_ONE
    fac = M_ONE

    if (mabs > M_ZERO) then
      do i = 1, mabs
        pmm = (-pmm)*fac*xysize
        fac = fac + M_TWO
      end do
    end if

    if (li == mabs) then
      plgndr = pmm
      dplg = (-li)*Rz*pmm/(xysize**2)
    else
      pmmp1 = Rz*(2*mabs + 1)*pmm
      if (li == mabs + 1) then
        plgndr = pmmp1
        dplg = -((li*Rz*pmmp1 - (mabs + li)*pmm)/(xysize**2))
      else
        do l = mabs + 2, li
          pll = (Rz*(2*l - 1)*pmmp1 - (l + mabs - 1)*pmm)/(l - mabs)
          pmm = pmmp1
          pmmp1 = pll
        end do
        plgndr = pll
        dplg = -((li*Rz*pll - (l + mabs - 1)*pmm)/(xysize**2))
      end if
    end if

    ilm0 = li*li + li
    cmi = c(ilm0 + mi)
    ylm = cmi*plgndr*phase

  end subroutine ylmr_real


  ! ---------------------------------------------------------
  !> Compute the weights for finite-difference calculations:
  !!
  !!  N -> highest order of the derivative to be approximated
  !!  M -> number of grid points to be used in the approximation.
  !!
  !!  c(j,k,i) -> ith order derivative at kth-order approximation
  !!              j=0,k: the coefficients acting of each point
  !!
  !!  side -> -1 left-sided, +1 right-sided, 0 centered (default)
  subroutine weights(N, M, cc, side)
    integer,           intent(in)  :: N, M
    FLOAT,             intent(out) :: cc(0:,0:,0:) !< (0:M, 0:M, 0:N)
    integer, optional, intent(in)  :: side

    integer :: i, j, k, mn, side_
    FLOAT :: c1, c2, c3, c4, c5, xi
    FLOAT, allocatable :: x(:)

    PUSH_SUB(weights)

    SAFE_ALLOCATE(x(0:M))

    if (present(side)) then
      side_ = side
    else
      side_ = 0
    end if

    select case (side_)
    case (-1)
      ! grid-points for left-side finite-difference formulas on an equi-spaced grid
      mn = M
      x(:) = (/(-i,i=0,mn)/)
    case (+1)
      ! grid-points for right-side finite-difference formulas on an equi-spaced grid
      mn = M
      x(:) = (/(i,i=0,mn)/)
    case default
      ! grid-points for centered finite-difference formulas on an equi-spaced grid
      mn = M/2
      x(:) = (/0,(-i,i,i=1,mn)/)
    end select



    xi = M_ZERO  ! point at which the approx. are to be accurate

    cc = M_ZERO
    cc(0,0,0) = M_ONE

    c1 = M_ONE
    c4 = x(0) - xi

    do j = 1, M
      mn = min(j,N)
      c2 = M_ONE
      c5 = c4
      c4 = x(j) - xi

      do k = 0, j - 1
        c3 = x(j) - x(k)
        c2 = c2*c3

        if (j <= N) cc(k, j - 1, j) = M_ZERO
        cc(k, j, 0) = c4*cc(k, j - 1, 0)/c3

        do i = 1, mn
          cc(k, j, i) = (c4*cc(k, j - 1, i) - i*cc(k, j - 1, i - 1))/c3
        end do

        cc(j, j, 0) = -c1*c5*cc(j - 1, j - 1, 0) / c2
      end do

      do i = 1, mn
        cc(j, j, i) = c1*(i*cc(j - 1, j - 1, i - 1) - c5*cc(j - 1, j - 1, i))/c2
      end do

      c1 = c2
    end do

    SAFE_DEALLOCATE_A(x)
    POP_SUB(weights)
  end subroutine weights

  ! ---------------------------------------------------------
  FLOAT pure function ddelta(i, j)
    integer, intent(in) :: i
    integer, intent(in) :: j

    ! no push_sub in pure function

    if (i == j) then
      ddelta = M_ONE
    else
      ddelta = M_ZERO
    end if

  end function ddelta


  ! ---------------------------------------------------------
  !> Construct out(1:length) = (/1, ..., n/) if in is not present,
  !! out(1:length) = in otherwise.
  subroutine make_idx_set(n, out, length, in)
    integer,              intent(in)  :: n
    integer, allocatable, intent(out) :: out(:)
    integer,              intent(out) :: length
    integer, optional,    intent(in)  :: in(:)

    integer :: i

    PUSH_SUB(make_idx_set)

    if (present(in)) then
      length = ubound(in, 1)
      SAFE_ALLOCATE(out(1:length))
      out = in
    else
      length = n
      SAFE_ALLOCATE(out(1:length))
      do i = 1, length
        out(i) = i
      end do
    end if

    POP_SUB(make_idx_set)
  end subroutine make_idx_set


  ! ---------------------------------------------------------
  !> Considers a(1:ubound(a, 1)) as an integer set and checks
  !! if n is a member of it.
  logical function member(n, a)
    integer, intent(in) :: n
    integer, intent(in) :: a(:)

    integer :: i

    ! no push_sub, called too frequently

    member = .false.

    do i = 1, ubound(a, 1)
      if (a(i) == n) then
        member = .true.
        exit
      end if
    end do

  end function member


  ! ---------------------------------------------------------
  subroutine interpolation_coefficients(nn, xa, xx, cc)
    integer, intent(in)  :: nn    !< the number of points and coefficients
    FLOAT,   intent(in)  :: xa(:) !< the nn points where we know the function
    FLOAT,   intent(in)  :: xx    !< the point where we want the function
    FLOAT,   intent(out) :: cc(:) !< the coefficients

    integer :: ii, kk

    ! no push_sub, called too frequently

    do ii = 1, nn
      cc(ii) = M_ONE
      do kk = 1, nn
        if (kk == ii) cycle
        cc(ii) = cc(ii)*(xx - xa(kk))/(xa(ii) - xa(kk))
      end do
    end do

  end subroutine interpolation_coefficients


  ! ---------------------------------------------------------
  !> Returns if n is even.
  logical pure function even(n)
    integer, intent(in) :: n

    even = (mod(n, 2) == 0)

  end function even


  ! ---------------------------------------------------------
  !> Returns if n is odd.
  logical pure function odd(n)
    integer, intent(in) :: n

    odd = .not. even(n)

  end function odd

  ! ---------------------------------------------------------
  !> Performs a transformation of an n-dimensional vector
  !! from Cartesian coordinates to hyperspherical coordinates
  subroutine cartesian2hyperspherical(x, u)
    FLOAT, intent(in)  :: x(:)
    FLOAT, intent(out) :: u(:)

    integer :: n, k, j
    FLOAT :: sumx2
    FLOAT, allocatable :: xx(:)

    PUSH_SUB(cartesian2hyperspherical)

    n = size(x)
    ASSERT(n>1)
    ASSERT(size(u) == n-1)

    ! These lines make the code less machine-dependent.
    SAFE_ALLOCATE(xx(1:n))
    do j = 1, n
      if (abs(x(j))<CNST(1.0e-8)) then
        xx(j) = M_ZERO
      else
        xx(j) = x(j)
      end if
    end do

    u = M_ZERO
    do k = 1, n-1
      sumx2 = M_ZERO
      do j = k+1, n
        sumx2 = sumx2 + xx(j)**2
      end do
      if (abs(sumx2) <= M_EPSILON .and. abs(xx(k)) <= M_EPSILON) exit
      if (k < n - 1) then
        u(k) = atan2(sqrt(sumx2), xx(k))
      else
        u(n-1) = atan2(xx(n), xx(n-1))
      end if
    end do


    POP_SUB(cartesian2hyperspherical)
  end subroutine cartesian2hyperspherical
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Performs the inverse transformation of
  !! cartesian2hyperspherical
  subroutine hyperspherical2cartesian(u, x)
    FLOAT, intent(in)  :: u(:)
    FLOAT, intent(out) :: x(:)

    integer :: n, j, k

    PUSH_SUB(hyperspherical2cartesian)

    n = size(x)
    ASSERT(n>1)
    ASSERT(size(u) == n-1)

    if (n == 2) then
      x(1) = cos(u(1))
      x(2) = sin(u(1))
    elseif (n == 3) then
      x(1) = cos(u(1))
      x(2) = sin(u(1))*cos(u(2))
      x(3) = sin(u(1))*sin(u(2))
    else
      x(1) = cos(u(1))
      x(2) = sin(u(1))*cos(u(2))
      x(3) = sin(u(1))*sin(u(2))*cos(u(3))
      do j = 4, n - 1
        x(j) = M_ONE
        do k = 1, j - 1
          x(j) = x(j) * sin(u(k))
        end do
        x(j) = x(j) * cos(u(j))
      end do
      x(n) = M_ONE
      do k = 1, n - 2
        x(n) = x(n) * sin(u(k))
      end do
      x(n) = x(n) * sin(u(n-1))
    end if

    POP_SUB(hyperspherical2cartesian)
  end subroutine hyperspherical2cartesian
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Gives the hyperspherical gradient matrix, which contains
  !! the derivatives of the Cartesian coordinates with respect
  !! to the hyperspherical angles.
  subroutine  hypersphere_grad_matrix(grad_matrix, r, x)
    FLOAT, intent(out)  :: grad_matrix(:,:)
    FLOAT, intent(in)   :: r     !< radius of hypersphere
    FLOAT, intent(in)   :: x(:)  !< array of hyperspherical angles

    integer :: n, l, m

    PUSH_SUB(hypersphere_grad_matrix)

    n = size(x)+1  ! the dimension of the matrix is (n-1)x(n)

    ! --- l=1 ---
    grad_matrix = M_ONE
    grad_matrix(1,1) = -r*sin(x(1))
    do m = 2, n-1
      grad_matrix(m,1) = M_ZERO
    end do

    ! --- l=2..(n-1) ---
    do l=2, n-1
      do m=1, n-1
        if (m == l) then
          grad_matrix(m,l) = -r*grad_matrix(m,l)*product(sin(x(1:l)))
        elseif (m < l) then
          grad_matrix(m,l) = grad_matrix(m,l)*r*cos(x(m))*cos(x(l))*product(sin(x(1:m-1)))*product(sin(x(m+1:l-1)))
        else
          grad_matrix(m,l) = M_ZERO
        end if
      end do
    end do

    ! --- l=n ---
    do m=1, n-1
      grad_matrix(m,n) = r*cos(x(m))*grad_matrix(m,n)*product(sin(x(1:m-1)))*product(sin(x(m+1:n-1)))
    end do

    POP_SUB(hypersphere_grad_matrix)
  end subroutine  hypersphere_grad_matrix

  ! ---------------------------------------------------------
  integer(i8) pure function pad88(size, blk)
    integer(i8), intent(in) :: size
    integer(i8), intent(in) :: blk

    integer(i8) :: mm

    mm = mod(size, blk)
    if (mm == 0) then
      pad88 = size
    else
      pad88 = size + blk - mm
    end if
  end function pad88

  ! ---------------------------------------------------------
  integer(i8) pure function pad48(size, blk)
    integer,     intent(in) :: size
    integer(i8), intent(in) :: blk

    pad48 = pad88(int(size, i8), blk)
  end function pad48

  ! ---------------------------------------------------------
  integer(i8) pure function pad8(size, blk)
    integer(i8), intent(in) :: size
    integer, intent(in) :: blk

    pad8 = pad88(size, int(blk, i8))
  end function pad8

  ! ---------------------------------------------------------
  integer pure function pad4(size, blk)
    integer, intent(in) :: size
    integer, intent(in) :: blk

    pad4 = int(pad88(int(size, i8), int(blk, i8)), i4)
  end function pad4

  ! ---------------------------------------------------------

  integer pure function pad_pow2(size)
    integer, intent(in) :: size

    integer :: mm, mm2

    mm = size
    pad_pow2 = 1

    ! loop below never terminates otherwise! just pick 1 as value.
    if (size == 0) return

    ! first we divide by two all the times we can, so we catch exactly
    ! the case when size is already a power of two
    do
      mm2 = mm/2
      if (mm - 2*mm2 /= 0) exit
      pad_pow2 = pad_pow2*2
      mm = mm2
    end do

    ! the rest is handled by log2
    if (mm /= 1) then
      pad_pow2 = pad_pow2*2**(ceiling(log(TOFLOAT(mm))/log(2.0_8)))
    end if

  end function pad_pow2

  ! -------------------------------------------------------

  FLOAT pure function dlog2(xx)
    FLOAT, intent(in) :: xx

    dlog2 = log(xx)/log(M_TWO)
  end function dlog2

  ! -------------------------------------------------------

  integer pure function ilog2(xx)
    integer, intent(in) :: xx

    ilog2 = nint(log2(TOFLOAT(xx)))
  end function ilog2

  ! -------------------------------------------------------

  integer(i8) pure function llog2(xx)
    integer(i8), intent(in) :: xx

    llog2 = nint(log2(TOFLOAT(xx)), kind=i8)
  end function llog2

  ! -------------------------------------------------------

  logical function is_prime(n)
    integer, intent(in) :: n

    integer :: i, root

    PUSH_SUB(is_prime)

    if (n < 1) then
      message(1) = "Internal error: is_prime does not take negative numbers."
      call messages_fatal(1)
    end if
    if (n == 1) then
      is_prime = .false.
      POP_SUB(is_prime)
      return
    end if

    root = nint(sqrt(TOFLOAT(n)))
    do i = 2, root
      if (mod(n,i) == 0) then
        is_prime = .false.
        POP_SUB(is_prime)
        return
      end if
    end do

    is_prime = .true.
    POP_SUB(is_prime)
  end function is_prime

  ! ---------------------------------------------------------
  !>  Generates a rotation matrix R to rotate a vector f to t.
  !!
  !!	T. Möller and J. F. Hughes, Journal of Graphics Tools 4, 1 (1999)
  !!
  subroutine generate_rotation_matrix(R, ff, tt)
    FLOAT,   intent(out)  :: R(:,:)
    FLOAT,   intent(in)   :: ff(:)
    FLOAT,   intent(in)   :: tt(:)

    integer            :: dim, i, j
    FLOAT              :: th, uv, uu, vv, ft
    FLOAT, allocatable :: axis(:), u(:), v(:), f(:), t(:), p(:)

    PUSH_SUB(generate_rotation_matrix)

    dim = size(ff,1)

    ASSERT((dim < 3) .or. (dim > 2))
    ASSERT(size(tt,1) == dim)
    ASSERT((size(R,1) == dim) .and. (size(R,2) == dim))

    SAFE_ALLOCATE(u(1:dim))
    SAFE_ALLOCATE(v(1:dim))
    SAFE_ALLOCATE(f(1:dim))
    SAFE_ALLOCATE(t(1:dim))


    !normalize
    f = ff / norm2(ff)
    t = tt / norm2(tt)

    ft = dot_product(f,t)

    if (abs(ft) < M_ONE) then
      select case (dim)
      case (2)
        th = acos(ft)
        R(1,1) = cos(th)
        R(1,2) = -sin(th)

        R(2,1) = sin(th)
        R(2,2) = cos(th)

      case (3)
        if (.false.) then
          !Old implementation
          SAFE_ALLOCATE(axis(1:dim))
          th = acos(ft)

          u = f / dot_product(f,f)
          v = t /dot_product(t,t)

          axis = dcross_product(u,v)
          axis = axis / norm2(axis)

          R(1,1) = cos(th) + axis(1)**2 * (1 - cos(th))
          R(1,2) = axis(1)*axis(2)*(1-cos(th)) + axis(3)*sin(th)
          R(1,3) = axis(1)*axis(3)*(1-cos(th)) - axis(2)*sin(th)

          R(2,1) = axis(2)*axis(1)*(1-cos(th)) - axis(3)*sin(th)
          R(2,2) = cos(th) + axis(2)**2 * (1 - cos(th))
          R(2,3) = axis(2)*axis(3)*(1-cos(th)) + axis(1)*sin(th)

          R(3,1) = axis(3)*axis(1)*(1-cos(th)) + axis(2)*sin(th)
          R(3,2) = axis(3)*axis(2)*(1-cos(th)) - axis(1)*sin(th)
          R(3,3) = cos(th) + axis(3)**2 * (1 - cos(th))

          SAFE_DEALLOCATE_A(axis)
        end if

        if (.true.) then
          !Naive implementation
          th = acos(ft)
          u = dcross_product(f,t)
          u = u / norm2(u)

          R(1,1) = u(1)**2 + (1-u(1)**2)*cos(th)
          R(1,2) = u(1)*u(2)*(1-cos(th)) - u(3)*sin(th)
          R(1,3) = u(1)*u(3) + u(2)*sin(th)

          R(2,1) = u(1)*u(2)*(1-cos(th)) + u(3)*sin(th)
          R(2,2) = u(2)**2 + (1-u(2)**2)*cos(th)
          R(2,3) = u(2)*u(3)*(1-cos(th)) - u(1)*sin(th)

          R(3,1) = u(1)*u(3)*(1-cos(th)) - u(2)*sin(th)
          R(3,2) = u(2)*u(3)*(1-cos(th)) + u(1)*sin(th)
          R(3,3) = u(3)**2 + (1-u(3)**2)*cos(th)
        end if

        if (.false.) then
          !Fast
          SAFE_ALLOCATE(p(1:dim))

          if (abs(f(1)) <= abs(f(2)) .and. abs(f(1)) < abs(f(3))) then
            p = (/M_ONE, M_ZERO, M_ZERO/)
          else if (abs(f(2)) < abs(f(1)) .and. abs(f(2)) <= abs(f(3))) then
            p = (/M_ZERO, M_ONE, M_ZERO/)
          else if (abs(f(3)) <= abs(f(1)) .and. abs(f(3)) < abs(f(2))) then
            p = (/M_ZERO, M_ZERO, M_ONE/)
          end if

          u = p - f
          v = p - t

          uu = dot_product(u,u)
          vv = dot_product(v,v)
          uv = dot_product(u,v)

          do i = 1,3
            do j = 1,3

              R(i,j) = ddelta(i,j) - M_TWO * u(i)*u(j)/uu - M_TWO * v(i)*v(j)/vv &
                + CNST(4)*uv * v(i)*u(j) /(uu*vv)

            end do
          end do

          SAFE_DEALLOCATE_A(p)
        end if


      end select

    else

      R = M_ZERO
      do i = 1,dim
        R(i,i) = M_ONE
      end do

    end if


    SAFE_DEALLOCATE_A(u)
    SAFE_DEALLOCATE_A(v)

    POP_SUB(generate_rotation_matrix)
  end subroutine generate_rotation_matrix

  ! ---------------------------------------------------------
  !> Numerical derivative (Ridder`s algorithm).
  !!
  !! This is an alternative to "loct_numerical_derivative" (which
  !! is just an interface to the GSL numerical derivative). This version
  !! is an implementation of Ridders algorithm [C. J. F. Ridders, Adv.
  !! Eng. Software 4, 75 (1982); also described in Numerical Recipes].
  !! It is more precise, but also typically more expensive, than  the
  !! simpler 4-point algorithm implemented in the GSL library.
  subroutine numder_ridders(x, h, res, err, f)
    implicit none
    real(r8), intent(in)  :: x, h
    real(r8), intent(out) :: res, err
    interface
      subroutine f(x, fx)
        use kind_oct_m
        implicit none
        real(r8), intent(in)    :: x
        real(r8), intent(inout) :: fx
      end subroutine f
    end interface

    real(r8) :: con = CNST(1.4), big = CNST(1.0e30), safe = M_TWO
    integer :: ntab = 20, i, j
    real(r8) :: errt, fac, hh, fx1, fx2
    real(r8), allocatable :: a(:, :)

    PUSH_SUB(numder_ridders)

    if (abs(h) <= M_EPSILON) then
      message(1) = "h must be nonzero in numder_ridders"
      call messages_fatal(1)
    end if

    SAFE_ALLOCATE(a(1:ntab, 1:ntab))

    hh = h
    call f(x+hh, fx1)
    call f(x-hh, fx2)
    a(1,1) = (fx1-fx2) / (M_TWO*hh)
    err = big
    do i = 2, ntab
      hh = hh / con
      call f(x+hh, fx1)
      call f(x-hh, fx2)
      a(1,i) = (fx1-fx2) / (M_TWO*hh)
      fac = con**2
      do j = 2, i
        a(j,i) = (a(j-1,i)*fac-a(j-1,i-1)) / (fac-M_ONE)
        fac = con**2*fac
        errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
        if (errt .le. err) then
          err = errt
          res = a(j,i)
        end if
      end do
      if (abs(a(i,i)-a(i-1,i-1)) .ge. safe*err) return
    end do

    SAFE_DEALLOCATE_A(a)
    POP_SUB(numder_ridders)
  end subroutine numder_ridders

  ! ---------------------------------------------------------
  pure function dzcross_product(a, b) result(c)
    FLOAT, intent(in) :: a(:) !< (3)
    CMPLX, intent(in) :: b(:) !< (3)

    CMPLX :: c(1:3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function dzcross_product

  ! ---------------------------------------------------------
  pure function zdcross_product(a, b) result(c)
    CMPLX, intent(in) :: a(:) !< (3)
    FLOAT, intent(in) :: b(:) !< (3)

    CMPLX :: c(1:3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function zdcross_product


!*****************************************************************************80
!
!! LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
!
!  First terms:
!
!    M = 0
!
!    Lm(0,0,X) =   1
!    Lm(1,0,X) =  -X   +  1
!    Lm(2,0,X) =   X^2 -  4 X   +  2
!    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
!    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
!    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
!    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
!
!    M = 1
!
!    Lm(0,1,X) =    0
!    Lm(1,1,X) =   -1,
!    Lm(2,1,X) =    2 X - 4,
!    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
!    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
!
!    M = 2
!
!    Lm(0,2,X) =    0
!    Lm(1,2,X) =    0,
!    Lm(2,2,X) =    2,
!    Lm(3,2,X) =   -6 X + 18,
!    Lm(4,2,X) =   12 X^2 - 96 X + 144
!
!    M = 3
!
!    Lm(0,3,X) =    0
!    Lm(1,3,X) =    0,
!    Lm(2,3,X) =    0,
!    Lm(3,3,X) =   -6,
!    Lm(4,3,X) =   24 X - 96
!
!    M = 4
!
!    Lm(0,4,X) =    0
!    Lm(1,4,X) =    0
!    Lm(2,4,X) =    0
!    Lm(3,4,X) =    0
!    Lm(4,4,X) =   24
!
!  Recursion:
!
!    Lm(0,M,X)   = 1
!    Lm(1,M,X)   = (M+1-X)
!
!    if 2 <= N:
!
!      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X)
!                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
!
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal
!    to the Laguerre polynomials L(N,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
!
!    Input, real ( kind = rk ) X(MM), the evaluation points.
!
!    Output, real ( kind = rk ) CX(MM,0:N), the associated Laguerre polynomials
!    of degrees 0 through N evaluated at the evaluation points.
!
! Taken from
! https://people.sc.fsu.edu/~jburkardt/f_src/laguerre_polynomial/laguerre_polynomial.html
! and adapted to Octopus by N. Tancogne-Dejean
  subroutine generalized_laguerre_polynomial ( np, nn, mm, xx, cx )
    integer, intent(in)  :: np
    integer, intent(in)  :: nn, mm
    FLOAT,   intent(in)  :: xx(np)
    FLOAT,   intent(out) :: cx(np)

    integer ii
    FLOAT, allocatable :: cx_tmp(:,:)

    ASSERT(mm >= 0)

    SAFE_ALLOCATE(cx_tmp(1:np, 0:nn))

    if (nn < 0) then
      cx = M_ZERO
      return
    end if

    if (nn == 0) then
      cx(1:np) = M_ONE
      return
    end if

    cx_tmp(1:np,0) = M_ONE
    cx_tmp(1:np,1) = TOFLOAT(mm + 1) - xx(1:np)

    do ii = 2, nn
      cx_tmp(1:np, ii) = &
        (( TOFLOAT(mm + 2*ii - 1) - xx(1:np)) * cx_tmp(1:np,ii-1)   &
        + TOFLOAT(-mm    -ii + 1)             * cx_tmp(1:np,ii-2)) / TOFLOAT(ii)
    end do

    cx(1:np) = cx_tmp(1:np, nn)

    SAFE_DEALLOCATE_A(cx_tmp)

  end subroutine generalized_laguerre_polynomial


#include "undef.F90"
#include "complex.F90"
#include "math_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "math_inc.F90"

end module math_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
