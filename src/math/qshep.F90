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

module qshep_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::            &
    qshep_t,           &
    qshep_init,        &
    qshep_interpolate, &
    qshep_end

  type qshepr_t
    private
    integer(i8) :: npoints, nq, nw, nr, dim
    integer(i8), allocatable :: lcell(:, :, :), lnext(:)
    real(r8), allocatable :: rsq(:), a(:, :)
    real(r8) :: xmin, ymin, dx, dy, rmax, xyzmin(3), xyzdel(3)

    real(r8), allocatable :: x(:), y(:), z(:)
  end type qshepr_t

  type qshep_t
    private
    integer(i8) :: kind !< 0 for real functions (im is not used);  for complex ones.
    type(qshepr_t) :: re, im
  end type qshep_t

  interface qshep_init
    module procedure dqshep_init, zqshep_init
  end interface

  interface qshep_interpolate
    module procedure dqshep_interpolate, zqshep_interpolate
  end interface qshep_interpolate

  interface
    real(r8) function qs2val(px, py, n, x, y, f, nr, lcell, lnext, xmin, &
      ymin, dx, dy, rmax, rsq, a)
      use kind_oct_m
      implicit none
      real(r8) :: px
      real(r8) :: py
      integer(i8) :: n
      real(r8) :: x(n)
      real(r8) :: y(n)
      real(r8) :: f(n)
      integer(i8) :: nr
      integer(i8) :: lcell(nr, nr)
      integer(i8) :: lnext(n)
      real(r8) :: xmin
      real(r8) :: ymin
      real(r8) :: dx
      real(r8) :: dy
      real(r8) :: rmax
      real(r8) :: rsq(n)
      real(r8) :: a(5, n)
    end function qs2val

    subroutine qs2grd(px, py, n, x, y, f, nr, lcell, lnext, xmin, &
      ymin, dx, dy, rmax, rsq, a, q, qx, qy, ier)
      use kind_oct_m
      implicit none
      real(r8) :: px
      real(r8) :: py
      integer(i8) :: n
      real(r8) :: x(n)
      real(r8) :: y(n)
      real(r8) :: f(n)
      integer(i8) :: nr
      integer(i8) :: lcell(nr, nr)
      integer(i8) :: lnext(n)
      real(r8) :: xmin
      real(r8) :: ymin
      real(r8) :: dx
      real(r8) :: dy
      real(r8) :: rmax
      real(r8) :: rsq(n)
      real(r8) :: a(5, n)
      real(r8) :: q
      real(r8) :: qx
      real(r8) :: qy
      integer(i8) :: ier
    end subroutine qs2grd

    real(r8) function qs3val(px, py, pz, n, x, y, z, f, nr, lcell, lnext, xyzmin, &
      xyzdel, rmax, rsq, a)
      use kind_oct_m
      implicit none
      real(r8) :: px
      real(r8) :: py
      real(r8) :: pz
      integer(i8) :: n
      real(r8) :: x(n)
      real(r8) :: y(n)
      real(r8) :: z(n)
      real(r8) :: f(n)
      integer(i8) :: nr
      integer(i8) :: lcell(nr, nr, nr)
      integer(i8) :: lnext(n)
      real(r8) :: xyzmin(3)
      real(r8) :: xyzdel(3)
      real(r8) :: rmax
      real(r8) :: rsq(n)
      real(r8) :: a(9, n)
    end function qs3val

    subroutine qs3grd(px, py, pz, n, x, y, z, f, nr, lcell, lnext, xyzmin, &
      xyzdel, rmax, rsq, a, q, qx, qy, qz, ier)
      use kind_oct_m
      implicit none
      real(r8) :: px
      real(r8) :: py
      real(r8) :: pz
      integer(i8) :: n
      real(r8) :: x(n)
      real(r8) :: y(n)
      real(r8) :: z(n)
      real(r8) :: f(n)
      integer(i8) :: nr
      integer(i8) :: lcell(nr, nr, nr)
      integer(i8) :: lnext(n)
      real(r8) :: xyzmin(3)
      real(r8) :: xyzdel(3)
      real(r8) :: rmax
      real(r8) :: rsq(n)
      real(r8) :: a(9, n)
      real(r8) :: q
      real(r8) :: qx
      real(r8) :: qy
      real(r8) :: qz
      integer(i8) :: ier
    end subroutine qs3grd
  end interface

contains

  ! ------------------------------------------------------------------------------
  subroutine dqshep_init(interp, npoints, f, x, y, z)
    type(qshep_t), intent(out) :: interp
    integer(i8),    intent(in)  :: npoints
    real(r8),       intent(in)  :: f(:)
    real(r8)                    :: x(:), y(:)
    real(r8),          optional :: z(:)

    PUSH_SUB(dqshep_init)

    interp%kind = 0
    call qshepr_init(interp%re, npoints, f, x, y, z = z)

    POP_SUB(dqshep_init)
  end subroutine dqshep_init


  ! ------------------------------------------------------------------------------
  subroutine zqshep_init(interp, npoints, f, x, y, z)
    type(qshep_t), intent(out) :: interp
    integer(i8),    intent(in)  :: npoints
    complex(r8),    intent(in)  :: f(:)
    real(r8)                    :: x(:), y(:)
    real(r8),       optional    :: z(:)

    PUSH_SUB(zqshep_init)

    interp%kind = 1
    call qshepr_init(interp%re, npoints, real(f),  x, y, z = z)
    call qshepr_init(interp%im, npoints, aimag(f), x, y, z = z)

    POP_SUB(zqshep_init)
  end subroutine zqshep_init


  ! ------------------------------------------------------------------------------
  subroutine qshepr_init(interp, npoints, f, x, y, z)
    type(qshepr_t), intent(out) :: interp
    integer(i8),     intent(in) :: npoints
    real(r8),        intent(in) :: f(:)
    real(r8),        target     :: x(:), y(:)
    real(r8),  target, optional :: z(:)

    integer(i8) :: ier

    PUSH_SUB(qshepr_init)

    interp%npoints = npoints

    if (present(z)) then
      interp%dim = 3
    else
      interp%dim = 2
    end if

    select case (interp%dim)
    case (2)
      interp%nq = 13
      interp%nw = 19
      interp%nr = nint(sqrt(interp%npoints/3.0_8))
      SAFE_ALLOCATE(interp%lcell(1:interp%nr, 1:interp%nr, 1))
      SAFE_ALLOCATE(interp%a(1:5, 1:npoints))
    case (3)
      interp%nq = 17 ! This is the recommended value in qshep3d.f90
      interp%nw = 16 ! The recommended value in qshep3d.f90 is 32, but this speeds up things.
      interp%nr = nint((interp%npoints/3.0_8)**(1.0_8/3.0_8))
      SAFE_ALLOCATE(interp%lcell(1:interp%nr, 1:interp%nr, 1:interp%nr))
      SAFE_ALLOCATE(interp%a(1:9, 1:interp%npoints))
    end select

    SAFE_ALLOCATE(interp%lnext(1:npoints))
    SAFE_ALLOCATE(interp%rsq(1:npoints))

    select case (interp%dim)
    case (2)
      call qshep2(npoints, x, y, f, interp%nq, interp%nw, interp%nr, interp%lcell(:, :, 1), &
        interp%lnext, interp%xmin, interp%ymin, interp%dx, interp%dy, &
        interp%rmax, interp%rsq, interp%a, ier)
      SAFE_ALLOCATE(interp%x(1:npoints))
      SAFE_ALLOCATE(interp%y(1:npoints))
      interp%x(1:npoints) = x(1:npoints)
      interp%y(1:npoints) = y(1:npoints)
    case (3)
      call qshep3(npoints, x, y, z, f, interp%nq, interp%nw, interp%nr, interp%lcell, &
        interp%lnext, interp%xyzmin, interp%xyzdel, interp%rmax, interp%rsq, &
        interp%a, ier)
      SAFE_ALLOCATE(interp%x(1:npoints))
      SAFE_ALLOCATE(interp%y(1:npoints))
      SAFE_ALLOCATE(interp%z(1:npoints))
      interp%x(1:npoints) = x(1:npoints)
      interp%y(1:npoints) = y(1:npoints)
      interp%z(1:npoints) = z(1:npoints)
    end select

    POP_SUB(qshepr_init)
  end subroutine qshepr_init


  ! ------------------------------------------------------------------------------
  real(r8) function qshep_interpolater(interp, f, p, gf) result(v)
    type(qshepr_t),    intent(in)    :: interp
    real(r8),           intent(in)    :: f(:)
    real(r8),           intent(in)    :: p(:)
    real(r8), optional, intent(inout) :: gf(:)

    integer(i8) :: ier

    PUSH_SUB(qshep_interpolater)

    select case (interp%dim)
    case (2)
      if (present(gf)) then
        call qs2grd( p(1), p(2), interp%npoints, interp%x, interp%y, &
          f, interp%nr, interp%lcell(:, :, 1), interp%lnext, interp%xmin, &
          interp%ymin, interp%dx, interp%dy, interp%rmax, interp%rsq, interp%a, &
          v, gf(1), gf(2), ier)
      else
        v = qs2val ( p(1), p(2), interp%npoints, interp%x, interp%y, &
          f, interp%nr, interp%lcell(:, :, 1), interp%lnext, interp%xmin, &
          interp%ymin, interp%dx, interp%dy, interp%rmax, interp%rsq, interp%a)
      end if
    case (3)
      if (present(gf)) then
        call qs3grd( p(1), p(2), p(3), interp%npoints, interp%x, interp%y, interp%z, &
          f, interp%nr, interp%lcell, interp%lnext, interp%xyzmin, &
          interp%xyzdel, interp%rmax, interp%rsq, interp%a , &
          v, gf(1), gf(2), gf(3), ier)
      else
        v = qs3val (p(1), p(2), p(3), interp%npoints, interp%x, interp%y, interp%z, &
          f, interp%nr, interp%lcell, interp%lnext, interp%xyzmin, &
          interp%xyzdel, interp%rmax, interp%rsq, interp%a)
      end if
    end select

    POP_SUB(qshep_interpolater)
  end function qshep_interpolater


  ! ------------------------------------------------------------------------------
  real(r8) function dqshep_interpolate(interp, f, p, gf) result(v)
    type(qshep_t),     intent(in)    :: interp
    real(r8),           intent(in)    :: f(:)
    real(r8),           intent(in)    :: p(:)
    real(r8), optional, intent(inout) :: gf(:)

    PUSH_SUB(dqshep_interpolate)

    v = qshep_interpolater(interp%re, f, p, gf = gf)

    POP_SUB(dqshep_interpolate)
  end function dqshep_interpolate


  ! ------------------------------------------------------------------------------
  complex(r8) function zqshep_interpolate(interp, f, p, gf) result(v)
    type(qshep_t),        intent(in)    :: interp
    complex(r8),           intent(in)    :: f(:)
    real(r8),              intent(in)    :: p(:)
    complex(r8), optional, intent(inout) :: gf(:)

    integer(i8) :: i
    real(r8), allocatable :: rgf(:), igf(:)

    PUSH_SUB(zqshep_interpolate)

    if (present(gf)) then
      SAFE_ALLOCATE(rgf(1:size(gf)))
      SAFE_ALLOCATE(igf(1:size(gf)))
      v = cmplx(qshep_interpolater(interp%re, real(f), p, rgf), &
        qshep_interpolater(interp%im, aimag(f), p, igf), 8)
      do i = 1, size(gf)
        gf(i) = cmplx( rgf(i), igf(i), 8)
      end do
      SAFE_DEALLOCATE_A(rgf)
      SAFE_DEALLOCATE_A(igf)
    else
      v = cmplx(qshep_interpolater(interp%re, real(f), p), &
        qshep_interpolater(interp%im, aimag(f), p), 8)
    end if

    POP_SUB(zqshep_interpolate)
  end function zqshep_interpolate


  ! ------------------------------------------------------------------------------
  subroutine qshep_end(interp)
    type(qshep_t), intent(inout) :: interp

    PUSH_SUB(qshep_end)

    call qshepr_end(interp%re)
    if (interp%kind == 1) call qshepr_end(interp%im)

    POP_SUB(qshep_end)
  end subroutine qshep_end


  ! ------------------------------------------------------------------------------
  subroutine qshepr_end(interp)
    type(qshepr_t), intent(inout) :: interp

    PUSH_SUB(qshepr_end)

    SAFE_DEALLOCATE_A(interp%lcell)
    SAFE_DEALLOCATE_A(interp%lnext)
    SAFE_DEALLOCATE_A(interp%rsq)
    SAFE_DEALLOCATE_A(interp%a)
    SAFE_DEALLOCATE_A(interp%x)
    SAFE_DEALLOCATE_A(interp%y)
    SAFE_DEALLOCATE_A(interp%z)

    POP_SUB(qshepr_end)
  end subroutine qshepr_end

end module qshep_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
