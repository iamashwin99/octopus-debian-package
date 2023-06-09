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

module poisson_cutoff_oct_m
  use global_oct_m
  use loct_math_oct_m

  implicit none

  private
  public ::                       &
    poisson_cutoff_3D_0D,         &
    poisson_cutoff_3D_1D,         &
    poisson_cutoff_3D_1D_finite,  &
    poisson_cutoff_3D_2D,         &
    poisson_cutoff_2D_0D,         &
    poisson_cutoff_2D_1D,         &
    poisson_cutoff_1D_0D,         &
    poisson_cutoff_intcoslog


  interface poisson_cutoff_3D_1D_finite
    real(r8) function c_poisson_cutoff_3d_1d_finite(gx, gperp, xsize, rsize)
      use kind_oct_m
      implicit none
      real(r8), intent(in) :: gx, gperp, rsize, xsize
    end function c_poisson_cutoff_3d_1d_finite
  end interface poisson_cutoff_3D_1D_finite

  interface poisson_cutoff_2D_0D
    real(r8) function c_poisson_cutoff_2d_0d(x, y)
      use kind_oct_m
      implicit none
      real(r8), intent(in) :: x, y
    end function c_poisson_cutoff_2d_0d
  end interface poisson_cutoff_2D_0D

  interface poisson_cutoff_2D_1D
    real(r8) function c_poisson_cutoff_2d_1d(gy, gx, r_c)
      use kind_oct_m
      implicit none
      real(r8), intent(in) :: gy, gx, r_c
    end function c_poisson_cutoff_2d_1d
  end interface poisson_cutoff_2D_1D

  interface poisson_cutoff_1D_0D
    real(r8) function c_poisson_cutoff_1d_0d(g, a, r_c)
      use kind_oct_m
      implicit none
      real(r8), intent(in) :: g, a, r_c
    end function c_poisson_cutoff_1d_0d
  end interface poisson_cutoff_1D_0D

  interface poisson_cutoff_intcoslog
    real(r8) function intcoslog(mu, gy, gx)
      use kind_oct_m
      implicit none
      real(r8), intent(in) :: mu, gy, gx
    end function intcoslog
  end interface poisson_cutoff_intcoslog


contains


  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_3D_0D(x, r) result(cutoff)
    FLOAT, intent(in) ::  x, r

    !no PUSH_SUB, called too often

    cutoff = M_ONE - cos(x*r)

  end function poisson_cutoff_3D_0D
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_3D_1D(x, p, rmax) result(cutoff)
    FLOAT, intent(in) ::  x, p, rmax

    integer :: j
    FLOAT :: dr, r, sum

    integer :: nr = CNST(1000)

    !no PUSH_SUB, called too often

    if (abs(x) < M_EPSILON) then
      ! Simpson rule for the G_x = 0 contribution -log(r)
      dr = rmax/TOFLOAT(nr)
      sum = M_ZERO
      do j = 1, nr - 1, 2
        r = j*dr
        sum = sum - M_FOUR*r*loct_bessel_j0(p*r)*log(r)
        r = (j+1)*dr
        sum = sum - M_TWO*r*loct_bessel_j0(p*r)*log(r)
      end do
      sum = sum - rmax*loct_bessel_j0(p*rmax)*log(rmax)
      cutoff = (p**2)*M_THIRD*sum*dr
    else
      cutoff = M_ONE + p*rmax*loct_bessel_j1(p*rmax)*loct_bessel_k0(x*rmax) &
        - x*rmax*loct_bessel_j0(p*rmax)*loct_bessel_k1(x*rmax)
    end if

  end function poisson_cutoff_3D_1D

  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_3D_2D(p, z, r) result(cutoff)
    FLOAT, intent(in) ::  p, z, r

    !no PUSH_SUB, called too often

    if (abs(p) < M_EPSILON) then
      cutoff = M_ONE - cos(z*r) - z*r*sin(z*r)
    else
      cutoff = M_ONE + exp(-p*r)*(z*sin(z*r)/p-cos(z*r))
    end if

  end function poisson_cutoff_3D_2D
  ! ---------------------------------------------------------

end module poisson_cutoff_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
