 !
 ! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
 ! This file is distributed under the terms of the
 ! GNU General Public License. See the file "License"
 ! in the root directory of the present distribution,
 ! or http://www.gnu.org/copyleft/gpl.txt .
 !
 ! routine adapted from the original Pwscf interface PP/pw2wannier90.f90
 !
#include "global.h"

module ylm_wannier_oct_m
  use global_oct_m
  use loct_math_oct_m


  implicit none

  private

  public :: ylm_wannier

contains
  subroutine ylm_wannier(ylm, l, mr, rr, xx, nr)
    ! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
    ! of the spherical harmonic identified  by indices (l,mr)
    ! in table 3.1 of the wannierf90 specification.
    !
    ! No reference to the particular ylm ordering internal to octopus
    ! is assumed.
    !
    ! If ordering in wannier90 code is changed or extended this should be the
    ! only place to be modified accordingly
    !

    ! I/O variables
    !
    FLOAT,   intent(inout) :: ylm(nr)
    integer, intent(in)    :: l, mr, nr
    FLOAT,   intent(inout) :: rr(nr), xx(3, nr)
    !
    ! local variables
    !
    !     FLOAT, EXTERNAL ::  s, pz,px,py, dz2, dxz, dyz, dx2my2, dxy
    !     FLOAT, EXTERNAL :: fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
    FLOAT :: cost, phi
    integer :: ir
    FLOAT :: bs2, bs3, bs6, bs12
    bs2  = M_ONE / sqrt(M_TWO)
    bs3  = M_ONE / sqrt(M_THREE)
    bs6  = M_ONE / sqrt(CNST(6.0))
    bs12 = M_ONE / sqrt(CNST(12.0))

    if (l > 3 .or. l < -5) stop 'ylm_wannier l out of range '
    if (l >= 0) then
      if (mr < 1 .or. mr > 2*l + 1) stop 'ylm_wannier mr out of range'
    else
      if (mr < 1 .or. mr > abs(l) + 1) stop 'ylm_wannier mr out of range'
    end if

    do ir = 1, nr
      ! if r=0, direction is undefined => make ylm=0 except for l=0
      if (rr(ir) < CNST(1e-12)) then
        if (l == 0) then
          ylm(ir) = s()
        else
          ylm(ir) = M_ZERO
        end if
        cycle
      end if

      cost =  xx(3, ir) / rr(ir)
      !
      !  beware the arc tan, it is defined modulo M_PI
      !
      if (xx(1, ir) > CNST(1e-12)) then
        phi = atan(xx(2, ir) / xx(1, ir))
      else if (xx(1,ir) < -CNST(1e-12)) then
        phi = atan(xx(2, ir) / xx(1, ir)) + M_PI
      else
        phi = sign(M_PI / M_TWO, xx(2, ir))
      end if

      select case (l)
      case (0)   ! s orbital
        ylm(ir) = s()
      case (1)   ! p orbitals
        if (mr == 1) ylm(ir) = pz(cost)
        if (mr == 2) ylm(ir) = px(cost, phi)
        if (mr == 3) ylm(ir) = py(cost, phi)
      case (2)   ! d orbitals
        if (mr == 1) ylm(ir) = dz2(cost)
        if (mr == 2) ylm(ir) = dxz(cost, phi)
        if (mr == 3) ylm(ir) = dyz(cost, phi)
        if (mr == 4) ylm(ir) = dx2my2(cost, phi)
        if (mr == 5) ylm(ir) = dxy(cost, phi)
      case (3)   ! f orbitals
        if (mr == 1) ylm(ir) = fz3(cost)
        if (mr == 2) ylm(ir) = fxz2(cost, phi)
        if (mr == 3) ylm(ir) = fyz2(cost, phi)
        if (mr == 4) ylm(ir) = fzx2my2(cost, phi)
        if (mr == 5) ylm(ir) = fxyz(cost, phi)
        if (mr == 6) ylm(ir) = fxx2m3y2(cost, phi)
        if (mr == 7) ylm(ir) = fy3x2my2(cost, phi)
      case (-1)  !  sp hybrids
        if (mr == 1) ylm(ir) = bs2 * (s() + px(cost, phi))
        if (mr == 2) ylm(ir) = bs2 * (s() - px(cost, phi))
      case (-2)  !  sp2 hybrids
        if (mr == 1) ylm(ir) = bs3 * s() - bs6 * px(cost, phi) + bs2 * py(cost, phi)
        if (mr == 2) ylm(ir) = bs3 * s() - bs6 * px(cost, phi) - bs2 * py(cost, phi)
        if (mr == 3) ylm(ir) = bs3 * s() + M_TWO * bs6 * px(cost, phi)
      case (-3)  !  sp3 hybrids
        if (mr == 1) ylm(ir) = M_HALF*(s() + px(cost, phi) + py(cost, phi) + pz(cost))
        if (mr == 2) ylm(ir) = M_HALF*(s() + px(cost, phi) - py(cost, phi) - pz(cost))
        if (mr == 3) ylm(ir) = M_HALF*(s() - px(cost, phi) + py(cost, phi) - pz(cost))
        if (mr == 4) ylm(ir) = M_HALF*(s() - px(cost, phi) - py(cost, phi) + pz(cost))
      case (-4)  !  sp3d hybrids
        if (mr == 1) ylm(ir) = bs3 * s() - bs6 * px(cost, phi) + bs2 * py(cost, phi)
        if (mr == 2) ylm(ir) = bs3 * s() - bs6 * px(cost, phi) - bs2 * py(cost, phi)
        if (mr == 3) ylm(ir) = bs3 * s() + M_TWO * bs6 * px(cost, phi)
        if (mr == 4) ylm(ir) = bs2 * pz(cost) + bs2 * dz2(cost)
        if (mr == 5) ylm(ir) =-bs2 * pz(cost) + bs2 * dz2(cost)
      case (-5)  ! sp3d2 hybrids
        if (mr == 1) ylm(ir) = bs6 * s() - bs2 * px(cost, phi) -bs12 * dz2(cost) + M_HALF * dx2my2(cost, phi)
        if (mr == 2) ylm(ir) = bs6 * s() + bs2 * px(cost, phi) -bs12 * dz2(cost) + M_HALF * dx2my2(cost, phi)
        if (mr == 3) ylm(ir) = bs6 * s() - bs2 * py(cost, phi) -bs12 * dz2(cost) - M_HALF * dx2my2(cost, phi)
        if (mr == 4) ylm(ir) = bs6 * s() + bs2 * py(cost, phi) -bs12 * dz2(cost) - M_HALF * dx2my2(cost, phi)
        if (mr == 5) ylm(ir) = bs6 * s() - bs2 * pz(cost) +bs3 * dz2(cost)
        if (mr == 6) ylm(ir) = bs6 * s() + bs2 * pz(cost) +bs3 * dz2(cost)
      end select

    end do
  end subroutine ylm_wannier

  !======== l = 0 =====================================================================
  function s()
    FLOAT :: s
    s = M_ONE / sqrt((M_FOUR * M_PI))
  end function s

  !======== l = 1 =====================================================================
  function pz(cost)
    FLOAT ::pz, cost
    pz =  sqrt(M_THREE / (M_FOUR * M_PI)) * cost
  end function pz

  function px(cost, phi)
    FLOAT ::px, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    px =  sqrt(M_THREE / (M_FOUR * M_PI)) * sint * cos(phi)
  end function px

  function py(cost, phi)
    FLOAT ::py, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    py =  sqrt(M_THREE / (M_FOUR * M_PI)) * sint * sin(phi)
  end function py

  !======== l = 2 =====================================================================
  function dz2(cost)
    FLOAT ::dz2, cost
    dz2 =  sqrt(CNST(1.25) / (M_FOUR * M_PI)) * (M_THREE* cost * cost - M_ONE)
  end function dz2

  function dxz(cost, phi)
    FLOAT ::dxz, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    dxz =  sqrt(CNST(15.0) / (M_FOUR * M_PI)) * sint * cost * cos(phi)
  end function dxz

  function dyz(cost, phi)
    FLOAT ::dyz, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    dyz =  sqrt(CNST(15.0) / (M_FOUR * M_PI)) * sint * cost * sin(phi)
  end function dyz

  function dx2my2(cost, phi)
    FLOAT ::dx2my2, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    dx2my2 =  sqrt(CNST(3.750) / (M_FOUR * M_PI)) * sint * sint * cos(M_TWO * phi)
  end function dx2my2

  function dxy(cost, phi)
    FLOAT ::dxy, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    dxy =  sqrt(CNST(3.750) / (M_FOUR * M_PI)) * sint * sint * sin(M_TWO * phi)
  end function dxy

  !======== l = 3 =====================================================================
  function fz3(cost)
    FLOAT ::fz3, cost
    fz3 =  CNST(0.25) * sqrt(CNST(7.0) / M_PI) * (CNST(5.0) * cost * cost - M_THREE) * cost
  end function fz3

  function fxz2(cost, phi)
    FLOAT ::fxz2, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    fxz2 =  CNST(0.25) * sqrt(CNST(10.5) / M_PI) * (CNST(5.0) * cost * cost - M_ONE) * sint * cos(phi)
  end function fxz2

  function fyz2(cost, phi)
    FLOAT ::fyz2, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    fyz2 =  CNST(0.25) * sqrt(CNST(10.5) / M_PI) * (CNST(5.0) * cost * cost - M_ONE) * sint * sin(phi)
  end function fyz2

  function fzx2my2(cost, phi)
    FLOAT ::fzx2my2, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    fzx2my2 =  CNST(0.25) * sqrt(CNST(105.0) / M_PI) * sint * sint * cost * cos(M_TWO * phi)
  end function fzx2my2

  function fxyz(cost, phi)
    FLOAT ::fxyz, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    fxyz =  CNST(0.25)*sqrt(CNST(105.0) / M_PI) * sint * sint * cost * sin(M_TWO * phi)
  end function fxyz

  function fxx2m3y2(cost, phi)
    FLOAT ::fxx2m3y2, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    fxx2m3y2 =  CNST(0.25) * sqrt(CNST(17.5) / M_PI) * sint * sint * sint * cos(M_THREE * phi)
  end function fxx2m3y2

  function fy3x2my2(cost, phi)
    FLOAT ::fy3x2my2, cost, phi, sint
    sint = sqrt(abs(M_ONE - cost * cost))
    fy3x2my2 =  CNST(0.25) * sqrt(CNST(17.5) / M_PI) * sint * sint * sint * sin(M_THREE * phi)
  end function fy3x2my2

end module ylm_wannier_oct_m
