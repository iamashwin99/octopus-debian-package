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

module ps_cpi_file_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use ps_in_grid_oct_m

  implicit none

  private

  public ::                &
    ps_cpi_file_t,          &
    ps_cpi_file_read,       &
    ps_cpi_file_end

  !> First, the contents of the file.
  type ps_cpi_file_t
    ! Components are public by default
    FLOAT              :: zval          !< valence charge
    integer            :: no_l_channels !< number of pseudo components (lmax+1)

    integer            :: nr            !< number of mesh points
    FLOAT              :: a             !< mesh multiplicative increment

    FLOAT, allocatable :: rofi(:)       !< radial mesh
    FLOAT, allocatable :: vps(:,:)      !< pseudopotential
    FLOAT, allocatable :: rphi(:,:)     !< r times the pseudowavefunctions

    logical            :: core_corrections
    FLOAT, allocatable :: chcore(:)     !< r times the core charge
    FLOAT, allocatable, private :: d1chcore(:)   !< first  derivative of chcore
    FLOAT, allocatable, private :: d2chcore(:)   !< second derivative of chcore
  end type ps_cpi_file_t

contains

  ! ---------------------------------------------------------
  subroutine ps_cpi_file_read(unit, psf)
    integer,             intent(in)    :: unit
    type(ps_cpi_file_t), intent(inout) :: psf

    integer  :: i, l, ios, idummy
    FLOAT    :: a, b, c, d

    PUSH_SUB(ps_cpi_file_read)

    read(unit, *) psf%zval, psf%no_l_channels
    ! skip 10 lines
    do i = 1, 10
      read(unit, *)
    end do

    read(unit, *) psf%nr, psf%a

    ! add extra point for the zero
    psf%nr = psf%nr + 1

    SAFE_ALLOCATE(psf%rofi   (1:psf%nr))
    SAFE_ALLOCATE(psf%vps    (1:psf%nr, 1:psf%no_l_channels))
    SAFE_ALLOCATE(psf%rphi   (1:psf%nr, 1:psf%no_l_channels))

    do l = 1, psf%no_l_channels
      if (l /= 1) read(unit, *)

      do i = 2, psf%nr
        read(unit, *) idummy, psf%rofi(i), psf%rphi(i, l), psf%vps(i, l)
      end do
    end do

    ! read core charge (if present)
    read(unit, *, iostat=ios) a, b, c, d
    if (ios == 0) then
      psf%core_corrections = .true.

      SAFE_ALLOCATE(psf%chcore  (1:psf%nr))
      SAFE_ALLOCATE(psf%d1chcore(1:psf%nr))
      SAFE_ALLOCATE(psf%d2chcore(1:psf%nr))

      psf%  chcore(2) = b
      psf%d1chcore(2) = c
      psf%d2chcore(2) = d

      do i = 3, psf%nr
        read(unit, *) a, psf%chcore(i), psf%d1chcore(i), psf%d2chcore(i)
      end do
    else
      psf%core_corrections = .false.
    end if

    ! add extra point at zero
    psf%rofi(1) = M_ZERO
    do l = 1, psf%no_l_channels
      psf%vps(1,  l) = first_point_extrapolate(psf%rofi, psf%vps(:, l))

      psf%rphi(1, l) = M_ZERO
    end do

    if (psf%core_corrections) then
      ! At this point, we use the normalization of the siesta format, where
      ! psf%chcore(:) = 4*pi*\tilde{rho} r**2. As in the Fritz-Haber file we
      ! have written 4*pi*\tilde{rho}, we multiply by r**2
      psf%chcore(:) = psf%chcore(:) * psf%rofi(:)**2

      psf%chcore(1) = first_point_extrapolate(psf%rofi, psf%chcore)
      psf%d1chcore(1) = first_point_extrapolate(psf%rofi, psf%d1chcore)
      psf%d2chcore(1) = first_point_extrapolate(psf%rofi, psf%d2chcore)
    end if

    ! WARNING: This should go away
    psf%vps(:,:) = psf%vps(:,:)*M_TWO ! convert to Rydbergs

    POP_SUB(ps_cpi_file_read)
  end subroutine ps_cpi_file_read


  ! ---------------------------------------------------------
  subroutine ps_cpi_file_end(psf)
    type(ps_cpi_file_t), intent(inout) :: psf

    PUSH_SUB(ps_cpi_file_end)

    SAFE_DEALLOCATE_A(psf%rofi)
    SAFE_DEALLOCATE_A(psf%vps)
    SAFE_DEALLOCATE_A(psf%rphi)

    if (psf%core_corrections) then
      SAFE_DEALLOCATE_A(psf%chcore)
      SAFE_DEALLOCATE_A(psf%d1chcore)
      SAFE_DEALLOCATE_A(psf%d2chcore)
    end if

    POP_SUB(ps_cpi_file_end)
  end subroutine ps_cpi_file_end

end module ps_cpi_file_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
