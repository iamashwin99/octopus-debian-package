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

module logrid_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::         &
    logrid_t,       &
    logrid_init,    &
    logrid_end,     &
    logrid_copy,    &
    logrid_index,   &
    logrid_radius,  &
    logrid_find_parameters, &
    derivative_in_log_grid

  integer, parameter, public :: &
    LOGRID_PSF  = 1, & !< log grid used in Troullier-Martins code
    LOGRID_CPI  = 2    !< log grid used in FHI code

  type logrid_t
    ! Components are public by default
    integer  :: flavor

    FLOAT    :: a, b
    integer  :: nrval

    FLOAT, allocatable :: rofi(:)  !< r value of the point i
    FLOAT, allocatable :: r2ofi(:) !< r value of the point i
    FLOAT, allocatable :: drdi(:)  !< Jacobian, i.e., the derivative of r in terms of i
    FLOAT, allocatable :: s(:)     !< sqrt of drdi
  end type logrid_t

contains

  ! ---------------------------------------------------------
  subroutine logrid_init(grid, flavor, aa, bb, nrval)
    type(logrid_t), intent(out) :: grid
    integer,        intent(in)  :: flavor
    FLOAT,          intent(in)  :: aa, bb
    integer,        intent(in)  :: nrval

    FLOAT :: rpb, ea
    integer  :: ir

    PUSH_SUB(logrid_init)

    ASSERT(flavor == LOGRID_PSF .or. flavor == LOGRID_CPI)

    grid%flavor = flavor
    grid%a = aa
    grid%b = bb
    grid%nrval = nrval

    SAFE_ALLOCATE(grid%rofi(1:nrval))
    SAFE_ALLOCATE(grid%r2ofi(1:nrval))
    SAFE_ALLOCATE(grid%drdi(1:nrval))
    SAFE_ALLOCATE(grid%s(1:nrval))

    select case (grid%flavor)
    case (LOGRID_PSF)
      rpb = bb
      ea  = exp(aa)
      do ir = 1, nrval
        grid%drdi(ir) = aa*rpb
        rpb           = rpb*ea
        grid%rofi(ir) = bb*(exp(aa*(ir-1)) - M_ONE)
      end do

    case (LOGRID_CPI)
      grid%rofi(1) = M_ZERO
      grid%drdi(1) = M_ZERO

      rpb = log(aa)
      grid%rofi(2) = bb
      grid%drdi(2) = bb*rpb
      do ir = 3, grid%nrval
        grid%rofi(ir) = grid%rofi(ir-1)*aa
        grid%drdi(ir) = grid%rofi(ir)*rpb
      end do
    end select

    ! calculate sqrt(drdi)
    do ir = 1, grid%nrval
      grid%s(ir)     = sqrt(grid%drdi(ir))
      grid%r2ofi(ir) = grid%rofi(ir)**2
    end do

    POP_SUB(logrid_init)
  end subroutine logrid_init


  ! ---------------------------------------------------------
  subroutine logrid_end(grid)
    type(logrid_t), intent(inout) :: grid

    PUSH_SUB(logrid_end)

    SAFE_DEALLOCATE_A(grid%rofi)
    SAFE_DEALLOCATE_A(grid%r2ofi)
    SAFE_DEALLOCATE_A(grid%drdi)
    SAFE_DEALLOCATE_A(grid%s)

    POP_SUB(logrid_end)
  end subroutine logrid_end


  ! ---------------------------------------------------------
  subroutine logrid_copy(grid_in, grid_out)
    type(logrid_t), intent(in)    :: grid_in
    type(logrid_t), intent(inout) :: grid_out

    PUSH_SUB(logrid_copy)

    call logrid_end(grid_out)

    grid_out%flavor = grid_in%flavor
    grid_out%a      = grid_in%a
    grid_out%b      = grid_in%b
    grid_out%nrval  = grid_in%nrval

    SAFE_ALLOCATE(grid_out%rofi (1:grid_out%nrval))
    SAFE_ALLOCATE(grid_out%r2ofi(1:grid_out%nrval))
    SAFE_ALLOCATE(grid_out%drdi (1:grid_out%nrval))
    SAFE_ALLOCATE(grid_out%s    (1:grid_out%nrval))

    grid_out%rofi(:)  = grid_in%rofi(:)
    grid_out%r2ofi(:) = grid_in%r2ofi(:)
    grid_out%drdi(:)  = grid_in%drdi(:)
    grid_out%s(:)     = grid_in%s(:)

    POP_SUB(logrid_copy)
  end subroutine logrid_copy


  ! ---------------------------------------------------------
  integer function logrid_index(grid, rofi) result(ii)
    type(logrid_t), intent(in) :: grid
    FLOAT,          intent(in) :: rofi

    integer :: ir

    PUSH_SUB(logrid_index)

    ii = 0
    do ir = 1, grid%nrval-1

      if (rofi >= grid%rofi(ir).and.rofi < grid%rofi(ir+1)) then
        if (abs(rofi-grid%rofi(ir)) < abs(rofi-grid%rofi(ir+1))) then
          ii = ir
        else
          ii = ir + 1
        end if
        exit
      end if

    end do

    POP_SUB(logrid_index)
  end function logrid_index


  ! ---------------------------------------------------------
  subroutine derivative_in_log_grid(grid, ff, dfdr)
    type(logrid_t), intent(in)   :: grid
    FLOAT,          intent(in)   :: ff(:)
    FLOAT,          intent(out)  :: dfdr(:)

    integer :: ii

    PUSH_SUB(derivative_in_log_grid)

    dfdr(1) = (ff(2) - ff(1))/(grid%rofi(2) - grid%rofi(1))
    do ii = 2, grid%nrval-1
      dfdr(ii) = (ff(ii+1) - ff(ii-1))/(grid%rofi(ii+1) - grid%rofi(ii-1))
    end do
    dfdr(grid%nrval) = (ff(grid%nrval) - ff(grid%nrval-1))/(grid%rofi(grid%nrval) - grid%rofi(grid%nrval-1))

    POP_SUB(derivative_in_log_grid)
  end subroutine derivative_in_log_grid

  ! ----------------------------------------------------------
  FLOAT pure function logrid_radius(grid) result(radius)
    type(logrid_t), intent(in)   :: grid

    radius = grid%rofi(grid%nrval)
  end function logrid_radius


  subroutine logrid_find_parameters(namespace, zz, aa, bb, np)
    type(namespace_t),  intent(in)  :: namespace
    integer,            intent(in)  :: zz
    FLOAT,              intent(out) :: aa, bb
    integer,            intent(out) :: np

    FLOAT :: xmin, xmax, a1, a2, f1, fm

    PUSH_SUB(logrid_find_parameters)

    ! Initializes the logarithmic grid.
    ! Parameters are obtained using the default values for the first non-zero point xmin,
    ! the last point xmax, and the number of points np
    ! These values have a default value obtained from the atomic number
    ! Adapted from APE
    xmin = sqrt(TOFLOAT(zz))*CNST(1e-5)
    xmax = sqrt(TOFLOAT(zz))*CNST(30.0)
    np = floor(sqrt(TOFLOAT(zz))*CNST(200))
    ! The code wants np to be an odd number
    np = floor(np/M_TWO)*2+1

    a1 = CNST(1e-8)
    f1 = func(xmin, xmax, TOFLOAT(np), a1)
    a2 = M_ONE
    do
      aa = (a2 + a1)*M_HALF
      fm = func(xmin, xmax, TOFLOAT(np), aa)
      if (M_HALF*abs(a1 - a2) < CNST(1.0e-16)) exit
      if (fm*f1 > M_ZERO) then
        a1 = aa
        f1 = fm
      else
        a2 = aa
      end if
    end do

    bb = xmin/(exp(aa)-M_ONE)

    if (debug%info) then
      write(message(1), '(a,es13.6,a,es13.6,a,i4)') 'Debug: Log grid parameters: a = ', aa, &
        ' b = ', bb, ' np = ', np
      call messages_info(1, namespace=namespace)
    end if

    POP_SUB(logrid_find_parameters)
  contains
    FLOAT function func(r1, rn, n, a)
      FLOAT, intent(in) :: r1, rn, a, n
      if((n-M_ONE)*a < M_MAX_EXP_ARG) then ! To avoid FPE
        func = exp((n-M_ONE)*a)*r1 - M_ONE*r1 - rn*exp(a) + rn*M_ONE
      else
        func = M_HUGE*r1 - M_ONE*r1 - rn*exp(a) + rn*M_ONE
      end if
    end function func
  end subroutine logrid_find_parameters

end module logrid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
