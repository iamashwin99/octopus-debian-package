!! Copyright (C) 2010 X. Andrade
!! Copyright (C) 2021 S. Ohlmann
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

module checksum_interface_oct_m

  public ::                 &
    checksum_calculate

  interface
    subroutine checksum_calculate(algorithm, narray, array, checksum)
      use kind_oct_m
      implicit none
      integer,     intent(in)  :: algorithm
      integer(i8), intent(in)  :: narray
      integer(i8), intent(in)  :: array
      integer(i8), intent(out) :: checksum
    end subroutine checksum_calculate
  end interface

end module checksum_interface_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
