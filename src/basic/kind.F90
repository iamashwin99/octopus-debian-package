!! Copyright (C) 2022 S. Ohlmann
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

module kind_oct_m
  use iso_fortran_env

  implicit none

  !> define kinds
  integer, parameter :: i4 = int32
  integer, parameter :: i8 = int64
  integer, parameter :: r4 = real32
  integer, parameter :: r8 = real64
end module kind_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
