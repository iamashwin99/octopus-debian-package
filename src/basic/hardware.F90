!! Copyright (C) 2008 X. Andrade
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
module hardware_oct_m

  implicit none

  private
  public ::             &
    hardware,           &
    cache_t,            &
    hardware_t,         &
    hardware_init,      &
    hardware_end

  type cache_t
    ! Components are public by default
    integer :: size
    integer :: line_size
  end type cache_t

  type hardware_t
    ! Components are public by default
    type(cache_t) :: l1
    type(cache_t) :: l2
    integer :: dblock_size
    integer :: zblock_size
  end type hardware_t

  type(hardware_t) :: hardware

contains

  subroutine hardware_init

    !for the moment we will use fixed values

    hardware%l1%size = 32*1024

    hardware%l1%line_size = 64

    hardware%l2%size = 4096*1024

    hardware%l2%line_size = 64

    ! set the block_size so each block fits in the l1 cache
    ! the block_size should be a multiple of the cache line (minus 2 lines to avoid powers of 2)

    hardware%dblock_size = hardware%l1%size / (4*8)
    hardware%dblock_size = hardware%dblock_size - mod(hardware%dblock_size, hardware%l1%line_size) - 2*hardware%l1%line_size

    hardware%zblock_size = hardware%l1%size / (4*16)
    hardware%zblock_size = hardware%zblock_size - mod(hardware%zblock_size, hardware%l1%line_size) - 2*hardware%l1%line_size

  end subroutine hardware_init

  subroutine hardware_end

  end subroutine hardware_end

end module hardware_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
