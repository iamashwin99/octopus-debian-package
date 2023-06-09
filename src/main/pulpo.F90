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

module pulpo_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m

  implicit none

  private
  public :: pulpo_print

contains

  subroutine pulpo_print()

    character(len=MAX_PATH_LEN) :: filename

    PUSH_SUB(pulpo_print)

    ! some white space
    message(1) = ''
    message(2) = ''
    call messages_info(2)

    call loct_printrecipe(trim(conf%share), filename)
    call io_dump_file(stdout, filename)
    call messages_info(2)
    call io_dump_file(stdout, trim(conf%share)//"/recipes/disclaimer.txt")
    call messages_info(2)

    POP_SUB(pulpo_print)

  end subroutine pulpo_print

end module pulpo_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
