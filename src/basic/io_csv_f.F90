!! Copyright (C) 2009 X. Andrade
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
#include "io_binary.h"

module io_csv_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m

  implicit none

  private

  public ::        &
    dread_csv,     &
    io_csv_get_info

contains

  subroutine dread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer(i8),         intent(in)    :: np
    real(r8),             intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

    interface
      subroutine read_csv(np, f, output_type, ierr, fname)
        use iso_c_binding
        implicit none
        integer(C_LONG),  intent(in)  :: np
        real(C_DOUBLE),   intent(in)  :: f
        integer(C_INT),   intent(in)  :: output_type
        integer(C_INT),   intent(out) :: ierr
        character(len=*), intent(in)  :: fname
      end subroutine read_csv
    end interface

    PUSH_SUB(dread_csv)

    call read_csv(np, ff(1), TYPE_DOUBLE, ierr, trim(fname))

    POP_SUB(dread_csv)
  end subroutine dread_csv

  subroutine io_csv_get_info(fname, dims, ierr)
    character(len=*),    intent(in)    :: fname
    integer(i8),         intent(inout) :: dims(:)
    integer,             intent(out)   :: ierr

    interface
      subroutine get_info_csv(dims, ierr, fname)
        use iso_c_binding
        implicit none
        integer(C_LONG),  intent(inout) :: dims(:)
        integer(C_INT),   intent(out)   :: ierr
        character(len=*), intent(in)    :: fname
      end subroutine get_info_csv
    end interface

    PUSH_SUB(io_csv_get_info)

    call get_info_csv(dims, ierr, trim(fname))

    POP_SUB(io_csv_get_info)
  end subroutine io_csv_get_info

end module io_csv_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
