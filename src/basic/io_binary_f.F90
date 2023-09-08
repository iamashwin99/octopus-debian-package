!! Copyright (C) 2009 X. Andrade
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

#include "global.h"
#include "io_binary.h"

module io_binary_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use mpi_oct_m
  use string_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                     &
    io_binary_is_little_endian, &
    io_binary_write,            &
    io_binary_write_parallel,   &
    io_binary_read,             &
    io_binary_read_parallel,    &
    io_binary_get_info,         &
    dwrite_header,              &
    zwrite_header,              &
    iwrite_header,              &
    lwrite_header

  interface io_binary_write
    module procedure dwrite_binary, zwrite_binary, iwrite_binary, lwrite_binary
    module procedure zwrite_binary2, dwrite_binary2, iwrite_binary2, lwrite_binary2
    module procedure zwrite_binary3, dwrite_binary3, iwrite_binary3, lwrite_binary3
    module procedure zwrite_binary4, dwrite_binary4, iwrite_binary4, lwrite_binary4
    module procedure zwrite_binary5, dwrite_binary5, iwrite_binary5, lwrite_binary5
  end interface io_binary_write

  interface io_binary_write_parallel
    module procedure dwrite_parallel, zwrite_parallel, iwrite_parallel, lwrite_parallel
  end interface io_binary_write_parallel

  interface io_binary_read
    module procedure dread_binary, zread_binary, iread_binary, lread_binary
    module procedure zread_binary2, dread_binary2, iread_binary2, lread_binary2
    module procedure zread_binary3, dread_binary3, iread_binary3, lread_binary3
    module procedure zread_binary4, dread_binary4, iread_binary4, lread_binary4
    module procedure zread_binary5, dread_binary5, iread_binary5, lread_binary5
  end interface io_binary_read

  interface io_binary_read_parallel
    module procedure dread_parallel, zread_parallel, iread_parallel, lread_parallel
  end interface io_binary_read_parallel

  interface
    subroutine get_info_binary(np, type, file_size, ierr, iio, fname) bind(c)
      use iso_c_binding
      integer(c_int64_t),     intent(out)   :: np        !< Number of points of the mesh, written in the header
      integer(c_int),         intent(out)   :: type      !< Type of number
      integer(c_int64_t),     intent(out)   :: file_size !< The actual size of the file
      integer(c_int),         intent(out)   :: ierr
      integer(c_int),         intent(inout) :: iio
      character(kind=c_char), intent(in)    :: fname(*)
    end subroutine get_info_binary

    subroutine write_header(np, type, ierr, iio, fname) bind(c, name="io_write_header")
      use iso_c_binding
      integer(c_int64_t),     intent(in)    :: np
      integer(c_int),         intent(in)    :: type
      integer(c_int),         intent(out)   :: ierr
      integer(c_int),         intent(inout) :: iio
      character(kind=c_char), intent(in)    :: fname(*)
    end subroutine write_header

    subroutine write_binary(np, ff, type, ierr, iio, nhd, flpe, fname) bind(c, name="write_binary")
      use iso_c_binding
      integer(c_int64_t),     intent(in)    :: np
      type(c_ptr),            value         :: ff
      integer(c_int),         intent(in)    :: type
      integer(c_int),         intent(out)   :: ierr
      integer(c_int),         intent(inout) :: iio
      integer(c_int),         intent(in)    :: nhd
      integer(c_int),         intent(in)    :: flpe
      character(kind=c_char), intent(in)    :: fname(*)
    end subroutine write_binary

    subroutine read_binary(np, offset, ff, output_type, ierr, iio, fname) bind(c, name="read_binary")
      use iso_c_binding
      integer(c_int64_t),     intent(in)    :: np
      integer(c_int64_t),     intent(in)    :: offset
      type(c_ptr),            value         :: ff
      integer(c_int),         intent(in)    :: output_type
      integer(c_int),         intent(out)   :: ierr
      integer(c_int),         intent(inout) :: iio
      character(kind=c_char), intent(in)    :: fname(*)
    end subroutine read_binary

  end interface

contains

  !> check endianness
  !> Logical output: true is the running architecture uses little endian ordering, false otherwise.
  logical pure function io_binary_is_little_endian() result(is_little)
    implicit none
    integer, parameter:: I4P  = selected_int_kind(9)
    integer, parameter:: I1P  = selected_int_kind(2)
    integer(I1P) :: int1(1:4) !< One byte integer array for casting 4 bytes integer.

    int1 = transfer(1_I4P, int1)
    is_little = (int1(1) == 1_I1P)

  end function io_binary_is_little_endian

  ! ------------------------------------------------------
  subroutine io_binary_parallel_start(fname, file_handle, comm, xlocal, np, sizeof_ff, is_write, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(out)   :: file_handle
    integer,             intent(in)    :: comm
    integer(i8),         intent(in)    :: xlocal
    integer,             intent(in)    :: np
    integer,             intent(in)    :: sizeof_ff
    logical,             intent(in)    :: is_write !< if false, is read.
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI
    integer(MPI_OFFSET_KIND) :: offset
    integer :: amode
#endif

    PUSH_SUB(io_binary_parallel_start)

    ASSERT(np > 0)

#ifdef HAVE_MPI
    offset = (xlocal-1)*sizeof_ff+64

    if (is_write) then
      amode = IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND)
    else
      amode = MPI_MODE_RDONLY
    end if
    call MPI_File_open(comm, fname, amode, MPI_INFO_NULL, file_handle, mpi_err)
    call io_incr_open_count()

    if (mpi_err == 0) then
      call MPI_File_set_atomicity(file_handle, .true., mpi_err)
      call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    end if
    ierr = mpi_err
#endif

    POP_SUB(io_binary_parallel_start)
  end subroutine io_binary_parallel_start

  ! ------------------------------------------------------

  subroutine io_binary_parallel_end(file_handle)
    integer, intent(inout) :: file_handle

    PUSH_SUB(io_binary_parallel_end)

#ifdef HAVE_MPI
    call MPI_File_close(file_handle, mpi_err)
    call io_incr_close_count()
#endif

    POP_SUB(io_binary_parallel_end)
  end subroutine io_binary_parallel_end


  ! ------------------------------------------------------

  subroutine try_dread_binary(fname, np, ff, ierr, offset)
    character(len=*),      intent(in)  :: fname
    integer(i8),           intent(in)  :: np
    complex(r8),           intent(out) :: ff(:)
    integer,               intent(out) :: ierr
    integer(i8), optional, intent(in)  :: offset

    integer(i8) :: read_np, file_size
    integer :: number_type, iio
    real(r8), allocatable :: read_ff(:)
    character(kind=c_char), dimension(c_str_len(fname)) :: cname


    PUSH_SUB(try_dread_binary)

    iio = 0
    cname = string_f_to_c(fname)
    call get_info_binary(read_np, number_type, file_size, ierr, iio, cname)
    call io_incr_counters(iio)

    ! if the type of the file is real, then read real numbers and convert to complex
    if (number_type /= TYPE_DOUBLE_COMPLEX) then
      if (debug%info) then
        write(message(1),'(a,i2,a,i2)') "Debug: Found type = ", number_type, " instead of ", TYPE_DOUBLE_COMPLEX
        call messages_info(1)
      end if

      SAFE_ALLOCATE(read_ff(1:np))
      call dread_binary(fname, np, read_ff, ierr, offset)
      ff = read_ff
      SAFE_DEALLOCATE_A(read_ff)
    else
      ierr = -1
    end if
    ! ierr will be 0 if dread_binary succeeded

    POP_SUB(try_dread_binary)
  end subroutine try_dread_binary

  !------------------------------------------------------

  subroutine try_dread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer(i8),         intent(in)    :: xlocal
    integer,             intent(in)    :: np
    complex(r8),         intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

    integer(i8) :: read_np, file_size
    integer :: number_type, iio
    real(r8), allocatable :: read_ff(:)

    PUSH_SUB(try_dread_parallel)

    iio = 0
    call get_info_binary(read_np, number_type, file_size, ierr, iio, string_f_to_c(fname))
    call io_incr_counters(iio)
    ! if the type of the file is real, then read real numbers and convert to complex
    if (number_type /= TYPE_DOUBLE_COMPLEX) then
      if (debug%info) then
        write(message(1),'(a,i2,a,i2)') "Debug: Found type = ", number_type, " instead of ", TYPE_DOUBLE_COMPLEX
        call messages_info(1)
      end if
      SAFE_ALLOCATE(read_ff(1:np))
      call dread_parallel(fname, comm, xlocal, np, read_ff, ierr)
      ff = read_ff
      SAFE_DEALLOCATE_A(read_ff)
    else
      ierr = -1
    end if
    ! ierr will be 0 if dread_parallel succeeded

    POP_SUB(try_dread_parallel)
  end subroutine try_dread_parallel

  !------------------------------------------------------

  subroutine io_binary_get_info(fname, np, file_size, ierr)
    character(len=*),    intent(in)    :: fname
    integer(i8),         intent(out)   :: np
    integer(i8),         intent(out)   :: file_size
    integer,             intent(out)   :: ierr

    integer :: type, iio

    PUSH_SUB(io_binary_get_info)

    iio = 0
    call get_info_binary(np, type, file_size, ierr, iio, string_f_to_c(fname))
    call io_incr_counters(iio)

    POP_SUB(io_binary_get_info)
  end subroutine io_binary_get_info

  ! ------------------------------------------------------
  integer pure function logical_to_integer(flag) result(iflag)
    logical, intent(in) :: flag
    iflag = 0
    if (flag) iflag = 1
  end function logical_to_integer

#include "complex.F90"
#include "io_binary_f_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "io_binary_f_inc.F90"

#include "undef.F90"

#include "integer.F90"
#include "io_binary_f_inc.F90"

#include "undef.F90"

#include "integer8.F90"
#include "io_binary_f_inc.F90"

end module io_binary_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
