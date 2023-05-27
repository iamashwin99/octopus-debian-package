!! Copyright (C) 2005-2006 Heiko Appel, Florian Lorenzen
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

! Routines to support MPI debugging.

#include "global.h"

module mpi_debug_oct_m
#if defined(MPI_MOD)
  use mpi
#endif
  use kind_oct_m
  use loct_oct_m

  implicit none

  ! some machines do not have a mpi module, but a mpi input file
#if defined(MPI_H)
  include "mpif.h"
#endif

  private

  public ::               &
    mpi_debug_statistics, &
    mpi_debug_init,       &
    mpi_debug_in,         &
    mpi_debug_out

  integer, parameter :: C_NUM_MPI_ROUTINES = 12

  integer, public, parameter ::  &
    C_MPI_BARRIER    = 1,        &
    C_MPI_SCATTERV   = 2,        &
    C_MPI_GATHERV    = 3,        &
    C_MPI_GATHER     = 4,        &
    C_MPI_ALLTOALLV  = 5,        &
    C_MPI_ALLGATHERV = 6,        &
    C_MPI_BCAST      = 7,        &
    C_MPI_ALLREDUCE  = 8,        &
    C_MPI_ALLTOALL   = 9,        &
    C_MPI_ALLGATHER  = 10,        &
    C_MPI_FILE_READ  = 11,       &
    C_MPI_FILE_WRITE = 12

  character(len=14), dimension(C_NUM_MPI_ROUTINES), public :: mpi_rlabel = &
    (/                           &
    'MPI_BARRIER   ',            &
    'MPI_SCATTERV  ',            &
    'MPI_GATHERV   ',            &
    'MPI_GATHER    ',            &
    'MPI_ALLTOALLV ',            &
    'MPI_ALLGATHERV',            &
    'MPI_BCAST     ',            &
    'MPI_ALLREDUCE ',            &
    'MPI_ALLTOALL  ',            &
    'MPI_ALLGATHER ',            &
    'MPI_FILE_READ ',            &
    'MPI_FILE_WRITE'             &
    /)

  integer, public :: call_counter(C_NUM_MPI_ROUTINES) = 0
  real(r8), public :: sec_accum(C_NUM_MPI_ROUTINES)   = 0_r8

  real(r8) :: sec_in

  logical :: debug_info
  integer :: mpi_rank

  !> max_lun is currently 99, i.e. we can hardwire unit_offset above 1000
  integer, parameter :: unit_offset = 1000

contains

  ! ---------------------------------------------------------
  subroutine mpi_debug_init(rank, info)
    integer, intent(in) :: rank
    logical, intent(in) :: info

    mpi_rank = rank
    debug_info = info
  end subroutine mpi_debug_init

  ! ---------------------------------------------------------
  subroutine mpi_debug_open_trace(iunit)
    integer, intent(out) :: iunit

    character(len=6) :: filenum

    iunit = mpi_rank + unit_offset
    write(filenum, '(i6.6)') iunit - unit_offset
    call loct_mkdir('debug')
    open(iunit, file = 'debug/debug_trace.node.'//filenum, &
      action='write', status='unknown', position='append')

  end subroutine mpi_debug_open_trace

  ! ---------------------------------------------------------
  subroutine mpi_debug_statistics()
    integer :: j, iunit
    real(r8) :: usec_call(C_NUM_MPI_ROUTINES)

    if (.not. debug_info) return

    call mpi_debug_open_trace(iunit)

    write(iunit,*)
    write(iunit,'(A)') '--------------------------------------------------------------------'
    write(iunit,*)
    write(iunit, '(23x,a,4x,a,8x,a)') 'total time', 'calls', 'usec/call'
    do j = 1, C_NUM_MPI_ROUTINES
      if (call_counter(j) <= 0) then
        usec_call(j) = 0
      else
        usec_call(j) = (sec_accum(j)*1000000)/call_counter(j)
      end if

      write(iunit,'(a,f13.6,6x,i4,6x,f13.0)') &
        mpi_rlabel(j)//' : ', sec_accum(j),          &
        call_counter(j), usec_call(j)
    end do
    write(iunit,*)
    write(iunit,'(A)') '--------------------------------------------------------------------'

    close(iunit)

  end subroutine mpi_debug_statistics


  ! ---------------------------------------------------------
  subroutine mpi_debug_in(comm, index)
    integer, intent(in) :: comm, index

    integer :: iunit

    if (.not. debug_info) return

    call mpi_debug_open_trace(iunit)

    call_counter(index) = call_counter(index) + 1
#if defined(HAVE_MPI)
    sec_in              = MPI_Wtime()
#endif
    write(iunit,'(a,f18.6,a,z8.8,a,i6.6,a,f13.6)') '* MPI_I ', &
      sec_in, ' '//mpi_rlabel(index)//' : 0x', comm, ' | ',  &
      call_counter(index), ' - ', sec_accum(index)

    close(iunit)

  end subroutine mpi_debug_in


  ! ---------------------------------------------------------
  subroutine mpi_debug_out(comm, index)
    integer, intent(in) :: comm, index

    integer :: iunit
    real(r8) :: sec_out, sec_diff

    if (.not. debug_info) return

    call mpi_debug_open_trace(iunit)

#if defined(HAVE_MPI)
    sec_out = MPI_Wtime()
#endif
    call mpi_time_accum(index, sec_out, sec_diff)
    write(iunit,'(a,f18.6,a,z8.8,a,i6.6,a,f13.6,a,f13.6)')         &
      '* MPI_O ', sec_out, ' '//mpi_rlabel(index)//' : 0x', comm, ' | ', &
      call_counter(index), ' - ', sec_accum(index), ' - ', sec_diff

    close(iunit)

  end subroutine mpi_debug_out


  ! ---------------------------------------------------------
  subroutine mpi_time_accum(index, sec, sec_diff)
    integer, intent(in)  :: index
    real(r8), intent(in)  :: sec
    real(r8), intent(out) :: sec_diff

    sec_diff         = sec - sec_in
    sec_accum(index) = sec_accum(index) + sec_diff

  end subroutine mpi_time_accum

end module mpi_debug_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
