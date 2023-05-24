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

#include "global.h"

module mpi_test_oct_m
  use global_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  implicit none
  private

  public :: test_mpiwrappers

contains

  subroutine test_mpiwrappers
    if (test_scatterv()) then
      write(message(1), *) "scatterv: PASS"
    else
      write(message(1), *) "scatterv: FAIL"
    end if
    call messages_info(1)

    if (test_gatherv()) then
      write(message(1), *) "gatherv: PASS"
    else
      write(message(1), *) "gatherv: FAIL"
    end if
    call messages_info(1)

    if (test_allgatherv()) then
      write(message(1), *) "allgatherv: PASS"
    else
      write(message(1), *) "allgatherv: FAIL"
    end if
    call messages_info(1)

    if (test_alltoallv()) then
      write(message(1), *) "alltoallv: PASS"
    else
      write(message(1), *) "alltoallv: FAIL"
    end if
    call messages_info(1)
  end subroutine test_mpiwrappers

  logical function test_scatterv()
    FLOAT, allocatable :: sendbuf(:), irecvbuf(:), lrecvbuf(:)
    integer, allocatable :: sendcnts(:), recvcounts(:)
    integer(i4), allocatable :: isdispls(:)
    integer(i8), allocatable :: lsdispls(:)
    integer, parameter :: N = 10003
    integer :: ii
    logical :: equal, allequal

    SAFE_ALLOCATE(sendcnts(1:mpi_world%size))
    SAFE_ALLOCATE(isdispls(1:mpi_world%size))
    SAFE_ALLOCATE(lsdispls(1:mpi_world%size))
    SAFE_ALLOCATE(recvcounts(1:1))
    do ii = 1, mpi_world%size
      isdispls(ii) = (ii-1) * N / mpi_world%size
      lsdispls(ii) = (ii-1) * N / mpi_world%size
      sendcnts(ii) = ii * N / mpi_world%size - isdispls(ii)
    end do
    if (mpi_world%rank == 0) then
      SAFE_ALLOCATE(sendbuf(1:N))
      do ii = 1, N
        sendbuf(ii) = ii
      end do
    end if
    recvcounts(1) = sendcnts(mpi_world%rank+1)
    SAFE_ALLOCATE(irecvbuf(1:recvcounts(1)))
    SAFE_ALLOCATE(lrecvbuf(1:recvcounts(1)))

    call mpi_world%scatterv(sendbuf, sendcnts, isdispls, MPI_FLOAT, irecvbuf, recvcounts(1), MPI_FLOAT, 0)
    call mpi_world%scatterv(sendbuf, sendcnts, lsdispls, MPI_FLOAT, lrecvbuf, recvcounts(1), MPI_FLOAT, 0)

    equal = all(irecvbuf == lrecvbuf)
    call mpi_world%allreduce(equal, allequal, 1, MPI_LOGICAL, MPI_LAND)
    test_scatterv = allequal

    SAFE_DEALLOCATE_A(sendcnts)
    SAFE_DEALLOCATE_A(isdispls)
    SAFE_DEALLOCATE_A(lsdispls)
    SAFE_DEALLOCATE_A(recvcounts)
    if (mpi_world%rank == 0) then
      SAFE_DEALLOCATE_A(sendbuf)
    end if
    SAFE_DEALLOCATE_A(irecvbuf)
    SAFE_DEALLOCATE_A(lrecvbuf)
  end function test_scatterv
    
  logical function test_gatherv()
    FLOAT, allocatable :: sendbuf(:), irecvbuf(:), lrecvbuf(:)
    integer, allocatable :: sendcnts(:), recvcounts(:)
    integer(i4), allocatable :: irdispls(:)
    integer(i8), allocatable :: lrdispls(:)
    integer, parameter :: N = 10003
    integer :: ii
    logical :: equal

    SAFE_ALLOCATE(recvcounts(1:mpi_world%size))
    SAFE_ALLOCATE(irdispls(1:mpi_world%size))
    SAFE_ALLOCATE(lrdispls(1:mpi_world%size))
    SAFE_ALLOCATE(sendcnts(1:1))
    do ii = 1, mpi_world%size
      irdispls(ii) = (ii-1) * N / mpi_world%size
      lrdispls(ii) = (ii-1) * N / mpi_world%size
      recvcounts(ii) = ii * N / mpi_world%size - irdispls(ii)
    end do
    if (mpi_world%rank == 0) then
      SAFE_ALLOCATE(irecvbuf(1:N))
      SAFE_ALLOCATE(lrecvbuf(1:N))
    end if
    sendcnts(1) = recvcounts(mpi_world%rank+1)
    SAFE_ALLOCATE(sendbuf(1:sendcnts(1)))
    do ii = 1, sendcnts(1)
      sendbuf(ii) = ii
    end do

    call mpi_world%gatherv(sendbuf, sendcnts(1), MPI_FLOAT, irecvbuf, recvcounts, irdispls, MPI_FLOAT, 0)
    call mpi_world%gatherv(sendbuf, sendcnts(1), MPI_FLOAT, lrecvbuf, recvcounts, lrdispls, MPI_FLOAT, 0)

    if (mpi_world%rank == 0) then
      equal = all(irecvbuf == lrecvbuf)
    end if
    test_gatherv = equal

    SAFE_DEALLOCATE_A(recvcounts)
    SAFE_DEALLOCATE_A(irdispls)
    SAFE_DEALLOCATE_A(lrdispls)
    SAFE_DEALLOCATE_A(sendcnts)
    if (mpi_world%rank == 0) then
      SAFE_DEALLOCATE_A(irecvbuf)
      SAFE_DEALLOCATE_A(lrecvbuf)
    end if
    SAFE_DEALLOCATE_A(sendbuf)
  end function test_gatherv

  logical function test_allgatherv()
    FLOAT, allocatable :: sendbuf(:), irecvbuf(:), lrecvbuf(:)
    integer, allocatable :: sendcnts(:), recvcounts(:)
    integer(i4), allocatable :: irdispls(:)
    integer(i8), allocatable :: lrdispls(:)
    integer, parameter :: N = 10003
    integer :: ii
    logical :: equal, allequal

    SAFE_ALLOCATE(recvcounts(1:mpi_world%size))
    SAFE_ALLOCATE(irdispls(1:mpi_world%size))
    SAFE_ALLOCATE(lrdispls(1:mpi_world%size))
    SAFE_ALLOCATE(sendcnts(1:1))
    do ii = 1, mpi_world%size
      irdispls(ii) = (ii-1) * N / mpi_world%size
      lrdispls(ii) = (ii-1) * N / mpi_world%size
      recvcounts(ii) = ii * N / mpi_world%size - irdispls(ii)
    end do
    SAFE_ALLOCATE(irecvbuf(1:N))
    SAFE_ALLOCATE(lrecvbuf(1:N))
    sendcnts(1) = recvcounts(mpi_world%rank+1)
    SAFE_ALLOCATE(sendbuf(1:sendcnts(1)))
    do ii = 1, sendcnts(1)
      sendbuf(ii) = ii
    end do

    call mpi_world%allgatherv(sendbuf, sendcnts(1), MPI_FLOAT, irecvbuf, recvcounts, irdispls, MPI_FLOAT)
    call mpi_world%allgatherv(sendbuf, sendcnts(1), MPI_FLOAT, lrecvbuf, recvcounts, lrdispls, MPI_FLOAT)

    equal = all(irecvbuf == lrecvbuf)
    call mpi_world%allreduce(equal, allequal, 1, MPI_LOGICAL, MPI_LAND)
    test_allgatherv = allequal

    SAFE_DEALLOCATE_A(recvcounts)
    SAFE_DEALLOCATE_A(irdispls)
    SAFE_DEALLOCATE_A(lrdispls)
    SAFE_DEALLOCATE_A(sendcnts)
    SAFE_DEALLOCATE_A(irecvbuf)
    SAFE_DEALLOCATE_A(lrecvbuf)
    SAFE_DEALLOCATE_A(sendbuf)
  end function test_allgatherv

  logical function test_alltoallv()
    FLOAT, allocatable :: sendbuf(:), irecvbuf(:), lrecvbuf(:)
    integer, allocatable :: sendcnts(:), recvcounts(:)
    integer(i4), allocatable :: isdispls(:), irdispls(:)
    integer(i8), allocatable :: lsdispls(:), lrdispls(:)
    integer, parameter :: N = 10003
    integer :: ii
    logical :: equal, allequal

    SAFE_ALLOCATE(sendcnts(1:mpi_world%size))
    SAFE_ALLOCATE(isdispls(1:mpi_world%size))
    SAFE_ALLOCATE(lsdispls(1:mpi_world%size))
    SAFE_ALLOCATE(irdispls(1:mpi_world%size))
    SAFE_ALLOCATE(lrdispls(1:mpi_world%size))
    SAFE_ALLOCATE(recvcounts(1:mpi_world%size))
    do ii = 1, mpi_world%size
      isdispls(ii) = (ii-1) * N / mpi_world%size
      lsdispls(ii) = (ii-1) * N / mpi_world%size
      sendcnts(ii) = ii * N / mpi_world%size - isdispls(ii)
    end do
    SAFE_ALLOCATE(sendbuf(1:N))
    do ii = 1, N
      sendbuf(ii) = ii
    end do
    call mpi_world%alltoall(sendcnts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER)
    SAFE_ALLOCATE(irecvbuf(1:sum(recvcounts)))
    SAFE_ALLOCATE(lrecvbuf(1:sum(recvcounts)))
    irdispls(1) = 0
    lrdispls(1) = 0
    do ii = 2, mpi_world%size
      irdispls(ii) = irdispls(ii-1) + recvcounts(ii-1)
      lrdispls(ii) = lrdispls(ii-1) + recvcounts(ii-1)
    end do

    call mpi_world%alltoallv(sendbuf, sendcnts, isdispls, MPI_FLOAT, &
      irecvbuf, recvcounts, irdispls, MPI_FLOAT)
    call mpi_world%alltoallv(sendbuf, sendcnts, lsdispls, MPI_FLOAT, &
      lrecvbuf, recvcounts, lrdispls, MPI_FLOAT)

    equal = all(irecvbuf == lrecvbuf)
    call mpi_world%allreduce(equal, allequal, 1, MPI_LOGICAL, MPI_LAND)
    test_alltoallv = allequal

    SAFE_DEALLOCATE_A(sendcnts)
    SAFE_DEALLOCATE_A(isdispls)
    SAFE_DEALLOCATE_A(lsdispls)
    SAFE_DEALLOCATE_A(irdispls)
    SAFE_DEALLOCATE_A(lrdispls)
    SAFE_DEALLOCATE_A(recvcounts)
    SAFE_DEALLOCATE_A(sendbuf)
    SAFE_DEALLOCATE_A(irecvbuf)
    SAFE_DEALLOCATE_A(lrecvbuf)
  end function test_alltoallv
end module mpi_test_oct_m
