!! Copyright (C) 2013 M. Oliveira
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

module partition_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use io_binary_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                         &
    partition_t,                    &
    partition_init,                 &
    partition_end,                  &
    partition_set,                  &
    partition_dump,                 &
    partition_load,                 &
    partition_get_local_size,       &
    partition_get_global,           &
    partition_get_partition_number, &
    partition_get_np_local,         &
    partition_get_npart,            &
    partition_get_part,             &
    partition_get_local


  !> The partition is an array that contains the mapping between some global index
  !! and a process, such that point ip will be stored in process partition%part(ip).
  !!
  !! In this module this array is distributed among the processes, such that each process
  !! only stores a portion of the full array. Because each process needs to know in a
  !! straighforward way (i.e. without having to perform any kind of communication or
  !! lengthy operations) which process stores the partition corresponding to any giving
  !! point, the distribution of the points is done using a block data decomposition:
  !! - partition i: starts at floor((i-1) * np_global/npart) + 1
  !! - partition i: ends at   floor( i * np_global/npart)
  !! - element j stored on partition: floor(npart*(j+1)/np_global) + 1
  !!
  !! Note 1: this module can be a bit confusing as they are in fact two partitions. One is the
  !! partition of some array (in the case of Octopus, this is typically the mesh functions),
  !! which is the main information stored in the partition_t object, and then there is the
  !! partition of the partition itself, as this is also distributed.
  !!
  !! Note 2: in principle, the mpi group used by the processes for the partition distribution
  !! does not need to be the same as the mpi group of the processes the partition refers to.
  type partition_t
    private

    !> The following components are the same for all processes:
    type(mpi_grp_t) :: mpi_grp   !< The mpi group use for distributing the partition data.
    integer(i8) ::     np_global !< The total number of points in the partition.
    integer     ::     npart     !< The number of partitions.
    integer, allocatable     :: np_local_vec(:)  !< The number of points for each partition
    integer(i8), allocatable :: istart_vec(:)    !< The position of the first point for each partition

    !> The following components are process-dependent:
    integer :: partno            !< local partition number (i.e. rank + 1)
    integer :: np_local          !< The number of points of the partition stored in this process.
    integer(i8) :: istart        !< The position of the first point stored in this process.
    integer, allocatable :: part(:)  !< The local portion of the partition.

  end type partition_t

  interface partition_get_partition_number
    module procedure partition_get_partition_number_4, partition_get_partition_number_8
  end interface partition_get_partition_number


contains

  !---------------------------------------------------------
  subroutine partition_init(partition, np_global, mpi_grp)
    type(partition_t), intent(out) :: partition
    integer(i8),       intent(in)  :: np_global
    type(mpi_grp_t),   intent(in)  :: mpi_grp

    integer :: ipart
    integer(i8) :: iend

    PUSH_SUB(partition_init)

    !Global variables
    partition%mpi_grp = mpi_grp
    partition%np_global = np_global
    partition%npart = mpi_grp%size
    partition%partno = mpi_grp%rank + 1

    SAFE_ALLOCATE(partition%np_local_vec(1:partition%npart))
    SAFE_ALLOCATE(partition%istart_vec(1:partition%npart))

    ! use block data decomposition
    do ipart = 1, partition%npart
      ! according to Fortran 2008 standard, Sect. 7.1.9.3, point 4, the kind parameter of an
      ! operation is the one with a larger decimal range, i.e., everything is promoted to i8
      partition%istart_vec(ipart) = (ipart-1) * np_global/partition%npart + 1
      iend  = ipart * np_global/partition%npart
      partition%np_local_vec(ipart) = i8_to_i4(iend - partition%istart_vec(ipart) + 1)
    end do
    partition%istart = partition%istart_vec(partition%partno)
    partition%np_local = partition%np_local_vec(partition%partno)

    !Allocate memory for the partition
    SAFE_ALLOCATE(partition%part(1:partition%np_local))

    POP_SUB(partition_init)
  end subroutine partition_init

  ! ---------------------------------------------------------
  subroutine partition_end(partition)
    type(partition_t), intent(inout) :: partition

    PUSH_SUB(partition_end)

    SAFE_DEALLOCATE_A(partition%part)

    POP_SUB(partition_end)
  end subroutine partition_end

  ! ---------------------------------------------------------
  subroutine partition_set(partition, part)
    type(partition_t), intent(inout) :: partition
    integer,           intent(in)    :: part(:) !< The local portion of the partition.

    PUSH_SUB(partition_set)

    partition%part(1:partition%np_local) = part(1:partition%np_local)

    POP_SUB(partition_set)
  end subroutine partition_set

  ! ---------------------------------------------------------
  !> Partition is written in parallel if MPI2 is available.
  !! Otherwise, is gathered in root and then written.
  subroutine partition_dump(partition, dir, filename, ierr)
    type(partition_t), intent(in)  :: partition
    character(len=*),  intent(in)  :: dir
    character(len=*),  intent(in)  :: filename
    integer,           intent(out) :: ierr

    integer :: err
    character(len=MAX_PATH_LEN) :: full_filename

    PUSH_SUB(partition_dump)

    ierr = 0
    full_filename = trim(dir)//'/'//trim(filename)

    ! Write the header (root only) and wait
    if (mpi_grp_is_root(partition%mpi_grp)) then
      call iwrite_header(full_filename, partition%np_global, err)
      if (err /= 0) ierr = ierr + 1
    end if
    call partition%mpi_grp%bcast(ierr, 1, MPI_INTEGER, 0)

    ASSERT(all(partition%part(:) > 0))

    ! Each process writes a portion of the partition
    if (ierr == 0) then
      call mpi_debug_in(partition%mpi_grp%comm, C_MPI_FILE_WRITE)
      ! Only one rank per partition group should write the partition restart information
      ! Otherwise, more than once is trying to be written data
      if (mod(mpi_world%rank, mpi_world%size/partition%mpi_grp%size) == 0) then
        call io_binary_write_parallel(full_filename, partition%mpi_grp%comm, partition%istart, &
          partition%np_local, partition%part, err)
        call mpi_debug_out(partition%mpi_grp%comm, C_MPI_FILE_WRITE)
        if (err /= 0) ierr = ierr + 2
      end if
    end if

    POP_SUB(partition_dump)
  end subroutine partition_dump

  ! ---------------------------------------------------------
  !> Partition is read in parallel if MPI2 is available.
  !! Otherwise, is read by the root and then scattered
  subroutine partition_load(partition, dir, filename, ierr)
    type(partition_t), intent(inout) :: partition
    character(len=*),  intent(in)    :: dir
    character(len=*),  intent(in)    :: filename
    integer,           intent(out)   :: ierr

    integer :: err
    integer(i8) :: np, file_size
    integer, allocatable :: scounts(:)
    integer(i8), allocatable :: sdispls(:)
    character(len=MAX_PATH_LEN) :: full_filename

    PUSH_SUB(partition_load)

    ierr = 0
    full_filename = trim(dir)//'/'//trim(filename)

    ! This is a writing to avoid an optimization of gfortran with -O3
    write(message(1),'(a,i8)') "Info: number of points in the partition (in root process) =", size(partition%part)
    call messages_info(1)

    ! Check if the file exists and has the proper size (only world root)
    if (mpi_world%rank == 0) then
      call io_binary_get_info(full_filename, np, file_size, err)
      ASSERT(np == partition%np_global)
    end if

    ! All nodes need to know the result
    call mpi_world%bcast(err, 1, MPI_INTEGER, 0)
    call mpi_world%bcast(file_size, 1, MPI_INTEGER8, 0)

    if (err /= 0) then
      ierr = ierr + 1
      POP_SUB(partition_load)
      return
    end if
    ! The size of the file is not as big as np_global
    if (file_size - 64 /= partition%np_global * FC_INTEGER_SIZE) then
      ierr = ierr + 2
      POP_SUB(partition_load)
      return
    end if

    ! Calculate displacements for reading
    SAFE_ALLOCATE(scounts(1:partition%npart))
    SAFE_ALLOCATE(sdispls(1:partition%npart))

    scounts = partition%np_local_vec
    sdispls = partition%istart_vec - 1

    ASSERT(sum(scounts(:)) == partition%np_global)

#ifdef HAVE_MPI
    call mpi_debug_in(partition%mpi_grp%comm, C_MPI_FILE_READ)
    call io_binary_read_parallel(full_filename, partition%mpi_grp%comm, partition%istart, &
      partition%np_local, partition%part, err)
    call mpi_debug_out(partition%mpi_grp%comm, C_MPI_FILE_READ)
    if (err /= 0) ierr = ierr + 4
#endif

    if (any(partition%part(:) <= 0)) then
      write(message(1),'(a)') 'Internal error: some elements of partition are <= 0.'
      write(message(2),*) 'filename = ', full_filename
      write(message(3),*) 'scounts = ', scounts(:)
      write(message(4),*) 'sdispls = ', sdispls(:)
      call messages_fatal(6)
    end if

    SAFE_DEALLOCATE_A(scounts)
    SAFE_DEALLOCATE_A(sdispls)

    POP_SUB(partition_load)
  end subroutine partition_load

  ! ---------------------------------------------------------
  subroutine partition_get_local_size(partition, istart, np_local)
    type(partition_t), intent(in)  :: partition
    integer(i8),       intent(out) :: istart   !< The number of points of the partition stored in this process.
    integer,           intent(out) :: np_local !< The position of the first point stored in this process.

    PUSH_SUB(partition_get_local_size)

    istart = partition%istart
    np_local = partition%np_local

    POP_SUB(partition_get_local_size)
  end subroutine partition_get_local_size

  ! ---------------------------------------------------------
  !> Returns the global partition. If root is present, the partition is
  !! gathered only in that node. Otherwise it is gathered in all nodes.
  subroutine partition_get_global(partition, part_global, root)
    type(partition_t), intent(in)  :: partition
    integer,           intent(out) :: part_global(:)
    integer, optional, intent(in)  :: root

    integer(i8), allocatable :: rdispls(:)
    integer, allocatable :: rcounts(:)

    PUSH_SUB(partition_get_global)

    SAFE_ALLOCATE(rdispls(1:partition%npart))
    SAFE_ALLOCATE(rcounts(1:partition%npart))

    rcounts = partition%np_local_vec
    rdispls = partition%istart_vec - 1

    ASSERT(all(partition%part(1:partition%np_local) > 0))

    if (present(root)) then
      call partition%mpi_grp%gatherv(partition%part, partition%np_local, MPI_INTEGER, &
        part_global, rcounts, rdispls, MPI_INTEGER, root)
    else
      call partition%mpi_grp%allgatherv(partition%part, partition%np_local, MPI_INTEGER, &
        part_global, rcounts, rdispls, MPI_INTEGER)
    end if

    if (.not. present(root) .or. partition%mpi_grp%rank == 0) then
      ASSERT(all(part_global(:) > 0))
    end if

#ifndef HAVE_MPI
    part_global = partition%part
#endif

    SAFE_DEALLOCATE_A(rdispls)
    SAFE_DEALLOCATE_A(rcounts)

    POP_SUB(partition_get_global)
  end subroutine partition_get_global

  ! ---------------------------------------------------------
  !> Given a list of _global_ indices, return the partition number
  !! where those points are stored.
  !! Note that this routine will accept global indices equal to 0. In that
  !! case it will return 0 as a partition number.
  subroutine partition_get_partition_number_4(partition, np, points, partno)
    type(partition_t), intent(in)  :: partition
    integer,           intent(in)  :: np
    integer(i8),       intent(in)  :: points(:)
    integer,           intent(out) :: partno(:)

    integer :: ip, nproc, rnp
    integer(i8), allocatable :: sbuffer(:), rbuffer(:)
    integer, allocatable :: scounts(:), rcounts(:)
    integer, allocatable :: sdispls(:), rdispls(:)
    integer, allocatable :: ipos(:), order(:)

    PUSH_SUB(partition_get_partition_number_4)

    SAFE_ALLOCATE(scounts(1:partition%npart))
    SAFE_ALLOCATE(rcounts(1:partition%npart))
    SAFE_ALLOCATE(sdispls(1:partition%npart))
    SAFE_ALLOCATE(rdispls(1:partition%npart))

    ! How many points will we have to send/receive from each process?
    scounts = 0
    do ip = 1, np
      ! Who knows where points(ip) is stored?
      nproc = partition_get_number(partition, points(ip))
      ! We increase the respective counter
      scounts(nproc) = scounts(nproc) + 1
    end do


    !Tell each process how many points we will need from it
    call partition%mpi_grp%alltoall(scounts, 1, MPI_INTEGER, &
      rcounts, 1, MPI_INTEGER)

    !Build displacement arrays
    sdispls(1) = 0
    rdispls(1) = 0
    do ip = 2, partition%npart
      sdispls(ip) = sdispls(ip-1) + scounts(ip-1)
      rdispls(ip) = rdispls(ip-1) + rcounts(ip-1)
    end do


    rnp = sum(rcounts)
    SAFE_ALLOCATE(sbuffer(1:np))
    SAFE_ALLOCATE(rbuffer(1:rnp))

    !Put points in correct order for sending
    SAFE_ALLOCATE(ipos(1:partition%npart))
    SAFE_ALLOCATE(order(1:np))
    ipos = 0
    do ip = 1, np
      ! Who knows where points(ip) is stored?
      nproc = partition_get_number(partition, points(ip))

      !We increase the respective counter
      ipos(nproc) = ipos(nproc) + 1

      !Put the point in the correct place
      order(ip) = sdispls(nproc) + ipos(nproc)
      sbuffer(order(ip)) = points(ip) ! global index of the point
    end do
    SAFE_DEALLOCATE_A(ipos)

    !Send the global indices of the points to the process that knows what is the corresponding partition
    call partition%mpi_grp%alltoallv(sbuffer, scounts, sdispls, MPI_INTEGER8, &
      rbuffer, rcounts, rdispls, MPI_INTEGER8)

    !We get the partition number from the global index. This will be send back.
    do ip = 1, rnp
      if (rbuffer(ip) == 0) cycle
      rbuffer(ip) = partition%part(rbuffer(ip) - partition%istart + 1)
    end do

    !Now we send the information backwards
    call partition%mpi_grp%alltoallv(rbuffer, rcounts, rdispls, MPI_INTEGER8, &
      sbuffer, scounts, sdispls, MPI_INTEGER8)

    !Reorder the points
    do ip = 1, np
      partno(ip) = i8_to_i4(sbuffer(order(ip)))
    end do

    !Deallocate memory
    SAFE_DEALLOCATE_A(order)
    SAFE_DEALLOCATE_A(sbuffer)
    SAFE_DEALLOCATE_A(scounts)
    SAFE_DEALLOCATE_A(sdispls)
    SAFE_DEALLOCATE_A(rbuffer)
    SAFE_DEALLOCATE_A(rcounts)
    SAFE_DEALLOCATE_A(rdispls)

    POP_SUB(partition_get_partition_number_4)
  end subroutine partition_get_partition_number_4

  !---------------------------------------------------------
  subroutine partition_get_partition_number_8(partition, np, points, partno)
    type(partition_t), intent(in)  :: partition
    integer(i8),       intent(in)  :: np
    integer(i8),       intent(in)  :: points(:)
    integer,           intent(out) :: partno(:)

    integer(i8) :: rounds, offset, iround

    PUSH_SUB(partition_get_partition_number_8)

    rounds = np/huge(0_i4)
    do iround = 1, rounds
      offset = (iround - 1)*huge(0_i4) + 1
      call partition_get_partition_number(partition, huge(0_i4), points(offset:offset+huge(0_i4)), &
        partno(offset:offset+huge(0_i4)))
    end do
    offset = rounds*huge(0_i4) + 1
    call partition_get_partition_number(partition, i8_to_i4(np-offset+1), points(offset:np), &
      partno(offset:np))

    POP_SUB(partition_get_partition_number_8)
  end subroutine partition_get_partition_number_8

  !---------------------------------------------------------
  !> Given the partition, returns the corresponding number of local
  !! points that each partition has.
  subroutine partition_get_np_local(partition, np_local_vec)
    type(partition_t),    intent(in)  :: partition       !< Current partition
    integer,              intent(inout) :: np_local_vec(:) !< Vector of local points (np_local)

    integer, allocatable :: np_local_vec_tmp(:)
    integer :: ip

    PUSH_SUB(partition_get_np_local)

    ASSERT(ubound(np_local_vec, 1) >= partition%npart)
    ASSERT(partition%npart > 0)
    ASSERT(all(partition%part(:) > 0))
    SAFE_ALLOCATE(np_local_vec_tmp(1:partition%npart))
    np_local_vec_tmp = 0

    ! Calculate locally the local points of each partition
    do ip = 1, partition%np_local
      np_local_vec_tmp(partition%part(ip)) = np_local_vec_tmp(partition%part(ip)) + 1
    end do

    ! Collect all the local points
    call partition%mpi_grp%allreduce(np_local_vec_tmp, np_local_vec, partition%npart, MPI_INTEGER, MPI_SUM)
    SAFE_DEALLOCATE_A(np_local_vec_tmp)

    POP_SUB(partition_get_np_local)
  end subroutine partition_get_np_local

  !---------------------------------------------------------
  !> Returns the total number of partitions
  pure integer function partition_get_npart(partition) result(npart)
    type(partition_t), intent(in) :: partition
    npart = partition%npart
  end function partition_get_npart

  !---------------------------------------------------------
  !> Returns the partition of the local point
  pure integer function partition_get_part(partition, local_point) result(part)
    type(partition_t), intent(in) :: partition
    integer,           intent(in) :: local_point
    part = partition%part(local_point)
  end function partition_get_part

  !---------------------------------------------------------
  !> Returns the partition number for a given global index
  !> If the index is zero, return local partition
  pure integer function partition_get_number(partition, global_point) result(part)
    type(partition_t), intent(in) :: partition
    integer(i8),       intent(in) :: global_point

    if (global_point == 0) then
      part = partition%partno
    else
      part = i8_to_i4((partition%npart*global_point - 1)/partition%np_global + 1)
    end if
  end function partition_get_number

  !---------------------------------------------------------
  !> Calculates the local vector of all partitions in parallel.
  !! Local vector stores the global point indices that each partition
  !! has.
  subroutine partition_get_local(partition, rbuffer, np_local)
    type(partition_t),    intent(in)    :: partition
    integer(i8),          intent(inout) :: rbuffer(:) !< The actual result, the local vector from 1 to np_local
    integer,              intent(out)   :: np_local   !< Number of elements, might be less than partition%np_local

    integer :: ip, ipart
    integer(i8) :: istart
    integer(i8), allocatable :: sdispls(:), rdispls(:), sbuffer(:)
    integer, allocatable :: scounts(:), rcounts(:)

    PUSH_SUB(partition_get_local)

    SAFE_ALLOCATE(sdispls(1:partition%npart))
    SAFE_ALLOCATE(scounts(1:partition%npart))
    SAFE_ALLOCATE(rcounts(1:partition%npart))

    ! Calculate the starting point of the running process
    istart = partition%istart - 1

    scounts = 0
    ! Count and store the local points for each partition
    do ip = 1, partition%np_local
      ipart = partition%part(ip)
      scounts(ipart) = scounts(ipart) + 1
    end do

    ! Create displacements
    sdispls(1) = 0
    do ipart = 2, partition%npart
      sdispls(ipart) = sdispls(ipart-1) + scounts(ipart-1)
    end do

    ! Allocate and fill the send buffer
    np_local = sum(scounts)
    scounts = 0
    SAFE_ALLOCATE(sbuffer(1:np_local))
    do ip = 1, np_local
      ipart = partition%part(ip)
      scounts(ipart) = scounts(ipart) + 1
      sbuffer(sdispls(ipart) + scounts(ipart)) = ip + istart
    end do

    ! Tell each process how many points we will need from it
    call partition%mpi_grp%alltoall(scounts, 1, MPI_INTEGER, &
      rcounts, 1, MPI_INTEGER)

    ! Create the displacement vector from the counts vector
    np_local = sum(rcounts)
    SAFE_ALLOCATE(rdispls(1:partition%npart))
    ASSERT(ubound(rbuffer, 1) >= np_local)

    rdispls(1) = 0
    do ipart = 2, partition%npart
      rdispls(ipart) = rcounts(ipart-1) + rdispls(ipart-1)
    end do

    ! Collect the corresponding points to each process
    call partition%mpi_grp%alltoallv(sbuffer, scounts, sdispls, MPI_INTEGER8, &
      rbuffer, rcounts, rdispls, MPI_INTEGER8)

    SAFE_DEALLOCATE_A(sdispls)
    SAFE_DEALLOCATE_A(scounts)
    SAFE_DEALLOCATE_A(sbuffer)
    SAFE_DEALLOCATE_A(rcounts)
    SAFE_DEALLOCATE_A(rdispls)

    POP_SUB(partition_get_local)
  end subroutine partition_get_local
end module partition_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
