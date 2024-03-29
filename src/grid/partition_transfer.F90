!! Copyright (C) 2011-2012 M. Oliveira
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

module partition_transfer_oct_m
  use debug_oct_m
  use global_oct_m
  use iihash_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                      &
    partition_transfer_t,        &
    partition_transfer_init,     &
    partition_transfer_end,      &
    dpartition_transfer,         &
    zpartition_transfer

  type partition_transfer_t
    private
    type(mpi_grp_t) :: mpi_grp

    integer, allocatable :: rdispls(:)
    integer, allocatable :: sdispls(:)
    integer, allocatable :: rcounts(:)
    integer, allocatable :: scounts(:)
  end type partition_transfer_t

  type(profile_t), save :: prof_transfer

contains

  ! -----------------------------------------------------------------
  !> \warning  Input and output groups may have a different number of
  !! processes. In that case the transfer will only work if some
  !! further constraints are met (see sanity checks below). One of the
  !! cases where it should work is when one of the groups is a subgroup
  !! of a Cartesian topology created from the other group. This is the
  !! case when one of the groups is mpi_world, the other is the
  !! parallelization in domains group, and there are no slaves.
  !! If the optional argument inverse is set to .true., the direction of
  !! the transfer is inverse: part_out is assumed to specify the partition
  !! on the incoming mpi group. This is useful if one needs to get points
  !! defined on the output group.
  subroutine partition_transfer_init(this, np, global_index, mpi_grp_in, mpi_grp_out, &
    part_out, nsend, nrec, order_in, order_out, inverse)
    type(partition_transfer_t), intent(out) :: this
    integer,                    intent(in)  :: np !< the number of local points in the input partition
    integer(i8),                intent(in)  :: global_index(:) !< the global indices of the points of the input partition
    type(mpi_grp_t), target,    intent(in)  :: mpi_grp_in
    type(mpi_grp_t), target,    intent(in)  :: mpi_grp_out
    integer,                    intent(in)  :: part_out(:) !< point -> partition
    integer,                    intent(out) :: nsend
    integer,                    intent(out) :: nrec
    integer(i8), allocatable,   intent(out) :: order_in(:)
    integer(i8), allocatable,   intent(out) :: order_out(:)
    logical, optional,          intent(in)  :: inverse

    logical :: found, inverse_
    integer :: n12, tmp_partno(2), ipart, opart, ip, pcount, mycolumn, irec, isend, ipos
    type(iihash_t) :: map_out
    type(mpi_grp_t), pointer :: grp1, grp2
    integer, allocatable :: partno_list(:,:), part_map(:,:)
    type(profile_t), save :: prof

    PUSH_SUB(partition_transfer_init)
    call profiling_in(prof,"P_TRANS_INIT")

    inverse_ = optional_default(inverse, .false.)

    ! In order to avoid unnecessary communications, all the data
    ! transfer is going to be made from the point of view of the group
    ! that has more processes.
    if (mpi_grp_in%size >= mpi_grp_out%size) then
      grp1 => mpi_grp_in
      grp2 => mpi_grp_out
    else
      grp1 => mpi_grp_out
      grp2 => mpi_grp_in
    end if
    call mpi_grp_copy(this%mpi_grp, grp1)

    ! The number of partitions in group 1 should be a multiple of the
    ! number of partitions in group 2.
    if (mod(grp1%size, grp2%size) /= 0) then
      message(1) = "Incompatible size of mpi groups in partition_transfer_init"
      call messages_fatal(1)
    end if
    n12 = grp1%size/grp2%size

    ! We need to know the partition number of all the processes in
    ! both groups
    SAFE_ALLOCATE(partno_list(1:2, 1:grp1%size))
    tmp_partno(1) = grp1%rank + 1
    tmp_partno(2) = grp2%rank + 1
    call this%mpi_grp%allgather(tmp_partno(1), 2, MPI_INTEGER, partno_list(1, 1), 2, MPI_INTEGER)

    ! Build partition map. This is a matrix with n12 columns and
    ! grp2%size lines. Each line contains the partition numbers of the
    ! processes of group 1 that also store the partition of group 2
    ! with the same number as the line number. The number of columns
    ! for each line should be exactly equal to grp1%size/grp2%size in
    ! all cases and there should be no repeated values.
    SAFE_ALLOCATE(part_map(1:grp2%size, 1:n12))
    part_map = 0
    do ipart = 1, grp2%size
      pcount = 0
      do ip = 1, grp1%size
        if (partno_list(2, ip) == ipart) then
          pcount = pcount + 1

          if (pcount > n12 .or. any(partno_list(1, ip) == part_map(1:ipart,:))) then
            message(1) = "Incompatible mpi groups in partition_transfer_init"
            call messages_fatal(1)
          end if
          part_map(ipart, pcount) = partno_list(1, ip)
          if (ip == grp1%rank + 1) mycolumn = pcount
        end if
      end do
      if (pcount /= n12) then
        message(1) = "Incompatible mpi groups in partition_transfer_init"
        call messages_fatal(1)
      end if
    end do

    ! Build mapping between all the possible receivers and the ouput
    ! group. This map is a hash table, where the keys are the possible
    ! receivers and the values are the output partition these
    ! receivers are responsible for.
    ! If group 1 is the input group, then all members of group 1 are
    ! possible receivers.
    ! If group 1 is the output group, then, in order to avoid
    ! unnecessary communications, each process will only send data to
    ! a subset of all the possible receivers.  We will choose the
    ! processes that are on the same column of the partition map than
    ! the local process. Note that this implies that there are always
    ! mpi_grp_in%size possible receivers.
    call iihash_init(map_out)
    if (.not. inverse_) then
      do ipart = 1, grp2%size
        if (mpi_grp_in%size >= mpi_grp_out%size) then
          do ip = 1, n12
            call iihash_insert(map_out, part_map(ipart, ip), ipart)
          end do
        else
          call iihash_insert(map_out, part_map(ipart, mycolumn), part_map(ipart, mycolumn))
        end if
      end do
    else
      do ipart = 1, grp2%size
        if (mpi_grp_in%size >= mpi_grp_out%size) then
          do ip = 1, n12
            call iihash_insert(map_out, part_map(ipart, ip), part_map(ipart, ip))
          end do
        else
          call iihash_insert(map_out, part_map(ipart, mycolumn), ipart)
        end if
      end do
    end if

    if (.not. inverse_) then
      ! Total number of points to be sent
      nsend = 0
      do irec = 1, grp1%size
        opart = iihash_lookup(map_out, irec, found)
        if (.not. found) cycle
        nsend = nsend + count(part_out(1:np) == opart)
      end do
    else
      ! Total number of points to be received
      nrec = 0
      do isend = 1, grp1%size
        opart = iihash_lookup(map_out, isend, found)
        if (.not. found) cycle
        nrec = nrec + count(part_out(1:np) == opart)
      end do
    end if

    ! List of points to be send
    SAFE_ALLOCATE(this%sdispls(1:grp1%size))
    SAFE_ALLOCATE(this%scounts(1:grp1%size))
    ! Displacements and number of points to be received
    SAFE_ALLOCATE(this%rdispls(1:grp1%size))
    SAFE_ALLOCATE(this%rcounts(1:grp1%size))

    if (.not. inverse_) then
      SAFE_ALLOCATE(order_in(1:max(1,nsend)))

      ipos = 0
      ! Loop over all possible receivers
      do irec = 1, grp1%size
        this%scounts(irec) = 0
        this%sdispls(irec) = ipos

        opart = iihash_lookup(map_out, irec, found)
        if (.not. found) cycle

        do ip = 1, np
          ! Should point ip be sent to partition opart?
          if (part_out(ip) == opart) then
            ipos = ipos + 1
            order_in(ipos) = global_index(ip)
            this%scounts(irec) = this%scounts(irec) + 1
          end if
        end do

      end do
    else
      SAFE_ALLOCATE(order_out(1:max(1,nrec)))

      ipos = 0
      ! Loop over all possible senders
      do isend = 1, grp1%size
        this%rcounts(isend) = 0
        this%rdispls(isend) = ipos

        opart = iihash_lookup(map_out, isend, found)
        if (.not. found) cycle

        do ip = 1, np
          ! Should point ip be received from partition opart?
          if (part_out(ip) == opart) then
            ipos = ipos + 1
            order_out(ipos) = global_index(ip)
            this%rcounts(isend) = this%rcounts(isend) + 1
          end if
        end do

      end do
    end if

    SAFE_DEALLOCATE_A(part_map)
    SAFE_DEALLOCATE_A(partno_list)
    call iihash_end(map_out)

    if (.not. inverse_) then
      ! assemble the receive information
      call this%mpi_grp%alltoall(this%scounts, 1, MPI_INTEGER, &
        this%rcounts, 1, MPI_INTEGER)
      nrec = 0
      do isend = 1, grp1%size
        this%rdispls(isend) = nrec
        nrec = nrec + this%rcounts(isend)
      end do
    else
      ! assemble the send information
      call this%mpi_grp%alltoall(this%rcounts, 1, MPI_INTEGER, &
        this%scounts, 1, MPI_INTEGER)
      nsend = 0
      do irec = 1, grp1%size
        this%sdispls(irec) = nsend
        nsend = nsend + this%scounts(irec)
      end do
    end if

    if (.not. inverse_) then
      ! Ordering of the points at output
      SAFE_ALLOCATE(order_out(1:max(1,nrec)))
      call this%mpi_grp%alltoallv(order_in, this%scounts, this%sdispls, MPI_INTEGER8, &
        order_out, this%rcounts, this%rdispls, MPI_INTEGER8)
    else
      ! Ordering of the points at input
      SAFE_ALLOCATE(order_in(1:max(1,nsend)))
      call this%mpi_grp%alltoallv(order_out, this%rcounts, this%rdispls, MPI_INTEGER8, &
        order_in, this%scounts, this%sdispls, MPI_INTEGER8)
    end if

    call profiling_out(prof)
    POP_SUB(partition_transfer_init)
  end subroutine partition_transfer_init

  ! -----------------------------------------------------------------
  subroutine partition_transfer_end(this)
    type(partition_transfer_t), intent(inout) :: this

    PUSH_SUB(partition_transfer_end)

    SAFE_DEALLOCATE_A(this%rdispls)
    SAFE_DEALLOCATE_A(this%sdispls)
    SAFE_DEALLOCATE_A(this%rcounts)
    SAFE_DEALLOCATE_A(this%scounts)

    POP_SUB(partition_transfer_end)
  end subroutine partition_transfer_end

#include "undef.F90"
#include "real.F90"
#include "partition_transfer_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "partition_transfer_inc.F90"

end module partition_transfer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
