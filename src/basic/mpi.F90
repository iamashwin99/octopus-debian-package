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

#include "global.h"

module mpi_oct_m
#if defined(MPI_MOD)
  use mpi
#endif
  use blacs_oct_m
  use mpi_debug_oct_m
  use kind_oct_m

  implicit none

  ! I do not make this module private on purpose, so that the symbols defined either in
  ! module mpi, or in mpif.h are exported

  ! some machines do not have a mpi module, but a mpi input file
#if defined(MPI_H)
  include "mpif.h"
#endif

  !> This is defined even when running serial
  type mpi_grp_t
    ! Components are public by default
    integer :: comm = -1 !< copy of the mpi communicator
    integer :: size = 0  !< size of comm (defined also in serial mode)
    integer :: rank = 0  !< rank of comm (defined also in serial mode)
  contains
    ! Wrapper functions for common MPI calls
    ! We do not check the error code in any of those wrappers because the behavior of
    ! an application is undefined after an MPI error according to the standard. The
    ! default is to let the application crash in such a case with an error message
    ! from the MPI runtime.
    procedure :: barrier => mpi_grp_barrier
    procedure :: dmpi_grp_scatterv, zmpi_grp_scatterv, impi_grp_scatterv, lmpi_grp_scatterv
    generic :: scatterv => dmpi_grp_scatterv, zmpi_grp_scatterv, impi_grp_scatterv, lmpi_grp_scatterv
    procedure :: dmpi_grp_scatterv_i8, zmpi_grp_scatterv_i8, impi_grp_scatterv_i8, lmpi_grp_scatterv_i8
    generic :: scatterv => dmpi_grp_scatterv_i8, zmpi_grp_scatterv_i8, impi_grp_scatterv_i8, lmpi_grp_scatterv_i8
    procedure :: dmpi_grp_gatherv, zmpi_grp_gatherv, impi_grp_gatherv, lmpi_grp_gatherv
    generic :: gatherv => dmpi_grp_gatherv, zmpi_grp_gatherv, impi_grp_gatherv, lmpi_grp_gatherv
    procedure :: dmpi_grp_gather_0, zmpi_grp_gather_0, impi_grp_gather_0, lmpi_grp_gather_0
    generic :: gather => dmpi_grp_gather_0, zmpi_grp_gather_0, impi_grp_gather_0, lmpi_grp_gather_0
    procedure :: dmpi_grp_gatherv_i8, zmpi_grp_gatherv_i8, impi_grp_gatherv_i8, lmpi_grp_gatherv_i8
    generic :: gatherv => dmpi_grp_gatherv_i8, zmpi_grp_gatherv_i8, impi_grp_gatherv_i8, lmpi_grp_gatherv_i8
    procedure :: dmpi_grp_alltoallv, zmpi_grp_alltoallv, impi_grp_alltoallv, lmpi_grp_alltoallv
    generic :: alltoallv => dmpi_grp_alltoallv, zmpi_grp_alltoallv, impi_grp_alltoallv, lmpi_grp_alltoallv
    procedure :: dmpi_grp_alltoallv_2, zmpi_grp_alltoallv_2, impi_grp_alltoallv_2, lmpi_grp_alltoallv_2
    generic :: alltoallv => dmpi_grp_alltoallv_2, zmpi_grp_alltoallv_2, impi_grp_alltoallv_2, lmpi_grp_alltoallv_2
    procedure :: dmpi_grp_alltoallv_3, zmpi_grp_alltoallv_3, impi_grp_alltoallv_3, lmpi_grp_alltoallv_3
    generic :: alltoallv => dmpi_grp_alltoallv_3, zmpi_grp_alltoallv_3, impi_grp_alltoallv_3, lmpi_grp_alltoallv_3
    procedure :: dmpi_grp_alltoallv_i8, zmpi_grp_alltoallv_i8, impi_grp_alltoallv_i8, lmpi_grp_alltoallv_i8
    generic :: alltoallv => dmpi_grp_alltoallv_i8, zmpi_grp_alltoallv_i8, impi_grp_alltoallv_i8, lmpi_grp_alltoallv_i8
    procedure :: dmpi_grp_alltoall, zmpi_grp_alltoall, impi_grp_alltoall, lmpi_grp_alltoall
    generic :: alltoall => dmpi_grp_alltoall, zmpi_grp_alltoall, impi_grp_alltoall, lmpi_grp_alltoall
    procedure :: dmpi_grp_allgatherv, zmpi_grp_allgatherv, impi_grp_allgatherv, lmpi_grp_allgatherv
    generic :: allgatherv => dmpi_grp_allgatherv, zmpi_grp_allgatherv, impi_grp_allgatherv, lmpi_grp_allgatherv
    procedure :: dmpi_grp_allgatherv_2, zmpi_grp_allgatherv_2, impi_grp_allgatherv_2, lmpi_grp_allgatherv_2
    generic :: allgatherv => dmpi_grp_allgatherv_2, zmpi_grp_allgatherv_2, impi_grp_allgatherv_2, lmpi_grp_allgatherv_2
    procedure :: dmpi_grp_allgatherv_3, zmpi_grp_allgatherv_3, impi_grp_allgatherv_3, lmpi_grp_allgatherv_3
    generic :: allgatherv => dmpi_grp_allgatherv_3, zmpi_grp_allgatherv_3, impi_grp_allgatherv_3, lmpi_grp_allgatherv_3
    procedure :: dmpi_grp_allgatherv_3_1, zmpi_grp_allgatherv_3_1, impi_grp_allgatherv_3_1, lmpi_grp_allgatherv_3_1
    generic :: allgatherv => dmpi_grp_allgatherv_3_1, zmpi_grp_allgatherv_3_1, impi_grp_allgatherv_3_1, lmpi_grp_allgatherv_3_1
    procedure :: dmpi_grp_allgatherv_i8, zmpi_grp_allgatherv_i8, impi_grp_allgatherv_i8, lmpi_grp_allgatherv_i8
    generic :: allgatherv => dmpi_grp_allgatherv_i8, zmpi_grp_allgatherv_i8, impi_grp_allgatherv_i8, lmpi_grp_allgatherv_i8
    procedure :: dmpi_grp_bcast, zmpi_grp_bcast, impi_grp_bcast, lmpi_grp_bcast
    generic :: bcast => dmpi_grp_bcast, zmpi_grp_bcast, impi_grp_bcast, lmpi_grp_bcast
    procedure :: dmpi_grp_bcast_0, zmpi_grp_bcast_0, impi_grp_bcast_0, lmpi_grp_bcast_0
    generic :: bcast => dmpi_grp_bcast_0, zmpi_grp_bcast_0, impi_grp_bcast_0, lmpi_grp_bcast_0
    procedure :: dmpi_grp_bcast_2, zmpi_grp_bcast_2, impi_grp_bcast_2, lmpi_grp_bcast_2
    generic :: bcast => dmpi_grp_bcast_2, zmpi_grp_bcast_2, impi_grp_bcast_2, lmpi_grp_bcast_2
    procedure :: dmpi_grp_bcast_3, zmpi_grp_bcast_3, impi_grp_bcast_3, lmpi_grp_bcast_3
    generic :: bcast => dmpi_grp_bcast_3, zmpi_grp_bcast_3, impi_grp_bcast_3, lmpi_grp_bcast_3
    procedure :: chmpi_grp_bcast_0, lompi_grp_bcast_0
    generic :: bcast => chmpi_grp_bcast_0, lompi_grp_bcast_0
    procedure :: dmpi_grp_bcast_0_l, zmpi_grp_bcast_0_l, impi_grp_bcast_0_l, lmpi_grp_bcast_0_l
    generic :: bcast => dmpi_grp_bcast_0_l, zmpi_grp_bcast_0_l, impi_grp_bcast_0_l, lmpi_grp_bcast_0_l
    procedure :: dmpi_grp_allreduce, zmpi_grp_allreduce, impi_grp_allreduce, lmpi_grp_allreduce
    generic :: allreduce => dmpi_grp_allreduce, zmpi_grp_allreduce, impi_grp_allreduce, lmpi_grp_allreduce
    procedure :: dmpi_grp_allreduce_2, zmpi_grp_allreduce_2, impi_grp_allreduce_2, lmpi_grp_allreduce_2
    generic :: allreduce => dmpi_grp_allreduce_2, zmpi_grp_allreduce_2, impi_grp_allreduce_2, lmpi_grp_allreduce_2
    procedure :: dmpi_grp_allreduce_3, zmpi_grp_allreduce_3, impi_grp_allreduce_3, lmpi_grp_allreduce_3
    generic :: allreduce => dmpi_grp_allreduce_3, zmpi_grp_allreduce_3, impi_grp_allreduce_3, lmpi_grp_allreduce_3
    procedure :: dmpi_grp_allreduce_0, zmpi_grp_allreduce_0, impi_grp_allreduce_0, lmpi_grp_allreduce_0
    generic :: allreduce => dmpi_grp_allreduce_0, zmpi_grp_allreduce_0, impi_grp_allreduce_0, lmpi_grp_allreduce_0
    procedure :: lompi_grp_allreduce_0
    generic :: allreduce => lompi_grp_allreduce_0
    procedure :: dmpi_grp_allreduce_inplace_0, zmpi_grp_allreduce_inplace_0
    procedure :: impi_grp_allreduce_inplace_0, lmpi_grp_allreduce_inplace_0
    procedure :: lompi_grp_allreduce_inplace_0
    generic :: allreduce_inplace => dmpi_grp_allreduce_inplace_0, zmpi_grp_allreduce_inplace_0
    generic :: allreduce_inplace => impi_grp_allreduce_inplace_0, lmpi_grp_allreduce_inplace_0
    generic :: allreduce_inplace => lompi_grp_allreduce_inplace_0
    procedure :: dmpi_grp_allgather, zmpi_grp_allgather, impi_grp_allgather, lmpi_grp_allgather
    generic :: allgather => dmpi_grp_allgather, zmpi_grp_allgather, impi_grp_allgather, lmpi_grp_allgather
    procedure :: dmpi_grp_allgather_0, zmpi_grp_allgather_0, impi_grp_allgather_0, lmpi_grp_allgather_0
    generic :: allgather => dmpi_grp_allgather_0, zmpi_grp_allgather_0, impi_grp_allgather_0, lmpi_grp_allgather_0
    procedure :: dmpi_grp_recv, zmpi_grp_recv, impi_grp_recv, lmpi_grp_recv
    generic :: recv => dmpi_grp_recv, zmpi_grp_recv, impi_grp_recv, lmpi_grp_recv
    procedure :: dmpi_grp_recv_0, zmpi_grp_recv_0, impi_grp_recv_0, lmpi_grp_recv_0
    generic :: recv => dmpi_grp_recv_0, zmpi_grp_recv_0, impi_grp_recv_0, lmpi_grp_recv_0
    procedure :: dmpi_grp_recv_2, zmpi_grp_recv_2, impi_grp_recv_2, lmpi_grp_recv_2
    generic :: recv => dmpi_grp_recv_2, zmpi_grp_recv_2, impi_grp_recv_2, lmpi_grp_recv_2
    procedure :: dmpi_grp_recv_3, zmpi_grp_recv_3, impi_grp_recv_3, lmpi_grp_recv_3
    generic :: recv => dmpi_grp_recv_3, zmpi_grp_recv_3, impi_grp_recv_3, lmpi_grp_recv_3
    procedure :: lompi_grp_recv_0
    generic :: recv => lompi_grp_recv_0
    procedure :: dmpi_grp_send, zmpi_grp_send, impi_grp_send, lmpi_grp_send
    generic :: send => dmpi_grp_send, zmpi_grp_send, impi_grp_send, lmpi_grp_send
    procedure :: dmpi_grp_send_0, zmpi_grp_send_0, impi_grp_send_0, lmpi_grp_send_0
    generic :: send => dmpi_grp_send_0, zmpi_grp_send_0, impi_grp_send_0, lmpi_grp_send_0
    procedure :: dmpi_grp_send_2, zmpi_grp_send_2, impi_grp_send_2, lmpi_grp_send_2
    generic :: send => dmpi_grp_send_2, zmpi_grp_send_2, impi_grp_send_2, lmpi_grp_send_2
    procedure :: dmpi_grp_send_3, zmpi_grp_send_3, impi_grp_send_3, lmpi_grp_send_3
    generic :: send => dmpi_grp_send_3, zmpi_grp_send_3, impi_grp_send_3, lmpi_grp_send_3
    procedure :: lompi_grp_send_0
    generic :: send => lompi_grp_send_0
    procedure :: dmpi_grp_irecv, zmpi_grp_irecv, impi_grp_irecv, lmpi_grp_irecv
    generic :: irecv => dmpi_grp_irecv, zmpi_grp_irecv, impi_grp_irecv, lmpi_grp_irecv
    procedure :: dmpi_grp_irecv_0, zmpi_grp_irecv_0, impi_grp_irecv_0, lmpi_grp_irecv_0
    generic :: irecv => dmpi_grp_irecv_0, zmpi_grp_irecv_0, impi_grp_irecv_0, lmpi_grp_irecv_0
    procedure :: dmpi_grp_irecv_2, zmpi_grp_irecv_2, impi_grp_irecv_2, lmpi_grp_irecv_2
    generic :: irecv => dmpi_grp_irecv_2, zmpi_grp_irecv_2, impi_grp_irecv_2, lmpi_grp_irecv_2
    procedure :: dmpi_grp_irecv_3, zmpi_grp_irecv_3, impi_grp_irecv_3, lmpi_grp_irecv_3
    generic :: irecv => dmpi_grp_irecv_3, zmpi_grp_irecv_3, impi_grp_irecv_3, lmpi_grp_irecv_3
    procedure :: dmpi_grp_irecv_0_i8, zmpi_grp_irecv_0_i8, impi_grp_irecv_0_i8, lmpi_grp_irecv_0_i8
    generic :: irecv => dmpi_grp_irecv_0_i8, zmpi_grp_irecv_0_i8, impi_grp_irecv_0_i8, lmpi_grp_irecv_0_i8
    procedure :: dmpi_grp_isend, zmpi_grp_isend, impi_grp_isend, lmpi_grp_isend
    generic :: isend => dmpi_grp_isend, zmpi_grp_isend, impi_grp_isend, lmpi_grp_isend
    procedure :: dmpi_grp_isend_0, zmpi_grp_isend_0, impi_grp_isend_0, lmpi_grp_isend_0
    generic :: isend => dmpi_grp_isend_0, zmpi_grp_isend_0, impi_grp_isend_0, lmpi_grp_isend_0
    procedure :: dmpi_grp_isend_2, zmpi_grp_isend_2, impi_grp_isend_2, lmpi_grp_isend_2
    generic :: isend => dmpi_grp_isend_2, zmpi_grp_isend_2, impi_grp_isend_2, lmpi_grp_isend_2
    procedure :: dmpi_grp_isend_3, zmpi_grp_isend_3, impi_grp_isend_3, lmpi_grp_isend_3
    generic :: isend => dmpi_grp_isend_3, zmpi_grp_isend_3, impi_grp_isend_3, lmpi_grp_isend_3
    procedure :: dmpi_grp_isend_0_i8, zmpi_grp_isend_0_i8, impi_grp_isend_0_i8, lmpi_grp_isend_0_i8
    generic :: isend => dmpi_grp_isend_0_i8, zmpi_grp_isend_0_i8, impi_grp_isend_0_i8, lmpi_grp_isend_0_i8
    procedure :: mpi_grp_wait, mpi_grp_waitall
    generic :: wait => mpi_grp_wait, mpi_grp_waitall
  end type mpi_grp_t

  type(mpi_grp_t), public :: mpi_world

  !> used to store return values of mpi calls
  integer, public :: mpi_err

#if !defined(HAVE_MPI)
  ! constants needed for compilation without MPI
  ! this makes it possible to get rid of lots of #ifdef clauses
  integer, parameter, public :: &
    MPI_INTEGER = 0,            &
    MPI_INTEGER8 = 0,           &
    MPI_DOUBLE_PRECISION = 0,   &
    MPI_DOUBLE_COMPLEX = 0,     &
    MPI_2FLOAT = 0,             &
    MPI_CHARACTER = 0,          &
    MPI_LOGICAL = 0,            &
    MPI_SUM = 0,                &
    MPI_MINLOC = 0,             &
    MPI_LOR = 0,                &
    MPI_LAND = 0,               &
    MPI_MAX = 0,                &
    MPI_MIN = 0,                &
    MPI_IN_PLACE = 0
#endif

contains

  ! ---------------------------------------------------------
  subroutine mpi_mod_init(is_serial)
    logical, intent(in) :: is_serial

#if defined(HAVE_MPI)
#ifdef HAVE_OPENMP
    integer :: provided
#endif
#ifdef HAVE_SCALAPACK
    integer :: iam, nprocs
    integer :: blacs_default_system_context !< for blacs/openmpi bug workaround
#endif
#endif

    if (is_serial) then
      call mpi_grp_init(mpi_world, -1)
      return
    end if

#if defined(HAVE_MPI)
    ! initialize MPI
#if defined(HAVE_OPENMP)
    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, mpi_err)
#else
    call MPI_Init(mpi_err)
#endif
    call mpi_grp_init(mpi_world, MPI_COMM_WORLD)
    call MPI_Barrier(mpi_world%comm, mpi_err)
#else
    call mpi_grp_init(mpi_world, -1)
#endif

#ifdef HAVE_SCALAPACK

    ! Initialize Blacs to be able to use ScaLAPACK
    ! Determine my process number and the number of processes in machine
    call blacs_pinfo(iam, nprocs)

    ! If machine needs additional set up, do it now
    if (nprocs < 1) then
      call blacs_setup(iam, mpi_world%size)
    end if

    !> ensure there was at least one call to blacs_gridinit() or blacs_gridmap()
    !> without it, blacs_exit() triggers an error
    !>
    !> *** An error occurred in MPI_Type_free
    !> *** MPI_ERR_TYPE: invalid datatype
    !>
    !> in openmpi
    call blacs_get(-1,0, blacs_default_system_context)
    call blacs_gridinit(blacs_default_system_context, 'R', 1, 1)

#endif
  end subroutine mpi_mod_init


  ! ---------------------------------------------------------
  subroutine mpi_mod_end()

#ifdef HAVE_SCALAPACK
    if (mpi_world%comm /= -1) call blacs_exit(1)
#endif

#if defined(HAVE_MPI)
    ! end MPI, if we started it
    if (mpi_world%comm /= -1) call MPI_Finalize(mpi_err)
#endif

  end subroutine mpi_mod_end


  ! ---------------------------------------------------------
  subroutine mpi_grp_init(grp, comm)
    type(mpi_grp_t), intent(out)  :: grp   !< information about this MPI group
    integer,         intent(in)   :: comm  !< the communicator that defined the group

    grp%comm = comm
#if defined(HAVE_MPI)
    if (grp%comm == MPI_COMM_NULL) grp%comm = -1
#endif

    if (grp%comm == -1) then
      grp%rank = 0
      grp%size = 1
#if defined(HAVE_MPI)
    else
      call MPI_Comm_rank(grp%comm, grp%rank, mpi_err)
      call MPI_Comm_size(grp%comm, grp%size, mpi_err)
#endif
    end if

  end subroutine mpi_grp_init

  ! ---------------------------------------------------------
  subroutine mpi_grp_copy(mpi_grp_out, mpi_grp_in)
    type(mpi_grp_t), intent(out) :: mpi_grp_out
    type(mpi_grp_t), intent(in)  :: mpi_grp_in

    mpi_grp_out%comm = mpi_grp_in%comm
    mpi_grp_out%size = mpi_grp_in%size
    mpi_grp_out%rank = mpi_grp_in%rank
  end subroutine mpi_grp_copy

  ! ---------------------------------------------------------
  subroutine mpi_grp_duplicate(mpi_grp_out, mpi_grp_in)
    type(mpi_grp_t), intent(out) :: mpi_grp_out
    type(mpi_grp_t), intent(in)  :: mpi_grp_in

#if defined(HAVE_MPI)
    call MPI_Comm_dup(mpi_grp_in%comm, mpi_grp_out%comm, mpi_err)
    call MPI_Comm_rank(mpi_grp_out%comm, mpi_grp_out%rank, mpi_err)
    call MPI_Comm_size(mpi_grp_out%comm, mpi_grp_out%size, mpi_err)
#else
    call mpi_grp_copy(mpi_grp_out, mpi_grp_in)
#endif
  end subroutine mpi_grp_duplicate

  ! ---------------------------------------------------------
  logical function mpi_grp_is_root(grp)
    type(mpi_grp_t), intent(in) :: grp

    mpi_grp_is_root = (grp%rank == 0)
  end function mpi_grp_is_root

  ! ---------------------------------------------------------
  subroutine mpi_grp_barrier(mpi_grp)
    class(mpi_grp_t), intent(in) :: mpi_grp

    if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
    call mpi_debug_in(mpi_grp%comm, C_MPI_BARRIER)
    call MPI_Barrier(mpi_grp%comm, mpi_err)
    call mpi_debug_out(mpi_grp%comm, C_MPI_BARRIER)
#endif
  end subroutine mpi_grp_barrier

  ! ---------------------------------------------------------
  subroutine chmpi_grp_bcast_0(mpi_grp, buf, cnt, sendtype, root)
    class(mpi_grp_t), intent(in)    :: mpi_grp
    character(len=*), intent(inout) :: buf
    integer,          intent(in)    :: cnt, sendtype, root

#if defined(HAVE_MPI)
    call mpi_debug_in(mpi_grp%comm, C_MPI_BCAST)
    if (mpi_grp%comm /= -1) then
      call MPI_Bcast(buf, cnt, sendtype, root, mpi_grp%comm, mpi_err)
    end if
    call mpi_debug_out(mpi_grp%comm, C_MPI_BCAST)
#endif
  end subroutine chmpi_grp_bcast_0

  ! ---------------------------------------------------------
  subroutine lompi_grp_bcast_0(mpi_grp, buf, cnt, sendtype, root)
    class(mpi_grp_t), intent(in)    :: mpi_grp
    logical,          intent(inout) :: buf
    integer,          intent(in)    :: cnt, sendtype, root

#if defined(HAVE_MPI)
    call mpi_debug_in(mpi_grp%comm, C_MPI_BCAST)
    if (mpi_grp%comm /= -1) then
      call MPI_Bcast(buf, cnt, sendtype, root, mpi_grp%comm, mpi_err)
    end if
    call mpi_debug_out(mpi_grp%comm, C_MPI_BCAST)
#endif
  end subroutine lompi_grp_bcast_0

  ! ---------------------------------------------------------
  ! copy routine for serial case
  subroutine lompi_grp_copy_0(sendbuf, recvbuf, count)
    use iso_c_binding
    logical, target,  intent(in)  :: sendbuf
    logical, target,  intent(out) :: recvbuf
    integer,          intent(in)  :: count
    integer :: ii
    logical, pointer :: send(:), recv(:)

    call c_f_pointer(c_loc(sendbuf), send, [count])
    call c_f_pointer(c_loc(recvbuf), recv, [count])
    do ii = 1, count
      recv(ii) = send(ii)
    end do
  end subroutine lompi_grp_copy_0

  ! ---------------------------------------------------------
  subroutine lompi_grp_allreduce_0(mpi_grp, sendbuf, recvbuf, count, datatype, op)
    class(mpi_grp_t), intent(in)  :: mpi_grp
    logical,          intent(in)  :: sendbuf
    logical,          intent(out) :: recvbuf
    integer,          intent(in)  :: count, datatype, op

#if defined(HAVE_MPI)
    call mpi_debug_in(mpi_grp%comm, C_MPI_ALLREDUCE)
    if (mpi_grp%comm /= -1) then
      call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
        mpi_grp%comm, mpi_err)
    else
      call lompi_grp_copy_0(sendbuf, recvbuf, count)
    end if
    call mpi_debug_out(mpi_grp%comm, C_MPI_ALLREDUCE)
#else
    call lompi_grp_copy_0(sendbuf, recvbuf, count)
#endif
  end subroutine lompi_grp_allreduce_0

  ! ---------------------------------------------------------
  subroutine lompi_grp_allreduce_inplace_0(mpi_grp, recvbuf, count, datatype, op)
    class(mpi_grp_t), intent(in)    :: mpi_grp
    logical,          intent(inout) :: recvbuf
    integer,          intent(in)    :: count, datatype, op

#if defined(HAVE_MPI)
    call mpi_debug_in(mpi_grp%comm, C_MPI_ALLREDUCE)
    if (mpi_grp%comm /= -1) then
      call MPI_Allreduce(MPI_IN_PLACE, recvbuf, count, datatype, op, &
        mpi_grp%comm, mpi_err)
    end if
    call mpi_debug_out(mpi_grp%comm, C_MPI_ALLREDUCE)
#endif
  end subroutine lompi_grp_allreduce_inplace_0

  ! ---------------------------------------------------------
  subroutine lompi_grp_recv_0(mpi_grp, recvbuf, recvcount, recvtype, source, tag)
    class(mpi_grp_t),  intent(in)  :: mpi_grp
    logical,           intent(out) :: recvbuf
    integer,           intent(in)  :: recvcount, recvtype
    integer,           intent(in)  :: source
    integer, optional, intent(in)  :: tag

    integer :: tag_

    tag_ = 0
    if (present(tag)) tag_ = tag
    if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
    call MPI_Recv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, MPI_STATUS_IGNORE, mpi_err)
#endif
  end subroutine lompi_grp_recv_0

  ! ---------------------------------------------------------
  subroutine lompi_grp_send_0(mpi_grp, sendbuf, sendcount, sendtype, dest, tag)
    class(mpi_grp_t),  intent(in)  :: mpi_grp
    logical,           intent(out) :: sendbuf
    integer,           intent(in)  :: sendcount, sendtype
    integer,           intent(in)  :: dest
    integer, optional, intent(in)  :: tag

    integer :: tag_

    tag_ = 0
    if (present(tag)) tag_ = tag
    if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
    call MPI_Send(sendbuf, sendcount, sendtype, dest, tag_, mpi_grp%comm, mpi_err)
#endif
  end subroutine lompi_grp_send_0

  ! ---------------------------------------------------------
  subroutine mpi_grp_wait(mpi_grp, request)
    class(mpi_grp_t),  intent(in)    :: mpi_grp
    integer,           intent(inout) :: request

    if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
    call MPI_Wait(request, MPI_STATUS_IGNORE, mpi_err)
#endif
  end subroutine mpi_grp_wait

  ! ---------------------------------------------------------
  subroutine mpi_grp_waitall(mpi_grp, count, requests)
    class(mpi_grp_t),  intent(in)    :: mpi_grp
    integer,           intent(in)    :: count
    integer,           intent(inout) :: requests(:)

    if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
    call MPI_Waitall(count, requests, MPI_STATUSES_IGNORE, mpi_err)
#endif
  end subroutine mpi_grp_waitall

#include "undef.F90"
#include "real.F90"
#include "mpi_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mpi_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "mpi_inc.F90"

#include "undef.F90"
#include "integer8.F90"
#include "mpi_inc.F90"

end module mpi_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
