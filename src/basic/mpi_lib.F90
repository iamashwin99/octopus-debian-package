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

!> This module contains some common usage patterns of MPI routines.
module mpi_lib_oct_m
  use debug_oct_m
  use global_oct_m
  use iso_c_binding
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use profiling_oct_m

  implicit none

  private

#if !defined(HAVE_MPI)
  integer, public :: mpi_lib_dummy !< this avoids compilers complaining about empty module
#else
  public ::                            &
    create_intranode_communicator,     &
    lmpi_create_shared_memory_window,  &
    lmpi_destroy_shared_memory_window, &
    lmpi_sync_shared_memory_window,    &
    lmpi_gen_allgatherv,               &
    lmpi_translate_rank

  interface lmpi_create_shared_memory_window
    module procedure dlmpi_create_shared_memory_window, zlmpi_create_shared_memory_window
    module procedure ilmpi_create_shared_memory_window, llmpi_create_shared_memory_window
  end interface lmpi_create_shared_memory_window

  interface lmpi_gen_allgatherv
    module procedure dlmpi_gen_allgatherv, zlmpi_gen_allgatherv
    module procedure ilmpi_gen_allgatherv, llmpi_gen_allgatherv
  end interface lmpi_gen_allgatherv

  type(profile_t), save :: prof_allgatherv
#endif

contains
#if defined(HAVE_MPI)
  ! ---------------------------------------------------------
  !> Returns the rank number of the node rank in from_comm for the
  !! to_comm communicator.
  integer function lmpi_translate_rank(from_comm, to_comm, rank)
    integer, intent(in) :: from_comm
    integer, intent(in) :: to_comm
    integer, intent(in) :: rank

    integer :: from_group, to_group, from_rank(1), to_rank(1)

    PUSH_SUB(lmpi_translate_rank)

    call MPI_Comm_group(from_comm, from_group, mpi_err)
    call MPI_Comm_group(to_comm, to_group, mpi_err)

    from_rank(1) = rank
    call MPI_Group_translate_ranks(from_group, 1, from_rank, to_group, to_rank, mpi_err)

    lmpi_translate_rank = to_rank(1)

    POP_SUB(lmpi_translate_rank)
  end function lmpi_translate_rank

  ! destroy a shared memory segment as an MPI window
  subroutine lmpi_destroy_shared_memory_window(window)
    integer, intent(inout) :: window

    PUSH_SUB(lmpi_destroy_shared_memory_window)

    ! end access epoch
    call MPI_Win_unlock_all(window, mpi_err)
    ! destroy window
    call MPI_Win_free(window, mpi_err)

    POP_SUB(lmpi_destroy_shared_memory_window)
  end subroutine lmpi_destroy_shared_memory_window

  ! synchronize a shared memory segment
  subroutine lmpi_sync_shared_memory_window(window, intranode_grp)
    integer,         intent(in) :: window
    type(mpi_grp_t), intent(in) :: intranode_grp

    PUSH_SUB(lmpi_sync_shared_memory_window)

    call MPI_Win_sync(window, mpi_err)
    call intranode_grp%barrier()

    POP_SUB(lmpi_sync_shared_memory_window)
  end subroutine lmpi_sync_shared_memory_window

  subroutine create_intranode_communicator(base_grp, intranode_grp, internode_grp)
    type(mpi_grp_t), intent(in)  :: base_grp
    type(mpi_grp_t), intent(out) :: intranode_grp
    type(mpi_grp_t), intent(out) :: internode_grp
    integer :: comm

    PUSH_SUB(create_intranode_communicator)

    ! create communicator for intranode communication
    ! this is used for MPI-3 shared memory windows
    call MPI_Comm_split_type(base_grp%comm, MPI_COMM_TYPE_SHARED, base_grp%rank, MPI_INFO_NULL, comm, mpi_err)
    call mpi_grp_init(intranode_grp, comm)

    ! create communicator for internode communication
    ! we have one communicator for each local rank
    call MPI_Comm_split(base_grp%comm, intranode_grp%rank, base_grp%rank, comm, mpi_err)
    call mpi_grp_init(internode_grp, comm)

    POP_SUB(create_intranode_communicator)
  end subroutine create_intranode_communicator

#include "undef.F90"
#include "real.F90"
#include "mpi_lib_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mpi_lib_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "mpi_lib_inc.F90"

#include "undef.F90"
#include "integer8.F90"
#include "mpi_lib_inc.F90"
#endif
end module mpi_lib_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
