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

module blacs_proc_grid_oct_m
  use debug_oct_m
  use global_oct_m
  use blacs_oct_m
  use mpi_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                      &
    blacs_proc_grid_t,           &
    blacs_proc_grid_init,        &
    blacs_proc_grid_end,         &
    blacs_proc_grid_copy,        &
    blacs_proc_grid_null

  type blacs_proc_grid_t
    ! Components are public by default
    integer          :: context = -1  !< The blacs context, -1 is object is null.
    integer          :: nprocs        !< Number of processors.
    integer          :: nprow         !< Number of processors per row.
    integer          :: npcol         !< Number of processors per column.
    integer, private :: iam           !< Process indentifier.
    integer          :: myrow         !< The row of the processor in the processor grid.
    integer          :: mycol         !< The column of the processor in the processor grid.
    integer, allocatable :: usermap(:, :) !< The index of each processor in the grid.
  end type blacs_proc_grid_t

contains

  ! -----------------------------------------------------------------------
  !> Initializes a blacs context from an MPI communicator with
  !! topological information.
  !!
  !! \Warning: For the moment this function only works if mpi_grp holds
  !! all the nodes of mpi_world.
  subroutine blacs_proc_grid_init(this, mpi_grp, procdim)
    type(blacs_proc_grid_t),           intent(inout) :: this
    type(mpi_grp_t),                   intent(in)    :: mpi_grp
    integer,                 optional, intent(in)    :: procdim(:)

#ifdef HAVE_SCALAPACK

    integer, parameter :: maxdims = 2
    integer :: dims(1:2), topo, coords(1:2), ix, iy, id, xy(2)
    logical :: periods(1:2)
    integer :: comm
    logical :: reorder
    integer, allocatable :: procmap(:)

    PUSH_SUB(blacs_proc_grid_init)

    call MPI_Topo_test(mpi_grp%comm, topo, mpi_err)

    if (topo /= MPI_CART .or. present(procdim)) then
      ! We create a new communicator with Cartesian topology
      if (present(procdim)) then
        dims(1) = procdim(1)
        dims(2) = procdim(2)
      else
        dims(1) = mpi_grp%size
        dims(2) = 1
      end if
      periods = .false.
      reorder = .false.
      call MPI_Cart_create(mpi_grp%comm, 2, dims, periods, reorder, comm, mpi_err)
    else
      comm = mpi_grp%comm
    end if

    call blacs_pinfo(this%iam, this%nprocs)

    ! The process ID from ScaLAPACK is not always the
    ! same as MPI, so we need to construct a map.
    SAFE_ALLOCATE(procmap(0:mpi_grp%size - 1))
    call MPI_Allgather(this%iam, 1, MPI_INTEGER, procmap(0), 1, MPI_INTEGER, comm, mpi_err)

    ASSERT(this%iam == procmap(mpi_grp%rank))

    dims = 1
    coords = 0

    call MPI_Cart_get(comm, maxdims, dims, periods, coords, mpi_err)

    SAFE_ALLOCATE(this%usermap(1:dims(1), 1:dims(2)))

    do ix = 1, dims(1)
      xy(1) = ix - 1
      do iy = 1, dims(2)
        xy(2) = iy - 1
        call MPI_Cart_rank(comm, xy, id, mpi_err)
        this%usermap(ix, iy) = procmap(id)
      end do
    end do

    ! get the default system context
    call blacs_get(-1, what = 0, val = this%context)

    ! now get the context associated with the map
    call blacs_gridmap(this%context, this%usermap(1, 1), dims(1), dims(1), dims(2))

    ! and fill the rest of the structure
    call blacs_gridinfo(this%context, this%nprow, this%npcol, this%myrow, this%mycol)

    !check that Blacs and MPI are consistent
    ASSERT(this%nprow == dims(1))
    ASSERT(this%npcol == dims(2))
    ASSERT(this%myrow == coords(1))
    ASSERT(this%mycol == coords(2))

    if (topo /= MPI_CART) then
      call MPI_Comm_free(comm, mpi_err)
    end if

    SAFE_DEALLOCATE_A(procmap)

    POP_SUB(blacs_proc_grid_init)
#endif
  end subroutine blacs_proc_grid_init

  ! ----------------------------------------------------

  subroutine blacs_proc_grid_end(this)
    type(blacs_proc_grid_t), intent(inout) :: this

    PUSH_SUB(blacs_proc_grid_end)

    if (this%context /= -1) then
#ifdef HAVE_SCALAPACK
      call blacs_gridexit(this%context)
#endif
      SAFE_DEALLOCATE_A(this%usermap)
    end if

    this%context = -1

    POP_SUB(blacs_proc_grid_end)
  end subroutine blacs_proc_grid_end

  ! ----------------------------------------------------

  subroutine blacs_proc_grid_copy(cin, cout)
    type(blacs_proc_grid_t), intent(in)    :: cin
    type(blacs_proc_grid_t), intent(inout) :: cout

    PUSH_SUB(blacs_proc_grid_copy)

    call blacs_proc_grid_end(cout)

    cout%context = cin%context

#ifdef HAVE_SCALAPACK
    cout%nprocs  = cin%nprocs
    cout%nprow   = cin%nprow
    cout%npcol   = cin%npcol
    cout%iam     = cin%iam
    cout%myrow   = cin%myrow
    cout%mycol   = cin%mycol

    if (cout%context /= -1) then
      ! we have to create a new context
      call blacs_get(-1, what = 0, val = cout%context)
      SAFE_ALLOCATE_SOURCE_A(cout%usermap, cin%usermap)
      call blacs_gridmap(cout%context, cout%usermap(1, 1), cout%nprow, cout%nprow, cout%npcol)
    end if

#endif
    POP_SUB(blacs_proc_grid_copy)
  end subroutine blacs_proc_grid_copy

  ! ----------------------------------------------------

  logical pure function blacs_proc_grid_null(this)
    type(blacs_proc_grid_t), intent(in) :: this

    blacs_proc_grid_null = this%context == -1
  end function blacs_proc_grid_null

end module blacs_proc_grid_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
