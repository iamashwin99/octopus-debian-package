!! Copyright (C) 2005-2006 Florian Lorenzen, Heiko Appel
!! Copyright (C) 2021 Sebastian Ohlmann
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

!> Generally:
!! Xpar_vec_gather and Xpar_vec_scatter only consider inner points.
!! Xpar_vec_scatter_bndry takes care of boundary points (there is
!! no Xpar_vec_gather_bndry as they are only written and not read).
!! Xpar_vec_scatter_all is Xpar_vec_scatter followd by Xpar_vec_scatter_bndry.

!! ---------------------------------------------------------
!! Scatters a vector v to all nodes in pv with respect to
!! to point -> node mapping in pv.
!! v_local has at least to be of size pv%np_local.
subroutine X(par_vec_scatter)(pv, root, v_local, v)
  type(par_vec_t), intent(in)  :: pv
  integer,         intent(in)  :: root
  R_TYPE,          intent(out) :: v_local(:)
  R_TYPE,          intent(in)  :: v(:)

  integer(i8)              :: ii        !< Counter.
  integer(i8), allocatable :: displs(:) !< Displacements for scatter.
  integer(i8), allocatable :: local_vec(:) !< mapping of points
  R_TYPE,  allocatable :: v_tmp(:)  !< Send buffer.
  type(profile_t), save :: prof_scatter

  PUSH_SUB(X(par_vec_scatter))
  call profiling_in(prof_scatter, TOSTRING(X(VEC_SCATTER)))

  ! Skip the MPI call if domain parallelization is not used.
  if (pv%npart < 2) then
    v_local(1:pv%np_global) = v(1:pv%np_global)
    POP_SUB(X(par_vec_scatter))
    return
  end if

  ! Unfortunately, pv%xlocal_vec ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:pv%npart))
  displs = pv%xlocal_vec - 1

  if (root == pv%rank) then
    SAFE_ALLOCATE(local_vec(1:pv%np_global))
  else
    SAFE_ALLOCATE(local_vec(0))
  end if
  call gather_local_vec(pv, root, local_vec)

  SAFE_ALLOCATE(v_tmp(1))
  if (pv%rank == root) then
    ! Fill send buffer.
    SAFE_DEALLOCATE_A(v_tmp)
    SAFE_ALLOCATE(v_tmp(1:pv%np_global))

    ! Rearrange copy of v. All points of node r are in
    ! v_tmp(xlocal_vec(r):xlocal_vec(r)+np_local_vec(r)-1).
    do ii = 1, pv%np_global
      v_tmp(ii) = v(local_vec(ii))
    end do
  end if
  SAFE_DEALLOCATE_A(local_vec)

  ! Careful: MPI rank numbers range from 0 to mpiv%numprocs-1
  ! But partition numbers from 1 to pv%npart with usually
  ! pv%npart = mpiv%numprocs.
  call pv%mpi_grp%scatterv(v_tmp, pv%np_local_vec, displs, R_MPITYPE, &
    v_local, pv%np_local, R_MPITYPE, root)

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call profiling_out(prof_scatter)

  POP_SUB(X(par_vec_scatter))
end subroutine X(par_vec_scatter)


! ---------------------------------------------------------
!> Reverse operation of Xpar_vec_scatter.
!! All v_locals from the nodes are packed together
!! into v on node root in correct order.
subroutine X(par_vec_gather)(pv, root, v_local, v)
  type(par_vec_t),      intent(in)  :: pv
  integer,              intent(in)  :: root
  R_TYPE,               intent(in)  :: v_local(:)
  R_TYPE,     optional, intent(out) :: v(:) !< in order to prevent unassociated pointer errors,
  !                                         !< this is optional, so that mpi ranks not expecting an output
  !                                         !< do not have to pass a null pointer.

  integer(i8)              :: ii        !< Counter.
  integer(i8), allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Receive buffer.
  integer(i8), allocatable :: local_vec(:) !< mapping of points

  PUSH_SUB(X(par_vec_gather))

  ! Skip the MPI call if domain parallelization is not used.
  if (pv%npart < 2) then
    v(1:pv%np_global) = v_local(1:pv%np_global)
    POP_SUB(X(par_vec_gather))
    return
  end if

  ! Unfortunately, pv%xlocal_vec ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:pv%npart))
  displs = pv%xlocal_vec - 1

  if (pv%rank == root) then
    SAFE_ALLOCATE(v_tmp(1:pv%np_global))
  else
    SAFE_ALLOCATE(v_tmp(1:1))
  end if

  call pv%mpi_grp%gatherv(v_local, pv%np_local, R_MPITYPE, v_tmp, &
    pv%np_local_vec, displs, R_MPITYPE, root)

  if (root == pv%rank) then
    SAFE_ALLOCATE(local_vec(1:pv%np_global))
  else
    SAFE_ALLOCATE(local_vec(1:1))
  end if
  call gather_local_vec(pv, root, local_vec)

  ! Copy values from v_tmp to their original position in v.
  if (pv%rank == root) then
    do ii = 1, pv%np_global
      v(local_vec(ii)) = v_tmp(ii)
    end do
  end if
  SAFE_DEALLOCATE_A(local_vec)

  SAFE_DEALLOCATE_A(local_vec)
  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  POP_SUB(X(par_vec_gather))
end subroutine X(par_vec_gather)

! ---------------------------------------------------------
!> Like Xpar_vec_gather but the result is gathered
!! on all nodes, i. e. v has to be a properly
!! allocated array on all nodes.
subroutine X(par_vec_allgather)(pv, v, v_local)
  type(par_vec_t), intent(in)  :: pv
  R_TYPE,          intent(out) :: v(:)
  R_TYPE,          intent(in)  :: v_local(:)

  type(profile_t), save :: prof_allgather
  integer, parameter :: root = 0

  PUSH_SUB(X(par_vec_allgather))
  call profiling_in(prof_allgather, TOSTRING(X(VEC_ALLGATHER)))

  call par_vec_gather(pv, root, v_local, v)
  if (pv%npart > 1) then
    call pv%mpi_grp%bcast(v(1), pv%np_global, R_MPITYPE, root)
  end if

  call profiling_out(prof_allgather)

  POP_SUB(X(par_vec_allgather))
end subroutine X(par_vec_allgather)

!--------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
