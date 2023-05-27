!! Copyright (C) 2005-2020 Florian Lorenzen, Heiko Appel, Martin Lueders
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

! ---------------------------------------------------------
!> Updates ghost points of every node. A vector suitable
!! for non-local operations contains local values and
!! ghost point values.
!! Length of v_local must be
!! pv%np_local+pv%np_ghost
subroutine X(par_vec_ghost_update)(pv, v_local)
  type(par_vec_t), intent(in)    :: pv
  R_TYPE,          intent(inout) :: v_local(:)

  R_TYPE,  allocatable :: ghost_send(:)
  integer :: ip
  type(profile_t), save :: prof_update

  call profiling_in(prof_update, TOSTRING(X(GHOST_UPDATE)))

  PUSH_SUB(X(par_vec_ghost_update))

  SAFE_ALLOCATE(ghost_send(1:pv%ghost_scount))

  ! pack data for sending
  do ip = 1, pv%ghost_scount
    ghost_send(ip) = v_local(pv%ghost_sendmap(ip))
  end do

  call pv%mpi_grp%alltoallv(ghost_send, pv%ghost_scounts, pv%ghost_sdispls, R_MPITYPE, &
    v_local(pv%np_local+1:), pv%ghost_rcounts, pv%ghost_rdispls, R_MPITYPE)

  SAFE_DEALLOCATE_A(ghost_send)

  POP_SUB(X(par_vec_ghost_update))

  call profiling_out(prof_update)
end subroutine X(par_vec_ghost_update)

! ---------------------------------------------------------

subroutine X(ghost_update_batch_start)(pv, v_local, handle)
  type(par_vec_t), target,      intent(in)    :: pv
  class(batch_t),  target,      intent(inout) :: v_local
  type(par_vec_handle_batch_t), intent(out)   :: handle

  integer :: ipart, pos, ii, nn, ip, ipart2
  integer :: offset, dim2, dim3, localsize
  type(profile_t), save :: prof_start, prof_irecv, prof_isend

  call profiling_in(prof_start, TOSTRING(X(GHOST_UPDATE_START)))
  PUSH_SUB(X(ghost_update_batch_start))

  ASSERT(v_local%nst_linear > 0)

  handle%nnb = 0
  handle%v_local => v_local
  handle%pv => pv

  SAFE_ALLOCATE(handle%requests(1:2*pv%npart*v_local%nst_linear))

  call profiling_in(prof_irecv, TOSTRING(X(GHOST_UPDATE_IRECV)))

  ! first post the receptions
  ! the communication scheme is in principle a sparse alltoallv:
  ! we use a ring scheme to post the receives and the sends which has
  ! the advantage that matching messages are posted at the same time,
  ! facilitating the matching of those messages
  select case (v_local%status())
  case (BATCH_DEVICE_PACKED)
    if (.not. accel%cuda_mpi) then
      SAFE_ALLOCATE(handle%X(recv_buffer)(1:v_local%pack_size(1)*pv%np_ghost))
      offset = 0
    else
      ! get device pointer for CUDA-aware MPI
      call accel_get_device_pointer(handle%X(recv_buffer), handle%v_local%ff_device, &
        [product(v_local%pack_size)])
      offset = pv%np_local*v_local%pack_size(1)
    end if

    ! ring scheme: count upwards from local rank for receiving
    do ipart2 = pv%partno, pv%partno + pv%npart - 1
      ipart = ipart2
      if (ipart > pv%npart) ipart = ipart - pv%npart
      if (pv%ghost_rcounts(ipart) == 0) cycle

      handle%nnb = handle%nnb + 1
      pos = 1 + pv%ghost_rdispls(ipart)*v_local%pack_size(1) + offset
      call pv%mpi_grp%irecv(handle%X(recv_buffer)(pos), pv%ghost_rcounts(ipart)*v_local%pack_size(1), R_MPITYPE, &
        ipart - 1, handle%requests(handle%nnb))
    end do

  case (BATCH_PACKED)
    !In this case, data from different vectors is contiguous. So we can use one message per partition.
    do ipart2 = pv%partno, pv%partno + pv%npart - 1
      ipart = ipart2
      if (ipart > pv%npart) ipart = ipart - pv%npart
      if (pv%ghost_rcounts(ipart) == 0) cycle

      handle%nnb = handle%nnb + 1
      pos = pv%np_local + 1 + pv%ghost_rdispls(ipart)
      call pv%mpi_grp%irecv(v_local%X(ff_pack)(1, pos), pv%ghost_rcounts(ipart)*v_local%pack_size(1), R_MPITYPE, &
        ipart - 1, handle%requests(handle%nnb))
    end do

  case (BATCH_NOT_PACKED)
    do ii = 1, v_local%nst_linear
      do ipart2 = pv%partno, pv%partno + pv%npart - 1
        ipart = ipart2
        if (ipart > pv%npart) ipart = ipart - pv%npart
        if (pv%ghost_rcounts(ipart) == 0) cycle

        handle%nnb = handle%nnb + 1
        pos = pv%np_local + 1 + pv%ghost_rdispls(ipart)
        call pv%mpi_grp%irecv(v_local%X(ff_linear)(pos, ii), pv%ghost_rcounts(ipart), R_MPITYPE, &
          ipart - 1, handle%requests(handle%nnb), tag=ii)
      end do
    end do

  end select
  call profiling_out(prof_irecv)

  call X(batch_init)(handle%ghost_send, v_local%dim, 1, v_local%nst, pv%ghost_scount, &
    packed=v_local%status()==BATCH_PACKED)

  if (v_local%status() == BATCH_DEVICE_PACKED) call handle%ghost_send%do_pack(copy = .false.)

  ! now pack the data for sending
  select case (handle%ghost_send%status())
  case (BATCH_PACKED)
    do ip = 1, pv%ghost_scount
      do ii = 1, handle%ghost_send%nst_linear
        handle%ghost_send%X(ff_pack)(ii, ip) = v_local%X(ff_pack)(ii, pv%ghost_sendmap(ip))
      end do
    end do
  case (BATCH_NOT_PACKED)
    do ii = 1, handle%ghost_send%nst_linear
      do ip = 1, pv%ghost_scount
        handle%ghost_send%X(ff_linear)(ip, ii) = v_local%X(ff_linear)(pv%ghost_sendmap(ip), ii)
      end do
    end do
  case (BATCH_DEVICE_PACKED)
    offset = 0
    call accel_set_kernel_arg(kernel_ghost_reorder, 0, handle%pv%ghost_scount)
    call accel_set_kernel_arg(kernel_ghost_reorder, 1, offset)
    call accel_set_kernel_arg(kernel_ghost_reorder, 2, pv%buff_sendmap)
    call accel_set_kernel_arg(kernel_ghost_reorder, 3, handle%v_local%ff_device)
    call accel_set_kernel_arg(kernel_ghost_reorder, 4, log2(handle%v_local%pack_size_real(1)))
    call accel_set_kernel_arg(kernel_ghost_reorder, 5, handle%ghost_send%ff_device)
    call accel_set_kernel_arg(kernel_ghost_reorder, 6, log2(handle%ghost_send%pack_size_real(1)))

    localsize = accel_kernel_workgroup_size(kernel_ghost_reorder)/handle%ghost_send%pack_size_real(1)

    dim3 = handle%pv%ghost_scount/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(handle%pv%ghost_scount, localsize))

    call accel_kernel_run(kernel_ghost_reorder, &
      (/handle%ghost_send%pack_size_real(1), dim2, dim3/), &
      (/handle%ghost_send%pack_size_real(1), localsize, 1/))

    call accel_finish()
  end select

  if (v_local%status() == BATCH_DEVICE_PACKED) then
    nn = product(handle%ghost_send%pack_size(1:2))
    if (.not. accel%cuda_mpi) then
      SAFE_ALLOCATE(handle%X(send_buffer)(1:nn))
      call accel_read_buffer(handle%ghost_send%ff_device, nn, handle%X(send_buffer))
    else
      call accel_get_device_pointer(handle%X(send_buffer), handle%ghost_send%ff_device, [nn])
    end if
  end if

  call profiling_in(prof_isend, TOSTRING(X(GHOST_UPDATE_ISEND)))
  select case (v_local%status())
  case (BATCH_DEVICE_PACKED)
    ! ring scheme: count downwards from local rank for sending
    do ipart2 = pv%partno, pv%partno - pv%npart + 1, -1
      ipart = ipart2
      if (ipart < 1) ipart = ipart + pv%npart
      if (pv%ghost_scounts(ipart) == 0) cycle
      handle%nnb = handle%nnb + 1
      call pv%mpi_grp%isend(handle%X(send_buffer)(1 + pv%ghost_sdispls(ipart)*v_local%pack_size(1)), &
        pv%ghost_scounts(ipart)*v_local%pack_size(1), R_MPITYPE, ipart - 1, handle%requests(handle%nnb))
    end do

  case (BATCH_PACKED)
    do ipart2 = pv%partno, pv%partno - pv%npart + 1, -1
      ipart = ipart2
      if (ipart < 1) ipart = ipart + pv%npart
      if (pv%ghost_scounts(ipart) == 0) cycle
      handle%nnb = handle%nnb + 1
      call pv%mpi_grp%isend(handle%ghost_send%X(ff_pack)(1, pv%ghost_sdispls(ipart)+1), &
        pv%ghost_scounts(ipart)*v_local%pack_size(1), &
        R_MPITYPE, ipart - 1, handle%requests(handle%nnb))
    end do

  case (BATCH_NOT_PACKED)
    do ii = 1, v_local%nst_linear
      do ipart2 = pv%partno, pv%partno - pv%npart + 1, -1
        ipart = ipart2
        if (ipart < 1) ipart = ipart + pv%npart
        if (pv%ghost_scounts(ipart) == 0) cycle
        handle%nnb = handle%nnb + 1
        call pv%mpi_grp%isend(handle%ghost_send%X(ff_linear)(pv%ghost_sdispls(ipart)+1, ii), &
          pv%ghost_scounts(ipart), R_MPITYPE, ipart - 1, handle%requests(handle%nnb), tag=ii)
      end do
    end do
  end select
  call profiling_out(prof_isend)

  POP_SUB(X(ghost_update_batch_start))
  call profiling_out(prof_start)

end subroutine X(ghost_update_batch_start)

! ---------------------------------------------------------

subroutine X(ghost_update_batch_finish)(handle)
  type(par_vec_handle_batch_t),  intent(inout)   :: handle

  type(profile_t), save :: prof_wait

  call profiling_in(prof_wait, TOSTRING(X(GHOST_UPDATE_WAIT)))
  PUSH_SUB(X(ghost_update_batch_finish))

  ASSERT(handle%nnb > 0)

  call handle%pv%mpi_grp%wait(handle%nnb, handle%requests)
  SAFE_DEALLOCATE_A(handle%requests)

  if (handle%v_local%status() == BATCH_DEVICE_PACKED) then
    ! First wait for the transfer to finish, then call accel_finish to
    ! synchronize the operate_map kernel for the inner points
    call accel_finish()

    ! copy to GPU if not using CUDA aware MPI
    if (.not. accel%cuda_mpi) then
      call accel_write_buffer(handle%v_local%ff_device, handle%v_local%pack_size(1)*handle%pv%np_ghost, &
        handle%X(recv_buffer), handle%v_local%pack_size(1)*handle%pv%np_local)
      SAFE_DEALLOCATE_P(handle%X(send_buffer))
      SAFE_DEALLOCATE_P(handle%X(recv_buffer))
    else
      nullify(handle%X(send_buffer))
      nullify(handle%X(recv_buffer))
    end if
  end if

  call handle%ghost_send%end()

  call profiling_out(prof_wait)
  POP_SUB(X(ghost_update_batch_finish))
end subroutine X(ghost_update_batch_finish)

! ---------------------------------------------------------
!> Set all boundary points in ffb to zero to implement zero
!! boundary conditions for the derivatives, in finite system;
!! or set according to periodic boundary conditions.
subroutine X(boundaries_set_batch)(boundaries, mesh, ffb, phase_correction)
  type(boundaries_t),     intent(in)    :: boundaries
  type(mesh_t),           intent(in)    :: mesh
  class(batch_t), target, intent(inout) :: ffb
  CMPLX,  optional,       intent(in)    :: phase_correction(:)

  integer :: bndry_start, bndry_end
  type(profile_t), save :: set_bc_prof
  type(profile_t), save :: set_bc_comm_prof
  type(profile_t), save :: set_bc_precomm_prof
  type(profile_t), save :: set_bc_postcomm_prof

  PUSH_SUB(X(boundaries_set_batch))
  call profiling_in(set_bc_prof, TOSTRING(X(SET_BC)))

  ASSERT(ffb%type() == R_TYPE_VAL)

  ! If the batch has a phase, this would be wrong to set the boundary conditions
  ! without the phase correction. On the other side, setting the phase correction for a
  ! batch that does not have a phase is also not correct.
  select type(ffb)
  class is(wfs_elec_t)
    if (present(phase_correction)) then
      ASSERT(ffb%has_phase)
    else
      ASSERT(.not. ffb%has_phase)
    end if
  class default
    !In this case we can only assume that the user knows what he is doing
  end select

  ! The boundary points are at different locations depending on the presence
  ! of ghost points due to domain parallelization.
  bndry_start = mesh%np + 1
  bndry_end   = mesh%np_part
  if (mesh%parallel_in_domains) then
    bndry_start = bndry_start + mesh%pv%np_ghost
  end if

  if (.not. boundaries%fully_periodic) call zero_boundaries()
  if (boundaries%periodic) then
    call periodic()
  end if

  call profiling_out(set_bc_prof)
  POP_SUB(X(boundaries_set_batch))

contains

  ! ---------------------------------------------------------
  subroutine zero_boundaries()
    integer :: ist, ip
    integer :: np

    PUSH_SUB(X(boundaries_set_batch).zero_boundaries)

    select case (ffb%status())
    case (BATCH_DEVICE_PACKED)
      np = ffb%pack_size(1)*(bndry_end - bndry_start + 1)
      call accel_set_buffer_to_zero(ffb%ff_device, ffb%type(), np, offset = ffb%pack_size(1)*(bndry_start - 1))
      call accel_finish()

    case (BATCH_PACKED)
      !$omp parallel do private(ist) schedule(static)
      do ip = bndry_start, bndry_end
        !$omp simd
        do ist = 1, ffb%nst_linear
          ffb%X(ff_pack)(ist, ip) = R_TOTYPE(M_ZERO)
        end do
      end do

    case (BATCH_NOT_PACKED)
      do ist = 1, ffb%nst_linear
        !$omp parallel do simd schedule(static)
        do ip = bndry_start, bndry_end
          ffb%X(ff_linear)(ip, ist) = R_TOTYPE(M_ZERO)
        end do
      end do

    end select

    POP_SUB(X(boundaries_set_batch).zero_boundaries)
  end subroutine zero_boundaries

  ! ---------------------------------------------------------
  subroutine periodic()
    integer :: ip, ist, ip_bnd, ip_inn

    R_TYPE, pointer:: sendbuffer(:, :, :)
    R_TYPE, pointer :: recvbuffer(:, :, :)
    integer, allocatable :: send_disp(:), send_count(:)
    integer, allocatable :: recv_disp(:), recv_count(:)
    integer :: ipart, npart, maxsend, maxrecv, ldbuffer, ip2
    type(accel_kernel_t), save :: kernel_send, kernel_recv, kernel_recv_corr, kernel, kernel_corr
    integer :: wgsize
    type(accel_mem_t) :: buff_send
    type(accel_mem_t) :: buff_recv
    type(accel_mem_t) :: buff_phase_corr

    PUSH_SUB(X(boundaries_set_batch).periodic)

    if (mesh%parallel_in_domains) then

      call profiling_in(set_bc_precomm_prof, TOSTRING(X(SET_BC_PRECOMM)))

      npart = mesh%pv%npart
      maxsend = maxval(boundaries%nsend(1:npart))
      maxrecv = maxval(boundaries%nrecv(1:npart))

      ldbuffer = ffb%nst_linear
      if (ffb%status() == BATCH_DEVICE_PACKED) ldbuffer = ffb%pack_size(1)

      if (ffb%status() /= BATCH_DEVICE_PACKED .or. .not. accel%cuda_mpi) then
        SAFE_ALLOCATE(sendbuffer(1:ldbuffer, 1:max(maxsend, 1), 1:npart))
      end if

      select case (ffb%status())

      case (BATCH_NOT_PACKED)

        do ipart = 1, npart
          !$omp parallel do private(ip2, ist)
          do ip = 1, boundaries%nsend(ipart)
            ip2 = boundaries%per_send(ip, ipart)
            do ist = 1, ffb%nst_linear
              sendbuffer(ist, ip, ipart) = ffb%X(ff_linear)(ip2, ist)
            end do
          end do
        end do

      case (BATCH_PACKED)

        do ipart = 1, npart
          !$omp parallel do private(ip2, ist)
          do ip = 1, boundaries%nsend(ipart)
            ip2 = boundaries%per_send(ip, ipart)
            do ist = 1, ffb%nst_linear
              sendbuffer(ist, ip, ipart) = ffb%X(ff_pack)(ist, ip2)
            end do
          end do
        end do

      case (BATCH_DEVICE_PACKED)
        call accel_create_buffer(buff_send, ACCEL_MEM_WRITE_ONLY, R_TYPE_VAL, max(ffb%pack_size(1)*maxsend*npart, 1))

        call accel_kernel_start_call(kernel_send, 'boundaries.cl', 'boundaries_periodic_send')

        call accel_set_kernel_arg(kernel_send, 0, maxsend)
        call accel_set_kernel_arg(kernel_send, 1, boundaries%buff_nsend)
        call accel_set_kernel_arg(kernel_send, 2, boundaries%buff_per_send)
        call accel_set_kernel_arg(kernel_send, 3, ffb%ff_device)
        call accel_set_kernel_arg(kernel_send, 4, log2(ffb%pack_size_real(1)))
        call accel_set_kernel_arg(kernel_send, 5, buff_send)

        wgsize = accel_kernel_workgroup_size(kernel_send)/ffb%pack_size_real(1)

        call accel_kernel_run(kernel_send, (/ffb%pack_size_real(1), pad(maxsend, wgsize), npart/), &
          (/ffb%pack_size_real(1), wgsize, 1/))

        call accel_finish()

        if (.not. accel%cuda_mpi) then
          call accel_read_buffer(buff_send, ffb%pack_size(1)*maxsend*npart, sendbuffer)
          call accel_release_buffer(buff_send)
        else
          call accel_get_device_pointer(sendbuffer, buff_send, [ffb%pack_size(1), max(maxsend, 1), npart])
        end if
      end select


      SAFE_ALLOCATE(send_count(1:npart))
      SAFE_ALLOCATE(send_disp(1:npart))
      SAFE_ALLOCATE(recv_count(1:npart))
      SAFE_ALLOCATE(recv_disp(1:npart))

      do ipart = 1, npart
        send_count(ipart) = ldbuffer*boundaries%nsend(ipart)
        send_disp(ipart)  = ldbuffer*maxsend*(ipart - 1)
        recv_count(ipart) = ldbuffer*boundaries%nrecv(ipart)
        recv_disp(ipart)  = ldbuffer*maxrecv*(ipart - 1)
      end do

      ASSERT(send_count(mesh%pv%partno) == 0)
      ASSERT(recv_count(mesh%pv%partno) == 0)

      if (ffb%status() /= BATCH_DEVICE_PACKED .or. .not. accel%cuda_mpi) then
        SAFE_ALLOCATE(recvbuffer(1:ldbuffer, 1:max(maxrecv, 1), 1:npart))
      else
        call accel_create_buffer(buff_recv, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, max(ffb%pack_size(1)*maxrecv*npart, 1)) 
        ! get device pointer for CUDA-aware MPI
        call accel_get_device_pointer(recvbuffer, buff_recv, [ffb%pack_size(1), max(maxrecv, 1), npart])
      end if

      call profiling_out(set_bc_precomm_prof)

      call profiling_in(set_bc_comm_prof, TOSTRING(X(SET_BC_COMM)))

      call mesh%mpi_grp%alltoallv(sendbuffer, send_count, send_disp, R_MPITYPE, &
        recvbuffer, recv_count, recv_disp, R_MPITYPE)

      ! Only release the buffer after the MPI call
      if(accel%cuda_mpi .and. ffb%status() == BATCH_DEVICE_PACKED) call accel_release_buffer(buff_send)

      call profiling_count_transfers(sum(boundaries%nsend(1:npart) + boundaries%nrecv(1:npart))*ffb%nst_linear, &
        R_TOTYPE(M_ONE))

      call profiling_out(set_bc_comm_prof)

      call profiling_in(set_bc_postcomm_prof, TOSTRING(X(SET_BC_POSTCOMM)))

      SAFE_DEALLOCATE_A(send_count)
      SAFE_DEALLOCATE_A(send_disp)
      SAFE_DEALLOCATE_A(recv_count)
      SAFE_DEALLOCATE_A(recv_disp)
      if(ffb%status() /= BATCH_DEVICE_PACKED .or. .not. accel%cuda_mpi) then
        SAFE_DEALLOCATE_P(sendbuffer)
      else
        nullify(sendbuffer)
      end if

      select case (ffb%status())

      case (BATCH_NOT_PACKED)

        if (.not. present(phase_correction)) then
          ! do not apply phase correction; phase is set in another step
          do ipart = 1, npart
            !$omp parallel do private(ip2, ist)
            do ip = 1, boundaries%nrecv(ipart)
              ip2 = boundaries%per_recv(ip, ipart)
              do ist = 1, ffb%nst_linear
                ffb%X(ff_linear)(ip2, ist) = recvbuffer(ist, ip, ipart)
              end do
            end do
          end do
        else
          ! apply phase correction when setting the BCs -> avoids unnecessary memory access
          ASSERT(lbound(phase_correction, 1) == 1)
          ASSERT(ubound(phase_correction, 1) == mesh%np_part - mesh%np)
          do ipart = 1, npart
            !$omp parallel do private(ip2, ist)
            do ip = 1, boundaries%nrecv(ipart)
              ip2 = boundaries%per_recv(ip, ipart)
              do ist = 1, ffb%nst_linear
#ifdef R_TCOMPLEX
                ffb%zff_linear(ip2, ist) = recvbuffer(ist, ip, ipart) * &
                  phase_correction(ip2-mesh%np)
#else
                ! No phase correction for real batches
                ASSERT(.false.)
#endif
              end do
            end do
          end do
        end if

      case (BATCH_PACKED)

        if (.not. present(phase_correction)) then
          ! do not apply phase correction; phase is set in another step
          do ipart = 1, npart
            !$omp parallel do private(ip2, ist)
            do ip = 1, boundaries%nrecv(ipart)
              ip2 = boundaries%per_recv(ip, ipart)
              do ist = 1, ffb%nst_linear
                ffb%X(ff_pack)(ist, ip2) = recvbuffer(ist, ip, ipart)
              end do
            end do
          end do
        else
          ! apply phase correction when setting the BCs -> avoids unnecessary memory access
          ASSERT(lbound(phase_correction, 1) == 1)
          ASSERT(ubound(phase_correction, 1) == mesh%np_part - mesh%np)
          do ipart = 1, npart
            !$omp parallel do private(ip2, ist)
            do ip = 1, boundaries%nrecv(ipart)
              ip2 = boundaries%per_recv(ip, ipart)
              do ist = 1, ffb%nst_linear
#ifdef R_TCOMPLEX
                ffb%zff_pack(ist, ip2) = recvbuffer(ist, ip, ipart) * &
                  phase_correction(ip2-mesh%np)
#else
                ! No phase correction for real batches
                ASSERT(.false.)
#endif
              end do
            end do
          end do
        end if

      case (BATCH_DEVICE_PACKED)
        if (.not. present(phase_correction)) then
          if (.not. accel%cuda_mpi) then
            call accel_create_buffer(buff_recv, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, max(ffb%pack_size(1)*maxrecv*npart, 1))
            call accel_write_buffer(buff_recv, ffb%pack_size(1)*maxrecv*npart, recvbuffer)
          end if

          call accel_kernel_start_call(kernel_recv, 'boundaries.cl', 'boundaries_periodic_recv')

          call accel_set_kernel_arg(kernel_recv, 0, maxrecv)
          call accel_set_kernel_arg(kernel_recv, 1, boundaries%buff_nrecv)
          call accel_set_kernel_arg(kernel_recv, 2, boundaries%buff_per_recv)
          call accel_set_kernel_arg(kernel_recv, 3, ubound(boundaries%per_recv, dim = 1))
          call accel_set_kernel_arg(kernel_recv, 4, buff_recv)
          call accel_set_kernel_arg(kernel_recv, 5, ffb%ff_device)
          call accel_set_kernel_arg(kernel_recv, 6, log2(ffb%pack_size_real(1)))

          wgsize = accel_kernel_workgroup_size(kernel_recv)/ffb%pack_size_real(1)

          call accel_kernel_run(kernel_recv, (/ffb%pack_size_real(1), pad(maxrecv, wgsize), npart/), &
            (/ffb%pack_size_real(1), wgsize, 1/))

          call accel_finish()

          call accel_release_buffer(buff_recv)
        else
          ASSERT(lbound(phase_correction, 1) == 1)
          ASSERT(ubound(phase_correction, 1) == mesh%np_part - mesh%np)
          ASSERT(R_TYPE_VAL == TYPE_CMPLX)

          if (.not. accel%cuda_mpi) then
            call accel_create_buffer(buff_recv, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, max(ffb%pack_size(1)*maxrecv*npart, 1))
            call accel_write_buffer(buff_recv, ffb%pack_size(1)*maxrecv*npart, recvbuffer)
          end if

          call accel_create_buffer(buff_phase_corr, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, mesh%np_part - mesh%np)
          call accel_write_buffer(buff_phase_corr, mesh%np_part - mesh%np, phase_correction)

          call accel_kernel_start_call(kernel_recv_corr, 'boundaries.cl', 'boundaries_periodic_recv_corr')

          call accel_set_kernel_arg(kernel_recv_corr, 0, maxrecv)
          call accel_set_kernel_arg(kernel_recv_corr, 1, boundaries%buff_nrecv)
          call accel_set_kernel_arg(kernel_recv_corr, 2, boundaries%buff_per_recv)
          call accel_set_kernel_arg(kernel_recv_corr, 3, ubound(boundaries%per_recv, dim = 1))
          call accel_set_kernel_arg(kernel_recv_corr, 4, buff_recv)
          call accel_set_kernel_arg(kernel_recv_corr, 5, ffb%ff_device)
          call accel_set_kernel_arg(kernel_recv_corr, 6, log2(ffb%pack_size(1)))
          call accel_set_kernel_arg(kernel_recv_corr, 7, buff_phase_corr)
          call accel_set_kernel_arg(kernel_recv_corr, 8, mesh%np)

          wgsize = accel_kernel_workgroup_size(kernel_recv_corr)/ffb%pack_size(1)

          call accel_kernel_run(kernel_recv_corr, (/ffb%pack_size(1), pad(maxrecv, wgsize), npart/), &
            (/ffb%pack_size(1), wgsize, 1/))

          call accel_finish()

          call accel_release_buffer(buff_recv)
          call accel_release_buffer(buff_phase_corr)
        end if
      end select

      if(ffb%status() /= BATCH_DEVICE_PACKED .or. .not. accel%cuda_mpi) then
        SAFE_DEALLOCATE_P(recvbuffer)
      else
        nullify(recvbuffer)
      end if

      call profiling_out(set_bc_postcomm_prof)

    end if

    select case (ffb%status())

    case (BATCH_NOT_PACKED)

      if (.not. present(phase_correction)) then
        ! do not apply phase correction; phase is set in another step
        do ist = 1, ffb%nst_linear
          do ip = 1, boundaries%nper
            ffb%X(ff_linear)(boundaries%per_points(POINT_BOUNDARY, ip), ist) = &
              ffb%X(ff_linear)(boundaries%per_points(POINT_INNER, ip), ist)
          end do
        end do
      else
        ! apply phase correction when setting the BCs -> avoids unnecessary memory access
        ASSERT(lbound(phase_correction, 1) == 1)
        ASSERT(ubound(phase_correction, 1) == mesh%np_part - mesh%np)
        do ist = 1, ffb%nst_linear
          do ip = 1, boundaries%nper
#ifdef R_TCOMPLEX
            ffb%X(ff_linear)(boundaries%per_points(POINT_BOUNDARY, ip), ist) = &
              ffb%X(ff_linear)(boundaries%per_points(POINT_INNER, ip), ist) * &
              phase_correction(boundaries%per_points(POINT_BOUNDARY, ip)-mesh%np)
#else
            ! No phase correction for real batches
            ASSERT(.false.)
#endif
          end do
        end do
      end if

    case (BATCH_PACKED)

      if (.not. present(phase_correction)) then
        ! do not apply phase correction; phase is set in another step
        !$omp parallel do private(ip_bnd, ip_inn, ist)
        do ip = 1, boundaries%nper
          ip_bnd = boundaries%per_points(POINT_BOUNDARY, ip)
          ip_inn = boundaries%per_points(POINT_INNER, ip)
          !$omp simd
          do ist = 1, ffb%nst_linear
            ffb%X(ff_pack)(ist, ip_bnd) = ffb%X(ff_pack)(ist, ip_inn)
          end do
        end do
      else
        ! apply phase correction when setting the BCs -> avoids unnecessary memory access
        ASSERT(lbound(phase_correction, 1) == 1)
        ASSERT(ubound(phase_correction, 1) == mesh%np_part - mesh%np)
        !$omp parallel do private(ip_bnd, ip_inn, ist)
        do ip = 1, boundaries%nper
          ip_bnd = boundaries%per_points(POINT_BOUNDARY, ip)
          ip_inn = boundaries%per_points(POINT_INNER, ip)
          !$omp simd
          do ist = 1, ffb%nst_linear
#ifdef R_TCOMPLEX
            ffb%X(ff_pack)(ist, ip_bnd) = ffb%X(ff_pack)(ist, ip_inn) * phase_correction(ip_bnd-mesh%np)
#else
            ! No phase correction for real batches
            ASSERT(.false.)
#endif
          end do
        end do
      end if

    case (BATCH_DEVICE_PACKED)
      if (.not. present(phase_correction)) then
        if (boundaries%nper > 0) then
          call accel_kernel_start_call(kernel, 'boundaries.cl', 'boundaries_periodic')

          call accel_set_kernel_arg(kernel, 0, boundaries%nper)
          call accel_set_kernel_arg(kernel, 1, boundaries%buff_per_points)
          call accel_set_kernel_arg(kernel, 2, ffb%ff_device)
          call accel_set_kernel_arg(kernel, 3, log2(ffb%pack_size_real(1)))

          wgsize = accel_kernel_workgroup_size(kernel)/ffb%pack_size_real(1)

          call accel_kernel_run(kernel, (/ffb%pack_size_real(1), pad(boundaries%nper, wgsize)/), &
            (/ffb%pack_size_real(1), wgsize/))

          call accel_finish()
        end if
      else

        ASSERT(lbound(phase_correction, 1) == 1)
        ASSERT(ubound(phase_correction, 1) == mesh%np_part - mesh%np)
        ASSERT(R_TYPE_VAL == TYPE_CMPLX)

        call accel_create_buffer(buff_phase_corr, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, mesh%np_part - mesh%np)
        call accel_write_buffer(buff_phase_corr, mesh%np_part - mesh%np, phase_correction)

        call accel_kernel_start_call(kernel_corr, 'boundaries.cl', 'boundaries_periodic_corr')

        call accel_set_kernel_arg(kernel_corr, 0, boundaries%nper)
        call accel_set_kernel_arg(kernel_corr, 1, boundaries%buff_per_points)
        call accel_set_kernel_arg(kernel_corr, 2, ffb%ff_device)
        call accel_set_kernel_arg(kernel_corr, 3, log2(ffb%pack_size(1)))
        call accel_set_kernel_arg(kernel_corr, 4, buff_phase_corr)
        call accel_set_kernel_arg(kernel_corr, 5, mesh%np)

        wgsize = accel_kernel_workgroup_size(kernel_corr)/ffb%pack_size(1)

        call accel_kernel_run(kernel_corr, (/ffb%pack_size(1), pad(boundaries%nper, wgsize)/), &
          (/ffb%pack_size(1), wgsize/))

        call accel_finish()

        call accel_release_buffer(buff_phase_corr)

      end if
    end select

    POP_SUB(X(boundaries_set_batch).periodic)
  end subroutine periodic

end subroutine X(boundaries_set_batch)

! ---------------------------------------------------------

subroutine X(boundaries_set_single)(boundaries, mesh, ff, phase_correction)
  type(boundaries_t),         intent(in)    :: boundaries
  type(mesh_t),               intent(in)    :: mesh
  R_TYPE, target, contiguous, intent(inout) :: ff(:)
  CMPLX, optional,            intent(in)    :: phase_correction(:)

  type(batch_t) :: batch_ff

  PUSH_SUB(X(boundaries_set_single))

  call batch_init(batch_ff, ff)

  call X(boundaries_set_batch)(boundaries, mesh, batch_ff, phase_correction=phase_correction)

  call batch_ff%end()
  POP_SUB(X(boundaries_set_single))

end subroutine X(boundaries_set_single)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
