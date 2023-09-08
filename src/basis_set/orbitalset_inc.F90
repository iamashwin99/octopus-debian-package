!! Copyright (C) 2017 N. Tancogne-Dejean
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

! ---------------------------------------------------------
subroutine X(orbitalset_get_coefficients)(os, ndim, psi, ik, has_phase, dot)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  R_TYPE,               intent(in) :: psi(:,:)
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,            intent(inout) :: dot(:,:)

  integer :: im, ip, idim, idim_orb
  type(profile_t), save :: prof, prof_reduce
  R_TYPE, allocatable :: spsi(:,:)
#ifdef R_TCOMPLEX
  CMPLX, allocatable ::tmp(:)
#endif
  logical :: use_submesh

  call profiling_in(prof, TOSTRING(X(ORBSET_GET_COEFFICIENTS)))

  PUSH_SUB(X(orbitalset_get_coefficients))

  ASSERT(ubound(dot, dim=2) >= os%norbs)

  use_submesh = os%use_submesh
  ! Because of possible phase corrections at the border, the array X(orb) is always stored
  ! on the submesh for complex wavefunctions
  ! This does only apply to X(orb). eorb_mesh/eorb_submesh are
  ! still stored according to the user choice given by basis%submesh.
  ! Hence, only if we do not have phases but have a complex wavefunctions we will access X(orb)
  ! always on the submesh
#ifdef R_TCOMPLEX
  if (.not. has_phase) use_submesh = .true.
#endif

  if (use_submesh) then
    SAFE_ALLOCATE(spsi(1:os%sphere%np, 1:ndim))
    do idim = 1, ndim
      !$omp parallel do
      do ip = 1, os%sphere%np
        spsi(ip,idim) = psi(os%sphere%map(ip), idim)
      end do
    end do
  end if

  !If we need to add the phase, we explicitly do the operation using the sphere
  if (has_phase) then
#ifdef R_TCOMPLEX
    if(os%ndim == 1 .and. .not. os%sphere%mesh%use_curvilinear) then
      SAFE_ALLOCATE(tmp(os%norbs))
      do idim = 1, ndim
        if (.not. os%use_submesh) then
          call blas_gemv('C', os%sphere%mesh%np, os%norbs, R_TOTYPE(os%sphere%mesh%volume_element), &
            os%eorb_mesh(1, 1, 1, ik), os%sphere%mesh%np, psi(1, idim), 1, R_TOTYPE(M_ZERO), tmp(1), 1)
        else
          if(os%sphere%np > 0) then
            call blas_gemv('C', os%sphere%np, os%norbs, R_TOTYPE(os%sphere%mesh%volume_element), &
              os%eorb_submesh(1, 1, 1, ik), os%sphere%np, spsi(1, idim), 1, R_TOTYPE(M_ZERO), tmp(1), 1)
          else
            tmp = M_ZERO
          end if
        end if
        dot(idim, 1:os%norbs) = tmp
      end do
      SAFE_DEALLOCATE_A(tmp)
    else
      do im = 1, os%norbs
        do idim = 1, ndim
          idim_orb = min(idim,os%ndim)
          if (.not. os%use_submesh) then
            dot(idim,im) = zmf_dotp(os%sphere%mesh, os%eorb_mesh(1:os%sphere%mesh%np, im, idim_orb, ik),&
              psi(1:os%sphere%mesh%np, idim), reduce=.false.)
          else
            dot(idim, im) = zmf_dotp(os%sphere%mesh, os%eorb_submesh(1:os%sphere%np, idim_orb, im, ik),&
              spsi(1:os%sphere%np, idim), reduce=.false., np=os%sphere%np)
          end if
        end do
      end do
    end if
#endif
  else
    do im = 1, os%norbs
      do idim = 1, ndim
        idim_orb = min(idim,os%ndim)
        if (.not. use_submesh) then
          dot(idim,im) = X(mf_dotp)(os%sphere%mesh, os%X(orb)(1:os%sphere%mesh%np, idim_orb, im),&
            psi(1:os%sphere%mesh%np, idim), reduce=.false.)
        else
          dot(idim,im) = X(mf_dotp)(os%sphere%mesh, os%X(orb)(1:os%sphere%np, idim_orb, im),&
            spsi(1:os%sphere%np, idim), reduce=.false., np=os%sphere%np)
        end if
      end do
    end do
  end if

  if (os%sphere%mesh%parallel_in_domains) then
    call profiling_in(prof_reduce, TOSTRING(X(ORBSET_GET_COEFF_REDUCE)))
    call os%sphere%mesh%allreduce(dot)
    call profiling_out(prof_reduce)
  end if

  SAFE_DEALLOCATE_A(spsi)

  POP_SUB(X(orbitalset_get_coefficients))
  call profiling_out(prof)
end subroutine X(orbitalset_get_coefficients)

! ---------------------------------------------------------
subroutine X(orbitalset_get_coeff_batch)(os, ndim, psib, dot)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  type(wfs_elec_t),     intent(in) :: psib
  R_TYPE,            intent(inout) :: dot(:,:,:) !< idim, iorb, ist

  integer :: ist
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psi(:,:)

  PUSH_SUB(X(orbitalset_get_coeff_batch))

  call profiling_in(prof, TOSTRING(X(ORBSET_GET_COEFF_BATCH)))

  if(accel_is_enabled() .and. ndim == 1 .and. .not. os%sphere%mesh%use_curvilinear) then
    call X(orbitalset_get_coeff_batch_accel)(os, ndim, psib, dot(1:ndim, 1:os%norbs, 1:psib%nst))
  else
    SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:ndim))
    do ist = 1, psib%nst
      call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
      call X(orbitalset_get_coefficients)(os, ndim, psi, psib%ik, psib%has_phase, dot(1:ndim,1:os%norbs,ist))
    end do
    SAFE_DEALLOCATE_A(psi)
  end if


  call profiling_out(prof)
  POP_SUB(X(orbitalset_get_coeff_batch))
end subroutine X(orbitalset_get_coeff_batch)

! ---------------------------------------------------------
subroutine X(orbitalset_get_coeff_batch_accel)(os, ndim, psib, dot)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  type(wfs_elec_t),     intent(in) :: psib
  R_TYPE,            intent(inout) :: dot(:,:,:) !< idim, iorb, ist

  integer :: ist, iorb, np, size
  integer :: global_sizes(3), local_sizes(3)
  type(profile_t), save :: prof, prof_reduce
  type(accel_mem_t) :: buff_dot
#ifdef R_TCOMPLEX
  type(accel_kernel_t), save, target :: ker_proj_bra_cmplx, ker_proj_bra_cmplx_submesh
#else
  type(accel_kernel_t), save, target :: ker_proj_bra, ker_proj_bra_submesh
#endif
  type(accel_kernel_t), pointer :: kernel
  R_TYPE, allocatable :: tmp_dot(:,:)
  logical :: use_submesh

  PUSH_SUB(X(orbitalset_get_coeff_batch_accel))

  ASSERT(ndim == 1)
  ASSERT(psib%status() == BATCH_DEVICE_PACKED)

  call profiling_in(prof, TOSTRING(X(ORBSET_GET_COEFF_ACCEL)))

  use_submesh = os%use_submesh
  ! Because of possible phase corrections at the border, the array X(orb) is always stored
  ! on the submesh for complex wavefunctions
  ! This does only apply to X(orb). eorb_mesh/eorb_submesh are
  ! still stored according to the user choice given by basis%submesh.
  ! Hence, only if we do not have phases but have a complex wavefunctions we will access X(orb)
  ! always on the submesh
#ifdef R_TCOMPLEX
  if (.not. psib%has_phase) use_submesh = .true.
#endif

  size = int(psib%pack_size(1), i4)
  call accel_create_buffer(buff_dot, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, size*os%norbs)

#ifdef R_TCOMPLEX
  if(.not. os%use_submesh) then
    call accel_kernel_start_call(ker_proj_bra_cmplx, 'dftu_projector.cl', 'dftu_projector_bra_cmplx')
    kernel => ker_proj_bra_cmplx
    np = os%sphere%mesh%np
  else
    call accel_kernel_start_call(ker_proj_bra_cmplx_submesh, 'dftu_projector.cl', &
      'dftu_projector_bra_cmplx_submesh')
    kernel => ker_proj_bra_cmplx_submesh
    np = os%sphere%np
  end if
#else
  if(os%use_submesh) then
    call accel_kernel_start_call(ker_proj_bra_submesh, 'dftu_projector.cl', 'dftu_projector_bra_submesh')
    kernel => ker_proj_bra_submesh
    np = os%sphere%np
  else
    call accel_kernel_start_call(ker_proj_bra, 'dftu_projector.cl', 'dftu_projector_bra')
    kernel => ker_proj_bra
    np = os%sphere%mesh%np
  end if
#endif

  call accel_set_kernel_arg(kernel, 0, np)
  call accel_set_kernel_arg(kernel, 1, psib%nst)
  call accel_set_kernel_arg(kernel, 2, os%norbs)
  if (psib%has_phase) then
    call accel_set_kernel_arg(kernel, 3, os%buff_eorb(psib%ik))
  else
    call accel_set_kernel_arg(kernel, 3, os%X(buff_orb))
  end if
  call accel_set_kernel_arg(kernel, 4, log2(os%ldorbs))
  call accel_set_kernel_arg(kernel, 5, psib%ff_device)
  call accel_set_kernel_arg(kernel, 6, log2(size))
  call accel_set_kernel_arg(kernel, 7, buff_dot)
#ifdef R_TCOMPLEX
  if(os%use_submesh) then
#else
  if(use_submesh) then
#endif
    call accel_set_kernel_arg(kernel, 8, os%sphere%buff_map)
  end if

  ! We use an optimized kernel, in which the loop over np is broken
  ! further into chunks, in order to parallelize over the threads within a warp.
  ! Therefore we need to launch warp_size * size kernels. The size of each block needs to
  ! have multiples of warp_size as x-dimension.
  global_sizes = (/ size * accel%warp_size, pad_pow2(os%norbs), 1 /)
  local_sizes  = (/ accel%warp_size,  1, 1 /)

  call accel_kernel_run(kernel, global_sizes, local_sizes)

  call accel_finish()

  SAFE_ALLOCATE(tmp_dot(1:size, 1:os%norbs))
  call accel_read_buffer(buff_dot, size*os%norbs, tmp_dot)
  call accel_release_buffer(buff_dot)

  if (os%sphere%mesh%parallel_in_domains) then
    call profiling_in(prof_reduce, TOSTRING(X(ORBSET_GET_COEFF_REDUCE)))
    call os%sphere%mesh%allreduce(tmp_dot, dim = (/psib%nst, os%norbs/))
    call profiling_out(prof_reduce)
  end if

  do ist = 1, psib%nst
    do iorb = 1, os%norbs
      dot(1, iorb, ist) = tmp_dot(ist, iorb)*os%sphere%mesh%volume_element
    end do
  end do

  SAFE_DEALLOCATE_A(tmp_dot)

  call profiling_out(prof)
  POP_SUB(X(orbitalset_get_coeff_batch_accel))
end subroutine X(orbitalset_get_coeff_batch_accel)

! ---------------------------------------------------------
subroutine X(orbitalset_add_to_batch)(os, ndim, psib, weight)
  type(orbitalset_t),   intent(in)    :: os
  integer,              intent(in)    :: ndim
  type(wfs_elec_t),     intent(inout) :: psib
  R_TYPE,               intent(in)    :: weight(:,:) !(os%norbs, psib%nst_linear)

  integer :: ip, iorb, ist, idim, bind, idim_orb
  integer :: idim1, idim2, idim3, idim4
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psi(:,:), sorb(:), tmp_weights(:,:)
  integer :: block_size, size, sp, ep
  logical :: use_submesh
  integer :: dim2, dim3, localsize
  type(accel_mem_t) :: buff_weight
#ifdef R_TCOMPLEX
  type(accel_kernel_t), save, target :: ker_proj_ket_cmplx, ker_proj_ket_cmplx_submesh
#else
  type(accel_kernel_t), save, target :: ker_proj_ket, ker_proj_ket_submesh
#endif
  type(accel_kernel_t), pointer :: kernel

  PUSH_SUB(X(orbitalset_add_to_batch))

  call profiling_in(prof, TOSTRING(X(ORBSET_ADD_TO_BATCH)))

  use_submesh = os%use_submesh
  ! Because of possible phase corrections at the border, the array X(orb) is always stored
  ! on the submesh for complex wavefunctions
  ! This does only apply to X(orb). eorb_mesh/eorb_submesh are
  ! still stored according to the user choice given by basis%submesh.
  ! Hence, only if we do not have phases but have a complex wavefunctions we will access X(orb)
  ! always on the submesh
#ifdef R_TCOMPLEX
  if (.not. psib%has_phase) use_submesh = .true.
#endif

  ! This routine uses blocking to optimize cache usage.
  block_size = hardware%X(block_size)

  if (os%sphere%mesh%use_curvilinear .or. (psib%status() == BATCH_DEVICE_PACKED .and. ndim > 1)) then
    !
    SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:ndim))
    do ist = 1, psib%nst
      call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
      !In case of phase, we have to apply the conjugate of the phase here
      if (psib%has_phase) then
#ifdef R_TCOMPLEX
        do idim = 1, ndim
          idim_orb = min(idim,os%ndim)
          bind = psib%ist_idim_to_linear((/ist, idim/))
          do iorb = 1, os%norbs
            if (.not. os%use_submesh) then
              call lalg_axpy(os%sphere%mesh%np, weight(iorb, bind), os%eorb_mesh(1:os%sphere%mesh%np, iorb, idim_orb, psib%ik), &
                psi(1:os%sphere%mesh%np, idim))
            else
              call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np, idim_orb, iorb, psib%ik), &
                psi(1:os%sphere%mesh%np,idim), weight(iorb,bind))
            end if
          end do
        end do
#endif
      else
        do iorb = 1, os%norbs
          do idim = 1, ndim
            idim_orb = min(idim,os%ndim)
            bind = psib%ist_idim_to_linear((/ist, idim/))
            if (.not. use_submesh) then
              call lalg_axpy(os%sphere%mesh%np, weight(iorb, bind), os%X(orb)(1:os%sphere%mesh%np,idim_orb,iorb), &
                psi(1:os%sphere%mesh%np,idim))
            else
              call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np, idim_orb, iorb), &
                psi(1:os%sphere%mesh%np,idim), weight(iorb, bind))
            end if
          end do
        end do
      end if
      call batch_set_state(psib, ist, os%sphere%mesh%np, psi)
    end do
    SAFE_DEALLOCATE_A(psi)
    !
  else
    !
    select case (psib%status())
    case (BATCH_NOT_PACKED)
      !
      if (psib%has_phase) then
#ifdef R_TCOMPLEX
        if (.not. os%use_submesh) then
          do sp = 1, os%sphere%mesh%np, block_size
            size = min(block_size, os%sphere%mesh%np - sp + 1)
            do ist = 1, psib%nst_linear
              idim = min(psib%linear_to_idim(ist), os%ndim)
              call blas_gemv('N', size, os%norbs, R_TOTYPE(M_ONE), os%eorb_mesh(sp, 1, idim, psib%ik), &
                os%sphere%mesh%np, weight(1, ist), 1, R_TOTYPE(M_ONE), psib%zff_linear(sp, ist), 1)
            end do
          end do
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist),os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do sp = 1, os%sphere%np, block_size
              size = min(block_size, os%sphere%np - sp + 1)
              ep = sp - 1 + size
              sorb(sp:ep) = R_TOTYPE(M_ZERO)
              do iorb = 1, os%norbs
                call blas_axpy(size, weight(iorb, ist), os%eorb_submesh(sp, idim, iorb, psib%ik), 1, sorb(sp), 1)
              end do
            end do
            do ip = 1,os%sphere%np
              psib%zff_linear(os%sphere%map(ip), ist) = &
                psib%zff_linear(os%sphere%map(ip), ist) + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
#endif
      else
        if (.not. use_submesh) then
          do sp = 1, os%sphere%mesh%np, block_size
            size = min(block_size, os%sphere%mesh%np - sp + 1)
            do ist = 1, psib%nst_linear
              idim = min(psib%linear_to_idim(ist), os%ndim)
              call blas_gemv('N', size, os%norbs, R_TOTYPE(M_ONE), os%X(orb)(sp, idim, 1), &
                os%sphere%mesh%np * os%ndim, weight(1, ist), 1, R_TOTYPE(M_ONE),            &
                psib%X(ff_linear)(sp, ist), 1)
            end do
          end do
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist), os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              call lalg_axpy(os%sphere%np, weight(iorb, ist), os%X(orb)(:, idim, iorb), sorb)
            end do
            call submesh_add_to_mesh(os%sphere, sorb, psib%X(ff_linear)(:,ist))
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
      end if

    case (BATCH_PACKED)
      !
      if (psib%has_phase) then
#ifdef R_TCOMPLEX
        if (.not. os%use_submesh) then
          !$omp parallel private(sp, ep, iorb, ist, idim, idim1, idim2, idim3, idim4, ip)
          do sp = 1, os%sphere%mesh%np, block_size
            ep = sp - 1 + min(block_size, os%sphere%mesh%np - sp + 1)
            do iorb = 1, os%norbs
              do ist = 1, psib%nst_linear - 4 + 1, 4
                idim1 = min(psib%linear_to_idim(ist),   os%ndim)
                idim2 = min(psib%linear_to_idim(ist+1), os%ndim)
                idim3 = min(psib%linear_to_idim(ist+2), os%ndim)
                idim4 = min(psib%linear_to_idim(ist+3), os%ndim)

                !$omp do
                do ip = sp, ep
                  psib%zff_pack(ist, ip) = &
                    psib%zff_pack(ist, ip) + weight(iorb, ist)*os%eorb_mesh(ip, iorb, idim1, psib%ik)
                  psib%zff_pack(ist + 1, ip) = &
                    psib%zff_pack(ist + 1, ip) + weight(iorb, ist + 1)*os%eorb_mesh(ip, iorb, idim2, psib%ik)
                  psib%zff_pack(ist + 2, ip) = &
                    psib%zff_pack(ist + 2, ip) + weight(iorb, ist + 2)*os%eorb_mesh(ip, iorb, idim3, psib%ik)
                  psib%zff_pack(ist + 3, ip) = &
                    psib%zff_pack(ist + 3, ip) + weight(iorb, ist + 3)*os%eorb_mesh(ip, iorb, idim4, psib%ik)
                end do
              end do

              do ist = ist, psib%nst_linear
                idim = min(psib%linear_to_idim(ist), os%ndim)
                !$omp do
                do ip = sp, ep
                  psib%zff_pack(ist, ip) = &
                    psib%zff_pack(ist, ip) + weight(iorb, ist)*os%eorb_mesh(ip, iorb, idim, psib%ik)
                end do
              end do
            end do
          end do
          !$omp end parallel
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist), os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              call lalg_axpy(os%sphere%np, weight(iorb, ist), os%eorb_submesh(:, idim, iorb, psib%ik), sorb)
            end do
            do ip = 1, os%sphere%np
              psib%zff_pack(ist, os%sphere%map(ip)) = psib%zff_pack(ist, os%sphere%map(ip)) &
                + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
#endif
      else
        if (.not. use_submesh) then
          !$omp parallel private(iorb, sp, ep, ist, idim, idim1, idim2, idim3, idim4, ip)
          do iorb = 1, os%norbs
            do sp = 1, os%sphere%mesh%np, block_size
              ep = sp - 1 + min(block_size, os%sphere%mesh%np - sp + 1)
              do ist = 1, psib%nst_linear - 4 + 1, 4
                idim1 = min(psib%linear_to_idim(ist)  , os%ndim)
                idim2 = min(psib%linear_to_idim(ist+1), os%ndim)
                idim3 = min(psib%linear_to_idim(ist+2), os%ndim)
                idim4 = min(psib%linear_to_idim(ist+3), os%ndim)

                !$omp do
                do ip = sp, ep
                  psib%X(ff_pack)(ist, ip) = &
                    psib%X(ff_pack)(ist, ip)   + weight(iorb, ist  ) * os%X(orb)(ip, idim1, iorb)
                  psib%X(ff_pack)(ist+1, ip) = &
                    psib%X(ff_pack)(ist+1, ip) + weight(iorb, ist+1) * os%X(orb)(ip, idim2, iorb)
                  psib%X(ff_pack)(ist+2, ip) = &
                    psib%X(ff_pack)(ist+2, ip) + weight(iorb, ist+2) * os%X(orb)(ip, idim3, iorb)
                  psib%X(ff_pack)(ist+3, ip) = &
                    psib%X(ff_pack)(ist+3, ip) + weight(iorb, ist+3) * os%X(orb)(ip, idim4, iorb)
                end do
              end do

              do ist = ist, psib%nst_linear
                idim = min(psib%linear_to_idim(ist), os%ndim)
                !$omp do
                do ip = sp, ep
                  psib%X(ff_pack)(ist, ip) = &
                    psib%X(ff_pack)(ist, ip) + weight(iorb, ist) * os%X(orb)(ip, idim, iorb)
                end do
              end do
            end do
          end do
          !$omp end parallel
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist), os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              call lalg_axpy(os%sphere%np, weight(iorb, ist), os%X(orb)(:, idim, iorb), sorb)
            end do
            do ip = 1, os%sphere%np
              psib%X(ff_pack)(ist,os%sphere%map(ip)) = psib%X(ff_pack)(ist,os%sphere%map(ip)) &
                + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
      end if

    case(BATCH_DEVICE_PACKED)

      size = int(psib%pack_size(1), i4)

      ! write the weights to the GPU
      call accel_create_buffer(buff_weight, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, os%norbs*size)
      SAFE_ALLOCATE(tmp_weights(1:size, 1:os%norbs))
      tmp_weights = M_ZERO
      do ist = 1, psib%nst_linear
        do iorb = 1, os%norbs
          tmp_weights(ist, iorb) = weight(iorb, ist)
        end do
      end do
      call accel_write_buffer(buff_weight, size*os%norbs, tmp_weights)
      SAFE_DEALLOCATE_A(tmp_weights)

#ifdef R_TCOMPLEX
      if(.not. use_submesh) then
        call accel_kernel_start_call(ker_proj_ket_cmplx, 'dftu_projector.cl', 'dftu_projector_ket_cmplx')
        kernel => ker_proj_ket_cmplx
        call zaccel_add_to_batch_mesh()
      else
        call accel_kernel_start_call(ker_proj_ket_cmplx_submesh, 'dftu_projector.cl', &
          'dftu_projector_ket_cmplx_submesh')
        kernel => ker_proj_ket_cmplx_submesh
        call zaccel_add_to_batch_submesh()
      end if
#else
      if(use_submesh) then
        call accel_kernel_start_call(ker_proj_ket_submesh, 'dftu_projector.cl', 'dftu_projector_ket_submesh')
        kernel => ker_proj_ket_submesh
        call daccel_add_to_batch_submesh()
      else
        call accel_kernel_start_call(ker_proj_ket, 'dftu_projector.cl', 'dftu_projector_ket')
        kernel => ker_proj_ket
        call daccel_add_to_batch_mesh()
      end if
#endif

      call accel_release_buffer(buff_weight)

    end select
  end if

  call profiling_out(prof)

  POP_SUB(X(orbitalset_add_to_batch))

contains
  subroutine X(accel_add_to_batch_mesh)()
    PUSH_SUB(X(orbitalset_add_to_batch.accel_mesh))

    call accel_set_kernel_arg(kernel, 0, os%norbs)
    call accel_set_kernel_arg(kernel, 1, os%sphere%mesh%np)
    call accel_set_kernel_arg(kernel, 2, buff_weight)
    if (psib%has_phase) then
      call accel_set_kernel_arg(kernel, 3, os%buff_eorb(psib%ik))
    else
      call accel_set_kernel_arg(kernel, 3, os%X(buff_orb))
    end if
    call accel_set_kernel_arg(kernel, 4, log2(os%ldorbs))
    call accel_set_kernel_arg(kernel, 5, psib%ff_device)
    call accel_set_kernel_arg(kernel, 6, log2(size))

    localsize = accel_kernel_workgroup_size(kernel)/size
    dim3 = os%sphere%mesh%np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(os%sphere%mesh%np, localsize))

    call accel_kernel_run(kernel, (/size, dim2, dim3/), (/size, localsize, 1/))

    call accel_finish()

    POP_SUB(X(orbitalset_add_to_batch.accel_mesh))
  end subroutine X(accel_add_to_batch_mesh)

  subroutine X(accel_add_to_batch_submesh)()
    integer :: ir

    PUSH_SUB(X(orbitalset_add_to_batch.accel_submesh))

    do ir = 1, os%sphere%num_regions
      call accel_set_kernel_arg(kernel, 0, os%norbs)
      call accel_set_kernel_arg(kernel, 1, os%sphere%regions(ir)-1)
      call accel_set_kernel_arg(kernel, 2, os%sphere%regions(ir+1)-1)
      call accel_set_kernel_arg(kernel, 3, buff_weight)
      if (psib%has_phase) then
        call accel_set_kernel_arg(kernel, 4, os%buff_eorb(psib%ik))
      else
        call accel_set_kernel_arg(kernel, 4, os%X(buff_orb))
      end if
      call accel_set_kernel_arg(kernel, 5, log2(os%ldorbs))
      call accel_set_kernel_arg(kernel, 6, psib%ff_device)
      call accel_set_kernel_arg(kernel, 7, log2(size))
      call accel_set_kernel_arg(kernel, 8, os%sphere%buff_map)

      localsize = accel_kernel_workgroup_size(kernel)/size
      dim3 = os%sphere%np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(os%sphere%np, localsize))

      call accel_kernel_run(kernel, (/size, dim2, dim3/), (/size, localsize, 1/))

      call accel_finish()
    end do

    POP_SUB(X(orbitalset_add_to_batch.accel_submesh))
  end subroutine X(accel_add_to_batch_submesh)
end subroutine X(orbitalset_add_to_batch)


! ---------------------------------------------------------
! This routine is used in the periodic case, as r|\psi> is not periodic
! The matrix element of <phi_m|r|\psi> are therefore computed on the submesh
subroutine X(orbitalset_get_position_matrix_elem)(os, ndim, psib, idir, dot)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  type(wfs_elec_t),     intent(in) :: psib
  integer,              intent(in) :: idir
  R_TYPE,            intent(inout) :: dot(:,:,:)

  integer :: im, ip, idim, idim_orb, ist, indb
  type(profile_t), save :: prof, prof_reduce
  R_TYPE, allocatable :: spsi(:,:)
  logical :: use_submesh
  R_TYPE  :: local_dot
  integer :: global_sizes(3), local_sizes(3), size
  type(accel_mem_t) :: buff_dot, buff_xx
#ifdef R_TCOMPLEX
  type(accel_kernel_t), save, target :: ker_pos_mat_elem_cmplx, ker_pos_mat_elem_cmplx_submesh
  type(accel_kernel_t), save, target :: ker_pos_mat_elem_phase
  type(accel_mem_t) :: buff_phase
#else
  type(accel_kernel_t), save, target :: ker_pos_mat_elem, ker_pos_mat_elem_submesh
#endif
  type(accel_kernel_t), pointer :: kernel
  R_TYPE, allocatable :: tmp_dot(:,:)
  FLOAT, allocatable :: xx(:)


  call profiling_in(prof, TOSTRING(X(ORBSET_GET_POS_MAT_ELEM)))

  PUSH_SUB(X(orbitalset_get_pos_mat_elem))

  ASSERT(ubound(dot, dim=2) >= os%norbs)

  use_submesh = os%use_submesh
  ! Because of possible phase corrections at the border, the array X(orb) is always stored
  ! on the submesh for complex wavefunctions
#ifdef R_TCOMPLEX
  if (.not. psib%has_phase) use_submesh = .true.
#endif

  SAFE_ALLOCATE(xx(os%sphere%np))
  !$omp parallel do
  do ip = 1, os%sphere%np
    xx(ip) = os%sphere%rel_x(idir, ip)+os%sphere%center(idir)
  end do

  if(psib%status() /= BATCH_DEVICE_PACKED) then
    SAFE_ALLOCATE(spsi(1:os%sphere%np, 1:ndim))
    do ist = 1, psib%nst
      do idim = 1, ndim
        indb = psib%ist_idim_to_linear((/ist, idim/))
        if(psib%status() == BATCH_PACKED) then
          !$omp parallel do
          do ip = 1, os%sphere%np
            spsi(ip,idim) = xx(ip)*psib%X(ff_pack)(indb, os%sphere%map(ip))
          end do
        else
          !$omp parallel do
          do ip = 1, os%sphere%np
            spsi(ip,idim) = xx(ip)*psib%X(ff_linear)(os%sphere%map(ip), indb)
          end do
        end if
      end do

      if (psib%has_phase) then
#ifdef R_TCOMPLEX
        do im = 1, os%norbs
          do idim = 1, ndim
            idim_orb = min(idim,os%ndim)
            if (.not. use_submesh) then
              local_dot = M_ZERO
              ! Here we must compute the dot product on the full mesh and hence reconstruct the
              ! the full orbital (with phase correction) on the fly
              !$omp parallel do reduction(+:local_dot)
              do ip = 1, os%sphere%np
                local_dot = local_dot &
                  + conjg(os%zorb(ip,idim_orb,im)*os%phase(ip, psib%ik)) * spsi(ip, idim)
              end do
              dot(idim, im, ist) = local_dot * os%sphere%mesh%volume_element
            else
              dot(idim, im, ist) = X(mf_dotp)(os%sphere%mesh, os%eorb_submesh(1:os%sphere%np, idim_orb, im, psib%ik),&
                spsi(1:os%sphere%np, idim), reduce = .false., np = os%sphere%np)
            end if
          end do
        end do
#endif
      else
        do im = 1, os%norbs
          do idim = 1, ndim
            idim_orb = min(idim,os%ndim)
            if (.not. use_submesh) then
              local_dot = M_ZERO
              !$omp parallel do reduction(+:local_dot)
              do ip = 1, os%sphere%np
                local_dot = local_dot &
                  + R_CONJ(os%X(orb)(os%sphere%map(ip), idim_orb, im)) * spsi(ip, idim)
              end do
              dot(idim, im, ist) = local_dot * os%sphere%mesh%volume_element
            else
              dot(idim, im, ist) = X(mf_dotp)(os%sphere%mesh, os%X(orb)(1:os%sphere%np, idim_orb, im),&
                spsi(1:os%sphere%np, idim), reduce = .false., np = os%sphere%np)
            end if
          end do !idim
        end do !im
      end if

    end do !ist

  else

    size = int(psib%pack_size(1), i4)
    call accel_create_buffer(buff_dot, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, size*os%norbs)

    call accel_create_buffer(buff_xx, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, os%sphere%np)
    call accel_write_buffer(buff_xx, os%sphere%np, xx)

#ifdef R_TCOMPLEX
    if(.not. os%use_submesh) then
      if(.not. psib%has_phase) then
        call accel_kernel_start_call(ker_pos_mat_elem_cmplx, 'dftu_projector.cl', 'dftu_pos_mat_elem_cmplx')
        kernel => ker_pos_mat_elem_cmplx
      else
        call accel_kernel_start_call(ker_pos_mat_elem_phase, 'dftu_projector.cl', 'dftu_pos_mat_elem_phase')
        kernel => ker_pos_mat_elem_phase
        call accel_create_buffer(buff_phase, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, os%sphere%np)
        call accel_write_buffer(buff_phase, os%sphere%np, os%phase(:,psib%ik))
      end if
    else
      call accel_kernel_start_call(ker_pos_mat_elem_cmplx_submesh, 'dftu_projector.cl', &
        'dftu_pos_mat_elem_cmplx_submesh')
      kernel => ker_pos_mat_elem_cmplx_submesh
    end if
#else
    if(os%use_submesh) then
      call accel_kernel_start_call(ker_pos_mat_elem_submesh, 'dftu_projector.cl', 'dftu_pos_mat_elem_submesh')
      kernel => ker_pos_mat_elem_submesh
    else
      call accel_kernel_start_call(ker_pos_mat_elem, 'dftu_projector.cl', 'dftu_pos_mat_elem')
      kernel => ker_pos_mat_elem
    end if
#endif

    call accel_set_kernel_arg(kernel, 0, os%sphere%np)
    call accel_set_kernel_arg(kernel, 1, psib%nst)
    call accel_set_kernel_arg(kernel, 2, os%norbs)
    if (psib%has_phase) then
      call accel_set_kernel_arg(kernel, 3, os%buff_eorb(psib%ik))
    else
      call accel_set_kernel_arg(kernel, 3, os%X(buff_orb))
    end if
    call accel_set_kernel_arg(kernel, 4, log2(os%ldorbs))
    call accel_set_kernel_arg(kernel, 5, psib%ff_device)
    call accel_set_kernel_arg(kernel, 6, log2(size))
    call accel_set_kernel_arg(kernel, 7, buff_xx)
    call accel_set_kernel_arg(kernel, 8, os%sphere%buff_map)
    call accel_set_kernel_arg(kernel, 9, buff_dot)
#ifdef R_TCOMPLEX
    if(.not. os%use_submesh .and. psib%has_phase) then
      call accel_set_kernel_arg(kernel, 10, buff_phase)
    end if
#endif

    ! We use an optimized kernel, in which the loop over np is broken
    ! further into chunks, in order to parallelize over the threads within a warp.
    ! Therefore we need to launch warp_size * size kernels. The size of each block needs to
    ! have multiples of warp_size as x-dimension.
    global_sizes = (/ size * accel%warp_size, pad_pow2(os%norbs), 1 /)
    local_sizes  = (/ accel%warp_size,  1, 1 /)

    call accel_kernel_run(kernel, global_sizes, local_sizes)

    call accel_finish()

    SAFE_ALLOCATE(tmp_dot(1:size, 1:os%norbs))
    call accel_read_buffer(buff_dot, size*os%norbs, tmp_dot)
    call accel_release_buffer(buff_dot)
    tmp_dot = tmp_dot * os%sphere%mesh%volume_element

    do ist = 1, psib%nst
      do im = 1, os%norbs
        dot(1, im, ist) = tmp_dot(ist, im)
      end do
    end do

    call accel_release_buffer(buff_xx)
#ifdef R_TCOMPLEX
    call accel_release_buffer(buff_phase)
#endif
  end if

  SAFE_DEALLOCATE_A(xx)

  if (os%sphere%mesh%parallel_in_domains) then
    call profiling_in(prof_reduce, TOSTRING(X(ORBSET_POS_MAT_REDUCE)))
    call os%sphere%mesh%allreduce(dot)
    call profiling_out(prof_reduce)
  end if

  SAFE_DEALLOCATE_A(spsi)

  POP_SUB(X(orbitalset_get_pos_mat_elem))
  call profiling_out(prof)
end subroutine X(orbitalset_get_position_matrix_elem)
