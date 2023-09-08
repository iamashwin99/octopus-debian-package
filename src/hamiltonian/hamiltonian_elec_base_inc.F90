!! Copyright (C) 2009-2020 X. Andrade, N. Tancogne-Dejean, M. Lueders
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

subroutine X(hamiltonian_elec_base_local)(this, mesh, std, ispin, psib, vpsib)
  type(hamiltonian_elec_base_t),  intent(in)    :: this
  class(mesh_t),                  intent(in)    :: mesh
  type(states_elec_dim_t),        intent(in)    :: std
  integer,                        intent(in)    :: ispin
  type(wfs_elec_t),               intent(in)    :: psib
  type(wfs_elec_t),               intent(inout) :: vpsib

  PUSH_SUB(X(hamiltonian_elec_base_local))

  if (psib%status() == BATCH_DEVICE_PACKED) then
    if (allocated(this%Impotential)) then
      call X(hamiltonian_elec_base_local_sub)(this%potential, mesh, std, ispin, &
        psib, vpsib, potential_accel=this%potential_accel, &
        impotential_accel=this%impotential_accel)
    else
      call X(hamiltonian_elec_base_local_sub)(this%potential, mesh, std, ispin, &
        psib, vpsib, potential_accel = this%potential_accel)
    end if
  else
    if (allocated(this%Impotential)) then
      call X(hamiltonian_elec_base_local_sub)(this%potential, mesh, std, ispin, &
        psib, vpsib, Impotential = this%Impotential)
    else
      call X(hamiltonian_elec_base_local_sub)(this%potential, mesh, std, ispin, psib, vpsib)
    end if
  end if

  POP_SUB(X(hamiltonian_elec_base_local))
end subroutine X(hamiltonian_elec_base_local)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_local_sub)(potential, mesh, std, ispin, psib, vpsib, &
  Impotential, potential_accel, impotential_accel)
  FLOAT,                             intent(in)    :: potential(:,:)
  class(mesh_t),                     intent(in)    :: mesh
  type(states_elec_dim_t),           intent(in)    :: std
  integer,                           intent(in)    :: ispin
  class(batch_t), target,            intent(in)    :: psib
  class(batch_t), target,            intent(inout) :: vpsib
  FLOAT, optional,                   intent(in)    :: impotential(:,:)
  type(accel_mem_t),  optional, target, intent(in) :: potential_accel
  type(accel_mem_t),  optional, target, intent(in) :: impotential_accel

  integer :: ist, ip, is
#ifdef R_TCOMPLEX
  R_TYPE :: psi1, psi2
  FLOAT  :: Imvv
  CMPLX  :: pot(1:4)
  CMPLX, pointer :: psi(:, :), vpsi(:, :)
#endif
  FLOAT   :: vv
  logical :: pot_is_cmplx
  integer :: pnp
  integer(i8) :: localsize, dim2, dim3
  type(accel_mem_t), pointer :: potential_accel_, impotential_accel_

  call profiling_in(prof_vlpsi, TOSTRING(X(VLPSI)))
  PUSH_SUB(X(hamiltonian_elec_base_local_sub))

  pot_is_cmplx = .false.
  if (present(Impotential)) pot_is_cmplx = .true.
  if (present(impotential_accel)) pot_is_cmplx = .true.

  call psib%check_compatibility_with(vpsib)

  select case (psib%status())
  case (BATCH_DEVICE_PACKED)
    pnp = accel_padded_size(mesh%np)

    if(.not. present(potential_accel)) then
      SAFE_ALLOCATE(potential_accel_)
      call accel_create_buffer(potential_accel_, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, pnp*std%nspin)
      if(std%ispin /= SPINORS) then
        call accel_write_buffer(potential_accel_, mesh%np, potential(:,ispin), offset=pnp*(ispin - 1))
      else
        do is = 1, std%nspin
          call accel_write_buffer(potential_accel_, mesh%np, potential(:,is), offset=pnp*(is - 1))
        end do
      end if
    else
      potential_accel_ => potential_accel
    end if

    if(.not. present(impotential_accel) .and. present(impotential)) then
      SAFE_ALLOCATE(impotential_accel_)
      call accel_create_buffer(impotential_accel_, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, pnp*std%nspin)
      if(std%ispin /= SPINORS) then
        call accel_write_buffer(impotential_accel_, mesh%np, impotential(:,ispin), offset=pnp*(ispin - 1))
      else
        do is = 1, std%nspin
          call accel_write_buffer(impotential_accel_, mesh%np, impotential(:,is), offset=pnp*(is - 1))
        end do
      end if
    else
      impotential_accel_ => impotential_accel
    end if

    if (.not. pot_is_cmplx) then
      select case (std%ispin)

      case (UNPOLARIZED, SPIN_POLARIZED)
        call accel_set_kernel_arg(kernel_vpsi, 0, pnp*(ispin - 1))
        call accel_set_kernel_arg(kernel_vpsi, 1, mesh%np)
        call accel_set_kernel_arg(kernel_vpsi, 2, potential_accel_)
        call accel_set_kernel_arg(kernel_vpsi, 3, psib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi, 4, int(log2(psib%pack_size_real(1)), i4))
        call accel_set_kernel_arg(kernel_vpsi, 5, vpsib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi, 6, int(log2(vpsib%pack_size_real(1)), i4))

        localsize = accel_kernel_workgroup_size(kernel_vpsi)/psib%pack_size_real(1)

        dim3 = mesh%np/(accel_max_size_per_dim(2)*localsize) + 1
        dim2 = min(accel_max_size_per_dim(2)*localsize, pad(mesh%np, localsize))

        call accel_kernel_run(kernel_vpsi, (/psib%pack_size_real(1), dim2, dim3/), (/psib%pack_size_real(1), localsize, 1_i8/))

      case (SPINORS)
        call accel_set_kernel_arg(kernel_vpsi_spinors, 0, mesh%np)
        call accel_set_kernel_arg(kernel_vpsi_spinors, 1, potential_accel_)
        call accel_set_kernel_arg(kernel_vpsi_spinors, 2, pnp)
        call accel_set_kernel_arg(kernel_vpsi_spinors, 3, psib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi_spinors, 4, int(psib%pack_size(1), i4))
        call accel_set_kernel_arg(kernel_vpsi_spinors, 5, vpsib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi_spinors, 6, int(vpsib%pack_size(1), i4))

        localsize = accel_kernel_workgroup_size(kernel_vpsi_spinors)/(psib%pack_size(1)/2)

        dim3 = pnp/(accel_max_size_per_dim(2)*localsize) + 1
        dim2 = min(accel_max_size_per_dim(2)*localsize, pad(pnp, localsize))

        call accel_kernel_run(kernel_vpsi_spinors, (/psib%pack_size(1)/2, dim2, dim3/), &
          (/psib%pack_size(1)/2, localsize, 1_i8/))

      end select
    else
      ! complex potentials
      select case (std%ispin)

      case (UNPOLARIZED, SPIN_POLARIZED)
        call accel_set_kernel_arg(kernel_vpsi_complex, 0, pnp*(ispin - 1))
        call accel_set_kernel_arg(kernel_vpsi_complex, 1, mesh%np)
        call accel_set_kernel_arg(kernel_vpsi_complex, 2, potential_accel_)
        call accel_set_kernel_arg(kernel_vpsi_complex, 3, impotential_accel)
        call accel_set_kernel_arg(kernel_vpsi_complex, 4, psib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi_complex, 5, int(log2(psib%pack_size(1)), i4))
        call accel_set_kernel_arg(kernel_vpsi_complex, 6, vpsib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi_complex, 7, int(log2(vpsib%pack_size(1)), i4))

        localsize = accel_kernel_workgroup_size(kernel_vpsi_complex)/psib%pack_size(1)

        dim3 = mesh%np/(accel_max_size_per_dim(2)*localsize) + 1
        dim2 = min(accel_max_size_per_dim(2)*localsize, pad(mesh%np, localsize))

        call accel_kernel_run(kernel_vpsi_complex, (/psib%pack_size(1), dim2, dim3/), (/psib%pack_size(1), localsize, 1_i8/))

      case (SPINORS)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 0, mesh%np)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 1, potential_accel_)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 2, pnp)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 3, impotential_accel)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 4, pnp)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 5, psib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 6, int(psib%pack_size(1), i4))
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 7, vpsib%ff_device)
        call accel_set_kernel_arg(kernel_vpsi_spinors_complex, 8, int(vpsib%pack_size(1), i4))

        localsize = accel_kernel_workgroup_size(kernel_vpsi_spinors_complex)/(psib%pack_size(1)/2)

        dim3 = pnp/(accel_max_size_per_dim(2)*localsize) + 1
        dim2 = min(accel_max_size_per_dim(2)*localsize, pad(pnp, localsize))

        call accel_kernel_run(kernel_vpsi_spinors_complex, (/psib%pack_size(1)/2, dim2, dim3/), &
          (/psib%pack_size(1)/2, localsize, 1_i8/))

      end select
    end if

    call accel_finish()

    call profiling_count_operations((R_MUL*psib%nst_linear)*mesh%np)
    call profiling_count_transfers(mesh%np, M_ONE)
    call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

    if(.not. present(potential_accel)) then
      call accel_release_buffer(potential_accel_)
      SAFE_DEALLOCATE_P(potential_accel_)
    end if
    if(.not. present(impotential_accel) .and. present(impotential)) then
      call accel_release_buffer(impotential_accel_)
      SAFE_DEALLOCATE_P(impotential_accel_)
    end if

  case (BATCH_PACKED)

    select case (std%ispin)
    case (UNPOLARIZED, SPIN_POLARIZED)
      if (pot_is_cmplx) then
#ifdef R_TCOMPLEX
        !$omp parallel do private(vv, Imvv, ist)
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          Imvv = Impotential(ip, ispin)
          !$omp simd
          do ist = 1, psib%nst_linear
            vpsib%zff_pack(ist, ip) = vpsib%zff_pack(ist, ip) + (vv + M_zI*Imvv)*psib%zff_pack(ist, ip)
          end do
        end do
        call profiling_count_operations(2*((R_ADD+R_MUL)*psib%nst_linear)*mesh%np)
#else
        ! Complex potential can only be applied to complex batches
        ASSERT(.false.)
#endif
      else
        !$omp parallel do private(vv, ist)
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          !$omp simd
          do ist = 1, psib%nst_linear
            vpsib%X(ff_pack)(ist, ip) = vpsib%X(ff_pack)(ist, ip) + vv*psib%X(ff_pack)(ist, ip)
          end do
        end do
        call profiling_count_operations(((R_ADD+R_MUL)*psib%nst_linear)*mesh%np)
      end if
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst_linear, R_TOTYPE(M_ONE))

    case (SPINORS)
#ifdef R_TCOMPLEX
      ASSERT(mod(psib%nst_linear, 2) == 0)
      !the spinor case is more complicated since it mixes the two components.
      if (pot_is_cmplx) then
        !$omp parallel do private(psi1, psi2, ist, pot)
        do ip = 1, mesh%np
          do ist = 1, psib%nst_linear, 2
            psi1 = psib%zff_pack(ist    , ip)
            psi2 = psib%zff_pack(ist + 1, ip)
            pot(1:4) = potential(ip, 1:4) + M_zI * Impotential(ip, 1:4)
            vpsib%zff_pack(ist    , ip) = vpsib%zff_pack(ist    , ip) + &
              pot(1)*psi1 + (pot(3) + M_zI*pot(4))*psi2
            vpsib%zff_pack(ist + 1, ip) = vpsib%zff_pack(ist + 1, ip) + &
              pot(2)*psi2 + (pot(3) - M_zI*pot(4))*psi1
          end do
        end do
        !$omp end parallel do
        call profiling_count_operations((7*R_ADD + 7*R_MUL)*mesh%np*psib%nst)
      else
        !$omp parallel do private(psi1, psi2, ist)
        do ip = 1, mesh%np
          do ist = 1, psib%nst_linear, 2
            psi1 = psib%zff_pack(ist    , ip)
            psi2 = psib%zff_pack(ist + 1, ip)
            vpsib%zff_pack(ist    , ip) = vpsib%zff_pack(ist    , ip) + &
              potential(ip, 1)*psi1 + (potential(ip, 3) + M_zI*potential(ip, 4))*psi2
            vpsib%zff_pack(ist + 1, ip) = vpsib%zff_pack(ist + 1, ip) + &
              potential(ip, 2)*psi2 + (potential(ip, 3) - M_zI*potential(ip, 4))*psi1
          end do
        end do
        !$omp end parallel do
        call profiling_count_operations((6*R_ADD + 6*R_MUL)*mesh%np*psib%nst)
      end if
#else
      ! Spinors always imply complex batches
      ASSERT(.false.)
#endif
    end select

  case (BATCH_NOT_PACKED)

    select case (std%ispin)
    case (UNPOLARIZED, SPIN_POLARIZED)
      if (pot_is_cmplx) then
#ifdef R_TCOMPLEX
        !$omp parallel private(ip, ist)
        do ist = 1, psib%nst
          !$omp do simd
          do ip = 1, mesh%np
            vpsib%X(ff)(ip, 1, ist) = vpsib%X(ff)(ip, 1, ist) + &
              (potential(ip, ispin)+ M_zI*Impotential(ip, ispin)) * psib%X(ff)(ip, 1, ist)
          end do
          !$omp end do simd nowait
        end do
        !$omp end parallel
        call profiling_count_operations(2*((R_ADD+R_MUL)*psib%nst)*mesh%np)
#else
        ! Complex potential can only be applied to complex batches
        ASSERT(.false.)
#endif
      else
        !$omp parallel private(ip, ist)
        do ist = 1, psib%nst
          !$omp do simd
          do ip = 1, mesh%np
            vpsib%X(ff)(ip, 1, ist) = vpsib%X(ff)(ip, 1, ist) + &
              potential(ip, ispin) * psib%X(ff)(ip, 1, ist)
          end do
          !$omp end do simd nowait
        end do
        !$omp end parallel
        call profiling_count_operations(((R_ADD+R_MUL)*psib%nst)*mesh%np)
      end if

      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

    case (SPINORS)
#ifdef R_TCOMPLEX
      !the spinor case is more complicated since it mixes the two components.
      if (pot_is_cmplx) then
        do ist = 1, psib%nst
          psi  => psib%zff(:, :, ist)
          vpsi => vpsib%zff(:, :, ist)

          do ip = 1, mesh%np
            pot(1:4) = potential(ip, 1:4) + M_zI * Impotential(ip, 1:4)
            vpsi(ip, 1) = vpsi(ip, 1) + pot(1)*psi(ip, 1) + &
              (pot(3) + M_zI*pot(4))*psi(ip, 2)
            vpsi(ip, 2) = vpsi(ip, 2) + pot(2)*psi(ip, 2) + &
              (pot(3) - M_zI*pot(4))*psi(ip, 1)
          end do
        end do
        call profiling_count_operations((7*R_ADD + 7*R_MUL)*mesh%np*psib%nst)

      else
        do ist = 1, psib%nst
          psi  => psib%zff(:, :, ist)
          vpsi => vpsib%zff(:, :, ist)

          do ip = 1, mesh%np
            vpsi(ip, 1) = vpsi(ip, 1) + potential(ip, 1)*psi(ip, 1) + &
              (potential(ip, 3) + M_zI*potential(ip, 4))*psi(ip, 2)
            vpsi(ip, 2) = vpsi(ip, 2) + potential(ip, 2)*psi(ip, 2) + &
              (potential(ip, 3) - M_zI*potential(ip, 4))*psi(ip, 1)
          end do
        end do
        call profiling_count_operations((6*R_ADD + 6*R_MUL)*mesh%np*psib%nst)
      end if
#else
      ! Spinors always imply complex batches
      ASSERT(.false.)
#endif
    end select

  end select

  call profiling_out(prof_vlpsi)
  POP_SUB(X(hamiltonian_elec_base_local_sub))

end subroutine X(hamiltonian_elec_base_local_sub)

! -----------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_magnetic)(this, mesh, der, std, ep, ispin, psib, vpsib)
  type(hamiltonian_elec_base_t),  intent(in)    :: this
  class(mesh_t),                  intent(in)    :: mesh
  type(derivatives_t),            intent(in)    :: der
  type(states_elec_dim_t),        intent(in)    :: std
  type(epot_t),                   intent(in)    :: ep
  integer,                        intent(in)    :: ispin
  type(wfs_elec_t), target,       intent(inout) :: psib
  type(wfs_elec_t), target,       intent(inout) :: vpsib

  integer :: idir, ip
  class(batch_t), allocatable :: gradb(:)
  type(wfs_elec_t) :: adotpb
  FLOAT, allocatable :: a2(:,:)
  FLOAT, allocatable :: avec(:,:,:)

  if (.not. hamiltonian_elec_base_has_magnetic(this)) return

  if (.not. allocated(this%vector_potential)) return

#ifndef R_TCOMPLEX
  ! Vector potential not allowed with real wavefunctions
  ASSERT(.false.)
#endif

  call profiling_in(prof_magnetic, TOSTRING(X(MAGNETIC)))
  PUSH_SUB(X(hamiltonian_elec_base_magnetic))

  call psib%clone_to_array(gradb, mesh%box%dim)
  call psib%copy_to(adotpb, copy_data = .false.)

  call X(derivatives_batch_grad)(der, psib, gradb)

  SAFE_ALLOCATE(a2(1:mesh%np,1:std%nspin))
  SAFE_ALLOCATE(avec(1:mesh%np_part, 1:std%nspin, mesh%box%dim))
  a2 = M_ZERO
  avec = M_ZERO

  !$omp parallel do
  do ip = 1, mesh%np
    a2(ip,1:std%spin_channels) = (M_HALF / this%mass) * sum(this%vector_potential(1:mesh%box%dim, ip)**2)
    do idir = 1, mesh%box%dim
      avec(ip, 1:std%spin_channels, idir) =  (M_ONE / this%mass) * this%vector_potential(idir, ip)
    end do
  end do

  ! Adding the A^2 term
  call zhamiltonian_elec_base_local_sub(a2, mesh, std, ispin, psib, vpsib)

  ! Adding the A.p term
  call batch_set_zero(adotpb)
  do idir = 1, mesh%box%dim
    call zhamiltonian_elec_base_local_sub(avec(:,:,idir), mesh, std, ispin, gradb(idir), adotpb)
  end do
  call batch_axpy(mesh%np, M_zI, adotpb, vpsib)

  SAFE_DEALLOCATE_A(a2)
  SAFE_DEALLOCATE_A(avec)

  do idir = 1, mesh%box%dim
    call gradb(idir)%end()
  end do
  call adotpb%end()
  SAFE_DEALLOCATE_A(gradb)


  POP_SUB(X(hamiltonian_elec_base_magnetic))
  call profiling_out(prof_magnetic)
end subroutine X(hamiltonian_elec_base_magnetic)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_start)(this, mesh, std, spiral_bnd, psib, projection)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  class(mesh_t),                    intent(in)    :: mesh
  type(states_elec_dim_t),          intent(in)    :: std
  logical,                          intent(in)    :: spiral_bnd
  type(wfs_elec_t),                 intent(in)    :: psib
  type(projection_t),               intent(out)   :: projection

  integer :: ist, ip, iproj, imat, nreal, iprojection
  integer :: npoints, nprojs, nst_linear, maxnpoints
  integer, allocatable :: ind(:)
  type(projector_matrix_t), pointer :: pmat
  integer(i8) :: padnprojs, lnprojs, size, size_unfolded
  integer :: nphase
  type(profile_t), save :: cl_prof
  type(accel_kernel_t), save, target :: zker_proj_bra, dker_proj_bra, dker_proj_bra_phase_spiral
  type(accel_kernel_t), save, target :: zker_proj_bra_phase, dker_proj_bra_phase, zker_proj_bra_phase_spiral
  type(accel_kernel_t), pointer :: kernel
  integer, allocatable :: spin_to_phase(:)
  R_TYPE, allocatable :: lpsi(:, :)
#ifdef R_TCOMPLEX
  integer :: iphase, map_ip
  CMPLX, allocatable :: tmp_proj(:, :)
#endif

  integer :: block_size

  if (.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)

  call profiling_in(prof_vnlpsi_start, TOSTRING(X(VNLPSI_MAT_BRA)))
  PUSH_SUB(X(hamiltonian_elec_base_nlocal_start))

  nst_linear = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst_linear
#else
  nreal = nst_linear
#endif
  nphase = 1

  if (psib%has_phase) then
    ASSERT(allocated(this%projector_phases))
  end if

  if (psib%status() == BATCH_DEVICE_PACKED) then

    call profiling_in(cl_prof, TOSTRING(X(CL_PROJ_BRA)))
    ! only do this if we have some points of projector matrices
    if (this%max_npoints > 0) then

      call accel_create_buffer(projection%buff_projection, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, &
        this%full_projection_size*psib%pack_size(1))

      if (allocated(this%projector_phases)) then
        ASSERT(R_TYPE_VAL == TYPE_CMPLX)

        if (spiral_bnd) then

          nphase = 3

          SAFE_ALLOCATE(spin_to_phase(1:psib%pack_size(1)))
          call accel_create_buffer(projection%buff_spin_to_phase, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, psib%pack_size(1))

          do ist = 1, nst_linear
            if (this%spin(3, psib%linear_to_ist(ist), psib%ik) > 0 .and. psib%linear_to_idim(ist) == 2) then
              spin_to_phase(ist) = 1
            else if (this%spin(3, psib%linear_to_ist(ist), psib%ik) < 0 .and. psib%linear_to_idim(ist) == 1) then
              spin_to_phase(ist) = 2
            else
              spin_to_phase(ist) = 0
            end if
          end do
          ! This might not be necessary:
          do ist = nst_linear+1, int(psib%pack_size(1), i4)
            spin_to_phase(ist) = 0
          end do

          call accel_write_buffer(projection%buff_spin_to_phase, psib%pack_size(1), spin_to_phase)

          if (this%projector_matrices(1)%is_cmplx) then
            call accel_kernel_start_call(zker_proj_bra_phase_spiral, 'projector.cl', 'zprojector_bra_phase_spiral',&
              flags = '-DRTYPE_COMPLEX')
            kernel => zker_proj_bra_phase_spiral
          else
            call accel_kernel_start_call(dker_proj_bra_phase_spiral, 'projector.cl', 'dprojector_bra_phase_spiral',&
              flags = '-DRTYPE_DOUBLE')
            kernel => dker_proj_bra_phase_spiral
          end if
          SAFE_DEALLOCATE_A(spin_to_phase)
        else
          if (this%projector_matrices(1)%is_cmplx) then
            call accel_kernel_start_call(zker_proj_bra_phase, 'projector.cl', 'zprojector_bra_phase',&
              flags = '-DRTYPE_COMPLEX')
            kernel => zker_proj_bra_phase
          else
            call accel_kernel_start_call(dker_proj_bra_phase, 'projector.cl', 'dprojector_bra_phase',&
              flags = '-DRTYPE_DOUBLE')
            kernel => dker_proj_bra_phase
          end if
        end if
        size = psib%pack_size(1)

      else
        if (this%projector_matrices(1)%is_cmplx) then
          call accel_kernel_start_call(zker_proj_bra, 'projector.cl', 'zprojector_bra', flags = '-DRTYPE_COMPLEX')
          kernel => zker_proj_bra
          size = psib%pack_size(1)
        else
          call accel_kernel_start_call(dker_proj_bra, 'projector.cl', 'dprojector_bra', flags = '-DRTYPE_DOUBLE')
          kernel => dker_proj_bra
          size = psib%pack_size_real(1)
        end if
      end if

      ! NOTE: For the projector_bra_* kernels, there is no problem with  self-overlapping spheres.

      call accel_set_kernel_arg(kernel, 0, this%nprojector_matrices)
      call accel_set_kernel_arg(kernel, 1, this%buff_offsets)
      call accel_set_kernel_arg(kernel, 2, this%buff_matrices)
      call accel_set_kernel_arg(kernel, 3, this%buff_maps)
      call accel_set_kernel_arg(kernel, 4, this%buff_scals)
      call accel_set_kernel_arg(kernel, 5, psib%ff_device)
      call accel_set_kernel_arg(kernel, 6, int(log2(size), i4))
      call accel_set_kernel_arg(kernel, 7, projection%buff_projection)
      call accel_set_kernel_arg(kernel, 8, int(log2(size), i4))

      if (allocated(this%projector_phases)) then
        call accel_set_kernel_arg(kernel, 9, this%buff_projector_phases)
        ! Note: we need to use this%nphase, as the kernel might be called with spiral=false, but
        !       the phases been built with spiralBC=true
        call accel_set_kernel_arg(kernel, 10, (psib%ik - std%kpt%start)*this%total_points*this%nphase)
        if (spiral_bnd) then
          call accel_set_kernel_arg(kernel, 11, projection%buff_spin_to_phase)
          call accel_set_kernel_arg(kernel, 12, this%nphase)
        end if
      end if

      ! In case of CUDA we use an optimized kernel, in which the loop over npoints is broken
      ! further into chunks, in order to parallelize over the threads within a warp.
      ! Therefore we need to launch warp_size * size kernels. The size of each block needs to
      ! have multiples of warp_size as x-dimension.

      size_unfolded = size * accel%warp_size
      padnprojs = pad_pow2(this%max_nprojs)
      lnprojs = min(accel_kernel_workgroup_size(kernel)/accel%warp_size, padnprojs)

      call accel_kernel_run(kernel, &
        (/size_unfolded, padnprojs, int(this%nprojector_matrices, i8)/), (/int(accel%warp_size, i8), lnprojs, 1_i8/))

      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)

        npoints = pmat%npoints
        nprojs = pmat%nprojs

        !! update number of operations for nphase !!
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst_linear*nprojs)
      end do

      call accel_finish()
    end if

    if (mesh%parallel_in_domains) then
      SAFE_ALLOCATE(projection%X(projection)(1:psib%pack_size(1), 1:this%full_projection_size))
      projection%X(projection) = M_ZERO
      if (this%max_npoints > 0) then
        call accel_read_buffer(projection%buff_projection, &
          this%full_projection_size*psib%pack_size(1), projection%X(projection))
      end if
    end if
    call profiling_out(cl_prof)

  else

    ! This routine uses blocking to optimize cache usage. One block of
    ! |phi> is loaded in cache L1 and then then we calculate the dot
    ! product of it with the corresponding blocks of |psi_k>, next we
    ! load another block and do the same. This way we only have to load
    ! |psi> from the L2 or memory.
    block_size = hardware%X(block_size)


    SAFE_ALLOCATE(projection%X(projection)(1:nst_linear, 1:this%full_projection_size))
    projection%X(projection) = M_ZERO

    SAFE_ALLOCATE(ind(1:this%nprojector_matrices))

    iprojection = 0
    maxnpoints = 0
    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)
      npoints = pmat%npoints
      maxnpoints = max(maxnpoints, npoints)
      nprojs = pmat%nprojs
      ind(imat) = iprojection
      iprojection = iprojection + nprojs
      call profiling_count_operations(nprojs*(R_ADD + R_MUL)*npoints + nst_linear*nprojs)
      if (allocated(this%projector_phases)) then
        call profiling_count_operations(R_MUL*npoints*nst_linear)
      end if
    end do

    SAFE_ALLOCATE(lpsi(1:nst_linear, 1:maxnpoints))

    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)
      iprojection = ind(imat)
      npoints = pmat%npoints
      nprojs = pmat%nprojs

      if (npoints == 0) cycle

      if (.not. allocated(this%projector_phases)) then
        if (psib%status() == BATCH_PACKED) then

          !$omp parallel do private(ist, ip)
          do ip = 1, npoints
            !$omp simd
            do ist= 1, nst_linear
              lpsi(ist, ip) = psib%X(ff_pack)(ist, pmat%map(ip))
            end do
          end do

        else

          do ist = 1, nst_linear
            !$omp parallel do
            do ip = 1, npoints
              lpsi(ist, ip) = psib%X(ff_linear)(pmat%map(ip), ist)
            end do
          end do

        end if

      else
        if (.not. spiral_bnd) then
          if (psib%status() == BATCH_PACKED) then
#ifdef R_TCOMPLEX
            !$omp parallel do private(ist, map_ip)
            do ip = 1, npoints
              map_ip = pmat%map(ip)
              !$omp simd
              do ist = 1, nst_linear
                lpsi(ist, ip) = psib%zff_pack(ist, map_ip)*this%projector_phases(ip, 1, imat, psib%ik)
              end do
            end do
#else
            ! Phases not allowed for real batches
            ASSERT(.false.)
#endif
          else
#ifdef R_TCOMPLEX
            do ist = 1, nst_linear
              !$omp parallel do
              do ip = 1, npoints
                lpsi(ist, ip) = psib%zff_linear(pmat%map(ip), ist)*this%projector_phases(ip, 1, imat, psib%ik)
              end do
            end do
#else
            ! Phases not allowed for real batches
            ASSERT(.false.)
#endif
          end if
        else
          if (psib%status() == BATCH_PACKED) then
#ifdef R_TCOMPLEX
            !$omp parallel do private(ist)
            do ip = 1, npoints
              do ist = 1, nst_linear, 2
                if (this%spin(3,psib%linear_to_ist(ist), psib%ik)>0) then
                  lpsi(ist, ip)   = psib%zff_pack(ist,   pmat%map(ip))*this%projector_phases(ip, 1, imat, psib%ik)
                  lpsi(ist+1, ip) = psib%zff_pack(ist+1, pmat%map(ip))*this%projector_phases(ip, 2, imat, psib%ik)
                else
                  lpsi(ist, ip)   = psib%zff_pack(ist,   pmat%map(ip))*this%projector_phases(ip, 3, imat, psib%ik)
                  lpsi(ist+1, ip) = psib%zff_pack(ist+1, pmat%map(ip))*this%projector_phases(ip, 1, imat, psib%ik)
                end if
              end do
            end do
#else
            ! Phases not allowed for real batches
            ASSERT(.false.)
#endif
          else
#ifdef R_TCOMPLEX
            do ist = 1, nst_linear
              if (this%spin(3, psib%linear_to_ist(ist), psib%ik) > 0 .and. psib%linear_to_idim(ist) == 2) then
                iphase = 2
              else if (this%spin(3, psib%linear_to_ist(ist), psib%ik) < 0 .and. psib%linear_to_idim(ist) == 1) then
                iphase = 3
              else
                iphase = 1
              end if
              !$omp parallel do
              do ip = 1, npoints
                lpsi(ist, ip) = psib%zff_linear(pmat%map(ip), ist)*this%projector_phases(ip, iphase, imat, psib%ik)
              end do
            end do
#else
            ! Phases not allowed for real batches
            ASSERT(.false.)
#endif
          end if

        end if
      end if

      if (pmat%is_cmplx) then
#ifdef R_TCOMPLEX
        SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:nst_linear))
        call blas_gemm('C', 'T', nprojs, nst_linear, npoints, &
          M_z1, pmat%zprojectors(1, 1), npoints, lpsi(1, 1), nst_linear, M_z0, tmp_proj(1,1), nprojs)
        !$omp parallel do private(iproj, ist)
        do iproj = 1, nprojs
          do ist = 1, nst_linear
            projection%X(projection)(ist, iprojection + iproj) = tmp_proj(iproj, ist)
          end do
        end do
        SAFE_DEALLOCATE_A(tmp_proj)
        call profiling_count_operations(nst_linear*nprojs*M_TWO*npoints)
#else
        ! Complex projection matrix not allowed for real batches
        ASSERT(.false.)
#endif
      else
        call blas_gemm('N', 'N', nreal, nprojs, npoints, &
          M_ONE, lpsi(1, 1), nreal, pmat%dprojectors(1, 1), npoints, M_ZERO,  projection%X(projection)(1, iprojection + 1), nreal)
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
      end if

      !$omp parallel do private(iproj, ist)
      do iproj = 1, nprojs
        do ist = 1, nst_linear
          projection%X(projection)(ist, iprojection + iproj) = projection%X(projection)(ist, iprojection + iproj)*pmat%scal(iproj)
        end do
      end do

    end do

    SAFE_DEALLOCATE_A(ind)
    SAFE_DEALLOCATE_A(lpsi)
  endif

  POP_SUB(X(hamiltonian_elec_base_nlocal_start))
  call profiling_out(prof_vnlpsi_start)
end subroutine X(hamiltonian_elec_base_nlocal_start)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_finish)(this, mesh, spiral_bnd, std, projection, vpsib)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  class(mesh_t),                    intent(in)    :: mesh
  type(states_elec_dim_t),          intent(in)    :: std
  logical,                          intent(in)    :: spiral_bnd
  type(projection_t),       target, intent(inout) :: projection
  class(wfs_elec_t),                intent(inout) :: vpsib

  integer :: ist, ip, imat, nreal, iprojection
  integer :: npoints, nprojs, nst_linear, nphase
  R_TYPE, allocatable :: psi(:, :)
  type(projector_matrix_t), pointer :: pmat
  type(profile_t), save :: reduce_prof
#ifdef R_TCOMPLEX
  integer :: iproj, idim, iphase
  CMPLX  :: phase, phase_pq, phase_mq
  CMPLX, allocatable :: tmp_proj(:, :, :)
#endif

  if (.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)

  call profiling_in(prof_vnlpsi_finish, TOSTRING(X(VNLPSI_MAT_KET)))
  PUSH_SUB(X(hamiltonian_elec_base_nlocal_finish))

  nst_linear = vpsib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst_linear
#else
  nreal = nst_linear
#endif
  nphase = 1
  if (spiral_bnd) nphase = 3

  if (vpsib%has_phase) then
    ASSERT(allocated(this%projector_phases))
  end if

  ! reduce the projections
  if (mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(VNLPSI_MAT_REDUCE)))
    call mesh%allreduce(projection%X(projection))
    call profiling_out(reduce_prof)
  end if

  if (vpsib%status() == BATCH_DEVICE_PACKED) then

    if (mesh%parallel_in_domains) then
      ! only do this if we have points of some projector matrices
      if (this%max_npoints > 0) then
        call accel_write_buffer(projection%buff_projection, &
          this%full_projection_size*vpsib%pack_size(1), projection%X(projection))
      end if
      SAFE_DEALLOCATE_A(projection%X(projection))
    end if

    call finish_accel()
    call accel_release_buffer(projection%buff_projection)
    if (spiral_bnd) then
      call accel_release_buffer(projection%buff_spin_to_phase)
    end if

  else

    ASSERT(allocated(projection%X(projection)))

    iprojection = 0
    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)

      npoints = pmat%npoints
      nprojs = pmat%nprojs

      if (allocated(pmat%zmix)) then
#ifdef R_TCOMPLEX
        SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:vpsib%nst, 1:std%dim))

        do ist = 1, vpsib%nst

          tmp_proj(1:nprojs, ist, 1) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 1), &
            projection%X(projection)((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs)) &
            + matmul(pmat%zmix(1:nprojs, 1:nprojs, 3), &
            projection%X(projection)((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs))

          tmp_proj(1:nprojs, ist, 2) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 2), &
            projection%X(projection)((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs)) &
            + matmul(pmat%zmix(1:nprojs, 1:nprojs, 4), &
            projection%X(projection)((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs))
        end do

        do ist = 1, vpsib%nst
          do idim = 1, std%dim
            do iproj = 1, nprojs
              projection%X(projection)((ist-1)*std%dim+idim, iprojection + iproj) = tmp_proj(iproj, ist, idim)
            end do
          end do
        end do

        SAFE_DEALLOCATE_A(tmp_proj)
#else
        ! Complex projection matrix not allowed for real batches
        ASSERT(.false.)
#endif
      else if (allocated(pmat%dmix)) then

        if (allocated(pmat%dmix)) then
          do ist = 1, nst_linear
            projection%X(projection)(ist, iprojection + 1:iprojection + nprojs) = &
              matmul(pmat%dmix(1:nprojs, 1:nprojs), projection%X(projection)(ist, iprojection + 1:iprojection + nprojs))
          end do
        end if

      end if

      if (npoints /=  0) then

        SAFE_ALLOCATE(psi(1:nst_linear, 1:npoints))

        ! Matrix-multiply again.
        ! the line below does: psi = matmul(projection, transpose(pmat%projectors))
        if (.not. pmat%is_cmplx) then
          call blas_gemm('N', 'T', nreal, npoints, nprojs, &
            M_ONE, projection%X(projection)(1, iprojection + 1), nreal, pmat%dprojectors(1, 1), npoints, &
            M_ZERO, psi(1, 1), nreal)
          call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
        else
#ifdef R_TCOMPLEX
          call blas_gemm('N', 'T', nst_linear, npoints, nprojs, &
            M_z1, projection%X(projection)(1, iprojection + 1), nst_linear, pmat%zprojectors(1, 1), npoints, &
            M_z0, psi(1, 1), nst_linear)
          call profiling_count_operations(nst_linear*nprojs*(R_ADD+R_MUL)*npoints)
#else
          ! Complex projection matrix not allowed for real batches
          ASSERT(.false.)
#endif
        end if

        call profiling_in(prof_scatter, TOSTRING(X(PROJ_MAT_SCATTER)))

        if (.not. allocated(this%projector_phases)) then
          ! and copy the points from the local buffer to its position
          if (vpsib%status() == BATCH_PACKED) then
            !$omp parallel do private(ip, ist) if (.not. this%projector_self_overlap)
            do ip = 1, npoints
              do ist = 1, nst_linear
                vpsib%X(ff_pack)(ist, pmat%map(ip)) = vpsib%X(ff_pack)(ist, pmat%map(ip)) + psi(ist, ip)
              end do
            end do
            !$omp end parallel do
          else
            do ist = 1, nst_linear
              !$omp parallel do if (.not. this%projector_self_overlap)
              do ip = 1, npoints
                vpsib%X(ff_linear)(pmat%map(ip), ist) = vpsib%X(ff_linear)(pmat%map(ip), ist) + psi(ist, ip)
              end do
              !$omp end parallel do
            end do
          end if
          call profiling_count_operations(nst_linear*npoints*R_ADD)
        else
#ifdef R_TCOMPLEX
          if (.not. spiral_bnd) then
            ! and copy the points from the local buffer to its position
            if (vpsib%status() == BATCH_PACKED) then
              !$omp parallel do private(ip, ist, phase) if (.not. this%projector_self_overlap)
              do ip = 1, npoints
                phase = conjg(this%projector_phases(ip, 1, imat, vpsib%ik))
                do ist = 1, nst_linear
                  vpsib%zff_pack(ist, pmat%map(ip)) = vpsib%zff_pack(ist, pmat%map(ip)) &
                    + psi(ist, ip)*phase
                end do
              end do
              !$omp end parallel do
            else
              do ist = 1, nst_linear
                !$omp parallel do if (.not. this%projector_self_overlap)
                do ip = 1, npoints
                  vpsib%zff_linear(pmat%map(ip), ist) = vpsib%zff_linear(pmat%map(ip), ist) &
                    + psi(ist, ip)*conjg(this%projector_phases(ip, 1, imat, vpsib%ik))
                end do
                !$omp end parallel do
              end do
            end if
            call profiling_count_operations(nst_linear*npoints*(R_ADD+R_MUL))
          else
            ! and copy the points from the local buffer to its position
            if (vpsib%status() == BATCH_PACKED) then
              !$omp parallel do private(ip, ist, phase, phase_pq, phase_mq) if (.not. this%projector_self_overlap)
              do ip = 1, npoints
                phase = conjg(this%projector_phases(ip, 1, imat, vpsib%ik))
                phase_pq = conjg(this%projector_phases(ip, 2, imat, vpsib%ik))
                phase_mq = conjg(this%projector_phases(ip, 3, imat, vpsib%ik))
                do ist = 1, nst_linear, 2
                  if (this%spin(3, vpsib%linear_to_ist(ist), vpsib%ik) > 0) then
                    vpsib%zff_pack(ist, pmat%map(ip)) = vpsib%zff_pack(ist, pmat%map(ip)) &
                      + psi(ist, ip)*phase
                    vpsib%zff_pack(ist+1, pmat%map(ip)) = vpsib%zff_pack(ist+1, pmat%map(ip)) &
                      + psi(ist+1, ip)*phase_pq
                  else
                    vpsib%zff_pack(ist, pmat%map(ip)) = vpsib%zff_pack(ist, pmat%map(ip)) &
                      + psi(ist, ip)*phase_mq
                    vpsib%zff_pack(ist+1, pmat%map(ip)) = vpsib%zff_pack(ist+1, pmat%map(ip)) &
                      + psi(ist+1, ip)*phase
                  end if
                end do
              end do
              !$omp end parallel do
            else
              do ist = 1, nst_linear
                if (this%spin(3, vpsib%linear_to_ist(ist), vpsib%ik) > 0 .and. vpsib%linear_to_idim(ist) == 2) then
                  iphase = 2
                else if (this%spin(3, vpsib%linear_to_ist(ist), vpsib%ik) < 0 .and. vpsib%linear_to_idim(ist) == 1) then
                  iphase = 3
                else
                  iphase = 1
                end if
                !$omp parallel do if (.not. this%projector_self_overlap)
                do ip = 1, npoints
                  vpsib%zff_linear(pmat%map(ip), ist) = vpsib%zff_linear(pmat%map(ip), ist) &
                    + psi(ist, ip)*conjg(this%projector_phases(ip, iphase, imat, vpsib%ik))
                end do
                !$omp end parallel do
              end do
            end if
            call profiling_count_operations(nst_linear*npoints*(R_ADD+R_MUL))
          end if
#else
          ! Phases not allowed for real batches
          ASSERT(.false.)
#endif
        end if
        call profiling_out(prof_scatter)
      end if

      SAFE_DEALLOCATE_A(psi)

      iprojection = iprojection + nprojs
    end do
  end if
  SAFE_DEALLOCATE_A(projection%X(projection))

  POP_SUB(X(hamiltonian_elec_base_nlocal_finish))
  call profiling_out(prof_vnlpsi_finish)

contains

  subroutine finish_accel()

    integer :: imat, iregion, nregions_self_overlap, iregion_self_overlap
    integer(i8) :: padnprojs, lnprojs, size, size2, wgsize
    type(profile_t), save :: cl_prof
    type(accel_kernel_t), save, target :: ker_proj_ket, dker_mix, zker_mix
    type(accel_kernel_t), save, target :: zker_proj_ket_phase, dker_proj_ket_phase
    type(accel_kernel_t), save, target :: dker_proj_ket_phase_spiral, zker_proj_ket_phase_spiral
    type(accel_kernel_t), pointer :: kernel, ker_mix
    type(accel_mem_t), pointer :: buff_proj

    PUSH_SUB(X(hamiltonian_elec_base_nlocal_finish).finish_accel)

    ! In this case we run one kernel per projector, since all write to
    ! the wave-function. Otherwise we would need to do atomic
    ! operations.

    ! only do this if we have points of some projector matrices
    if (this%max_npoints > 0) then

      if (this%projector_mix) then
        call profiling_in(cl_prof, TOSTRING(X(CL_PROJ_MIX)))

        SAFE_ALLOCATE(buff_proj)
        if (allocated(this%projector_matrices(1)%zmix)) then
          size = vpsib%pack_size(1)
          size2 = size/2
          call accel_kernel_start_call(zker_mix, 'projector.cl', 'zprojector_mix', flags = '-DRTYPE_COMPLEX')
          call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, this%full_projection_size*size)
          ker_mix => zker_mix
        else
          size = vpsib%pack_size_real(1)
          size2 = size
          call accel_kernel_start_call(dker_mix, 'projector.cl', 'dprojector_mix', flags = '-DRTYPE_DOUBLE')
          call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, this%full_projection_size*size)
          ker_mix => dker_mix
        end if

        call accel_set_kernel_arg(ker_mix, 0, this%nprojector_matrices)
        call accel_set_kernel_arg(ker_mix, 1, this%buff_offsets)
        call accel_set_kernel_arg(ker_mix, 2, this%buff_mix)
        call accel_set_kernel_arg(ker_mix, 3, projection%buff_projection)
        call accel_set_kernel_arg(ker_mix, 4, int(log2(size), i4))
        call accel_set_kernel_arg(ker_mix, 5, buff_proj)

        padnprojs = pad_pow2(this%max_nprojs)
        lnprojs = min(accel_kernel_workgroup_size(ker_mix)/size2, padnprojs)

        call accel_kernel_run(ker_mix, &
          (/size2, padnprojs, int(this%nprojector_matrices, i8)/), (/size2, lnprojs, 1_i8/))

        call accel_finish()

        call profiling_out(cl_prof)
      else

        buff_proj => projection%buff_projection

      end if

      call profiling_in(cl_prof, TOSTRING(X(CL_PROJ_KET)))

      if (allocated(this%projector_phases)) then
        ASSERT(R_TYPE_VAL == TYPE_CMPLX)
        if (spiral_bnd) then
          if (this%projector_matrices(1)%is_cmplx) then
            call accel_kernel_start_call(zker_proj_ket_phase_spiral, 'projector.cl', 'zprojector_ket_phase_spiral',&
              flags = '-DRTYPE_COMPLEX')
            kernel => zker_proj_ket_phase_spiral
          else
            call accel_kernel_start_call(dker_proj_ket_phase_spiral, 'projector.cl', 'dprojector_ket_phase_spiral',&
              flags = '-DRTYPE_DOUBLE')
            kernel => dker_proj_ket_phase_spiral
          end if
        else
          if (this%projector_matrices(1)%is_cmplx) then
            call accel_kernel_start_call(zker_proj_ket_phase, 'projector.cl', 'zprojector_ket_phase',&
              flags = '-DRTYPE_COMPLEX')
            kernel => zker_proj_ket_phase
          else
            call accel_kernel_start_call(dker_proj_ket_phase, 'projector.cl', 'dprojector_ket_phase',&
              flags = '-DRTYPE_DOUBLE')
            kernel => dker_proj_ket_phase
          end if
        end if
        size = vpsib%pack_size(1)
      else
        if (this%projector_matrices(1)%is_cmplx) then
          call accel_kernel_start_call(ker_proj_ket, 'projector.cl', 'zprojector_ket', flags = '-DRTYPE_COMPLEX')
          size = vpsib%pack_size(1)
        else
          call accel_kernel_start_call(ker_proj_ket, 'projector.cl', 'dprojector_ket', flags = '-DRTYPE_DOUBLE')
          size = vpsib%pack_size_real(1)
        end if
        kernel => ker_proj_ket
      end if

      do iregion = 1, this%nregions

        do imat = this%regions(iregion), this%regions(iregion+1)-1
          nregions_self_overlap = this%projector_matrices(imat)%nregions

          do iregion_self_overlap = 1, nregions_self_overlap

            call accel_set_kernel_arg(kernel, 0, this%nprojector_matrices)
            call accel_set_kernel_arg(kernel, 1, imat - 1)
            call accel_set_kernel_arg(kernel, 2, this%buff_offsets)
            call accel_set_kernel_arg(kernel, 3, this%buff_matrices)
            call accel_set_kernel_arg(kernel, 4, this%buff_maps)
            call accel_set_kernel_arg(kernel, 5, buff_proj)
            call accel_set_kernel_arg(kernel, 6, int(log2(size), i4))
            call accel_set_kernel_arg(kernel, 7, vpsib%ff_device)
            call accel_set_kernel_arg(kernel, 8, int(log2(size), i4))

            call accel_set_kernel_arg(kernel, 9, this%projector_matrices(imat)%regions(iregion_self_overlap) - 1)
            call accel_set_kernel_arg(kernel,10, this%projector_matrices(imat)%regions(iregion_self_overlap+1) - 1)

            if (allocated(this%projector_phases)) then
              call accel_set_kernel_arg(kernel, 11, this%buff_projector_phases)
              ! Note: we need to use this%nphase, as the kernel might be called with spiral=false, but
              !       the phases been built with spiralBC=true
              call accel_set_kernel_arg(kernel, 12, (vpsib%ik - std%kpt%start)*this%total_points*this%nphase)
              if (spiral_bnd) then
                call accel_set_kernel_arg(kernel, 13, projection%buff_spin_to_phase)
                call accel_set_kernel_arg(kernel, 14, this%nphase)
              end if
            end if

            wgsize = accel_kernel_workgroup_size(kernel)/size

            ! implicit loops over ist, ip
            call accel_kernel_run(kernel, &
              (/size, int(pad(this%max_npoints, wgsize), i8), 1_i8/), &
              (/size, wgsize, 1_i8/))

            call accel_finish()
          end do ! iregion_self_overlap
        end do ! iatom

      end do ! iregion

      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)
        npoints = pmat%npoints
        nprojs = pmat%nprojs
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
        call profiling_count_operations(nst_linear*npoints*R_ADD)
      end do

      if (this%projector_mix) then
        call accel_release_buffer(buff_proj)
        SAFE_DEALLOCATE_P(buff_proj)
      end if
      call profiling_out(cl_prof)
    end if

    POP_SUB(X(hamiltonian_elec_base_nlocal_finish).finish_accel)
  end subroutine finish_accel

end subroutine X(hamiltonian_elec_base_nlocal_finish)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_force)(this, mesh, st, spiral_bnd, iqn, ndim, psi1b, psi2b, force)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  class(mesh_t),                         intent(in)    :: mesh
  type(states_elec_t),                   intent(in)    :: st
  logical,                               intent(in)    :: spiral_bnd
  integer,                               intent(in)    :: iqn
  integer,                               intent(in)    :: ndim
  type(wfs_elec_t),                      intent(in)    :: psi1b
  type(wfs_elec_t),                      intent(in)    :: psi2b(:)
  FLOAT,                                 intent(inout) :: force(:, :)

  integer :: ii, ist, ip, iproj, imat, nreal, iprojection, iatom, idir
  integer :: npoints, nprojs, nst_linear
  FLOAT, allocatable :: ff(:)
  R_TYPE, allocatable :: psi(:, :, :), projs(:, :, :)
  type(projector_matrix_t), pointer :: pmat
  integer(i8) :: padnprojs, lnprojs, size, size_unfolded
  type(accel_mem_t) :: buff_projs
  type(profile_t), save :: cl_prof
  type(accel_kernel_t), save, target :: dker_proj_bra_force, zker_proj_bra_force
  type(accel_kernel_t), save, target :: zker_proj_bra_force_phase, dker_proj_bra_force_phase
  type(accel_kernel_t), pointer :: kernel
#ifdef R_TCOMPLEX
  integer :: idim
  CMPLX, allocatable :: tmp_proj(:, :, :)
#endif

  if (.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)

  call profiling_in(prof_matelement, TOSTRING(X(VNLPSI_MAT_ELEM)))
  PUSH_SUB(X(hamiltonian_elec_base_nlocal_force))

  ASSERT(psi1b%nst_linear == psi2b(1)%nst_linear)
  ASSERT(psi1b%status() == psi2b(1)%status())

  ASSERT(.not. spiral_bnd)

  nst_linear = psi1b%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst_linear
#else
  nreal = nst_linear
#endif

  if( .not. psi1b%status() == BATCH_DEVICE_PACKED) then
    SAFE_ALLOCATE(projs(0:ndim, 1:nst_linear, 1:this%full_projection_size))
    projs = CNST(0.0)

    iprojection = 0
    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)

      npoints = pmat%npoints
      nprojs = pmat%nprojs

      if (npoints /= 0) then

        SAFE_ALLOCATE(psi(0:ndim, 1:nst_linear, 1:npoints))

        call profiling_in(prof_matelement_gather, TOSTRING(X(PROJ_MAT_ELEM_GATHER)))

        ! collect all the points we need in a continuous array
        if (psi1b%status() == BATCH_PACKED) then
          do ip = 1, npoints
            do ist = 1, nst_linear
              psi(0, ist, ip) = psi1b%X(ff_pack)(ist, pmat%map(ip))
              do idir = 1, ndim
                psi(idir, ist, ip) = psi2b(idir)%X(ff_pack)(ist, pmat%map(ip))
              end do
            end do
          end do
        else
          do ip = 1, npoints
            do ist = 1, nst_linear
              psi(0, ist, ip) = psi1b%X(ff_linear)(pmat%map(ip), ist)
              do idir = 1, ndim
                psi(idir, ist, ip) = psi2b(idir)%X(ff_linear)(pmat%map(ip), ist)
              end do
            end do
          end do
        end if

        if (allocated(this%projector_phases)) then
#ifdef R_TCOMPLEX
          do ip = 1, npoints
            do ist = 1, nst_linear
              do idir = 0, ndim
                psi(idir, ist, ip) = this%projector_phases(ip, 1, imat, psi1b%ik)*psi(idir, ist, ip)
              end do
            end do
          end do
#else
          ! Phases not allowed for real batches
          ASSERT(.false.)
#endif
        end if

        call profiling_out(prof_matelement_gather)

        ! Now matrix-multiply to calculate the projections. We can do all the matrix multiplications at once
        if (.not. pmat%is_cmplx) then
          call blas_gemm('N', 'N', (ndim + 1)*nreal, nprojs, npoints, M_ONE, &
            psi(0, 1, 1), (ndim + 1)*nreal, pmat%dprojectors(1, 1), npoints, &
            M_ZERO, projs(0, 1, iprojection + 1), (ndim + 1)*nreal)

          call profiling_count_operations(nreal*(ndim + 1)*nprojs*M_TWO*npoints)
        else
#ifdef R_TCOMPLEX
          SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:nst_linear*(ndim + 1), 1))
          call blas_gemm('C', 'T', nprojs, (ndim + 1)*nst_linear, npoints, &
            M_z1, pmat%zprojectors(1, 1), npoints, psi(0, 1, 1), (ndim + 1)*nst_linear, &
            M_z0, tmp_proj(1,1,1), nprojs)
          do iproj = 1, nprojs
            do ist = 1, nst_linear
              do idir = 0, ndim
                projs(idir , ist, iprojection + iproj) = tmp_proj(iproj, (ist-1)*(ndim+1)+idir+1, 1)
              end do
            end do
          end do
          SAFE_DEALLOCATE_A(tmp_proj)

          call profiling_count_operations(nst_linear*(ndim + 1)*nprojs*(R_ADD+R_MUL)*npoints)
#endif
        end if

      else

        projs(0:ndim, 1:nst_linear, iprojection + 1:iprojection + nprojs) = CNST(0.0)

      end if

      SAFE_DEALLOCATE_A(psi)

      iprojection = iprojection + nprojs

    end do


  else
    ! Mostly Cargo Cult programming, copying accelerated projector application
    ! from hamiltonian_elec_base_nlocal_start, but modified to deal with
    ! the gradient of psi, psi2b, having ndim directions
    if (allocated(this%projector_phases)) then
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_proj_bra_force_phase, 'projector.cl', 'zprojector_bra_force_phase',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_proj_bra_force_phase
      else
        call accel_kernel_start_call(dker_proj_bra_force_phase, 'projector.cl', 'dprojector_bra_force_phase',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_proj_bra_force_phase
      end if
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
      size = psi1b%pack_size(1)
    else
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_proj_bra_force, 'projector.cl', 'zprojector_bra_force',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_proj_bra_force
        size = psi1b%pack_size(1)
      else
        call accel_kernel_start_call(dker_proj_bra_force, 'projector.cl', 'dprojector_bra_force',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_proj_bra_force
        size = psi1b%pack_size_real(1)
      end if
    end if

    SAFE_ALLOCATE(projs(0:ndim, 1:psi1b%pack_size(1), 1:this%full_projection_size))
    projs = CNST(0.0)

    if (this%max_npoints > 0) then
      call accel_create_buffer(buff_projs, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, &
        (ndim+1)*this%full_projection_size*psi1b%pack_size(1))
      call profiling_in(cl_prof, TOSTRING(X(CL_PROJ_FORCE)))

      do idir=0,ndim
        call accel_set_kernel_arg(kernel, 0, idir)
        call accel_set_kernel_arg(kernel, 1, ndim+1)
        call accel_set_kernel_arg(kernel, 2, int(size/psi1b%pack_size(1), i4))
        call accel_set_kernel_arg(kernel, 3, this%nprojector_matrices)
        call accel_set_kernel_arg(kernel, 4, this%buff_offsets)
        call accel_set_kernel_arg(kernel, 5, this%buff_matrices)
        call accel_set_kernel_arg(kernel, 6, this%buff_maps)
        if ( idir == 0 ) then
          call accel_set_kernel_arg(kernel, 7, psi1b%ff_device)
        else
          call accel_set_kernel_arg(kernel, 7, psi2b(idir)%ff_device)
        end if
        call accel_set_kernel_arg(kernel, 8, int(log2(size), i4))
        call accel_set_kernel_arg(kernel, 9, buff_projs)
        call accel_set_kernel_arg(kernel,10, int(log2(size), i4))

        if (allocated(this%projector_phases)) then
          call accel_set_kernel_arg(kernel, 11, this%buff_projector_phases)
          ! Note: we need to use this%nphase, as the kernel might be called with spiral=false, but
          !       the phases been built with spiralBC=true
          !> spiral BC ik should be the same for psi2b
          call accel_set_kernel_arg(kernel, 12, (psi1b%ik - st%d%kpt%start)*this%total_points*this%nphase)
        end if
        ! In case of CUDA we use an optimized kernel, in which the loop over npoints is broken
        ! further into chunks, in order to parallelize over the threads within a warp.
        ! Therefore we need to launch warp_size * size kernels. The size of each block needs to
        ! have multiples of warp_size as x-dimension.

        size_unfolded = size * accel%warp_size
        padnprojs = pad_pow2(this%max_nprojs)
        lnprojs = min(accel_kernel_workgroup_size(kernel)/accel%warp_size, padnprojs)

        call accel_kernel_run(kernel, &
          (/size_unfolded, padnprojs, int(this%nprojector_matrices, i8)/), (/int(accel%warp_size, i8), lnprojs, 1_i8/))

        do imat = 1, this%nprojector_matrices
          pmat => this%projector_matrices(imat)

          npoints = pmat%npoints
          nprojs = pmat%nprojs

          !! update number of operations for nphase !!
          call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst_linear*nprojs)
        end do

        call accel_finish()
      end do
      call accel_read_buffer(buff_projs,(ndim+1)*this%full_projection_size*psi1b%pack_size(1), projs)
      call accel_release_buffer(buff_projs)

      call profiling_out(cl_prof)
    end if

  end if

  if (mesh%parallel_in_domains) then
    call profiling_in(prof_matelement_reduce, TOSTRING(X(VNLPSI_MAT_ELEM_REDUCE)))
    call mesh%allreduce(projs)
    call profiling_out(prof_matelement_reduce)
  end if

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    iatom = this%projector_to_atom(imat)

    if (allocated(pmat%zmix)) then
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:psi1b%nst, 1:st%d%dim))

      do idir = 1, ndim
        do ist = 1, psi1b%nst
          tmp_proj(1:nprojs, ist, 1) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 1), &
            projs(idir, (ist-1)*st%d%dim+1, iprojection + 1:iprojection + nprojs)) &
            + matmul(pmat%zmix(1:nprojs, 1:nprojs, 3), &
            projs(idir, (ist-1)*st%d%dim+2, iprojection + 1:iprojection + nprojs))
          tmp_proj(1:nprojs, ist, 2) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 2), &
            projs(idir, (ist-1)*st%d%dim+2, iprojection + 1:iprojection + nprojs)) &
            + matmul(pmat%zmix(1:nprojs, 1:nprojs, 4), &
            projs(idir, (ist-1)*st%d%dim+1, iprojection + 1:iprojection + nprojs))
        end do

        do ist = 1, psi1b%nst
          do idim = 1, st%d%dim
            do iproj = 1, nprojs
              projs(idir, (ist-1)*st%d%dim+idim, iprojection + iproj) = tmp_proj(iproj, ist, idim)
            end do
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(tmp_proj)
#else
      ! Complex projector matrix not allowed for real batches
      ASSERT(.false.)
#endif
    else if (allocated(pmat%dmix)) then

      do idir = 1, ndim
        do ist = 1, nst_linear
          projs(idir, ist, iprojection + 1:iprojection + nprojs) = &
            matmul(pmat%dmix(1:nprojs, 1:nprojs), projs(idir, ist, iprojection + 1:iprojection + nprojs))
        end do
      end do
    end if

    SAFE_ALLOCATE(ff(1:ndim))

    ff(1:ndim) = CNST(0.0)

    do ii = 1, psi1b%nst_linear
      ist = psi1b%linear_to_ist(ii)
      if (st%d%kweights(iqn)*abs(st%occ(ist, iqn)) <= M_EPSILON) cycle
      do iproj = 1, nprojs
        do idir = 1, ndim
          ff(idir) = ff(idir) - M_TWO*st%d%kweights(iqn)*st%occ(ist, iqn)*pmat%scal(iproj)*mesh%volume_element* &
            R_REAL(R_CONJ(projs(0, ii, iprojection + iproj))*projs(idir, ii, iprojection + iproj))
        end do
      end do
    end do

    force(1:ndim, iatom) = force(1:ndim, iatom) + ff(1:ndim)

    call profiling_count_operations((R_ADD + 2*R_MUL)*nst_linear*ndim*nprojs)

    SAFE_DEALLOCATE_A(ff)

    iprojection = iprojection + nprojs

  end do

  SAFE_DEALLOCATE_A(projs)
  call profiling_out(prof_matelement)

  POP_SUB(X(hamiltonian_elec_base_nlocal_force))
end subroutine X(hamiltonian_elec_base_nlocal_force)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_position_commutator)(this, mesh, std, spiral_bnd, psib, commpsib)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  class(mesh_t),                         intent(in)    :: mesh
  type(states_elec_dim_t),               intent(in)    :: std
  logical,                               intent(in)    :: spiral_bnd
  type(wfs_elec_t),                      intent(in)    :: psib
  class(wfs_elec_t), target,             intent(inout) :: commpsib(:)

  integer :: ist, ip, iproj, imat, nreal, iprojection, idir
  integer :: npoints, nprojs, nst
  integer, allocatable :: ind(:)
  R_TYPE :: aa, bb, cc, dd
  R_TYPE, allocatable :: projections(:, :, :)
  R_TYPE, allocatable :: psi(:, :, :), lpsi(:,:)
  type(projector_matrix_t), pointer :: pmat
  type(profile_t), save :: prof, reduce_prof
  integer(i8) :: wgsize, size
  class(wfs_elec_t), pointer :: commpsib_(:)
#ifdef R_TCOMPLEX
  integer :: idim
  CMPLX :: phase
  CMPLX, allocatable :: tmp_proj(:, :, :)
#endif

  if (.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)

  PUSH_SUB(X(hamiltonian_elec_base_nlocal_position_commutator))
  call profiling_in(prof, TOSTRING(X(COMMUTATOR)))

  ASSERT(psib%is_packed())
  ASSERT(.not. spiral_bnd)

  nst = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif
  ! this pointer assigment is needed as a workaround for gcc 12
  ! otherwise it throws errors related to openmp and polymorphic arrays
  commpsib_ => commpsib

  if (psib%status() == BATCH_DEVICE_PACKED) then
    call X(commutator_accel)()
    call profiling_out(prof)
    POP_SUB(X(hamiltonian_elec_base_nlocal_position_commutator))
    return
  end if

  SAFE_ALLOCATE(projections(1:nst, 1:this%full_projection_size, 0:3))
  projections = M_ZERO

  SAFE_ALLOCATE(ind(1:this%nprojector_matrices))

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs
    ind(imat) = iprojection
    iprojection = iprojection + nprojs
    !    call profiling_count_operations(nprojs*(R_ADD + R_MUL)*npoints + nst*nprojs)
  end do

  !$omp parallel do private(imat, pmat, iprojection, npoints, nprojs, iproj, ist, aa, bb, cc, dd, ip, lpsi)
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    iprojection = ind(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if (npoints == 0) cycle

    SAFE_ALLOCATE(lpsi(1:npoints, 1:nst))
    if (.not. allocated(this%projector_phases)) then
      do ist = 1, nst
        !$omp simd
        do ip = 1, npoints
          lpsi(ip, ist) = psib%X(ff_pack)(ist, pmat%map(ip))
        end do
      end do
    else
#ifdef R_TCOMPLEX
      do ip = 1, npoints
        !$omp simd
        do ist = 1, nst
          lpsi(ip, ist) = psib%zff_pack(ist, pmat%map(ip)) &
            *this%projector_phases(ip, 1, imat, psib%ik)
        end do
      end do
#else
      ! Phases not allowed for real batches
      ASSERT(.false.)
#endif
    end if

    do iproj = 1, nprojs

      if (pmat%is_cmplx) then
#ifdef R_TCOMPLEX
        do ist = 1, nst
          aa = CNST(0.0)
          bb = CNST(0.0)
          cc = CNST(0.0)
          dd = CNST(0.0)

          !$omp simd reduction(+:aa, bb, cc, dd)
          do ip = 1, npoints
            aa = aa + R_CONJ(pmat%zprojectors(ip, iproj))*lpsi(ip, ist)
            bb = bb + R_CONJ(pmat%zprojectors(ip, iproj))*pmat%position(1, ip)*lpsi(ip, ist)
            cc = cc + R_CONJ(pmat%zprojectors(ip, iproj))*pmat%position(2, ip)*lpsi(ip, ist)
            dd = dd + R_CONJ(pmat%zprojectors(ip, iproj))*pmat%position(3, ip)*lpsi(ip, ist)
          end do
          projections(ist, iprojection + iproj, 0) = pmat%scal(iproj)*aa
          projections(ist, iprojection + iproj, 1) = pmat%scal(iproj)*bb
          projections(ist, iprojection + iproj, 2) = pmat%scal(iproj)*cc
          projections(ist, iprojection + iproj, 3) = pmat%scal(iproj)*dd
        end do
#else
        ! Complex projection matrix not allowed for real batches
        ASSERT(.false.)
#endif
      else
        do ist = 1, nst
          aa = CNST(0.0)
          bb = CNST(0.0)
          cc = CNST(0.0)
          dd = CNST(0.0)
          !$omp simd reduction(+:aa, bb, cc, dd)
          do ip = 1, npoints
            aa = aa + pmat%dprojectors(ip, iproj)*lpsi(ip, ist)
            bb = bb + pmat%dprojectors(ip, iproj)*pmat%position(1, ip)*lpsi(ip, ist)
            cc = cc + pmat%dprojectors(ip, iproj)*pmat%position(2, ip)*lpsi(ip, ist)
            dd = dd + pmat%dprojectors(ip, iproj)*pmat%position(3, ip)*lpsi(ip, ist)
          end do
          projections(ist, iprojection + iproj, 0) = pmat%scal(iproj)*aa
          projections(ist, iprojection + iproj, 1) = pmat%scal(iproj)*bb
          projections(ist, iprojection + iproj, 2) = pmat%scal(iproj)*cc
          projections(ist, iprojection + iproj, 3) = pmat%scal(iproj)*dd
        end do

      end if
    end do

    SAFE_DEALLOCATE_A(lpsi)
  end do
  !$omp end parallel do

  ! reduce the projections
  if (mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(COMMUTATOR_REDUCE)))
    call mesh%allreduce(projections)
    call profiling_out(reduce_prof)
  end if

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if (allocated(pmat%zmix)) then
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:psib%nst, 1:std%dim))

      do idir = 0, 3
        do ist = 1, psib%nst
          tmp_proj(1:nprojs, ist, 1) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 1), &
            projections((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs, idir)) &
            + matmul(pmat%zmix(1:nprojs, 1:nprojs, 3), &
            projections((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs, idir))
          tmp_proj(1:nprojs, ist, 2) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 2), &
            projections((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs, idir)) &
            + matmul(pmat%zmix(1:nprojs, 1:nprojs, 4), &
            projections((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs, idir))
        end do

        do ist = 1, psib%nst
          do idim = 1, std%dim
            do iproj = 1, nprojs
              projections((ist-1)*std%dim+idim, iprojection + iproj, idir) = tmp_proj(iproj, ist, idim)
            end do
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(tmp_proj)
#else
      ! Complex projection matrix not allowed for real batches
      ASSERT(.false.)
#endif
    else if (allocated(pmat%dmix)) then
      do idir = 0, 3
        do ist = 1, nst
          projections(ist, iprojection + 1:iprojection + nprojs, idir) = &
            matmul(pmat%dmix(1:nprojs, 1:nprojs), projections(ist, iprojection + 1:iprojection + nprojs, idir))
        end do
      end do
    end if

    if (npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints, 0:3))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))

      if (.not. pmat%is_cmplx) then
        do idir = 0, 3
          call blas_gemm('N', 'T', nreal, npoints, nprojs, &
            M_ONE, projections(1, iprojection + 1, idir), nreal, pmat%dprojectors(1, 1), npoints, &
            M_ZERO, psi(1, 1, idir), nreal)
        end do
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints*4)

      else
#ifdef R_TCOMPLEX
        do idir = 0, 3
          call blas_gemm('N', 'T', nst, npoints, nprojs, &
            M_z1, projections(1, iprojection + 1, idir), nst, pmat%zprojectors(1, 1), npoints, &
            M_z0, psi(1, 1, idir), nst)
        end do
#endif
        call profiling_count_operations(nst*nprojs*(R_ADD+R_MUL)*npoints*4)
      end if

      if (allocated(this%projector_phases)) then
#ifdef R_TCOMPLEX
        !$omp parallel private(ip, ist, phase)
        do idir = 0, 3
          !$omp do
          do ip = 1, npoints
            phase = conjg(this%projector_phases(ip, 1, imat, psib%ik))
            !$omp simd
            do ist = 1, nst
              psi(ist, ip, idir) = phase*psi(ist, ip, idir)
            end do
          end do
          !$omp end do nowait
        end do
        !$omp end parallel
        call profiling_count_operations(nst*npoints*3*R_MUL)
#else
        ! Phases not allowed for real batches
        ASSERT(.false.)
#endif
      end if

      !$omp parallel private(ip, ist)  if (.not. this%projector_self_overlap)
      do idir = 1, 3
        !$omp do
        do ip = 1, npoints
          do ist = 1, nst
            commpsib_(idir)%X(ff_pack)(ist, pmat%map(ip)) = commpsib_(idir)%X(ff_pack)(ist, pmat%map(ip)) &
              - psi(ist, ip, idir) + pmat%position(idir, ip)*psi(ist, ip, 0)
          end do
        end do
        !$omp end do nowait
      end do
      !$omp end parallel

      call profiling_count_operations(nst*npoints*3*(2*R_ADD+R_MUL))
    end if

    SAFE_DEALLOCATE_A(psi)

    iprojection = iprojection + nprojs
  end do

  SAFE_DEALLOCATE_A(ind)

  call profiling_out(prof)
  POP_SUB(X(hamiltonian_elec_base_nlocal_position_commutator))

contains

  subroutine X(commutator_accel)()
    type(accel_kernel_t), target, save :: dker_commutator_bra, dker_commutator_bra_phase, dker_mix
    type(accel_kernel_t), target, save :: dker_commutator_ket, dker_commutator_ket_phase
    type(accel_kernel_t), target, save :: zker_commutator_bra, zker_commutator_bra_phase, zker_mix
    type(accel_kernel_t), target, save :: zker_commutator_ket, zker_commutator_ket_phase
    type(accel_kernel_t), pointer :: kernel, ker_mix
    type(accel_mem_t), target :: buff_proj
    type(accel_mem_t), pointer :: buff_proj_copy
    integer(i8) :: padnprojs, lnprojs, size2
    integer :: iregion, nregions_self_overlap, iregion_self_overlap
    R_TYPE, allocatable :: proj(:)

    padnprojs = pad_pow2(this%max_nprojs)


    if (allocated(this%projector_phases)) then
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_bra_phase, 'projector_commutator.cl', 'zprojector_commutator_bra_phase',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_commutator_bra_phase
      else
        call accel_kernel_start_call(dker_commutator_bra_phase, 'projector_commutator.cl', 'dprojector_commutator_bra_phase',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_commutator_bra_phase
      end if
      size = psib%pack_size(1)
    else
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_bra, 'projector_commutator.cl', 'zprojector_commutator_bra',&
          flags = '-DRTYPE_COMPLEX')
        size = psib%pack_size(1)
        kernel => zker_commutator_bra
      else
        call accel_kernel_start_call(dker_commutator_bra, 'projector_commutator.cl', 'dprojector_commutator_bra',&
          flags = '-DRTYPE_DOUBLE')
        size = psib%pack_size_real(1)
        kernel => dker_commutator_bra
      end if
    end if

    call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, 4*this%full_projection_size*psib%pack_size(1))

    call accel_set_kernel_arg(kernel,  0, this%nprojector_matrices)
    call accel_set_kernel_arg(kernel,  1, this%buff_offsets)
    call accel_set_kernel_arg(kernel,  2, this%buff_matrices)
    call accel_set_kernel_arg(kernel,  3, this%buff_maps)
    call accel_set_kernel_arg(kernel,  4, this%buff_scals)
    call accel_set_kernel_arg(kernel,  5, this%buff_position)
    call accel_set_kernel_arg(kernel,  6, psib%ff_device)
    call accel_set_kernel_arg(kernel,  7, log2(int(size, i4)))
    call accel_set_kernel_arg(kernel,  8, buff_proj)
    call accel_set_kernel_arg(kernel,  9, log2(int(size, i4)))

    if (allocated(this%projector_phases)) then
      call accel_set_kernel_arg(kernel, 10, this%buff_projector_phases)
      call accel_set_kernel_arg(kernel, 11, (psib%ik - std%kpt%start)*this%total_points)
    end if

    lnprojs = min(accel_kernel_workgroup_size(kernel)/size, padnprojs)

    call accel_kernel_run(kernel, (/size, padnprojs, int(this%nprojector_matrices, i8)/), (/size, lnprojs, 1_i8/))

    call accel_finish()

    if (mesh%parallel_in_domains) then
      SAFE_ALLOCATE(proj(1:4*this%full_projection_size*psib%pack_size(1)))
      call accel_read_buffer(buff_proj, 4*this%full_projection_size*psib%pack_size(1), proj)
      call mesh%allreduce(proj)
      call accel_write_buffer(buff_proj, 4*this%full_projection_size*psib%pack_size(1), proj)
      SAFE_DEALLOCATE_A(proj)
    end if

    if (this%projector_mix) then

      SAFE_ALLOCATE(buff_proj_copy)

      if (allocated(this%projector_matrices(1)%zmix)) then
        call accel_kernel_start_call(zker_mix, 'projector_commutator.cl', 'zprojector_mix_commutator', flags = '-DRTYPE_COMPLEX')
        ker_mix => zker_mix
        size = psib%pack_size(1)
        size2 = size/2
        call accel_create_buffer(buff_proj_copy, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, &
          this%full_projection_size*size)
      else
        call accel_kernel_start_call(dker_mix, 'projector_commutator.cl', 'dprojector_mix_commutator', flags = '-DRTYPE_DOUBLE')
        ker_mix => dker_mix
        size = psib%pack_size_real(1)
        size2 = size
        call accel_create_buffer(buff_proj_copy, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, &
          this%full_projection_size*size)
      end if
      call accel_set_kernel_arg(ker_mix, 0, this%nprojector_matrices)
      call accel_set_kernel_arg(ker_mix, 1, this%buff_offsets)
      call accel_set_kernel_arg(ker_mix, 2, this%buff_mix)
      call accel_set_kernel_arg(ker_mix, 3, buff_proj)
      call accel_set_kernel_arg(ker_mix, 4, log2(int(size, i4)))
      call accel_set_kernel_arg(ker_mix, 5, buff_proj_copy)

      padnprojs = pad_pow2(this%max_nprojs)
      lnprojs = min(accel_kernel_workgroup_size(ker_mix)/size2, padnprojs)

      call accel_kernel_run(ker_mix, (/size2, padnprojs, int(this%nprojector_matrices, i8)/), (/size2, lnprojs, 1_i8/))

      call accel_finish()

    else

      buff_proj_copy => buff_proj

    end if

    if (allocated(this%projector_phases)) then
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)

      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_ket_phase, 'projector_commutator.cl', 'zprojector_commutator_ket_phase',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_commutator_ket_phase
      else
        call accel_kernel_start_call(dker_commutator_ket_phase, 'projector_commutator.cl', 'dprojector_commutator_ket_phase',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_commutator_ket_phase
      end if
      size = psib%pack_size(1)
    else
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_ket, 'projector_commutator.cl', 'zprojector_commutator_ket',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_commutator_ket
        size = psib%pack_size(1)
      else
        call accel_kernel_start_call(dker_commutator_ket, 'projector_commutator.cl', 'dprojector_commutator_ket',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_commutator_ket
        size = psib%pack_size_real(1)
      end if
    end if

    do iregion = 1, this%nregions

      do imat = this%regions(iregion), this%regions(iregion+1)-1

        nregions_self_overlap = this%projector_matrices(imat)%nregions

        do iregion_self_overlap = 1, nregions_self_overlap

          call accel_set_kernel_arg(kernel,  0, this%nprojector_matrices)
          call accel_set_kernel_arg(kernel,  1, imat - 1)
          call accel_set_kernel_arg(kernel,  2, this%buff_offsets)
          call accel_set_kernel_arg(kernel,  3, this%buff_matrices)
          call accel_set_kernel_arg(kernel,  4, this%buff_maps)
          call accel_set_kernel_arg(kernel,  5, this%buff_position)
          call accel_set_kernel_arg(kernel,  6, buff_proj_copy)
          call accel_set_kernel_arg(kernel,  7, log2(int(size, i4)))
          call accel_set_kernel_arg(kernel,  8, commpsib(1)%ff_device)
          call accel_set_kernel_arg(kernel,  9, commpsib(2)%ff_device)
          call accel_set_kernel_arg(kernel, 10, commpsib(3)%ff_device)
          call accel_set_kernel_arg(kernel, 11, log2(int(size, i4)))

          call accel_set_kernel_arg(kernel, 12, this%projector_matrices(imat)%regions(iregion_self_overlap) - 1)
          call accel_set_kernel_arg(kernel, 13, this%projector_matrices(imat)%regions(iregion_self_overlap+1) - 1)

          if (allocated(this%projector_phases)) then
            call accel_set_kernel_arg(kernel, 14, this%buff_projector_phases)
            call accel_set_kernel_arg(kernel, 15, (psib%ik - std%kpt%start)*this%total_points)
          end if

          wgsize = accel_kernel_workgroup_size(kernel)/size

          call accel_kernel_run(kernel, &
            (/size, pad(this%max_npoints, wgsize), 1_i8 /), &
            (/size, wgsize, 1_i8/))

          call accel_finish()

        end do

      end do

    end do

    if (this%projector_mix) then
      call accel_release_buffer(buff_proj_copy)
      SAFE_ALLOCATE(buff_proj_copy)
    end if

    call accel_release_buffer(buff_proj)

  end subroutine X(commutator_accel)

end subroutine X(hamiltonian_elec_base_nlocal_position_commutator)

! ---------------------------------------------------------------------------------------
! Accumulates to commpsib the result of x Vnl | \psib >
subroutine X(hamiltonian_elec_base_r_vnlocal)(this, mesh, std, spiral_bnd, psib, commpsib)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  class(mesh_t),                         intent(in)    :: mesh
  type(states_elec_dim_t),               intent(in)    :: std
  logical,                               intent(in)    :: spiral_bnd
  type(wfs_elec_t),                      intent(in)    :: psib
  class(wfs_elec_t), target,             intent(inout) :: commpsib(1:3)

  integer :: ist, ip, iproj, imat, nreal, iprojection, idir
  integer :: npoints, nprojs, nst
  integer, allocatable :: ind(:)
  R_TYPE :: aa
  R_TYPE, allocatable :: projections(:, :)
  R_TYPE, allocatable :: psi(:, :), lpsi(:,:)
  type(projector_matrix_t), pointer :: pmat
  type(profile_t), save :: prof, reduce_prof
  integer(i8) :: wgsize, size
  class(wfs_elec_t), pointer :: commpsib_(:)
#ifdef R_TCOMPLEX
  integer :: idim
  CMPLX :: phase
  CMPLX, allocatable :: tmp_proj(:, :, :)
#endif

  if (.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)

  PUSH_SUB(X(hamiltonian_elec_base_r_vnlocal))
  call profiling_in(prof, TOSTRING(X(R_VNL)))

  ASSERT(psib%is_packed())
  ASSERT(.not. spiral_bnd)

  nst = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif
  ! this pointer assigment is needed as a workaround for gcc 12
  ! otherwise it throws errors related to openmp and polymorphic arrays
  commpsib_ => commpsib

  if (psib%status() == BATCH_DEVICE_PACKED) then
    call X(commutator_accel)()
    call profiling_out(prof)
    POP_SUB(X(hamiltonian_elec_base_r_nlocal))
    return
  end if

  SAFE_ALLOCATE(projections(1:nst, 1:this%full_projection_size))
  projections = M_ZERO

  SAFE_ALLOCATE(ind(1:this%nprojector_matrices))

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs
    ind(imat) = iprojection
    iprojection = iprojection + nprojs
  end do

  !$omp parallel do private(imat, pmat, iprojection, npoints, nprojs, iproj, ist, aa, ip, lpsi)
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    iprojection = ind(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if (npoints == 0) cycle

    SAFE_ALLOCATE(lpsi(1:npoints, 1:nst))
    if (.not. allocated(this%projector_phases)) then
      do ist = 1, nst
        !$omp simd
        do ip = 1, npoints
          lpsi(ip, ist) = psib%X(ff_pack)(ist, pmat%map(ip))
        end do
      end do
    else
#ifdef R_TCOMPLEX
      do ip = 1, npoints
        !$omp simd
        do ist = 1, nst
          lpsi(ip, ist) = psib%zff_pack(ist, pmat%map(ip)) &
            *this%projector_phases(ip, 1, imat, psib%ik)
        end do
      end do
#else
      ! Phases not allowed for real batches
      ASSERT(.false.)
#endif
    end if

    do iproj = 1, nprojs

      if (pmat%is_cmplx) then
#ifdef R_TCOMPLEX
        do ist = 1, nst
          aa = CNST(0.0)

          !$omp simd reduction(+:aa)
          do ip = 1, npoints
            aa = aa + R_CONJ(pmat%zprojectors(ip, iproj))*lpsi(ip, ist)
          end do
          projections(ist, iprojection + iproj) = pmat%scal(iproj)*aa
        end do
#else
        ! Complex projection matrix not allowed for real batches
        ASSERT(.false.)
#endif
      else
        do ist = 1, nst
          aa = CNST(0.0)
          !$omp simd reduction(+:aa)
          do ip = 1, npoints
            aa = aa + pmat%dprojectors(ip, iproj)*lpsi(ip, ist)
          end do
          projections(ist, iprojection + iproj) = pmat%scal(iproj)*aa
        end do

      end if
    end do

    SAFE_DEALLOCATE_A(lpsi)
  end do
  !$omp end parallel do

  ! reduce the projections
  if (mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(R_VNL_REDUCE)))
    call mesh%allreduce(projections)
    call profiling_out(reduce_prof)
  end if

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if (allocated(pmat%zmix)) then
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:psib%nst, 1:std%dim))

      do ist = 1, psib%nst
        tmp_proj(1:nprojs, ist, 1) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 1), &
          projections((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs)) &
          + matmul(pmat%zmix(1:nprojs, 1:nprojs, 3), &
          projections((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs))
        tmp_proj(1:nprojs, ist, 2) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 2), &
          projections((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs)) &
          + matmul(pmat%zmix(1:nprojs, 1:nprojs, 4), &
          projections((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs))
      end do

      do ist = 1, psib%nst
        do idim = 1, std%dim
          do iproj = 1, nprojs
            projections((ist-1)*std%dim+idim, iprojection + iproj) = tmp_proj(iproj, ist, idim)
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(tmp_proj)
#else
      ! Complex projection matrix not allowed for real batches
      ASSERT(.false.)
#endif
    else if (allocated(pmat%dmix)) then
      do ist = 1, nst
        projections(ist, iprojection + 1:iprojection + nprojs) = &
          matmul(pmat%dmix(1:nprojs, 1:nprojs), projections(ist, iprojection + 1:iprojection + nprojs))
      end do
    end if

    if (npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))

      if (.not. pmat%is_cmplx) then
        call blas_gemm('N', 'T', nreal, npoints, nprojs, &
          M_ONE, projections(1, iprojection + 1), nreal, pmat%dprojectors(1, 1), npoints, &
          M_ZERO, psi(1, 1), nreal)
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints*4)

      else
#ifdef R_TCOMPLEX
        call blas_gemm('N', 'T', nst, npoints, nprojs, &
          M_z1, projections(1, iprojection + 1), nst, pmat%zprojectors(1, 1), npoints, &
          M_z0, psi(1, 1), nst)
#endif
        call profiling_count_operations(nst*nprojs*(R_ADD+R_MUL)*npoints*4)
      end if

      if (allocated(this%projector_phases)) then
#ifdef R_TCOMPLEX
        !$omp parallel do private(ip, ist, phase)
        do ip = 1, npoints
          phase = conjg(this%projector_phases(ip, 1, imat, psib%ik))
          !$omp simd
          do ist = 1, nst
            psi(ist, ip) = phase*psi(ist, ip)
          end do
        end do
        !$omp end parallel do
        call profiling_count_operations(nst*npoints*3*R_MUL)
#else
        ! Phases not allowed for real batches
        ASSERT(.false.)
#endif
      end if

      !$omp parallel private(ip, ist)  if (.not. this%projector_self_overlap)
      do idir = 1, 3
        !$omp do
        do ip = 1, npoints
          do ist = 1, nst
            commpsib_(idir)%X(ff_pack)(ist, pmat%map(ip)) = &
              commpsib_(idir)%X(ff_pack)(ist, pmat%map(ip)) + pmat%position(idir, ip)*psi(ist, ip)
          end do
        end do
        !$omp end do nowait
      end do
      !$omp end parallel

      call profiling_count_operations(nst*npoints*3*(2*R_ADD+R_MUL))
    end if

    SAFE_DEALLOCATE_A(psi)

    iprojection = iprojection + nprojs
  end do

  SAFE_DEALLOCATE_A(ind)

  call profiling_out(prof)
  POP_SUB(X(hamiltonian_elec_base_r_vnlocal))

contains

  subroutine X(commutator_accel)()
    type(accel_kernel_t), target, save :: dker_commutator_bra, dker_commutator_bra_phase, dker_mix
    type(accel_kernel_t), target, save :: dker_commutator_ket, dker_commutator_ket_phase
    type(accel_kernel_t), target, save :: zker_commutator_bra, zker_commutator_bra_phase, zker_mix
    type(accel_kernel_t), target, save :: zker_commutator_ket, zker_commutator_ket_phase
    type(accel_kernel_t), pointer :: kernel, ker_mix
    type(accel_mem_t), target :: buff_proj
    type(accel_mem_t), pointer :: buff_proj_copy
    integer(i8) :: padnprojs, lnprojs, size2
    integer :: iregion, nregions_self_overlap, iregion_self_overlap
    R_TYPE, allocatable :: proj(:)

    padnprojs = pad_pow2(this%max_nprojs)


    if (allocated(this%projector_phases)) then
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_bra_phase, 'projector.cl', 'zprojector_r_vnl_bra_phase',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_commutator_bra_phase
      else
        call accel_kernel_start_call(dker_commutator_bra_phase, 'projector.cl', 'dprojector_r_vnl_bra_phase',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_commutator_bra_phase
      end if
      size = psib%pack_size(1)
    else
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_bra, 'projector.cl', 'zprojector_r_vnl_bra',&
          flags = '-DRTYPE_COMPLEX')
        size = psib%pack_size(1)
        kernel => zker_commutator_bra
      else
        call accel_kernel_start_call(dker_commutator_bra, 'projector.cl', 'dprojector_r_`vnl_bra',&
          flags = '-DRTYPE_DOUBLE')
        size = psib%pack_size_real(1)
        kernel => dker_commutator_bra
      end if
    end if

    call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, this%full_projection_size*psib%pack_size(1))

    call accel_set_kernel_arg(kernel,  0, this%nprojector_matrices)
    call accel_set_kernel_arg(kernel,  1, this%buff_offsets)
    call accel_set_kernel_arg(kernel,  2, this%buff_matrices)
    call accel_set_kernel_arg(kernel,  3, this%buff_maps)
    call accel_set_kernel_arg(kernel,  4, this%buff_scals)
    call accel_set_kernel_arg(kernel,  5, psib%ff_device)
    call accel_set_kernel_arg(kernel,  6, log2(int(size, i4)))
    call accel_set_kernel_arg(kernel,  7, buff_proj)
    call accel_set_kernel_arg(kernel,  8, log2(int(size, i4)))

    if (allocated(this%projector_phases)) then
      call accel_set_kernel_arg(kernel, 9, this%buff_projector_phases)
      call accel_set_kernel_arg(kernel, 10, (psib%ik - std%kpt%start)*this%total_points)
    end if

    lnprojs = min(accel_kernel_workgroup_size(kernel)/size, padnprojs)

    call accel_kernel_run(kernel, (/size, padnprojs, int(this%nprojector_matrices, i8)/), (/size, lnprojs, 1_i8/))

    call accel_finish()

    if (mesh%parallel_in_domains) then
      SAFE_ALLOCATE(proj(1:this%full_projection_size*psib%pack_size(1)))
      call accel_read_buffer(buff_proj, this%full_projection_size*psib%pack_size(1), proj)
      call mesh%allreduce(proj)
      call accel_write_buffer(buff_proj, this%full_projection_size*psib%pack_size(1), proj)
      SAFE_DEALLOCATE_A(proj)
    end if

    if (this%projector_mix) then

      SAFE_ALLOCATE(buff_proj_copy)

      if (allocated(this%projector_matrices(1)%zmix)) then
        call accel_kernel_start_call(zker_mix, 'projector.cl', 'zprojector_mix', flags = '-DRTYPE_COMPLEX')
        ker_mix => zker_mix
        size = psib%pack_size(1)
        ! In the case of spinors, we need to distribute states in pairs, else access random memory
        size2 = size/2
        call accel_create_buffer(buff_proj_copy, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, &
          this%full_projection_size*size)
      else
        call accel_kernel_start_call(dker_mix, 'projector.cl', 'dprojector_mix', flags = '-DRTYPE_DOUBLE')
        ker_mix => dker_mix
        size = psib%pack_size_real(1)
        size2 = size
        call accel_create_buffer(buff_proj_copy, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, &
          this%full_projection_size*size)
      end if
      call accel_set_kernel_arg(ker_mix, 0, this%nprojector_matrices)
      call accel_set_kernel_arg(ker_mix, 1, this%buff_offsets)
      call accel_set_kernel_arg(ker_mix, 2, this%buff_mix)
      call accel_set_kernel_arg(ker_mix, 3, buff_proj)
      call accel_set_kernel_arg(ker_mix, 4, log2(int(size, i4)))
      call accel_set_kernel_arg(ker_mix, 5, buff_proj_copy)

      padnprojs = pad_pow2(this%max_nprojs)
      lnprojs = min(accel_kernel_workgroup_size(ker_mix)/size2, padnprojs)

      call accel_kernel_run(ker_mix, (/size2, padnprojs, int(this%nprojector_matrices, i8)/), (/size2, lnprojs, 1_i8/))

      call accel_finish()

    else

      buff_proj_copy => buff_proj

    end if

    if (allocated(this%projector_phases)) then
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)

      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_ket_phase, 'projector.cl', 'zprojector_r_vnl_ket_phase',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_commutator_ket_phase
      else
        call accel_kernel_start_call(dker_commutator_ket_phase, 'projector.cl', 'dprojector_r_vnl_ket_phase',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_commutator_ket_phase
      end if
      size = psib%pack_size(1)
    else
      if (this%projector_matrices(1)%is_cmplx) then
        call accel_kernel_start_call(zker_commutator_ket, 'projector.cl', 'zprojector_r_vnl_ket',&
          flags = '-DRTYPE_COMPLEX')
        kernel => zker_commutator_ket
        size = psib%pack_size(1)
      else
        call accel_kernel_start_call(dker_commutator_ket, 'projector.cl', 'dprojector_r_vnl_ket',&
          flags = '-DRTYPE_DOUBLE')
        kernel => dker_commutator_ket
        size = psib%pack_size_real(1)
      end if
    end if

    do iregion = 1, this%nregions

      do imat = this%regions(iregion), this%regions(iregion+1)-1

        nregions_self_overlap = this%projector_matrices(imat)%nregions

        do iregion_self_overlap = 1, nregions_self_overlap

          call accel_set_kernel_arg(kernel,  0, this%nprojector_matrices)
          call accel_set_kernel_arg(kernel,  1, imat - 1)
          call accel_set_kernel_arg(kernel,  2, this%buff_offsets)
          call accel_set_kernel_arg(kernel,  3, this%buff_matrices)
          call accel_set_kernel_arg(kernel,  4, this%buff_maps)
          call accel_set_kernel_arg(kernel,  5, this%buff_position)
          call accel_set_kernel_arg(kernel,  6, buff_proj_copy)
          call accel_set_kernel_arg(kernel,  7, log2(int(size, i4)))
          call accel_set_kernel_arg(kernel,  8, commpsib(1)%ff_device)
          call accel_set_kernel_arg(kernel,  9, commpsib(2)%ff_device)
          call accel_set_kernel_arg(kernel, 10, commpsib(3)%ff_device)
          call accel_set_kernel_arg(kernel, 11, log2(int(size, i4)))

          call accel_set_kernel_arg(kernel, 12, this%projector_matrices(imat)%regions(iregion_self_overlap) - 1)
          call accel_set_kernel_arg(kernel, 13, this%projector_matrices(imat)%regions(iregion_self_overlap+1) - 1)

          if (allocated(this%projector_phases)) then
            call accel_set_kernel_arg(kernel, 14, this%buff_projector_phases)
            call accel_set_kernel_arg(kernel, 15, (psib%ik - std%kpt%start)*this%total_points)
          end if

          wgsize = accel_kernel_workgroup_size(kernel)/size

          call accel_kernel_run(kernel, &
            (/size, pad(this%max_npoints, wgsize), 1_i8 /), &
            (/size, wgsize, 1_i8/))

          call accel_finish()

        end do

      end do

    end do

    if (this%projector_mix) then
      call accel_release_buffer(buff_proj_copy)
      SAFE_ALLOCATE(buff_proj_copy)
    end if

    call accel_release_buffer(buff_proj)

  end subroutine X(commutator_accel)

end subroutine X(hamiltonian_elec_base_r_vnlocal)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
