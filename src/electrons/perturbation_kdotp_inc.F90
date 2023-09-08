!! Copyright (C) 2007 X. Andrade
!! Copyright (C) 2021 N. Tancogne-Dejean
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

! --------------------------------------------------------------------------
!> Returns f_out = H' f_in, where H' is perturbation Hamiltonian
!! Note that e^ikr phase is applied to f_in, then is removed afterward
subroutine X(perturbation_kdotp_apply)(this, namespace, space, gr, hm, ik, f_in, f_out, set_bc)
  class(perturbation_kdotp_t), intent(in)    :: this
  type(namespace_t),           intent(in)    :: namespace
  type(space_t),               intent(in)    :: space
  type(grid_t),                intent(in)    :: gr
  type(hamiltonian_elec_t),    intent(in)    :: hm
  integer,                     intent(in)    :: ik
  R_TYPE,                      intent(in)    :: f_in(:, :)
  R_TYPE,                      intent(out)   :: f_out(:, :)
  logical,           optional, intent(in)    :: set_bc

  R_TYPE, allocatable :: f_in_copy(:, :)
  logical :: apply_kpoint, set_bc_
  integer :: ip, idim
  type(profile_t), save :: prof
  R_TYPE, allocatable :: grad(:, :, :), Hxpsi(:,:)
  integer :: iatom

  PUSH_SUB(X(perturbation_kdotp_apply))

  call profiling_in(prof, TOSTRING(X(PERT_KDOTP_APPLY)))

  ASSERT(this%dir /= -1)

  set_bc_ = optional_default(set_bc, .true.)

  SAFE_ALLOCATE(f_in_copy(1:gr%np_part, 1:hm%d%dim))
  if (set_bc_) then
    call lalg_copy(gr%np, hm%d%dim, f_in, f_in_copy)
    do idim = 1, hm%d%dim
      call boundaries_set(gr%der%boundaries, gr, f_in_copy(:,idim))
    end do
  else
    call lalg_copy(gr%np_part, hm%d%dim, f_in, f_in_copy)
  end if

  apply_kpoint = allocated(hm%hm_base%phase)
  if (this%vel_method == OPTION__KDOTPVELMETHOD__HCOM_VEL) then
    apply_kpoint = .false.
  end if

  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_in_copy, hm%hm_base%phase(1:gr%np_part, ik), gr%np_part, .false.)
#endif
  end if

  ! perturbation is grad + [V,r]
  if (this%vel_method /= OPTION__KDOTPVELMETHOD__HCOM_VEL) then
    SAFE_ALLOCATE(grad(1:gr%np, 1:space%dim, 1:hm%d%dim))

    do idim = 1, hm%d%dim
      call X(derivatives_grad) (gr%der, f_in_copy(:, idim), grad(:, :, idim), set_bc = .false.)
      ! set_bc done already separately
    end do

    do idim = 1, hm%d%dim
      call lalg_copy(gr%np, grad(:, this%dir, idim), f_out(:, idim))
    end do

    SAFE_DEALLOCATE_A(grad)

    ! i delta_H_k = i (-i*grad + k) . delta_k
    ! representation on psi is just grad . delta_k
    ! note that second-order term is left out
    if (this%use_nonlocalpps) then
      do iatom = 1, this%ions%natoms
        if (species_is_ps(this%ions%atom(iatom)%species)) then
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, gr%der%boundaries, hm%d%dim, this%dir, ik, f_in_copy, f_out)
        end if
      end do
    end if

  else

    SAFE_ALLOCATE(Hxpsi(1:gr%np,1:hm%d%dim))
    call X(hamiltonian_elec_apply_single)(hm, namespace, gr, f_in_copy(:,:), Hxpsi(:,:), 1, ik, set_bc = .false.)

    do idim = 1, hm%d%dim
      !$omp parallel do
      do ip = 1, gr%np
        f_out(ip,idim) = gr%x(ip,this%dir)*Hxpsi(ip,idim)
      end do

      !$omp parallel do
      do ip = 1, gr%np_part
        f_in_copy(ip,idim) = gr%x(ip,this%dir)*f_in_copy(ip,idim)
      end do
    end do

    call X(hamiltonian_elec_apply_single)(hm, namespace, gr, f_in_copy(:,:), Hxpsi(:,:), 1, ik, set_bc = .false.)

    call lalg_axpy(gr%np, hm%d%dim, -M_ONE, Hxpsi, f_out)
    SAFE_DEALLOCATE_A(Hxpsi)
  end if

  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_out, hm%hm_base%phase(1:gr%np, ik), gr%np, .true.)
#endif
  end if

  SAFE_DEALLOCATE_A(f_in_copy)

  call profiling_out(prof)
  POP_SUB(X(perturbation_kdotp_apply))
end subroutine X(perturbation_kdotp_apply)

! --------------------------------------------------------------------------
!> d^2/dki dkj (-(1/2) ki kj [ri,[rj,H]]) =
!! for i  = j : 1 - [ri,[rj,Vnl]]
!! for i != j : -(1/2) [ri,[rj,Vnl]]
!! Ref: Eq. 3 from M Cardona and FH Pollak, Phys. Rev. 142, 530-543 (1966)
!! Except everything is times minus one, since our kdotp perturbation is d/d(ik)
subroutine X(perturbation_kdotp_apply_order_2) (this, namespace, space, gr, hm, ik, f_in, f_out)
  class(perturbation_kdotp_t), intent(in)    :: this
  type(namespace_t),           intent(in)    :: namespace
  type(space_t),               intent(in)    :: space
  type(grid_t),                intent(in)    :: gr
  type(hamiltonian_elec_t),    intent(in)    :: hm
  integer,                     intent(in)    :: ik
  R_TYPE,                      intent(in)    :: f_in(:, :)
  R_TYPE,                      intent(out)   :: f_out(:, :)

  integer :: ip, idim
  R_TYPE, allocatable :: f_in_copy(:,:)
  integer :: iatom
  R_TYPE, allocatable :: cpsi(:,:)
  type(perturbation_kdotp_t), pointer :: perturbation_kdotp

  PUSH_SUB(X(perturbation_kdotp_apply_order_2))

  ASSERT(this%dir2 /= -1)

  SAFE_ALLOCATE(f_in_copy(1:gr%np_part, 1:hm%d%dim))
  call lalg_copy(gr%np, hm%d%dim, f_in, f_in_copy)

  ! kdotp has the perturbation written in terms of the periodic part with the phase

  if (this%vel_method /= OPTION__KDOTPVELMETHOD__HCOM_VEL) then
    f_out(1:gr%np, 1:hm%d%dim) = M_ZERO
    SAFE_ALLOCATE(cpsi(1:gr%np, 1:hm%d%dim))
    cpsi(1:gr%np, 1:hm%d%dim) = M_ZERO

    if (this%use_nonlocalpps) then

      do iatom = 1, this%ions%natoms
        if (species_is_ps(this%ions%atom(iatom)%species)) then
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, gr%der%boundaries, hm%d%dim, this%dir, ik, f_in_copy, cpsi(:, :))
        end if
      end do

      do idim = 1, hm%d%dim
        !$omp parallel do
        do ip = 1, gr%np
          f_out(ip, idim) = f_out(ip, idim) + gr%x(ip, this%dir2) * cpsi(ip, idim)
          f_in_copy(ip,idim) = gr%x(ip,this%dir2)*f_in_copy(ip,idim)
        end do
      end do

      cpsi(1:gr%np, 1:hm%d%dim) = M_ZERO
      do iatom = 1, this%ions%natoms
        if (species_is_ps(this%ions%atom(iatom)%species)) then
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, gr%der%boundaries, hm%d%dim, this%dir, ik, f_in_copy, cpsi(:, :))
        end if
      end do

      call lalg_axpy(gr%np, hm%d%dim, -M_ONE, cpsi, f_out)

    end if

    if (this%dir == this%dir2) then
      ! add delta_ij
      call lalg_axpy(gr%np, hm%d%dim, -M_ONE, f_in, f_out)
    end if

  else

    do idim = 1, hm%d%dim
      call boundaries_set(gr%der%boundaries, gr, f_in_copy(:,idim))
    end do

    SAFE_ALLOCATE(cpsi(1:gr%np,1:hm%d%dim))

    perturbation_kdotp => perturbation_kdotp_t(namespace, this%ions)
    call perturbation_kdotp%setup_dir(this%dir)
    call perturbation_kdotp%X(apply)(namespace, space, gr, hm, ik, f_in_copy, cpsi, set_bc=.false.)

    do idim = 1, hm%d%dim
      !$omp parallel do
      do ip = 1, gr%np
        f_out(ip,idim) = gr%x(ip,this%dir2)*cpsi(ip,idim)
      end do
      !$omp parallel do
      do ip = 1, gr%np_part
        f_in_copy(ip,idim) = gr%x(ip,this%dir2)*f_in_copy(ip,idim)
      end do
    end do

    call perturbation_kdotp%X(apply)(namespace, space, gr, hm, ik, f_in_copy, cpsi, set_bc=.false.)
    call lalg_axpy(gr%np, hm%d%dim, -M_ONE, cpsi, f_out)

    SAFE_DEALLOCATE_P(perturbation_kdotp)
  end if

  if (this%dir /= this%dir2) then
    call lalg_scal(gr%np, hm%d%dim, M_HALF, f_out)
  end if

  SAFE_DEALLOCATE_A(cpsi)

  SAFE_DEALLOCATE_A(f_in_copy)

  POP_SUB(X(perturbation_kdotp_apply_order_2))
end subroutine X(perturbation_kdotp_apply_order_2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
