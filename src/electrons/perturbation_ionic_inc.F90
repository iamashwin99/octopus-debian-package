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
subroutine X(perturbation_ionic_apply)(this, namespace, space, gr, hm, ik, f_in, f_out, set_bc)
  class(perturbation_ionic_t), intent(in)    :: this
  type(namespace_t),           intent(in)    :: namespace
  type(space_t),               intent(in)    :: space
  type(grid_t),                intent(in)    :: gr
  type(hamiltonian_elec_t),    intent(in)    :: hm
  integer,                     intent(in)    :: ik
  R_TYPE,                      intent(in)    :: f_in(:, :)
  R_TYPE,                      intent(out)   :: f_out(:, :)
  logical,        optional,    intent(in)    :: set_bc

  R_TYPE, allocatable :: f_in_copy(:, :)
  logical :: apply_kpoint, set_bc_
  integer :: idim, iatom, idir
  type(profile_t), save :: prof
  R_TYPE, allocatable  :: tmp(:)

  PUSH_SUB(X(perturbation_ionic_apply))

  call profiling_in(prof, TOSTRING(X(PERT_ION_APPLY)))

  ASSERT(this%dir /= -1)
  ASSERT(hm%d%dim == 1)

  set_bc_ = optional_default(set_bc, .true.)

  SAFE_ALLOCATE(f_in_copy(1:gr%np_part, 1:hm%d%dim))
  if (set_bc_) then
    call lalg_copy(gr%np, hm%d%dim, f_in, f_in_copy)
    do idim = 1, hm%d%dim
      call boundaries_set(gr%der%boundaries, gr, f_in_copy(:, idim))
    end do
  else
    call lalg_copy(gr%np_part, hm%d%dim, f_in, f_in_copy)
  end if

  apply_kpoint = allocated(hm%hm_base%phase)
  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_in_copy, hm%hm_base%phase(1:gr%np_part, ik), gr%np_part, .false.)
#endif
  end if

  f_out(1:gr%np, 1) = M_ZERO

  SAFE_ALLOCATE(tmp(1:gr%np))
  do iatom = 1, this%ions%natoms
    do idir = 1, this%ions%space%dim

      if (this%pure_dir .and. iatom /= this%atom1 .and. idir /= this%dir) cycle

      call X(ionic_perturbation)(gr, namespace, this%ions, hm, ik, f_in_copy(:, 1), tmp, iatom, idir)

      call lalg_axpy(gr%np, this%mix1(iatom, idir), tmp, f_out(:, 1))

    end do
  end do
  SAFE_DEALLOCATE_A(tmp)

  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_out, hm%hm_base%phase(1:gr%np, ik), gr%np, .true.)
#endif
  end if

  SAFE_DEALLOCATE_A(f_in_copy)

  call profiling_out(prof)
  POP_SUB(X(perturbation_ionic_apply))
end subroutine X(perturbation_ionic_apply)

 ! --------------------------------------------------------------------------
subroutine X(ionic_perturbation)(gr, namespace, ions, hm, ik, f_in, f_out, iatom, idir)
  type(grid_t),              intent(in)    :: gr
  type(namespace_t),         intent(in)    :: namespace
  type(ions_t),              intent(in)    :: ions
  type(hamiltonian_elec_t),  intent(in)    :: hm
  integer,                   intent(in)    :: ik
  R_TYPE,                    intent(in)    :: f_in(:)
  R_TYPE,                    intent(out)   :: f_out(:)
  integer,                   intent(in)    :: iatom, idir

  ! FIX ME: may need to tell derivatives_perform not to apply boundary conditions
  ! more things about ghost points may need to be done

  R_TYPE, allocatable :: grad(:,:), fin(:, :), fout(:, :)
  FLOAT,  allocatable :: vloc(:)
  integer :: ip

  PUSH_SUB(X(ionic_perturbation))

  ! The above derivatives are only valid for orthogonal cells
  ASSERT(.not. ions%latt%nonorthogonal)

  ! Formula: grad(V_nl) psi = grad(V_nl psi) - V_nl (grad psi)

  SAFE_ALLOCATE(vloc(1:gr%np))
  vloc(1:gr%np) = M_ZERO
  call epot_local_potential(hm%ep, namespace, ions%space, ions%latt, gr, ions%atom(iatom)%species, &
    ions%pos(:, iatom), iatom, vloc)

  SAFE_ALLOCATE(fin(1:gr%np_part, 1))
  call lalg_copy(gr%np_part, f_in, fin(:, 1))

  !d^T v |f>
  SAFE_ALLOCATE(fout(1:gr%np_part, 1))
  !$omp parallel do
  do ip = 1, gr%np
    fout(ip, 1) = vloc(ip)*fin(ip, 1)
  end do
  call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, 1, fin, fout, ik)
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fout(:,1), f_out)

  !v d |f>
  SAFE_ALLOCATE(grad(1:gr%np, 1))
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fin(:,1), grad(:,1))
  !$omp parallel do
  do ip = 1, gr%np
    fout(ip, 1) = vloc(ip)*grad(ip, 1)
  end do
  call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, 1, grad, fout, ik)
  !$omp parallel do
  do ip = 1, gr%np
    f_out(ip) = -f_out(ip) + fout(ip, 1)
  end do

  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(fin)
  SAFE_DEALLOCATE_A(fout)
  SAFE_DEALLOCATE_A(vloc)
  POP_SUB(X(ionic_perturbation))

end subroutine X(ionic_perturbation)

! --------------------------------------------------------------------------
subroutine X(perturbation_ionic_apply_order_2) (this, namespace, space, gr, hm, ik, f_in, f_out)
  class(perturbation_ionic_t), intent(in)    :: this
  type(namespace_t),           intent(in)    :: namespace
  type(space_t),               intent(in)    :: space
  type(grid_t),                intent(in)    :: gr
  type(hamiltonian_elec_t),    intent(in)    :: hm
  integer,                     intent(in)    :: ik
  R_TYPE,                      intent(in)    :: f_in(:, :)
  R_TYPE,                      intent(out)   :: f_out(:, :)

  integer :: idim
  R_TYPE, allocatable :: f_in_copy(:,:)
  logical :: apply_kpoint
  integer :: iatom, idir, jdir
  R_TYPE, allocatable  :: tmp(:)

  PUSH_SUB(X(perturbation_ionic_apply_order_2))

  ASSERT(this%dir2 /= -1)

  SAFE_ALLOCATE(f_in_copy(1:gr%np_part, 1:hm%d%dim))
  call lalg_copy(gr%np, hm%d%dim, f_in, f_in_copy)
  do idim = 1, hm%d%dim
    call boundaries_set(gr%der%boundaries, gr, f_in_copy(:, idim))
  end do

  apply_kpoint = allocated(hm%hm_base%phase)
  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_in_copy, hm%hm_base%phase(1:gr%np_part, ik), gr%np_part, .false.)
#endif
  end if

  ASSERT(hm%d%dim == 1)

  SAFE_ALLOCATE(tmp(1:gr%np))

  f_out(1:gr%np, 1) = M_ZERO

  do iatom = 1, this%ions%natoms
    do idir = 1, this%ions%space%dim
      do jdir = 1, this%ions%space%dim

        if (this%pure_dir &
          .and. iatom /= this%atom1 .and. idir /= this%dir &
          .and. iatom /= this%atom2 .and. jdir /= this%dir2) cycle

        call X(ionic_perturbation_order_2)(gr, namespace, this%ions, hm, ik, f_in_copy(:, 1), tmp, iatom, idir, jdir)

        call lalg_axpy(gr%np, this%mix1(iatom, idir)*this%mix2(iatom, jdir), tmp, f_out(:, 1))

      end do
    end do
  end do

  SAFE_DEALLOCATE_A(tmp)

  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_out, hm%hm_base%phase(1:gr%np, ik), gr%np, .true.)
#endif
  end if

  SAFE_DEALLOCATE_A(f_in_copy)

  POP_SUB(X(perturbation_ionic_apply_order_2))
end subroutine X(perturbation_ionic_apply_order_2)


! --------------------------------------------------------------------------
subroutine X(ionic_perturbation_order_2) (gr, namespace, ions, hm, ik, f_in, f_out, iatom, idir, jdir)
  type(grid_t),        intent(in)    :: gr
  type(namespace_t),   intent(in)    :: namespace
  type(ions_t),        intent(in)    :: ions
  type(hamiltonian_elec_t), intent(in) :: hm
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: f_in(:)
  R_TYPE,              intent(out)   :: f_out(:)
  integer,             intent(in)    :: iatom, idir, jdir

  ! FIXME: may need to tell derivatives_oper not to apply boundary conditions

  R_TYPE, allocatable :: fin(:, :)
  R_TYPE, allocatable :: tmp1(:, :), tmp2(:,:)
  FLOAT,  allocatable :: vloc(:)
  integer :: ip

  PUSH_SUB(X(ionic_perturbation_order_2))

  ! The above derivatives are only valid for orthogonal cells
  ASSERT(.not. ions%latt%nonorthogonal)

  SAFE_ALLOCATE( fin(1:gr%np_part, 1))
  SAFE_ALLOCATE(tmp1(1:gr%np_part, 1))
  SAFE_ALLOCATE(tmp2(1:gr%np_part, 1))
  SAFE_ALLOCATE(vloc(1:gr%np))

  !$omp parallel do
  do ip = 1, gr%np
    vloc(ip) = M_ZERO
  end do
  call epot_local_potential(hm%ep, namespace, ions%space, ions%latt, gr, ions%atom(iatom)%species, &
    ions%pos(:, iatom), iatom, vloc)

  call lalg_copy(gr%np_part, f_in, fin(:, 1))

  !di^T dj^T v |f>
  !$omp parallel do
  do ip = 1, gr%np
    tmp1(ip, 1) = vloc(ip)*fin(ip, 1)
  end do
  call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, 1, fin, tmp1, ik)
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, tmp1(:,1), tmp2(:,1))
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, tmp2(:,1), f_out)

  !di^T v dj |f>
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, fin(:,1), tmp1(:,1))
  !$omp parallel do
  do ip = 1, gr%np
    tmp2(ip, 1) = vloc(ip)*tmp1(ip, 1)
  end do
  call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, tmp2(:,1), tmp1(:,1))
  call lalg_axpy(gr%np, -M_ONE, tmp1(:,1), f_out)

  !dj^T v di |f>
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fin(:,1), tmp1(:,1))
  !$omp parallel do
  do ip = 1, gr%np
    tmp2(ip, 1) = vloc(ip)*tmp1(ip, 1)
  end do
  call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, tmp2(:,1), tmp1(:,1))
  call lalg_axpy(gr%np, -M_ONE, tmp1(:,1), f_out)

  !v di dj |f>
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fin(:,1), tmp1(:,1))
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, tmp1(:,1), tmp2(:,1))
  !$omp parallel do
  do ip = 1, gr%np
    tmp1(ip, 1) = vloc(ip)*tmp2(ip, 1)
  end do
  call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, 1, tmp2, tmp1, ik)
  call lalg_axpy(gr%np, M_ONE, tmp1(:,1), f_out)

  POP_SUB(X(ionic_perturbation_order_2))

end subroutine X(ionic_perturbation_order_2)

! --------------------------------------------------------------------------
subroutine X(ionic_pert_matrix_elements_2)(gr, namespace, space, ions, hm, ik, st, vib, factor, matrix)
  type(grid_t),        intent(in)    :: gr
  type(namespace_t),   intent(in)    :: namespace
  type(space_t),       intent(in)    :: space
  type(ions_t),        intent(in)    :: ions
  type(hamiltonian_elec_t), intent(in) :: hm
  integer,             intent(in)    :: ik
  type(states_elec_t), intent(in)    :: st
  type(vibrations_t),  intent(in)    :: vib
  FLOAT,               intent(in)    :: factor
  FLOAT,               intent(inout) :: matrix(:, :) !< this is an expectation value of a Hermitian operator

  integer :: ist, idim, ip
  integer :: imat, jmat, iatom, idir, jdir
  FLOAT, allocatable :: vloc(:)
  R_TYPE, allocatable :: gpsi(:, :, :), g2psi(:, :, :, :), tmp1(:, :), psi(:, :)
  FLOAT :: dot

  PUSH_SUB(X(ionic_pert_matrix_elements_2))

  ASSERT(.not. st%parallel_in_states)

  SAFE_ALLOCATE( vloc(1:gr%np))
  SAFE_ALLOCATE(psi(1:gr%np_part, 1:st%d%dim))
  SAFE_ALLOCATE( gpsi(1:gr%np_part, 1:st%d%dim, 1:space%dim))
  SAFE_ALLOCATE(g2psi(1:gr%np_part, 1:st%d%dim, 1:space%dim, 1:space%dim))
  SAFE_ALLOCATE( tmp1(1:gr%np, 1:st%d%dim))

  do ist = 1, st%nst

    call states_elec_get_state(st, gr, ist, ik, psi)

    do idim = 1, st%d%dim
      call X(derivatives_grad)(gr%der, psi(:, idim), gpsi(:, idim, :))
      do idir = 1, space%dim
        call X(derivatives_grad)(gr%der, gpsi(:, idim, idir), g2psi(:, idim, idir, :))
      end do
    end do

    ! This term applies only to matrix elements (iatom, idir; iatom, jdir)
    do imat = 1, vib%num_modes
      iatom = vibrations_get_atom(vib, imat)
      idir  = vibrations_get_dir (vib, imat)

      !$omp parallel do
      do ip = 1, gr%np
        vloc(ip) = M_ZERO
      end do
      call epot_local_potential(hm%ep, namespace, ions%space, ions%latt, gr, ions%atom(iatom)%species, &
        ions%pos(:, iatom), iatom, vloc)

      do jdir = 1, space%dim
        jmat = vibrations_get_index(vib, iatom, jdir)

        dot = M_ZERO

        !2<f|dj^T v di |f>
        do idim = 1, st%d%dim
          !$omp parallel do
          do ip = 1, gr%np
            tmp1(ip, idim) = vloc(ip)*gpsi(ip, idim, idir)
          end do
        end do
        call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, st%d%dim, gpsi(:, :, idir), tmp1, ik)
        dot = dot + M_TWO*R_REAL(X(mf_dotp)(gr, st%d%dim, gpsi(:, :, jdir), tmp1))

        !2<f|di^T dj^T v |f>
        do idim = 1, st%d%dim
          !$omp parallel do
          do ip = 1, gr%np
            tmp1(ip, idim) = vloc(ip)*psi(ip, idim)
          end do
        end do
        call X(project_psi)(gr, gr%der%boundaries, hm%ep%proj(iatom:iatom), 1, st%d%dim, psi, tmp1, ik)
        dot = dot + M_TWO*R_REAL(X(mf_dotp)(gr, st%d%dim, g2psi(:, :, idir, jdir), tmp1))

        matrix(jmat, imat) = matrix(jmat, imat) + dot*st%occ(ist, ik)*factor

      end do

    end do
  end do

  POP_SUB(X(ionic_pert_matrix_elements_2))
end subroutine X(ionic_pert_matrix_elements_2)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
