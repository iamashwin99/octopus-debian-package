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

! ---------------------------------------------------------
!> conjugate-gradients method.
subroutine X(eigensolver_cg) (namespace, mesh, st, hm, xc, pre, tol, niter, converged, ik, diff, &
  energy_change_threshold, orthogonalize_to_all, conjugate_direction, additional_terms, shift)
  type(namespace_t),        intent(in)    :: namespace
  class(mesh_t),            intent(in)    :: mesh
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(preconditioner_t),   intent(in)    :: pre
  type(xc_t),               intent(in)    :: xc
  FLOAT,                    intent(in)    :: tol
  integer,                  intent(inout) :: niter
  integer,                  intent(inout) :: converged
  integer,                  intent(in)    :: ik
  FLOAT,                    intent(out)   :: diff(:) !< (1:st%nst)
  FLOAT,                    intent(in)    :: energy_change_threshold
  logical,                  intent(in)    :: orthogonalize_to_all
  integer,                  intent(in)    :: conjugate_direction
  logical,                  intent(in)    :: additional_terms
  FLOAT, pointer, optional, intent(in)    :: shift(:,:)

  R_TYPE, allocatable :: h_psi(:,:), sd(:,:), sd_precond(:,:), cg(:,:), h_cg(:,:)
  R_TYPE, allocatable :: psi(:, :), psi2(:, :), sd_previous(:,:), psi_j(:,:)
  R_TYPE   :: sd_product, sd_product_previous, sd_product_mixed, gamma, overlap, dot
  FLOAT    :: cg_norm, sd_norm, alpha, beta, theta, old_energy, lam, lam_conj, cg_phi
  FLOAT    :: residue, residue_previous
  FLOAT    :: delta_e, first_delta_e
  FLOAT    :: stheta, ctheta
  FLOAT, allocatable :: chi(:, :), omega(:, :), fxc(:, :, :), lam_sym(:)
  FLOAT    :: integral_hartree, integral_xc
  integer  :: ist, jst, iter, maxter, idim, ip, isp, ixc, ib
  R_TYPE   :: sb(2)
  FLOAT    :: a0, b0, dsb(3)
  logical  :: fold_ ! use folded spectrum operator (H-shift)^2
  logical  :: add_xc_term
  logical  :: small_residual_change
  type(states_elec_group_t) :: hpsi_j

  PUSH_SUB(X(eigensolver_cg))

  ! if the optional shift argument is present, assume we are computing a folded spectrum
  fold_ =  present(shift)

  ! make sure the passed optional pointer is allocated
  if(fold_) then
    ASSERT(associated(shift))
  end if

  residue_previous = M_HUGE

  ! do we add the XC term? needs derivatives of the XC functional
  add_xc_term = additional_terms
  if(st%d%ispin == UNPOLARIZED) then
    isp = 1
  else
    isp = 2
  end if
  do ixc = 1, 2
    if(bitand(xc%kernel(ixc, isp)%flags, XC_FLAGS_HAVE_FXC) == 0) then
      add_xc_term = .false.
    end if
  end do
  if(bitand(xc%kernel_family, XC_FAMILY_LDA) == 0) then
    add_xc_term = .false.
  end if
  ! TODO: extend to spinors
  if(st%d%ispin == SPINORS) then
    add_xc_term = .false.
  end if

  maxter = niter
  niter = 0

  SAFE_ALLOCATE(psi(1:mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(cg(1:mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(sd(1:mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(sd_precond(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(sd_previous(1:mesh%np, 1:st%d%dim))
  if(additional_terms) then
    SAFE_ALLOCATE(chi(1:mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(omega(1:mesh%np_part, 1:st%d%dim))
    if(st%d%ispin == UNPOLARIZED) then
      SAFE_ALLOCATE(fxc(1:mesh%np, 1:1, 1:1))
    else if(st%d%ispin == SPIN_POLARIZED) then
      SAFE_ALLOCATE(fxc(1:mesh%np, 1:2, 1:2))
    end if
  end if
  if(fold_) then
    SAFE_ALLOCATE( psi2(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(h_psi(1:mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE( h_cg(1:mesh%np_part, 1:st%d%dim))
  else
    SAFE_ALLOCATE(h_psi(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE( h_cg(1:mesh%np, 1:st%d%dim))
  end if

  if(hm%theory_level == RDMFT) then
    SAFE_ALLOCATE(psi_j(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(lam_sym(1:st%nst))
    call states_elec_group_copy(st%d, st%group, hpsi_j, copy_data=.false.)
    do ib = hpsi_j%block_start, hpsi_j%block_end
      call X(hamiltonian_elec_apply_batch) (hm, namespace, mesh, st%group%psib(ib, ik), hpsi_j%psib(ib, ik))
    end do
  end if

  h_psi = R_TOTYPE(M_ZERO)
  cg = R_TOTYPE(M_ZERO)
  sd = R_TOTYPE(M_ZERO)
  sd_precond = R_TOTYPE(M_ZERO)
  h_cg = R_TOTYPE(M_ZERO)
  sd_previous = R_TOTYPE(M_ZERO)
  residue_previous = M_ZERO

  ! get derivative once here -> the density does not change in the loop
  if(add_xc_term) then
    fxc = M_ZERO
    call xc_get_fxc(xc, mesh, namespace, st%rho, st%d%ispin, fxc)
  end if

  ! Set the diff to zero, since it is intent(out)
  diff(1:st%nst) = M_ZERO

  ! Start of main loop, which runs over all the eigenvectors searched
  ASSERT(converged >= 0)

  ! The steps in this loop follow closely the algorithm from
  ! Payne et al. (1992), Rev. Mod. Phys. 64, 4, section V.B
  eigenfunction_loop : do ist = converged + 1, st%nst
    sd_product_mixed = R_TOTYPE(M_ZERO)

    call states_elec_get_state(st, mesh, ist, ik, psi)
    ! Orthogonalize starting eigenfunctions to those already calculated...
    if(ist > 1) call X(states_elec_orthogonalize_single)(st, mesh, ist - 1, ik, psi, normalize = .true.)
    !if (ist > 1) call X(states_elec_orthogonalize_single_batch)(st, mesh, ist - 1, ik, psi, normalize = .true.)

    ! Calculate starting gradient: |hpsi> = H|psi>
    call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, psi, h_psi, ist, ik)

    if(fold_) then
      call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, h_psi, psi2, ist, ik)
      ! h_psi = (H-shift)^2 psi
      do idim = 1, st%d%dim
        !$omp parallel do simd schedule(static)
        do ip = 1, mesh%np
          h_psi(ip, idim) = psi2(ip, idim) - M_TWO*shift(ist,ik)*h_psi(ip, idim) + shift(ist,ik)**2*psi(ip, idim)
        end do
      end do
    end if

    ! Calculates starting eigenvalue: e(p) = <psi(p)|H|psi>
    st%eigenval(ist, ik) = R_REAL(X(mf_dotp) (mesh, st%d%dim, psi, h_psi))
    old_energy = st%eigenval(ist, ik)

    small_residual_change = .false.
    ! Starts iteration for this band
    iter_loop: do iter = 1, maxter
      ! need to save sd from previous iteration for Polak-Ribiere method
      if(conjugate_direction == OPTION__CGDIRECTION__POLAK) then
        if(iter /= 1) then
          sd_previous = sd
        else
          sd_previous = M_ZERO
        end if
      end if

      ! PTA92, eq. 5.10 (get steepest descent vector)
      do idim = 1, st%d%dim
        !$omp parallel do simd schedule(static)
        do ip = 1, mesh%np
          sd(ip, idim) = -h_psi(ip, idim) + st%eigenval(ist, ik)*psi(ip, idim)
        end do
      end do

      if (hm%theory_level == RDMFT) then
        ! For RDMFT, the gradient of the total energy functional differs from
        ! the DFT and HF cases, as the lagrange multiplier matrix lambda cannot
        ! be diagonalized together with the Hamiltonian. This is because the
        ! orbitals of the minimization are not the eigenstates of the
        ! single-body Hamiltonian, but of the systems 1RDM.
        ! The functional to be minimized here is:
        !  F = E[psi_i]-sum_ij lam_ij (<psi_i|psi_j> - delta_ij) + const.
        ! The respective gradient reads:
        !   dF/dpsi_i= dE/dphi_i - sum_j lam_ij |psi_j>= H|psi_i> - sum_j lam_ij |psi_j>
        ! We get the expression for lam_ij from the gradient with respect to
        ! psi*: lam_ij = <psi_i|dE/dpsi_j^*> = <psi_i|H|psi_j>
        ! NB: lam_ij != lam_ji until SCF convergence!
        do jst = 1, st%nst
          if (jst == ist) then
            lam_sym(jst) = M_TWO*st%eigenval(jst, ik)
          else
            call states_elec_get_state(st, mesh, jst, ik, psi_j)
            do idim = 1, st%d%dim
              call batch_get_state(hpsi_j%psib(hpsi_j%iblock(jst, ik), ik), (/jst, idim/), mesh%np, h_cg(:, idim))
            end do

            ! calculate <phi_j|H|phi_i> = lam_ji and <phi_i|H|phi_j> = lam_ij
            lam = R_REAL(X(mf_dotp) (mesh, st%d%dim, psi_j, h_psi))
            lam_conj = R_REAL(X(mf_dotp) (mesh, st%d%dim, psi, h_cg))
            lam_sym(jst) = lam + lam_conj

            do idim = 1, st%d%dim
              call lalg_axpy(mesh%np, lam_conj, psi_j(:, idim), sd(:, idim))
            end do
          end if
        end do
      end if

      ! PTA92, eq. 5.12
      ! Orthogonalize to all states -> seems not to be needed
      !call X(states_elec_orthogonalize_single_batch)(st, mesh, ist - 1, ik, sd, normalize = .true., &
      !    against_all=orthogonalize_to_all)

      ! PTA92, eq. 5.17
      ! Approximate inverse preconditioner
      call  X(preconditioner_apply)(pre, namespace, mesh, hm, sd(:,:), sd_precond(:,:), ik)

      ! PTA92, eq. 5.18
      dot = X(mf_dotp) (mesh, st%d%dim, psi, sd_precond)
      ! This needs to be done before the orthogonalization_single call, as psi is not guaranted
      ! to be orthogonal to the other bands here
      call lalg_axpy(mesh%np, st%d%dim, -dot, psi, sd_precond)

      ! orthogonalize against previous or all states, depending on the optional argument orthogonalize_to_all
      call X(states_elec_orthogonalize_single_batch)(st, mesh, ist - 1, ik, sd_precond, normalize = .false., &
        against_all=orthogonalize_to_all)

      sd_norm = X(mf_nrm2)(mesh, st%d%dim, sd_precond)
      if(sd_norm < CNST(1e-150)) then
        if (debug%info) then
          write(message(1), '(a,i8,a,i4,a,i4,a,es13.6)') 'Debug: CG Eigensolver - zero norm - ik', ik, &
            ' ist ', ist, ' iter ', iter, ' res ', residue
          call messages_info(1, namespace=namespace)
        end if
        exit
      else
        call lalg_scal(mesh%np, st%d%dim, M_ONE/sd_norm, sd_precond)
      end if

      ! dot products needed for conjugate gradient
      sd_product = X(mf_dotp) (mesh, st%d%dim, sd_precond, sd, reduce = .false.)
      if(iter /= 1 .and. conjugate_direction == OPTION__CGDIRECTION__POLAK) then
        ! only needed for Polak-Ribiere
        sd_product_mixed = X(mf_dotp) (mesh, st%d%dim, sd_precond, sd_previous, reduce = .false.)
      end if

      if(mesh%parallel_in_domains) then
        sb(1) = sd_product_mixed
        sb(2) = sd_product
        call mesh%allreduce(sb, dim = 2)
        sd_product_mixed = sb(1)
        sd_product  = sb(2)
      end if

      if(iter  ==  1) then
        ! simply use the preconditioned steepest descent vector in the first iteration
        sd_product_previous = sd_product
        call lalg_copy(mesh%np, st%d%dim, sd_precond, cg)
      else
        ! compute conjugate gradient vector from current and previous iteration
        select case (conjugate_direction)
        case (OPTION__CGDIRECTION__FLETCHER)
          ! PTA eq. 5.20
          gamma = sd_product/sd_product_previous        ! (Fletcher-Reeves)
        case (OPTION__CGDIRECTION__POLAK)
          gamma = (sd_product - sd_product_mixed)/sd_product_previous   ! (Polack-Ribiere)
        case default
          call messages_input_error(namespace, 'Conjugate Direction')
        end select
        ! save for next iteration
        sd_product_previous = sd_product

        ! PTA92, eq. 5.19
        do idim = 1, st%d%dim
          !$omp parallel do simd schedule(static)
          do ip = 1, mesh%np
            cg(ip, idim) = sd_precond(ip, idim) + gamma*cg(ip, idim)
          end do
        end do
        call profiling_count_operations(st%d%dim*mesh%np*(2*R_ADD + 2*R_MUL))
      end if

      ! PTA92, eq. 5.21
      overlap =  X(mf_dotp) (mesh, st%d%dim, psi, cg)
      call lalg_axpy(mesh%np, st%d%dim, -overlap, psi, cg)

      ! normalize cg here (PTA92, eq. 5.22)
      cg_norm = X(mf_nrm2) (mesh, st%d%dim, cg)
      if(cg_norm > CNST(1e-150)) then
        call lalg_scal(mesh%np, st%d%dim, M_ONE/cg_norm, cg)
      else
        if (debug%info) then
          write(message(1), '(a,i8,a,i4,a,i4,a,es13.6)') 'Debug: CG Eigensolver - zero norm - ik', ik, &
            ' ist ', ist, ' iter ', iter, ' res ', residue
          call messages_info(1, namespace=namespace)
        end if
        exit
      end if

      ! cg contains now the conjugate gradient
      call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, cg, h_cg, ist, ik)

      if(fold_) then
        call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, h_cg, psi2, ist, ik)
        ! h_psi = (H-shift)^2 psi
        do idim = 1, st%d%dim
          h_cg(1:mesh%np, idim) = psi2(1:mesh%np, idim) - M_TWO*shift(ist,ik)*h_cg(1:mesh%np, idim) &
            + shift(ist,ik)**2*cg(1:mesh%np, idim)
        end do
      end if

      ! Line minimization (eq. 5.23 to 5.38)
      a0 = M_TWO*R_REAL(X(mf_dotp) (mesh, st%d%dim, psi, h_cg, reduce = .false.)) !Eq. 5.26
      b0  = R_REAL(X(mf_dotp) (mesh, st%d%dim, cg, h_cg, reduce = .false.))

      if(mesh%parallel_in_domains) then
        dsb(1) = a0
        dsb(2) = b0
        call mesh%allreduce(dsb, dim=2)
        a0 = dsb(1)
        b0 = dsb(2)
      end if

      ! compare eq. 5.31
      alpha = M_TWO * (st%eigenval(ist, ik) - b0)

      !No contribution from unoccupied states
      if (additional_terms .and. abs(st%d%kweights(ik)*st%occ(ist, ik)) > M_EPSILON) then
        ! more terms here, see PTA92 eqs 5.31, 5.32, 5.33, 5.36
        ! Hartree term
        do idim = 1, st%d%dim
          !$omp parallel do simd schedule(static)
          do ip = 1, mesh%np
            chi(ip, idim) = M_TWO * R_REAL(R_CONJ(cg(ip, idim)) * psi(ip, idim))
          end do
        end do
        call dpoisson_solve(hm%psolver, namespace, omega(:, 1), chi(:, 1), all_nodes = .false.)
        integral_hartree = dmf_dotp(mesh, st%d%dim, chi, omega)

        ! exchange term
        ! TODO: adapt to different spin cases
        if(add_xc_term) then
          integral_xc = dmf_dotp(mesh, st%d%dim, fxc(:, :, 1), chi(:, :)**2)
        else
          integral_xc = M_ZERO
        end if

        ! add additional terms to alpha (alpha is -d2e/dtheta2 from eq. 5.31)
        alpha = alpha - (st%d%kweights(ik)*st%occ(ist, ik))**2 * (integral_hartree + integral_xc)
      end if

      beta = a0 * M_TWO

      ! For RDMFT, we get a different formula for the line minimization, which turns out to
      ! only change the beta of the original expression.
      ! beta -> beta + beta_rdmft, with beta_rdmft= - sum_j (lam_ji <cg_i|phi_k> + c.c.)
      if(hm%theory_level == RDMFT) then
        do jst = 1, st%nst
          call states_elec_get_state(st, mesh, jst, ik, psi_j)
          cg_phi = M_TWO*R_REAL(X(mf_dotp) (mesh, st%d%dim, psi_j, cg))
          beta = beta - cg_phi * lam_sym(jst)
        end do
      end if

      ! Eq. 5.37
      theta = atan(beta/alpha)*M_HALF
      ! Choose the minimum solutions.
      ! theta is in the range [-pi/2:pi/2] and we want the solution
      ! which is in between 0 and pi/2
      ! The sign of theta is given by the sign of alpha
      ! However, sometimes theta is slightly positive, in this case we still want the
      ! same solution.
      if (alpha/abs(st%eigenval(ist, ik)) > CNST(0.1)) then
        ctheta = -sin(theta)
        stheta = cos(theta)
      else
        ctheta = cos(theta)
        stheta = sin(theta)
      end if

      ! PTA92, eq. 5.38
      do idim = 1, st%d%dim
        !$omp parallel do simd schedule(static)
        do ip = 1, mesh%np
          psi(ip, idim) = ctheta*psi(ip, idim) + stheta*cg(ip, idim)
          h_psi(ip, idim) = ctheta*h_psi(ip, idim) + stheta*h_cg(ip, idim)
        end do
      end do

      call profiling_count_operations(st%d%dim*mesh%np*(2*R_ADD + 4*R_MUL))

      ! New eigenvalue
      st%eigenval(ist, ik) = st%eigenval(ist, ik)*ctheta**2 + stheta**2*b0 + ctheta*stheta*a0

      residue = X(states_elec_residue)(mesh, st%d%dim, h_psi, st%eigenval(ist, ik), psi)
      ! We compute it analytically, as st%eigenval(ist, ik) - old_energy
      ! The direct evaluation can lead to numerically zero if the difference is the eigenvalues
      ! are very close (theta close to zero).
      ! This is avoided by computing analytically the difference
      delta_e = stheta*(-st%eigenval(ist, ik)*stheta + stheta*b0 + ctheta*a0)
      if(iter == 1) first_delta_e = delta_e

      if (hm%theory_level == RDMFT) then
        do idim = 1, st%d%dim
          call batch_set_state(hpsi_j%psib(hpsi_j%iblock(ist, ik), ik), (/ist, idim/), mesh%np, h_psi(:, idim))
        end do
      end if

      if(debug%info) then
        write(message(1), '(a,i4,a,i4,a,i4,a,i4,a,es12.5,a,es12.5,a,es12.5,a,es12.5,a,2es12.5)') &
          'Debug: CG Eigensolver - ik', ik, ' ist ', ist, &
          ' iter ', iter, ' max ', maxter, &
          ' deltae ', abs(delta_e)/old_energy, &
          ' theta ', theta, &
          ' alpha ', alpha, &
          ' beta ', beta, &
          ' residue ', residue
        call messages_info(1)
      end if

      ! consider change in energy
      if (iter == 1) then
        first_delta_e = abs(delta_e)
      end if

      if (iter > 1) then
        ! This criterion is discussed in Sec. V.B.6
        if (abs(delta_e) < first_delta_e*energy_change_threshold) then
          exit iter_loop
        end if
      end if
      old_energy = st%eigenval(ist, ik)

      ! Test convergence.
      if (residue < tol) then
        ! require residue below tolerance for two consecutive steps
        if (iter > 1 .and. residue_previous < tol) then
          if (converged == ist - 1) converged = ist ! only consider the first converged eigenvectors
          exit iter_loop
        end if
      end if
      residue_previous = residue
    end do iter_loop

    ! if the folded operator was used, compute the actual eigenvalue
    if(fold_) then
      call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, psi, h_psi, ist, ik)
      st%eigenval(ist, ik) = R_REAL(X(mf_dotp) (mesh, st%d%dim, psi, h_psi, reduce = .true.))
      residue = X(states_elec_residue)(mesh, st%d%dim, h_psi, st%eigenval(ist, ik), psi)
    end if

    if (.not. orthogonalize_to_all) then
      ! do the orthogonalization here only if not done during the iteration
      call X(states_elec_orthogonalize_single_batch)(st, mesh, ist - 1, ik, psi, &
        normalize = .true., against_all=.false.)
    end if
    call states_elec_set_state(st, mesh, ist, ik, psi)

    niter = niter + iter + 1

    diff(ist) = residue

    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      call loct_progress_bar(st%lnst*(ik - 1) +  ist, st%lnst*st%d%kpt%nlocal)
    end if

  end do eigenfunction_loop

  ! Deallocation of variables
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(h_psi)
  SAFE_DEALLOCATE_A(sd)
  SAFE_DEALLOCATE_A(sd_precond)
  SAFE_DEALLOCATE_A(sd_previous)
  SAFE_DEALLOCATE_A(cg)
  SAFE_DEALLOCATE_A(h_cg)
  SAFE_DEALLOCATE_A(chi)
  SAFE_DEALLOCATE_A(omega)
  SAFE_DEALLOCATE_A(fxc)
  SAFE_DEALLOCATE_A(psi2)
  if(hm%theory_level == RDMFT) then
    SAFE_DEALLOCATE_A(psi_j)
    SAFE_DEALLOCATE_A(lam_sym)
    call states_elec_group_end(hpsi_j, st%d)
  end if

  POP_SUB(X(eigensolver_cg))
end subroutine X(eigensolver_cg)


! ---------------------------------------------------------
!> The algorithm is essentially taken from Jiang et al. Phys. Rev. B 68, 165337 (2003).
subroutine X(eigensolver_cg_jiang) (namespace, mesh, st, hm, tol, niter, converged, ik, diff)
  type(namespace_t),        intent(in)    :: namespace
  class(mesh_t),            intent(in)    :: mesh
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  FLOAT,                    intent(in)    :: tol
  integer,                  intent(inout) :: niter
  integer,                  intent(inout) :: converged
  integer,                  intent(in)    :: ik
  FLOAT,          optional, intent(out)   :: diff(:) !< (1:st%nst)

  integer :: nst, dim, ist, maxter, i, conv, ip, idim
  R_TYPE, allocatable :: psi(:,:), phi(:, :), hcgp(:, :), cg(:, :), sd(:, :), cgp(:, :)
  FLOAT :: ctheta, stheta, ctheta2, stheta2, mu, lambda, dump, &
    gamma, sol(2), alpha, beta, theta, theta2, res, norm
  R_TYPE :: dot

  PUSH_SUB(X(eigensolver_cg_jiang))

  dim = st%d%dim
  nst = st%nst

  maxter = niter
  niter = 0

  SAFE_ALLOCATE( phi(1:mesh%np     , 1:dim))
  SAFE_ALLOCATE( psi(1:mesh%np_part, 1:dim))
  SAFE_ALLOCATE(  cg(1:mesh%np     , 1:dim))
  SAFE_ALLOCATE(hcgp(1:mesh%np     , 1:dim))
  SAFE_ALLOCATE(  sd(1:mesh%np     , 1:dim))
  SAFE_ALLOCATE( cgp(1:mesh%np_part, 1:dim))

  phi(1:mesh%np, 1:dim) = R_TOTYPE(M_ZERO)
  psi(1:mesh%np, 1:dim) = R_TOTYPE(M_ZERO)
  cgp(1:mesh%np, 1:dim) = R_TOTYPE(M_ZERO)

  ! Set the diff to zero, since it is intent(out)
  if(present(diff)) diff(1:st%nst) = M_ZERO

  conv = converged
  states: do ist = conv + 1, nst

    call states_elec_get_state(st, mesh, ist, ik, psi)

    ! Orthogonalize starting eigenfunctions to those already calculated...
    if (ist > 1) call X(states_elec_orthogonalize_single)(st, mesh, ist - 1, ik, psi, normalize = .true.)

    ! Calculate starting gradient: |hpsi> = H|psi>
    call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, psi, phi, ist, ik)
    niter = niter + 1

    ! Initial settings for scalar variables.
    ctheta = M_ONE
    stheta = M_ZERO
    mu     = M_ONE

    ! Initialize to zero the vector variables.
    hcgp = R_TOTYPE(M_ZERO)
    cg   = R_TOTYPE(M_ZERO)


    band: do i = 1, maxter - 1 ! One operation has already been made.

      if( i >1 ) then ! Get H|psi> (through the linear formula)
        do idim = 1, st%d%dim
          do ip = 1, mesh%np
            phi(ip, idim) = ctheta*phi(ip, idim) + stheta*hcgp(ip, idim)
          end do
        end do
      end if

      ! lambda = <psi|H|psi> = <psi|phi>
      lambda = R_REAL(X(mf_dotp)(mesh, dim, psi, phi))

      ! Check convergence
      res = X(states_elec_residue)(mesh, dim, phi, lambda, psi)

      if(debug%info) then
        norm = X(mf_nrm2)(mesh, dim, phi)
        write(message(1), '(a,i4,a,i4,a,i4,a,es13.6,a,es13.6)') 'Debug: CG New Eigensolver - ik', ik, &
          ' ist ', ist, ' iter ', i + 1, ' res ', res, ' ', res/norm
        call messages_info(1, namespace=namespace)
      end if

      if(present(diff)) diff(ist) = res
      if(res < tol) then
        if(conv == ist - 1) conv = ist
        exit band
      end if

      ! Get steepest descent vector
      do idim = 1, st%d%dim
        do ip = 1, mesh%np
          sd(ip, idim) = lambda*psi(ip, idim) - phi(ip, idim)
        end do
      end do

      if (ist > 1) call X(states_elec_orthogonalize_single)(st, mesh, ist - 1, ik, sd, normalize = .false.)

      ! Get conjugate-gradient vector
      dump = X(mf_nrm2)(mesh, dim, sd)**2
      gamma = dump/mu
      mu    = dump

      do idim = 1, st%d%dim
        do ip = 1, mesh%np
          cg(ip, idim) = sd(ip, idim) + gamma*cg(ip, idim)
        end do
      end do

      dot = X(mf_dotp)(mesh, dim, psi, cg)

      do idim = 1, st%d%dim
        do ip = 1, mesh%np
          cgp(ip, idim) = cg(ip, idim) - dot*psi(ip, idim)
        end do
      end do

      norm = X(mf_nrm2)(mesh, dim, cgp)

      call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, cgp, hcgp, ist, ik)

      niter = niter + 1

      alpha = -lambda + R_REAL(X(mf_dotp)(mesh, dim, cgp, hcgp))/norm**2
      beta  = M_TWO*R_REAL(X(mf_dotp)(mesh, dim, cgp, phi))/norm
      theta = M_HALF*atan(-beta/alpha)
      ctheta = cos(theta)
      stheta = sin(theta)

      ! This checks whether we are picking the maximum or the minimum.
      theta2 = theta + M_PI/M_TWO
      ctheta2 = cos(theta2)
      stheta2 = sin(theta2)
      sol(1) = lambda + stheta**2*alpha + beta*stheta*ctheta
      sol(2) = lambda + stheta2**2*alpha + beta*stheta2*ctheta2

      if(sol(2) < sol(1)) then
        theta = theta2
        stheta = stheta2
        ctheta = ctheta2
      end if
      stheta = stheta/norm

      do idim = 1, st%d%dim
        do ip = 1, mesh%np
          psi(ip, idim) = ctheta*psi(ip, idim) + stheta*cgp(ip, idim)
        end do
      end do

    end do band

    call states_elec_set_state(st, mesh, ist, ik, psi)

    st%eigenval(ist, ik) = lambda

    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      call loct_progress_bar(st%lnst*(ik - 1) +  ist, st%lnst*st%d%kpt%nlocal)
    end if

  end do states

  converged = conv

  SAFE_DEALLOCATE_A(phi)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(cg)
  SAFE_DEALLOCATE_A(hcgp)
  SAFE_DEALLOCATE_A(sd)
  SAFE_DEALLOCATE_A(cgp)

  POP_SUB(X(eigensolver_cg_jiang))
end subroutine X(eigensolver_cg_jiang)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
