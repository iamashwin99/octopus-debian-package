!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012-2013 M. Gruning, P. Melo, M. Oliveira
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
!!
! ---------------------------------------------------------
!> This file handles the evaluation of the OEP potential, in the KLI or full OEP
!! as described in S. Kuemmel and J. Perdew, PRL 90, 043004 (2003)
!!
!! This file has to be outside the module xc, for it requires the Hpsi.
!! This is why it needs the xc_functional module. I prefer to put it here since
!! the rest of the Hamiltonian module does not know about the gory details
!! of how xc is defined and calculated.
subroutine X(xc_oep_calc)(oep, namespace, xcs, gr, hm, st, space, rcell_volume, ex, ec, vxc)
  type(xc_oep_t),           intent(inout) :: oep
  type(namespace_t),        intent(in)    :: namespace
  type(xc_t),               intent(inout) :: xcs
  type(grid_t),             intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(inout) :: st
  type(space_t),            intent(in)    :: space
  FLOAT,                    intent(in)    :: rcell_volume
  FLOAT,                    intent(inout) :: ex, ec
  FLOAT, optional,          intent(inout) :: vxc(:,:) !< vxc(mesh%np, st%d%nspin)

  FLOAT :: eig
  integer :: is, ist, ixc, nspin_, isp, idm, ib, ik
  logical, save :: first = .true.
  R_TYPE, allocatable :: psi(:), xpsi(:)
  type(states_elec_t) :: xst
  logical :: exx

  if (oep%level == OEP_LEVEL_NONE) return

  call profiling_in(C_PROFILING_XC_OEP, TOSTRING(X(XC_OEP)))
  PUSH_SUB(X(xc_oep_calc))

  ! initialize oep structure
  nspin_ = min(st%d%nspin, 2)
  SAFE_ALLOCATE(oep%eigen_type (1:st%nst))
  SAFE_ALLOCATE(oep%eigen_index(1:st%nst))
  SAFE_ALLOCATE(oep%X(lxc)(1:gr%np, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
  oep%X(lxc) = M_ZERO
  SAFE_ALLOCATE(oep%uxc_bar(1:st%d%dim, 1:st%nst, st%d%kpt%start:st%d%kpt%end))

  !We first apply the exchange operator to all the states
  call xst%nullify()

  exx = .false.
  select case(oep%type)
  case(OEP_TYPE_EXX, OEP_TYPE_PHOTONS)
    functl_loop: do ixc = 1, 2
      if (xcs%functional(ixc, 1)%family /= XC_FAMILY_OEP) cycle
      select case (xcs%functional(ixc,1)%id)
      case (XC_OEP_X)
        call X(exchange_operator_compute_potentials)(hm%exxop, namespace, space, gr, st, xst, hm%kpoints, eig)
        ex = ex + eig
        exx = .true.
      end select
    end do functl_loop

    ! MGGA from OEP
  case(OEP_TYPE_MGGA)
    if (.not. xst%group%block_initialized) then
      call states_elec_copy(xst, st)
      call states_elec_set_zero(xst)
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        ! Here we apply the term by calling hamiltonian_elec_apply_batch, which takes care of setting the
        ! phases properly, as well as upating boundary points/ghost points when needed
        call X(hamiltonian_elec_apply_batch)(hm, namespace, gr, st%group%psib(ib, ik), &
          xst%group%psib(ib, ik), terms = TERM_MGGA)
      end do
    end do
    exx = .true.

    ! this part handles the (pure) orbital functionals
    ! SIC a la PZ is handled here
  case(OEP_TYPE_SIC)
    ! this routine is only prepared for finite systems. (Why not?)
    if (st%d%nik > st%d%ispin) then
      call messages_not_implemented("OEP-SIC for periodic systems", namespace=namespace)
    end if
    spin: do is = 1, nspin_
      ! distinguish between 'is' being the spin_channel index (collinear)
      ! and being the spinor (noncollinear)
      if (st%d%ispin == SPINORS) then
        isp = 1
        idm = is
      else
        isp = is
        idm = 1
      end if

      call X(oep_sic) (xcs, gr, hm%psolver, namespace, space, rcell_volume, st, hm%kpoints, is, oep, ex, ec)
    end do spin

  case default
    ASSERT(.false.)
  end select

  ! calculate uxc_bar for the occupied states

  SAFE_ALLOCATE(psi(1:gr%np))
  SAFE_ALLOCATE(xpsi(1:gr%np))

  oep%uxc_bar(:, :, :) = M_ZERO
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idm = 1, st%d%dim
        call states_elec_get_state(st, gr, idm, ist, ik, psi)
        if (exx) then
          ! Here we copy the state from xst to X(lxc).
          ! This will be removed in the future, but it allows to keep both EXX and PZ-SIC in the code
          call states_elec_get_state(xst, gr, idm, ist, ik, xpsi)
          ! There is a complex conjugate here, as the lxc is defined as <\psi|X and
          ! exchange_operator_compute_potentials returns X|\psi>
#ifndef R_TREAL
          xpsi = R_CONJ(xpsi)
#endif
          call lalg_axpy(gr%np, M_ONE, xpsi, oep%X(lxc)(1:gr%np, idm, ist, ik))
        end if
        oep%uxc_bar(idm, ist, ik) = R_REAL(X(mf_dotp)(gr, psi, oep%X(lxc)(1:gr%np, idm, ist, ik), reduce = .false., dotu = .true.))
      end do
    end do
  end do
  if (gr%parallel_in_domains) then
    call gr%allreduce(oep%uxc_bar(1:st%d%dim, 1:st%st_end, st%d%kpt%start:st%d%kpt%end))
  end if

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(xpsi)

  call states_elec_end(xst)

  if (st%parallel_in_states) then
    call st%mpi_grp%barrier()
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = 1, st%nst
        call st%mpi_grp%bcast(oep%uxc_bar(1, ist, ik), st%d%dim, MPI_FLOAT, st%node(ist))
      end do
    end do
  end if

  if (st%d%ispin == SPINORS) then
    call xc_oep_AnalyzeEigen(oep, st, 1)
    call xc_KLI_Pauli_solve(gr, st, oep)
    if (present(vxc)) then
      vxc(1:gr%np, 1:4) = vxc(1:gr%np,1:4) + oep%vxc(1:gr%np,1:4)
    end if
    ! full OEP not implemented!
  else
    if (oep%level == OEP_LEVEL_FULL .and. (.not. first)) then
      do is = 1, nspin_
        ! get the HOMO state
        call xc_oep_AnalyzeEigen(oep, st, is)

        if (present(vxc)) then
          if (oep%type == OEP_TYPE_PHOTONS) then
            call X(xc_oep_solve_photon)(namespace, gr, hm, st, is, vxc(:,is), oep)
          else
            call X(xc_oep_solve)(namespace, gr, hm, st, is, vxc(:,is), oep)
          end if

          call lalg_axpy(gr%np, M_ONE, oep%vxc(1:gr%np, is), vxc(1:gr%np, is))
        end if
      end do

    else !KLI
      if (oep%type == OEP_TYPE_PHOTONS) then
        oep%vxc = M_ZERO
        do is = 1, nspin_
          call xc_oep_AnalyzeEigen(oep, st, is)
          call X(xc_KLI_solve_photon) (namespace, gr, hm, st, is, oep, first)
          call lalg_axpy(gr%np, M_ONE, oep%vxc(1:gr%np, is), vxc(1:gr%np, is))
        end do
      else
        call X(xc_KLI_solve) (space, gr, st, oep, rcell_volume)
        call lalg_axpy(gr%np, st%d%nspin, M_ONE, oep%vxc, vxc)
      end if

    end if
    first = .false.
  end if
  SAFE_DEALLOCATE_A(oep%eigen_type)
  SAFE_DEALLOCATE_A(oep%eigen_index)
  SAFE_DEALLOCATE_A(oep%X(lxc))
  SAFE_DEALLOCATE_A(oep%uxc_bar)

  POP_SUB(X(xc_oep_calc))
  call profiling_out(C_PROFILING_XC_OEP)
end subroutine X(xc_OEP_calc)


! ---------------------------------------------------------
!> This routine follows closely the one of PRB 68, 035103 (2003)
!> Below we refer to the equation number of this paper
subroutine X(xc_oep_solve) (namespace, mesh, hm, st, is, vxc, oep)
  type(namespace_t),        intent(in)    :: namespace
  class(mesh_t),            intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer,                  intent(in)    :: is
  FLOAT,                    intent(inout) :: vxc(:) !< (mesh%np, given for the spin is)
  type(xc_oep_t),           intent(inout) :: oep

  integer :: iter, ist, iter_used
  FLOAT :: vxc_bar, ff, residue
  FLOAT, allocatable :: ss(:), vxc_old(:)
  FLOAT, allocatable :: psi2(:,:)
  R_TYPE, allocatable :: bb(:,:), psi(:, :)
  logical, allocatable :: orthogonal(:)

  call profiling_in(C_PROFILING_OEP_LEVEL_FULL, TOSTRING(X(OEP_LEVEL_FULL)))
  PUSH_SUB(X(xc_oep_solve))

  if (st%parallel_in_states) then
    call messages_not_implemented("Full OEP parallel in states", namespace=namespace)
  end if

  SAFE_ALLOCATE(     bb(1:mesh%np, 1))
  SAFE_ALLOCATE(     ss(1:mesh%np))
  SAFE_ALLOCATE(vxc_old(1:mesh%np))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(orthogonal(1:st%nst))

  call lalg_copy(mesh%np, vxc, vxc_old)

  if (.not. lr_is_allocated(oep%lr)) then
    call lr_allocate(oep%lr, st, mesh)
    oep%lr%X(dl_psi)(:,:, :, :) = M_ZERO
  end if

  ! fix xc potential (needed for Hpsi)
  call lalg_axpy(mesh%np, M_ONE, vxc_old(:), oep%vxc(:, is))

  do iter = 1, oep%scftol%max_iter
    ! iteration over all states
    ss = M_ZERO
    do ist = 1, st%nst

      if (abs(st%occ(ist,is)) <= M_EPSILON) cycle  !only over occupied states
      call states_elec_get_state(st, mesh, ist, is, psi)
      psi2(:, 1) = TOFLOAT(R_CONJ(psi(:, 1))*psi(:,1))

      ! evaluate right-hand side
      vxc_bar = dmf_integrate(mesh, psi2(:, 1)*oep%vxc(1:mesh%np, is))

      ! This the right-hand side of Eq. 21
      bb(1:mesh%np, 1) = -(oep%vxc(1:mesh%np, is) - (vxc_bar - oep%uxc_bar(1, ist, is)))* &
        R_CONJ(psi(:, 1)) + oep%X(lxc)(1:mesh%np, 1, ist, is)

      call X(lr_orth_vector) (mesh, st, bb, ist, is, R_TOTYPE(M_ZERO))

      ! Sternheimer equation [H-E_i]psi_i = bb_i, where psi_i the orbital shift, see Eq. 21
      call X(linear_solver_solve_HXeY)(oep%solver, namespace, hm, mesh, st, ist, is, &
        oep%lr%X(dl_psi)(:,:, ist, is), bb, &
        R_TOTYPE(-st%eigenval(ist, is)), oep%scftol%final_tol, residue, iter_used)

      if (debug%info) then
        write(message(1),'(a,i3,a,es14.6,a,es14.6,a,i4)') "Debug: OEP - iter ", iter, &
          " linear solver residual ", residue, " tol ", &
          oep%scftol%final_tol, " iter_used ", iter_used
        call messages_info(1, namespace=namespace)
      end if

      !We project out the occupied bands
      call X(lr_orth_vector) (mesh, st, oep%lr%X(dl_psi)(:,:, ist, is), ist, is, R_TOTYPE(M_ZERO))

      ! calculate this funny function ss
      ! ss = ss + 2*dl_psi*psi
      ! This is Eq. 25
      call lalg_axpy(mesh%np, M_TWO, R_REAL(oep%lr%X(dl_psi)(1:mesh%np, 1, ist, is)*psi(:, 1)), ss(:))
    end do

    !Here we mix the xc potential
    call X(xc_oep_mix)(oep, mesh, ss, st%rho(:,is), is)

    !Here we enforce Eq. (24), see the discussion below Eq. 26
    do ist = 1, st%nst
      if (oep%eigen_type(ist) == 2) then
        call states_elec_get_state(st, mesh, ist, is, psi)
        psi2(:, 1) = TOFLOAT(R_CONJ(psi(:, 1))*psi(:,1))
        vxc_bar = dmf_integrate(mesh, psi2(:, 1)*oep%vxc(1:mesh%np, is))
        oep%vxc(1:mesh%np,is) = oep%vxc(1:mesh%np,is) - (vxc_bar - oep%uxc_bar(1, ist,is))
      end if
    end do

    ff = dmf_nrm2(mesh, ss)

    if (debug%info) then
      write(message(1),'(a,i3,a,es14.6,a,i4)') "Debug: OEP - iter ", iter, " residual ", ff, " max ", oep%scftol%max_iter
      call messages_info(1, namespace=namespace)
    end if

    if (ff < oep%scftol%conv_abs_dens) exit
  end do

  if (is == 1) then
    oep%norm2ss = ff
  else
    oep%norm2ss = oep%norm2ss + ff !adding up spin up and spin down component
  end if

  if (ff > oep%scftol%conv_abs_dens) then
    write(message(1), '(a)') "OEP did not converge."
    call messages_warning(1, namespace=namespace)

    ! otherwise the number below will be one too high
    iter = iter - 1
  end if

  write(message(1), '(a,i4,a,es14.6)') "Info: After ", iter, " iterations, the OEP residual = ", ff
  message(2) = ''
  call messages_info(2, namespace=namespace)

  call lalg_copy(mesh%np, vxc_old, vxc)

  SAFE_DEALLOCATE_A(bb)
  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(vxc_old)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi2)
  SAFE_DEALLOCATE_A(orthogonal)

  POP_SUB(X(xc_oep_solve))
  call profiling_out(C_PROFILING_OEP_LEVEL_FULL)
end subroutine X(xc_oep_solve)

! ---------------------------------------------------------
!> This is the photon version of the xc_oep_solve routine
subroutine X(xc_oep_solve_photon) (namespace, mesh, hm, st, is, vxc, oep)
  type(namespace_t),        intent(in)    :: namespace
  class(mesh_t),            intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer,                  intent(in)    :: is
  FLOAT,                    intent(inout) :: vxc(:) !< (mesh%np, given for the spin is)
  type(xc_oep_t),           intent(inout) :: oep

  integer :: iter, ist, iter_used
  FLOAT :: vxc_bar, ff, residue
  FLOAT, allocatable :: ss(:), vxc_old(:)
  FLOAT, allocatable :: psi2(:,:)
  R_TYPE, allocatable :: bb(:,:), psi(:, :), phi1(:,:,:)
  logical, allocatable :: orthogonal(:)

  call profiling_in(C_PROFILING_OEP_LEVEL_FULL, TOSTRING(X(OEP_LEVEL_FULL_PHOTON)))
  PUSH_SUB(X(xc_oep_solve_photon))

  if (st%parallel_in_states) then
    call messages_not_implemented("Full OEP parallel in states", namespace=namespace)
  end if

#ifndef R_TREAL
  ! Photons with OEP are only implemented for real states
  ASSERT(.false.)
#endif

  SAFE_ALLOCATE(     bb(1:mesh%np, 1))
  SAFE_ALLOCATE(     ss(1:mesh%np))
  SAFE_ALLOCATE(vxc_old(1:mesh%np))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(orthogonal(1:st%nst))
  SAFE_ALLOCATE(phi1(1:mesh%np, 1:st%d%dim, 1:oep%noccst))

  call lalg_copy(mesh%np, vxc, vxc_old)

  if (.not. lr_is_allocated(oep%lr)) then
    call lr_allocate(oep%lr, st, mesh)
    oep%lr%X(dl_psi)(:,:, :, :) = M_ZERO
  end if

  if (.not. lr_is_allocated(oep%photon_lr)) then
    call lr_allocate(oep%photon_lr, st, mesh)
    oep%photon_lr%ddl_psi(:, :, :, :) = M_ZERO
  end if
#ifdef R_TREAL
  call xc_oep_pt_phi(namespace, mesh, hm, st, is, oep, phi1)
#endif

  ! fix xc potential (needed for Hpsi)
  call lalg_axpy(mesh%np, M_ONE, vxc_old(:), oep%vxc(:, is))

  do iter = 1, oep%scftol%max_iter
    ! iteration over all states
    ss = M_ZERO
    do ist = 1, st%nst

      if (abs(st%occ(ist,is))<= M_EPSILON) cycle  !only over occupied states
      call states_elec_get_state(st, mesh, ist, is, psi)
      psi2(:, 1) = TOFLOAT(R_CONJ(psi(:, 1))*psi(:,1))

      ! evaluate right-hand side
      vxc_bar = dmf_integrate(mesh, psi2(:, 1)*oep%vxc(1:mesh%np, is))

      ! This the right-hand side of Eq. 21
      bb(1:mesh%np, 1) = -(oep%vxc(1:mesh%np, is) - (vxc_bar - oep%uxc_bar(1, ist, is)))* &
        R_CONJ(psi(:, 1)) + oep%X(lxc)(1:mesh%np, 1, ist, is)

#ifdef R_TREAL
      call xc_oep_pt_rhs(mesh, st, is, oep, phi1, ist, bb)
#endif

      orthogonal = .true.
      orthogonal(ist) = .false.
      call X(states_elec_orthogonalize_single)(st, mesh, st%nst, is, bb, normalize = .false., mask = orthogonal)

      ! Sternheimer equation [H-E_i]psi_i = bb_i, where psi_i the orbital shift, see Eq. 21
      call X(linear_solver_solve_HXeY)(oep%solver, namespace, hm, mesh, st, ist, is, &
        oep%lr%X(dl_psi)(:,:, ist, is), bb, &
        R_TOTYPE(-st%eigenval(ist, is)), oep%scftol%final_tol, residue, iter_used)

      if (debug%info) then
        write(message(1),'(a,i3,a,es14.6,a,es14.6,a,i4)') "Debug: OEP - iter ", iter, &
          " linear solver residual ", residue, " tol ", &
          oep%scftol%final_tol, " iter_used ", iter_used
        call messages_info(1, namespace=namespace)
      end if

      !We project out the occupied bands
      orthogonal = .true.
      orthogonal(ist) = .false.
      call X(states_elec_orthogonalize_single)(st, mesh, st%nst, is, &
        oep%lr%X(dl_psi)(:,:, ist, is), normalize = .false., mask = orthogonal)

      ! calculate this funny function ss
      ! ss = ss + 2*dl_psi*psi
      ! This is Eq. 25
      call lalg_axpy(mesh%np, M_TWO, R_REAL(oep%lr%X(dl_psi)(1:mesh%np, 1, ist, is)*psi(:, 1)), ss(:))

#ifdef R_TREAL
      call xc_oep_pt_inhomog(mesh, st, is, phi1, ist, ss)
#endif
    end do

    !Here we mix the xc potential
    call X(xc_oep_mix)(oep, mesh, ss, st%rho(:,is), is)

    !Here we enforce Eq. (24), see the discussion below Eq. 26
    do ist = 1, st%nst
      if (oep%eigen_type(ist) == 2) then
        call states_elec_get_state(st, mesh, ist, is, psi)
        psi2(:, 1) = TOFLOAT(R_CONJ(psi(:, 1))*psi(:,1))
        vxc_bar = dmf_integrate(mesh, psi2(:, 1)*oep%vxc(1:mesh%np, is))
#ifdef R_TREAL
        call xc_oep_pt_uxcbar(mesh, st, is, oep, phi1, ist, vxc_bar)
#endif
        oep%vxc(1:mesh%np,is) = oep%vxc(1:mesh%np,is) - (vxc_bar - oep%uxc_bar(1, ist,is))
      end if
    end do

    ff = dmf_nrm2(mesh, ss)

    if (debug%info) then
      write(message(1),'(a,i3,a,es14.6,a,i4)') "Debug: OEP - iter ", iter, " residual ", ff, " max ", oep%scftol%max_iter
      call messages_info(1, namespace=namespace)
    end if

    if (ff < oep%scftol%conv_abs_dens) exit
  end do

  if (is == 1) then
    oep%norm2ss = ff
  else
    oep%norm2ss = oep%norm2ss + ff !adding up spin up and spin down component
  end if

  if (ff > oep%scftol%conv_abs_dens) then
    write(message(1), '(a)') "OEP did not converge."
    call messages_warning(1, namespace=namespace)

    ! otherwise the number below will be one too high
    iter = iter - 1
  end if

  write(message(1), '(a,i4,a,es14.6)') "Info: After ", iter, " iterations, the OEP residual = ", ff
  message(2) = ''
  call messages_info(2, namespace=namespace)

  call lalg_copy(mesh%np, vxc_old, vxc)

  SAFE_DEALLOCATE_A(bb)
  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(vxc_old)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi2)
  SAFE_DEALLOCATE_A(orthogonal)
  SAFE_DEALLOCATE_A(phi1)

  POP_SUB(X(xc_oep_solve_photon))
  call profiling_out(C_PROFILING_OEP_LEVEL_FULL)
end subroutine X(xc_oep_solve_photon)

!----------------------------------------------------------------------
!> A routine that takes care of mixing the potential
subroutine X(xc_oep_mix)(oep, mesh, ss, rho, is)
  type(xc_oep_t),           intent(inout) :: oep
  type(mesh_t),             intent(in)    :: mesh
  FLOAT,                    intent(in)    :: ss(:)
  FLOAT,                    intent(in)    :: rho(:)
  integer,                  intent(in)    :: is

  integer :: ip
  FLOAT, allocatable :: mix(:)

  PUSH_SUB(X(xc_oep_mix))

  !Here we mix the xc potential
  select case (oep%mixing_scheme)
  case (OEP_MIXING_SCHEME_CONST)
    !This is Eq. 26
    call lalg_axpy(mesh%np, oep%mixing, ss(:), oep%vxc(:, is))

  case (OEP_MIXING_SCHEME_DENS)

    ! See Eq. 28 of the Kuemmel paper
    SAFE_ALLOCATE(mix(1:mesh%np))
    do ip = 1, mesh%np
      mix(ip) = - M_HALF * oep%vxc(ip, is) / (rho(ip) + M_TINY)
      ! To avoid nonsense local mixing, we put a maximum value here
      ! This tipically occurs when the wavefunctions are not well converged
      if (abs(mix(ip)) > CNST(1e3)) mix(ip) = sign(CNST(1e3),mix(ip))
      mix(ip) = mix(ip) * ss(ip)
    end do

    call lalg_axpy(mesh%np, M_ONE, mix, oep%vxc(:, is))
    SAFE_DEALLOCATE_A(mix)

  case (OEP_MIXING_SCHEME_BB)
    !This is the Barzilai-Borwein scheme, as explained in
    !Hollins, et al. PRB 85, 235126 (2012)
    if (dmf_nrm2(mesh, oep%vxc_old(1:mesh%np,is)) > M_EPSILON) then ! do not do it for the first run
      oep%mixing = -dmf_dotp(mesh, oep%vxc(1:mesh%np,is) - oep%vxc_old(1:mesh%np,is), ss - oep%ss_old(:, is)) &
        / dmf_dotp(mesh, ss - oep%ss_old(:, is), ss - oep%ss_old(:, is))
    end if

    if (debug%info) then
      write(message(1), '(a,es14.6,a,es15.8)') "Info: oep%mixing:", oep%mixing, " norm2ss: ", dmf_nrm2(mesh, ss)
      call messages_info(1)
    end if

    call lalg_copy(mesh%np, oep%vxc(:,is), oep%vxc_old(:,is))
    call lalg_copy(mesh%np, ss, oep%ss_old(:, is))
    call lalg_axpy(mesh%np, oep%mixing, ss(:), oep%vxc(:, is))

  end select

  POP_SUB(X(xc_oep_mix))
end subroutine X(xc_oep_mix)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
