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
!> calculates the eigenvalues of the orbitals
subroutine X(calculate_eigenvalues)(namespace, hm, der, st)
  type(namespace_t),        intent(in)    :: namespace
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(derivatives_t),      intent(in)    :: der
  type(states_elec_t),      intent(inout) :: st

  R_TYPE, allocatable :: eigen(:, :)

  PUSH_SUB(X(calculate_eigenvalues))

  if (debug%info) then
    write(message(1), '(a)') 'Debug: Calculating eigenvalues.'
    call messages_info(1, namespace=namespace)
  end if

  st%eigenval = M_ZERO

  SAFE_ALLOCATE(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
  call X(calculate_expectation_values)(namespace, hm, der, st, eigen)

  st%eigenval(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end) = &
    TOFLOAT(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

  call comm_allreduce(st%st_kpt_mpi_grp, st%eigenval)

  SAFE_DEALLOCATE_A(eigen)

  POP_SUB(X(calculate_eigenvalues))
end subroutine X(calculate_eigenvalues)

subroutine X(calculate_expectation_values)(namespace, hm, der, st, eigen, terms)
  type(namespace_t),        intent(in)    :: namespace
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(derivatives_t),      intent(in)    :: der
  type(states_elec_t),      intent(inout) :: st
  R_TYPE,                   intent(out)   :: eigen(st%st_start:, st%d%kpt%start:) !< (:st%st_end, :st%d%kpt%end)
  integer, optional,        intent(in)    :: terms

  integer :: ik, minst, maxst, ib
  type(wfs_elec_t) :: hpsib
  type(profile_t), save :: prof

  PUSH_SUB(X(calculate_expectation_values))

  call profiling_in(prof, TOSTRING(X(EIGENVALUE_CALC)))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = st%group%block_start, st%group%block_end

      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)

      if (hm%apply_packed()) then
        call st%group%psib(ib, ik)%do_pack()
      end if

      call st%group%psib(ib, ik)%copy_to(hpsib)

      call X(hamiltonian_elec_apply_batch)(hm, namespace, der%mesh, st%group%psib(ib, ik), hpsib, terms = terms)
      call X(mesh_batch_dotp_vector)(der%mesh, st%group%psib(ib, ik), hpsib, eigen(minst:maxst, ik), reduce = .false.)

      if (hm%apply_packed()) then
        call st%group%psib(ib, ik)%do_unpack(copy = .false.)
      end if

      call hpsib%end()
    end do
  end do

  if (der%mesh%parallel_in_domains) call der%mesh%allreduce(&
    eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

  call profiling_out(prof)
  POP_SUB(X(calculate_expectation_values))
end subroutine X(calculate_expectation_values)

! ---------------------------------------------------------
FLOAT function X(energy_calc_electronic)(namespace, hm, der, st, terms) result(energy)
  type(namespace_t),        intent(in)    :: namespace
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(derivatives_t),      intent(in)    :: der
  type(states_elec_t),      intent(inout) :: st
  integer,                  intent(in)    :: terms

  R_TYPE, allocatable  :: tt(:, :)

  PUSH_SUB(X(energy_calc_electronic))

  SAFE_ALLOCATE(tt(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

  call X(calculate_expectation_values)(namespace, hm, der, st, tt, terms = terms)

  energy = states_elec_eigenvalues_sum(st, TOFLOAT(tt))

  SAFE_DEALLOCATE_A(tt)
  POP_SUB(X(energy_calc_electronic))
end function X(energy_calc_electronic)

subroutine X(calculate_expectation_values_matrix)(namespace, hm, der, st, eigen, terms, diagonal_states)
  type(namespace_t),        intent(in)    :: namespace
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(derivatives_t),      intent(in)    :: der
  type(states_elec_t),      intent(inout) :: st
  R_TYPE,                   intent(out)   :: eigen(st%st_start:,st%st_start:,st%d%kpt%start:) !:st%st_end,:st%st_end,:st%d%kpt%end
  integer, optional,        intent(in)    :: terms
  logical, optional,        intent(in)    :: diagonal_states ! if true (default), computes only <i|H|i>, otherwise <i|H|j>

  integer :: ik, ib, jb
  type(wfs_elec_t) :: hpsib
  type(profile_t), save :: prof

  PUSH_SUB(X(calculate_expectation_values_matrix))

  call profiling_in(prof, TOSTRING(X(calc_exp_values_matrix)))

  eigen = M_ZERO


  do ik = st%d%kpt%start, st%d%kpt%end
    ! States are groupped in blocks. The dimension of each block depends on the hardware.
    ! The following do loop goes through all the groups of states (the parallelization here is DOMAINS or KPT)
    do ib = st%group%block_start, st%group%block_end

      if (hm%apply_packed()) call st%group%psib(ib, ik)%do_pack()

      call st%group%psib(ib, ik)%copy_to(hpsib)

      ! Perform H|psi_i> of the group of states
      call X(hamiltonian_elec_apply_batch)(hm, namespace, der%mesh, st%group%psib(ib, ik), hpsib, terms = terms)

      ! Now loop again over the groups of states, paying attention to avoid the combinations <i|H|j> == conj(<j|H|i>)
      do jb = ib, st%group%block_end
        if (optional_default(diagonal_states, .true.) .and. jb > ib) exit

        ! Perform the operation <psi_j|hpsi_i>, where hpsi_i = H|psi_i>. Normally it would be a vector multiplication,
        ! but since we are considering a block of states it is actually a matrix multiplication
        call X(mesh_batch_dotp_matrix)(der%mesh, st%group%psib(jb, ik), hpsib, eigen(:, :, ik), reduce = .false.)
        ! Remark: despite it seems that we are filling the upper part of the matrix,
        ! mesh_batch_dotp_matrix is actually computing only the lower triangular part
      end do

      if (hm%apply_packed()) call st%group%psib(ib, ik)%do_unpack(copy = .false.)
      call hpsib%end()
    end do
  end do

  if (der%mesh%parallel_in_domains) then
    call der%mesh%allreduce(eigen(st%st_start:st%st_end, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
  end if

  call profiling_out(prof)
  POP_SUB(X(calculate_expectation_values_matrix))
end subroutine X(calculate_expectation_values_matrix)

! ---------------------------------------------------------
subroutine X(one_body_matrix_elements)(st, namespace, gr, hm, nint, iindex, jindex, oneint)
  type(states_elec_t), intent(inout) :: st
  type(namespace_t),   intent(in)    :: namespace
  type(grid_t),        intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)  :: hm
  integer,             intent(in)    :: nint
  integer,             intent(out)   :: iindex(:)
  integer,             intent(out)   :: jindex(:)
  R_TYPE,              intent(out)   :: oneint(:)

  integer :: ist, jst, np, iint
  R_TYPE  :: me
  R_TYPE, allocatable :: psii(:, :), psij(:, :)
  R_TYPE, allocatable :: one_bodies(:, :, :)

  PUSH_SUB(X(one_body_matrix_elements))

  SAFE_ALLOCATE(psii(1:gr%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%np_part, 1:st%d%dim))

  if (st%d%ispin == SPINORS) then
    call messages_not_implemented("One-body integrals with spinors", namespace=namespace)
  end if
  ASSERT(.not. st%parallel_in_states) ! TODO: Add states parallelization support

  iint = 1

  ! Select case on theory level
  select case(hm%theory_level)
  case (HARTREE_FOCK)
    ! For Hartree Fock, the one body me should include only the kinetic energy and the external potential
    SAFE_ALLOCATE(one_bodies(st%st_start:st%st_end, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

    call X(calculate_expectation_values_matrix)(namespace, hm, gr%der, st, one_bodies, &
      terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL, diagonal_states = .false.)

    ! Now we have to fill the arrays iindex/jindex/oneint
    do ist = st%st_start, st%st_end
      do jst = st%st_start, ist
        iindex(iint) = ist
        jindex(iint) = jst
        oneint(iint) = one_bodies(ist, jst, 1)
        iint = iint + 1
      end do
    end do
    SAFE_DEALLOCATE_A(one_bodies)

  case (KOHN_SHAM_DFT)
    np = gr%np

    do ist = 1, st%nst

      call states_elec_get_state(st, gr, ist, 1, psii)

      do jst = 1, st%nst
        if (jst > ist) cycle

        call states_elec_get_state(st, gr, jst, 1, psij)

        psij(1:np, 1) = R_CONJ(psii(1:np, 1))*hm%vhxc(1:np, 1)*psij(1:np, 1)

        me = - X(mf_integrate)(gr, psij(:, 1))

        if (ist == jst) me = me + st%eigenval(ist,1)

        iindex(iint) = ist
        jindex(iint) = jst
        oneint(iint) = me
        iint = iint + 1
      end do
    end do
  case default
    call messages_not_implemented("One-body integrals with TheoryLevel not DFT or Hartree-Fock", namespace=namespace)
  end select

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)

  POP_SUB(X(one_body_matrix_elements))
end subroutine X(one_body_matrix_elements)

! ---------------------------------------------------------
subroutine X(two_body_matrix_elements) (st, namespace, space, gr, kpoints, psolver, st_min, st_max, iindex, &
  jindex, kindex, lindex, twoint, phase, singularity, exc_k)
  type(states_elec_t), target,   intent(inout) :: st
  type(namespace_t),             intent(in)    :: namespace
  type(space_t),                 intent(in)    :: space
  type(grid_t),                  intent(in)    :: gr
  type(kpoints_t),               intent(in)    :: kpoints
  type(poisson_t),               intent(inout) :: psolver
  integer,                       intent(in)    :: st_min, st_max
  integer,                       intent(out)   :: iindex(:,:)
  integer,                       intent(out)   :: jindex(:,:)
  integer,                       intent(out)   :: kindex(:,:)
  integer,                       intent(out)   :: lindex(:,:)
  R_TYPE,                        intent(out)   :: twoint(:)  !
  CMPLX,               optional, intent(in)    :: phase(:,st%d%kpt%start:)
  type(singularity_t), optional, intent(in)    :: singularity
  logical,             optional, intent(in)    :: exc_k

  integer :: ist, jst, kst, lst, ijst, klst, ikpt, jkpt, kkpt, lkpt
  integer :: ist_global, jst_global, kst_global, lst_global, nst, nst_tot
  integer :: iint, ikpoint, jkpoint, ip, ibind, npath
  R_TYPE  :: me
  R_TYPE, allocatable :: nn(:), vv(:), two_body_int(:), tmp(:)
  R_TYPE, pointer :: psii(:), psij(:), psik(:), psil(:)
  FLOAT :: qq(space%dim)
  logical :: exc_k_
  class(wfs_elec_t), pointer :: wfs
  type(fourier_space_op_t) :: coulb

  PUSH_SUB(X(two_body_matrix_elements))

  SAFE_ALLOCATE(nn(1:gr%np))
  SAFE_ALLOCATE(vv(1:gr%np))
  SAFE_ALLOCATE(tmp(1:gr%np))
  SAFE_ALLOCATE(two_body_int(1:gr%np))

  if (st%d%ispin == SPINORS) then
    call messages_not_implemented("Two-body integrals with spinors", namespace=namespace)
  end if

  ASSERT(present(phase) .eqv. present(singularity))
#ifdef R_TCOMPLEX
  ASSERT(present(phase))
#endif

  npath = kpoints%nkpt_in_path()

  if (st%are_packed()) call st%unpack()

  ijst = 0
  iint = 1

  nst_tot = (st_max-st_min+1)*st%d%nik
  nst = (st_max-st_min+1)

  exc_k_ = .false.
  if (present(exc_k)) exc_k_ = exc_k

  if (present(singularity)) then
    qq = M_ZERO
    call poisson_build_kernel(psolver, namespace, space, coulb, qq, M_ZERO)
  end if

  do ist_global = 1, nst_tot
    ist = mod(ist_global - 1, nst) + 1
    ikpt = (ist_global - ist) / nst + 1
    ikpoint = st%d%get_kpoint_index(ikpt)

    wfs => st%group%psib(st%group%iblock(ist+st_min-1, ikpt), ikpt)
    ASSERT(wfs%status() /= BATCH_DEVICE_PACKED)
    ibind = wfs%inv_index((/ist+st_min-1, 1/))
    if (wfs%status() == BATCH_NOT_PACKED) then
      psii => wfs%X(ff_linear)(:, ibind)
    else if (wfs%status() == BATCH_PACKED) then
      psii => wfs%X(ff_pack)(ibind, :)
    else
      ASSERT(.false.) ! TODO: Add GPU support
    end if

    do jst_global = 1, nst_tot
      jst = mod(jst_global - 1, nst) + 1
      jkpt = (jst_global - jst) / nst + 1
      jkpoint = st%d%get_kpoint_index(jkpt)

      if (exc_k_ .and. ist /= jst) cycle

      if (present(singularity)) then
        qq(:) = kpoints%get_point(ikpoint, absolute_coordinates=.false.) &
          - kpoints%get_point(jkpoint, absolute_coordinates=.false.)
        ! In case of k-points, the poisson solver must contains k-q
        ! in the Coulomb potential, and must be changed for each q point
        call poisson_build_kernel(psolver, namespace, space, coulb, qq, M_ZERO, &
          -(kpoints%full%npoints-npath)*kpoints%latt%rcell_volume*(singularity%Fk(jkpoint)-singularity%FF))
      end if

#ifndef R_TCOMPLEX
      if (jst_global > ist_global) cycle
#endif
      ijst=ijst+1

      wfs => st%group%psib(st%group%iblock(jst+st_min-1, jkpt), jkpt)
      ibind = wfs%inv_index((/jst+st_min-1, 1/))
      if (wfs%status() == BATCH_NOT_PACKED) then
        psij => wfs%X(ff_linear)(:, ibind)
      else if (wfs%status() == BATCH_PACKED) then
        psij => wfs%X(ff_pack)(ibind, :)
      else
        ASSERT(.false.) ! TODO: Add GPU support
      end if

      nn(1:gr%np) = R_CONJ(psii(1:gr%np))*psij(1:gr%np)
      if (present(singularity)) then
        call X(poisson_solve)(psolver, namespace, vv, nn, all_nodes=.false., kernel=coulb)
      else
        call X(poisson_solve)(psolver, namespace, vv, nn, all_nodes=.false.)
      end if

      !We now put back the phase that we treated analytically using the Poisson solver
#ifdef R_TCOMPLEX
      do ip = 1, gr%np
        vv(ip) = vv(ip) * exp(M_zI*sum(qq(:)*gr%x(ip, :)))
      end do
#endif

      klst=0
      do kst_global = 1, nst_tot
        kst = mod(kst_global - 1, nst) + 1
        kkpt = (kst_global - kst) / nst + 1

        if (exc_k_ .and. kkpt /= jkpt) cycle

        wfs => st%group%psib(st%group%iblock(kst + st_min - 1, kkpt), kkpt)
        ibind = wfs%inv_index((/kst + st_min - 1, 1/))
        if (wfs%status() == BATCH_NOT_PACKED) then
          psik => wfs%X(ff_linear)(:, ibind)
        else if (wfs%status() == BATCH_PACKED) then
          psik => wfs%X(ff_pack)(ibind, :)
        else
          ASSERT(.false.) ! TODO: Add GPU support
        end if

        if (present(phase)) then
#ifdef R_TCOMPLEX
          !$omp parallel do
          do ip = 1, gr%np
            tmp(ip) = vv(ip) * R_CONJ(psik(ip) * phase(ip, kkpt))
          end do
          !$omp end parallel do
#endif
        else
          !$omp parallel do
          do ip = 1, gr%np
            tmp(ip) = vv(ip)*R_CONJ(psik(ip))
          end do
          !$omp end parallel do
        end if

        do lst_global = 1, nst_tot
          lst = mod(lst_global - 1, nst) + 1
          lkpt = (lst_global - lst)/nst + 1

#ifndef R_TCOMPLEX
          if (lst_global > kst_global) cycle
          klst=klst+1
          if (klst > ijst) cycle
#endif

          if (exc_k_ .and. kst /= lst) cycle
          if (exc_k_ .and. lkpt /= ikpt) cycle
          wfs => st%group%psib(st%group%iblock(lst+st_min-1, lkpt), lkpt)
          ibind = wfs%inv_index((/lst+st_min-1, 1/))
          if (wfs%status() == BATCH_NOT_PACKED) then
            psil => wfs%X(ff_linear)(:, ibind)
          else if (wfs%status() == BATCH_PACKED) then
            psil => wfs%X(ff_pack)(ibind, :)
          else
            ASSERT(.false.) ! TODO: Add GPU support
          end if

          if (present(phase)) then
#ifdef R_TCOMPLEX
            !$omp parallel do
            do ip = 1, gr%np
              two_body_int(ip) = tmp(ip)*psil(ip)*phase(ip, lkpt)
            end do
            !$omp end parallel do
#endif
          else
            !$omp parallel do
            do ip = 1, gr%np
              two_body_int(ip) = tmp(ip)*psil(ip)
            end do
            !$omp end parallel do
          end if

          me = X(mf_integrate)(gr, two_body_int(:), reduce = .false.)

          iindex(1,iint) = ist + st_min - 1
          iindex(2,iint) = ikpt
          jindex(1,iint) = jst + st_min - 1
          jindex(2,iint) = jkpt
          kindex(1,iint) = kst + st_min - 1
          kindex(2,iint) = kkpt
          lindex(1,iint) = lst + st_min - 1
          lindex(2,iint) = lkpt
          twoint(iint) = me
          iint = iint + 1

        end do
      end do
    end do
  end do

  if (gr%parallel_in_domains) then
    call gr%allreduce(twoint)
  end if

  if (present(singularity)) then
    call fourier_space_op_end(coulb)
  end if

  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(two_body_int)

  POP_SUB(X(two_body_matrix_elements))
end subroutine X(two_body_matrix_elements)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
