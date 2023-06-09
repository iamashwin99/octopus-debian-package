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

!> This subroutine implements the preconditioned Lanczos eigensolver as
!! described in the paper:
!!
!! Y. Saad, A. Stathopoulos, J. Chelikowsky, K. Wu and S. Ogut,
!! "Solution of Large Eigenvalue Problems in Electronic Structure Calculations",
!! BIT 36 563-578 (1996) doi:10.1007/BF01731934 .
!!
!! We also implement the "smoothing" preconditioning described in that paper.
subroutine X(eigensolver_plan) (namespace, mesh, st, hm, pre, tol, niter, converged, ik, diff)
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(states_elec_t),         intent(inout) :: st
  type(hamiltonian_elec_t),    intent(in)    :: hm
  type(preconditioner_t),      intent(in)    :: pre
  FLOAT,                       intent(in)    :: tol
  integer,                     intent(inout) :: niter
  integer,                     intent(out)   :: converged
  integer,                     intent(in)    :: ik
  FLOAT,             optional, intent(out)   :: diff(:) !< (1:st%nst)

  !  integer :: n       ! Dimension of the problem.
  integer :: ned        ! Number of smallest eigenpairs desired
  integer :: nec        ! number of eigen-pairs converged, if initially
  ! nec > 0, the first nec elements of eigenval, res and
  ! first nec columns of eigenvec are assumed to have converged
  ! eigen-pairs and corresponding residual norms.
  integer :: maxmatvecs ! Maximum number of matrix-vectors applications allowed.
  ! On exit reset to actual number of MATVECs used
  integer :: me         ! array size of eigenval, res and number of columns in eigenvec.
  FLOAT,  allocatable :: eigenval(:)     ! The eigenvalues
  R_TYPE, allocatable :: eigenvec(:,:,:) ! The eigenvectors
  FLOAT,  allocatable :: res(:)          ! The residuals
  R_TYPE, allocatable :: vv(:,:,:)       ! The Krylov subspace basis vectors
  R_TYPE, allocatable :: av(:,:,:)       ! Workspace: W = A V
  FLOAT,  allocatable :: tmp(:)          ! Workspace.
  R_TYPE, allocatable :: ham(:,:)        ! Projection of the Hamiltonian onto Krylov subspace.
  R_TYPE, allocatable :: hevec(:,:)
  R_TYPE, allocatable :: aux(:,:)
  type(wfs_elec_t) :: vvb, avb
  integer  :: blk, ist, ii, dim, jst, d1, d2, matvec, nconv
  FLOAT :: xx

  ! Some hard-coded parameters.
  integer, parameter  :: winsiz = 5  ! window size, number of eigenvalues computed simultaneously
  integer, parameter  :: krylov = 15 ! The Krylov subspace size.
  integer, parameter  :: krylov_half = 7 ! Half the Krylov subspace size (rounded down).

  PUSH_SUB(X(eigensolver_plan))

  !  n          = mesh%np*st%d%dim
  dim        = st%d%dim
  ned        = st%nst
  nec        = 0
  maxmatvecs = niter*st%d%nik*st%nst
  me         = ned + winsiz - 1

  ! Allocate memory
  ! Careful: aux has to range from 1 to mesh%np_part because it is input to
  ! hpsi. In parallel the space NP+1:mesh%np_part is needed for ghost points
  ! in the non local operator.
  SAFE_ALLOCATE(eigenvec(1:mesh%np, 1:dim, 1:me))
  SAFE_ALLOCATE(aux(1:mesh%np_part, 1:dim))
  SAFE_ALLOCATE(tmp(1:krylov))
  SAFE_ALLOCATE(vv(1:mesh%np, 1:dim, 1:krylov))
  SAFE_ALLOCATE(ham(1:krylov, 1:krylov))
  SAFE_ALLOCATE(eigenval(1:me))
  SAFE_ALLOCATE(av(1:mesh%np, 1:dim, 1:krylov))
  SAFE_ALLOCATE(hevec(1:krylov, 1:krylov))
  SAFE_ALLOCATE(res(1:me))

  eigenval = M_ZERO
  eigenvec = R_TOTYPE(M_ZERO)
  res      = M_ZERO
  vv       = R_TOTYPE(M_ZERO)
  av       = R_TOTYPE(M_ZERO)
  tmp      = M_ZERO
  ham      = R_TOTYPE(M_ZERO)
  hevec    = R_TOTYPE(M_ZERO)
  aux      = R_TOTYPE(M_ZERO)

  niter = 0 ! Initialize the total matrix-vector multiplication counter.

  ! First of all, copy the initial estimates.
  do ist = 1, st%nst
    call states_elec_get_state(st, mesh, ist, ik, eigenvec(:, :, ist))
    eigenval(ist) = st%eigenval(ist, ik)
  end do

  ! Initialization of counters...
  matvec = 0 ! Set the matrix-multiplication counter to zero.
  nec    = 0 ! Sets the converged vectors counter to zero.
  d1     = 0 ! index for inner loop
  nconv  = 0 ! number of eigen-pairs converged.

  ! Sets the projected Hamiltonian matrix to zero.
  ham = R_TOTYPE(M_ZERO)

  ! Beginning of the outer loop; start/restart
  outer_loop : do

    if (nec >= ned)           exit outer_loop ! :)   Already converged!
    if (matvec >= maxmatvecs) exit outer_loop ! :(   Maximum number of mat-vec operation surpassed...

    if (d1 <= winsiz) then !start from beginning
      blk = winsiz
    else                   !restart to work on another set of eigen-pairs
      blk = min(krylov_half, d1)
    end if

    !copy next set of Ritz vector/initial guesses to vv
    do ist = 1, winsiz
      call lalg_copy(mesh%np, dim, eigenvec(:, :, nec+ist), vv(:, :, ist))
    end do

    ! Beginning of the inner loop.
    d1 = 0
    inner_loop: do
      d2 = d1 + blk

      ! Orthonormalization. The vectors in vv are orthonormalized against the converged
      ! eigenvectors, and amongst themselves. A check is done for the case of linear dependence,
      ! and random vectors are created in that case.
      ist = d1 + 1
      ortho: do
        if (ist>d2) exit ortho
        do ii = 1, nec
          av(ii, 1, d1 + 1) = X(mf_dotp)(mesh, dim, eigenvec(:, :, ii), vv(:, :, ist))
          call lalg_axpy(mesh%np, dim, -av(ii, 1, d1 + 1), eigenvec(:, :, ii), vv(:, :, ist))
        end do
        do ii = 1, ist - 1
          av(ii, 1, d1 + 1) = X(mf_dotp)(mesh, dim, vv(:, :, ii), vv(:, :, ist))
          call lalg_axpy(mesh%np, dim, -av(ii, 1, d1 + 1), vv(:, :, ii), vv(:, :, ist))
        end do
        xx = X(mf_nrm2)(mesh, dim, vv(:, :, ist))
        if (xx  <=  M_EPSILON) then
          if (st%randomization == PAR_INDEPENDENT) then
            call X(mf_random)(mesh, vv(:, 1, ist), &
              pre_shift = mesh%pv%xlocal-1, &
              post_shift = mesh%pv%np_global - mesh%pv%xlocal - mesh%np + 1)
          else
            call X(mf_random)(mesh, vv(:, 1, ist))
          end if
        else
          call lalg_scal(mesh%np, dim, R_TOTYPE(M_ONE/xx), vv(:, :, ist))
          ist = ist + 1
        end if
      end do ortho

      call X(wfs_elec_init)(vvb, st%d%dim, 1, blk, mesh%np_part, ik)
      call wfs_elec_init(avb, st%d%dim, d1 + 1, d1 + blk, av(:, :, d1 + 1:), ik)

      ! we need to copy to mesh%np_part size array
      do ist = 1, blk
        call batch_set_state(vvb, ist, mesh%np, vv(:, :, d1 + ist))
      end do

      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, vvb, avb)
      matvec = matvec + blk

      call vvb%end()
      call avb%end()

      ! Here we calculate the last blk columns of H = V^T A V. We do not need the lower
      ! part of the matrix since it is symmetric (LAPACK routine only needs the upper triangle).
      do ist = d1 + 1, d2
        do ii = 1, ist
          ham(ii, ist) = X(mf_dotp)(mesh, dim, vv(:, :, ii), av(:, :, ist))
        end do
      end do

      ! Diagonalization in the subspace, using LAPACK.
      hevec(1:d2, 1:d2) = ham(1:d2, 1:d2)
      call lalg_eigensolve(d2, hevec, tmp)

      ! Store the Ritz values as approximate eigenvalues.
      call lalg_copy(winsiz, tmp, eigenval(nec + 1:nec + winsiz))

      if (d2 + 1 <= krylov .and. matvec < maxmatvecs) then
        ! In this case, compute only the lowest Ritz eigenpair.
        call lalg_gemv(mesh%np, dim, d2, R_TOTYPE(M_ONE), vv(:, :, 1:d2), hevec(1:d2, 1), &
          R_TOTYPE(M_ZERO), eigenvec(:, :, nec + 1))
        call lalg_gemv(mesh%np, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, 1), &
          R_TOTYPE(M_ZERO), av(:, :, d2 + 1))
        call residual(av(:, :, d2+1), eigenvec(:, :, nec+1), tmp(1), av(:, :, d2+1), res(nec+1))

        ! If the first Ritz eigen-pair converged, compute all
        ! Ritz vectors and the residual norms.
        if (res(nec + 1) < tol) then
          do ist = 2, winsiz
            call lalg_gemv(mesh%np, dim, d2, R_TOTYPE(M_ONE), vv(:, :, 1:d2), hevec(1:d2, ist), &
              R_TOTYPE(M_ZERO), eigenvec(:, :, nec+ist))
          end do
          do ist = 2, winsiz
            call lalg_gemv(mesh%np, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, ist), &
              R_TOTYPE(M_ZERO), vv(:, :, ist))
          end do
          do ist = 2, winsiz
            call residual(vv(:, :, ist), eigenvec(:, :, nec+ist), tmp(ist), av(:, :, ist), res(nec+ist))
          end do
        end if
        d1 = d2
      else
        do ist = 1, winsiz
          call lalg_gemv(mesh%np, dim, d2, R_TOTYPE(M_ONE), vv(:, :, 1:d2), hevec(1:d2, ist), &
            R_TOTYPE(M_ZERO), eigenvec(:, :, nec+ist))
        end do
        do ist = 1, winsiz
          call lalg_gemv(mesh%np, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, ist), &
            R_TOTYPE(M_ZERO), vv(:, :, ist))
        end do
        do ist = 1, winsiz
          call lalg_copy(mesh%np, dim, vv(:, :, ist), av(:, :, ist))
          call lalg_copy(mesh%np, dim, eigenvec(:, :, nec + ist), vv(:, :, ist))
          call residual(av(:, :, ist), vv(:, :, ist), tmp(ist), av(:, :, winsiz+ist), res(nec+ist))
        end do

        ! Forms the first winsiz rows of H = V^T A V
        do ist = 1, winsiz
          do ii = 1, ist
            ham(ii, ist) = X(mf_dotp)(mesh, dim, vv(:, :, ii), av(:, :, ist))
          end do
        end do
        d1 = winsiz
      end if
      blk = 1

      ! Convergence test, and reordering of the eigenpairs. Starts checking
      ! the convergence of the eigenpairs of the window, and stops checking
      ! whenever it finds one not converged. Then, for each converged eigenpair,
      ! compares its eigenvalue to the previous one, swapping them if
      ! necessary.
      nconv = 0
      ordering: do ist = nec + 1, nec + winsiz - 1
        if (res(ist) >= tol) exit ordering
        nconv = nconv + 1
        do jst = ist, 2, -1
          if (eigenval(jst-1) <= eigenval(jst)) exit

          xx = eigenval(jst-1)
          eigenval(jst-1) = eigenval(jst)
          eigenval(jst) = xx

          xx = res(jst-1)
          res(jst-1) = res(jst)
          res(jst) = xx

          call lalg_swap(mesh%np, dim, eigenvec(:, :, jst), eigenvec(:, :, jst-1))
        end do
      end do ordering

      if (debug%info) then
        do ist = 1, st%nst
          ! there do not seem to be counted iterations here
          write(message(1), '(a,i4,a,i4,a,es13.6)') &
            'Debug: PLAN Eigensolver - ik', ik, ' ist ', ist, ' res ', res(ist)
          call messages_info(1, namespace=namespace)
        end do
      end if

      ! If the maximum mat-vecs is exceeded, get out of here.
      if (matvec > maxmatvecs) exit outer_loop

      ! Restart if any eigenpair is converged.
      if (nconv > 0) then
        nec = nec + nconv
        if (d2+1 > krylov) d1 = d2
        cycle outer_loop
      end if

      ! Preconditioning
      call lalg_copy(mesh%np, dim, av(:, :, d1 + 1), aux)
      call X(preconditioner_apply)(pre, namespace, mesh, hm, aux(:,:), vv(:,:, d1+1), ik)

    end do inner_loop
  end do outer_loop

  do ist = 1, st%nst
    call states_elec_set_state(st, mesh, ist, ik, eigenvec(:, :, ist))
    st%eigenval(ist, ik) = eigenval(ist)
    diff(ist) = res(ist)
  end do

  converged = nec
  niter = niter + matvec

  SAFE_DEALLOCATE_A(eigenval)
  SAFE_DEALLOCATE_A(eigenvec)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(av)
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(ham)
  SAFE_DEALLOCATE_A(hevec)
  SAFE_DEALLOCATE_A(aux)
  POP_SUB(X(eigensolver_plan))

contains

  ! ---------------------------------------------------------
  subroutine residual(hv, vv, ee, res, rr)
    R_TYPE,  intent(inout) :: hv(:,:)
    R_TYPE,  intent(inout) :: vv(:,:)
    FLOAT,   intent(in)    :: ee
    R_TYPE,  intent(inout) :: res(:,:)
    FLOAT,   intent(out)   :: rr

    PUSH_SUB(X(eigensolver_plan).residual)

    res = hv - ee*vv
    rr = X(mf_nrm2)(mesh, dim, res)

    POP_SUB(X(eigensolver_plan).residual)
  end subroutine residual

end subroutine X(eigensolver_plan)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
