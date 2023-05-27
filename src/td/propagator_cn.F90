!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2022 N. Tancogne-Dejean
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

#include "global.h"

module propagator_cn_oct_m
  use accel_oct_m
  use batch_ops_oct_m
  use debug_oct_m
  use density_oct_m
  use exponential_oct_m
  use ext_partner_list_oct_m
  use gauge_field_oct_m
  use grid_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use interaction_partner_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_base_oct_m
  use propagation_ops_elec_oct_m
  use solvers_oct_m
  use space_oct_m
  use sparskit_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use wfs_elec_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                    &
    td_crank_nicolson

  type(namespace_t),        pointer, private :: namespace_p
  type(mesh_t),             pointer, private :: mesh_p
  type(hamiltonian_elec_t), pointer, private :: hm_p
  type(propagator_base_t),  pointer, private :: tr_p
  type(states_elec_t),      pointer, private :: st_p
  integer,                           private :: ik_op, ist_op, dim_op
  FLOAT,                             private :: t_op, dt_op

contains

  ! ---------------------------------------------------------
  !> Crank-Nicolson propagator
  subroutine td_crank_nicolson(hm, namespace, space, gr, st, tr, time, dt, ionic_scale, ions_dyn, ions, ext_partners, use_sparskit)
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(namespace_t),        target, intent(in)    :: namespace
    type(space_t),                    intent(in)    :: space
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_base_t),  target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions_dyn
    type(ions_t),                     intent(inout) :: ions
    type(partner_list_t),             intent(in)    :: ext_partners
    logical,                          intent(in)    :: use_sparskit

    CMPLX, allocatable :: zpsi(:), rhs(:)
    integer :: ik, ist, idim, np, max_iter
    FLOAT :: cgtol = CNST(1.0e-12)
    integer :: ib, minst, maxst
    integer, allocatable :: iter_used(:)
    FLOAT, allocatable :: residue(:)
    type(density_calc_t) :: dens_calc
    type(wfs_elec_t) :: psib_rhs
    type(gauge_field_t), pointer :: gfield

    PUSH_SUB(propagator_dt.td_crank_nicolson)

    !TODO: Add gauge field support
    gfield => list_get_gauge_field(ext_partners)
    if(associated(gfield)) then
      ASSERT(gauge_field_is_propagated(gfield) .eqv. .false.)
    end if


#ifndef HAVE_SPARSKIT
    if (use_sparskit) then
      message(1) = "Cannot use SPARSKIT in Crank-Nicolson propagator: not compiled with SPARSKIT support."
      call messages_fatal(1, namespace=namespace)
    end if
#endif

    np = gr%mesh%np

    ! define pointer and variables for usage in td_zop, td_zopt routines
    namespace_p => namespace
    mesh_p      => gr%mesh
    hm_p        => hm
    tr_p        => tr
    st_p        => st
    dt_op = dt
    t_op  = time - dt/M_TWO
    dim_op = st%d%dim

    ! we (ab)use exponential_apply to compute (1-i\delta t/2 H_n)\psi^n
    ! exponential order needs to be only 1
    tr%te%exp_method = EXP_TAYLOR
    tr%te%exp_order  = 1

    SAFE_ALLOCATE(zpsi(1:np*st%d%dim))
    SAFE_ALLOCATE(rhs(1:np*st%d%dim))

    !move the ions to time 'time - dt/2', and save the current status to return to it later.
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, space, ions_dyn, ions, &
      ext_partners, time - M_HALF*dt, ionic_scale*M_HALF*dt, save_pos = .true.)

    if (family_is_mgga_with_exc(hm%xc)) then
      call potential_interpolation_interpolate(tr%vksold, 3, &
        time, dt, time -dt/M_TWO, hm%vhxc, vtau = hm%vtau)
    else
      call potential_interpolation_interpolate(tr%vksold, 3, &
        time, dt, time -dt/M_TWO, hm%vhxc)
    end if

    call propagation_ops_elec_update_hamiltonian(namespace, space, st, gr%mesh, hm, ext_partners, time - dt*M_HALF)

    call density_calc_init(dens_calc, st, gr, st%rho)

    ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
    do ik = st%d%kpt%start, st%d%kpt%end
      call propagation_ops_do_pack(st, hm, st%group%block_start, ik)
      do ib = st%group%block_start, st%group%block_end
        if (ib + 1 <= st%group%block_end) call propagation_ops_do_pack(st, hm, ib+1, ik)
        call accel_set_stream(ib)

        call st%group%psib(ib, ik)%copy_to(psib_rhs)
        dt_op = -dt
        call propagator_qmr_op_batch (st%group%psib(ib, ik), psib_rhs)
        dt_op = dt

        if (hamiltonian_elec_inh_term(hm)) then
          call batch_axpy(np, dt, hm%inh_st%group%psib(ib, ik), psib_rhs)
        end if

        if(use_sparskit) then !Unbatchified path
          minst = states_elec_block_min(st, ib)
          maxst = states_elec_block_max(st, ib)
          do ist = minst, maxst
            ! put the values in a continuous array
            do idim = 1, st%d%dim
              call batch_get_state(st%group%psib(ib, ik), (/ist, idim/), gr%mesh%np, zpsi((idim - 1)*np+1:idim*np))
              call batch_get_state(psib_rhs, (/ist, idim/), np, rhs((idim - 1)*np+1:idim*np))
            end do

            ist_op = ist
            ik_op = ik

            call zsparskit_solver_run(namespace, tr%tdsk, td_zop, td_zopt, zpsi, rhs)

            do idim = 1, st%d%dim
              call batch_set_state(st%group%psib(ib, ik), (/ist, idim/), gr%mesh%np, zpsi((idim - 1)*np+1:idim*np))
            end do

          end do

        else

          SAFE_ALLOCATE(iter_used(psib_rhs%nst))
          SAFE_ALLOCATE(residue(psib_rhs%nst))
          max_iter = 2000
          call zbatch_qmr_dotu(namespace, gr%mesh, st, st%group%psib(ib, ik), psib_rhs, &
            propagator_qmr_op_batch, max_iter, iter_used, residue, cgtol, .false.)

          if (any(iter_used == max_iter)) then
            write(message(1),'(a)') 'The linear solver used for the Crank-Nicolson'
            write(message(2),'(a)') 'propagator did not converge: '
            do ist=1, psib_rhs%nst
              write(message(ist+2),'(a,i2,a,es14.4)') 'State = ', ist, ' Residual = ', residue(ist)
            end do
            call messages_warning(psib_rhs%nst+2, namespace=namespace)
          end if

          SAFE_DEALLOCATE_A(iter_used)
          SAFE_DEALLOCATE_A(residue)
        end if

        !use the dt propagation to calculate the density
        call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))

        call psib_rhs%end()

        call propagation_ops_do_unpack(st, hm, ib, ik)
        if (ib-1 >= st%group%block_start) call propagation_ops_finish_unpack(st, hm, ib-1, ik)
      end do
      call propagation_ops_finish_unpack(st, hm, st%group%block_end, ik)
    end do

    call density_calc_end(dens_calc)

    !restore to time 'time - dt'
    call propagation_ops_elec_restore_ions(tr%propagation_ops_elec, ions_dyn, ions)

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(rhs)

    POP_SUB(propagator_dt.td_crank_nicolson)
  end subroutine td_crank_nicolson
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> operators for Crank-Nicolson scheme
  subroutine td_zop(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(td_zop)

    SAFE_ALLOCATE(zpsi(1:mesh_p%np_part, 1:dim_op))
    zpsi = M_z0
    do idim = 1, dim_op
      zpsi(1:mesh_p%np, idim) = xre((idim-1)*mesh_p%np+1:idim*mesh_p%np) + &
        M_zI * xim((idim-1)*mesh_p%np+1:idim*mesh_p%np)
    end do

    call exponential_apply(tr_p%te, namespace_p, mesh_p, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO)

    do idim = 1, dim_op
      yre((idim-1)*mesh_p%np+1:idim*mesh_p%np) = TOFLOAT(zpsi(1:mesh_p%np, idim))
      yim((idim-1)*mesh_p%np+1:idim*mesh_p%np) = aimag(zpsi(1:mesh_p%np, idim))
    end do

    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(td_zop)
  end subroutine td_zop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Transpose of H (called e.g. by bi-conjugate gradient solver)
  subroutine td_zopt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(td_zopt)

    ! To act with the transpose of H on the wfn we apply H to the conjugate of psi
    ! and conjugate the resulting hpsi (note that H is not a purely real operator
    ! for scattering wavefunctions anymore).
    SAFE_ALLOCATE(zpsi(1:mesh_p%np_part, 1:dim_op))
    zpsi = M_z0
    do idim = 1, dim_op
      zpsi(1:mesh_p%np, idim) = xre((idim-1)*mesh_p%np+1:idim*mesh_p%np) - &
        M_zI * xim((idim-1)*mesh_p%np+1:idim*mesh_p%np)
    end do

    call exponential_apply(tr_p%te, namespace_p, mesh_p, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO)

    do idim = 1, dim_op
      yre((idim-1)*mesh_p%np+1:idim*mesh_p%np) =   TOFLOAT(zpsi(1:mesh_p%np, idim))
      yim((idim-1)*mesh_p%np+1:idim*mesh_p%np) = - aimag(zpsi(1:mesh_p%np, idim))
    end do

    SAFE_DEALLOCATE_A(zpsi)
    POP_SUB(td_zopt)
  end subroutine td_zopt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> operators for Crank-Nicolson scheme
  subroutine propagator_qmr_op_batch(xxb, yyb)
    type(wfs_elec_t), intent(inout)  :: xxb
    type(wfs_elec_t), intent(inout) :: yyb

    PUSH_SUB(propagator_qmr_op_batch)

    call xxb%copy_data_to(mesh_p%np, yyb)
    call exponential_apply_batch(tr_p%te, namespace_p, mesh_p, hm_p, yyb, -dt_op/M_TWO)

    POP_SUB(propagator_qmr_op_batch)
  end subroutine propagator_qmr_op_batch
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> for complex symmetric matrices
  !! W Chen and B Poirier, J Comput Phys 219, 198-209 (2006)
  !!
  !! Adapted from linear_solver_qmr_dotp
  !! This modified version does not have the preconditioner, as it is not used for the Crank-Nicolson,
  !! which allows for improved performances
  subroutine zbatch_qmr_dotu(namespace, mesh, st, xb, bb, op, max_iter, iter_used, &
    residue, threshold, use_initial_guess)
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(states_elec_t),      intent(in)    :: st
    type(wfs_elec_t),         intent(inout) :: xb
    type(wfs_elec_t),         intent(in)    :: bb
    interface
      subroutine op(xxb, yyb)         !< the preconditioner
        import wfs_elec_t
        implicit none
        type(wfs_elec_t), intent(inout) :: xxb
        type(wfs_elec_t), intent(inout) :: yyb
      end subroutine op
    end interface
    integer,                  intent(in)    :: max_iter
    integer,                  intent(out)   :: iter_used(:)
    FLOAT,                    intent(out)   :: residue(:)   !< the residue = abs(Ax-b)
    FLOAT,                    intent(in)    :: threshold    !< convergence threshold
    logical, optional,        intent(in)    :: use_initial_guess


    type(wfs_elec_t) :: vvb, res, qqb, ppb, deltax, deltar
    integer             :: ii, iter
    FLOAT, allocatable  :: rho(:), oldrho(:), norm_b(:), xsi(:), gamma(:), alpha(:), theta(:), oldtheta(:), saved_res(:)
    FLOAT, allocatable  :: oldgamma(:)
    CMPLX, allocatable :: eta(:), beta(:), delta(:), eps(:), exception_saved(:, :, :)
    integer, allocatable :: status(:), saved_iter(:)

    integer, parameter ::        &
      QMR_NOT_CONVERGED    = 0,  &
      QMR_CONVERGED        = 1,  &
      QMR_RES_ZERO         = 2,  &
      QMR_B_ZERO           = 3,  &
      QMR_BREAKDOWN_PB     = 4,  &
      QMR_BREAKDOWN_VZ     = 5,  &
      QMR_BREAKDOWN_QP     = 6,  &
      QMR_BREAKDOWN_GAMMA  = 7

    PUSH_SUB(zbatch_qmr_dotu)

    SAFE_ALLOCATE(rho(1:xb%nst))
    SAFE_ALLOCATE(oldrho(1:xb%nst))
    SAFE_ALLOCATE(norm_b(1:xb%nst))
    SAFE_ALLOCATE(xsi(1:xb%nst))
    SAFE_ALLOCATE(gamma(1:xb%nst))
    SAFE_ALLOCATE(oldgamma(1:xb%nst))
    SAFE_ALLOCATE(alpha(1:xb%nst))
    SAFE_ALLOCATE(eta(1:xb%nst))
    SAFE_ALLOCATE(theta(1:xb%nst))
    SAFE_ALLOCATE(oldtheta(1:xb%nst))
    SAFE_ALLOCATE(beta(1:xb%nst))
    SAFE_ALLOCATE(delta(1:xb%nst))
    SAFE_ALLOCATE(eps(1:xb%nst))
    SAFE_ALLOCATE(saved_res(1:xb%nst))

    SAFE_ALLOCATE(status(1:xb%nst))
    SAFE_ALLOCATE(saved_iter(1:xb%nst))

    SAFE_ALLOCATE(exception_saved(1:mesh%np, 1:st%d%dim, 1:xb%nst))

    call xb%copy_to(vvb)
    call xb%copy_to(res)
    call xb%copy_to(qqb)
    call xb%copy_to(ppb)
    call xb%copy_to(deltax)
    call xb%copy_to(deltar)

    if (optional_default(use_initial_guess, .true.)) then
      ! Compared to the original algorithm, we assume that we have an initial guess
      ! This means that instead of setting x^(0)=0, we have x^(0)=xb
      ! TODO: We should implement here the proper recursion for x^(0) /= 0
      ! as published in "Preconditioning of Symmetric, but Highly Indefinite Linear Systems"
      ! R. W. Freund, 15th IMACS World Congress on Scientific Computation, Modelling and Applied Mathematics
      ! 2, 551 (1997)
      call op(xb, vvb)

      call batch_xpay(mesh%np, bb, CNST(-1.0), vvb)
      call vvb%copy_data_to(mesh%np, res)

      ! Norm of the right-hand side
      call mesh_batch_nrm2(mesh, vvb, rho)
      call mesh_batch_nrm2(mesh, bb, norm_b)

      do ii = 1, xb%nst
        ! If the initial guess is a good enough solution
        if (abs(rho(ii)) <= threshold) then
          status(ii) = QMR_RES_ZERO
          residue(ii) = rho(ii)
          call batch_get_state(xb, ii, mesh%np, exception_saved(:, :, ii))
          saved_iter(ii) = 0
          saved_res(ii) = residue(ii)
        end if
      end do

    else ! If we don't know any guess, let's stick to the original algorithm

      call batch_set_zero(xb)
      call bb%copy_data_to(mesh%np, vvb)
      call vvb%copy_data_to(mesh%np, res)
      call mesh_batch_nrm2(mesh, bb, norm_b)
      rho = norm_b

    end if

    status = QMR_NOT_CONVERGED

    iter = 0

    do ii = 1, xb%nst

      residue(ii) = rho(ii)
      ! if b is zero, the solution is trivial
      if (status(ii) == QMR_NOT_CONVERGED .and. abs(norm_b(ii)) <= M_TINY) then
        exception_saved = M_ZERO
        status(ii) = QMR_B_ZERO
        residue(ii) = norm_b(ii)
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if

    end do

    xsi = rho
    gamma = M_ONE
    oldgamma = gamma
    eta   = CNST(-1.0)
    alpha = M_ONE
    theta = M_ZERO

    do while(iter < max_iter)
      iter = iter + 1

      ! Exit condition
      if (all(status /= QMR_NOT_CONVERGED)) exit

      ! Failure of the algorithm
      do ii = 1, xb%nst
        if (status(ii) == QMR_NOT_CONVERGED .and. (abs(rho(ii)) < M_TINY .or. abs(xsi(ii)) < M_TINY)) then
          call batch_get_state(xb, ii, mesh%np, exception_saved(:, :, ii))
          status(ii) = QMR_BREAKDOWN_PB
          saved_iter(ii) = iter
          saved_res(ii) = residue(ii)
        end if

        alpha(ii) = alpha(ii)*xsi(ii)/rho(ii)
      end do

      ! v^(i) = v^(i) / \rho_i
      call batch_scal(mesh%np, M_ONE/rho, vvb, a_full = .false.)
      ! \delta_i = v^(i)\ldotu z^(i)
      call zmesh_batch_dotp_vector(mesh, vvb, vvb, delta, cproduct=.true.)

      delta = delta / alpha


      !If \delta_i = 0, method fails
      do ii = 1, xb%nst
        if (status(ii) == QMR_NOT_CONVERGED .and. abs(delta(ii)) < M_TINY) then
          call batch_get_state(xb, ii, mesh%np, exception_saved(:, :, ii))
          status(ii) = QMR_BREAKDOWN_VZ
          saved_iter(ii) = iter
          saved_res(ii) = residue(ii)
        end if
      end do

      if (iter == 1) then
        ! q^(1) = z^(1)
        call vvb%copy_data_to(mesh%np, qqb)
      else
        ! q^(i) = z^(i) - (\rho_i\delta_i)/(\eps_{i-1})q^(i-1)
        call batch_xpay(mesh%np, vvb, -rho*delta/eps, qqb, a_full = .false.)
      end if

      ! ppb = H*qqb
      call op(qqb, ppb)
      ! p^(i) = \alpha_{i+1} (H-shift)*q^(i)
      call batch_scal(mesh%np, alpha, ppb, a_full = .false.)

      ! \eps_i = q^{(i)}\ldotu p^{(i)}
      call zmesh_batch_dotp_vector(mesh, qqb, ppb, eps, cproduct=.true.)

      ! If \eps_i == 0, method fails
      do ii = 1, xb%nst
        if (status(ii) == QMR_NOT_CONVERGED .and. abs(eps(ii)) < M_TINY) then
          call batch_get_state(xb, ii, mesh%np, exception_saved(:, :, ii))
          status(ii) = QMR_BREAKDOWN_QP
          saved_iter(ii) = iter
          saved_res(ii) = residue(ii)
        end if

        beta(ii) = eps(ii)/delta(ii)
      end do

      ! v^(i+1) = p^(i) - \beta_i v^(i)
      call batch_xpay(mesh%np, ppb, -beta, vvb, a_full = .false.)

      do ii = 1, xb%nst
        oldrho(ii) = rho(ii)
      end do

      ! \rho_{i+1} = ||v^{i+1}||_2
      call mesh_batch_nrm2(mesh, vvb, rho)

      ! \xsi_{i+1} = ||z^{i+1}||_2
      xsi = rho / (alpha**2)

      do ii = 1, xb%nst

        oldtheta(ii) = theta(ii)
        ! \theta_i = \rho_{i+1}/(\gamma_{i-1} |\beta_i|)
        theta(ii) = rho(ii)/(gamma(ii)*abs(beta(ii)))

        oldgamma(ii) = gamma(ii)
        ! \gamma_i = 1/sqrt(1+\theta_i^2)
        gamma(ii) = M_ONE/sqrt(M_ONE + theta(ii)**2)

        ! If \gamma_i == 0, method fails
        if (status(ii) == QMR_NOT_CONVERGED .and. abs(gamma(ii)) < M_TINY) then
          call batch_get_state(xb, ii, mesh%np, exception_saved(:, :, ii))
          status(ii) = QMR_BREAKDOWN_GAMMA
          saved_iter(ii) = iter
          saved_res(ii) = residue(ii)
        end if

        ! \eta_i = -\eta_{i-1}\rho_i \gamma_i^2/ (\beta_i \gamma_{i-1}^2)
        eta(ii) = -eta(ii)*oldrho(ii)*gamma(ii)**2/(beta(ii)*oldgamma(ii)**2)
      end do

      if (iter == 1) then

        ! \delta_x^(1) = \eta_1 \alpha_2 q^{(1)}
        call qqb%copy_data_to(mesh%np, deltax)
        call batch_scal(mesh%np, eta*alpha, deltax, a_full = .false.)

        ! \delta_r^(1) = \eta_1 p^1
        call ppb%copy_data_to(mesh%np, deltar)
        call batch_scal(mesh%np, eta, deltar, a_full = .false.)

      else

        ! \delta_x^{i} = (\theta_{i-1}\gamma_i)^2 \delta_x^{i-1} + \eta_i\alpha_{i+1} q^i
        call batch_scal(mesh%np, (oldtheta*gamma)**2, deltax, a_full = .false.)
        call batch_axpy(mesh%np, eta*alpha, qqb, deltax, a_full = .false.)

        ! \delta_r^{i} = (\theta_{i-1}\gamma_i)^2 \delta_r^{i-1} + \eta_i p^i
        call batch_scal(mesh%np, (oldtheta*gamma)**2, deltar, a_full = .false.)
        call batch_axpy(mesh%np, eta, ppb, deltar, a_full = .false.)

      end if

      ! FIXME: if the states are converged, why changing them here
      ! x^{i} = x^{i-1} + \delta_x^{i}
      call batch_axpy(mesh%np, M_ONE, deltax, xb)
      ! r^{i} = r^{i-1} - \delta_r^i
      ! This is given by r^{i} = b - Hx^{i}
      call batch_axpy(mesh%np, CNST(-1.0), deltar, res)

      ! We compute the norm of the residual
      call mesh_batch_nrm2(mesh, res, residue)
      do ii = 1, xb%nst
        residue(ii) = residue(ii)/norm_b(ii)
      end do

      ! Convergence is reached once the residues are smaller than the threshold
      do ii = 1, xb%nst
        if (status(ii) == QMR_NOT_CONVERGED .and. residue(ii) < threshold) then
          status(ii) = QMR_CONVERGED
          if (debug%info) then
            write(message(1),*) 'Debug: State ', xb%ist(ii), ' converged in ', iter, ' iterations.'
            call messages_info(1, namespace=namespace)
          end if
        end if
      end do

    end do

    do ii = 1, xb%nst
      if (status(ii) == QMR_NOT_CONVERGED .or. status(ii) == QMR_CONVERGED) then
        ! We stop at the entrance of the next iteraction, so we substract one here
        iter_used(ii) = iter -1
      else
        call batch_set_state(xb, ii, mesh%np, exception_saved(:, :, ii))
        iter_used(ii) = saved_iter(ii)
        residue(ii) = saved_res(ii)
      end if

      select case (status(ii))
      case (QMR_NOT_CONVERGED)
        write(message(1), '(a)') "QMR solver not converged!"
        write(message(2), '(a)') "Try increasing the maximum number of iterations or the tolerance."
        call messages_warning(2, namespace=namespace)
      case (QMR_BREAKDOWN_PB)
        write(message(1), '(a)') "QMR breakdown, cannot continue: b or P*b is the zero vector!"
        call messages_warning(1, namespace=namespace)
      case (QMR_BREAKDOWN_VZ)
        write(message(1), '(a)') "QMR breakdown, cannot continue: v^T*z is zero!"
        call messages_warning(1, namespace=namespace)
      case (QMR_BREAKDOWN_QP)
        write(message(1), '(a)') "QMR breakdown, cannot continue: q^T*p is zero!"
        call messages_warning(1, namespace=namespace)
      case (QMR_BREAKDOWN_GAMMA)
        write(message(1), '(a)') "QMR breakdown, cannot continue: gamma is zero!"
        call messages_warning(1, namespace=namespace)
      end select

    end do

    call vvb%end()
    call res%end()
    call qqb%end()
    call ppb%end()
    call deltax%end()
    call deltar%end()

    SAFE_DEALLOCATE_A(exception_saved)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(oldrho)
    SAFE_DEALLOCATE_A(norm_b)
    SAFE_DEALLOCATE_A(xsi)
    SAFE_DEALLOCATE_A(gamma)
    SAFE_DEALLOCATE_A(oldgamma)
    SAFE_DEALLOCATE_A(alpha)
    SAFE_DEALLOCATE_A(eta)
    SAFE_DEALLOCATE_A(theta)
    SAFE_DEALLOCATE_A(oldtheta)
    SAFE_DEALLOCATE_A(beta)
    SAFE_DEALLOCATE_A(delta)
    SAFE_DEALLOCATE_A(eps)

    SAFE_DEALLOCATE_A(saved_res)
    SAFE_DEALLOCATE_A(status)
    SAFE_DEALLOCATE_A(saved_iter)

    POP_SUB(zbatch_qmr_dotu)
  end subroutine zbatch_qmr_dotu

end module propagator_cn_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
