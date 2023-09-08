!! Copyright (C) 2023. A Buccheri.
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

!> @brief Driver for Chebyshev filter-based solver.
!!
!! This routine is implemented according to Algorithm 4.1 in the paper:
!! ""Chebyshev-filtered subspace iteration method free of sparse diagonalization for
!!  solving the Kohn–Sham equation"". http://dx.doi.org/10.1016/j.jcp.2014.06.056
!!
!! The scaled Chebyshev algorithm is always utilised, as we get an
!! estimate of the lowest eigenvalue of H from the lowest Ritz value of
!! the prior step. The reason for the scaling is to prevent a potential overflow,
!! which may happen if the filter degree is large or if the smallest eigenvalue is
!! mapped far away from −1.
!!
!! For the first SCF step, a Chebyshev polynomial of filter_params%degree is applied
!! filter_params%n_iter times to a set of search vectors found from the initial guess
!! for the density.
subroutine X(chebyshev_filter_solver)(namespace, sdiag, mesh, st, hm, ik, subspace_tol, &
  filter_params, scf_iter)
  type(namespace_t),             intent(in)    :: namespace        !< Calling namespace
  type(subspace_t),              intent(in)    :: sdiag            !< Subspace diagonalisation choice
  class(mesh_t),                 intent(in)    :: mesh             !< Real-space mesh
  type(states_elec_t),           intent(inout) :: st               !< Eigenstates
  type(hamiltonian_elec_t),      intent(in)    :: hm               !< Hamiltonian
  integer,                       intent(in)    :: ik               !< k-point index
  FLOAT,                         intent(in)    :: subspace_tol     !< Subspace iterative solver tolerance
  type(eigen_chebyshev_t),       intent(in)    :: filter_params    !< Chebyshev filter parameters
  integer,                       intent(in)    :: scf_iter         !< SCF iteration

  type(chebyshev_filter_bounds_t), pointer :: bounds               !< Filter bounds
  FLOAT :: a_l, lower_bound, upper_bound                           !< Limits

  PUSH_SUB(X(chebyshev_filter_solver))

  if (scf_iter == 1 .and. filter_params%n_iter > 0) then
    ! Estimate for first SCF step, applying the filter iteratively
    call X(firstscf_iterative_chebyshev_filter)(namespace, sdiag, mesh, st, hm, ik, subspace_tol, filter_params)
  else
    ! Update lower bound and min eigenvalue from largest and smallest prior Ritz values
    lower_bound = maxval(st%eigenval(:, ik)) + CNST(1e-3)
    a_l = minval(st%eigenval(:, ik))
    ! Estimate upper bound with Lanczos estimator
    upper_bound = X(upper_bound_estimator)(namespace, mesh, st, hm, ik, filter_params%n_lanczos)
    bounds => chebyshev_filter_bounds_t(lower_bound, upper_bound, a_l=a_l)

    call X(chebyshev_filter)(namespace, mesh, st, hm, filter_params%degree, bounds, ik)

    call X(states_elec_orthogonalization_full)(st, namespace, mesh, ik)
    SAFE_DEALLOCATE_P(bounds)
  endif

  POP_SUB(X(chebyshev_filter_solver))
end subroutine X(chebyshev_filter_solver)


!> @brief Iterative application of Chebyshev filter, for use with the first SCF step.
!!
!! The initial search vectors are found using the initial for the density (sum of atomic
!! densities LCAO, etc), and returned in st: st%group%psib.
subroutine X(firstscf_iterative_chebyshev_filter)(namespace, sdiag, mesh, st, hm, ik, tolerance, filter_params)
  type(namespace_t),        intent(in)    :: namespace      !< Calling namespace
  type(subspace_t),         intent(in)    :: sdiag          !< Subspace diagonalisation choice
  class(mesh_t),            intent(in)    :: mesh           !< Real-space mesh
  type(states_elec_t),      intent(inout) :: st             !< Initial guess at search vectors, and much more
  type(hamiltonian_elec_t), intent(in)    :: hm             !< Hamiltonian
  integer,                  intent(in)    :: ik             !< k-point index
  FLOAT,                    intent(in)    :: tolerance      !< Subspace iterative solver tolerance
  type(eigen_chebyshev_t),  intent(in)    :: filter_params  !< Chebyshev filter parameters

  FLOAT,   allocatable :: e_diff(:)                         !< Differences in the eigenvalues
  logical, allocatable :: converged(:)                      !< Has each eigenvalue converged
  R_TYPE,  allocatable :: search_vector(:, :)               !< Random search vector
  type(chebyshev_filter_bounds_t), pointer :: bounds        !< Filter bounds
  FLOAT   :: upper_bound                                    !< Upper bound
  integer :: ifilter, ist                                   !< Loop indices

  PUSH_SUB(X(firstscf_iterative_chebyshev_filter))

  SAFE_ALLOCATE(search_vector(1:mesh%np_part, 1:hm%d%dim))
  call states_elec_generate_random_vector(mesh, st, search_vector, normalized=.true.)

  ! Initial bounds estimate
  call X(filter_bounds_estimator)(namespace, mesh, hm, ik, filter_params%n_lanczos, &
    filter_params%bound_mixing, search_vector, bounds)
  SAFE_DEALLOCATE_A(search_vector)

  ! Diagonalise subspace to get Ritz pairs
  SAFE_ALLOCATE(e_diff(1:st%nst))
  call X(subspace_diag)(sdiag, namespace, mesh, st, hm, ik, st%eigenval(:, ik), e_diff)

  converged = e_diff < tolerance
  if ((all(converged))) return

  upper_bound = bounds%upper
  SAFE_DEALLOCATE_P(bounds)
  bounds => chebyshev_filter_bounds_t(maxval(st%eigenval(:, ik)), &    ! Largest Ritz value
    upper_bound, &                       ! Unchanged
    a_l=minval(st%eigenval(:, ik)))  ! Smallest Ritz value

  do ifilter = 1, filter_params%n_iter
    ! Apply Chebyshev filter
    call X(chebyshev_filter)(namespace, mesh, st, hm, filter_params%degree, bounds, ik)

    ! Orthogonalise returned vectors
    call X(states_elec_orthogonalization_full)(st, namespace, mesh, ik)

    ! Diagonalise subspace to get Ritz pairs
    call X(subspace_diag)(sdiag, namespace, mesh, st, hm, ik, st%eigenval(:, ik), e_diff)

    converged = e_diff < tolerance
    if ((all(converged))) return

    SAFE_DEALLOCATE_P(bounds)
    bounds => chebyshev_filter_bounds_t(maxval(st%eigenval(:, ik)), &    ! Largest Ritz value
      upper_bound, &                       ! Unchanged
      a_l=minval(st%eigenval(:, ik)))  ! Smallest Ritz value
  enddo

  if(debug%info) then
    write(message(1), '(A, X, I7)') 'Debug: Chebyshev 1st iterative step state convergence for ik = ', ik
    call messages_info(1, namespace=namespace)

    do ist = 1, size(converged)
      write(message(1), '(I8, X, f16.12, X, L)') ist, e_diff(ist), converged(ist)
      call messages_info(1, namespace=namespace)
    enddo
  end if

  SAFE_DEALLOCATE_A(e_diff)
  SAFE_DEALLOCATE_P(bounds)
  POP_SUB(X(firstscf_iterative_chebyshev_filter))
end subroutine X(firstscf_iterative_chebyshev_filter)


!> @brief Chebyshev Filter.
!!
!! Filter an eigenspectrum by an m-degree Chebyshev polynomial that dampens on the interval [a,b].
!! Based on algorithm 3.2 of ""Chebyshev-filtered subspace iteration method free of sparse diagonalization
!! for solving the Kohn–Sham equation"". http://dx.doi.org/10.1016/j.jcp.2014.06.056
!!
!! Application of the simple or scaled filter depends upon the choice of bounds.
!! In the simple case, sigma will reduce to 1.
subroutine X(chebyshev_filter) (namespace, mesh, st, hm, degree, bounds, ik)
  type(namespace_t),               intent(in)    :: namespace  !< Calling namespace
  type(mesh_t),                    intent(in)    :: mesh       !< Real-space mesh
  class(hamiltonian_elec_t),       intent(in)    :: hm         !< Hamiltonian
  integer,                         intent(in)    :: degree     !< Chebyshev polynomial degree
  !                                                            !< (filter applied using recursivel definition)
  type(chebyshev_filter_bounds_t), intent(in)    :: bounds     !< Polynomial filter bounds of the subspace to dampen
  integer,                         intent(in)    :: ik         !< k-point index
  type(states_elec_t), target,  intent(inout)    :: st         !< KS containing input eigenvectors and
  !                                                            !< returned with updated eigenvectors

  FLOAT :: c                                                   !< Centre of the filter limits
  FLOAT :: hw                                                  !< Half-width of the filter limits
  FLOAT :: sigma, sigma_new, tau                               !< Filter scaling
  type(wfs_elec_t), pointer :: Y_old, Y, Y_new                 !< Pointers, used to avoid memory transfer
  integer :: idegree, iblock                                   !< Loop indices
  type(profile_t), save :: prof                                !< Profiling data
  type(batch_pointer_t) :: physical_memory(3)                  !< Physical memory
  integer :: iperm(3)                                          !< Indices for the circular permutation of the pointers

  PUSH_SUB(X(chebyshev_filter))

  call profiling_in(prof, TOSTRING(X(CHEBY)))

  hw = bounds%half_width()
  c = bounds%center()
  sigma = bounds%sigma()
  ! Tau is not updated w.r.t. application  of the filter
  tau = M_TWO / sigma

  SAFE_ALLOCATE(physical_memory(2)%batch)
  SAFE_ALLOCATE(physical_memory(3)%batch)
  SAFE_ALLOCATE(Y_old)
  SAFE_ALLOCATE(Y)
  SAFE_ALLOCATE(Y_new)

  ! Block index
  do iblock = st%group%block_start, st%group%block_end

    physical_memory(1)%batch => st%group%psib(iblock, ik)
    call hamiltonian_elec_base_set_phase_corr(hm%hm_base, mesh, st%group%psib(iblock, ik))
    call physical_memory(1)%batch%copy_to(physical_memory(2)%batch)
    call physical_memory(1)%batch%copy_to(physical_memory(3)%batch)

    ! Define Y = (H X - c X) * (sigma / hw)  in stages:
    ! Y = H X
    call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, physical_memory(1)%batch, physical_memory(2)%batch)
    ! Y -> (Y - cX)
    call batch_axpy(mesh%np, -c, physical_memory(1)%batch, physical_memory(2)%batch)
    ! Y -> Y * sigma / hw
    call batch_scal(mesh%np, sigma / hw, physical_memory(2)%batch)

    do idegree = 2, degree
      sigma_new = M_ONE / (tau - sigma)

      ! Having st%group%psib(iblock, ik), Yb, and HY been allocated, we have enough memory to
      ! avoid memory transfer. To do this, we employ pointers that we swap during the execution
      ! In order to avoid memory transfer, we do a circular permutation of the pointers
      iperm = [mod(idegree - 2, 3) + 1, mod(idegree - 1, 3) + 1, mod(idegree, 3) + 1]
      Y_old => physical_memory(iperm(1))%batch
      Y     => physical_memory(iperm(2))%batch
      Y_new => physical_memory(iperm(3))%batch

      ! Construct Y_new = (HY - c * Y) * (2 * sigma_new / hw) - (sigma * sigma_new) * X
      ! in stages:
      ! HY = H Y
      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, Y, Y_new)
      ! HY -> (HY - cY)
      call batch_axpy(mesh%np, -c, Y, Y_new)
      ! HY -> (2 * sigma_new / hw) * HY
      call batch_scal(mesh%np, M_TWO * sigma_new / hw, Y_new)
      ! HY -> HY - (sigma * sigma_new) * Xb
      call batch_axpy(mesh%np, -sigma * sigma_new, Y_old, Y_new)

      sigma = sigma_new
    end do

    ! Copy the final data
    call Y_new%copy_data_to(mesh%np, st%group%psib(iblock, ik))
    nullify(Y_old, Y, Y_new)
    nullify(physical_memory(1)%batch)
    call physical_memory(2)%batch%end()
    call physical_memory(3)%batch%end()

    ! Remove the phase correction
    call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, mesh, st%group%psib(iblock, ik))
  end do

  SAFE_DEALLOCATE_P(physical_memory(2)%batch)
  SAFE_DEALLOCATE_P(physical_memory(3)%batch)

  call profiling_out(prof)

  POP_SUB(X(chebyshev_filter))
end subroutine X(chebyshev_filter)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
