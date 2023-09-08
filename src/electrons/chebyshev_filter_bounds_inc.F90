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

!> @brief Compute a tridiagonal matrix with Lanczos.
!!
!! A k-step Lanczos leads to a Lanczos-decomposition:
!!
!! \f[ H V_k=V_k T_k+f_k e_k^{\mathrm{T}} \f]
!!
!! where \f$ V_k \f$ contains the k Lanczos basis, \(T_k\) is a size-k tridiagonal matrix,
!! \f$ f_k \f$ is a residual vector, and \(e_k\) is a length k unit vector with only the
!! first element nonzero.
!!
!! It's the caller's responsibility to allocate return arrays:
!! ```
!!  allocate(tridiagonal(size_t, size_t))
!!  allocate(f_residual(mesh%np, hm%d%dim))
!! ```
!! Note, `hamiltonian_elec_apply_single` is always initialised with state = 1.
!! If the Hamiltonian is different for each state (orbital-dependent functionals)
!! then the state index in principle needs to be provided - the bounds of the Hamiltonian should change.
!! This can be neglected if the maximum cutoff energy is used as an estimate of the upper bound.
subroutine X(lanczos_tridiagonal)(namespace, mesh, hm, ik, vs, tridiagonal, f_residual, size_t)
  type(namespace_t),         intent(in)    :: namespace          !< Calling namespace
  class(mesh_t),             intent(in)    :: mesh               !< Real-space mesh
  class(hamiltonian_elec_t), intent(in)    :: hm                 !< Hamiltonian
  integer,                   intent(in)    :: ik                 !< k-point index
  R_TYPE, contiguous,        intent(inout) :: vs(:, :)           !< Random search vector
  FLOAT,                     intent(out)   :: tridiagonal(:, :)  !< Tridiagonal matrix
  R_TYPE, contiguous,        intent(out)   :: f_residual(:, :)   !< Residual vector
  integer,                   intent(inout) :: size_t             !< Size of tridiagonal matrix. If beta is found to be ~ 0
  !                                                              !< for i < size_t, size_t will be set equal to i.

  integer, parameter :: dummy_i = 1                              !< Initialise a batch with a specific state
  R_TYPE, allocatable :: v0(:, :)                                !< Random search vector
  FLOAT :: alpha                                                 !< The dot product <f_residual, vs>
  FLOAT :: beta                                                  !< Norm of the residual vector
  integer :: n_iterations                                        !< Number of Lanczos iterations
  integer :: i                                                   !< Loop index
  type(profile_t), save :: prof

  PUSH_SUB(X(lanczos_tridiagonal))

  call profiling_in(prof, TOSTRING(X(BOUND_ESTIMATE)))

  ! Check random vector dimensions are correct
  ASSERT(size(vs, 1) == mesh%np_part)
  ASSERT(size(vs, 2) == hm%d%dim)

  ASSERT(hm%is_hermitian())
  n_iterations = size(tridiagonal, 1)
  ASSERT(size(tridiagonal, 2) == n_iterations)
  tridiagonal = M_ZERO

  ! v0 does not operate on H, hence it does not include boundaries in size allocation
  SAFE_ALLOCATE(v0(1:mesh%np, 1:hm%d%dim))

  ! f_residual = Hv
  call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, vs, f_residual, dummy_i, ik)
  ! alpha = <f_residual, vs>
  alpha = real((X(mf_dotp)(mesh, hm%d%dim, f_residual, vs)))
  tridiagonal(1, 1) = alpha
  ! f_residual -> f_residual - alpha * vs
  call lalg_axpy(mesh%np, hm%d%dim, -alpha, vs, f_residual)

  do i = 2, n_iterations
    ! beta = |f_residual|
    beta = X(mf_nrm2)(mesh, hm%d%dim, f_residual)
    if (beta < M_EPSILON) then
      size_t = i
      return
    endif

    ! v0 = copy(vs)
    call lalg_copy(mesh%np, hm%d%dim, vs, v0)

    ! vs = f_residual/beta
    call lalg_copy(mesh%np, hm%d%dim, f_residual, vs)
    call lalg_scal(mesh%np, hm%d%dim, M_ONE / beta, vs)

    ! f_residual = np.matmul(H, f_residual) / beta
    call X(hamiltonian_elec_apply_single)(hm, namespace, mesh, vs, f_residual, dummy_i, ik)

    ! f_residual -> f_residual  - beta * v0
    call lalg_axpy(mesh%np, hm%d%dim, -beta, v0, f_residual)

    ! alpha = <f_residual, vs>
    alpha = real(X(mf_dotp)(mesh, hm%d%dim, f_residual, vs))

    ! f_residual -> f_residual - alpha * vs
    call lalg_axpy(mesh%np, hm%d%dim, -alpha, vs, f_residual)

    tridiagonal(i, i) = alpha
    tridiagonal(i - 1, i) = beta
    tridiagonal(i, i - 1) = beta
  end do

  SAFE_DEALLOCATE_A(v0)

  call profiling_out(prof)

  POP_SUB(X(lanczos_tridiagonal))
end subroutine X(lanczos_tridiagonal)


!> @brief Estimate the upper bound of sigma(H_k) by k-step Lanczos.
!!
!! The Lanczos iteration quickly approximates the outermost eigenvalues, such that:
!!
!! \f[
!!    \left\|H V_k\right\|_2=\left\|V_k T_k+f_k e_k^{\mathrm{T}}\right\|_2
!!         \leqslant\left\|T_k\right\|_2+\left\|f_k\right\|_2
!! \f]
!!
!! We use \f$ left\|T_k\right\|_2+\left\|f_k\right\|_2 \f$ as an approximation
!! for the upper bound of H i.e. the maximal eigenvalue.
!!
!! Based on algorithm 4.4. in "Self-consistent-field calculations using "
!! "Chebyshev-filtered subspace iteration".
function X(upper_bound_estimator)(namespace, mesh, st, hm, ik, n_iterations) result(upper_bound)
  type(namespace_t),         intent(in) :: namespace      !< Calling namespace
  class(mesh_t),             intent(in) :: mesh           !< Real-space mesh
  type(states_elec_t),       intent(in) :: st             !< Basis information
  class(hamiltonian_elec_t), intent(in) :: hm             !< Hamiltonian
  integer,                   intent(in) :: ik             !< k-point index
  integer,                   intent(in) :: n_iterations   !< Number of Lanczos iterations, and consequentally the dimensions
  !                                                       !< of the tesulting tridiagonal matrix
  FLOAT :: upper_bound                                    !< Approximate max eigenvalue of H

  integer, parameter   :: MAX_LANCZOS_STEPS = 10          !< Maximum number of Lanczos steps allowed when determining upper bound
  R_TYPE,  allocatable :: vs(:, :)                        !< Random search vector
  FLOAT,   allocatable :: tridiagonal(:, :)               !< Tridiagonal matrix
  R_TYPE,  allocatable :: f_residual(:, :)                !< Residual vector
  FLOAT                :: norm_tri                        !< Norm of tridiagonal vector
  integer              :: size_t                          !< Size of tridiagonal matrix

  PUSH_SUB(X(upper_bound_estimator))

  ASSERT(n_iterations > 0)
  ASSERT(n_iterations <= MAX_LANCZOS_STEPS)
  SAFE_ALLOCATE(vs(1:mesh%np_part, 1:hm%d%dim))
  call states_elec_generate_random_vector(mesh, st, vs, normalized=.true.)

  size_t = n_iterations
  SAFE_ALLOCATE(tridiagonal(1:size_t, 1:size_t))
  SAFE_ALLOCATE(f_residual(1:mesh%np, 1:hm%d%dim))

  call X(lanczos_tridiagonal)(namespace, mesh, hm, ik, vs, tridiagonal, f_residual, size_t)
  call lalg_matrix_norm2(size_t, size_t, tridiagonal, norm_tri)
  upper_bound = norm_tri + X(mf_nrm2)(mesh, hm%d%dim, f_residual)

  SAFE_DEALLOCATE_A(vs)
  SAFE_DEALLOCATE_A(tridiagonal)
  SAFE_DEALLOCATE_A(f_residual)
  POP_SUB(X(upper_bound_estimator))
end function X(upper_bound_estimator)


!> @brief Estimate the lower and upper bounds of filter interval to discard.
!!
!! Use a k-Lanczos algorithm to compute a tridiagonal matrix, T.
!! Diagonalise T and use min/max eigenvalues to get an approximation
!! of the min/max eigenvalues of H.
!!
!! NOTE: Not extended to work with orbital-dependent functionals.
subroutine X(filter_bounds_estimator)(namespace, mesh, hm, ik, n_iterations, beta_mixing, search_vector, bounds)
  type(namespace_t),         intent(in)    :: namespace               !< Calling namespace
  class(mesh_t),             intent(in)    :: mesh                    !< Real-space mesh
  class(hamiltonian_elec_t), intent(in)    :: hm                      !< Hamiltonian
  integer,                   intent(in)    :: ik                      !< k-point index
  integer,                   intent(in)    :: n_iterations            !< Number of Lanczos iterations, and consequentally the dimensions
  !                                                                   !< of the resulting tridiagonal matrix
  FLOAT,                     intent(in)    :: beta_mixing             !< Linear mixing parameter for the definition of the
  !                                                                   !< lower bound, defined as \( \beta E_{min} + (1 - \beta)  E_{max} \)
  R_TYPE, contiguous,        intent(inout) :: search_vector(:, :)     !< Random search vector
  type(chebyshev_filter_bounds_t), pointer, intent(out) :: bounds     !< Bound estimates for the Chebyshev filter

  FLOAT,    allocatable :: tridiagonal(:, :)                          !< Tridiagonal matrix
  R_TYPE,   allocatable :: f_residual(:, :)                           !< Residual vector

  FLOAT,    allocatable :: eigenvalues(:)                             !< Eigenvalues of tridiagonal matrix
  FLOAT                 :: e_min, e_max                               !< min and max eigenvalues of tridiagonal matrix
  integer               :: size_t                                     !< Size of tridiagonal matrix
  FLOAT                 :: norm_tri                                   !< Norm of tridiagonal vectors
  FLOAT                 :: upper_bound                                !< Estimated upper bound of H

  PUSH_SUB(X(filter_bounds_estimator))
  ASSERT(beta_mixing >= M_ZERO)
  ASSERT(beta_mixing <= M_ONE)

  ! Ensure random vector is normalised
  call X(mf_normalize)(mesh, hm%d%dim, search_vector)

  size_t = n_iterations
  SAFE_ALLOCATE(tridiagonal(1:size_t, 1:size_t))
  SAFE_ALLOCATE(f_residual(1:mesh%np, 1:hm%d%dim))
  call X(lanczos_tridiagonal)(namespace, mesh, hm, ik, search_vector, tridiagonal, f_residual, size_t)
  call lalg_matrix_norm2(size_t, size_t, tridiagonal, norm_tri)
  upper_bound = norm_tri + X(mf_nrm2)(mesh, hm%d%dim, f_residual)
  SAFE_DEALLOCATE_A(f_residual)

  SAFE_ALLOCATE(eigenvalues(1:size_t))
  call lalg_eigensolve(size_t, tridiagonal, eigenvalues)
  e_min = minval(eigenvalues, dim=1)
  e_max = maxval(eigenvalues, dim=1)
  bounds => chebyshev_filter_bounds_t( &
    beta_mixing * e_min + (M_ONE - beta_mixing) * e_max, &
    max(e_max, upper_bound), &
    a_l = e_min)

  SAFE_DEALLOCATE_A(tridiagonal)
  SAFE_DEALLOCATE_A(eigenvalues)

  POP_SUB(X(filter_bounds_estimator))
end subroutine X(filter_bounds_estimator)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
