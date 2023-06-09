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
subroutine controlfunction_basis_to_theta(par)
  type(controlfunction_t), intent(inout) :: par
  integer :: j, n, dof
  FLOAT, allocatable :: ep(:), e(:), y(:), x(:)

  PUSH_SUB(controlfunction_basis_to_theta)

  ASSERT(par%current_representation /= ctr_internal)

  select case (par%current_representation)
  case (ctr_rt)
    do j = 1, par%dim
      par%theta(j) = tdf(par%f(1), j)
    end do

  case (ctr_fourier_series_h)
    n = par%dim
    dof = par%dof
    SAFE_ALLOCATE( e(1:n))
    SAFE_ALLOCATE(ep(1:n))
    SAFE_ALLOCATE(x(1:dof))

    do j = 1, n
      ep(j) = tdf(par%f(1), j)
    end do
    e = matmul(par%utransf, ep)
    call cartesian2hyperspherical(e, x(1:n-1))

    par%theta = x
    SAFE_DEALLOCATE_A(e)
    SAFE_DEALLOCATE_A(ep)
    SAFE_DEALLOCATE_A(x)

  case (ctr_zero_fourier_series_h)
    n = par%dim
    dof = par%dof
    SAFE_ALLOCATE(ep(1:n))
    SAFE_ALLOCATE(x(1:dof))

    do j = 1, n
      ep(j) = tdf(par%f(1), j)
    end do
    SAFE_ALLOCATE(y(1:n-1))
    y = matmul(par%utransf, ep(2:n))
    call cartesian2hyperspherical(y, x(1:n-2))
    SAFE_DEALLOCATE_A(y)

    par%theta = x
    SAFE_DEALLOCATE_A(ep)
    SAFE_DEALLOCATE_A(x)

  case (ctr_fourier_series)
    do j = 1, par%dim
      par%theta(j) = tdf(par%f(1), j)
    end do

  case (ctr_zero_fourier_series)
    ! In this case, the transformation is (n = par%dim):
    ! theta(1) = a(2)
    ! theta(2) = a(3)
    ! ...      = ...
    ! theta(n/2-1)   = a(n/2)
    ! theta(n/2) = b(1)
    ! ...      = ...
    ! theta(n-1) = b(n/2)
    ! where a are the coefficients of the cosines in the Fourier series, and b are
    ! the coefficients of the sines

    SAFE_ALLOCATE(e(1:par%dim))
    do j = 1, par%dim
      e(j) = tdf(par%f(1), j)
    end do

    do j = 2, par%dim
      par%theta(j-1) = e(j)
    end do

    SAFE_DEALLOCATE_A(e)

  end select

  POP_SUB(controlfunction_basis_to_theta)
end subroutine controlfunction_basis_to_theta
 ! ---------------------------------------------------------


 ! ---------------------------------------------------------
subroutine controlfunction_theta_to_basis(par)
  type(controlfunction_t), intent(inout) :: par

  FLOAT, allocatable :: a(:), e(:), ep(:), x(:)
  integer :: n, dof, j

  PUSH_SUB(controlfunction_theta_to_basis)

  ASSERT(par%current_representation /= ctr_internal)

  select case (par%current_representation)
  case (ctr_rt)
    call tdf_set_numerical(par%f(1), par%theta)

  case (ctr_fourier_series_h)

    n = par%dim
    dof = par%dof
    SAFE_ALLOCATE( e(1:n))
    SAFE_ALLOCATE(ep(1:n))
    SAFE_ALLOCATE(x(1:dof))
    x = par%theta

    call hyperspherical2cartesian(x(1:n-1), e)
    e = sqrt(cf_common%targetfluence) * e
    ep = matmul(par%utransfi, e)
    call tdf_set_numerical(par%f(1), ep)

    SAFE_DEALLOCATE_A(e)
    SAFE_DEALLOCATE_A(ep)
    SAFE_DEALLOCATE_A(x)

  case (ctr_zero_fourier_series_h)

    n = par%dim
    dof = par%dof
    SAFE_ALLOCATE( e(1:n))
    SAFE_ALLOCATE( ep(1:n))
    SAFE_ALLOCATE(x(1:dof))
    x = par%theta

    call hyperspherical2cartesian(x, ep(2:n))
    e(2:n) = matmul(par%utransfi, ep(2:n))
    e(1) = -sum(e(2:n/2))
    e = sqrt(cf_common%targetfluence) * e
    call tdf_set_numerical(par%f(1), e)

    SAFE_DEALLOCATE_A(ep)
    SAFE_DEALLOCATE_A(e)
    SAFE_DEALLOCATE_A(x)

  case (ctr_fourier_series)

    call tdf_set_numerical(par%f(1), par%theta)

  case (ctr_zero_fourier_series)

    ! In this case, the transformation is (n = par%dim):
    ! a(1) = -sum(a(2)...a(n/2))  represents the constraint
    ! theta(1) = a(2)
    ! theta(2) = a(3)
    ! ...      = ...
    ! theta(n/2-1)   = a(n/2)
    ! theta(n/2) = b(1)
    ! ...      = ...
    ! theta(n-1) = b(n/2)
    ! where a are the coefficients of the cosines in the Fourier series, and b are
    ! the coefficients of the sines (Note that in the next lines "a" are the coefficients
    ! of both sines and cosines.

    n = par%dim
    SAFE_ALLOCATE(a(1:n))

    do j = 2, n
      a(j) = par%theta(j-1)
    end do
    a(1) = -sum(a(2:n/2))
    call tdf_set_numerical(par%f(1), a)

    SAFE_DEALLOCATE_A(a)

  end select

  POP_SUB(controlfunction_theta_to_basis)
end subroutine controlfunction_theta_to_basis
 ! ---------------------------------------------------------


 ! ---------------------------------------------------------
subroutine controlfunction_get_theta(par, theta)
  type(controlfunction_t), intent(in) :: par
  FLOAT, intent(inout) :: theta(:)

  PUSH_SUB(controlfunction_get_theta)
  theta = par%theta

  POP_SUB(controlfunction_get_theta)
end subroutine controlfunction_get_theta
 ! ---------------------------------------------------------


 ! ---------------------------------------------------------
subroutine controlfunction_set_theta(par, theta)
  type(controlfunction_t), intent(inout) :: par
  FLOAT, intent(in) :: theta(:)

  PUSH_SUB(controlfunction_set_theta)
  par%theta = theta

  POP_SUB(controlfunction_set_theta)
end subroutine controlfunction_set_theta
 ! ---------------------------------------------------------


 ! ---------------------------------------------------------
subroutine controlfunction_trans_matrix(par)
  type(controlfunction_t), intent(inout) :: par

  integer :: i, mm, nn, n, j, k
  FLOAT :: t, w1, dt
  FLOAT, allocatable :: neigenvec(:, :), eigenvec(:, :), eigenval(:)

  type(tdf_t), allocatable :: fnn(:)

  PUSH_SUB(controlfunction_trans_matrix)

  ! First, we will construct the matrix u, from which the fluence is computed
  ! as dot_product(x, matmul(u, x)).
  select case (cf_common%representation)
  case (ctr_rt)

    ! Since it is diagonal, instead of par%u being a matrix, we only store the diagonal.
    SAFE_ALLOCATE(par%u(1:par%dim, 1))
    par%u = M_ZERO

    dt = tdf_dt(par%f(1))
    par%u(1, 1) = M_HALF * dt
    do mm = 2, par%dim - 1
      par%u(mm, 1) = dt
    end do
    par%u(par%dim, 1) = M_HALF * dt

  case (ctr_fourier_series, ctr_fourier_series_h)

    SAFE_ALLOCATE(par%u(1:par%dim, 1:par%dim))
    par%u = diagonal_matrix(par%dim, M_ONE)

    if (cf_common%mode == controlfunction_mode_f) then

      SAFE_ALLOCATE(fnn(1:par%dim))
      do mm = 1, par%dim
        call tdf_init_numerical(fnn(mm), tdf_niter(par%f(1)), tdf_dt(par%f(1)), cf_common%omegamax, rep = TDF_FOURIER_SERIES)
        call tdf_set_numerical(fnn(mm), mm, M_ONE)
        call tdf_fourier_to_numerical(fnn(mm))
        do i = 1, tdf_niter(fnn(mm)) + 1
          t = (i-1)*tdf_dt(fnn(mm))
          call tdf_set_numerical(fnn(mm), i, tdf(fnn(mm), i)*cos(par%w0*t))
        end do
      end do

      do mm = 1, par%dim
        do nn = mm, par%dim
          par%u(mm, nn) = tdf_dot_product(fnn(mm), fnn(nn))
        end do
      end do

      do mm = 1, par%dim
        call tdf_end(fnn(mm))
        do nn = 1, mm - 1
          par%u(mm, nn) = par%u(nn, mm)
        end do
      end do

      SAFE_DEALLOCATE_A(fnn)

    end if

  case (ctr_zero_fourier_series, ctr_zero_fourier_series_h)

    n = par%dim - 1
    SAFE_ALLOCATE(par%u (1:n, 1:n))
    par%u = M_ZERO
    w1 = (M_TWO*M_PI/(tdf_dt(par%f(1))*tdf_niter(par%f(1))))

    SAFE_ALLOCATE(fnn(1:n))
    do mm = 1, n
      call tdf_init_numerical(fnn(mm), tdf_niter(par%f(1)), tdf_dt(par%f(1)), cf_common%omegamax, rep = TDF_ZERO_FOURIER)
      call tdf_set_numerical(fnn(mm), mm+1, M_ONE)
      call tdf_fourier_to_numerical(fnn(mm))
      if (mm <= n/2) then
        do i = 1, tdf_niter(fnn(mm)) + 1
          t = (i-1)*tdf_dt(fnn(mm))
          call tdf_set_numerical(fnn(mm), i, &
            (tdf(fnn(mm), i)-sqrt(M_TWO/(tdf_dt(fnn(mm))*tdf_niter(fnn(mm))))*cos(w1*t))*cos(par%w0*t))
        end do
      else
        do i = 1, tdf_niter(fnn(mm)) + 1
          t = (i-1)*tdf_dt(fnn(mm))
          call tdf_set_numerical(fnn(mm), i, tdf(fnn(mm), i)*cos(par%w0*t))
        end do
      end if
    end do

    do mm = 1, n
      do nn = mm, n
        par%u(mm, nn) = tdf_dot_product(fnn(mm), fnn(nn))
      end do
    end do

    do mm = 1, n
      call tdf_end(fnn(mm))
      do nn = 1, mm - 1
        par%u(mm, nn) = par%u(nn, mm)
      end do
    end do

    SAFE_DEALLOCATE_A(fnn)

  end select

  ! Now the transformation matrices, in case they are needed.
  select case (cf_common%representation)

  case (ctr_fourier_series_h)

    n = par%dim
    SAFE_ALLOCATE(par%utransf(1:n, 1:n))
    SAFE_ALLOCATE(par%utransfi(1:n, 1:n))
    SAFE_ALLOCATE(eigenvec(1:n, 1:n))
    SAFE_ALLOCATE(eigenval(1:n))

    eigenvec = par%u
    call lalg_eigensolve(par%dim, eigenvec, eigenval)

    ! We need to make sure that eigenvectors have the same sign on all machines, which is not guaranteed
    ! by LAPACK. So, we will use the following criterion: the sign of the first non-null component should be
    ! positive.
    do nn = 1, par%dim
      do mm = 1, par%dim
        if (eigenvec(mm, nn)*eigenvec(mm, nn) > CNST(1.0e-20)) then
          !eigenvec(1:par%dim, nn) = sign(eigenvec(mm, nn), M_ONE) * eigenvec(1:par%dim, nn)
          if (eigenvec(mm, nn) < M_ZERO) eigenvec(1:par%dim, nn) = - eigenvec(1:par%dim, nn)
          exit
        end if
      end do
    end do

    do mm = 1, par%dim
      do nn = 1, par%dim
        eigenvec(mm, nn) = eigenvec(mm, nn) * sqrt(eigenval(nn))
      end do
    end do

    par%utransf = transpose(eigenvec)
    par%utransfi = par%utransf
    call lalg_inverter(par%dim, par%utransfi)

    SAFE_DEALLOCATE_A(eigenvec)
    SAFE_DEALLOCATE_A(eigenval)


  case (ctr_zero_fourier_series_h)

    n = par%dim
    SAFE_ALLOCATE(eigenvec(1:n-1, 1:n-1))
    SAFE_ALLOCATE(neigenvec(1:n-1, 1:n-1))
    SAFE_ALLOCATE(eigenval(1:n-1))
    SAFE_ALLOCATE(par%utransf (1:n-1, 1:n-1))
    SAFE_ALLOCATE(par%utransfi(1:n-1, 1:n-1))
    par%utransf  = diagonal_matrix(n-1, M_ONE)
    par%utransfi = diagonal_matrix(n-1, M_ONE)

    if (cf_common%mode == controlfunction_mode_f) then

      eigenvec = par%u
      call lalg_eigensolve(n-1, eigenvec, eigenval)
      ! We need to make sure that eigenvectors have the same sign on all machines, which is not guaranteed
      ! by LAPACK. So, we will use the following criterion: the sign of the first non-null component should be
      ! positive.
      do nn = 1, n-1
        do mm = 1, n-1
          if (eigenvec(mm, nn)*eigenvec(mm, nn) > CNST(1.0e-20)) then
            !eigenvec(1:par%dim, nn) = sign(eigenvec(mm, nn), M_ONE) * eigenvec(1:par%dim, nn)
            if (eigenvec(mm, nn) < M_ZERO) eigenvec(1:n-1, nn) = - eigenvec(1:n-1, nn)
            exit
          end if
        end do
      end do

    else

      ! In this case, we can write explicitly the eigenvalues. It is better to do it this way, because there
      ! are degenerate eigenvalues, and different machines will give different eigenvectors (all of them valid). This
      ! would imply different results in different machines, which is not very nice.
      ! k = 1
      eigenvec = M_ZERO
      eigenvec(1:n/2-1, 1) = M_ONE
      eigenvec(n/2:n-1, 1) = M_ZERO
      eigenval(1) = n/2
      ! k = 2, ...., n/2-1
      do k = 2, n/2-1
        eigenval(k) = M_ONE
        eigenvec(:, k) = M_ZERO
        eigenvec(1, k) = -M_ONE
        eigenvec(k, k) = M_ONE
      end do
      ! k = n/2, ..., n-1
      do k = n/2, n-1
        eigenval(k) = M_ONE
        eigenvec(:, k) = M_ZERO
        eigenvec(k, k) = M_ONE
      end do

      do k = 1, n-1
        neigenvec(:, k) = eigenvec(:, k)
        do j = 1, k-1
          neigenvec(:, k) = neigenvec(:, k) - &
            (dot_product(eigenvec(:, k), neigenvec(:, j)) / dot_product(neigenvec(:, j), neigenvec(:, j))) * neigenvec(:, j)
        end do
      end do
      do k = 1, n-1
        neigenvec(:, k) = neigenvec(:, k)/norm2(neigenvec(:, k))
      end do
      eigenvec = neigenvec
    end if


    do mm = 1, n-1
      do nn = 1, n-1
        eigenvec(mm, nn) = eigenvec(mm, nn) * sqrt(eigenval(nn))
      end do
    end do
    par%utransf = transpose(eigenvec)
    par%utransfi = par%utransf
    call lalg_inverter(n-1, par%utransfi)

    SAFE_DEALLOCATE_A(eigenvec)
    SAFE_DEALLOCATE_A(neigenvec)
    SAFE_DEALLOCATE_A(eigenval)

  end select

  POP_SUB(controlfunction_trans_matrix)
end subroutine controlfunction_trans_matrix
 ! ---------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
