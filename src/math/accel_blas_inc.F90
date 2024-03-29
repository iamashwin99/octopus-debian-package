!! Copyright (C) 2016 X. Andrade
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

subroutine X(accel_herk)(uplo, trans, n, k, alpha, a, offa, lda, beta, c, offc, ldc)
  integer,           intent(in)    :: uplo
  integer,           intent(in)    :: trans
  integer(i8),       intent(in)    :: n
  integer(i8),       intent(in)    :: k
  real(r8),           intent(in)    :: alpha
  type(accel_mem_t), intent(in)    :: a
  integer(i8),       intent(in)    :: offa
  integer(i8),       intent(in)    :: lda
  real(r8),           intent(in)    :: beta
  type(accel_mem_t), intent(inout) :: c
  integer(i8),       intent(in)    :: offc
  integer(i8),       intent(in)    :: ldc

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  integer :: ierr
#endif
#ifdef HAVE_CUDA
  type(accel_mem_t) :: alpha_buffer, beta_buffer
  type(c_ptr) :: A_with_offset, C_with_offset
#endif

  PUSH_SUB(X(accel_herk))

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
#ifdef R_TREAL
  call clblasDsyrkEx(order = clblasColumnMajor, uplo = uplo, transA = trans, N = n, K = k, &
    alpha = alpha, A = a%mem, offA = offa, lda = lda, &
    beta = beta, C = c%mem, offC = offc, ldc = ldc, &
    CommandQueue = accel%command_queue, status = ierr)
  if (ierr /= clblasSuccess) call clblas_print_error(ierr, 'clblasDsyrkEx')
#else
  call clblasZherkEx(order = clblasColumnMajor, uplo = uplo, transA = trans, N = n, K = k, &
    alpha = alpha, A = a%mem, offA = offa, lda = lda, &
    beta = beta, C = c%mem, offC = offc, ldc = ldc, &
    CommandQueue = accel%command_queue, status = ierr)
  if (ierr /= clblasSuccess) call clblas_print_error(ierr, 'clblasZherkEx')
#endif
  call accel_finish()
#else

#endif

#ifdef HAVE_CUDA

  A_with_offset = X(accel_get_pointer_with_offset)(a%mem, offa)
  C_with_offset = X(accel_get_pointer_with_offset)(c%mem, offc)

  call accel_create_buffer(alpha_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 1)
  call accel_create_buffer(beta_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 1)

  call accel_write_buffer(alpha_buffer, alpha)
  call accel_write_buffer(beta_buffer, beta)

#ifdef R_TREAL
  call cuda_blas_dsyrk(handle = accel%cublas_handle, uplo = uplo, trans = trans, &
    n = n, k = k, &
    alpha = alpha_buffer%mem, &
    A = A_with_offset, lda = lda, &
    beta = beta_buffer%mem, &
    C = C_with_offset, ldc = ldc)
#else
  call cuda_blas_zherk(handle = accel%cublas_handle, uplo = uplo, trans = trans, &
    n = n, k = k, &
    alpha = alpha_buffer%mem, &
    A = A_with_offset, lda = lda, &
    beta = beta_buffer%mem, &
    C = C_with_offset, ldc = ldc)
#endif

  call accel_finish()

  call accel_release_buffer(alpha_buffer)
  call accel_release_buffer(beta_buffer)

  call accel_clean_pointer(A_with_offset)
  call accel_clean_pointer(C_with_offset)
#endif

  POP_SUB(X(accel_herk))
end subroutine X(accel_herk)

! -----------------------------------------------------------------------------------

subroutine X(accel_trsm)(side, uplo, trans, diag, m, n, alpha, a, offa, lda, b, offb, ldb)
  integer,           intent(in)    :: side
  integer,           intent(in)    :: uplo
  integer,           intent(in)    :: trans
  integer,           intent(in)    :: diag
  integer(i8),       intent(in)    :: m
  integer(i8),       intent(in)    :: n
  R_TYPE,            intent(in)    :: alpha
  type(accel_mem_t), intent(inout) :: A
  integer(i8),       intent(in)    :: offa
  integer(i8),       intent(in)    :: lda
  type(accel_mem_t), intent(inout) :: B
  integer(i8),       intent(in)    :: offb
  integer(i8),       intent(in)    :: ldb

#ifdef HAVE_CUDA
  type(accel_mem_t) :: alpha_buffer
  type(c_ptr) :: A_with_offset, B_with_offset
#endif
#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  integer :: ierr
#endif

  PUSH_SUB(X(accel_trsm))

#ifdef HAVE_CUDA

  A_with_offset = X(accel_get_pointer_with_offset)(a%mem, offa)
  B_with_offset = X(accel_get_pointer_with_offset)(b%mem, offb)

  call accel_create_buffer(alpha_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, 1)
  call accel_write_buffer(alpha_buffer, alpha)

  call aX(cuda_blas_, trsm)(handle = accel%cublas_handle, side = side, uplo = uplo, trans = trans, diag = diag, &
    m = m, n = n, alpha = alpha_buffer%mem, A = A_with_offset, lda = lda, B = B_with_offset, ldb = ldb)
#endif

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  call aX(clblas,trsmEx)(order = clblasColumnMajor, side = side, uplo = uplo, transA = trans, diag = diag, &
    M = m, N = n, alpha = alpha, A = a%mem, offA = offa, lda = lda, &
    B = b%mem, offB = offb, ldb = ldb, &
    CommandQueue = accel%command_queue, status = ierr)
  if (ierr /= clblasSuccess) call clblas_print_error(ierr, TOSTRING(aX(clblas,trsmEx)))
#endif

  call accel_finish()

#ifdef HAVE_CUDA
  call accel_release_buffer(alpha_buffer)
  call accel_clean_pointer(A_with_offset)
  call accel_clean_pointer(B_with_offset)
#endif

  POP_SUB(X(accel_trsm))
end subroutine X(accel_trsm)

! -----------------------------------------------------------------------------------

subroutine X(accel_gemm)(transa, transb, m, n, k, alpha, A, offa, lda, B, offb, ldb, beta, C, offc, ldc)
  integer,            intent(in)    :: transa
  integer,            intent(in)    :: transb
  integer(i8),        intent(in)    :: m
  integer(i8),        intent(in)    :: n
  integer(i8),        intent(in)    :: k
  R_TYPE,             intent(in)    :: alpha
  type(accel_mem_t),  intent(in)    :: A
  integer(i8),        intent(in)    :: offa
  integer(i8),        intent(in)    :: lda
  type(accel_mem_t),  intent(in)    :: B
  integer(i8),        intent(in)    :: offb
  integer(i8),        intent(in)    :: ldb
  R_TYPE,             intent(in)    :: beta
  type(accel_mem_t),  intent(inout) :: C
  integer(i8),        intent(in)    :: offc
  integer(i8),        intent(in)    :: ldc

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  integer :: ierr
#endif
#ifdef HAVE_CUDA
  type(accel_mem_t) :: alpha_buffer, beta_buffer
  type(c_ptr) :: A_with_offset, B_with_offset, C_with_offset
#endif

  PUSH_SUB(X(accel_gemm))

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)

  call aX(clblas,gemmEx)(order = clblasColumnMajor, transA = transa, transB = transb, &
    M = m, N = n, K = k, &
    alpha = alpha, A = a%mem, offa = offa, lda = lda, &
    B = b%mem, offB = offb, ldb = ldb, &
    beta = beta, C = c%mem, offC = offc, ldc = ldc, &
    CommandQueue = accel%command_queue, status = ierr)
  if (ierr /= clblasSuccess) call clblas_print_error(ierr, 'clblasXgemmEx')

#endif
#ifdef HAVE_CUDA

  A_with_offset = X(accel_get_pointer_with_offset)(a%mem, offa)
  B_with_offset = X(accel_get_pointer_with_offset)(b%mem, offb)
  C_with_offset = X(accel_get_pointer_with_offset)(c%mem, offc)

  call accel_create_buffer(alpha_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, 1)
  call accel_create_buffer(beta_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, 1)

  call accel_write_buffer(alpha_buffer, alpha)
  call accel_write_buffer(beta_buffer, beta)

  call aX(cuda_blas_,gemm)(handle = accel%cublas_handle, transa = transa, transb = transb, &
    m = m, n = n, k = k, &
    alpha = alpha_buffer%mem, a = A_with_offset, lda = lda, &
    b = B_with_offset, ldb = ldb, &
    beta = beta_buffer%mem, c = C_with_offset, ldc = ldc)

  call accel_finish()

  call accel_release_buffer(alpha_buffer)
  call accel_release_buffer(beta_buffer)

  call accel_clean_pointer(A_with_offset)
  call accel_clean_pointer(B_with_offset)
  call accel_clean_pointer(C_with_offset)
#endif

  POP_SUB(X(accel_gemm))
end subroutine X(accel_gemm)

! -----------------------------------------------------------------

subroutine X(accel_dot)(n, x, offx, incx, y, offy, incy, res, offres)
  integer(i8),       intent(in)    :: n
  type(accel_mem_t), intent(in)    :: x
  integer(i8),       intent(in)    :: offx
  integer(i8),       intent(in)    :: incx
  type(accel_mem_t), intent(in)    :: y
  integer(i8),       intent(in)    :: offy
  integer(i8),       intent(in)    :: incy
  type(accel_mem_t), intent(inout) :: res
  integer(i8),       intent(in)    :: offres

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  integer :: status
  type(accel_mem_t)  :: scratch_buffer
#endif

  PUSH_SUB(X(accel_dot))

#ifdef HAVE_CUDA
#ifdef R_TREAL
  call cuda_blas_ddot &
#else
    call cuda_blas_zdotc &
#endif
    (handle = accel%cublas_handle, n = n, x = x%mem, offx = offx, incx = incx, &
    y = y%mem, offy = offy, incy = incy, res = res%mem, offres = offres)
#endif

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  call accel_create_buffer(scratch_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, n)

#ifdef R_TREAL
  call clblasDdot &
#else
    call clblasZdotc &
#endif
    (N = n, dotProduct = res%mem, offDP = offres, &
    X = x%mem, offx = offx, incx = int(incx), &
    Y = y%mem, offy = offy, incy = int(incy), &
    scratchBuff = scratch_buffer%mem, CommandQueue = accel%command_queue, status = status)
  if (status /= clblasSuccess) call clblas_print_error(status, 'clblasXdot')

  call accel_finish()

  call accel_release_buffer(scratch_buffer)
#endif

  POP_SUB(X(accel_dot))
end subroutine X(accel_dot)

! -----------------------------------------------------------------

subroutine X(accel_dotu)(n, x, offx, incx, y, offy, incy, res, offres)
  integer(i8),       intent(in)    :: n
  type(accel_mem_t), intent(in)    :: x
  integer(i8),       intent(in)    :: offx
  integer(i8),       intent(in)    :: incx
  type(accel_mem_t), intent(in)    :: y
  integer(i8),       intent(in)    :: offy
  integer(i8),       intent(in)    :: incy
  type(accel_mem_t), intent(inout) :: res
  integer(i8),       intent(in)    :: offres

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  integer :: status
  type(accel_mem_t)  :: scratch_buffer
#endif

  PUSH_SUB(X(accel_dotu))

#ifdef HAVE_CUDA
#ifdef R_TREAL
  call cuda_blas_ddot &
#else
    call cuda_blas_zdotu &
#endif
    (handle = accel%cublas_handle, n = n, x = x%mem, offx = offx, incx = incx, &
    y = y%mem, offy = offy, incy = incy, res = res%mem, offres = offres)
#endif

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  call accel_create_buffer(scratch_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, n)

#ifdef R_TREAL
  call clblasDdot &
#else
    call clblasZdotu &
#endif
    (N = n, dotProduct = res%mem, offDP = offres, &
    X = x%mem, offx = offx, incx = int(incx), &
    Y = y%mem, offy = offy, incy = int(incy), &
    scratchBuff = scratch_buffer%mem, CommandQueue = accel%command_queue, status = status)
  if (status /= clblasSuccess) call clblas_print_error(status, 'clblasXdot')

  call accel_finish()

  call accel_release_buffer(scratch_buffer)
#endif

  POP_SUB(X(accel_dotu))
end subroutine X(accel_dotu)


! -----------------------------------------------------------------

subroutine X(accel_nrm2)(n, x, offx, incx, res, offres)
  integer(i8),       intent(in)    :: n
  type(accel_mem_t), intent(in)    :: x
  integer(i8),       intent(in)    :: offx
  integer(i8),       intent(in)    :: incx
  type(accel_mem_t), intent(inout) :: res
  integer(i8),       intent(in)    :: offres

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  integer :: status
  type(accel_mem_t)  :: scratch_buffer
#endif

  PUSH_SUB(X(accel_nrm2))

#ifdef HAVE_CUDA
  call aX(cuda_blas_,nrm2)(handle = accel%cublas_handle, n = n, x = x%mem, offx = offx, incx = incx, &
    res = res%mem, offres = offres)
#endif


#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  call accel_create_buffer(scratch_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, n)

#ifdef R_TREAL
  call clblasDnrm2 &
#else
    call clblasDznrm2 &
#endif
    (N = n, NRM2 = res%mem, offNRM2 = offres, &
    X = x%mem, offx = offx, incx = int(incx), &
    scratchBuff = scratch_buffer%mem, CommandQueue = accel%command_queue, status = status)
  if (status /= clblasSuccess) call clblas_print_error(status, 'clblasXnrm2')

  call accel_finish()

  call accel_release_buffer(scratch_buffer)
#endif

  POP_SUB(X(accel_nrm2))
end subroutine X(accel_nrm2)

! -----------------------------------------------------------------------------------

subroutine X(accel_gemv)(transa, m, n, alpha, A, lda, x, incx, beta, y, incy)
  integer,            intent(in)    :: transa
  integer(i8),        intent(in)    :: m
  integer(i8),        intent(in)    :: n
  R_TYPE,             intent(in)    :: alpha
  type(accel_mem_t),  intent(in)    :: A
  integer(i8),        intent(in)    :: lda
  type(accel_mem_t),  intent(in)    :: x
  integer(i8),        intent(in)    :: incx
  R_TYPE,             intent(in)    :: beta
  type(accel_mem_t),  intent(inout) :: y
  integer(i8),        intent(in)    :: incy

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)
  integer :: ierr
#endif
#ifdef HAVE_CUDA
  type(accel_mem_t) :: alpha_buffer, beta_buffer
#endif

  PUSH_SUB(X(accel_gemv))

#if defined(HAVE_CLBLAS) || defined(HAVE_CLBLAST)

  call aX(clblas,gemvEx)(order = clblasColumnMajor, transA = transa, &
    M = m, N = n,  &
    alpha = alpha, A = a%mem, offa = 0_i8, lda = lda, &
    x = x%mem, offx = 0_i8, incx = incx, &
    beta = beta, y = y%mem, offy = 0_i8, incy = incy, &
    CommandQueue = accel%command_queue, status = ierr)
  if (ierr /= clblasSuccess) call clblas_print_error(ierr, 'clblasXgemvEx')

#endif
#ifdef HAVE_CUDA

  call accel_create_buffer(alpha_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, 1)
  call accel_create_buffer(beta_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, 1)

  call accel_write_buffer(alpha_buffer, alpha)
  call accel_write_buffer(beta_buffer, beta)

  call aX(cuda_blas_,gemv)(handle = accel%cublas_handle, transa = transa, &
    m = m, n = n,  &
    alpha = alpha_buffer%mem, &
    a = a%mem, lda = lda, &
    x = x%mem, incx = incx, &
    beta = beta_buffer%mem, &
    y = y%mem, incy = incy)

  call accel_finish()

  call accel_release_buffer(alpha_buffer)
  call accel_release_buffer(beta_buffer)

#endif

  POP_SUB(X(accel_gemv))
end subroutine X(accel_gemv)
