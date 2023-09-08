/* Copyright (C) 2022 N. Tancogne-Dejean
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301, USA.
*/

#include <config.h>

#ifdef HAVE_CLBLAST

#include <clblast_c.h>

void FC_FUNC_(clblasdtrsmex_low, CLBLASDTRSMEX_LOW)(int * order,
							  int * side,
							  int * uplo,
							  int * transA,
							  int * diag,
							  cl_long * M,
							  cl_long * N,
							  double * alpha,
							  const cl_mem * A,
							  size_t * offA,
							  size_t * lda,
							  cl_mem * B,
							  size_t * offB,
							  size_t * ldb, 
							  cl_command_queue * CommandQueue, 
							  int * status){


  *status = CLBlastDtrsm((CLBlastLayout) *order, (CLBlastSide) *side, (CLBlastTriangle) *uplo, 
			(CLBlastTranspose) *transA, (CLBlastDiagonal) *diag,
			(size_t) *M, (size_t) *N, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, 
			CommandQueue, NULL);
}


void FC_FUNC_(clblasztrsmex_low, CLBLASZTRSMEX_LOW)(int * order,
							  int * side,
							  int * uplo,
							  int * transA,
							  int * diag,
							  cl_long * M,
							  cl_long * N,
							  cl_double2 * alpha,
							  const cl_mem * A,
							  size_t * offA,
							  size_t * lda,
							  cl_mem * B,
							  size_t * offB,
							  size_t * ldb, 
							  cl_command_queue * CommandQueue, 
							  int * status){


  *status = CLBlastZtrsm((CLBlastLayout) *order, (CLBlastSide) *side, (CLBlastTriangle) *uplo, 
			(CLBlastTranspose) *transA, (CLBlastDiagonal) *diag,
			(size_t) *M, (size_t) *N, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, 
			CommandQueue, NULL);
}

void FC_FUNC_(clblasdgemvex_low, CLBLASDGEMVEX_LOW)(int * order,
                                                          int * transA,
                                                          cl_long * M,
                                                          cl_long * N,
                                                          double * alpha,
                                                          const cl_mem * A,
                                                          cl_long * offA,
                                                          cl_long * lda,
                                                          const cl_mem * X,
                                                          cl_long * offX,
                                                          cl_long * incx,
                                                          double * beta,
                                                          cl_mem * Y,
                                                          cl_long * offY,
                                                          cl_long * incy,
                                                          cl_command_queue * CommandQueue,
                                                          int * status){

  *status = CLBlastDgemv((CLBlastLayout) *order, (CLBlastTranspose) *transA,
                        (size_t) *M, (size_t) *N, *alpha,
                        *A, (size_t) *offA, (size_t) *lda,
                        *X, (size_t) *offX, (size_t) *incx, *beta,
                        *Y, (size_t) *offY, (size_t) *incy,
                        CommandQueue, NULL);
}

void FC_FUNC_(clblaszgemvex_low, CLBLASZGEMVEX_LOW)(int * order,
                                                          int * transA,
                                                          cl_long * M,
                                                          cl_long * N,
                                                          cl_double2 * alpha,
                                                          const cl_mem * A,
                                                          cl_long * offA,
                                                          cl_long * lda,
                                                          const cl_mem * X,
                                                          cl_long * offX,
                                                          cl_long * incx,
                                                          cl_double2 * beta,
                                                          cl_mem * Y,
                                                          cl_long * offY,
                                                          cl_long * incy,
                                                          cl_command_queue * CommandQueue,
                                                          int * status){

  *status = CLBlastZgemv((CLBlastLayout) *order, (CLBlastTranspose) *transA,
                        (size_t) *M, (size_t) *N, *alpha,
                        *A, (size_t) *offA, (size_t) *lda,
                        *X, (size_t) *offX, (size_t) *incx, *beta,
                        *Y, (size_t) *offY, (size_t) *incy,
                        CommandQueue, NULL);
}



void FC_FUNC_(clblasdgemmex_low, CLBLASDGEMMEX_LOW)(int * order,
							  int * transA, 
							  int * transB, 
							  cl_long * M,
							  cl_long * N,
							  cl_long * K,
							  double * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  const cl_mem * B,
							  cl_long * offB,
							  cl_long * ldb, 
							  double * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = CLBlastDgemm((CLBlastLayout) *order, (CLBlastTranspose) *transA, (CLBlastTranspose) *transB, 
			(size_t) *M, (size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			CommandQueue, NULL);
}

void FC_FUNC_(clblaszgemmex_low, CLBLASDGEMMEX_LOW)(int * order,
							  int * transA, 
							  int * transB, 
							  cl_long * M,
							  cl_long * N,
							  cl_long * K,
							  cl_double2 * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  const cl_mem * B,
							  cl_long * offB,
							  cl_long * ldb, 
							  cl_double2 * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = CLBlastZgemm((CLBlastLayout) *order, (CLBlastTranspose) *transA, (CLBlastTranspose) *transB, 
			(size_t) *M, (size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			CommandQueue, NULL);
}

void FC_FUNC_(clblasdsyrkex_low, CLBLASDSYRKEX_LOW)(int * order,
							  int * uplo, 
							  int * transA, 
							  cl_long * N,
							  cl_long * K,
							  double * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  double * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = CLBlastDsyrk((CLBlastLayout) *order, (CLBlastTriangle) *uplo, (CLBlastTranspose) *transA, 
			(size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			CommandQueue, NULL);
}

void FC_FUNC_(clblaszherkex_low, CLBLASZHERKEX_LOW)(int * order,
							  int * uplo, 
							  int * transA, 
							  cl_long * N,
							  cl_long * K,
							  double * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  double * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = CLBlastZherk((CLBlastLayout) *order, (CLBlastTriangle) *uplo, (CLBlastTranspose) *transA, 
			(size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			CommandQueue, NULL);
}


CLBlastStatusCode FC_FUNC_(clblasddot_low, CLBLASDDOT_LOW)(cl_long * N,
						      cl_mem * dotProduct,
						      cl_long * offDP,
						      const cl_mem *X,
						      cl_long * offx,
						      int *incx,
						      const cl_mem * Y,
						      cl_long * offy,
						      int * incy,
						      cl_mem * scratchBuff,
						      cl_command_queue * CommandQueue,
						      int * status){

  
  *status = CLBlastDdot((size_t) *N, *dotProduct, (size_t) *offDP, *X, (size_t) *offx, *incx,
		       *Y, (size_t) *offy, *incy, CommandQueue, NULL);

}


CLBlastStatusCode FC_FUNC_(clblaszdotc_low, CLBLASZDOTc_LOW)(cl_long * N,
                                                        cl_mem * dotProduct,
                                                        cl_long * offDP,
                                                        const cl_mem *X,
                                                        cl_long * offx,
                                                        int *incx,
                                                        const cl_mem * Y,
                                                        cl_long * offy,
                                                        int * incy,
                                                        cl_mem * scratchBuff,
                                                        cl_command_queue * CommandQueue,
                                                        int * status){


  *status = CLBlastZdotc((size_t) *N, *dotProduct, (size_t) *offDP, *X, (size_t) *offx, *incx,
                        *Y, (size_t) *offy, *incy, CommandQueue, NULL);

}


CLBlastStatusCode FC_FUNC_(clblaszdotu_low, CLBLASZDOTU_LOW)(cl_long * N,
							cl_mem * dotProduct,
							cl_long * offDP,
							const cl_mem *X,
							cl_long * offx,
							int *incx,
							const cl_mem * Y,
							cl_long * offy,
							int * incy,
							cl_mem * scratchBuff,
							cl_command_queue * CommandQueue,
							int * status){

  
  *status = CLBlastZdotu((size_t) *N, *dotProduct, (size_t) *offDP, *X, (size_t) *offx, *incx,
			*Y, (size_t) *offy, *incy, CommandQueue, NULL);

}


CLBlastStatusCode FC_FUNC_(clblasdnrm2_low, CLBLASDNRM2_LOW)(cl_long * N,
							cl_mem * NRM2,
							cl_long * offNRM2,
							const cl_mem *X,
							cl_long * offx,
							int *incx,
							cl_mem * scratchBuff,
							cl_command_queue * CommandQueue,
							int * status){

  
  *status = CLBlastDnrm2((size_t) *N, *NRM2, (size_t) *offNRM2, *X, (size_t) *offx, *incx,
			CommandQueue, NULL);

}

CLBlastStatusCode FC_FUNC_(clblasdznrm2_low, CLBLASDZNRM2_LOW)(cl_long * N,
							  cl_mem * NRM2,
							  cl_long * offNRM2,
							  const cl_mem *X,
							  cl_long * offx,
							  int *incx,
							  cl_mem * scratchBuff,
							  cl_command_queue * CommandQueue,
							  int * status){
  
  
  *status = CLBlastDznrm2((size_t) *N, *NRM2, (size_t) *offNRM2, *X, (size_t) *offx, *incx,
			 CommandQueue, NULL);
  
}

#endif

