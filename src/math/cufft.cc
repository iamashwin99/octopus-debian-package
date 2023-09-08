/*
 Copyright (C) 2016 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <fortran_types.h>

#ifdef HAVE_CUDA
#include <cuda.h>
#include <cufft.h>
#else
typedef int cufftHandle;
typedef int CUdeviceptr;
typedef int cudaStream_t;
#endif

#define CUFFT_SAFE_CALL(x)                                                     \
  do {                                                                         \
    cufftResult safe_call_result = x;                                          \
    if (safe_call_result != CUFFT_SUCCESS) {                                   \
      std::cerr << "\nerror: " #x " failed with error " << safe_call_result    \
                << std::endl;                                                  \
      if (safe_call_result == CUFFT_ALLOC_FAILED) {                            \
        std::cerr << "cuFFT failed to allocate GPU or CPU memory."             \
                  << std::endl;                                                \
        std::cerr << "Please try to use more nodes for the parallelization."   \
                  << std::endl;                                                \
      }                                                                        \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

extern "C" void FC_FUNC_(cuda_fft_plan3d,
                         CUDA_FFT_PLAN3D)(cufftHandle **plan, fint *nx,
                                          fint *ny, fint *nz, fint *type,
                                          cudaStream_t **stream) {
  *plan = new cufftHandle;
#ifdef HAVE_CUDA
  CUFFT_SAFE_CALL(cufftPlan3d(*plan, *nx, *ny, *nz, (cufftType)*type));
  CUFFT_SAFE_CALL(cufftSetStream(**plan, **stream));
#endif
}

extern "C" void FC_FUNC_(cuda_fft_destroy,
                         CUDA_FFT_DESTROY)(cufftHandle **plan) {
#ifdef HAVE_CUDA
  CUFFT_SAFE_CALL(cufftDestroy(**plan));
#endif
  delete *plan;
}

extern "C" void FC_FUNC_(cuda_fft_execute_d2z,
                         CUDA_FFT_EXECUTE_D2Z)(cufftHandle **plan,
                                               CUdeviceptr **idata,
                                               CUdeviceptr **odata,
                                               cudaStream_t **stream) {
#ifdef HAVE_CUDA
  CUFFT_SAFE_CALL(cufftSetStream(**plan, **stream));
  CUFFT_SAFE_CALL(cufftExecD2Z(**plan, (cufftDoubleReal *)**idata,
                               (cufftDoubleComplex *)**odata));
#endif
}

extern "C" void FC_FUNC_(cuda_fft_execute_z2d,
                         CUDA_FFT_EXECUTE_Z2D)(cufftHandle **plan,
                                               CUdeviceptr **idata,
                                               CUdeviceptr **odata,
                                               cudaStream_t **stream) {
#ifdef HAVE_CUDA
  CUFFT_SAFE_CALL(cufftSetStream(**plan, **stream));
  CUFFT_SAFE_CALL(cufftExecZ2D(**plan, (cufftDoubleComplex *)**idata,
                               (cufftDoubleReal *)**odata));
#endif
}

extern "C" void FC_FUNC_(cuda_fft_execute_z2z_forward,
                         CUDA_FFT_EXECUTE_Z2Z_FORWARD)(cufftHandle **plan,
                                                       CUdeviceptr **idata,
                                                       CUdeviceptr **odata,
                                                       cudaStream_t **stream) {
#ifdef HAVE_CUDA
  CUFFT_SAFE_CALL(cufftSetStream(**plan, **stream));
  CUFFT_SAFE_CALL(cufftExecZ2Z(**plan, (cufftDoubleComplex *)**idata,
                               (cufftDoubleComplex *)**odata, CUFFT_FORWARD));
#endif
}

extern "C" void FC_FUNC_(cuda_fft_execute_z2z_backward,
                         CUDA_FFT_EXECUTE_Z2Z_BACKWARD)(cufftHandle **plan,
                                                        CUdeviceptr **idata,
                                                        CUdeviceptr **odata,
                                                        cudaStream_t **stream) {
#ifdef HAVE_CUDA
  CUFFT_SAFE_CALL(cufftSetStream(**plan, **stream));
  CUFFT_SAFE_CALL(cufftExecZ2Z(**plan, (cufftDoubleComplex *)**idata,
                               (cufftDoubleComplex *)**odata, CUFFT_INVERSE));
#endif
}
