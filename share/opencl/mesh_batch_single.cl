/*
 Copyright (C) 2012-2020 X. Andrade, M. Lueders
 Copyright (C) 2023 N. Tancogne-Dejean

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

#include <cl_global.h>
#include <cl_complex.h>
#include <cl_reduce.h>


// The X(batch_mf_dotp) kernels should be called on a global grid of (nst, 1, 1)
__kernel void zbatch_mf_dotp(
      const int np,    //< number of mesh points
      const int nst,   //< number of states
      const int ndim,  //< number of spin components (2 for spinors)
      __global double2 const * restrict xx_buffer, const int ldxx,   //< batch of states
      __global double2 const * restrict psi_buffer, const int ldpsi, //< single state
      __global double2 * restrict  dot_buffer                         //< vector of dot products
) {
  int ist = get_global_id(0);

  double2 tmp_dot = complex_number(0.0, 0.0);

  if(ist >= nst) return;

  for(int ip=0; ip<np; ip++) {
    for(int idim=0; idim< ndim; idim++) {
      tmp_dot += complex_dotp(xx_buffer[idim + (ndim)*ist + (ip<<ldxx)], psi_buffer[ip + (idim<<ldpsi)]);
    }
  }
  dot_buffer[ist] = tmp_dot;
}

__kernel void dbatch_mf_dotp(
      const int np,    //< number of mesh points
      const int nst,   //< number of states
      const int ndim,  //< number of spin components (2 for spinors)
      __global double const * restrict xx_buffer, const int ldxx,   //< batch of states
      __global double const * restrict psi_buffer, const int ldpsi, //< single state
      __global double * restrict  dot_buffer                         //< vector of dot products
) {
  int ist = get_global_id(0);

  double tmp_dot = 0.0;

  if(ist >= nst) return;

  for(int ip=0; ip<np; ip++) {
    for(int idim=0; idim< ndim; idim++) {
      tmp_dot += xx_buffer[idim + (ndim)*ist + (ip<<ldxx)] * psi_buffer[ip + (idim<<ldpsi)];
    }
  }
  dot_buffer[ist] = tmp_dot;
}



// See for instance 
// https://www.nvidia.com/content/GTC-2010/pdfs/2131_GTC2010.pdf
// https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
__kernel void dbatch_dotpv(
      const int np, const int nst_linear,
      const int ndim, const double vol_elem,
      __global double const * restrict aa_buffer, const int ldaa,  
      __global double const * restrict bb_buffer, const int ldbb, 
      __global double * restrict dot_buffer                    
#if !defined(CUDA)
      , __local double * cache
#endif
     ){

#if defined(CUDA)
  extern __shared__ double cache[];
#endif

  const unsigned int ip = get_global_id(0);
  const unsigned int ist_linear = get_global_id(1);
  const unsigned int ist = (int) (ist_linear/ndim);
  const unsigned int localIdx = get_local_id(0);
  const unsigned int tid = localIdx + get_local_id(1) * get_local_size(0);

  cache[tid] = 0.0;
  if (ip < np && ist_linear < nst_linear) {
    cache[tid] = aa_buffer[ist_linear + (ip<<ldaa)] * bb_buffer[ist_linear + (ip<<ldbb)];
  }

#ifdef CUDA
  int my_warp_size = warpSize;
  if(get_local_size(0)/2 <= warpSize) my_warp_size=0;
#else
  const int my_warp_size=0;
#endif

  barrier(CLK_LOCAL_MEM_FENCE);

  for (unsigned int i = get_local_size(0)/2; i>my_warp_size; i>>=1) {
    if (localIdx < i) {
      cache[tid] += cache[tid + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

#ifdef CUDA
  if (localIdx < warpSize && get_local_size(0)/2 > warpSize) {
    dwarpReduce_shared(cache, tid);
  }
#endif

  if (localIdx == 0 && ist_linear < nst_linear && ip < np) {
    atomicAdd(&dot_buffer[ist], vol_elem*cache[tid]);
  }

}

__kernel void zbatch_dotpv(
      const int np, const int nst_linear, 
      const int ndim, const double vol_elem,
      __global double2 const * restrict aa_buffer, const int ldaa,
      __global double2 const * restrict bb_buffer, const int ldbb,
      __global double2 * restrict  dot_buffer
#if !defined(CUDA)
      , __local double * cache
#endif
     ){

#if defined(CUDA)
  extern __shared__ double cache[];
#endif

  const unsigned int ip = get_global_id(0);
  const unsigned int ist_linear = get_global_id(1);
  const unsigned int ist = (int) (ist_linear/ndim);
  const unsigned int localIdx = get_local_id(0);
  const unsigned int tidx = localIdx + get_local_id(1) * get_local_size(0);
  const unsigned int tidy = tidx + get_local_size(0)*get_local_size(1);

  cache[tidx] = 0.0;
  cache[tidy] = 0.0;
  if (ip < np && ist_linear < nst_linear) {
    double2 tmp = complex_dotp(aa_buffer[ist_linear + (ip<<ldaa)], bb_buffer[ist_linear + (ip<<ldbb)]);
    cache[tidx] = tmp.x;
    cache[tidy] = tmp.y;
  }

#ifdef CUDA
  int my_warp_size = warpSize;
  if(get_local_size(0)/2 <= warpSize) my_warp_size=0;
#else
  const int my_warp_size=0;
#endif

  barrier(CLK_LOCAL_MEM_FENCE);

  for (unsigned int i = get_local_size(0)/2; i>my_warp_size; i>>=1) {
    if (localIdx < i) {
      cache[tidx] += cache[tidx + i];
      cache[tidy] += cache[tidy + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

#ifdef CUDA
  if (localIdx < warpSize && get_local_size(0)/2 > warpSize) {
    dwarpReduce_shared(cache, tidx);
    dwarpReduce_shared(cache, tidy);
  }
#endif

  if (localIdx == 0 && ist_linear < nst_linear && ip < np) {
#ifdef CUDA
    atomicAdd(&(dot_buffer[ist].x), vol_elem*cache[tidx]);
    atomicAdd(&(dot_buffer[ist].y), vol_elem*cache[tidy]);
#else
    atomicAdd(&(dot_buffer[ist]), vol_elem*cache[tidx]);
    atomicAdd(&(dot_buffer[ist])+1, vol_elem*cache[tidy]);
#endif
  }
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
