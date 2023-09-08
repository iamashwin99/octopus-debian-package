/*
 Copyright (C) 2022 N. Tancogne-Dejean

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

//------------------------------------------------------------------------------------
__kernel void dftu_projector_bra(const int np, const int nst, const int norbs,
      __global const double* __restrict orbs, const int ldorbs,
      __global const double* __restrict psi, const int ldpsi,
      __global double* __restrict dot_buffer 
) {
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  const int iorb = get_global_id(1);
  if(ist >= nst) return;
  if(iorb >= norbs) return;

  const unsigned int tid = get_local_id(0)%my_warp_size;

#ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * tid ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  double tmp_dot = 0.0;
  for(int ip=start; ip<end; ip++) {
    tmp_dot += psi[ist + (ip<<ldpsi)] * orbs[ip + (iorb<<ldorbs)];
  }
#ifdef CUDA
  tmp_dot = dwarpReduce(tmp_dot);
  if(tid == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = tmp_dot;
}

//------------------------------------------------------------------------------------
__kernel void dftu_projector_bra_submesh(const int np, const int nst, const int norbs,
      __global const double* __restrict orbs, const int ldorbs,
      __global const double* __restrict psi, const int ldpsi,
      __global double* __restrict dot_buffer,
      __global int const * restrict map
      ) {
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const unsigned int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  const unsigned int iorb = get_global_id(1);
  if(ist >= nst) return;
  if(iorb >= norbs) return;

  const unsigned int tid = get_local_id(0)%my_warp_size;

  //We slice the loop over np by slices of size warp_size
#ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * tid ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  double sum = 0.0;
  //do reduction for locally
  for(unsigned int ip=start; ip<end; ip++) {
    sum += psi[ist + ((map[ip]-1)<<ldpsi)] * orbs[ip + (iorb<<ldorbs)];
  }

  //Do reduction over one warp of threads and write the result to global memory
#ifdef CUDA
  sum = dwarpReduce(sum);
  if(tid == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = sum;
}

//------------------------------------------------------------------------------------
__kernel void dftu_projector_bra_cmplx(const int np, const int nst, const int norbs,
      __global const double2* __restrict orbs, const int ldorbs,
      __global const double2* __restrict psi, const int ldpsi,
      __global double2* __restrict dot_buffer
) {

#ifdef CUDA
  const int ist = get_global_id(0)/warpSize; // [0:nst-1]
  const int iorb = get_global_id(1);
  if(ist >= nst) return;
  if(iorb >= norbs) return;

  const unsigned int tid = get_local_id(0);

  //We slice the loop over np by slices of size warpSize
  const unsigned int slice = np%warpSize==0 ? np/warpSize : np/warpSize+1;
  const unsigned int start = slice * tid;
  const unsigned int end   = min( start + slice , np );

  double2 sum = 0.0;
  for(unsigned int ip = start; ip < end; ip++) {
    sum += complex_dotp(orbs[ip + (iorb<<ldorbs)], psi[ist + (ip<<ldpsi)]);
  }

  //Do reduction over one warp of threads and write the result to global memory
  sum = zwarpReduce(sum);
  if(tid == 0)
    dot_buffer[ist + (iorb<<ldpsi)] = sum;

#else
  const int ist = get_global_id(0);
  const int iorb = get_global_id(1);
  if(ist >= nst) return;
  if(iorb >= norbs) return;

  double2 sum = 0.0;
  for(int ip = 0; ip < np; ip++) {
    sum += complex_mul(psi[ist + (ip<<ldpsi)], complex_conj(orbs[ip + (iorb<<ldorbs)]));
  }
  dot_buffer[ist + (iorb<<ldpsi)] = sum;
#endif
}

//------------------------------------------------------------------------------------
__kernel void dftu_projector_bra_cmplx_submesh(const int np, const int nst, const int norbs,
      __global const double2* __restrict orbs, const int ldorbs,
      __global const double2* __restrict psi, const int ldpsi,
      __global double2* __restrict dot_buffer,
      __global int const * restrict map
      ) {
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  if(ist >= nst) return;

  const int iorb = get_global_id(1);
  if(iorb >= norbs) return;

  double2 tmp_dot = 0.0;

#ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  for(int ip=start; ip<end; ip++) {
    tmp_dot += complex_dotp(orbs[ip + (iorb<<ldorbs)], psi[ist + ((map[ip]-1)<<ldpsi)]);
  }

#ifdef CUDA
  tmp_dot = zwarpReduce(tmp_dot);
  if(get_local_id(0) == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = tmp_dot;

}

//------------------------------------------------------------------------------------
__kernel void dftu_projector_ket_submesh(const int norbs, const int ip_start, const int ip_end,
          __global double const * restrict weights,
          __global double const * restrict orbs, const int ldorbs,
          __global double * restrict psi, const int ldpsi,
          __global int const * restrict map
          ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if( (ip < ip_start) || (ip >= ip_end) ) return;

  double aa = 0.0;
  for(int iorb = 0; iorb < norbs; iorb++){
    aa += weights[ist + (iorb<<ldpsi)] * orbs[ip + (iorb<<ldorbs)];
  }

  psi[ist + ((map[ip] - 1)<<ldpsi)] += aa;

}

//------------------------------------------------------------------------------------
__kernel void dftu_projector_ket(const int norbs, const int npoints,
          __global double const * restrict weights,
          __global double const * restrict orbs, const int ldorbs,
          __global double * restrict psi, const int ldpsi
          ){

  const int ist = get_global_id(0);
  const int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if(ip >= npoints) return;

  double aa = 0.0;
  for(int iorb = 0; iorb < norbs; iorb++){
    aa += weights[ist + (iorb<<ldpsi)] * orbs[ip + (iorb<<ldorbs)];
  }

  psi[ist + (ip<<ldpsi)] += aa;

}

//------------------------------------------------------------------------------------
__kernel void dftu_projector_ket_cmplx_submesh(const int norbs, const int ip_start, const int ip_end,
          __global double2 const * restrict weights,
          __global double2 const * restrict orbs, const int ldorbs,
          __global double2 * restrict psi, const int ldpsi,
          __global int const * restrict map
          ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if( (ip < ip_start) || (ip >= ip_end) ) return;

  double2 aa = 0.0;
  for(int iorb = 0; iorb < norbs; iorb++){
    aa += complex_mul(weights[ist + (iorb<<ldpsi)], orbs[ip + (iorb<<ldorbs)]);
  }

  psi[ist + ((map[ip] - 1)<<ldpsi)] += aa;

}

//------------------------------------------------------------------------------------
__kernel void dftu_projector_ket_cmplx(const int norbs, const int npoints,
          __global double2 const * restrict weights,
          __global double2 const * restrict orbs, const int ldorbs,
          __global double2 * restrict psi, const int ldpsi
          ){

  const int ist = get_global_id(0);
  const int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if(ip >= npoints) return;

  double2 aa = 0.0;
  for(int iorb = 0; iorb < norbs; iorb++){
    aa += complex_mul(weights[ist + (iorb<<ldpsi)], orbs[ip + (iorb<<ldorbs)]);
  }

  psi[ist + (ip<<ldpsi)] += aa;

}

//------------------------------------------------------------------------------------
__kernel void dftu_pos_mat_elem(const int np, const int nst, const int norbs,
      __global const double* __restrict orbs, const int ldorbs,
      __global const double* __restrict psi, const int ldpsi,
      __global const double* __restrict xx,
      __global int const * restrict map,
      __global double* __restrict dot_buffer
) {
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  if(ist >= nst) return;

  const int iorb = get_global_id(1);
  if(iorb >= norbs) return;

  double tmp_dot = 0.0;

#ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  for(int ip=start; ip<end; ip++) {
    tmp_dot += psi[ist + ((map[ip]-1)<<ldpsi)] * xx[ip] * orbs[(map[ip]-1) + (iorb<<ldorbs)];
  }
#ifdef CUDA
  tmp_dot = dwarpReduce(tmp_dot);
  if(get_local_id(0)%my_warp_size == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = tmp_dot;
}

//------------------------------------------------------------------------------------
__kernel void dftu_pos_mat_elem_submesh(const int np, const int nst, const int norbs,
      __global const double* __restrict orbs, const int ldorbs,
      __global const double* __restrict psi, const int ldpsi,
      __global const double* __restrict xx,
      __global int const * restrict map,
      __global double* __restrict dot_buffer
      ) {
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  if(ist >= nst) return;

  const int iorb = get_global_id(1);
  if(iorb >= norbs) return;

  double tmp_dot = 0.0;

#ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  for(int ip=start; ip<end; ip++) {
    tmp_dot += psi[ist + ((map[ip]-1)<<ldpsi)] * xx[ip] * orbs[ip + (iorb<<ldorbs)];
  }

#ifdef CUDA
  tmp_dot = dwarpReduce(tmp_dot);
  if(get_local_id(0) == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = tmp_dot;
}

//------------------------------------------------------------------------------------
__kernel void dftu_pos_mat_elem_cmplx(const int np, const int nst, const int norbs,
      __global const double2* __restrict orbs, const int ldorbs,
      __global const double2* __restrict psi, const int ldpsi,
      __global const double* __restrict xx,
      __global int const * restrict map,
      __global double2* __restrict dot_buffer
) {

#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  if(ist >= nst) return;

  const int iorb = get_global_id(1);
  if(iorb >= norbs) return;

  double2 tmp_dot = 0.0;

#ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  for(int ip=start; ip<end; ip++) {
    tmp_dot += complex_dotp(orbs[(map[ip]-1) + (iorb<<ldorbs)], complex_mul(psi[ist + ((map[ip]-1)<<ldpsi)], xx[ip]));
  }

#ifdef CUDA
  tmp_dot = zwarpReduce(tmp_dot);
  if(get_local_id(0) == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = tmp_dot;
}

//------------------------------------------------------------------------------------
__kernel void dftu_pos_mat_elem_phase(const int np, const int nst, const int norbs,
      __global const double2* __restrict orbs, const int ldorbs,
      __global const double2* __restrict psi, const int ldpsi,
      __global const double* __restrict xx,
      __global int const * restrict map,
      __global double2* __restrict dot_buffer,
      __global const double2* __restrict phase
) {

#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  if(ist >= nst) return;

  const int iorb = get_global_id(1);
  if(iorb >= norbs) return;

  double2 tmp_dot = 0.0;

#ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  for(int ip=start; ip<end; ip++) {
    double2 phaseorb = complex_mul(orbs[(map[ip]-1) + (iorb<<ldorbs)], phase[ip]);
    tmp_dot += complex_dotp(phaseorb, complex_mul(psi[ist + ((map[ip]-1)<<ldpsi)], xx[ip]));
  }

#ifdef CUDA
  tmp_dot = zwarpReduce(tmp_dot);
  if(get_local_id(0) == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = tmp_dot;
}


//------------------------------------------------------------------------------------
__kernel void dftu_pos_mat_elem_cmplx_submesh(const int np, const int nst, const int norbs,
      __global const double2* __restrict orbs, const int ldorbs,
      __global const double2* __restrict psi, const int ldpsi,
      __global const double* __restrict xx,
      __global int const * restrict map,
      __global double2* __restrict dot_buffer
      ) {
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size; // [0:nst-1]
  if(ist >= nst) return;

  const int iorb = get_global_id(1);
  if(iorb >= norbs) return;

  double2 tmp_dot = 0.0;

  #ifdef CUDA
  const int slice = np%my_warp_size==0 ? np/my_warp_size : np/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice , np );
#else
  const int start = 0;
  const int end = np;
#endif

  for(int ip=start; ip<end; ip++) {
    tmp_dot += complex_dotp(orbs[ip + (iorb<<ldorbs)], complex_mul(psi[ist + ((map[ip]-1)<<ldpsi)], xx[ip]));
  }

#ifdef CUDA
  tmp_dot = zwarpReduce(tmp_dot);
  if(get_local_id(0) == 0)
#endif
    dot_buffer[ist + (iorb<<ldpsi)] = tmp_dot;

}



/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

