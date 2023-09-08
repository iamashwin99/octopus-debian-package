/*
 Copyright (C) 2010 X. Andrade

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
#include <cl_rtype.h>

#define OFFSET_SIZE 6 /* defined in src/hamiltonian/hamiltonian_base.F90 */


__kernel void X(projector_bra)(const int nmat,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global double const * restrict scal,
          __global rtype const * restrict psi, const int ldpsi,
          __global rtype * restrict projection, const int ldprojection
          ){
  
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size = 1;
#endif

  const int ist  = get_global_id(0) / my_warp_size;
  const int ipj  = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints * ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size == 0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * (get_local_id(0)%my_warp_size) ;
  const int end   = min(start + slice, npoints);
  const int step  = 1;
#else
  const int start = 0;
  const int end   = npoints;
  const int step  = 1;
#endif

  rtype aa = 0.0;
  for(int ip = start; ip < end; ip += step){
    aa += MUL(CONJ(matrix[matrix_offset + ip + nppj]), psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
  }

#ifdef CUDA
  aa = X(warpReduce)(aa);
  if(get_local_id(0)%my_warp_size == 0)
#endif
    projection[ist + ((scal_offset + ipj)<<ldprojection)] = MUL(scal[scal_offset + ipj], aa);

}

__kernel void X(projector_bra_phase)(const int nmat,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global double const * restrict scal,
          __global double2 const * restrict psi, const int ldpsi,
          __global double2 * restrict projection, const int ldprojection,
          __global double2 const * restrict phases, const int phases_offset
          ){
  

#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size;  // the kernel is to be called for (at least) all ist<nst_linear.
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size==0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice, npoints );
  const int step  = 1;
#else
  const int start = 0;
  const int end = npoints;
  const int step = 1;
#endif

  double2 aa = 0.0;
  for(int ip = start; ip < end; ip+=step){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa += MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
  }

#ifdef CUDA
  aa = zwarpReduce(aa);
  if(get_local_id(0)%my_warp_size==0) 
#endif
    projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa;

}

// The *_phase_spiral kernels also could be used for the *_phase case, if a spin_to_phase array with all zeros is passed.
// This could make the old kernel obsolete, at the expense of slightly more obscure code.

__kernel void X(projector_bra_phase_spiral)(const int nmat,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global double const * restrict scal,
          __global double2 const * restrict psi, const int ldpsi,
          __global double2 * restrict projection, const int ldprojection,
          __global double2 const * restrict phases, const int phases_offset, 
          __global int const * restrict spin_to_phase,
	  const int nphase
          ){
  

#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size;
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0]; // number of points in projector_matrix(imat)
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1]; // number of projectors
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2]; // cumulative of pmat%npoints * pmat%nprojs
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3]; // cumulative of pmat%npoints for each imap
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4]; // cumulative of pmat%nprojs

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size==0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice, npoints );
  const int step  = 1;
#else
  const int start = 0;
  const int end = npoints;
  const int step = 1;
#endif

  // phase_offset = (psib%ik - std%kpt%start)*this%total_points)
  // phases[phase_offset] points to the starting address of projector_phases(1,1,1,ik)
  // phases[phase_offset + map_offset + ip] points to projector_phases(ip,1,imat,ik) (for nphase = 1!!)

  // We need to access phases[] at: projector_phases(ip,iphase,imat,ik)
  // The starting addresses for successive iphase blocks are shifted by npoints (which is imap dependent).
  // The correct iphase can be obtained from spin_to_phase[ist], where ist ranged from 1:nst_linear

  const int phases_total_offset = phases_offset + map_offset*nphase + npoints * spin_to_phase[ist];

  double2 aa = 0.0;
  for(int ip = start; ip < end; ip+=step){
    double2 phasepsi = complex_mul(phases[phases_total_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa += MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
  }

#ifdef CUDA
  aa = zwarpReduce(aa);
  if(get_local_id(0)%my_warp_size==0) 
#endif
    projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa;

}

__kernel void X(projector_ket)(const int nmat,
          const int imat_offset,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global rtype const * restrict projection, const int ldprojection,
          __global rtype * restrict psi, const int ldpsi,
          const int ip_start,
          const int ip_end
          ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1);
  const int imat = get_global_id(2) + imat_offset;

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  // these need to be determined at runtime from imat and i_self_region
  if( (ip < ip_start) || (ip >= ip_end) ) return;

  rtype aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += MUL(matrix[matrix_offset + ip + npoints*ipj], projection[ist + ((scal_offset + ipj)<<ldprojection)]);
  }

  psi[((map[map_offset + ip] - 1)<<ldpsi) + ist] += aa;

}

__kernel void X(projector_ket_phase)(const int nmat,
          const int imat_offset,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global double2 const * restrict projection, const int ldprojection,
          __global double2 * restrict psi, const int ldpsi,
          const int ip_start,
          const int ip_end,
          __global double2 const * restrict phases, const int phases_offset
          ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1);
  const int imat = get_global_id(2) + imat_offset;

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ip < ip_start || ip >= ip_end) return;

  double2 aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[ist + ((scal_offset + ipj)<<ldprojection)]);
  }

  psi[((map[map_offset + ip] - 1)<<ldpsi) + ist] += complex_mul(complex_conj(phases[phases_offset + map_offset + ip]), aa);

}


__kernel void X(projector_ket_phase_spiral)(const int nmat,
          const int imat_offset,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global double2 const * restrict projection, const int ldprojection,
          __global double2 * restrict psi, const int ldpsi,
          const int ip_start,
          const int ip_end,
          __global double2 const * restrict phases, const int phases_offset,
          __global int const * restrict spin_to_phase,
          const int nphase
          ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1);
  const int imat = get_global_id(2) + imat_offset;

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ip < ip_start || ip >= ip_end) return;

  double2 aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[ist + ((scal_offset + ipj)<<ldprojection)]);
  }

  const int phases_total_offset = phases_offset + map_offset*nphase + npoints * spin_to_phase[ist];

  psi[((map[map_offset + ip] - 1)<<ldpsi) + ist] += complex_mul(complex_conj(phases[phases_total_offset + ip]), aa);

}

__kernel void dprojector_mix(const int nmat,
          __global int const * restrict offsets,
          __global double const * restrict mix,
          __global double * const restrict projection, const int ldprojection,
          __global double * restrict mixprojection
          ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];
  const int mix_offset    = offsets[OFFSET_SIZE*imat + 5];

  if(ipj >= nprojs) return;
  
  double aa = 0.0;
  if(mix_offset != -1) {
    for(int jpj = 0; jpj < nprojs; jpj++){
      aa += mix[mix_offset + nprojs*jpj + ipj] * projection[ist + ((scal_offset + jpj)<<ldprojection)];
    }
  }
  else {
    aa = projection[ist + ((scal_offset + ipj)<<ldprojection)];
  }

  mixprojection[ist + ((scal_offset + ipj)<<ldprojection)] = aa;

}

__kernel void zprojector_mix(const int nmat,
          __global int const * restrict offsets,
          __global double2 const * restrict mix,
          __global double2 * const restrict projection, const int ldprojection,
          __global double2 * restrict mixprojection
          ){
  
  const int ist  = get_global_id(0);
  const int ipj  = get_global_id(1);
  const int imat = get_global_id(2);

  const int nprojs       = offsets[OFFSET_SIZE*imat + 1];
  const int scal_offset  = offsets[OFFSET_SIZE*imat + 4];
  const int mix_offset_0 = offsets[OFFSET_SIZE*imat + 5];
  const int mix_offset_1 = mix_offset_0 + nprojs*nprojs;
  const int mix_offset_2 = mix_offset_1 + nprojs*nprojs;
  const int mix_offset_3 = mix_offset_2 + nprojs*nprojs;

  if(ipj >= nprojs) return;

  double2 aa0 = 0.0, aa1 = 0.0, aa2 = 0.0, aa3 = 0.0;
  if (mix_offset_0 != -1) {
    for(int jpj = 0; jpj < nprojs; jpj++){
      aa0 += complex_mul(mix[mix_offset_0 + nprojs*jpj + ipj], projection[ist*2   + ((scal_offset + jpj)<<ldprojection)]);
      aa1 += complex_mul(mix[mix_offset_1 + nprojs*jpj + ipj], projection[ist*2+1 + ((scal_offset + jpj)<<ldprojection)]);
      aa2 += complex_mul(mix[mix_offset_2 + nprojs*jpj + ipj], projection[ist*2+1 + ((scal_offset + jpj)<<ldprojection)]);
      aa3 += complex_mul(mix[mix_offset_3 + nprojs*jpj + ipj], projection[ist*2   + ((scal_offset + jpj)<<ldprojection)]);
    }
  }
  else {
    aa0 = projection[ist*2   + ((scal_offset + ipj)<<ldprojection)];
    aa1 = projection[ist*2+1 + ((scal_offset + ipj)<<ldprojection)];
  }

  mixprojection[ist*2   + ((scal_offset + ipj)<<ldprojection)] = aa0 + aa2;
  mixprojection[ist*2+1 + ((scal_offset + ipj)<<ldprojection)] = aa1 + aa3;

}

__kernel void dprojector_mix_commutator(const int nmat,
          __global int const * restrict offsets,
          __global double const * restrict mix,
          __global double * const restrict projection, const int ldprojection,
          __global double * restrict mixprojection
          ){

  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];
  const int mix_offset    = offsets[OFFSET_SIZE*imat + 5];

  if(ipj >= nprojs) return;

  double aa = 0.0;
  double bb = 0.0;
  double cc = 0.0;
  double dd = 0.0;
  if(mix_offset != -1) {
    for(int jpj = 0; jpj < nprojs; jpj++){
      aa += mix[mix_offset + nprojs*jpj + ipj] * projection[4*(ist + ((scal_offset + jpj)<<ldprojection))+0];
      bb += mix[mix_offset + nprojs*jpj + ipj] * projection[4*(ist + ((scal_offset + jpj)<<ldprojection))+1];
      cc += mix[mix_offset + nprojs*jpj + ipj] * projection[4*(ist + ((scal_offset + jpj)<<ldprojection))+2];
      dd += mix[mix_offset + nprojs*jpj + ipj] * projection[4*(ist + ((scal_offset + jpj)<<ldprojection))+3];
    }
  }
  else {
    aa = projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0];
    bb = projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1];
    cc = projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2];
    dd = projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3];
  }

  mixprojection[ist + ((scal_offset + ipj)<<ldprojection)] = aa;
  mixprojection[ist + ((scal_offset + ipj)<<ldprojection)] = bb;
  mixprojection[ist + ((scal_offset + ipj)<<ldprojection)] = cc;
  mixprojection[ist + ((scal_offset + ipj)<<ldprojection)] = dd;

}

__kernel void zprojector_mix_commutator(const int nmat,
          __global int const * restrict offsets,
          __global double2 const * restrict mix,
          __global double2 * const restrict projection, const int ldprojection,
          __global double2 * restrict mixprojection
          ){

  const int ist  = get_global_id(0);
  const int ipj  = get_global_id(1);
  const int imat = get_global_id(2);

  const int nprojs       = offsets[OFFSET_SIZE*imat + 1];
  const int scal_offset  = offsets[OFFSET_SIZE*imat + 4];
  const int mix_offset_0 = offsets[OFFSET_SIZE*imat + 5];
  const int mix_offset_1 = mix_offset_0 + nprojs*nprojs;
  const int mix_offset_2 = mix_offset_1 + nprojs*nprojs;
  const int mix_offset_3 = mix_offset_2 + nprojs*nprojs;

  if(ipj >= nprojs) return;

  double2 aa0 = 0.0, aa1 = 0.0, aa2 = 0.0, aa3 = 0.0;
  if (mix_offset_0 != -1) {
    for(int jpj = 0; jpj < nprojs; jpj++){
      aa0 += complex_mul(mix[mix_offset_0 + nprojs*jpj + ipj], projection[ist*2   + ((scal_offset + jpj)<<ldprojection)]);
      aa1 += complex_mul(mix[mix_offset_1 + nprojs*jpj + ipj], projection[ist*2+1 + ((scal_offset + jpj)<<ldprojection)]);
      aa2 += complex_mul(mix[mix_offset_2 + nprojs*jpj + ipj], projection[ist*2+1 + ((scal_offset + jpj)<<ldprojection)]);
      aa3 += complex_mul(mix[mix_offset_3 + nprojs*jpj + ipj], projection[ist*2   + ((scal_offset + jpj)<<ldprojection)]);
    }
  }
  else {
    aa0 = projection[ist*2   + ((scal_offset + ipj)<<ldprojection)];
    aa1 = projection[ist*2+1 + ((scal_offset + ipj)<<ldprojection)];
  }

  mixprojection[ist*2   + ((scal_offset + ipj)<<ldprojection)] = aa0 + aa2;
  mixprojection[ist*2+1 + ((scal_offset + ipj)<<ldprojection)] = aa1 + aa3;

}


__kernel void X(projector_bra_force)(const int idir,
            const int ndim, const int ratio, 
            const int nmat,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global rtype const * restrict psi, const int ldpsi,
          __global rtype * restrict projection, const int ldprojection
          ){
  
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size;
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size==0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * get_local_id(0) ;
  const int end   = min( start + slice , npoints );
#else
  const int start = 0;
  const int end = npoints;
#endif

  rtype aa = 0.0;
  for(int ip = start; ip < end; ip++){
    aa += MUL(CONJ(matrix[matrix_offset + ip + nppj]),psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
  }

  // See comment below
  const int ist_scaled = ((int) (ist/ratio))*ratio;

#ifdef CUDA
  aa = X(warpReduce)(aa);
  if(get_local_id(0) == 0) 
#endif
    // In the complex case, we double the number of states and use real-real multiplications.
    // The line below takes care of ordering properly the real values, such that reading of pack_size(1) 
    // complex numbers on the CPU side produces the expected result.
    // In this case, ratio would be 2, instead of 1 in the other cases
    // Note that ldprojection already contains the multiplication by ratio
    projection[idir*ratio + ist%ratio + ndim*(ist_scaled + ((scal_offset + ipj)<<ldprojection))] = aa;
}

__kernel void X(projector_bra_force_phase)(const int idir,
            const int ndim, const int ratio, 
            const int nmat,
          __global int const * restrict offsets,
          __global rtype const * restrict matrix,
          __global int const * restrict map,
          __global double2 const * restrict psi, const int ldpsi,
          __global double2 * restrict projection, const int ldprojection,
          __global double2 const * restrict phases, const int phases_offset
          ){
  

#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size;  // the kernel is to be called for (at least) all ist<nst_linear.
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size==0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( start + slice, npoints );
#else
  const int start = 0;
  const int end = npoints;
#endif

  double2 aa = 0.0;
  for(int ip = start; ip < end; ip++){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa += MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
  }

#ifdef CUDA
  aa = zwarpReduce(aa);
  if(get_local_id(0)%my_warp_size==0) 
#endif
    projection[idir + ndim*(ist + ((scal_offset + ipj)<<ldprojection))] = aa;
}

__kernel void X(projector_r_vnl_bra)(const int nmat,
               __global int const * restrict offsets,
               __global rtype const * restrict matrix,
               __global int const * restrict map,
               __global double const * restrict scal,
               __global rtype const * restrict psi, const int ldpsi,
               __global rtype * restrict projection, const int ldprojection
               ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

  rtype aa0 = 0.0;
  
  for(int ip = 0; ip < npoints; ip++){
    aa0 += MUL(CONJ(matrix[matrix_offset + ip + nppj]),psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
  }

  projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa0;
  
}


__kernel void X(projector_r_vnl_bra_phase)(const int nmat,
               __global int const * restrict offsets,
               __global rtype const * restrict matrix,
               __global int const * restrict map,
               __global double const * restrict scal,
               __global double2 const * restrict psi, const int ldpsi,
               __global double2 * restrict projection, const int ldprojection,
               __global double2 const * restrict phases, const int phases_offset
               ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

  double2 aa0 = 0.0;
  
  for(int ip = 0; ip < npoints; ip++){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa0 += MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
  }

  projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa0;
}

__kernel void X(projector_r_vnl_ket)(const int nmat,
               const int imat_offset,
               __global int const * restrict offsets,
               __global rtype const * restrict matrix,
               __global int const * restrict map,
               __global double const * restrict position,
               __global rtype const * restrict projection, const int ldprojection,
               __global rtype * restrict cpsi1,
               __global rtype * restrict cpsi2,
               __global rtype * restrict cpsi3, const int ldpsi,
               const int ip_start,
               const int ip_end
               ){

  const int ist = get_global_id(0);
  const int ip = get_global_id(1);
  const int imat = get_global_id(2) + imat_offset;

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ip < ip_start || ip >= ip_end) return;
  
  rtype aa0 = 0.0;
    
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa0 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0]);
  }

  cpsi1[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 0]*aa0;
  cpsi2[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 1]*aa0;
  cpsi3[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 2]*aa0;

}

__kernel void X(projector_r_vnl_ket_phase)(const int nmat,
               const int imat_offset,
               __global int const * restrict offsets,
               __global rtype const * restrict matrix,
               __global int const * restrict map,
               __global double const * restrict position,
               __global double2 const * restrict projection, const int ldprojection,
               __global double2 * restrict cpsi1,
               __global double2 * restrict cpsi2,
               __global double2 * restrict cpsi3, const int ldpsi,
               const int ip_start,
               const int ip_end,
               __global double2 const * restrict phases, const int phases_offset
               ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1);
  const int imat = get_global_id(2) + imat_offset;

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ip < ip_start || ip >= ip_end) return;

  double2 aa0 = 0.0;
    
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa0 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[ist + ((scal_offset + ipj)<<ldprojection)]);
  }

  double2 phase = complex_conj(phases[phases_offset + map_offset + ip]);
  aa0 = complex_mul(phase, aa0);
  
  cpsi1[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 0]*aa0;
  cpsi2[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 1]*aa0;
  cpsi3[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 2]*aa0;

}


/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
