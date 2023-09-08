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

__kernel void X(projector_commutator_bra)(const int nmat,
               __global int const * restrict offsets,
               __global rtype const * restrict matrix,
               __global int const * restrict map,
               __global double const * restrict scal,
               __global double const * restrict position,
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
  rtype aa1 = 0.0;
  rtype aa2 = 0.0;
  rtype aa3 = 0.0;
  
  for(int ip = 0; ip < npoints; ip++){
    aa0 += MUL(CONJ(matrix[matrix_offset + ip + nppj]),psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa1 += position[(map_offset + ip)*3 + 0]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa2 += position[(map_offset + ip)*3 + 1]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa3 += position[(map_offset + ip)*3 + 2]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
  }

  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0] = scal[scal_offset + ipj]*aa0;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1] = scal[scal_offset + ipj]*aa1;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2] = scal[scal_offset + ipj]*aa2;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3] = scal[scal_offset + ipj]*aa3;
  
}


__kernel void X(projector_commutator_bra_phase)(const int nmat,
               __global int const * restrict offsets,
               __global rtype const * restrict matrix,
               __global int const * restrict map,
               __global double const * restrict scal,
               __global double const * restrict position,
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
  double2 aa1 = 0.0;
  double2 aa2 = 0.0;
  double2 aa3 = 0.0;
  
  for(int ip = 0; ip < npoints; ip++){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa0 += MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
    aa1 += position[(map_offset + ip)*3 + 0]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
    aa2 += position[(map_offset + ip)*3 + 1]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
    aa3 += position[(map_offset + ip)*3 + 2]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
  }

  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0] = scal[scal_offset + ipj]*aa0;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1] = scal[scal_offset + ipj]*aa1;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2] = scal[scal_offset + ipj]*aa2;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3] = scal[scal_offset + ipj]*aa3;
  
}


__kernel void X(projector_commutator_bra_phase_spiral)(const int nmat,
               __global int const * restrict offsets,
               __global rtype const * restrict matrix,
               __global int const * restrict map,
               __global double const * restrict scal,
               __global double const * restrict position,
               __global double2 const * restrict psi, const int ldpsi,
               __global double2 * restrict projection, const int ldprojection,
               __global double2 const * restrict phases, const int phases_offset,
               __global int const * restrict spin_to_phase,
	       const int nphase
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
  double2 aa1 = 0.0;
  double2 aa2 = 0.0;
  double2 aa3 = 0.0;

  const int phases_total_offset = phases_offset + map_offset*nphase + npoints * spin_to_phase[ist];

  for(int ip = 0; ip < npoints; ip++){
    double2 phasepsi = complex_mul(phases[phases_total_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa0 += MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
    aa1 += position[(map_offset + ip)*3 + 0]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
    aa2 += position[(map_offset + ip)*3 + 1]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
    aa3 += position[(map_offset + ip)*3 + 2]*MUL(CONJ(matrix[matrix_offset + ip + nppj]),phasepsi);
  }

  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0] = scal[scal_offset + ipj]*aa0;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1] = scal[scal_offset + ipj]*aa1;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2] = scal[scal_offset + ipj]*aa2;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3] = scal[scal_offset + ipj]*aa3;
  
}

__kernel void X(projector_commutator_ket)(const int nmat,
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
  rtype aa1 = 0.0;
  rtype aa2 = 0.0;
  rtype aa3 = 0.0;
    
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa0 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0]);
    aa1 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1]);
    aa2 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2]);
    aa3 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3]);
  }

  cpsi1[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 0]*aa0 - aa1;
  cpsi2[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 1]*aa0 - aa2;
  cpsi3[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 2]*aa0 - aa3;

}


__kernel void X(projector_commutator_ket_phase)(const int nmat,
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
  double2 aa1 = 0.0;
  double2 aa2 = 0.0;
  double2 aa3 = 0.0;
    
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa0 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0]);
    aa1 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1]);
    aa2 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2]);
    aa3 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3]);
  }

  double2 phase = complex_conj(phases[phases_offset + map_offset + ip]);
  aa0 = complex_mul(phase, aa0);
  aa1 = complex_mul(phase, aa1);
  aa2 = complex_mul(phase, aa2);
  aa3 = complex_mul(phase, aa3);
  
  cpsi1[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 0]*aa0 - aa1;
  cpsi2[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 1]*aa0 - aa2;
  cpsi3[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 2]*aa0 - aa3;

}


__kernel void X(projector_commutator_ket_phase_spiral)(const int nmat,
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

  double2 aa0 = 0.0;
  double2 aa1 = 0.0;
  double2 aa2 = 0.0;
  double2 aa3 = 0.0;
    
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa0 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0]);
    aa1 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1]);
    aa2 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2]);
    aa3 += MUL(matrix[matrix_offset + ip + npoints*ipj],projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3]);
  }

  const int phases_total_offset = phases_offset + map_offset*nphase + npoints * spin_to_phase[ist];

  double2 phase = complex_conj(phases[phases_total_offset  + ip]);
  aa0 = complex_mul(phase, aa0);
  aa1 = complex_mul(phase, aa1);
  aa2 = complex_mul(phase, aa2);
  aa3 = complex_mul(phase, aa3);
  
  cpsi1[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 0]*aa0 - aa1;
  cpsi2[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 1]*aa0 - aa2;
  cpsi3[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 2]*aa0 - aa3;

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

  mixprojection[4*(ist + ((scal_offset + ipj)<<ldprojection))+0] = aa;
  mixprojection[4*(ist + ((scal_offset + ipj)<<ldprojection))+1] = bb;
  mixprojection[4*(ist + ((scal_offset + ipj)<<ldprojection))+2] = cc;
  mixprojection[4*(ist + ((scal_offset + ipj)<<ldprojection))+3] = dd;

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
  double2 bb0 = 0.0, bb1 = 0.0, bb2 = 0.0, bb3 = 0.0;
  double2 cc0 = 0.0, cc1 = 0.0, cc2 = 0.0, cc3 = 0.0;
  double2 dd0 = 0.0, dd1 = 0.0, dd2 = 0.0, dd3 = 0.0;
  if (mix_offset_0 != -1) {
    for(int jpj = 0; jpj < nprojs; jpj++){
      aa0 += complex_mul(mix[mix_offset_0 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+0]);
      aa1 += complex_mul(mix[mix_offset_1 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+0]);
      aa2 += complex_mul(mix[mix_offset_2 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+0]);
      aa3 += complex_mul(mix[mix_offset_3 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+0]);
      bb0 += complex_mul(mix[mix_offset_0 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+1]);
      bb1 += complex_mul(mix[mix_offset_1 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+1]);
      bb2 += complex_mul(mix[mix_offset_2 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+1]);
      bb3 += complex_mul(mix[mix_offset_3 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+1]);
      cc0 += complex_mul(mix[mix_offset_0 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+2]);
      cc1 += complex_mul(mix[mix_offset_1 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+2]);
      cc2 += complex_mul(mix[mix_offset_2 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+2]);
      cc3 += complex_mul(mix[mix_offset_3 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+2]);
      dd0 += complex_mul(mix[mix_offset_0 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+3]);
      dd1 += complex_mul(mix[mix_offset_1 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+3]);
      dd2 += complex_mul(mix[mix_offset_2 + nprojs*jpj + ipj], projection[4*(ist*2+1 + ((scal_offset + jpj)<<ldprojection))+3]);
      dd3 += complex_mul(mix[mix_offset_3 + nprojs*jpj + ipj], projection[4*(ist*2   + ((scal_offset + jpj)<<ldprojection))+3]);
    }
  }
  else {
    aa0 = projection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+0];
    aa1 = projection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+0];
    bb0 = projection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+1];
    bb1 = projection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+1];
    cc0 = projection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+2];
    cc1 = projection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+2];
    dd0 = projection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+3];
    dd1 = projection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+3];

  }

  mixprojection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+0] = aa0 + aa2;
  mixprojection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+0] = aa1 + aa3;
  mixprojection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+1] = bb0 + bb2;
  mixprojection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+1] = bb1 + bb3;
  mixprojection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+2] = cc0 + cc2;
  mixprojection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+2] = cc1 + cc3;
  mixprojection[4*(ist*2   + ((scal_offset + ipj)<<ldprojection))+3] = dd0 + dd2;
  mixprojection[4*(ist*2+1 + ((scal_offset + ipj)<<ldprojection))+3] = dd1 + dd3;
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
