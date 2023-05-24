/*
 Copyright (C) 2010 X. Andrade
 Copyright (C) 2022 S. Ohlmann

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

__kernel void vpsi(const int offset,
		   const int np,
		   __global double const * restrict vv,
		   __global double const * restrict psi, const int ldpsi,
		   __global double * restrict vpsi, const int ldvpsi){

  const int ist = get_global_id(0);
  const int ip  = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if(ip < np){
    vpsi[(ip<<ldvpsi) + ist] += vv[offset + ip]*psi[(ip<<ldpsi) + ist];
  }

}

__kernel void vpsi_complex(const int offset,
		           const int np,
		           __global double const * restrict vv,
		           __global double const * restrict imvv,
		           __global double2 const * restrict psi, const int ldpsi,
		           __global double2 * restrict vpsi, const int ldvpsi){

  const int ist = get_global_id(0);
  const int ip  = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if(ip < np){
    vpsi[(ip<<ldvpsi) + ist] += complex_mul(complex_number(vv[offset + ip], imvv[offset + ip]), psi[(ip<<ldpsi) + ist]);
  }

}

__kernel void vpsi_spinors(const int np,
			   __global double const * restrict vv, const int ldvv,
			   __global double2 const * restrict psi, const int ldpsi,
			   __global double2 * restrict vpsi, const int ldvpsi){
  const int ist = 2*get_global_id(0);
  const int ip  = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if(ip < np){

    const double vi1 = vv[         ip];
    const double vi2 = vv[ldvv   + ip];
    const double vi3 = vv[2*ldvv + ip];
    const double vi4 = vv[3*ldvv + ip];

    const double2 psi1 = psi[ip*ldpsi + ist];
    const double2 psi2 = psi[ip*ldpsi + ist + 1];

#ifdef CUDA
    vpsi[ip*ldvpsi + ist] += vi1*psi1 + double2(vi3*psi2.x - vi4*psi2.y, vi3*psi2.y + vi4*psi2.x);
    vpsi[ip*ldvpsi + ist + 1] += vi2*psi2 + double2(vi3*psi1.x + vi4*psi1.y, vi3*psi1.y - vi4*psi1.x);
#else
    vpsi[ip*ldvpsi + ist] += vi1*psi1 + (double2)(vi3*psi2.x - vi4*psi2.y, vi3*psi2.y + vi4*psi2.x);
    vpsi[ip*ldvpsi + ist + 1] += vi2*psi2 + (double2)(vi3*psi1.x + vi4*psi1.y, vi3*psi1.y - vi4*psi1.x);
#endif

  }
}

__kernel void vpsi_spinors_complex(const int np,
			           __global double const * restrict vv, const int ldvv,
			           __global double const * restrict imvv, const int ldimvv,
			           __global double2 const * restrict psi, const int ldpsi,
			           __global double2 * restrict vpsi, const int ldvpsi){
  const int ist = 2*get_global_id(0);
  const int ip  = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if(ip < np){

    const double2 vi1 = complex_number(vv[         ip], imvv[           ip]);
    const double2 vi2 = complex_number(vv[ldvv   + ip], imvv[ldimvv   + ip]);
    const double2 vi3 = complex_number(vv[2*ldvv + ip], imvv[2*ldimvv + ip]);
    const double2 vi4 = complex_number(vv[3*ldvv + ip], imvv[3*ldimvv + ip]);
    const double2 imag = complex_number(0, 1);

    const double2 psi1 = psi[ip*ldpsi + ist];
    const double2 psi2 = psi[ip*ldpsi + ist + 1];

    vpsi[ip*ldvpsi + ist] += complex_mul(vi1, psi1) + complex_mul(vi3 + complex_mul(imag, vi4), psi2);
    vpsi[ip*ldvpsi + ist + 1] += complex_mul(vi2, psi2) + complex_mul(vi3 - complex_mul(imag, vi4), psi1);
  }
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

