/*
 Copyright (C) 2011 X. Andrade

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

__kernel void get_points(const int sp,
			 const int ep,
			 __global const int * linear_to_ist,
			 __global const int * linear_to_idim,
			 const int nst,
			 __global double const * restrict psi, const int ldpsi,
			 __global double * restrict points,
                         const int ldpoints1, const int ldpoints2){
  const int ist_linear = get_global_id(0);
  const int ip  = get_global_id(1);
  const int sip = ip + sp - 1;

  if(ist_linear < nst) {
    int ist = linear_to_ist[ist_linear];
    int idim = linear_to_idim[ist_linear];

    points[ldpoints1*ldpoints2*ip + idim*ldpoints1 + ist] = psi[ldpsi*sip + ist_linear];
  }

}

__kernel void set_points(const int sp,
			 const int ep,
			 __global const int * linear_to_ist,
			 __global const int * linear_to_idim,
			 const int nst,
			 __global double const * restrict points,
                         const int ldpoints1, const int ldpoints2,
			 __global double * restrict psi, const int ldpsi){
  const int ist_linear = get_global_id(0);
  const int ip  = get_global_id(1);
  const int sip = ip + sp - 1;

  if(ist_linear < nst) {
    int ist = linear_to_ist[ist_linear];
    int idim = linear_to_idim[ist_linear];

    psi[ldpsi*sip + ist_linear] = points[ldpoints1*ldpoints2*ip + idim*ldpoints1 + ist];
  }

}
/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
