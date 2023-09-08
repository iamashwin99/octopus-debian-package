/*
 Copyright (C) 2022,2023 S. Ohlmann

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

__kernel void pml_apply(const int with_medium,
    const int npoints,
    const int pml_dir,
    const double P_c,
    const int rs_sign,
    __global int * restrict pml_map,
    __global double2 * restrict grad,
    const int ldgrad,
    const __global double * restrict pml_a,
    const __global double * restrict pml_b,
    const __global double * restrict pml_c,
    __global double2 * restrict pml_conv_plus,
    __global double2 * restrict pml_conv_minus)
{
  const int field_dirs[3][2] = {{2, 3}, {1, 3}, {1, 2}};
  const int ip_in = get_global_id(0);
  if (ip_in >= npoints) return;

  const int ip = pml_map[ip_in]-1;
  const long long np = npoints;

  for (int ifield = 0; ifield <= 1; ifield++) {
    int field_dir = field_dirs[pml_dir][ifield]-1;
    grad[(ip<<ldgrad) + field_dir] = complex_number(
        // real part
        pml_c[ip_in + pml_dir*np]*(
          (1.0+pml_a[ip_in + pml_dir*np])/P_c * grad[(ip<<ldgrad) + field_dir].x +
          rs_sign * pml_b[ip_in + pml_dir*np] *
            pml_conv_plus[ip_in + np*(pml_dir + 3*field_dir)].x),
         // imaginary part
        pml_c[ip_in + pml_dir*np]*(
          (1.0+pml_a[ip_in + pml_dir*np])/P_c * grad[(ip<<ldgrad) + field_dir].y +
          rs_sign * pml_b[ip_in + pml_dir*np] *
            pml_conv_plus[ip_in + np*(pml_dir + 3*field_dir)].y));
    if (with_medium) {
      grad[(ip<<ldgrad) + field_dir + 3] = complex_number(
          // real part
          pml_c[ip_in + pml_dir*np]*(
            (1.0+pml_a[ip_in + pml_dir*np])/P_c * grad[(ip<<ldgrad) + field_dir + 3].x +
            rs_sign * pml_b[ip_in + pml_dir*np] *
              pml_conv_minus[ip_in + np*(pml_dir + 3*field_dir)].x),
          // imaginary part
          pml_c[ip_in + pml_dir*np]*(
            (1.0+pml_a[ip_in + pml_dir*np])/P_c * grad[(ip<<ldgrad) + field_dir + 3].y +
            rs_sign * pml_b[ip_in + pml_dir*np] *
              pml_conv_minus[ip_in + np*(pml_dir + 3*field_dir)].y));
    }
  }
}

__kernel void pml_update_conv(const int with_medium,
    const int npoints,
    const int pml_dir,
    __global int * restrict pml_map,
    __global double2 * restrict grad,
    const int ldgrad,
    const __global double * restrict pml_a,
    const __global double * restrict pml_b,
    __global double2 * restrict pml_conv_plus,
    __global double2 * restrict pml_conv_minus)
{
  const int field_dirs[3][2] = {{2, 3}, {1, 3}, {1, 2}};
  const int ip_in = get_global_id(0);
  if (ip_in >= npoints) return;

  const int ip = pml_map[ip_in]-1;
  const long long np = npoints;

  for (int ifield = 0; ifield <= 1; ifield++) {
    int field_dir = field_dirs[pml_dir][ifield]-1;
    pml_conv_plus[ip_in + np*(pml_dir + 3*field_dir)] = complex_number(
        // real part
        pml_a[ip_in + pml_dir*np] * grad[(ip<<ldgrad) + field_dir].x +
        pml_b[ip_in + pml_dir*np] * pml_conv_plus[ip_in + np*(pml_dir + 3*field_dir)].x,
        // imaginary part
        pml_a[ip_in + pml_dir*np] * grad[(ip<<ldgrad) + field_dir].y +
        pml_b[ip_in + pml_dir*np] * pml_conv_plus[ip_in + np*(pml_dir + 3*field_dir)].y);
    if (with_medium) {
      pml_conv_minus[ip_in + np*(pml_dir + 3*field_dir)] = complex_number(
          // real part
          pml_a[ip_in + pml_dir*np] * grad[(ip<<ldgrad) + field_dir + 3].x +
          pml_b[ip_in + pml_dir*np] * pml_conv_minus[ip_in + np*(pml_dir + 3*field_dir)].x,
          // imaginary part
          pml_a[ip_in + pml_dir*np] * grad[(ip<<ldgrad) + field_dir + 3].y +
          pml_b[ip_in + pml_dir*np] * pml_conv_minus[ip_in + np*(pml_dir + 3*field_dir)].y);
    }
  }
}

__kernel void pml_apply_new(const int npoints,
    const int jdir,
    __global int * restrict pml_map,
    __global double2 * restrict grad,
    const int ldgrad,
    const __global double * restrict pml_b,
    const __global double * restrict pml_c,
    __global double2 * restrict pml_conv_plus,
    __global double2 * restrict pml_conv_plus_old)
{
  const int idirs[3][2] = {{2, 3}, {1, 3}, {1, 2}};
  const int ip_in = get_global_id(0);
  if (ip_in >= npoints) return;

  const int ip = pml_map[ip_in]-1;
  const long long np = npoints;

  for (int i= 0; i <= 1; i++) {
    int idir = idirs[jdir][i]-1;
    pml_conv_plus[ip_in + np*(jdir + 3*idir)] = pml_c[ip_in + jdir*np] * grad[(ip<<ldgrad) + idir]
      + pml_b[ip_in + jdir*np] * pml_conv_plus_old[ip_in + np*(jdir + 3*idir)];

    grad[(ip<<ldgrad) + idir] += pml_conv_plus[ip_in + np*(jdir + 3*idir)];
  }
}

__kernel void pml_copy(const int npoints,
    __global double2 * restrict pml_conv_plus,
    __global double2 * restrict pml_conv_plus_old)
{
  const int ip_in = get_global_id(0);
  if (ip_in >= npoints) return;

  pml_conv_plus_old[ip_in] = pml_conv_plus[ip_in];
}
