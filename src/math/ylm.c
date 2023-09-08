/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

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
#include <float.h>
#include <fortran_types.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_sf_legendre.h>

/* normalization constants for the spherical harmonics */
static double sph_cnsts[16] = {
    /* l = 0 */
    0.28209479177387814347403972578,
    /* l = 1 */
    0.48860251190291992158638462283, 0.48860251190291992158638462283,
    0.48860251190291992158638462283,
    /* l = 2 */
    1.09254843059207907054338570580, 1.09254843059207907054338570580,
    0.31539156525252000603089369029, 1.09254843059207907054338570580,
    0.54627421529603953527169285290,
    /* l = 3 */
    0.59004358992664351034561027754, 2.89061144264055405538938015439,
    0.45704579946446573615802069691, 0.37317633259011539141439591319,
    0.45704579946446573615802069691, 1.44530572132027702769469007719,
    0.59004358992664351034561027754};

void oct_ylm_base(const double *r, const double *x, const double *y,
                  const double *z, const fint *l, const fint *m, double *ylm);

/* Computes real spherical harmonics ylm in the direction of vector r:
         ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
         ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
         with (theta,phi) the polar angles of r, c a positive normalization
         constant and plm associated legendre polynomials. */

void FC_FUNC_(oct_ylm, OCT_YLM)(const fint *np, const double *x,
                                const double *y, const double *z, const fint *l,
                                const fint *m, double *restrict ylm) {
  if (l[0] == 0) {
    for (int ip = 0; ip < np[0]; ip++)
      ylm[ip] = sph_cnsts[0];
    return;
  }

  for (int ip = 0; ip < np[0]; ip++) {

    double rr = sqrt(x[ip] * x[ip] + y[ip] * y[ip] + z[ip] * z[ip]);

    oct_ylm_base(&rr, &x[ip], &y[ip], &z[ip], l, m, &ylm[ip]);
  }
}

void FC_FUNC_(oct_ylm_vector, OCT_YLM_VECTOR)(const fint *np, const double *xx,
                                              const double *rr, const fint *l,
                                              const fint *m,
                                              double *restrict ylm) {
  if (l[0] == 0) {
    for (int ip = 0; ip < np[0]; ip++)
      ylm[ip] = sph_cnsts[0];
    return;
  }

  for (int ip = 0; ip < np[0]; ip++) {

    oct_ylm_base(&rr[ip], &xx[ip * 3], &xx[ip * 3 + 1], &xx[ip * 3 + 2], l, m,
                 &ylm[ip]);
  }
}

void oct_ylm_base(const double *r, const double *x, const double *y,
                  const double *z, const fint *l, const fint *m, double *ylm) {

  double r2, r3, rr, rx, ry, rz, cosphi, sinphi, cosm, sinm, phase;
  int i;

  r2 = r[0] * r[0];

  /* if r=0, direction is undefined => make ylm=0 except for l=0 */
  if (r2 < 1.0e-30) {
    ylm[0] = 0.0;
    return;
  }

  switch (l[0]) {
  case 1:
    rr = 1.0 / r[0];
    switch (m[0]) {
    case -1:
      ylm[0] = -sph_cnsts[1] * rr * y[0];
      return;
    case 0:
      ylm[0] = sph_cnsts[2] * rr * z[0];
      return;
    case 1:
      ylm[0] = -sph_cnsts[3] * rr * x[0];
      return;
    }
  case 2:
    switch (m[0]) {
    case -2:
      ylm[0] = sph_cnsts[4] * x[0] * y[0] / r2;
      return;
    case -1:
      ylm[0] = -sph_cnsts[5] * y[0] * z[0] / r2;
      return;
    case 0:
      ylm[0] = sph_cnsts[6] * (3.0 * z[0] * z[0] / r2 - 1.0);
      return;
    case 1:
      ylm[0] = -sph_cnsts[7] * x[0] * z[0] / r2;
      return;
    case 2:
      ylm[0] = sph_cnsts[8] * (x[0] * x[0] - y[0] * y[0]) / r2;
      return;
    }
  case 3:
    r3 = r2 * r[0];
    switch (m[0]) {
    case -3:
      ylm[0] = sph_cnsts[9] * ((3.0 * x[0] * x[0] - y[0] * y[0]) * y[0]) / r3;
      return;
    case -2:
      ylm[0] = sph_cnsts[10] * (x[0] * y[0] * z[0]) / r3;
      return;
    case -1:
      ylm[0] = sph_cnsts[11] * (y[0] * (5.0 * z[0] * z[0] - r2)) / r3;
      return;
    case 0:
      ylm[0] = sph_cnsts[12] * (z[0] * (5.0 * z[0] * z[0] - 3.0 * r2)) / r3;
      return;
    case 1:
      ylm[0] = sph_cnsts[13] * (x[0] * (5.0 * z[0] * z[0] - r2)) / r3;
      return;
    case 2:
      ylm[0] = sph_cnsts[14] * ((x[0] * x[0] - y[0] * y[0]) * z[0]) / r3;
      return;
    case 3:
      ylm[0] = sph_cnsts[15] * ((x[0] * x[0] - 3.0 * y[0] * y[0]) * x[0]) / r3;
      return;
    }
  }

  /* get phase */
  rr = 1.0 / r[0];

  rx = x[0] * rr;
  ry = y[0] * rr;
  rz = z[0] * rr;

  rr = hypot(rx, ry);
  if (rr < 1e-20) {
    cosphi = 0.0;
    sinphi = 0.0;
  } else {
    cosphi = rx / rr;
    sinphi = ry / rr;
  }

  /* compute sin(mphi) and cos(mphi) by adding cos/sin */
  cosm = 1.;
  sinm = 0.;
  for (i = 0; i < abs(m[0]); i++) {
    double a = cosm, b = sinm;
    cosm = a * cosphi - b * sinphi;
    sinm = a * sinphi + b * cosphi;
  }
  phase = m[0] < 0 ? sinm : cosm;
  phase = m[0] == 0 ? phase : sqrt(2.0) * phase;

  /* lower bound with a small number (~= 10^-308) to avoid floating invalids */
  rz = rz < DBL_MIN ? DBL_MIN : rz;

  rr = gsl_sf_legendre_sphPlm(l[0], abs(m[0]), rz);

  /* I am not sure whether we are including the Condon-Shortley factor (-1)^m */
  ylm[0] = rr * phase;
}
