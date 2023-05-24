!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module mesh_function_oct_m
  use blas_oct_m
  use comm_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use loct_math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use qshep_oct_m
  use quickrnd_oct_m
  use splines_oct_m

  implicit none

  private
  public ::                &
    dmf_integrate,         &
    zmf_integrate,         &
    dmf_dotp,              &
    zmf_dotp,              &
    dmf_nrm2,              &
    zmf_nrm2,              &
    dmf_moment,            &
    zmf_moment,            &
    dmf_random,            &
    zmf_random,            &
    dmf_interpolate_points,&
    zmf_interpolate_points,&
    mf_surface_integral,   &
    mf_line_integral,      &
    dmf_dotp_aux,          &
    zmf_dotp_aux,          &
    dmf_multipoles,        &
    zmf_multipoles,        &
    dmf_dipole,            &
    zmf_dipole,            &
    mesh_init_mesh_aux,    &
    dmf_dotu_aux,          &
    zmf_dotu_aux,          &
    dmf_nrm2_aux,          &
    zmf_nrm2_aux,          &
    dmf_normalize,         &
    zmf_normalize,         &
    zmf_fix_phase

  ! These variables are to be used by the "distdot" function, that is outside the module
  ! but inside this file.
  ! FIXME: This is very ugly, at least these values should be set by a function.
  public :: mesh_aux
  logical, public :: sp_parallel
  integer, public :: sp_np, sp_dim, sp_st1, sp_st2, sp_kp1, sp_kp2
  integer, public :: sp_distdot_mode
  type(mpi_grp_t), public :: sp_grp

  interface mf_surface_integral
    module procedure dmf_surface_integral_scalar, dmf_surface_integral_vector, &
      zmf_surface_integral_scalar, zmf_surface_integral_vector
  end interface mf_surface_integral

  interface mf_line_integral
    module procedure dmf_line_integral_scalar, dmf_line_integral_vector, &
      zmf_line_integral_scalar, zmf_line_integral_vector
  end interface mf_line_integral

  interface dmf_dotp
    module procedure dmf_dotp_1, dmf_dotp_2
  end interface dmf_dotp

  interface zmf_dotp
    module procedure zmf_dotp_1, zmf_dotp_2
  end interface zmf_dotp

  interface dmf_nrm2
    module procedure dmf_nrm2_1, dmf_nrm2_2
  end interface dmf_nrm2

  interface zmf_nrm2
    module procedure zmf_nrm2_1, zmf_nrm2_2
  end interface zmf_nrm2

  type(mesh_t), pointer :: mesh_aux => null()

  type(profile_t), save ::        &
    dPROFILING_MF_INTEGRATE,      &
    zPROFILING_MF_INTEGRATE,      &
    dPROFILING_MF_DOTP,           &
    zPROFILING_MF_DOTP,           &
    dPROFILING_MF_REDUCE,         &
    zPROFILING_MF_REDUCE,         &
    dPROFILING_MF_NRM2,           &
    zPROFILING_MF_NRM2

contains

  subroutine mesh_init_mesh_aux(mesh)
    type(mesh_t), target, intent(in) :: mesh

    PUSH_SUB(mesh_init_mesh_aux)
    mesh_aux => mesh

    POP_SUB(mesh_init_mesh_aux)
  end subroutine mesh_init_mesh_aux

  ! Given a complex function on the grid, we want to minimize its imaginary part
  ! For this, we determine the angle that will minimize the imaginary part after rotation
  ! in the complex plane.
  ! We then choose the solution that maximizes the real part
  subroutine zmf_fix_phase(mesh, ff)
    type(mesh_t),      intent(in)     :: mesh
    CMPLX,             intent(inout)  :: ff(:)

    FLOAT, allocatable :: ff_re(:), ff_im(:)
    FLOAT :: dot_re, dot_im, dot_reim, sol1, sol2, maxv, maxvloc
    FLOAT :: theta, ctheta, stheta
    integer :: ip, ipmax
    CMPLX   :: zval

    PUSH_SUB(zmf_fix_phase)

    SAFE_ALLOCATE(ff_re(1:mesh%np))
    SAFE_ALLOCATE(ff_im(1:mesh%np))

    ff_re = TOFLOAT(ff)
    ff_im = aimag(ff)

    dot_re = dmf_nrm2(mesh, ff_re)**2
    dot_im = dmf_nrm2(mesh, ff_im)**2
    dot_reim = dmf_dotp(mesh, ff_re, ff_im)

    ! We have tan(2\theta) = 2 (\int \psi_R \psi_I)/( \int \psi_R^2-\psi_I^2)
    if(dot_re+dot_im > CNST(1e-8)) then
      if(M_TWO*dot_reim/(dot_re-dot_im) > CNST(1e-8)) then
        theta = M_HALF * atan(M_TWO*dot_reim/(dot_re-dot_im))

        ! Get the maximum for the real part
        sol1 = cos(theta)**2 * dot_re + sin(theta)**2 * dot_im + M_TWO*sin(theta)*cos(theta)*dot_reim
        sol2 = sin(theta)**2 * dot_re + cos(theta)**2 * dot_re - M_TWO*sin(theta)*cos(theta)*dot_reim
        if(sol2>sol1) theta = theta + M_HALF*M_PI

        ctheta = cos(theta)
        stheta = sin(theta)

        !We do the rotation
        !$omp parallel do
        do ip = 1, mesh%np
          ff(ip) = TOCMPLX(ctheta*ff_re(ip)-stheta*ff_im(ip),stheta*ff_re(ip)+ctheta*ff_im(ip))
        end do

      else ! We do not know what is the direction here
        !    we set it using the maximum value of the function from rank 0
        ipmax = 0
        maxv = M_ZERO
        do ip = 1, mesh%np
          maxvloc = TOFLOAT(ff(ip)*conjg(ff(ip)))
          if(maxvloc > maxv) then
            ipmax = ip
            maxv = maxvloc
          end if
        end do
        zval = ff(ipmax)

        !We want to use the same phase on all processors
        if(mesh%parallel_in_domains) then
          call mesh%mpi_grp%bcast(zval, 1, MPI_FLOAT, 0)
        end if

        call lalg_scal(mesh%np, abs(zval)/zval, ff)

      end if
    else
      message(1) = "Internal error fixing phase: mesh function has zero norm."
      call messages_fatal(1)
    end if

    SAFE_DEALLOCATE_A(ff_re)
    SAFE_DEALLOCATE_A(ff_im)

    POP_SUB(zmf_fix_phase)
  end subroutine zmf_fix_phase

#include "undef.F90"
#include "real.F90"
#include "mesh_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_function_inc.F90"

end module mesh_function_oct_m

#ifdef HAVE_SPARSKIT

!> This function will be linked by SPARSKIT to perform dot products.
!! It expects complex numbers as an array with first real parts, then imaginary parts.
! ---------------------------------------------------------
REAL_DOUBLE function distdot(n, x, ix, y, iy)
use comm_oct_m
use global_oct_m
use messages_oct_m
use mesh_function_oct_m
use profiling_oct_m

implicit none

integer,     intent(in) :: n
REAL_DOUBLE, intent(in) :: x(n)
integer,     intent(in) :: ix
REAL_DOUBLE, intent(in) :: y(n)
integer,     intent(in) :: iy

integer :: j, ik, ist, idim, k

 ! SPARSKIT only calls this function with ix, iy = 1 i.e. no stride.
ASSERT(ix == 1)
ASSERT(iy == 1)

select case (sp_distdot_mode)
case (1)
  distdot = dmf_dotp_aux(x(1:n/2), y(1:n/2)) + dmf_dotp_aux(x(n/2+1:n), y(n/2+1:n))

case (2)
  distdot = M_ZERO
  j = 1
  k = n/2+1
  do ik = sp_kp1, sp_kp2
    do ist = sp_st1, sp_st2
      do idim = 1, sp_dim
        distdot = distdot + dmf_dotp_aux(x(j: j+sp_np-1), y(j:j+sp_np-1))
        distdot = distdot + dmf_dotp_aux(x(k: k+sp_np-1), y(k:k+sp_np-1))
        j = j + sp_np
        k = k + sp_np
      end do
    end do
  end do
  if (sp_parallel) call comm_allreduce(sp_grp, distdot)

case (3)
  distdot = M_ZERO
  j = 1
  k = n/2+1
  do ik = sp_kp1, sp_kp2
    do ist = sp_st1, sp_st2
      do idim = 1, sp_dim
        distdot = distdot + dmf_dotp_aux(x(j: j+sp_np-1), y(j:j+sp_np-1))
        distdot = distdot + dmf_dotp_aux(x(k: k+sp_np-1), y(k:k+sp_np-1))
        j = j + sp_np
        k = k + sp_np
      end do
    end do
  end do
  do ik = sp_kp1, sp_kp2
    do ist = sp_st1, sp_st2
      do idim = 1, sp_dim
        distdot = distdot + dmf_dotp_aux(x(j: j+sp_np-1), y(j:j+sp_np-1))
        distdot = distdot + dmf_dotp_aux(x(k: k+sp_np-1), y(k:k+sp_np-1))
        j = j + sp_np
        k = k + sp_np
      end do
    end do
  end do
  if (sp_parallel) call comm_allreduce(sp_grp, distdot)

end select

end function distdot

#endif /* HAVE_SPARSKIT */

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
