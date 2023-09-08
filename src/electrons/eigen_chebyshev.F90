!! Copyright (C) 2023. A Buccheri.
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

module eigen_chebyshev_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use debug_oct_m
  use chebyshev_filter_bounds_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use kind_oct_m, only: dp => r8
  use lalg_adv_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use subspace_oct_m
  use wfs_elec_oct_m

  implicit none
  private

  !> @class Chebyshev input parameters.
  type eigen_chebyshev_t
    integer  :: n_lanczos              !< Number of Lanczos iterations used to determine the  estimate upper bound of H
    integer :: degree                  !< Degree of Chebyshev polynomial
    FLOAT :: bound_mixing              !< Coefficient for linear mixing of min and max eigenvalues,
    !                                  !< for approximation of the spectral filter''s lower bound:
    !                                  !< lower_bound = bm * min(e_approx) + (bm - 1) * max(eig_approx)
    !                                  !< such that bm = 0 => lower_bound = max(eig_approx)
    !                                  !<           bm = 1 => lower_bound = min(eig_approx)
    integer :: n_iter                  !< Number of iterations used for the first SCF step
  end type eigen_chebyshev_t

  !> Default Chebyshev input parameters
  !> Arguments 1 and 2 taken from 10.1016/j.jcp.2006.03.017
  !> Argument 3 taken from 10.1016/j.jcp.2014.06.056
  !> Argument 4 set empirically. A value > 2 results in a more accurate density.
  !> Values > 6 show minimal to no improvement in the number of SCF steps.
  type(eigen_chebyshev_t), protected :: default_chebyshev_params = eigen_chebyshev_t(5, 10, M_HALF, 5)

  type batch_pointer_t
    private
    type(wfs_elec_t), pointer :: batch
  end type batch_pointer_t


  public :: &
    eigen_chebyshev_t, &
    default_chebyshev_params, &
    dchebyshev_filter_solver, zchebyshev_filter_solver

contains

#include "real.F90"
#include "eigen_chebyshev_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "eigen_chebyshev_inc.F90"
#include "undef.F90"

end module eigen_chebyshev_oct_m
