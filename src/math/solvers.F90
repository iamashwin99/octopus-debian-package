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

!> This module is intended to contain "only mathematical" functions
!! and procedures.

module solvers_oct_m
  use blas_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                     &
    dconjugate_gradients,       &
    zconjugate_gradients,       &
    zqmr_sym_gen_dotu,          &
    dqmr_sym_gen_dotu,          &
    zqmr_gen_dotu,              &
    dqmr_gen_dotu,              &
    didrs,                      &
    zidrs


  !> ---------------------------------------------------------
  !! QMR (quasi-minimal residual) algorithm for complex symmetric matrices
  !! algorithm taken from:
  !! Parallel implementation of efficient preconditioned linear solver for
  !! grid-based applications in chemical physics. II: QMR linear solver
  !! Appendix A. Simplified QMR algorithm
  !! W Chen and B Poirier, J Comput Phys 219, 198-209 (2006)

  !> ---------------------------------------------------------
  !! QMR (quasi-minimal residual) algorithm for complex matrices
  !! algorithm taken from: An Implementation of the QMR Method based on
  !! Coupled Two-Term Recurrences by R. W. Freund and N. M. Nachtigal (page 25)
  !! http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950017192_1995117192.pdf

  interface dconjugate_gradients
    module procedure dsym_conjugate_gradients, dbi_conjugate_gradients
  end interface dconjugate_gradients

  interface zconjugate_gradients
    module procedure zsym_conjugate_gradients, zbi_conjugate_gradients
  end interface zconjugate_gradients

contains

#include "undef.F90"
#include "complex.F90"
#include "solvers_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "solvers_inc.F90"

end module solvers_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
