!! Copyright (C) 2008 X. Andrade
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

module eigen_rmmdiis_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use debug_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use preconditioners_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                     &
    deigensolver_rmmdiis,       &
    zeigensolver_rmmdiis,       &
    deigensolver_rmmdiis_min,   &
    zeigensolver_rmmdiis_min

  type batch_pointer_t
    private
    type(wfs_elec_t), pointer :: batch
  end type batch_pointer_t

  type(profile_t), save :: prof, prof_iter
  type(profile_t), save :: prof_lc

contains

  subroutine find_lambda(ca, cb, cc, lambda, ik, ist)
    FLOAT,   intent(in)  :: ca, cb, cc
    FLOAT,   intent(out) :: lambda
    integer, intent(in)  :: ik, ist

    FLOAT :: qq

    PUSH_SUB(find_lambda)


    ! We follow numerical recipes here to get the two stable solutions
    qq = -M_HALF * (cb + sign(sqrt(cb**2 - M_FOUR*ca*cc), cb))
    lambda = qq/ca

    ! To find the minimum, we want that the derivative changes from negative to positive
    ! Else, we pick the other solution
    if(ca*(M_TWO*lambda)+cb < M_ZERO) then ! This is a maximum, so we pick the other solution
      lambda = cc / qq
    end if

    if(debug%info) then
      write(message(1), '(2(a,i4),4(a,es12.5))') &
      'Debug: RMMDIIS Eigensolver - ik', ik, ' ist ', ist, &
        ' lambda= ', lambda, ' a= ', ca, ' b= ', cb, ' c= ', cc
      call messages_info(1)
    end if

    ! restrict the value of lambda to be between 0.1 and 1000
    ! This is not the same range as in the original RMMDIIS paper,
    ! but this is ont important as our preconditioner is different, so the norm of the preconditioned
    ! residue is also changed
    if (abs(lambda) > CNST(1000.0)) lambda = sign(CNST(1000.0),lambda)
    if (abs(lambda) < CNST(0.1)) lambda = sign(CNST(0.1),lambda)
 
    POP_SUB(find_lambda)
  end subroutine find_lambda

#include "real.F90"
#include "eigen_rmmdiis_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "eigen_rmmdiis_inc.F90"
#include "undef.F90"

end module eigen_rmmdiis_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
