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

module physics_op_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::            &
    dphysics_op_L,     &
    zphysics_op_L,     &
    dphysics_op_L2,    &
    zphysics_op_L2

contains

#include "undef.F90"
#include "real.F90"
#include "physics_op_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "physics_op_inc.F90"


end module physics_op_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
