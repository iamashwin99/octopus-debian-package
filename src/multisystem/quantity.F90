!! Copyright (C) 2020 N. Tancogne-Dejean
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

!> This module defines the quantity_t class and the IDs for quantities, which can be exposed by a system,
!! and used by an interaction.
module quantity_oct_m
  use clock_oct_m
  implicit none

  private
  public ::                   &
    quantity_t

  !# doc_start quantity
  integer, public, parameter ::         &
    POSITION                     =  1,  &
    VELOCITY                     =  2,  &
    CURRENT                      =  3,  &
    DENSITY                      =  4,  &
    SCALAR_POTENTIAL             =  5,  &
    VECTOR_POTENTIAL             =  6,  &
    E_FIELD                      =  7,  &
    B_FIELD                      =  8,  &
    MASS                         =  9,  &
    CHARGE                       = 10,  &
    PERMITTIVITY                 = 11,  &
    PERMEABILITY                 = 12,  &
    E_CONDUCTIVITY               = 13,  &
    M_CONDUCTIVITY               = 14,  &
    MAX_QUANTITIES               = 14
  !# doc_end

  character(len=17), public, parameter :: QUANTITY_LABEL(MAX_QUANTITIES) = (/ &
    "position        ", &
    "velocity        ", &
    "current         ", &
    "density         ", &
    "scalar potential", &
    "vector potential", &
    "E field         ", &
    "B field         ", &
    "mass            ", &
    "charge          ", &
    "permittivity    ", &
    "permeability    ", &
    "e_conductivity  ", &
    "m_conductivity  "  &
    /)

  !> Systems (system_t) can expose quantities that can be used to calculate interactions
  !! with other systems.
  !!
  !! Some quantities are dynamical variables of the system. Such quantities are
  !! usually updated by the propagation algorithm and cannot be calculated
  !! on-demand. Such quantities must be marked as "protected".
  type quantity_t
    private
    type(clock_t), public :: clock               !< Clock storing the time at which the quantity was last updated.
    logical,       public :: required = .false.  !< Should this quantities be calculated?
    logical,       public :: available_at_any_time = .false. !< Can we use this quantity at any requested time? (e.g., this will be true for a static quantity, but false for a quantity that is only updated at specific time-steps)
    logical,       public :: updated_on_demand = .true. !< If true, the quantity is only updated when requested. If false, the quantity is updated automatically during the execution of an algorithm.
  end type quantity_t

contains

end module quantity_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
