!! Copyright (C) 2019 N. Tancogne-Dejean
!! Copyright (C) 2020 M. Oliveira, Heiko Appel
!! Copyright (C) 2021 S. Ohlmann
!! Copyright (C) 2023 M. Oliveira
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

module propagator_factory_oct_m
  use algorithm_factory_oct_m
  use algorithm_oct_m
  use debug_oct_m
  use global_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use propagator_aetrs_oct_m
  use propagator_beeman_oct_m
  use propagator_bomd_oct_m
  use propagator_exp_mid_oct_m
  use propagator_exp_mid_2step_oct_m
  use propagator_leapfrog_oct_m
  use propagator_oct_m
  use propagator_rk4_oct_m
  use propagator_static_oct_m
  use propagator_verlet_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use system_oct_m
  implicit none

  private
  public :: propagator_factory_t

  ! Known multisystem propagators
  integer, public, parameter ::       &
    PROP_STATIC                 = 0,  &
    PROP_VERLET                 = 1,  &
    PROP_BEEMAN                 = 2,  &
    PROP_BEEMAN_SCF             = 3,  &
    PROP_EXPMID_2STEP           = 4,  &
    PROP_EXPMID_2STEP_SCF       = 5,  &
    PROP_AETRS_MS               = 6,  &
    PROP_RK4                    = 7,  &
    PROP_EXPMID                 = 8,  &
    PROP_LEAPFROG               = 9,  &
    PROP_BOMD                   = 10

  type, extends(algorithm_factory_t) :: propagator_factory_t
    private
    FLOAT :: final_time !< Final propagation time
  contains
    procedure :: create => propagator_factory_create
    procedure :: create_static => propagator_factory_create_static
  end type propagator_factory_t

  interface propagator_factory_t
    module procedure propagator_factory_constructor
  end interface propagator_factory_t


contains
  ! ---------------------------------------------------------------------------------------
  function propagator_factory_constructor(namespace) result(factory)
    type(namespace_t), intent(in) :: namespace
    type(propagator_factory_t) :: factory

    PUSH_SUB(propagator_factory_constructor)

    ! Get final propagation time from input
    ! This variable is also defined (and properly documented) in td/td.F90.
    ! This is temporary, until all the propagators are moved to the new framework.
    call parse_variable(namespace, 'TDPropagationTime', CNST(-1.0), factory%final_time, unit = units_inp%time)
    if (factory%final_time <= M_ZERO) then
      call messages_input_error(namespace, 'TDPropagationTime', 'must be greater than zero')
    end if
    call messages_print_var_value('TDPropagationTime', factory%final_time)

    POP_SUB(propagator_factory_constructor)
  end function propagator_factory_constructor

  ! ---------------------------------------------------------------------------------------
  function propagator_factory_create(this, system) result(algorithm)
    class(propagator_factory_t),  intent(in) :: this
    class(interaction_partner_t), intent(in), target :: system
    class(algorithm_t), pointer :: algorithm

    integer :: prop_type
    FLOAT :: dt
    class(propagator_t), pointer :: propagator

    PUSH_SUB(propagator_factory_create)

    !%Variable TDSystemPropagator
    !%Type integer
    !%Default static
    !%Section Time-Dependent::Propagation
    !%Description
    !% A variable to set the propagator in the multisystem framework.
    !% This is a temporary solution, and should be replaced by the
    !% TDPropagator variable.
    !%Option static 0
    !% (Experimental) Do not propagate the system in time.
    !%Option verlet 1
    !% (Experimental) Verlet propagator.
    !%Option beeman 2
    !% (Experimental) Beeman propagator without predictor-corrector.
    !%Option beeman_scf 3
    !% (Experimental) Beeman propagator with predictor-corrector scheme.
    !%Option exp_mid_2step 4
    !% (Experimental) Exponential midpoint propagator without predictor-corrector.
    !%Option exp_mid_2step_scf 5
    !% (Experimental) Exponential midpoint propagator with predictor-corrector scheme.
    !%Option prop_aetrs 6
    !% (Experimental) Approximate ETRS propagator
    !%Option prop_rk4 7
    !% (Experimental) RK4 propagator
    !%Option prop_expmid 8
    !% (Experimental) Exponential midpoint propagator with extrapolation.
    !%Option prop_leapfrog 9
    !% (Experimental) Leap frog algorithm
    !%Option prop_bomd 10
    !% (Experimental) Born-Oppenheimer MD propagator for the matter system.
    !%End
    call parse_variable(system%namespace, 'TDSystemPropagator', PROP_STATIC, prop_type)
    if (.not. varinfo_valid_option('TDSystemPropagator', prop_type)) then
      call messages_input_error(system%namespace, 'TDSystemPropagator')
    end if
    call messages_print_var_option('TDSystemPropagator', prop_type, namespace=system%namespace)

    dt = propagator_factory_read_dt(this, system%namespace)

    select case (prop_type)
    case (PROP_STATIC)
      propagator => propagator_static_t(dt, 1)
    case (PROP_VERLET)
      propagator => propagator_verlet_t(dt)
    case (PROP_BEEMAN)
      propagator => propagator_beeman_t(dt, predictor_corrector=.false.)
    case (PROP_BEEMAN_SCF)
      propagator => propagator_beeman_t(dt, predictor_corrector=.true.)
    case (PROP_EXPMID_2STEP)
      propagator => propagator_exp_mid_2step_t(dt, predictor_corrector=.false.)
    case (PROP_EXPMID_2STEP_SCF)
      propagator => propagator_exp_mid_2step_t(dt, predictor_corrector=.true.)
    case (PROP_AETRS_MS)
      propagator => propagator_aetrs_t(dt)
    case (PROP_RK4)
      propagator => propagator_rk4_t(dt)
    case (PROP_EXPMID)
      propagator => propagator_exp_mid_t(dt)
    case (PROP_LEAPFROG)
      propagator => propagator_leapfrog_t(dt)
    case (PROP_BOMD)
      propagator => propagator_bomd_t(dt)
    case default
      call messages_input_error(system%namespace, 'TDSystemPropagator')
    end select

    propagator%final_time = this%final_time

    select type (system)
    class is (system_t)
      propagator%system => system
    class default
      ASSERT(.false.)
    end select

    algorithm => propagator

    POP_SUB(propagator_factory_create)
  end function propagator_factory_create

  ! ---------------------------------------------------------------------------------------
  function propagator_factory_create_static(this, system) result(algorithm)
    class(propagator_factory_t),  intent(in) :: this
    class(interaction_partner_t), intent(in), target :: system
    class(algorithm_t), pointer :: algorithm

    FLOAT :: largest_dt
    integer :: nsteps
    class(propagator_t), pointer :: propagator

    PUSH_SUB(propagator_factory_create_static)

    select type (system)
    type is (multisystem_basic_t)
      ! The multisystem container is an exception, as its time-step and the
      ! number of algorithmic steps is determined by the subsystems.
      largest_dt = system%largest_dt()
      nsteps = int(largest_dt/system%smallest_algo_dt())
      propagator => propagator_static_t(largest_dt, nsteps)

    class default
      algorithm => propagator_static_t(propagator_factory_read_dt(this, system%namespace), 1)
    end select

    propagator%final_time = this%final_time

    select type (system)
    class is (system_t)
      propagator%system => system
    class default
      ASSERT(.false.)
    end select

    algorithm => propagator

    POP_SUB(propagator_factory_create_static)
  end function propagator_factory_create_static

  ! ---------------------------------------------------------------------------------------
  FLOAT function propagator_factory_read_dt(this, namespace) result(dt)
    class(propagator_factory_t),  intent(in) :: this
    type(namespace_t),            intent(in) :: namespace

    ! This variable is also defined (and properly documented) in td/td.F90.
    ! This is temporary, until all the propagators are moved to the new framework.
    call parse_variable(namespace, 'TDTimeStep', CNST(10.0), dt)
    if (dt <= M_ZERO) then
      call messages_input_error(namespace, 'TDTimeStep', "must be greater than zero")
    end if
    call messages_print_var_value('TDTimeStep', dt, namespace=namespace)

  end function propagator_factory_read_dt

end module propagator_factory_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
