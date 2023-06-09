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


 ! ---------------------------------------------------------
 ! Tries to avoid ill-defined combinations of run modes.
 ! ---------------------------------------------------------
subroutine check_faulty_runmodes(sys, tr)
  type(electrons_t),       intent(in) :: sys
  type(propagator_base_t), intent(in) :: tr

  integer :: no_electrons, n_filled, n_partially_filled, n_half_filled

  PUSH_SUB(check_faulty_runmodes)

  ! No QOCT runs with periodic boundary conditions.
  if (sys%space%is_periodic()) then
    write(message(1), '(a)') 'No QOCT runs with periodic boundary conditions. '
    call messages_fatal(1, namespace=sys%namespace)
  end if

  ! This should check that we really have occupation one for
  ! one of the spin-orbitals, and occupation zero for all the others.
  ! Otherwise the algorithms are bound to fail.
  select case (sys%st%d%ispin)
  case (UNPOLARIZED)
    call occupied_states(sys%st, sys%namespace, 1, n_filled, n_partially_filled, n_half_filled)
    no_electrons = 2*n_filled + n_half_filled
    if (n_partially_filled > 0) then
      write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
      call messages_fatal(1, namespace=sys%namespace)
    end if
  case (SPIN_POLARIZED)
    call occupied_states(sys%st, sys%namespace, 1, n_filled, n_partially_filled, n_half_filled)
    if (n_partially_filled > 0 .or. n_half_filled > 0) then
      write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
      call messages_fatal(1, namespace=sys%namespace)
    end if
    no_electrons = n_filled
    call occupied_states(sys%st, sys%namespace, 2, n_filled, n_partially_filled, n_half_filled)
    no_electrons = n_filled + no_electrons
    if (n_partially_filled > 0 .or. n_half_filled > 0) then
      write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
      call messages_fatal(1, namespace=sys%namespace)
    end if
  case (SPINORS)
    call occupied_states(sys%st, sys%namespace, 1, n_filled, n_partially_filled, n_half_filled)
    no_electrons = n_filled
    if (n_partially_filled > 0 .or. n_half_filled > 0) then
      write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
      call messages_fatal(1, namespace=sys%namespace)
    end if
  end select

  if (abs(sys%st%qtot - TOFLOAT(no_electrons)) > CNST(1.0e-8)) then
    write(message(1), '(a)') 'Error in check_faulty_runmodes'
    call messages_fatal(1, namespace=sys%namespace)
  end if

  if (oct%algorithm == OPTION__OCTSCHEME__OCT_ZBR98) then
    select case (target_type(oct_target))
    case (oct_tg_groundstate, oct_tg_gstransformation, &
      oct_tg_userdefined)
    case default
      write(message(1), '(a)') 'The scheme "OCTScheme = oct_zbr98 can only be used if'
      write(message(2), '(a)') 'the target state is "OCTTargetOperator = oct_tg_gstransformation"'
      write(message(3), '(a)') 'or "OCTTargetOperator = oct_tg_groundstate"'
      write(message(4), '(a)') 'or "OCTTargetOperator = oct_tg_userdefined".'
      call messages_fatal(4, namespace=sys%namespace)
    end select
  end if

  ! Filters only with the WG05 scheme.
  if (filter_number(filter) /= 0) then
    if (oct%algorithm /= OPTION__OCTSCHEME__OCT_WG05) then
      write(message(1), '(a)') 'Filters can only be used with the WG05 QOCT algorithm.'
      call messages_fatal(1, namespace=sys%namespace)
    end if
  end if

  ! local targets only in ZR98 and WG05
  if (target_type(oct_target) == oct_tg_local .or. &
    target_type(oct_target) == oct_tg_jdensity .or. &
    target_type(oct_target) == oct_tg_td_local) then
    if (oct%algorithm == OPTION__OCTSCHEME__OCT_ZBR98) then
      write(message(1), '(a)') 'Cannot use ZBR98 OCT scheme if the target is oct_tg_jdensity,'
      write(message(2), '(a)') 'oct_tg_local or oct_tg_td_local.'
      call messages_fatal(2, namespace=sys%namespace)
    end if
  end if

  ! the inh term in the bwd evolution of chi is taken into
  ! consideration only for certain propagators
  if (.not. oct_algorithm_is_direct(oct)) then
    if (target_mode(oct_target) == oct_targetmode_td) then
      select case (tr%method)
      case (PROP_CRANK_NICOLSON)
      case (PROP_QOCT_TDDFT_PROPAGATOR)
        select case (tr%te%exp_method)
        case (EXP_TAYLOR)
        case default
          write(message(1), '(a)') 'If you use time-dependent target, and you set'
          write(message(2), '(a)') '"TDPropagator = qoct_tddft_propagator", '
          write(message(3), '(a)') 'then you must set "TDExponentialMethod = taylor".'
          call messages_fatal(3, namespace=sys%namespace)
        end select
      case (PROP_EXPONENTIAL_MIDPOINT)
        select case (tr%te%exp_method)
        case (EXP_LANCZOS)
        case default
          write(message(1), '(a)') 'If you use time-dependent target, and you set'
          write(message(2), '(a)') '"TDPropagator = exp_mid", '
          write(message(3), '(a)') 'then you must set "TDExponentialMethod = lanczos".'
          call messages_fatal(3, namespace=sys%namespace)
        end select
      case default
        write(message(1), '(a)') 'If you use time-dependent target, then you must set'
        write(message(2), '(a)') '"TDPropagator = crank_nicolson", '
        write(message(3), '(a)') '"TDPropagator = qoct_tddft_propagator", or'
        write(message(4), '(a)') '"TDPropagator = exp_mid".'
        call messages_fatal(4, namespace=sys%namespace)
      end select
    end if
  end if


  if (target_type(oct_target) == oct_tg_excited) then
    if (sys%st%d%ispin == UNPOLARIZED) then
      write(message(1), '(a)') 'If OCTTargetMode = oct_tg_excited, then you must run either with'
      write(message(1), '(a)') 'SpinComponents = spin_polarized or SpinComponents = spinors.'
      call messages_fatal(2, namespace=sys%namespace)
    end if
  end if

  if (sys%hm%theory_level /= INDEPENDENT_PARTICLES) then
    if (sys%hm%theory_level /= KOHN_SHAM_DFT) then
      write(message(1), '(a)') 'In optimal control theory mode, you can only use either independent'
      write(message(2), '(a)') 'particles "TheoryLevel = independent_particles", or Kohn-Sham DFT'
      write(message(3), '(a)') '"TheoryLevel = kohn_sham".'
      call messages_fatal(3, namespace=sys%namespace)
    end if
    if ((tr%method /= PROP_QOCT_TDDFT_PROPAGATOR) .and. &
      (tr%method /= PROP_EXPLICIT_RUNGE_KUTTA4)   .and. &
      (tr%method /= PROP_RUNGE_KUTTA2)) then
      if (.not. oct_algorithm_is_direct(oct)) then
        write(message(1), '(a)') 'When doing QOCT with interacting electrons, then you must set'
        write(message(2), '(a)') 'TDPropagator = qoct_tddft_propagator'
        call messages_fatal(2, namespace=sys%namespace)
      end if
    end if
  end if

  if (sys%hm%abs_boundaries%abtype == MASK_ABSORBING) then
    if ((oct%algorithm /= OPTION__OCTSCHEME__OCT_DIRECT) .and. &
      (oct%algorithm /= OPTION__OCTSCHEME__OCT_NLOPT_BOBYQA)) then
      write(message(1), '(a)') 'Cannot do QOCT with mask absorbing boundaries. Use either'
      write(message(2), '(a)') '"AbsorbingBoundaries = cap" or "AbsorbingBoundaries = no".'
      call messages_fatal(2, namespace=sys%namespace)
    end if
  end if

  if (target_type(oct_target) == oct_tg_exclude_state) then
    if (no_electrons > 1) then
      write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
      write(message(2), '(a)') 'one-electron runs.'
      call messages_fatal(2, namespace=sys%namespace)
    end if
    if (sys%st%d%ispin == SPIN_POLARIZED) then
      write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
      write(message(2), '(a)') 'runs in spin restricted, or in spinors mode (spin-polarized is'
      write(message(3), '(a)') 'is not allowed.'
      call messages_fatal(3, namespace=sys%namespace)
    end if
  end if

  if (target_type(oct_target) == oct_tg_velocity) then
    if ((oct%algorithm /= OPTION__OCTSCHEME__OCT_DIRECT) .and. &
      (oct%algorithm /= OPTION__OCTSCHEME__OCT_NLOPT_BOBYQA) .and. &
      (oct%algorithm /= OPTION__OCTSCHEME__OCT_BFGS) .and. &
      (oct%algorithm /= OPTION__OCTSCHEME__OCT_CG)) then
      write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_velocity", you can only use'
      write(message(2), '(a)') '"OCTScheme = oct_direct" or'
      write(message(3), '(a)') '"OCTScheme = oct_bobyqa" or'
      write(message(4), '(a)') '"OCTScheme = oct_cg" for the optimization.'
      call messages_fatal(4, namespace=sys%namespace)
    end if
    if (((oct%algorithm == OPTION__OCTSCHEME__OCT_CG) .or. (oct%algorithm == OPTION__OCTSCHEME__OCT_BFGS)) &
      .and. target_move_ions(oct_target)) then
      write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_velocity", and'
      write(message(2), '(a)') '"OCTScheme = oct_cg", or "OCTScheme = oct_bfgs",'
      write(message(3), '(a)') 'then you have to set "MoveIons = false"'
      call messages_fatal(3, namespace=sys%namespace)
    end if
  end if

  if (target_type(oct_target) == oct_tg_hhgnew) then
    if ((oct%algorithm /= OPTION__OCTSCHEME__OCT_DIRECT) .and. &
      (oct%algorithm /= OPTION__OCTSCHEME__OCT_NLOPT_BOBYQA) .and. &
      (oct%algorithm /= OPTION__OCTSCHEME__OCT_BFGS) .and. &
      (oct%algorithm /= OPTION__OCTSCHEME__OCT_CG)) then
      write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_hhgnew", you can only use'
      write(message(2), '(a)') '"OCTScheme = oct_direct" or'
      write(message(3), '(a)') '"OCTScheme = oct_bobyqa" or'
      write(message(4), '(a)') '"OCTScheme = oct_cg" for the optimization.'
      call messages_fatal(4, namespace=sys%namespace)
    end if
    if (((oct%algorithm == OPTION__OCTSCHEME__OCT_CG) .or. (oct%algorithm == OPTION__OCTSCHEME__OCT_BFGS)) &
      .and. target_move_ions(oct_target)) then
      write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_hhgnew", and'
      write(message(2), '(a)') '"OCTScheme = oct_cg", or "OCTScheme = oct_bfgs",'
      write(message(3), '(a)') 'then you have to set "MoveIons = false"'
      call messages_fatal(3, namespace=sys%namespace)
    end if
  end if

  if (target_curr_functional(oct_target) /= oct_no_curr) then
    select case (sys%st%d%ispin)
    case (UNPOLARIZED)
    case (SPIN_POLARIZED)
      message(1) = 'Spin_polarized! Do not use OCT current functionals.'
      call messages_fatal(1, namespace=sys%namespace)
    case (SPINORS)
      message(1) = 'Spinors! Do not use OCT current functionals.'
      call messages_fatal(1, namespace=sys%namespace)
    end select
  end if

  select case (controlfunction_mode())
  case (controlfunction_mode_f)
    if (.not. oct_algorithm_is_direct(oct)) then
      if (oct%algorithm /= OPTION__OCTSCHEME__OCT_CG .and. oct%algorithm /= OPTION__OCTSCHEME__OCT_BFGS) then
        message(1) = 'If you attempt an envelope-only or phase-only optimization, then'
        message(2) = 'you must use either a gradient-free algorithm, oct_algorithm_cg, or'
        message(3) = 'oct_algorithm_bfgs algorithm.'
        call messages_fatal(3, namespace=sys%namespace)
      end if
    end if
  end select

  POP_SUB(check_faulty_runmodes)
end subroutine check_faulty_runmodes


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
