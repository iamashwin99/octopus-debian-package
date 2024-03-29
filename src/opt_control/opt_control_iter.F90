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

module opt_control_iter_oct_m
  use controlfunction_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use electrons_oct_m

  implicit none

  private
  public ::                   &
    oct_iterator_t,           &
    oct_iterator_init,        &
    oct_iterator_end,         &
    iteration_manager,        &
    iteration_manager_direct, &
    oct_iterator_bestpar,     &
    oct_iterator_current,     &
    oct_iterator_maxiter,     &
    oct_iterator_tolerance,   &
    iteration_manager_main,   &
    velocities_write


  type oct_iterator_t
    private
    FLOAT              :: eps
    integer            :: ctr_iter_max
    integer            :: ctr_iter
    integer            :: ctr_iter_main
    FLOAT              :: bestJ1
    FLOAT              :: bestJ1_fluence
    FLOAT              :: bestJ1_J
    integer            :: bestJ1_ctr_iter
    type(controlfunction_t) :: best_par
    integer            :: convergence_iunit
    integer            :: velocities_iunit
    logical            :: dump_intermediate
  end type oct_iterator_t

contains


  ! ---------------------------------------------------------
  subroutine oct_iterator_init(iterator, namespace, par)
    type(oct_iterator_t),      intent(inout) :: iterator
    type(namespace_t),         intent(in)    :: namespace
    type(controlfunction_t),   intent(in)    :: par

    PUSH_SUB(oct_iterator_init)

    !%Variable OCTEps
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 1.0e-6
    !%Description
    !% Define the convergence threshold. It computes the difference between the "input"
    !% field in the iterative procedure, and the "output" field. If this difference is
    !% less than <tt>OCTEps</tt> the iteration is stopped. This difference is defined as:
    !%
    !% <math>
    !% D[\varepsilon^{in},\varepsilon^{out}] = \int_0^T dt \left| \varepsilon^{in}(t)-\varepsilon^{out}(t)\right|^2
    !% </math>
    !%
    !% (If there are several control fields, this difference is defined as the sum over
    !% all the individual differences.)
    !%
    !% Whenever this condition is satisfied, it means that we have reached a solution point
    !% of the QOCT equations, <i>i.e.</i> a critical point of the QOCT functional (not
    !% necessarily a maximum, and not necessarily the global maximum).
    !%End
    call parse_variable(namespace, 'OCTEps', CNST(1.0e-6), iterator%eps)
    if (iterator%eps < M_ZERO) iterator%eps = tiny(M_ONE)

    !%Variable OCTMaxIter
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 10
    !%Description
    !% The maximum number of iterations.
    !% Typical values range from 10-100.
    !%End
    call parse_variable(namespace, 'OCTMaxIter', 10, iterator%ctr_iter_max)

    if (iterator%ctr_iter_max < 0 .and. iterator%eps < M_ZERO) then
      message(1) = "OCTMaxIter and OCTEps cannot be both < 0."
      call messages_fatal(1, namespace=namespace)
    end if
    if (iterator%ctr_iter_max < 0) iterator%ctr_iter_max = huge(iterator%ctr_iter_max)

    !%Variable OCTDumpIntermediate
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default true
    !%Description
    !% Writes to disk the laser pulse data during the OCT algorithm at intermediate steps.
    !% These are files called <tt>opt_control/laser.xxxx</tt>, where <tt>xxxx</tt> is the iteration number.
    !%End
    call parse_variable(namespace, 'OCTDumpIntermediate', .false., iterator%dump_intermediate)
    call messages_print_var_value("OCTDumpIntermediate", iterator%dump_intermediate, namespace=namespace)

    iterator%ctr_iter = 0
    iterator%ctr_iter_main = 0
    iterator%bestJ1          = -HUGE(iterator%bestJ1)
    iterator%bestJ1_fluence  = M_ZERO
    iterator%bestJ1_J        = M_ZERO
    iterator%bestJ1_ctr_iter = 0

    call controlfunction_copy(iterator%best_par, par)

    iterator%convergence_iunit = io_open(OCT_DIR//'convergence', namespace, action='write')

    write(iterator%convergence_iunit, '(91(''#''))')
    write(iterator%convergence_iunit, '(5(a))') '# iteration', '  J[Psi,chi,epsilon]', &
      '            J_1[Psi]', &
      '        J_2[epsilon]', &
      '               Delta'
    write(iterator%convergence_iunit, '(91(''#''))')

    if (parse_is_defined(namespace, 'OCTVelocityTarget')) then
      iterator%velocities_iunit = io_open(OCT_DIR//'velocities', namespace, action='write')
    end if

    POP_SUB(oct_iterator_init)
  end subroutine oct_iterator_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_iterator_end(iterator, namespace)
    type(oct_iterator_t), intent(inout) :: iterator
    type(namespace_t),    intent(in)    :: namespace

    PUSH_SUB(oct_iterator_end)

    call controlfunction_write(OCT_DIR//'laser.bestJ1', iterator%best_par, namespace)

    call controlfunction_end(iterator%best_par)
    write(iterator%convergence_iunit, '(91("#"))')
    call io_close(iterator%convergence_iunit)

    if (parse_is_defined(namespace, 'OCTVelocityTarget')) then
      call io_close(iterator%velocities_iunit)
    end if

    POP_SUB(oct_iterator_end)
  end subroutine oct_iterator_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical function iteration_manager(namespace, j1, par, par_prev, iterator) result(stoploop)
    type(namespace_t),       intent(in)    :: namespace
    FLOAT,                   intent(in)    :: j1
    type(controlfunction_t), intent(in)    :: par
    type(controlfunction_t), intent(in)    :: par_prev
    type(oct_iterator_t),    intent(inout) :: iterator

    FLOAT :: fluence, jfunctional, j2, delta
    logical :: bestj1

    PUSH_SUB(iteration_manager)

    if (iterator%dump_intermediate) call iterator_write(namespace, iterator, par)

    stoploop = .false.

    fluence = controlfunction_fluence(par)
    j2 = controlfunction_j2(par)
    jfunctional = j1 + j2
    delta = controlfunction_diff(par, par_prev)

    if (iterator%ctr_iter == iterator%ctr_iter_max) then
      message(1) = "Info: Maximum number of iterations reached."
      call messages_info(1, namespace=namespace)
      stoploop = .true.
    end if

    if ((iterator%eps > M_ZERO) .and. &
      (delta < iterator%eps) .and. &
      (iterator%ctr_iter > 0)) then
      message(1) = "Info: Convergence threshold reached."
      call messages_info(1, namespace=namespace)
      stoploop = .true.
    end if

    write(message(1), '(a,i5)') 'Optimal control iteration #', iterator%ctr_iter
    call messages_print_with_emphasis(msg=trim(message(1)), namespace=namespace)

    write(message(1), '(6x,a,f12.5)')  " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)')  " => J        = ", jfunctional
    write(message(3), '(6x,a,f12.5)')  " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)')  " => Fluence  = ", fluence
    write(message(5), '(6x,a,f12.5)')  " => Penalty  = ", controlfunction_alpha(par, 1)
    write(message(6), '(6x,a,es12.2)') " => D[e,e']  = ", delta
    if (iterator%ctr_iter /= 0) then
      call messages_info(6, namespace=namespace)
    else
      call messages_info(5, namespace=namespace)
    end if
    call messages_print_with_emphasis(namespace=namespace)

    bestj1 = (j1 > iterator%bestj1)
    ! store field with best J1
    if (bestj1) then
      iterator%bestJ1          = j1
      iterator%bestJ1_J        = jfunctional
      iterator%bestJ1_fluence  = fluence
      iterator%bestJ1_ctr_iter = iterator%ctr_iter
      call controlfunction_end(iterator%best_par)
      call controlfunction_copy(iterator%best_par, par)
      if (.not. iterator%dump_intermediate) call controlfunction_write(OCT_DIR//'laser.bestJ1', &
        iterator%best_par, namespace)
    end if

    write(iterator%convergence_iunit, '(i11,4f20.8)')                &
      iterator%ctr_iter, jfunctional, j1, j2, controlfunction_diff(par, par_prev)

    iterator%ctr_iter = iterator%ctr_iter + 1

    POP_SUB(iteration_manager)
  end function iteration_manager
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine iteration_manager_direct(j, par, iterator, sys, dx)
    FLOAT,                   intent(in)    :: j
    type(controlfunction_t), intent(in)    :: par
    type(oct_iterator_t),    intent(inout) :: iterator
    type(electrons_t),       intent(in)    :: sys
    FLOAT,         optional, intent(in)    :: dx

    FLOAT :: j1, j2, fluence, delta
    PUSH_SUB(iteration_manager_direct)

    if (iterator%dump_intermediate) call iterator_write(sys%namespace, iterator, par)

    fluence = controlfunction_fluence(par)
    j2 = controlfunction_j2(par)
    j1 = j - j2

    if (present(dx)) then
      delta = dx
    else
      delta = M_ZERO
    end if

    if (iterator%ctr_iter == 0) then
      call messages_print_with_emphasis(msg='Initial-guess field', namespace=sys%namespace)
    else
      write(message(1), '(a,i5)') 'Function evaluation #', iterator%ctr_iter
      call messages_print_with_emphasis(msg=trim(message(1)), namespace=sys%namespace)
    end if

    write(message(1), '(6x,a,f12.5)')    " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)')    " => J        = ", j
    write(message(3), '(6x,a,f12.5)')    " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)')    " => Fluence  = ", fluence
    call messages_info(4, namespace=sys%namespace)
    if (present(dx)) then
      write(message(1), '(6x,a,f12.5)')  " => Delta    = ", dx
      call messages_info(1, namespace=sys%namespace)
    end if
    call messages_print_with_emphasis(namespace=sys%namespace)

    ! store field with best J1
    if (j1 > iterator%bestJ1) then
      iterator%bestJ1          = j1
      iterator%bestJ1_J        = j
      iterator%bestJ1_fluence  = fluence
      iterator%bestJ1_ctr_iter = iterator%ctr_iter
      call controlfunction_end(iterator%best_par)
      call controlfunction_copy(iterator%best_par, par)
      if (.not. iterator%dump_intermediate) call controlfunction_write(OCT_DIR//'laser.bestJ1', &
        iterator%best_par, sys%namespace)
    end if

    write(iterator%convergence_iunit, '(i11,4f20.8)')                &
      iterator%ctr_iter, j, j1, j2, delta

    if (parse_is_defined(sys%namespace, 'OCTVelocityTarget')) then
      call velocities_write(iterator, sys)
    end if

    iterator%ctr_iter = iterator%ctr_iter + 1

    POP_SUB(iteration_manager_direct)
  end subroutine iteration_manager_direct
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine iteration_manager_main(iterator, j, j1, j2, delta)
    type(oct_iterator_t), intent(inout) :: iterator

    FLOAT, intent(in) :: j, j1, j2, delta

    PUSH_SUB(iteration_manager_main)

    iterator%ctr_iter_main = iterator%ctr_iter_main + 1
    write(iterator%convergence_iunit, '("### MAIN ITERATION")')
    write(iterator%convergence_iunit, '(a2,i9,4f20.8)')                &
      '##', iterator%ctr_iter_main, j, j1, j2, delta
    write(iterator%convergence_iunit, '("###")')

    POP_SUB(iteration_manager_main)
  end subroutine iteration_manager_main
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine iterator_write(namespace, iterator, par)
    type(namespace_t),              intent(in) :: namespace
    type(oct_iterator_t),           intent(in) :: iterator
    type(controlfunction_t),        intent(in) :: par

    character(len=80)  :: filename

    PUSH_SUB(iterator_write)

    write(filename,'(a,i4.4)') OCT_DIR//'laser.', iterator%ctr_iter
    call controlfunction_write(filename, par, namespace)

    POP_SUB(iterator_write)
  end subroutine iterator_write
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_iterator_bestpar(par, iterator)
    type(oct_iterator_t), target, intent(inout) :: iterator
    type(controlfunction_t), pointer    :: par

    PUSH_SUB(oct_iterator_bestpar)
    par => iterator%best_par

    POP_SUB(oct_iterator_bestpar)
  end subroutine oct_iterator_bestpar
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function oct_iterator_current(iterator)
    type(oct_iterator_t), intent(in)     :: iterator
    oct_iterator_current = iterator%ctr_iter
  end function oct_iterator_current
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function oct_iterator_maxiter(iterator)
    type(oct_iterator_t), intent(in)     :: iterator
    oct_iterator_maxiter = iterator%ctr_iter_max
  end function oct_iterator_maxiter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function oct_iterator_tolerance(iterator)
    type(oct_iterator_t), intent(in)     :: iterator
    oct_iterator_tolerance = iterator%eps
  end function oct_iterator_tolerance
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine velocities_write(iterator, sys)
    type(oct_iterator_t),    intent(in) :: iterator
    type(electrons_t),       intent(in) :: sys

    character (len=100) :: temp_str
    character (len=2) :: atoms_str
    character (len=1) :: dim_str
    integer :: i, j, n_atoms, dim

    PUSH_SUB(velocities_write)

    n_atoms = sys%ions%natoms
    dim = sys%space%dim

    ! write header of the velocities output file
    if (iterator%ctr_iter == 0) then
      write(iterator%velocities_iunit,'(100("#"))')
      write(iterator%velocities_iunit,'("#  iter")',advance='no')
      do i = 1, n_atoms
        write(atoms_str,'(i2.2)') i
        do j = 1, dim
          write(dim_str,'(i1)') j
          temp_str = "v[" // atoms_str // "," // dim_str // "]"
          write(iterator%velocities_iunit,'(a16)',advance='no') trim(temp_str)
        end do
      end do
      write(iterator%velocities_iunit,'("")')
      write(iterator%velocities_iunit,'(100("#"))')
    end if

    ! write data
    write(iterator%velocities_iunit,'(i7)',advance='no') iterator%ctr_iter
    do i = 1, n_atoms
      do j = 1, dim
        write(iterator%velocities_iunit,'(4(" "),(f12.10))',advance='no') sys%ions%vel(j, i)
      end do
    end do
    write(iterator%velocities_iunit,'("")')

    POP_SUB(velocities_write)
  end subroutine velocities_write
  ! ---------------------------------------------------------


end module opt_control_iter_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
