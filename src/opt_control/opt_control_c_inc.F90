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
#if defined(HAVE_NLOPT)
subroutine opt_control_nlopt_func(val, n, x, grad, need_gradient, f_data)
  real(c_double), intent(out) :: val
  integer(c_int), intent(in)  :: n
  real(c_double), intent(in)  :: x(*)
  real(c_double), intent(out) :: grad(*)
  integer(c_int), intent(in)  :: need_gradient
  type(c_ptr),    intent(in)  :: f_data

  integer :: getgrad
  PUSH_SUB(opt_control_nlopt_func)

  getgrad = 0
  if (need_gradient.ne.0) getgrad = 1
  call opt_control_cg_calc(n, x, val, getgrad, grad)

  POP_SUB(opt_control_nlopt_func)
end subroutine opt_control_nlopt_func
#endif

 !> ---------------------------------------------------------
 !! The following routines are to be called by C routines, which in turn
 !! are called by the main procedure of this module, opt_control_run, which
 !! is below.
 !! ---------------------------------------------------------
subroutine opt_control_function_forward(x, f)
  REAL_DOUBLE, intent(in)    :: x
  REAL_DOUBLE, intent(inout) :: f

  FLOAT, allocatable :: theta(:), y(:)
  integer :: dof
  type(opt_control_state_t) :: qcpsi

  PUSH_SUB(opt_control_function_forward)

  dof = controlfunction_dof(par_)
  SAFE_ALLOCATE(theta(1:dof))
  SAFE_ALLOCATE(y(1:dof))

  theta = x_
  theta(index_) = x

  call controlfunction_set_theta(par_, theta)
  call opt_control_state_null(qcpsi)
  call opt_control_state_copy(qcpsi, initial_st)
  call propagate_forward(sys_, td_, par_, oct_target, qcpsi)
  f = - target_j1(oct_target, sys_%namespace, sys_%gr, sys_%kpoints, qcpsi, sys_%ions) - controlfunction_j2(par_)

  SAFE_DEALLOCATE_A(theta)
  SAFE_DEALLOCATE_A(y)
  POP_SUB(opt_control_function_forward)
end subroutine opt_control_function_forward


 ! ---------------------------------------------------------
subroutine opt_control_cg_calc(n, x, f, getgrad, df)
  integer,        intent(in)    :: n
  REAL_DOUBLE,    intent(in)    :: x(n)
  REAL_DOUBLE,    intent(inout) :: f
  integer,        intent(in)    :: getgrad
  REAL_DOUBLE,    intent(inout) :: df(n)

  integer :: j
  type(controlfunction_t) :: par_new
  FLOAT :: j1, dx
  FLOAT, allocatable :: theta(:), abserr(:), dfn(:), dff(:)
  type(opt_control_state_t) :: qcpsi


  PUSH_SUB(opt_control_cg_calc)

  SAFE_ALLOCATE(theta(1:n))
  if (getgrad == 1) then
    theta = x
    call controlfunction_set_theta(par_, theta)
    call controlfunction_copy(par_new, par_)
    call f_striter(sys_, td_, par_new, j1)
    f = - j1 - controlfunction_j2(par_)
    call iteration_manager_direct(TOFLOAT(-f), par_, iterator, sys_)
    SAFE_ALLOCATE(dff(1:n))
    dff = df
    call controlfunction_gradient(par_, par_new, dff)
    df = dff

    ! Check if the gradient has been computed properly... This should be done only
    ! for debugging purposes.
    if (abs(oct%check_gradient) > M_ZERO) then
      dx = oct%check_gradient
      SAFE_ALLOCATE(dfn(1:n))
      SAFE_ALLOCATE(x_(1:n))
      SAFE_ALLOCATE(abserr(1:n))

      x_ = x
      do j = 1, n
        index_ = j
        call numder_ridders(x(j), dx, dfn(j), abserr(j), opt_control_function_forward)
      end do

      write(message(1), '(70(''#''))')
      write(message(2), *) &
        'GRADIENT (FORWARD-BACKWARD) |         GRADIENT (NUMERICAL)          |'
      call messages_info(2)
      do j = 1, n
        write(message(1), '(4x,es18.8,7x,a,3x,es18.8,a4,es8.1,6x,a)') &
          df(j), '|', dfn(j), ' +/-', abserr(j), '|'
        call messages_info(1)
      end do
      write(message(1), '(70(''-''))')
      write(message(2), '(a,es18.8,''                                        |'')') 'ABS DIFF = ', &
        norm2(df-dfn)
      write(message(3), '(a,es18.8,''                                        |'')') 'REL DIFF = ', &
        norm2(df-dfn)/norm2(dfn)
      write(message(4), '(70(''#''))')
      call messages_info(4)

      SAFE_DEALLOCATE_A(dfn)
      SAFE_DEALLOCATE_A(x_)
      SAFE_DEALLOCATE_A(abserr)
    end if

    call controlfunction_end(par_new)

  else
    theta = x
    call controlfunction_set_theta(par_, theta)
    call opt_control_state_null(qcpsi)
    call opt_control_state_copy(qcpsi, initial_st)
    call propagate_forward(sys_, td_, par_, oct_target, qcpsi)
    f = - target_j1(oct_target, sys_%namespace, sys_%gr, sys_%kpoints, qcpsi, sys_%ions) - controlfunction_j2(par_)
    call opt_control_state_end(qcpsi)
    call iteration_manager_direct(TOFLOAT(-f), par_, iterator, sys_)
  end if

  SAFE_DEALLOCATE_A(theta)
  POP_SUB(opt_control_cg_calc)
end subroutine opt_control_cg_calc
 ! ---------------------------------------------------------


 ! ---------------------------------------------------------
 !> interface is required by its being passed as dummy routine to minimize_multidim
subroutine opt_control_cg_write_info(iter, n, val, maxdx, maxdf, x)
  integer,     intent(in) :: iter, n
  REAL_DOUBLE, intent(in) :: val, maxdx, maxdf
  REAL_DOUBLE, intent(in) :: x(n)

  FLOAT :: fluence, j1, j2, j

  PUSH_SUB(opt_control_cg_write_info)

  j = - val
  fluence = controlfunction_fluence(par_)
  j2 = controlfunction_j2(par_)
  j1 = j - j2

  write(message(1), '(a,i5)') 'CG optimization iteration #', iter
  call messages_print_with_emphasis(msg=trim(message(1)), namespace=sys_%namespace)

  write(message(1), '(6x,a,f12.5)')    " => J1       = ", j1
  write(message(2), '(6x,a,f12.5)')    " => J        = ", j
  write(message(3), '(6x,a,f12.5)')    " => J2       = ", j2
  write(message(4), '(6x,a,f12.5)')    " => Fluence  = ", fluence
  write(message(5), '(6x,a,f12.5)')    " => Delta    = ", maxdx
  call messages_info(5)
  call messages_print_with_emphasis(namespace=sys_%namespace)

  call iteration_manager_main(iterator, j, j1, j2, TOFLOAT(maxdx))

  POP_SUB(opt_control_cg_write_info)
end subroutine opt_control_cg_write_info
 ! ---------------------------------------------------------


 ! ---------------------------------------------------------
 !> No intents here is unfortunately required because this will be passed to newuoa
 !! routines as a dummy function, whose interface has no intents.
subroutine opt_control_direct_calc(n, x, f)
  integer      :: n
  REAL_DOUBLE  :: x(n)
  REAL_DOUBLE  :: f

  FLOAT :: j1, delta
  FLOAT, allocatable :: theta(:)
  type(opt_control_state_t) :: qcpsi
  type(controlfunction_t) :: par_new

  PUSH_SUB(opt_control_direct_calc)

  SAFE_ALLOCATE(theta(1:n))
  theta = x
  call controlfunction_set_theta(par_, theta)

  if (abs(oct%delta) <= M_EPSILON) then
    ! We only need the value of the target functional.
    call opt_control_state_null(qcpsi)
    call opt_control_state_copy(qcpsi, initial_st)
    call propagate_forward(sys_, td_, par_, oct_target, qcpsi)
    f = - target_j1(oct_target, sys_%namespace, sys_%gr, sys_%kpoints, qcpsi, sys_%ions) - controlfunction_j2(par_)
    call opt_control_state_end(qcpsi)
    call iteration_manager_direct(TOFLOAT(-f), par_, iterator, sys_)
  else
    call controlfunction_copy(par_new, par_)
    call f_striter(sys_, td_, par_new, j1)
    delta = controlfunction_diff(par_, par_new)
    f = - oct%eta * j1 + oct%delta * delta
    call iteration_manager_direct(TOFLOAT(-f), par_, iterator, sys_, delta)
    call controlfunction_end(par_new)
  end if

  SAFE_DEALLOCATE_A(theta)
  POP_SUB(opt_control_direct_calc)
end subroutine opt_control_direct_calc
 ! ---------------------------------------------------------


 ! ---------------------------------------------------------
subroutine opt_control_direct_message_info(iter, n, val, maxdx, x)
  integer,     intent(in) :: iter, n
  REAL_DOUBLE, intent(in) :: val, maxdx
  REAL_DOUBLE, intent(in) :: x(n)

  FLOAT :: fluence, j1, j2, j
  FLOAT, allocatable :: theta(:)

  PUSH_SUB(opt_control_direct_message_info)

  SAFE_ALLOCATE(theta(1:n))
  theta = x
  call controlfunction_set_theta(par_, theta)
  SAFE_DEALLOCATE_A(theta)

  j = - val
  fluence = controlfunction_fluence(par_)
  j2 = controlfunction_j2(par_)
  j1 = j - j2

  write(message(1), '(a,i5)') 'Direct optimization iteration #', iter
  call messages_print_with_emphasis(msg=trim(message(1)), namespace=sys_%namespace)

  write(message(1), '(6x,a,f12.5)')    " => J1       = ", j1
  write(message(2), '(6x,a,f12.5)')    " => J        = ", j
  write(message(3), '(6x,a,f12.5)')    " => J2       = ", j2
  write(message(4), '(6x,a,f12.5)')    " => Fluence  = ", fluence
  write(message(5), '(6x,a,f12.5)')    " => Delta    = ", maxdx
  call messages_info(5)
  call messages_print_with_emphasis(namespace=sys_%namespace)

  call iteration_manager_main(iterator, j, j1, j2, TOFLOAT(maxdx))

  POP_SUB(opt_control_direct_message_info)
end subroutine opt_control_direct_message_info
 ! ---------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
