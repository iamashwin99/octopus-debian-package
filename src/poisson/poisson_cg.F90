!! Cropyright (C) 2004-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module poisson_cg_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use poisson_corrections_oct_m
  use profiling_oct_m
  use solvers_oct_m

  implicit none

  private
  public ::           &
    poisson_cg_init,  &
    poisson_cg_end,   &
    poisson_cg1,      &
    poisson_cg2

  FLOAT, public :: threshold
  integer, public :: maxiter

contains


  ! ---------------------------------------------------------
  subroutine poisson_cg_init(thr, itr)
    integer,      intent(in) :: itr
    FLOAT,        intent(in) :: thr

    PUSH_SUB(poisson_cg_init)
    threshold = thr
    maxiter = itr
    POP_SUB(poisson_cg_init)
  end subroutine poisson_cg_init


  ! ---------------------------------------------------------
  subroutine poisson_cg_end

  end subroutine poisson_cg_end


  ! ---------------------------------------------------------
  subroutine poisson_cg1(namespace, der, corrector, pot, rho)
    type(namespace_t),            intent(in)    :: namespace
    type(derivatives_t),  target, intent(in)    :: der
    type(poisson_corr_t),         intent(in)    :: corrector
    FLOAT,                        intent(inout) :: pot(:)
    FLOAT,                        intent(in)    :: rho(:)

    integer :: iter
    FLOAT :: res
    FLOAT, allocatable :: wk(:), lwk(:), zk(:), pk(:)

    PUSH_SUB(poisson_cg1)

    SAFE_ALLOCATE( wk(1:der%mesh%np_part))
    SAFE_ALLOCATE(lwk(1:der%mesh%np_part))
    SAFE_ALLOCATE( zk(1:der%mesh%np_part))
    SAFE_ALLOCATE( pk(1:der%mesh%np_part))

    ! build initial guess for the potential
    wk(1:der%mesh%np) = pot(1:der%mesh%np)
    call poisson_boundary_conditions(corrector, der%mesh, rho, wk)
    call dderivatives_lapl(der, wk, lwk, ghost_update = .true., set_bc = .false.)

    zk(1:der%mesh%np) = -M_FOUR*M_PI*rho(1:der%mesh%np) - lwk(1:der%mesh%np)
    SAFE_DEALLOCATE_A(wk)
    SAFE_DEALLOCATE_A(lwk) ! they are no longer needed

    der_pointer  => der
    mesh_pointer => der%mesh
    call lalg_copy(der%mesh%np, zk, pk)
    pk(der%mesh%np + 1:) = M_ZERO
    iter = maxiter
    call dconjugate_gradients(der%mesh%np, pk, zk, internal_laplacian_op, internal_dotp, iter, res, threshold)
    if (res >= threshold) then
      message(1) = 'Conjugate-gradients Poisson solver did not converge.'
      write(message(2), '(a,i8)')    '  Iter = ',iter
      write(message(3), '(a,e14.6)') '  Res = ', res
      call messages_warning(3, namespace=namespace)
    end if
    nullify(der_pointer, mesh_pointer)
    pot(1:der%mesh%np) = pot(1:der%mesh%np) + pk(1:der%mesh%np)

    SAFE_DEALLOCATE_A(zk)
    SAFE_DEALLOCATE_A(pk)
    POP_SUB(poisson_cg1)
  end subroutine poisson_cg1


  ! ---------------------------------------------------------
  subroutine poisson_cg2(namespace, der, pot, rho)
    type(namespace_t),           intent(in)    :: namespace
    type(derivatives_t), target, intent(in)    :: der
    FLOAT,                       intent(inout) :: pot(:)
    FLOAT,                       intent(in)    :: rho(:)

    integer :: iter, ip
    FLOAT, allocatable :: potc(:), rhs(:)
    FLOAT :: res

    PUSH_SUB(poisson_cg2)

    iter = maxiter
    der_pointer  => der
    mesh_pointer => der%mesh

    SAFE_ALLOCATE(rhs(1:der%mesh%np))
    SAFE_ALLOCATE(potc(1:der%mesh%np_part))

    do ip = 1, der%mesh%np
      rhs(ip) = CNST(-4.0)*M_PI*rho(ip)
    end do
    call lalg_copy(der%mesh%np, pot, potc)

    call dconjugate_gradients(der%mesh%np, potc, rhs, internal_laplacian_op, internal_dotp, iter, res, threshold)

    if (res >= threshold) then
      message(1) = 'Conjugate-gradients Poisson solver did not converge.'
      write(message(2), '(a,i8)')    '  Iter = ', iter
      write(message(3), '(a,e14.6)') '  Res = ', res
      call messages_warning(3, namespace=namespace)
    end if

    call lalg_copy(der%mesh%np, potc, pot)

    nullify(der_pointer, mesh_pointer)

    SAFE_DEALLOCATE_A(rhs)
    SAFE_DEALLOCATE_A(potc)

    POP_SUB(poisson_cg2)
  end subroutine poisson_cg2

end module poisson_cg_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
