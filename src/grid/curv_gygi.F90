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

!> This module implements the curvilinear coordinates given in
!! F. Gygi and G. Galli, PRB 52 R2229 (1996).

module curv_gygi_oct_m
  use coordinate_system_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_adv_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use root_solver_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                   &
    curv_gygi_t,              &
    curv_gygi_copy

  type, extends(coordinate_system_t) :: curv_gygi_t
    private
    FLOAT, public :: A             !< local reduction in grid spacing is 1/(1+A)
    FLOAT, public :: alpha         !< range of enhancement of the resolution
    FLOAT, public :: beta          !< distance over which Euclidian coordinates are recovered
    FLOAT, allocatable :: pos(:, :)
    integer :: npos
    type(root_solver_t) :: rs
  contains
    procedure :: to_cartesian => curv_gygi_to_cartesian
    procedure :: from_cartesian => curv_gygi_from_cartesian
    procedure :: det_jac => curv_gygi_det_jac
    procedure :: write_info => curv_gygi_write_info
    procedure :: surface_element => curv_gygi_surface_element
    final :: curv_gygi_finalize
  end type curv_gygi_t

  interface curv_gygi_t
    procedure curv_gygi_constructor
  end interface curv_gygi_t

  ! Auxiliary variables for the root solver.
  class(curv_gygi_t), pointer  :: gygi_p
  integer :: i_p
  FLOAT, allocatable :: chi_p(:)

contains

  ! ---------------------------------------------------------
  function curv_gygi_constructor(namespace, dim, npos, pos) result(gygi)
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(in)  :: dim
    integer,           intent(in)  :: npos
    FLOAT,             intent(in)  :: pos(1:dim,1:npos)
    class(curv_gygi_t), pointer :: gygi

    PUSH_SUB(curv_gygi_constructor)

    SAFE_ALLOCATE(gygi)

    gygi%dim = dim
    gygi%local_basis = .true.
    gygi%orthogonal = .true. ! This needs to be checked.

    gygi%npos = npos
    SAFE_ALLOCATE(gygi%pos(1:dim, 1:gygi%npos))
    gygi%pos = pos

    !%Variable CurvGygiA
    !%Type float
    !%Default 0.5
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% The grid spacing is reduced locally around each atom, and the reduction is
    !% given by 1/(1+<i>A</i>), where <i>A</i> is specified by this variable. So, if
    !% <i>A</i>=1/2 (the default), the grid spacing is reduced to two thirds = 1/(1+1/2).
    !% [This is the <math>A_{\alpha}</math> variable in Eq. 2 of F. Gygi and G. Galli, <i>Phys.
    !% Rev. B</i> <b>52</b>, R2229 (1995)]. It must be larger than zero.
    !%End
    call parse_variable(namespace, 'CurvGygiA', M_HALF, gygi%A)

    !%Variable CurvGygiAlpha
    !%Type float
    !%Default 2.0 a.u.
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% This number determines the region over which the grid is enhanced (range of
    !% enhancement of the resolution). That is, the grid is enhanced on a sphere
    !% around each atom, whose radius is given by this variable. [This is the <math>a_{\alpha}</math>
    !% variable in Eq. 2 of F. Gygi and G. Galli, <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)].
    !% It must be larger than zero.
    !%End

    call parse_variable(namespace, 'CurvGygiAlpha', M_TWO, gygi%alpha, units_inp%length)
    !%Variable CurvGygiBeta
    !%Type float
    !%Default 4.0 a.u.
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% This number determines the distance over which Euclidean coordinates are
    !% recovered. [This is the <math>b_{\alpha}</math> variable in Eq. 2 of F. Gygi and G. Galli,
    !% <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)]. It must be larger than zero.
    !%End
    call parse_variable(namespace, 'CurvGygiBeta', M_FOUR, gygi%beta, units_inp%length)

    if (gygi%a <= M_ZERO)     call messages_input_error(namespace, 'CurvGygiA')
    if (gygi%alpha <= M_ZERO) call messages_input_error(namespace, 'CurvGygiAlpha')
    if (gygi%beta <= M_ZERO)  call messages_input_error(namespace, 'CurvGygiBeta')

    gygi%min_mesh_scaling_product = (M_ONE / (M_ONE + gygi%A))**gygi%dim

    ! initialize root solver
    call root_solver_init(gygi%rs, namespace, dim, solver_type = ROOT_NEWTON, maxiter = 500, abs_tolerance = CNST(1.0e-10))

    POP_SUB(curv_gygi_constructor)
  end function curv_gygi_constructor

  ! ---------------------------------------------------------
  subroutine curv_gygi_copy(this_out, this_in)
    type(curv_gygi_t), intent(inout) :: this_out
    type(curv_gygi_t), intent(in)    :: this_in

    PUSH_SUB(curv_gygi_copy)

    this_out%A = this_in%A
    this_out%alpha = this_in%alpha
    this_out%beta = this_in%beta
    SAFE_ALLOCATE_SOURCE_A(this_out%pos, this_in%pos)
    this_out%npos = this_in%npos

    POP_SUB(curv_gygi_copy)
  end subroutine curv_gygi_copy

  ! ---------------------------------------------------------
  subroutine curv_gygi_finalize(this)
    type(curv_gygi_t), intent(inout) :: this

    PUSH_SUB(curv_gygi_finalize)

    SAFE_DEALLOCATE_A(this%pos)

    POP_SUB(curv_gygi_finalize)
  end subroutine curv_gygi_finalize

  ! ---------------------------------------------------------
  function curv_gygi_to_cartesian(this, chi) result(xx)
    class(curv_gygi_t), target, intent(in)  :: this
    FLOAT,                      intent(in)  :: chi(:)  !< chi(dim)
    FLOAT :: xx(1:this%dim)

    integer :: i
    logical :: conv

    ! no PUSH_SUB, called too often

    gygi_p => this
    i_p =  this%npos
    SAFE_ALLOCATE(chi_p(1:this%dim))
    chi_p(:) = chi(:)

    call droot_solver_run(this%rs, getf, xx, conv, startval = chi)

    if (.not. conv) then
      do i = 1, this%npos
        conv = .false.
        i_p = i
        call droot_solver_run(this%rs, getf, xx, conv, startval = xx(1:this%dim))
      end do
    end if

    nullify(gygi_p)
    SAFE_DEALLOCATE_A(chi_p)

    if (.not. conv) then
      message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
      message(2) = "method did not converge for point:"
      write(message(3),'(9f14.6)') xx(1:this%dim)
      message(4) = "Try varying the Gygi parameters -- usually reducing CurvGygiA or"
      message(5) = "CurvGygiAlpha (or both) solves the problem."
      call messages_fatal(5)
    end if

  end function curv_gygi_to_cartesian

  ! ---------------------------------------------------------
  pure function curv_gygi_from_cartesian(this, xx) result(chi)
    class(curv_gygi_t), target, intent(in)  :: this
    FLOAT,                      intent(in)  :: xx(:)    !< xx(dim)
    FLOAT :: chi(1:this%dim)

    integer :: i, ia
    FLOAT   :: r, ar, th, ex

    ! no PUSH_SUB, called too often

    chi(1:this%dim) = xx(1:this%dim)
    do ia = 1, this%npos
      r = max(norm2(xx(1:this%dim) - this%pos(1:this%dim, ia)), CNST(1e-6))
      ar = this%A*this%alpha/r
      th = tanh(r/this%alpha)
      ex = exp(-(r/this%beta)**2)
      do i = 1, this%dim
        chi(i) = chi(i) + (xx(i) - this%pos(i, ia)) * this%a * ar * th * ex
      end do
    end do

  end function curv_gygi_from_cartesian

  ! ---------------------------------------------------------
  FLOAT function curv_gygi_det_jac(this, xx, chi) result(jdet)
    class(curv_gygi_t), intent(in)  :: this
    FLOAT,              intent(in)  :: xx(:)   !< xx(dim)
    FLOAT,              intent(in)  :: chi(:)  !< chi(dim)

    FLOAT :: dummy(this%dim)
    FLOAT :: jac(1:this%dim, 1:this%dim)

    ! no PUSH_SUB, called too often

    call curv_gygi_jacobian(this, xx, dummy, jac)
    jdet = M_ONE/lalg_determinant(this%dim, jac, preserve_mat = .false.)

  end function curv_gygi_det_jac

  ! ---------------------------------------------------------
  subroutine curv_gygi_write_info(this, iunit, namespace)
    class(curv_gygi_t),           intent(in) :: this
    integer,            optional, intent(in) :: iunit
    type(namespace_t),  optional, intent(in) :: namespace

    PUSH_SUB(curv_gygi_write_info)

    write(message(1), '(a)')  '  Curvilinear Method = gygi'
    write(message(2), '(a)')  '  Gygi Parameters:'
    write(message(3), '(4x,a,f6.3)')  'A = ', this%a
    write(message(4), '(4x,3a,f6.3)') 'alpha [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%alpha)
    write(message(5), '(4x,3a,f6.3)') 'beta  [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%beta)
    call messages_info(5, iunit=iunit, namespace=namespace)

    POP_SUB(curv_gygi_write_info)
  end subroutine curv_gygi_write_info

  ! ---------------------------------------------------------
  FLOAT function curv_gygi_surface_element(this, idir) result(ds)
    class(curv_gygi_t), intent(in) :: this
    integer,            intent(in) :: idir

    ds = M_ZERO
    message(1) = 'Surface element with gygi curvilinear coordinates not implemented'
    call messages_fatal(1)

  end function curv_gygi_surface_element

  ! ---------------------------------------------------------
  pure subroutine curv_gygi_jacobian(this, xx, chi, jac, natoms)
    class(curv_gygi_t), intent(in)  :: this
    FLOAT,              intent(in)  :: xx(:)       !< x(dim)
    FLOAT,              intent(out) :: chi(:)     !< chi(dim)
    FLOAT,              intent(out) :: jac(:, :)  !< jac(dim,dim), the Jacobian
    integer,  optional, intent(in)  :: natoms

    integer :: i, ix, iy, natoms_
    FLOAT :: r, f_alpha, df_alpha
    FLOAT :: th, ex, ar

    ! no PUSH_SUB, called too often

    jac(1:this%dim, 1:this%dim) = diagonal_matrix(this%dim, M_ONE)
    chi(1:this%dim) = xx(1:this%dim)

    natoms_ = this%npos
    if (present(natoms)) natoms_ = natoms

    do i = 1, natoms_
      r = max(norm2(xx(1:this%dim) - this%pos(1:this%dim, i)), CNST(1e-6))

      ar = this%A*this%alpha/r
      th = tanh(r/this%alpha)
      ex = exp(-(r/this%beta)**2)

      f_alpha  = ar * th * ex
      df_alpha = ar*(-th*ex/r + ex/(this%alpha*cosh(r/this%alpha)**2) - th*M_TWO*r*ex/this%beta**2)

      do ix = 1, this%dim
        chi(ix) = chi(ix) + f_alpha*(xx(ix) - this%pos(ix, i))

        jac(ix, ix) = jac(ix, ix) + f_alpha
        do iy = 1, this%dim
          jac(ix, iy) = jac(ix, iy) + (xx(ix) - this%pos(ix, i))*(xx(iy) - this%pos(iy, i))/r*df_alpha
        end do
      end do
    end do

  end subroutine curv_gygi_jacobian

  ! ---------------------------------------------------------
  pure subroutine getf(y, f, jf)
    FLOAT, intent(in)    :: y(:)
    FLOAT, intent(out)   :: f(:), jf(:, :)

    ! no PUSH_SUB, called too often

    call curv_gygi_jacobian(gygi_p, y, f, jf, i_p)
    f(1:gygi_p%dim) = f(1:gygi_p%dim) - chi_p(1:gygi_p%dim)

  end subroutine getf

end module curv_gygi_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
