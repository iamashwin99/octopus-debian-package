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
!! E.L. Briggs, D.J. Sullivan, and J. Bernholc, PRB 54 14362 (1996)
!!
!! It assumes that the Oxygen atom is located at x0=0 (see Eq. (12))

module curv_briggs_oct_m
  use coordinate_system_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                     &
    curv_briggs_t,              &
    curv_briggs_copy

  type, extends(coordinate_system_t) :: curv_briggs_t
    private
    FLOAT, allocatable :: lsize(:) !< size of the box
    FLOAT :: beta                  !< adjustable parameter between 0 and 1 that controls the degree of scaling
  contains
    procedure :: to_cartesian => curv_briggs_to_cartesian
    procedure :: from_cartesian => curv_briggs_from_cartesian
    procedure :: det_jac => curv_briggs_det_jac
    procedure :: write_info => curv_briggs_write_info
    procedure :: surface_element => curv_briggs_surface_element
    final :: curv_briggs_finalize
  end type curv_briggs_t

  interface curv_briggs_t
    procedure curv_briggs_constructor
  end interface curv_briggs_t

contains

  ! ---------------------------------------------------------
  function curv_briggs_constructor(namespace, dim, lsize, spacing) result(briggs)
    type(namespace_t),   intent(in)  :: namespace
    integer,             intent(in)  :: dim
    FLOAT,               intent(in)  :: lsize(1:dim)
    FLOAT,               intent(in)  :: spacing(1:dim)
    class(curv_briggs_t), pointer :: briggs

    integer :: idim

    PUSH_SUB(curv_briggs_constructor)

    SAFE_ALLOCATE(briggs)

    briggs%local_basis = .true.
    briggs%orthogonal = .true. ! This needs to be checked
    briggs%dim = dim
    SAFE_ALLOCATE(briggs%lsize(1:dim))
    briggs%lsize(1:dim) = lsize(1:dim)

    call parse_variable(namespace, 'CurvBriggsBeta', M_HALF, briggs%beta)

    if (briggs%beta < M_ZERO .or. briggs%beta > M_ONE) then
      message(1) = 'The parameter "CurvBriggsBeta" must lie between 0 and 1.'
      call messages_fatal(1, namespace=namespace)
    end if

    briggs%min_mesh_scaling_product = M_ONE
    do idim = 1, briggs%dim
      ! corresponds to the distance of grid points at [+spacing/2,-spacing/2]
      briggs%min_mesh_scaling_product = briggs%min_mesh_scaling_product * (M_ONE / &
        (M_ONE - briggs%lsize(idim) * briggs%beta / (M_PI * spacing(idim)) * sin(M_PI * spacing(idim) / briggs%lsize(idim))))
    end do

    POP_SUB(curv_briggs_constructor)
  end function curv_briggs_constructor

  ! ---------------------------------------------------------
  subroutine curv_briggs_finalize(this)
    type(curv_briggs_t), intent(inout) :: this

    PUSH_SUB(curv_briggs_finalize)

    SAFE_DEALLOCATE_A(this%lsize)

    POP_SUB(curv_briggs_finalize)
  end subroutine curv_briggs_finalize

  ! ---------------------------------------------------------
  subroutine curv_briggs_copy(this_out, this_in)
    type(curv_briggs_t), intent(inout) :: this_out
    type(curv_briggs_t), intent(in)    :: this_in

    PUSH_SUB(curv_briggs_copy)

    SAFE_ALLOCATE_SOURCE_A(this_out%lsize, this_in%lsize)
    this_out%beta = this_in%beta

    POP_SUB(curv_briggs_copy)
  end subroutine curv_briggs_copy

  ! ---------------------------------------------------------
  function curv_briggs_to_cartesian(this, chi) result(xx)
    class(curv_briggs_t), target, intent(in)  :: this
    FLOAT,                        intent(in)  :: chi(:)  !< chi(dim)
    FLOAT :: xx(1:this%dim)

    ! no PUSH_SUB, called too often

    xx = chi - this%lsize*this%beta/(M_TWO*M_PI)*sin(M_TWO*M_PI*chi/this%lsize)

  end function curv_briggs_to_cartesian

  ! ---------------------------------------------------------
  function curv_briggs_from_cartesian(this, xx) result(chi)
    class(curv_briggs_t), target, intent(in)  :: this
    FLOAT,                        intent(in)  :: xx(:)   !< xx(dim)
    FLOAT :: chi(1:this%dim)

    ! no PUSH_SUB, called too often

    chi = M_ZERO
    message(1) = "Internal error in curv_briggs_from_cartesian"
    call messages_fatal(1)

  end function curv_briggs_from_cartesian

  ! ---------------------------------------------------------
  FLOAT function curv_briggs_det_jac(this, xx, chi) result(jdet)
    class(curv_briggs_t), intent(in)  :: this
    FLOAT,                intent(in)  :: xx(:)
    FLOAT,                intent(in)  :: chi(:)

    integer :: idim
    FLOAT :: jac(1:this%dim)

    ! no PUSH_SUB, called too often

    ! Jacobian is diagonal in this method
    do idim = 1, this%dim
      jac(idim) = M_ONE - this%beta*cos(M_TWO*M_PI*chi(idim)/this%lsize(idim))
    end do
    jdet = product(jac)

  end function curv_briggs_det_jac

  ! ---------------------------------------------------------
  subroutine curv_briggs_write_info(this, iunit, namespace)
    class(curv_briggs_t),         intent(in) :: this
    integer,            optional, intent(in) :: iunit
    type(namespace_t),  optional, intent(in) :: namespace

    PUSH_SUB(curv_briggs_write_info)

    write(message(1), '(a)') '  Curvilinear Method = briggs'
    call messages_info(1, iunit=iunit, namespace=namespace)

    POP_SUB(curv_briggs_write_info)
  end subroutine curv_briggs_write_info

  ! ---------------------------------------------------------
  FLOAT function curv_briggs_surface_element(this, idir) result(ds)
    class(curv_briggs_t), intent(in) :: this
    integer,              intent(in) :: idir

    ds = M_ZERO
    message(1) = 'Surface element with briggs curvilinear coordinates not implemented'
    call messages_fatal(1)

  end function curv_briggs_surface_element

end module curv_briggs_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
