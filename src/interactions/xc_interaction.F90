!! Copyright (C) 2021 Nicolas Tancogne-Dejean
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

module xc_interaction_oct_m
  use comm_oct_m
  use debug_oct_m
  use distributed_oct_m
  use global_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use density_interaction_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use space_oct_m
  use states_elec_dim_oct_m
  use xc_f03_lib_m
  use xc_functional_oct_m

  implicit none

  private
  public ::                &
    xc_interaction_t,      &
    xc_interaction_compute,&
    calc_tb09_c,           &
    calc_mvorb_alpha 

  type, extends(density_interaction_t) :: xc_interaction_t
    private

  contains
    procedure :: init => xc_interaction_init
    procedure :: calculate => xc_interaction_calculate
    procedure :: calculate_energy => xc_interaction_calculate_energy
    procedure :: end => xc_interaction_end
    final :: xc_interaction_finalize
  end type xc_interaction_t

  interface xc_interaction_t
    module procedure xc_interaction_constructor
  end interface xc_interaction_t

  integer, public, parameter :: &
    FUNC_X = 1,         &
    FUNC_C = 2


contains

  ! ---------------------------------------------------------
  function xc_interaction_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(xc_interaction_t),               pointer       :: this

    PUSH_SUB(xc_interaction_constructor)

    SAFE_ALLOCATE(this)

    this%label = "exchange-correlation interaction"

    this%partner => partner

    !We do not need any system quantity here
    this%n_system_quantities = 0

    POP_SUB(xc_interaction_constructor)
  end function xc_interaction_constructor

  ! ---------------------------------------------------------
  subroutine xc_interaction_init(this)
    class(xc_interaction_t),    intent(inout) :: this

    PUSH_SUB(xc_interaction_init)

    POP_SUB(xc_interaction_init)
  end subroutine xc_interaction_init

  ! ---------------------------------------------------------
  subroutine xc_interaction_calculate(this)
    class(xc_interaction_t),             intent(inout) :: this

    PUSH_SUB(xc_interaction_calculate)

    this%density = M_ZERO

    POP_SUB(xc_interaction_calculate)
  end subroutine xc_interaction_calculate

  ! ---------------------------------------------------------
  subroutine xc_interaction_calculate_energy(this)
    class(xc_interaction_t),             intent(inout) :: this

    PUSH_SUB(xc_interaction_calculate_energy)

    this%energy = M_ZERO

    POP_SUB(xc_interaction_calculate_energy)
  end subroutine xc_interaction_calculate_energy


  ! ---------------------------------------------------------
  subroutine xc_interaction_compute(this) 
    class(xc_interaction_t),             intent(inout) :: this
    PUSH_SUB(xc_interaction_compute)

    this%density = M_ZERO

    POP_SUB(xc_interaction_compute)
  end subroutine xc_interaction_compute

  ! ---------------------------------------------------------
  subroutine xc_interaction_end(this)
    class(xc_interaction_t), intent(inout) :: this

    PUSH_SUB(xc_interaction_end)

    SAFE_DEALLOCATE_A(this%density)

    call interaction_with_partner_end(this)

    POP_SUB(xc_interaction_end)
  end subroutine xc_interaction_end


  ! ---------------------------------------------------------
  subroutine xc_interaction_finalize(this)
    type(xc_interaction_t), intent(inout) :: this

    PUSH_SUB(xc_interaction_finalize)

    call xc_interaction_end(this)

    POP_SUB(xc_interaction_finalize)
  end subroutine xc_interaction_finalize

  ! -----------------------------------------------------
  subroutine calc_tb09_c(mesh, space, functl, dens, gdens, ispin, rcell_volume)
    type(mesh_t),          intent(in) :: mesh
    type(space_t),         intent(in) :: space
    type(xc_functional_t), intent(inout) :: functl(:)
    FLOAT,                 intent(in) :: dens(:,:)
    FLOAT,                 intent(in) :: gdens(:,:,:)
    integer,               intent(in) :: ispin
    FLOAT,                 intent(in) :: rcell_volume
  
    FLOAT, allocatable :: gnon(:)
    FLOAT :: gn(MAX_DIM), n, parameters(1)
    integer :: ii
  
    PUSH_SUB(calc_tb09_c)
  
    SAFE_ALLOCATE(gnon(1:mesh%np))
  
    do ii = 1, mesh%np
      if (ispin == UNPOLARIZED) then
        n = dens(ii, 1)
        gn(1:space%dim) = gdens(ii, 1:space%dim, 1)
      else
        n = dens(ii, 1) + dens(ii, 2)
        gn(1:space%dim) = gdens(ii, 1:space%dim, 1) + gdens(ii, 1:space%dim, 2)
      end if
  
      if (n <= CNST(1e-7)) then
        gnon(ii) = M_ZERO
      else
        gnon(ii) = sqrt(sum((gn(1:space%dim)/n)**2))
      end if
    end do
  
    parameters(1) =  -CNST(0.012) + CNST(1.023)*sqrt(dmf_integrate(mesh, gnon)/rcell_volume)
  
    call xc_f03_func_set_ext_params(functl(1)%conf, parameters)
  
    SAFE_DEALLOCATE_A(gnon)
  
    POP_SUB(calc_tb09_c)
  end subroutine calc_tb09_c
  
  ! ---------------------------------------------------------
  subroutine calc_mvorb_alpha(mesh, namespace, space, functl, dens, gdens, ispin, rcell_volume, &
    cam_alpha, cam_beta, cam_omega)
    type(mesh_t),          intent(in) :: mesh
    type(namespace_t),     intent(in) :: namespace
    type(space_t),         intent(in) :: space
    type(xc_functional_t), intent(inout) :: functl(:)
    FLOAT,                 intent(in) :: dens(:,:)
    FLOAT,                 intent(in) :: gdens(:,:,:)
    integer,               intent(in) :: ispin
    FLOAT,                 intent(in) :: rcell_volume
    FLOAT,                 intent(inout) :: cam_alpha, cam_beta, cam_omega
  
    FLOAT, allocatable :: gnon(:)
    FLOAT :: tb09_c, alpha
    FLOAT :: gn(MAX_DIM), n
    integer :: ii
    FLOAT :: parameters(3)
  
    PUSH_SUB(calc_mvorb_alpha)
  
    SAFE_ALLOCATE(gnon(1:mesh%np))
  
    do ii = 1, mesh%np
      if (ispin == UNPOLARIZED) then
        n = dens(ii, 1)
        gn(1:space%dim) = gdens(ii, 1:space%dim, 1)
      else
        n = dens(ii, 1) + dens(ii, 2)
        gn(1:space%dim) = gdens(ii, 1:space%dim, 1) + gdens(ii, 1:space%dim, 2)
      end if
  
      if (n <= CNST(1e-7)) then
        gnon(ii) = M_ZERO
      else
        gnon(ii) = sqrt(sum((gn(1:space%dim)/n)**2))
        gnon(ii) = sqrt(gnon(ii))
      end if
    end do
  
    tb09_c =  dmf_integrate(mesh, gnon)/rcell_volume
  
    SAFE_DEALLOCATE_A(gnon)
  
    select case (functl(FUNC_C)%id)
    case (XC_HYB_GGA_XC_MVORB_HSE06)
      alpha = CNST(0.121983)+CNST(0.130711)*tb09_c**4
  
      if (alpha > 1) then
        write(message(1), '(a,f6.3,a)') 'MVORB mixing parameter bigger than one (' , alpha ,').'
        call messages_warning(1, namespace=namespace)
        alpha = CNST(0.25)
      end if
  
  
      parameters(1) = alpha
      parameters(2) = cam_omega
      parameters(3) = cam_omega
      call xc_f03_func_set_ext_params(functl(FUNC_C)%conf, parameters)
      !The name is confusing. Here alpha is the beta of hybrids in functionals,
      !but is called alpha in the original paper.
      cam_beta = alpha
  
    case (XC_HYB_GGA_XC_MVORB_PBEH)
      alpha = -CNST(1.00778)+CNST(1.10507)*tb09_c
  
      if (alpha > 1) then
        write(message(1), '(a,f6.3,a)') 'MVORB mixing parameter bigger than one (' , alpha ,').'
        call messages_warning(1, namespace=namespace)
        alpha = CNST(0.25)
      end if
      if (alpha < 0) then
        write(message(1), '(a,f6.3,a)') 'MVORB mixing parameter smaller than zero (' , alpha ,').'
        call messages_warning(1, namespace=namespace)
        alpha = CNST(0.25)
      end if
  
  #if defined HAVE_LIBXC5
      parameters(1) = alpha
      call xc_f03_func_set_ext_params(functl(FUNC_C)%conf, parameters)
  #else
      call messages_not_implemented("MVORB with PBE0 requires libxc 5", namespace=namespace)
  #endif
      cam_alpha = alpha
    case default
      call messages_not_implemented("MVORB density-based mixing for functionals other than PBE0 and HSE06", namespace=namespace)
    end select
  
    POP_SUB(calc_mvorb_alpha)
  end subroutine calc_mvorb_alpha


end module xc_interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
