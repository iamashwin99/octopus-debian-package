!! Copyright (C) 2007 X. Andrade
!! Copyright (C) 2021 N. Tancogne-Dejean
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

module perturbation_kdotp_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use perturbation_oct_m
  use profiling_oct_m
  use projector_oct_m
  use space_oct_m
  use species_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use varinfo_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                           &
    perturbation_kdotp_t

  type, extends(perturbation_t) :: perturbation_kdotp_t
    private
    integer               :: vel_method
    logical               :: use_nonlocalpps
    type(ions_t), pointer :: ions => null()
  contains
    procedure :: copy_to => perturbation_kdotp_copy
    generic   :: assignment(=) => copy_to
    procedure :: info => perturbation_kdotp_info
    procedure :: dapply => dperturbation_kdotp_apply
    procedure :: zapply => zperturbation_kdotp_apply
    procedure :: dapply_order_2 => dperturbation_kdotp_apply_order_2
    procedure :: zapply_order_2 => zperturbation_kdotp_apply_order_2
    final :: perturbation_kdotp_finalize
  end type perturbation_kdotp_t

  interface perturbation_kdotp_t
    procedure perturbation_kdotp_constructor
  end interface perturbation_kdotp_t


contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function perturbation_kdotp_constructor(namespace, ions) result(pert)
    class(perturbation_kdotp_t), pointer    :: pert
    type(namespace_t),           intent(in) :: namespace
    type(ions_t), target,        intent(in) :: ions

    PUSH_SUB(perturbation_kdotp_constructor)

    SAFE_ALLOCATE(pert)

    call perturbation_kdotp_init(pert, namespace, ions)

    POP_SUB(perturbation_kdotp_constructor)
  end function perturbation_kdotp_constructor


  ! --------------------------------------------------------------------
  subroutine perturbation_kdotp_init(this, namespace, ions)
    type(perturbation_kdotp_t), intent(out) :: this
    type(namespace_t),          intent(in)  :: namespace
    type(ions_t), target,       intent(in)  :: ions

    PUSH_SUB(perturbation_kdotp_init)

    this%dir = -1
    this%dir2 = -1

    this%ions => ions

    !%Variable KdotPUseNonLocalPseudopotential
    !%Type logical
    !%Default true
    !%Section Linear Response::KdotP
    !%Description
    !% For testing purposes, set to false to ignore the term <math>-i \left[\vec{r}, V\right]</math> in
    !% the <math>\vec{k} \cdot \vec{p}</math> perturbation, which is due to non-local pseudopotentials.
    !%End
    call messages_obsolete_variable(namespace, 'KdotP_UseNonLocalPseudopotential', 'KdotPUseNonLocalPseudopotential')
    call parse_variable(namespace, 'KdotPUseNonLocalPseudopotential', .true., this%use_nonlocalpps)

    !%Variable KdotPVelMethod
    !%Type integer
    !%Default grad_vel
    !%Section Linear Response::KdotP
    !%Description
    !% Method of velocity calculation.
    !%Option grad_vel 0
    !% <math>-i \left(\nabla + \left[r, V_{\rm nl} \right] \right)</math>
    !%Option hcom_vel 1
    !% As a commutator of the position operator and Hamiltonian, <math>-i \left[ r, H \right]</math>.
    !%End
    call parse_variable(namespace, 'KdotPVelMethod', OPTION__KDOTPVELMETHOD__GRAD_VEL, this%vel_method)

    POP_SUB(perturbation_kdotp_init)
  end subroutine perturbation_kdotp_init

  ! --------------------------------------------------------------------
  subroutine perturbation_kdotp_finalize(this)
    type(perturbation_kdotp_t), intent(inout) :: this

    PUSH_SUB(perturbation_kdotp_finalize)

    POP_SUB(perturbation_kdotp_finalize)
  end subroutine perturbation_kdotp_finalize

  ! --------------------------------------------------------------------
  subroutine perturbation_kdotp_copy(this, source)
    class(perturbation_kdotp_t), intent(out) :: this
    class(perturbation_kdotp_t), intent(in)  :: source

    PUSH_SUB(perturbation_kdotp_copy)

    call perturbation_copy(this, source)
    this%ions => source%ions

    this%use_nonlocalpps = source%use_nonlocalpps
    this%vel_method = source%vel_method

    POP_SUB(perturbation_kdotp_copy)
  end subroutine perturbation_kdotp_copy

  ! --------------------------------------------------------------------
  subroutine perturbation_kdotp_info(this)
    class(perturbation_kdotp_t), intent(in) :: this

    PUSH_SUB(perturbation_kdotp_info)

    if (.not. this%use_nonlocalpps) then
      write(message(1), '(a)') 'Ignoring non-local pseudopotential term.'
      call messages_info(1)
    end if

    POP_SUB(perturbation_kdotp_info)
  end subroutine perturbation_kdotp_info

#include "undef.F90"
#include "real.F90"
#include "perturbation_kdotp_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "perturbation_kdotp_inc.F90"

end module perturbation_kdotp_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
