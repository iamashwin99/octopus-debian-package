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

module perturbation_none_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use perturbation_oct_m
  use space_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                           &
    perturbation_none_t

  type, extends(perturbation_t) ::  perturbation_none_t
    private
  contains
    procedure :: copy => perturbation_none_copy
    generic   :: assignment(=) => copy
    procedure :: info => perturbation_none_info
    procedure :: dapply => dperturbation_none_apply
    procedure :: zapply => zperturbation_none_apply
    procedure :: apply_batch => perturbation_none_apply_batch
    procedure :: dapply_order_2 => dperturbation_none_apply_order_2
    procedure :: zapply_order_2 => zperturbation_none_apply_order_2
    final :: perturbation_none_finalize
  end type perturbation_none_t

  interface perturbation_none_t
    procedure perturbation_none_constructor
  end interface perturbation_none_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function perturbation_none_constructor(namespace) result(pert)
    class(perturbation_none_t), pointer    :: pert
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(perturbation_none_constructor)

    SAFE_ALLOCATE(pert)

    call perturbation_none_init(pert, namespace)

    POP_SUB(perturbation_none_constructor)
  end function perturbation_none_constructor

  ! --------------------------------------------------------------------
  subroutine perturbation_none_init(this, namespace)
    type(perturbation_none_t),   intent(out) :: this
    type(namespace_t),               intent(in)  :: namespace

    PUSH_SUB(perturbation_none_init)

    this%dir = -1
    this%dir2 = -1

    POP_SUB(perturbation_none_init)
  end subroutine perturbation_none_init

  ! --------------------------------------------------------------------
  subroutine perturbation_none_finalize(this)
    type(perturbation_none_t), intent(inout) :: this

    PUSH_SUB(perturbation_none_finalize)

    POP_SUB(perturbation_none_finalize)
  end subroutine perturbation_none_finalize

  ! --------------------------------------------------------------------
  subroutine perturbation_none_copy(this, source)
    class(perturbation_none_t), intent(out) :: this
    class(perturbation_none_t), intent(in)  :: source

    PUSH_SUB(perturbation_none_copy)

    call perturbation_copy(this, source)

    POP_SUB(perturbation_none_copy)
  end subroutine perturbation_none_copy


  ! --------------------------------------------------------------------
  subroutine perturbation_none_info(this)
    class(perturbation_none_t), intent(in) :: this

    PUSH_SUB(perturbation_none_info)

    POP_SUB(perturbation_none_info)
  end subroutine perturbation_none_info

  ! --------------------------------------------------------------------------
  subroutine perturbation_none_apply_batch(this, namespace, space, gr, hm, f_in, f_out)
    class(perturbation_none_t),        intent(in)    :: this
    type(namespace_t),                 intent(in)    :: namespace
    type(space_t),                     intent(in)    :: space
    type(grid_t),                      intent(in)    :: gr
    type(hamiltonian_elec_t),          intent(in)    :: hm
    type(wfs_elec_t),                  intent(in)    :: f_in
    type(wfs_elec_t),                  intent(inout) :: f_out

    PUSH_SUB(perturbation_apply_batch)

    ASSERT(f_in%status() == f_out%status())

    call batch_set_zero(f_out)

    POP_SUB(perturbation_none_apply_batch)
  end subroutine perturbation_none_apply_batch



#include "undef.F90"
#include "real.F90"
#include "perturbation_none_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "perturbation_none_inc.F90"

end module perturbation_none_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
