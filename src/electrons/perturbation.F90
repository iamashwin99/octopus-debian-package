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

module perturbation_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use lalg_basic_oct_m
  use mpi_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use types_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                           &
    perturbation_t,                           &
    perturbation_copy

  type, abstract :: perturbation_t
    integer            :: dir
    integer            :: dir2
  contains
    procedure :: setup_dir => perturbation_setup_dir
    procedure :: apply_batch => perturbation_apply_batch
    procedure :: dexpectation_value => dperturbation_expectation_value
    procedure :: zexpectation_value => zperturbation_expectation_value
    procedure :: dstates_elec_expectation_value => dperturbation_states_elec_expectation_value
    procedure :: zstates_elec_expectation_value => zperturbation_states_elec_expectation_value
    procedure :: dexpectation_density => dperturbation_expectation_density
    procedure :: zexpectation_density => zperturbation_expectation_density
    procedure(perturbation_info),           deferred :: info
    procedure(dperturbation_apply),         deferred :: dapply
    procedure(zperturbation_apply),         deferred :: zapply
    procedure(dperturbation_apply_order_2), deferred :: dapply_order_2
    procedure(zperturbation_apply_order_2), deferred :: zapply_order_2
  end type perturbation_t

  abstract interface
    ! ---------------------------------------------------------
    subroutine dperturbation_apply(this, namespace, space, gr, hm, ik, f_in, f_out, set_bc)
      import perturbation_t
      import namespace_t
      import space_t
      import grid_t
      import hamiltonian_elec_t
      class(perturbation_t),    intent(in)    :: this
      type(namespace_t),        intent(in)    :: namespace
      type(space_t),            intent(in)    :: space
      type(grid_t),             intent(in)    :: gr
      type(hamiltonian_elec_t), intent(in)    :: hm
      integer,                  intent(in)    :: ik
      FLOAT,                    intent(in)    :: f_in(:, :)
      FLOAT,                    intent(out)   :: f_out(:, :)
      logical,        optional, intent(in)    :: set_bc
    end subroutine dperturbation_apply

    ! ---------------------------------------------------------
    subroutine zperturbation_apply(this, namespace, space, gr, hm, ik, f_in, f_out, set_bc)
      import perturbation_t
      import namespace_t
      import space_t
      import grid_t
      import hamiltonian_elec_t
      class(perturbation_t),    intent(in)    :: this
      type(namespace_t),        intent(in)    :: namespace
      type(space_t),            intent(in)    :: space
      type(grid_t),             intent(in)    :: gr
      type(hamiltonian_elec_t), intent(in)    :: hm
      integer,                  intent(in)    :: ik
      CMPLX,                    intent(in)    :: f_in(:, :)
      CMPLX,                    intent(out)   :: f_out(:, :)
      logical,        optional, intent(in)    :: set_bc
    end subroutine zperturbation_apply

    ! ---------------------------------------------------------
    subroutine dperturbation_apply_order_2(this, namespace, space, gr, hm, ik, f_in, f_out)
      import perturbation_t
      import namespace_t
      import space_t
      import grid_t
      import hamiltonian_elec_t
      class(perturbation_t),    intent(in)    :: this
      type(namespace_t),        intent(in)    :: namespace
      type(space_t),            intent(in)    :: space
      type(grid_t),             intent(in)    :: gr
      type(hamiltonian_elec_t), intent(in)    :: hm
      integer,                  intent(in)    :: ik
      FLOAT,                    intent(in)    :: f_in(:, :)
      FLOAT,                    intent(out)   :: f_out(:, :)
    end subroutine dperturbation_apply_order_2

    ! ---------------------------------------------------------
    subroutine zperturbation_apply_order_2(this, namespace, space, gr, hm, ik, f_in, f_out)
      import perturbation_t
      import namespace_t
      import space_t
      import grid_t
      import hamiltonian_elec_t
      class(perturbation_t),    intent(in)    :: this
      type(namespace_t),        intent(in)    :: namespace
      type(space_t),            intent(in)    :: space
      type(grid_t),             intent(in)    :: gr
      type(hamiltonian_elec_t), intent(in)    :: hm
      integer,                  intent(in)    :: ik
      CMPLX,                    intent(in)    :: f_in(:, :)
      CMPLX,                    intent(out)   :: f_out(:, :)
    end subroutine zperturbation_apply_order_2


    ! ---------------------------------------------------------
    subroutine perturbation_info(this)
      import perturbation_t
      class(perturbation_t),   intent(in) :: this
    end subroutine perturbation_info
  end interface

contains

  ! --------------------------------------------------------------------
  subroutine perturbation_copy(this, source)
    class(perturbation_t), intent(out) :: this
    class(perturbation_t), intent(in)  :: source

    PUSH_SUB(perturbation_copy)

    this%dir = source%dir
    this%dir2 = source%dir2

    POP_SUB(perturbation_copy)
  end subroutine perturbation_copy

  ! --------------------------------------------------------------------
  subroutine perturbation_setup_dir(this, dir, dir2)
    class(perturbation_t), intent(inout) :: this
    integer,               intent(in)    :: dir
    integer, optional,     intent(in)    :: dir2

    PUSH_SUB(perturbation_setup_dir)

    this%dir = dir
    if (present(dir2)) this%dir2 = dir2

    POP_SUB(perturbation_setup_dir)
  end subroutine perturbation_setup_dir

  ! --------------------------------------------------------------------------
  subroutine perturbation_apply_batch(this, namespace, space, gr, hm, f_in, f_out)
    class(perturbation_t),       intent(in)    :: this
    type(namespace_t),           intent(in)    :: namespace
    type(space_t),               intent(in)    :: space
    type(grid_t),                intent(in)    :: gr
    type(hamiltonian_elec_t),    intent(in)    :: hm
    type(wfs_elec_t),            intent(in)    :: f_in
    type(wfs_elec_t),            intent(inout) :: f_out

    integer :: ist
    FLOAT, allocatable :: dfi(:, :), dfo(:, :)
    CMPLX, allocatable :: zfi(:, :), zfo(:, :)

    PUSH_SUB(perturbation_apply_batch)

    ASSERT(f_in%status() == f_out%status())

    if (f_in%type() == TYPE_FLOAT) then
      SAFE_ALLOCATE(dfi(1:gr%np, 1:hm%d%dim))
      SAFE_ALLOCATE(dfo(1:gr%np, 1:hm%d%dim))

      do ist = 1, f_in%nst
        call batch_get_state(f_in, ist, gr%np, dfi)
        call this%dapply(namespace, space, gr, hm, f_in%ik, dfi, dfo)
        call batch_set_state(f_out, ist, gr%np, dfo)
      end do

      SAFE_DEALLOCATE_A(dfi)
      SAFE_DEALLOCATE_A(dfo)

    else

      SAFE_ALLOCATE(zfi(1:gr%np, 1:hm%d%dim))
      SAFE_ALLOCATE(zfo(1:gr%np, 1:hm%d%dim))

      do ist = 1, f_in%nst
        call batch_get_state(f_in, ist, gr%np, zfi)
        call this%zapply(namespace, space, gr, hm, f_in%ik, zfi, zfo)
        call batch_set_state(f_out, ist, gr%np, zfo)
      end do

      SAFE_DEALLOCATE_A(zfi)
      SAFE_DEALLOCATE_A(zfo)

    end if

    POP_SUB(perturbation_apply_batch)
  end subroutine perturbation_apply_batch


#include "undef.F90"
#include "real.F90"
#include "perturbation_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "perturbation_inc.F90"

end module perturbation_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
