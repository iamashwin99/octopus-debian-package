!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
!! Copyright (C) 2021 S. Ohlmann
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

module multigrid_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use index_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_init_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use par_vec_oct_m
  use space_oct_m
  use stencil_oct_m
  use transfer_table_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                         &
    multigrid_level_t,              &
    multigrid_t,                    &
    multigrid_init,                 &
    multigrid_end,                  &
    multigrid_mesh_half,            &
    multigrid_mesh_double,          &
    multigrid_number_of_levels,     &
    dmultigrid_fine2coarse,         &
    zmultigrid_fine2coarse,         &
    dmultigrid_coarse2fine,         &
    zmultigrid_coarse2fine,         &
    dmultigrid_fine2coarse_batch,   &
    zmultigrid_fine2coarse_batch,   &
    dmultigrid_coarse2fine_batch,   &
    zmultigrid_coarse2fine_batch,   &
    multigrid_get_transfer_tables

  integer, parameter, public :: &
    INJECTION  = 1,             &
    FULLWEIGHT = 2

  type multigrid_level_t
    ! Components are public by default
    type(transfer_table_t)          :: tt
    type(mesh_t),          pointer  :: mesh => NULL()
    type(derivatives_t),   pointer  :: der  => NULL()
  end type multigrid_level_t

  type multigrid_t
    private
    integer                                      :: n_levels
    type(multigrid_level_t), allocatable, public :: level(:)

    integer              :: tp
    integer, allocatable :: sp(:)
    integer, allocatable :: ep(:)
    integer, allocatable :: ep_part(:)
  end type multigrid_t

  type(profile_t), save :: interp_prof, injection_prof, restrict_prof
  type(profile_t), save :: interp_prof_batch, restrict_prof_batch

contains

  ! ---------------------------------------------------------
  subroutine multigrid_init(mgrid, namespace, space, mesh, der, stencil, mc, used_for_preconditioner)
    type(multigrid_t),     target, intent(out) :: mgrid
    type(namespace_t),             intent(in)  :: namespace
    type(space_t),                 intent(in)  :: space
    class(mesh_t),         target, intent(in)  :: mesh
    type(derivatives_t),   target, intent(in)  :: der
    type(stencil_t),               intent(in)  :: stencil
    type(multicomm_t),             intent(in)  :: mc
    logical, optional,             intent(in)  :: used_for_preconditioner

    integer :: i, n_levels, np, order

    PUSH_SUB(multigrid_init)

    !%Variable MultigridLevels
    !%Type integer
    !%Default max_levels
    !%Section Mesh
    !%Description
    !% Number of levels in the grid hierarchy used for multigrid. Positive
    !% numbers indicate an absolute number of levels, negative
    !% numbers are subtracted from the maximum number of levels possible.
    !%Option max_levels 0
    !% Calculate the optimal number of levels for the grid.
    !%End

    call parse_variable(namespace, 'MultigridLevels', 0, n_levels)

    ! default:
    order = der%order
    if (optional_default(used_for_preconditioner, .false.)) then
      n_levels = 3
      write(message(1), '(a)') "Set number of multigrid levels to 3 for preconditioner. This ignores the value of MultigridLevels."
      call messages_info(1, namespace=namespace)

      !%Variable MultigridDerivativesOrder
      !%Type integer
      !%Default 1
      !%Section Mesh::Derivatives
      !%Description
      !% This variable gives the discretization order for the approximation of
      !% the differential operators on the different levels of the multigrid.
      !% For more details, see the variable DerivativesOrder.
      !%End
      call parse_variable(namespace, 'MultigridDerivativesOrder', 1, order)
      ! set order to a minimum of 2 for general star stencil, fails otherwise
      ! the parameter DER_STARGENERAL is private to the derivatives module
      if (der%stencil_type == 5) then
        order = max(2, order)
      end if
    end if

    if (n_levels <= 0) then
      n_levels = n_levels - 3
      np = mesh%np

      do while (np > 1)
        np = np / 8
        n_levels = n_levels + 1
      end do
    else
      n_levels = n_levels - 1
    end if

    mgrid%n_levels = n_levels

    SAFE_ALLOCATE(mgrid%level(0:n_levels))

    mgrid%level(0)%mesh => mesh
    mgrid%level(0)%der => der

    mgrid%level(0)%tt%n_fine = mesh%np
    SAFE_ALLOCATE(mgrid%level(0)%tt%fine_i(1:mesh%np))

    write(message(1), '(a,i3)') "Multigrid levels:", n_levels + 1
    call messages_info(1, namespace=namespace)

    do i = 1, mgrid%n_levels
      SAFE_ALLOCATE(mgrid%level(i)%mesh)
      SAFE_ALLOCATE(mgrid%level(i)%der)

      call multigrid_mesh_half(space, namespace, mgrid%level(i-1)%mesh, mgrid%level(i)%mesh, stencil)

      call derivatives_init(mgrid%level(i)%der, namespace, space, mesh%coord_system, order=order)

      call mesh_init_stage_3(mgrid%level(i)%mesh, namespace, space, stencil, mc, parent = mgrid%level(i - 1)%mesh)

      call multigrid_get_transfer_tables(mgrid%level(i)%tt, mgrid%level(i-1)%mesh, mgrid%level(i)%mesh)

      call derivatives_build(mgrid%level(i)%der, namespace, space, mgrid%level(i)%mesh)

      call mesh_write_info(mgrid%level(i)%mesh, namespace=namespace)

      mgrid%level(i)%der%finer => mgrid%level(i - 1)%der
      mgrid%level(i - 1)%der%coarser => mgrid%level(i)%der
      mgrid%level(i)%der%to_finer => mgrid%level(i)%tt
      mgrid%level(i - 1)%der%to_coarser => mgrid%level(i)%tt
    end do

    SAFE_ALLOCATE(mgrid%sp(0:mgrid%n_levels))
    SAFE_ALLOCATE(mgrid%ep(0:mgrid%n_levels))
    SAFE_ALLOCATE(mgrid%ep_part(0:mgrid%n_levels))

    mgrid%tp = 0
    do i = 0, mgrid%n_levels
      mgrid%sp(i) = 1 + mgrid%tp
      mgrid%ep(i) = mgrid%tp + mgrid%level(i)%mesh%np
      mgrid%tp = mgrid%tp + mgrid%level(i)%mesh%np_part
      mgrid%ep_part(i) = mgrid%tp
    end do

    POP_SUB(multigrid_init)
  end subroutine multigrid_init

  ! ---------------------------------------------------------
  !> creates the lookup tables to go between the coarse and fine meshes
  subroutine multigrid_get_transfer_tables(tt, fine, coarse)
    type(transfer_table_t), intent(inout) :: tt
    type(mesh_t),           intent(in)    :: fine, coarse

    integer :: i, i1, i2, i4, i8, pt
    integer :: x(MAX_DIM), mod2(MAX_DIM), idx(MAX_DIM)

    PUSH_SUB(multigrid_get_transfer_tables)

    tt%n_coarse = coarse%np
    SAFE_ALLOCATE(tt%to_coarse(1:tt%n_coarse))

    ! GENERATE THE TABLE TO MAP FROM THE FINE TO THE COARSE GRID
    do i = 1, tt%n_coarse
      ! locate the equivalent fine grid point
      call mesh_local_index_to_coords(coarse, i, idx)
      tt%to_coarse(i) = mesh_local_index_from_coords(fine, 2*idx)
    end do

    ! count
    tt%n_fine = fine%np
    SAFE_ALLOCATE(tt%fine_i(1:tt%n_fine))

    tt%n_fine1 = 0
    tt%n_fine2 = 0
    tt%n_fine4 = 0
    tt%n_fine8 = 0
    do i = 1, tt%n_fine
      call mesh_local_index_to_coords(fine, i, idx)
      mod2 = mod(idx, 2)

      pt = sum(abs(mod2(1:3)))

      select case (pt)
      case (0)
        tt%n_fine1 = tt%n_fine1 + 1
        tt%fine_i(i) = 1
      case (1)
        tt%n_fine2 = tt%n_fine2 + 1
        tt%fine_i(i) = 2
      case (2)
        tt%n_fine4 = tt%n_fine4 + 1
        tt%fine_i(i) = 4
      case (3)
        tt%n_fine8 = tt%n_fine8 + 1
        tt%fine_i(i) = 8
      end select
    end do

    ASSERT(tt%n_fine1 + tt%n_fine2 + tt%n_fine4 + tt%n_fine8 == tt%n_fine)

    SAFE_ALLOCATE(tt%to_fine1(1, 1:tt%n_fine1))
    SAFE_ALLOCATE(tt%to_fine2(1:2, 1:tt%n_fine2))
    SAFE_ALLOCATE(tt%to_fine4(1:4, 1:tt%n_fine4))
    SAFE_ALLOCATE(tt%to_fine8(1:8, 1:tt%n_fine8))

    ! and now build the tables
    i1 = 0
    i2 = 0
    i4 = 0
    i8 = 0
    do i = 1, fine%np
      call mesh_local_index_to_coords(fine, i, idx)
      x(1:3)    = idx(1:3)/2
      mod2(1:3) = mod(idx(1:3), 2)

      pt = sum(abs(mod2(1:3)))

      select case (pt)
      case (0)
        i1 = i1 + 1
        tt%to_fine1(1, i1) = mesh_local_index_from_coords(coarse, [x(1), x(2), x(3)])

      case (1)
        i2 = i2 + 1
        tt%to_fine2(1, i2) = mesh_local_index_from_coords(coarse, [x(1)          , x(2)          , x(3)          ])
        tt%to_fine2(2, i2) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2) + mod2(2), x(3) + mod2(3)])

      case (2)
        i4 = i4 + 1
        tt%to_fine4(1, i4) = mesh_local_index_from_coords(coarse, [x(1)          , x(2) + mod2(2), x(3) + mod2(3)])
        tt%to_fine4(2, i4) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2)          , x(3) + mod2(3)])
        tt%to_fine4(3, i4) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2) + mod2(2), x(3)          ])
        tt%to_fine4(4, i4) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2) + mod2(2), x(3) + mod2(3)])

      case (3)
        i8 = i8 + 1
        tt%to_fine8(1, i8) = mesh_local_index_from_coords(coarse, [x(1)          , x(2)          , x(3)          ])
        tt%to_fine8(2, i8) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2)          , x(3)          ])
        tt%to_fine8(3, i8) = mesh_local_index_from_coords(coarse, [x(1)          , x(2) + mod2(2), x(3)          ])
        tt%to_fine8(4, i8) = mesh_local_index_from_coords(coarse, [x(1)          , x(2)          , x(3) + mod2(3)])
        tt%to_fine8(5, i8) = mesh_local_index_from_coords(coarse, [x(1)          , x(2) + mod2(2), x(3) + mod2(3)])
        tt%to_fine8(6, i8) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2)          , x(3) + mod2(3)])
        tt%to_fine8(7, i8) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2) + mod2(2), x(3)          ])
        tt%to_fine8(8, i8) = mesh_local_index_from_coords(coarse, [x(1) + mod2(1), x(2) + mod2(2), x(3) + mod2(3)])

      end select

    end do

    ASSERT(i1 == tt%n_fine1 .and. i2 == tt%n_fine2 .and. i4 == tt%n_fine4 .and. i8 == tt%n_fine8)


    POP_SUB(multigrid_get_transfer_tables)
  end subroutine multigrid_get_transfer_tables

  !---------------------------------------------------------------------------------
  !> Creates a mesh that has twice the spacing betwen the points than the in mesh.
  !! This is used in the multi-grid routines
  !---------------------------------------------------------------------------------
  subroutine multigrid_mesh_half(space, namespace, mesh_in, mesh_out, stencil)
    type(space_t),              intent(in)    :: space
    type(namespace_t),          intent(in)    :: namespace
    type(mesh_t),       target, intent(in)    :: mesh_in
    type(mesh_t),               intent(inout) :: mesh_out
    type(stencil_t),            intent(in)    :: stencil

    integer :: idim

    PUSH_SUB(multigrid_mesh_half)

    mesh_out%box              => mesh_in%box
    mesh_out%idx%dim          =  mesh_in%idx%dim
    mesh_out%use_curvilinear  =  mesh_in%use_curvilinear
    mesh_out%masked_periodic_boundaries = mesh_in%masked_periodic_boundaries
    mesh_out%coord_system     => mesh_in%coord_system

    mesh_out%spacing(:) = 2*mesh_in%spacing(:)
    mesh_out%idx%enlarge(:) = mesh_in%idx%enlarge(:)
    mesh_out%idx%nr(1,:) = (mesh_in%idx%nr(1,:)+mesh_in%idx%enlarge(:))/2
    mesh_out%idx%nr(2,:) = (mesh_in%idx%nr(2,:)-mesh_in%idx%enlarge(:))/2
    mesh_out%idx%ll(:) = mesh_out%idx%nr(2, :) - mesh_out%idx%nr(1, :) + 1

    mesh_out%idx%stride(1) = 1
    do idim = 2, MAX_DIM+1
      mesh_out%idx%stride(idim) = mesh_out%idx%stride(idim-1) *  &
        (mesh_out%idx%ll(idim-1) + 2*mesh_out%idx%enlarge(idim-1))
    end do

    call mesh_init_stage_2(mesh_out, namespace, space, mesh_out%box, stencil)

    POP_SUB(multigrid_mesh_half)
  end subroutine multigrid_mesh_half

  !---------------------------------------------------------------------------------
  subroutine multigrid_mesh_double(space, namespace, mesh_in, mesh_out, stencil)
    type(space_t),              intent(in)    :: space
    type(namespace_t),          intent(in)    :: namespace
    type(mesh_t),       target, intent(in)    :: mesh_in
    type(mesh_t),               intent(inout) :: mesh_out
    type(stencil_t),            intent(in)    :: stencil

    integer :: idim
    PUSH_SUB(multigrid_mesh_double)

    mesh_out%box              => mesh_in%box
    mesh_out%idx%dim          =  mesh_in%idx%dim
    mesh_out%use_curvilinear =  mesh_in%use_curvilinear
    mesh_out%masked_periodic_boundaries = mesh_in%masked_periodic_boundaries
    mesh_out%coord_system    => mesh_in%coord_system

    mesh_out%spacing(:) = M_HALF*mesh_in%spacing(:)
    mesh_out%idx%enlarge(:) = mesh_in%idx%enlarge(:)
    mesh_out%idx%nr(1,:) = (mesh_in%idx%nr(1,:)+mesh_in%idx%enlarge(:))*2
    mesh_out%idx%nr(2,:) = (mesh_in%idx%nr(2,:)-mesh_in%idx%enlarge(:))*2
    ! We need to make the possible number of grid points larger by one for
    ! the non-periodic dimensions because the spacing is only half of the
    ! original mesh and thus we could get points in the new boundary
    ! that are still inside the simulation box.
    ! For the periodic dimensions, we are anyway commensurate with the size
    ! of the box, so we are still commensurate when taking twice the number
    ! of points.
    do idim = space%periodic_dim + 1, space%dim
      mesh_out%idx%nr(1, idim) = mesh_out%idx%nr(1, idim) - 1
      mesh_out%idx%nr(2, idim) = mesh_out%idx%nr(2, idim) + 1
    end do
    mesh_out%idx%ll(:) = mesh_out%idx%nr(2, :) - mesh_out%idx%nr(1, :) + 1

    mesh_out%idx%stride(1) = 1
    do idim = 2, MAX_DIM+1
      mesh_out%idx%stride(idim) = mesh_out%idx%stride(idim-1) *  &
        (mesh_out%idx%ll(idim-1) + 2*mesh_out%idx%enlarge(idim-1))
    end do

    call mesh_init_stage_2(mesh_out, namespace, space, mesh_out%box, stencil)

    POP_SUB(multigrid_mesh_double)
  end subroutine multigrid_mesh_double

  ! ---------------------------------------------------------
  subroutine multigrid_end(mgrid)
    type(multigrid_t), target, intent(inout) :: mgrid

    integer :: i
    type(multigrid_level_t), pointer :: level

    PUSH_SUB(multigrid_end)

    SAFE_DEALLOCATE_A(mgrid%sp)
    SAFE_DEALLOCATE_A(mgrid%ep)
    SAFE_DEALLOCATE_A(mgrid%ep_part)

    SAFE_DEALLOCATE_A(mgrid%level(0)%tt%fine_i)

    do i = 1, mgrid%n_levels
      level => mgrid%level(i)

      call derivatives_end(level%der)
      call mesh_end(level%mesh)
      SAFE_DEALLOCATE_P(level%mesh)
      SAFE_DEALLOCATE_P(level%der)

      SAFE_DEALLOCATE_A(level%tt%to_coarse)
      SAFE_DEALLOCATE_A(level%tt%to_fine1)
      SAFE_DEALLOCATE_A(level%tt%to_fine2)
      SAFE_DEALLOCATE_A(level%tt%to_fine4)
      SAFE_DEALLOCATE_A(level%tt%to_fine8)
      SAFE_DEALLOCATE_A(level%tt%fine_i)
    end do

    SAFE_DEALLOCATE_A(mgrid%level)

    POP_SUB(multigrid_end)
  end subroutine multigrid_end

  !---------------------------------------------------------------------------------
  integer function multigrid_number_of_levels(base_der) result(number)
    type(derivatives_t), target, intent(in)  :: base_der

    type(derivatives_t), pointer :: next_der

    next_der => base_der%coarser

    number = 0
    do
      number = number + 1
      next_der => next_der%coarser
      if (.not. associated(next_der)) exit
    end do

  end function multigrid_number_of_levels


#include "undef.F90"
#include "real.F90"
#include "multigrid_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "multigrid_inc.F90"

end module multigrid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
