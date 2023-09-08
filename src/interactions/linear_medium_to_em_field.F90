!! Copyright (C) 2021 F. Bonafé
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

module linear_medium_to_em_field_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quantity_oct_m

  implicit none

  private
  public ::                    &
    linear_medium_to_em_field_t,  &
    single_medium_box_t,       &
    single_medium_box_allocate,&
    single_medium_box_end

  type single_medium_box_t
    FLOAT, allocatable            :: ep(:) !< permitivity of the linear media
    FLOAT, allocatable            :: mu(:) !< permeability of the linear media
    FLOAT, allocatable            :: c(:) !< speed of light in the linear media
    FLOAT, allocatable            :: sigma_e(:) !< electric conductivy of (lossy) medium
    FLOAT, allocatable            :: sigma_m(:) !< magnetic conductivy of (lossy) medium
    integer                       :: points_number
    integer                       :: global_points_number
    integer, allocatable          :: points_map(:)
    logical                       :: has_mapping = .true.
    integer                       :: bdry_number
    integer, allocatable          :: bdry_map(:)
    FLOAT, allocatable            :: aux_ep(:,:) !< auxiliary array for the epsilon derivative profile
    FLOAT, allocatable            :: aux_mu(:,:) !< auxiliary array for the softened mu profile
    integer                       :: box_shape
    FLOAT                         :: center(3) !< center of a box
    FLOAT                         :: lsize(3)  !< length in each direction of a box
    character(len=256)            :: filename
    logical                       :: check_medium_points = .false.
  contains
    procedure :: to_grid => single_medium_box_to_grid
  end type single_medium_box_t

  type, extends(interaction_with_partner_t) :: linear_medium_to_em_field_t
    private

    type(grid_t), pointer, public    :: system_gr !< pointer to grid of the Maxwell system
    type(single_medium_box_t), public :: partner_medium_box
    type(single_medium_box_t), public :: medium_box

  contains
    procedure :: init => linear_medium_to_em_field_init
    procedure :: calculate => linear_medium_to_em_field_calculate
    procedure :: calculate_energy => linear_medium_to_em_field_calculate_energy
    final :: linear_medium_to_em_field_finalize
  end type linear_medium_to_em_field_t


  interface linear_medium_to_em_field_t
    module procedure linear_medium_to_em_field_constructor
  end interface linear_medium_to_em_field_t

contains

  ! ---------------------------------------------------------
  function linear_medium_to_em_field_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(linear_medium_to_em_field_t),      pointer       :: this

    PUSH_SUB(linear_medium_to_em_field_constructor)

    SAFE_ALLOCATE(this)

    this%label = "linear_medium_to_em_field"
    this%partner => partner

    this%n_system_quantities = 0
    nullify(this%system_gr)

    this%n_partner_quantities = 4
    SAFE_ALLOCATE(this%partner_quantities(1:this%n_partner_quantities))
    this%partner_quantities(1) = PERMITTIVITY
    this%partner_quantities(2) = PERMEABILITY
    this%partner_quantities(3) = E_CONDUCTIVITY
    this%partner_quantities(4) = M_CONDUCTIVITY
    this%intra_interaction = .false. ! This interaction does not support intra-interactions

    POP_SUB(linear_medium_to_em_field_constructor)
  end function linear_medium_to_em_field_constructor


  subroutine linear_medium_to_em_field_init(this, gr)
    class(linear_medium_to_em_field_t), intent(inout) :: this
    type(grid_t), target, intent(in)         :: gr

    PUSH_SUB(linear_medium_to_em_field_init)

    this%system_gr => gr
    ! allocate medium box
    call single_medium_box_allocate(this%medium_box, gr%np)
    this%medium_box%has_mapping = .false.

    POP_SUB(linear_medium_to_em_field_init)
  end subroutine linear_medium_to_em_field_init

  ! ---------------------------------------------------------
  subroutine linear_medium_to_em_field_finalize(this)
    type(linear_medium_to_em_field_t), intent(inout) :: this

    PUSH_SUB(linear_medium_to_em_field_finalize)

    POP_SUB(linear_medium_to_em_field_finalize)

  end subroutine linear_medium_to_em_field_finalize


  ! ---------------------------------------------------------
  subroutine linear_medium_to_em_field_calculate(this)
    class(linear_medium_to_em_field_t), intent(inout) :: this

    PUSH_SUB(linear_medium_to_em_field_calculate)

    POP_SUB(linear_medium_to_em_field_calculate)
  end subroutine linear_medium_to_em_field_calculate

  ! ---------------------------------------------------------
  subroutine linear_medium_to_em_field_calculate_energy(this)
    class(linear_medium_to_em_field_t), intent(inout) :: this

    PUSH_SUB(linear_medium_to_em_field_calculate_energy)

    this%energy = M_ZERO

    POP_SUB(linear_medium_to_em_field_calculate_energy)
  end subroutine linear_medium_to_em_field_calculate_energy


  ! ---------------------------------------------------------
  !> Allocation of medium_box components
  subroutine single_medium_box_allocate(medium_box, n_points)
    type(single_medium_box_t),   intent(inout)    :: medium_box
    integer,                     intent(in)       :: n_points

    type(profile_t), save :: prof

    PUSH_SUB(medium_box_allocate)

    call profiling_in(prof, 'MEDIUM_BOX_ALLOC')
    SAFE_ALLOCATE(medium_box%aux_ep(1:n_points,1:3))
    SAFE_ALLOCATE(medium_box%aux_mu(1:n_points,1:3))
    SAFE_ALLOCATE(medium_box%c(1:n_points))
    SAFE_ALLOCATE(medium_box%ep(1:n_points))
    SAFE_ALLOCATE(medium_box%mu(1:n_points))
    SAFE_ALLOCATE(medium_box%sigma_e(1:n_points))
    SAFE_ALLOCATE(medium_box%sigma_m(1:n_points))
    SAFE_ALLOCATE(medium_box%points_map(1:n_points))
    medium_box%points_map = 0
    medium_box%aux_ep(:,1:3) = M_ZERO
    medium_box%aux_mu(:,1:3) = M_ZERO
    medium_box%ep(:) = M_ZERO
    medium_box%mu(:) = M_ZERO
    medium_box%c(:) = M_ZERO
    medium_box%sigma_e(:) = M_ZERO
    medium_box%sigma_m(:) = M_ZERO
    medium_box%points_number = n_points
    call profiling_out(prof)

    POP_SUB(medium_box_allocate)

  end subroutine single_medium_box_allocate

  ! ---------------------------------------------------------
  !> Deallocation of medium_box components
  subroutine single_medium_box_end(medium_box)
    type(single_medium_box_t),   intent(inout)    :: medium_box

    type(profile_t), save :: prof

    PUSH_SUB(medium_box_end)

    call profiling_in(prof, 'MEDIUM_BOX_END')

    SAFE_DEALLOCATE_A(medium_box%points_map)
    SAFE_DEALLOCATE_A(medium_box%bdry_map)
    SAFE_DEALLOCATE_A(medium_box%aux_ep)
    SAFE_DEALLOCATE_A(medium_box%aux_mu)
    SAFE_DEALLOCATE_A(medium_box%c)
    SAFE_DEALLOCATE_A(medium_box%ep)
    SAFE_DEALLOCATE_A(medium_box%mu)
    SAFE_DEALLOCATE_A(medium_box%sigma_e)
    SAFE_DEALLOCATE_A(medium_box%sigma_m)

    call profiling_out(prof)

    POP_SUB(medium_box_end)

  end subroutine single_medium_box_end

  ! ---------------------------------------------------------
  ! return a medium box with all functions mapped to the parent grid
  function single_medium_box_to_grid(medium_box, grid_out) result(medium_box_out)
    class(single_medium_box_t), intent(in) :: medium_box
    type(grid_t),               intent(in) :: grid_out
    class(single_medium_box_t), pointer :: medium_box_out

    type(profile_t), save :: prof
    integer :: ip, ip_out

    PUSH_SUB(medium_box_to_grid)
    call profiling_in(prof, 'MEDIUM_BOX_TO_GRID')

    SAFE_ALLOCATE(medium_box_out)
    call single_medium_box_allocate(medium_box_out, grid_out%np)
    medium_box_out%has_mapping = .false.

    do ip = 1, medium_box%points_number
      ip_out = medium_box%points_map(ip)
      medium_box_out%aux_ep(ip_out,1:3) = medium_box%aux_ep(ip,1:3)
      medium_box_out%aux_mu(ip_out,1:3) = medium_box%aux_mu(ip,1:3)
      medium_box_out%ep(ip_out) = medium_box%ep(ip)
      medium_box_out%mu(ip_out) = medium_box%mu(ip)
      medium_box_out%c(ip_out) = medium_box%c(ip)
      medium_box_out%sigma_e(ip_out) = medium_box%sigma_e(ip)
      medium_box_out%sigma_m(ip_out) = medium_box%sigma_m(ip)
    end do

    call profiling_out(prof)
    POP_SUB(medium_box_to_grid)
  end function single_medium_box_to_grid

end module linear_medium_to_em_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
