!! Copyright (C) 2023 S. Ohlmann
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

module field_transfer_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use regridding_oct_m
  use restart_oct_m
  use states_mxll_oct_m
  use time_interpolation_oct_m

  implicit none

  private
  public :: &
    field_transfer_t

  type, extends(interaction_with_partner_t), abstract :: field_transfer_t
    private
    FLOAT, allocatable, public     :: partner_field(:,:) !< field from partner
    FLOAT, allocatable, public     :: system_field(:,:)  !< field transferred to system grid
    type(grid_t), pointer, public  :: system_gr !< pointer to grid of the system
    type(regridding_t), pointer, public :: regridding => NULL()
    type(time_interpolation_t), pointer, public :: interpolation => NULL()
    integer, public :: ndim
    logical, public :: interpolation_initialized = .false.

  contains
    procedure :: init => field_transfer_init
    procedure :: init_from_partner => field_transfer_init_from_partner
    procedure :: init_interpolation => field_transfer_init_interpolation
    procedure :: calculate => field_transfer_calculate
    procedure :: dfield_transfer_interpolate, zfield_transfer_interpolate
    generic :: interpolate => dfield_transfer_interpolate, zfield_transfer_interpolate
    procedure :: calculate_energy => field_transfer_calculate_energy
    procedure :: read_restart => field_transfer_read_restart
    procedure :: write_restart => field_transfer_write_restart
    procedure :: end => field_transfer_end
  end type field_transfer_t

contains

  subroutine field_transfer_init(this, gr, ndim)
    class(field_transfer_t), intent(inout) :: this
    type(grid_t), target,    intent(in)    :: gr
    integer,                 intent(in)    :: ndim

    PUSH_SUB(field_transfer_init)

    this%system_gr => gr
    this%ndim = ndim
    SAFE_ALLOCATE(this%system_field(1:gr%np, 1:ndim))
    this%system_field(:,:) = M_zero

    POP_SUB(field_transfer_init)
  end subroutine field_transfer_init

  subroutine field_transfer_init_from_partner(this, partner_gr, partner_space, partner_namespace)
    class(field_transfer_t), intent(inout) :: this
    type(grid_t),            intent(in)    :: partner_gr
    type(space_t),           intent(in)    :: partner_space
    type(namespace_t),       intent(in)    :: partner_namespace

    PUSH_SUB(field_transfer_init_from_partner)

    SAFE_ALLOCATE(this%partner_field(1:partner_gr%np, 1:this%ndim))
    this%partner_field(:,:) = M_zero
    this%regridding => regridding_t(this%system_gr, partner_gr, partner_space, partner_namespace)

    POP_SUB(field_transfer_init_from_partner)
  end subroutine field_transfer_init_from_partner

  subroutine field_transfer_init_interpolation(this, depth, label, cmplx)
    class(field_transfer_t), intent(inout) :: this
    integer,                 intent(in)    :: depth
    character(len=*),        intent(in)    :: label
    logical, optional,       intent(in)    :: cmplx

    PUSH_SUB(field_transfer_init_interpolation)

    this%interpolation => time_interpolation_t(this%system_gr%np, this%ndim, depth, &
      optional_default(cmplx, .false.), label)
    this%interpolation_initialized = .true.

    POP_SUB(field_transfer_init_interpolation)
  end subroutine field_transfer_init_interpolation

  subroutine field_transfer_end(this)
    class(field_transfer_t), intent(inout) :: this

    PUSH_SUB(field_transfer_end)

    call interaction_with_partner_end(this)
    SAFE_DEALLOCATE_A(this%partner_field)
    SAFE_DEALLOCATE_A(this%system_field)
    SAFE_DEALLOCATE_P(this%regridding)
    SAFE_DEALLOCATE_P(this%interpolation)

    POP_SUB(field_transfer_end)
  end subroutine field_transfer_end

  subroutine field_transfer_calculate(this)
    class(field_transfer_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(field_transfer_calculate)

    call profiling_in(prof,"FIELD_TRANSFER_REGRIDDING")

    ASSERT(this%interpolation_initialized)

    if (allocated(this%partner_field)) then
      call this%regridding%do_transfer(this%system_field, this%partner_field)
    end if

    call this%interpolation%add_time(this%partner%quantities(this%partner_quantities(1))%clock%time(), &
      this%system_field)

    call profiling_out(prof)
    POP_SUB(field_transfer_calculate)
  end subroutine field_transfer_calculate

  subroutine dfield_transfer_interpolate(this, time, field)
    class(field_transfer_t), intent(in)  :: this
    FLOAT,                   intent(in)  :: time
    FLOAT,                   intent(out) :: field(:, :)

    type(profile_t), save :: prof

    PUSH_SUB(dfield_transfer_interpolate)

    call profiling_in(prof,"FIELD_TRANSFER_INTERPOLATE")

    call this%interpolation%interpolate(time, field)

    call profiling_out(prof)
    POP_SUB(dfield_transfer_interpolate)
  end subroutine dfield_transfer_interpolate

  subroutine zfield_transfer_interpolate(this, time, field)
    class(field_transfer_t), intent(in)  :: this
    FLOAT,                   intent(in)  :: time
    CMPLX,                   intent(out) :: field(:, :)

    type(profile_t), save :: prof

    PUSH_SUB(zfield_transfer_interpolate)

    call profiling_in(prof,"FIELD_TRANSFER_INTERPOLATE")

    call this%interpolation%interpolate(time, field)

    call profiling_out(prof)
    POP_SUB(zfield_transfer_interpolate)
  end subroutine zfield_transfer_interpolate

  subroutine field_transfer_calculate_energy(this)
    class(field_transfer_t),    intent(inout) :: this

    PUSH_SUB(field_transfer_calculate_energy)

    ! interaction energy is zero, since it is only re-gridding the quantities of one system
    ! on the mesh of the other
    this%energy = M_ZERO

    POP_SUB(field_transfer_calculate_energy)
  end subroutine field_transfer_calculate_energy

  subroutine field_transfer_read_restart(this, mesh, space, restart, err)
    class(field_transfer_t), intent(inout) :: this
    class(mesh_t),           intent(in)    :: mesh
    type(space_t),           intent(in)    :: space
    type(restart_t),         intent(in)    :: restart
    integer,                 intent(out)   :: err

    PUSH_SUB(field_transfer_read_restart)

    call this%interpolation%read_restart(mesh, space, restart, err)

    POP_SUB(field_transfer_read_restart)
  end subroutine field_transfer_read_restart

  subroutine field_transfer_write_restart(this, mesh, space, restart, err)
    class(field_transfer_t), intent(inout) :: this
    class(mesh_t),           intent(in)    :: mesh
    type(space_t),           intent(in)    :: space
    type(restart_t),         intent(in)    :: restart
    integer,                 intent(out)   :: err

    PUSH_SUB(field_transfer_write_restart)

    call this%interpolation%write_restart(mesh, space, restart, err)

    POP_SUB(field_transfer_write_restart)
  end subroutine field_transfer_write_restart
end module field_transfer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
