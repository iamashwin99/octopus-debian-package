!! Copyright (C) 2022 F. Bonafé
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

module current_to_mxll_field_oct_m
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
  use states_mxll_oct_m

  implicit none

  private
  public ::                    &
    current_to_mxll_field_t
  
  type, extends(interaction_with_partner_t) :: current_to_mxll_field_t
    private

    type(grid_t), pointer, public  :: system_gr !< pointer to grid of the Maxwell system
    integer, allocatable, public   :: partner_points_map(:)
    FLOAT, allocatable, public     :: partner_current_p(:,:) !< polarization current, size given by number of points on partner grid
    FLOAT, allocatable             :: system_current_p(:,:) !< polarization current, size given by number of points on system grid
    CMPLX, allocatable, public     :: rs_current_p(:,:) !< polarization current density, size given by number of points on system grid
    integer, allocatable, public   :: partner_to_system_map(:)
    integer, public                :: partner_points_number

  contains
    procedure :: init => current_to_mxll_field_init
    procedure :: calculate => current_to_mxll_field_calculate
    procedure :: calculate_energy => current_to_mxll_field_calculate_energy
    final :: current_to_mxll_field_finalize
  end type current_to_mxll_field_t


  interface current_to_mxll_field_t
    module procedure current_to_mxll_field_constructor
  end interface current_to_mxll_field_t

contains

  ! ---------------------------------------------------------
  function current_to_mxll_field_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(current_to_mxll_field_t),               pointer       :: this

    PUSH_SUB(current_to_mxll_field_constructor)

    SAFE_ALLOCATE(this)

    this%label = "current_to_mxll_field"
    this%partner => partner

    this%n_system_quantities = 0
    SAFE_ALLOCATE(this%system_quantities(1:this%n_system_quantities))
    nullify(this%system_gr)

    this%n_partner_quantities = 1
    SAFE_ALLOCATE(this%partner_quantities(1:this%n_partner_quantities))
    this%partner_quantities(1) = CURRENT

    POP_SUB(current_to_mxll_field_constructor)
  end function current_to_mxll_field_constructor


  subroutine current_to_mxll_field_init(this, gr)
    class(current_to_mxll_field_t), intent(inout) :: this
    type(grid_t), target, intent(in)         :: gr

    PUSH_SUB(current_to_mxll_field_init)

    this%system_gr => gr
    SAFE_ALLOCATE(this%rs_current_p(gr%mesh%np, gr%box%dim))
    SAFE_ALLOCATE(this%system_current_p(gr%mesh%np, gr%box%dim))
    this%rs_current_p = M_z0
    this%system_current_p = M_ZERO

    POP_SUB(current_to_mxll_field_init)
  end subroutine current_to_mxll_field_init

  ! ---------------------------------------------------------
  subroutine current_to_mxll_field_finalize(this)
    type(current_to_mxll_field_t), intent(inout) :: this

    PUSH_SUB(current_to_mxll_field_finalize)

    call interaction_with_partner_end(this)
    SAFE_DEALLOCATE_A(this%rs_current_p)
    SAFE_DEALLOCATE_A(this%system_current_p)
    SAFE_DEALLOCATE_A(this%partner_current_p)
    SAFE_DEALLOCATE_A(this%partner_points_map)
    SAFE_DEALLOCATE_A(this%partner_to_system_map)

    POP_SUB(current_to_mxll_field_finalize)
  end subroutine current_to_mxll_field_finalize

  ! ---------------------------------------------------------
  subroutine current_to_mxll_field_calculate(this)
    class(current_to_mxll_field_t), intent(inout) :: this

    integer :: ip, ip_mxll, ip_in

    type(profile_t), save :: prof

    PUSH_SUB(current_to_mxll_field_calculate)

    call profiling_in(prof,"CURRENT_TO_MXLL_FIELD_CALCULATE")

    do ip_in = 1, this%partner_points_number
      ip = this%partner_points_map(ip_in)
      ip_mxll = this%partner_to_system_map(ip)
      if (ip_mxll > 0) then
        this%system_current_p(ip_mxll,:) = this%partner_current_p(ip,:)
      end if
    end do
    call build_rs_current_state(this%system_current_p, this%system_gr%mesh, this%rs_current_p)

    call profiling_out(prof)

    POP_SUB(current_to_mxll_field_calculate)
  end subroutine current_to_mxll_field_calculate

  ! ---------------------------------------------------------
  subroutine current_to_mxll_field_calculate_energy(this)
    class(current_to_mxll_field_t),    intent(inout) :: this

    PUSH_SUB(current_to_mxll_field_calculate_energy)

    ! interaction energy is zero, since it is only re-gridding the quantities of one system
    ! on the mesh of the other
    this%energy = M_ZERO

    POP_SUB(current_to_mxll_field_calculate_energy)
  end subroutine current_to_mxll_field_calculate_energy

end module current_to_mxll_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
