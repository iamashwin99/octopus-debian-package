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

module mxll_field_to_medium_oct_m
  use clock_oct_m
  use debug_oct_m
  use field_transfer_oct_m
  use global_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use plane_wave_oct_m
  use profiling_oct_m
  use quantity_oct_m

  implicit none

  private
  public ::                    &
    mxll_field_to_medium_t

  type, extends(field_transfer_t) :: mxll_field_to_medium_t
    private
    integer, public                :: type
    type(plane_wave_t), public     :: plane_wave
    logical, public :: ext_source_flag = .false.
    logical, public :: grid_based_partner = .true.
  contains
    procedure :: calculate => mxll_field_to_medium_calculate
    final :: mxll_field_to_medium_finalize
  end type mxll_field_to_medium_t


  interface mxll_field_to_medium_t
    module procedure mxll_field_to_medium_constructor
  end interface mxll_field_to_medium_t

  integer, public, parameter :: &
    MXLL_FIELD_NONE  = -1,       &
    MXLL_FIELD_TOTAL = 0,       &
    MXLL_FIELD_TRANS = 1,       &
    MXLL_FIELD_LONG  = 2,       &
    MXLL_VEC_POT_TRANS = 3

contains

  function mxll_field_to_medium_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(mxll_field_to_medium_t), pointer :: this

    PUSH_SUB(mxll_field_to_medium_constructor)

    SAFE_ALLOCATE(this)

    this%label = "mxll_field_to_medium"
    this%partner => partner

    this%n_system_quantities = 0
    SAFE_ALLOCATE(this%system_quantities(1:this%n_system_quantities))

    this%n_partner_quantities = 1
    SAFE_ALLOCATE(this%partner_quantities(1:this%n_partner_quantities))
    this%partner_quantities(1) = E_FIELD
    this%type = MXLL_FIELD_NONE

    this%intra_interaction = .false.

    POP_SUB(mxll_field_to_medium_constructor)
  end function mxll_field_to_medium_constructor

  ! ---------------------------------------------------------
  subroutine mxll_field_to_medium_calculate(this)
    class(mxll_field_to_medium_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(mxll_field_to_medium_calculate)

    call profiling_in(prof,"MXLL_FIELD_TO_MEDIUM_CALC")

    ! In some cases partner evaluates system_e_field directly e.g. external_source
    ! If not, we have to do the regridding

    if(this%grid_based_partner) then
      call this%regridding%do_transfer(this%system_field, this%partner_field)
    end if

    if (this%ext_source_flag) then
      this%system_field = M_ZERO
      ! Loop over plane waves occurrences inside plane_waves_eval
      call plane_waves_eval(this%plane_wave, this%partner%quantities(this%partner_quantities(1))%clock%time(), &
        this%system_gr, this%system_field)
    end if

    call this%interpolation%add_time(this%partner%quantities(this%partner_quantities(1))%clock%time(), this%system_field)

    call profiling_out(prof)
    POP_SUB(mxll_field_to_medium_calculate)

  end subroutine mxll_field_to_medium_calculate

! ---------------------------------------------------------
  subroutine mxll_field_to_medium_finalize(this)
    type(mxll_field_to_medium_t), intent(inout) :: this

    PUSH_SUB(mxll_field_to_medium_finalize)

    call this%end()

    POP_SUB(mxll_field_to_medium_finalize)
  end subroutine mxll_field_to_medium_finalize

end module mxll_field_to_medium_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
