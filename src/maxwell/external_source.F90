!! Copyright (C) 2023 E.I. Albar, H. Appel and F. Bonafe
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

module external_source_oct_m
  use clock_oct_m
  use debug_oct_m
  use ghost_interaction_oct_m
  use global_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use interaction_with_partner_oct_m
  use interactions_factory_oct_m
  use lorentz_force_oct_m
  use maxwell_boundary_op_oct_m
  use maxwell_function_oct_m
  use mxll_field_to_medium_oct_m
  use mpi_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use plane_wave_oct_m
  use propagator_mxll_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use tdfunction_oct_m

  implicit none

  private
  public ::                       &
    external_source_t,            &
    load_external_source

  integer, public, parameter ::     &
    E_FIELD_NONE             =  0,  &
    E_FIELD_ELECTRIC         =  1

  !> @brief External source is an analytically described electromagnetic field calculated in the box.
  !! This feature couples an electromagnetic field via interactions to other systems.
  !! This field does not propagate, hence does not interfere with boundaries and does not have scattering.
  !! It solely consists of a formula evaluated on the grid and timestep of the partner system.
  type, extends(interaction_partner_t) :: external_source_t
    private

    integer, public      :: no_external_source         !< number of external sources
    type(plane_wave_t)   :: plane_wave
    logical :: ext_source_flag = .false.

  contains
    procedure :: update_exposed_quantities => external_source_update_exposed_quantities
    procedure :: init_interaction_as_partner => external_source_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => external_source_copy_quantities_to_interaction
    final :: external_source_finalize
  end type

  interface external_source_t
    module procedure external_source_constructor
  end interface external_source_t


contains

  function external_source_constructor(namespace) result(this)
    class(external_source_t), pointer :: this
    type(namespace_t), intent(in) :: namespace

    PUSH_SUB(external_source_constructor)

    SAFE_ALLOCATE(this)

    this%namespace = namespace_t("ExternalSource", parent=namespace)

    !%Variable AnalyticalExternalSource
    !%Type logical
    !%Default no
    !%Section Maxwell
    !%Description
    !% This means the analytical evaluation of formula will be used, Maxwell propagation will not be used.
    !%End
    call parse_variable(namespace, 'AnalyticalExternalSource', .false., this%ext_source_flag)
    this%space%dim = 3
    this%space%periodic_dim = 0
    this%no_external_source = 0

    if (this%ext_source_flag) then
      message(1) = 'External Source is currently always 3D and non-periodic.'
      call messages_warning(1)
      call plane_wave_init(this%plane_wave, this%namespace)
      this%no_external_source = this%plane_wave%number

      this%quantities(E_FIELD)%available_at_any_time = .true.
      this%quantities(E_FIELD)%required = .true.
      this%quantities(E_FIELD)%updated_on_demand = .false.

      ! call partners%supported_interactions_as_partner%add(LORENTZ_FORCE)
      call this%supported_interactions_as_partner%add(MXLL_FIELD_TO_MEDIUM)
      ! Initialize clock without a time-step, as the source will not be propagated
      this%clock = clock_t()
      this%quantities(E_FIELD)%clock = clock_t()
    end if

    POP_SUB(external_source_constructor)
  end function external_source_constructor

  ! ---------------------------------------------------------
  subroutine external_source_finalize(this)
    type(external_source_t), intent(inout) :: this

    PUSH_SUB(external_source_finalize)
    call plane_wave_end(this%plane_wave)

    POP_SUB(external_source_finalize)
  end subroutine external_source_finalize

  ! ---------------------------------------------------------
  logical function external_source_update_exposed_quantities(partner, requested_time, interaction) &
    result(allowed_to_update)
    class(external_source_t), intent(inout) :: partner
    type(clock_t),               intent(in)    :: requested_time
    class(interaction_t),        intent(inout) :: interaction

    integer :: iq

    PUSH_SUB(external_source_update_exposed_quantities)

    ! Always allowed to update, as the external source is not propagated
    allowed_to_update = .true.

    call partner%clock%set_time(requested_time)
    select type(interaction)
    class is (interaction_with_partner_t)
      do iq=1, interaction%n_partner_quantities
        ASSERT(partner%quantities(interaction%partner_quantities(iq))%required)
        call partner%quantities(interaction%partner_quantities(iq))%clock%set_time(requested_time)
      end do
    class default
      message(1) = "Incompatible interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(external_source_update_exposed_quantities)
  end function external_source_update_exposed_quantities

  ! ---------------------------------------------------------
  subroutine external_source_init_interaction_as_partner(partner, interaction)
    class(external_source_t),           intent(in)    :: partner
    class(interaction_t),       intent(inout) :: interaction

    PUSH_SUB(external_source_init_interaction_as_partner)

    select type (interaction)
    type is (lorentz_force_t)
      ! Nothing to be initialized for the Lorentz force.
    type is (mxll_field_to_medium_t)
      interaction%grid_based_partner = .false.
      interaction%ext_source_flag = partner%ext_source_flag
      interaction%plane_wave = partner%plane_wave

      ! Nothing to be done
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select
    POP_SUB(external_source_init_interaction_as_partner)
  end subroutine external_source_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine external_source_copy_quantities_to_interaction(partner, interaction)
    class(external_source_t),      intent(inout) :: partner
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(external_source_copy_quantities_to_interaction)

    POP_SUB(external_source_copy_quantities_to_interaction)
  end subroutine external_source_copy_quantities_to_interaction


  ! ---------------------------------------------------------
  ! Load the external source for the multisystem framework
  subroutine load_external_source(partners, namespace)
    class(partner_list_t), intent(inout)  :: partners
    type(namespace_t),    intent(in)     :: namespace

    class(external_source_t), pointer :: ext_source

    PUSH_SUB(load_external_source)

    ext_source => external_source_t(namespace)

    if(ext_source%no_external_source > 0) then
      call partners%add(ext_source)
    else
      SAFE_DEALLOCATE_P(ext_source)
    end if

    POP_SUB(load_external_source)
  end subroutine load_external_source


end module external_source_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
