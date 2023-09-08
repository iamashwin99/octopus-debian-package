!! Copyright (C) 2020 H. Appel, S. Ohlmann, M. Oliveira, N. Tancogne-Dejean
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

module charged_particle_oct_m
  use algorithm_oct_m
  use classical_particle_oct_m
  use clock_oct_m
  use coulomb_force_oct_m
  use current_to_mxll_field_oct_m
  use debug_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use lorentz_force_oct_m
  use messages_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_oct_m

  implicit none

  private
  public ::               &
    charged_particle_t,   &
    charged_particle_init

  !> @brief class for a charged classical particle
  !!
  type, extends(classical_particle_t) :: charged_particle_t
    private

    FLOAT :: charge(1)

  contains
    procedure :: init_interaction => charged_particle_init_interaction
    procedure :: initial_conditions => charged_particle_initial_conditions
    procedure :: is_tolerance_reached => charged_particle_is_tolerance_reached
    procedure :: update_quantity => charged_particle_update_quantity
    procedure :: update_exposed_quantity => charged_particle_update_exposed_quantity
    procedure :: init_interaction_as_partner => charged_particle_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => charged_particle_copy_quantities_to_interaction
  end type charged_particle_t

  interface charged_particle_t
    procedure charged_particle_constructor
  end interface charged_particle_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function charged_particle_constructor(namespace) result(sys)
    class(charged_particle_t), pointer  :: sys
    type(namespace_t),       intent(in) :: namespace

    PUSH_SUB(charged_particle_constructor)

    SAFE_ALLOCATE(sys)

    call charged_particle_init(sys, namespace)

    POP_SUB(charged_particle_constructor)
  end function charged_particle_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  subroutine charged_particle_init(this, namespace)
    class(charged_particle_t), intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace

    PUSH_SUB(charged_particle_init)

    call classical_particle_init(this%classical_particle_t, namespace)

    !%Variable ParticleCharge
    !%Type float
    !%Section ClassicalParticles
    !%Description
    !% Charge of classical particle
    !%End
    call parse_variable(namespace, 'ParticleCharge', M_ONE, this%charge(1))
    call messages_print_var_value('ParticleCharge', this%charge(1), namespace=namespace)

    call this%supported_interactions%add(LORENTZ_FORCE)
    call this%supported_interactions%add(COULOMB_FORCE)
    call this%supported_interactions_as_partner%add(COULOMB_FORCE)
    call this%supported_interactions_as_partner%add(CURRENT_TO_MXLL_FIELD)

    this%quantities(CURRENT)%required = .true.
    this%quantities(CHARGE)%updated_on_demand = .false.
    this%quantities(CHARGE)%available_at_any_time = .true.


    POP_SUB(charged_particle_init)
  end subroutine charged_particle_init

  ! ---------------------------------------------------------
  subroutine charged_particle_init_interaction(this, interaction)
    class(charged_particle_t), target, intent(inout) :: this
    class(interaction_t),              intent(inout) :: interaction

    PUSH_SUB(charged_particle_init_interaction)

    select type (interaction)
    type is (coulomb_force_t)
      call interaction%init(this%space%dim, 1, this%quantities, this%charge, this%pos)
    type is (lorentz_force_t)
      call interaction%init(this%space%dim, 1, this%quantities, this%charge, this%pos, this%vel, this%namespace)
    class default
      call this%classical_particle_t%init_interaction(interaction)
    end select

    POP_SUB(charged_particle_init_interaction)
  end subroutine charged_particle_init_interaction

  ! ---------------------------------------------------------
  subroutine charged_particle_initial_conditions(this)
    class(charged_particle_t), intent(inout) :: this

    PUSH_SUB(charged_particle_initial_conditions)

    call this%classical_particle_t%initial_conditions()

    POP_SUB(charged_particle_initial_conditions)
  end subroutine charged_particle_initial_conditions

  ! ---------------------------------------------------------
  logical function charged_particle_is_tolerance_reached(this, tol) result(converged)
    class(charged_particle_t), intent(in) :: this
    FLOAT,                     intent(in) :: tol

    PUSH_SUB(charged_particle_is_tolerance_reached)

    converged = this%classical_particle_t%is_tolerance_reached(tol)

    POP_SUB(charged_particle_is_tolerance_reached)
  end function charged_particle_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine charged_particle_update_quantity(this, iq)
    class(charged_particle_t), intent(inout) :: this
    integer,                   intent(in)    :: iq

    PUSH_SUB(charged_particle_update_quantity)

    select case (iq)
    case(CURRENT)
      ! The charged particle has a velocity, giving rise to the current.
      ! Therefore we set it to be in sync with the particle velocity
      call this%quantities(iq)%clock%set_time(this%quantities(VELOCITY)%clock)
    case default
      ! Other quantities should be handled by the parent class
      call this%classical_particle_t%update_quantity(iq)
    end select

    POP_SUB(charged_particle_update_quantity)
  end subroutine charged_particle_update_quantity

  ! ---------------------------------------------------------
  subroutine charged_particle_update_exposed_quantity(partner, iq)
    class(charged_particle_t), intent(inout) :: partner
    integer,                   intent(in)    :: iq

    PUSH_SUB(charged_particle_update_exposed_quantity)

    select case (iq)
    case (CURRENT)
      ! We still need to set its clock, so we set it to be in sync with the particle velocity.
      call partner%quantities(iq)%clock%set_time(partner%quantities(VELOCITY)%clock)
      call multisystem_debug_write_marker(partner%namespace, &
        event_clock_update_t("quantity", QUANTITY_LABEL(iq), partner%quantities(iq)%clock, "set"))
    case default
      ! Other quantities should be handled by the parent class
      call partner%classical_particle_t%update_exposed_quantity(iq)
    end select

    POP_SUB(charged_particle_update_exposed_quantity)
  end subroutine charged_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine charged_particle_init_interaction_as_partner(partner, interaction)
    class(charged_particle_t), intent(in)    :: partner
    class(interaction_t),      intent(inout) :: interaction

    PUSH_SUB(charged_particle_init_interaction_as_partner)

    select type (interaction)
    type is (coulomb_force_t)
      interaction%partner_np = 1
      SAFE_ALLOCATE(interaction%partner_charge(1))
      SAFE_ALLOCATE(interaction%partner_pos(1:partner%space%dim, 1))
    type is (current_to_mxll_field_t)
      interaction%partner_np = 1
      interaction%grid_based_partner = .false.
      SAFE_ALLOCATE(interaction%partner_charge(1))
      SAFE_ALLOCATE(interaction%partner_pos(1:partner%space%dim, 1))
      SAFE_ALLOCATE(interaction%partner_vel(1:partner%space%dim, 1))
    class default
      ! Other interactions should be handled by the parent class
      call partner%classical_particle_t%init_interaction_as_partner(interaction)
    end select

    POP_SUB(charged_particle_init_interaction_as_partner)
  end subroutine charged_particle_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine charged_particle_copy_quantities_to_interaction(partner, interaction)
    class(charged_particle_t), intent(inout) :: partner
    class(interaction_t),      intent(inout) :: interaction

    PUSH_SUB(charged_particle_copy_quantities_to_interaction)

    select type (interaction)
    type is (coulomb_force_t)
      interaction%partner_charge(1) = partner%charge(1)
      interaction%partner_pos(:,1) = partner%pos(:, 1)

    type is (current_to_mxll_field_t)
      interaction%partner_charge(1) = partner%charge(1)
      interaction%partner_pos(:,1) = partner%pos(:, 1)
      interaction%partner_vel(:,1) = partner%vel(:, 1)

    class default
      call partner%classical_particle_t%copy_quantities_to_interaction(interaction)
    end select

    POP_SUB(charged_particle_copy_quantities_to_interaction)
  end subroutine charged_particle_copy_quantities_to_interaction

end module charged_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
