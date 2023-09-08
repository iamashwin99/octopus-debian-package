!! Copyright (C) 2020 M. Oliveira, Heiko Appel
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

!> @brief This module defines the abstract interaction_t class, and some auxiliary classes for interactions.
!!
!! Usually, two systems are interacting with each other in both directions. Each system affects the other system.
!! Of course, these effects can be of different nature. For instance, electrons affect a Maxwell system through their
!! current, while the Maxwell system acts on the electrons throught the electro-magnetic fields. Symmetric interactions can
!! occur between two systems of the same type.
!!
!! In Octopus, physical interactions between two systems are split into uni-directional interactions,
!! having a well defined source (causing the interaction) and a system, which is affected by the interaction.
!!
!! \dot
!!    digraph interaction {
!!     "partner A" -> "partner B" [label=" A -> B   "];
!!     "partner B" -> "partner A" [label=" B -> A   "];
!!    }
!! \enddot
!!
!! In some approximations, only one of them are used.
!!
!! The ''partners'' are instances of (interaction_partner_oct_m::interaction_partner_t) derived types.
!! The interactions (from now on uni-directional) are implemented as derived types of the abstract
!! (interaction_oct_m::interaction_t) class.
!!
!! \dot
!!   digraph uni {
!!     "partner A" -> "partner B" [label=" interaction "];
!!   }
!! \enddot
!!
!! When processing the interactions of a system, e.g. B, the interaction uses the
!! quantities, exposed by the partner system, in this example, A.
!!

module interaction_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  implicit none

  private
  public ::                &
    interaction_t,         &
    interaction_end,       &
    interaction_list_t,    &
    interaction_iterator_t

  !> @brief abstract interaction class defining only minimal structures
  !!
  !! This is the minimal data structure for interactions, which only contains information relating to partner A,
  !! who owns the interacion. In particular, the interaction has a list of quantities, which partner A needs to expose,
  !! so that the interaction can be calculated. The IDs for the exposed quantities are defined in the section 'Exposed Quantities'.
  !!
  type, abstract :: interaction_t
    private
    ! The interaction requires access to some quantities from a system to be evaluated.
    logical,              public :: intra_interaction    !< Is this an interaction of a system with itself?
    integer,              public :: n_system_quantities  !< Number of quantities needed from the system
    integer, allocatable, public :: system_quantities(:) !< Identifiers of the quantities needed from the system
    type(clock_t), public :: clock !< Clock storing the time at which the interaction was last updated.
    character(len=:), public, allocatable :: label !< label of an interaction, printed for debug purposes
    FLOAT, public :: energy !< Energy associated with the interaction.
  contains
    procedure(interaction_update),    deferred :: update
    procedure(interaction_calculate), deferred :: calculate
    procedure(interaction_calculate_energy), deferred :: calculate_energy
    procedure :: restart_read => interaction_restart_read   !< @copydoc interaction_oct_m::interaction_restart_read
    procedure :: restart_write => interaction_restart_write !< @copydoc interaction_oct_m::interaction_restart_write
  end type interaction_t

  abstract interface
    logical function interaction_update(this, requested_time)
      import interaction_t
      import clock_t
      class(interaction_t),      intent(inout) :: this
      class(clock_t),            intent(in)    :: requested_time
    end function interaction_update

    subroutine interaction_calculate(this)
      import interaction_t
      class(interaction_t),      intent(inout) :: this
    end subroutine interaction_calculate

    subroutine interaction_calculate_energy(this)
      import interaction_t
      class(interaction_t),      intent(inout) :: this
    end subroutine interaction_calculate_energy
  end interface

  !> These classes extend the list and list iterator to make an interaction list
  type, extends(linked_list_t) :: interaction_list_t
    private
  contains
    procedure :: add => interaction_list_add_node
  end type interaction_list_t

  !> These class extend the list and list iterator to make an interaction list
  type, extends(linked_list_iterator_t) :: interaction_iterator_t
    private
  contains
    procedure :: get_next => interaction_iterator_get_next
  end type interaction_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine interaction_end(this)
    class(interaction_t), intent(inout) :: this

    PUSH_SUB(interaction_end)

    SAFE_DEALLOCATE_A(this%system_quantities)
    if (allocated(this%label)) then
      deallocate(this%label)
    end if

    POP_SUB(interaction_end)
  end subroutine interaction_end

  ! ---------------------------------------------------------
  subroutine interaction_list_add_node(this, interaction)
    class(interaction_list_t)         :: this
    class(interaction_t),      target :: interaction

    PUSH_SUB(interaction_list_add_node)

    call this%add_ptr(interaction)

    POP_SUB(interaction_list_add_node)
  end subroutine interaction_list_add_node

  ! ---------------------------------------------------------
  function interaction_iterator_get_next(this) result(interaction)
    class(interaction_iterator_t), intent(inout) :: this
    class(interaction_t),          pointer       :: interaction

    PUSH_SUB(interaction_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (interaction_t)
      interaction => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(interaction_iterator_get_next)
  end function interaction_iterator_get_next

  ! ---------------------------------------------------------
  !> @brief read restart information
  !!
  !! return .true. on success.
  logical function interaction_restart_read(this, namespace)
    class(interaction_t), intent(inout) :: this
    type(namespace_t),    intent(in)    :: namespace

    PUSH_SUB(interaction_restart_read)

    interaction_restart_read = this%clock%restart_read('restart_clock_interaction_'//trim(this%label), &
      namespace)

    POP_SUB(interaction_restart_read)
  end function interaction_restart_read

  ! ---------------------------------------------------------
  subroutine interaction_restart_write(this, namespace)
    class(interaction_t), intent(inout) :: this
    type(namespace_t),    intent(in)    :: namespace

    PUSH_SUB(interaction_restart_write)

    call this%clock%restart_write('restart_clock_interaction_'//trim(this%label), namespace)

    POP_SUB(interaction_restart_write)
  end subroutine interaction_restart_write
end module interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
