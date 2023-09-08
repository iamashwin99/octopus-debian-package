!! Copyright (C) 2020 M. Oliveira
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

!> @brief This module defines classes and functions for interaction partners.
!!
!! Interaction partners are general objects, which define the "source" of an interaction (interaction_oct_m::interaction_t).
module interaction_partner_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use interaction_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use quantity_oct_m
  use space_oct_m
  implicit none

  private
  public ::                &
    interaction_partner_t, &
    partner_list_t,        &
    partner_iterator_t

  !> @brief abstract class for general interaction partners
  !!
  !! Some interactions require a partner. This is usually a system, but it could
  !! also be some external entity, like an external field.
  !! An interaction partner must expose some quantities that the interaction can use.
  type, abstract :: interaction_partner_t
    private
    type(namespace_t), public :: namespace
    type(clock_t),     public :: clock
    type(space_t),     public :: space

    type(integer_list_t), public :: supported_interactions_as_partner !< list of interactions, which support
    !!                                                                   this interaction_partner_t as partner

    type(quantity_t),  public :: quantities(MAX_QUANTITIES) !< Array of all possible quantities.
    !!                                                         The elements of the array are accessed using the
    !!                                                         quantity`s identifiers.
  contains
    procedure(interaction_partner_update_exposed_quantities),      deferred :: update_exposed_quantities
    !< @copydoc interaction_partner_oct_m::interaction_partner_update_exposed_quantities
    procedure(interaction_partner_init_interaction_as_partner),    deferred :: init_interaction_as_partner
    !< @copydoc interaction_partner_oct_m::interaction_partner_init_interaction_as_partner
    procedure(interaction_partner_copy_quantities_to_interaction), deferred :: copy_quantities_to_interaction
    !< @copydoc interaction_partner_oct_m::interaction_partner_copy_quantities_to_interaction
  end type interaction_partner_t

  abstract interface
    ! ---------------------------------------------------------
    !> @brief update all exposed quantities as interaction partner to a requested time
    !!
    !! @result update was sucessful
    logical function interaction_partner_update_exposed_quantities(partner, requested_time, interaction)
      import interaction_partner_t
      import clock_t
      import interaction_t
      class(interaction_partner_t), intent(inout) :: partner         !< the current interacion partner (this)
      type(clock_t),                intent(in)    :: requested_time  !< the time at which we need the quantities
      class(interaction_t),         intent(inout) :: interaction     !< the interaction
    end function interaction_partner_update_exposed_quantities

    ! ---------------------------------------------------------
    subroutine interaction_partner_init_interaction_as_partner(partner, interaction)
      import interaction_partner_t
      import interaction_t
      class(interaction_partner_t),     intent(in)    :: partner     !< the current interacion partner (this)
      class(interaction_t),             intent(inout) :: interaction !< the interaction
    end subroutine interaction_partner_init_interaction_as_partner

    ! ---------------------------------------------------------
    subroutine interaction_partner_copy_quantities_to_interaction(partner, interaction)
      import interaction_partner_t
      import interaction_t
      class(interaction_partner_t),     intent(inout) :: partner     !< the current interacion partner (this)
      class(interaction_t),             intent(inout) :: interaction !< the interaction
    end subroutine interaction_partner_copy_quantities_to_interaction


  end interface

  !> @brief the list of partners
  !!
  !! This class extends the list to create a partner list
  type, extends(linked_list_t) :: partner_list_t
    private
  contains
    procedure :: add => partner_list_add_node !< @copydoc interaction_partner_oct_m::partner_list_add_node
  end type partner_list_t

  !> @brief iterator for the list of partners
  !!
  !! This class extends the list iterator to create a partner list iterator
  type, extends(linked_list_iterator_t) :: partner_iterator_t
    private
  contains
    procedure :: get_next => partner_iterator_get_next !< @copydoc interaction_partner_oct_m::partner_iterator_get_next
  end type partner_iterator_t

contains

  ! ---------------------------------------------------------
  !> @brief add a partner to the list
  !!
  subroutine partner_list_add_node(this, partner)
    class(partner_list_t)                :: this     !< the partner list
    class(interaction_partner_t), target :: partner  !< the partner to add

    PUSH_SUB(partner_list_add_node)

    call this%add_ptr(partner)

    POP_SUB(partner_list_add_node)
  end subroutine partner_list_add_node

  ! ---------------------------------------------------------
  !> @brief get next partner from the list
  !!
  function partner_iterator_get_next(this) result(partner)
    class(partner_iterator_t),    intent(inout) :: this     !< the partner list
    class(interaction_partner_t), pointer       :: partner  !< the next element of the list

    PUSH_SUB(partner_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (interaction_partner_t)
      partner => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(partner_iterator_get_next)
  end function partner_iterator_get_next

end module interaction_partner_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
