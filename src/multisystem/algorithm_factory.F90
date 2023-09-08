!! Copyright (C) 2023 M. Oliveira
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

module algorithm_factory_oct_m
  use algorithm_oct_m
  use interaction_partner_oct_m
  implicit none

  private
  public ::                        &
    algorithm_factory_t

  !> @brief Abstract class for the algoritm factories
  !!
  !! The abstract algorithm factory defines the abstract interface for the algorithm factories and needs
  !! to be specialized for actual use.
  type, abstract :: algorithm_factory_t
  contains
    procedure(algorithm_factory_create),        deferred :: create
    procedure(algorithm_factory_create_static), deferred :: create_static
  end type algorithm_factory_t

  abstract interface
    function algorithm_factory_create(this, system) result(algorithm)
      import :: algorithm_factory_t
      import interaction_partner_t
      import algorithm_t
      class(algorithm_factory_t),   intent(in)          :: this
      class(interaction_partner_t), intent(in),  target :: system
      class(algorithm_t),           pointer             :: algorithm
    end function algorithm_factory_create

    function algorithm_factory_create_static(this, system) result(algorithm)
      import :: algorithm_factory_t
      import interaction_partner_t
      import algorithm_t
      class(algorithm_factory_t),   intent(in)          :: this
      class(interaction_partner_t), intent(in),  target :: system
      class(algorithm_t),           pointer             :: algorithm
    end function algorithm_factory_create_static
  end interface

end module algorithm_factory_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
