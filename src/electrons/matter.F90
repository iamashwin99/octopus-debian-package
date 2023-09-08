!! Copyright (C) 2021 Micael Oliveira
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

!> This module defines a container system for electrons and ions
!!
module matter_oct_m
  use debug_oct_m
  use electrons_oct_m
  use global_oct_m
  use ions_oct_m
  use messages_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use profiling_oct_m
  implicit none

  private
  public ::               &
    matter_t

  !> @brief container class for for electrons and ions
  !!
  !! Currently, this system type is not yet used, as the electrons_oct_m::electrons_t type still contains the ions.
  !! Once the ions are fully disentangled from the electrons, this class will be used for electron-ion systems.
  !!
  type, extends(multisystem_t) :: matter_t
    class(electrons_t), pointer :: electrons => NULL()
    class(ions_t),      pointer :: ions => NULL()
  contains
    final :: matter_finalizer
  end type matter_t

  interface matter_t
    procedure matter_constructor
  end interface matter_t

contains

  ! ---------------------------------------------------------------------------------------
  function matter_constructor(namespace) result(matter)
    type(namespace_t), intent(in) :: namespace
    class(matter_t),   pointer    :: matter

    PUSH_SUB(matter_constructor)

    SAFE_ALLOCATE(matter)

    matter%namespace = namespace

    matter%ions => ions_t(namespace_t("ions", parent=matter%namespace))
    matter%electrons => electrons_t(namespace_t("electrons", parent=matter%namespace))

    call matter%list%add(matter%ions)
    call matter%list%add(matter%electrons)

    POP_SUB(matter_constructor)
  end function matter_constructor

  ! ---------------------------------------------------------
  subroutine matter_finalizer(this)
    type(matter_t), intent(inout) :: this

    PUSH_SUB(matter_finalizer)

    call multisystem_end(this)
    nullify(this%electrons)
    nullify(this%ions)

    POP_SUB(matter_finalizer)
  end subroutine matter_finalizer

end module matter_oct_m
