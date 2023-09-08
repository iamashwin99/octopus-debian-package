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

module system_factory_oct_m
  use charged_particle_oct_m
  use classical_particle_oct_m
  use debug_oct_m
  use dftb_oct_m
  use dispersive_medium_oct_m
  use global_oct_m
  use linear_medium_oct_m
  use matter_oct_m
  use maxwell_oct_m
  use messages_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use electrons_oct_m
  use system_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::                         &
    system_factory_t

  !# doc_start system_types
  integer, parameter ::             &
    SYSTEM_ELECTRONIC         = 1,  & !< electronic system (electrons_oct_m::electrons_t)
    SYSTEM_MAXWELL            = 2,  & !< maxwell system, (maxwell_oct_m::maxwell_t)
    SYSTEM_CLASSICAL_PARTICLE = 3,  & !< single classical particle (classical_particle_oct_m::classical_particle_t)
    SYSTEM_CHARGED_PARTICLE   = 4,  & !< single charged classical particle (charged_particle_oct_m::charged_particle_t)
    SYSTEM_DFTBPLUS           = 5,  & !< tight binding system (dftb_oct_m::dftb_t)
    SYSTEM_LINEAR_MEDIUM      = 6,  & !< linear medium for Maxwell calculations (linear_medium_oct_m::linear_medium_t)
    SYSTEM_MATTER             = 7,  & !< electrons including ions (matter_oct_m::matter_t)
    SYSTEM_DISPERSIVE_MEDIUM  = 8,  & !< dispersive medium for classical electrodynamics (dispersive_medium_oct_m::dispersive_medium_t)
    SYSTEM_MULTISYSTEM        = 9     !< container system. (multisystem_basic_oct_m::multisystem_basic_t)
  !# doc_end

  !> @brief factory for classes, derived from the abstract system_cot_m::system_t class
  !!
  !! The systems factory performs the book-keeping when creating systems.
  !! Systems are actually created by their own ''create'' function, to which
  !! an instance of the systems_factory_t has to be passed.
  type, extends(system_factory_abst_t) :: system_factory_t
  contains
    procedure :: create => system_factory_create         !< @copydoc system_factory_oct_m::system_factory_create
    procedure :: block_name => system_factory_block_name !< @copydoc system_factory_oct_m::system_factory_block_name
  end type system_factory_t

contains

  ! ---------------------------------------------------------------------------------------
  !> @brief create a new system.
  !!
  !! This function is called by multisystem_create_system() and run().
  !! It calls the constructors of specific systems.
  recursive function system_factory_create(this, namespace, name, type) result(system)
    class(system_factory_t), intent(in) :: this       !< the system factory
    type(namespace_t),       intent(in) :: namespace  !< namespace of the parent
    character(len=*),        intent(in) :: name       !< name of the system (will define namespace)
    integer,                 intent(in) :: type       !< type of the system to create
    class(system_t),         pointer    :: system     !< pointer to newly created system

    PUSH_SUB(system_factory_create)

    !%Variable Systems
    !%Type block
    !%Section System
    !%Description
    !% List of systems that will be treated in the calculation.
    !% The first column should be a string containing the system name.
    !% The second column should be the system type. See below for a list of
    !% available system types.
    !%Option electronic 1
    !% An electronic system. (not fully implemented yet)
    !%Option maxwell 2
    !% A maxwell system.
    !%Option classical_particle 3
    !% A classical particle. Used for testing purposes only.
    !%Option charged_particle 4
    !% A charged classical particle.
    !%Option dftbplus 5
    !% A DFTB+ system
    !%Option linear_medium 6
    !% A linear medium for classical electrodynamics.
    !%Option matter 7
    !% A matter system containing electrons and classical ions.
    !%Option dispersive_medium 8
    !% (Experimental) A dispersive medium for classical electrodynamics.
    !%Option multisystem 9
    !% A system containing other systems.
    !%End
    select case (type)
    case (SYSTEM_MULTISYSTEM)
      system => multisystem_basic_t(namespace_t(name, parent=namespace), this)
    case (SYSTEM_ELECTRONIC)
      system => electrons_t(namespace_t(name, parent=namespace))
    case (SYSTEM_MAXWELL)
      system => maxwell_t(namespace_t(name, parent=namespace))
    case (SYSTEM_CLASSICAL_PARTICLE)
      system => classical_particle_t(namespace_t(name, parent=namespace))
    case (SYSTEM_CHARGED_PARTICLE)
      system => charged_particle_t(namespace_t(name, parent=namespace))
    case (SYSTEM_DFTBPLUS)
      system => dftb_t(namespace_t(name, parent=namespace))
    case (SYSTEM_LINEAR_MEDIUM)
      system => linear_medium_t(namespace_t(name, parent=namespace))
    case (SYSTEM_MATTER)
      system => matter_t(namespace_t(name, parent=namespace))
    case (SYSTEM_DISPERSIVE_MEDIUM)
      system => dispersive_medium_t(namespace_t(name, parent=namespace))
      call messages_experimental('dispersive_medium', namespace=namespace)
    case default
      system => null()
    end select

    POP_SUB(system_factory_create)
  end function system_factory_create

  ! ---------------------------------------------------------------------------------------
  !> @brief block the name for new systems to avoid clashes.
  character(len=80) function system_factory_block_name(this) result(name)
    class(system_factory_t), intent(in)    :: this

    PUSH_SUB(system_factory_block_name)

    name = "Systems"

    POP_SUB(system_factory_block_name)
  end function system_factory_block_name

end module system_factory_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
