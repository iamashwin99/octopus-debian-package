!! Copyright (C) 2019 N. Tancogne-Dejean
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

module classical_particle_oct_m
  use algorithm_oct_m
  use classical_particles_oct_m
  use clock_oct_m
  use debug_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use gravity_oct_m
  use io_oct_m
  use iso_c_binding
  use lennard_jones_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::                  &
    classical_particle_t,    &
    classical_particle_init

  integer, parameter ::     &
    OUTPUT_COORDINATES = 1, &
    OUTPUT_ENERGY      = 2

  !> @brief class for a neutral classical particle
  !!
  type, extends(classical_particles_t) :: classical_particle_t
    type(c_ptr) :: output_handle(2)

  contains
    procedure :: init_interaction => classical_particle_init_interaction
    procedure :: initial_conditions => classical_particle_initial_conditions
    procedure :: output_start => classical_particle_output_start
    procedure :: output_write => classical_particle_output_write
    procedure :: output_finish => classical_particle_output_finish
    procedure :: restart_write_data => classical_particle_restart_write_data
    procedure :: update_quantity => classical_particle_update_quantity
    procedure :: update_exposed_quantity => classical_particle_update_exposed_quantity
    procedure :: init_interaction_as_partner => classical_particle_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => classical_particle_copy_quantities_to_interaction
    final :: classical_particle_finalize
  end type classical_particle_t

  interface classical_particle_t
    procedure classical_particle_constructor
  end interface classical_particle_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function classical_particle_constructor(namespace) result(sys)
    class(classical_particle_t), pointer    :: sys
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(classical_particle_constructor)

    SAFE_ALLOCATE(sys)

    call classical_particle_init(sys, namespace)

    POP_SUB(classical_particle_constructor)
  end function classical_particle_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine classical_particle_init(this, namespace)
    class(classical_particle_t), intent(inout) :: this
    type(namespace_t),           intent(in)    :: namespace

    PUSH_SUB(classical_particle_init)

    this%namespace = namespace

    call space_init(this%space, namespace)
    if (this%space%is_periodic()) then
      call messages_not_implemented('Classical particle for periodic systems', namespace=namespace)
    end if

    call messages_print_with_emphasis(msg="Classical Particle", namespace=namespace)

    call classical_particles_init(this, 1)

    !%Variable ParticleMass
    !%Type float
    !%Section ClassicalParticles
    !%Description
    !% Mass of classical particle in Kg.
    !%End
    call parse_variable(namespace, 'ParticleMass', M_ONE, this%mass(1))
    call messages_print_var_value('ParticleMass', this%mass(1), namespace=namespace)

    call this%supported_interactions%add(GRAVITY)
    call this%supported_interactions_as_partner%add(GRAVITY)
    call this%supported_interactions%add(LENNARD_JONES)
    call this%supported_interactions_as_partner%add(LENNARD_JONES)

    !%Variable LennardJonesEpsilon
    !%Type float
    !%Section ClassicalParticles
    !%Description
    !% Epsilon parameter (dispersion energy) of Lennard-Jones interaction for this species.
    !% In case two particles have a different epsilon, the combination rule will be computed
    !% <math>\epsilon_{12} = \sqrt{\epsilon_1 + \epsilon_2}</math>.
    !%End
    call parse_variable(namespace, 'LennardJonesEpsilon', M_ONE, this%lj_epsilon(1))
    call messages_print_var_value('LennardJonesEpsilon', this%lj_epsilon(1), namespace=namespace)

    !%Variable LennardJonesSigma
    !%Type float
    !%Section ClassicalParticles
    !%Description
    !% Sigma parameter (particle size) of Lennard-Jones interaction for this species.
    !% In case two particles have a different sigma, the combination rule will be computed
    !% <math>\sigma_{12} = (\sigma_1 + \sigma_2) / 2.
    !%End
    call parse_variable(namespace, 'LennardJonesSigma', M_ONE, this%lj_sigma(1))
    call messages_print_var_value('LennardJonesSigma', this%lj_sigma(1), namespace=namespace)

    call messages_print_with_emphasis(namespace=namespace)

    POP_SUB(classical_particle_init)
  end subroutine classical_particle_init

  ! ---------------------------------------------------------
  subroutine classical_particle_init_interaction(this, interaction)
    class(classical_particle_t), target, intent(inout) :: this
    class(interaction_t),                intent(inout) :: interaction

    PUSH_SUB(classical_particle_init_interaction)

    select type (interaction)
    type is (gravity_t)
      call interaction%init(this%space%dim, 1, this%quantities, this%mass, this%pos)
    type is (lennard_jones_t)
      if (.not. (parse_is_defined(this%namespace, 'LennardJonesSigma') .and. &
        parse_is_defined(this%namespace, 'LennardJonesEpsilon') )) then
        write(message(1),'(a,es9.2)') 'Using default value for Lennard-Jones parameter.'
        call messages_warning(1, namespace=this%namespace)
      end if

      call interaction%init(this%space%dim, 1, this%quantities, this%pos, this%lj_epsilon(1), &
        this%lj_sigma(1))
    class default
      call classical_particles_init_interaction(this, interaction)
    end select

    POP_SUB(classical_particle_init_interaction)
  end subroutine classical_particle_init_interaction

  ! ---------------------------------------------------------
  subroutine classical_particle_initial_conditions(this)
    class(classical_particle_t), intent(inout) :: this

    integer :: n_rows, idir
    type(block_t) :: blk

    PUSH_SUB(classical_particle_initial_conditions)

    !%Variable ParticleInitialPosition
    !%Type block
    !%Section ClassicalParticles
    !%Description
    !% Initial position of classical particle, in Km.
    !%End
    this%pos = M_ZERO
    if (parse_block(this%namespace, 'ParticleInitialPosition', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call messages_input_error(this%namespace, 'ParticleInitialPosition')

      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%pos(idir, 1))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value('ParticleInitialPosition', this%pos(1:this%space%dim, 1), namespace=this%namespace)

    !%Variable ParticleInitialVelocity
    !%Type block
    !%Section ClassicalParticles
    !%Description
    !% Initial velocity of classical particle in Km/s.
    !%End
    this%vel = M_ZERO
    if (parse_block(this%namespace, 'ParticleInitialVelocity', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call messages_input_error(this%namespace, 'ParticleInitialVelocity')
      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%vel(idir, 1))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value('ParticleInitialVelocity', this%vel(1:this%space%dim, 1), namespace=this%namespace)

    POP_SUB(classical_particle_initial_conditions)
  end subroutine classical_particle_initial_conditions

  ! ---------------------------------------------------------
  subroutine classical_particle_output_start(this)
    class(classical_particle_t), intent(inout) :: this

    integer :: iteration

    PUSH_SUB(classical_particle_output_start)

    iteration = this%clock%get_tick()
    if (iteration > 0) then
      ! add one if restarting because the output routine is only called at the end of the timestep
      iteration = iteration + CLOCK_TICK
    end if
    ! Create output handle
    call io_mkdir('td.general', this%namespace)
    if (mpi_grp_is_root(this%grp)) then
      call write_iter_init(this%output_handle(OUTPUT_COORDINATES), iteration, this%algo%dt, &
        trim(io_workpath("td.general/coordinates", this%namespace)))
      call write_iter_init(this%output_handle(OUTPUT_ENERGY), iteration, this%algo%dt, &
        trim(io_workpath("td.general/energy", this%namespace)))
    end if

    ! Output info for first iteration
    if (iteration == 0) then
      call this%output_write()
    end if

    POP_SUB(classical_particle_output_start)
  end subroutine classical_particle_output_start

  ! ---------------------------------------------------------
  subroutine classical_particle_output_finish(this)
    class(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_output_finish)

    if (mpi_grp_is_root(this%grp)) then
      call write_iter_end(this%output_handle(OUTPUT_COORDINATES))
      call write_iter_end(this%output_handle(OUTPUT_ENERGY))
    end if

    POP_SUB(classical_particle_output_finish)
  end subroutine classical_particle_output_finish

  ! ---------------------------------------------------------
  subroutine classical_particle_output_write(this)
    class(classical_particle_t), intent(inout) :: this

    integer :: idir, ii
    character(len=50) :: aux
    FLOAT :: tmp(this%space%dim)

    if (.not. mpi_grp_is_root(this%grp)) return ! only first node outputs

    PUSH_SUB(classical_particle_output_write)

    if (this%clock%get_tick() == 0) then
      ! header
      call write_iter_clear(this%output_handle(OUTPUT_COORDINATES))
      call write_iter_string(this%output_handle(OUTPUT_COORDINATES), &
        '############################################################################')
      call write_iter_nl(this%output_handle(OUTPUT_COORDINATES))
      call write_iter_string(this%output_handle(OUTPUT_COORDINATES),'# HEADER')
      call write_iter_nl(this%output_handle(OUTPUT_COORDINATES))

      ! first line: column names
      call write_iter_header_start(this%output_handle(OUTPUT_COORDINATES))
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'x(', idir, ')'
        call write_iter_header(this%output_handle(OUTPUT_COORDINATES), aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'v(', idir, ')'
        call write_iter_header(this%output_handle(OUTPUT_COORDINATES), aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'f(', idir, ')'
        call write_iter_header(this%output_handle(OUTPUT_COORDINATES), aux)
      end do
      call write_iter_nl(this%output_handle(OUTPUT_COORDINATES))

      ! second line: units
      call write_iter_string(this%output_handle(OUTPUT_COORDINATES), '#[Iter n.]')
      call write_iter_header(this%output_handle(OUTPUT_COORDINATES), '[' // trim(units_abbrev(units_out%time)) // ']')
      do idir = 1, this%space%dim
        call write_iter_header(this%output_handle(OUTPUT_COORDINATES), '[' // trim(units_abbrev(units_out%length)) // ']')
      end do
      do idir = 1, this%space%dim
        call write_iter_header(this%output_handle(OUTPUT_COORDINATES), '[' // trim(units_abbrev(units_out%velocity)) // ']')
      end do
      do idir = 1, this%space%dim
        call write_iter_header(this%output_handle(OUTPUT_COORDINATES), '[' // trim(units_abbrev(units_out%force)) // ']')
      end do
      call write_iter_nl(this%output_handle(OUTPUT_COORDINATES))

      call write_iter_string(this%output_handle(OUTPUT_COORDINATES), &
        '############################################################################')
      call write_iter_nl(this%output_handle(OUTPUT_COORDINATES))
    end if

    call write_iter_start(this%output_handle(OUTPUT_COORDINATES))

    ! Position
    tmp(1:this%space%dim) = units_from_atomic(units_out%length, this%pos(1:this%space%dim, 1))
    call write_iter_double(this%output_handle(OUTPUT_COORDINATES), tmp, this%space%dim)
    ! Velocity
    tmp(1:this%space%dim) = units_from_atomic(units_out%velocity, this%vel(1:this%space%dim, 1))
    call write_iter_double(this%output_handle(OUTPUT_COORDINATES), tmp, this%space%dim)
    ! Force
    tmp(1:this%space%dim) = units_from_atomic(units_out%force, this%tot_force(1:this%space%dim, 1))
    call write_iter_double(this%output_handle(OUTPUT_COORDINATES), tmp, this%space%dim)
    call write_iter_nl(this%output_handle(OUTPUT_COORDINATES))

    ! Energies
    if (this%clock%get_tick() == 0) then
      ! header
      call write_iter_clear(this%output_handle(OUTPUT_ENERGY))
      call write_iter_string(this%output_handle(OUTPUT_ENERGY), &
        '############################################################################')
      call write_iter_nl(this%output_handle(OUTPUT_ENERGY))
      call write_iter_string(this%output_handle(OUTPUT_ENERGY),'# HEADER')
      call write_iter_nl(this%output_handle(OUTPUT_ENERGY))

      ! first line: column names
      call write_iter_header_start(this%output_handle(OUTPUT_ENERGY))
      call write_iter_header(this%output_handle(OUTPUT_ENERGY), 'Total')
      call write_iter_header(this%output_handle(OUTPUT_ENERGY), 'Kinetic')
      call write_iter_header(this%output_handle(OUTPUT_ENERGY), 'Potential')
      call write_iter_header(this%output_handle(OUTPUT_ENERGY), 'Internal')
      call write_iter_nl(this%output_handle(OUTPUT_ENERGY))

      ! second line: units
      call write_iter_string(this%output_handle(OUTPUT_ENERGY), '#[Iter n.]')
      call write_iter_header(this%output_handle(OUTPUT_ENERGY), '[' // trim(units_abbrev(units_out%time)) // ']')
      do ii = 1, 4
        call write_iter_header(this%output_handle(OUTPUT_ENERGY), '[' // trim(units_abbrev(units_out%energy)) // ']')
      end do
      call write_iter_nl(this%output_handle(OUTPUT_ENERGY))

      call write_iter_string(this%output_handle(OUTPUT_ENERGY), &
        '############################################################################')
      call write_iter_nl(this%output_handle(OUTPUT_ENERGY))
    end if
    call write_iter_start(this%output_handle(OUTPUT_ENERGY))

    call write_iter_double(this%output_handle(OUTPUT_ENERGY), units_from_atomic(units_out%energy, this%total_energy), 1)
    call write_iter_double(this%output_handle(OUTPUT_ENERGY), units_from_atomic(units_out%energy, this%kinetic_energy), 1)
    call write_iter_double(this%output_handle(OUTPUT_ENERGY), units_from_atomic(units_out%energy, this%potential_energy), 1)
    call write_iter_double(this%output_handle(OUTPUT_ENERGY), units_from_atomic(units_out%energy, this%internal_energy), 1)
    call write_iter_nl(this%output_handle(OUTPUT_ENERGY))

    POP_SUB(classical_particle_output_write)
  end subroutine classical_particle_output_write

  ! ---------------------------------------------------------
  subroutine classical_particle_restart_write_data(this)
    class(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_restart_write_data)

    call classical_particles_restart_write_data(this)

    if (mpi_grp_is_root(this%grp)) then
      call write_iter_flush(this%output_handle(OUTPUT_COORDINATES))
      call write_iter_flush(this%output_handle(OUTPUT_ENERGY))
    end if

    POP_SUB(classical_particle_restart_write_data)
  end subroutine classical_particle_restart_write_data


  ! ---------------------------------------------------------
  subroutine classical_particle_update_quantity(this, iq)
    class(classical_particle_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particle_update_quantity)

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(this%quantities(iq)%updated_on_demand)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_quantity(this, iq)
    end select

    POP_SUB(classical_particle_update_quantity)
  end subroutine classical_particle_update_quantity

  ! ---------------------------------------------------------
  subroutine classical_particle_update_exposed_quantity(partner, iq)
    class(classical_particle_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particle_update_exposed_quantity)

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(partner%quantities(iq)%updated_on_demand)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_exposed_quantity(partner, iq)
    end select

    POP_SUB(classical_particle_update_exposed_quantity)
  end subroutine classical_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine classical_particle_init_interaction_as_partner(partner, interaction)
    class(classical_particle_t), intent(in)    :: partner
    class(interaction_t),        intent(inout) :: interaction

    PUSH_SUB(classical_particle_init_interaction_as_partner)

    select type (interaction)
    type is (gravity_t)
      interaction%partner_np = 1
      SAFE_ALLOCATE(interaction%partner_mass(1))
      SAFE_ALLOCATE(interaction%partner_pos(1:partner%space%dim, 1))

    type is (lennard_jones_t)
      interaction%partner_np = 1
      SAFE_ALLOCATE(interaction%partner_pos(1:partner%space%dim, 1))
      ! in case the LennardJones epsilon and sigma of system and partner are different,
      ! we compute the combination rules with arithmetic and geometric means:
      ! (they give back the original parameters if they happen to be equal):
      ! sigma_12 = (sigma_1 + sigma_2)/2,    epsilon_12 = sqrt(epsilon_1 * epsilon_2)
      interaction%lj_sigma = M_HALF * (partner%lj_sigma(1) + interaction%lj_sigma)
      interaction%lj_epsilon = sqrt(partner%lj_epsilon(1) * interaction%lj_epsilon)

    class default
      ! Other interactions should be handled by the parent class
      call classical_particles_init_interaction_as_partner(partner, interaction)
    end select

    POP_SUB(classical_particle_init_interaction_as_partner)
  end subroutine classical_particle_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine classical_particle_copy_quantities_to_interaction(partner, interaction)
    class(classical_particle_t),          intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(classical_particle_copy_quantities_to_interaction)

    select type (interaction)
    type is (gravity_t)
      interaction%partner_mass(1) = partner%mass(1)
      interaction%partner_pos(:,1) = partner%pos(:,1)

    type is (lennard_jones_t)
      interaction%partner_pos(:,1) = partner%pos(:,1)

    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(classical_particle_copy_quantities_to_interaction)
  end subroutine classical_particle_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine classical_particle_finalize(this)
    type(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_finalize)

    call classical_particles_end(this)

    POP_SUB(classical_particle_finalize)
  end subroutine classical_particle_finalize

end module classical_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
