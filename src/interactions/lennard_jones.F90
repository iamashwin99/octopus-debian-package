!! Copyright (C) 2023 F. Bonafe
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

module lennard_jones_oct_m
  use debug_oct_m
  use force_interaction_oct_m
  use global_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quantity_oct_m

  implicit none

  private
  public ::                &
    lennard_jones_t

  !> Lennard-Jones interaction between two systems of particles.
  type, extends(force_interaction_t) :: lennard_jones_t
    private
    FLOAT, pointer :: system_pos(:,:) !< pointer to array storing the positions of the particles
    FLOAT, public  :: lj_epsilon
    FLOAT, public  :: lj_sigma

    integer, public :: partner_np = 0 !< number of particles in the partner system
    FLOAT, allocatable, public :: partner_pos(:,:) !< array storing a copy of the positions of the partner particles

  contains
    procedure :: init => lennard_jones_init
    procedure :: calculate => lennard_jones_calculate
    procedure :: calculate_energy => lennard_jones_calculate_energy
    final :: lennard_jones_finalize
  end type lennard_jones_t

  interface lennard_jones_t
    module procedure lennard_jones_constructor
  end interface lennard_jones_t

contains

  ! ---------------------------------------------------------
  function lennard_jones_constructor(partner, intra_interaction) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(lennard_jones_t),               pointer       :: this
    logical,                              intent(in)    :: intra_interaction

    PUSH_SUB(lennard_jones_constructor)

    SAFE_ALLOCATE(this)

    this%label = "lennard_jones"

    this%partner => partner

    ! The Lennard Jones interaction needs the position of the particles and the parameters
    this%n_system_quantities = 1
    SAFE_ALLOCATE(this%system_quantities(1:this%n_system_quantities))
    this%system_quantities(1) = POSITION
    nullify(this%system_pos)

    ! The Lennard Jones interaction needs the position of the partner
    this%n_partner_quantities = 1
    SAFE_ALLOCATE(this%partner_quantities(1:this%n_partner_quantities))
    this%partner_quantities(1) = POSITION
    this%intra_interaction = intra_interaction

    POP_SUB(lennard_jones_constructor)
  end function lennard_jones_constructor

  ! ---------------------------------------------------------
  subroutine lennard_jones_init(this, dim, system_np, system_quantities, system_pos, system_eps, system_sigma)
    class(lennard_jones_t),             intent(inout) :: this
    integer,                            intent(in)    :: dim !< number of dimensions in space
    integer,                            intent(in)    :: system_np  !< number of particles in the system that owns this interaction
    type(quantity_t),                   intent(inout) :: system_quantities(:)
    FLOAT,                      target, intent(in)    :: system_pos(:,:)
    FLOAT,                              intent(in)    :: system_eps
    FLOAT,                              intent(in)    :: system_sigma

    PUSH_SUB(lennard_jones_init)

    this%dim = dim
    this%system_np = system_np
    SAFE_ALLOCATE(this%force(1:dim, 1:system_np))

    this%lj_epsilon = system_eps
    this%lj_sigma = system_sigma

    this%system_pos => system_pos

    POP_SUB(lennard_jones_init)
  end subroutine lennard_jones_init

  ! ---------------------------------------------------------
  subroutine lennard_jones_calculate(this)
    class(lennard_jones_t),             intent(inout) :: this

    integer :: ip, jp
    FLOAT :: dist, rr(1:this%dim), lj_force

    PUSH_SUB(lennard_jones_calculate)

    ASSERT(allocated(this%partner_pos))

    do ip = 1, this%system_np
      do jp = 1, this%partner_np
        if (this%intra_interaction .and. ip == jp ) cycle

        ! r_ij = r_i - r_j
        rr(1:this%dim) = this%system_pos(1:this%dim, ip) - this%partner_pos(1:this%dim, jp)
        dist = sqrt(sum(rr(1:this%dim)**2))

        ! lj_force = -d U(|r_ij|) / d |r_ij|
        lj_force = CNST(48.0) * this%lj_epsilon * (this%lj_sigma**12 / dist**13 - &
          M_HALF * this%lj_sigma**6 / dist**7)

        ! F_i = - d U(r) / d r_i = - (d U (r) / d |r_ij|) * (d r_ij / d r_i) =
        ! = - d U (r) / d r_ij, because d r_ij / d r_i = 1
        this%force(1:this%dim, ip) = rr(1:this%dim) / dist * lj_force
      end do
    end do

    POP_SUB(lennard_jones_calculate)
  end subroutine lennard_jones_calculate

  ! ---------------------------------------------------------
  subroutine lennard_jones_calculate_energy(this)
    class(lennard_jones_t),             intent(inout) :: this

    integer :: ip, jp
    FLOAT :: dist

    PUSH_SUB(lennard_jones_calculate_energy)

    ASSERT(allocated(this%partner_pos))

    this%energy = M_ZERO
    do ip = 1, this%system_np
      do jp = 1, this%partner_np
        if (this%intra_interaction .and. ip == jp ) cycle

        dist = sqrt(sum((this%system_pos(1:this%dim, ip) - this%partner_pos(1:this%dim, jp))**2)) + M_EPSILON

        ! this%energy = this%energy + 0.5 * U(r_ij)
        this%energy = this%energy + M_TWO * this%lj_epsilon * ( (this%lj_sigma / dist)**12 - &
          (this%lj_sigma / dist)**6 )
      end do
    end do

    POP_SUB(lennard_jones_calculate_energy)
  end subroutine lennard_jones_calculate_energy


  ! ---------------------------------------------------------
  subroutine lennard_jones_finalize(this)
    type(lennard_jones_t), intent(inout) :: this

    PUSH_SUB(lennard_jones_finalize)

    this%force = M_ZERO
    nullify(this%system_pos)
    SAFE_DEALLOCATE_A(this%partner_pos)
    SAFE_DEALLOCATE_A(this%force)

    call interaction_with_partner_end(this)

    POP_SUB(lennard_jones_finalize)
  end subroutine lennard_jones_finalize

end module lennard_jones_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
