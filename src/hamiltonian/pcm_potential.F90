!! Copyright (C) 2014 Alain Delgado Gran, Carlo Andrea Rozzi, Stefano Corni, Gabriel Gil
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

module pcm_potential_oct_m
  use box_minimum_oct_m
  use comm_oct_m
  use debug_oct_m
  use ext_partner_list_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use index_oct_m
  use interaction_partner_oct_m
  use io_oct_m
  use ions_oct_m
  use kick_oct_m
  use lasers_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mesh_interpolation_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use pcm_oct_m
  use pcm_eom_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use species_oct_m
  use space_oct_m
  use varinfo_oct_m

  ! to output debug info
  use unit_oct_m
  use unit_system_oct_m
  use io_function_oct_m
  use mesh_function_oct_m

  implicit none

  private

  public ::                 &
    pcm_hartree_potential

contains

  ! -----------------------------------------------------------------------------
  !> PCM reaction field due to the electronic density
  subroutine pcm_hartree_potential(pcm, space, mesh, psolver, ext_partners, vhartree, density, pcm_corr, &
    kick, time)
    type(pcm_t),         intent(inout) :: pcm
    type(space_t),       intent(in)    :: space
    class(mesh_t),       intent(in)    :: mesh
    type(poisson_t),     intent(inout) :: psolver
    type(partner_list_t),intent(in)    :: ext_partners
    FLOAT,               intent(in)    :: vhartree(:)
    FLOAT,               intent(in)    :: density(:)
    FLOAT,               intent(out)   :: pcm_corr
    type(kick_t), optional, intent(in) :: kick
    FLOAT, optional,     intent(in)    :: time

    FLOAT, allocatable :: potx(:)
    CMPLX, allocatable :: kick_eval(:)
    FLOAT, allocatable :: kick_real(:)
    integer :: ii

    logical :: kick_time
    type(lasers_t), pointer :: lasers

    PUSH_SUB(pcm_hartree_potential)

    if (.not. pcm%run_pcm .or. .not. pcm_update(pcm)) then
      pcm_corr = M_ZERO
      POP_SUB(pcm_hartree_potential)
      return
    end if

    !> Generates the real-space PCM potential due to electrons during the SCF calculation.
    if (pcm%solute) then
      call pcm_calc_pot_rs(pcm, mesh, psolver, v_h = vhartree, time_present = present(time))
    end if

    !> Local field effects due to the applied electrostatic potential representing the laser and the kick (if they were).
    !! For the laser, the latter is only valid in the long-wavelength limit.
    !! Static potentials are included in subroutine hamiltonian_elec_epot_generate (module hamiltonian).
    !! The sign convention for typical potentials and kick are different...
    if (pcm%localf .and. present(time)) then
      lasers => list_get_lasers(ext_partners)
      if (associated(lasers) .and. present(kick)) then !< external potential and kick
        SAFE_ALLOCATE(potx(1:mesh%np_part))
        SAFE_ALLOCATE(kick_eval(1:mesh%np_part))
        SAFE_ALLOCATE(kick_real(1:mesh%np_part))
        potx = M_ZERO
        kick_eval = M_ZERO
        do ii = 1, lasers%no_lasers
          call laser_potential(lasers%lasers(ii), mesh, potx, time)
        end do
        kick_real = M_ZERO
        kick_time = ((pcm%iter-1)*pcm%dt <= kick%time) .and. (pcm%iter*pcm%dt > kick%time)
        if (kick_time) then
          call kick_function_get(space, mesh, kick, kick_eval, 1, to_interpolate = .true.)
          kick_eval = kick%delta_strength * kick_eval
          kick_real = TOFLOAT(kick_eval)
        end if
        call pcm_calc_pot_rs(pcm, mesh, psolver, v_ext = potx, kick = -kick_real, &
          time_present = present(time), kick_time = kick_time)
        SAFE_DEALLOCATE_A(potx)
        SAFE_DEALLOCATE_A(kick_eval)
        SAFE_DEALLOCATE_A(kick_real)
      else if (associated(lasers) .and. .not.present(kick)) then !< just external potential
        SAFE_ALLOCATE(potx(1:mesh%np_part))
        potx = M_ZERO
        do ii = 1, lasers%no_lasers
          call laser_potential(lasers%lasers(ii), mesh, potx, time)
        end do
        call pcm_calc_pot_rs(pcm, mesh, psolver, v_ext = potx, time_present = present(time))
        SAFE_DEALLOCATE_A(potx)
      else if (.not.associated(lasers) .and. present(kick)) then !< just kick
        SAFE_ALLOCATE(kick_eval(1:mesh%np_part))
        SAFE_ALLOCATE(kick_real(1:mesh%np_part))
        kick_eval = M_ZERO
        kick_real = M_ZERO
        kick_time =((pcm%iter-1)*pcm%dt <= kick%time) .and. (pcm%iter*pcm%dt > kick%time)
        if (kick_time) then
          call kick_function_get(space, mesh, kick, kick_eval, 1, to_interpolate = .true.)
          kick_eval = kick%delta_strength * kick_eval
          kick_real = TOFLOAT(kick_eval)
        end if
        call pcm_calc_pot_rs(pcm, mesh, psolver, kick = -kick_real, &
          time_present = present(time), kick_time = kick_time)
        SAFE_DEALLOCATE_A(kick_eval)
        SAFE_DEALLOCATE_A(kick_real)
      end if

      ! Calculating the PCM term renormalizing the sum of the single-particle energies
      ! to keep the idea of pcm_corr... but it will be added later on
      pcm_corr = dmf_dotp( mesh, density, pcm%v_e_rs + pcm%v_n_rs + pcm%v_ext_rs)
    else
      ! Calculating the PCM term renormalizing the sum of the single-particle energies
      pcm_corr = dmf_dotp( mesh, density, pcm%v_e_rs + pcm%v_n_rs)
    end if

    if (debug%info) then
      call messages_write(' PCM potential updated')
      call messages_new_line()
      call messages_write(' PCM update iteration counter: ')
      call messages_write(pcm%iter)
      call messages_info()
    end if

    POP_SUB(pcm_hartree_potential)

  end subroutine pcm_hartree_potential

end module pcm_potential_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
