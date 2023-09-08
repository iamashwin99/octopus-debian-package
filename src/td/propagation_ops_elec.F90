!! Copyright (C) 2019 N. Tancogne-Dejean
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

module propagation_ops_elec_oct_m
  use accel_oct_m
  use batch_oct_m
  use debug_oct_m
  use density_oct_m
  use exponential_oct_m
  use ext_partner_list_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use interaction_partner_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_verlet_oct_m
  use space_oct_m
  use states_elec_oct_m
  use varinfo_oct_m
  use wfs_elec_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                             &
    propagation_ops_elec_t,                      &
    propagation_ops_elec_update_hamiltonian,     &
    propagation_ops_elec_exp_apply,              &
    propagation_ops_elec_fuse_density_exp_apply, &
    propagation_ops_elec_move_ions,              &
    propagation_ops_elec_restore_ions,           &
    propagation_ops_elec_propagate_gauge_field,  &
    propagation_ops_elec_restore_gauge_field,    &
    propagation_ops_elec_interpolate_get,        &
    propagation_ops_do_pack,                     &
    propagation_ops_do_unpack,                   &
    propagation_ops_finish_unpack

  type :: propagation_ops_elec_t
    private
    type(ion_state_t) :: ions_state
    FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)
  end type propagation_ops_elec_t


contains

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_update_hamiltonian(namespace, space, st, mesh, hm, ext_partners, time)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(states_elec_t),      intent(inout) :: st
    class(mesh_t),            intent(in)    :: mesh
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(partner_list_t),     intent(in)    :: ext_partners
    FLOAT,                    intent(in)    :: time

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_update_hamiltonian)

    call profiling_in(prof, 'ELEC_UPDATE_H')

    call calculate_mxll_dipole_field(hm, mesh, st)
    call hm%update(mesh, namespace, space, ext_partners, time = time)
    call lda_u_update_occ_matrices(hm%lda_u, namespace, mesh, st, hm%hm_base, hm%energy)

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_update_hamiltonian)
  end subroutine propagation_ops_elec_update_hamiltonian


  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_move_ions(wo, gr, hm, st, namespace, space, ions_dyn, ions, &
    ext_partners, time, dt, save_pos)
    class(propagation_ops_elec_t), intent(inout) :: wo
    type(grid_t),                  intent(in)    :: gr
    type(hamiltonian_elec_t),      intent(inout) :: hm
    type(states_elec_t),           intent(inout) :: st
    type(namespace_t),             intent(in)    :: namespace
    type(space_t),                 intent(in)    :: space
    type(ion_dynamics_t),          intent(inout) :: ions_dyn
    type(ions_t),                  intent(inout) :: ions
    type(partner_list_t),          intent(in)    :: ext_partners
    FLOAT,                         intent(in)    :: time
    FLOAT,                         intent(in)    :: dt
    logical, optional,             intent(in)    :: save_pos

    type(profile_t), save :: prof
    FLOAT :: dt_ions

    PUSH_SUB(propagation_ops_elec_move_ions)

    call profiling_in(prof, 'ELEC_MOVE_IONS')


    if (ion_dynamics_ions_move(ions_dyn)) then
      dt_ions = dt * ions_dyn%ionic_scale
      if (optional_default(save_pos, .false.)) then
        call ion_dynamics_save_state(ions_dyn, ions, wo%ions_state)
      end if
      call ion_dynamics_propagate(ions_dyn, ions, time, dt_ions, namespace)
      call hamiltonian_elec_epot_generate(hm, namespace, space, gr, ions, ext_partners, st, time = time)
    end if

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_move_ions)
  end subroutine propagation_ops_elec_move_ions

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_restore_ions(wo, ions_dyn, ions)
    class(propagation_ops_elec_t),    intent(inout) :: wo
    type(ion_dynamics_t),    intent(inout) :: ions_dyn
    type(ions_t),            intent(inout) :: ions

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_restore_ions)

    call profiling_in(prof, 'ELEC_RESTORE_IONS')

    if (ion_dynamics_ions_move(ions_dyn)) then
      call ion_dynamics_restore_state(ions_dyn, ions, wo%ions_state)
    end if

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_restore_ions)
  end subroutine propagation_ops_elec_restore_ions

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_propagate_gauge_field(wo, gfield, dt, time, save_gf)
    class(propagation_ops_elec_t), intent(inout) :: wo
    type(gauge_field_t),           intent(inout) :: gfield
    FLOAT,                         intent(in)    :: dt
    FLOAT,                         intent(in)    :: time
    logical,  optional,            intent(in)    :: save_gf

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_propagate_gauge_field)

    call profiling_in(prof, 'ELEC_MOVE_GAUGE')

    if (gauge_field_is_propagated(gfield)) then
      if (optional_default(save_gf, .false.)) then
        call gauge_field_get_vec_pot(gfield, wo%vecpot)
        call gauge_field_get_vec_pot_vel(gfield, wo%vecpot_vel)
      end if
      call gauge_field_do_algorithmic_operation(gfield, OP_VERLET_COMPUTE_ACC, dt, time)
    end if

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_propagate_gauge_field)
  end subroutine propagation_ops_elec_propagate_gauge_field

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_restore_gauge_field(wo, namespace, space, hm, mesh, ext_partners)
    class(propagation_ops_elec_t), intent(in)    :: wo
    type(namespace_t),             intent(in)    :: namespace
    type(space_t),                 intent(in)    :: space
    type(hamiltonian_elec_t),      intent(inout) :: hm
    class(mesh_t),                 intent(in)    :: mesh
    type(partner_list_t),          intent(in)    :: ext_partners

    type(profile_t), save :: prof
    type(gauge_field_t), pointer :: gfield

    PUSH_SUB(propagation_ops_elec_restore_gauge_field)

    call profiling_in(prof, 'ELEC_RESTORE_GAUGE')

    gfield => list_get_gauge_field(ext_partners)
    if (associated(gfield)) then
      call gauge_field_set_vec_pot(gfield, wo%vecpot)
      call gauge_field_set_vec_pot_vel(gfield, wo%vecpot_vel)
      call hm%update(mesh, namespace, space, ext_partners)
    end if

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_restore_gauge_field)
  end subroutine propagation_ops_elec_restore_gauge_field

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_exp_apply(te, namespace, st, mesh, hm, dt)
    type(exponential_t),      intent(inout) :: te
    type(namespace_t),        intent(in)    :: namespace
    type(states_elec_t),      intent(inout) :: st
    class(mesh_t),            intent(in)    :: mesh
    type(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT,                    intent(in)    :: dt

    integer :: ik, ib
    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_exp_apply)

    call profiling_in(prof, 'ELEC_EXP_APPLY')

    do ik = st%d%kpt%start, st%d%kpt%end
      call propagation_ops_do_pack(st, hm, st%group%block_start, ik)
      do ib = st%group%block_start, st%group%block_end
        if (ib + 1 <= st%group%block_end) call propagation_ops_do_pack(st, hm, ib+1, ik)
        call accel_set_stream(ib)

        call hamiltonian_elec_base_set_phase_corr(hm%hm_base, mesh, st%group%psib(ib, ik))
        if (hamiltonian_elec_inh_term(hm)) then
          call te%apply_batch(namespace, mesh, hm, st%group%psib(ib, ik), dt, &
            inh_psib = hm%inh_st%group%psib(ib, ik))
        else
          call te%apply_batch(namespace, mesh, hm, st%group%psib(ib, ik), dt)
        end if
        call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, mesh, st%group%psib(ib, ik))

        call propagation_ops_do_unpack(st, hm, ib, ik)
        if (ib-1 >= st%group%block_start) call propagation_ops_finish_unpack(st, hm, ib-1, ik)
      end do
      call propagation_ops_finish_unpack(st, hm, st%group%block_end, ik)
    end do

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_exp_apply)

  end subroutine propagation_ops_elec_exp_apply

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_fuse_density_exp_apply(te, namespace, st, gr, hm, dt, dt2, vmagnus)
    type(exponential_t),      intent(inout) :: te
    type(namespace_t),        intent(in)    :: namespace
    type(states_elec_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT,                    intent(in)    :: dt
    FLOAT, optional,          intent(in)    :: dt2
    FLOAT, optional,          intent(in)    :: vmagnus(:,:,:)

    integer :: ik, ib
    type(wfs_elec_t) :: zpsib_dt
    type(density_calc_t) :: dens_calc
    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_fuse_density_exp_apply)

    call profiling_in(prof, 'ELEC_FUSE_DENS_EXP_APPLY')

    call density_calc_init(dens_calc, st, gr, st%rho)

    do ik = st%d%kpt%start, st%d%kpt%end
      call propagation_ops_do_pack(st, hm, st%group%block_start, ik)
      do ib = st%group%block_start, st%group%block_end
        if (ib + 1 <= st%group%block_end) call propagation_ops_do_pack(st, hm, ib+1, ik)
        call accel_set_stream(ib)

        call hamiltonian_elec_base_set_phase_corr(hm%hm_base, gr, st%group%psib(ib, ik))
        if (present(dt2)) then
          call st%group%psib(ib, ik)%copy_to(zpsib_dt)
          if (st%group%psib(ib, ik)%is_packed()) call zpsib_dt%do_pack(copy = .false.)

          !propagate the state to dt/2 and dt, simultaneously, with H(time - dt)
          if (hamiltonian_elec_inh_term(hm)) then
            call te%apply_batch(namespace, gr, hm, st%group%psib(ib, ik), dt, psib2 = zpsib_dt, &
              deltat2 = dt2, inh_psib = hm%inh_st%group%psib(ib, ik))
          else
            call te%apply_batch(namespace, gr, hm, st%group%psib(ib, ik), dt, psib2 = zpsib_dt, &
              deltat2 = dt2)
          end if
          call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, gr, st%group%psib(ib, ik))
          call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, gr, zpsib_dt)

          !use the dt propagation to calculate the density
          call density_calc_accumulate(dens_calc, zpsib_dt)

          call zpsib_dt%end()
        else
          !propagate the state to dt with H(time - dt)
          if (hamiltonian_elec_inh_term(hm)) then
            call te%apply_batch(namespace, gr, hm, st%group%psib(ib, ik), dt, vmagnus=vmagnus, &
              inh_psib = hm%inh_st%group%psib(ib, ik))
          else
            call te%apply_batch(namespace, gr, hm, st%group%psib(ib, ik), dt, vmagnus=vmagnus)
          end if
          call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, gr, st%group%psib(ib, ik))

          !use the dt propagation to calculate the density
          call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))
        end if

        call propagation_ops_do_unpack(st, hm, ib, ik)
        if (ib-1 >= st%group%block_start) call propagation_ops_finish_unpack(st, hm, ib-1, ik)
      end do
      call propagation_ops_finish_unpack(st, hm, st%group%block_end, ik)
    end do

    call density_calc_end(dens_calc)

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_fuse_density_exp_apply)
  end subroutine propagation_ops_elec_fuse_density_exp_apply

  ! ---------------------------------------------------------
  subroutine propagation_ops_do_pack(st, hm, ib, ik)
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,                  intent(in)    :: ib
    integer,                  intent(in)    :: ik

    PUSH_SUB(propagation_ops_do_pack)
    if (hm%apply_packed()) then
      call accel_set_stream(ib)
      call st%group%psib(ib, ik)%do_pack(async=.true.)
      if (hamiltonian_elec_inh_term(hm)) call hm%inh_st%group%psib(ib, ik)%do_pack(async=.true.)
    end if
    POP_SUB(propagation_ops_do_pack)
  end subroutine propagation_ops_do_pack

  ! ---------------------------------------------------------
  subroutine propagation_ops_do_unpack(st, hm, ib, ik)
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,                  intent(in)    :: ib
    integer,                  intent(in)    :: ik

    PUSH_SUB(propagation_ops_do_unpack)
    if (hm%apply_packed()) then
      call accel_set_stream(ib)
      call st%group%psib(ib, ik)%do_unpack(async=.true.)
      if (hamiltonian_elec_inh_term(hm)) call hm%inh_st%group%psib(ib, ik)%do_unpack(async=.true.)
    end if
    POP_SUB(propagation_ops_do_unpack)
  end subroutine propagation_ops_do_unpack

  ! ---------------------------------------------------------
  subroutine propagation_ops_finish_unpack(st, hm, ib, ik)
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,                  intent(in)    :: ib
    integer,                  intent(in)    :: ik

    PUSH_SUB(propagation_ops_finish_unpack)
    if (hm%apply_packed()) then
      call accel_set_stream(ib)
      call st%group%psib(ib, ik)%finish_unpack()
      if (hamiltonian_elec_inh_term(hm)) call hm%inh_st%group%psib(ib, ik)%finish_unpack()
    end if
    POP_SUB(propagation_ops_finish_unpack)
  end subroutine propagation_ops_finish_unpack

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_interpolate_get(mesh, hm, interp)
    class(mesh_t),                   intent(in)    :: mesh
    type(hamiltonian_elec_t),        intent(inout) :: hm
    type(potential_interpolation_t), intent(inout) :: interp

    PUSH_SUB(propagation_ops_elec_interpolate_get)

    call potential_interpolation_get(interp, mesh%np, hm%d%nspin, 0, hm%vhxc, vtau = hm%vtau)

    POP_SUB(propagation_ops_elec_interpolate_get)

  end subroutine propagation_ops_elec_interpolate_get

  ! ---------------------------------------------------------
  subroutine calculate_mxll_dipole_field(hm, mesh, st)
    type(hamiltonian_elec_t),        intent(inout) :: hm
    class(mesh_t),                   intent(in)    :: mesh
    type(states_elec_t),             intent(in)    :: st

    FLOAT, allocatable :: density(:,:), total_density(:), mask_density(:)
    FLOAT :: integral_mask
    FLOAT, parameter :: density_threshold = CNST(1.0e-8)
    integer :: idir

    if (hm%mxll_coupling_mode == LENGTH_GAUGE_DIPOLE .or. &
      hm%mxll_coupling_mode == VELOCITY_GAUGE_DIPOLE) then

      if (hm%mxll_dipole_field == DIPOLE_AVERAGE) then
        SAFE_ALLOCATE(density(1:mesh%np,1:st%d%nspin))
        SAFE_ALLOCATE(total_density(1:mesh%np))
        SAFE_ALLOCATE(mask_density(size(total_density)))
        call states_elec_total_density(st, mesh, density)
        total_density = sum(density, dim=2)
        mask_density = merge(M_ONE, M_ZERO, total_density > density_threshold)
        integral_mask = dmf_integrate(mesh, mask_density)
        do idir = 1, mesh%box%dim
          ! field_mxll_dip will be E field or A field, depending on the mxll_coupling_mode
          hm%field_mxll_dip(idir) = dmf_integrate(mesh, mask_density*hm%field_mxll(:,idir))/integral_mask
        end do
        SAFE_DEALLOCATE_A(total_density)
        SAFE_DEALLOCATE_A(density)
        SAFE_DEALLOCATE_A(mask_density)

      elseif (hm%mxll_dipole_field == DIPOLE_AT_COM) then
        hm%field_mxll_dip(:) = hm%field_mxll(hm%center_of_mass_ip,:)
      end if
    end if

  end subroutine calculate_mxll_dipole_field

end module propagation_ops_elec_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
