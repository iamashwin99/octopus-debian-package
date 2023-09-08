!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module energy_calc_oct_m
  use batch_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use exchange_operator_oct_m
  use ext_partner_list_oct_m
  use fourier_space_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lda_u_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use interaction_partner_oct_m
  use pcm_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use singularity_oct_m
  use smear_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use wfs_elec_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                               &
    energy_calc_total,                    &
    denergy_calc_electronic,              &
    zenergy_calc_electronic,              &
    done_body_matrix_elements,            &
    zone_body_matrix_elements,            &
    dtwo_body_matrix_elements,            &
    ztwo_body_matrix_elements,            &
    energy_calc_eigenvalues,              &
    energy_calc_virial_ex

contains

  ! ---------------------------------------------------------
  !> This subroutine calculates the total energy of the system. Basically, it
  !! adds up the KS eigenvalues, and then it subtracts whatever double
  !! counts exist (see TDDFT theory for details).
  subroutine energy_calc_total(namespace, space, hm, gr, st, ext_partners, iunit, full)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(inout) :: st
    type(partner_list_t),     intent(in)    :: ext_partners
    integer, optional,        intent(in)    :: iunit
    logical, optional,        intent(in)    :: full

    logical :: full_
    type(states_elec_t) :: xst !< The states after the application of the Fock operator
    !!                            This is needed to construct the ACE operator
    type(gauge_field_t), pointer :: gfield

    PUSH_SUB(energy_calc_total)

    full_ = .false.
    if (present(full)) full_ = full

    hm%energy%eigenvalues = states_elec_eigenvalues_sum(st)

    if (full_ .or. hm%theory_level == HARTREE .or. hm%theory_level == HARTREE_FOCK &
      .or. hm%theory_level == GENERALIZED_KOHN_SHAM_DFT) then
      if (states_are_real(st)) then
        hm%energy%kinetic  = denergy_calc_electronic(namespace, hm, gr%der, st, terms = TERM_KINETIC)
        hm%energy%extern_local = denergy_calc_electronic(namespace, hm, gr%der, st, terms = TERM_LOCAL_EXTERNAL)
        hm%energy%extern_non_local   = denergy_calc_electronic(namespace, hm, gr%der, st, terms = TERM_NON_LOCAL_POTENTIAL)
        hm%energy%extern = hm%energy%extern_local + hm%energy%extern_non_local
      else
        hm%energy%kinetic = zenergy_calc_electronic(namespace, hm, gr%der, st, terms = TERM_KINETIC)

        hm%energy%extern_local = zenergy_calc_electronic(namespace, hm, gr%der, st, terms = TERM_LOCAL_EXTERNAL)

        hm%energy%extern_non_local = zenergy_calc_electronic(namespace, hm, gr%der, st, terms = TERM_NON_LOCAL_POTENTIAL)
        hm%energy%extern   = hm%energy%extern_local   + hm%energy%extern_non_local
      end if
    end if

    ! Computes the HF energy when not doing the ACE
    if( full_ .and. (hm%theory_level == HARTREE_FOCK .or. &
      (hm%theory_level == GENERALIZED_KOHN_SHAM_DFT .and. family_is_hybrid(hm%xc))) &
      .and. .not. hm%exxop%useACE  ) then
      call xst%nullify()
      if (states_are_real(st)) then
        call dexchange_operator_compute_potentials(hm%exxop, namespace, space, gr, st, &
          xst, hm%kpoints, hm%energy%exchange_hf)
      else
        call zexchange_operator_compute_potentials(hm%exxop, namespace, space, gr, st, &
          xst, hm%kpoints, hm%energy%exchange_hf)
      end if
      call states_elec_end(xst)
    else
      if(.not. hm%exxop%useACE) hm%energy%exchange_hf = M_ZERO
    end if

    if (hm%pcm%run_pcm) then
      hm%pcm%counter = hm%pcm%counter + 1
      if (hm%pcm%localf) then
        call pcm_elect_energy(hm%ions, hm%pcm, hm%energy%int_ee_pcm, hm%energy%int_en_pcm, &
          hm%energy%int_ne_pcm, hm%energy%int_nn_pcm, &
          E_int_e_ext = hm%energy%int_e_ext_pcm,      &
          E_int_n_ext = hm%energy%int_n_ext_pcm       )
      else
        call pcm_elect_energy(hm%ions, hm%pcm, hm%energy%int_ee_pcm, hm%energy%int_en_pcm, &
          hm%energy%int_ne_pcm, hm%energy%int_nn_pcm  )
      end if
    end if

    select case (hm%theory_level)
    case (INDEPENDENT_PARTICLES)
      hm%energy%total = hm%ep%eii + hm%energy%eigenvalues

    case (HARTREE)
      hm%energy%total = hm%ep%eii + M_HALF * (hm%energy%eigenvalues + hm%energy%kinetic + hm%energy%extern)

    case (HARTREE_FOCK, GENERALIZED_KOHN_SHAM_DFT)
      hm%energy%total = hm%ep%eii + &
        M_HALF*(hm%energy%eigenvalues + hm%energy%kinetic + hm%energy%extern - hm%energy%intnvxc &
        - hm%energy%int_dft_u) &
        + hm%energy%exchange + hm%energy%exchange_hf + hm%energy%correlation &
        + hm%energy%vdw - hm%energy%intnvstatic + hm%energy%dft_u

      ! FIXME: pcm terms are only added to total energy in DFT case

    case (KOHN_SHAM_DFT)
      hm%energy%total = hm%ep%eii + hm%energy%eigenvalues &
        - hm%energy%hartree + hm%energy%exchange + hm%energy%correlation + hm%energy%vdw - hm%energy%intnvxc &
        - hm%energy%pcm_corr + hm%energy%int_ee_pcm + hm%energy%int_en_pcm &
        + hm%energy%int_nn_pcm + hm%energy%int_ne_pcm &
        + hm%energy%int_e_ext_pcm + hm%energy%int_n_ext_pcm &
        + hm%energy%dft_u -  hm%energy%int_dft_u - hm%energy%intnvstatic &
        + hm%energy%photon_exchange

    end select

    hm%energy%entropy = smear_calc_entropy(st%smear, st%eigenval, st%d%nik, st%nst, st%d%kweights, st%occ)
    if (st%smear%method == SMEAR_FIXED_OCC) then ! no temperature available
      hm%energy%TS = M_ZERO
    else
      hm%energy%TS = st%smear%dsmear * hm%energy%entropy
    end if

    gfield => list_get_gauge_field(ext_partners)
    if(associated(gfield)) hm%energy%total = hm%energy%total + gauge_field_get_energy(gfield)

    if (allocated(hm%vberry)) then
      hm%energy%total = hm%energy%total + hm%energy%berry
    else
      hm%energy%berry = M_ZERO
    end if

    if (optional_default(iunit, -1) > 0) then
      write(message(1), '(6x,a, f18.8)')'Total       = ', units_from_atomic(units_out%energy, hm%energy%total)
      write(message(2), '(6x,a, f18.8)')'Free        = ', units_from_atomic(units_out%energy, hm%energy%total - hm%energy%TS)
      write(message(3), '(6x,a)') '-----------'
      call messages_info(3, iunit)

      write(message(1), '(6x,a, f18.8)')'Ion-ion     = ', units_from_atomic(units_out%energy, hm%ep%eii)
      write(message(2), '(6x,a, f18.8)')'Eigenvalues = ', units_from_atomic(units_out%energy, hm%energy%eigenvalues)
      write(message(3), '(6x,a, f18.8)')'Hartree     = ', units_from_atomic(units_out%energy, hm%energy%hartree)
      write(message(4), '(6x,a, f18.8)')'Int[n*v_xc] = ', units_from_atomic(units_out%energy, hm%energy%intnvxc)
      write(message(5), '(6x,a, f18.8)')'Exchange    = ', &
        units_from_atomic(units_out%energy, hm%energy%exchange+hm%energy%exchange_hf)
      write(message(6), '(6x,a, f18.8)')'Correlation = ', units_from_atomic(units_out%energy, hm%energy%correlation)
      write(message(7), '(6x,a, f18.8)')'vanderWaals = ', units_from_atomic(units_out%energy, hm%energy%vdw)
      write(message(8), '(6x,a, f18.8)')'Delta XC    = ', units_from_atomic(units_out%energy, hm%energy%delta_xc)
      write(message(9), '(6x,a, f18.8)')'Entropy     = ', hm%energy%entropy ! the dimensionless sigma of Kittel&Kroemer
      write(message(10), '(6x,a, f18.8)')'-TS         = ', -units_from_atomic(units_out%energy, hm%energy%TS)
      write(message(11), '(6x,a, f18.8)')'Photon ex.  = ', units_from_atomic(units_out%energy, hm%energy%photon_exchange)
      call messages_info(11, iunit)

      if (hm%pcm%run_pcm) then
        write(message(1),'(6x,a, f18.8)')'E_e-solvent = ',  units_from_atomic(units_out%energy, hm%energy%int_ee_pcm + &
          hm%energy%int_en_pcm + &
          hm%energy%int_e_ext_pcm)
        write(message(2),'(6x,a, f18.8)')'E_n-solvent = ',  units_from_atomic(units_out%energy, hm%energy%int_nn_pcm + &
          hm%energy%int_ne_pcm + &
          hm%energy%int_n_ext_pcm)
        write(message(3),'(6x,a, f18.8)')'E_M-solvent = ',  units_from_atomic(units_out%energy, &
          hm%energy%int_ee_pcm + hm%energy%int_en_pcm + &
          hm%energy%int_nn_pcm + hm%energy%int_ne_pcm + &
          hm%energy%int_e_ext_pcm + hm%energy%int_n_ext_pcm)
        call messages_info(3, iunit)
      end if

      if (full_) then
        write(message(1), '(6x,a, f18.8)')'Kinetic     = ', units_from_atomic(units_out%energy, hm%energy%kinetic)
        write(message(2), '(6x,a, f18.8)')'External    = ', units_from_atomic(units_out%energy, hm%energy%extern)
        write(message(3), '(6x,a, f18.8)')'Non-local   = ', units_from_atomic(units_out%energy, hm%energy%extern_non_local)
        write(message(4), '(6x,a, f18.8)')'Int[n*v_E]  = ', units_from_atomic(units_out%energy, hm%energy%intnvstatic)
        call messages_info(4, iunit)
      end if
      if (allocated(hm%vberry) .and. space%is_periodic()) then
        write(message(1), '(6x,a, f18.8)')'Berry       = ', units_from_atomic(units_out%energy, hm%energy%berry)
        call messages_info(1, iunit)
      end if
      if (hm%lda_u_level /= DFT_U_NONE) then
        write(message(1), '(6x,a, f18.8)')'Hubbard     = ', units_from_atomic(units_out%energy, hm%energy%dft_u)
        write(message(2), '(6x,a, f18.8)')'Int[n*v_U]  = ', units_from_atomic(units_out%energy, hm%energy%int_dft_u)
        call messages_info(2, iunit)
      end if
    end if

    POP_SUB(energy_calc_total)
  end subroutine energy_calc_total

  ! --------------------------------------------------------------------

  subroutine energy_calc_eigenvalues(namespace, hm, der, st)
    type(namespace_t),        intent(in)    :: namespace
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(derivatives_t),      intent(in)    :: der
    type(states_elec_t),      intent(inout) :: st

    PUSH_SUB(energy_calc_eigenvalues)

    if (states_are_real(st)) then
      call dcalculate_eigenvalues(namespace, hm, der, st)
    else
      call zcalculate_eigenvalues(namespace, hm, der, st)
    end if

    POP_SUB(energy_calc_eigenvalues)
  end subroutine energy_calc_eigenvalues

  ! ---------------------------------------------------------
  subroutine energy_calc_virial_ex(der, vxc, st, ex)
    type(derivatives_t),      intent(in) :: der
    FLOAT,                    intent(in) :: vxc(:,:)
    type(states_elec_t),      intent(in) :: st
    FLOAT,                   intent(out) :: ex

    integer :: idir, ip, isp
    FLOAT, allocatable :: gradvx(:,:), nrgradvx(:)
    FLOAT :: rr, xx(3)

    PUSH_SUB(energy_calc_virial_ex)

    ASSERT(st%d%ispin /= SPINORS)

    SAFE_ALLOCATE(nrgradvx(1:der%mesh%np_part))
    SAFE_ALLOCATE(gradvx(1:der%mesh%np, 1:der%dim))

    ex = M_ZERO
    do isp = 1, st%d%nspin
      !We output the energy from the virial relation
      nrgradvx(1:der%mesh%np) = vxc(1:der%mesh%np,isp)
      call dderivatives_grad(der, nrgradvx, gradvx)
      nrgradvx = M_ZERO
      do idir = 1, der%dim
        do ip = 1, der%mesh%np
          call mesh_r(der%mesh, ip, rr, coords=xx)
          nrgradvx(ip) = nrgradvx(ip) - gradvx(ip, idir) * st%rho(ip, isp) * xx(idir)
        end do
      end do
      ex = ex + dmf_integrate(der%mesh, nrgradvx)
    end do

    SAFE_DEALLOCATE_A(gradvx)
    SAFE_DEALLOCATE_A(nrgradvx)

    POP_SUB(energy_calc_virial_ex)
  end subroutine energy_calc_virial_ex


#include "undef.F90"
#include "real.F90"
#include "energy_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "energy_calc_inc.F90"

end module energy_calc_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
