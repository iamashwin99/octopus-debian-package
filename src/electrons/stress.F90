!! Copyright (C) 2002-2016 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2023 N. Tancogne-Dejean
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

! ---------------------------------------------------------
!> @brief This module implements the calculation of the stress tensor
!!
module stress_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use debug_oct_m
  use density_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use energy_calc_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use interaction_partner_oct_m
  use hamiltonian_elec_base_oct_m
  use ions_oct_m
  use ion_interaction_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lda_u_oct_m
  use loct_math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use ps_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use splines_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use wfs_elec_oct_m
  use xc_f03_lib_m

  implicit none

  private
  public ::                    &
    stress_calculate,          &
    output_stress,             &
    output_pressure

contains

  ! ---------------------------------------------------------
  !> @brief This computes the total stress on the lattice
  subroutine stress_calculate(namespace, gr, hm, st, ions, ks, ext_partners)
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(inout) :: gr !< grid
    type(hamiltonian_elec_t), intent(inout) :: hm !< the Hamiltonian
    type(states_elec_t),      intent(inout) :: st !< the electronic states
    type(ions_t),             intent(inout) :: ions !< geometry
    type(v_ks_t),             intent(in)    :: ks   !< the Kohn-Sham system
    type(partner_list_t),     intent(in)    :: ext_partners !< external interaction partners

    FLOAT, allocatable    :: rho_total(:)
    type(profile_t), save :: stress_prof
    FLOAT :: stress(3,3) ! stress tensor in Cartesian coordinate
    FLOAT :: stress_kin(3,3), stress_Hartree(3,3), stress_xc(3,3), stress_xc_nlcc(3,3)
    FLOAT :: stress_ps(3,3), stress_ps_nl(3,3), stress_ps_local(3,3), stress_ii(3,3)
    integer :: ip
    FLOAT, allocatable :: vh(:) !< Hartree potential
    FLOAT, allocatable :: grad_vh(:,:) !< Gradient of the Hartree potential
    FLOAT :: ehartree !< Hartree energy


    call profiling_in(stress_prof, "STRESS_CALCULATE")
    PUSH_SUB(stress_calculate)

    if (st%wfs_type /= TYPE_CMPLX) then
      write(message(1),'(a)') 'The stress tensors for real wavefunctions has not been implemented!'

      if (hm%kpoints%full%npoints == 1) then
        write(message(2),'(a)') 'For testing this feature, you can add ForceComplex=yes to the input file'
        call messages_fatal(2, namespace=namespace)
      end if

      call messages_fatal(1, namespace=namespace)
    end if

    if (ions%space%periodic_dim /= 3) then
      call messages_not_implemented("Stress tensors for periodicity different from 3D", namespace=namespace)
    end if

    if ( .not. (ks%theory_level == KOHN_SHAM_DFT .and. bitand(hm%xc%family, XC_FAMILY_LDA) /= 0) &
      .and. .not. (ks%theory_level == INDEPENDENT_PARTICLES)) then
      write(message(1),'(a)') 'The stress tensor computation is currently only possible at the Kohn-Sham DFT level'
      write(message(2),'(a)') 'with LDA functionals or for independent particles.'
      call messages_fatal(2, namespace=namespace)
    end if

    if (ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
      write(message(1),'(a)') 'The stress tensor is currently not properly computed with vdW corrections'
      call messages_fatal(1, namespace=namespace)
    end if

    if (hm%pcm%run_pcm) then
      call messages_not_implemented('Stress tensor with PCM')
    end if

    if (hm%lda_u_level /= DFT_U_NONE) then
      call messages_not_implemented('Stress tensor with DFT+U')
    end if

    if (allocated(hm%v_static)) then
      call messages_not_implemented('Stress tensor with static electric fields')
    end if

    if (ks%has_photons) then
      call messages_not_implemented('Stress tensor with photon modes')
    end if

    if (.not. hm%hm_base%apply_projector_matrices) then
      call messages_not_implemented('Stress tensor with relativistic Kleinman-Bylander pseudopotential')
    end if

    stress(:,:) = M_ZERO

    SAFE_ALLOCATE(rho_total(1:gr%np_part))
    do ip = 1, gr%np
      rho_total(ip) = sum(st%rho(ip, 1:st%d%nspin))
    end do

    ! As we rely on some of the full energy components, we need to recompute it first
    ! TODO: We should restrict the components of the energy needed to be computed
    call energy_calc_total(namespace, ions%space, hm, gr, st, ext_partners, iunit = 0, full = .true.)

    ! In order to get the electrostatic part (Hartree and local pseudopotential part),
    ! we need to get the Hartree potential and its gradient
    SAFE_ALLOCATE(vh(1:gr%np_part))
    SAFE_ALLOCATE(grad_vh(1:gr%np, 1:gr%der%dim))
    call lalg_copy(gr%np, hm%vhartree, vh)
    ehartree = hm%energy%hartree

    ! We also compute the gradient here
    call dderivatives_grad(gr%der, vh, grad_vh)

    ! We now compute the various contributions to the stress tensor

    ! Stress from kinetic energy of electrons
    call stress_from_kinetic(gr%der, ions%space, hm, st, gr%symm, ions%latt%rcell_volume, stress_kin)
    stress = stress + stress_kin

    if (ks%theory_level == INDEPENDENT_PARTICLES) then
      stress_Hartree = M_ZERO
      stress_xc = M_ZERO
    else
      call stress_from_Hartree(gr, ions%space, ions%latt%rcell_volume, vh, grad_vh, ehartree, stress_Hartree)
      stress = stress + stress_Hartree

      call stress_from_xc(hm%energy, ions%latt%rcell_volume, stress_xc)

      ! Nonlinear core correction contribution
      if (allocated(st%rho_core)) then
        call stress_from_xc_nlcc(ions%latt%rcell_volume, gr, st, ions, hm%vxc, stress_xc_nlcc)
        stress_xc = stress_xc + stress_xc_nlcc
      end if
      stress = stress + stress_xc

    end if

    call stress_from_pseudo_local(gr, st, hm, hm%ep%proj, ions, rho_total, vh, grad_vh, stress_ps_local)
    stress_ps = stress_ps_local
    stress = stress + stress_ps_local

    SAFE_DEALLOCATE_A(vh)
    SAFE_DEALLOCATE_A(grad_vh)

    call stress_from_pseudo_nonloc(gr, st, hm, hm%ep%proj, ions, stress_ps_nl)
    stress_ps = stress_ps + stress_ps_nl
    stress = stress + stress_ps_nl

    call ion_interaction_stress(ions%ion_interaction, ions%space, ions%latt, ions%atom, ions%natoms, ions%pos, stress_ii)
    stress = stress + stress_ii

    ! Stress from kinetic energy of ion
    ! Stress from ion-field interaction

    ! Sign changed to fit conventional definition
    stress = - stress

    st%stress_tensors%total(1:3,1:3) = stress(1:3,1:3)
    st%stress_tensors%kinetic(1:3,1:3) = stress_kin(1:3,1:3)
    st%stress_tensors%Hartree(1:3,1:3) = stress_Hartree(1:3,1:3)
    st%stress_tensors%xc(1:3,1:3) = stress_xc(1:3,1:3)
    st%stress_tensors%pseudopotential(1:3,1:3) = stress_ps(1:3,1:3)
    st%stress_tensors%ion_ion(1:3,1:3) = stress_ii(1:3,1:3)

    SAFE_DEALLOCATE_A(rho_total)

    POP_SUB(stress_calculate)
    call profiling_out(stress_prof)
  end subroutine stress_calculate

  ! -------------------------------------------------------
  !> @brief Computes the contribution to the stress tensor from the kinetic energy
  !!
  !! We use the real space formula from Sharma and Suryanarayana
  !! On the calculation of the stress tensor in real-space Kohn-Sham density functional theory
  !! J. Chem. Phys. 149, 194104 (2018)
  !!
  !! More precisely, this routines computes
  !! \f[
  !! \sigma_{ij} = \frac{1}{V}\sum_n\sum_k^{BZ} w_kf_{n,k}\int d^3r \partial_i \psi^*_{n,k}(r) \partial_j \psi_{n,k}(r)\,
  !! \f]
  !!
  !! where \f$V\f$ is the cell volume, \f$ w_k\f$ is the weight of the k-point k,
  !! \f$ f_{n,k}\f$ is the occupation number of the band  with a k-point index k,
  !! and \f$ \psi_{n,k}\f$ is the corresponding Bloch state.
  subroutine stress_from_kinetic(der, space, hm, st, symm, rcell_volume, stress_kin)
    type(derivatives_t),            intent(in)    :: der
    type(space_t),                  intent(in)    :: space
    type(hamiltonian_elec_t),       intent(in)    :: hm
    type(states_elec_t),            intent(inout) :: st
    type(symmetries_t),             intent(in)    :: symm
    FLOAT,                          intent(in)    :: rcell_volume
    FLOAT,                          intent(out)   :: stress_kin(3, 3)

    integer :: ik, ist, idir, jdir, ib, minst, maxst
    CMPLX, allocatable :: stress_l_block(:)
    type(profile_t), save :: prof
    type(wfs_elec_t) :: psib, gpsib(space%dim)

    call profiling_in(prof, "STRESS_FROM_KINETIC")
    PUSH_SUB(stress_from_kinetic)

    stress_kin(:,:) = M_ZERO

    SAFE_ALLOCATE(stress_l_block(1:st%d%block_size))

    do ik = st%d%kpt%start, st%d%kpt%end
      if (st%d%kweights(ik) <= M_EPSILON) cycle

      do ib = st%group%block_start, st%group%block_end
        minst = states_elec_block_min(st, ib)
        maxst = states_elec_block_max(st, ib)

        call st%group%psib(ib, ik)%do_pack(copy = .true.)
        call st%group%psib(ib, ik)%copy_to(psib)
        ! set the boundary conditions
        call boundaries_set(der%boundaries, der%mesh, st%group%psib(ib, ik))

        ! set the phase for periodic systems
        if (allocated(hm%hm_base%phase)) then
          call hamiltonian_elec_base_phase(hm%hm_base, der%mesh, der%mesh%np_part, &
            conjugate = .false., psib = psib, src = st%group%psib(ib, ik))
        else
          call st%group%psib(ib, ik)%copy_data_to(der%mesh%np_part, psib)
        end if

        ! calculate the gradient
        do idir = 1, space%dim
          call psib%copy_to(gpsib(idir))
        end do
        call zderivatives_batch_grad(der, psib, gpsib, set_bc=.false.)

        ! Accumulate the result
        do idir = 1, der%dim
          do jdir = idir, der%dim
            call zmesh_batch_dotp_vector(der%mesh, gpsib(idir), gpsib(jdir), stress_l_block)

            do ist = minst, maxst
              stress_kin(idir,jdir) = stress_kin(idir,jdir) &
                + st%d%kweights(ik) * st%occ(ist, ik) &
                * TOFLOAT(stress_l_block(ist - minst + 1))
            end do
          end do
        end do

        do idir = 1, space%dim
          call gpsib(idir)%end()
        end do
        call psib%end()
        call st%group%psib(ib, ik)%do_unpack(copy = .false.)

      end do
    end do

    if (st%parallel_in_states .or. st%d%kpt%parallel) then
      call comm_allreduce(st%st_kpt_mpi_grp, stress_kin)
    end if


    ! Symmetrize the kinetic stress tensor
    do idir = 1, der%dim
      do jdir = idir+1, der%dim
        stress_kin(jdir,idir) = stress_kin(idir,jdir)
      end do
    end do

    ! Symmetrize the stress tensor if we use k-point symmetries
    if (hm%kpoints%use_symmetries) then
      call dsymmetrize_tensor_cart(symm, stress_kin)
    end if

    stress_kin = stress_kin / rcell_volume

    call profiling_out(prof)
    POP_SUB(stress_from_kinetic)
  end subroutine stress_from_kinetic

  ! -------------------------------------------------------
  !> @brief Computes the contribution to the stress tensor from the Hartree energy
  !!
  !! We use the real space formula from Sharma and Suryanarayana
  !! On the calculation of the stress tensor in real-space Kohn-Sham density functional theory
  !! J. Chem. Phys. 149, 194104 (2018)
  !!
  !! More precisely, this routines computes
  !! \f[
  !! \sigma_{ij} = \frac{1}{4\pi V}\int d^3r \partial_i v_{\rm H}(r) \partial_j v_{\rm H}(r) - \delta_{ij} \frac{1}{V}E_{\rm H}\,
  !! \f]
  !!
  !! where \f$V\f$ is the cell volume, \f$ v_{\rm H}(r) \f$ is the Hartree potential, and \f$ E_{\rm H} \f$ is the Hartree energy
  !! defined as
  !! \f[
  !! E_{\rm H} = \frac{1}{2} \int d^3r n(r) v_{\rm H}(r) \,.
  !! \f]
  subroutine stress_from_Hartree(gr, space, volume, vh, grad_vh, ehartree, stress_Hartree)
    type(grid_t),       intent(in)    :: gr
    type(space_t),      intent(in)    :: space
    FLOAT,              intent(in)    :: volume
    FLOAT,              intent(in)    :: vh(:) !< Hartree potential
    FLOAT,              intent(in)    :: grad_vh(:,:) !< Gradient of the Hartree potential
    FLOAT,              intent(in)    :: ehartree !< Hartree      U = (1/2)*Int [n v_Hartree]
    FLOAT,              intent(out)   :: stress_Hartree(3, 3)

    integer :: idir, jdir
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_HARTREE")
    PUSH_SUB(stress_from_Hartree)

    stress_Hartree(:,:) = M_ZERO

    do idir = 1,3
      do jdir = idir,3
        stress_Hartree(idir, jdir) = -dmf_dotp(gr, grad_vh(:,idir), grad_vh(:, jdir))/M_FOUR/M_PI
        stress_Hartree(jdir, idir) = stress_Hartree(idir, jdir)
      end do
      stress_Hartree(idir, idir) = stress_Hartree(idir, idir) + ehartree
    end do

    stress_Hartree =  stress_Hartree/volume

    call profiling_out(prof)
    POP_SUB(stress_from_Hartree)
  end subroutine stress_from_Hartree


  ! -------------------------------------------------------
  !> @brief Computes the contribution to the stress tensor from the xc energy
  !!
  !! This routines computes the xc stress tensor assuming an LDA functional
  !! \f[
  !! \sigma_{ij} = \frac{1}{V}\delta_{ij}(-\int d^3r n(r) v_{\rm xc}(r) + E_{\rm xc})\,
  !! \f]
  !!
  !! where \f$V\f$ is the cell volume, \f$n(r) \f$ is the electronic density,
  !! \f$ v_{\rm xc}(r) \f$ is the exchange-correlation potential,
  !! and \f$ E_{\rm xc} \f$ is the exchange-correlation energy.
  !!
  ! Note: We assume hm%energy%echange, correlation, and intnvxc
  ! have already been calculated before.
  subroutine stress_from_xc(energy, rcell_volume, stress_xc)
    type(energy_t),           intent(in)    :: energy
    FLOAT,                    intent(in)    :: rcell_volume
    FLOAT,                    intent(out)   :: stress_xc(3, 3)

    integer :: idir
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_XC")
    PUSH_SUB(stress_from_xc)

    stress_xc = M_ZERO
    do idir = 1, 3
      stress_xc(idir, idir) = - energy%exchange - energy%correlation + energy%intnvxc
    end do
    stress_xc(:,:) = stress_xc(:,:) / rcell_volume

    call profiling_out(prof)
    POP_SUB(stress_from_xc)
  end subroutine stress_from_xc


  ! -------------------------------------------------------
  !> @brief Computes the NLCC contribution to the stress tensor from the xc energy
  !!
  !! The nonlinear core correction term is given by
  !! \f[
  !! \sigma_{ij}^{\rm NLCC} = \frac{1}{V}\int d^3r v_{xc}(r) \frac{\partial \rho_{\rm NLCC}(\epsilon r)}{\partial \epsilon_{ij}}\Bigg|_{\epsilon=I}\,.
  !! \f]
  subroutine stress_from_xc_nlcc(rcell_volume, gr, st, ions, vxc, stress_xc_nlcc)
    FLOAT,                    intent(in)    :: rcell_volume
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(in)    :: st
    type(ions_t),             intent(in)    :: ions
    FLOAT,                    intent(in)    :: vxc(:,:)
    FLOAT,                    intent(out)   :: stress_xc_nlcc(3, 3)

    integer :: idir, jdir, iat
    FLOAT, allocatable :: gnlcc(:,:), gnlcc_x(:,:,:), vxc_tot(:)
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_XC_NLCC")
    PUSH_SUB(stress_from_xc_nlcc)

    ASSERT(allocated(st%rho_core))

    stress_xc_nlcc = M_ZERO

    ! We first accumulate the contribution from all the pseudo-ions
    SAFE_ALLOCATE(gnlcc(gr%np, gr%der%dim))
    SAFE_ALLOCATE(gnlcc_x(gr%np, gr%der%dim, gr%der%dim))
    gnlcc_x = M_ZERO
    do iat = ions%atoms_dist%start, ions%atoms_dist%end
      ASSERT(species_is_ps(ions%atom(iat)%species))
      call species_get_nlcc_grad(ions%atom(iat)%species, ions%space, ions%latt, &
        ions%pos(:,iat), gr, gnlcc, gnlcc_x)
    end do
    SAFE_DEALLOCATE_A(gnlcc)

    if (ions%atoms_dist%parallel) then
      call comm_allreduce(ions%atoms_dist%mpi_grp, gnlcc_x)
    end if

    ! Sum over spin of the xc potential
    SAFE_ALLOCATE(vxc_tot(1:gr%np))
    vxc_tot = vxc(1:gr%np, 1)
    if(st%d%nspin > 1) vxc_tot = vxc_tot + vxc(1:gr%np, 2)

    do idir = 1, 3
      do jdir = idir, 3
        stress_xc_nlcc(idir, jdir) = dmf_dotp(gr, vxc_tot, gnlcc_x(:,idir, jdir))
        stress_xc_nlcc(jdir, idir) = stress_xc_nlcc(idir, jdir)
      end do
    end do
    SAFE_DEALLOCATE_A(vxc_tot)
    SAFE_DEALLOCATE_A(gnlcc_x)

    stress_xc_nlcc(:,:) = stress_xc_nlcc(:,:) / rcell_volume

    call profiling_out(prof)
    POP_SUB(stress_from_xc_nlcc)
  end subroutine stress_from_xc_nlcc

  ! -------------------------------------------------------
  !> @brief Computes the contribution to the stress tensor from the nonlocal part of the pseudopotentials
  !!
  !! More precisely, this routines computes
  !! \f[
  !! \sigma_{ij} = \frac{2}{V}\sum_n\sum_k^{BZ} \sum_I w_kf_{n,k} \langle \partial_i \psi_{n,k}| (r_j-R_j) \hat{V}_{\rm NL,I}| \psi_{n,k}\rangle
  !! + \delta_{ij} \frac{1}{V}E_{\rm NL}
  !! \f]
  !!
  !! where \f$V\f$ is the cell volume, \f$ w_k\f$ is the weight of the k-point k,
  !! \f$ f_{n,k}\f$ is the occupation number of the band  with a k-point index k,
  !! \f$ \psi_{n,k}\f$ is the corresponding Bloch state, and \f$ \hat{V}_{\rm NL, I}\f$ is the
  !! pseudopotential non-local operator from atom I, centered on the position \f$R_I\f$.
  !! \f$ E_{\rm NL} \f$ is the nonlocal energy from the nonlocal pseudopotential.
  !!
  !! See Sharma and Suryanarayana,
  !! On the calculation of the stress tensor in real-space Kohn-Sham density functional theory,
  !! J. Chem. Phys. 149, 194104 (2018) for more details
  !!
  subroutine stress_from_pseudo_nonloc(gr, st, hm, proj, ions, stress_ps_nl)
    type(grid_t),      target,           intent(in) :: gr !< grid
    type(states_elec_t),              intent(inout) :: st
    type(hamiltonian_elec_t),            intent(in) :: hm
    type(projector_t), target,           intent(in) :: proj(:)
    type(ions_t),                        intent(in) :: ions !< ions
    FLOAT,                              intent(out) :: stress_ps_nl(3, 3)

    integer :: ik, ist, idir, jdir
    integer :: ib, minst, maxst
    type(profile_t), save :: prof
    type(wfs_elec_t) :: psib, rvnl_psib(3), gpsib(3)
    CMPLX, allocatable :: stress_tmp(:)

    call profiling_in(prof, "STRESS_FROM_PSEUDO_NL")
    PUSH_SUB(stress_from_pseudo_nonloc)

    ASSERT(st%wfs_type == TYPE_CMPLX)

    SAFE_ALLOCATE(stress_tmp(1:st%d%block_size))

    stress_ps_nl = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end

      if (st%d%kweights(ik) <= M_EPSILON) cycle

      do ib = st%group%block_start, st%group%block_end
        minst = states_elec_block_min(st, ib)
        maxst = states_elec_block_max(st, ib)

        call st%group%psib(ib, ik)%do_pack(copy = .true.)
        call st%group%psib(ib, ik)%copy_to(psib)
        ! set the boundary conditions
        call boundaries_set(gr%der%boundaries, gr, st%group%psib(ib, ik))

        ! set the phase for periodic systems
        if (allocated(hm%hm_base%phase)) then
          call hamiltonian_elec_base_phase(hm%hm_base, gr, gr%np_part, &
            conjugate = .false., psib = psib, src = st%group%psib(ib, ik))
        else
          call st%group%psib(ib, ik)%copy_data_to(gr%np_part, psib)
        end if

        ! calculate the gradient
        do idir = 1, gr%der%dim
          call psib%copy_to(gpsib(idir))
        end do
        call zderivatives_batch_grad(gr%der, psib, gpsib, set_bc=.false.)


        ! Get rV_NL |\psi> for all atoms
        do idir = 1, gr%der%dim
          call psib%copy_to(rvnl_psib(idir))
          call batch_set_zero(rvnl_psib(idir))
        end do
        call zhamiltonian_elec_base_r_vnlocal(hm%hm_base, gr, st%d, &
          gr%der%boundaries%spiral, psib, rvnl_psib)

        do idir = 1,3
          do jdir = idir,3
            call zmesh_batch_dotp_vector(gr, gpsib(idir), rvnl_psib(jdir), stress_tmp)

            do ist = minst, maxst
              stress_ps_nl(idir, jdir) = stress_ps_nl(idir, jdir) &
                + M_TWO * st%d%kweights(ik) * st%occ(ist, ik) * TOFLOAT(stress_tmp(ist-minst+1))
            end do

          end do
        end do

        do idir = 1, gr%der%dim
          call rvnl_psib(idir)%end()
          call gpsib(idir)%end()
        end do
        call psib%end()
        call st%group%psib(ib, ik)%do_unpack(copy = .false.)
      end do
    end do

    SAFE_DEALLOCATE_A(stress_tmp)

    if (st%parallel_in_states .or. st%d%kpt%parallel) then
      call comm_allreduce(st%st_kpt_mpi_grp, stress_ps_nl)
    end if

    ! Symmetrize the kinetic stress tensor
    do idir = 1, gr%der%dim
      do jdir = idir+1, gr%der%dim
        stress_ps_nl(jdir,idir) = stress_ps_nl(idir,jdir)
      end do
    end do

    ! Symmetrize the stress tensor if we use k-point symmetries
    if (hm%kpoints%use_symmetries) then
      call dsymmetrize_tensor_cart(gr%symm, stress_ps_nl)
    end if

    ! Add the nonlocal energy
    do idir = 1,3
      stress_ps_nl(idir, idir) = stress_ps_nl(idir, idir) + hm%energy%extern_non_local
    end do

    stress_ps_nl = stress_ps_nl/ions%latt%rcell_volume

    call profiling_out(prof)
    POP_SUB(stress_from_pseudo_nonloc)

  end subroutine stress_from_pseudo_nonloc


  ! -------------------------------------------------------
  !> @brief Computes the contribution from the local part of the pseudopotential
  !!
  !! We use a real space formulation, which computes two parts. One short range (SR)
  !! \f[
  !! \sigma_{ij} =\frac{1}{V} \sum_I\int d^3r [\partial_i \rho(r)] (x_j-R_{I,j}) v_{\rm SR, I}(r) + \delta_{ij} \frac{1}{V}E_{\rm loc, SR}\,
  !! \f]
  !!
  !! and a long range part
  !!
  !! \f[
  !! \sigma_{ij} = \frac{1}{V}\int d\vec{x} \sum_I n_{\rm LR}^{I}(|\vec{x}-\vec{R}_I|) (x_j-R_{I,j}) [\partial_i v_H(\vec{x})]
  !! + \delta_{ij} \frac{1}{V} E_{\rm loc, LR} -\frac{2}{4\pi V}\int d\vec{r} [\partial_i v_{\rm H}(\vec{r})][\partial_j v_{\rm LR}(\vec{r})]
  !! \f]
  !!
  !! where \f$V\f$ is the cell volume, \f$ v_{\rm H}(r) \f$ is the Hartree potential,
  !! \f$\vec{R}_I\f$ refers to the position of the atom \f$I\f$, and \f$n_{\rm LR}^{I}\f$ is the long-range density associated with the long-range potential \f$v_{\rm LR}^I\f$.
  subroutine stress_from_pseudo_local(gr, st, hm, proj, ions, rho_total, vh, grad_vh, stress_ps_local)
    type(grid_t),      target,           intent(in) :: gr !< grid
    type(states_elec_t),              intent(inout) :: st
    type(hamiltonian_elec_t),            intent(in) :: hm
    type(projector_t), target,           intent(in) :: proj(:)
    type(ions_t),                        intent(in) :: ions !< ions
    FLOAT, contiguous,                intent(inout) :: rho_total(:)
    FLOAT,                               intent(in) :: vh(:) !< Hartree potential
    FLOAT,                               intent(in) :: grad_vh(:,:) !< Gradient of the Hartree potential
    FLOAT,                              intent(out) :: stress_ps_local(3, 3)


    FLOAT :: stress_SR(3, 3), stress_LR(3, 3)
    FLOAT :: energy_ps_SR
    FLOAT,  allocatable :: vloc(:), rvloc(:,:), rho_local_lr(:), rho_lr(:)
    FLOAT,  allocatable :: grad_rho(:,:), rho_lr_x(:,:), vlr(:), grad_vlr(:,:)
    integer :: idir, jdir, iatom
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_PSEUDO_LOC")
    PUSH_SUB(stress_from_pseudo_local)

    ! calculate stress from short-range local pseudopotentials
    stress_SR = M_ZERO

    SAFE_ALLOCATE(vloc(1:gr%np))
    vloc = M_ZERO
    SAFE_ALLOCATE(rvloc(1:gr%np, 1:gr%der%dim))
    rvloc = M_ZERO
    do iatom = 1, ions%natoms
      call epot_local_pseudopotential_SR(gr, ions, iatom, vloc, rvloc)
    end do
    SAFE_DEALLOCATE_A(vloc)

    SAFE_ALLOCATE(grad_rho(1:gr%np,1:gr%der%dim))
    call dderivatives_grad(gr%der, rho_total, grad_rho)

    energy_ps_SR = hm%energy%extern_local
    do idir = 1,3
      do jdir = idir,3
        stress_SR(idir, jdir) = stress_SR(idir, jdir) &
          +dmf_dotp(gr, rvloc(:, jdir), grad_rho(:, idir))
        stress_SR(jdir, idir) = stress_SR(idir, jdir)
      end do
      stress_SR(idir,idir) = stress_SR(idir,idir) + energy_ps_SR
    end do

    stress_SR = stress_SR/ions%latt%rcell_volume

    SAFE_DEALLOCATE_A(rvloc)
    SAFE_DEALLOCATE_A(grad_rho)


    ! calculate stress from long-range local pseudopotentials
    stress_LR = M_ZERO

    ! We treat the long-range part of the local potential as the Hartree term
    ! We first sum the long range densities from atoms
    SAFE_ALLOCATE(rho_lr(1:gr%np_part))
    SAFE_ALLOCATE(rho_lr_x(1:gr%np, 1:gr%der%dim))
    rho_lr = M_ZERO
    rho_lr_x = M_ZERO
    SAFE_ALLOCATE(rho_local_lr(1:gr%np))
    do iatom = ions%atoms_dist%start, ions%atoms_dist%end
      ASSERT(species_is_ps(ions%atom(iatom)%species))
      call species_get_long_range_density(ions%atom(iatom)%species, ions%namespace, ions%space, ions%latt, &
        ions%pos(:, iatom), gr, rho_local_lr, nlr_x=rho_lr_x)

      call lalg_axpy(gr%np, M_ONE, rho_local_lr, rho_lr)
    end do
    SAFE_DEALLOCATE_A(rho_local_lr)

    if (ions%atoms_dist%parallel) then
      call comm_allreduce(ions%atoms_dist%mpi_grp, rho_lr)
      call comm_allreduce(ions%atoms_dist%mpi_grp, rho_lr_x)
    end if

    do idir = 1, 3
      do jdir = idir, 3
        stress_LR(idir, jdir) = stress_LR(idir, jdir) + dmf_dotp(gr, rho_lr_x(:,jdir), grad_vh(:, idir))
      end do
    end do
    SAFE_DEALLOCATE_A(rho_lr_x)

    SAFE_ALLOCATE(vlr(1:gr%np_part))
    if (poisson_solver_is_iterative(hm%psolver)) then
      ! vl has to be initialized before entering routine
      ! and our best guess for the potential is zero
      vlr(1:gr%np) = M_ZERO
    end if
    call dpoisson_solve(hm%psolver, ions%namespace, vlr, rho_lr, all_nodes = .true.)
    SAFE_DEALLOCATE_A(rho_lr)

    SAFE_ALLOCATE(grad_vlr(1:gr%np, 1:gr%der%dim))
    call dderivatives_grad(gr%der, vlr, grad_vlr)
    SAFE_DEALLOCATE_A(vlr)

    do idir = 1,3
      do jdir = idir, 3
        stress_LR(idir, jdir) = stress_LR(idir, jdir) - dmf_dotp(gr, grad_vh(:,idir), grad_vlr(:, jdir))/M_TWO/M_PI
        stress_LR(jdir, idir) = stress_LR(idir, jdir)
      end do
    end do

    SAFE_DEALLOCATE_A(grad_vlr)

    stress_LR = stress_LR/ions%latt%rcell_volume

    stress_ps_local = stress_SR + stress_LR

!!! NOTE!! This part is moved to Ewald contribution
!! Contribution from G=0 component of the long-range part
!    charge = M_ZERO
!    do iatom = 1, ions%natoms
!       zi = species_zval(ions%atom(iatom)%species)
!       charge = charge + zi
!    end do
!
!    do idir = 1,3
!       stress_ps(idir,idir) = stress_ps(idir,idir) &
!            + 2d0*M_PI*sigma_erf**2*charge**2 /ions%latt%rcell_volume**2
!    end do

    call profiling_out(prof)
    POP_SUB(stress_from_pseudo_local)

  end subroutine stress_from_pseudo_local

  ! -------------------------------------------------------
  subroutine epot_local_pseudopotential_SR(mesh, ions, iatom, vpsl, rvpsl)
    class(mesh_t),            intent(in)    :: mesh
    type(ions_t),             intent(in)    :: ions
    integer,                  intent(in)    :: iatom
    FLOAT,                    intent(inout) :: vpsl(:)
    FLOAT,                    intent(inout) :: rvpsl(:,:)

    integer :: ip
    FLOAT :: radius, vl_ip
    type(submesh_t)  :: sphere
    type(profile_t), save :: prof
    type(ps_t), pointer :: ps

    PUSH_SUB(epot_local_pseudopotential_sr)
    call profiling_in(prof, "EPOT_LOCAL_PS_SR")

    !the localized part

    if (species_is_ps(ions%atom(iatom)%species)) then

      ps => species_ps(ions%atom(iatom)%species)

      radius = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold)*CNST(1.05)

      call submesh_init(sphere, ions%space, mesh, ions%latt, ions%pos(:, iatom), radius)

      ! Cannot be written (correctly) as a vector expression since for periodic systems,
      ! there can be values ip, jp such that sphere%map(ip) == sphere%map(jp).
      do ip = 1, sphere%np
        vl_ip = spline_eval(ps%vl, sphere%r(ip))
        vpsl(sphere%map(ip)) = vpsl(sphere%map(ip)) + vl_ip
        rvpsl(sphere%map(ip) ,1) = rvpsl(sphere%map(ip), 1) + sphere%rel_x(1, ip) * vl_ip
        rvpsl(sphere%map(ip) ,2) = rvpsl(sphere%map(ip), 2) + sphere%rel_x(2, ip) * vl_ip
        rvpsl(sphere%map(ip) ,3) = rvpsl(sphere%map(ip), 3) + sphere%rel_x(3, ip) * vl_ip
      end do

      call submesh_end(sphere)

      nullify(ps)

    end if


    call profiling_out(prof)
    POP_SUB(epot_local_pseudopotential_sr)
  end subroutine epot_local_pseudopotential_SR


  ! -------------------------------------------------------
  subroutine output_stress(iunit, namespace, space_dim, stress_tensors, all_terms)
    integer,           intent(in) :: iunit
    type(namespace_t), intent(in) :: namespace
    integer,           intent(in) :: space_dim
    type(stress_t),    intent(in) :: stress_tensors
    logical, optional, intent(in) :: all_terms  ! write each contributing term separately

    logical :: write_all_terms
    character(len=16) :: stress_unit

    if (present(all_terms)) then
      write_all_terms = all_terms
    else
      write_all_terms = .true.
    end if

    write(stress_unit, '(4a,i1)') trim(units_abbrev(units_out%energy)), '/', &
      trim(units_abbrev(units_out%length)), '^', space_dim

    if (mpi_grp_is_root(mpi_world)) then

      if (write_all_terms) then
        write(iunit, '(3a)') 'Kinetic stress tensor [', trim(stress_unit), '] ='
        call print_stress_tensor(iunit, space_dim, stress_tensors%kinetic)

        write(iunit, '(3a)') 'Hartree stress tensor [', trim(stress_unit), '] ='
        call print_stress_tensor(iunit, space_dim, stress_tensors%Hartree)

        write(iunit, '(3a)') 'XC stress tensor [', trim(stress_unit), '] ='
        call print_stress_tensor(iunit, space_dim, stress_tensors%xc)

        write(iunit, '(3a)') 'Pseudopotential stress tensor [', trim(stress_unit), '] ='
        call print_stress_tensor(iunit, space_dim, stress_tensors%pseudopotential)

        write(iunit, '(3a)') 'Ion-ion stress tensor [', trim(stress_unit), '] ='
        call print_stress_tensor(iunit, space_dim, stress_tensors%ion_ion)
      end if

      write(iunit, '(3a)') 'Total stress tensor [', trim(stress_unit), '] ='
      call print_stress_tensor(iunit, space_dim, stress_tensors%total)

    end if
  end subroutine output_stress


  subroutine output_pressure(iunit, space_dim, total_stress_tensor)
    integer, intent(in) :: iunit
    integer, intent(in) :: space_dim
    FLOAT,   intent(in) :: total_stress_tensor(3,3)

    integer :: idim
    FLOAT :: pressure = M_ZERO
    character(len=16) :: stress_unit

    write(stress_unit, '(4a,i1)') trim(units_abbrev(units_out%energy)), '/', &
      trim(units_abbrev(units_out%length)), '^', space_dim

    do idim = 1, space_dim
      pressure = pressure - total_stress_tensor(idim, idim) / TOFLOAT(space_dim)
    end do

    write(iunit,'(3a,es16.8)', advance="no") 'Pressure [', trim(stress_unit), '] = ', &
      units_from_atomic(units_out%energy/units_out%length**space_dim, pressure)
    if (space_dim == 3) then
      write(iunit,'(2x,a,f16.8)') 'Pressure [GPa] = ', pressure * CNST(29421.02648438959)
    else
      write(iunit,*)
    end if

  end subroutine output_pressure

  subroutine print_stress_tensor(ounit, space_dim, tensor)
    integer, intent(in) :: ounit
    integer, intent(in) :: space_dim
    FLOAT,   intent(in) :: tensor(3,3)

    FLOAT   :: tensor_with_unit(3,3)
    integer :: idim, jdim

    tensor_with_unit = units_from_atomic(units_out%energy/units_out%length**space_dim, tensor)

    write(ounit,'(a9,2x)', advance="no")"T_{ij}"
    do jdim = 1, space_dim
      write(ounit,'(i18)', advance="no") jdim
    end do
    write(ounit,*)
    do idim = 1, space_dim
      write(ounit,'(i9,2x)', advance="no") idim
      do jdim = 1, space_dim
        write(ounit,'(es18.9)', advance="no") tensor_with_unit(idim, jdim)
      end do
      write(ounit,*)
    end do
    write(ounit,*)

  end subroutine print_stress_tensor


end module stress_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
