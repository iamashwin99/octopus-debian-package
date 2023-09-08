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

module output_berkeleygw_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use debug_oct_m
  use exchange_operator_oct_m
  use fft_oct_m
  use fourier_shell_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use ions_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use symm_op_oct_m
  use symmetries_oct_m
  use v_ks_oct_m
#if defined(HAVE_BERKELEYGW)
  use wfn_rho_vxc_io_m
#endif
  use xc_oct_m

  implicit none

  private
  public ::                 &
    output_bgw_t,           &
    output_berkeleygw_init, &
    output_berkeleygw

  type output_bgw_t
    private
    integer           :: nbands
    integer           :: vxc_diag_nmin
    integer           :: vxc_diag_nmax
    integer           :: vxc_offdiag_nmin
    integer           :: vxc_offdiag_nmax
    logical           :: complex
    character(len=80) :: wfn_filename
    logical           :: calc_exchange
    logical           :: calc_vmtxel
    integer           :: vmtxel_ncband
    integer           :: vmtxel_nvband
    FLOAT             :: vmtxel_polarization(3)
  end type output_bgw_t

contains

  ! ---------------------------------------------------------
  subroutine output_berkeleygw_init(nst, namespace, bgw, periodic_dim)
    integer,            intent(in)  :: nst
    type(namespace_t),  intent(in)  :: namespace
    type(output_bgw_t), intent(out) :: bgw
    integer,            intent(in)  :: periodic_dim

    integer :: idir
    FLOAT :: norm
    type(block_t) :: blk

    PUSH_SUB(output_berkeleygw_init)

    call messages_experimental("BerkeleyGW output", namespace=namespace)

#ifndef HAVE_BERKELEYGW
    message(1) = "Cannot do BerkeleyGW output: the library was not linked."
    call messages_fatal(1, namespace=namespace)
#endif

    !%Variable BerkeleyGW_NumberBands
    !%Type integer
    !%Default all states
    !%Section Output::BerkeleyGW
    !%Description
    !% Wavefunctions for bands up to this number will be output. Must be between <= number of states.
    !% If < 1, no wavefunction file will be output.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_NumberBands', nst, bgw%nbands)

    ! these cannot be checked earlier, since output is initialized before unocc determines nst
    if (bgw%nbands > nst) then
      message(1) = "BerkeleyGW_NumberBands must be <= number of states."
      call messages_fatal(1, only_root_writes = .true., namespace=namespace)
    end if

    !%Variable BerkeleyGW_Vxc_diag_nmin
    !%Type integer
    !%Default 1
    !%Section Output::BerkeleyGW
    !%Description
    !% Lowest band for which to write diagonal exchange-correlation matrix elements. Must be <= number of states.
    !% If < 1, diagonals will be skipped.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_Vxc_diag_nmin', 1, bgw%vxc_diag_nmin)

    if (bgw%vxc_diag_nmin > nst) then
      message(1) = "BerkeleyGW_Vxc_diag_nmin must be <= number of states."
      call messages_fatal(1, only_root_writes = .true., namespace=namespace)
    end if

    !%Variable BerkeleyGW_Vxc_diag_nmax
    !%Type integer
    !%Default nst
    !%Section Output::BerkeleyGW
    !%Description
    !% Highest band for which to write diagonal exchange-correlation matrix elements. Must be between <= number of states.
    !% If < 1, diagonals will be skipped.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_Vxc_diag_nmax', nst, bgw%vxc_diag_nmax)

    if (bgw%vxc_diag_nmax > nst) then
      message(1) = "BerkeleyGW_Vxc_diag_nmax must be <= number of states."
      call messages_fatal(1, only_root_writes = .true., namespace=namespace)
    end if

    if (bgw%vxc_diag_nmin <= 0 .or. bgw%vxc_diag_nmax <= 0) then
      bgw%vxc_diag_nmin = 0
      bgw%vxc_diag_nmax = 0
    end if

    !%Variable BerkeleyGW_Vxc_offdiag_nmin
    !%Type integer
    !%Default 1
    !%Section Output::BerkeleyGW
    !%Description
    !% Lowest band for which to write off-diagonal exchange-correlation matrix elements. Must be <= number of states.
    !% If < 1, off-diagonals will be skipped.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_Vxc_offdiag_nmin', 1, bgw%vxc_offdiag_nmin)

    if (bgw%vxc_offdiag_nmin > nst) then
      message(1) = "BerkeleyGW_Vxc_offdiag_nmin must be <= number of states."
      call messages_fatal(1, only_root_writes = .true., namespace=namespace)
    end if

    !%Variable BerkeleyGW_Vxc_offdiag_nmax
    !%Type integer
    !%Default nst
    !%Section Output::BerkeleyGW
    !%Description
    !% Highest band for which to write off-diagonal exchange-correlation matrix elements. Must be <= number of states.
    !% If < 1, off-diagonals will be skipped.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_Vxc_offdiag_nmax', nst, bgw%vxc_offdiag_nmax)

    if (bgw%vxc_offdiag_nmax > nst) then
      message(1) = "BerkeleyGW_Vxc_offdiag_nmax must be <= number of states."
      call messages_fatal(1, only_root_writes = .true., namespace=namespace)
    end if

    if (bgw%vxc_offdiag_nmin <= 0 .or. bgw%vxc_offdiag_nmax <= 0) then
      bgw%vxc_offdiag_nmin = 0
      bgw%vxc_offdiag_nmax = 0
    end if

    !!%Variable BerkeleyGW_Complex
    !!%Type logical
    !!%Default false
    !!%Section Output::BerkeleyGW
    !!%Description
    !!% Even when wavefunctions, density, and XC potential could be real in reciprocal space,
    !!% they will be output as complex.
    !!%End
    !call parse_variable(namespace, 'BerkeleyGW_Complex', .false., bgw%complex)

    bgw%complex = .true.
    ! real output not implemented, so currently this is always true

    !%Variable BerkeleyGW_WFN_filename
    !%Type string
    !%Default WFN
    !%Section Output::BerkeleyGW
    !%Description
    !% Filename for the wavefunctions.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_WFN_filename', 'WFN', bgw%wfn_filename)

    !%Variable BerkeleyGW_CalcExchange
    !%Type logical
    !%Default false
    !%Section Output::BerkeleyGW
    !%Description
    !% Whether to calculate exchange matrix elements, to be written in <tt>x.dat</tt>.
    !% These will be calculated anyway by BerkeleyGW <tt>Sigma</tt>, so this is useful
    !% mainly for comparison and testing.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_CalcExchange', .false., bgw%calc_exchange)

    !%Variable BerkeleyGW_CalcDipoleMtxels
    !%Type logical
    !%Default false
    !%Section Output::BerkeleyGW
    !%Description
    !% Whether to calculate dipole matrix elements, to be written in <tt>vmtxel</tt>.
    !% This should be done when calculating <tt>WFN_fi</tt> for Bethe-Salpeter calculations
    !% with light polarization in a finite direction. In that case, a shifted grid
    !% <tt>WFNq_fi</tt> cannot be calculated, but we can instead use matrix elements of
    !% <math>r</math> in a more exact scheme. In <tt>absorption.inp</tt>, set <tt>read_vmtxel</tt>
    !% and <tt>use_momentum</tt>. Specify the number of conduction and valence bands you will
    !% use in BSE here with <tt>BerkeleyGW_VmtxelNumCondBands</tt> and <tt>BerkeleyGW_VmtxelNumValBands</tt>.
    !%End
    call parse_variable(namespace, 'BerkeleyGW_CalcDipoleMtxels', .false., bgw%calc_vmtxel)

    !%Variable BerkeleyGW_VmtxelPolarization
    !%Type block
    !%Default (1, 0, 0)
    !%Section Output::BerkeleyGW
    !%Description
    !% Polarization, <i>i.e.</i> direction vector, for which to calculate <tt>vmtxel</tt>, if you have set
    !% <tt>BerkeleyGW_CalcDipoleMtxels = yes</tt>. May not have any component in a periodic direction.
    !% The vector will be normalized.
    !%End

    bgw%vmtxel_polarization(1:3) = M_ZERO
    bgw%vmtxel_polarization(1) = M_ONE

    if (bgw%calc_vmtxel .and. parse_block(namespace, 'BerkeleyGW_VmtxelPolarization', blk) == 0) then
      do idir = 1, 3
        call parse_block_float(blk, 0, idir - 1, bgw%vmtxel_polarization(idir))

        if (idir <= periodic_dim .and. abs(bgw%vmtxel_polarization(idir)) > M_EPSILON) then
          message(1) = "You cannot calculate vmtxel with polarization in a periodic direction. Use WFNq_fi instead."
          call messages_fatal(1, only_root_writes = .true., namespace=namespace)
        end if
      end do
      call parse_block_end(blk)
      norm = sum(abs(bgw%vmtxel_polarization(1:3))**2)
      if (norm < M_EPSILON) then
        message(1) = "A non-zero value must be set for BerkeleyGW_VmtxelPolarization when BerkeleyGW_CalcDipoleMtxels = yes."
        call messages_fatal(1, namespace=namespace)
      end if
      bgw%vmtxel_polarization(1:3) = bgw%vmtxel_polarization(1:3) / sqrt(norm)
    end if

    !%Variable BerkeleyGW_VmtxelNumCondBands
    !%Type integer
    !%Default 0
    !%Section Output::BerkeleyGW
    !%Description
    !% Number of conduction bands for which to calculate <tt>vmtxel</tt>, if you have set
    !% <tt>BerkeleyGW_CalcDipoleMtxels = yes</tt>. This should be equal to the number to be
    !% used in BSE.
    !%End
    if (bgw%calc_vmtxel) call parse_variable(namespace, 'BerkeleyGW_VmtxelNumCondBands', 0, bgw%vmtxel_ncband)
    ! The default should be the minimum number of occupied states on any k-point or spin.

    !%Variable BerkeleyGW_VmtxelNumValBands
    !%Type integer
    !%Default 0
    !%Section Output::BerkeleyGW
    !%Description
    !% Number of valence bands for which to calculate <tt>vmtxel</tt>, if you have set
    !% <tt>BerkeleyGW_CalcDipoleMtxels = yes</tt>. This should be equal to the number to be
    !% used in BSE.
    !%End
    if (bgw%calc_vmtxel) call parse_variable(namespace, 'BerkeleyGW_VmtxelNumValBands', 0, bgw%vmtxel_nvband)
    ! The default should be the minimum number of unoccupied states on any k-point or spin.

    POP_SUB(output_berkeleygw_init)
  end subroutine output_berkeleygw_init


  ! ---------------------------------------------------------
  subroutine output_berkeleygw(bgw, namespace, space, dir, st, gr, ks, hm, ions)
    type(output_bgw_t),       intent(in)    :: bgw
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(states_elec_t),      intent(in)    :: st
    type(grid_t), target,     intent(in)    :: gr
    type(v_ks_t),             intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(ions_t),             intent(in)    :: ions

#ifdef HAVE_BERKELEYGW
    integer :: ik, is, ikk, ist, itran, iunit, iatom, mtrx(3, 3, 48), FFTgrid(3), ngkmax
    integer, pointer :: ifmin(:,:), ifmax(:,:), atyp(:), ngk(:)
    character(len=3) :: sheader
    FLOAT :: adot(3,3), bdot(3,3), recvol, tnp(3, 48), ecutrho, ecutwfc
    FLOAT, pointer :: energies(:,:,:), occupations(:,:,:), apos(:,:)
    FLOAT, allocatable :: vxc(:,:), dpsi(:,:)
    CMPLX, allocatable :: field_g(:,:), zpsi(:,:)
    type(cube_t) :: cube
    type(cube_function_t) :: cf
    type(fourier_shell_t) :: shell_density, shell_wfn
#endif

    PUSH_SUB(output_berkeleygw)

    if (space%dim /= 3) then
      message(1) = "BerkeleyGW output only available in 3D."
      call messages_fatal(1, namespace=namespace)
    end if

    if (st%d%ispin == SPINORS) call messages_not_implemented("BerkeleyGW output for spinors", namespace=namespace)

    if (st%parallel_in_states) call messages_not_implemented("BerkeleyGW output parallel in states", namespace=namespace)

    if (st%d%kpt%parallel) call messages_not_implemented("BerkeleyGW output parallel in k-points", namespace=namespace)

    if (ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK .or. xc_is_orbital_dependent(ks%xc)) then
      call messages_not_implemented("BerkeleyGW output with orbital-dependent functionals", namespace=namespace)
    end if

    if (hm%ep%nlcc) call messages_not_implemented("BerkeleyGW output with NLCC", namespace=namespace)

#ifdef HAVE_BERKELEYGW

    SAFE_ALLOCATE(vxc(1:gr%np, 1:st%d%nspin))
    vxc(:,:) = M_ZERO
    ! we should not include core rho here. that is why we do not just use hm%vxc
    call xc_get_vxc(gr, ks%xc, st, hm%kpoints, hm%psolver, namespace, space, st%rho, st%d%ispin, &
      hm%ions%latt%rcell_volume, vxc)

    message(1) = "BerkeleyGW output: vxc.dat"
    if (bgw%calc_exchange) message(1) = trim(message(1)) // ", x.dat"
    call messages_info(1, namespace=namespace)

    if (states_are_real(st)) then
      call dbgw_vxc_dat(bgw, namespace, space, dir, st, gr, hm, vxc)
    else
      call zbgw_vxc_dat(bgw, namespace, space, dir, st, gr, hm, vxc)
    end if

    call cube_init(cube, gr%idx%ll, namespace, space, gr%spacing, gr%coord_system, &
      fft_type = FFT_COMPLEX, dont_optimize = .true., nn_out = FFTgrid)
    call cube_init_cube_map(cube, gr)
    if (any(gr%idx%ll(1:3) /= FFTgrid(1:3))) then ! paranoia check
      message(1) = "Cannot do BerkeleyGW output: FFT grid has been modified."
      call messages_fatal(1, namespace=namespace)
    end if
    call zcube_function_alloc_rs(cube, cf)
    call cube_function_alloc_fs(cube, cf)

    ! NOTE: in BerkeleyGW, no G-vector may have coordinate equal to the half the FFT grid size.
    call fourier_shell_init(shell_density, namespace, space, cube, gr)
    ecutrho = shell_density%ekin_cutoff
    SAFE_ALLOCATE(field_g(1:shell_density%ngvectors, 1:st%d%nspin))

    call bgw_setup_header()


    if (bgw%calc_vmtxel) then
      write(message(1),'(a,3f12.6)') "BerkeleyGW output: vmtxel. Polarization = ", bgw%vmtxel_polarization(1:3)
      call messages_info(1, namespace=namespace)

      if (states_are_real(st)) then
        call dbgw_vmtxel(bgw, namespace, dir, st, gr, ifmax)
      else
        call zbgw_vmtxel(bgw, namespace, dir, st, gr, ifmax)
      end if
    end if

    message(1) = "BerkeleyGW output: VXC"
    call messages_info(1, namespace=namespace)

    sheader = 'VXC'
    if (mpi_grp_is_root(mpi_world)) then
      iunit = io_open(trim(dir) // 'VXC', namespace, form = 'unformatted', action = 'write')
      call bgw_write_header(sheader, iunit)
    end if
    ! convert from Ha to Ry, make usable with same processing as RHO
    vxc(:,:) = vxc(:,:) * M_TWO / (product(cube%rs_n_global(1:3)) * gr%volume_element)
    call dbgw_write_FS(namespace, iunit, vxc, field_g, shell_density, st%d%nspin, gr, cube, cf, is_wfn = .false.)
    if (mpi_grp_is_root(mpi_world)) call io_close(iunit)
    SAFE_DEALLOCATE_A(vxc)


    message(1) = "BerkeleyGW output: RHO"
    call messages_info(1, namespace=namespace)

    sheader = 'RHO'
    if (mpi_grp_is_root(mpi_world)) then
      iunit = io_open(trim(dir) // 'RHO', namespace, form = 'unformatted', action = 'write')
      call bgw_write_header(sheader, iunit)
    end if
    call dbgw_write_FS(namespace, iunit, st%rho, field_g, shell_density, st%d%nspin, gr, cube, cf, is_wfn = .false.)
    if (mpi_grp_is_root(mpi_world)) call io_close(iunit)

    message(1) = "BerkeleyGW output: WFN"
    write(message(2),'(a,f12.6,a)') "Wavefunction cutoff for BerkeleyGW: ", &
      fourier_shell_cutoff(space, cube, gr, .true.) * M_TWO, " Ry"
    call messages_info(2, namespace=namespace)

    if (states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(1:gr%np, 1:st%d%nspin))
    else
      SAFE_ALLOCATE(zpsi(1:gr%np, 1:st%d%nspin))
    end if

    sheader = 'WFN'
    if (mpi_grp_is_root(mpi_world)) then
      iunit = io_open(trim(dir) // bgw%wfn_filename, namespace, form = 'unformatted', action = 'write')
      call bgw_write_header(sheader, iunit)
    end if

    call fourier_shell_end(shell_density)

    ! FIXME: is parallelization over k-points possible?
    do ik = st%d%kpt%start, st%d%kpt%end, st%d%nspin
      call fourier_shell_init(shell_wfn, namespace, space, cube, gr, kk = hm%kpoints%reduced%red_point(:, ik))

      if (mpi_grp_is_root(mpi_world)) then
        call write_binary_gvectors(iunit, shell_wfn%ngvectors, shell_wfn%ngvectors, shell_wfn%red_gvec)
      end if
      do ist = 1, st%nst
        do is = 1, st%d%nspin
          ikk = ik + is - 1
          if (states_are_real(st)) then
            call states_elec_get_state(st, gr, 1, ist, ikk, dpsi(:, is))
          else
            call states_elec_get_state(st, gr, 1, ist, ikk, zpsi(:, is))
          end if
        end do
        if (states_are_real(st)) then
          call dbgw_write_FS(namespace, iunit, dpsi, field_g, shell_wfn, st%d%nspin, gr, cube, cf, is_wfn = .true.)
        else
          call zbgw_write_FS(namespace, iunit, zpsi, field_g, shell_wfn, st%d%nspin, gr, cube, cf, is_wfn = .true.)
        end if
      end do
      call fourier_shell_end(shell_wfn)
    end do

    if (mpi_grp_is_root(mpi_world)) call io_close(iunit)

    ! deallocate everything
    call cube_function_free_fs(cube, cf)
    call zcube_function_free_rs(cube, cf)
    call cube_end(cube)

    if (states_are_real(st)) then
      SAFE_DEALLOCATE_A(dpsi)
    else
      SAFE_DEALLOCATE_A(zpsi)
    end if
    SAFE_DEALLOCATE_A(vxc)
    SAFE_DEALLOCATE_A(field_g)
    SAFE_DEALLOCATE_P(ifmin)
    SAFE_DEALLOCATE_P(ifmax)
    SAFE_DEALLOCATE_P(ngk)
    SAFE_DEALLOCATE_P(energies)
    SAFE_DEALLOCATE_P(occupations)
    SAFE_DEALLOCATE_P(atyp)
    SAFE_DEALLOCATE_P(apos)

#else
    message(1) = "Cannot do BerkeleyGW output: the library was not linked."
    call messages_fatal(1, namespace=namespace)
#endif

    POP_SUB(output_berkeleygw)

#ifdef HAVE_BERKELEYGW
  contains

    subroutine bgw_setup_header()
      PUSH_SUB(output_berkeleygw.bgw_setup_header)

      if (space%periodic_dim /= 3) then
        message(1) = "BerkeleyGW for mixed-periodicity is currently not implemented."
        call messages_fatal(1, namespace=namespace)
      end if

      ! The rlattice, klattice and rcell_volume used here are not correct for
      ! mixid periodicity. Note also that the BerkeleyGW treats the z direction
      ! as periodic for wires, while in Octopus it is the x direction that is
      ! periodic.
      adot(1:3, 1:3) = matmul(ions%latt%rlattice(1:3, 1:3), ions%latt%rlattice(1:3, 1:3))
      bdot(1:3, 1:3) = matmul(ions%latt%klattice(1:3, 1:3), ions%latt%klattice(1:3, 1:3))
      recvol = (M_TWO * M_PI)**3 / ions%latt%rcell_volume

      ! symmetry is not analyzed by Octopus for finite systems, but we only need it for periodic ones
      do itran = 1, symmetries_number(gr%symm)
        mtrx(:,:, itran) = symm_op_rotation_matrix_red(gr%symm%ops(itran))
        tnp(:, itran) = symm_op_translation_vector_red(gr%symm%ops(itran))
      end do
      ! some further work on conventions of mtrx and tnp is required!

      SAFE_ALLOCATE(ifmin(1:hm%kpoints%reduced%npoints, 1:st%d%nspin))
      SAFE_ALLOCATE(ifmax(1:hm%kpoints%reduced%npoints, 1:st%d%nspin))
      SAFE_ALLOCATE(energies(1:st%nst, 1:hm%kpoints%reduced%npoints, 1:st%d%nspin))
      SAFE_ALLOCATE(occupations(1:st%nst, 1:hm%kpoints%reduced%npoints, 1:st%d%nspin))

      ifmin(:,:) = 1
!     This is how semiconducting smearing "should" work, but not in our implementation.
!      if (smear_is_semiconducting(st%smear)) then
!        ifmax(:,:) = nint(st%qtot / st%smear%el_per_state)
!      end if
      do ik = 1, st%d%nik
        is = st%d%get_spin_index(ik)
        ikk = st%d%get_kpoint_index(ik)
        energies(1:st%nst, ikk, is) = st%eigenval(1:st%nst,ik) * M_TWO
        occupations(1:st%nst, ikk, is) = st%occ(1:st%nst, ik) / st%smear%el_per_state
        do ist = 1, st%nst
          ! M_EPSILON needed since e_fermi is top of valence band for fixed_occ and semiconducting smearing
          if (st%eigenval(ist, ik) < st%smear%e_fermi + M_EPSILON) then
            ifmax(ikk, is) = ist
          else
            exit
          end if
        end do
      end do

      SAFE_ALLOCATE(ngk(1:hm%kpoints%reduced%npoints))
      do ik = 1, st%d%nik, st%d%nspin
        call fourier_shell_init(shell_wfn, namespace, space, cube, gr, kk = hm%kpoints%reduced%red_point(:, ik))
        if (ik == 1) ecutwfc = shell_wfn%ekin_cutoff ! should be the same for all, anyway
        ngk(ik) = shell_wfn%ngvectors
        call fourier_shell_end(shell_wfn)
      end do
      ngkmax = maxval(ngk)

      SAFE_ALLOCATE(atyp(1:ions%natoms))
      SAFE_ALLOCATE(apos(1:3, 1:ions%natoms))
      do iatom = 1, ions%natoms
        atyp(iatom) = species_index(ions%atom(iatom)%species)
        apos(1:3, iatom) = ions%pos(1:3, iatom)
      end do

      if (any(hm%kpoints%nik_axis(1:3) == 0)) then
        message(1) = "KPointsGrid has a zero component. Set KPointsGrid appropriately,"
        message(2) = "or this WFN will only be usable in BerkeleyGW's inteqp."
        call messages_warning(1, namespace=namespace)
      end if

      POP_SUB(output_berkeleygw.bgw_setup_header)
    end subroutine bgw_setup_header

    ! ---------------------------------------------------------
    subroutine bgw_write_header(sheader, iunit)
      character(len=3), intent(inout) :: sheader
      integer,          intent(in)    :: iunit

      FLOAT, pointer :: weight(:), red_point(:,:)

      PUSH_SUB(output_berkeleygw.bgw_write_header)

      weight => hm%kpoints%reduced%weight
      red_point => hm%kpoints%reduced%red_point

      call write_binary_header(iunit, sheader, 2, st%d%nspin, shell_density%ngvectors, &
        symmetries_number(gr%symm), 0, ions%natoms, &
        hm%kpoints%reduced%npoints, st%nst, ngkmax, ecutrho * M_TWO,  &
        ecutwfc * M_TWO, FFTgrid, hm%kpoints%nik_axis, hm%kpoints%full%shifts, &
        ions%latt%rcell_volume, M_ONE, ions%latt%rlattice, adot, recvol, &
        M_ONE, ions%latt%klattice, bdot, mtrx, tnp, atyp, &
        apos, ngk, weight, red_point, &
        ifmin, ifmax, energies, occupations, warn = .false.)

      call write_binary_gvectors(iunit, shell_density%ngvectors, shell_density%ngvectors, shell_density%red_gvec)

      POP_SUB(output_berkeleygw.bgw_write_header)
    end subroutine bgw_write_header

#endif

  end subroutine output_berkeleygw

#include "undef.F90"
#include "complex.F90"
#ifdef HAVE_BERKELEYGW
#include "output_berkeleygw_inc.F90"
#endif

#include "undef.F90"
#include "real.F90"
#ifdef HAVE_BERKELEYGW
#include "output_berkeleygw_inc.F90"
#endif

end module output_berkeleygw_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
