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

module v_ks_oct_m
  use accel_oct_m
  use comm_oct_m
  use current_oct_m
  use debug_oct_m
  use density_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use energy_calc_oct_m
  use exchange_operator_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use interaction_partner_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use lda_u_oct_m
  use libvdwxc_oct_m
  use magnetic_oct_m
  use magnetic_constrain_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use photon_mode_mf_oct_m
  use photon_mode_oct_m
  use profiling_oct_m
  use pseudo_oct_m
  use pcm_potential_oct_m
  use sort_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_parallel_oct_m
  use varinfo_oct_m
  use vdw_ts_oct_m
  use xc_oct_m
  use xc_f03_lib_m
  use xc_functional_oct_m
  use xc_interaction_oct_m
  use xc_ks_inversion_oct_m
  use xc_oep_oct_m
  use xc_sic_oct_m

  ! from the dftd3 library
  use dftd3_api

  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    v_ks_h_setup,       &
    v_ks_calc,          &
    v_ks_calc_t,        &
    v_ks_calc_start,    &
    v_ks_calc_finish,   &
    v_ks_freeze_hxc,    &
    v_ks_calculate_current

  type v_ks_calc_t
    private
    logical                           :: calculating
    logical                           :: time_present
    FLOAT                             :: time
    FLOAT,                allocatable :: density(:, :)
    logical                           :: total_density_alloc
    FLOAT,                pointer     :: total_density(:)
    type(energy_t),       allocatable :: energy
    type(states_elec_t),  pointer     :: hf_st
    FLOAT,                allocatable :: vxc(:, :)
    FLOAT,                allocatable :: vtau(:, :)
    FLOAT,                allocatable :: axc(:, :, :)
    FLOAT,                allocatable :: a_ind(:, :)
    FLOAT,                allocatable :: b_ind(:, :)
    logical                           :: calc_energy

    FLOAT,                allocatable :: vdw_forces(:, :)
  end type v_ks_calc_t

  type v_ks_t
    private
    integer,                  public :: theory_level = -1

    logical,                  public :: frozen_hxc = .false. !< For RPA and SAE calculations.

    integer,                  public :: xc_family = 0  !< the XC stuff
    integer,                  public :: xc_flags = 0  !< the XC flags
    type(xc_t),               public :: xc
    type(xc_oep_t),           public :: oep
    type(xc_ks_inversion_t),  public :: ks_inversion
    type(xc_sic_t),           public :: sic
    type(grid_t), pointer,    public :: gr
    type(v_ks_calc_t)                :: calc
    logical                          :: calculate_current = .false.
    type(current_t)                  :: current_calculator
    integer,                  public :: vdw_correction = -1
    logical                          :: vdw_self_consistent = .false.
    type(vdw_ts_t),           public :: vdw_ts
    type(dftd3_calc)                 :: vdw_d3
    logical                          :: include_td_field = .false.
    logical,                  public :: has_photons = .false.
    type(photon_mode_t),      public :: pt
    type(mf_t),               public :: pt_mx
  end type v_ks_t

contains

  ! ---------------------------------------------------------
  subroutine v_ks_init(ks, namespace, gr, st, ions, mc, space, kpoints)
    type(v_ks_t),            intent(inout) :: ks
    type(namespace_t),       intent(in)    :: namespace
    type(grid_t),    target, intent(inout) :: gr
    type(states_elec_t),     intent(in)    :: st
    type(ions_t),            intent(inout) :: ions
    type(multicomm_t),       intent(in)    :: mc
    type(space_t),           intent(in)    :: space
    type(kpoints_t),         intent(in)    :: kpoints

    integer :: x_id, c_id, xk_id, ck_id, default, val, iatom
    logical :: parsed_theory_level, using_hartree_fock
    type(dftd3_input) :: d3_input
    character(len=20) :: d3func_def, d3func
    integer :: pseudo_x_functional, pseudo_c_functional
    integer :: oep_type

    PUSH_SUB(v_ks_init)

    ! We need to parse TheoryLevel and XCFunctional, this is
    ! complicated because they are interdependent.

    !%Variable TheoryLevel
    !%Type integer
    !%Section Hamiltonian
    !%Description
    !% The calculations can be run with different "theory levels" that
    !% control how electrons are simulated. The default is
    !% <tt>dft</tt>. When hybrid functionals are requested, through
    !% the <tt>XCFunctional</tt> variable, the default is
    !% <tt>hartree_fock</tt>.
    !%Option independent_particles 2
    !% Particles will be considered as independent, <i>i.e.</i> as non-interacting.
    !% This mode is mainly used for testing purposes, as the code is usually
    !% much faster with <tt>independent_particles</tt>.
    !%Option hartree 1
    !% Calculation within the Hartree method (experimental). Note that, contrary to popular
    !% belief, the Hartree potential is self-interaction-free. Therefore, this run
    !% mode will not yield the same result as <tt>kohn-sham</tt> without exchange-correlation.
    !%Option hartree_fock 3
    !% This is the traditional Hartree-Fock scheme. Like the Hartree scheme, it is fully
    !% self-interaction-free.
    !%Option kohn_sham 4
    !% This is the default density-functional theory scheme. Note that you can also use
    !% hybrid functionals in this scheme, but they will be handled the "DFT" way, <i>i.e.</i>,
    !% solving the OEP equation.
    !%Option generalized_kohn_sham 5
    !% This is similar to the <tt>kohn-sham</tt> scheme, except that this allows for nonlocal operators.
    !% This is the default mode to run hybrid functionals, meta-GGA functionals, or DFT+U.
    !% It can be more convenient to use <tt>kohn-sham</tt> DFT within the OEP scheme to get similar (but not the same) results.
    !% Note that within this scheme you can use a correlation functional, or a hybrid
    !% functional (see <tt>XCFunctional</tt>). In the latter case, you will be following the
    !% quantum-chemistry recipe to use hybrids.
    !%Option rdmft 7
    !% (Experimental) Reduced Density Matrix functional theory.
    !%End

    ks%xc_family = XC_FAMILY_NONE
    ks%sic%level = SIC_NONE
    ks%oep%level = OEP_LEVEL_NONE

    ks%theory_level = KOHN_SHAM_DFT
    parsed_theory_level = .false.

    ! the user knows what he wants, give her that
    if (parse_is_defined(namespace, 'TheoryLevel')) then
      call parse_variable(namespace, 'TheoryLevel', KOHN_SHAM_DFT, ks%theory_level)
      if (.not. varinfo_valid_option('TheoryLevel', ks%theory_level)) call messages_input_error(namespace, 'TheoryLevel')

      parsed_theory_level = .true.
    end if

    ! parse the XC functional
    default = 0

    call get_functional_from_pseudos(pseudo_x_functional, pseudo_c_functional)

    if (ks%theory_level == KOHN_SHAM_DFT .or. ks%theory_level == GENERALIZED_KOHN_SHAM_DFT) then
      if (pseudo_x_functional /= PSEUDO_EXCHANGE_ANY) then
        default = pseudo_x_functional
      else
        select case (space%dim)
        case (3)
          default = XC_LDA_X
        case (2)
          default = XC_LDA_X_2D
        case (1)
          ! TODO: When the support of libxc4 is removed, this should be changed to XC_LDA_X_1D_SOFT
          default = XC_LDA_X_1D
        end select
      end if
    end if

    ASSERT(default >= 0)

    if (ks%theory_level == KOHN_SHAM_DFT .or. ks%theory_level == GENERALIZED_KOHN_SHAM_DFT) then
      if (pseudo_c_functional /= PSEUDO_CORRELATION_ANY) then
        default = default + 1000*pseudo_c_functional
      else
        select case (space%dim)
        case (3)
          default = default + 1000*XC_LDA_C_PZ_MOD
        case (2)
          default = default + 1000*XC_LDA_C_2D_AMGB
        case (1)
          default = default + 1000*XC_LDA_C_1D_CSC
        end select
      end if
    end if

    ASSERT(default >= 0)

    if (.not. parse_is_defined(namespace, 'XCFunctional') &
      .and. (pseudo_x_functional /= PSEUDO_EXCHANGE_ANY .or. pseudo_c_functional /= PSEUDO_CORRELATION_ANY)) then
      call messages_write('Info: the XCFunctional has been selected to match the pseudopotentials', new_line = .true.)
      call messages_write('      used in the calculation.')
      call messages_info()
    end if

    ! The description of this variable can be found in file src/xc/functionals_list.F90
    call parse_variable(namespace, 'XCFunctional', default, val)

    ! the first 3 digits of the number indicate the X functional and
    ! the next 3 the C functional.
    c_id = val / 1000
    x_id = val - c_id*1000

    if ((x_id /= pseudo_x_functional .and. pseudo_x_functional /= PSEUDO_EXCHANGE_ANY) .or. &
      (c_id /= pseudo_c_functional .and. pseudo_c_functional /= PSEUDO_EXCHANGE_ANY)) then
      call messages_write('The XCFunctional that you selected does not match the one used', new_line = .true.)
      call messages_write('to generate the pseudopotentials.')
      call messages_warning(namespace=namespace)
    end if

    ! FIXME: we rarely need this. We should only parse when necessary.

    !%Variable XCKernel
    !%Type integer
    !%Section Hamiltonian::XC
    !%Description
    !% Defines the exchange-correlation kernel. Only LDA kernels are available currently.
    !% The options are the same as <tt>XCFunctional</tt>.
    !% Note: the kernel is only needed for Casida, Sternheimer, or optimal-control calculations.
    !% Defaults:
    !% <br>1D: <tt>lda_x_1d + lda_c_1d_csc</tt>
    !% <br>2D: <tt>lda_x_2d + lda_c_2d_amgb</tt>
    !% <br>3D: <tt>lda_x + lda_c_pz_mod</tt>
    !%Option xc_functional -1
    !% The same functional defined by <tt>XCFunctional</tt>.
    !%End
    call parse_variable(namespace, 'XCKernel', default, val)
    if (-1 == val) then
      ck_id = c_id
      xk_id = x_id
    else
      ck_id = val / 1000
      xk_id = val - ck_id*1000
    end if

    call messages_obsolete_variable(namespace, 'XFunctional', 'XCFunctional')
    call messages_obsolete_variable(namespace, 'CFunctional', 'XCFunctional')


    !%Variable EnablePhotons
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% This variable can be used to enable photons in several types of calculations.
    !% It can be used to activate the one-photon OEP formalism.
    !% In the case of CalculationMode = casida, it enables photon modes as
    !% described in ACS Photonics 2019, 6, 11, 2757-2778.
    !% Finally, if set to yes when solving the ferquency-dependent Sternheimer
    !% equation, the photons are coupled to the electronic subsystem.
    !%End
    call messages_obsolete_variable(namespace, 'OEPPtX', 'EnablePhotons')
    call parse_variable(namespace, 'EnablePhotons', .false., ks%has_photons)

    ! initialize XC modules

    ! This is a bit ugly, theory_level might not be generalized KS or HF now
    ! but it might become generalized KS or HF later. This is safe because it
    ! becomes generalized KS in the cases where the functional is hybrid
    ! and the ifs inside check for both conditions.
    using_hartree_fock = (ks%theory_level == HARTREE_FOCK) &
      .or. (ks%theory_level == GENERALIZED_KOHN_SHAM_DFT .and. family_is_hybrid(ks%xc))
    call xc_init(ks%xc, namespace, space%dim, space%periodic_dim, st%qtot, &
      x_id, c_id, xk_id, ck_id, hartree_fock = using_hartree_fock)

    if (bitand(ks%xc%family, XC_FAMILY_LIBVDWXC) /= 0) then
      call libvdwxc_set_geometry(ks%xc%functional(FUNC_C,1)%libvdwxc, namespace, space, gr)
    end if

    ks%xc_family = ks%xc%family
    ks%xc_flags  = ks%xc%flags

    if (.not. parsed_theory_level) then
      default = KOHN_SHAM_DFT

      ! the functional is a hybrid, use Hartree-Fock as theory level by default
      if (family_is_hybrid(ks%xc) .or. family_is_mgga_with_exc(ks%xc)) then
        default = GENERALIZED_KOHN_SHAM_DFT
      end if

      ! In principle we do not need to parse. However we do it for consistency
      call parse_variable(namespace, 'TheoryLevel', default, ks%theory_level)
      if (.not. varinfo_valid_option('TheoryLevel', ks%theory_level)) call messages_input_error(namespace, 'TheoryLevel')

    end if

    ! In case we need OEP, we need to find which type of OEP it is
    oep_type = -1
    if (family_is_mgga_with_exc(ks%xc)) then
      call messages_experimental('MGGA energy functionals')

      if (accel_is_enabled() .and. (gr%parallel_in_domains .or. st%parallel_in_states .or. st%d%kpt%parallel)) then
        !At the moment this combination produces wrong results
        call messages_not_implemented("MGGA with energy functionals and CUDA+MPI")
      end if

      if (ks%theory_level == KOHN_SHAM_DFT) then
        call messages_experimental("MGGA within the Kohn-Sham scheme")
        ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)
        oep_type = OEP_TYPE_MGGA
      end if
    end if

    call messages_obsolete_variable(namespace, 'NonInteractingElectrons', 'TheoryLevel')
    call messages_obsolete_variable(namespace, 'HartreeFock', 'TheoryLevel')

    ! Due to how the code is made, we need to set this to have theory level other than DFT
    ! correct...
    ks%sic%amaldi_factor = M_ONE

    select case (ks%theory_level)
    case (INDEPENDENT_PARTICLES)

    case (HARTREE)
      call messages_experimental("Hartree theory level")
      if (space%periodic_dim == space%dim) then
        call messages_experimental("Hartree in fully periodic system")
      end if
      if (kpoints%full%npoints > 1) then
        call messages_not_implemented("Hartree with k-points", namespace=namespace)
      end if

    case (HARTREE_FOCK)
      if (kpoints%full%npoints > 1) then
        call messages_experimental("Hartree-Fock with k-points")
      end if

    case (GENERALIZED_KOHN_SHAM_DFT)
      if (kpoints%full%npoints > 1 .and. family_is_hybrid(ks%xc)) then
        call messages_experimental("Hybrid functionals with k-points")
      end if

    case (RDMFT)
      call messages_experimental('RDMFT theory level')

    case (KOHN_SHAM_DFT)

      ! check for SIC
      if (bitand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA) /= 0) then
        call xc_sic_init(ks%sic, namespace, gr, st, mc, space)
      end if

      if (bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        select case (ks%xc%functional(FUNC_X,1)%id)
        case (XC_OEP_X_SLATER)
          if (kpoints%reduced%npoints > 1) then
            call messages_not_implemented("Slater with k-points", namespace=namespace)
          end if
          ks%oep%level = OEP_LEVEL_NONE
        case (XC_OEP_X_FBE)
          if (kpoints%reduced%npoints > 1) then
            call messages_not_implemented("FBE functional with k-points", namespace=namespace)
          end if
          ks%oep%level = OEP_LEVEL_NONE
        case default
          if(oep_type == -1) then ! Else we have a MGGA
            if(ks%has_photons) then
              oep_type = OEP_TYPE_PHOTONS
            else
              oep_type = OEP_TYPE_EXX
            end if
          end if
          call xc_oep_init(ks%oep, namespace, gr, st, mc, space, oep_type)
        end select
      else
        ks%oep%level = OEP_LEVEL_NONE
      end if

      if (bitand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_init(ks%ks_inversion, namespace, gr, ions, st, ks%xc, mc, space, kpoints)
      end if

    end select

    if (st%d%ispin == SPINORS) then
      if (bitand(ks%xc_family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA) /= 0) then
        call messages_not_implemented("GGA with spinors", namespace=namespace)
      end if
      if (bitand(ks%xc_family, XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0) then
        call messages_not_implemented("MGGA with spinors", namespace=namespace)
      end if
    end if

    ks%frozen_hxc = .false.

    call v_ks_write_info(ks, namespace=namespace)

    ks%gr => gr
    ks%calc%calculating = .false.

    !The value of ks%calculate_current is set to false or true by Output
    call current_init(ks%current_calculator, namespace)

    !%Variable VDWCorrection
    !%Type integer
    !%Default no
    !%Section Hamiltonian::XC
    !%Description
    !% (Experimental) This variable selects which van der Waals
    !% correction to apply to the correlation functional.
    !%Option none 0
    !% No correction is applied.
    !%Option vdw_ts 1
    !% The scheme of Tkatchenko and Scheffler, Phys. Rev. Lett. 102
    !% 073005 (2009).
    !%Option vdw_d3 3
    !% The DFT-D3 scheme of S. Grimme, J. Antony, S. Ehrlich, and
    !% S. Krieg, J. Chem. Phys. 132, 154104 (2010).
    !%End
    call parse_variable(namespace, 'VDWCorrection', OPTION__VDWCORRECTION__NONE, ks%vdw_correction)

    if (ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
      call messages_experimental('VDWCorrection')

      select case (ks%vdw_correction)
      case (OPTION__VDWCORRECTION__VDW_TS)

        !%Variable VDWSelfConsistent
        !%Type logical
        !%Default yes
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls whether the VDW correction is applied
        !% self-consistently, the default, or just as a correction to
        !% the total energy. This option only works with vdw_ts.
        !%End
        call parse_variable(namespace, 'VDWSelfConsistent', .true., ks%vdw_self_consistent)

        call vdw_ts_init(ks%vdw_ts, namespace, ions)

      case (OPTION__VDWCORRECTION__VDW_D3)
        ks%vdw_self_consistent = .false.

        if (space%dim /= 3) then
          call messages_write('vdw_d3 can only be used in 3-dimensional systems')
          call messages_fatal(namespace=namespace)
        end if

        do iatom = 1, ions%natoms
          if (.not. species_represents_real_atom(ions%atom(iatom)%species)) then
            call messages_write('vdw_d3 is not implemented when non-atomic species are present')
            call messages_fatal(namespace=namespace)
          end if
        end do

        d3func_def = ''

        ! The list of valid values can be found in 'external_libs/dftd3/core.f90'.
        ! For the moment I include the most common ones.
        if (x_id == OPTION__XCFUNCTIONAL__GGA_X_B88 .and. c_id*1000 == OPTION__XCFUNCTIONAL__GGA_C_LYP) then
          d3func_def = 'b-lyp'
        end if
        if (x_id == OPTION__XCFUNCTIONAL__GGA_X_PBE .and. c_id*1000 == OPTION__XCFUNCTIONAL__GGA_C_PBE) then
          d3func_def = 'pbe'
        end if
        if (x_id == OPTION__XCFUNCTIONAL__GGA_X_PBE_SOL .and. c_id*1000 == OPTION__XCFUNCTIONAL__GGA_C_PBE_SOL) then
          d3func_def = 'pbesol'
        end if
        if (c_id*1000 == OPTION__XCFUNCTIONAL__HYB_GGA_XC_B3LYP) then
          d3func_def = 'b3-lyp'
        end if
        if (c_id*1000 == OPTION__XCFUNCTIONAL__HYB_GGA_XC_PBEH) then
          d3func_def = 'pbe0'
        end if

        !%Variable VDWD3Functional
        !%Type string
        !%Section Hamiltonian::XC
        !%Description
        !% (Experimental) You can use this variable to override the
        !% parametrization used by the DFT-D3 van deer Waals
        !% correction. Normally you need not set this variable, as the
        !% proper value will be selected by Octopus (if available).
        !%
        !% This variable takes a string value, the valid values can
        !% be found in the source file 'external_libs/dftd3/core.f90'.
        !% For example you can use:
        !%
        !%  VDWD3Functional = 'pbe'
        !%
        !%End
        if (parse_is_defined(namespace, 'VDWD3Functional')) call messages_experimental('VDWD3Functional')
        call parse_variable(namespace, 'VDWD3Functional', d3func_def, d3func)

        if (d3func == '') then
          call messages_write('Cannot find  a matching parametrization  of DFT-D3 for the current')
          call messages_new_line()
          call messages_write('XCFunctional.  Please select a different XCFunctional, or select a')
          call messages_new_line()
          call messages_write('functional for DFT-D3 using the <tt>VDWD3Functional</tt> variable.')
          call messages_fatal(namespace=namespace)
        end if

        if (space%periodic_dim /= 0 .and. space%periodic_dim /= 3) then
          call messages_write('For partially periodic systems,  the vdw_d3 interaction is assumed')
          call messages_new_line()
          call messages_write('to be periodic in three dimensions.')
          call messages_warning(namespace=namespace)
        end if

        call dftd3_init(ks%vdw_d3, d3_input, trim(conf%share)//'/dftd3/pars.dat')
        call dftd3_set_functional(ks%vdw_d3, func = d3func, version = 4, tz = .false.)

      case default
        ASSERT(.false.)
      end select

    else
      ks%vdw_self_consistent = .false.
    end if

    if (ks%has_photons) then
      call messages_experimental('EnablePhotons = yes')
      call photon_mode_init(ks%pt, namespace, gr, space%dim, st%qtot)
      write(message(1), '(a,i5,a)') 'Happy to have ', ks%pt%nmodes, ' photon modes with us.'
      call messages_info(1, namespace=namespace)
      call mf_init(ks%pt_mx, ks%gr, st, ions, ks%pt)
    end if


    POP_SUB(v_ks_init)

  contains

    subroutine get_functional_from_pseudos(x_functional, c_functional)
      integer, intent(out) :: x_functional
      integer, intent(out) :: c_functional

      integer :: xf, cf, ispecies
      logical :: warned_inconsistent

      x_functional = PSEUDO_EXCHANGE_ANY
      c_functional = PSEUDO_CORRELATION_ANY

      warned_inconsistent = .false.
      do ispecies = 1, ions%nspecies
        xf = species_x_functional(ions%species(ispecies))
        cf = species_c_functional(ions%species(ispecies))

        if (xf == PSEUDO_EXCHANGE_UNKNOWN .or. cf == PSEUDO_CORRELATION_UNKNOWN) then
          call messages_write("Unknown XC functional for species '"//trim(species_label(ions%species(ispecies)))//"'")
          call messages_warning(namespace=namespace)
          cycle
        end if

        if (x_functional == PSEUDO_EXCHANGE_ANY) then
          x_functional = xf
        else
          if (xf /= x_functional .and. .not. warned_inconsistent) then
            call messages_write('Inconsistent XC functional detected between species')
            call messages_warning(namespace=namespace)
            warned_inconsistent = .true.
          end if
        end if

        if (c_functional == PSEUDO_CORRELATION_ANY) then
          c_functional = cf
        else
          if (cf /= c_functional .and. .not. warned_inconsistent) then
            call messages_write('Inconsistent XC functional detected between species')
            call messages_warning(namespace=namespace)
            warned_inconsistent = .true.
          end if
        end if

      end do

      ASSERT(x_functional /= PSEUDO_EXCHANGE_UNKNOWN)
      ASSERT(c_functional /= PSEUDO_CORRELATION_UNKNOWN)

    end subroutine get_functional_from_pseudos
  end subroutine v_ks_init
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine v_ks_end(ks)
    type(v_ks_t),     intent(inout) :: ks

    PUSH_SUB(v_ks_end)

    select case (ks%vdw_correction)
    case (OPTION__VDWCORRECTION__VDW_TS)
      call vdw_ts_end(ks%vdw_ts)
    end select

    select case (ks%theory_level)
    case (KOHN_SHAM_DFT)
      if (bitand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_end(ks%ks_inversion)
      end if
      if (bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        call xc_oep_end(ks%oep)
      end if
      call xc_end(ks%xc)
    case (HARTREE_FOCK, GENERALIZED_KOHN_SHAM_DFT)
      call xc_end(ks%xc)
    end select

    call xc_sic_end(ks%sic)

    if (ks%has_photons) then
      call photon_mode_end(ks%pt)
      call mf_end(ks%pt_mx)
    end if



    POP_SUB(v_ks_end)
  end subroutine v_ks_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit, namespace)
    type(v_ks_t),                intent(in) :: ks
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(v_ks_write_info)

    call messages_print_with_emphasis(msg="Theory Level", iunit=iunit, namespace=namespace)
    call messages_print_var_option("TheoryLevel", ks%theory_level, iunit=iunit, namespace=namespace)

    select case (ks%theory_level)
    case (HARTREE_FOCK, GENERALIZED_KOHN_SHAM_DFT)
      call messages_info(iunit=iunit, namespace=namespace)
      call xc_write_info(ks%xc, iunit, namespace)

    case (KOHN_SHAM_DFT)
      call messages_info(iunit=iunit, namespace=namespace)
      call xc_write_info(ks%xc, iunit, namespace)

      call messages_info(iunit=iunit, namespace=namespace)

      call xc_sic_write_info(ks%sic, iunit, namespace)
      call xc_oep_write_info(ks%oep, iunit, namespace)
      call xc_ks_inversion_write_info(ks%ks_inversion, iunit, namespace)

    end select

    call messages_print_with_emphasis(iunit=iunit, namespace=namespace)

    POP_SUB(v_ks_write_info)
  end subroutine v_ks_write_info
  ! ---------------------------------------------------------


  !----------------------------------------------------------
  subroutine v_ks_h_setup(namespace, space, gr, ions, ext_partners, st, ks, hm, calc_eigenval, calc_current)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    type(partner_list_t),     intent(in)    :: ext_partners
    type(states_elec_t),      intent(inout) :: st
    type(v_ks_t),             intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    logical,        optional, intent(in)    :: calc_eigenval !< default is true
    logical,        optional, intent(in)    :: calc_current !< default is true

    integer, allocatable :: ind(:)
    integer :: ist, ik
    FLOAT, allocatable :: copy_occ(:)
    logical :: calc_eigenval_
    logical :: calc_current_

    PUSH_SUB(v_ks_h_setup)

    calc_eigenval_ = optional_default(calc_eigenval, .true.)
    calc_current_ = optional_default(calc_current, .true.)
    call states_elec_fermi(st, namespace, gr)
    call density_calc(st, gr, st%rho)
    call v_ks_calc(ks, namespace, space, hm, st, ions, ext_partners, calc_eigenval = calc_eigenval_, calc_current = calc_current_) ! get potentials

    if (st%restart_reorder_occs .and. .not. st%fromScratch) then
      message(1) = "Reordering occupations for restart."
      call messages_info(1, namespace=namespace)

      SAFE_ALLOCATE(ind(1:st%nst))
      SAFE_ALLOCATE(copy_occ(1:st%nst))

      do ik = 1, st%d%nik
        call sort(st%eigenval(:, ik), ind)
        copy_occ(1:st%nst) = st%occ(1:st%nst, ik)
        do ist = 1, st%nst
          st%occ(ist, ik) = copy_occ(ind(ist))
        end do
      end do

      SAFE_DEALLOCATE_A(ind)
      SAFE_DEALLOCATE_A(copy_occ)
    end if

    if (calc_eigenval_) call states_elec_fermi(st, namespace, gr) ! occupations
    call energy_calc_total(namespace, space, hm, gr, st, ext_partners)

    POP_SUB(v_ks_h_setup)
  end subroutine v_ks_h_setup

  ! ---------------------------------------------------------
  subroutine v_ks_calc(ks, namespace, space, hm, st, ions, ext_partners, &
    calc_eigenval, time, calc_energy, calc_current, use_vxc)
    type(v_ks_t),               intent(inout) :: ks
    type(namespace_t),          intent(in)    :: namespace
    type(space_t),              intent(in)    :: space
    type(hamiltonian_elec_t),   intent(inout) :: hm
    type(states_elec_t),        intent(inout) :: st
    type(ions_t),               intent(in)    :: ions
    type(partner_list_t),       intent(in)    :: ext_partners
    logical,          optional, intent(in)    :: calc_eigenval
    FLOAT,            optional, intent(in)    :: time
    logical,          optional, intent(in)    :: calc_energy
    logical,          optional, intent(in)    :: calc_current
    logical,          optional, intent(in)    :: use_vxc

    logical :: calc_current_

    PUSH_SUB(v_ks_calc)

    calc_current_ = optional_default(calc_current, .true.)

    call v_ks_calc_start(ks, namespace, space, hm, st, ions, ext_partners, time, &
      calc_energy, calc_current_, use_vxc)
    call v_ks_calc_finish(ks, hm, namespace, space, st, ext_partners, use_vxc)

    if (optional_default(calc_eigenval, .false.)) then
      call energy_calc_eigenvalues(namespace, hm, ks%gr%der, st)
    end if

    !Update the magnetic constrain
    call magnetic_constrain_update(hm%magnetic_constrain, ks%gr, st%d, ions%space, ions%latt, ions%pos)

    POP_SUB(v_ks_calc)
  end subroutine v_ks_calc

  ! ---------------------------------------------------------

  !> This routine starts the calculation of the Kohn-Sham
  !! potential. The routine v_ks_calc_finish must be called to finish
  !! the calculation. The argument hm is not modified. The argument st
  !! can be modified after the function have been used.
  subroutine v_ks_calc_start(ks, namespace, space, hm, st, ions, ext_partners, time, &
    calc_energy, calc_current, use_vxc)
    type(v_ks_t),              target, intent(inout) :: ks
    type(namespace_t),                 intent(in)    :: namespace
    type(space_t),                     intent(in)    :: space
    type(hamiltonian_elec_t),  target, intent(in)    :: hm !< This MUST be intent(in), changes to hm are done in v_ks_calc_finish.
    type(states_elec_t),               intent(inout) :: st
    type(ions_t),                      intent(in)    :: ions
    type(partner_list_t),              intent(in)    :: ext_partners
    FLOAT,                   optional, intent(in)    :: time
    logical,                 optional, intent(in)    :: calc_energy
    logical,                 optional, intent(in)    :: calc_current
    logical,                 optional, intent(in)    :: use_vxc

    type(profile_t), save :: prof
    logical :: calc_current_

    PUSH_SUB(v_ks_calc_start)

    calc_current_ = optional_default(calc_current, .true.)  &
      .and. ks%calculate_current &
      .and. states_are_complex(st) &
      .or. hamiltonian_elec_needs_current(hm, states_are_real(st))

    call profiling_in(prof, "KOHN_SHAM_CALC")

    ASSERT(.not. ks%calc%calculating)
    ks%calc%calculating = .true.

    if (debug%info) then
      write(message(1), '(a)') 'Debug: Calculating Kohn-Sham potential.'
      call messages_info(1, namespace=namespace)
    end if

    ks%calc%time_present = present(time)
    ks%calc%time = optional_default(time, M_ZERO)

    ks%calc%calc_energy = optional_default(calc_energy, .true.)

    ! If the Hxc term is frozen, there is nothing more to do (WARNING: MISSING ks%calc%energy%intnvxc)
    if (ks%frozen_hxc) then
      if (calc_current_) then
        call states_elec_allocate_current(st, space, ks%gr)
        call current_calculate(ks%current_calculator, namespace, ks%gr, hm, space, st)
      end if

      call profiling_out(prof)
      POP_SUB(v_ks_calc_start)
      return
    end if

    SAFE_ALLOCATE(ks%calc%energy)

    call energy_copy(hm%energy, ks%calc%energy)

    ks%calc%energy%intnvxc = M_ZERO

    nullify(ks%calc%total_density)

    if (ks%theory_level /= INDEPENDENT_PARTICLES .and. ks%sic%amaldi_factor /= M_ZERO) then

      call calculate_density()

      if (poisson_is_async(hm%psolver)) then
        call dpoisson_solve_start(hm%psolver, ks%calc%total_density)
      end if

      if (ks%theory_level /= HARTREE .and. ks%theory_level /= RDMFT) call v_a_xc(hm, use_vxc)
    else
      ks%calc%total_density_alloc = .false.
    end if

    if (calc_current_) then
      call states_elec_allocate_current(st, space, ks%gr)
      call current_calculate(ks%current_calculator, namespace, ks%gr, hm, space, st)
    end if

    nullify(ks%calc%hf_st)
    if (ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK &
      .or. ks%theory_level == RDMFT .or. (ks%theory_level == GENERALIZED_KOHN_SHAM_DFT &
      .and. family_is_hybrid(ks%xc))) then
      SAFE_ALLOCATE(ks%calc%hf_st)
      call states_elec_copy(ks%calc%hf_st, st)

      if (st%parallel_in_states) then
        if (accel_is_enabled()) then
          call messages_write('State parallelization of Hartree-Fock exchange is not supported')
          call messages_new_line()
          call messages_write('when running with OpenCL/CUDA. Please use domain parallelization')
          call messages_new_line()
          call messages_write("or disable acceleration using 'DisableAccel = yes'.")
          call messages_fatal(namespace=namespace)
        end if
        call states_elec_parallel_remote_access_start(ks%calc%hf_st)
      end if
    end if


    ! Calculate the vector potential induced by the electronic current.
    ! WARNING: calculating the self-induced magnetic field here only makes
    ! sense if it is going to be used in the Hamiltonian, which does not happen
    ! now. Otherwise one could just calculate it at the end of the calculation.
    if (hm%self_induced_magnetic) then
      SAFE_ALLOCATE(ks%calc%a_ind(1:ks%gr%np_part, 1:space%dim))
      SAFE_ALLOCATE(ks%calc%b_ind(1:ks%gr%np_part, 1:space%dim))
      call magnetic_induced(namespace, ks%gr, st, hm%psolver, hm%kpoints, ks%calc%a_ind, ks%calc%b_ind)
    end if

    if ((ks%has_photons).and.(ks%calc%time_present)) then
      call mf_calc(ks%pt_mx, ks%gr, st, hm%ions, ks%pt, time)
    end if

    ! if (ks%has_vibrations) then
    !   call vibrations_eph_coup(ks%vib, ks%gr, hm, ions, st)
    ! end if

    call profiling_out(prof)
    POP_SUB(v_ks_calc_start)

  contains

    subroutine calculate_density()
      integer :: ip

      PUSH_SUB(v_ks_calc_start.calculate_density)

      ! get density taking into account non-linear core corrections
      SAFE_ALLOCATE(ks%calc%density(1:ks%gr%np, 1:st%d%nspin))
      call states_elec_total_density(st, ks%gr, ks%calc%density)

      ! Amaldi correction
      if (ks%sic%level == SIC_AMALDI) then
        call lalg_scal(ks%gr%np, st%d%nspin, ks%sic%amaldi_factor, ks%calc%density)
      end if

      nullify(ks%calc%total_density)
      if (allocated(st%rho_core) .or. hm%d%spin_channels > 1) then
        ks%calc%total_density_alloc = .true.

        SAFE_ALLOCATE(ks%calc%total_density(1:ks%gr%np))

        do ip = 1, ks%gr%np
          ks%calc%total_density(ip) = sum(ks%calc%density(ip, 1:hm%d%spin_channels))
        end do

        ! remove non-local core corrections
        if (allocated(st%rho_core)) then
          call lalg_axpy(ks%gr%np, -ks%sic%amaldi_factor, st%rho_core,  ks%calc%total_density)
        end if
      else
        ks%calc%total_density_alloc = .false.
        ks%calc%total_density => ks%calc%density(:, 1)
      end if

      POP_SUB(v_ks_calc_start.calculate_density)
    end subroutine calculate_density

    ! ---------------------------------------------------------
    subroutine v_a_xc(hm, use_vxc)
      type(hamiltonian_elec_t),  intent(in) :: hm
      logical, optional,         intent(in) :: use_vxc

      type(profile_t), save :: prof
      FLOAT :: factor
      integer :: ispin, iatom, idir
      FLOAT, allocatable :: vvdw(:)
      FLOAT :: vdw_stress(1:3, 1:3)
      integer, allocatable :: atnum(:)
      FLOAT :: rlattice(3, 3)

      PUSH_SUB(v_ks_calc_start.v_a_xc)
      call profiling_in(prof, "XC")

      ks%calc%energy%exchange = M_ZERO
      ks%calc%energy%correlation = M_ZERO
      ks%calc%energy%xc_j = M_ZERO
      ks%calc%energy%vdw = M_ZERO

      SAFE_ALLOCATE(ks%calc%vxc(1:ks%gr%np, 1:st%d%nspin))
      ks%calc%vxc = M_ZERO

      if (family_is_mgga_with_exc(hm%xc)) then
        SAFE_ALLOCATE(ks%calc%vtau(1:ks%gr%np, 1:st%d%nspin))
        ks%calc%vtau = M_ZERO
      end if

      if (.not. optional_default(use_vxc, .true.)) then
        call profiling_out(prof)
        POP_SUB(v_ks_calc_start.v_a_xc)
        return
      end if

      ! Get the *local* XC term
      if (ks%calc%calc_energy) then
        if (family_is_mgga_with_exc(hm%xc)) then
          call xc_get_vxc(ks%gr, ks%xc, st, hm%kpoints, hm%psolver, namespace, space, ks%calc%density, st%d%ispin, &
            hm%ions%latt%rcell_volume, ks%calc%vxc, ex = ks%calc%energy%exchange, ec = ks%calc%energy%correlation, &
            deltaxc = ks%calc%energy%delta_xc, vtau = ks%calc%vtau)
        else
          call xc_get_vxc(ks%gr, ks%xc, st, hm%kpoints, hm%psolver, namespace, space, ks%calc%density, st%d%ispin, &
            hm%ions%latt%rcell_volume, ks%calc%vxc, ex = ks%calc%energy%exchange, ec = ks%calc%energy%correlation, &
            deltaxc = ks%calc%energy%delta_xc)
        end if
      else
        if (family_is_mgga_with_exc(hm%xc)) then
          call xc_get_vxc(ks%gr, ks%xc, st, hm%kpoints, hm%psolver, namespace, space, ks%calc%density, &
            st%d%ispin, hm%ions%latt%rcell_volume, ks%calc%vxc, vtau = ks%calc%vtau)
        else
          call xc_get_vxc(ks%gr, ks%xc, st, hm%kpoints, hm%psolver, namespace, space, ks%calc%density, &
            st%d%ispin, hm%ions%latt%rcell_volume, ks%calc%vxc)
        end if
      end if

      ! ADSIC correction
      if (ks%sic%level == SIC_ADSIC) then
        if (family_is_mgga(hm%xc%family)) then
          call messages_not_implemented('ADSIC with MGGAs', namespace=namespace)
        end if
        if (ks%calc%calc_energy) then
          call xc_sic_calc_adsic(ks%sic, namespace, space, ks%gr, st, hm, ks%xc, ks%calc%density, &
            ks%calc%total_density, ks%calc%vxc, ex = ks%calc%energy%exchange, ec = ks%calc%energy%correlation)
        else
          call xc_sic_calc_adsic(ks%sic, namespace, space, ks%gr, st, hm, ks%xc, ks%calc%density, &
            ks%calc%total_density, ks%calc%vxc)
        end if
      end if
      if(ks%sic%level == SIC_PZ_OEP) then
        if (states_are_real(st)) then
          call dxc_oep_calc(ks%sic%oep, namespace, ks%xc, ks%gr, hm, st, space, &
            hm%ions%latt%rcell_volume, ks%calc%energy%exchange, ks%calc%energy%correlation, vxc = ks%calc%vxc)
        else
          call zxc_oep_calc(ks%sic%oep, namespace, ks%xc, ks%gr, hm, st, space, &
            hm%ions%latt%rcell_volume, ks%calc%energy%exchange, ks%calc%energy%correlation, vxc = ks%calc%vxc)
        end if
      end if

      if (ks%theory_level == KOHN_SHAM_DFT) then
        ! The OEP family has to be handled specially
        if (bitand(ks%xc_family, XC_FAMILY_OEP) /= 0 .or. family_is_mgga_with_exc(ks%xc)) then

          if (ks%xc%functional(FUNC_X,1)%id == XC_OEP_X_SLATER) then
            if (states_are_real(st)) then
              call dslater_calc(namespace, ks%gr, space, hm%exxop, st, hm%kpoints, ks%calc%energy%exchange, &
                vxc = ks%calc%vxc)
            else
              call zslater_calc(namespace, ks%gr, space, hm%exxop, st, hm%kpoints, ks%calc%energy%exchange, &
                vxc = ks%calc%vxc)
            end if
          else if (ks%xc%functional(FUNC_X,1)%id == XC_OEP_X_FBE) then
            if (states_are_real(st)) then
              call dx_fbe_calc(namespace, hm%psolver, ks%gr, ks%gr%der, st, ks%calc%energy%exchange, vxc = ks%calc%vxc)
            else
              call zx_fbe_calc(namespace, hm%psolver, ks%gr, ks%gr%der, st, ks%calc%energy%exchange, vxc = ks%calc%vxc)
            end if

          else

            if (states_are_real(st)) then
              call dxc_oep_calc(ks%oep, namespace, ks%xc, ks%gr, hm, st, space, &
                hm%ions%latt%rcell_volume, ks%calc%energy%exchange, ks%calc%energy%correlation, vxc = ks%calc%vxc)
            else
              call zxc_oep_calc(ks%oep, namespace, ks%xc, ks%gr, hm, st, space, &
                hm%ions%latt%rcell_volume, ks%calc%energy%exchange, ks%calc%energy%correlation, vxc = ks%calc%vxc)
            end if
            if (ks%oep%type == OEP_TYPE_PHOTONS) then
              ks%calc%energy%photon_exchange = ks%oep%pt%ex
            end if
          end if

        end if

        if (bitand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
          ! Also treat KS inversion separately (not part of libxc)
          call xc_ks_inversion_calc(ks%ks_inversion, namespace, space, ks%gr, hm, ext_partners, st, vxc = ks%calc%vxc, &
            time = ks%calc%time)
        end if
      end if

      if (ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
        ASSERT(ions%space%dim == 3)

        SAFE_ALLOCATE(vvdw(1:ks%gr%np))
        SAFE_ALLOCATE(ks%calc%vdw_forces(1:ions%space%dim, 1:ions%natoms))

        select case (ks%vdw_correction)

        case (OPTION__VDWCORRECTION__VDW_TS)
          vvdw = CNST(0.0)
          call vdw_ts_calculate(ks%vdw_ts, namespace, ions, ks%gr, st%d%nspin, st%rho, &
            ks%calc%energy%vdw, vvdw, ks%calc%vdw_forces)

        case (OPTION__VDWCORRECTION__VDW_D3)

          SAFE_ALLOCATE(atnum(1:ions%natoms))

          do iatom = 1, ions%natoms
            atnum(iatom) = nint(species_z(ions%atom(iatom)%species))
          end do

          if (space%is_periodic()) then
            ! DFTD3 treats interactions as 3D periodic. In the case of
            ! mixed-periodicity, we set the lenght of the periodic lattice along
            ! the non-periodic dimensions to be much larger than the system
            ! size, such that for all practical purposes the system is treated
            ! as isolated.
            rlattice = ions%latt%rlattice(1:3, 1:3)
            if (space%periodic_dim < 3) then
              message(1) = "Info: Using DFT-D3 van der Walls corrections for a system with mixed-periodicity,"
              message(2) = "      but DFTD3 treats system as fully periodic. Octopus will set the DFT-D3"
              message(3) = "      lattice vectors along non-periodic dimensions to a suitably large value."
              call messages_info(3, namespace=namespace)

              do idir = space%periodic_dim + 1, 3
                rlattice(idir,idir) = (maxval(abs(ions%pos(idir,:))) + M_ONE)*CNST(1000)
              end do
            end if
            call dftd3_pbc_dispersion(ks%vdw_d3, ions%pos, atnum, rlattice, ks%calc%energy%vdw, ks%calc%vdw_forces, &
              vdw_stress)
          else
            call dftd3_dispersion(ks%vdw_d3, ions%pos, atnum, ks%calc%energy%vdw, ks%calc%vdw_forces)
          end if

          SAFE_DEALLOCATE_A(atnum)

        case default
          ASSERT(.false.)

        end select

        if (ks%vdw_self_consistent) then
          do ispin = 1, hm%d%nspin
            ks%calc%vxc(1:ks%gr%np, ispin) = ks%calc%vxc(1:ks%gr%np, ispin) + vvdw(1:ks%gr%np)
          end do
        end if

        SAFE_DEALLOCATE_A(vvdw)

      end if

      if (ks%calc%calc_energy) then
        ! Now we calculate Int[n vxc] = energy%intnvxc
        ks%calc%energy%intnvxc = M_ZERO

        do ispin = 1, hm%d%nspin
          if (ispin <= 2) then
            factor = M_ONE
          else
            factor = M_TWO
          end if
          ks%calc%energy%intnvxc = ks%calc%energy%intnvxc + &
            factor*dmf_dotp(ks%gr, st%rho(:, ispin), ks%calc%vxc(:, ispin), reduce = .false.)
        end do
        if (ks%gr%parallel_in_domains) call ks%gr%allreduce(ks%calc%energy%intnvxc)

        ! MGGA vtau contribution is done after copying vtau to hm%vtau

        if (hm%lda_u_level /= DFT_U_NONE) then
          if (states_are_real(st)) then
            ks%calc%energy%int_dft_u = denergy_calc_electronic(namespace, hm, ks%gr%der, st, terms = TERM_DFT_U)
          else
            ks%calc%energy%int_dft_u = zenergy_calc_electronic(namespace, hm, ks%gr%der, st, terms = TERM_DFT_U)
          end if
        end if
      end if

      call profiling_out(prof)
      POP_SUB(v_ks_calc_start.v_a_xc)
    end subroutine v_a_xc

  end subroutine v_ks_calc_start
  ! ---------------------------------------------------------

  subroutine v_ks_calc_finish(ks, hm, namespace, space, st, ext_partners, use_vxc)
    type(v_ks_t),     target, intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(states_elec_t),      intent(inout) :: st
    type(partner_list_t),     intent(in)    :: ext_partners
    logical,       optional,  intent(in)    :: use_vxc

    integer                           :: ip, ispin
    type(states_elec_t) :: xst !< The states after the application of the Fock operator
    !!                            This is needed to construct the ACE operator
    FLOAT :: exx_energy

    PUSH_SUB(v_ks_calc_finish)

    ASSERT(ks%calc%calculating)
    ks%calc%calculating = .false.

    if (ks%frozen_hxc) then
      POP_SUB(v_ks_calc_finish)
      return
    end if

    !change the pointer to the energy object
    SAFE_DEALLOCATE_A(hm%energy)
    call move_alloc(ks%calc%energy, hm%energy)

    if (hm%self_induced_magnetic) then
      hm%a_ind(1:ks%gr%np, 1:space%dim) = ks%calc%a_ind(1:ks%gr%np, 1:space%dim)
      hm%b_ind(1:ks%gr%np, 1:space%dim) = ks%calc%b_ind(1:ks%gr%np, 1:space%dim)

      SAFE_DEALLOCATE_A(ks%calc%a_ind)
      SAFE_DEALLOCATE_A(ks%calc%b_ind)
    end if

    if (allocated(hm%v_static)) then
      hm%energy%intnvstatic = dmf_dotp(ks%gr, ks%calc%total_density, hm%v_static)
    else
      hm%energy%intnvstatic = M_ZERO
    end if

    if (ks%theory_level == INDEPENDENT_PARTICLES .or. abs(ks%sic%amaldi_factor) <= M_EPSILON) then

      hm%vhxc = M_ZERO
      hm%energy%intnvxc     = M_ZERO
      hm%energy%hartree     = M_ZERO
      hm%energy%exchange    = M_ZERO
      hm%energy%exchange_hf = M_ZERO
      hm%energy%correlation = M_ZERO
    else

      if (ks%theory_level /= HARTREE .and. ks%theory_level /= RDMFT) then
        ! move allocation of vxc from ks%calc to hm
        SAFE_DEALLOCATE_A(hm%vxc)
        call move_alloc(ks%calc%vxc, hm%vxc)

        if (family_is_mgga_with_exc(hm%xc)) then
          call lalg_copy(ks%gr%np, hm%d%nspin, ks%calc%vtau, hm%vtau)
          SAFE_DEALLOCATE_A(ks%calc%vtau)

          ! We need to evaluate the energy after copying vtau to hm%vtau
          if (ks%theory_level == GENERALIZED_KOHN_SHAM_DFT .and. ks%calc%calc_energy) then
            ! MGGA vtau contribution
            if (states_are_real(st)) then
              hm%energy%intnvxc = hm%energy%intnvxc &
                + denergy_calc_electronic(namespace, hm, ks%gr%der, st, terms = TERM_MGGA)
            else
              hm%energy%intnvxc = hm%energy%intnvxc &
                + zenergy_calc_electronic(namespace, hm, ks%gr%der, st, terms = TERM_MGGA)
            end if
          end if
        end if

      else
        hm%vxc = M_ZERO
      end if

      hm%energy%hartree = M_ZERO
      call v_ks_hartree(namespace, ks, space, hm, ext_partners)


      ! Build Hartree + XC potential

      do ip = 1, ks%gr%np
        hm%vhxc(ip, 1) = hm%vxc(ip, 1) + hm%vhartree(ip)
      end do
      if (allocated(hm%vberry)) then
        do ip = 1, ks%gr%np
          hm%vhxc(ip, 1) = hm%vhxc(ip, 1) + hm%vberry(ip, 1)
        end do
      end if

      if (hm%d%ispin > UNPOLARIZED) then
        do ip = 1, ks%gr%np
          hm%vhxc(ip, 2) = hm%vxc(ip, 2) + hm%vhartree(ip)
        end do
        if (allocated(hm%vberry)) then
          do ip = 1, ks%gr%np
            hm%vhxc(ip, 2) = hm%vhxc(ip, 2) + hm%vberry(ip, 2)
          end do
        end if
      end if

      if (hm%d%ispin == SPINORS) then
        do ispin=3, 4
          do ip = 1, ks%gr%np
            hm%vhxc(ip, ispin) = hm%vxc(ip, ispin)
          end do
        end do
      end if

      ! Note: this includes hybrids calculated with the Fock operator instead of OEP
      hm%energy%exchange_hf = M_ZERO
      if ((ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK &
        .or. ks%theory_level == RDMFT .or. ks%theory_level == GENERALIZED_KOHN_SHAM_DFT)) then

        ! swap the states object
        if (associated(hm%exxop%st)) then
          if (hm%exxop%st%parallel_in_states) call states_elec_parallel_remote_access_stop(hm%exxop%st)
          call states_elec_end(hm%exxop%st)
          SAFE_DEALLOCATE_P(hm%exxop%st)
        end if
        if (associated(ks%calc%hf_st) .and. hm%exxop%useACE) then
          if (ks%calc%hf_st%parallel_in_states) call states_elec_parallel_remote_access_stop(ks%calc%hf_st)
        end if

        !At the moment this block is called before the reinit call. This way the LCAO does not call the
        !exchange operator.
        !This should be changed and the CAM parameters should also be obtained from the restart information
        !Maybe the parameters should be mixed too.
        exx_energy = M_ZERO
        if (optional_default(use_vxc, .true.)) then
          if ((ks%theory_level == HARTREE_FOCK .or. ks%theory_level == RDMFT &
            .or. (ks%theory_level == GENERALIZED_KOHN_SHAM_DFT &
            .and. family_is_hybrid(ks%xc))) .and. hm%exxop%useACE) then
            call xst%nullify()
            if (states_are_real(ks%calc%hf_st)) then
              call dexchange_operator_compute_potentials(hm%exxop, namespace, space, ks%gr, &
                ks%calc%hf_st, xst, hm%kpoints, exx_energy)
              call dexchange_operator_ACE(hm%exxop, namespace, ks%gr, ks%calc%hf_st, xst)
            else
              call zexchange_operator_compute_potentials(hm%exxop, namespace, space, ks%gr, &
                ks%calc%hf_st, xst, hm%kpoints, exx_energy)
              if (allocated(hm%hm_base%phase)) then
                call zexchange_operator_ACE(hm%exxop, namespace, ks%gr, ks%calc%hf_st, xst, &
                  hm%hm_base%phase(1:ks%gr%np, ks%calc%hf_st%d%kpt%start:ks%calc%hf_st%d%kpt%end))
              else
                call zexchange_operator_ACE(hm%exxop, namespace, ks%gr, ks%calc%hf_st, xst)
              end if
            end if
            call states_elec_end(xst)
            exx_energy = exx_energy + hm%exxop%singul%energy
          end if
        end if

        select case (ks%theory_level)
        case (GENERALIZED_KOHN_SHAM_DFT)
          if (family_is_hybrid(ks%xc)) then
            ! Add the energy only the ACE case. In the non-ACE case, the exchange energy is not correct
            hm%energy%exchange_hf = hm%energy%exchange_hf + ks%xc%cam_alpha*exx_energy
            call exchange_operator_reinit(hm%exxop, ks%xc%cam_omega, ks%xc%cam_alpha, ks%xc%cam_beta, ks%calc%hf_st)
          end if
        case (HARTREE_FOCK)
          ! Add the energy only the ACE case. In the non-ACE case, the exchange energy is not correct
          hm%energy%exchange_hf = hm%energy%exchange_hf + ks%xc%cam_alpha*exx_energy
          call exchange_operator_reinit(hm%exxop, ks%xc%cam_omega, ks%xc%cam_alpha, ks%xc%cam_beta, ks%calc%hf_st)
        case (HARTREE)
          call exchange_operator_reinit(hm%exxop, M_ZERO, M_ONE, M_ZERO, ks%calc%hf_st)
        case (RDMFT)
          call exchange_operator_reinit(hm%exxop, M_ZERO, M_ONE, M_ZERO, ks%calc%hf_st)
        end select
      end if

    end if

    ! Because of the intent(in) in v_ks_calc_start, we need to update the parameters of hybrids for OEP
    ! here
    if (ks%theory_level == KOHN_SHAM_DFT .and. bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
      if (ks%xc%functional(FUNC_X,1)%id /= XC_OEP_X_SLATER .and. ks%xc%functional(FUNC_X,1)%id /= XC_OEP_X_FBE) then
        call exchange_operator_reinit(hm%exxop, ks%xc%cam_omega, ks%xc%cam_alpha, ks%xc%cam_beta)
      end if
    end if

    if (ks%has_photons) then
      if (associated(ks%pt_mx%vmf)) then
        forall(ip = 1:ks%gr%np) hm%vhxc(ip, 1) = hm%vhxc(ip, 1) + ks%pt_mx%vmf(ip)
        if (hm%d%ispin > UNPOLARIZED) then
          forall(ip = 1:ks%gr%np) hm%vhxc(ip, 2) = hm%vhxc(ip, 2) + ks%pt_mx%vmf(ip)
        end if
      end if
      hm%ep%photon_forces(1:space%dim) = ks%pt_mx%fmf(1:space%dim)
    end if

    if (ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
      hm%ep%vdw_forces = ks%calc%vdw_forces
      SAFE_DEALLOCATE_A(ks%calc%vdw_forces)
    else
      hm%ep%vdw_forces = CNST(0.0)
    end if

    if (ks%calc%time_present .or. hm%time_zero) then
      call hm%update(ks%gr, namespace, space, ext_partners, time = ks%calc%time)
    else
      call hamiltonian_elec_update_pot(hm, ks%gr, accel_copy=.true.)
    end if


    SAFE_DEALLOCATE_A(ks%calc%density)
    if (ks%calc%total_density_alloc) then
      SAFE_DEALLOCATE_P(ks%calc%total_density)
    end if
    nullify(ks%calc%total_density)

    POP_SUB(v_ks_calc_finish)
  end subroutine v_ks_calc_finish

  ! ---------------------------------------------------------
  !
  !> Hartree contribution to the KS potential. This function is
  !! designed to be used by v_ks_calc_finish and it cannot be called
  !! directly.
  !
  subroutine v_ks_hartree(namespace, ks, space, hm, ext_partners)
    type(namespace_t),                intent(in)    :: namespace
    type(v_ks_t),                     intent(inout) :: ks
    type(space_t),                    intent(in)    :: space
    type(hamiltonian_elec_t),         intent(inout) :: hm
    type(partner_list_t),             intent(in)    :: ext_partners

    PUSH_SUB(v_ks_hartree)

    if (.not. poisson_is_async(hm%psolver)) then
      ! solve the Poisson equation
      call dpoisson_solve(hm%psolver, namespace, hm%vhartree, ks%calc%total_density)
    else
      ! The calculation was started by v_ks_calc_start.
      call dpoisson_solve_finish(hm%psolver, hm%vhartree)
    end if

    if (ks%calc%calc_energy) then
      ! Get the Hartree energy
      hm%energy%hartree = M_HALF*dmf_dotp(ks%gr, ks%calc%total_density, hm%vhartree)
    end if

    !> PCM reaction field due to the electronic density
    if(ks%calc%time_present) then
      if(hamiltonian_elec_has_kick(hm)) then
        call pcm_hartree_potential(hm%pcm, space, ks%gr, hm%psolver, ext_partners, hm%vhartree, &
          ks%calc%total_density, hm%energy%pcm_corr, kick=hm%kick, time=ks%calc%time)
      else
        call pcm_hartree_potential(hm%pcm, space, ks%gr, hm%psolver, ext_partners, hm%vhartree, &
          ks%calc%total_density, hm%energy%pcm_corr, time=ks%calc%time)
      end if
    else
      if(hamiltonian_elec_has_kick(hm)) then
        call pcm_hartree_potential(hm%pcm, space, ks%gr, hm%psolver, ext_partners, hm%vhartree, &
          ks%calc%total_density, hm%energy%pcm_corr, kick=hm%kick)
      else
        call pcm_hartree_potential(hm%pcm, space, ks%gr, hm%psolver, ext_partners, hm%vhartree, &
          ks%calc%total_density, hm%energy%pcm_corr)
      end if
    end if

    POP_SUB(v_ks_hartree)
  end subroutine v_ks_hartree
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_freeze_hxc(ks)
    type(v_ks_t), intent(inout) :: ks

    PUSH_SUB(v_ks_freeze_hxc)

    ks%frozen_hxc = .true.

    POP_SUB(v_ks_freeze_hxc)
  end subroutine v_ks_freeze_hxc
  ! ---------------------------------------------------------

  subroutine v_ks_calculate_current(this, calc_cur)
    type(v_ks_t), intent(inout) :: this
    logical,      intent(in)    :: calc_cur

    PUSH_SUB(v_ks_calculate_current)

    this%calculate_current = calc_cur

    POP_SUB(v_ks_calculate_current)
  end subroutine v_ks_calculate_current

  subroutine get_rotation_matrix(dens, alpha, betar, betai)
    FLOAT,  intent(in)  :: dens(:)
    FLOAT,  intent(out) :: alpha, betar, betai

    FLOAT :: mz, mm

    mz = dens(1) - dens(2)

    mm = sqrt(mz**2 + M_FOUR*(dens(3)**2 + dens(4)**2))

    !Fully spin unpolarized system
    if (mm < CNST(1.0e-12)) then
      alpha = M_ONE
      betar = M_ZERO
      betai = M_ZERO
      return
    end if

    alpha = sqrt((mm + abs(mz))/(M_TWO * mm))
    !We find the absolute values of real and imaginary parts of beta
    betar = M_TWO * dens(3) / sqrt(M_TWO * mm * (mm + abs(mz)))
    betai = M_TWO * dens(4) / sqrt(M_TWO * mm * (mm + abs(mz)))

    if (mz < M_ZERO) then
      betar = -betar
      betai = -betai
    end if

  end subroutine get_rotation_matrix

  !Given a matrix in spin space, this routine rotates is according to the rotation
  !matrix R defined by the alpha and beta coefficients
  !rotmat = R mat R^T
  subroutine rotate_to_local(mat, alpha, betar, betai, alpha2, beta2, rot_mat)
    FLOAT,  intent(in)  :: mat(:)
    FLOAT,  intent(in)  :: alpha, betar, betai, alpha2, beta2
    FLOAT,  intent(out) :: rot_mat(:)

    CMPLX :: cross

    rot_mat(1) = alpha2 * mat(1) + beta2 * mat(2) + M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    rot_mat(2) = alpha2 * mat(2) + beta2 * mat(1) - M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    cross = (TOCMPLX(betar, betai))**2 * TOCMPLX(mat(3), -mat(4))
    rot_mat(3) = alpha2 * mat(3) + alpha * betar * (mat(2)-mat(1)) - real(cross)
    rot_mat(4) = alpha2 * mat(4) + alpha * betai * (mat(2)-mat(1)) - aimag(cross)

  end subroutine rotate_to_local

  !Given a matrix in spin space, this routine rotates is according to the rotation
  !matrix R defined by the alpha and beta coefficients
  !rotmat = R^T mat R
  subroutine rotate_to_global(mat, alpha, betar, betai, alpha2, beta2, rot_mat)
    FLOAT,  intent(in)  :: mat(:)
    FLOAT,  intent(in)  :: alpha, betar, betai, alpha2, beta2
    FLOAT,  intent(out) :: rot_mat(:)

    CMPLX :: cross

    rot_mat(1) = alpha2 * mat(1) + beta2 * mat(2) - M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    rot_mat(2) = alpha2 * mat(2) + beta2 * mat(1) + M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    cross = (TOCMPLX(betar, betai))**2 * TOCMPLX(mat(3), -mat(4))
    rot_mat(3) = alpha2 * mat(3) - alpha * betar * (mat(2)-mat(1)) - real(cross)
    rot_mat(4) = alpha2 * mat(4) - alpha * betai * (mat(2)-mat(1)) - aimag(cross)

  end subroutine rotate_to_global



#include "undef.F90"
#include "real.F90"
#include "xc_slater_inc.F90"
#include "x_fbe_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_slater_inc.F90"
#include "x_fbe_inc.F90"

end module v_ks_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
