!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012-2013 M. Gruning, P. Melo, M. Oliveira
!! Copyright (C) 2022 N. Tancogne-Dejean
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

module xc_oep_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use exchange_operator_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lalg_adv_oct_m
  use lattice_vectors_oct_m
  use linear_response_oct_m
  use linear_solver_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use photon_mode_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use scf_tol_oct_m
  use varinfo_oct_m
  use xc_oct_m
  use xc_f03_lib_m
  use xc_functional_oct_m

  implicit none

  private
  public ::                     &
    xc_oep_t,                   &
    xc_oep_init,                &
    xc_oep_end,                 &
    xc_oep_write_info,          &
    dxc_oep_calc,               &
    zxc_oep_calc

  !> the OEP levels
  integer, public, parameter ::  &
    OEP_LEVEL_NONE   = 1,           &
    OEP_LEVEL_KLI    = 3,           &
    OEP_LEVEL_FULL   = 5

  !> Mixing schemes
  integer, public, parameter ::  &
    OEP_MIXING_SCHEME_CONST = 1, &
    OEP_MIXING_SCHEME_BB    = 2, &
    OEP_MIXING_SCHEME_DENS  = 3

  !> The different types of OEP that we can work with
  integer, public, parameter ::  &
    OEP_TYPE_EXX     = 1,          &
    OEP_TYPE_MGGA    = 2,          &
    OEP_TYPE_SIC     = 3,          &
    OEP_TYPE_PHOTONS = 4             ! one-photon OEP

  type xc_oep_t
    private
    integer,              public :: level      !< 0 = no oep, 1 = Slater, 2 = KLI, 4 = full OEP
    FLOAT                        :: mixing     !< how much of the function S(r) to add to vxc in every iteration
    type(lr_t)                   :: lr         !< to solve the equation H psi = b
    type(linear_solver_t)        :: solver
    type(scf_tol_t)              :: scftol
    integer                      :: eigen_n
    integer, allocatable         :: eigen_type(:), eigen_index(:)
    FLOAT                        :: socc, sfact
    FLOAT,   allocatable, public :: vxc(:,:), uxc_bar(:,:,:)
    FLOAT,   allocatable         :: dlxc(:, :, :, :)
    CMPLX,   allocatable         :: zlxc(:, :, :, :)
    integer                      :: mixing_scheme
    type(photon_mode_t),  public :: pt
    type(lr_t)                   :: photon_lr     !< to solve the equation H psi = b
    FLOAT,                public :: norm2ss
    FLOAT,   allocatable         :: vxc_old(:,:), ss_old(:,:)
    integer                      :: noccst
    logical                      :: coc_translation

    integer, public              :: type = -1
  end type xc_oep_t

  type(profile_t), save ::      &
    C_PROFILING_XC_OEP,         &
    C_PROFILING_XC_SIC,         &
    C_PROFILING_OEP_LEVEL_FULL,    &
    C_PROFILING_XC_KLI

contains

  ! ---------------------------------------------------------
  subroutine xc_oep_init(oep, namespace, gr, st, mc, space, oep_type)
    type(xc_oep_t),      intent(out)   :: oep
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(inout) :: gr
    type(states_elec_t), intent(in)    :: st
    type(multicomm_t),   intent(in)    :: mc
    type(space_t),       intent(in)    :: space
    integer,             intent(in)    :: oep_type

    PUSH_SUB(xc_oep_init)

    !%Variable OEPLevel
    !%Type integer
    !%Default oep_kli
    !%Section Hamiltonian::XC
    !%Description
    !% At what level shall <tt>Octopus</tt> handle the optimized effective potential (OEP) equation.
    !%Option oep_none 1
    !% Do not solve OEP equation.
    !%Option oep_kli 3
    !% Krieger-Li-Iafrate (KLI) approximation.
    !% Ref: JB Krieger, Y Li, GJ Iafrate, <i>Phys. Lett. A</i> <b>146</b>, 256 (1990).
    !%Option oep_full 5
    !% (Experimental) Full solution of OEP equation using the Sternheimer approach.
    !% The linear solver will be controlled by the variables in section <tt>Linear Response::Solver</tt>,
    !% and the iterations for OEP by <tt>Linear Response::SCF in LR calculations</tt> and variable
    !% <tt>OEPMixing</tt>. Note that default for <tt>LRMaximumIter</tt> is set to 10.
    !% Ref: S. Kuemmel and J. Perdew, <i>Phys. Rev. Lett.</i> <b>90</b>, 043004 (2003).
    !%End
    call messages_obsolete_variable(namespace, 'OEP_Level', 'OEPLevel')
    call parse_variable(namespace, 'OEPLevel', OEP_LEVEL_KLI, oep%level)
    if (.not. varinfo_valid_option('OEPLevel', oep%level)) call messages_input_error(namespace, 'OEPLevel')

    if (oep%level == OEP_LEVEL_NONE) then
      POP_SUB(xc_oep_init)
      return
    end if

    oep%type = oep_type

    if (oep%type == OEP_TYPE_PHOTONS) then
      call messages_experimental("EnablePhotons = yes")
      call photon_mode_init(oep%pt, namespace, gr, space%dim, st%qtot)
      if (oep%pt%nmodes > 1) then
        call messages_not_implemented('Photon OEP for more than one photon mode')
      end if
      if (oep%level == OEP_LEVEL_FULL .and. st%d%nspin /= UNPOLARIZED) then
        call messages_not_implemented('Spin-polarized calculations with photon OEP')
      end if

      if (states_are_complex(st)) then
        call messages_not_implemented('Photon OEP with complex wavefunctions')
      end if

      if (st%d%nik > st%d%ispin .and. oep%level == OEP_LEVEL_FULL) then
        call messages_not_implemented("Full OEP for periodic systems", namespace=namespace)
      end if
      if (st%d%nik > st%d%ispin .and. st%d%ispin==SPINORS) then
        call messages_not_implemented("OEP for periodic systems with spinors", namespace=namespace)
      end if
    end if

    if (oep%level == OEP_LEVEL_FULL) then

      if (st%d%nspin == SPINORS) then
        call messages_not_implemented("Full OEP with spinors", namespace=namespace)
      end if

      call messages_experimental("Full OEP")

      !%Variable OEPMixing
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::XC
      !%Description
      !% The linear mixing factor used to solve the Sternheimer
      !% equation in the full OEP procedure.
      !%End
      call messages_obsolete_variable(namespace, 'OEP_Mixing', 'OEPMixing')
      call parse_variable(namespace, 'OEPMixing', M_ONE, oep%mixing)

      !%Variable OEPMixingScheme
      !%Type integer
      !%Default oep_mixing_scheme_const
      !%Section Hamiltonian::XC
      !%Description
      !% Different Mixing Schemes are possible
      !%Option oep_mixing_scheme_const 1
      !% Use a constant
      !% Reference: S. Kuemmel and J. Perdew, <i>Phys. Rev. Lett.</i> <b>90</b>, 4, 043004 (2003)
      !%Option oep_mixing_scheme_bb 2
      !% Use the Barzilai-Borwein (BB) Method
      !% Reference: T. W. Hollins, S. J. Clark, K. Refson, and N. I. Gidopoulos,
      !% <i>Phys. Rev. B</i> <b>85</b>, 235126 (2012)
      !%Option oep_mixing_scheme_dens 3
      !% Use the inverse of the electron density
      !% Reference: S. Kuemmel and J. Perdew, <i>Phys. Rev. B</i> <b>68</b>, 035103 (2003)
      !%End
      call parse_variable(namespace, 'OEPMixingScheme', OEP_MIXING_SCHEME_CONST, oep%mixing_scheme)

      if (oep%mixing_scheme == OEP_MIXING_SCHEME_BB) then
        SAFE_ALLOCATE(oep%vxc_old(1:gr%np,st%d%ispin))
        SAFE_ALLOCATE(oep%ss_old(1:gr%np,st%d%ispin))
        oep%vxc_old = M_ZERO
        oep%ss_old = M_ZERO
      end if

      oep%norm2ss = M_ZERO
    end if

    ! obtain the spin factors
    call xc_oep_SpinFactor(oep, st%d%nspin)

    ! This variable will keep vxc across iterations
    if ((st%d%ispin == SPINORS) .or. oep%level == OEP_LEVEL_FULL) then
      SAFE_ALLOCATE(oep%vxc(1:gr%np,st%d%nspin))
    else
      SAFE_ALLOCATE(oep%vxc(1:gr%np,1:min(st%d%nspin, 2)))
    end if
    oep%vxc = M_ZERO

    !%Variable KLIPhotonCOC
    !%Type logical
    !%Default .false.
    !%Section Hamiltonian::XC
    !%Description
    !% Activate the center of charge translation of the electric dipole operator which should avoid the dependence of the photon KLI on an permanent dipole.
    !%End

    ! when performing full OEP, we need to solve a linear equation
    if ((oep%level == OEP_LEVEL_FULL).or.(oep%type == OEP_TYPE_PHOTONS)) then
      call scf_tol_init(oep%scftol, namespace, st%qtot, def_maximumiter=10)
      call linear_solver_init(oep%solver, namespace, gr, states_are_real(st), mc, space)
      call lr_init(oep%lr)
      if (oep%type == OEP_TYPE_PHOTONS) then
        call lr_init(oep%photon_lr)
        call parse_variable(namespace, 'KLIPhotonCOC', .false., oep%coc_translation)
      end if
    end if

    ! the linear equation has to be more converged if we are to attain the required precision
    !oep%lr%conv_abs_dens = oep%lr%conv_abs_dens / (oep%mixing)

    if (st%d%nspin == SPINORS) then
      call messages_experimental("OEP with spinors")
    end if

    if (st%d%kpt%parallel .and. oep%level == OEP_LEVEL_FULL) then
      call messages_experimental("OEP parallel in spin/k-points", namespace=namespace)
    end if

    POP_SUB(xc_oep_init)
  end subroutine xc_oep_init


  ! ---------------------------------------------------------
  subroutine xc_oep_end(oep)
    type(xc_oep_t), intent(inout) :: oep

    PUSH_SUB(xc_oep_end)

    if (oep%level /= OEP_LEVEL_NONE) then
      SAFE_DEALLOCATE_A(oep%vxc)
      if (oep%level == OEP_LEVEL_FULL .or. oep%type == OEP_TYPE_PHOTONS) then
        call lr_dealloc(oep%lr)
        call linear_solver_end(oep%solver)
      end if
      if (oep%type == OEP_TYPE_PHOTONS) then
        call lr_dealloc(oep%photon_lr)
        call photon_mode_end(oep%pt)
      end if
      if (oep%level == OEP_LEVEL_FULL .and. oep%mixing_scheme == OEP_MIXING_SCHEME_BB) then
        SAFE_DEALLOCATE_A(oep%vxc_old)
        SAFE_DEALLOCATE_A(oep%ss_old)
      end if
    end if

    POP_SUB(xc_oep_end)
  end subroutine xc_oep_end


  ! ---------------------------------------------------------
  subroutine xc_oep_write_info(oep, iunit, namespace)
    type(xc_oep_t),              intent(in) :: oep
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    if (oep%level == OEP_LEVEL_NONE) return

    PUSH_SUB(xc_oep_write_info)

    call messages_print_var_option('OEPLevel', oep%level, iunit=iunit, namespace=namespace)

    POP_SUB(xc_oep_write_info)
  end subroutine xc_oep_write_info


  ! ---------------------------------------------------------
  !> A couple of auxiliary functions for oep
  ! ---------------------------------------------------------
  subroutine xc_oep_SpinFactor(oep, nspin)
    type(xc_oep_t), intent(inout) :: oep
    integer,        intent(in)    :: nspin

    PUSH_SUB(xc_oep_SpinFactor)

    select case (nspin)
    case (1) ! we need to correct or the spin occupancies
      oep%socc  = M_HALF
      oep%sfact = M_TWO
    case (2, 4)
      oep%socc  = M_ONE
      oep%sfact = M_ONE
    case default ! cannot handle any other case
      ASSERT(.false.)
    end select

    POP_SUB(xc_oep_SpinFactor)
  end subroutine xc_oep_SpinFactor


  ! ---------------------------------------------------------
  subroutine xc_oep_AnalyzeEigen(oep, st, is)
    type(xc_oep_t),       intent(inout) :: oep
    type(states_elec_t),  intent(in)    :: st
    integer,              intent(in)    :: is

    integer  :: ist
    FLOAT :: max_eigen
    FLOAT, allocatable :: eigenval(:), occ(:)

    PUSH_SUB(xc_oep_AnalyzeEigen)

    SAFE_ALLOCATE(eigenval(1:st%nst))
    SAFE_ALLOCATE     (occ(1:st%nst))
    eigenval = M_ZERO
    occ = M_ZERO

    do ist = st%st_start, st%st_end
      eigenval(ist) = st%eigenval(ist, is)
      occ(ist) = st%occ(ist, is)
    end do

    if (st%parallel_in_states) then
      call st%mpi_grp%barrier()
      do ist = 1, st%nst
        call st%mpi_grp%bcast(eigenval(ist), 1, MPI_FLOAT, st%node(ist))
        call st%mpi_grp%bcast(occ(ist), 1, MPI_FLOAT, st%node(ist))
      end do
    end if

    ! find out the top occupied state, to correct for the asymptotics
    ! of the potential
    max_eigen = CNST(-1e30)
    do ist = 1, st%nst
      if ((occ(ist) > M_EPSILON).and.(eigenval(ist) > max_eigen)) then
        max_eigen = eigenval(ist)
      end if
    end do

    oep%eigen_n = 1
    do ist = 1, st%nst
      if (occ(ist) > M_EPSILON) then
        ! criterion for degeneracy
        if (abs(eigenval(ist)-max_eigen) <= CNST(1e-3)) then
          oep%eigen_type(ist) = 2
        else
          oep%eigen_type(ist) = 1
          oep%eigen_index(oep%eigen_n) = ist
          oep%eigen_n = oep%eigen_n + 1
        end if
      else
        oep%eigen_type(ist) = 0
      end if
    end do
    oep%eigen_n = oep%eigen_n - 1

    ! find how many states are occupied.
    oep%noccst = 0
    do ist = 1, st%nst
      if (st%occ(ist, is) > M_EPSILON) oep%noccst = ist
    end do

    SAFE_DEALLOCATE_A(eigenval)
    SAFE_DEALLOCATE_A(occ)
    POP_SUB(xc_oep_AnalyzeEigen)
  end subroutine xc_oep_AnalyzeEigen


#include "xc_kli_pauli_inc.F90"
#include "xc_oep_qed_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "xc_kli_inc.F90"
#include "xc_oep_sic_inc.F90"
#include "xc_oep_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_kli_inc.F90"
#include "xc_oep_sic_inc.F90"
#include "xc_oep_inc.F90"

end module xc_oep_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
