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

!> @brief This module implements a unit-test like runmode for Octopus.
!!
#include "global.h"

module test_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use box_oct_m
  use box_factory_oct_m
  use calc_mode_par_oct_m
  use cgal_polyhedra_oct_m
  use current_oct_m
  use clock_oct_m
  use debug_oct_m
  use density_oct_m
  use derivatives_oct_m
  use exponential_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use helmholtz_decomposition_m
  use helmholtz_decomposition_test_m
  use iihash_oct_m
  use ion_interaction_oct_m
  use iso_c_binding
  use io_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use mesh_interpolation_oct_m
  use messages_oct_m
  use maxwell_oct_m
  use mpi_oct_m
  use mpi_test_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use orbitalbasis_oct_m
  use orbitalset_oct_m
  use parser_oct_m
  use poisson_oct_m
  use preconditioners_oct_m
  use profiling_oct_m
  use projector_oct_m
  use regridding_oct_m
  use sihash_oct_m
  use solvers_oct_m
  use space_oct_m
  use sphash_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use states_mxll_oct_m
  use submesh_oct_m
  use subspace_oct_m
  use electrons_oct_m
  use types_oct_m
  use v_ks_oct_m
  use wfs_elec_oct_m
  use unit_system_oct_m
  use io_function_oct_m

  implicit none

  type test_parameters_t
    private
    integer :: type
    integer :: repetitions
    integer :: min_blocksize
    integer :: max_blocksize
  end type test_parameters_t

  !Auxiliary quantities needed by the linear solver
  FLOAT :: shift_aux
  type(derivatives_t), pointer :: der_aux  => null()
  type(preconditioner_t)       :: prec_aux

  public :: test_run

contains

  ! ---------------------------------------------------------
  subroutine test_run(namespace)
    type(namespace_t),       intent(in)    :: namespace

    type(test_parameters_t) :: param
    integer :: test_mode

    PUSH_SUB(test_run)

    call messages_obsolete_variable(namespace, 'WhichTest', 'TestMode')

    !%Variable TestMode
    !%Type integer
    !%Default hartree
    !%Section Calculation Modes::Test
    !%Description
    !% Decides what kind of test should be performed.
    !%Option hartree 1
    !% Tests the Poisson solvers used to calculate the Hartree potential.
    !%Option derivatives 2
    !% Tests and benchmarks the implementation of the finite-difference operators, used to calculate derivatives.
    !%Option orthogonalization 3
    !% Tests the implementation of the orthogonalization routines.
    !%Option interpolation 4
    !% Test the interpolation routines.
    !%Option ion_interaction 5
    !% Tests the ion-ion interaction routines.
    !%Option projector 6
    !% Tests the code that applies the nonlocal part of the pseudopotentials
    !% in case of spin-orbit coupling
    !%Option dft_u 7
    !% Tests the DFT+U part of the code for projections on the basis.
    !%Option hamiltonian_apply 8
    !% Tests the application of the Hamiltonian, or a part of it
    !%Option density_calc 9
    !% Calculation of the density.
    !%Option exp_apply 10
    !% Tests the exponential of the Hamiltonian
    !%Option boundaries 11
    !% Tests the boundaries conditions
    !%Option subspace_diag 12
    !% Tests the subspace diagonalization
    !%Option batch_ops 13
    !% Tests the batch operations
    !%Option clock 18
    !% Tests for clock
    !%Option linear_solver 19
    !% Tests the linear solvers
    !%Option cgal 20
    !% Tests for cgal interface
    !%Option dense_eigensolver 21
    !% Tests for dense eigensolvers (especially parallel ones)
    !%Option grid_interpolation 22
    !% Tests for grid interpolation and multigrid methods.
    !%Option iihash 23
    !% Tests for the integer-integer hash table.
    !%Option sihash 24
    !% Tests for the string-integer hash table.
    !%Option sphash 25
    !% Tests for the string-polymorphic hash table.
    !%Option mpiwrappers 26
    !% Tests for the MPI wrappers with large integer displacements.
    !%Option regridding 27
    !% Tests the regridding between two different grids.
    !%Option helmholtz_decomposition 28
    !% Test for the Helmholtz decomposition subroutines
    !%Option vecpot_analytical 29
    !% Tests analytically the vector potential from B field.
    !%Option current_density 30
    !% Tests the different contributions to the total electronic current density
    !%End
    call parse_variable(namespace, 'TestMode', OPTION__TESTMODE__HARTREE, test_mode)

    call messages_obsolete_variable(namespace, 'TestDerivatives', 'TestType')
    call messages_obsolete_variable(namespace, 'TestOrthogonalization', 'TestType')

    !%Variable TestType
    !%Type integer
    !%Default all
    !%Section Calculation Modes::Test
    !%Description
    !% Decides on what type of values the test should be performed.
    !%Option real 1
    !% Test for double-precision real functions.
    !%Option complex 2
    !%Option all 3
    !% Tests for double-precision real and complex functions.
    !%End
    call parse_variable(namespace, 'TestType', OPTION__TESTTYPE__ALL, param%type)
    if (param%type < 1 .or. param%type > 5) then
      message(1) = "Invalid option for TestType."
      call messages_fatal(1, only_root_writes = .true., namespace=namespace)
    end if

    !%Variable TestRepetitions
    !%Type integer
    !%Default 1
    !%Section Calculation Modes::Test
    !%Description
    !% This variable controls the behavior of oct-test for performance
    !% benchmarking purposes. It sets the number of times the
    !% computational kernel of a test will be executed, in order to
    !% provide more accurate timings.
    !%
    !% Currently this variable is used by the <tt>hartree_test</tt>,
    !% <tt>derivatives</tt>, and <tt>projector</tt> tests.
    !%End
    call parse_variable(namespace, 'TestRepetitions', 1, param%repetitions)

    !%Variable TestMinBlockSize
    !%Type integer
    !%Default 1
    !%Section Calculation Modes::Test
    !%Description
    !% Some tests can work with multiple blocksizes, in this case of
    !% range of blocksizes will be tested. This variable sets the lower
    !% bound of that range.
    !%
    !% Currently this variable is only used by the derivatives test.
    !%End
    call parse_variable(namespace, 'TestMinBlockSize', 1, param%min_blocksize)

    !%Variable TestMaxBlockSize
    !%Type integer
    !%Default 128
    !%Section Calculation Modes::Test
    !%Description
    !% Some tests can work with multiple blocksizes, in this case of
    !% range of blocksizes will be tested. This variable sets the lower
    !% bound of that range.
    !%
    !% Currently this variable is only used by the derivatives test.
    !%End
    call parse_variable(namespace, 'TestMaxBlockSize', 128, param%max_blocksize)

    call messages_print_with_emphasis(msg="Test mode", namespace=namespace)
    call messages_print_var_option("TestMode", test_mode, namespace=namespace)
    call messages_print_var_option("TestType", param%type, namespace=namespace)
    call messages_print_var_value("TestRepetitions", param%repetitions, namespace=namespace)
    call messages_print_var_value("TestMinBlockSize", param%min_blocksize, namespace=namespace)
    call messages_print_var_value("TestMaxBlockSize", param%max_blocksize, namespace=namespace)
    call messages_print_with_emphasis(namespace=namespace)

    select case (test_mode)
    case (OPTION__TESTMODE__HARTREE)
      call test_hartree(param, namespace)
    case (OPTION__TESTMODE__DERIVATIVES)
      call test_derivatives(param, namespace)
    case (OPTION__TESTMODE__ORTHOGONALIZATION)
      call test_orthogonalization(param, namespace)
    case (OPTION__TESTMODE__INTERPOLATION)
      call test_interpolation(param, namespace)
    case (OPTION__TESTMODE__ION_INTERACTION)
      call test_ion_interaction(namespace)
    case (OPTION__TESTMODE__PROJECTOR)
      call test_projector(param, namespace)
    case (OPTION__TESTMODE__DFT_U)
      call test_dft_u(param, namespace)
    case (OPTION__TESTMODE__HAMILTONIAN_APPLY)
      call test_hamiltonian(param, namespace)
    case (OPTION__TESTMODE__DENSITY_CALC)
      call test_density_calc(param, namespace)
    case (OPTION__TESTMODE__EXP_APPLY)
      call test_exponential(param, namespace)
    case (OPTION__TESTMODE__BOUNDARIES)
      call test_boundaries(param, namespace)
    case (OPTION__TESTMODE__SUBSPACE_DIAG)
      call test_subspace_diagonalization(param, namespace)
    case (OPTION__TESTMODE__BATCH_OPS)
      call test_batch_ops(param, namespace)
    case (OPTION__TESTMODE__CLOCK)
      call test_clock()
    case (OPTION__TESTMODE__LINEAR_SOLVER)
      call test_linear_solver(namespace)
    case (OPTION__TESTMODE__CGAL)
      call test_cgal()
    case (OPTION__TESTMODE__DENSE_EIGENSOLVER)
      call test_dense_eigensolver()
    case (OPTION__TESTMODE__GRID_INTERPOLATION)
      call test_grid_interpolation()
    case (OPTION__TESTMODE__IIHASH)
      call test_iihash()
    case (OPTION__TESTMODE__SIHASH)
      call test_sihash()
    case (OPTION__TESTMODE__SPHASH)
      call test_sphash(namespace)
    case (OPTION__TESTMODE__MPIWRAPPERS)
      call test_mpiwrappers()
    case (OPTION__TESTMODE__REGRIDDING)
      call test_regridding(namespace)
    case (OPTION__TESTMODE__HELMHOLTZ_DECOMPOSITION)
      call test_helmholtz_decomposition(namespace)
    case (OPTION__TESTMODE__VECPOT_ANALYTICAL)
      call test_vecpot_analytical(namespace)
    case (OPTION__TESTMODE__CURRENT_DENSITY)
      call test_current_density(namespace)
    end select

    POP_SUB(test_run)
  end subroutine test_run

  ! ---------------------------------------------------------
  subroutine test_hartree(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_hartree)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)
    call poisson_test(sys%hm%psolver, sys%space, sys%gr, sys%ions%latt, namespace, param%repetitions)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_hartree)
  end subroutine test_hartree

  ! ---------------------------------------------------------
  subroutine test_helmholtz_decomposition(namespace)
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t),      pointer :: sys
    type(helmholtz_decomposition_t) :: helmholtz

    PUSH_SUB(test_helmholtz_decomposition)

    ! First of all we have to initialize the grid and the poisson solver
    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call helmholtz%init(namespace, sys%gr, sys%mc, sys%space)

    ! Then we have to initialize the exact fields
    write(message(1),'(a)') "Helmholtz decomposition: beginning Hertzian dipole test"
    call messages_info(1, namespace=namespace)
    call hertzian_dipole_test(helmholtz, sys%gr, sys%namespace, sys%space)

    write(message(1),'(a)') "Helmholtz decomposition: beginning Gaussian test"
    call messages_info(1, namespace=namespace)
    call gaussian_test(helmholtz, sys%gr, sys%namespace, sys%space)

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_helmholtz_decomposition)
  end subroutine test_helmholtz_decomposition

  ! ---------------------------------------------------------
  subroutine test_linear_solver(namespace)
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    FLOAT, allocatable :: rho(:), x(:)
    FLOAT :: center(MAX_DIM)
    FLOAT :: rr, alpha, beta, res
    integer :: ip, iter

    FLOAT, parameter :: threshold = CNST(1e-7)

    PUSH_SUB(test_linear_solver)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    ! We need to set up some auxiliary quantities called by the linear solver
    call mesh_init_mesh_aux(sys%gr)
    ! Shift of the linear equation
    shift_aux = CNST(0.25)
    ! Preconditioner used for the QMR algorithm
    call preconditioner_init(prec_aux, namespace, sys%gr, sys%mc, sys%space)
    ! Derivative object needed
    call set_der_aux(sys%gr%der)

    ! Here we put a Gaussian as the right-hand side of the linear solver
    ! Values are taken from the poisson_test routine
    alpha = M_FOUR * sys%gr%spacing(1)
    beta = M_ONE / (alpha**sys%space%dim * sqrt(M_PI)**sys%space%dim)
    ! The Gaussian is centered around the origin
    center = M_ZERO

    SAFE_ALLOCATE(rho(1:sys%gr%np))
    rho = M_ZERO
    do ip = 1, sys%gr%np
      call mesh_r(sys%gr, ip, rr, origin = center(:))
      rho(ip) = beta*exp(-(rr/alpha)**2)
    end do

    SAFE_ALLOCATE(x(1:sys%gr%np))

    !Test the CG linear solver
    x = M_ZERO
    iter = 1000
    call dconjugate_gradients(sys%gr%np, x, rho, laplacian_op, dmf_dotp_aux, iter, res, threshold)
    write(message(1),'(a,i6,a)')  "Info: CG converged with ", iter, " iterations."
    write(message(2),'(a,e14.6)')    "Info: The residue is ", res
    write(message(3),'(a,e14.6)') "Info: Norm solution CG ", dmf_nrm2(sys%gr, x)
    call messages_info(3, namespace=namespace)

    call preconditioner_end(prec_aux)
    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_linear_solver)
  contains

    subroutine set_der_aux(der)
      type(derivatives_t), target, intent(in) :: der
      PUSH_SUB(test_linear_solver.set_der_aux)
      der_aux => der
      POP_SUB(test_linear_solver.set_der_aux)
    end subroutine set_der_aux

    ! ---------------------------------------------------------
    !> Computes \f$ Hx = (-\Laplacian + shift) x \f$
    subroutine laplacian_op(x, hx)
      FLOAT,            intent(in)    :: x(:)
      FLOAT, contiguous,intent(out)   :: Hx(:)

      FLOAT, allocatable :: tmpx(:)

      ASSERT(associated(mesh_aux))

      SAFE_ALLOCATE(tmpx(1:mesh_aux%np_part))
      call lalg_copy(mesh_aux%np, x, tmpx)
      call dderivatives_lapl(der_aux, tmpx, Hx)
      call lalg_scal(mesh_aux%np, -M_ONE, hx)
      call lalg_axpy(mesh_aux%np, shift_aux, x, hx)
      SAFE_DEALLOCATE_A(tmpx)

    end subroutine laplacian_op

  end subroutine test_linear_solver


  ! ---------------------------------------------------------
  subroutine test_projector(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(wfs_elec_t) :: epsib
    integer :: itime
    CMPLX, allocatable :: psi(:, :)

    PUSH_SUB(test_projector)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the nonlocal part of the pseudopotential with SOC')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr, wfs_type = TYPE_CMPLX)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)

    ! Initialize external potential
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%space, sys%gr, sys%ions, sys%ext_partners, sys%st)

    call sys%st%group%psib(1, 1)%copy_to(epsib)

    call batch_set_zero(epsib)

    do itime = 1, param%repetitions
      call zproject_psi_batch(sys%gr, sys%gr%der%boundaries, sys%hm%ep%proj,  &
        sys%hm%ep%natoms, 2, sys%st%group%psib(1, 1), epsib)
    end do

    SAFE_ALLOCATE(psi(1:sys%gr%np, 1:sys%st%d%dim))
    do itime = 1, epsib%nst
      call batch_get_state(epsib, itime, sys%gr%np, psi)
      write(message(1),'(a,i1,3x, f12.6)') "Norm state  ", itime, zmf_nrm2(sys%gr, 2, psi)
      call messages_info(1, namespace=sys%namespace)
    end do
    SAFE_DEALLOCATE_A(psi)

    call epsib%end()
    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_projector)
  end subroutine test_projector

  ! ---------------------------------------------------------
  subroutine test_dft_u(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(wfs_elec_t) :: epsib, epsib2
    integer :: itime, ist
    type(orbitalbasis_t) :: basis
    FLOAT, allocatable :: ddot(:,:,:), dweight(:,:)
    CMPLX, allocatable :: zdot(:,:,:), zweight(:,:)

    PUSH_SUB(test_dft_u)

    call calc_mode_par_unset_parallelization(P_STRATEGY_STATES)
    call calc_mode_par_unset_parallelization(P_STRATEGY_KPOINTS)
    call calc_mode_par_set_parallelization(P_STRATEGY_DOMAINS, default = .true.)

    call messages_write('Info: Testing some DFT+U routines')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)
    if (sys%st%d%pack_states) call sys%st%pack()

    call sys%st%group%psib(1, 1)%copy_to(epsib2, copy_data = .true.)

    ! We set the phase of the batch if needed
    if (.not. allocated(sys%hm%hm_base%phase)) then
      call sys%st%group%psib(1, 1)%copy_to(epsib, copy_data = .true.)
    else
      call sys%st%group%psib(1, 1)%copy_to(epsib)
      call hamiltonian_elec_base_phase(sys%hm%hm_base, sys%gr, sys%gr%np, &
        .false., epsib, src=sys%st%group%psib(1, 1))
      epsib2%has_phase = .true.
    end if

    ! Initialize the orbital basis
    call orbitalbasis_init(basis, sys%namespace, sys%space%periodic_dim)
    if (states_are_real(sys%st)) then
      call dorbitalbasis_build(basis, sys%namespace, sys%ions, sys%gr, sys%st%d%kpt, sys%st%d%dim, &
        allocated(sys%hm%hm_base%phase), .false., .false.)
      SAFE_ALLOCATE(dweight(1:basis%orbsets(1)%norbs, 1:epsib%nst_linear))
      SAFE_ALLOCATE(ddot(1:sys%st%d%dim, 1:basis%orbsets(1)%norbs, 1:epsib%nst))
    else
      call zorbitalbasis_build(basis, sys%namespace, sys%ions, sys%gr, sys%st%d%kpt, sys%st%d%dim, &
        allocated(sys%hm%hm_base%phase), .false., .false.)
      call orbitalset_update_phase(basis%orbsets(1), sys%space%dim, sys%st%d%kpt, sys%kpoints, (sys%st%d%ispin==SPIN_POLARIZED))
      SAFE_ALLOCATE(zweight(1:basis%orbsets(1)%norbs, 1:epsib%nst_linear))
      SAFE_ALLOCATE(zdot(1:sys%st%d%dim, 1:basis%orbsets(1)%norbs, 1:epsib%nst))

      !We set the phase of the orbitals if needed
      if (allocated(sys%hm%hm_base%phase)) then
        call orbitalset_update_phase(basis%orbsets(1), sys%space%dim, sys%st%d%kpt, sys%kpoints, &
          (sys%st%d%ispin==SPIN_POLARIZED))
      end if
    end if

    do itime = 1, param%repetitions
      call batch_set_zero(epsib2)
      if (states_are_real(sys%st)) then
        dweight = M_ONE
        ddot = M_ZERO
        call dorbitalset_get_coeff_batch(basis%orbsets(1), sys%st%d%dim, sys%st%group%psib(1, 1), ddot)
        call dorbitalset_add_to_batch(basis%orbsets(1), sys%st%d%dim, epsib2, dweight)
      else
        zweight = M_ONE
        zdot = M_ZERO
        call zorbitalset_get_coeff_batch(basis%orbsets(1), sys%st%d%dim, epsib, zdot)
        call zorbitalset_add_to_batch(basis%orbsets(1), sys%st%d%dim, epsib2, zweight)
      end if
    end do

    if (epsib%is_packed()) then
      call epsib%do_unpack(force = .true.)
    end if

    do ist = 1, epsib%nst
      if (states_are_real(sys%st)) then
        write(message(1),'(a,i2,3x,e13.6)') "Dotp state ", ist, ddot(1,1,ist)
      else
        write(message(1),'(a,i2,2(3x,e13.6))') "Dotp state ", ist, zdot(1,1,ist)
      end if
      call messages_info(1)
    end do


    call test_prints_info_batch(sys%st, sys%gr, epsib2)

    SAFE_DEALLOCATE_A(dweight)
    SAFE_DEALLOCATE_A(zweight)
    SAFE_DEALLOCATE_A(ddot)
    SAFE_DEALLOCATE_A(zdot)

    call epsib%end()
    call epsib2%end()
    call orbitalbasis_end(basis)
    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_dft_u)
  end subroutine test_dft_u

  ! ---------------------------------------------------------
  subroutine test_hamiltonian(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(wfs_elec_t) :: hpsib
    integer :: itime, terms

    PUSH_SUB(test_hamiltonian)

    !%Variable TestHamiltonianApply
    !%Type integer
    !%Default term_all
    !%Section Calculation Modes::Test
    !%Description
    !% Decides which part of the Hamiltonian is applied.
    !%Option term_all 0
    !% Apply the full Hamiltonian.
    !%Option term_kinetic 1
    !% Apply only the kinetic operator
    !%Option term_local_potential 2
    !% Apply only the local potential.
    !%Option term_non_local_potential 4
    !% Apply only the non_local potential.
    !%End
    call parse_variable(namespace, 'TestHamiltonianApply', OPTION__TESTHAMILTONIANAPPLY__TERM_ALL, terms)
    if (terms == 0) terms = huge(1)


    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the application of the Hamiltonian')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)

    ! Initialize external potential
    if (sys%st%d%pack_states .and. sys%hm%apply_packed()) call sys%st%pack()
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%space, sys%gr, sys%ions, sys%ext_partners, sys%st)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, sys%st, sys%ions, sys%ext_partners)

    call boundaries_set(sys%gr%der%boundaries, sys%gr, sys%st%group%psib(1, 1))

    call sys%st%group%psib(1, 1)%copy_to(hpsib)

    if (sys%hm%apply_packed()) then
      call sys%st%group%psib(1, 1)%do_pack()
      call hpsib%do_pack(copy = .false.)
    end if

    do itime = 1, param%repetitions
      if (states_are_real(sys%st)) then
        call dhamiltonian_elec_apply_batch(sys%hm, sys%namespace, sys%gr, sys%st%group%psib(1, 1), hpsib, terms = terms, &
          set_bc = .false.)
      else
        call zhamiltonian_elec_apply_batch(sys%hm, sys%namespace, sys%gr, sys%st%group%psib(1, 1), hpsib, terms = terms, &
          set_bc = .false.)
      end if
    end do

    if (hpsib%is_packed()) then
      call hpsib%do_unpack(force = .true.)
    end if

    call test_prints_info_batch(sys%st, sys%gr, hpsib)

    call hpsib%end(copy = .false.)
    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_hamiltonian)
  end subroutine test_hamiltonian


  ! ---------------------------------------------------------
  subroutine test_density_calc(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime

    PUSH_SUB(test_density_calc)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing density calculation')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)
    if (sys%st%d%pack_states) call sys%st%pack()

    do itime = 1, param%repetitions
      call density_calc(sys%st, sys%gr, sys%st%rho)
    end do

    write(message(1),'(a,3x, f12.6)') "Norm density  ", dmf_nrm2(sys%gr, sys%st%rho(:,1))
    call messages_info(1, namespace=namespace)

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_density_calc)
  end subroutine test_density_calc


  ! ---------------------------------------------------------
  subroutine test_boundaries(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime

    PUSH_SUB(test_density_calc)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing boundary conditions')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)
    if (sys%st%d%pack_states) call sys%st%pack()

    do itime = 1, param%repetitions
      call boundaries_set(sys%gr%der%boundaries, sys%gr, sys%st%group%psib(1, 1))
    end do

    call test_prints_info_batch(sys%st, sys%gr, sys%st%group%psib(1, 1))

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_density_calc)
  end subroutine test_boundaries


  ! ---------------------------------------------------------
  subroutine test_exponential(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(exponential_t) :: te
    integer :: itime

    PUSH_SUB(test_exponential)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing exponential')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr, wfs_type=TYPE_CMPLX)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)

    ! Initialize external potential
    if (sys%st%d%pack_states .and. sys%hm%apply_packed()) call sys%st%pack()
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%space, sys%gr, sys%ions, sys%ext_partners, sys%st)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, sys%st, sys%ions, sys%ext_partners)

    call exponential_init(te, namespace)

    if (sys%hm%apply_packed()) then
      call sys%st%group%psib(1, 1)%do_pack()
    end if

    do itime = 1, param%repetitions
      call te%apply_batch(sys%namespace, sys%gr, sys%hm, sys%st%group%psib(1, 1), M_ONE)
    end do

    call test_prints_info_batch(sys%st, sys%gr, sys%st%group%psib(1, 1))

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_exponential)
  end subroutine test_exponential


  ! ---------------------------------------------------------
  subroutine test_subspace_diagonalization(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime
    type(subspace_t) :: sdiag
    FLOAT, allocatable :: diff(:)

    PUSH_SUB(test_subspace_diagonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing boundary conditions')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)

    if (sys%st%d%pack_states .and. sys%hm%apply_packed()) call sys%st%pack()
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%space, sys%gr, sys%ions, sys%ext_partners, sys%st)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, sys%st, sys%ions, sys%ext_partners)

    call subspace_init(sdiag, sys%namespace, sys%st)

    SAFE_ALLOCATE(diff(1:sys%st%nst))

    do itime = 1, param%repetitions
      if (states_are_real(sys%st)) then
        call dsubspace_diag(sdiag, sys%namespace, sys%gr, sys%st, sys%hm, 1, sys%st%eigenval(:, 1), diff)
      else
        call zsubspace_diag(sdiag, sys%namespace, sys%gr, sys%st, sys%hm, 1, sys%st%eigenval(:, 1), diff)
      end if
    end do

    SAFE_DEALLOCATE_A(diff)

    call test_prints_info_batch(sys%st, sys%gr, sys%st%group%psib(1, 1))

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_subspace_diagonalization)
  end subroutine test_subspace_diagonalization


  ! ---------------------------------------------------------
  subroutine test_batch_ops(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime, ops, ops_default, ist, jst, nst
    type(wfs_elec_t) :: xx, yy
    FLOAT, allocatable :: tmp(:)
    FLOAT, allocatable :: ddotv(:)
    CMPLX, allocatable :: zdotv(:)
    FLOAT, allocatable :: ddot(:,:)
    CMPLX, allocatable :: zdot(:,:)


    PUSH_SUB(test_batch_ops)

    !%Variable TestBatchOps
    !%Type flag
    !%Default ops_axpy + ops_scal + ops_nrm2
    !%Section Calculation Modes::Test
    !%Description
    !% Decides which part of the Hamiltonian is applied.
    !%Option ops_axpy bit(1)
    !% Tests batch_axpy operation
    !%Option ops_scal bit(2)
    !% Tests batch_scal operation
    !%Option ops_nrm2 bit(3)
    !% Tests batch_nrm2 operation
    !%Option ops_dotp_matrix bit(4)
    !% Tests X(mesh_batch_dotp_matrix)
    !%Option ops_dotp_self bit(5)
    !% Tests X(mesh_batch_dotp_self)
    !%Option ops_dotp_vector bit(6)
    !% Tests X(mesh_batch_dotp_vector)
    !%End
    ops_default = &
      OPTION__TESTBATCHOPS__OPS_AXPY + &
      OPTION__TESTBATCHOPS__OPS_SCAL + &
      OPTION__TESTBATCHOPS__OPS_NRM2 + &
      OPTION__TESTBATCHOPS__OPS_DOTP_MATRIX + &
      OPTION__TESTBATCHOPS__OPS_DOTP_SELF + &
      OPTION__TESTBATCHOPS__OPS_DOTP_VECTOR

    call parse_variable(namespace, 'TestBatchOps', ops_default, ops)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing batch operations')
    call messages_new_line()
    call messages_new_line()
    call messages_info(namespace=namespace)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call states_elec_allocate_wfns(sys%st, sys%gr)
    call test_batch_set_gaussian(sys%st%group%psib(1, 1), sys%gr)
    if (sys%st%d%pack_states) call sys%st%pack()

    if (bitand(ops, OPTION__TESTBATCHOPS__OPS_AXPY) /= 0) then
      message(1) = 'Info: Testing axpy'
      call messages_info(1, namespace=sys%namespace)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      do itime = 1, param%repetitions
        call batch_axpy(sys%gr%np, CNST(0.1), xx, yy)
      end do
      call test_prints_info_batch(sys%st, sys%gr, yy, string = "axpy")

      call xx%end()
      call yy%end()
    end if

    if (bitand(ops, OPTION__TESTBATCHOPS__OPS_SCAL) /= 0) then
      message(1) = 'Info: Testing scal'
      call messages_info(1, namespace=sys%namespace)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      do itime = 1, param%repetitions
        call batch_scal(sys%gr%np, CNST(0.1), yy)
      end do
      call test_prints_info_batch(sys%st, sys%gr, yy, string="scal")

      call xx%end()
      call yy%end()
    end if

    if (bitand(ops, OPTION__TESTBATCHOPS__OPS_NRM2) /= 0) then
      message(1) = 'Info: Testing nrm2'
      call messages_info(1, namespace=sys%namespace)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      SAFE_ALLOCATE(tmp(1:xx%nst))

      do itime = 1, param%repetitions
        call mesh_batch_nrm2(sys%gr, yy, tmp)
      end do
      do itime = 1, xx%nst
        write(message(1),'(a,i1,3x,e13.6)') "Nrm2 norm state  ", itime, tmp(itime)
        call messages_info(1, namespace=sys%namespace)
      end do

      SAFE_DEALLOCATE_A(tmp)

      call xx%end()
      call yy%end()
    end if

    if (bitand(ops, OPTION__TESTBATCHOPS__OPS_DOTP_MATRIX) /= 0) then

      message(1) = 'Info: Testing dotp_matrix'
      call messages_info(1, namespace=sys%namespace)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      nst = sys%st%group%psib(1, 1)%nst

      if (states_are_real(sys%st)) then
        SAFE_ALLOCATE(ddot(1:nst, 1:nst))
        call dmesh_batch_dotp_matrix(sys%gr, xx, yy, ddot)

        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,e13.6)') 'Dotp_matrix states', ist, jst, ddot(ist,jst)
          end do
          call messages_info(nst, namespace=sys%namespace)
        end do
        SAFE_DEALLOCATE_A(ddot)
      else
        SAFE_ALLOCATE(zdot(1:nst, 1:nst))
        call zmesh_batch_dotp_matrix(sys%gr, xx, yy, zdot)

        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,2e14.6)') 'Dotp_matrix states', ist, jst, zdot(ist,jst)
          end do
          call messages_info(nst, namespace=sys%namespace)
        end do
        SAFE_DEALLOCATE_A(zdot)
      end if

      call xx%end()
      call yy%end()
    end if

    if (bitand(ops, OPTION__TESTBATCHOPS__OPS_DOTP_VECTOR) /= 0) then

      message(1) = 'Info: Testing dotp_vector'
      call messages_info(1, namespace=sys%namespace)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      nst = sys%st%group%psib(1, 1)%nst

      if (states_are_real(sys%st)) then
        SAFE_ALLOCATE(ddotv(1:nst))
        call dmesh_batch_dotp_vector(sys%gr, xx, yy, ddotv)

        do ist = 1, nst
          write(message(ist), '(a,i3,3x,e13.6)') 'Dotp_vector state', ist, ddotv(ist)
        end do
        call messages_info(nst, namespace=sys%namespace)
        SAFE_DEALLOCATE_A(ddotv)
      else
        SAFE_ALLOCATE(zdotv(1:nst))
        call zmesh_batch_dotp_vector(sys%gr, xx, yy, zdotv)
        do ist = 1, nst
          write(message(ist), '(a,i3,3x,2e14.6)') 'Dotp_vector state', ist, zdotv(ist)
        end do
        call messages_info(nst, namespace=sys%namespace)
        SAFE_DEALLOCATE_A(zdotv)
      end if

      call xx%end()
      call yy%end()
    end if

    if (bitand(ops, OPTION__TESTBATCHOPS__OPS_DOTP_SELF) /= 0) then

      message(1) = 'Info: Testing dotp_self'
      call messages_info(1, namespace=sys%namespace)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)

      nst = sys%st%group%psib(1, 1)%nst

      if (states_are_real(sys%st)) then
        SAFE_ALLOCATE(ddot(1:nst, 1:nst))
        call dmesh_batch_dotp_self(sys%gr, xx, ddot)

        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,e13.6)') 'Dotp_self states', ist, jst, ddot(ist,jst)
          end do
          call messages_info(nst, namespace=sys%namespace)
        end do
        SAFE_DEALLOCATE_A(ddot)
      else
        SAFE_ALLOCATE(zdot(1:nst, 1:nst))
        call zmesh_batch_dotp_self(sys%gr, xx, zdot)
        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,2e14.6)') 'Dotp_self states', ist, jst, zdot(ist,jst)
          end do
          call messages_info(nst, namespace=sys%namespace)
        end do
        SAFE_DEALLOCATE_A(zdot)
      end if

      call xx%end()
    end if

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_batch_ops)
  end subroutine test_batch_ops

! ---------------------------------------------------------
  subroutine test_derivatives(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_derivatives)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    message(1) = 'Info: Testing the finite-differences derivatives.'
    message(2) = ''
    call messages_info(2, namespace=namespace)

    if (param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      call dderivatives_test(sys%gr%der, sys%namespace, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    if (param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      call zderivatives_test(sys%gr%der, sys%namespace, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_derivatives)
  end subroutine test_derivatives

  ! ---------------------------------------------------------

  subroutine test_orthogonalization(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime

    PUSH_SUB(test_orthogonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
    call calc_mode_par_set_scalapack_compat()

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    message(1) = 'Info: Testing orthogonalization.'
    message(2) = ''
    call messages_info(2, namespace=namespace)

    if (param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      message(1) = 'Info: Real wave-functions.'
      call messages_info(1, namespace=namespace)
      do itime = 1, param%repetitions
        call dstates_elec_calc_orth_test(sys%st, sys%namespace, sys%gr, sys%kpoints)
      end do
    end if

    if (param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      message(1) = 'Info: Complex wave-functions.'
      call messages_info(1, namespace=namespace)
      do itime = 1, param%repetitions
        call zstates_elec_calc_orth_test(sys%st, sys%namespace, sys%gr, sys%kpoints)
      end do
    end if

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_orthogonalization)
  end subroutine test_orthogonalization

  ! ---------------------------------------------------------

  subroutine test_interpolation(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_interpolation)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    if (param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      call messages_write('Info: Testing real interpolation routines')
      call messages_new_line()
      call messages_new_line()
      call messages_info(namespace=namespace)

      call dmesh_interpolation_test(sys%gr)
    end if

    if (param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      call messages_new_line()
      call messages_write('Info: Testing complex interpolation routines')
      call messages_new_line()
      call messages_new_line()
      call messages_info(namespace=namespace)

      call zmesh_interpolation_test(sys%gr)
    end if

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_interpolation)
  end subroutine test_interpolation


  ! ---------------------------------------------------------

  subroutine test_ion_interaction(namespace)
    type(namespace_t),        intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_ion_interaction)

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call ion_interaction_test(sys%space, sys%ions%latt, sys%ions%atom, sys%ions%natoms, sys%ions%pos, &
      sys%gr%box%bounding_box_l, namespace, sys%mc)

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_ion_interaction)
  end subroutine test_ion_interaction

  ! ---------------------------------------------------------

  subroutine test_prints_info_batch(st, gr, psib, string)
    type(states_elec_t), intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    class(batch_t),      intent(inout) :: psib
    character(*), optional,  intent(in)    :: string

    integer :: itime
    CMPLX, allocatable :: zpsi(:, :)
    FLOAT, allocatable :: dpsi(:, :)

    character(80)      :: string_

    string_ = optional_default(string, "")

    PUSH_SUB(test_prints_info_batch)

    if (states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(1:gr%np, 1:st%d%dim))
    else
      SAFE_ALLOCATE(zpsi(1:gr%np, 1:st%d%dim))
    end if

    do itime = 1, psib%nst
      if (states_are_real(st)) then
        call batch_get_state(psib, itime, gr%np, dpsi)
        write(message(1),'(a,i2,3x,e13.6)') "Norm state "//trim(string_)//" ", itime, dmf_nrm2(gr, st%d%dim, dpsi)
      else
        call batch_get_state(psib, itime, gr%np, zpsi)
        write(message(1),'(a,i2,3x,e13.6)') "Norm state "//trim(string_)//" ", itime, zmf_nrm2(gr, st%d%dim, zpsi)
      end if
      call messages_info(1)
    end do

    if (states_are_real(st)) then
      SAFE_DEALLOCATE_A(dpsi)
    else
      SAFE_DEALLOCATE_A(zpsi)
    end if

    POP_SUB(test_prints_info_batch)

  end subroutine test_prints_info_batch


  ! ---------------------------------------------------------
  subroutine test_clock()

    type(clock_t) :: test_clock_a, test_clock_b

    PUSH_SUB(test_clock)

    test_clock_a = clock_t(time_step=M_TWO, initial_tick=100)
    test_clock_b = clock_t(time_step=M_ONE)
    call test_clock_a%print()
    call test_clock_b%print()

    call test_clock_a%set_time(test_clock_b)
    call test_clock_a%print()
    test_clock_a = test_clock_a + CLOCK_TICK
    call test_clock_a%print()
    test_clock_a = test_clock_a - CLOCK_TICK
    call test_clock_a%print()
    test_clock_a = test_clock_a + CLOCK_TICK
    call test_clock_a%print()
    call test_clock_a%reset()
    call test_clock_a%print()
    test_clock_a = test_clock_a + 3*CLOCK_TICK
    call test_clock_a%print()
    test_clock_a = test_clock_a - 2*CLOCK_TICK
    call test_clock_a%print()
    message(1) = test_clock_a%print_str()
    call messages_info(1)

    write(message(1),'(A,x,I10.10)') &
      'clock_get_tick', test_clock_a%get_tick()
    write(message(2),'(A,x,F15.10)') &
      'clock_time', test_clock_a%time()
    write(message(3),'(A,x,I1)')     &
      'clock_is_earlier', abs(transfer(test_clock_a .lt. test_clock_b, 0))
    write(message(4),'(A,x,I1)')     &
      'clock_is_equal_or_earlier', abs(transfer(test_clock_a .le. test_clock_b, 0))
    write(message(5),'(A,x,I1)')     &
      'clock_is_later', abs(transfer(test_clock_a .gt. test_clock_b, 0))
    write(message(6),'(A,x,I1)')     &
      'clock_is_equal_or_later', abs(transfer(test_clock_a .ge. test_clock_b, 0))
    write(message(7),'(A,x,I1)')     &
      'clock_is_equal', abs(transfer(test_clock_a .eq. test_clock_b, 0))
    call messages_info(7)

    POP_SUB(test_clock)
  end subroutine test_clock


  ! ---------------------------------------------------------
  subroutine test_cgal()

    type(cgal_polyhedra_t) :: cgal_poly

    PUSH_SUB(test_cgal)

    call cgal_polyhedron_init(cgal_poly, "28-cgal.02-X.off", verbose = .true.)

    if (cgal_polyhedron_point_inside(cgal_poly, CNST(30.), CNST(10.), CNST(30.))) then
      message(1) = "cgal_polyhedron_point_inside"
      call messages_info(1)
    end if

    call cgal_polyhedron_end(cgal_poly)

    POP_SUB(test_cgal)
  end subroutine test_cgal


  ! ---------------------------------------------------------
  subroutine test_dense_eigensolver()
    integer :: N, ii, jj, N_list(4), i_N
    FLOAT, allocatable :: matrix(:, :), eigenvectors(:, :), eigenvalues(:), test(:)
    FLOAT, allocatable :: differences(:)

    PUSH_SUB(test_dense_eigensolver)

    N_list = [15, 32, 100, 500]

    do i_N = 1, 4
      N = N_list(i_N)
      SAFE_ALLOCATE(matrix(1:N, 1:N))
      SAFE_ALLOCATE(eigenvectors(1:N, 1:N))
      SAFE_ALLOCATE(eigenvalues(1:N))
      SAFE_ALLOCATE(test(1:N))
      SAFE_ALLOCATE(differences(1:N))


      ! set up symmetrix matrix
      do jj = 1, N
        do ii = 1, N
          matrix(ii, jj) = ii * jj
        end do
      end do

      ! parallel solver
      eigenvectors(1:N, 1:N) = matrix(1:N, 1:N)
      call lalg_eigensolve_parallel(N, eigenvectors, eigenvalues)

      do ii = 1, N
        test(:) = matmul(matrix, eigenvectors(:, ii)) - eigenvalues(ii) * eigenvectors(:, ii)
        differences(ii) = sum(abs(test)) / sum(abs(eigenvectors(:, ii)))
      end do
      write(message(1), "(A, I3, A, E13.6)") "Parallel solver - N: ", N, &
        ", average difference: ", sum(differences)/N
      call messages_info(1)

      ! serial solver
      eigenvectors(1:N, 1:N) = matrix(1:N, 1:N)
      call lalg_eigensolve(N, eigenvectors, eigenvalues)

      do ii = 1, N
        test(:) = matmul(matrix, eigenvectors(:, ii)) - eigenvalues(ii) * eigenvectors(:, ii)
        differences(ii) = sum(abs(test)) / sum(abs(eigenvectors(:, ii)))
      end do
      write(message(1), "(A, I3, A, E13.6)") "Serial solver   - N: ", N, &
        ", average difference: ", sum(differences)/N
      call messages_info(1)

      SAFE_DEALLOCATE_A(matrix)
      SAFE_DEALLOCATE_A(eigenvectors)
      SAFE_DEALLOCATE_A(eigenvalues)
      SAFE_DEALLOCATE_A(test)
      SAFE_DEALLOCATE_A(differences)
    end do

    POP_SUB(test_dense_eigensolver)
  end subroutine test_dense_eigensolver

  subroutine test_batch_set_gaussian(psib, mesh)
    class(batch_t), intent(inout) :: psib
    class(mesh_t),  intent(in)    :: mesh

    FLOAT, allocatable :: dff(:)
    CMPLX, allocatable :: zff(:)
    integer :: ist, ip
    FLOAT :: da, db, dc
    CMPLX :: za, zb, zc

    PUSH_SUB(test_batch_set_gaussian)

    ! use a similar function as in the derivatives test
    da = M_ONE/mesh%box%bounding_box_l(1)
    db = CNST(10.0)
    dc = CNST(100.0)

    if (type_is_complex(psib%type())) then
      ! we make things more "complex"
      za = da + M_ZI*CNST(0.01)
      zb = db*exp(M_ZI*CNST(0.345))
      zc = dc - M_ZI*CNST(50.0)

      SAFE_ALLOCATE(zff(1:mesh%np))
      do ist = 1, psib%nst_linear
        za = za * ist
        zb = zb / ist
        do ip = 1, mesh%np
          zff(ip) = zb*exp(-za*sum(mesh%x(ip, :)**2)) + zc
        end do
        call batch_set_state(psib, ist, mesh%np, zff)
      end do
      call zmesh_batch_normalize(mesh, psib)
      SAFE_DEALLOCATE_A(zff)
    else
      SAFE_ALLOCATE(dff(1:mesh%np))
      do ist = 1, psib%nst_linear
        da = da * ist
        db = db / ist
        do ip = 1, mesh%np
          dff(ip) = db*exp(-da*sum(mesh%x(ip, :)**2)) + dc
        end do
        call batch_set_state(psib, ist, mesh%np, dff)
      end do
      call dmesh_batch_normalize(mesh, psib)
      SAFE_DEALLOCATE_A(dff)
    end if

    POP_SUB(test_batch_set_gaussian)
  end subroutine test_batch_set_gaussian

  ! ---------------------------------------------------------
  subroutine test_grid_interpolation()
    type(electrons_t), pointer :: sys
    type(multigrid_t) :: mgrid

    PUSH_SUB(test_grid_interpolation)

    sys => electrons_t(global_namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    call multigrid_init(mgrid, global_namespace, sys%space, sys%gr, sys%gr%der, &
      sys%gr%stencil, sys%mc, used_for_preconditioner = .true.)

    call multigrid_test_interpolation(mgrid, sys%space)

    call multigrid_end(mgrid)

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_grid_interpolation)
  end subroutine test_grid_interpolation

  subroutine test_iihash()

    type(iihash_t) :: h
    integer :: value
    logical :: found

    call iihash_init(h)
    call iihash_insert(h, 1,5)
    call iihash_insert(h, 3,7)

    value = iihash_lookup(h, 1, found)
    write(message(1),*) 'hash[1] :', found, value
    call messages_info()

    value = iihash_lookup(h, 2, found)
    write(message(1),*) 'hash[2] :', found, value
    call messages_info()

    value = iihash_lookup(h, 3, found)
    write(message(1),*) 'hash[3] :', found, value
    call messages_info()

    call iihash_end(h)

  end subroutine test_iihash

  subroutine test_sihash()
    type(sihash_t)          :: h
    type(sihash_iterator_t) :: it

    integer :: counter
    integer :: value, sum
    logical :: found

    call sihash_init(h)
    call sihash_insert(h, "one",5)
    call sihash_insert(h, "three",7)
    call sihash_insert(h, "the answer", 42)

    value = sihash_lookup(h, "one", found)
    write(message(1),*) 'hash["one"]:   ', found, value
    call messages_info()

    value = sihash_lookup(h, "two", found)
    write(message(1),*) 'hash["two"]:   ', found, value
    call messages_info()

    value = sihash_lookup(h, "three", found)
    write(message(1),*) 'hash["three"]: ', found, value
    call messages_info()

    sum = 0
    counter = 1
    call it%start(h)

    do while (it%has_next())
      value = it%get_next()
      sum = sum + value
      write(message(1),'(I3,A,I5)') counter,': hash[...] = ',value
      call messages_info()
      counter = counter + 1
    end do
    write(message(1),*) 'counter = ', counter
    write(message(2),*) 'sum = ', sum
    call messages_info(2)


    call sihash_end(h)

  end subroutine test_sihash


  subroutine test_sphash(namespace)
    type(namespace_t), intent(in)  :: namespace

    type(sphash_t)          :: h
    type(sphash_iterator_t) :: it

    logical :: found


    type(clock_t)              :: clock_1
    type(clock_t), allocatable :: clock_2
    type(space_t) :: space_1

    class(*), pointer :: value

    integer :: count_clock, count_space

    SAFE_ALLOCATE(clock_2)

    clock_1 = clock_t(1.d-5, 1)
    clock_2 = clock_t(2.d-5, 2)
    call space_init(space_1, namespace)

    call sphash_init(h)
    call sphash_insert(h, "one",   clock_1)
    call sphash_insert(h, "two",   space_1)
    call sphash_insert(h, "three", clock_2, clone=.true.)

    value => sphash_lookup(h, "one", found)
    select type(value)
    type is (clock_t)
      write(message(1),*) 'hash["one"]:   ', found, value%get_tick()
      call messages_info()
    type is (space_t)
      write(message(1),*) 'hash["one"]:   ', found, value%short_info()
      call messages_info()
    class default
      write(message(1),*) 'wrong type. found = ', found
      call messages_info()
    end select

    value => sphash_lookup(h, "two", found)
    select type(value)
    type is (clock_t)
      write(message(1),*) 'hash["two"]:   ', found, value%get_tick()
      call messages_info()
    type is (space_t)
      write(message(1),*) 'hash["two"]:   ', found, value%short_info()
      call messages_info()
    class default
      write(message(1),*) 'wrong type. found = ',found
      call messages_info()
    end select

    SAFE_DEALLOCATE_A(clock_2)

    value => sphash_lookup(h, "three", found)
    select type(value)
    type is (clock_t)
      write(message(1),*) 'hash["three"]:   ', found, value%get_tick()
      call messages_info()
    type is (space_t)
      write(message(1),*) 'hash["three"]:   ', found, value%short_info()
      call messages_info()
    class default
      write(message(1),*) 'wrong type. found = ',found
      call messages_info()
    end select

    count_clock = 0
    count_space = 0

    call it%start(h)

    do while (it%has_next())
      value => it%get_next()
      select type(value)
      type is (clock_t)
        count_clock = count_clock + 1
      type is (space_t)
        count_space = count_space + 1
      end select
    end do

    write(message(1), *) 'Count_clock = ', count_clock
    write(message(2), *) 'Count_space = ', count_space
    call messages_info(2)

    call sphash_end(h)

  end subroutine test_sphash

  ! ---------------------------------------------------------
  subroutine test_regridding(namespace)
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sysA, sysB
    type(regridding_t), pointer :: regridding
    FLOAT, allocatable :: ff_A(:), ff_A_reference(:), ff_B(:), ff_B_reference(:), diff_A(:), diff_B(:)
    FLOAT :: norm_ff, norm_diff
    integer :: ip, ierr

    PUSH_SUB(test_regridding)

    sysA => electrons_t(namespace_t("A", namespace), generate_epot=.false.)
    sysB => electrons_t(namespace_t("B", namespace), generate_epot=.false.)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default=.true.)
    call sysA%init_parallelization(mpi_world)
    call sysB%init_parallelization(mpi_world)


    SAFE_ALLOCATE(ff_A(1:sysA%gr%np))
    SAFE_ALLOCATE(ff_A_reference(1:sysA%gr%np))
    SAFE_ALLOCATE(diff_A(1:sysA%gr%np))
    SAFE_ALLOCATE(ff_B(1:sysB%gr%np))
    SAFE_ALLOCATE(ff_B_reference(1:sysB%gr%np))
    SAFE_ALLOCATE(diff_B(1:sysB%gr%np))

    do ip = 1, sysA%gr%np
      ff_A_reference(ip) = values(sysA%gr%x(ip, :))
    end do
    do ip = 1, sysB%gr%np
      ff_B_reference(ip) = values(sysB%gr%x(ip, :))
    end do

    ! forward mapping A->B
    regridding => regridding_t(sysB%gr, sysA%gr, sysA%space, sysA%namespace)
    call regridding%do_transfer(ff_B, ff_A_reference)
    SAFE_DEALLOCATE_P(regridding)

    ! check that mapped function is exactly the same for the points that are available on both grids
    do ip = 1, sysB%gr%np
      if (ff_B(ip) == M_ZERO) then
        ff_B_reference(ip) = M_ZERO
        diff_B(ip) = M_ZERO
      else
        diff_B(ip) = abs(ff_B_reference(ip) - ff_B(ip))
      end if
    end do
    norm_ff = dmf_nrm2(sysB%gr, ff_B_reference)
    norm_diff = dmf_nrm2(sysB%gr, diff_B)

    write(message(1),'(a, E14.6)') "Forward: difference of reference to mapped function (rel.): ", &
      norm_diff/norm_ff
    call messages_info(1, namespace=namespace)

    call dio_function_output(io_function_fill_how('AxisX'), ".", "forward_reference", namespace, sysB%space, &
      sysB%gr, ff_B_reference, unit_one, ierr)
    call dio_function_output(io_function_fill_how('AxisX'), ".", "forward_mapped", namespace, sysB%space, &
      sysB%gr, ff_B, unit_one, ierr)
    call dio_function_output(io_function_fill_how('AxisX'), ".", "forward_original", namespace, sysA%space, &
      sysA%gr, ff_A_reference, unit_one, ierr)

    ! backward mapping B->A
    regridding => regridding_t(sysA%gr, sysB%gr, sysB%space, sysB%namespace)
    call regridding%do_transfer(ff_A, ff_B_reference)
    SAFE_DEALLOCATE_P(regridding)

    do ip = 1, sysA%gr%np
      if (ff_A(ip) == M_ZERO) then
        ff_A_reference(ip) = M_ZERO
        diff_A(ip) = M_ZERO
      else
        diff_A(ip) = abs(ff_A_reference(ip) - ff_A(ip))
      end if
    end do
    norm_ff = dmf_nrm2(sysA%gr, ff_A_reference)
    norm_diff = dmf_nrm2(sysA%gr, diff_A)

    write(message(1),'(a, E14.6)') "Backward: difference of reference to mapped function (rel.): ", &
      norm_diff/norm_ff
    call messages_info(1, namespace=namespace)

    call dio_function_output(io_function_fill_how('AxisX'), ".", "backward_reference", namespace, sysA%space, &
      sysA%gr, ff_A_reference, unit_one, ierr)
    call dio_function_output(io_function_fill_how('AxisX'), ".", "backward_mapped", namespace, sysA%space, &
      sysA%gr, ff_A, unit_one, ierr)
    call dio_function_output(io_function_fill_how('AxisX'), ".", "backward_original", namespace, sysB%space, &
      sysB%gr, ff_B_reference, unit_one, ierr)

    SAFE_DEALLOCATE_A(ff_A)
    SAFE_DEALLOCATE_A(ff_A_reference)
    SAFE_DEALLOCATE_A(ff_B)
    SAFE_DEALLOCATE_A(ff_B_reference)
    SAFE_DEALLOCATE_A(diff_A)
    SAFE_DEALLOCATE_A(diff_B)
    SAFE_DEALLOCATE_P(sysA)
    SAFE_DEALLOCATE_P(sysB)

    call messages_info(1, namespace=namespace)

    POP_SUB(test_regridding)
  contains
    FLOAT function values(xx)
      FLOAT, intent(in) :: xx(:)
      FLOAT :: xx0(1:size(xx, dim=1))
      FLOAT, parameter :: aa = M_HALF
      FLOAT, parameter :: bb = M_FOUR

      ! no push_sub/pop_sub because this function is called often
      xx0(:) = M_ONE
      values = bb * exp(-aa*sum((xx-xx0)**2))
    end function values
  end subroutine test_regridding

  !> Here, analytical formulation for vector potential and B field are used.
  !! Ref: Sangita Sen and Erik I. Tellgren, <i>J. Chem. Theory Comput.</i> <i>17</i>, <b>3</b> (2021).
  !! Analytical input for vector potential /f$ A_{r}= 1/3[ -xz, yz, x^2 - y^2] \f$
  !! When bounded, above expression is multiplied with gaussian envelope
  !! \f$ (1/box_size)*exp^(-x^2-y^2-z^2)] \f$
  subroutine test_vecpot_analytical(namespace)
    type(namespace_t), intent(in)  :: namespace

    class(maxwell_t), pointer :: maxwell_system

    FLOAT, allocatable :: magnetic_field(:,:)
    FLOAT, allocatable :: vector_potential_mag(:,:)
    FLOAT, allocatable :: vector_potential_analytical(:,:)
    FLOAT, allocatable :: delta(:,:)
    FLOAT              :: exp_factor
    FLOAT              :: xx
    FLOAT              :: yy
    FLOAT              :: zz
    FLOAT              :: sigma
    integer :: ip, j, ierr, nn
    integer(i8) :: out_how
    character(len=MAX_PATH_LEN) :: fname, fname2, fname3

    out_how = 32
    maxwell_system => maxwell_t(namespace)
    sigma = maxwell_system%gr%box%bounding_box_l(1)/CNST(10)  ! this is exponential width
    call maxwell_system%init_parallelization(mpi_world)

    SAFE_ALLOCATE(magnetic_field(1:maxwell_system%gr%np_part, 1:3))
    SAFE_ALLOCATE(vector_potential_mag(1:maxwell_system%gr%np_part, 1:3))
    SAFE_ALLOCATE(vector_potential_analytical(1:maxwell_system%gr%np_part, 1:3))
    SAFE_ALLOCATE(delta(1:maxwell_system%gr%np, 1:3))

    !%Variable TestVectorPotentialType
    !%Type integer
    !%Default bounded
    !%Section Calculation Modes::Test
    !%Description
    !% Select whether bounded or unbounded type will be used for vector potential tests
    !%Option bounded 1
    !% Analytical Vector Potential formulation is bounded by spatial gaussian
    !%Option unbounded 2
    !% Analytical Vector Potential is not bounded
    !%End
    call parse_variable(namespace, 'TestVectorPotentialType', OPTION__TESTVECTORPOTENTIALTYPE__BOUNDED, nn)

    select case (nn)
    case (OPTION__TESTVECTORPOTENTIALTYPE__BOUNDED)  ! bounded input
      do ip = 1, maxwell_system%gr%np_part
        xx = maxwell_system%gr%x(ip, 1)
        yy = maxwell_system%gr%x(ip, 2)
        zz = maxwell_system%gr%x(ip, 3)
        exp_factor = exp((-xx**2 - yy**2 - zz**2)*1/(2*sigma**2))
        magnetic_field(ip, 1) = exp_factor*yy*(1 - (-xx**2 + yy**2)/(3*sigma**2) - zz**2/(3*sigma**2))
        magnetic_field(ip, 2) = exp_factor * xx * (1 + (-xx**2 + yy**2)/(3*sigma**2) - zz**2/(3*sigma**2))
        magnetic_field(ip, 3) = exp_factor * 2 * xx * yy * zz * 1/(3*sigma**2)

        vector_potential_analytical(ip, 1) = M_THIRD * xx * zz * exp_factor
        vector_potential_analytical(ip, 2) = - M_THIRD * yy * zz * exp_factor
        vector_potential_analytical(ip, 3) = M_THIRD * (-xx**2 + yy**2) * exp_factor
      end do
    case (OPTION__TESTVECTORPOTENTIALTYPE__UNBOUNDED)  ! unbounded input, TODO this unit test requires implementation of BCs for Helmholtz decomposition
      do ip = 1, maxwell_system%gr%np_part
        magnetic_field(ip, 1) = maxwell_system%gr%x(ip, 2)
        magnetic_field(ip, 2) = maxwell_system%gr%x(ip, 1)
        magnetic_field(ip, 3) = M_ZERO

        vector_potential_analytical(ip, 1) =   M_THIRD * maxwell_system%gr%x(ip, 1) * maxwell_system%gr%x(ip, 3)
        vector_potential_analytical(ip, 2) = - M_THIRD * maxwell_system%gr%x(ip, 2) * maxwell_system%gr%x(ip, 3)
        vector_potential_analytical(ip, 3) = - M_THIRD * (maxwell_system%gr%x(ip, 1)**2 - maxwell_system%gr%x(ip, 2)**2)
      end do
    end select
    call maxwell_system%helmholtz%get_vector_potential(namespace, vector_potential_mag, magnetic_field)

    do j = 1, 3
      delta(:,:) = M_ZERO
      do ip = 1, maxwell_system%gr%np
        delta(ip,j) = vector_potential_analytical(ip, j) - vector_potential_mag(ip, j)
      end do
    end do

    do j = 1, 3
      write(message(j),*) 'j, norm2(delta)', j, norm2(delta(:,j))
    end do
    call messages_info(3)

    write(fname, '(a)') 'deviation_from_analytical_formulation' ! use messages later
    call io_function_output_vector(out_how , './' , trim(fname), namespace, maxwell_system%space, maxwell_system%gr, &
      delta, unit_one, ierr)
    write(fname2, '(a)') 'vector_potential_analytical'
    call io_function_output_vector(out_how , './' , trim(fname2), namespace, maxwell_system%space, maxwell_system%gr, &
      vector_potential_analytical, unit_one, ierr)
    write(fname3, '(a)') 'vector_potential_mag'
    call io_function_output_vector(out_how , './' , trim(fname3), namespace, maxwell_system%space, maxwell_system%gr, &
      vector_potential_mag, unit_one, ierr)

    SAFE_DEALLOCATE_A(magnetic_field)
    SAFE_DEALLOCATE_A(vector_potential_mag)
    SAFE_DEALLOCATE_A(vector_potential_analytical)
    SAFE_DEALLOCATE_A(delta)
    SAFE_DEALLOCATE_P(maxwell_system)

  end subroutine test_vecpot_analytical

  ! ---------------------------------------------------------
  subroutine multigrid_test_interpolation(mgrid, space)
    type(multigrid_t), intent(in) :: mgrid
    type(space_t),     intent(in) :: space

    FLOAT, allocatable :: guess0(:), res0(:), guess1(:)
    type(mesh_t), pointer :: mesh0, mesh1
    FLOAT :: delta, xx(3,2), alpha, beta, rr
    integer :: nn, ip, ierr, order

    PUSH_SUB(multigrid_test_interpolation)

    !%Variable InterpolationTestOrder
    !%Type integer
    !%Default 1
    !%Section Calculation Modes::Test
    !%Description
    !% This variable controls the order of the grid interpolation
    !% used in the corresponding unit test.
    !%End
    call parse_variable(global_namespace, 'InterpolationTestOrder', 1, order)

    message(1) = 'Info: Testing the grid interpolation.'
    message(2) = ''
    call messages_info(2)

    mesh0 => mgrid%level(0)%mesh
    mesh1 => mgrid%level(1)%mesh

    SAFE_ALLOCATE(guess0(1:mesh0%np_part))
    SAFE_ALLOCATE(res0(1:mesh0%np))
    SAFE_ALLOCATE(guess1(1:mesh1%np_part))

    alpha = M_FOUR*mesh0%spacing(1)
    beta = M_ONE / (alpha**space%dim * sqrt(M_PI)**space%dim)

    ! Set the centers of the Gaussians by hand
    xx(1, 1) = M_ONE
    xx(2, 1) = -M_HALF
    xx(3, 1) = M_TWO
    xx(1, 2) = -M_TWO
    xx(2, 2) = M_ZERO
    xx(3, 2) = -M_ONE
    xx = xx * alpha

    ! Density as sum of Gaussians
    guess0 = M_ZERO
    do nn = 1, 2
      do ip = 1, mesh0%np
        call mesh_r(mesh0, ip, rr, origin = xx(:, nn))
        guess0(ip) = guess0(ip) + (-1)**nn * beta*exp(-(rr/alpha)**2)
      end do
    end do

    call dio_function_output (io_function_fill_how('AxisX'), ".", "interpolation_target", global_namespace, &
      space, mesh0, guess0, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisZ'), ".", "interpolation_target", global_namespace, &
      space, mesh0, guess0, unit_one, ierr)
    call dio_function_output (io_function_fill_how('PlaneZ'), ".", "interpolation_target", global_namespace, &
      space, mesh0, guess0, unit_one, ierr)

    ! We start by testing the interpolation scheme. For this, we generate a function on the fine grid
    ! and we inject it on the coarse grid. Then we interpolate it back to the fine grid and we compare
    ! This allows for testing the quality of the interpolation scheme

    ! move to level  1
    call dmultigrid_fine2coarse(mgrid%level(1)%tt, mgrid%level(0)%der, mesh1, guess0, guess1, INJECTION)
    ! back to level 0
    call dmultigrid_coarse2fine(mgrid%level(1)%tt, mgrid%level(1)%der, mesh0, guess1, res0, order = order)

    call dio_function_output (io_function_fill_how('AxisX'), ".", "interpolation_result", global_namespace, &
      space, mesh0, res0, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisZ'), ".", "interpolation_result", global_namespace, &
      space, mesh0, res0, unit_one, ierr)
    call dio_function_output (io_function_fill_how('PlaneZ'), ".", "interpolation_result", global_namespace, &
      space, mesh0, res0, unit_one, ierr)

    delta = dmf_nrm2(mesh0, guess0(1:mesh0%np)-res0)
    write(message(1),'(a,e13.6)') 'Interpolation test (abs.) = ', delta

    ! Now we test if the restriction+interpolation combination returns the original result or not
    ! This allows to test if restriction and interpolation are adjoint operators or not.

    ! move to level  1
    call dmultigrid_fine2coarse(mgrid%level(1)%tt, mgrid%level(0)%der, mesh1, guess0, guess1, FULLWEIGHT)
    ! back to level 0
    call dmultigrid_coarse2fine(mgrid%level(1)%tt, mgrid%level(1)%der, mesh0, guess1, res0, order = order)

    call dio_function_output (io_function_fill_how('AxisX'), ".", "restriction_result", global_namespace, &
      space, mesh0, res0, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisZ'), ".", "restriction_result", global_namespace, &
      space, mesh0, res0, unit_one, ierr)
    call dio_function_output (io_function_fill_how('PlaneZ'), ".", "restriction_result", global_namespace, &
      space, mesh0, res0, unit_one, ierr)

    delta = dmf_nrm2(mesh0, guess0(1:mesh0%np)-res0)
    write(message(2),'(a,e13.6)')  'Restriction test (abs.) = ', delta
    call messages_info(2)

    SAFE_DEALLOCATE_A(guess0)
    SAFE_DEALLOCATE_A(res0)
    SAFE_DEALLOCATE_A(guess1)

    POP_SUB(multigrid_test_interpolation)
  end subroutine multigrid_test_interpolation


  !> Here we test the different contributions to the total electronic current density
  subroutine test_current_density(namespace)
    type(namespace_t), intent(in) :: namespace
    type(electrons_t), pointer :: sys

    type(current_t) :: current
    character(len=MAX_PATH_LEN) :: fname
    integer :: ierr, ip, idir
    integer(i8) :: out_how
    FLOAT, allocatable :: current_para_ref(:,:), current_dia_ref(:,:), current_mag_ref(:,:), delta(:)
    FLOAT :: xx(3), rr, a0, mag_curr, sin_thet, sin_phi, cos_phi, vec_pot_slope
    CMPLX :: alpha

    sys => electrons_t(namespace, generate_epot=.false.)
    call sys%init_parallelization(mpi_world)

    alpha = (CNST(0.0), CNST(0.5))
    a0 = M_ONE
    vec_pot_slope = CNST(0.4) ! units of B field

    call states_elec_allocate_wfns(sys%st, sys%gr, wfs_type = TYPE_CMPLX)
    call set_hydrogen_states(sys%st%group%psib(1, 1), sys%gr, namespace, alpha, a0)

    call current_init(current, namespace)

    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%space, sys%gr, sys%ions, sys%ext_partners, sys%st)

    SAFE_ALLOCATE(sys%hm%hm_base%vector_potential(1:3, 1:sys%gr%np))

    sys%hm%hm_base%vector_potential = M_ZERO
    do ip = 1, sys%gr%np
      xx = sys%gr%x(ip, 1:3)
      sys%hm%hm_base%vector_potential(2, ip) = vec_pot_slope * xx(1) / P_c  ! vector potential is here devided by c_0
    end do

    call states_elec_allocate_current(sys%st, sys%space, sys%gr)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    ! paramagnetic + diamagnetic current densities
    call current_calculate(current, namespace, sys%gr, sys%hm, sys%space, sys%st)

    SAFE_ALLOCATE(current_para_ref(1:sys%gr%np,1:3))
    SAFE_ALLOCATE(current_dia_ref(1:sys%gr%np,1:3))
    SAFE_ALLOCATE(current_mag_ref(1:sys%gr%np,1:3))
    SAFE_ALLOCATE(delta(1:sys%gr%np))

    ! analytic paramagnetic
    current_para_ref(:,:) = M_ZERO
    do ip = 1, sys%gr%np
      call mesh_r(sys%gr, ip, rr)
      xx = sys%gr%x(ip, 1:3)
      if (rr > R_SMALL) then
        current_para_ref(ip,1:3) = - (psi_2s(rr, a0) * dr_psi_1s(rr, a0)  - &
          psi_1s(rr, a0) * dr_psi_2s(rr, a0) ) * imag(alpha) / (1 + abs(alpha)**2) * xx(1:3) / rr
      end if
    end do

    write(fname, '(a)') 'current_para'
    out_how = io_function_fill_how("PlaneZ")
    call io_function_output_vector(out_how , './' , trim(fname), namespace, sys%space, sys%gr, &
      sys%st%current_para(:,:,1), unit_one, ierr)

    write(fname, '(a)') 'current_para-ref'
    out_how = io_function_fill_how("PlaneZ")
    call io_function_output_vector(out_how , './' , trim(fname), namespace, sys%space, sys%gr, &
      current_para_ref(:,:), unit_one, ierr)

    do idir = 1, 3
      delta = M_ZERO
      delta(:) = current_para_ref(1:sys%gr%np, idir) - sys%st%current_para(1:sys%gr%np, idir, 1)
      write(message(idir),*) 'idir =',idir,', norm2(delta paramagnetic)',norm2(delta)
    end do
    call messages_info(3)

    ! analytic diamagnetic
    current_dia_ref(:,:) = M_ZERO
    do ip = 1, sys%gr%np
      call mesh_r(sys%gr, ip, rr)
      current_dia_ref(ip,1:3) = - sys%hm%hm_base%vector_potential(1:3,ip) *&
        M_ONE / (1 + abs(alpha)**2) * abs(lc_hydrogen_state(rr, alpha, a0))**2
    end do

    write(fname, '(a)') 'current_dia'
    out_how = io_function_fill_how("PlaneZ")
    call io_function_output_vector(out_how , './' , trim(fname), namespace, sys%space, sys%gr, &
      sys%st%current_dia(:,:,1), unit_one, ierr)

    write(fname, '(a)') 'current_dia-ref'
    out_how = io_function_fill_how("PlaneZ")
    call io_function_output_vector(out_how , './' , trim(fname), namespace, sys%space, sys%gr, &
      current_dia_ref(:,:), unit_one, ierr)

    do idir = 1, 3
      delta = M_ZERO
      delta(:) = current_dia_ref(1:sys%gr%np, idir) - sys%st%current_dia(1:sys%gr%np, idir, 1)
      write(message(idir),*) 'idir =',idir,', norm2(delta diamagnetic)',norm2(delta)
    end do
    call messages_info(3)

    ! magnetization current
    call current_calculate_mag(current, namespace, sys%gr%der, sys%hm, sys%st)

    ! analytic magnetization
    current_mag_ref(:,:) = M_ZERO
    do ip = 1, sys%gr%np
      call mesh_r(sys%gr, ip, rr)
      xx = sys%gr%x(ip, 1:3)
      if (norm2(xx(1:2)) > R_SMALL .and. rr > R_SMALL) then
        sin_thet = norm2(xx(1:2)) / rr
        sin_phi = xx(2) / norm2(xx(1:2))
        cos_phi = xx(1) / norm2(xx(1:2))
        mag_curr = M_TWO * (psi_1s(rr,a0)*dr_psi_1s(rr,a0) + real(alpha)*(psi_1s(rr,a0)*dr_psi_2s(rr,a0) &
          + dr_psi_1s(rr,a0)*psi_2s(rr,a0)) + abs(alpha)**2 * psi_2s(rr,a0)*dr_psi_2s(rr,a0))
        ! minus signs are reversed because of the charge of the electron
        current_mag_ref(ip,1) = M_HALF * mag_curr * sin_thet * sin_phi / (1+abs(alpha)**2)
        current_mag_ref(ip,2) = -M_HALF * mag_curr * sin_thet * cos_phi / (1+abs(alpha)**2)
      endif
    end do

    write(fname, '(a)') 'current_mag'
    out_how = io_function_fill_how("PlaneZ")
    call io_function_output_vector(out_how , './' , trim(fname), namespace, sys%space, sys%gr, &
      sys%st%current_mag(:,:,1), unit_one, ierr)

    write(fname, '(a)') 'current_mag-ref'
    out_how = io_function_fill_how("PlaneZ")
    call io_function_output_vector(out_how , './' , trim(fname), namespace, sys%space, sys%gr, &
      current_mag_ref(:,:), unit_one, ierr)

    do idir = 1, 3
      delta = M_ZERO
      delta(:) = current_mag_ref(1:sys%gr%np, idir) - sys%st%current_mag(1:sys%gr%np, idir, 1)
      write(message(idir),*) 'idir =',idir,', norm2(delta magnetization)',norm2(delta)
    end do
    call messages_info(3)

    SAFE_DEALLOCATE_A(current_para_ref)
    SAFE_DEALLOCATE_A(current_dia_ref)
    SAFE_DEALLOCATE_A(current_mag_ref)
    SAFE_DEALLOCATE_A(delta)
    SAFE_DEALLOCATE_P(sys)

  end subroutine test_current_density


  subroutine set_hydrogen_states(psib, mesh, namespace, alpha, a0)
    class(batch_t), intent(inout) :: psib
    class(mesh_t),  intent(in)    :: mesh
    type(namespace_t), intent(in) :: namespace
    CMPLX,             intent(in) :: alpha
    FLOAT,             intent(in) :: a0

    CMPLX, allocatable :: zff(:)
    integer :: ip
    FLOAT :: rr

    PUSH_SUB(set_hydrogen_states)

    ! psi = phi_{1s} + \alpha \phi_{2s}
    SAFE_ALLOCATE(zff(1:mesh%np))
    if (type_is_complex(psib%type())) then
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr)
        zff(ip) = M_ONE / sqrt(1 + abs(alpha)**2) * lc_hydrogen_state(rr, alpha, a0)
        call batch_set_state(psib, 1, mesh%np, zff)
      end do
      SAFE_DEALLOCATE_A(zff)
    else
      write(message(1),*) "States should be complex for the linear combination of hydrogenic states to work"
      call messages_info(1, namespace=namespace)
    end if

    POP_SUB(set_hydrogen_states)
  end subroutine set_hydrogen_states

  CMPLX function lc_hydrogen_state(rr, alpha, a0)
    FLOAT, intent(in) :: a0, rr
    CMPLX, intent(in) :: alpha

    lc_hydrogen_state = psi_1s(rr, a0) + alpha * psi_2s(rr, a0)
  end function lc_hydrogen_state

  ! phi_{1s} = R_{10}(r)Y_{00}(\theta, \phi) = 1/(2\sqrt(\pi)) * 2 a0^{-3/2} * \exp(-r/a0)
  FLOAT function psi_1s(rr, a0)
    FLOAT, intent(in) :: a0, rr

    psi_1s = a0**(-M_THREE/M_TWO) * exp(-rr / a0) / sqrt(M_PI)
  end function psi_1s

  ! phi_{2s} = R_{20}(r)Y_{00}(\theta, \phi) = 1/(2\sqrt(\pi)) * \sqrt(2) a0^{-3/2}  (1 - r/(2a0)) \exp(-r/(2a0))
  FLOAT function psi_2s(rr, a0)
    FLOAT, intent(in) :: a0, rr

    psi_2s = sqrt(M_TWO) * a0**(-M_THREE/M_TWO) * (M_ONE - rr/(M_TWO * a0)) * exp(-rr/(M_TWO * a0)) &
      / (M_TWO * sqrt(M_PI))
  end function psi_2s

  FLOAT function dr_psi_1s(rr, a0)
    FLOAT, intent(in) :: a0, rr

    dr_psi_1s = -(M_ONE / a0) * psi_1s(rr, a0)
  end function dr_psi_1s

  FLOAT function dr_psi_2s(rr, a0)
    FLOAT, intent(in) :: a0, rr

    dr_psi_2s = -(M_HALF / a0) * psi_2s(rr, a0) - &
      a0**(-M_FIVE/M_TWO) * exp(-rr/(M_TWO * a0)) / (sqrt(M_TWO) * M_TWO * sqrt(M_PI))
  end function dr_psi_2s

end module test_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
