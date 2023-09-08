!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2009 X. Andrade
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


module electrons_oct_m
  use accel_oct_m
  use absorbing_boundaries_oct_m
  use algorithm_oct_m
  use algorithm_factory_oct_m
  use calc_mode_par_oct_m
  use classical_particles_oct_m
  use clock_oct_m
  use current_oct_m
  use current_to_mxll_field_oct_m
  use debug_oct_m
  use density_oct_m
  use elf_oct_m
  use energy_calc_oct_m
  use ext_partner_list_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use interactions_factory_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use kick_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use loct_oct_m
  use mesh_oct_m
  use messages_oct_m
  use modelmb_particles_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use mxll_field_to_medium_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use pes_oct_m
  use poisson_oct_m
  use potential_interpolation_oct_m
  use propagator_oct_m
  use propagator_base_oct_m
  use propagator_bomd_oct_m
  use propagator_aetrs_oct_m
  use propagator_elec_oct_m
  use propagator_exp_mid_oct_m
  use propagator_verlet_oct_m
  use propagation_ops_elec_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use regridding_oct_m
  use scf_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use stress_oct_m
  use sort_oct_m
  use system_oct_m
  use td_oct_m
  use td_write_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use xc_oct_m
  use xc_oep_oct_m
  use xc_interaction_oct_m

  implicit none

  private
  public ::               &
    electrons_t


  !> @brief Class describing the electron system
  !!
  !! This class describes a system of electrons and ions.
  !!
  !! \todo move the ions into their own ions_t class.
  type, extends(system_t) :: electrons_t
    ! Components are public by default
    type(ions_t),     pointer    :: ions => NULL() !< the ion component of the system
    type(grid_t)                 :: gr    !< the mesh
    type(states_elec_t)          :: st    !< the states
    type(v_ks_t)                 :: ks    !< the Kohn-Sham potentials
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators
    type(hamiltonian_elec_t)     :: hm    !< the Hamiltonian
    type(td_t)                   :: td    !< everything related to time propagation
    type(current_t)              :: current_calculator
    type(scf_t)                  :: scf   !< SCF for BOMD

    type(kpoints_t) :: kpoints                   !< the k-points

    logical :: generate_epot

    type(states_elec_t)          :: st_copy  !< copy of the states

    ! At the moment this is not treated as an external potential
    class(lasers_t), pointer :: lasers => null()      !< lasers
    class(gauge_field_t), pointer :: gfield => null()      !< gauge field

    ! List with all the external partners
    ! This will become a list of interactions in the future
    type(partner_list_t) :: ext_partners

    !TODO: have a list of self interactions
    type(xc_interaction_t), pointer   :: xc_interaction => null()

    logical :: ions_propagated = .false.

    logical :: needs_current = .false.   !< If current is needed by an interaction
  contains
    procedure :: init_interaction => electrons_init_interaction
    procedure :: init_parallelization => electrons_init_parallelization
    procedure :: init_algorithm => electrons_init_algorithm
    procedure :: initial_conditions => electrons_initial_conditions
    procedure :: do_algorithmic_operation => electrons_do_algorithmic_operation
    procedure :: is_tolerance_reached => electrons_is_tolerance_reached
    procedure :: update_quantity => electrons_update_quantity
    procedure :: update_exposed_quantity => electrons_update_exposed_quantity
    procedure :: init_interaction_as_partner => electrons_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => electrons_copy_quantities_to_interaction
    procedure :: output_start => electrons_output_start
    procedure :: output_write => electrons_output_write
    procedure :: output_finish => electrons_output_finish
    procedure :: process_is_slave  => electrons_process_is_slave
    procedure :: restart_write_data => electrons_restart_write_data
    procedure :: restart_read_data => electrons_restart_read_data
    procedure :: update_kinetic_energy => electrons_update_kinetic_energy
    procedure :: propagation_start => electrons_propagation_start
    final :: electrons_finalize
  end type electrons_t

  interface electrons_t
    procedure electrons_constructor
  end interface electrons_t

contains

  !----------------------------------------------------------
  function electrons_constructor(namespace, generate_epot) result(sys)
    class(electrons_t), pointer    :: sys
    type(namespace_t),  intent(in) :: namespace
    logical,  optional, intent(in) :: generate_epot

    integer :: iatom
    type(lattice_vectors_t) :: latt_inp
    type(profile_t), save :: prof

    PUSH_SUB(electrons_constructor)
    call profiling_in(prof,"ELECTRONS_CONSTRUCTOR")

    SAFE_ALLOCATE(sys)

    sys%namespace = namespace

    call messages_obsolete_variable(sys%namespace, 'SystemName')

    call space_init(sys%space, sys%namespace)
    call sys%space%write_info(sys%namespace)
    if (sys%space%has_mixed_periodicity()) then
      call messages_experimental('Support for mixed periodicity systems')
    end if

    sys%ions => ions_t(sys%namespace, latt_inp=latt_inp)

    call grid_init_stage_1(sys%gr, sys%namespace, sys%space, sys%ions%symm, latt_inp, sys%ions%natoms, sys%ions%pos)

    if (sys%space%is_periodic()) then
      call sys%ions%latt%write_info(sys%namespace)
    end if

    ! Sanity check for atomic coordinates
    do iatom = 1, sys%ions%natoms
      if (.not. sys%gr%box%contains_point(sys%ions%pos(:, iatom))) then
        if (sys%space%periodic_dim /= sys%space%dim) then
          ! FIXME: This could fail for partial periodicity systems
          ! because contains_point is too strict with atoms close to
          ! the upper boundary to the cell.
          write(message(1), '(a,i5,a)') "Atom ", iatom, " is outside the box."
          call messages_warning(1, namespace=sys%namespace)
        end if
      end if
    end do

    ! we need k-points for periodic systems
    call kpoints_init(sys%kpoints, sys%namespace, sys%gr%symm, sys%space%dim, sys%space%periodic_dim, sys%ions%latt)

    call states_elec_init(sys%st, sys%namespace, sys%space, sys%ions%val_charge(), sys%kpoints)
    call sys%st%write_info(sys%namespace)
    ! if independent particles in N dimensions are being used, need to initialize them
    !  after masses are set to 1 in grid_init_stage_1 -> derivatives_init
    call modelmb_copy_masses (sys%st%modelmbparticles, sys%gr%der%masses)
    call elf_init(sys%namespace)

    sys%generate_epot = optional_default(generate_epot, .true.)

    call sys%supported_interactions_as_partner%add(CURRENT_TO_MXLL_FIELD)
    call sys%supported_interactions%add(MXLL_FIELD_TO_MEDIUM)
    sys%quantities(CURRENT)%required = .true.
    sys%quantities(CURRENT)%updated_on_demand = .false.
    call current_init(sys%current_calculator, sys%namespace)

    call profiling_out(prof)
    POP_SUB(electrons_constructor)
  end function electrons_constructor

  ! ---------------------------------------------------------
  subroutine electrons_init_interaction(this, interaction)
    class(electrons_t), target, intent(inout) :: this
    class(interaction_t),       intent(inout) :: interaction

    FLOAT :: dmin
    integer :: rankmin

    PUSH_SUB(electrons_init_interactions)

    select type (interaction)
    type is (mxll_field_to_medium_t)
      call interaction%init(this%gr, 3)
      select case (this%hm%mxll_coupling_mode)
      case (LENGTH_GAUGE_DIPOLE)
        ASSERT(this%space%periodic_dim == 0)
        interaction%type = MXLL_FIELD_TRANS
      case (VELOCITY_GAUGE_DIPOLE, FULL_MINIMAL_COUPLING)
        interaction%type = MXLL_VEC_POT_TRANS
      case default
        message(1) = "Unknown Maxwell-matter coupling level."
        message(2) = "You should set the MaxwellCouplingMode variable."
        call messages_fatal(2, namespace=this%namespace)
      end select
      this%hm%center_of_mass(1:3) = this%ions%center_of_mass()
      this%hm%center_of_mass_ip = mesh_nearest_point(this%gr, this%hm%center_of_mass, dmin, rankmin)
    class default
      message(1) = "Trying to initialize an unsupported interaction by the electrons."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(electrons_init_interactions)
  end subroutine electrons_init_interaction

  ! ---------------------------------------------------------
  subroutine electrons_init_parallelization(this, grp)
    class(electrons_t), intent(inout) :: this
    type(mpi_grp_t),    intent(in)    :: grp

    integer(i8) :: index_range(4)
    FLOAT :: mesh_global, mesh_local, wfns
    integer :: idir
    FLOAT :: spiral_q(3), spiral_q_red(3)
    type(block_t) :: blk

    PUSH_SUB(electrons_init_parallelization)

    call mpi_grp_copy(this%grp, grp)

    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%np_global  ! Number of points in mesh
    index_range(2) = this%st%nst             ! Number of states
    index_range(3) = this%st%d%nik           ! Number of k-points
    index_range(4) = 100000                 ! Some large number

    ! create index and domain communicators
    call multicomm_init(this%mc, this%namespace, this%grp, calc_mode_par_parallel_mask(), calc_mode_par_default_parallel_mask(), &
      mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

    call this%ions%partition(this%mc)
    call kpoints_distribute(this%st%d, this%mc)
    call states_elec_distribute_nodes(this%st, this%namespace, this%mc)


    if (parse_is_defined(this%namespace, 'TDMomentumTransfer') .or. &
      parse_is_defined(this%namespace, 'TDReducedMomentumTransfer')) then
      if (parse_block(this%namespace, 'TDMomentumTransfer', blk) == 0) then
        do idir = 1, this%space%dim
          call parse_block_float(blk, 0, idir - 1, spiral_q(idir))
        end do
      else if(parse_block(this%namespace, 'TDReducedMomentumTransfer', blk) == 0) then
        do idir = 1, this%space%dim
          call parse_block_float(blk, 0, idir - 1, spiral_q_red(idir))
        end do
        call kpoints_to_absolute(this%kpoints%latt, spiral_q_red, spiral_q)
      end if
      call parse_block_end(blk)
      call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc, spiral_q)
    else
      call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc)
    end if

    if (this%st%symmetrize_density) then
      call mesh_check_symmetries(this%gr, this%gr%symm, this%ions%space%periodic_dim)
    end if

    call output_init(this%outp, this%namespace, this%space, this%st, this%st%nst, this%ks)
    call states_elec_densities_init(this%st, this%gr)
    call states_elec_exec_init(this%st, this%namespace, this%mc)

    call v_ks_init(this%ks, this%namespace, this%gr, this%st, this%ions, this%mc, this%space, this%kpoints)
    if (this%ks%theory_level == KOHN_SHAM_DFT .or. this%ks%theory_level == GENERALIZED_KOHN_SHAM_DFT) then
      this%xc_interaction => xc_interaction_t(this)
    end if

    ! Temporary place for the initialization of the lasers
    this%lasers => lasers_t(this%namespace, this%gr, this%space, this%ions%latt)
    if(this%lasers%no_lasers > 0) then
      call this%ext_partners%add(this%lasers)
      call lasers_check_symmetries(this%lasers, this%kpoints)
    else
      deallocate(this%lasers)
    end if

    ! Temporary place for the initialization of the gauge field
    this%gfield => gauge_field_t(this%namespace, this%ions%latt%rcell_volume)
    if(gauge_field_is_used(this%gfield)) then
      call this%ext_partners%add(this%gfield)
      call gauge_field_check_symmetries(this%gfield, this%kpoints)
    else
      deallocate(this%gfield)
    end if

    call hamiltonian_elec_init(this%hm, this%namespace, this%space, this%gr, this%ions, this%ext_partners, &
      this%st, this%ks%theory_level, this%ks%xc, this%mc, this%kpoints, &
      need_exchange = output_need_exchange(this%outp) .or. this%ks%oep%level /= OEP_LEVEL_NONE)

    if (this%hm%pcm%run_pcm .and. this%mc%par_strategy /= P_STRATEGY_SERIAL .and. this%mc%par_strategy /= P_STRATEGY_STATES) then
      call messages_experimental('Parallel in domain calculations with PCM')
    end if

    ! Print memory requirements
    call messages_print_with_emphasis(msg='Approximate memory requirements', namespace=this%namespace)

    mesh_global = mesh_global_memory(this%gr)
    mesh_local  = mesh_local_memory(this%gr)

    call messages_write('Mesh')
    call messages_new_line()
    call messages_write('  global  :')
    call messages_write(mesh_global, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_write('  local   :')
    call messages_write(mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_write('  total   :')
    call messages_write(mesh_global + mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_info()

    wfns = states_elec_wfns_memory(this%st, this%gr)
    call messages_write('States')
    call messages_new_line()
    call messages_write('  real    :')
    call messages_write(wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()
    call messages_write('  complex :')
    call messages_write(2.0_8*wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()
    call messages_info()

    call messages_print_with_emphasis(namespace=this%namespace)

    if (this%generate_epot) then
      message(1) = "Info: Generating external potential"
      call messages_info(1, namespace=this%namespace)
      call hamiltonian_elec_epot_generate(this%hm, this%namespace, this%space, this%gr, this%ions, &
        this%ext_partners, this%st)
      message(1) = "      done."
      call messages_info(1, namespace=this%namespace)
    end if

    if (this%ks%theory_level /= INDEPENDENT_PARTICLES) then
      call poisson_async_init(this%hm%psolver, this%mc)
      ! slave nodes do not call the calculation routine
      if (multicomm_is_slave(this%mc)) then
        !for the moment we only have one type of slave
        call poisson_slave_work(this%hm%psolver, this%namespace)
      end if
    end if


    POP_SUB(electrons_init_parallelization)
  end subroutine electrons_init_parallelization

  ! ---------------------------------------------------------
  subroutine electrons_init_algorithm(this, factory)
    class(electrons_t),         intent(inout) :: this
    class(algorithm_factory_t), intent(in)    :: factory

    integer :: depth, i_interaction
    type(interaction_iterator_t) :: iter
    class(interaction_t), pointer :: interaction
    character(len=256) :: label

    PUSH_SUB(electrons_init_algorithm)

    call system_init_algorithm(this, factory)
    ! TODO: once GS algorithms are supported, this should only be executed for TD algorithms
    call td_init(this%td, this%namespace, this%space, this%gr, this%ions, this%st, this%ks, &
      this%hm, this%ext_partners, this%outp)

    ! this corresponds to the first part of td_init_run
    call td_allocate_wavefunctions(this%td, this%namespace, this%mc, this%gr, this%ions, this%st, &
      this%ks, this%hm, this%space)
    call td_init_gaugefield(this%td, this%namespace, this%gr, this%st, this%ks, this%hm, &
      this%ext_partners, this%space)

    ! interpolation depth depends on the propagator
    select type (prop => this%algo)
    type is (propagator_exp_mid_t)
      depth = 3
    type is (propagator_aetrs_t)
      depth = 3
    type is (propagator_bomd_t)
      depth = 1
    class default
      message(1) = "Propagator does not yet support initializing of interaction interpolation"
      call messages_fatal(1)
    end select
    ! set interpolation depth for interactions
    i_interaction = 0
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      select type (interaction)
      class is (mxll_field_to_medium_t)
        write(label, "(A, I5.5)") "mxll_field_to_medium_", i_interaction
        call interaction%init_interpolation(depth, trim(label))
      end select
      i_interaction = i_interaction + 1
    end do

    POP_SUB(electrons_init_algorithm)
  end subroutine electrons_init_algorithm

  ! ---------------------------------------------------------
  subroutine electrons_initial_conditions(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_initial_conditions)

    call td_set_from_scratch(this%td, .true.)
    call td_load_restart_from_gs(this%td, this%namespace, this%space, this%mc, this%gr, this%ions, &
      this%ext_partners, this%st, this%ks, this%hm)

    POP_SUB(electrons_initial_conditions)
  end subroutine electrons_initial_conditions

  ! ---------------------------------------------------------
  subroutine electrons_propagation_start(this)
    class(electrons_t),      intent(inout) :: this

    PUSH_SUB(electrons_propagation_start)

    call system_propagation_start(this)

    ! additional initialization needed for electrons
    call td_init_with_wavefunctions(this%td, this%namespace, this%space, this%mc, this%gr, this%ions, &
      this%ext_partners, this%st, this%ks, this%hm, this%outp, td_get_from_scratch(this%td))

    POP_SUB(electrons_propagation_start)
  end subroutine electrons_propagation_start

  ! ---------------------------------------------------------
  logical function electrons_do_algorithmic_operation(this, operation) result(done)
    class(electrons_t),             intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    logical :: update_energy_
    type(gauge_field_t), pointer :: gfield

    PUSH_SUB(electrons_do_algorithmic_operation)

    update_energy_ = .true.

    ! kick at t > 0 still missing!

    done = .true.
    select case (operation%id)
    case (AETRS_FIRST_HALF)
      ! propagate half of the time step with H(time)
      call get_fields_from_interaction(this, this%algo%clock%time(), this%hm%field_mxll)
      call propagation_ops_elec_update_hamiltonian(this%namespace, this%space, this%st, this%gr, &
        this%hm, this%ext_partners, this%algo%clock%time())
      call propagation_ops_elec_exp_apply(this%td%tr%te, this%namespace, this%st, this%gr, this%hm, M_HALF*this%algo%dt)

    case (AETRS_EXTRAPOLATE)
      ! Do the extrapolation of the Hamiltonian
      ! First the extrapolation of the potential
      call potential_interpolation_new(this%td%tr%vksold, this%gr%np, this%st%d%nspin, &
        this%algo%clock%time()+this%algo%dt, this%algo%dt, this%hm%vhxc, vtau=this%hm%vtau)

      !Get the potentials from the interpolator
      call propagation_ops_elec_interpolate_get(this%gr, this%hm, this%td%tr%vksold)

      ! Second the ions and gauge field which later on will be treated as extrapolation
      ! of interactions, but this is not yet possible

      ! move the ions to time t
      call propagation_ops_elec_move_ions(this%td%tr%propagation_ops_elec, this%gr, this%hm, &
        this%st, this%namespace, this%space, this%td%ions_dyn, this%ions, this%ext_partners, &
        this%algo%clock%time()+this%algo%dt, this%algo%dt)

      !Propagate gauge field
      gfield => list_get_gauge_field(this%ext_partners)
      if(associated(gfield)) then
        call propagation_ops_elec_propagate_gauge_field(this%td%tr%propagation_ops_elec, gfield, &
          this%algo%dt, this%algo%clock%time()+this%algo%dt)
      end if

      !Update Hamiltonian and current
      call get_fields_from_interaction(this, this%algo%clock%time()+this%algo%dt, this%hm%field_mxll)
      call propagation_ops_elec_update_hamiltonian(this%namespace, this%space, this%st, this%gr, &
        this%hm, this%ext_partners, this%algo%clock%time()+this%algo%dt)

    case (AETRS_SECOND_HALF)
      !Do the time propagation for the second half of the time step
      call propagation_ops_elec_fuse_density_exp_apply(this%td%tr%te, this%namespace, this%st, &
        this%gr, this%hm, M_HALF*this%algo%dt)
      call calc_current_for_interaction(this)

    case (EXPMID_EXTRAPOLATE)
      ! the half step of this propagator screws with the gauge field kick
      gfield => list_get_gauge_field(this%ext_partners)
      if(associated(gfield)) then
        ASSERT(gauge_field_is_propagated(gfield) .eqv. .false.)
      end if

      ! Do the extrapolation of the Hamiltonian
      ! First the extrapolation of the potential
      call potential_interpolation_new(this%td%tr%vksold, this%gr%np, this%st%d%nspin, &
        this%algo%clock%time()+this%algo%dt, this%algo%dt, this%hm%vhxc, vtau=this%hm%vtau)

      ! get the potentials from the interpolator
      if (this%hm%theory_level /= INDEPENDENT_PARTICLES) then
        call potential_interpolation_interpolate(this%td%tr%vksold, 3, &
          this%algo%clock%time()+this%algo%dt, this%algo%dt, this%algo%clock%time() + this%algo%dt/M_TWO, &
          this%hm%vhxc, vtau=this%hm%vtau)
      end if

      ! Second the ions which later on will be treated as extrapolation of interactions,
      ! but this is not yet possible
      ! move the ions to the half step
      call propagation_ops_elec_move_ions(this%td%tr%propagation_ops_elec, this%gr, this%hm, this%st, &
        this%namespace, this%space, this%td%ions_dyn, this%ions, this%ext_partners, &
        this%algo%clock%time() + M_HALF*this%algo%dt, M_HALF*this%algo%dt, save_pos=.true.)

      call get_fields_from_interaction(this, this%algo%clock%time() + M_HALF*this%algo%dt, this%hm%field_mxll)
      call propagation_ops_elec_update_hamiltonian(this%namespace, this%space, this%st, this%gr, &
        this%hm, this%ext_partners, this%algo%clock%time() + M_HALF*this%algo%dt)

    case (EXPMID_PROPAGATE)
      ! Do the actual propagation step
      call propagation_ops_elec_fuse_density_exp_apply(this%td%tr%te, this%namespace, this%st, &
        this%gr, this%hm, this%algo%dt)
      call calc_current_for_interaction(this)

      ! restore to previous time
      call propagation_ops_elec_restore_ions(this%td%tr%propagation_ops_elec, this%td%ions_dyn, this%ions)

    case (BOMD_START)
      call scf_init(this%scf, this%namespace, this%gr, this%ions, this%st, this%mc, this%hm, this%ks, this%space)
      ! the ions are propagated inside the propagation step already, so no need to do it at the end
      this%ions_propagated = .true.

    case (VERLET_UPDATE_POS)
      ! move the ions to time t
      call ion_dynamics_propagate(this%td%ions_dyn, this%ions, this%algo%clock%time()+this%algo%dt, &
        this%algo%dt, this%namespace)

    case (BOMD_ELEC_SCF)
      call hamiltonian_elec_epot_generate(this%hm, this%namespace, this%space, this%gr, this%ions, &
        this%ext_partners, this%st, time = this%algo%clock%time()+this%algo%dt)
      ! now calculate the eigenfunctions
      call scf_run(this%scf, this%namespace, this%space, this%mc, this%gr, this%ions, &
        this%ext_partners, this%st, this%ks, this%hm, this%outp, gs_run = .false., verbosity = VERB_COMPACT)
      ! TODO: Check if this call is realy needed. - NTD
      call hamiltonian_elec_epot_generate(this%hm, this%namespace, this%space, this%gr, this%ions, &
        this%ext_partners, this%st, time = this%algo%clock%time()+this%algo%dt)

      ! update Hamiltonian and eigenvalues (fermi is *not* called)
      call v_ks_calc(this%ks, this%namespace, this%space, this%hm, this%st, this%ions, this%ext_partners, &
        calc_eigenval = .true., time = this%algo%clock%time()+this%algo%dt, calc_energy = .true.)

      ! Get the energies.
      call energy_calc_total(this%namespace, this%space, this%hm, this%gr, this%st, this%ext_partners, iunit = -1)

    case (VERLET_COMPUTE_ACC)
      ! Do nothing, forces are computing in scf_run

    case (VERLET_COMPUTE_VEL)
      call ion_dynamics_propagate_vel(this%td%ions_dyn, this%ions)
      ! TODO: Check if this call is realy needed. - NTD
      call hamiltonian_elec_epot_generate(this%hm, this%namespace, this%space, this%gr, this%ions, &
        this%ext_partners, this%st, time = this%algo%clock%time()+this%algo%dt)
      call this%ions%update_kinetic_energy()

    case (EXPMID_START)
      this%ions_propagated = .false.

    case (AETRS_START)
      ! the ions are propagated inside the propagation step already, so no need to do it at the end
      this%ions_propagated = .true.

    case (STEP_DONE)
      call electrons_exec_end_of_timestep_tasks(this)
      done = .false.

    case (BOMD_FINISH)
      call scf_end(this%scf)

    case (EXPMID_FINISH, AETRS_FINISH)
    case default
      done = .false.
    end select

    POP_SUB(electrons_do_algorithmic_operation)
  end function electrons_do_algorithmic_operation

  ! ---------------------------------------------------------
  logical function electrons_is_tolerance_reached(this, tol) result(converged)
    class(electrons_t), intent(in) :: this
    FLOAT,              intent(in) :: tol

    PUSH_SUB(electrons_is_tolerance_reached)

    converged = .false.

    POP_SUB(electrons_is_tolerance_reached)
  end function electrons_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine electrons_update_quantity(this, iq)
    class(electrons_t),   intent(inout) :: this
    integer,              intent(in)    :: iq

    PUSH_SUB(electrons_update_quantity)

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(this%quantities(iq)%updated_on_demand)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(electrons_update_quantity)
  end subroutine electrons_update_quantity

  ! ---------------------------------------------------------
  subroutine electrons_update_exposed_quantity(partner, iq)
    class(electrons_t), intent(inout) :: partner
    integer,            intent(in)    :: iq

    PUSH_SUB(electrons_update_exposed_quantity)

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(partner%quantities(iq)%updated_on_demand)

    select case (iq)
    case (CURRENT)
      ! Nothing to do, it is now done in calc_current_for_interaction
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(electrons_update_exposed_quantity)
  end subroutine electrons_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine electrons_init_interaction_as_partner(partner, interaction)
    class(electrons_t),   intent(in)    :: partner
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(electrons_init_interaction_as_partner)

    select type (interaction)
    type is (current_to_mxll_field_t)
      call interaction%init_from_partner(partner%gr, partner%space, partner%namespace)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(electrons_init_interaction_as_partner)
  end subroutine electrons_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine electrons_copy_quantities_to_interaction(partner, interaction)
    class(electrons_t),   intent(inout) :: partner
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(electrons_copy_quantities_to_interaction)

    select type (interaction)
    type is (current_to_mxll_field_t)
      partner%needs_current = .true.
      if (allocated(partner%st%current)) then
        interaction%partner_field(:,:) = partner%st%current(1:partner%gr%np,:,1)
      end if
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(electrons_copy_quantities_to_interaction)
  end subroutine electrons_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine electrons_output_start(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_output_start)

    POP_SUB(electrons_output_start)
  end subroutine electrons_output_start

  ! ---------------------------------------------------------
  subroutine electrons_output_write(this)
    class(electrons_t), intent(inout) :: this

    integer :: iter

    PUSH_SUB(electrons_output_write)

    iter = this%clock%get_tick()

    call td_write_iter(this%td%write_handler, this%namespace, this%space, this%outp, this%gr, &
      this%st, this%hm, this%ions, this%ext_partners, this%hm%kick, this%ks, this%td%dt, iter)

    if (this%outp%anything_now(iter)) then ! output
      call td_write_output(this%namespace, this%space, this%gr, this%st, this%hm, this%ks, &
        this%outp, this%ions, this%ext_partners, iter, this%td%dt)
    end if

    POP_SUB(electrons_output_write)
  end subroutine electrons_output_write

  ! ---------------------------------------------------------
  subroutine electrons_output_finish(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_output_finish)

    POP_SUB(electrons_output_finish)
  end subroutine electrons_output_finish

  ! ---------------------------------------------------------
  logical function electrons_process_is_slave(this) result(is_slave)
    class(electrons_t), intent(in) :: this

    PUSH_SUB(electrons_process_is_slave)

    is_slave = multicomm_is_slave(this%mc)

    POP_SUB(electrons_process_is_slave)
  end function electrons_process_is_slave

  ! ---------------------------------------------------------
  subroutine electrons_exec_end_of_timestep_tasks(this)
    class(electrons_t), intent(inout) :: this
    logical :: stopping
    logical :: generate
    logical :: update_energy_
    integer :: nt
    FLOAT :: dt, time
    type(gauge_field_t), pointer :: gfield

    PUSH_SUB(electrons_exec_end_of_timestep_tasks)

    stopping = .false.

    dt = this%algo%dt
    nt = this%td%iter
    ! this is the time at the end of the timestep, as required in all routines here
    time = dt*nt
    update_energy_ = .true.

    !Apply mask absorbing boundaries
    if (this%hm%abs_boundaries%abtype == MASK_ABSORBING) call zvmask(this%gr, this%hm, this%st)

    !Photoelectron stuff
    if (this%td%pesv%calc_spm .or. this%td%pesv%calc_mask .or. this%td%pesv%calc_flux) then
      call pes_calc(this%td%pesv, this%namespace, this%space, this%gr, this%st, &
        this%td%dt, nt, this%gr%der, this%hm%kpoints, this%ext_partners, stopping)
    end if

    ! For BOMD, we do not want the lines below to be executed
    select type(prop => this%algo)
    type is(propagator_bomd_t)
      POP_SUB(electrons_exec_end_of_timestep_tasks)
      return
    end select

    ! The propagation of the ions and the gauge field is currently done here.
    ! TODO: this code is to be moved to their own systems at some point
    generate = .false.
    if (ion_dynamics_ions_move(this%td%ions_dyn)) then
      if (.not. this%ions_propagated) then
        call ion_dynamics_propagate(this%td%ions_dyn, this%ions, abs(nt*dt), this%td%ions_dyn%ionic_scale*dt, &
          this%namespace)
        generate = .true.
      end if
    end if

    gfield => list_get_gauge_field(this%ext_partners)
    if(associated(gfield)) then
      if (gauge_field_is_propagated(gfield) .and. .not. this%ions_propagated) then
        call gauge_field_do_algorithmic_operation(gfield, OP_VERLET_COMPUTE_ACC, dt, time)
      end if
    end if

    if (generate .or. this%ions%has_time_dependent_species()) then
      call hamiltonian_elec_epot_generate(this%hm, this%namespace, this%space, this%gr, this%ions, &
        this%ext_partners, this%st, time = abs(nt*dt))
    end if

    call v_ks_calc(this%ks, this%namespace, this%space, this%hm, this%st, this%ions, this%ext_partners, &
      calc_eigenval = update_energy_, time = abs(nt*dt), calc_energy = update_energy_)

    if (update_energy_) then
      call energy_calc_total(this%namespace, this%space, this%hm, this%gr, this%st, this%ext_partners, iunit = -1)
    end if

    ! Recalculate forces, update velocities...
    if (ion_dynamics_ions_move(this%td%ions_dyn)) then
      call forces_calculate(this%gr, this%namespace, this%ions, this%hm, this%ext_partners, this%st, &
        this%ks, t = abs(nt*dt), dt = dt)
      call ion_dynamics_propagate_vel(this%td%ions_dyn, this%ions, atoms_moved = generate)
      call this%ions%update_kinetic_energy()
    else
      if (this%outp%what(OPTION__OUTPUT__FORCES) .or. this%td%write_handler%out(OUT_SEPARATE_FORCES)%write) then
        call forces_calculate(this%gr, this%namespace, this%ions, this%hm, this%ext_partners, &
          this%st, this%ks, t = abs(nt*dt), dt = dt)
      end if
    end if

    if (this%outp%what(OPTION__OUTPUT__STRESS)) then
      call stress_calculate(this%namespace, this%gr, this%hm, this%st, this%ions, this%ks, this%ext_partners)
    end if

    if(associated(gfield)) then
      if(gauge_field_is_propagated(gfield)) then
        if(this%ks%xc%kernel_lrc_alpha>M_EPSILON) then
          call gauge_field_get_force(gfield, this%gr, this%st%d%spin_channels, this%st%current, this%ks%xc%kernel_lrc_alpha)
        else
          call gauge_field_get_force(gfield, this%gr, this%st%d%spin_channels, this%st%current)
        endif
        call gauge_field_do_algorithmic_operation(gfield, OP_VERLET_COMPUTE_VEL, dt, time)
      end if
    end if

    !We update the occupation matrices
    call lda_u_update_occ_matrices(this%hm%lda_u, this%namespace, this%gr, this%st, this%hm%hm_base, this%hm%energy)

    ! this is needed to be compatible with the code in td_*
    this%td%iter = this%td%iter + 1

    POP_SUB(electrons_exec_end_of_timestep_tasks)
  end subroutine electrons_exec_end_of_timestep_tasks

  ! ---------------------------------------------------------
  subroutine electrons_restart_write_data(this)
    class(electrons_t), intent(inout) :: this

    integer :: ierr

    PUSH_SUB(electrons_restart_write_data)

    call td_write_data(this%td%write_handler)
    call td_dump(this%td, this%namespace, this%space, this%gr, this%st, this%hm, &
      this%ks, this%ext_partners, this%clock%get_tick(), ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write time-dependent restart information."
      call messages_warning(1, namespace=this%namespace)
    end if

    ! TODO: this is here because of legacy reasons and should be moved to the output framework
    call pes_output(this%td%pesv, this%namespace, this%space, this%gr, this%st, this%clock%get_tick(), &
      this%outp, this%td%dt, this%ions)

    POP_SUB(electrons_restart_write_data)
  end subroutine electrons_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function electrons_restart_read_data(this)
    class(electrons_t), intent(inout) :: this

    logical :: from_scratch

    PUSH_SUB(electrons_restart_read_data)

    from_scratch = .false.
    call td_load_restart_from_td(this%td, this%namespace, this%space, this%mc, this%gr, this%ions, &
      this%ext_partners, this%st, this%ks, this%hm, from_scratch)
    call td_set_from_scratch(this%td, from_scratch)
    if (from_scratch) then
      ! restart data could not be loaded
      electrons_restart_read_data = .false.
    else
      ! restart data could be loaded
      electrons_restart_read_data = .true.
    end if

    POP_SUB(electrons_restart_read_data)
  end function electrons_restart_read_data

  !----------------------------------------------------------
  subroutine electrons_update_kinetic_energy(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_update_kinetic_energy)

    if (states_are_real(this%st)) then
      this%kinetic_energy = denergy_calc_electronic(this%namespace, this%hm, this%gr%der, this%st, terms = TERM_KINETIC)
    else
      this%kinetic_energy = zenergy_calc_electronic(this%namespace, this%hm, this%gr%der, this%st, terms = TERM_KINETIC)
    end if

    POP_SUB(electrons_update_kinetic_energy)

  end subroutine electrons_update_kinetic_energy

  ! ---------------------------------------------------------
  subroutine get_fields_from_interaction(this, time, fields)
    class(electrons_t), intent(inout) :: this
    FLOAT,              intent(in)    :: time
    FLOAT,              intent(inout) :: fields(:, :)

    type(interaction_iterator_t) :: iter
    FLOAT, allocatable :: field_tmp(:, :)

    PUSH_SUB(get_fields_from_interaction)
    SAFE_ALLOCATE(field_tmp(1:this%gr%np, 1:this%gr%box%dim))
    fields = M_ZERO

    ! interpolate field from interaction
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (mxll_field_to_medium_t)
        call interaction%interpolate(time, field_tmp)
        call lalg_axpy(this%gr%np, 3, M_ONE, field_tmp, fields)
      end select
    end do

    SAFE_DEALLOCATE_A(field_tmp)
    POP_SUB(get_fields_from_interaction)

  end subroutine get_fields_from_interaction

  ! ---------------------------------------------------------
  subroutine calc_current_for_interaction(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(calc_current_for_interaction)

    if (this%needs_current) then
      call states_elec_allocate_current(this%st, this%space, this%gr)
      call current_calculate(this%current_calculator, this%namespace, this%gr,&
        this%hm, this%space, this%st)
      this%quantities(CURRENT)%clock = this%quantities(CURRENT)%clock + CLOCK_TICK
    end if

    POP_SUB(calc_current_for_interaction)
  end subroutine calc_current_for_interaction

  !----------------------------------------------------------
  subroutine electrons_finalize(sys)
    type(electrons_t), intent(inout) :: sys

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner

    PUSH_SUB(electrons_finalize)

    if (associated(sys%algo)) then
      call td_end_run(sys%td, sys%st, sys%hm)
      call td_end(sys%td)
    end if

    if (sys%ks%theory_level /= INDEPENDENT_PARTICLES) then
      call poisson_async_end(sys%hm%psolver, sys%mc)
    end if

    call iter%start(sys%ext_partners)
    do while (iter%has_next())
      partner => iter%get_next()
      SAFE_DEALLOCATE_P(partner)
    end do
    call sys%ext_partners%empty()

    SAFE_DEALLOCATE_P(sys%xc_interaction)

    call hamiltonian_elec_end(sys%hm)

    nullify(sys%gfield)
    nullify(sys%lasers)

    call multicomm_end(sys%mc)

    call v_ks_end(sys%ks)

    call states_elec_end(sys%st)

    SAFE_DEALLOCATE_P(sys%ions)

    call kpoints_end(sys%kpoints)

    call grid_end(sys%gr)

    call system_end(sys)

    POP_SUB(electrons_finalize)
  end subroutine electrons_finalize

end module electrons_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
