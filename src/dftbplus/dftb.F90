!! Copyright (C) 2020 F. Bonafé, H. Appel
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

module dftb_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
#ifdef HAVE_DFTBPLUS
  use dftbplus
#endif
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use iso_c_binding
  use lasers_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use species_oct_m
  use system_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::           &
    dftb_t,    &
    dftb_init

  integer, parameter ::     &
    OUTPUT_COORDINATES = 1, &
    OUTPUT_FORCES      = 2

  !> @brief class for a tight binding
  !!
  type, extends(system_t) :: dftb_t
    integer :: n_atom
    FLOAT, allocatable :: coords(:,:), gradients(:,:)
    FLOAT, allocatable :: acc(:,:)
    FLOAT, allocatable :: tot_force(:,:)
    FLOAT, allocatable :: vel(:,:)
    FLOAT, allocatable :: prev_tot_force(:,:) !< Used for the SCF convergence criterium
    integer, allocatable :: species(:)
    integer              :: dynamics
    FLOAT, allocatable :: mass(:)
    FLOAT, allocatable :: atom_charges(:,:) !< shape is (n_atoms, n_spin)
    character(len=LABEL_LEN), allocatable  :: labels(:)
    FLOAT, allocatable :: prev_acc(:,:,:) !< A storage of the prior times.
    FLOAT :: scc_tolerance
    type(ions_t),     pointer :: ions => NULL()
    type(c_ptr) :: output_handle(2)
    type(ion_dynamics_t) :: ions_dyn
    class(lasers_t), pointer :: lasers => null()
    logical :: laser_field
    FLOAT :: field(3)
    FLOAT :: energy
#ifdef HAVE_DFTBPLUS
    type(TDftbPlus) :: dftbp
#endif
  contains
    procedure :: init_interaction => dftb_init_interaction
    procedure :: initial_conditions => dftb_initial_conditions
    procedure :: do_algorithmic_operation => dftb_do_algorithmic_operation
    procedure :: output_start => dftb_output_start
    procedure :: output_write => dftb_output_write
    procedure :: output_finish => dftb_output_finish
    procedure :: is_tolerance_reached => dftb_is_tolerance_reached
    procedure :: update_quantity => dftb_update_quantity
    procedure :: update_exposed_quantity => dftb_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => dftb_copy_quantities_to_interaction
    procedure :: init_interaction_as_partner => dftb_init_interaction_as_partner
    procedure :: update_interactions_start => dftb_update_interactions_start
    procedure :: update_interactions_finish => dftb_update_interactions_finish
    procedure :: restart_write_data => dftb_restart_write_data
    procedure :: restart_read_data => dftb_restart_read_data
    procedure :: update_kinetic_energy => dftb_update_kinetic_energy
    final :: dftb_finalize
  end type dftb_t

  interface dftb_t
    procedure dftb_constructor
  end interface dftb_t

  !> Parameters.
  integer, parameter :: &
    EHRENFEST = 1,      &
    BO        = 2

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function dftb_constructor(namespace) result(sys)
    class(dftb_t),       pointer    :: sys
    type(namespace_t),   intent(in) :: namespace

    PUSH_SUB(dftb_constructor)

    SAFE_ALLOCATE(sys)

    call dftb_init(sys, namespace)

    POP_SUB(dftb_constructor)
  end function dftb_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine dftb_init(this, namespace)
    class(dftb_t),     target, intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace

    integer :: ii, jj, ispec, il, ierr
    character(len=MAX_PATH_LEN) :: slako_dir
    character(len=1), allocatable  :: max_ang_mom(:)
    character(len=LABEL_LEN) :: this_max_ang_mom, this_label
    integer :: n_maxang_block
    type(block_t) :: blk
    character(len=200) :: envelope_expression, phase_expression
    FLOAT :: omega0, initial_temp
    CMPLX :: pol(3)
    type(tdf_t) :: ff, phi

#ifdef HAVE_DFTBPLUS
    type(TDftbPlusInput) :: input
    type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pAnalysis
    type(fnode), pointer :: pParserOpts
#endif
#ifdef HAVE_DFTBPLUS_DEVEL
    type(fnode), pointer :: pElecDyn, pPerturb, pLaser
#endif

    PUSH_SUB(dftb_init)

    this%namespace = namespace

    call messages_print_with_emphasis(msg="DFTB+ System", namespace=namespace)

    call space_init(this%space, namespace)
    if (this%space%is_periodic()) then
      call messages_not_implemented('DFTB+ for periodic systems', namespace=namespace)
    end if
    if (this%space%dim /= 3) then
      call messages_not_implemented('DFTB+ for Dimensions /= 3', namespace=namespace)
    end if

    this%ions => ions_t(namespace)
    this%n_atom = this%ions%natoms
    SAFE_ALLOCATE(this%coords(1:3, 1:this%n_atom))
    SAFE_ALLOCATE(this%acc(1:3, 1:this%n_atom))
    SAFE_ALLOCATE(this%vel(1:3, 1:this%n_atom))
    SAFE_ALLOCATE(this%tot_force(1:3, 1:this%n_atom))
    SAFE_ALLOCATE(this%prev_tot_force(1:3, 1:this%n_atom))
    SAFE_ALLOCATE(this%gradients(1:3, 1:this%n_atom))
    SAFE_ALLOCATE(this%species(1:this%n_atom))
    SAFE_ALLOCATE(this%mass(1:this%n_atom))
    SAFE_ALLOCATE(this%atom_charges(1:this%n_atom, 1))
    SAFE_ALLOCATE(this%labels(1:this%ions%nspecies))
    SAFE_ALLOCATE(max_ang_mom(1:this%ions%nspecies))

    ispec = 1
    this%species(1) = 1
    this%labels(1) = trim(this%ions%atom(1)%label)

    this%coords = this%ions%pos
    ! mass is read from the default pseudopotential files
    this%mass = this%ions%mass
    do ii = 1, this%n_atom
      if ((ii > 1) .and. .not. (any(this%labels(1:ispec) == this%ions%atom(ii)%label))) then
        ispec = ispec + 1
        this%labels(ispec) = trim(this%ions%atom(ii)%label)
      end if
      do jj = 1, ispec
        if (trim(this%ions%atom(ii)%label) == trim(this%labels(jj))) then
          this%species(ii) = jj
        end if
      end do
    end do
    this%vel = M_ZERO
    this%tot_force = M_ZERO

    !%Variable MaxAngularMomentum
    !%Type block
    !%Section DFTBPlusInterface
    !%Description
    !% Specifies the highest angular momentum for each atom type. All orbitals up
    !% to that angular momentum will be included in the calculation.
    !% Possible values for the angular momenta are s, p, d, f.
    !% These are examples:
    !%
    !% <tt>%MaxAngularMomentum
    !% <br>&nbsp;&nbsp;'O'   | 'p'
    !% <br>&nbsp;&nbsp;'H'   | 's'
    !% <br>%</tt>
    !%End
    n_maxang_block = 0
    if (parse_block(namespace, 'MaxAngularMomentum', blk) == 0) then
      n_maxang_block = parse_block_n(blk)
      if (n_maxang_block /= this%ions%nspecies) then
        call messages_input_error(namespace, "MaxAngularMomentum", "Wrong number of species.")
      end if

      do ii = 1, n_maxang_block
        call parse_block_string(blk, ii-1, 0, this_label)
        call parse_block_string(blk, ii-1, 1, this_max_ang_mom)
        if (any(["s","p","d","f"] == trim(this_max_ang_mom))) then
          call messages_input_error(namespace, "MaxAngularMomentum", "Wrong maximum angular momentum for element"//trim(this_label))
        end if
        do jj = 1, this%ions%nspecies
          if (trim(adjustl(this_label)) == trim(adjustl(this%labels(jj)))) then
            max_ang_mom(jj) = trim(adjustl(this_max_ang_mom))
          end if
        end do
      end do
    end if
    call parse_block_end(blk)

    !%Variable SlakoDir
    !%Type string
    !%Default "./"
    !%Section Execution::IO
    !%Description
    !% Folder containing the Slako files
    !%End
    call parse_variable(namespace, 'SlakoDir', './', slako_dir)


    ! Dynamics variables

    call ion_dynamics_init(this%ions_dyn, namespace, this%ions)

    call parse_variable(namespace, 'TDDynamics', BO, this%dynamics)
    call messages_print_var_option('TDDynamics', this%dynamics, namespace=namespace)
    if (this%dynamics == BO) then
      call ion_dynamics_unfreeze(this%ions_dyn)
    end if

    !%Variable InitialIonicTemperature
    !%Type float
    !%Default 0.0
    !%Section DFTBPlusInterface
    !%Description
    !% If this variable is present, the ions will have initial velocities
    !% velocities to the atoms following a Boltzmann distribution with
    !% this temperature (in Kelvin). Used only if <tt>TDDynamics = Ehrenfest</tt>
    !% and  <tt>MoveIons = yes</tt>.
    !%End
    call parse_variable(namespace, 'InitialIonicTemperature', M_zero, initial_temp, unit = unit_kelvin)

    allocate(this%lasers)
    this%lasers%no_lasers = 0
    if (parse_block(namespace, 'TDExternalFields', blk) == 0) then
      this%laser_field = .true.
      ! No call to safe_deallocate macro here, as it gives an ICE with gfortran
      this%lasers%no_lasers = parse_block_n(blk)
      SAFE_ALLOCATE(this%lasers%lasers(1:this%lasers%no_lasers))

      do il = 1, this%lasers%no_lasers

        call parse_block_cmplx(blk, il-1, 0, pol(1))
        call parse_block_cmplx(blk, il-1, 1, pol(2))
        call parse_block_cmplx(blk, il-1, 2, pol(3))
        call parse_block_float(blk, il-1, 3, omega0)
        omega0 = units_to_atomic(units_inp%energy, omega0)

        call laser_set_frequency(this%lasers%lasers(il), omega0)
        pol(1:3) = pol(1:3)/norm2(abs(pol(1:3)))
        call laser_set_polarization(this%lasers%lasers(il), pol)

        call parse_block_string(blk, il-1, 4, envelope_expression)
        call tdf_read(ff, namespace, trim(envelope_expression), ierr)
        call laser_set_f(this%lasers%lasers(il), ff)

        ! Check if there is a phase.
        if (parse_block_cols(blk, il-1) > 5) then
          call parse_block_string(blk, il-1, 5, phase_expression)
          call tdf_read(phi, namespace, trim(phase_expression), ierr)
          call laser_set_phi(this%lasers%lasers(il), phi)
          if (ierr /= 0) then
            write(message(1),'(3A)') 'Error in the "', trim(envelope_expression), '" field defined in the TDExternalFields block:'
            write(message(2),'(3A)') 'Time-dependent phase function "', trim(phase_expression), '" not found.'
            call messages_warning(2, namespace=namespace)
          end if
        else
          call laser_set_empty_phi(this%lasers%lasers(il))
        end if
      end do

      call parse_block_end(blk)
    else
      this%laser_field = .false.
    end if

#ifdef HAVE_DFTBPLUS
    call TDftbPlus_init(this%dftbp, mpicomm=mpi_world%comm)

    call this%dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call setChild(pRoot, "Geometry", pGeo)
    call setChildValue(pGeo, "Periodic", .false.)
    call setChildValue(pGeo, "TypeNames", this%labels(1:this%ions%nspecies))
    call setChildValue(pGeo, "TypesAndCoordinates", reshape(this%species, [1, size(this%species)]), this%coords)
    call setChild(pRoot, "Hamiltonian", pHam)
    call setChild(pHam, "Dftb", pDftb)
    call setChildValue(pDftb, "Scc", .true.)

    !%Variable SccTolerance
    !%Type float
    !%Section DFTBPlusInterface
    !%Description
    !% Self-consistent-charges convergence tolerance. Once this
    !% tolerance has been achieved the SCC cycle will stop.
    !%End
    call parse_variable(namespace, 'SccTolerance', CNST(1e-9), this%scc_tolerance)
    call messages_print_var_value('SccTolerance', this%scc_tolerance, namespace=namespace)
    call setChildValue(pDftb, "SccTolerance", this%scc_tolerance)

    ! sub-block inside hamiltonian for the maximum angular momenta
    call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
    ! explicitly set the maximum angular momenta for the species
    do ii = 1, this%ions%nspecies
      call setChildValue(pMaxAng, this%labels(ii), max_ang_mom(ii))
    end do

    ! get the SK data
    ! You should provide the skfiles as found in the external/slakos/origin/mio-1-1/ folder. These can
    ! be downloaded with the utils/get_opt_externals script
    call setChild(pDftb, "SlaterKosterFiles", pSlakos)
    call setChild(pSlakos, "Type2FileNames", pType2Files)
    call setChildValue(pType2Files, "Prefix", slako_dir)
    call setChildValue(pType2Files, "Separator", "-")
    call setChildValue(pType2Files, "Suffix", ".skf")

    !  set up analysis options
    call setChild(pRoot, "Analysis", pAnalysis)
    call setChildValue(pAnalysis, "CalculateForces", .true.)

    call setChild(pRoot, "ParserOptions", pParserOpts)
    call setChildValue(pParserOpts, "ParserVersion", 5)
#endif

    if (this%dynamics == EHRENFEST) then
#ifdef HAVE_DFTBPLUS_DEVEL
      call setChild(pRoot, "ElectronDynamics", pElecDyn)
      call setChildValue(pElecDyn, "IonDynamics", ion_dynamics_ions_move(this%ions_dyn))
      if (ion_dynamics_ions_move(this%ions_dyn)) then
        call setChildValue(pElecDyn, "InitialTemperature", initial_temp)
      end if

      ! initialize with wrong arguments for the moment, will be overriden later
      call setChildValue(pElecDyn, "Steps", 1)
      call setChildValue(pElecDyn, "TimeStep", M_ONE)
      call setChild(pElecDyn, "Perturbation", pPerturb)
      if (this%laser_field) then
        call setChild(pPerturb, "Laser", pLaser)
        call setChildValue(pLaser, "PolarizationDirection", [ M_ONE , M_ZERO , M_ZERO ])
        call setChildValue(pLaser, "LaserEnergy", M_ONE)
        call setChildValue(pElecDyn, "FieldStrength", M_ONE)
      else
        call setChild(pPerturb, "None", pLaser)
      end if
#else
      message(1) = "DFTB Ehrenfest dynamics enabled only in DFTB development library"
      call messages_fatal(1, namespace=namespace)
#endif
    end if

#ifdef HAVE_DFTBPLUS
    message(1) = 'Input tree in HSD format:'
    call messages_info(1, namespace=namespace)
    call dumpHsd(input%hsdTree, stdout)

    ! initialise the DFTB+ calculator
    call this%dftbp%setupCalculator(input)
    call this%dftbp%setGeometry(this%coords)
#endif

    POP_SUB(dftb_init)
  end subroutine dftb_init

  ! ---------------------------------------------------------
  subroutine dftb_init_interaction(this, interaction)
    class(dftb_t),        target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(dftb_init_interaction)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by DFTB+."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(dftb_init_interaction)
  end subroutine dftb_init_interaction

  ! ---------------------------------------------------------
  subroutine dftb_initial_conditions(this)
    class(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_initial_conditions)

    select case (this%dynamics)
    case (BO)
#ifdef HAVE_DFTBPLUS
      call this%dftbp%getGradients(this%gradients)
      this%tot_force = -this%gradients
#endif
    case (EHRENFEST)
#ifdef HAVE_DFTBPLUS_DEVEL
      call this%dftbp%getEnergy(this%energy)
      call this%dftbp%initializeTimeProp(this%algo%dt, this%laser_field)
#endif
    end select

    POP_SUB(dftb_initial_conditions)
  end subroutine dftb_initial_conditions

  ! ---------------------------------------------------------
  logical function dftb_do_algorithmic_operation(this, operation) result(done)
    class(dftb_t),                  intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    integer :: ii, jj, il
    type(tdf_t) :: ff, phi
    CMPLX :: amp, pol(this%space%dim)
    FLOAT :: time, omega

    PUSH_SUB(dftb_do_algorithmic_operation)

    done = .true.
    select case (this%dynamics)
    case (BO)
      ! Born-Oppenheimer dynamics
      select case (operation%id)
      case (STORE_CURRENT_STATUS)
        ! Do nothing

      case (VERLET_START)
        SAFE_ALLOCATE(this%prev_acc(1:this%space%dim, 1:this%n_atom, 1))
        do jj = 1, this%n_atom
          this%acc(1:this%space%dim, jj) = this%tot_force(1:this%space%dim, jj) / this%mass(jj)
        end do

      case (VERLET_FINISH)
        SAFE_DEALLOCATE_A(this%prev_acc)

      case (VERLET_UPDATE_POS)
        do jj = 1, this%n_atom
          this%coords(1:this%space%dim, jj) = this%coords(1:this%space%dim, jj) + this%algo%dt * this%vel(1:this%space%dim, jj) &
            + M_HALF * this%algo%dt**2 * this%acc(1:this%space%dim, jj)
        end do
        this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", &
          QUANTITY_LABEL(POSITION), this%quantities(POSITION)%clock, "tick"))

      case (VERLET_COMPUTE_ACC)
        do ii = size(this%prev_acc, dim=3) - 1, 1, -1
          this%prev_acc(1:this%space%dim, 1:this%n_atom, ii + 1) = this%prev_acc(1:this%space%dim, 1:this%n_atom, ii)
        end do
        this%prev_acc(1:this%space%dim, 1:this%n_atom, 1) = this%acc(1:this%space%dim, 1:this%n_atom)
#ifdef HAVE_DFTBPLUS
        call this%dftbp%setGeometry(this%coords)
        call this%dftbp%getGradients(this%gradients)
        this%tot_force = -this%gradients
#endif
        do jj = 1, this%n_atom
          this%acc(1:this%space%dim, jj) = this%tot_force(1:this%space%dim, jj) / this%mass(jj)
        end do

      case (VERLET_COMPUTE_VEL)
        this%vel(1:this%space%dim, 1:this%n_atom) = this%vel(1:this%space%dim, 1:this%n_atom) &
          + M_HALF * this%algo%dt * (this%prev_acc(1:this%space%dim, 1:this%n_atom, 1) + &
          this%acc(1:this%space%dim, 1:this%n_atom))
        this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t("quantity", &
          QUANTITY_LABEL(VELOCITY), this%quantities(VELOCITY)%clock, "tick"))

      case default
        done = .false.
      end select

    case (EHRENFEST)
      ! Ehrenfest dynamics
      select case (operation%id)
      case (STORE_CURRENT_STATUS)
        ! Do nothing
      case (VERLET_START)
        !Do nothing
      case (VERLET_FINISH)
        !Do nothing
      case (VERLET_UPDATE_POS)
        this%field = M_zero
        time = this%clock%time()
        do il = 1, this%lasers%no_lasers
          ! get properties of laser
          call laser_get_f(this%lasers%lasers(il), ff)
          call laser_get_phi(this%lasers%lasers(il), phi)
          omega = laser_carrier_frequency(this%lasers%lasers(il))
          pol = laser_polarization(this%lasers%lasers(il))
          ! calculate electric field from laser
          amp = tdf(ff, time) * exp(M_zI * (omega*time + tdf(phi, time)))
          this%field(1:3) = this%field(1:3) + TOFLOAT(amp*pol(1:3))
        end do
#ifdef HAVE_DFTBPLUS_DEVEL
        call this%dftbp%setTdElectricField(this%field)
        call this%dftbp%doOneTdStep(this%clock%get_tick(), atomNetCharges=this%atom_charges, coord=this%coords,&
          force=this%tot_force, energy=this%energy)
#endif
      case (VERLET_COMPUTE_ACC)
        !Do nothing
      case (VERLET_COMPUTE_VEL)
        !Do nothing
      case default
        done = .false.
      end select

    end select

    POP_SUB(dftb_do_algorithmic_operation)
  end function dftb_do_algorithmic_operation

  ! ---------------------------------------------------------
  logical function dftb_is_tolerance_reached(this, tol) result(converged)
    class(dftb_t),        intent(in)    :: this
    FLOAT,                intent(in)    :: tol

    PUSH_SUB(dftb_is_tolerance_reached)

    ! this routine is never called at present, no reason to be here
    ASSERT(.false.)
    converged = .false.

    POP_SUB(dftb_is_tolerance_reached)
  end function dftb_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine dftb_output_start(this)
    class(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_output_start)

    ! Create output handle
    call io_mkdir('td.general', this%namespace)
    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_init(this%output_handle(OUTPUT_COORDINATES), 0, this%algo%dt, &
        trim(io_workpath("td.general/coordinates", this%namespace)))
      call write_iter_init(this%output_handle(OUTPUT_FORCES), 0, this%algo%dt, &
        trim(io_workpath("td.general/forces", this%namespace)))
    end if

    ! Output info for first iteration
    call this%output_write()

    POP_SUB(dftb_output_start)
  end subroutine dftb_output_start

  ! ---------------------------------------------------------
  subroutine dftb_output_finish(this)
    class(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_output_finish)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle(OUTPUT_COORDINATES))
      call write_iter_end(this%output_handle(OUTPUT_FORCES))
    end if

    POP_SUB(dftb_output_finish)
  end subroutine dftb_output_finish

  ! ---------------------------------------------------------
  subroutine dftb_output_write(this)
    class(dftb_t), intent(inout) :: this

    integer :: idir, iat, iout
    character(len=50) :: aux
    character(1) :: out_label(2)
    FLOAT :: tmp(this%space%dim)

    if (.not. mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(dftb_output_write)

    out_label(1) = "x"
    out_label(2) = "f"

    if (this%clock%get_tick() == 0) then
      ! header
      do iout = 1, 2
        call write_iter_clear(this%output_handle(iout))
        call write_iter_string(this%output_handle(iout),'#####################################################################')
        call write_iter_nl(this%output_handle(iout))
        call write_iter_string(this%output_handle(iout),'# HEADER')
        call write_iter_nl(this%output_handle(iout))

        ! first line: column names
        call write_iter_header_start(this%output_handle(iout))

        do iat = 1, this%n_atom
          do idir = 1, this%space%dim
            write(aux, '(a1,a1,i3,a1,i3,a1)') out_label(iout),'(', iat, ',', idir, ')'
            call write_iter_header(this%output_handle(iout), aux)
          end do
        end do
        call write_iter_nl(this%output_handle(iout))

        ! second line: units
        call write_iter_string(this%output_handle(iout), '#[Iter n.]')
        call write_iter_header(this%output_handle(iout), '[' // trim(units_abbrev(units_out%time)) // ']')
      end do

      do iat = 1, this%n_atom
        do idir = 1, this%space%dim
          call write_iter_header(this%output_handle(OUTPUT_COORDINATES), '[' // trim(units_abbrev(units_out%length)) // ']')
          call write_iter_header(this%output_handle(OUTPUT_FORCES), '[' // trim(units_abbrev(units_out%force)) // ']')
        end do
      end do

      do iout = 1, 2
        call write_iter_nl(this%output_handle(iout))
        call write_iter_string(this%output_handle(iout),'#######################################################################')
        call write_iter_nl(this%output_handle(iout))
      end do
    end if

    call write_iter_start(this%output_handle(OUTPUT_COORDINATES))
    call write_iter_start(this%output_handle(OUTPUT_FORCES))

    do iat = 1, this%n_atom
      ! Position
      tmp(1:this%space%dim) = units_from_atomic(units_out%length, this%coords(1:this%space%dim, iat))
      call write_iter_double(this%output_handle(OUTPUT_COORDINATES), tmp, this%space%dim)
      ! Force
      tmp(1:this%space%dim) = units_from_atomic(units_out%force, this%tot_force(1:this%space%dim, iat))
      call write_iter_double(this%output_handle(OUTPUT_FORCES), tmp, this%space%dim)
    end do

    call write_iter_nl(this%output_handle(OUTPUT_COORDINATES))
    call write_iter_nl(this%output_handle(OUTPUT_FORCES))

    POP_SUB(dftb_output_write)
  end subroutine dftb_output_write

  ! ---------------------------------------------------------
  subroutine dftb_update_quantity(this, iq)
    class(dftb_t),        intent(inout) :: this
    integer,              intent(in)    :: iq

    PUSH_SUB(dftb_update_quantity)

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(this%quantities(iq)%updated_on_demand)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(dftb_update_quantity)
  end subroutine dftb_update_quantity

  ! ---------------------------------------------------------
  subroutine dftb_update_exposed_quantity(partner, iq)
    class(dftb_t),        intent(inout) :: partner
    integer,              intent(in)    :: iq

    PUSH_SUB(dftb_update_exposed_quantity)

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(partner%quantities(iq)%updated_on_demand)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(dftb_update_exposed_quantity)
  end subroutine dftb_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine dftb_init_interaction_as_partner(partner, interaction)
    class(dftb_t),             intent(in)    :: partner
    class(interaction_t),      intent(inout) :: interaction

    PUSH_SUB(dftb_init_interaction_as_partner)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(dftb_init_interaction_as_partner)
  end subroutine dftb_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine dftb_copy_quantities_to_interaction(partner, interaction)
    class(dftb_t),          intent(inout) :: partner
    class(interaction_t),   intent(inout) :: interaction

    PUSH_SUB(dftb_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(dftb_copy_quantities_to_interaction)
  end subroutine dftb_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine dftb_update_interactions_start(this)
    class(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_update_interactions_start)

    ! Store previous force, as it is used as SCF criterium
    this%prev_tot_force(1:this%space%dim, 1:this%n_atom) = this%tot_force(1:this%space%dim, 1:this%n_atom)

    POP_SUB(dftb_update_interactions_start)
  end subroutine dftb_update_interactions_start

  ! ---------------------------------------------------------
  subroutine dftb_update_kinetic_energy(this)
    class(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_update_kinetic_energy)

    this%kinetic_energy = M_ZERO

    POP_SUB(dftb_update_kinetic_energy)
  end subroutine dftb_update_kinetic_energy

  ! ---------------------------------------------------------
  subroutine dftb_update_interactions_finish(this)
    class(dftb_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(dftb_update_interactions_finish)

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class default
        message(1) = "Interactions not implemented for DFTB+ systems."
        call messages_fatal(1, namespace=this%namespace)
      end select
    end do

    POP_SUB(dftb_update_interactions_finish)
  end subroutine dftb_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine dftb_restart_write_data(this)
    class(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_restart_write_data)

    message(1) = "DFTB system "//trim(this%namespace%get())//" cannot write restart data."
    call messages_warning(1, namespace=this%namespace)

    POP_SUB(dftb_restart_write_data)
  end subroutine dftb_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function dftb_restart_read_data(this)
    class(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_restart_read_data)

    ! restarting not yet supported
    dftb_restart_read_data = .false.

    POP_SUB(dftb_restart_read_data)
  end function dftb_restart_read_data

  ! ---------------------------------------------------------
  subroutine dftb_finalize(this)
    type(dftb_t), intent(inout) :: this

    PUSH_SUB(dftb_finalize)

    SAFE_DEALLOCATE_A(this%coords)
    SAFE_DEALLOCATE_A(this%acc)
    SAFE_DEALLOCATE_A(this%vel)
    SAFE_DEALLOCATE_A(this%tot_force)
    SAFE_DEALLOCATE_A(this%prev_tot_force)
    SAFE_DEALLOCATE_A(this%gradients)
    SAFE_DEALLOCATE_A(this%species)
    SAFE_DEALLOCATE_A(this%mass)
    SAFE_DEALLOCATE_P(this%ions)
    call ion_dynamics_end(this%ions_dyn)

    ! No call to safe_deallocate macro here, as it gives an ICE with gfortran
    if (associated(this%lasers)) then
      deallocate(this%lasers)
    end if

    call system_end(this)

    POP_SUB(dftb_finalize)
  end subroutine dftb_finalize

end module dftb_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
