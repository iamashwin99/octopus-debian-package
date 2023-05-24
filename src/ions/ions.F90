!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 M. Oliveira
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

module ions_oct_m
  use atom_oct_m
  use classical_particles_oct_m
  use iso_c_binding
  use debug_oct_m
  use distributed_oct_m
  use global_oct_m
  use interaction_oct_m
  use io_binary_oct_m
  use io_oct_m
  use ion_interaction_oct_m
  use lattice_vectors_oct_m
  use math_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use read_coords_oct_m
  use space_oct_m
  use species_oct_m
  use symmetries_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public :: ions_t

  type, extends(classical_particles_t) :: ions_t
    ! Components are public by default

    type(lattice_vectors_t) :: latt

    integer                   :: natoms
    type(atom_t), allocatable :: atom(:)

    integer                             :: ncatoms  !< For QM+MM calculations
    type(atom_classical_t), allocatable :: catom(:)

    type(symmetries_t) :: symm

    type(distributed_t) :: atoms_dist

    !> Information about the species
    integer                      :: nspecies
    type(species_t), allocatable :: species(:)
    logical                      :: only_user_def          !< Do we want to treat only user-defined species?
    logical,         private     :: species_time_dependent !< For time-dependent user defined species

    logical                 :: force_total_enforce
    type(ion_interaction_t) :: ion_interaction

    !> variables for external forces over the ions
    logical,     private :: apply_global_force
    type(tdf_t), private :: global_force_function
  contains
    procedure :: copy => ions_copy
    generic   :: assignment(=) => copy
    procedure :: partition => ions_partition
    procedure :: init_interaction => ions_init_interaction
    procedure :: initial_conditions => ions_initial_conditions
    procedure :: update_quantity => ions_update_quantity
    procedure :: update_exposed_quantity => ions_update_exposed_quantity
    procedure :: init_interaction_as_partner => ions_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => ions_copy_quantities_to_interaction
    procedure :: fold_atoms_into_cell => ions_fold_atoms_into_cell
    procedure :: min_distance => ions_min_distance
    procedure :: has_time_dependent_species => ions_has_time_dependent_species
    procedure :: val_charge => ions_val_charge
    procedure :: dipole => ions_dipole
    procedure :: translate => ions_translate
    procedure :: rotate => ions_rotate
    procedure :: global_force => ions_global_force
    procedure :: write_xyz => ions_write_xyz
    procedure :: read_xyz => ions_read_xyz
    procedure :: write_crystal => ions_write_crystal
    procedure :: write_bild_forces_file => ions_write_bild_forces_file
    procedure :: write_vtk_geometry => ions_write_vtk_geometry
    final :: ions_finalize
  end type ions_t

  interface ions_t
    procedure ions_constructor
  end interface ions_t

contains

  ! ---------------------------------------------------------
  function ions_constructor(namespace, print_info, latt_inp) result(ions)
    type(namespace_t),                 intent(in)    :: namespace
    logical,                 optional, intent(in)    :: print_info
    type(lattice_vectors_t), optional, intent(out)   :: latt_inp
    class(ions_t), pointer :: ions

    type(read_coords_info) :: xyz
    integer :: ia, ierr, idir
    character(len=100)  :: function_name
    FLOAT :: mindist
    FLOAT, allocatable :: factor(:)
    integer, allocatable :: site_type(:)
    logical, allocatable :: spherical_site(:)
    FLOAT, parameter :: threshold = CNST(1e-5)

    PUSH_SUB(ions_constructor)

    SAFE_ALLOCATE(ions)

    ions%namespace = namespace

    call space_init(ions%space, namespace)

    call species_init_global(namespace)

    ! initialize geometry
    call read_coords_init(xyz)

    ! load positions of the atoms
    call read_coords_read('Coordinates', xyz, ions%space, namespace)

    if (xyz%n < 1) then
      message(1) = "Coordinates have not been defined."
      call messages_fatal(1, namespace=namespace)
    end if

    ! Initialize parent class
    call classical_particles_init(ions, xyz%n)

    ! copy information from xyz to ions
    ions%natoms = xyz%n
    SAFE_ALLOCATE(ions%atom(1:ions%natoms))
    do ia = 1, ions%natoms
      call atom_init(ions%atom(ia), xyz%atom(ia)%label)
      ions%pos(:,ia) = xyz%atom(ia)%x(1:ions%space%dim)
      if (bitand(xyz%flags, XYZ_FLAGS_MOVE) /= 0) then
        ions%fixed(ia) = .not. xyz%atom(ia)%move
      end if
    end do

    if (allocated(xyz%latvec)) then
      ! Build lattice vectors from the XSF input
      ions%latt = lattice_vectors_t(namespace, ions%space, xyz%latvec)
    else
      ! Build lattice vectors from input file
      ions%latt = lattice_vectors_t(namespace, ions%space)
    end if

    ! Convert coordinates to Cartesian in case we have reduced coordinates
    if (xyz%source == READ_COORDS_REDUCED) then
      do ia = 1, ions%natoms
        ions%pos(:, ia) = ions%latt%red_to_cart(ions%pos(:, ia))
      end do
    end if

    call read_coords_end(xyz)

    ! load positions of the classical atoms, if any
    call read_coords_init(xyz)
    ions%ncatoms = 0
    call read_coords_read('Classical', xyz, ions%space, namespace)
    if (xyz%source /= READ_COORDS_ERR) then ! found classical atoms
      if (.not. bitand(xyz%flags, XYZ_FLAGS_CHARGE) /= 0) then
        message(1) = "Need to know charge for the classical atoms."
        message(2) = "Please use a .pdb"
        call messages_fatal(2, namespace=namespace)
      end if
      ions%ncatoms = xyz%n
      write(message(1), '(a,i8)') 'Info: Number of classical atoms = ', ions%ncatoms
      call messages_info(1, namespace=namespace)
      if (ions%ncatoms>0) then
        SAFE_ALLOCATE(ions%catom(1:ions%ncatoms))
        do ia = 1, ions%ncatoms
          call atom_classical_init(ions%catom(ia), xyz%atom(ia)%label, xyz%atom(ia)%x, xyz%atom(ia)%charge)
        end do
      end if
      call read_coords_end(xyz)
    end if


    call ions_fold_atoms_into_cell(ions)
    call ions_init_species(ions, print_info=print_info)
    call distributed_nullify(ions%atoms_dist, ions%natoms)

    if (present(latt_inp)) then
      ! The lattice as read from the input might be needed by some other part of the code, so we save it
      latt_inp = ions%latt
    end if

    ! Now that we have processed the atomic coordinates, we renormalize the
    ! lattice parameters along the non-periodic dimensions
    if (ions%space%has_mixed_periodicity()) then
      SAFE_ALLOCATE(factor(ions%space%dim))
      do idir = 1, ions%space%periodic_dim
        factor(idir) = M_ONE
      end do
      do idir = ions%space%periodic_dim + 1, ions%space%dim
        factor(idir) = M_ONE/norm2(ions%latt%rlattice(1:ions%space%dim, idir))
      end do
      call ions%latt%scale(factor)
      SAFE_DEALLOCATE_A(factor)
    end if

    ! Set the masses. This needs to be done after initializing the species.
    do ia = 1, ions%natoms
      ions%mass(ia) = species_mass(ions%atom(ia)%species)
    end do

    ! Check that atoms are not too close
    if (ions%natoms > 1) then
      mindist = ions_min_distance(ions, real_atoms_only = .false.)
      if (mindist < threshold) then
        write(message(1), '(a)') "Some of the atoms seem to sit too close to each other."
        write(message(2), '(a)') "Please review your input files and the output geometry (in 'static/')."
        write(message(3), '(a, f12.6, 1x, a)') "Minimum distance = ", &
          units_from_atomic(units_out%length, mindist), trim(units_abbrev(units_out%length))
        call messages_warning(3, namespace=namespace)

        ! then write out the geometry, whether asked for or not in Output variable
        call io_mkdir(STATIC_DIR, namespace)
        call ions%write_xyz(trim(STATIC_DIR)//'/geometry')
      end if

      if (ions_min_distance(ions, real_atoms_only = .true.) < threshold) then
        message(1) = "It cannot be correct to run with physical atoms so close."
        call messages_fatal(1, namespace=namespace)
      end if
    end if

    !Initialize symmetries
    SAFE_ALLOCATE(spherical_site(1:ions%natoms))
    SAFE_ALLOCATE(site_type(1:ions%natoms))
    do ia = 1, ions%natoms
      spherical_site(ia) = .not. (&
        species_type(ions%atom(ia)%species) == SPECIES_USDEF          .or. &
        species_type(ions%atom(ia)%species) == SPECIES_JELLIUM_SLAB   .or. &
        species_type(ions%atom(ia)%species) == SPECIES_CHARGE_DENSITY .or. &
        species_type(ions%atom(ia)%species) == SPECIES_FROM_FILE)

      site_type(ia) = species_index(ions%atom(ia)%species)
    end do

    ions%symm = symmetries_t(ions%namespace, ions%space, ions%latt, ions%natoms, ions%pos, site_type, spherical_site)

    SAFE_DEALLOCATE_A(spherical_site)
    SAFE_DEALLOCATE_A(site_type)


    call ion_interaction_init(ions%ion_interaction, namespace, ions%space, ions%natoms)

    !%Variable ForceTotalEnforce
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% (Experimental) If this variable is set to "yes", then the sum
    !% of the total forces will be enforced to be zero.
    !%End
    call parse_variable(namespace, 'ForceTotalEnforce', .false., ions%force_total_enforce)
    if (ions%force_total_enforce) call messages_experimental('ForceTotalEnforce', namespace=namespace)

    !%Variable TDGlobalForce
    !%Type string
    !%Section Time-Dependent
    !%Description
    !% If this variable is set, a global time-dependent force will be
    !% applied to the ions in the x direction during a time-dependent
    !% run. This variable defines the base name of the force, that
    !% should be defined in the <tt>TDFunctions</tt> block. This force
    !% does not affect the electrons directly.
    !%End

    if (parse_is_defined(namespace, 'TDGlobalForce')) then

      ions%apply_global_force = .true.

      call parse_variable(namespace, 'TDGlobalForce', 'none', function_name)
      call tdf_read(ions%global_force_function, namespace, trim(function_name), ierr)

      if (ierr /= 0) then
        call messages_write("You have enabled the GlobalForce option but Octopus could not find")
        call messages_write("the '"//trim(function_name)//"' function in the TDFunctions block.")
        call messages_fatal(namespace=namespace)
      end if

    else

      ions%apply_global_force = .false.

    end if

    POP_SUB(ions_constructor)
  end function ions_constructor

  ! ---------------------------------------------------------
  subroutine ions_init_species(ions, print_info)
    type(ions_t),      intent(inout) :: ions
    logical, optional, intent(in)    :: print_info

    logical :: print_info_, spec_user_defined
    integer :: i, j, k, ispin

    PUSH_SUB(ions_init_species)

    print_info_ = .true.
    if (present(print_info)) then
      print_info_ = print_info
    end if
    ! First, count the species
    ions%nspecies = 0
    atoms1:  do i = 1, ions%natoms
      do j = 1, i - 1
        if (atom_same_species(ions%atom(j), ions%atom(i))) cycle atoms1
      end do
      ions%nspecies = ions%nspecies + 1
    end do atoms1

    ! Allocate the species structure.
    SAFE_ALLOCATE(ions%species(1:ions%nspecies))

    ! Now, read the data.
    k = 0
    ions%only_user_def = .true.
    atoms2: do i = 1, ions%natoms
      do j = 1, i - 1
        if (atom_same_species(ions%atom(j), ions%atom(i))) cycle atoms2
      end do
      k = k + 1
      call species_init(ions%species(k), atom_get_label(ions%atom(j)), k)
      call species_read(ions%species(k), ions%namespace)
      ! these are the species which do not represent atoms
      ions%only_user_def = ions%only_user_def .and. .not. species_represents_real_atom(ions%species(k))

      if (species_is_ps(ions%species(k)) .and. ions%space%dim /= 3) then
        message(1) = "Pseudopotentials may only be used with Dimensions = 3."
        call messages_fatal(1, namespace=ions%namespace)
      end if

      if (species_type(ions%species(k)) == SPECIES_JELLIUM_SLAB) then
        if (ions%space%is_periodic() .and. ions%space%periodic_dim /= 2) then
          message(1) = "Periodic jelium slab can only be used if PeriodicDim = 2"
          call messages_fatal(1, namespace=ions%namespace)
        end if
      end if

    end do atoms2

    ! Reads the spin components. This is read here, as well as in states_init,
    ! to be able to pass it to the pseudopotential initializations subroutine.
    call parse_variable(ions%namespace, 'SpinComponents', 1, ispin)
    if (.not. varinfo_valid_option('SpinComponents', ispin)) call messages_input_error(ions%namespace, 'SpinComponents')
    ispin = min(2, ispin)

    if (print_info_) then
      call messages_print_stress(msg="Species", namespace=ions%namespace)
    end if
    do i = 1, ions%nspecies
      call species_build(ions%species(i), ions%namespace, ispin, ions%space%dim, print_info=print_info_)
    end do
    if (print_info_) then
      call messages_print_stress(namespace=ions%namespace)
    end if

    !%Variable SpeciesTimeDependent
    !%Type logical
    !%Default no
    !%Section System::Species
    !%Description
    !% When this variable is set, the potential defined in the block <tt>Species</tt> is calculated
    !% and applied to the Hamiltonian at each time step. You must have at least one <tt>species_user_defined</tt>
    !% type of species to use this.
    !%End
    call parse_variable(ions%namespace, 'SpeciesTimeDependent', .false., ions%species_time_dependent)
    ! we must have at least one user defined species in order to have time dependency
    do i = 1,ions%nspecies
      if (species_type(ions%species(i)) == SPECIES_USDEF) then
        spec_user_defined = .true.
      end if
    end do
    if (ions%species_time_dependent .and. .not. spec_user_defined) then
      call messages_input_error(ions%namespace, 'SpeciesTimeDependent')
    end if

    !  assign species
    do i = 1, ions%natoms
      do j = 1, ions%nspecies
        if (atom_same_species(ions%atom(i), ions%species(j))) then
          call atom_set_species(ions%atom(i), ions%species(j))
          exit
        end if
      end do
    end do

    POP_SUB(ions_init_species)
  end subroutine ions_init_species

  !--------------------------------------------------------------
  subroutine ions_copy(ions_out, ions_in)
    class(ions_t),     intent(out) :: ions_out
    class(ions_t),     intent(in)  :: ions_in

    PUSH_SUB(ions_copy)

    call classical_particles_copy(ions_out, ions_in)

    ions_out%latt = ions_in%latt

    ions_out%natoms = ions_in%natoms
    SAFE_ALLOCATE(ions_out%atom(1:ions_out%natoms))
    ions_out%atom = ions_in%atom

    ions_out%ncatoms = ions_in%ncatoms
    SAFE_ALLOCATE(ions_out%catom(1:ions_out%ncatoms))
    if (ions_in%ncatoms > 0) then
      ions_out%catom(1:ions_out%ncatoms) = ions_in%catom(1:ions_in%ncatoms)
    end if

    ions_out%nspecies = ions_in%nspecies
    SAFE_ALLOCATE(ions_out%species(1:ions_out%nspecies))
    ions_out%species = ions_in%species

    ions_out%only_user_def     = ions_in%only_user_def

    call distributed_copy(ions_in%atoms_dist, ions_out%atoms_dist)

    POP_SUB(ions_copy)
  end subroutine ions_copy

  ! ---------------------------------------------------------
  subroutine ions_partition(this, mc)
    class(ions_t),               intent(inout) :: this
    type(multicomm_t),           intent(in)    :: mc

    PUSH_SUB(ions_partition)

    call distributed_init(this%atoms_dist, this%natoms, mc%group_comm(P_STRATEGY_STATES), "atoms")

    call ion_interaction_init_parallelization(this%ion_interaction, this%natoms, mc)

    POP_SUB(ions_partition)
  end subroutine ions_partition

  ! ---------------------------------------------------------
  subroutine ions_init_interaction(this, interaction)
    class(ions_t),        target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(ions_init_interaction)

    select type (interaction)
    class default
      call classical_particles_init_interaction(this, interaction)
    end select

    POP_SUB(ions_init_interaction)
  end subroutine ions_init_interaction

  ! ---------------------------------------------------------
  subroutine ions_initial_conditions(this)
    class(ions_t), intent(inout) :: this

    PUSH_SUB(ions_initial_conditions)

    POP_SUB(ions_initial_conditions)
  end subroutine ions_initial_conditions

  ! ---------------------------------------------------------
  subroutine ions_update_quantity(this, iq)
    class(ions_t), intent(inout) :: this
    integer,       intent(in)    :: iq

    PUSH_SUB(ions_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_quantity(this, iq)
    end select

    POP_SUB(ions_update_quantity)
  end subroutine ions_update_quantity

  ! ---------------------------------------------------------
  subroutine ions_update_exposed_quantity(partner, iq)
    class(ions_t), intent(inout) :: partner
    integer,       intent(in)    :: iq

    PUSH_SUB(ions_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_exposed_quantity(partner, iq)
    end select

    POP_SUB(ions_update_exposed_quantity)
  end subroutine ions_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine ions_init_interaction_as_partner(partner, interaction)
    class(ions_t),                intent(in)    :: partner
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(ions_init_interaction_as_partner)

    select type (interaction)
    class default
      call classical_particles_init_interaction_as_partner(partner, interaction)
    end select

    POP_SUB(ions_init_interaction_as_partner)
  end subroutine ions_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine ions_copy_quantities_to_interaction(partner, interaction)
    class(ions_t),          intent(inout) :: partner
    class(interaction_t),   intent(inout) :: interaction

    PUSH_SUB(ions_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(ions_copy_quantities_to_interaction)
  end subroutine ions_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine ions_fold_atoms_into_cell(this)
    class(ions_t),       intent(inout) :: this

    integer :: iatom

    PUSH_SUB(ions_fold_atoms_into_cell)

    do iatom = 1, this%natoms
      this%pos(:, iatom) = this%latt%fold_into_cell(this%pos(:, iatom))
    end do

    POP_SUB(ions_fold_atoms_into_cell)
  end subroutine ions_fold_atoms_into_cell

  ! ---------------------------------------------------------
  FLOAT function ions_min_distance(this, real_atoms_only) result(rmin)
    class(ions_t),      intent(in) :: this
    logical,  optional, intent(in) :: real_atoms_only

    integer :: iatom, jatom, idir
    FLOAT   :: xx(this%space%dim)
    logical :: real_atoms_only_
    type(species_t), pointer :: species

    PUSH_SUB(ions_min_distance)

    real_atoms_only_ = optional_default(real_atoms_only, .false.)

    rmin = huge(rmin)
    do iatom = 1, this%natoms
      call atom_get_species(this%atom(iatom), species)
      if (real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
      do jatom = iatom + 1, this%natoms
        call atom_get_species(this%atom(iatom), species)
        if (real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
        xx = abs(this%pos(:, iatom) - this%pos(:, jatom))
        if (this%space%is_periodic()) then
          xx = this%latt%cart_to_red(xx)
          do idir = 1, this%space%periodic_dim
            xx(idir) = xx(idir) - floor(xx(idir) + M_HALF)
          end do
          xx = this%latt%red_to_cart(xx)
        end if
        rmin = min(norm2(xx), rmin)
      end do
    end do

    if (.not. (this%only_user_def .and. real_atoms_only_)) then
      ! what if the nearest neighbors are periodic images?
      do idir = 1, this%space%periodic_dim
        rmin = min(rmin, norm2(this%latt%rlattice(:,idir)))
      end do
    end if

    POP_SUB(ions_min_distance)
  end function ions_min_distance

  ! ---------------------------------------------------------
  logical function ions_has_time_dependent_species(this) result(time_dependent)
    class(ions_t),     intent(in) :: this

    PUSH_SUB(ions_has_time_dependent_species)

    time_dependent = this%species_time_dependent

    POP_SUB(ions_has_time_dependent_species)
  end function ions_has_time_dependent_species

  ! ---------------------------------------------------------
  FLOAT function ions_val_charge(this, mask) result(val_charge)
    class(ions_t),              intent(in) :: this
    logical,          optional, intent(in) :: mask(:)

    integer :: iatom

    PUSH_SUB(ions_val_charge)

    val_charge = M_ZERO
    do iatom = 1, this%natoms
      if (present(mask)) then
        if (.not. mask(iatom)) cycle
      end if
      val_charge = val_charge - species_zval(this%atom(iatom)%species)
    end do

    POP_SUB(ions_val_charge)
  end function ions_val_charge

  ! ---------------------------------------------------------
  function ions_dipole(this, mask) result(dipole)
    class(ions_t),               intent(in) :: this
    logical,           optional, intent(in) :: mask(:)
    FLOAT :: dipole(this%space%dim)

    integer :: ia

    PUSH_SUB(ions_dipole)

    dipole = M_ZERO
    do ia = 1, this%natoms
      if (present(mask)) then
        if (.not. mask(ia)) cycle
      end if
      dipole = dipole + species_zval(this%atom(ia)%species)*this%pos(:, ia)
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(ions_dipole)
  end function ions_dipole

  ! ---------------------------------------------------------
  subroutine ions_translate(this, xx)
    class(ions_t),     intent(inout) :: this
    FLOAT,             intent(in)    :: xx(this%space%dim)

    integer  :: iatom

    PUSH_SUB(ions_translate)

    do iatom = 1, this%natoms
      this%pos(:, iatom) = this%pos(:, iatom) - xx
    end do
    do iatom = 1, this%ncatoms
      this%catom(iatom)%x(1:this%space%dim) = this%catom(iatom)%x(1:this%space%dim) - xx
    end do

    POP_SUB(ions_translate)
  end subroutine ions_translate

  ! ---------------------------------------------------------
  subroutine ions_rotate(this, from, from2, to)
    class(ions_t),     intent(inout) :: this
    FLOAT,             intent(in)    :: from(this%space%dim)   !< assumed to be normalized
    FLOAT,             intent(in)    :: from2(this%space%dim)  !< assumed to be normalized
    FLOAT,             intent(in)    :: to(this%space%dim)     !< assumed to be normalized

    integer :: iatom
    FLOAT :: m1(3, 3), m2(3, 3)
    FLOAT :: m3(3, 3), f2(3), per(3)
    FLOAT :: alpha, r

    PUSH_SUB(ions_rotate)

    if (this%space%dim /= 3) then
      call messages_not_implemented("ions_rotate in other than 3 dimensions", namespace=this%namespace)
    end if

    ! initialize matrices
    m1 = diagonal_matrix(3, M_ONE)

    ! rotate the to-axis to the z-axis
    if (to(2) /= M_ZERO) then
      alpha = atan2(to(2), to(1))
      call rotate(m1, alpha, 3)
    end if
    alpha = atan2(norm2(to(1:2)), to(3))
    call rotate(m1, -alpha, 2)

    ! get perpendicular to z and from
    f2 = matmul(m1, from)
    per(1) = -f2(2)
    per(2) =  f2(1)
    per(3) = M_ZERO
    r = norm2(per)
    if (r > M_ZERO) then
      per = per/r
    else
      per(2) = M_ONE
    end if

    ! rotate perpendicular axis to the y-axis
    m2 = diagonal_matrix(3, M_ONE)
    alpha = atan2(per(1), per(2))
    call rotate(m2, -alpha, 3)

    ! rotate from => to (around the y-axis)
    m3 = diagonal_matrix(3, M_ONE)
    alpha = acos(sum(from*to))
    call rotate(m3, -alpha, 2)

    ! join matrices
    m2 = matmul(transpose(m2), matmul(m3, m2))

    ! rotate around the z-axis to get the second axis
    per = matmul(m2, matmul(m1, from2))
    alpha = atan2(per(1), per(2))
    call rotate(m2, -alpha, 3) ! second axis is now y

    ! get combined transformation
    m1 = matmul(transpose(m1), matmul(m2, m1))

    ! now transform the coordinates
    ! it is written in this way to avoid what I consider a bug in the Intel compiler
    do iatom = 1, this%natoms
      f2 = this%pos(:, iatom)
      this%pos(:, iatom) = matmul(m1, f2)
    end do

    do iatom = 1, this%ncatoms
      f2 = this%catom(iatom)%x(1:this%space%dim)
      this%catom(iatom)%x(1:this%space%dim) = matmul(m1, f2)
    end do

    POP_SUB(ions_rotate)
  contains

    ! ---------------------------------------------------------
    subroutine rotate(m, angle, dir)
      FLOAT,   intent(inout) :: m(3, 3)
      FLOAT,   intent(in)    :: angle
      integer, intent(in)    :: dir

      FLOAT :: aux(3, 3), ca, sa

      PUSH_SUB(ions_rotate.rotate)

      ca = cos(angle)
      sa = sin(angle)

      aux = M_ZERO
      select case (dir)
      case (1)
        aux(1, 1) = M_ONE
        aux(2, 2) = ca
        aux(3, 3) = ca
        aux(2, 3) = sa
        aux(3, 2) = -sa
      case (2)
        aux(2, 2) = M_ONE
        aux(1, 1) = ca
        aux(3, 3) = ca
        aux(1, 3) = sa
        aux(3, 1) = -sa
      case (3)
        aux(3, 3) = M_ONE
        aux(1, 1) = ca
        aux(2, 2) = ca
        aux(1, 2) = sa
        aux(2, 1) = -sa
      end select

      m = matmul(aux, m)

      POP_SUB(ions_rotate.rotate)
    end subroutine rotate

  end subroutine ions_rotate

  ! ---------------------------------------------------------
  function ions_global_force(this, time) result(force)
    class(ions_t),        intent(in)    :: this
    FLOAT,                intent(in)    :: time
    FLOAT :: force(this%space%dim)

    PUSH_SUB(ions_global_force)

    force = M_ZERO

    if (this%apply_global_force) then
      force(1) = units_to_atomic(units_inp%force, tdf(this%global_force_function, time))
    end if

    POP_SUB(ions_global_force)
  end function ions_global_force

  ! ---------------------------------------------------------
  subroutine ions_write_xyz(this, fname, append, comment)
    class(ions_t),              intent(in) :: this
    character(len=*),           intent(in) :: fname
    logical,          optional, intent(in) :: append
    character(len=*), optional, intent(in) :: comment

    integer :: iatom, idim, iunit
    character(len=6) position
    character(len=19) :: frmt

    if (.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(ions_write_xyz)

    position = 'asis'
    if (present(append)) then
      if (append) position = 'append'
    end if
    iunit = io_open(trim(fname)//'.xyz', this%namespace, action='write', position=position)

    write(iunit, '(i4)') this%natoms
    if (present(comment)) then
      write(iunit, '(1x,a)') comment
    else
      write(iunit, '(1x,a,a)') 'units: ', trim(units_abbrev(units_out%length_xyz_file))
    end if

    write(unit=frmt, fmt="(a5,i2.2,a4,i2.2,a6)") "(6x,a", LABEL_LEN, ",2x,", this%space%dim,"f12.6)"
    do iatom = 1, this%natoms
      write(unit=iunit, fmt=frmt) this%atom(iatom)%label, &
        (units_from_atomic(units_out%length_xyz_file, this%pos(idim, iatom)), idim=1, this%space%dim)
    end do
    call io_close(iunit)

    if (this%ncatoms > 0) then
      iunit = io_open(trim(fname)//'_classical.xyz', this%namespace, action='write', position=position)
      write(iunit, '(i4)') this%ncatoms
      write(iunit, '(1x)')
      do iatom = 1, this%ncatoms
        call atom_classical_write_xyz(this%catom(iatom), this%space%dim, iunit)
      end do
      call io_close(iunit)
    end if

    POP_SUB(ions_write_xyz)
  end subroutine ions_write_xyz

  ! ---------------------------------------------------------
  subroutine ions_read_xyz(this, fname, comment)
    class(ions_t),              intent(inout) :: this
    character(len=*),           intent(in)    :: fname
    character(len=*), optional, intent(in)    :: comment

    integer :: iatom, idir, iunit
    character(len=19) :: frmt, dum
    FLOAT :: tmp(this%space%dim)

    PUSH_SUB(ions_read_xyz)

    iunit = io_open(trim(fname)//'.xyz', this%namespace, action='read', position='rewind')

    read(iunit, '(i4)') this%natoms
    if (present(comment)) then
      read(iunit, *)
    else
      read(iunit, *)
    end if
    write(unit=frmt, fmt="(a5,i2.2,a4,i2.2,a6)") "(6x,a", LABEL_LEN, ",2x,", this%space%dim, "f12.6)"
    do iatom = 1, this%natoms
      read(unit=iunit, fmt=frmt) dum, (tmp(idir), idir=1, this%space%dim)

      this%pos(:, iatom) = units_to_atomic(units_out%length_xyz_file, tmp)
    end do
    call io_close(iunit)

    if (this%ncatoms > 0) then
      iunit = io_open(trim(fname)//'_classical.xyz', this%namespace, action='read', position='rewind')
      read(iunit, '(i4)') this%ncatoms
      read(iunit, *)
      do iatom = 1, this%ncatoms
        call atom_classical_read_xyz(this%catom(iatom), this%space%dim, iunit)
      end do
      call io_close(iunit)
    end if

    POP_SUB(ions_read_xyz)
  end subroutine ions_read_xyz

  ! ----------------------------------------------------------------
  !> This subroutine creates a crystal by replicating the geometry and
  !! writes the result to dir//'crystal.xyz'
  subroutine ions_write_crystal(this, dir)
    class(ions_t),           intent(in) :: this
    character(len=*),        intent(in) :: dir

    type(lattice_iterator_t) :: latt_iter
    FLOAT :: radius, pos(this%space%dim)
    integer :: iatom, icopy, iunit

    PUSH_SUB(ions_write_crystal)

    radius = maxval(M_HALF*norm2(this%latt%rlattice(:,1:this%space%periodic_dim), dim=1))*(M_ONE + M_EPSILON)
    latt_iter = lattice_iterator_t(this%latt, radius)

    if (mpi_grp_is_root(mpi_world)) then

      iunit = io_open(trim(dir)//'/crystal.xyz', this%namespace, action='write')

      write(iunit, '(i9)') this%natoms*latt_iter%n_cells
      write(iunit, '(a)') '#generated by Octopus'

      do iatom = 1, this%natoms
        do icopy = 1, latt_iter%n_cells
          pos = units_from_atomic(units_out%length, this%pos(:, iatom) + latt_iter%get(icopy))
          write(iunit, '(a, 99f12.6)') this%atom(iatom)%label, pos
        end do
      end do

      call io_close(iunit)
    end if

    POP_SUB(ions_write_crystal)
  end subroutine ions_write_crystal

  ! ---------------------------------------------------------
  subroutine ions_write_bild_forces_file(this, dir, fname)
    class(ions_t),      intent(in) :: this
    character(len=*),   intent(in) :: dir, fname

    integer :: iunit, iatom, idir
    FLOAT :: force(this%space%dim), center(this%space%dim)
    character(len=20) frmt

    if (.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(write_bild_forces_file)

    call io_mkdir(dir, this%namespace)
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.bild', this%namespace, action='write', &
      position='asis')

    write(frmt,'(a,i0,a)')'(a,2(', this%space%dim,'f16.6,1x))'

    write(iunit, '(a)')'.comment : force vectors in ['//trim(units_abbrev(units_out%force))//']'
    write(iunit, *)
    write(iunit, '(a)')'.color red'
    write(iunit, *)
    do iatom = 1, this%natoms
      center = units_from_atomic(units_out%length, this%pos(:, iatom))
      force = units_from_atomic(units_out%force, this%tot_force(:, iatom))
      write(iunit, '(a,1x,i4,1x,a2,1x,a6,1x,f10.6,a)') '.comment :', iatom, trim(this%atom(iatom)%label), &
        'force:', norm2(force), '['//trim(units_abbrev(units_out%force))//']'
      write(iunit,fmt=trim(frmt)) '.arrow', (center(idir), idir = 1, this%space%dim), &
        (center(idir) + force(idir), idir = 1, this%space%dim)
      write(iunit,*)
    end do

    call io_close(iunit)

    POP_SUB(ions_write_bild_forces_file)
  end subroutine ions_write_bild_forces_file

  ! -----------------------------------------------------
  subroutine ions_write_vtk_geometry(this, filename, ascii)
    class(ions_t),              intent(in) :: this
    character(len=*),           intent(in) :: filename
    logical,          optional, intent(in) :: ascii

    integer :: iunit, iatom, ierr
    logical :: ascii_
    FLOAT, allocatable :: data(:, :)
    integer, allocatable :: idata(:, :)
    character(len=MAX_PATH_LEN) :: fullname

    PUSH_SUB(ions_write_vtk_geometry)

    ASSERT(this%space%dim == 3)

    ascii_ = optional_default(ascii, .true.)

    fullname = trim(filename)//".vtk"

    iunit = io_open(trim(fullname), this%namespace, action='write')

    write(iunit, '(1a)') '# vtk DataFile Version 3.0 '
    write(iunit, '(6a)') 'Generated by octopus ', trim(conf%version), ' -  git: ', &
      trim(conf%git_commit), " build: ",  trim(conf%build_time)

    if (ascii_) then
      write(iunit, '(1a)') 'ASCII'
    else
      write(iunit, '(1a)') 'BINARY'
    end if

    write(iunit, '(1a)') 'DATASET POLYDATA'

    write(iunit, '(a,i9,a)') 'POINTS ', this%natoms, ' double'

    if (ascii_) then
      do iatom = 1, this%natoms
        write(iunit, '(3f12.6)') this%pos(1:3, iatom)
      end do
    else
      call io_close(iunit)
      SAFE_ALLOCATE(data(1:3, 1:this%natoms))
      do iatom = 1, this%natoms
        data(1:3, iatom) = this%pos(1:3, iatom)
      end do
      call io_binary_write(io_workpath(fullname, this%namespace), i4_to_i8(3*this%natoms), data, &
        ierr, nohead = .true., fendian = io_binary_is_little_endian())
      SAFE_DEALLOCATE_A(data)
      iunit = io_open(trim(fullname), this%namespace, action='write', position = 'append')
      write(iunit, '(1a)') ''
    end if

    write(iunit, '(a,2i9)') 'VERTICES ', this%natoms, 2*this%natoms

    if (ascii_) then
      do iatom = 1, this%natoms
        write(iunit, '(2i9)') 1, iatom - 1
      end do
    else
      call io_close(iunit)
      SAFE_ALLOCATE(idata(1:2, 1:this%natoms))
      do iatom = 1, this%natoms
        idata(1, iatom) = 1
        idata(2, iatom) = iatom - 1
      end do
      call io_binary_write(io_workpath(fullname, this%namespace), i4_to_i8(2*this%natoms), idata, &
        ierr, nohead = .true., fendian = io_binary_is_little_endian())
      SAFE_DEALLOCATE_A(idata)
      iunit = io_open(trim(fullname), this%namespace, action='write', position = 'append')
      write(iunit, '(1a)') ''
    end if

    write(iunit, '(a,i9)') 'POINT_DATA', this%natoms
    write(iunit, '(a)') 'SCALARS element integer'
    write(iunit, '(a)') 'LOOKUP_TABLE default'

    if (ascii_) then
      do iatom = 1, this%natoms
        write(iunit, '(i9)') nint(species_z(this%atom(iatom)%species))
      end do
    else
      call io_close(iunit)

      SAFE_ALLOCATE(idata(1:this%natoms, 1))

      do iatom = 1, this%natoms
        idata(iatom, 1) = nint(species_z(this%atom(iatom)%species))
      end do

      call io_binary_write(io_workpath(fullname, this%namespace), i4_to_i8(this%natoms), idata, &
        ierr, nohead = .true., fendian = io_binary_is_little_endian())

      SAFE_DEALLOCATE_A(idata)

      iunit = io_open(trim(fullname), this%namespace, action='write', position = 'append')
      write(iunit, '(1a)') ''
    end if

    call io_close(iunit)

    POP_SUB(ions_write_vtk_geometry)
  end subroutine ions_write_vtk_geometry

  ! ---------------------------------------------------------
  subroutine ions_finalize(ions)
    type(ions_t),     intent(inout) :: ions

    PUSH_SUB(ions_finalize)

    call classical_particles_end(ions)

    call distributed_end(ions%atoms_dist)

    call ion_interaction_end(ions%ion_interaction)

    SAFE_DEALLOCATE_A(ions%atom)
    ions%natoms=0
    SAFE_DEALLOCATE_A(ions%catom)
    ions%ncatoms=0

    call species_end(ions%nspecies, ions%species)
    SAFE_DEALLOCATE_A(ions%species)
    ions%nspecies=0

    call species_end_global()

    POP_SUB(ions_finalize)
  end subroutine ions_finalize

end module ions_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
