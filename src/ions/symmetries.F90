!! Copyright (C) 2009 X. Andrade
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

module symmetries_oct_m
  use debug_oct_m
  use global_oct_m
  use iso_c_binding
  use lattice_vectors_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use string_oct_m
  use space_oct_m
  use spglib_f08
  use symm_op_oct_m

  implicit none

  private

  public ::                        &
    symmetries_t,                  &
    symmetries_number,             &
    symmetries_apply_kpoint_red,   &
    symmetries_space_group_number, &
    symmetries_have_break_dir,     &
    symmetries_identity_index,     &
    symmetries_write_info

  type symmetries_t
    private
    type(symm_op_t), allocatable, public :: ops(:)
    integer, public          :: nops = 0
    FLOAT                    :: breakdir(1:3)
    integer, public          :: periodic_dim
    integer                  :: space_group
    logical                  :: any_non_spherical
    logical                  :: symmetries_compute
    character(len=6)         :: group_name
    character(len=30)        :: group_elements
    character(len=11)        :: symbol
    character(len=7)         :: schoenflies
  contains
    procedure :: copy => symmetries_copy
    generic   :: assignment(=) => copy
    final :: symmetries_finalizer
  end type symmetries_t

  FLOAT, public :: symprec = CNST(1e-5)

  !> NOTE: unfortunately, these routines use global variables shared among them
  interface
    subroutine symmetries_finite_init(natoms, types, positions, verbosity, point_group)
      integer, intent(in)  :: natoms
      integer, intent(in)  :: types !< (natoms)
      FLOAT,   intent(in)  :: positions !< (3, natoms)
      integer, intent(in)  :: verbosity
      integer, intent(out) :: point_group
    end subroutine symmetries_finite_init

    subroutine symmetries_finite_get_group_name(point_group, name)
      integer,          intent(in)  :: point_group
      character(len=*), intent(out) :: name
    end subroutine symmetries_finite_get_group_name

    subroutine symmetries_finite_get_group_elements(point_group, elements)
      integer,          intent(in)  :: point_group
      character(len=*), intent(out) :: elements
    end subroutine symmetries_finite_get_group_elements

    subroutine symmetries_finite_end()
    end subroutine symmetries_finite_end
  end interface

  interface symmetries_t
    module procedure symmetries_constructor
  end interface symmetries_t

contains

  function symmetries_constructor(namespace, space, latt, n_sites, site_pos, site_type, spherical_site) result(this)
    type(namespace_t),       intent(in)  :: namespace
    type(space_t),           intent(in)  :: space
    type(lattice_vectors_t), intent(in)  :: latt
    integer,                 intent(in)  :: n_sites
    FLOAT,                   intent(in)  :: site_pos(:,:)
    integer,                 intent(in)  :: site_type(1:n_sites)
    logical,                 intent(in)  :: spherical_site(:)
    type(symmetries_t) :: this

    integer :: max_size, dim4syms
    integer :: idir, isite, isite_symm, iop, verbosity, point_group
    FLOAT   :: lattice(1:3, 1:3)
    FLOAT, allocatable :: position(:, :)
    type(block_t) :: blk
    type(symm_op_t) :: tmpop
    integer :: identity(3,3)
    logical :: found_identity, is_supercell
    integer                  :: fullnops
    integer, allocatable     :: rotation(:, :, :)
    FLOAT,   allocatable     :: translation(:, :)
    character(kind=c_char) :: c_symbol(11), c_schoenflies(7)
    logical :: def_sym_comp
    FLOAT :: rsite(space%dim)

    PUSH_SUB(symmetries_constructor)

    ! if someone cares, they could try to analyze the symmetry point group of the individual species too
    this%any_non_spherical = any(.not. spherical_site)

    this%periodic_dim = space%periodic_dim

    dim4syms = min(3, space%dim)

    def_sym_comp = (n_sites < 100) .or. space%is_periodic()
    def_sym_comp = def_sym_comp .and. space%dim == 3

    !%Variable SymmetriesCompute
    !%Type logical
    !%Section Execution::Symmetries
    !%Description
    !% If disabled, <tt>Octopus</tt> will not compute
    !% nor print the symmetries.
    !%
    !% By default, symmetries are computed when running in 3
    !% dimensions for systems with less than 100 atoms.
    !% For periodic systems, the default is always true, irrespective of the number of atoms.
    !%End
    call parse_variable(namespace, 'SymmetriesCompute', def_sym_comp, this%symmetries_compute)

    !%Variable SymmetriesTolerance
    !%Type float
    !%Section Execution::Symmetries
    !%Description
    !% For periodic systems, this variable controls the tolerance used by the symmetry finder
    !% (spglib) to find the spacegroup and symmetries of the crystal.
    !%End
    call parse_variable(namespace, 'SymmetriesTolerance', CNST(1e-5), symprec)


    if (this%symmetries_compute .and. space%dim /= 3) then
      call messages_experimental('symmetries for non 3D systems', namespace=namespace)
    end if

    if (this%any_non_spherical .or. .not. this%symmetries_compute) then
      call init_identity()

      POP_SUB(symmetries_constructor)
      return
    end if

    ! In all cases, we must check that the grid respects the symmetries. --DAS

    if (space%periodic_dim == 0) then

      call init_identity()

      ! for the moment symmetries are only used for information, so we compute them only on one node.
      if (mpi_grp_is_root(mpi_world)) then
        verbosity = -1

        SAFE_ALLOCATE(position(1:3, 1:n_sites))

        do isite = 1, n_sites
          position(1:3, isite) = M_ZERO
          position(1:dim4syms, isite) = site_pos(1:dim4syms, isite)
        end do

        if (this%symmetries_compute) then
          call symmetries_finite_init(n_sites, site_type(1), position(1, 1), verbosity, point_group)
          if (point_group > -1) then
            call symmetries_finite_get_group_name(point_group, this%group_name)
            call symmetries_finite_get_group_elements(point_group, this%group_elements)
          end if
          call symmetries_finite_end()
        end if
        SAFE_DEALLOCATE_A(position)
      end if

    else

      SAFE_ALLOCATE(position(1:3, 1:n_sites))

      do isite = 1, n_sites
        position(1:3,isite) = M_ZERO

        ! Transform atomic positions to reduced coordinates
        position(1:dim4syms, isite) = latt%cart_to_red(site_pos(1:dim4syms, isite))
        position(1:dim4syms, isite) = position(1:dim4syms, isite) - M_HALF
        do idir = 1, dim4syms
          position(idir, isite) = position(idir, isite) - anint(position(idir, isite))
        end do
        position(1:dim4syms, isite) = position(1:dim4syms,isite) + M_HALF
      end do

      lattice = M_ZERO
      !NOTE: Why "inverse matrix" ? (NTD)
      ! get inverse matrix to extract reduced coordinates for spglib
      lattice(1:space%periodic_dim, 1:space%periodic_dim) = latt%rlattice(1:space%periodic_dim, 1:space%periodic_dim)
      ! transpose the lattice vectors for use in spglib as row-major matrix
      lattice(:,:) = transpose(lattice(:,:))
      ! spglib always assume 3D periodicity, so we need to set the lattice vectors along the non-periodic dimensions.
      ! Here we choose a value that should guarantee that no incorrect symmetries are found along the non-periodic dimensions.
      do idir = space%periodic_dim + 1, 3
        if (idir <= space%dim) then
          lattice(idir, idir) = M_TWO*(maxval(site_pos(idir, :)) - minval(site_pos(idir, :))) + M_ONE
        else
          lattice(idir, idir) = M_ONE
        end if
      end do

      this%space_group = spg_get_international(c_symbol, lattice(1,1), position(1,1), site_type(1), n_sites, symprec)
      this%space_group = spg_get_schoenflies(c_schoenflies, lattice(1, 1), position(1, 1), site_type(1), n_sites, symprec)

      if (this%space_group == 0) then
        message(1) = "Symmetry analysis failed in spglib. Disabling symmetries."
        call messages_warning(1, namespace=namespace)

        do isite = 1, n_sites
          write(message(1),'(a,i6,a,3f12.6,a,3f12.6)') 'type ', site_type(isite), &
            ' reduced coords ', position(:, isite), ' cartesian coords ', site_pos(:, isite)
          call messages_info(1, namespace=namespace)
        end do

        call init_identity()
        POP_SUB(symmetries_init)
        return
      end if

      call string_c_to_f(c_symbol, this%symbol)
      call string_c_to_f(c_schoenflies, this%schoenflies)

      max_size = spg_get_multiplicity(lattice(1, 1), position(1, 1), site_type(1), n_sites, symprec)

      ! spglib returns row-major not column-major matrices!!! --DAS
      SAFE_ALLOCATE(rotation(1:3, 1:3, 1:max_size))
      SAFE_ALLOCATE(translation(1:3, 1:max_size))

      fullnops = spg_get_symmetry(rotation(1, 1, 1), translation(1, 1), &
        max_size, lattice(1, 1), position(1, 1), site_type(1), n_sites, symprec)

      do iop = 1, fullnops
        ! transpose due to array order difference between C and fortran
        rotation(:,:,iop) = transpose(rotation(:,:,iop))
        ! sometimes spglib may return lattice vectors as 'fractional' translations
        translation(:, iop) = translation(:, iop) - anint(translation(:, iop) + M_HALF*symprec)
      end do

      ! we need to check that it is not a supercell, as in the QE routine (sgam_at)
      ! they disable fractional translations if the identity has one, because the sym ops might not form a group.
      ! spglib may return duplicate operations in this case!
      is_supercell = (fullnops > 48)
      found_identity = .false.
      identity = diagonal_matrix(3, 1)
      do iop = 1, fullnops
        if (all(rotation(1:3, 1:3, iop) == identity(1:3, 1:3))) then
          found_identity = .true.
          if (any(abs(translation(1:3, iop)) > TOFLOAT(symprec))) then
            is_supercell = .true.
            write(message(1),'(a,3f12.6)') 'Identity has a fractional translation ', translation(1:3, iop)
            call messages_info(1, namespace=namespace)
          end if
        end if
      end do
      if (.not. found_identity) then
        message(1) = "Symmetries internal error: Identity is missing from symmetry operations."
        call messages_fatal(1, namespace=namespace)
      end if

      if (is_supercell) then
        message(1) = "Disabling fractional translations. System appears to be a supercell."
        call messages_info(1, namespace=namespace)
      end if
      ! actually, we do not use fractional translations regardless currently

      ! this is a hack to get things working, this variable should be
      ! eliminated and the direction calculated automatically from the
      ! perturbations.

      !%Variable SymmetryBreakDir
      !%Type block
      !%Section Mesh::Simulation Box
      !%Description
      !% This variable specifies a direction in which the symmetry of
      !% the system will be broken. This is useful for generating <i>k</i>-point
      !% grids when an external perturbation is applied.
      !%End

      this%breakdir(1:3) = M_ZERO

      if (parse_block(namespace, 'SymmetryBreakDir', blk) == 0) then

        do idir = 1, dim4syms
          call parse_block_float(blk, 0, idir - 1, this%breakdir(idir))
        end do

        call parse_block_end(blk)

      end if

      SAFE_ALLOCATE(this%ops(1:fullnops))


      ! check all operations and leave those that kept the symmetry-breaking
      ! direction invariant and (for the moment) that do not have a translation
      this%nops = 0
      do iop = 1, fullnops
        call symm_op_init(tmpop, rotation(1:3, 1:3, iop), latt, dim4syms, TOFLOAT(translation(1:3, iop)))

        if (symm_op_invariant_cart(tmpop, this%breakdir, TOFLOAT(symprec)) &
          .and. .not. symm_op_has_translation(tmpop, TOFLOAT(symprec))) then
          this%nops = this%nops + 1
          call symm_op_copy(tmpop, this%ops(this%nops))
        end if
      end do

      SAFE_DEALLOCATE_A(position)
      SAFE_DEALLOCATE_A(rotation)
      SAFE_DEALLOCATE_A(translation)

    end if


    ! Checks if that the atomic coordinates are compatible with the symmetries
    !
    ! We want to use for instance that
    !
    ! \int dr f(Rr) V_isite(r) \nabla f(R(v)) = R\int dr f(r) V_isite(R*r) f(r)
    !
    ! and that the operator R should map the position of atom
    ! isite to the position of some other atom isite_symm, so that
    !
    ! V_isite(R*r) = V_isite_symm(r)
    !
    do iop = 1, symmetries_number(this)
      if (iop == symmetries_identity_index(this)) cycle

      do isite = 1, n_sites
        rsite = symm_op_apply_cart(this%ops(iop), site_pos(:, isite))

        rsite = latt%fold_into_cell(rsite)

        ! find isite_symm
        do isite_symm = 1, n_sites
          if (all(abs(rsite - site_pos(:, isite_symm)) < CNST(1.0e-5))) exit
        end do

        if (isite_symm > n_sites) then
          write(message(1),'(a,i6)') 'Internal error: could not find symetric partner for atom number', isite
          write(message(2),'(a,i3,a)') 'with symmetry operation number ', iop, '.'
          call messages_fatal(2, namespace=namespace)
        end if

      end do
    end do


    call symmetries_write_info(this, space, namespace=namespace)

    POP_SUB(symmetries_constructor)

  contains

    subroutine init_identity()

      PUSH_SUB(symmetries_init.init_identity)

      SAFE_ALLOCATE(this%ops(1))
      this%nops = 1
      call symm_op_init(this%ops(1), diagonal_matrix(3, 1), latt, dim4syms)
      this%breakdir = M_ZERO
      this%space_group = 1

      POP_SUB(symmetries_init.init_identity)

    end subroutine init_identity

  end function symmetries_constructor

  ! -------------------------------------------------------------------------------
  subroutine symmetries_copy(lhs, rhs)
    class(symmetries_t),  intent(out) :: lhs
    class(symmetries_t),  intent(in)  :: rhs

    integer :: iop

    PUSH_SUB(symmetries_copy)

    SAFE_ALLOCATE(lhs%ops(1:rhs%nops))
    do iop = 1, rhs%nops
      call symm_op_copy(rhs%ops(iop), lhs%ops(iop))
    end do
    lhs%nops = rhs%nops
    lhs%breakdir = rhs%breakdir
    lhs%periodic_dim = rhs%periodic_dim
    lhs%space_group = rhs%space_group
    lhs%any_non_spherical = rhs%any_non_spherical
    lhs%symmetries_compute = rhs%symmetries_compute
    lhs%group_name = rhs%group_name
    lhs%group_elements = rhs%group_elements
    lhs%symbol = rhs%symbol
    lhs%schoenflies = rhs%schoenflies

    POP_SUB(symmetries_copy)
  end subroutine symmetries_copy

  ! -------------------------------------------------------------------------------

  subroutine symmetries_finalizer(this)
    type(symmetries_t),  intent(inout) :: this

    PUSH_SUB(symmetries_finalizer)

    SAFE_DEALLOCATE_A(this%ops)

    POP_SUB(symmetries_finalizer)
  end subroutine symmetries_finalizer

  ! -------------------------------------------------------------------------------

  integer pure function symmetries_number(this) result(number)
    type(symmetries_t),  intent(in) :: this

    number = this%nops
  end function symmetries_number

  ! -------------------------------------------------------------------------------

  subroutine symmetries_apply_kpoint_red(this, iop, aa, bb)
    type(symmetries_t),  intent(in)  :: this
    integer,             intent(in)  :: iop
    FLOAT,               intent(in)  :: aa(1:3)
    FLOAT,               intent(out) :: bb(1:3)

    PUSH_SUB(symmetries_apply_kpoint_red)

    ASSERT(iop /= 0 .and. abs(iop) <= this%nops)

    !Time-reversal symmetry
    if (iop < 0) then
      bb(1:3) = symm_op_apply_transpose_red(this%ops(abs(iop)), -aa(1:3))
    else
      bb(1:3) = symm_op_apply_transpose_red(this%ops(abs(iop)), aa(1:3))
    end if

    POP_SUB(symmetries_apply_kpoint_red)
  end subroutine symmetries_apply_kpoint_red

  ! -------------------------------------------------------------------------------

  integer pure function symmetries_space_group_number(this) result(number)
    type(symmetries_t),  intent(in) :: this

    number = this%space_group
  end function symmetries_space_group_number

  ! -------------------------------------------------------------------------------

  logical pure function symmetries_have_break_dir(this) result(have)
    type(symmetries_t),  intent(in)  :: this

    have = any(abs(this%breakdir(1:3)) > M_EPSILON)
  end function symmetries_have_break_dir

  ! -------------------------------------------------------------------------------

  integer pure function symmetries_identity_index(this) result(index)
    type(symmetries_t),  intent(in)  :: this

    integer :: iop

    do iop = 1, this%nops
      if (symm_op_is_identity(this%ops(iop))) then
        index = iop
        cycle
      end if
    end do

  end function symmetries_identity_index

  ! ---------------------------------------------------------
  subroutine symmetries_write_info(this, space, iunit, namespace)
    type(symmetries_t),          intent(in) :: this
    type(space_t),               intent(in) :: space
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer :: iop

    PUSH_SUB(symmetries_write_info)

    call messages_print_stress(msg='Symmetries', iunit=iunit, namespace=namespace)

    if (this%any_non_spherical) then
      message(1) = "Symmetries are disabled since non-spherically symmetric species may be present."
      call messages_info(1,iunit = iunit, namespace=namespace)
      call messages_print_stress(iunit=iunit, namespace=namespace)
    end if

    if (.not. this%symmetries_compute) then
      message(1) = "Symmetries have been disabled by SymmetriesCompute = false."
      call messages_info(1, iunit=iunit, namespace=namespace)
      call messages_print_stress(iunit=iunit, namespace=namespace)
      POP_SUB(symmetries_write_info)
      return
    end if

    if (space%periodic_dim == 0) then
      ! At the moment only the root node has information about symetries of finite systems.
      if (mpi_grp_is_root(mpi_world)) then
        if (this%symmetries_compute) then
          call messages_write('Symmetry elements : '//trim(this%group_elements), new_line = .true.)
          call messages_write('Symmetry group    : '//trim(this%group_name))
          call messages_info(iunit=iunit, namespace=namespace)
        end if
      end if
    else
      write(message(1),'(a, i4)') 'Space group No. ', this%space_group
      write(message(2),'(2a)') 'International: ', trim(this%symbol)
      write(message(3),'(2a)') 'Schoenflies: ', trim(this%schoenflies)
      call messages_info(3, iunit=iunit, namespace=namespace)

      write(message(1),'(a7,a31,12x,a33)') 'Index', 'Rotation matrix', 'Fractional translations'
      call messages_info(1, iunit=iunit, namespace=namespace)
      do iop = 1, this%nops
        ! list all operations and leave those that kept the symmetry-breaking
        ! direction invariant and (for the moment) that do not have a translation
        if (space%dim == 1) then
          write(message(1),'(i5,1x,a,2x,1(1i4,2x),1f12.6)') iop, ':', symm_op_rotation_matrix_red(this%ops(iop)), &
            symm_op_translation_vector_red(this%ops(iop))
        end if
        if (space%dim == 2) then
          write(message(1),'(i5,1x,a,2x,2(2i4,2x),2f12.6)') iop, ':', symm_op_rotation_matrix_red(this%ops(iop)), &
            symm_op_translation_vector_red(this%ops(iop))
        end if
        if (space%dim == 3) then
          write(message(1),'(i5,1x,a,2x,3(3i4,2x),3f12.6)') iop, ':', symm_op_rotation_matrix_red(this%ops(iop)), &
            symm_op_translation_vector_red(this%ops(iop))
        end if
        call messages_info(1, iunit=iunit, namespace=namespace)
      end do
      write(message(1), '(a,i5,a)') 'Info: The system has ', this%nops, ' symmetries that can be used.'
      call messages_info(iunit=iunit, namespace=namespace)
    end if
    call messages_print_stress(iunit=iunit, namespace=namespace)

    POP_SUB(symmetries_write_info)
  end subroutine symmetries_write_info


end module symmetries_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
