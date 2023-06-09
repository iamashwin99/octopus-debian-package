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

module read_coords_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                     &
    read_coords_atom,              &
    read_coords_info,              &
    read_coords_init,              &
    read_coords_end,               &
    read_coords_read

  !> for read_coords_info::file_type
  integer, public, parameter :: &
    READ_COORDS_ERR      = 0,      &
    READ_COORDS_PDB      = 1,      &
    READ_COORDS_XYZ      = 2,      &
    READ_COORDS_INP      = 3,      &
    READ_COORDS_REDUCED  = 4,      &
    READ_COORDS_XSF      = 5

  !> for read_coords_info::flags
  integer, public, parameter :: &
    XYZ_FLAGS_RESIDUE = 1,      &
    XYZ_FLAGS_CHARGE  = 2,      &
    XYZ_FLAGS_MOVE    = 4

  type read_coords_atom
    character(len=15) :: label  !< stuff that is always known
    FLOAT             :: x(MAX_DIM)

    FLOAT             :: charge !< stuff specific to PDB files
    character(len=3)  :: residue
    logical           :: move   !< stuff specific to the inp file
  end type read_coords_atom

  type read_coords_info
    integer :: source
    integer :: flags

    integer :: n                !< number of atoms in file
    type(read_coords_atom), allocatable :: atom(:)

    FLOAT, allocatable :: latvec(:,:)
  end type read_coords_info


contains

  ! ---------------------------------------------------------
  subroutine read_coords_init(gf)
    type(read_coords_info), intent(out) :: gf

    PUSH_SUB(read_coords_init)

    gf%source = READ_COORDS_ERR
    gf%flags     = 0
    gf%n         = 0

    POP_SUB(read_coords_init)
  end subroutine read_coords_init


  ! ---------------------------------------------------------
  subroutine read_coords_end(gf)
    type(read_coords_info), intent(inout) :: gf

    PUSH_SUB(read_coords_end)

    SAFE_DEALLOCATE_A(gf%atom)
    call read_coords_init(gf)
    SAFE_DEALLOCATE_A(gf%latvec)

    POP_SUB(read_coords_end)
  end subroutine read_coords_end


  ! ---------------------------------------------------------
  subroutine read_coords_read(what, gf, space, namespace)
    character(len=*),       intent(in)    :: what
    type(read_coords_info), intent(inout) :: gf
    type(space_t),          intent(in)    :: space
    type(namespace_t),      intent(in)    :: namespace

    integer :: ia, ncol, iunit, jdir, int_one, nsteps, istep, step_to_use, periodic_dim
    type(block_t) :: blk
    character(len=256) :: str
    logical :: done

    PUSH_SUB(read_coords_read)

    done = .false.

    !%Variable PDBCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If this variable is present, the program tries to read the atomic coordinates
    !% from the file specified by its value. The PDB (<a href=http://www.rcsb.org/pdb>Protein Data Bank</a>)
    !% format is quite complicated, and it goes
    !% well beyond the scope of this manual. You can find a comprehensive
    !% description <a href=http://www.wwpdb.org/docs.html>here</a>.
    !% From the plethora of instructions defined in the PDB standard, <tt>Octopus</tt>
    !% only reads two, <tt>ATOM</tt> and <tt>HETATOM</tt>. From these fields, it reads:
    !% <ul>
    !% <li> columns 13-16: The species; in fact <tt>Octopus</tt> only cares about the
    !% first letter - <tt>CA</tt> and <tt>CB</tt> will both refer to carbon - so elements whose
    !% chemical symbol has more than one letter cannot be represented in this way.
    !% So, if you want to run mercury (Hg), please use one of the other methods
    !% to input the coordinates.
    !% <li> columns 18-21: The residue. Ignored.
    !% <li> columns 31-54: The Cartesian coordinates. The Fortran format is <tt>(3f8.3)</tt>.</li>
    !% <li> columns 61-65: Classical charge of the atom. Required if reading classical atoms, ignored otherwise.
    !% The Fortran format is <tt>(f6.2)</tt>.</li>
    !% </ul>
    !% NOTE: The coordinates are treated in the units specified by <tt>Units</tt> and/or <tt>UnitsInput</tt>.
    !%End

    !%Variable PDBClassical
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If this variable is present, the program tries to read the atomic coordinates for classical atoms.
    !% from the file specified by its value. The same as <tt>PDBCoordinates</tt>, except that the
    !% classical charge colum must be present. The interaction from the
    !% classical atoms is specified by <tt>ClassicalPotential</tt>, for QM/MM calculations.
    !% Not available in periodic systems.
    !%End

    if (parse_is_defined(namespace, 'PDB'//trim(what))) then
      call check_duplicated(done)

      gf%source = READ_COORDS_PDB
      gf%flags = ior(gf%flags, XYZ_FLAGS_RESIDUE)
      gf%flags = ior(gf%flags, XYZ_FLAGS_CHARGE)

      ! no default, since we do not do this unless the input tag is present
      call parse_variable(namespace, 'PDB'//trim(what), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1, namespace=namespace)

      iunit = io_open(str, namespace, action='read')
      call read_coords_read_PDB(what, iunit, gf)
      call io_close(iunit)
    end if

    ! PDB is the only acceptable format for classical atoms.
    if (trim(what) == "Classical") then
      POP_SUB(read_coords_read)
      return
    end if

    !%Variable XYZCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If <tt>PDBCoordinates</tt> is not present, the program reads the atomic coordinates from
    !% the XYZ file specified by the variable <tt>XYZCoordinates</tt> -- in case this variable
    !% is present. The XYZ format is very simple: The first line of the file has an integer
    !% indicating the number of atoms. The second can contain comments that are simply ignored by
    !% <tt>Octopus</tt>. Then there follows one line per atom, containing the chemical species and
    !% the Cartesian coordinates of the atom.
    !%
    !% If you want to specify the unit of the XYZ file, you can use the variable <tt>UnitsXYZFiles</tt>.
    !%End

    if (parse_is_defined(namespace, 'XYZ'//trim(what))) then ! read an xyz file
      call check_duplicated(done)

      gf%source = READ_COORDS_XYZ
      ! no default, since we do not do this unless the input tag is present
      call parse_variable(namespace, 'XYZ'//trim(what), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1, namespace=namespace)

      iunit = io_open(str, namespace, action='read', status='old')
      read(iunit, *) gf%n

      if (gf%n <= 0) then
        write(message(1),'(a,i6)') "Invalid number of atoms ", gf%n
        call messages_fatal(1, namespace=namespace)
      end if

      read(iunit, *) ! skip comment line

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        read(iunit,*) gf%atom(ia)%label, gf%atom(ia)%x(1:space%dim)
      end do

      call io_close(iunit)
    end if

    !%Variable XSFCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% Another option besides PDB and XYZ coordinates formats is XSF, as <a href=http://www.xcrysden.org/doc/XSF.html>defined</a>
    !% by the XCrySDen visualization program. Specify the filename with this variable.
    !% The lattice vectors will also be read from this file and the value of
    !% <tt>PeriodicDimensions</tt> needs to be compatible with the first line
    !% (<tt>CRYSTAL</tt>, <tt>SLAB</tt>, <tt>POLYMER</tt>, or <tt>MOLECULE</tt>).
    !% The file should not contain <tt>ATOMS</tt>, <tt>CONVVEC</tt>, or <tt>PRIMCOORD</tt>.
    !% NOTE: The coordinates are treated in the units specified by <tt>Units</tt> and/or <tt>UnitsInput</tt>.
    !%End

    if (parse_is_defined(namespace, 'XSF'//trim(what))) then ! read an xsf file
      call check_duplicated(done)

      gf%source = READ_COORDS_XSF
      ! no default, since we do not do this unless the input tag is present
      call parse_variable(namespace, 'XSF'//trim(what), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1, namespace=namespace)

      iunit = io_open(str, namespace, action='read', status='old')

      read(iunit, '(a256)') str
      if (str(1:9) == 'ANIMSTEPS') then
        read(str(10:), *) nsteps
        read(iunit, *) str

        !%Variable XSFCoordinatesAnimStep
        !%Type integer
        !%Default 1
        !%Section System::Coordinates
        !%Description
        !% If an animated file is given with <tt>XSFCoordinates</tt>, this variable selects which animation step
        !% will be used. The <tt>PRIMVEC</tt> block must be written for each step.
        !%End
        call parse_variable(namespace, 'XSFCoordinatesAnimStep', 1, step_to_use)
        if (step_to_use < 1) then
          message(1) = "XSFCoordinatesAnimStep must be > 0."
          call messages_fatal(1, namespace=namespace)
        else if (step_to_use > nsteps) then
          write(message(1),'(a,i9)') "XSFCoordinatesAnimStep must be <= available number of steps ", nsteps
          call messages_fatal(1, namespace=namespace)
        end if
        write(message(1),'(a,i9,a,i9)') 'Using animation step ', step_to_use, ' out of ', nsteps
        call messages_info(1, namespace=namespace)
      else
        nsteps = 0
        step_to_use = 1
      end if

      ! periodicity = 'CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE'
      select case (trim(str))
      case ('CRYSTAL')
        periodic_dim = 3
      case ('SLAB')
        periodic_dim = 2
      case ('POLYMER')
        periodic_dim = 1
      case ('MOLECULE')
        periodic_dim = 0
      case ('ATOMS')
        call messages_not_implemented("Input from XSF file beginning with ATOMS", namespace=namespace)
      case default
        write(message(1),'(3a)') 'Line in file was "', trim(str), '" instead of CRYSTAL/SLAB/POLYMER/MOLECULE.'
        call messages_fatal(1, namespace=namespace)
      end select
      if (periodic_dim /= space%periodic_dim) then
        message(1) = "Periodicity in XSF input is incompatible with the value of PeriodicDimensions."
        call messages_fatal(1, namespace=namespace)
      end if

      do istep = 1, step_to_use - 1
        read(iunit, *) ! PRIMVEC
        do jdir = 1, space%dim
          read(iunit, *) ! lattice vectors
        end do
        read(iunit, *) ! PRIMCOORD istep
        read(iunit, *) gf%n, int_one
        do ia = 1, gf%n
          read(iunit, *) ! atoms
        end do
      end do

      read(iunit, '(a256)') str
      if (str(1:7) /= 'PRIMVEC') then
        write(message(1),'(3a)') 'Line in file was "', trim(str), '" instead of "PRIMVEC".'
        call messages_warning(1, namespace=namespace)
      end if
      if (nsteps > 0) then
        read(str(8:), *) istep
        if (istep /= step_to_use) then
          write(message(1), '(a,i9,a,i9)') 'Found PRIMVEC index ', istep, ' instead of ', step_to_use
          call messages_fatal(1, namespace=namespace)
        end if
      end if

      SAFE_ALLOCATE(gf%latvec(1:space%dim, 1:space%dim))
      gf%latvec = M_ZERO
      do jdir = 1, space%dim
        read(iunit, *) gf%latvec(1:space%dim, jdir)
      end do
      gf%latvec = units_to_atomic(units_inp%length, gf%latvec)

      read(iunit, '(a256)') str
      if (str(1:9) /= 'PRIMCOORD') then
        write(message(1),'(3a)') 'Line in file was "', trim(str), '" instead of "PRIMCOORD".'
        call messages_warning(1, namespace=namespace)
      end if
      if (nsteps > 0) then
        read(str(10:), *) istep
        if (istep /= step_to_use) then
          write(message(1), '(a,i9,a,i9)') 'Found PRIMCOORD index ', istep, ' instead of ', step_to_use
          call messages_fatal(1, namespace=namespace)
        end if
      end if

      read(iunit, *) gf%n, int_one
      if (gf%n <= 0) then
        write(message(1),'(a,i6)') "Invalid number of atoms ", gf%n
        call messages_fatal(1, namespace=namespace)
      end if
      if (int_one /= 1) then
        write(message(1),'(a,i6,a)') 'Number in file was ', int_one, ' instead of 1.'
        call messages_warning(1, namespace=namespace)
      end if
      SAFE_ALLOCATE(gf%atom(1:gf%n))

      ! TODO: add support for velocities as vectors here?
      do ia = 1, gf%n
        read(iunit,*) gf%atom(ia)%label, gf%atom(ia)%x(1:space%dim)
      end do

      call io_close(iunit)
    end if

    !%Variable Coordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% If <tt>XYZCoordinates</tt>, <tt>PDBCoordinates</tt>, and <tt>XSFCoordinates</tt> were not found,
    !% <tt>Octopus</tt> tries to read the coordinates for the atoms from the block <tt>Coordinates</tt>. The
    !% format is quite straightforward:
    !%
    !% <tt>%Coordinates
    !% <br>&nbsp;&nbsp;'C' |      -0.56415 | 0.0 | 0.0 | no
    !% <br>&nbsp;&nbsp;'O' | &nbsp;0.56415 | 0.0 | 0.0 | no
    !% <br>%</tt>
    !%
    !% The first line defines a carbon atom at coordinates (-0.56415, 0.0, 0.0),
    !% that is <b>not</b> allowed to move during dynamical simulations. The second line has
    !% a similar meaning. This block obviously defines a carbon monoxide molecule, if the
    !% input units are <tt>eV_Angstrom</tt>. The number of coordinates for each species
    !% must be equal to the dimension of your space (generally 3).
    !% Note that in this way it is possible to fix some of the atoms (this
    !% is not possible when specifying the coordinates through a <tt>PDBCoordinates</tt> or
    !% <tt>XYZCoordinates</tt> file). The last column is optional, and the default is yes.
    !% It is always possible to fix <b>all</b> atoms using the <tt>MoveIons</tt> directive.
    !%End

    if (parse_block(namespace, trim(what), blk) == 0) then
      call check_duplicated(done)

      gf%n = parse_block_n(blk)

      gf%source = READ_COORDS_INP
      gf%flags = ior(gf%flags, XYZ_FLAGS_MOVE)

      message(1) = "Reading " // trim(what) // " from " // trim(what) // " block"
      call messages_info(1, namespace=namespace)

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        ncol = parse_block_cols(blk, ia - 1)
        if ((ncol  <  space%dim + 1) .or. (ncol > space%dim + 2)) then
          write(message(1), '(3a,i2)') 'Error in block ', what, ' line #', ia
          call messages_fatal(1, namespace=namespace)
        end if
        call parse_block_string (blk, ia - 1, 0, gf%atom(ia)%label)
        do jdir = 1, space%dim
          call parse_block_float  (blk, ia - 1, jdir, gf%atom(ia)%x(jdir))
        end do
        if (ncol == space%dim + 2) then
          call parse_block_logical(blk, ia - 1, space%dim + 1, gf%atom(ia)%move)
        else
          gf%atom(ia)%move = .true.
        end if
      end do

      call parse_block_end(blk)
    end if

    !%Variable ReducedCoordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% This block gives the atomic coordinates relative to the real
    !% space unit cell. The format is the same as the
    !% <tt>Coordinates</tt> block.
    !%
    !% Note that in Octopus the origin of coordinates is in the center
    !% of the cell, so the coordinates inside the cell are in the
    !% range [-0.5, 0.5).
    !%
    !% This block cannot be used with the <tt>minimum</tt> box shape.
    !%End

    ! This is valid only for Coordinates.
    if (trim(what) == 'Coordinates') then
      if (parse_block(namespace, 'Reduced'//trim(what), blk) == 0) then
        call check_duplicated(done)

        gf%n = parse_block_n(blk)

        gf%source = READ_COORDS_REDUCED
        gf%flags = ior(gf%flags, XYZ_FLAGS_MOVE)

        message(1) = "Reading " // trim(what) // " from Reduced" // trim(what) // " block"
        call messages_info(1, namespace=namespace)

        SAFE_ALLOCATE(gf%atom(1:gf%n))

        do ia = 1, gf%n
          ncol = parse_block_cols(blk, ia - 1)
          if ((ncol  <  space%dim + 1) .or. (ncol > space%dim + 2)) then
            write(message(1), '(3a,i2)') 'Error in block ', what, ' line #', ia
            call messages_fatal(1, namespace=namespace)
          end if
          call parse_block_string (blk, ia - 1, 0, gf%atom(ia)%label)
          do jdir = 1, space%dim
            call parse_block_float  (blk, ia - 1, jdir, gf%atom(ia)%x(jdir))
          end do
          do jdir = space%dim + 1, MAX_DIM
            gf%atom(ia)%x(jdir) = M_ZERO
          end do
          if (ncol == space%dim + 2) then
            call parse_block_logical(blk, ia - 1, space%dim + 1, gf%atom(ia)%move)
          else
            gf%atom(ia)%move = .true.
          end if
        end do

        call parse_block_end(blk)
      end if
    end if

    if (gf%source /= READ_COORDS_REDUCED) then
      ! adjust units
      do ia = 1, gf%n
        do jdir = space%dim + 1, MAX_DIM
          gf%atom(ia)%x(jdir) = M_ZERO
        end do

        if (gf%source == READ_COORDS_XYZ) then
          gf%atom(ia)%x = units_to_atomic(units_inp%length_xyz_file, gf%atom(ia)%x)
        else
          gf%atom(ia)%x = units_to_atomic(units_inp%length, gf%atom(ia)%x)
        end if

      end do
    end if

    POP_SUB(read_coords_read)

  contains

    subroutine check_duplicated(done)
      logical, intent(inout) :: done

      PUSH_SUB(read_coords_read.check_duplicated)

      if (.not. done) then
        done = .true.
      else
        message(1) = 'Multiple definitions of '//trim(what)//' in the input file.'
        call messages_fatal(1, namespace=namespace)
      end if

      POP_SUB(read_coords_read.check_duplicated)
    end subroutine check_duplicated

  end subroutine read_coords_read


  ! ---------------------------------------------------------
  subroutine read_coords_read_PDB(what, iunit, gf)
    character(len=*),    intent(in)    :: what
    integer,             intent(in)    :: iunit
    type(read_coords_info), intent(inout) :: gf

    character(len=80) :: record
    character(len=6)  :: record_name
    integer :: na

    PUSH_SUB(read_coords_read_PDB)

    ! TODO: read the specification of the crystal structure.

    ! First count number of atoms
    rewind(iunit)
    do
      read(iunit, '(a80)', err=990, end=990) record
      read(record, '(a6)') record_name
      if (trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        gf%n = gf%n + 1
      end if
    end do
990 continue

    SAFE_ALLOCATE(gf%atom(1:gf%n))

    ! read in the data
    rewind(iunit)
    na = 1
    do
      read(iunit, '(a80)', err=991, end=991) record
      read(record, '(a6)') record_name
      if (trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(12x,a4,1x,a3)') gf%atom(na)%label, gf%atom(na)%residue
        call str_trim(gf%atom(na)%label)
        gf%atom(na)%label = gf%atom(na)%label(1:1)
        call str_trim(gf%atom(na)%residue)

        if (trim(what) == 'Classical') then
          read(record, '(30x,3f8.3,6x,f5.2)') gf%atom(na)%x(1:3), gf%atom(na)%charge
        else
          read(record, '(30x,3f8.3)') gf%atom(na)%x(1:3)
        end if

        na = na + 1
      end if
    end do
991 continue

    POP_SUB(read_coords_read_PDB)
  end subroutine read_coords_read_PDB

end module read_coords_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
