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
!
!> Atomic weights should be read in "atomic mass units" (u) (not to
!! be confused with mass in "atomic units"), that is, it should be given
!! the relative atomic weight). 1 u is roughly the mass of the proton,
!! and exactly one twelfth of mass of the ^{12}C isotope. The relation of the
!! atomic mass unit and the atomic unit of mass, au_[mass], is:
!!
!! 1 au_[mass] = 5.485799110e-4 u
!!
module unit_system_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use unit_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                  &
    unit_system_t,           &
    unit_system_init,        &
    unit_system_get,         &
    unit_system_from_file

  type unit_system_t
    ! Components are public by default
    type(unit_t) :: length
    type(unit_t) :: length_xyz_file
    type(unit_t) :: energy
    type(unit_t) :: time
    type(unit_t) :: velocity
    type(unit_t) :: mass
    type(unit_t) :: force
    type(unit_t) :: acceleration
    type(unit_t) :: polarizability
    type(unit_t) :: hyperpolarizability
  end type unit_system_t

  !> the units systems for reading and writing
  type(unit_system_t), public :: units_inp, units_out

  !> some special units required for particular quantities
  type(unit_t),        public :: unit_one           !< For unitless quantities and arithmetics with units.
  type(unit_t),        public :: unit_angstrom      !< For XYZ files.
  type(unit_t),        public :: unit_ppm           !< Parts per million.
  type(unit_t),        public :: unit_debye         !< For dipoles.
  type(unit_t),        public :: unit_invcm         !< For vibrational frequencies.
  type(unit_t),        public :: unit_susc_ppm_cgs  !< Some magnetic stuff.
  type(unit_t),        public :: unit_kelvin        !< For converting energies into temperatures.
  type(unit_t),        public :: unit_femtosecond   !< Time in femtoseconds.
  type(unit_t),        public :: unit_amu           !< Mass in atomic mass units (AKA Dalton).
  type(unit_t),        public :: unit_kilobytes     !< For small amounts of data (natural code units are bytes)
  type(unit_t),        public :: unit_megabytes     !< For large amounts of data (natural code units are bytes)
  type(unit_t),        public :: unit_gigabytes     !< For larger amounts of data (natural code units are bytes)
  type(unit_t),        public :: unit_eV            !< For output energies in eV.

  integer, parameter, public :: UNITS_ATOMIC = 0, UNITS_EVA = 1, UNITS_FS = 2

contains


  ! ---------------------------------------------------------
  subroutine unit_system_init(namespace)
    type(namespace_t), intent(in) :: namespace

    integer :: cc, cinp, cout, xyz_units

    PUSH_SUB(unit_system_init)

    !%Variable Units
    !%Type virtual
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% (Virtual) These are the units that can be used in the input file.
    !%
    !%Option angstrom        1.8897261328856432
    !%Option pm              0.018897261328856432
    !%Option picometer       0.018897261328856432
    !%Option nm              18.897261328856432
    !%Option nanometer       18.897261328856432
    !%Option ry              0.5
    !%Option rydberg         0.5
    !%Option ev              0.03674932539796232
    !%Option electronvolt    0.03674932539796232
    !%Option invcm           4.5563353e-06
    !%Option kelvin          3.1668105e-06
    !%Option kjoule_mol      0.00038087988
    !%Option kcal_mol        0.0015936014
    !%Option as              0.0413413737896
    !%Option attosecond      0.0413413737896
    !%Option fs              41.3413737896
    !%Option femtosecond     41.3413737896
    !%Option ps              41341.3737896
    !%Option picosecond      41341.3737896
    !%Option c               137.035999139
    !%End

    !%Variable UnitsOutput
    !%Type integer
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% This variable selects the units that Octopus use for output.
    !%
    !% Atomic units seem to be the preferred system in the atomic and
    !% molecular physics community. Internally, the code works in
    !% atomic units. However, for output, some people like
    !% to use a system based on electron-Volts (eV) for energies
    !% and Angstroms (&Aring;) for length.
    !%
    !% Normally time units are derived from energy and length units,
    !% so it is measured in <math>\hbar</math>/Hartree or
    !% <math>\hbar</math>/eV.
    !%
    !% Warning 1: All files read on input will also be treated using
    !% these units, including XYZ geometry files.
    !%
    !% Warning 2: Some values are treated in their most common units,
    !% for example atomic masses (a.m.u.), electron effective masses
    !% (electron mass), vibrational frequencies
    !% (cm<sup>-1</sup>) or temperatures (Kelvin). The unit of charge is always
    !% the electronic charge <i>e</i>.
    !%
    !%Option atomic        0
    !% Atomic units.
    !%Option ev_angstrom   1
    !% Electronvolts for energy, Angstroms for length, the rest of the
    !% units are derived from these and <math>\hbar=1</math>.
    !%End

    if (parse_is_defined(namespace, 'Units') .or. parse_is_defined(namespace, 'Units')) then
      call messages_write("The 'Units' variable is obsolete. Now Octopus always works in atomic", new_line = .true.)
      call messages_write("units. For different units you can use values like 'angstrom', 'eV' ", new_line = .true.)
      call messages_write("and others in the input file.")
      call messages_fatal()
    end if

    call messages_obsolete_variable(namespace, 'Units')
    call messages_obsolete_variable(namespace, 'UnitsInput')

    cinp = UNITS_ATOMIC

    call parse_variable(namespace, 'UnitsOutput', UNITS_ATOMIC, cc)
    if (.not. varinfo_valid_option('Units', cc, is_flag = .true.)) call messages_input_error(namespace, 'UnitsOutput')
    cout = cc

    unit_one%factor = M_ONE
    unit_one%abbrev = '1'
    unit_one%name   = 'one'

    unit_angstrom%factor = P_Ang  ! 1 a.u. = 0.529 A
    unit_angstrom%abbrev = "A"
    unit_angstrom%name   = "Angstrom"

    unit_ppm%factor = CNST(1e-6)
    unit_ppm%abbrev = 'ppm a.u.'
    unit_ppm%name   = 'parts per million'

    unit_susc_ppm_cgs%factor = CNST(1e-6)/CNST(8.9238878e-2)
    unit_susc_ppm_cgs%abbrev = 'ppm cgs/mol'
    unit_susc_ppm_cgs%name   = 'magnetic susceptibility parts per million cgs'

    unit_debye%factor = M_ONE/CNST(2.5417462)
    unit_debye%abbrev = 'Debye'
    unit_debye%name   = 'Debye'

    ! this is the factor to convert 1 Ha to hc/cm in energy units
    unit_invcm%factor = M_ONE/CNST(219474.63)
    unit_invcm%abbrev = 'cm^-1'
    unit_invcm%name   = 'h times c over centimeters'

    unit_kelvin%factor = P_KB
    unit_kelvin%abbrev = 'K'
    unit_kelvin%name   = 'degrees Kelvin'

    !atomic mass units
    unit_amu%factor = M_ONE/CNST(5.485799110e-4)
    unit_amu%abbrev = 'u'
    unit_amu%name   = '1/12 of the mass of C^12'

    unit_femtosecond%factor = CNST(1.0)/CNST(0.024188843)
    unit_femtosecond%abbrev = 'fs'
    unit_femtosecond%name   = 'femtoseconds'

    unit_kilobytes%factor = CNST(2.0)**10
    unit_kilobytes%abbrev = 'KiB'
    unit_kilobytes%name   = 'kibibytes'

    unit_megabytes%factor = CNST(2.0)**20
    unit_megabytes%abbrev = 'MiB'
    unit_megabytes%name   = 'mebibytes'

    unit_gigabytes%factor = CNST(2.0)**30
    unit_gigabytes%abbrev = 'GiB'
    unit_gigabytes%name   = 'gibibytes'

    unit_eV%abbrev = "eV"
    unit_eV%name   = "electronvolt"
    unit_eV%factor = M_ONE/(M_TWO*P_Ry)   ! 1 a.u. = 27.2 eV

    call unit_system_get(units_inp, cinp)
    call unit_system_get(units_out, cout)

    !%Variable UnitsXYZFiles
    !%Type integer
    !%Default angstrom_units
    !%Section Execution::Units
    !%Description
    !% This variable selects in which units I/O of XYZ files should be
    !% performed.
    !%Option bohr_units       0
    !% The XYZ will be assumed to be in Bohr atomic units.
    !%Option angstrom_units   1
    !% XYZ files will be assumed to be always in Angstrom,
    !% independently of the units used by Octopus. This ensures
    !% compatibility with most programs, that assume XYZ files have
    !% coordinates in Angstrom.
    !%End

    call parse_variable(namespace, 'UnitsXYZFiles', OPTION__UNITSXYZFILES__ANGSTROM_UNITS, xyz_units)

    if (.not. varinfo_valid_option('UnitsXYZFiles', xyz_units)) then
      call messages_input_error(namespace, 'UnitsXYZFiles', 'Invalid option')
    end if

    select case (xyz_units)

    case (OPTION__UNITSXYZFILES__BOHR_UNITS)
      ! Use units_inp%length for initialization, as units_inp are always in atomic units.
      units_inp%length_xyz_file = units_inp%length
      units_out%length_xyz_file = units_inp%length

    case (OPTION__UNITSXYZFILES__ANGSTROM_UNITS)
      units_inp%length_xyz_file = unit_angstrom
      units_out%length_xyz_file = unit_angstrom

    case default
      ASSERT(.false.)
    end select

    POP_SUB(unit_system_init)

  end subroutine unit_system_init

  ! ---------------------------------------------------------
  subroutine unit_system_get(uu, cc)
    type(unit_system_t), intent(out) :: uu
    integer,             intent(in)  :: cc

    PUSH_SUB(unit_system_get)

    select case (cc)
    case (UNITS_ATOMIC)
      call unit_system_init_atomic(uu)
    case (UNITS_EVA)
      call unit_system_init_eV_Ang(uu)
    case default
      message(1) = "Unknown units in unit_system_get"
      call messages_fatal(1)
    end select

    POP_SUB(unit_system_get)
  end subroutine unit_system_get


  ! ---------------------------------------------------------
  !> These routines output the unit-conversion factors, defined by
  !! [a.u.] = input*u.unit, output = [a.u.]/u.unit
  ! ---------------------------------------------------------
  subroutine unit_system_init_atomic(uu)
    type(unit_system_t), intent(out) :: uu

    PUSH_SUB(unit_system_init_atomic)

    uu%length%abbrev = "b"
    uu%length%name   = "Bohr"
    uu%length%factor = M_ONE

    uu%energy%abbrev = "H"
    uu%energy%name   = "Hartree"
    uu%energy%factor = M_ONE

    uu%time%abbrev = "hbar/H"
    uu%time%name   = "hbar/Hartree"
    uu%time%factor = M_ONE/uu%energy%factor

    uu%velocity%abbrev = "bH(2pi/h)"
    uu%velocity%name   = "Bohr times Hartree over hbar"
    uu%velocity%factor = M_ONE

    uu%mass%abbrev   = "me"
    uu%mass%name     = "electron mass"
    uu%mass%factor   = M_ONE

    uu%force%abbrev  = "H/b"
    uu%force%name    = "Hartree/Bohr"
    uu%force%factor  = M_ONE

    uu%acceleration%abbrev = "bH(2pi/h)^2"
    uu%acceleration%name   = "Bohr times (Hartree over h bar) squared"
    uu%acceleration%factor = M_ONE

    uu%polarizability%abbrev  = "b^3"
    uu%polarizability%name    = "Bohr^3"
    uu%polarizability%factor  = M_ONE
    ! By convention, this unit appears more commonly than the
    ! equivalent b^2/H. It does not depend on the dimensionality
    ! of the system, despite analogies to a volume.

    uu%hyperpolarizability%abbrev  = "b^5"
    uu%hyperpolarizability%name    = "Bohr^5"
    uu%hyperpolarizability%factor  = M_ONE

    POP_SUB(unit_system_init_atomic)
  end subroutine unit_system_init_atomic


  ! ---------------------------------------------------------
  subroutine unit_system_init_eV_Ang(uu)
    type(unit_system_t), intent(out) :: uu

    PUSH_SUB(unit_system_init_eV_Ang)

    uu%length%abbrev = "A"
    uu%length%name   = "Angstrom"
    uu%length%factor = P_Ang  ! 1 a.u. = 0.529 A

    uu%energy%abbrev = "eV"
    uu%energy%name   = "electronvolt"
    uu%energy%factor = M_ONE/(M_TWO*P_Ry)   ! 1 a.u. = 27.2 eV

    uu%time%abbrev = "hbar/eV"
    uu%time%name   = "hbar/electronvolt"
    uu%time%factor = M_ONE/uu%energy%factor

    uu%velocity%abbrev = "A eV/hbar"
    uu%velocity%name   = "Angstrom times electronvolts over hbar"
    uu%velocity%factor = uu%length%factor*uu%energy%factor

    uu%mass%abbrev   = "hbar^2/(eV A^2)"
    uu%mass%name     = "hbar^2/(electronvolt * Angstrom^2)"
    uu%mass%factor   = M_ONE/(uu%energy%factor *  uu%length%factor**2)

    uu%force%abbrev  = "eV/A"
    uu%force%name    = "electronvolt/Angstrom"
    uu%force%factor  = uu%energy%factor/uu%length%factor

    uu%acceleration%abbrev = "A (eV/hbar)^2"
    uu%acceleration%name   = "Angstrom times (electronvolt over hbar) squared"
    uu%acceleration%factor = uu%length%factor/uu%time%factor**2

    uu%polarizability      = uu%length**3
    uu%hyperpolarizability = uu%length**5

    POP_SUB(unit_system_init_eV_Ang)
  end subroutine unit_system_init_eV_Ang


  ! ---------------------------------------------------------
  !> This is a very primitive procedure that attempts to find out
  !! which units were used to write an octopus file, whether
  !! "multipoles", "cross_section_tensor", etc.
  !! \todo  Although it seems to work in most cases, it is obviously
  !! a very weak code.
  ! ---------------------------------------------------------
  subroutine unit_system_from_file(uu, fname, namespace, ierr)
    type(unit_system_t), intent(inout) :: uu
    character(len=*),    intent(in)    :: fname
    type(namespace_t),   intent(in)    :: namespace
    integer,             intent(inout) :: ierr

    integer            :: iunit, ios
    character(len=256) :: line

    PUSH_SUB(unit_system_from_file)

    iunit = io_open(trim(fname), namespace, action='read', status='old', die=.false.)
    if (iunit < 0) then
      ierr = -2
      POP_SUB(unit_system_from_file)
      return
    end if

    ierr = 0
    do
      read(iunit, '(a)', iostat = ios) line
      if (ios /= 0) exit
      if (index(line,'[A]') /= 0 .or. index(line,'eV') /= 0) then
        call unit_system_get(uu, UNITS_EVA)
        POP_SUB(unit_system_from_file)
        return
      else if (index(line,'[b]') /= 0) then
        call unit_system_get(uu, UNITS_ATOMIC)
        POP_SUB(unit_system_from_file)
        return
      end if
    end do

    ierr = -1

    POP_SUB(unit_system_from_file)
  end subroutine unit_system_from_file

end module unit_system_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
