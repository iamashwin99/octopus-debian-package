!! Copyright (C) 2009 D. Strubbe
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

module born_charges_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use ions_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use species_oct_m
  use unit_system_oct_m
  use utils_oct_m

  implicit none

  private
  public ::                &
    born_charges_t,        &
    born_charges_init,     &
    born_charges_end,      &
    born_output_charges

  type born_charges_t
    private
    CMPLX, allocatable, public :: charge(:, :, :)    !< i, j, atom: Z*(i,j) = dF(j)/dE(i) = dP(i) / dR(j)
    CMPLX :: sum_ideal(MAX_DIM, MAX_DIM) !< the sum of born charges according to acoustic sum rule
    CMPLX :: delta(MAX_DIM, MAX_DIM)     !< discrepancy of sum of born charge tensors from sum rule, per atom
    logical :: correct                   !< correct according to sum rule?
  end type born_charges_t

contains

  ! ---------------------------------------------------------
  subroutine born_charges_init(this, namespace, ions, val_charge, qtot, dim)
    type(born_charges_t), intent(out) :: this
    type(namespace_t),    intent(in)  :: namespace
    type(ions_t),         intent(in)  :: ions
    FLOAT,                intent(in)  :: val_charge
    FLOAT,                intent(in)  :: qtot
    integer,              intent(in)  :: dim

    integer :: idir

    PUSH_SUB(born_charges_init)

    SAFE_ALLOCATE(this%charge(1:dim, 1:dim, 1:ions%natoms))
    this%charge(1:dim, 1:dim, 1:ions%natoms) = M_ZERO
    this%delta(1:dim, 1:dim) = M_ZERO

    this%sum_ideal(1:dim, 1:dim) = M_ZERO
    do idir = 1, dim
      this%sum_ideal(idir, idir) = -(val_charge + qtot) ! total charge
    end do

    !%Variable BornChargeSumRuleCorrection
    !%Type logical
    !%Default true
    !%Section Linear Response::Polarizabilities
    !%Description
    !% Enforce the acoustic sum rule by distributing the excess sum of Born charges equally among the atoms.
    !% Sum rule: <math>\sum_{\alpha} Z^{*}_{\alpha, i, j} = Z_{\rm tot} \delta_{ij}</math>.
    !% Violation of the sum rule may be caused by inadequate spacing, box size (in finite directions),
    !% or <i>k</i>-point sampling (in periodic directions).
    !%End

    call parse_variable(namespace, 'BornChargeSumRuleCorrection', .true., this%correct)

    POP_SUB(born_charges_init)
  end subroutine born_charges_init

  ! ---------------------------------------------------------
  subroutine born_charges_end(this)
    type(born_charges_t), intent(inout) :: this

    PUSH_SUB(born_charges_end)

    SAFE_DEALLOCATE_A(this%charge)

    POP_SUB(born_charges_end)
  end subroutine born_charges_end

  ! ---------------------------------------------------------
  !> The sum over atoms of a given tensor component of the born charges
  !!  should be Z delta_ij to satisfy the acoustic sum rule, where Z is total charge of system
  subroutine correct_born_charges(this, ions, dim)
    type(born_charges_t), intent(inout) :: this
    type(ions_t),         intent(in)    :: ions
    integer,              intent(in)    :: dim

    CMPLX :: born_sum(dim, dim)        ! the sum of born charges from the calculation
    integer :: iatom

    PUSH_SUB(correct_born_charges)

    born_sum(1:dim, 1:dim) = M_ZERO

    do iatom = 1, ions%natoms
      born_sum(1:dim, 1:dim) = born_sum(1:dim, 1:dim) + this%charge(1:dim, 1:dim, iatom)
    end do

    this%delta(1:dim, 1:dim) = (born_sum(1:dim, 1:dim) - this%sum_ideal(1:dim, 1:dim)) / ions%natoms

    if (this%correct) then
      do iatom = 1, ions%natoms
        this%charge(1:dim, 1:dim, iatom) = &
          this%charge(1:dim, 1:dim, iatom) - this%delta(1:dim, 1:dim)
      end do
    end if

    POP_SUB(correct_born_charges)
  end subroutine correct_born_charges

  ! ---------------------------------------------------------
  subroutine born_output_charges(this, ions, namespace, dim, dirname, write_real)
    type(born_charges_t), intent(inout) :: this
    type(ions_t),         intent(in)    :: ions
    type(namespace_t),    intent(in)    :: namespace
    integer,              intent(in)    :: dim
    character(len=*),     intent(in)    :: dirname
    logical,              intent(in)    :: write_real
    !< set write_real to true if they are all real, to suppress writing imaginary part and phase

    integer iatom, iunit
    FLOAT :: phase(dim, dim)

    PUSH_SUB(born_output_charges)

    call correct_born_charges(this, ions, dim)

    if (mpi_grp_is_root(mpi_world)) then ! only first node outputs
      iunit = io_open(trim(dirname)//'/born_charges', namespace, action='write')
      write(iunit,'(a)') '# (Frequency-dependent) Born effective charge tensors'
      if (.not. write_real) write(iunit,'(a)') '# Real and imaginary parts'
      do iatom = 1, ions%natoms
        write(iunit,'(a,i5,a,a5,a,f10.4)') 'Index: ', iatom, '   Label: ', trim(species_label(ions%atom(iatom)%species)), &
          '   Ionic charge: ', species_zval(ions%atom(iatom)%species)

        if (.not. write_real) write(iunit,'(a)') 'Real:'
        call output_tensor(TOFLOAT(this%charge(:, :, iatom)), dim, unit_one, iunit=iunit)

        if (.not. write_real) then
          write(iunit,'(a)') 'Imaginary:'
          call output_tensor(aimag(this%charge(:, :, iatom)), dim, unit_one, iunit=iunit)
        end if

        write(iunit,'(a)')
      end do

      if (.not. write_real) then
        write(iunit,'(a)') '# Magnitude and phase'
        do iatom = 1, ions%natoms
          write(iunit,'(a,i5,a,a5,a,f10.4)') 'Index: ', iatom, '   Label: ', trim(species_label(ions%atom(iatom)%species)), &
            '   Ionic charge: ', species_zval(ions%atom(iatom)%species)

          write(iunit,'(a)') 'Magnitude:'
          call output_tensor(TOFLOAT(abs(this%charge(:, :, iatom))), dim, unit_one, iunit=iunit)

          write(iunit,'(a)') 'Phase:'

          where (abs(this%charge(1:dim, 1:dim, iatom)) > M_EPSILON)
            phase(1:dim, 1:dim) = atan2(aimag(this%charge(1:dim, 1:dim, iatom)), TOFLOAT(this%charge(1:dim, 1:dim, iatom)))
          else where
            phase(1:dim, 1:dim) = M_ZERO
          end where
          call output_tensor(phase(:, :), dim, unit_one, write_average = .false., iunit=iunit)
          write(iunit,'(a)')
        end do
      end if

      write(iunit,'(a)') '# Discrepancy of Born effective charges from acoustic sum rule before correction, per atom'
      if (.not. write_real) write(iunit,'(a)') 'Real:'
      call output_tensor(TOFLOAT(this%delta(:, :)), dim, unit_one, iunit=iunit)
      if (.not. write_real) then
        write(iunit,'(a)') 'Imaginary:'
        call output_tensor(aimag(this%delta(:, :)), dim, unit_one, iunit=iunit)
      end if

      call io_close(iunit)
    end if

    POP_SUB(born_output_charges)
  end subroutine born_output_charges

end module born_charges_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
