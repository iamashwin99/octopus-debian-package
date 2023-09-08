!! Copyright (C) 2007 X. Andrade
!! Copyright (C) 2021 N. Tancogne-Dejean
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

module perturbation_ionic_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use perturbation_oct_m
  use profiling_oct_m
  use projector_oct_m
  use space_oct_m
  use species_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use vibrations_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                           &
    perturbation_ionic_t,                           &
    dionic_pert_matrix_elements_2,    &
    zionic_pert_matrix_elements_2

  type, extends(perturbation_t) :: perturbation_ionic_t
    private

    type(ions_t), pointer :: ions => null()

    !> if pure_dir is .false. then the perturbation is a combination of
    !! displacements of atoms.
    !! If pure_dir is .true., next mix1 and mix2 arrays are allocated
    !! If pure_dir is .false., atom, dir, atom2 and dir2 are used
    logical :: pure_dir
    integer            :: atom1, atom2
    FLOAT, allocatable :: mix1(:,:) !< mix1(natoms, ndim)
    FLOAT, allocatable :: mix2(:,:)

  contains
    procedure :: copy_to => perturbation_ionic_copy
    generic   :: assignment(=) => copy_to
    procedure :: info => perturbation_ionic_info
    procedure :: dapply => dperturbation_ionic_apply
    procedure :: zapply => zperturbation_ionic_apply
    procedure :: dapply_order_2 => dperturbation_ionic_apply_order_2
    procedure :: zapply_order_2 => zperturbation_ionic_apply_order_2
    procedure :: setup_dir => perturbation_ionic_setup_dir
    procedure :: setup_atom => perturbation_ionic_setup_atom
    procedure :: setup_mix_dir => perturbation_ionic_setup_mixed_dir
    final :: perturbation_ionic_finalize
  end type perturbation_ionic_t

  interface perturbation_ionic_t
    procedure perturbation_ionic_constructor
  end interface perturbation_ionic_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function perturbation_ionic_constructor(namespace, ions) result(pert)
    class(perturbation_ionic_t), pointer    :: pert
    type(namespace_t),           intent(in) :: namespace
    type(ions_t),  target,       intent(in) :: ions

    PUSH_SUB(perturbation_ionic_constructor)

    SAFE_ALLOCATE(pert)

    call perturbation_ionic_init(pert, namespace, ions)

    POP_SUB(perturbation_ionic_constructor)
  end function perturbation_ionic_constructor


  ! --------------------------------------------------------------------
  subroutine perturbation_ionic_init(this, namespace, ions)
    type(perturbation_ionic_t),  intent(out) :: this
    type(namespace_t),           intent(in)  :: namespace
    type(ions_t), target,        intent(in)  :: ions

    PUSH_SUB(perturbation_ionic_init)

    this%dir = -1
    this%dir2 = -1
    this%atom1 = -1
    this%atom2 = -1

    this%ions => ions

    this%pure_dir = .false.

    SAFE_ALLOCATE(this%mix1(1:ions%natoms, 1:ions%space%dim))
    SAFE_ALLOCATE(this%mix2(1:ions%natoms, 1:ions%space%dim))

    POP_SUB(perturbation_ionic_init)
  end subroutine perturbation_ionic_init

  ! --------------------------------------------------------------------
  subroutine perturbation_ionic_copy(this, source)
    class(perturbation_ionic_t), intent(out) :: this
    class(perturbation_ionic_t), intent(in)  :: source

    PUSH_SUB(perturbation_ionic_copy)

    call perturbation_copy(this, source)
    this%atom1 = source%atom1
    this%atom2 = source%atom2
    this%ions => source%ions

    this%pure_dir = source%pure_dir

    SAFE_ALLOCATE(this%mix1(1:source%ions%natoms, 1:source%ions%space%dim))
    SAFE_ALLOCATE(this%mix2(1:source%ions%natoms, 1:source%ions%space%dim))
    this%mix1 = source%mix1
    this%mix2 = source%mix2

    POP_SUB(perturbation_ionic_copy)
  end subroutine perturbation_ionic_copy


  ! --------------------------------------------------------------------
  subroutine perturbation_ionic_finalize(this)
    type(perturbation_ionic_t), intent(inout) :: this

    PUSH_SUB(perturbation_ionic_finalize)

    SAFE_DEALLOCATE_A(this%mix1)
    SAFE_DEALLOCATE_A(this%mix2)

    POP_SUB(perturbation_ionic_finalize)
  end subroutine perturbation_ionic_finalize

  ! --------------------------------------------------------------------
  subroutine perturbation_ionic_info(this)
    class(perturbation_ionic_t), intent(in) :: this

    PUSH_SUB(perturbation_ionic_info)

    POP_SUB(perturbation_ionic_info)
  end subroutine perturbation_ionic_info

  ! --------------------------------------------------------------------
  subroutine perturbation_ionic_setup_dir(this, dir, dir2)
    class(perturbation_ionic_t),   intent(inout) :: this
    integer,                       intent(in)    :: dir
    integer, optional,             intent(in)    :: dir2

    PUSH_SUB(perturbation_ionic_setup_dir)

    this%dir = dir
    if (present(dir2)) this%dir2 = dir2

    this%pure_dir = .true.

    this%mix1 = M_ZERO
    this%mix2 = M_ZERO

    if (this%dir  > 0 .and. this%atom1 > 0) this%mix1(this%atom1, this%dir ) = M_ONE
    if (this%dir2 > 0 .and. this%atom2 > 0) this%mix2(this%atom2, this%dir2) = M_ONE

    POP_SUB(perturbation_ionic_setup_dir)
  end subroutine perturbation_ionic_setup_dir

  ! --------------------------------------------------------------------
  subroutine perturbation_ionic_setup_atom(this, iatom, iatom2)
    class(perturbation_ionic_t),   intent(inout) :: this
    integer,                       intent(in)    :: iatom
    integer, optional,             intent(in)    :: iatom2

    PUSH_SUB(perturbation_ionic_setup_atom)

    this%atom1 = iatom
    if (present(iatom2)) this%atom2 = iatom2

    this%pure_dir = .true.

    this%mix1 = M_ZERO
    this%mix2 = M_ZERO

    if (this%dir  > 0 .and. this%atom1 > 0) this%mix1(this%atom1, this%dir ) = M_ONE
    if (this%dir2 > 0 .and. this%atom2 > 0) this%mix2(this%atom2, this%dir2) = M_ONE

    POP_SUB(perturbation_ionic_setup_atom)
  end subroutine perturbation_ionic_setup_atom

  ! --------------------------------------------------------------------
  subroutine perturbation_ionic_setup_mixed_dir(this, iatom, idir, val, jatom, jdir, valuej)
    class(perturbation_ionic_t),      intent(inout) :: this
    integer,           intent(in)    :: iatom
    integer,           intent(in)    :: idir
    FLOAT,             intent(in)    :: val
    integer, optional, intent(in)    :: jatom
    integer, optional, intent(in)    :: jdir
    FLOAT,   optional, intent(in)    :: valuej

    logical :: have_dir_2

    PUSH_SUB(perturbation_ionic_setup_mixed_dir)

    this%pure_dir = .false.

    this%mix1(iatom, idir) = val

    have_dir_2 = present(jatom) .and. present(jdir) .and. present(jatom)

    if (have_dir_2) then
      this%mix1(jatom, jdir) = valuej
    else
      ASSERT(.not. present(jatom) .and. .not. present(jdir) .and. .not. present(jatom))
    end if

    POP_SUB(perturbation_ionic_setup_mixed_dir)
  end subroutine perturbation_ionic_setup_mixed_dir



#include "undef.F90"
#include "real.F90"
#include "perturbation_ionic_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "perturbation_ionic_inc.F90"

end module perturbation_ionic_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
