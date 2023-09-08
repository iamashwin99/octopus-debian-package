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

module perturbation_magnetic_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
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
  use parser_oct_m
  use perturbation_oct_m
  use physics_op_oct_m
  use profiling_oct_m
  use projector_oct_m
  use space_oct_m
  use species_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use varinfo_oct_m
  use vibrations_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                           &
    perturbation_magnetic_t

  integer, public, parameter :: &
    GAUGE_GIPAW  = 1, &
    GAUGE_ICL    = 2

  type, extends(perturbation_t) :: perturbation_magnetic_t
    private
    integer            :: gauge
    type(ions_t), pointer :: ions => null()
  contains
    procedure :: copy_to => perturbation_magnetic_copy
    generic   :: assignment(=) => copy_to
    procedure :: info => perturbation_magnetic_info
    procedure :: dapply => dperturbation_magnetic_apply
    procedure :: zapply => zperturbation_magnetic_apply
    procedure :: dapply_order_2 => dperturbation_magnetic_apply_order_2
    procedure :: zapply_order_2 => zperturbation_magnetic_apply_order_2
    final :: perturbation_magnetic_finalize
  end type perturbation_magnetic_t

  interface perturbation_magnetic_t
    procedure perturbation_magnetic_constructor
  end interface perturbation_magnetic_t


contains
  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function perturbation_magnetic_constructor(namespace, ions) result(pert)
    class(perturbation_magnetic_t), pointer    :: pert
    type(namespace_t),              intent(in) :: namespace
    type(ions_t),  target,          intent(in) :: ions

    PUSH_SUB(perturbation_magnetic_constructor)

    SAFE_ALLOCATE(pert)

    call perturbation_magnetic_init(pert, namespace, ions)

    POP_SUB(perturbation_magnetic_constructor)
  end function perturbation_magnetic_constructor


  ! --------------------------------------------------------------------
  subroutine perturbation_magnetic_init(this, namespace, ions)
    type(perturbation_magnetic_t),    intent(out) :: this
    type(namespace_t),                intent(in)  :: namespace
    type(ions_t), target,              intent(in) :: ions

    PUSH_SUB(perturbation_magnetic_init)

    this%dir = -1
    this%dir2 = -1

    this%ions => ions

    !%Variable MagneticGaugeCorrection
    !%Type integer
    !%Default gipaw
    !%Section Linear Response
    !%Description
    !% For magnetic linear response: how to handle gauge-invariance in the description
    !% of the coupling of electrons to the magnetic field.
    !%Option none 0
    !% No correction.
    !%Option gipaw 1
    !% GIPAW correction: C Pickard and F Mauri, <i>Phys. Rev. Lett.</i> <b>91</b>, 196401 (2003).
    !%Option icl 2
    !% ICL correction: S Ismail-Beigi, EK Chang, and SG Louie, <i>Phys. Rev. Lett.</i> <b>87</b>, 087402 (2001).
    !%End

    call parse_variable(namespace, 'MagneticGaugeCorrection', GAUGE_GIPAW, this%gauge)
    if (.not. varinfo_valid_option('MagneticGaugeCorrection', this%gauge)) then
      call messages_input_error(namespace, 'MagneticGaugeCorrection')
    end if

    POP_SUB(perturbation_magnetic_init)
  end subroutine perturbation_magnetic_init

  ! --------------------------------------------------------------------
  subroutine perturbation_magnetic_finalize(this)
    type(perturbation_magnetic_t), intent(inout) :: this

    PUSH_SUB(perturbation_magnetic_finalize)

    POP_SUB(perturbation_magnetic_finalize)
  end subroutine perturbation_magnetic_finalize

  ! --------------------------------------------------------------------
  subroutine perturbation_magnetic_copy(this, source)
    class(perturbation_magnetic_t), intent(out) :: this
    class(perturbation_magnetic_t), intent(in)  :: source

    PUSH_SUB(perturbation_magnetic_copy)

    call perturbation_copy(this, source)
    this%ions => source%ions

    this%gauge = source%gauge

    POP_SUB(perturbation_magnetic_copy)
  end subroutine perturbation_magnetic_copy

  ! --------------------------------------------------------------------
  subroutine perturbation_magnetic_info(this)
    class(perturbation_magnetic_t), intent(in) :: this

    PUSH_SUB(perturbation_magnetic_info)

    POP_SUB(perturbation_magnetic_info)
  end subroutine perturbation_magnetic_info


#include "undef.F90"
#include "real.F90"
#include "perturbation_magnetic_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "perturbation_magnetic_inc.F90"

end module perturbation_magnetic_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
