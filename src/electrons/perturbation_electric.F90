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

module perturbation_electric_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use perturbation_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use types_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                           &
    perturbation_electric_t,          &
    perturbation_electric_init

  type, extends(perturbation_t) :: perturbation_electric_t
    private
  contains
    procedure :: copy => perturbation_electric_copy
    generic   :: assignment(=) => copy
    procedure :: info => perturbation_electric_info
    procedure :: dapply => dperturbation_electric_apply
    procedure :: zapply => zperturbation_electric_apply
    procedure :: apply_batch => perturbation_electric_apply_batch
    procedure :: dapply_order_2 => dperturbation_electric_apply_order_2
    procedure :: zapply_order_2 => zperturbation_electric_apply_order_2
    final :: perturbation_electric_finalize
  end type perturbation_electric_t

  interface perturbation_electric_t
    procedure perturbation_electric_constructor
  end interface perturbation_electric_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function perturbation_electric_constructor(namespace) result(pert)
    class(perturbation_electric_t), pointer    :: pert
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(perturbation_electric_constructor)

    SAFE_ALLOCATE(pert)

    call perturbation_electric_init(pert, namespace)

    POP_SUB(perturbation_electric_constructor)
  end function perturbation_electric_constructor

  ! --------------------------------------------------------------------
  subroutine perturbation_electric_init(this, namespace)
    type(perturbation_electric_t),   intent(out) :: this
    type(namespace_t),               intent(in)  :: namespace

    PUSH_SUB(perturbation_electric_init)

    this%dir = -1
    this%dir2 = -1

    POP_SUB(perturbation_electric_init)
  end subroutine perturbation_electric_init

  ! --------------------------------------------------------------------
  subroutine perturbation_electric_finalize(this)
    type(perturbation_electric_t), intent(inout) :: this

    PUSH_SUB(perturbation_electric_finalize)

    POP_SUB(perturbation_electric_finalize)
  end subroutine perturbation_electric_finalize

  ! --------------------------------------------------------------------
  subroutine perturbation_electric_copy(this, source)
    class(perturbation_electric_t), intent(out) :: this
    class(perturbation_electric_t), intent(in)  :: source

    PUSH_SUB(perturbation_electric_copy)

    call perturbation_copy(this, source)

    POP_SUB(perturbation_electric_copy)
  end subroutine perturbation_electric_copy


  ! --------------------------------------------------------------------
  subroutine perturbation_electric_info(this)
    class(perturbation_electric_t), intent(in) :: this

    PUSH_SUB(perturbation_electric_info)

    POP_SUB(perturbation_electric_info)
  end subroutine perturbation_electric_info

  ! --------------------------------------------------------------------------
  subroutine perturbation_electric_apply_batch(this, namespace, space, gr, hm, f_in, f_out)
    class(perturbation_electric_t),    intent(in)    :: this
    type(namespace_t),                 intent(in)    :: namespace
    type(space_t),                     intent(in)    :: space
    type(grid_t),                      intent(in)    :: gr
    type(hamiltonian_elec_t),          intent(in)    :: hm
    type(wfs_elec_t),                  intent(in)    :: f_in
    type(wfs_elec_t),                  intent(inout) :: f_out

    integer :: ii, ip

    PUSH_SUB(perturbation_electric_apply_batch)

    ! electric does not need it since (e^-ikr)r(e^ikr) = r

    ASSERT(f_in%status() == f_out%status())

    select case (f_in%status())

    case (BATCH_NOT_PACKED)

      if (f_in%type() == TYPE_FLOAT) then

        do ii = 1, f_in%nst_linear
          !$omp parallel do private(ip)
          do ip = 1, gr%np
            f_out%dff_linear(ip, ii) = gr%x(ip, this%dir)*f_in%dff_linear(ip, ii)
          end do
        end do

      else

        do ii = 1, f_in%nst_linear
          !$omp parallel do private(ip)
          do ip = 1, gr%np
            f_out%zff_linear(ip, ii) = gr%x(ip, this%dir)*f_in%zff_linear(ip, ii)
          end do
        end do

      end if

    case (BATCH_PACKED)

      if (f_in%type() == TYPE_FLOAT) then

        !$omp parallel do private(ip, ii)
        do ip = 1, gr%np
          do ii = 1, f_in%nst_linear
            f_out%dff_pack(ii, ip) = gr%x(ip, this%dir)*f_in%dff_pack(ii, ip)
          end do
        end do

      else

        !$omp parallel do private(ip, ii)
        do ip = 1, gr%np
          do ii = 1, f_in%nst_linear
            f_out%zff_pack(ii, ip) = gr%x(ip, this%dir)*f_in%zff_pack(ii, ip)
          end do
        end do

      end if

    case (BATCH_DEVICE_PACKED)

      call messages_not_implemented("perturbation_electric_apply_batch on GPUs", namespace = namespace)

    end select

    POP_SUB(perturbation_electric_apply_batch)
  end subroutine perturbation_electric_apply_batch


#include "undef.F90"
#include "real.F90"
#include "perturbation_electric_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "perturbation_electric_inc.F90"

end module perturbation_electric_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
