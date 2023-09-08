!! Copyright (C) 2008 X. Andrade
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

module batch_ops_oct_m
  use accel_oct_m
  use batch_oct_m
  use blas_oct_m
  use debug_oct_m
  use iso_c_binding
  use global_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use messages_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private
  public ::                         &
    batch_set_zero,                 &
    batch_axpy,                     &
    batch_scal,                     &
    batch_xpay,                     &
    batch_set_state,                &
    batch_get_state,                &
    batch_get_points,               &
    batch_set_points,               &
    batch_points_block_size,        &
    batch_mul,                      &
    batch_add_with_map,             &
    batch_copy_with_map,            &
    dbatch_axpy_function,           &
    zbatch_axpy_function,           &
    dbatch_ax_function_py,          &
    zbatch_ax_function_py,          &
    batch_split_complex

  interface batch_axpy
    module procedure dbatch_axpy_const
    module procedure zbatch_axpy_const
    module procedure dbatch_axpy_vec
    module procedure zbatch_axpy_vec
  end interface batch_axpy

  interface batch_scal
    module procedure dbatch_scal_const
    module procedure zbatch_scal_const
    module procedure dbatch_scal_vec
    module procedure zbatch_scal_vec
  end interface batch_scal

  interface batch_xpay
    module procedure dbatch_xpay_vec
    module procedure zbatch_xpay_vec
    module procedure dbatch_xpay_const
    module procedure zbatch_xpay_const
  end interface batch_xpay

  interface batch_add_with_map
    module procedure batch_add_with_map_cpu
    module procedure batch_add_with_map_cl
  end interface batch_add_with_map

  interface batch_copy_with_map
    module procedure batch_copy_with_map_cpu
    module procedure batch_copy_with_map_cl
  end interface batch_copy_with_map

  ! There are several ways how to call batch_set_state and batch_get_state:
  ! 1. With a 1d array and a single index: batch_get_state(psib, ist, np, psi)
  !    In this case, ist is between 1 and nst_linear of the batch (i.e. it is the
  !    linear index in the batch) and this call will fetch the first np points of
  !    that part of the batch.
  ! 2. With a 1d array and an index tuple: batch_get_state(psib, (/ist, idim/), np, psi)
  !    In this case, ist corresponds to the real state index in the system, i.e., it
  !    should be between states_elec_block_min(st, ib) and states_elec_block_max(st, ib)
  !    for the batch with index ib. idim gives the spin space index.
  ! 3. With a 2d array and a single index: batch_get_state(psib, ist, np, psi(:, :))
  !    In this case, ist is between 1 and nst of the batch (i.e. the state index of
  !    the batch) and this call will fetch the wavefunctions for all spin indices for
  !    that state index. This is more efficient than looping over spin space and
  !    getting the wavefunctions for different spin indices individually.
  interface batch_set_state
    module procedure dbatch_set_state1
    module procedure zbatch_set_state1
    module procedure dbatch_set_state2
    module procedure zbatch_set_state2
    module procedure dbatch_set_state3
    module procedure zbatch_set_state3
  end interface batch_set_state

  interface batch_get_state
    module procedure dbatch_get_state1
    module procedure zbatch_get_state1
    module procedure dbatch_get_state2
    module procedure zbatch_get_state2
    module procedure dbatch_get_state3
    module procedure zbatch_get_state3
  end interface batch_get_state

  interface batch_get_points
    module procedure dbatch_get_points
    module procedure zbatch_get_points
    module procedure batch_get_points_cl
  end interface batch_get_points

  interface batch_set_points
    module procedure dbatch_set_points
    module procedure zbatch_set_points
    module procedure batch_set_points_cl
  end interface batch_set_points

  interface batch_mul
    module procedure dbatch_mul
    module procedure zbatch_mul
  end interface batch_mul

  type(profile_t), save :: get_points_prof, set_points_prof

contains

  !--------------------------------------------------------------

  subroutine batch_set_zero(this)
    class(batch_t),     intent(inout) :: this

    type(profile_t), save :: prof
    integer :: ist_linear, ist, ip

    PUSH_SUB(batch_set_zero)

    ASSERT(not_in_openmp())

    call profiling_in(prof, "BATCH_SET_ZERO")

    select case (this%status())
    case (BATCH_DEVICE_PACKED)
      call accel_set_buffer_to_zero(this%ff_device, this%type(), product(this%pack_size))

    case (BATCH_PACKED)
      if (this%type() == TYPE_FLOAT) then
        !$omp parallel do private(ist) schedule(static)
        do ip = 1, int(this%pack_size(2), i4)
          !$omp simd
          do ist = 1, int(this%pack_size(1), i4)
            this%dff_pack(ist, ip) = M_ZERO
          end do
        end do
      else
        !$omp parallel do private(ist) schedule(static)
        do ip = 1, int(this%pack_size(2), i4)
          !$omp simd
          do ist = 1, int(this%pack_size(1), i4)
            this%zff_pack(ist, ip) = M_z0
          end do
        end do
      end if

    case (BATCH_NOT_PACKED)
      if (this%type() == TYPE_FLOAT) then
        do ist_linear = 1, this%nst_linear
          !$omp parallel do schedule(static)
          do ip = 1, ubound(this%dff_linear, dim=1)
            this%dff_linear(ip, ist_linear) = M_ZERO
          end do
        end do
      else
        do ist_linear = 1, this%nst_linear
          !$omp parallel do schedule(static)
          do ip = 1, ubound(this%zff_linear, dim=1)
            this%zff_linear(ip, ist_linear) = M_z0
          end do
        end do
      end if

    case default
      ASSERT(.false.)

    end select

    call profiling_out(prof)

    POP_SUB(batch_set_zero)
  end subroutine batch_set_zero

! --------------------------------------------------------------

  subroutine batch_get_points_cl(this, sp, ep, psi, ldpsi1, ldpsi2)
    class(batch_t),      intent(in)    :: this
    integer,             intent(in)    :: sp
    integer,             intent(in)    :: ep
    type(accel_mem_t),   intent(inout) :: psi
    integer,             intent(in)    :: ldpsi1
    integer,             intent(in)    :: ldpsi2

    integer :: tsize, ii, it
    type(accel_kernel_t), save :: kernel
    integer, allocatable :: linear_to_ist(:), linear_to_idim(:)
    type(accel_mem_t) :: buff_linear_to_ist, buff_linear_to_idim

    PUSH_SUB(batch_get_points_cl)
    call profiling_in(get_points_prof, "GET_POINTS")

    select case (this%status())
    case (BATCH_NOT_PACKED, BATCH_PACKED)
      call messages_not_implemented('batch_get_points_cl for non-CL batches')

    case (BATCH_DEVICE_PACKED)

      tsize = types_get_size(this%type())/types_get_size(TYPE_FLOAT)
      SAFE_ALLOCATE(linear_to_ist(1:this%nst_linear*tsize))
      SAFE_ALLOCATE(linear_to_idim(1:this%nst_linear*tsize))
      do ii = 1, this%nst_linear
        do it = 1, tsize
          linear_to_ist(tsize*(ii-1)+it) = tsize*(this%linear_to_ist(ii) - 1) + it - 1
          linear_to_idim(tsize*(ii-1)+it) = this%linear_to_idim(ii) - 1
        end do
      end do

      call accel_create_buffer(buff_linear_to_ist, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nst_linear*tsize)
      call accel_write_buffer(buff_linear_to_ist, this%nst_linear*tsize, linear_to_ist)
      call accel_create_buffer(buff_linear_to_idim, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nst_linear*tsize)
      call accel_write_buffer(buff_linear_to_idim, this%nst_linear*tsize, linear_to_idim)

      call accel_kernel_start_call(kernel, 'points.cl', 'get_points')

      call accel_set_kernel_arg(kernel, 0, sp)
      call accel_set_kernel_arg(kernel, 1, ep)
      call accel_set_kernel_arg(kernel, 2, buff_linear_to_ist)
      call accel_set_kernel_arg(kernel, 3, buff_linear_to_idim)
      call accel_set_kernel_arg(kernel, 4, this%nst_linear*tsize)
      call accel_set_kernel_arg(kernel, 5, this%ff_device)
      call accel_set_kernel_arg(kernel, 6, int(this%pack_size_real(1), i4))
      call accel_set_kernel_arg(kernel, 7, psi)
      call accel_set_kernel_arg(kernel, 8, ldpsi1*tsize)
      call accel_set_kernel_arg(kernel, 9, ldpsi2)

      call accel_kernel_run(kernel, (/this%pack_size_real(1), int(ep - sp + 1, i8)/), (/this%pack_size_real(1), 1_i8/))

      call accel_release_buffer(buff_linear_to_ist)
      call accel_release_buffer(buff_linear_to_idim)
      SAFE_DEALLOCATE_A(linear_to_ist)
      SAFE_DEALLOCATE_A(linear_to_idim)

    end select

    call profiling_out(get_points_prof)

    POP_SUB(batch_get_points_cl)
  end subroutine batch_get_points_cl

! --------------------------------------------------------------

  subroutine batch_set_points_cl(this, sp, ep, psi, ldpsi1, ldpsi2)
    class(batch_t),      intent(inout) :: this
    integer,             intent(in)    :: sp
    integer,             intent(in)    :: ep
    type(accel_mem_t),   intent(in)    :: psi
    integer,             intent(in)    :: ldpsi1
    integer,             intent(in)    :: ldpsi2

    integer :: tsize, ii, it
    type(accel_kernel_t), save :: kernel
    integer, allocatable :: linear_to_ist(:), linear_to_idim(:)
    type(accel_mem_t) :: buff_linear_to_ist, buff_linear_to_idim

    PUSH_SUB(batch_set_points_cl)
    call profiling_in(set_points_prof, "SET_POINTS")

    select case (this%status())
    case (BATCH_NOT_PACKED, BATCH_PACKED)
      call messages_not_implemented('batch_set_points_cl for non-CL batches')

    case (BATCH_DEVICE_PACKED)

      tsize = types_get_size(this%type())/types_get_size(TYPE_FLOAT)
      SAFE_ALLOCATE(linear_to_ist(1:this%nst_linear*tsize))
      SAFE_ALLOCATE(linear_to_idim(1:this%nst_linear*tsize))
      do ii = 1, this%nst_linear
        do it = 1, tsize
          linear_to_ist(tsize*(ii-1)+it) = tsize*(this%linear_to_ist(ii) - 1) + it - 1
          linear_to_idim(tsize*(ii-1)+it) = this%linear_to_idim(ii) - 1
        end do
      end do

      call accel_create_buffer(buff_linear_to_ist, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nst_linear*tsize)
      call accel_write_buffer(buff_linear_to_ist, this%nst_linear*tsize, linear_to_ist)
      call accel_create_buffer(buff_linear_to_idim, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nst_linear*tsize)
      call accel_write_buffer(buff_linear_to_idim, this%nst_linear*tsize, linear_to_idim)

      call accel_kernel_start_call(kernel, 'points.cl', 'set_points')

      call accel_set_kernel_arg(kernel, 0, sp)
      call accel_set_kernel_arg(kernel, 1, ep)
      call accel_set_kernel_arg(kernel, 2, buff_linear_to_ist)
      call accel_set_kernel_arg(kernel, 3, buff_linear_to_idim)
      call accel_set_kernel_arg(kernel, 4, this%nst_linear*tsize)
      call accel_set_kernel_arg(kernel, 5, psi)
      call accel_set_kernel_arg(kernel, 6, ldpsi1*tsize)
      call accel_set_kernel_arg(kernel, 7, ldpsi2)
      call accel_set_kernel_arg(kernel, 8, this%ff_device)
      call accel_set_kernel_arg(kernel, 9, int(this%pack_size_real(1), i4))

      call accel_kernel_run(kernel, (/this%pack_size_real(1), int(ep - sp + 1, i8)/), (/this%pack_size_real(1), 1_i8/))

      call accel_release_buffer(buff_linear_to_ist)
      call accel_release_buffer(buff_linear_to_idim)
      SAFE_DEALLOCATE_A(linear_to_ist)
      SAFE_DEALLOCATE_A(linear_to_idim)

    end select

    call profiling_out(set_points_prof)

    POP_SUB(batch_set_points_cl)
  end subroutine batch_set_points_cl

! -------------------------

  integer pure function batch_points_block_size() result(block_size)

    block_size = 61440

  end function batch_points_block_size

! -------------------------
  subroutine batch_add_with_map_cpu(np, map, xx, yy, zz)
    integer,           intent(in)    :: np
    integer,           intent(in)    :: map(:)
    class(batch_t),    intent(in)    :: xx
    class(batch_t),    intent(in)    :: yy
    class(batch_t),    intent(inout) :: zz
    type(accel_mem_t) :: buff_map

    PUSH_SUB(batch_add_with_map_cpu)

    if (xx%status() /= BATCH_DEVICE_PACKED) then
      if (xx%type() == TYPE_FLOAT) then
        call dbatch_add_with_map(np, map, xx, yy, zz)
      else
        call zbatch_add_with_map(np, map, xx, yy, zz)
      end if
    else
      ! copy map to GPU if not already there
      call accel_create_buffer(buff_map, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, np)
      call accel_write_buffer(buff_map, np, map)
      call batch_add_with_map_cl(np, buff_map, xx, yy, zz)
      call accel_release_buffer(buff_map)
    end if

    POP_SUB(batch_add_with_map_cpu)
  end subroutine batch_add_with_map_cpu

! -------------------------
  subroutine batch_add_with_map_cl(np, map, xx, yy, zz)
    integer,            intent(in)    :: np
    class(accel_mem_t), intent(in)    :: map
    class(batch_t),     intent(in)    :: xx
    class(batch_t),     intent(in)    :: yy
    class(batch_t),     intent(inout) :: zz

    type(accel_kernel_t), save :: kernel
    integer(i8) :: localsize, dim3, dim2

    PUSH_SUB(batch_add_with_map_cl)

    call accel_kernel_start_call(kernel, 'copy.cl', 'add_with_map')

    call accel_set_kernel_arg(kernel, 0, np)
    call accel_set_kernel_arg(kernel, 1, map)
    call accel_set_kernel_arg(kernel, 2, xx%ff_device)
    call accel_set_kernel_arg(kernel, 3, log2(int(xx%pack_size_real(1), i4)))
    call accel_set_kernel_arg(kernel, 4, yy%ff_device)
    call accel_set_kernel_arg(kernel, 5, log2(int(yy%pack_size_real(1), i4)))
    call accel_set_kernel_arg(kernel, 6, zz%ff_device)
    call accel_set_kernel_arg(kernel, 7, log2(int(zz%pack_size_real(1), i4)))

    localsize = accel_kernel_workgroup_size(kernel)/xx%pack_size_real(1)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))

    call accel_kernel_run(kernel, (/xx%pack_size_real(1), dim2, dim3/), (/xx%pack_size_real(1), localsize, 1_i8/))

    POP_SUB(batch_add_with_map_cl)
  end subroutine batch_add_with_map_cl

! -------------------------
  subroutine batch_copy_with_map_cpu(np, map, xx, yy)
    integer,           intent(in)    :: np
    integer,           intent(in)    :: map(:)
    class(batch_t),    intent(in)    :: xx
    class(batch_t),    intent(inout) :: yy
    type(accel_mem_t) :: buff_map

    PUSH_SUB(batch_copy_with_map_cpu)

    if (xx%status() /= BATCH_DEVICE_PACKED) then
      if (xx%type() == TYPE_FLOAT) then
        call dbatch_copy_with_map(np, map, xx, yy)
      else
        call zbatch_copy_with_map(np, map, xx, yy)
      end if
    else
      ! copy map to GPU if not already there
      call accel_create_buffer(buff_map, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, np)
      call accel_write_buffer(buff_map, np, map)
      call batch_copy_with_map_cl(np, buff_map, xx, yy)
      call accel_release_buffer(buff_map)
    end if

    POP_SUB(batch_copy_with_map_cpu)
  end subroutine batch_copy_with_map_cpu

! -------------------------
  subroutine batch_copy_with_map_cl(np, map, xx, yy)
    integer,            intent(in)    :: np
    class(accel_mem_t), intent(in)    :: map
    class(batch_t),     intent(in)    :: xx
    class(batch_t),     intent(inout) :: yy

    type(accel_kernel_t), save :: kernel
    integer(i8) :: localsize, dim3, dim2

    PUSH_SUB(batch_copy_with_map_cl)

    call accel_kernel_start_call(kernel, 'copy.cl', 'copy_with_map')

    call accel_set_kernel_arg(kernel, 0, np)
    call accel_set_kernel_arg(kernel, 1, map)
    call accel_set_kernel_arg(kernel, 2, xx%ff_device)
    call accel_set_kernel_arg(kernel, 3, log2(int(xx%pack_size_real(1), i4)))
    call accel_set_kernel_arg(kernel, 4, yy%ff_device)
    call accel_set_kernel_arg(kernel, 5, log2(int(yy%pack_size_real(1), i4)))

    localsize = accel_kernel_workgroup_size(kernel)/xx%pack_size_real(1)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))

    call accel_kernel_run(kernel, (/xx%pack_size_real(1), dim2, dim3/), (/xx%pack_size_real(1), localsize, 1_i8/))

    POP_SUB(batch_copy_with_map_cl)
  end subroutine batch_copy_with_map_cl

! -------------------------
  subroutine batch_split_complex(np, xx, yy, zz)
    integer,           intent(in)    :: np
    class(batch_t),    intent(in)    :: xx
    class(batch_t),    intent(inout) :: yy
    class(batch_t),    intent(inout) :: zz

    integer :: ist_linear, ip
    type(accel_kernel_t), save :: kernel
    integer(i8) :: localsize, dim3, dim2

    PUSH_SUB(batch_split_complex)

    ASSERT(xx%type() == TYPE_CMPLX)
    ASSERT(yy%type() == TYPE_FLOAT)
    ASSERT(zz%type() == TYPE_FLOAT)
    ASSERT(xx%status() == yy%status())
    ASSERT(xx%status() == zz%status())

    select case (xx%status())
    case (BATCH_NOT_PACKED)
      do ist_linear = 1, xx%nst_linear
        !$omp parallel do schedule(static)
        do ip = 1, np
          yy%dff_linear(ip, ist_linear) = TOFLOAT(xx%zff_linear(ip, ist_linear))
          zz%dff_linear(ip, ist_linear) = aimag(xx%zff_linear(ip, ist_linear))
        end do
      end do
    case (BATCH_PACKED)
      !$omp parallel do private(ist_linear) schedule(static)
      do ip = 1, np
        do ist_linear = 1, xx%nst_linear
          yy%dff_pack(ist_linear, ip) = TOFLOAT(xx%zff_pack(ist_linear, ip))
          zz%dff_pack(ist_linear, ip) = aimag(xx%zff_pack(ist_linear, ip))
        end do
      end do
    case (BATCH_DEVICE_PACKED)
      call accel_kernel_start_call(kernel, 'split.cl', 'split_complex')

      call accel_set_kernel_arg(kernel, 0, int(xx%pack_size(2), i4))
      call accel_set_kernel_arg(kernel, 1, xx%ff_device)
      call accel_set_kernel_arg(kernel, 2, log2(int(xx%pack_size(1), i4)))
      call accel_set_kernel_arg(kernel, 3, yy%ff_device)
      call accel_set_kernel_arg(kernel, 4, log2(int(yy%pack_size(1), i4)))
      call accel_set_kernel_arg(kernel, 5, zz%ff_device)
      call accel_set_kernel_arg(kernel, 6, log2(int(zz%pack_size(1), i4)))

      localsize = accel_kernel_workgroup_size(kernel)/xx%pack_size(1)

      dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))

      call accel_kernel_run(kernel, (/xx%pack_size(1), dim2, dim3/), (/xx%pack_size(1), localsize, 1_i8/))
    end select

    POP_SUB(batch_split_complex)
  end subroutine batch_split_complex

#include "real.F90"
#include "batch_ops_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "batch_ops_inc.F90"
#include "undef.F90"

end module batch_ops_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
