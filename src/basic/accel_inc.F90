!! Copyright (C) 2010 X. Andrade
!! Copyright (C) 2022 S. Ohlmann
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

subroutine X(accel_write_buffer_single)(this, data, async)
  type(accel_mem_t),               intent(inout) :: this
  R_TYPE,                          intent(in)    :: data
  logical,               optional, intent(in)    :: async

  PUSH_SUB(X(accel_write_buffer_single))

  call X(accel_write_buffer_0)(this, 1_i8, data, async=async)

  POP_SUB(X(accel_write_buffer_single))
end subroutine X(accel_write_buffer_single)


subroutine X(accel_write_buffer_0)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE,                           intent(in)    :: data
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(i8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif

  PUSH_SUB(X(accel_write_buffer_0))
  call profiling_in(prof_write, TOSTRING(X(CL_WRITE_BUFFER)))

  ! it does not make sense to write a buffer that the kernels cannot read
  ASSERT(this%flags /= ACCEL_MEM_WRITE_ONLY)

  fsize = size*R_SIZEOF
  offset_ = 0_i8
  if (present(offset)) offset_ = offset*R_SIZEOF

  async_ = optional_default(async, .false.)

  if (fsize > 0) then
#ifdef HAVE_OPENCL
    call clEnqueueWriteBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data, ierr)
    if (ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueWriteBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_htod(this%mem, data, fsize, offset_)
#endif

    call profiling_count_transfers(size, data)
    if (.not. async_) call accel_finish()
  end if

  call profiling_out(prof_write)
  POP_SUB(X(accel_write_buffer_0))
end subroutine X(accel_write_buffer_0)


subroutine X(accel_write_buffer_1)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE, contiguous,               intent(in)    :: data(:)
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  if (size == 0) return

  PUSH_SUB(X(accel_write_buffer_1))

  ASSERT(ubound(data, dim=1) >= size)

  call X(accel_write_buffer_0)(this, size, data(1), offset, async)

  POP_SUB(X(accel_write_buffer_1))
end subroutine X(accel_write_buffer_1)


subroutine X(accel_write_buffer_2)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE, contiguous,               intent(in)    :: data(:, :)
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  if (size == 0) return

  PUSH_SUB(X(accel_write_buffer_2))

  ASSERT(ubound(data, dim=1) >= 1)
  ASSERT(ubound(data, dim=2) >= 1)

  call X(accel_write_buffer_0)(this, size, data(1, 1), offset, async)

  POP_SUB(X(accel_write_buffer_2))
end subroutine X(accel_write_buffer_2)


subroutine X(accel_write_buffer_3)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE, contiguous,               intent(in)    :: data(:, :, :)
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  if (size == 0) return

  PUSH_SUB(X(accel_write_buffer_3))

  ASSERT(ubound(data, dim=1) >= 1)
  ASSERT(ubound(data, dim=2) >= 1)
  ASSERT(ubound(data, dim=3) >= 1)

  call X(accel_write_buffer_0)(this, size, data(1, 1, 1), offset, async)

  POP_SUB(X(accel_write_buffer_3))
end subroutine X(accel_write_buffer_3)


subroutine X(accel_write_buffer_0_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(in)    :: data
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_write_buffer_0_i4))

  if (present(offset)) then
    call accel_write_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_write_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_write_buffer_0_i4))
end subroutine X(accel_write_buffer_0_i4)


subroutine X(accel_write_buffer_1_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE, contiguous,               intent(in)    :: data(:)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_write_buffer_1_i4))

  if (present(offset)) then
    call accel_write_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_write_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_write_buffer_1_i4))
end subroutine X(accel_write_buffer_1_i4)


subroutine X(accel_write_buffer_2_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE, contiguous,               intent(in)    :: data(:, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_write_buffer_2_i4))

  if (present(offset)) then
    call accel_write_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_write_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_write_buffer_2_i4))
end subroutine X(accel_write_buffer_2_i4)


subroutine X(accel_write_buffer_3_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE, contiguous,               intent(in)    :: data(:, :, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_write_buffer_3_i4))

  if (present(offset)) then
    call accel_write_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_write_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_write_buffer_3_i4))
end subroutine X(accel_write_buffer_3_i4)


subroutine X(accel_read_buffer_0)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE,                           intent(out)   :: data
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(i8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif

  PUSH_SUB(X(accel_read_buffer_0))
  call profiling_in(prof_read, TOSTRING(X(CL_READ_BUFFER)))

  ! it does not make sense to read a buffer that the kernels cannot write
  ASSERT(this%flags /= ACCEL_MEM_READ_ONLY)

  fsize = size*R_SIZEOF
  offset_ = 0
  if (present(offset)) offset_ = offset*R_SIZEOF

  async_ = optional_default(async, .false.)

  if (fsize > 0) then
#ifdef HAVE_OPENCL
    call clEnqueueReadBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data, ierr)
    if (ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueReadBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_dtoh(this%mem, data, fsize, offset_)
#endif

    call profiling_count_transfers(size, data)
    if (.not. async_) call accel_finish()
  end if

  call profiling_out(prof_read)
  POP_SUB(X(accel_read_buffer_0))
end subroutine X(accel_read_buffer_0)


subroutine X(accel_read_buffer_1)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE, contiguous,               intent(out)   :: data(:)
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  if(size == 0) return

  PUSH_SUB(X(accel_read_buffer_1))

  ASSERT(ubound(data, dim=1) >= size)

  call X(accel_read_buffer_0)(this, size, data(1), offset, async)

  POP_SUB(X(accel_read_buffer_1))
end subroutine X(accel_read_buffer_1)


subroutine X(accel_read_buffer_2)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE, contiguous,               intent(out)   :: data(:, :)
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  if(size == 0) return

  PUSH_SUB(X(accel_read_buffer_2))

  ASSERT(ubound(data, dim=1) >= 1)
  ASSERT(ubound(data, dim=2) >= 1)

  call X(accel_read_buffer_0)(this, size, data(1, 1), offset, async)

  POP_SUB(X(accel_read_buffer_2))
end subroutine X(accel_read_buffer_2)


subroutine X(accel_read_buffer_3)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer(i8),                      intent(in)    :: size
  R_TYPE, contiguous,               intent(out)   :: data(:, :, :)
  integer(i8),            optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  if(size == 0) return

  PUSH_SUB(X(accel_read_buffer_3))

  ASSERT(ubound(data, dim=1) >= 1)
  ASSERT(ubound(data, dim=2) >= 1)
  ASSERT(ubound(data, dim=3) >= 1)

  call X(accel_read_buffer_0)(this, size, data(1, 1, 1), offset, async)

  POP_SUB(X(accel_read_buffer_3))
end subroutine X(accel_read_buffer_3)


subroutine X(accel_read_buffer_0_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(out)   :: data
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_read_buffer_0_i4))

  if (present(offset)) then
    call accel_read_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_read_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_read_buffer_0_i4))
end subroutine X(accel_read_buffer_0_i4)


subroutine X(accel_read_buffer_1_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer,                          intent(in)    :: size
  R_TYPE, contiguous,               intent(out)   :: data(:)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_read_buffer_1_i4))

  if (present(offset)) then
    call accel_read_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_read_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_read_buffer_1_i4))
end subroutine X(accel_read_buffer_1_i4)


subroutine X(accel_read_buffer_2_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer,                          intent(in)    :: size
  R_TYPE, contiguous,               intent(out)   :: data(:, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_read_buffer_2_i4))

  if (present(offset)) then
    call accel_read_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_read_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_read_buffer_2_i4))
end subroutine X(accel_read_buffer_2_i4)


subroutine X(accel_read_buffer_3_i4)(this, size, data, offset, async)
  type(accel_mem_t),                intent(in)    :: this
  integer,                          intent(in)    :: size
  R_TYPE, contiguous,               intent(out)   :: data(:, :, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  PUSH_SUB(X(accel_read_buffer_3_i4))

  if (present(offset)) then
    call accel_read_buffer(this, int(size, i8), data, int(offset, i8), async)
  else
    call accel_read_buffer(this, int(size, i8), data, async=async)
  end if

  POP_SUB(X(accel_read_buffer_3_i4))
end subroutine X(accel_read_buffer_3_i4)


subroutine X(accel_set_kernel_arg_data)(kernel, narg, data)
  type(accel_kernel_t), intent(inout) :: kernel
  integer,              intent(in)    :: narg
  R_TYPE,               intent(in)    :: data

#ifdef HAVE_OPENCL
  integer :: ierr
#endif

  PUSH_SUB(X(accel_set_kernel_arg_data))

  ! no push_sub, called too frequently
#ifdef HAVE_CUDA
  call cuda_kernel_set_arg_value(kernel%arguments, data, narg, types_get_size(R_TYPE_VAL))
#endif

#ifdef HAVE_OPENCL
  call clSetKernelArg(kernel%kernel, narg, data, ierr)
  if (ierr /= CL_SUCCESS) call opencl_print_error(ierr, "set_kernel_arg_data")
#endif

  POP_SUB(X(accel_set_kernel_arg_data))
end subroutine X(accel_set_kernel_arg_data)


subroutine X(accel_get_device_pointer_1)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer,   intent(inout) :: host_pointer(:)
  type(accel_mem_t), intent(in)    :: device_pointer
  integer,           intent(in)    :: dimensions(:)

#ifdef HAVE_CUDA
  type(c_ptr) :: tmp_pointer
#endif

  PUSH_SUB(X(accel_get_device_pointer_1))

#ifdef HAVE_CUDA
  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)
#else
  ASSERT(.false.)
#endif


  POP_SUB(X(accel_get_device_pointer_1))
end subroutine X(accel_get_device_pointer_1)


subroutine X(accel_get_device_pointer_2)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer,   intent(inout) :: host_pointer(:, :)
  type(accel_mem_t), intent(in)    :: device_pointer
  integer,           intent(in)    :: dimensions(:)

#ifdef HAVE_CUDA
  type(c_ptr) :: tmp_pointer
#endif

  PUSH_SUB(X(accel_get_device_pointer_2))

#ifdef HAVE_CUDA
  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)
#else
  ASSERT(.false.)
#endif

  POP_SUB(X(accel_get_device_pointer_2))
end subroutine X(accel_get_device_pointer_2)

subroutine X(accel_get_device_pointer_3)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer,   intent(inout) :: host_pointer(:, :, :)
  type(accel_mem_t), intent(in)    :: device_pointer
  integer,           intent(in)    :: dimensions(:)

#ifdef HAVE_CUDA
  type(c_ptr) :: tmp_pointer
#endif

  PUSH_SUB(X(accel_get_device_pointer_3))

#ifdef HAVE_CUDA
  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)
#else
  ASSERT(.false.)
#endif

  POP_SUB(X(accel_get_device_pointer_3))
end subroutine X(accel_get_device_pointer_3)

subroutine X(accel_get_device_pointer_1l)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer,   intent(inout) :: host_pointer(:)
  type(accel_mem_t), intent(in)    :: device_pointer
  integer(i8),       intent(in)    :: dimensions(:)

#ifdef HAVE_CUDA
  type(c_ptr) :: tmp_pointer
#endif

  PUSH_SUB(X(accel_get_device_pointer_1l))

#ifdef HAVE_CUDA
  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)
#else
  ASSERT(.false.)
#endif

  POP_SUB(X(accel_get_device_pointer_1l))
end subroutine X(accel_get_device_pointer_1l)


subroutine X(accel_get_device_pointer_2l)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer,   intent(inout) :: host_pointer(:, :)
  type(accel_mem_t), intent(in)    :: device_pointer
  integer(i8),       intent(in)    :: dimensions(:)

#ifdef HAVE_CUDA
  type(c_ptr) :: tmp_pointer
#endif

  PUSH_SUB(X(accel_get_device_pointer_2l))

#ifdef HAVE_CUDA
  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)
#else
  ASSERT(.false.)
#endif

  POP_SUB(X(accel_get_device_pointer_2l))
end subroutine X(accel_get_device_pointer_2l)

subroutine X(accel_get_device_pointer_3l)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer,   intent(inout) :: host_pointer(:, :, :)
  type(accel_mem_t), intent(in)    :: device_pointer
  integer(i8),       intent(in)    :: dimensions(:)

#ifdef HAVE_CUDA
  type(c_ptr) :: tmp_pointer
#endif

  PUSH_SUB(X(accel_get_device_pointer_3l))

#ifdef HAVE_CUDA
  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)
#else
  ASSERT(.false.)
#endif

  POP_SUB(X(accel_get_device_pointer_3l))
end subroutine X(accel_get_device_pointer_3l)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
