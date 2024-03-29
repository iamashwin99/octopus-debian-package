!! Copyright (C) 2008 X. Andrade, 2020 S. Ohlmann
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

!> \ingroup Fortran_Module
module batch_oct_m
  use accel_oct_m
  use allocate_hardware_aware_oct_m
  use blas_oct_m
  use debug_oct_m
  use global_oct_m
  use hardware_oct_m
  use iso_c_binding
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private
  public ::                         &
    batch_t,                        &
    batch_init,                     &
    dbatch_init,                    &
    zbatch_init

  !> \ingroup Fortran_Class
  !! @brief Class defining batches of mesh functions
  type batch_t
    private
    integer,                        public :: nst
    integer,                        public :: dim
    integer                                :: np

    integer                                :: ndims
    integer,           allocatable         :: ist_idim_index(:, :)
    integer,           allocatable, public :: ist(:)

    logical                                :: is_allocated
    logical                                :: own_memory !< does the batch own the memory or is it foreign memory?
    !> We also need a linear array with the states in order to calculate derivatives, etc.
    integer,                        public :: nst_linear !< nst_linear = nst*std%dim

    integer                                :: status_of
    integer                                :: status_host
    type(type_t)                           :: type_of !< either TYPE_FLOAT or TYPE_COMPLEX
    integer                                :: device_buffer_count !< whether there is a copy in the opencl buffer
    integer                                :: host_buffer_count !< whether the batch was packed on the cpu
    logical                                :: special_memory
    logical                                :: needs_finish_unpack


    !> unpacked variables; linear variables are pointers with different shapes
    FLOAT, pointer, contiguous,     public :: dff(:, :, :)
    CMPLX, pointer, contiguous,     public :: zff(:, :, :)
    FLOAT, pointer, contiguous,     public :: dff_linear(:, :)
    CMPLX, pointer, contiguous,     public :: zff_linear(:, :)
    !> packed variables; only rank-2 arrays due to padding to powers of 2
    FLOAT, pointer, contiguous,     public :: dff_pack(:, :)
    CMPLX, pointer, contiguous,     public :: zff_pack(:, :)

    integer(i8),                    public :: pack_size(1:2)
    integer(i8),                    public :: pack_size_real(1:2) !< pack_size_real = pack_size; if batch type is complex, then pack_size_real(1) = 2*pack_size(1)

    type(accel_mem_t),              public :: ff_device

  contains
    procedure :: check_compatibility_with => batch_check_compatibility_with
    procedure :: clone_to => batch_clone_to
    procedure :: clone_to_array => batch_clone_to_array
    procedure :: copy_to => batch_copy_to
    procedure :: copy_data_to => batch_copy_data_to
    procedure :: do_pack => batch_do_pack
    procedure :: do_unpack => batch_do_unpack
    procedure :: finish_unpack => batch_finish_unpack
    procedure :: end => batch_end
    procedure :: inv_index => batch_inv_index
    procedure :: is_packed => batch_is_packed
    procedure :: ist_idim_to_linear => batch_ist_idim_to_linear
    procedure :: linear_to_idim => batch_linear_to_idim
    procedure :: linear_to_ist => batch_linear_to_ist
    procedure :: pack_total_size => batch_pack_total_size
    procedure :: remote_access_start => batch_remote_access_start
    procedure :: remote_access_stop => batch_remote_access_stop
    procedure :: status => batch_status
    procedure :: type => batch_type
    procedure :: type_as_int => batch_type_as_integer
    procedure, private :: dallocate_unpacked_host => dbatch_allocate_unpacked_host
    procedure, private :: zallocate_unpacked_host => zbatch_allocate_unpacked_host
    procedure, private :: allocate_unpacked_host => batch_allocate_unpacked_host
    procedure, private :: dallocate_packed_host => dbatch_allocate_packed_host
    procedure, private :: zallocate_packed_host => zbatch_allocate_packed_host
    procedure, private :: allocate_packed_host => batch_allocate_packed_host
    procedure, private :: allocate_packed_device => batch_allocate_packed_device
    procedure, private :: deallocate_unpacked_host => batch_deallocate_unpacked_host
    procedure, private :: deallocate_packed_host => batch_deallocate_packed_host
    procedure, private :: deallocate_packed_device => batch_deallocate_packed_device
  end type batch_t

  !--------------------------------------------------------------
  interface batch_init
    module procedure dbatch_init_with_memory_3
    module procedure zbatch_init_with_memory_3
    module procedure dbatch_init_with_memory_2
    module procedure zbatch_init_with_memory_2
    module procedure dbatch_init_with_memory_1
    module procedure zbatch_init_with_memory_1
  end interface batch_init

  integer, public, parameter :: &
    BATCH_NOT_PACKED     = 0,   &
    BATCH_PACKED         = 1,   &
    BATCH_DEVICE_PACKED  = 2

  integer, parameter :: CL_PACK_MAX_BUFFER_SIZE = 4 !< this value controls the size (in number of wave-functions)
  !!                                                   of the buffer used to copy states to the opencl device.

contains

  !--------------------------------------------------------------
  subroutine batch_end(this, copy)
    class(batch_t),          intent(inout) :: this
    logical,       optional, intent(in)    :: copy

    PUSH_SUB(batch_end)

    if (this%own_memory .and. this%is_packed()) then
      !deallocate directly to avoid unnecessary copies
      if (this%status() == BATCH_DEVICE_PACKED) then
        call this%deallocate_packed_device()
      end if
      if (this%status() == BATCH_PACKED .or. this%status_host == BATCH_PACKED) then
        call this%deallocate_packed_host()
      end if
      this%status_of = BATCH_NOT_PACKED
      this%status_host = BATCH_NOT_PACKED
      this%host_buffer_count = 0
      this%device_buffer_count = 0
    end if
    if (this%status() == BATCH_DEVICE_PACKED) call this%do_unpack(copy, force = .true.)
    if (this%status() == BATCH_PACKED) call this%do_unpack(copy, force = .true.)

    if (this%is_allocated) then
      call this%deallocate_unpacked_host()
    end if

    SAFE_DEALLOCATE_A(this%ist_idim_index)
    SAFE_DEALLOCATE_A(this%ist)

    POP_SUB(batch_end)
  end subroutine batch_end

  !--------------------------------------------------------------
  subroutine batch_deallocate_unpacked_host(this)
    class(batch_t),  intent(inout) :: this

    PUSH_SUB(batch_deallocate_unpacked_host)

    this%is_allocated = .false.

    if (this%special_memory) then
      if (associated(this%dff)) then
        call deallocate_hardware_aware(c_loc(this%dff(1,1,1)), int(this%np, i8)*this%dim*this%nst*8)
      end if
      if (associated(this%zff)) then
        call deallocate_hardware_aware(c_loc(this%zff(1,1,1)), int(this%np, i8)*this%dim*this%nst*16)
      end if
    else
      SAFE_DEALLOCATE_P(this%dff)
      SAFE_DEALLOCATE_P(this%zff)
    end if
    nullify(this%dff)
    nullify(this%dff_linear)
    nullify(this%zff)
    nullify(this%zff_linear)

    POP_SUB(batch_deallocate_unpacked_host)
  end subroutine batch_deallocate_unpacked_host

  !--------------------------------------------------------------
  subroutine batch_deallocate_packed_host(this)
    class(batch_t),  intent(inout) :: this

    PUSH_SUB(batch_deallocate_packed_host)

    if (this%special_memory) then
      if (associated(this%dff_pack)) then
        call deallocate_hardware_aware(c_loc(this%dff_pack(1,1)), int(this%pack_size(1), i8)*this%pack_size(2)*8)
      end if
      if (associated(this%zff_pack)) then
        call deallocate_hardware_aware(c_loc(this%zff_pack(1,1)), int(this%pack_size(1), i8)*this%pack_size(2)*16)
      end if
    else
      SAFE_DEALLOCATE_P(this%dff_pack)
      SAFE_DEALLOCATE_P(this%zff_pack)
    end if
    nullify(this%dff_pack)
    nullify(this%zff_pack)

    POP_SUB(batch_deallocate_packed_host)
  end subroutine batch_deallocate_packed_host

  !--------------------------------------------------------------
  subroutine batch_deallocate_packed_device(this)
    class(batch_t),  intent(inout) :: this

    PUSH_SUB(batch_deallocate_packed_device)

    call accel_release_buffer(this%ff_device)

    POP_SUB(batch_deallocate_packed_device)
  end subroutine batch_deallocate_packed_device

  !--------------------------------------------------------------
  subroutine batch_allocate_unpacked_host(this)
    class(batch_t),  intent(inout) :: this

    PUSH_SUB(batch_allocate_unpacked_host)

    if (this%type() == TYPE_FLOAT) then
      call this%dallocate_unpacked_host()
    else if (this%type() == TYPE_CMPLX) then
      call this%zallocate_unpacked_host()
    end if

    POP_SUB(batch_allocate_unpacked_host)
  end subroutine batch_allocate_unpacked_host

  !--------------------------------------------------------------
  subroutine batch_allocate_packed_host(this)
    class(batch_t),  intent(inout) :: this

    PUSH_SUB(batch_allocate_packed_host)

    if (this%type() == TYPE_FLOAT) then
      call this%dallocate_packed_host()
    else if (this%type() == TYPE_CMPLX) then
      call this%zallocate_packed_host()
    end if

    POP_SUB(batch_allocate_packed_host)
  end subroutine batch_allocate_packed_host

  !--------------------------------------------------------------
  subroutine batch_allocate_packed_device(this)
    class(batch_t),  intent(inout) :: this

    PUSH_SUB(batch_allocate_packed_device)

    call accel_create_buffer(this%ff_device, ACCEL_MEM_READ_WRITE, this%type(), &
      product(this%pack_size))

    POP_SUB(batch_allocate_packed_device)
  end subroutine batch_allocate_packed_device

  !--------------------------------------------------------------
  subroutine batch_init_empty (this, dim, nst, np)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: dim
    integer,       intent(in)    :: nst
    integer,       intent(in)    :: np

    PUSH_SUB(batch_init_empty)

    this%is_allocated = .false.
    this%own_memory = .false.
    this%special_memory = .false.
    this%needs_finish_unpack = .false.
    this%nst = nst
    this%dim = dim
    this%type_of = TYPE_NONE

    this%nst_linear = nst*dim

    this%np = np
    this%device_buffer_count = 0
    this%host_buffer_count = 0
    this%status_of = BATCH_NOT_PACKED
    this%status_host = BATCH_NOT_PACKED

    this%ndims = 2
    SAFE_ALLOCATE(this%ist_idim_index(1:this%nst_linear, 1:this%ndims))
    SAFE_ALLOCATE(this%ist(1:this%nst))

    nullify(this%dff, this%zff, this%dff_linear, this%zff_linear)
    nullify(this%dff_pack, this%zff_pack)

    POP_SUB(batch_init_empty)
  end subroutine batch_init_empty

  !--------------------------------------------------------------

  subroutine batch_clone_to(this, dest, pack, copy_data, new_np)
    class(batch_t),              intent(in)    :: this
    class(batch_t), allocatable, intent(out)   :: dest
    logical,        optional,    intent(in)    :: pack       !< If .false. the new batch will not be packed. Default: batch_is_packed(this)
    logical,        optional,    intent(in)    :: copy_data  !< If .true. the batch data will be copied to the destination batch. Default: .false.
    integer,        optional,    intent(in)    :: new_np

    PUSH_SUB(batch_clone_to)

    if (.not. allocated(dest)) then
      SAFE_ALLOCATE_TYPE(batch_t, dest)
    else
      message(1) = "Internal error: destination batch in batch_clone_to has been previously allocated."
      call messages_fatal(1)
    end if

    call this%copy_to(dest, pack, copy_data, new_np)

    POP_SUB(batch_clone_to)
  end subroutine batch_clone_to

  !--------------------------------------------------------------

  subroutine batch_clone_to_array(this, dest, n_batches, pack, copy_data)
    class(batch_t),              intent(in)    :: this
    class(batch_t), allocatable, intent(out)   :: dest(:)
    integer,                     intent(in)    :: n_batches
    logical,        optional,    intent(in)    :: pack       !< If .false. the new batch will not be packed. Default: batch_is_packed(this)
    logical,        optional,    intent(in)    :: copy_data  !< If .true. the batch data will be copied to the destination batch. Default: .false.

    integer :: ib

    PUSH_SUB(batch_clone_to_array)

    if (.not. allocated(dest)) then
      SAFE_ALLOCATE_TYPE_ARRAY(batch_t, dest, (1:n_batches))
    else
      message(1) = "Internal error: destination batch in batch_clone_to_array has been previously allocated."
      call messages_fatal(1)
    end if

    do ib = 1, n_batches
      call this%copy_to(dest(ib), pack, copy_data)
    end do

    POP_SUB(batch_clone_to_array)
  end subroutine batch_clone_to_array

  !--------------------------------------------------------------

  subroutine batch_copy_to(this, dest, pack, copy_data, new_np)
    class(batch_t),          intent(in)    :: this
    class(batch_t),          intent(out)   :: dest
    logical,       optional, intent(in)    :: pack       !< If .false. the new batch will not be packed. Default: batch_is_packed(this)
    logical,       optional, intent(in)    :: copy_data  !< If .true. the batch data will be copied to the destination batch. Default: .false.
    integer,       optional, intent(in)    :: new_np !< If present, this replaces this%np in the initialization

    logical :: host_packed, special
    integer :: np_

    PUSH_SUB(batch_copy_to)

    np_ = optional_default(new_np, this%np)

    host_packed = this%host_buffer_count > 0
    ! use special memory here only for batches not on the GPU to avoid allocating
    ! pinned memory for temporary batches because that leads to a severe performance
    ! decrease for GPU runs (up to 20x)
    special = this%special_memory .and. .not. this%device_buffer_count > 0
    if (this%type() == TYPE_FLOAT) then
      call dbatch_init(dest, this%dim, 1, this%nst, np_, packed=host_packed, special=special)
    else if (this%type() == TYPE_CMPLX) then
      call zbatch_init(dest, this%dim, 1, this%nst, np_, packed=host_packed, special=special)
    else
      message(1) = "Internal error: unknown batch type in batch_copy_to."
      call messages_fatal(1)
    end if

    if (this%status() /= dest%status() .and. optional_default(pack, this%is_packed())) call dest%do_pack(copy = .false.)

    dest%ist_idim_index(1:this%nst_linear, 1:this%ndims) = this%ist_idim_index(1:this%nst_linear, 1:this%ndims)
    dest%ist(1:this%nst) = this%ist(1:this%nst)

    if (optional_default(copy_data, .false.)) then
      ASSERT(np_ == this%np)
      call this%copy_data_to(min(this%np, np_), dest)
    end if

    POP_SUB(batch_copy_to)
  end subroutine batch_copy_to

  ! ----------------------------------------------------
  !> THREADSAFE
  type(type_t) pure function batch_type(this) result(btype)
    class(batch_t),      intent(in)    :: this

    btype = this%type_of

  end function batch_type

  ! ----------------------------------------------------
  !> For debuging purpose only
  integer pure function batch_type_as_integer(this) result(itype)
    class(batch_t),      intent(in)    :: this

    type(type_t) :: btype

    itype = 0
    btype = this%type()
    if (btype == TYPE_FLOAT) itype = 1
    if (btype == TYPE_CMPLX) itype = 2

  end function batch_type_as_integer

  ! ----------------------------------------------------
  !> THREADSAFE
  integer pure function batch_status(this) result(bstatus)
    class(batch_t),      intent(in)    :: this

    bstatus = this%status_of
  end function batch_status

  ! ----------------------------------------------------

  logical pure function batch_is_packed(this) result(in_buffer)
    class(batch_t),      intent(in)    :: this

    in_buffer = (this%device_buffer_count > 0) .or. (this%host_buffer_count > 0)
  end function batch_is_packed

  ! ----------------------------------------------------

  integer(i8) function batch_pack_total_size(this) result(size)
    class(batch_t),      intent(inout) :: this

    size = this%np
    if (accel_is_enabled()) size = accel_padded_size(size)
    size = size*pad_pow2(this%nst_linear)*types_get_size(this%type())

  end function batch_pack_total_size

  ! ----------------------------------------------------

  subroutine batch_do_pack(this, copy, async)
    class(batch_t),      intent(inout) :: this
    logical,   optional, intent(in)    :: copy
    logical,   optional, intent(in)    :: async

    logical               :: copy_
    logical               :: async_
    type(profile_t), save :: prof
    integer               :: source, target

    ! no push_sub, called too frequently

    call profiling_in(prof, "BATCH_DO_PACK")

    copy_ = optional_default(copy, .true.)

    async_ = optional_default(async, .false.)

    ! get source and target states for this batch
    source = this%status()
    select case (source)
    case (BATCH_NOT_PACKED, BATCH_PACKED)
      if (accel_is_enabled()) then
        target = BATCH_DEVICE_PACKED
      else
        target = BATCH_PACKED
      end if
    case (BATCH_DEVICE_PACKED)
      target = BATCH_DEVICE_PACKED
    end select

    ! only do something if target is different from source
    if (source /= target) then
      select case (target)
      case (BATCH_DEVICE_PACKED)
        call this%allocate_packed_device()
        this%status_of = BATCH_DEVICE_PACKED

        if (copy_) then
          select case (source)
          case (BATCH_NOT_PACKED)
            ! copy from unpacked host array to device
            call batch_write_unpacked_to_device(this)
          case (BATCH_PACKED)
            ! copy from packed host array to device
            call batch_write_packed_to_device(this, async_)
          end select
        end if
      case (BATCH_PACKED)
        call this%allocate_packed_host()
        this%status_of = BATCH_PACKED
        this%status_host = BATCH_PACKED

        if (copy_) then
          if (this%type() == TYPE_FLOAT) then
            call dbatch_pack_copy(this)
          else if (this%type() == TYPE_CMPLX) then
            call zbatch_pack_copy(this)
          end if
        end if
        if (this%own_memory) call this%deallocate_unpacked_host()
      end select
    end if

    select case (target)
    case (BATCH_DEVICE_PACKED)
      this%device_buffer_count = this%device_buffer_count + 1
    case (BATCH_PACKED)
      this%host_buffer_count = this%host_buffer_count + 1
    end select

    call profiling_out(prof)
  end subroutine batch_do_pack

  ! ----------------------------------------------------

  subroutine batch_do_unpack(this, copy, force, async)
    class(batch_t),     intent(inout) :: this
    logical, optional,  intent(in)    :: copy
    logical, optional,  intent(in)    :: force  !< if force = .true., unpack independently of the counter
    logical, optional,  intent(in)    :: async

    logical :: copy_, force_, async_
    type(profile_t), save :: prof
    integer               :: source, target

    PUSH_SUB(batch_do_unpack)

    call profiling_in(prof, "BATCH_DO_UNPACK")

    copy_ = optional_default(copy, .true.)

    force_ = optional_default(force, .false.)

    async_ = optional_default(async, .false.)

    ! get source and target states for this batch
    source = this%status()
    select case (source)
    case (BATCH_NOT_PACKED)
      target = source
    case (BATCH_PACKED)
      target = BATCH_NOT_PACKED
    case (BATCH_DEVICE_PACKED)
      target = this%status_host
    end select

    ! only do something if target is different from source
    if (source /= target) then
      select case (source)
      case (BATCH_PACKED)
        if (this%host_buffer_count == 1 .or. force_) then
          if (this%own_memory) call this%allocate_unpacked_host()
          ! unpack from packed_host to unpacked_host
          if (copy_ .or. this%own_memory) then
            if (this%type() == TYPE_FLOAT) then
              call dbatch_unpack_copy(this)
            else if (this%type() == TYPE_CMPLX) then
              call zbatch_unpack_copy(this)
            end if
          end if
          call this%deallocate_packed_host()
          this%status_host = target
          this%status_of = target
          this%host_buffer_count = 1
        end if
        this%host_buffer_count = this%host_buffer_count - 1
      case (BATCH_DEVICE_PACKED)
        if (this%device_buffer_count == 1 .or. force_) then
          if (copy_) then
            select case (target)
              ! unpack from packed_device to unpacked_host
            case (BATCH_NOT_PACKED)
              call batch_read_device_to_unpacked(this)
              ! unpack from packed_device to packed_host
            case (BATCH_PACKED)
              call batch_read_device_to_packed(this, async_)
            end select
          end if
          if (async_) then
            this%needs_finish_unpack = .true.
          else
            call this%deallocate_packed_device()
          end if
          this%status_of = target
          this%device_buffer_count = 1
        end if
        this%device_buffer_count = this%device_buffer_count - 1
      end select
    end if

    call profiling_out(prof)

    POP_SUB(batch_do_unpack)
  end subroutine batch_do_unpack

  ! ----------------------------------------------------
  subroutine batch_finish_unpack(this)
    class(batch_t),      intent(inout)  :: this

    PUSH_SUB(batch_finish_unpack)
    if (this%needs_finish_unpack) then
      call accel_finish()
      call this%deallocate_packed_device()
      this%needs_finish_unpack = .false.
    end if
    POP_SUB(batch_finish_unpack)
  end subroutine batch_finish_unpack

  ! ----------------------------------------------------

  subroutine batch_write_unpacked_to_device(this)
    class(batch_t),      intent(inout)  :: this

    integer :: ist, ist2
    integer(i8) :: unroll
    type(accel_mem_t) :: tmp
    type(profile_t), save :: prof, prof_pack
    type(accel_kernel_t), pointer :: kernel

    PUSH_SUB(batch_write_unpacked_to_device)

    call profiling_in(prof, "BATCH_WRT_UNPACK_ACCEL")
    if (this%nst_linear == 1) then
      ! we can copy directly
      if (this%type() == TYPE_FLOAT) then
        call accel_write_buffer(this%ff_device, ubound(this%dff_linear, dim=1), this%dff_linear(:, 1))
      else if (this%type() == TYPE_CMPLX) then
        call accel_write_buffer(this%ff_device, ubound(this%zff_linear, dim=1), this%zff_linear(:, 1))
      else
        ASSERT(.false.)
      end if

    else
      ! we copy to a temporary array and then we re-arrange data

      if (this%type() == TYPE_FLOAT) then
        kernel => dpack
      else
        kernel => zpack
      end if

      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack_size(1))

      call accel_create_buffer(tmp, ACCEL_MEM_READ_ONLY, this%type(), unroll*this%pack_size(2))

      do ist = 1, this%nst_linear, int(unroll, i4)

        ! copy a number 'unroll' of states to the buffer
        do ist2 = ist, min(ist + int(unroll, i4) - 1, this%nst_linear)

          if (this%type() == TYPE_FLOAT) then
            call accel_write_buffer(tmp, ubound(this%dff_linear, dim=1, kind=i8), this%dff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
          else
            call accel_write_buffer(tmp, ubound(this%zff_linear, dim=1, kind=i8), this%zff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
          end if
        end do

        ! now call an opencl kernel to rearrange the data
        call accel_set_kernel_arg(kernel, 0, int(this%pack_size(1), i4))
        call accel_set_kernel_arg(kernel, 1, int(this%pack_size(2), i4))
        call accel_set_kernel_arg(kernel, 2, ist - 1)
        call accel_set_kernel_arg(kernel, 3, tmp)
        call accel_set_kernel_arg(kernel, 4, this%ff_device)

        call profiling_in(prof_pack, "CL_PACK")
        call accel_kernel_run(kernel, (/this%pack_size(2), unroll/), (/accel_max_workgroup_size()/unroll, unroll/))

        if (this%type() == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack_size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack_size(2), M_ZI)
        end if

        call accel_finish()
        call profiling_out(prof_pack)

      end do

      call accel_release_buffer(tmp)

    end if

    call profiling_out(prof)
    POP_SUB(batch_write_unpacked_to_device)
  end subroutine batch_write_unpacked_to_device

  ! ------------------------------------------------------------------

  subroutine batch_read_device_to_unpacked(this)
    class(batch_t),      intent(inout) :: this

    integer :: ist, ist2
    integer(i8) :: unroll
    type(accel_mem_t) :: tmp
    type(accel_kernel_t), pointer :: kernel
    type(profile_t), save :: prof, prof_unpack

    PUSH_SUB(batch_read_device_to_unpacked)
    call profiling_in(prof, "BATCH_READ_UNPACKED_ACCEL")

    if (this%nst_linear == 1) then
      ! we can copy directly
      if (this%type() == TYPE_FLOAT) then
        call accel_read_buffer(this%ff_device, ubound(this%dff_linear, dim=1), this%dff_linear(:, 1))
      else
        call accel_read_buffer(this%ff_device, ubound(this%zff_linear, dim=1), this%zff_linear(:, 1))
      end if
    else

      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack_size(1))

      ! we use a kernel to move to a temporary array and then we read
      call accel_create_buffer(tmp, ACCEL_MEM_WRITE_ONLY, this%type(), unroll*this%pack_size(2))

      if (this%type() == TYPE_FLOAT) then
        kernel => dunpack
      else
        kernel => zunpack
      end if

      do ist = 1, this%nst_linear, int(unroll, i4)
        call accel_set_kernel_arg(kernel, 0, int(this%pack_size(1), i4))
        call accel_set_kernel_arg(kernel, 1, int(this%pack_size(2), i4))
        call accel_set_kernel_arg(kernel, 2, ist - 1)
        call accel_set_kernel_arg(kernel, 3, this%ff_device)
        call accel_set_kernel_arg(kernel, 4, tmp)

        call profiling_in(prof_unpack, "CL_UNPACK")
        call accel_kernel_run(kernel, (/unroll, this%pack_size(2)/), (/unroll, accel_max_workgroup_size()/unroll/))

        if (this%type() == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack_size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack_size(2), M_ZI)
        end if

        call accel_finish()
        call profiling_out(prof_unpack)

        ! copy a number 'unroll' of states from the buffer
        do ist2 = ist, min(ist + int(unroll, i4) - 1, this%nst_linear)

          if (this%type() == TYPE_FLOAT) then
            call accel_read_buffer(tmp, ubound(this%dff_linear, dim=1, kind=i8), this%dff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
          else
            call accel_read_buffer(tmp, ubound(this%zff_linear, dim=1, kind=i8), this%zff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
          end if
        end do

      end do

      call accel_release_buffer(tmp)
    end if

    call profiling_out(prof)
    POP_SUB(batch_read_device_to_unpacked)
  end subroutine batch_read_device_to_unpacked

  ! ------------------------------------------------------------------
  subroutine batch_write_packed_to_device(this, async)
    class(batch_t),      intent(inout)  :: this
    logical,   optional, intent(in)     :: async

    type(profile_t), save :: prof_pack

    PUSH_SUB(batch_write_packed_to_device)

    call profiling_in(prof_pack, "BATCH_WRITE_PACKED_ACCEL")
    if (this%type() == TYPE_FLOAT) then
      call accel_write_buffer(this%ff_device, product(this%pack_size), this%dff_pack, async=async)
    else
      call accel_write_buffer(this%ff_device, product(this%pack_size), this%zff_pack, async=async)
    end if
    call profiling_out(prof_pack)

    POP_SUB(batch_write_packed_to_device)
  end subroutine batch_write_packed_to_device

  ! ------------------------------------------------------------------
  subroutine batch_read_device_to_packed(this, async)
    class(batch_t),      intent(inout) :: this
    logical,   optional, intent(in)    :: async

    type(profile_t), save :: prof_unpack

    PUSH_SUB(batch_read_device_to_packed)

    call profiling_in(prof_unpack, "BATCH_READ_PACKED_ACCEL")
    if (this%type() == TYPE_FLOAT) then
      call accel_read_buffer(this%ff_device, product(this%pack_size), this%dff_pack, async=async)
    else
      call accel_read_buffer(this%ff_device, product(this%pack_size), this%zff_pack, async=async)
    end if
    call profiling_out(prof_unpack)

    POP_SUB(batch_read_device_to_packed)
  end subroutine batch_read_device_to_packed

! ------------------------------------------------------
  integer function batch_inv_index(this, cind) result(index)
    class(batch_t),     intent(in)    :: this
    integer,            intent(in)    :: cind(:)

    do index = 1, this%nst_linear
      if (all(cind(1:this%ndims) == this%ist_idim_index(index, 1:this%ndims))) exit
    end do

    ASSERT(index <= this%nst_linear)

  end function batch_inv_index

! ------------------------------------------------------

  integer pure function batch_ist_idim_to_linear(this, cind) result(index)
    class(batch_t),     intent(in)    :: this
    integer,            intent(in)    :: cind(:)

    if (ubound(cind, dim = 1) == 1) then
      index = cind(1)
    else
      index = (cind(1) - 1)*this%dim + cind(2)
    end if

  end function batch_ist_idim_to_linear

! ------------------------------------------------------

  integer pure function batch_linear_to_ist(this, linear_index) result(ist)
    class(batch_t),     intent(in)    :: this
    integer,            intent(in)    :: linear_index

    ist = this%ist_idim_index(linear_index, 1)

  end function batch_linear_to_ist

! ------------------------------------------------------

  integer pure function batch_linear_to_idim(this, linear_index) result(idim)
    class(batch_t),     intent(in)    :: this
    integer,            intent(in)    :: linear_index

    idim = this%ist_idim_index(linear_index, 2)

  end function batch_linear_to_idim

! ------------------------------------------------------

  subroutine batch_remote_access_start(this, mpi_grp, rma_win)
    class(batch_t),  intent(inout) :: this
    type(mpi_grp_t), intent(in)    :: mpi_grp
    integer,         intent(out)   :: rma_win

    PUSH_SUB(batch_remote_access_start)

    ASSERT(.not. accel_is_enabled())

    if (mpi_grp%size > 1) then
      call this%do_pack()

      if (this%type() == TYPE_CMPLX) then
#ifdef HAVE_MPI
        call MPI_Win_create(this%zff_pack(1, 1), int(product(this%pack_size)*types_get_size(this%type()), MPI_ADDRESS_KIND), &
          types_get_size(this%type()), MPI_INFO_NULL, mpi_grp%comm, rma_win, mpi_err)
#endif
      else if (this%type() == TYPE_FLOAT) then
#ifdef HAVE_MPI
        call MPI_Win_create(this%dff_pack(1, 1), int(product(this%pack_size)*types_get_size(this%type()), MPI_ADDRESS_KIND), &
          types_get_size(this%type()), MPI_INFO_NULL, mpi_grp%comm, rma_win, mpi_err)
#endif
      else
        message(1) = "Internal error: unknown batch type in batch_remote_access_start."
        call messages_fatal(1)
      end if

    else
      rma_win = -1
    end if

    POP_SUB(batch_remote_access_start)
  end subroutine batch_remote_access_start

! ------------------------------------------------------

  subroutine batch_remote_access_stop(this, rma_win)
    class(batch_t),     intent(inout) :: this
    integer,            intent(inout) :: rma_win

    PUSH_SUB(batch_remote_access_stop)

    if (rma_win /= -1) then
#ifdef HAVE_MPI
      call MPI_Win_free(rma_win, mpi_err)
#endif
      call this%do_unpack()
    end if

    POP_SUB(batch_remote_access_stop)
  end subroutine batch_remote_access_stop

! --------------------------------------------------------------

  subroutine batch_copy_data_to(this, np, dest)
    class(batch_t),    intent(in)    :: this
    integer,           intent(in)    :: np
    class(batch_t),    intent(inout) :: dest

    integer(i8) :: localsize, dim2, dim3
    type(profile_t), save :: prof
    integer :: ist, ip

    PUSH_SUB(batch_copy_data_to)
    call profiling_in(prof, "BATCH_COPY_DATA_TO")

    call this%check_compatibility_with(dest)

    select case (this%status())
    case (BATCH_DEVICE_PACKED)
      call accel_set_kernel_arg(kernel_copy, 0, np)
      call accel_set_kernel_arg(kernel_copy, 1, this%ff_device)
      call accel_set_kernel_arg(kernel_copy, 2, log2(int(this%pack_size_real(1), i4)))
      call accel_set_kernel_arg(kernel_copy, 3, dest%ff_device)
      call accel_set_kernel_arg(kernel_copy, 4, log2(int(dest%pack_size_real(1), i4)))

      localsize = accel_kernel_workgroup_size(kernel_copy)/dest%pack_size_real(1)

      dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(int(np, i8), localsize))

      call accel_kernel_run(kernel_copy, (/dest%pack_size_real(1), dim2, dim3/), (/dest%pack_size_real(1), localsize, 1_i8/))

      call accel_finish()

    case (BATCH_PACKED)
      if (np*this%pack_size(1) > huge(0_i4)) then
        ! BLAS cannot handle 8-byte integers, so we need a special version here
        do ip = 1, np
          if (dest%type() == TYPE_FLOAT) then
            call blas_copy(int(this%pack_size(1), i4), this%dff_pack(1, ip), 1, dest%dff_pack(1, ip), 1)
          else
            call blas_copy(int(this%pack_size(1), i4), this%zff_pack(1, ip), 1, dest%zff_pack(1, ip), 1)
          end if
        end do
      else
        if (dest%type() == TYPE_FLOAT) then
          !$omp parallel do private(ist)
          do ip = 1, np
            !$omp simd
            do ist = 1, int(this%pack_size(1), i4)
              dest%dff_pack(ist, ip) = this%dff_pack(ist, ip)
            end do
          end do
        else
          !$omp parallel do private(ist)
          do ip = 1, np
            !$omp simd
            do ist = 1, int(this%pack_size(1), i4)
              dest%zff_pack(ist, ip) = this%zff_pack(ist, ip)
            end do
          end do
        end if
      end if

    case (BATCH_NOT_PACKED)
      do ist = 1, dest%nst_linear
        if (dest%type() == TYPE_CMPLX) then
          call blas_copy(np, this%zff_linear(1, ist), 1, dest%zff_linear(1, ist), 1)
        else
          call blas_copy(np, this%dff_linear(1, ist), 1, dest%dff_linear(1, ist), 1)
        end if
      end do

    end select

    call profiling_out(prof)
    POP_SUB(batch_copy_data_to)
  end subroutine batch_copy_data_to

! --------------------------------------------------------------

  subroutine batch_check_compatibility_with(this, target, only_check_dim)
    class(batch_t),    intent(in) :: this
    class(batch_t),    intent(in) :: target
    logical, optional, intent(in) :: only_check_dim

    PUSH_SUB(batch_check_compatibility_with)

    ASSERT(this%type() == target%type())
    if (.not. optional_default(only_check_dim, .false.)) then
      ASSERT(this%nst_linear == target%nst_linear)
    end if
    ASSERT(this%status() == target%status())
    ASSERT(this%dim == target%dim)

    POP_SUB(batch_check_compatibility_with)

  end subroutine batch_check_compatibility_with

!--------------------------------------------------------------
  subroutine batch_build_indices(this, st_start, st_end)
    class(batch_t), intent(inout) :: this
    integer,        intent(in)    :: st_start
    integer,        intent(in)    :: st_end

    integer :: idim, ii, ist

    PUSH_SUB(batch_build_indices)

    do ist = st_start, st_end
      ! now we also populate the linear array
      do idim = 1, this%dim
        ii = this%dim*(ist - st_start) + idim
        this%ist_idim_index(ii, 1) = ist
        this%ist_idim_index(ii, 2) = idim
      end do
      this%ist(ist - st_start + 1) = ist
    end do

    ! compute packed sizes
    this%pack_size(1) = pad_pow2(this%nst_linear)
    this%pack_size(2) = this%np
    if (accel_is_enabled()) this%pack_size(2) = accel_padded_size(this%pack_size(2))

    this%pack_size_real = this%pack_size
    if (type_is_complex(this%type())) this%pack_size_real(1) = 2*this%pack_size_real(1)

    POP_SUB(batch_build_indices)
  end subroutine batch_build_indices


#include "real.F90"
#include "batch_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "batch_inc.F90"
#include "undef.F90"

end module batch_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
