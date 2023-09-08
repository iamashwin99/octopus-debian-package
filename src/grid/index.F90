!! Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
!! Copyright (C) 2021 S. Ohlmann
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

module index_oct_m
  use debug_oct_m
  use global_oct_m
  use iihash_oct_m
  use io_oct_m
  use io_binary_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                 &
    index_t,                &
    index_from_coords,      &
    index_from_coords_vec,  &
    index_to_coords,        &
    index_point_to_spatial, &
    index_spatial_to_point, &
    index_dump,             &
    index_load,             &
    get_blocked_index

  type index_t
    ! Components are public by default
    integer              :: dim              !< the dimension
    integer              :: nr(2, MAX_DIM)   !< dimensions of the box where the points are contained
    integer              :: ll(MAX_DIM)      !< literally nr(2,:) - nr(1,:) + 1 - 2*enlarge(:)
    integer              :: enlarge(MAX_DIM) !< number of points to add for boundary conditions
    integer(i8)          :: checksum
    integer(i8), pointer  :: grid_to_spatial_global(:) !< map: global grid index -> spatial index
    integer(i8), pointer  :: spatial_to_grid_global(:) !< inverse map: spatial index -> global grid index
    integer              :: type               !< index type
    integer              :: bits               !< bits per dimension for Hilbert index
    integer              :: offset(MAX_DIM)    !< offset for getting the indices from the spatial index
    integer              :: stride(MAX_DIM+1)
    integer              :: window_grid_to_spatial !< handle to shared memory window for map
    integer              :: window_spatial_to_grid !< handle to shared memory window for inverse map
  end type index_t

  integer, parameter, public :: &
    IDX_CUBIC   = 1, &
    IDX_HILBERT = 2


  interface
    subroutine hilbert_index_to_point(dim, nbits, index, point)
      use kind_oct_m
      implicit none

      integer,     intent(in)       :: dim
      integer,     intent(in)       :: nbits
      integer(i8), intent(in)       :: index
      integer,     intent(out)      :: point(*) !< (1:3)
    end subroutine hilbert_index_to_point

    subroutine hilbert_point_to_index(dim, nbits, index, point)
      use kind_oct_m
      implicit none

      integer,     intent(in)       :: dim
      integer,     intent(in)       :: nbits
      integer(i8), intent(out)      :: index
      integer,     intent(in)       :: point(*) !< (1:3)
    end subroutine hilbert_point_to_index

    subroutine hilbert_index_to_point_int(dim, nbits, index, point)
      use kind_oct_m
      implicit none

      integer,     intent(in)       :: dim
      integer,     intent(in)       :: nbits
      integer(i4), intent(in)       :: index
      integer,     intent(out)      :: point(*) !< (1:3)
    end subroutine hilbert_index_to_point_int

    subroutine hilbert_point_to_index_int(dim, nbits, index, point)
      use kind_oct_m
      implicit none

      integer,     intent(in)       :: dim
      integer,     intent(in)       :: nbits
      integer(i4), intent(out)      :: index
      integer,     intent(in)       :: point(*) !< (1:3)
    end subroutine hilbert_point_to_index_int
  end interface

contains

  !> This function takes care of the boundary conditions for a given
  !! vector of integer coordinates it returns the true _global_ index
  !! of the point.
  integer(i8) function index_from_coords(idx, ix) result(index)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: ix(1:idx%dim)

    integer(i8) :: ispatial

    ! No PUSH SUB, called too often
    call index_point_to_spatial(idx, idx%dim, ispatial, ix)
    if (ispatial < 0 .or. ispatial > ubound(idx%spatial_to_grid_global, dim=1, kind=i8)) then
      index = 0
    else
      index = idx%spatial_to_grid_global(ispatial)
    end if
  end function index_from_coords

  ! -----------------------------------------------------------------

  subroutine index_from_coords_vec(idx, npoints, ix, index)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: npoints
    integer,       intent(in)    :: ix(1:idx%dim, 1:npoints)
    integer(i8),   intent(out)   :: index(1:npoints)

    integer :: ip
    integer(i8) :: ispatial

    ! No PUSH SUB, called too often
    do ip = 1, npoints
      call index_point_to_spatial(idx, idx%dim, ispatial, ix(1:idx%dim, ip))
      if (ispatial < 0 .or. ispatial > ubound(idx%spatial_to_grid_global, dim=1, kind=i8)) then
        index(ip) = 0
      else
        index(ip) = idx%spatial_to_grid_global(ispatial)
      end if
    end do

  end subroutine index_from_coords_vec


  !> Given a _global_ point index, this function returns the set of
  !! integer coordinates of the point.
  subroutine index_to_coords(idx, ip, ix)
    type(index_t), intent(in)    :: idx
    integer(i8),   intent(in)    :: ip
    integer,       intent(out)   :: ix(:)

    ! We set all ix to zero first (otherwise the non-existent dimensions would be
    ! undefined on exit).
    ix = 0
    call index_spatial_to_point(idx, idx%dim, idx%grid_to_spatial_global(ip), ix)
  end subroutine index_to_coords

  ! --------------------------------------------------------------
  subroutine index_dump(idx, np, dir, mpi_grp, namespace, ierr)
    type(index_t),    intent(in)  :: idx
    integer(i8),      intent(in)  :: np
    character(len=*), intent(in)  :: dir
    type(mpi_grp_t),  intent(in)  :: mpi_grp
    type(namespace_t),intent(in)  :: namespace
    integer,          intent(out) :: ierr

    integer :: err

    PUSH_SUB(index_dump)

    ierr = 0

    if (mpi_grp_is_root(mpi_grp)) then
      ! the index array is a global function and only root will write
      ASSERT(associated(idx%grid_to_spatial_global))
      call io_binary_write(trim(io_workpath(dir, namespace))//"/indices.obf", np, &
        idx%grid_to_spatial_global, err)
      if (err /= 0) then
        ierr = ierr + 1
        message(1) = "Unable to write index function to '"//trim(dir)//"/indices.obf'."
        call messages_warning(1, namespace=namespace)
      end if
    end if

    call mpi_grp%bcast(ierr, 1, MPI_INTEGER, 0)

    POP_SUB(index_dump)
  end subroutine index_dump


  ! --------------------------------------------------------------
  !> Load the index arrays from a file
  subroutine index_load(idx, np, dir, mpi_grp, namespace, ierr)
    type(index_t),     intent(inout) :: idx
    integer(i8),       intent(in)    :: np
    character(len=*),  intent(in)    :: dir
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(out)   :: ierr

    integer :: err
    integer(i8) :: ip, global_size
    logical :: exists

    PUSH_SUB(index_load)

    ierr = 0

    ASSERT(associated(idx%grid_to_spatial_global))

    if (mpi_grp_is_root(mpi_grp)) then
      ! check for existence of lxyz.obf, print error message if found
      inquire(file=trim(trim(io_workpath(dir, namespace))//"/lxyz.obf"), exist=exists)
      if (exists) then
        message(1) = "Found lxyz.obf file. This means you created the restart files with an old version of the code."
        message(2) = "Please generate the restart files again with the current version of the code"
        message(3) = "because the internal format has changed."
        call messages_fatal(3)
      end if
      ! the index array is a global function and only root will write
      call io_binary_read(trim(io_workpath(dir, namespace))//"/indices.obf", np, &
        idx%grid_to_spatial_global, err)
      if (err /= 0) then
        ierr = ierr + 1
        message(1) = "Unable to read index function from '"//trim(dir)//"/indices.obf'."
        call messages_warning(1, namespace=namespace)
      end if
    end if

    ! Broadcast the results and synchronize
    call mpi_grp%bcast(ierr, 1, MPI_INTEGER, 0)
    if (ierr == 0) then
      call mpi_grp%bcast(idx%grid_to_spatial_global(1), np, MPI_INTEGER8, 0)
    end if

    ! fill global hash map
    global_size = product(idx%nr(2, 1:MAX_DIM) - idx%nr(1, 1:MAX_DIM) + 1)
    SAFE_ALLOCATE(idx%spatial_to_grid_global(0:global_size))
    ! TODO: remove this restriction
    ASSERT(np <= huge(0_i4))
    do ip = 1, np
      idx%spatial_to_grid_global(idx%grid_to_spatial_global(ip)) = i8_to_i4(ip)
    end do

    POP_SUB(index_load)
  end subroutine index_load

  subroutine index_spatial_to_point(idx, dim, ispatial, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(i8),   intent(in)  :: ispatial
    integer,       intent(out) :: point(1:dim)

    ! no push_sub/pop_sub, called too often
    select case (idx%type)
    case (IDX_CUBIC)
      call index_cubic_to_point(idx, dim, ispatial, point)
    case (IDX_HILBERT)
      call index_hilbert_to_point(idx, dim, ispatial, point)
    case default
      message(1) = "Unknown index type."
      call messages_fatal(1)
    end select
  end subroutine index_spatial_to_point

  subroutine index_point_to_spatial(idx, dim, ispatial, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(i8),   intent(out) :: ispatial
    integer,       intent(in)  :: point(1:dim)

    ! no push_sub/pop_sub, called too often
    select case (idx%type)
    case (IDX_CUBIC)
      call index_point_to_cubic(idx, dim, ispatial, point)
    case (IDX_HILBERT)
      call index_point_to_hilbert(idx, dim, ispatial, point)
    case default
      message(1) = "Unknown index type."
      call messages_fatal(1)
    end select
  end subroutine index_point_to_spatial

  subroutine index_cubic_to_point(idx, dim, icubic, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(i8),   intent(in)  :: icubic
    integer,       intent(out) :: point(1:dim)

    integer :: ii
    integer(i8) :: tmp

    ! no push_sub/pop_sub, called too often
    tmp = icubic
    do ii = dim, 2, -1
      point(ii) = int(tmp/idx%stride(ii) - idx%offset(ii), i4)
      tmp = mod(tmp, idx%stride(ii))
    end do
    point(1) = int(tmp - idx%offset(1), i4)
  end subroutine index_cubic_to_point

  subroutine index_point_to_cubic(idx, dim, icubic, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(i8),   intent(out) :: icubic
    integer,       intent(in)  :: point(1:dim)

    integer :: ii

    ! no push_sub/pop_sub, called too often
    icubic = 0
    do ii = 1, dim
      if (point(ii) < idx%nr(1, ii) .or. point(ii) > idx%nr(2, ii)) then
        icubic = -1
        return
      end if
      icubic = icubic + int(point(ii)+idx%offset(ii), i8) * idx%stride(ii)
    end do
  end subroutine index_point_to_cubic

  subroutine index_hilbert_to_point(idx, dim, ihilbert, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(i8),   intent(in)  :: ihilbert
    integer,       intent(out) :: point(:)

    ! no push_sub/pop_sub, called too often
    call hilbert_index_to_point(dim, idx%bits, ihilbert, point)
    point(1:dim) = point(1:dim) - idx%offset(1:dim)
  end subroutine index_hilbert_to_point

  subroutine index_point_to_hilbert(idx, dim, ihilbert, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(i8),   intent(out) :: ihilbert
    integer,       intent(in)  :: point(1:dim)

    integer :: point_copy(1:dim)

    ! no push_sub/pop_sub, called too often
    point_copy(1:dim) = point(1:dim) + idx%offset(1:dim)
    call hilbert_point_to_index(dim, idx%bits, ihilbert, point_copy)
  end subroutine index_point_to_hilbert

  ! get index along a curve that follows small parallelepipeds
  ! this corresponds to blocked loops over n-dimensional space
  integer(i8) function get_blocked_index(dim, point, bsize, number_of_blocks, increase_with_dimensions)
    integer,           intent(in) :: dim
    integer,           intent(in) :: point(1:dim)
    integer,           intent(in) :: bsize(1:dim)
    integer,           intent(in) :: number_of_blocks(1:dim)
    logical, optional, intent(in) :: increase_with_dimensions
    integer(i8) :: idim, j_local, j_block, stride_local, stride_block
    integer :: start, stop, step

    ! no push_sub/pop_sub, called too often
    if (optional_default(increase_with_dimensions, .true.)) then
      ! fastest increasing index is in first dimension
      start = 1
      stop = dim
      step = 1
    else
      ! fastest increasing index is in last dimension
      start = dim
      stop = 1
      step = -1
    end if

    ! j_local: index in local block
    ! j_block: block index
    j_local = 0
    j_block = 0
    stride_local = 1
    stride_block = 1
    do idim = start, stop, step
      j_local = j_local + mod(point(idim), bsize(idim)) * stride_local
      stride_local = stride_local * bsize(idim)
      j_block = j_block + point(idim) / bsize(idim) * stride_block
      stride_block = stride_block * number_of_blocks(idim)
    end do
    ! total index along the curve
    get_blocked_index = j_local + j_block * product(int(bsize(1:dim), i8))
  end function get_blocked_index
end module index_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
