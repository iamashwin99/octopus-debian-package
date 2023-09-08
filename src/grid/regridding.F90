!! Copyright (C) 2022 F. Bonafé, S. Ohlmann
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


! Implementation details for regridding
!
! Restriction and prolongation are chosen such that they are adjoint operators.
!
! 1. Simple mapping
! - find all points of output mesh that are also on the local input mesh
! - get their partition number on the output mesh
! - use the partition_transfer class to communicate those points
! - use the corresponding mappings for the input and output buffers
!
! 2. Restriction: mapping from fine to coarse meshes.
!
! The strategy is as follows:
! - find all points of the coarse mesh that correspond to local points of
!   the fine mesh
! - find all points of the coarse mesh around those that have been found
!   that are connected by an order-1 cube stencil - the idea is that those
!   points have contributions in the restriction operation from points on
!   the fine mesh
! - create the stencil for the restriction and compute the weights; the
!   size is 2*(eta-1)+1 in each dimension, where eta is the ratio of the
!   grid spacings
! - save the indices of all of these points for communication
! - when doing the transfer:
!   - apply the restriction by looping over the locally matching coarse
!     points and applying the restriction stencil, adding up the
!     contribution from the points on the fine mesh
!   - communicate those values
!   - do a reduction of the received values; each rank can receive
!     contributions from different ranks, depending on the parallelization
!
! 3. Prolongation: mapping from coarse to fine meshes using linear interpolation
! in nd space
!
! The strategy is slightly different from the restriction operator:
! - first, get all local coarse points corresponding to local fine points
!   on the output mesh, thus loop over the output mesh here
! - get all points on the coarse input mesh reachable by a cube stencil
! - these are now the points needed for the output mesh
! - use the inverse partition transfer to communicate the points; this is
!   needed because the normal partition transfer can only communicate
!   points that are locally available, but we also potentially need some
!   from other cores
! - when doing the transfer, loop over the coarse points and apply the
!   transfer stencil (which corresponds to linear interpolation) to the
!   adjacent fine mesh points


module regridding_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use iihash_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use partition_oct_m
  use partition_transfer_oct_m
  use profiling_oct_m
  use space_oct_m
  use stencil_oct_m
  use stencil_cube_oct_m
  use utils_oct_m

  implicit none
  public :: regridding_t

  type :: regridding_t
    private
    class(mesh_t), pointer :: mesh_in, mesh_out
    type(partition_transfer_t) :: partition_transfer
    integer :: nsend, nrecv, dim
    integer, allocatable :: order_in(:), order_out(:)
    integer(i8), allocatable :: order_in_global(:), order_out_global(:)
    logical :: do_restriction, do_prolongation
    type(stencil_t) :: transfer_stencil
    FLOAT, allocatable :: weights(:)
    integer, allocatable :: eta(:)
  contains
    procedure :: dregridding_do_transfer_1, zregridding_do_transfer_1
    procedure :: dregridding_do_transfer_2, zregridding_do_transfer_2
    generic :: do_transfer => dregridding_do_transfer_1, zregridding_do_transfer_1
    generic :: do_transfer => dregridding_do_transfer_2, zregridding_do_transfer_2
    final :: regridding_finalize
  end type regridding_t

  interface regridding_t
    procedure regridding_init
  end interface regridding_t

contains

  ! ---------------------------------------------------------
  !> Generate a re-mapping of points from mesh_in to mesh_out
  function regridding_init(mesh_out, mesh_in, space_in, namespace) result(this)
    class(mesh_t), target, intent(in) :: mesh_out
    class(mesh_t), target, intent(in) :: mesh_in
    class(space_t),        intent(in) :: space_in
    type(namespace_t),     intent(in) :: namespace
    class(regridding_t), pointer :: this

    integer :: ip_in, ip_out, index(1:space_in%dim), idim, ii, is, size_array
    integer(i8), allocatable :: global_indices_in(:), order_in(:), order_out(:)
    integer(i8), allocatable :: global_indices_in_tmp(:)
    integer, allocatable :: partition_in(:)
    integer(i8) :: ipg_out, ipg_in
    type(profile_t), save :: prof
    FLOAT :: spacing_ratio(1:space_in%dim)
    integer :: index_stencil(1:space_in%dim)
    logical :: same_eta, on_coarse_grid, found
    type(stencil_t) :: cube_stencil
    type(lihash_t) :: global_indices

    PUSH_SUB(regridding_init)

    call profiling_in(prof,"REGRIDDING_INIT")

    SAFE_ALLOCATE(this)

    ! check some conditions which are not yet supported
    if (mesh_out%coord_system%local_basis .or. &
      mesh_in%coord_system%local_basis .or. &
      .not. mesh_out%coord_system%orthogonal .or. &
      .not. mesh_in%coord_system%orthogonal .or. &
      .not. (mesh_out%coord_system%dim == space_in%dim) .or. &
      space_in%is_periodic()) then
      message(1) = "Currently, regridding is only possible for orthogonal, non-periodic systems."
      call messages_fatal(1, namespace=namespace)
    end if
    this%dim = space_in%dim
    spacing_ratio(:) = mesh_out%spacing(1:space_in%dim)/mesh_in%spacing(1:space_in%dim)

    this%do_restriction = all(spacing_ratio > M_ONE)
    this%do_prolongation = all(spacing_ratio < M_ONE)
    ! invert spacing ratio for prolongations
    if (this%do_prolongation) spacing_ratio = M_ONE/spacing_ratio
    ! get the integer ratio of the spacings
    SAFE_ALLOCATE(this%eta(1:space_in%dim))
    do idim = 1, space_in%dim
      this%eta(idim) = nint(spacing_ratio(idim))
    end do
    if (any(abs((spacing_ratio - this%eta)/spacing_ratio) > M_EPSILON)) then
      message(1) = 'Only commensurate grids allowed for regridding.'
      call messages_fatal(1, namespace=namespace)
    end if
    same_eta = .true.
    do idim = 2, space_in%dim
      same_eta = same_eta .and. this%eta(idim) == this%eta(idim-1)
    end do
    if (.not. same_eta) then
      message(1) = 'Commensurate grids need to have same ratio in all dimensions for regridding.'
      call messages_fatal(1, namespace=namespace)
    end if

    this%mesh_in => mesh_in
    this%mesh_out => mesh_out

    ! collect all locally available points in mesh_in that are also in mesh_out
    call lihash_init(global_indices)
    if (.not. this%do_prolongation) then
      ! same spacing or restriction
      size_array = mesh_in%np
      SAFE_ALLOCATE(global_indices_in_tmp(size_array))

      ii = 0
      do ip_in = 1, mesh_in%np
        call mesh_local_index_to_coords(mesh_in, ip_in, index)
        if (this%do_restriction) then
          on_coarse_grid = .true.
          do idim = 1, space_in%dim
            on_coarse_grid = on_coarse_grid .and. mod(index(idim), this%eta(idim)) == 0
          end do
          if (.not. on_coarse_grid) cycle
          ! translate between fine and coarse grid
          index(:) = index(:) / this%eta(:)
        end if
        ipg_out = mesh_global_index_from_coords(mesh_out, index)
        if (ipg_out > 0 .and. ipg_out <= mesh_out%np_global) then
          ii = ii + 1
          ! store global indices of mesh_out here
          global_indices_in_tmp(ii) = ipg_out
          call lihash_insert(global_indices, ipg_out, ii)
        end if
      end do
      if (this%do_restriction) then
        call stencil_cube_get_lapl(cube_stencil, space_in%dim, order=1)

        ! now get all coarse points that are reachable by a cube stencil of order 1
        ! these are needed for the restriction
        do ip_out = 1, ii
          call mesh_global_index_to_coords(mesh_out, global_indices_in_tmp(ip_out), index)
          do is = 1, cube_stencil%size
            if (cube_stencil%center == is) cycle
            index_stencil(:) = index(:) + cube_stencil%points(1:space_in%dim, is)
            ipg_out = mesh_global_index_from_coords(mesh_out, index_stencil)
            if (ipg_out > 0 .and. ipg_out <= mesh_out%np_global) then
              ip_in = lihash_lookup(global_indices, ipg_out, found)
              if (found) cycle
              ii = ii + 1
              ! enlarge array if necessary
              if (ii >= size_array) then
                size_array = size_array * 2
                call make_array_larger(global_indices_in_tmp, size_array)
              end if
              global_indices_in_tmp(ii) = ipg_out
              call lihash_insert(global_indices, ipg_out, ii)
            end if
          end do
        end do
        call stencil_end(cube_stencil)
      end if
      call lihash_end(global_indices)
    else
      ! collect points for prolongation
      ! here the logic is the other way round: we need to get all points
      ! we need for the prolongation to the finer output mesh
      size_array = mesh_out%np
      SAFE_ALLOCATE(global_indices_in_tmp(size_array))
      ii = 0
      do ip_out = 1, mesh_out%np
        call mesh_local_index_to_coords(mesh_out, ip_out, index)
        ! translate between fine and coarse grid
        on_coarse_grid = .true.
        do idim = 1, space_in%dim
          on_coarse_grid = on_coarse_grid .and. mod(index(idim), this%eta(idim)) == 0
        end do
        if (.not. on_coarse_grid) cycle
        ! translate between fine and coarse grid
        index(:) = index(:) / this%eta(:)
        ipg_in = mesh_global_index_from_coords(mesh_in, index)
        if (ipg_in > 0 .and. ipg_in <= mesh_in%np_global) then
          ii = ii + 1
          ! store global indices of mesh_in here
          global_indices_in_tmp(ii) = ipg_in
          call lihash_insert(global_indices, ipg_in, ii)
        end if
      end do
      call stencil_cube_get_lapl(cube_stencil, space_in%dim, order=1)

      ! now get all coarse points that are reachable by a cube stencil of order 1 on mesh_in
      ! these are needed for the prolongation
      do ip_in = 1, ii
        call mesh_global_index_to_coords(mesh_in, global_indices_in_tmp(ip_in), index)
        do is = 1, cube_stencil%size
          if (cube_stencil%center == is) cycle
          index_stencil(:) = index(:) + cube_stencil%points(1:space_in%dim, is)
          ipg_in = mesh_global_index_from_coords(mesh_in, index_stencil)
          if (ipg_in > 0 .and. ipg_in <= mesh_in%np_global) then
            ip_out = lihash_lookup(global_indices, ipg_in, found)
            if (found) cycle
            ii = ii + 1
            ! enlarge array if necessary
            if (ii >= size_array) then
              size_array = size_array * 2
              call make_array_larger(global_indices_in_tmp, size_array)
            end if
            global_indices_in_tmp(ii) = ipg_in
            call lihash_insert(global_indices, ipg_in, ii)
          end if
        end do
      end do
      call stencil_end(cube_stencil)
      call lihash_end(global_indices)
    end if

    SAFE_ALLOCATE(global_indices_in(ii))
    SAFE_ALLOCATE(partition_in(ii))
    global_indices_in(1:ii) = global_indices_in_tmp(1:ii)
    SAFE_DEALLOCATE_A(global_indices_in_tmp)

    if (.not. this%do_prolongation) then
      if (mesh_out%parallel_in_domains) then
        ! determine where the points of mesh_out are stored
        call partition_get_partition_number(mesh_out%partition, ii, global_indices_in, partition_in)
      else
        partition_in = 1
      end if
    else
      if (mesh_in%parallel_in_domains) then
        ! determine where the points of mesh_in are stored
        call partition_get_partition_number(mesh_in%partition, ii, global_indices_in, partition_in)
      else
        partition_in = 1
      end if
    end if

    ! Init partition transfer
    ! need inverse direction for prolongations
    call partition_transfer_init(this%partition_transfer, ii, global_indices_in, &
      mesh_in%mpi_grp, mesh_out%mpi_grp, partition_in, &
      this%nsend, this%nrecv, order_in, order_out, inverse=this%do_prolongation)
    ! convert from global to local indices
    SAFE_ALLOCATE(this%order_in(1:this%nsend))
    SAFE_ALLOCATE(this%order_in_global(1:this%nsend))
    SAFE_ALLOCATE(this%order_out(1:this%nrecv))
    SAFE_ALLOCATE(this%order_out_global(1:this%nrecv))

    ! get the mapping for mesh_in in the order given by the global indices of mesh_out
    do ip_in = 1, this%nsend
      if (this%do_restriction) then
        ! in this case, we have an index on mesh_out
        this%order_in_global(ip_in) = order_in(ip_in)
        this%order_in(ip_in) = 0
      else if (this%do_prolongation) then
        ! in this case, we have an index on mesh_in
        this%order_in_global(ip_in) = order_in(ip_in)
        this%order_in(ip_in) = mesh_global2local(mesh_in, this%order_in_global(ip_in))
      else
        ! convert back from the global index of mesh_out to a local index of mesh_in
        call mesh_global_index_to_coords(mesh_out, order_in(ip_in), index)
        this%order_in_global(ip_in) = mesh_global_index_from_coords(mesh_in, index)
        this%order_in(ip_in) = mesh_global2local(mesh_in, this%order_in_global(ip_in))
      end if
      if (.not. this%do_restriction) then
        if (this%order_in(ip_in) == 0 .or. this%order_in(ip_in) > mesh_in%np) then
          write(message(1),'(a,i10,a,i10)') "Error in regridding part 1: mesh point ", &
            this%order_in(ip_in), " is not stored in partition ", mesh_in%pv%partno
          call messages_fatal(1, namespace=namespace)
        end if
      end if
    end do

    ! for mapping back to the global grid after the transfer, convert to local indices of mesh_out
    do ip_out = 1, this%nrecv
      if (.not. this%do_prolongation) then
        this%order_out_global(ip_out) = order_out(ip_out)
        this%order_out(ip_out) = mesh_global2local(mesh_out, order_out(ip_out))
      else
        ! store the global index of mesh_in
        this%order_out_global(ip_out) = order_out(ip_out)
        this%order_out(ip_out) = 0
      end if
      if (.not. this%do_prolongation) then
        if (this%order_out(ip_out) == 0 .or. this%order_out(ip_out) > mesh_out%np) then
          write(message(1),'(a,i10,a,i10)') "Error in regridding part 2: mesh point ", &
            this%order_out(ip_out), " is not stored in partition ", mesh_out%pv%partno
          call messages_fatal(1, namespace=namespace)
        end if
      end if
    end do

    ! create transfer stencil
    if (this%do_restriction .or. this%do_prolongation) then
      call stencil_cube_get_lapl(this%transfer_stencil, space_in%dim, order=this%eta(1)-1)
      SAFE_ALLOCATE(this%weights(1:this%transfer_stencil%size))
      this%weights(:) = M_ONE
      do ii = 1, this%transfer_stencil%size
        do idim = 1, space_in%dim
          ! the weights correspond to linear interpolation in d dimensions
          this%weights(ii) = this%weights(ii) * &
            (M_ONE - abs(this%transfer_stencil%points(idim, ii))/TOFLOAT(this%eta(idim)))
          ! for the restriction, we need to divide by eta for each dimension
          ! like this, restriction is adjoint to prolongation
          if (this%do_restriction) then
            this%weights(ii) = this%weights(ii) / TOFLOAT(this%eta(idim))
          end if
        end do
      end do
    end if

    SAFE_DEALLOCATE_A(partition_in)
    SAFE_DEALLOCATE_A(global_indices_in)
    SAFE_DEALLOCATE_A(order_in)
    SAFE_DEALLOCATE_A(order_out)

    call profiling_out(prof)
    POP_SUB(regridding_init)
  end function regridding_init

  subroutine regridding_finalize(this)
    type(regridding_t), intent(inout) :: this

    PUSH_SUB(regridding_finalize)

    call partition_transfer_end(this%partition_transfer)
    SAFE_DEALLOCATE_A(this%eta)
    SAFE_DEALLOCATE_A(this%order_in)
    SAFE_DEALLOCATE_A(this%order_in_global)
    SAFE_DEALLOCATE_A(this%order_out)
    SAFE_DEALLOCATE_A(this%order_out_global)
    SAFE_DEALLOCATE_A(this%weights)

    POP_SUB(regridding_finalize)
  end subroutine regridding_finalize

#include "real.F90"
#include "regridding_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "regridding_inc.F90"
#include "undef.F90"

end module regridding_oct_m
