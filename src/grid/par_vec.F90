!! Copyright (C) 2005-2006 Florian Lorenzen, Heiko Appel, J. Alberdi-Rodriguez
!! Copyright (C) 2021 Sebastian Ohlmann
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

 !> Some general things and nomenclature:
 !!
 !! - Points that are stored only on one process are
 !!   called local points.
 !! - Local points that are stored redundantly on
 !!   another process because of the partitioning are
 !!   called ghost points.
 !! - Boundary points are stored locally such that each
 !!   process has all points it needs for the finite differences
 !! - np is the total number of inner points.
 !!
 !! A globally defined vector v has two parts:
 !! - v(1:np) are the inner points
 !! - v(np+1:np_part) are the boundary points
 !! In the typical case of zero boundary conditions
 !! v(np+1:np_part) is 0.
 !! The two parts are split according to the partitions.
 !! The result of this split are local vectors vl on each process
 !! which consist of three parts:
 !! - vl(1:np_local)                                     local points.
 !! - vl(np_local+1:np_local+np_ghost)                   ghost points.
 !! - vl(np_local+np_ghost+1:np_local+np_ghost+np_bndry) boundary points.
 !!
 !! Usage example for par_vec routines.
 !! \verbatim
 !!
 !! ! Initialize parallelization with mesh and operator op
 !! ! initialized and given.
 !! ! mesh       = sys%gr%mesh
 !! ! stencil    = op%stencil
 !!
 !! FLOAT              :: uu(np_global), vv(np_global)
 !! FLOAT, allocatable :: ul(:), vl(:), wl(:)
 !! type(mesh_t)       :: mesh
 !!
 !! ! Fill uu, vv with sensible values.
 !! ! ...
 !!
 !! ! Allocate space for local vectors.
 !! allocate(ul(np_part))
 !! allocate(vl(np_part))
 !! allocate(wl(np_part))
 !!
 !! ! Distribute vectors.
 !! call vec_scatter(vp, uu, ul)
 !! call vec_scatter(vp, vv, vl)
 !!
 !! ! Compute some operator op: vl = op ul
 !! call X(vec_ghost_update)(vp, ul)
 !! call X(nl_operator_operate)(op, ul, vl)
 !! !! Gather result of op in one vector vv.
 !! call vec_gather(vp, vv, vl)
 !!
 !! ! Clean up.
 !! deallocate(ul, vl, wl)
 !! \endverbatim
module par_vec_oct_m
  use accel_oct_m
  use debug_oct_m
  use global_oct_m
  use iihash_oct_m
  use index_oct_m
  use io_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use partition_oct_m
  use profiling_oct_m
  use sort_oct_m
  use space_oct_m
  use stencil_oct_m
  use types_oct_m
  use utils_oct_m

  implicit none

  private

  public ::                &
    par_vec_t,             & !< parallel information
    par_vec_init,          &
    par_vec_end,           &
    par_vec_scatter,       &
    par_vec_gather,        &
    par_vec_allgather,     &
    par_vec_global2local,  &
    par_vec_local2global

  !> Parallel information
  type par_vec_t
    ! Components are public by default

    ! The content of these members is process-dependent.
    integer              :: rank                 !< Our rank in the communicator.
    !> Partition number of the
    !! current process
    integer              :: partno
    integer, allocatable :: ghost_sendpos(:)  !< The positions of the points for the ghost communication

    integer, allocatable :: ghost_rdispls(:)  !< Ghost points receive displacements
    integer, allocatable :: ghost_sdispls(:)  !< Ghost points send displacements
    integer, allocatable :: ghost_rcounts(:)  !< Number of ghost points to receive
    integer, allocatable :: ghost_scounts(:)  !< Number of ghost points to send
    integer              :: ghost_scount      !< Total number of ghost points to send
    integer, allocatable :: ghost_sendmap(:)  !< map for packing ghost points
    integer, allocatable :: ghost_recvmap(:)  !< map for unpacking ghost points
    type(accel_mem_t)    :: buff_sendmap      !< buffer for send map on GPUs
    type(accel_mem_t)    :: buff_recvmap      !< buffer for recv map on GPUs

    ! The following members are set independent of the processs.
    type(mpi_grp_t)         :: mpi_grp              !< MPI group to use
    integer                 :: npart                !< Number of partitions.
    integer(i8)             :: np_global            !< Number of points in mesh.

    integer, allocatable    :: np_local_vec(:)      !< How many points has partition r?
    !!                                                 Global vector; npart elements.
    integer                 :: np_local             !< How many points has running partition?
    !!                                                 Local value.
    integer(i8), allocatable, private :: xlocal_vec(:)  !< Points of partition r start at
    !!                                                 xlocal_vec(r) in local. Global start point
    !!                                                 of the local index.
    !!                                                 Global vector; npart elements.
    integer(i8)             :: xlocal               !< Starting index of running process in local(:) vector.
    !!                                                 Local value.

    integer(i8), allocatable, private    :: local(:)     !< Local points of running process
    !!                                                      Local vector; np_local elements
    integer, allocatable    :: recv_count(:)        !< Number of points to receive from all the other processes
    integer, allocatable    :: send_count(:)        !< Number of points to send to all the other processes
    !!                                                 in a MPI_Alltoallv.
    integer, allocatable    :: recv_disp(:)         !< Displacement of points to receive from all the other processes
    integer, allocatable    :: send_disp(:)         !< Displacement of points to send to all the other processes
    integer(i8), allocatable:: sendmap(:)           !< map for packing initial global points
    integer(i8), allocatable:: recvmap(:)           !< map for unpacking initial global points
    !!                                                 in a MPI_Alltoallv.
    integer                 :: np_bndry             !< Number of local boundary points.

    integer(i8), allocatable, private :: bndry(:)       !< local to global mapping of boundary points, np_bndry elements

    type(lihash_t), private :: global               !< global contains the global -> local mapping

    integer                 :: np_ghost             !< number of local ghost points
    integer(i8), allocatable, private :: ghost(:)       !< Global indices of ghost points, np_ghost elements
  end type par_vec_t

  interface par_vec_scatter
    module procedure dpar_vec_scatter
    module procedure zpar_vec_scatter
    module procedure ipar_vec_scatter
  end interface par_vec_scatter

  interface par_vec_gather
    module procedure dpar_vec_gather
    module procedure zpar_vec_gather
    module procedure ipar_vec_gather
  end interface par_vec_gather

  interface par_vec_allgather
    module procedure dpar_vec_allgather
    module procedure zpar_vec_allgather
    module procedure ipar_vec_allgather
  end interface par_vec_allgather

contains

  !> Initializes a par_vec_type object (parallel vector).
  !! It computes the local-to-global and global-to-local index tables
  !! and the ghost point exchange.
  !! \warning The naming scheme for the np_ variables is different
  !! from how it is in the rest of the code (for historical reasons
  !! and also because the vec_init has more a global than local point
  !! of view on the mesh): See the comments in the parameter list.
  subroutine par_vec_init(mpi_grp, np_global, idx, stencil, space, partition, pv, namespace)
    type(mpi_grp_t),       intent(in)    :: mpi_grp        !< MPI group to use.
    !> The next seven entries come from the mesh.
    integer(i8),           intent(in)    :: np_global      !< mesh%np_global
    type(index_t),         intent(in)    :: idx
    type(stencil_t),       intent(in)    :: stencil        !< The stencil for which to calculate ghost points.
    type(space_t),         intent(in)    :: space
    type(partition_t),     intent(in)    :: partition
    type(par_vec_t),       intent(inout) :: pv             !< Description of partition.
    type(namespace_t),     intent(in)    :: namespace

    ! Careful: MPI counts process ranks from 0 to numproc-1.
    ! Partition numbers from METIS range from 1 to numproc.
    ! For this reason, all ranks are incremented by one.
    integer                     :: npart            !< Number of partitions.
    integer(i8)                 :: gip, index
    integer                     :: ip, jp, jj, inode
    integer                     :: p1(MAX_DIM)      !< Points.
    type(lihash_t)              :: boundary
    type(lihash_t)              :: ghost
    integer(i8), allocatable    :: boundary_inv(:), ghost_inv(:)
    integer                     :: iunit            !< For debug output to files.
    character(len=6)            :: filenum
    logical                     :: found

    integer                     :: idir, ipart, itmp
    integer(i8)                 :: tmp
    integer, allocatable        :: points(:), part_ghost(:), part_ghost_tmp(:)
    integer(i8), allocatable    :: indices(:), ghost_tmp(:)

    PUSH_SUB(par_vec_init)

    ! Shortcuts.
    npart = mpi_grp%size

    ! Store partition number and rank for later reference.
    ! Having both variables is a bit redundant but makes the code readable.
    pv%rank = mpi_grp%rank
    pv%partno = pv%rank + 1

    call mpi_grp_copy(pv%mpi_grp, mpi_grp)
    pv%np_global = np_global
    pv%npart     = npart


    SAFE_ALLOCATE(pv%np_local_vec(1:npart))
    SAFE_ALLOCATE(pv%xlocal_vec(1:npart))
    SAFE_ALLOCATE(pv%ghost_rcounts(1:npart))
    SAFE_ALLOCATE(pv%ghost_scounts(1:npart))

    ! Count number of points for each process.
    ! Local points.
    call partition_get_np_local(partition, pv%np_local_vec)
    pv%np_local = pv%np_local_vec(pv%partno)


    ! Set up local-to-global index table for local points (xlocal_vec, local)
    pv%xlocal_vec(1) = 1
    ! Set the starting point of local and boundary points
    do inode = 2, npart
      pv%xlocal_vec(inode) = pv%xlocal_vec(inode - 1) + pv%np_local_vec(inode - 1)
    end do
    pv%xlocal = pv%xlocal_vec(pv%partno)

    SAFE_ALLOCATE(pv%local(pv%xlocal:pv%xlocal + pv%np_local - 1))
    ! Calculate the local vector in parallel
    call partition_get_local(partition, pv%local(pv%xlocal:), itmp)

    ! Create hash table.
    call lihash_init(pv%global)
    ! Insert local points.
    do ip = 1, pv%np_local
      call lihash_insert(pv%global, pv%local(pv%xlocal + ip - 1), ip)
    end do

    ! reorder points if needed
    call reorder_points()

    call lihash_init(boundary)
    call lihash_init(ghost)
    SAFE_ALLOCATE(boundary_inv(1:pv%np_local))
    SAFE_ALLOCATE(ghost_inv(1:pv%np_local))
    pv%np_ghost = 0
    pv%np_bndry = 0
    do gip = pv%xlocal, pv%xlocal + pv%np_local - 1
      ! Get coordinates of current point.
      call index_to_coords(idx, pv%local(gip), p1)
      ! For all points in stencil.
      do jj = 1, stencil%size
        ! Get point number of possible ghost point.
        index = index_from_coords(idx, p1(:) + stencil%points(:, jj))
        ASSERT(index /= 0)
        ! check if this point is a local point
        tmp = lihash_lookup(pv%global, index, found)
        if (found) cycle
        ! now check if the point is a potential ghost or boundary point
        if (index > np_global) then
          tmp = lihash_lookup(boundary, index, found)
          if (found) cycle
          pv%np_bndry = pv%np_bndry + 1
          call lihash_insert(boundary, index, pv%np_bndry)
          if (pv%np_bndry >= ubound(boundary_inv, 1)) call make_array_larger(boundary_inv, pv%np_bndry*2)
          boundary_inv(pv%np_bndry) = index
        else
          tmp = lihash_lookup(ghost, index, found)
          if (found) cycle
          pv%np_ghost = pv%np_ghost + 1
          call lihash_insert(ghost, index, pv%np_ghost)
          if (pv%np_ghost >= ubound(ghost_inv, 1)) call make_array_larger(ghost_inv, pv%np_ghost*2)
          ghost_inv(pv%np_ghost) = index
        end if
      end do
    end do

    SAFE_ALLOCATE(pv%bndry(1:pv%np_bndry))
    do ip = 1, pv%np_bndry
      pv%bndry(ip) = boundary_inv(ip)
    end do

    ! first get the temporary array of ghost points, will be later reorder by partition
    SAFE_ALLOCATE(ghost_tmp(1:pv%np_ghost))
    do ip = 1, pv%np_ghost
      ghost_tmp(ip) = ghost_inv(ip)
    end do
    call lihash_end(ghost)
    call lihash_end(boundary)
    SAFE_DEALLOCATE_A(ghost_inv)
    SAFE_DEALLOCATE_A(boundary_inv)

    SAFE_ALLOCATE(part_ghost_tmp(1:pv%np_ghost))
    call partition_get_partition_number(partition, pv%np_ghost, &
      ghost_tmp, part_ghost_tmp)

    ! determine parallel distribution (counts, displacements)
    pv%ghost_rcounts(:) = 0
    do ip = 1, pv%np_ghost
      ipart = part_ghost_tmp(ip)
      pv%ghost_rcounts(ipart) = pv%ghost_rcounts(ipart)+1
    end do
    ASSERT(sum(pv%ghost_rcounts) == pv%np_ghost)

    SAFE_ALLOCATE(pv%ghost_rdispls(1:pv%npart))
    pv%ghost_rdispls(1) = 0
    do ipart = 2, pv%npart
      pv%ghost_rdispls(ipart) = pv%ghost_rdispls(ipart - 1) + pv%ghost_rcounts(ipart - 1)
    end do

    ! reorder points by partition
    SAFE_ALLOCATE(pv%ghost(1:pv%np_ghost))
    SAFE_ALLOCATE(part_ghost(1:pv%np_ghost))
    SAFE_ALLOCATE(points(1:pv%npart))
    points = 0
    do ip = 1, pv%np_ghost
      ipart = part_ghost_tmp(ip)
      points(ipart) = points(ipart)+1
      ! jp is the new index, sorted according to partitions
      jp = pv%ghost_rdispls(ipart) + points(ipart)
      pv%ghost(jp) = ghost_tmp(ip)
      part_ghost(jp) = part_ghost_tmp(ip)
    end do
    SAFE_DEALLOCATE_A(points)
    SAFE_DEALLOCATE_A(ghost_tmp)
    SAFE_DEALLOCATE_A(part_ghost_tmp)

    call pv%mpi_grp%alltoall(pv%ghost_rcounts, 1, MPI_INTEGER, &
      pv%ghost_scounts, 1, MPI_INTEGER)

    SAFE_ALLOCATE(pv%ghost_sdispls(1:pv%npart))

    pv%ghost_sdispls(1) = 0
    do ipart = 2, pv%npart
      pv%ghost_sdispls(ipart) = pv%ghost_sdispls(ipart - 1) + pv%ghost_scounts(ipart - 1)
    end do
    pv%ghost_scount = sum(pv%ghost_scounts)

    SAFE_ALLOCATE(pv%ghost_recvmap(1:pv%np_ghost))
    SAFE_ALLOCATE(points(1:pv%npart))
    points = 0
    do ip = 1, pv%np_ghost
      ipart = part_ghost(ip)
      points(ipart) =  points(ipart) + 1
      pv%ghost_recvmap(ip) = pv%ghost_rdispls(ipart) + points(ipart)
    end do
    SAFE_DEALLOCATE_A(points)

    SAFE_ALLOCATE(indices(1:pv%np_ghost))
    do ip = 1, pv%np_ghost
      indices(pv%ghost_recvmap(ip)) = pv%ghost(ip)
    end do
    SAFE_ALLOCATE(pv%ghost_sendmap(1:pv%ghost_scount))
    SAFE_ALLOCATE(ghost_tmp(1:pv%ghost_scount))
    call pv%mpi_grp%alltoallv(indices, pv%ghost_rcounts, pv%ghost_rdispls, MPI_INTEGER8, &
      ghost_tmp, pv%ghost_scounts, pv%ghost_sdispls, MPI_INTEGER8)
    do ip = 1, pv%ghost_scount
      ! get local index
      jp = lihash_lookup(pv%global, ghost_tmp(ip), found)
      ASSERT(found)
      pv%ghost_sendmap(ip) = jp
    end do
    SAFE_DEALLOCATE_A(indices)
    SAFE_DEALLOCATE_A(ghost_tmp)

    if (accel_is_enabled()) then
      ! copy maps to GPU
      call accel_create_buffer(pv%buff_sendmap, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, pv%ghost_scount)
      call accel_write_buffer(pv%buff_sendmap, pv%ghost_scount, pv%ghost_sendmap)

      call accel_create_buffer(pv%buff_recvmap, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, pv%np_ghost)
      call accel_write_buffer(pv%buff_recvmap, pv%np_ghost, pv%ghost_recvmap)
    end if

    if (debug%info) then
      ! Write numbers and coordinates of each process` ghost points
      ! to a single file (like in mesh_partition_init) called
      ! debug/mesh_partition/ghost_points.###.
      call io_mkdir('debug/mesh_partition', namespace, parents=.true.)

      write(filenum, '(i6.6)') pv%partno
      iunit = io_open('debug/mesh_partition/ghost_points.'//filenum, namespace, action='write')
      do ip = 1, pv%np_ghost
        call index_to_coords(idx, pv%ghost(ip), p1)
        write(iunit, '(99i8)') pv%ghost(ip), (p1(idir), idir = 1, space%dim)
      end do

      call io_close(iunit)
    end if

    call init_MPI_Alltoall()

    ! Insert ghost points.
    do ip = 1, pv%np_ghost
      call lihash_insert(pv%global, pv%ghost(ip), ip + pv%np_local)
    end do
    ! Insert boundary points.
    do ip = 1, pv%np_bndry
      call lihash_insert(pv%global, pv%bndry(ip), ip + pv%np_local + pv%np_ghost)
    end do

    SAFE_DEALLOCATE_A(part_ghost)

    POP_SUB(par_vec_init)

  contains
    subroutine reorder_points()
      integer :: ip, nn, direction
      integer(i8) :: ipg
      integer :: bsize(space%dim), order, point(space%dim), number_of_blocks(space%dim)
      integer(i8), allocatable :: reorder_indices(:)
      integer(i8), allocatable :: reordered(:)
      integer, allocatable :: global_indices(:)
      type(block_t) :: blk
      logical :: increase_with_dimension
      integer, parameter :: &
        ORDER_BLOCKS     =  1, &
        ORDER_CUBE       =  3, &
        ORDER_GLOBAL     =  4

      PUSH_SUB(par_vec_init.reorder_points)

      !%Variable MeshLocalOrder
      !%Default blocks
      !%Type integer
      !%Section Execution::Optimization
      !%Description
      !% This variable controls how the grid points are mapped to a
      !% linear array. This influences the performance of the code.
      !%Option order_blocks 1
      !% The grid is mapped using small parallelepipedic grids. The size
      !% of the blocks is controlled by <tt>MeshBlockSize</tt>.
      !%Option order_cube 3
      !% The grid is mapped using a full cube, i.e. without blocking.
      !%Option order_global 4
      !% Use the ordering from the global mesh
      !%End
      call parse_variable(namespace, 'MeshLocalOrder', ORDER_BLOCKS, order)
      ! no reordering in 1D necessary
      if (space%dim == 1) then
        order = ORDER_GLOBAL
      end if

      !%Variable MeshLocalBlockDirection
      !%Type integer
      !%Section Execution::Optimization
      !%Description
      !% Determines the direction in which the dimensions are chosen to compute
      !% the blocked index for sorting the mesh points (see MeshLocalBlockSize).
      !% The default is increase_with_dimensions, corresponding to xyz ordering
      !% in 3D.
      !%Option increase_with_dimension 1
      !% The fastest changing index is in the first dimension, i.e., in 3D this
      !% corresponds to ordering in xyz directions.
      !%Option decrease_with_dimension 2
      !% The fastest changing index is in the last dimension, i.e., in 3D this
      !% corresponds to ordering in zyx directions.
      !%End
      call parse_variable(namespace, 'MeshLocalBlockDirection', 1, direction)
      increase_with_dimension = direction == 1
      if (direction /= 1 .and. direction /= 2) then
        call messages_input_error(namespace, 'MeshLocalBlockDirection')
      end if

      select case (order)
      case (ORDER_GLOBAL)
        ! nothing to do, points are ordered along the global distribution by default
      case (ORDER_BLOCKS, ORDER_CUBE)
        if (order == ORDER_CUBE) then
          bsize(1:space%dim) = idx%nr(2, 1:space%dim) - idx%nr(1, 1:space%dim) + 1
        else
          !%Variable MeshLocalBlockSize
          !%Type block
          !%Section Execution::Optimization
          !%Description
          !% To improve memory-access locality when calculating derivatives,
          !% <tt>Octopus</tt> arranges mesh points in blocks. This variable
          !% controls the size of this blocks in the different
          !% directions. The default is selected according to the value of
          !% the <tt>StatesBlockSize</tt> variable. (This variable only affects the
          !% performance of <tt>Octopus</tt> and not the results.)
          !%End
          if (conf%target_states_block_size < 16) then
            bsize(1) = 80 * 4 / abs(conf%target_states_block_size)
            if (space%dim > 1) bsize(2) = 4
            if (space%dim > 2) bsize(3:) = 10
          else
            bsize(1) = max(4 * 16 / abs(conf%target_states_block_size), 1)
            if (space%dim > 1) bsize(2) = 15
            if (space%dim > 2) bsize(3:) = 15
          end if

          if (parse_block(namespace, 'MeshLocalBlockSize', blk) == 0) then
            nn = parse_block_cols(blk, 0)
            if (nn /= space%dim) then
              message(1) = "Error: number of entries in MeshLocalBlockSize must match the  umber of dimensions."
              call messages_fatal(1, namespace=namespace)
            end if
            do idir = 1, nn
              call parse_block_integer(blk, 0, idir - 1, bsize(idir))
            end do
          end if
        end if

        number_of_blocks(1:space%dim) = idx%ll(1:space%dim)/bsize(1:space%dim) + 1

        ! compute new indices locally
        SAFE_ALLOCATE(reorder_indices(1:pv%np_local))
        SAFE_ALLOCATE(global_indices(1:pv%np_local))
        do ip = 1, pv%np_local
          ipg = pv%local(ip + pv%xlocal - 1)
          call index_spatial_to_point(idx, space%dim, idx%grid_to_spatial_global(ipg), point)
          point(1:space%dim) = point(1:space%dim) + idx%offset(1:space%dim) + idx%enlarge(1:space%dim)
          reorder_indices(ip) = get_blocked_index(space%dim, point, bsize, number_of_blocks, increase_with_dimension)
        end do
        ! sort the local array
        call sort(reorder_indices, global_indices)
        SAFE_DEALLOCATE_A(reorder_indices)

        ! reorder according to new order
        SAFE_ALLOCATE(reordered(1:pv%np_local))
        do ip = 1, pv%np_local
          reordered(ip) = pv%local(global_indices(ip) + pv%xlocal - 1)
        end do
        SAFE_DEALLOCATE_A(global_indices)

        ! assign the reordered points
        do ip = 1, pv%np_local
          pv%local(ip + pv%xlocal - 1) = reordered(ip)
        end do
        SAFE_DEALLOCATE_A(reordered)

        ! Recreate hash table.
        call lihash_end(pv%global)
        call lihash_init(pv%global)
        ! Insert local points.
        do ip = 1, pv%np_local
          call lihash_insert(pv%global, pv%local(pv%xlocal + ip - 1), ip)
        end do
      end select

      POP_SUB(par_vec_init.reorder_points)
    end subroutine reorder_points

    subroutine init_MPI_Alltoall()
      integer, allocatable :: part_local(:)

      PUSH_SUB(par_vec_init.init_MPI_Alltoall)

      SAFE_ALLOCATE(part_local(1:pv%np_local))
      SAFE_ALLOCATE(indices(1:pv%np_local))
      do ip = 1, pv%np_local
        indices(ip) = pv%xlocal + ip - 1
      end do
      call partition_get_partition_number(partition, pv%np_local, &
        indices, part_local)
      SAFE_DEALLOCATE_A(indices)

      SAFE_ALLOCATE(pv%send_count(1:npart))
      SAFE_ALLOCATE(pv%recv_count(1:npart))
      SAFE_ALLOCATE(pv%send_disp(1:npart))
      SAFE_ALLOCATE(pv%recv_disp(1:npart))

      pv%send_count = 0
      do ip = 1, pv%np_local
        ipart = part_local(ip)
        pv%send_count(ipart) = pv%send_count(ipart) + 1
      end do
      pv%send_disp(1) = 0
      do ipart = 2, npart
        pv%send_disp(ipart) = pv%send_disp(ipart - 1) + pv%send_count(ipart - 1)
      end do

      call pv%mpi_grp%alltoall(pv%send_count, 1, MPI_INTEGER, &
        pv%recv_count, 1, MPI_INTEGER)

      pv%recv_disp(1) = 0
      do ipart = 2, npart
        pv%recv_disp(ipart) = pv%recv_disp(ipart - 1) + pv%recv_count(ipart - 1)
      end do

      ! create maps
      SAFE_ALLOCATE(pv%sendmap(1:pv%np_local))
      SAFE_ALLOCATE(points(1:pv%npart))
      points = 0
      do ip = 1, pv%np_local
        ipart = part_local(ip)
        points(ipart) =  points(ipart) + 1
        pv%sendmap(ip) = pv%send_disp(ipart) + points(ipart)
      end do
      SAFE_DEALLOCATE_A(points)

      SAFE_ALLOCATE(indices(1:pv%np_local))
      do ip = 1, pv%np_local
        indices(pv%sendmap(ip)) = pv%xlocal + ip - 1
      end do
      SAFE_ALLOCATE(pv%recvmap(1:sum(pv%recv_count)))
      call pv%mpi_grp%alltoallv(indices, pv%send_count, pv%send_disp, MPI_INTEGER8, &
        pv%recvmap, pv%recv_count, pv%recv_disp, MPI_INTEGER8)
      do ip = 1, sum(pv%recv_count)
        ! get local index
        index = lihash_lookup(pv%global, pv%recvmap(ip), found)
        ASSERT(found)
        pv%recvmap(ip) = index
      end do
      SAFE_DEALLOCATE_A(indices)

      POP_SUB(par_vec_init.init_MPI_Alltoall)
    end subroutine init_MPI_Alltoall
  end subroutine par_vec_init

  ! ---------------------------------------------------------
  !> Deallocate memory used by pv.
  subroutine par_vec_end(pv)
    type(par_vec_t), intent(inout) :: pv

    PUSH_SUB(par_vec_end)

    SAFE_DEALLOCATE_A(pv%ghost_rdispls)
    SAFE_DEALLOCATE_A(pv%ghost_sdispls)
    SAFE_DEALLOCATE_A(pv%ghost_rcounts)
    SAFE_DEALLOCATE_A(pv%ghost_scounts)
    SAFE_DEALLOCATE_A(pv%ghost_sendpos)
    SAFE_DEALLOCATE_A(pv%send_disp)
    SAFE_DEALLOCATE_A(pv%recv_disp)
    SAFE_DEALLOCATE_A(pv%np_local_vec)
    SAFE_DEALLOCATE_A(pv%xlocal_vec)
    SAFE_DEALLOCATE_A(pv%local)
    SAFE_DEALLOCATE_A(pv%send_count)
    SAFE_DEALLOCATE_A(pv%recv_count)
    SAFE_DEALLOCATE_A(pv%bndry)
    SAFE_DEALLOCATE_A(pv%ghost)

    call lihash_end(pv%global)

    if (accel_is_enabled()) then
      call accel_release_buffer(pv%buff_recvmap)
      call accel_release_buffer(pv%buff_sendmap)
    end if

    POP_SUB(par_vec_end)
  end subroutine par_vec_end


  ! ---------------------------------------------------------
  !> Returns local number of global point ip on the local node
  !! If the result is zero, the point is not available on the local node
  integer function par_vec_global2local(pv, ipg) result(ip)
    type(par_vec_t), intent(in) :: pv
    integer(i8),     intent(in) :: ipg

    integer :: nn
    logical :: found

! no push_sub because called too frequently

    if (pv%npart > 1) then
      ip = 0
      nn = lihash_lookup(pv%global, ipg, found)
      if (found) ip = nn
    else
      ASSERT(ipg < huge(ip))
      ip = int(ipg, i4)
    end if
  end function par_vec_global2local

  ! ---------------------------------------------------------
  !> Returns global index of local point ip
  integer(i8) function par_vec_local2global(pv, ip) result(ipg)
    type(par_vec_t), intent(in) :: pv
    integer,         intent(in) :: ip

! no push_sub because called too frequently

    if (pv%npart == 1) then
      ipg = ip
    else
      if (ip <= pv%np_local) then
        ipg = pv%local(pv%xlocal + ip - 1)
      else if (ip <= pv%np_local + pv%np_ghost) then
        ipg = pv%ghost(ip - pv%np_local)
      else if (ip <= pv%np_local + pv%np_ghost + pv%np_bndry) then
        ipg = pv%bndry(ip - pv%np_local - pv%np_ghost)
      else
        ipg = 0_i8
      end if
    end if
  end function par_vec_local2global

  ! gather all local arrays into a global one on rank root
  ! this gives the global mapping of the index in the partition to the global index
  subroutine gather_local_vec(pv, root, local_vec)
    type(par_vec_t),          intent(in)    :: pv
    integer,                  intent(in)    :: root
    integer(i8), allocatable, intent(inout) :: local_vec(:)

    integer(i8), allocatable :: xlocal_tmp(:)

    PUSH_SUB(gather_local_vec)

    SAFE_ALLOCATE(xlocal_tmp(1:pv%npart))
    xlocal_tmp = pv%xlocal_vec - 1
    ! Gather all the local vectors in a unique big one
    call pv%mpi_grp%gatherv(pv%local(pv%xlocal:), pv%np_local, MPI_INTEGER8, &
      local_vec, pv%np_local_vec, xlocal_tmp, MPI_INTEGER8, root)
    SAFE_DEALLOCATE_A(xlocal_tmp)

    POP_SUB(gather_local_vec)
  end subroutine gather_local_vec

#include "undef.F90"
#include "complex.F90"
#include "par_vec_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "par_vec_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "par_vec_inc.F90"

end module par_vec_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
