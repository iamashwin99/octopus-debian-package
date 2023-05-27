!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module mesh_init_oct_m
  use affine_coordinates_oct_m
  use box_oct_m
  use checksum_interface_oct_m
  use coordinate_system_oct_m
  use debug_oct_m
  use global_oct_m
  use iihash_oct_m
  use index_oct_m
  use math_oct_m
  use merge_sorted_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use mesh_partition_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use partition_oct_m
  use profiling_oct_m
  use restart_oct_m
  use sort_oct_m
  use space_oct_m
  use stencil_oct_m
  use utils_oct_m

  implicit none

  private
  public ::                    &
    mesh_init_stage_1,         &
    mesh_init_stage_2,         &
    mesh_init_stage_3

  type(profile_t), save :: mesh_init_prof

  integer, parameter :: INNER_POINT = 1
  integer, parameter :: ENLARGEMENT_POINT = 2
  integer, parameter :: BOUNDARY = -1

contains

! ---------------------------------------------------------
  subroutine mesh_init_stage_1(mesh, namespace, space, box, coord_system, spacing, enlarge)
    type(mesh_t),                intent(inout) :: mesh
    type(namespace_t),           intent(in)    :: namespace
    type(space_t),               intent(in)    :: space
    class(box_t), target,        intent(in)    :: box
    class(coordinate_system_t), target, intent(in) :: coord_system
    FLOAT,                       intent(in)    :: spacing(1:MAX_DIM)
    integer,                     intent(in)    :: enlarge(MAX_DIM)

    integer :: idir, jj, delta
    FLOAT   :: x(MAX_DIM), chi(MAX_DIM), spacing_new(-1:1), box_bounds(2, space%dim)
    logical :: out

    PUSH_SUB(mesh_init_stage_1)
    call profiling_in(mesh_init_prof, "MESH_INIT")

    mesh%box => box

    mesh%spacing = spacing ! this number can change in the following
    mesh%use_curvilinear = coord_system%local_basis
    mesh%coord_system => coord_system

    mesh%idx%dim = space%dim
    mesh%idx%enlarge = enlarge

    ! get box bounds along the axes that generate the grid points
    select type (coord_system)
    class is (affine_coordinates_t)
      box_bounds = box%bounds(coord_system%basis)
    class default
      box_bounds = box%bounds()
    end select

    ! adjust nr
    mesh%idx%nr = 0
    do idir = 1, space%dim
      chi = M_ZERO
      ! the upper border
      jj = 0
      out = .false.
      do while(.not.out)
        jj = jj + 1
        chi(idir) = TOFLOAT(jj) * mesh%spacing(idir)
        if (mesh%use_curvilinear) then
          x(1:space%dim) = coord_system%to_cartesian(chi(1:space%dim))
          out = x(idir) > maxval(abs(box_bounds(:, idir))) + BOX_BOUNDARY_DELTA
        else
          ! do the same comparison here as in simul_box_contains_points
          out = chi(idir) > maxval(abs(box_bounds(:, idir))) + BOX_BOUNDARY_DELTA
        end if
      end do
      mesh%idx%nr(2, idir) = jj - 1
    end do

    ! we have a symmetric mesh (for now)
    mesh%idx%nr(1,:) = -mesh%idx%nr(2,:)

    ! we have to adjust a couple of things for the periodic directions
    do idir = 1, space%periodic_dim
      if (mesh%idx%nr(2, idir) == 0) then
        ! this happens if Spacing > box size
        mesh%idx%nr(2, idir) =  1
        mesh%idx%nr(1, idir) = -1
      end if

      ! We have to adjust the spacing to be commensurate with the box,
      ! for this we scan the possible values of the grid size around the
      ! one we selected. We choose the size that has the spacing closest
      ! to the requested one.
      do delta = -1, 1
        spacing_new(delta) = M_TWO*maxval(abs(box_bounds(:, idir))) / TOFLOAT(2 * mesh%idx%nr(2, idir) + 1 - delta)
        spacing_new(delta) = abs(spacing_new(delta) - spacing(idir))
      end do

      delta = minloc(spacing_new, dim = 1) - 2

      ASSERT(delta >= -1)
      ASSERT(delta <=  1)

      mesh%spacing(idir) = M_TWO*maxval(abs(box_bounds(:, idir))) / TOFLOAT(2 * mesh%idx%nr(2, idir) + 1 - delta)

      ! we need to adjust the grid by adding or removing one point
      if (delta == -1) then
        mesh%idx%nr(1, idir) = mesh%idx%nr(1, idir) - 1
      else if (delta == 1) then
        mesh%idx%nr(2, idir) = mesh%idx%nr(2, idir) - 1
      end if

    end do

    if ( any(abs(mesh%spacing(1:space%periodic_dim) - spacing(1:space%periodic_dim)) > CNST(1e-6)) ) then
      call messages_write('The spacing has been modified to make it commensurate with the periodicity of the system.')
      call messages_warning()
    end if

    do idir = space%periodic_dim + 1, space%dim
      if (mesh%idx%nr(2, idir) == 0) then
        write(message(1),'(a,i2)') 'Spacing > box size in direction ', idir
        call messages_fatal(1, namespace=namespace)
      end if
    end do

    mesh%idx%ll(1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) - mesh%idx%nr(1, 1:MAX_DIM) + 1
    ! compute strides for cubic indices
    mesh%idx%stride(:) = 1
    do idir = 2, space%dim
      mesh%idx%stride(idir) = mesh%idx%stride(idir-1) *  &
        (mesh%idx%ll(idir-1) + 2*mesh%idx%enlarge(idir-1))
    end do

    call profiling_out(mesh_init_prof)
    POP_SUB(mesh_init_stage_1)
  end subroutine mesh_init_stage_1

! ---------------------------------------------------------
!> This subroutine creates the global array of spatial indices
!! and the inverse mapping.
  subroutine mesh_init_stage_2(mesh, namespace, space, box, stencil)
    type(mesh_t),        intent(inout) :: mesh
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    class(box_t),        intent(in)    :: box
    type(stencil_t),     intent(in)    :: stencil

    integer :: is
    FLOAT   :: chi(1:space%dim)
    FLOAT :: pos(1:space%dim)
    integer :: point(1:MAX_DIM), point_stencil(1:MAX_DIM), grid_sizes(1:MAX_DIM)
    integer(i8) :: global_size
    integer(i4) :: local_size
    integer(i8) :: ispatial, ispatialb, istart, iend, spatial_size, ipg
    integer :: ip, ib, ib2, np, np_boundary, ii
    logical :: found
    type(lihash_t) :: spatial_to_boundary
    integer(i8), allocatable :: boundary_to_spatial(:), boundary_to_spatial_reordered(:)
    integer(i8), allocatable :: grid_to_spatial(:), grid_to_spatial_initial(:), grid_to_spatial_reordered(:)
    integer(i8), allocatable :: spatial_to_grid(:)
    integer, allocatable :: sizes(:)
    integer(i8), allocatable :: offsets(:)
    integer :: size_boundary
#ifdef HAVE_MPI
    integer(i8), pointer :: ptr(:)
    type(mpi_grp_t) :: internode_grp, intranode_grp
#endif

    PUSH_SUB(mesh_init_stage_2)
    call profiling_in(mesh_init_prof, "MESH_INIT")

    ! enlarge mesh for boundary points
    mesh%idx%nr(1, 1:MAX_DIM) = mesh%idx%nr(1, 1:MAX_DIM) - mesh%idx%enlarge(1:MAX_DIM)
    mesh%idx%nr(2, 1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) + mesh%idx%enlarge(1:MAX_DIM)

    !%Variable MeshIndexType
    !%Type integer
    !%Default idx_cubic
    !%Section Mesh
    !%Description
    !% Determine index type. Must be the same for restarting a calculation.
    !%Option idx_cubic 1
    !% Cubic indices are used to map the spatial information to the grid points.
    !%Option idx_hilbert 2
    !% A Hilbert space-filling curve is used to map the spatial information to
    !% the grid points.
    !%End
    call parse_variable(namespace, 'MeshIndexType', IDX_CUBIC, mesh%idx%type)

    grid_sizes(1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) - mesh%idx%nr(1, 1:MAX_DIM) + 1
    mesh%idx%offset(1:MAX_DIM) = grid_sizes(1:MAX_DIM)/2
    if (space%dim > 1 .and. any(grid_sizes > 2**(int(63/space%dim, i8)))) then
      write(message(1), '(A, I10, A, I2, A)') "Error: grid too large, more than ", 2**(int(63/space%dim, i8)), &
        " points in one direction for ", space%dim, " dimensions. This is not supported."
      call messages_fatal(1, namespace=namespace)
    end if
    global_size = product(int(grid_sizes, i8))
    ! compute the bits per dimension: grid_sizes(i) <= 2**bits
    mesh%idx%bits = maxval(ceiling(log(TOFLOAT(grid_sizes))/log(2.)))

    if (mesh%idx%type == IDX_CUBIC) then
      spatial_size = global_size
    else if (mesh%idx%type == IDX_HILBERT) then
      spatial_size = 2**(space%dim*mesh%idx%bits)
    end if

    ! use block data decomposition of spatial indices
    istart = spatial_size * mpi_world%rank/mpi_world%size
    iend = spatial_size * (mpi_world%rank+1)/mpi_world%size - 1
    ASSERT(iend - istart + 1 < huge(0_i4))
    local_size = int(iend - istart + 1, i4)

    SAFE_ALLOCATE(grid_to_spatial_initial(1:local_size))

    ! get inner grid indices
    ip = 1
    do ispatial = istart, iend
      call index_spatial_to_point(mesh%idx, space%dim, ispatial, point)
      ! first check if point is outside bounding box
      if (any(point(1:space%dim) < mesh%idx%nr(1, 1:space%dim) + mesh%idx%enlarge(1:space%dim))) cycle
      if (any(point(1:space%dim) > mesh%idx%nr(2, 1:space%dim) - mesh%idx%enlarge(1:space%dim))) cycle
      ! then check if point is inside simulation box
      chi(1:space%dim) = TOFLOAT(point(1:space%dim)) * mesh%spacing(1:space%dim)
      pos(1:space%dim) = mesh%coord_system%to_cartesian(chi(1:space%dim))
      if (.not. box%contains_point(pos)) cycle
      grid_to_spatial_initial(ip) = ispatial
      ASSERT(ip + 1 < huge(ip))
      ip = ip + 1
    end do
    np = ip - 1

    call rebalance_array(grid_to_spatial_initial(1:np), grid_to_spatial, sizes)
    np = sizes(mpi_world%rank)

    SAFE_DEALLOCATE_A(grid_to_spatial_initial)

    SAFE_ALLOCATE(spatial_to_grid(grid_to_spatial(1):grid_to_spatial(np)))
    SAFE_DEALLOCATE_A(sizes)

    !$omp parallel do
    do ispatial = grid_to_spatial(1), grid_to_spatial(np)
      spatial_to_grid(ispatial) = -1
    end do
    !$omp parallel do
    do ip = 1, np
      spatial_to_grid(grid_to_spatial(ip)) = ip
    end do

    ! get local boundary indices
    call lihash_init(spatial_to_boundary)
    size_boundary = np
    SAFE_ALLOCATE(boundary_to_spatial(1:size_boundary))
    ib = 1
    do ip = 1, np
      call index_spatial_to_point(mesh%idx, space%dim, grid_to_spatial(ip), point)
      do is = 1, stencil%size
        if (stencil%center == is) cycle
        point_stencil(1:space%dim) = point(1:space%dim) + stencil%points(1:space%dim, is)
        ! check if point is in inner part
        call index_point_to_spatial(mesh%idx, space%dim, ispatialb, point_stencil)
        if (ispatialb >= lbound(spatial_to_grid, dim=1, kind=i8) .and. &
          ispatialb <= ubound(spatial_to_grid, dim=1, kind=i8)) then
          if (spatial_to_grid(ispatialb) > 0) cycle
        end if
        ! then check if point is inside simulation box
        chi(1:space%dim) = TOFLOAT(point_stencil(1:space%dim)) * mesh%spacing(1:space%dim)
        pos(1:space%dim) = mesh%coord_system%to_cartesian(chi(1:space%dim))
        if (box%contains_point(pos)) cycle
        ! it has to be a boundary point now
        ! check if already counted
        ib2 = lihash_lookup(spatial_to_boundary, ispatialb, found)
        if (found) cycle
        boundary_to_spatial(ib) = ispatialb
        call lihash_insert(spatial_to_boundary, ispatialb, ib)
        ib = ib + 1
        ! enlarge array
        if (ib >= size_boundary) then
          size_boundary = size_boundary * 2
          call make_array_larger(boundary_to_spatial, size_boundary)
        end if
      end do
    end do
    np_boundary = ib - 1
    call lihash_end(spatial_to_boundary)
    SAFE_DEALLOCATE_A(spatial_to_grid)

    ! reorder inner points
    call reorder_points(namespace, space, mesh%idx, grid_to_spatial, grid_to_spatial_reordered)
    SAFE_DEALLOCATE_A(grid_to_spatial)

    call rebalance_array(grid_to_spatial_reordered, grid_to_spatial, sizes)
    np = sizes(mpi_world%rank)
    mesh%np_global = sizes(0)
    do ii = 1, mpi_world%size - 1
      mesh%np_global = mesh%np_global + sizes(ii)
    end do
    SAFE_DEALLOCATE_A(sizes)
    SAFE_DEALLOCATE_A(grid_to_spatial_reordered)

    ! reorder boundary points
    call make_array_larger(boundary_to_spatial, np_boundary)
    call reorder_points(namespace, space, mesh%idx, boundary_to_spatial, boundary_to_spatial_reordered)
    SAFE_DEALLOCATE_A(boundary_to_spatial)

    call rebalance_array(boundary_to_spatial_reordered, boundary_to_spatial, sizes)
    SAFE_DEALLOCATE_A(boundary_to_spatial_reordered)

    ! global grid size
    np_boundary = sizes(mpi_world%rank)
    mesh%np_part_global = mesh%np_global + sizes(0)
    do ii = 1, mpi_world%size - 1
      mesh%np_part_global = mesh%np_part_global + sizes(ii)
    end do
    SAFE_DEALLOCATE_A(sizes)


    ! get global indices
#ifdef HAVE_MPI
    ! create shared memory window and fill it only on root
    call create_intranode_communicator(mpi_world, intranode_grp, internode_grp)
    call lmpi_create_shared_memory_window(mesh%np_part_global, intranode_grp, &
      mesh%idx%window_grid_to_spatial, mesh%idx%grid_to_spatial_global)
#else
    SAFE_ALLOCATE(mesh%idx%grid_to_spatial_global(1:mesh%np_part_global))
#endif
    ! inner grid
    call get_sizes_offsets(np, sizes, offsets)
    call mpi_world%gatherv(grid_to_spatial, np, MPI_INTEGER8, &
      mesh%idx%grid_to_spatial_global, sizes, offsets, MPI_INTEGER8, 0)

    ! boundary indices
    call get_sizes_offsets(np_boundary, sizes, offsets)
    call mpi_world%gatherv(boundary_to_spatial, np_boundary, MPI_INTEGER8, &
      mesh%idx%grid_to_spatial_global(mesh%np_global+1:), sizes, offsets, MPI_INTEGER8, 0)

    ! fill global hash map
#ifdef HAVE_MPI
    ! create shared memory window and fill it only on root
    call lmpi_create_shared_memory_window(spatial_size, intranode_grp, &
      mesh%idx%window_spatial_to_grid, ptr)
    mesh%idx%spatial_to_grid_global(0:spatial_size-1) => ptr(1:spatial_size)
#else
    SAFE_ALLOCATE(mesh%idx%spatial_to_grid_global(0:spatial_size-1))
#endif
    if (mpi_grp_is_root(mpi_world)) then
      ! fill only on root, then broadcast
      !$omp parallel do
      do ispatial = 0, spatial_size-1
        mesh%idx%spatial_to_grid_global(ispatial) = 0
      end do
      !$omp parallel do
      do ipg = 1, mesh%np_part_global
        mesh%idx%spatial_to_grid_global(mesh%idx%grid_to_spatial_global(ipg)) = ipg
      end do
    end if

#ifdef HAVE_MPI
    ! now broadcast the global arrays to local rank 0 on each node
    if (intranode_grp%rank == 0) then
      call internode_grp%bcast(mesh%idx%grid_to_spatial_global(1), mesh%np_part_global, MPI_INTEGER8, 0)
      call internode_grp%bcast(mesh%idx%spatial_to_grid_global(0), spatial_size, MPI_INTEGER8, 0)
    end if
    call lmpi_sync_shared_memory_window(mesh%idx%window_grid_to_spatial, intranode_grp)
    call lmpi_sync_shared_memory_window(mesh%idx%window_spatial_to_grid, intranode_grp)
#endif

    SAFE_DEALLOCATE_A(offsets)
    SAFE_DEALLOCATE_A(sizes)

    SAFE_DEALLOCATE_A(boundary_to_spatial)
    SAFE_DEALLOCATE_A(grid_to_spatial)

    call profiling_out(mesh_init_prof)
    POP_SUB(mesh_init_stage_2)
  end subroutine mesh_init_stage_2

! ---------------------------------------------------------
!> When running parallel in domains, stencil and np_stencil
!! are needed to compute the ghost points.
!! mpi_grp is the communicator group that will be used for
!! this mesh.
! ---------------------------------------------------------
  subroutine mesh_init_stage_3(mesh, namespace, space, stencil, mc, parent)
    type(mesh_t),              intent(inout) :: mesh
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    type(stencil_t),           intent(in)    :: stencil
    type(multicomm_t),         intent(in)    :: mc
    type(mesh_t),    optional, intent(in)    :: parent

    integer :: ip

    PUSH_SUB(mesh_init_stage_3)
    call profiling_in(mesh_init_prof, "MESH_INIT")

    call mpi_grp_init(mesh%mpi_grp, mc%group_comm(P_STRATEGY_DOMAINS))

    ! check if we are running in parallel in domains
    mesh%parallel_in_domains = (mesh%mpi_grp%size > 1)

    call checksum_calculate(1, mesh%np_part_global, mesh%idx%grid_to_spatial_global(1), &
      mesh%idx%checksum)

    if (mesh%parallel_in_domains) then
      call do_partition()
    else
      ! When running serially those two are the same.
      ASSERT(mesh%np_part_global < huge(mesh%np_part))
      mesh%np      = i8_to_i4(mesh%np_global)
      mesh%np_part = i8_to_i4(mesh%np_part_global)

      ! These must be initialized for par_vec_gather, par_vec_scatter to work
      ! as copy operations when running without domain parallelization.
      mesh%pv%np_global = mesh%np_global
      mesh%pv%np_ghost = 0
      mesh%pv%np_bndry = mesh%np_part - mesh%np
      mesh%pv%npart = 1
      mesh%pv%xlocal = 1
    end if

    ! Compute mesh%x
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part, 1:space%dim))
    mesh%x(:, :) = M_ZERO
    do ip = 1, mesh%np_part
      mesh%x(ip, 1:space%dim) = mesh_x_global(mesh, mesh_local2global(mesh, ip))
    end do

    call mesh_cube_map_init(mesh%cube_map, mesh%idx, mesh%np_global)

    call mesh_get_vol_pp()

    call profiling_out(mesh_init_prof)
    POP_SUB(mesh_init_stage_3)

  contains
    ! ---------------------------------------------------------
    subroutine do_partition()
#ifdef HAVE_MPI
      integer :: jj, ipart, jpart
      integer(i8) :: ipg
      integer, allocatable :: gindex(:), gedges(:)
      logical, allocatable :: nb(:, :)
      integer              :: idx(1:MAX_DIM), jx(1:MAX_DIM)
      integer              :: graph_comm, iedge, nnb
      logical              :: use_topo, reorder, partition_print
      integer              :: ierr

      logical :: has_virtual_partition = .false.
      integer :: vsize !< 'virtual' partition size
      type(restart_t) :: restart_load, restart_dump
      integer, allocatable :: part_vec(:)

      PUSH_SUB(mesh_init_stage_3.do_partition)

      !Try to load the partition from the restart files
      call restart_init(restart_load, namespace, RESTART_PARTITION, RESTART_TYPE_LOAD, mc, ierr, mesh=mesh, exact=.true.)
      if (ierr == 0) call mesh_partition_load(restart_load, mesh, ierr)
      call restart_end(restart_load)

      if (ierr /= 0) then

        !%Variable MeshPartitionVirtualSize
        !%Type integer
        !%Default mesh mpi_grp size
        !%Section Execution::Parallelization
        !%Description
        !% Gives the possibility to change the partition nodes.
        !% Afterward, it crashes.
        !%End
        call parse_variable(namespace, 'MeshPartitionVirtualSize', mesh%mpi_grp%size, vsize)

        if (vsize /= mesh%mpi_grp%size) then
          write(message(1),'(a,I7)') "Changing the partition size to", vsize
          write(message(2),'(a)') "The execution will crash."
          call messages_warning(2, namespace=namespace)
          has_virtual_partition = .true.
        else
          has_virtual_partition = .false.
        end if

        if (.not. present(parent)) then
          call mesh_partition(mesh, namespace, stencil, vsize)
        else
          ! if there is a parent grid, use its partition
          call mesh_partition_from_parent(mesh, parent)
        end if

        !Now that we have the partitions, we save them
        call restart_init(restart_dump, namespace, RESTART_PARTITION, RESTART_TYPE_DUMP, mc, ierr, mesh=mesh)
        call mesh_partition_dump(restart_dump, mesh, vsize, ierr)
        call restart_end(restart_dump)
      end if

      if (has_virtual_partition) then
        call profiling_end(namespace)
        call print_date("Calculation ended on ")
        write(message(1),'(a)') "Execution has ended."
        write(message(2),'(a)') "If you want to run your system, do not use MeshPartitionVirtualSize."
        call messages_warning(2, namespace=namespace)
        call messages_end()
        call global_end()
        stop
      end if

      !%Variable MeshUseTopology
      !%Type logical
      !%Default false
      !%Section Execution::Parallelization
      !%Description
      !% (experimental) If enabled, <tt>Octopus</tt> will use an MPI virtual
      !% topology to map the processors. This can improve performance
      !% for certain interconnection systems.
      !%End
      call parse_variable(namespace, 'MeshUseTopology', .false., use_topo)

      if (use_topo) then
        ! At the moment we still need the global partition. This will be removed in near future.
        SAFE_ALLOCATE(part_vec(1:mesh%np_part_global))
        call partition_get_global(mesh%partition, part_vec(1:mesh%np_global))


        ! generate a table of neighbours

        SAFE_ALLOCATE(nb(1:mesh%mpi_grp%size, 1:mesh%mpi_grp%size))
        nb = .false.

        do ipg = 1, mesh%np_global
          ipart = part_vec(ipg)
          call mesh_global_index_to_coords(mesh, ipg, idx)
          do jj = 1, stencil%size
            jx(1:MAX_DIM) = idx(1:MAX_DIM) + stencil%points(1:MAX_DIM, jj)
            if (all(jx(1:MAX_DIM) >= mesh%idx%nr(1, 1:MAX_DIM)) .and. all(jx(1:MAX_DIM) <= mesh%idx%nr(2, 1:MAX_DIM))) then
              jpart = part_vec(mesh_global_index_from_coords(mesh, jx))
              if (ipart /= jpart ) nb(ipart, jpart) = .true.
            end if
          end do
        end do
        SAFE_DEALLOCATE_A(part_vec)

        ! now generate the information of the graph

        SAFE_ALLOCATE(gindex(1:mesh%mpi_grp%size))
        SAFE_ALLOCATE(gedges(1:count(nb)))

        ! and now generate it
        iedge = 0
        do ipart = 1, mesh%mpi_grp%size
          do jpart = 1, mesh%mpi_grp%size
            if (nb(ipart, jpart)) then
              iedge = iedge + 1
              gedges(iedge) = jpart - 1
            end if
          end do
          gindex(ipart) = iedge
        end do

        ASSERT(iedge == count(nb))

        reorder = .true.
        call MPI_Graph_create(mesh%mpi_grp%comm, mesh%mpi_grp%size, gindex, gedges, reorder, graph_comm, mpi_err)

        ! we have a new communicator
        call mpi_grp_init(mesh%mpi_grp, graph_comm)

        SAFE_DEALLOCATE_A(nb)
        SAFE_DEALLOCATE_A(gindex)
        SAFE_DEALLOCATE_A(gedges)

      end if

      call par_vec_init(mesh%mpi_grp, mesh%np_global, mesh%idx, stencil,&
        space, mesh%partition, mesh%pv, namespace)

      ! check the number of ghost neighbours in parallel
      nnb = 0
      jpart =  mesh%pv%partno
      do ipart = 1, mesh%pv%npart
        if (ipart == jpart) cycle
        if (mesh%pv%ghost_scounts(ipart) /= 0) nnb = nnb + 1
      end do
      ASSERT(nnb >= 0 .and. nnb < mesh%pv%npart)

      ! Set local point numbers.
      mesh%np      = mesh%pv%np_local
      mesh%np_part = mesh%np + mesh%pv%np_ghost + mesh%pv%np_bndry

      !%Variable PartitionPrint
      !%Type logical
      !%Default true
      !%Section Execution::Parallelization
      !%Description
      !% (experimental) If disabled, <tt>Octopus</tt> will not compute
      !% nor print the partition information, such as local points,
      !% no. of neighbours, ghost points and boundary points.
      !%End
      call parse_variable(namespace, 'PartitionPrint', .true., partition_print)

      if (partition_print) then
        call mesh_partition_write_info(mesh, namespace=namespace)
        call mesh_partition_messages_debug(mesh, namespace)
      end if
#endif

      POP_SUB(mesh_init_stage_3.do_partition)
    end subroutine do_partition


    ! ---------------------------------------------------------
    !> calculate the volume of integration
    subroutine mesh_get_vol_pp()

      integer :: jj(1:MAX_DIM), ip, np
      FLOAT   :: chi(MAX_DIM)

      PUSH_SUB(mesh_init_stage_3.mesh_get_vol_pp)

      np = 1
      if (mesh%use_curvilinear) np = mesh%np_part

      SAFE_ALLOCATE(mesh%vol_pp(1:np))

      do ip = 1, np
        mesh%vol_pp(ip) = product(mesh%spacing(1:space%dim))
      end do

      jj(space%dim + 1:MAX_DIM) = 0

      do ip = 1, np
        call mesh_local_index_to_coords(mesh, ip, jj)
        chi(1:space%dim) = jj(1:space%dim)*mesh%spacing(1:space%dim)
        mesh%vol_pp(ip) = mesh%vol_pp(ip)*mesh%coord_system%det_Jac(mesh%x(ip, 1:space%dim), chi(1:space%dim))
      end do

      if (mesh%use_curvilinear) then
        mesh%volume_element = M_ONE
      else
        mesh%volume_element = mesh%vol_pp(1)
      end if

      POP_SUB(mesh_init_stage_3.mesh_get_vol_pp)
    end subroutine mesh_get_vol_pp

  end subroutine mesh_init_stage_3

  subroutine rebalance_array(data_input, data_output, output_sizes)
    integer(i8),                    intent(in)  :: data_input(:)
    integer(i8), allocatable,       intent(out) :: data_output(:)
    integer, allocatable, optional, intent(out) :: output_sizes(:)

    integer, allocatable :: initial_sizes(:), final_sizes(:)
    integer(i8), allocatable :: initial_offsets(:), final_offsets(:)
    integer, allocatable :: scounts(:), sdispls(:)
    integer, allocatable :: rcounts(:), rdispls(:)
    integer :: irank
    integer(i8) :: itmp

    PUSH_SUB(rebalance_array)

    ! collect current sizes of distributed array
    SAFE_ALLOCATE(initial_sizes(0:mpi_world%size-1))
    call mpi_world%allgather(size(data_input), 1, MPI_INTEGER, initial_sizes(0), 1, MPI_INTEGER)
    SAFE_ALLOCATE(initial_offsets(0:mpi_world%size))
    initial_offsets(0) = 0
    do irank = 1, mpi_world%size
      initial_offsets(irank) = initial_offsets(irank-1) + initial_sizes(irank-1)
    end do

    ! now redistribute the arrays
    ! use block data decomposition of grid indices
    SAFE_ALLOCATE(final_offsets(0:mpi_world%size))
    SAFE_ALLOCATE(final_sizes(0:mpi_world%size-1))

    do irank = 0, mpi_world%size
      final_offsets(irank) = sum(int(initial_sizes, i8)) * irank/mpi_world%size
    end do
    do irank = 0, mpi_world%size - 1
      ASSERT(final_offsets(irank + 1) - final_offsets(irank) < huge(0_i4))
      final_sizes(irank) = int(final_offsets(irank + 1) - final_offsets(irank), i4)
    end do

    SAFE_ALLOCATE(scounts(0:mpi_world%size-1))
    SAFE_ALLOCATE(sdispls(0:mpi_world%size-1))
    SAFE_ALLOCATE(rcounts(0:mpi_world%size-1))
    SAFE_ALLOCATE(rdispls(0:mpi_world%size-1))
    ! determine communication pattern
    scounts = 0
    do irank = 0, mpi_world%size - 1
      ! get overlap of initial and final distribution
      itmp = min(final_offsets(irank+1), initial_offsets(mpi_world%rank+1)) - &
        max(final_offsets(irank), initial_offsets(mpi_world%rank))
      ASSERT(itmp < huge(0_i4))
      if (itmp < 0) then
        scounts(irank) = 0
      else
        scounts(irank) = int(itmp, i4)
      end if
    end do
    sdispls(0) = 0
    do irank = 1, mpi_world%size - 1
      sdispls(irank) = sdispls(irank - 1) + scounts(irank - 1)
    end do
    ASSERT(sum(int(scounts, i8)) < huge(0_i4))
    ASSERT(sum(scounts) == initial_sizes(mpi_world%rank))

    rcounts = 0
    do irank = 0, mpi_world%size - 1
      ! get overlap of initial and final distribution
      itmp = min(final_offsets(mpi_world%rank+1), initial_offsets(irank+1)) - &
        max(final_offsets(mpi_world%rank), initial_offsets(irank))
      ASSERT(itmp < huge(0_i4))
      if (itmp < 0) then
        rcounts(irank) = 0
      else
        rcounts(irank) = int(itmp, i4)
      end if
    end do
    rdispls(0) = 0
    do irank = 1, mpi_world%size - 1
      rdispls(irank) = rdispls(irank - 1) + rcounts(irank - 1)
    end do
    ! check for consistency between sending and receiving
    ASSERT(sum(rcounts) == final_sizes(mpi_world%rank))

    SAFE_ALLOCATE(data_output(1:final_sizes(mpi_world%rank)))
    call mpi_world%alltoallv(data_input, scounts, sdispls, MPI_INTEGER8, &
      data_output, rcounts, rdispls, MPI_INTEGER8)

    ! save final sizes of array if optional argument present
    if (present(output_sizes)) then
      SAFE_ALLOCATE(output_sizes(0:mpi_world%size-1))
      output_sizes(:) = final_sizes(:)
    end if

    POP_SUB(rebalance_array)
  end subroutine rebalance_array

  subroutine reorder_points(namespace, space, idx, grid_to_spatial, grid_to_spatial_reordered)
    type(namespace_t),        intent(in)  :: namespace
    type(space_t),            intent(in)  :: space
    type(index_t),            intent(in)  :: idx
    integer(i8),              intent(in)  :: grid_to_spatial(:)
    integer(i8), allocatable, intent(out) :: grid_to_spatial_reordered(:)

    integer :: bsize(space%dim), order, default
    integer :: nn, idir, ipg, ip, number_of_blocks(space%dim)
    type(block_t) :: blk
    integer, parameter :: &
      ORDER_BLOCKS     =  1, &
      ORDER_ORIGINAL    =  2, &
      ORDER_CUBE       =  3
    integer :: point(1:space%dim)
    integer(i8), allocatable :: reorder_indices(:), reorder_recv(:)
    integer, allocatable :: index_map(:), indices(:)
    integer(i8), allocatable :: grid_to_spatial_recv(:)
    integer, allocatable :: initial_sizes(:)
    integer(i8), allocatable :: initial_offsets(:)
    integer(i8) :: istart, iend, indstart, indend, spatial_size
    integer :: irank, local_size, num_recv
    integer :: iunique, nunique
    integer :: direction
    logical :: increase_with_dimension

    integer, allocatable :: scounts(:), sdispls(:), rcounts(:), rdispls(:)
    integer(i8), allocatable :: spatial_cutoff(:)

    PUSH_SUB(reorder_points)

    !%Variable MeshOrder
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% This variable controls how the grid points are mapped to a
    !% linear array for global arrays. For runs that are parallel
    !% in domains, the local mesh order may be different (see
    !% <tt>MeshLocalOrder</tt>).
    !% The default is blocks when serial in domains and cube when
    !% parallel in domains with the local mesh order set to blocks.
    !%Option order_blocks 1
    !% The grid is mapped using small parallelepipedic grids. The size
    !% of the blocks is controlled by <tt>MeshBlockSize</tt>.
    !%Option order_original 2
    !% The original order of the indices is used to map the grid.
    !%Option order_cube 3
    !% The grid is mapped using a full cube, i.e. without blocking.
    !%End
    default = ORDER_BLOCKS
    call parse_variable(namespace, 'MeshOrder', default, order)
    ! no reordering in 1D necessary
    if (space%dim == 1) then
      order = ORDER_ORIGINAL
    end if

    !%Variable MeshBlockDirection
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% Determines the direction in which the dimensions are chosen to compute
    !% the blocked index for sorting the mesh points (see MeshBlockSize).
    !% The default is increase_with_dimensions, corresponding to xyz ordering
    !% in 3D.
    !%Option increase_with_dimension 1
    !% The fastest changing index is in the first dimension, i.e., in 3D this
    !% corresponds to ordering in xyz directions.
    !%Option decrease_with_dimension 2
    !% The fastest changing index is in the last dimension, i.e., in 3D this
    !% corresponds to ordering in zyx directions.
    !%End
    call parse_variable(namespace, 'MeshBlockDirection', 1, direction)
    increase_with_dimension = direction == 1
    if (direction /= 1 .and. direction /= 2) then
      call messages_input_error(namespace, 'MeshBlockDirection')
    end if

    select case (order)
    case (ORDER_ORIGINAL)
      ! only copy points, they stay in their original ordering
      SAFE_ALLOCATE(grid_to_spatial_reordered(1:size(grid_to_spatial)))
      grid_to_spatial_reordered(1:size(grid_to_spatial)) = grid_to_spatial(1:size(grid_to_spatial))
    case (ORDER_BLOCKS, ORDER_CUBE)
      if (order == ORDER_CUBE) then
        bsize(1:space%dim) = idx%nr(2, 1:space%dim) - idx%nr(1, 1:space%dim) + 1
      else
        !%Variable MeshBlockSize
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

        if (parse_block(namespace, 'MeshBlockSize', blk) == 0) then
          nn = parse_block_cols(blk, 0)
          if (nn /= space%dim) then
            message(1) = "Error: number of entries in MeshBlockSize must match the number of dimensions."
            call messages_fatal(1, namespace=namespace)
          end if
          do idir = 1, nn
            call parse_block_integer(blk, 0, idir - 1, bsize(idir))
          end do
        end if
      end if

      number_of_blocks(1:space%dim) = (idx%nr(2, 1:space%dim) - idx%nr(1, 1:space%dim) + 1) &
        /bsize(1:space%dim) + 1


      ! do the global reordering in parallel, use block data decomposition of global indices
      ! reorder indices along blocked parallelepiped curve

      ! collect current sizes of distributed array
      SAFE_ALLOCATE(initial_sizes(0:mpi_world%size-1))
      call mpi_world%allgather(size(grid_to_spatial), 1, MPI_INTEGER, initial_sizes(0), 1, MPI_INTEGER)
      SAFE_ALLOCATE(initial_offsets(0:mpi_world%size))
      initial_offsets(0) = 0
      do irank = 1, mpi_world%size
        initial_offsets(irank) = initial_offsets(irank-1) + initial_sizes(irank-1)
      end do

      ! get local range and size
      istart = initial_offsets(mpi_world%rank)
      iend = initial_offsets(mpi_world%rank + 1) - 1
      ASSERT(iend - istart + 1 < huge(0_i4))
      local_size = int(iend - istart + 1, i4)
      ASSERT(local_size == initial_sizes(mpi_world%rank))

      ! compute new indices locally
      SAFE_ALLOCATE(reorder_indices(1:local_size))
      SAFE_ALLOCATE(indices(1:local_size))
      SAFE_ALLOCATE(grid_to_spatial_reordered(1:local_size))
      !$omp parallel do private(point)
      do ip = 1, local_size
        call index_spatial_to_point(idx, space%dim, grid_to_spatial(ip), point)
        point(1:space%dim) = point(1:space%dim) + idx%offset(1:space%dim)
        reorder_indices(ip) = get_blocked_index(space%dim, point, bsize, number_of_blocks, increase_with_dimension)
      end do
      ! parallel sort according to the new indices
      ! sort the local array
      call sort(reorder_indices, indices)
      ! save reordered indices to send to other processes
      !$omp parallel do
      do ip = 1, local_size
        grid_to_spatial_reordered(ip) = grid_to_spatial(indices(ip))
      end do

      ! get minimum and maximum
      indstart = reorder_indices(1)
      indend = reorder_indices(local_size)
      call mpi_world%allreduce_inplace(indstart, 1, MPI_INTEGER8, MPI_MIN)
      call mpi_world%allreduce_inplace(indend, 1, MPI_INTEGER8, MPI_MAX)
      spatial_size = indend - indstart + 1

      ! get index ranges for each rank
      SAFE_ALLOCATE(spatial_cutoff(0:mpi_world%size-1))
      do irank = 0, mpi_world%size - 1
        spatial_cutoff(irank) = spatial_size * (irank+1)/mpi_world%size + indstart
      end do

      SAFE_ALLOCATE(scounts(0:mpi_world%size-1))
      SAFE_ALLOCATE(sdispls(0:mpi_world%size-1))
      SAFE_ALLOCATE(rcounts(0:mpi_world%size-1))
      SAFE_ALLOCATE(rdispls(0:mpi_world%size-1))
      ! get send counts
      scounts = 0
      irank = 0
      ! the indices are ordered, so we can go through them and increase
      ! the rank to which they are associated to when we cross a cutoff
      do ip = 1, local_size
        if (reorder_indices(ip) >= spatial_cutoff(irank)) then
          ! this do loop is needed in case some ranks do not have any points
          do while (reorder_indices(ip) >= spatial_cutoff(irank))
            irank = irank + 1
          end do
          ASSERT(irank < mpi_world%size)
        end if
        scounts(irank) = scounts(irank) + 1
      end do
      SAFE_DEALLOCATE_A(spatial_cutoff)
      ASSERT(sum(scounts) == local_size)

      ! compute communication pattern (sdispls, rcounts, rdispls)
      sdispls(0) = 0
      do irank = 1, mpi_world%size - 1
        sdispls(irank) = sdispls(irank - 1) + scounts(irank - 1)
      end do

      call mpi_world%alltoall(scounts, 1, MPI_INTEGER, &
        rcounts, 1, MPI_INTEGER)

      rdispls(0) = 0
      do irank = 1, mpi_world%size - 1
        rdispls(irank) = rdispls(irank - 1) + rcounts(irank - 1)
      end do

      ! make sure the arrays get allocated also if we do not receive anything
      num_recv = max(sum(rcounts), 1)
      ! communicate the locally sorted indices
      SAFE_ALLOCATE(reorder_recv(1:num_recv))
      call mpi_world%alltoallv(reorder_indices, scounts, sdispls, MPI_INTEGER8, &
        reorder_recv, rcounts, rdispls, MPI_INTEGER8)
      SAFE_DEALLOCATE_A(reorder_indices)

      ! communicate the corresponding spatial indices
      SAFE_ALLOCATE(grid_to_spatial_recv(1:num_recv))
      call mpi_world%alltoallv(grid_to_spatial_reordered, scounts, sdispls, MPI_INTEGER8, &
        grid_to_spatial_recv, rcounts, rdispls, MPI_INTEGER8)
      SAFE_DEALLOCATE_A(grid_to_spatial_reordered)

      ! do k-way merge of sorted indices
      SAFE_ALLOCATE(reorder_indices(1:num_recv))
      SAFE_ALLOCATE(index_map(1:num_recv))
      if (sum(rcounts) > 0) then
        call merge_sorted_arrays(reorder_recv, rcounts, reorder_indices, index_map)

        ! get number of unique indices, needed for boundary
        nunique = 1
        do ipg = 2, sum(rcounts)
          if (reorder_indices(ipg) /= reorder_indices(ipg-1)) then
            nunique = nunique + 1
          end if
        end do

        ! reorder according to new order, but remove duplicate entries
        SAFE_ALLOCATE(grid_to_spatial_reordered(1:nunique))
        iunique = 1
        grid_to_spatial_reordered(iunique) = grid_to_spatial_recv(index_map(1))
        do ipg = 2, sum(rcounts)
          if (reorder_indices(ipg) /= reorder_indices(ipg-1)) then
            iunique = iunique + 1
            grid_to_spatial_reordered(iunique) = grid_to_spatial_recv(index_map(ipg))
          end if
        end do
      else
        SAFE_ALLOCATE(grid_to_spatial_reordered(1:0))
      end if

      SAFE_DEALLOCATE_A(initial_offsets)
      SAFE_DEALLOCATE_A(initial_sizes)

      SAFE_DEALLOCATE_A(reorder_indices)
      SAFE_DEALLOCATE_A(reorder_recv)

      SAFE_DEALLOCATE_A(grid_to_spatial_recv)
      SAFE_DEALLOCATE_A(index_map)
      SAFE_DEALLOCATE_A(indices)

      SAFE_DEALLOCATE_A(scounts)
      SAFE_DEALLOCATE_A(sdispls)
      SAFE_DEALLOCATE_A(rcounts)
      SAFE_DEALLOCATE_A(rdispls)

    end select
    POP_SUB(reorder_points)
  end subroutine reorder_points

  subroutine get_sizes_offsets(local_size, sizes, offsets, mpi_grp)
    integer,                   intent(in)  :: local_size
    integer, allocatable,      intent(out) :: sizes(:)
    integer(i8), allocatable,  intent(out) :: offsets(:)
    type(mpi_grp_t), optional, intent(in)  :: mpi_grp

    integer :: irank
    type(mpi_grp_t) :: mpi_grp_

    PUSH_SUB(get_sizes_offsets)

    if (present(mpi_grp)) then
      mpi_grp_ = mpi_grp
    else
      mpi_grp_ = mpi_world
    end if

    SAFE_ALLOCATE(sizes(0:mpi_grp_%size-1))
    call mpi_grp_%allgather(local_size, 1, MPI_INTEGER, sizes(0), 1, MPI_INTEGER)
    SAFE_ALLOCATE(offsets(0:mpi_grp_%size))
    offsets(0) = 0
    do irank = 1, mpi_grp_%size
      offsets(irank) = offsets(irank-1) + sizes(irank-1)
    end do

    POP_SUB(get_sizes_offsets)
  end subroutine get_sizes_offsets

end module mesh_init_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
