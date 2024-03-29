!! Copyright (C) 2002-2013 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module mesh_partition_oct_m
  use debug_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use mesh_oct_m
  use metis_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use partition_oct_m
  use profiling_oct_m
  use restart_oct_m
  use stencil_oct_m
  use stencil_star_oct_m

  implicit none

  private
  public ::                      &
    mesh_partition,              &
    mesh_partition_from_parent,  &
    mesh_partition_dump,         &
    mesh_partition_load,         &
    mesh_partition_write_info,   &
    mesh_partition_messages_debug

  integer, parameter :: &
    METIS    = 1,       &
    PARMETIS = 2,       &
    HILBERT  = 3

  integer, parameter :: &
    RCB    = 1,         &
    GRAPH  = 2

contains

  ! ---------------------------------------------------------------
  !> Converts the mesh given by grid points into a graph. Each
  !! point is a vertex in the graph and closest neighbours are
  !! connected by an edge (at most 6 in 3D and 4 in 2D, 2 in
  !! 1D, fewer at the boundaries).
  !! Then calls METIS to get npart partitions.
  !! Stored the mapping point no. -> partition no. into part,
  !! which has to be allocated beforehand.
  !! (mesh_partition_end should be called later.)
  ! ---------------------------------------------------------------
  subroutine mesh_partition(mesh, namespace, lapl_stencil, vsize)
    type(mesh_t),      intent(inout)  :: mesh
    type(namespace_t), intent(in)     :: namespace
    type(stencil_t),   intent(in)     :: lapl_stencil
    integer,           intent(in)     :: vsize        !< number of partitions to be created. Might
    !!                                                   be different from the actual number of
    !!                                                   domains

    integer              :: jp, jpart
    integer(imetis)      :: iv
    integer(i8)          :: ivtmp, inb
    integer              :: ix(1:MAX_DIM), jx(1:MAX_DIM)

    integer              :: npart          !< Number of partitions.
    integer              :: ipart          !< number of the current partition

    integer(imetis)              :: ne_global !< Global number of edges.
    integer(imetis)              :: nv_global !< Global number of vertices.
    !! The global number of vertices (nv_global) is equal to the number of
    !! points np_global and the maximum global number of edges (ne_global)
    !! is 2*mesh%box%dim*np_global (there are a little fewer because points
    !! on the border have fewer than two neighbours per dimension).
    !! xadj_global has nv_global+1 entries because last entry contains the total
    !! number of edges.
    integer(imetis)              :: ne        !< Local number of edges.
    integer                      :: nv        !< Local number of vertices.

    integer(imetis), allocatable :: xadj_global(:)   !< Indices of adjacency list in adjncy_global.
    integer(imetis), allocatable :: adjncy_global(:) !< Adjacency lists.
    integer(imetis), allocatable :: adjncy(:)        !< Local part of adjacency list
    integer(imetis), allocatable :: xadj(:)          !< Local part of xadj

    integer(imetis), allocatable :: options(:)     !< Options to (Par)METIS.
#ifdef HAVE_METIS
    integer(imetis)              :: edgecut        !< Number of edges cut by partitioning.
#endif
    REAL_SINGLE, allocatable :: tpwgts(:)  !< The fraction of vertex weight that should be distributed

    integer                      :: iunit          !< For debug output to files.
    integer, allocatable :: rcounts(:)
    integer(imetis)      :: tmp
    integer(imetis), allocatable :: rdispls(:)
    integer(imetis), allocatable :: part(:)        !< The local partition
    integer(imetis), allocatable :: part_global(:) !< The global partition (only used with METIS)
    integer(imetis), allocatable :: vtxdist(:)     !< Initial distribution of the points

    type(stencil_t) :: stencil
    integer :: stencil_to_use, default_method, method
    integer :: library
    integer, parameter   :: STAR = 1, LAPLACIAN = 2
    integer(imetis), allocatable :: istart(:)
    integer, allocatable :: lsize(:)

    type(profile_t), save :: prof
    integer :: default, ierr

    call profiling_in(prof, "MESH_PARTITION")
    PUSH_SUB(mesh_partition)

    if (mesh%np_global == 0) then
      message(1) = 'The mesh is empty and cannot be partitioned.'
      call messages_fatal(1, namespace=namespace)
    end if

    !%Variable MeshPartitionPackage
    !%Type integer
    !%Section Execution::Parallelization
    !%Description
    !% Decides which library to use to perform the mesh partition.
    !% By default ParMETIS is used when available, otherwise METIS is used.
    !%Option metis 1
    !% METIS library.
    !%Option parmetis 2
    !% (Experimental) Use ParMETIS library to perform the mesh partition.
    !% Only available if the code was compiled with ParMETIS support.
    !%Option part_hilbert 3
    !% Use the ordering along the Hilbert curve for partitioning.
    !%End
    default = METIS
#ifdef HAVE_PARMETIS
    default = PARMETIS
#endif
    call parse_variable(namespace, 'MeshPartitionPackage', default, library)

#if !defined(HAVE_METIS)
    if (library == METIS) then
      message(1) = 'METIS was requested, but Octopus was compiled without it.'
      call messages_fatal(1, only_root_writes = .true.)
    end if
#endif
#if !defined(HAVE_PARMETIS)
    if (library == PARMETIS) then
      message(1) = 'PARMETIS was requested, but Octopus was compiled without it.'
      call messages_fatal(1, only_root_writes = .true.)
    end if
#endif

    !%Variable MeshPartitionStencil
    !%Type integer
    !%Default stencil_star
    !%Section Execution::Parallelization
    !%Description
    !% To partition the mesh, it is necessary to calculate the connection
    !% graph connecting the points. This variable selects which stencil
    !% is used to do this.
    !%Option stencil_star 1
    !% An order-one star stencil.
    !%Option laplacian 2
    !% The stencil used for the Laplacian is used to calculate the
    !% partition. This in principle should give a better partition, but
    !% it is slower and requires more memory.
    !%End
    call parse_variable(namespace, 'MeshPartitionStencil', STAR, stencil_to_use)

    if (stencil_to_use == STAR) then
      call stencil_star_get_lapl(stencil, mesh%box%dim, order = 1)
    else if (stencil_to_use == LAPLACIAN) then
      call stencil_copy(lapl_stencil, stencil)
    else
      call messages_input_error(namespace, 'MeshPartitionStencil')
    end if

    ! Shortcut to the global number of vertices
    if (mesh%np_global > huge(nv_global)) then
      message(1) = "Error: too many grid points for this version of metis/parmetis."
      message(2) = "Please use the 64-bit version."
      call messages_fatal(2, namespace=namespace)
    end if
    nv_global = i8_to_imetis(mesh%np_global)

    call partition_init(mesh%partition, imetis_to_i8(nv_global), mesh%mpi_grp)

    ! Get number of partitions.
    npart = mesh%mpi_grp%size
    ipart = mesh%mpi_grp%rank + 1

    SAFE_ALLOCATE(istart(1:npart))
    SAFE_ALLOCATE(lsize(1:npart))
    SAFE_ALLOCATE(vtxdist(1:npart+1))
    call partition_get_local_size(mesh%partition, ivtmp, nv)
    ASSERT(ivtmp <= huge(0_imetis))
    iv = i8_to_imetis(ivtmp)

    ! Allocate local output matrix
    SAFE_ALLOCATE(part(1:nv))
    if (library == HILBERT) then
      part = ipart
    else

      call mesh%mpi_grp%allgather(iv, 1, MPI_METIS_INT, istart(1), 1, MPI_METIS_INT)
      call mesh%mpi_grp%allgather(nv, 1, MPI_INTEGER, lsize(1), 1, MPI_INTEGER)

      vtxdist(1:npart) = istart(1:npart)
      vtxdist(npart+1) = nv_global + 1

      ! Create graph with each point being represented by a
      ! vertex and edges between neighbouring points.
      SAFE_ALLOCATE(xadj(1:nv+1))
      SAFE_ALLOCATE(adjncy(1:i4_to_imetis(stencil%size - 1)*nv))
      ne = 1
      ! Iterate over number of vertices.
      do iv = 1, lsize(ipart)
        ! Get coordinates of point iv (vertex iv).
        call mesh_global_index_to_coords(mesh, imetis_to_i8(istart(ipart)-1 + iv), ix)

        ! Set entry in index table.
        xadj(iv) = ne
        ! Check all possible neighbours.
        do jp = 1, stencil%size
          if (jp == stencil%center) cycle

          ! Store coordinates of possible neighbors, they
          ! are needed several times in the check below.
          jx(1:MAX_DIM) = ix(1:MAX_DIM) + stencil%points(1:MAX_DIM, jp)

          if (all(jx(1:MAX_DIM) >= mesh%idx%nr(1, 1:MAX_DIM)) .and. all(jx(1:MAX_DIM) <= mesh%idx%nr(2, 1:MAX_DIM))) then
            ! Only points inside the mesh or its enlargement
            ! are included in the graph.
            inb = mesh_global_index_from_coords(mesh, jx)
            if (inb /= 0 .and. inb <= nv_global) then
              ! Store a new edge and increment edge counter.
              adjncy(ne) = i8_to_imetis(inb)
              ne = ne + 1
            end if
          end if
        end do
      end do
      ! We started with ne=1 for simplicity, so ne is off by one in the end.
      ne = ne - 1

      ! Set the total number of edges plus one as last index.
      ! (NOTE: the plus one is because we are using Fortran numbering)
      ! The reason is: neighbours of node i are stored in adjncy(xadj(i):xadj(i+1)-1).
      ! Setting the last index as mentioned makes special handling of last element
      ! unnecessary (this indexing is a METIS requirement).
      xadj(nv+1) = ne + 1


      if (debug%info .or. library == METIS) then
        !Gather the global xadj and adjncy arrays
        SAFE_ALLOCATE(rcounts(1:npart))
        SAFE_ALLOCATE(rdispls(1:npart))

        SAFE_ALLOCATE(xadj_global(1:nv_global + 1))
        xadj_global(1) = 1
        rcounts(1:npart) = lsize(1:npart)
        rdispls(1:npart) = vtxdist(1:npart) - 1
        call mesh%mpi_grp%allgatherv(xadj(2:), lsize(ipart), MPI_METIS_INT, &
          xadj_global(2:), rcounts, rdispls, MPI_METIS_INT)
        do jpart = 2, npart
          do iv = 1, lsize(jpart)
            xadj_global(istart(jpart) + iv) = xadj_global(istart(jpart) + iv) + xadj_global(istart(jpart)) - 1
          end do
        end do

        ne_global = xadj_global(nv_global + 1)
        SAFE_ALLOCATE(adjncy_global(1:ne_global))
        do jpart = 1, npart
          rdispls(jpart) = xadj_global(vtxdist(jpart)) - 1
          tmp = xadj_global(vtxdist(jpart+1)) - 1 - rdispls(jpart)
          ASSERT(tmp < huge(0_i4))
          rcounts(jpart) = imetis_to_i4(tmp)
        end do
        ASSERT(xadj(nv+1)-1 < huge(0_i4))
        call mesh%mpi_grp%allgatherv(adjncy, imetis_to_i4(xadj(nv+1)-1), MPI_METIS_INT, &
          adjncy_global, rcounts, rdispls, MPI_METIS_INT)

        SAFE_DEALLOCATE_A(rcounts)
        SAFE_DEALLOCATE_A(rdispls)
      end if


      if (debug%info) then
        ! DEBUG output. Write graph to file mesh_graph.txt.
        message(1) = 'Info: Adjacency lists of the graph representing the grid'
        message(2) = 'Info: are stored in debug/mesh_partition/mesh_graph.txt.'
        message(3) = 'Info: Compatible with METIS programs pmetis and kmetis.'
        message(4) = 'Info: First line contains number of vertices and edges.'
        message(5) = 'Info: Edges are not directed and appear twice in the lists.'
        call messages_info(5)
        if (mpi_grp_is_root(mpi_world)) then
          call io_mkdir('debug/mesh_partition', namespace)
          iunit = io_open('debug/mesh_partition/mesh_graph.txt', namespace, action='write')
          write(iunit, *) nv_global, ne_global/2
          do iv = 1, nv
            write(iunit, *) adjncy_global(xadj_global(iv):xadj_global(iv+1) - 1)
          end do
          call io_close(iunit)
        end if
        message(1) = "Info: Done writing mesh_graph.txt."
        call messages_info(1, namespace=namespace)
      end if

      ! The sum of all tpwgts elements has to be 1 and
      ! we don`t care about the weights. So; 1/npart
      SAFE_ALLOCATE(tpwgts(1:vsize))
      tpwgts(1:vsize) = real(1.0, 4)/real(vsize, 4)

      select case (library)
      case (METIS)
        SAFE_ALLOCATE(options(1:40))
        options = 0
#ifdef HAVE_METIS
        call oct_metis_setdefaultoptions(options(1)) ! is equal to: options = -1
#endif
        options(METIS_OPTION_NUMBERING) = 1 ! Fortran style: start counting from 1

        if (vsize  <  8) then
          default_method = RCB
        else
          default_method = GRAPH
        end if

        !%Variable MeshPartition
        !%Type integer
        !%Section Execution::Parallelization
        !%Description
        !% When using METIS to perform the mesh partitioning, decides which
        !% algorithm is used. By default, <tt>graph</tt> partitioning
        !% is used for 8 or more partitions, and <tt>rcb</tt> for fewer.
        !%Option rcb 1
        !% Recursive coordinate bisection partitioning.
        !%Option graph 2
        !% Graph partitioning (called 'k-way' by METIS).
        !%End
        call parse_variable(namespace, 'MeshPartition', default_method, method)

        SAFE_ALLOCATE(part_global(1:nv_global))

        !Now we can call METIS
        select case (method)
        case (RCB)
          message(1) = 'Info: Using METIS 5 multilevel recursive bisection to partition the mesh.'
          call messages_info(1, namespace=namespace)
#ifdef HAVE_METIS
          ierr = oct_metis_partgraphrecursive(nv_global, 1_imetis, xadj_global(1), adjncy_global(1), &
            i4_to_imetis(vsize), tpwgts(1), 1.01_4, options(1), edgecut, part_global(1))
#endif
          call metis_error_code(ierr)
        case (GRAPH)
          message(1) = 'Info: Using METIS 5 multilevel k-way algorithm to partition the mesh.'
          call messages_info(1, namespace=namespace)
#ifdef HAVE_METIS
          ierr = oct_metis_partgraphkway(nv_global, 1_imetis, xadj_global(1), adjncy_global(1), &
            i4_to_imetis(vsize), tpwgts(1), 1.01_4, options(1), edgecut, part_global(1))
#endif
          call metis_error_code(ierr)
        case default
          message(1) = 'Selected partition method is not available in METIS 5.'
          call messages_fatal(1, namespace=namespace)
        end select

        part(1:nv) = part_global(istart(ipart):istart(ipart) + nv - 1)

        SAFE_DEALLOCATE_A(options)
        SAFE_DEALLOCATE_A(part_global)

      case (PARMETIS)

        SAFE_ALLOCATE(options(1:3))
        options = 0 ! For the moment we use default options

        write(message(1),'(a)') 'Info: Using ParMETIS multilevel k-way algorithm to partition the mesh.'
        call messages_info(1, namespace=namespace)

        ! Call to ParMETIS with no imbalance tolerance
#ifdef HAVE_PARMETIS
        call oct_parmetis_v3_partkway(vtxdist(1), xadj(1), adjncy(1), &
          1_imetis, i4_to_imetis(mesh%mpi_grp%size), tpwgts(1), 1.05_4, &
          options(1), edgecut, part(1), mesh%mpi_grp%comm)
#endif

        ! Deallocate matrices
        SAFE_DEALLOCATE_A(options)

      end select

      SAFE_DEALLOCATE_A(vtxdist)

      if (debug%info .or. library == METIS) then
        SAFE_DEALLOCATE_A(xadj_global)
        SAFE_DEALLOCATE_A(adjncy_global)
      end if
      SAFE_DEALLOCATE_A(tpwgts)
      SAFE_DEALLOCATE_A(xadj)
      SAFE_DEALLOCATE_A(adjncy)
      SAFE_DEALLOCATE_A(istart)
      SAFE_DEALLOCATE_A(lsize)
    end if

    ASSERT(all(part(1:nv) > 0))
    ASSERT(all(part(1:nv) <= vsize))
    call partition_set(mesh%partition, imetis_to_i4(part))

    SAFE_DEALLOCATE_A(part)

    call stencil_end(stencil)
    POP_SUB(mesh_partition)
    call profiling_out(prof)

  contains

    subroutine metis_error_code(ierr)
      integer, intent(in) :: ierr

      PUSH_SUB(mesh_partition.metis_error_code)

      select case (ierr)
      case (METIS_OK)
        !Everything OK
      case (METIS_ERROR_INPUT)
        message(1) = "Metis: Input error."
        call messages_fatal(1, namespace=namespace)
      case (METIS_ERROR_MEMORY)
        message(1) = "Metis: Unable to allocate required memory."
        call messages_fatal(1, namespace=namespace)
      case (METIS_ERROR)
        message(1) = "Metis: Some type of error."
        call messages_fatal(1, namespace=namespace)
      end select

      POP_SUB(mesh_partition.metis_error_code)
    end subroutine metis_error_code

  end subroutine mesh_partition

  ! ----------------------------------------------------
  subroutine mesh_partition_from_parent(mesh, parent)
    type(mesh_t), intent(inout) :: mesh
    type(mesh_t), intent(in)    :: parent

    integer :: np, ip_local
    integer(i8) :: istart, ip_global
    integer(i8), allocatable :: points(:)
    integer, allocatable :: part(:)
    integer :: idx(1:3)

    PUSH_SUB(mesh_partition_from_parent)

    !Initialize partition
    call partition_init(mesh%partition, mesh%np_global, mesh%mpi_grp)
    call partition_get_local_size(mesh%partition, istart, np)

    !Get the partition number of each point from the parent mesh
    SAFE_ALLOCATE(points(1:np))
    SAFE_ALLOCATE(part(1:np))
    do ip_local = 1, np
      ip_global = istart + ip_local - 1
      call mesh_global_index_to_coords(mesh, ip_global, idx)
      points(ip_local) = mesh_global_index_from_coords(parent, 2*idx)
    end do
    call partition_get_partition_number(parent%partition, np, points, part)

    !Set the partition
    call partition_set(mesh%partition, part)

    !Free memory
    SAFE_DEALLOCATE_A(points)
    SAFE_DEALLOCATE_A(part)

    POP_SUB(mesh_partition_from_parent)
  end subroutine mesh_partition_from_parent

  ! ----------------------------------------------------
  subroutine mesh_partition_dump(restart, mesh, vsize, ierr)
    type(restart_t), intent(in)  :: restart
    type(mesh_t),    intent(in)  :: mesh
    integer,         intent(in)  :: vsize
    integer,         intent(out) :: ierr

    integer :: err
    character(len=6) :: numstring
    type(profile_t), save :: prof

    PUSH_SUB(mesh_partition_dump)

    call profiling_in(prof, "PARTITION_WRITE")

    ierr = 0

    if (restart_skip(restart)) then
      call profiling_out(prof)
      POP_SUB(mesh_partition_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing mesh partition restart."
      call messages_info(1)
    end if

    write(numstring, '(i6.6)') vsize

    call partition_dump(mesh%partition, restart_dir(restart), 'inner_partition_'//trim(numstring)//'.obf', err)
    if (err /= 0) then
      message(1) = "Unable to write inner mesh partition to 'inner_partition_"//trim(numstring)//".obf'"
      call messages_warning(1)
      ierr = ierr + 1
    end if

    if (debug%info) then
      message(1) = "Debug: Writing mesh partition restart done."
      call messages_info(1)
    end if

    call profiling_out(prof)

    POP_SUB(mesh_partition_dump)
  end subroutine mesh_partition_dump

  ! ----------------------------------------------------
  subroutine mesh_partition_load(restart, mesh, ierr)
    type(restart_t), intent(in)    :: restart
    type(mesh_t),    intent(inout) :: mesh
    integer,         intent(out)   :: ierr

    integer :: err
    character(len=6) :: numstring
    character(len=MAX_PATH_LEN) :: filename, dir

    PUSH_SUB(mesh_partition_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(mesh_partition_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading mesh partition restart."
      call messages_info(1)
    end if

    call partition_init(mesh%partition, mesh%np_global, mesh%mpi_grp)

    write(numstring, '(i6.6)') mesh%mpi_grp%size
    dir = restart_dir(restart)

    !Read inner partition
    filename = 'inner_partition_'//numstring//'.obf'
    if (.not. io_file_exists(trim(dir)//"/"//filename)) then
      message(1) = "Unable to read inner mesh partition, file '"//trim(filename)//"' does not exist."
      call messages_warning(1)
      ierr = ierr + 1
    else
      call partition_load(mesh%partition, restart_dir(restart), filename, err)
      if (err /= 0) then
        message(1) = "Unable to read inner mesh partition from '"//trim(filename)//"'."
        call messages_warning(1)
        ierr = ierr + 2
      end if
    end if

    ! Free the memory in case we were unable to read the partitions
    if (ierr /= 0) then
      call partition_end(mesh%partition)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading mesh partition restart done."
      call messages_info(1)
    end if

    POP_SUB(mesh_partition_load)
  end subroutine mesh_partition_load


  ! ----------------------------------------------------------------------
  subroutine mesh_partition_write_info(mesh, iunit, namespace)
    type(mesh_t),                intent(in) :: mesh
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer(i8) :: npoints
    integer :: npart
    integer, allocatable :: nghost(:), nbound(:), nlocal(:), nneigh(:)
    integer :: num_ghost, num_bound, num_local, num_neigh
    FLOAT :: quality

    integer :: ipart
    type(profile_t), save :: prof
    logical, allocatable :: is_a_neigh(:)
    FLOAT :: scal

    PUSH_SUB(mesh_partition_write_info)

    call profiling_in(prof, "MESH_PARTITION_WRITE_INFO")

    ! Build information about the mesh partition
    npart = mesh%mpi_grp%size
    npoints = mesh%np_part_global

    SAFE_ALLOCATE(nghost(1:npart))
    SAFE_ALLOCATE(nbound(1:npart))
    SAFE_ALLOCATE(nlocal(1:npart))
    SAFE_ALLOCATE(nneigh(1:npart))
    SAFE_ALLOCATE(is_a_neigh(1:npart))

    ipart = mesh%mpi_grp%rank + 1

    num_local = mesh%np
    num_ghost = mesh%pv%np_ghost
    num_bound = mesh%pv%np_bndry

    is_a_neigh = mesh%pv%ghost_rcounts > 0
    num_neigh = count(is_a_neigh(1:npart))

    SAFE_DEALLOCATE_A(is_a_neigh)

    call mesh%mpi_grp%gather(num_neigh, 1, MPI_INTEGER, nneigh(1), 1, MPI_INTEGER, 0)
    call mesh%mpi_grp%gather(num_local, 1, MPI_INTEGER, nlocal(1), 1, MPI_INTEGER, 0)
    call mesh%mpi_grp%gather(num_ghost, 1, MPI_INTEGER, nghost(1), 1, MPI_INTEGER, 0)
    call mesh%mpi_grp%gather(num_bound, 1, MPI_INTEGER, nbound(1), 1, MPI_INTEGER, 0)

    ! Calculate partition quality
    scal = TOFLOAT(npart)/npoints

    quality = M_ZERO

    quality = quality + (TOFLOAT(maxval(nlocal)) - TOFLOAT(minval(nlocal)))**3
    quality = quality + (sum(TOFLOAT(nghost)**2))

    quality = M_ONE/(M_ONE + quality)


    ! Only the root has the information (nlocal, nghost, ...)
    if (mpi_grp_is_root(mesh%mpi_grp)) then
      ! Write information about the partition
      message(1) = &
        'Info: Mesh partition:'
      message(2) = ''
      call messages_info(2, iunit=iunit, namespace=namespace)

      write(message(1),'(a,e16.6)') &
        '      Partition quality:', quality
      message(2) = ''
      call messages_info(2, iunit=iunit, namespace=namespace)

      write(message(1),'(a)') &
        '                 Neighbours         Ghost points'
      write(message(2),'(a,i5,a,i10)') &
        '      Average  :      ', nint(sum(TOFLOAT(nneigh))/npart), '           ', nint(sum(TOFLOAT(nghost))/npart)
      write(message(3),'(a,i5,a,i10)') &
        '      Minimum  :      ', minval(nneigh),    '           ', minval(nghost)
      write(message(4),'(a,i5,a,i10)') &
        '      Maximum  :      ', maxval(nneigh),    '           ', maxval(nghost)
      message(5) = ''
      call messages_info(5, iunit=iunit, namespace=namespace)

      do ipart = 1, npart
        write(message(1),'(a,i5)')  &
          '      Nodes in domain-group  ', ipart
        write(message(2),'(a,i10,a,i10)') &
          '        Neighbours     :', nneigh(ipart), &
          '        Local points    :', nlocal(ipart)
        write(message(3),'(a,i10,a,i10)') &
          '        Ghost points   :', nghost(ipart), &
          '        Boundary points :', nbound(ipart)
        call messages_info(3, iunit=iunit, namespace=namespace)
      end do

      message(1) = ''
      call messages_info(1, iunit=iunit, namespace=namespace)
    end if

    SAFE_DEALLOCATE_A(nghost)
    SAFE_DEALLOCATE_A(nbound)
    SAFE_DEALLOCATE_A(nlocal)
    SAFE_DEALLOCATE_A(nneigh)

    call profiling_out(prof)

    POP_SUB(mesh_partition_write_info)
  end subroutine mesh_partition_write_info

  ! ----------------------------------------------------
  subroutine mesh_partition_messages_debug(mesh, namespace)
    type(mesh_t),      intent(in)    :: mesh
    type(namespace_t), intent(in)    :: namespace

    integer              :: ip
    integer(i8)          :: ipg
    integer              :: iunit          ! For debug output to files.
    character(len=6)     :: filenum

    if (.not. debug%info) return

    PUSH_SUB(mesh_partition_messages_debug)

    call io_mkdir('debug/mesh_partition', namespace)

    ! Debug output. Write points of each partition in a different file.
    write(filenum, '(i6.6)') mesh%mpi_grp%rank+1

    ! without boundary
    iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
      namespace, action='write')
    do ip = 1, mesh%np
      ipg = mesh_local2global(mesh, ip)
      write(iunit, '(i8,99f18.8)') ipg, mesh_x_global(mesh, ipg)
    end do
    call io_close(iunit)

    ! with boundary included
    iunit = io_open('debug/mesh_partition/mesh_partition_all.'//filenum, &
      namespace, action='write')
    do ip = 1, mesh%np_part
      ipg = mesh_local2global(mesh, ip)
      write(iunit, '(i8,99f18.8)') ipg, mesh_x_global(mesh, ipg)
    end do
    call io_close(iunit)

    ! points from enlargement
    if (mpi_grp_is_root(mpi_world)) then
      iunit = io_open('debug/mesh_partition/mesh_partition_boundary', &
        namespace, action='write')
      do ipg = mesh%np_global+1, mesh%np_part_global
        write(iunit, '(i8,99f18.8)') ipg, mesh_x_global(mesh, ipg)
      end do
      call io_close(iunit)
    end if

    POP_SUB(mesh_partition_messages_debug)

  end subroutine mesh_partition_messages_debug

end module mesh_partition_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
