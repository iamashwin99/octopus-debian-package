!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

!> @brief This module defines the meshes, which are used in Octopus
!!
module mesh_oct_m
  use basis_set_abst_oct_m
  use box_oct_m
  use comm_oct_m
  use coordinate_system_oct_m
  use debug_oct_m
  use global_oct_m
  use iihash_oct_m
  use index_oct_m
  use io_oct_m
  use io_binary_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use partition_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use symmetries_oct_m
  use symm_op_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                        &
    mesh_t,                        &
    mesh_plane_t,                  &
    mesh_line_t,                   &
    mesh_check_dump_compatibility, &
    mesh_end,                      &
    mesh_double_box,               &
    mesh_r,                        &
    mesh_gcutoff,                  &
    mesh_write_info,               &
    mesh_nearest_point,            &
    mesh_periodic_point,           &
    mesh_global_memory,            &
    mesh_local_memory,             &
    mesh_x_global,                 &
    mesh_write_fingerprint,        &
    mesh_read_fingerprint,         &
    mesh_compact_boundaries,       &
    mesh_check_symmetries,         &
    mesh_global_index_to_coords,   &
    mesh_global_index_from_coords, &
    mesh_local_index_to_coords,    &
    mesh_local_index_from_coords,  &
    mesh_local2global,             &
    mesh_global2local

  !> @brief Describes mesh distribution to nodes.
  !!
  !! Some general things:
  !! All members of type(mesh_t) are equal on all
  !! nodes when running parallel except
  !! - np, np_part
  !! - x, vol_pp
  !! These four are defined for all the points the node is responsible for.
  type, extends(basis_set_abst_t) :: mesh_t
    ! Components are public by default
    class(box_t),   pointer :: box  !< simulation box of box_t
    class(coordinate_system_t), pointer :: coord_system
    type(index_t)                :: idx
    logical :: use_curvilinear

    FLOAT :: spacing(MAX_DIM)         !< the (constant) spacing between the points

    ! When running serially, the local number of points is
    ! equal to the global number of points.
    ! Otherwise, the next two are different on each node.
    integer      :: np               !< Local number of points in mesh
    integer      :: np_part          !< Local points plus ghost points plus boundary points.
    integer(i8)  :: np_global        !< Global number of points in mesh.
    integer(i8)  :: np_part_global   !< Global number of inner points and boundary points.
    logical         :: parallel_in_domains !< will I run parallel in domains?
    type(mpi_grp_t) :: mpi_grp             !< the mpi group describing parallelization in domains
    type(par_vec_t) :: pv                  !< describes parallel vectors defined on the mesh.
    type(partition_t) :: partition         !< describes how the inner points are assigned to the domains

    FLOAT,   allocatable :: x(:,:)            !< The (local) \b points
    FLOAT                :: volume_element    !< The global volume element.
    FLOAT,   allocatable :: vol_pp(:)         !< Element of volume for curvilinear coordinates.

    logical :: masked_periodic_boundaries
    character(len=256) :: periodic_boundary_mask
  contains
    procedure :: end => mesh_end
    procedure :: init => mesh_init
    procedure :: write_info => mesh_write_info
    procedure :: dmesh_allreduce_0, zmesh_allreduce_0, imesh_allreduce_0
    procedure :: dmesh_allreduce_1, zmesh_allreduce_1, imesh_allreduce_1
    procedure :: dmesh_allreduce_2, zmesh_allreduce_2, imesh_allreduce_2
    procedure :: dmesh_allreduce_3, zmesh_allreduce_3, imesh_allreduce_3
    procedure :: dmesh_allreduce_4, zmesh_allreduce_4, imesh_allreduce_4
    procedure :: dmesh_allreduce_5, zmesh_allreduce_5, imesh_allreduce_5
    generic :: allreduce => dmesh_allreduce_0, zmesh_allreduce_0, imesh_allreduce_0
    generic :: allreduce => dmesh_allreduce_1, zmesh_allreduce_1, imesh_allreduce_1
    generic :: allreduce => dmesh_allreduce_2, zmesh_allreduce_2, imesh_allreduce_2
    generic :: allreduce => dmesh_allreduce_3, zmesh_allreduce_3, imesh_allreduce_3
    generic :: allreduce => dmesh_allreduce_4, zmesh_allreduce_4, imesh_allreduce_4
    generic :: allreduce => dmesh_allreduce_5, zmesh_allreduce_5, imesh_allreduce_5
  end type mesh_t

  !> @brief define a grid on a plane.
  !!
  !! This data type defines a plane, and a regular grid defined on
  !! this plane (or, rather, on a portion of this plane)
  !! n should be a unit vector, that determines the normal of the plane.
  !! Origin is a point belonging to the plane
  !! u and v are unit orthogonal vectors belonging to the plane
  !! The grid is generated by the vectors u and v:
  !!   x_{i,j} = origin + i*spacing*u + j*spacing*v,
  !! for nu <= i <= mu and nv <= j <= mv
  type mesh_plane_t
    ! Components are public by default
    FLOAT :: n(MAX_DIM)             !< normal vector (unit vector)
    FLOAT :: u(MAX_DIM), v(MAX_DIM) !< unit vectors, spanning the plane
    FLOAT :: origin(MAX_DIM)        !< origin of the plane
    FLOAT :: spacing
    integer :: nu, mu, nv, mv
  end type mesh_plane_t

  !> This data type defines a line, and a regular grid defined on this
  !! line (or rather, on a portion of this line).
  type mesh_line_t
    ! Components are public by default
    FLOAT :: n(MAX_DIM)          !< normal vector (unit vector)
    FLOAT :: u(MAX_DIM)          !< direction of the line
    FLOAT :: origin(MAX_DIM)     !< origin of the line
    FLOAT :: spacing
    integer :: nu, mu
  end type mesh_line_t

contains

  subroutine mesh_init(this)
    class(mesh_t), intent(inout) :: this

    PUSH_SUB(mesh_init)

    call this%set_time_dependent(.false.)

    POP_SUB(mesh_init)
  end subroutine mesh_init

! ---------------------------------------------------------
  !> finds the dimension of a box doubled in the non-periodic dimensions
  subroutine mesh_double_box(space, mesh, alpha, db)
    type(space_t),     intent(in)  :: space
    type(mesh_t),      intent(in)  :: mesh
    FLOAT,             intent(in)  :: alpha !< enlargement factor for double box
    integer,           intent(out) :: db(:)

    integer :: idir

    PUSH_SUB(mesh_double_box)

    db = 1

    ! double mesh with 2n points
    do idir = 1, space%periodic_dim
      db(idir) = mesh%idx%ll(idir)
    end do
    do idir = space%periodic_dim + 1, space%dim
      db(idir) = nint(alpha * (mesh%idx%ll(idir) - 1)) + 1
    end do

    POP_SUB(mesh_double_box)
  end subroutine mesh_double_box


  ! ---------------------------------------------------------
  subroutine mesh_write_info(this, iunit, namespace)
    class(mesh_t),               intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer :: ii
    FLOAT :: cutoff

    PUSH_SUB(mesh_write_info)

    write(message(1),'(3a)') '  Spacing [', trim(units_abbrev(units_out%length)), '] = ('
    do ii = 1, this%box%dim
      if (ii > 1) write(message(1), '(2a)') trim(message(1)), ','
      write(message(1), '(a,f6.3)') trim(message(1)), units_from_atomic(units_out%length, this%spacing(ii))
    end do
    write(message(1), '(5a,f12.5)') trim(message(1)), ') ', &
      '   volume/point [', trim(units_abbrev(units_out%length**this%box%dim)), '] = ',      &
      units_from_atomic(units_out%length**this%box%dim, this%vol_pp(1))

    write(message(2),'(a, i10)') '  # inner mesh = ', this%np_global
    write(message(3),'(a, i10)') '  # total mesh = ', this%np_part_global

    cutoff = mesh_gcutoff(this)**2 / M_TWO
    write(message(4),'(3a,f12.6,a,f12.6)') '  Grid Cutoff [', trim(units_abbrev(units_out%energy)),'] = ', &
      units_from_atomic(units_out%energy, cutoff), '    Grid Cutoff [Ry] = ', cutoff * M_TWO
    call messages_info(4, iunit=iunit, namespace=namespace)

    POP_SUB(mesh_write_info)
  end subroutine mesh_write_info


  ! ---------------------------------------------------------
  subroutine mesh_r(mesh, ip, rr, origin, coords)
    class(mesh_t), intent(in)  :: mesh
    integer,       intent(in)  :: ip
    FLOAT,         intent(out) :: rr
    FLOAT,         intent(in),  optional :: origin(:) !< origin(sb%dim)
    FLOAT,         intent(out), optional :: coords(:) !< coords(sb%dim)

    FLOAT :: xx(1:mesh%box%dim)

    ! no push_sub because it is called too frequently

    xx(1:mesh%box%dim) = mesh%x(ip, 1:mesh%box%dim)
    if (present(origin)) xx(1:mesh%box%dim) = xx(1:mesh%box%dim) - origin(1:mesh%box%dim)
    rr = norm2(xx(1:mesh%box%dim))

    if (present(coords)) then
      coords(1:mesh%box%dim) = xx(1:mesh%box%dim)
    end if

  end subroutine mesh_r

  !---------------------------------------------------------------------
  !> Returns the index of the point which is nearest to a given vector
  !! position pos. Variable dmin will hold, on exit, the distance between
  !! pos and this nearest mesh point. rankmin will be zero, if the mesh is
  !! not partitioned, and the rank of the processor which holds the point
  !! ind if the mesh is partitioned.
  ! ----------------------------------------------------------------------
  integer function mesh_nearest_point(mesh, pos, dmin, rankmin) result(ind)
    class(mesh_t),intent(in)  :: mesh
    FLOAT,        intent(in)  :: pos(:)
    FLOAT,        intent(out) :: dmin
    integer,      intent(out) :: rankmin

    FLOAT :: dd
    integer :: imin, ip
    FLOAT :: min_loc_in(2), min_loc_out(2)

    PUSH_SUB(mesh_nearest_point)

    ! find the point of the grid that is closer to the atom
    dmin = M_ZERO
    do ip = 1, mesh%np
      dd = sum((pos(1:mesh%box%dim) - mesh%x(ip, 1:mesh%box%dim))**2)
      if ((dd < dmin) .or. (ip == 1)) then
        imin = ip
        dmin = dd
      end if
    end do

    rankmin = 0
    if (mesh%parallel_in_domains) then
      min_loc_in(1) = dmin
      min_loc_in(2) = mesh%mpi_grp%rank
      call mesh%mpi_grp%allreduce(min_loc_in, min_loc_out, 1, MPI_2FLOAT, MPI_MINLOC)
      dmin = min_loc_out(1)
      rankmin = nint(min_loc_out(2))
      call mesh%mpi_grp%bcast(imin, 1, MPI_INTEGER, rankmin)
    end if

    ind = imin
    POP_SUB(mesh_nearest_point)
  end function mesh_nearest_point


  ! --------------------------------------------------------------
  !> mesh_gcutoff returns the "natural" band limitation of the
  !! grid mesh, in terms of the maximum G vector. For a cubic regular
  !! grid, it is M_PI/spacing.
  ! --------------------------------------------------------------
  FLOAT function mesh_gcutoff(mesh) result(gmax)
    class(mesh_t), intent(in) :: mesh

    PUSH_SUB(mesh_gcutoff)
    gmax = M_PI / (maxval(mesh%spacing))

    POP_SUB(mesh_gcutoff)
  end function mesh_gcutoff

  ! --------------------------------------------------------------
  subroutine mesh_write_fingerprint(mesh, dir, filename, mpi_grp, namespace, ierr)
    type(mesh_t),     intent(in)  :: mesh
    character(len=*), intent(in)  :: dir
    character(len=*), intent(in)  :: filename
    type(mpi_grp_t),  intent(in)  :: mpi_grp
    type(namespace_t),intent(in)  :: namespace
    integer,          intent(out) :: ierr

    integer :: iunit, ii

    PUSH_SUB(mesh_write_fingerprint)

    ierr = 0

    iunit = io_open(trim(dir)//"/"//trim(filename), namespace, action='write', &
      die=.false., grp=mpi_grp)
    if (iunit <= 0) then
      message(1) = "Unable to open file '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1, namespace=namespace)
      ierr = ierr + 1
    else
      if (mpi_grp_is_root(mpi_grp)) then
        write(iunit, '(a20,i21)')  'np_part_global=     ', mesh%np_part_global
        write(iunit, '(a20,i21)')  'np_global=          ', mesh%np_global
        write(iunit, '(a20,i21)')  'algorithm=          ', 1
        write(iunit, '(a20,i21)')  'checksum=           ', mesh%idx%checksum
        write(iunit, '(a20,i21)')  'bits=               ', mesh%idx%bits
        write(iunit, '(a20,i21)')  'dim=                ', mesh%idx%dim
        write(iunit, '(a20,i21)')  'type=               ', mesh%idx%type
        do ii = 1, mesh%idx%dim
          write(iunit, '(a7,i2,a11,i21)') 'offset(',ii,')=         ', mesh%idx%offset(ii)
        end do
        do ii = 1, mesh%idx%dim
          write(iunit, '(a7,i2,a11,i21)') 'nn(',ii,')=             ', mesh%idx%nr(2, ii) - mesh%idx%nr(1, ii) + 1
        end do
      end if
      call io_close(iunit, grp=mpi_grp)
    end if

    POP_SUB(mesh_write_fingerprint)
  end subroutine mesh_write_fingerprint


  ! -----------------------------------------------------------------------
  !> This function reads the fingerprint of a mesh written in
  !! filename. If the meshes are equal (same fingerprint) return values
  !! are 0, otherwise it returns the size of the mesh stored.
  !! fingerprint cannot be read, it returns ierr /= 0.
  subroutine mesh_read_fingerprint(mesh, dir, filename, mpi_grp, namespace, &
    read_np_part, read_np, bits, type, offset, nn, ierr)
    type(mesh_t),     intent(in)  :: mesh
    character(len=*), intent(in)  :: dir
    character(len=*), intent(in)  :: filename
    type(mpi_grp_t),  intent(in)  :: mpi_grp
    type(namespace_t),intent(in)  :: namespace
    integer(i8),      intent(out) :: read_np_part
    integer(i8),      intent(out) :: read_np
    integer,          intent(out) :: bits
    integer,          intent(out) :: type
    integer,          intent(out) :: offset(1:mesh%idx%dim)
    integer,          intent(out) :: nn(1:mesh%idx%dim)
    integer,          intent(out) :: ierr

    character(len=20)  :: str
    character(len=100) :: lines(7)
    integer :: iunit, algorithm, dim, err, ii
    integer(i8) :: checksum

    PUSH_SUB(mesh_read_fingerprint)

    ierr = 0

    read_np_part = 0_i8
    read_np = 0_i8

    iunit = io_open(trim(dir)//"/"//trim(filename), namespace, action='read', &
      status='old', die=.false., grp=mpi_grp)
    if (iunit <= 0) then
      ierr = ierr + 1
      message(1) = "Unable to open file '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1, namespace=namespace)
    else
      call iopar_read(mpi_grp, iunit, lines, 7, err)
      if (err /= 0) then
        ierr = ierr + 4
      else
        read(lines(1), '(a20,i21)')  str, read_np_part
        read(lines(2), '(a20,i21)')  str, read_np
        read(lines(3), '(a20,i21)')  str, algorithm
        read(lines(4), '(a20,i21)')  str, checksum
        read(lines(5), '(a20,i21)')  str, bits
        read(lines(6), '(a20,i21)')  str, dim
        read(lines(7), '(a20,i21)')  str, type
        ! only allow restarting simulations with the same dimensions
        if (dim /= mesh%idx%dim) then
          ierr = ierr + 8
        else
          ! read offset, has dim lines
          call iopar_read(mpi_grp, iunit, lines, dim, err)
          if (err /= 0) then
            ierr = ierr + 4
          else
            do ii = 1, dim
              read(lines(ii), '(a20,i21)')  str, offset(ii)
            end do
          end if

          ! read nn, has dim lines
          call iopar_read(mpi_grp, iunit, lines, dim, err)
          if (err /= 0) then
            ierr = ierr + 4
          else
            do ii = 1, dim
              read(lines(ii), '(a20,i21)')  str, nn(ii)
            end do
          end if
        end if

        ASSERT(read_np_part >= read_np)

        if (read_np_part == mesh%np_part_global &
          .and. read_np == mesh%np_global &
          .and. algorithm == 1 &
          .and. checksum == mesh%idx%checksum) then
          read_np_part = 0
          read_np = 0
        end if
      end if

      call io_close(iunit, grp=mpi_grp)
    end if

    POP_SUB(mesh_read_fingerprint)
  end subroutine mesh_read_fingerprint

  ! ---------------------------------------------------------
  subroutine mesh_check_dump_compatibility(mesh, dir, filename, namespace, mpi_grp, grid_changed, grid_reordered, map, ierr)
    type(mesh_t),             intent(in)  :: mesh
    character(len=*),         intent(in)  :: dir
    character(len=*),         intent(in)  :: filename
    type(namespace_t),        intent(in)  :: namespace
    type(mpi_grp_t),          intent(in)  :: mpi_grp
    logical,                  intent(out) :: grid_changed
    logical,                  intent(out) :: grid_reordered
    integer(i8), allocatable, intent(out) :: map(:)
    integer,                  intent(out) :: ierr

    integer(i8) :: ipg, ipg_new, read_np_part, read_np
    integer :: err, idir
    integer :: bits, type, offset(mesh%idx%dim), point(mesh%idx%dim), nn(mesh%idx%dim)
    integer(i8), allocatable :: read_indices(:)
    type(index_t) :: idx_old

    PUSH_SUB(mesh_check_dump_compatibility)

    ierr = 0

    grid_changed = .false.
    grid_reordered = .false.

    ! Read the mesh fingerprint
    call mesh_read_fingerprint(mesh, dir, filename, mpi_grp, namespace, read_np_part, read_np, &
      bits, type, offset, nn, err)
    if (err /= 0) then
      ierr = ierr + 1
      message(1) = "Unable to read mesh fingerprint from '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1, namespace=namespace)

    else if (read_np > 0) then
      if (.not. associated(mesh%box)) then
        ! We can only check the compatibility of two meshes that have different fingerprints if we also
        ! have the simulation box. In the case we do not, we will assume that the fingerprint is enough.
        ierr = ierr + 2
      else
        grid_changed = .true.

        ! perhaps only the order of the points changed, this can only
        ! happen if the number of points is the same and no points maps
        ! to zero (this is checked below)
        grid_reordered = (read_np == mesh%np_global)

        ! the grid is different, so we read the coordinates.
        SAFE_ALLOCATE(read_indices(1:read_np_part))
        call io_binary_read(trim(io_workpath(dir, namespace))//"/indices.obf", read_np_part, &
          read_indices, err)
        if (err /= 0) then
          ierr = ierr + 4
          message(1) = "Unable to read index map from '"//trim(dir)//"'."
          call messages_warning(1, namespace=namespace)
        else
          ! dummy index object
          idx_old%type = type
          idx_old%bits = bits
          idx_old%nr(1, 1:mesh%idx%dim) = -offset(1:mesh%idx%dim)
          idx_old%nr(2, 1:mesh%idx%dim) = -offset(1:mesh%idx%dim) + nn(1:mesh%idx%dim) - 1
          idx_old%offset(1:mesh%idx%dim) = offset(1:mesh%idx%dim)
          idx_old%stride(1) = 1
          do idir = 2, mesh%idx%dim
            idx_old%stride(idir) = idx_old%stride(idir-1) * nn(idir-1)
          end do
          ! generate the map
          SAFE_ALLOCATE(map(1:read_np))
          do ipg = 1, read_np
            ! get nd-index from old 1d index
            call index_spatial_to_point(idx_old, mesh%idx%dim, read_indices(ipg), point)
            ! get new global index
            ipg_new = mesh_global_index_from_coords(mesh, point)
            map(ipg) = ipg_new
            ! ignore boundary points
            if (map(ipg) > mesh%np_global) map(ipg) = 0
            ! if the map is zero for one point, it is not a simple reordering
            if (map(ipg) == 0) grid_reordered = .false.
          end do
        end if

        SAFE_DEALLOCATE_A(read_indices)
      end if
    end if

    POP_SUB(mesh_check_dump_compatibility)
  end subroutine mesh_check_dump_compatibility


  ! --------------------------------------------------------------
  recursive subroutine mesh_end(this)
    class(mesh_t), intent(inout)   :: this

    PUSH_SUB(mesh_end)

#ifdef HAVE_MPI
    call lmpi_destroy_shared_memory_window(this%idx%window_grid_to_spatial)
    call lmpi_destroy_shared_memory_window(this%idx%window_spatial_to_grid)
    nullify(this%idx%grid_to_spatial_global)
    nullify(this%idx%spatial_to_grid_global)
#else
    SAFE_DEALLOCATE_P(this%idx%grid_to_spatial_global)
    SAFE_DEALLOCATE_P(this%idx%spatial_to_grid_global)
#endif

    SAFE_DEALLOCATE_A(this%x)
    SAFE_DEALLOCATE_A(this%vol_pp)

    if (this%parallel_in_domains) then
      call par_vec_end(this%pv)
      call partition_end(this%partition)
    end if

    POP_SUB(mesh_end)
  end subroutine mesh_end


  !> This function returns the point inside the grid corresponding to
  !! a boundary point when PBCs are used. In case the point does not
  !! have a correspondence (i.e. other BCs are used in that direction),
  !! the same point is returned. Note that this function returns a
  !! global point number when parallelization in domains is used.
  ! ---------------------------------------------------------
  integer(i8) function mesh_periodic_point(mesh, space, ip) result(ipg)
    class(mesh_t), intent(in)    :: mesh
    type(space_t), intent(in)    :: space
    integer,       intent(in)    :: ip     !< local point for which periodic copy is searched

    integer :: ix(MAX_DIM), nr(2, MAX_DIM), idim
    FLOAT :: xx(MAX_DIM), rr, ufn_re, ufn_im

    ! no push_sub, called too frequently

    call mesh_local_index_to_coords(mesh, ip, ix)
    nr(1, :) = mesh%idx%nr(1, :) + mesh%idx%enlarge(:)
    nr(2, :) = mesh%idx%nr(2, :) - mesh%idx%enlarge(:)

    do idim = 1, space%periodic_dim
      if (ix(idim) < nr(1, idim)) ix(idim) = ix(idim) + mesh%idx%ll(idim)
      if (ix(idim) > nr(2, idim)) ix(idim) = ix(idim) - mesh%idx%ll(idim)
    end do

    ipg = mesh_global_index_from_coords(mesh, ix)
    ASSERT(ipg > 0)

    if (mesh%masked_periodic_boundaries) then
      call mesh_r(mesh, ip, rr, coords = xx)
      call parse_expression(ufn_re, ufn_im, space%dim, xx, rr, M_ZERO, mesh%periodic_boundary_mask)
      if (int(ufn_re) == 0) ipg = mesh_local2global(mesh, ip) ! Nothing will be done
    end if

  end function mesh_periodic_point


  ! ---------------------------------------------------------
  FLOAT pure function mesh_global_memory(mesh) result(memory)
    class(mesh_t), intent(in) :: mesh

    ! 2 global index arrays
    memory = SIZEOF_UNSIGNED_LONG_LONG * TOFLOAT(mesh%np_part_global) * 2

  end function mesh_global_memory


  ! ---------------------------------------------------------
  FLOAT pure function mesh_local_memory(mesh) result(memory)
    class(mesh_t), intent(in) :: mesh

    memory = M_ZERO

    ! x
    memory = memory + REAL_PRECISION * TOFLOAT(mesh%np_part) * MAX_DIM
    ! local index arrays
    memory = memory + SIZEOF_UNSIGNED_LONG_LONG * TOFLOAT(mesh%np_part) * 2
  end function mesh_local_memory


  ! ---------------------------------------------------------
  function mesh_x_global(mesh, ipg) result(xx)
    class(mesh_t),      intent(in) :: mesh
    integer(i8),        intent(in) :: ipg
    FLOAT                          :: xx(1:mesh%box%dim)

    FLOAT :: chi(1:mesh%box%dim)
    integer :: ix(1:mesh%box%dim)

! no push_sub because function is called too frequently

    call mesh_global_index_to_coords(mesh, ipg, ix)
    chi(1:mesh%box%dim) = ix(1:mesh%box%dim) * mesh%spacing(1:mesh%box%dim)
    xx = mesh%coord_system%to_cartesian(chi(1:mesh%box%dim))

  end function mesh_x_global


  ! ---------------------------------------------------------
  logical pure function mesh_compact_boundaries(mesh) result(cb)
    type(mesh_t),       intent(in) :: mesh

    cb = .not. mesh%use_curvilinear .and. &
      .not. mesh%parallel_in_domains

  end function mesh_compact_boundaries


  ! ---------------------------------------------------------
  subroutine mesh_check_symmetries(mesh, symm, periodic_dim)
    class(mesh_t),       intent(in) :: mesh
    type(symmetries_t),  intent(in) :: symm
    integer,             intent(in) :: periodic_dim

    integer :: iop, ip, idim, nops, ix(1:3)
    FLOAT :: destpoint(1:3), srcpoint(1:3), lsize(1:3), offset(1:3)

    !If all the axis have the same spacing and the same length
    !the grid is by obviously symmetric
    !Indeed, reduced coordinates are proportional to the point index
    !and the reduced rotation are integer matrices
    !The result of the product is also proportional to an integer
    !and therefore belong to the grid.
    if (mesh%idx%ll(1) == mesh%idx%ll(2)  .and. &
      mesh%idx%ll(2)   == mesh%idx%ll(3)  .and. &
      mesh%spacing(1)  == mesh%spacing(2) .and. &
      mesh%spacing(2)  == mesh%spacing(3)) return

    PUSH_SUB(mesh_check_symmetries)

    message(1) = "Checking if the real-space grid is symmetric"
    call messages_info(1)

    lsize(1:3) = TOFLOAT(mesh%idx%ll(1:3))
    offset(1:3) = TOFLOAT(mesh%idx%nr(1, 1:3) + mesh%idx%enlarge(1:3))

    nops = symmetries_number(symm)

    do ip = 1, mesh%np
      !We use floating point coordinates to check if the symmetric point
      !belong to the grid.
      !If yes, it should have integer reduced coordinates
      call mesh_local_index_to_coords(mesh, ip, ix)
      destpoint(1:3) = TOFLOAT(ix(1:3)) - offset(1:3)
      ! offset moves corner of cell to origin, in integer mesh coordinates
      ASSERT(all(destpoint >= 0))
      ASSERT(all(destpoint < lsize))

      ! move to center of cell in real coordinates
      destpoint = destpoint - TOFLOAT(int(lsize)/2)

      !convert to proper reduced coordinates
      do idim = 1, 3
        destpoint(idim) = destpoint(idim)/lsize(idim)
      end do

      ! iterate over all points that go to this point by a symmetry operation
      do iop = 1, nops
        srcpoint = symm_op_apply_red(symm%ops(iop), destpoint)

        !We now come back to what should be an integer, if the symmetric point beloings to the grid
        do idim = 1, 3
          srcpoint(idim) = srcpoint(idim)*lsize(idim)
        end do

        ! move back to reference to origin at corner of cell
        srcpoint = srcpoint + TOFLOAT(int(lsize)/2)

        ! apply periodic boundary conditions in periodic directions
        do idim = 1, periodic_dim
          if (nint(srcpoint(idim)) < 0 .or. nint(srcpoint(idim)) >= lsize(idim)) then
            srcpoint(idim) = modulo(srcpoint(idim)+M_HALF*symprec, lsize(idim))
          end if
        end do
        ASSERT(all(srcpoint >= -symprec))
        ASSERT(all(srcpoint < lsize))

        srcpoint(1:3) = srcpoint(1:3) + offset(1:3)

        if (any(srcpoint-anint(srcpoint)> symprec*M_TWO)) then
          message(1) = "The real-space grid breaks at least one of the symmetries of the system."
          message(2) = "Change your spacing or use SymmetrizeDensity=no."
          call messages_fatal(2)
        end if
      end do
    end do

    POP_SUB(mesh_check_symmetries)
  end subroutine

  !> This function returns the true _global_ index of the point for a given
  !! vector of integer coordinates.
  integer(i8) function mesh_global_index_from_coords(mesh, ix) result(index)
    class(mesh_t), intent(in)    :: mesh
    integer,       intent(in)    :: ix(:)

    index = index_from_coords(mesh%idx, ix)
  end function mesh_global_index_from_coords

  !> Given a _global_ point index, this function returns the set of
  !! integer coordinates of the point.
  subroutine mesh_global_index_to_coords(mesh, ipg, ix)
    type(mesh_t),  intent(in)    :: mesh
    integer(i8),   intent(in)    :: ipg
    integer,       intent(out)   :: ix(:)

    call index_to_coords(mesh%idx, ipg, ix)
  end subroutine mesh_global_index_to_coords

  !> This function returns the _local_ index of the point for a given
  !! vector of integer coordinates.
  integer function mesh_local_index_from_coords(mesh, ix) result(ip)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ix(:)

    integer(i8) :: ipg

    ipg = index_from_coords(mesh%idx, ix)
    ip = mesh_global2local(mesh, ipg)
  end function mesh_local_index_from_coords

  !> Given a _local_ point index, this function returns the set of
  !! integer coordinates of the point.
  subroutine mesh_local_index_to_coords(mesh, ip, ix)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ip
    integer,       intent(out)   :: ix(:)

    integer(i8) :: ipg

    ipg = mesh_local2global(mesh, ip)
    call index_to_coords(mesh%idx, ipg, ix)
  end subroutine mesh_local_index_to_coords

  !> This function returns the global mesh index for a given local index
  integer(i8) function mesh_local2global(mesh, ip) result(ipg)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ip

    ipg = par_vec_local2global(mesh%pv, ip)
  end function mesh_local2global

  !> This function returns the local mesh index for a given global index
  integer function mesh_global2local(mesh, ipg) result(ip)
    type(mesh_t),  intent(in)    :: mesh
    integer(i8),   intent(in)    :: ipg

    ip = par_vec_global2local(mesh%pv, ipg)
  end function mesh_global2local

#include "undef.F90"
#include "real.F90"
#include "mesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "mesh_inc.F90"

end module mesh_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
