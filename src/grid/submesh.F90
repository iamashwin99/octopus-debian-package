!! Copyright (C) 2007 X. Andrade
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

module submesh_oct_m
  use accel_oct_m
  use batch_oct_m
  use boundaries_oct_m
  use box_oct_m
  use comm_oct_m
  use debug_oct_m
  use global_oct_m
  use index_oct_m
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use messages_oct_m
  use sort_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use mpi_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use space_oct_m
  use types_oct_m

  implicit none
  private

  public ::                      &
    submesh_t,                   &
    submesh_init,                &
    submesh_init_box,            &
    submesh_merge,               &
    submesh_shift_center,        &
    submesh_broadcast,           &
    submesh_get_inv,             &
    submesh_build_global,        &
    submesh_end_global,          &
    dsm_integrate,               &
    zsm_integrate,               &
    dsm_integrate_frommesh,      &
    zsm_integrate_frommesh,      &
    dsm_nrm2,                    &
    zsm_nrm2,                    &
    submesh_add_to_mesh,         &
    dsubmesh_batch_add,          &
    zsubmesh_batch_add,          &
    dzsubmesh_batch_add,         &
    submesh_to_mesh_dotp,        &
    dsubmesh_batch_add_matrix,   &
    zsubmesh_batch_add_matrix,   &
    dsubmesh_batch_dotp_matrix,  &
    zsubmesh_batch_dotp_matrix,  &
    submesh_overlap,             &
    dsubmesh_copy_from_mesh,     &
    zsubmesh_copy_from_mesh,     &
    dsubmesh_copy_from_mesh_batch,     &
    zsubmesh_copy_from_mesh_batch,     &
    submesh_compatible,          &
    submesh_end,                 &
    submesh_get_cube_dim,        &
    submesh_init_cube_map,       &
    submesh_end_cube_map

  !> A submesh is a type of mesh, used for the projectors in the pseudopotentials
  !! It contains points on a regular mesh confined to a sphere of a given radius.
  type submesh_t
    ! Components are public by default
    FLOAT, allocatable    :: center(:)
    FLOAT                 :: radius = M_ZERO
    class(box_t), pointer :: box => NULL()  !< simulation box of box_t
    integer               :: np = -1        !< number of points inside the submesh
    integer,  allocatable :: map(:)         !< maps point inside the submesh to a point inside the underlying mesh
    integer               :: num_regions    !< number of injective regions of the map
    integer,  allocatable :: regions(:)     !< offsets for regions of the map
    type(accel_mem_t)     :: buff_map
    FLOAT,    allocatable :: rel_x(:,:)     !< Relative position to the center (1:space%dim, 1:np_part)
    FLOAT,    allocatable :: r(:)           !< distance from centre of the submesh.
    type(mesh_t), pointer :: mesh => NULL() !< pointer to the underlying mesh
    logical               :: overlap        !< .true. if the submesh has more than one point that is mapped to a mesh point,
    !!                                         i.e. the submesh overlaps with itself (as can happen in periodic systems)
    integer               :: np_global = -1 !< total number of points in the entire mesh
    FLOAT,    allocatable :: rel_x_global(:,:)
    integer,  allocatable :: part_v(:)
    integer,  allocatable :: global2local(:)

    type(mesh_cube_map_t) :: cube_map
  end type submesh_t

  interface submesh_add_to_mesh
    module procedure ddsubmesh_add_to_mesh, zdsubmesh_add_to_mesh, zzsubmesh_add_to_mesh
  end interface submesh_add_to_mesh

  interface submesh_to_mesh_dotp
    module procedure ddsubmesh_to_mesh_dotp, zdsubmesh_to_mesh_dotp, zzsubmesh_to_mesh_dotp
  end interface submesh_to_mesh_dotp

contains

  ! -------------------------------------------------------------
  ! Multipliers for recursive formulation of n-ellipsoid volume
  ! simplifying the Gamma function
  ! f(n) = 2f(n-2)/n, f(0)=1, f(1)=2
  recursive FLOAT function f_n(dims) result(fn)
    integer :: dims

    if (dims == 0) then
      fn = M_ONE
    else if (dims == 1) then
      fn = M_TWO
    else
      fn = M_TWO * f_n(dims - 2) / dims
    end if

  end function f_n

  !> This subroutine creates a submesh which has the same box of the one that the user passes
  subroutine submesh_init_box(this, space, mesh, box, center)
    type(submesh_t),         intent(inout) :: this
    type(space_t),           intent(in)    :: space
    class(mesh_t), target,   intent(in)    :: mesh
    class(box_t), target,    intent(in)    :: box   !< Simulation box of the submesh. Must be contained in the box of the mesh!
    FLOAT,                   intent(in)    :: center(:)

    integer :: original_mesh_index, submesh_index
    logical :: submesh_points(1:mesh%np)
    FLOAT, allocatable :: xtmp(:, :), rtmp(:)

    PUSH_SUB(submesh_init_box)

    ASSERT(.not. space%is_periodic())

    ! Assign the mesh and the box
    this%mesh => mesh
    this%box => box

    ! Get the total number of points inside the submesh. First of all get the number of points inside the submesh box
    submesh_points = box%contains_points(mesh%np, mesh%x)
    this%np = count(submesh_points)
    ! Verify that the box of the submesh is contained into the box of the mesh
    ASSERT(this%np <= mesh%np)

    SAFE_ALLOCATE(this%map(1:this%np))
    SAFE_ALLOCATE(xtmp(1:space%dim, 1:this%np))
    SAFE_ALLOCATE(rtmp(1:this%np))
    xtmp = M_ZERO
    rtmp = M_ZERO

    ! Now we generate a map between the points of the submesh and the points of the original mesh
    submesh_index = 0
    do original_mesh_index = 1, mesh%np
      ! check that the point is contained in the new box
      if (submesh_points(original_mesh_index)) then
        submesh_index = submesh_index + 1
        this%map(submesh_index) = original_mesh_index
        ! Assign the coordinates to this new point
        xtmp(:, submesh_index) = mesh%x(original_mesh_index, :) - center(:)
        rtmp(submesh_index) = norm2(xtmp(:, submesh_index))
      end if
    end do

    call submesh_reorder_points(this, space, xtmp, rtmp)

    SAFE_DEALLOCATE_A(xtmp)
    SAFE_DEALLOCATE_A(rtmp)

    POP_SUB(submesh_init_box)
  end subroutine submesh_init_box

  ! -------------------------------------------------------------
  subroutine submesh_init(this, space, mesh, latt, center, rc)
    type(submesh_t),         intent(inout)  :: this
    type(space_t),           intent(in)     :: space
    class(mesh_t), target,   intent(in)     :: mesh
    type(lattice_vectors_t), intent(in)     :: latt
    FLOAT,                   intent(in)     :: center(1:space%dim)
    FLOAT,                   intent(in)     :: rc

    FLOAT :: r2, rc2, xx(space%dim), rc_norm_n
    FLOAT, allocatable :: center_copies(:,:), xtmp(:, :), rtmp(:), x(:,:)
    integer :: icell, is, ip, ix, iy, iz
    integer(i8) :: max_elements_count
    type(profile_t), save :: submesh_init_prof, prof1, prof2
    type(lattice_iterator_t) :: latt_iter
    integer, allocatable :: map_inv(:), map_temp(:)
    integer :: nmax(3), nmin(3)


    PUSH_SUB(submesh_init)
    call profiling_in(submesh_init_prof, "SUBMESH_INIT")

    ASSERT(space%dim <= 3)

    this%mesh => mesh

    SAFE_ALLOCATE(this%center(1:space%dim))
    this%center = center

    this%radius = rc
    rc2 = rc**2

    ! The spheres are generated differently for periodic coordinates,
    ! mainly for performance reasons.
    if (.not. space%is_periodic()) then

      call profiling_in(prof1, "SUBMESH_INIT_MAP_INV")
      SAFE_ALLOCATE(map_inv(0:this%mesh%np))
      map_inv(0:this%mesh%np) = 0

      nmin = 0
      nmax = 0

      ! get a cube of points that contains the sphere
      nmin(1:space%dim) = int((center(1:space%dim) - abs(rc))/mesh%spacing(1:space%dim)) - 1
      nmax(1:space%dim) = int((center(1:space%dim) + abs(rc))/mesh%spacing(1:space%dim)) + 1

      ! make sure that the cube is inside the grid
      ! parts of the cube which would fall outside the simulation box are chopped off.
      nmin(1:space%dim) = max(mesh%idx%nr(1, 1:space%dim), nmin(1:space%dim))
      nmax(1:space%dim) = min(mesh%idx%nr(2, 1:space%dim), nmax(1:space%dim))

      ! Get the total number of points inside the sphere
      is = 0   ! this index counts inner points
      do iz = nmin(3), nmax(3)
        do iy = nmin(2), nmax(2)
          do ix = nmin(1), nmax(1)
            ip = mesh_local_index_from_coords(mesh, [ix, iy, iz])
            if (ip == 0 .or. ip > mesh%np) cycle
            r2 = sum((mesh%x(ip, :) - center)**2)
            if (r2 <= rc2) then
              is = is + 1
              map_inv(ip) = is
            end if
          end do
        end do
      end do
      this%np = is
      call profiling_out(prof1)

      call profiling_in(prof2, "SUBMESH_INIT_RTMP")
      SAFE_ALLOCATE(this%map(1:this%np))
      SAFE_ALLOCATE(xtmp(1:space%dim, 1:this%np))
      SAFE_ALLOCATE(rtmp(1:this%np))

      ! Generate the table and the positions
      do iz = nmin(3), nmax(3)
        do iy = nmin(2), nmax(2)
          do ix = nmin(1), nmax(1)
            ip = mesh_local_index_from_coords(mesh, [ix, iy, iz])
            if (ip == 0 .or. ip > mesh%np) cycle
            is = map_inv(ip)
            if (is == 0) cycle
            this%map(is) = ip
            xtmp(:, is) = mesh%x(ip, :) - center
            rtmp(is) = norm2(xtmp(:,is))
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(map_inv)
      call profiling_out(prof2)

      ! This is the case for a periodic system
    else

      ! Get the total number of points inside the sphere considering
      ! replicas along PBCs

      ! this requires some optimization

      latt_iter = lattice_iterator_t(latt, rc)

      SAFE_ALLOCATE(center_copies(1:space%dim, 1:latt_iter%n_cells))
      do icell = 1, latt_iter%n_cells
        center_copies(:, icell) = center + latt_iter%get(icell)
      end do

      !Recursive formulation for the volume of n-ellipsoid
      !Garry Tee, NZ J. Mathematics Vol. 34 (2005) p. 165 eqs. 53,55
      rc_norm_n = product(ceiling(rc / mesh%spacing(1:space%dim), i8) + M_ONE)
      if (mesh%use_curvilinear) rc_norm_n = rc_norm_n / mesh%coord_system%min_mesh_scaling_product
      max_elements_count = 3**space%dim * int(M_PI**floor(0.5 * space%dim) * rc_norm_n * f_n(space%dim), i8)

      SAFE_ALLOCATE(x(1:space%dim,1:mesh%np))
      do ip = 1, mesh%np
        x(:,ip) = mesh%x(ip,:)
      end do
      call profiling_in(prof1, "SUBMESH_INIT_PERIODIC_R")
      SAFE_ALLOCATE(map_temp(1:max_elements_count))
      SAFE_ALLOCATE(xtmp(1:space%dim, 1:max_elements_count))
      SAFE_ALLOCATE(rtmp(1:max_elements_count))

      is = 0
      do ip = 1, mesh%np
        do icell = 1, latt_iter%n_cells
          xx = x(:,ip) - center_copies(:, icell)
          if(any(abs(xx)>rc)) cycle
          r2 = sum(xx**2)
          if (r2 > rc2) cycle
          is = is + 1
          map_temp(is) = ip
          rtmp(is) = sqrt(r2)
          xtmp(:, is) = xx
          ! Note that xx can be outside the unit cell
        end do
      end do
      ASSERT(is < huge(is))
      this%np = is

      SAFE_ALLOCATE(this%map(1:this%np))
      this%map(1:this%np) = map_temp(1:this%np)
      call profiling_out(prof1)

      SAFE_DEALLOCATE_A(map_temp)
      SAFE_DEALLOCATE_A(center_copies)
      SAFE_DEALLOCATE_A(x)

    end if

    call submesh_reorder_points(this, space, xtmp, rtmp)
    call profiling_out(submesh_init_prof)

    SAFE_DEALLOCATE_A(xtmp)
    SAFE_DEALLOCATE_A(rtmp)

    POP_SUB(submesh_init)
  end subroutine submesh_init

  subroutine submesh_reorder_points(this, space, xtmp, rtmp)
    type(submesh_t), intent(inout) :: this
    type(space_t),   intent(in)    :: space
    FLOAT,           intent(in)    :: xtmp(:, :), rtmp(:)

    integer :: ip, i_region, offset
    integer, allocatable :: order(:), order_new(:)
    integer, allocatable :: map_new(:)
    integer, allocatable :: np_region(:), tmp_array(:)
    type(profile_t), save :: profiler_1, profiler_2

    PUSH_SUB(submesh_reorder_points)

    ! now order points for better locality

    call profiling_in(profiler_1, "SUBMESH_INIT_ORDER")
    SAFE_ALLOCATE(order(1:this%np))
    SAFE_ALLOCATE(this%rel_x(1:space%dim, 1:this%np))
    SAFE_ALLOCATE(this%r(1:this%np))

    do ip = 1, this%np
      order(ip) = ip
    end do

    ! First we just reorder in order to determine overlap:
    call sort(this%map, order)

    !check whether points overlap (i.e. whetehr a submesh contains the same point more than once)
    this%overlap = .false.
    do ip = 1, this%np - 1
      if (this%map(ip) == this%map(ip + 1)) then
        ! this simplified test works, as the points are ordered.
        this%overlap = .true.
        exit
      end if
    end do

    this%num_regions = 1
    call profiling_out(profiler_1)

    if(this%overlap) then
      call profiling_in(profiler_2, "SUBMESH_INIT_OVERLAP")
      !disentangle the map into injective regions

      SAFE_ALLOCATE(tmp_array(1:this%np))
      SAFE_ALLOCATE(order_new(1:this%np))
      SAFE_ALLOCATE(np_region(1:this%np))
      SAFE_ALLOCATE(map_new(  1:this%np))

      np_region(1) = 1
      tmp_array(1) = 1
      i_region = 1

      do ip = 2, this%np
        if (this%map(ip) == this%map(ip - 1)) then
          i_region = i_region + 1
          if (i_region > this%num_regions) then
            this%num_regions = i_region
            np_region(i_region) = 0
          end if
        else
          i_region = 1
        end if
        tmp_array(ip) = i_region                        ! which region does ip belong to
        np_region(i_region) = np_region(i_region) + 1   ! increase number of points in i_region
      end do

      ASSERT( .not. allocated(this%regions))

      ! construct array of offsets
      SAFE_ALLOCATE(this%regions(1:this%num_regions+1))

      this%regions(1) = 1

      if(this%num_regions > 1) then
        do i_region = 1, this%num_regions
          this%regions(i_region + 1) = this%regions(i_region) + np_region(i_region)
        end do
      else
        this%regions(2) = this%np + 1
      end if

      np_region(1:this%np) = 0
      order_new(1:this%np) = -1
      map_new(1:this%np) = -1


      !reassemble regions into global map array
      do ip = 1, this%np
        i_region = tmp_array(ip)
        np_region(i_region) = np_region(i_region) + 1
        offset = this%regions(i_region) - 1
        map_new(   offset + np_region(i_region) ) = this%map(ip)
        order_new( offset + np_region(i_region) ) = order(ip)
      end do

      order(1:this%np) = order_new(1:this%np)
      this%map(1:this%np) = map_new(1:this%np)

      SAFE_DEALLOCATE_A(tmp_array)
      SAFE_DEALLOCATE_A(order_new)
      SAFE_DEALLOCATE_A(np_region)
      SAFE_DEALLOCATE_A(map_new)
      call profiling_out(profiler_2)

    else
      this%num_regions = 1
      SAFE_ALLOCATE(this%regions(1:2))
      this%regions(1) = 1
      this%regions(2) = this%np + 1
    end if

    ! Lastly, reorder the points according to the new scheme
    do ip = 1, this%np
      this%rel_x(:, ip) = xtmp(:, order(ip))
      this%r(ip) = rtmp(order(ip))
    end do

    SAFE_DEALLOCATE_A(order)

    POP_SUB(submesh_reorder_points)
  end subroutine submesh_reorder_points


  ! --------------------------------------------------------------
  !This routine takes two submeshes and merge them into a bigger submesh
  !The grid is centered on the first center
  subroutine submesh_merge(this, space, mesh, sm1, sm2, shift)
    type(submesh_t),       intent(inout)  :: this !< valgrind objects to intent(out) due to the initializations above
    type(space_t),         intent(in)     :: space
    class(mesh_t), target, intent(in)     :: mesh
    type(submesh_t),       intent(in)     :: sm1
    type(submesh_t),       intent(in)     :: sm2
    FLOAT, optional,       intent(in)     :: shift(:) !< If present, shifts the center of sm2

    FLOAT :: r2
    integer :: ip, is
    type(profile_t), save :: prof
    FLOAT :: xx(space%dim), diff_centers(space%dim)

    PUSH_SUB(submesh_merge)
    call profiling_in(prof, "SUBMESH_MERGE")

    this%mesh => mesh

    SAFE_ALLOCATE(this%center(1:space%dim))
    this%center  = sm1%center
    this%radius = sm1%radius

    ! This is a quick fix to prevent uninitialized variables. To properly check the self-overlap,
    ! a similar approach as in submesh_init should be taken with respect to the merged map.
    this%overlap = sm1%overlap .or. sm2%overlap

    diff_centers = sm1%center - sm2%center
    if (present(shift)) diff_centers = diff_centers - shift

    !As we take the union of the two submeshes, we know that we have all the points from the first one included.
    !The extra points from the second submesh are those which are not included in the first one
    is = sm1%np
    do ip = 1, sm2%np
      !sm2%x contains points coordinates defined with respect to sm2%center
      xx = sm2%rel_x(:, ip) - diff_centers
      !If the point is not in sm1, we add it
      if (sum(xx**2) > sm1%radius**2) is = is + 1
    end do

    this%np = is

    SAFE_ALLOCATE(this%map(1:this%np))
    SAFE_ALLOCATE(this%rel_x(1:space%dim, 1:this%np))
    SAFE_ALLOCATE(this%r(1:this%np))
    this%map(1:sm1%np) = sm1%map(1:sm1%np)
    this%rel_x(:, 1:sm1%np) = sm1%rel_x(:, 1:sm1%np)
    this%r(1:sm1%np) = sm1%r(1:sm1%np)

    !iterate again to fill the tables
    is = sm1%np
    do ip = 1, sm2%np
      xx = sm2%rel_x(:, ip) - diff_centers
      r2 = sum(xx**2)
      if (r2 > sm1%radius**2) then
        is = is + 1
        this%map(is) = sm2%map(ip)
        this%r(is) = sqrt(r2)
        this%rel_x(:, is) = xx
      end if
    end do

    call profiling_out(prof)
    POP_SUB(submesh_merge)
  end subroutine submesh_merge

  ! --------------------------------------------------------------
  !This routine shifts the center of a submesh, without changing the grid points
  subroutine submesh_shift_center(this, space, newcenter)
    type(submesh_t),      intent(inout)  :: this
    type(space_t),        intent(in)     :: space
    FLOAT,                intent(in)     :: newcenter(:)

    FLOAT :: xx(space%dim), diff_centers(space%dim), oldcenter(space%dim)
    integer :: ip
    type(profile_t), save :: prof

    PUSH_SUB(submesh_shift_center)
    call profiling_in(prof, "SUBMESH_SHIFT")

    oldcenter = this%center
    this%center  = newcenter

    diff_centers = newcenter - oldcenter

    do ip = 1, this%np
      xx = this%rel_x(:, ip) - diff_centers
      this%r(ip) = norm2(xx)
      this%rel_x(:, ip) = xx
    end do

    call profiling_out(prof)
    POP_SUB(submesh_shift_center)
  end subroutine submesh_shift_center

  ! --------------------------------------------------------------
  subroutine submesh_broadcast(this, space, mesh, center, radius, root, mpi_grp)
    type(submesh_t),      intent(inout)  :: this
    type(space_t),        intent(in)     :: space
    type(mesh_t), target, intent(in)     :: mesh
    FLOAT,                intent(in)     :: center(1:space%dim)
    FLOAT,                intent(in)     :: radius
    integer,              intent(in)     :: root
    type(mpi_grp_t),      intent(in)     :: mpi_grp

    integer :: nparray(1:3)
    type(profile_t), save :: prof

    PUSH_SUB(submesh_broadcast)
    call profiling_in(prof, 'SUBMESH_BCAST')

    if (root /= mpi_grp%rank) then
      this%mesh => mesh
      this%center = center
      this%radius = radius
    end if

    if (mpi_grp%size > 1) then

      if (root == mpi_grp%rank) then
        nparray(1) = this%np
        nparray(2) = this%num_regions
        if (this%overlap) then
          nparray(3) = 1
        else
          nparray(3) = 0
        end if
      end if

      call mpi_grp%bcast(nparray, 2, MPI_INTEGER, root)
      this%np = nparray(1)
      this%num_regions = nparray(2)
      this%overlap = (nparray(3) == 1)

      if (root /= mpi_grp%rank) then
        SAFE_ALLOCATE(this%map(1:this%np))
        SAFE_ALLOCATE(this%rel_x(1:space%dim, 1:this%np))
        SAFE_ALLOCATE(this%r(1:this%np))
        SAFE_ALLOCATE(this%regions(1:this%num_regions+1))
      end if

      call mpi_grp%bcast(this%regions(1), this%num_regions+1, MPI_INTEGER, root)

      if (this%np > 0) then
        call mpi_grp%bcast(this%map(1), this%np, MPI_INTEGER, root)
        call mpi_grp%bcast(this%rel_x(1, 1), this%np*space%dim, MPI_FLOAT, root)
        call mpi_grp%bcast(this%r(1), this%np, MPI_FLOAT, root)
      end if

    end if

    call profiling_out(prof)
    POP_SUB(submesh_broadcast)
  end subroutine submesh_broadcast

  ! --------------------------------------------------------------
  logical function submesh_compatible(this, radius, center, dx) result(compatible)
    type(submesh_t),   intent(in) :: this
    FLOAT,             intent(in) :: radius
    FLOAT,             intent(in) :: center(:)
    FLOAT,             intent(in) :: dx

    compatible =.false.
    if (allocated(this%center)) then
      !> At the moment the center check doesnt matter because everywhere submeshes are recycled at the moment,
      !> between uses, the position doesnt change
      if (radius <= this%radius .and. all(abs(this%center - center) < 0.25*dx)) then
        compatible = .true.
      end if
    end if
    return

  end function submesh_compatible

  ! --------------------------------------------------------------
  subroutine submesh_end(this)
    type(submesh_t),   intent(inout)  :: this

    PUSH_SUB(submesh_end)

    if (this%np /= -1) then
      nullify(this%mesh)
      this%np = -1
      SAFE_DEALLOCATE_A(this%center)
      SAFE_DEALLOCATE_A(this%map)
      SAFE_DEALLOCATE_A(this%rel_x)
      SAFE_DEALLOCATE_A(this%r)
      SAFE_DEALLOCATE_A(this%regions)
    end if

    if (accel_is_enabled()) then
      call accel_release_buffer(this%buff_map)
    end if

    POP_SUB(submesh_end)
  end subroutine submesh_end

  ! --------------------------------------------------------------

  subroutine submesh_get_inv(this, map_inv)
    type(submesh_t),      intent(in)   :: this
    integer,              intent(out)  :: map_inv(:)

    integer :: is

    PUSH_SUB(submesh_get_inv)

    map_inv(1:this%mesh%np_part) = 0
    do is = 1, this%np
      map_inv(this%map(is)) = is
    end do

    POP_SUB(submesh_get_inv)
  end subroutine submesh_get_inv

  ! --------------------------------------------------------------
  logical function submesh_overlap(sm1, sm2, space) result(overlap)
    type(submesh_t),      intent(in)   :: sm1
    type(submesh_t),      intent(in)   :: sm2
    type(space_t),        intent(in)   :: space

    integer :: ii, jj, dd
    FLOAT :: distance

    !no PUSH_SUB, called too often

    if (.not. space%is_periodic()) then
      !first check the distance
      distance = sum((sm1%center - sm2%center)**2)
      overlap = distance <= (CNST(1.5)*(sm1%radius + sm2%radius))**2

      ! if they are very far, no need to check in detail
      if (.not. overlap) return
    end if

    ! Otherwise check whether they have the some point in common. We
    ! can make the comparison faster using that the arrays are sorted.
    overlap = .false.
    ii = 1
    jj = 1
    do while(ii <= sm1%np .and. jj <= sm2%np)
      dd = sm1%map(ii) - sm2%map(jj)
      if (dd < 0) then
        ii = ii + 1
      else if (dd > 0) then
        jj = jj + 1
      else
        overlap = .true.
        exit
      end if
    end do

    if (sm1%mesh%parallel_in_domains) then
      call sm1%mesh%mpi_grp%allreduce_inplace(overlap, 1, MPI_LOGICAL, MPI_LOR)
    end if

  end function submesh_overlap

  ! -------------------------------------------------------------
  subroutine submesh_build_global(this, space)
    type(submesh_t),      intent(inout) :: this
    type(space_t),        intent(in)    :: space

    integer, allocatable :: part_np(:)
    integer :: ipart, ind, ip

    PUSH_SUB(submesh_build_global)

    if (.not. this%mesh%parallel_in_domains) then
      this%np_global = this%np
      POP_SUB(submesh_build_global)
      return
    end if

    SAFE_ALLOCATE(part_np(this%mesh%pv%npart))
    part_np = 0
    part_np(this%mesh%pv%partno) = this%np

    call this%mesh%allreduce(part_np)
    this%np_global = sum(part_np)

    SAFE_ALLOCATE(this%rel_x_global(1:space%dim, 1:this%np_global))
    SAFE_ALLOCATE(this%part_v(1:this%np_global))
    SAFE_ALLOCATE(this%global2local(1:this%np_global))
    this%rel_x_global(1:space%dim, 1:this%np_global) = M_ZERO
    this%part_v(1:this%np_global) = 0
    this%global2local(1:this%np_global) = 0

    ind = 0
    do ipart = 1, this%mesh%pv%npart
      if (ipart == this%mesh%pv%partno) then
        do ip = 1, this%np
          this%rel_x_global(:, ind + ip) = this%rel_x(:, ip)
          this%part_v(ind + ip) = this%mesh%pv%partno
          this%global2local(ind + ip) = ip
        end do
      end if
      ind = ind + part_np(ipart)
    end do

    call this%mesh%allreduce(this%rel_x_global)
    call this%mesh%allreduce(this%part_v)
    call this%mesh%allreduce(this%global2local)

    SAFE_DEALLOCATE_A(part_np)

    POP_SUB(submesh_build_global)
  end subroutine submesh_build_global

  ! -----------------------------------------------------------
  subroutine submesh_end_global(this)
    type(submesh_t),      intent(inout)   :: this

    PUSH_SUB(submesh_end_global)

    SAFE_DEALLOCATE_A(this%rel_x_global)
    this%np_global = -1
    SAFE_DEALLOCATE_A(this%part_v)
    SAFE_DEALLOCATE_A(this%global2local)

    POP_SUB(submesh_end_global)
  end subroutine submesh_end_global


  ! -----------------------------------------------------------
  subroutine zzsubmesh_add_to_mesh(this, sphi, phi, factor)
    type(submesh_t),  intent(in)    :: this
    CMPLX,            intent(in)    :: sphi(:)
    CMPLX,            intent(inout) :: phi(:)
    CMPLX,  optional, intent(in)    :: factor

    integer :: ip, m

    PUSH_SUB(zzdsubmesh_add_to_mesh)

    if (present(factor)) then
      !Loop unrolling inspired by BLAS axpy routine
      m = mod(this%np,4)
      do ip = 1, m
        phi(this%map(ip)) = phi(this%map(ip)) + factor*sphi(ip)
      end do
      if (this%np.ge.4) then
        do ip = m+1, this%np, 4
          phi(this%map(ip))   = phi(this%map(ip))   + factor*sphi(ip)
          phi(this%map(ip+1)) = phi(this%map(ip+1)) + factor*sphi(ip+1)
          phi(this%map(ip+2)) = phi(this%map(ip+2)) + factor*sphi(ip+2)
          phi(this%map(ip+3)) = phi(this%map(ip+3)) + factor*sphi(ip+3)
        end do
      end if
    else
      m = mod(this%np,4)
      do ip = 1, m
        phi(this%map(ip)) = phi(this%map(ip)) + sphi(ip)
      end do
      if (this%np.ge.4) then
        do ip = m+1, this%np, 4
          phi(this%map(ip))   = phi(this%map(ip))   + sphi(ip)
          phi(this%map(ip+1)) = phi(this%map(ip+1)) + sphi(ip+1)
          phi(this%map(ip+2)) = phi(this%map(ip+2)) + sphi(ip+2)
          phi(this%map(ip+3)) = phi(this%map(ip+3)) + sphi(ip+3)
        end do
      end if
    end if

    POP_SUB(zzdsubmesh_add_to_mesh)
  end subroutine zzsubmesh_add_to_mesh

  !------------------------------------------------------------
  CMPLX function zzsubmesh_to_mesh_dotp(this, sphi, phi, reduce) result(dotp)
    type(submesh_t),   intent(in) :: this
    CMPLX,             intent(in) :: sphi(:)
    CMPLX,             intent(in) :: phi(:)
    logical, optional, intent(in) :: reduce

    integer :: is, m, ip
    type(profile_t), save :: prof_sm_reduce

    PUSH_SUB(zzsubmesh_to_mesh_dotp)

    dotp = M_z0

    if (this%mesh%use_curvilinear) then
      do is = 1, this%np
        dotp = dotp + this%mesh%vol_pp(this%map(is))*phi(this%map(is))*conjg(sphi(is))
      end do
    else
      m = mod(this%np,4)
      do ip = 1, m
        dotp = dotp + phi(this%map(ip))*conjg(sphi(ip))
      end do
      if (this%np.ge.4) then
        do ip = m+1, this%np, 4
          dotp = dotp + phi(this%map(ip))*conjg(sphi(ip))     &
            + phi(this%map(ip+1))*conjg(sphi(ip+1)) &
            + phi(this%map(ip+2))*conjg(sphi(ip+2)) &
            + phi(this%map(ip+3))*conjg(sphi(ip+3))
        end do
      end if
      dotp = dotp*this%mesh%vol_pp(1)
    end if

    if (optional_default(reduce, .true.) .and. this%mesh%parallel_in_domains) then
      call profiling_in(prof_sm_reduce, "SM_REDUCE_DOTP")
      call this%mesh%allreduce(dotp)
      call profiling_out(prof_sm_reduce)
    end if

    POP_SUB(zzsubmesh_to_mesh_dotp)
  end function zzsubmesh_to_mesh_dotp

  !------------------------------------------------------------
  !> finds the dimension of a box containing the submesh
  subroutine submesh_get_cube_dim(sm, space, db)
    type(submesh_t), target,  intent(in)  :: sm
    type(space_t),            intent(in)  :: space
    integer,                  intent(out) :: db(1:space%dim)

    integer :: ip, idir
    FLOAT :: chi(space%dim)
#ifdef HAVE_MPI
    integer :: db_red(1:space%dim)
#endif

    PUSH_SUB(submesh_get_cube_dim)

    db = 1

    do ip = 1, sm%np
      chi = sm%mesh%coord_system%from_cartesian(sm%rel_x(:, ip))
      do idir = 1, space%dim
        db(idir) = max(db(idir), nint(abs(chi(idir))/sm%mesh%spacing(idir) + M_HALF))
      end do
    end do

    if(sm%mesh%parallel_in_domains) then
#ifdef HAVE_MPI
      call sm%mesh%mpi_grp%allreduce(db(1), db_red(1), space%dim, MPI_INTEGER, MPI_MAX)
      db(1:space%dim) = db_red(1:space%dim)
#endif
    end if

    db = 2 * db + 1

    POP_SUB(submesh_get_cube_dim)
  end subroutine submesh_get_cube_dim

  !------------------------------------------------------------
  subroutine submesh_init_cube_map(sm, space)
    type(submesh_t), target,  intent(inout) :: sm
    type(space_t),            intent(in)    :: space

    integer(i8) :: ip
    integer ::idir
    FLOAT :: chi(space%dim), shift(space%dim)

    PUSH_SUB(submesh_init_cube_map)

    sm%cube_map%nmap = sm%np

    SAFE_ALLOCATE(sm%cube_map%map(1:space%dim, 1:sm%cube_map%nmap))

    !The center of the submesh does not belong to the mesh
    !So we first need to find the closest grid point, and center the cube to it
    chi = sm%mesh%coord_system%from_cartesian(sm%center)
    do idir = 1, space%dim
      shift(idir) = nint(chi(idir)/sm%mesh%spacing(idir))*sm%mesh%spacing(idir)
    end do
    shift = sm%mesh%coord_system%to_cartesian(shift)
    shift = shift - sm%center

    do ip = 1, sm%cube_map%nmap
      chi = sm%mesh%coord_system%from_cartesian(sm%rel_x(:,ip) - shift)
      do idir = 1, space%dim
        sm%cube_map%map(idir, ip) = nint(chi(idir)/sm%mesh%spacing(idir))
      end do
    end do

    if (accel_is_enabled()) then
      call accel_create_buffer(sm%cube_map%map_buffer, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, sm%cube_map%nmap*space%dim)
      call accel_write_buffer(sm%cube_map%map_buffer, sm%cube_map%nmap*space%dim, sm%cube_map%map)
    end if

    POP_SUB(submesh_init_cube_map)
  end subroutine submesh_init_cube_map

  !------------------------------------------------------------
  subroutine submesh_end_cube_map(sm)
    type(submesh_t),   intent(inout)  :: sm

    PUSH_SUB(submesh_end_cube_map)

    call mesh_cube_map_end(sm%cube_map)

    POP_SUB(submesh_end_cube_map)
  end subroutine submesh_end_cube_map

  !> The following function takes a batch of functions defined in
  !! submesh (ss) and adds one of them to each of the mesh functions in
  !! other batch (mm).  Each one is multiplied by a factor given by the
  !! array factor.
  !! In this version of the routine, the submesh batch is real
  !! and the mesh batch is complex
  subroutine dzsubmesh_batch_add(this, ss, mm)
    type(submesh_t),  intent(in)    :: this
    class(batch_t),   intent(in)    :: ss
    class(batch_t),   intent(inout) :: mm

    integer :: ist, idim, jdim, is

    PUSH_SUB(dzsubmesh_batch_add)

    ASSERT(.not. mm%is_packed())
    ASSERT(ss%type() == TYPE_FLOAT)
    ASSERT(mm%type() == TYPE_CMPLX)
    ASSERT(ss%nst_linear == mm%nst_linear)
    ASSERT(ss%status() == mm%status())
    ASSERT(ss%dim == mm%dim)

    ASSERT(mm%nst == ss%nst)

    !$omp parallel do private(ist, idim, jdim, is) if(.not. this%overlap)
    do ist =  1, mm%nst
      do idim = 1, mm%dim
        jdim = min(idim, ss%dim)

        do is = 1, this%np
          mm%zff(this%map(is), idim, ist) = &
            mm%zff(this%map(is), idim, ist) + ss%dff(is, jdim, ist)
        end do

      end do
    end do
    !$omp end parallel do

    POP_SUB(dzsubmesh_batch_add)
  end subroutine dzsubmesh_batch_add


#include "undef.F90"
#include "real.F90"
#include "submesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "submesh_inc.F90"

end module submesh_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
