!! Copyright (C) 2012-2014 M. Oliveira, J. Alberdi-Rodriguez
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

module mesh_cube_parallel_map_oct_m
  use cube_oct_m
  use debug_oct_m
  use global_oct_m
  use index_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use messages_oct_m
  use mpi_oct_m
  use partition_oct_m
  use partition_transfer_oct_m
  use par_vec_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                       &
    mesh_cube_parallel_map_t,     &
    mesh_cube_parallel_map_init,  &
    mesh_cube_parallel_map_end

  type mesh_cube_parallel_map_t
    ! Components are public by default
    logical :: initialized = .false.

    ! Mesh to cube:
    type(partition_transfer_t) :: m2c
    integer :: m2c_nsend !< How many points will this process send to the cube partition
    integer :: m2c_nrec  !< How many points will this process receive from the mesh partition
    integer, allocatable :: m2c_mf_order(:)   !< How the points of the mesh function should be ordered to execute
    !!                                           the transfer from the mesh partition to the cube partition.
    !!                                           These are *local* indices.
    integer, allocatable :: m2c_cf_order(:,:) !< How the points of the mesh function are ordered after executing
    !!                                           the transfer from the mesh partition to the cube partition with
    !!                                           respect to the *local* cube indices.

    ! Cube to mesh:
    type(partition_transfer_t) :: c2m
    integer :: c2m_nsend !< How many points will this process send to the mesh partition
    integer :: c2m_nrec  !< How many points will this process receive from the cube partition
    integer, allocatable :: c2m_cf_order(:,:) !< How the points of the mesh function should be ordered to execute
    !!                                           the transfer from the cube partition to the mesh partition with
    !!                                           respect to the *local* cube indices.
    integer, allocatable :: c2m_mf_order(:)   !< How the points of the mesh function are ordered after executing
    !!                                           the transfer from the cube partition to the mesh partition.
    !!                                           These are *local* indices.
  end type mesh_cube_parallel_map_t

contains

  ! ---------------------------------------------------------
  subroutine mesh_cube_parallel_map_init(this, mesh, cube)
    type(mesh_cube_parallel_map_t), intent(out) :: this
    type(mesh_t),                   intent(in)  :: mesh
    type(cube_t),                   intent(in)  :: cube

    integer :: ip, ixyz(3), lxyz(3), ipos, cube_np
    integer :: ix, iy, iz
    integer, allocatable :: cube_part_local(:)
    integer(i8), allocatable :: global_index(:)
    integer(i8), allocatable :: mf_order(:), cf_order(:)
    integer(i8) :: ipg
    type(dimensions_t), allocatable :: part(:)

    type(profile_t), save :: prof

    PUSH_SUB(mesh_cube_parallel_map_init)
    call profiling_in(prof,"MC_PAR_INIT")

    !Get the cube partition on the mesh and the number of local cube points that are also on the mesh
    SAFE_ALLOCATE(part(1:cube%mpi_grp%size))
    call cube_partition(cube, part)

    ! Mesh to cube
    ! We will work only with the local mesh points and we need to know the global index of those points.
    SAFE_ALLOCATE(cube_part_local(1:mesh%np))
    SAFE_ALLOCATE(global_index(1:mesh%np))
    do ip = 1, mesh%np
      global_index(ip) = mesh_local2global(mesh, ip)
      call mesh_local_index_to_coords(mesh, ip, ixyz)
      ixyz = ixyz + cube%center
      cube_part_local(ip) = cube_point_to_process(ixyz, part)
    end do

    ! Init partition transfer
    call partition_transfer_init(this%m2c, mesh%np, global_index, mesh%mpi_grp, &
      cube%mpi_grp, cube_part_local, &
      this%m2c_nsend, this%m2c_nrec, mf_order, cf_order)
    ! Convert ordering of mesh and cube points from global mesh index to local mesh and cube indices
    SAFE_ALLOCATE(this%m2c_mf_order(1:this%m2c_nsend))
    SAFE_ALLOCATE(this%m2c_cf_order(1:this%m2c_nrec, 1:3))

    ! Initialize all to 0, to detect possible errors
    this%m2c_mf_order = 0

    do ip = 1, this%m2c_nsend
      this%m2c_mf_order(ip) = mesh_global2local(mesh, mf_order(ip))
      if (this%m2c_mf_order(ip) == 0) then
        write(message(1),'(a,i4,a,i4)') "Error in mesh_cube_parallel_map_init (m2c): mesh point ", &
          mf_order(ip), " is not stored in partition ", mesh%pv%partno
        call messages_fatal(1)
      end if
    end do

    do ip = 1, this%m2c_nrec
      call mesh_global_index_to_coords(mesh, cf_order(ip), ixyz)
      ixyz = ixyz + cube%center

      if (.not. cube_global2local(cube, ixyz, lxyz)) then
        write(message(1),'(a,3i4,a,i4)') "Error in mesh_cube_parallel_map_init (m2c): cube point ", &
          lxyz(1:3), " is not stored in partition ", cube%mpi_grp%rank + 1
        call messages_fatal(1)
      end if

      this%m2c_cf_order(ip, 1:3) = lxyz(1:3)
    end do
    SAFE_DEALLOCATE_A(mf_order)
    SAFE_DEALLOCATE_A(cf_order)
    SAFE_DEALLOCATE_A(cube_part_local)
    SAFE_DEALLOCATE_A(global_index)


    ! Cube to mesh
    ASSERT(product(i4_to_i8(cube%rs_n)) < huge(0_i4))
    cube_np = product(cube%rs_n)

    ! We will work only with the local cube points and we need to know the global index of those points.
    SAFE_ALLOCATE(cube_part_local(1:cube_np))
    SAFE_ALLOCATE(global_index(1:cube_np))
    ipos = 0
    do iz = part(cube%mpi_grp%rank+1)%start_xyz(3), part(cube%mpi_grp%rank+1)%end_xyz(3)
      do iy = part(cube%mpi_grp%rank+1)%start_xyz(2), part(cube%mpi_grp%rank+1)%end_xyz(2)
        do ix = part(cube%mpi_grp%rank+1)%start_xyz(1), part(cube%mpi_grp%rank+1)%end_xyz(1)
          ixyz(1) = ix
          ixyz(2) = iy
          ixyz(3) = iz
          ixyz = ixyz - cube%center
          ipg = mesh_global_index_from_coords(mesh, ixyz)
          ! do not map boundary points
          if (ipg > 0 .and. ipg <= mesh%np_global) then
            ipos = ipos + 1
            global_index(ipos) = ipg
          end if
        end do
      end do
    end do
    cube_np = ipos

    if (mesh%parallel_in_domains) then
      call partition_get_partition_number(mesh%partition, cube_np, global_index, cube_part_local)
    else
      cube_part_local = 1
    end if

    SAFE_DEALLOCATE_A(part)

    ! Init partition transfer
    call partition_transfer_init(this%c2m, cube_np, global_index, cube%mpi_grp, &
      mesh%mpi_grp, cube_part_local, &
      this%c2m_nsend, this%c2m_nrec, cf_order, mf_order)

    ! Convert ordering of mesh and cube points from global mesh index to local mesh and cube indices
    SAFE_ALLOCATE(this%c2m_cf_order(1:this%c2m_nsend, 1:3))
    SAFE_ALLOCATE(this%c2m_mf_order(1:this%c2m_nrec))
    do ip = 1, this%c2m_nsend
      call mesh_global_index_to_coords(mesh, cf_order(ip), ixyz)
      ixyz = ixyz + cube%center

      if (.not. cube_global2local(cube, ixyz, lxyz)) then
        write(message(1),'(a,3i4,a,i4)') "Error in mesh_cube_parallel_map_init (c2m): cube point ", &
          lxyz(1:3), " is not stored in partition ", cube%mpi_grp%rank + 1
        call messages_fatal(1)
      end if
      this%c2m_cf_order(ip, 1:3) = lxyz(1:3)
    end do
    do ip = 1, this%c2m_nrec
      this%c2m_mf_order(ip) = mesh_global2local(mesh, mf_order(ip))
      if (this%c2m_mf_order(ip) == 0) then
        write(message(1),'(a,i3,a,i3)') "Error in mesh_cube_parallel_map_init (c2m): mesh point ", &
          mf_order(ip), " is not stored in partition ", mesh%pv%partno
        call messages_fatal(1)
      end if
    end do

    SAFE_DEALLOCATE_A(mf_order)
    SAFE_DEALLOCATE_A(cf_order)
    SAFE_DEALLOCATE_A(cube_part_local)
    SAFE_DEALLOCATE_A(global_index)

    this%initialized = .true.

    call profiling_out(prof)
    POP_SUB(mesh_cube_parallel_map_init)
  end subroutine mesh_cube_parallel_map_init

  ! ---------------------------------------------------------
  subroutine mesh_cube_parallel_map_end(this)
    type(mesh_cube_parallel_map_t), intent(inout) :: this

    PUSH_SUB(mesh_cube_parallel_map_end)

    if (this%initialized) then
      SAFE_DEALLOCATE_A(this%m2c_mf_order)
      SAFE_DEALLOCATE_A(this%m2c_cf_order)
      SAFE_DEALLOCATE_A(this%c2m_mf_order)
      SAFE_DEALLOCATE_A(this%c2m_cf_order)
      call partition_transfer_end(this%m2c)
      call partition_transfer_end(this%c2m)
      this%initialized = .false.
    end if

    POP_SUB(mesh_cube_parallel_map_end)
  end subroutine mesh_cube_parallel_map_end

end module mesh_cube_parallel_map_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
