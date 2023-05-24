!! Copyright (C) 2022 F. Bonafé
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

module regridding_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none
  public :: grid_transfer_mapping

contains

  ! ---------------------------------------------------------
  !> Generate a re-mapping of points from grid_in to grid_out
  !> Works only on 3D and orthogonal axes
  subroutine grid_transfer_mapping(grid_out, grid_in, regridding_map, namespace)
    type(grid_t), intent(in)    :: grid_out
    type(grid_t), intent(in)    :: grid_in
    integer, allocatable, intent(out) :: regridding_map(:) ! it will be of size grid_in%mesh%np
    type(namespace_t),  intent(in) :: namespace

    integer :: ip, index(1:grid_in%coord_system%dim), ip_system
    FLOAT :: chi(1:grid_in%coord_system%dim)
    type(profile_t), save :: prof

    PUSH_SUB(grid_transfer_mapping)

    call profiling_in(prof,"GRID_TRANSFER_MAPPING")

    if (all(grid_out%mesh%spacing == grid_in%mesh%spacing) .and. &
        (.not.grid_out%mesh%parallel_in_domains .or. &
        grid_out%mesh%idx%checksum == grid_in%mesh%idx%checksum)) then
       SAFE_ALLOCATE(regridding_map(grid_in%mesh%np))
       regridding_map = 0

       do ip = 1, grid_in%mesh%np
          chi = grid_in%mesh%coord_system%from_cartesian(grid_in%mesh%x(ip,:))
          index = int(chi/grid_in%mesh%spacing(1:grid_in%coord_system%dim))
          ip_system = mesh_local_index_from_coords(grid_out%mesh, index)
          if (ip_system > grid_out%mesh%np) ip_system = 0
          regridding_map(ip) = ip_system
       end do
    else
      message(1) = 'Only meshes with equal spacing (in serial) or exactly equal meshes (in parallel)'
      message(2) = 'are supported for regridding'
      call messages_fatal(2, namespace=namespace)
    end if

    call profiling_out(prof)
    POP_SUB(grid_transfer_mapping)
  end subroutine grid_transfer_mapping

end module regridding_oct_m
