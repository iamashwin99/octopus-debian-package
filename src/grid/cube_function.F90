!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module cube_function_oct_m
  use accel_oct_m
  use comm_oct_m
  use cube_oct_m
  use debug_oct_m
  use fft_oct_m
  use global_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use mesh_cube_parallel_map_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use partition_transfer_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use submesh_oct_m
  use types_oct_m

  implicit none
  private
  public ::                        &
    cube_function_t,               &
    dcube_function_alloc_RS,       &
    zcube_function_alloc_RS,       &
    dcube_function_free_RS,        &
    zcube_function_free_RS,        &
    dcube_function_surface_average,&
    zcube_function_surface_average,&
    dmesh_to_cube_parallel,        &
    zmesh_to_cube_parallel,        &
    dcube_to_mesh_parallel,        &
    zcube_to_mesh_parallel,        &
    dmesh_to_cube,                 &
    zmesh_to_cube,                 &
    dcube_to_mesh,                 &
    zcube_to_mesh,                 &
    dsubmesh_to_cube,              &
    zsubmesh_to_cube,              &
    dcube_to_submesh,              &
    zcube_to_submesh,              &
    dcube_function_allgather,      &
    zcube_function_allgather

  type cube_function_t
    ! Components are public by default
    FLOAT, pointer :: dRS(:, :, :) => NULL() !< real-space grid
    CMPLX, pointer :: zRS(:, :, :) => NULL() !< real-space grid, complex numbers
    CMPLX, pointer :: FS(:, :, :)  => NULL() !< Fourier-space grid
    logical            :: forced_alloc = .false. !< Forced to be allocated even when PFFT is associated with the cube
    logical            :: in_device_memory = .false.
    type(accel_mem_t) :: real_space_buffer
    type(accel_mem_t) :: fourier_space_buffer
  end type cube_function_t

contains

#include "undef.F90"
#include "real.F90"
#include "cube_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "cube_function_inc.F90"

end module cube_function_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
