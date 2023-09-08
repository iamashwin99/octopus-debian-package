!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio,
!! G. Bertsch, M. Oliveira, J. Alberdi-Rodriguez, P. Garcia Risueño
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


! ---------------------------------------------------------
!> Allocates locally the real space grid, if PFFT library is not used.
!! Otherwise, it assigns the PFFT real space grid to the cube real space grid,
!! via pointer.
subroutine X(cube_function_alloc_rs)(cube, cf, in_device, force_alloc)
  type(cube_t), target,  intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf
  logical, optional,     intent(in)    :: in_device
  logical, optional,     intent(in)    :: force_alloc

  logical :: is_allocated

  PUSH_SUB(X(cube_function_alloc_rs))

  ASSERT(.not. associated(cf%X(rs)))

  is_allocated = .false.

  cf%forced_alloc = optional_default(force_alloc, .false.)

  if (allocated(cube%fft)) then
    select case (cube%fft%library)
    case (FFTLIB_PFFT)

      if (.not. cf%forced_alloc) then
        is_allocated = .true.
        cf%X(rs) => cube%fft%X(rs_data)(1:cube%rs_n(1), 1:cube%rs_n(2), 1:cube%rs_n(3))
      end if
    case (FFTLIB_ACCEL)
      if (optional_default(in_device, .true.)) then
        is_allocated = .true.
        cf%in_device_memory = .true.
        ! conversion to i8 needed to avoid overflow
        call accel_create_buffer(cf%real_space_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, product(int(cube%rs_n(1:3), i8)))
      end if
      !We use aligned memory for FFTW
    case (FFTLIB_FFTW)
      if (.not. cf%forced_alloc) then
        ASSERT(cube%fft%aligned_memory)
        is_allocated = .true.
        cf%X(rs) => cube%fft%X(rs_data)(1:cube%rs_n(1), 1:cube%rs_n(2), 1:cube%rs_n(3))
      end if
    end select
  end if

  if (.not. is_allocated) then
    SAFE_ALLOCATE(cf%X(rs)(1:cube%rs_n(1), 1:cube%rs_n(2), 1:cube%rs_n(3)))
  end if

  POP_SUB(X(cube_function_alloc_rs))
end subroutine X(cube_function_alloc_rs)


! ---------------------------------------------------------
!> Deallocates the real space grid
subroutine X(cube_function_free_rs)(cube, cf)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  logical :: deallocated

  PUSH_SUB(X(cube_function_free_rs))

  deallocated = .false.

  if (allocated(cube%fft)) then
    select case (cube%fft%library)
    case (FFTLIB_PFFT)
      if (.not. cf%forced_alloc) then
        deallocated = .true.
        nullify(cf%X(rs))
      end if
    case (FFTLIB_ACCEL)
      if (cf%in_device_memory) then
        deallocated = .true.
        ASSERT(cf%in_device_memory)
        call accel_release_buffer(cf%real_space_buffer)
        cf%in_device_memory = .false.
      end if
    case (FFTLIB_FFTW)
      if (.not. cf%forced_alloc) then
        deallocated = .true.
        nullify(cf%X(rs))
      end if
    end select
  end if

  if (.not. deallocated) then
    SAFE_DEALLOCATE_P(cf%X(rs))
  end if

  POP_SUB(X(cube_function_free_rs))
end subroutine X(cube_function_free_rs)

! ---------------------------------------------------------
subroutine X(cube_function_allgather)(cube, cf, cf_local, order, gatherfs)
  type(cube_t),      intent(in) :: cube
  R_TYPE,            intent(out):: cf(:,:,:)
  R_TYPE,            intent(in) :: cf_local(:,:,:)
  integer, optional, intent(in) :: order(3)
  logical, optional, intent(in) :: gatherfs

#ifdef HAVE_MPI
  integer :: ix, iy, iz, index, order_(1:3)
  integer(i8) :: number_points
  R_TYPE, allocatable :: cf_tmp(:)
  type(profile_t), save :: prof_allgather
#endif

  PUSH_SUB(X(cube_function_allgather))

#ifdef HAVE_MPI
  call profiling_in(prof_allgather, TOSTRING(X(CF_ALLGATHER)))

  if (cube_getFFTLibrary(cube) == FFTLIB_PFFT .or. &
    (cube_getFFTLibrary(cube) == FFTLIB_PNFFT .and. .not. optional_default(gatherfs, .false.))) then

    ! make sure we do not run into integer overflow here
    number_points = cube%rs_n_global(1) * cube%rs_n_global(2)
    number_points = number_points * cube%rs_n_global(3)
    if (number_points >= HUGE(0)) then
      message(1) = "Error: too many points for the normal cube. Please try to use a distributed FFT."
      call messages_fatal(1)
    end if
    SAFE_ALLOCATE(cf_tmp(1:cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3)))
    ! Warning: in the next line we have to pass the full cf_local array, not just the first element.
    ! This is because cf_local might be a pointer to a subarray when using PFFT, such that
    ! memory will not be contiguous (see cube_function_alloc_rs). In that case the
    ! Fortran compiler should do a temporary copy.
    call cube%mpi_grp%allgatherv(cf_local, cube%np_local(cube%mpi_grp%rank+1), R_MPITYPE, &
      cf_tmp, cube%np_local, cube%xlocal - 1, R_MPITYPE)

    if (optional_default(gatherfs,.false.)) then
      order_ = (/2,3,1/)
    else
      order_ = (/1,2,3/)
    end if
    if (present(order)) order_ = order


    do index = 1, cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3)
      ix = cube%local(index, order_(1))
      iy = cube%local(index, order_(2))
      iz = cube%local(index, order_(3))
      cf(ix, iy, iz) = cf_tmp(index)
    end do

  else
    order_ = (/1,2,3/)
    if (present(order)) order_ = order

    ! make sure we do not run into integer overflow here
    number_points = cube%fs_n_global(1) * cube%fs_n_global(2)
    number_points = number_points * cube%fs_n_global(3)
    if (number_points >= HUGE(0)) then
      message(1) = "Error: too many points for the normal cube. Please try to use a distributed FFT."
      call messages_fatal(1)
    end if
    SAFE_ALLOCATE(cf_tmp(1:cube%fs_n_global(1)*cube%fs_n_global(2)*cube%fs_n_global(3)))
    ! Warning: in the next line we have to pass the full cf_local array, not just the first element.
    ! This is because cf_local might be a pointer to a subarray when using PFFT, such that
    ! memory will not be contiguous (see cube_function_alloc_rs). In that case the
    ! Fortran compiler should do a temporary copy.
    call cube%mpi_grp%allgatherv(cf_local, cube%np_local_fs(cube%mpi_grp%rank+1), R_MPITYPE, &
      cf_tmp, cube%np_local_fs, cube%xlocal_fs - 1, R_MPITYPE)

    do index = 1, cube%fs_n_global(1)*cube%fs_n_global(2)*cube%fs_n_global(3)
      ix = cube%local_fs(index, order_(1))
      iy = cube%local_fs(index, order_(2))
      iz = cube%local_fs(index, order_(3))
      cf(ix, iy, iz) = cf_tmp(index)
    end do

  end if

  SAFE_DEALLOCATE_A(cf_tmp)

  call profiling_out(prof_allgather)

#else
  cf = cf_local
#endif
  POP_SUB(X(cube_function_allgather))
end subroutine X(cube_function_allgather)


! ---------------------------------------------------------
!> The next two subroutines convert a function between the normal
!! mesh and the cube.
!!
!! Note that the function in the mesh should be defined
!! globally, not just in a partition (when running in
!! parallel in real-space domains).
! ---------------------------------------------------------

subroutine X(mesh_to_cube)(mesh, mf, cube, cf)
  class(mesh_t),         intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:) !< function defined on the mesh, mesh%np points
  type(cube_t), target,  intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  integer :: ip, ix, iy, iz, iproc
  integer :: im, ii, nn
  integer :: nps(0:mesh%mpi_grp%size-1), nmaps(0:mesh%mpi_grp%size-1)
  integer :: bsize
  type(accel_mem_t) :: mf_buffer, recv_buffer, map_buffer
  R_TYPE, pointer :: recv(:), mf_pointer(:)
  integer, pointer :: map_pointer(:, :), map(:, :)
  integer :: requests(1:4), recv_proc, send_proc
  type(accel_kernel_t), save :: X(kernel)
  type(profile_t), save :: prof_m2c

  PUSH_SUB(X(mesh_to_cube))
  call profiling_in(prof_m2c, TOSTRING(X(MESH_TO_CUBE)))

  ASSERT(ubound(mf, dim=1) == mesh%np .or. ubound(mf, dim=1) == mesh%np_part)

  if (.not. cube%cube_map_present) then
    message(1) = "Internal error: cube_map needs to be initialized."
    message(2) = "Please call cube_init_cube_map."
    call messages_fatal(2)
  end if

  ! get mesh and map sizes
  call mesh%mpi_grp%allgather(mesh%np, 1, MPI_INTEGER, nps(0), 1, MPI_INTEGER)
  call mesh%mpi_grp%allgather(cube%cube_map%nmap, 1, MPI_INTEGER, nmaps(0), 1, MPI_INTEGER)

  if (.not. cf%in_device_memory) then

    ASSERT(associated(cf%X(rs)))

    !$omp parallel workshare
    cf%X(rs) = M_ZERO
    !$omp end parallel workshare

    ASSERT(allocated(cube%cube_map%map))
    ASSERT(mesh%box%dim <= 3)

    if (mesh%mpi_grp%size == 1) then
      recv => mf
      map => cube%cube_map%map
    else
      SAFE_ALLOCATE(recv(1:maxval(nps)))
      SAFE_ALLOCATE(map(1:5, 1:maxval(nmaps)))
    end if

    ! We basically get the mesh data and mapping from each process and execute
    ! the mapping locally. This is faster than using collective operations.

    ! communicate in ring pattern
    do iproc = 0, mesh%mpi_grp%size - 1
      ! first communicate mesh function and map
      recv_proc = mesh%mpi_grp%rank - iproc
      if (recv_proc < 0) recv_proc = recv_proc + mesh%mpi_grp%size
      send_proc = mesh%mpi_grp%rank + iproc
      if (send_proc >= mesh%mpi_grp%size) send_proc = send_proc - mesh%mpi_grp%size
      if (mesh%mpi_grp%size > 1) then
        call mesh%mpi_grp%irecv(recv, nps(recv_proc), R_MPITYPE, recv_proc, requests(1))
        call mesh%mpi_grp%isend(mf, mesh%np, R_MPITYPE, send_proc, requests(2))
        call mesh%mpi_grp%irecv(map, nmaps(recv_proc)*5, MPI_INTEGER, recv_proc, requests(3))
        call mesh%mpi_grp%isend(cube%cube_map%map, cube%cube_map%nmap*5, MPI_INTEGER, send_proc, requests(4))
        call mesh%mpi_grp%wait(4, requests)
      end if

      ! now do the mapping
      !$omp parallel do private(ix, iy, iz, ip, nn, ii)
      do im = 1, nmaps(recv_proc)
        ix = map(1, im) + cube%center(1)
        iy = map(2, im) + cube%center(2)
        iz = map(3, im) + cube%center(3)

        ip = map(MCM_POINT, im)
        nn = map(MCM_COUNT, im)
        do ii = 0, nn - 1
          cf%X(rs)(ix + ii, iy, iz) = recv(ip + ii)
        end do
      end do
      !$omp end parallel do
    end do

    if (mesh%mpi_grp%size == 1) then
      nullify(recv)
      nullify(map)
    else
      SAFE_DEALLOCATE_P(recv)
      SAFE_DEALLOCATE_P(map)
    end if

  else

    ! conversion to i8 needed to avoid overflow
    call accel_set_buffer_to_zero(cf%real_space_buffer, R_TYPE_VAL, product(int(cube%rs_n(1:3), i8)))

    call accel_create_buffer(mf_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, mesh%np)
    call accel_write_buffer(mf_buffer, mesh%np, mf)
    if (accel%cuda_mpi) then
      call accel_get_device_pointer(mf_pointer, mf_buffer, [mesh%np])
    else
      mf_pointer => mf
    end if

    if (accel%cuda_mpi) then
      call accel_get_device_pointer(map_pointer, cube%cube_map%map_buffer, [5, cube%cube_map%nmap])
    else
      map_pointer => cube%cube_map%map
    end if

    call accel_create_buffer(recv_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, maxval(nps))
    if (accel%cuda_mpi) then
      call accel_get_device_pointer(recv, recv_buffer, [maxval(nps)])
    else
      SAFE_ALLOCATE(recv(1:maxval(nps)))
    end if

    call accel_create_buffer(map_buffer, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, maxval(nmaps)*5)
    if (accel%cuda_mpi) then
      call accel_get_device_pointer(map, map_buffer, [5, maxval(nmaps)])
    else
      SAFE_ALLOCATE(map(1:5, 1:maxval(nmaps)))
    end if

    ! The following could also be done in a pipelined fashion to decrease the
    ! copy latency, but seems to be fast enough as it is right now.
    ! We basically get the mesh data and mapping from each process and execute
    ! the mapping locally. Especially with GPU-aware MPI, this keeps the data
    ! on the GPU and reduces communication overhead.

    ! communicate in ring pattern
    do iproc = 0, mesh%mpi_grp%size - 1
      ! first communicate mesh function and map
      recv_proc = mesh%mpi_grp%rank - iproc
      if (recv_proc < 0) recv_proc = recv_proc + mesh%mpi_grp%size
      send_proc = mesh%mpi_grp%rank + iproc
      if (send_proc >= mesh%mpi_grp%size) send_proc = send_proc - mesh%mpi_grp%size
      if (mesh%mpi_grp%size > 1) then
        call mesh%mpi_grp%irecv(recv, nps(recv_proc), R_MPITYPE, recv_proc, requests(1))
        call mesh%mpi_grp%isend(mf_pointer, mesh%np, R_MPITYPE, send_proc, requests(2))
        call mesh%mpi_grp%irecv(map, nmaps(recv_proc)*5, MPI_INTEGER, recv_proc, requests(3))
        call mesh%mpi_grp%isend(map_pointer, cube%cube_map%nmap*5, MPI_INTEGER, send_proc, requests(4))
        call mesh%mpi_grp%wait(4, requests)

        ! copy data to GPU if we cannot use CUDA-aware MPI
        if (.not. accel%cuda_mpi) then
          call accel_write_buffer(recv_buffer, nps(recv_proc), recv)
          call accel_write_buffer(map_buffer, nmaps(recv_proc)*5, map)
        end if
      end if

      ! then do the actual mapping
      call accel_kernel_start_call(X(kernel), 'mesh_to_cube.cl', TOSTRING(X(mesh_to_cube)), &
        flags = '-D' + R_TYPE_CL)

      call accel_set_kernel_arg(X(kernel), 0, nmaps(recv_proc))
      call accel_set_kernel_arg(X(kernel), 1, cube%fft%stride_rs(1))
      call accel_set_kernel_arg(X(kernel), 2, cube%fft%stride_rs(2))
      call accel_set_kernel_arg(X(kernel), 3, cube%fft%stride_rs(3))
      call accel_set_kernel_arg(X(kernel), 4, cube%center(1))
      call accel_set_kernel_arg(X(kernel), 5, cube%center(2))
      call accel_set_kernel_arg(X(kernel), 6, cube%center(3))
      if (mesh%mpi_grp%size > 1) then
        call accel_set_kernel_arg(X(kernel), 7, map_buffer)
        call accel_set_kernel_arg(X(kernel), 8, recv_buffer)
      else
        ! serial version
        call accel_set_kernel_arg(X(kernel), 7, cube%cube_map%map_buffer)
        call accel_set_kernel_arg(X(kernel), 8, mf_buffer)
      end if
      call accel_set_kernel_arg(X(kernel), 9, cf%real_space_buffer)

      bsize = accel_kernel_workgroup_size(X(kernel))

      call accel_kernel_run(X(kernel), (/pad(nmaps(recv_proc), bsize)/), (/bsize/))

      call accel_finish()
    end do
    call accel_release_buffer(mf_buffer)
    call accel_release_buffer(map_buffer)
    call accel_release_buffer(recv_buffer)

    nullify(mf_pointer)
    nullify(map_pointer)
    if (accel%cuda_mpi) then
      nullify(recv)
      nullify(map)
    else
      SAFE_DEALLOCATE_P(recv)
      SAFE_DEALLOCATE_P(map)
    end if
  end if

  call profiling_count_transfers(mesh%np, mf(1))

  call profiling_out(prof_m2c)
  POP_SUB(X(mesh_to_cube))
end subroutine X(mesh_to_cube)

! ---------------------------------------------------------
subroutine X(cube_to_mesh) (cube, cf, mesh, mf)
  type(cube_t),          intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  class(mesh_t),         intent(in)  :: mesh
  R_TYPE,  target,       intent(out) :: mf(:) !< function defined on the mesh, mesh%np points

  integer :: ip, ix, iy, iz
  integer :: im, ii, nn
  integer                    :: bsize
  type(accel_mem_t)         :: mf_buffer
  type(accel_kernel_t), save :: X(kernel)
  type(profile_t), save :: prof_c2m

  PUSH_SUB(X(cube_to_mesh))

  call profiling_in(prof_c2m, TOSTRING(X(CUBE_TO_MESH)))

  ASSERT(ubound(mf, dim=1) == mesh%np .or. ubound(mf, dim=1) == mesh%np_part)
  if (.not. cube%cube_map_present) then
    message(1) = "Internal error: cube_map needs to be initialized."
    message(2) = "Please call cube_init_cube_map."
    call messages_fatal(2)
  end if

  if (.not. cf%in_device_memory) then

    ASSERT(associated(cf%X(rs)))
    ASSERT(allocated(cube%cube_map%map))

    !$omp parallel do private(ix, iy, iz, ip, nn, ii)
    do im = 1, cube%cube_map%nmap
      ix = cube%cube_map%map(1, im) + cube%center(1)
      iy = cube%cube_map%map(2, im) + cube%center(2)
      iz = cube%cube_map%map(3, im) + cube%center(3)

      ip = cube%cube_map%map(MCM_POINT, im)
      nn = cube%cube_map%map(MCM_COUNT, im)

      do ii = 0, nn - 1
        mf(ip + ii) = cf%X(rs)(ix + ii, iy, iz)
      end do
    end do
    !$omp end parallel do

  else

    call accel_create_buffer(mf_buffer, ACCEL_MEM_WRITE_ONLY, R_TYPE_VAL, mesh%np)

    call accel_kernel_start_call(X(kernel), 'mesh_to_cube.cl', TOSTRING(X(cube_to_mesh)), flags = '-D' + R_TYPE_CL)

    call accel_set_kernel_arg(X(kernel), 0, cube%cube_map%nmap)
    call accel_set_kernel_arg(X(kernel), 1, cube%fft%stride_rs(1))
    call accel_set_kernel_arg(X(kernel), 2, cube%fft%stride_rs(2))
    call accel_set_kernel_arg(X(kernel), 3, cube%fft%stride_rs(3))
    call accel_set_kernel_arg(X(kernel), 4, cube%center(1))
    call accel_set_kernel_arg(X(kernel), 5, cube%center(2))
    call accel_set_kernel_arg(X(kernel), 6, cube%center(3))
    call accel_set_kernel_arg(X(kernel), 7, cube%cube_map%map_buffer)
    call accel_set_kernel_arg(X(kernel), 8, cf%real_space_buffer)
    call accel_set_kernel_arg(X(kernel), 9, mf_buffer)

    bsize = accel_kernel_workgroup_size(X(kernel))

    call accel_kernel_run(X(kernel), (/pad(cube%cube_map%nmap, bsize)/), (/bsize/))
    call accel_finish()

    call accel_read_buffer(mf_buffer, mesh%np, mf)
    call accel_release_buffer(mf_buffer)

  end if

  call profiling_count_transfers(mesh%np, mf(1))

  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh))

end subroutine X(cube_to_mesh)

! ---------------------------------------------------------
!> The next two subroutines convert a function between the normal
!! mesh and the cube in parallel.
! ---------------------------------------------------------
subroutine X(mesh_to_cube_parallel)(mesh, mf, cube, cf, map)
  class(mesh_t),         intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:)  !< mf(mesh%np)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf
  type(mesh_cube_parallel_map_t), intent(in) :: map

  integer :: ip, ix, iy, iz
  integer :: im, ii, nn
  integer :: min_x, min_y, min_z, max_x, max_y, max_z
  R_TYPE, allocatable :: in(:), out(:)
  type(profile_t), save :: prof_m2c

  PUSH_SUB(X(mesh_to_cube_parallel))
  call profiling_in(prof_m2c, TOSTRING(X(MESH_TO_CUBE_PARALLEL)))

  ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)
  if (.not. cube%cube_map_present) then
    message(1) = "Internal error: cube_map needs to be initialized."
    message(2) = "Please call cube_init_cube_map."
    call messages_fatal(2)
  end if

  if (mesh%parallel_in_domains) then
    SAFE_ALLOCATE(in(1:map%m2c_nsend))
    SAFE_ALLOCATE(out(1:map%m2c_nrec))

    !Put the mesh function data in correct order for transfer
    do ip = 1, map%m2c_nsend
      in(ip) = mf(map%m2c_mf_order(ip))
    end do

    !Transfer the mesh function from the mesh partition to the cube partition
    call X(partition_transfer)(map%m2c, in, out)

    !Put the transfered mesh function in the cube
    cf%X(rs) = R_TOTYPE(M_ZERO)
    do ip = 1, map%m2c_nrec
      ix = map%m2c_cf_order(ip, 1)
      iy = map%m2c_cf_order(ip, 2)
      iz = map%m2c_cf_order(ip, 3)
      cf%X(rs)(ix, iy, iz) = out(ip)
    end do

    SAFE_DEALLOCATE_A(in)
    SAFE_DEALLOCATE_A(out)

  else
    ! Save the limit values
    min_x = cube%rs_istart(1)
    min_y = cube%rs_istart(2)
    min_z = cube%rs_istart(3)
    max_x = cube%rs_istart(1) + cube%rs_n(1)
    max_y = cube%rs_istart(2) + cube%rs_n(2)
    max_z = cube%rs_istart(3) + cube%rs_n(3)

    ! Initialize to zero the input matrix
    do iz = 1, cube%rs_n(3)
      do iy = 1, cube%rs_n(2)
        do ix = 1, cube%rs_n(1)
          cf%X(rs)(ix, iy, iz) = M_ZERO
        end do
      end do
    end do

    ! Do the actual transform, only for the output values
    do im = 1, cube%cube_map%nmap
      ip = cube%cube_map%map(MCM_POINT, im)
      nn = cube%cube_map%map(MCM_COUNT, im)

      iz = cube%cube_map%map(3, im) + cube%center(3)
      if (iz >= min_z .and. iz < max_z) then
        iy = cube%cube_map%map(2, im) + cube%center(2)
        if (iy >= min_y .and. iy < max_y) then
          ix = cube%cube_map%map(1, im) + cube%center(1)
          do ii = 0, nn - 1
            if (ix+ii >= min_x .and. ix+ii < max_x) then
              cf%X(rs)(ix+ii-min_x+1, iy-min_y+1, iz-min_z+1) = mf(ip + ii)
            end if
          end do
        end if
      end if

    end do

  end if


  call profiling_count_transfers(mesh%np_global, mf(1))
  call profiling_out(prof_m2c)

  POP_SUB(X(mesh_to_cube_parallel))
end subroutine X(mesh_to_cube_parallel)

! ---------------------------------------------------------
subroutine X(cube_to_mesh_parallel) (cube, cf, mesh, mf, map)
  type(cube_t),          intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  class(mesh_t),         intent(in)  :: mesh
  R_TYPE,                intent(out) :: mf(:)  !< mf(mesh%np)
  type(mesh_cube_parallel_map_t), intent(in) :: map

  integer :: ip, ix, iy, iz, im, ii, nn, ixyz(3)
  R_TYPE, allocatable :: gcf(:,:,:)
  R_TYPE, allocatable :: in(:), out(:)
  type(profile_t), save :: prof_c2m

  PUSH_SUB(X(cube_to_mesh_parallel))
  call profiling_in(prof_c2m, TOSTRING(X(CUBE_TO_MESH_PARALLEL)))

  ASSERT(.not. cf%in_device_memory)
  ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)
  if (.not. cube%cube_map_present) then
    message(1) = "Internal error: cube_map needs to be initialized."
    message(2) = "Please call cube_init_cube_map."
    call messages_fatal(2)
  end if

  if (cube%parallel_in_domains) then
    SAFE_ALLOCATE(in(1:map%c2m_nsend))
    SAFE_ALLOCATE(out(1:map%c2m_nrec))

    !Put the cube function in the mesh in the correct order for transfer
    do ip = 1, map%c2m_nsend
      ix = map%c2m_cf_order(ip, 1)
      iy = map%c2m_cf_order(ip, 2)
      iz = map%c2m_cf_order(ip, 3)
      in(ip) = cf%X(rs)(ix, iy, iz)
    end do

    !Transfer the mesh function from the cube partition to the mesh partition
    call X(partition_transfer)(map%c2m, in, out)

    !Put the transfered mesh function data in correct order
    do ip = 1, map%c2m_nrec
      mf(map%c2m_mf_order(ip)) = out(ip)
    end do

    SAFE_DEALLOCATE_A(in)
    SAFE_DEALLOCATE_A(out)

  else
    !collect the data in all processes
    SAFE_ALLOCATE(gcf(1:cube%rs_n_global(1), 1:cube%rs_n_global(2), 1:cube%rs_n_global(3)))
    if (mpi_world%size > 1) then
#ifdef HAVE_MPI
      call X(cube_function_allgather)(cube, gcf, cf%X(rs))
#endif
    else
      gcf = cf%X(rs)
    end if

    ! cube to mesh
    do im = 1, cube%cube_map%nmap
      ip = cube%cube_map%map(MCM_POINT, im)
      nn = cube%cube_map%map(MCM_COUNT, im)

      ixyz(1:3) = cube%cube_map%map(1:3, im) + cube%center(1:3)

      do ii = 0, nn - 1
        mf(ip + ii) = gcf(ixyz(1)+ii, ixyz(2), ixyz(3))
      end do
    end do

    SAFE_DEALLOCATE_A(gcf)
  end if

  call profiling_count_transfers(mesh%np, mf(1))
  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh_parallel))
end subroutine X(cube_to_mesh_parallel)

! ---------------------------------------------------------
!> This function calculates the surface average of any function.
!! \warning Some more careful testing should be done on this.
R_TYPE function X(cube_function_surface_average)(cube, cf) result(x)
  type(cube_t),          intent(in) :: cube
  type(cube_function_t), intent(in) :: cf

  integer :: ii, jj, kk, ix, iy, iz, npoints
  R_TYPE :: tmp_x

  ASSERT(.not. cf%in_device_memory)

  PUSH_SUB(X(cube_function_surface_average))

  tmp_x = M_ZERO
  do ii = 1, cube%rs_n(1)
    do jj = 1, cube%rs_n(2)
      do kk = 1, cube%rs_n(3)
        ix = ii + cube%rs_istart(1) - 1
        iy = jj + cube%rs_istart(2) - 1
        iz = kk + cube%rs_istart(3) - 1
        if ((ix == 1 .or. ix == cube%rs_n_global(1)) .or. &
          ((iy == 1 .or. iy == cube%rs_n_global(2)) .and. (ix /= 1 .and. ix /= cube%rs_n_global(1))) .or. &
          ((iz == 1 .or. iz == cube%rs_n_global(3)) .and. (ix /= 1 .and. ix /= cube%rs_n_global(1) .and. &
          iy /= 1 .and. iy /= cube%rs_n_global(2)))) then
          tmp_x = tmp_x + cf%X(RS)(ii, jj, kk)
        end if
      end do
    end do
  end do


  if (cube%parallel_in_domains) then
    call cube%mpi_grp%allreduce(tmp_x, x, 1, R_MPITYPE, MPI_SUM)
  else
    x = tmp_x
  end if

  npoints = 2*(cube%rs_n_global(1)-2)**2 + 4*(cube%rs_n_global(1)-2) + &
    2*(cube%rs_n_global(2)-2)**2 + 4*(cube%rs_n_global(2)-2) + &
    2*(cube%rs_n_global(3)-2)**2 + 4*(cube%rs_n_global(3)-2) + 8
  x = x/npoints

  POP_SUB(X(cube_function_surface_average))
end function X(cube_function_surface_average)


!> The next two subroutines convert a function between a
!! submesh and the cube.
! ---------------------------------------------------------

subroutine X(submesh_to_cube)(sm, mf, cube, cf)
  type(submesh_t),       intent(in)    :: sm
  R_TYPE,  target,       intent(in)    :: mf(:) !< function defined on the submesh.
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  integer(i8) :: ix, iy, iz, im
  type(profile_t), save :: prof_sm2c

  PUSH_SUB(X(submesh_to_cube))
  call profiling_in(prof_sm2c, TOSTRING(X(SUBMESH_TO_CUBE)))

  ASSERT(ubound(mf, dim=1) == sm%np)

  !Not implemented, but not reachable at the moment
  ASSERT(.not. cf%in_device_memory)

  ASSERT(associated(cf%X(rs)))

  !$omp parallel workshare
  cf%X(rs) = M_ZERO
  !$omp end parallel workshare

  ASSERT(allocated(sm%cube_map%map))
  ASSERT(sm%mesh%box%dim <= 3)

  do im = 1, sm%np
    ix = sm%cube_map%map(1, im) + cube%center(1)
    iy = sm%cube_map%map(2, im) + cube%center(2)
    iz = sm%cube_map%map(3, im) + cube%center(3)
    cf%X(rs)(ix, iy, iz) = mf(im)
  end do

  if (sm%mesh%parallel_in_domains) then
    call sm%mesh%allreduce(cf%X(rs))
    call profiling_count_transfers(sm%np_global, mf(1))
  else
    call profiling_count_transfers(sm%np, mf(1))
  end if

  call profiling_out(prof_sm2c)
  POP_SUB(X(submesh_to_cube))
end subroutine X(submesh_to_cube)

! ---------------------------------------------------------
subroutine X(cube_to_submesh) (cube, cf, sm, mf)
  type(cube_t),          intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(submesh_t),       intent(in)  :: sm
  R_TYPE,  target,       intent(out) :: mf(:) !< function defined on the submesh.

  integer(i8) :: ix, iy, iz, im
  type(profile_t), save :: prof_c2sm

  PUSH_SUB(X(cube_to_submesh))

  call profiling_in(prof_c2sm, TOSTRING(X(CUBE_TO_SUBMESH)))

  ASSERT(ubound(mf, dim=1) == sm%np)

  !Not implemented, but not reachable at the moment
  ASSERT(.not. cf%in_device_memory)

  ASSERT(associated(cf%X(rs)))

  ASSERT(allocated(sm%cube_map%map))

  do im = 1, sm%np
    ix = sm%cube_map%map(1, im) + cube%center(1)
    iy = sm%cube_map%map(2, im) + cube%center(2)
    iz = sm%cube_map%map(3, im) + cube%center(3)
    mf(im) = cf%X(rs)(ix, iy, iz)
  end do

  call profiling_count_transfers(sm%np_global, mf(1))

  call profiling_out(prof_c2sm)

  POP_SUB(X(cube_to_submesh))

end subroutine X(cube_to_submesh)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
