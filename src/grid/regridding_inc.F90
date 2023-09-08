!! Copyright (C) 2022 S. Ohlmann
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

subroutine X(regridding_do_transfer_1)(this, ff_out, ff_in)
  class(regridding_t), intent(in)    :: this
  R_TYPE,              intent(inout) :: ff_out(:)
  R_TYPE,              intent(in)    :: ff_in(:)

  type(profile_t), save :: prof
  integer :: ip_in, ip_out, ip, is
  R_TYPE, allocatable :: out_buffer(:), in_buffer(:)
  integer :: index(1:this%dim), index_stencil(1:this%dim)

  PUSH_SUB(X(regridding_do_transfer_1))

  call profiling_in(prof,TOSTRING(X(REGRIDDING_DO_TRANSFER_1)))

  ASSERT(size(ff_out, dim=1) >= this%mesh_out%np)
  ASSERT(size(ff_in, dim=1) >= this%mesh_in%np)

  ff_out(:) = M_ZERO

  SAFE_ALLOCATE(out_buffer(1:this%nrecv))
  SAFE_ALLOCATE(in_buffer(1:this%nsend))
  in_buffer = M_ZERO
  if (this%do_restriction) then
    ! apply restriction operator
    do ip_in = 1, this%nsend
      ! get the index on the coarse mesh
      call mesh_global_index_to_coords(this%mesh_out, this%order_in_global(ip_in), index)
      index = index * this%eta
      ! now apply the stencil over the adjacent fine mesh points
      do is = 1, this%transfer_stencil%size
        index_stencil(:) = index(:) + this%transfer_stencil%points(1:this%dim, is)
        ip = mesh_local_index_from_coords(this%mesh_in, index_stencil)
        if (ip == 0 .or. ip > this%mesh_in%np) cycle
        ! save to coarse mesh point, to be transferred
        in_buffer(ip_in) = in_buffer(ip_in) + ff_in(ip) * this%weights(is)
      end do
    end do
  else
    ! copy in the right order to the buffer
    do ip_in = 1, this%nsend
      in_buffer(ip_in) = ff_in(this%order_in(ip_in))
    end do
  end if
  ! do the transfer
  call X(partition_transfer)(this%partition_transfer, in_buffer, out_buffer)
  if (this%do_prolongation) then
    ! apply prolongation operator
    do ip_out = 1, this%nrecv
      ! get the index on the coarse mesh
      call mesh_global_index_to_coords(this%mesh_in, this%order_out_global(ip_out), index)
      ! translate to fine mesh
      index = index * this%eta
      ! now apply stencil to adjacent fine mesh points
      do is = 1, this%transfer_stencil%size
        index_stencil(:) = index(:) + this%transfer_stencil%points(1:this%dim, is)
        ip = mesh_local_index_from_coords(this%mesh_out, index_stencil)
        if (ip == 0 .or. ip > this%mesh_out%np) cycle
        ! save to fine mesh point
        ff_out(ip) = ff_out(ip) + out_buffer(ip_out) * this%weights(is)
      end do
    end do
  else
    ! copy back in the correct order, do a reduction in case of restriction
    do ip_out = 1, this%nrecv
      ff_out(this%order_out(ip_out)) = ff_out(this%order_out(ip_out)) + out_buffer(ip_out)
    end do
  end if

  SAFE_DEALLOCATE_A(out_buffer)
  SAFE_DEALLOCATE_A(in_buffer)

  call profiling_out(prof)
  POP_SUB(X(regridding_do_transfer_1))
end subroutine X(regridding_do_transfer_1)

subroutine X(regridding_do_transfer_2)(this, ff_out, ff_in)
  class(regridding_t), intent(in)    :: this
  R_TYPE,              intent(inout) :: ff_out(:, :)
  R_TYPE,              intent(in)    :: ff_in(:, :)

  type(profile_t), save :: prof
  integer :: idim

  PUSH_SUB(X(regridding_do_transfer_2))

  call profiling_in(prof,TOSTRING(X(REGRIDDING_DO_TRANSFER_2)))

  ASSERT(size(ff_out, dim=2) == size(ff_in, dim=2))

  do idim = 1, size(ff_out, dim=2)
    call this%do_transfer(ff_out(:, idim), ff_in(:, idim))
  end do

  call profiling_out(prof)
  POP_SUB(X(regridding_do_transfer_2))
end subroutine X(regridding_do_transfer_2)
