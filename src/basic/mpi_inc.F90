!! Copyright (C) 2005-2006 Heiko Appel, Florian Lorenzen
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


! ---------------------------------------------------------
! these copy routines are for the serial case
subroutine X(mpi_grp_copy_0)(sendbuf, recvbuf, count)
  use iso_c_binding
  R_TYPE, target,   intent(in)  :: sendbuf
  R_TYPE, target,   intent(out) :: recvbuf
  integer,          intent(in)  :: count
  integer :: ii
  R_TYPE, pointer :: send(:), recv(:)

  call c_f_pointer(c_loc(sendbuf), send, [count])
  call c_f_pointer(c_loc(recvbuf), recv, [count])
  do ii = 1, count
    recv(ii) = send(ii)
  end do
end subroutine X(mpi_grp_copy_0)

subroutine X(mpi_grp_copy_1)(sendbuf, recvbuf, count)
  use iso_c_binding
  R_TYPE, target,   intent(in)  :: sendbuf(:)
  R_TYPE, target,   intent(out) :: recvbuf(:)
  integer,          intent(in)  :: count
  integer :: ii
  R_TYPE, pointer :: send(:), recv(:)

  call c_f_pointer(c_loc(sendbuf), send, [count])
  call c_f_pointer(c_loc(recvbuf), recv, [count])
  do ii = 1, count
    recv(ii) = send(ii)
  end do
end subroutine X(mpi_grp_copy_1)

subroutine X(mpi_grp_copy_2)(sendbuf, recvbuf, count)
  use iso_c_binding
  R_TYPE, target,   intent(in)  :: sendbuf(:, :)
  R_TYPE, target,   intent(out) :: recvbuf(:, :)
  integer,          intent(in)  :: count
  integer :: ii
  R_TYPE, pointer :: send(:), recv(:)

  call c_f_pointer(c_loc(sendbuf), send, [count])
  call c_f_pointer(c_loc(recvbuf), recv, [count])
  do ii = 1, count
    recv(ii) = send(ii)
  end do
end subroutine X(mpi_grp_copy_2)

subroutine X(mpi_grp_copy_3)(sendbuf, recvbuf, count)
  use iso_c_binding
  R_TYPE, target,   intent(in)  :: sendbuf(:, :, :)
  R_TYPE, target,   intent(out) :: recvbuf(:, :, :)
  integer,          intent(in)  :: count
  integer :: ii
  R_TYPE, pointer :: send(:), recv(:)

  call c_f_pointer(c_loc(sendbuf), send, [count])
  call c_f_pointer(c_loc(recvbuf), recv, [count])
  do ii = 1, count
    recv(ii) = send(ii)
  end do
end subroutine X(mpi_grp_copy_3)

subroutine X(mpi_grp_copy_3_1)(sendbuf, recvbuf, count)
  use iso_c_binding
  R_TYPE, target,   intent(in)  :: sendbuf(:, :, :)
  R_TYPE, target,   intent(out) :: recvbuf(:)
  integer,          intent(in)  :: count
  integer :: ii
  R_TYPE, pointer :: send(:), recv(:)

  call c_f_pointer(c_loc(sendbuf), send, [count])
  call c_f_pointer(c_loc(recvbuf), recv, [count])
  do ii = 1, count
    recv(ii) = send(ii)
  end do
end subroutine X(mpi_grp_copy_3_1)

! ---------------------------------------------------------
subroutine X(mpi_grp_scatterv)(mpi_grp, sendbuf, sendcnts, displs, sendtype, recvbuf, &
  recvcount, recvtype, root)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts(:), displs(:), sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount
  integer,          intent(in)  :: recvtype, root

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_SCATTERV)
  if (mpi_grp%comm /= -1) then
    call MPI_Scatterv(sendbuf, sendcnts, displs, sendtype, recvbuf, &
      recvcount, recvtype, root, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_SCATTERV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
#endif
end subroutine X(mpi_grp_scatterv)

! ---------------------------------------------------------
subroutine X(mpi_grp_gatherv)(mpi_grp, sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype, root)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts, sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount(:), displs(:)
  integer,          intent(in)  :: recvtype, root

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_GATHERV)
  if (mpi_grp%comm /= -1) then
    call MPI_Gatherv(sendbuf, sendcnts, sendtype, recvbuf, &
      recvcount, displs, recvtype, root, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_GATHERV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_gatherv)

! ---------------------------------------------------------
subroutine X(mpi_grp_gather_0)(mpi_grp, sendbuf, sendcount, sendtype, recvbuf, &
  recvcount, recvtype, root)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf
  integer,          intent(in)  :: sendcount, sendtype
  R_TYPE,           intent(out) :: recvbuf
  integer,          intent(in)  :: recvcount, recvtype, root

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_GATHER)
  if (mpi_grp%comm /= -1) then
    call MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, &
      recvcount, recvtype, root, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_0)(sendbuf, recvbuf, recvcount)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_GATHER)
#else
  call X(mpi_grp_copy_0)(sendbuf, recvbuf, recvcount)
#endif
end subroutine X(mpi_grp_gather_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_alltoallv)(mpi_grp, sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
  recvcount, rdispls, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts(:), sdispls(:), sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount(:), rdispls(:), recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLTOALLV)
  if (mpi_grp%comm /= -1) then
    call MPI_Alltoallv(sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
      recvcount, rdispls, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLTOALLV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_alltoallv)

! ---------------------------------------------------------
subroutine X(mpi_grp_alltoallv_2)(mpi_grp, sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
  recvcount, rdispls, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:, :)
  integer,          intent(in)  :: sendcnts(:), sdispls(:), sendtype
  R_TYPE,           intent(out) :: recvbuf(:, :)
  integer,          intent(in)  :: recvcount(:), rdispls(:), recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLTOALLV)
  if (mpi_grp%comm /= -1) then
    call MPI_Alltoallv(sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
      recvcount, rdispls, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_2)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLTOALLV)
#else
  call X(mpi_grp_copy_2)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_alltoallv_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_alltoallv_3)(mpi_grp, sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
  recvcount, rdispls, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:, :, :)
  integer,          intent(in)  :: sendcnts(:), sdispls(:), sendtype
  R_TYPE,           intent(out) :: recvbuf(:, :, :)
  integer,          intent(in)  :: recvcount(:), rdispls(:), recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLTOALLV)
  if (mpi_grp%comm /= -1) then
    call MPI_Alltoallv(sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
      recvcount, rdispls, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_3)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLTOALLV)
#else
  call X(mpi_grp_copy_3)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_alltoallv_3)

! ---------------------------------------------------------
subroutine X(mpi_grp_alltoall)(mpi_grp, sendbuf, sendcount, sendtype, recvbuf, &
  recvcount, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcount, sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount, recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLTOALLV)
  if (mpi_grp%comm /= -1) then
    call MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, &
      recvcount, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLTOALLV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
#endif
end subroutine X(mpi_grp_alltoall)

! ---------------------------------------------------------
subroutine X(mpi_grp_allgatherv)(mpi_grp, sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts, sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount(:), displs(:)
  integer,          intent(in)  :: recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHERV)
  if (mpi_grp%comm /= -1) then
    call MPI_Allgatherv(sendbuf, sendcnts, sendtype, recvbuf, &
      recvcount, displs, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHERV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_allgatherv)

! ---------------------------------------------------------
subroutine X(mpi_grp_allgatherv_2)(mpi_grp, sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:, :)
  integer,          intent(in)  :: sendcnts, sendtype
  R_TYPE,           intent(out) :: recvbuf(:, :)
  integer,          intent(in)  :: recvcount(:), displs(:)
  integer,          intent(in)  :: recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHERV)
  if (mpi_grp%comm /= -1) then
    call MPI_Allgatherv(sendbuf, sendcnts, sendtype, recvbuf, &
      recvcount, displs, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_2)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHERV)
#else
  call X(mpi_grp_copy_2)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_allgatherv_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_allgatherv_3)(mpi_grp, sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:, :, :)
  integer,          intent(in)  :: sendcnts, sendtype
  R_TYPE,           intent(out) :: recvbuf(:, :, :)
  integer,          intent(in)  :: recvcount(:), displs(:)
  integer,          intent(in)  :: recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHERV)
  if (mpi_grp%comm /= -1) then
    call MPI_Allgatherv(sendbuf, sendcnts, sendtype, recvbuf, &
      recvcount, displs, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_3)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHERV)
#else
  call X(mpi_grp_copy_3)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_allgatherv_3)

! ---------------------------------------------------------
subroutine X(mpi_grp_allgatherv_3_1)(mpi_grp, sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:, :, :)
  integer,          intent(in)  :: sendcnts, sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount(:), displs(:)
  integer,          intent(in)  :: recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHERV)
  if (mpi_grp%comm /= -1) then
    call MPI_Allgatherv(sendbuf, sendcnts, sendtype, recvbuf, &
      recvcount, displs, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_3_1)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHERV)
#else
  call X(mpi_grp_copy_3_1)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_allgatherv_3_1)

! ---------------------------------------------------------
subroutine X(mpi_grp_bcast)(mpi_grp, buf, cnt, sendtype, root)
  class(mpi_grp_t), intent(in)    :: mpi_grp
  R_TYPE,           intent(inout) :: buf(:)
  integer,          intent(in)    :: cnt, sendtype, root

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_BCAST)
  if (mpi_grp%comm /= -1) then
    call MPI_Bcast(buf, cnt, sendtype, root, mpi_grp%comm, mpi_err)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_BCAST)
#endif
end subroutine X(mpi_grp_bcast)

! ---------------------------------------------------------
subroutine X(mpi_grp_bcast_0)(mpi_grp, buf, cnt, sendtype, root)
  class(mpi_grp_t), intent(in)    :: mpi_grp
  R_TYPE,           intent(inout) :: buf
  integer,          intent(in)    :: cnt, sendtype, root

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_BCAST)
  if (mpi_grp%comm /= -1) then
    call MPI_Bcast(buf, cnt, sendtype, root, mpi_grp%comm, mpi_err)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_BCAST)
#endif
end subroutine X(mpi_grp_bcast_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_bcast_2)(mpi_grp, buf, cnt, sendtype, root)
  class(mpi_grp_t), intent(in)    :: mpi_grp
  R_TYPE,           intent(inout) :: buf(:, :)
  integer,          intent(in)    :: cnt, sendtype, root

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_BCAST)
  if (mpi_grp%comm /= -1) then
    call MPI_Bcast(buf, cnt, sendtype, root, mpi_grp%comm, mpi_err)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_BCAST)
#endif
end subroutine X(mpi_grp_bcast_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_bcast_3)(mpi_grp, buf, cnt, sendtype, root)
  class(mpi_grp_t), intent(in)    :: mpi_grp
  R_TYPE,           intent(inout) :: buf(:, :, :)
  integer,          intent(in)    :: cnt, sendtype, root

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_BCAST)
  if (mpi_grp%comm /= -1) then
    call MPI_Bcast(buf, cnt, sendtype, root, mpi_grp%comm, mpi_err)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_BCAST)
#endif
end subroutine X(mpi_grp_bcast_3)

! ---------------------------------------------------------
subroutine X(mpi_grp_bcast_0_l)(mpi_grp, buf, cnt, sendtype, root)
  use iso_c_binding
  class(mpi_grp_t), intent(in)    :: mpi_grp
  R_TYPE, target,   intent(inout) :: buf
  integer(i8),      intent(in)    :: cnt
  integer,          intent(in)    :: sendtype, root

#if defined(HAVE_MPI)
  integer :: rounds, iround, size
  integer(i8) :: offset
  R_TYPE, pointer :: bufptr(:)

  call mpi_debug_in(mpi_grp%comm, C_MPI_BCAST)
  if (mpi_grp%comm /= -1) then
    ! need to do the broadcast in rounds that fit into i4 integers
    call c_f_pointer(c_loc(buf), bufptr, [cnt])
    rounds = int(cnt/huge(0_i4), i4)
    do iround = 1, rounds
      offset = int(huge(0_i4), i8) * (iround - 1) + 1
      call MPI_Bcast(bufptr(offset), huge(0_i4), sendtype, root, mpi_grp%comm, mpi_err)
    end do
    ! broadcast the remainder
    offset = int(huge(0_i4), i8) * rounds + 1
    size = int(mod(cnt, huge(0_i4)), i4)
    call MPI_Bcast(bufptr(offset), size, sendtype, root, mpi_grp%comm, mpi_err)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_BCAST)
#endif
end subroutine X(mpi_grp_bcast_0_l)

! ---------------------------------------------------------
subroutine X(mpi_grp_allreduce)(mpi_grp, sendbuf, recvbuf, count, datatype, op)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: count, datatype, op

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLREDUCE)
  if (mpi_grp%comm /= -1) then
    call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
      mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, count)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLREDUCE)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, count)
#endif
end subroutine X(mpi_grp_allreduce)

! ---------------------------------------------------------
subroutine X(mpi_grp_allreduce_2)(mpi_grp, sendbuf, recvbuf, count, datatype, op)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:, :)
  R_TYPE,           intent(out) :: recvbuf(:, :)
  integer,          intent(in)  :: count, datatype, op

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLREDUCE)
  if (mpi_grp%comm /= -1) then
    call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
      mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_2)(sendbuf, recvbuf, count)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLREDUCE)
#else
  call X(mpi_grp_copy_2)(sendbuf, recvbuf, count)
#endif
end subroutine X(mpi_grp_allreduce_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_allreduce_3)(mpi_grp, sendbuf, recvbuf, count, datatype, op)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:, :, :)
  R_TYPE,           intent(out) :: recvbuf(:, :, :)
  integer,          intent(in)  :: count, datatype, op

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLREDUCE)
  if (mpi_grp%comm /= -1) then
    call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
      mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_3)(sendbuf, recvbuf, count)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLREDUCE)
#else
  call X(mpi_grp_copy_3)(sendbuf, recvbuf, count)
#endif
end subroutine X(mpi_grp_allreduce_3)

! ---------------------------------------------------------
subroutine X(mpi_grp_allreduce_0)(mpi_grp, sendbuf, recvbuf, count, datatype, op)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf
  R_TYPE,           intent(out) :: recvbuf
  integer,          intent(in)  :: count, datatype, op

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLREDUCE)
  if (mpi_grp%comm /= -1) then
    call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
      mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_0)(sendbuf, recvbuf, count)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLREDUCE)
#else
  call X(mpi_grp_copy_0)(sendbuf, recvbuf, count)
#endif
end subroutine X(mpi_grp_allreduce_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_allreduce_inplace_0)(mpi_grp, recvbuf, count, datatype, op)
  class(mpi_grp_t), intent(in)    :: mpi_grp
  R_TYPE,           intent(inout) :: recvbuf
  integer,          intent(in)    :: count, datatype, op

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLREDUCE)
  if (mpi_grp%comm /= -1) then
    call MPI_Allreduce(MPI_IN_PLACE, recvbuf, count, datatype, op, &
      mpi_grp%comm, mpi_err)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLREDUCE)
#endif
end subroutine X(mpi_grp_allreduce_inplace_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_allgather)(mpi_grp, sendbuf, sendcount, sendtype, recvbuf, &
  recvcount, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcount, sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount, recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHER)
  if (mpi_grp%comm /= -1) then
    call MPI_Allgather(sendbuf, sendcount, sendtype,&
      recvbuf, recvcount, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHER)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
#endif
end subroutine X(mpi_grp_allgather)

! ---------------------------------------------------------
subroutine X(mpi_grp_allgather_0)(mpi_grp, sendbuf, sendcount, sendtype, recvbuf, &
  recvcount, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf
  integer,          intent(in)  :: sendcount, sendtype
  R_TYPE,           intent(out) :: recvbuf
  integer,          intent(in)  :: recvcount, recvtype

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHER)
  if (mpi_grp%comm /= -1) then
    call MPI_Allgather(sendbuf, sendcount, sendtype,&
      recvbuf, recvcount, recvtype, mpi_grp%comm, mpi_err)
  else
    call X(mpi_grp_copy_0)(sendbuf, recvbuf, sendcount)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHER)
#else
  call X(mpi_grp_copy_0)(sendbuf, recvbuf, sendcount)
#endif
end subroutine X(mpi_grp_allgather_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_recv_0)(mpi_grp, recvbuf, recvcount, recvtype, source, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: recvbuf
  integer,           intent(in)  :: recvcount, recvtype
  integer,           intent(in)  :: source
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Recv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, MPI_STATUS_IGNORE, mpi_err)
#endif
end subroutine X(mpi_grp_recv_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_recv)(mpi_grp, recvbuf, recvcount, recvtype, source, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: recvbuf(:)
  integer,           intent(in)  :: recvcount, recvtype
  integer,           intent(in)  :: source
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Recv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, MPI_STATUS_IGNORE, mpi_err)
#endif
end subroutine X(mpi_grp_recv)

! ---------------------------------------------------------
subroutine X(mpi_grp_recv_2)(mpi_grp, recvbuf, recvcount, recvtype, source, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: recvbuf(:, :)
  integer,           intent(in)  :: recvcount, recvtype
  integer,           intent(in)  :: source
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Recv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, MPI_STATUS_IGNORE, mpi_err)
#endif
end subroutine X(mpi_grp_recv_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_recv_3)(mpi_grp, recvbuf, recvcount, recvtype, source, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: recvbuf(:, :, :)
  integer,           intent(in)  :: recvcount, recvtype
  integer,           intent(in)  :: source
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Recv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, MPI_STATUS_IGNORE, mpi_err)
#endif
end subroutine X(mpi_grp_recv_3)

! ---------------------------------------------------------
subroutine X(mpi_grp_send_0)(mpi_grp, sendbuf, sendcount, sendtype, dest, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: sendbuf
  integer,           intent(in)  :: sendcount, sendtype
  integer,           intent(in)  :: dest
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Send(sendbuf, sendcount, sendtype, dest, tag_, mpi_grp%comm, mpi_err)
#endif
end subroutine X(mpi_grp_send_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_send)(mpi_grp, sendbuf, sendcount, sendtype, dest, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: sendbuf(:)
  integer,           intent(in)  :: sendcount, sendtype
  integer,           intent(in)  :: dest
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Send(sendbuf, sendcount, sendtype, dest, tag_, mpi_grp%comm, mpi_err)
#endif
end subroutine X(mpi_grp_send)

! ---------------------------------------------------------
subroutine X(mpi_grp_send_2)(mpi_grp, sendbuf, sendcount, sendtype, dest, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: sendbuf(:, :)
  integer,           intent(in)  :: sendcount, sendtype
  integer,           intent(in)  :: dest
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Send(sendbuf, sendcount, sendtype, dest, tag_, mpi_grp%comm, mpi_err)
#endif
end subroutine X(mpi_grp_send_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_send_3)(mpi_grp, sendbuf, sendcount, sendtype, dest, tag)
  class(mpi_grp_t),  intent(in)  :: mpi_grp
  R_TYPE,            intent(out) :: sendbuf(:, :, :)
  integer,           intent(in)  :: sendcount, sendtype
  integer,           intent(in)  :: dest
  integer, optional, intent(in)  :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Send(sendbuf, sendcount, sendtype, dest, tag_, mpi_grp%comm, mpi_err)
#endif
end subroutine X(mpi_grp_send_3)

! ---------------------------------------------------------
subroutine X(mpi_grp_irecv_0_i8)(mpi_grp, recvbuf, recvcount, recvtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(out)   :: recvbuf
  integer(i8),       intent(in)    :: recvcount
  integer,           intent(in)    :: recvtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  ASSERT(recvcount < huge(0_i4))
  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Irecv(recvbuf, int(recvcount, i4), recvtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_irecv_0_i8)

! ---------------------------------------------------------
subroutine X(mpi_grp_irecv_0)(mpi_grp, recvbuf, recvcount, recvtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(out)   :: recvbuf
  integer,           intent(in)    :: recvcount, recvtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Irecv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_irecv_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_irecv)(mpi_grp, recvbuf, recvcount, recvtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(out)   :: recvbuf(:)
  integer,           intent(in)    :: recvcount, recvtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Irecv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_irecv)

! ---------------------------------------------------------
subroutine X(mpi_grp_irecv_2)(mpi_grp, recvbuf, recvcount, recvtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(out)   :: recvbuf(:, :)
  integer,           intent(in)    :: recvcount, recvtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Irecv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_irecv_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_irecv_3)(mpi_grp, recvbuf, recvcount, recvtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(out)   :: recvbuf(:, :, :)
  integer,           intent(in)    :: recvcount, recvtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
#if defined(HAVE_MPI)
  call MPI_Irecv(recvbuf, recvcount, recvtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_irecv_3)

! ---------------------------------------------------------
subroutine X(mpi_grp_isend_0_i8)(mpi_grp, sendbuf, sendcount, sendtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(in)    :: sendbuf
  integer(i8),       intent(in)    :: sendcount
  integer,           intent(in)    :: sendtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  ASSERT(sendcount < huge(0_i4))
  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Isend(sendbuf, int(sendcount, i4), sendtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_isend_0_i8)

! ---------------------------------------------------------
subroutine X(mpi_grp_isend_0)(mpi_grp, sendbuf, sendcount, sendtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(in)    :: sendbuf
  integer,           intent(in)    :: sendcount, sendtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Isend(sendbuf, sendcount, sendtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_isend_0)

! ---------------------------------------------------------
subroutine X(mpi_grp_isend)(mpi_grp, sendbuf, sendcount, sendtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(in)    :: sendbuf(:)
  integer,           intent(in)    :: sendcount, sendtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Isend(sendbuf, sendcount, sendtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_isend)

! ---------------------------------------------------------
subroutine X(mpi_grp_isend_2)(mpi_grp, sendbuf, sendcount, sendtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(in)    :: sendbuf(:, :)
  integer,           intent(in)    :: sendcount, sendtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Isend(sendbuf, sendcount, sendtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_isend_2)

! ---------------------------------------------------------
subroutine X(mpi_grp_isend_3)(mpi_grp, sendbuf, sendcount, sendtype, source, request, tag)
  class(mpi_grp_t),  intent(in)    :: mpi_grp
  R_TYPE,            intent(in)    :: sendbuf(:, :, :)
  integer,           intent(in)    :: sendcount, sendtype
  integer,           intent(in)    :: source
  integer,           intent(inout) :: request
  integer, optional, intent(in)    :: tag

  integer :: tag_

  tag_ = 0
  if (present(tag)) tag_ = tag
  if (mpi_grp%comm == -1) return
#if defined(HAVE_MPI)
  call MPI_Isend(sendbuf, sendcount, sendtype, source, tag_, mpi_grp%comm, request, mpi_err)
#endif
end subroutine X(mpi_grp_isend_3)

! Below are the interfaces for functions with large offsets.
! We need to implement them ourselves with send/recv patterns
! because the MPI 3.1 standard does not support them. Although
! the MPI 4.0 standard supports this, it is not yet implemented
! by MPI libraries as of today (beginning of 2022).

! ---------------------------------------------------------
subroutine X(mpi_grp_scatterv_i8)(mpi_grp, sendbuf, sendcnts, displs, sendtype, recvbuf, &
  recvcount, recvtype, root)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts(:)
  integer(i8),      intent(in)  :: displs(:)
  integer,          intent(in)  :: sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount
  integer,          intent(in)  :: recvtype, root

#if defined(HAVE_MPI)
  integer, parameter :: tag = 10
  integer :: request, irank, sendrequests(0:mpi_grp%size-1)
#endif

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_SCATTERV)
  if (mpi_grp%comm /= -1) then
    if (recvcount == 0) then
      request = MPI_REQUEST_NULL
    else
      call MPI_Irecv(recvbuf, recvcount, recvtype, root, tag, mpi_grp%comm, request, mpi_err)
    end if
    if (mpi_grp%rank == root) then
      do irank = 0, mpi_grp%size - 1
        if (sendcnts(irank+1) == 0) then
          sendrequests(irank) = MPI_REQUEST_NULL
        else
          call MPI_Isend(sendbuf(displs(irank+1)+1), sendcnts(irank+1), sendtype, irank, tag, mpi_grp%comm, &
            sendrequests(irank), mpi_err)
        end if
      end do
      call MPI_Waitall(mpi_grp%size, sendrequests, MPI_STATUSES_IGNORE, mpi_err)
    end if
    call MPI_Wait(request, MPI_STATUS_IGNORE, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_SCATTERV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount)
#endif
end subroutine X(mpi_grp_scatterv_i8)

! ---------------------------------------------------------
subroutine X(mpi_grp_gatherv_i8)(mpi_grp, sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype, root)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts, sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount(:)
  integer(i8),      intent(in)  :: displs(:)
  integer,          intent(in)  :: recvtype, root

#if defined(HAVE_MPI)
  integer, parameter :: tag = 11
  integer :: request, irank, recvrequests(0:mpi_grp%size-1)
#endif

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_GATHERV)
  if (mpi_grp%comm /= -1) then
    if (sendcnts == 0) then
      request = MPI_REQUEST_NULL
    else
      call MPI_Isend(sendbuf, sendcnts, sendtype, root, tag, mpi_grp%comm, request, mpi_err)
    end if
    if (mpi_grp%rank == root) then
      do irank = 0, mpi_grp%size - 1
        if (recvcount(irank+1) == 0) then
          recvrequests(irank) = MPI_REQUEST_NULL
        else
          call MPI_Irecv(recvbuf(displs(irank+1)+1), recvcount(irank+1), recvtype, irank, tag, mpi_grp%comm, &
            recvrequests(irank), mpi_err)
        end if
      end do
      call MPI_Waitall(mpi_grp%size, recvrequests, MPI_STATUSES_IGNORE, mpi_err)
    end if
    call MPI_Wait(request, MPI_STATUS_IGNORE, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_GATHERV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_gatherv_i8)

! ---------------------------------------------------------
subroutine X(mpi_grp_alltoallv_i8)(mpi_grp, sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
  recvcount, rdispls, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts(:)
  integer(i8),      intent(in)  :: sdispls(:)
  integer,          intent(in)  :: sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount(:)
  integer(i8),      intent(in)  :: rdispls(:)
  integer,          intent(in)  :: recvtype

#if defined(HAVE_MPI)
  integer, parameter :: tag = 11
  integer :: irank, irank_, recvrequests(0:mpi_grp%size-1), sendrequests(0:mpi_grp%size-1)
#endif

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLTOALLV)
  if (mpi_grp%comm /= -1) then
    ! post matching recvs and sends; first neighbors, then next-neighbors, and so on
    do irank = mpi_grp%rank, mpi_grp%rank + mpi_grp%size - 1
      if (irank > mpi_grp%size - 1) then
        irank_ = irank - mpi_grp%size
      else
        irank_ = irank
      end if
      if (recvcount(irank_+1) == 0) then
        recvrequests(irank_) = MPI_REQUEST_NULL
      else
        call MPI_Irecv(recvbuf(rdispls(irank_+1)+1), recvcount(irank_+1), recvtype, irank_, tag, &
          mpi_grp%comm, recvrequests(irank_), mpi_err)
      end if
    end do
    do irank = mpi_grp%rank, mpi_grp%rank - mpi_grp%size + 1, - 1
      if (irank < 0) then
        irank_ = irank + mpi_grp%size
      else
        irank_ = irank
      end if
      if (sendcnts(irank_+1) == 0) then
        sendrequests(irank_) = MPI_REQUEST_NULL
      else
        call MPI_Isend(sendbuf(sdispls(irank_+1)+1), sendcnts(irank_+1), sendtype, irank_, tag, &
          mpi_grp%comm, sendrequests(irank_), mpi_err)
      end if
    end do
    call MPI_Waitall(mpi_grp%size, recvrequests, MPI_STATUSES_IGNORE, mpi_err)
    call MPI_Waitall(mpi_grp%size, sendrequests, MPI_STATUSES_IGNORE, mpi_err)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
  end if
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLTOALLV)
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_alltoallv_i8)

! ---------------------------------------------------------
subroutine X(mpi_grp_allgatherv_i8)(mpi_grp, sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype)
  class(mpi_grp_t), intent(in)  :: mpi_grp
  R_TYPE,           intent(in)  :: sendbuf(:)
  integer,          intent(in)  :: sendcnts, sendtype
  R_TYPE,           intent(out) :: recvbuf(:)
  integer,          intent(in)  :: recvcount(:)
  integer(i8),      intent(in)  :: displs(:)
  integer,          intent(in)  :: recvtype

#if defined(HAVE_MPI)
  integer, parameter :: tag = 13
  integer :: irank, irank_, recvrequests(1:mpi_grp%size), sendrequests(1:mpi_grp%size)
  integer :: nrecvreqs, nsendreqs
#endif

#if defined(HAVE_MPI)
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHERV)
  if (mpi_grp%comm /= -1) then
    ASSERT(all(displs < sum(recvcount)))
    ASSERT(all(displs >= 0))
    ! post matching recvs and sends; first neighbors, then next-neighbors, and so on
    nrecvreqs = 0
    do irank = mpi_grp%rank, mpi_grp%rank + mpi_grp%size - 1
      if (irank > mpi_grp%size - 1) then
        irank_ = irank - mpi_grp%size
      else
        irank_ = irank
      end if
      if (recvcount(irank_+1) /= 0) then
        nrecvreqs = nrecvreqs + 1
        call MPI_Irecv(recvbuf(displs(irank_+1)+1_i8), recvcount(irank_+1), recvtype, irank_, tag, &
          mpi_grp%comm, recvrequests(nrecvreqs), mpi_err)
      end if
    end do
    nsendreqs = 0
    do irank = mpi_grp%rank, mpi_grp%rank - mpi_grp%size + 1, - 1
      if (irank < 0) then
        irank_ = irank + mpi_grp%size
      else
        irank_ = irank
      end if
      if (sendcnts /= 0) then
        nsendreqs = nsendreqs + 1
        call MPI_Isend(sendbuf, sendcnts, sendtype, irank_, tag, &
          mpi_grp%comm, sendrequests(nsendreqs), mpi_err)
      end if
    end do
    call MPI_Waitall(nrecvreqs, recvrequests, MPI_STATUSES_IGNORE, mpi_err)
    call MPI_Waitall(nsendreqs, sendrequests, MPI_STATUSES_IGNORE, mpi_err)
    call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHERV)
  else
    call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
  end if
#else
  call X(mpi_grp_copy_1)(sendbuf, recvbuf, recvcount(1))
#endif
end subroutine X(mpi_grp_allgatherv_i8)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
