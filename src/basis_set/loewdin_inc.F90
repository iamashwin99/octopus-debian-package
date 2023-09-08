!! Copyright (C) 2018 N. Tancogne-Dejean
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

subroutine X(loewdin_orthogonalize)(basis, kpt, namespace)
  type(orbitalbasis_t), target, intent(inout) :: basis
  type(distributed_t),          intent(in)    :: kpt
  type(namespace_t),            intent(in)    :: namespace

  R_TYPE, allocatable :: overlap(:,:), overlap2(:,:)
  FLOAT,  allocatable :: eigenval(:)
  integer :: ik, is, np
  integer :: ind, ind2, ios, ios2, iorb, iorb2, ns, idim
  type(orbitalset_t), pointer :: os, os2
  type(profile_t), save :: prof

  PUSH_SUB(X(loewdin_orthogonalize))

  call profiling_in(prof, TOSTRING(X(LOEWDIN)))

  SAFE_ALLOCATE(overlap(1:basis%size,1:basis%size))
  SAFE_ALLOCATE(overlap2(1:basis%size,1:basis%size))
  SAFE_ALLOCATE(eigenval(1:basis%size))

  if (debug%info) then
    write(message(1), '(a)') 'Debug: Orthogonalizing the atomic orbital basis.'
    call messages_info(1, namespace=namespace)
  end if

  do ik = kpt%start, kpt%end

    call X(loewdin_overlap)(basis, overlap, ik)
    if (debug%info) call X(print_matrix)(basis, 'overlap', namespace, overlap, ik)

    !We compute S^{-1/2} from S
    call lalg_eigensolve(basis%size, overlap, eigenval)
    if (debug%info) call X(print_matrix)(basis, 'eigenvec', namespace, overlap, ik)
    overlap2 = R_CONJ(transpose(overlap))
    do is = 1, basis%size
      eigenval(is) = M_ONE/sqrt(eigenval(is))
      overlap2(is, 1:basis%size) = eigenval(is) * overlap2(is, 1:basis%size)
    end do
    overlap = matmul(overlap,overlap2)
    if (debug%info) call X(print_matrix)(basis, 'loewdin', namespace, overlap, ik)

    ! In the case of no phase, the eorb_mesh is allocated and used as a temporary array
    do ios = 1, basis%norbsets
      os => basis%orbsets(ios)
      if (.not. allocated(os%phase)) then
        if (os%use_submesh) then
          SAFE_ALLOCATE(os%eorb_mesh(1:os%sphere%np, 1:os%norbs, 1:os%ndim, 1:1))
        else
          SAFE_ALLOCATE(os%eorb_mesh(1:os%sphere%mesh%np, 1:os%norbs, 1:os%ndim, 1:1))
        end if
      end if
    end do

    ! We now contruct the orthogonalized basis
    do ind = 1, basis%size
      ios  = basis%global2os(1, ind)
      iorb = basis%global2os(2, ind)
      os => basis%orbsets(ios)
      os%eorb_mesh(:, iorb, :, ik) = R_TOTYPE(M_ZERO)
      do idim = 1, os%ndim
        do ind2 = 1, basis%size
          ios2  = basis%global2os(1, ind2)
          iorb2 = basis%global2os(2, ind2)
          os2 => basis%orbsets(ios2)
          ns = os2%sphere%np

          if (allocated(os%phase)) then
#ifdef R_TCOMPLEX
            ! if (abs(overlap(ind2, ind)) < CNST(1e-6)) cycle
            do is = 1, ns
              os%eorb_mesh(os2%sphere%map(is), iorb, idim, ik) &
                = os%eorb_mesh(os2%sphere%map(is), iorb, idim, ik) &
                + overlap(ind2, ind)*os2%zorb(is, idim, iorb2)*os2%phase(is, ik)
            end do
#else
            ASSERT(.false.)
#endif
          else
            if (.not. os%use_submesh) then
              !$omp parallel do
              do is = 1, os2%sphere%mesh%np
                os%eorb_mesh(is, iorb, idim, ik) &
                  = os%eorb_mesh(is, iorb, idim, ik) &
                  + overlap(ind2, ind)*os2%X(orb)(is, idim, iorb2)
              end do
            end if
          end if
        end do
      end do
    end do

    ! We now replace the old, nonorthogonalized basis by the new one
    do ind = 1, basis%size
      ios  = basis%global2os(1, ind)
      iorb = basis%global2os(2, ind)
      os => basis%orbsets(ios)
      if (.not. allocated(os%phase)) then
#ifdef R_TCOMPLEX
        os%X(orb)(:,:,iorb) = os%eorb_mesh(:,iorb,:,1)
#else
        os%X(orb)(:,:,iorb) = TOFLOAT(os%eorb_mesh(:,iorb,:,1))
#endif
      end if
    end do
    do ios = 1, basis%norbsets
      os => basis%orbsets(ios)
      if (.not. allocated(os%phase)) then
        SAFE_DEALLOCATE_A(os%eorb_mesh)
      end if
    end do

    ! If we use GPUs, we need to transfer the orbitals on the device
    if (accel_is_enabled() .and. os%ndim == 1) then
      do ios = 1, basis%norbsets
        os => basis%orbsets(ios)
        np = os%sphere%np
        if(.not. os%use_submesh) np = os%sphere%mesh%np

        if (.not. allocated(os%phase)) then
          do iorb = 1, os%norbs
            call accel_write_buffer(os%X(buff_orb), np, os%X(orb)(:, 1, iorb), offset = (iorb - 1)*os%ldorbs)
          end do
        else
          if(os%use_submesh) then
            do iorb = 1, os%norbs
              call accel_write_buffer(os%buff_eorb(ik), np, &
                os%eorb_submesh(:, 1, iorb, ik), offset = (iorb - 1)*os%ldorbs)
            end do
          else
            do iorb = 1, os%norbs
              call accel_write_buffer(os%buff_eorb(ik), np, &
                os%eorb_mesh(:, iorb, 1, ik), offset = (iorb - 1)*os%ldorbs)
            end do
          end if
        end if
      end do
    end if


    ! For debugging, we want to control what is the overlap matrix after
    ! orthogonalization
    if (debug%info) then
      call X(loewdin_overlap)(basis, overlap, ik)
      call X(print_matrix)(basis, 'overlap_after', namespace, overlap, ik)
    end if

  end do

  SAFE_DEALLOCATE_A(overlap)
  SAFE_DEALLOCATE_A(eigenval)

  if (debug%info) then
    write(message(1), '(a)') 'Debug: Orthogonalization completed.'
    call messages_info(1, namespace=namespace)
  end if

  call profiling_out(prof)

  POP_SUB(X(loewdin_orthogonalize))

end subroutine X(loewdin_orthogonalize)

subroutine X(loewdin_overlap)(basis, overlap, ik)
  type(orbitalbasis_t), target, intent(in)   :: basis
  R_TYPE,                       intent(inout):: overlap(:,:) !overlap matrix (norb, norb)
  integer,                      intent(in)   :: ik

  integer :: ios, ios2, iorb, iorb2
  integer :: ind, ind2
  integer :: idim
  type(orbitalset_t), pointer :: os, os2

  PUSH_SUB(X(loewdin_overlap))

  overlap(1:basis%size, 1:basis%size) = R_TOTYPE(M_ZERO)

  do ind = 1, basis%size
    ios  = basis%global2os(1, ind)
    iorb = basis%global2os(2, ind)
    os => basis%orbsets(ios)
    do ind2 = ind, basis%size
      ios2  = basis%global2os(1, ind2)
      iorb2 = basis%global2os(2, ind2)
      os2 => basis%orbsets(ios2)

      if (allocated(os%phase) .and. .not. os%use_submesh) then
#ifdef R_TCOMPLEX
        overlap(ind,ind2) = M_Z0
        do idim = 1, os%ndim
          overlap(ind, ind2) = overlap(ind, ind2) + zmf_dotp(os%sphere%mesh, &
            os%eorb_mesh(:, iorb, idim, ik), os2%eorb_mesh(:, iorb2, idim, ik))
        end do
#endif
      else
        if (os%use_submesh) then
          call messages_not_implemented("Lowdin orthogonalization with submeshes")
        else
          overlap(ind,ind2) = M_Z0
          do idim = 1, os%ndim
            overlap(ind, ind2) = overlap(ind, ind2) + X(mf_dotp)(os%sphere%mesh, &
              os%X(orb)(:, idim, iorb), os2%X(orb)(:, idim, iorb2))
          end do
        end if
      end if
    end do !ind2
  end do !ind

  !The overlap matrix is Hermitian
  do ind = 1, basis%size
    do ind2 = 1, ind-1
      overlap(ind, ind2) = R_CONJ(overlap(ind2, ind))
    end do
  end do


  POP_SUB(X(loewdin_overlap))
end subroutine X(loewdin_overlap)

subroutine X(loewdin_info)(basis, kpt, namespace)
  type(orbitalbasis_t),    intent(inout):: basis
  type(distributed_t),     intent(in)   :: kpt
  type(namespace_t),       intent(in)   :: namespace

  R_TYPE, allocatable :: overlap(:,:)
  integer :: ik

  PUSH_SUB(X(loewdin_info))

  SAFE_ALLOCATE(overlap(1:basis%size,1:basis%size))

  do ik = kpt%start, kpt%end
    call X(loewdin_overlap)(basis, overlap, ik)
    if (debug%info) call X(print_matrix)(basis, 'overlap', namespace, overlap, ik)
  end do

  SAFE_DEALLOCATE_A(overlap)

  POP_SUB(X(loewdin_info))
end subroutine X(loewdin_info)

subroutine X(print_matrix)(basis, label, namespace, overlap, ik)
  type(orbitalbasis_t),   intent(in) :: basis
  character(len=*),       intent(in) :: label
  type(namespace_t),      intent(in) :: namespace
  R_TYPE,                 intent(in) :: overlap(:,:) !< (basis%size,basis%size)
  integer,                intent(in) :: ik

  integer :: is, is2, iunit
  character(len=5) :: filename

  if (.not. mpi_grp_is_root(mpi_world)) return

  PUSH_SUB(X(print_matrix))

  write(filename, '(i5.5)') ik
  iunit = io_open(trim(basis%debugdir) // "/" //trim(label) // "_nk" // filename, namespace, &
    action='write')
  write(iunit,'(a)') ' Orthogonalization matrix '
#ifdef R_TCOMPLEX
  write(iunit,'(a)') ' Real part '
#endif
  do is = 1, basis%size
    do is2 = 1, basis%size-1
      write(iunit,'(es19.12,1x)',advance='no') TOFLOAT(overlap(is,is2))
    end do
    write(iunit,'(es19.12)') TOFLOAT(overlap(is,basis%size))
  end do
#ifdef R_TCOMPLEX
  write(iunit,'(a)') ' Imaginary part '
  do is = 1, basis%size
    do is2 = 1, basis%size-1
      write(iunit,'(es19.12,1x)',advance='no') aimag(overlap(is,is2))
    end do
    write(iunit,'(es19.12)') aimag(overlap(is,basis%size))
  end do
#endif
  call io_close(iunit)


  POP_SUB(X(print_matrix))
end subroutine
