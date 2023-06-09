!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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


!------------------------------------------------------------------------------
!> X(kb_project) calculates the action of the projector kb_p on the psi
!! wavefunction. The action of the projector kb_p is defined as:
!! \f[
!! \hat{kb_p} |psi> = \sum_{i}^kb_p\%n_c p\%e(i) |kb_p\%p(:, i)><kb_p\%p(:, i)|psi>
!! \f]
!! The result is summed up to ppsi.
!------------------------------------------------------------------------------
subroutine X(kb_project)(mesh, sm, kb_p, dim, psi, ppsi)
  type(mesh_t),         intent(in)    :: mesh
  type(submesh_t),      intent(in)    :: sm
  type(kb_projector_t), intent(in)    :: kb_p
  integer,              intent(in)    :: dim
  R_TYPE,               intent(in)    :: psi(:, :)  !< (kb%n_s, dim)
  R_TYPE,               intent(inout) :: ppsi(:, :) !< (kb%n_s, dim)

  R_TYPE, allocatable :: uvpsi(:,:)

  PUSH_SUB(X(kb_project))

  SAFE_ALLOCATE(uvpsi(1:dim, 1:kb_p%n_c))
  call X(kb_project_bra)(mesh, sm, kb_p, dim, psi, uvpsi)

  if (mesh%parallel_in_domains) call mesh%allreduce(uvpsi)

  call X(kb_project_ket)(kb_p, dim, uvpsi, ppsi)
  SAFE_DEALLOCATE_A(uvpsi)

  POP_SUB(X(kb_project))

end subroutine X(kb_project)

!--------------------------------------------------------------
!> THREADSAFE
subroutine X(kb_project_bra)(mesh, sm, kb_p, dim, psi, uvpsi)
  type(mesh_t),         intent(in)  :: mesh
  type(submesh_t),      intent(in)  :: sm
  type(kb_projector_t), intent(in)  :: kb_p
  integer,              intent(in)  :: dim
  R_TYPE,               intent(in)  :: psi(:,:)   !< (1:ns, 1:dim)
  R_TYPE,               intent(out) :: uvpsi(:,:) !< (1:dim, 1:kb_p%n_c)

  integer :: ic, idim, ns, is

  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(KB_PROJECT_BRA)))


#ifndef HAVE_OPENMP
  PUSH_SUB(X(kb_project_bra))
#endif

  ns = kb_p%n_s

  uvpsi(1:dim, 1:kb_p%n_c) = M_ZERO

  if (mesh%use_curvilinear) then

    do idim = 1, dim
      do ic = 1, kb_p%n_c
        do is = 1, ns
          uvpsi(idim, ic) = uvpsi(idim, ic) + (kb_p%p(is, ic)*psi(is, idim))*mesh%vol_pp(sm%map(is))
        end do
      end do
    end do

  else

    do idim = 1, dim
      do ic = 1, kb_p%n_c
        do is = 1, ns
          uvpsi(idim, ic) = uvpsi(idim, ic) + psi(is, idim)*kb_p%p(is, ic)
        end do
      end do
    end do

  end if

  uvpsi(1:dim, 1:kb_p%n_c) = uvpsi(1:dim, 1:kb_p%n_c)*mesh%volume_element

#ifndef HAVE_OPENMP
  POP_SUB(X(kb_project_bra))
#endif

  call profiling_out(prof)

end subroutine X(kb_project_bra)

!--------------------------------------------------------------
!> THREADSAFE
subroutine X(kb_project_ket)(kb_p, dim, uvpsi, psi)
  type(kb_projector_t), intent(in)    :: kb_p
  integer,              intent(in)    :: dim
  R_TYPE,               intent(inout) :: uvpsi(:,:) !< (1:dim, 1:kb_p%n_c)
  R_TYPE,               intent(inout) :: psi(:, :) !< (1:ns, 1:dim)

  integer :: ic, idim, ns, is
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(KB_PROJECT_KET)))


#ifndef HAVE_OPENMP
  PUSH_SUB(X(kb_project_ket))
#endif

  ns = kb_p%n_s

  call X(kb_mul_energies)(kb_p, dim, uvpsi)

  do idim = 1, dim
    do ic = 1, kb_p%n_c
      do is = 1, ns
        psi(is, idim) = psi(is, idim) + uvpsi(idim, ic)*kb_p%p(is, ic)
      end do
    end do
  end do

#ifndef HAVE_OPENMP
  POP_SUB(X(kb_project_ket))
#endif

  call profiling_out(prof)

end subroutine X(kb_project_ket)

!--------------------------------------------------------------
!> THREADSAFE
subroutine X(kb_mul_energies)(kb_p, dim, uvpsi)
  type(kb_projector_t), intent(in)    :: kb_p
  integer,              intent(in)    :: dim
  R_TYPE,               intent(inout) :: uvpsi(:,:) !< (1:dim, 1:kb_p%n_c)

  integer :: idim

  do idim = 1, dim
    uvpsi(idim, 1:kb_p%n_c) = uvpsi(idim, 1:kb_p%n_c)*kb_p%e(1:kb_p%n_c)
  end do
end subroutine X(kb_mul_energies)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
