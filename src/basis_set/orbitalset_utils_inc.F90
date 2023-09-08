!! Copyright (C) 2017 N. Tancogne-Dejean
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


 !At the present time this routine can only return atomic orbitals, but this could be generalized
subroutine X(orbitalset_utils_getorbitals)(namespace, os, ions, mesh, use_mesh, normalize)
  type(namespace_t),       intent(in)    :: namespace
  type(orbitalset_t),      intent(inout) :: os
  type(ions_t),            intent(in)    :: ions
  type(mesh_t),            intent(in)    :: mesh
  logical,                 intent(in)    :: use_mesh
  logical,                 intent(in)    :: normalize

  integer :: iorb

  PUSH_SUB(X(orbitalset_utils_getorbitals))

  do iorb = 1, os%norbs
    if (debug%info) then
      write(message(1),'(a,i3,1x,i1,1x,i1,1x,i1,1x,f3.1)')  'get_atomic_orbital ', os%iatom, &
        iorb, os%ii, os%ll, os%jj
      call messages_info(1, namespace=namespace)
    end if
    ! We obtain the orbital
    call X(get_atomic_orbital)(namespace, ions%space, ions%latt, ions%pos(:,os%iatom), ions%atom(os%iatom)%species, mesh, &
      os%sphere, os%ii, os%ll, os%jj, os, iorb, os%radius, os%ndim, use_mesh, normalize)
  end do !iorb

  POP_SUB(X(orbitalset_utils_getorbitals))

end subroutine X(orbitalset_utils_getorbitals)


! ----------------------------------------------------------------------------
! For periodic systems we employ the method described here
! https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
subroutine X(orbitalset_get_center_of_mass)(os, space, mesh, latt)
  type(orbitalset_t),      intent(inout) :: os
  type(space_t),           intent(in)    :: space
  type(mesh_t),            intent(in)    :: mesh
  type(lattice_vectors_t), intent(in)    :: latt

  integer :: ip, idim, iorb, idir
  FLOAT :: den, rr, xx_red(space%dim), angle
  FLOAT :: av_eta(space%periodic_dim), av_chi(space%periodic_dim)
  FLOAT, allocatable :: xx(:,:), eta(:,:), chi(:,:)

  PUSH_SUB(X(orbitalset_get_center_of_mass))

  SAFE_ALLOCATE(xx(1:space%dim, 1:mesh%np))

  if(space%is_periodic()) then
    SAFE_ALLOCATE(eta(1:space%periodic_dim, 1:mesh%np))
    SAFE_ALLOCATE(chi(1:space%periodic_dim, 1:mesh%np))
  end if
  if (mesh%use_curvilinear) then
    do ip = 1, mesh%np
      call mesh_r(mesh, ip, rr, coords=xx(:,ip))
      xx_red = latt%cart_to_red(xx(:,ip))
      ! Working with reduced angle for periodic directions
      do idir = 1, space%periodic_dim
        angle = xx_red(idir)*M_TWO*M_PI + M_PI
        eta(idir, ip) = cos(angle) * mesh%vol_pp(ip)
        chi(idir, ip) = sin(angle) * mesh%vol_pp(ip)
      end do
      xx(space%periodic_dim+1:space%dim, ip) = xx(space%periodic_dim+1:space%dim, ip) * mesh%vol_pp(ip)
    end do
  else
    do ip = 1, mesh%np
      call mesh_r(mesh, ip, rr, coords=xx(:,ip))
      xx_red = latt%cart_to_red(xx(:,ip))
      ! Working with reduced angle for periodic directions
      do idir = 1, space%periodic_dim
        angle = xx_red(idir)*M_TWO*M_PI + M_PI
        eta(idir, ip) = cos(angle)
        chi(idir, ip) = sin(angle)
      end do
    end do
  end if

  SAFE_ALLOCATE(os%sphere%center(1:space%dim))
  os%sphere%center = M_ZERO
  av_eta = M_ZERO
  av_chi = M_ZERO
  do iorb = 1, os%norbs
    do idim = 1, os%ndim
      do ip = 1, mesh%np
        den = TOFLOAT(R_CONJ(os%X(orb)(ip, idim, iorb))*os%X(orb)(ip, idim, iorb))
        ! Periodic directions
        av_eta = av_eta + eta(:,ip) * den
        av_chi = av_chi + chi(:,ip) * den
        ! Non-periodic directions
        os%sphere%center(space%periodic_dim+1:space%dim) = os%sphere%center(space%periodic_dim+1:space%dim) &
          + xx(space%periodic_dim+1:space%dim,ip)*den
      end do
    end do
  end do

  if (mesh%parallel_in_domains) then
    call mesh%allreduce(av_eta)
    call mesh%allreduce(av_chi)
    call mesh%allreduce(os%sphere%center)
  end if

  os%sphere%center = os%sphere%center * mesh%volume_element / os%norbs
  av_eta = av_eta * mesh%volume_element / os%norbs
  av_chi = av_chi * mesh%volume_element / os%norbs

  ! Returning to absolute coordinates
  do idir = 1, space%periodic_dim
    angle = atan2(-av_chi(idir), -av_eta(idir))
    os%sphere%center(idir) = angle/(M_TWO*M_PI)
  end do
  xx_red = latt%red_to_cart(os%sphere%center)
  os%sphere%center(1:space%periodic_dim) = xx_red(1:space%periodic_dim)

  ! We draw a sphere within the unit cell
  ! A larger sphere would contain periodic copies
  ! Note that here we assume the orbitals to be spherical, which might not be a good approximation
  os%radius = minval(norm2(latt%rlattice(1:space%dim, 1:space%periodic_dim), dim=1))*M_HALF

  if(debug%info .and. space%dim == 3) then
    write(message(1), '(a,3(f6.3,a))') 'Info : Center of mass of the orbital at (', &
      os%sphere%center(1), ', ', os%sphere%center(2), ', ', os%sphere%center(3), ')'
    write(message(2), '(a,f6.3,a)') 'Info : Radius of the orbital is ', os%radius, ' [B]'
    call messages_info(2)
  end if

  SAFE_DEALLOCATE_A(xx)
  SAFE_DEALLOCATE_A(eta)
  SAFE_DEALLOCATE_A(chi)

  POP_SUB(X(orbitalset_get_center_of_mass))
end subroutine X(orbitalset_get_center_of_mass)

