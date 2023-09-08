!! Copyright (C) 2016 N. Tancogne-Dejean
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

module orbitalset_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use comm_oct_m
  use debug_oct_m
  use distributed_oct_m
  use global_oct_m
  use hardware_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use species_oct_m
  use submesh_oct_m
  use types_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                         &
    orbitalset_t,                   &
    orbitalset_init,                &
    orbitalset_end,                 &
    orbitalset_update_phase,        &
    orbitalset_update_phase_shift,  &
    dorbitalset_get_coefficients,   &
    zorbitalset_get_coefficients,   &
    dorbitalset_get_coeff_batch,    &
    zorbitalset_get_coeff_batch,    &
    dorbitalset_get_coeff_batch_accel, &
    zorbitalset_get_coeff_batch_accel, &
    dorbitalset_get_position_matrix_elem, &
    zorbitalset_get_position_matrix_elem, &
    dorbitalset_add_to_batch,       &
    zorbitalset_add_to_batch,       &
    orbitalset_set_jln

  type orbitalset_t
    ! Components are public by default
    integer             :: nn, ll, ii
    FLOAT               :: jj
    integer             :: norbs
    integer             :: ndim
    integer             :: iatom
    type(submesh_t)     :: sphere             !< The submesh of the orbital
    CMPLX, allocatable  :: phase(:,:)         !< Correction to the global phase
    !                                         !< if the sphere cross the border of the box
    FLOAT               :: Ueff               !< The effective U of the simplified rotational invariant form
    FLOAT               :: Ubar, Jbar
    FLOAT               :: alpha              !< A potential used to constrained occupations, as defined in PRB 71, 035105 (2005)
    integer             :: nneighbors         !< Number of neighbouring atoms on which the intersite
    !                                         !< interaction is considered
    FLOAT, allocatable  :: V_ij(:,:)          !< The list of intersite interaction parameters
    FLOAT, allocatable  :: coulomb_IIJJ(:,:,:,:,:) !< Coulomb integrales with neighboring atoms
    integer, allocatable:: map_os(:)
    CMPLX, allocatable  :: phase_shift(:,:)

    FLOAT               :: radius
    type(species_t), pointer :: spec
    integer             :: spec_index

    FLOAT, allocatable  :: dorb(:,:,:) !< The orbital, if real, on the submesh
    CMPLX, allocatable  :: zorb(:,:,:) !< The orbital, if complex, on the submesh
    CMPLX, allocatable  :: eorb_submesh(:,:,:,:) !< Orbitals with its phase factor, on the submesh (for isolated system with TD phase)
    CMPLX, allocatable  :: eorb_mesh(:,:,:,:) !< Orbitals with its phase factor, on the mesh (for periodic systems GS and TD)
    integer             :: ldorbs
    type(accel_mem_t)   :: dbuff_orb, zbuff_orb !< The accel buffers containing the orbitals
    type(accel_mem_t), allocatable   :: buff_eorb (:)

    logical             :: use_submesh        !< Do we use or not submeshes for the orbitals

    type(poisson_t)     :: poisson            !< For computing the Coulomb integrals
  end type orbitalset_t

contains

  subroutine orbitalset_init(this)
    type(orbitalset_t),             intent(inout) :: this

    PUSH_SUB(orbitalset_init)

    this%iatom = -1
    this%nneighbors = 0
    this%nn = 0
    this%ll = 0
    this%jj = M_ONE
    this%ii = 0
    this%iatom = 0
    this%ndim = 1
    this%spec_index = 0

    this%Ueff = M_ZERO
    this%Ubar = M_ZERO
    this%Jbar = M_ZERO
    this%alpha = M_ZERO
    this%radius = M_ZERO

    POP_SUB(orbitalset_init)
  end subroutine orbitalset_init


  subroutine orbitalset_end(this)
    type(orbitalset_t), intent(inout) :: this

    integer :: ik

    PUSH_SUB(orbitalset_end)

    SAFE_DEALLOCATE_A(this%phase)
    SAFE_DEALLOCATE_A(this%dorb)
    SAFE_DEALLOCATE_A(this%zorb)
    SAFE_DEALLOCATE_A(this%eorb_submesh)
    SAFE_DEALLOCATE_A(this%eorb_mesh)
    nullify(this%spec)
    call submesh_end(this%sphere)

    SAFE_DEALLOCATE_A(this%V_ij)
    SAFE_DEALLOCATE_A(this%coulomb_IIJJ)
    SAFE_DEALLOCATE_A(this%map_os)
    SAFE_DEALLOCATE_A(this%phase_shift)

    if(accel_is_enabled()) then
      call accel_release_buffer(this%dbuff_orb)
      call accel_release_buffer(this%zbuff_orb)
      if(allocated(this%buff_eorb)) then
        do ik = lbound(this%buff_eorb, dim=1), ubound(this%buff_eorb, dim=1)
          call accel_release_buffer(this%buff_eorb(ik))
        end do
        SAFE_DEALLOCATE_A(this%buff_eorb)
      end if
    end if

    POP_SUB(orbitalset_end)
  end subroutine orbitalset_end

  subroutine orbitalset_set_jln(this, jj, ll, nn)
    type(orbitalset_t), intent(inout) :: this
    FLOAT,              intent(in)    :: jj
    integer,            intent(in)    :: ll, nn

    PUSH_SUB(orbitalset_set_jln)

    this%jj = jj
    this%ll = ll
    this%nn = nn

    POP_SUB(orbitalset_set_jln)
  end subroutine orbitalset_set_jln


  !> Build the phase correction to the global phase in case the orbital crosses the border of the simulaton box
  subroutine orbitalset_update_phase(os, dim, kpt, kpoints, spin_polarized, vec_pot, vec_pot_var, kpt_max)
    type(orbitalset_t),            intent(inout) :: os
    integer,                       intent(in)    :: dim
    type(distributed_t),           intent(in)    :: kpt
    type(kpoints_t),               intent(in)    :: kpoints
    logical,                       intent(in)    :: spin_polarized
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot(:) !< (sb%dim)
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot_var(:, :) !< (1:sb%dim, 1:ns)
    integer, optional,             intent(in)    :: kpt_max

    integer :: ns, iq, is, ikpoint, im, idim, kpt_end
    FLOAT   :: kr, kpoint(1:dim), dx(1:dim)
    integer :: iorb

    PUSH_SUB(orbitalset_update_phase)

    ns = os%sphere%np

    kpt_end = kpt%end
    if (present(kpt_max)) kpt_end = min(kpt_max, kpt_end)

    do iq = kpt%start, kpt_end
      !This is durty but avoids to refer to states_get_kpoint_index
      if (spin_polarized) then
        ikpoint = 1 + (iq - 1)/2
      else
        ikpoint = iq
      end if

      ! if this fails, it probably means that sb is not compatible with std
      ASSERT(ikpoint <= kpoints_number(kpoints))

      kpoint(1:dim) = kpoints%get_point(ikpoint)

      do is = 1, ns
        ! this is only the correction to the global phase, that can
        ! appear if the sphere crossed the boundary of the cell.
        dx(1:dim) = os%sphere%rel_x(1:dim, is) - os%sphere%mesh%x(os%sphere%map(is), 1:dim) + os%sphere%center(1:dim)
        kr = sum(kpoint(1:dim)*dx(1:dim))
        if (present(vec_pot)) then
          if (allocated(vec_pot)) kr = kr + sum(vec_pot(1:dim)*dx(1:dim))
        end if

        if (present(vec_pot_var)) then
          if (allocated(vec_pot_var)) kr = kr + sum(vec_pot_var(1:dim, os%sphere%map(is)) &
            *(os%sphere%rel_x(1:dim, is)+os%sphere%center))
        end if

        os%phase(is, iq) = exp(M_zI*kr)
      end do

      if (.not. os%use_submesh) then
        !We now compute the so-called Bloch sum of the localized orbitals
        os%eorb_mesh(:,:,:,iq) = M_Z0
        do idim = 1, os%ndim
          do im = 1, os%norbs
            do is = 1, ns
              os%eorb_mesh(os%sphere%map(is),im,idim,iq) = os%eorb_mesh(os%sphere%map(is),im,idim,iq) &
                + os%zorb(is,idim,im)*os%phase(is, iq)
            end do
          end do
        end do
      else !In the case of the isolated system, we still use the submesh
        do im = 1, os%norbs
          do idim = 1, os%ndim
            do is = 1, ns
              os%eorb_submesh(is,idim,im,iq) = os%zorb(is,idim,im)*os%phase(is, iq)
            end do
          end do
        end do
      end if

      if(accel_is_enabled() .and. os%ndim == 1) then
        if(os%use_submesh) then
          do iorb = 1, os%norbs
            call accel_write_buffer(os%buff_eorb(iq), ns, os%eorb_submesh(:, 1, iorb, iq), offset = (iorb - 1)*os%ldorbs)
          end do
        else
          do iorb = 1, os%norbs
            call accel_write_buffer(os%buff_eorb(iq), os%sphere%mesh%np, &
              os%eorb_mesh(:, iorb, 1, iq), offset = (iorb - 1)*os%ldorbs)
          end do
        end if
      end if
    end do

    POP_SUB(orbitalset_update_phase)
  end subroutine orbitalset_update_phase


  !> Build the phase shift for the intersite interaction
  subroutine orbitalset_update_phase_shift(os, dim, kpt, kpoints, spin_polarized, vec_pot, vec_pot_var, kpt_max)
    type(orbitalset_t),            intent(inout) :: os
    integer,                       intent(in)    :: dim
    type(distributed_t),           intent(in)    :: kpt
    type(kpoints_t),               intent(in)    :: kpoints
    logical,                       intent(in)    :: spin_polarized
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot(:) !< (sb%dim)
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot_var(:, :) !< (1:sb%dim, 1:ns)
    integer, optional,             intent(in)    :: kpt_max

    integer :: iq, ikpoint
    FLOAT   :: kr, kpoint(dim), dx(dim)
    integer :: inn, kpt_end

    PUSH_SUB(orbitalset_update_phase_shift)

    kpt_end = kpt%end
    if(present(kpt_max)) kpt_end = min(kpt_max, kpt_end)

    do iq = kpt%start, kpt_end
      !This is durty but avoids to refer to states_get_kpoint_index
      if(spin_polarized) then
        ikpoint = 1 + (iq - 1)/2
      else
        ikpoint = iq
      end if

      ! if this fails, it probably means that sb is not compatible with std
      ASSERT(ikpoint <= kpoints_number(kpoints))

      kpoint(1:dim) = kpoints%get_point(ikpoint)

      if (os%nneighbors > 0) then
        do inn = 1, os%nneighbors
          dx(1:dim) = os%V_ij(inn,1:dim)
          kr = sum(kpoint(1:dim)*dx(1:dim))
          if (present(vec_pot)) then
            if (allocated(vec_pot)) kr = kr + sum(vec_pot(1:dim)*dx(1:dim))
          end if

          !At the moment the uniform vector potential is in vec_pot_var
          if (present(vec_pot_var)) then
            if (allocated(vec_pot_var))  kr = kr + sum(vec_pot_var(1:dim, 1)*dx(1:dim))
          end if

          !The sign is different as this is applied on the wavefunction and not the orbitals
          os%phase_shift(inn, iq) = exp(-M_zI*kr)
        end do
      end if
    end do


    POP_SUB(orbitalset_update_phase_shift)
  end subroutine orbitalset_update_phase_shift

#include "undef.F90"
#include "real.F90"
#include "orbitalset_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "orbitalset_inc.F90"

end module orbitalset_oct_m
