!! Copyright (C) 2007 X. Andrade
!! Copyright (C) 2021 N. Tancogne-Dejean
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

! --------------------------------------------------------------------------
!> Returns f_out = H' f_in, where H' is perturbation Hamiltonian
!! Note that e^ikr phase is applied to f_in, then is removed afterward
subroutine X(perturbation_magnetic_apply)(this, namespace, space, gr, hm, ik, f_in, f_out, set_bc)
  class(perturbation_magnetic_t), intent(in)    :: this
  type(namespace_t),              intent(in)    :: namespace
  type(space_t),                  intent(in)    :: space
  type(grid_t),                   intent(in)    :: gr
  type(hamiltonian_elec_t),       intent(in)    :: hm
  integer,                        intent(in)    :: ik
  R_TYPE,                         intent(in)    :: f_in(:, :)
  R_TYPE,                         intent(out)   :: f_out(:, :)
  logical,        optional,       intent(in)    :: set_bc

  R_TYPE, allocatable :: f_in_copy(:, :), lf(:,:), vrnl(:,:,:)
  logical :: apply_kpoint, set_bc_
  integer :: ip, idim, iatom, idir
  type(profile_t), save :: prof
  R_TYPE :: cross(1:space%dim)

  PUSH_SUB(X(perturbation_magnetic_apply))

  call profiling_in(prof, TOSTRING(X(PERT_MAG_APPLY)))

  ASSERT(this%dir /= -1)

  set_bc_ = optional_default(set_bc, .true.)

  SAFE_ALLOCATE(f_in_copy(1:gr%np_part, 1:hm%d%dim))
  if (set_bc_) then
    call lalg_copy(gr%np, hm%d%dim, f_in, f_in_copy)
    do idim = 1, hm%d%dim
      call boundaries_set(gr%der%boundaries, gr, f_in_copy(:, idim))
    end do
  else
    call lalg_copy(gr%np_part, hm%d%dim, f_in, f_in_copy)
  end if

  apply_kpoint = allocated(hm%hm_base%phase)
  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_in_copy, hm%hm_base%phase(1:gr%np_part, ik), gr%np_part, .false.)
#endif
  end if

  SAFE_ALLOCATE(lf(1:gr%np, 1:space%dim))

  do idim = 1, hm%d%dim
    ! Note that we leave out the term 1/P_c
    call X(physics_op_L)(gr%der, f_in_copy(:, idim), lf, set_bc = .false.)
    f_out(1:gr%np, idim) = M_HALF*lf(1:gr%np, this%dir)
  end do

  SAFE_DEALLOCATE_A(lf)

  if (this%gauge == GAUGE_GIPAW .or. this%gauge == GAUGE_ICL) then
    SAFE_ALLOCATE(vrnl(1:gr%np, 1:hm%d%dim, 1:space%dim))

    do iatom = 1, this%ions%natoms

      vrnl = M_ZERO
      do idir = 1, space%dim
        if (this%dir == idir) cycle ! this direction is not used in the cross product
        call X(projector_commute_r)(hm%ep%proj(iatom), gr, gr%der%boundaries, hm%d%dim, idir, ik, f_in_copy, vrnl(:, :, idir))
      end do

      do idim = 1, hm%d%dim
        if (this%gauge == GAUGE_ICL) then
          do ip = 1, gr%np
#if !defined(R_TCOMPLEX)
            cross = dcross_product(gr%x(ip, 1:space%dim), vrnl(ip, idim, 1:space%dim))
            f_out(ip, idim) = f_out(ip, idim) + M_HALF*cross(this%dir)
#else
            cross = dzcross_product(gr%x(ip, 1:space%dim), vrnl(ip, idim, 1:space%dim))
            f_out(ip, idim) = f_out(ip, idim) - M_zI*M_HALF*cross(this%dir)
#endif
          end do
        else
          do ip = 1, gr%np
#if !defined(R_TCOMPLEX)
            cross = dcross_product(this%ions%pos(:,iatom), vrnl(ip, idim, 1:space%dim))
            f_out(ip, idim) = f_out(ip, idim) + M_HALF*cross(this%dir)
#else
            cross = dzcross_product(this%ions%pos(:,iatom), vrnl(ip, idim, 1:space%dim))
            f_out(ip, idim) = f_out(ip, idim) - M_zI*M_HALF*cross(this%dir)
#endif
          end do
        end if
      end do

    end do

    SAFE_DEALLOCATE_A(vrnl)
  end if

  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_out, hm%hm_base%phase(1:gr%np, ik), gr%np, .true.)
#endif
  end if

  SAFE_DEALLOCATE_A(f_in_copy)

  call profiling_out(prof)
  POP_SUB(X(perturbation_magnetic_apply))
end subroutine X(perturbation_magnetic_apply)

! --------------------------------------------------------------------------
subroutine X(perturbation_magnetic_apply_order_2) (this, namespace, space, gr, hm, ik, f_in, f_out)
  class(perturbation_magnetic_t), intent(in)    :: this
  type(namespace_t),              intent(in)    :: namespace
  type(space_t),                  intent(in)    :: space
  type(grid_t),                   intent(in)    :: gr
  type(hamiltonian_elec_t),       intent(in)    :: hm
  integer,                        intent(in)    :: ik
  R_TYPE,                         intent(in)    :: f_in(:, :)
  R_TYPE,                         intent(out)   :: f_out(:, :)

  integer :: ip, idim
  R_TYPE, allocatable :: f_in_copy(:,:)
  logical :: apply_kpoint
  R_TYPE, allocatable :: f_in2(:, :, :), dnl(:, :, :), vrnl(:,:), xf(:, :)
  FLOAT :: cross1(1:3), bdir(1:space%dim, 2)
  FLOAT :: rdelta, delta
  R_TYPE :: contr
  integer :: iatom, idir, idir2

  PUSH_SUB(X(perturbation_magnetic_apply_order_2))

  ASSERT(this%dir2 /= -1)

  SAFE_ALLOCATE(f_in_copy(1:gr%np, 1:hm%d%dim))
  call lalg_copy(gr%np, hm%d%dim, f_in, f_in_copy)

  apply_kpoint = allocated(hm%hm_base%phase)

  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_in_copy, hm%hm_base%phase(1:gr%np, ik), gr%np, .false.)
#endif
  end if

  delta = ddelta(this%dir, this%dir2)
  do idim = 1, hm%d%dim
    !$omp parallel do private(rdelta)
    do ip = 1, gr%np
      rdelta = sum(gr%x(ip, 1:space%dim)**2)*delta
      f_out(ip, idim) = M_FOURTH*(rdelta - gr%x(ip, this%dir)*gr%x(ip, this%dir2))*f_in_copy(ip, idim)
    end do
  end do

  ! gauge correction
  apply_gauge: if (this%gauge == GAUGE_GIPAW .or. this%gauge == GAUGE_ICL) then
    bdir(:, :) = M_ZERO
    bdir(this%dir,  1)   = M_ONE
    bdir(this%dir2, 2)   = M_ONE

    SAFE_ALLOCATE(f_in2(1:gr%np, 1:hm%d%dim, 1:space%dim))
    SAFE_ALLOCATE( vrnl(1:gr%np, 1:hm%d%dim))
    SAFE_ALLOCATE(  dnl(1:gr%np, 1:hm%d%dim, 1:space%dim))
    SAFE_ALLOCATE(   xf(1:gr%np, 1:hm%d%dim))

    f_in2 = R_TOTYPE(M_ZERO)
    atoms: do iatom = 1, this%ions%natoms

      ! This calculates f_in2 = (B x r) f_in_copy
      if(this%gauge == GAUGE_GIPAW) then
        cross1 = dcross_product(bdir(:, 2), this%ions%pos(:, iatom))
      end if

      do ip = 1, gr%np
        if(this%gauge == GAUGE_ICL) then
          cross1 = dcross_product(bdir(:, 2), gr%x(ip, :))
        end if

        do idim = 1,hm%d%dim
          f_in2(ip, idim, 1:space%dim) = cross1(1:space%dim)*f_in_copy(ip, idim)
        end do
      end do

      ! let us now get sum_beta Dnl f_in2
      dnl = R_TOTYPE(M_ZERO)
      do idir = 1, space%dim
        do idir2 = 1, space%dim
          vrnl = M_ZERO
          !calculate dnl |f_in2> = -[x,vnl] |f_in2>
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, gr%der%boundaries, hm%d%dim, idir2, ik, f_in2(:, :, idir2), vrnl)

          do idim = 1, hm%d%dim
            !$omp parallel do
            do ip = 1, gr%np
              ! -x vnl |f>
              dnl(ip, idim, idir) = dnl(ip, idim, idir) - gr%x(ip, idir)*vrnl(ip, idim)
              ! vnl x |f>
              xf(ip, idim) = gr%x(ip, idir) * f_in2(ip, idim, idir2)
            end do
          end do

          vrnl = M_ZERO
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, gr%der%boundaries, hm%d%dim, idir2, ik, xf, vrnl)

          call lalg_axpy(gr%np, hm%d%dim, M_ONE, vrnl, dnl(:,:,idir))
        end do
      end do

      if(this%gauge == GAUGE_GIPAW) then
        cross1 = dcross_product(bdir(:, 1), this%ions%pos(:, iatom))
      end if

      do ip = 1, gr%np
        if (this%gauge == GAUGE_ICL) then
          cross1 = dcross_product(bdir(:, 1), gr%x(ip, :))
        end if

        do idim = 1, hm%d%dim
          contr = M_ZERO
          do idir = 1, space%dim
            contr = contr + cross1(idir)*dnl(ip, idim, idir)
          end do
          f_out(ip, idim) = f_out(ip, idim) + M_FOURTH * contr
        end do
      end do
    end do atoms

    SAFE_DEALLOCATE_A(f_in2)
    SAFE_DEALLOCATE_A(vrnl)
    SAFE_DEALLOCATE_A(dnl)
    SAFE_DEALLOCATE_A(xf)
  end if apply_gauge

  if (apply_kpoint) then
#ifdef R_TCOMPLEX
    call states_elec_set_phase(hm%d, f_out, hm%hm_base%phase(1:gr%np, ik), gr%np, .true.)
#endif
  end if

  SAFE_DEALLOCATE_A(f_in_copy)

  POP_SUB(X(perturbation_magnetic_apply_order_2))
end subroutine X(perturbation_magnetic_apply_order_2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
