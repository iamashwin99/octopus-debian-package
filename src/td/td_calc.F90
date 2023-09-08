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

#include "global.h"

module td_calc_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use ext_partner_list_oct_m
  use forces_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_base_oct_m
  use hamiltonian_elec_oct_m
  use interaction_partner_oct_m
  use ions_oct_m
  use iso_c_binding
  use kpoints_oct_m
  use lasers_oct_m
  use loct_math_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_calc_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m

  implicit none

  private
  public ::       &
    td_calc_tacc, &
    td_calc_tvel, &
    td_calc_ionch

contains

! ---------------------------------------------------------
!> Electronic acceleration (to calculate harmonic spectrum...)
!! It is calculated as:
!!
!! \f[
!! d2<x>/dt2 = d<p>/dt + i<[H,[V_nl,x]]> =
!!           = i<[V_l,p]> + i<[V_nl,p]> - E(t)N + i<[H,[V_nl,x]]>
!! \f]
!! \warning This subroutine only works if ions are not
!!          allowed to move
! ---------------------------------------------------------
  subroutine td_calc_tacc(namespace, space, gr, ions, ext_partners, st, hm, acc, time)
    type(namespace_t),        intent(in)  :: namespace
    type(space_t),            intent(in)  :: space
    type(grid_t),             intent(in)  :: gr
    type(ions_t),             intent(in)  :: ions
    type(partner_list_t),     intent(in)  :: ext_partners
    type(states_elec_t),      intent(in)  :: st
    type(hamiltonian_elec_t), intent(in)  :: hm
    FLOAT,                    intent(in)  :: time
    FLOAT,                    intent(out) :: acc(:)

    FLOAT :: field(space%dim), x(space%dim)
    CMPLX, allocatable :: zpsi(:, :), hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:)
    integer  :: j, k, ik, ist, idim
    type(lasers_t), pointer :: lasers

    PUSH_SUB(td_calc_tacc)

    ! The term i<[V_l,p]> + i<[V_nl,p]> may be considered as equal but opposite to the
    ! force exerted by the electrons on the ions. COMMENT: This has to be thought about.
    ! Maybe we are forgetting something....
    call total_force_calculate(space, gr, ions, hm%ep, st, hm%kpoints, acc, hm%lda_u_level)

    ! Adds the laser contribution : i<[V_laser, p]>
    ! WARNING: this ignores the possibility of non-electric td external fields.
    field = M_ZERO
    lasers => list_get_lasers(ext_partners)
    if(associated(lasers)) then
      do j = 1, lasers%no_lasers
        call laser_electric_field(lasers%lasers(j), field(1:space%dim), time, CNST(0.001))
        acc(1:space%dim) = acc(1:space%dim) - st%qtot*field(1:space%dim)
      end do
    end if

    if (.not. hm%ep%non_local) then
      POP_SUB(td_calc_tacc)
      return
    end if

    ! And now, i<[H,[V_nl,x]]>
    x = M_ZERO
    SAFE_ALLOCATE(zpsi(1:gr%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hzpsi (1:gr%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hhzpsi(1:3, 1:gr%np_part))

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end

        call states_elec_get_state(st, gr, ist, ik, zpsi)

        call zhamiltonian_elec_apply_single(hm, namespace, gr, zpsi, hzpsi, ist, ik)

        SAFE_ALLOCATE(xzpsi    (1:gr%np_part, 1:st%d%dim, 1:3))
        SAFE_ALLOCATE(vnl_xzpsi(1:gr%np_part, 1:st%d%dim))
        xzpsi = M_z0
        do k = 1, gr%np
          do j = 1, space%dim
            xzpsi(k, 1:st%d%dim, j) = gr%x(k, j)*zpsi(k, 1:st%d%dim)
          end do
        end do

        do j = 1, space%dim
          call zhamiltonian_elec_apply_single(hm, namespace, gr, xzpsi(:, :, j), vnl_xzpsi, ist, ik, &
            terms = TERM_NON_LOCAL_POTENTIAL)

          do idim = 1, st%d%dim
            x(j) = x(j) - 2*st%occ(ist, ik)*TOFLOAT(zmf_dotp(gr, hzpsi(1:gr%np, idim), vnl_xzpsi(:, idim)))
          end do
        end do

        xzpsi = M_z0
        do k = 1, gr%np
          do j = 1, space%dim
            xzpsi(k, 1:st%d%dim, j) = gr%x(k, j)*hzpsi(k, 1:st%d%dim)
          end do
        end do

        do j = 1, space%dim
          call zhamiltonian_elec_apply_single(hm, namespace, gr, xzpsi(:, :, j), vnl_xzpsi, ist, ik, &
            terms = TERM_NON_LOCAL_POTENTIAL)

          do idim = 1, st%d%dim
            x(j) = x(j) + 2*st%occ(ist, ik)*TOFLOAT(zmf_dotp(gr, zpsi(:, idim), vnl_xzpsi(:, idim)))
          end do
        end do
        SAFE_DEALLOCATE_A(xzpsi)
        SAFE_DEALLOCATE_A(vnl_xzpsi)

      end do
    end do
    SAFE_DEALLOCATE_A(hzpsi)
    SAFE_DEALLOCATE_A(hhzpsi)

    if (st%parallel_in_states) then
      call comm_allreduce(st%mpi_grp, x)
    end if
    acc = acc + x

    POP_SUB(td_calc_tacc)
  end subroutine td_calc_tacc

! ---------------------------------------------------------
!> Electronic velocity (to calculate harmonic spectrum...)
!! It is calculated as:
!! \f[
!! d<x>/dt = <p>
!! \f]
! ---------------------------------------------------------
  subroutine td_calc_tvel(space, der, st, kpoints, vel)
    type(space_t),       intent(in)  :: space
    type(derivatives_t), intent(in)  :: der
    type(states_elec_t), intent(in)  :: st
    type(kpoints_t),     intent(in)  :: kpoints
    FLOAT,               intent(out) :: vel(:)

    FLOAT, allocatable :: momentum(:,:,:)

    PUSH_SUB(td_calc_tvel)

    SAFE_ALLOCATE(momentum(1:space%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
    call states_elec_calc_momentum(st, space, der, kpoints, momentum)

    momentum(1:space%dim, st%st_start:st%st_end, 1) = &
      sum(momentum(1:space%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end), 3)
    momentum(1:space%dim, 1, 1) = sum(momentum(1:space%dim, st%st_start:st%st_end, 1), 2)
    vel = momentum(:, 1, 1)

    SAFE_DEALLOCATE_A(momentum)

    POP_SUB(td_calc_tvel)
  end subroutine td_calc_tvel

! ---------------------------------------------------------
!> Multiple ionization probabilities calculated form the KS orbital densities
!! C. Ullrich, Journal of Molecular Structure: THEOCHEM 501, 315 (2000).
! ---------------------------------------------------------
  subroutine td_calc_ionch(mesh, st, ch, Nch)
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
    integer,             intent(in)    :: Nch
    FLOAT,               intent(out)   :: ch(0:Nch)

    integer :: ik, ist, ii, jj, idim, Nid
    FLOAT   :: prod, prod0
    FLOAT, allocatable :: N(:), Nnot(:)
    CMPLX,  allocatable :: zpsi(:)
    !combinations
    integer :: next
    type(c_ptr) :: c
    integer, allocatable :: idx0(:), idx(:), idxref(:)

    PUSH_SUB(td_calc_ionch)

    SAFE_ALLOCATE(   N(1: Nch))
    SAFE_ALLOCATE(Nnot(1: Nch))
    SAFE_ALLOCATE(zpsi(1:mesh%np))

    N(:)   = M_ZERO
    Nnot(:)= M_ZERO


    ii = 1
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do idim = 1, st%d%dim

          if (st%st_start <= ist .and. ist <= st%st_end .and. &
            st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
            call states_elec_get_state(st, mesh, idim, ist, ik, zpsi)
            N(ii) = zmf_nrm2(mesh, zpsi)
            Nnot(ii) = M_ONE - N(ii)**2
          end if

          ii = ii + 1

        end do
      end do
    end do

    if (st%parallel_in_states) then
      call comm_allreduce(st%mpi_grp, N)
      call comm_allreduce(st%mpi_grp, Nnot)
    end if

! print * ,mpi_world%rank, "N    =", N(:)
! print * ,mpi_world%rank, "Nnot =", Nnot(:)

    ch(:) = M_ZERO

    Nid = Nch - 1
    SAFE_ALLOCATE(idx(0:Nid))
    SAFE_ALLOCATE(idxref(0:Nid))
    idxref = (/ (ii, ii = 0, Nid) /)
    do ii = 0, Nch
!     print *, "Nch= ", Nch, "ii", ii
      call loct_combination_init(c, Nch, ii)
      if (ii == 0) then
        SAFE_ALLOCATE(idx0(0:0))
      else
        SAFE_ALLOCATE(idx0(0:ii-1))
      end if
!     print *,"size(idx0,1)=",size(idx0,1)
      if (debug%info) then
        call messages_write("P(")
        call messages_write(ii)
        call messages_write(") = 0")
        call messages_new_line()
      end if
      ! Loop over combinations
      do
        prod  = M_ONE
        prod0 = M_ONE
        if (debug%info) then
          call messages_write("P(")
          call messages_write(ii)
          call messages_write(") += ")
        end if
        call loct_get_combination(c, idx0)
        idx(:) = idxref(:)
        do jj = 0, ii-1
          idx(idx0(jj))= -1
          if (debug%info) then
            call messages_write(" No(")
            call messages_write(idx0(jj)+1)
            call messages_write(") ")
          end if
          prod0 = prod0 * Nnot(idx0(jj)+1)
        end do

        do jj = 0, Nid
          if (idx(jj) >= 0) then
            if (debug%info) then
              call messages_write(" N(")
              call messages_write(idx(jj)+1)
              call messages_write(") ")
            end if
            prod = prod * N(idx(jj)+1)
          end if
        end do

        if (debug%info) call messages_new_line()

        ch(ii) = ch(ii) + prod*prod0

        call loct_combination_next(c, next)
        if (next /= 0) exit
      end do
      SAFE_DEALLOCATE_A(idx0)
      call loct_combination_end(c)

      if (debug%info) call messages_info()
    end do

    SAFE_DEALLOCATE_A(idx)
    SAFE_DEALLOCATE_A(idxref)


    SAFE_DEALLOCATE_A(N)
    SAFE_DEALLOCATE_A(Nnot)
    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(td_calc_ionch)
  end subroutine td_calc_ionch


end module td_calc_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
