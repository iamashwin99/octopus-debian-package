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
!> This routine includes occupations for psib if perturbation_order == 2, correct if used
!! as \f$ <\psi(0)|H(2)|\psi(0)> \f$. It does not include occupations if perturbation_order == 1,
!! correct if used as <psi(0)|H(1)|psi(1)> since the LR wavefunctions include the
!! occupation. This routine must be modified if used differently than these two
!! ways.
subroutine X(perturbation_expectation_density) (this, namespace, space, gr, hm, st, psia, psib, density, perturbation_order)
  class(perturbation_t),    intent(in)    :: this
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(grid_t),             intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  R_TYPE,                   intent(in)    :: psia(:, :, :, st%d%kpt%start:)
  R_TYPE,                   intent(in)    :: psib(:, :, :, st%d%kpt%start:)
  R_TYPE,                   intent(out)   :: density(:)
  integer, optional,        intent(in)    :: perturbation_order

  R_TYPE, allocatable :: pertpsib(:, :)
  integer :: ik, ist, idim, order
  FLOAT   :: ikweight

  PUSH_SUB(X(perturbation_expectation_density))

  SAFE_ALLOCATE(pertpsib(1:gr%np, 1:st%d%dim))

  order = optional_default(perturbation_order, 1)
  ASSERT(order == 1 .or. order == 2)

  density(1:gr%np) = R_TOTYPE(M_ZERO)

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      if (order == 1) then
        call this%X(apply)(namespace, space, gr, hm, ik, psib(:, :, ist, ik), pertpsib)
        ikweight = st%d%kweights(ik) * st%smear%el_per_state
      else
        call this%X(apply_order_2)(namespace, space, gr, hm, ik, psib(:, :, ist, ik), pertpsib)
        ikweight = st%d%kweights(ik) * st%occ(ist, ik)
      end if

      do idim = 1, st%d%dim
        density(1:gr%np) = density(1:gr%np) + ikweight * &
          R_CONJ(psia(1:gr%np, idim, ist, ik)) * pertpsib(1:gr%np, idim)
      end do

    end do
  end do

  SAFE_DEALLOCATE_A(pertpsib)
  POP_SUB(X(perturbation_expectation_density))

end subroutine X(perturbation_expectation_density)

! --------------------------------------------------------------------------
R_TYPE function X(perturbation_expectation_value) (this, namespace, space, gr, hm, st, psia, psib, &
  perturbation_order) result(expval)
  class(perturbation_t),    intent(in)    :: this
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(grid_t),             intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  R_TYPE,                   intent(in)    :: psia(:, :, :, st%d%kpt%start:)
  R_TYPE,                   intent(in)    :: psib(:, :, :, st%d%kpt%start:)
  integer, optional,        intent(in)    :: perturbation_order

  R_TYPE, allocatable :: density(:)
  integer :: order

  PUSH_SUB(X(perturbation_expectation_value))

  order = optional_default(perturbation_order, 1)

  ASSERT(order == 1 .or. order == 2)

  SAFE_ALLOCATE(density(1:gr%np))

  call X(perturbation_expectation_density)(this, namespace, space, gr, hm, st, psia, psib, density, perturbation_order = order)

  expval = X(mf_integrate)(gr, density)
  SAFE_DEALLOCATE_A(density)

  if (st%parallel_in_states .or. st%d%kpt%parallel) then
    call comm_allreduce(st%st_kpt_mpi_grp, expval)
  end if

  POP_SUB(X(perturbation_expectation_value))

end function X(perturbation_expectation_value)

! --------------------------------------------------------------------------

R_TYPE function X(perturbation_states_elec_expectation_value)(this, namespace, space, gr, hm, st, &
  perturbation_order) result(expval)
  class(perturbation_t),    intent(in)    :: this
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(grid_t),             intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer, optional,        intent(in)    :: perturbation_order

  integer :: order, ik, ib, minst, maxst, ist
  R_TYPE, allocatable :: tt(:)
  type(wfs_elec_t) :: hpsib

  PUSH_SUB(X(perturbation_states_elec_expectation_value))

  order = optional_default(perturbation_order, 1)

  ASSERT(order == 1)

  SAFE_ALLOCATE(tt(st%st_start:st%st_end))

  expval = M_ZERO
  do ik = st%d%kpt%start, st%d%kpt%end
    tt = M_ZERO

    do ib = st%group%block_start, st%group%block_end
      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)

      call st%group%psib(ib, ik)%copy_to(hpsib)

      call this%apply_batch(namespace, space, gr, hm, st%group%psib(ib, ik), hpsib)
      call X(mesh_batch_dotp_vector)(gr, st%group%psib(ib, ik), hpsib, tt(minst:maxst))

      call hpsib%end(copy = .false.)

    end do

    do ist = st%st_start, st%st_end
      expval = expval + st%d%kweights(ik)*st%smear%el_per_state*tt(ist)
    end do

  end do

  if (st%parallel_in_states .or. st%d%kpt%parallel) then
    call comm_allreduce(st%st_kpt_mpi_grp, expval)
  end if


  POP_SUB(X(perturbation_states_elec_expectation_value))

end function X(perturbation_states_elec_expectation_value)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
