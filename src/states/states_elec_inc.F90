!! Copyright (C) 2011 X. Andrade
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

subroutine X(states_elec_get_state2)(st, mesh, ist, iqn, psi)
  type(states_elec_t), intent(in)    :: st
  class(mesh_t),       intent(in)    :: mesh
  integer,             intent(in)    :: ist       !< current state
  integer,             intent(in)    :: iqn       !< current k-point
  R_TYPE,              intent(inout) :: psi(:, :)

  integer :: idim

  PUSH_SUB(X(states_elec_get_state2))

  do idim =  1, st%d%dim
    call X(states_elec_get_state1)(st, mesh, idim, ist, iqn, psi(:, idim))
  end do

  POP_SUB(X(states_elec_get_state2))
end subroutine X(states_elec_get_state2)

! ------------------------------------------------------------

subroutine X(states_elec_get_state1)(st, mesh, idim, ist, iqn, psi)
  type(states_elec_t), intent(in)    :: st
  class(mesh_t),       intent(in)    :: mesh
  integer,             intent(in)    :: idim   !< current dimension
  integer,             intent(in)    :: ist
  integer,             intent(in)    :: iqn    !< current k-point
  R_TYPE,              intent(inout) :: psi(:)

  PUSH_SUB(X(states_elec_get_state1))

  call batch_get_state(st%group%psib(st%group%iblock(ist, iqn), iqn), (/ist, idim/), mesh%np, psi)

  POP_SUB(X(states_elec_get_state1))
end subroutine X(states_elec_get_state1)

! ------------------------------------------------------------

subroutine X(states_elec_get_state4)(st, mesh, psi)
  type(states_elec_t), intent(in)    :: st
  class(mesh_t),       intent(in)    :: mesh
  R_TYPE,              intent(inout) :: psi(:, :, st%st_start:, st%d%kpt%start:)

  integer :: iqn, ist

  PUSH_SUB(X(states_elec_get_state4))

  do iqn = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      call states_elec_get_state(st, mesh, ist, iqn, psi(:, :, ist, iqn))
    end do
  end do

  POP_SUB(X(states_elec_get_state4))
end subroutine X(states_elec_get_state4)

! ------------------------------------------------------------

subroutine X(states_elec_get_state3)(st, mesh, iqn, psi)
  type(states_elec_t), intent(in)    :: st
  class(mesh_t),       intent(in)    :: mesh
  integer,             intent(in)    :: iqn
  R_TYPE,              intent(inout) :: psi(:, :, st%st_start:)

  integer :: ist

  PUSH_SUB(X(states_elec_get_state3))

  do ist = st%st_start, st%st_end
    call states_elec_get_state(st, mesh, ist, iqn, psi(:, :, ist))
  end do

  POP_SUB(X(states_elec_get_state3))
end subroutine X(states_elec_get_state3)

! ------------------------------------------------------------

subroutine X(states_elec_set_state2)(st, mesh, ist, iqn, psi)
  type(states_elec_t), intent(inout) :: st
  class(mesh_t),       intent(in)    :: mesh
  integer,             intent(in)    :: ist       !< current dimension
  integer,             intent(in)    :: iqn       !< current k-point
  R_TYPE,              intent(in)    :: psi(:, :)

  integer :: idim

  PUSH_SUB(X(states_elec_set_state2))

  do idim =  1, st%d%dim
    call X(states_elec_set_state1)(st, mesh, idim, ist, iqn, psi(:, idim))
  end do

  POP_SUB(X(states_elec_set_state2))
end subroutine X(states_elec_set_state2)

! ------------------------------------------------------------

subroutine X(states_elec_set_state1)(st, mesh, idim, ist, iqn, psi)
  type(states_elec_t), intent(inout) :: st
  class(mesh_t),       intent(in)    :: mesh
  integer,             intent(in)    :: idim   !< current dimension
  integer,             intent(in)    :: ist    !< current state
  integer,             intent(in)    :: iqn    !< current k-point
  R_TYPE,              intent(in)    :: psi(:)

  PUSH_SUB(X(states_elec_set_state1))

  call batch_set_state(st%group%psib(st%group%iblock(ist, iqn), iqn), (/ist, idim/), mesh%np, psi)

  POP_SUB(X(states_elec_set_state1))
end subroutine X(states_elec_set_state1)

! ------------------------------------------------------------

subroutine X(states_elec_set_state3)(st, mesh, iqn, psi)
  type(states_elec_t), intent(inout) :: st
  class(mesh_t),       intent(in)    :: mesh
  integer,             intent(in)    :: iqn
  R_TYPE,              intent(in)    :: psi(:, :, st%st_start:)

  integer :: ist

  PUSH_SUB(X(states_elec_set_state3))

  do ist = st%st_start, st%st_end
    call states_elec_set_state(st, mesh, ist, iqn, psi(:, :, ist))
  end do

  POP_SUB(X(states_elec_set_state3))
end subroutine X(states_elec_set_state3)

! ------------------------------------------------------------

subroutine X(states_elec_set_state4)(st, mesh, psi)
  type(states_elec_t), intent(inout) :: st
  class(mesh_t),       intent(in)    :: mesh
  R_TYPE,              intent(in)    :: psi(:, :, st%st_start:, st%d%kpt%start:)

  integer :: iqn, ist

  PUSH_SUB(X(states_elec_set_state4))

  do iqn = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      call states_elec_set_state(st, mesh, ist, iqn, psi(:, :, ist, iqn))
    end do
  end do

  POP_SUB(X(states_elec_set_state4))
end subroutine X(states_elec_set_state4)

! ------------------------------------------------------------

!> Returns the value of all the states in the range of points
!> [start_point:end_point].

subroutine X(states_elec_get_points1)(st, start_point, end_point, iqn, psi)
  type(states_elec_t),  intent(in)    :: st
  integer,              intent(in)    :: start_point
  integer,              intent(in)    :: end_point
  integer,              intent(in)    :: iqn
  R_TYPE,               intent(out)   :: psi(:, :, :)

  integer :: ib

  PUSH_SUB(X(states_elec_get_points1))

  do ib = st%group%block_start, st%group%block_end
    call batch_get_points(st%group%psib(ib, iqn), start_point, end_point, psi)
  end do

  POP_SUB(X(states_elec_get_points1))
end subroutine X(states_elec_get_points1)

! ------------------------------------------------------------
! ------------------------------------------------------------

!> Returns the value of all the states in the range of points
!> [start_point:end_point].

subroutine X(states_elec_get_points2)(st, start_point, end_point, psi)
  type(states_elec_t), intent(in)    :: st
  integer,             intent(in)    :: start_point
  integer,             intent(in)    :: end_point
  R_TYPE,              intent(out)   :: psi(:, :, :, :)

  integer :: iqn

  PUSH_SUB(X(states_elec_get_points2))

  do iqn = st%d%kpt%start, st%d%kpt%end
    call X(states_elec_get_points1)(st, start_point, end_point, iqn, psi(:, :, :, iqn))
  end do

  POP_SUB(X(states_elec_get_points2))
end subroutine X(states_elec_get_points2)


!> @brief Generate a random vector.
!!
!! @warning This will not work for SPINORS when the spins are not fixed.
subroutine X(states_elec_generate_random_vector)(mesh, st, vector, normalized)
  class(mesh_t),          intent(in)    :: mesh           !< System grid.
  type(states_elec_t),    intent(in)    :: st             !< Provide information on the basis.
  R_TYPE,                 intent(out)   :: vector(:,  :)  !< Random vector.
  logical, optional,      intent(in)    :: normalized     !< Normalize the random vector.

  type(batch_t) :: ffb              !< A batch which contains the mesh functions whose points will be exchanged.
  logical       :: normalized_      !< Local instance of normalized
  integer       :: n_spin           !< Number of spin channels
  integer       :: ispin            !< Spin index

  PUSH_SUB(X(states_elec_generate_random_vector))

  if ((st%d%ispin == SPINORS) .and. st%fixed_spins) then
    call messages_not_implemented('`generate_random_vector` with a fixed spinor direction')
  endif

  normalized_ = optional_default(normalized, .true.)
  n_spin = size(vector, dim=2)

  do ispin = 1, n_spin
    call X(mf_random)(mesh, vector(1:mesh%np, ispin), &
      pre_shift = mesh%pv%xlocal-1, &
      post_shift = mesh%pv%np_global - mesh%pv%xlocal - mesh%np + 1, &
      normalized = .false.)
  enddo

  ! Ensures that the grid points are properly distributed in the domain parallel case
  if(mesh%parallel_in_domains) then
    do ispin = 1, n_spin
      call batch_init(ffb, vector(1:mesh%np, ispin))
      call X(mesh_batch_exchange_points)(mesh, ffb, backward_map = .true.)
      call ffb%end()
    enddo
  end if

  if (normalized_) call X(mf_normalize)(mesh, n_spin, vector)

  POP_SUB(X(states_elec_generate_random_vector))
end subroutine X(states_elec_generate_random_vector)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
