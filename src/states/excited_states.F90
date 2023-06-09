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

module excited_states_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m

  implicit none

  private
  public ::                         &
    states_pair_t,                  &
    excited_states_t,               &
    excited_states_init,            &
    excited_states_kill,            &
    excited_states_output,          &
    dstates_elec_mpdotp,                 &
    zstates_elec_mpdotp,                 &
    zstates_elec_mpmatrixelement,        &
    dstates_elec_mpmatrixelement,        &
    dstates_elec_matrix_swap,            &
    zstates_elec_matrix_swap



  interface dstates_elec_mpdotp
    module procedure dstates_elec_mpdotp_g, dstates_elec_mpdotp_x
  end interface dstates_elec_mpdotp

  interface zstates_elec_mpdotp
    module procedure zstates_elec_mpdotp_g, zstates_elec_mpdotp_x
  end interface zstates_elec_mpdotp

  interface dstates_elec_mpmatrixelement
    module procedure dstates_elec_mpmatrixelement_g
  end interface dstates_elec_mpmatrixelement

  interface zstates_elec_mpmatrixelement
    module procedure zstates_elec_mpmatrixelement_g
  end interface zstates_elec_mpmatrixelement

  type states_pair_t
    ! Components are public by default
    integer :: i
    integer :: a
    integer :: kk
  end type states_pair_t

  type excited_states_t
    ! Components are public by default
    type(states_elec_t),     pointer :: st
    integer                          :: n_pairs
    type(states_pair_t), allocatable :: pair(:)
    FLOAT,               allocatable :: weight(:)
  end type excited_states_t

contains

  ! ---------------------------------------------------------
  !> Fills in an excited_state structure, by reading a file called
  !! "filename". This file describes the "promotions" from occupied
  !! to unoccupied levels that change the initial Slater determinant
  !! structure specified in ground_state. These promotions are a set
  !! of electron-hole pairs. The structure of the file is thus four
  !! columns:
  !!
  !! i  a  k  weight
  !!
  !! where i should be an occupied state, a an unoccupied one, and k
  !! the k-point (including spin) of the corresponding orbital. This pair is then associated with a
  !! creation-annihilation pair a^t_{a,k} a_{i,k}, so that the
  !! excited state will be a linear combination in the form:
  !!
  !! |ExcitedState> = Sum [ weight(i,a,k) a^t_{a,k} a_{i,k} |GroundState> ]
  !!
  !! where weight is the number in the fourth column.
  !! These weights should be normalized to one; otherwise the routine
  !! will normalize them, and write a warning.
  !!
  !! This file structure is the one written by the casida run mode, in the files
  !! in the directory "excitations".
  !! ---------------------------------------------------------
  subroutine excited_states_init(excited_state, ground_state, filename, namespace)
    type(excited_states_t),      intent(inout) :: excited_state
    type(states_elec_t), target, intent(in)    :: ground_state
    character(len=*),            intent(in)    :: filename
    type(namespace_t),           intent(in)    :: namespace

    integer :: iunit, nst, ispin, nik, &
      n_possible_pairs, ipair, jpair, ist, ios, nspin, trash
    integer, allocatable :: n_filled(:), n_partially_filled(:), n_half_filled(:), n_empty(:), &
      filled(:, :), partially_filled(:, :), half_filled(:, :)
    FLOAT :: dump
    logical :: ok

    PUSH_SUB(excited_states_init)

    ! This is just to make the code more readable.
    nst   = ground_state%nst
    nik   = ground_state%d%nik
    ispin = ground_state%d%ispin
    nspin = ground_state%d%nspin

    if (nik > 2 .or. ( (nik == 2) .and. (ispin /= SPIN_POLARIZED))) then
      message(1) = 'Cannot calculate projections onto excited states for periodic systems.'
      call messages_fatal(1, namespace=namespace)
    end if

    SAFE_ALLOCATE(          n_filled(1:nspin))
    SAFE_ALLOCATE(           n_empty(1:nspin))
    SAFE_ALLOCATE(n_partially_filled(1:nspin))
    SAFE_ALLOCATE(     n_half_filled(1:nspin))
    SAFE_ALLOCATE(          filled(1:nst, 1:nspin))
    SAFE_ALLOCATE(partially_filled(1:nst, 1:nspin))
    SAFE_ALLOCATE(     half_filled(1:nst, 1:nspin))

    select case (ispin)
    case (UNPOLARIZED)
      call occupied_states(ground_state, namespace, 1, n_filled(1), n_partially_filled(1), n_half_filled(1), &
        filled(:, 1), partially_filled(:, 1), half_filled(:, 1))
      if (n_partially_filled(1) > 0) then
        message(1) = 'Cannot calculate projections onto excited states if there are partially filled orbitals.'
        call messages_fatal(1, namespace=namespace)
      end if
      ! We will not accept, for the time being, constructing excited states in spin-restricted mode if
      ! there are single-particle states that are half-filled, *unless* there is only one state and it is
      ! half-filled (single-particle calculation).
      if ((n_half_filled(1) /= 0 .and. n_filled(1) > 0) .or. (n_half_filled(1) > 1)) then
        message(1) = 'Cannot construct excited states from ground states that contain half-filled'
        message(2) = 'orbitals - unless they are just one-particle states with only one half-filled'
        message(3) = 'orbital and no doubly occupied ones. Try using the spin-unrestricted mode.'
        call messages_fatal(3, namespace=namespace)
      end if
      if (n_half_filled(1) == 0) then
        n_empty(1) = nst - n_filled(1)
        n_possible_pairs = n_filled(1) * n_empty(1)
      else ! This is for the one-electron case.
        n_empty(1) = nst - 1
        n_possible_pairs = n_empty(1)
      end if

    case (SPIN_POLARIZED)
      call occupied_states(ground_state, namespace, 1, n_filled(1), n_partially_filled(1), n_half_filled(1), &
        filled(:, 1), partially_filled(:, 1), half_filled(:, 1))
      call occupied_states(ground_state, namespace, 2, n_filled(2), n_partially_filled(2), n_half_filled(2), &
        filled(:, 2), partially_filled(:, 2), half_filled(:, 2))
      if (n_partially_filled(1) * n_partially_filled(2) > 0) then
        message(1) = 'Cannot calculate projections onto excited states if there are partially filled orbitals.'
        call messages_fatal(1, namespace=namespace)
      end if
      n_empty(1) = nst - n_filled(1)
      n_empty(2) = nst - n_filled(2)
      n_possible_pairs = n_filled(1) * n_empty(1) + n_filled(2) * n_empty(2)

    case (SPINORS)
      call occupied_states(ground_state, namespace, 1, n_filled(1), n_partially_filled(1), n_half_filled(1), &
        filled(:, 1), partially_filled(:, 1), half_filled(:, 1))
      if (n_partially_filled(1) > 0) then
        message(1) = 'Cannot calculate projections onto excited states if there are partially filled orbitals.'
        call messages_fatal(1, namespace=namespace)
      end if
      n_empty(1) = nst - n_filled(1)
      n_possible_pairs = n_filled(1) * n_empty(1)
    end select

    iunit = io_open(trim(filename), namespace, action = 'read', status = 'old', die = .true.)
    call io_skip_header(iunit)

    ! Now we count the number of pairs in the file
    ipair = 0
    do
      read(iunit, *, end = 101)
      backspace(iunit)
      read(iunit, *, iostat = ios) trash, trash, trash, dump
      if (ios /= 0) then
        message(1) = 'Error attempting to read the electron-hole pairs in file "'//trim(filename)//'"'
        call messages_fatal(1, namespace=namespace)
      end if
      ipair = ipair + 1
    end do
101 continue
    if (ipair == 0) then
      message(1) = 'File "'//trim(filename)//'" is empty?'
      call messages_fatal(1, namespace=namespace)
    elseif (ipair > n_possible_pairs) then
      message(1) = 'File "'//trim(filename)//'" contains too many electron-hole pairs.'
      call messages_fatal(1, namespace=namespace)
    end if

    excited_state%n_pairs = ipair
    SAFE_ALLOCATE(excited_state%pair(1:ipair))
    SAFE_ALLOCATE(excited_state%weight(1:ipair))

    rewind(iunit)
    call io_skip_header(iunit)
    do ipair = 1, excited_state%n_pairs
      read(iunit, *) excited_state%pair(ipair)%i, excited_state%pair(ipair)%a, &
        excited_state%pair(ipair)%kk, excited_state%weight(ipair)
      if (((ispin == UNPOLARIZED) .or. (ispin == SPINORS)) .and. (excited_state%pair(ipair)%kk /= 1)) then
        message(1) = 'Error reading excited state in file "'//trim(filename)//'":'
        message(2) = 'Cannot treat a electron-hole pair for "down" spin when not working in spin-polarized mode.'
        call messages_fatal(2, namespace=namespace)
      end if
      ! Check whether it is a legitimate electron-hole swap.

      ! First, whether the occupied state belongs to the list of occupied states.
      ok = .false.
      do ist = 1, n_filled(excited_state%pair(ipair)%kk)
        ok = excited_state%pair(ipair)%i == filled(ist, excited_state%pair(ipair)%kk)
        if (ok) exit
      end do

      ! Treat differently the one-electron case in unpolarized mode
      if (ispin == UNPOLARIZED .and. (n_half_filled(1) == 1)) then
        ok = excited_state%pair(ipair)%i == half_filled(1, excited_state%pair(ipair)%kk)
      end if

      if (.not. ok) then
        write(message(1),'(a6,i3,a1,i3,a8,i1,a)') 'Pair (', excited_state%pair(ipair)%i, ',', &
          excited_state%pair(ipair)%a, '; k =', excited_state%pair(ipair)%kk, ') is not valid.'
        call messages_fatal(1, namespace=namespace)
      end if

      ! Then, whether the unoccupied state is really unoccupied.
      ok = .true.
      do ist = 1, n_filled(excited_state%pair(ipair)%kk)
        ok = .not. (excited_state%pair(ipair)%a == filled(ist, excited_state%pair(ipair)%kk))
      end do
      ! Treat differently the one-electron case in unpolarized mode
      if (ispin == UNPOLARIZED .and. (n_half_filled(1) == 1)) then
        ok = .not. (excited_state%pair(ipair)%i == half_filled(1, excited_state%pair(ipair)%kk))
      end if
      ok = .not. (excited_state%pair(ipair)%a > nst)
      if (.not. ok) then
        write(message(1),'(a6,i3,a1,i3,a8,i1,a)') 'Pair (', excited_state%pair(ipair)%i, ',', &
          excited_state%pair(ipair)%a, '; k =', excited_state%pair(ipair)%kk, ') is not valid.'
        call messages_fatal(1, namespace=namespace)
      end if

      ! Now, we check that there are no repetitions:
      do jpair = 1, ipair - 1
        ok = .not. pair_is_eq(excited_state%pair(ipair), excited_state%pair(jpair))
        if (.not. ok) then
          write(message(1),'(a6,i3,a1,i3,a8,i1,a)') 'Pair (', excited_state%pair(ipair)%i, ',', &
            excited_state%pair(ipair)%a, '; k =', excited_state%pair(ipair)%kk, ') is repeated in the file.'
          call messages_fatal(1, namespace=namespace)
        end if
      end do

    end do

    ! Now we point to the ground state from which the excited state is defined.
    excited_state%st => ground_state

    ! Check the normalization.
    dump = sum(excited_state%weight(1:excited_state%n_pairs)**2)
    if (.not. abs(dump - M_ONE) < CNST(1.0e-5)) then
      excited_state%weight(1:excited_state%n_pairs) = excited_state%weight(1:excited_state%n_pairs) / sqrt(dump)
      message(1) = 'The excited state in file "'//trim(filename)//'" was not normalized.'
      call messages_warning(1, namespace=namespace)
    end if

    call io_close(iunit)

    SAFE_DEALLOCATE_A(n_filled)
    SAFE_DEALLOCATE_A(n_partially_filled)
    SAFE_DEALLOCATE_A(n_half_filled)
    SAFE_DEALLOCATE_A(n_empty)
    SAFE_DEALLOCATE_A(filled)
    SAFE_DEALLOCATE_A(partially_filled)
    SAFE_DEALLOCATE_A(half_filled)
    POP_SUB(excited_states_init)
  end subroutine excited_states_init


  ! ---------------------------------------------------------
  !> Kills an excited_state structure.
  subroutine excited_states_kill(excited_state)
    type(excited_states_t), intent(inout) :: excited_state

    PUSH_SUB(excited_states_kill)

    nullify(excited_state%st)
    SAFE_DEALLOCATE_A(excited_state%pair)
    SAFE_DEALLOCATE_A(excited_state%weight)

    POP_SUB(excited_states_kill)
  end subroutine excited_states_kill


  ! ---------------------------------------------------------
  subroutine excited_states_output(excited_state, dirname, namespace)
    type(excited_states_t), intent(in) :: excited_state
    character(len=*),       intent(in) :: dirname
    type(namespace_t),      intent(in) :: namespace

    integer :: iunit, ipair

    PUSH_SUB(excited_states_output)

    iunit = io_open(trim(dirname)//'/excitations', namespace, action = 'write', status = 'replace')
    do ipair = 1, excited_state%n_pairs
      write(iunit, '(3i5,es20.12)') excited_state%pair(ipair)%i, excited_state%pair(ipair)%a, &
        excited_state%pair(ipair)%kk, excited_state%weight(ipair)
    end do

    call io_close(iunit)
    POP_SUB(excited_states_output)
  end subroutine excited_states_output


  ! ---------------------------------------------------------
  logical function pair_is_eq(pair1, pair2) result(res)
    type(states_pair_t), intent(in) :: pair1, pair2

    res = (pair1%i == pair2%i) .and. (pair1%a == pair2%a) .and. (pair1%kk == pair2%kk)
  end function pair_is_eq

#include "undef.F90"
#include "real.F90"
#include "excited_states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "excited_states_inc.F90"
#include "undef.F90"

end module excited_states_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
