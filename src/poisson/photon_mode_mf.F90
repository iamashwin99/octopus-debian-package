!! Copyright (C) 2017 Johannes Flick
!! Copyright (C) 2021 Davis Welakuh (merged)
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

!!   v_mxc = - omega*q_alpha*(lambda*r) + (lambda*dipole)*(lambda*r)
!!   see PRL 115, 093001 (2015) and PRL 121, 113002 (2018)
#include "global.h"

module photon_mode_mf_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use io_function_oct_m
  use io_oct_m
  use ions_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use photon_mode_oct_m
  use profiling_oct_m
  use ps_oct_m
  use restart_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m


  implicit none

  private

  public ::                &
    mf_t,                  &
    mf_init,               &
    mf_end,                &
    mf_calc,               &
    mf_photons_dump,       &
    mf_photons_load

  type mf_t
    FLOAT,   pointer     :: vmf(:)
    FLOAT,   allocatable :: dipole(:, :)
    FLOAT,   allocatable :: dipole_former(:, :)
    CMPLX,   allocatable :: integral(:)
    FLOAT,   allocatable :: pt_q(:), pt_p(:)
    FLOAT,   allocatable :: pt_q_former(:)
    FLOAT                :: time_former
    FLOAT,   pointer     :: fmf(:)  !! meanfield force
    logical              :: has_restart

  end type mf_t

contains

  subroutine mf_init(this, gr, st, ions, pt_mode)
    type(mf_t),          intent(out)  :: this
    type(grid_t),        intent(in)   :: gr
    type(states_elec_t), intent(in)   :: st
    type(ions_t),        intent(in)   :: ions
    type(photon_mode_t), intent(in)   :: pt_mode

    FLOAT, allocatable :: e_dip(:, :), n_dip(:)
    integer :: jj, ispin, ions_dim

    PUSH_SUB(mf_init)

    ions_dim = gr%mesh%box%dim

    SAFE_ALLOCATE(e_dip(1:ions_dim+1, 1:st%d%nspin))
    SAFE_ALLOCATE(n_dip(1:ions_dim))

    SAFE_ALLOCATE(this%vmf(1:gr%mesh%np))
    SAFE_ALLOCATE(this%dipole(1:ions_dim, 1:st%d%nspin))
    SAFE_ALLOCATE(this%dipole_former(1:ions_dim, 1:st%d%nspin))

    SAFE_ALLOCATE(this%integral(1:pt_mode%nmodes))
    SAFE_ALLOCATE(this%pt_q(1:pt_mode%nmodes))
    SAFE_ALLOCATE(this%pt_p(1:pt_mode%nmodes))
    SAFE_ALLOCATE(this%pt_q_former(1:pt_mode%nmodes))
    SAFE_ALLOCATE(this%fmf(1:ions_dim))

    this%vmf = M_ZERO
    this%has_restart = .false.
    !call dmf_multipoles(gr%mesh, st%rho(:, 1), 1, e_dip(:, 1))
    do ispin = 1, st%d%nspin
      call dmf_multipoles(gr%mesh, st%rho(:, ispin), 1, e_dip(:, ispin))
    end do

    n_dip = ions%dipole()
    do jj = 1, ions_dim
      e_dip(jj+1, 1) = sum(e_dip(jj+1, :))
      this%dipole(jj,1) = - n_dip(jj) - e_dip(jj+1, 1)  ! dipole moment <mu_el> = \sum_i -e <x_i>
    end do

    this%dipole_former = M_ZERO
    this%integral = M_ZERO

    if (pt_mode%has_q0_p0) then
      this%pt_q(1:pt_mode%nmodes) = pt_mode%pt_coord_q0(1:pt_mode%nmodes)
      this%pt_p(1:pt_mode%nmodes) = pt_mode%pt_momen_p0(1:pt_mode%nmodes)
    else
      this%pt_q = M_ZERO
      this%pt_p = M_ZERO
    end if
    
    this%pt_q_former = M_ZERO
    this%time_former = M_ZERO

    this%fmf(1:ions_dim) = M_ZERO

    SAFE_DEALLOCATE_A(e_dip)
    SAFE_DEALLOCATE_A(n_dip)

    POP_SUB(mf_init)
  end subroutine mf_init

!------------------------------------------

  subroutine mf_end(this)
    type(mf_t), intent(inout) :: this

    PUSH_SUB(mf_end)

    SAFE_DEALLOCATE_P(this%vmf)
    SAFE_DEALLOCATE_A(this%dipole_former)
    SAFE_DEALLOCATE_A(this%dipole)

    SAFE_DEALLOCATE_A(this%integral)
    SAFE_DEALLOCATE_A(this%pt_p)
    SAFE_DEALLOCATE_A(this%pt_q)
    SAFE_DEALLOCATE_A(this%pt_q)
    SAFE_DEALLOCATE_A(this%pt_q_former)

    SAFE_DEALLOCATE_P(this%fmf)

    POP_SUB(mf_end)
  end subroutine mf_end

!------------------------------------------

  subroutine mf_calc(this, gr, st, ions, pt_mode, time)
    type(mf_t),          intent(inout)    :: this
    type(grid_t),        intent(inout)    :: gr
    type(states_elec_t), intent(inout)    :: st
    type(ions_t),        intent(in)       :: ions
    type(photon_mode_t), intent(in)       :: pt_mode
    FLOAT,               intent(in)       :: time

    FLOAT :: lambda_pol_dipole, lambda_pol_dipole_former
    CMPLX :: integrand
    FLOAT, allocatable :: e_dip(:, :), n_dip(:)
    FLOAT :: q0, p0
    integer :: ii, jj, ispin, ions_dim
    logical, save :: first = .true.

    PUSH_SUB(mf_calc)

    ions_dim = gr%mesh%box%dim

    SAFE_ALLOCATE(e_dip(1:ions_dim+1, 1:st%d%nspin))
    SAFE_ALLOCATE(n_dip(1:ions_dim))

    if (.not. (first .and. this%has_restart)) then
      do ii = 1, pt_mode%nmodes
        this%pt_q_former(ii) = this%pt_q(ii)
      end do
      this%dipole_former(:, :) = this%dipole(:, :)

      do ispin = 1, st%d%nspin
        call dmf_multipoles(gr%mesh, st%rho(:, ispin), 1, e_dip(:, ispin))
      end do

      n_dip = ions%dipole()
      do jj = 1, ions_dim
        e_dip(jj+1, 1) = sum(e_dip(jj+1, :))
        this%dipole(jj,1) = - n_dip(jj) - e_dip(jj+1, 1)  ! dipole moment <mu_el> = \sum_i -e <x_i>
      end do
    end if

    this%vmf(1:gr%mesh%np) = M_ZERO
    this%fmf(1:ions_dim) = M_ZERO

    do ii = 1, pt_mode%nmodes
      lambda_pol_dipole = M_ZERO
      lambda_pol_dipole_former = M_ZERO
      do jj = 1, ions_dim
        lambda_pol_dipole = lambda_pol_dipole +  &
          pt_mode%lambda(ii)*pt_mode%pol(ii, jj)*this%dipole(jj, 1)
        lambda_pol_dipole_former = lambda_pol_dipole_former + &
          pt_mode%lambda(ii)*pt_mode%pol(ii, jj)*this%dipole_former(jj, 1)
      end do

      if (.not.(first .and. this%has_restart)) then
        if (time == M_ZERO) then
          this%integral(ii) = M_ZERO
        else if (this%integral(ii) == M_ZERO) then
          this%integral(ii) = M_HALF*lambda_pol_dipole*exp(-M_zI*pt_mode%omega(ii)*(this%time_former))
        else
          this%integral(ii) = this%integral(ii) + lambda_pol_dipole*exp(-M_zI*pt_mode%omega(ii)*(this%time_former))
        end if
      end if

      integrand = -(this%integral(ii)*exp(M_zI*pt_mode%omega(ii)*(time)*pt_mode%mu))* &
        (time*pt_mode%mu - this%time_former)

      this%pt_q(ii) = aimag(integrand)
      this%pt_p(ii) = pt_mode%omega(ii)*real(integrand)

      if (pt_mode%has_q0_p0) then
        q0 = pt_mode%pt_coord_q0(ii)*cos(pt_mode%omega(ii)*(time)*pt_mode%mu)
        q0 = q0 + pt_mode%pt_momen_p0(ii)/pt_mode%omega(ii)*sin(pt_mode%omega(ii)*(time)*pt_mode%mu)
        this%pt_q(ii) = this%pt_q(ii) + q0

        p0 = -pt_mode%pt_coord_q0(ii)*pt_mode%omega(ii)*sin(pt_mode%omega(ii)*(time)*pt_mode%mu)
        p0 = p0 + pt_mode%pt_momen_p0(ii)*cos(pt_mode%omega(ii)*(time)*pt_mode%mu)
        this%pt_p(ii) = this%pt_p(ii) + p0
      end if

      ! we need the negative sign due to the electric pol_dipole
      this%vmf(1:gr%mesh%np) = this%vmf(1:gr%mesh%np) - &
        M_HALF*((lambda_pol_dipole + lambda_pol_dipole_former) + pt_mode%omega(ii)* &
        (this%pt_q(ii) + this%pt_q_former(ii)))*(pt_mode%lambda(ii)*pt_mode%pol_dipole(1:gr%mesh%np,ii))
      do jj = 1, ions_dim
        this%fmf(jj) = this%fmf(jj) - pt_mode%omega(ii)*pt_mode%lambda(ii)* &
          pt_mode%pol(ii, jj)*(this%pt_q(ii) + lambda_pol_dipole/pt_mode%omega(ii)) !minus?
      end do
    end do

    if (first .and. this%has_restart) then
      first = .false.
    end if

    this%time_former = time*pt_mode%mu

    SAFE_DEALLOCATE_A(e_dip)
    SAFE_DEALLOCATE_A(n_dip)

    POP_SUB(mf_calc)
  end subroutine mf_calc

! ---------------------------------------------------------

  subroutine mf_photons_dump(restart, this, gr, dt, pt_mode, ierr)
    type(restart_t), intent(in)  :: restart
    type(mf_t),      intent(in)  :: this
    type(grid_t),        intent(in)    :: gr
    FLOAT,               intent(in)    :: dt
    type(photon_mode_t), intent(in)    :: pt_mode
    integer,         intent(out) :: ierr

    character(len=80), allocatable :: lines(:)
    integer :: iunit, err, jj, ions_dim

    PUSH_SUB(mf_photons_dump)
    ions_dim = gr%mesh%box%dim

    SAFE_ALLOCATE(lines(1:2*ions_dim + 4))

    ierr = 0

    iunit = restart_open(restart, 'photon_mf')
    write(lines(1), '(a10,2x,es19.12)') 'pt_integral_real', real(this%integral(1))
    write(lines(2), '(a10,2x,es19.12)') 'pt_integral_aimag', aimag(this%integral(1))
    write(lines(3), '(a10,2x,es19.12)') 'pt_q_former', this%pt_q_former(1)
    write(lines(4), '(a10,2x,es19.12)') 'pt_time_former', this%time_former - dt*pt_mode%mu
    do jj = 1, ions_dim
      write(lines(4 + jj), '(a10,2x,es19.12)') 'dipole', this%dipole(jj, 1)
    end do
    do jj = 1, ions_dim
      write(lines(4 + ions_dim + jj), '(a10,2x,es19.12)') 'dipole_former', this%dipole_former(jj, 1)
    end do
    call restart_write(restart, iunit, lines, 2*ions_dim + 4, err)
    if (err /= 0) ierr = ierr + 1
    call restart_close(restart, iunit)

    SAFE_DEALLOCATE_A(lines)

    POP_SUB(mf_photons_dump)
  end subroutine mf_photons_dump

! ---------------------------------------------------------

  subroutine mf_photons_load(restart, this, gr, ierr)
    type(restart_t),   intent(in)    :: restart
    type(mf_t),        intent(inout) :: this
    type(grid_t),      intent(in)    :: gr
    integer,           intent(out)   :: ierr

    integer :: err, iunit, jj, ions_dim
    character(len=128), allocatable :: lines(:)
    character(len=7) :: dummy
    FLOAT, allocatable :: rr(:)

    PUSH_SUB(mf_photons_load)

    ierr = 0
    ions_dim = gr%mesh%box%dim

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(mf_photons_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading Photons restart."
      call messages_info(1, namespace=restart%namespace)
    end if

    SAFE_ALLOCATE(rr(1:2*ions_dim + 4))
    SAFE_ALLOCATE(lines(1:2*ions_dim + 4))
    iunit = restart_open(restart, 'photon_mf')
    call restart_read(restart, iunit, lines, 2*ions_dim + 4, err)
    if (err /= 0) then
      ierr = ierr + 1
    else
      do jj = 1, 2*ions_dim + 4
        read(lines(jj),'(a10,2x,es19.12)') dummy, rr(jj)
      end do

      this%integral(1) = rr(1) + M_zI*rr(2)
      this%pt_q_former(1) = rr(3)
      this%time_former = rr(4)
      do jj = 1, ions_dim
        this%dipole(jj, 1) = rr(jj + 4)
      end do
      do jj = 1, ions_dim
        this%dipole_former(jj, 1) = rr(jj + ions_dim + 4)
      end do
      this%has_restart = .true.

    end if
    call restart_close(restart, iunit)

    if (debug%info) then
      message(1) = "Debug: Reading Photons restart done."
      call messages_info(1, namespace=restart%namespace)
    end if

    SAFE_DEALLOCATE_A(rr)
    SAFE_DEALLOCATE_A(lines)

    POP_SUB(mf_photons_load)
  end subroutine mf_photons_load

end module photon_mode_mf_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
