!! Copyright (C) 2008-2019 X. Andrade, F. Bonafe, R. Jestaedt, H. Appel
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

module current_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use exchange_operator_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_base_oct_m
  use hamiltonian_elec_oct_m
  use lalg_basic_oct_m
  use lda_u_oct_m
  use math_oct_m
  use magnetic_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use projector_oct_m
  use scissor_oct_m
  use space_oct_m
  use states_elec_dim_oct_m
  use states_elec_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use wfs_elec_oct_m
  use xc_oct_m

  implicit none

  private

  type current_t
    private
    integer :: method
  end type current_t


  public ::                               &
    current_t,                            &
    current_init,                         &
    current_calculate,                    &
    current_calculate_mag,                &
    current_heat_calculate,               &
    current_calculate_mel

  integer, parameter, public ::           &
    CURRENT_GRADIENT           = 1,       &
    CURRENT_GRADIENT_CORR      = 2,       &
    CURRENT_HAMILTONIAN        = 3

contains

  subroutine current_init(this, namespace)
    type(current_t),   intent(out)   :: this
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(current_init)

    !%Variable CurrentDensity
    !%Default gradient_corrected
    !%Type integer
    !%Section Hamiltonian
    !%Description
    !% This variable selects the method used to
    !% calculate the current density. For the moment this variable is
    !% for development purposes and users should not need to use
    !% it.
    !%Option gradient 1
    !% The calculation of current is done using the gradient operator. (Experimental)
    !%Option gradient_corrected 2
    !% The calculation of current is done using the gradient operator
    !% with additional corrections for the total current from non-local operators.
    !%Option hamiltonian 3
    !% The current density is obtained from the commutator of the
    !% Hamiltonian with the position operator. (Experimental)
    !%End

    call parse_variable(namespace, 'CurrentDensity', CURRENT_GRADIENT_CORR, this%method)
    if (.not. varinfo_valid_option('CurrentDensity', this%method)) call messages_input_error(namespace, 'CurrentDensity')
    if (this%method /= CURRENT_GRADIENT_CORR) then
      call messages_experimental("CurrentDensity /= gradient_corrected")
    end if

    POP_SUB(current_init)
  end subroutine current_init

  ! ---------------------------------------------------------

  subroutine current_batch_accumulate(st, der, ik, ib, psib, gpsib)
    type(states_elec_t), intent(inout) :: st
    type(derivatives_t), intent(inout) :: der
    integer,             intent(in)    :: ik
    integer,             intent(in)    :: ib
    type(wfs_elec_t),    intent(in)    :: psib
    class(wfs_elec_t),   intent(in)    :: gpsib(:)

    integer :: ist, idir, ii, ip, idim, wgsize
    CMPLX, allocatable :: psi(:, :), gpsi(:, :)
    FLOAT, allocatable :: current_tmp(:, :)
    CMPLX :: c_tmp
    FLOAT :: ww
    FLOAT, allocatable :: weight(:)
    type(accel_mem_t) :: buff_weight, buff_current
    type(accel_kernel_t), save :: kernel

    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np_part, 1:st%d%dim))

    if (st%d%ispin == SPINORS .or. (psib%status() == BATCH_DEVICE_PACKED .and. der%dim /= 3)) then

      do idir = 1, der%dim
        do ist = states_elec_block_min(st, ib), states_elec_block_max(st, ib)

          ww = st%d%kweights(ik)*st%occ(ist, ik)
          if (abs(ww) <= M_EPSILON) cycle

          do idim = 1, st%d%dim
            ii = st%group%psib(ib, ik)%inv_index((/ist, idim/))
            call batch_get_state(psib, ii, der%mesh%np, psi(:, idim))
            call batch_get_state(gpsib(idir), ii, der%mesh%np, gpsi(:, idim))
          end do

          if (st%d%ispin /= SPINORS) then
            !$omp parallel do
            do ip = 1, der%mesh%np
              st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) + ww*aimag(conjg(psi(ip, 1))*gpsi(ip, 1))
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(c_tmp)
            do ip = 1, der%mesh%np
              st%current_para(ip, idir, 1) = st%current_para(ip, idir, 1) + ww*aimag(conjg(psi(ip, 1))*gpsi(ip, 1))
              st%current_para(ip, idir, 2) = st%current_para(ip, idir, 2) + ww*aimag(conjg(psi(ip, 2))*gpsi(ip, 2))
              c_tmp = -M_HALF*M_ZI*(conjg(psi(ip, 2))*gpsi(ip, 1) - psi(ip, 1)*conjg(gpsi(ip, 2)))
              st%current_para(ip, idir, 3) = st%current_para(ip, idir, 3) + ww*TOFLOAT(c_tmp)
              st%current_para(ip, idir, 4) = st%current_para(ip, idir, 4) + ww*aimag(c_tmp)
            end do
            !$omp end parallel do
          end if

        end do
      end do

    else if (psib%status() == BATCH_DEVICE_PACKED) then

      ASSERT(der%dim == 3)

      SAFE_ALLOCATE(weight(1:psib%nst))
      do ist = 1, psib%nst
        weight(ist) = st%d%kweights(ik)*st%occ(psib%ist(ist), ik)
      end do


      call accel_create_buffer(buff_weight, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, psib%nst)
      call accel_write_buffer(buff_weight, psib%nst, weight)

      call accel_create_buffer(buff_current, ACCEL_MEM_WRITE_ONLY, TYPE_FLOAT, der%mesh%np*3)

      call accel_kernel_start_call(kernel, 'density.cl', 'current_accumulate')

      call accel_set_kernel_arg(kernel, 0, psib%nst)
      call accel_set_kernel_arg(kernel, 1, der%mesh%np)
      call accel_set_kernel_arg(kernel, 2, buff_weight)
      call accel_set_kernel_arg(kernel, 3, psib%ff_device)
      call accel_set_kernel_arg(kernel, 4, log2(int(psib%pack_size(1), i4)))
      call accel_set_kernel_arg(kernel, 5, gpsib(1)%ff_device)
      call accel_set_kernel_arg(kernel, 6, gpsib(2)%ff_device)
      call accel_set_kernel_arg(kernel, 7, gpsib(3)%ff_device)
      call accel_set_kernel_arg(kernel, 8, log2(int(gpsib(1)%pack_size(1), i4)))
      call accel_set_kernel_arg(kernel, 9, buff_current)

      wgsize = accel_kernel_workgroup_size(kernel)

      call accel_kernel_run(kernel, (/pad(der%mesh%np, wgsize)/), (/wgsize/))

      SAFE_ALLOCATE(current_tmp(1:der%dim, 1:der%mesh%np))

      call accel_finish()

      call accel_read_buffer(buff_current, der%mesh%np*3, current_tmp)

      !$omp parallel private(ip)
      do idir = 1, der%dim
        !$omp do
        do ip = 1, der%mesh%np
          st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) + current_tmp(idir, ip)
        end do
        !$omp end do nowait
      end do
      !$omp end parallel

      SAFE_DEALLOCATE_A(current_tmp)

      call accel_release_buffer(buff_weight)
      call accel_release_buffer(buff_current)

      SAFE_DEALLOCATE_A(weight)

    else

      ASSERT(psib%is_packed() .eqv. gpsib(1)%is_packed())

      !$omp parallel private(ip, ist, ww, idir)
      do ii = 1, psib%nst
        ist = states_elec_block_min(st, ib) + ii - 1
        ww = st%d%kweights(ik)*st%occ(ist, ik)
        if (abs(ww) <= M_EPSILON) cycle

        if (psib%is_packed()) then
          do idir = 1, der%dim
            !$omp do
            do ip = 1, der%mesh%np
              st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) &
                + ww*aimag(conjg(psib%zff_pack(ii, ip))*gpsib(idir)%zff_pack(ii, ip))
            end do
            !$omp end do nowait
          end do
        else
          do idir = 1, der%dim
            !$omp do
            do ip = 1, der%mesh%np
              st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) &
                + ww*aimag(conjg(psib%zff(ip, 1, ii))*gpsib(idir)%zff(ip, 1, ii))
            end do
            !$omp end do nowait
          end do
        end if
      end do
      !$omp end parallel

    end if

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(gpsi)

  end subroutine current_batch_accumulate


  ! ---------------------------------------------------------
  !> Compute total electronic current density.
  subroutine current_calculate(this, namespace, gr, hm, space, st)
    type(current_t),          intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(inout) :: gr
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(space_t),            intent(in)    :: space
    type(states_elec_t),      intent(inout) :: st

    type(profile_t), save :: prof

    call profiling_in(prof, "CURRENT_TOTAL")
    PUSH_SUB(current_calculate)

    call current_calculate_para(this, namespace, gr, hm, space, st)
    call current_calculate_dia_non_unif_vec_pot(this, namespace, gr%der, hm, st)

    ! Diamagnetic current from non-uniform vector potential and magnetization
    ! current not included by default in total current
    st%current = st%current_para

    call profiling_out(prof)
    POP_SUB(current_calculate)

  end subroutine current_calculate

  ! ---------------------------------------------------------
  !> Compute diamagnetic current density from non-uniform vector potential
  !> (the part coming from the uniform vector potential is already included in the paramagnetic term)
  subroutine current_calculate_dia_non_unif_vec_pot(this, namespace, der, hm, st)
    type(current_t),          intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(derivatives_t),      intent(inout) :: der
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(states_elec_t),      intent(inout) :: st

    type(profile_t), save :: prof
    integer :: ispin, idir, ip

    call profiling_in(prof, "CURRENT_DIA_NON_UNIF_A")
    PUSH_SUB(current_calculate_dia_non_unif_vec_pot)

    st%current_dia = M_ZERO

    if(allocated(hm%hm_base%vector_potential)) then
      do ispin = 1, st%d%nspin
        do idir = 1, der%dim
          !$omp parallel do
          do ip = 1, der%mesh%np
            ! the vector potential is assumed to be devided by c_0 already
            st%current_dia(ip, idir, ispin) = st%current_dia(ip, idir, ispin) - &
              st%rho(ip, ispin)*hm%hm_base%vector_potential(idir, ip)
          end do
          !$omp end parallel do
        end do
      end do
    end if

    call profiling_out(prof)
    POP_SUB(current_calculate_dia_non_unif_vec_pot)

  end subroutine current_calculate_dia_non_unif_vec_pot

  ! ---------------------------------------------------------
  !> Compute magnetization current
  !> Note: due to the the numerical curl, the magnetization current could deviate from the analytical
  !> one near the origin (see test for "current_density").
  subroutine current_calculate_mag(this, namespace, der, hm, st)
    type(current_t),          intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(derivatives_t),      intent(inout) :: der
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(states_elec_t),      intent(inout) :: st

    type(profile_t), save :: prof
    integer :: idir, ip
    FLOAT, allocatable :: magnetization_density(:, :), curl_mag(:, :)

    call profiling_in(prof, "CURRENT_MAG")
    PUSH_SUB(current_calculate_mag)

    st%current_mag = M_ZERO

    if (st%d%ispin /= UNPOLARIZED) then
      SAFE_ALLOCATE(magnetization_density(1:der%mesh%np_part, 1:der%dim))
      SAFE_ALLOCATE(curl_mag(1:der%mesh%np_part, 1:der%dim))

      call magnetic_density(der%mesh, st%d, st%rho, magnetization_density)
      call dderivatives_curl(der, magnetization_density, curl_mag)

      do idir = 1, der%dim
        !$omp parallel do
        do ip = 1, der%mesh%np
          st%current_mag(ip, idir, 1) = st%current_mag(ip, idir, 1) + M_HALF * curl_mag(ip, idir)
        end do
        !$omp end parallel do
      end do

      SAFE_DEALLOCATE_A(magnetization_density)
      SAFE_DEALLOCATE_A(curl_mag)
    end if

    call profiling_out(prof)
    POP_SUB(current_calculate_mag)

  end subroutine current_calculate_mag

  ! ---------------------------------------------------------
  !> Compute paramagnetic current density (including full diamagnetic term if method = Hamiltonian
  !> us used, and including diamagnetic from uniform vector potential otherwise)
  subroutine current_calculate_para(this, namespace, gr, hm, space, st)
    type(current_t),          intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(inout) :: gr
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(space_t),            intent(in)    :: space
    type(states_elec_t),      intent(inout) :: st

    integer :: ik, ist, idir, idim, ip, ib, ii, ispin
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), hpsi(:, :), rhpsi(:, :), rpsi(:, :), hrpsi(:, :)
    type(profile_t), save :: prof
    type(wfs_elec_t) :: hpsib, rhpsib, rpsib, hrpsib, epsib
    class(wfs_elec_t), allocatable :: commpsib(:)
    FLOAT :: ww
    CMPLX :: c_tmp

    call profiling_in(prof, "CURRENT_PARA")
    PUSH_SUB(current_calculate_para)

    ! spin not implemented or tested
    ASSERT(all(ubound(st%current_para) == (/gr%np_part, gr%der%dim, st%d%nspin/)))
    ASSERT(all(ubound(st%current_kpt) == (/gr%np, gr%der%dim, st%d%kpt%end/)))
    ASSERT(all(lbound(st%current_kpt) == (/1, 1, st%d%kpt%start/)))

    SAFE_ALLOCATE(psi(1:gr%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:gr%np, 1:gr%der%dim, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:gr%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rhpsi(1:gr%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rpsi(1:gr%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hrpsi(1:gr%np_part, 1:st%d%dim))
    SAFE_ALLOCATE_TYPE_ARRAY(wfs_elec_t, commpsib, (1:gr%der%dim))

    !$omp parallel private(idir, ip, ispin)
    do ik = st%d%kpt%start, st%d%kpt%end
      do idir = 1, gr%der%dim
        !$omp do
        do ip = 1, gr%np
          st%current_kpt(ip, idir, ik) = M_ZERO
        end do
        !$omp end do nowait
      end do
    end do
    do ispin = 1, st%d%nspin
      do idir = 1, gr%der%dim
        !$omp do
        do ip = 1, gr%np
          st%current_para(ip, idir, ispin) = M_ZERO
        end do
        !$omp end do nowait
      end do
    end do
    !$omp end parallel

    select case (this%method)

    case (CURRENT_HAMILTONIAN)

      do ik = st%d%kpt%start, st%d%kpt%end
        ispin = st%d%get_spin_index(ik)
        do ib = st%group%block_start, st%group%block_end

          call st%group%psib(ib, ik)%do_pack(copy = .true.)

          call st%group%psib(ib, ik)%copy_to(hpsib)
          call st%group%psib(ib, ik)%copy_to(rhpsib)
          call st%group%psib(ib, ik)%copy_to(rpsib)
          call st%group%psib(ib, ik)%copy_to(hrpsib)

          call boundaries_set(gr%der%boundaries, gr, st%group%psib(ib, ik))
          call zhamiltonian_elec_apply_batch(hm, namespace, gr, st%group%psib(ib, ik), hpsib, set_bc = .false.)

          do idir = 1, gr%der%dim

            call batch_mul(gr%np, gr%x(:, idir), hpsib, rhpsib)
            call batch_mul(gr%np_part, gr%x(:, idir), st%group%psib(ib, ik), rpsib)

            call zhamiltonian_elec_apply_batch(hm, namespace, gr, rpsib, hrpsib, set_bc = .false.)

            do ist = states_elec_block_min(st, ib), states_elec_block_max(st, ib)
              ww = st%d%kweights(ik)*st%occ(ist, ik)
              if (ww <= M_EPSILON) cycle

              do idim = 1, st%d%dim
                ii = st%group%psib(ib, ik)%inv_index((/ist, idim/))
                call batch_get_state(st%group%psib(ib, ik), ii, gr%np, psi(:, idim))
                call batch_get_state(hrpsib, ii, gr%np, hrpsi(:, idim))
                call batch_get_state(rhpsib, ii, gr%np, rhpsi(:, idim))
              end do

              if (st%d%ispin /= SPINORS) then
                !$omp parallel do
                do ip = 1, gr%np
                  st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) &
                    - ww*aimag(conjg(psi(ip, 1))*hrpsi(ip, 1) - conjg(psi(ip, 1))*rhpsi(ip, 1))
                end do
                !$omp end parallel do
              else
                !$omp parallel do  private(c_tmp)
                do ip = 1, gr%np
                  st%current_para(ip, idir, 1) = st%current_para(ip, idir, 1) + &
                    ww*aimag(conjg(psi(ip, 1))*hrpsi(ip, 1) - conjg(psi(ip, 1))*rhpsi(ip, 1))
                  st%current_para(ip, idir, 2) = st%current_para(ip, idir, 2) + &
                    ww*aimag(conjg(psi(ip, 2))*hrpsi(ip, 2) - conjg(psi(ip, 2))*rhpsi(ip, 2))
                  c_tmp = -M_HALF*M_ZI*(conjg(psi(ip, 2))*hrpsi(ip, 1) - conjg(psi(ip, 2))*rhpsi(ip, 1) &
                    -psi(ip, 1)*conjg(hrpsi(ip, 2)) - psi(ip, 1)*conjg(rhpsi(ip, 2)))
                  st%current_para(ip, idir, 3) = st%current_para(ip, idir, 3) + ww*TOFLOAT(c_tmp)
                  st%current_para(ip, idir, 4) = st%current_para(ip, idir, 4) + ww*aimag(c_tmp)
                end do
                !$omp end parallel do
              end if

            end do

          end do

          call st%group%psib(ib, ik)%do_unpack(copy = .false.)

          call hpsib%end()
          call rhpsib%end()
          call rpsib%end()
          call hrpsib%end()

        end do
      end do

    case (CURRENT_GRADIENT, CURRENT_GRADIENT_CORR)

      if (this%method == CURRENT_GRADIENT_CORR .and. .not. family_is_mgga_with_exc(hm%xc) &
        .and. hm%theory_level /= HARTREE_FOCK &
        .and. hm%theory_level /= GENERALIZED_KOHN_SHAM_DFT &
        .and. hm%theory_level /= RDMFT) then

        ! we can use the packed version

        do ik = st%d%kpt%start, st%d%kpt%end
          ispin = st%d%get_spin_index(ik)
          do ib = st%group%block_start, st%group%block_end

            call st%group%psib(ib, ik)%do_pack(copy = .true.)
            call st%group%psib(ib, ik)%copy_to(epsib)
            call boundaries_set(gr%der%boundaries, gr, st%group%psib(ib, ik))

            if (allocated(hm%hm_base%phase)) then
              call hamiltonian_elec_base_phase(hm%hm_base, gr, gr%np_part, &
                conjugate = .false., psib = epsib, src = st%group%psib(ib, ik))
            else
              call st%group%psib(ib, ik)%copy_data_to(gr%np_part, epsib)
            end if

            ! this now takes non-orthogonal axis into account
            do idir = 1, gr%der%dim
              call epsib%copy_to(commpsib(idir))
            end do
            call zderivatives_batch_grad(gr%der, epsib, commpsib, set_bc=.false.)

            call zhamiltonian_elec_base_nlocal_position_commutator(hm%hm_base, gr, st%d, &
              gr%der%boundaries%spiral, epsib, commpsib)

            call zlda_u_commute_r(hm%lda_u, gr, space, st%d, namespace, epsib, commpsib)

            call current_batch_accumulate(st, gr%der, ik, ib, epsib, commpsib)

            do idir = 1, gr%der%dim
              call commpsib(idir)%end()
            end do

            call epsib%end()
            call st%group%psib(ib, ik)%do_unpack(copy = .false.)

          end do
        end do

      else

        ! use the slow non-packed version

        do ik = st%d%kpt%start, st%d%kpt%end
          ispin = st%d%get_spin_index(ik)
          do ist = st%st_start, st%st_end

            ww = st%d%kweights(ik)*st%occ(ist, ik)
            if (abs(ww) <= M_EPSILON) cycle

            call states_elec_get_state(st, gr, ist, ik, psi)

            do idim = 1, st%d%dim
              call boundaries_set(gr%der%boundaries, gr, psi(:, idim))
            end do

            if (allocated(hm%hm_base%phase)) then
              call states_elec_set_phase(st%d, psi, hm%hm_base%phase(1:gr%np_part, ik), gr%np_part, .false.)
            end if

            do idim = 1, st%d%dim
              call zderivatives_grad(gr%der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
            end do

            if (this%method == CURRENT_GRADIENT_CORR) then
              !A nonlocal contribution from the MGGA potential must be included
              !This must be done first, as this is like a position-dependent mass
              if (family_is_mgga_with_exc(hm%xc)) then
                do idim = 1, st%d%dim
                  do idir = 1, gr%der%dim
                    !$omp parallel do
                    do ip = 1, gr%np
                      gpsi(ip, idir, idim) = (M_ONE+M_TWO*hm%vtau(ip,ispin))*gpsi(ip, idir, idim)
                    end do
                    !$omp end parallel do
                  end do
                end do
              end if

              !A nonlocal contribution from the pseudopotential must be included
              call zprojector_commute_r_allatoms_alldir(hm%ep%proj, hm%ions, gr, st%d%dim, &
                gr%der%boundaries, ik, psi, gpsi)
              !A nonlocal contribution from the scissor must be included
              if (hm%scissor%apply) then
                call scissor_commute_r(hm%scissor, gr, ik, psi, gpsi)
              end if

              call zlda_u_commute_r_single(hm%lda_u, gr, space, st%d, namespace, ist, ik, &
                psi, gpsi, allocated(hm%hm_base%phase))

              call zexchange_operator_commute_r(hm%exxop, namespace, gr, st%d, ik, psi, gpsi)

            end if

            if (st%d%ispin /= SPINORS) then
              do idir = 1, gr%der%dim
                !$omp parallel do
                do ip = 1, gr%np
                  st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) + &
                    ww*aimag(conjg(psi(ip, 1))*gpsi(ip, idir, 1))
                end do
                !$omp end parallel do
              end do
            else
              do idir = 1, gr%der%dim
                !$omp parallel do  private(c_tmp)
                do ip = 1, gr%np
                  st%current_para(ip, idir, 1) = st%current_para(ip, idir, 1) + &
                    ww*aimag(conjg(psi(ip, 1))*gpsi(ip, idir, 1))
                  st%current_para(ip, idir, 2) = st%current_para(ip, idir, 2) + &
                    ww*aimag(conjg(psi(ip, 2))*gpsi(ip, idir, 2))
                  c_tmp = -M_HALF*M_ZI*(conjg(psi(ip, 2))*gpsi(ip, idir, 1) - psi(ip, 1)*conjg(gpsi(ip, idir, 2)))
                  st%current_para(ip, idir, 3) = st%current_para(ip, idir, 3) + ww*TOFLOAT(c_tmp)
                  st%current_para(ip, idir, 4) = st%current_para(ip, idir, 4) + ww*aimag(c_tmp)
                end do
                !$omp end parallel do
              end do
            end if

          end do
        end do

      end if

    case default

      ASSERT(.false.)

    end select

    if (st%d%ispin /= SPINORS) then
      !We sum the current over k-points
      do ik = st%d%kpt%start, st%d%kpt%end
        ispin = st%d%get_spin_index(ik)
        call lalg_axpy(gr%np, gr%der%dim, M_ONE, st%current_kpt(:, :, ik), st%current_para(:, :, ispin))
      end do
    end if

    if (st%parallel_in_states .or. st%d%kpt%parallel) then
      call comm_allreduce(st%st_kpt_mpi_grp, st%current_para, dim = (/gr%np, gr%der%dim, st%d%nspin/))
    end if

    if (st%symmetrize_density) then
      do ispin = 1, st%d%nspin
        call dgrid_symmetrize_vector_field(gr, st%current_para(:, :, ispin), suppress_warning = .true.)
      end do
    end if

    SAFE_DEALLOCATE_A(gpsi)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(rhpsi)
    SAFE_DEALLOCATE_A(rpsi)
    SAFE_DEALLOCATE_A(hrpsi)
    SAFE_DEALLOCATE_A(commpsib)

    call profiling_out(prof)
    POP_SUB(current_calculate_para)

  end subroutine current_calculate_para


  ! ---------------------------------------------------------
  ! Calculate the current matrix element between two states
  ! I_{ij}(t) = <i| J(t) |j>
  ! This is used only in the floquet_observables utility and
  ! is highly experimental

  subroutine current_calculate_mel(der, hm, psi_i, psi_j, ik,  cmel)
    type(derivatives_t),  intent(inout)    :: der
    type(hamiltonian_elec_t),  intent(in)  :: hm
    CMPLX,                intent(in)       :: psi_i(:,:)
    CMPLX,                intent(in)       :: psi_j(:,:)
    integer,              intent(in)       :: ik
    CMPLX,                intent(out)      :: cmel(:,:) ! the current vector cmel(1:der%dim, 1:st%d%nspin)

    integer ::  idir, idim, ip, ispin
    CMPLX, allocatable :: gpsi_j(:, :, :), ppsi_j(:,:),  gpsi_i(:, :, :), ppsi_i(:,:)

    PUSH_SUB(current_calculate_mel)

    SAFE_ALLOCATE(gpsi_i(1:der%mesh%np, 1:der%dim, 1:hm%d%dim))
    SAFE_ALLOCATE(ppsi_i(1:der%mesh%np_part,1:hm%d%dim))
    SAFE_ALLOCATE(gpsi_j(1:der%mesh%np, 1:der%dim, 1:hm%d%dim))
    SAFE_ALLOCATE(ppsi_j(1:der%mesh%np_part,1:hm%d%dim))

    cmel = M_z0

    ispin = hm%d%get_spin_index(ik)
    ppsi_i(:,:) = M_z0
    ppsi_i(1:der%mesh%np,:) = psi_i(1:der%mesh%np,:)
    ppsi_j(:,:) = M_z0
    ppsi_j(1:der%mesh%np,:) = psi_j(1:der%mesh%np,:)


    do idim = 1, hm%d%dim
      call boundaries_set(der%boundaries, der%mesh, ppsi_i(:, idim))
      call boundaries_set(der%boundaries, der%mesh, ppsi_j(:, idim))
    end do

    if (allocated(hm%hm_base%phase)) then
      ! Apply the phase that contains both the k-point and vector-potential terms.
      do idim = 1, hm%d%dim
        !$omp parallel do
        do ip = 1, der%mesh%np_part
          ppsi_i(ip, idim) = hm%hm_base%phase(ip, ik)*ppsi_i(ip, idim)
          ppsi_j(ip, idim) = hm%hm_base%phase(ip, ik)*ppsi_j(ip, idim)
        end do
        !$omp end parallel do
      end do
    end if

    do idim = 1, hm%d%dim
      call zderivatives_grad(der, ppsi_i(:, idim), gpsi_i(:, :, idim), set_bc = .false.)
      call zderivatives_grad(der, ppsi_j(:, idim), gpsi_j(:, :, idim), set_bc = .false.)
    end do

    !A nonlocal contribution from the MGGA potential must be included
    !This must be done first, as this is like a position-dependent mass
    if (family_is_mgga_with_exc(hm%xc)) then
      do idim = 1, hm%d%dim
        do idir = 1, der%dim
          !$omp parallel do
          do ip = 1, der%mesh%np
            gpsi_i(ip, idir, idim) = (M_ONE+M_TWO*hm%vtau(ip,ispin))*gpsi_i(ip, idir, idim)
            gpsi_j(ip, idir, idim) = (M_ONE+M_TWO*hm%vtau(ip,ispin))*gpsi_j(ip, idir, idim)
          end do
          !$omp end parallel do
        end do
      end do

      !A nonlocal contribution from the pseudopotential must be included
      call zprojector_commute_r_allatoms_alldir(hm%ep%proj, hm%ions, der%mesh, hm%d%dim, &
        der%boundaries, ik, ppsi_i, gpsi_i)
      call zprojector_commute_r_allatoms_alldir(hm%ep%proj, hm%ions, der%mesh, hm%d%dim, &
        der%boundaries, ik, ppsi_j, gpsi_j)
      !A nonlocal contribution from the scissor must be included
      if (hm%scissor%apply) then
        call scissor_commute_r(hm%scissor, der%mesh, ik, ppsi_i, gpsi_i)
        call scissor_commute_r(hm%scissor, der%mesh, ik, ppsi_j, gpsi_j)
      end if

    end if


    do idir = 1, der%dim

      do idim = 1, hm%d%dim

        cmel(idir,ispin) = M_zI * zmf_dotp(der%mesh, psi_i(:, idim), gpsi_j(:, idir,idim), reduce = .false.)
        cmel(idir,ispin) = cmel(idir,ispin) - M_zI * zmf_dotp(der%mesh, gpsi_i(:, idir, idim), psi_j(:, idim), reduce = .false.)

      end do
    end do

    if (der%mesh%parallel_in_domains) call der%mesh%allreduce(cmel)



    SAFE_DEALLOCATE_A(gpsi_i)
    SAFE_DEALLOCATE_A(ppsi_i)
    SAFE_DEALLOCATE_A(gpsi_j)
    SAFE_DEALLOCATE_A(ppsi_j)

    POP_SUB(current_calculate_mel)

  end subroutine current_calculate_mel

  ! ---------------------------------------------------------
  subroutine current_heat_calculate(space, der, hm, st, current)
    type(space_t),            intent(in)    :: space
    type(derivatives_t),      intent(in)    :: der
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(states_elec_t),      intent(in)    :: st
    FLOAT,                    intent(out)   :: current(:, :, :)

    integer :: ik, ist, idir, idim, ip, ispin, ndim
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), g2psi(:, :, :, :)
    CMPLX :: tmp

    PUSH_SUB(current_heat_calculate)

    ASSERT(space%is_periodic())
    ASSERT(st%d%dim == 1)

    ndim = space%dim

    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np_part, 1:ndim, 1:st%d%dim))
    SAFE_ALLOCATE(g2psi(1:der%mesh%np, 1:ndim, 1:ndim, 1:st%d%dim))

    do ip = 1, der%mesh%np
      current(ip, 1:ndim, 1:st%d%nspin) = st%current(ip, 1:ndim, 1:st%d%nspin)*hm%ep%vpsl(ip)
    end do


    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = st%d%get_spin_index(ik)
      do ist = st%st_start, st%st_end

        if (abs(st%d%kweights(ik)*st%occ(ist, ik)) <= M_EPSILON) cycle

        call states_elec_get_state(st, der%mesh, ist, ik, psi)
        do idim = 1, st%d%dim
          call boundaries_set(der%boundaries, der%mesh, psi(:, idim))
        end do

        if (allocated(hm%hm_base%phase)) then
          call states_elec_set_phase(st%d, psi, hm%hm_base%phase(1:der%mesh%np_part, ik), der%mesh%np_part,  conjugate = .false.)
        end if

        do idim = 1, st%d%dim
          call zderivatives_grad(der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
        end do
        do idir = 1, ndim
          if (allocated(hm%hm_base%phase)) then
            call states_elec_set_phase(st%d, gpsi(:, idir, :), hm%hm_base%phase(1:der%mesh%np_part, ik), der%mesh%np, &
              conjugate = .true.)
          end if

          do idim = 1, st%d%dim
            call boundaries_set(der%boundaries, der%mesh, gpsi(:,idir, idim))
          end do

          if (allocated(hm%hm_base%phase)) then
            call states_elec_set_phase(st%d, gpsi(:, idir, :), hm%hm_base%phase(1:der%mesh%np_part, ik), &
              der%mesh%np_part,  conjugate = .false.)
          end if

          do idim = 1, st%d%dim
            call zderivatives_grad(der, gpsi(:, idir, idim), g2psi(:, :, idir, idim), set_bc = .false.)
          end do
        end do
        idim = 1
        do ip = 1, der%mesh%np
          do idir = 1, ndim
            !tmp = sum(conjg(g2psi(ip, idir, 1:ndim, idim))*gpsi(ip, idir, idim)) - sum(conjg(gpsi(ip, 1:ndim, idim))*g2psi(ip, idir, 1:ndim, idim))
            tmp = sum(conjg(g2psi(ip, 1:ndim, idir, idim))*gpsi(ip, 1:ndim, idim)) - &
              sum(conjg(gpsi(ip, 1:ndim, idim))*g2psi(ip, 1:ndim, idir, idim))
            tmp = tmp - conjg(gpsi(ip, idir, idim))*sum(g2psi(ip, 1:ndim, 1:ndim, idim)) + &
              sum(conjg(g2psi(ip, 1:ndim, 1:ndim, idim)))*gpsi(ip, idir, idim)
            current(ip, idir, ispin) = current(ip, idir, ispin) + st%d%kweights(ik)*st%occ(ist, ik)*aimag(tmp)/CNST(8.0)
          end do
        end do
      end do
    end do


    POP_SUB(current_heat_calculate)

  end subroutine current_heat_calculate

end module current_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
