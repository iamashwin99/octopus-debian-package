!! Copyright (C) 2002-2010 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

module density_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use iso_c_binding
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use kpoints_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use smear_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use types_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                           &
    density_calc_t,                   &
    density_calc_init,                &
    density_calc_accumulate,          &
    density_calc_end,                 &
    density_calc,                     &
    states_elec_freeze_orbitals,      &
    states_elec_total_density,        &
    ddensity_accumulate_grad,         &
    zdensity_accumulate_grad,         &
    states_elec_freeze_redistribute_states,&
    states_elec_freeze_adjust_qtot

  type density_calc_t
    private
    FLOAT,                pointer :: density(:, :)
    type(states_elec_t),  pointer :: st
    type(grid_t),         pointer :: gr
    type(accel_mem_t)            :: buff_density
    integer(i8)                   :: pnp
    logical                       :: packed
  end type density_calc_t

contains

  subroutine density_calc_init(this, st, gr, density)
    type(density_calc_t),           intent(out)   :: this
    type(states_elec_t),  target,   intent(in)    :: st
    type(grid_t),         target,   intent(in)    :: gr
    FLOAT,                target,   intent(out)   :: density(:, :)

    integer :: ip, ispin
    logical :: correct_size

    PUSH_SUB(density_calc_init)

    this%st => st
    this%gr => gr

    this%density => density
    !$omp parallel private(ip)
    do ispin = 1, st%d%nspin
      !$omp do
      do ip = 1, gr%np
        this%density(ip, ispin) = M_ZERO
      end do
      !$omp end do nowait
    end do
    !$omp end parallel

    this%packed = .false.

    correct_size = ubound(this%density, dim = 1) == this%gr%np .or. &
      ubound(this%density, dim = 1) == this%gr%np_part
    ASSERT(correct_size)

    POP_SUB(density_calc_init)
  end subroutine density_calc_init

  ! ---------------------------------------------------

  subroutine density_calc_pack(this)
    type(density_calc_t),           intent(inout)   :: this

    PUSH_SUB(density_calc_pack)

    this%packed = .true.
    this%pnp = accel_padded_size(this%gr%np)
    call accel_create_buffer(this%buff_density, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, this%pnp*this%st%d%nspin)

    ! set to zero
    call accel_set_buffer_to_zero(this%buff_density, TYPE_FLOAT, this%pnp*this%st%d%nspin)

    POP_SUB(density_calc_pack)
  end subroutine density_calc_pack

  ! ---------------------------------------------------

  subroutine density_calc_state(this, psib, istin)
    type(density_calc_t),         intent(inout) :: this
    type(wfs_elec_t),             intent(in)    :: psib
    integer,                      intent(in)    :: istin

    integer :: ist, ip, ispin, istin_, wgsize
    CMPLX   :: term, psi1, psi2
    FLOAT, allocatable :: weight(:)
    type(profile_t), save :: prof
    type(accel_mem_t) :: buff_weight
    type(accel_kernel_t), pointer :: kernel

    PUSH_SUB(density_calc_state)
    call profiling_in(prof, "CALC_STATE_DENSITY")

    ispin = this%st%d%get_spin_index(psib%ik)

    istin_= 0

    SAFE_ALLOCATE(weight(1:psib%nst))
    do ist = 1, psib%nst
      weight(ist) = this%st%d%kweights(psib%ik)*this%st%occ(psib%ist(ist), psib%ik)
      if (psib%ist(ist) == istin) istin_ = ist
    end do

    if (abs(weight(istin_)) <= M_EPSILON) then
      return
      POP_SUB(density_calc_state)
    end if

    !fix ist for the remaining code
    ist = istin_


    select case (psib%status())
    case (BATCH_NOT_PACKED)
      select case (this%st%d%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        if (states_are_real(this%st)) then
          !$omp parallel do simd schedule(static)
          do ip = 1, this%gr%np
            this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%dff(ip, 1, ist)**2
          end do
        else
          !$omp parallel do schedule(static)
          do ip = 1, this%gr%np
            this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
              TOFLOAT(conjg(psib%zff(ip, 1, ist))*psib%zff(ip, 1, ist))
          end do
        end if
      case (SPINORS)
        !$omp parallel do schedule(static) private(psi1, psi2, term)
        do ip = 1, this%gr%np
          psi1 = psib%zff(ip, 1, ist)
          psi2 = psib%zff(ip, 2, ist)
          this%density(ip, 1) = this%density(ip, 1) + weight(ist)*TOFLOAT(conjg(psi1)*psi1)
          this%density(ip, 2) = this%density(ip, 2) + weight(ist)*TOFLOAT(conjg(psi2)*psi2)
          term = weight(ist)*psi1*conjg(psi2)
          this%density(ip, 3) = this%density(ip, 3) + TOFLOAT(term)
          this%density(ip, 4) = this%density(ip, 4) + aimag(term)
        end do
      end select

    case (BATCH_PACKED)

      select case (this%st%d%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        if (states_are_real(this%st)) then
          !$omp parallel do schedule(static)
          do ip = 1, this%gr%np
            this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%dff_pack(ist, ip)**2
          end do
        else

          !$omp parallel do schedule(static)
          do ip = 1, this%gr%np
            this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
              TOFLOAT(conjg(psib%zff_pack(ist, ip))*psib%zff_pack(ist, ip))
          end do
        end if
      case (SPINORS)
        ASSERT(mod(psib%nst_linear, 2) == 0)
        !$omp parallel do schedule(static) private(ist, psi1, psi2, term)
        do ip = 1, this%gr%np
          psi1 = psib%zff_pack(2*ist - 1, ip)
          psi2 = psib%zff_pack(2*ist,     ip)
          term = weight(ist)*psi1*conjg(psi2)

          this%density(ip, 1) = this%density(ip, 1) + weight(ist)*TOFLOAT(conjg(psi1)*psi1)
          this%density(ip, 2) = this%density(ip, 2) + weight(ist)*TOFLOAT(conjg(psi2)*psi2)
          this%density(ip, 3) = this%density(ip, 3) + TOFLOAT(term)
          this%density(ip, 4) = this%density(ip, 4) + aimag(term)
        end do
      end select

    case (BATCH_DEVICE_PACKED)
      if (.not. this%packed) call density_calc_pack(this)

      call accel_create_buffer(buff_weight, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, psib%nst)
      call accel_write_buffer(buff_weight, psib%nst, weight)

      select case (this%st%d%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        if (states_are_real(this%st)) then
          kernel => kernel_density_real
        else
          kernel => kernel_density_complex
        end if

        call accel_set_kernel_arg(kernel, 0, psib%nst)
        call accel_set_kernel_arg(kernel, 1, this%gr%np)
        call accel_set_kernel_arg(kernel, 2, int(this%pnp*(ispin - 1), i4))
        call accel_set_kernel_arg(kernel, 3, buff_weight)
        call accel_set_kernel_arg(kernel, 4, psib%ff_device)
        call accel_set_kernel_arg(kernel, 5, int(log2(psib%pack_size(1)), i4))
        call accel_set_kernel_arg(kernel, 6, this%buff_density)

      case (SPINORS)
        kernel => kernel_density_spinors

        call accel_set_kernel_arg(kernel, 0, psib%nst)
        call accel_set_kernel_arg(kernel, 1, this%gr%np)
        call accel_set_kernel_arg(kernel, 2, int(this%pnp, i4))
        call accel_set_kernel_arg(kernel, 3, buff_weight)
        call accel_set_kernel_arg(kernel, 4, psib%ff_device)
        call accel_set_kernel_arg(kernel, 5, int(log2(psib%pack_size(1)), i4))
        call accel_set_kernel_arg(kernel, 6, this%buff_density)
      end select

      wgsize = accel_kernel_workgroup_size(kernel)

      call accel_kernel_run(kernel, (/pad(this%gr%np, wgsize)/), (/wgsize/))

      call accel_finish()

      call accel_release_buffer(buff_weight)

    end select

    SAFE_DEALLOCATE_A(weight)

    call profiling_out(prof)

    POP_SUB(density_calc_state)
  end subroutine density_calc_state

  ! ---------------------------------------------------

  subroutine density_calc_accumulate(this, psib)
    type(density_calc_t),         intent(inout) :: this
    type(wfs_elec_t),             intent(in)    :: psib

    integer :: ist, ip, ispin, wgsize
    CMPLX   :: term, psi1, psi2
    FLOAT, allocatable :: weight(:)
    type(profile_t), save :: prof
    type(accel_mem_t) :: buff_weight
    type(accel_kernel_t), pointer :: kernel

    PUSH_SUB(density_calc_accumulate)
    call profiling_in(prof, "CALC_DENSITY")

    ispin = this%st%d%get_spin_index(psib%ik)

    SAFE_ALLOCATE(weight(1:psib%nst))
    do ist = 1, psib%nst
      weight(ist) = this%st%d%kweights(psib%ik)*this%st%occ(psib%ist(ist), psib%ik)
    end do

    select case (psib%status())
    case (BATCH_NOT_PACKED)
      select case (this%st%d%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        if (states_are_real(this%st)) then
          do ist = 1, psib%nst
            if (abs(weight(ist)) <= M_EPSILON) cycle
            !$omp parallel do simd schedule(static)
            do ip = 1, this%gr%np
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%dff(ip, 1, ist)**2
            end do
          end do
        else
          do ist = 1, psib%nst
            if (abs(weight(ist)) <= M_EPSILON) cycle
            !$omp parallel do schedule(static)
            do ip = 1, this%gr%np
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
                TOFLOAT(conjg(psib%zff(ip, 1, ist))*psib%zff(ip, 1, ist))
            end do
          end do
        end if
      case (SPINORS)
        do ist = 1, psib%nst
          if (abs(weight(ist)) <= M_EPSILON) cycle
          !$omp parallel do schedule(static) private(psi1, psi2, term)
          do ip = 1, this%gr%np
            psi1 = psib%zff(ip, 1, ist)
            psi2 = psib%zff(ip, 2, ist)
            this%density(ip, 1) = this%density(ip, 1) + weight(ist)*TOFLOAT(conjg(psi1)*psi1)
            this%density(ip, 2) = this%density(ip, 2) + weight(ist)*TOFLOAT(conjg(psi2)*psi2)
            term = weight(ist)*psi1*conjg(psi2)
            this%density(ip, 3) = this%density(ip, 3) + TOFLOAT(term)
            this%density(ip, 4) = this%density(ip, 4) + aimag(term)
          end do
        end do
      end select

    case (BATCH_PACKED)

      select case (this%st%d%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        if (states_are_real(this%st)) then
          !$omp parallel do schedule(static)
          do ip = 1, this%gr%np
            do ist = 1, psib%nst
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%dff_pack(ist, ip)**2
            end do
          end do
        else
          !$omp parallel do schedule(static)
          do ip = 1, this%gr%np
            do ist = 1, psib%nst
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
                TOFLOAT(conjg(psib%zff_pack(ist, ip))*psib%zff_pack(ist, ip))
            end do
          end do
        end if
      case (SPINORS)
        ASSERT(mod(psib%nst_linear, 2) == 0)
        !$omp parallel do schedule(static) private(ist, psi1, psi2, term)
        do ip = 1, this%gr%np
          do ist = 1, psib%nst
            psi1 = psib%zff_pack(2*ist - 1, ip)
            psi2 = psib%zff_pack(2*ist,     ip)
            term = weight(ist)*psi1*conjg(psi2)

            this%density(ip, 1) = this%density(ip, 1) + weight(ist)*TOFLOAT(conjg(psi1)*psi1)
            this%density(ip, 2) = this%density(ip, 2) + weight(ist)*TOFLOAT(conjg(psi2)*psi2)
            this%density(ip, 3) = this%density(ip, 3) + TOFLOAT(term)
            this%density(ip, 4) = this%density(ip, 4) + aimag(term)
          end do
        end do
      end select

    case (BATCH_DEVICE_PACKED)
      if (.not. this%packed) call density_calc_pack(this)

      call accel_create_buffer(buff_weight, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, psib%nst)
      call accel_write_buffer(buff_weight, psib%nst, weight)

      select case (this%st%d%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        if (states_are_real(this%st)) then
          kernel => kernel_density_real
        else
          kernel => kernel_density_complex
        end if

        call accel_set_kernel_arg(kernel, 0, psib%nst)
        call accel_set_kernel_arg(kernel, 1, this%gr%np)
        call accel_set_kernel_arg(kernel, 2, int(this%pnp, i4)*(ispin - 1))
        call accel_set_kernel_arg(kernel, 3, buff_weight)
        call accel_set_kernel_arg(kernel, 4, psib%ff_device)
        call accel_set_kernel_arg(kernel, 5, int(log2(psib%pack_size(1)), i4))
        call accel_set_kernel_arg(kernel, 6, this%buff_density)

      case (SPINORS)
        kernel => kernel_density_spinors

        call accel_set_kernel_arg(kernel, 0, psib%nst)
        call accel_set_kernel_arg(kernel, 1, this%gr%np)
        call accel_set_kernel_arg(kernel, 2, int(this%pnp, i4))
        call accel_set_kernel_arg(kernel, 3, buff_weight)
        call accel_set_kernel_arg(kernel, 4, psib%ff_device)
        call accel_set_kernel_arg(kernel, 5, int(log2(psib%pack_size(1)), i4))
        call accel_set_kernel_arg(kernel, 6, this%buff_density)
      end select

      wgsize = accel_kernel_workgroup_size(kernel)

      call accel_kernel_run(kernel, (/pad(this%gr%np, wgsize)/), (/wgsize/))

      call accel_finish()

      call accel_release_buffer(buff_weight)

    end select

    SAFE_DEALLOCATE_A(weight)

    call profiling_out(prof)

    POP_SUB(density_calc_accumulate)
  end subroutine density_calc_accumulate

  ! ---------------------------------------------------

  subroutine density_calc_end(this, allreduce, symmetrize)
    type(density_calc_t), intent(inout) :: this
    logical, optional,    intent(in)    :: allreduce
    logical, optional,    intent(in)    :: symmetrize

    integer :: ispin, ip
    type(profile_t), save :: reduce_prof
    FLOAT, allocatable :: tmpdensity(:)

    PUSH_SUB(density_calc_end)

    if (this%packed) then
      SAFE_ALLOCATE(tmpdensity(1:this%gr%np_part))

      ! the density is in device memory
      do ispin = 1, this%st%d%nspin
        call accel_read_buffer(this%buff_density, int(this%gr%np, i8), tmpdensity, offset = (ispin - 1)*this%pnp)

        !$omp parallel do
        do ip = 1, this%gr%np
          this%density(ip, ispin) = this%density(ip, ispin) + tmpdensity(ip)
        end do
      end do

      this%packed = .false.
      call accel_release_buffer(this%buff_density)
      SAFE_DEALLOCATE_A(tmpdensity)
    end if

    ! reduce over states and k-points
    if ((this%st%parallel_in_states .or. this%st%d%kpt%parallel) .and. optional_default(allreduce, .true.)) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      call comm_allreduce(this%st%st_kpt_mpi_grp, this%density, dim = (/this%gr%np, this%st%d%nspin/))
      call profiling_out(reduce_prof)
    end if

    if (this%st%symmetrize_density .and. optional_default(symmetrize, .true.)) then
      do ispin = 1, this%st%d%nspin
        call dgrid_symmetrize_scalar_field(this%gr, this%density(:, ispin))
      end do
    end if

    POP_SUB(density_calc_end)
  end subroutine density_calc_end


  ! ---------------------------------------------------------
  !> Computes the density from the orbitals in st.
  subroutine density_calc(st, gr, density, istin)
    type(states_elec_t),     intent(in)  :: st
    type(grid_t),            intent(in)  :: gr
    FLOAT,                   intent(out) :: density(:, :)
    integer, optional,       intent(in)  :: istin

    integer :: ik, ib
    type(density_calc_t) :: dens_calc

    PUSH_SUB(density_calc)

    ASSERT(ubound(density, dim = 1) == gr%np .or. ubound(density, dim = 1) == gr%np_part)

    call density_calc_init(dens_calc, st, gr, density)

    if (.not. present(istin)) then
      ! ordinary density accumulate
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))
        end do
      end do
    else
      ! state-resolved density
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          if (st%group%psib(ib, ik)%ist(1) <= istin .and. st%group%psib(ib, ik)%ist(st%group%psib(ib, ik)%nst) >= istin) then
            call density_calc_state(dens_calc, st%group%psib(ib, ik), istin)
          end if
        end do
      end do

    end if


    call density_calc_end(dens_calc, allreduce = .not. present(istin))

    POP_SUB(density_calc)
  end subroutine density_calc

  ! ---------------------------------------------------------

  subroutine states_elec_freeze_orbitals(st, namespace, space, gr, mc, kpoints, n, family_is_mgga)
    type(states_elec_t), intent(inout) :: st
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    type(grid_t),        intent(in)    :: gr
    type(multicomm_t),   intent(in)    :: mc
    type(kpoints_t),     intent(in)    :: kpoints
    integer,             intent(in)    :: n
    logical,             intent(in)    :: family_is_mgga

    integer :: ist, ik, ib, nblock
    integer :: nodeto, nodefr
    type(states_elec_t) :: staux
    CMPLX, allocatable :: psi(:, :, :), rec_buffer(:,:)
    type(wfs_elec_t)  :: psib
    type(density_calc_t) :: dens_calc

    PUSH_SUB(states_elec_freeze_orbitals)

    if (n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call messages_fatal(2, namespace=namespace)
    end if

    ASSERT(states_are_complex(st))

    if (.not. allocated(st%frozen_rho)) then
      SAFE_ALLOCATE(st%frozen_rho(1:gr%np, 1:st%d%nspin))
    end if

    call density_calc_init(dens_calc, st, gr, st%frozen_rho)

    do ik = st%d%kpt%start, st%d%kpt%end
      if (n < st%st_start) cycle

      do ib =  st%group%block_start, st%group%block_end
        !We can use the full batch
        if (states_elec_block_max(st, ib) <= n) then

          call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))
          if (states_elec_block_max(st, ib) == n) exit

        else !Here we only use a part of this batch

          nblock = n - states_elec_block_min(st, ib) + 1

          SAFE_ALLOCATE(psi(1:gr%np, 1:st%d%dim, 1:nblock))

          do ist = 1, nblock
            call states_elec_get_state(st, gr, states_elec_block_min(st, ib) + ist - 1, ik, psi(:, :, ist))
          end do

          call wfs_elec_init(psib, st%d%dim, states_elec_block_min(st, ib), n, psi, ik)
          call density_calc_accumulate(dens_calc, psib)
          call psib%end()
          SAFE_DEALLOCATE_A(psi)

          exit

        end if

      end do
    end do

    call density_calc_end(dens_calc)

    if (family_is_mgga) then
      if (.not. allocated(st%frozen_tau)) then
        SAFE_ALLOCATE(st%frozen_tau(1:gr%np, 1:st%d%nspin))
      end if
      if (.not. allocated(st%frozen_gdens)) then
        SAFE_ALLOCATE(st%frozen_gdens(1:gr%np, 1:space%dim, 1:st%d%nspin))
      end if
      if (.not. allocated(st%frozen_ldens)) then
        SAFE_ALLOCATE(st%frozen_ldens(1:gr%np, 1:st%d%nspin))
      end if

      call states_elec_calc_quantities(gr, st, kpoints, .true., kinetic_energy_density = st%frozen_tau, &
        density_gradient = st%frozen_gdens, density_laplacian = st%frozen_ldens, st_end = n)
    end if


    call states_elec_copy(staux, st)

    call states_elec_freeze_redistribute_states(st, namespace, gr, mc, n)

    SAFE_ALLOCATE(psi(1:gr%np, 1:st%d%dim, 1))
    SAFE_ALLOCATE(rec_buffer(1:gr%np, 1:st%d%dim))

    do ik = st%d%kpt%start, st%d%kpt%end
      ! loop over all frozen states
      do ist = 1, st%nst
        ! new distribution: frozen states
        nodeto = st%node(ist)
        ! old distribution: full states
        nodefr = staux%node(ist+n)

        if (st%mpi_grp%rank == nodeto .and. st%mpi_grp%rank == nodefr) then
          ! do a simple copy on the same rank, also for serial version
          call states_elec_get_state(staux, gr, ist+n, ik, psi(:, :, 1))
          call states_elec_set_state(st, gr, ist, ik, psi(:, :, 1))
        else
          ! we can use blocking send/recv calls here because we always have different ranks
          if (st%mpi_grp%rank == nodefr) then
            call states_elec_get_state(staux, gr, ist+n, ik, psi(:, :, 1))
            call st%mpi_grp%send(psi(1, 1, 1), gr%np*st%d%dim, MPI_CMPLX, nodeto, ist)
          end if
          if (st%mpi_grp%rank == nodeto) then
            call st%mpi_grp%recv(rec_buffer(1, 1), gr%np*st%d%dim, MPI_CMPLX, nodefr, ist)
            call states_elec_set_state(st, gr, ist, ik, rec_buffer(:, :))
          end if
        end if
      end do
    end do

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(rec_buffer)

    do ik = 1, st%d%nik
      do ist = 1, st%nst
        st%occ(ist, ik) = staux%occ(n+ist, ik)
        st%eigenval(ist, ik) = staux%eigenval(n+ist, ik)
      end do
    end do

    call states_elec_end(staux)
    POP_SUB(states_elec_freeze_orbitals)
  end subroutine states_elec_freeze_orbitals

  ! ---------------------------------------------------------
  subroutine states_elec_freeze_redistribute_states(st, namespace, mesh, mc, nn)
    type(states_elec_t), intent(inout) :: st
    type(namespace_t),   intent(in)    :: namespace
    class(mesh_t),       intent(in)    :: mesh
    type(multicomm_t),   intent(in)    :: mc
    integer,             intent(in)    :: nn

    PUSH_SUB(states_elec_freeze_redistribute_states)

    st%nst = st%nst - nn

    call states_elec_deallocate_wfns(st)
    call states_elec_distribute_nodes(st, namespace, mc)
    call states_elec_allocate_wfns(st, mesh, TYPE_CMPLX)

    SAFE_DEALLOCATE_A(st%eigenval)
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)

    SAFE_DEALLOCATE_A(st%occ)
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    st%occ      = M_ZERO


    POP_SUB(states_elec_freeze_redistribute_states)
  end subroutine states_elec_freeze_redistribute_states

  ! ---------------------------------------------------------
  subroutine states_elec_freeze_adjust_qtot(st)
    type(states_elec_t), intent(inout) :: st

    integer :: ik, ist

    PUSH_SUB(states_elec_freeze_adjust_qtot)

    ! Change the smearing method by fixing the occupations to
    ! that of the ground-state such that the unfrozen states inherit
    ! those values.
    st%smear%method = SMEAR_FIXED_OCC

    ! Set total charge
    st%qtot = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%qtot = st%qtot + st%occ(ist, ik) * st%d%kweights(ik)
      end do
    end do

    if (st%parallel_in_states .or. st%d%kpt%parallel) then
      call comm_allreduce(st%st_kpt_mpi_grp, st%qtot)
    end if


    POP_SUB(states_elec_freeze_adjust_qtot)
  end subroutine states_elec_freeze_adjust_qtot


  ! ---------------------------------------------------------
  !> this routine calculates the total electronic density,
  !! which is the sum of the part coming from the orbitals, the
  !! non-linear core corrections and the frozen orbitals
  subroutine states_elec_total_density(st, mesh, total_rho)
    type(states_elec_t),  intent(in)  :: st
    class(mesh_t),        intent(in)  :: mesh
    FLOAT,                intent(out) :: total_rho(:,:)

    integer :: is, ip

    PUSH_SUB(states_elec_total_density)

    do is = 1, st%d%nspin
      do ip = 1, mesh%np
        total_rho(ip, is) = st%rho(ip, is)
      end do
    end do

    if (allocated(st%rho_core)) then
      !$omp parallel private(ip)
      do is = 1, st%d%spin_channels
        !$omp do
        do ip = 1, mesh%np
          total_rho(ip, is) = total_rho(ip, is) + st%rho_core(ip)/st%d%spin_channels
        end do
        !$omp end do nowait
      end do
      !$omp end parallel
    end if

    ! Add, if it exists, the frozen density from the inner orbitals.
    if (allocated(st%frozen_rho)) then
      !$omp parallel private(ip)
      do is = 1, st%d%nspin
        !$omp do
        do ip = 1, mesh%np
          total_rho(ip, is) = total_rho(ip, is) + st%frozen_rho(ip, is)
        end do
        !$omp end do nowait
      end do
      !$omp end parallel
    end if

    POP_SUB(states_elec_total_density)
  end subroutine states_elec_total_density

#include "undef.F90"
#include "real.F90"
#include "density_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "density_inc.F90"

end module density_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
