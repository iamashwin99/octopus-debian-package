!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2022 N. Tancogne-Dejean
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

module xc_sic_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use varinfo_oct_m
  use xc_oct_m
  use xc_f03_lib_m
  use xc_oep_oct_m

  implicit none

  private
  public ::                     &
    xc_sic_t,                   &
    xc_sic_init,                &
    xc_sic_calc_adsic,          &
    xc_sic_write_info,          &
    xc_sic_end

  !> the SIC levels
  integer, parameter, public :: &
    SIC_NONE   = 1,     &  !< no self-interaction correction
    SIC_PZ_OEP = 2,     &  !< Perdew-Zunger SIC (OEP way)
    SIC_AMALDI = 3,     &  !< Amaldi correction term
    SIC_ADSIC  = 4         !< Averaged density SIC

  type xc_sic_t
    private
    integer,        public :: level = SIC_NONE  !< what kind of self-interaction correction to apply
    FLOAT,          public :: amaldi_factor
    type(xc_oep_t), public :: oep
  end type xc_sic_t

contains

  ! ---------------------------------------------------------
  subroutine xc_sic_init(sic, namespace, gr, st, mc, space)
    type(xc_sic_t),      intent(out)   :: sic
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(inout) :: gr
    type(states_elec_t), intent(in)    :: st
    type(multicomm_t),   intent(in)    :: mc
    type(space_t),       intent(in)    :: space


    PUSH_SUB(xc_sic_init)

    !%Variable SICCorrection
    !%Type integer
    !%Default sic_none
    !%Section Hamiltonian::XC
    !%Description
    !% This variable controls which form of self-interaction correction to use. Note that
    !% this correction will be applied to the functional chosen by <tt>XCFunctional</tt>.
    !%Option sic_none 1
    !% No self-interaction correction.
    !%Option sic_pz 2
    !% Perdew-Zunger SIC, handled by the OEP technique.
    !%Option sic_amaldi 3
    !% Amaldi correction term.
    !%Option sic_adsic 4
    !% Average-density SIC.
    !% C. Legrand <i>et al.</i>, <i>J. Phys. B</i> <b>35</b>, 1115 (2002).
    !%End
    call parse_variable(namespace, 'SICCorrection', SIC_NONE, sic%level)
    if (.not. varinfo_valid_option('SICCorrection', sic%level)) call messages_input_error(namespace, 'SICCorrection')

    ! check whether we should introduce the Amaldi SIC correction
    sic%amaldi_factor = M_ONE
    if (sic%level == SIC_AMALDI) sic%amaldi_factor = (st%qtot - M_ONE)/st%qtot

    if(sic%level == SIC_PZ_OEP) then
      call xc_oep_init(sic%oep, namespace, gr, st, mc, space, oep_type = OEP_TYPE_SIC)
    end if

    POP_SUB(xc_sic_init)
  end subroutine xc_sic_init

  ! ---------------------------------------------------------
  subroutine xc_sic_end(sic)
    type(xc_sic_t),  intent(inout) :: sic

    if (sic%level == SIC_NONE) return

    PUSH_SUB(xc_sic_end)

    if(sic%level == SIC_PZ_OEP) call xc_oep_end(sic%oep)

    POP_SUB(xc_sic_end)
  end subroutine xc_sic_end


  ! ---------------------------------------------------------
  subroutine xc_sic_write_info(sic, iunit, namespace)
    type(xc_sic_t),              intent(in) :: sic
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    if (sic%level == SIC_NONE) return

    PUSH_SUB(xc_sic_write_info)

    call messages_print_var_option('SICCorrection', sic%level, iunit=iunit, namespace=namespace)

    POP_SUB(xc_sic_write_info)
  end subroutine xc_sic_write_info

  ! ---------------------------------------------------------
  ! ADSIC potential is:
  ! V_ADSIC[n] = V_ks[n] - (V_h[n/N] + V_xc[n_{up}/N_{up},0] + Vxc(0, n_{dn}/N_{dn}))
  !
  ! E_ADSIC[n] = E - [N E_H[n/N] + N_{up} E_xc[n_{up}/N_{up},0] + N_{dn} Exc(0, n_{dn}/N_{dn})
  !
  ! C. Legrand et al., J. Phys. B: At. Mol. Opt. Phys. 35 (2002) 1115–1128
  subroutine xc_sic_calc_adsic(sic, namespace, space, gr, st, hm, xc, density, total_density, vxc, ex, ec)
    type(xc_sic_t),              intent(in) :: sic
    type(namespace_t),           intent(in) :: namespace
    type(space_t),               intent(in) :: space
    type(grid_t),                intent(in) :: gr
    type(states_elec_t),         intent(in) :: st
    type(hamiltonian_elec_t),    intent(in) :: hm
    type(xc_t),               intent(inout) :: xc
    FLOAT,                       intent(in) :: density(:,:)
    FLOAT,                       intent(in) :: total_density(:)
    FLOAT,                    intent(inout) :: vxc(:,:)
    FLOAT, optional,          intent(inout) :: ex, ec

    integer            :: ip, ispin, ist, ik
    FLOAT, allocatable :: vxc_sic(:,:),  vh_sic(:), rho(:, :), qsp(:)
    FLOAT :: ex_sic, ec_sic

    PUSH_SUB(xc_sic_calc_adsic)

    ASSERT(sic%level == SIC_ADSIC)

    if (st%d%ispin == SPINORS) then
      call messages_not_implemented('ADSIC with non-collinear spin', namespace=namespace)
    end if

    ! We compute here the number of electrons per spin channel
    SAFE_ALLOCATE(qsp(1:2))
    qsp = M_ZERO
    select case (st%d%ispin)
    case (UNPOLARIZED, SPIN_POLARIZED)
      do ist = 1, st%nst
        do ik = 1, st%d%nik
          ispin = st%d%get_spin_index(ik)
          qsp(ispin) = qsp(ispin) + st%occ(ist, ik) * st%d%kweights(ik)
        end do
      end do
    case (SPINORS)
      do ist = 1, st%nst
        do ik = 1, st%d%nik
          qsp(1) = qsp(1) + st%occ(ist, ik) * st%d%kweights(ik)
        end do
      end do
    end select

    SAFE_ALLOCATE(vxc_sic(1:gr%np, 1:2))
    vxc_sic = M_ZERO
    SAFE_ALLOCATE(rho(1:gr%np, 1:2))
    ! We first compute the average xc self-interction error and we substract it
    select case (st%d%ispin)
    case (UNPOLARIZED, SPIN_POLARIZED)
      do ispin = 1, st%d%nspin
        if (abs(qsp(ispin)) <= M_EPSILON) cycle

        rho = M_ZERO
        vxc_sic = M_ZERO

        rho(:, ispin) = density(:, ispin) / qsp(ispin)
        ! This needs always to be called for the spin-polarized case
        if(present(ex) .and. present(ec)) then
          ex_sic = M_ZERO
          ec_sic = M_ZERO
          call xc_get_vxc(gr, xc, st, hm%kpoints, hm%psolver, namespace, space, &
            rho, SPIN_POLARIZED, hm%ions%latt%rcell_volume, vxc_sic, ex = ex_sic, ec = ec_sic)
          ex = ex - ex_sic * qsp(ispin)
          ec = ec - ec_sic * qsp(ispin)
        else
          call xc_get_vxc(gr, xc, st, hm%kpoints, hm%psolver, namespace, space, &
            rho, SPIN_POLARIZED, hm%ions%latt%rcell_volume, vxc_sic)
        end if

        vxc(:, ispin) = vxc(:, ispin) - vxc_sic(:, ispin)
      end do

    case (SPINORS)
      !TODO
    end select
    SAFE_DEALLOCATE_A(qsp)
    SAFE_DEALLOCATE_A(vxc_sic)

    SAFE_ALLOCATE(vh_sic(1:gr%np))
    vh_sic = M_ZERO
    ! We now substract the averaged Hartree self-interaction error
    rho(:, 1) = total_density / st%qtot
    call dpoisson_solve(hm%psolver, namespace, vh_sic, rho(:,1))
    do ip = 1, gr%np
      vxc(ip,:) = vxc(ip,:) - vh_sic(ip)
    end do
    ! Compute the corresponding energy contribution
    if(present(ex)) then
      ex = ex - dmf_dotp(gr, total_density, vh_sic)
    end if


    SAFE_DEALLOCATE_A(vh_sic)
    SAFE_DEALLOCATE_A(rho)

    POP_SUB(xc_sic_calc_adsic)
  end subroutine xc_sic_calc_adsic

end module xc_sic_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
