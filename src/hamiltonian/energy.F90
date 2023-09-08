!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012 M. Oliveira
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

module energy_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m

  implicit none

  private

  public ::         &
    energy_t,       &
    energy_copy

  type energy_t
    ! Components are public by default
    ! Energies
    FLOAT :: total            = M_ZERO !< Total energy
    !!                                    \f[
    !!                                      E = E_{ii} + \sum[{\rm Eigenvalues}] - U + E_x + E_c - \int n v_{xc}
    !!                                    - 1/2 \int n^e v_{pcm}  + 1/2 \int n^n v_{pcm}  - \int n v_U
    !!                                    \f]
    FLOAT :: eigenvalues      = M_ZERO !< \f$ \sum[{\rm Eigenvalues}] \f$
    FLOAT :: exchange         = M_ZERO
    FLOAT :: exchange_hf      = M_ZERO !< Exchange energy for the Hartree-Fock case only
    FLOAT :: correlation      = M_ZERO
    FLOAT :: vdw              = M_ZERO
    FLOAT :: xc_j             = M_ZERO
    FLOAT :: intnvxc          = M_ZERO !< \f$ \int n v_{\rm xc} \f$
    FLOAT :: hartree          = M_ZERO !< Hartree     \f$ U = (1/2)  \int n v_{\rm Hartree} \f$
    FLOAT :: int_ee_pcm       = M_ZERO !< \f$ 1/2 v_{\rm Hartree} q_{pcm_e} \f$ dot product of vectors of dimension n_tesserae
    FLOAT :: int_en_pcm       = M_ZERO !< \f$ 1/2 v_{\rm Hartree} q_{pcm_n} \f$
    FLOAT :: int_ne_pcm       = M_ZERO !< \f$ 1/2 v_n q_{pcm_e} \f$
    FLOAT :: int_nn_pcm       = M_ZERO !< \f$ 1/2 v_n q_{pcm_n} \f$
    FLOAT :: int_e_ext_pcm    = M_ZERO !< \f$ v_{\rm Hartree} *  q_{pcm_ext} \f$
    FLOAT :: int_n_ext_pcm    = M_ZERO !< \f$ v_n * q_{pcm_ext} \f$
    FLOAT :: pcm_corr         = M_ZERO !< \f$ \int [n (v_{e_rs} + v_{n_rs})] \f$
    FLOAT :: kinetic          = M_ZERO !< Kinetic energy of the non-interacting (KS) system of electrons
    FLOAT :: extern           = M_ZERO !< External     \f$ V = <\Phi|V|\Phi> = \int n v  \f$ (if no non-local pseudos exist)
    FLOAT :: extern_local     = M_ZERO !< The local part of the external energy ( \f$ \int n v \f$ )
    FLOAT :: extern_non_local = M_ZERO !< The part of the external energy coming from the non-local part of the pseudos
    FLOAT :: entropy          = M_ZERO
    FLOAT :: ts               = M_ZERO !< TS
    FLOAT :: berry            = M_ZERO !< Berry energy correction = \f$ -\mu  E - <V_{\rm berry}> \f$
    FLOAT :: delta_xc         = M_ZERO !< the XC derivative discontinuity
    FLOAT :: dft_u            = M_ZERO !< DFT+U contribution
    FLOAT :: int_dft_u        = M_ZERO !< \f$ \int n v_U \f$
    FLOAT :: intnvstatic      = M_ZERO !< \f$ \int n v_{\rm static} \f$ (static electric field)
    FLOAT :: photon_exchange  = M_ZERO
  end type energy_t

contains

  subroutine energy_copy(ein, eout)
    type(energy_t), intent(in)  :: ein
    type(energy_t), intent(out) :: eout

    PUSH_SUB(energy_copy)

    eout%total        = ein%total
    eout%eigenvalues  = ein%eigenvalues
    eout%exchange     = ein%exchange
    eout%exchange_hf  = ein%exchange_hf
    eout%correlation  = ein%correlation
    eout%vdw          = ein%vdw
    eout%xc_j         = ein%xc_j
    eout%intnvxc      = ein%intnvxc
    eout%hartree      = ein%hartree
    eout%int_ee_pcm   = ein%int_ee_pcm
    eout%int_en_pcm   = ein%int_en_pcm
    eout%int_nn_pcm   = ein%int_nn_pcm
    eout%int_ne_pcm   = ein%int_ne_pcm
    eout%int_e_ext_pcm   = ein%int_e_ext_pcm
    eout%int_n_ext_pcm   = ein%int_n_ext_pcm
    eout%pcm_corr     = ein%pcm_corr
    eout%kinetic      = ein%kinetic
    eout%extern       = ein%extern
    eout%extern_local = ein%extern_local
    eout%extern_non_local = ein%extern_non_local
    eout%entropy      = ein%entropy
    eout%ts           = ein%ts
    eout%berry        = ein%berry
    eout%delta_xc     = ein%delta_xc
    eout%dft_u        = ein%dft_u
    eout%int_dft_u    = ein%int_dft_u
    eout%intnvstatic  = ein%intnvstatic
    eout%photon_exchange  = ein%photon_exchange

    POP_SUB(energy_copy)
  end subroutine energy_copy

end module energy_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
