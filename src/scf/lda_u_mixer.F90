!! Copyright (C) 2016 N. Tancogne-Dejean
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

module lda_u_mixer_oct_m
  use debug_oct_m
  use global_oct_m
  use lda_u_oct_m
  use messages_oct_m
  use mix_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use types_oct_m

  implicit none

  private

  public ::                     &
    lda_u_mixer_t,              &
    lda_u_mixer_init,           &
    lda_u_mixer_init_auxmixer,  &
    lda_u_mixer_clear,          &
    lda_u_mixer_end,            &
    lda_u_mixer_set_vout,       &
    lda_u_mixer_set_vin,        &
    lda_u_mixer_get_vnew

  type lda_u_mixer_t
    private
    integer :: occsize
    logical :: realstates
    logical :: apply = .false.
    logical :: mixU = .false.

    FLOAT, allocatable :: dtmp_occ(:,:), tmpU(:,:)
    CMPLX, allocatable :: ztmp_occ(:,:)

    type(mixfield_t)  :: mixfield_occ, mixfield_U
  end type lda_u_mixer_t

contains

  ! ---------------------------------------------------------
  subroutine lda_u_mixer_init_auxmixer(this, namespace, mixer, smix, st)
    type(lda_u_t),       intent(in)    :: this
    type(namespace_t),   intent(in)    :: namespace
    type(lda_u_mixer_t), intent(inout) :: mixer
    type(mix_t),         intent(inout) :: smix
    type(states_elec_t), intent(in)    :: st

    integer :: dim1

    if (this%level == DFT_U_NONE) return
    PUSH_SUB(lda_u_mixer_init_auxmixer)

    dim1 = this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets
    if (this%level == DFT_U_ACBN0) then
      dim1 = dim1*2
      if (this%intersite) then
        dim1 = dim1 + 2*this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets*this%maxneighbors
      end if
    end if

    if (states_are_real(st)) then
      call mixfield_init(smix, mixer%mixfield_occ, dim1, 1, 1, mix_d4(smix), TYPE_FLOAT)
      mixer%realstates = .true.
    else
      call mixfield_init(smix, mixer%mixfield_occ, dim1, 1, 1, mix_d4(smix), TYPE_CMPLX)
      mixer%realstates = .false.
    end if
    call mixfield_clear(mix_scheme(smix), mixer%mixfield_occ)
    call mix_add_auxmixfield(namespace, smix, mixer%mixfield_occ)

    if (this%level == DFT_U_ACBN0) then
      call mixfield_init(smix, mixer%mixfield_U, this%norbsets, 1, 1,  mix_d4(smix), TYPE_FLOAT)
      call mixfield_clear(mix_scheme(smix), mixer%mixfield_U)
      call mix_add_auxmixfield(namespace, smix, mixer%mixfield_U)
      mixer%mixU = .true.
    else
      mixer%mixU = .false.
    end if

    POP_SUB(lda_u_mixer_init_auxmixer)
  end subroutine lda_u_mixer_init_auxmixer


  ! ---------------------------------------------------------
  subroutine lda_u_mixer_init(this, mixer, st)
    type(lda_u_t),       intent(in)    :: this
    type(lda_u_mixer_t), intent(inout) :: mixer
    type(states_elec_t), intent(in)    :: st

    if (this%level == DFT_U_NONE) return
    PUSH_SUB(lda_u_mixer_init)

    mixer%apply = .true.

    mixer%occsize = this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets
    if (this%level == DFT_U_ACBN0) then
      mixer%occsize = mixer%occsize*2
      if (this%intersite) then
        mixer%occsize = mixer%occsize + 2*this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets*this%maxneighbors
      end if
    end if

    if (states_are_real(st)) then
      SAFE_ALLOCATE(mixer%dtmp_occ(1:mixer%occsize, 1))
    else
      SAFE_ALLOCATE(mixer%ztmp_occ(1:mixer%occsize, 1))
    end if

    if (this%level == DFT_U_ACBN0) then
      SAFE_ALLOCATE(mixer%tmpU(1:this%norbsets, 1))
    end if

    POP_SUB(lda_u_mixer_init)
  end subroutine lda_u_mixer_init

  ! ---------------------------------------------------------
  subroutine lda_u_mixer_clear(mixer, smix)
    type(lda_u_mixer_t), intent(inout) :: mixer
    type(mix_t),         intent(inout)   :: smix

    if (.not. mixer%apply) return
    PUSH_SUB(lda_u_mixer_clear)

    call mixfield_clear(mix_scheme(smix), mixer%mixfield_occ)
    if (mixer%mixU) call mixfield_clear(mix_scheme(smix), mixer%mixfield_U)

    POP_SUB(lda_u_mixer_clear)
  end subroutine lda_u_mixer_clear

  ! ---------------------------------------------------------
  subroutine lda_u_mixer_end(mixer, smix)
    type(lda_u_mixer_t), intent(inout) :: mixer
    type(mix_t),         intent(inout) :: smix

    if (.not. mixer%apply) return
    PUSH_SUB(lda_u_mixer_end)

    SAFE_DEALLOCATE_A(mixer%dtmp_occ)
    SAFE_DEALLOCATE_A(mixer%ztmp_occ)
    SAFE_DEALLOCATE_A(mixer%tmpU)

    call mixfield_end(smix,mixer%mixfield_occ)
    if (mixer%mixU) call mixfield_end(smix, mixer%mixfield_U)

    POP_SUB(lda_u_mixer_end)
  end subroutine lda_u_mixer_end

  ! ---------------------------------------------------------
  subroutine lda_u_mixer_set_vout(this, mixer)
    type(lda_u_t),       intent(in)    :: this
    type(lda_u_mixer_t), intent(inout) :: mixer

    if (this%level == DFT_U_NONE) return
    PUSH_SUB(lda_u_mixer_set_vout)

    if (mixer%realstates) then
      call dlda_u_get_occupations(this, mixer%dtmp_occ(1:mixer%occsize, 1))
      call mixfield_set_vout(mixer%mixfield_occ, mixer%dtmp_occ)
    else
      call zlda_u_get_occupations(this, mixer%ztmp_occ(1:mixer%occsize, 1))
      call mixfield_set_vout(mixer%mixfield_occ, mixer%ztmp_occ)
    end if

    if (this%level == DFT_U_ACBN0) then
      call lda_u_get_effectiveU(this, mixer%tmpU(1:this%norbsets, 1))
      call mixfield_set_vout(mixer%mixfield_U, mixer%tmpU)
    end if

    POP_SUB(lda_u_mixer_set_vout)
  end subroutine lda_u_mixer_set_vout

  ! ---------------------------------------------------------
  subroutine lda_u_mixer_set_vin(this, mixer)
    type(lda_u_t),       intent(in)    :: this
    type(lda_u_mixer_t), intent(inout) :: mixer

    if (.not. mixer%apply) return
    PUSH_SUB(lda_u_mixer_set_vin)

    if (mixer%realstates) then
      call dlda_u_get_occupations(this, mixer%dtmp_occ(1:mixer%occsize, 1))
      call mixfield_set_vin(mixer%mixfield_occ, mixer%dtmp_occ)
    else
      call zlda_u_get_occupations(this, mixer%ztmp_occ(1:mixer%occsize, 1))
      call mixfield_set_vin(mixer%mixfield_occ, mixer%ztmp_occ)
    end if

    if (this%level == DFT_U_ACBN0) then
      call lda_u_get_effectiveU(this, mixer%tmpU(1:this%norbsets, 1))
      call mixfield_set_vin(mixer%mixfield_U, mixer%tmpU)
    end if

    POP_SUB(lda_u_mixer_set_vin)
  end subroutine lda_u_mixer_set_vin

  ! ---------------------------------------------------------
  subroutine lda_u_mixer_get_vnew(this, mixer, st)
    type(lda_u_t),       intent(inout) :: this
    type(lda_u_mixer_t), intent(inout) :: mixer
    type(states_elec_t), intent(in)    :: st

    if (.not. mixer%apply) return
    PUSH_SUB(lda_u_mixer_get_vnew)

    if (this%level == DFT_U_ACBN0) then
      call mixfield_get_vnew(mixer%mixfield_U, mixer%tmpU)
      call lda_u_set_effectiveU(this, mixer%tmpU(1:this%norbsets, 1))
    end if


    if (mixer%realstates) then
      call mixfield_get_vnew(mixer%mixfield_occ, mixer%dtmp_occ)
      call dlda_u_set_occupations(this, mixer%dtmp_occ(1:mixer%occsize, 1))
      call dlda_u_update_potential(this, st)
    else
      call mixfield_get_vnew(mixer%mixfield_occ, mixer%ztmp_occ)
      call zlda_u_set_occupations(this, mixer%ztmp_occ(1:mixer%occsize, 1))
      call zlda_u_update_potential(this, st)
    end if

    POP_SUB(lda_u_mixer_get_vnew)
  end subroutine lda_u_mixer_get_vnew


end module lda_u_mixer_oct_m
