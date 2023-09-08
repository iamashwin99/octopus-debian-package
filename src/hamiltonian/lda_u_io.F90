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

module lda_u_io_oct_m
  use atomic_orbital_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lalg_basic_oct_m
  use lda_u_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::                             &
    lda_u_write_occupation_matrices, &
    lda_u_write_effectiveU,          &
    lda_u_write_kanamoriU,           &
    lda_u_write_U,                   &
    lda_u_write_V,                   &
    lda_u_write_magnetization,       &
    lda_u_load,                      &
    lda_u_dump

contains

  !> Prints the occupation matrices at the end of the scf calculation.
  subroutine lda_u_write_occupation_matrices(dir, this, st, namespace)
    type(lda_u_t),       intent(in)    :: this
    character(len=*),    intent(in)    :: dir
    type(states_elec_t), intent(in)    :: st
    type(namespace_t),   intent(in)    :: namespace

    integer :: iunit, ios, ispin, im, imp

    PUSH_SUB(lda_u_write_occupation_matrices)

    if (mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
      iunit = io_open(trim(dir) // "/occ_matrices", namespace, action='write')
      write(iunit,'(a)') ' Occupation matrices '

      do ios = 1, this%norbsets
        do ispin = 1,st%d%nspin
          write(iunit,'(a, i4, a, i4)') ' Orbital set ', ios, ' spin ', ispin
          do im = 1, this%orbsets(ios)%norbs
            write(iunit,'(1x)',advance='no')

            if (states_are_real(st)) then
              do imp = 1, this%orbsets(ios)%norbs-1
                write(iunit,'(f14.8)',advance='no') this%dn(im,imp,ispin,ios)
              end do
              write(iunit,'(f14.8)') this%dn(im,this%orbsets(ios)%norbs,ispin,ios)
            else
              do imp = 1, this%orbsets(ios)%norbs-1
                write(iunit,'(f14.8,f14.8)',advance='no') this%zn(im,imp,ispin,ios)
              end do
              write(iunit,'(f14.8,f14.8)') this%zn(im,this%orbsets(ios)%norbs,ispin,ios)
            end if
          end do
        end do !ispin
      end do !iatom
      call io_close(iunit)

      if (this%level == DFT_U_ACBN0) then
        iunit = io_open(trim(dir) // "/renorm_occ_matrices", namespace, action='write')
        write(iunit,'(a)') ' Renormalized occupation matrices '

        do ios = 1, this%norbsets
          do ispin = 1,st%d%nspin
            write(iunit,'(a, i4, a, i4)') ' Orbital set ', ios, ' spin ', ispin
            do im = 1, this%orbsets(ios)%norbs
              write(iunit,'(1x)',advance='no')

              if (states_are_real(st)) then
                do imp = 1, this%orbsets(ios)%norbs-1
                  write(iunit,'(f14.8)',advance='no') this%dn_alt(im,imp,ispin,ios)
                end do
                write(iunit,'(f14.8)') this%dn_alt(im,this%orbsets(ios)%norbs,ispin,ios)
              else
                do imp = 1, this%orbsets(ios)%norbs-1
                  write(iunit,'(f14.8,f14.8)',advance='no') this%zn_alt(im,imp,ispin,ios)
                end do
                write(iunit,'(f14.8,f14.8)') this%zn_alt(im,this%orbsets(ios)%norbs,ispin,ios)
              end if
            end do
          end do !ispin
        end do !iatom
        call io_close(iunit)
      end if

    end if

    POP_SUB(lda_u_write_occupation_matrices)
  end subroutine lda_u_write_occupation_matrices

  !---------------------------------------------------------
  subroutine lda_u_write_effectiveU(dir, this, namespace)
    type(lda_u_t),     intent(in)    :: this
    character(len=*),  intent(in)    :: dir
    type(namespace_t), intent(in)    :: namespace

    integer :: iunit, ios

    PUSH_SUB(lda_u_write_effectiveU)

    if (mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
      iunit = io_open(trim(dir) // "/effectiveU", namespace, action='write')
      call lda_u_write_U(this, iunit)

      write(iunit, '(a,a,a,f7.3,a)') 'Hubbard U [', &
        trim(units_abbrev(units_out%energy)),']:'
      write(iunit,'(a,6x,14x,a)') ' Orbital',  'U'
      do ios = 1, this%norbsets
        if (.not. this%basisfromstates) then
          if (this%orbsets(ios)%ndim == 1) then
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
            else
              write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
            end if
          else
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
            else
              write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
            end if
          end if
        else
          write(iunit,'(i4,a10, 3x, f15.6)') ios, 'states', units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
        end if
      end do


      write(iunit, '(a,a,a,f7.3,a)') 'Hund J [', &
        trim(units_abbrev(units_out%energy)),']:'
      write(iunit,'(a,6x,14x,a)') ' Orbital',  'J'
      do ios = 1, this%norbsets
        if (.not. this%basisfromstates) then
          if (this%orbsets(ios)%ndim == 1) then
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
            else
              write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
            end if
          else
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
            else
              write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
            end if
          end if
        else
          write(iunit,'(i4,a10, f15.6)') ios, 'states', units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
        end if
      end do

      call io_close(iunit)
    end if

    POP_SUB(lda_u_write_effectiveU)
  end subroutine lda_u_write_effectiveU

  !---------------------------------------------------------
  subroutine lda_u_write_kanamoriU(dir, st, this, namespace)
    type(lda_u_t),       intent(in)    :: this
    type(states_elec_t), intent(in)    :: st
    character(len=*),    intent(in)    :: dir
    type(namespace_t),   intent(in)    :: namespace

    integer :: iunit, ios
    FLOAT, allocatable :: kanamori(:,:)

    PUSH_SUB(lda_u_write_kanamoriU)

    if (mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
      SAFE_ALLOCATE(kanamori(1:3,1:this%norbsets))

      call compute_ACBNO_U_kanamori(this, st, kanamori)

      iunit = io_open(trim(dir) // "/kanamoriU", namespace, action='write')

      write(iunit, '(a,a,a,f7.3,a)') 'Intraorbital U [', &
        trim(units_abbrev(units_out%energy)),']:'
      write(iunit,'(a,6x,14x,a)') ' Orbital',  'U'
      do ios = 1, this%norbsets
        if (.not. this%basisfromstates) then
          if (this%orbsets(ios)%ndim == 1) then
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, kanamori(1,ios))
            else
              write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, kanamori(1,ios))
            end if
          else
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, kanamori(1,ios))
            else
              write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, kanamori(1,ios))
            end if
          end if
        else
          write(iunit,'(i4,a10, 3x, f15.6)') ios, 'states', units_from_atomic(units_out%energy, kanamori(1,ios))
        end if
      end do


      write(iunit, '(a,a,a,f7.3,a)') 'Interorbital Up [', &
        trim(units_abbrev(units_out%energy)),']:'
      write(iunit,'(a,6x,14x,a)') ' Orbital',  'Up'
      do ios = 1, this%norbsets
        if (.not. this%basisfromstates) then
          if (this%orbsets(ios)%ndim == 1) then
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, kanamori(2,ios))
            else
              write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, kanamori(2,ios))
            end if
          else
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, kanamori(2,ios))
            else
              write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, kanamori(2,ios))
            end if
          end if
        else
          write(iunit,'(i4,a10, f15.6)') ios, 'states', units_from_atomic(units_out%energy, kanamori(2,ios))
        end if
      end do

      write(iunit, '(a,a,a,f7.3,a)') 'Hund J [', &
        trim(units_abbrev(units_out%energy)),']:'
      write(iunit,'(a,6x,14x,a)') ' Orbital',  'J'
      do ios = 1, this%norbsets
        if (.not. this%basisfromstates) then
          if (this%orbsets(ios)%ndim == 1) then
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, kanamori(3,ios))
            else
              write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                units_from_atomic(units_out%energy, kanamori(3,ios))
            end if
          else
            if (this%orbsets(ios)%nn /= 0) then
              write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, kanamori(3,ios))
            else
              write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
                units_from_atomic(units_out%energy, kanamori(3,ios))
            end if
          end if
        else
          write(iunit,'(i4,a10, f15.6)') ios, 'states', units_from_atomic(units_out%energy, kanamori(3,ios))
        end if
      end do


      call io_close(iunit)

      SAFE_DEALLOCATE_A(kanamori)
    end if


    POP_SUB(lda_u_write_kanamoriU)
  end subroutine lda_u_write_kanamoriU



  !---------------------------------------------------------
  subroutine lda_u_write_magnetization(dir, this, ions, mesh, st, namespace)
    type(lda_u_t),       intent(in)    :: this
    character(len=*),    intent(in)    :: dir
    type(ions_t),        intent(in)    :: ions
    class(mesh_t),       intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
    type(namespace_t),   intent(in)    :: namespace

    integer :: iunit, ia, ios, im
    FLOAT, allocatable :: mm(:,:)

    if (.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(lda_u_write_magnetization)

    call io_mkdir(dir, namespace)
    iunit = io_open(trim(dir)//"/magnetization.xsf", namespace, action='write', position='asis')

    if (this%nspins > 1) then
      SAFE_ALLOCATE(mm(1:mesh%box%dim, 1:ions%natoms))
      mm = M_ZERO
      !We compute the magnetization vector for each orbital set
      do ios = 1, this%norbsets
        ia = this%orbsets(ios)%iatom
        do im = 1, this%orbsets(ios)%norbs
          if (states_are_real(st)) then
            mm(3, ia) = mm(3, ia) + this%dn(im,im,1,ios) - this%dn(im,im,2,ios)
          else
            mm(3, ia) = mm(3, ia) + TOFLOAT(this%zn(im,im,1,ios) - this%zn(im,im,2,ios))
            !Spinors
            if (this%nspins /= this%spin_channels) then
              mm(1, ia) = mm(1, ia) + 2*TOFLOAT(this%zn(im,im,3,ios))
              mm(2, ia) = mm(2, ia) - 2*aimag(this%zn(im,im,3,ios))
            end if
          end if
        end do !im
      end do ! ios
      call write_xsf_geometry(iunit,  ions%space, ions%latt, ions%pos, ions%atom, mesh, forces = mm)
      SAFE_DEALLOCATE_A(mm)
    end if

    call io_close(iunit)

    POP_SUB(lda_u_write_magnetization)
  end subroutine lda_u_write_magnetization

  !---------------------------------------------------------
  subroutine lda_u_write_U(this, iunit, namespace)
    type(lda_u_t),               intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer :: ios

    PUSH_SUB(lda_u_write_U)

    write(message(1), '(a,a,a,f7.3,a)') 'Effective Hubbard U [', &
      trim(units_abbrev(units_out%energy)),']:'
    write(message(2),'(a,6x,14x,a)') ' Orbital',  'U'
    call messages_info(2, iunit=iunit, namespace=namespace)

    do ios = 1, this%norbsets
      if (.not. this%basisfromstates) then
        if (this%orbsets(ios)%ndim == 1) then
          if (this%orbsets(ios)%nn /= 0) then
            write(message(1),'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
              this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
              units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
          else
            write(message(1),'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
              l_notation(this%orbsets(ios)%ll), &
              units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
          end if
        else
          if (this%orbsets(ios)%nn /= 0) then
            write(message(1),'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
              this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
              int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
              units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
          else
            write(message(1),'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
              l_notation(this%orbsets(ios)%ll), &
              int(M_TWO*(this%orbsets(ios)%jj)), '/2', &
              units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
          end if
        end if
      else
        write(message(1),'(i4,a10, f15.6)') ios, 'states', units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
      end if
      call messages_info(1, iunit=iunit, namespace=namespace)
    end do

    POP_SUB(lda_u_write_U)
  end subroutine lda_u_write_U

  !---------------------------------------------------------
  subroutine lda_u_write_V(this, iunit, namespace)
    type(lda_u_t),               intent(in) :: this
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer :: ios, icopies, ios2

    if (.not. this%intersite) return

    PUSH_SUB(lda_u_write_V)

    write(message(1), '(a,a,a,f7.3,a)') 'Effective intersite V [', &
      trim(units_abbrev(units_out%energy)),']:'
    write(message(2),'(a,14x,a)') ' Orbital',  'V'
    call messages_info(2, iunit=iunit, namespace=namespace)

    do ios = 1, this%norbsets
      do icopies = 1, this%orbsets(ios)%nneighbors
        ios2 = this%orbsets(ios)%map_os(icopies)
        if (.not.this%basisfromstates) then
          if (this%orbsets(ios)%ndim == 1) then
            if (this%orbsets(ios)%nn /= 0) then
              write(message(1),'(i4,a10, 2x, i1, a1, i4, 1x, i1, a1, f7.3, f15.6)') ios, &
                trim(species_label(this%orbsets(ios)%spec)), this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), ios2, &
                this%orbsets(ios2)%nn, l_notation(this%orbsets(ios2)%ll), &
                units_from_atomic(units_out%length, this%orbsets(ios)%V_ij(icopies,3+1)), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%V_ij(icopies,0))
              call messages_info(1, iunit=iunit, namespace=namespace)
            else
              write(message(1),'(i4,a10, 3x, a1, i4, 1x, a1, f7.3, f15.6)') ios, &
                trim(species_label(this%orbsets(ios)%spec)), l_notation(this%orbsets(ios)%ll), ios2, &
                l_notation(this%orbsets(ios2)%ll), units_from_atomic(units_out%length, this%orbsets(ios)%V_ij(icopies,3+1)), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%V_ij(icopies,0))
              call messages_info(1, iunit=iunit, namespace=namespace)
            end if
          else
            if (this%orbsets(ios)%nn /= 0) then
              write(message(1),'(i4,a10, 2x, i1, a1, i1, a2, i4, f7.3, f15.6)') ios, &
                trim(species_label(this%orbsets(ios)%spec)), this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2',  ios2,      &
                units_from_atomic(units_out%length, this%orbsets(ios)%V_ij(icopies,3+1)), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%V_ij(icopies,0))
              call messages_info(1, iunit=iunit, namespace=namespace)
            else
              write(message(1),'(i4,a10, 3x, a1, i1, a2, i4, f7.3, f15.6)') ios, &
                trim(species_label(this%orbsets(ios)%spec)), l_notation(this%orbsets(ios)%ll), &
                int(M_TWO*(this%orbsets(ios)%jj)), '/2',  ios2,  &
                units_from_atomic(units_out%length, this%orbsets(ios)%V_ij(icopies,3+1)), &
                units_from_atomic(units_out%energy, this%orbsets(ios)%V_ij(icopies,0))
              call messages_info(1, iunit=iunit, namespace=namespace)
            end if ! this%orbsets(ios)%nn /= 0
          end if ! this%orbsets(ios)%ndim == 1
        else
          write(message(1),'(i4,a10, i4, f7.3, f15.6)') ios, 'states', ios2, &
            units_from_atomic(units_out%length, this%orbsets(ios)%V_ij(icopies,3+1)), &
            units_from_atomic(units_out%energy, this%orbsets(ios)%V_ij(icopies,0))
          call messages_info(1, iunit=iunit, namespace=namespace)
        end if

      end do ! icopies
    end do ! ios

    POP_SUB(lda_u_write_V)
  end subroutine lda_u_write_V


  ! ---------------------------------------------------------
  subroutine lda_u_dump(restart, namespace, this, st, mesh, ierr)
    type(restart_t),      intent(in)  :: restart
    type(namespace_t),    intent(in)  :: namespace
    type(lda_u_t),        intent(in)  :: this
    type(states_elec_t),  intent(in)  :: st
    class(mesh_t),        intent(in)  :: mesh
    integer,              intent(out) :: ierr

    integer :: err, occsize, ios, ncount
    FLOAT, allocatable :: Ueff(:), docc(:), Veff(:)
    CMPLX, allocatable :: zocc(:)

    PUSH_SUB(lda_u_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(lda_u_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing DFT+U restart."
      call messages_info(1, namespace=namespace)
    end if

    occsize = this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets
    if (this%level == DFT_U_ACBN0) then
      occsize = occsize*2
      if (this%intersite) then
        occsize = occsize + 2*this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets*this%maxneighbors
      end if
    end if


    if (states_are_real(st)) then
      SAFE_ALLOCATE(docc(1:occsize))
      docc = M_ZERO
      call dlda_u_get_occupations(this, docc)
      call drestart_write_binary(restart, "lda_u_occ", occsize, docc, err)
      if (err /= 0) ierr = ierr + 1
      SAFE_DEALLOCATE_A(docc)
    else
      SAFE_ALLOCATE(zocc(1:occsize))
      zocc = M_ZERO
      call zlda_u_get_occupations(this, zocc)
      call zrestart_write_binary(restart, "lda_u_occ", occsize, zocc, err)
      if (err /= 0) ierr = ierr + 1
      SAFE_DEALLOCATE_A(zocc)
    end if


    if (this%level == DFT_U_ACBN0) then
      SAFE_ALLOCATE(Ueff(1:this%norbsets))
      Ueff = M_ZERO
      call lda_u_get_effectiveU(this, Ueff(:))
      call drestart_write_binary(restart, "lda_u_Ueff", this%norbsets, Ueff, err)
      SAFE_DEALLOCATE_A(Ueff)
      if (err /= 0) ierr = ierr + 1

      if (this%intersite .and. this%maxneighbors > 0) then
        ncount = 0
        do ios = 1, this%norbsets
          ncount = ncount + this%orbsets(ios)%nneighbors
        end do
        SAFE_ALLOCATE(Veff(1:ncount))
        call lda_u_get_effectiveV(this, Veff(:))
        call drestart_write_binary(restart, "lda_u_Veff", ncount, Veff, err)
        SAFE_DEALLOCATE_A(Veff)
        if (err /= 0) ierr = ierr + 1
      end if
    end if

    call lda_u_dump_coulomb_integrals(this, namespace, restart, st, mesh, ierr)

    if (debug%info) then
      message(1) = "Debug: Writing DFT+U restart done."
      call messages_info(1, namespace=namespace)
    end if

    POP_SUB(lda_u_dump)
  end subroutine lda_u_dump


  ! ---------------------------------------------------------
  subroutine lda_u_load(restart, this, st, dftu_energy, ierr, occ_only, u_only)
    type(restart_t),      intent(in)    :: restart
    type(lda_u_t),        intent(inout) :: this
    type(states_elec_t),  intent(in)    :: st
    FLOAT,                intent(out)   :: dftu_energy
    integer,              intent(out)   :: ierr
    logical, optional,    intent(in)    :: occ_only
    logical, optional,    intent(in)    :: u_only

    integer :: err, occsize, ncount, ios
    FLOAT, allocatable :: Ueff(:), docc(:), Veff(:)
    CMPLX, allocatable :: zocc(:)

    PUSH_SUB(lda_u_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(lda_u_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading DFT+U restart."
      call messages_info(1)
    end if

    !We have to read the effective U first, as we call lda_u_uptade_potential latter
    if (this%level == DFT_U_ACBN0 .and. .not. optional_default(occ_only, .false.)) then
      SAFE_ALLOCATE(Ueff(1:this%norbsets))
      call drestart_read_binary(restart, "lda_u_Ueff", this%norbsets, Ueff, err)
      if (err /= 0) then
        ierr = ierr + 1
        Ueff = M_ZERO
      end if
      call lda_u_set_effectiveU(this, Ueff)
      SAFE_DEALLOCATE_A(Ueff)

      if (this%intersite .and. this%maxneighbors > 0) then
        ncount = 0
        do ios = 1, this%norbsets
          ncount = ncount + this%orbsets(ios)%nneighbors
        end do
        SAFE_ALLOCATE(Veff(1:ncount))
        call drestart_read_binary(restart, "lda_u_Veff", ncount, Veff, err)
        if (err /= 0) then
          ierr = ierr + 1
          Veff = M_ZERO
        end if
        call lda_u_set_effectiveV(this, Veff)
        SAFE_DEALLOCATE_A(Veff)
      end if
    end if


    if (.not. optional_default(u_only, .false.)) then
      occsize = this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets
      if (this%level == DFT_U_ACBN0) then
        occsize = occsize*2
        if (this%intersite) then
          occsize = occsize + 2*this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets*this%maxneighbors
        end if
      end if


      if (states_are_real(st)) then
        SAFE_ALLOCATE(docc(1:occsize))
        call drestart_read_binary(restart, "lda_u_occ", occsize, docc, err)
        if (err /= 0) then
          ierr = ierr + 1
          docc = M_ZERO
        end if
        call dlda_u_set_occupations(this, docc)
        SAFE_DEALLOCATE_A(docc)
      else
        SAFE_ALLOCATE(zocc(1:occsize))
        call zrestart_read_binary(restart, "lda_u_occ", occsize, zocc, err)
        if (err /= 0) then
          ierr = ierr + 1
          zocc = M_z0
        end if
        call zlda_u_set_occupations(this, zocc)
        SAFE_DEALLOCATE_A(zocc)
      end if
    end if

    if (states_are_real(st)) then
      call dcompute_dftu_energy(this, dftu_energy, st)
      call dlda_u_update_potential(this, st)
    else
      call zcompute_dftu_energy(this, dftu_energy, st)
      call zlda_u_update_potential(this, st)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading DFT+U restart done."
      call messages_info(1)
    end if

    POP_SUB(lda_u_load)
  end subroutine lda_u_load

  ! ---------------------------------------------------------
  subroutine lda_u_dump_coulomb_integrals(this, namespace, restart, st, mesh, ierr)
    type(lda_u_t),        intent(in)    :: this
    type(namespace_t),    intent(in)    :: namespace
    type(restart_t),      intent(in)    :: restart
    type(states_elec_t),  intent(in)    :: st
    class(mesh_t),        intent(in)    :: mesh
    integer,              intent(out)   :: ierr

    integer              :: coulomb_int_file, idim1, idim2, err
    integer              :: ios, im, imp, impp, imppp
    character(len=256)   :: lines(3)
    logical              :: complex_coulomb_integrals

    PUSH_SUB(lda_u_dump_coulomb_integrals)

    ierr = 0

    if (restart_skip(restart) .or. this%level == DFT_U_EMPIRICAL) then
      ierr = -1
      POP_SUB(lda_u_dump_coulomb_integrals)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing Coulomb integrals restart."
      call messages_info(1)
    end if

    complex_coulomb_integrals = .false.
    do ios = 1, this%norbsets
      if (this%orbsets(ios)%ndim  > 1) complex_coulomb_integrals = .true.
    end do

    coulomb_int_file  = restart_open(restart, 'coulomb_integrals')
    ! sanity checks. Example file 'coulomb_int_file':
    ! norb=                        2
    ! dim=                         1
    ! checksum=                    xxxxxxxxxxx
    write(lines(1), '(a20,i21)') "norb=", this%norbsets
    write(lines(2), '(a20,i21)') "dim=", st%d%dim
    write(lines(3), '(a20,i21)') "checksum=", mesh%idx%checksum
    call restart_write(restart, coulomb_int_file, lines, 3, err)
    if (err /= 0) ierr = ierr - 2

    do ios = 1, this%norbsets
      do im = 1, this%orbsets(ios)%norbs
        do imp = 1, this%orbsets(ios)%norbs
          do impp = 1, this%orbsets(ios)%norbs
            do imppp = 1, this%orbsets(ios)%norbs
              if(.not. complex_coulomb_integrals) then
                write(lines(1), '(i4,i4,i4,i4,i4,e20.12)') ios, im, imp, impp, imppp, this%coulomb(im, imp, impp, imppp, ios)
              else
                do idim1 = 1, st%d%dim
                  do idim2 = 1, st%d%dim
                    write(lines(1), '(i4,i4,i4,i4,i4,2e20.12)') ios, im, imp, impp, imppp, idim1, idim2, &
                      TOFLOAT(this%zcoulomb(im, imp, impp, imppp, idim1, idim2, ios)), &
                      aimag(this%zcoulomb(im, imp, impp, imppp, idim1, idim2, ios))
                  end do
                end do
              end if
              call restart_write(restart, coulomb_int_file, lines, 1, err)
              if (err /=0) then
                ierr = ierr - 2**2
                exit
              end if
            end do
          end do
        end do
      end do
    end do
    call restart_close(restart, coulomb_int_file)

    if (debug%info) then
      message(1) = "Debug: Writing Coulomb integrals restart done."
      call messages_info(1)
    end if

    POP_SUB(lda_u_dump_coulomb_integrals)
  end subroutine lda_u_dump_coulomb_integrals

end module lda_u_io_oct_m
