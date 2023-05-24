!! Copyright (C) 2021 S. Ohlmann, I. Albar, A. Obzhirov, A. Geondzhian, M. Lawan
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

module propagator_data_classical_particles_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_beeman_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m

  implicit none

  private
  public ::                                      &
    propagator_data_t

  type :: propagator_data_t
    !> The following variables are work arrays used by the different propagators:
    FLOAT, allocatable :: acc(:,:)        !< Acceleration of the particles
    FLOAT, allocatable :: prev_acc(:,:,:) !< A storage of the prior times.
    FLOAT, allocatable :: save_pos(:,:)   !< A storage for the SCF loops
    FLOAT, allocatable :: save_vel(:,:)   !< A storage for the SCF loops
    FLOAT, allocatable :: prev_tot_force(:,:) !< Used for the SCF convergence criterium
    FLOAT, allocatable :: prev_pos(:,:,:) !< Used for extrapolation
    FLOAT, allocatable :: prev_vel(:,:,:) !< Used for extrapolation
    FLOAT, allocatable :: hamiltonian_elements(:,:)
    logical :: initialized = .false.
  contains
    procedure :: restart_write => propagator_data_restart_write
    procedure :: restart_read => propagator_data_restart_read
    procedure :: initialize => propagator_data_initialize
    procedure :: end => propagator_data_end
    procedure :: propagator_data_copy
    generic   :: assignment(=) => propagator_data_copy
  end type propagator_data_t

contains
  ! ---------------------------------------------------------
  subroutine propagator_data_restart_write(this, namespace, prop)
    class(propagator_data_t), intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    class(propagator_t),      intent(in)    :: prop

    integer                                 :: restart_file_unit

    PUSH_SUB(propagator_data_restart_write)
    call io_mkdir('restart/'//TD_DIR, namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_classical_particles_propagation', &
      namespace, action='write')

    select type(prop)
    type is (propagator_verlet_t)
      write(restart_file_unit,*) this%acc(:,:)
      write(restart_file_unit,*) this%prev_acc(:,:,:)
    type is (propagator_beeman_t)
      write(restart_file_unit,*) this%acc(:,:)
      write(restart_file_unit,*) this%prev_acc(:,:,:)
      if (prop%predictor_corrector) then
        write(restart_file_unit,*) this%prev_tot_force(:,:)
        write(restart_file_unit,*) this%save_vel(:,:)
        write(restart_file_unit,*) this%save_pos(:,:)
      end if
    type is (propagator_exp_mid_t)
      write(restart_file_unit,*) this%prev_vel(:,:,:)
      write(restart_file_unit,*) this%prev_pos(:,:,:)
      write(restart_file_unit,*) this%save_vel(:,:)
      write(restart_file_unit,*) this%save_pos(:,:)
    end select

    call io_close(restart_file_unit)

    POP_SUB(propagator_data_restart_write)
  end subroutine propagator_data_restart_write

  ! ---------------------------------------------------------
  logical function propagator_data_restart_read(this, namespace, prop)
    class(propagator_data_t), intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    class(propagator_t),      intent(in)    :: prop

    integer                                 :: restart_file_unit

    PUSH_SUB(propagator_data_restart_read)

    call io_mkdir('restart/'//TD_DIR, namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_classical_particles_propagation', namespace, &
      action='read', die=.false.)
    if (restart_file_unit > 0) then
      select type(prop)
      type is (propagator_verlet_t)
        read(restart_file_unit,*) this%acc(:,:)
        read(restart_file_unit,*) this%prev_acc(:,:,:)

      type is (propagator_beeman_t)
        read(restart_file_unit,*) this%acc(:,:)
        read(restart_file_unit,*) this%prev_acc(:,:,:)
        if (prop%predictor_corrector) then
          read(restart_file_unit,*) this%prev_tot_force(:,:)
          read(restart_file_unit,*) this%save_vel(:,:)
          read(restart_file_unit,*) this%save_pos(:,:)
        end if
      type is (propagator_exp_mid_t)
        read(restart_file_unit,*) this%prev_vel(:,:,:)
        read(restart_file_unit,*) this%prev_pos(:,:,:)
        read(restart_file_unit,*) this%save_vel(:,:)
        read(restart_file_unit,*) this%save_pos(:,:)
      end select

      call io_close(restart_file_unit)
      propagator_data_restart_read = .true.
    else
      ! error opening file
      propagator_data_restart_read = .false.
    end if

    POP_SUB(propagator_data_restart_read)
  end function propagator_data_restart_read

  ! ---------------------------------------------------------
  subroutine propagator_data_initialize(this, prop, dim, np)
    class(propagator_data_t), intent(inout) :: this
    class(propagator_t),      intent(in)    :: prop
    integer,                  intent(in)    :: dim
    integer,                  intent(in)    :: np

    PUSH_SUB(propagator_data_initialize)

    if (.not. this%initialized) then
      select type(prop)
      type is (propagator_verlet_t)
        SAFE_ALLOCATE(this%acc(1:dim, 1:np))
        SAFE_ALLOCATE(this%prev_acc(1:dim, 1:np, 1))
      type is (propagator_beeman_t)
        if (prop%predictor_corrector) then
          SAFE_ALLOCATE(this%save_pos(1:dim, 1:np))
          SAFE_ALLOCATE(this%save_vel(1:dim, 1:np))
          SAFE_ALLOCATE(this%prev_tot_force(1:dim, 1:np))
        end if
        SAFE_ALLOCATE(this%acc(1:dim, 1:np))
        SAFE_ALLOCATE(this%prev_acc(1:dim, 1:np, 1:2))
      type is (propagator_exp_mid_t)
        SAFE_ALLOCATE(this%save_pos(1:dim, 1:np))
        SAFE_ALLOCATE(this%save_vel(1:dim, 1:np))
        SAFE_ALLOCATE(this%hamiltonian_elements(1:dim, 1:np))
        SAFE_ALLOCATE(this%prev_pos(1:dim, 1:np, 1))
        SAFE_ALLOCATE(this%prev_vel(1:dim, 1:np, 1))
      end select
      this%initialized = .true.
    end if

    POP_SUB(propagator_data_initialize)
  end subroutine propagator_data_initialize

  ! ---------------------------------------------------------
  subroutine propagator_data_end(this)
    class(propagator_data_t),     intent(inout) :: this

    PUSH_SUB(propagator_data_end)

    SAFE_DEALLOCATE_A(this%acc)
    SAFE_DEALLOCATE_A(this%prev_acc)
    SAFE_DEALLOCATE_A(this%prev_tot_force)
    SAFE_DEALLOCATE_A(this%save_pos)
    SAFE_DEALLOCATE_A(this%save_vel)
    SAFE_DEALLOCATE_A(this%hamiltonian_elements)
    SAFE_DEALLOCATE_A(this%prev_pos)
    SAFE_DEALLOCATE_A(this%prev_vel)

    POP_SUB(propagator_data_end)
  end subroutine propagator_data_end

  ! ---------------------------------------------------------
  subroutine propagator_data_copy(this, prop_data_in)
    class(propagator_data_t), intent(inout) :: this
    class(propagator_data_t), intent(in)    :: prop_data_in

    PUSH_SUB(propagator_data_copy)

    SAFE_ALLOCATE_SOURCE_A(this%acc,                  prop_data_in%acc)
    SAFE_ALLOCATE_SOURCE_A(this%prev_acc,             prop_data_in%prev_acc)
    SAFE_ALLOCATE_SOURCE_A(this%save_pos,             prop_data_in%save_pos)
    SAFE_ALLOCATE_SOURCE_A(this%save_vel,             prop_data_in%save_vel)
    SAFE_ALLOCATE_SOURCE_A(this%prev_tot_force,       prop_data_in%prev_tot_force)
    SAFE_ALLOCATE_SOURCE_A(this%prev_pos,             prop_data_in%prev_pos)
    SAFE_ALLOCATE_SOURCE_A(this%prev_vel,             prop_data_in%prev_vel)
    SAFE_ALLOCATE_SOURCE_A(this%hamiltonian_elements, prop_data_in%hamiltonian_elements)

    POP_SUB(propagator_data_copy)
  end subroutine propagator_data_copy
end module propagator_data_classical_particles_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
