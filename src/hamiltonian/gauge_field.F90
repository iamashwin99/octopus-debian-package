!! Copyright (C) 2008 X. Andrade
!! Copyright (C) 2020 N. Tancogne-Dejean
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

module gauge_field_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use debug_oct_m
  use ghost_interaction_oct_m
  use global_oct_m
  use grid_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use iso_c_binding
  use kpoints_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use propagator_verlet_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use symmetries_oct_m
  use symm_op_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private

  public ::                               &
    gauge_field_t,                        &
    gauge_field_init,                     &
    gauge_field_init_vec_pot,             &
    gauge_field_is_propagated,            &
    gauge_field_is_used,                  &
    gauge_field_set_vec_pot,              &
    gauge_field_set_vec_pot_vel,          &
    gauge_field_get_vec_pot,              &
    gauge_field_get_vec_pot_vel,          &
    gauge_field_get_vec_pot_acc,          &
    gauge_field_get_energy,               &
    gauge_field_dump,                     &
    gauge_field_load,                     &
    gauge_field_end,                      &
    gauge_field_get_force,                &
    gauge_field_do_algorithmic_operation,                    &
    gauge_field_output_write,             &
    gauge_field_check_symmetries

  type, extends(interaction_partner_t) :: gauge_field_t
    private
    FLOAT, allocatable :: vecpot(:)
    FLOAT, allocatable :: vecpot_vel(:)
    FLOAT, allocatable :: vecpot_acc(:)
    FLOAT, allocatable :: vecpot_kick(:)
    FLOAT, allocatable :: force(:)
    FLOAT   :: wp2
    logical :: with_gauge_field = .false.
    integer :: dynamics
    FLOAT   :: kicktime

    FLOAT   :: volume

  contains
    procedure :: update_exposed_quantities => gauge_field_update_exposed_quantities
    procedure :: update_exposed_quantity => gauge_field_update_exposed_quantity
    procedure :: init_interaction_as_partner => gauge_field_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => gauge_field_copy_quantities_to_interaction
    final :: gauge_field_finalize
  end type

  interface gauge_field_t
    module procedure gauge_field_init
  end interface gauge_field_t


contains

  ! ---------------------------------------------------------
  function gauge_field_init(namespace, volume) result(this)
    class(gauge_field_t), pointer :: this
    type(namespace_t), intent(in) :: namespace
    FLOAT,             intent(in) :: volume

    integer :: ii
    type(block_t) :: blk

    PUSH_SUB(gauge_field_init)

    SAFE_ALLOCATE(this)

    this%namespace = namespace_t("GaugeField", parent=namespace)
    call space_init(this%space, this%namespace)

    ! Initialize clock without a time-step, as the gauge-field will not be propagated directly
    ! for the moment
    this%clock = clock_t()

    this%with_gauge_field = .false.

    this%volume = volume

    SAFE_ALLOCATE(this%vecpot(1:this%space%dim))
    SAFE_ALLOCATE(this%vecpot_vel(1:this%space%dim))
    SAFE_ALLOCATE(this%vecpot_acc(1:this%space%dim))
    SAFE_ALLOCATE(this%vecpot_kick(1:this%space%dim))
    SAFE_ALLOCATE(this%force(1:this%space%dim))
    this%vecpot = M_ZERO
    this%vecpot_vel = M_ZERO
    this%vecpot_acc = M_ZERO
    this%vecpot_kick = M_ZERO
    this%force = M_ZERO

    !%Variable GaugeFieldDynamics
    !%Type integer
    !%Default polarization
    !%Section Hamiltonian
    !%Description
    !% This variable select the dynamics of the gauge field used to
    !% apply a finite electric field to periodic systems in
    !% time-dependent runs.
    !%Option none 0
    !% The gauge field does not have dynamics. The induced polarization field is zero.
    !%Option polarization 1
    !% The gauge field follows the dynamic described in
    !% Bertsch et al, Phys. Rev. B 62 7998 (2000).
    !%End

    call parse_variable(namespace, 'GaugeFieldDynamics', OPTION__GAUGEFIELDDYNAMICS__POLARIZATION, this%dynamics)

    !%Variable GaugeFieldPropagate
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% Propagate the gauge field with initial condition set by GaugeVectorField or zero if not specified
    !%End

    call parse_variable(namespace, 'GaugeFieldPropagate', .false., this%with_gauge_field)

    !%Variable GaugeVectorField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% The gauge vector field is used to include a uniform (but time-dependent)
    !% external electric field in a time-dependent run for
    !% a periodic system. An optional second row specifies the initial
    !% value for the time derivative of the gauge field (which is set
    !% to zero by default). By default this field is not included.
    !% If <tt>KPointsUseSymmetries = yes</tt>, then <tt>SymmetryBreakDir</tt>
    !% must be set in the same direction.
    !% This is used with utility <tt>oct-dielectric_function</tt>
    !% according to GF Bertsch, J-I Iwata, A Rubio, and K Yabana,
    !% <i>Phys. Rev. B</i> <b>62</b>, 7998-8002 (2000).
    !%End
    ! Read the initial gauge vector field

    if (parse_block(namespace, 'GaugeVectorField', blk) == 0) then

      this%with_gauge_field = .true.

      do ii = 1, this%space%dim
        call parse_block_float(blk, 0, ii - 1, this%vecpot_kick(ii))
      end do

      call parse_block_end(blk)
      if (.not. this%space%is_periodic()) then
        message(1) = "GaugeVectorField is intended for periodic systems."
        call messages_warning(1, namespace=namespace)
      end if

    end if

    !%Variable GaugeFieldDelay
    !%Type float
    !%Default 0.
    !%Section Hamiltonian
    !%Description
    !% The application of the gauge field acts as a probe of the system. For dynamical
    !% systems one can apply this probe with a delay relative to the start of the simulation.
    !%End

    call parse_variable(namespace, 'GaugeFieldDelay', M_ZERO, this%kicktime)

    if (abs(this%kicktime) <= M_EPSILON) then
      this%vecpot(1:this%space%dim) = this%vecpot_kick(1:this%space%dim)
    end if

    POP_SUB(gauge_field_init)
  end function gauge_field_init

  ! ---------------------------------------------------------
  subroutine gauge_field_check_symmetries(this, kpoints)
    type(gauge_field_t),  intent(in) :: this
    type(kpoints_t),      intent(in) :: kpoints

    integer :: iop

    PUSH_SUB(gauge_field_check_symmetries)

    if (kpoints%use_symmetries) then
      do iop = 1, symmetries_number(kpoints%symm)
        if (iop == symmetries_identity_index(kpoints%symm)) cycle
        if (.not. symm_op_invariant_cart(kpoints%symm%ops(iop), this%vecpot_kick, CNST(1e-5))) then
          message(1) = "The GaugeVectorField breaks (at least) one of the symmetries used to reduce the k-points."
          message(2) = "Set SymmetryBreakDir equal to GaugeVectorField."
          call messages_fatal(2, namespace=this%namespace)
        end if
      end do
    end if

    POP_SUB(gauge_field_check_symmetries)
  end subroutine gauge_field_check_symmetries


  ! ---------------------------------------------------------
  subroutine gauge_field_end(this)
    type(gauge_field_t),     intent(inout) :: this

    PUSH_SUB(gauge_field_end)

    this%with_gauge_field = .false.
    SAFE_DEALLOCATE_A(this%vecpot)
    SAFE_DEALLOCATE_A(this%vecpot_vel)
    SAFE_DEALLOCATE_A(this%vecpot_acc)
    SAFE_DEALLOCATE_A(this%vecpot_kick)
    SAFE_DEALLOCATE_A(this%force)

    POP_SUB(gauge_field_end)
  end subroutine gauge_field_end


  ! ---------------------------------------------------------
  logical pure function gauge_field_is_propagated(this) result(is_propagated)
    type(gauge_field_t),  intent(in) :: this

    is_propagated = this%with_gauge_field
  end function gauge_field_is_propagated

  ! ---------------------------------------------------------
  logical pure function gauge_field_is_used(this) result(is_used)
    type(gauge_field_t),  intent(in) :: this

    is_used = this%with_gauge_field .or. (norm2(this%vecpot_kick(1:this%space%dim)) > M_EPSILON)
  end function gauge_field_is_used



  ! ---------------------------------------------------------
  subroutine gauge_field_set_vec_pot(this, vec_pot)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot(:) !< (this%space%dim)

    PUSH_SUB(gauge_field_set_vec_pot)
    this%vecpot(1:this%space%dim) = vec_pot(1:this%space%dim)

    POP_SUB(gauge_field_set_vec_pot)
  end subroutine gauge_field_set_vec_pot


  ! ---------------------------------------------------------
  subroutine gauge_field_set_vec_pot_vel(this, vec_pot_vel)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot_vel(:) !< (this%space%dim)

    PUSH_SUB(gauge_field_set_vec_pot_vel)
    this%vecpot_vel(1:this%space%dim) = vec_pot_vel(1:this%space%dim)

    POP_SUB(gauge_field_set_vec_pot_vel)
  end subroutine gauge_field_set_vec_pot_vel


  ! ---------------------------------------------------------
  subroutine gauge_field_get_vec_pot(this, vec_pot)
    type(gauge_field_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot(:) !< (this%space%dim)

    PUSH_SUB(gauge_field_get_vec_pot)
    vec_pot(1:this%space%dim) = this%vecpot(1:this%space%dim)

    POP_SUB(gauge_field_get_vec_pot)
  end subroutine gauge_field_get_vec_pot


  ! ---------------------------------------------------------
  subroutine gauge_field_get_vec_pot_vel(this, vec_pot_vel)
    type(gauge_field_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot_vel(:) !< (this%space%dim)

    PUSH_SUB(gauge_field_get_vec_pot_vel)
    vec_pot_vel(1:this%space%dim) = this%vecpot_vel(1:this%space%dim)

    POP_SUB(gauge_field_get_vec_pot_vel)
  end subroutine gauge_field_get_vec_pot_vel


  ! ---------------------------------------------------------
  subroutine gauge_field_get_vec_pot_acc(this, vec_pot_acc)
    type(gauge_field_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot_acc(:) !< (this%space%dim)

    PUSH_SUB(gauge_field_get_vec_pot_acc)
    vec_pot_acc(1:this%space%dim) = this%vecpot_acc(1:this%space%dim)

    POP_SUB(gauge_field_get_vec_pot_acc)
  end subroutine gauge_field_get_vec_pot_acc

  ! ---------------------------------------------------------
  subroutine gauge_field_propagate(this, dt, time)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: dt
    FLOAT,                intent(in)    :: time

    logical, save :: warning_shown = .false.
    integer :: idim

    PUSH_SUB(gauge_field_propagate)

    this%vecpot_acc(1:this%space%dim) = this%force(1:this%space%dim)

    ! apply kick, in case kicktime=0 the kick has already been applied
    if (this%kicktime > M_ZERO .and. time-dt <= this%kicktime .and. time > this%kicktime) then
      this%vecpot(1:this%space%dim) = this%vecpot(1:this%space%dim) +  this%vecpot_kick(1:this%space%dim)
      call messages_write('     ----------------  Applying gauge kick  ----------------')
      call messages_info(namespace=this%namespace)
    end if

    this%vecpot(1:this%space%dim) = this%vecpot(1:this%space%dim) + dt * this%vecpot_vel(1:this%space%dim) + &
      M_HALF * dt**2 * this%force(1:this%space%dim)

    !In the case of a kick, the induced field could not be higher than the initial kick
    do idim = 1, this%space%dim
      if (.not. warning_shown .and. this%vecpot_kick(idim) /= M_ZERO .and.  &
        abs(this%vecpot(idim))> abs(this%vecpot_kick(idim))*1.01 .and. .not. this%kicktime > M_ZERO) then

        warning_shown = .true.

        write(message(1),'(a)') 'It seems that the gauge-field might be diverging. You should probably check'
        write(message(2),'(a)') 'the simulation parameters, in particular the number of k-points.'
        call messages_warning(2, namespace=this%namespace)
      end if
    end do
    POP_SUB(gauge_field_propagate)
  end subroutine gauge_field_propagate

  ! ---------------------------------------------------------
  subroutine gauge_field_init_vec_pot(this, qtot)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: qtot

    PUSH_SUB(gauge_field_init_vec_pot)

    this%wp2 = M_FOUR*M_PI*qtot/this%volume

    write (message(1), '(a,f12.6,a)') "Info: Electron-gas plasmon frequency", &
      units_from_atomic(units_out%energy, sqrt(this%wp2)), " ["//trim(units_abbrev(units_out%energy))//"]"
    call messages_info(1, namespace=this%namespace)

    POP_SUB(gauge_field_init_vec_pot)
  end subroutine gauge_field_init_vec_pot

  ! ---------------------------------------------------------
  FLOAT function gauge_field_get_energy(this) result(energy)
    type(gauge_field_t),  intent(in)    :: this

    PUSH_SUB(gauge_field_get_energy)

    if(allocated(this%vecpot_vel)) then
      energy = this%volume / (CNST(8.0) * M_PI * P_c**2) * sum(this%vecpot_vel(1:this%space%dim)**2)
    else
      energy = M_ZERO
    end if

    POP_SUB(gauge_field_get_energy)
  end function gauge_field_get_energy


  ! ---------------------------------------------------------
  subroutine gauge_field_dump(restart, gfield, ierr)
    type(restart_t),      intent(in)  :: restart
    type(gauge_field_t),  intent(in)  :: gfield
    integer,              intent(out) :: ierr

    integer :: err
    FLOAT, allocatable :: vecpot(:,:)

    PUSH_SUB(gauge_field_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(gauge_field_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing gauge field restart."
      call messages_info(1, namespace=gfield%namespace)
    end if

    SAFE_ALLOCATE(vecpot(1:gfield%space%dim, 1:2))
    vecpot = M_ZERO
    call gauge_field_get_vec_pot(gfield, vecpot(:, 1))
    call gauge_field_get_vec_pot_vel(gfield, vecpot(:, 2))

    call drestart_write_binary(restart, "gauge_field", 2*gfield%space%dim, vecpot, err)
    SAFE_DEALLOCATE_A(vecpot)
    if (err /= 0) ierr = ierr + 1

    if (debug%info) then
      message(1) = "Debug: Writing gauge field restart done."
      call messages_info(1, namespace=gfield%namespace)
    end if

    POP_SUB(gauge_field_dump)
  end subroutine gauge_field_dump


  ! ---------------------------------------------------------
  subroutine gauge_field_load(restart, gfield, ierr)
    type(restart_t),      intent(in)    :: restart
    type(gauge_field_t),  intent(inout) :: gfield
    integer,              intent(out)   :: ierr

    integer :: err
    FLOAT, allocatable :: vecpot(:,:)

    PUSH_SUB(gauge_field_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(gauge_field_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading gauge field restart."
      call messages_info(1, namespace=gfield%namespace)
    end if

    SAFE_ALLOCATE(vecpot(1:gfield%space%dim, 1:2))
    call drestart_read_binary(restart, "gauge_field", 2*gfield%space%dim, vecpot, err)
    if (err /= 0) ierr = ierr + 1

    call gauge_field_set_vec_pot(gfield, vecpot(:,1))
    call gauge_field_set_vec_pot_vel(gfield, vecpot(:,2))
    SAFE_DEALLOCATE_A(vecpot)

    if (debug%info) then
      message(1) = "Debug: Reading gauge field restart done."
      call messages_info(1, namespace=gfield%namespace)
    end if

    POP_SUB(gauge_field_load)
  end subroutine gauge_field_load

  ! ---------------------------------------------------------

  subroutine gauge_field_get_force(this, gr, spin_channels, current, lrc_alpha)
    type(gauge_field_t),     intent(inout) :: this
    type(grid_t),            intent(in)    :: gr
    integer,                 intent(in)    :: spin_channels
    FLOAT,                   intent(in)    :: current(:,:,:)
    FLOAT,optional,          intent(in)    :: lrc_alpha

    integer :: idir, ispin
    FLOAT :: lrc_alpha_

    PUSH_SUB(gauge_field_get_force)

    lrc_alpha_ = optional_default(lrc_alpha, M_ZERO)

    select case (this%dynamics)
    case (OPTION__GAUGEFIELDDYNAMICS__NONE)
      this%force(1:this%space%dim) = M_ZERO

    case (OPTION__GAUGEFIELDDYNAMICS__POLARIZATION)
      do idir = 1, this%space%periodic_dim
        this%force(idir) = M_ZERO
        do ispin = 1, spin_channels
          if(lrc_alpha_ > M_EPSILON) then
            this%force(idir) = this%force(idir) + &
              lrc_alpha*P_c/this%volume*dmf_integrate(gr, current(:, idir, ispin))
          else
            this%force(idir) = this%force(idir) - &
              M_FOUR*M_PI*P_c/this%volume*dmf_integrate(gr, current(:, idir, ispin))
          endif
        end do
      end do

    case default
      ASSERT(.false.)
    end select

    POP_SUB(gauge_field_get_force)
  end subroutine gauge_field_get_force

  ! ---------------------------------------------------------

  subroutine gauge_field_do_algorithmic_operation(this, operation, dt, time)
    class(gauge_field_t),          intent(inout) :: this
    type(algorithmic_operation_t), intent(in)    :: operation
    FLOAT,                         intent(in)    :: dt
    FLOAT,                         intent(in)    :: time

    PUSH_SUB(gauge_field_do_algorithmic_operation)

    select case (operation%id)
    case (VERLET_START)
      !Does nothing at the moment
    case (VERLET_FINISH)
      !Does nothing at the moment
    case (VERLET_UPDATE_POS)
      !This is inside the gauge_field_propagate routine at the moment
    case (VERLET_COMPUTE_ACC)
      call gauge_field_propagate(this, dt, time)

    case (VERLET_COMPUTE_VEL)
      this%vecpot_vel(1:this%space%dim) = this%vecpot_vel(1:this%space%dim) + &
        M_HALF * dt * (this%vecpot_acc(1:this%space%dim) + this%force(1:this%space%dim))
    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(gauge_field_do_algorithmic_operation)
  end subroutine gauge_field_do_algorithmic_operation


  ! ---------------------------------------------------------
  subroutine gauge_field_output_write(this, out_gauge, iter)
    type(gauge_field_t), intent(in)    :: this
    type(c_ptr),         intent(inout) :: out_gauge
    integer,             intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux
    FLOAT :: temp(this%space%dim)

    if (.not. mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(gauge_field_output_write)

    if (iter == 0) then
      call write_iter_clear(out_gauge)
      call write_iter_string(out_gauge,'################################################################################')
      call write_iter_nl(out_gauge)
      call write_iter_string(out_gauge,'# HEADER')
      call write_iter_nl(out_gauge)

      ! first line: column names
      call write_iter_header_start(out_gauge)

      do idir = 1, this%space%dim
        write(aux, '(a2,i1,a1)') 'A(', idir, ')'
        call write_iter_header(out_gauge, aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a6,i1,a1)') 'dA/dt(', idir, ')'
        call write_iter_header(out_gauge, aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a10,i1,a1)') 'd^2A/dt^2(', idir, ')'
        call write_iter_header(out_gauge, aux)
      end do
      call write_iter_nl(out_gauge)

      ! second line: units
      !call write_iter_string(out_gauge, '#[Iter n.]')
      !call write_iter_header(out_gauge, '[' // trim(units_abbrev(units_out%time)) // ']')
      !call write_iter_string(out_gauge, &
      !  'A Vector potential in '   // trim(units_abbrev(units_out%length)) &
      !  'A dot in '                // trim(units_abbrev(units_out%length)) &
      !  'A dot dot in '            // trim(units_abbrev(units_out%length))
      !call write_iter_nl(out_gauge)

      call write_iter_string(out_gauge,'################################################################################')
      call write_iter_nl(out_gauge)

    end if

    call write_iter_start(out_gauge)

    do idir = 1, this%space%dim
      temp(idir) = units_from_atomic(units_out%energy, this%vecpot(idir))
    end do
    call write_iter_double(out_gauge, temp, this%space%dim)

    do idir = 1, this%space%dim
      temp(idir) = units_from_atomic(units_out%energy / units_out%time, this%vecpot_vel(idir))
    end do
    call write_iter_double(out_gauge, temp, this%space%dim)

    do idir = 1, this%space%dim
      temp(idir) = units_from_atomic(units_out%energy / units_out%time**2, this%vecpot_acc(idir))
    end do
    call write_iter_double(out_gauge, temp, this%space%dim)

    call write_iter_nl(out_gauge)

    POP_SUB(gauge_field_output_write)
  end subroutine gauge_field_output_write

  ! ---------------------------------------------------------
  subroutine gauge_field_finalize(this)
    type(gauge_field_t), intent(inout) :: this

    PUSH_SUB(gauge_field_finalize)

    call gauge_field_end(this)

    POP_SUB(gauge_field_finalize)
  end subroutine gauge_field_finalize

  ! ---------------------------------------------------------
  logical function gauge_field_update_exposed_quantities(partner, requested_time, interaction) &
    result(allowed_to_update)
    class(gauge_field_t), intent(inout) :: partner
    type(clock_t),        intent(in)    :: requested_time
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(gauge_field_update_exposed_quantities)

    ! Always allowed to update, as the external potentials are not propagated
    allowed_to_update = .true.

    call partner%clock%set_time(requested_time)

    select type (interaction)
    type is (ghost_interaction_t)
      ! Nothing to copy. We still need to check that we are at the right
      ! time for the update though!
    class default
      call partner%copy_quantities_to_interaction(interaction)
    end select

    POP_SUB(gauge_field_update_exposed_quantities)
  end function gauge_field_update_exposed_quantities

  ! ---------------------------------------------------------
  subroutine gauge_field_update_exposed_quantity(partner, iq)
    class(gauge_field_t),      intent(inout) :: partner
    integer,                   intent(in)    :: iq

    PUSH_SUB(gauge_field_update_exposed_quantities)

    ! We are only allowed to update quantities that can be updated on demand
    ASSERT(partner%quantities(iq)%updated_on_demand)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(gauge_field_update_exposed_quantities)
  end subroutine gauge_field_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine gauge_field_init_interaction_as_partner(partner, interaction)
    class(gauge_field_t), intent(in)    :: partner
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(gauge_field_init_interaction_as_partner)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(gauge_field_init_interaction_as_partner)
  end subroutine gauge_field_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine gauge_field_copy_quantities_to_interaction(partner, interaction)
    class(gauge_field_t), intent(inout) :: partner
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(gauge_field_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(gauge_field_copy_quantities_to_interaction)
  end subroutine gauge_field_copy_quantities_to_interaction

end module gauge_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
