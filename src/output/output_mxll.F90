!! Copyright (C) 2019 F. Bonafe, R. Jestaedt, H. Appel
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

module output_mxll_oct_m
  use debug_oct_m
  use external_densities_oct_m
  use hamiltonian_mxll_oct_m
  use helmholtz_decomposition_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_mxll_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m


  implicit none

  private
  public ::              &
    output_mxll_init,    &
    output_mxll

contains

  ! ---------------------------------------------------------
  subroutine output_mxll_init(outp, namespace, space)
    type(output_t),            intent(out) :: outp
    type(namespace_t),         intent(in)  :: namespace
    type(space_t),             intent(in)  :: space
    integer :: what_i
    PUSH_SUB(output_mxll_init)

    !%Variable MaxwellOutput
    !%Type block
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what to print. The output files are written at the end of the run into the output directory for the
    !% Maxwell run.
    !% Time-dependent simulations print only per iteration, including always the last. The frequency of output per iteration
    !% is set by <tt>OutputInterval</tt> and the directory is set by <tt>OutputIterDir</tt>.
    !% Each option must be in a separate row. Optionally individual output formats and output intervals can be defined
    !% for each row or they can be read separately from <tt>OutputFormat</tt> and <tt>MaxwellOutputInterval</tt> variables
    !% in the input file.
    !%
    !% Example:
    !% <br><br><tt>%MaxwellOutput
    !% <br>&nbsp;&nbsp;electric_field
    !% <br>&nbsp;&nbsp;magnetic_field
    !% <br>%<br></tt>
    !% This block supports all the formats of the <tt>Output</tt> block.
    !% See <tt>Output</tt>.
    !%Option electric_field 1
    !% Output of the electric field
    !%Option magnetic_field 2
    !% Output of the magnetic field
    !%Option trans_electric_field 3
    !% Output of the transversal electric field
    !%Option trans_magnetic_field 4
    !% Output of the transversal magnetic field
    !%Option long_electric_field 5
    !% Output of the longitudinal electric field
    !%Option long_magnetic_field 6
    !% Output of the longitudinal magnetic field
    !%Option div_electric_field 7
    !% Output of the divergence of the electric field
    !%Option div_magnetic_field 8
    !% Output of the divergence of the magnetic field
    !%Option poynting_vector 9
    !% Output of the Maxwell Poynting vector
    !%Option maxwell_energy_density 10
    !% Output of the electromagnetic density
    !%Option external_current 11
    !% Output of the external Maxwell current
    !%Option charge_density 12
    !% Output of the charge density calculated by the divergence of the electric field.
    !%Option orbital_angular_momentum 13
    !% Output of the orbital angular momentum
    !%Option vector_potential_mag 14
    !% Output of the vector potential from magnetic field
    !%Option magnetic_field_diff 15
    !% Output of the magnetic field difference
    !%Option total_current_mxll 16
    !% Output of the total current density
    !%End

    !%Variable MaxwellOutputInterval
    !%Type integer
    !%Default 50
    !%Section Output
    !%Description
    !% The output requested by variable <tt>MaxwellOutput</tt> is written
    !% to the directory <tt>MaxwellOutputIterDir</tt>
    !% when the iteration number is a multiple of the <tt>MaxwellOutputInterval</tt> variable.
    !% Subdirectories are named Y.X, where Y is <tt>td</tt>, <tt>scf</tt>, or <tt>unocc</tt>, and
    !% X is the iteration number. To use the working directory, specify <tt>"."</tt>
    !% (Output of restart files is instead controlled by <tt>MaxwellRestartWriteInterval</tt>.)
    !% Must be >= 0. If it is 0, then no output is written.
    !% This variable can also be defined inside the <tt>MaxwellOutput</tt> block.
    !% See <tt>MaxwellOutput</tt>.
    !%End

    outp%what = .false.
    call io_function_read_what_how_when(namespace, space, outp%what, outp%how, outp%output_interval, &
      'MaxwellOutput', 'OutputFormat', 'MaxwellOutputInterval')

    do what_i = lbound(outp%what, 1), ubound(outp%what, 1)
      ! xyz format is not available for Maxwell
      if (bitand(outp%how(what_i), OPTION__OUTPUTFORMAT__XYZ) /= 0) then
        message(1) = "OutputFormat = xyz is not compatible with Maxwell systems"
        call messages_warning(1)
      end if
    end do



    !%Variable MaxwellOutputIterDir
    !%Default "output_iter"
    !%Type string
    !%Section Output
    !%Description
    !% The name of the directory where <tt>Octopus</tt> stores information
    !% such as the density, forces, etc. requested by variable <tt>MaxwellOutput</tt>
    !% in the format specified by <tt>OutputHow</tt>.
    !% This information is written while iterating <tt>CalculationMode = maxwell</tt>
    !% according to <tt>OutputInterval</tt>, and has nothing to do with the restart information.
    !%End
    call parse_variable(namespace, 'MaxwellOutputIterDir', "output_iter", outp%iter_dir)
    if (any(outp%what) .and. maxval(outp%output_interval) > 0) then
      call io_mkdir(outp%iter_dir, namespace)
    end if
    call add_last_slash(outp%iter_dir)

    !%Variable MaxwellRestartWriteInterval
    !%Type integer
    !%Default 50
    !%Section Execution::IO
    !%Description
    !% Restart data is written when the iteration number is a multiple of the
    !% <tt>MaxwellRestartWriteInterval</tt> variable. (Other output is controlled by <tt>MaxwellOutputInterval</tt>.)
    !%End
    call parse_variable(namespace, 'MaxwellRestartWriteInterval', 50, outp%restart_write_interval)
    if (outp%restart_write_interval <= 0) then
      message(1) = "MaxwellRestartWriteInterval must be > 0."
      call messages_fatal(1, namespace=namespace)
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__ELECTRIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_ELECTRIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_MAGNETIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if


    POP_SUB(output_mxll_init)
  end subroutine output_mxll_init


  ! ---------------------------------------------------------
  subroutine output_mxll(outp, namespace, space, gr_mxll, st_mxll, hm_mxll, helmholtz, time, dir, gr_elec, st_elec, hm_elec)
    type(output_t),                     intent(in)    :: outp
    type(namespace_t),                  intent(in)    :: namespace
    type(space_t),                      intent(in)    :: space
    type(states_mxll_t),                intent(inout) :: st_mxll
    type(hamiltonian_mxll_t),           intent(in)    :: hm_mxll
    type(grid_t),                       intent(inout) :: gr_mxll
    type(helmholtz_decomposition_t),    intent(inout) :: helmholtz
    FLOAT,                              intent(in)    :: time
    character(len=*),                   intent(in)    :: dir
    type(grid_t),             optional, intent(inout) :: gr_elec
    type(states_elec_t),      optional, intent(inout) :: st_elec
    type(hamiltonian_elec_t), optional, intent(in)    :: hm_elec

    PUSH_SUB(output_mxll)

    if (any(outp%what)) then
      message(1) = "Info: Writing output to " // trim(dir)
      call messages_info(1, namespace=namespace)
      call io_mkdir(dir, namespace)
    end if

    call output_states_mxll(outp, namespace, space, dir, st_mxll, gr_mxll, hm_mxll)
    call output_energy_density_mxll(outp, namespace, space, dir, hm_mxll, gr_mxll)
    call output_poynting_vector_orbital_angular_momentum(outp, namespace, space, dir, st_mxll, gr_mxll, hm_mxll)
    call output_transverse_rs_state(helmholtz, outp, namespace, space, dir, st_mxll, gr_mxll)
    call output_longitudinal_rs_state(helmholtz, outp, namespace, space, dir, st_mxll, gr_mxll)
    call output_vector_potential_mag(outp, helmholtz, namespace, gr_mxll, space, dir, st_mxll)
    call output_divergence_rs_state(outp, namespace, space, dir, st_mxll, gr_mxll)
    call output_external_current_density(outp, namespace, space, dir, st_mxll, gr_mxll, hm_mxll, time)
    call output_total_current_mxll(outp, namespace, space, dir, st_mxll, gr_mxll, hm_mxll, time, st_mxll%ep, gr_mxll%np)
    call output_charge_density_mxll(outp, namespace, space, dir, st_mxll, gr_mxll, hm_mxll)

    if (present(hm_elec) .and. present(gr_elec) .and. present(st_elec)) then
      call output_coupling_potentials(outp, namespace, dir, hm_elec, gr_elec)
      call output_current_density(outp, namespace, dir, st_mxll, gr_mxll, hm_mxll, st_elec, gr_elec, hm_elec, time)
    end if

    POP_SUB(output_mxll)
  end subroutine output_mxll


  ! ---------------------------------------------------------
  subroutine output_states_mxll(outp, namespace, space, dir, st, mesh, hm)
    type(output_t),           intent(in)    :: outp
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(states_mxll_t),      intent(inout) :: st
    class(mesh_t),            intent(in)    :: mesh
    type(hamiltonian_mxll_t), intent(in)    :: hm

    integer :: ierr
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:)

    PUSH_SUB(output_states_mxll)

    ! Electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__ELECTRIC_FIELD)) then
      fn_unit = units_out%energy/units_out%length
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:mesh%box%dim))
      call get_electric_field_state(st%rs_state, mesh, dtmp, st%ep, mesh%np)
      fname = 'e_field'
      call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__ELECTRIC_FIELD), dir, fname, namespace, space, &
        mesh, dtmp, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! Magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD)) then
      fn_unit = unit_one/units_out%length**2
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:mesh%box%dim))
      call get_magnetic_field_state(st%rs_state, mesh, st%rs_sign, dtmp, st%mu(1:mesh%np), mesh%np)
      fname = 'b_field'
      call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD), dir, fname, namespace, space, &
        mesh, dtmp, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_states_mxll)
  end subroutine output_states_mxll


  !----------------------------------------------------------
  subroutine output_energy_density_mxll(outp, namespace, space, dir, hm, mesh)   !< have to set unit output correctly
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    class(mesh_t),             intent(in)    :: mesh
    type(output_t),            intent(in)    :: outp

    integer :: ierr

    PUSH_SUB(output_energy_density_mxll)

    ! Maxwell energy density
    if (outp%what(OPTION__MAXWELLOUTPUT__MAXWELL_ENERGY_DENSITY)) then
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MAXWELL_ENERGY_DENSITY), dir, "maxwell_energy_density", &
        namespace, space, mesh, hm%energy%energy_density(:), units_out%energy/units_out%length**3, ierr)
    end if

    POP_SUB(output_energy_density_mxll)
  end subroutine output_energy_density_mxll


  !----------------------------------------------------------
  subroutine output_poynting_vector_orbital_angular_momentum(outp, namespace, space, dir, st, mesh, hm)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st
    class(mesh_t),             intent(inout) :: mesh
    type(hamiltonian_mxll_t),  intent(in)    :: hm

    integer :: ierr
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: poynting_vector(:,:), oam(:,:)

    PUSH_SUB(output_poynting_vector_orbital_angular_momentum)

    if (outp%what(OPTION__MAXWELLOUTPUT__POYNTING_VECTOR).or.outp%what(OPTION__MAXWELLOUTPUT__ORBITAL_ANGULAR_MOMENTUM)) then

      ! the Poynting vector is needed in both cases, so we compute it here first
      fn_unit = units_out%energy/(units_out%time*units_out%length**2)
      SAFE_ALLOCATE(poynting_vector(1:mesh%np, 1:mesh%box%dim))
      call get_poynting_vector(mesh, st, st%rs_state, st%rs_sign, poynting_vector, st%ep, st%mu)

      if(outp%what(OPTION__MAXWELLOUTPUT__POYNTING_VECTOR)) then
        fname = 'poynting_vector'
        call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__POYNTING_VECTOR), dir, fname, namespace, space, &
          mesh, poynting_vector, fn_unit, ierr)
      end if

      if(outp%what(OPTION__MAXWELLOUTPUT__ORBITAL_ANGULAR_MOMENTUM)) then
        SAFE_ALLOCATE(oam(1:mesh%np, 1:mesh%box%dim))
        call get_orbital_angular_momentum(mesh, st, poynting_vector, oam)
        fname = 'orbital_angular_momentum'
        call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__ORBITAL_ANGULAR_MOMENTUM), dir, fname, namespace, space, &
          mesh, oam, fn_unit, ierr)
        SAFE_DEALLOCATE_A(oam)
      end if

      SAFE_DEALLOCATE_A(poynting_vector)
    end if

    POP_SUB(output_poynting_vector_orbital_angular_momentum)
  end subroutine output_poynting_vector_orbital_angular_momentum


  ! add documentation to magnetic field subtraction order
  subroutine output_vector_potential_mag(outp, helmholtz, namespace, gr, space, dir, st)    !< have to set unit output correctly
    type(output_t),                  intent(in)    :: outp
    type(helmholtz_decomposition_t), intent(inout) :: helmholtz
    type(namespace_t),               intent(in)    :: namespace
    type(grid_t),                    intent(inout) :: gr
    type(space_t),                   intent(in)    :: space
    character(len=*),                intent(in)    :: dir
    type(states_mxll_t),             intent(in)    :: st

    character(len=MAX_PATH_LEN) :: fname
    integer :: ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:),vtmp(:,:), mag_fromvec(:,:) , delta(:,:)

    PUSH_SUB(output_vector_potential_mag)

    ! vector potential from magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__VECTOR_POTENTIAL_MAG) .or. &
      outp%what(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD_DIFF)) then

      fn_unit = unit_one/units_out%length
      SAFE_ALLOCATE(dtmp(1:gr%np_part, 1:gr%box%dim))
      SAFE_ALLOCATE(vtmp(1:gr%np_part, 1:gr%box%dim))
      call get_magnetic_field_state(st%rs_state, gr, st%rs_sign, dtmp, st%mu(1:gr%np), gr%np)
      ! Get vector potential for magnetic field
      call helmholtz%get_vector_potential(namespace, vtmp, dtmp)

      if (outp%what(OPTION__MAXWELLOUTPUT__VECTOR_POTENTIAL_MAG)) then
        message(1) = 'Vector potential is currently missing surface contributions'
        message(2) = 'If in doubt, use magnetic_field_diff output which shows deviation of B field'
        message(3) = 'Large B field deviations mean the calculated vector potential is unreliable'
        call messages_warning(3)

        fname = 'vector_potential_mag'
        call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__VECTOR_POTENTIAL_MAG), dir, fname, namespace, &
          space, gr, vtmp, fn_unit, ierr)
      end if

      if (outp%what(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD_DIFF)) then
        fn_unit = unit_one/units_out%length**2
        SAFE_ALLOCATE(mag_fromvec(1:gr%np_part, 1:gr%box%dim))
        SAFE_ALLOCATE(delta(1:gr%np, 1:gr%box%dim))
        ! Get transverse field from vector potential
        call helmholtz%get_trans_field(namespace, mag_fromvec, vector_potential=vtmp)
        ! mag_fromvec is of size np_part, but should be used up to np
        ! Due to how the Helmholtz decomposition works (see its documentation), we have to remove a stencil from dtmp
        call helmholtz%apply_inner_stencil_mask(dtmp)
        delta = dtmp(1:gr%np, 1:gr%box%dim) - mag_fromvec(1:gr%np, 1:gr%box%dim)
        fname = 'magnetic_field_diff'
        call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD_DIFF), dir, fname, namespace, &
          space, gr, delta, fn_unit, ierr)
        SAFE_DEALLOCATE_A(mag_fromvec)
        SAFE_DEALLOCATE_A(delta)
      end if

      SAFE_DEALLOCATE_A(dtmp)
      SAFE_DEALLOCATE_A(vtmp)
    end if

    POP_SUB(output_vector_potential_mag)
  end subroutine output_vector_potential_mag


  !----------------------------------------------------------
  subroutine output_transverse_rs_state(helmholtz, outp, namespace, space, dir, st, gr)    !< have to set unit output correctly
    type(helmholtz_decomposition_t), intent(inout) :: helmholtz
    type(output_t),                  intent(in)    :: outp
    type(namespace_t),               intent(in)    :: namespace
    type(space_t),                   intent(in)    :: space
    character(len=*),                intent(in)    :: dir
    type(states_mxll_t),             intent(inout) :: st
    type(grid_t),                    intent(in)    :: gr

    character(len=MAX_PATH_LEN) :: fname
    integer :: ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:)

    PUSH_SUB(output_transverse_rs_state)

    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_ELECTRIC_FIELD) .or. outp%what(OPTION__MAXWELLOUTPUT__TRANS_MAGNETIC_FIELD)) then
      call helmholtz%get_trans_field(namespace, st%rs_state_trans, total_field=st%rs_state)
    end if

    ! transverse component of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_ELECTRIC_FIELD)) then
      fn_unit = units_out%energy/units_out%length
      SAFE_ALLOCATE(dtmp(1:gr%np, 1:gr%box%dim))
      call get_electric_field_state(st%rs_state_trans, gr, dtmp, st%ep(1:gr%np), gr%np)
      fname = 'e_field_trans'
      call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__TRANS_ELECTRIC_FIELD), dir, fname, namespace, space, &
        gr, dtmp, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! transverse component of the magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_MAGNETIC_FIELD)) then
      fn_unit = unit_one/units_out%length**2
      SAFE_ALLOCATE(dtmp(1:gr%np, 1:gr%box%dim))
      call get_magnetic_field_state(st%rs_state_trans, gr, st%rs_sign, dtmp, st%mu(1:gr%np), gr%np)
      fname = 'b_field_trans'
      call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__TRANS_MAGNETIC_FIELD), dir, fname, namespace, space, &
        gr, dtmp, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_transverse_rs_state)
  end subroutine output_transverse_rs_state


  !----------------------------------------------------------
  !This subroutine is left here, will be corrected and called in near future
  subroutine output_longitudinal_rs_state(helmholtz, outp, namespace, space, dir, st, gr)    !< have to set unit output correctly
    type(helmholtz_decomposition_t), intent(inout) :: helmholtz
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st
    type(grid_t),              intent(in)    :: gr

    character(len=MAX_PATH_LEN) :: fname
    integer :: ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:)

    PUSH_SUB(output_longitudinal_rs_state)

    if (outp%what(OPTION__MAXWELLOUTPUT__LONG_ELECTRIC_FIELD) .or. outp%what(OPTION__MAXWELLOUTPUT__LONG_MAGNETIC_FIELD)) then
      call helmholtz%get_long_field(namespace, st%rs_state_long, total_field=st%rs_state)
    end if

    ! longitudinal component of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__LONG_ELECTRIC_FIELD)) then
      fn_unit = units_out%energy/units_out%length
      SAFE_ALLOCATE(dtmp(1:gr%np, 1:gr%box%dim))
      call get_electric_field_state(st%rs_state_long, gr, dtmp, st%ep(1:gr%np), gr%np)
      fname = 'e_field_long'
      call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__LONG_ELECTRIC_FIELD), dir, fname, namespace, space, &
        gr, dtmp, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! longitudinal component of the magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__LONG_MAGNETIC_FIELD)) then
      fn_unit = unit_one/units_out%length**2
      SAFE_ALLOCATE(dtmp(1:gr%np, 1:gr%box%dim))
      call get_magnetic_field_state(st%rs_state_long, gr, st%rs_sign, dtmp, st%mu(1:gr%np), gr%np)
      fname = 'b_field_long'
      call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__LONG_MAGNETIC_FIELD), dir, fname, namespace, space, &
        gr, dtmp, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_longitudinal_rs_state)
  end subroutine output_longitudinal_rs_state


  !----------------------------------------------------------
  subroutine output_divergence_rs_state(outp, namespace, space, dir, st, gr)    !< have to set unit output correctly
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(in)    :: st
    type(grid_t),              intent(in)    :: gr

    integer :: ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp_1(:,:), dtmp_2(:)

    PUSH_SUB(output_divergence_rs_state)

    ! divergence of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__DIV_ELECTRIC_FIELD)) then
      fn_unit = units_out%length**(-space%dim)
      SAFE_ALLOCATE(dtmp_1(1:gr%np_part, 1:gr%box%dim))
      SAFE_ALLOCATE(dtmp_2(1:gr%np))
      dtmp_1 = M_ZERO
      call get_electric_field_state(st%rs_state, gr, dtmp_1(1:gr%np,:), st%ep(1:gr%np), gr%np)
      call get_divergence_field(gr, dtmp_1, dtmp_2, .false.)
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__DIV_ELECTRIC_FIELD), dir, "e_field_div", namespace, space, &
        gr, dtmp_2, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp_1)
      SAFE_DEALLOCATE_A(dtmp_2)
    end if

    ! divergence of the magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__DIV_MAGNETIC_FIELD)) then
      ! unit does not matter, should be zero
      fn_unit = unit_one
      SAFE_ALLOCATE(dtmp_1(1:gr%np_part, 1:gr%box%dim))
      SAFE_ALLOCATE(dtmp_2(1:gr%np))
      dtmp_1 = M_ZERO
      call get_magnetic_field_state(st%rs_state, gr, st%rs_sign, dtmp_1(1:gr%np,:), st%mu(1:gr%np), gr%np)
      call get_divergence_field(gr, dtmp_1, dtmp_2, .false.)
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__DIV_MAGNETIC_FIELD), dir, "b_field_div", namespace, space, &
        gr, dtmp_2, fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp_1)
      SAFE_DEALLOCATE_A(dtmp_2)
    end if

    POP_SUB(output_divergence_rs_state)
  end subroutine output_divergence_rs_state


  !------------------------------------------------------------
  subroutine output_coupling_potentials(outp, namespace, dir, hm, gr)    !< have to set unit output correctly
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    character(len=*),          intent(in)    :: dir
    type(hamiltonian_elec_t),  intent(in)    :: hm
    type(grid_t),              intent(in)    :: gr

    message(1) = "Maxwell-matter coupling potentials not implemented yet."
    call messages_fatal(1, namespace=namespace)

  end subroutine output_coupling_potentials


  !----------------------------------------------------------
  subroutine output_current_density(outp, namespace, dir, st_mxll, gr_mxll, hm_mxll, st_elec, gr_elec, hm_elec, time)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st_mxll
    type(grid_t),              intent(inout) :: gr_mxll
    type(hamiltonian_mxll_t),  intent(in)    :: hm_mxll
    type(states_elec_t),       intent(inout) :: st_elec
    type(grid_t),              intent(inout) :: gr_elec
    type(hamiltonian_elec_t),  intent(in)    :: hm_elec
    FLOAT,                     intent(in)    :: time

    message(1) = "Current density can not yet be calculated in a Maxwell propagation."
    call messages_fatal(1, namespace=namespace)

  end subroutine output_current_density


  !----------------------------------------------------------
  subroutine output_external_current_density(outp, namespace, space, dir, st, mesh, hm, time)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st
    class(mesh_t),             intent(in)    :: mesh
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    FLOAT,                     intent(in)    :: time

    character(len=MAX_PATH_LEN) :: fname
    integer                     :: ierr
    FLOAT, allocatable          :: dtmp(:,:)
    type(unit_t)                :: fn_unit

    PUSH_SUB(output_external_current_density)

    ! Maxwell current density
    if (outp%what(OPTION__MAXWELLOUTPUT__EXTERNAL_CURRENT)) then
      if (hm%current_density_ext_flag) then
        fn_unit = (unit_one/units_out%time)/(units_out%length**2)      !< test both if its the same
        SAFE_ALLOCATE(dtmp(1:mesh%np, 1:mesh%box%dim))
        call external_current_calculation(st, mesh, time, dtmp)
        fname = 'external_current'
        call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__EXTERNAL_CURRENT), dir, fname, namespace, space, &
          mesh, dtmp, fn_unit, ierr)
        SAFE_DEALLOCATE_A(dtmp)
      end if
    end if

    POP_SUB(output_external_current_density)
  end subroutine output_external_current_density

  subroutine output_total_current_mxll(outp, namespace, space, dir, st, mesh, hm, time, ep_field, np)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st
    class(mesh_t),             intent(in)    :: mesh
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    FLOAT,                     intent(in)    :: time
    FLOAT,   optional, intent(in)    :: ep_field(:)
    integer, optional, intent(in)    :: np

    FLOAT, allocatable          :: ctmp(:,:)
    character(len=MAX_PATH_LEN) :: fname
    integer                     :: ierr, ip, idim, np_, ff_dim
    type(unit_t)                :: fn_unit

    PUSH_SUB(output_total_current_mxll)
    SAFE_ALLOCATE(ctmp(1:mesh%np, 1:mesh%box%dim))
    if (outp%what(OPTION__MAXWELLOUTPUT__TOTAL_CURRENT_MXLL)) then
      np_ = optional_default(np, mesh%np)
      ff_dim = mesh%box%dim
      if (present(ep_field)) then
        do idim = 1, ff_dim
          do ip = 1, np_
            ctmp(ip, idim) = sqrt(M_TWO*ep_field(ip)) * real(st%rs_current_density_t1(ip, idim))
          end do
        end do
      else
        do idim = 1, ff_dim
          do ip = 1, np_
            ctmp(ip, idim) = sqrt(M_TWO*P_ep) * real(st%rs_current_density_t1(ip, idim))
          end do
        end do
      end if

      fn_unit = (unit_one/units_out%time)/(units_out%length**2)      !< test both if its the same
      fname = 'total_current_mxll'
      call io_function_output_vector(outp%how(OPTION__MAXWELLOUTPUT__TOTAL_CURRENT_MXLL), dir, fname, namespace, space, &
        mesh, ctmp, fn_unit, ierr)
    end if
    SAFE_DEALLOCATE_A(ctmp)
    POP_SUB(output_total_current_mxll)
  end subroutine output_total_current_mxll


  !----------------------------------------------------------
  subroutine output_charge_density_mxll(outp, namespace, space, dir, st, gr, hm)    !< have to set unit output correctly
    type(output_t),                intent(in)    :: outp
    type(namespace_t),             intent(in)    :: namespace
    type(space_t),                 intent(in)    :: space
    character(len=*),              intent(in)    :: dir
    type(states_mxll_t),           intent(in)    :: st
    type(grid_t),                  intent(in)    :: gr
    type(hamiltonian_mxll_t),      intent(in)    :: hm

    integer :: ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp_1(:,:), dtmp_2(:)

    PUSH_SUB(output_charge_density_mxll)

    ! charge density calculated by the divergence of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__CHARGE_DENSITY)) then
      fn_unit = units_out%length**(-space%dim)
      SAFE_ALLOCATE(dtmp_1(1:gr%np_part,1:gr%box%dim))
      SAFE_ALLOCATE(dtmp_2(1:gr%np))
      call get_electric_field_state(st%rs_state, gr, dtmp_1, st%ep, gr%np)
      call get_divergence_field(gr, dtmp_1, dtmp_2, .true.)
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__CHARGE_DENSITY), dir, "charge_density", namespace, space, &
        gr, dtmp_2(:), fn_unit, ierr)
      SAFE_DEALLOCATE_A(dtmp_1)
      SAFE_DEALLOCATE_A(dtmp_2)
    end if

    POP_SUB(output_charge_density_mxll)
  end subroutine output_charge_density_mxll

end module output_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
