!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!> this module contains the output system
module output_oct_m
  use accel_oct_m
  use basins_oct_m
  use box_oct_m
  use comm_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use current_oct_m
  use debug_oct_m
  use density_oct_m
  use derivatives_oct_m
  use dos_oct_m
  use elf_oct_m
#if defined(HAVE_ETSF_IO)
  use etsf_io
  use etsf_io_tools
#endif
  use exchange_operator_oct_m
  use ext_partner_list_oct_m
  use fft_oct_m
  use fourier_shell_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use interaction_partner_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use lda_u_io_oct_m
  use linear_response_oct_m
  use loct_oct_m
  use magnetic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use orbitalset_oct_m
  use output_me_oct_m
  use output_berkeleygw_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use smear_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_io_oct_m
  use stress_oct_m
  use string_oct_m
  use submesh_oct_m
  use symm_op_oct_m
  use symmetries_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use vtk_oct_m
  use xc_oct_m
  use xc_oep_oct_m
  use xc_f03_lib_m

  implicit none

  private
  public ::              &
    output_t,            &
    output_bgw_t,        &
    output_init,         &
    output_states,       &
    output_hamiltonian,  &
    output_all,          &
    output_current_flow, &
    doutput_lr,          &
    zoutput_lr,          &
    output_scalar_pot,   &
    output_needs_current, &
    output_need_exchange

  !> @brief output handler class
  !!
  type output_t
    private
    !> General output variables:
    logical, public    :: what(MAX_OUTPUT_TYPES)             !< what to output
    integer(i8), public :: how(0:MAX_OUTPUT_TYPES)           !< how to output

    type(output_me_t) :: me        !< this handles the output of matrix elements

    ! These variables fine-tune the output for some of the possible output options:
    integer, public :: output_interval(0:MAX_OUTPUT_TYPES)     !< output every iter
    integer, public :: restart_write_interval
    logical, public :: duringscf
    character(len=80), public :: wfs_list  !< If output_wfs, this list decides which wavefunctions to print.
    character(len=MAX_PATH_LEN), public :: iter_dir  !< The folder name, if information will be output while iterating.

    type(mesh_plane_t) :: plane    !< This is to calculate the current flow across a plane
    type(mesh_line_t)  :: line     !< or through a line (in 2D)

    type(output_bgw_t) :: bgw      !< parameters for BerkeleyGW output

  contains
    procedure :: what_now => output_what_now
    procedure :: anything_now => output_anything_now
  end type output_t

contains

  subroutine output_init(outp, namespace, space, st, nst, ks)
    type(output_t),            intent(out)   :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    type(states_elec_t),       intent(in)    :: st
    integer,                   intent(in)    :: nst
    type(v_ks_t),              intent(inout) :: ks

    type(block_t) :: blk
    FLOAT :: norm
    character(len=80) :: nst_string, default

    PUSH_SUB(output_init)
    outp%what = .false.

    call io_function_read_what_how_when(namespace, space, outp%what, outp%how, outp%output_interval)

    if (outp%what(OPTION__OUTPUT__WFS_FOURIER)) then
      if (accel_is_enabled()) then
        message(1) = "Wave functions in Fourier space not supported on GPUs."
        call messages_fatal(1, namespace=namespace)
      end if
      call messages_experimental("Wave-functions in Fourier space", namespace=namespace)
    end if

    ! cannot calculate the ELF in 1D
    if (outp%what(OPTION__OUTPUT__ELF) .or. outp%what(OPTION__OUTPUT__ELF_BASINS)) then
      if (space%dim /= 2 .and. space%dim /= 3) then
        outp%what(OPTION__OUTPUT__ELF) = .false.
        outp%what(OPTION__OUTPUT__ELF_BASINS) = .false.
        write(message(1), '(a)') 'Cannot calculate ELF except in 2D and 3D.'
        call messages_warning(1, namespace=namespace)
      end if
    end if


    if (outp%what(OPTION__OUTPUT__MMB_WFS)) then
      call messages_experimental("Model many-body wfs", namespace=namespace)
    end if

    if (outp%what(OPTION__OUTPUT__XC_TORQUE)) then
      if (st%d%ispin /= SPINORS) then
        write(message(1), '(a)') 'The output xc_torque can only be computed for spinors.'
        call messages_fatal(1, namespace=namespace)
      end if
      if (space%dim /= 3) then
        write(message(1), '(a)') 'The output xc_torque can only be computed in the 3D case.'
        call messages_fatal(1, namespace=namespace)
      end if
    end if
    if (outp%what(OPTION__OUTPUT__MMB_DEN)) then
      call messages_experimental("Model many-body density matrix", namespace=namespace)
      ! NOTES:
      !   could be made into block to be able to specify which dimensions to trace
      !   in principle all combinations are interesting, but this means we need to
      !   be able to output density matrices for multiple particles or multiple
      !   dimensions. The current 1D 1-particle case is simple.
    end if

    if (outp%what(OPTION__OUTPUT__ENERGY_DENSITY)) call messages_experimental("'Output = energy_density'", namespace=namespace)
    if (outp%what(OPTION__OUTPUT__HEAT_CURRENT)) call messages_experimental("'Output = heat_current'", namespace=namespace)

    if (outp%what(OPTION__OUTPUT__WFS) .or. outp%what(OPTION__OUTPUT__WFS_SQMOD)) then

      !%Variable OutputWfsNumber
      !%Type string
      !%Default all states
      !%Section Output
      !%Description
      !% Which wavefunctions to print, in list form: <i>i.e.</i>, "1-5" to print the first
      !% five states, "2,3" to print the second and the third state, etc.
      !% If more states are specified than available, extra ones will be ignored.
      !%End

      write(nst_string,'(i6)') nst
      write(default,'(a,a)') "1-", trim(adjustl(nst_string))
      call parse_variable(namespace, 'OutputWfsNumber', default, outp%wfs_list)
    end if

    if (parse_block(namespace, 'CurrentThroughPlane', blk) == 0) then
      if (.not. outp%what(OPTION__OUTPUT__J_FLOW)) then
        outp%what(OPTION__OUTPUT__J_FLOW) = .true.
        call parse_variable(namespace, 'OutputInterval', 50, outp%output_interval(OPTION__OUTPUT__J_FLOW))
      end if

      !%Variable CurrentThroughPlane
      !%Type block
      !%Section Output
      !%Description
      !% The code can calculate current
      !% traversing a user-defined portion of a plane, as specified by this block.
      !% A small plain-text file <tt>current-flow</tt> will be written containing this information.
      !% Only available for 1D, 2D, or 3D.
      !% In the format below, <tt>origin</tt> is a point in the plane.
      !% <tt>u</tt> and <tt>v</tt> are the (dimensionless) vectors defining the plane;
      !% they will be normalized. <tt>spacing</tt> is the fineness of the mesh
      !% on the plane. Integers <tt>nu</tt> and <tt>mu</tt> are the length and
      !% width of the portion of the plane, in units of <tt>spacing</tt>.
      !% Thus, the grid points included in the plane are
      !% <tt>x_ij = origin + i*spacing*u + j*spacing*v</tt>,
      !% for <tt>nu <= i <= mu</tt> and <tt>nv <= j <= mv</tt>.
      !% Analogously, in the 2D case, the current flow is calculated through a line;
      !% in the 1D case, the current flow is calculated through a point. Note that the spacing
      !% can differ from the one used in the main calculation; an interpolation will be performed.
      !%
      !% Example (3D):
      !%
      !% <tt>%CurrentThroughPlane
      !% <br>&nbsp;&nbsp; 0.0 | 0.0 | 0.0  # origin
      !% <br>&nbsp;&nbsp; 0.0 | 1.0 | 0.0  # u
      !% <br>&nbsp;&nbsp; 0.0 | 0.0 | 1.0  # v
      !% <br>&nbsp;&nbsp; 0.2              # spacing
      !% <br>&nbsp;&nbsp; 0 | 50           # nu | mu
      !% <br>&nbsp;&nbsp; -50 | 50         # nv | mv
      !% <br>%</tt>
      !%
      !% Example (2D):
      !%
      !% <tt>%CurrentThroughPlane
      !% <br>&nbsp;&nbsp; 0.0 | 0.0        # origin
      !% <br>&nbsp;&nbsp; 1.0 | 0.0        # u
      !% <br>&nbsp;&nbsp; 0.2              # spacing
      !% <br>&nbsp;&nbsp; 0 | 50           # nu | mu
      !% <br>%</tt>
      !%
      !% Example (1D):
      !%
      !% <tt>%CurrentThroughPlane
      !% <br>&nbsp;&nbsp; 0.0              # origin
      !% <br>%</tt>
      !%
      !%End

      select case (space%dim)
      case (3)

        call parse_block_float(blk, 0, 0, outp%plane%origin(1), units_inp%length)
        call parse_block_float(blk, 0, 1, outp%plane%origin(2), units_inp%length)
        call parse_block_float(blk, 0, 2, outp%plane%origin(3), units_inp%length)
        call parse_block_float(blk, 1, 0, outp%plane%u(1))
        call parse_block_float(blk, 1, 1, outp%plane%u(2))
        call parse_block_float(blk, 1, 2, outp%plane%u(3))
        call parse_block_float(blk, 2, 0, outp%plane%v(1))
        call parse_block_float(blk, 2, 1, outp%plane%v(2))
        call parse_block_float(blk, 2, 2, outp%plane%v(3))
        call parse_block_float(blk, 3, 0, outp%plane%spacing, units_inp%length)
        call parse_block_integer(blk, 4, 0, outp%plane%nu)
        call parse_block_integer(blk, 4, 1, outp%plane%mu)
        call parse_block_integer(blk, 5, 0, outp%plane%nv)
        call parse_block_integer(blk, 5, 1, outp%plane%mv)

        norm = norm2(outp%plane%u(1:3))
        if (norm < M_EPSILON) then
          write(message(1), '(a)') 'u-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1, namespace=namespace)
        end if
        outp%plane%u(1:3) = outp%plane%u(1:3) / norm

        norm = norm2(outp%plane%v(1:3))
        if (norm < M_EPSILON) then
          write(message(1), '(a)') 'v-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1, namespace=namespace)
        end if
        outp%plane%v(1:3) = outp%plane%v(1:3) / norm

        outp%plane%n(1) = outp%plane%u(2)*outp%plane%v(3) - outp%plane%u(3)*outp%plane%v(2)
        outp%plane%n(2) = outp%plane%u(3)*outp%plane%v(1) - outp%plane%u(1)*outp%plane%v(3)
        outp%plane%n(3) = outp%plane%u(1)*outp%plane%v(2) - outp%plane%u(2)*outp%plane%v(1)

      case (2)

        call parse_block_float(blk, 0, 0, outp%line%origin(1), units_inp%length)
        call parse_block_float(blk, 0, 1, outp%line%origin(2), units_inp%length)
        call parse_block_float(blk, 1, 0, outp%line%u(1))
        call parse_block_float(blk, 1, 1, outp%line%u(2))
        call parse_block_float(blk, 2, 0, outp%line%spacing, units_inp%length)
        call parse_block_integer(blk, 3, 0, outp%line%nu)
        call parse_block_integer(blk, 3, 1, outp%line%mu)

        norm = norm2(outp%line%u(1:2))
        if (norm < M_EPSILON) then
          write(message(1), '(a)') 'u-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1, namespace=namespace)
        end if
        outp%line%u(1:2) = outp%line%u(1:2) / norm

        outp%line%n(1) = -outp%line%u(2)
        outp%line%n(2) =  outp%line%u(1)

      case (1)

        call parse_block_float(blk, 0, 0, outp%line%origin(1), units_inp%length)

      case default

        call messages_not_implemented("CurrentThroughPlane for 4D or higher", namespace=namespace)

      end select
      call parse_block_end(blk)
    end if

    if (outp%what(OPTION__OUTPUT__MATRIX_ELEMENTS)) then
      call output_me_init(outp%me, namespace, space, st, nst)
    else
      outp%me%what = .false.
    end if

    if (outp%what(OPTION__OUTPUT__BERKELEYGW)) then
      if (accel_is_enabled()) then
        message(1) = "BerkeleyGW is not compatible with GPUs."
        call messages_fatal(1, namespace=namespace)
      end if
      call output_berkeleygw_init(nst, namespace, outp%bgw, space%periodic_dim)
    end if

    ! required for output_hamiltonian()
    if (outp%what(OPTION__OUTPUT__POTENTIAL_GRADIENT) .and. .not. outp%what(OPTION__OUTPUT__POTENTIAL)) then
      outp%what(OPTION__OUTPUT__POTENTIAL) = .true.
      outp%output_interval(OPTION__OUTPUT__POTENTIAL) = outp%output_interval(OPTION__OUTPUT__POTENTIAL_GRADIENT)
    end if


    !%Variable OutputDuringSCF
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% During <tt>gs</tt> and <tt>unocc</tt> runs, if this variable is set to yes,
    !% output will be written after every <tt>OutputInterval</tt> iterations.
    !%End
    call parse_variable(namespace, 'OutputDuringSCF', .false., outp%duringscf)

    !%Variable RestartWriteInterval
    !%Type integer
    !%Default 50
    !%Section Execution::IO
    !%Description
    !% Restart data is written when the iteration number is a multiple
    !% of the <tt>RestartWriteInterval</tt> variable. For
    !% time-dependent runs this includes the update of the output
    !% controlled by the <tt>TDOutput</tt> variable. (Other output is
    !% controlled by <tt>OutputInterval</tt>.)
    !%End
    call parse_variable(namespace, 'RestartWriteInterval', 50, outp%restart_write_interval)
    if (outp%restart_write_interval <= 0) then
      message(1) = "RestartWriteInterval must be > 0."
      call messages_fatal(1, namespace=namespace)
    end if

    !%Variable OutputIterDir
    !%Default "output_iter"
    !%Type string
    !%Section Output
    !%Description
    !% The name of the directory where <tt>Octopus</tt> stores information
    !% such as the density, forces, etc. requested by variable <tt>Output</tt>
    !% in the format specified by <tt>OutputFormat</tt>.
    !% This information is written while iterating <tt>CalculationMode = gs</tt>, <tt>unocc</tt>, or <tt>td</tt>,
    !% according to <tt>OutputInterval</tt>, and has nothing to do with the restart information.
    !%End
    call parse_variable(namespace, 'OutputIterDir', "output_iter", outp%iter_dir)
    if (any(outp%what) .and. any(outp%output_interval > 0)) then
      call io_mkdir(outp%iter_dir, namespace)
    end if
    call add_last_slash(outp%iter_dir)

    ! At this point, we don`t know whether the states will be real or complex.
    ! We therefore pass .false. to states_are_real, and need to check for real states later.

    if (output_needs_current(outp, .false.)) then
      call v_ks_calculate_current(ks, .true.)
    else
      call v_ks_calculate_current(ks, .false.)
    end if


    POP_SUB(output_init)
  end subroutine output_init

  ! ---------------------------------------------------------
  subroutine output_all(outp, namespace, space, dir, gr, ions, iter, st, hm, ks)
    type(output_t),           intent(in)    :: outp
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    integer,                  intent(in)    :: iter
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(v_ks_t),             intent(inout) :: ks

    integer :: idir, ierr, iout, iunit
    character(len=80) :: fname
    type(profile_t), save :: prof

    PUSH_SUB(output_all)
    call profiling_in(prof, "OUTPUT_ALL")

    if (any(outp%what)) then
      message(1) = "Info: Writing output to " // trim(dir)
      call messages_info(1, namespace=namespace)
      call io_mkdir(dir, namespace)
    end if

    if (outp%what_now(OPTION__OUTPUT__MESH_R, iter)) then
      do idir = 1, space%dim
        write(fname, '(a,a)') 'mesh_r-', index2axis(idir)
        call dio_function_output(outp%how(OPTION__OUTPUT__MESH_R), dir, fname, namespace, space, &
          gr, gr%x(:,idir), units_out%length, ierr, pos=ions%pos, atoms=ions%atom)
      end do
    end if

    call output_states(outp, namespace, space, dir, st, gr, ions, hm, iter)
    call output_hamiltonian(outp, namespace, space, dir, hm, st, gr%der, ions, gr, iter, st%st_kpt_mpi_grp)

    if (outp%what_now(OPTION__OUTPUT__EL_PRESSURE, iter)) then
      if(ks%theory_level /= KOHN_SHAM_DFT) then
        call messages_not_implemented("el_pressure for TheoryLevel different from kohn_sham", namespace=namespace)
      end if
      if (st%d%spin_channels > 1) then
        call messages_not_implemented("el_pressure for spin-polarized or spinors", namespace=namespace)
      end if
    end if

    call output_localization_funct(outp, namespace, space, dir, st, hm, gr, ions, iter)

    if (outp%what_now(OPTION__OUTPUT__J_FLOW, iter)) then
      call output_current_flow(outp, namespace, space, dir, gr, st, hm%kpoints)
    end if

    if (outp%what_now(OPTION__OUTPUT__GEOMETRY, iter)) then
      if (bitand(outp%how(OPTION__OUTPUT__GEOMETRY), OPTION__OUTPUTFORMAT__XCRYSDEN) /= 0) then
        call write_xsf_geometry_file(dir, "geometry", space, ions%latt, ions%pos, ions%atom, gr, namespace)
      end if
      if (bitand(outp%how(OPTION__OUTPUT__GEOMETRY), OPTION__OUTPUTFORMAT__XYZ) /= 0) then
        call ions%write_xyz(trim(dir)//'/geometry')
        if (ions%space%is_periodic()) then
          call ions%write_crystal(dir)
        end if
      end if
      if (bitand(outp%how(OPTION__OUTPUT__GEOMETRY), OPTION__OUTPUTFORMAT__VTK) /= 0) then
        call ions%write_vtk_geometry(trim(dir)//'/geometry')
      end if
    end if

    if (outp%what_now(OPTION__OUTPUT__FORCES, iter)) then
      if (bitand(outp%how(OPTION__OUTPUT__FORCES), OPTION__OUTPUTFORMAT__BILD) /= 0) then
        call ions%write_bild_forces_file(dir, "forces")
      else
        call write_xsf_geometry_file(dir, "forces", space, ions%latt, ions%pos, ions%atom, &
          gr, namespace, total_forces = ions%tot_force)
      end if
    end if

    if (outp%what_now(OPTION__OUTPUT__MATRIX_ELEMENTS, iter)) then
      call output_me(outp%me, namespace, space, dir, st, gr, ions, hm)
    end if

    do iout = lbound(outp%how, 1), ubound(outp%how, 1)
      if (bitand(outp%how(iout), OPTION__OUTPUTFORMAT__ETSF) /= 0) then
        call output_etsf(outp, namespace, space, dir, st, gr, hm%kpoints, ions, iter)
        exit
      end if
    end do

    if (outp%what_now(OPTION__OUTPUT__BERKELEYGW, iter)) then
      call output_berkeleygw(outp%bgw, namespace, space, dir, st, gr, ks, hm, ions)
    end if

    if (outp%what_now(OPTION__OUTPUT__ENERGY_DENSITY, iter)) then
      call output_energy_density(outp, namespace, space, dir, hm, ks, st, ions, gr)
    end if

    if (outp%what_now(OPTION__OUTPUT__STRESS, iter)) then
      call io_mkdir(dir, namespace)
      iunit = io_open(trim(dir)//'/stress', namespace, action='write')
      call output_stress(iunit, namespace, space%dim, st%stress_tensors)
      call io_close(iunit)
    end if

    if (hm%lda_u_level /= DFT_U_NONE) then
      if (outp%what_now(OPTION__OUTPUT__OCC_MATRICES, iter))&
        call lda_u_write_occupation_matrices(dir, hm%lda_u, st, namespace)

      if (outp%what_now(OPTION__OUTPUT__EFFECTIVEU, iter))&
        call lda_u_write_effectiveU(dir, hm%lda_u, namespace)

      if (outp%what_now(OPTION__OUTPUT__MAGNETIZATION, iter))&
        call lda_u_write_magnetization(dir, hm%lda_u, ions, gr, st, namespace)

      if (outp%what_now(OPTION__OUTPUT__LOCAL_ORBITALS, iter))&
        call output_dftu_orbitals(outp, dir, namespace, space, hm%lda_u, st, gr, ions, allocated(hm%hm_base%phase))

      if (outp%what_now(OPTION__OUTPUT__KANAMORIU, iter))&
        call lda_u_write_kanamoriU(dir, st, hm%lda_u, namespace)
    end if

    if (bitand(ks%xc_family, XC_FAMILY_OEP) /= 0 .and. ks%theory_level /= HARTREE_FOCK &
      .and. ks%theory_level /= GENERALIZED_KOHN_SHAM_DFT) then
      if (ks%oep%level == OEP_LEVEL_FULL .and. ks%oep%type == OEP_TYPE_PHOTONS) then
        if (outp%what_now(OPTION__OUTPUT__PHOTON_CORRELATOR, iter)) then
          write(fname, '(a)') 'photon_correlator'
          call dio_function_output(outp%how(OPTION__OUTPUT__PHOTON_CORRELATOR), dir, trim(fname), namespace, space, &
            gr, ks%oep%pt%correlator(:,1), units_out%length, ierr, pos=ions%pos, atoms=ions%atom)
        end if
      end if
    end if

    call output_xc_torque(outp, namespace, dir, gr, hm, st, ions, ions%space)

    call profiling_out(prof)
    POP_SUB(output_all)
  end subroutine output_all


  ! ---------------------------------------------------------
  subroutine output_localization_funct(outp, namespace, space, dir, st, hm, gr, ions, iter)
    type(output_t),           intent(in)    :: outp
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    integer,                  intent(in)    :: iter

    FLOAT, allocatable :: f_loc(:,:)
    character(len=256) :: fname
    integer :: is, ierr, imax
    type(mpi_grp_t) :: mpi_grp

    PUSH_SUB(output_localization_funct)

    mpi_grp = st%dom_st_kpt_mpi_grp

    ! if SPIN_POLARIZED, the ELF contains one extra channel: the total ELF
    imax = st%d%nspin
    if (st%d%ispin == SPIN_POLARIZED) imax = 3

    SAFE_ALLOCATE(f_loc(1:gr%np, 1:imax))

    ! First the ELF in real space
    if (outp%what_now(OPTION__OUTPUT__ELF, iter) .or. outp%what_now(OPTION__OUTPUT__ELF_BASINS, iter)) then
      ASSERT(space%dim /= 1)

      call elf_calc(space, st, gr, hm%kpoints, f_loc)

      ! output ELF in real space
      if (outp%what_now(OPTION__OUTPUT__ELF, iter)) then
        write(fname, '(a)') 'elf_rs'
        call dio_function_output(outp%how(OPTION__OUTPUT__ELF), dir, trim(fname), namespace, space, gr, &
          f_loc(:,imax), unit_one, ierr, pos=ions%pos, atoms=ions%atom, grp = mpi_grp)
        ! this quantity is dimensionless

        if (st%d%ispin /= UNPOLARIZED) then
          do is = 1, 2
            write(fname, '(a,i1)') 'elf_rs-sp', is
            call dio_function_output(outp%how(OPTION__OUTPUT__ELF), dir, trim(fname), namespace, space, gr, &
              f_loc(:, is), unit_one, ierr, pos=ions%pos, atoms=ions%atom, grp = mpi_grp)
            ! this quantity is dimensionless
          end do
        end if
      end if

      if (outp%what_now(OPTION__OUTPUT__ELF_BASINS, iter)) then
        call out_basins(f_loc(:,1), "elf_rs_basins", outp%how(OPTION__OUTPUT__ELF_BASINS))
      end if
    end if

    ! Now Bader analysis
    if (outp%what_now(OPTION__OUTPUT__BADER, iter)) then
      do is = 1, st%d%nspin
        call dderivatives_lapl(gr%der, st%rho(:,is), f_loc(:,is))
        if (st%d%nspin == 1) then
          write(fname, '(a)') 'bader'
        else
          write(fname, '(a,i1)') 'bader-sp', is
        end if
        call dio_function_output(outp%how(OPTION__OUTPUT__BADER), dir, trim(fname), namespace, space, gr, &
          f_loc(:,is), units_out%length**(-2 - space%dim), ierr, &
          pos=ions%pos, atoms=ions%atom, grp = mpi_grp)

        if (st%d%nspin == 1) then
          write(fname, '(a)') 'bader_basins'
        else
          write(fname, '(a,i1)') 'bader_basins-sp', is
        end if
        call out_basins(f_loc(:,1), fname, outp%how(OPTION__OUTPUT__BADER))
      end do
    end if

    ! Now the pressure
    if (outp%what_now(OPTION__OUTPUT__EL_PRESSURE, iter)) then
      call calc_electronic_pressure(namespace, st, hm, gr, f_loc(:,1))
      call dio_function_output(outp%how(OPTION__OUTPUT__EL_PRESSURE), dir, "el_pressure", namespace, space, gr, &
        f_loc(:,1), unit_one, ierr, pos=ions%pos, atoms=ions%atom, grp = mpi_grp)
      ! this quantity is dimensionless
    end if

    SAFE_DEALLOCATE_A(f_loc)

    POP_SUB(output_localization_funct)

  contains
    ! ---------------------------------------------------------
    subroutine out_basins(ff, filename, output_how)
      FLOAT,            intent(in)    :: ff(:)
      character(len=*), intent(in)    :: filename
      integer(i8),      intent(in)    :: output_how

      character(len=256) :: fname
      type(basins_t)     :: basins
      integer            :: iunit

      PUSH_SUB(output_localization_funct.out_basins)

      call basins_init(basins, namespace, gr)
      call basins_analyze(basins, namespace, gr, ff(:), st%rho, CNST(0.01))

      call dio_function_output(output_how, dir, trim(filename), namespace, space, gr, &
        TOFLOAT(basins%map), unit_one, ierr, pos=ions%pos, atoms=ions%atom, grp = mpi_grp)
      ! this quantity is dimensionless

      write(fname,'(4a)') trim(dir), '/', trim(filename), '.info'
      iunit = io_open(trim(fname), namespace, action = 'write')
      call basins_write(basins, gr, iunit)
      call io_close(iunit)

      call basins_end(basins)

      POP_SUB(output_localization_funct.out_basins)
    end subroutine out_basins

  end subroutine output_localization_funct


  ! ---------------------------------------------------------
  subroutine calc_electronic_pressure(namespace, st, hm, gr, pressure)
    type(namespace_t),        intent(in)    :: namespace
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(grid_t),             intent(in)    :: gr
    FLOAT,                    intent(out)   :: pressure(:)

    FLOAT, allocatable :: rho(:,:), lrho(:), tau(:,:)
    FLOAT   :: p_tf, dens
    integer :: is, ii

    PUSH_SUB(calc_electronic_pressure)

    SAFE_ALLOCATE( rho(1:gr%np_part, 1:st%d%nspin))
    SAFE_ALLOCATE(lrho(1:gr%np))
    SAFE_ALLOCATE( tau(1:gr%np, 1:st%d%nspin))

    rho = M_ZERO
    call density_calc(st, gr, rho)
    call states_elec_calc_quantities(gr, st, hm%kpoints, .false., kinetic_energy_density = tau)

    pressure = M_ZERO
    do is = 1, st%d%spin_channels
      lrho = M_ZERO
      call dderivatives_lapl(gr%der, rho(:, is), lrho)

      pressure(:) = pressure(:) + &
        tau(:, is)/M_THREE - lrho(:)/M_FOUR
    end do

    do ii = 1, gr%np
      dens = sum(rho(ii,1:st%d%spin_channels))

      p_tf = M_TWO/M_FIVE*(M_THREE*M_PI**2)**(M_TWO/M_THREE)* &
        dens**(M_FIVE/M_THREE)

      ! add XC pressure
      ! FIXME: Not correct for spinors and spin-polarized here
      pressure(ii) = pressure(ii) + (dens*hm%vxc(ii,1) - hm%energy%exchange - hm%energy%correlation)

      pressure(ii) = pressure(ii)/p_tf
      pressure(ii) = M_HALF*(M_ONE + pressure(ii)/sqrt(M_ONE + pressure(ii)**2))
    end do

    POP_SUB(calc_electronic_pressure)
  end subroutine calc_electronic_pressure


  ! ---------------------------------------------------------
  subroutine output_energy_density(outp, namespace, space, dir, hm, ks, st, ions, gr)
    type(output_t),            intent(in) :: outp
    type(namespace_t),         intent(in) :: namespace
    type(space_t),             intent(in) :: space
    character(len=*),          intent(in) :: dir
    type(hamiltonian_elec_t),  intent(in) :: hm
    type(v_ks_t),              intent(inout) :: ks
    type(states_elec_t),       intent(in) :: st
    type(ions_t),              intent(in) :: ions
    type(grid_t),              intent(in) :: gr

    integer :: is, ierr, ip
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: energy_density(:, :)
    FLOAT, allocatable :: ex_density(:)
    FLOAT, allocatable :: ec_density(:)

    PUSH_SUB(output_energy_density)

    fn_unit = units_out%energy*units_out%length**(-space%dim)
    SAFE_ALLOCATE(energy_density(1:gr%np, 1:st%d%nspin))

    ! the kinetic energy density
    call states_elec_calc_quantities(gr, st, hm%kpoints, .true., kinetic_energy_density = energy_density)

    ! the external potential energy density
    do is = 1, st%d%nspin
      do ip = 1, gr%np
        energy_density(ip, is) = energy_density(ip, is) + st%rho(ip, is)*hm%ep%vpsl(ip)
      end do
    end do

    ! the hartree energy density
    do is = 1, st%d%nspin
      do ip = 1, gr%np
        energy_density(ip, is) = energy_density(ip, is) + M_HALF*st%rho(ip, is)*hm%vhartree(ip)
      end do
    end do

    ! the XC energy density
    SAFE_ALLOCATE(ex_density(1:gr%np))
    SAFE_ALLOCATE(ec_density(1:gr%np))

    call xc_get_vxc(gr, ks%xc, st, hm%kpoints, hm%psolver, namespace, space, st%rho, st%d%ispin, &
      hm%ions%latt%rcell_volume, ex_density = ex_density, ec_density = ec_density)
    do is = 1, st%d%nspin
      do ip = 1, gr%np
        energy_density(ip, is) = energy_density(ip, is) + ex_density(ip) + ec_density(ip)
      end do
    end do

    SAFE_DEALLOCATE_A(ex_density)
    SAFE_DEALLOCATE_A(ec_density)

    select case (st%d%ispin)
    case (UNPOLARIZED)
      write(fname, '(a)') 'energy_density'
      call dio_function_output(outp%how(OPTION__OUTPUT__ENERGY_DENSITY), dir, trim(fname), namespace, space, gr, &
        energy_density(:,1), unit_one, ierr, pos=ions%pos, atoms=ions%atom, grp = st%dom_st_kpt_mpi_grp)
    case (SPIN_POLARIZED, SPINORS)
      do is = 1, 2
        write(fname, '(a,i1)') 'energy_density-sp', is
        call dio_function_output(outp%how(OPTION__OUTPUT__ENERGY_DENSITY), dir, trim(fname), namespace, space, gr, &
          energy_density(:, is), unit_one, ierr, pos=ions%pos, atoms=ions%atom, grp = st%dom_st_kpt_mpi_grp)
      end do
    end select
    SAFE_DEALLOCATE_A(energy_density)

    POP_SUB(output_energy_density)
  end subroutine output_energy_density

  !--------------------------------------------------------------

  logical function output_need_exchange(outp) result(need_exx)
    type(output_t),         intent(in)    :: outp

    need_exx =(outp%what(OPTION__OUTPUT__BERKELEYGW) &
      .or. outp%me%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY) &
      .or. outp%me%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY_EXC_K))
  end function output_need_exchange


  ! ---------------------------------------------------------
  subroutine output_dftu_orbitals(outp, dir, namespace, space, this, st, mesh, ions, has_phase)
    type(output_t),      intent(in) :: outp
    character(len=*),    intent(in) :: dir
    type(namespace_t),   intent(in) :: namespace
    type(space_t),       intent(in) :: space
    type(lda_u_t),       intent(in) :: this
    type(states_elec_t), intent(in) :: st
    class(mesh_t),       intent(in) :: mesh
    type(ions_t),        intent(in) :: ions
    logical,             intent(in) :: has_phase

    integer :: ios, im, ik, idim, ierr
    CMPLX, allocatable :: tmp(:)
    FLOAT, allocatable :: dtmp(:)
    type(orbitalset_t), pointer :: os
    type(unit_t) :: fn_unit
    character(len=32) :: fname

    PUSH_SUB(output_dftu_orbitals)

    fn_unit = sqrt(units_out%length**(-space%dim))

    if (this%basis%use_submesh) then
      if (states_are_real(st)) then
        SAFE_ALLOCATE(dtmp(1:mesh%np))
      else
        SAFE_ALLOCATE(tmp(1:mesh%np))
      end if
    end if

    do ios = 1, this%norbsets
      os => this%orbsets(ios)
      do ik = st%d%kpt%start, st%d%kpt%end
        do im = 1, this%orbsets(ios)%norbs
          do idim = 1, min(os%ndim, st%d%dim)
            if (st%d%nik > 1) then
              if (min(os%ndim, st%d%dim) > 1) then
                write(fname, '(a,i1,a,i3.3,a,i8.8,a,i1)') 'orb', im, '-os', ios, '-k', ik, '-sp', idim
              else
                write(fname, '(a,i1,a,i3.3,a,i8.8)') 'orb', im, '-os', ios, '-k', ik
              end if
            else
              if (min(os%ndim, st%d%dim) > 1) then
                write(fname, '(a,i1,a,i3.3,a,i1)') 'orb', im, '-os', ios, '-sp', idim
              else
                write(fname, '(a,i1,a,i3.3)') 'orb', im, '-os', ios
              end if
            end if
            if (has_phase) then
              if (.not. this%basis%use_submesh) then
                call zio_function_output(outp%how(OPTION__OUTPUT__LOCAL_ORBITALS), dir, fname, namespace, space, &
                  mesh, os%eorb_mesh(1:mesh%np,im,idim,ik), fn_unit, ierr, pos=ions%pos, atoms=ions%atom)
              else
                tmp = M_Z0
                call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np,idim,im,ik), tmp)
                call zio_function_output(outp%how(OPTION__OUTPUT__LOCAL_ORBITALS), dir, fname, namespace, space, &
                  mesh, tmp, fn_unit, ierr, pos=ions%pos, atoms=ions%atom)
              end if
            else
              if (.not. this%basis%use_submesh) then
                if (states_are_real(st)) then
                  call dio_function_output(outp%how(OPTION__OUTPUT__LOCAL_ORBITALS), dir, fname, namespace, space, mesh, &
                    os%dorb(1:mesh%np,idim,im), fn_unit, ierr, pos=ions%pos, atoms=ions%atom)
                else
                  call zio_function_output(outp%how(OPTION__OUTPUT__LOCAL_ORBITALS), dir, fname, namespace, space, mesh, &
                    os%zorb(1:mesh%np,idim,im), fn_unit, ierr, pos=ions%pos, atoms=ions%atom)
                end if
              else
                if (states_are_real(st)) then
                  dtmp = M_Z0
                  call submesh_add_to_mesh(os%sphere, os%dorb(1:os%sphere%np,idim,im), dtmp)
                  call dio_function_output(outp%how(OPTION__OUTPUT__LOCAL_ORBITALS), dir, fname, namespace, space, &
                    mesh, dtmp, fn_unit, ierr, pos=ions%pos, atoms=ions%atom)
                else
                  tmp = M_Z0
                  call submesh_add_to_mesh(os%sphere, os%zorb(1:os%sphere%np,idim,im), tmp)
                  call zio_function_output(outp%how(OPTION__OUTPUT__LOCAL_ORBITALS), dir, fname, namespace, space, &
                    mesh, tmp, fn_unit, ierr, pos=ions%pos, atoms=ions%atom)
                end if
              end if
            end if
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(dtmp)

    POP_SUB(output_dftu_orbitals)
  end subroutine output_dftu_orbitals

  ! ---------------------------------------------------------
  logical function output_needs_current(outp, states_are_real)
    type(output_t),      intent(in) :: outp
    logical,             intent(in) :: states_are_real

    output_needs_current = .false.

    if (outp%what(OPTION__OUTPUT__CURRENT) &
      .or. outp%what(OPTION__OUTPUT__HEAT_CURRENT) &
      .or. outp%what(OPTION__OUTPUT__CURRENT_KPT)) then
      if (.not. states_are_real) then
        output_needs_current = .true.
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1)
      end if
    end if


  end function

  ! ---------------------------------------------------------

  function output_what_now(this, what_id, iter) result(write_now)
    class(output_t), intent(in) :: this
    integer(i8),     intent(in) :: what_id
    integer,         intent(in) :: iter

    logical :: write_now

    write_now = .false.
    if ((what_id > 0) .and. (this%output_interval(what_id) > 0)) then
      if (this%what(what_id) .and. (iter == -1 .or. mod(iter, this%output_interval(what_id)) == 0)) then
        write_now = .true.
      end if
    end if

  end function output_what_now

  ! ---------------------------------------------------------
  !> return true if any output is to be written now
  logical function output_anything_now(this, iter)
    class(output_t), intent(in) :: this
    integer,         intent(in) :: iter

    integer(i8) :: what_it

    PUSH_SUB(output_anything_now)

    output_anything_now = .false.
    do what_it = lbound(this%output_interval, 1), ubound(this%output_interval, 1)
      if (this%what_now(what_it, iter)) then
        output_anything_now = .true.
        exit
      end if
    end do

    POP_SUB(output_anything_now)
  end function output_anything_now

#include "output_etsf_inc.F90"

#include "output_states_inc.F90"

#include "output_h_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "output_linear_response_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "output_linear_response_inc.F90"

end module output_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
