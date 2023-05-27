!! Copyright (C) 2017 Johannes Flick
!! Copyright (C) 2019 F. Buchholz, M. Oliveira
!! Copyright (C) 2020 S. Ohlmann
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

module photon_mode_oct_m
  use comm_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                       &
    photon_mode_t,                &
    photon_mode_init,             &
    photon_mode_end,              &
    photon_mode_write_info,       &
    photon_mode_add_poisson_terms

  type photon_mode_t
    ! All components are public by default
    integer               :: nmodes             !< Number of photon modes
    integer               :: dim                !< Dimensionality of the electronic system
    FLOAT, allocatable    :: omega(:)           !< Mode frequencies
    FLOAT, allocatable    :: lambda(:)          !< Interaction strength
    FLOAT, allocatable    :: pol(:,:)           !< Polarization of the photon field
    FLOAT, allocatable    :: pol_dipole(:,:)    !< Polarization*dipole operator
    FLOAT                 :: ex                 !< Photon exchange energy
    FLOAT, allocatable    :: number(:)          !< Number of photons in mode
    FLOAT, allocatable    :: correlator(:,:)    !< Correlation function <n(r)(ad+a)>
    FLOAT                 :: n_electrons        !< Number of electrons
    FLOAT, pointer        :: pt_coord_q0(:) => null()   !< Photon coordinates, initial value or gs result
    FLOAT, pointer        :: pt_momen_p0(:) => null()   !< Photon momenta, initial value or gs result
    FLOAT                 :: mu
    logical               :: has_q0_p0
  end type photon_mode_t

contains

  ! ---------------------------------------------------------
  subroutine photon_mode_init(this, namespace, mesh, dim, n_electrons)
    type(photon_mode_t),  intent(out) :: this
    type(namespace_t),    intent(in)  :: namespace
    type(mesh_t),         intent(in)  :: mesh
    integer,              intent(in)  :: dim
    FLOAT,                intent(in)  :: n_electrons

    type(block_t)         :: blk
    integer               :: ii, ip, idir, iunit, ncols
    logical               :: file_exists
    character(MAX_PATH_LEN) :: filename

    PUSH_SUB(photon_mode_init)

    this%nmodes = 0

    this%dim = dim
    this%has_q0_p0 = .false.
    this%n_electrons = n_electrons

    !%Variable PhotonmodesFilename
    !%Type string
    !%Default "photonmodes"
    !%Section Linear Response::Casida
    !%Description
    !% Filename for photon modes in text format
    !%  - first line contains 2 integers: number of photon modes and number of
    !%    columns
    !%  - each further line contains the given number of floats for one photon
    !%    mode
    !%End
    call parse_variable(namespace, 'PhotonmodesFilename', 'photonmodes', filename)
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
      if (mpi_grp_is_root(mpi_world)) then
        message(1) = 'Opening '//trim(filename)
        call messages_info(1, namespace=namespace)
        ! open file on root
        iunit = io_open(trim(filename), namespace, action='read', form='formatted')

        ! get dimensions from first line
        read(iunit, *) this%nmodes, ncols

        write(message(1), '(3a,i7,a,i3,a)') 'Reading file ', trim(filename), ' with ', &
          this%nmodes, ' photon modes and ', ncols, ' columns.'
        call messages_info(1, namespace=namespace)

        SAFE_ALLOCATE(this%omega(1:this%nmodes))
        SAFE_ALLOCATE(this%lambda(1:this%nmodes))
        SAFE_ALLOCATE(this%pol(1:this%nmodes,1:3))
        SAFE_ALLOCATE(this%pol_dipole(1:mesh%np, 1:this%nmodes))

        ! now read in all modes
        do ii = 1, this%nmodes
          if (ncols == 5) then
            read(iunit, *) this%omega(ii), this%lambda(ii), &
              this%pol(ii,1), this%pol(ii,2), this%pol(ii,3)
          else if (ncols == 7) then
            this%has_q0_p0 = .true.
            if (.not. associated(this%pt_coord_q0)) then
              SAFE_ALLOCATE(this%pt_coord_q0(1:this%nmodes))
            end if
            if (.not. associated(this%pt_momen_p0)) then
              SAFE_ALLOCATE(this%pt_momen_p0(1:this%nmodes))
            end if
            read(iunit, *) this%omega(ii), this%lambda(ii), &
              this%pol(ii,1), this%pol(ii,2), this%pol(ii,3), &
              this%pt_coord_q0(ii), this%pt_momen_p0(ii)
          else
            ! error if not 5 columns
            message(1) = 'Error: unexpected number of columns in file:'
            message(2) = filename
            call messages_fatal(2, namespace=namespace)
          end if

          ! Normalize polarization vector
          this%pol(ii,:) = this%pol(ii,:)/norm2(this%pol(ii,:))

          ! Calculate polarization dipole
          do ip = 1, mesh%np
            this%pol_dipole(ip, ii) = dot_product(this%pol(ii, 1:this%dim), mesh%x(ip, 1:this%dim))
          end do

        end do
        call io_close(iunit)
      end if
      ! broadcast first array dimensions, then allocate and broadcast arrays
      call mpi_world%bcast(this%nmodes, 1, MPI_INTEGER, 0)
      call mpi_world%bcast(ncols, 1, MPI_INTEGER, 0)
      if (.not. mpi_grp_is_root(mpi_world)) then
        SAFE_ALLOCATE(this%omega(1:this%nmodes))
        SAFE_ALLOCATE(this%lambda(1:this%nmodes))
        SAFE_ALLOCATE(this%pol(1:this%nmodes,1:3))
        SAFE_ALLOCATE(this%pol_dipole(1:mesh%np, 1:this%nmodes))
        if (ncols == 7) then
          SAFE_ALLOCATE(this%pt_coord_q0(1:this%nmodes))
          SAFE_ALLOCATE(this%pt_momen_p0(1:this%nmodes))
        end if
      end if
      call mpi_world%bcast(this%omega(1), this%nmodes, MPI_FLOAT, 0)
      call mpi_world%bcast(this%lambda(1), this%nmodes, MPI_FLOAT, 0)
      call mpi_world%bcast(this%pol(1,1), this%nmodes*3, MPI_FLOAT, 0)
      call mpi_world%bcast(this%pol_dipole(1,1), mesh%np*this%nmodes, MPI_FLOAT, 0)
      if (ncols == 7) then
        call mpi_world%bcast(this%pt_coord_q0(1), this%nmodes, MPI_FLOAT, 0)
        call mpi_world%bcast(this%pt_momen_p0(1), this%nmodes, MPI_FLOAT, 0)
      end if
    else
      if (.not. parse_is_defined(namespace, 'PhotonModes')) then
        call messages_write('You need to specify the correct external photon modes file')
        call messages_write('or define the PhotonModes variable!')
        call messages_fatal(namespace=namespace)
      end if
    end if

    !%Variable PhotonModes
    !%Type block
    !%Section Hamiltonian::XC
    !%Description
    !% Each line of the block should specify one photon mode. The syntax is the following:
    !%
    !% %PhotonModes
    !%  omega1 | lambda1| PolX1 | PolY1 | PolZ1
    !%  ...
    !% %
    !%
    !% The first column is the mode frequency, in units of energy.
    !% The second column is the coupling strength, in units of energy.
    !% The remaining columns specify the polarization direction of the mode.
    !% If the polarization vector should be normalized to one. If that is not the case
    !% the code will normalize it.
    !%End

    if (.not. file_exists) then
      if (parse_block(namespace, 'PhotonModes', blk) == 0) then

        this%nmodes = parse_block_n(blk)
        SAFE_ALLOCATE(this%omega(1:this%nmodes))
        SAFE_ALLOCATE(this%lambda(1:this%nmodes))
        SAFE_ALLOCATE(this%pol(1:this%nmodes, 1:3))
        SAFE_ALLOCATE(this%pol_dipole(1:mesh%np, 1:this%nmodes))

        this%pol = M_ZERO

        do ii = 1, this%nmodes
          ncols = parse_block_cols(blk, ii-1)

          ! Read line
          call parse_block_float(blk, ii-1, 0, this%omega(ii), units_inp%energy)  ! frequency
          call parse_block_float(blk, ii-1, 1, this%lambda(ii), units_inp%energy) ! coupling strength
          do idir = 1, this%dim
            call parse_block_float(blk, ii-1, idir + 1, this%pol(ii, idir)) ! polarization vector components
          end do

          if (ncols > 5) then
            this%has_q0_p0 = .true.
            if (.not. associated(this%pt_coord_q0)) then
              SAFE_ALLOCATE(this%pt_coord_q0(1:this%nmodes))
            end if
            if (.not. associated(this%pt_momen_p0)) then
              SAFE_ALLOCATE(this%pt_momen_p0(1:this%nmodes))
            end if
            call parse_block_float(blk, ii-1, 5, this%pt_coord_q0(ii))   !row, column
            call parse_block_float(blk, ii-1, 6, this%pt_momen_p0(ii))   !row, column
          end if
          ! Sanity check
          if ((ncols /= this%dim + 2) .and. (ncols /= this%dim + 4)) then
            call messages_input_error(namespace, 'PhotonModes', 'Incorrect number of columns')
          end if

          ! Normalize polarization vector
          this%pol(ii,:) = this%pol(ii,:)/norm2(this%pol(ii,:))

          ! Calculate polarization dipole
          do ip = 1, mesh%np
            this%pol_dipole(ip, ii) = dot_product(this%pol(ii, 1:this%dim), mesh%x(ip, 1:this%dim))
          end do

        end do
        call parse_block_end(blk)
      else
        call messages_write('You need to specify the photon modes!')
        call messages_fatal(namespace=namespace)
      end if
    end if

    !%Variable TDPhotonicTimeScale
    !%Type float
    !%Default 1.0
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable defines the factor between the timescale of photonic
    !% and electronic movement.
    !% for more details see the documentation of TDIonicTimeScale
    !% If you also use TDIonicTimeScale, we advise to set
    !% TDPhotonicTimeScale = TDIonicTimeScale, in the case the
    !% photon frequency is in a vibrational energy range.
    !% Important: The electronic time step will be the value of
    !% <tt>TDTimeStep</tt> divided by this variable, so if you have determined an
    !% optimal electronic time step (that we can call <i>dte</i>), it is
    !% recommended that you define your time step as:
    !%
    !% <tt>TDTimeStep</tt> = <i>dte</i> * <tt>TDPhotonicTimeScale</tt>
    !%
    !% so you will always use the optimal electronic time step
    !% (<a href=http://arxiv.org/abs/0710.3321>more details</a>).
    !%End

    call parse_variable(namespace, 'TDPhotonicTimeScale', CNST(1.0), this%mu)

    this%ex = M_ZERO
    SAFE_ALLOCATE(this%number(1:this%nmodes))
    this%number = M_ZERO

    SAFE_ALLOCATE(this%correlator(1:mesh%np, 1:this%nmodes))
    this%correlator = M_ZERO

    POP_SUB(photon_mode_init)
  end subroutine photon_mode_init

  ! ---------------------------------------------------------
  subroutine photon_mode_end(this)
    type(photon_mode_t), intent(inout) :: this

    PUSH_SUB(photon_mode_end)

    SAFE_DEALLOCATE_A(this%correlator)

    SAFE_DEALLOCATE_A(this%omega)
    SAFE_DEALLOCATE_A(this%lambda)
    SAFE_DEALLOCATE_A(this%number)

    SAFE_DEALLOCATE_A(this%pol)
    SAFE_DEALLOCATE_A(this%pol_dipole)

    if (this%has_q0_p0) then
      SAFE_DEALLOCATE_P(this%pt_coord_q0)
      SAFE_DEALLOCATE_P(this%pt_momen_p0)
    end if

    POP_SUB(photon_mode_end)
  end subroutine photon_mode_end

  !-----------------------------------------------------------------
  subroutine photon_mode_write_info(this, iunit, namespace)
    type(photon_mode_t),           intent(in) :: this
    integer,             optional, intent(in) :: iunit
    type(namespace_t),   optional, intent(in) :: namespace

    integer :: im, idir

    PUSH_SUB(photon_mode_write_info)

    call messages_print_stress(msg="Photon Modes", iunit=iunit, namespace=namespace)
    write(iunit, '(6x,a,10x,a,3x)', advance='no') 'Omega', 'Lambda'
    write(iunit, '(1x,a,i1,a)') ('Pol.(', idir, ')', idir = 1, this%dim)
    do im = 1, this%nmodes
      write(iunit, '(1x,f14.12)', advance='no') this%omega(im)
      write(iunit, '(1x,f14.12)', advance='no') this%lambda(im)
      write(iunit, '(2X,f5.3,1X)') (this%pol(im, idir), idir = 1, this%dim)
    end do
    call messages_print_stress(iunit=iunit, namespace=namespace)

    POP_SUB(photon_mode_write_info)
  end subroutine photon_mode_write_info

  !-----------------------------------------------------------------
  subroutine photon_mode_add_poisson_terms(this, mesh, rho, pot)
    type(photon_mode_t), intent(in)    :: this
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: rho(:)
    FLOAT,               intent(inout) :: pot(:)

    integer :: ip
    FLOAT :: lx, ld, dipole(mesh%box%dim)

    PUSH_SUB(photon_mode_add_poisson_terms)

    ! Currently this only works with one photon mode
    ASSERT(this%nmodes == 1)

    call dmf_dipole(mesh, rho, dipole)
    ld = dot_product(this%pol(1, 1:this%dim), dipole(1:this%dim))*this%lambda(1)

    do ip = 1, mesh%np
      lx = this%pol_dipole(ip, 1)*this%lambda(1)
      pot(ip) = pot(ip) - this%omega(1)/sqrt(this%n_electrons)*(mesh%x(ip, this%dim + 1)*ld + lx*dipole(this%dim + 1)) + lx*ld
    end do

    POP_SUB(photon_mode_add_poisson_terms)
  end subroutine photon_mode_add_poisson_terms

end module photon_mode_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
