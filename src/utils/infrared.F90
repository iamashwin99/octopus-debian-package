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

program infrared
  use command_line_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  integer :: iunit, ierr, ii, jj, iter, read_iter, max_iter, ini_iter, end_iter
  FLOAT :: start_time, end_time
  FLOAT, allocatable :: time(:), dipole(:,:)
  CMPLX, allocatable :: ftdipole(:,:)
  type(space_t)     :: space
  FLOAT :: ww, irtotal
  FLOAT :: dw, max_energy
  integer :: ifreq, idir
  integer, parameter :: max_freq = 10000

  ! Initialize stuff
  call global_init(is_serial = .true.)

  call getopt_init(ierr)

  call getopt_end()

  call parser_init()

  call messages_init()

  call io_init()
  call profiling_init(global_namespace)

  call unit_system_init(global_namespace)

  !These variables are documented in src/td/spectrum.F90
  call parse_variable(global_namespace, 'TDMaxSteps', 1500, max_iter)
  call parse_variable(global_namespace, 'PropagationSpectrumStartTime',  M_ZERO, start_time, units_inp%time)
  call parse_variable(global_namespace, 'PropagationSpectrumEndTime',  -M_ONE, end_time, units_inp%time)
  call parse_variable(global_namespace, 'PropagationSpectrumMaxEnergy', &
    units_from_atomic(units_inp%energy, units_to_atomic(unit_invcm, CNST(10000.0))), max_energy, units_inp%energy)

  dw = max_energy/(max_freq-M_ONE) !Initializes the wavevector step dw

  if (end_time < M_ZERO) end_time = huge(end_time)

  call space_init(space, global_namespace)

  SAFE_ALLOCATE(dipole(0:max_iter+1, 1:3))

  call read_dipole(dipole)

  SAFE_ALLOCATE(ftdipole(1:max_freq, 1:3))

  do idir = 1, 3
    call fourier(dipole(:, idir), ftdipole(:, idir))
  end do

  !and print the spectrum
  iunit = io_open('td.general/infrared', global_namespace, action='write')

100 FORMAT(100('#'))

  write(unit = iunit, iostat = ierr, fmt = 100)
  write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
  write(unit = iunit, iostat = ierr, fmt = '(a25,a1,a1,a10)') &
    '# all absorptions are in ',units_out%length%abbrev,'*',units_out%time%abbrev
  write(unit = iunit, iostat = ierr, fmt = '(a1)') '#'
  write(unit = iunit, iostat = ierr, fmt = '(a19,41x,a17)') '#        Energy    ', 'absorption'
  write(unit = iunit, iostat = ierr, fmt = '(a15,13x,a5,15x,a7,13x,a7,13x,a7)') &
    '#        [1/cm]','total','FT(<x>)','FT(<y>)','FT(<z>)'
  write(unit = iunit, iostat = ierr, fmt = 0100)

  do ifreq = 1, max_freq
    ww = dw * ifreq
    irtotal = norm2(abs(ftdipole(ifreq, 1:3)))
    write(unit = iunit, iostat = ierr, fmt = '(5e20.10)') &
      units_from_atomic(unit_invcm, ww), units_from_atomic(units_out%length*units_out%time, ww*irtotal), &
      (units_from_atomic(units_out%length*units_out%time,  ww*abs(ftdipole(ifreq, idir))), idir = 1, 3)
  end do
  call io_close(iunit)

  call profiling_end(global_namespace)
  call io_end()
  call messages_end()

  call parser_end()
  call global_end()

contains

  subroutine read_dipole(dipole)
    FLOAT,   intent(out)   :: dipole(0:, :)

    FLOAT :: charge

    PUSH_SUB(read_dipole)

    ! Opens the coordinates files.
    iunit = io_open('td.general/multipoles', global_namespace, action='read')

    call io_skip_header(iunit)

    ini_iter = -1
    iter = 0

    SAFE_ALLOCATE(time(1:max_iter))

    do while(.true.)

      read(unit = iunit, iostat = ierr, fmt = *) read_iter, time(iter), &
        charge, dipole(iter, 1), dipole(iter, 2), dipole(iter, 3)

      time(iter) =  units_to_atomic(units_out%time, time(iter))

      !dipole moment has unit of charge*length, charge has the same unit in both systems
      do ii = 1, 3
        dipole(iter, ii) = units_to_atomic(units_out%length, dipole(iter, ii))
      end do

      if (ierr /= 0) then
        iter = iter - 1 !last iteration is not valid
        exit
      end if

      ASSERT(iter == read_iter)

      if (time(iter) >= end_time) exit

      if (time(iter) >= start_time .and. ini_iter == -1) then
        ini_iter = iter
      end if

      iter = iter + 1
    end do

    call io_close(iunit)

    if (ini_iter == 0) ini_iter = 1
    end_iter = iter - 1

    write (message(1), '(a)') "Read dipole moment from '"// &
      trim(io_workpath('td.general/multipoles', global_namespace))//"'."
    call messages_info(1)

    POP_SUB(read_dipole)
  end subroutine read_dipole

  ! -------------------------------------------------

  subroutine fourier(fi, ftfi)
    implicit none
    FLOAT, intent(inout)  :: fi(:)
    CMPLX, intent(out)    :: ftfi(:)

    FLOAT :: ww, av, dt
    integer :: ifreq, count

    PUSH_SUB(fourier)

    !apply an envelope
    do jj = ini_iter, end_iter
      fi(jj) = fi(jj) * sin((time(jj)-time(ini_iter+1))*M_PI/(time(end_iter)-time(ini_iter)))
    end do

    dt = time(2) - time(1)

    !remove the DC component
    av = M_ZERO
    count = 0
    do jj = ini_iter, end_iter
      av = av + fi(jj)*dt
      count = count + 1
    end do

    do jj = ini_iter, end_iter
      fi(jj) = fi(jj) - av/(dt*count)
    end do

    write (message(1), '(a)') "Taking the Fourier transform."
    call messages_info(1)

    !now calculate the FT
    !$omp parallel do private(ww, jj)
    do ifreq = 1, max_freq
      ww = dw * ifreq
      ftfi(ifreq) = M_ZERO
      do jj = ini_iter, end_iter
        ftfi(ifreq) = ftfi(ifreq) + &
          exp(M_zI * ww * time(jj))*fi(jj)*dt
      end do
    end do
    !$omp end parallel do

    write (message(1), '(a)') "Done."
    call messages_info(1)

    POP_SUB(fourier)
  end subroutine fourier

end program infrared

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
