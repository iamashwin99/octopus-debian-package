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

program harmonic_spectrum
  use command_line_oct_m
  use fft_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use spectrum_oct_m
  use unit_system_oct_m

  implicit none

  integer :: mode, ierr, ar
  FLOAT :: w0, vec(1:3)
  type(spectrum_t) :: spectrum
  character :: pol
  logical :: get_maxima

  integer, parameter :: &
    HS_FROM_MULT = 1,   &
    HS_FROM_ACC  = 2,   &
    HS_FROM_CURR = 3

  call getopt_init(ierr)
  if (ierr /= 0) then
    call messages_write("Your Fortran compiler doesn't support command-line arguments;")
    call messages_new_line()
    call messages_write("the oct-harmonic-spectrum command is not available.")
    call messages_fatal()
  end if
  ! These are the default values.
  get_maxima = .true.
  w0 = M_ZERO
  pol = 'x'
  mode = 1
  ar = 0
  vec = (/M_ONE, M_ZERO, M_ONE /)
  call getopt_harmonic_spectrum(w0, mode, ar, vec(1),vec(2),vec(3),pol)
  if (w0 <= M_ZERO) get_maxima = .false.
  call getopt_end()


  ! Initialize stuff
  call global_init(is_serial = .true.)

  call parser_init()

  call messages_init()

  call io_init()
  call profiling_init(global_namespace)
  call unit_system_init(global_namespace)
  call fft_all_init(global_namespace)

  call spectrum_init(spectrum, global_namespace)

  call messages_obsolete_variable(global_namespace, 'HarmonicSpectrumPolarization')
  call messages_obsolete_variable(global_namespace, 'HarmonicSpectrumMode')

  if ((pol /= 'x') .and. &
    (pol /= 'y') .and. &
    (pol /= 'z') .and. &
    (pol /= '+') .and. &
    (pol /= '-') .and. &
    (pol /= 'v')) then
    message(1) = 'The polarization direction given in the command line is not valid.'
    call messages_fatal(1)
  end if

  select case (mode)
  case (HS_FROM_MULT)
    if (get_maxima) then
      call spectrum_hs_from_mult(spectrum, global_namespace, 'hs-mult-maxima', pol, vec, w0)
    else
      if (ar == 1) then
        message(1)= "Calculating angle-resolved hs from multipoles."
        call messages_info(1)
        call spectrum_hs_ar_from_mult(spectrum, global_namespace, 'hs-mult', vec)
      else
        call spectrum_hs_from_mult(spectrum, global_namespace, 'hs-mult', pol, vec)
      end if
    end if
  case (HS_FROM_ACC)
    if (get_maxima) then
      call spectrum_hs_from_acc(spectrum, global_namespace, 'hs-acc-maxima', pol, vec, w0)
    else
      if (ar == 1) then
        message(1)= "Calculating angle-resolved hs from acceleration."
        call messages_info(1)
        call spectrum_hs_ar_from_acc(spectrum, global_namespace, 'hs-acc', vec)
      else
        call spectrum_hs_from_acc(spectrum, global_namespace, 'hs-acc', pol, vec)
      end if
    end if
  case (HS_FROM_CURR)
    if (get_maxima) then
      call spectrum_hs_from_mult(spectrum, global_namespace, 'hs-curr-maxima', pol, vec, w0)
    else
      call spectrum_hs_from_current(spectrum, global_namespace, 'hs-curr', pol, vec)
    end if
  case default
    message(1) = 'The harmonic-spectrum mode given in the command line is not valid.'
    call messages_fatal(1)
  end select

  call profiling_end(global_namespace)
  call io_end()
  call messages_end()

  call parser_end()

  call global_end()
end program harmonic_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
