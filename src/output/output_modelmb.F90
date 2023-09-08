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

module output_modelmb_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use messages_oct_m
  use mesh_oct_m
  use modelmb_density_matrix_oct_m
  use namespace_oct_m
  use output_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use young_oct_m


  implicit none

  private
  public ::              &
    output_modelmb

contains

  ! ---------------------------------------------------------
  subroutine output_modelmb(outp, namespace, space, dir, gr, ions, iter, st)
    type(output_t),           intent(in)    :: outp
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    integer,                  intent(in)    :: iter
    type(states_elec_t),      intent(inout) :: st

    type(profile_t), save :: prof

    PUSH_SUB(output_modelmb)
    call profiling_in(prof, "OUTPUT_MODELMB")

    if (outp%what_now(OPTION__OUTPUT__MMB_DEN, iter) .or. outp%what_now(OPTION__OUTPUT__MMB_WFS, iter)) then
      if (states_are_real(st)) then
        call doutput_modelmb(outp, namespace, space, trim(dir), gr, st, ions)
      else
        call zoutput_modelmb(outp, namespace, space, trim(dir), gr, st, ions)
      end if
    end if

    call profiling_out(prof)
    POP_SUB(output_modelmb)
  end subroutine output_modelmb

#include "undef.F90"
#include "complex.F90"
#include "output_modelmb_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "output_modelmb_inc.F90"

end module output_modelmb_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
