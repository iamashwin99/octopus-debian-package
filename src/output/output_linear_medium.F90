!! Copyright (C) 2021 F. Bonafe
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

module output_linear_medium_oct_m
  use debug_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use io_function_oct_m
  use messages_oct_m
  use mesh_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                    &
    output_linear_medium_init, &
    output_linear_medium

contains

  ! ---------------------------------------------------------
  subroutine output_linear_medium_init(outp, namespace, space)
    type(output_t),            intent(out) :: outp
    type(namespace_t),         intent(in)  :: namespace
    type(space_t),             intent(in)  :: space

    PUSH_SUB(output_linear_medium_init)

    !%Variable LinearMediumOutput
    !%Type block
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what to print. The output files are written at the beginning of the run into the output directory for the
    !% linear medium.
    !% Each option must be in a separate row. Optionally individual output formats can be defined
    !% for each row (VTK format is supported) or they can be read separately from <tt>OutputFormat</tt> in the input file.
    !%
    !% Example:
    !% <br><br><tt>%LinearMediumOutput
    !% <br>&nbsp;&nbsp;permittivity
    !% <br>&nbsp;&nbsp;permeability
    !% <br>%<br></tt>
    !% This block supports all the formats of the <tt>Output</tt> block. See <tt>Output</tt>.
    !%Option points 1
    !% Outputs 1 if the a given point is inside the medium, and 0 otherwise. This can be used to check the grid points of the medium region.
    !%Option permittivity 2
    !% Output of the (static) space-dependent relative permittivity
    !%Option permeability 3
    !% Output of the (static) space-dependent relative permeability
    !%Option speed_of_light 3
    !% Output of the speed of light in atomic units
    !%End

    outp%what = .false.
    call io_function_read_what_how_when(namespace, space, outp%what, outp%how, outp%output_interval, &
      'LinearMediumOutput', 'OutputFormat', 'OutputInterval')


    !%Variable LinearMediumOutputDir
    !%Default "output_iter"
    !%Type string
    !%Section Output
    !%Description
    !% The name of the directory where <tt>Octopus</tt> stores the information
    !% about the linear medium system, as required by the <tt>LinearMediumOutput</tt> variable.
    !%End
    call parse_variable(namespace, 'LinearMediumOutputDir', "static", outp%iter_dir)
    if (any(outp%what) .and. maxval(outp%output_interval) > 0) then
      call io_mkdir(outp%iter_dir, namespace)
    end if
    call add_last_slash(outp%iter_dir)

    POP_SUB(output_linear_medium_init)
  end subroutine output_linear_medium_init


  ! ---------------------------------------------------------
  subroutine output_linear_medium(outp, namespace, space, mesh, dir, points_map, ep, mu, cc)
    type(output_t),                     intent(in)    :: outp
    type(namespace_t),                  intent(in)    :: namespace
    type(space_t),                      intent(in)    :: space
    class(mesh_t),                      intent(in)    :: mesh
    character(len=*),                   intent(in)    :: dir
    FLOAT, intent(in)            :: ep(:)
    FLOAT, intent(in)            :: mu(:)
    FLOAT, intent(in)            :: cc(:)
    integer, intent(in)          :: points_map(:)

    integer :: ierr
    FLOAT, allocatable :: dtmp(:), dtmp2(:)
    character(len=MAX_PATH_LEN) :: fname

    PUSH_SUB(output_linear_medium)

    if (any(outp%what)) then
      message(1) = "Info: Writing output to " // trim(dir)
      call messages_info(1, namespace=namespace)
      call io_mkdir(dir, namespace)
    endif

    ! Permittivity
    if (outp%what(OPTION__LINEARMEDIUMOUTPUT__PERMITTIVITY)) then
      write(fname, '(1a)') 'medium-permittivity'
      SAFE_ALLOCATE(dtmp(1:mesh%np))
      dtmp(:) = P_ep
      call get_medium_property(ep, points_map, dtmp) ! ep already has the P_ep factor
      dtmp(:) = dtmp / P_ep ! to print relative permittivity
      call dio_function_output(outp%how(OPTION__LINEARMEDIUMOUTPUT__PERMITTIVITY), dir, fname, namespace, space, &
        mesh, dtmp(:), unit_one, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! Permeability
    if (outp%what(OPTION__LINEARMEDIUMOUTPUT__PERMEABILITY)) then
      write(fname, '(1a)') 'medium-permeability'
      SAFE_ALLOCATE(dtmp(1:mesh%np))
      dtmp(:) = P_mu
      call get_medium_property(mu, points_map, dtmp) ! mu alredy has the P_mu factor
      dtmp(:) = dtmp / P_mu ! to print relative permability
      call dio_function_output(outp%how(OPTION__LINEARMEDIUMOUTPUT__PERMEABILITY), dir, fname, namespace, space, &
        mesh, dtmp(:), unit_one, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! Speed of light
    if (outp%what(OPTION__LINEARMEDIUMOUTPUT__SPEED_OF_LIGHT)) then
      write(fname, '(1a)') 'medium-speed-of-light'
      SAFE_ALLOCATE(dtmp(1:mesh%np))
      dtmp(:) = P_c
      call get_medium_property(cc, points_map, dtmp)
      call dio_function_output(outp%how(OPTION__LINEARMEDIUMOUTPUT__SPEED_OF_LIGHT), dir, fname, namespace, space, &
        mesh, dtmp(:), unit_one, ierr)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! Only points
    if (outp%what(OPTION__LINEARMEDIUMOUTPUT__POINTS)) then
      write(fname, '(1a)') 'medium-points'
      SAFE_ALLOCATE(dtmp(1:mesh%np))
      SAFE_ALLOCATE(dtmp2(1:mesh%np))
      dtmp(:) = M_ZERO
      dtmp2(:) = M_ONE
      call get_medium_property(dtmp2, points_map, dtmp) ! dtmp will have 1 only at the medium grid points
      call dio_function_output(outp%how(OPTION__LINEARMEDIUMOUTPUT__POINTS), dir, fname, namespace, space, &
        mesh, dtmp(:), unit_one, ierr)
      SAFE_DEALLOCATE_A(dtmp)
      SAFE_DEALLOCATE_A(dtmp2)
    end if


    POP_SUB(output_linear_medium)

  contains

    subroutine get_medium_property(medium_func, points_map, io_func)
      FLOAT,              intent(in)    :: medium_func(:)
      integer,            intent(in)    :: points_map(:)
      FLOAT,              intent(out)   :: io_func(:)

      integer :: ip, ip_in

      do ip_in = 1, size(points_map)
        ip          = points_map(ip_in)
        io_func(ip) = medium_func(ip_in)
      end do
    end subroutine get_medium_property

  end subroutine output_linear_medium

end module output_linear_medium_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
