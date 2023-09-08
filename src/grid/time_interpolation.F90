!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2023 S. Ohlmann
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

module time_interpolation_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m

  implicit none

  private
  public :: &
    time_interpolation_t

  type time_interpolation_t
    private
    CMPLX, allocatable :: zfield(:, :, :) !< complex fields for interpolation
    FLOAT, allocatable :: dfield(:, :, :) !< real fields for interpolation
    FLOAT, allocatable :: times(:)        !< corresponding times
    integer :: max_depth                  !< maximum interpolation depth
    integer :: depth                      !< current interpolation depth (allow for smaller depth in beginning)
    integer :: np                         !< first array dimension of fields
    integer :: dim                        !< second array dimension of fields
    logical :: cmplx                      !< do we have real or complex values
    character(len=MAX_PATH_LEN) :: label  !< a label needed for restart
  contains
    procedure :: dtime_interpolation_add_time, ztime_interpolation_add_time
    generic :: add_time => dtime_interpolation_add_time, ztime_interpolation_add_time
    procedure :: dtime_interpolation_interpolate, ztime_interpolation_interpolate
    generic :: interpolate => dtime_interpolation_interpolate, ztime_interpolation_interpolate
    procedure :: read_restart => time_interpolation_read_restart
    procedure :: write_restart => time_interpolation_write_restart
    final :: time_interpolation_finalize
  end type time_interpolation_t

  interface time_interpolation_t
    procedure time_interpolation_constructor
  end interface time_interpolation_t

contains
  function time_interpolation_constructor(np, dim, depth, cmplx, label) result(this)
    integer,                     intent(in) :: np
    integer,                     intent(in) :: dim
    integer,                     intent(in) :: depth
    logical,                     intent(in) :: cmplx
    character(len=*),            intent(in) :: label
    class(time_interpolation_t), pointer :: this

    PUSH_SUB(time_interpolation_constructor)

    SAFE_ALLOCATE(this)

    this%np = np
    this%dim = dim
    this%max_depth = depth
    this%depth = 0
    this%cmplx = cmplx
    this%label = trim(label)

    if (this%cmplx) then
      SAFE_ALLOCATE(this%zfield(1:np, 1:dim, 1:depth))
    else
      SAFE_ALLOCATE(this%dfield(1:np, 1:dim, 1:depth))
    end if
    SAFE_ALLOCATE(this%times(1:depth))

    POP_SUB(time_interpolation_constructor)
  end function time_interpolation_constructor

  subroutine time_interpolation_finalize(this)
    type(time_interpolation_t), intent(inout) :: this

    PUSH_SUB(time_interpolation_finalize)

    if (this%cmplx) then
      SAFE_DEALLOCATE_A(this%zfield)
    else
      SAFE_DEALLOCATE_A(this%dfield)
    end if
    SAFE_DEALLOCATE_A(this%times)

    POP_SUB(time_interpolation_finalize)
  end subroutine time_interpolation_finalize

  subroutine time_interpolation_write_restart(this, mesh, space, restart, err)
    class(time_interpolation_t), intent(in)    :: this
    class(mesh_t),                intent(in)    :: mesh
    type(space_t),                intent(in)    :: space
    type(restart_t),              intent(in)    :: restart
    integer,                      intent(out)   :: err

    integer :: itime, idim, err_restart, iunit
    character(len=MAX_PATH_LEN) :: filename, lines(1)

    PUSH_SUB(time_interpolation_write_restart)

    ! write fields
    err = 0
    do itime = 1, this%depth
      do idim = 1, this%dim
        write(filename, '(a1,i2.2,a1,i3.3)') '_', itime, '_',  idim
        filename = "field_" // trim(this%label) // trim(filename)
        if (this%cmplx) then
          call zrestart_write_mesh_function(restart, space, filename, mesh, &
            this%zfield(1:this%np, idim, itime), err_restart)
        else
          call drestart_write_mesh_function(restart, space, filename, mesh, &
            this%dfield(1:this%np, idim, itime), err_restart)
        end if
        if (err_restart /= 0) err = err + 1
      end do
    end do

    ! write times
    call drestart_write_binary(restart, "field_times_"//trim(this%label), this%depth, this%times(1:this%depth), err_restart)
    if (err_restart /= 0) err = err + 1

    ! write depth
    iunit = restart_open(restart, "field_"//trim(this%label))
    write(lines(1), '(i2.2)') this%depth
    call restart_write(restart, iunit, lines, 1, err_restart)
    if (err_restart /= 0) err = err + 1
    call restart_close(restart, iunit)

    POP_SUB(time_interpolation_write_restart)
  end subroutine time_interpolation_write_restart

  subroutine time_interpolation_read_restart(this, mesh, space, restart, err)
    class(time_interpolation_t),  intent(inout) :: this
    class(mesh_t),                intent(in)    :: mesh
    type(space_t),                intent(in)    :: space
    type(restart_t),              intent(in)    :: restart
    integer,                      intent(out)   :: err

    integer :: itime, idim, err_restart, iunit
    character(len=MAX_PATH_LEN) :: filename, lines(1)

    PUSH_SUB(time_interpolation_read_restart)

    ! read depth
    iunit = restart_open(restart, "field_"//trim(this%label))
    call restart_read(restart, iunit, lines, 1, err_restart)
    if (err_restart /= 0) then
      err = err + 1
      POP_SUB(time_interpolation_read_restart)
      return
    else
      read(lines(1), '(i2.2)') this%depth
    end if
    call restart_close(restart, iunit)

    ! read fields
    err = 0
    do itime = 1, this%depth
      do idim = 1, this%dim
        write(filename, '(a1,i2.2,a1,i3.3)') '_', itime, '_',  idim
        filename = "field_" // trim(this%label) // trim(filename)
        if (this%cmplx) then
          call zrestart_read_mesh_function(restart, space, filename, mesh, &
            this%zfield(1:this%np, idim, itime), err_restart)
        else
          call drestart_read_mesh_function(restart, space, filename, mesh, &
            this%dfield(1:this%np, idim, itime), err_restart)
        end if
        if (err_restart /= 0) err = err + 1
      end do
    end do

    ! read times
    call drestart_read_binary(restart, "field_times_"//trim(this%label), this%depth, this%times(1:this%depth), err_restart)
    if (err_restart /= 0) err = err + 1

    POP_SUB(time_interpolation_read_restart)
  end subroutine time_interpolation_read_restart

#include "real.F90"
#include "time_interpolation_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "time_interpolation_inc.F90"
#include "undef.F90"

end module time_interpolation_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
