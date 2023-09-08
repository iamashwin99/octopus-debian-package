!! Copyright (C) 2019 M. Lueders
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

!> \ingroup Fortran_Module

module basis_set_abst_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public :: basis_set_abst_t

  !> @brief abstract class for basis sets
  !!
  !! Abstract class for basis sets. In our context, basis sets include meshes, atomic orbitals, plane waves,
  !! possibly wavelets, etc. For this reason, the abstract class is very minimalistic, and only contains
  !! essential quantities and functions
  type, abstract :: basis_set_abst_t
    private
    logical :: time_dependent   !< flag for time-dependent basis sets
  contains
    procedure(init),       deferred        :: init
    procedure(end),        deferred        :: end
    procedure(write_info), deferred        :: write_info
    procedure,             non_overridable :: is_time_dependent
    procedure,             non_overridable :: set_time_dependent
  end type basis_set_abst_t

  abstract interface
    subroutine init(this)
      import basis_set_abst_t
      class(basis_set_abst_t), intent(inout) :: this
    end subroutine init

    subroutine end(this)
      import basis_set_abst_t
      class(basis_set_abst_t), intent(inout) :: this
    end subroutine end

    subroutine write_info(this, iunit, namespace)
      import basis_set_abst_t
      import namespace_t
      class(basis_set_abst_t),     intent(in) :: this
      integer,           optional, intent(in) :: iunit
      type(namespace_t), optional, intent(in) :: namespace
    end subroutine write_info
  end interface

contains

  function is_time_dependent(this) result(td_flag)
    class(basis_set_abst_t), intent(in) :: this
    logical :: td_flag

    PUSH_SUB(is_time_dependent)

    td_flag = this%time_dependent

    POP_SUB(is_time_dependent)
  end function is_time_dependent

  subroutine set_time_dependent(this, td_flag)
    class(basis_set_abst_t), intent(inout) :: this
    logical,                 intent(in)    :: td_flag

    PUSH_SUB(set_time_dependent)

    this%time_dependent = td_flag

    POP_SUB(set_time_dependent)
  end subroutine set_time_dependent

end module basis_set_abst_oct_m
