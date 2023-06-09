!! Copyright (C) 2010 X. Andrade
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

module types_oct_m

  implicit none

  private

  public ::            &
    type_t,            &
    types_get_size,    &
    operator(==),      &
    operator(/=),      &
    type_is_complex

  type type_t
    private
    integer :: itype
  end type type_t

  type(type_t), public :: TYPE_NONE         = type_t(0)
  type(type_t), public :: TYPE_FLOAT        = type_t(1)
  type(type_t), public :: TYPE_CMPLX        = type_t(2)
  type(type_t), public :: TYPE_INTEGER      = type_t(3)
  type(type_t), public :: TYPE_BYTE         = type_t(4)
  type(type_t), public :: TYPE_INTEGER8     = type_t(5)

  interface operator(==)
    module procedure types_equal
  end interface operator(==)

  interface operator(/=)
    module procedure types_not_equal
  end interface operator(/=)

  integer :: sizes(5) = (/8, 16, 4, 1, 8/)

contains

  integer pure function types_get_size(this) result(size)
    type(type_t), intent(in) :: this

    size = sizes(this%itype)
  end function types_get_size

  ! -----------------------------------------------------

  logical pure function types_equal(ta, tb) result(equal)
    type(type_t), intent(in) :: ta
    type(type_t), intent(in) :: tb

    equal = ta%itype == tb%itype

  end function types_equal

  ! -----------------------------------------------------

  logical pure function types_not_equal(ta, tb) result(equal)
    type(type_t), intent(in) :: ta
    type(type_t), intent(in) :: tb

    equal = ta%itype /= tb%itype

  end function types_not_equal

  ! -----------------------------------------------------

  logical pure function type_is_complex(this) result(is_complex)
    type(type_t), intent(in) :: this

    is_complex = this == TYPE_CMPLX

  end function type_is_complex

end module types_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
