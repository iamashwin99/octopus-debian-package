!! Copyright (C) 2002-2021 M. Marques, A. Castro, A. Rubio, G. Bertsch,
!!                         M. Lueders
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

!> This module implements a simple hash table for string valued keys
!! and integer values using the C++ STL unordered_map container.

module sihash_oct_m
  use global_oct_m
  use string_oct_m

  use iso_c_binding

  implicit none

  private
  public ::        &
    sihash_t,      &
    sihash_init,   &
    sihash_end,    &
    sihash_insert, &
    sihash_lookup, &
    sihash_iterator_t

  type sihash_t
    private

    type(c_ptr) :: map
  end type sihash_t

  type sihash_iterator_t
    private

    type(c_ptr)  :: iterator
    type(c_ptr)  :: end

  contains
    procedure :: start    => sihash_iterator_start
    procedure :: has_next => sihash_iterator_has_next
    procedure :: get_next => sihash_iterator_get_next
  end type sihash_iterator_t

contains

  ! ---------------------------------------------------------
  !> Initialize a hash table h with size entries. Since we use separate
  !! chaining, the number of entries in the hash table is, in
  !! principle, unlimited. We take the smallest prime number as table
  !! size that is greater or equal than the requested size to reduce
  !! collisions.
  subroutine sihash_init(h)
    type(sihash_t), intent(out) :: h

    interface
      subroutine sihash_map_init(map) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr) :: map
      end subroutine sihash_map_init
    end interface


    call sihash_map_init(h%map)

  end subroutine sihash_init


  ! ---------------------------------------------------------
  !> Free a hash table.
  subroutine sihash_end(h)
    type(sihash_t), intent(inout) :: h

    interface
      subroutine sihash_map_end(map) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr) :: map
      end subroutine sihash_map_end
    end interface

    call sihash_map_end(h%map)

  end subroutine sihash_end


  ! ---------------------------------------------------------
  !> Insert a (key, val) pair into the hash table h.
  subroutine sihash_insert(h, key, val)
    type(sihash_t),    intent(inout) :: h
    character(len=*),  intent(in)    :: key
    integer,           intent(in)    :: val

    interface
      subroutine sihash_map_insert(map, key, val) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr),            value :: map
        character(kind=c_char), intent(in)    :: key(*)
        integer(kind=c_int),    value   :: val
      end subroutine sihash_map_insert
    end interface

    call sihash_map_insert(h%map, string_f_to_c(key), val)
  end subroutine sihash_insert


  ! ---------------------------------------------------------
  !> Look up a value in the hash table h. If found is present, it
  !! indicates if key could be found in the table. If found = .false.,
  !! the return value of iihash_lookup is meaningless (and essentially
  !! undefined).
  integer function sihash_lookup(h, key, found)
    type(sihash_t),    intent(in)  :: h
    character(len=*),  intent(in)  :: key
    logical, optional, intent(out) :: found

    interface
      subroutine sihash_map_lookup(map, key, ifound, val) bind(c)
        use iso_c_binding
        implicit none

        type(c_ptr),            value         :: map
        character(kind=c_char), intent(in)    :: key(*)
        integer(kind=c_int),    intent(out)   :: ifound
        integer(kind=c_int),    intent(out)   :: val
      end subroutine sihash_map_lookup
    end interface

    integer :: ifound, val

    call sihash_map_lookup(h%map, string_f_to_c(key), ifound, val)

    found = (ifound == 1)

    sihash_lookup = -1
    if(found) sihash_lookup = val

  end function sihash_lookup

  ! ---------------------------------------------------------
  subroutine sihash_iterator_start(this, h)
    class(sihash_iterator_t),  intent(inout) :: this
    class(sihash_t),           intent(in)    :: h

    interface
      subroutine sihash_iterator_low_start(iterator, end, map) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr)            :: iterator
        type(c_ptr)            :: end
        type(c_ptr), value     :: map
      end subroutine sihash_iterator_low_start
    end interface

    call sihash_iterator_low_start(this%iterator, this%end, h%map)

  end subroutine sihash_iterator_start

  ! ---------------------------------------------------------
  logical function sihash_iterator_has_next(this)
    class(sihash_iterator_t), intent(in) :: this

    integer :: value

    interface
      subroutine sihash_iterator_low_has_next(iterator, end, value) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr),  value, intent(in)  :: iterator
        type(c_ptr),  value, intent(in)  :: end
        integer(kind=c_int), intent(out) :: value

      end subroutine sihash_iterator_low_has_next
    end interface

    call sihash_iterator_low_has_next(this%iterator, this%end, value)

    sihash_iterator_has_next = (value /= 0)

  end function sihash_iterator_has_next

  ! ---------------------------------------------------------
  integer function sihash_iterator_get_next(this) result(value)
    class(sihash_iterator_t), intent(inout) :: this

    interface
      subroutine sihash_iterator_low_get(iterator, value) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr),         intent(in)  :: iterator
        integer(kind=c_int), intent(out) :: value

      end subroutine sihash_iterator_low_get
    end interface

    call sihash_iterator_low_get(this%iterator, value)

  end function sihash_iterator_get_next


end module sihash_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
