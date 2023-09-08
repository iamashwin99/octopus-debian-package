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

module sphash_oct_m
  use global_oct_m
  use string_oct_m

  use iso_c_binding

  implicit none

  private
  public ::        &
    sphash_t,      &
    sphash_init,   &
    sphash_end,    &
    sphash_insert, &
    sphash_lookup, &
    sphash_iterator_t

  ! ---------------------------------------------------------

  type sphash_value_t
    private
    class(*), pointer :: value
    logical           :: clone
  contains
    procedure :: get => sphash_value_get
    final :: sphash_value_finalize
  end type sphash_value_t

  interface sphash_value_t
    procedure :: sphash_value_constructor
  end interface sphash_value_t

  ! ---------------------------------------------------------

  type sphash_t
    private

    type(c_ptr) :: map
  end type sphash_t

  ! ---------------------------------------------------------

  type sphash_iterator_t
    private

    type(c_ptr)  :: iterator
    type(c_ptr)  :: end

  contains
    procedure :: start    => sphash_iterator_start
    procedure :: has_next => sphash_iterator_has_next
    procedure :: get_next => sphash_iterator_get_next
  end type sphash_iterator_t

contains


  function sphash_value_constructor(value, clone) result(constructor)

    class(*),              target :: value
    logical,             optional :: clone
    type(sphash_value_t), pointer :: constructor

    allocate(constructor)
    constructor%clone = optional_default(clone, .false.)

    if (constructor%clone) then
      allocate(constructor%value, source=value)
    else
      constructor%value => value
    end if

  end function sphash_value_constructor


  subroutine sphash_value_finalize(this)

    type(sphash_value_t) :: this

    if (associated(this%value)) then
      if (this%clone) then
        deallocate(this%value)
      else
        nullify(this%value)
      end if
    end if

  end subroutine sphash_value_finalize


  function sphash_value_get(this) result(value)
    class(sphash_value_t), intent(in)  :: this
    class(*),              pointer     :: value

    value => this%value

  end function sphash_value_get

  ! ---------------------------------------------------------
  !> Initialize a hash table h with size entries. Since we use separate
  !! chaining, the number of entries in the hash table is, in
  !! principle, unlimited. We take the smallest prime number as table
  !! size that is greater or equal than the requested size to reduce
  !! collisions.
  subroutine sphash_init(h)
    type(sphash_t), intent(out) :: h

    interface
      subroutine sphash_map_init(map) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr) :: map
      end subroutine sphash_map_init
    end interface


    call sphash_map_init(h%map)

  end subroutine sphash_init


  ! ---------------------------------------------------------
  !> Free a hash table.
  subroutine sphash_end(h)
    type(sphash_t), intent(inout) :: h

    type(sphash_iterator_t) :: it
    type(c_ptr) :: tmp_ptr
    type(sphash_value_t), pointer :: tmp_value

    interface
      subroutine sphash_iterator_low_get(iterator, value_ptr) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr),    intent(in)  :: iterator
        type(c_ptr),    intent(out) :: value_ptr

      end subroutine sphash_iterator_low_get

      subroutine sphash_map_end(map) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr) :: map
      end subroutine sphash_map_end
    end interface

    call it%start(h)

    do while (it%has_next())
      call sphash_iterator_low_get(it%iterator, tmp_ptr)
      call c_f_pointer(tmp_ptr, tmp_value)
      deallocate(tmp_value)
    end do

    call sphash_map_end(h%map)

  end subroutine sphash_end


  ! ---------------------------------------------------------
  !> Insert a (key, val) pair into the hash table h.
  !> If clone=.true., the object will be copied.
  subroutine sphash_insert(h, key, val, clone)
    type(sphash_t),    intent(inout) :: h
    character(len=*),  intent(in)    :: key
    class(*),  target, intent(in)    :: val
    logical, optional, intent(in)    :: clone

    type(sphash_value_t), pointer    :: value

    interface
      subroutine sphash_map_insert(map, key, val) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr),            value      :: map
        character(kind=c_char), intent(in) :: key(*)
        type(c_ptr),            value      :: val
      end subroutine sphash_map_insert
    end interface

    value => sphash_value_t(val, clone)

    call sphash_map_insert(h%map, string_f_to_c(key), c_loc(value))

  end subroutine sphash_insert


  ! ---------------------------------------------------------
  !> Look up a value in the hash table h. If found is present, it
  !! indicates if key could be found in the table. If found = .false.,
  !! the return value of iihash_lookup is meaningless (and essentially
  !! undefined).
  function sphash_lookup(h, key, found) result(value)
    type(sphash_t),    intent(in)  :: h
    character(len=*),  intent(in)  :: key
    logical, optional, intent(out) :: found
    class(*),  pointer             :: value

    type(sphash_value_t), pointer :: tmp_value
    character(kind=c_char), dimension(c_str_len(key)) :: c_key

    interface
      subroutine sphash_map_lookup(map, key, ifound, val) bind(c)
        use iso_c_binding
        implicit none

        type(c_ptr),                  value         :: map
        character(kind=c_char),       intent(in)    :: key(*)
        integer(kind=c_int),          intent(out)   :: ifound
        type(c_ptr),                  intent(out)   :: val
      end subroutine sphash_map_lookup
    end interface

    integer :: ifound
    type(c_ptr) :: val

    c_key = string_f_to_c(key)
    !Passing string_f_to_c as an arg creates a temporary array (ifort)
    call sphash_map_lookup(h%map, c_key, ifound, val)

    found = (ifound == 1)

    nullify(value)
    if (found) then
      call c_f_pointer(val, tmp_value)
      value => tmp_value%get()
    end if

  end function sphash_lookup

  ! ---------------------------------------------------------
  subroutine sphash_iterator_start(this, h)
    class(sphash_iterator_t),  intent(inout) :: this
    class(sphash_t),           intent(in)    :: h

    interface
      subroutine sphash_iterator_low_start(iterator, end, map) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr)            :: iterator
        type(c_ptr)            :: end
        type(c_ptr), value     :: map
      end subroutine sphash_iterator_low_start
    end interface

    call sphash_iterator_low_start(this%iterator, this%end, h%map)

  end subroutine sphash_iterator_start

  ! ---------------------------------------------------------
  logical function sphash_iterator_has_next(this)
    class(sphash_iterator_t), intent(in) :: this

    integer :: value

    interface
      subroutine sphash_iterator_low_has_next(iterator, end, value) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr),  value, intent(in)  :: iterator
        type(c_ptr),  value, intent(in)  :: end
        integer(kind=c_int), intent(out) :: value

      end subroutine sphash_iterator_low_has_next
    end interface

    call sphash_iterator_low_has_next(this%iterator, this%end, value)

    sphash_iterator_has_next = (value /= 0)

  end function sphash_iterator_has_next

  ! ---------------------------------------------------------
  function sphash_iterator_get_next(this) result(value)
    class(sphash_iterator_t), intent(inout) :: this
    class(*),                 pointer       :: value

    type(c_ptr) :: tmp_ptr
    type(sphash_value_t), pointer :: tmp_value

    interface
      subroutine sphash_iterator_low_get(iterator, value_ptr) bind(c)
        use iso_c_binding
        import
        implicit none

        type(c_ptr),    intent(in)  :: iterator
        type(c_ptr),    intent(out) :: value_ptr

      end subroutine sphash_iterator_low_get
    end interface

    call sphash_iterator_low_get(this%iterator, tmp_ptr)

    call c_f_pointer(tmp_ptr, tmp_value)
    value => tmp_value%get()

  end function sphash_iterator_get_next


end module sphash_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
