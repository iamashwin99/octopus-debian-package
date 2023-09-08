!! Copyright (C) 2022  N. Tancogne-Dejean
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

module ext_partner_list_oct_m
  use debug_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use interaction_partner_oct_m
  use lasers_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                          &
    list_has_lasers,                      &
    list_has_gauge_field,                 &
    list_get_lasers,                      &
    list_get_gauge_field

contains

  ! ---------------------------------------------------------
  logical function list_has_lasers(partners)
    type(partner_list_t), intent(in)  :: partners

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner

    PUSH_SUB(list_has_lasers)

    list_has_lasers = .false.
    call iter%start(partners)
    do while (iter%has_next() .and. .not. list_has_lasers)
      partner => iter%get_next()
      select type(partner)
      type is(lasers_t)
        list_has_lasers = .true.
      end select
    end do

    POP_SUB(list_has_lasers)
  end function list_has_lasers

  ! ---------------------------------------------------------
  logical function list_has_gauge_field(partners)
    type(partner_list_t), intent(in)  :: partners

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner

    PUSH_SUB(list_has_gauge_field)

    list_has_gauge_field = .false.
    call iter%start(partners)
    do while (iter%has_next() .and. .not. list_has_gauge_field)
      partner => iter%get_next()
      select type(partner)
      type is(gauge_field_t)
        list_has_gauge_field = .true.
      end select
    end do

    POP_SUB(list_has_gauge_field)
  end function list_has_gauge_field

  ! ---------------------------------------------------------
  function list_get_gauge_field(partners) result(value)
    type(partner_list_t), intent(in)  :: partners
    type(gauge_field_t),  pointer     :: value

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner

    PUSH_SUB(list_get_gauge_field)

    value => null()
    call iter%start(partners)
    do while (iter%has_next() .and. .not. associated(value))
      partner => iter%get_next()
      select type(partner)
      type is(gauge_field_t)
        value => partner
      end select
    end do

    POP_SUB(list_get_gauge_field)
  end function list_get_gauge_field


  ! ---------------------------------------------------------
  function list_get_lasers(partners) result(value)
    type(partner_list_t), intent(in)  :: partners
    type(lasers_t),       pointer     :: value

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner

    PUSH_SUB(list_get_lasers)

    value => null()
    call iter%start(partners)
    do while (iter%has_next() .and. .not. associated(value))
      partner => iter%get_next()
      select type(partner)
      type is(lasers_t)
        value => partner
      end select
    end do

    POP_SUB(list_get_lasers)
  end function list_get_lasers

end module ext_partner_list_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
