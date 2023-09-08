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

subroutine X(time_interpolation_add_time)(this, time, field)
  class(time_interpolation_t), intent(inout) :: this
  FLOAT,                       intent(in)    :: time
  R_TYPE,                      intent(in)    :: field(:, :)

  integer :: itime

  PUSH_SUB(X(time_interpolation_add_time))

#ifdef R_TCOMPLEX
  ASSERT(this%cmplx)
#else
  ASSERT(.not.this%cmplx)
#endif

  ! only add the field if the time is not already the last time added
  if (.not.(this%depth > 0 .and. abs(time - this%times(1)) < M_TINY)) then
    if (this%depth < this%max_depth) then
      this%depth = this%depth + 1
    end if
    ! shift old entries by one
    do itime = this%depth, 2, -1
      this%times(itime) = this%times(itime-1)
      call lalg_copy(this%np, this%dim, this%X(field)(:, :, itime-1), this%X(field)(:, :, itime))
    end do
    ! set most recent entry
    this%times(1) = time
    call lalg_copy(this%np, this%dim, field(:, :), this%X(field)(:, :, 1))
  end if

  POP_SUB(X(time_interpolation_add_time))
end subroutine X(time_interpolation_add_time)

subroutine X(time_interpolation_interpolate)(this, time, field)
  class(time_interpolation_t), intent(inout) :: this
  FLOAT,                       intent(in)    :: time
  R_TYPE,                      intent(inout) :: field(:, :)

  PUSH_SUB(X(time_interpolation_interpolate))

#ifdef R_TCOMPLEX
  ASSERT(this%cmplx)
#else
  ASSERT(.not.this%cmplx)
#endif

  ! at least one timestep needs to have been added
  ASSERT(this%depth > 0)

  if (abs(time - this%times(1)) < M_TINY) then
    ! special case: exact time requested, does not need to be interpolated
    call lalg_copy(this%np, this%dim, this%X(field)(:, :, 1), field(:, :))
  else
    call interpolate(this%times(1:this%depth), this%X(field)(:, :, 1:this%depth), time, field(:, :))
  end if

  POP_SUB(X(time_interpolation_interpolate))
end subroutine X(time_interpolation_interpolate)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
