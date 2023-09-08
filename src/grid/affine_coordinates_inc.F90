!! Copyright (C) 2021 M. Oliveira
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

! -----------------------------------------------------------------------------
pure subroutine X(affine_coordinates_vector_from_cartesian)(this, xx, vv, src)
  class(affine_coordinates_t), intent(in)    :: this
  FLOAT,                       intent(in)    :: xx(:)
  R_TYPE,                      intent(inout) :: vv(:)
  R_TYPE, optional,            intent(in)    :: src(:)

  ! no PUSH_SUB, called too often

  ! Contravariant vectors transform as V_uvw = B V_xyz
  if (present(src)) then
    vv(1:this%dim) = matmul(this%basis%change_of_basis_matrix(1:this%dim, 1:this%dim), src(1:this%dim))
  else
    vv(1:this%dim) = matmul(this%basis%change_of_basis_matrix(1:this%dim, 1:this%dim), vv(1:this%dim))
  end if

end subroutine X(affine_coordinates_vector_from_cartesian)

! -----------------------------------------------------------------------------
pure subroutine X(affine_coordinates_covector_to_cartesian)(this, xx, cv, src)
  class(affine_coordinates_t), intent(in)    :: this
  FLOAT,                       intent(in)    :: xx(:)
  R_TYPE,                      intent(inout) :: cv(:)
  R_TYPE, optional,            intent(in)    :: src(:)

  ! no PUSH_SUB, called too often

  ! Covectors, like the gradient, transform as V_xyw = Bt V_uvw (see Chelikowsky after Eq. 10).
  if (present(src)) then
    cv(1:this%dim) = matmul(src(1:this%dim), this%basis%change_of_basis_matrix(1:this%dim, 1:this%dim))
  else
    cv(1:this%dim) = matmul(cv(1:this%dim), this%basis%change_of_basis_matrix(1:this%dim, 1:this%dim))
  end if

end subroutine X(affine_coordinates_covector_to_cartesian)
