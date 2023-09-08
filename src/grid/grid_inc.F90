!! Copyright (C) 2022 N. Tancogne-Dejean
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

!-------------------------------------------------------------------
subroutine X(grid_symmetrize_scalar_field)(gr, field, suppress_warning)
  type(grid_t),      intent(in)    :: gr
  R_TYPE,            intent(inout) :: field(:)
  logical, optional, intent(in)    :: suppress_warning

  R_TYPE, allocatable :: tmp(:)
  logical :: suppress_warning_

  PUSH_SUB(X(grid_symmetrize_scalar_field))

  ASSERT(ubound(field, dim=1) >= gr%np)

  suppress_warning_ = optional_default(suppress_warning, .false.)

  SAFE_ALLOCATE(tmp(1:gr%np))

  call X(symmetrizer_apply)(gr%symmetrizer, gr, field = field, &
    symmfield = tmp, suppress_warning = suppress_warning_)
  field(1:gr%np) = tmp(1:gr%np)

  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(X(grid_symmetrize_scalar_field))
end subroutine X(grid_symmetrize_scalar_field)

!-------------------------------------------------------------------
subroutine X(grid_symmetrize_vector_field)(gr, field, suppress_warning)
  type(grid_t),      intent(in)    :: gr
  R_TYPE,            intent(inout) :: field(:,:)
  logical, optional, intent(in)    :: suppress_warning

  R_TYPE, allocatable :: tmp(:,:)
  logical :: suppress_warning_

  PUSH_SUB(X(grid_symmetrize_vector_field))

  ASSERT(ubound(field, dim=1)>= gr%np)
  ASSERT(ubound(field, dim=2)>= gr%der%dim)

  suppress_warning_ = optional_default(suppress_warning, .false.)

  SAFE_ALLOCATE(tmp(1:gr%np, 1:gr%der%dim))

  call X(symmetrizer_apply)(gr%symmetrizer, gr, field_vector = field, &
    symmfield_vector = tmp, suppress_warning = suppress_warning_)
  field(1:gr%np, 1:gr%der%dim) = tmp(1:gr%np, 1:gr%der%dim)

  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(X(grid_symmetrize_vector_field))
end subroutine X(grid_symmetrize_vector_field)


!-------------------------------------------------------------------
subroutine X(grid_symmetrize_single)(gr, iop, field, symm_field, suppress_warning)
  type(grid_t),      intent(in)    :: gr
  integer,           intent(in)    :: iop
  R_TYPE,            intent(in)    :: field(:)
  R_TYPE,            intent(out)   :: symm_field(:)
  logical, optional, intent(in)    :: suppress_warning

  logical :: suppress_warning_

  PUSH_SUB(X(grid_symmetrize_single))

  ASSERT(ubound(field, dim=1) >= gr%np)
  ASSERT(ubound(symm_field, dim=1) >= gr%np)

  suppress_warning_ = optional_default(suppress_warning, .false.)

  call X(symmetrizer_apply_single)(gr%symmetrizer, gr, iop, field, symm_field)

  POP_SUB(X(grid_symmetrize_single))
end subroutine X(grid_symmetrize_single)
