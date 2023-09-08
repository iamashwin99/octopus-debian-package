!! Copyright (C) 2007 X. Andrade
!! Copyright (C) 2021 N. Tancogne-Dejean
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

! --------------------------------------------------------------------------
!> Returns f_out = H' f_in, where H' is perturbation Hamiltonian
!! Note that e^ikr phase is applied to f_in, then is removed afterward
subroutine X(perturbation_electric_apply)(this, namespace, space, gr, hm, ik, f_in, f_out, set_bc)
  class(perturbation_electric_t), intent(in)    :: this
  type(namespace_t),              intent(in)    :: namespace
  type(space_t),                  intent(in)    :: space
  type(grid_t),                   intent(in)    :: gr
  type(hamiltonian_elec_t),       intent(in)    :: hm
  integer,                        intent(in)    :: ik
  R_TYPE,                         intent(in)    :: f_in(:, :)
  R_TYPE,                         intent(out)   :: f_out(:, :)
  logical,              optional, intent(in)    :: set_bc

  integer :: ip, idim
  type(profile_t), save :: prof

  PUSH_SUB(X(perturbation_electric_apply))
  call profiling_in(prof, TOSTRING(X(PERT_ELEC_APPLY)))

  ASSERT(this%dir /= -1)

  ! electric perturbation does not need phases since (e^-ikr)r(e^ikr) = r

  do idim = 1, hm%d%dim
    do ip = 1, gr%np
      f_out(ip, idim) = f_in(ip, idim) * gr%x(ip, this%dir)
    end do
  end do

  call profiling_out(prof)
  POP_SUB(X(perturbation_electric_apply))
end subroutine X(perturbation_electric_apply)

! --------------------------------------------------------------------------
subroutine X(perturbation_electric_apply_order_2) (this, namespace, space, gr, hm, ik, f_in, f_out)
  class(perturbation_electric_t), intent(in)    :: this
  type(namespace_t),              intent(in)    :: namespace
  type(space_t),                  intent(in)    :: space
  type(grid_t),                   intent(in)    :: gr
  type(hamiltonian_elec_t),       intent(in)    :: hm
  integer,                        intent(in)    :: ik
  R_TYPE,                         intent(in)    :: f_in(:, :)
  R_TYPE,                         intent(out)   :: f_out(:, :)

  PUSH_SUB(X(perturbation_electric_apply_order_2))

  f_out(1:gr%np, 1:hm%d%dim) = R_TOTYPE(M_ZERO)

  POP_SUB(X(perturbation_electric_apply_order_2))
end subroutine X(perturbation_electric_apply_order_2)
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
