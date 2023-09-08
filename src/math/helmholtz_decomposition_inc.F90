!! Copyright (C) 2019 R. Jestaedt, H. Appel, F. Bonafe, M. Oliveira, N. Tancogne-Dejean
!! Copyright (C) 2022-2023 F. Troisi
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

subroutine X(apply_inner_stencil_mask)(this, field)
  class(helmholtz_decomposition_t), intent(inout) :: this
  R_TYPE,                           intent(inout) :: field(:,:)

  integer :: ip, ii

  PUSH_SUB(X(apply_inner_stencil_mask))

  !$omp parallel do private(ii)
  do ip = 1, SIZE(this%inner_stencil, dim = 1)
    ii = this%inner_stencil(ip)
    field(ii, :) = M_ZERO
  end do

  POP_SUB(X(apply_inner_stencil_mask))
end subroutine X(apply_inner_stencil_mask)

subroutine X(get_vector_potential)(this, namespace, vector_potential, total_field, apply_boundary)
  class(helmholtz_decomposition_t), intent(inout) :: this                  !< Helmholtz object
  type(namespace_t),               intent(in)     :: namespace
  R_TYPE,                          intent(out)    :: vector_potential(:,:) !< Vector potential obtained from the total_field.
  !!                                                                          The field will be defined on the grid of the
  !!                                                                          Helmholtz type (this%grid), UP UNTIL THIS%GRID%NP
  !!                                                                          vector_potential(1:this%grid%np_part,
  !!                                                                                           1:this%grid%box%dim)
  R_TYPE,                          intent(inout)  :: total_field(:,:)      !< Total field to decompose. To ensure the correct
  !!                                                                          behaviour this field must be correctly defined on the
  !!                                                                          grid of the system that calls Helmholtz (sys_gr),
  !!                                                                          UP UNTIL SYS_GR%NP_PART
  !!                                                                          total_field(1:sys_gr%np_part, 1:sys_gr%box%dim)
  logical, optional,               intent(in)     :: apply_boundary        !< Should the curl apply boundary conditions?

  integer             :: idim
  R_TYPE, allocatable :: support_field(:,:)
  logical             :: apply_boundary_

  PUSH_SUB(X(get_vector_potential))

  SAFE_ALLOCATE(support_field(1:this%sys_grid%np_part, 1:this%sys_grid%box%dim))
  support_field = M_ZERO
  ! if not specified, we apply the boundary condition when computing the curl
  apply_boundary_ = optional_default(apply_boundary, .true.)

  ! First of all, we have to compute the curl of the field. We do this on the system box
  call X(derivatives_curl)(this%sys_grid%der, total_field, support_field, set_bc = apply_boundary_)

  ! Then we solve poisson equation to compute helmholtz decomposition integral.
  ! The curl is a vector, so we apply poisson in all directions
  vector_potential = M_ZERO
  do idim = 1, this%sys_grid%box%dim
    call X(poisson_solve)(this%poisson_solver, namespace, vector_potential(:, idim), support_field(:, idim))
  end do

  if (this%compute_surface_correction) then
    support_field = M_ZERO
    ! Compute the surface correction
    call X(compute_surface_correction_trans_field)(this, namespace, total_field, support_field)
    ! Add the surface correction to the vector potential
    vector_potential(:, :) = vector_potential(:, :) - support_field(:, :)
  end if
  ! Add the prefactor 1/(4*pi)
  vector_potential = R_TOTYPE(M_ONE/(M_FOUR * M_PI)) * vector_potential

  SAFE_DEALLOCATE_A(support_field)

  POP_SUB(X(get_vector_potential))
end subroutine X(get_vector_potential)

subroutine X(get_trans_field)(this, namespace, transverse_field, total_field, vector_potential, apply_boundary)
  class(helmholtz_decomposition_t), intent(inout) :: this                  !< Helmholtz object
  type(namespace_t),                intent(in)    :: namespace
  R_TYPE,                           intent(out)   :: transverse_field(:,:) !< Transverse field obtained from the total_field.
  !!                                                                          It will be defined on the inner_grid of the Helmholtz
  !!                                                                          type (this%inner_grid), UP UNTIL INNER_GR%NP
  !!                                                                          transverse_field(1:this%inner_grid%np_part,
  !!                                                                                          1:this%inner_grid%box%dim)
  R_TYPE, optional,                 intent(inout) :: total_field(:,:)      !< Total field to decompose. To ensure the correct
  !!                                                                          behaviour this field must be correctly defined on the
  !!                                                                          grid of the system that calls Helmholtz (sys_gr),
  !!                                                                          UP UNTIL SYS_GR%NP_PART
  !!                                                                          total_field(1:sys_gr%np_part, 1:sys_gr%box%dim)
  R_TYPE, optional,                 intent(inout) :: vector_potential(:,:) !< Vector potential from which trans field whould be
  !!                                                                          computed
  logical, optional,                intent(in)    :: apply_boundary        !< Should the curl apply boundary conditions?

  R_TYPE, allocatable :: vec_pot(:,:)

  PUSH_SUB(X(get_trans_field))

  if (present(vector_potential)) then
    ! Compute transverse field from vector potential
    call X(derivatives_curl)(this%sys_grid%der, vector_potential, transverse_field, set_bc = .true.)
  else if (present(total_field)) then
    ! Compute transverse field from total field
    SAFE_ALLOCATE(vec_pot(1:this%sys_grid%np_part, 1:this%sys_grid%box%dim))
    ! First, get the vector potential
    call this%get_vector_potential(namespace, vec_pot, total_field, apply_boundary)

    ! Then we compute the curl again to retrieve the divergence-free term in the helmholtz decomposition
    call X(derivatives_curl)(this%sys_grid%der, vec_pot, transverse_field, set_bc = .true.)
    SAFE_DEALLOCATE_A(vec_pot)
  else
    message(1) = "Helmholtz decomposition error - you must pass either the total_field or the vector_potential"
    call messages_fatal(1)
  end if

  ! Finally, apply the mask to take care of the spikes generated by the curl
  call this%apply_inner_stencil_mask(transverse_field)

  POP_SUB(X(get_trans_field))
end subroutine X(get_trans_field)

! TODO: Issue 705 (ftroisi) - Check that inclusion of reduction over mesh does not affect the results. Add openmp support
subroutine X(compute_surface_correction_trans_field)(this, namespace, field, surface_correction)
  class(helmholtz_decomposition_t), intent(in)  :: this                     !< Helmholtz object
  type(namespace_t),                intent(in)  :: namespace
  R_TYPE,                           intent(in)  :: field(:,:)
  R_TYPE,                           intent(out) :: surface_correction(:,:)

  integer :: sp, vp, ii
  FLOAT   :: rr, surface_element
  FLOAT   :: r_vp(this%sys_grid%box%dim), r_sp(this%sys_grid%box%dim), normal_vector(this%sys_grid%box%dim)
  R_TYPE  :: integrand(this%sys_grid%box%dim)

  PUSH_SUB(X(compute_surface_correction_trans_field))
  ASSERT(.not. this%sys_grid%use_curvilinear)
  surface_correction = M_ZERO

  ! First of all cycle over all surface points
  do ii = 1, SIZE(this%surface_points, dim = 1)
    ! Check that the point belongs to the surface and that is not the same of the volume point (avoid 0 at denominator)
    sp = this%surface_points(ii)
    ! Get coordinates of surface point
    call mesh_r(this%sys_grid, sp, rr, coords = r_sp)
    ! Then get the normal vector in that point and the Jacobian to account for the fact that the spacing between
    ! surface points is different than the spacing on the normal grid. Since surface points are always the same,
    ! it is enough to compute the normal vector and the Jacobian once
    call this%sys_grid%box%get_surface_point_info(r_sp, this%sys_grid%spacing(:), normal_vector(:), surface_element)

    ! Compute the cross prod between the normal vec and the field, times the surface_element (does not depend on the volume point)
    integrand(:) = X(cross_product)(R_TOTYPE(normal_vector(:)), field(sp, :)) * surface_element

    ! Then cycle over the volume points of the submesh
    do vp = 1, this%sys_grid%np
      ! If vp == sp we are dealing with the singullar point, which must be treated differently
      if (vp /= sp) then
        ! Get coordinates of volume point
        call mesh_r(this%sys_grid, vp, rr, coords = r_vp)
        ! Perform the integral
        surface_correction(vp, :) = surface_correction(vp, :) + integrand(:) / norm2(r_sp - r_vp)
      else
        surface_correction(vp, :) = surface_correction(vp, :) + this%poisson_prefactor * integrand(:)
      end if
    end do
  end do

  if (this%sys_grid%parallel_in_domains) then
    call this%sys_grid%allreduce(surface_correction)
  end if

  POP_SUB(X(compute_surface_correction_trans_field))
end subroutine X(compute_surface_correction_trans_field)

subroutine X(get_scalar_potential)(this, namespace, scalar_potential, total_field, apply_boundary)
  class(helmholtz_decomposition_t), intent(inout) :: this                !< Helmholtz object
  type(namespace_t),                intent(in)    :: namespace
  R_TYPE,                           intent(out)   :: scalar_potential(:) !< Vector potential obtained from the total_field.
  !!                                                                        It will be defined on the grid of the Helmholtz type
  !!                                                                        (this%grid), UP UNTIL THIS%GRID%NP:
  !!                                                                        scalar_potential(1:this%grid%np_part)
  R_TYPE,                           intent(inout) :: total_field(:,:)    !< Total field to decompose. To ensure the correct
  !!                                                                        behaviour this field must be correctly defined on the
  !!                                                                        grid of the system that calls Helmholtz (sys_gr),
  !!                                                                        UP UNTIL SYS_GR%NP_PART:
  !!                                                                        total_field(1:sys_gr%np_part, 1:sys_gr%box%dim)
  logical, optional,                intent(in)    :: apply_boundary      !< Should the divergence apply boundary conditions?

  logical              :: apply_boundary_
  R_TYPE, allocatable  :: support_field(:)

  PUSH_SUB(X(get_scalar_potential))

  SAFE_ALLOCATE(support_field(1:this%sys_grid%np_part))
  support_field = M_ZERO
  ! if not specified, we apply the boundary condition when computing the div
  apply_boundary_ = optional_default(apply_boundary, .true.)

  ! First of all, we have to compute the divergence of the field
  call X(derivatives_div)(this%sys_grid%der, total_field, support_field, set_bc = apply_boundary_)

  ! Then we solve poisson equation to compute helmholtz decomposition integral.
  ! The divergence is scalar, so we apply the poisson only in one direction
  scalar_potential = M_ZERO
  call X(poisson_solve)(this%poisson_solver, namespace, scalar_potential(:), support_field(:))

  if (this%compute_surface_correction) then
    support_field = M_ZERO
    ! Compute the surface correction
    call X(compute_surface_correction_long_field)(this, namespace, total_field, support_field)
    ! Add the surface correction to the scalar potential
    scalar_potential(:) = scalar_potential(:) - support_field(:)
  end if
  ! Add the prefactor 1/(4*pi)
  scalar_potential = R_TOTYPE(M_ONE/(M_FOUR * M_PI)) * scalar_potential

  SAFE_DEALLOCATE_A(support_field)

  POP_SUB(X(get_scalar_potential))
end subroutine X(get_scalar_potential)

!----------------------------------------------------------
subroutine X(get_long_field)(this, namespace, longitudinal_field, total_field, scalar_potential, apply_boundary)
  class(helmholtz_decomposition_t), intent(inout) :: this                    !< Helmholtz object
  type(namespace_t),                intent(in)    :: namespace
  R_TYPE,                           intent(out)   :: longitudinal_field(:,:) !< Transverse field obtained from the total_field. It
  !!                                                                            will be defined on the inner_grid of the Helmholtz
  !!                                                                            type (this%inner_grid), UP UNTIL INNER_GR%NP
  !!                                                                            longitudinal_field(1:this%inner_grid%np_part,
  !!                                                                                               1:this%inner_grid%box%dim)
  R_TYPE, optional,                 intent(inout) :: total_field(:,:)        !< Total field to decompose. To ensure the correct
  !!                                                                            behaviour this field must be correctly defined on
  !!                                                                            the grid of the system that calls Helmholtz (sys_gr)
  !!                                                                            UP UNTIL SYS_GR%NP_PART
  !!                                                                            total_field(1:sys_gr%np_part, 1:sys_gr%box%dim)
  R_TYPE, optional,                 intent(inout) :: scalar_potential(:)     !< Scalar potential from which long field whould be
  !!                                                                            computed
  logical, optional,                intent(in)    :: apply_boundary          !< Should the divergence apply boundary conditions?

  R_TYPE, allocatable  :: scal_pot(:)

  PUSH_SUB(X(get_long_field))

  if (present(scalar_potential)) then
    ! Compute longitudinal field from scalar potential
    call X(derivatives_grad)(this%sys_grid%der, scalar_potential, longitudinal_field, set_bc = .true.)
  else if (present(total_field)) then
    ! Compute longitudinal field from total field
    SAFE_ALLOCATE(scal_pot(1:this%sys_grid%np_part))
    ! First, get the scalar potential
    call this%get_scalar_potential(namespace, scal_pot, total_field, apply_boundary)

    ! Finally we compute the gradient to retrieve the curl-free term in the helmholtz decomposition
    call X(derivatives_grad)(this%sys_grid%der, scal_pot, longitudinal_field, set_bc = .true.)
    SAFE_DEALLOCATE_A(scal_pot)
  else
    message(1) = "Helmholtz decomposition error - you must pass either the total_field or the vector_potential"
    call messages_fatal(1)
  end if
  ! Add the - in front
  longitudinal_field = -longitudinal_field

  ! Finally, apply the mask to take care of the spikes generated by the curl
  call this%apply_inner_stencil_mask(longitudinal_field)

  POP_SUB(X(get_long_field))
end subroutine X(get_long_field)

! TODO: Issue 705 (ftroisi) - Check that inclusion of reduction over mesh does not affect the results. Add openmp support
subroutine X(compute_surface_correction_long_field)(this, namespace, field, surface_correction)
  class(helmholtz_decomposition_t), intent(in)  :: this                  !< Helmholtz object
  type(namespace_t),                intent(in)  :: namespace
  R_TYPE,                           intent(in)  :: field(:,:)            !< (1:this%sys_grid%np, 1:this%sys_grid%box%dim)
  R_TYPE,                           intent(out) :: surface_correction(:) !< (1:this%sys_grid%np)

  integer :: sp, vp, ii
  FLOAT   :: rr, surface_element
  FLOAT   :: r_vp(this%sys_grid%box%dim), r_sp(this%sys_grid%box%dim), normal_vector(this%sys_grid%box%dim)
  R_TYPE  :: integrand

  PUSH_SUB(X(compute_surface_correction_long_field))
  ASSERT(.not. this%sys_grid%use_curvilinear)
  surface_correction = M_ZERO

  ! First of all cycle over all surface points
  do ii = 1, SIZE(this%surface_points, dim = 1)
    ! Check that the point belongs to the surface and that is not the same of the volume point (avoid 0 at denominator)
    sp = this%surface_points(ii)
    ! Get coordinates of surface point
    call mesh_r(this%sys_grid, sp, rr, coords = r_sp)
    ! Then get the normal vector in that point and the Jacobian to account for the fact that the spacing between
    ! surface points is different than the spacing on the normal grid. Since surface points are always the same,
    ! it is enough to compute the normal vector and the Jacobian once
    call this%sys_grid%box%get_surface_point_info(r_sp, this%sys_grid%spacing(:), normal_vector(:), surface_element)

    ! Compute the cross prod between the normal vec and the field, times the surface_element (does not depend on the volume point)
    integrand = dot_product(normal_vector(:), field(sp, :)) * surface_element

    ! Then cycle over the volume points
    do vp = 1, this%sys_grid%np
      ! If vp == sp we are dealing with the singular point, which must be treated differently
      if (vp /= sp) then
        ! Get coordinates of volume point
        call mesh_r(this%sys_grid, vp, rr, coords = r_vp)
        ! Perform the integral
        surface_correction(vp) = surface_correction(vp) + integrand / norm2(r_sp - r_vp)
      else
        surface_correction(vp) = surface_correction(vp) + this%poisson_prefactor * integrand
      end if
    end do
  end do

  if (this%sys_grid%parallel_in_domains) then
    call this%sys_grid%allreduce(surface_correction)
  end if

  POP_SUB(X(compute_surface_correction_long_field))
end subroutine X(compute_surface_correction_long_field)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
