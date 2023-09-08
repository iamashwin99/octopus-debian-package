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

#include "global.h"

!> @brief Test suit for the Helmholtz decomposition module

module helmholtz_decomposition_test_m
  use blas_oct_m
  use box_oct_m
  use box_factory_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use helmholtz_decomposition_m
  use io_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::               &
    helmholtz_test_t,     &
    hertzian_dipole_test, &
    gaussian_test

  type helmholtz_test_t
    private
    FLOAT, allocatable :: total_field(:,:)
    FLOAT, allocatable :: trans_field(:,:)
    FLOAT, allocatable :: long_field(:,:)
    FLOAT, allocatable :: vector_potential(:,:)
    FLOAT, allocatable :: scalar_potential(:)
    FLOAT              :: abs_norm_self_consistency, abs_norm_long_field, abs_norm_trans_field, abs_norm_scal_pot, abs_norm_vec_pot
    FLOAT              :: rel_norm_self_consistency, rel_norm_long_field, rel_norm_trans_field, rel_norm_scal_pot, rel_norm_vec_pot
  contains
    procedure :: init => helmholtz_test_init
    final     :: helmholtz_test_finalize
  end type helmholtz_test_t

contains

  subroutine helmholtz_test_init(this, np_part, dim)
    class(helmholtz_test_t), intent(inout) :: this
    integer,                 intent(in)    :: np_part
    integer,                 intent(in)    :: dim

    PUSH_SUB(helmholtz_test_init)

    ! Scalar potential
    SAFE_ALLOCATE(this%scalar_potential(1:np_part))
    this%scalar_potential = M_ZERO
    ! Vector potential
    SAFE_ALLOCATE(this%vector_potential(1:np_part, 1:dim))
    this%vector_potential = M_ZERO
    ! Total Field
    SAFE_ALLOCATE(this%total_field(1:np_part, 1:dim))
    this%total_field = M_ZERO
    ! Transverse field
    SAFE_ALLOCATE(this%trans_field(1:np_part, 1:dim))
    this%trans_field = M_ZERO
    ! Longitudinal field
    SAFE_ALLOCATE(this%long_field(1:np_part, 1:dim))
    this%long_field = M_ZERO

    POP_SUB(helmholtz_test_init)
  end subroutine helmholtz_test_init

  subroutine helmholtz_test_finalize(this)
    type(helmholtz_test_t), intent(inout) :: this

    PUSH_SUB(helmholtz_test_finalize)

    SAFE_DEALLOCATE_A(this%total_field)
    SAFE_DEALLOCATE_A(this%trans_field)
    SAFE_DEALLOCATE_A(this%long_field)
    SAFE_DEALLOCATE_A(this%vector_potential)
    SAFE_DEALLOCATE_A(this%scalar_potential)

    POP_SUB(helmholtz_test_finalize)
  end subroutine helmholtz_test_finalize

  subroutine hertzian_dipole_test(helmholtz, sys_grid, namespace, space)
    class(helmholtz_decomposition_t), intent(inout) :: helmholtz
    type(grid_t),                     intent(in)    :: sys_grid
    type(namespace_t),                intent(in)    :: namespace
    type(space_t),                    intent(in)    :: space

    ! Fields
    type(helmholtz_test_t) :: helm_analytical
    type(helmholtz_test_t) :: helm_comp_test_1
    type(helmholtz_test_t) :: helm_comp_test_2

    character(len=15)      :: test_case

    PUSH_SUB(hertzian_dipole_test)
    ASSERT(space%dim == 3)

    ! Allocate and initialize the fields
    call helm_analytical%init(sys_grid%np_part, sys_grid%box%dim)
    call helm_comp_test_1%init(sys_grid%np_part, sys_grid%box%dim)
    ! Define the field to be tested (Hertzian dipole).
    ! Reference: The Helmholtz Decomposition and the Coulomb Gauge, Kirk T. McDonald, Joseph Henry Laboratories, Princeton (2008)
    ! http://kirkmcd.princeton.edu/examples/helmholtz.pdf
    call init_hertzian_dipole_field(helm_analytical, sys_grid)
    call io_mkdir("hertzian_dipole_test", namespace=namespace)
    ! Output analytical fields
    call io_mkdir("hertzian_dipole_test/exact_fields", namespace=namespace)
    call output_fields(sys_grid, helm_analytical, space, namespace, "hertzian_dipole_test/exact_fields")

    ! TEST 1
    ! Compute the transverse and longitudinal fields numerically
    call helmholtz%get_trans_field(namespace, helm_comp_test_1%trans_field, total_field=helm_analytical%total_field, &
      apply_boundary=.false.)
    call helmholtz%get_long_field(namespace, helm_comp_test_1%long_field, total_field=helm_analytical%total_field, &
      apply_boundary=.false.)
    ! Compute the scalar and vector potential numerically
    call helmholtz%get_vector_potential(namespace, helm_comp_test_1%vector_potential, helm_analytical%total_field, &
      apply_boundary=.false.)
    call helmholtz%get_scalar_potential(namespace, helm_comp_test_1%scalar_potential, helm_analytical%total_field, &
      apply_boundary=.false.)
    ! Output fields
    test_case = "no_surf_corr"
    if (helmholtz%compute_surface_correction) test_case = "with_surf_corr"
    call io_mkdir("hertzian_dipole_test/" // trim(test_case), namespace=namespace)
    call output_fields(sys_grid,  helm_comp_test_1, space, namespace, "hertzian_dipole_test/" // trim(test_case))
    call compute_norms(helmholtz,sys_grid, helm_analytical, helm_comp_test_1, namespace, "hertzian_dipole_test/"//trim(test_case))

    ! TEST 2
    ! If the TEST 1 was done with surface correction, then perform a second one without surface correction. This allows to test
    ! that the surface correction actualy brings some improvements
    if (helmholtz%compute_surface_correction) then
      helmholtz%compute_surface_correction = .false.
      call helm_comp_test_2%init(sys_grid%np_part, sys_grid%box%dim)
      ! Compute the transverse and longitudinal fields numerically
      call helmholtz%get_trans_field(namespace, helm_comp_test_2%trans_field, total_field=helm_analytical%total_field, &
        apply_boundary=.false.)
      call helmholtz%get_long_field(namespace, helm_comp_test_2%long_field, total_field=helm_analytical%total_field, &
        apply_boundary=.false.)
      ! Compute the scalar and vector potential numerically
      call helmholtz%get_vector_potential(namespace, helm_comp_test_2%vector_potential, helm_analytical%total_field, &
        apply_boundary=.false.)
      call helmholtz%get_scalar_potential(namespace, helm_comp_test_2%scalar_potential, helm_analytical%total_field, &
        apply_boundary=.false.)
      ! Output fields
      test_case = "no_surf_corr"
      call io_mkdir("hertzian_dipole_test/" // trim(test_case), namespace=namespace)
      call output_fields(sys_grid, helm_comp_test_2, space, namespace, "hertzian_dipole_test/" // trim(test_case))
      call compute_norms(helmholtz,sys_grid,helm_analytical, helm_comp_test_2, namespace, "hertzian_dipole_test/"//trim(test_case))
      ! After the test is done, reset the compute_surface_correction flag
      helmholtz%compute_surface_correction = .true.
      ! If we are here, we did a test with surface correction and one without. So we should test that with surface correction
      ! we got better results
      call compare_norms(helm_comp_test_1, helm_comp_test_2, namespace, "hertzian_dipole_test/")
    end if

    POP_SUB(hertzian_dipole_test)
  end subroutine hertzian_dipole_test

  subroutine init_hertzian_dipole_field(this, sys_grid)
    class(helmholtz_test_t), intent(inout) :: this
    type(grid_t),            intent(in)    :: sys_grid

    FLOAT              :: lambda, omega, kk, rr, time, norm
    FLOAT              :: xx(sys_grid%box%dim), xx_origin(sys_grid%box%dim), pp(sys_grid%box%dim), cross_p(sys_grid%box%dim)
    integer            :: ip
    FLOAT, parameter   :: nearly_zero = CNST(1e-14)

    PUSH_SUB(init_hertzian_dipole_field)

    ! Parameters that define the dipole
    lambda = CNST(1.55)
    omega = M_TWO * M_PI * (P_c / lambda)
    kk = M_TWO * M_PI / lambda
    time = M_ZERO
    pp = cos(omega * time)
    xx_origin(2:sys_grid%box%dim) = M_ZERO
    xx_origin(1) = CNST(3.11)

    do ip = 1, sys_grid%np_part
      call mesh_r(sys_grid, ip, rr, coords = xx, origin = xx_origin)
      if (rr >= sys_grid%spacing(1)) then
        xx = xx / rr
        ! Longitudinal field - equation 31 of reference paper
        this%long_field(ip, :) = (M_THREE * dot_product(pp, xx) * xx - pp) * cos(omega*time)/(rr**3)
        ! Transverse field - equation 32 of reference paper
        cross_p = dcross_product(xx, pp)
        this%trans_field(ip, :) = (kk**2) * cross_p * (cos(kk*rr - omega*time)/rr - sin(kk*rr - omega*time)/(kk * rr**2))
        ! Vector potential - equation 53 of reference paper
        this%vector_potential(ip, :) = kk * (pp - dot_product(pp, xx)*xx) * sin(kk*rr - omega*time)/rr
        this%vector_potential(ip, :) = this%vector_potential(ip, :) * (pp - M_THREE*dot_product(pp, xx)*xx) * &
          (cos(kk*rr - omega*time)/(rr**2) - (sin(kk*rr - omega*time) - sin(omega*time))/(kk * rr**3))
        ! Scalar potential - equation 51 of reference paper
        this%scalar_potential(ip) = dot_product(pp, xx) * cos(omega*time) / (rr**2)
        ! Total field
        this%total_field = this%trans_field + this%long_field
      end if
    end do
    ! Check that the exact fields ere defined correctly
    norm = dmf_nrm2(sys_grid, sys_grid%box%dim, this%total_field - this%trans_field - this%long_field)
    ASSERT(norm < nearly_zero)

    POP_SUB(init_hertzian_dipole_field)
  end subroutine init_hertzian_dipole_field

  subroutine gaussian_test(helmholtz, sys_grid, namespace, space)
    class(helmholtz_decomposition_t), intent(inout) :: helmholtz
    type(grid_t),                     intent(in)    :: sys_grid
    type(namespace_t),                intent(in)    :: namespace
    type(space_t),                    intent(in)    :: space

    ! Fields
    type(helmholtz_test_t) :: helm_analytical
    type(helmholtz_test_t) :: helm_comp_test_1
    type(helmholtz_test_t) :: helm_comp_test_2

    character(len=15)      :: test_case

    PUSH_SUB(gaussian_test)
    ASSERT(space%dim == 3)

    ! Allocate and initialize the fields
    call helm_analytical%init(sys_grid%np_part, sys_grid%box%dim)
    call helm_comp_test_1%init(sys_grid%np_part, sys_grid%box%dim)

    ! Define the field to be tested
    ! https://gitlab.com/octopus-code/octopus-hackathons/-/blob/main/01-full-minimal-coupling-2022-07-04-2022-07-08/documentation/francesco/Helmholtz_Decomposition_Coulomb_Gauge.pdf
    call init_gaussian_field(helm_analytical, sys_grid)
    call io_mkdir("gaussian_field_test", namespace=namespace)
    ! Output analytical fields
    call io_mkdir("gaussian_field_test/exact_fields", namespace=namespace)
    call output_fields(sys_grid, helm_analytical, space, namespace, "gaussian_field_test/exact_fields")

    ! TEST 1
    ! Compute the transverse and longitudinal fields numerically
    call helmholtz%get_trans_field(namespace, helm_comp_test_1%trans_field, total_field=helm_analytical%total_field, &
      apply_boundary=.false.)
    call helmholtz%get_long_field(namespace, helm_comp_test_1%long_field, total_field=helm_analytical%total_field, &
      apply_boundary=.false.)
    ! Compute the scalar and vector potential numerically
    call helmholtz%get_vector_potential(namespace, helm_comp_test_1%vector_potential, helm_analytical%total_field, &
      apply_boundary=.false.)
    call helmholtz%get_scalar_potential(namespace, helm_comp_test_1%scalar_potential, helm_analytical%total_field, &
      apply_boundary=.false.)
    ! Output fields
    test_case = "no_surf_corr"
    if (helmholtz%compute_surface_correction) test_case = "with_surf_corr"
    call io_mkdir("gaussian_field_test/" // trim(test_case), namespace=namespace)
    call output_fields(sys_grid,  helm_comp_test_1, space, namespace, "gaussian_field_test/" // trim(test_case))
    call compute_norms(helmholtz, sys_grid, helm_analytical, helm_comp_test_1, namespace, "gaussian_field_test/"//trim(test_case))

    ! TEST 2
    ! If the TEST 1 was done with surface correction, then perform a second one without surface correction. This allows to test
    ! that the surface correction actualy brings some improvements
    if (helmholtz%compute_surface_correction) then
      helmholtz%compute_surface_correction = .false.
      call helm_comp_test_2%init(sys_grid%np_part, sys_grid%box%dim)
      ! Compute the transverse and longitudinal fields numerically
      call helmholtz%get_trans_field(namespace, helm_comp_test_2%trans_field, total_field=helm_analytical%total_field, &
        apply_boundary=.false.)
      call helmholtz%get_long_field(namespace, helm_comp_test_2%long_field, total_field=helm_analytical%total_field, &
        apply_boundary=.false.)
      ! Compute the scalar and vector potential numerically
      call helmholtz%get_vector_potential(namespace, helm_comp_test_2%vector_potential, helm_analytical%total_field, &
        apply_boundary=.false.)
      call helmholtz%get_scalar_potential(namespace, helm_comp_test_2%scalar_potential, helm_analytical%total_field, &
        apply_boundary=.false.)
      ! Output fields
      test_case = "no_surf_corr"
      call io_mkdir("gaussian_field_test/" // trim(test_case), namespace=namespace)
      call output_fields(sys_grid, helm_comp_test_2, space, namespace, "gaussian_field_test/" // trim(test_case))
      call compute_norms(helmholtz,sys_grid, helm_analytical, helm_comp_test_2, namespace, "gaussian_field_test/"//trim(test_case))
      ! After the test is done, reset the compute_surface_correction flag
      helmholtz%compute_surface_correction = .true.
      ! If we are here, we did a test with surface correction and one without. So we should test that with surface correction
      ! we got better results
      call compare_norms(helm_comp_test_1, helm_comp_test_2, namespace, "gaussian_field_test/")
    end if

    POP_SUB(gaussian_test)
  end subroutine gaussian_test

  subroutine init_gaussian_field(this, sys_grid)
    class(helmholtz_test_t), intent(inout) :: this
    type(grid_t),            intent(in)    :: sys_grid

    FLOAT              :: alpha, beta, a_prefactor, b_prefactor, ra, rb, norm
    FLOAT              :: aa(sys_grid%box%dim), a_origin(sys_grid%box%dim), bb(sys_grid%box%dim), b_origin(sys_grid%box%dim)
    integer            :: ip

    PUSH_SUB(init_gaussian_field)

    ! Parameters that define the gaussian
    alpha = CNST(0.3)
    a_prefactor = CNST(1.5)
    a_origin = M_ONE / 10
    beta = CNST(0.3)
    b_prefactor = CNST(1.7)
    b_origin = -M_ONE / 10

    ! Define the field to be tested
    do ip = 1, sys_grid%np_part
      call mesh_r(sys_grid, ip, ra, coords = aa, origin = a_origin)
      ! Longitudinal field
      this%long_field(ip, :) = (2 * a_prefactor / (alpha**2)) * aa(:) * exp(-(ra / alpha)**2)
      ! Analytical scalar potential
      this%scalar_potential(ip) = a_prefactor * exp(-(ra / alpha)**2)

      call mesh_r(sys_grid, ip, rb, coords = bb, origin = b_origin)
      ! Transverse field
      this%trans_field(ip, 1) = - (2 * b_prefactor / (beta**2)) * bb(1) * bb(3) * exp(-(rb / beta)**2)
      this%trans_field(ip, 2) = - (2 * b_prefactor / (beta**2)) * bb(2) * bb(3) * exp(-(rb / beta)**2)
      this%trans_field(ip, 3) = (2 * b_prefactor / (beta**2)) * (bb(1)**2 + bb(2)**2 - beta**2) * exp(-(rb / beta)**2)
      ! Analytical vector potential
      this%vector_potential(ip, 1) = b_prefactor * bb(2) * exp(-(rb / beta)**2)
      this%vector_potential(ip, 2) = - b_prefactor * bb(1) * exp(-(rb / beta)**2)

      ! Total field
      this%total_field(ip, :) = this%long_field(ip, :) + this%trans_field(ip, :)
    end do

    ! Check that the exact fields ere defined correctly
    norm = dmf_nrm2(sys_grid, sys_grid%box%dim, this%total_field - this%trans_field - this%long_field)
    ASSERT(norm < M_EPSILON)

    POP_SUB(init_gaussian_field)
  end subroutine init_gaussian_field

  subroutine output_fields(sys_grid, test, space, namespace, path)
    type(grid_t),            intent(in) :: sys_grid
    type(helmholtz_test_t),  intent(in) :: test
    type(space_t),           intent(in) :: space
    type(namespace_t),       intent(in) :: namespace
    character(len=*),        intent(in) :: path

    integer            :: ierr

    PUSH_SUB(output_fields)

    ! Output the values of the fields on the grid
    call io_function_output_vector (io_function_fill_how('PlaneZ'), ".", trim(path) // "/total_field", namespace, space, &
      sys_grid, test%total_field, unit_one, ierr)
    call io_function_output_vector (io_function_fill_how('PlaneZ'), ".", trim(path) // "/trans_field", namespace, space, &
      sys_grid, test%trans_field, unit_one, ierr)
    call io_function_output_vector (io_function_fill_how('PlaneZ'), ".", trim(path) // "/long_field", namespace, space, &
      sys_grid, test%long_field, unit_one, ierr)
    call io_function_output_vector (io_function_fill_how('PlaneZ'), ".", trim(path) // "/vector_potential", namespace, space, &
      sys_grid, test%vector_potential, unit_one, ierr)
    call dio_function_output (io_function_fill_how('PlaneZ'), ".", trim(path) // "/scalar_potential", namespace, space, &
      sys_grid, test%scalar_potential, unit_one, ierr)

    POP_SUB(output_fields)
  end subroutine output_fields

  subroutine compute_norms(helmholtz, sys_grid, analytical, test, namespace, path)
    type(helmholtz_decomposition_t), intent(inout) :: helmholtz
    type(grid_t),                    intent(in)    :: sys_grid
    type(helmholtz_test_t),          intent(in)    :: analytical
    type(helmholtz_test_t),          intent(inout) :: test
    type(namespace_t),               intent(in)    :: namespace
    character(len=*),                intent(in)    :: path

    integer            :: iunit, dim, np
    FLOAT, allocatable :: local_field(:,:)

    PUSH_SUB(compute_norms)
    SAFE_ALLOCATE(local_field(1:sys_grid%np, 1:sys_grid%box%dim))

    dim = sys_grid%box%dim ! All boxes have the same dimension
    np = sys_grid%np

    ! First, compute the norms
    ! Transverse field
    local_field = analytical%trans_field(1:np,:)
    call helmholtz%apply_inner_stencil_mask(local_field)
    test%abs_norm_trans_field = dmf_nrm2(sys_grid, dim, test%trans_field(1:np,:) - local_field(1:np,:))
    test%rel_norm_trans_field = test%abs_norm_trans_field / dmf_nrm2(sys_grid, dim, local_field(1:np, :))
    ! Longitudinal field
    local_field = analytical%long_field(1:np,:)
    call helmholtz%apply_inner_stencil_mask(local_field)
    test%abs_norm_long_field = dmf_nrm2(sys_grid, dim, test%long_field(1:np, :) - local_field(1:np, :))
    test%rel_norm_long_field = test%abs_norm_long_field / dmf_nrm2(sys_grid, dim, local_field(1:np, :))
    ! Vector potential
    test%abs_norm_vec_pot = dmf_nrm2(sys_grid, dim, test%vector_potential(1:np,:) - analytical%vector_potential(1:np,:))
    test%rel_norm_vec_pot = test%abs_norm_vec_pot / dmf_nrm2(sys_grid, dim, analytical%vector_potential(1:np, :))
    ! Scalar potential
    test%abs_norm_scal_pot = dmf_nrm2(sys_grid, test%scalar_potential(1:np) - analytical%scalar_potential(1:np))
    test%rel_norm_scal_pot = test%abs_norm_scal_pot / dmf_nrm2(sys_grid, analytical%scalar_potential(1:np))
    ! Self consistency
    local_field = analytical%total_field(1:np,:)
    call helmholtz%apply_inner_stencil_mask(local_field)
    test%abs_norm_self_consistency = dmf_nrm2(sys_grid,dim,local_field(1:np,:)-test%trans_field(1:np,:)-test%long_field(1:np,:))
    test%rel_norm_self_consistency = test%abs_norm_self_consistency / dmf_nrm2(sys_grid, dim, local_field(1:np, :))

    SAFE_DEALLOCATE_A(local_field)

    ! Output the norms to file
    iunit = io_open(path // "/norms", namespace, action='write')
    write(iunit, '(a,f19.13)') 'Transverse field test (abs.). = ', test%abs_norm_trans_field
    write(iunit, '(a,f19.13)') 'Transverse field test (rel.). = ', test%rel_norm_trans_field
    write(iunit, '(a,f19.13)') 'Longitudinal field test (abs.). = ', test%abs_norm_long_field
    write(iunit, '(a,f19.13)') 'Longitudinal field test (rel.). = ', test%rel_norm_long_field
    write(iunit, '(a,f19.13)') 'Vector potential test (abs.). = ', test%abs_norm_vec_pot
    write(iunit, '(a,f19.13)') 'Vector potential test (rel.). = ', test%rel_norm_vec_pot
    write(iunit, '(a,f19.13)') 'Scalar potential test (abs.). = ', test%abs_norm_scal_pot
    write(iunit, '(a,f19.13)') 'Scalar potential test (rel.). = ', test%rel_norm_scal_pot
    write(iunit, '(a,f19.13)') 'Self consistency test (abs.). = ', test%abs_norm_self_consistency
    write(iunit, '(a,f19.13)') 'Self consistency test (rel.). = ', test%rel_norm_self_consistency

    call io_close(iunit)

    POP_SUB(compute_norms)
  end subroutine compute_norms

  subroutine compare_norms(test_with_sc, test_no_sc, namespace, path)
    type(helmholtz_test_t),  intent(in) :: test_with_sc
    type(helmholtz_test_t),  intent(in) :: test_no_sc
    type(namespace_t),       intent(in) :: namespace
    character(len=*),        intent(in) :: path

    integer                      :: iunit
    character(len=24), parameter :: form = '(a,f19.13,a,f19.13)'

    PUSH_SUB(compare_norms)

    ! Output the norms to file
    iunit = io_open(path // "norm_comparison", namespace, action='write')
    ! Trans field
    write(iunit, form) 'Transverse field test (abs.). No surf correction = ', test_no_sc%abs_norm_trans_field, &
      " ;With surf correction = ", test_with_sc%abs_norm_trans_field
    write(iunit, form) 'Transverse field test (rel.). No surf correction = ', test_no_sc%rel_norm_trans_field, &
      " ;With surf correction = ", test_with_sc%rel_norm_trans_field
    ! Long field
    write(iunit, form) 'Longitudinal field test (abs.). No surf correction = ', test_no_sc%abs_norm_long_field, &
      " ;With surf correction = ", test_with_sc%abs_norm_long_field
    write(iunit, form) 'Longitudinal field test (rel.). No surf correction = ', test_no_sc%rel_norm_long_field, &
      " ;With surf correction = ", test_with_sc%rel_norm_long_field
    ! Vector potential
    write(iunit, form) 'Vector potential test (abs.). No surf correction = ', test_no_sc%abs_norm_vec_pot, &
      " ;With surf correction = ", test_with_sc%abs_norm_vec_pot
    write(iunit, form) 'Vector potential test (rel.). No surf correction = ', test_no_sc%rel_norm_vec_pot, &
      " ;With surf correction = ", test_with_sc%rel_norm_vec_pot
    ! Scalar potential
    write(iunit, form) 'Scalar potential test (abs.). No surf correction = ', test_no_sc%abs_norm_scal_pot, &
      " ;With surf correction = ", test_with_sc%abs_norm_scal_pot
    write(iunit, form) 'Scalar potential test (rel.). No surf correction = ', test_no_sc%rel_norm_scal_pot, &
      " ;With surf correction = ", test_with_sc%rel_norm_scal_pot
    ! Self consistency
    write(iunit, form) 'Self consistency test (abs.). No surf correction = ', test_no_sc%abs_norm_self_consistency, &
      " ;With surf correction = ", test_with_sc%abs_norm_self_consistency
    write(iunit, form) 'Self consistency test (rel.). No surf correction = ', test_no_sc%rel_norm_self_consistency, &
      " ;With surf correction = ", test_with_sc%rel_norm_self_consistency
    call io_close(iunit)

    ! Trans field
    ASSERT(test_with_sc%abs_norm_trans_field < test_no_sc%abs_norm_trans_field)
    ASSERT(test_with_sc%rel_norm_trans_field < test_no_sc%rel_norm_trans_field)
    ! Long field
    ASSERT(test_with_sc%abs_norm_long_field < test_no_sc%abs_norm_long_field)
    ASSERT(test_with_sc%rel_norm_long_field < test_no_sc%rel_norm_long_field)
    ! Vector potential
    ASSERT(test_with_sc%abs_norm_vec_pot < test_no_sc%abs_norm_vec_pot)
    ASSERT(test_with_sc%rel_norm_vec_pot < test_no_sc%rel_norm_vec_pot)
    ! Scalar potential
    ASSERT(test_with_sc%abs_norm_scal_pot < test_no_sc%abs_norm_scal_pot)
    ASSERT(test_with_sc%rel_norm_scal_pot < test_no_sc%rel_norm_scal_pot)
    ! Self Consistency
    ASSERT(test_with_sc%abs_norm_self_consistency < test_no_sc%abs_norm_self_consistency)
    ASSERT(test_with_sc%rel_norm_self_consistency < test_no_sc%rel_norm_self_consistency)

    POP_SUB(compare_norms)
  end subroutine compare_norms

#include "undef.F90"
#include "complex.F90"

#include "undef.F90"
#include "real.F90"

end module helmholtz_decomposition_test_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
