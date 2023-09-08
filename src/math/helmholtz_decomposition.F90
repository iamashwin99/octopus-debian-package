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

!> @brief The Helmholtz decomposition is intended to contain "only mathematical" functions and procedures to compute the
!! Helmholtz decomposition of a generic field.
!!
!! Given a generic vector field F, it is possibile to define a curl-free component div(phi) and a divergence free component
!! curl(A) such that \f$ \mathbf{F} = -div(phi) + curl(A) \f$, where
!! \f[
!!   \phi = \frac{1}{4\pi} \left( \int \frac{\nabla \cdot \mathbf{F}}{|\mathbf{r} -
!!   \mathbf{r_1}|} \,dr_1 - \oint \mathbf{n} \cdot \frac{\mathbf{F}}{|\mathbf{r} - \mathbf{r_1}|} \,dS_1  \right)
!! \f]
!!
!! \f[\mathbf{A} = \frac{1}{4\pi} \left( \int \frac{\nabla \times \mathbf{F}}{|\mathbf{r} - \mathbf{r_1}|} \,dr_1 - \oint \mathbf{n} \times \frac{\mathbf{F}}{|\mathbf{r} - \mathbf{r_1}|} \,dS_1  \right) \f]
!!
!! We start from a field defined on the system grid, i.e. the grid used for the Maxwell simulation.
!!
!! ## VECTOR and SCALAR POTENTIAL
!! If the user asks to compute the vector potential A or the scalar potential phi, then they are computed on the system grid.
!! All the points up to system_grid%np will have physical meaning
!! ## TRANSVERSE and LONGITUDINAL FIELDS
!! If the user asks to compute the transverse or the longitudinal field, then they are computed on the system grid. However,
!! the points the points contained in a layer of the width of the stencil between the innermost layer of points and last row
!! of points in the grid will be set to zero. This is a visual representation for stencil = 2:
!! 1 1 1 1 1 1 1
!! 1 1 1 1 1 1 1
!! 1 1 0 0 0 1 1
!! 1 1 0 0 0 1 1
!! 1 1 0 0 0 1 1
!! 1 1 1 1 1 1 1
!! 1 1 1 1 1 1 1
!! To understand this, one has to remember that when we get the scalar or vector potential, we know the solution only up to
!! system_grid%np. When we take the final divergence or curl, also the points from system_grid%np+1 to system_grid%np_part
!! matter (due to the way in which the finite differences are computed). Since potentials are set to zero in this non-physical
!! region, there might be spikes at the border. Therefore, we set to zero the points of the fields in the aforementioned mask,
!! so that we do not have spikes
!! ## SURFACE CORRECTION
!! In the equations above, there are two integrals, the volume one (always computed) and the surface one. The surface integral
!! is important when the total field is significantly different from zero at the boundaries of the box

! TODO: Issue 705 (ftroisi): Add a test for non Coulomb Gauge
module helmholtz_decomposition_m
  use blas_oct_m
  use box_oct_m
  use box_factory_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
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
  public  :: helmholtz_decomposition_t

  type helmholtz_decomposition_t
    private
    ! Helmholtz solution handling
    type(grid_t), pointer :: sys_grid                   !< Pointer to the system grid
    logical, public       :: compute_surface_correction !< Should the surface correction be computed?
    integer, allocatable  :: surface_points(:) !< Array containing the list of surface points of the box
    integer, allocatable  :: inner_stencil(:)  !< Mask for points that are less than a stencil away from sys_grid boundary
    ! Poisson Solver
    FLOAT                 :: poisson_prefactor !< Prefactor to handle singularities in the surface correction
    type(poisson_t)       :: poisson_solver    !< The poisson solver for the Helmholtz decomposition

  contains
    ! Init
    procedure :: init => helmholtz_decomposition_init
    ! Apply mask
    procedure :: dapply_inner_stencil_mask, zapply_inner_stencil_mask
    generic   :: apply_inner_stencil_mask => dapply_inner_stencil_mask, zapply_inner_stencil_mask
    ! Vector potential
    procedure :: dget_vector_potential, zget_vector_potential
    generic   :: get_vector_potential => dget_vector_potential, zget_vector_potential
    ! Scalar potential
    procedure :: dget_scalar_potential, zget_scalar_potential
    generic   :: get_scalar_potential => dget_scalar_potential, zget_scalar_potential
    ! Transverse field
    procedure :: dget_trans_field, zget_trans_field
    generic   :: get_trans_field => dget_trans_field, zget_trans_field
    ! Longitudinal field
    procedure :: dget_long_field, zget_long_field
    generic   :: get_long_field => dget_long_field, zget_long_field
    ! Finalize
    final     :: helmholtz_finalize
  end type helmholtz_decomposition_t

contains

  !> Initialize Helmholtz decomposition object
  subroutine helmholtz_decomposition_init(this, namespace, sys_grid, system_mc, space)
    class(helmholtz_decomposition_t), intent(inout) :: this         !< Helmholtz decomposition type to initialize
    type(namespace_t),                intent(in)    :: namespace
    type(grid_t), target,             intent(in)    :: sys_grid     !< Grid of the system calling Helmholtz
    type(multicomm_t),                intent(in)    :: system_mc    !< MultiCommunicator of the system calling Helmholtz
    type(space_t),                    intent(in)    :: space        !< Space of the system calling Helmholtz

    logical, allocatable  :: mask(:)
    logical               :: visualize_boxes
    integer               :: number_of_points

    PUSH_SUB(helmholtz_decomposition_init)

    ! Allocate pointer to system grid
    this%sys_grid => sys_grid

    ! 1. Initialize poisson solver
    call poisson_init(this%poisson_solver, namespace, space, sys_grid%der, system_mc, label = "Helmholtz")

    ! 2. Create a mask for the points contained in a layer of the width of the stencil between the innermost layer of points and
    !    last row of points in the grid.
    SAFE_ALLOCATE(mask(1:sys_grid%np))
    ! Get the mask
    mask = derivatives_get_inner_boundary_mask(sys_grid%der)
    number_of_points = count(mask)
    ! Store the indices of the mask points
    SAFE_ALLOCATE(this%inner_stencil(1:number_of_points))
    call get_indices_from_mask(sys_grid%np, mask, this%inner_stencil)

    ! 3. Surface correction TODO: See Issue 705 (@ftroisi)

    !%Variable SurfaceCorrection
    !%Type logical
    !%Default no
    !%Section Calculation Modes::Test
    !%Description
    !% Compute the surface correction for Helmholtz decomposition?
    !%End
    call parse_variable(namespace, 'SurfaceCorrection', .false., this%compute_surface_correction)
    if (this%compute_surface_correction) then
      call messages_not_implemented("Surface correction for Helmholtz decomposition")
    end if

    if (this%compute_surface_correction) then
      ! First of all we have to get the surface points
      mask = sys_grid%box%get_surface_points(namespace, sys_grid%spacing, sys_grid%np, sys_grid%x)
      number_of_points = count(mask)

      ! Allocate array surface_points and assign those points
      SAFE_ALLOCATE(this%surface_points(1:number_of_points))
      call get_indices_from_mask(sys_grid%np, mask, this%surface_points)

      ! Initialize the poisson prefactor to treat the singular points (for the surface correction)
      select case (sys_grid%box%dim)
      case (3)
        this%poisson_prefactor = M_TWO * M_PI * (M_THREE / (M_PI * M_FOUR))**(M_TWOTHIRD)
      case (2)
        this%poisson_prefactor = M_TWO * sqrt(M_PI)
      case (1)
        this%poisson_prefactor = M_ONE
      case default
        message(1) = "Internal error: surface correction for helmholtz decomposition can only be called for 1D, 2D or 3D."
        call messages_fatal(1, namespace = namespace)
      end select
    end if
    SAFE_DEALLOCATE_A(mask)

    ! 4. Optionally, the user may want to visualize the relevant regions used to compute the Helmholtz decomposition
    !%Variable HelmholtzVisualizeBoxes
    !%Type logical
    !%Default no
    !%Section Calculation Modes::Test
    !%Description
    !% If true, output the volume points for the three boxes of the Helmholtz surface correction.
    !% 1) The volume points of the system box
    !% 2) The inner mask for the system box. This region has the thickness of the stencil and it is used to set to zero
    !%    the longitudinal or transverse field after computing the final divergence or curl (to avoid spikes)
    !% 3) The surface points of the system box
    !%
    !%End
    call parse_variable(namespace, 'HelmholtzVisualizeBoxes', .false., visualize_boxes)
    if (visualize_boxes) then
      call helmholtz_visualize_boxes(this, namespace, space)
    end if

    POP_SUB(helmholtz_decomposition_init)
  end subroutine helmholtz_decomposition_init

  subroutine get_indices_from_mask(np, mask, indices)
    integer, intent(in)  :: np
    logical, intent(in)  :: mask(:)
    integer, intent(out) :: indices(:)

    integer :: ip, j
    PUSH_SUB(get_indices_from_mask)

    j = 1
    do ip = 1, np
      if (mask(ip)) then
        indices(j) = ip
        j = j + 1
      end if
    end do

    POP_SUB(get_indices_from_mask)
  end subroutine get_indices_from_mask

  subroutine helmholtz_finalize(this)
    type(helmholtz_decomposition_t), intent(inout) :: this
    PUSH_SUB(helmholtz_finalize)

    call poisson_end(this%poisson_solver)

    SAFE_DEALLOCATE_A(this%inner_stencil)
    SAFE_DEALLOCATE_A(this%surface_points)

    POP_SUB(helmholtz_finalize)
  end subroutine helmholtz_finalize

  !> Visualise boxes for use in Helmholtz Decomposition
  subroutine helmholtz_visualize_boxes(this, namespace, space)
    type(helmholtz_decomposition_t), intent(inout) :: this
    type(namespace_t), intent(in)                  :: namespace
    type(space_t),     intent(in)                  :: space

    FLOAT                :: coords(space%dim), normal_vector(space%dim), surface_element, ra
    integer              :: ip, ii, iunit

    PUSH_SUB(helmholtz_visualize_boxes)

    ! SYSTEM BOX - Volume points
    iunit = io_open("system_box_volume_points", namespace, action='write')
    write(iunit, '(a6, 5X, a1, 12X, a1, 12X, a1)') "#point", "x", "y", "z"
    do ip = 1, this%sys_grid%np
      call mesh_r(this%sys_grid, ip, ra, coords = coords)
      write(iunit, '(i7, 3f12.7)') ip, coords(:)
    end do
    call io_close(iunit)

    ! SYSTEM BOX - Mask for inner stencil
    iunit = io_open("system_box_inner_stencil", namespace, action='write')
    write(iunit, '(a6, 5X, a1, 12X, a1, 12X, a1)') "#point", "x", "y", "z"
    do ip = 1, SIZE(this%inner_stencil, dim = 1)
      ii = this%inner_stencil(ip)
      call mesh_r(this%sys_grid, ii, ra, coords = coords)
      write(iunit, '(i7, 3f12.7)') ip, coords(:)
    end do
    call io_close(iunit)

    ! SYSTEM BOX - Surface Points
    if (this%compute_surface_correction) then
      iunit = io_open("system_box_surface_points", namespace, action='write')
      write(iunit, '(a,12X,a,12X,a,12X,a,12X,a,12X,a,12X,a)') "x", "y", "z", "norm_x", "norm_y", "norm_z", "surf_element"
      do ip = 1, SIZE(this%surface_points, dim = 1)
        ii = this%surface_points(ip)
        call mesh_r(this%sys_grid, ii, ra, coords = coords)
        call this%sys_grid%box%get_surface_point_info(coords, this%sys_grid%spacing, normal_vector, surface_element)
        write(iunit, '(7f12.7)') coords(:), normal_vector(:), surface_element
      end do
      call io_close(iunit)
    end if

    POP_SUB(helmholtz_visualize_boxes)
  end subroutine helmholtz_visualize_boxes

#include "undef.F90"
#include "complex.F90"
#include "helmholtz_decomposition_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "helmholtz_decomposition_inc.F90"

end module helmholtz_decomposition_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
