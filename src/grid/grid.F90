!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 S. Ohlmann
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

!> \ingroup Fortran_Module
module grid_oct_m
  use affine_coordinates_oct_m
  use box_oct_m
  use box_factory_oct_m
  use box_image_oct_m
  use cartesian_oct_m
  use coordinate_system_oct_m
  use cube_oct_m
  use curv_briggs_oct_m
  use curv_gygi_oct_m
  use curv_modine_oct_m
  use derivatives_oct_m
  use global_oct_m
  use lattice_vectors_oct_m
  use mesh_oct_m
  use mesh_init_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use nl_operator_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use stencil_oct_m
  use stencil_cube_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                        &
    grid_t,                        &
    grid_init_stage_1,             &
    grid_init_from_grid_stage_1,   &
    grid_init_stage_2,             &
    grid_end,                      &
    grid_write_info,               &
    dgrid_symmetrize_scalar_field, &
    zgrid_symmetrize_scalar_field, &
    dgrid_symmetrize_vector_field, &
    zgrid_symmetrize_vector_field, &
    dgrid_symmetrize_single,       &
    zgrid_symmetrize_single

  !> \ingroup Fortran_Class
  !! @brief Description of the grid, containing information on derivatives, stencil, and symmetries
  type, extends(mesh_t) :: grid_t
    ! Components are public by default
    type(derivatives_t)         :: der
    type(stencil_t)             :: stencil
    type(symmetries_t)          :: symm
    type(symmetrizer_t)         :: symmetrizer
  end type grid_t

  integer, parameter :: &
    CURV_AFFINE  = 1,   &
    CURV_GYGI    = 2,   &
    CURV_BRIGGS  = 3,   &
    CURV_MODINE  = 4

contains

  !-------------------------------------------------------------------
  subroutine grid_init_stage_1(gr, namespace, space, symm, latt, n_sites, site_position)
    type(grid_t),                      intent(inout) :: gr
    type(namespace_t),                 intent(in)    :: namespace
    type(space_t),                     intent(in)    :: space
    type(symmetries_t),      optional, intent(in)    :: symm
    type(lattice_vectors_t), optional, intent(in)    :: latt
    integer,                 optional, intent(in)    :: n_sites
    FLOAT,                   optional, intent(in)    :: site_position(:,:)

    type(stencil_t) :: cube
    integer :: enlarge(1:MAX_DIM)
    type(block_t) :: blk
    integer :: idir
    FLOAT :: grid_spacing(1:MAX_DIM)

    PUSH_SUB(grid_init_stage_1)

    gr%box => box_factory_create(namespace, space, latt=latt, n_sites=n_sites, site_position=site_position)

    if (present(symm)) then
      gr%symm = symm
    end if

    if (.not. parse_is_defined(namespace, 'Spacing')) then
      call messages_input_error(namespace, 'Spacing', 'Spacing needs to be set in the input file.')
    end if

    !%Variable Spacing
    !%Type float
    !%Section Mesh
    !%Description
    !% The spacing between the points in the mesh. This controls the
    !% quality of the discretization: smaller spacing gives more
    !% precise results but increased computational cost.
    !%
    !% When using curvilinear coordinates, this is a canonical spacing
    !% that will be changed locally by the transformation. In periodic
    !% directions, your spacing may be slightly different than what
    !% you request here, since the box size must be an integer
    !% multiple of the spacing.
    !%
    !% The default value is defined by the species if only default pseudopotentials are used
    !% or by the image resolution if <tt>BoxShape = box_image</tt>. Otherwise, there is
    !% no default.
    !%
    !% It is possible to have a different spacing in each one of the Cartesian directions
    !% if we define <tt>Spacing</tt> as block of the form
    !%
    !% <tt>%Spacing
    !% <br>&nbsp;&nbsp;spacing_x | spacing_y | spacing_z
    !% <br>%</tt>
    !%End

    grid_spacing(space%dim+1:MAX_DIM) = -M_ONE
    if(parse_block(namespace, 'Spacing', blk) == 0) then
      if(parse_block_cols(blk,0) < space%dim) call messages_input_error(namespace, 'Spacing')
      do idir = 1, space%dim
        call parse_block_float(blk, 0, idir - 1, grid_spacing(idir), units_inp%length)
      end do
      call parse_block_end(blk)
    else
      call parse_variable(namespace, 'Spacing', -M_ONE, grid_spacing(1), units_inp%length)
      grid_spacing(1:space%dim) = grid_spacing(1)
    end if

    if (associated(gr%box)) then
      select type (box => gr%box)
      type is (box_image_t)
        do idir = 1, space%dim
          ! default grid_spacing is determined from the pixel size such that one grid point = one pixel.
          if (grid_spacing(idir) < M_ZERO) then
            grid_spacing(idir) = box%pixel_size(idir)
          end if
        end do
      end select
    end if

    if (any(grid_spacing(1:space%dim) < M_EPSILON)) then
      message(1) = " Your input for 'Spacing' is negative or zero."
      call messages_fatal(1, namespace=namespace)
    end if

    !%Variable PeriodicBoundaryMask
    !%Type block
    !%Section Mesh
    !%Description
    !% (Experimental) Defines a mask for which periodic boundaries are replaced by zero boundary conditions.
    !%End
    if (parse_block(namespace, 'PeriodicBoundaryMask', blk) < 0) then
      gr%masked_periodic_boundaries = .false.
    else
      gr%masked_periodic_boundaries = .true.
      call parse_block_string(blk, 0, 0, gr%periodic_boundary_mask)
      call messages_experimental('PeriodicBoundaryMask')
    end if

    ! Initialize coordinate system
    call initialize_coordinate_system(gr, namespace, space, latt, n_sites, site_position)

    ! initialize derivatives
    call derivatives_init(gr%der, namespace, space, gr%coord_system)
    ! the stencil used to generate the grid is a union of a cube (for
    ! multigrid) and the Laplacian
    call stencil_cube_get_lapl(cube, space%dim, order = 2)
    call stencil_union(space%dim, cube, gr%der%lapl%stencil, gr%stencil)
    call stencil_end(cube)

    enlarge = 0
    enlarge(1:space%dim) = 2
    enlarge = max(enlarge, gr%der%n_ghost)

    call mesh_init_stage_1(gr, namespace, space, gr%box, gr%coord_system, grid_spacing, enlarge)
    call mesh_init_stage_2(gr, namespace, space, gr%box, gr%stencil)

    POP_SUB(grid_init_stage_1)
  end subroutine grid_init_stage_1

  !--------------------------------------------------------------------
  !> this subroutine allows to create a grid from an existing grid
  subroutine grid_init_from_grid_stage_1(gr, original_gr, namespace, space, box, spacing_prefactor, latt, n_sites, site_position)
    type(grid_t),                      intent(inout) :: gr !< The grid to initialize
    type(grid_t),                      intent(in)    :: original_gr !< The original grid from where to copy the data
    type(namespace_t),                 intent(in)    :: namespace
    type(space_t),                     intent(in)    :: space
    class(box_t), target,    optional, intent(in)    :: box !< the box that contains the grid
    FLOAT,                   optional, intent(in)    :: spacing_prefactor(:)
    type(lattice_vectors_t), optional, intent(in)    :: latt
    integer,                 optional, intent(in)    :: n_sites
    FLOAT,                   optional, intent(in)    :: site_position(:,:)

    type(stencil_t) :: cube
    integer :: enlarge(1:MAX_DIM)

    PUSH_SUB(grid_init_from_grid_stage_1)

    ! First of all, create and associate the box that contains the grid
    if (.not. present(box)) then
      gr%box => box_factory_create(namespace, space)
    else
      gr%box => box
    end if

    !!! Copy all data !!!
    ! Symmetries
    gr%symm = original_gr%symm
    ! Spacing
    if (present(spacing_prefactor)) then
      gr%spacing(:) = spacing_prefactor(:) * original_gr%spacing(:)
      if (any(gr%spacing /= M_ONE)) then
        write(message(1), '(a)') "The two grids have different spacings, it will not be possible to establish a map between them"
        call messages_warning(1)
      end if
    else
      gr%spacing(:) = original_gr%spacing(:)
    end if
    ! Periodic boundaries
    gr%masked_periodic_boundaries = original_gr%masked_periodic_boundaries
    if (gr%masked_periodic_boundaries) gr%periodic_boundary_mask = original_gr%periodic_boundary_mask
    ! Coordinate System
    call initialize_coordinate_system(gr, namespace, space, latt, n_sites, site_position)

    ! initialize derivatives
    call derivatives_init(gr%der, namespace, space, gr%coord_system)
    ! the stencil used to generate the grid is a union of a cube (for
    ! multigrid) and the Laplacian
    call stencil_cube_get_lapl(cube, space%dim, order = 2)
    call stencil_union(space%dim, cube, gr%der%lapl%stencil, gr%stencil)
    call stencil_end(cube)

    enlarge = 0
    enlarge(1:space%dim) = 2
    enlarge = max(enlarge, gr%der%n_ghost)

    call mesh_init_stage_1(gr, namespace, space, gr%box, gr%coord_system, gr%spacing, enlarge)
    call mesh_init_stage_2(gr, namespace, space, gr%box, gr%stencil)

    POP_SUB(grid_init_from_grid_stage_1)
  end subroutine grid_init_from_grid_stage_1

  subroutine initialize_coordinate_system(gr, namespace, space, latt, n_sites, site_position)
    type(grid_t),                      intent(inout) :: gr
    type(namespace_t),                 intent(in)    :: namespace
    type(space_t),                     intent(in)    :: space
    type(lattice_vectors_t), optional, intent(in)    :: latt
    integer,                 optional, intent(in)    :: n_sites
    FLOAT,                   optional, intent(in)    :: site_position(:,:)

    integer :: cv_method

    PUSH_SUB(initialize_coordinate_system)

    !%Variable CurvMethod
    !%Type integer
    !%Default curv_uniform
    !%Section Mesh::Curvilinear
    !%Description
    !% The relevant functions in octopus are represented on a mesh in real space.
    !% This mesh may be an evenly spaced regular rectangular grid (standard mode),
    !% or else an adaptive or curvilinear grid. We have implemented
    !% three kinds of adaptive meshes, although only one is currently working,
    !% the one invented by F. Gygi (<tt>curv_gygi</tt>). The code will stop if any of
    !% the other two is invoked. All are experimental with domain parallelization.
    !%Option curv_affine 1
    !% Regular, uniform rectangular grid.
    !%Option curv_gygi 2
    !% The deformation of the grid is done according to the scheme described by
    !% F. Gygi [F. Gygi and G. Galli, <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)].
    !%Option curv_briggs 3
    !% The deformation of the grid is done according to the scheme described by
    !% Briggs [E.L. Briggs, D.J. Sullivan, and J. Bernholc, <i>Phys. Rev. B</i> <b>54</b> 14362 (1996)]
    !% (NOT WORKING).
    !%Option curv_modine 4
    !% The deformation of the grid is done according to the scheme described by
    !% Modine [N.A. Modine, G. Zumbach and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289 (1997)]
    !% (NOT WORKING).
    !%End
    call parse_variable(namespace, 'CurvMethod', CURV_AFFINE, cv_method)
    if (.not. varinfo_valid_option('CurvMethod', cv_method)) call messages_input_error(namespace, 'CurvMethod')
    if (present(latt)) then
      if (cv_method /= CURV_AFFINE .and. latt%nonorthogonal) then
        call messages_input_error(namespace, 'CurvMethod', 'Curvilinear coordinates with non-orthogonal cells are not implemented')
      end if
    end if
    call messages_print_var_option("CurvMethod", cv_method, namespace=namespace)

    ! FIXME: The other two methods are apparently not working
    if (cv_method > CURV_GYGI) then
      call messages_experimental('Selected curvilinear coordinates method')
    end if

    select case (cv_method)
    case (CURV_BRIGGS)
      gr%coord_system => curv_briggs_t(namespace, space%dim, &
        gr%box%bounding_box_l(1:space%dim), gr%spacing(1:space%dim))
    case (CURV_GYGI)
      if (present(site_position) .and. present(n_sites)) then
        gr%coord_system => curv_gygi_t(namespace, space%dim, n_sites, site_position)
      else
        message(1) = "Option CurvMethod = curv_gygi is not currently implemented without ions present."
        call messages_fatal(1, namespace=namespace)
      end if
    case (CURV_MODINE)
      if (present(site_position) .and. present(n_sites)) then
        gr%coord_system => curv_modine_t(namespace, space%dim, n_sites, site_position, &
          gr%box%bounding_box_l(1:space%dim), gr%spacing(1:space%dim))
      else
        message(1) = "Option CurvMethod = curv_modine is not currently implemented without ions present."
        call messages_fatal(1, namespace=namespace)
      end if
    case (CURV_AFFINE)
      if (present(latt)) then
        if (latt%nonorthogonal) then
          gr%coord_system => affine_coordinates_t(namespace, space%dim, latt%rlattice_primitive)
        else
          gr%coord_system => cartesian_t(namespace, space%dim)
        end if
      else
        gr%coord_system => cartesian_t(namespace, space%dim)
      end if
    end select

    POP_SUB(initialize_coordinate_system)
  end subroutine initialize_coordinate_system


  !-------------------------------------------------------------------
  subroutine grid_init_stage_2(gr, namespace, space, mc, qvector)
    type(grid_t), target, intent(inout) :: gr
    type(namespace_t),    intent(in)    :: namespace
    type(space_t),        intent(in)    :: space
    type(multicomm_t),    intent(in)    :: mc
    FLOAT, optional,      intent(in)    :: qvector(:)

    PUSH_SUB(grid_init_stage_2)

    call mesh_init_stage_3(gr, namespace, space, gr%stencil, mc)

    call nl_operator_global_init(namespace)
    call derivatives_build(gr%der, namespace, space, gr, qvector)

    ! print info concerning the grid
    call grid_write_info(gr, namespace=namespace)

    if(space%dim == 3) then
      call symmetrizer_init(gr%symmetrizer, gr, gr%symm)
    end if

    POP_SUB(grid_init_stage_2)
  end subroutine grid_init_stage_2


  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_t), intent(inout) :: gr

    class(coordinate_system_t), pointer :: coord_system
    class(box_t), pointer :: box

    PUSH_SUB(grid_end)

    call nl_operator_global_end()

    call derivatives_end(gr%der)

    ! We need to take a pointer here, otherwise we run into a gfortran bug.
    coord_system => gr%coord_system
    SAFE_DEALLOCATE_P(coord_system)
    box => gr%box
    SAFE_DEALLOCATE_P(box)

    call mesh_end(gr)

    call stencil_end(gr%stencil)

    call symmetrizer_end(gr%symmetrizer)

    POP_SUB(grid_end)
  end subroutine grid_end


  !-------------------------------------------------------------------
  subroutine grid_write_info(gr, iunit, namespace)
    type(grid_t),                intent(in) :: gr
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(grid_write_info)

    call messages_print_with_emphasis(msg="Grid", iunit=iunit, namespace=namespace)

    message(1) = 'Simulation Box:'
    call messages_info(1, iunit=iunit, namespace=namespace)
    call gr%box%write_info(iunit, namespace)

    message(1) = "Main mesh:"
    call messages_info(1, iunit=iunit, namespace=namespace)
    call mesh_write_info(gr, iunit, namespace)

    if (gr%use_curvilinear) then
      call gr%coord_system%write_info(iunit, namespace)
    end if

    call messages_print_with_emphasis(iunit=iunit, namespace=namespace)

    POP_SUB(grid_write_info)
  end subroutine grid_write_info

#include "undef.F90"
#include "real.F90"
#include "grid_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "grid_inc.F90"

end module grid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
