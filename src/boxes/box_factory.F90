!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2013,2021 M. Oliveira
!! Copyright (C) 2021 K. Lively, A. Obzhirov, I. Albar
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

!> for the moment, just have a subroutine which returns a box
!> could be generalized into a box_factory_t class
module box_factory_oct_m
  use box_oct_m
  use box_cgal_oct_m
  use box_cylinder_oct_m
  use box_image_oct_m
  use box_minimum_oct_m
  use box_parallelepiped_oct_m
  use box_sphere_oct_m
  use box_user_defined_oct_m
  use global_oct_m
  use lattice_vectors_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public  ::        &
    box_factory_create

  integer, parameter, public :: &
    SPHERE         = 1,         &
    CYLINDER       = 2,         &
    MINIMUM        = 3,         &
    PARALLELEPIPED = 4,         &
    BOX_IMAGE      = 5,         &
    BOX_CGAL       = 6,         &
    BOX_USDEF      = 77

contains

  !> initialize a box of any type
  !> if the system is periodic, the lattice vectors should be passed
  !> if ions are present, ions should be passed
  function box_factory_create(namespace, space, latt, n_sites, site_position) result(box)
    type(namespace_t),                 intent(in)    :: namespace
    type(space_t),                     intent(in)    :: space
    type(lattice_vectors_t), optional, intent(in)    :: latt
    integer,                 optional, intent(in)    :: n_sites
    FLOAT,                   optional, intent(in)    :: site_position(:,:)
    class(box_t), pointer ::  box

    type(block_t) :: blk
    integer       :: default_boxshape, idir, box_shape
    FLOAT         :: center(space%dim), axes(space%dim, space%dim), rsize, xsize, lsize(space%dim)
    character(len=1024) :: filename
    character(len=1024) :: user_def

    PUSH_SUB(box_factory_create)

    ! We must have boht n_sites and site_position or none.
    ASSERT(present(n_sites) .eqv. present(site_position))

    !%Variable BoxShape
    !%Type integer
    !%Section Mesh::Simulation Box
    !%Description
    !% This variable decides the shape of the simulation box.
    !% The default is <tt>minimum</tt> for finite systems and <tt>parallelepiped</tt> for periodic systems.
    !% Note that some incompatibilities apply:
    !% <ul><li>Spherical or minimum mesh is not allowed for periodic systems.
    !% <li>Cylindrical mesh is not allowed for systems that are periodic in more than one dimension.
    !% <li><tt>box_image</tt> is only allowed in 2D.</ul>
    !%Option sphere 1
    !% The simulation box will be a sphere of radius <tt>Radius</tt>. (In 2D, this is a circle.)
    !%Option cylinder 2
    !% The simulation box will be a cylinder with radius <tt>Radius</tt> and height (in the <i>x</i>-direction)
    !% of 2 <tt>Xlength</tt>.
    !%Option minimum 3
    !% The simulation box will be constructed by adding spheres created around each
    !% atom (or user-defined potential), of radius <tt>Radius</tt>.
    !%Option parallelepiped 4
    !% The simulation box will be a parallelepiped whose dimensions are taken from
    !% the variable <tt>Lsize</tt>.
    !%Option box_image 5
    !% The simulation box will be defined through an image, specified with <tt>BoxShapeImage</tt>.
    !% White (RGB = 255,255,255) means that the point
    !% is contained in the simulation box, while any other color means that the point is out.
    !% The image will be scaled to fit <tt>Lsize</tt>, while its resolution will define the default <tt>Spacing</tt>.
    !% The actual box may be slightly larger than <tt>Lsize</tt> to ensure one grid point = one pixel for
    !% default <tt>Spacing</tt>.
    !%Option box_cgal 6
    !% The simulation box will be defined by a file read using the CGAL library.
    !% The file name needs to be specified with <tt>BoxCgalFile</tt>.
    !% <tt>Lsize</tt> needs to be large enough to contain the shape defined in the file.
    !%Option user_defined 77
    !% The shape of the simulation box will be read from the variable <tt>BoxShapeUsDef</tt>.
    !%End

    if (space%is_periodic()) then
      default_boxshape = PARALLELEPIPED
      if(.not. present(latt)) then
        message(1) = "Periodic system box is trying to be generated without lattice vectors given."
        call messages_fatal(1, namespace=namespace)
      endif
    else
      if (present(site_position)) then
        default_boxshape = MINIMUM !< only well defined if we have a list of sites and their positions
      else
        default_boxshape = PARALLELEPIPED
      endif
    end if
    call parse_variable(namespace, 'BoxShape', default_boxshape, box_shape)
    if (.not. varinfo_valid_option('BoxShape', box_shape)) call messages_input_error(namespace, 'BoxShape')
    select case (box_shape)
    case (SPHERE, MINIMUM, BOX_USDEF)
      if (space%dim > 1 .and. space%is_periodic()) call messages_input_error(namespace, 'BoxShape')
    case (CYLINDER)
      if (space%dim == 2) then
        message(1) = "BoxShape = cylinder is not meaningful in 2D. Use sphere if you want a circle."
        call messages_fatal(1, namespace=namespace)
      end if
      if (space%periodic_dim > 1) call messages_input_error(namespace, 'BoxShape')
    end select

    ! ignore box_shape in 1D
    if (space%dim == 1 .and. box_shape /= PARALLELEPIPED) then
      box_shape = SPHERE
    end if

    ! Cannot use images in 1D or 3D
    if (space%dim /= 2 .and. box_shape == BOX_IMAGE) call messages_input_error(namespace, 'BoxShape')

    if (space%dim > 3 .and. box_shape /= PARALLELEPIPED) then
      message(1) = "For more than 3 dimensions, you can only use the parallelepiped box."
      call messages_fatal(1, namespace=namespace)
      ! FIXME: why not a hypersphere as another option?
    end if

    !%Variable Radius
    !%Type float
    !%Section Mesh::Simulation Box
    !%Description
    !% Defines the radius for <tt>BoxShape</tt> = <tt>sphere</tt>,
    !% <tt>cylinder</tt>, or <tt>minimum</tt>. Must be a positive
    !% number.
    !%End
    if (box_shape == SPHERE .or. box_shape == CYLINDER .or. box_shape == MINIMUM) then
      if (.not. parse_is_defined(namespace, 'Radius')) then
        message(1) = "BoxShape = sphere, cylinder and minimum require the Radius variable to be defined."
        call messages_fatal(1,namespace=namespace)
      else
        call parse_variable(namespace, 'Radius', -M_ONE, rsize, units_inp%length)
      end if
      if (rsize < M_ZERO) call messages_input_error(namespace, 'Radius')
    end if

    if (box_shape == MINIMUM .and. .not. present(site_position)) then
      !> We should only be here if ions are present, something has gone wrong
      message(1) = "BoxShape = minimum only makes sense when the calculation includes atomic sites."
      call messages_fatal(1, namespace=namespace)
    end if

    if (box_shape == CYLINDER) then
      !%Variable Xlength
      !%Default <tt>Radius</tt>
      !%Type float
      !%Section Mesh::Simulation Box
      !%Description
      !% If <tt>BoxShape</tt> is <tt>cylinder</tt>, the total length of the cylinder is twice <tt>Xlength</tt>.
      !% Note that when PeriodicDimensions = 1, then the length of the cylinder is determined from the lattice vectors.
      !%End
      if (space%is_periodic()) then
        xsize = norm2(latt%rlattice(1:space%periodic_dim, 1))/M_TWO
      else
        call parse_variable(namespace, 'Xlength', rsize, xsize, units_inp%length)
      end if
    end if

    lsize = M_ZERO
    if (box_shape == PARALLELEPIPED .or. box_shape == BOX_IMAGE .or. &
      box_shape == BOX_USDEF .or. box_shape == BOX_CGAL) then

      !%Variable Lsize
      !%Type block
      !%Section Mesh::Simulation Box
      !%Description
      !% If <tt>BoxShape</tt> is <tt>parallelepiped</tt>, <tt>box_image</tt>,
      !% or <tt>user_defined</tt>, this is a block of the form:
      !%
      !% <tt>%Lsize
      !% <br>&nbsp;&nbsp;sizex | sizey | sizez | ...
      !% <br>%</tt>
      !%
      !% where the <tt>size*</tt> are half the lengths of the box in each direction.
      !%
      !% The number of columns must match the dimensionality of the
      !% calculation. If you want a cube you can also set <tt>Lsize</tt> as a
      !% single variable.
      !%End
      if (present(latt)) then
        !> lattice should only be passed if periodic
        ! lsize along the periodic dimensions must always be set from the norm of the lattice vectors
        do idir = 1, space%periodic_dim
          lsize(idir) = norm2(latt%rlattice(1:space%dim, idir))/M_TWO
        end do
      endif

      if (space%is_periodic()) then
        ! For mixed-periodicity, lsize along the non-periodic dimensions is
        ! by default set from the lattice parameters (this can still be
        ! overriden by setting Lsize, see bellow).
        do idir = space%periodic_dim + 1, space%dim
          lsize(idir) = norm2(latt%rlattice(1:space%dim, idir))/M_TWO
        end do
      else
        ! Lsize must be set for finite systems, as in that case we do not have the lattice parameters
        if (.not. parse_is_defined(namespace, 'Lsize')) then
          call messages_input_error(namespace, 'Lsize', 'Lsize is required for finite systems')
        end if
      end if

      ! Note that for cases with mixed-periodicidy, the user still has the
      ! option to set Lsize to override the size of the box along the
      ! non-periodic dimensions given by the lattice parameters. This
      ! requires the user to also set Lsize for the periodic dimensions,
      ! which at the moment must match exactly the corresponding values
      ! given by the lattice vectors.
      if (parse_block(namespace, 'Lsize', blk) == 0) then
        ! Lsize is specified as a block
        if (parse_block_cols(blk, 0) < space%dim) then
          call messages_input_error(namespace, 'Lsize')
        end if

        do idir = 1, space%dim
          call parse_block_float(blk, 0, idir - 1, lsize(idir), units_inp%length)
        end do
        call parse_block_end(blk)

      else if (parse_is_defined(namespace, 'Lsize')) then
        ! Lsize is specified as a scalar
        call parse_variable(namespace, 'Lsize', -M_ONE, lsize(1), units_inp%length)
        if (abs(lsize(1) + M_ONE) <= M_EPSILON) then
          call messages_input_error(namespace, 'Lsize')
        end if
        do idir = 2, space%dim
          lsize(idir) = lsize(1)
        end do

      end if

      ! Check that lsize is consistent with the lattice vectors along the periodic dimensions
      do idir = 1, space%periodic_dim
        if (abs(M_TWO*lsize(idir) - norm2(latt%rlattice(1:space%dim, idir))) > M_EPSILON) then
          call messages_input_error(namespace, 'Lsize', &
            'Lsize must be exactly half the length of the lattice vectors along periodic dimensions')
        end if
      end do
    end if

    ! read in image for box_image
    if (box_shape == BOX_IMAGE) then

      !%Variable BoxShapeImage
      !%Type string
      !%Section Mesh::Simulation Box
      !%Description
      !% Name of the file that contains the image that defines the simulation box
      !% when <tt>BoxShape = box_image</tt>. No default. Will search in current
      !% directory and <tt>OCTOPUS-HOME/share/</tt>.
      !%End
#if defined(HAVE_GDLIB)
      call parse_variable(namespace, 'BoxShapeImage', '', filename)
      if (trim(filename) == "") then
        message(1) = "Must specify BoxShapeImage if BoxShape = box_image."
        call messages_fatal(1, namespace=namespace)
      end if
#else
      message(1) = "To use 'BoxShape = box_image', you have to compile Octopus"
      message(2) = "with GD library support."
      call messages_fatal(2, namespace=namespace)
#endif
    end if

    ! read in box shape for user-defined boxes
    if (box_shape == BOX_USDEF) then

      !%Variable BoxShapeUsDef
      !%Type string
      !%Section Mesh::Simulation Box
      !%Description
      !% Boolean expression that defines the interior of the simulation box. For example,
      !% <tt>BoxShapeUsDef = "(sqrt(x^2+y^2) <= 4) && z>-2 && z<2"</tt> defines a cylinder
      !% with axis parallel to the <i>z</i>-axis.
      !%End

      call parse_variable(namespace, 'BoxShapeUsDef', 'x^2+y^2+z^2 < 4', user_def)
    end if

    ! get filename for cgal boxes
    if (box_shape == BOX_CGAL) then
#ifndef HAVE_CGAL
      message(1) = "To use 'BoxShape = box_cgal', you have to compile Octopus"
      message(2) = "with CGAL library support."
      call messages_fatal(2, namespace=namespace)
#endif
      !%Variable BoxCgalFile
      !%Type string
      !%Section Mesh::Simulation Box
      !%Description
      !% Filename to be read in by the cgal library. It should describe a shape that
      !% is used for the simulation box
      !%End
      if (.not. parse_is_defined(namespace, 'BoxCgalFile')) then
        message(1) = "Must specify BoxCgalFile if BoxShape = box_cgal."
        call messages_fatal(1, namespace=namespace)
      end if
      call parse_variable(namespace, 'BoxCgalFile', '', filename)
    end if

    call messages_obsolete_variable(namespace, 'BoxOffset')

    !%Variable BoxCenter
    !%Type float
    !%Section Mesh::Simulation Box
    !%Description
    !% This block defines the coordinate center of the simulation box
    !%End
    if(parse_block(namespace, 'BoxCenter', blk) == 0) then
      if(parse_block_cols(blk, 0) < space%dim) call messages_input_error(namespace, 'BoxCenter')
      do idir = 1, space%dim
        call parse_block_float(blk, 0, idir - 1, center(idir), units_inp%length)
      end do
      call parse_block_end(blk)
    else
      center = M_ZERO  ! In case no BoxCenter is explicitly defined in the input file the default is the origin
    end if

    if (present(latt)) then
      ! Use the lattice vectors
      axes = latt%rlattice
    else
      ! Use Cartesian axes (instead, we could read this from the input file here)
      axes = diagonal_matrix(space%dim, M_ONE)
    end if

    ! Now generate the box
    select case (box_shape)
    case (SPHERE)
      box => box_sphere_t(space%dim, center, rsize, namespace)
    case (CYLINDER)
      box => box_cylinder_t(space%dim, center, axes, rsize, M_TWO * xsize, namespace, periodic_boundaries=space%is_periodic())
    case (PARALLELEPIPED)
      box => box_parallelepiped_t(space%dim, center, axes, M_TWO * lsize(1:space%dim), namespace, &
        n_periodic_boundaries=space%periodic_dim)
    case (BOX_USDEF)
      box => box_user_defined_t(space%dim, center, axes, user_def, M_TWO*lsize(1:space%dim), namespace)
    case (BOX_CGAL)
      box => box_cgal_t(space%dim, center, filename, M_TWO*lsize(1:space%dim), namespace)
    case (MINIMUM)
      box => box_minimum_t(space%dim, rsize, n_sites, site_position, namespace)
    case (BOX_IMAGE)
      box => box_image_t(center, axes, lsize, filename, space%periodic_dim, namespace)
    end select

    POP_SUB(box_factory_create)
  end function box_factory_create

end module box_factory_oct_m
