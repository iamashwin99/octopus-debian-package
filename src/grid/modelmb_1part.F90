!! Copyright (C) 2009 N. Helbig and M. Verstraete
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

!>  general module for modelmb particles (eg 4 electrons in 1D equiv to
!!  1 in 4D). Also calculate different densities on request.
module modelmb_1part_oct_m
  use debug_oct_m
  use global_oct_m
  use hypercube_oct_m
  use messages_oct_m
  use mesh_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::             &
    modelmb_1part_init, &
    modelmb_1part_end,  &
    modelmb_1part_t

  !>  container type for the position and dimensions for 1 particle (out of MAX_DIM/dims)
  type modelmb_1part_t
    ! All components are public by default
    integer :: ndim1part
    integer, private :: npt_part
    integer :: npt
    FLOAT   :: vol_elem_1part
    FLOAT,   allocatable :: origin(:)
    integer, allocatable :: enlarge_1part(:)
    integer, allocatable :: nr_1part(:,:)
    integer, allocatable, private :: ll(:)
    FLOAT,   allocatable :: h_1part(:)
    type(hypercube_t) :: hypercube_1part
  end type modelmb_1part_t

contains

  subroutine modelmb_1part_init(this, mesh, ikeeppart, ndim1part)
    type(modelmb_1part_t), intent(out) :: this
    type(mesh_t),          intent(in)  :: mesh
    integer,               intent(in)  :: ikeeppart
    integer,               intent(in)  :: ndim1part

    integer :: idir, irealdir

    PUSH_SUB(modelmb_1part_init)

    this%ndim1part = ndim1part

!   get full size of arrays for 1 particle only in ndim_modelmb dimensions
    this%npt_part = 1
    do idir = 1, ndim1part
      this%npt_part = this%npt_part*(mesh%idx%nr(2,(ikeeppart - 1)*ndim1part + idir) &
        - mesh%idx%nr(1,(ikeeppart - 1)*ndim1part + idir) + 1)
    end do

!real bounds for indices in
    this%npt = 1
    SAFE_ALLOCATE(this%ll(1:ndim1part))
    do idir = 1, ndim1part
      this%ll(idir) = mesh%idx%ll((ikeeppart-1)*ndim1part+idir)
      this%npt = this%npt*this%ll(idir)
    end do

!   volume element for the chosen particle
    SAFE_ALLOCATE(this%h_1part(1:ndim1part))
    this%vol_elem_1part = M_ONE
    do idir = 1,ndim1part
      irealdir = (ikeeppart-1)*ndim1part + idir
      this%vol_elem_1part = this%vol_elem_1part*mesh%spacing(irealdir)
      this%h_1part(idir) = mesh%spacing(irealdir)
    end do

!   store start and end positions for the relevant dimensions for this particle
    SAFE_ALLOCATE(this%nr_1part(1:2, 1:ndim1part))
    this%nr_1part(:,:) = mesh%idx%nr(:,(ikeeppart-1)*ndim1part+1:ikeeppart*ndim1part)

!   initialize a hypercube for just this particle
!   NB: hypercube_* presume that enlarge is the same for all dimensions!
    SAFE_ALLOCATE(this%enlarge_1part(1:ndim1part))
    this%enlarge_1part = mesh%idx%enlarge((ikeeppart-1)*ndim1part+1:ikeeppart*ndim1part)
    call hypercube_init(this%hypercube_1part, ndim1part, this%nr_1part, this%enlarge_1part(1))

    ! not always the real origin if the box is shifted, no?
    !  which happens to be my case...
    !  only important for printout, so it is ok
    SAFE_ALLOCATE(this%origin(1:ndim1part))
    do idir = 1,ndim1part
      irealdir = (ikeeppart-1)*ndim1part + idir
      !origin(idir) = (npoints(irealdir)/2)*gr%mesh%spacing(irealdir)
      this%origin(idir) = M_ZERO
    end do

    POP_SUB(modelmb_1part_init)
  end subroutine modelmb_1part_init


  subroutine modelmb_1part_end(this)
    type(modelmb_1part_t), intent(inout) :: this

    PUSH_SUB(modelmb_1part_end)

    SAFE_DEALLOCATE_A(this%origin)
    SAFE_DEALLOCATE_A(this%enlarge_1part)
    SAFE_DEALLOCATE_A(this%nr_1part)
    SAFE_DEALLOCATE_A(this%ll)
    SAFE_DEALLOCATE_A(this%h_1part)

    call hypercube_end(this%hypercube_1part)

    POP_SUB(modelmb_1part_end)
  end subroutine modelmb_1part_end


end module modelmb_1part_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
