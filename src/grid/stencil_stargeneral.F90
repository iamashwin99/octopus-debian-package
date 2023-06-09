!! Copyright (C) 2016 U. De Giovannini, H Huebener; 2019 N. Tancogne-Dejean
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

module stencil_stargeneral_oct_m
  use affine_coordinates_oct_m
  use coordinate_system_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use space_oct_m
  use stencil_oct_m

  implicit none

  private
  public ::                     &
    stencil_stargeneral_size_lapl, &
    stencil_stargeneral_extent,    &
    stencil_stargeneral_get_lapl,  &
    stencil_stargeneral_pol_lapl,  &
    stencil_stargeneral_get_arms



contains

  ! ---------------------------------------------------------
  subroutine stencil_stargeneral_get_arms(this, space, coord_system)
    type(stencil_t),            intent(inout) :: this
    type(space_t),              intent(in)    :: space
    class(coordinate_system_t), intent(in)    :: coord_system

    integer :: dim
    FLOAT   :: theta(space%dim), basis_vectors(space%dim, space%dim)
    integer :: arm(space%dim)
    integer :: ii, jj, kk
    FLOAT   :: norm, min_norm(space%dim)

    PUSH_SUB(stencil_stargeneral_get_arms)

    ASSERT(space%dim <= 3)

    dim = space%dim

    this%stargeneral%narms = 0
    this%stargeneral%arms = 0

    if (dim == 1) then
      !we are done
      POP_SUB(stencil_stargeneral_get_arms)
      return
    end if

    select type (coord_system)
    class is (affine_coordinates_t)
      basis_vectors = coord_system%basis%vectors
    class default
      message(1) = "Stencil stargeneral can currently only be used with affine coordinate systems."
      call messages_fatal(1)
    end select

    if (dim == 2) then
      !get the angle between the primitive vectors
      theta(1) = acos(dot_product(basis_vectors(:, 1), basis_vectors(:, 2)))

      if (theta(1) < M_PI*M_HALF) then
        this%stargeneral%narms = this%stargeneral%narms + 1
        arm = (/1,-1/)
        this%stargeneral%arms(this%stargeneral%narms, 1:dim) = arm
      else if (theta(1) > M_PI*M_HALF) then
        this%stargeneral%narms = this%stargeneral%narms + 1
        arm = (/1,+1/)
        this%stargeneral%arms(this%stargeneral%narms, 1:dim) = arm
      end if
      !if theta == pi/2 we do not need additional arms

      ! NB: you have supposed the axis of the 2D system is along z.
      !we are done
      POP_SUB(stencil_stargeneral_get_arms)
      return
    end if

    ! dim>2

    !We first count how many arms we need
    theta(1) = acos(dot_product(basis_vectors(:, 1), basis_vectors(:, 2)))
    if (abs(theta(1) - M_PI*M_HALF) > CNST(1e-8)) this%stargeneral%narms = this%stargeneral%narms + 1

    theta(2) = acos(dot_product(basis_vectors(:, 2), basis_vectors(:, 3)))
    if (abs(theta(2)-M_PI*M_HALF) > CNST(1e-8)) this%stargeneral%narms = this%stargeneral%narms + 1

    theta(3) = acos(dot_product(basis_vectors(:, 3), basis_vectors(:, 1)))
    if (abs(theta(3) - M_PI*M_HALF) > CNST(1e-8)) this%stargeneral%narms = this%stargeneral%narms + 1


    !Cubic cell, no extra arms
    if (this%stargeneral%narms == 0) then
      POP_SUB(stencil_stargeneral_get_arms)
      return
    end if


    !We find the three closest neighbours, to define the diretion of the arms
    min_norm = M_HUGE
    do ii= -1,1
      do jj = -1,1
        do kk = -1,1
          !The nearest points could not be along the primitive axis
          if (ii == 0 .and. jj == 0) cycle
          if (ii == 0 .and. kk == 0) cycle
          if (kk == 0 .and. jj == 0) cycle

          !Nor in a plane formed by orthogonal vectors
          if (abs(theta(1) - M_PI*M_HALF) <= CNST(1e-8) .and. kk == 0) cycle
          if (abs(theta(2) - M_PI*M_HALF) <= CNST(1e-8) .and. ii == 0) cycle
          if (abs(theta(3) - M_PI*M_HALF) <= CNST(1e-8) .and. jj == 0) cycle

          norm = sum((ii*basis_vectors(:, 1) + jj*basis_vectors(:, 2) + kk*basis_vectors(:, 3))**2)

          if (norm < min_norm(1)) then
            if (min_norm(1) < min_norm(2) .and. this%stargeneral%narms == 2) then
              if (min_norm(2) < min_norm(3) .and. this%stargeneral%narms == 3) then
                this%stargeneral%arms(3, 1:3) = this%stargeneral%arms(2, 1:3)
                min_norm(3) =  min_norm(2)
              end if
              this%stargeneral%arms(2, 1:3) = this%stargeneral%arms(1, 1:3)
              min_norm(2) = min_norm(1)
            end if
            min_norm(1) = norm
            this%stargeneral%arms(1, 1:3) = (/ii, jj, kk/)
            cycle
          end if

          if (this%stargeneral%arms(1, 1) == -ii &
            .and. this%stargeneral%arms(1, 2) == -jj &
            .and. this%stargeneral%arms(1, 3) == -kk) cycle

          if (this%stargeneral%narms == 1) cycle

          if (norm < min_norm(2)) then
            if (min_norm(2)< min_norm(3) .and. this%stargeneral%narms == 3) then
              this%stargeneral%arms(3, 1:3) = this%stargeneral%arms(2, 1:3)
              min_norm(3) =  min_norm(2)
            end if
            min_norm(2) = norm
            this%stargeneral%arms(2, 1:3) = (/ii, jj, kk/)
            cycle
          end if

          if (this%stargeneral%narms == 2) cycle

          if (this%stargeneral%arms(2, 1) == -ii &
            .and. this%stargeneral%arms(2, 2) == -jj &
            .and. this%stargeneral%arms(2, 3) == -kk) cycle

          if (norm < min_norm(3)) then
            min_norm(3) = norm
            this%stargeneral%arms(3, 1:3) = (/ii, jj, kk/)
          end if
        end do
      end do
    end do

    POP_SUB(stencil_stargeneral_get_arms)
  end subroutine stencil_stargeneral_get_arms


  ! ---------------------------------------------------------
  integer function stencil_stargeneral_size_lapl(this, dim, order) result(n)
    type(stencil_t),     intent(inout) :: this
    integer, intent(in) :: dim
    integer, intent(in) :: order

    PUSH_SUB(stencil_stargeneral_size_lapl)

    !normal star
    n = 2*dim*order + 1

    ! star general
    n = n + 2 * order * this%stargeneral%narms


    POP_SUB(stencil_stargeneral_size_lapl)
  end function stencil_stargeneral_size_lapl


  ! ---------------------------------------------------------
  !> Returns maximum extension of the stencil in spatial direction
  !! dir = 1, 2, 3 for a given discretization order.
  integer function stencil_stargeneral_extent(dir, order)
    integer, intent(in) :: dir
    integer, intent(in) :: order

    integer :: extent

    PUSH_SUB(stencil_stargeneral_extent)

    extent = 0
    if (dir >= 1 .or. dir <= 3) then
      if (order <= 2) then
        extent = 2
      else
        extent = order
      end if
    end if
    stencil_stargeneral_extent = extent

    POP_SUB(stencil_stargeneral_extent)
  end function stencil_stargeneral_extent



  ! ---------------------------------------------------------
  subroutine stencil_stargeneral_get_lapl(this, dim, order)
    type(stencil_t), intent(inout) :: this
    integer,         intent(in)    :: dim
    integer,         intent(in)    :: order

    integer :: i, j, n
    logical :: got_center

    PUSH_SUB(stencil_stargeneral_get_lapl)

    call stencil_allocate(this, stencil_stargeneral_size_lapl(this, dim, order))

    n = 1
    select case (dim)
    case (1)
      n = 1
      do i = 1, dim
        do j = -order, order
          if (j == 0) cycle
          n = n + 1
          this%points(i, n) = j
        end do
      end do
    case (2)
      n = 1
      do i = 1, dim
        do j = -order, order
          if (j == 0) cycle
          n = n + 1
          this%points(i, n) = j
        end do
      end do

      do j = -order, order
        if (j == 0) cycle
        do i = 1, this%stargeneral%narms
          n = n + 1
          this%points(1:2, n) = this%stargeneral%arms(i, 1:2)*j
        end do
      end do

    case (3)
      got_center = .false.

      n = 0
      do i = 1, dim
        do j = -order, order

          ! count center only once
          if (j == 0) then
            if (got_center) then
              cycle
            else
              got_center = .true.
            end if

          end if
          n = n + 1
          this%points(i, n) = j
        end do
      end do

      do j = -order, order
        if (j == 0) cycle
        do i = 1, this%stargeneral%narms
          n = n + 1
          this%points(1:3, n) = this%stargeneral%arms(i, 1:3)*j
        end do
      end do

    end select

    call stencil_init_center(this)

    POP_SUB(stencil_stargeneral_get_lapl)
  end subroutine stencil_stargeneral_get_lapl




  ! ---------------------------------------------------------
  subroutine stencil_stargeneral_pol_lapl(this, dim, order, pol)
    type(stencil_t), intent(in) :: this
    integer, intent(in)          :: dim
    integer, intent(in)          :: order
    integer, intent(out)         :: pol(:,:) !< pol(dim, order)

    integer :: i, j, n, j1, j2

    PUSH_SUB(stencil_stargeneral_pol_lapl)

    n = 1
    select case (dim)
    case (1)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do
    case (2)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do

      do j = 1, 2*order
        do i = 1, this%stargeneral%narms
          n = n + 1
          if (sum(this%stargeneral%arms(i,1:dim)) == 0) then
            pol(1:2, n) = (/j,1/)
          else
            pol(1:2, n) = (/1,j/)
          end if
        end do
      end do

    case (3)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do

      do i = 1, this%stargeneral%narms

        if (this%stargeneral%arms(i, 1) == 0) then
          ! sum(this%stargeneral%arms(i, 1:dim)) == 0 just checks whether we have a -1 in the arm vector or not
          if (sum(this%stargeneral%arms(i, 1:dim)) == 0) then
            do j1 = 1, 2*order
              do j2 = 1, 2*order - j1
                n = n + 1
                pol(1:3, n) = (/0, j1, j2/)
              end do
            end do
          else
            do j1 = 1, 2*order
              do j2 = 1, 2*order - j1
                n = n + 1
                pol(1:3, n) = (/0, j2, j1/)
              end do
            end do
          end if
        end if

        if (this%stargeneral%arms(i, 2) == 0) then
          ! sum(this%stargeneral%arms(i, 1:dim)) == 0 just checks whether we have a -1 in the arm vector or not
          if (sum(this%stargeneral%arms(i, 1:dim)) == 0) then
            do j1 = 1, 2*order
              do j2 = 1, 2*order - j1
                n = n + 1
                pol(1:3, n) = (/j2, 0, j1/)
              end do
            end do
          else
            do j1 = 1, 2*order
              do j2 = 1, 2*order - j1
                n = n + 1
                pol(1:3, n) = (/j1, 0, j2/)
              end do
            end do
          end if
        end if

        if (this%stargeneral%arms(i, 3) == 0) then
          ! sum(this%stargeneral%arms(i, 1:dim)) == 0 just checks whether we have a -1 in the arm vector or not
          if (sum(this%stargeneral%arms(i,1:dim)) == 0) then
            do j1 = 1, 2*order
              do j2 = 1, 2*order - j1
                n = n + 1
                pol(1:3, n) = (/j1, j2, 0/)
              end do
            end do
          else
            do j1 = 1, 2*order
              do j2 = 1, 2*order - j1
                n = n + 1
                pol(1:3, n) = (/j2, j1, 0/)
              end do
            end do
          end if
        end if
      end do ! end loop on number of arms

!       !FCC
!       do j = 1, 2*order
!         n = n + 1
!         pol(1:3, n) = (/j,1,0/)
!         n = n + 1
!         pol(1:3, n) = (/1,0,j/)
!         n = n + 1
!         pol(1:3, n) = (/0,j,1/)
!       end do

!       !HEX
!       do j = 1, 2*order
!         n = n + 1
!         pol(1:3, n) = (/j, 1, 0/)
!       end do


    end select

    POP_SUB(stencil_stargeneral_pol_lapl)
  end subroutine stencil_stargeneral_pol_lapl


end module stencil_stargeneral_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
