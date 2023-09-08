!! Copyright (C) 2022 F. Bonafé
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

module current_to_mxll_field_oct_m
  use clock_oct_m
  use debug_oct_m
  use field_transfer_oct_m
  use global_oct_m
  use grid_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use lattice_vectors_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use regridding_oct_m
  use space_oct_m
  use states_mxll_oct_m
  use submesh_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  integer, parameter :: COS2 = 1
  public  ::                    &
    current_to_mxll_field_t

  type, extends(field_transfer_t) :: current_to_mxll_field_t
    private
    logical, public                :: grid_based_partner = .true.

    !> For partners that are point-wise charges
    integer, public :: partner_np = 0 !< number of particles in the partner system
    FLOAT, allocatable, public :: partner_charge(:) !< array storing a copy of the masses of the partner particles
    FLOAT, allocatable, public :: partner_pos(:,:) !< array storing a copy of the positions of the partner particles
    FLOAT, allocatable, public :: partner_vel(:,:) !< array storing a copy of the positions of the partner particles

    integer :: reg_type !< regularization function type
    FLOAT   :: reg_width
    FLOAT, allocatable :: reg(:)

    type(space_t), pointer :: system_space => null()
    type(lattice_vectors_t), pointer :: system_latt => null()


  contains
    procedure :: init_space_latt => current_to_mxll_field_init_space_latt
    procedure :: calculate => current_to_mxll_field_calculate
    procedure :: calculate_energy => current_to_mxll_field_calculate_energy
    final :: current_to_mxll_field_finalize
  end type current_to_mxll_field_t


  interface current_to_mxll_field_t
    module procedure current_to_mxll_field_constructor
  end interface current_to_mxll_field_t

contains

  ! ---------------------------------------------------------
  function current_to_mxll_field_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(current_to_mxll_field_t),               pointer       :: this

    PUSH_SUB(current_to_mxll_field_constructor)

    SAFE_ALLOCATE(this)

    this%label = "current_to_mxll_field"
    this%partner => partner

    this%n_system_quantities = 0
    SAFE_ALLOCATE(this%system_quantities(1:this%n_system_quantities))
    nullify(this%system_gr)

    this%n_partner_quantities = 1
    SAFE_ALLOCATE(this%partner_quantities(1:this%n_partner_quantities))
    ! TODO: For point-wise systems, this should be different, as we
    ! need there the charge and velocity to do the charge regularization.
    ! Once each partner has a basis class, we can do a derived type
    ! TODO: look for user input option of the regularization function for point partner

    !%Variable RegularizationFunction
    !%Type integer
    !%Default COS2
    !%Section Maxwell
    !%Description
    !% The current arising from charged point particles must be mapped onto the Maxwell
    !% propagation grid. This requires a smearing or regularization function $\phi(\mathbf{r})$ attached to
    !% each particle position $\mathbf{r}_i$ with user defined cutoff width, $\sigma$
    !%Option COS2 1
    !% $\phi(r)=\text{cos}^2(\frac{\pi}{2}\frac{|\mathbf{r}-\mathbf{r}_i|}{\sigma})$
    !% if $|\mahtbf{r}-\mathbf{r}_i|<\sigma$, and 0 otherwise.
    !%End

    call parse_variable(partner%namespace, 'RegularizationFunction', COS2, this%reg_type)
    if (.not. varinfo_valid_option('RegularizationFunction', this%reg_type)) &
      call messages_input_error(partner%namespace, 'RegularizationFunction')

    !%Variable RegularizationFunctionWidth
    !%Type float
    !%Default 2
    !%Section Maxwell
    !%Description
    !% The current arising from charged point particles must be mapped onto the Maxwell
    !% propagation grid. This requires a smearing or regularization function $\phi(\mathbf{r})$ attached to
    !% each particle position $\mathbf{r}_i$ with user defined cutoff width, $\sigma$.
    !% Default 2 bohrradii
    !%End

    call parse_variable(partner%namespace, 'RegularizationFunctionWidth', CNST(2.0), &
      this%reg_width, units_inp%length)

    this%partner_quantities(1) = CURRENT

    this%intra_interaction = .false.

    POP_SUB(current_to_mxll_field_constructor)
  end function current_to_mxll_field_constructor


  subroutine current_to_mxll_field_init_space_latt(this, space, latt)
    class(current_to_mxll_field_t), intent(inout) :: this
    type(space_t), target, intent(in)             :: space
    type(lattice_vectors_t), target, intent(in)   :: latt

    PUSH_SUB(current_to_mxll_field_init_space_latt)

    this%system_space => space
    this%system_latt => latt

    POP_SUB(current_to_mxll_field_init_space_latt)
  end subroutine current_to_mxll_field_init_space_latt

  ! ---------------------------------------------------------
  subroutine current_to_mxll_field_finalize(this)
    type(current_to_mxll_field_t), intent(inout) :: this

    PUSH_SUB(current_to_mxll_field_finalize)

    call this%end()
    nullify(this%system_space)
    nullify(this%system_latt)

    POP_SUB(current_to_mxll_field_finalize)
  end subroutine current_to_mxll_field_finalize

  ! ---------------------------------------------------------
  subroutine current_to_mxll_field_calculate(this)
    class(current_to_mxll_field_t), intent(inout) :: this

    integer :: part_ind, i_dim, ip
    type(submesh_t) :: submesh
    FLOAT :: norm
    CMPLX, allocatable :: rs_current(:, :)

    type(profile_t), save :: prof

    PUSH_SUB(current_to_mxll_field_calculate)

    call profiling_in(prof,"CURRENT_MXLL_FIELD_CALC")

    if(this%grid_based_partner) then
      call this%regridding%do_transfer(this%system_field, this%partner_field)
    else
      this%system_field = M_ZERO
      do part_ind = 1, this%partner_np
        call submesh_init(submesh, this%system_space, this%system_gr, &
          this%system_latt, this%partner_pos(:,part_ind), this%reg_width)

        SAFE_ALLOCATE(this%reg(1:submesh%np))
        if (this%reg_type == COS2) then
          do ip = 1, submesh%np
            this%reg(ip) = cos( (M_PI/M_TWO) * (submesh%r(ip)/this%reg_width) )**2
          end do
        end if

        norm = dsm_integrate(this%system_gr, submesh, this%reg)

        do i_dim = 1, this%system_space%dim
          call submesh_add_to_mesh(submesh, this%reg, this%system_field(:,i_dim), &
            this%partner_charge(part_ind)*this%partner_vel(i_dim,part_ind)/norm)
        end do
        SAFE_DEALLOCATE_A(this%reg)

        call submesh_end(submesh)
      end do

    end if
    ! add the RS current to the interpolation
    SAFE_ALLOCATE(rs_current(1:this%system_gr%np, 1:this%ndim))
    call build_rs_current_state(this%system_field, this%system_gr, rs_current)
    ASSERT(this%interpolation_initialized)
    call this%interpolation%add_time(this%partner%quantities(this%partner_quantities(1))%clock%time(), &
      rs_current)
    SAFE_DEALLOCATE_A(rs_current)

    call profiling_out(prof)

    POP_SUB(current_to_mxll_field_calculate)
  end subroutine current_to_mxll_field_calculate

  ! ---------------------------------------------------------
  subroutine current_to_mxll_field_calculate_energy(this)
    class(current_to_mxll_field_t),    intent(inout) :: this

    PUSH_SUB(current_to_mxll_field_calculate_energy)

    ! interaction energy is zero, since it is only re-gridding the quantities of one system
    ! on the mesh of the other
    this%energy = M_ZERO

    POP_SUB(current_to_mxll_field_calculate_energy)
  end subroutine current_to_mxll_field_calculate_energy

end module current_to_mxll_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
