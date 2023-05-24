!! Copyright (C) 2020 F. Bonafé, H. Appel, R. Jestädt
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

module linear_medium_oct_m
  use algorithm_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use cgal_polyhedra_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use iso_c_binding
  use io_oct_m
  use linear_medium_to_em_field_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use comm_oct_m
  use grid_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use regridding_oct_m
  use space_oct_m
  use system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                    &
    linear_medium_t,           &
    linear_medium_init,        &
    get_medium_box_points_map, &
    medium_box_init

  integer, parameter :: &
    MEDIUM_PARALLELEPIPED = 1,         &
    MEDIUM_BOX_FILE       = 2

  type, extends(system_t) :: linear_medium_t
    FLOAT                 :: ep_factor !< permitivity before applying edge profile
    FLOAT                 :: mu_factor !< permeability before applying edge profile
    FLOAT                 :: sigma_e_factor !< electric conductivy before applying edge profile
    FLOAT                 :: sigma_m_factor !< magnetic conductivity before applying edge4 profile
    integer               :: edge_profile  !< edge shape profile (smooth or steep)
    type(output_t)        :: outp
    type(grid_t)          :: gr    !< the mesh
    type(multicomm_t)     :: mc    !< index and domain communicators
    type(single_medium_box_t) :: medium_box

  contains
    procedure :: init_interaction => linear_medium_init_interaction
    procedure :: init_interaction_as_partner => linear_medium_init_interaction_as_partner
    procedure :: initial_conditions => linear_medium_initial_conditions
    procedure :: do_td_operation => linear_medium_do_td
    procedure :: is_tolerance_reached => linear_medium_is_tolerance_reached
    procedure :: update_quantity => linear_medium_update_quantity
    procedure :: update_exposed_quantity => linear_medium_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => linear_medium_copy_quantities_to_interaction
    procedure :: restart_write_data => linear_medium_restart_write_data
    procedure :: restart_read_data => linear_medium_restart_read_data
    procedure :: update_kinetic_energy => linear_medium_update_kinetic_energy
    final :: linear_medium_finalize
  end type linear_medium_t

  interface linear_medium_t
    procedure linear_medium_constructor
  end interface linear_medium_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function linear_medium_constructor(namespace) result(sys)
    class(linear_medium_t), pointer    :: sys
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(linear_medium_constructor)

    SAFE_ALLOCATE(sys)

    call linear_medium_init(sys, namespace)

    POP_SUB(linear_medium_constructor)
  end function linear_medium_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine linear_medium_init(this, namespace)
    class(linear_medium_t), target, intent(inout) :: this
    type(namespace_t),            intent(in)      :: namespace

    integer :: nlines, ncols,n_points, n_global_points
    integer(i8) :: index_range(4)
    integer, allocatable :: tmp(:)
    type(block_t) :: blk
    type(profile_t), save :: prof

    PUSH_SUB(linear_medium_init)

    call profiling_in(prof, 'MEDIUM_BOX_INIT')

    this%namespace = namespace
    call space_init(this%space, this%namespace)
    if (this%space%is_periodic()) then
      call messages_not_implemented('Linear medium for periodic systems', namespace=namespace)
    end if
    call grid_init_stage_1(this%gr, this%namespace, this%space)
    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%mesh%np_global  ! Number of points in mesh
    index_range(2) = 1                      ! Number of states
    index_range(3) = 1                      ! Number of k-points
    index_range(4) = 100000                 ! Some large number

    ! create index and domain communicators
    call multicomm_init(this%mc, this%namespace, mpi_world, calc_mode_par_parallel_mask(), &
         &calc_mode_par_default_parallel_mask(), mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))
    call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc)

    call medium_box_init(this%medium_box, this%namespace)

    !%Variable LinearMediumProperties
    !%Type block
    !%Section Maxwell
    !%Description
    !% Defines electromagnetic parameters for a linear medium box.
    !%
    !% Example:
    !%
    !% <tt>%LinearMediumProperties
    !% <br>&nbsp;&nbsp;   epsilon_factor | mu_factor | sigma_e | sigma_m
    !% <br>%</tt>
    !%
    !% Permittivity factor, permeability factor, electric conductivity and magnetic conductivity of the medium box.
    !%End

    if (parse_block(namespace, 'LinearMediumProperties', blk) == 0) then

      nlines = parse_block_n(blk)
      if (nlines /=  1) then
        call messages_input_error(namespace, 'LinearMediumProperties', 'should consist of one line', nlines)
      end if
      ncols = parse_block_cols(blk, 0)
      if (ncols /= 4) then
        call messages_input_error(namespace, 'LinearMediumProperties', 'should consist of four columns')
      end if
      call parse_block_float(blk, 0, 0, this%ep_factor)
      call parse_block_float(blk, 0, 1, this%mu_factor)
      call parse_block_float(blk, 0, 2, this%sigma_e_factor)
      call parse_block_float(blk, 0, 3, this%sigma_m_factor)
      write(message(1),'(a,es9.2)') 'Box epsilon factor: ', this%ep_factor
      write(message(2),'(a,es9.2)') 'Box mu factor:      ', this%mu_factor
      write(message(3),'(a,es9.2)') 'Box electric sigma: ', this%sigma_e_factor
      write(message(4),'(a,es9.2)') 'Box magnetic sigma: ', this%sigma_m_factor
      write(message(5),'(a)') ""
      call messages_info(5, namespace=namespace)
      call parse_block_end(blk)

      call messages_print_stress(namespace=namespace)
    else
      message(1) = 'You must specify the properties of your linear medium through the LinearMediumProperties block.'
      call messages_fatal(1, namespace=namespace)
    end if

    !%Variable LinearMediumEdgeProfile
    !%Type integer
    !%Section Maxwell
    !%Description
    !% Defines the type of numerical approximation used for the derivatives at the edges of the medium box.
    !% Default is edged. When the box shape is read from file, only the edged profile is supported.
    !%
    !%Option edged 1
    !% Medium box edges are considered steep for derivatives.
    !%Option smooth 2
    !% Medium box edged and softened for derivatives.
    !%End
    call parse_variable(namespace, 'LinearMediumEdgeProfile', OPTION__LINEARMEDIUMEDGEPROFILE__EDGED, this%edge_profile)
    if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__EDGED) then
      write(message(1),'(a,a)')   'Box shape:          ', 'edged'
    else if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__SMOOTH) then
      write(message(1),'(a,a)')   'Box shape:          ', 'smooth'
    end if
    call messages_info(1, namespace=namespace)

    call this%supported_interactions_as_partner%add(LINEAR_MEDIUM_TO_EM_FIELD)
    this%quantities(PERMITTIVITY)%required = .true.
    this%quantities(PERMEABILITY)%required = .true.
    this%quantities(E_CONDUCTIVITY)%required = .true.
    this%quantities(M_CONDUCTIVITY)%required = .true.

    call output_linear_medium_init(this%outp, this%namespace, this%space)

    call get_medium_box_points_map(this%medium_box, this%gr)
    call get_linear_medium_em_properties(this, this%medium_box, this%gr)
    call output_linear_medium(this%outp, this%namespace, this%space, this%gr%mesh, &
         this%outp%iter_dir, this%medium_box%points_map, this%medium_box%ep, &
         this%medium_box%mu, this%medium_box%c)

    if (this%medium_box%check_medium_points) then
      SAFE_ALLOCATE(tmp(1:this%gr%mesh%np))
      n_global_points = 0
      write(message(1),'(a, a, a)')   'Check of points inside surface of medium ', trim(this%medium_box%filename), ":"
      call messages_info(1, namespace=this%namespace)
      call get_points_map_from_file(this%medium_box%filename, this%gr%mesh, n_points, n_global_points, tmp, CNST(0.99))
      write(message(1),'(a, I8)')'Number of points inside medium (normal coordinates):', &
           this%medium_box%global_points_number
      write(message(2),'(a, I8)')'Number of points inside medium (rescaled coordinates):', n_global_points
      write(message(3), '(a)') ""
      call messages_info(3, namespace=this%namespace)
      SAFE_DEALLOCATE_A(tmp)
    end if

    call profiling_out(prof)

    POP_SUB(linear_medium_init)
  end subroutine linear_medium_init

  ! ---------------------------------------------------------
  subroutine linear_medium_init_interaction(this, interaction)
    class(linear_medium_t), target, intent(inout) :: this
    class(interaction_t),                intent(inout) :: interaction

    PUSH_SUB(linear_medium_init_interaction)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by a linear medium."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(linear_medium_init_interaction)
  end subroutine linear_medium_init_interaction

  ! ---------------------------------------------------------
  subroutine linear_medium_init_interaction_as_partner(partner, interaction)
    class(linear_medium_t),   intent(in)    :: partner
    class(interaction_t),     intent(inout) :: interaction

    integer  :: ip_in, ip_med, ip_mxll, ip_med_gr, ip_mxll_gr
    integer, allocatable:: aux_map(:),  partner_to_system_map(:), medium_inverse_mapping(:)

    PUSH_SUB(linear_medium_init_interaction_as_partner)

    select type (interaction)
    type is (linear_medium_to_em_field_t)
      call grid_transfer_mapping(interaction%system_gr, partner%gr, partner_to_system_map, partner%namespace)

      ! find which medium box points are part of Maxwell grid
      ip_mxll = 0
      allocate(aux_map(partner%medium_box%points_number)) ! at the most it will have these many points
      allocate(medium_inverse_mapping(partner%gr%mesh%np))
      aux_map = 0
      medium_inverse_mapping = 0
      do ip_in = 1, partner%medium_box%points_number
        ip_med_gr = partner%medium_box%points_map(ip_in)
        ip_mxll_gr = partner_to_system_map(ip_med_gr)
        if (ip_mxll_gr > 0) then
          ip_mxll = ip_mxll + 1
          aux_map(ip_mxll) = ip_mxll_gr
          medium_inverse_mapping(ip_mxll) = ip_in
        end if
      end do

      ! allocate and fill medium box on the Maxwell side
      call single_medium_box_allocate(interaction%medium_box, ip_mxll)
      interaction%medium_box%points_map(1:ip_mxll) = aux_map(1:ip_mxll)
      interaction%medium_box%points_number = ip_mxll

      do ip_in = 1, interaction%medium_box%points_number
        ip_med = medium_inverse_mapping(ip_in)
        interaction%medium_box%aux_ep(ip_in,1:3) = partner%medium_box%aux_ep(ip_med,1:3)
        interaction%medium_box%aux_mu(ip_in,1:3) = partner%medium_box%aux_mu(ip_med,1:3)
        interaction%medium_box%ep(ip_in) = partner%medium_box%ep(ip_med)
        interaction%medium_box%mu(ip_in) = partner%medium_box%mu(ip_med)
        interaction%medium_box%c(ip_in) = partner%medium_box%c(ip_med)
        interaction%medium_box%sigma_e(ip_in) = partner%medium_box%sigma_e(ip_med)
        interaction%medium_box%sigma_m(ip_in) = partner%medium_box%sigma_m(ip_med)
      end do
      SAFE_DEALLOCATE_A(aux_map)
      SAFE_DEALLOCATE_A(medium_inverse_mapping)
      SAFE_DEALLOCATE_A(partner_to_system_map)
    class default
      message(1) = "Trying to initialize an unsupported interaction by a linear medium."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(linear_medium_init_interaction_as_partner)
  end subroutine linear_medium_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine linear_medium_initial_conditions(this)
    class(linear_medium_t), intent(inout) :: this

    PUSH_SUB(linear_medium_initial_conditions)

    POP_SUB(linear_medium_initial_conditions)
  end subroutine linear_medium_initial_conditions

  ! ---------------------------------------------------------
  subroutine linear_medium_do_td(this, operation)
    class(linear_medium_t),    intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    PUSH_SUB(linear_medium_do_td)

    select case (operation%id)
    case (SKIP)
      ! Do nothing
    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(linear_medium_do_td)
  end subroutine linear_medium_do_td

  ! ---------------------------------------------------------
  logical function linear_medium_is_tolerance_reached(this, tol) result(converged)
    class(linear_medium_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(linear_medium_is_tolerance_reached)

    ! this routine is never called at present, no reason to be here
    ASSERT(.false.)
    converged = .false.

    POP_SUB(linear_medium_is_tolerance_reached)
  end function linear_medium_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine linear_medium_update_quantity(this, iq)
    class(linear_medium_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(linear_medium_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(linear_medium_update_quantity)
  end subroutine linear_medium_update_quantity

  ! ---------------------------------------------------------
  subroutine linear_medium_update_exposed_quantity(partner, iq)
    class(linear_medium_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(linear_medium_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case (PERMITTIVITY, PERMEABILITY, E_CONDUCTIVITY, M_CONDUCTIVITY)
      call partner%quantities(iq)%clock%set_time(partner%prop%clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(linear_medium_update_exposed_quantity)
  end subroutine linear_medium_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine linear_medium_copy_quantities_to_interaction(partner, interaction)
    class(linear_medium_t),          intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(linear_medium_copy_quantities_to_interaction)

    select type (interaction)
    type is (linear_medium_to_em_field_t)
      ! No need to copy anything
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(linear_medium_copy_quantities_to_interaction)
  end subroutine linear_medium_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine linear_medium_restart_write_data(this)
    class(linear_medium_t), intent(inout) :: this

    integer :: restart_file_unit

    PUSH_SUB(linear_medium_restart_write_data)

    message(1) = "Linear medium system "//trim(this%namespace%get())//" will only write a dummy restart file."
    call messages_warning(1, namespace=this%namespace)

    call io_mkdir('restart/'//TD_DIR, this%namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_linear_medium', this%namespace, action='write')
    write(restart_file_unit, *) 'Linear medium restart file'
    call io_close(restart_file_unit)

    message(1) = "Successfully wrote restart data for system "//trim(this%namespace%get())
    call messages_info(1, namespace=this%namespace)

    POP_SUB(linear_medium_restart_write_data)
  end subroutine linear_medium_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function linear_medium_restart_read_data(this)
    class(linear_medium_t), intent(inout) :: this

    integer :: restart_file_unit

    PUSH_SUB(linear_medium_restart_read_data)

    call io_mkdir('restart/'//TD_DIR, this%namespace, parents=.true.)
    restart_file_unit = io_open('restart/'//TD_DIR// 'restart_linear_medium', this%namespace, action='read', die=.false.)
    if (restart_file_unit > 0) then
      call io_close(restart_file_unit)
      linear_medium_restart_read_data = .true.
    else
      ! could not open file
      linear_medium_restart_read_data = .false.
    end if

    if (linear_medium_restart_read_data) then
      message(1) = "Successfully read restart data for system "//trim(this%namespace%get())
      call messages_info(1, namespace=this%namespace)
    end if

    POP_SUB(linear_medium_restart_read_data)
  end function linear_medium_restart_read_data

  subroutine linear_medium_update_kinetic_energy(this)
    class(linear_medium_t), intent(inout) :: this

    PUSH_SUB(linear_medium_update_kinetic_energy)

    this%kinetic_energy = M_ZERO

    POP_SUB(linear_medium_update_kinetic_energy)

  end subroutine linear_medium_update_kinetic_energy

  ! ---------------------------------------------------------
  subroutine linear_medium_finalize(this)
    type(linear_medium_t), intent(inout) :: this

    PUSH_SUB(linear_medium_finalize)

    call single_medium_box_end(this%medium_box)
    call system_end(this)
    call multicomm_end(this%mc)
    call grid_end(this%gr)

    POP_SUB(linear_medium_finalize)
  end subroutine linear_medium_finalize


  ! Specific routines for this system:

  ! ---------------------------------------------------------
  subroutine get_medium_box_points_map(medium_box, gr)
    type(single_medium_box_t),   intent(inout) :: medium_box
    type(grid_t),                intent(in)    :: gr

    integer :: ip, ip_in, ip_bd, idim
    integer, allocatable :: tmp_points_map(:), tmp_bdry_map(:)
    FLOAT   :: bounds(2,gr%box%dim), xx(gr%box%dim)
    logical :: inside
    type(profile_t), save :: prof

    PUSH_SUB(get_medium_box_points_map)

    call profiling_in(prof, 'GET_MEDIUM_BOX_POINTS_MAP')

    SAFE_ALLOCATE(tmp_points_map(1:gr%mesh%np))
    SAFE_ALLOCATE(tmp_bdry_map(1:gr%mesh%np))
    tmp_points_map = 0
    tmp_bdry_map = 0

    if (medium_box%box_shape == MEDIUM_BOX_FILE) then

      call get_points_map_from_file(medium_box%filename, gr%mesh, medium_box%points_number,&
        medium_box%global_points_number, tmp_points_map)

    else

      do idim = 1, 3
        bounds(1,idim) = medium_box%center(idim) - medium_box%lsize(idim)/M_TWO
        bounds(2,idim) = medium_box%center(idim) + medium_box%lsize(idim)/M_TWO
      end do
      ip_in = 0
      ip_bd = 0
      do ip = 1, gr%mesh%np
        xx(1:3) = gr%mesh%x(ip, 1:3)
        inside = check_point_in_bounds(xx, bounds(:,:))
        if (check_point_in_bounds(xx, bounds(:,:))) then
          ip_in = ip_in + 1
          tmp_points_map(ip_in) = ip
        end if
        if (check_point_on_bounds(xx, bounds(:,:))) then
          ip_bd = ip_bd + 1
          tmp_bdry_map(ip_bd) = ip
        end if
      end do

      medium_box%points_number = ip_in
      medium_box%bdry_number = ip_bd

      SAFE_ALLOCATE(medium_box%bdry_map(1:medium_box%bdry_number))
      medium_box%bdry_map = 0
      medium_box%bdry_map = tmp_bdry_map(1:medium_box%bdry_number)

    end if
    call single_medium_box_allocate(medium_box, medium_box%points_number)
    medium_box%points_map(:) = tmp_points_map(1:medium_box%points_number)

    ! TODO: add some check that medium boxes do not overlap (now they are different systems)
    SAFE_DEALLOCATE_A(tmp_points_map)
    SAFE_DEALLOCATE_A(tmp_bdry_map)

    call profiling_out(prof)

    POP_SUB(get_medium_box_points_map)
  contains

    logical pure function check_point_in_bounds(xx, bounds) result (check)
      FLOAT, intent(in) :: xx(:)
      FLOAT, intent(in) :: bounds(:,:)

      check = .false.
      if ((xx(1) >= bounds(1,1)) .and. (xx(1) <= bounds(2,1)) .and. &
        (xx(2) >= bounds(1,2)) .and. (xx(2) <= bounds(2,2)) .and. &
        (xx(3) >= bounds(1,3)) .and. (xx(3) <= bounds(2,3)) ) then
        check = .true.
      end if

    end function check_point_in_bounds

    logical pure function check_point_on_bounds(xx, bounds) result (check)
      FLOAT, intent(in) :: xx(:)
      FLOAT, intent(in) :: bounds(:,:)

      check = .false.
      if (xx(1) == bounds(1,1) .and. (xx(2) >= bounds(1,2) .and. xx(3) >= bounds(1,3)) &
        .and. (xx(2) <= bounds(2,2) .and. xx(3) <= bounds(2,3)) .or. &
        xx(2) == bounds(1,2) .and. (xx(1) >= bounds(1,1) .and. xx(3) >= bounds(1,3)) &
        .and. (xx(1) <= bounds(2,1) .and. xx(3) <= bounds(2,3)) .or. &
        xx(3) == bounds(1,3) .and. (xx(1) >= bounds(1,1) .and. xx(2) >= bounds(1,2)) &
        .and. (xx(1) <= bounds(2,1) .and. xx(2) <= bounds(2,2)) .or. &
        xx(1) == bounds(2,1) .and. (xx(2) >= bounds(1,2) .and. xx(3) >= bounds(1,3)) &
        .and. (xx(2) <= bounds(2,2) .and. xx(3) <= bounds(2,3)) .or. &
        xx(2) == bounds(2,2) .and. (xx(1) >= bounds(1,1) .and. xx(3) >= bounds(1,3)) &
        .and. (xx(1) <= bounds(2,1) .and. xx(3) <= bounds(2,3)) .or. &
        xx(3) == bounds(2,3) .and. (xx(1) >= bounds(1,1) .and. xx(2) >= bounds(1,2)) &
        .and. (xx(1) <= bounds(2,1) .and. xx(2) <= bounds(2,2)) ) then
        check = .true.
      end if

    end function check_point_on_bounds

  end subroutine get_medium_box_points_map

  ! ----------------------------------------------------------
  !> Evaluate electromagnetic properties of linear medium
  subroutine get_linear_medium_em_properties(this, medium_box, gr)
    type(linear_medium_t), intent(in)          :: this
    type(single_medium_box_t),   intent(inout) :: medium_box
    type(grid_t),                intent(in)    :: gr

    integer :: ip, ip_in, ip_bd, ipp
    FLOAT   :: xx(gr%box%dim), xxp(gr%box%dim), dd, dd_max, dd_min
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)
    type(profile_t), save :: prof

    PUSH_SUB(get_linear_medium_em_properties)

    call profiling_in(prof, 'GET_LINEAR_MEDIUM_EM_PROPERTIES')

    SAFE_ALLOCATE(tmp(1:gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(1:gr%mesh%np_part,1:gr%box%dim))
    dd_max = max(2*gr%mesh%spacing(1), 2*gr%mesh%spacing(2), 2*gr%mesh%spacing(3))

    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__SMOOTH) then
        ASSERT(allocated(medium_box%bdry_map))

        xx(1:3) = gr%mesh%x(ip,1:3)
        dd_min = M_HUGE

        do ip_bd = 1, medium_box%bdry_number
          ipp = medium_box%bdry_map(ip_bd)
          xxp(1:3) = gr%mesh%x(ipp,1:3)
          dd = norm2(xx(1:3) - xxp(1:3))
          if (dd < dd_min) dd_min = dd
        end do

        medium_box%ep(ip_in) = P_ep + ((P_ep * this%ep_factor - P_ep)  &
          * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max))))
        medium_box%mu(ip_in) = P_mu + ((P_mu * this%mu_factor - P_mu) &
          * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max))))
        medium_box%c(ip_in) = M_ONE/sqrt(medium_box%ep(ip_in)*medium_box%mu(ip_in))
        medium_box%sigma_e(ip_in) = this%sigma_e_factor &
          * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max)))
        medium_box%sigma_m(ip_in) = this%sigma_m_factor &
          * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max)))

      else if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__EDGED) then

        medium_box%ep(ip_in) = P_ep * this%ep_factor
        medium_box%mu(ip_in) = P_mu * this%mu_factor
        medium_box%c(ip_in) = M_ONE/sqrt(medium_box%ep(ip_in)*medium_box%mu(ip_in))
        medium_box%sigma_e(ip_in) = this%sigma_e_factor
        medium_box%sigma_m(ip_in) = this%sigma_m_factor

      end if
    end do

    tmp(:) = P_ep
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      tmp(ip)= medium_box%ep(ip_in)
    end do
    call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      medium_box%aux_ep(ip_in, :) = tmp_grad(ip, :)/(M_FOUR * medium_box%ep(ip_in))
    end do

    tmp(:) = P_mu
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      tmp(ip) = medium_box%mu(ip_in)
    end do
    call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      medium_box%aux_mu(ip_in, :) = tmp_grad(ip, :)/(M_FOUR * medium_box%mu(ip_in))
    end do

    !TODO: add print information about the medium box

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_grad)

    call profiling_out(prof)

    POP_SUB(get_linear_medium_em_properties)
  end subroutine get_linear_medium_em_properties

  ! ----------------------------------------------------------
  !> Populate list of point indices for points inside the polyhedron
  subroutine get_points_map_from_file(filename, mesh, n_points, global_points_number, tmp_map, scale_factor)
    character(len=256),       intent(in)    :: filename
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(out)   :: n_points
    integer,                  intent(out)   :: global_points_number
    integer,                  intent(inout) :: tmp_map(:)
    FLOAT, optional,          intent(in)    :: scale_factor

    integer :: ip_in, ip
    FLOAT   :: xx(3)
    type(cgal_polyhedra_t) :: cgal_poly
    type(profile_t), save :: prof

    PUSH_SUB(get_points_map_from_file)

    call profiling_in(prof, 'GET_POINTS_MAP_FROM_FILE')

    call cgal_polyhedron_init(cgal_poly, trim(filename), verbose = .false.)

    ip_in = 0
    do ip = 1, mesh%np
      if (present(scale_factor)) then
        xx(1:3) = scale_factor * mesh%x(ip, 1:3)
      else
        xx(1:3) = mesh%x(ip, 1:3)
      end if
      if (cgal_polyhedron_point_inside(cgal_poly, xx(1), xx(2), xx(3))) then
        ip_in = ip_in + 1
        tmp_map(ip_in) = ip
      end if
    end do
    n_points = ip_in
    call cgal_polyhedron_end(cgal_poly)

    call mpi_world%allreduce(ip_in, global_points_number, 1, MPI_INTEGER, MPI_SUM)

    call profiling_out(prof)

    POP_SUB(get_points_map_from_file)

  end subroutine get_points_map_from_file

  ! ----------------------------------------------------------
  !> Parse and store geometry of medium box
  subroutine medium_box_init(medium_box, namespace)
    type(single_medium_box_t), intent(inout) :: medium_box
    type(namespace_t),      intent(in)       :: namespace

    integer :: nlines, ncols, idim
    type(block_t) :: blk

    PUSH_SUB(medium_box_init)

    !%Variable LinearMediumBoxShape
    !%Type integer
    !%Section Maxwell
    !%Description
    !% This variable defines the shape of the linear medium box.
    !% The default is <tt>medium_parallelepiped</tt>.
    !%Option medium_parallelepiped 1
    !% The medium box will be a parallelepiped whose center and dimensions are taken from
    !% the variable <tt>LinearMediumBoxSize</tt>.
    !%Option medium_box_file 2
    !% The simulation box will be read from an external file in OFF format, defined by the variable <tt>LinearMediumBoxFile</tt>.
    !%End
    call parse_variable(namespace, 'LinearMediumBoxShape', MEDIUM_PARALLELEPIPED, medium_box%box_shape)

    if (.not. varinfo_valid_option('LinearMediumBoxShape', medium_box%box_shape)) then
      call messages_input_error(namespace, 'LinearMediumBoxShape')
    end if

    select case (medium_box%box_shape)
    case (MEDIUM_PARALLELEPIPED)

      !%Variable LinearMediumBoxSize
      !%Type block
      !%Section Maxwell
      !%Description
      !% Defines center and size of a parallelepiped linear medium box.
      !%
      !% Example:
      !%
      !% <tt>%LinearMediumBoxSize
      !% <br>&nbsp;&nbsp;   center_x | center_y | center_z | x_length | y_length | z_length
      !% <br>%</tt>
      !%End

      if (parse_block(namespace, 'LinearMediumBoxSize', blk) == 0) then
        call messages_print_stress(msg=trim('Linear Medium box center and size:'), namespace=namespace)

        nlines = parse_block_n(blk)
        if (nlines /=  1) then
          call messages_input_error(namespace, 'LinearMediumBoxSize', 'should consist of one line')
        end if
        ncols = parse_block_cols(blk, 0)
        if (ncols /= 6) then
          call messages_input_error(namespace, 'LinearMediumBoxSize', 'should consist of six columns')
        end if
        do idim = 1, 3
          call parse_block_float(blk, 0, idim-1, medium_box%center(idim))
          call parse_block_float(blk, 0, idim+2, medium_box%lsize(idim))
        end do
        write(message(1),'(a,es9.2,a,es9.2,a,es9.2)') 'Box center:         ', medium_box%center(1), ' | ',&
          medium_box%center(2), ' | ', medium_box%center(3)
        write(message(2),'(a,es9.2,a,es9.2,a,es9.2)') 'Box size  :         ', medium_box%lsize(1), ' | ', &
          medium_box%lsize(2), ' | ', medium_box%lsize(3)
        write(message(3),'(a)') ""
        call messages_info(3)
        call parse_block_end(blk)

        call messages_print_stress(namespace=namespace)
      else
        message(1) = "For parallelepiped box shapes, you must provide a LinearMediumBoxSize block."
        call messages_fatal(1, namespace=namespace)
      end if

    case (MEDIUM_BOX_FILE)
      !%Variable LinearMediumBoxFile
      !%Type string
      !%Section Maxwell
      !%Description
      !% File in OFF format with the shape of the linear medium.
      !%End
      if (parse_is_defined(namespace, 'LinearMediumBoxFile')) then
        call parse_variable(namespace, 'LinearMediumBoxFile', 'mediumboxfile', medium_box%filename)
      else
        message(1) = "When using box_file as the box shape, you must provide a filename through the LinearMediumBoxFile variable."
        call messages_fatal(1, namespace=namespace)
      end if

      !%Variable CheckPointsMediumFromFile
      !%Type logical
      !%Default no
      !%Section Maxwell
      !%Description
      !% Whether to re-calculate the points map by artificially shrinking the coordinate system by a factor of
      !% 0.99 to check if the points inside the medium surface are properly detected. This works for only one
      !% medium surface which is centered in the origin of the coordinate system.
      !%End
      call parse_variable(namespace, 'CheckPointsMediumFromFile', .false., medium_box%check_medium_points)

    end select

    POP_SUB(medium_box_init)

  end subroutine medium_box_init

end module linear_medium_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
