!! Copyright (C) 2023 E.I. Albar, F. Bonafe and Heiko Appel
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WAx_normANTY; without even the implied wax_normanty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.

#include "global.h"

module plane_wave_oct_m
  use accel_oct_m
  use box_sphere_oct_m
  use box_parallelepiped_oct_m
  use cube_function_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use index_oct_m
  use io_oct_m
  use io_function_oct_m
  use maxwell_function_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use string_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use space_oct_m
  use states_mxll_oct_m

  implicit none

  private
  public ::                    &
    plane_wave_t,              &
    plane_wave_init,           &
    plane_waves_eval,          &
    plane_wave_end

  type plane_wave_t
    integer                          :: points_number  !< number of points of plane wave boundary
    integer,             allocatable :: points_map(:) !< points map for plane waves boundary
    integer                          :: number !< number of plane waves given by user
    integer,             allocatable :: modus(:) !< input file modus, either parser or Maxwell function
    character(len=1024), allocatable :: e_field_string(:,:) !< string in case of parser
    FLOAT,               allocatable :: k_vector(:,:) !< k vector for each plane wave
    FLOAT,               allocatable :: v_vector(:,:) !< velocity vector for each plane wave
    CMPLX,               allocatable :: e_field(:,:) !< field amplitude for each plane wave
    type(mxf_t),         allocatable :: mx_function(:) !< Maxwell function for each plane wave
    type(space_t), private :: space
    logical                          :: evaluate_on_one_side = .false.
    integer                          :: side_of_the_box
    type(accel_mem_t)                :: buff_map
  end type plane_wave_t

contains
  !> @brief Here, plane wave is evaluated from analytical formulae on grid.
  ! ---------------------------------------------------------
  subroutine plane_wave_init(plane_wave, namespace)
    type(plane_wave_t),        intent(inout) :: plane_wave
    type(namespace_t),      intent(in)    :: namespace
    type(block_t)        :: blk
    integer              :: il, nlines, ncols, iex_norm, idim
    FLOAT                :: k_vector(3), velocity(3), xx(3), x_norm, dummy(3), test, test_limit!, angle, sigma
    CMPLX                :: e_field(3)
    character(len=1024)  :: k_string(3)
    character(len=1024)  :: mxf_expression
    type(profile_t), save :: prof

    PUSH_SUB(plane_wave_init)

    call profiling_in(prof, 'PLANE_WAVE_INIT')

    test_limit = CNST(10.0e-9)
    plane_wave%space%dim = 3

    !%Variable MaxwellIncidentWaves
    !%Type block
    !%Section MaxwellStates
    !%Description
    !% The initial electromagnetic fields can be set by the user
    !% with the <tt>MaxwellIncidentWaves</tt> block variable.
    !% The electromagnetic fields have to fulfill the
    !% Maxwells equations in vacuum.
    !%
    !% Example:
    !%
    !% <tt>%MaxwellIncidentWaves
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | "k1x" | "k1y" | "k1z" | "E1x" | "E1z" | "E1x"
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | "k2x" | "k2y" | "k2z" | "E2x" | "E2y" | "E2z"
    !% <br>&nbsp;&nbsp;   plane_wave_gauss       | "k3x" | "k3y" | "k3z" | "E3x" | "E3y" | "E3z" | "width" | "shift"
    !% <br>&nbsp;&nbsp;   plane_wave_mx_function | "E4x" | "E4y" | "E4z" | mx_envelope_name
    !% <br>%</tt>
    !%
    !% Description about MaxwellIncidentWaves follows
    !%
    !%Option plane_wave_parser 0
    !% Parser input modus
    !%Option plane_wave_mx_function 1
    !% The incident wave envelope is defined by an mx_function
    !%End

    if (parse_block(namespace, 'MaxwellIncidentWaves', blk) == 0) then

      call messages_print_with_emphasis(msg='Substitution of the electromagnetic incident waves', namespace=namespace)

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)

      plane_wave%number = nlines
      SAFE_ALLOCATE(plane_wave%modus(1:nlines))
      SAFE_ALLOCATE(plane_wave%e_field_string(1:plane_wave%space%dim, 1:nlines))
      SAFE_ALLOCATE(plane_wave%e_field(1:plane_wave%space%dim, 1:nlines))
      SAFE_ALLOCATE(plane_wave%k_vector(1:plane_wave%space%dim, 1:nlines))
      SAFE_ALLOCATE(plane_wave%v_vector(1:plane_wave%space%dim, 1:nlines))
      SAFE_ALLOCATE(plane_wave%mx_function(1:nlines))

      ! read all lines
      do il = 1, nlines
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, il - 1)
        if ((ncols /= 5) .and. (ncols /= 7) .and. (ncols /= 9)) then
          message(1) = 'Each line in the MaxwellIncidentWaves block must have five, seven or nine columns.'
          call messages_fatal(1, namespace=namespace)
        end if

        ! check input modus e.g. parser of defined functions
        call parse_block_integer(blk, il - 1, 0, plane_wave%modus(il))

        ! parse formula string
        if (plane_wave%modus(il) == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_PARSER) then

          call parse_block_string( blk, il - 1, 1, k_string(1))
          call parse_block_string( blk, il - 1, 2, k_string(2))
          call parse_block_string( blk, il - 1, 3, k_string(3))
          call parse_block_string( blk, il - 1, 4, plane_wave%e_field_string(1, il))
          call parse_block_string( blk, il - 1, 5, plane_wave%e_field_string(2, il))
          call parse_block_string( blk, il - 1, 6, plane_wave%e_field_string(3, il))

          write(message(1), '(a,i2,a) ') 'Substituting electromagnetic incident wave ', il, ' with the expressions: '
          call messages_info(1, namespace=namespace)
          write(message(1), '(6a)')     '  Wave vector k(x)   = ', trim(k_string(1))
          write(message(2), '(2a)')     '  Wave vector k(y)   = ', trim(k_string(2))
          write(message(3), '(2a)')     '  Wave vector k(z)   = ', trim(k_string(3))
          write(message(4), '(2a)')     '  E-field(x) for t_0 = ', trim(plane_wave%e_field_string(1, il))
          write(message(5), '(2a)')     '  E-field(y) for t_0 = ', trim(plane_wave%e_field_string(2, il))
          write(message(6), '(2a)')     '  E-field(z) for t_0 = ', trim(plane_wave%e_field_string(3, il))
          call messages_info(6, namespace=namespace)

          call conv_to_C_string(k_string(1))
          call conv_to_C_string(k_string(2))
          call conv_to_C_string(k_string(3))
          call conv_to_C_string(plane_wave%e_field_string(1, il))
          call conv_to_C_string(plane_wave%e_field_string(2, il))
          call conv_to_C_string(plane_wave%e_field_string(3, il))

          xx(:) = M_ZERO
          x_norm    = M_ZERO
          call parse_expression(k_vector(1), dummy(1), 1, xx, x_norm, M_ZERO, k_string(1))
          call parse_expression(k_vector(2), dummy(2), 2, xx, x_norm, M_ZERO, k_string(2))
          call parse_expression(k_vector(3), dummy(3), 3, xx, x_norm, M_ZERO, k_string(3))
          k_vector = units_to_atomic(unit_one/units_inp%length, k_vector)

          velocity(:)    = k_vector(:) / norm2(k_vector) * P_c
          plane_wave%k_vector(:,il) = k_vector(:)
          plane_wave%v_vector(:,il) = velocity(:)

        else if (plane_wave%modus(il) == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
          do idim = 1, size(e_field)
            call parse_block_cmplx( blk, il - 1, idim, e_field(idim))
          end do
          call parse_block_string( blk, il - 1, size(e_field) + 1, mxf_expression)

          write(message(1), '(a,i2) ') 'Substituting electromagnetic incident wave ', il
          write(message(2), '(a)'    ) 'with the expression: '
          call messages_info(2, namespace=namespace)
          write(message(1), '(a,f9.4,sp,f9.4,"i")') '  E-field(x) complex amplitude  = ', real(e_field(1)), aimag(e_field(1))
          write(message(2), '(a,f9.4,sp,f9.4,"i")') '  E-field(y) complex amplitude  = ', real(e_field(2)), aimag(e_field(2))
          write(message(3), '(a,f9.4,sp,f9.4,"i")') '  E-field(z) complex amplitude  = ', real(e_field(3)), aimag(e_field(3))
          write(message(4), '(2a)'    )      '  Maxwell wave function name = ', trim(mxf_expression)
          call messages_info(4, namespace=namespace)
          call mxf_read(plane_wave%mx_function(il), namespace, trim(mxf_expression), iex_norm)
          if (iex_norm /= 0) then
            write(message(1),'(3A)') 'Ex_normor in the ""', trim(mxf_expression), &
              '"" field defined in the MaxwellIncidentWaves block'
            call messages_fatal(1, namespace=namespace)
          end if
          e_field  = units_to_atomic(units_inp%energy/units_inp%length, e_field)
          k_vector(1:3) = plane_wave%mx_function(il)%k_vector(1:3)

          test = TOFLOAT(dot_product(k_vector, e_field))
          if (abs(test) > test_limit) then
            message(1) = 'The wave vector k or its electric field E-field '
            message(2) = 'is not perpendicular enough.'
            call messages_fatal(2, namespace=namespace)
          end if
          if (norm2(k_vector) < 1e-10) then
            message(1) = 'The k vector is not defined set correctly.'
            call messages_fatal(1, namespace=namespace)
          end if

          plane_wave%e_field(:,il)  = e_field(:)
          plane_wave%k_vector(:,il) = k_vector(:)
          plane_wave%v_vector(:,il) = k_vector(:) / norm2(k_vector) * P_c

        end if
      end do

      call parse_block_end(blk)

      call messages_print_with_emphasis(namespace=namespace)
    else
      plane_wave%number = 0

    end if

    call profiling_out(prof)

    POP_SUB(plane_wave_init)
  end subroutine plane_wave_init

  ! ---------------------------------------------------------
  subroutine plane_wave_end(plane_wave)
    type(plane_wave_t),   intent(inout) :: plane_wave

    PUSH_SUB(plane_wave_end)

    SAFE_DEALLOCATE_A(plane_wave%points_map)
    SAFE_DEALLOCATE_A(plane_wave%modus)
    SAFE_DEALLOCATE_A(plane_wave%e_field_string)
    SAFE_DEALLOCATE_A(plane_wave%k_vector)
    SAFE_DEALLOCATE_A(plane_wave%v_vector)
    SAFE_DEALLOCATE_A(plane_wave%e_field)
    SAFE_DEALLOCATE_A(plane_wave%mx_function)
    if (accel_is_enabled()) then
      call accel_release_buffer(plane_wave%buff_map)
    end if

    POP_SUB(plane_wave_end)
  end subroutine plane_wave_end

  ! ---------------------------------------------------------

  !> Calculation of plane waves from parsed formula
  subroutine plane_waves_eval(plane_wave, time, mesh, e_field_total, b_field_total)
    type(plane_wave_t),        intent(inout) :: plane_wave
    FLOAT,                     intent(in)    :: time
    class(mesh_t),             intent(in)    :: mesh
    FLOAT,                     intent(inout) :: e_field_total(:, :)
    FLOAT, optional,           intent(inout) :: b_field_total(:, :)

    integer              :: ip, wn, idim, i
    FLOAT                :: x_prop(plane_wave%space%dim), x_norm                    !< Propagated position
    FLOAT                :: velocity(plane_wave%space%dim)                          !< PW velocity
    FLOAT                :: velocity_time(plane_wave%space%dim)                     !< Velocity times time
    FLOAT                :: k_vector(plane_wave%space%dim), k_vector_abs            !< PW k-vector and its norm
    FLOAT                :: e_field(plane_wave%space%dim)                           !< E and B field
    FLOAT                :: dummy(plane_wave%space%dim)                             !< Dummy array required for `parse_expression` call
    CMPLX                :: e0(plane_wave%space%dim)                                !< Modulus of E field
    type(profile_t), save :: prof

    integer, allocatable :: indices_pw_parser(:)      !< PW modus indices that correspond to PLANE_WAVE_PARSER option
    integer, allocatable :: indices_mx_ftc(:)         !< PW modus indices that correspond to PLANE_WAVE_MX_FUNCTION option

    PUSH_SUB(plane_waves_eval)

    call profiling_in(prof, 'PLANE_WAVES_EVAL')

    e_field_total(:,:) = M_ZERO
    ! TODO(Ilke): Issue #708 Refactor plane wave and vector potential
    ! TODO(Ilke): Issue 708 Swap the indices of e_field_total, it will speed up the calculation.
    ! TODO(Ilke): Issue 708 Fix double access to memory in e_field_total.

    indices_pw_parser = pack([(wn, wn = 1,plane_wave%number)], &
      plane_wave%modus == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_PARSER)

    do i = 1, size(indices_pw_parser)
      wn = indices_pw_parser(i)
      velocity_time(:) = plane_wave%v_vector(1:plane_wave%space%dim, wn) * time
      do ip = 1, mesh%np
        x_prop = mesh%x(ip, :) - velocity_time
        x_norm = norm2(x_prop(1:plane_wave%space%dim))
        do idim = 1, plane_wave%space%dim
          call parse_expression(e_field(idim), dummy(idim), plane_wave%space%dim, x_prop, x_norm, M_ZERO, &
            plane_wave%e_field_string(idim, wn))
          e_field_total(ip, idim) = e_field_total(ip, idim) + units_to_atomic(units_inp%energy/units_inp%length, e_field(idim))
        end do
      end do
    end do

    indices_mx_ftc = pack([(wn, wn = 1,plane_wave%number)], &
      plane_wave%modus == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION)

    do i = 1, size(indices_mx_ftc)
      wn = indices_mx_ftc(i)
      velocity_time(:) = plane_wave%v_vector(1:plane_wave%space%dim, wn) * time
      e0(:) = plane_wave%e_field(1:plane_wave%space%dim, wn)
      do ip = 1, mesh%np
        x_prop = mesh%x(ip, :) - velocity_time
        x_norm = norm2(x_prop(1:plane_wave%space%dim))
        e_field_total(ip, :) = e_field_total(ip, :) + &
          TOFLOAT(e0(1:plane_wave%space%dim) * mxf(plane_wave%mx_function(wn), x_prop(1:plane_wave%space%dim)))
      end do
    end do

    if (present(b_field_total)) then
      if (plane_wave%number == 0)  then
        call profiling_out(prof)
        POP_SUB(plane_waves_eval)
        return
      end if
      wn = 1
      velocity = plane_wave%v_vector(1:plane_wave%space%dim, wn)
      k_vector = plane_wave%k_vector(1:plane_wave%space%dim, wn)
      k_vector_abs = norm2(k_vector(1:plane_wave%space%dim))
      e0 = plane_wave%e_field(1:plane_wave%space%dim, wn)

      do ip = 1, mesh%np
        b_field_total(ip, :) = M_ONE/(P_c * k_vector_abs) * dcross_product(k_vector, e_field_total(ip, :))
      end do

      do wn = 2, plane_wave%number
        velocity = plane_wave%v_vector(1:plane_wave%space%dim, wn)
        k_vector = plane_wave%k_vector(1:plane_wave%space%dim, wn)
        k_vector_abs = norm2(k_vector(1:plane_wave%space%dim))
        e0 = plane_wave%e_field(1:plane_wave%space%dim, wn)
        do ip = 1, mesh%np
          b_field_total(ip, :) = b_field_total(ip, :) + M_ONE/(P_c * k_vector_abs) * dcross_product(k_vector, &
            e_field_total(ip, :))
        end do
      end do
    end if

    call profiling_out(prof)

    POP_SUB(plane_waves_eval)
  end subroutine plane_waves_eval

end module plane_wave_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
