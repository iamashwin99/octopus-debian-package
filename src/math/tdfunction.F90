!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!>--------------------------------------------------------------
!! This module defines "time-dependent functions", to be used by
!! the lasers module, or in the future in order to define time-dependent
!! magnetic fields.
!!--------------------------------------------------------------
module tdfunction_oct_m
  use iso_c_binding
  use debug_oct_m
  use fft_oct_m
  use global_oct_m
  use io_oct_m
  use loct_math_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use splines_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                      &
    tdf_t,                       &
    tdf_init,                    &
    tdf_init_cw,                 &
    tdf_init_gaussian,           &
    tdf_init_cosinoidal,         &
    tdf_init_trapezoidal,        &
    tdf_init_fromfile,           &
    tdf_init_fromexpr,           &
    tdf_init_numerical,          &
    tdf_set_numerical,           &
    tdf_set_random,              &
    tdf_to_numerical,            &
    tdf,                         &
    tdf_dot_product,             &
    tdf_diff,                    &
    tdf_scalar_multiply,         &
    tdf_cosine_multiply,         &
    tdf_numerical_to_fourier,    &
    tdf_fourier_to_numerical,    &
    tdf_numerical_to_zerofourier,&
    tdf_zerofourier_to_numerical,&
    tdf_fourier_grid,            &
    tdf_write,                   &
    tdf_niter,                   &
    tdf_nfreqs,                  &
    tdf_dt,                      &
    tdf_copy,                    &
    tdf_read,                    &
    tdf_is_empty,                &
    tdf_end


  integer, public, parameter ::  &
    TDF_EMPTY         =  10001,  &
    TDF_CW            =  10002,  &
    TDF_GAUSSIAN      =  10003,  &
    TDF_COSINOIDAL    =  10004,  &
    TDF_TRAPEZOIDAL   =  10005,  &
    TDF_FROM_FILE     =  10006,  &
    TDF_NUMERICAL     =  10007,  &
    TDF_FROM_EXPR     =  10008,  &
    TDF_FOURIER_SERIES=  10010,  &
    TDF_ZERO_FOURIER  =  10011

  type tdf_t
    private
    integer :: mode        = TDF_EMPTY
    FLOAT   :: t0          = M_ZERO  !< the time at the maximum of the pulse
    FLOAT   :: tau0        = M_ZERO  !< the width of the pulse
    FLOAT   :: tau1        = M_ZERO  !< for the ramped shape, the length of the "ramping" intervals
    FLOAT   :: a0          = M_ZERO
    FLOAT   :: omega0      = M_ZERO
    FLOAT   :: dt          = M_ZERO !< the time-discretization value.
    FLOAT   :: init_time   = M_ZERO
    FLOAT   :: final_time  = M_ZERO
    integer :: niter       = 0
    integer :: nfreqs      = 0

    type(spline_t)         :: amplitude
    character(len=1024)     :: expression
    FLOAT, allocatable :: val(:)
    FLOAT, allocatable :: valww(:)
    type(fft_t) :: fft_handler
  end type tdf_t

  interface tdf_set_numerical
    module procedure tdf_set_numericalr, tdf_set_numericalr1
  end interface tdf_set_numerical

  interface tdf
    module procedure tdfi, tdft
  end interface tdf

contains

  !------------------------------------------------------------
  !> This function initializes "f" from the TDFunctions block.
  subroutine tdf_read(f, namespace, function_name, ierr)
    type(tdf_t),       intent(inout) :: f
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: function_name
    integer,           intent(out)   :: ierr  !< Error code, 0 on success.

    type(block_t) :: blk
    integer :: nrows, i, function_type
    character(len=1024) :: row_name, filename, function_expression
    FLOAT :: a0, tau0, t0, tau1

    PUSH_SUB(tdf_read)

    !%Variable TDFunctions
    !%Type block
    !%Section Time-Dependent
    !%Description
    !% This block specifies the shape of a "time-dependent function", such as the
    !% envelope needed when using the <tt>TDExternalFields</tt> block. Each line in the block
    !% specifies one function. The first element of each line will be a string
    !% that defines the name of the function. The second element specifies which type
    !% of function we are using; in the following we provide an example for each of the
    !% possible types:
    !%
    !%Option tdf_cw 10002
    !%
    !% <tt>%TDFunctions
    !% <br>&nbsp;&nbsp; "function-name" | tdf_cw | amplitude
    !% <br>%</tt>
    !%
    !% The function is just a constant of value <tt>amplitude</tt>: <math> f(t) </math> = amplitude
    !%
    !%Option tdf_gaussian 10003
    !%
    !% <tt>%TDFunctions
    !% <br>&nbsp;&nbsp; "function-name" | tdf_gaussian | amplitude | tau0 | t0
    !% <br>%</tt>
    !%
    !% The function is a Gaussian, <math> f(t) = F_0 \exp( - (t-t_0)^2/(2\tau_0^2) ) </math>,
    !% where <math>F_0</math> = amplitude.
    !%
    !%Option tdf_cosinoidal 10004
    !%
    !% <tt>%TDFunctions
    !% <br>&nbsp;&nbsp; "function-name" | tdf_cosinoidal | amplitude | tau0 | t0
    !% <br>%</tt>
    !%
    !% <math> f(t) =  F_0 \cos( \frac{\pi}{2} \frac{t-2\tau_0-t_0}{\tau0} )  </math>
    !%
    !% If <math> | t - t_0 | > \tau_0 </math>, then <math> f(t) = 0 </math>.
    !%
    !%Option tdf_trapezoidal 10005
    !%
    !% <tt>%TDFunctions
    !% <br>&nbsp;&nbsp; "function-name" | tdf_trapezoidal | amplitude | tau0 | t0 | tau1
    !% <br>%</tt>
    !%
    !% This function is a trapezoidal centered around <tt>t0</tt>. The
    !% shape is determined by <tt>tau0</tt> and <tt>tau1</tt>. The
    !% function ramps linearly for <tt>tau1</tt> time units, stays
    !% constant for <tt>tau0</tt> time units, and then decays to zero
    !% linearly again for <tt>tau1</tt> time units.
    !%
    !%Option tdf_from_file 10006
    !%
    !% <tt>%TDFunctions
    !% <br>&nbsp;&nbsp; "function-name" | tdf_from_file | "filename"
    !% <br>%</tt>
    !%
    !% The temporal shape of the function is contained in a file called <tt>filename</tt>. This file
    !% should contain three columns: first column is time, second and third column are the
    !% real part and the imaginary part of the temporal function <i>f</i>(<i>t</i>).
    !%
    !%Option tdf_from_expr 10008
    !%
    !% <tt>%TDFunctions
    !% <br>&nbsp;&nbsp; "function-name" | tdf_from_expr | "expression"
    !% <br>%</tt>
    !%
    !% The temporal shape of the field is given as an expression (e.g., <tt>cos(2.0*t)</tt>. The
    !% letter <i>t</i> means time, obviously. The expression is used to construct the function <i>f</i>
    !% that defines the field.
    !%End
    ierr = -3
    if (parse_block(namespace, 'TDFunctions', blk) /= 0) then
      ierr = -1
      POP_SUB(tdf_read)
      return
    end if

    nrows = parse_block_n(blk)
    row_loop: do i = 1, nrows
      call parse_block_string(blk, i-1, 0, row_name)
      if (trim(row_name) == trim(function_name)) then

        call parse_block_integer(blk, i-1, 1, function_type)

        a0   = M_ZERO
        tau0 = M_ZERO
        t0   = M_ZERO
        tau1 = M_ZERO
        select case (function_type)
        case (TDF_CW)
          call parse_block_float(blk, i-1, 2, a0, units_inp%energy/units_inp%length)
          call tdf_init_cw(f, a0, M_ZERO)
        case (TDF_GAUSSIAN)
          call parse_block_float(blk, i-1, 2, a0, units_inp%energy/units_inp%length)
          call parse_block_float(blk, i-1, 3, tau0, units_inp%time)
          call parse_block_float(blk, i-1, 4, t0, units_inp%time)
          call tdf_init_gaussian(f, a0, M_ZERO, t0, tau0)
        case (TDF_COSINOIDAL)
          call parse_block_float(blk, i-1, 2, a0, units_inp%energy/units_inp%length)
          call parse_block_float(blk, i-1, 3, tau0, units_inp%time)
          call parse_block_float(blk, i-1, 4, t0, units_inp%time)
          call tdf_init_cosinoidal(f, a0, M_ZERO, t0, tau0)
        case (TDF_TRAPEZOIDAL)
          call parse_block_float(blk, i-1, 2, a0, units_inp%energy/units_inp%length)
          call parse_block_float(blk, i-1, 3, tau0, units_inp%time)
          call parse_block_float(blk, i-1, 4, t0, units_inp%time)
          call parse_block_float(blk, i-1, 5, tau1, units_inp%time)
          call tdf_init_trapezoidal(f, a0, M_ZERO, t0, tau0, tau1)
        case (TDF_FROM_FILE)
          call parse_block_string(blk, i-1, 2, filename)
          call tdf_init_fromfile(f, trim(filename), namespace, ierr)
        case (TDF_FROM_EXPR)
          call parse_block_string(blk, i-1, 2, function_expression)
          call tdf_init_fromexpr(f, trim(function_expression))
        case default
          ierr = -2
          call parse_block_end(blk)
          POP_SUB(tdf_read)
          return
        end select

        ierr = 0
        exit row_loop
      end if
    end do row_loop

    call parse_block_end(blk)
    POP_SUB(tdf_read)
  end subroutine tdf_read
  !------------------------------------------------------------


  !------------------------------------------------------------
  integer pure function tdf_niter(f)
    type(tdf_t), intent(in) :: f
    tdf_niter = f%niter
  end function tdf_niter
  !------------------------------------------------------------


  !------------------------------------------------------------
  integer pure function tdf_nfreqs(f)
    type(tdf_t), intent(in) :: f
    tdf_nfreqs = f%nfreqs
  end function tdf_nfreqs
  !------------------------------------------------------------


  !------------------------------------------------------------
  FLOAT pure function tdf_dt(f)
    type(tdf_t), intent(in) :: f
    tdf_dt = f%dt
  end function tdf_dt
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init(f)
    type(tdf_t), intent(inout) :: f

    PUSH_SUB(tdf_init)

    f%mode = TDF_EMPTY
    f%niter = 0
    f%dt = M_ZERO

    POP_SUB(tdf_init)
  end subroutine tdf_init
  !------------------------------------------------------------


  !------------------------------------------------------------
  logical function tdf_is_empty(f)
    type(tdf_t), intent(in) :: f

    PUSH_SUB(tdf_is_empty)
    tdf_is_empty = (f%mode == TDF_EMPTY)

    POP_SUB(tdf_is_empty)
  end function tdf_is_empty
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_cw(f, a0, omega0)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, omega0

    PUSH_SUB(tdf_init_cw)

    f%mode = TDF_CW
    f%a0 = a0
    f%omega0 = omega0

    POP_SUB(tdf_init_cw)
  end subroutine tdf_init_cw
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_gaussian(f, a0, omega0, t0, tau0)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, omega0, t0, tau0

    PUSH_SUB(tdf_init_gaussian)

    f%mode = TDF_GAUSSIAN
    f%a0 = a0
    f%omega0 = omega0
    f%t0 = t0
    f%tau0 = tau0

    POP_SUB(tdf_init_gaussian)
  end subroutine tdf_init_gaussian
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_cosinoidal(f, a0, omega0, t0, tau0)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, omega0, t0, tau0

    PUSH_SUB(tdf_init_cosinoidal)

    f%mode = TDF_COSINOIDAL
    f%a0 = a0
    f%omega0 = omega0
    f%t0 = t0
    f%tau0 = tau0

    POP_SUB(tdf_init_cosinoidal)
  end subroutine tdf_init_cosinoidal
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_trapezoidal(f, a0, omega0, t0, tau0, tau1)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, omega0, t0, tau0, tau1

    PUSH_SUB(tdf_init_trapezoidal)

    f%mode = TDF_TRAPEZOIDAL
    f%a0 = a0
    f%omega0 = omega0
    f%t0 = t0
    f%tau0 = tau0
    f%tau1 = tau1

    POP_SUB(tdf_init_trapezoidal)
  end subroutine tdf_init_trapezoidal
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_fromexpr(f, expression)
    type(tdf_t),      intent(inout) :: f
    character(len=*), intent(in)    :: expression

    PUSH_SUB(tdf_init_fromexpr)

    f%mode = TDF_FROM_EXPR
    f%expression = trim(expression)

    POP_SUB(tdf_init_fromexpr)
  end subroutine tdf_init_fromexpr
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_fromfile(f, filename, namespace, ierr)
    type(tdf_t),      intent(inout) :: f
    character(len=*), intent(in)    :: filename
    type(namespace_t),intent(in)    :: namespace
    integer,          intent(out)   :: ierr

    integer :: iunit, lines, j
    FLOAT :: dummy
    FLOAT, allocatable :: t(:), am(:)

    PUSH_SUB(tdf_init_fromfile)

    f%mode = TDF_FROM_FILE
    ierr = 0

    iunit = io_open(trim(filename), namespace, action='read', status='old')

    ! count lines in file
    call io_skip_header(iunit)
    lines = 0
    do
      read(iunit, *, err=100, end=100) dummy, dummy
      lines = lines + 1
    end do
100 continue
    rewind(iunit)
    call io_skip_header(iunit)

    ! allocate and read info
    SAFE_ALLOCATE( t(1:lines))
    SAFE_ALLOCATE(am(1:lines))
    do j = 1, lines
      read(iunit, *) t(j), am(j)
    end do
    call io_close(iunit)

    f%init_time  = t(1)
    f%final_time = t(lines)

    call spline_init(f%amplitude)
    call spline_fit(lines, t, am, f%amplitude)

    SAFE_DEALLOCATE_A(t)
    SAFE_DEALLOCATE_A(am)

    POP_SUB(tdf_init_fromfile)
  end subroutine tdf_init_fromfile
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_numerical(f, niter, dt, omegamax, initval, rep)
    type(tdf_t),       intent(inout) :: f
    integer,           intent(in)    :: niter
    FLOAT,             intent(in)    :: dt
    FLOAT,             intent(in)    :: omegamax
    FLOAT,   optional, intent(in)    :: initval
    integer, optional, intent(in)    :: rep

    integer :: n(3), optimize_parity(3)
    logical :: optimize(3)
    FLOAT :: bigt

    PUSH_SUB(tdf_init_numerical)

    f%mode = TDF_NUMERICAL
    f%niter = niter
    SAFE_ALLOCATE(f%val(1:niter+1))
    if (present(initval)) then
      f%val = initval
    else
      f%val = M_ZERO
    end if
    f%dt = dt

    f%init_time = M_ZERO
    f%final_time = f%dt * f%niter

    if (omegamax > M_ZERO) then
      bigt = f%final_time - f%init_time
      f%nfreqs = int(bigt * omegamax / (M_TWO * M_PI)) + 1
    else
      f%nfreqs = f%niter/2+1
    end if

    n(1:3) = (/ f%niter, 1, 1 /)
    optimize(1:3) = .false.
    optimize_parity(1:3) = -1
    call fft_init(f%fft_handler, n, 1, FFT_REAL, FFTLIB_FFTW, optimize, optimize_parity)

    if (present(rep)) then
      select case (rep)
      case (TDF_FOURIER_SERIES,TDF_ZERO_FOURIER)
        SAFE_ALLOCATE(f%valww(1:2*(f%nfreqs-1)+1))
        f%valww = M_ZERO
        SAFE_DEALLOCATE_A(f%val)
        f%mode = rep
      end select
    end if

    POP_SUB(tdf_init_numerical)
  end subroutine tdf_init_numerical
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_fourier_grid(f, wgrid)
    type(tdf_t), intent(in)    :: f
    FLOAT,       intent(inout) :: wgrid(:)

    integer :: i
    FLOAT   :: df

    PUSH_SUB(tdf_fourier_grid)

    wgrid = M_ZERO
    select case (f%mode)
    case (TDF_FOURIER_SERIES, TDF_ZERO_FOURIER)
      df = M_TWO * M_PI / (f%final_time-f%init_time)
      do i = 1, f%nfreqs
        wgrid(i) = (i-1)*df
      end do
    case default
      message(1) = "Illegal mode in tdf_fourier_grid."
      call messages_fatal(1)
    end select

    POP_SUB(tdf_fourier_grid)
  end subroutine tdf_fourier_grid
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_numerical_to_fourier(f)
    type(tdf_t), intent(inout) :: f

    integer :: j
    CMPLX, allocatable :: tmp(:)

    PUSH_SUB(tdf_numerical_to_fourier)

    SAFE_ALLOCATE(tmp(1:f%niter/2+1))
    ! Moving to Fourier space implies periodic functions. However, the function on entry
    ! may not be periodic (f%val(f%niter +1) may not be equal to f%val(1)), so it is
    ! better if we take the average of those two values.
    f%val(1) = M_HALF*(f%val(1)+f%val(f%niter+1))
    f%val(f%niter+1) = f%val(1)
    call dfft_forward(f%fft_handler, f%val(1:f%niter), tmp)
    tmp = tmp * f%dt * sqrt(M_ONE/(f%final_time-f%init_time))
    f%mode = TDF_FOURIER_SERIES
    SAFE_ALLOCATE(f%valww(1:2*(f%nfreqs-1)+1))
    f%valww(1) = TOFLOAT(tmp(1))
    do j = 2, f%nfreqs
      f%valww(j) = (sqrt(M_TWO)) * TOFLOAT(tmp(j))
    end do
    do j = f%nfreqs+1, 2*f%nfreqs-1
      f%valww(j) = - (sqrt(M_TWO)) * aimag(tmp(j-f%nfreqs+1))
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(f%val)

    POP_SUB(tdf_numerical_to_fourier)
  end subroutine tdf_numerical_to_fourier
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_fourier_to_numerical(f)
    type(tdf_t), intent(inout) :: f

    integer :: j
    CMPLX, allocatable :: tmp(:)

    PUSH_SUB(tdf_fourier_to_numerical)

    SAFE_ALLOCATE(tmp(1:f%niter/2+1))
    tmp = M_z0
    tmp(1) = f%valww(1)
    do j = 2, f%nfreqs
      tmp(j) = TOCMPLX((sqrt(M_TWO)/M_TWO)*f%valww(j), -(sqrt(M_TWO)/M_TWO)*f%valww(j+f%nfreqs-1))
    end do
    SAFE_ALLOCATE(f%val(1:f%niter+1))
    call dfft_backward(f%fft_handler, tmp, f%val(1:f%niter))
    f%val(f%niter+1) = f%val(1)
    f%val = f%val * f%niter * sqrt(M_ONE/(f%final_time-f%init_time))
    f%mode = TDF_NUMERICAL

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(f%valww)

    POP_SUB(tdf_fourier_to_numerical)
  end subroutine tdf_fourier_to_numerical
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_numerical_to_zerofourier(f)
    type(tdf_t), intent(inout) :: f

    FLOAT :: s
    PUSH_SUB(tdf_numerical_to_zerofourier)

    call tdf_numerical_to_fourier(f)
    f%valww(1) = M_ZERO
    s = sum(f%valww(2:f%nfreqs))
    f%valww(2:f%nfreqs) = f%valww(2:f%nfreqs) - s/f%nfreqs
    f%mode = TDF_ZERO_FOURIER

    POP_SUB(tdf_numerical_to_zerofourier)
  end subroutine tdf_numerical_to_zerofourier
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_zerofourier_to_numerical(f)
    type(tdf_t), intent(inout) :: f

    PUSH_SUB(tdf_zerofourier_to_numerical)

    ASSERT(abs(f%valww(1)) <= M_EPSILON)
    call tdf_fourier_to_numerical(f)

    POP_SUB(tdf_zerofourier_to_numerical)
  end subroutine tdf_zerofourier_to_numerical
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_set_numericalr(f, values)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in) :: values(:)

    ! no push_sub because it is called too frequently

    select case (f%mode)
    case (TDF_NUMERICAL)
      f%val(1:f%niter+1) = values(1:f%niter+1)
    case (TDF_FOURIER_SERIES)
      f%valww(1:2*f%nfreqs-1) = values(1:2*f%nfreqs-1)
    case (TDF_ZERO_FOURIER)
      f%valww(1) = M_ZERO
      f%valww(2:2*f%nfreqs-1) = values(1:2*f%nfreqs-2)
    end select

  end subroutine tdf_set_numericalr
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_set_numericalr1(f, index, val)
    type(tdf_t), intent(inout) :: f
    integer,     intent(in)    :: index
    FLOAT,       intent(in)    :: val

    ! no push_sub because it is called too frequently

    select case (f%mode)
    case (TDF_NUMERICAL)
      f%val(index) = val
    case (TDF_FOURIER_SERIES)
      f%valww(index) = val
    case (TDF_ZERO_FOURIER)
      f%valww(index+1) = val
    end select

  end subroutine tdf_set_numericalr1
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_set_random(f, fdotf)
    type(tdf_t),     intent(inout) :: f
    FLOAT, optional, intent(in)    :: fdotf

    type(c_ptr) :: random_gen_pointer
    integer :: i, n
    FLOAT :: fdotf_, nrm
    FLOAT, allocatable :: e(:)

    PUSH_SUB(tdf_set_random)

    select case (f%mode)
    case (TDF_FOURIER_SERIES)
      n = 2*f%nfreqs-1
    case (TDF_ZERO_FOURIER)
      n = 2*f%nfreqs-2
    case default
      message(1) = "Illegal value for f%mode in tdf_set_random."
      call messages_fatal(1)
    end select
    SAFE_ALLOCATE(e(1:n))

    if (mpi_grp_is_root(mpi_world)) then
      if (present(fdotf)) then
        fdotf_ = fdotf
      else
        fdotf_ = tdf_dot_product(f, f)
      end if

      call loct_ran_init(random_gen_pointer)

      do i = 1, n
        e(i) = loct_ran_gaussian(random_gen_pointer, M_ONE)
      end do
      nrm = norm2(e)
      e = sqrt(fdotf_) * e/ nrm

      if (f%mode == TDF_ZERO_FOURIER) then
        e(1:f%nfreqs-1) = e(1:f%nfreqs-1) - sum(e(1:f%nfreqs-1))/(f%nfreqs-1)
        nrm = norm2(e)
        e = sqrt(fdotf_) * e/ nrm
      end if

      call loct_ran_end(random_gen_pointer)
    end if

    call mpi_world%bcast(e(1), n, MPI_FLOAT, 0)

    call tdf_set_numerical(f, e)
    SAFE_DEALLOCATE_A(e)
    POP_SUB(tdf_set_random)
  end subroutine tdf_set_random
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_to_numerical(f, niter, dt, omegamax)
    type(tdf_t),       intent(inout) :: f
    integer, optional, intent(in)    :: niter
    FLOAT,   optional, intent(in)    :: dt
    FLOAT,   optional, intent(in)    :: omegamax

    FLOAT :: t
    integer :: j
    FLOAT, allocatable :: val(:)

    if (f%mode == TDF_NUMERICAL) return
    PUSH_SUB(tdf_to_numerical)

    select case (f%mode)
    case (TDF_ZERO_FOURIER)
      call tdf_zerofourier_to_numerical(f)
    case (TDF_FOURIER_SERIES)
      call tdf_fourier_to_numerical(f)
    case default
      SAFE_ALLOCATE(val(1:niter+1))
      do j = 1, niter + 1
        t = (j-1)*dt
        val(j) = tdf(f, t)
      end do
      call tdf_end(f)
      ASSERT(present(niter))
      ASSERT(present(dt))
      ASSERT(present(omegamax))
      call tdf_init_numerical(f, niter, dt, omegamax)
      call tdf_set_numerical(f, val)
      SAFE_DEALLOCATE_A(val)
    end select

    POP_SUB(tdf_to_numerical)
  end subroutine tdf_to_numerical
  !------------------------------------------------------------


  !------------------------------------------------------------
  FLOAT pure function tdfi(f, i) result(y)
    type(tdf_t), intent(in) :: f
    integer,     intent(in) :: i

    ! Maybe there should be a grid for any kind of function, so
    ! that a meaningful number is produced in any case.
    y = M_ZERO
    select case (f%mode)
    case (TDF_NUMERICAL)
      y = f%val(i)
    case (TDF_FOURIER_SERIES)
      y = f%valww(i)
    case (TDF_ZERO_FOURIER)
      y = f%valww(i+1)
    end select

  end function tdfi
  !------------------------------------------------------------


  !------------------------------------------------------------
  FLOAT function tdft(f, t) result(y)
    type(tdf_t), intent(in) :: f
    FLOAT,       intent(in) :: t

    FLOAT, allocatable :: timearray(:), valarray(:)
    FLOAT :: r, fre, fim, tcu
    integer :: il, iu

    ! no push_sub because it is called too frequently

    select case (f%mode)

    case (TDF_CW)

      y = f%a0 * cos(f%omega0*t)

    case (TDF_GAUSSIAN)

      r = exp(-(t - f%t0)**2 / (M_TWO*f%tau0**2))
      y = f%a0 * r * cos(f%omega0 * t)

    case (TDF_COSINOIDAL)

      r = M_ZERO
      if (abs(t - f%t0) <= f%tau0) then
        r = cos((M_Pi / 2) * ((t - 2 * f%tau0 - f%t0) / f%tau0))
      end if
      y = f%a0 * r * cos(f%omega0 * t)

    case (TDF_TRAPEZOIDAL)
      if (t > f%t0-f%tau0/M_TWO-f%tau1 .and. t <= f%t0-f%tau0 / M_TWO) then
        r = (t - (f%t0 - f%tau0/M_TWO - f%tau1)) / f%tau1
      else if (t>f%t0-f%tau0/M_TWO .and. t <= f%t0+f%tau0 / M_TWO) then
        r = M_ONE
      else if (t>f%t0+f%tau0/M_TWO .and. t <= f%t0+f%tau0 / M_TWO+f%tau1) then
        r = (f%t0 + f%tau0/M_TWO + f%tau1 - t) / f%tau1
      else
        r = M_ZERO
      end if
      y = f%a0 * r * cos(f%omega0 * t)
    case (TDF_FROM_FILE)

      if (t >= f%init_time .and. t <= f%final_time) then
        y = spline_eval(f%amplitude, t)
      else
        y = M_ZERO
      end if

    case (TDF_NUMERICAL)

      il = int(t/f%dt)+1; iu = il+1

      SAFE_ALLOCATE(timearray(1:4))
      SAFE_ALLOCATE(valarray(1:4))

      if (il <= 1) then
        timearray = (/ M_ZERO, f%dt, M_TWO*f%dt, M_THREE*f%dt  /)
        valarray =  (/ f%val(1), f%val(2), f%val(3), f%val(4) /)
      elseif (il >= f%niter) then
        timearray = (/ (f%niter-3)*f%dt, (f%niter-2)*f%dt, (f%niter-1)*f%dt, f%niter*f%dt  /)
        valarray  = (/ f%val(f%niter-2), f%val(f%niter-1), f%val(f%niter), f%val(f%niter+1) /)
      else
        timearray = (/ (il-2)*f%dt, (il-1)*f%dt, il*f%dt, (il+1)*f%dt  /)
        valarray =  (/ f%val(il-1), f%val(il), f%val(il+1), f%val(il+2) /)
      end if

      call interpolate(timearray, valarray, t, y)

      SAFE_DEALLOCATE_A(valarray)
      SAFE_DEALLOCATE_A(timearray)

    case (TDF_FROM_EXPR)
      tcu = units_from_atomic(units_inp%time, t)
      call parse_expression(fre, fim, 't', tcu, f%expression)
      y = fre

    case default
      y = M_ZERO

    end select

  end function tdft
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_end(f)
    type(tdf_t), intent(inout) :: f

    PUSH_SUB(tdf_end)

    select case (f%mode)
    case (TDF_FROM_FILE)
      call spline_end(f%amplitude)
    case (TDF_NUMERICAL)
      call fft_end(f%fft_handler)
    case (TDF_FOURIER_SERIES, TDF_ZERO_FOURIER)
      SAFE_DEALLOCATE_A(f%valww)
      call fft_end(f%fft_handler)
    end select
    f%mode = TDF_EMPTY
    SAFE_DEALLOCATE_A(f%val)

    POP_SUB(tdf_end)
  end subroutine tdf_end
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_copy(fout, fin)
    type(tdf_t), intent(inout) :: fout
    type(tdf_t), intent(in)    :: fin

    PUSH_SUB(tdf_copy)

    ASSERT((fin%mode >= TDF_EMPTY) .and. (fin%mode <= TDF_ZERO_FOURIER))

    call tdf_end(fout)
    call tdf_init(fout)

    fout%t0     = fin%t0
    fout%tau0   = fin%tau0
    fout%tau1   = fin%tau1
    fout%dt     = fin%dt
    fout%a0     = fin%a0
    fout%omega0 = fin%omega0
    fout%niter  = fin%niter
    fout%final_time = fin%final_time
    fout%init_time  = fin%init_time
    fout%expression = fin%expression
    fout%nfreqs = fin%nfreqs
    if (fin%mode == TDF_FROM_FILE) then
      fout%amplitude = fin%amplitude
    end if
    if (fin%mode == TDF_NUMERICAL) then
      SAFE_ALLOCATE(fout%val(1:fout%niter+1))
      fout%val  = fin%val
      call fft_copy(fin%fft_handler, fout%fft_handler)
    end if
    if (fin%mode == TDF_FOURIER_SERIES .or. fin%mode == TDF_ZERO_FOURIER) then
      SAFE_ALLOCATE(fout%valww(1:2*fout%nfreqs-1))
      fout%valww  = fin%valww
      call fft_copy(fin%fft_handler, fout%fft_handler)
    end if
    fout%mode   = fin%mode

    POP_SUB(tdf_copy)
  end subroutine tdf_copy
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_scalar_multiply(alpha, f)
    FLOAT,       intent(in)    :: alpha
    type(tdf_t), intent(inout) :: f

    PUSH_SUB(tdf_scalar_multiply)

    select case (f%mode)
    case (TDF_CW, TDF_GAUSSIAN, TDF_COSINOIDAL, TDF_TRAPEZOIDAL)
      f%a0 = alpha*f%a0
    case (TDF_NUMERICAL)
      f%val = alpha*f%val
    case (TDF_FOURIER_SERIES,TDF_ZERO_FOURIER)
      f%valww = alpha*f%valww
    case (TDF_FROM_FILE)
      call spline_times(alpha, f%amplitude)
    end select

    POP_SUB(tdf_scalar_multiply)
  end subroutine tdf_scalar_multiply
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_cosine_multiply(omega, f)
    FLOAT,       intent(in)    :: omega
    type(tdf_t), intent(inout) :: f

    integer :: j
    FLOAT :: t

    PUSH_SUB(tdf_cosine_multiply)

    ! For the moment, we will just assume that f and g are of the same type.
    ASSERT(f%mode == TDF_NUMERICAL)

    do j = 1, f%niter + 1
      t = f%init_time + (j-1)*f%dt
      f%val(j) = f%val(j) * cos(omega*t)
    end do

    POP_SUB(tdf_cosine_multiply)
  end subroutine tdf_cosine_multiply
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_write(f, iunit, namespace)
    type(tdf_t),                 intent(in) :: f
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer :: n_msg

    PUSH_SUB(tdf_write)

    select case (f%mode)
    case (TDF_CW)
      write(message(1),'(6x,a)')          'Mode: continuous wave.'
      write(message(2),'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      n_msg = 2
    case (TDF_GAUSSIAN)
      write(message(1),'(6x,a)')          'Mode: Gaussian envelope.'
      write(message(2),'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(message(3),'(6x,a,f10.4,3a)') 'Width:     ', units_from_atomic(units_out%time, f%tau0), &
        ' [', trim(units_abbrev(units_out%time)), ']'
      write(message(4),'(6x,a,f10.4,3a)') 'Middle t:  ', units_from_atomic(units_out%time, f%t0), &
        ' [', trim(units_abbrev(units_out%time)), ']'
      n_msg = 4
    case (TDF_COSINOIDAL)
      write(message(1),'(6x,a)') 'Mode: cosinoidal envelope.'
      write(message(2),'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(message(3),'(6x,a,f10.4,3a)') 'Width:     ', units_from_atomic(units_out%time, f%tau0), &
        ' [', trim(units_abbrev(units_out%time)), ']'
      write(message(4),'(6x,a,f10.4,3a)') 'Middle t:  ', units_from_atomic(units_out%time, f%t0), &
        ' [', trim(units_abbrev(units_out%time)), ']'
      n_msg = 4
    case (TDF_TRAPEZOIDAL)
      write(message(1),'(6x,a)') 'Mode: trapezoidal envelope.'
      write(message(2),'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(message(3),'(6x,a,f10.4,3a)') 'Width:     ', units_from_atomic(units_out%time, f%tau0), &
        ' [', trim(units_abbrev(units_out%time)), ']'
      write(message(4),'(6x,a,f10.4,3a)') 'Middle t:  ', units_from_atomic(units_out%time, f%t0), &
        ' [', trim(units_abbrev(units_out%time)), ']'
      write(message(5),'(6x,a,f10.4,3a)') 'Ramp time: ', units_from_atomic(units_out%time, f%tau1), &
        ' [', trim(units_abbrev(units_out%time)), ']'
      n_msg = 5
    case (TDF_FROM_FILE)
      write(message(1),'(6x,a)') 'Mode: time-dependent function read from file.'
      n_msg = 1
    case (TDF_NUMERICAL)
      write(message(1),'(6x,a)') 'Mode: time-dependent function stored in a numerical array.'
      n_msg = 1
    case (TDF_FROM_EXPR)
      write(message(1),'(6x,a)') 'Mode: time-dependent function parsed from the expression:'
      write(message(2),'(6x,a)') '      f(t) = '//trim(f%expression)
      n_msg = 2
    end select

    if (f%omega0 /= M_ZERO) then
      n_msg = n_msg + 1
      write(message(n_msg),'(6x,a,f10.4,3a)') 'Frequency: ', units_from_atomic(units_out%energy, f%omega0), &
        ' [', trim(units_abbrev(units_out%energy)), ']'
    end if

    call messages_info(n_msg, iunit=iunit, namespace=namespace)

    POP_SUB(tdf_write)
  end subroutine tdf_write
  !------------------------------------------------------------


  !------------------------------------------------------------
  ! Returns the dot product of f and g, defined as:
  !    < f | g > = \int_0^T dt f^*(t) g(t)
  ! It assumes that both f and m are in the same mode, otherwise
  ! it will fail and stop the code.
  !------------------------------------------------------------
  FLOAT function tdf_dot_product(f, g) result (fg)
    type(tdf_t), intent(in) :: f, g

    integer :: i
    FLOAT :: t

    PUSH_SUB(tdf_dot_product)

    fg = M_ZERO

    ! For the moment, we will just assume that f and g are of the same type.
    ASSERT(f%mode == g%mode)

    select case (f%mode)
    case (TDF_NUMERICAL)
      ! We assume that the grid is the same for both functions.
      ! We will apply Simpson`s rule. However, note that this rule is only valid if there is an even
      ! number of orbitals. So if there is an odd number, we will apply a correction term (that will
      ! reduce the error order from 4 to 2, similar to the simple trapezoidal rule).
      fg = M_ZERO
      do i = 1, f%niter/2
        fg = fg + f%val(2*i-2+1)*g%val(2*i-2+1) + M_FOUR*f%val(2*i-1+1)*g%val(2*i-1+1) + f%val(2*i+1)*g%val(2*i+1)
      end do
      fg = fg * f%dt / M_THREE
      ! This is the correction term.
      if (mod(f%niter, 2).eq.1) then
        fg = fg + M_HALF * (f%val(f%niter)*g%val(f%niter) + f%val(f%niter+1) * g%val(f%niter+1)) * f%dt
      end if

    case (TDF_FOURIER_SERIES)
      fg = dot_product(f%valww, g%valww)
    case (TDF_ZERO_FOURIER)
      ASSERT(abs(f%valww(1)) <= M_EPSILON)
      fg = dot_product(f%valww, g%valww)
    case default
      do i = 1, f%niter + 1
        t = (i-1) * f%dt
        fg = fg + tdf(f, i) * tdf(g, i)
      end do

    end select

    POP_SUB(tdf_dot_product)
  end function tdf_dot_product
  !------------------------------------------------------------


  !------------------------------------------------------------
  ! Returns the difference of f and g, defined as:
  !    < f-g | f-g > = \int_0^T dt (f^*(t)-g^*(t))*(f(t)-g(t))
  ! It assumes that both f and m are in the same mode, otherwise
  ! it will fail and stop the code.
  !------------------------------------------------------------
  FLOAT function tdf_diff(f, g) result (fg)
    type(tdf_t), intent(in) :: f, g

    integer :: i
    type(tdf_t) :: fminusg

    PUSH_SUB(tdf_diff)

    ! For the moment, we will just assume that f and g are of the same type.
    ASSERT(f%mode == g%mode)

    call tdf_copy(fminusg, f)

    select case (f%mode)
    case (TDF_NUMERICAL)
      do i = 1, f%niter
        fminusg%val(i) = fminusg%val(i) - g%val(i)
      end do
    case (TDF_FOURIER_SERIES, TDF_ZERO_FOURIER)
      do i = 1, 2*(f%nfreqs-1)+1
        fminusg%valww(i) = fminusg%valww(i) - g%valww(i)
      end do
    case default
      message(1) = "Illegal value for f%mode in tdf_diff"
      call messages_fatal(1)
    end select

    fg = tdf_dot_product(fminusg, fminusg)

    call tdf_end(fminusg)

    POP_SUB(tdf_diff)
  end function tdf_diff
  !------------------------------------------------------------

end module tdfunction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
