!! Copyright (C) 2023. A. Buccheri
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

module chebyshev_filter_bounds_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use debug_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use lalg_basic_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use wfs_elec_oct_m

  implicit none
  private
  public :: &
    dupper_bound_estimator, &
    zupper_bound_estimator, &
    dfilter_bounds_estimator, &
    zfilter_bounds_estimator, &
    chebyshev_filter_bounds_t

  !> @class Chebyshev filter bounds
  !! @brief Chebyshev filter bounds
  type chebyshev_filter_bounds_t
    FLOAT :: lower         !< Lower bound of filter
    FLOAT :: upper         !< Upper bound of filter
    FLOAT :: lower_scaled  !< Scaled lower bound, which should approximate the lowest eigenvalue of H
  contains
    procedure :: center
    procedure :: half_width
    procedure :: sigma
    final :: finalize
  end type chebyshev_filter_bounds_t

  !> Overload default constructor
  interface chebyshev_filter_bounds_t
    procedure constructor
  end interface chebyshev_filter_bounds_t

contains

  !> @brief Create an instance of chebyshev_filter_bounds_t
  !!
  !! If \f$a_l\f$ is not passed, its value will default to
  !! the lower bound of the filter, which results in
  !! application of the simple Chebyshev filter.
  function constructor(lower, upper, a_l) result(this)
    FLOAT,           intent(in) :: lower, upper  !< Lower and upper filter bounds
    FLOAT, optional, intent(in) :: a_l           !< Smallest eigenvalue in spectrum, for scaling the filter
    class(chebyshev_filter_bounds_t), pointer :: this

    if(upper <= lower) then
      message(1) = "Chebychev filtering lower bound cannot be >= the upper bound."
      call messages_fatal(1)
    end if

    SAFE_ALLOCATE(this)
    this%lower = lower
    this%upper = upper

    this%lower_scaled = lower
    if (present(a_l)) then
      ASSERT(a_l < this%lower)
      this%lower_scaled = a_l
    endif

  end function constructor

  !> @brief Finalizer
  !!
  !! Implemented zeroing to avoid warning of stub routine
  subroutine finalize(this)
    type(chebyshev_filter_bounds_t) :: this

    this%lower = M_ZERO
    this%upper = M_ZERO
    this%lower_scaled = M_ZERO

  end subroutine finalize

  !> @brief Center of the filter interval.
  pure FLOAT function center(this)
    class(chebyshev_filter_bounds_t), intent(in) :: this

    center = M_HALF * (this%lower + this%upper)

  end function center

  !> @brief Half-width of the filter interval.
  !!
  !! Denoted as `e` in Algorithm 3.2 of
  !! [Zhou et. al.](http://dx.doi.org/10.1016/j.jcp.2014.06.056)
  pure FLOAT function half_width(this)
    class(chebyshev_filter_bounds_t), intent(in) :: this

    half_width = M_HALF * (this%upper - this%lower)

  end function half_width

  !> @brief Sigma scaling function, arising from application of the 3-term recurrence
  !> relationship of the Chebyshev polynomial.
  !!
  !! If \f$a_l\f$ is not defined in the constructor, sigma is set to one, and
  !! a simple Chebyshev filter \f$p(H)\f$ is used.
  !!
  !! If \f$a_l\f$ is defined, one applies a scaled Chebyshev filter,
  !! \f$ \tilde{p}(H)=p(H) / p \left( a_L \right) \f$.
  !!
  !! See equations 8 - 10 in [Zhou et. al.](http://dx.doi.org/10.1016/j.jcp.2014.06.056)
  !! for more details.
  FLOAT function sigma(this)
    class(chebyshev_filter_bounds_t), intent(in) :: this

    ASSERT(this%lower_scaled <= this%lower)
    if (is_close(this%lower, this%lower_scaled)) then
      sigma = M_ONE
    endif
    sigma = this%half_width() / (this%center() - this%lower_scaled)

  end function sigma

#include "real.F90"
#include "chebyshev_filter_bounds_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "chebyshev_filter_bounds_inc.F90"
#include "undef.F90"

end module chebyshev_filter_bounds_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
