!! Copyright (C) 2013 Umberto De Giovannini
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


!> The includes for the PNFFT
module pnfft_params_oct_m
  use,intrinsic :: iso_c_binding
  use pfft_params_oct_m
  implicit none
#ifdef HAVE_PNFFT
  include "pnfft.f03"
#endif
end module pnfft_params_oct_m


!> The low level module to work with the PNFFT library.
!! http://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
module pnfft_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use loct_math_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pfft_oct_m
  use pnfft_params_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::               &
    pnfft_t,              &
    pnfft_write_info,     &
    pnfft_init_plan,      &
    pnfft_copy_params,    &
    pnfft_init_procmesh,  &
    pnfft_end,            &
    pnfft_set_sp_nodes,   &
    pnfft_guru_options,   &
    zpnfft_forward,       &
    zpnfft_backward,      &
    dpnfft_forward,       &
    dpnfft_backward


  type pnfft_t
    private

! Parameters
    integer                       :: np(3)           !> Processes
    integer(C_INTPTR_T)           :: N_local(3)     !> Number of Fourier coefficients
    integer(C_INTPTR_T)           :: N(3)            !> Number of Fourier coefficients local
    integer(C_INTPTR_T)           :: Nos(3)          !> FFT grid size
    integer(C_INTPTR_T), public   :: M(3)
    integer                       :: M_istart(3)
    integer(C_INTPTR_T)           :: local_M         !> Local number of nodes per process
    integer                       :: mm              !> Real space cut-off
    FLOAT,               public   :: sigma           !> oversampling factor
    integer                       :: flags           !> PNFFT initialization options
    logical,             public   :: set_defaults = .false. !> set default values from the code

    FLOAT,               public   :: norm       !> Normalization


    integer               :: comm
! Data
    type(C_PTR)           :: plan            !> pnfft plan

    complex(C_DOUBLE_COMPLEX), pointer :: f_lin(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: f(:,:,:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: f_hat(:,:,:) => NULL()
    real(C_DOUBLE),            pointer :: x_lin(:,:) => NULL()
    real(C_DOUBLE),            pointer :: x(:,:,:,:) => NULL()

    real(C_DOUBLE)        :: lower_border(3) !> contains the real-space nodes local borders
    real(C_DOUBLE)        :: upper_border(3) !> parallelization
    FLOAT                 :: lo_global(3)
    FLOAT                 :: up_global(3)

  end type pnfft_t

contains

  ! ---------------------------------------------------------
  subroutine pnfft_guru_options(pnfft, namespace)
    type(pnfft_t),     intent(inout) :: pnfft
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(pnfft_guru_options)


    !%Variable PNFFTCutoff
    !%Type integer
    !%Default 6
    !%Section Mesh::FFTs
    !%Description
    !% Cut-off parameter of the window function.
    !%End
    call parse_variable(namespace, 'PNFFTCutoff', 6, pnfft%mm)

    !%Variable PNFFTOversampling
    !%Type float
    !%Default 2.0
    !%Section Mesh::FFTs
    !%Description
    !% PNFFT oversampling factor (sigma). This will rule the size of the FFT under the hood.
    !%End
    call parse_variable(namespace, 'PNFFTOversampling', M_TWO, pnfft%sigma)

    POP_SUB(pnfft_guru_options)
  end subroutine pnfft_guru_options

  ! ---------------------------------------------------------
  subroutine pnfft_init_params(pnfft, pnfft_options, nn, optimize)
    type(pnfft_t),     intent(inout) :: pnfft
    type(pnfft_t),     intent(in)    :: pnfft_options
    integer,           intent(in)    :: nn(3) !> pnfft bandwidths
    logical, optional, intent(in)    :: optimize

    integer :: ii, my_nn(3)
    logical :: optimize_

    PUSH_SUB(pnfft_init_params)

    optimize_ = optional_default(optimize, .true.)

    if (.not. pnfft%set_defaults) then
      !Set defaults
      pnfft%mm = pnfft_options%mm
      pnfft%sigma = pnfft_options%sigma
    end if

    my_nn = 0
    do ii = 1, 3
      my_nn(ii) = nn(ii)*int(pnfft%sigma)
      if (optimize_) call loct_fft_optimize(my_nn(ii), 1) ! ask for an odd number
    end do

    pnfft%Nos(1:3) = my_nn(1:3)

#ifdef HAVE_PNFFT
    pnfft%flags = PNFFT_MALLOC_X + PNFFT_MALLOC_F_HAT + PNFFT_MALLOC_F + &
      PNFFT_WINDOW_KAISER_BESSEL
#endif

    POP_SUB(pnfft_init_params)
  end subroutine pnfft_init_params


  ! ---------------------------------------------------------
  subroutine pnfft_init_procmesh(pnfft, mpi_grp, comm)
    type(pnfft_t), intent(inout)  :: pnfft
    type(mpi_grp_t),   intent(in) :: mpi_grp
    integer,           intent(out):: comm

    integer :: ierror

    PUSH_SUB(pnfft_init_procmesh)

#ifdef HAVE_PNFFT
    call pnfft_init()
#endif

    pnfft%np(1:3) = 1

    call pfft_decompose(mpi_grp%size, pnfft%np(1), pnfft%np(2))

#ifdef HAVE_PNFFT
    ierror = pnfft_create_procmesh(2, mpi_grp%comm,  pnfft%np, comm)
#else
    ierror = 0
    comm = 0
#endif

    if (ierror /= 0) then
      message(1) = "The number of rows and columns in PNFFT processor grid is not equal to "
      message(2) = "the number of processor in the MPI communicator."
      message(3) = "Please check it."
      call messages_fatal(3)
    end if

    POP_SUB(pnfft_init_procmesh)
  end subroutine pnfft_init_procmesh


  ! ---------------------------------------------------------
  subroutine pnfft_copy_params(in, out)
    type(pnfft_t), intent(in)  :: in
    type(pnfft_t), intent(out) :: out


    PUSH_SUB(pnfft_copy_params)

    out%mm = in%mm
    out%sigma = in%sigma
    out%set_defaults = in%set_defaults

    POP_SUB(pnfft_copy_params)
  end subroutine pnfft_copy_params

  ! ---------------------------------------------------------
  subroutine pnfft_write_info(pnfft)
    type(pnfft_t), intent(in) :: pnfft

    integer :: idir

    PUSH_SUB(pnfft_write_info)



    call messages_write("Info: PNFFT parameters")
    call messages_new_line()
    call messages_write("      Fourier coefficients      N = ")
    do idir = 1, 3
      call messages_write(pnfft%N(idir))
      if (idir < 3) call messages_write(" x ")
    end do
    call messages_new_line()
    call messages_write("      Spatial nodes per process   = ")
    call messages_write(pnfft%local_M)
    call messages_new_line()
    call messages_write("      Oversampling factor   sigma = ")
    call messages_write(pnfft%sigma)
    call messages_new_line()
    call messages_write("      FFT grid size             n = ")
    do idir = 1, 3
      call messages_write(pnfft%Nos(idir))
      if (idir < 3) call messages_write(" x ")
    end do
    call messages_new_line()
    call messages_write("      Real Space cutoff           = ")
    call messages_write(pnfft%mm)
    call messages_new_line()
    call messages_write("      Process mesh             np = ")
    do idir = 1, 3
      call messages_write(pnfft%np(idir))
      if (idir < 3) call messages_write(" x ")
    end do
    call messages_info()

    POP_SUB(pnfft_write_info)
  end subroutine pnfft_write_info

  ! ---------------------------------------------------------
  subroutine pnfft_init_plan(pnfft, pnfft_options, mpi_comm, fs_n_global, fs_n, fs_istart, rs_n, rs_istart)
    type(pnfft_t),   intent(inout) :: pnfft
    type(pnfft_t),   intent(inout) :: pnfft_options
    integer,         intent(in)    :: mpi_comm         !< MPI comunicator
    integer,         intent(in)    :: fs_n_global(1:3) !< The general number of elements in each dimension in Fourier space
    integer,         intent(out)   :: fs_n(1:3)        !< Local number of elements in each direction in Fourier space
    integer,         intent(out)   :: fs_istart(1:3)   !< Where does the local portion of the function start in Fourier space
    integer,         intent(out)   :: rs_n(1:3)        !< Local number of elements in each direction in real space
    integer,         intent(out)   :: rs_istart(1:3)   !< Where does the local portion of the function start in real space

    real(C_DOUBLE) :: lower_border(3), upper_border(3)
    real(C_DOUBLE) :: x_max(3)
    integer(C_INTPTR_T) :: local_N(3), local_N_start(3), d=3, local_M
    type(C_PTR) :: cf_hat, cf, cx


    PUSH_SUB(pnfft_init_plan)

    pnfft%N(1:3) = fs_n_global(1:3)

    call pnfft_init_params(pnfft, pnfft_options, fs_n_global(1:3), optimize = .true.)

    x_max(:) = CNST(0.4)

#ifdef HAVE_PNFFT
    call pnfft_local_size_guru(3, pnfft%N, pnfft%Nos, x_max, pnfft%mm, mpi_comm, &
      PNFFT_TRANSPOSED_F_HAT, local_N, local_N_start, lower_border, upper_border)
#else
    local_N = 0
    local_N_start = 0
    lower_border = M_ZERO
    upper_border = M_ZERO
#endif

    pnfft%comm = mpi_comm

    pnfft%lower_border = lower_border
    pnfft%upper_border = upper_border
    pnfft%N_local(1:3)   = local_N(1:3)

    pnfft%M(1)   = pnfft%N_local(2)
    pnfft%M(2)   = pnfft%N_local(3)
    pnfft%M(3)   = pnfft%N_local(1)

    local_M = pnfft%M(1) * pnfft%M(2) * pnfft%M(3)

    fs_n(1)      = int(local_N(1))
    fs_n(2)      = int(local_N(3))
    fs_n(3)      = int(local_N(2))
    fs_istart(1) = int(pnfft%N(1)/2 + local_N_start(1) + 1)
    fs_istart(2) = int(pnfft%N(3)/2 + local_N_start(3) + 1)
    fs_istart(3) = int(pnfft%N(2)/2 + local_N_start(2) + 1)


    rs_n(1:3) = int(pnfft%M(1:3))

    rs_istart(1) = fs_istart(3)
    rs_istart(2) = fs_istart(2)
    rs_istart(3) = fs_istart(1)

    pnfft%M_istart(:) = rs_istart(:)

#ifdef HAVE_PNFFT
    pnfft%plan = pnfft_init_guru(3, pnfft%N, pnfft%Nos, x_max, local_M, pnfft%mm, &
      pnfft%flags, PFFT_ESTIMATE, mpi_comm)
#endif

    pnfft%local_M=local_M

#ifdef HAVE_PNFFT
    ! Get data pointers in C format
    cf_hat = pnfft_get_f_hat(pnfft%plan)
    cf     = pnfft_get_f(pnfft%plan)
    cx     = pnfft_get_x(pnfft%plan)
#else
    cf_hat = C_NULL_PTR
    cf     = C_NULL_PTR
    cx     = C_NULL_PTR
#endif

    ! Convert data pointers to Fortran format
    call c_f_pointer(cf_hat, pnfft%f_hat, [local_N(1),local_N(3),local_N(2)])
    call c_f_pointer(cf,     pnfft%f_lin, [pnfft%local_M])
    call c_f_pointer(cf,     pnfft%f,     [pnfft%M(1),pnfft%M(2),pnfft%M(3)])
    call c_f_pointer(cx,     pnfft%x_lin, [d, pnfft%local_M])
    call c_f_pointer(cx,     pnfft%x,     [d, pnfft%M(1),pnfft%M(2),pnfft%M(3)])



    write(6,*) mpi_world%rank, "local_N(3)       ", local_N
    write(6,*) mpi_world%rank, "local_N_start(3) ", local_N_start
    write(6,*) mpi_world%rank, "fs_istart(1:3)   ", fs_istart
    write(6,*) mpi_world%rank, "fs_n(1:3)        ", fs_n
    write(6,*) mpi_world%rank, "rs_istart(1:3)   ", rs_istart
    write(6,*) mpi_world%rank, "rs_n(1:3)        ", rs_n
    write(6,*) mpi_world%rank, "lower_border     ", lower_border
    write(6,*) mpi_world%rank, "upper_border     ", upper_border
    write(6,*) mpi_world%rank, "rs_range         ", upper_border(:) - lower_border(:)
    write(6,*) mpi_world%rank, "local_M          ", local_M
    write(6,*) mpi_world%rank, "pnfft%N_local    ", pnfft%N_local
    write(6,*) mpi_world%rank, "pnfft%M          ", pnfft%M
    write(6,*) mpi_world%rank, "pnfft%M_istart   ", pnfft%M_istart
    write(6,*) mpi_world%rank, "size(pnfft%f_hat)", size(pnfft%f_hat,1), size(pnfft%f_hat,2), size(pnfft%f_hat, 3)
    write(6,*) mpi_world%rank, "size(pnfft%f)    ", size(pnfft%f,1), size(pnfft%f,2), size(pnfft%f,3)

    POP_SUB(pnfft_init_plan)
  end subroutine pnfft_init_plan

  ! ---------------------------------------------------------
  subroutine pnfft_end(pnfft)
    type(pnfft_t), intent(inout) :: pnfft

    PUSH_SUB(pnfft_end)

#ifdef HAVE_PNFFT
    call pnfft_finalize(pnfft%plan, PNFFT_FREE_X + PNFFT_FREE_F_HAT + PNFFT_FREE_F)
    call pnfft_cleanup()
#endif

    nullify(pnfft%f_lin)
    nullify(pnfft%f)
    nullify(pnfft%f_hat)
    nullify(pnfft%x)
    nullify(pnfft%x_lin)

    POP_SUB(pnfft_end)
  end subroutine pnfft_end


  ! ---------------------------------------------------------
  ! Set the coordinates for the spatial nodes rescaled to
  ! the 3D torus [-0.5,0.5)
  ! ---------------------------------------------------------
  subroutine pnfft_set_sp_nodes(pnfft, namespace, X)
    type(pnfft_t),    intent(inout) :: pnfft
    type(namespace_t),intent(in)    :: namespace
    FLOAT,            intent(in)    :: X(:,:) !X(i, dim)

    FLOAT   :: len(3), cc(3), eps,temp, lo(3), up(3)
    integer :: ii, idir, i1, i2, i3
    FLOAT, allocatable ::  dX(:,:)
!    integer :: j,t

    PUSH_SUB(pnfft_set_sp_nodes)

    eps = CNST(1.25) ! the sample nodes must be in [0.5,0.5)

    lo = pnfft%lower_border
    up = pnfft%upper_border



    do idir = 1,3
      len(:) = (maxval(X(:,idir))-minval(X(:,idir)))*eps
      cc(:) = (minval(X(:,idir))+maxval(X(:,idir)))/M_TWO
    end do




    do i1 = 1, int(pnfft%M(1))
      do i2 = 1, int(pnfft%M(2))
        do i3 = 1, int(pnfft%M(3))
          pnfft%x(1, i1,i2,i3) = (X(pnfft%M_istart(1)+i1-1, 1)  - cc(1))/len(1)
          pnfft%x(2, i1,i2,i3) = (X(pnfft%M_istart(2)+i2-1, 2)  - cc(2))/len(2)
          pnfft%x(3, i1,i2,i3) = (X(pnfft%M_istart(3)+i3-1, 3)  - cc(3))/len(3)

!           pnfft%x_lin(1, pnfft_idx_3to1(pnfft,i1,i2,i3)) = real((X(rs_istart(1)+i1-1, 1)  - cc(1))/len(1), C_DOUBLE)
!           pnfft%x_lin(2, pnfft_idx_3to1(pnfft,i1,i2,i3)) = real((X(rs_istart(2)+i2-1, 2)  - cc(2))/len(2), C_DOUBLE)
!           pnfft%x_lin(3, pnfft_idx_3to1(pnfft,i1,i2,i3)) = real((X(rs_istart(3)+i3-1, 3)  - cc(3))/len(3), C_DOUBLE)

!           pnfft%x_lin(1, pnfft_idx_3to1(pnfft,i1,i2,i3)) = &
!               (pnfft%upper_border(1) - pnfft%lower_border(1)) * rand(0) + pnfft%lower_border(1)
!           pnfft%x_lin(2, pnfft_idx_3to1(pnfft,i1,i2,i3)) = &
!               (pnfft%upper_border(2) - pnfft%lower_border(2)) * rand(0) + pnfft%lower_border(2)
!           pnfft%x_lin(3, pnfft_idx_3to1(pnfft,i1,i2,i3)) = &
!               (pnfft%upper_border(3) - pnfft%lower_border(3)) * rand(0) + pnfft%lower_border(3)

!           pnfft%x(1, i1,i2,i3) = (pnfft%upper_border(1) - pnfft%lower_border(1)) * rand(0) + pnfft%lower_border(1)
!           pnfft%x(2, i1,i2,i3) = (pnfft%upper_border(2) - pnfft%lower_border(2)) * rand(0) + pnfft%lower_border(2)
!           pnfft%x(3, i1,i2,i3) = (pnfft%upper_border(3) - pnfft%lower_border(3)) * rand(0) + pnfft%lower_border(3)

          temp = (X(pnfft%M_istart(3)+i3-1, 3)  - cc(3))/len(3)
          if (temp > pnfft%upper_border(3) .or. temp < pnfft%lower_border(3)) then
            write(6,*) mpi_world%rank, "out of bounds x3 = ", temp,"-- ", pnfft%lower_border(3), pnfft%upper_border(3)
          end if


        end do

        temp = (X(pnfft%M_istart(2)+i2-1, 2)  - cc(2))/len(2)
        if (temp > pnfft%upper_border(2) .or. temp < pnfft%lower_border(2)) then
          write(6,*) mpi_world%rank, "out of bounds x2 = ", temp,"-- ", pnfft%lower_border(2), pnfft%upper_border(2)
        end if

      end do

      temp = (X(pnfft%M_istart(1)+i1-1, 1)  - cc(1))/len(1)
      if (temp > pnfft%upper_border(1) .or. temp < pnfft%lower_border(1)) then
        write(6,*) mpi_world%rank, "out of bounds x1 = ", temp,"-- ", pnfft%lower_border(1), pnfft%upper_border(1)
      end if

    end do


!     do j = 1,pnfft%local_M
!       do t=1,3
!         pnfft%x_lin(t,j) = (pnfft%upper_border(t) - pnfft%lower_border(t)) * rand(0) + pnfft%lower_border(t)
!       end do
!     end do


    call pnfft_messages_debug(pnfft, namespace)


    SAFE_ALLOCATE( dX(1:maxval(pnfft%M(:))-1, 1:3))


    ! Set the normalization factor
    do idir = 1,3
      do ii = 1, size(X(:,idir))-1
        dX(ii,idir)= abs(X(ii+1, idir)-X(ii, idir))
!         dX(ii,2)= abs(x2_(ii+1)-x2_(ii))
!         dX(ii,3)= abs(x3_(ii+1)-x3_(ii))
      end do
    end do

!     pnfft%norm = M_ONE/(minval(dX(:,1)) * minval(dX(:,2)) * minval(dX(:,3)))
    pnfft%norm = M_ONE/(pnfft%N(1)*pnfft%N(2)*pnfft%N(3))

    write(6,*) mpi_world%rank, "pnfft%norm", pnfft%norm

    POP_SUB(pnfft_set_sp_nodes)
  end subroutine pnfft_set_sp_nodes

  ! ---------------------------------------------------------
  subroutine pnfft_messages_debug(pnfft, namespace)
    type(pnfft_t),     intent(in) :: pnfft
    type(namespace_t), intent(in) :: namespace

    integer          :: nn, i1, i2, i3
    integer          :: iunit          !< For debug output to files.
    character(len=3) :: filenum
#ifdef HAVE_MPI
    integer          :: ierr
#endif

    PUSH_SUB(pnfft_messages_debug)

    if (debug%info) then

      if (mpi_grp_is_root(mpi_world)) then
        call io_mkdir('debug/PNFFT', namespace)
      end if
#ifdef HAVE_MPI
      call MPI_Barrier(pnfft%comm, ierr)
#endif

      nn = mpi_world%rank
      write(filenum, '(i3.3)') nn

      iunit = io_open('debug/PNFFT/rs_partition.'//filenum, &
        namespace, action='write')

      do i1 = 1, int(pnfft%M(1))
        do i2 = 1, int(pnfft%M(2))
          do i3 = 1, int(pnfft%M(3))
            write(iunit, '(3f18.8)') pnfft%x(1, i1,i2,i3), pnfft%x(2, i1,i2,i3), pnfft%x(3, i1,i2,i3)
          end do
        end do
      end do
      call io_close(iunit)

    end if

    POP_SUB(pnfft_messages_debug)

  end subroutine pnfft_messages_debug

  !---------------------------------------------------------------------------------
!  integer function pnfft_idx_3to1(pnfft, ix , iy, iz) result(idx)
!    type(pnfft_t),  intent(in) :: pnfft
!    integer,        intent(in) :: ix
!    integer,        intent(in) :: iy
!    integer,        intent(in) :: iz
!
!    idx =  (ix-1)*pnfft%M(2)*pnfft%M(3) + (iy-1)*pnfft%M(3) + (iz-1) + 1
!
!  end function pnfft_idx_3to1

#include "undef.F90"
#include "real.F90"
#include "pnfft_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "pnfft_inc.F90"

end module pnfft_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
