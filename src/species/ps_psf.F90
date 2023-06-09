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

module ps_psf_oct_m
  use atomic_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use logrid_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use ps_in_grid_oct_m
  use ps_psf_file_oct_m

  implicit none

  private
  public ::     &
    ps_psf_t,       &
    ps_psf_init,    &
    ps_psf_end,     &
    ps_psf_process, &
    ps_psf_get_eigen

  type ps_psf_t
    ! Components are public by default
    type(ps_psf_file_t), private :: psf_file
    type(ps_in_grid_t)           :: ps_grid

    type(valconf_t)              :: conf
    FLOAT, allocatable,  private :: eigen(:, :)
    integer,             private :: ispin
  end type ps_psf_t

contains

  ! ---------------------------------------------------------
  subroutine ps_psf_init(pstm, ispin, filename, namespace)
    type(ps_psf_t),    intent(inout) :: pstm
    integer,           intent(in)    :: ispin
    character(len=*),  intent(in)    :: filename
    type(namespace_t), intent(in)    :: namespace

    character(len=MAX_PATH_LEN) :: fullpath
    integer :: iunit
    logical :: found, ascii

    PUSH_SUB(ps_psf_init)

    ! Sets the spin components
    pstm%ispin = ispin

    ! Find out where the hell the file is.

    found = .false.

    ASSERT(trim(filename) /= '')

    ascii = .true.
    fullpath = trim(filename)
    inquire(file = fullpath, exist = found)

    if (.not. found) then
      message(1) = "Pseudopotential file '" // trim(fullpath) // " not found"
      call messages_fatal(1, namespace=namespace)
    end if

    if (ascii) then
      iunit = io_open(fullpath, action='read', form='formatted', status='old')
    else
      iunit = io_open(fullpath, action='read', form='unformatted', status='old')
    end if
    call ps_psf_file_read(iunit, ascii, pstm%psf_file, namespace)
    call io_close(iunit)

    ! Fills the valence configuration data.
    call build_valconf(pstm%psf_file, ispin, pstm%conf, namespace)

    ! Hack
    if (mod(pstm%psf_file%nr, 2) == 0) pstm%psf_file%nr = pstm%psf_file%nr - 1

    call file_to_grid(pstm%psf_file, pstm%ps_grid)

    POP_SUB(ps_psf_init)
  end subroutine ps_psf_init


  ! ---------------------------------------------------------
  subroutine ps_psf_end(ps_psf)
    type(ps_psf_t), intent(inout) :: ps_psf

    PUSH_SUB(ps_psf_end)

    SAFE_DEALLOCATE_A(ps_psf%eigen)
    call ps_in_grid_end(ps_psf%ps_grid)
    call ps_psf_file_end(ps_psf%psf_file)

    POP_SUB(ps_psf_end)
  end subroutine ps_psf_end


  ! ---------------------------------------------------------
  subroutine build_valconf(psf_file, ispin, conf, namespace)
    type(ps_psf_file_t), intent(in)  :: psf_file
    integer,             intent(in)  :: ispin
    type(valconf_t),     intent(out) :: conf
    type(namespace_t),   intent(in)  :: namespace

    character(len=1)   :: char1(6), char2
    character(len=256) :: r_fmt
    integer :: l

    PUSH_SUB(build_valconf)

    conf%symbol = psf_file%namatm
    conf%p = psf_file%npotd
    write(char2,'(i1)') conf%p

    select case (psf_file%irel)
    case ('nrl')
      write(r_fmt, '(3a)') '(', char2, '(i1,a1,f5.2,10x))'
      read(psf_file%title, r_fmt)        &
        (conf%n(l), char1(l), conf%occ(l,1), l = 1, conf%p)
    case ('isp')
      write(r_fmt, '(3a)') '(', char2, '(i1,a1,f4.2,1x,f4.2,6x))'
      read(psf_file%title, r_fmt)        &
        (conf%n(l), char1(l), conf%occ(l,1), conf%occ(l,2), l = 1, conf%p)
    case ('rel')
      write(r_fmt, '(3a)') '(', char2, '(i1,a1,f5.2,10x))'
      read(psf_file%title, r_fmt)        &
        (conf%n(l), char1(l), conf%occ(l,1), l = 1, conf%p)
    end select

    do l = 1, conf%p
      select case (char1(l))
      case ('s')
        conf%l(l) = 0
      case ('p')
        conf%l(l) = 1
      case ('d')
        conf%l(l) = 2
      case ('f')
        conf%l(l) = 3
      case default
        message(1) = 'Error reading pseudopotential file.'
        call messages_fatal(1, namespace=namespace)
      end select
    end do

    if (ispin == 2 .and. psf_file%irel /= 'isp') then
      call valconf_unpolarized_to_polarized(conf)
    end if

    POP_SUB(build_valconf)
  end subroutine build_valconf

  !----------------------------------------------------------------
  subroutine file_to_grid(psf_file, ps_grid)
    type(ps_psf_file_t), intent(in)  :: psf_file
    type(ps_in_grid_t),  intent(out) :: ps_grid

    integer :: nrval

    PUSH_SUB(file_to_grid)

    ! Initializes the pseudo in the logarithmic grid.
    call ps_in_grid_init(ps_grid,                      &
      LOGRID_PSF, psf_file%a, psf_file%b, psf_file%nr,  &
      psf_file%npotd, psf_file%npotu)

    nrval = ps_grid%g%nrval

    ps_grid%zval      = psf_file%zval
    ps_grid%vps(1:nrval,:)  = psf_file%vps(1:nrval,:)
    ps_grid%chcore(1:nrval) = psf_file%chcore(1:nrval)
    if (ps_grid%so_no_l_channels > 0) then
      ps_grid%so_vps(1:nrval,:) = psf_file%vso(1:nrval,:)
    end if

    ps_grid%core_corrections = .true.
    if (trim(psf_file%icore) == 'nc') ps_grid%core_corrections = .false.

    POP_SUB(file_to_grid)
  end subroutine file_to_grid


  ! ---------------------------------------------------------
  subroutine ps_psf_process(ps_psf, namespace, lmax, lloc)
    type(ps_psf_t),    intent(inout) :: ps_psf
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: lmax, lloc

    PUSH_SUB(psf_process)

    ! get the pseudoatomic eigenfunctions
    SAFE_ALLOCATE(ps_psf%eigen(1:ps_psf%psf_file%npotd, 1:3))
    call solve_schroedinger(ps_psf%psf_file, namespace, ps_psf%ps_grid%g, &
      ps_psf%conf, ps_psf%ispin, ps_psf%ps_grid%rphi, ps_psf%eigen)

    ! check norm of rphi
    call ps_in_grid_check_rphi(ps_psf%ps_grid, namespace)

    ! Fix the local potential. Final argument is the core radius
    call ps_in_grid_vlocal(ps_psf%ps_grid, lloc, M_THREE, namespace)

    ! Calculate kb cosines and norms
    call ps_in_grid_kb_cosines(ps_psf%ps_grid, lloc)

    ! Ghost analysis.
    call ghost_analysis(ps_psf%psf_file, ps_psf%ps_grid, ps_psf%ps_grid%g, namespace, ps_psf%eigen, lmax)

    ! Define the KB-projector cut-off radii
    call ps_in_grid_cutoff_radii(ps_psf%ps_grid, lloc)

    ! Calculate KB-projectors
    call ps_in_grid_kb_projectors(ps_psf%ps_grid)

    POP_SUB(psf_process)
  end subroutine ps_psf_process

  ! ---------------------------------------------------------
  subroutine ps_psf_get_eigen(ps_psf, eigen)
    type(ps_psf_t), intent(in)  :: ps_psf
    FLOAT,          intent(out) :: eigen(:,:)

    integer :: i

    PUSH_SUB(ps_psf_get_eigen)

    do i = 1, ps_psf%conf%p
      if (ps_psf%ispin == 1) then
        eigen(i, 1) = ps_psf%eigen(i, 1)
      else
        eigen(i, 1:2) = ps_psf%eigen(i, 2:3)
      end if
    end do

    POP_SUB(ps_psf_get_eigen)
  end subroutine ps_psf_get_eigen

  ! ---------------------------------------------------------
  subroutine solve_schroedinger(psf_file, namespace, g, conf, ispin, rphi, eigen)
    type(ps_psf_file_t), intent(inout) :: psf_file
    type(namespace_t),   intent(in)    :: namespace
    type(logrid_t),      intent(in)    :: g
    type(valconf_t),     intent(in)    :: conf
    integer,             intent(in)    :: ispin
    FLOAT,               intent(out)   :: rphi(:,:,:), eigen(:,:)


    character(len=3) :: functl
    integer :: iter, ir, is, l, nnode, nprin, ierr
    FLOAT :: vtot, diff, a2b4
    FLOAT, allocatable :: ve(:, :), rho(:, :), prev(:, :)

    ! These variables are in double precision, even if the single-precision version of
    ! octopus is compiled, because they are passed to egofv.
    REAL_DOUBLE :: e, z, dr, rmax
    REAL_DOUBLE, allocatable :: s(:), hato(:), gg(:), y(:)

    PUSH_SUB(solve_schroedinger)

    ! Allocation.
    SAFE_ALLOCATE(s   (1:g%nrval))
    SAFE_ALLOCATE(hato(1:g%nrval))
    SAFE_ALLOCATE(gg  (1:g%nrval))
    SAFE_ALLOCATE(y   (1:g%nrval))
    SAFE_ALLOCATE(ve  (1:g%nrval, 1:ispin))
    SAFE_ALLOCATE(rho (1:g%nrval, 1:ispin))
    SAFE_ALLOCATE(prev(1:g%nrval, 1:ispin))
    s = M_ZERO
    hato = M_ZERO
    gg = M_ZERO
    y = M_ZERO
    ve = M_ZERO
    rho = M_ZERO
    prev = M_ZERO

    ! These numerical parameters have to be fixed for egofv to work.
    s(:) = g%drdi(:)**2
    a2b4 = M_FOURTH*psf_file%a**2

    !  ionic pseudopotential if core-correction for Hartree
    if ((psf_file%icore == 'pche') .or. (psf_file%icore == 'fche')) then
      call vhrtre(psf_file%chcore, ve(:, 1), g%rofi, g%drdi, g%s, g%nrval, g%a)
      do l = 1, psf_file%npotd
        psf_file%vps(:, l) = psf_file%vps(:, l) + ve(:, 1)
      end do
    end if

    ! Calculation of the valence screening potential from the density:
    !       ve(1:nrval) is the hartree+xc potential created by the pseudo -
    !               valence charge distribution (everything in Rydbergs, and Bohrs)
    rho(1:g%nrval, 1) = psf_file%rho_val(1:g%nrval)
    select case (psf_file%icorr)
    case ('pb')
      functl = 'GGA'
    case default
      functl = 'LDA'
    end select

    call atomhxc(namespace, functl, g, 1, rho(:, 1:1), ve(:, 1:1), psf_file%chcore)

    ! Calculation of the pseudo-wavefunctions.
    !       rphi(1:nrval, 1:conf%p, :) : radial pseudo-wavefunctions. They are normalized so that
    !           int dr {rphi^2} = 1. Thus its units are Bohr^(-1/2).
    !       eigen(1:conf%p, :)        : eigenvalues, in Rydbergs.
    do l = 1, psf_file%npotd
      do ir = 2, g%nrval
        vtot = psf_file%vps(ir, l) + ve(ir, 1) + TOFLOAT((l-1)*l)/(g%r2ofi(ir))
        hato(ir) = vtot*s(ir) + a2b4
      end do
      hato(1) = first_point_extrapolate(g%rofi, hato)

      nnode = 1
      nprin = l
      e     = -((psf_file%zval/TOFLOAT(nprin))**2)
      z     = psf_file%zval
      dr    = CNST(-1.0e5)
      rmax = g%rofi(g%nrval)

      call egofv(namespace, hato, s, g%nrval, e, gg, y, l-1, z, &
        TOFLOAT(g%a), TOFLOAT(g%b), rmax, nprin, nnode, dr, ierr)

      if (ierr /= 0) then
        write(message(1),'(a)') 'The algorithm that calculates atomic wavefunctions could not'
        write(message(2),'(a)') 'do its job. The program will terminate, since the wavefunctions'
        write(message(3),'(a)') 'are needed. Change the pseudopotential or improve the code.'
        call messages_fatal(3, namespace=namespace)
      end if
      eigen(l, 1) = e

      rphi(:, l, 1) = gg(:) * sqrt(g%drdi(:))

    end do

    ! Make a copy of the wavefunctions.
    rphi(:, :, 2) = rphi(:, :, 1)

    ! And now, for the spin-polarized case...
    spin_polarized: if (ispin == 2) then

      rho = M_ZERO    ! Here information from previous calculation could be used, but
      prev = M_ZERO   ! to save code lines, let us start from scratch.
      diff = CNST(1.0e5)
      iter = 0
      self_consistent: do
        prev = rho
        iter = iter + 1

        spin: do is = 1, ispin
          ang: do l = 1, psf_file%npotd
            do ir = 2, g%nrval
              vtot = psf_file%vps(ir, l) + ve(ir, is) + TOFLOAT((l-1)*l)/g%r2ofi(ir)
              hato(ir) = vtot*s(ir) + a2b4
            end do
            hato(1) = first_point_extrapolate(g%rofi, hato)

            nnode = 1
            nprin = l
            e     = -((psf_file%zval/TOFLOAT(nprin))**2)
            z     = psf_file%zval
            dr    = -CNST(1.0e5)
            rmax = g%rofi(g%nrval)

            call egofv(namespace, hato, s, g%nrval, e, gg, y, l-1, z, &
              TOFLOAT(g%a), TOFLOAT(g%b), rmax, nprin, nnode, dr, ierr)

            if (ierr /= 0) then
              write(message(1),'(a)') 'The algorithm that calculates atomic wavefunctions could not'
              write(message(2),'(a)') 'do its job. The program will terminate, since the wavefunctions'
              write(message(3),'(a)') 'are needed. Change the pseudopotential or improve the code.'
              call messages_fatal(3, namespace=namespace)
            end if
            eigen(l, 1 + is) = e

            rphi(:, l, 1 + is) = gg(:) * sqrt(g%drdi(:))
          end do ang
        end do spin

        rho = M_ZERO
        do is = 1, ispin
          do l = 1, psf_file%npotd
            rho(:, is) = rho(:, is) + conf%occ(l, is)*rphi(:, l, 1 + is)**2
          end do
        end do

        diff = M_ZERO
        do is = 1, ispin
          diff = diff + sqrt(sum(g%drdi(:) * (rho(:, is) - prev(:, is))**2))
        end do

        if (diff < M_EPSILON*CNST(1e2)) exit self_consistent
        if (iter > 1) rho = M_HALF * rho + M_HALF * prev

        !write(message(1),'(a,i4,a,e10.2)') '      Iter =', iter, '; Diff =', diff
        !call messages_info(1, namespace=namespace)

        call atomhxc(namespace, functl, g, 2, rho, ve, psf_file%chcore)
      end do self_consistent

    end if spin_polarized

    ! Exit this...
    SAFE_DEALLOCATE_A(s)
    SAFE_DEALLOCATE_A(hato)
    SAFE_DEALLOCATE_A(gg)
    SAFE_DEALLOCATE_A(y)
    SAFE_DEALLOCATE_A(ve)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(prev)

    POP_SUB(solve_schroedinger)
  end subroutine solve_schroedinger


  ! ---------------------------------------------------------
  subroutine ghost_analysis(psf_file, ps_grid, g, namespace, eigen, lmax)
    type(ps_psf_file_t), intent(in) :: psf_file
    type(ps_in_grid_t), intent(in) :: ps_grid
    type(logrid_t),     intent(in) :: g
    type(namespace_t),  intent(in) :: namespace
    FLOAT,              intent(in) :: eigen(:,:)
    integer,            intent(in) :: lmax

    character(len=3) :: functl
    integer :: ir, l, nnode, nprin, ighost, ierr
    FLOAT :: vtot, a2b4
    FLOAT, allocatable :: ve(:), elocal(:,:)

    REAL_DOUBLE :: z, e, dr, rmax
    REAL_DOUBLE, allocatable :: hato(:), s(:), gg(:), y(:)

    PUSH_SUB(ghost_analysis)

    SAFE_ALLOCATE(ve    (1:g%nrval))
    SAFE_ALLOCATE(s     (1:g%nrval))
    SAFE_ALLOCATE(hato  (1:g%nrval))
    SAFE_ALLOCATE(gg    (1:g%nrval))
    SAFE_ALLOCATE(y     (1:g%nrval))
    SAFE_ALLOCATE(elocal(1:2, 1:lmax+1))

    select case (psf_file%icorr)
    case ('pb')
      functl = 'GGA'
    case default
      functl = 'LDA'
    end select

    call atomhxc(namespace, functl, g, 1, psf_file%rho_val(:), ve(:), ps_grid%chcore(:))

    ! calculate eigenvalues of the local potential for ghost analysis
    s(:) = g%drdi(:)**2
    a2b4 = M_FOURTH*g%a**2

    do l = 1, lmax+1
      do ir = 2, g%nrval
        vtot = ps_grid%vlocal(ir) + ve(ir) + TOFLOAT((l-1)*l)/(g%r2ofi(ir))
        hato(ir) = vtot*s(ir) + a2b4
      end do
      hato(1) = first_point_extrapolate(g%rofi, hato)

      do nnode = 1, 2
        nprin = l
        e     = -(ps_grid%zval/TOFLOAT(nprin))**2
        z     = ps_grid%zval
        dr    = CNST(-1.0e5)
        rmax  = g%rofi(g%nrval)

        call egofv(namespace, hato, s, g%nrval, e, gg, y, l, z, &
          TOFLOAT(g%a), TOFLOAT(g%b), rmax, nprin, nnode, dr, ierr)

        elocal(nnode, l) = e
      end do
    end do

    ! Ghost analysis
    do l = 1, lmax+1
      ighost = -1

      if (ps_grid%dkbcos(l) > M_ZERO) then
        if (eigen(l, 1) > elocal(2, l)) then
          ighost = 1
        end if
      else if (ps_grid%dkbcos(l) < M_ZERO) then
        if (eigen(l, 1) > elocal(1, l)) then
          ighost = 1
        end if
      end if

      if (ighost >= 0) then
        write(message(1), '(a,i2)') "Ghost state found for l = ", l
        call messages_warning(1, namespace=namespace)
      end if
    end do

    SAFE_DEALLOCATE_A(hato)
    SAFE_DEALLOCATE_A(gg)
    SAFE_DEALLOCATE_A(y)
    SAFE_DEALLOCATE_A(elocal)
    SAFE_DEALLOCATE_A(s)
    SAFE_DEALLOCATE_A(ve)

    POP_SUB(ghost_analysis)
  end subroutine ghost_analysis

end module ps_psf_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
