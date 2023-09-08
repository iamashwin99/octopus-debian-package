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

!> This module implements the curvilinear coordinates given in
!! N. A. Modine, G. Zumbach, and E. Kaxiras, Phys. Rev. B 55, 10289-10301 (1997)
!!
!! The local refinement was changed for a simple exponential.
!! I believe that the recipe given by the authors is too complicated
!! for me to sort out.

module curv_modine_oct_m
  use coordinate_system_oct_m
  use debug_oct_m
  use global_oct_m
  use lalg_adv_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use root_solver_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                     &
    curv_modine_t,              &
    curv_modine_copy

  type, extends(coordinate_system_t) :: curv_modine_t
    private
    FLOAT, allocatable :: lsize(:)   !< size of the box
    FLOAT              :: xbar       !< size of central flat region (in units of L)
    FLOAT              :: Jbar       !< increase in density of points is 1/J

    FLOAT, allocatable :: Jlocal(:)  !< local (around the atoms) refinement
    FLOAT, allocatable :: Jrange(:)  !< local refinement range

    FLOAT, allocatable :: chi_atoms(:,:)
    FLOAT, allocatable :: csi(:,:)

    integer :: npos
    type(root_solver_t) :: rs
  contains
    procedure :: to_cartesian => curv_modine_to_cartesian
    procedure :: from_cartesian => curv_modine_from_cartesian
    procedure :: det_jac => curv_modine_det_jac
    procedure :: write_info => curv_modine_write_info
    procedure :: surface_element => curv_modine_surface_element
    final :: curv_modine_finalize
  end type curv_modine_t

  interface curv_modine_t
    procedure curv_modine_constructor
  end interface curv_modine_t

  ! Auxiliary variables for the root solver.
  class(curv_modine_t), pointer :: modine_p
  FLOAT, allocatable :: x_p(:)

contains

  ! ---------------------------------------------------------
  function curv_modine_constructor(namespace, dim, npos, pos, lsize, spacing) result(modine)
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(in)  :: dim
    integer,           intent(in)  :: npos
    FLOAT,             intent(in)  :: pos(1:dim,1:npos)
    FLOAT,             intent(in)  :: lsize(1:dim)
    FLOAT,             intent(in)  :: spacing(1:dim)
    class(curv_modine_t), pointer :: modine

    PUSH_SUB(curv_modine_constructor)

    SAFE_ALLOCATE(modine)

    modine%dim = dim
    modine%local_basis = .true.
    modine%orthogonal = .true. ! This needs to be checked.

    modine%npos = npos

    !%Variable CurvModineXBar
    !%Type float
    !%Default 1/3
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Size of central flat region (in units of <tt>Lsize</tt>). Must be between 0 and 1.
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineXBar', M_ONE/M_THREE, modine%xbar)

    !%Variable CurvModineJBar
    !%Type float
    !%Default 1/2
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Increase in density of points is inverse of this parameter.
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineJBar', M_HALF, modine%Jbar)

    SAFE_ALLOCATE(modine%lsize(1:dim))
    modine%lsize = lsize / modine%Jbar

    if (modine%xbar < M_ZERO .or. modine%xbar > M_ONE) then
      message(1) = 'The parameter "CurvModineXBar" must lie between 0 and 1.'
      call messages_fatal(1, namespace=namespace)
    end if

    SAFE_ALLOCATE(modine%Jlocal(1:npos))
    SAFE_ALLOCATE(modine%Jrange(1:npos))

    ! \warning: the reading has to be done for each atom kind

    !%Variable CurvModineJlocal
    !%Type float
    !%Default 0.25
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Local refinement around the atoms. Must be between 0 and 1.
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineJlocal', CNST(0.25), modine%Jlocal(1))

    !%Variable CurvModineJrange
    !%Type float
    !%Default 2 b
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Local refinement range (a length).
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineJrange', M_TWO, modine%Jrange(1), units_inp%length)

    if (modine%Jlocal(1)<M_ZERO.or.modine%Jlocal(1)>M_ONE) then
      message(1) = 'The parameter "CurvModineJlocal" must lie between 0 and 1.'
      call messages_fatal(1, namespace=namespace)
    end if

    modine%Jlocal(:) = modine%Jlocal(1)
    modine%Jrange(:) = modine%Jrange(1)

    ! initialize root solver for the optimization
    call root_solver_init(modine%rs, namespace, dim, &
      solver_type = ROOT_NEWTON, maxiter = 500, abs_tolerance = CNST(1.0e-10))

    call find_atom_points()
    call optimize()

    modine%min_mesh_scaling_product = modine%Jbar**modine%dim

    POP_SUB(curv_modine_constructor)
  contains

    subroutine find_atom_points()
      integer :: iatom, jj

      PUSH_SUB(curv_modine_constructor.find_atom_points)

      ! Initialize csi
      SAFE_ALLOCATE(modine%csi(1:modine%dim, 1:modine%npos))
      modine%csi = npos

      ! get first estimate for chi_atoms
      SAFE_ALLOCATE(modine%chi_atoms(1:modine%dim, 1:modine%npos))
      do jj = 1, 10  ! \warning: make something better
        do iatom = 1, modine%npos
          modine%chi_atoms(:, iatom) = modine%from_cartesian(pos(:, iatom))
        end do
        modine%csi(:,:) = modine%chi_atoms(:,:)
      end do

      do iatom = 1, modine%npos
        ! These are the chi positions where we want the atoms.
        modine%chi_atoms(:, iatom) = nint(modine%chi_atoms(:, iatom) / spacing(:)) * spacing(:)
      end do

      POP_SUB(curv_modine_constructor.find_atom_points)
    end subroutine find_atom_points

    subroutine optimize()
      logical :: conv
      integer :: iatom, idim, index
      FLOAT, allocatable :: my_csi(:), start_csi(:)

      PUSH_SUB(curv_modine_constructor.optimize)

      modine_p  => modine

      SAFE_ALLOCATE(x_p(1:modine%dim*modine%npos))
      SAFE_ALLOCATE(my_csi(1:modine%dim*modine%npos))
      SAFE_ALLOCATE(start_csi(1:modine%dim*modine%npos))

      do iatom = 1, modine%npos
        do idim = 1, modine%dim
          index = (iatom-1)*dim + idim
          x_p(index)       = pos(idim, iatom)
          start_csi(index) = modine%chi_atoms(idim, iatom)
        end do
      end do

      call droot_solver_run(modine%rs, getf2, my_csi, conv, startval=start_csi)

      if (.not. conv) then
        message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
        message(2) = "method did not converge."
        call messages_fatal(2, namespace=namespace)
      end if

      ! Now set csi to the new values
      do iatom = 1, modine%npos
        do idim = 1, modine%dim
          index = (iatom-1)*modine%dim + idim
          modine_p%csi(idim, iatom) = my_csi(index)
        end do
      end do

      SAFE_DEALLOCATE_A(x_p)
      SAFE_DEALLOCATE_A(my_csi)
      SAFE_DEALLOCATE_A(start_csi)

      nullify(modine_p)

      POP_SUB(curv_modine_constructor.optimize)
    end subroutine optimize

  end function curv_modine_constructor

  ! ---------------------------------------------------------
  subroutine curv_modine_copy(this_out, this_in)
    type(curv_modine_t), intent(inout) :: this_out
    type(curv_modine_t), intent(in)    :: this_in

    PUSH_SUB(curv_modine_copy)

    this_out%dim = this_in%dim
    this_out%local_basis = this_in%local_basis
    SAFE_ALLOCATE_SOURCE_A(this_out%lsize, this_in%lsize)
    this_out%xbar = this_in%xbar
    this_out%Jbar = this_in%Jbar
    SAFE_ALLOCATE_SOURCE_A(this_out%Jlocal, this_in%Jlocal)
    SAFE_ALLOCATE_SOURCE_A(this_out%Jrange, this_in%Jrange)
    SAFE_ALLOCATE_SOURCE_A(this_out%chi_atoms, this_in%chi_atoms)
    SAFE_ALLOCATE_SOURCE_A(this_out%csi, this_in%csi)
    this_out%npos = this_in%npos

    POP_SUB(curv_modine_copy)
  end subroutine curv_modine_copy

  ! ---------------------------------------------------------
  subroutine curv_modine_finalize(this)
    type(curv_modine_t), intent(inout) :: this

    PUSH_SUB(curv_modine_finalize)

    SAFE_DEALLOCATE_A(this%lsize)
    SAFE_DEALLOCATE_A(this%Jlocal)
    SAFE_DEALLOCATE_A(this%Jrange)
    SAFE_DEALLOCATE_A(this%chi_atoms)
    SAFE_DEALLOCATE_A(this%csi)

    POP_SUB(curv_modine_finalize)
  end subroutine curv_modine_finalize

  ! ---------------------------------------------------------
  function curv_modine_to_cartesian(this, chi) result(xx)
    class(curv_modine_t), target, intent(in)  :: this
    FLOAT,                        intent(in)  :: chi(:)
    FLOAT :: xx(1:this%dim)

    FLOAT :: chi2(this%dim), rr2, dd
    integer :: iatom

    ! no PUSH_SUB, called too often

    call curv_modine_chi2chi2(this, chi, chi2)

    xx(:) = chi2(:)
    do iatom = 1, this%npos
      rr2 = sum((chi2(:) - this%csi(:, iatom))**2)
      dd = exp(-rr2/(M_TWO*this%Jrange(iatom)**2))

      xx(:) = xx(:) - this%Jlocal(iatom)*(chi2(:) - this%csi(:, iatom)) * dd
    end do

  end function curv_modine_to_cartesian

  ! ---------------------------------------------------------
  function curv_modine_from_cartesian(this, xx) result(chi)
    class(curv_modine_t), target, intent(in)  :: this
    FLOAT,                        intent(in)  :: xx(:)
    FLOAT :: chi(1:this%dim)

    logical :: conv

    ! no PUSH_SUB, called too often

    modine_p  => this
    SAFE_ALLOCATE(x_p(1:this%dim))
    x_p(:) = xx(:)

    call droot_solver_run(this%rs, getf, chi, conv, startval = xx)

    SAFE_DEALLOCATE_A(x_p)
    nullify(modine_p)

    if (.not. conv) then
      message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
      message(2) = "method did not converge for point:"
      write(message(3),'(3f14.6)') xx(1:this%dim)
      call messages_fatal(3)
    end if

  end function curv_modine_from_cartesian

  ! ---------------------------------------------------------
  FLOAT function curv_modine_det_jac(this, xx, chi) result(jdet)
    class(curv_modine_t), intent(in)  :: this
    FLOAT,                intent(in)  :: xx(:)
    FLOAT,                intent(in)  :: chi(:)

    FLOAT :: dummy(this%dim), jac(1:this%dim, 1:this%dim)

    ! no PUSH_SUB, called too often

    call curv_modine_jacobian_inv(this, chi, dummy, Jac)
    jdet = lalg_determinant(this%dim, jac, preserve_mat = .false.)

  end function curv_modine_det_jac

  ! ---------------------------------------------------------
  subroutine curv_modine_write_info(this, iunit, namespace)
    class(curv_modine_t),         intent(in) :: this
    integer,            optional, intent(in) :: iunit
    type(namespace_t),  optional, intent(in) :: namespace

    PUSH_SUB(curv_modine_write_info)

    write(message(1), '(a)') '  Curvilinear Method = modine'
    call messages_info(1, iunit=iunit, namespace=namespace)

    POP_SUB(curv_modine_write_info)
  end subroutine curv_modine_write_info

  ! ---------------------------------------------------------
  FLOAT function curv_modine_surface_element(this, idir) result(ds)
    class(curv_modine_t), intent(in) :: this
    integer,              intent(in) :: idir

    ds = M_ZERO
    message(1) = 'Surface element with modine curvilinear coordinates not implemented'
    call messages_fatal(1)

  end function curv_modine_surface_element

  ! ---------------------------------------------------------
  pure subroutine curv_modine_chi2chi2(this, chi_, chi2, Jac)
    class(curv_modine_t), intent(in)  :: this
    FLOAT,                intent(in)  :: chi_(:)
    FLOAT,                intent(out) :: chi2(:)
    FLOAT,      optional, intent(out) :: Jac(:)   !< the Jacobian of this transformation is diagonal

    integer, parameter :: qq = 3
    FLOAT :: chibar(this%dim), rr, chi
    logical :: neg
    integer :: i

    ! no PUSH_SUB, called too often

    chibar = this%xbar*this%lsize

    do i = 1, this%dim
      neg = (chi_(i) < 0)
      chi = abs(chi_(i))

      chi2(i)  = this%Jbar * chi
      if (present(Jac)) Jac(i) = this%Jbar

      if (chi > chibar(i)) then
        rr = (chi - chibar(i))/(this%lsize(i) - chibar(i))

        chi2(i)  = chi2(i) + this%lsize(i)/M_TWO*(1 - this%Jbar) * rr**qq * &
          (qq + M_ONE - (qq - M_ONE)*rr)

        if (present(Jac)) then
          Jac(i) = Jac(i) + this%lsize(i)/M_TWO*(M_ONE - this%Jbar) * rr**(qq - 1)/(this%lsize(i) - chibar(i)) * &
            (qq*(qq + M_ONE) - (qq**2 - M_ONE)*rr)
        end if
      end if

      if (neg) chi2(i) = -chi2(i)
      ! CHECK if Jacobian does not have to be negated!
    end do

  end subroutine curv_modine_chi2chi2

  ! ---------------------------------------------------------
  pure subroutine curv_modine_jacobian_inv(this, chi, xx, Jac)
    type(curv_modine_t), intent(in)  :: this
    FLOAT,               intent(in)  :: chi(:)
    FLOAT,               intent(out) :: xx(:)
    FLOAT,               intent(out) :: Jac(:, :) !< the Jacobian

    FLOAT :: chi2(this%dim), rr2, dd, J2(this%dim)
    integer :: iatom, idim, idim2

    ! no PUSH_SUB, called too often

    call curv_modine_chi2chi2(this, chi, chi2, J2)

    ! initialize both xx and the Jacobian
    xx(:) = chi2(:)
    Jac(:,:) = M_ZERO
    do idim = 1, this%dim
      Jac(idim, idim) = M_ONE
    end do

    do iatom = 1, this%npos
      rr2 = sum((chi2(:) - this%csi(:, iatom))**2)
      dd = exp(-rr2/(M_TWO*this%Jrange(iatom)**2))

      xx(:) = xx(:) -  this%Jlocal(iatom)*(chi2(:) - this%csi(:, iatom)) * dd

      do idim = 1, this%dim
        Jac(idim, idim) = Jac(idim, idim) - this%Jlocal(iatom) * dd
        do idim2 = 1, this%dim
          Jac(idim, idim2) = Jac(idim, idim2) + &
            this%Jlocal(iatom)*(chi2(idim) - this%csi(idim, iatom))*(chi2(idim2) - this%csi(idim2, iatom)) * &
            M_TWO/(M_TWO*this%Jrange(iatom)**2) * dd
        end do
      end do
    end do

    do idim = 1, this%dim
      Jac(idim, :) = Jac(idim, :) * J2(:)
    end do

  end subroutine curv_modine_jacobian_inv

  ! ---------------------------------------------------------
  pure subroutine getf(yy, ff, jf)
    FLOAT, intent(in)  :: yy(:)
    FLOAT, intent(out) :: ff(:), jf(:, :)

    ! no PUSH_SUB, called too often

    call curv_modine_jacobian_inv(modine_p, yy, ff, jf)
    ff(:) = ff(:) - x_p(:)

  end subroutine getf

  ! ---------------------------------------------------------
  subroutine getf2(csi, ff, jf)
    FLOAT, intent(in)  :: csi(:)
    FLOAT, intent(out) :: ff(:), jf(:, :)

    integer :: i1, j1, i2, j2, index1, index2
    FLOAT :: rr2, dd, dd2
    FLOAT, allocatable :: xx(:), chi2(:)

    ! no PUSH_SUB, called too often

    SAFE_ALLOCATE(xx(1:modine_p%dim))
    SAFE_ALLOCATE(chi2(1:modine_p%dim))

    ! first we fill in coord_system%csi with the values we have
    index1 = 1
    do i1 = 1, modine_p%npos
      do j1 = 1, modine_p%dim
        modine_p%csi(j1, i1) = csi(index1)
        index1 = index1 + 1
      end do
    end do

    ! get ff and jf
    jf(:,:) = M_ZERO
    do i1 = 1, modine_p%npos
      call curv_modine_chi2chi2(modine_p, modine_p%chi_atoms(:,i1), chi2)
      xx = chi2

      do i2 = 1, modine_p%npos
        rr2 = sum((chi2 - modine_p%csi(:,i2))**2)
        dd = exp(-rr2/(M_TWO*modine_p%Jrange(i2)**2))

        xx = xx - modine_p%Jlocal(i2)*(chi2 - modine_p%csi(:,i2)) * dd
      end do

      do j1 = 1, modine_p%dim
        index1 = (i1 - 1)*modine_p%dim + j1
        ff(index1) = xx(j1) - x_p(index1)

        do i2 = 1, modine_p%npos
          rr2  = sum((chi2 - modine_p%csi(:,i2))**2)
          dd  = exp(-rr2/(M_TWO*modine_p%Jrange(i2)**2))
          dd2 = -M_TWO/(M_TWO*modine_p%Jrange(i2)**2)*dd

          index2 = (i2 - 1)*modine_p%dim + j1
          jf(index1, index2) = modine_p%Jlocal(i2) * dd

          do j2 = 1, modine_p%dim
            index2 = (i2 - 1)*modine_p%dim + j2

            jf(index1, index2) =  jf(index1, index2) + modine_p%Jlocal(i2) * dd2 * &
              (chi2(j1) - modine_p%csi(j1,i2))*(chi2(j2) - modine_p%csi(j2,i2))
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(xx)
    SAFE_DEALLOCATE_A(chi2)

  end subroutine getf2

end module curv_modine_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
