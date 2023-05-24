!! Copyright (C) 2021 N. Tancogne-Dejean
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

module magnetic_constrain_oct_m
  use debug_oct_m
  use global_oct_m
  use ions_oct_m
  use magnetic_oct_m
  use messages_oct_m
  use mesh_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m

  implicit none

  private
  public ::                          &
    magnetic_constrain_t,            &
    magnetic_constrain_init,         &
    magnetic_constrain_update,       &
    magnetic_constrain_end

  integer, public, parameter ::        &
    CONSTRAIN_NONE                = 0, &
    CONSTRAIN_DIR                 = 1, &
    CONSTRAIN_FULL                = 2


  type magnetic_constrain_t
    private
    integer, public    :: level
    FLOAT              :: lambda !< Lagrange multiplier
    FLOAT, allocatable :: constrain(:,:) !< Constrained moments (3, ions%natoms)
    FLOAT, allocatable, public :: pot(:,:) !< Potential (mesh%np, d%nspin)
    FLOAT,              public :: energy   !< Energy

    FLOAT              :: lmm_r

    FLOAT, pointer     :: rho(:,:)
  end type magnetic_constrain_t

contains

  subroutine magnetic_constrain_init(this, namespace, mesh, std, rho, ions)
    type(magnetic_constrain_t), intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace
    type(mesh_t),              intent(in)    :: mesh
    type(states_elec_dim_t),   intent(in)    :: std
    FLOAT, target,             intent(in)    :: rho(:,:)
    type(ions_t),              intent(in)    :: ions

    integer :: ia, idir
    FLOAT   :: rmin, norm
    type(block_t) :: blk

    PUSH_SUB(magnetic_constrain_init)

    !%Variable MagneticConstrain
    !%Type integer
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% This variable selects which magnetic constrain expression is added to the Hamiltonian.
    !%Option constrain_none 0
    !% No constrain is added to the Hamiltonian.
    !%Option constrain_dir 1
    !% We are adding a constrain for the direction of the magnetic moments only,
    !% see PRB 91, 054420 (2015).
    !%Option constrain_full 2
    !% We are adding a constrain for the direction and norm of the magnetic moments only,
    !% see PRB 91, 054420 (2015).
    !%End
    call parse_variable(namespace, 'MagneticConstrain', CONSTRAIN_NONE, this%level)
    call messages_print_var_value('MagneticConstrain', this%level, namespace=namespace)
    if (this%level == CONSTRAIN_NONE) then
      POP_SUB(magnetic_constrain_init)
      return
    end if

    call messages_experimental('MagneticConstrain', namespace=namespace)

    !%Variable MagneticConstrainStrength
    !%Type float
    !%Default 0.01
    !%Section Hamiltonian
    !%Description
    !% This variable determines the value of the Lagrange multiplier used for the constrain term.
    !%End
    call parse_variable(namespace, 'MagneticConstrainStrength', CNST(0.01), this%lambda)
    call messages_print_var_value('MagneticConstrainStrength', this%lambda, namespace=namespace)


    !We are using the AtomsMagnetDirection block for the contrain
    if (parse_block(namespace, 'AtomsMagnetDirection', blk) < 0) then
      message(1) = "AtomsMagnetDirection block is not defined."
      message(2) = "Magnetic constrained DFT cannot be used without it."
      call messages_fatal(2, namespace=namespace)
    end if

    if (parse_block_n(blk) /= ions%natoms) then
      message(1) = "AtomsMagnetDirection block has the wrong number of rows."
      call messages_fatal(1, namespace=namespace)
    end if

    if(std%ispin == UNPOLARIZED) then
      message(1) = "Magnetic constrains can only be used for spin-polized and spinor calculations."
      call messages_fatal(1, namespace=namespace)
    end if

    SAFE_ALLOCATE(this%constrain(1:3, ions%natoms))
    do ia = 1, ions%natoms
      !Read from AtomsMagnetDirection block
      if (std%ispin == SPIN_POLARIZED) then
        call parse_block_float(blk, ia-1, 0, this%constrain(3, ia))
        this%constrain(1:2, ia) = M_ZERO
      elseif (std%ispin == SPINORS) then
        do idir = 1, 3
          call parse_block_float(blk, ia-1, idir-1, this%constrain(idir, ia))
          if (abs(this%constrain(idir, ia)) < CNST(1.0e-20)) this%constrain(idir, ia) = M_ZERO
        end do
      end if

      if(this%level == CONSTRAIN_DIR) then
        !Renormalization for constrained directions
        norm = norm2(this%constrain(1:3, ia))
        this%constrain(1:3, ia) = this%constrain(1:3, ia) / norm
      end if
    end do
    call parse_block_end(blk)

    SAFE_ALLOCATE(this%pot(1:mesh%np, std%nspin))
    this%pot = M_ZERO

    rmin = ions%min_distance()
    call parse_variable(namespace, 'LocalMagneticMomentsSphereRadius', min(M_HALF*rmin, CNST(100.0)), &
      this%lmm_r)

    this%rho => rho

    POP_SUB(magnetic_constrain_init)
  end subroutine magnetic_constrain_init

  subroutine magnetic_constrain_end(this)
    type(magnetic_constrain_t), intent(inout) :: this

    PUSH_SUB(magnetic_constrain_end)

    SAFE_DEALLOCATE_A(this%constrain)
    SAFE_DEALLOCATE_A(this%pot)
    nullify(this%rho)

    POP_SUB(magnetic_constrain_end)
  end subroutine magnetic_constrain_end

! ---------------------------------------------------------
  subroutine magnetic_constrain_update(this, mesh, std, ions)
    type(magnetic_constrain_t), intent(inout) :: this
    type(mesh_t),               intent(in)    :: mesh
    type(states_elec_dim_t),    intent(in)    :: std
    type(ions_t),               intent(in)    :: ions

    integer :: ia, idir, ip
    FLOAT :: bb(3), b2, lmm(3), dotp, norm, xx
    FLOAT, allocatable :: md(:,:), mdf(:), mask(:)
    type(submesh_t) :: sphere
    type(profile_t), save :: prof

    if (this%level == CONSTRAIN_NONE) return

    PUSH_SUB(magnetic_constrain_update)

    call profiling_in(prof, TOSTRING(MAGNETIC_CONSTRAIN))

    this%pot = M_ZERO
    this%energy = M_ZERO

    SAFE_ALLOCATE(md(1:mesh%np, 1:max(mesh%box%dim, 3)))
    SAFE_ALLOCATE(mdf(1:mesh%np))
    call magnetic_density(mesh, std, this%rho, md)

    do ia = 1, ions%natoms
      call submesh_init(sphere, ions%space, mesh, ions%latt, ions%pos(:, ia), this%lmm_r)

      SAFE_ALLOCATE(mask(1:sphere%np))
      mask = M_ZERO
      do ip = 1, sphere%np
        xx = sphere%r(ip) * M_PI / this%lmm_r
        !if(xx < 0.1) then
        !  mask(ip) = M_ONE-xx*xx/CNST(6.0)
        !else
        !  mask(ip) = sin(xx)/xx
        !end if
        mask(ip) = M_ONE
      end do

      ! We compute the local moment here
      ! We multiply here by a function that decreases monotonically to zero
      !  at the boundary of the atomic sphere.
      mdf = M_ZERO

      do idir = 1, max(mesh%box%dim, 3)
        do ip = 1, sphere%np
          mdf(sphere%map(ip)) = md(sphere%map(ip), idir) * mask(ip)
        end do
        lmm(idir) = dsm_integrate_frommesh(mesh, sphere, mdf)
      end do

      ! See PRB 91, 054420 (2015)
      ! See for instance https://www.vasp.at/wiki/index.php/I_CONSTRAINED_M
      select case(this%level)
      case(CONSTRAIN_DIR)
        norm = max(norm2(lmm),CNST(1e-10))
        dotp = dot_product(lmm, this%constrain(1:3, ia))
        this%energy = this%energy + this%lambda * (norm - dotp)
        bb = lmm/norm - this%constrain(:,ia)
      case(CONSTRAIN_FULL)
        this%energy = this%energy + this%lambda * sum((lmm - this%constrain(:,ia))**2)
        bb = M_TWO * (lmm - this%constrain(:,ia))
      end select

      ! Here we need to be carefull because we might end up in the case where
      ! the constrained is fullfilled in one direction, say along z
      ! and small noise in the x-y plane would then break the symmetry and
      ! if we have independent particles, this determines the magnetization direction
      !
      ! To prevent problems with the indenpendent particles,
      ! we just set a small fraction of the constrain instead of zero, to break
      ! the symmetries
      do idir = 1, max(mesh%box%dim, 3)
        if(abs(bb(idir)) < CNST(1e-3)) then
          bb(idir) = -this%constrain(idir, ia)*CNST(1e-2)
        else
          bb(idir) = bb(idir) * this%lambda
        end if
      end do

      ! We are adding an effective Zeeman term within the atomic sphere
      select case(std%ispin)
      case (SPIN_POLARIZED)
        b2 = norm2(bb(1:max(mesh%box%dim, 3)))
        do ip = 1, sphere%np
          this%pot(sphere%map(ip), 1) = this%pot(sphere%map(ip), 1) + b2 * mask(ip)
          this%pot(sphere%map(ip), 2) = this%pot(sphere%map(ip), 2) - b2 * mask(ip)
        end do

      case (SPINORS)
        do ip = 1, sphere%np
          this%pot(sphere%map(ip), 1) = this%pot(sphere%map(ip), 1) + bb(3) * mask(ip)
          this%pot(sphere%map(ip), 2) = this%pot(sphere%map(ip), 2) - bb(3) * mask(ip)
          this%pot(sphere%map(ip), 3) = this%pot(sphere%map(ip), 3) + bb(1) * mask(ip)
          this%pot(sphere%map(ip), 4) = this%pot(sphere%map(ip), 4) - bb(2) * mask(ip)
        end do

      end select

      SAFE_DEALLOCATE_A(mask)
      call submesh_end(sphere)
    end do

    SAFE_DEALLOCATE_A(md)
    SAFE_DEALLOCATE_A(mdf)

    call profiling_out(prof)

    POP_SUB(magnetic_constrain_update)
  end subroutine magnetic_constrain_update

end module magnetic_constrain_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
