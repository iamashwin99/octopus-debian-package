!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Verstraete
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


! ---------------------------------------------------------
!> integrates a function
R_TYPE function X(mf_integrate) (mesh, ff, mask, reduce) result(dd)
  class(mesh_t), intent(in) :: mesh
  R_TYPE,       intent(in) :: ff(:)  !< (mesh%np)
  logical, optional, intent(in) :: mask(:)
  logical, optional, intent(in) :: reduce

  integer :: ip

  call profiling_in(X(PROFILING_MF_INTEGRATE), TOSTRING(X(MF_INTEGRATE)))
  PUSH_SUB(X(mf_integrate))

  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)
  ASSERT(not_in_openmp())

  !TODO: This case need to be implemented
  if (present(mask)) then
    ASSERT(.not. mesh%use_curvilinear)
  end if

  dd = R_TOTYPE(M_ZERO)
  if (mesh%use_curvilinear) then
    !$omp parallel do reduction(+:dd)
    do ip = 1, mesh%np
      dd = dd + ff(ip)*mesh%vol_pp(ip)
    end do
  else if (present(mask)) then
    dd = sum(ff(1:mesh%np), mask=mask(1:mesh%np))
  else
    !$omp parallel do reduction(+:dd)
    do ip = 1, mesh%np
      dd = dd + ff(ip)
    end do
  end if

  dd = dd*mesh%volume_element

  if (mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(X(PROFILING_MF_REDUCE), TOSTRING(X(MF_REDUCE)))
    call mesh%allreduce(dd)
    call profiling_out(X(PROFILING_MF_REDUCE))
  end if

  POP_SUB(X(mf_integrate))
  call profiling_out(X(PROFILING_MF_INTEGRATE))

end function X(mf_integrate)

! ---------------------------------------------------------
subroutine X(mf_normalize)(mesh, dim, psi, norm)
  class(mesh_t),   intent(in)    :: mesh
  integer,         intent(in)    :: dim
  R_TYPE,          intent(inout) :: psi(:,:)
  FLOAT, optional, intent(out)   :: norm

  FLOAT   :: norm_

  PUSH_SUB(X(mf_normalize))

  norm_ = X(mf_nrm2) (mesh, dim, psi)
  if (abs(norm_) <= CNST(1e-150)) then
    message(1) = "Mesh function has zero norm; cannot normalize."
    call messages_fatal(1)
  end if

  call lalg_scal(mesh%np, dim, M_ONE / norm_, psi)

  if (present(norm)) then
    norm = norm_
  end if

  POP_SUB(X(mf_normalize))
end subroutine X(mf_normalize)



!> ---------------------------------------------------------
!! This function returns the dot product between two vectors,
!! but using the mesh_aux defined as a global object in this
!! module. This way it can be called by external libraries,
!! passing only the two vectors. First, one has to
!! make sure that mesh_aux is pointing to some defined
!! mesh data structure, by calling mesh_init_mesh_aux.
!! ---------------------------------------------------------
R_TYPE function X(mf_dotp_aux)(f1, f2) result(dotp)
  R_TYPE,            intent(in) :: f1(:), f2(:)

  PUSH_SUB(X(mf_dotp_aux))

  ASSERT(associated(mesh_aux))
  dotp = X(mf_dotp)(mesh_aux, f1, f2)

  POP_SUB(X(mf_dotp_aux))
end function X(mf_dotp_aux)

!> Same as above, but no conjugation.
!! ---------------------------------------------------------
R_TYPE function X(mf_dotu_aux)(f1, f2) result(dotu)
  R_TYPE,            intent(in) :: f1(:), f2(:)

  PUSH_SUB(X(mf_dotu_aux))

  ASSERT(associated(mesh_aux))
  dotu = X(mf_dotp)(mesh_aux, f1, f2, dotu = .true.)

  POP_SUB(X(mf_dotu_aux))
end function X(mf_dotu_aux)

!> Same as above, but for norm.
!! ---------------------------------------------------------
FLOAT function X(mf_nrm2_aux)(ff) result(norm)
  R_TYPE,            intent(in) :: ff(:)

  PUSH_SUB(X(mf_nrm2_aux))

  ASSERT(associated(mesh_aux))
  norm = X(mf_nrm2)(mesh_aux, ff)

  POP_SUB(X(mf_nrm2_aux))
end function X(mf_nrm2_aux)


! ---------------------------------------------------------
!> this function returns the dot product between two vectors
R_TYPE function X(mf_dotp_1)(mesh, f1, f2, reduce, dotu, np) result(dotp)
  class(mesh_t),     intent(in) :: mesh
  R_TYPE,            intent(in) :: f1(:), f2(:)
  logical, optional, intent(in) :: reduce
  logical, optional, intent(in) :: dotu
  !< if true, use blas_dotu instead of blas_dot;
  !! no complex conjugation.  Default is false.
  !! has no effect if working with real version
  integer, optional, intent(in) :: np

#ifdef R_TCOMPLEX
  logical             :: dotu_
#endif
  integer             :: ip, np_

  PUSH_SUB(X(mf_dotp_1))
  call profiling_in(X(PROFILING_MF_DOTP), TOSTRING(X(MF_DOTP)))

  ASSERT(not_in_openmp())

  np_ = optional_default(np, mesh%np)

  if (np_ == 0) then
    POP_SUB(X(mf_dotp_1))
    return
  end if

  ASSERT(ubound(f1, dim = 1) == np_ .or. ubound(f1, dim = 1) == mesh%np_part)
  ASSERT(ubound(f2, dim = 1) == np_ .or. ubound(f2, dim = 1) == mesh%np_part)

#ifdef R_TCOMPLEX
  dotu_ = optional_default(dotu, .false.)
#endif

  if (mesh%use_curvilinear) then
    dotp = R_TOTYPE(M_ZERO)
    ! preprocessor conditionals necessary since blas_dotu only exists for complex input
#ifdef R_TCOMPLEX
    if (.not. dotu_) then
#endif
      !$omp parallel do reduction(+:dotp)
      do ip = 1, np_
        dotp = dotp + mesh%vol_pp(ip)*R_CONJ(f1(ip))*f2(ip)
      end do
#ifdef R_TCOMPLEX
    else
      !$omp parallel do reduction(+:dotp)
      do ip = 1, np_
        dotp = dotp + mesh%vol_pp(ip)*f1(ip)*f2(ip)
      end do
    end if
#endif
    call profiling_count_operations(np_*(2*R_ADD + R_MUL))
  else
#ifdef R_TCOMPLEX
    if (.not. dotu_) then
#endif
      dotp = blas_dot(np_, f1(1), 1, f2(1), 1)
#ifdef R_TCOMPLEX
    else
      dotp = blas_dotu(np_, f1(1), 1, f2(1), 1)
    end if
#endif
    call profiling_count_operations(np_*(R_ADD + R_MUL))

  end if

  dotp = dotp*mesh%volume_element

  if (mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(X(PROFILING_MF_REDUCE), TOSTRING(X(MF_REDUCE)))
    call mesh%allreduce(dotp)
    call profiling_out(X(PROFILING_MF_REDUCE))
  end if

  call profiling_out(X(PROFILING_MF_DOTP))
  POP_SUB(X(mf_dotp_1))
end function X(mf_dotp_1)


! ---------------------------------------------------------
R_TYPE function X(mf_dotp_2)(mesh, dim, f1, f2, reduce, dotu, np) result(dotp)
  class(mesh_t),     intent(in) :: mesh
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f1(:,:), f2(:,:)
  logical, optional, intent(in) :: reduce
  logical, optional, intent(in) :: dotu
  !< if true, use lalg_dotu instead of lalg_dot;
  !! no complex conjugation.  Default is false.
  integer, optional, intent(in) :: np

  integer :: idim

  PUSH_SUB(X(mf_dotp_2))

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp_1)(mesh, f1(:, idim), f2(:, idim), reduce = .false., dotu = dotu, np = np)
  end do

  if (mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(X(PROFILING_MF_REDUCE), TOSTRING(X(MF_REDUCE)))
    call mesh%allreduce(dotp)
    call profiling_out(X(PROFILING_MF_REDUCE))
  end if

  POP_SUB(X(mf_dotp_2))

end function X(mf_dotp_2)


! ---------------------------------------------------------
!> this function returns the the norm of a vector
FLOAT function X(mf_nrm2_1)(mesh, ff, reduce) result(nrm2)
  class(mesh_t),     intent(in) :: mesh
  R_TYPE,            intent(in) :: ff(:)
  logical, optional, intent(in) :: reduce

  R_TYPE, allocatable :: ll(:)

  call profiling_in(X(PROFILING_MF_NRM2), TOSTRING(X(MF_NRM2)))
  PUSH_SUB(X(mf_nrm2_1))

  if (mesh%use_curvilinear) then
    SAFE_ALLOCATE(ll(1:mesh%np))
    ll(1:mesh%np) = ff(1:mesh%np)*sqrt(mesh%vol_pp(1:mesh%np))
    nrm2 = lalg_nrm2(mesh%np, ll)
    SAFE_DEALLOCATE_A(ll)
  else
    nrm2 = lalg_nrm2(mesh%np, ff)
  end if

  nrm2 = nrm2*sqrt(mesh%volume_element)

  if (mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(X(PROFILING_MF_REDUCE), TOSTRING(X(MF_REDUCE)))
    nrm2 = nrm2**2
    call mesh%allreduce(nrm2)
    nrm2 = sqrt(nrm2)
    call profiling_out(X(PROFILING_MF_REDUCE))
  end if

  POP_SUB(X(mf_nrm2_1))
  call profiling_out(X(PROFILING_MF_NRM2))

end function X(mf_nrm2_1)

! ---------------------------------------------------------
FLOAT function X(mf_nrm2_2)(mesh, dim, ff, reduce) result(nrm2)
  class(mesh_t),     intent(in) :: mesh
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: ff(:,:)
  logical, optional, intent(in) :: reduce

  integer :: idim

  PUSH_SUB(X(mf_nrm2_2))

  nrm2 = M_ZERO

  do idim = 1, dim
    nrm2 = hypot(nrm2, X(mf_nrm2)(mesh, ff(:, idim), reduce = reduce))
  end do

  POP_SUB(X(mf_nrm2_2))

end function X(mf_nrm2_2)


! ---------------------------------------------------------
!> This function calculates the "order" moment of the function ff
R_TYPE function X(mf_moment) (mesh, ff, idir, order) result(rr)
  class(mesh_t),intent(in) :: mesh
  R_TYPE,       intent(in) :: ff(:)
  integer,      intent(in) :: idir
  integer,      intent(in) :: order

  R_TYPE, allocatable :: fxn(:)
  integer :: ip

  PUSH_SUB(X(mf_moment))

  ASSERT(not_in_openmp())

  SAFE_ALLOCATE(fxn(1:mesh%np))

  !$omp parallel do
  do ip = 1, mesh%np
    fxn(ip) = ff(ip)*mesh%x(ip, idir)**order
  end do
  rr = X(mf_integrate)(mesh, fxn)

  SAFE_DEALLOCATE_A(fxn)

  POP_SUB(X(mf_moment))
end function X(mf_moment)

! ---------------------------------------------------------
!> This subroutine fills a function with randon values.
subroutine X(mf_random)(mesh, ff, pre_shift, post_shift, seed, normalized)
  class(mesh_t),         intent(in)  :: mesh
  R_TYPE,                intent(out) :: ff(:)
  integer(i8), optional, intent(in)  :: pre_shift
  integer(i8), optional, intent(in)  :: post_shift
  integer, optional,     intent(in)  :: seed
  logical, optional,     intent(in)  :: normalized !< whether generate states should have norm 1, true by default

  integer, save :: iseed = 123
  R_BASE  :: rr
  type(profile_t), save :: prof

  PUSH_SUB(X(mf_random))

  call profiling_in(prof, TOSTRING(X(RANDOMIZE)))

  if (present(seed)) then
    iseed = iseed + seed
  end if

  if (present(pre_shift)) then
    !We skip shift times the seed to compensate for MPI tasks dealing with previous mesh points
    call shiftseed(iseed, pre_shift)
#if defined(R_TCOMPLEX)
    ! for complex wave functions we need to shift twice (real and imag part).
    call shiftseed(iseed, pre_shift)
#endif
  end if

  call quickrnd(iseed, mesh%np, ff(1:mesh%np))

  if (present(post_shift)) then
    !We skip shift times the seed to compensate for MPI tasks dealing with posteriour mesh points
    call shiftseed(iseed, post_shift)
#if defined(R_TCOMPLEX)
    ! for complex wave functions we need to shift twice (real and imag part).
    call shiftseed(iseed, post_shift)
#endif
  end if


  if (optional_default(normalized, .true.)) then
    rr = X(mf_nrm2)(mesh, ff)
    call lalg_scal(mesh%np, R_TOTYPE(1.0)/rr, ff)
  end if

  call profiling_out(prof)

  POP_SUB(X(mf_random))
end subroutine X(mf_random)

! ---------------------------------------------------------
!> This function receives a function f_in defined in a mesh, and returns
!! the interpolated values of the function over the npoints_in defined
!! by x_in.
subroutine X(mf_interpolate_points) (ndim, npoints_in, x_in, f_in, npoints_out, x_out, f_out)
  integer, intent(in)  :: ndim, npoints_in, npoints_out
  R_TYPE,  intent(in)  :: f_in(:)    !< (npoints_in)
  FLOAT,   intent(in)  :: x_in(:, :)
  FLOAT,   intent(in)  :: x_out(:,:)
  R_TYPE,  intent(out) :: f_out(:)   !< (npoints_out)

  FLOAT :: pp(ndim)
  integer :: ip
  type(qshep_t) :: interp
#ifndef R_TCOMPLEX
  type(spline_t) :: interp1d
#endif

  PUSH_SUB(X(mf_interpolate_points))

  select case (ndim)
  case (2)
    call qshep_init(interp, i4_to_i8(npoints_in), f_in, x_in(:, 1), x_in(:, 2))
    do ip = 1, npoints_out
      pp(1:2)   = x_out(ip, 1:2)
      f_out(ip) = qshep_interpolate(interp, f_in, pp(1:2))
    end do
    call qshep_end(interp)

  case (3)
    call qshep_init(interp, i4_to_i8(npoints_in), f_in, x_in(:, 1), x_in(:, 2), x_in(:, 3))
    do ip = 1, npoints_out
      pp(1:3)   = x_out(ip, 1:3)
      f_out(ip) = qshep_interpolate(interp, f_in, pp(1:3))
    end do
    call qshep_end(interp)

  case (1)
#ifdef R_TCOMPLEX
    message(1) = 'Believe it or not: cannot do 1D complex interpolation, only 2D or 3D.'
    call messages_fatal(1)
#else
    call spline_init(interp1d)
    call spline_fit(npoints_in, x_in(:, 1), f_in, interp1d)
    do ip = 1, npoints_out
      f_out(ip) = spline_eval(interp1d, x_out(ip, 1))
    end do
    call spline_end(interp1d)
#endif
  end select

  POP_SUB(X(mf_interpolate_points))
end subroutine X(mf_interpolate_points)

! ---------------------------------------------------------
!> Given a function ff defined on mesh, and a plane, it gives
!! back the values of ff on the plane, by doing the appropriate
!! interpolation.
subroutine X(mf_interpolate_on_plane)(mesh, plane, ff, f_in_plane)
  class(mesh_t),      intent(in)  :: mesh
  type(mesh_plane_t), intent(in)  :: plane
  R_TYPE, target,     intent(in)  :: ff(:)
  R_TYPE,             intent(out) :: f_in_plane(plane%nu:plane%mu, plane%nv:plane%mv)

  integer :: iu, iv
  integer(i8) :: ip
  R_TYPE, pointer :: f_global(:)
  FLOAT :: pp(3)
  type(qshep_t) :: interp
  FLOAT, allocatable :: xglobal(:, :)

  PUSH_SUB(X(mf_interpolate_on_plane))

  ASSERT(not_in_openmp())

  SAFE_ALLOCATE(xglobal(1:mesh%np_part_global, 1:mesh%box%dim))
  !$omp parallel do
  do ip = 1, mesh%np_part_global
    xglobal(ip, 1:mesh%box%dim) = mesh_x_global(mesh, ip)
  end do

  if (mesh%parallel_in_domains) then
    SAFE_ALLOCATE(f_global(1:mesh%np_global))
    call par_vec_allgather(mesh%pv, f_global, ff)
  else
    f_global => ff
  end if

  call qshep_init(interp, mesh%np_global, f_global, xglobal(:, 1), xglobal(:, 2), xglobal(:, 3))

  do iu = plane%nu, plane%mu
    do iv = plane%nv, plane%mv
      pp(1) = plane%origin(1) + iu*plane%spacing * plane%u(1) + iv * plane%spacing * plane%v(1)
      pp(2) = plane%origin(2) + iu*plane%spacing * plane%u(2) + iv * plane%spacing * plane%v(2)
      pp(3) = plane%origin(3) + iu*plane%spacing * plane%u(3) + iv * plane%spacing * plane%v(3)
      f_in_plane(iu, iv) = qshep_interpolate(interp, f_global, pp(1:3))
    end do
  end do

  call qshep_end(interp)

  SAFE_DEALLOCATE_A(xglobal)
  if (mesh%parallel_in_domains) then
    SAFE_DEALLOCATE_P(f_global)
  end if

  POP_SUB(X(mf_interpolate_on_plane))
end subroutine X(mf_interpolate_on_plane)

! ---------------------------------------------------------
!> Given a function ff defined on mesh, and a line, it gives
!! back the values of ff on the line, by doing the appropriate
!! interpolation.
subroutine X(mf_interpolate_on_line)(mesh, line, ff, f_in_line)
  class(mesh_t),      intent(in)  :: mesh
  type(mesh_line_t),  intent(in)  :: line
  R_TYPE, target,     intent(in)  :: ff(:)
  R_TYPE,             intent(out) :: f_in_line(line%nu:line%mu)

  integer :: iu
  integer(i8) :: ip
  R_TYPE, pointer :: f_global(:)
  FLOAT :: pp(2)
  type(qshep_t) :: interp
  FLOAT , allocatable :: xglobal(:, :)

  PUSH_SUB(X(mf_interpolate_on_line))

  SAFE_ALLOCATE(xglobal(1:mesh%np_part_global, 1:mesh%box%dim))
  !$omp parallel do
  do ip = 1, mesh%np_part_global
    xglobal(ip, 1:mesh%box%dim) = mesh_x_global(mesh, ip)
  end do

  if (mesh%parallel_in_domains) then
    SAFE_ALLOCATE(f_global(1:mesh%np_global))
    call par_vec_allgather(mesh%pv, f_global, ff)
  else
    f_global => ff
  end if

  call qshep_init(interp, mesh%np_global, f_global, xglobal(:, 1), xglobal(:, 2))
  do iu = line%nu, line%mu
    pp(1) = line%origin(1) + iu * line%spacing * line%u(1)
    pp(2) = line%origin(2) + iu * line%spacing * line%u(2)
    f_in_line(iu) = qshep_interpolate(interp, f_global, pp(1:2))
  end do
  call qshep_end(interp)

  SAFE_DEALLOCATE_A(xglobal)
  if (mesh%parallel_in_domains) then
    SAFE_DEALLOCATE_P(f_global)
  end if

  POP_SUB(X(mf_interpolate_on_line))
end subroutine X(mf_interpolate_on_line)

! ---------------------------------------------------------
!> This subroutine calculates the surface integral of a scalar
!! function on a given plane.
R_TYPE function X(mf_surface_integral_scalar) (mesh, ff, plane) result(dd)
  class(mesh_t),      intent(in) :: mesh
  R_TYPE,             intent(in) :: ff(:)  !< (mesh%np)
  type(mesh_plane_t), intent(in) :: plane

  R_TYPE, allocatable :: f_in_plane(:, :)

  PUSH_SUB(X(mf_surface_integral_scalar))

  if (mesh%box%dim /= 3) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality.'
    call messages_fatal(1)
  end if

  SAFE_ALLOCATE(f_in_plane(plane%nu:plane%mu, plane%nv:plane%mv))

  call X(mf_interpolate_on_plane)(mesh, plane, ff, f_in_plane)

  dd = sum(f_in_plane(:, :) * plane%spacing**2)

  SAFE_DEALLOCATE_A(f_in_plane)
  POP_SUB(X(mf_surface_integral_scalar))
end function X(mf_surface_integral_scalar)


! ---------------------------------------------------------
!> This subroutine calculates the surface integral of a vector
!! function on a given plane.
R_TYPE function X(mf_surface_integral_vector) (mesh, ff, plane) result(dd)
  class(mesh_t),intent(in)       :: mesh
  R_TYPE,       intent(in)       :: ff(:, :)  !< (mesh%np, mesh%box%dim)
  type(mesh_plane_t), intent(in) :: plane

  R_TYPE, allocatable :: fn(:)
  integer :: ip

  PUSH_SUB(X(mf_surface_integral_vector))

  ASSERT(not_in_openmp())

  SAFE_ALLOCATE(fn(1:mesh%np))
  !$omp parallel do
  do ip = 1, mesh%np
    fn(ip) = sum(ff(ip, 1:mesh%box%dim) * plane%n(1:mesh%box%dim))
  end do

  dd =  X(mf_surface_integral_scalar)(mesh, fn, plane)

  SAFE_DEALLOCATE_A(fn)

  POP_SUB(X(mf_surface_integral_vector))
end function X(mf_surface_integral_vector)


! ---------------------------------------------------------
!> This subroutine calculates the line integral of a scalar
!! function on a given line.
R_TYPE function X(mf_line_integral_scalar) (mesh, ff, line) result(dd)
  class(mesh_t),     intent(in) :: mesh
  R_TYPE,            intent(in) :: ff(:)  !< (mesh%np)
  type(mesh_line_t), intent(in) :: line

  R_TYPE, allocatable :: f_in_line(:)

  PUSH_SUB(X(mf_line_integral_scalar))

  if (mesh%box%dim /= 2) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality.'
    call messages_fatal(1)
  end if

  SAFE_ALLOCATE(f_in_line(line%nu:line%mu))

  call X(mf_interpolate_on_line)(mesh, line, ff, f_in_line)

  dd = sum(f_in_line(:) * line%spacing)

  SAFE_DEALLOCATE_A(f_in_line)
  POP_SUB(X(mf_line_integral_scalar))
end function X(mf_line_integral_scalar)


! ---------------------------------------------------------
!> This subroutine calculates the line integral of a vector
!! function on a given line.
R_TYPE function X(mf_line_integral_vector) (mesh, ff, line) result(dd)
  class(mesh_t),     intent(in) :: mesh
  R_TYPE,            intent(in) :: ff(:, :)  !< (mesh%np, mesh%box%dim)
  type(mesh_line_t), intent(in) :: line

  R_TYPE, allocatable :: fn(:)
  integer :: ip

  PUSH_SUB(X(mf_line_integral_vector))

  ASSERT(not_in_openmp())

  SAFE_ALLOCATE(fn(1:mesh%np))
  !$omp parallel do
  do ip = 1, mesh%np
    fn(ip) = sum(ff(ip, 1:mesh%box%dim) * line%n(1:mesh%box%dim))
  end do

  dd = X(mf_line_integral_scalar)(mesh, fn, line)

  SAFE_DEALLOCATE_A(fn)
  POP_SUB(X(mf_line_integral_vector))
end function X(mf_line_integral_vector)


! -----------------------------------------------------------------------------
!> This routine calculates the multipoles of a function ff,
!! defined in the following way:
!! multipole(1) is the trace of ff (defined to be positive; integral
!!   of ff).
!! multipole(2:4) contains the dipole: integral of ff times x, y or z.
!! multipole(5:9, is) contains the quadrupole, defined in the usual way using
!!   the spherical harmonics: multipole(5) = Integral [ ff * Y_{2,-2} ],
!!   multipole(6, is) = Integral [ f * Y_{2, -1} ].
!! And so on.
!! -----------------------------------------------------------------------------
subroutine X(mf_multipoles) (mesh, ff, lmax, multipole, mask)
  class(mesh_t),     intent(in)  :: mesh
  R_TYPE,            intent(in)  :: ff(:)
  integer,           intent(in)  :: lmax
  R_TYPE,            intent(out) :: multipole(:) !< ((lmax + 1)**2)
  logical, optional, intent(in)  :: mask(:) !< (mesh%np)

  integer :: idim, ip, ll, lm, add_lm
  FLOAT   :: xx(mesh%box%dim), rr, ylm
  R_TYPE, allocatable :: ff2(:)

  PUSH_SUB(X(mf_multipoles))

  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)

  SAFE_ALLOCATE(ff2(1:mesh%np))

  ff2(1:mesh%np) = ff(1:mesh%np)
  multipole(1) = X(mf_integrate)(mesh, ff2, mask = mask)

  if (lmax > 0) then
    do idim = 1, mesh%box%dim
      ff2(1:mesh%np) = ff(1:mesh%np) * mesh%x(1:mesh%np, idim)
      multipole(1 + idim) = X(mf_integrate)(mesh, ff2, mask = mask)
    end do
  end if

  if (lmax > 1) then
    if (mesh%box%dim /= 3) then
      message(1) = "multipoles for l > 1 are only available in 3D."
      call messages_fatal(1)
    end if
    add_lm = 5
    do ll = 2, lmax
      do lm = -ll, ll
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, coords=xx)
          call loct_ylm(1, xx(1), xx(2), xx(3), ll, lm, ylm)
          ff2(ip) = ff(ip) * ylm * rr**ll
        end do
        multipole(add_lm) = X(mf_integrate)(mesh, ff2, mask = mask)
        add_lm = add_lm + 1
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(ff2)
  POP_SUB(X(mf_multipoles))
end subroutine X(mf_multipoles)

! -----------------------------------------------------------------------------
!> This routine calculates the dipole of a function ff, for arbitrary dimensions
subroutine X(mf_dipole) (mesh, ff, dipole, mask)
  class(mesh_t),     intent(in)  :: mesh
  R_TYPE,            intent(in)  :: ff(:)
  R_TYPE,            intent(out) :: dipole(:) !< (mesh%box%dim)
  logical, optional, intent(in)  :: mask(:)   !< (mesh%np)

  integer :: idim
  R_TYPE, allocatable :: ff2(:)

  PUSH_SUB(X(mf_dipole))

  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)

  SAFE_ALLOCATE(ff2(1:mesh%np))

  do idim = 1, mesh%box%dim
    ff2(1:mesh%np) = ff(1:mesh%np) * mesh%x(1:mesh%np, idim)
    dipole(idim) = X(mf_integrate)(mesh, ff2, mask = mask)
  end do

  SAFE_DEALLOCATE_A(ff2)
  POP_SUB(X(mf_dipole))
end subroutine X(mf_dipole)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
