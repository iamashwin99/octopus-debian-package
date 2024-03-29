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

module epot_oct_m
  use comm_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use ion_interaction_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use ps_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use splines_oct_m
  use spline_filter_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                        &
    epot_t,                        &
    epot_init,                     &
    epot_end,                      &
    epot_generate,                 &
    epot_local_potential,          &
    epot_precalc_local_potential,  &
    epot_have_external_potentials, &
    local_potential_has_density

  integer, public, parameter :: &
    NOREL      = 0,             &
    SPIN_ORBIT = 1

  type epot_t
    ! Components are public by default

    ! Ions
    FLOAT,             allocatable :: vpsl(:)       !< the local part of the pseudopotentials
    !                                               !< plus the potential from static electric fields
    type(projector_t), allocatable :: proj(:)       !< non-local projectors
    logical                        :: non_local
    integer                        :: natoms

    ! External e-m fields
    FLOAT, allocatable     :: E_field(:)           !< static electric field
    FLOAT, allocatable     :: v_ext(:)             !< static scalar potential - 1:gr%np_part
    FLOAT, allocatable     :: B_field(:)           !< static magnetic field
    FLOAT, allocatable     :: A_static(:,:)        !< static vector potential
    integer                :: reltype              !< type of relativistic correction to use

    !> The gyromagnetic ratio (-2.0 for the electron, but different if we treat
    !! *effective* electrons in a quantum dot. It affects the spin Zeeman term.)
    FLOAT :: gyromagnetic_ratio

    !> SO prefactor (1.0 = normal SO, 0.0 = no SO)
    FLOAT, private :: so_strength

    !> the ion-ion energy and force
    FLOAT              :: eii
    FLOAT, allocatable :: fii(:, :)
    FLOAT, allocatable :: vdw_forces(:, :)
    FLOAT, allocatable :: photon_forces(:)

    FLOAT, allocatable, private :: local_potential(:,:)
    logical,            private :: local_potential_precalculated

    logical,                  private :: have_density
    type(poisson_t), pointer, private :: poisson_solver

    logical              :: nlcc = .false.   !< does any species have non-local core corrections?
  end type epot_t

contains

  ! ---------------------------------------------------------
  subroutine epot_init(ep, namespace, space, gr, ions, psolver, ispin, xc_family, kpoints)
    type(epot_t),                       intent(out)   :: ep
    type(namespace_t),                  intent(in)    :: namespace
    type(space_t),                      intent(in)    :: space
    type(grid_t),                       intent(in)    :: gr
    type(ions_t),                       intent(inout) :: ions
    type(poisson_t),  target,           intent(in)    :: psolver
    integer,                            intent(in)    :: ispin
    integer,                            intent(in)    :: xc_family
    type(kpoints_t),                    intent(in)    :: kpoints


    integer :: ispec, ia
    integer :: filter

    PUSH_SUB(epot_init)

    !%Variable FilterPotentials
    !%Type integer
    !%Default filter_ts
    !%Section Hamiltonian
    !%Description
    !% <tt>Octopus</tt> can filter the pseudopotentials so that they no
    !% longer contain Fourier components larger than the mesh itself. This is
    !% very useful to decrease the egg-box effect, and so should be used in
    !% all instances where atoms move (<i>e.g.</i> geometry optimization,
    !% molecular dynamics, and vibrational modes).
    !%Option filter_none 0
    !% Do not filter.
    !%Option filter_TS 2
    !% The filter of M. Tafipolsky and R. Schmid, <i>J. Chem. Phys.</i> <b>124</b>, 174102 (2006).
    !%Option filter_BSB 3
    !% The filter of E. L. Briggs, D. J. Sullivan, and J. Bernholc, <i>Phys. Rev. B</i> <b>54</b>, 14362 (1996).
    !%End
    call parse_variable(namespace, 'FilterPotentials', PS_FILTER_TS, filter)
    if (.not. varinfo_valid_option('FilterPotentials', filter)) call messages_input_error(namespace, 'FilterPotentials')
    call messages_print_var_option("FilterPotentials", filter, namespace=namespace)

    if (family_is_mgga(xc_family) .and. filter /= PS_FILTER_NONE) then
      call messages_not_implemented("FilterPotentials different from filter_none with MGGA", namespace=namespace)
    end if

    if (filter == PS_FILTER_TS) call spline_filter_mask_init()
    do ispec = 1, ions%nspecies
      call species_pot_init(ions%species(ispec), namespace, mesh_gcutoff(gr), filter)
    end do

    SAFE_ALLOCATE(ep%vpsl(1:gr%np))

    ep%vpsl(1:gr%np) = M_ZERO

    ! No more "UserDefinedTDPotential" from this version on.
    call messages_obsolete_variable(namespace, 'UserDefinedTDPotential', 'TDExternalFields')

    call messages_obsolete_variable(namespace, 'ClassicalPotential')

    !%Variable GyromagneticRatio
    !%Type float
    !%Default 2.0023193043768
    !%Section Hamiltonian
    !%Description
    !% The gyromagnetic ratio of the electron. This is of course a physical
    !% constant, and the default value is the exact one that you should not
    !% touch, unless:
    !% (i)  You want to disconnect the anomalous Zeeman term in the Hamiltonian
    !% (then set it to zero; this number only affects that term);
    !% (ii) You are using an effective Hamiltonian, as is the case when
    !% you calculate a 2D electron gas, in which case you have an effective
    !% gyromagnetic factor that depends on the material.
    !%End
    call parse_variable(namespace, 'GyromagneticRatio', P_g, ep%gyromagnetic_ratio)

    !%Variable RelativisticCorrection
    !%Type integer
    !%Default non_relativistic
    !%Section Hamiltonian
    !%Description
    !% The default value means that <i>no</i> relativistic correction is used. To
    !% include spin-orbit coupling turn <tt>RelativisticCorrection</tt> to <tt>spin_orbit</tt>
    !% (this will only work if <tt>SpinComponents</tt> has been set to <tt>non_collinear</tt>, which ensures
    !% the use of spinors).
    !%Option non_relativistic 0
    !% No relativistic corrections.
    !%Option spin_orbit 1
    !% Spin-orbit.
    !%End
    call parse_variable(namespace, 'RelativisticCorrection', NOREL, ep%reltype)
    if (.not. varinfo_valid_option('RelativisticCorrection', ep%reltype)) then
      call messages_input_error(namespace, 'RelativisticCorrection')
    end if
    if (ispin /= SPINORS .and. ep%reltype == SPIN_ORBIT) then
      message(1) = "The spin-orbit term can only be applied when using spinors."
      call messages_fatal(1, namespace=namespace)
    end if

    if(ep%reltype == SPIN_ORBIT .and. kpoints%use_symmetries) then
      call messages_not_implemented("Spin-orbit coupling and k-point symmetries", namespace=namespace)
    end if

    call messages_print_var_option("RelativisticCorrection", ep%reltype, namespace=namespace)

    !%Variable SOStrength
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian
    !%Description
    !% Tuning of the spin-orbit coupling strength: setting this value to zero turns off spin-orbit terms in
    !% the Hamiltonian, and setting it to one corresponds to full spin-orbit.
    !%End
    if (ep%reltype == SPIN_ORBIT) then
      call parse_variable(namespace, 'SOStrength', M_ONE, ep%so_strength)
    else
      ep%so_strength = M_ONE
    end if

    SAFE_ALLOCATE(ep%proj(1:ions%natoms))

    ep%natoms = ions%natoms
    ep%non_local = .false.

    ep%eii = M_ZERO
    SAFE_ALLOCATE(ep%fii(1:ions%space%dim, 1:ions%natoms))
    ep%fii = M_ZERO

    SAFE_ALLOCATE(ep%vdw_forces(1:ions%space%dim, 1:ions%natoms))
    ep%vdw_forces = M_ZERO

    SAFE_ALLOCATE(ep%photon_forces(1:ions%space%dim))
    ep%photon_forces = M_ZERO

    ep%local_potential_precalculated = .false.


    ep%have_density = .false.
    do ia = 1, ions%nspecies
      if (local_potential_has_density(ions%space, ions%species(ia))) then
        ep%have_density = .true.
        exit
      end if
    end do

    if (ep%have_density) then
      ep%poisson_solver => psolver
    else
      nullify(ep%poisson_solver)
    end if

    ! find out if we need non-local core corrections
    ep%nlcc = .false.
    do ia = 1, ions%nspecies
      ep%nlcc = (ep%nlcc.or.species_has_nlcc(ions%species(ia)))
    end do

    POP_SUB(epot_init)
  end subroutine epot_init

  ! ---------------------------------------------------------
  subroutine epot_end(ep)
    type(epot_t), intent(inout) :: ep

    integer :: iproj

    PUSH_SUB(epot_end)

    if (ep%have_density) then
      nullify(ep%poisson_solver)
    end if

    SAFE_DEALLOCATE_A(ep%local_potential)
    SAFE_DEALLOCATE_A(ep%fii)
    SAFE_DEALLOCATE_A(ep%vdw_forces)
    SAFE_DEALLOCATE_A(ep%vpsl)
    SAFE_DEALLOCATE_A(ep%photon_forces)

    ! the macroscopic fields
    SAFE_DEALLOCATE_A(ep%E_field)
    SAFE_DEALLOCATE_A(ep%v_ext)
    SAFE_DEALLOCATE_A(ep%B_field)
    SAFE_DEALLOCATE_A(ep%A_static)

    do iproj = 1, ep%natoms
      if (projector_is_null(ep%proj(iproj))) cycle
      call projector_end(ep%proj(iproj))
    end do

    ASSERT(allocated(ep%proj))
    SAFE_DEALLOCATE_A(ep%proj)

    POP_SUB(epot_end)

  end subroutine epot_end

  ! ---------------------------------------------------------
  subroutine epot_generate(ep, namespace, mesh, ions, st_d)
    type(epot_t),             intent(inout) :: ep
    type(namespace_t),        intent(in)    :: namespace
    class(mesh_t),    target, intent(in)    :: mesh
    type(ions_t),     target, intent(inout) :: ions
    type(states_elec_dim_t),  intent(inout) :: st_d

    integer :: ia
    type(profile_t), save :: epot_generate_prof
    type(ps_t), pointer :: ps

    call profiling_in(epot_generate_prof, "EPOT_GENERATE")
    PUSH_SUB(epot_generate)

    ! Local part
    ep%vpsl = M_ZERO

    ! we assume that we need to recalculate the ion-ion energy
    call ion_interaction_calculate(ions%ion_interaction, ions%space, ions%latt, ions%atom, &
      ions%natoms, ions%pos, mesh%box%bounding_box_l, ep%eii, ep%fii)

    ! the pseudopotential part.
    do ia = 1, ions%natoms
      if (.not. species_is_ps(ions%atom(ia)%species)) cycle
      call projector_end(ep%proj(ia))
      call projector_init(ep%proj(ia), ions%atom(ia), namespace, st_d%dim, ep%reltype)
    end do

    do ia = ions%atoms_dist%start, ions%atoms_dist%end
      if (ep%proj(ia)%type == PROJ_NONE) cycle
      ps => species_ps(ions%atom(ia)%species)
      call submesh_init(ep%proj(ia)%sphere, ions%space, mesh, ions%latt, ions%pos(:, ia), ps%rc_max)
    end do

    if (ions%atoms_dist%parallel) then
      do ia = 1, ions%natoms
        if (ep%proj(ia)%type == PROJ_NONE) cycle
        ps => species_ps(ions%atom(ia)%species)
        call submesh_broadcast(ep%proj(ia)%sphere, ions%space, mesh, ions%pos(:, ia), ps%rc_max, &
          ions%atoms_dist%node(ia), ions%atoms_dist%mpi_grp)
      end do
    end if

    do ia = 1, ions%natoms
      call projector_build(ep%proj(ia), ions%atom(ia), ep%so_strength)
      if (.not. projector_is(ep%proj(ia), PROJ_NONE)) ep%non_local = .true.
    end do

    POP_SUB(epot_generate)
    call profiling_out(epot_generate_prof)
  end subroutine epot_generate

  ! ---------------------------------------------------------

  logical pure function local_potential_has_density(space, species) result(has_density)
    type(space_t),        intent(in)    :: space
    type(species_t),      intent(in)    :: species

    has_density = species_has_density(species) .or. (species_is_ps(species) .and. space%is_periodic())

  end function local_potential_has_density

  ! ---------------------------------------------------------
  subroutine epot_local_potential(ep, namespace, space, latt, mesh, species, pos, iatom, vpsl)
    type(epot_t),             intent(in)    :: ep
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(lattice_vectors_t),  intent(in)    :: latt
    class(mesh_t),            intent(in)    :: mesh
    type(species_t),          intent(in)    :: species
    FLOAT,                    intent(in)    :: pos(1:space%dim)
    integer,                  intent(in)    :: iatom
    FLOAT,                    intent(inout) :: vpsl(:)

    integer :: ip
    FLOAT :: radius
    FLOAT, allocatable :: vl(:), rho(:)
    type(submesh_t)  :: sphere
    type(profile_t), save :: prof
    type(ps_t), pointer :: ps

    PUSH_SUB(epot_local_potential)
    call profiling_in(prof, "EPOT_LOCAL")

    if (ep%local_potential_precalculated) then

      do ip = 1, mesh%np
        vpsl(ip) = vpsl(ip) + ep%local_potential(ip, iatom)
      end do
    else

      !Local potential, we can get it by solving the Poisson equation
      !(for all-electron species or pseudopotentials in periodic
      !systems) or by applying it directly to the grid
      SAFE_ALLOCATE(vl(1:mesh%np))

      if (local_potential_has_density(space, species)) then
        SAFE_ALLOCATE(rho(1:mesh%np))

        call species_get_long_range_density(species, namespace, space, latt, pos, mesh, rho, sphere)

        if (poisson_solver_is_iterative(ep%poisson_solver)) then
          ! vl has to be initialized before entering routine
          ! and our best guess for the potential is zero
          vl(1:mesh%np) = M_ZERO
        end if

        call dpoisson_solve(ep%poisson_solver, namespace, vl, rho, all_nodes = .false.)

        SAFE_DEALLOCATE_A(rho)

      else

        call species_get_local(species, namespace, space, latt, pos, mesh, vl)

      end if

      call lalg_axpy(mesh%np, M_ONE, vl, vpsl)
      SAFE_DEALLOCATE_A(vl)

      !the localized part
      if (species_is_ps(species)) then

        ps => species_ps(species)

        radius = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold)*CNST(1.05)
        if (.not. submesh_compatible(sphere, radius, pos, minval(mesh%spacing(1:space%dim)))) then
          call submesh_end(sphere)
          call submesh_init(sphere, space, mesh, latt, pos, radius)
        end if
        SAFE_ALLOCATE(vl(1:sphere%np))
        vl = M_ZERO

        do ip = 1, sphere%np
          if(sphere%r(ip) <= radius) then
            vl(ip) = spline_eval(ps%vl, sphere%r(ip))
          end if
        end do

        call submesh_add_to_mesh(sphere, vl, vpsl)

        SAFE_DEALLOCATE_A(vl)
        nullify(ps)

      end if
      call submesh_end(sphere)

    end if

    call profiling_out(prof)
    POP_SUB(epot_local_potential)
  end subroutine epot_local_potential

  ! ---------------------------------------------------------
  subroutine epot_precalc_local_potential(ep, namespace, gr, ions)
    type(epot_t),      intent(inout) :: ep
    type(namespace_t), intent(in)    :: namespace
    type(grid_t),      intent(in)    :: gr
    type(ions_t),      intent(in)    :: ions

    integer :: iatom

    PUSH_SUB(epot_precalc_local_potential)

    if (.not. allocated(ep%local_potential)) then
      SAFE_ALLOCATE(ep%local_potential(1:gr%np, 1:ions%natoms))
    end if

    ep%local_potential_precalculated = .false.

    do iatom = 1, ions%natoms
      ep%local_potential(1:gr%np, iatom) = M_ZERO
      call epot_local_potential(ep, namespace, ions%space, ions%latt, gr, ions%atom(iatom)%species, &
        ions%pos(:, iatom), iatom, ep%local_potential(1:gr%np, iatom))!, time)
    end do
    ep%local_potential_precalculated = .true.

    POP_SUB(epot_precalc_local_potential)
  end subroutine epot_precalc_local_potential

  ! ---------------------------------------------------------

  logical function epot_have_external_potentials(ep)
    type(epot_t), intent(in)  :: ep

    PUSH_SUB(epot_have_external_potentials)

    epot_have_external_potentials = allocated(ep%E_field)

    POP_SUB(epot_have_external_potentials)

  end function epot_have_external_potentials

end module epot_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
