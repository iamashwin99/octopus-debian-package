#include "global.h"

module atom_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                               &
    atom_init,                            &
    atom_end,                             &
    atom_get_label,                       &
    atom_set_species,                     &
    atom_get_species,                     &
    atom_same_species

  type, public :: atom_t
    !private
    character(len=LABEL_LEN)  :: label = ""
    type(species_t), pointer  :: species => null() !< pointer to species
    integer, dimension(MAX_DIM) :: c   = 0      !< Constrain on te atom (0 or 1)

    !Components of the force
    FLOAT, dimension(MAX_DIM) :: f_ii     = M_ZERO !< Ion-Ion part
    FLOAT, dimension(MAX_DIM) :: f_vdw    = M_ZERO !< Van der Waals part
    FLOAT, dimension(MAX_DIM) :: f_loc    = M_ZERO !< Local electronic part
    FLOAT, dimension(MAX_DIM) :: f_nl     = M_ZERO !< NL electronic part
    FLOAT, dimension(MAX_DIM) :: f_fields = M_ZERO !< Lasers
    FLOAT, dimension(MAX_DIM) :: f_u      = M_ZERO !< Hubbard forces
    FLOAT, dimension(MAX_DIM) :: f_scf    = M_ZERO !< SCF forces
    FLOAT, dimension(MAX_DIM) :: f_nlcc   = M_ZERO !< NLCC forces
    FLOAT, dimension(MAX_DIM) :: f_photons= M_ZERO !< Photons forces
  end type atom_t

  interface atom_same_species
    module procedure atom_same_species_aa
    module procedure atom_same_species_as
  end interface atom_same_species

contains

  ! ---------------------------------------------------------
  subroutine atom_init(this, label, species)
    type(atom_t),                      intent(out) :: this
    character(len=*),                  intent(in)  :: label
    type(species_t), target, optional, intent(in)  :: species

    PUSH_SUB(atom_init)

    this%label = trim(adjustl(label))
    this%species => null()
    if (present(species)) this%species => species

    this%f_ii      = M_ZERO
    this%f_vdw     = M_ZERO
    this%f_loc     = M_ZERO
    this%f_nl      = M_ZERO
    this%f_fields  = M_ZERO
    this%f_u       = M_ZERO
    this%f_photons = M_ZERO

    POP_SUB(atom_init)
  end subroutine atom_init

  ! ---------------------------------------------------------
  elemental subroutine atom_end(this)
    type(atom_t), intent(inout) :: this

    this%label = ""
    this%species => null()

    this%f_ii      = M_ZERO
    this%f_vdw     = M_ZERO
    this%f_loc     = M_ZERO
    this%f_nl      = M_ZERO
    this%f_fields  = M_ZERO
    this%f_u       = M_ZERO
    this%f_photons = M_ZERO

  end subroutine atom_end

  ! ---------------------------------------------------------

  pure function atom_get_label(this) result(label)
    type(atom_t), intent(in) :: this

    character(len=len_trim(adjustl(this%label))) :: label

    label=trim(adjustl(this%label))

  end function atom_get_label

  ! ---------------------------------------------------------
  subroutine atom_set_species(this, species)
    type(atom_t),            intent(inout) :: this
    type(species_t), target, intent(in)    :: species

    PUSH_SUB(atom_set_species)

    this%species => species
    POP_SUB(atom_set_species)

  end subroutine atom_set_species

  ! ---------------------------------------------------------
  subroutine atom_get_species(this, species)
    type(atom_t),    target,  intent(in)  :: this
    type(species_t), pointer, intent(out) :: species

    ! NO PUSH_SUB, called too often

    species => null()
    if (associated(this%species)) species => this%species

  end subroutine atom_get_species

  ! ---------------------------------------------------------
  elemental function atom_same_species_aa(this, that) result(is)
    type(atom_t), intent(in) :: this
    type(atom_t), intent(in) :: that

    logical :: is

    is = (atom_get_label(this) == atom_get_label(that))

  end function atom_same_species_aa

  ! ---------------------------------------------------------
  elemental function atom_same_species_as(this, species) result(is)
    type(atom_t),    intent(in) :: this
    type(species_t), intent(in) :: species

    logical :: is

    is = (atom_get_label(this) == species_label(species))

  end function atom_same_species_as

end module atom_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
