!! Copyright (C) 2009 N. Helbig and M. Verstraete
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

module modelmb_density_matrix_oct_m
  use batch_oct_m
  use comm_oct_m
  use debug_oct_m
  use global_oct_m
  use hypercube_oct_m
  use io_oct_m
  use index_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use modelmb_1part_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_oct_m

  implicit none

  private

  public ::                        &
    zmodelmb_density_matrix_write, &
    dmodelmb_density_matrix_write, &
    modelmb_density_matrix_init,   &
    modelmb_density_matrix_end,    &
    modelmb_denmat_t

  type modelmb_denmat_t
    private
    integer :: ndensmat_to_calculate
    character(len=200) :: dirname
    character(80), allocatable :: labels(:)
    integer, allocatable :: particle_kept(:)
    integer, allocatable :: nnatorb_prt(:)
  end type modelmb_denmat_t

contains

  subroutine modelmb_density_matrix_init(dir, namespace, st, denmat)
    character(len=*),       intent(in)  :: dir
    type(namespace_t),      intent(in)  :: namespace
    type(states_elec_t),    intent(in)  :: st
    type(modelmb_denmat_t), intent(out) :: denmat

    integer :: ncols, ipart
    type(block_t) :: blk

    PUSH_SUB(modelmb_density_matrix_init)

    !%Variable DensitytoCalc
    !%Type block
    !%Section States::ModelMB
    !%Description
    !% Choice of which particle density (event. matrices) will be calculated and output, in the
    !% modelmb particles scheme.
    !%
    !% <tt>%DensitytoCalc
    !% <br>&nbsp;&nbsp; "proton"   | 1 | 10
    !% <br>&nbsp;&nbsp; "electron" | 2 | 15
    !% <br>%</tt>
    !%
    !% would ask octopus to calculate the density matrix corresponding to the 1st
    !% particle (whose coordinates correspond to dimensions 1 to ndim_modelmb),
    !% which is an proton, then that corresponding to the 2nd particle
    !% (electron with dimensions ndim_modelmb+1 to 2*ndim_modelmb), printing
    !% 10 natural orbitals for the first and 15 for the second.
    !%
    !% <tt>%DensitytoCalc
    !% <br>&nbsp;&nbsp; "proton"   | 1 | -1
    !% <br>&nbsp;&nbsp; "electron" | 2 | -1
    !% <br>%</tt>
    !%
    !% would ask octopus to print out just the densities for particles 1 and 2
    !% without any density matrix output.
    !%
    !%End

    call messages_obsolete_variable(namespace, 'DensityMatrixtoCalc', 'DensitytoCalc')
    call messages_obsolete_variable(namespace, 'DensitiestoCalc', 'DensitytoCalc')

    if (parse_block(namespace, 'DensitytoCalc', blk) /= 0) then
      message(1) = 'To print out density (matrices), you must specify the DensitytoCalc block in input'
      call messages_fatal(1, namespace=namespace)
    end if

    ncols = parse_block_cols(blk, 0)
    if (ncols /= 3) then
      call messages_input_error(namespace, "DensitytoCalc")
    end if
    denmat%ndensmat_to_calculate=parse_block_n(blk)
    if (denmat%ndensmat_to_calculate < 0 .or. &
      denmat%ndensmat_to_calculate > st%modelmbparticles%nparticle) then
      call messages_input_error(namespace, "DensitytoCalc")
    end if

    SAFE_ALLOCATE(denmat%labels(1:denmat%ndensmat_to_calculate))
    SAFE_ALLOCATE(denmat%particle_kept(1:denmat%ndensmat_to_calculate))
    SAFE_ALLOCATE(denmat%nnatorb_prt(1:denmat%ndensmat_to_calculate))

    do ipart=1,denmat%ndensmat_to_calculate
      call parse_block_string(blk, ipart-1, 0, denmat%labels(ipart))
      call parse_block_integer(blk, ipart-1, 1, denmat%particle_kept(ipart))
      call parse_block_integer(blk, ipart-1, 2, denmat%nnatorb_prt(ipart))

      write (message(1),'(a,a)') 'labels_densmat = ', denmat%labels(ipart)
      write (message(2),'(a,i6)') 'particle_kept_densmat = ', denmat%particle_kept(ipart)
      write (message(3),'(a,i6)') 'nnatorb_prt_densmat = ', denmat%nnatorb_prt(ipart)
      call messages_info(3, namespace=namespace)
    end do
    call parse_block_end(blk)
    ! END reading in of input var block DensitytoCalc

    denmat%dirname = trim(dir)

    POP_SUB(modelmb_density_matrix_init)
  end subroutine modelmb_density_matrix_init

  ! ---------------------------------------------------------
  subroutine modelmb_density_matrix_end(this)
    type(modelmb_denmat_t), intent(inout) :: this

    PUSH_SUB(modelmb_density_matrix_end)

    SAFE_DEALLOCATE_A(this%labels)
    SAFE_DEALLOCATE_A(this%particle_kept)
    SAFE_DEALLOCATE_A(this%nnatorb_prt)

    POP_SUB(modelmb_density_matrix_end)
  end subroutine modelmb_density_matrix_end


#include "undef.F90"
#include "real.F90"
#include "modelmb_density_matrix_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "modelmb_density_matrix_inc.F90"
#include "undef.F90"

end module modelmb_density_matrix_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
