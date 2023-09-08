!! Copyright (C) 2013 M. Oliveira
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

! -----------------------------------------------------------------------
!> This module contains interfaces for METIS and PARMETIS routines
! -----------------------------------------------------------------------

module metis_oct_m
  use kind_oct_m
  use mpi_oct_m
  implicit none

  public

  ! Options codes (copied from metis.h)
  integer, parameter ::          &
    METIS_OPTION_PTYPE     = 1,  &
    METIS_OPTION_OBJTYPE   = 2,  &
    METIS_OPTION_CTYPE     = 3,  &
    METIS_OPTION_IPTYPE    = 4,  &
    METIS_OPTION_RTYPE     = 5,  &
    METIS_OPTION_DBGLVL    = 6,  &
    METIS_OPTION_NITER     = 7,  &
    METIS_OPTION_NCUTS     = 8,  &
    METIS_OPTION_SEED      = 9,  &
    METIS_OPTION_NO2HOP    = 10, &
    METIS_OPTION_MINCONN   = 11, &
    METIS_OPTION_CONTIG    = 12, &
    METIS_OPTION_COMPRESS  = 13, &
    METIS_OPTION_CCORDER   = 14, &
    METIS_OPTION_PFACTOR   = 15, &
    METIS_OPTION_NSEPS     = 16, &
    METIS_OPTION_UFACTOR   = 17, &
    METIS_OPTION_NUMBERING = 18, &
    METIS_OPTION_HELP      = 19, &
    METIS_OPTION_TPWGTS    = 20, &
    METIS_OPTION_NCOMMON   = 21, &
    METIS_OPTION_NOOUTPUT  = 22, &
    METIS_OPTION_BALANCE   = 23, &
    METIS_OPTION_GTYPE     = 24, &
    METIS_OPTION_UBVEC     = 25

  !Error code from metis.h
  integer, parameter ::         &
    METIS_OK              = 1,  & !< Returned normally
    METIS_ERROR_INPUT     = -2, & !< Returned due to erroneous inputs and/or options
    METIS_ERROR_MEMORY    = -3, & !< Returned due to insufficient memory
    METIS_ERROR           = -4    !< Some other errors

  ! define the integer kind as needed for the metis installation
  integer, parameter :: &
#if METIS_IDXTYPEWIDTH == 64
    imetis = i8, &
    MPI_METIS_INT = MPI_INTEGER8
#else
  imetis = i4, &
    MPI_METIS_INT = MPI_INTEGER
#endif

  interface

#if defined(HAVE_PARMETIS) || defined(HAVE_METIS)
    subroutine oct_metis_setdefaultoptions(options)
      import :: imetis
      implicit none
      integer(imetis), intent(inout) :: options
    end subroutine oct_metis_setdefaultoptions

    integer function oct_metis_partgraphrecursive(nvtxs, ncon, xadj, adjncy, nparts, tpwgts, ubvec, options, objval, part)
      import :: imetis
      implicit none
      integer(imetis), intent(in)  :: nvtxs  !< The number of vertices in the graph.
      integer(imetis), intent(in)  :: ncon   !< The number of balancing constraints. It should be at least 1.
      integer(imetis), intent(in)  :: xadj   !< The adjacency structure of the graph.
      integer(imetis), intent(in)  :: adjncy !< The adjacency structure of the graph.
      integer(imetis), intent(in)  :: nparts !< The number of parts to partition the graph.
      REAL_SINGLE,     intent(in)  :: tpwgts !< This is an array of size nparts x ncon that specifies the desired
      !!                                        weight for each partition and constraint.
      REAL_SINGLE,     intent(in)  :: ubvec  !< This is an array of size ncon that specifies the allowed load imbalance
      !!                                        tolerance for each constraint.
      integer(imetis), intent(in)  :: options!< This is the array of options.
      integer(imetis), intent(out) :: objval !< Upon successful completion, this variable stores the edge-cut or the total
      !!                                        communication volume of the partitioning solution.
      integer(imetis), intent(out) :: part   !< This is a vector of size nvtxs that upon successful completion stores the
      !!                                        partition vector of the graph.
    end function oct_metis_partgraphrecursive

    integer function oct_metis_partgraphkway(nvtxs, ncon, xadj, adjncy, nparts, tpwgts, ubvec, options, objval, part)
      import :: imetis
      implicit none
      integer(imetis), intent(in)  :: nvtxs  !< The number of vertices in the graph.
      integer(imetis), intent(in)  :: ncon   !< The number of balancing constraints. It should be at least 1.
      integer(imetis), intent(in)  :: xadj   !< The adjacency structure of the graph.
      integer(imetis), intent(in)  :: adjncy !< The adjacency structure of the graph.
      integer(imetis), intent(in)  :: nparts !< The number of parts to partition the graph.
      REAL_SINGLE,     intent(in)  :: tpwgts !< This is an array of size nparts x ncon that specifies the desired weight for
      !!                                        each partition and constraint.
      REAL_SINGLE,     intent(in)  :: ubvec  !< This is an array of size ncon that specifies the allowed load imbalance
      !!                                        tolerance for each constraint.
      integer(imetis), intent(in)  :: options!< This is the array of options.
      integer(imetis), intent(out) :: objval !< Upon successful completion, this variable stores the edge-cut or the total
      !!                                        communication volume of the partitioning solution.
      integer(imetis), intent(out) :: part   !< This is a vector of size nvtxs that upon successful completion stores the
      !!                                        partition vector of the graph.
    end function oct_metis_partgraphkway

#endif
#if defined(HAVE_PARMETIS)

    subroutine oct_parmetis_v3_partkway(vtxdist, xadj, adjncy, ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm)
      import :: imetis
      implicit none
      integer(imetis), intent(in)  :: vtxdist!< This array describes how the vertices of the graph are distributed among
      !!                                        the processors. Its contents are identical for every processor.
      integer(imetis), intent(in)  :: xadj   !< These store the (local) adjacency structure of the graph at each processor.
      integer(imetis), intent(in)  :: adjncy !< These store the (local) adjacency structure of the graph at each processor.
      integer(imetis), intent(in)  :: ncon   !< This is used to specify the number of weights that each vertex has. It is
      !!                                        also the number of balance constraints that must be satisfied.
      integer(imetis), intent(in)  :: nparts !< This is used to specify the number of sub-domains that are desired. Note
      !!                                        that the number of sub-domains is independent of the number of processors
      !!                                        that call this routine.
      REAL_SINGLE,     intent(in)  :: tpwgts !< An array of size ncon x nparts that is used to specify the fraction of vertex
      !!                                        weight that should be distributed to each sub-domain for each balance
      !!                                        constraint. If all of the sub-domains are to be of the same size for every
      !!                                        vertex weight, then each of the ncon x nparts elements should be set to
      !!                                        a value of 1/nparts. If ncon is greater than 1, the target sub-domain weights
      !!                                        for each sub-domain are stored contiguously (similar to the vwgt array).
      !!                                        Note that the sum of all of the tpwgts for a given vertex weight should be one.
      REAL_SINGLE,     intent(in)  :: ubvec  !< An array of size ncon that is used to specify the imbalance tolerance for each
      !!                                        vertex weight, with 1 being perfect balance and nparts being perfect imbalance.
      !!                                        A value of 1.05 for each of the ncon weights is recommended.
      integer(imetis), intent(in)  :: options!< This is an array of integers that is used to pass additional parameters for
      !!                                        the routine. The first element (i.e., options[0]) can take either the value of
      !!                                        0 or 1. If it is 0, then the default values are used.
      integer(imetis), intent(out) :: edgecut!< Upon successful completion, the number of edges that are cut by the
      !!                                        partitioning is written to this parameter.
      integer(imetis), intent(out) :: part   !< This is an array of size equal to the number of locally-stored vertices. Upon
      !!                                        successful completion the partition vector of the locally-stored vertices is
      !!                                        written to this array.
      integer,         intent(in)  :: comm   !< This is a pointer to the MPI communicator of the processes that call PARMETIS.
    end subroutine oct_parmetis_v3_partkway

#endif
  end interface

  interface i4_to_imetis
    module procedure i4_to_imetis_0, i4_to_imetis_1
  end interface i4_to_imetis

  interface imetis_to_i4
    module procedure imetis_to_i4_0, imetis_to_i4_1
  end interface imetis_to_i4

  interface i8_to_imetis
    module procedure i8_to_imetis_0, i8_to_imetis_1
  end interface i8_to_imetis

  interface imetis_to_i8
    module procedure imetis_to_i8_0, imetis_to_i8_1
  end interface imetis_to_i8

contains

  integer(imetis) pure function i4_to_imetis_0(ii)
    integer(i4), intent(in) :: ii

    i4_to_imetis_0 = int(ii, imetis)
  end function i4_to_imetis_0

  integer(i4) pure function imetis_to_i4_0(ii)
    integer(imetis), intent(in) :: ii

    imetis_to_i4_0 = int(ii, i4)
  end function imetis_to_i4_0

  pure function i4_to_imetis_1(ii)
    integer(i4), intent(in) :: ii(:)
    integer(imetis) :: i4_to_imetis_1(lbound(ii, 1):ubound(ii, 1))

    i4_to_imetis_1 = int(ii, imetis)
  end function i4_to_imetis_1

  pure function imetis_to_i4_1(ii)
    integer(imetis), intent(in) :: ii(:)
    integer(i4) :: imetis_to_i4_1(lbound(ii, 1, kind=imetis):ubound(ii, 1, kind=imetis))

    imetis_to_i4_1 = int(ii, i4)
  end function imetis_to_i4_1

  integer(imetis) pure function i8_to_imetis_0(ii)
    integer(i8), intent(in) :: ii

    i8_to_imetis_0 = int(ii, imetis)
  end function i8_to_imetis_0

  integer(i8) pure function imetis_to_i8_0(ii)
    integer(imetis), intent(in) :: ii

    imetis_to_i8_0 = int(ii, i8)
  end function imetis_to_i8_0

  pure function i8_to_imetis_1(ii)
    integer(i8), intent(in) :: ii(:)
    integer(imetis) :: i8_to_imetis_1(lbound(ii, 1, kind=i8):ubound(ii, 1, kind=i8))

    i8_to_imetis_1 = int(ii, imetis)
  end function i8_to_imetis_1

  pure function imetis_to_i8_1(ii)
    integer(imetis), intent(in) :: ii(:)
    integer(i8) :: imetis_to_i8_1(lbound(ii, 1, kind=imetis):ubound(ii, 1, kind=imetis))

    imetis_to_i8_1 = int(ii, i8)
  end function imetis_to_i8_1
end module metis_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
