!! Copyright (C) 2002-2015 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

!> @brief This module handles the communicators for the various parallelization strategies.
!!
!! Given a number n_index of indices and their ranges index_range(1:n_index),
!! we divide the n_nodes in groups and create communicators for each group.
!! Each group is associated with one index. The min_range indicates the minimum
!! number of elements in each processor. For example, given min_range(1) = 25000,
!! the algorithm will try to put at least 25000 points in each processor
!!
!! Example:
!! \verbatim
!! Given 3 indices with ranges
!!   index_range(1) = 100000, (number of points in the mesh)
!!   index_range(2) = 15,     (number of states)
!!   index_range(3) = 2,      (number of k-points)
!! \endverbatim
!! and 12 processors, we could get
!! \verbatim
!!   mc%group_sizes = (2, 3, 2)
!! \endverbatim
!! which means that
!! - space is divided in 2 domains per state,
!! - the states are divided in 3 groups, i.e. 5 states per processor, and
!! - the whole setting is duplicated because of the 2 k-points.
!!
!! To perform collective operations (like a reduce), you can use the communicators
!! provided in mc%group_comm(:). For example, to sum over states, the communicator
!! to use is mc%group_comm(P_STRATEGY_STATES)
!!
!! You can use the routine multicomm_strategy_is_parallel to know if a certain
!! index is parallelized.

module multicomm_oct_m
  use blacs_oct_m
  use debug_oct_m
  use global_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use namespace_oct_m
#if defined(HAVE_OPENMP)
  use omp_lib
#endif
  use parser_oct_m
  use profiling_oct_m
  use utils_oct_m

  implicit none

  private

  public ::                          &
    multicomm_divide_range,          &
    multicomm_divide_range_omp,      &
#if defined(HAVE_MPI)
    multicomm_create_all_pairs,      &
#endif
    multicomm_t,                     &
    multicomm_all_pairs_t,           &
    multicomm_init,                  &
    multicomm_end,                   &
    multicomm_all_pairs_copy,        &
    multicomm_strategy_is_parallel,  &
    multicomm_is_slave,              &
    multicomm_have_slaves

  !> possible parallelization strategies
  integer, public, parameter ::      &
    P_STRATEGY_SERIAL  = 0,          & !< single domain, all states, k-points on a single processor
    P_STRATEGY_DOMAINS = 1,          & !< parallelization in domains
    P_STRATEGY_STATES  = 2,          & !< parallelization in states
    P_STRATEGY_KPOINTS = 3,          & !< parallelization in k-points
    P_STRATEGY_OTHER   = 4,          & !< something else like e-h pairs
    P_STRATEGY_MAX     = 4

  integer, public, parameter ::      &
    P_MASTER           = 1,          &
    P_SLAVE            = 2

  integer, public, parameter ::      &
    PAR_AUTO            = -1,        &
    PAR_NO              =  0

  integer,           parameter :: n_par_types = 4
  character(len=11), parameter :: par_types(0:n_par_types) = &
    (/                              &
    "serial    ",                   &
    "ParDomains",                   &
    "ParStates ",                   &
    "ParKPoints",                   &
    "ParOther  "                    &
    /)

  integer, parameter :: MAX_INDEX = 5

  !> Stores all communicators and groups
  type multicomm_t
    private
    integer          :: n_node           !< Total number of nodes.

    integer, public  :: par_strategy     !< What kind of parallelization strategy should we use?

    integer, allocatable         :: group_sizes(:)     !< Number of processors in each group.
    integer, allocatable, public :: who_am_i(:)        !< Rank in the "line"-communicators.
    integer, allocatable, public :: group_comm(:)      !< "Line"-communicators I belong to.
    integer, public              :: dom_st_comm        !< States-domain plane communicator.
    integer, public              :: st_kpt_comm        !< Kpoints-states plane communicator.
    integer, public              :: dom_st_kpt_comm    !< Kpoints-states-domain cube communicator.
    type(mpi_grp_t), public      :: intranode_grp      !< group for ranks on a compute node
    type(mpi_grp_t), public      :: internode_grp      !< group for communicating between compute nodes

    integer          :: nthreads         !< Number of OMP threads
    integer          :: node_type        !< Is this node a P_MASTER or a P_SLAVE?
    logical          :: have_slaves      !< are slaves available?

    integer          :: full_comm        !< The base communicator.
    integer          :: full_comm_rank   !< The rank in the base communicator.
    integer, public  :: master_comm      !< The communicator without slaves.
    integer          :: master_comm_rank !< The rank in the communicator without slaves.
    integer, public  :: slave_intercomm  !< the intercomm to communicate with slaves

    logical          :: reorder_ranks    !< do we reorder ranks in a more compact way?
  end type multicomm_t

  !> An all-pairs communication schedule for a given group.
  type multicomm_all_pairs_t
    private
    type(mpi_grp_t)  :: grp                    !< Schedule for this group.
    integer          :: rounds                 !< This many comm. rounds.
    integer, allocatable, public :: schedule(:, :) !< This is the schedule.
  end type multicomm_all_pairs_t

contains

  ! ---------------------------------------------------------
  subroutine multicomm_all_pairs_copy(apout, apin)
    type(multicomm_all_pairs_t), intent(inout) :: apout
    type(multicomm_all_pairs_t), intent(in)    :: apin

    PUSH_SUB(multicomm_all_pairs_copy)

    call mpi_grp_copy(apout%grp, apin%grp)
    apout%rounds = apin%rounds
    if (allocated(apin%schedule)) then
      SAFE_ALLOCATE(apout%schedule(1:size(apin%schedule, 1), 1:size(apin%schedule, 2)))
      apout%schedule = apin%schedule
    end if

    POP_SUB(multicomm_all_pairs_copy)
  end subroutine multicomm_all_pairs_copy

  ! ---------------------------------------------------------
  !> create index and domain communicators
  subroutine multicomm_init(mc, namespace, base_grp, parallel_mask, default_mask, n_node, index_range, min_range)
    type(multicomm_t), intent(out)   :: mc
    type(namespace_t), intent(in)    :: namespace
    type(mpi_grp_t),   intent(in)    :: base_grp
    integer,           intent(in)    :: parallel_mask
    integer,           intent(in)    :: default_mask
    integer,           intent(in)    :: n_node
    integer(i8),       intent(inout) :: index_range(:)
    integer,           intent(in)    :: min_range(:)

    integer :: ii, num_slaves, slave_level, ipar
    integer :: parse(1:P_STRATEGY_MAX), default(1:P_STRATEGY_MAX)

    PUSH_SUB(multicomm_init)

    mc%n_node  = n_node

    call messages_print_with_emphasis(msg="Parallelization", namespace=namespace)

    !%Variable ReorderRanks
    !%Default no
    !%Type logical
    !%Section Execution::Parallelization
    !%Description
    !% This variable controls whether the ranks are reorganized to have a more
    !% compact distribution with respect to domain parallelization which needs
    !% to communicate most often. Depending on the system, this can improve
    !% communication speeds.
    !%End
    call parse_variable(namespace, 'ReorderRanks', .false., mc%reorder_ranks)

    call messages_obsolete_variable(namespace, 'ParallelizationStrategy')
    call messages_obsolete_variable(namespace, 'ParallelizationGroupRanks')

    do ipar = 1, P_STRATEGY_MAX
      default(ipar) = PAR_NO
      if (bitand(default_mask, ibset(0, ipar - 1)) /= 0) then
        default(ipar) = PAR_AUTO
      end if
    end do

    !%Variable ParDomains
    !%Type integer
    !%Default auto
    !%Section Execution::Parallelization
    !%Description
    !% This variable controls the number of processors used for the
    !% parallelization in domains.
    !% The special value <tt>auto</tt>, the default, lets Octopus
    !% decide how many processors will be assigned for this
    !% strategy. To disable parallelization in domains, you can use
    !% <tt>ParDomains = no</tt> (or set the number of processors to
    !% 1).
    !%
    !% The total number of processors required is the multiplication
    !% of the processors assigned to each parallelization strategy.
    !%Option auto -1
    !% The number of processors is assigned automatically.
    !%Option no 0
    !% This parallelization strategy is not used.
    !%End
    call parse_variable(namespace, 'ParDomains', default(P_STRATEGY_DOMAINS), parse(P_STRATEGY_DOMAINS))

    !%Variable ParStates
    !%Type integer
    !%Section Execution::Parallelization
    !%Description
    !% This variable controls the number of processors used for the
    !% parallelization in states. The special value <tt>auto</tt> lets
    !% Octopus decide how many processors will be assigned for this
    !% strategy. To disable parallelization in states, you can use
    !% <tt>ParStates = no</tt> (or set the number of processors to 1).
    !%
    !% The default value depends on the <tt>CalculationMode</tt>. For
    !% <tt>CalculationMode = td</tt> the default is <tt>auto</tt>, while
    !% for for other modes the default is <tt>no</tt>.
    !%
    !% The total number of processors required is the multiplication
    !% of the processors assigned to each parallelization strategy.
    !%Option auto -1
    !% The number of processors is assigned automatically.
    !%Option no 0
    !% This parallelization strategy is not used.
    !%End
    call parse_variable(namespace, 'ParStates', default(P_STRATEGY_STATES), parse(P_STRATEGY_STATES))

    !%Variable ParKPoints
    !%Type integer
    !%Default auto
    !%Section Execution::Parallelization
    !%Description
    !% This variable controls the number of processors used for the
    !% parallelization in K-Points and/or spin.
    !% The special value <tt>auto</tt> lets Octopus decide how many processors will be
    !% assigned for this strategy. To disable parallelization in
    !% KPoints, you can use <tt>ParKPoints = no</tt> (or set the
    !% number of processors to 1).
    !%
    !% The total number of processors required is the multiplication
    !% of the processors assigned to each parallelization strategy.
    !%Option auto -1
    !% The number of processors is assigned automatically.
    !%Option no 0
    !% This parallelization strategy is not used.
    !%End
    call parse_variable(namespace, 'ParKPoints', default(P_STRATEGY_KPOINTS), parse(P_STRATEGY_KPOINTS))

    !%Variable ParOther
    !%Type integer
    !%Default auto
    !%Section Execution::Parallelization
    !%Description
    !% This variable controls the number of processors used for the
    !% 'other' parallelization mode, that is CalculatioMode
    !% dependent. For <tt>CalculationMode = casida</tt>, it means
    !% parallelization in electron-hole pairs.
    !%
    !% The special value <tt>auto</tt>,
    !% the default, lets Octopus decide how many processors will be
    !% assigned for this strategy. To disable parallelization in
    !% Other, you can use <tt>ParOther = no</tt> (or set the
    !% number of processors to 1).
    !%
    !% The total number of processors required is the multiplication
    !% of the processors assigned to each parallelization strategy.
    !%Option auto -1
    !% The number of processors is assigned automatically.
    !%Option no 0
    !% This parallelization strategy is not used.
    !%End
    call parse_variable(namespace, 'ParOther', default(P_STRATEGY_OTHER), parse(P_STRATEGY_OTHER))

    do ipar = 1, P_STRATEGY_MAX
      if (parse(ipar) == PAR_NO) parse(ipar) = 1
    end do

    call strategy()

    mc%have_slaves = .false.

    if (mc%par_strategy /= P_STRATEGY_SERIAL) then
      SAFE_ALLOCATE(mc%group_sizes(1:P_STRATEGY_MAX))

      mc%group_sizes = 1

      do ipar = 1, P_STRATEGY_MAX
        if (multicomm_strategy_is_parallel(mc, ipar)) then
          mc%group_sizes(ipar) = parse(ipar)
        else if (parse(ipar) /= 1) then
          call messages_write('Ignoring specification for ' // par_types(ipar))
          call messages_new_line()
          call messages_write('This parallelization strategy is not available.')
          call messages_warning()
        end if
      end do

      call assign_nodes()


      !%Variable ParallelizationNumberSlaves
      !%Type integer
      !%Default 0
      !%Section Execution::Parallelization
      !%Description
      !% Slaves are nodes used for task parallelization. The number of
      !% such nodes is given by this variable multiplied by the number
      !% of domains used in domain parallelization.
      !%End
      call parse_variable(namespace, 'ParallelizationNumberSlaves', 0, num_slaves)

      ! the slaves must be defined at a certain parallelization level, for the moment this is state parallelization.
      slave_level = P_STRATEGY_STATES
      mc%have_slaves = (num_slaves > 0)

      if (mc%have_slaves) then
        call messages_experimental('Task parallelization')
      end if

      ! clear parallel strategies that were available but will not be used
      do ii = 1, P_STRATEGY_MAX
        if (mc%group_sizes(ii) == 1) mc%par_strategy = ibclr(mc%par_strategy, ii - 1)
      end do

      ! reset
      call sanity_check()
    end if

    call group_comm_create()

    call messages_print_with_emphasis(namespace=namespace)

    POP_SUB(multicomm_init)

  contains

    ! ---------------------------------------------------------
    subroutine strategy()
      integer :: jj, ipar

      PUSH_SUB(multicomm_init.strategy)

      if (base_grp%size > 1) then

        mc%par_strategy = 0

        do ipar = 1, P_STRATEGY_MAX
          if (parse(ipar) == PAR_AUTO .or. parse(ipar) > 1) then
            mc%par_strategy = ibset(mc%par_strategy, ipar - 1)
          end if
        end do

        if (mc%par_strategy /= bitand(mc%par_strategy, parallel_mask)) then
          call messages_write('Parallelization strategies unavailable for this run mode are being discarded.')
          call messages_warning()
        end if

        mc%par_strategy = bitand(mc%par_strategy, parallel_mask)

        if (mc%par_strategy == P_STRATEGY_SERIAL) then
          message(1) = "More than one node is available, but this run mode cannot run with the requested parallelization."
          message(2) = "Please select a parallelization strategy compatible with"
          jj = 2
          do ii = 1, n_par_types
            if (bitand(parallel_mask, 2**(ii - 1)) /= 0) then
              jj = jj + 1
              write(message(jj), '(2a)') "  -> ", par_types(ii)
            end if
          end do
          jj=jj+1
          write(message(jj),'(a,i6)') "mc%par_strategy is : ",mc%par_strategy
          call messages_fatal(jj, only_root_writes = .true.)
        end if
      else
        mc%par_strategy = P_STRATEGY_SERIAL
      end if

      mc%nthreads = 1
#if defined(HAVE_OPENMP)
      !$omp parallel
      !$omp master
      mc%nthreads = omp_get_num_threads()
      !$omp end master
      !$omp end parallel
#endif

      if (mc%par_strategy == P_STRATEGY_SERIAL .and. mc%nthreads == 1) then
        message(1) = "Info: Octopus will run in *serial*"
        call messages_info(1, namespace=namespace)
      else
        write(message(1),'(a)')     'Info: Octopus will run in *parallel*'
        write(message(2),'(a)')     ''
        write(message(3),'(a, i8)') '      Number of processes           :', base_grp%size
        write(message(4),'(a, i8)') '      Number of threads per process :', mc%nthreads
        write(message(5),'(a)')     ''
        call messages_info(5, namespace=namespace)
      end if

      POP_SUB(multicomm_init.strategy)
    end subroutine strategy

    ! ---------------------------------------------------------

    subroutine assign_nodes()
      integer :: ii, nn, kk, n_divisors, divisors(1:50)
      integer(i8) :: n_group_max(1:P_STRATEGY_MAX)

      PUSH_SUB(multicomm_init.assign_nodes)

      ! this is the maximum number of processors in each group
      n_group_max(1:P_STRATEGY_MAX) = max(index_range(1:P_STRATEGY_MAX), 1)
      do kk = 1, P_STRATEGY_MAX
        if (.not. multicomm_strategy_is_parallel(mc, kk)) n_group_max(kk) = 1
      end do

      if (debug%info) then
        call messages_write('Debug info: Allowable group ranks:', new_line = .true.)
        do kk = 1, P_STRATEGY_MAX
          call messages_write(par_types(kk), fmt = '2x,a12,":",1x')
          call messages_write(n_group_max(kk), new_line = .true.)
        end do
        call messages_info()
      end if

      nn = mc%n_node

      ! first loop, check the processors assigned by the user
      do ipar = P_STRATEGY_MAX, 1, -1

        if (mc%group_sizes(ipar) == PAR_AUTO) cycle

        if (mc%group_sizes(ipar) > n_group_max(ipar)) then
          call messages_write('The number of processors specified for '//par_types(ipar)//'(')
          call messages_write(mc%group_sizes(ipar))
          call messages_write(')', new_line = .true.)
          call messages_write('is larger than the degrees of freedom for that level (')
          call messages_write(n_group_max(ipar))
          call messages_write(').')
          call messages_warning()
        end if

        if (mod(nn, mc%group_sizes(ipar)) /= 0) then
          call messages_write('The number of processors specified for '//par_types(ipar)//'(')
          call messages_write(mc%group_sizes(ipar))
          call messages_write(')', new_line = .true.)
          call messages_write('is not a divisor of the number of processors (')
          call messages_write(mc%n_node)
          call messages_write(').')
          call messages_fatal()
        end if

        nn = nn/mc%group_sizes(ipar)

      end do

      ! second loop, now assign the rest automatically
      do ipar = P_STRATEGY_MAX, 1, -1

        if (mc%group_sizes(ipar) /= PAR_AUTO) cycle

        n_divisors = ubound(divisors, dim = 1)
        call get_divisors(nn, n_divisors, divisors)

        mc%group_sizes(ipar) = nn
        do ii = 2, n_divisors
          if (divisors(ii) > n_group_max(ipar)) then
            mc%group_sizes(ipar) = divisors(ii - 1)
            exit
          end if
        end do

        nn = nn/mc%group_sizes(ipar)

      end do

      POP_SUB(multicomm_init.assign_nodes)
    end subroutine assign_nodes


    ! ---------------------------------------------------------
    !> check if a balanced distribution of nodes will be used
    subroutine sanity_check()
      FLOAT :: frac
      integer :: ii, kk
      integer(i8) :: jj, n_max
      integer :: real_group_sizes(1:MAX_INDEX)

      PUSH_SUB(multicomm_init.sanity_check)

      if (num_slaves > 0) then

        if (mc%group_sizes(slave_level) < num_slaves + 1) then
          message(1) = 'Too many nodes assigned to task parallelization.'
          call messages_fatal(1)
        end if

        write(message(1),'(a,i6)') 'Info: Number of slaves nodes              :', &
          num_slaves*product(mc%group_sizes(1:slave_level - 1))
        call messages_info(1)

      end if

      ! print out some info
      ii = 0
      do kk = P_STRATEGY_MAX, 1, -1
        real_group_sizes(kk) = mc%group_sizes(kk)
        if (.not. multicomm_strategy_is_parallel(mc, kk)) cycle
        ii = ii + 1
        if (kk == slave_level) real_group_sizes(kk) = real_group_sizes(kk) - num_slaves
        write(message(ii),'(3a,i6,a,i12,a)') 'Info: Number of nodes in ', &
          par_types(kk), ' group:', real_group_sizes(kk), ' (', index_range(kk), ')'
      end do
      call messages_info(ii)

      ! do we have the correct number of processors
      if (product(mc%group_sizes(1:P_STRATEGY_MAX)) /= base_grp%size) then
        write(message(1),'(a)') 'Inconsistent number of processors:'
        write(message(2),'(a,i6)') '  MPI processes      = ', base_grp%size
        write(message(3),'(a,i6)') '  Required processes = ', product(mc%group_sizes(1:P_STRATEGY_MAX))
        message(4) = ''
        message(5) = 'You probably have a problem in the ParDomains, ParStates, ParKPoints or ParOther.'
        call messages_fatal(5, only_root_writes = .true.)
      end if

      if (any(real_group_sizes(1:P_STRATEGY_MAX) > index_range(1:P_STRATEGY_MAX))) then
        message(1) = "Could not distribute nodes in parallel job. Most likely you are trying to"
        message(2) = "use too many nodes for the job."
        call messages_fatal(2, only_root_writes = .true.)
      end if

      if (any(index_range(1:P_STRATEGY_MAX) / real_group_sizes(1:P_STRATEGY_MAX) < min_range(1:P_STRATEGY_MAX) .and. &
        real_group_sizes(1:P_STRATEGY_MAX) >  1)) then
        message(1) = "I have fewer elements in a parallel group than recommended."
        message(2) = "Maybe you should reduce the number of nodes."
        call messages_warning(2)
      end if

      ! calculate fraction of idle time
      frac = M_ONE
      do ii = 1, P_STRATEGY_MAX
        n_max = ceiling(TOFLOAT(index_range(ii)) / TOFLOAT(real_group_sizes(ii)))
        jj = n_max*real_group_sizes(ii)
        frac = frac*(M_ONE - TOFLOAT(jj - index_range(ii)) / TOFLOAT(jj))
      end do

      write(message(1), '(a,f5.2,a)') "Info: Octopus will waste at least ", &
        (M_ONE - frac)*CNST(100.0), "% of computer time."
      if (frac < CNST(0.8)) then
        message(2) = "Usually a number of processors which is a multiple of small primes is best."
        call messages_warning(2)
      else
        call messages_info(1)
      end if

      POP_SUB(multicomm_init.sanity_check)
    end subroutine sanity_check

    ! ---------------------------------------------------------
    subroutine group_comm_create()
#if defined(HAVE_MPI)
      logical :: dim_mask(MAX_INDEX)
      integer :: i_strategy, irank
      logical :: reorder, periodic_mask(MAX_INDEX)
      integer :: coords(MAX_INDEX)
      integer :: new_comm, new_comm_size
      character(len=6) :: node_type
      type(mpi_grp_t) :: reorder_grp
      integer :: base_group, reorder_group, ranks(base_grp%size)
      integer :: ii, jj, kk, ll, nn, reorder_comm
#endif

      PUSH_SUB(multicomm_init.group_comm_create)

      mc%node_type = P_MASTER

      SAFE_ALLOCATE(mc%group_comm(1:P_STRATEGY_MAX))
      SAFE_ALLOCATE(mc%who_am_i(1:P_STRATEGY_MAX))

#if defined(HAVE_MPI)
      mc%full_comm = MPI_COMM_NULL
      mc%slave_intercomm = MPI_COMM_NULL
      if (mc%par_strategy /= P_STRATEGY_SERIAL) then
        if (mc%reorder_ranks) then
          ! first, reorder the ranks
          ! this is done to get a column-major ordering of the ranks in the
          ! Cartesian communicator, since they a ordered row-major otherwise
          call MPI_Comm_group(base_grp%comm, base_group, mpi_err)
          if (mpi_err /= MPI_SUCCESS) then
            message(1) = "Error in getting MPI group!"
            call messages_fatal(1)
          end if
          ! now transpose the hypercube => get rank numbers in column-major order
          nn = 1
          do ii = 1, mc%group_sizes(1)
            do jj = 1, mc%group_sizes(2)
              do kk = 1, mc%group_sizes(3)
                do ll = 1, mc%group_sizes(4)
                  ranks(nn) = (ll-1)*mc%group_sizes(3)*mc%group_sizes(2)*mc%group_sizes(1) &
                    + (kk-1)*mc%group_sizes(2)*mc%group_sizes(1) &
                    + (jj-1)*mc%group_sizes(1) + ii - 1
                  nn = nn + 1
                end do
              end do
            end do
          end do
          call MPI_Group_incl(base_group, base_grp%size, ranks, reorder_group, mpi_err)
          if (mpi_err /= MPI_SUCCESS) then
            message(1) = "Error in creating MPI group!"
            call messages_fatal(1)
          end if
          ! now get the reordered communicator
          call MPI_Comm_create(base_grp%comm, reorder_group, reorder_comm, mpi_err)
          if (mpi_err /= MPI_SUCCESS) then
            message(1) = "Error in creating reordered communicator!"
            call messages_fatal(1)
          end if
          call mpi_grp_init(reorder_grp, reorder_comm)
        else
          call mpi_grp_copy(reorder_grp, base_grp)
        end if

        ! Multilevel parallelization is organized in a hypercube. We
        ! use an MPI Cartesian topology to generate the communicators
        ! that correspond to each level.

        ! create the topology
        periodic_mask = .false.
        reorder = .true.

        ! The domain and states dimensions have to be periodic (2D torus)
        ! in order to circulate matrix blocks.
        periodic_mask(P_STRATEGY_DOMAINS) = multicomm_strategy_is_parallel(mc, P_STRATEGY_DOMAINS)
        periodic_mask(P_STRATEGY_STATES)  = multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)

        ! We allow reordering of ranks.
        call MPI_Cart_create(reorder_grp%comm, P_STRATEGY_MAX, mc%group_sizes, periodic_mask, reorder, mc%full_comm, mpi_err)

        call MPI_Comm_rank(mc%full_comm, mc%full_comm_rank, mpi_err)

        ! get the coordinates of the current processor
        call MPI_Cart_coords(mc%full_comm, mc%full_comm_rank, P_STRATEGY_MAX, coords, mpi_err)

        ! find out what type of node this is
        if (coords(slave_level) >= mc%group_sizes(slave_level) - num_slaves) then
          mc%node_type = P_SLAVE
        end if

        if (mc%node_type == P_MASTER) then
          mc%group_sizes(slave_level) = mc%group_sizes(slave_level) - num_slaves
        else
          mc%group_sizes(slave_level) = num_slaves
        end if

        call MPI_Comm_split(mc%full_comm, mc%node_type, mc%full_comm_rank, new_comm, mpi_err)
        ASSERT(new_comm /= MPI_COMM_NULL)
        call MPI_Comm_size(new_comm, new_comm_size, mpi_err)

        reorder = .false.
        if (product(mc%group_sizes(:)) /= new_comm_size) then
          write(stderr,*) 'node ', mpi_world%rank, ': mc%group_sizes = ', mc%group_sizes, ' new_comm_size = ', new_comm_size
          call mpi_world%barrier()
          ASSERT(product(mc%group_sizes(:)) == new_comm_size)
        end if
        call MPI_Cart_create(new_comm, P_STRATEGY_MAX, mc%group_sizes, periodic_mask, reorder, mc%master_comm, mpi_err)
        ASSERT(mc%master_comm /= MPI_COMM_NULL)

        call MPI_Comm_free(new_comm, mpi_err)

        call MPI_Comm_rank(mc%master_comm, mc%master_comm_rank, mpi_err)

        ! The "lines" of the Cartesian grid.
        ! Initialize all the communicators, even if they are not parallelized
        do i_strategy = 1, P_STRATEGY_MAX
          dim_mask             = .false.
          dim_mask(i_strategy) = .true.
          call MPI_Cart_sub(mc%master_comm, dim_mask, mc%group_comm(i_strategy), mpi_err)
          call MPI_Comm_rank(mc%group_comm(i_strategy), mc%who_am_i(i_strategy), mpi_err)
        end do

        ! The domain-state "planes" of the grid (the ones with periodic dimensions).
        dim_mask                     = .false.
        dim_mask(P_STRATEGY_DOMAINS) = .true.
        dim_mask(P_STRATEGY_STATES)  = .true.
        call MPI_Cart_sub(mc%master_comm, dim_mask, mc%dom_st_comm, mpi_err)

        ! The state-kpoints "planes" of the grid
        dim_mask                     = .false.
        dim_mask(P_STRATEGY_STATES)  = .true.
        dim_mask(P_STRATEGY_KPOINTS) = .true.
        call MPI_Cart_sub(mc%master_comm, dim_mask, mc%st_kpt_comm, mpi_err)

        ! The domains-states-kpoints "cubes" of the grid
        dim_mask                     = .false.
        dim_mask(P_STRATEGY_DOMAINS) = .true.
        dim_mask(P_STRATEGY_STATES)  = .true.
        dim_mask(P_STRATEGY_KPOINTS) = .true.
        call MPI_Cart_sub(mc%master_comm, dim_mask, mc%dom_st_kpt_comm, mpi_err)

        if (num_slaves > 0) call create_slave_intercommunicators()

        call create_intranode_communicator(base_grp, mc%intranode_grp, mc%internode_grp)
      else
        ! we initialize these communicators so we can use them even in serial
        mc%group_comm = base_grp%comm
        mc%who_am_i   = 0
        mc%master_comm = base_grp%comm
        mc%dom_st_comm = base_grp%comm
        mc%st_kpt_comm = base_grp%comm
        mc%dom_st_kpt_comm = base_grp%comm
        call mpi_grp_copy(mc%intranode_grp, base_grp)
        call mpi_grp_copy(mc%internode_grp, base_grp)
      end if

      ! This is temporary debugging information.
      if (debug%info .and. mc%par_strategy /= P_STRATEGY_SERIAL) then
        write(message(1),'(a)') 'Debug: MPI Task Assignment to MPI Groups'
        write(message(2),'(5a10)') 'World', 'Domains', 'States', 'K-Points', 'Other'
        call messages_info(1)

        if (mc%node_type == P_SLAVE) then
          node_type = "slave"
        else
          node_type = "master"
        end if
        do irank = 0, mpi_world%size - 1
          if (mpi_world%rank == irank) then
            write(message(1),'(5i10,5x,a)') mpi_world%rank, mc%who_am_i(P_STRATEGY_DOMAINS), mc%who_am_i(P_STRATEGY_STATES), &
              mc%who_am_i(P_STRATEGY_KPOINTS), mc%who_am_i(P_STRATEGY_OTHER), trim(node_type)
            call messages_info(1, all_nodes = .true.)
          end if
          call mpi_world%barrier()
        end do
      end if

#else
      mc%group_comm = -1
      mc%who_am_i   = 0
      mc%full_comm = -1
      mc%master_comm = -1
      mc%dom_st_comm = -1
      mc%st_kpt_comm = -1
      mc%dom_st_kpt_comm = -1
      mc%slave_intercomm = -1
#endif

      POP_SUB(multicomm_init.group_comm_create)
    end subroutine group_comm_create

    ! -----------------------------------------------------
#ifdef HAVE_MPI
    subroutine create_slave_intercommunicators()
      integer :: remote_leader
      integer :: tag
      integer :: coords(MAX_INDEX)

      PUSH_SUB(multicomm_init.create_slave_intercommunicators)

      ! create the intercommunicators to communicate with slaves

      ! get the coordinates of the current processor
      call MPI_Cart_coords(mc%full_comm, mc%full_comm_rank, P_STRATEGY_MAX, coords, mpi_err)

      !now get the rank of the remote_leader
      if (mc%node_type == P_SLAVE) then
        coords(slave_level) = 0
      else
        coords(slave_level) = mc%group_sizes(slave_level)
      end if
      call MPI_Cart_rank(mc%full_comm, coords, remote_leader, mpi_err)

      ! now create the intercommunicator
      tag = coords(P_STRATEGY_DOMAINS)
      call MPI_Intercomm_create(mc%group_comm(slave_level), 0, base_grp%comm, remote_leader, tag, mc%slave_intercomm, mpi_err)

      POP_SUB(multicomm_init.create_slave_intercommunicators)
    end subroutine create_slave_intercommunicators
#endif
  end subroutine multicomm_init


  ! ---------------------------------------------------------
  subroutine multicomm_end(mc)
    type(multicomm_t), intent(inout) :: mc

#if defined(HAVE_MPI)
    integer :: ii
#endif

    PUSH_SUB(multicomm_end)

    if (mc%par_strategy /= P_STRATEGY_SERIAL) then
#if defined(HAVE_MPI)
      ! Delete communicators.
      do ii = 1, P_STRATEGY_MAX
        ! initialized even if not parallelized
        call MPI_Comm_free(mc%group_comm(ii), mpi_err)
      end do
      call MPI_Comm_free(mc%dom_st_comm, mpi_err)
      call MPI_Comm_free(mc%st_kpt_comm, mpi_err)
      call MPI_Comm_free(mc%dom_st_kpt_comm, mpi_err)
      call MPI_Comm_free(mc%full_comm, mpi_err)
      call MPI_Comm_free(mc%master_comm, mpi_err)

      if (multicomm_have_slaves(mc)) call MPI_Comm_free(mc%slave_intercomm, mpi_err)

#endif
    end if

    SAFE_DEALLOCATE_A(mc%group_sizes)
    SAFE_DEALLOCATE_A(mc%group_comm)
    SAFE_DEALLOCATE_A(mc%who_am_i)

    POP_SUB(multicomm_end)
  end subroutine multicomm_end


  ! ---------------------------------------------------------
  logical pure function multicomm_strategy_is_parallel(mc, level) result(rr)
    type(multicomm_t), intent(in) :: mc
    integer,           intent(in) :: level

    rr = bitand(mc%par_strategy, 2**(level - 1)) /= 0

  end function multicomm_strategy_is_parallel


  ! ---------------------------------------------------------
  !> This routine uses the one-factorization (or near-one-factorization
  !! of a complete graph to construct an all-pair communication
  !! schedule (cf. Wang, X., Blum, E. K., Parker, D. S., and Massey,
  !! D. 1997. The dance-party problem and its application to collective
  !! communication in computer networks. Parallel Comput. 23, 8
  !! (Aug. 1997), 1141-1156.
#if defined(HAVE_MPI)

  subroutine multicomm_create_all_pairs(mpi_grp, ap)
    type(mpi_grp_t),             intent(in)  :: mpi_grp
    type(multicomm_all_pairs_t), intent(out) :: ap

    integer :: grp_size, rounds, ir, in

    PUSH_SUB(create_all_pairs)

    ap%grp = mpi_grp
    grp_size = mpi_grp%size

    ! Number of rounds.
    if (mod(grp_size, 2) == 0) then
      rounds = grp_size - 1
    else
      rounds = grp_size
    end if
    ap%rounds = rounds

    ! Calculate schedule.
    SAFE_ALLOCATE(ap%schedule(0:grp_size - 1, 1:rounds))
    do ir = 1, rounds
      do in = 0, grp_size - 1
        ap%schedule(in, ir) = get_partner(in + 1, ir) - 1
      end do
    end do

    POP_SUB(create_all_pairs)

  contains

    ! ---------------------------------------------------------
    !> Those are from the paper cited above.
    integer pure function get_partner(in, ir)
      integer, intent(in) :: in, ir

      ! No PUSH SUB, called too often.

      if (mod(grp_size, 2) == 0) then
        get_partner = get_partner_even(grp_size, in - 1, ir - 1) + 1
      else
        get_partner = get_partner_odd(grp_size, in - 1, ir - 1) + 1
      end if

    end function get_partner

    ! ---------------------------------------------------------
    integer pure function get_partner_even(grp_size, ii, rr) result(pp)
      integer, intent(in) :: grp_size, ii, rr

      integer :: mm

      ! No PUSH SUB, called too often.

      mm = grp_size / 2

      if (ii == 0) then
        pp = rr + 1
      elseif (ii == rr + 1) then
        pp = 0
      else
        ! I never know when to use which remainder function, but here
        ! it has to be the modulo one. Do not change that!
        pp = modulo(2 * rr - ii + 1, 2 * mm - 1) + 1
      end if

    end function get_partner_even

    ! ---------------------------------------------------------
    integer pure function get_partner_odd(grp_size, ii, rr) result(pp)
      integer, intent(in) :: grp_size, ii, rr

      integer :: mm

      ! No PUSH SUB, called too often.

      mm = (grp_size + 1) / 2

      pp = get_partner_even(grp_size + 1, ii, rr)

      if (pp == 2 * mm - 1) then
        pp = ii
      end if

    end function get_partner_odd

  end subroutine multicomm_create_all_pairs
#endif

  !---------------------------------------------------
  !> Function to divide the range of numbers from 1 to nobjs
  !! between nprocs processors.
  !! THREADSAFE
  subroutine multicomm_divide_range(nobjs, nprocs, istart, ifinal, lsize, scalapack_compat)
    integer,           intent(in)  :: nobjs !< Number of points to divide
    integer,           intent(in)  :: nprocs !< Number of processors
    integer,           intent(out) :: istart(:)
    integer,           intent(out) :: ifinal(:)
    integer, optional, intent(out) :: lsize(:) !< Number of objects in each partition
    logical, optional, intent(in)  :: scalapack_compat

    integer :: ii, jj, rank
    logical :: scalapack_compat_
    integer :: nbl, size

    ! no push_sub, threadsafe

    scalapack_compat_ = optional_default(scalapack_compat, .false.)
#ifndef HAVE_SCALAPACK
    scalapack_compat_ = .false.
#endif

    if (scalapack_compat_) then
      nbl = nobjs/nprocs
      if (mod(nobjs, nprocs) /= 0) nbl = nbl + 1

      istart(1) = 1
      do rank = 1, nprocs
#ifdef HAVE_SCALAPACK
        size = numroc(nobjs, nbl, rank - 1, 0, nprocs)
#endif
        if (size > 0) then
          if (rank > 1) istart(rank) = ifinal(rank - 1) + 1
          ifinal(rank) = istart(rank) + size - 1
        else
          istart(rank) = 1
          ifinal(rank) = 0
        end if
      end do
    else

      if (nprocs <= nobjs) then

        ! procs are assigned to groups by round robin
        do rank = 0, nprocs - 1
          jj = nobjs / nprocs
          ii = nobjs - jj*nprocs
          if (ii > 0 .and. rank < ii) then
            jj = jj + 1
            istart(rank + 1) = rank*jj + 1
            ifinal(rank + 1) = istart(rank + 1) + jj - 1
          else
            ifinal(rank + 1) = nobjs - (nprocs - rank - 1)*jj
            istart(rank + 1) = ifinal(rank + 1) - jj + 1
          end if
        end do

      else
        do ii = 1, nprocs
          if (ii <= nobjs) then
            istart(ii) = ii
            ifinal(ii) = ii
          else
            istart(ii) = 1
            ifinal(ii) = 0
          end if
        end do
      end if
    end if

    if (present(lsize)) then
      lsize(1:nprocs) = ifinal(1:nprocs) - istart(1:nprocs) + 1
      ASSERT(sum(lsize(1:nprocs)) == nobjs)
    end if

  end subroutine multicomm_divide_range

  ! ---------------------------------------------------------
  !> Function to divide the range of numbers from 1 to nobjs
  !! between all available threads with OpenMP.
  ! THREADSAFE
  subroutine multicomm_divide_range_omp(nobjs, ini, nobjs_loc)
    integer, intent(in)    :: nobjs       !< Number of points to divide
    integer, intent(out)   :: ini         !< Start point of the partition
    integer, intent(out)   :: nobjs_loc   !< Number of objects in each partition

    integer :: rank
#ifdef HAVE_OPENMP
    integer, allocatable :: istart(:), ifinal(:), lsize(:)
    integer :: nthreads
#endif

    ! no push_sub, threadsafe
    rank = 1
#ifdef HAVE_OPENMP
    nthreads = omp_get_num_threads()
    SAFE_ALLOCATE(istart(1:nthreads))
    SAFE_ALLOCATE(ifinal(1:nthreads))
    SAFE_ALLOCATE(lsize(1:nthreads))
    call multicomm_divide_range(nobjs, nthreads, istart, ifinal, lsize)
    rank   = 1 + omp_get_thread_num()
    ini    = istart(rank)
    nobjs_loc = lsize(rank)
    SAFE_DEALLOCATE_A(istart)
    SAFE_DEALLOCATE_A(ifinal)
    SAFE_DEALLOCATE_A(lsize)
#else
    ini = 1
    nobjs_loc = nobjs
#endif

  end subroutine multicomm_divide_range_omp

  ! ---------------------------------------------------------

  logical pure function multicomm_is_slave(this) result(slave)
    type(multicomm_t), intent(in) :: this

    slave = this%node_type == P_SLAVE
  end function multicomm_is_slave

  ! ---------------------------------------------------------

  logical pure function multicomm_have_slaves(this) result(have_slaves)
    type(multicomm_t), intent(in) :: this

    have_slaves = this%have_slaves
  end function multicomm_have_slaves

end module multicomm_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
