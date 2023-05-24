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

module io_oct_m
  use debug_oct_m
  use global_oct_m
  use loct_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m

  implicit none

  private
  public ::              &
    io_workpath,         &
    io_open,             &
    io_mkdir,            &
    io_rm,               &
    io_init,             &
    io_end,              &
    io_status,           &
    io_dump_file,        &
    io_free,             &
    io_close,            &
    io_assign,           &
    io_get_extension,    &
    io_debug_on_the_fly, &
    iopar_read,          &
    iopar_backspace,     &
    iopar_find_line,     &
    io_skip_header,      &
    io_file_exists,      &
    io_dir_exists,       &
    io_get_open_count,   &
    io_get_close_count,  &
    io_incr_open_count,  &
    io_incr_close_count, &
    io_incr_counters

  integer, parameter, public :: iunit_out = 8
  integer, parameter, public :: iunit_err = 9
  integer, parameter :: min_lun=10, max_lun=99
  logical            :: lun_is_free(min_lun:max_lun)
  character(len=MAX_PATH_LEN) :: work_dir    !< name of the output directory
  integer(i8), save   :: io_open_count
  integer(i8), save   :: io_close_count

contains

  ! ---------------------------------------------------------
  !> If the argument defaults is present and set to true, then the routine
  !! will not try to read anything from the inp file, but set everything
  !! to the default values.
  subroutine io_init(defaults)
    logical, optional, intent(in)    :: defaults

    character(len=MAX_PATH_LEN) :: filename

    io_open_count = 0
    io_close_count = 0

    ! cannot use push/pop before initializing io

    if (present(defaults)) then
      if (defaults) then
        lun_is_free(min_lun:max_lun)=.true.
        stdin  = 5
        stdout = 6
        stderr = 0
        work_dir = '.'
        return
      end if
    end if

    lun_is_free(min_lun:max_lun)=.true.
    stdin = 5

    !%Variable stdout
    !%Type string
    !%Default "-"
    !%Section Execution::IO
    !%Description
    !% The standard output by default goes to, well, to standard output. This can
    !% be changed by setting this variable: if you give it a name (other than "-")
    !% the output stream is printed in that file instead.
    !%End
    call parse_variable(global_namespace, 'stdout', '-', filename)
    stdout = 6
    if (trim(filename) /= '-') then
      close(stdout)
      open(stdout, file=filename, status='unknown')
    end if

    !%Variable stderr
    !%Type string
    !%Default "-"
    !%Section Execution::IO
    !%Description
    !% The standard error by default goes to, well, to standard error. This can
    !% be changed by setting this variable: if you give it a name (other than "-")
    !% the output stream is printed in that file instead.
    !%End
    call parse_variable(global_namespace, 'stderr', '-', filename)
    stderr = 0
    if (trim(filename) /= '-') then
      close(stderr)
      open(stderr, file=filename, status='unknown')
    end if

    !%Variable WorkDir
    !%Type string
    !%Default "."
    !%Section Execution::IO
    !%Description
    !% By default, all files are written and read from the working directory,
    !% <i>i.e.</i> the directory from which the executable was launched. This behavior can
    !% be changed by setting this variable. If you set <tt>WorkDir</tt> to a name other than ".",
    !% the following directories are written and read in that directory:
    !%<ul>
    !% <li>"casida/"</li>
    !% <li>"em_resp_fd/"</li>
    !% <li>"em_resp/"</li>
    !% <li>"geom/"</li>
    !% <li>"kdotp/"</li>
    !% <li>"local.general"</li>
    !% <li>"pcm/"</li>
    !% <li>"profiling/"</li>
    !% <li>"restart/"</li>
    !% <li>"static/"</li>
    !% <li>"td.general/"</li>
    !% <li>"vdw/"</li>
    !% <li>"vib_modes/"</li>
    !%</ul>
    !% Furthermore, some of the debug information (see <tt>Debug</tt>) is also written to <tt>WorkDir</tt> and
    !% the non-absolute paths defined in <tt>OutputIterDir</tt> are relative to <tt>WorkDir</tt>.
    !%End
    call parse_variable(global_namespace, 'WorkDir', '.', work_dir)
    ! ... and if necessary create workdir (will not harm if work_dir is already there)
    if (work_dir /= '.') call loct_mkdir(trim(work_dir))

    if (debug%info .or. debug%interaction_graph .or. debug%propagation_graph) then
      call io_mkdir('debug', global_namespace)
    end if

    if (debug%trace_file) then
      !wipe out debug trace files from previous runs to start fresh rather than appending
      call debug_delete_trace()
    end if

  end subroutine io_init

  ! ---------------------------------------------------------
  subroutine io_end()

    ! no PUSH/POP, because the POP would write to stderr after it was closed.

    if (stderr /= 0) call io_close(stderr)
    if (stdin  /= 5) call io_close(stdin)
    if (stdout /= 6) call io_close(stdout)

  end subroutine io_end


  ! ---------------------------------------------------------
  subroutine io_assign(got_lun)
    integer, intent(out) :: got_lun

    integer :: iostat, lun
    logical :: used

    got_lun = -1

    ! Looks for a free unit and assigns it to lun
    do lun = min_lun, max_lun
      if (lun_is_free(lun)) then
        inquire(unit=lun, opened=used, iostat=iostat)

        if (iostat /= 0) used = .true.
        lun_is_free(lun) = .false.
        if (.not. used) then
          got_lun = lun
          exit
        end if
      end if
    end do

  end subroutine io_assign


  ! ---------------------------------------------------------
  subroutine io_free(lun)
    integer, intent(in) :: lun

    if (lun >= min_lun .and. lun  <=  max_lun) &
      lun_is_free(lun) = .true.

  end subroutine io_free


  ! ---------------------------------------------------------
  character(len=MAX_PATH_LEN) function io_workpath(path, namespace) result(wpath)
    character(len=*),            intent(in) :: path
    type(namespace_t), optional, intent(in) :: namespace

    logical :: absolute_path
    integer :: total_len

    ! use the logical to avoid problems with the string length
    absolute_path = .false.
    if (len_trim(path) > 0) then
      absolute_path = path(1:1) == '/'
    end if

    ! check that the path is not longer than the maximum allowed
    total_len = len_trim(path)
    if (.not. absolute_path) then
      total_len = total_len + len_trim(work_dir) + 1
      if (present(namespace)) then
        if (namespace%len() > 0) total_len = total_len + namespace%len() + 1
      end if
    end if
    if (total_len > MAX_PATH_LEN) then
      write(stderr,"(A,I5)") "Path is longer than the maximum path length of ", MAX_PATH_LEN
    end if

    if (absolute_path) then
      ! we do not change absolute path names
      wpath = trim(path)
    else
      wpath = trim(work_dir)
      if (present(namespace)) then
        ! insert namespace into path
        if (namespace%len() > 0) wpath = trim(wpath) + "/" + trim(namespace%get('/'))
      end if
      wpath = trim(wpath) + "/" + trim(path)
    end if

  end function io_workpath


  ! ---------------------------------------------------------
  subroutine io_mkdir(fname, namespace, parents)
    character(len=*),            intent(in) :: fname
    type(namespace_t), optional, intent(in) :: namespace
    logical,           optional, intent(in) :: parents

    logical :: parents_
    integer :: last_slash, pos, length

    parents_ = .false.
    if (present(parents)) parents_ = parents

    if (.not. parents_) then
      call loct_mkdir(trim(io_workpath("", namespace=namespace)))
      call loct_mkdir(trim(io_workpath(fname, namespace=namespace)))
    else
      last_slash = max(index(fname, "/", .true.), len_trim(fname))
      pos = 1
      length = index(fname, '/') - 1
      do while (pos < last_slash)
        call loct_mkdir(trim(io_workpath(fname(1:pos+length-1), namespace=namespace)))
        pos = pos + length + 1
        length = index(fname(pos:), "/") - 1
        if (length < 1) length = len_trim(fname(pos:))
      end do

    end if

  end subroutine io_mkdir


  ! ---------------------------------------------------------
  subroutine io_rm(fname, namespace)
    character(len=*),            intent(in) :: fname
    type(namespace_t), optional, intent(in) :: namespace

    call loct_rm(trim(io_workpath(fname, namespace=namespace)))

  end subroutine io_rm


  ! ---------------------------------------------------------
  integer function io_open(file, namespace, action, status, form, position, die, recl, grp) result(iunit)
    character(len=*), intent(in)           :: file, action
    type(namespace_t),intent(in), optional :: namespace
    character(len=*), intent(in), optional :: status, form, position
    logical,          intent(in), optional :: die
    integer,          intent(in), optional :: recl
    type(mpi_grp_t),  intent(in), optional :: grp

    character(len=20)  :: status_, form_, position_
    character(len=MAX_PATH_LEN) :: file_
    logical            :: die_
    integer            :: iostat
    character(len=100) :: io_emsg
    type(mpi_grp_t)    :: grp_

    if (present(grp)) then
      grp_%comm = grp%comm
      grp_%rank = grp%rank
      grp_%size = grp%size
    else
      call mpi_grp_init(grp_, -1)
    end if


    if (mpi_grp_is_root(grp_)) then

      status_ = 'unknown'
      if (present(status)) status_ = status
      form_ = 'formatted'
      if (present(form)) form_  = form
      position_ = 'asis'
      if (present(position)) position_ = position
      die_ = .true.
      if (present(die)) die_ = die

      call io_assign(iunit)
      if (iunit < 0) then
        if (die_) then
          write(stderr, '(a)') '*** IO Error: Too many files open.'
        end if
        return
      end if

      file_ = io_workpath(file, namespace=namespace)

      if (present(recl)) then
        open(unit=iunit, file=trim(file_), status=trim(status_), form=trim(form_), &
          recl=recl, action=trim(action), position=trim(position_), iostat=iostat, iomsg=io_emsg)
      else
        open(unit=iunit, file=trim(file_), status=trim(status_), form=trim(form_), &
          action=trim(action), position=trim(position_), iostat=iostat, iomsg=io_emsg)
      end if

      io_open_count = io_open_count + 1

      if (iostat /= 0) then
        call io_free(iunit)
        iunit = -1
        if (die_) then
          write(stderr, '(a,a)') '*** IO Error: ', trim(io_emsg)
        end if
      end if

    end if

    if (grp_%size > 1) then
      call grp_%bcast(iunit, 1, MPI_INTEGER, 0)
    end if

  end function io_open


  ! ---------------------------------------------------------
  subroutine io_close(iunit, grp)
    integer, intent(inout) :: iunit
    type(mpi_grp_t),  intent(in), optional :: grp

    type(mpi_grp_t)    :: grp_

    if (present(grp)) then
      grp_%comm = grp%comm
      grp_%rank = grp%rank
      grp_%size = grp%size
    else
      call mpi_grp_init(grp_, -1)
    end if

    if (mpi_grp_is_root(grp_)) then
      close(iunit)
      io_close_count = io_close_count + 1
      call io_free(iunit)
    end if

    iunit = -1

  end subroutine io_close


  ! ---------------------------------------------------------
  !> Prints a list of the connected logical units and the names of
  !! the associated files
  ! ---------------------------------------------------------
  subroutine io_status(iunit)
    integer, intent(in) :: iunit

    integer :: ii, iostat
    logical :: opened, named
    character(len=MAX_PATH_LEN) :: filename
    character(len=11) :: form

    write(iunit, '(a)') '******** io_status ********'
    do ii = 0, max_lun
      inquire(ii, opened=opened, named=named, name=filename, form=form, iostat=iostat)
      if (iostat == 0) then
        if (opened) then
          if (.not. named) filename = 'No name available'
          write(iunit, '(i4,5x,a,5x,a)') ii, form, filename
        end if
      else
        write(iunit, '(i4,5x,a)') ii, 'Iostat error'
      end if
    end do
    write(iunit,'(a)') '********           ********'

  end subroutine io_status


  ! ---------------------------------------------------------
  subroutine io_dump_file(ounit, filename)
    integer,          intent(in) :: ounit
    character(len=*), intent(in) :: filename

    integer :: iunit, err
    character(len=80) :: line

    if (.not. mpi_grp_is_root(mpi_world)) return

    call io_assign(iunit)
    open(unit=iunit, file=filename, iostat=err, action='read', status='old')

    do while(err == 0)
      read(iunit, fmt='(a80)', iostat=err) line
      if (err == 0) then
        write(ounit, '(a)') trim(line)
      end if
    end do

    call io_close(iunit)

  end subroutine io_dump_file


  ! ---------------------------------------------------------
  !> Given a path, it returns the extension (if it exists) of the file
  !! (that is, the part of the name that comes after its last point).
  !! If the filename does not have an extension, it returns the empty string.
  character(len=8) function io_get_extension(path) result(ext)
    character(len = *), intent(in)  :: path
    integer :: i, j

    i = index(path, ".", back = .true.)
    j = index(path(i+1:), "/")
    if (i == 0 .or. j /= 0) then
      ext = ""
    else
      ext = path(i+1:)
    end if

  end function io_get_extension


  ! ---------------------------------------------------------
  !> check if debug mode should be enabled or disabled on the fly
  subroutine io_debug_on_the_fly(namespace)
    type(namespace_t), intent(in) :: namespace

    ! only root node performs the check
    if (mpi_grp_is_root(mpi_world)) then
      if (io_file_exists('enable_debug_mode', msg='Enabling DebugMode')) then
        call debug_enable(debug)
        ! this call does not hurt if the directory is already there
        ! but is otherwise required
        call io_mkdir('debug', namespace)
        ! we have been notified by the user, so we can cleanup the file
        call loct_rm('enable_debug_mode')
        ! artificially increase sub stack to avoid underflow
        no_sub_stack = no_sub_stack + 8
      end if

      if (io_file_exists('disable_debug_mode', msg='Disabling DebugMode')) then
        call debug_disable(debug)
        ! we have been notified by the user, so we can cleanup the file
        call loct_rm('disable_debug_mode')
      end if

    end if

  end subroutine io_debug_on_the_fly


  !> Returns true if a file with name 'filename' exists
  !! and issues a reminder.
  ! ---------------------------------------------------------
  logical function io_file_exists(filename, msg) result(file_exists)
    character(len=*),           intent(in)  :: filename
    character(len=*), optional, intent(in)  :: msg

    file_exists = .false.
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists .and. present(msg)) then
      write(stderr, '(a)') trim(msg)
    end if

  end function io_file_exists

  !> Returns true if a dir with name 'dir' exists
  ! ---------------------------------------------------------
  logical function io_dir_exists(dir, namespace)
    character(len=*), intent(in)  :: dir
    type(namespace_t),   intent(in) :: namespace

    io_dir_exists = loct_dir_exists(trim(io_workpath(dir, namespace)))

  end function io_dir_exists

  ! ---------------------------------------------------------
  subroutine iopar_read(grp, iunit, lines, n_lines, ierr)
    type(mpi_grp_t),  intent(in)  :: grp
    integer,          intent(in)  :: iunit
    character(len=*), intent(out) :: lines(:)
    integer,          intent(in)  :: n_lines
    integer,          intent(out) :: ierr

    integer :: il

    ASSERT(n_lines <= size(lines))

    if (mpi_grp_is_root(grp)) then
      do il = 1, n_lines
        read(iunit, '(a)', iostat=ierr) lines(il)
        if (ierr /= 0) exit
      end do
    end if

    call grp%bcast(ierr, 1, MPI_INTEGER, 0)
    call grp%bcast(lines(1), len(lines(1))*n_lines, MPI_CHARACTER, 0)

  end subroutine iopar_read

  ! ---------------------------------------------------------
  subroutine iopar_backspace(grp, iunit)
    type(mpi_grp_t), intent(in)  :: grp
    integer,         intent(in)  :: iunit

    if (mpi_grp_is_root(grp)) then
      backspace(iunit)
    end if

  end subroutine iopar_backspace


  ! ---------------------------------------------------------
  subroutine iopar_find_line(grp, iunit, line, ierr)
    type(mpi_grp_t),  intent(in)  :: grp
    integer,          intent(in)  :: iunit
    character(len=*), intent(in)  :: line
    integer,          intent(out) :: ierr

    character(len=80) :: read_line

    if (mpi_grp_is_root(grp)) then
      rewind(iunit)
      do
        read(iunit, '(a)', iostat=ierr) read_line
        if (ierr /= 0 .or. trim(line) == trim(read_line)) exit
      end do
    end if

    call grp%bcast(ierr, 1, MPI_INTEGER, 0)

  end subroutine iopar_find_line


  ! ---------------------------------------------------------
  subroutine io_skip_header(iunit)
    integer, intent(in) :: iunit

    character(len=1) :: a

    rewind(iunit)
    read(iunit,'(a)') a
    do while(a == '#')
      read(iunit,'(a)') a
    end do
    backspace(iunit)

  end subroutine io_skip_header

  ! ---------------------------------------------------------
  integer(i8) pure function io_get_open_count() result(count)

    count = io_open_count

  end function io_get_open_count

  ! ---------------------------------------------------------
  integer(i8) pure function io_get_close_count() result(count)

    count = io_close_count

  end function io_get_close_count

  ! ---------------------------------------------------------
  subroutine io_incr_open_count()

    io_open_count = io_open_count + 1

  end subroutine io_incr_open_count

  ! ---------------------------------------------------------
  subroutine io_incr_close_count()

    io_close_count = io_close_count + 1

  end subroutine io_incr_close_count

  ! ---------------------------------------------------------
  subroutine io_incr_counters(iio)
    integer, intent(in) :: iio

    integer :: open_count

    open_count = int(iio/100)
    io_open_count = io_open_count + open_count
    io_close_count = io_close_count + iio - open_count * 100

  end subroutine io_incr_counters


end module io_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
