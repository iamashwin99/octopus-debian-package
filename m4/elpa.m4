## Copyright (C) 2016 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
##

AC_DEFUN([ACX_ELPA],
[

if test x"$acx_scalapack_ok" = xyes; then
  acx_elpa_ok=no
  test_elpa_with_openmp=False
  test_elpa_without_openmp=False

  dnl BACKUP LIBS AND FCFLAGS
  acx_elpa_save_LIBS="$LIBS"
  acx_elpa_save_FCFLAGS="$FCFLAGS"

  dnl Check if the library was given in the command line
  AC_ARG_WITH(elpa-prefix, [AS_HELP_STRING([--with-elpa-prefix=DIR], [Directory where elpa was installed.])])

  FCFLAGS_ELPA="${FCFLAGS_ELPA} ${ELPA_FCFLAGS}"

  # Set FCFLAGS_ELPA only if not set from environment
  if [ "$FCFLAGS_ELPA" = " " ] ; then
    case $with_elpa_prefix in
      "") FCFLAGS_ELPA="-I/usr/include -I/usr/include/elpa/modules" ;;
      *)  FCFLAGS_ELPA="-I$with_elpa_prefix/include" ;;
    esac
  fi

  AC_MSG_CHECKING([for elpa])


  elpa_program="AC_LANG_PROGRAM([],[
    use :: elpa
    implicit none
    class(elpa_t), pointer :: e
    integer :: status

    status = elpa_init(20170403)
    e => elpa_allocate()

  ])"

  FCFLAGS="$FCFLAGS_ELPA $acx_elpa_save_FCFLAGS"

  # If octopus is to be compiled without enable_openmp,
  # then skip check with openMP and test without openMP only
  # Here is the table for the different cases
  # |                | LIBS_ELPA="-lelpa_openmp"      | LIBS_ELPA="-lelpa"   | LIBS_ELPA not set                     |
  # | -------------- | ------------------------------ | -------------------- | ------------------------------------- |
  # | Octopus+openmp | Test Elpa_openmp               | Test Elpa sequential | Test Elpa_openmp then Elpa sequential |
  # | Octopus Serial | Fail with incompatible message | Test Elpa sequential | Test Elpa sequential                  |


  if test -z "$enable_openmp" -o x"$enable_openmp" = x"no"; then
    if test "$LIBS_ELPA" = "-lelpa_openmp"; then
      AC_MSG_WARN([Octopus is compiled without openMP but elpa_openmp library is requested, skipping elpa inorder to avoid bringing in openMP])
    elif test -z "$LIBS_ELPA" -o "$LIBS_ELPA" = "-lelpa"; then
      # LIBS_ELPA is not set/ set to serial elpa, check only for elpa without openmp
      test_elpa_without_openmp=True
    else
      AC_MSG_WARN([Unkonw value for LIBS_ELPA: $LIBS_ELPA])
    fi
  else
    if test -z "$LIBS_ELPA"; then
      # LIBS_ELPA is not set, check for elpa with openmp first and then without openmp
      test_elpa_with_openmp=True
      test_elpa_without_openmp=True
    elif test "$LIBS_ELPA" = "-lelpa_openmp"; then
      # LIBS_ELPA is set to lelpa_openmp
      test_elpa_with_openmp=True
    elif test "$LIBS_ELPA" = "-lelpa"; then
      # LIBS_ELPA is set to lelpa
      test_elpa_without_openmp=True
    else
      AC_MSG_WARN([Unkonw value for LIBS_ELPA: $LIBS_ELPA])
    fi
  fi

  # Section to test elpa with openMP support
  if test x$test_elpa_with_openmp = xTrue; then
    AC_MSG_CHECKING([for elpa with openMP support])
    # Set LIBS_ELPA 
    if test ! -z "$with_elpa_prefix"; then
      LIBS_ELPA="-L$with_elpa_prefix/lib -lelpa_openmp"
    else
      LIBS_ELPA="-lelpa_openmp"
    fi

    LIBS="$LIBS_ELPA $acx_elpa_save_LIBS $LIBS_LAPACK $LIBS_BLAS"
    AC_LINK_IFELSE($elpa_program, [acx_elpa_ok=yes], [acx_elpa_ok=no])

    AC_MSG_RESULT([$acx_elpa_ok ($FCFLAGS_ELPA $LIBS_ELPA)])

    if test x$acx_elpa_ok != xyes; then

      AC_MSG_WARN([Could not find the elpa library with openMP support])

      if test x"$test_elpa_without_openmp" = xTrue; then
        AC_MSG_WARN([Trying to find elpa library without openMP support])
      else
        FCFLAGS_ELPA=""
        LIBS_ELPA=""
      fi

    else

      AC_DEFINE(HAVE_ELPA, 1, [Define if ELPA is available])
      # Elpa with openmp is found, no need to test elpa without openmp
      test_elpa_without_openmp=False

    fi
  fi

  # Section to test elpa without openMP support
  if test x$test_elpa_without_openmp = xTrue; then

    AC_MSG_CHECKING([for elpa without openMP support])
      # Set LIBS_ELPA 
      if test ! -z "$with_elpa_prefix"; then
        LIBS_ELPA="-L$with_elpa_prefix/lib -lelpa"
      else
        LIBS_ELPA="-lelpa"
      fi

    LIBS="$LIBS_ELPA $acx_elpa_save_LIBS $LIBS_LAPACK $LIBS_BLAS"
    AC_LINK_IFELSE($elpa_program, [acx_elpa_ok=yes], [acx_elpa_ok=no])

    AC_MSG_RESULT([$acx_elpa_ok ($FCFLAGS_ELPA $LIBS_ELPA)])

    if test x$acx_elpa_ok = xyes; then
      AC_DEFINE(HAVE_ELPA, 1, [Define if ELPA is available with openMP support])
    else
      AC_MSG_WARN([Could not find the elpa library, compiling without elpa support])
      FCFLAGS_ELPA=""
      LIBS_ELPA=""
    fi

  fi

  AC_SUBST(FCFLAGS_ELPA)
  AC_SUBST(LIBS_ELPA)

  FCFLAGS="$acx_elpa_save_FCFLAGS"
  LIBS="$acx_elpa_save_LIBS"
else
  dnl Scalapack is needed to use elpa
  AC_MSG_WARN([elpa library is not used because scalapack is not found which is a prerequisite for elpa])
fi

])
