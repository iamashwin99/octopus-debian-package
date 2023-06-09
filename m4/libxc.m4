## Copyright (C) 2010-2015 M. Marques, X. Andrade, D. Strubbe, M. Oliveira
## Copyright (C) 2020-2021 M. Oliveira
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

AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no
acx_libxc_v5=no
acx_libxc_kxc=yes
acx_libxc_fxc=yes

dnl Check if the library was given in the command line
dnl if not, use environment variables or defaults
AC_ARG_WITH(libxc-prefix, [AS_HELP_STRING([--with-libxc-prefix=DIR], [Directory where libxc was installed.])])

# Set FCFLAGS_LIBXC only if not set from environment
if test x"$FCFLAGS_LIBXC" = x; then
  case $with_libxc_prefix in
    "") FCFLAGS_LIBXC="$ax_cv_f90_modflag/usr/include" ;;
    *)  FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_prefix/include" ;;
  esac
fi

AC_ARG_WITH(libxc-include, [AS_HELP_STRING([--with-libxc-include=DIR], [Directory where libxc Fortran headers were installed.])])
case $with_libxc_include in
  "") ;;
  *)  FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_libxc_save_LIBS="$LIBS"
acx_libxc_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for libxc])

testprog4="AC_LANG_PROGRAM([],[
  use xc_f03_lib_m
  implicit none
  integer :: major
  integer :: minor
  integer :: micro
  integer :: flags = XC_FLAGS_NEEDS_LAPLACIAN
  call xc_f03_version(major, minor, micro)]
  write(*,*) flags)"

dnl This will test to see if the libxc (version 5.1 or greater) has been compiled with the fxc support
testprog5_fxc="AC_LANG_PROGRAM([],[
  use xc_f03_lib_m
  implicit none
  integer :: xc_flags
  type(xc_f03_func_t) :: xc_func
  type(xc_f03_func_info_t) :: xc_info

  call xc_f03_func_init(xc_func, XC_GGA_X_B88, XC_UNPOLARIZED)
  xc_info = xc_f03_func_get_info(xc_func)

  xc_flags = xc_f03_func_info_get_flags(xc_info)
  if(IAND(xc_flags, XC_FLAGS_HAVE_FXC) == 0) then
    stop 1
  endif
  call xc_f03_func_end(xc_func)
  ])"

dnl This will test to see if the libxc (version 5.1 or greater) has been compiled with the kxc support
testprog5_kxc="AC_LANG_PROGRAM([],[
  use xc_f03_lib_m
  implicit none
  integer :: xc_flags
  type(xc_f03_func_t) :: xc_func
  type(xc_f03_func_info_t) :: xc_info

  call xc_f03_func_init(xc_func, XC_GGA_X_B88, XC_UNPOLARIZED)
  xc_info = xc_f03_func_get_info(xc_func)

  xc_flags = xc_f03_func_info_get_flags(xc_info)
  if(IAND(xc_flags, XC_FLAGS_HAVE_KXC) == 0) then
    stop 1
  endif
  call xc_f03_func_end(xc_func)
  ])"

dnl Note that the f03 suffix of the Fortran 2003 interface has been changed to f90 in Libxc 5.0.
dnl Then an interface with the f03 suffix identical to the one with the f90 suffix was again added in Libxc 5.1.
dnl In Libxc 6, only the f03 suffix will be kept.
testprog5="AC_LANG_PROGRAM([],[
  use xc_f03_lib_m
  implicit none
  integer :: major
  integer :: minor
  integer :: micro
  integer :: flags = XC_FLAGS_HAVE_ALL
  call xc_f03_version(major, minor, micro)]
  write(*,*) flags)"


FCFLAGS="$FCFLAGS_LIBXC $acx_libxc_save_FCFLAGS"

# set from environment variable, if not blank
if test ! -z "$LIBS_LIBXC"; then
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog4, [acx_libxc_ok=yes], [])
fi

# dynamic linkage, version 5.1 or greater
if test x"$acx_libxc_ok" = xno; then
  if test ! -z "$with_libxc_prefix"; then
    LIBS_LIBXC="-L$with_libxc_prefix/lib"
  else
    LIBS_LIBXC=""
  fi
  LIBS_LIBXC="$LIBS_LIBXC -lxcf03 -lxc"
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog5, [acx_libxc_ok=yes; acx_libxc_v5=yes], [])
  AC_RUN_IFELSE($testprog5_fxc, [], [acx_libxc_fxc=no])
  AC_RUN_IFELSE($testprog5_kxc, [], [acx_libxc_kxc=no])
fi

# dynamic linkage, version 4
if test x"$acx_libxc_ok" = xno; then
  if test ! -z "$with_libxc_prefix"; then
    LIBS_LIBXC="-L$with_libxc_prefix/lib"
  else
    LIBS_LIBXC=""
  fi
  LIBS_LIBXC="$LIBS_LIBXC -lxcf03 -lxc"
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog4, [acx_libxc_ok=yes], [])
fi

if test ! -z "$with_libxc_prefix"; then
  # static linkage, version 5.1 or greater
  if test x"$acx_libxc_ok" = xno; then
    LIBS_LIBXC="$with_libxc_prefix/lib/libxcf03.a $with_libxc_prefix/lib/libxc.a"
    LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
    AC_LINK_IFELSE($testprog5, [acx_libxc_ok=yes; acx_libxc_v5=yes], [])
    AC_RUN_IFELSE($testprog5_fxc, [], [acx_libxc_fxc=no])
    AC_RUN_IFELSE($testprog5_kxc, [], [acx_libxc_kxc=no])
  fi

  # static linkage, version 4
  if test x"$acx_libxc_ok" = xno; then
    LIBS_LIBXC="$with_libxc_prefix/lib/libxcf03.a $with_libxc_prefix/lib/libxc.a"
    LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
    AC_LINK_IFELSE($testprog4, [acx_libxc_ok=yes], [])
  fi
fi

AC_MSG_RESULT([$acx_libxc_ok ($FCFLAGS_LIBXC $LIBS_LIBXC)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libxc_ok" = xno; then
  AC_MSG_ERROR([Could not find required libxc library ( >= v 4.0.0 and != 5.0.0 ).])
fi

AC_MSG_CHECKING([whether libxc version is 5.1 or greater])
AC_MSG_RESULT([$acx_libxc_v5])

AC_MSG_CHECKING([whether libxc is compiled with fxc support])
AC_MSG_RESULT([$acx_libxc_fxc])

AC_MSG_CHECKING([whether libxc is compiled with kxc support])
AC_MSG_RESULT([$acx_libxc_kxc])

if test x"$acx_libxc_v5" = xyes; then
  AC_DEFINE(HAVE_LIBXC5, 1, [Defined if you have version 5.1 or greater of the LIBXC library.])
fi

if test x"$acx_libxc_fxc" = xyes; then
  AC_DEFINE(HAVE_LIBXC_FXC, 1, [Defined if you have compiled the LIBXC library with the FXC support.])
fi

if test x"$acx_libxc_kxc" = xyes; then
  AC_DEFINE(HAVE_LIBXC_KXC, 1, [Defined if you have compiled the LIBXC library with the KXC support.])
fi

AC_SUBST(FCFLAGS_LIBXC)
AC_SUBST(LIBS_LIBXC)
FCFLAGS="$acx_libxc_save_FCFLAGS"
LIBS="$acx_libxc_save_LIBS"
])dnl ACX_LIBXC
