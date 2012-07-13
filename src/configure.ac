# Copyright (C) 2012-2012 Katsuhiko Nishimra
#
# This file is part of MolDS.
#
# MolDS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MolDS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MolDS.  If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.59])
AC_INIT([MolDS], [0.1], [molds-dev@lists.sourceforge.jp])

dnl Check if pwd is source directory
AS_IF([test -x ./configure],
      AC_MSG_FAILURE([Don't execute configure script in source directory!]))

AM_INIT_AUTOMAKE([MolDS], [0.1])
AC_CONFIG_SRCDIR([Main.cpp])
AC_CONFIG_HEADERS([config.h])

dnl Setup CXXFLAGS
AC_MSG_CHECKING([release mode])
AC_ARG_ENABLE([release],
              AS_HELP_STRING([--enable-release],
                             [Enable optimize option and disable debug option])
              ,[AS_IF([test "$enableval" = "no"],
                      [release_enabled=no],
                      [release_enabled=yes])]
              ,[release_enabled=no])
AC_MSG_RESULT([$release_enabled])
AS_IF([test -z "$CXXFLAGS"],
      [AS_IF([test "$release_enabled" = "yes"],
             [CXXFLAGS="-O2"],
             [CXXFLAGS="-g -O0"])])

dnl Checks for programs.
AC_PROG_CXX([icpc g++])

dnl Use C++
AC_LANG_PUSH([C++])

dnl Check for pthread
AC_CHECK_LIB([pthread], [main])

dnl Check for OpenMP
CXXFLAGS_SAVED="$CXXFLAGS"
AX_OPENMP(,[AC_MSG_FAILURE([Cannot use OpenMP!])])
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
AC_CHECK_HEADERS([omp.h],,[AC_MSG_FAILURE([Cannot find omp.h!])])
AC_CHECK_LIB([iomp5], [omp_set_num_threads],,[
              AC_CHECK_FUNC([omp_set_num_threads],,
                            [AC_MSG_FAILURE([Cannot find omp_set_num_threads!])
                            ])
              ])
AC_CHECK_FUNC([omp_get_wtime],,
              [AC_MSG_FAILURE([Cannot find omp_get_wtime!])])

dnl OpenMP debug options
AS_IF([test "$release_enabled" = "no" ],
      [CXXFLAGS="$CXXFLAGS -openmp-report2"
       AC_MSG_CHECKING([for openmp diagnosis option])
       AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(){return 0;}])],
      [AS_IF([test -f penmp-report2],
             [rm penmp-report2
              AC_MSG_RESULT([no])],
             [OPENMP_CXXFLAGS="$OPENMP_CXXFLAGS -openmp-report2"
              AC_MSG_RESULT([-openmp-report2])])],
      [AC_MSG_RESULT([no])
       ])])
CXXFLAGS="$CXXFLAGS_SAVED"
AC_SUBST([OPENMP_CXXFLAGS])

dnl Check for LAPACK
AX_BLAS(,[AC_MSG_FAILURE([Cannot find BLAS library!])])
AX_LAPACK(,[AC_MSG_FAILURE([Cannot find LAPACK library!])])
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
AC_CHECK_HEADERS([mkl.h lapacke.h],[molds_lapack_header_found=yes],)
AS_IF([test "x$molds_lapack_header_found" != "xyes"],
      [AC_MSG_FAILURE([Cannot find mkl.h or lapacke.h!])])

dnl Boost
AX_BOOST_BASE([1.48],,[AC_MSG_FAILURE([boost library is too old!])])

CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS -Wl,-rpath=$BOOST_LIBDIR $BOOST_LDFLAGS"

AC_CHECK_HEADERS([boost/random.hpp])

dnl Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

dnl Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
