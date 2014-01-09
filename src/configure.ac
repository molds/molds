# Copyright (C) 2012-2014 Katsuhiko Nishimra
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

AC_PREREQ([2.60])
AC_INIT([MolDS], [0.2], [molds-dev@lists.sourceforge.jp], [MolDS])

dnl Check if pwd is source directory
AS_IF([test -x ./configure],
      AC_MSG_FAILURE([Don't execute configure script in source directory!]))

AC_CONFIG_SRCDIR([Main.cpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([no-define])

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
AC_PROG_CC([icc gcc])
AC_PROG_CXX([mpiicpc mpicxx])
AC_PROG_F77([ifort gfortran])

dnl Check if any MPI compiler is availbale
AX_MPI([],
       [AC_MSG_FAILURE([Cannot find any MPI compiler!])])

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
CXXFLAGS="$CXXFLAGS_SAVED"
AC_SUBST([OPENMP_CXXFLAGS])

dnl --with-openblas option
AC_ARG_WITH([openblas],
            [AC_HELP_STRING([--with-openblas[=DIR]],
                           [Use OpenBLAS library.])],
            [AS_CASE([$withval],
                     [yes],[BLAS_LIBS=-lopenblas
                            OPENBLAS_ROOT=],
                     [no],[OPENBLAS_ROOT=],
                     [AS_IF([test -d $withval],
                            [AC_MSG_NOTICE([Using OpenBLAS under `$withval'.])
                             BLAS_LIBS=-lopenblas
                             OPENBLAS_ROOT=$withval
                             CPPFLAGS="$CPPFLAGS -I$OPENBLAS_ROOT/include"
                             LDFLAGS="$LDFLAGS -L$OPENBLAS_ROOT/lib"],
                            [AC_MSG_FAILURE([Cannot find directory `$withval'!])])])],
            [OPENBLAS_ROOT=])

dnl Check for BLAS
AX_BLAS(,[AC_MSG_FAILURE([Cannot find BLAS library!])])
LIBS="$BLAS_LIBS $LIBS"
AC_CHECK_HEADERS([mkl.h cblas.h],[molds_blas_header_found=yes],)
AS_IF([test "x$molds_blas_header_found" != "xyes"],
      [AC_MSG_FAILURE([Cannot find mkl.h or blas.h!])])

dnl Check for CBLAS
AC_CHECK_FUNCS([cblas_dcopy cblas_daxpy cblas_ddot cblas_dgemv cblas_dgemm cblas_dsymv cblas_dsyr cblas_dsyrk],
                [molds_blas_cinterface_found=yes],
                [molds_blas_cinterface_found=no;break])

AS_IF([test x$molds_blas_cinterface_found != xyes],
      [AC_MSG_FAILURE([Cannot find C interface functions of BLAS!])])

dnl Check for LAPACK
AX_LAPACK(,[AC_MSG_FAILURE([Cannot find LAPACK library!])])
LIBS="$LAPACK_LIBS $LIBS"
AC_CHECK_HEADERS([mkl_lapacke.h lapacke.h],[molds_lapack_header_found=yes],)
AS_IF([test "x$molds_lapack_header_found" != "xyes"],
      [AC_MSG_FAILURE([Cannot find mkl_lapacke.h or lapacke.h!])])

dnl Check for LAPACK C interface
AC_CHECK_FUNCS([dgetrf dgetrs dsyevd dsysv],
               [molds_lapack_cinterface_found=yes],
               [molds_lapack_cinterface_found=no;break])

AS_IF([test x$molds_lapack_cinterface_found != xyes],
      [AC_SEARCH_LIBS([LAPACKE_dgetrf_work], [lapacke])])

AS_IF([test x$molds_lapack_cinterface_found != xyes],
      [AC_CHECK_FUNCS([LAPACKE_dgetrf_work LAPACKE_dgetrs_work LAPACKE_dsyevd_work LAPACKE_dsysv_work],
                      [molds_lapack_cinterface_found=yes],
                      [molds_lapack_cinterface_found=no;break])])

AS_IF([test x$molds_lapack_cinterface_found != xyes],
      [AC_MSG_FAILURE([Cannot find C interface functions of LAPACK!])])

dnl Check size of BLAS int
AC_LANG_PUSH([C])
AC_CACHE_CHECK([size of integer variable used in BLAS library],[molds_cv_sizeof_blasint],
               [AC_RUN_IFELSE([AC_LANG_PROGRAM([
                                                #include <stdint.h>
                                                #ifdef HAVE_DAXPY
                                                void daxpy(      int64_t*, double*, double*, int64_t*, double*, int64_t*);
                                                #elif defined(HAVE_CBLAS_DAXPY)
                                                void cblas_daxpy(int64_t , double , double*, int64_t , double*, int64_t);
                                                #endif
                                               ],[
                                                double alpha = 1,
                                                       x     = 1,
                                                       y     = 0;
                                                int64_t N    = 0x80000000,
                                                        incx = 0,
                                                        incy = 0;
                                                #ifdef HAVE_DAXPY
                                                daxpy(&N, &alpha, &x, &incx, &y, &incy);
                                                #elif defined(HAVE_CBLAS_DAXPY)
                                                cblas_daxpy(N, alpha, &x, incx, &y, incy);
                                                #else
                                                #error Cannot find either daxpy or cblas_daxpy!
                                                #endif
                                                return (int64_t)y == 0x80000000 ? 0 : 1;
                                               ])
                              ],[molds_cv_sizeof_blasint=64], [molds_cv_sizeof_blasint=32])
               ])
AC_LANG_POP([C])
AC_DEFINE_UNQUOTED([SIZEOF_BLASINT], [($molds_cv_sizeof_blasint)], [The size of integer variables used by BLAS library])

dnl Size of BLAS int
AC_DEFINE([SIZEOF_LAPACKINT], [SIZEOF_BLASINT], [The size of integer variables used by LAPACK library])
AS_IF([test x$molds_cv_sizeof_blasint = x64],
      [AC_DEFINE([MKL_ILP64], [],
                 [Whether ILP64 version of Intel Math Kernel Library])])

dnl Boost
AX_BOOST_BASE([1.46],,[AC_MSG_FAILURE([boost library is too old!])])
AX_BOOST_SERIALIZATION
AX_BOOST_MPI

dnl Check whether ld accepts -rpath option
AC_CACHE_CHECK([whether ld accepts -rpath option],[molds_cv_ld_accepts_rpath],
               [LDFLAGS_save=$LDFLAGS]
               [LIBS_save=$LIBS]
               [LDFLAGS="-Wl,-rpath=$BOOST_LIBDIR"]
               [LIBS=""]
               [AC_LINK_IFELSE([AC_LANG_PROGRAM(,[return 0;])],
                               [molds_cv_ld_accepts_rpath=yes],
                               [molds_cv_ld_accepts_rpath=no])]
               [LDFLAGS=$LDFLAGS_save]
               [LIBS=$LIBS_save])

AS_IF([test x$molds_cv_ld_accepts_rpath == xyes],
      [BOOST_LDFLAGS="$BOOST_LDFLAGS -Wl,-rpath=$BOOST_LIBDIR"])

dnl Check for Boost.Thread compiled library
AX_BOOST_THREAD

dnl Check for boost headers
CPPFLAGS_save=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
AC_CHECK_HEADERS([boost/format.hpp],,[AC_MSG_FAILURE([Cannot find boost/format.hpp!])])
AC_CHECK_HEADERS([boost/random.hpp],,[AC_MSG_FAILURE([Cannot find boost/random.hpp!])])
AC_CHECK_HEADERS([boost/shared_ptr.hpp],,[AC_MSG_FAILURE([Cannot find boost/shared_ptr.hpp!])])

AC_CHECK_HEADERS([boost/math/special_functions/factorials.hpp])
CPPFLAGS=$CPPFLAGS_save

dnl Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

dnl Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
