#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([MolDS], [0.1], [molds-dev@lists.sourceforge.jp])
AM_INIT_AUTOMAKE([MolDS], [0.1])
AC_CONFIG_SRCDIR([Main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX([icpc g++])

# Use C++
AC_LANG_PUSH([C++])

#Check for pthread
AC_CHECK_LIB([pthread], [main])

#Check for OpenMP
AX_OPENMP(,[AC_MSG_FAILURE([Cannot use OpenMP!])])
export CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
AC_CHECK_HEADERS([omp.h],,[AC_MSG_FAILURE([Cannot find omp.h!])])
AC_CHECK_LIB([iomp5], [omp_set_num_threads],,[AC_MSG_FAILURE([Cannot find libiomp5!])])

# Check for mkl
export LIBS="-lmkl_intel_thread -lmkl_core $LIBS"
AC_CHECK_LIB([mkl_intel_lp64], [main],,
						 [
							AC_CHECK_LIB([mkl_intel], [main],,
													 [AC_MSG_FAILURE([Cannot find mkl!])
													 ])
							])

# Boost
AX_BOOST_BASE([1.48],,[AC_MSG_FAILURE([boost library is too old!])])

export CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
export LDFLAGS="$LDFLAGS -Wl,-rpath=$BOOST_LIBDIR $BOOST_LDFLAGS"

AC_CHECK_HEADERS([boost/random.hpp])
AX_CXX_CHECK_LIB([boost_random], [boost::random::random_device],,[AC_MSG_ERROR(Cannot find -lboost_random!)])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lmkl_core':
# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
