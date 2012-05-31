#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([MolDS], [0.1], [molds-dev@lists.sourceforge.jp])
AM_INIT_AUTOMAKE([MolDS], [0.1])
AC_CONFIG_SRCDIR([Main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lboost_random':
AC_CHECK_LIB([boost_random], [main])
# FIXME: Replace `main' with a function in `-liomp5':
AC_CHECK_LIB([iomp5], [main])
# FIXME: Replace `main' with a function in `-lmkl_core':
AC_CHECK_LIB([mkl_core], [main])
# FIXME: Replace `main' with a function in `-lmkl_intel':
AC_CHECK_LIB([mkl_intel], [main])
# FIXME: Replace `main' with a function in `-lmkl_intel_lp64':
AC_CHECK_LIB([mkl_intel_lp64], [main])
# FIXME: Replace `main' with a function in `-lmkl_intel_thread':
AC_CHECK_LIB([mkl_intel_thread], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
