# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                            RBGen Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - package name, version number, and e-mail address below
AC_INIT(rbgen, 0.1d, hkthorn@sandia.gov)

# Hello World!
#  #auto np# - package name below
echo "----------------------------------------"
echo "Running RBGen Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #auto np# - file name below
AC_CONFIG_SRCDIR([src/RBGen_AnasaziPOD.cpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# #np# This a a subdirectory of Trilinos/packages/rbgen that contains
# #np# a number of useful m4 files.  Any package that uses these m4 files
# #np# will have to add them to the repository.  The m4 files perform various
# #np# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
#  #auto np# - Change file names in next line 
# Configure should create src/RBGen_config.h from src/RBGen_config.h.in

AM_CONFIG_HEADER(src/RBGen_config.h:src/RBGen_config.h.in)

# Allow users to specify their own "install" command.  If none is specified,
# the default is install-sh found in the config subdirectory.
                                                                                
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)
                                                                                
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake

#  - Required version of automake.
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check.  Please note that if you use
#	  F77_FUNC to determine Fortran name mangling, you should not remove
#	  the Fortran compiler check or the check for Fortran flags.  Doing
#	  so will prevent the detection of the proper name mangling in some
#	  cases.
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package.  Please 
#	  note that if you use F77_FUNC to determine Fortran name mangling, 
#	  you should not remove the Fortran compiler check or the check for 
#	  Fortran flags.  Doing so will prevent the detection of the proper 
#	  name mangling in some cases.
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

#  #np# - Customize the checks for headers needed in package being added.
# ------------------------------------------------------------------------
# Checks for header files - used in RBGen_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([vector vector.h], [break],)
AC_CHECK_HEADERS([stdexcept stdexcept.h], [break],)
AC_CHECK_HEADERS([algorithm algo.h algorithm.h], [break],)
AC_CHECK_HEADERS([iterator iterator.h], [break],)
AC_CHECK_HEADERS([functional function.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  RBGen requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  RBGen requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  RBGen requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  RBGen requires at least one or the other.]])
fi

if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cmath or math.h.  RBGen requires at least one or the other.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find vector or vector.h.  RBGen requires at least one or the other.]])
fi

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# #np# newp_swahili is an example of a package feature which requires both
# a cpp #define command and one or more source or header files  - ie 
# conditional compilation which requires a change to the makefile.
# If configure is called with "--enable-newp_swahili", it will place 
# a #define HAVE_NEWP_SWAHILI in RBGen_config.h and add source 
# and header files as specified in Makefile.am 
#
# To add an option to your package, you will need to modify this file
# (i.e. configure.ac), and src/Makefile.am.  A bootstrap is also required.
# Look for SWAHILI in both files listed above.
#
# TAC_ARG_ENABLE_OPTION performs several tasks, including:
#   1)  Setting a help string
#   2)  Printing a message to standard out
#   3)  Defining HAVE_NEWP_SWAHILI (For use in RBGen_config.h)
#
# AM_CONDITIONAL makes HAVE_NEWP_SWAHILI available for use in Makefile.am
# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#
#  #np# - These options are specific to Epetra and are left only as examples.
#         If uncommented, the following lines would define HAVE_FATAL_MESSAGES
#         in the (package)_config.h file.  Developers should use a 
#         more specific name such as 'rbgen-fatal_messages', rather
#         than 'fatal_messages'.
#TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
#                   FATAL_MESSAGES, yes)
#TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
#                   WARNING_MESSAGES, no)
#TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
#                   FORMAT_IO, yes)

#  #np# - These are specific to Epetra and are left only as examples.
#         If uncommented, the following line would define HAVE_ZOLTAN in
#         the (package)_config.h file.  Developers should use a more
#         specific name such as 'rbgen-zoltan', rather than 'zoltan'.
# TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
# TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)

# ------------------------------------------------------------------------
# Checks if the rbgen drivers should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE_SUB( rbgen, drivers, [Build RBGen drivers (binaries)], RBGEN_DRIVERS, yes)
AM_CONDITIONAL(BUILD_RBGEN_DRIVERS, [test "X$ac_cv_use_rbgen_drivers" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( rbgen, netcdf, [Enable NetCDF support], RBGEN_NETCDF, no)
AM_CONDITIONAL(BUILD_RBGEN_NETCDF, [test "X$ac_cv_use_rbgen_netcdf" != "Xno"])

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  #np# - These options can disable the tests and examples of a package.
#  #np# - Packages that do not have tests or examples should #-out the 
#  #np# - option(s) that does (do) not apply.

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( rbgen, tests, [Build RBGen tests], RBGEN_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_rbgen_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( rbgen, examples, [Build RBGen examples], RBGEN_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_rbgen_examples" != "Xno")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Trilinos Package Dependencies
# Required: Epetra, Teuchos
# Optional: EpetraExt, Anasazi
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(rbgen, epetra, RBGEN_EPETRA, BUILD_RBGEN_EPETRA, yes, [Build source code / examples which use Epetra.], [Build source code / examples which use Epetra. Can be overridden with --disable-rbgen-epetra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(rbgen, epetraext, RBGEN_EPETRAEXT, BUILD_RBGEN_EPETRAEXT, yes, [Build source code / examples which use EpetraExt.], [Build source code / examples which use EpetraExt. Can be overridden with --disable-rbgen-epetraext.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(rbgen, anasazi, RBGEN_ANASAZI, BUILD_RBGEN_ANASAZI, yes, [Build source code / examples which use Anasazi], [Build source code / examples which uses Anasazi. Can be overridden with --disable-rbgen-anasazi.])

# ------------------------------------------------------------------------
# Enable optional python wrappers
# ------------------------------------------------------------------------
#TAC_ARG_ENABLE_FEATURE_SUB( rbgen, python, [Build RBGen Python wrappers (requires --enable-python, which is NOT set by default)], RBGEN_PYTHON, yes)
#AM_CONDITIONAL(BUILD_RBGEN_PYTHON, test "X$ac_cv_use_rbgen_python" != "Xno")
#TAC_ENABLE_PYTHON

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK - trash WORKGXX 
AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

if test "X$ac_cv_use_rbgen_netcdf" != "Xno"; then
#
#  Check for NetCDF headers.
#
     AC_MSG_CHECKING(for netCDF headers (netcdf.h) ) 
     AC_TRY_CPP([#include "netcdf.h"],
	[AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find NetCDF headers."
     	 echo "Use --with-incdirs to specify NetCDF include directory."
     	 echo "-----"
     	 AC_MSG_ERROR(NetCDF headers cannot be found)
     	])
#
#  Check for NetCDF library.
#
     AC_MSG_CHECKING(for netDF library (libnetcdf.a) )
     AC_LANG_SAVE
     AC_LANG_CPLUSPLUS
     AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "netcdf.h"]], 
					[[int ncid=0; int status = nc_open("temp.file",NC_NOWRITE,&ncid);]]) ],
	[AC_MSG_RESULT(yes)],
	[AC_MSG_RESULT(no)
		 echo "-----"
	     	 echo "Cannot find NetCDF library."
     		 echo "Use --with-ldflags to specify NetCDF library directory."
	     	 echo "Use --with-libs to specify NetCDF library (-lnetcdf)."
	     	 echo "-----"
	     	 AC_MSG_ERROR(NetCDF test cannot link)])
     AC_LANG_RESTORE

fi



# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.

# #np# -
# If a package does not have tests or examples, the corresponding check(s)
# should be pulled out of the "if" statement below.
#if test "X$ac_cv_use_rbgen_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then
if test "X$ac_cv_use_rbgen_tests" != "Xno" || test "X$ac_cv_use_rbgen_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

# #np# If Lapack is also needed, see the  
# configure.ac file for epetra (Trilinos/packages/epetra)

ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

# -- DEC CXX --
#if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
#  case $target in
#    alpha*)
#      if test $CXX = cxx; then
#        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
#      fi
#    ;;
#  esac
#fi
                                                                  
#AC_MSG_CHECKING([for template-specific object files to be linked])if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
# AC_MSG_RESULT(none)
#else
# AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
#fi
                                                                  
#AC_SUBST(TRILINOS_TEMPLATE_OBJS)

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
##
#  #np# - Most packages should have an example and a full suite 
#  of tests.  New package comes only with src and example.
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
		Makefile
		Makefile.export.rbgen
		src/Makefile
		example/Makefile
		example/Drivers/Makefile
		example/Drivers/Epetra/Makefile
		example/Drivers/Epetra/Lena/Makefile
		example/Drivers/Epetra/Moon/Makefile
		example/Drivers/Epetra/sampleXML/Makefile
		])

AC_OUTPUT()

# Bye World!
#  #auto np# change package name
echo "----------------------------------------"
echo "Finished Running RBGen Configure Script"
echo "----------------------------------------"
