
# Note: we should really be probing for these values, but will hard-code them
# in the short-term.  Most of these should be fine for all but ridiculously
# out-of-date C++ compilers, but you never know.

SET(HAVE_ALGORITHM TRUE)
SET(HAVE_CASSERT TRUE)
SET(HAVE_CCTYPE TRUE)
SET(HAVE_CERRNO TRUE)
SET(HAVE_CLIMITS TRUE)
SET(HAVE_CMATH TRUE)
SET(HAVE_COMPLEX TRUE)
SET(HAVE_CSTDARG TRUE)
SET(HAVE_CSTDIO TRUE)
SET(HAVE_CSTDLIB TRUE)
SET(HAVE_CSTRING TRUE)
SET(HAVE_IOMANIP TRUE)
SET(HAVE_IOSTREAM TRUE)
SET(HAVE_ITERATOR TRUE)
SET(HAVE_LIST TRUE)
SET(HAVE_MAP TRUE)
SET(HAVE_MEMORY TRUE)
SET(HAVE_MUTABLE TRUE)
SET(HAVE_NAMESPACES TRUE)
SET(HAVE_NEW_FOR_SCOPING TRUE)
SET(HAVE_NUMERIC TRUE)
SET(HAVE_NUMERIC_LIMITS TRUE)
SET(HAVE_POW TRUE)
SET(HAVE_SET TRUE)
SET(HAVE_SSTREAM TRUE)
SET(HAVE_FSTREAM TRUE)
SET(HAVE_STDEXCEPT TRUE)
SET(HAVE_STRING TRUE)
SET(HAVE_VECTOR TRUE)

# Enable compilers

IF(Trilinos_ENABLE_C)
  ENABLE_LANGUAGE(C)
ENDIF()

IF(Trilinos_ENABLE_CXX)
  ENABLE_LANGUAGE(CXX)
ENDIF()

IF(Trilinos_ENABLE_FORTRAN)
  ENABLE_LANGUAGE(Fortran)
ENDIF()

# Find Perl

FIND_PACKAGE(Perl)

# Find Python

FIND_PACKAGE(PythonInterp)

# Find the hostname used in selecting or deselecting tests by the
# TRILINOS_ADD_TEST(...) function.

SITE_NAME(Trilinos_HOSTNAME)
MARK_AS_ADVANCED(Trilinos_HOSTNAME)

# Set up MPI

IF(TPL_ENABLE_MPI)
  INCLUDE(ConfigureMPI)
ENDIF()

# Set some other options (NOTE: These should be handed in XXX_config.h files!

IF(CMAKE_SIZEOF_VOID_P GREATER 4)
  SET(EPETRA_ADDRESS64BIT ON)
ENDIF()

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
ENDIF()

# Do Fortran stuff

INCLUDE(TrilinosFortranMangling)

# Get BLAS name mangling
 
INCLUDE(TrilinosBLASMangling)

# Determine compiler type

INCLUDE(CMakeDetermineCXXCompiler)
PRINT_VAR(CMAKE_CXX_COMPILER_ID)
# See CMake/Modules/CMakeCXXCompilerId.cpp.in in the CMake source
# directory for a listing of know compiler types.
