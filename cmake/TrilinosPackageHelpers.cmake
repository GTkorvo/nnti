
INCLUDE(Parse_Variable_Arguments)

#
# Macro called at the very beginning of a Trilinos package's top-level
# CMakeLists.txt file
#
MACRO(TRILINOS_PACKAGE PACKAGE_NAME_IN)

  IF (DEFINED PACKAGE_NAME_GLOBAL)
    IF (NOT ${PACKAGE_NAME_IN} STREQUAL ${PACKAGE_NAME_GLOBAL})
      MESSAGE(FATAL_ERROR "Error, the pacakge-defined package name '${PACKAGE_NAME_IN}' is not the same as the package name defined at the global level '${PACKAGE_NAME_GLOBAL}'")
    ENDIF()
  ENDIF()

   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS( PARSE "LANGUAGES" "" ${ARGN} )

  #
  # B) Set up the CMake "project" for this Trilinos package and define some
  # top-level varaibles.
  #

  SET(PACKAGE_NAME "${PACKAGE_NAME_IN}")
  MESSAGE(STATUS "Processing enabled Trilinos package: ${PACKAGE_NAME}")

  # Set up a new CMake "project" for this Trilinos "package"
  IF (PARSE_LANGUAGES)
    PROJECT(${PACKAGE_NAME} ${PARSE_LANGUAGES})
  ELSE()
    PROJECT(${PACKAGE_NAME})
  ENDIF()

  # Write PACKAGE versions of common variables
  SET(PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  SET(PACKAGE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # Get the name of the directory this Trilinos package is in
  FILE(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR} STANDARD_PACKAGE_SOURCE_DIR)
  STRING(REGEX REPLACE "/.+/(.+)" "\\1" PACKAGE_DIR_NAME "${STANDARD_PACKAGE_SOURCE_DIR}")

  #
  # C) Define package linkage varaibles
  #

  SET(${PACKAGE_NAME}_INCLUDE_DIRS "" CACHE INTERNAL "")
  SET(${PACKAGE_NAME}_LIBRARY_DIRS "" CACHE INTERNAL "")
  SET(${PACKAGE_NAME}_LIBRARIES "" CACHE INTERNAL "")

  #
  # D) Set up package dependancy-related variables
  #

  # ToDo: Implement!

  #
  # E) Define standard runtests targets for home-grown perl-based test harness
  #
  
  ADD_CUSTOM_TARGET(
    runtests-serial-${PACKAGE_DIR_NAME}
     ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
    --trilinos-dir=${TRILINOS_HOME_DIR}
    --comm=serial
    --build-dir=${TRILINOS_BUILD_DIR}
    --category=${TRILINOS_TEST_CATEGORY}
    --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
    --verbosity=1
    --packages=${PACKAGE_DIR_NAME}
    )
  
  ADD_CUSTOM_TARGET(
    runtests-mpi-${PACKAGE_DIR_NAME}
     ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
    --trilinos-dir=${TRILINOS_HOME_DIR}
    --comm=mpi
    --mpi-go="${TRILINOS_MPI_GO}"
    --max-proc=${MPIEXEC_MAX_NUMPROCS}
    --build-dir=${TRILINOS_BUILD_DIR}
    --category=${TRILINOS_TEST_CATEGORY}
    --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
    --verbosity=1
    --packages=${PACKAGE_DIR_NAME}
    )

ENDMACRO()


#
# Macro called at the very end of a Trilinos package's top-level
# CMakeLists.txt file
#
MACRO(TRILINOS_PACKAGE_POSTPROCESS)

  LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)
  LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARY_DIRS)
  LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARIES)

ENDMACRO()
