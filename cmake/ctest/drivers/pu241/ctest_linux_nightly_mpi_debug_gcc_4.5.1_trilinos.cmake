#
# Build all Secondary Stable Trilinos packages in core Trilinos with GCC 4.5.1 compiler
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.pu241.gcc.4.5.1.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/SubmitToTrilinos.cmake")

SET(COMM_TYPE MPI)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME MPI_DEBUG_GCC)
SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_TEST_TIMEOUT 900)
SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
SET(EXTRA_CONFIGURE_OPTIONS
  )

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
