INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.pu241.gcc.4.6.1.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/SubmitToTrilinos.cmake")

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_RELEASE_GCC)
SET(CTEST_TEST_TYPE Nightly)
#SET(CTEST_TEST_TIMEOUT 900)
SET(EXTRA_CONFIGURE_OPTIONS
  ${EXTRA_CONFIGURE_OPTIONS}
  -DSTK_ENABLE_BoostLib:BOOL=OFF
  )

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
