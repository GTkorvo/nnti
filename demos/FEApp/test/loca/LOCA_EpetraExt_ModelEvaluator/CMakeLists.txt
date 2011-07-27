
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_EpetraExt_ModelEvaluator_Test
    SOURCES 
    tcubed.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils 
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )
