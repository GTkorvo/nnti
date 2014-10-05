
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Brusselator_Example
    SOURCES brusselator.cpp
    COMM serial mpi
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Bratu_Example
    SOURCES bratu.cpp
    COMM serial mpi
    ARGS -v
    NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "All tests passed"
    )
