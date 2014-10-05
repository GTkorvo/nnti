

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NOX_Brusselator_Example
    SOURCES brusselator.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS -v
    PASS_REGULAR_EXPRESSION "Test Passed"
    )
