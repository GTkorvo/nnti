
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(Stokhos_ENABLE_ML)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Linear_Diffusion_PCE_Example
    SOURCES linear_diffusion_pce.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

ENDIF()
