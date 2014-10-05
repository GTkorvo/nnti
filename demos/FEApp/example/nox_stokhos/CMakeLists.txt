
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(Stokhos_ENABLE_ML)
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NOX_PCE_Example
    SOURCES heat_nonlinear_source_pce.cpp     
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NOX_Linear_Diffusion_PCE_Example
    SOURCES linear_diffusion_pce_nox.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NOX_Linear_Diffusion_PCE_MPNI_Example
    SOURCES linear_diffusion_pce_mpni.cpp 
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

ENDIF()

IF( Stokhos_ENABLE_Ifpack AND Stokhos_ENABLE_Anasazi )

  TRIBITS_ADD_EXECUTABLE(
    NOX_Jacobian_KL_Example
    SOURCES hns_jacobian_KL.cpp
    COMM serial mpi
    )

ENDIF()
