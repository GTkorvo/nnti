
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(FEApp_ENABLE_Stokhos AND Stokhos_ENABLE_ML)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Linear_Diffusion_PCE_Example
    SOURCES linear_diffusion_pce.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Linear_Diffusion_PCE_Jacobi_Example
      SOURCES linear_diffusion_pce_jacobi.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      ) 

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Linear_Diffusion_PCE_GS_Example
      SOURCES linear_diffusion_pce_gs.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      )
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Linear_Diffusion_PCE_Jacobi_MRHS
      SOURCES linear_diffusion_pce_jacobi_mrhs.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      )
   
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Linear_Diffusion_PCE_Jacobi_Belos
      SOURCES linear_diffusion_pce_jacobi_belos.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      ) 
      
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Linear_Diffusion_PCE_Jacobi_Kkk_Example
      SOURCES linear_diffusion_pce_jacobi_Kkk.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      )

ENDIF()
