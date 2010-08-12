
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_EpetraExt AND 
  FEApp_ENABLE_NOX AND FEApp_ENABLE_Stokhos AND Stokhos_ENABLE_ML)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    NOX_PCE_Example
    SOURCES heat_nonlinear_source_pce.cpp ENAT_NOXSolver.hpp ENAT_NOXSolver.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    NOX_Linear_Diffusion_PCE_Example
    SOURCES linear_diffusion_pce_nox.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS
    PASS_REGULAR_EXPRESSION "Test Passed"
    )
 
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      NOX_Linear_Diffusion_PCE_SGJacobi
      SOURCES linear_diffusion_pce_nox_jacobi.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      )
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      NOX_Linear_Diffusion_PCE_SGGS
      SOURCES linear_diffusion_pce_nox_gs.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      )

  IF(Stokhos_ENABLE_Anasazi)	
    PACKAGE_ADD_EXECUTABLE(
      NOX_Jacobian_KL_Example
      SOURCES hns_jacobian_KL.cpp ENAT_NOXSolver.hpp ENAT_NOXSolver.cpp
      COMM serial mpi
      )
  ENDIF()

ENDIF()
