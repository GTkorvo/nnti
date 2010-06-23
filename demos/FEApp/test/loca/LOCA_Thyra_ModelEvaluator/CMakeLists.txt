
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA AND NOX_ENABLE_LOCA AND
   NOX_ENABLE_EpetraExt AND FEApp_ENABLE_NOX)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Thyra_ModelEvaluator_Test
    SOURCES 
    tcubed.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils 
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Thyra_ModelEvaluator_TurningPoint_Test
    SOURCES 
    tcubedTP.cpp
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()
