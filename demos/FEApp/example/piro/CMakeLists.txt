
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  Solve
  SOURCES FEApp_Solve.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  SGSolve
  SOURCES FEApp_SGSolve.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  SGAnalysis
  SOURCES FEApp_SGAnalysis.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
  DEST_FILES  input_nox_brusselator.xml
              input_loca_brusselator.xml 
              input_loca_bratu.xml
  EXEDEPS Solve
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(copySGTestInputFiles
  DEST_FILES  input_sg_nox_linear.xml
              input_mpni_nox_linear.xml
  	      input_sg_nox_hns.xml input_nox_hns.xml
  EXEDEPS SGSolve
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(copySGAnalysisTestInputFiles
  DEST_FILES  input_sg_analysis_linear.xml
              dakota_sg_linear.in dakota_sg_linear_target.txt
  EXEDEPS SGSolve
  )

IF (Piro_ENABLE_TriKota)
  TRIBITS_ADD_EXECUTABLE(
    Dakota
    SOURCES FEApp_Dakota.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyDakotaTestInputFiles
    DEST_FILES  input_dakota_pce_mp_linear.xml dakota_pce_mp_linear.in
    EXEDEPS Dakota
    )
ENDIF()

TRIBITS_ADD_TEST(
  Solve
  NAME NOX_Brusselator
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_nox_brusselator.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  ) 

TRIBITS_ADD_TEST(
  Solve
  NAME LOCA_Brusselator
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_loca_brusselator.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  )

TRIBITS_ADD_TEST(
  Solve
  NAME LOCA_Bratu
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_loca_bratu.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  )

IF(Stokhos_ENABLE_ML)

  TRIBITS_ADD_TEST(
    SGSolve
    NAME SG_NOX_Linear
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS "--input=input_sg_nox_linear.xml"
    PASS_REGULAR_EXPRESSION "Test Passed"
    )
  
  TRIBITS_ADD_TEST(
    SGSolve
    NAME MPNI_NOX_Linear
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS "--input=input_mpni_nox_linear.xml"
    PASS_REGULAR_EXPRESSION "Test Passed"
    )
  
  TRIBITS_ADD_TEST(
    SGSolve
    NAME SG_NOX_HNS
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS "--input=input_sg_nox_hns.xml --input_ig=input_nox_hns.xml"
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

ENDIF()