
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
  Solve
  SOURCES FEApp_Solve.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  SolveSG
  SOURCES FEApp_SGSolve.cpp
  COMM serial mpi
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
  DEST_FILES  input_nox_brusselator.xml
              input_loca_brusselator.xml 
              input_loca_bratu.xml
  EXEDEPS Solve
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(copySGTestInputFiles
  DEST_FILES  input_sg_nox_linear.xml
              input_mpni_nox_linear.xml
  	      input_sg_nox_hns.xml input_nox_hns.xml
  EXEDEPS SolveSG
  )


PACKAGE_ADD_TEST(
  Solve
  NAME NOX_Brusselator
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_nox_brusselator.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  ) 

PACKAGE_ADD_TEST(
  Solve
  NAME LOCA_Brusselator
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_loca_brusselator.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  )

PACKAGE_ADD_TEST(
  Solve
  NAME LOCA_Bratu
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_loca_bratu.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  )

PACKAGE_ADD_TEST(
  SolveSG
  NAME SG_NOX_Linear
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_sg_nox_linear.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  )

PACKAGE_ADD_TEST(
  SolveSG
  NAME MPNI_NOX_Linear
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_mpni_nox_linear.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  )

PACKAGE_ADD_TEST(
  SolveSG
  NAME SG_NOX_HNS
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS "--input=input_sg_nox_hns.xml --input_ig=input_nox_hns.xml"
  PASS_REGULAR_EXPRESSION "Test Passed"
  )
