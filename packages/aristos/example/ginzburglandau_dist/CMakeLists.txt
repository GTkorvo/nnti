INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/example/ginzburglandau_dist)

SET(${PACKAGE_NAME}_TESTGLDIST_SOURCES "")

SET(${PACKAGE_NAME}_TESTGLDIST_SOURCES ${${PACKAGE_NAME}_TESTGLDIST_SOURCES}
    GLdistApp_GLdistYUEpetraConstraints.cpp
    GLdistApp_GLdistYUEpetraDataPool.cpp
    GLdistApp_GLdistYUEpetraFeasStep.cpp
    GLdistApp_GLdistYUEpetraHessVec.cpp
    GLdistApp_GLdistYUEpetraLagMult.cpp
    GLdistApp_GLdistYUEpetraObjective.cpp
    GLdistApp_SchurOp.cpp
    assemble_bdry.cpp
    assemble.cpp
    lassembly.cpp
    meshreader.cpp
    nonlinhessvec.cpp
    nonlinjac.cpp
    nonlinvec.cpp
    print2MATLAB.cpp
    quadrature.cpp
    usr_par.cpp
    testGLdist.cpp
  )

IF(Trilinos_ENABLE_Amesos AND Trilinos_ENABLE_Ifpack AND Trilinos_ENABLE_AztecOO AND Trilinos_ENABLE_Belos AND Trilinos_ENABLE_EpetraExt AND Trilinos_ENABLE_Stratimikos AND Trilinos_ENABLE_Thyra)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestGLdistFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/data
    SOURCE_FILES
      square.2.000
      square.2.001
      square.2.002
      square.2.003
      stratimikos_aztecoo_ifpack.xml
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestGLdistInstructions
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    SOURCE_FILES HOW_TO_RUN
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    testGLdist
    SOURCES ${${PACKAGE_NAME}_TESTGLDIST_SOURCES}
    ARGS "--geomfile=data/square.2"
    COMM mpi
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )


  TRIBITS_ADD_TEST(
    testGLdist
    POSTFIX_AND_ARGS_0 stratimikos
      --geomfile=data/square.2 --use-stratimikos
      --stratimikos-params-file=data/stratimikos_aztecoo_ifpack.xml
    COMM mpi
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

ENDIF()
