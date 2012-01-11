INCLUDE(TribitsAddExecutableAndTest)

SET(${PACKAGE_NAME}_TESTSIMPLE_SOURCES "")

SET(${PACKAGE_NAME}_TESTSIMPLE_SOURCES ${${PACKAGE_NAME}_TESTSIMPLE_SOURCES}
    SimpleApp_SimpleEpetraConstraints.cpp
    SimpleApp_SimpleEpetraDataPool.cpp
    SimpleApp_SimpleEpetraFeasStep.cpp
    SimpleApp_SimpleEpetraHessVec.cpp
    SimpleApp_SimpleEpetraLagMult.cpp
    SimpleApp_SimpleEpetraObjective.cpp
    testSimple.cpp
  )

IF(Trilinos_ENABLE_Amesos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    testSimple
    SOURCES ${${PACKAGE_NAME}_TESTSIMPLE_SOURCES}
    ARGS
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    )
ENDIF()
