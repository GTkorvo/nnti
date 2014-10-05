
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

APPEND_SET(HEADERS
  Aristos_ConfigDefs.hpp
  )

#
# Shared headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
APPEND_SET(HEADERS
  Aristos_ConfigDefs.h
  Aristos_Constraints.hpp
  Aristos_DataPool.hpp
  Aristos_EpetraVector.hpp
  Aristos_FeasStep.hpp
  Aristos_HessVec.hpp
  Aristos_LagMult.hpp
  Aristos_Objective.hpp
  Aristos_SQPAlgoDecl.hpp
  Aristos_SQPAlgo.hpp
  Aristos_Vector.hpp
  Aristos_YUEpetraVector.hpp
  )

#
# All sources
#

APPEND_SET(SOURCES
  Aristos_EpetraVector.cpp
  Aristos_YUEpetraVector.cpp
  Aristos_Version.cpp
  )


IF (TPL_ENABLE_MPI)
ENDIF()


#
# C) Define the targets for package's library(s)
# 

  TRIBITS_ADD_LIBRARY(
    aristos
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
