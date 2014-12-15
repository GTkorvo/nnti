cmake_minimum_required(VERSION 2.8.3)
cmake_policy(VERSION 2.8.3)
project(nnti)
include(CheckLibraryExists)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddOptionAndDefine)

include(TriosProbeTimers)
include(TriosProbeGettid)
include(TriosProbeNetwork)
include(TriosProbeSemaphores)
include(TriosProbeThreads)
include(TriosProbeXDR)

SET(PACKAGE_NAME nnti)
SET(PARENT_PACKAGE_NAME CMAKE)

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

CHECK_LIBRARY_EXISTS (ibverbs ibv_create_qp "" HAVE_IBVERBS)
CHECK_LIBRARY_EXISTS (ugni gni_err_str "" HAVE_GNI)

SET (DEPLIBS "")

SET(TRANSPORT_FOUND FALSE)

IF (HAVE_IBVERBS)
   SET(HAVE_TRIOS_INFINIBAND 1)
   SET(nnti_ENABLE_InfiniBand 1)
   SET(TRANSPORT_FOUND TRUE)
   INCLUDE_DIRECTORIES(/usr/include/infiniband)
   list (APPEND DEPLIBS ibverbs)
ENDIF (HAVE_IBVERBS)

IF (HAVE_GNI)
   SET(HAVE_TRIOS_GEMINI 1)
   SET(nnti_ENABLE_Gemini 1)
   SET(TRANSPORT_FOUND TRUE)
   list (APPEND DEPLIBS ugni)
ENDIF (HAVE_GNI)

IF (NOT ${TRANSPORT_FOUND})
   MESSAGE (FATAL_ERROR "No supported NNTI transport identified.  Cmake will exit.")
ENDIF (NOT ${TRANSPORT_FOUND})
    
ADD_SUBDIRECTORY( support )

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/support/xdr)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/support/logger)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/support/threads)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/support/timer)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/support/signal)

ADD_SUBDIRECTORY( nnti )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/Trios_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/Trios_config.h )
