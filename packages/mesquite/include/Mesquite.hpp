#ifndef MESQUITE_HPP
#define MESQUITE_HPP

#ifdef WIN32
#pragma warning ( 4 : 4786)
#endif


#ifdef USE_STD_INCLUDES
#include <iostream>
#include <stdexcept>
#else
#include <iostream.h>
#include <stdexcept.h>
#endif

#ifdef USE_C_PREFIX_INCLUDES
#include <cmath>
#include <cstddef>
#include <cstring>
#else
#include <math.h>
#include <stddef.h>
#include <string.h>
#endif


#include "TSTT_Base.h"
/*! \file Mesquite.hpp
 */

/*!
  \namespace Mesquite
  
         Copyright 2002 Sandia Corporation. Under the terms of 
         Contract DE-AC04-94AL85000, there is a non-exclusive
         license for use of this work by or on behalf of the U.S.
         Government. Export of this program may require a license
         from the United States Government.
*/
namespace Mesquite
{
   
  typedef int StatusCode;

  typedef double real;

  enum StatusCodeValues
  {
    MSQ_FAILURE = 0,
    MSQ_SUCCESS
  };

  enum EntityTopology
  {
    POLYGON =7,
    TRIANGLE =8,
    QUADRILATERAL =9,
    POLYHEDRON =10,
    TETRAHEDRON =11,
    HEXAHEDRON =12,
    PRISM =13,
    PYRAMID =14,
    SEPTAHEDRON =15,
    MIXED
  };

    /*!Enumeration telling where the geometry can be found.
      NONE: No geometry is known.
      TSTT_DEFINED:  The geometry is defined through the TSTT interface.
      MSQ_DEFINED:  The geometry is stored in a SimplifiedGeometryEngine.
      */
  enum GeometryEngine
  {
    NONE,
    TSTT_DEFINED,
    MSQ_DEFINED
  };
  
  inline TSTT::EntityTopology mesquite_to_tstt(EntityTopology type)
  {
    switch (type)
    {
      case POLYGON:
        return TSTT::POLYGON;
      case TRIANGLE:
        return TSTT::TRIANGLE;
      case QUADRILATERAL:
        return TSTT::QUADRILATERAL;
      case POLYHEDRON:
        return TSTT::POLYHEDRON;
      case TETRAHEDRON:
        return TSTT::TETRAHEDRON;
      case HEXAHEDRON:
        return TSTT::HEXAHEDRON;
      case PRISM:
        return TSTT::PRISM;
      case PYRAMID:
        return TSTT::PYRAMID;
      case SEPTAHEDRON:
        return TSTT::SEPTAHEDRON;
      case MIXED:
      default:
        return TSTT::UNDEFINED;
    }
  }
  
  inline EntityTopology tstt_to_mesquite(TSTT::EntityTopology type)
  {
    switch (type)
    {
      case TSTT::POLYGON:
        return Mesquite::POLYGON;
      case TSTT::TRIANGLE:
        return Mesquite::TRIANGLE;
      case TSTT::QUADRILATERAL:
        return Mesquite::QUADRILATERAL;
      case TSTT::POLYHEDRON:
        return Mesquite::POLYHEDRON;
      case TSTT::TETRAHEDRON:
        return Mesquite::TETRAHEDRON;
      case TSTT::HEXAHEDRON:
        return Mesquite::HEXAHEDRON;
      case TSTT::PRISM:
        return Mesquite::PRISM;
      case TSTT::PYRAMID:
        return Mesquite::PYRAMID;
      case TSTT::SEPTAHEDRON:
        return Mesquite::SEPTAHEDRON;
      case TSTT::UNDEFINED:
      default:
        return Mesquite::MIXED;
    }
  }

  //GLOBALS variables
  const int MSQ_MAX_NUM_VERT_PER_ENT=8;
  const int MSQ_HIST_SIZE=7;//number of division in histogram
  const double MSQ_MIN=1.e-12;
  const double MSQ_MAX_CAP=1.e6;//Upper cap on metric values	
  static const double MSQ_SQRT_TWO = sqrt(2.0);
  static const double MSQ_SQRT_THREE = sqrt(3.0);
  static const double MSQ_SQRT_THREE_DIV_TWO=MSQ_SQRT_THREE/2.0;
  static const double MSQ_SQRT_THREE_INV=1.0/MSQ_SQRT_THREE;
  static const double MSQ_SQRT_TWO_INV=1.0/MSQ_SQRT_TWO;
  static const double MSQ_SQRT_TWO_DIV_SQRT_THREE=MSQ_SQRT_TWO/MSQ_SQRT_THREE;
  static const double MSQ_TWO_THIRDS = 2.0 / 3.0;
}

#define MSQ_DEBUG
#define MSQ_DEBUG4 // notify of all constructor / copy / destructor uses


/* the Mesquite Debugging system
     Level 0 provides no information and the debug macros are empty
     Level 1 provides user function information only
             e.g. threshold set, function used, etc.
     Level 2 provides basic algorithmic information for each local
             submesh
     Level 3 provides more information, data structures, and details
             than most users would want to know about 

     The default is Level 0
*/
//#define MSQ_DBG3

#ifdef MSQ_DBG3
#define MSQ_DEBUG_LEVEL 3
#elif  MSQ_DBG2
#define MSQ_DEBUG_LEVEL 2
#elif  MSQ_DBG1
#define MSQ_DEBUG_LEVEL 1
#elif  MSQ_DBG0
#define MSQ_DEBUG_LEVEL 0
#else
#define MSQ_DEBUG_LEVEL 0
#endif

#if MSQ_DEBUG_LEVEL != 0
#define MSQ_DEBUG_PRINT(level, statement)\
{\
   if ((level <= MSQ_DEBUG_LEVEL))\
     {\
     fprintf(stdout,statement);\
     fflush(stdout);\
     }\
}
#define MSQ_DEBUG_ACTION(level, action)\
{\
   if ((level <= MSQ_DEBUG_LEVEL))\
     {\
     action\
     fflush(stdout);\
     }\
}
#else
#define MSQ_DEBUG_PRINT(level, statement)\
{\
}
#define MSQ_DEBUG_ACTION(level, statement)\
{\
}
#endif

#endif
