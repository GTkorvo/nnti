#ifndef MesquiteUtilities_hpp
#define MesquiteUtilities_hpp

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

#include "Mesquite.hpp"

#include "TSTT_Base.h"

#include "MesquiteError.hpp"

namespace Mesquite {

//   void writeVtkMesh(const char filebase[128], TSTT::cMesh_Handle mesh_h,
//                            MsqError &err);
//   void writeTSTTShowMeMesh(const char filebase[128], TSTT::cMesh_Handle mesh_h,
//                            MsqError &err);
//   void writeTSTTFacetMesh(const char filebase[128], TSTT::cMesh_Handle mesh_h,
//                            MsqError &err);

} // namespace

#endif // MesquiteUtilities_hpp
