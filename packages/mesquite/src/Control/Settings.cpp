/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Settings.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Settings.hpp"
#include "LinearTriangle.hpp"
#include "LinearQuadrilateral.hpp"
#include "LinearTetrahedron.hpp"
#include "LinearPyramid.hpp"
#include "LinearPrism.hpp"
#include "LinearHexahedron.hpp"

namespace MESQUITE_NS {

#ifdef MSQ_TRAP_FPE
const bool IQ_TRAP_FPE_DEFAULT = true;
#else
const bool IQ_TRAP_FPE_DEFAULT = false;
#endif

struct SettingData {
  SettingData(); //!< Initialize to default settings.
  SettingData( const SettingData& other );
  SettingData& operator=( const SettingData& other );
  bool trapFPE;
  Settings::FixedVertexMode fixedMode;
  Settings::HigherOrderSlaveMode slaveMode;
  std::vector<const MappingFunction*> mapArray;
  std::vector<const MappingFunction2D*> mapArray2D;
  std::vector<const MappingFunction3D*> mapArray3D;
  
  LinearTriangle      linTriFunc;
  LinearQuadrilateral linQuadFunc;
  LinearTetrahedron   linTetFunc;
  LinearPyramid       linPyrFunc;
  LinearPrism         linPriFunc;
  LinearHexahedron    linHexFunc;
  
private:
  void fix_copy( const SettingData& other);
};

SettingData::SettingData()
  : trapFPE( IQ_TRAP_FPE_DEFAULT ),
    fixedMode( Settings::FIXED_FLAG ),
    slaveMode( Settings::SLAVE_ALL ),
    mapArray( MIXED, 0 ),
    mapArray2D( MIXED, 0 ),
    mapArray3D( MIXED, 0 )
{
  mapArray[TRIANGLE     ] = &linTriFunc;
  mapArray[QUADRILATERAL] = &linQuadFunc;
  mapArray[TETRAHEDRON  ] = &linTetFunc;
  mapArray[PYRAMID      ] = &linPyrFunc;
  mapArray[PRISM        ] = &linPriFunc;
  mapArray[HEXAHEDRON   ] = &linHexFunc;
  mapArray2D[TRIANGLE     ] = &linTriFunc;
  mapArray2D[QUADRILATERAL] = &linQuadFunc;
  mapArray3D[TETRAHEDRON  ] = &linTetFunc;
  mapArray3D[PYRAMID      ] = &linPyrFunc;
  mapArray3D[PRISM        ] = &linPriFunc;
  mapArray3D[HEXAHEDRON   ] = &linHexFunc;
}


SettingData::SettingData( const SettingData& other )
  : trapFPE( other.trapFPE ),
    fixedMode( other.fixedMode ),
    slaveMode( other.slaveMode ),
    mapArray( other.mapArray ),
    mapArray2D( other.mapArray2D ),
    mapArray3D( other.mapArray3D )
{
  fix_copy( other );
}

SettingData& SettingData::operator=( const SettingData& other )
{
  trapFPE = other.trapFPE;
  fixedMode = other.fixedMode;
  slaveMode = other.slaveMode;
  mapArray = other.mapArray;
  mapArray2D = other.mapArray2D;
  mapArray3D = other.mapArray3D;
  fix_copy( other );
  return *this;
}

void SettingData::fix_copy( const SettingData& other )
{
  if (mapArray[TRIANGLE] == &other.linTriFunc)
    mapArray[TRIANGLE] = &linTriFunc;
  if (mapArray[QUADRILATERAL] == &other.linQuadFunc)
    mapArray[QUADRILATERAL] = &linQuadFunc;
  if (mapArray[TETRAHEDRON] == &other.linTetFunc)
    mapArray[TETRAHEDRON] = &linTetFunc;
  if (mapArray[PYRAMID] == &other.linPyrFunc)
    mapArray[PYRAMID] = &linPyrFunc;
  if (mapArray[PRISM] == &other.linPriFunc)
    mapArray[PRISM] = &linPriFunc;
  if (mapArray[HEXAHEDRON] == &other.linHexFunc)
    mapArray[HEXAHEDRON] = &linHexFunc;
  if (mapArray2D[TRIANGLE] == &other.linTriFunc)
    mapArray2D[TRIANGLE] = &linTriFunc;
  if (mapArray2D[QUADRILATERAL] == &other.linQuadFunc)
    mapArray2D[QUADRILATERAL] = &linQuadFunc;
  if (mapArray3D[TETRAHEDRON] == &other.linTetFunc)
    mapArray3D[TETRAHEDRON] = &linTetFunc;
  if (mapArray3D[PYRAMID] == &other.linPyrFunc)
    mapArray3D[PYRAMID] = &linPyrFunc;
  if (mapArray3D[PRISM] == &other.linPriFunc)
    mapArray3D[PRISM] = &linPriFunc;
  if (mapArray3D[HEXAHEDRON] == &other.linHexFunc)
    mapArray3D[HEXAHEDRON] = &linHexFunc;
}

Settings::Settings() 
  : mData( new SettingData ) 
  { }
Settings::Settings( const Settings& other )
  : mData( new SettingData(other->mData) ) 
  { }
Settings::~Settings() 
  { delete mData; }
Settings& Settings::operator=( const Settings& other )
  { *mData = *(other.mData); return *this; }

void Settings::set_mapping_function( const MappingFunction* func )
{ 
  EntityTopology type = func->element_topology();
  if (mData->mapArray.size() <= (size_t)type)
    mData->mapArray.resize( type+1, 0 );
  mData->mapArray[type] = func;
  if (TopologyInfo::dimension(type) == 2 && mData->mapArray2D.size() > (size_t)type)
    mData->mapArray2D[type] = 0;
  else if (TopologyInfo::dimension(type) == 3 && mData->mapArray3D.size() > (size_t)type)
    mData->mapArray3D[type] = 0;
}

void Settings::set_mapping_function( const MappingFunction2D* func )
{ 
  unsigned type = func->element_topology();
  if (mData->mapArray.size() <= type)
    mData->mapArray.resize(type+1, 0);
  mData->mapArray[type] = func;
  if (mData->mapArray2D.size() <= type)
    mData->mapArray2D.resize( type+1, 0 );
  mData->mapArray2D[type] = func;
}

void Settings::set_mapping_function( const MappingFunction3D* func )
{ 
  unsigned type = func->element_topology();
  if (mData->mapArray.size() <= type)
    mData->mapArray.resize(type+1, 0);
  mData->mapArray[type] = func;
  if (mData->mapArray3D.size() <= type)
    mData->mapArray3D.resize( type+1, 0 );
  mData->mapArray3D[type] = func;
}

void Settings::set_mapping_functions( const MappingFunction* const* array,
                                      size_t array_len )
{
  for (size_t i = 0; i < array_len; ++i)
    set_mapping_function( array[i] );
}
void Settings::set_mapping_functions( const MappingFunction2D* const* array,
                                      size_t array_len )
{
  for (size_t i = 0; i < array_len; ++i)
    set_mapping_function( array[i] );
}
void Settings::set_mapping_functions( const MappingFunction3D* const* array,
                                      size_t array_len )
{
  for (size_t i = 0; i < array_len; ++i)
    set_mapping_function( array[i] );
}

const MappingFunction* Settings::get_mapping_function( EntityTopology type ) const
  { return (size_t)type < mData->mapArray.size() ? mData->mapArray[type] : 0; }

const MappingFunction2D* Settings::get_mapping_function_2D( EntityTopology type ) const
  { return (size_t)type < mData->mapArray2D.size() ? mData->mapArray2D[type] : 0; }

const MappingFunction3D* Settings::get_mapping_function_3D( EntityTopology type ) const
  { return (size_t)type < mData->mapArray3D.size() ? mData->mapArray3D[type] : 0; }

void Settings::set_fixed_vertex_mode( Settings::FixedVertexMode mode )
  { mData->fixedMode = mode; }

Settings::FixedVertexMode Settings::get_fixed_vertex_mode() const
  { return mData->fixedMode; }

void Settings::set_slaved_ho_node_mode( Settings::HigherOrderSlaveMode mode ) 
  { mData->slaveMode = mode; } 

Settings::HigherOrderSlaveMode Settings::get_slaved_ho_node_mode() const
  { return mData->slaveMode; }

void Settings::trap_floating_point_exception( bool enable )
  { mData->trapFPE = enable; }
bool Settings::trap_floating_point_exception() const
  { return mData->trapFPE; }



} // namespace Mesquite
