/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#ifndef phdmesh_WedgeTopology_hpp
#define phdmesh_WedgeTopology_hpp

#include <element/LocalTopology.hpp>
#include <element/TriangleTopology.hpp>
#include <element/QuadrilateralTopology.hpp>

namespace phdmesh {

template< unsigned NumNode = 0 > struct Wedge ;        // 3D only

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 > ,
                IndexList< 1 , 2 > ,
                IndexList< 2 , 0 > ,
                IndexList< 3 , 4 > ,
                IndexList< 4 , 5 > ,
                IndexList< 5 , 3 > ,
                IndexList< 0 , 3 > ,
                IndexList< 1 , 4 > ,
                IndexList< 2 , 5 > >::type WedgeEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 , 3 > ,
                IndexList< 1 , 2 , 5 , 4 > ,
                IndexList< 0 , 3 , 5 , 2 > ,
                IndexList< 0 , 2 , 1 > ,
                IndexList< 3 , 4 , 5 > >::type WedgeSideNodeMap ;

template<>
struct Wedge<0>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                6 , 0 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      WedgeEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<> , 
                    BoundaryQuadrilateral<> ,
                    BoundaryQuadrilateral<> ,
                    BoundaryTriangle<> ,
                    BoundaryTriangle<> >::type ,
      WedgeSideNodeMap >
{
  typedef Wedge<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Wedge<6>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                6 , 6 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      WedgeEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<4> , 
                    BoundaryQuadrilateral<4> ,
                    BoundaryQuadrilateral<4> ,
                    BoundaryTriangle<3> ,
                    BoundaryTriangle<3> >::type ,
      WedgeSideNodeMap >
{
  typedef Wedge<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

