/*------------------------------------------------------------------------*/
/*                  shards : Shared Discretization Tools                  */
/*                Copyright (2008) Sandia Corporation                     */
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

#include <iostream>

#include <Shards_BasicTopologies.hpp>

extern "C" {
typedef CellTopologyData::Subcell  Subcell ;
}

namespace shards {

enum { MAXIMUM_INDICES = 256 };

const unsigned * index_identity_array()
{
  static unsigned self[ MAXIMUM_INDICES ];
  static int init = 1 ;

  if ( init ) {
    for ( unsigned i = 0 ; i < MAXIMUM_INDICES ; ++i ) { self[i] = i ; }
    init = 0 ;
  }

  return self ;
}

const Subcell * subcell_nodes_array()
{
  static Subcell self[ MAXIMUM_INDICES ];
  static int init = 1 ;

  if ( init ) {

    const CellTopologyData * const top = getCellTopologyData<Node>();
    const unsigned     * const ident = index_identity_array();

    for ( int i = 0 ; i < MAXIMUM_INDICES ; ++i ) {
      self[i].topology = top ;
      self[i].node     = ident + i ;
    }

    init = 0 ;
  }

  return self ;
}

namespace {

template< class IList >
const unsigned * index_list( const IList & )
{
  static const unsigned self[] = {
    IndexListAt< IList ,  0 >::value ,
    IndexListAt< IList ,  1 >::value ,
    IndexListAt< IList ,  2 >::value ,
    IndexListAt< IList ,  3 >::value ,
    IndexListAt< IList ,  4 >::value ,
    IndexListAt< IList ,  5 >::value ,
    IndexListAt< IList ,  6 >::value ,
    IndexListAt< IList ,  7 >::value ,
    IndexListAt< IList ,  8 >::value ,
    IndexListAt< IList ,  9 >::value ,
    IndexListAt< IList , 10 >::value ,
    IndexListAt< IList , 11 >::value ,
    IndexListAt< IList , 12 >::value ,
    IndexListAt< IList , 13 >::value ,
    IndexListAt< IList , 14 >::value ,
    IndexListAt< IList , 15 >::value ,
    IndexListAt< IList , 16 >::value ,
    IndexListAt< IList , 17 >::value ,
    IndexListAt< IList , 18 >::value ,
    IndexListAt< IList , 19 >::value ,
    IndexListAt< IList , 20 >::value ,
    IndexListAt< IList , 21 >::value ,
    IndexListAt< IList , 22 >::value ,
    IndexListAt< IList , 23 >::value ,
    IndexListAt< IList , 24 >::value ,
    IndexListAt< IList , 25 >::value ,
    IndexListAt< IList , 26 >::value ,
    IndexListAt< IList , 27 >::value ,
    IndexListAt< IList , 28 >::value ,
    IndexListAt< IList , 29 >::value ,
    IndexListAt< IList , 30 >::value ,
    IndexListAt< IList , 31 >::value
  };

  return self ;
}

//----------------------------------------------------------------------

template< class TList , class IList , unsigned N >
struct SubcellValue ;

template< class TList , class IList >
struct SubcellValue<TList,IList,0>
{ static void assign( Subcell * ) {} };

template< class TList , class IList , unsigned N >
struct SubcellValue {
  static void assign( Subcell * s )
    {
      enum { I = N - 1 };
      SubcellValue<TList,IList,I>::assign( s );
      s[I].topology = getCellTopologyData< typename TypeListAt<TList,I>::type >();
      s[I].node     = index_list( typename TypeListAt<IList,I>::type() );
    }
};

template< class TList , class IList >
struct SubcellArray {

  enum { N = TypeListLength<TList>::value };

  Subcell array[ N ];

  SubcellArray() { SubcellValue<TList,IList,N>::assign( array ); }
};
  
//----------------------------------------------------------------------

template< class Traits > struct Descriptor ;

template< unsigned Number_Vertex ,
          unsigned Number_Node ,
          class    EdgeList ,
          class    EdgeMaps ,
          class    FaceList ,
          class    FaceMaps >
struct Descriptor<
  CellTopologyTraits< 3 , Number_Vertex , Number_Node ,
                      EdgeList , EdgeMaps , FaceList , FaceMaps > >
{
  typedef CellTopologyTraits< 3 , Number_Vertex , Number_Node ,
                              EdgeList , EdgeMaps ,
                              FaceList , FaceMaps > Traits ;

  typedef SubcellArray< EdgeList , EdgeMaps > EdgeArray ;
  typedef SubcellArray< FaceList , FaceMaps > FaceArray ;

  EdgeArray edges ;
  FaceArray faces ;

  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name )
    : edges(), faces()
    {
      typedef typename Traits::template subcell<0> subcell_0 ;
      typedef typename Traits::template subcell<1> subcell_1 ;
      typedef typename Traits::template subcell<2> subcell_2 ;
      typedef typename Traits::template subcell<3> subcell_3 ;

      self.topology = & top ;
      self.node     = index_identity_array();

      top.base             = base ? base : & top ;
      top.name             = name ;
      top.key              = Traits::key ;
      top.dimension        = 3 ;
      top.vertex_count     = Number_Vertex ;
      top.node_count       = Number_Node ;
      top.edge_count       = EdgeArray::N ;
      top.side_count       = Traits::side_count ;
      top.subcell_homogeneity[0] = subcell_0::homogeneity ;
      top.subcell_homogeneity[1] = subcell_1::homogeneity ;
      top.subcell_homogeneity[2] = subcell_2::homogeneity ;
      top.subcell_homogeneity[3] = subcell_3::homogeneity ;
      top.subcell_count[0] = Number_Node ;
      top.subcell_count[1] = EdgeArray::N ;
      top.subcell_count[2] = FaceArray::N ;
      top.subcell_count[3] = 1 ;
      top.subcell[0]       = subcell_nodes_array();
      top.subcell[1]       = edges.array ;
      top.subcell[2]       = faces.array ;
      top.subcell[3]       = & self ;
      top.side             = faces.array ;
      top.edge             = edges.array ;
    };
};

template< unsigned Number_Vertex ,
          unsigned Number_Node ,
          class    EdgeList ,
          class    EdgeMaps >
struct Descriptor<
  CellTopologyTraits< 2 , Number_Vertex , Number_Node ,
                      EdgeList , EdgeMaps , TypeListEnd , TypeListEnd > >
{
  typedef CellTopologyTraits< 2 , Number_Vertex , Number_Node ,
                              EdgeList , EdgeMaps ,
                              TypeListEnd , TypeListEnd > Traits ;
  
  typedef SubcellArray< EdgeList , EdgeMaps > EdgeArray ;

  EdgeArray edges ;

  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name ) : edges()
    {
      typedef typename Traits::template subcell<0> subcell_0 ;
      typedef typename Traits::template subcell<1> subcell_1 ;
      typedef typename Traits::template subcell<2> subcell_2 ;
      typedef typename Traits::template subcell<3> subcell_3 ;

      self.topology = & top ;
      self.node     = index_identity_array();

      top.base             = base ? base : & top ;
      top.name             = name ;
      top.key              = Traits::key ;
      top.dimension        = 2 ;
      top.vertex_count     = Number_Vertex ;
      top.node_count       = Number_Node ;
      top.edge_count       = EdgeArray::N ;
      top.side_count       = Traits::side_count ;
      top.subcell_homogeneity[0] = subcell_0::homogeneity ;
      top.subcell_homogeneity[1] = subcell_1::homogeneity ;
      top.subcell_homogeneity[2] = subcell_2::homogeneity ;
      top.subcell_homogeneity[3] = subcell_3::homogeneity ;
      top.subcell_count[0] = Number_Node ;
      top.subcell_count[1] = EdgeArray::N ;
      top.subcell_count[2] = 1 ;
      top.subcell_count[3] = 0 ;
      top.subcell[0]       = subcell_nodes_array();
      top.subcell[1]       = edges.array ;
      top.subcell[2]       = & self ;
      top.subcell[3]       = NULL ;
      top.side             = edges.array ;
      top.edge             = edges.array ;
    };
};

template< unsigned Number_Node , unsigned Number_Vertex >
struct Descriptor<
  CellTopologyTraits< 1 , Number_Vertex , Number_Node ,
                      TypeListEnd , TypeListEnd ,
                      TypeListEnd , TypeListEnd > >
{
  typedef CellTopologyTraits< 1 , Number_Vertex , Number_Node ,
                              TypeListEnd , TypeListEnd ,
                              TypeListEnd , TypeListEnd > Traits ;
  
  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name )
    {
      typedef typename Traits::template subcell<0> subcell_0 ;
      typedef typename Traits::template subcell<1> subcell_1 ;
      typedef typename Traits::template subcell<2> subcell_2 ;
      typedef typename Traits::template subcell<3> subcell_3 ;

      self.topology = & top ;
      self.node     = index_identity_array();

      top.base             = base ? base : & top ;
      top.name             = name ;
      top.key              = Traits::key ;
      top.dimension        = 1 ;
      top.vertex_count     = Number_Vertex ;
      top.node_count       = Number_Node ;
      top.edge_count       = 0 ;
      top.side_count       = 0 ;
      top.subcell_homogeneity[0] = subcell_0::homogeneity ;
      top.subcell_homogeneity[1] = subcell_1::homogeneity ;
      top.subcell_homogeneity[2] = subcell_2::homogeneity ;
      top.subcell_homogeneity[3] = subcell_3::homogeneity ;
      top.subcell_count[0] = Number_Node ;
      top.subcell_count[1] = 1 ;
      top.subcell_count[2] = 0 ;
      top.subcell_count[3] = 0 ;
      top.subcell[0]       = subcell_nodes_array();
      top.subcell[1]       = & self ;
      top.subcell[2]       = NULL ;
      top.subcell[3]       = NULL ;
      top.side             = NULL ;
      top.edge             = NULL ;
    };
};

template<>
struct Descriptor<
  CellTopologyTraits< 0 , 0 , 0 ,
                      TypeListEnd , TypeListEnd ,
                      TypeListEnd , TypeListEnd > >
{
  typedef CellTopologyTraits< 0 , 0 , 0 ,
                              TypeListEnd , TypeListEnd ,
                              TypeListEnd , TypeListEnd > Traits ;
  
  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name )
    {
      self.topology = & top ;
      self.node     = index_identity_array();

      top.base             = base ? base : & top ;
      top.name             = name ;
      top.key              = Traits::key ;
      top.dimension        = 0 ;
      top.vertex_count     = 0 ;
      top.node_count       = 0 ;
      top.edge_count       = 0 ;
      top.side_count       = Traits::side_count ;
      top.subcell_homogeneity[0] = true ;
      top.subcell_homogeneity[1] = false ;
      top.subcell_homogeneity[2] = false ;
      top.subcell_homogeneity[3] = false ;
      top.subcell_count[0] = 1 ;
      top.subcell_count[1] = 0 ;
      top.subcell_count[2] = 0 ;
      top.subcell_count[3] = 0 ;
      top.subcell[0]       = & self ;
      top.subcell[1]       = NULL ;
      top.subcell[2]       = NULL ;
      top.subcell[3]       = NULL ;
      top.side             = NULL ;
      top.edge             = NULL ;
    };
};

}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const CellTopologyData & t )
{
  s << t.name ;
  s << " { D = " << t.dimension ;
  s << " , NV = " << t.vertex_count ;
  s << " , K = 0x" << std::hex << t.key << std::dec ;
  s << std::endl ;

  for ( unsigned d = 0 ; d < 4 ; ++d ) {
    for ( unsigned i = 0 ; i < t.subcell_count[d] ; ++i ) {

      const Subcell & sub = t.subcell[d][i] ;

      s << "  subcell[" << d << "][" << i << "] = { " ;

      s << sub.topology->name ;
      s << " ," ;
      for ( unsigned j = 0 ; j < sub.topology->node_count ; ++j ) {
        s << " " << sub.node[j] ;
      }
      s << " }" << std::endl ;
    }
  }

  s << "}" << std::endl ;
  return s ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Node>()
{
  static const char name[] = "Node" ;
  static const Descriptor< Node::Traits > self( NULL , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData< Particle >()
{
  static const char name[] = "Particle" ;
  static const Descriptor< Particle::Traits > self( NULL , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Line<2> >()
{
  static const char name[] = "Line_2" ;
  static const Descriptor< Line<2>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< Line<3> >()
{
  static const char name[] = "Line_3" ;
  static const Descriptor< Line<3>::Traits >
    self( getCellTopologyData<Line<2> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< Beam<2> >()
{
  static const char name[] = "Beam_2" ;
  static const Descriptor< Beam<2>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< Beam<3> >()
{
  static const char name[] = "Beam_3" ;
  static const Descriptor< Beam<3>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< ShellLine<2> >()
{
  static const char name[] = "ShellLine_2" ;
  static const Descriptor< ShellLine<2>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< ShellLine<3> >()
{
  static const char name[] = "ShellLine_3" ;
  static const Descriptor< ShellLine<3>::Traits >
    self( getCellTopologyData< ShellLine<2> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Triangle<3> >()
{
  static const char name[] = "Triangle_3" ;
  static const Descriptor< Triangle<3>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Triangle<4> >()
{
  static const char name[] = "Triangle_4" ;
  static const Descriptor< Triangle<4>::Traits >
    self( getCellTopologyData< Triangle<3> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Triangle<6> >()
{
  static const char name[] = "Triangle_6" ;
  static const Descriptor< Triangle<6>::Traits >
    self( getCellTopologyData< Triangle<3> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<ShellTriangle<3> >()
{
  static const char name[] = "ShellTriangle_3" ;
  static const Descriptor< ShellTriangle<3>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<ShellTriangle<6> >()
{
  static const char name[] = "ShellTriangle_6" ;
  static const Descriptor< ShellTriangle<6>::Traits >
    self( getCellTopologyData< ShellTriangle<3> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Quadrilateral<4> >()
{
  static const char name[] = "Quadrilateral_4" ;
  static const Descriptor< Quadrilateral<4>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Quadrilateral<8> >()
{
  static const char name[] = "Quadrilateral_8" ;
  static const Descriptor< Quadrilateral<8>::Traits >
    self( getCellTopologyData<Quadrilateral<4> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Quadrilateral<9> >()
{
  static const char name[] = "Quadrilateral_9" ;
  static const Descriptor< Quadrilateral<9>::Traits >
    self( getCellTopologyData<Quadrilateral<4> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<ShellQuadrilateral<4> >()
{
  static const char name[] = "ShellQuadrilateral_4" ;
  static const Descriptor< ShellQuadrilateral<4>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<ShellQuadrilateral<8> >()
{
  static const char name[] = "ShellQuadrilateral_8" ;
  static const Descriptor< ShellQuadrilateral<8>::Traits >
    self( getCellTopologyData<ShellQuadrilateral<4> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<ShellQuadrilateral<9> >()
{
  static const char name[] = "ShellQuadrilateral_9" ;
  static const Descriptor< ShellQuadrilateral<9>::Traits >
    self( getCellTopologyData<ShellQuadrilateral<4> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Hexahedron<8> >()
{
  static const char name[] = "Hexahedron_8" ;
  static const Descriptor< Hexahedron<8>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Hexahedron<20> >()
{
  static const char name[] = "Hexahedron_20" ;
  static const Descriptor< Hexahedron<20>::Traits >
    self( getCellTopologyData<Hexahedron<8> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Hexahedron<27> >()
{
  static const char name[] = "Hexahedron_27" ;
  static const Descriptor< Hexahedron<27>::Traits >
    self( getCellTopologyData<Hexahedron<8> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Tetrahedron<4> >()
{
  static const char name[] = "Tetrahedron_4" ;
  static const Descriptor< Tetrahedron<4>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Tetrahedron<10> >()
{
  static const char name[] = "Tetrahedron_10" ;
  static const Descriptor< Tetrahedron<10>::Traits >
    self( getCellTopologyData<Tetrahedron<4> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Tetrahedron<8> >()
{
  static const char name[] = "Tetrahedron_8" ;
  static const Descriptor< Tetrahedron<8>::Traits >
    self( getCellTopologyData<Tetrahedron<4> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Pyramid<5> >()
{
  static const char name[] = "Pyramid_5" ;
  static const Descriptor< Pyramid<5>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Pyramid<13> >()
{
  static const char name[] = "Pyramid_13" ;
  static const Descriptor< Pyramid<13>::Traits >
    self( getCellTopologyData<Pyramid<5> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Pyramid<14> >()
{
  static const char name[] = "Pyramid_14" ;
  static const Descriptor< Pyramid<14>::Traits >
    self( getCellTopologyData<Pyramid<5> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Wedge<6> >()
{
  static const char name[] = "Wedge_6" ;
  static const Descriptor< Wedge<6>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Wedge<15> >()
{
  static const char name[] = "Wedge_15" ;
  static const Descriptor< Wedge<15>::Traits >
    self( getCellTopologyData<Wedge<6> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Wedge<18> >()
{
  static const char name[] = "Wedge_18" ;
  static const Descriptor< Wedge<18>::Traits >
    self( getCellTopologyData<Wedge<6> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Pentagon<5> >()
{
  static const char name[] = "Pentagon_5" ;
  static const Descriptor< Pentagon<5>::Traits > self( NULL , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Hexagon<6> >()
{
  static const char name[] = "Hexagon_6" ;
  static const Descriptor< Hexagon<6>::Traits > self( NULL , name );
  return & self.top ;
}


}//namespace shards

