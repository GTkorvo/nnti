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

#ifndef Shards_CellTopologyData_h
#define Shards_CellTopologyData_h

#if defined( __cplusplus )
extern "C" {
#endif

/** \addtogroup shards_package_cell_topology
 *  \{
 */

/** \brief A simple 'C' struct of cell topology attributes.
 *
 *  The topology may be extended such that the number of nodes
 *  (subcells of dimension zero) is greater than the number of
 *  vertices.  In this case the vertices must be ordered first.
 *
 * Nodes, edges, and sides are subcells with a particular dimension.
 * A cell has edges and sides only if its dimension is greater than one.
 *  - node has Dim == 0
 *  - edge has Dim == 1
 *  - side has Dim == dimension - 1.
 */
struct CellTopologyData {
  /** \brief  Base, a.k.a. not-extended, version of this topology
   *          where vertex_count == edge_count.
   */
  const struct CellTopologyData * base ;

  /** \brief Intuitive name for this topology */
  const char * name ;

  /** \brief Unique key for this topology */
  unsigned key ;

  /** \brief Topological dimension */
  unsigned dimension ;

  /** \brief Number of vertices. */
  unsigned vertex_count ;

  /** \brief Number of nodes (a.k.a. Cell^0 subcells).
   *
   *  A topology is <em> extended </em> if node_count > vertex_count
   */
  unsigned node_count ;

  /** \brief Number of edges (a.k.a. Cell^1 subcells). */
  unsigned edge_count ;

  /** \brief Number of sides (a.k.a. Cell^(D-1) subcells). */
  unsigned side_count ;

  /** \brief Flag if the subcells of a given dimension are homogeneous */
  unsigned subcell_homogeneity ;

  /** \brief Number of subcells of each dimension. */
  unsigned subcell_count[4] ;

  /** \brief Subcell information.
   *  - required: 0 <= Ord <= subcell_count[Dim]
   *  - required: 0 <= J   <  subcell[Dim][Ord]->subcell_count[0]
   *  - subcell[Dim][Ord].topology
   *  - subcell[Dim][Ord].node[J]
   */
  struct Subcell {
    /** \brief Subcell topology */
    const struct CellTopologyData * topology ;
  
    /** \brief Subcell indexing of Cell^0 with respect to parent cell. */
    const unsigned * node ;
  };

  /** \brief  Array of subcells of each dimension
   *
   *  The length of each subcell array is subcell_count[Dim]
   *  - <b> subcell[Dim][Ord].topology </b> topology of the subcell
   *  - <b> subcell[Dim][Ord].node[I]  </b> node ordinal of the subcell's node I
   */
  const struct Subcell * subcell[4] ;

  /** \brief  Array of side subcells of length side_count
   *
   *  The length of the side array is side_count
   *  - <b> side[Ord].topology </b> topology of the side
   *  - <b> side[Ord].node[I]  </b> node ordinal of the side's node I
   */
  const struct Subcell * side ;

  /** \brief  Array of edges subcells of length edge_count
   *
   *  The length of the edge array is edge_count
   *  - <b> edge[Ord].topology </b> topology of the edge
   *  - <b> edge[Ord].node[I]  </b> node ordinal of the edge's node I
   */
  const struct Subcell * edge ;
};

/** \brief  Self-typedef */
typedef struct CellTopologyData  CellTopologyData ;

/** \} */

#if defined( __cplusplus )
} /* extern "C" */
#endif

#endif /* Shards_CellTopologyData_h */

