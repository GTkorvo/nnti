/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

//----------------------------------------------------------------------------

namespace HybridFEM {

//----------------------------------------------------------------------------

template< typename ScalarType ,
          typename CoordScalarType >
struct GatherFill< 
  Kokkos::CrsMatrix< ScalarType , KOKKOS_MACRO_DEVICE > ,
  FEMesh< CoordScalarType , 8 , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE                              device_type ;
  typedef device_type::size_type                           size_type ;
  typedef Kokkos::CrsMatrix< ScalarType , device_type >    matrix_type ;
  typedef typename matrix_type::coefficients_type   coefficients_type ;
  typedef Kokkos::MultiVector< ScalarType , device_type >  vector_type ;
  typedef Kokkos::MDArray< ScalarType , device_type >      elem_contrib_type ;
  typedef Kokkos::MDArray< size_type , device_type >       elem_graph_type ;

  static const size_type ElemNodeCount = 8 ;

  typedef FEMesh< CoordScalarType , ElemNodeCount , device_type > mesh_type ;
  typedef typename mesh_type::node_elem_ids_type node_elem_ids_type ;

private:

  node_elem_ids_type node_elem_ids ;
  elem_graph_type    elem_graph ;
  elem_contrib_type  elem_matrices ;
  elem_contrib_type  elem_vectors ;
  coefficients_type  system_coeff ;
  vector_type        system_rhs ;

public:

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type irow ) const
  {
    const size_type node_elem_begin = node_elem_ids.row_entry_begin(irow);
    const size_type node_elem_end   = node_elem_ids.row_entry_end(irow);

    //  for each element that a node belongs to 

    for ( size_type i = node_elem_begin ; i < node_elem_end ; i++ ) {

      const size_type elem_id   = node_elem_ids( i, 0);
      const size_type row_index = node_elem_ids( i, 1);

      system_rhs(irow) += elem_vectors(elem_id, row_index);

      //  for each node in a particular related element  
      //  gather the contents of the element stiffness
      //  matrix that belong in irow

      for ( size_type j = 0 ; j < ElemNodeCount ; ++j ){
        const size_type A_index = elem_graph( elem_id , row_index , j );

        system_coeff( A_index ) += elem_matrices( elem_id, row_index, j );
      }
    }
  }


  static void apply( const matrix_type & matrix ,
                     const vector_type & rhs ,
                     const mesh_type   & mesh ,
                     const elem_graph_type   & elem_graph ,
                     const elem_contrib_type & elem_matrices ,
                     const elem_contrib_type & elem_vectors )
  {
    GatherFill op ;
    op.node_elem_ids = mesh.node_elem_ids ;
    op.elem_graph    = elem_graph ;
    op.elem_matrices = elem_matrices ;
    op.elem_vectors  = elem_vectors ;
    op.system_coeff  = matrix.coefficients ;
    op.system_rhs    = rhs ;

    parallel_for( matrix.graph.row_count() , op );
  }
};

//----------------------------------------------------------------------------

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


