/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

namespace HybridFEM {
namespace NonLinear {

//----------------------------------------------------------------------------

template< typename ScalarCoordType , unsigned ElemNode , typename ScalarType >
struct ElementComputation<
  FEMesh< ScalarCoordType , ElemNode , KOKKOSARRAY_MACRO_DEVICE > , ScalarType >
{
  typedef KOKKOSARRAY_MACRO_DEVICE  device_type;
  typedef ScalarType           scalar_type ;

  static const unsigned ElementNodeCount = ElemNode ;

  typedef FEMesh< ScalarCoordType , ElementNodeCount , device_type > mesh_type ;

  typedef HexElement_Data< ElementNodeCount > element_data_type ;

  static const unsigned SpatialDim       = element_data_type::spatial_dimension ;
  static const unsigned FunctionCount    = element_data_type::function_count ;
  static const unsigned IntegrationCount = element_data_type::integration_count ;
  static const unsigned TensorDim        = SpatialDim * SpatialDim ;

  typedef KokkosArray::View< scalar_type[][FunctionCount][FunctionCount] , device_type > elem_matrices_type ;
  typedef KokkosArray::View< scalar_type[][FunctionCount] , device_type > elem_vectors_type ;
  typedef KokkosArray::View< scalar_type[] , device_type > value_vector_type ;


private:

  const element_data_type                 elem_data ;
  typename mesh_type::elem_node_ids_type  elem_node_ids ;
  typename mesh_type::node_coords_type    node_coords ;
  value_vector_type                       nodal_values ;
  elem_matrices_type                      element_matrices ;
  elem_vectors_type                       element_vectors ;
  scalar_type                             coeff_K ;

public:

  ElementComputation( const mesh_type   & arg_mesh ,
                      const elem_matrices_type  & arg_element_matrices , 
                      const elem_vectors_type   & arg_element_vectors ,
                      const value_vector_type   & arg_nodal_values ,
	              const scalar_type   arg_coeff_K )
  : elem_data()
  , elem_node_ids( arg_mesh.elem_node_ids )
  , node_coords(   arg_mesh.node_coords )
  , nodal_values(   arg_nodal_values )
  , element_matrices( arg_element_matrices )
  , element_vectors( arg_element_vectors )
  , coeff_K( arg_coeff_K )
  {
    const size_t elem_count = arg_mesh.elem_node_ids.dimension(0);

    parallel_for( elem_count , *this );
  }

  //------------------------------------

  static const unsigned FLOPS_transform_gradients =
     /* Jacobian */           FunctionCount * TensorDim * 2 +
     /* Inverse jacobian */   TensorDim * 6 + 6 +
     /* Gradient transform */ FunctionCount * 15 ;

  KOKKOSARRAY_INLINE_FUNCTION
  float transform_gradients(
    const float grad[][ FunctionCount ] , // Gradient of bases master element
    const double x[] ,
    const double y[] ,
    const double z[] ,
    float dpsidx[] ,
    float dpsidy[] ,
    float dpsidz[] ) const
  {
    enum { j11 = 0 , j12 = 1 , j13 = 2 ,
           j21 = 3 , j22 = 4 , j23 = 5 ,
           j31 = 6 , j32 = 7 , j33 = 8 };

    // Jacobian accumulation:

    double J[ TensorDim ] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

    for( unsigned i = 0; i < FunctionCount ; ++i ) {
      const double x1 = x[i] ;
      const double x2 = y[i] ;
      const double x3 = z[i] ;

      const float g1 = grad[0][i] ;
      const float g2 = grad[1][i] ;
      const float g3 = grad[2][i] ;

      J[j11] += g1 * x1 ;
      J[j12] += g1 * x2 ;
      J[j13] += g1 * x3 ;

      J[j21] += g2 * x1 ;
      J[j22] += g2 * x2 ;
      J[j23] += g2 * x3 ;

      J[j31] += g3 * x1 ;
      J[j32] += g3 * x2 ;
      J[j33] += g3 * x3 ;
    }

    // Inverse jacobian:

    float invJ[ TensorDim ] = {
      ( J[j22] * J[j33] - J[j23] * J[j32] ) ,
      ( J[j13] * J[j32] - J[j12] * J[j33] ) ,
      ( J[j12] * J[j23] - J[j13] * J[j22] ) ,

      ( J[j23] * J[j31] - J[j21] * J[j33] ) ,
      ( J[j11] * J[j33] - J[j13] * J[j31] ) ,
      ( J[j13] * J[j21] - J[j11] * J[j23] ) ,

      ( J[j21] * J[j32] - J[j22] * J[j31] ) ,
      ( J[j12] * J[j31] - J[j11] * J[j32] ) ,
      ( J[j11] * J[j22] - J[j12] * J[j21] ) };

    const float detJ = J[j11] * invJ[j11] +
                       J[j21] * invJ[j12] +
                       J[j31] * invJ[j13] ;

    const float detJinv = 1.0 / detJ ;

    for ( unsigned i = 0 ; i < TensorDim ; ++i ) { invJ[i] *= detJinv ; }

    // Transform gradients:

    for( unsigned i = 0; i < FunctionCount ; ++i ) {
      const float g0 = grad[0][i];
      const float g1 = grad[1][i];
      const float g2 = grad[2][i];

      dpsidx[i] = g0 * invJ[j11] + g1 * invJ[j12] + g2 * invJ[j13];
      dpsidy[i] = g0 * invJ[j21] + g1 * invJ[j22] + g2 * invJ[j23];
      dpsidz[i] = g0 * invJ[j31] + g1 * invJ[j32] + g2 * invJ[j33];
    }

    return detJ ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  void contributeResidualJacobian(
    const float coeff_k ,
    const double dof_values[] ,
    const float dpsidx[] ,
    const float dpsidy[] ,
    const float dpsidz[] ,
    const float detJ ,
    const float integ_weight ,
    const float bases_vals[] ,
    double elem_res[] ,
    double elem_mat[][ FunctionCount ] ) const
  {
    double value_at_pt = 0 ;
    double gradx_at_pt = 0 ;
    double grady_at_pt = 0 ;
    double gradz_at_pt = 0 ;

    for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
      value_at_pt += dof_values[m] * bases_vals[m] ;
      gradx_at_pt += dof_values[m] * dpsidx[m] ;
      grady_at_pt += dof_values[m] * dpsidy[m] ;
      gradz_at_pt += dof_values[m] * dpsidz[m] ;
    }

    const scalar_type k_detJ_weight = coeff_k        * detJ * integ_weight ;
    const double res_val = value_at_pt * value_at_pt * detJ * integ_weight ;
    const double mat_val = 2.0 * value_at_pt         * detJ * integ_weight ; 

    // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$ 
    // $$ J_{i,j} = \frac{\partial R_i}{\partial T_j} = \int_{\Omega} k \nabla \phi_i \cdot \nabla \phi_j + 2 \phi_i \phi_j T d \Omega $$ 

    for ( unsigned m = 0; m < FunctionCount; m++) {
      double * const mat = elem_mat[m] ;
      const float bases_val_m = bases_vals[m];
      const float dpsidx_m    = dpsidx[m] ;
      const float dpsidy_m    = dpsidy[m] ;
      const float dpsidz_m    = dpsidz[m] ;

      elem_res[m] += k_detJ_weight * ( dpsidx_m * gradx_at_pt +
                                       dpsidy_m * grady_at_pt +
                                       dpsidz_m * gradz_at_pt ) +
                     res_val * bases_val_m ;

      for( unsigned n = 0; n < FunctionCount; n++) {
      
        mat[n] += k_detJ_weight * ( dpsidx_m * dpsidx[n] + 
                                    dpsidy_m * dpsidy[n] +
                                    dpsidz_m * dpsidz[n] ) +
                  mat_val * bases_val_m * bases_vals[n]; 
      }
    }
  }

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    // Gather nodal coordinates and solution vector:

    double x[ FunctionCount ] ;
    double y[ FunctionCount ] ;
    double z[ FunctionCount ] ;
    double val[ FunctionCount ] ;

    for ( unsigned i = 0 ; i < ElementNodeCount ; ++i ) {
      const unsigned node_index = elem_node_ids( ielem , i );

      x[i] = node_coords( node_index , 0 );
      y[i] = node_coords( node_index , 1 );
      z[i] = node_coords( node_index , 2 );

      val[i] = nodal_values( node_index );
    }

    double elem_vec[ FunctionCount ] ;
    double elem_mat[ FunctionCount ][ FunctionCount ] ;

    for( unsigned i = 0; i < FunctionCount ; i++ ) {
      elem_vec[i] = 0 ;
      for( unsigned j = 0; j < FunctionCount ; j++){
        elem_mat[i][j] = 0 ;
      }
    }

    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {
      float dpsidx[ FunctionCount ] ;
      float dpsidy[ FunctionCount ] ;
      float dpsidz[ FunctionCount ] ;

      const float detJ =
        transform_gradients( elem_data.gradients[i] , x , y , z ,
                             dpsidx , dpsidy , dpsidz );

      contributeResidualJacobian( coeff_K ,
                                  val , dpsidx , dpsidy , dpsidz ,
                                  detJ ,
                                  elem_data.weights[i] ,
                                  elem_data.values[i] , 
                                  elem_vec , elem_mat );
    }

    for( unsigned i = 0; i < FunctionCount ; i++){
      element_vectors(ielem, i) = elem_vec[i] ;
      for( unsigned j = 0; j < FunctionCount ; j++){
        element_matrices(ielem, i, j) = elem_mat[i][j] ;
      }
    }
  }

}; /* ElementComputation */

//----------------------------------------------------------------------------

} /* namespace NonLinear */
} /* namespace HybridFEM */

