//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"


/**  \brief  Returns the nodes and weights for the integration \n
             on the interval [0,1] (dim = 1) \n
             on the triangle with vertices (0,0), (1,0), (0,1) (if dim = 2) \n
             on the tetrahedron with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1) (if dim = 3).
             
  \param  dim     [in]   - spatial dimension (dim = 1, 2)
  \param  order   [in]   - required degree of polynomials that integrate exactly
  \param  nodes   [out]  - Matrix in which the i-th row of nodes gives coordinates of the i-th quadrature node
  \param  weights [out]  - quadrature weights

  \return 0                if successful
*/
int quadrature(const int dim, const int order,
               Epetra_SerialDenseMatrix & nodes,
               Epetra_SerialDenseVector & weights)
{
  
  if (dim == 1) {

    // Gauss quadrature nodes and weights on the interval [0,1]

    if (order == 1) {
      nodes.Shape(1,1);
      nodes(0,0) = 0.5;
      weights.Size(1);
      weights(0) = 1.0;
    }
    else if (order == 2) {
      nodes.Shape(2,1);
      nodes(0,0) = (1.0-1.0/sqrt(3.0))/2.0;
      nodes(1,0) = (1.0+1.0/sqrt(3.0))/2.0;
      weights.Size(2);
      weights(0) = 0.5;
      weights(1) = 0.5;
    }
    else if (order == 3) {
      nodes.Shape(3,1);
      nodes(0,0) = (1.0-sqrt(3.0/5.0))/2.0;
      nodes(1,0) = 0.5;
      nodes(2,0) = (1.0+sqrt(3.0/5.0))/2.0;
      weights.Size(3);
      weights(0) = 5.0/18.0;
      weights(1) = 4.0/9.0;
      weights(2) = 5.0/18.0;
    }
    else {
      cout << "Quadrature for dim = " << dim << " and order = ";
      cout << order << " not available.\n";
      exit(-1);
    }

  }
  else if (dim == 2) {
    
    // Quadrature nodes and weights on the unit simplex with
    // vertices (0,0), (1,0), and (0,1).

    if (order == 1) {
      nodes.Shape(1,2);
      nodes(0,0) = 1.0/3.0; nodes (0,1) = 1.0/3.0;
      weights.Size(1);
      weights(0) = 0.5;
    }
    else if (order == 2) {
      nodes.Shape(3,2);
      nodes(0,0) = 1.0/6.0; nodes (0,1) = 1.0/6.0;
      nodes(1,0) = 2.0/3.0; nodes (1,1) = 1.0/6.0;
      nodes(2,0) = 1.0/6.0; nodes (2,1) = 2.0/3.0;
      weights.Size(3);
      weights(0) = 1.0/6.0;
      weights(1) = 1.0/6.0;
      weights(2) = 1.0/6.0;
    }
    else if (order == 3) {
      nodes.Shape(4,2);
      nodes(0,0) = 1.0/3.0; nodes (0,1) = 1.0/3.0;
      nodes(1,0) = 3.0/5.0; nodes (1,1) = 1.0/5.0;
      nodes(2,0) = 1.0/5.0; nodes (2,1) = 3.0/5.0;
      nodes(3,0) = 1.0/5.0; nodes (3,1) = 1.0/5.0;
      weights.Size(4);
      weights(0) = -9.0/32.0;
      weights(1) = 25.0/96.0;
      weights(2) = 25.0/96.0;
      weights(3) = 25.0/96.0;
    }
    else {
      cout << "Quadrature for dim = " << dim << " and order = ";
      cout << order << " not available.\n";
      exit(-1);
    }

  }
  else {
    cout << "Quadrature for dim = " << dim << " not available.\n";
    exit(-1);
  }

  return(0);
}
