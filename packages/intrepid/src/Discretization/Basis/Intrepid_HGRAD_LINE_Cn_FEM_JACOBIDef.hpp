// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_TRI_C1_FEM_Def.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(grad) functions on TRI.
    \author Created by P. Bochev and D. Ridzal.
 */

namespace Intrepid {


template<class Scalar, class ArrayScalar>
Basis_HGRAD_LINE_Cn_FEM_JACOBI<Scalar,ArrayScalar>::Basis_HGRAD_LINE_Cn_FEM_JACOBI(int order, Scalar alpha, Scalar beta) {
    this -> basisCardinality_  = order+1;
    this -> basisDegree_       = order;    
    this -> jacobiAlpha_       = alpha;    
    this -> jacobiBeta_        = beta;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    this -> basisType_         = BASIS_FEM_HIERARCHICAL;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
}



template<class Scalar, class ArrayScalar> 
void Basis_HGRAD_LINE_Cn_FEM_JACOBI<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                                    const ArrayScalar &  inputPoints,
                                                                    const EOperator      operatorType) const {
  
  // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif
  
  // Number of evaluation points = dimension 0 of inputPoints
  int numPoints = inputPoints.dimension(0);  
  
  Teuchos::Array<Scalar> tmpPoints(numPoints);
  Teuchos::Array<Scalar> jacobiPolyAtPoints(numPoints);

  // Copy inputPoints into tmpPoints, to prepare for call to jacobfd
  for (int i=0; i<numPoints; i++) {
    tmpPoints[i] = inputPoints(i, 0);
  }

  switch (operatorType) {
    
    case OPERATOR_VALUE: {
      for (int ord = 0; ord < this -> basisCardinality_; ord++) {
        IntrepidPolylib::jacobfd(numPoints, &tmpPoints[0], &jacobiPolyAtPoints[0], (Scalar*)0, ord, jacobiAlpha_, jacobiBeta_);
        for (int pt = 0; pt < numPoints; pt++) {
          // outputValues is a rank-2 array with dimensions (basisCardinality_, numPoints)
          outputValues(ord, pt) = jacobiPolyAtPoints[pt];
        }
      }
    }
    break;

    case OPERATOR_GRAD:
    case OPERATOR_CURL:
    case OPERATOR_DIV:
    case OPERATOR_D1: {
      for (int ord = 0; ord < this -> basisCardinality_; ord++) {
        IntrepidPolylib::jacobd(numPoints, &tmpPoints[0], &jacobiPolyAtPoints[0], ord, jacobiAlpha_, jacobiBeta_);
        for (int pt = 0; pt < numPoints; pt++) {
          // outputValues is a rank-2 array with dimensions (basisCardinality_, numPoints)
          outputValues(ord, pt, 0) = jacobiPolyAtPoints[pt];
        }
      }
    }
    break;

    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10: {
      TEST_FOR_EXCEPTION( 1, std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM_JACOBI): Higher-order (>1st) derivatives not implemented.");
    }
    break;

    default:
      TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM_JACOBI): Invalid operator type.");
  }
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Cn_FEM_JACOBI<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                                    const ArrayScalar &    inputPoints,
                                                                    const ArrayScalar &    cellVertices,
                                                                    const EOperator        operatorType) const {
  TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM_JACOBI): FEM Basis calling an FVD member function");
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Cn_FEM_JACOBI<Scalar, ArrayScalar>::setBasisParameters(int n, Scalar alpha, Scalar beta) {
  this -> basisCardinality_  = n+1;
  this -> basisDegree_       = n;
  this -> jacobiAlpha_       = alpha;
  this -> jacobiBeta_        = beta;
  this -> initializeTags();
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Cn_FEM_JACOBI<Scalar, ArrayScalar>::initializeTags() {

  // Basis-dependent initializations

  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // poisition in the tag, counting from 0, of the subcell dim
  int posScOrd = 1;        // poisition in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  FieldContainer<int> tags(this->basisCardinality_, 4);

  for (int i=0; i < this->basisCardinality_; i++) {
    tags(i, 0) = 1;                        // these are all "internal" i.e. "volume" DoFs
    tags(i, 1) = 0;                        // there is only one line
    tags(i, 2) = i;                        // local DoF id 
    tags(i, 3) = this->basisCardinality_;  // total number of DoFs 
  }

  // Basis-independent function, sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              &tags[0],
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
}

}// namespace Intrepid
