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
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//                    FIAT (fiat-dev@fenics.org) (must join mailing list first, see www.fenics.org)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_F2_TET_I5_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 4 for 2-forms on TET cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F2_TET_I5_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 2 , 0 , 0 , 15 , 2 , 0 , 1 , 15 , 2 , 0 , 2 , 15 , 2 , 0 , 3 , 15 , 2 , 0 , 4 , 15 , 2 , 0 , 5 , 15 , 2 , 0 , 6 , 15 , 2 , 0 , 7 , 15 , 2 , 0 , 8 , 15 , 2 , 0 , 9 , 15 , 2 , 0 , 10 , 15 , 2 , 0 , 11 , 15 , 2 , 0 , 12 , 15 , 2 , 0 , 13 , 15 , 2 , 0 , 14 , 15 , 2 , 1 , 0 , 15 , 2 , 1 , 1 , 15 , 2 , 1 , 2 , 15 , 2 , 1 , 3 , 15 , 2 , 1 , 4 , 15 , 2 , 1 , 5 , 15 , 2 , 1 , 6 , 15 , 2 , 1 , 7 , 15 , 2 , 1 , 8 , 15 , 2 , 1 , 9 , 15 , 2 , 1 , 10 , 15 , 2 , 1 , 11 , 15 , 2 , 1 , 12 , 15 , 2 , 1 , 13 , 15 , 2 , 1 , 14 , 15 , 2 , 2 , 0 , 15 , 2 , 2 , 1 , 15 , 2 , 2 , 2 , 15 , 2 , 2 , 3 , 15 , 2 , 2 , 4 , 15 , 2 , 2 , 5 , 15 , 2 , 2 , 6 , 15 , 2 , 2 , 7 , 15 , 2 , 2 , 8 , 15 , 2 , 2 , 9 , 15 , 2 , 2 , 10 , 15 , 2 , 2 , 11 , 15 , 2 , 2 , 12 , 15 , 2 , 2 , 13 , 15 , 2 , 2 , 14 , 15 , 2 , 3 , 0 , 15 , 2 , 3 , 1 , 15 , 2 , 3 , 2 , 15 , 2 , 3 , 3 , 15 , 2 , 3 , 4 , 15 , 2 , 3 , 5 , 15 , 2 , 3 , 6 , 15 , 2 , 3 , 7 , 15 , 2 , 3 , 8 , 15 , 2 , 3 , 9 , 15 , 2 , 3 , 10 , 15 , 2 , 3 , 11 , 15 , 2 , 3 , 12 , 15 , 2 , 3 , 13 , 15 , 2 , 3 , 14 , 15 , 3 , 0 , 0 , 60 , 3 , 0 , 1 , 60 , 3 , 0 , 2 , 60 , 3 , 0 , 3 , 60 , 3 , 0 , 4 , 60 , 3 , 0 , 5 , 60 , 3 , 0 , 6 , 60 , 3 , 0 , 7 , 60 , 3 , 0 , 8 , 60 , 3 , 0 , 9 , 60 , 3 , 0 , 10 , 60 , 3 , 0 , 11 , 60 , 3 , 0 , 12 , 60 , 3 , 0 , 13 , 60 , 3 , 0 , 14 , 60 , 3 , 0 , 15 , 60 , 3 , 0 , 16 , 60 , 3 , 0 , 17 , 60 , 3 , 0 , 18 , 60 , 3 , 0 , 19 , 60 , 3 , 0 , 20 , 60 , 3 , 0 , 21 , 60 , 3 , 0 , 22 , 60 , 3 , 0 , 23 , 60 , 3 , 0 , 24 , 60 , 3 , 0 , 25 , 60 , 3 , 0 , 26 , 60 , 3 , 0 , 27 , 60 , 3 , 0 , 28 , 60 , 3 , 0 , 29 , 60 , 3 , 0 , 30 , 60 , 3 , 0 , 31 , 60 , 3 , 0 , 32 , 60 , 3 , 0 , 33 , 60 , 3 , 0 , 34 , 60 , 3 , 0 , 35 , 60 , 3 , 0 , 36 , 60 , 3 , 0 , 37 , 60 , 3 , 0 , 38 , 60 , 3 , 0 , 39 , 60 , 3 , 0 , 40 , 60 , 3 , 0 , 41 , 60 , 3 , 0 , 42 , 60 , 3 , 0 , 43 , 60 , 3 , 0 , 44 , 60 , 3 , 0 , 45 , 60 , 3 , 0 , 46 , 60 , 3 , 0 , 47 , 60 , 3 , 0 , 48 , 60 , 3 , 0 , 49 , 60 , 3 , 0 , 50 , 60 , 3 , 0 , 51 , 60 , 3 , 0 , 52 , 60 , 3 , 0 , 53 , 60 , 3 , 0 , 54 , 60 , 3 , 0 , 55 , 60 , 3 , 0 , 56 , 60 , 3 , 0 , 57 , 60 , 3 , 0 , 58 , 60 , 3 , 0 , 59 , 60 };
  
  // Basis-independent function sets tag and enum data in the static arrays:
  Intrepid::setEnumTagData(tagToEnum_,
                           enumToTag_,
                           tags,
                           numDof_,
                           tagSize,
                           posScDim,
                           posScId,
                           posBfId);
}

template<class Scalar> 
void Basis_F2_TET_I5_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 4 (I4) has 120 basis functions on a tetrahedron that are 2-forms in 3D
  int    numFields = 120;
  EField fieldType = FIELD_FORM_2;
  int    spaceDim  = 3;

  // temporaries
  int countPt  = 0;               // point counter
  Teuchos::Array<int> indexV(3);  // multi-index for values
  Teuchos::Array<int> indexD(2);  // multi-index for divergences

  // Shape the FieldContainer for the output values using these values:
  outputValues.resize(numPoints,
                     numFields,
                     fieldType,
                     operatorType,
                     spaceDim);

#ifdef HAVE_INTREPID_DEBUG
  for (countPt=0; countPt<numPoints; countPt++) {
    // Verify that all points are inside the TET reference cell
    TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_TET, inputPoints[countPt]),
                        std::invalid_argument,
                        ">>> ERROR (Basis_F2_TET_I5_FEM_FIAT): Evaluation point is outside the TET reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(56,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(120,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,5,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<120;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<120;i++) {
          outputValues(countPt,i,1) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm2_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<120;i++) {
          outputValues(countPt,i,2) = result(i,countPt);
        }
      }
    }
    break;

    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TET_I5_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TET_I5_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_DIV:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(56,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(56,56);
      Teuchos::SerialDenseMatrix<int,Scalar> div(120,numPoints);
      div.putScalar(0.0);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,5,inputPoints,expansions);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,1.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,1.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm2_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<120;i++) {
          outputValues(countPt,i) = div(i,countPt);
        }
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
                            (operatorType != OPERATOR_GRAD)  &&
                            (operatorType != OPERATOR_CURL)  &&
                            (operatorType != OPERATOR_DIV)   &&
                            (operatorType != OPERATOR_D1)    &&
                            (operatorType != OPERATOR_D2)    &&
                            (operatorType != OPERATOR_D3)    &&
                            (operatorType != OPERATOR_D4)    &&
                            (operatorType != OPERATOR_D5)    &&
                            (operatorType != OPERATOR_D6)    &&
                            (operatorType != OPERATOR_D7)    &&
                            (operatorType != OPERATOR_D8)    &&
                            (operatorType != OPERATOR_D9)    &&
                            (operatorType != OPERATOR_D10) ),
                          std::invalid_argument,
                          ">>> ERROR (Basis_F2_TET_I5_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F2_TET_I5_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F2_TET_I5_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F2_TET_I5_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F2_TET_I5_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F2_TET_I5_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F2_TET_I5_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F2_TET_I5_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TET;
}



template<class Scalar>
EBasis Basis_F2_TET_I5_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F2_TET_I5_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F2_TET_I5_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

