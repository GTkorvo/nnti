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
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_TRI_C1_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_TRI_Cn_FEM_ORTH class.
    \author Created by Robert Kirby
*/

#ifndef INTREPID_HGRAD_TRI_Cn_FEM_ORTHHPP
#define INTREPID_HGRAD_TRI_Cn_FEM_ORTHHPP

#include "Intrepid_Basis.hpp"
#include "Sacado.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_TRI_Cn_FEM_ORTH
    \brief  Implementation of the default H(grad)-compatible orthogonal basis (Dubiner) of
            arbitrary degree on triangle.
  
  \verbatim  
    \remarks
    \li     DefaultBasisFactory will select this class if the following parameters are specified:
  \verbatim
  |=======================|===================================|
  |  CellTopology         |  Triangle                         |
  |-----------------------|-----------------------------------|
  |  EFunctionSpace       |  FUNCTION_SPACE_HGRAD             |
  |-----------------------|-----------------------------------|
  |  EDiscreteSpace       |  DISCRETE_SPACE_COMPLETE          |
  |-----------------------|-----------------------------------|
  |  degree               |  1                                |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_HIERARCHICAL           |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

    \li   All degrees of freedom are considered to be internal (ie not assembled)
  \endverbatim
  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_TRI_Cn_FEM_ORTH: public Basis<Scalar, ArrayScalar> {
private:
  /** \brief  Initializes <var>tagToOrdinal_</var> and
      <var>ordinalToTag_</var> lookup arrays.   */
  void initializeTags();

public:
  
  /** \brief  Constructor.
   */
  Basis_HGRAD_TRI_Cn_FEM_ORTH( int degree );  
  
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Triangle</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Triangle</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec .

      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
      \param  operatorType      [in]  - the operator acting on the basis functions    
  */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const EOperator        operatorType) const;
  
  
  /**  \brief  FVD basis evaluation: invocation of this method throws an exception.
   */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const ArrayScalar &    cellVertices,
                 const EOperator        operatorType = OPERATOR_VALUE) const;

  template<unsigned deriv>
  static void getVals( ArrayScalar       &outputValues ,
		       const int deg ,
		       const ArrayScalar &inputPoints );
		    
};

template<typename Scalar,typename ArrayScalar, unsigned derivOrder>
class TabulatorTri
{
public:
  static void tabulate( ArrayScalar & outputValues ,
			const int deg ,
			const ArrayScalar &inputPoints );
};

template<typename Scalar,typename ArrayScalar>
class TabulatorTri<Scalar,ArrayScalar,0>
{
public:
  static void tabulate( ArrayScalar & outputValues ,
			const int deg ,
			const ArrayScalar &inputPoints );
    
};

template<typename Scalar,typename ArrayScalar>
class TabulatorTri<Scalar,ArrayScalar,1>
{
public:
  static void tabulate( ArrayScalar & outputValues ,
			const int deg ,
			const ArrayScalar &inputPoints );
    
};



}// namespace Intrepid

#include "Intrepid_HGRAD_TRI_Cn_FEM_ORTHDef.hpp"

#endif

