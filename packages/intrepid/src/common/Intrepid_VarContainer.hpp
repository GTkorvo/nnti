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

/** \file   Intrepid_VarContainer.hpp
\brief  Header file for utility class to provide variable containers. VarContainer is used to 
store values of vector and scalar fields and their derivatives at a set of points.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_VARCONTAINER_HPP
#define INTREPID_VARCONTAINER_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_LexContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \class Intrepid::VarContainer
  \brief Implementation of a templated variable container to store and access values of an operator 
  applied to a <strong>set of fields</strong> (scalar, vector or tensor), evaluated at a 
  <strong>set of points</strong>  on a <\strong>single</strong> cell. 
  
  Provides more specialized storage option than a LexContainer by allowing access to values by point 
  index, function index, function coordinate, and enumeration of derivatives of order k. 
  Values can also be accessed by their enumeration value using an overloaded [] operator. 
  Uses LexContainer for actual storage.     
  */
  template<class Scalar>
  class VarContainer : public LexContainer<Scalar> {
  private:
    
    /** \brief Type of the field whose values will be stored
    */
    EField fieldType_;
    
    /** \brief Type of the operator that is applied to the field
    */
    EOperator operatorType_;
    
    /** \brief The space dimension in which this container is defined
    */
    int spaceDim_;
    
  public:
      
    /** \brief Default destructor.
    */
    ~VarContainer() {};
    
    
    /** \brief Default constructor: LexContainer with zero capacity, dummy field, operator & dim.
    */
    VarContainer() : LexContainer<Scalar>::LexContainer() {
      fieldType_    = FIELD_FORM_0;
      operatorType_ = OPERATOR_VALUE;
      spaceDim_     = 1;
    };
    
    
    /** \brief Copy constructor.
    */
    VarContainer(const VarContainer& right);
    
    
    /** \brief Define VarContainer to store values based on number of points, number of fields
               field type, applied operator type and space dimension. 
      
      \param numPoints       [in]        - number of evaluation points
      \param numFields       [in]        - number of fields that will be evaluated
      \param fieldType       [in]        - type of the field that will be evaluated
      \param operatorType    [in]        - type of the operator that will be applied to the field
      \param spaceDim        [in]        - dimension of the ambient space
    */
    VarContainer(const int         numPoints,
                 const int         numFields,
                 const EField      fieldType,
                 const EOperator   operatorType,
                 const int         spaceDim);
    
    
    /** \brief Returns number of evaluation points. Same as getIndexBound(0).
    */
    int getNumPoints() const;
    
    
    /** \brief Returns number of fields that are being evaluated. Same as getIndexBound(1).
    */
    int getNumFields() const;
    
    
    /** \brief Returns type of the fields that are being evaluated.
    */
    EField getFieldType() const;
    
    
    /** \brief Returns type of the operator that is applied to the fields.
    */
    EOperator getOperatorType() const;
    
    
    /** \brief Returns dimension of the ambient space where fields are defined.
    */
    int getSpaceDim() const;
    
        
    /** \brief Returns the enumeration of a value based on its indices. The number, type  
      and bounds of the indices depend on the number of points, number of fields, the rank of the
      field and the rank of the operator applied to the field. They are summarized below:
      \verbatim
      |--------------------|-------------------|-------------------|------------------------|
      |operator/field rank |       rank 0      | rank 1 2D/3D      | rank 2 2D/3D           |
      |--------------------|-------------------|-------------------|------------------------|
      |       VALUE        | [P][F]            | [P][F][D]         | [P][F][D][D]           |
      |--------------------|-------------------|-------------------|------------------------|
      |     GRAD, D1       | [P][F][D]         | [P][F][D][D]      | [P][F][D][D][D]        |
      |--------------------|-------------------|-------------------|------------------------|
      |        CURL        | [P][F][D](undef3D)| [P][F]/[P][F][D]  | [P][F][D]/[P][F][D][D] |
      |--------------------|-------------------|-------------------|------------------------|
      |        DIV         | [P][F][D](only 1D)| [P][F]            | [P][F][D]              |
      |--------------------|-------------------|-------------------|------------------------|
      |    D1,D2,..,D10    | [P][F][K]         | [P][F][D][K]      | [P][F][D][D][K]        |
      |--------------------|-------------------|-------------------|------------------------|

      Legend:
        P -> point index            range 0 <= P < numPoints
        F -> field index            range 0 <= F < numFields
        D -> field component index  range 0 <= D < spaceDim
        K -> enumeration of Dk      range 0 <= K < DkCardinality
      \endverbatim
      
      \note If operator is one of the k-th derivative operators Dk, the enumeration index K of
      a particular derivative in the set of all kth order derivatives can be obtained from the
      multiplicities of dx, dy and dz using global function getDkEnumeration.
     
      \note The rank of the VarContainer equals the number of indices needed to tag its values 
      and ranges from 2 to 5. The input array must contain exactly the same number of 
      indices, and each index value must be within its admissible bounds as defined above.      
      
      \param multiIndex      [in]     - a set of indices whose enumeration we want
    */
    int getEnumeration(const Teuchos::Array<int>& multiIndex) const;
    
    
    /** \brief Resets VarContainer to store values of an operator applied to a set of numFields
      fields, evaluated at numPoints points.
      
      \param numPoints       [in]        - number of evaluation points
      \param numFields       [in]        - number of fields that will be evaluated
      \param fieldType       [in]        - type of the field whose basis will be evaluated
      \param operatorType    [in]        - type of the operator that will be applied to the basis
      \param spaceDim        [in]        - dimension of the ambient space
    */
    void reset(const int       numPoints,
               const int       numFields,
               const EField    fieldType,
               const EOperator operatorType,
               const int       spaceDim);
    
    
    /** \brief Assignment operator <var>*this = right</var>.
    */
    VarContainer& operator  = (const VarContainer& right);
    
    
    /** \brief Prints basic facts about the container: size, rank, index bounds, number of
      points, number of fields, operator type and space dimension
    */
    void print(std::ostream& os) const;

    
    
  }; // end class VarContainer
  
  //===========================================================================//
  //                                                                           //
  //                Function declarations related to VarContainer              //
  //                                                                           //
  //===========================================================================//
  
  /** \relates VarContainer
  Outputs a formated stream with VarContainer data. For debugging purposes.
  */
  template<class Scalar>
  std::ostream& operator <<  (std::ostream& os, const VarContainer<Scalar>& container);
  
  
  
  /** \relates VarContainer
    Field rank is defined as the number of indices needed to specify the value. Intrepid provides
    three basic types of fields: scalar, vectors and tensors with ranks 0,1, and 2 respectively. 
    The scalar field types are FIELD_FORM_0 and FIELD _FORM_3. The vector field types are 
    FIELD_FORM_1, FIELD_FORM_2 and FIELD_VECTOR. The tensor field type is FIELD_TENSOR.
    
    \param fieldType         [in]     - field type whose rank we want to know
    */
  int getFieldRank(const EField fieldType); 
  
  
  
  /**\relates varContainer
    Returns rank of an operator. When an operator acts on a field with a certain rank, the result
    is a field with another rank. Operator rank is defined as the difference between the ranks of
    the output field and the input field:
    
                  Rank(OPERATOR) = Rank(OPERATOR(FIELD)) - Rank(FIELD)
    
    This function is used by VarContainer class to determine how many indices are needed to store
    the values of OPERATOR(FIELD). By knowing the rank of the operator and the rank of the field we
    can figure out the rank of the result:
    
                  Rank(OPERATOR(FIELD)) = Rank(OPERATOR) + Rank(FIELD)
    
    Operator ranks are defined in the following table: (X denotes undefined, below slash means 3D).
    By default, in 1D any operator other than VALUE has rank 1, i.e., GRAD, CURL and DIV default to
    d/dx and Dk are the higher-order derivatives d^k/dx^k. Only scalar functions are allowed in 1D.
    
    \verbatim
    |--------|------|--------------------------|----------|----------|----------|----------|
    | field  | rank |  field type ENUM         |   VALUE  | GRAD, Dk |   CURL   |    DIV   |
    |--------|------|--------------------------|----------|----------|----------|----------|
    | scalar |   0  |  FORM_0, FORM_3          |     0    |     1    | 3-dim/X  |     X    |
    | vector |   1  |  FORM_1, FORM_2, VECTOR  |     0    |     1    | dim - 3  |    -1    |
    | tensor |   2  |  TENSOR                  |     0    |     1    | dim - 3  |    -1    |
    |--------|------|--------------------------|----------|----------|----------|----------|
    |   1D   |   0  |  FORM_0, FORM_3 only     |     0    |     1    |     1    |     1    |
    |--------|------|--------------------------|----------|----------|----------|----------|
    \endverbatim
      
      \param operatorType     [in]    - the operator acting on the field
      \param fieldRank        [in]    - rank of the field being acted upon (use getFieldRank)
      \param spaceDim         [in]    - dimension of the space where all happens (1,2 or 3)
    */
  int getOperatorRank(const EOperator operatorType,
                      const int       fieldRank,
                      const int       spaceDim);

  
  
  /**\relates VarContainer 
    Returns order of an operator: ranges from 0 (for OPERATOR_VALUE) to 10 (OPERATOR_D10)
    
    \param operatorType       [in]    - type of the operator whose order we want to know
    */
  int getOperatorOrder(const EOperator operatorType);
  
  
  
  /** \relates VarContainer
    Returns enumeration of a partial derivative of order k based on multiplicities of the partial
    derivatives dx, dy and dz. Enumeration follows the lexicographical order of the partial 
    derivatives multiplicities. For example, the 10 derivatives of order 3 in 3D are enumerated as:
    
    D3 = {(3,0,0), (2,1,0), (2,0,1), (1,2,0), (1,1,1), (1,0,2), (0,3,0), (0,2,1), (0,1,2), (0,0,3)}
  
    Enumeration formula for lexicographically ordered partial derivatives of order k is given by
          
    <table>   
      <tr> <td>\f$i(xMult)               = 0\f$</td>              <td>in 1D (only 1 derivative)</td> </tr>
      <tr> <td>\f$i(xMult, yMult)        = yMult\f$</td>                              <td>in 2D</td> </tr>
      <tr> <td>\f$i(xMult, yMult, zMult) = zMult + \displaystyle\sum_{r = 0}^{k - xMult} r\f$</td> <td>in 3D</td> </tr>
    </table>
  
    where the order k of Dk is implicitly defined by xMult + yMult + zMult. Space dimension is
    implicitly defined by the default values of the multiplicities of y and z derivatives.
  
    \param xMult            [in]    - multiplicity of the partial derivative in x
    \param yMult            [in]    - multiplicity of the partial derivative in y (default = -1)
    \param zMult            [in]    - multiplicity of the partial derivative in z (default = -1)
  */
  int getDkEnumeration(const int xMult,
                       const int yMult = -1,
                       const int zMult = -1);
  
  
  
  /** \relates VarContainer
    Returns multiplicities of the partials in x, y, and z based on the enumeration of the partial
    derivative, its order and the space dimension. Inverse of the getDkEnumeration() method.
    
    \param partialMult      [out]    - array with partial derivative multiplicities
    \param derivativeEnum   [in]     - enumeration of the partial derivative
    \param derivativeOrder  [in]     - order of Dk
    \param spaceDim         [in]     - space dimension
  */
  void getDkMultiplicities(Teuchos::Array<int>& partialMult,
                           const int derivativeEnum,
                           const int derivativeOrder,
                           const int spaceDim);
  
    
  
  /** \relates VarContainer
    Returns cardinality of Dk, i.e., the number of all derivatives of order k. The set of all 
    partial derivatives of order k is isomorphic to the set of all multisets of cardinality k with 
    elements taken from the sets {x}, {x,y}, and {x,y,z} in 1D, 2D, and 3D respectively. 
    For example, the partial derivative
    
    \f$\displaystyle D\{1,2,1\}f = \frac{d^4 f}{dx dy^2 dz}\f$   maps to multiset
    \f$\{x, y, z\}\f$ with multiplicities \f$\{1,2,1\}\f$ 
    
    The number of all such multisets is given by the binomial coefficient

    \f$
      \begin{pmatrix}
        spaceDim + k - 1 \\
        spaceDim - 1
      \end{pmatrix}
    \f$
    
    1D: cardinality = 1\n
    2D: cardinality = k + 1\n
    3D: cardinality = (k + 1)*(k + 2)/2
    
    \param derivativeOrder  [in]     - order of Dk
    \param spaceDim         [in]     - space dimension
    */
  int getDkCardinality(const int derivativeOrder,
                       const int spaceDim);
  
  //===========================================================================//
  //                                                                           //
  //              End Function declarations related to VarContainer            //
  //                                                                           //
  //===========================================================================//
  
} // end namespace Intrepid

// include templated definitions
#include <Intrepid_VarContainerDef.hpp>

#endif
