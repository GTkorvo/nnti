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

/** \file   Intrepid_LexContainer.hpp
    \brief  Header file for utility class to provide lexicographical containers.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_LEXCONTAINER_HPP
#define INTREPID_LEXCONTAINER_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \class Intrepid::LexContainer
  \brief Implementation of a templated lexicographical container for a multi-indexed quantity. 
  LexContainer object stores a multi-indexed value using the lexicographical index ordering: the
  rightmost index changes first and the leftmost index changes last. LexContainer can be viewed
  as a dynamic multidimensional array whose values can be accessed in two ways: by their
  multi-index or by their address, using an overloaded [] operator. The address gives the enumeration
  of the multi-indexed value with respect to the lexicographical order. 
  The number of indices is unlimited.
  */
template<class Scalar>
class LexContainer {
  private:
    
    /** \brief Array to store the multi-indexed quantity 
      */
    Teuchos::Array<Scalar> data_;
    
    /**\brief Array to store upper bounds for the multi-indices. Admissible range for the k-th index is
       0 <= i_k < indexRange_[k]. Size of this array definees the rank of the multi-indexed quantity, 
       i.e., number of its indices.  
      */
    Teuchos::Array<int> indexRange_;
    
  public:
      
    /** \brief Default destructor.
      */
    ~LexContainer() {};
    
  
    /** \brief Default constructor.
      */
    LexContainer() {
      data_.resize(0);
      indexRange_.resize(0);
    };
    
    
    /** \brief Copy constructor.
    */
    LexContainer(const LexContainer& right);

    
    /** \brief Create empty LexContainer using specified index ranges. The size of the argument
      implicitly defines the rank of the LexContainer. Capacity of the LexContainer is implicitly
      defined by the specified upper index range values.
      
      \param indexRange[in]           - array with upper index bounds
      */
    LexContainer(const Teuchos::Array<int>& indexRange);
    
    
    /** \brief Create LexContainer from Teuchos::Arrays of data and indices. If DEBUG mode is enabled,
      will check if size of data array matches the capacity required by upper index ranges.
      
      \param indexRange[in]           - array with upper index bounds (its size implicitly defines 
                                        the rank of the new LexContainer object)
      \param data[in]                 - array with container values
      */
    LexContainer(const Teuchos::Array<int>&    indexRange,
                 const Teuchos::Array<Scalar>& data);
    
        
    /** \brief Return rank of the LexContainer = number of indices used to tag the multi-indexed value
    */
    int getRank() const;

    
    /** \brief Compute size of the LexContainer using index ranges stored in the object.
      */
    int getSize() const;
    
    
    /** \brief Returns address of value based on its multi-index (in Teuchos::Array format). In DEBUG
      mode checks if number of multi-indices matches the rank of the LexContainer and 
      whether each index is within its admissible range.
      */
    int getAddress(Teuchos::Array<int> multiIndex) const;
    
    
    /** \brief Returns address of value based on its multi-index (in an int array format). 
      \warning Method does not check if number of multi-indices matches rank of the LexContainer, it
      will use as many index values as required by the current rank of the LexContainer.
      */
    int getAddress(int* multiIndexPtr) const;
      
    
    /** \brief Returns multi-index corresponding to a LexContainer address
      
      \param multiIndex   [out]       - array containg multi-index of the specified address
      \param valueAddress [in]        - address of the value in the LexContainer
     */
    void getMultiIndex(Teuchos::Array<int>& multiIndex,
                       const int valueAddress) const;
    
    
    /** \brief Retrieve value by its multi-index. In DEBUG mode checks if number of multi-indices 
      matches the rank of the LexContainer and  whether each index is within its admissible range. 
      To retrieve value by its address use overloaded [] operator.
      
      \param multiIndex [in]            - array containing multi-index of the desired value
     */
    Scalar getValue(const Teuchos::Array<int>& multiIndex) const;

    
    /** \brief Resets LexContainer to trivial container (one with rank = 0 and size = 0)
      */
    void emptyContainer();
    
    
    /** \brief Resets LexContainer to unit length and stores a single zero. Default for multi-indexed
      quantity whose values are all zero (such as for example, 3rd derivatives of a linear basis)
      */
    void storeZero(); 
    
    
    
    /** \brief General purpose resize method: resizes LexContainer based on a new array of index 
      ranges. The size of the input array implicitly defines the rank of the container. 
      The size of the container is computed from the new index ranges in the array argument.
      
      \param newIndexRange[in]             - new upper values for index ranges
      */
    void resizeContainer(const Teuchos::Array<int>& newIndexRange);
    
  
    
    /** \brief Specialized resize method, useful for resizing containers to store a multi-indexed
      quantity which represents the values of an operator applied to a set of scalar, vector or
      tensor fields, and evaluated at a set of points. Such quantities have three distinct groups
      of indices:
      group 1:  2 indices for number of points and number of fields
      group 2:  0,1, or 2 indices for scalar, vector or tensor field, respectively
      group 3:  0,1,2,.., or 10 indices: equal to the order of total derivative operator Dk
      The multi-index of such quantity looks like this: {np,nf; f0,..,f_{rank-1}; d0,...,d_{ord-1} }
      
      \param numPoints       [in]        - number of evaluation points
      \param numFields       [in]        - number of fields that will be evaluated
      \param fieldRank       [in]        - rank of the field: 0=scalar, 1=vector, 2=tensor
      \param operatorOrd     [in]        - order of the specified operator (0 - 10)
      \param spaceDim        [in]        - dimension of the ambient space
     */
    void resizeContainer(const int numPoints,
                         const int numFields,
                         const int fieldRank,
                         const int operatorOrd,
                         const int spaceDim);
    
    
    
    /** \brief Shapes LexContainer based on fieldType and operatorType arguments. This method is
      a wrapper for the specialized resizeContainer method, that translates its arguments into 
      a set of arguments for the specialized resizeContainer method.
      
      \param numPoints       [in]        - number of evaluation points
      \param numFields       [in]        - number of fields that will be evaluated
      \param fieldType       [in]        - type of the field whose basis will be evaluated
      \param operatorType    [in]        - type of the operator that will be applied to the basis
      \param spaceDim        [in]        - dimension of the ambient space
      */
    void shapeContainer(const int       numPoints,
                        const int       numFields,
                        const EField    fieldType,
                        const EOperator operatorType,
                        const int       spaceDim);
    
    
    /** \brief Assign value by its multi-index. Does not change rank or index ranges of the LexContainer.
      In DEBUG mode checks if number of multi-indices matches the rank of the LexContainer and  
      whether each index is within its admissible range.
      
      \param dataValue [in]             - value to be assigned
      \param multiIndex [in]            - multi-index of the value
      */
    void setValue(const Scalar               dataValue,
                  const Teuchos::Array<int>& multiIndex);
    
    
    /**\brief Assign values from Teuchos::Array without changing rank and index ranges of LexContainer. 
      In DEBUG mode checks if size of input argument matches capacity of LexContainer.
      
      \param dataArray[in]               - new values
      */
    void setValues(const Teuchos::Array<Scalar>& dataArray);
    
    
    /** \brief Assign values from array without changing rank and index ranges of LexContainer.
       
      \param dataPtr[in]                - pointer to array of data values
      \param dataSize[in]               - number of values in the array
      */
    void setValues(const Scalar* dataPtr,
                   const int     dataSize);
    
     
    /** \brief   Overloaded [] operator. returns value based on its address
    */
    const Scalar & operator [] (const int address) const;
    
    /** \brief Assignment operator <var>*this = right</var>.
    */
    LexContainer& operator  = (const LexContainer& right);
    
}; // end class LexContainer
  
//===========================================================================//
//                                                                           //
//                Function declarations related to LexContainer              //
//                                                                           //
//===========================================================================//

/** \relates LexContainer
    Outputs a formated stream with LexContainer data. For debugging purposes.
*/
template<class Scalar>
  std::ostream& operator << (std::ostream& os, const LexContainer<Scalar>& container);



/** \relates LexContainer
    Returns rank of the field: 0, 1, or 2 for scalar, vector, or tensor fields
*/
int getFieldRank(const EField fieldType); 



/**\relates LexContainer 
   Returns order of an operator: ranges from 0 (for OPERATOR_VALUE) to 10 (OPERATOR_D10)
*/
int getOperatorOrder(const EOperator operatorType);

//===========================================================================//
//                                                                           //
//              End Function declarations related to LexContainer            //
//                                                                           //
//===========================================================================//

  
} // end namespace Intrepid

// include templated definitions
#include <Intrepid_LexContainerDef.hpp>

#endif
