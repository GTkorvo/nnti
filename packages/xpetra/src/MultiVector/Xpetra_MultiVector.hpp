// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_MULTIVECTOR_HPP
#define XPETRA_MULTIVECTOR_HPP

/* this file is automatically generated - do not edit (see script/interfaces.py) */

#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DistObject.hpp"
#include "Xpetra_Map.hpp"

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Vector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class Vector;
#endif

  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MultiVector
    : public DistObject< Scalar, LocalOrdinal, GlobalOrdinal, Node >
  {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Destructor.
    virtual ~MultiVector() { }

   //@}

    //! @name Post-construction modification routines
    //@{

    //! Replace value, using global (row) index.
    virtual void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value)= 0;

    //! Add value to existing value, using global (row) index.
    virtual void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value)= 0;

    //! Replace value, using local (row) index.
    virtual void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value)= 0;

    //! Add value to existing value, using local (row) index.
    virtual void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value)= 0;

    //! Set all values in the multivector with the given value.
    virtual void putScalar(const Scalar &value)= 0;

    //@}

    //! @name Data Copy and View get methods
    //@{

    //! Const view of the local values in a particular vector of this multivector.
    virtual Teuchos::ArrayRCP< const Scalar > getData(size_t j) const = 0;

    //! View of the local values in a particular vector of this multivector.
    virtual Teuchos::ArrayRCP< Scalar > getDataNonConst(size_t j)= 0;

    //@}

    //! @name Mathematical methods
    //@{

    //! Compute dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i]).
    virtual void dot(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Teuchos::ArrayView< Scalar > &dots) const = 0;

    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this).
    virtual void abs(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A)= 0;

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    virtual void reciprocal(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A)= 0;

    //! Scale the current values of a multi-vector, this = alpha*this.
    virtual void scale(const Scalar &alpha)= 0;

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    virtual void update(const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Scalar &beta)= 0;

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    virtual void update(const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Scalar &beta, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, const Scalar &gamma)= 0;

    //! Compute 1-norm of each vector in multi-vector.
    virtual void norm1(const Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const = 0;

    //! Compute 2-norm of each vector in multi-vector.
    virtual void norm2(const Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const = 0;

    //! Compute Inf-norm of each vector in multi-vector.
    virtual void normInf(const Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const = 0;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    virtual void normWeighted(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &weights, const Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const = 0;

    //! Compute mean (average) value of each vector in multi-vector.
    virtual void meanValue(const Teuchos::ArrayView< Scalar > &means) const = 0;

    //! Matrix-matrix multiplication: this = beta*this + alpha*op(A)*op(B).
    virtual void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, const Scalar &beta)= 0;

    //! Element-wise multiply of a Vector A with a MultiVector B.
    virtual void elementWiseMultiply(Scalar scalarAB, const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, Scalar scalarThis)= 0;

    //@}

    //! @name Attribute access functions
    //@{

    //! Number of columns in the multivector.
    virtual size_t getNumVectors() const = 0;

    //! Local number of rows on the calling process.
    virtual size_t getLocalLength() const = 0;

    //! Global number of rows in the multivector.
    virtual global_size_t getGlobalLength() const = 0;

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    virtual std::string description() const = 0;

    //! Print the object with the given verbosity level to a FancyOStream.
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

    //@}

    //! @name Xpetra specific
    //@{
 
    //! Set seed for Random function.
    virtual void setSeed(unsigned int seed)= 0;


    virtual void randomize(bool bUseXpetraImplementation = false)= 0;

    //! Set multi-vector values to random numbers. XPetra implementation
    virtual void Xpetra_randomize()
    {
        typedef Teuchos::ScalarTraits<Scalar> SCT;

        const size_t numVectors = getNumVectors();
        for (size_t i = 0; i < numVectors; i++)
        {
            Teuchos::ArrayRCP< Scalar > datai = getDataNonConst(i);

            const size_t myLength = getLocalLength();
            for(size_t j=0; j<myLength; j++)
            {
                datai[j] = SCT::random();
            }
        }
    }


    //@}

  }; // MultiVector class

} // Xpetra namespace

#define XPETRA_MULTIVECTOR_SHORT
#endif // XPETRA_MULTIVECTOR_HPP
