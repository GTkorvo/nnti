// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_MULTIVECTOR_DECL_HPP
#define TPETRA_MULTIVECTOR_DECL_HPP

#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CombineMode.hpp"

// TODO: add principal use case instructions for memory management interfaces (view/copy extraction)
// TODO: expand user-visible documentation 

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of MultiVectorData, needed to prevent circular inclusions
  template<class Scalar> class MultiVectorData;
  // forward declaration of Vector, needed to prevent circular inclusions
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal> class Vector;
#endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal>
  class MultiVector : public DistObject<Scalar,LocalOrdinal,GlobalOrdinal> {

    public:

    /** \name Convenient typedefs */
    //@{ 

    /*! Non-const pointer-like typedef. 
        In a debug build (<tt>--enable-teuchos-debug</tt>), this is an 
        ArrayRCP<Scalar>. In a non-debug build, it is a <tt>Scalar *</tt>. In either case, 
        the syntax is the same: pointer arithmetic and indexing are supported. */
    typedef typename Teuchos::ArrayView<Scalar>::iterator pointer;

    /*! Const pointer-like typedef. 
        In a debug build (<tt>--enable-teuchos-debug</tt>), this is an 
        ArrayRCP<const Scalar>. In a non-debug build, it is a <tt>const Scalar *</tt>. In either case, 
        the syntax is the same: pointer arithmetic and indexing are supported. */
    typedef typename Teuchos::ArrayView<const Scalar>::iterator const_pointer;

    /*! Non-const double pointer-like typedef. 
        In a debug build (<tt>--enable-teuchos-debug</tt>), this is an 
        ArrayRCP<const ArrayRPC<Scalar> >. In a non-debug build, it is a <tt>Scalar * const *</tt>. In either case, 
        the syntax is the same, as a two-dimensional data structure. */
    typedef typename Teuchos::ArrayView<const typename Teuchos::ArrayView<Scalar>::iterator>::iterator double_pointer;

    /*! Const double pointer-like typedef. 
        In a debug build (<tt>--enable-teuchos-debug</tt>), this is an 
        ArrayRCP<const ArrayRCP<const Scalar> >. In a non-debug build, it is a <tt>const Scalar * const *</tt>. In either case, 
        the syntax is the same, as a two-dimensional data structure. */
    typedef typename Teuchos::ArrayView<const typename Teuchos::ArrayView<const Scalar>::iterator>::iterator const_double_pointer;

    //@}

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Basic MultiVector constuctor.
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos_Ordinal NumVectors, bool zeroOut=true);

    //! MultiVector copy constructor.
    MultiVector(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &source);

    //! \brief Set multi-vector values from array of pointers using C pointers.
    /*! \c CopyView indicates whether the data will be copied from the input array or if the MultiVector object will encapsulate the data throughout its existence.
        \c OwnsMem indicates whether the MultiVector object owns the memory and is therefore responsible for deleting it. This flag is ignored if \c CopyView is \c Teuchos::Copy.
      
       Post-condition: constantStride() == false
     */
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos::DataAccess CopyView, Scalar **ArrayOfPtrs, Teuchos_Ordinal NumVectors, bool OwnsMem = false);

    //! \brief Set multi-vector values from two-dimensional array using a C pointer.
    /*! \c CopyView indicates whether the data will be copied from the input array or if the MultiVector object will encapsulate the data throughout its existence.
        \c OwnsMem indicates whether the MultiVector object owns the memory and is therefore responsible for deleting it. This flag is ignored if \c CopyView is \c Teuchos::Copy.
      
       Post-condition: constantStride() == true
     */
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos::DataAccess CopyView, Scalar *A, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors, bool OwnsMem = false);

    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Scalar> &A, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors);

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, Teuchos_Ordinal NumVectors);

    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (view)
    /*! Post-condition: constantStride() == true */
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayRCP<Scalar> &A, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors);

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (view)
    /*! Post-condition: constantStride() == false */
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayRCP<Scalar> > &ArrayOfPtrs, Teuchos_Ordinal NumVectors);

    //! MultiVector destructor.
    virtual ~MultiVector();

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    void replaceGlobalValue(GlobalOrdinal globalRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    void sumIntoGlobalValue(GlobalOrdinal globalRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    void replaceMyValue(LocalOrdinal myRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    void sumIntoMyValue(LocalOrdinal myRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Initialize all values in a multi-vector with specified value.
    void putScalar(const Scalar &value);

    //! Set multi-vector values to random numbers.
    void random();

    //! Replace the underlying Map with a compatible one.
    void replaceMap(const Map<LocalOrdinal,GlobalOrdinal> &map);

    //! Instruct a local (non-distributed) MultiVector to sum values across all nodes.
    void reduce();

    //! = Operator.
    /*! \param In A - Multivector to copy
     */
    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>& operator=(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &source);

    //@}

    //! @name Data Copy and View extraction methods
    /** These methods are used to extract the data underlying the MultiVector. They return data in one of three forms: 
      - a MultiVector with a subset of the columns of the target MultiVector
      - a raw C pointer or array of raw C pointers
      - one of the Teuchos memory management classes
      Not all of these methods are valid for a particular MultiVector. For instance, calling a method that accesses a 
      view of the data in a 1-D format (i.e., extractView1D) requires that the target MultiVector has constant stride.
     */
    //@{

    //! Returns a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > subCopy(const Teuchos::Range1D &colRng) const;

    //! Returns a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > subCopy(const Teuchos::ArrayView<const Teuchos_Index> &cols) const;

    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > subView(const Teuchos::Range1D &colRng);

    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > subView(const Teuchos::ArrayView<const Teuchos_Index> &cols);

    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > subViewConst(const Teuchos::Range1D &colRng) const;

    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > subViewConst(const Teuchos::ArrayView<const Teuchos_Index> &cols) const;

    //! Vector access function.
    Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal> > operator()(Teuchos_Ordinal j);

    //! Vector access function.
    Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal> > operator() (Teuchos_Ordinal j) const;

    //! Local vector access function.
    //! Pointer to the local values in a particular vector of this multi-vector.
    inline pointer operator[](Teuchos_Ordinal j);

    //! Local vector access function.
    //! Pointer to the local values in a particular vector of this multi-vector.
    inline const_pointer operator[](Teuchos_Ordinal j) const;

    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    void extractCopy1D(typename Teuchos::ArrayView<Scalar> A, Teuchos_Ordinal LDA) const;

    //! Return multi-vector values in user-provided two-dimensional array (using C pointers).
    void extractCopy1D(Scalar *A, Teuchos_Ordinal LDA) const;

    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    void extractCopy2D(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const;

    //! Return multi-vector values in user-provided array of pointers (using C pointers).
    void extractCopy2D(Scalar * const * ArrayOfPtrs) const;

    //! Return non-const non-persisting view of values in a one-dimensional array (using C pointers). Throws std::runtime_error if the underlying data is non-contiguous.
    void extractView1D(Scalar * &A, Teuchos_Ordinal &MyLDA);

    //! Return non-const non-persisting view of values in a one-dimensional array (using Teuchos memory management classes). Throws std::runtime_error if the underlying data is non-contiguous.
    void extractView1D(Teuchos::ArrayView<Scalar> &A, Teuchos_Ordinal &MyLDA);

    //! Return non-const non-persisting pointers to values.
    inline double_pointer extractView2D();

    //! Return const non-persisting view of values in a one-dimensional array (using Teuchos memory management classes). Throws std::runtime_error if the underlying data is non-contiguous.
    void extractConstView1D(Teuchos::ArrayView<const Scalar> &A, Teuchos_Ordinal &MyLDA) const;

    //! Return const non-persisting view of values in a one-dimensional array (using C pointers). Throws std::runtime_error if the underlying data is non-contiguous.
    void extractConstView1D(const Scalar * &A, Teuchos_Ordinal &MyLDA) const;

    //! Return const non-persisting pointers to values.
    inline const_double_pointer extractConstView2D() const;

    //@}

    //! @name Mathematical methods
    //@{ 

    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A, const Teuchos::ArrayView<Scalar> &dots) const;

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A);

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const Scalar &alpha);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(Teuchos::ArrayView<const Scalar> alpha);

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    void scale(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A);

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A, const Scalar &beta);

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A, const Scalar &beta, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &B, const Scalar &gamma);

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute 2-norm of each vector in multi-vector.
    void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    void normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute mean (average) value of each vector in multi-vector.
    void meanValue(const Teuchos::ArrayView<Scalar> &means) const;

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &B, const Scalar &beta);

    //! Multiply a MultiVector with another, element-by-element: this(i,j) = beta*this(i,j) + alpha*A(i,j)*B(i,j)
    void multiply(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &B, const Scalar &beta);

    //! Multiply a MultiVector by the reciprocal of another, element-by-element. this(i,j) = beta*this(i,j) + alpha*B(i,j)/A(i,j)
    void reciprocalMultiply(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &B, const Scalar &beta);

    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    Teuchos_Ordinal numVectors() const;

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    LocalOrdinal myLength() const;

    //! Returns the global vector length of vectors in the multi-vector.
    GlobalOrdinal globalLength() const;

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    Teuchos_Ordinal stride() const;

    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    bool constantStride() const;

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

/*
    //! @name Expert-only unsupported methods
    //@{ 

    //! Reset the view of an existing multivector to point to new user data.
    void resetView(const Teuchos::ArrayRCP<const Teuchos::ArrayRCP<Scalar> > &ArrayOfPtrs);

    //! Get pointer to MultiVector values.
    const Teuchos::ArrayRCP<Scalar> & values();

    //! Get pointer to MultiVector values.
    const Teuchos::ArrayRCP<const Scalar> & valuesConst() const;

    //@}
*/

    protected:

    Teuchos::RCP<MultiVectorData<Scalar> > MVData_;

    // Advanced MultiVector constuctor for creating views.
    MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::RCP<MultiVectorData<Scalar> > &mvdata);

    // four functions needed for DistObject derivation
    bool checkSizes(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj, Teuchos_Ordinal &packetSize);

    void copyAndPermute(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj,
                              Teuchos_Ordinal numSameIDs,
                        const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                        const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs);

    void packAndPrepare(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj,
                        const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                        const Teuchos::ArrayView<Scalar> &exports,
                        Distributor &distor);

    void unpackAndCombine(const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                          const Teuchos::ArrayView<const Scalar> &imports,
                          Distributor &distor,
                          CombineMode CM);

  }; // class MultiVector

} // namespace Tpetra


#endif // TPETRA_MULTIVECTOR_DECL_HPP
