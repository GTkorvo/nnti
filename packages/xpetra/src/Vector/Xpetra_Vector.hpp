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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_VECTOR_HPP
#define XPETRA_VECTOR_HPP

/* this file is automatically generated - do not edit (see script/interfaces.py) */

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MultiVector.hpp"

namespace Xpetra {

  template <class Scalar = MultiVector<>::scalar_type,
            class LocalOrdinal =
              typename MultiVector<Scalar>::local_ordinal_type,
            class GlobalOrdinal =
              typename MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node =
              typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class Vector
    : public virtual MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >
  {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::dot;          // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::norm1;        // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::norm2;        // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::normInf;      // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::normWeighted; // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::meanValue;    // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::replaceGlobalValue; // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::sumIntoGlobalValue; // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::replaceLocalValue; // overloading, not hiding
    using MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::sumIntoLocalValue; // overloading, not hiding

    //! @name Constructor/Destructor Methods
    //@{

    //! Destructor.
    virtual ~Vector() { }

   //@}

    //! @name Post-construction modification routines
    //@{

    //! Replace current value at the specified location with specified value.
    virtual void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value)= 0;

    //! Adds specified value to existing value at the specified location.
    virtual void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value)= 0;

    //! Replace current value at the specified location with specified values.
    virtual void replaceLocalValue(LocalOrdinal myRow, const Scalar &value)= 0;

    //! Adds specified value to existing value at the specified location.
    virtual void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value)= 0;

    //@}

    //! @name Mathematical methods
    //@{

    //! Computes dot product of this Vector against input Vector x.
    virtual Scalar dot(const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &a) const = 0;

    //! Return 1-norm of this Vector.
    virtual typename Teuchos::ScalarTraits< Scalar >::magnitudeType norm1() const = 0;

    //! Compute 2-norm of this Vector.
    virtual typename Teuchos::ScalarTraits< Scalar >::magnitudeType norm2() const = 0;

    //! Compute Inf-norm of this Vector.
    virtual typename Teuchos::ScalarTraits< Scalar >::magnitudeType normInf() const = 0;

    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    virtual typename Teuchos::ScalarTraits< Scalar >::magnitudeType normWeighted(const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &weights) const = 0;

    //! Compute mean (average) value of this Vector.
    virtual Scalar meanValue() const = 0;

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    virtual std::string description() const = 0;

    //! Print the object with some verbosity level to an FancyOStream object.
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

    //@}

  }; // Vector class

} // Xpetra namespace

#define XPETRA_VECTOR_SHORT
#endif // XPETRA_VECTOR_HPP
