// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
//
#ifndef ANASAZI_MULTI_VEC_TRAITS_HPP
#define ANASAZI_MULTI_VEC_TRAITS_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Anasazi {

  template< class TYPE, class MV >
  struct UndefinedMultiVecTraits
  {
    //! This function should not compile if there is an attempt to instantiate!
    static inline TYPE notDefined() { return MV::this_type_is_missing_a_specialization(); };
  };
  
  template<class TYPE, class MV>
  class MultiVecTraits 
  {
  public:
    ///
    static Teuchos::RefCountPtr<MV> Clone( const MV& mv, const int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<MV> CloneCopy( const MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<MV> CloneCopy( const MV& mv, int index[], int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<MV> CloneView( MV& mv, int index[], int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<const MV> CloneView( const MV& mv, int index[], int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static int GetVecLength( const MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return 0; }     
    ///
    static int GetNumberVecs( const MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return 0; }     
    ///
    static void MvTimesMatAddMv( TYPE alpha, const MV& A, 
				 const Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta, MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvAddMv( TYPE alpha, const MV& A, TYPE beta, const MV& B, MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvTransMv( const MV& mv, TYPE alpha, const MV& A, Teuchos::SerialDenseMatrix<int,TYPE>& B )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvNorm( const MV& mv, TYPE *normvec )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void SetBlock( const MV& A, int index[], int numvecs, MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvRandom( MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvInit( MV& mv, TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvPrint( const MV& mv, ostream& os )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
  };
  
} // namespace Anasazi

#endif // ANASAZI_MULTI_VEC_TRAITS_HPP
