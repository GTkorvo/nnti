// File BelosTSFCoreInterface.hpp: interface for the BelosTSFCore class.
//
#ifndef BELOS_TSFCORE_HPP
#define BELOS_TSFCORE_HPP

/*! \file BelosTSFCoreInterface.hpp
  \brief Provides several interfaces between Belos virtual classes and the TSFCore virtual classes.
*/

// Belos files
#include "BelosOperator.hpp"
#include "BelosMultiVec.hpp"
#include "BelosReturnType.hpp"
#include "BelosConfigDefs.hpp"

// TSFCore files
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreLinearOp.hpp"

// Teuchos files
#include "Teuchos_RefCountPtr.hpp"

namespace Belos {
//
//--------template class Belos::TSFCoreMat-----------------------
//
template <class TYPE> 
class TSFCoreMat : public Operator<TYPE> {
public:
  TSFCoreMat( Teuchos::RefCountPtr< TSFCore::LinearOp<TYPE> > Op_in ) { Op = Op_in; }
  ReturnType Apply ( const MultiVec<TYPE>& x, 
		     MultiVec<TYPE>& y ) const;
private:
  Teuchos::RefCountPtr< TSFCore::LinearOp<TYPE> > Op;
};

//------------------------------------------------------------
//
//--------template class Belos::TSFCoreVec-------------------------------------
//

template <class TYPE>
class TSFCoreVec : public MultiVec<TYPE> {
public:
  //! Enumeration for accessing the data in the TSFCore multivectors.
  enum DataAccess { Copy, /*!< Deep Copy */
		    View  /*!< Shallow Copy */
  };

  friend class TSFCoreMat<TYPE>;
  // constructors
  TSFCoreVec( Teuchos::RefCountPtr<const TSFCore::VectorSpace<TYPE> > source_space, const int NumVecs );
  TSFCoreVec( TSFCore::MultiVector<TYPE>& source ) { TSFCoreMV = Teuchos::rcp( &source, false); }
  TSFCoreVec( const TSFCoreVec<TYPE>& source );
  TSFCoreVec( DataAccess type, TSFCoreVec<TYPE>& source, int index[], int NumVecs); 
  //
  //  member functions inherited from MultiVec
  //
  //  the following is a virtual copy constructor returning
  //  a pointer to the pure virtual class. vector values are
  //  not copied; instead a new MultiVec is created containing
  //  a non-zero amount of columns.
  //
  MultiVec<TYPE> * Clone ( const int );
  //
  //  the following is a virtual copy constructor returning
  //  a pointer to the pure virtual class. vector values are
  //  copied and a new stand-alone MultiVector is created.
  //  (deep copy).
  //
  MultiVec<TYPE> * CloneCopy ();
  //
  //  Selective deep copy (or copy) constructor.
  //
  MultiVec<TYPE> * CloneCopy ( int [], int );
  //
  //  the following is a virtual view constructor returning
  //  a pointer to the pure virtual class. vector values are 
  //  shared and hence no memory is allocated for the columns.
  //
  MultiVec<TYPE> * CloneView ( int [], int );
  //
  int GetNumberVecs () const;
  int GetVecLength () const;
  //
  //  set a block of this multivec with the multivecs specified by
  //  the index.
  //
  void SetBlock ( MultiVec<TYPE>& A, int index[], 
			  int NumVecs ); 
  //
  // *this <- alpha * A * B + beta * (*this)
  //
  void MvTimesMatAddMv ( TYPE alpha, MultiVec<TYPE>& A, 
				 Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta );
  //
  // *this <- alpha * A + beta * B
  //
  void MvAddMv ( TYPE alpha , MultiVec<TYPE>& A, TYPE beta,
			 MultiVec<TYPE>& B);
  //
  // B <- alpha * A^T * (*this)
  //
  void MvTransMv ( TYPE alpha, MultiVec<TYPE>& A, 
			   Teuchos::SerialDenseMatrix<int,TYPE>& B );
  //
  // alpha[i] = norm of i-th column of (*this)
  //
  ReturnType MvNorm ( TYPE* normvec, NormType norm_type = TwoNorm );
  //
  // random vectors in i-th column of (*this)
  //
  void MvRandom();
  //
  // initializes each element of (*this) with alpha
  //
  void MvInit ( TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() );
  //
private:
  // Data container
  Teuchos::RefCountPtr< TSFCore::MultiVector<TYPE> > TSFCoreMV;
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

template<class TYPE>
TSFCoreVec<TYPE>::TSFCoreVec( Teuchos::RefCountPtr<const TSFCore::VectorSpace<TYPE> > source_space,
			      const int NumVecs )
{
  TSFCoreMV = source_space->createMembers( NumVecs );
}
  
template<class TYPE>
TSFCoreVec<TYPE>::TSFCoreVec( const TSFCoreVec<TYPE>& source )
{
  TSFCoreMV = (source.TSFCoreMV)->clone_mv();
}
  
template<class TYPE>
TSFCoreVec<TYPE>::TSFCoreVec( DataAccess type, TSFCoreVec<TYPE>& source, 
			       	int index[], int NumVecs )
{
  // Alter indexing for one-based indexing in TSFCore.
  for (int i=0; i<NumVecs; i++) index[i]++;
  //
  if (type == Copy) {
    Teuchos::RefCountPtr< TSFCore::MultiVector<TYPE> > tempMV = (source.TSFCoreMV)->subView( NumVecs, index );
    TSFCoreMV = tempMV->clone_mv();
  } 
  else {
    TSFCoreMV = (source.TSFCoreMV)->subView( NumVecs, index );
  }
}

//
//  member functions inherited from MultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to TSFCoreVec<TYPE>
//  (the derived type) instead of a pointer to the pure virtual base class.
//
template<class TYPE>
MultiVec<TYPE>* TSFCoreVec<TYPE>::Clone ( const int NumVecs ) {
  TSFCoreVec *ptr_alv = new TSFCoreVec<TYPE>( TSFCoreMV->range(), NumVecs);
  return ptr_alv; // safe upcast.
}
//
//  the following is a virtual copy constructor returning
//  a pointer to the pure virtual class. vector values are
//  copied.
//
template<class TYPE>
MultiVec<TYPE>* TSFCoreVec<TYPE>::CloneCopy() {
  TSFCoreVec *ptr_alv = new TSFCoreVec<TYPE>(*this);
  return ptr_alv; // safe upcast
}

template<class TYPE>
MultiVec<TYPE>* TSFCoreVec<TYPE>::CloneCopy ( int index[], int NumVecs ) {
  
  TSFCoreVec *ptr_alv = new TSFCoreVec<TYPE>( Copy, *this, index, NumVecs );
  return ptr_alv; // safe upcast.
}

template<class TYPE>
MultiVec<TYPE>* TSFCoreVec<TYPE>::CloneView ( int index[], int NumVecs ) {
  
  TSFCoreVec *ptr_alv = new TSFCoreVec<TYPE>( View, *this, index, NumVecs );
  return ptr_alv; // safe upcast.
}

template<class TYPE>
int TSFCoreVec<TYPE>::GetNumberVecs () const {
  return (TSFCoreMV->domain())->dim();
}

template<class TYPE>
int TSFCoreVec<TYPE>::GetVecLength () const {
  return (TSFCoreMV->range())->dim();
}

template<class TYPE>
void TSFCoreVec<TYPE>::SetBlock( MultiVec<TYPE>& A, int index[], int NumVecs ) {
  
  TSFCoreVec<TYPE> *A_vec = dynamic_cast<TSFCoreVec<TYPE> *>(&A); assert(A_vec!=NULL);
  //
  // Alter indexing for one-based indexing in TSFCore.
  for (int i=0; i<NumVecs; i++) index[i]++;
  //
  Teuchos::RefCountPtr<TSFCore::MultiVector<TYPE> > tempMV = TSFCoreMV->subView( NumVecs, index );
  TSFCore::assign( tempMV.get(), *(A_vec->TSFCoreMV->subView( TSFCore::Range1D( 1, NumVecs ) )) );
}
//
// *this <- alpha * A * B + beta * (*this)
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvTimesMatAddMv ( TYPE alpha, MultiVec<TYPE>& A, 
					 Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta ) {
  
  TSFCoreVec<TYPE> *A_vec = dynamic_cast<TSFCoreVec<TYPE>*>(&A); assert(A_vec!=NULL);
  TSFCore::scale( beta, TSFCoreMV.get() );
  for (int j=0; j<B.numCols(); j++) {
    for (int i=0; i<B.numRows(); i++) {
      TSFCore::Vp_StV( (TSFCoreMV->col(j+1)).get(), alpha*B(i,j), *((A_vec->TSFCoreMV)->col(i+1)) );
    }
  }
}
//
// *this <- alpha * A + beta * B
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvAddMv ( TYPE alpha , MultiVec<TYPE>& A, 
				 TYPE beta, MultiVec<TYPE>& B) {
  TSFCoreVec<TYPE> *A_vec = dynamic_cast<TSFCoreVec<TYPE>*>(&A); assert(A_vec!=NULL);
  TSFCoreVec<TYPE> *B_vec = dynamic_cast<TSFCoreVec<TYPE>*>(&B); assert(B_vec!=NULL);
  //
  // We must be aware of any of the input MultiVec's that are the same as *this.
  // 
  if ( (A_vec->TSFCoreMV).get() == TSFCoreMV.get() ) {
    //
    // *this *= alpha
    TSFCore::scale( alpha, TSFCoreMV.get() );
    //
    // *this += beta * B
    TSFCore::update( beta, *(B_vec->TSFCoreMV), TSFCoreMV.get() );
    //
  } else if ( (B_vec->TSFCoreMV).get() == TSFCoreMV.get() ) { 
    //
    // *this *=beta
    TSFCore::scale( beta, TSFCoreMV.get() );
    //
    // *this += alpha * A
    TSFCore::update( alpha, *(A_vec->TSFCoreMV), TSFCoreMV.get() );
    //
  } else {
    // *this <- A
    TSFCore::assign( TSFCoreMV.get(), *(A_vec->TSFCoreMV) );
    //
    // *this *= alpha
    TSFCore::scale( alpha, TSFCoreMV.get() );
    //
    // *this += beta * B
    TSFCore::update( beta, *(B_vec->TSFCoreMV), TSFCoreMV.get() );
  }
}
//
// dense B <- alpha * A^T * (*this)
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvTransMv ( TYPE alpha, MultiVec<TYPE>& A,
				   Teuchos::SerialDenseMatrix<int,TYPE>& B) {
  
  TSFCoreVec<TYPE> *A_vec = dynamic_cast<TSFCoreVec<TYPE>*>(&A); assert(A_vec!=NULL);
  //
  // The indexing for the "col" method is 1-based, so the index is altered in the "dot" call.
  for (int j=0; j<B.numCols(); j++) {
    for (int i=0; i<B.numRows(); i++) {
      B(i,j) =  TSFCore::dot( *(TSFCoreMV->col(j+1)), *((A_vec->TSFCoreMV)->col(i+1)) );
    }
  }
  B.scale( alpha );
}
//
// array[i] = norm of i-th column of (*this)
//
template<class TYPE>
ReturnType TSFCoreVec<TYPE>::MvNorm ( TYPE * normvec, NormType norm_type ) 
{
  if (normvec) {
    switch( norm_type ) {
    case ( OneNorm ) :
      for (int i=0; i<GetNumberVecs(); i++)
	normvec[i] = TSFCore::norm_1( *(TSFCoreMV->col(i+1)) );
      return Ok;
    case ( TwoNorm ) :
      for (int i=0; i<GetNumberVecs(); i++)
	normvec[i] = TSFCore::norm_2( *(TSFCoreMV->col(i+1)) );
      return Ok;
    case ( InfNorm ) :
      for (int i=0; i<GetNumberVecs(); i++) 
	normvec[i] = TSFCore::norm_inf( *(TSFCoreMV->col(i+1)) );
      return Ok;
    default :
      return Undefined;
    }
  }
  return Undefined;
}
//
// random vectors in i-th column of (*this)
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvRandom () 
{
  TSFCore::randomize( 0.0, 1.0, TSFCoreMV.get() );
}
//
// initializes each element of (*this) with alpha
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvInit ( TYPE alpha ) 
{
  TSFCore::assign( TSFCoreMV.get(), alpha );
}

//--------------------------------------------------------------
//
// implementation of the Belos::TSFCoreMat class.
//
template <class TYPE>
ReturnType TSFCoreMat<TYPE>::Apply ( const MultiVec<TYPE>& x, 
				     MultiVec<TYPE>& y ) const 
{	
  MultiVec<TYPE> &temp_x = const_cast<MultiVec<TYPE> &>(x);
  TSFCoreVec<TYPE> *x_vec = dynamic_cast<TSFCoreVec<TYPE> *>(&temp_x); assert(x_vec!=NULL);
  TSFCoreVec<TYPE> *y_vec = dynamic_cast<TSFCoreVec<TYPE> *>(&y); assert(y_vec!=NULL);
  Op->apply( TSFCore::NOTRANS, *(x_vec->TSFCoreMV), (y_vec->TSFCoreMV).get() );
  return Ok;
}


} // end Belos namespace
#endif 
 // end of file BELOS_TSFCORE_HPP
