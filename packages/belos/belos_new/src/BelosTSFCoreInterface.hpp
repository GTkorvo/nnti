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
//#include "TSFCoreVectorSpaceDecl.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
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
  TSFCoreVec( Teuchos::RefCountPtr<const TSFCore::VectorSpace<TYPE> >& source_space, const int NumVecs );
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
  virtual MultiVec<TYPE> * Clone ( const int );
  //
  //  the following is a virtual copy constructor returning
  //  a pointer to the pure virtual class. vector values are
  //  copied and a new stand-alone MultiVector is created.
  //  (deep copy).
  //
  virtual MultiVec<TYPE> * CloneCopy ();
  //
  //  Selective deep copy (or copy) constructor.
  //
  virtual MultiVec<TYPE> * CloneCopy ( int [], int );
  //
  //  the following is a virtual view constructor returning
  //  a pointer to the pure virtual class. vector values are 
  //  shared and hence no memory is allocated for the columns.
  //
  virtual MultiVec<TYPE> * CloneView ( int [], int );
  //
  virtual int GetNumberVecs () const;
  virtual int GetVecLength () const;
  //
  //  set a block of this multivec with the multivecs specified by
  //  the index.
  //
  virtual void SetBlock ( MultiVec<TYPE>& A, int index[], 
			  int NumVecs ); 
  //
  // *this <- alpha * A * B + beta * (*this)
  //
  virtual void MvTimesMatAddMv ( TYPE alpha, MultiVec<TYPE>& A, 
				 Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta );
  //
  // *this <- alpha * A + beta * B
  //
  virtual void MvAddMv ( TYPE alpha , MultiVec<TYPE>& A, TYPE beta,
			 MultiVec<TYPE>& B);
  //
  // B <- alpha * A^T * (*this)
  //
  virtual void MvTransMv ( TYPE alpha, MultiVec<TYPE>& A, 
			   Teuchos::SerialDenseMatrix<int,TYPE>& B );
  //
  // alpha[i] = norm of i-th column of (*this)
  //
  virtual ReturnType MvNorm ( TYPE* normvec, NormType norm_type = TwoNorm );
  //
  // random vectors in i-th column of (*this)
  //
  virtual void MvRandom();
  //
  // initializes each element of (*this) with alpha
  //
  virtual void MvInit ( TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() );
  //
private:
  // Data container
  Teuchos::RefCountPtr< TSFCore::MultiVector<TYPE> > TSFCoreMV;
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

template<class TYPE>
TSFCoreVec<TYPE>::TSFCoreVec( Teuchos::RefCountPtr<const TSFCore::VectorSpace<TYPE> >& source_space,
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
	TSFCoreVec *ptr_alv = new TSFCoreVec( TSFCoreMV->domain(), NumVecs);
	return ptr_alv; // safe upcast.
}
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied.
	//
template<class TYPE>
MultiVec<TYPE>* TSFCoreVec<TYPE>::CloneCopy() {
	TSFCoreVec *ptr_alv = new TSFCoreVec(*this);
	return ptr_alv; // safe upcast
}

template<class TYPE>
MultiVec<TYPE>* TSFCoreVec<TYPE>::CloneCopy ( int index[], int NumVecs ) {
	
	TSFCoreVec *ptr_alv = new TSFCoreVec( Copy, *this, index, NumVecs );
	return ptr_alv; // safe upcast.
}

template<class TYPE>
MultiVec<TYPE>* TSFCoreVec<TYPE>::CloneView ( int index[], int NumVecs ) {
	
	TSFCoreVec *ptr_alv = new TSFCoreVec( View, *this, index, NumVecs );
	return ptr_alv; // safe upcast.
}

template<class TYPE>
int TSFCoreVec<TYPE>::GetNumberVecs () const {
        return (TSFCoreMV->range())->dim();
}

template<class TYPE>
int TSFCoreVec<TYPE>::GetVecLength () const {
	return (TSFCoreMV->domain())->dim();
}

template<class TYPE>
void TSFCoreVec<TYPE>::SetBlock( MultiVec<TYPE>& A, int index[], int NumVecs ) {

	TSFCoreVec *A_vec = dynamic_cast<TSFCoreVec *>(&A); assert(A_vec!=NULL);
	Teuchos::RefCountPtr<TSFCore::MultiVector<TYPE> > tempMV = TSFCoreMV->subView( NumVecs, index );
	TSFCore::assign( tempMV.get(), *(A_vec->TSFCoreMV->subView( TSFCore::Range1D( 1, NumVecs ) )) );
	delete &tempMV;
}
//
// *this <- alpha * A * B + beta * (*this)
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvTimesMatAddMv ( TYPE alpha, MultiVec<TYPE>& A, 
						   Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta ) 
{
	TSFCoreVec *A_vec = dynamic_cast<TSFCoreVec *>(&A); assert(A_vec!=NULL);
}
//
// *this <- alpha * A + beta * B
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvAddMv ( TYPE alpha , MultiVec<TYPE>& A, 
						   TYPE beta, MultiVec<TYPE>& B) {
	TSFCoreVec *A_vec = dynamic_cast<TSFCoreVec *>(&A); assert(A_vec!=NULL);
	TSFCoreVec *B_vec = dynamic_cast<TSFCoreVec *>(&B); assert(B_vec!=NULL);
}
//
// dense B <- alpha * A^T * (*this)
//
template<class TYPE>
void TSFCoreVec<TYPE>::MvTransMv ( TYPE alpha, MultiVec<TYPE>& A,
						   Teuchos::SerialDenseMatrix<int,TYPE>& B) {

	TSFCoreVec *A_vec = dynamic_cast<TSFCoreVec *>(&A); assert(A_vec!=NULL);
	TSFCore::dot( *(A_vec->TSFCoreMV), *TSFCoreMV, B.values() );
	B.scale( alpha );
}
//
// array[i] = norm of i-th column of (*this)
//
template<class TYPE>
ReturnType TSFCoreVec<TYPE>::MvNorm ( TYPE * normvec, NormType norm_type ) 
{
	if (normvec) {
		for (int i=0; i<(TSFCoreMV->range())->dim(); i++) {
			normvec[i] = TSFCore::norm_2( *(TSFCoreMV->col(i)) );
		}
	}
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
	(*Op)->apply( TSFCore::NOTRANS, x_vec->TSFCoreMV, *y_vec );
        return Ok;
}


///////////////////////////////////////////////////////////////
//--------template class Belos::PetraPrec--------------------

template <class TYPE>
class TSFCorePrec : public Operator<TYPE> {
public:
        TSFCorePrec(const TSFCore::LinearOp<TYPE>& Prec_in);
        ~TSFCorePrec();
        void Apply ( const MultiVec<TYPE>& x, MultiVec<TYPE>& y ) const;
private:
   	const TSFCore::LinearOp<TYPE>& TSFCore_Prec;
};
//--------------------------------------------------------------
//
// implementation of the Belos::TSFCorePrec class.
//
// Constructor.
//
template <class TYPE>
TSFCorePrec<TYPE>::TSFCorePrec(const TSFCore::LinearOp<TYPE>& Prec_in) : TSFCore_Prec(Prec_in) 
{
}
//
// Destructor.
//
template <class TYPE>
TSFCorePrec<TYPE>::~TSFCorePrec() 
{
}
/////////////////////////////////////////////////////////////
//
// Precondition Operator application
//
template <class TYPE>
void TSFCorePrec<TYPE>::Apply ( const MultiVec<TYPE>& x,
				MultiVec<TYPE>& y) const {
	int info=0;
	MultiVec<TYPE> & temp_x = const_cast<MultiVec<TYPE> &>(x);
	TSFCoreVec<TYPE>* vec_x = dynamic_cast<TSFCoreVec<TYPE>* >(&temp_x); assert(vec_x != NULL);
	TSFCoreVec<TYPE>* vec_y = dynamic_cast<TSFCoreVec<TYPE>* >(&y); assert(vec_y != NULL);
	Op.apply( TSFCore::NOTRANS, vec_x->TSFCoreMV, &vec_y->TSFCoreMV );
}


} // end Belos namespace
#endif 
 // end of file BELOS_TSFCORE_HPP
