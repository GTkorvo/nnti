/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/

#ifndef TSFEPETRAVECTOR_HPP
#define TSFEPETRAVECTOR_HPP

#include "SundanceDefs.hpp"
#include "SundancePrintable.hpp"
#include "TSFIndexableVector.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFRawDataAccessibleVector.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"
#include "TSFEpetraVectorSpace.hpp"


namespace TSFExtended
{
using Teuchos::Range1D;
using namespace Thyra;
using namespace Teuchos;
/**
 * TSF extension of Thyra::EpetraVector, implementing the LoadableVector
 * interface allowing an application to access elements. This class derives
 * from Thyra::VectorDefaultBase, so it can be used seamlessly in any 
 * Thyra-based code.
 */
class EpetraVector : public Thyra::VectorDefaultBase<double>,
                     public IndexableVector<double>,
                     public RawDataAccessibleVector<double>,
                     public SundanceUtils::Printable
{
public:

  /** Construct with a smart pointer to an Epetra vector space. */
  EpetraVector(const RefCountPtr<const VectorSpaceBase<double> >& vs);

  /** Construct with smart pointers to an Epetra vector space
      and an existing Epetra vector. */
  EpetraVector(const RefCountPtr<const VectorSpaceBase<double> >& vs,
    const RefCountPtr<Epetra_Vector>& vec);


  /** \name VectorBase interface */
  //@{
  /** */
   RefCountPtr< const VectorSpaceBase<double> > 
   space() const {return vecSpace_;}

#ifndef TRILINOS_8
  /** */
  void applyOpImpl(const RTOpPack::RTOpT< double >& op,
		const ArrayView< const Ptr< const VectorBase< double > > > &  	vecs,
		const ArrayView< const Ptr< VectorBase< double > > > &  	targ_vecs,
		const Ptr< RTOpPack::ReductTarget > &  	reduct_obj,
		const Index  	first_ele_offset,
		const Index  	sub_dim,
		const Index  	global_offset	 
    ) const ;
#else
  virtual void applyOp(
    const RTOpPack::RTOpT<double> &op,
    const int num_vecs,
    const VectorBase<double>*const vecs[],
    const int num_targ_vecs,
    VectorBase<double>*const targ_vecs[],
    RTOpPack::ReductTarget *reduct_obj,
    const Index first_ele_offset,
    const Index sub_dim,
    const Index global_offset
    ) const ;
#endif

  /** */
  void acquireDetachedVectorViewImpl(const Range1D& rng,
		RTOpPack::ConstSubVectorView<double>* sub_vec) const ;

  /** */
  void releaseDetachedVectorViewImpl(
    RTOpPack::ConstSubVectorView<double>* sub_vec) const ;

  /** */
  void acquireNonconstDetachedVectorViewImpl(const Range1D& rng,
		RTOpPack::SubVectorView<double> * sub_vec);	 


  /** */
  void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<double>* sub_vec);
  
  //@}

  /** \name IndexableVector interface */
  //@{
  /** read the element at the given global index */
  virtual const double& operator[](Index globalIndex) const 
    {return getElement(globalIndex);}

  /** writable access to the element at the given global index */
  virtual double& operator[](Index globalIndex) ;
  //@}

  /** \name Raw data access interface */
  //@{
  /** */
  virtual const double* dataPtr() const {return &(epetraVec_->operator[](0));}
  /** */
  virtual double* dataPtr() {return &(epetraVec_->operator[](0));}
  //@}

  /** \name LoadableVector interface */
  //@{
  /** set a single element */
  void setElement(Index globalIndex, const double& value);

  /** add to a single element */
  void addToElement(Index globalIndex, const double& value);

  /** set a group of elements */
  void setElements(size_t numElems, const int* globalIndices, 
    const double* values);


  /** add to a group of elements */
  void addToElements(size_t numElems, const int* globalIndices, 
    const double* values);

  /** */
  void finalizeAssembly();
  //@}

  /** \name AccessibleVector interface */
  //@{
  /** */
  const double& getElement(Index globalIndex) const ;

  /** */
  void getElements(const Index* globalIndices, int numElems,
    Teuchos::Array<double>& elems) const ;
  //@}

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  void print(std::ostream& os) const ;
  //@}
      

  /** */
  const RefCountPtr<Epetra_Vector>& epetraVec() const 
    {return epetraVec_;}

  /** */
  RefCountPtr<Epetra_Vector>& epetraVec() {return epetraVec_;}

  /** Get a read-only Epetra_Vector */
  static const Epetra_Vector& getConcrete(const TSFExtended::Vector<double>& tsfVec);
  /** Get a read-write Epetra_Vector */
  static Epetra_Vector& getConcrete(TSFExtended::Vector<double>& tsfVec);
  /** Get a read-write Epetra_Vector pointer */
  static Epetra_Vector* getConcretePtr(TSFExtended::Vector<double>& tsfVec);

  

    

protected:    
  /** */
  const RefCountPtr<const Epetra_Map>& epetraMap() const {return epetraMap_;}

  /** */
  Range1D validateRange(const Range1D& rng) const ;

private:

  RefCountPtr<Epetra_Vector> epetraVec_;

  RefCountPtr<const Thyra::VectorSpaceBase<double> > vecSpace_;

  RefCountPtr<const EpetraVectorSpace> epetraVecSpace_;

  RefCountPtr<const Epetra_Map> epetraMap_;

  int localOffset_;

  int localSubDim_;

  int globalDim_;

  mutable bool in_applyOpImpl_;
};
  
}


#endif
