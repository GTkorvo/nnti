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

#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraVector.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Thyra_DefaultSpmdVectorSpaceFactory.hpp"
#include "TSFOut.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif



using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;

EpetraVectorSpace::EpetraVectorSpace(const RefCountPtr<const Epetra_Map>& m)
  : ScalarProdVectorSpaceBase<double>(),
    SpmdVectorSpaceBase<double>(),
    smallVecSpcFactory_(rcp(new DefaultSpmdVectorSpaceFactory<double>())),
    epetraMap_(m),
    comm_(epetraCommToTeuchosComm(m->Comm())),
    localSubDim_(epetraMap_->NumMyElements()),
    localOffset_(epetraMap_->MinMyGID())
{}


Index EpetraVectorSpace::dim() const 
{
  return epetraMap_->NumGlobalElements();
}

bool EpetraVectorSpace::isCompatible(const VectorSpaceBase<double>& other) const 
{
  const EpetraVectorSpace* epvs = dynamic_cast<const EpetraVectorSpace*>(&other);
  if (epvs != 0)
  {
    return epetraMap_->SameAs(*(epvs->epetraMap_));
  }
  return false;
}

RefCountPtr<const VectorSpaceFactoryBase<double> > 
EpetraVectorSpace::smallVecSpcFcty() const 
{
  return smallVecSpcFactory_;
}



// Overridden from VectorSpace

Teuchos::RefCountPtr<VectorBase<double> >
EpetraVectorSpace::createMember() const
{
//  cout << "creating vector" << endl;
  return rcp(new EpetraVector(rcp(this, false)));
}



Teuchos::RefCountPtr<MultiVectorBase<double> >
EpetraVectorSpace::createMembers(int n) const
{
  RefCountPtr<const VectorSpaceBase<double> > self = rcp(this, false);
  RefCountPtr<const VectorSpaceBase<double> > small 
    = rcp(new DefaultSpmdVectorSpace<double>(n));
  Array<RefCountPtr<VectorBase<double> > > vecs(n);
  for (unsigned int i=0; i<vecs.size(); i++)
    {
      vecs[i] = createMember();
    }
  return rcp(
    new Thyra::DefaultColumnwiseMultiVector<double>(
      self, small,
#ifdef TRILINOS_DEV
      vecs
#else
      &(vecs[0])
#endif
      )
    );
}



Teuchos::RefCountPtr< const VectorSpaceBase<double> >
EpetraVectorSpace::clone() const
{
  return Teuchos::rcp(new EpetraVectorSpace(epetraMap_));
}



string EpetraVectorSpace::description() const
{
  return "EpetraVectorSpace[dim=" + Teuchos::toString(dim()) + ", local="
    + Teuchos::toString(localSubDim()) + "]";
}



Teuchos::RefCountPtr<const Teuchos::Comm<Index> > 
EpetraVectorSpace::epetraCommToTeuchosComm(const Epetra_Comm& epComm) const 
{
  RefCountPtr<const Comm<Index> > rtn;

#ifdef HAVE_MPI
  const Epetra_MpiComm* mpiComm 
    = dynamic_cast<const Epetra_MpiComm*>(&epComm);
#endif

  const Epetra_SerialComm* serialComm 
    = dynamic_cast<const Epetra_SerialComm*>(&epComm);

  if (serialComm != 0)
  {
    rtn  = rcp(new SerialComm<Index>());
  }
#ifdef HAVE_MPI
  else if (mpiComm != 0)
  {
    MPI_Comm rawMpiComm = mpiComm->GetMpiComm();
    RefCountPtr<const OpaqueWrapper<MPI_Comm> > ptr 
      = rcp(new OpaqueWrapper<MPI_Comm>(rawMpiComm));
    rtn  = rcp(new MpiComm<Index>(ptr));
  }
#endif
  else
  {
    TEST_FOR_EXCEPTION(true, std::runtime_error, "Epetra_Comm is neither "
      "a SerialComm or MpiComm");
  }
  return rtn;
}




