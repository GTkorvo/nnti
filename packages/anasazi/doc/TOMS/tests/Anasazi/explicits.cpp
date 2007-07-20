#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
template class Teuchos::RefCountPtr<Epetra_CrsMatrix>;
template class Teuchos::RefCountPtr<Epetra_Operator>;
template class Teuchos::RefCountPtr<Epetra_MultiVector>;
template class Teuchos::RefCountPtr<Anasazi::EpetraMultiVec>;
template Teuchos::RefCountPtr<Epetra_Operator>::RefCountPtr(Teuchos::RefCountPtr<Epetra_CrsMatrix> const&);
template class Teuchos::RefCountPtr<Anasazi::EpetraSymOp>;

#include "AnasaziMultiVec.hpp"
template class Anasazi::MultiVec<double>;
template class Teuchos::RefCountPtr<Anasazi::MultiVec<double> >;
template Teuchos::RefCountPtr<Anasazi::MultiVec<double> const>::RefCountPtr(Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> const&);

#include "AnasaziOperator.hpp"
template class Anasazi::Operator<double>;
template class Teuchos::RefCountPtr<Anasazi::Operator<double> >;
template class Teuchos::RefCountPtr<Anasazi::Operator<double> const>;
template Teuchos::RefCountPtr<Anasazi::Operator<double> const>::RefCountPtr(Teuchos::RefCountPtr<Anasazi::EpetraSymOp> const&);

#include "AnasaziBasicOutputManager.hpp"
template class Anasazi::BasicOutputManager<double>;
template class Teuchos::RefCountPtr<Anasazi::OutputManager<double> >;

#include "AnasaziBasicEigenproblem.hpp"
template class Anasazi::BasicEigenproblem<double,Epetra_MultiVector,Epetra_Operator>;
template class Teuchos::RefCountPtr<Anasazi::Eigenproblem<double, Epetra_MultiVector, Epetra_Operator> >;
template class Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> >;
template Teuchos::RefCountPtr<Anasazi::Eigenproblem<double, Epetra_MultiVector, Epetra_Operator> >::RefCountPtr(Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> > const&);

#include "AnasaziBlockDavidsonSolMgr.hpp"
template class Anasazi::BlockDavidsonSolMgr<double,Epetra_MultiVector,Epetra_Operator>;

#include "AnasaziBlockKrylovSchur.hpp"
template class Anasazi::BlockKrylovSchur<double,Epetra_MultiVector,Epetra_Operator>;

#include "AnasaziBlockKrylovSchurSolMgr.hpp"
template class Anasazi::BlockKrylovSchurSolMgr<double,Epetra_MultiVector,Epetra_Operator>;
// template void Anasazi::scaleRitzVectors<std::complex<double> >(std::vector<Teuchos::ScalarTraits<std::complex<double> >::magnitudeType, std::allocator<Teuchos::ScalarTraits<std::complex<double> >::magnitudeType> > const&, Teuchos::SerialDenseMatrix<int, std::complex<double> >*);



#include "AnasaziLOBPCGSolMgr.hpp"
template class Anasazi::LOBPCGSolMgr<double,Epetra_MultiVector,Epetra_Operator>;

// template class BasicEigenproblem<double,Thyra::MultiVectorBase<double>,Thyra::LinearOpBase<double> >;
// template class BlockDavidsonSolMgr<double,Thyra::MultiVectorBase<double>,Thyra::LinearOpBase<double> >;


#include "ModeLaplace1DQ1.h"
template class Teuchos::RefCountPtr<ModalProblem>;
template class Teuchos::RefCountPtr<ModeLaplace1DQ1>;
template Teuchos::RefCountPtr<ModalProblem>::RefCountPtr(Teuchos::RefCountPtr<ModeLaplace1DQ1> const&);

#include "ModeLaplace2DQ1.h"
template class Teuchos::RefCountPtr<ModeLaplace2DQ1>;
template Teuchos::RefCountPtr<ModalProblem>::RefCountPtr(Teuchos::RefCountPtr<ModeLaplace2DQ1> const&);

#include "AnasaziMVOPTester.hpp"
template bool Anasazi::TestMultiVecTraits<double, Anasazi::MultiVec<double> >(Teuchos::RefCountPtr<Anasazi::OutputManager<double> > const&, Teuchos::RefCountPtr<Anasazi::MultiVec<double> const> const&);
template bool Anasazi::TestOperatorTraits<double, Anasazi::MultiVec<double>, Anasazi::Operator<double> >(Teuchos::RefCountPtr<Anasazi::OutputManager<double> > const&, Teuchos::RefCountPtr<Anasazi::MultiVec<double> const> const&, Teuchos::RefCountPtr<Anasazi::Operator<double> const> const&);

