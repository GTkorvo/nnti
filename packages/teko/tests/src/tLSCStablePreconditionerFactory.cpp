#include "tLSCStablePreconditionerFactory.hpp"
#include "NS/PB_LSCPreconditionerFactory.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"

#include <vector>

// This whole test rig is based on inverting the matrix
// 
//      [  1  2  1 -1 ]
//  A = [  2  1 -3  1 ]
//      [  1 -3  0  0 ]
//      [ -1  1  0  0 ]
//
// see the matlab file

namespace PB {
namespace Test {

using namespace Teuchos;
using namespace Thyra;
using namespace PB::NS;

void tLSCStablePreconditionerFactory::initializeTest()
{
   std::vector<int> indicies(2);
   std::vector<double> row0(2),row1(2);

   tolerance_ = 1.0e-13;

   comm = rcp(new Epetra_SerialComm());
   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));

   const RCP<Epetra_CrsMatrix> ptrF  = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   const RCP<Epetra_CrsMatrix> ptrB  = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   const RCP<Epetra_CrsMatrix> ptrBt = rcp(new Epetra_CrsMatrix(Copy,*map,2));

   const RCP<Epetra_CrsMatrix> ptrInvF = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   const RCP<Epetra_CrsMatrix> ptrInvS = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   const RCP<Epetra_CrsMatrix> ptrInvMass = rcp(new Epetra_CrsMatrix(Copy,*map,2));

   indicies[0] = 0;
   indicies[1] = 1;

   // build F matrix
   row0[0] = 1.0; row0[1] = 2.0; 
   row1[0] = 2.0; row1[1] = 1.0; 
   ptrF->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   ptrF->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   ptrF->FillComplete();
   F_ = Thyra::epetraLinearOp(ptrF,"ptrF");
   
   // build B matrix
   row0[0] =  1.0; row0[1] = -3.0; 
   row1[0] = -1.0; row1[1] =  1.0; 
   ptrB->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   ptrB->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   ptrB->FillComplete();
   B_ = Thyra::epetraLinearOp(ptrB,"ptrB");
   
   // build Bt matrix
   row0[0] =  1.0; row0[1] = -1.0; 
   row1[0] = -3.0; row1[1] =  1.0; 
   ptrBt->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   ptrBt->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   ptrBt->FillComplete();
   Bt_ = Thyra::epetraLinearOp(ptrBt,"ptrBt");

   // build inv(F) matrix
   row0[0] = -1.0/3.0; row0[1] =  2.0/3.0;
   row1[0] =  2.0/3.0; row1[1] = -1.0/3.0;
   ptrInvF->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   ptrInvF->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   ptrInvF->FillComplete();
   invF_ = Thyra::epetraLinearOp(ptrInvF,"ptrInvF");

   // build inv(Mass) matrix
   row0[0] =  3.0; row0[1] =  0.0;
   row1[0] =  0.0; row1[1] =  2.0;
   ptrInvMass->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   ptrInvMass->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   ptrInvMass->FillComplete();
   invMass_ = Thyra::epetraLinearOp(ptrInvMass,"ptrInvMass");

   // build inv(Pschur) matrix
   row0[0] = 0.208333333333333; row0[1] = 0.375000000000000;
   row1[0] = 0.375000000000000; row1[1] = 0.875000000000000;
   ptrInvS->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   ptrInvS->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   ptrInvS->FillComplete();
   invBQBt_ = Thyra::scale<double>(-1.0,Thyra::epetraLinearOp(ptrInvS,"ptrInvS"));

   A_ = Thyra::block2x2<double>(F_,Bt_,B_,Thyra::zero(Bt_->range(),B_->domain()),"A");
}

int tLSCStablePreconditionerFactory::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLSCStablePreconditionerFactory";

   status = test_createPrec(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"createPrec\" ... PASSED","   \"createPrec\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_initializePrec(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"initializePrec\" ... PASSED","   \"initializePrec\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_uninitializePrec(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"uninitializePrec\" ... PASSED","   \"uninitializePrec\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_isCompatable(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"isCompatable\" ... PASSED","   \"isCompatable\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_identity(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"identity\" ... PASSED","   \"identity\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_diagonal(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"diagonal\" ... PASSED","   \"diagonal\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_result(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"result\" ... PASSED","   \"result\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tLSCStablePreconditionedFactory...PASSED","tLSCStablePreconditionedFactory...FAILED");
   }
   else {// Normal Operatoring Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tLSCStablePreconditionedFactory...FAILED");
   }

   return failcount;
}

bool tLSCStablePreconditionerFactory::test_createPrec(int verbosity,std::ostream & os)
{
   //RCP<LSCStablePreconditionerFactory> fact = rcp(new LSCStablePreconditionerFactory(invF_,invBQBt_));
   const RCP<const Thyra::PreconditionerFactoryBase<double> > fact = rcp(new LSCPreconditionerFactory(invF_,invBQBt_,Teuchos::null));

   try {
      // preconditioner factory should return a DefaultPreconditionerBase
      rcp_dynamic_cast<DefaultPreconditioner<double> >(fact->createPrec(),true);
   }
   catch(std::exception & e) {
      // if the dynamic cast fails...so does the test
      os << std::endl << "   test_createPrec: dynamic cast to \"DefaultPreconditioner\" FAILED" << std::endl;
      os << "   Descriptive exception \"" << e.what() << "\""<< std::endl;

      return false;
   }

   return true;
}

bool tLSCStablePreconditionerFactory::test_initializePrec(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   // Build block2x2 preconditioner
   //RCP<Thyra::PreconditionerFactoryBase<double> > precFactory 
   //      = rcp(new LSCStablePreconditionerFactory(invF_,invBQBt_));
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory 
         = rcp(new LSCPreconditionerFactory(invF_,invBQBt_,Teuchos::null));
   RCP<Thyra::PreconditionerBase<double> > prec = precFactory->createPrec();

   // initialize the preconditioner
   precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

   RCP<const Thyra::LinearOpBase<double> > op;

   op = prec->getUnspecifiedPrecOp();
   status = (op!=Teuchos::null);
   if(not status) {
      os << std::endl << "   tLSCStablePreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
      os << "      " << "Preconditioner \"getUnspecifiedPrecOp\" is null (it should not be!)" << std::endl;;
   }
   allPassed &= status;

   op = prec->getRightPrecOp();
   status = (op==Teuchos::null);
   if(not status) {
      os << std::endl << "   tLSCStablePreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
      os << "      " << "Preconditioner \"getRightPrecOp\" is not null (it should be!)" << std::endl;;
   }
   allPassed &= status;

   op = prec->getLeftPrecOp();
   status = (op==Teuchos::null);
   if(not status) {
      os << std::endl << "   tLSCStablePreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
      os << "      " << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;;
   }
   allPassed &= status;

   return allPassed;
}

bool tLSCStablePreconditionerFactory::test_uninitializePrec(int verbosity,std::ostream & os)
{
   return true;
}

bool tLSCStablePreconditionerFactory::test_isCompatable(int verbosity,std::ostream & os)
{
   return true;
}

bool tLSCStablePreconditionerFactory::test_identity(int verbosity,std::ostream & os)
{
   // make sure the preconditioner is working by testing against the identity matrix
   typedef RCP<const Thyra::VectorBase<double> > Vector;
   typedef RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;
   typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

   bool status = false;
   bool allPassed = true;

   LinearOp Iu = Thyra::identity<double>(invF_->range());
   LinearOp Ip = Thyra::identity<double>(invBQBt_->range());
   LinearOp Zp = Thyra::zero<double>(invBQBt_->range(),invF_->domain());
   LinearOp invBQBt = Ip;

   LinearOp A = Thyra::block2x2(Iu,Ip,Iu,Zp);
   //RCP<Thyra::PreconditionerFactoryBase<double> > precFactory 
   //      = rcp(new LSCStablePreconditionerFactory(Iu,invBQBt));
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory 
         = rcp(new LSCPreconditionerFactory(Iu,invBQBt,Teuchos::null));
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,A);

   // build linear operator
   RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));
   // construct a couple of vectors
   Epetra_Vector ea(*map),eb(*map);
   Epetra_Vector ef(*map),eg(*map);
   const RCP<const Thyra::VectorBase<double> > x = BlockVector(ea,eb,A->domain());
   const RCP<const Thyra::VectorBase<double> > z = BlockVector(ef,eg,A->domain());
   const RCP<Thyra::VectorBase<double> > y = Thyra::createMember(A->range()); 

   // test vector [0 1 1 3]
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 1.0; eb[1] = 3.0;
   ef[0] = eb[0]; ef[1] = eb[1]; eg[0] = ea[0]-eb[0]; eg[1] = ea[1]-eb[1];
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = (PB::Test::Difference(y,z)<tolerance_);
   if(not status || verbosity>=10) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_Identity " << toString(status) << ": (y=inv(A)*x) != z" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",y);
   }
   allPassed &= status;

   // test vector [-2 4 7 9]
   ea[0] =-2.0; ea[1] = 4.0; eb[0] = 7.0; eb[1] = 9.0;
   ef[0] = eb[0]; ef[1] = eb[1]; eg[0] = ea[0]-eb[0]; eg[1] = ea[1]-eb[1];
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = (PB::Test::Difference(y,z)<tolerance_);
   if(not status || verbosity>=10) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_Identity " << toString(status) << ": (y=inv(A)*x) != z" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",y);
   }
   allPassed &= status;

   // test vector [1 0 0 -5]
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] =-5.0;
   ef[0] = eb[0]; ef[1] = eb[1]; eg[0] = ea[0]-eb[0]; eg[1] = ea[1]-eb[1];
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = (PB::Test::Difference(y,z)<tolerance_);
   if(not status || verbosity>=10) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_Identity " << toString(status) << ": (y=inv(A)*x) != z" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",y);
   }
   allPassed &= status;

   // test vector [4 -4 6 12]
   ea[0] = 4.0; ea[1] =-4.0; eb[0] = 6.0; eb[1] =12.0;
   ef[0] = eb[0]; ef[1] = eb[1]; eg[0] = ea[0]-eb[0]; eg[1] = ea[1]-eb[1];
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = (PB::Test::Difference(y,z)<tolerance_);
   if(not status || verbosity>=10) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_Identity " << toString(status) << ": (y=inv(A)*x) != z" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",y);
   }
   allPassed &= status;

   return allPassed;
}

bool tLSCStablePreconditionerFactory::test_diagonal(int verbosity,std::ostream & os)
{
   // make sure the preconditioner is working by testing against the identity matrix
   typedef RCP<const Thyra::VectorBase<double> > Vector;
   typedef RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;
   typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

   bool status = false;
   bool allPassed = true;
   double vec[2];
   double diff = 0.0;

   // build 4x4 matrix with block 2x2 diagonal subblocks
   //
   //            [ 1 0 7 0 ]
   // [ F G ] =  [ 0 2 0 8 ]
   // [ D C ]    [ 5 0 0 0 ]
   //            [ 0 6 0 0 ]
   //

   vec[0] = 1.0; vec[1] = 2.0;
   LinearOp F = PB::Test::DiagMatrix(2,vec);

   vec[0] = 7.0; vec[1] = 8.0;
   LinearOp G = PB::Test::DiagMatrix(2,vec);

   vec[0] = 5.0; vec[1] = 6.0;
   LinearOp D = PB::Test::DiagMatrix(2,vec);

   vec[0] = 0.0; vec[1] = 0.0;
   LinearOp C = PB::Test::DiagMatrix(2,vec);

   vec[0] = 1.0; vec[1] = 0.5;
   LinearOp iF = PB::Test::DiagMatrix(2,vec);

   // S = -C+D*iF*G
   vec[0] = 0.028571428571429; vec[1] = 0.020833333333333; 
   LinearOp iBBt = PB::Test::DiagMatrix(2,vec);

   LinearOp A = Thyra::block2x2(F,G,D,C);
   //RCP<Thyra::PreconditionerFactoryBase<double> > precFactory 
   //      = rcp(new LSCStablePreconditionerFactory(iF,iBBt));
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory 
         = rcp(new LSCPreconditionerFactory(iF,iBBt,Teuchos::null));
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,A);

   // build linear operator
   RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));
   // construct a couple of vectors
   Epetra_Vector ea(*map),eb(*map);
   Epetra_Vector ef(*map),eg(*map);
   const RCP<const Thyra::VectorBase<double> > x = BlockVector(ea,eb,A->domain());
   const RCP<const Thyra::VectorBase<double> > z = BlockVector(ef,eg,A->domain());
   const RCP<Thyra::VectorBase<double> > y = Thyra::createMember(A->range()); 

   // now checks of the preconditioner (should be exact!)
   /////////////////////////////////////////////////////////////////////////

   // test vector [0 1 1 3]
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 1.0; eb[1] = 3.0;
   ef[0] =  0.200000000000000; ef[1] =  0.500000000000000;  
   eg[0] = -0.028571428571429; eg[1] =  0;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_diagonal " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [-2 4 7 9]
   ea[0] =-2.0; ea[1] = 4.0; eb[0] = 7.0; eb[1] = 9.0;
   ef[0] =  1.400000000000000; ef[1] =  1.500000000000000;
   eg[0] = -0.485714285714286; eg[1] =  0.125000000000000;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_diagonal " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [1 0 0 -5]
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] =-5.0;
   ef[0] =  0.000000000000000; ef[1] = -0.833333333333333;
   eg[0] =  0.142857142857143; eg[1] =  0.208333333333333;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_diagonal " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [4 -4 6 12]
   ea[0] = 4.0; ea[1] =-4.0; eb[0] = 6.0; eb[1] =12.0;
   ef[0] =  1.200000000000000; ef[1] =  2.000000000000000;
   eg[0] =  0.400000000000000; eg[1] = -1.000000000000000;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_diagonal " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   return allPassed;
}

bool tLSCStablePreconditionerFactory::test_result(int verbosity,std::ostream & os)
{
   typedef RCP<const Thyra::VectorBase<double> > Vector;
   typedef RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;

   bool status = false;
   bool allPassed = true;
   double diff;
 
   // Build block2x2 preconditioner
   //RCP<Thyra::PreconditionerFactoryBase<double> > precFactory 
   //      = rcp(new LSCStablePreconditionerFactory(invF_,invBQBt_,invMass_));
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory 
         = rcp(new LSCPreconditionerFactory(invF_,invBQBt_,invMass_));
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,A_);

   // build linear operator
   RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));
   // construct a couple of vectors
   Epetra_Vector ea(*map),eb(*map);
   Epetra_Vector ef(*map),eg(*map);
   
   const RCP<const Thyra::VectorBase<double> > x = BlockVector(ea,eb,A_->domain());
   const RCP<const Thyra::VectorBase<double> > z = BlockVector(ef,eg,A_->domain());
   const RCP<Thyra::VectorBase<double> > y = Thyra::createMember(A_->range()); 

   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);

   // now checks of the preconditioner (should be exact!)
   /////////////////////////////////////////////////////////////////////////

   // test vector [0 1 1 3]
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 1.0; eb[1] = 3.0;
   ef[0] = -5.0000000000000000e0;  ef[1] = -1.999999999999994e0;
   eg[0] = -1.0999999999999991e1; eg[1] = -1.9999999999999979e1;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_result " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [-2 4 7 9]
   ea[0] =-2.0; ea[1] = 4.0; eb[0] = 7.0; eb[1] = 9.0;
   ef[0] = -1.700000000000000e+01; ef[1] = -7.999999999999977e+00;
   eg[0] = -3.849999999999996e+01; eg[1] = -6.949999999999991e+01;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_result " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [1 0 0 -5]
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] =-5.0;
   ef[0] = 7.500000000000000e+00; ef[1] = 2.499999999999991e+00;
   eg[0] = 1.449999999999999e+01; eg[1] = 2.599999999999997e+01;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_result " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [4 -4 6 12]
   ea[0] = 4.0; ea[1] =-4.0; eb[0] = 6.0; eb[1] =12.0;
   ef[0] = -2.100000000000000e+01; ef[1] = -8.999999999999975e+00;
   eg[0] = -4.499999999999996e+01; eg[1] = -8.799999999999991e+01;
   Thyra::apply(*precOp,NONCONJ_ELE,*x,&*y);
   status = ((diff = PB::Test::Difference(y,z)/Thyra::norm_2(*z))<tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCStablePreconditionerFactory::test_result " << toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " 
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   return allPassed;
}

} // end namespace Test
} // end namespace PB
