
#include "Thyra_EpetraExtAddTransformer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Thyra::EpetraExtAddTransformer;
using Thyra::epetraExtAddTransformer;
using Thyra::VectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::LinearOpTester;
using Thyra::adjoint;
using Thyra::multiply;
using Thyra::diagonal;


std::string matrixFile = "";
std::string matrixFile2 = "";


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "matrix-file", &matrixFile,
    "Defines the Epetra_CrsMatrix to read in."  );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "matrix-file-2", &matrixFile2,
    "Defines the Epetra_CrsMatrix to read in."  );
}

const Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildAddOperator(int scenario,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A,
                               const Teuchos::RCP<const Thyra::LinearOpBase<double> > & B)
{
   // build operators for the various addition/adjoint scenarios
   RCP<const Thyra::LinearOpBase<double> > M;
   
   switch(scenario) {
   case 0:
      M = Thyra::add(A,B,"A+B");
      break;
   case 1:
      M = Thyra::add(A,B,"A+adj(B)");
      break;
   case 2:
      M = Thyra::add(A,B,"adj(A)+B");
      break;
   case 3:
      M = Thyra::add(A,B,"adb(A)+adb(B)");
      break;
   default:
      TEUCHOS_ASSERT(false);
      break;
   }

   return M;
}

TEUCHOS_UNIT_TEST( EpetraExtAddTransformer, basic_Add )
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\' ...\n";
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile2<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_A;
  RCP<Epetra_CrsMatrix> epetra_B;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A, NULL, NULL, NULL );
  EpetraExt::readEpetraLinearSystem( matrixFile2, comm, &epetra_B, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
  double scaleA=3.7;
  double scaleB=-2.9;
 
  const RCP<const Thyra::LinearOpBase<double> > A = Thyra::scale<double>(scaleA,Thyra::epetraLinearOp(epetra_B));
  const RCP<const Thyra::LinearOpBase<double> > B = Thyra::scale<double>(scaleB,Thyra::epetraLinearOp(epetra_B));

  out << "\nA = " << *A;
  out << "\nB = " << *B;

  for(int scenario=0;scenario<4;scenario++) {
     //
     // C) Create implicit A+B operator
     //
   
     const RCP<const Thyra::LinearOpBase<double> > M = buildAddOperator(scenario,A,B);
   
     //
     // D) Do the transformation
     //
   
     const RCP<EpetraExtAddTransformer> ApB_transformer = epetraExtAddTransformer();
   
     TEST_ASSERT(ApB_transformer != null);
   
     const RCP<LinearOpBase<double> > M_explicit = ApB_transformer->createOutputOp();
     ApB_transformer->transform( *M, M_explicit.ptr() );
   
     out << "\nM_explicit = " << *M_explicit;
   
     //
     // E) Check the explicit operator
     //
   
     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;
   
     const bool result = M_explicit_tester.compare( *M, *M_explicit, &out );
     if (!result) success = false;
  }
}

} // end namespace
