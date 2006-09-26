
#include "RBGen_Params.h"
#include "RBGen_Utils.h"
#include "RBGen_EpetraMVFileIOFactory.h"
#include "RBGen_EpetraMVMethodFactory.h"
#include "RBGen_EpetraMVPreprocessorFactory.h"
#include "RBGen_PODMethod.hpp"

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_LAPACK.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_Array.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"
	
int main( int argc, char* argv[] )
{

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::Array<Teuchos::RefCountPtr<Teuchos::Time> > timersRBGen;
  //
  // ---------------------------------------------------------------
  //  CREATE THE INITIAL PARAMETER LIST FROM THE COMMAND LINE OPTIONS
  // ---------------------------------------------------------------
  //
  Teuchos::RefCountPtr<Teuchos::ParameterList> BasisParams = RBGen::createParams( argc, argv );
  //
  // ---------------------------------------------------------------
  //  CREATE THE FILE I/O HANDLER
  // ---------------------------------------------------------------
  //
  //  - First create the abstract factory for the file i/o handler.
  //
  RBGen::EpetraMVFileIOFactory fio_factory;
  //
  //  - Then use the abstract factory to create the file i/o handler specified in the parameter list.
  //
  Teuchos::RefCountPtr<Teuchos::Time> timerFileIO = Teuchos::rcp( new Teuchos::Time("Create File I/O Handler") );
  timersRBGen.push_back( timerFileIO );
  //
  Teuchos::RefCountPtr< RBGen::FileIOHandler<Epetra_MultiVector> > fileio;
  {
    Teuchos::TimeMonitor lcltimer( *timerFileIO );
    fileio = fio_factory.create( *BasisParams );
    //					    
    // Initialize file IO handler
    //
    fileio->Initialize( BasisParams );
  }    
  //
  // ---------------------------------------------------------------
  //  RETRIEVE THE DATA SET / SNAPSHOT SET & PREPROCESS
  //  ( this will be a separate abstract class type )
  // ---------------------------------------------------------------
  //
  std::vector<std::string> filenames = Teuchos::getParameter<std::vector<std::string> >( *BasisParams, "Snapshot Files");
  Teuchos::RefCountPtr<Teuchos::Time> timerSnapshotIn = Teuchos::rcp( new Teuchos::Time("Reading in Snapshot Set") );
  timersRBGen.push_back( timerSnapshotIn );
  //
  Teuchos::RefCountPtr<Epetra_MultiVector> testMV;
  {
    Teuchos::TimeMonitor lcltimer( *timerSnapshotIn );
    testMV = fileio->Read( filenames );
  } 

  RBGen::EpetraMVPreprocessorFactory preprocess_factory;

  Teuchos::RefCountPtr<Teuchos::Time> timerCreatePreprocessor = Teuchos::rcp( new Teuchos::Time("Create Preprocessor") );
  timersRBGen.push_back( timerCreatePreprocessor );
  Teuchos::RefCountPtr<RBGen::Preprocessor<Epetra_MultiVector> > prep;
  {
    Teuchos::TimeMonitor lcltimer( *timerCreatePreprocessor );
    prep = preprocess_factory.create( *BasisParams );
    //
    // Initialize preprocessor.
    //
    prep->Initialize( BasisParams, fileio );
  }

  Teuchos::RefCountPtr<Teuchos::Time> timerPreprocess = Teuchos::rcp( new Teuchos::Time("Preprocess Snapshot Set") );  
  timersRBGen.push_back( timerPreprocess );
  {
    Teuchos::TimeMonitor lcltimer( *timerPreprocess );
    prep->Preprocess( testMV );
  }
  //
  // ---------------------------------------------------------------
  //  COMPUTE THE REDUCED BASIS
  // ---------------------------------------------------------------
  //
  //  - First create the abstract factory for the reduced basis methods.
  //
  RBGen::EpetraMVMethodFactory mthd_factory;
  //
  //  - Then use the abstract factory to create the method specified in the parameter list.
  //
  Teuchos::RefCountPtr<Teuchos::Time> timerCreateMethod = Teuchos::rcp( new Teuchos::Time("Create Reduced Basis Method") );
  timersRBGen.push_back( timerCreateMethod );
  Teuchos::RefCountPtr<RBGen::Method<Epetra_MultiVector> > method;
  {
    Teuchos::TimeMonitor lcltimer( *timerCreateMethod );  
    method = mthd_factory.create( *BasisParams );
    //
    // Initialize reduced basis method.
    //
    method->Initialize( BasisParams, testMV );
  }
  //
  //  - Call the computeBasis method on the reduced basis method object.
  //
  Teuchos::RefCountPtr<Teuchos::Time> timerComputeBasis = Teuchos::rcp( new Teuchos::Time("Reduced Basis Computation") );
  timersRBGen.push_back( timerComputeBasis );
  {
    Teuchos::TimeMonitor lcltimer( *timerComputeBasis );  
    method->computeBasis();
  }
  //
  //  - Retrieve the computed basis from the method object.
  //
  Teuchos::RefCountPtr<const Epetra_MultiVector> basisMV = method->getBasis();
  //
  //  Since we're using a POD method, we can dynamic cast to get the singular values.
  //
  Teuchos::RefCountPtr<RBGen::PODMethod<double> > pod_method = Teuchos::rcp_dynamic_cast<RBGen::PODMethod<double> >( method );
  const std::vector<double> sv = pod_method->getSingularValues();
  //
  /*if (Comm.MyPID() == 0) {
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"Computed Singular Values : "<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    for (unsigned int i=0; i<sv.size(); ++i) { cout << sv[i] << endl; }
  }      
  */
  if (Comm.MyPID() == 0) {
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"RBGen Computation Time Breakdown (seconds) : "<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    for (unsigned int i=0; i<timersRBGen.size(); ++i)
      cout << std::left << std::setw(40) << timersRBGen[i]->name() << " : "
	   << std::setw(15) << timersRBGen[i]->totalElapsedTime() << endl;
    cout<<"-------------------------------------------------------"<<endl;
  }
  //
  // ---------------------------------------------------------------
  //  POSTPROCESS BASIS (not necessary right now)
  // ---------------------------------------------------------------
  //
  //
  // ---------------------------------------------------------------
  //  WRITE OUT THE REDUCED BASIS
  // ---------------------------------------------------------------
  //
  if ( BasisParams->isParameter( "Reduced Basis Output File" ) ) {
    std::string outfile = Teuchos::getParameter<std::string>( *BasisParams, "Reduced Basis Output File" );
    fileio->Write( basisMV, outfile );
  }
  //
#ifdef EPETRA_MPI
  // Finalize MPI
  MPI_Finalize();
#endif

  return 0;
}


