// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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

#include <cstdlib>
#include <fstream>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <XpetraExt_MatrixMatrix.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_Operator.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#endif
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <MueLu_CreateEpetraPreconditioner.hpp>
#endif

#include <MueLu_UseDefaultTypes.hpp>

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"

#ifdef HAVE_MUELU_EPETRA
#include "BelosEpetraAdapter.hpp"
#endif
#endif

void run_sed(const std::string& pattern, const std::string& baseFile);

const std::string thickSeparator = "==========================================================================================================================";
const std::string thinSeparator  = "--------------------------------------------------------------------------------------------------------------------------";

const std::string prefSeparator = "=====================================";

namespace MueLuExamples {
#include <MueLu_UseShortNames.hpp>
// --------------------------------------------------------------------------------------
  void solve_system_list(Xpetra::UnderlyingLib & lib, Teuchos::RCP<Matrix> & A, Teuchos::RCP<Vector>&  X, Teuchos::RCP<Vector> & B, Teuchos::ParameterList & MueLuList, Teuchos::RCP<Teuchos::ParameterList> & SList, const std::string & fname) {
    using Teuchos::RCP;
    using Teuchos::rcp;
    int myRank  = A->getRowMap()->getComm()->getRank();

    std::filebuf    buffer;
    std::streambuf* oldbuffer = NULL;


#ifdef HAVE_MUELU_BELOS
#ifdef HAVE_MUELU_TPETRA
    typedef Tpetra::Operator<SC,LO,GO> Tpetra_Operator;
    typedef Tpetra::CrsMatrix<SC,LO,GO> Tpetra_CrsMatrix;
    typedef Tpetra::Vector<SC,LO,GO> Tpetra_Vector;
    typedef Tpetra::MultiVector<SC,LO,GO> Tpetra_MultiVector;
    if(lib==Xpetra::UseTpetra) {
      RCP<Tpetra_CrsMatrix>   At = Xpetra::MatrixMatrix::Op2NonConstTpetraCrs(A);

      // Redirect output
      if (myRank == 0) {
	buffer.open((fname + ".out").c_str(), std::ios::out);
	oldbuffer = std::cout.rdbuf(&buffer);
      }
      RCP<Tpetra_Operator>    Mt = MueLu::CreateTpetraPreconditioner(At,MueLuList);
      // Redirect output back
      if(myRank==0) {
	std::cout.rdbuf(oldbuffer);
	buffer.close();        
      }
      
      RCP<Tpetra_MultiVector> Xt = Xpetra::toTpetra(*X);
      RCP<Tpetra_MultiVector> Bt = Xpetra::toTpetra(*B);
      typedef Tpetra_MultiVector MV;
      typedef Tpetra_Operator OP;
      RCP<Belos::LinearProblem<SC,MV,OP> > belosProblem = rcp(new Belos::LinearProblem<SC,MV,OP>(At, Xt, Bt));
      belosProblem->setRightPrec(Mt);
      belosProblem->setProblem(Xt,Bt);

      Belos::SolverFactory<SC, MV, OP> BelosFactory;
      Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > BelosSolver = BelosFactory.create(std::string("CG"), SList);
      BelosSolver->setProblem(belosProblem);
      Belos::ReturnType result = BelosSolver->solve();
      if(result==Belos::Unconverged)
        throw std::runtime_error("Belos failed to converge");
    }
#endif
#ifdef HAVE_MUELU_EPETRA
    if(lib==Xpetra::UseEpetra) {
      RCP<Epetra_CrsMatrix>   Ae = Xpetra::MatrixMatrix::Op2NonConstEpetraCrs(A);

      // Redirect output
      if (myRank == 0) {
	buffer.open((fname + ".out").c_str(), std::ios::out);
	oldbuffer = std::cout.rdbuf(&buffer);
      }
      RCP<Epetra_Operator>    Me = MueLu::CreateEpetraPreconditioner(Ae,MueLuList);
      // Redirect output back
      if(myRank==0) {
	std::cout.rdbuf(oldbuffer);
	buffer.close();        
      }
      RCP<Epetra_MultiVector> Xe = rcp(&Xpetra::toEpetra(*X),false);
      RCP<Epetra_MultiVector> Be = rcp(&Xpetra::toEpetra(*B),false);
      typedef Epetra_MultiVector MV;
      typedef Epetra_Operator OP;
      RCP<Belos::LinearProblem<SC,MV,OP> > belosProblem = rcp(new Belos::LinearProblem<SC,MV,OP>(Ae, Xe, Be));
      Teuchos::RCP<Belos::EpetraPrecOp> PrecWrap = Teuchos::rcp(new Belos::EpetraPrecOp(Me));
      belosProblem->setRightPrec(PrecWrap);
      belosProblem->setProblem(Xe,Be);

      Belos::SolverFactory<SC, MV, OP> BelosFactory;
      Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > BelosSolver = BelosFactory.create(std::string("CG"), SList);
      BelosSolver->setProblem(belosProblem);
      Belos::ReturnType result = BelosSolver->solve();
      if(result==Belos::Unconverged)
        throw std::runtime_error("Belos failed to converge");
    }
#endif
#endif // #ifdef HAVE_MUELU_BELOS


  }

// --------------------------------------------------------------------------------------
void solve_system_hierarchy(Xpetra::UnderlyingLib & lib, Teuchos::RCP<Matrix> & A, Teuchos::RCP<Vector>&  X, Teuchos::RCP<Vector> & B, Teuchos::RCP<Hierarchy> & H, Teuchos::RCP<Teuchos::ParameterList> & SList) {
    using Teuchos::RCP;
    using Teuchos::rcp;
    
#ifdef HAVE_MUELU_BELOS
#ifdef HAVE_MUELU_TPETRA
    typedef Tpetra::Operator<SC,LO,GO> Tpetra_Operator;
    typedef Tpetra::CrsMatrix<SC,LO,GO> Tpetra_CrsMatrix;
    typedef Tpetra::Vector<SC,LO,GO> Tpetra_Vector;
    typedef Tpetra::MultiVector<SC,LO,GO> Tpetra_MultiVector;
    if(lib==Xpetra::UseTpetra) {
      RCP<Tpetra_CrsMatrix>   At = Xpetra::MatrixMatrix::Op2NonConstTpetraCrs(A);
      RCP<Tpetra_Operator>    Mt = rcp(new MueLu::TpetraOperator<SC,LO,GO>(H));
      RCP<Tpetra_MultiVector> Xt = Xpetra::toTpetra(*X);
      RCP<Tpetra_MultiVector> Bt = Xpetra::toTpetra(*B);
      typedef Tpetra_MultiVector MV;
      typedef Tpetra_Operator OP;
      RCP<Belos::LinearProblem<SC,MV,OP> > belosProblem = rcp(new Belos::LinearProblem<SC,MV,OP>(At, Xt, Bt));
      belosProblem->setRightPrec(Mt);
      belosProblem->setProblem(Xt,Bt);

      Belos::SolverFactory<SC, MV, OP> BelosFactory;
      Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > BelosSolver = BelosFactory.create(std::string("CG"), SList);
      BelosSolver->setProblem(belosProblem);
      Belos::ReturnType result = BelosSolver->solve();
      if(result==Belos::Unconverged)
        throw std::runtime_error("Belos failed to converge");
    }
#endif
#ifdef HAVE_MUELU_EPETRA
    if(lib==Xpetra::UseEpetra) {
      RCP<Epetra_CrsMatrix>      Ae = Xpetra::MatrixMatrix::Op2NonConstEpetraCrs(A);
      RCP<MueLu::EpetraOperator> Me = rcp(new MueLu::EpetraOperator(H));
      RCP<Epetra_MultiVector>    Xe = rcp(&Xpetra::toEpetra(*X),false);
      RCP<Epetra_MultiVector>    Be = rcp(&Xpetra::toEpetra(*B),false);
      typedef Epetra_MultiVector MV;
      typedef Epetra_Operator OP;
      RCP<Belos::LinearProblem<SC,MV,OP> > belosProblem = rcp(new Belos::LinearProblem<SC,MV,OP>(Ae, Xe, Be));
      Teuchos::RCP<Belos::EpetraPrecOp> PrecWrap = Teuchos::rcp(new Belos::EpetraPrecOp(Me));
      belosProblem->setRightPrec(PrecWrap);
      belosProblem->setProblem(Xe,Be);

      Belos::SolverFactory<SC, MV, OP> BelosFactory;
      Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > BelosSolver = BelosFactory.create(std::string("CG"), SList);
      BelosSolver->setProblem(belosProblem);
      Belos::ReturnType result = BelosSolver->solve();
      if(result==Belos::Unconverged)
        throw std::runtime_error("Belos failed to converge");
    }
#endif
#endif // #ifdef HAVE_MUELU_BELOS

  }

// --------------------------------------------------------------------------------------
// This routine generate's the user's original A matrix and nullspace
void generate_user_matrix_and_nullspace(std::string &matrixType,  Xpetra::UnderlyingLib & lib,Teuchos::ParameterList &galeriList,  Teuchos::RCP<const Teuchos::Comm<int> > &comm, Teuchos::RCP<Matrix> & A,Teuchos::RCP<MultiVector> & nullspace){
    using Teuchos::RCP;

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *fancy;

    RCP<const Map>   map;
    RCP<MultiVector> coordinates;
    if (matrixType == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" || matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D", map, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, (matrixType == "Elasticity2D" ? 2 : 3));

    out << "Processor subdomains in x direction: " << galeriList.get<int>("mx") << std::endl
        << "Processor subdomains in y direction: " << galeriList.get<int>("my") << std::endl
        << "Processor subdomains in z direction: " << galeriList.get<int>("mz") << std::endl
        << "========================================================" << std::endl;

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixType, map, galeriList);

    A = Pr->BuildMatrix();

    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      nullspace = Pr->BuildNullspace();
      A->SetFixedBlockSize((matrixType == "Elasticity2D") ? 2 : 3);
    }    
}



bool compare_to_gold(int myRank, const std::string & baseFile) { 
  bool failed=false;
  if (myRank == 0) {
    
    // Create a copy of outputs
    std::string cmd = "cp -f ";
    system((cmd + baseFile + ".res "  + baseFile + ".resorig").c_str());
    system((cmd + baseFile + ".out "  + baseFile + ".outorig").c_str());
    
    // Tpetra produces different eigenvalues in Chebyshev due to using
    // std::rand() for generating random vectors, which may be initialized
    // using different seed, and may have different algorithm from one
    // gcc version to another, or to anogther compiler (like clang)
    // This leads to us always failing this test.
    // NOTE1 : Epetra, on the other hand, rolls out its out random number
    // generator, which always produces same results
    
    // Ignore the value of "lambdaMax"
    run_sed("'s/lambdaMax: [0-9]*.[0-9]*/lambdaMax = <ignored>/'", baseFile);
    
    // Ignore the value of "lambdaMin"
    run_sed("'s/lambdaMin: [0-9]*.[0-9]*/lambdaMin = <ignored>/'", baseFile);
    
    // Ignore the value of "chebyshev: max eigenvalue"
    // NOTE: we skip lines with default value ([default])
    run_sed("'/[default]/! s/chebyshev: max eigenvalue = [0-9]*.[0-9]*/chebyshev: max eigenvalue = <ignored>/'", baseFile);
    
    // Ignore the exact type of direct solver (it is selected semi-automatically
    // depending on how Trilinos was configured
    run_sed("'s/Amesos\\([2]*\\)Smoother{type = .*}/Amesos\\1Smoother{type = <ignored>}/'", baseFile);
    run_sed("'s/SuperLU solver interface, direct solve/<Direct> solver interface/'", baseFile);
    run_sed("'s/KLU2 solver interface/<Direct> solver interface/'", baseFile);
    run_sed("'s/Basker solver interface/<Direct> solver interface/'", baseFile);
    
    // Nuke all pointers
    run_sed("'s/0x[0-9a-f]*//g'",baseFile);

    // Strip template args for some classes
    std::vector<std::string> classes;
    classes.push_back("Xpetra::Matrix");
    classes.push_back("MueLu::Constraint");
    for (size_t q = 0; q < classes.size(); q++)
      run_sed("'s/" + classes[q] + "<.*>/" + classes[q] + "<ignored> >/'", baseFile);
    
#ifdef __APPLE__
    // Some Macs print outs ptrs as 0x0 instead of 0, fix that
    run_sed("'/RCP/ s/=0x0/=0/g'", baseFile);
#endif
    
    // Run comparison (ignoring whitespaces)
    cmd = "diff -u -w -I\"^\\s*$\" " + baseFile + ".res " + baseFile + ".out";
    int ret = system(cmd.c_str());
    if (ret)
      failed = true;
    
    std::cout << baseFile << ": " << (ret ? "failed" : "passed") << std::endl;
  }

  return !failed;
}


void run_sed(const std::string& pattern, const std::string& baseFile) {
  // sed behaviour differs between Mac and Linux
  // You can run "sed -i 's//' " in Linux, but you always have to specify
  // "sed -i "<smth,could be empty>" 's//'" in Mac. Both, however, take '-i<extension>'
  std::string sed_pref = "sed -i ";
#ifdef __APPLE__
  sed_pref = sed_pref +  "\"\" ";
#endif

  system((sed_pref + pattern + " " + baseFile + ".res").c_str());
  system((sed_pref + pattern + " " + baseFile + ".out").c_str());
}


}//namespace

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = true;
  //  bool verbose = true;
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int numProc = comm->getSize();
  int myRank  = comm->getRank();
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  
  // =========================================================================
  // Parameters initialization
  // =========================================================================
  Teuchos::CommandLineProcessor clp(false);
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  ::Xpetra::Parameters xpetraParameters(clp);
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();
  ParameterList galeriList = galeriParameters.GetParameterList();
  out << thickSeparator << std::endl << xpetraParameters << galeriParameters;
  
  // =========================================================================
  // Problem construction
  // =========================================================================
  RCP<const Map>   map;
  RCP<Matrix> A,P,R, Ac;
  RCP<MultiVector> nullspace;
  std::string matrixType = galeriParameters.GetMatrixType();
  MueLuExamples::generate_user_matrix_and_nullspace(matrixType,lib,galeriList,comm,A,nullspace);
  map=A->getRowMap();
  
  // =========================================================================
  // Setups and solves
  // =========================================================================
  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);
  B->setSeed(846930886);
  B->randomize();
  RCP<TimeMonitor> tm;
  
#ifdef HAVE_MUELU_BELOS
  // Belos Options
  RCP<Teuchos::ParameterList> SList = rcp(new Teuchos::ParameterList );
  SList->set("Verbosity",Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  SList->set("Output Frequency",10);
  SList->set("Output Style",Belos::Brief);
  SList->set("Maximum Iterations",200);
  SList->set("Convergence Tolerance",1e-10);
#endif
  
  
  // =========================================================================
  // Solve #1 (standard MueLu)
  // =========================================================================
  out << thickSeparator << std::endl;
  out << prefSeparator << " Solve 1: Standard "<< prefSeparator <<std::endl;
  {
    std::string fname = "Output/operator_solve_1_np" + Teuchos::toString(numProc);
    std::srand(12345);

    Teuchos::ParameterList MueLuList;
    MueLuList.set("verbosity","high");
    MueLuList.set("coarse: max size",100);

    std::filebuf    buffer;
    std::streambuf* oldbuffer = NULL;
    if (myRank == 0) {
      // Redirect output
      buffer.open((fname + ".out").c_str(), std::ios::out);
      oldbuffer = std::cout.rdbuf(&buffer);
    }
    
    ParameterListInterpreter mueLuFactory(MueLuList);
    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    Teuchos::RCP<FactoryManagerBase> LevelFactory = mueLuFactory.GetFactoryManager(1);
    H->setlib(lib);
    H->AddNewLevel();
    H->GetLevel(1)->Keep("Nullspace",LevelFactory->GetFactory("Nullspace").get());
    H->GetLevel(0)->Set("A", A);
    mueLuFactory.SetupHierarchy(*H);
    
    // Redirect output back
    if(myRank==0) {
      std::cout.rdbuf(oldbuffer);
      buffer.close();        
    }

#ifdef HAVE_MUELU_BELOS
    // Solve
    MueLuExamples::solve_system_hierarchy(lib,A,X,B,H,SList);
    success = success && MueLuExamples::compare_to_gold(myRank,fname);
#endif
    
    // Extract R,P & Ac for LevelWrap Usage
    H->GetLevel(1)->Get("R",R);
    H->GetLevel(1)->Get("P",P);
    H->GetLevel(1)->Get("A",Ac);    
    nullspace = H->GetLevel(1)->Get<RCP<MultiVector> >("Nullspace",LevelFactory->GetFactory("Nullspace").get());
  }
  

  // =========================================================================
  // Solve #5 (level wrap, the fast way, P, R + Nullspace)
  // =========================================================================
  out << thickSeparator << std::endl;
  out << prefSeparator << " Solve 5: LevelWrap, Fast Way, P, R "<< prefSeparator <<std::endl;
  {
    std::string fname = "Output/operator_solve_5_np" + Teuchos::toString(numProc);
    std::srand(12345);

    Teuchos::ParameterList MueLuList, level1;
    level1.set("R",R);
    level1.set("P",P);
    level1.set("Nullspace",nullspace);
    MueLuList.set("level 1",level1);
    MueLuList.set("verbosity","high");
    MueLuList.set("coarse: max size",100);
#ifdef HAVE_MUELU_BELOS
    MueLuExamples::solve_system_list(lib,A,X,B,MueLuList,SList,fname);
    success = success && MueLuExamples::compare_to_gold(myRank,fname);
#endif
    }


  // =========================================================================
  // Solve #6 (level wrap, the fast way, P only, explicit transpose)
  // =========================================================================
  out << thickSeparator << std::endl;
  out << prefSeparator << " Solve 6: LevelWrap, Fast Way, P only, explicit transpose "<< prefSeparator <<std::endl;
  {
    std::string fname = "Output/operator_solve_6_np" + Teuchos::toString(numProc);
    std::srand(12345);

    Teuchos::ParameterList MueLuList, level1;
    level1.set("P",P);
    level1.set("Nullspace",nullspace);
    MueLuList.set("level 1",level1);
    MueLuList.set("verbosity","high");
    MueLuList.set("coarse: max size",100);
    MueLuList.set("transpose: use implicit",false);
    MueLuList.set("max levels",4);
#ifdef HAVE_MUELU_BELOS
    MueLuExamples::solve_system_list(lib,A,X,B,MueLuList,SList,fname);
    success = success && MueLuExamples::compare_to_gold(myRank,fname);
#endif
  }


  // =========================================================================
  // Solve #7 (level wrap, the fast way, P only, implicit transpose)
  // =========================================================================
  out << thickSeparator << std::endl;
  out << prefSeparator << " Solve 7: LevelWrap, Fast Way, P only, implicit transpose "<< prefSeparator <<std::endl;
  {
    std::string fname = "Output/operator_solve_7_np" + Teuchos::toString(numProc);
    std::srand(12345);

    Teuchos::ParameterList MueLuList, level1;
    level1.set("P",P);
    level1.set("Nullspace",nullspace);
    MueLuList.set("level 1",level1);
    MueLuList.set("verbosity","high");
    MueLuList.set("coarse: max size",100);
    MueLuList.set("transpose: use implicit",true);
    MueLuList.set("max levels",2);
#ifdef HAVE_MUELU_BELOS
    MueLuExamples::solve_system_list(lib,A,X,B,MueLuList,SList,fname);
    success = success && MueLuExamples::compare_to_gold(myRank,fname);
#endif
  }

  
  if (myRank == 0)
    std::cout << std::endl << "End Result: TEST " << (!success ? "FAILED" : "PASSED") << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

void run_sed(const std::string& pattern, const std::string& baseFile) {
  // sed behaviour differs between Mac and Linux
  // You can run "sed -i 's//' " in Linux, but you always have to specify
  // "sed -i "<smth,could be empty>" 's//'" in Mac. Both, however, take '-i<extension>'
  std::string sed_pref = "sed -i ";
#ifdef __APPLE__
  sed_pref = sed_pref +  "\"\" ";
#endif

  system((sed_pref + pattern + " " + baseFile + ".res").c_str());
  system((sed_pref + pattern + " " + baseFile + ".out").c_str());
}
