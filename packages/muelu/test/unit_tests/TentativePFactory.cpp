#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(TentativePFactory, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory> tentPFact = rcp(new TentativePFactory);
    TEST_EQUALITY(tentPFact != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(TentativePFactory, SetGetMethods)
  {
    out << "version: " << MueLu::Version() << std::endl;

    TentativePFactory tentPFact;

    bool flag = tentPFact.TentativeWithQR();
    TEST_EQUALITY(flag, false);
    tentPFact.TentativeWithQR(true);
    flag = tentPFact.TentativeWithQR();
    TEST_EQUALITY(flag, true);
  } //SetGetMethods

  //TODO test BuildP

  TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentative_LapackQR)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test QR with user-supplied nullspace" << std::endl;

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(199);
    fineLevel.Set("A",A);

      LO NSdim = 2;
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      nullSpace->randomize();
      fineLevel.Set("Nullspace",nullSpace);
      RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
      UCAggFact->SetMinNodesPerAggregate(3);
      UCAggFact->SetMaxNeighAlreadySelected(0);
      UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
      UCAggFact->SetPhase3AggCreation(0.5);

      RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory(UCAggFact));

      coarseLevel.Request("P",TentativePFact.get());         // request Ptent
      coarseLevel.Request("Nullspace",TentativePFact.get()); // request coarse nullspace (false, since no PR factory available)
      coarseLevel.Request(*TentativePFact);
      TentativePFact->Build(fineLevel,coarseLevel);

      RCP<Operator> Ptent; 
      coarseLevel.Get("P",Ptent,TentativePFact.get());

      RCP<MultiVector> coarseNullSpace; 
      coarseLevel.Get("Nullspace",coarseNullSpace, TentativePFact.get());

      //check interpolation
      RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(),NSdim);
      Ptent->apply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

      RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      diff->putScalar(0.0);

      coarseLevel.Release("P",TentativePFact.get()); // release Ptent
      coarseLevel.Release("Nullspace",TentativePFact.get());   // release coarse nullspace

      //diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
      diff->update(1.0,*nullSpace,-1.0,*PtN,0.0);

      Teuchos::Array<ST::magnitudeType> norms(NSdim);
      diff->norm2(norms);
      for (LO i=0; i<NSdim; ++i) {
        out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
        TEST_EQUALITY(norms[i]<1e-12, true);
      }

      Teuchos::ArrayRCP<const double> col1 = coarseNullSpace->getData(0);
      Teuchos::ArrayRCP<const double> col2 = coarseNullSpace->getData(1);
      TEST_EQUALITY(col1.size() == col2.size(), true);

  } //MakeTentative  Lapack QR

  TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentative)
  {
      out << "version: " << MueLu::Version() << std::endl;
      out << "Test QR with user-supplied nullspace" << std::endl;

      Level fineLevel, coarseLevel;
      TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

      RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(199);
      fineLevel.Request("A");
      fineLevel.Set("A",A);

      // only one NS vector -> exercises manual orthogonalization
      LO NSdim = 1;
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      nullSpace->randomize();
      fineLevel.Set("Nullspace",nullSpace);
      RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
      UCAggFact->SetMinNodesPerAggregate(3);
      UCAggFact->SetMaxNeighAlreadySelected(0);
      UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
      UCAggFact->SetPhase3AggCreation(0.5);

      //fineLevel.Request("Aggregates",&UCAggFact); //FIXME putting this in to avoid error until Merge needs business
      RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory(UCAggFact));

      coarseLevel.Request("P",TentativePFact.get());  // request Ptent
      coarseLevel.Request("Nullspace",TentativePFact.get());    // request coarse nullspace (false, since no PRFactory available)
      coarseLevel.Request(*TentativePFact);
      TentativePFact->Build(fineLevel,coarseLevel);

      RCP<Operator> Ptent;
      coarseLevel.Get("P",Ptent,TentativePFact.get());

      RCP<MultiVector> coarseNullSpace;
      coarseLevel.Get("Nullspace", coarseNullSpace, TentativePFact.get());

      //check interpolation
      RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(),NSdim);
      Ptent->apply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

      RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      diff->putScalar(0.0);

      coarseLevel.Release("P",TentativePFact.get()); // release Ptent
      coarseLevel.Release("Nullspace", TentativePFact.get());   // release coarse nullspace

      //diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
      diff->update(1.0,*nullSpace,-1.0,*PtN,0.0);

      Teuchos::Array<ST::magnitudeType> norms(NSdim);
      diff->norm2(norms);
      for (LO i=0; i<NSdim; ++i) {
        out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
        TEST_EQUALITY(norms[i]<1e-12, true);
      }

      // check normalization and orthogonality of prolongator columns
      Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > PtentTPtent = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(Ptent,true,Ptent,false);
      Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > diagVec = Xpetra::VectorFactory<Scalar,LO,GO>::Build(PtentTPtent->getRowMap());
      PtentTPtent->getLocalDiagCopy(*diagVec);
      //std::cout << diagVec->norm1() << " " << diagVec->normInf() << " " << diagVec->meanValue() << std::endl;
      TEST_EQUALITY(diagVec->norm1(), diagVec->getGlobalLength());
      TEST_EQUALITY(diagVec->normInf()-1 < 1e-12, true);
      TEST_EQUALITY(diagVec->meanValue(), 1.0);
      TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

  } //MakeTentative

  TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentativeUsingDefaultNullSpace)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test QR when nullspace isn't supplied by user" << std::endl;

    Level fineLevel, coarseLevel; TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(199);

    fineLevel.Set("A",A);

    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    UCAggFact->SetMinNodesPerAggregate(3);
    UCAggFact->SetMaxNeighAlreadySelected(0);
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    UCAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory(UCAggFact));

    coarseLevel.Request("P",tentativePFact.get());  // request Ptent
    coarseLevel.Request("Nullspace", tentativePFact.get());  // request coarse nullspace
    coarseLevel.Request(*tentativePFact);
    tentativePFact->Build(fineLevel,coarseLevel);

    RCP<Operator> Ptent; 
    coarseLevel.Get("P",Ptent,tentativePFact.get());

    RCP<MultiVector> coarseNullSpace; 
    coarseLevel.Get("Nullspace",coarseNullSpace, tentativePFact.get());

    coarseLevel.Release("P",tentativePFact.get()); // release Ptent
    coarseLevel.Release("Nullspace", tentativePFact.get());   // release coarse nullspace

    //grab default fine level nullspace (vector of all ones)
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), 1);
    nullSpace->putScalar(1.0);

    //check interpolation
    LO NSdim = 1;
    RCP<MultiVector> PtN = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    Ptent->apply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

    RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    diff->putScalar(0.0);

    //diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
    diff->update(1.0,*nullSpace,-1.0,*PtN,0.0);

    Teuchos::Array<ST::magnitudeType> norms(NSdim);
    diff->norm2(norms);
    for (LO i=0; i<NSdim; ++i) {
      out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
      TEST_EQUALITY(norms[i]<1e-12, true);
    }

  } //MakeTentative

  TEUCHOS_UNIT_TEST(TentativePFactory, TentativePFactory_NonStandardMaps)
  {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LO maxLevels = 3;
    //LO its=10;
    GO nEle = 63;
    GO nIndexBase = 10;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, nIndexBase, comm);


    RCP<CrsOperator> mtx = MueLu::Gallery::MatrixTraits<Map,CrsOperator>::Build(map, 3);

    LocalOrdinal NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
    GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
    assert(NumGlobalElements == nEle);

    GlobalOrdinal NumEntries;
    LocalOrdinal nnz=2;
    std::vector<Scalar> Values(nnz);
    std::vector<GlobalOrdinal> Indices(nnz);

    Scalar a = 2.0;
    Scalar b = -1.0;
    Scalar c = -1.0;

    for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i)
      {
        if (MyGlobalElements[i] == nIndexBase)
          {
            // off-diagonal for first row
            Indices[0] = nIndexBase;
            NumEntries = 1;
            Values[0] = c;
          }
        else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1)
          {
            // off-diagonal for last row
            Indices[0] = nIndexBase + NumGlobalElements - 2;
            NumEntries = 1;
            Values[0] = b;
          }
        else
          {
            // off-diagonal for internal row
            Indices[0] = MyGlobalElements[i] - 1;
            Values[1] = b;
            Indices[1] = MyGlobalElements[i] + 1;
            Values[0] = c;
            NumEntries = 2;
          }

        // put the off-diagonal entries
        // Xpetra wants ArrayViews (sigh)
        Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
        Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
        mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

        // Put in the diagonal entry
        mtx->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(a) );

      } //for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i)


    mtx->fillComplete(map,map);

    std::cout << map->getIndexBase() << std::endl;

    RCP<Operator> Op = Teuchos::rcp_dynamic_cast<Operator>(mtx);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    RCP<Level> Finest = H->GetLevel(); // first associate level with hierarchy (for defaultFactoryHandler!)

    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",Op);                      // set fine level matrix
    Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

    // define transfer operators
    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    UCAggFact->SetMinNodesPerAggregate(3);
    UCAggFact->SetMaxNeighAlreadySelected(0);
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    UCAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Pfact = rcp(new TentativePFactory(UCAggFact));
    RCP<RFactory>          Rfact = rcp( new TransPFactory() );
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    H->SetMaxCoarseSize(1);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LO) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother(lib, "RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    Teuchos::ParameterList status;
    status = H->FullPopulate(*Pfact,*Rfact, *Acfact,*SmooFact,0,maxLevels);

    RCP<Level> coarseLevel = H->GetLevel(1);
    TEST_EQUALITY(coarseLevel->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsKept("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsKept("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsKept("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsKept("PostSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsKept("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsKept("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsKept("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsKept("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsKept("R",Rfact.get()), false);
    RCP<Level> coarseLevel2 = H->GetLevel(2);
    TEST_EQUALITY(coarseLevel2->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsKept("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsKept("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsKept("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsKept("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsKept("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsKept("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsKept("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsKept("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsKept("R",Rfact.get()), false);

  }

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_TPETRA)
  TEUCHOS_UNIT_TEST(TentativePFactory, TentativePFactory_EpetraVsTpetra)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test QR when nullspace isn't supplied by user" << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::Array<ST::magnitudeType> results(2);

    // run test only on 1 proc
    if(comm->getSize() == 1)
    {
        Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

        // run Epetra and Tpetra test
        for (int run = 0; run < 2; run++)
        {
            if (run == 0) lib = Xpetra::UseEpetra;
            else lib = Xpetra::UseTpetra;

            // generate problem
            LO maxLevels = 3;
            LO its=10;
            LO nEle = 63;
            const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
            Teuchos::ParameterList matrixParameters;
            matrixParameters.set("nx",nEle);
            RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>("Laplace1D", map, matrixParameters);

            // build nullspace
            RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
            nullSpace->putScalar( (SC) 1.0);
            Teuchos::Array<ST::magnitudeType> norms(1);
            nullSpace->norm1(norms);
            if (comm->getRank() == 0)
              out << "||NS|| = " << norms[0] << std::endl;

            // fill hierarchy
            RCP<Hierarchy> H = rcp( new Hierarchy() );
            H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
            RCP<Level> Finest = H->GetLevel(); // first associate level with hierarchy (for defaultFactoryHandler!)

            Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
            Finest->Set("A",Op);                      // set fine level matrix
            Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

            // define transfer operators
            RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
            UCAggFact->SetMinNodesPerAggregate(3);
            UCAggFact->SetMaxNeighAlreadySelected(0);
            UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
            UCAggFact->SetPhase3AggCreation(0.5);

            RCP<TentativePFactory> Pfact = rcp(new TentativePFactory(UCAggFact));
            RCP<RFactory>          Rfact = rcp( new TransPFactory() );
            RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
            H->SetMaxCoarseSize(1);

            // setup smoothers
            Teuchos::ParameterList smootherParamList;
            smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
            smootherParamList.set("relaxation: sweeps", (LO) 1);
            smootherParamList.set("relaxation: damping factor", (SC) 1.0);
            RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother(lib, "RELAXATION", smootherParamList) );
            RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
            Acfact->setVerbLevel(Teuchos::VERB_HIGH);

            Teuchos::ParameterList status;
            status = H->FullPopulate(*Pfact,*Rfact, *Acfact,*SmooFact,0,maxLevels);

            SmootherFactory coarseSolveFact(smooProto);
            H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

            // test some basic multgrid data
            RCP<Level> coarseLevel = H->GetLevel(1);
            RCP<Operator> P1 = coarseLevel->Get< RCP<Operator> >("P");
            RCP<Operator> R1 = coarseLevel->Get< RCP<Operator> >("R");
            TEST_EQUALITY(P1->getGlobalNumRows(), 63);
            TEST_EQUALITY(P1->getGlobalNumCols(), 21);
            TEST_EQUALITY(R1->getGlobalNumRows(), 21);
            TEST_EQUALITY(R1->getGlobalNumCols(), 63);
            RCP<Level> coarseLevel2 = H->GetLevel(2);
            RCP<Operator> P2 = coarseLevel2->Get< RCP<Operator> >("P");
            RCP<Operator> R2 = coarseLevel2->Get< RCP<Operator> >("R");
            TEST_EQUALITY(P2->getGlobalNumRows(), 21);
            TEST_EQUALITY(P2->getGlobalNumCols(), 7);
            TEST_EQUALITY(R2->getGlobalNumRows(), 7);
            TEST_EQUALITY(R2->getGlobalNumCols(), 21);

            Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > PtentTPtent = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(P1,true,P1,false);
            Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > diagVec = Xpetra::VectorFactory<Scalar,LO,GO>::Build(PtentTPtent->getRowMap());
            PtentTPtent->getLocalDiagCopy(*diagVec);
            TEST_EQUALITY(diagVec->norm1()-diagVec->getGlobalLength() < 1e-12, true);
            TEST_EQUALITY(diagVec->normInf()-1 < 1e-12, true);
            TEST_EQUALITY(diagVec->meanValue()-1 < 1e-12, true);
            TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

            // Define RHS
            RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
            RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

            X->putScalar(1.0);
            X->norm2(norms);
            if (comm->getRank() == 0)
              out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

            Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

            // Use AMG directly as an iterative method
            {
              X->putScalar( (SC) 0.0);

              H->Iterate(*RHS,its,*X);

              X->norm2(norms);
              if (comm->getRank() == 0)
                out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
              results[run] = norms[0];
            }
        }

        TEST_EQUALITY(results[0], results[1]); // check results of EPETRA vs TPETRA
    } // comm->getSize == 1

  } // TentativePFactory_EpetraVsTpetra
#endif


} // namespace MueLuTests
