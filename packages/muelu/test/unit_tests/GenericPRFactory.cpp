#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(GenericPRFactory, Constructor_NoArgs)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<GenericPRFactory> PRFact = rcp(new GenericPRFactory());
    TEST_EQUALITY(PRFact != Teuchos::null, true);
  } //Constructor_NoArgs

  TEUCHOS_UNIT_TEST(GenericPRFactory, Constructor_TwoArgs)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory>  PFact = rcp(new SaPFactory());
    RCP<TransPFactory>  RFact = rcp(new TransPFactory());
    RCP<GenericPRFactory>  PRFact = rcp(new GenericPRFactory(PFact,RFact));
    TEST_EQUALITY(PRFact != Teuchos::null, true);
  } //Constructor_TwoArgs

  TEUCHOS_UNIT_TEST(GenericPRFactory, GetSetMethods)
  {

    out << "version: " << MueLu::Version() << std::endl;

    GenericPRFactory genericPR = GenericPRFactory();
    genericPR.ReUseAggregates(true);
    TEST_EQUALITY( genericPR.ReUseAggregates(), true);
    genericPR.ReUseAggregates(false);
    TEST_EQUALITY( genericPR.ReUseAggregates(), false);
    genericPR.ReUseGraph(true);
    TEST_EQUALITY( genericPR.ReUseGraph(), true);
    genericPR.ReUseGraph(false);
    TEST_EQUALITY( genericPR.ReUseGraph(), false);
  } //GetSetMethods

  TEUCHOS_UNIT_TEST(GenericPRFactory, TooCoarse_DoNotBuild)
  {
    bool r;

    out << "version: " << MueLu::Version() << std::endl;
    out << "Test that Build returns early if the coarsest matrix is smaller than specified MaxCoarseSize" << std::endl;

    // build test-specific default factory handler
    RCP<DefaultFactoryHandlerBase> defHandler = rcp(new DefaultFactoryHandlerBase());
    defHandler->SetDefaultFactory("A", rcp(MueLu::NoFactory::get(),false));         // dummy factory for A
    defHandler->SetDefaultFactory("Nullspace", rcp(new NullspaceFactory()));        // real null space factory for Ptent
    defHandler->SetDefaultFactory("Graph", rcp(new CoalesceDropFactory()));         // real graph factory for Ptent
    defHandler->SetDefaultFactory("Aggregates", rcp(new UCAggregationFactory()));   // real aggregation factory for Ptent

    Level fineLevel, coarseLevel; TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel, defHandler);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(240);
    fineLevel.Request("A",NULL);
    fineLevel.Set("A",A,NULL);

    GenericPRFactory genericPR = GenericPRFactory();
    genericPR.SetMaxCoarseSize(500);
    genericPR.DeclareInput(fineLevel,coarseLevel);
    r = genericPR.Build(fineLevel,coarseLevel);
    TEST_EQUALITY( r , false);

    out << "Test that Build completes if the coarsest matrix is larger than specified MaxCoarseSize" << std::endl;
    genericPR.SetMaxCoarseSize(239);
    r = genericPR.Build(fineLevel,coarseLevel);
    TEST_EQUALITY( r , true);
  } //TooCoarse_DoNotBuild

}//namespace MueLuTests

