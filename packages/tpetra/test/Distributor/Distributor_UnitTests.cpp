#include "Teuchos_UnitTestHarness.hpp"

// #include "../tpetra_test_util.hpp"
#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream 

#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"

#include "Tpetra_Distributor.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::Distributor;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  template <typename T>
  T generateValue(T const x, T const y) {
    const T two = as<T>(2);
    // formula for z(x,y) = 0.5(x^2 + y^2 + 3x + y) + xy
    return(((x*x + y*y + x+x+x + y) / two) + (x*y));
  }

  // puts the generated values for (x, 0) ... (x, length-1) into vector
  template <typename T>
  void generateColumn(std::vector<T>& vector, const int x, const int length) {
    vector.resize(length);
    for(int y = 0; y < length; y++) {
      vector[y] = generateValue(as<T>(x), as<T>(y));
    }
  }


  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  template<class Ordinal>
  RCP<const Platform<Ordinal> > getDefaultPlatform()
  {
    if (testMpi) {
      return DefaultPlatform<Ordinal>::getPlatform();
    }
    return rcp(new Tpetra::SerialPlatform<Ordinal>());
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, basic, Ordinal )
  {
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    out << "platform = " << *platform << std::endl;
    TEST_INEQUALITY_CONST( platform->createComm(), Teuchos::null );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSends, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    out << "platform = " << *platform << std::endl;
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const Ordinal numImages = comm->getSize();
    const Ordinal ZERO = OT::zero();

    // send data to each image, including myself
    const Ordinal numExportIDs = as<Ordinal>(numImages); 
    Ordinal numRemoteIDs = ZERO;

    // fill exportImageIDs with {0, 1, 2, ... numImages-1}
    vector<Ordinal> exportImageIDs; 
    exportImageIDs.reserve(numExportIDs);
    for(Ordinal i = ZERO; i < numExportIDs; ++i) {
      exportImageIDs.push_back(i);
    }

    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numRemoteIDs);
    TEST_EQUALITY(numRemoteIDs, numImages);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromReceives, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    out << "platform = " << *platform << std::endl;
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getSize();
    const Ordinal ZERO = OT::zero();
    const Ordinal length = as<Ordinal>(numImages);
    const Ordinal invalid = as<Ordinal>(-2);

    Ordinal numRemoteIDs = length;

    // fill remoteGIDs with row from generator
    std::vector<Ordinal> remoteGIDs;
    remoteGIDs.reserve(numRemoteIDs);
    for(Ordinal i = ZERO; i < length; i++) {
      remoteGIDs.push_back( generateValue(i, as<Ordinal>(myImageID)) );
    }

    // fill remoteImageIDs with {0, 1, 2, ... length-1}
    vector<Ordinal> remoteImageIDs;
    remoteImageIDs.reserve(numRemoteIDs);
    for(Ordinal i = ZERO; i < numRemoteIDs; ++i) {
      remoteImageIDs.push_back(i);
    }

    Ordinal numExportIDs = ZERO;
    vector<Ordinal> exportGIDs(length,invalid),
                exportImageIDs(length,invalid);

    Distributor<Ordinal> distributor(comm);
    distributor.createFromRecvs(remoteGIDs, remoteImageIDs, exportGIDs, exportImageIDs);

    std::vector<Ordinal> expectedGIDs;
    generateColumn(expectedGIDs, myImageID, length);

    TEST_EQUALITY(numExportIDs, numRemoteIDs);
    TEST_COMPARE_ARRAYS(exportGIDs, expectedGIDs);
    TEST_COMPARE_ARRAYS(exportImageIDs, remoteImageIDs);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  //
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
# define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, basic, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSends, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, basic, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSends, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL )

    UNIT_TEST_GROUP_ORDINAL(char)
    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    UNIT_TEST_GROUP_ORDINAL_WITH_PAIRS(int)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL_WITH_PAIRS(LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}



/*

//======================================================================
template <typename Ordinal, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
  // ========================================
  // test createFromRecvs
  // ========================================
  // ========================================
  // test doPosts
  // ========================================
  if(verbose) cout << "Testing doPostsAndWaits... ";
  comm->barrier();

  std::vector<ScalarType> imports;
  std::vector<ScalarType> exports;
  generateColumn(exports, myImageID, numImages);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "exports: " + Tpetra::toString(exports));
  }
  // FINISH
  distributorS.doPostsAndWaits(exports, imports);
  if(debug) {
    outputData(myImageID, numImages, "imports: " + Tpetra::toString(imports));
    if(verbose) cout << "doPostsAndWaits test: ";
  }

  if(ierr != 0) {
    ierr = 1;
    if(verbose) cout << "failed" << endl;
  }
  else {
    if(verbose) cout << "passed" << endl;
  }
  returnierr += ierr;
  ierr = 0;
#endif // HAVE_MPI

  // ======================================================================
  // finish up
  // ======================================================================

  comm->barrier();
  if(verbose) {
    if(returnierr == 0) {
      outputHeading("Unit tests for " + className + " passed.");
    }
    else {
      outputHeading("Unit tests for " + className + " failed.");
    }
  }
  return(returnierr);
}
*/
