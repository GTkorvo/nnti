#include "Teuchos_UnitTestHarness.hpp"


#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


namespace {


using Teuchos::as;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::GlobalMPISession;
using Teuchos::defaultSmallNumber;


bool testMpi = true;


double errorTolSlack = 1e+1;


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
RCP<const Comm<Ordinal> > getDefaultComm()
{
  if (testMpi) {
    return DefaultComm<Ordinal>::getComm();
  }
  return rcp(new Teuchos::SerialComm<Ordinal>);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMpiComm, basic, Ordinal )
{
  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  out << "comm = " << Teuchos::describe(*comm);
  TEST_EQUALITY( size(*comm), GlobalMPISession::getNProc() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, reduceAllAndScatter_1, Ordinal, Packet )
{

  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  //const Ordinal procRank = rank(*comm);

#ifdef TEUCHOS_MPI_COMM_DUMP
  Teuchos::MpiComm<Ordinal>::show_dump = true;
#endif

  Array<Packet> sendBuffer(as<Ordinal>(numProcs));
  for (Ordinal k = 0; k < numProcs; ++k) {
    sendBuffer[k] = as<Packet>(1);
  }

  Array<Ordinal> recvCounts(as<Ordinal>(numProcs), as<Ordinal>(1));

  Array<Packet> myGlobalReducts(1);

  Teuchos::reduceAllAndScatter<Ordinal,Packet>(
    *comm, Teuchos::REDUCE_SUM,
    as<Ordinal>(sendBuffer.size()), &sendBuffer[0],
    &recvCounts[0], &myGlobalReducts[0]
    );

  if (std::numeric_limits<Packet>::is_integer) {
    TEST_EQUALITY( myGlobalReducts[0], as<Packet>(numProcs) );
  }
  else {
    const PacketMag local_errorTolSlack = static_cast<PacketMag>(errorTolSlack);
    TEST_FLOATING_EQUALITY( myGlobalReducts[0], as<Packet>(numProcs),
      as<PacketMag>(defaultSmallNumber<PacketMag>() * local_errorTolSlack / numProcs)
      );
  }
  
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, reduceAllAndScatter_2, Ordinal, Packet )
{

  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  const Ordinal procRank = rank(*comm);

  Array<Packet> sendBuffer(as<Ordinal>(numProcs));
  for (Ordinal k = 0; k < numProcs; ++k) {
    sendBuffer[k] = as<Packet>(procRank + k);
  }

  Array<Ordinal> recvCounts(as<Ordinal>(numProcs), as<Ordinal>(1));

  Array<Packet> myGlobalReducts(1);

  Teuchos::reduceAllAndScatter<Ordinal,Packet>(
    *comm, Teuchos::REDUCE_SUM,
    as<Ordinal>(sendBuffer.size()), &sendBuffer[0],
    &recvCounts[0], &myGlobalReducts[0]
    );

  const Packet expectedMyGlobalReduct =
    numProcs * procRank + ((numProcs - 1) * numProcs)/2;

  if (std::numeric_limits<Packet>::is_integer) {
    TEST_EQUALITY( myGlobalReducts[0], expectedMyGlobalReduct );
  }
  else {
    const PacketMag local_errorTolSlack = static_cast<PacketMag>(errorTolSlack);
    TEST_FLOATING_EQUALITY( myGlobalReducts[0], expectedMyGlobalReduct,
      as<PacketMag>(defaultSmallNumber<PacketMag>() * local_errorTolSlack / numProcs)
      );
  }

}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, NonblockingSendReceive_1, Ordinal, Packet )
{

  using Teuchos::as;
  using Teuchos::rcpFromRef;
  using Teuchos::outArg;
  using Teuchos::isend;
  using Teuchos::ireceive;
  using Teuchos::wait;
  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  const Ordinal procRank = rank(*comm);

  if (numProcs == 1) {
    out << "\nCan't test send/receive on one process so successs!";
    return;
  }

  const Packet input_data = as<Packet>(1);
  Packet output_data = as<Packet>(-1);

  RCP<Teuchos::CommRequest> recvRequest;
  RCP<Teuchos::CommRequest> sendRequest;

  if (procRank == 0) {
    // Create copy of data to make sure that peristing relationship is
    // maintained!
    sendRequest = isend<Ordinal, Packet>(
      *comm, Teuchos::rcp(new Packet(input_data)), numProcs-1);
  }
  if (procRank == numProcs-1) {
    // We will need to read output_data after wait(...) below
    recvRequest = ireceive<Ordinal, Packet>(
      *comm, rcpFromRef(output_data), 0);
  }

  if (procRank == 0) {
    wait( *comm, outArg(sendRequest) );
  }
  if (procRank == numProcs-1) {
    wait<Ordinal>( *comm, outArg(recvRequest) );
  }
  
  TEST_EQUALITY_CONST( sendRequest, Teuchos::null );
  TEST_EQUALITY_CONST( recvRequest, Teuchos::null );

  if (procRank == numProcs-1) {
    TEST_EQUALITY( output_data, input_data );
  }

  int globalSuccess_int = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );

}


//
// Instantiations
//


#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME, ORDINAL)\
     typedef std::complex<float> ComplexFloat; \
     TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(TEST_GROUP, TEST_NAME, ORDINAL, ComplexFloat)
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME, ORDINAL)\
     typedef std::complex<double> ComplexDouble; \
     TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(TEST_GROUP, TEST_NAME, ORDINAL, ComplexDouble)
#else
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME, ORDINAL)
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME, ORDINAL)
#endif


#define UNIT_TEST_GROUP_ORDINAL_PACKET( ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, reduceAllAndScatter_1, ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, reduceAllAndScatter_2, ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, NonblockingSendReceive_1, ORDINAL, PACKET )


#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DefaultMpiComm, basic, ORDINAL ) \
  UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, char) \
  UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, int) \
  UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, float) \
  UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, double) \
  UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, reduceAllAndScatter_1, ORDINAL) \
  UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, reduceAllAndScatter_2, ORDINAL) \
  UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, NonblockingSendReceive_1, ORDINAL) \
  UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, reduceAllAndScatter_1, ORDINAL) \
  UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, reduceAllAndScatter_2, ORDINAL) \
  UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, NonblockingSendReceive_1, ORDINAL)


UNIT_TEST_GROUP_ORDINAL(char)
typedef short int ShortInt;
UNIT_TEST_GROUP_ORDINAL(ShortInt)
UNIT_TEST_GROUP_ORDINAL(int)
typedef long int LongInt;
UNIT_TEST_GROUP_ORDINAL(LongInt)

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int LongLongInt;
UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#endif


} // namespace
