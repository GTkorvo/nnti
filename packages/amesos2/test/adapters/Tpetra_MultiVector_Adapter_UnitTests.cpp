// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_DefaultPlatform.hpp>

#include "Amesos2_TpetraMultiVecAdapter.hpp"
#include "Amesos2_Util.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace {

  using std::cout;
  using std::endl;
  using std::string;

  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::ptrInArg;
  using Teuchos::outArg;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::ScalarTraits;
  using Teuchos::OrdinalTraits;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;

  using Tpetra::global_size_t;
  using Tpetra::MultiVector;
  using Tpetra::DefaultPlatform;
  using Tpetra::Map;

  using Amesos::MultiVecAdapter;
  using Amesos::createMultiVecAdapter;

  using Amesos::Util::is_same;
  using Amesos::Util::get_1d_copy_helper;
  using Amesos::Util::put_1d_data_helper;
  using Amesos::Util::Rooted;
  using Amesos::Util::Distributed;
  using Amesos::Util::Globally_Replicated;


  typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

  bool testMpi = false;

  // Where to look for input files
  string filedir;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption("filedir",&filedir,"Directory of matrix files.");
    clp.addOutputSetupOptions(true);
    clp.setOption("test-mpi", "test-serial", &testMpi,
		  "Test Serial by default (for now) or force MPI test.  In a serial build,"
		  " this option is ignored and a serial comm is always used." );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    RCP<const Comm<int> > ret;
    if( testMpi ){
      ret = DefaultPlatform::getDefaultPlatform().getComm();
    } else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  RCP<FancyOStream> getDefaultOStream()
  {
    return( VerboseObjectBase::getDefaultOStream() );
  }

  /*
   * UNIT TESTS
   */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Initialization, SCALAR, LO, GO )
  {
    /* Test correct initialization of the MultiVecAdapter
     *
     * - All Constructors
     * - Correct initialization of class members
     * - Correct typedefs ( using Amesos::is_same<> )
     */
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    // const size_t numprocs = comm->getSize();
    // const size_t rank     = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    RCP<MV> mv = rcp(new MV(map,11));
    mv->randomize();
    // RCP<FancyOStream> os = getDefaultOStream();
    // mv->describe(*os,Teuchos::VERB_EXTREME);

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);

    // The following should all pass at compile time
    TEST_ASSERT( (is_same<SCALAR,        typename ADAPT::scalar_t>::value) );
    TEST_ASSERT( (is_same<LO,            typename ADAPT::local_ordinal_t>::value) );
    TEST_ASSERT( (is_same<GO,            typename ADAPT::global_ordinal_t>::value) );
    TEST_ASSERT( (is_same<Node,          typename ADAPT::node_t>::value) );
    TEST_ASSERT( (is_same<global_size_t, typename ADAPT::global_size_t>::value) );
    TEST_ASSERT( (is_same<MV,            typename ADAPT::multivec_t>::value) );

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Dimensions, SCALAR, LO, GO )
  {
    // Test that the dimensions reported by the adapter match those as reported
    // by the Tpetra::MultiVector
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    // const size_t numprocs = comm->getSize();
    // const size_t rank     = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    RCP<MV> mv = rcp(new MV(map,11));
    mv->randomize();
    // RCP<FancyOStream> os = getDefaultOStream();
    // mv->describe(*os,Teuchos::VERB_EXTREME);

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);

    TEST_EQUALITY( mv->getLocalLength(),  adapter->getLocalLength()      );
    TEST_EQUALITY( mv->getNumVectors(),   adapter->getLocalNumVectors()  );
    TEST_EQUALITY( mv->getNumVectors(),   adapter->getGlobalNumVectors() );
    TEST_EQUALITY( mv->getGlobalLength(), adapter->getGlobalLength()     );
    TEST_EQUALITY( mv->getStride(),       adapter->getStride()           );
  
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Copy, SCALAR, LO, GO )
  {
    /* Test the get1dCopy() method of MultiVecAdapter.  We can check against a
     * known multivector and also check against what is returned by the
     * Tpetra::MultiVector.
     */
    typedef ScalarTraits<SCALAR> ST;
    typedef typename ScalarTraits<SCALAR>::magnitudeType MAG;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();

    // create a Map
    const size_t numVectors = 7;
    const size_t numLocal = 13;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    RCP<MV> mv = rcp(new MV(map,numVectors));
    mv->randomize();

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);
    Array<SCALAR> original(numVectors*numLocal*numprocs);
    Array<SCALAR> copy(numVectors*numLocal*numprocs);

    get_1d_copy_helper<ADAPT,SCALAR>::do_get(ptrInArg(*adapter), copy(), numLocal*numprocs, Rooted);

    // Only rank=0 process has global copy of the mv data, check against an import
    size_t my_num_elems = OrdinalTraits<size_t>::zero();
    if( rank == 0 ) my_num_elems = numLocal*numprocs;
    Map<LO,GO,Node> root_map(numLocal*numprocs, my_num_elems, 0, comm);
    MV root_mv(rcpFromRef(root_map), numVectors);
    Tpetra::Import<LO,GO,Node> importer(map,rcpFromRef(root_map));
    root_mv.doImport(*mv, importer, Tpetra::REPLACE);

    root_mv.get1dCopy(original(),numLocal*numprocs);

    TEST_EQUALITY( original, copy );

    // Check getting copy of just local data
    original.clear();
    original.resize(numVectors*numLocal);
    copy.clear();
    copy.resize(numVectors*numLocal);
    mv->randomize();
    
    mv->get1dCopy(original(),mv->getLocalLength());
    get_1d_copy_helper<ADAPT,SCALAR>::do_get(ptrInArg(*adapter), copy(), numLocal, Distributed);
    
    // Check that the values remain the same
    TEST_EQUALITY( original, copy );

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Globalize, SCALAR, LO, GO )
  {
    typedef ScalarTraits<SCALAR> ST;
    typedef typename ScalarTraits<SCALAR>::magnitudeType MAG;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();
    // create a Map
    const size_t numVectors = 7;
    const size_t numLocal = 13;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    RCP<MV> mv = rcp(new MV(map,numVectors));
    mv->randomize();

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);
    Array<SCALAR> original(numVectors*numLocal*numprocs);
    Array<SCALAR> copy(numVectors*numLocal*numprocs);

    if( rank == 0 ){
      std::fill(original.begin(), original.end(), 1.9);
    }

    // distribute rank 0's data
    put_1d_data_helper<ADAPT,SCALAR>::do_put(outArg(*adapter), original(),
					     numLocal*numprocs,
					     Rooted);

    // Send rank 0's array to everyone else
    Teuchos::broadcast(*comm, 0, original());

    // Now have everyone get a copy from the multivector adapter
    get_1d_copy_helper<ADAPT,SCALAR>::do_get(ptrInArg(*adapter), copy(),
					     numLocal*numprocs,
					     Globally_Replicated);

    TEST_EQUALITY( original, copy );
  }


  /*
   * Instantiations
   */

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)	\
  typedef std::complex<float> ComplexFloat;		\
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)	\
  typedef std::complex<double> ComplexDouble;			\
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexDouble) 
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR )		\
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Initialization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Dimensions, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Copy, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Globalize, SCALAR, LO, GO ) \
  
#define UNIT_TEST_GROUP_ORDINAL( ORDINAL )		\
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )	\
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double)	\
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT( LO, GO )
  UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )	\
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, float)		\
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, double)	\
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)		\
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO,GO)

  UNIT_TEST_GROUP_ORDINAL(int)

  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongInt )
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongLongInt )
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

} // end anonymous namespace
