#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_as.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

// FINISH: add test for MultiVector with a node containing zero local entries
// FINISH: add tests for local MultiVectors 

namespace Teuchos {
  template <>
    ScalarTraits<int>::magnitudeType
    relErr( const int &s1, const int &s2 )
    {
      typedef ScalarTraits<int> ST;
      return ST::magnitude(s1-s2);
    }

  template <>
    ScalarTraits<char>::magnitudeType
    relErr( const char &s1, const char &s2 )
    {
      typedef ScalarTraits<char> ST;
      return ST::magnitude(s1-s2);
    }
}

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::arrayView;
  using std::copy;
  using std::ostream_iterator;
  using std::string;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

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

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  Node& getDefaultNode()
  {
    return DefaultPlatform::getDefaultPlatform().getNode();
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, basic, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal numLocal = 10;
    const Teuchos_Ordinal numVecs  = 5;
    Map<Ordinal> map(INVALID,numLocal,0,comm);
    MV mvec(node,map,numVecs,true);
    TEST_EQUALITY( mvec.numVectors(), numVecs );
    TEST_EQUALITY( mvec.myLength(), numLocal );
    TEST_EQUALITY( mvec.globalLength(), numImages*numLocal );
    // we zeroed it out in the constructor; all norms should be zero
    Array<Magnitude> norms(numVecs), zeros(numVecs);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Magnitude>::zero());
    mvec.norm2(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    mvec.norm1(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    mvec.normInf(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    // print it
    out << mvec << endl;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstNumVecs, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal numLocal = 10;
    Map<Ordinal> map(INVALID,numLocal,0,comm);
    TEST_THROW(MV mvec(node,map,0), std::invalid_argument);
    TEST_THROW(MV mvec(node,map,-1), std::invalid_argument);
  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstLDA, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    // numlocal > LDA
//REFACTOR//    // ergo, the arrayview doesn't contain enough data to specify the entries
//REFACTOR//    // also, if bounds checking is enabled, check that bad bounds are caught
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 2;
//REFACTOR//    const Teuchos_Ordinal numVecs = 2;
//REFACTOR//    // multivector has two vectors, each proc having two values per vector
//REFACTOR//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//REFACTOR//    // we need 4 scalars to specify values on each proc
//REFACTOR//    Array<Scalar> values(4);
//REFACTOR//#ifdef HAVE_TPETRA_DEBUG
//REFACTOR//    // too small an ArrayView (less than 4 values) is met with an exception, if debugging is on
//REFACTOR//    TEST_THROW(MV mvec(node,map,values(0,3),2,numVecs), std::runtime_error);
//REFACTOR//    // it could also be too small for the given LDA: 
//REFACTOR//    TEST_THROW(MV mvec(node,map,values(),2+1,numVecs), std::runtime_error);
//REFACTOR//    // too small for number of entries in a Vector
//REFACTOR//    TEST_THROW(V   vec(node,map,values(0,1)), std::runtime_error);
//REFACTOR//#endif
//REFACTOR//    // LDA < numLocal throws an exception anytime
//REFACTOR//    TEST_THROW(MV mvec(node,map,values(0,4),1,numVecs), std::runtime_error);
//REFACTOR//  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, LabeledObject, LO, GO, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    Node &node = getDefaultNode();
    // create Map
    Map<LO,GO> map(INVALID,3,0,comm);
    // test labeling
    const string lbl("mvecA");
    MV mvecA(node,map,2);
    string desc1 = mvecA.description();
    if (myImageID==0) out << desc1 << endl;
    mvecA.setObjectLabel(lbl);
    string desc2 = mvecA.description();
    if (myImageID==0) out << desc2 << endl;
    if (myImageID==0) {
      TEST_EQUALITY( mvecA.getObjectLabel(), lbl );
    }
    // test describing at different verbosity levels
    if (myImageID==0) out << "Describing with verbosity VERB_DEFAULT..." << endl;
    mvecA.describe(out);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_NONE..." << endl;
    mvecA.describe(out,VERB_NONE);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_LOW..." << endl;
    mvecA.describe(out,VERB_LOW);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_MEDIUM..." << endl;
    mvecA.describe(out,VERB_MEDIUM);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_HIGH..." << endl;
    mvecA.describe(out,VERB_HIGH);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_EXTREME..." << endl;
    mvecA.describe(out,VERB_EXTREME);
    comm->barrier();
    comm->barrier();
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadMultiply, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    const Scalar S1 = ScalarTraits<Scalar>::one(),
                S0 = ScalarTraits<Scalar>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    {
      // create local Maps
      Map<Ordinal> map3l(3,0,comm,true),
                   map2l(2,0,comm,true);
      MV mvecA(node,map3l,2),
         mvecB(node,map2l,3),
         mvecD(node,map2l,2);
      // failures, 8 combinations:
      // [NTC],[NTC]: A,B don't match
      // [NTC],[NTC]: C doesn't match A,B
      TEST_THROW( mvecD.multiply(NO_TRANS  ,NO_TRANS  ,S1,mvecA,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(NO_TRANS  ,CONJ_TRANS,S1,mvecA,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(CONJ_TRANS,NO_TRANS  ,S1,mvecB,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(CONJ_TRANS,CONJ_TRANS,S1,mvecB,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(NO_TRANS  ,NO_TRANS  ,S1,mvecA,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(NO_TRANS  ,CONJ_TRANS,S1,mvecA,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(CONJ_TRANS,NO_TRANS  ,S1,mvecB,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(CONJ_TRANS,CONJ_TRANS,S1,mvecB,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
    }
    // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
    {
      Map<Ordinal> map3n(INVALID,3,0,comm),
                   map2n(INVALID,2,0,comm),
                   map2l(2,0,comm,true),
                   map3l(3,0,comm,true);
      MV mv3nx2(node,map3n,2),
         mv2nx2(node,map2n,2),
         mv2lx2(node,map2l,2),
         mv2lx3(node,map2l,3),
         mv3lx2(node,map3l,2),
         mv3lx3(node,map3l,3);
      // non-matching input lengths
      TEST_THROW( mv2lx2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv2nx2,S0), std::runtime_error);   // (2 x 3n) x (2n x 2) not compat
      TEST_THROW( mv2lx2.multiply(CONJ_TRANS,NO_TRANS,S1,mv2nx2,mv3nx2,S0), std::runtime_error);   // (2 x 2n) x (3n x 2) not compat
      // non-matching output size
      TEST_THROW( mv3lx3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 3x3
      TEST_THROW( mv3lx2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 3x2
      TEST_THROW( mv2lx3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 2x3
    }
    // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
    {
      Map<Ordinal> map3n(INVALID,3,0,comm),
                   map2n(INVALID,2,0,comm),
                   map2l(2,0,comm,true),
                   map3l(3,0,comm,true);
      MV mv3nx2(node,map3n,2),
         mv2nx2(node,map2n,2),
         mv2x3(node,map2l,3),
         mv3x2(node,map3l,2);
      // non-matching input lengths
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx2,mv2x3,S0), std::runtime_error);   // (3n x 2) x (3 x 2) (trans) not compat
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,S1,mv3nx2,mv3x2,S0), std::runtime_error);   // (3n x 2) x (3 x 2) (nontrans) not compat
      // non-matching output sizes
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx2,mv3x2,S0), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,S1,mv3nx2,mv2x3,S0), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Multiply, Ordinal, Scalar )
  {
    using Teuchos::View;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // const int numImages = comm->getSize();
    Node &node = getDefaultNode();
    // create a Map
    Map<Ordinal> map3n(INVALID,3,0,comm),
                 map2n(INVALID,2,0,comm),
                 lmap3(3,0,comm,true),
                 lmap2(2,0,comm,true);
    const Scalar S1 = ScalarTraits<Scalar>::one(),
                 S0 = ScalarTraits<Scalar>::zero();
    const Mag    M0 = ScalarTraits<Mag>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // deterministic input/output
    {
      MV mv3x2l(node,lmap3,2),
         mv2x3l(node,lmap2,3),
         mv2x2l(node,lmap2,2),
         mv3x3l(node,lmap3,3);
      // fill multivectors with ones
      mv3x2l.putScalar(ScalarTraits<Scalar>::one());
      mv2x3l.putScalar(ScalarTraits<Scalar>::one());
      // fill expected answers Array
      Teuchos::Array<Scalar> check2(4,3); // each entry (of four) is the product [1 1 1]*[1 1 1]' = 3
      Teuchos::Array<Scalar> check3(9,2); // each entry (of nine) is the product [1 1]*[1 1]' = 2
      // test
      Array<Scalar> tmpCopy(9);  // FINISH: after extractConstView1D is finished, use it instead of extractCopy1D
      mv3x3l.multiply(NO_TRANS  ,NO_TRANS  ,S1,mv3x2l,mv2x3l,S0);
      mv3x3l.extractCopy1D(tmpCopy(),3); TEST_COMPARE_FLOATING_ARRAYS(tmpCopy(0,9),check3,M0);
      mv2x2l.multiply(NO_TRANS  ,CONJ_TRANS,S1,mv2x3l,mv2x3l,S0);
      mv2x2l.extractCopy1D(tmpCopy(),2); TEST_COMPARE_FLOATING_ARRAYS(tmpCopy(0,4),check2,M0);
      mv2x2l.multiply(CONJ_TRANS,NO_TRANS  ,S1,mv3x2l,mv3x2l,S0);
      mv2x2l.extractCopy1D(tmpCopy(),2); TEST_COMPARE_FLOATING_ARRAYS(tmpCopy(0,4),check2,M0);
      mv3x3l.multiply(CONJ_TRANS,CONJ_TRANS,S1,mv2x3l,mv3x2l,S0);
      mv3x3l.extractCopy1D(tmpCopy(),3); TEST_COMPARE_FLOATING_ARRAYS(tmpCopy(0,9),check3,M0);
    }
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // random input/output
    {
      Array<Scalar>     tmvCopy1(6), tmvCopy2(6);
      ArrayView<Scalar> sdmView(Teuchos::null);
      MV tmv3x2(node,lmap3,2),
         tmv2x3(node,lmap2,3),
         tmv2x2(node,lmap2,2),
         tmv3x3(node,lmap3,3);
      // fill multivectors with random, get copy of contents
      tmv3x2.random();  tmv3x2.extractCopy1D(tmvCopy1(),3); 
      tmv2x3.random();  tmv2x3.extractCopy1D(tmvCopy2(),2);
      // point SerialDenseMatrices at copies
      SerialDenseMatrix<int,Scalar> sdm3x2(View,tmvCopy1.getRawPtr(),3,3,2);
      SerialDenseMatrix<int,Scalar> sdm2x3(View,tmvCopy2.getRawPtr(),2,2,3);
      // space for answers
      SerialDenseMatrix<int,Scalar> sdm2x2(2,2), sdm3x3(3,3);
      // test: perform local Tpetra::MultiVector multiply and Teuchos::SerialDenseMatrix multiply, then check that answers are equivalent
      {
        Array<Scalar> tmpcopy(9);
        tmv3x3.multiply(NO_TRANS,NO_TRANS,S1,tmv3x2,tmv2x3,S0);
        sdm3x3.multiply(NO_TRANS,NO_TRANS,S1,sdm3x2,sdm2x3,S0);
        tmv3x3.extractCopy1D(tmpcopy(),3); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpcopy(),sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
      {
        Array<Scalar> tmpcopy(4);
        tmv2x2.multiply(NO_TRANS,CONJ_TRANS,S1,tmv2x3,tmv2x3,S0);
        sdm2x2.multiply(NO_TRANS,CONJ_TRANS,S1,sdm2x3,sdm2x3,S0);
        tmv2x2.extractCopy1D(tmpcopy(),2); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpcopy(),sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
      {
        Array<Scalar> tmpcopy(4);
        tmv2x2.multiply(CONJ_TRANS,NO_TRANS,S1,tmv3x2,tmv3x2,S0);
        sdm2x2.multiply(CONJ_TRANS,NO_TRANS,S1,sdm3x2,sdm3x2,S0);
        tmv2x2.extractCopy1D(tmpcopy(),2); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpcopy(),sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
      {
        Array<Scalar> tmpcopy(9);
        tmv3x3.multiply(CONJ_TRANS,CONJ_TRANS,S1,tmv2x3,tmv3x2,S0);
        sdm3x3.multiply(CONJ_TRANS,CONJ_TRANS,S1,sdm2x3,sdm3x2,S0);
        tmv3x3.extractCopy1D(tmpcopy(),3); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpcopy(),sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
    }
    // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
    // {
    //   MV mv3nx2(node,map3n,2),
    //      mv3nx3(node,map3n,3),
    //      // locals
    //      mv2x2(node,lmap2,2),
    //      mv2x3(node,lmap2,3),
    //      mv3x2(node,lmap3,2),
    //      mv3x3(node,lmap3,3);
    //   // fill multivectors with ones
    //   mv3nx3.putScalar(ScalarTraits<Scalar>::one());
    //   mv3nx2.putScalar(ScalarTraits<Scalar>::one());
    //   // fill expected answers Array
    //   ArrayView<const Scalar> tmpView(Teuchos::null); Teuchos_Ordinal dummy;
    //   Teuchos::Array<Scalar> check(9,3*numImages);
    //   // test
    //   mv2x2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0); 
    //   mv2x2.extractConstView1D(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
    //   mv2x3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx3,S0);
    //   mv2x3.extractConstView1D(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
    //   mv3x2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx3,mv3nx2,S0);
    //   mv3x2.extractConstView1D(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
    //   mv3x3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx3,mv3nx3,S0);
    //   mv3x3.extractConstView1D(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
    // }
    // // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
    // {
    //   // FINISH
    // }
  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstAA, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    // constructor takes ArrayView<ArrayView<Scalar> A, NumVectors
//REFACTOR//    // A.size() == NumVectors
//REFACTOR//    // A[i].size() >= MyLength
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 2;
//REFACTOR//    // multivector has two vectors, each proc having two values per vector
//REFACTOR//    Map<Ordinal> map2(INVALID,numLocal  ,indexBase,comm),
//REFACTOR//                 map3(INVALID,numLocal+1,indexBase,comm);
//REFACTOR//    // we need 4 scalars to specify values on each proc
//REFACTOR//    Array<Scalar> values(4);
//REFACTOR//    Array<ArrayView<const Scalar> > arrOfarr(2,ArrayView<const Scalar>(Teuchos::null));
//REFACTOR//    Array<ArrayView<const Scalar> > emptyArr;
//REFACTOR//    arrOfarr[0] = values(0,2);
//REFACTOR//    arrOfarr[1] = values(2,2);
//REFACTOR//    // arrOfarr.size() == 0
//REFACTOR//    TEST_THROW(MV mvec(node,map2,emptyArr(),0), std::runtime_error);
//REFACTOR//#ifdef HAVE_TPETRA_DEBUG
//REFACTOR//    // individual ArrayViews could be too small
//REFACTOR//    TEST_THROW(MV mvec(node,map3,arrOfarr(),2), std::runtime_error);
//REFACTOR//#endif
//REFACTOR//  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    Map<Ordinal> map1(INVALID,1,indexBase,comm),
                 map2(INVALID,2,indexBase,comm);
    {
      MV mv12(node,map1,1),
         mv21(node,map2,1),
         mv22(node,map2,2);
      Array<Scalar> dots(2);
      // incompatible maps
      TEST_THROW(mv12.dot(mv21,dots()),std::runtime_error);
      // incompatible numvecs
      TEST_THROW(mv22.dot(mv21,dots()),std::runtime_error);
      // too small output array
#ifdef TEUCHOS_DEBUG
      TEST_THROW(mv22.dot(mv22,dots(0,1)),std::runtime_error);
#endif
    }
    {
      V v1(node,map1),
        v2(node,map2);
      // incompatible maps
      TEST_THROW(v1.dot(v2),std::runtime_error);
      TEST_THROW(v2.dot(v1),std::runtime_error);
      // wrong size output array through MultiVector interface
      Array<Scalar> dots(2);
#ifdef TEUCHOS_DEBUG
      TEST_THROW(v1.dot(v2,dots()),std::runtime_error);
      TEST_THROW(v2.dot(v1,dots()),std::runtime_error);
#endif
    }
  }


//  ////
//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, OrthoDot, Ordinal, Scalar )
//  {
//    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//    const Scalar S0 = ScalarTraits<Scalar>::zero();
//    const Mag M0 = ScalarTraits<Mag>::zero();
//    // get a comm and node
//    RCP<const Comm<int> > comm = getDefaultComm();
//    Node &node = getDefaultNode();
//    // create a Map
//    const Ordinal indexBase = 0;
//    const Ordinal numLocal = 2;
//    const Teuchos_Ordinal numVectors = 2;
//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//    Array<Scalar> values(5);
//    // values = {0, 1, 0, 1, 0}
//    // values(0,4) = {0, 1, 0, 1} = [0 0]
//    //                            = [1 1]
//    // values(1,4) = {1, 0, 1, 0} = [1 1]
//    //                            = [0 0]
//    // these should be numerical orthogonal even in finite arithmetic
//    values[0] = as<Scalar>(1);
//    values[1] = as<Scalar>(0);
//    values[2] = as<Scalar>(1);
//    values[3] = as<Scalar>(0);
//    values[4] = as<Scalar>(1);
//    MV mvec1(node,map,values(0,4),2,numVectors),
//       mvec2(node,map,values(1,4),2,numVectors);
//    Array<Scalar> dots1(numVectors), dots2(numVectors), zeros(numVectors);
//    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Scalar>::zero());
//    mvec1.dot(mvec2,dots1());
//    mvec2.dot(mvec1,dots2());
//    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,M0);
//    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
//    TEST_EQUALITY_CONST( mvec1(0)->dot(*mvec2(0)), S0);
//  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, ZeroScaleUpdate, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    const Mag M0 = ScalarTraits<Mag>::zero();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 2;
//REFACTOR//    const Teuchos_Ordinal numVectors = 2;
//REFACTOR//    const Teuchos_Ordinal LDA = 2;
//REFACTOR//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//REFACTOR//    Array<Scalar> values(6);
//REFACTOR//    // values = {1, 1, 2, 2, 4, 4}
//REFACTOR//    // values(0,4) = {1, 1, 2, 2} = [1 2]
//REFACTOR//    //                            = [1 2]
//REFACTOR//    // values(2,6) = {2, 2, 4, 4} = [2 4]
//REFACTOR//    //                            = [2 4]
//REFACTOR//    // a multivector A constructed from the first 
//REFACTOR//    // has values .5 of a multivector B constructed from the second
//REFACTOR//    // then 2*A - B = 0
//REFACTOR//    // we test both scale(), both update(), and norm()
//REFACTOR//    values[0] = as<Scalar>(1);
//REFACTOR//    values[1] = as<Scalar>(1);
//REFACTOR//    values[2] = as<Scalar>(2);
//REFACTOR//    values[3] = as<Scalar>(2);
//REFACTOR//    values[4] = as<Scalar>(4);
//REFACTOR//    values[5] = as<Scalar>(4);
//REFACTOR//    MV A(node,map,values(0,4),LDA,numVectors),
//REFACTOR//       B(node,map,values(2,4),LDA,numVectors);
//REFACTOR//    Array<Mag> norms(numVectors), zeros(numVectors);
//REFACTOR//    std::fill(zeros.begin(),zeros.end(),M0);
//REFACTOR//    //
//REFACTOR//    //      [.... ....]
//REFACTOR//    // A == [ones ones] 
//REFACTOR//    //      [.... ....]
//REFACTOR//    // 
//REFACTOR//    //      [.... ....]
//REFACTOR//    // B == [twos twos]
//REFACTOR//    //      [.... ....]
//REFACTOR//    //
//REFACTOR//    //   set A2 = A
//REFACTOR//    //   scale it by 2 in situ
//REFACTOR//    //   check that it equals B: subtraction in situ
//REFACTOR//    {
//REFACTOR//      MV A2(A);
//REFACTOR//      A2.scale(as<Scalar>(2));
//REFACTOR//      A2.update(as<Scalar>(-1),B,as<Scalar>(1));
//REFACTOR//      A2.norm2(norms);
//REFACTOR//      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
//REFACTOR//    }
//REFACTOR//    //   set A2 = A
//REFACTOR//    //   check that it equals B: scale,subtraction in situ
//REFACTOR//    {
//REFACTOR//      MV A2(A);
//REFACTOR//      A2.update(as<Scalar>(-1),B,as<Scalar>(2));
//REFACTOR//      A2.norm2(norms);
//REFACTOR//      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
//REFACTOR//    }
//REFACTOR//    //   set C random
//REFACTOR//    //   set it to zero by combination with A,B
//REFACTOR//    {
//REFACTOR//      MV C(node,map,numVectors);
//REFACTOR//      C.random();
//REFACTOR//      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
//REFACTOR//      C.norm2(norms);
//REFACTOR//      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
//REFACTOR//    }
//REFACTOR//    //   set C random
//REFACTOR//    //   scale it ex-situ
//REFACTOR//    //   check that it equals B: subtraction in situ
//REFACTOR//    {
//REFACTOR//      MV C(node,map,numVectors);
//REFACTOR//      C.scale(as<Scalar>(2),A);
//REFACTOR//      C.update(as<Scalar>(1),B,as<Scalar>(-1));
//REFACTOR//      C.norm2(norms);
//REFACTOR//      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
//REFACTOR//    }
//REFACTOR//  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, ZeroScaleUpdate, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
//REFACTOR//    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    const Mag M0 = ScalarTraits<Mag>::zero();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    Map<Ordinal> map(INVALID,2,0,comm);
//REFACTOR//    Array<Scalar> values(6);
//REFACTOR//    // values = {1, 1, 2, 2}
//REFACTOR//    // values(0,2) = {1, 1} = [1]
//REFACTOR//    //                      = [1]
//REFACTOR//    // values(2,2) = {2, 2} = [2]
//REFACTOR//    //                      = [2]
//REFACTOR//    // a vector A constructed from the first 
//REFACTOR//    // has values .5 of a vector B constructed from the second
//REFACTOR//    // thus 2*A - B = 0
//REFACTOR//    // we test both scale(), both update(), and norm()
//REFACTOR//    values[0] = as<Scalar>(1);
//REFACTOR//    values[1] = as<Scalar>(1);
//REFACTOR//    values[2] = as<Scalar>(2);
//REFACTOR//    values[3] = as<Scalar>(2);
//REFACTOR//    V A(node,map,values(0,2)),
//REFACTOR//      B(node,map,values(2,2));
//REFACTOR//    Mag norm;
//REFACTOR//    Array<Mag> norms(1);
//REFACTOR//    //
//REFACTOR//    //      [....]
//REFACTOR//    // A == [ones]
//REFACTOR//    //      [....]
//REFACTOR//    // 
//REFACTOR//    //      [....]
//REFACTOR//    // B == [twos]
//REFACTOR//    //      [....]
//REFACTOR//    //
//REFACTOR//    //   set A2 = A
//REFACTOR//    //   scale it by 2 in situ
//REFACTOR//    //   check that it equals B: subtraction in situ
//REFACTOR//    {
//REFACTOR//      V A2(A);
//REFACTOR//      A2.scale(as<Scalar>(2));
//REFACTOR//      A2.update(as<Scalar>(-1),B,as<Scalar>(1));
//REFACTOR//      norm = A2.norm2(); A2.norm2(norms());
//REFACTOR//      TEST_EQUALITY(norm,M0);
//REFACTOR//      TEST_EQUALITY(norm,norms[0]);
//REFACTOR//    }
//REFACTOR//    //   set A2 = A
//REFACTOR//    //   check that it equals B: scale,subtraction in situ
//REFACTOR//    {
//REFACTOR//      V A2(A);
//REFACTOR//      A2.update(as<Scalar>(-1),B,as<Scalar>(2));
//REFACTOR//      norm = A2.norm2(); A2.norm2(norms());
//REFACTOR//      TEST_EQUALITY(norm,M0);
//REFACTOR//      TEST_EQUALITY(norm,norms[0]);
//REFACTOR//    }
//REFACTOR//    //   set C random
//REFACTOR//    //   set it to zero by combination with A,B
//REFACTOR//    {
//REFACTOR//      V C(map);
//REFACTOR//      C.random();
//REFACTOR//      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
//REFACTOR//      norm = C.norm2(); C.norm2(norms());
//REFACTOR//      TEST_EQUALITY(norm,M0);
//REFACTOR//      TEST_EQUALITY(norm,norms[0]);
//REFACTOR//    }
//REFACTOR//    //   set C random
//REFACTOR//    //   scale it ex-situ
//REFACTOR//    //   check that it equals B: subtraction in situ
//REFACTOR//    {
//REFACTOR//      V C(map);
//REFACTOR//      C.random();
//REFACTOR//      C.scale(as<Scalar>(2),A);
//REFACTOR//      C.update(as<Scalar>(1),B,as<Scalar>(-1));
//REFACTOR//      norm = C.norm2(); C.norm2(norms());
//REFACTOR//      TEST_EQUALITY(norm,M0);
//REFACTOR//      TEST_EQUALITY(norm,norms[0]);
//REFACTOR//    }
//REFACTOR//  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    // create random MV
    MV morig(node,map,numVectors);
    morig.random();
    // copy it
    MV mcopy1(morig), mcopy2(morig);
    // verify that all three have identical values
    Array<Magnitude> norig(numVectors), ncopy1(numVectors), ncopy2(numVectors);
    morig.normInf(norig);
    mcopy1.normInf(ncopy1);
    mcopy2.normInf(ncopy2);
    TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy1,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy2,M0);
    TEST_COMPARE_FLOATING_ARRAYS(ncopy1,ncopy2,M0);
    // modify all three
    morig.putScalar(as<Scalar>(0));
    mcopy1.putScalar(as<Scalar>(1));
    mcopy2.putScalar(as<Scalar>(2));
    // compute norms
    morig.normInf(norig);
    mcopy1.normInf(ncopy1);
    mcopy2.normInf(ncopy2);
    // check them
    bool local_success = true;
    for (Teuchos_Ordinal i=0; i<numVectors; ++i) {
      TEST_ARRAY_ELE_EQUALITY( norig,  i, as<Scalar>(0) );
      TEST_ARRAY_ELE_EQUALITY( ncopy1, i, as<Scalar>(1) );
      TEST_ARRAY_ELE_EQUALITY( ncopy2, i, as<Scalar>(2) );
    }
    success &= local_success;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    Map<Ordinal> map(INVALID,2,0,comm);
    // create random MV
    V morig(node,map);
    morig.random();
    // copy it
    V mcopy1(morig), mcopy2(morig);
    // verify that all three have identical values
    Magnitude norig, ncopy1, ncopy2;
    norig = morig.normInf();
    ncopy1 = mcopy1.normInf();
    ncopy2 = mcopy2.normInf();
    TEST_EQUALITY(norig,ncopy1);
    TEST_EQUALITY(norig,ncopy2);
    TEST_EQUALITY(ncopy1,ncopy2);
    // modify all three
    morig.putScalar(as<Scalar>(0));
    mcopy1.putScalar(as<Scalar>(1));
    mcopy2.putScalar(as<Scalar>(2));
    // compute norms
    norig = morig.normInf();
    ncopy1 = mcopy1.normInf();
    ncopy2 = mcopy2.normInf();
    // check them
    TEST_EQUALITY(norig, as<Scalar>(0));
    TEST_EQUALITY(ncopy1,as<Scalar>(1));
    TEST_EQUALITY(ncopy2,as<Scalar>(2));
  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, SingleVecNormalize, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    // this documents a usage case in Anasazi::SVQBOrthoManager, which was failing
//REFACTOR//    // error turned out to be a neglected return in both implementations of update(), 
//REFACTOR//    // after passing the buck to scale() in the case of alpha==0 or beta==0 or gamma=0
//REFACTOR//    if (ScalarTraits<Scalar>::isOrdinal) return;
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    const Magnitude M1  = ScalarTraits<Magnitude>::one();
//REFACTOR//    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 10;
//REFACTOR//    const Teuchos_Ordinal numVectors = 6;
//REFACTOR//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//REFACTOR//    // create random MV
//REFACTOR//    MV mv(node,map,numVectors);
//REFACTOR//    mv.random();
//REFACTOR//    // compute the norms
//REFACTOR//    Array<Magnitude> norms(numVectors);
//REFACTOR//    mv.norm2(norms());
//REFACTOR//    for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
//REFACTOR//      // get a view of column j, normalize it using update()
//REFACTOR//      Array<Teuchos_Ordinal> ind(1,j);
//REFACTOR//      RCP<MV> mvj = mv.subView(ind());
//REFACTOR//      switch (j){
//REFACTOR//      case 0:
//REFACTOR//        mvj->scale( M1/norms[j] );
//REFACTOR//        break;
//REFACTOR//      case 1:
//REFACTOR//        mvj->update( M1/norms[j], *mvj, M0 );
//REFACTOR//        break;
//REFACTOR//      case 2:
//REFACTOR//        mvj->update( M0        , *mvj, M1/norms[j] );
//REFACTOR//        break;
//REFACTOR//      case 3:
//REFACTOR//        mvj->update( M0        , *mvj, M1/norms[j], *mvj, M0 );
//REFACTOR//        break;
//REFACTOR//      case 4:
//REFACTOR//        mvj->update( M1/norms[j], *mvj, M0        , *mvj, M0 );
//REFACTOR//        break;
//REFACTOR//      case 5:
//REFACTOR//        mvj->update( M0        , *mvj, M0        , *mvj, M1/norms[j] );
//REFACTOR//        break;
//REFACTOR//      }
//REFACTOR//    }
//REFACTOR//    mv.norm2(norms()); // should be all one now
//REFACTOR//    Array<Magnitude> ones(numVectors,M1);
//REFACTOR//    TEST_COMPARE_FLOATING_ARRAYS(norms,ones,ScalarTraits<Magnitude>::eps()*as<Magnitude>(10.));
//REFACTOR//  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDot, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    const int numImages = comm->getSize();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 2;
//REFACTOR//    const Teuchos_Ordinal numVectors = 3;
//REFACTOR//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//REFACTOR//    Array<Scalar> values(6);
//REFACTOR//    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
//REFACTOR//    //                               [0 1 2]
//REFACTOR//    // dot(values,values) = [0*0+0*0 1+1*1 + 2*2+2*2] = [0 2 8]
//REFACTOR//    // summed over all procs, this is [0 2*nprocs 8*nprocs]
//REFACTOR//    values[0] = as<Scalar>(0);
//REFACTOR//    values[1] = as<Scalar>(0);
//REFACTOR//    values[2] = as<Scalar>(1);
//REFACTOR//    values[3] = as<Scalar>(1);
//REFACTOR//    values[4] = as<Scalar>(2);
//REFACTOR//    values[5] = as<Scalar>(2);
//REFACTOR//    MV mvec1(node,map,values(),2,numVectors),
//REFACTOR//       mvec2(node,map,values(),2,numVectors);
//REFACTOR//    Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
//REFACTOR//    answer[0] = as<Scalar>(0);
//REFACTOR//    answer[1] = as<Scalar>(2*numImages);
//REFACTOR//    answer[2] = as<Scalar>(8*numImages);
//REFACTOR//    // do the dots
//REFACTOR//    mvec1.dot(mvec2,dots1());
//REFACTOR//    mvec2.dot(mvec1,dots2());
//REFACTOR//    // check the answers
//REFACTOR//    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,M0);
//REFACTOR//    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,M0);
//REFACTOR//  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDotNonTrivLDA, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    // same as CountDot, but the A,LDA has a non-trivial LDA (i.e., LDA != myLen)
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    const int numImages = comm->getSize();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 2;
//REFACTOR//    const Teuchos_Ordinal numVectors = 3;
//REFACTOR//    const Teuchos_Ordinal LDA = 3;
//REFACTOR//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//REFACTOR//    Array<Scalar> values(9);
//REFACTOR//    // A = {0, 0, -1, 1, 1, -1, 2, 2, -1} = [0   1  2]
//REFACTOR//    //                                      [0   1  2]
//REFACTOR//    //                                      [-1 -1 -1]
//REFACTOR//    // processed as a 2 x 3 with LDA==3, the result it
//REFACTOR//    //            values =       [0 1 2]
//REFACTOR//    //                           [0 1 2]
//REFACTOR//    // dot(values,values) = [0*0+0*0 1+1*1 + 2*2+2*2] = [0 2 8]
//REFACTOR//    // summed over all procs, this is [0 2*nprocs 8*nprocs]
//REFACTOR//    values[0] = as<Scalar>(0);
//REFACTOR//    values[1] = as<Scalar>(0);
//REFACTOR//    values[2] = as<Scalar>(-1);
//REFACTOR//    values[3] = as<Scalar>(1);
//REFACTOR//    values[4] = as<Scalar>(1);
//REFACTOR//    values[5] = as<Scalar>(-1);
//REFACTOR//    values[6] = as<Scalar>(2);
//REFACTOR//    values[7] = as<Scalar>(2);
//REFACTOR//    values[8] = as<Scalar>(-1);
//REFACTOR//    MV mvec1(node,map,values(),LDA,numVectors),
//REFACTOR//       mvec2(node,map,values(),LDA,numVectors);
//REFACTOR//    Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
//REFACTOR//    answer[0] = as<Scalar>(0);
//REFACTOR//    answer[1] = as<Scalar>(2*numImages);
//REFACTOR//    answer[2] = as<Scalar>(8*numImages);
//REFACTOR//    // do the dots
//REFACTOR//    mvec1.dot(mvec2,dots1());
//REFACTOR//    mvec2.dot(mvec1,dots2());
//REFACTOR//    // check the answers
//REFACTOR//    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,M0);
//REFACTOR//    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,M0);
//REFACTOR//  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNorm1, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    const MT M0 = ScalarTraits<MT>::zero();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    const int numImages = comm->getSize();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 2;
//REFACTOR//    const Ordinal numVectors = 3;
//REFACTOR//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//REFACTOR//    Array<Scalar> values(6);
//REFACTOR//    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
//REFACTOR//    //                               [0 1 2]
//REFACTOR//    // norm1(values) = [0 2 4]
//REFACTOR//    // over all procs, this is [0 2*nprocs 4*nprocs]
//REFACTOR//    values[0] = as<Scalar>(0);
//REFACTOR//    values[1] = as<Scalar>(0);
//REFACTOR//    values[2] = as<Scalar>(1);
//REFACTOR//    values[3] = as<Scalar>(1);
//REFACTOR//    values[4] = as<Scalar>(2);
//REFACTOR//    values[5] = as<Scalar>(2);
//REFACTOR//    MV mvec(node,map,values(),2,numVectors);
//REFACTOR//    Array<MT> norms(numVectors), answer(numVectors);
//REFACTOR//    answer[0] = as<MT>(0);
//REFACTOR//    answer[1] = as<MT>(2*numImages);
//REFACTOR//    answer[2] = as<MT>(4*numImages);
//REFACTOR//    // do the dots
//REFACTOR//    mvec.norm1(norms());
//REFACTOR//    // check the answers
//REFACTOR//    TEST_COMPARE_FLOATING_ARRAYS(norms,answer,M0);
//REFACTOR//  }


//REFACTOR//  ////
//REFACTOR//  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNormInf, Ordinal, Scalar )
//REFACTOR//  {
//REFACTOR//    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
//REFACTOR//    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
//REFACTOR//    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
//REFACTOR//    const MT M0 = ScalarTraits<MT>::zero();
//REFACTOR//    // get a comm and node
//REFACTOR//    RCP<const Comm<int> > comm = getDefaultComm();
//REFACTOR//    Node &node = getDefaultNode();
//REFACTOR//    // create a Map
//REFACTOR//    const Ordinal indexBase = 0;
//REFACTOR//    const Ordinal numLocal = 2;
//REFACTOR//    const Ordinal numVectors = 3;
//REFACTOR//    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
//REFACTOR//    Array<Scalar> values(6);
//REFACTOR//    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
//REFACTOR//    //                               [0 1 2]
//REFACTOR//    // normInf(values) = [0 1 2]
//REFACTOR//    // over all procs, this is [0 1 2]
//REFACTOR//    values[0] = as<Scalar>(0);
//REFACTOR//    values[1] = as<Scalar>(0);
//REFACTOR//    values[2] = as<Scalar>(1);
//REFACTOR//    values[3] = as<Scalar>(1);
//REFACTOR//    values[4] = as<Scalar>(2);
//REFACTOR//    values[5] = as<Scalar>(2);
//REFACTOR//    MV mvec(node,map,values(),2,numVectors);
//REFACTOR//    Array<MT> norms(numVectors), answer(numVectors);
//REFACTOR//    answer[0] = as<MT>(0);
//REFACTOR//    answer[1] = as<MT>(1);
//REFACTOR//    answer[2] = as<MT>(2);
//REFACTOR//    // do the dots
//REFACTOR//    mvec.normInf(norms());
//REFACTOR//    // check the answers
//REFACTOR//    TEST_COMPARE_FLOATING_ARRAYS(norms,answer,M0);
//REFACTOR//  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Norm2, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Ordinal numVectors = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    MV mvec(node,map,2,numVectors);
    // randomize the multivector
    mvec.random();
    // take norms; they should not be zero
    Array<MT> normsRand(numVectors), normsZero(numVectors);
    mvec.norm2(normsRand());
    // zero the vector
    mvec.putScalar(ScalarTraits<Scalar>::zero());
    // take norms; they should be zero
    mvec.norm2(normsZero());
    // check the answers
    bool local_success = true;
    for (Teuchos_Ordinal i=0; i<numVectors; ++i) {
      TEST_ARRAY_ELE_INEQUALITY(normsRand,i,M0);
      TEST_ARRAY_ELE_EQUALITY(normsZero,i,M0);
    }
    success &= local_success;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, MinMaxMean, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 2;
    const Teuchos_Ordinal LDA = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(4);
    // on proc i of n:
    // values = {i, i+1, n-i, n-i-1} = [ i   n-i ]
    //                                 [i+1 n-i-1]
    // min values are [0,0], but from different procs
    // max values are [numImages,numImages], again from different procs
    // mean values are 
    //  1  n-1
    // --  sum i+i+1 == [(n-1)*n + n]/2n == (n-1)/2 
    // 2n  i=0
    //               == n*n/2n == n/2
    values[0] = as<Scalar>(myImageID);
    values[1] = as<Scalar>(myImageID+1);
    values[2] = as<Scalar>(numImages-myImageID);
    values[3] = as<Scalar>(numImages-myImageID-1);
    MV mvec(node,map,values(),LDA,numVectors);
    // FINISH this test
  }


#define PRINT_TYPE_AND_VALUE(val) { out << std::setw(30) << #val << std::setw(30) << Teuchos::typeName(val) << ": " << val << endl; }



  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadCombinations, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Scalar rnd = ScalarTraits<Scalar>::random();
    // two maps: one has two entires per node, the other disagrees on node 0
    Map<Ordinal> map1(INVALID,2,indexBase,comm),
                 map2(INVALID,(myImageID == 0 ? 1 : 2),indexBase,comm);
    // multivectors from different maps are incompatible for all ops
    // multivectors from the same map are compatible only if they have the same number of
    //    columns
    MV m1n1(node,map1,1), m1n2(node,map1,2), m2n2(node,map2,2), m1n2_2(node,map1,2);
    Array<Scalar> dots(1);
    Array<Mag>    norms(1);
    // FINISH: test multiply (both), reciprocalMultiply
    TEST_THROW(m1n2.dot(m1n1,dots()), std::runtime_error); // dot
    TEST_THROW(m1n2.dot(m2n2,dots()), std::runtime_error);
    TEST_THROW(m1n2.abs(m1n1), std::runtime_error);       // abs
    TEST_THROW(m1n2.abs(m2n2), std::runtime_error);
    TEST_THROW(m1n2.abs(m1n1), std::runtime_error);       // abs
    TEST_THROW(m1n2.abs(m2n2), std::runtime_error);
    TEST_THROW(m1n2.scale(rnd,m1n1), std::runtime_error); // abs
    TEST_THROW(m1n2.scale(rnd,m2n2), std::runtime_error);
    TEST_THROW(m1n2.update(rnd,m1n1,rnd), std::runtime_error); // update(alpha,A,beta)
    TEST_THROW(m1n2.update(rnd,m2n2,rnd), std::runtime_error);
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m1n2_2,rnd), std::runtime_error); // update(alpha,A,beta,B,gamma) // A incompat
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m1n2_2,rnd), std::runtime_error); // incompt is length            // A incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m2n2  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m2n2  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m2n2  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m2n2  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n2_2,rnd), std::runtime_error); // incompt is numVecs           // A incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n2_2,rnd), std::runtime_error);                                 // A incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m1n1  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m1n1  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n1  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n1  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n1.normWeighted(m1n2,norms()), std::runtime_error);        // normWeighted
    TEST_THROW(m1n2.normWeighted(m2n2,norms()), std::runtime_error);
    TEST_THROW(m1n2.reciprocal(m1n1), std::runtime_error);                  // reciprocal
    TEST_THROW(m1n2.reciprocal(m2n2), std::runtime_error);
  }


  /* TODO 
     Many constructors left to test
     MultiVector (const Map< Ordinal > &map, const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &arrayOfArrays, Ordinal numVectors)

     MultiVector<Ordinal,Scalar> subCopy(const Teuchos::Range1D &colRng) const;
     MultiVector<Ordinal,Scalar> subCopy(const Teuchos::ArrayView<Teuchos_Index> &cols) const;
     MultiVector<Ordinal,Scalar> subView(const Teuchos::Range1D &colRng);
     MultiVector<Ordinal,Scalar> subView(const Teuchos::ArrayView<Teuchos_Index> &cols);
     const MultiVector<Ordinal,Scalar> subViewConst(const Teuchos::Range1D &colRng) const;
     const MultiVector<Ordinal,Scalar> subViewConst(const Teuchos::ArrayView<Teuchos_Index> &cols) const;

     Mod routines left to test
     void replaceGlobalValue (Ordinal globalRow, Ordinal vectorIndex, const Scalar &value)
     void sumIntoGlobalValue (Ordinal globalRow, Ordinal vectorIndex, const Scalar &value)
     void replaceMyValue (Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)
     void sumIntoMyValue (Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)

     Arithmetic methods left to test:
     void multiply (Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
     void multiply (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
     void reciprocalMultiply (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
  */

// 
// INSTANTIATIONS
//

#ifdef HAVE_TEUCHOS_BLASFLOAT
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT(ORDINAL)\
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
#else
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT(ORDINAL)
#endif

#ifdef HAVE_TEUCHOS_COMPLEX
#  ifdef HAVE_TEUCHOS_BLASFLOAT
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)\
       typedef std::complex<float> ComplexFloat; \
       UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexFloat)
#  else
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
#  endif
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)\
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexDouble)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, basic             , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstNumVecs   , ORDINAL, SCALAR ) \
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstLDA       , ORDINAL, SCALAR ) */\
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstAA        , ORDINAL, SCALAR ) */\
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CopyConst         , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(      Vector, CopyConst         , ORDINAL, SCALAR ) \
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, OrthoDot          , ORDINAL, SCALAR ) */\
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountDot          , ORDINAL, SCALAR ) */\
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountDotNonTrivLDA, ORDINAL, SCALAR ) */\
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadDot            , ORDINAL, SCALAR ) \
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountNorm1        , ORDINAL, SCALAR ) */\
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountNormInf      , ORDINAL, SCALAR ) */\
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, Norm2             , ORDINAL, SCALAR ) \
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, ZeroScaleUpdate   , ORDINAL, SCALAR ) */\
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(      Vector, ZeroScaleUpdate   , ORDINAL, SCALAR ) */\
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadCombinations   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadMultiply       , ORDINAL, SCALAR ) \
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, SingleVecNormalize, ORDINAL, SCALAR ) */\
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, Multiply          , ORDINAL, SCALAR ) */\
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, LabeledObject     , ORDINAL, ORDINAL, SCALAR ) 


#ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
     UNIT_TEST_GROUP_ORDINAL(int)

#else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_FLOAT(ORDINAL) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)  \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)
     typedef long int LongInt;   UNIT_TEST_GROUP_ORDINAL(LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt; UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#    endif

#endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
