#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Array.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Teuchos_as.hpp"

#include "Tpetra_SerialPlatform.hpp"
#ifdef HAVE_TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#endif

#include <Kokkos_SerialNode.hpp>

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

  using Teuchos::OrdinalTraits;
  using Teuchos::RCP;
  using Teuchos::Comm;
  using Tpetra::Platform;
  using Teuchos::rcp;
  using Tpetra::SerialPlatform;
#ifdef HAVE_TPETRA_MPI
  using Tpetra::MpiPlatform;
#endif
  using Kokkos::SerialNode;

  template <class PLAT>
  RCP<PLAT> getPlatform() {
    TEST_FOR_EXCEPTION(true,std::logic_error,"Platform type not defined.");
  }

  SerialNode snode;
  template <>
  RCP<SerialPlatform<float,int,int,SerialNode> > getPlatform() {
    return rcp(new SerialPlatform<float,int,int,SerialNode>(snode));
  }

#ifdef TPETRA_HAVE_MPI
  template <>
  RCP<MpiPlatform<float,int,int,SerialNode> > getPlatform() {
    return rcp(new MpiPlatform<float,int,int,SerialNode>(snode));
  }
#endif

#define PRINT_VECTOR(v) \
   { \
     out << #v << ": "; \
     copy(v.begin(), v.end(), ostream_iterator<Ordinal>(out," ")); \
     out << endl; \
   }

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Platform, basic, PlatformType )
  {
    out << "Testing " << Teuchos::TypeNameTraits<PlatformType>::name() << std::endl;
    typedef typename PlatformType::ScalarType          S;
    typedef typename PlatformType::LocalOrdinalType   LO;
    typedef typename PlatformType::GlobalOrdinalType  GO;
    typedef typename PlatformType::NodeType            N;
    // create a platform  
    RCP<Platform<S,LO,GO,N> > platform = getPlatform<PlatformType>();
    platform->setObjectLabel("not the default label");
    // get the comm for this platform
    RCP<Comm<int> > comm = platform->getComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    TEST_EQUALITY( myImageID < numImages, true );
    // clone the platform and get the new comm, test that it is different
    RCP<const Platform<S,LO,GO,N> > platform2 = platform->clone();
    RCP<Comm<int> > comm2 = platform2->getComm();
    TEST_EQUALITY_CONST( comm == comm2, false );
    N &node  = platform->getNode();
    N &node2 = platform2->getNode();
    TEST_EQUALITY( &node, &node2 );
    TEST_EQUALITY_CONST( platform->getObjectLabel() == platform2->getObjectLabel(), false );
  }


  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

typedef SerialPlatform<float,int,int,SerialNode> SP;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Platform, basic, SP)
#ifdef HAVE_TPETRA_MPI
typedef MpiPlatform<float,int,int,SerialNode> MP;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Platform, basic, MP )
#endif

}
