/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_Tuple.hpp>
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "example.hpp"

#define NUM_GLOBAL_ELEMENTS 100

using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
typedef int LO;
typedef long GO;
typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
typedef Tpetra::Map<int, long, Node> Map;
typedef Tpetra::Directory<LO, GO, Node> Directory;
    

namespace {
  TEUCHOS_UNIT_TEST( OneToOne, AlreadyOneToOne)
  {
    //Creates a map that is already one-to-one and tests to make sure 
    //it is not changed.
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();
  
    RCP<const Map> map = Tpetra::createUniformContigMapWithNode<LO, GO>(NUM_GLOBAL_ELEMENTS, comm, node);

    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);

    TEST_ASSERT(map->isSameAs(*new_map));
  }

  TEUCHOS_UNIT_TEST(OneToOne, LargeOverlap)
  {
    //Creates a map with large overlaps
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();
    
    int unit = (NUM_GLOBAL_ELEMENTS/(numProc+1));
    int num_loc_elems=0;
    Array<GO> elementList;

    if(myRank<numProc-1)
    {
      elementList = Array<GO>(2*unit);
      num_loc_elems=2*unit;
    }
    else //I'm the last proc and need to take the leftover elements.
    {
      elementList = Array<GO>((2*unit)+(NUM_GLOBAL_ELEMENTS%numProc));
      num_loc_elems=(2*unit)+(NUM_GLOBAL_ELEMENTS%numProc);
    }

    for(int i=(myRank*unit),j=0;i<(myRank*unit)+(2*unit);++i,++j)
    {
      elementList[j]=i;
    }
    //Don't forget to assign leftovers to the last proc!
    if(myRank==numProc-1 && NUM_GLOBAL_ELEMENTS%numProc != 0)
    {
      for(int i=(myRank*unit)+2*unit, j=2*unit;i<NUM_GLOBAL_ELEMENTS;++i,++j)
      {
        elementList[j]=i;
      }
    }

    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    //std::cout<<map->description();
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    //std::cout<<new_map->description();
    //Now we need to test if we lost anything.
    
    Array<int> nodeIDlist(num_loc_elems); //This is unimportant. It's only used in the following function.
    Tpetra::LookupStatus stat=new_map->getRemoteIndexList(elementList (), nodeIDlist ()); //this style from a tpetra test.
    //If we pass this we didn't lose any IDs.
    TEST_EQUALITY_CONST(stat, Tpetra::AllIDsPresent);
    TEST_EQUALITY(new_map->getGlobalNumElements(),NUM_GLOBAL_ELEMENTS);

    //Now we need to make sure they're in the right place. Keep in mind Tpetra
    //Directory gives precidence to the higher numbered proc.

    ArrayView<const GO> my_owned = new_map->getNodeElementList();

    for(int i=(myRank*unit),j=0;i<(myRank*unit)+unit; ++i,++j)
    {
      TEST_EQUALITY(my_owned[j], i);
    }
    //And the last proc should have all that it started with.
    if(myRank==numProc-1)
    {
      for(int i=0;i<num_loc_elems;++i)
      {
        TEST_EQUALITY(elementList[i],my_owned[i]);
      }
    }
  }

  TEUCHOS_UNIT_TEST(OneToOne, AllOnOneProc)
  {
    //Will create a non-contig map with all of the elements on a single
    //processor. Nothing should change.
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();

    Array<GO> elementList;
    if(myRank==0)
    {
      elementList = Array<GO>(NUM_GLOBAL_ELEMENTS);
      for(int i=0;i<NUM_GLOBAL_ELEMENTS;++i)
      {
        elementList[i]=i;
      }
    }
    else
    {
      elementList = Array<GO>(0);
    }
    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    TEST_ASSERT(map->isSameAs(*new_map));
  }
  
  TEUCHOS_UNIT_TEST(OneToOne, NoIDs)
  {
    //An empty map.
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();

    Array<GO> elementList (0);

    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    TEST_ASSERT(map->isSameAs(*new_map));


  }

  TEUCHOS_UNIT_TEST(OneToOne, AllOwnEvery)
  {
    //Every processor starts by owning all of them.
    //After one-to-one, only the last processor should own all of them.
    
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();
    
    Array<GO> elementList (NUM_GLOBAL_ELEMENTS);

    for(int i=0;i<NUM_GLOBAL_ELEMENTS;++i)
    {
      elementList[i]=i;
    }
    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    
    if(myRank<numProc-1)//I shouldn't have any elements.
    {
      TEST_EQUALITY(new_map->getNodeNumElements(),0);
    }
    else//I should have all of them.
    {
      TEST_EQUALITY(new_map->getNodeNumElements(),NUM_GLOBAL_ELEMENTS);
    }
  }
}
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
typedef int LO;
typedef long GO;
typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
typedef Tpetra::Map<int, long, Node> Map;
typedef Tpetra::Directory<LO, GO, Node> Directory;
    

namespace {
  TEUCHOS_UNIT_TEST( OneToOne, AlreadyOneToOne)
  {
    //Creates a map that is already one-to-one and tests to make sure 
    //it is not changed.
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();
  
    RCP<const Map> map = Tpetra::createUniformContigMapWithNode<LO, GO>(NUM_GLOBAL_ELEMENTS, comm, node);

    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);

    TEST_ASSERT(map->isSameAs(*new_map));
  }

  TEUCHOS_UNIT_TEST(OneToOne, LargeOverlap)
  {
    //Creates a map with large overlaps
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();
    
    int unit = (NUM_GLOBAL_ELEMENTS/(numProc+1));
    int num_loc_elems=0;
    Array<GO> elementList;

    if(myRank<numProc-1)
    {
      elementList = Array<GO>(2*unit);
      num_loc_elems=2*unit;
    }
    else //I'm the last proc and need to take the leftover elements.
    {
      elementList = Array<GO>((2*unit)+(NUM_GLOBAL_ELEMENTS%numProc));
      num_loc_elems=(2*unit)+(NUM_GLOBAL_ELEMENTS%numProc);
    }

    for(int i=(myRank*unit),j=0;i<(myRank*unit)+(2*unit);++i,++j)
    {
      elementList[j]=i;
    }
    //Don't forget to assign leftovers to the last proc!
    if(myRank==numProc-1 && NUM_GLOBAL_ELEMENTS%numProc != 0)
    {
      for(int i=(myRank*unit)+2*unit, j=2*unit;i<NUM_GLOBAL_ELEMENTS;++i,++j)
      {
        elementList[j]=i;
      }
    }

    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    //std::cout<<map->description();
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    //std::cout<<new_map->description();
    //Now we need to test if we lost anything.
    
    Array<int> nodeIDlist(num_loc_elems); //This is unimportant. It's only used in the following function.
    Tpetra::LookupStatus stat=new_map->getRemoteIndexList(elementList (), nodeIDlist ()); //this style from a tpetra test.
    //If we pass this we didn't lose any IDs.
    TEST_EQUALITY_CONST(stat, Tpetra::AllIDsPresent);
    TEST_EQUALITY(new_map->getGlobalNumElements(),NUM_GLOBAL_ELEMENTS);

    //Now we need to make sure they're in the right place. Keep in mind Tpetra
    //Directory gives precidence to the higher numbered proc.

    ArrayView<const GO> my_owned = new_map->getNodeElementList();

    for(int i=(myRank*unit),j=0;i<(myRank*unit)+unit; ++i,++j)
    {
      TEST_EQUALITY(my_owned[j], i);
    }
    //And the last proc should have all that it started with.
    if(myRank==numProc-1)
    {
      for(int i=0;i<num_loc_elems;++i)
      {
        TEST_EQUALITY(elementList[i],my_owned[i]);
      }
    }
  }

  TEUCHOS_UNIT_TEST(OneToOne, AllOnOneProc)
  {
    //Will create a non-contig map with all of the elements on a single
    //processor. Nothing should change.
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();

    Array<GO> elementList;
    if(myRank==0)
    {
      elementList = Array<GO>(NUM_GLOBAL_ELEMENTS);
      for(int i=0;i<NUM_GLOBAL_ELEMENTS;++i)
      {
        elementList[i]=i;
      }
    }
    else
    {
      elementList = Array<GO>(0);
    }
    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    TEST_ASSERT(map->isSameAs(*new_map));
  }
  
  TEUCHOS_UNIT_TEST(OneToOne, NoIDs)
  {
    //An empty map.
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();

    Array<GO> elementList (0);

    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    TEST_ASSERT(map->isSameAs(*new_map));


  }

  TEUCHOS_UNIT_TEST(OneToOne, AllOwnEvery)
  {
    //Every processor starts by owning all of them.
    //After one-to-one, only the last processor should own all of them.
    
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    RCP<Node> node = platform.getNode();
    const int myRank = comm->getRank();
    const int numProc = comm->getSize();
    
    Array<GO> elementList (NUM_GLOBAL_ELEMENTS);

    for(int i=0;i<NUM_GLOBAL_ELEMENTS;++i)
    {
      elementList[i]=i;
    }
    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO>(elementList,comm,node);
    RCP<const Map> new_map = createOneToOne<LO,GO,Node>(map);
    
    if(myRank<numProc-1)//I shouldn't have any elements.
    {
      TEST_EQUALITY(new_map->getNodeNumElements(),0);
    }
    else//I should have all of them.
    {
      TEST_EQUALITY(new_map->getNodeNumElements(),NUM_GLOBAL_ELEMENTS);
    }
  }
}
