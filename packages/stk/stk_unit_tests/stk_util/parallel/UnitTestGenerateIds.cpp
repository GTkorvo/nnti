// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>
#include <vector>
#include <mpi.h>
#include <algorithm>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <fstream>

#include <stdint.h>

namespace {

typedef int INTMPI;

class MpiInfo
{
public:
    MpiInfo(MPI_Comm comm) : mProcId(-1), mNumProcs(-1), mComm(comm)
    {
        MPI_Comm_rank(mComm, &mProcId);
        MPI_Comm_size(mComm, &mNumProcs);
    }
    ~MpiInfo() {}

    INTMPI getProcId() const { return mProcId; }
    INTMPI getNumProcs() const { return mNumProcs; }
    MPI_Comm getMpiComm() const { return mComm; }

private:
    INTMPI mProcId;
    INTMPI mNumProcs;
    MPI_Comm mComm;
};

INTMPI whichProcOwnsId(const uint64_t maxId, const uint64_t id, INTMPI numProcs);
void distributeIds(std::vector<uint64_t> &myIds, const MpiInfo& mpiInfo);
uint64_t getNumIdsPerProc(const uint64_t maxId, const INTMPI numProcs);
void getIdUsageAcrossAllProcs(std::vector< std::vector<uint64_t> > &idsToComm, std::vector<uint64_t> &idsInUseAcrossAllProcsInMyRange, const MpiInfo& mpiInfo);

void retrieveIds(const INTMPI root, uint64_t id, MPI_Comm comm, uint64_t numIdsToGetPerProc, std::vector<int>& areIdsBeingUsed);
void receiveIdAndGetIds(const int root, const uint64_t maxId, const std::vector<uint64_t> &idsInUse, MPI_Comm comm);

bool sendIdToCheck(const INTMPI root, uint64_t id, MPI_Comm comm);
void receiveIdAndCheck(const int root, const std::vector<uint64_t> &idsInUse, MPI_Comm comm);
void checkUniqueIds(const std::vector<uint64_t> &myIds, const std::vector<uint64_t> &uniqueIds, const MpiInfo& mpiInfo);
void writeIdsToFile(const std::string &filename, const INTMPI myProcId, const std::vector<uint64_t>& myIds, const std::vector<uint64_t> &uniqueIds);


void generate_ids(const uint64_t maxId, const std::vector<uint64_t> &idsInUse, std::vector<uint64_t> &uniqueIds, const MpiInfo& mpiInfo)
{
    std::vector< std::vector<uint64_t> > idsToComm(mpiInfo.getNumProcs());

    for (size_t i=0;i<idsInUse.size();i++)
    {
        INTMPI procOwner = whichProcOwnsId(maxId, idsInUse[i], mpiInfo.getNumProcs());
        ThrowRequireMsg(static_cast<int>(procOwner)<mpiInfo.getNumProcs(), "Id generation error. Please contact sierra-help. procOwner = " << procOwner << ", maxId = " << maxId << ", and id = " << idsInUse[i] << std::endl);
        idsToComm[procOwner].push_back(idsInUse[i]);
    }

    std::vector<uint64_t> idsInUseAcrossAllProcsInMyRange;
    getIdUsageAcrossAllProcs(idsToComm, idsInUseAcrossAllProcsInMyRange, mpiInfo);

    std::vector<uint64_t> uniqueIdsFound;

    uint64_t myIndexStart = getNumIdsPerProc(maxId, mpiInfo.getNumProcs())*mpiInfo.getProcId()+1;
    uint64_t myIndexEnd = getNumIdsPerProc(maxId, mpiInfo.getNumProcs())*(mpiInfo.getProcId()+1)+1;

    size_t numIdsNeeded = uniqueIds.size();
    size_t numIdsAvailableThisProc =  getNumIdsPerProc(maxId, mpiInfo.getNumProcs()) - idsInUseAcrossAllProcsInMyRange.size();
    ThrowRequireMsg(numIdsNeeded <= numIdsAvailableThisProc, "Not enough unique ids available (Id generation error). Plrease contact sierra-help for support. Number of ids needed: "
            << numIdsNeeded << " and num available ids: " << numIdsAvailableThisProc);

    if ( !uniqueIds.empty() )
    {
        for (uint64_t i=myIndexStart;i<myIndexEnd;i++)
        {
            if ( !std::binary_search(idsInUseAcrossAllProcsInMyRange.begin(), idsInUseAcrossAllProcsInMyRange.end(), i) )
            {
                uniqueIdsFound.push_back(i);
                if ( uniqueIdsFound.size() == uniqueIds.size() )
                {
                    break;
                }
            }
        }
    }

    ThrowRequireMsg(uniqueIdsFound.size() == uniqueIds.size(), "Id generation error. Could not obtain needed ids. Please contact sierra-help for support.");
    std::copy(uniqueIdsFound.begin(), uniqueIdsFound.end(), uniqueIds.begin());
}

void getAvailableIds(const std::vector<uint64_t> &myIds, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startId, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    std::vector<uint64_t> receivedInfo(mpiInfo.getNumProcs(),0);
    MPI_Allgather(&numIdsNeeded, 1, MPI_UINT64_T, &receivedInfo[0], 1, MPI_UINT64_T, mpiInfo.getMpiComm());

    std::vector<uint64_t> sortedIds(myIds.begin(), myIds.end());
    std::sort(sortedIds.begin(), sortedIds.end());

    for (INTMPI procIndex=0;procIndex<mpiInfo.getNumProcs();procIndex++)
    {
        if ( receivedInfo[procIndex] != 0 )
        {
            if ( procIndex == mpiInfo.getProcId() )
            {
                while ( startId < maxId && idsObtained.size() < numIdsNeeded )
                {
                    int requestNumIds = numIdsNeeded - idsObtained.size();
                    std::vector<int> areIdsBeingUsed(requestNumIds,0);
                    if ( !std::binary_search(sortedIds.begin(), sortedIds.end(), startId) )
                    {
                        retrieveIds(procIndex, startId, mpiInfo.getMpiComm(), requestNumIds, areIdsBeingUsed);
                        for (size_t i=0;i<areIdsBeingUsed.size();i++)
                        {
                            if ( areIdsBeingUsed[i] == 0 )
                            {
                                idsObtained.push_back(startId+i);
                            }
                        }
                        startId += areIdsBeingUsed.size();
                    }
                    else
                    {
                        startId++;
                    }
                }
                ThrowRequireMsg(idsObtained.size()==numIdsNeeded, "Id generation error. Ran out of ids. Please contact sierra-help for support.");
                sendIdToCheck(procIndex, 0, mpiInfo.getMpiComm()); // a zero terminates communication
            }
            else
            {
                receiveIdAndGetIds(procIndex, maxId, sortedIds, mpiInfo.getMpiComm());
            }
            // updated starting id across all procs
            MPI_Bcast(&startId, 1, MPI_UINT64_T, procIndex, mpiInfo.getMpiComm());
        }
    }
}

TEST(GeneratedIds, whichProcOwnsId)
{
    uint64_t maxId=10;
    INTMPI numProcs=3;

    uint64_t ids[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    INTMPI procOwner[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3};
    for (int i=0;i<10;i++)
    {
        EXPECT_EQ(procOwner[i], whichProcOwnsId(maxId, ids[i], numProcs));
    }
}

TEST(GeneratedIds, distributeIds)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    int numIdsPerProc = 10;
    std::vector<uint64_t> myIds(numIdsPerProc);
    distributeIds(myIds, mpiInfo);

    uint64_t offsetForProc = numIdsPerProc*mpiInfo.getProcId();
    for (size_t i=0;i<myIds.size();i++)
    {
        uint64_t oneBasedId = i+1;
        EXPECT_EQ(oneBasedId+offsetForProc, myIds[i]);
    }
}

TEST(GeneratedIds, getNumIdsPerProc)
{
    INTMPI numProcs = 10;
    uint64_t findIdsLessThan = 1000;

    uint64_t numIdsOneOneProc = getNumIdsPerProc(findIdsLessThan, numProcs);
    EXPECT_EQ(findIdsLessThan/numProcs, numIdsOneOneProc);
}

TEST(GeneratedIds, getIdUsageAcrossAllProcs)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() != 1 )
    {
        std::vector< std::vector<uint64_t> > idsToComm(mpiInfo.getNumProcs());
        INTMPI procToSendId = 0;
        uint64_t madeUpId = static_cast<uint64_t>(mpiInfo.getProcId());
        idsToComm[procToSendId].push_back(madeUpId);
        std::vector<uint64_t> idsInUseAcrossAllProcsInMyRange;
        getIdUsageAcrossAllProcs(idsToComm, idsInUseAcrossAllProcsInMyRange, mpiInfo);
        if ( mpiInfo.getProcId() == 0 )
        {
            for (uint64_t id=0;id<idsInUseAcrossAllProcsInMyRange.size();id++)
            {
                EXPECT_EQ(id, idsInUseAcrossAllProcsInMyRange[id]);
            }
        }
        else
        {
            EXPECT_EQ(0u, idsInUseAcrossAllProcsInMyRange.size());
        }
    }
}

TEST(GeneratedIds, findUniqueIdAcrossProcs)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    uint64_t totalIdsInUse = 362880;
    uint64_t numIdsThisProc = totalIdsInUse/mpiInfo.getNumProcs();;

    std::vector<uint64_t> myIds(numIdsThisProc,0);
    distributeIds(myIds, mpiInfo);

    uint64_t startId = 1;
    uint64_t numIdsNeeded = 100;
    std::vector<uint64_t> idsObtained;
    uint64_t maxId = 10000000;

    double startTime = stk::wall_time();

    getAvailableIds(myIds, numIdsNeeded, idsObtained, startId, maxId, mpiInfo);

    double endTime = stk::wall_time();

    if ( mpiInfo.getProcId() == 0 )
    {
        std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
    }

    checkUniqueIds(myIds, idsObtained, mpiInfo);
    writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, idsObtained);

    for (size_t i=0;i<idsObtained.size();i++)
    {
        uint64_t goldId = totalIdsInUse + numIdsNeeded*mpiInfo.getProcId()+1 + i;
        EXPECT_EQ(goldId, idsObtained[i]);
    }
}

TEST(GeneratedIds, findUniqueIdAcrossProcsVaryingNumIdsInUse)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    uint64_t totalIdsInUse = 10000*mpiInfo.getNumProcs();
    uint64_t numIdsThisProc = totalIdsInUse/mpiInfo.getNumProcs();;

    std::vector<uint64_t> myIds(numIdsThisProc,0);
    distributeIds(myIds, mpiInfo);

    uint64_t startId = 1;
    uint64_t numIdsNeeded = 100;
    std::vector<uint64_t> idsObtained;
    uint64_t maxId = 10000000;

    double startTime = stk::wall_time();

    getAvailableIds(myIds, numIdsNeeded, idsObtained, startId, maxId, mpiInfo);

    double endTime = stk::wall_time();

    if ( mpiInfo.getProcId() == 0 )
    {
        std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
    }

    checkUniqueIds(myIds, idsObtained, mpiInfo);
    writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, idsObtained);

    for (size_t i=0;i<idsObtained.size();i++)
    {
        uint64_t goldId = totalIdsInUse + numIdsNeeded*mpiInfo.getProcId()+1 + i;
        EXPECT_EQ(goldId, idsObtained[i]);
    }
}

TEST(GeneratedIds, numIdsInUseVariesWithNumberProcessorsNumIdsNeededVariesWithNumberProcessors)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() < 10 )
    {
        uint64_t maxId = 10000000;
        uint64_t numIdsThisProc = 100000;

        std::vector<uint64_t> myIds(numIdsThisProc,0);
        distributeIds(myIds, mpiInfo);

        uint64_t numIdsNeeded = 10000;
        std::vector<uint64_t> uniqueIds(numIdsNeeded,0);

        MPI_Barrier(mpiInfo.getMpiComm());

        double startTime = stk::wall_time();

        generate_ids(maxId, myIds, uniqueIds, mpiInfo);

        double endTime = stk::wall_time();

        if ( mpiInfo.getProcId() == 0 )
        {
            std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
        }

        EXPECT_EQ(numIdsNeeded, uniqueIds.size());
        checkUniqueIds(myIds, uniqueIds, mpiInfo);

        writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, uniqueIds);
    }
}

// 2^56 = 72, 057, 594, 037, 927, 936 --> 72 quadrillion, 57 trillion, 594 billion, 37 million, 927 thousand, 936 ids!

TEST(GeneratedIds, numIdsInUseConstantNumIdsNeededConstant)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() <= 10 )
    {
        uint64_t totalIdsInUse = 362880;
        uint64_t numIdsThisProc = totalIdsInUse/mpiInfo.getNumProcs();;

        std::vector<uint64_t> myIds(numIdsThisProc,0);
        distributeIds(myIds, mpiInfo);

        uint64_t maxId = 10000000;

        uint64_t numIdsNeeded = 10000;
        if ( mpiInfo.getProcId() != 0 )
        {
            numIdsNeeded = 0;
        }

        std::vector<uint64_t> uniqueIds(numIdsNeeded,0);

        MPI_Barrier(mpiInfo.getMpiComm());

        double startTime = stk::wall_time();

        generate_ids(maxId, myIds, uniqueIds, mpiInfo);

        double endTime = stk::wall_time();

        if ( mpiInfo.getProcId() == 0 )
        {
            std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
        }

        EXPECT_EQ(numIdsNeeded, uniqueIds.size());
        checkUniqueIds(myIds, uniqueIds, mpiInfo);

        writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, uniqueIds);
    }
}

TEST(GeneratedIds, numIdsInUseConstantNumIdsNeededVariesWithNumberProcessors)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() <= 10 )
    {
        uint64_t totalIdsInUse = 362880;
        uint64_t numIdsThisProc = totalIdsInUse/mpiInfo.getNumProcs();;

        std::vector<uint64_t> myIds(numIdsThisProc,0);
        distributeIds(myIds, mpiInfo);

        uint64_t maxId = 10000000;

        uint64_t numIdsNeeded = 10000*mpiInfo.getNumProcs();
        if ( mpiInfo.getProcId() != 0 )
        {
            numIdsNeeded = 0;
        }

        std::vector<uint64_t> uniqueIds(numIdsNeeded,0);

        MPI_Barrier(mpiInfo.getMpiComm());

        double startTime = stk::wall_time();

        generate_ids(maxId, myIds, uniqueIds, mpiInfo);

        double endTime = stk::wall_time();

        if ( mpiInfo.getProcId() == 0 )
        {
            std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
        }

        EXPECT_EQ(numIdsNeeded, uniqueIds.size());
        checkUniqueIds(myIds, uniqueIds, mpiInfo);

        writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, uniqueIds);
    }
}

////////////////////////////////////////////////////////////////////

void distributeIds(std::vector<uint64_t> &myIds, const MpiInfo &mpiInfo)
{
    uint64_t offset = myIds.size();

    for (uint64_t i=0;i<myIds.size();i++)
    {
        myIds[i] = offset*mpiInfo.getProcId() + i + 1;
    }
}

////////////////////////////////////////////////////////////////////

void getIdUsageAcrossAllProcs(std::vector< std::vector<uint64_t> > &idsToComm, std::vector<uint64_t> &idsInUseAcrossAllProcsInMyRange, const MpiInfo &mpiInfo)
{
    for (INTMPI i=0;i<mpiInfo.getNumProcs();i++)
    {
        if ( i == mpiInfo.getProcId() )
        {
            for (INTMPI j=0;j<mpiInfo.getNumProcs();j++)
            {
                if ( j != mpiInfo.getProcId() )
                {
                    uint64_t numItemsToComm = idsToComm[j].size();
                    MPI_Send(&numItemsToComm, 1, MPI_UINT64_T, j, mpiInfo.getNumProcs()*i+j, mpiInfo.getMpiComm());
                    if ( numItemsToComm > 0 )
                    {
                        MPI_Send(&idsToComm[j][0], numItemsToComm, MPI_UINT64_T, j, mpiInfo.getNumProcs()*i+j, mpiInfo.getMpiComm());
                    }
                }
            }
        }
        else
        {
            uint64_t numItemsToReceive=0;
            MPI_Status status;
            MPI_Recv(&numItemsToReceive, 1, MPI_UINT64_T, i, mpiInfo.getNumProcs()*i+mpiInfo.getProcId(), mpiInfo.getMpiComm(), &status);
            if ( numItemsToReceive > 0 )
            {
                std::vector<uint64_t> idsFromOtherProc(numItemsToReceive,0);
                MPI_Request request;
                MPI_Irecv(&idsFromOtherProc[0], numItemsToReceive, MPI_UINT64_T, i, mpiInfo.getNumProcs()*i+mpiInfo.getProcId(), mpiInfo.getMpiComm(), &request);
                MPI_Status status;
                MPI_Wait(&request, &status);
                idsInUseAcrossAllProcsInMyRange.insert(idsInUseAcrossAllProcsInMyRange.end(), idsFromOtherProc.begin(), idsFromOtherProc.end());
            }
        }
    }

    idsInUseAcrossAllProcsInMyRange.insert(idsInUseAcrossAllProcsInMyRange.end(), idsToComm[mpiInfo.getProcId()].begin(), idsToComm[mpiInfo.getProcId()].end());
    std::sort(idsInUseAcrossAllProcsInMyRange.begin(), idsInUseAcrossAllProcsInMyRange.end());
}

////////////////////////////////////////////////////////////////////

void checkUniqueIds(const std::vector<uint64_t> &myIds, const std::vector<uint64_t> &uniqueIds, const MpiInfo &mpiInfo)
{
    std::vector<uint64_t> sortedIds(uniqueIds.begin(), uniqueIds.end());
    std::sort(sortedIds.begin(), sortedIds.end());

    for (INTMPI i=0;i<mpiInfo.getNumProcs();i++)
    {
        if (mpiInfo.getNumProcs() == i)
        {
            for (size_t j=0;j<uniqueIds.size();j++)
            {
                ThrowRequireMsg(uniqueIds[j]>0, "Id generation error. Please contact sierra-help for support.");
                EXPECT_FALSE(std::binary_search(myIds.begin(), myIds.end(), uniqueIds[j]));
                EXPECT_EQ(true, sendIdToCheck(i, uniqueIds[j], mpiInfo.getMpiComm()));
            }
            sendIdToCheck(i, 0, mpiInfo.getMpiComm());
        }
        else
        {
            receiveIdAndCheck(i, sortedIds, mpiInfo.getMpiComm());
        }
    }
}

////////////////////////////////////////////////////////////////////

void writeIdsToFile(const std::string &filename, const INTMPI myProcId, const std::vector<uint64_t>& myIds, const std::vector<uint64_t> &uniqueIds)
{
    std::ostringstream os;
    os << "ids_" << myProcId << ".m";
    std::ofstream out(os.str().c_str());
//    for (size_t i=0;i<myIds.size();i++)
//    {
//        out << "myids(" << i+1 << ") = " << myIds[i] << ";" << std::endl;
//    }
    for (size_t i=0;i<uniqueIds.size();i++)
    {
        out << "unique_ids(" << i+1 << ") = " << uniqueIds[i] << ";" << std::endl;
    }
    out.close();
}

////////////////////////////////////////////////////////////////////

uint64_t getNumIdsPerProc(const uint64_t maxId, const INTMPI numProcs)
{
    return maxId/numProcs;
}

////////////////////////////////////////////////////////////////////

INTMPI whichProcOwnsId(const uint64_t maxId, const uint64_t id, INTMPI numProcs)
{
    ThrowRequireMsg(id>0, "Invalid Id. Contact sierra-help for support.");
    uint64_t numIdsPerProc = getNumIdsPerProc(maxId, numProcs);
    INTMPI procOwner = (id-1)/numIdsPerProc;
    return procOwner;
}

////////////////////////////////////////////////////////////////////

bool sendIdToCheck(const INTMPI root, uint64_t id, MPI_Comm comm)
{
    MPI_Bcast(&id, 1, MPI_UINT64_T, root, comm);
    bool goodId = true;
    if ( id != 0 )
    {
        uint64_t good = 0;
        uint64_t received = 0;
        MPI_Reduce(&good, &received, 1, MPI_UINT64_T, MPI_SUM, root, comm);

        if ( received > 0 )
        {
            goodId = false;
        }
    }
    return goodId;
}

////////////////////////////////////////////////////////////////////

void receiveIdAndCheck(const int root, const std::vector<uint64_t> &idsInUse, MPI_Comm comm)
{
    uint64_t id=0;

    while( true )
    {
        MPI_Bcast(&id, 1, MPI_UINT64_T, root, comm);
        if ( id == 0) break;

        bool found = std::binary_search(idsInUse.begin(), idsInUse.end(), id);
        uint64_t result = 0;
        if ( found )
        {
            result = 1;
        }
        MPI_Reduce(&result, &result, 1, MPI_UINT64_T, MPI_SUM, root, comm);
    }
}

////////////////////////////////////////////////////////////////////


void receiveIdAndGetIds(const int root, const uint64_t maxId, const std::vector<uint64_t> &idsInUse, MPI_Comm comm)
{
    uint64_t id=0;

    while( true )
    {
        MPI_Bcast(&id, 1, MPI_UINT64_T, root, comm);
        if ( id == 0) break;
        uint64_t numIdsToGet=0;
        MPI_Bcast(&numIdsToGet, 1, MPI_UINT64_T, root, comm);

        std::vector<int> areIdsBeingused(numIdsToGet,0);
        for (size_t i=0;i<areIdsBeingused.size();i++)
        {
            if ( std::binary_search(idsInUse.begin(), idsInUse.end(), id+i) )
            {
                areIdsBeingused[i] = 1;
            }
        }
        int *rbuff = 0;
        MPI_Reduce(&areIdsBeingused[0], rbuff, numIdsToGet, MPI_INT, MPI_SUM, root, comm);
    }
}

////////////////////////////////////////////////////////////////////

void retrieveIds(const INTMPI root, uint64_t id, MPI_Comm comm, uint64_t numIdsToGetPerProc, std::vector<int>& areIdsBeingUsed)
{
    MPI_Bcast(&id, 1, MPI_UINT64_T, root, comm);
    if ( id != 0 )
    {
        MPI_Bcast(&numIdsToGetPerProc, 1, MPI_UINT64_T, root, comm);
        std::vector<uint64_t> zeroids(numIdsToGetPerProc,0);
        MPI_Reduce(&zeroids[0], &areIdsBeingUsed[0], numIdsToGetPerProc, MPI_INT, MPI_SUM, root, comm);
    }
}
}
