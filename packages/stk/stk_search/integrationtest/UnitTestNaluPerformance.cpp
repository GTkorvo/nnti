#include <gtest/gtest.h>
#include <unit_tests/UnitTestUtils.hpp>
#include <unit_tests/MeshUtilsForBoundingVolumes.hpp>

namespace
{

void printGoldResults(const GtkBoxVector &domainBoxes, const std::vector< std::pair<Sphere, Ident> > &spheres)
{
    SearchResults boxIdPairResults;
    for (size_t i=0;i<domainBoxes.size();++i)
    {
        for (size_t j=0;j<spheres.size();++j)
        {
            if ( stk::search::intersects(domainBoxes[i].first, spheres[j].first) )
            {
                boxIdPairResults.push_back(std::make_pair(domainBoxes[i].second, spheres[j].second));
            }
        }
    }
    std::cerr << "Gold: Found " << boxIdPairResults.size() << " interactions.\n";
}

void printGoldResults(const GtkBoxVector &domainBoxes, const GtkBoxVector &spheres)
{
    SearchResults boxIdPairResults;
    for (size_t i=0;i<domainBoxes.size();++i)
    {
        for (size_t j=0;j<spheres.size();++j)
        {
            if ( stk::search::intersects(domainBoxes[i].first, spheres[j].first) )
            {
                boxIdPairResults.push_back(std::make_pair(domainBoxes[i].second, spheres[j].second));
            }
        }
    }
    std::cerr << "Gold: Found " << boxIdPairResults.size() << " interactions.\n";
}

struct Options
{
     std::string mSphereFile;
     std::string mVolumeFile;
     bool mCommunicateRangeBoxes;
     NewSearchMethod mSearchMethod;
     bool mSpheresFirstThenBoxes;
     bool mTestToGetGoldResults;

     void checkForRequiredFile(const std::string &option, const std::string &file)
     {
         ThrowRequireMsg(file != "NO_FILE_SPECIFIED", option << " required for this unit test.");
     }

     void setSphereFile()
     {
         std::string optionString = "-sphere";
         mSphereFile = getOption(optionString, "NO_FILE_SPECIFIED");
         checkForRequiredFile(optionString, mSphereFile);
     }

     void setVolumeFile()
     {
         std::string optionString = "-volume";
         mVolumeFile = getOption(optionString, "NO_FILE_SPECIFIED");
         checkForRequiredFile(optionString, mVolumeFile);
     }

     void setSearchMethod()
     {
         std::string optionString = "-method";
         mSearchMethod = BOOST_RTREE;
         std::string searchString = getOption(optionString, "boost");
         if ( searchString == "octree")
         {
             mSearchMethod = OCTREE;
         }
         else if ( searchString == "gtk" )
         {
             mSearchMethod = GTK;
         }
     }

     void setRangeBoxCommunication()
     {
         std::string optionString = "-rangeBoxComm";
         mCommunicateRangeBoxes = true;
         if ( getOption(optionString, "yes") == "no" )
         {
             mCommunicateRangeBoxes = false;
         }
     }

     void setSphereBoxes()
     {
         std::string optionString = "-sb";
         mSpheresFirstThenBoxes = false;
         if ( getOption(optionString, "no" ) == "yes" )
         {
             mSpheresFirstThenBoxes = true;
         }
     }

     void setIfTestIsGoldTestRun()
     {
         bool mTestToGetGoldResults = false;
         std::string optionString = "-getGold";
         if ( getOption(optionString, "no") == "yes" )
         {
             mTestToGetGoldResults = true;
         }
     }

};

Options getOptionsForTest()
{
    Options local;
    local.setSphereFile();
    local.setVolumeFile();
    local.setSearchMethod();
    local.setRangeBoxCommunication();
    local.setSphereBoxes();
    local.setIfTestIsGoldTestRun();

    return local;
}


TEST(NaluPerformance, BoxSphereIntersections)
{
    Options options = getOptionsForTest();

    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector< std::pair<Sphere, Ident> > spheres;
    fillBoundingVolumesUsingNodesFromFile(comm, options.mSphereFile, spheres);

    GtkBoxVector domainBoxes;
    fillBoxesUsingElementBlocksFromFile(comm, options.mVolumeFile, domainBoxes);

    SearchResults searchResults;

    double startTime = stk::wall_time();

    if ( options.mSpheresFirstThenBoxes )
    {
        stk::search::coarse_search(spheres, domainBoxes, mapSearchMethodToStk(options.mSearchMethod), comm, searchResults, options.mCommunicateRangeBoxes);
    }
    else
    {
        stk::search::coarse_search(domainBoxes, spheres, mapSearchMethodToStk(options.mSearchMethod), comm, searchResults, options.mCommunicateRangeBoxes);
    }

    double elapsedTime = stk::wall_time() - startTime;
    printPeformanceStats(elapsedTime, comm);

    if ( options.mTestToGetGoldResults )
    {
        int numProcs=0;
        MPI_Comm_size(comm, &numProcs);
        if ( numProcs != 1 )
        {
            std::cerr << "Gold results are available only on serial runs.\n";
        }
        else
        {
            printGoldResults(domainBoxes, spheres);
        }
    }
    else
    {
        gatherResultstoProcZero(comm, searchResults);

        int procId=-1;
        MPI_Comm_rank(comm, &procId);
        if ( procId == 0 )
        {
            std::vector< std::pair<int,int> > globalIdMapping(searchResults.size());
            for (size_t i=0; i<searchResults.size(); i++)
            {
                globalIdMapping[i] = std::make_pair(searchResults[i].first.id(), searchResults[i].second.id());
            }
            std::sort(globalIdMapping.begin(), globalIdMapping.end());
            std::vector< std::pair<int,int> >::iterator iter_end = std::unique(globalIdMapping.begin(), globalIdMapping.end());
            globalIdMapping.erase(iter_end, globalIdMapping.end());

            size_t numInteractions = getGoldValueForTest();
            EXPECT_EQ(numInteractions, globalIdMapping.size());
        }
    }
}


TEST(NaluPerformance, BoxBoxIntersections)
{
    Options options = getOptionsForTest();
    MPI_Comm comm = MPI_COMM_WORLD;

    GtkBoxVector spheres;
    fillBoundingVolumesUsingNodesFromFile(comm, options.mSphereFile, spheres);

    GtkBoxVector domainBoxes;
    fillBoxesUsingElementBlocksFromFile(comm, options.mVolumeFile, domainBoxes);

    SearchResults searchResults;

    double startTime = stk::wall_time();

    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    if ( options.mSpheresFirstThenBoxes )
    {
        coarse_search_new(spheres, domainBoxes, options.mSearchMethod, comm, searchResults);
    }
    else
    {
        coarse_search_new(domainBoxes, spheres, options.mSearchMethod, comm, searchResults);
    }

    double elapsedTime = stk::wall_time() - startTime;
    printPeformanceStats(elapsedTime, comm);

    if ( options.mTestToGetGoldResults )
    {
        int numProcs=0;
        MPI_Comm_size(comm, &numProcs);
        if ( numProcs != 1 )
        {
            std::cerr << "Gold results are available only on serial runs.\n";
        }
        else
        {
            printGoldResults(domainBoxes, spheres);
        }
    }
    else
    {
        gatherResultstoProcZero(comm, searchResults);

        if ( procId == 0 )
        {
            std::vector< std::pair<int,int> > globalIdMapping(searchResults.size());
            for (size_t i=0; i<searchResults.size(); i++)
            {
                globalIdMapping[i] = std::make_pair(searchResults[i].first.id(), searchResults[i].second.id());
            }
            std::sort(globalIdMapping.begin(), globalIdMapping.end());
            std::vector< std::pair<int,int> >::iterator iter_end = std::unique(globalIdMapping.begin(), globalIdMapping.end());
            globalIdMapping.erase(iter_end, globalIdMapping.end());

            size_t numInteractions = getGoldValueForTest();
            EXPECT_EQ(numInteractions, globalIdMapping.size());
        }
    }
}

TEST(stkSearch, boxSphereIntersection)
{
    GtkBox box(0,0,0,1,1,1);
    Sphere sphere(Point(2,2,2), 0.5);
    EXPECT_FALSE(stk::search::intersects(box, sphere));
    EXPECT_FALSE(stk::search::intersects(sphere, box));
    Sphere sphere1(Point(1.1, 1.1, 1.1), 0.2);
    EXPECT_TRUE(stk::search::intersects(box, sphere1));
    EXPECT_TRUE(stk::search::intersects(sphere1, box));
    Sphere sphere2(Point(1.1, 1.1, 1.1), 0.17321);
    EXPECT_TRUE(stk::search::intersects(box, sphere2));
    EXPECT_TRUE(stk::search::intersects(sphere2, box));
    Sphere sphere3(Point(0.5, 0.5, 0.5), 1);
    EXPECT_TRUE(stk::search::intersects(box, sphere3));
    EXPECT_TRUE(stk::search::intersects(sphere3, box));
    Sphere sphere4(Point(0.5, 0.5, 0.5), 0.1);
    EXPECT_TRUE(stk::search::intersects(box, sphere4));
    EXPECT_TRUE(stk::search::intersects(sphere4, box));
}

}
