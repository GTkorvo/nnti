
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_Directory.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"


//==============================================================================
// Epetra_Directory constructor for a Epetra_BlockMap object
Epetra_Directory::Epetra_Directory(Epetra_BlockMap* Map)
  : Map_(Map),
    DirectoryMap_(0),
    ProcList_(0),
    LocalIndexList_(0),
    AllMinGIDs_(0)
{

  // Test for simple cases

  // Uniprocessor and local map cases (Nothing to set up)

  if (!(Map->DistributedGlobal())) return;

  // Linear Map case

  else if (Map->LinearMap()) {

    // Build a list of the Minimum global ids for all processors on each processor.
    // Since the map is linear, we know that all GIDs are contiguous on each processor
    // and can be found using the MinGIDs.

    int NumProc = Map->Comm().NumProc();
    AllMinGIDs_ = new int[NumProc+1];
    int MinMyGID = Map->MinMyGID();
    Map->Comm().GatherAll(&MinMyGID, AllMinGIDs_, 1);
    AllMinGIDs_[NumProc] = 1 + Map->MaxAllGID(); // Set max cap
  }

  // General case.  Need to build a directory via calls to communication functions
  else {
    
    int flag = Generate();
    assert(flag==0);
  }
}
//==============================================================================
// Epetra_Directory copy constructor
Epetra_Directory::Epetra_Directory(const Epetra_Directory & Directory)
  : Map_(Directory.Map_),
    DirectoryMap_(0),
    ProcList_(0),
    LocalIndexList_(0),
    AllMinGIDs_(0)
{
  int i;

  if (Directory.DirectoryMap_!=0) DirectoryMap_ = new Epetra_Map(Directory.DirectoryMap());

  int Dir_NumMyElements = DirectoryMap_->NumMyElements();

  if (Directory.ProcList_!=0) {
    ProcList_ = new int[Dir_NumMyElements];
    for (i=0; i<Dir_NumMyElements; i++) ProcList_[i] = Directory.ProcList_[i];
  }
  if (Directory.LocalIndexList_!=0) {
    LocalIndexList_ = new int[Dir_NumMyElements];
    for (int i=0; i<Dir_NumMyElements; i++) LocalIndexList_[i] = Directory.LocalIndexList_[i];
    }
  if (Directory.AllMinGIDs_!=0) {
    int NumProc = Map_->Comm().NumProc();
    AllMinGIDs_ = new int[NumProc+1];
    for (int i=0; i<NumProc+1; i++) AllMinGIDs_[i] = Directory.AllMinGIDs_[i];
    }

}
//==============================================================================
// Epetra_Directory destructor 
Epetra_Directory::~Epetra_Directory()
{
  if( DirectoryMap_ != 0 ) delete DirectoryMap_;
  if( ProcList_ != 0 ) delete [] ProcList_;
  if( LocalIndexList_ != 0 ) delete [] LocalIndexList_;
  if( AllMinGIDs_ != 0 ) delete [] AllMinGIDs_;

  DirectoryMap_ = 0;
  ProcList_ = 0 ;
  LocalIndexList_ = 0;
  AllMinGIDs_ = 0;
}

//==============================================================================
// Generate: Generates Directory Tables
int Epetra_Directory::Generate()
{
  int i;

  int MinAllGID = Map_->MinAllGID();
  int MaxAllGID = Map_->MaxAllGID();
  // DirectoryMap will have a range of elements from the minimum to the maximum
  // GID of the user map, and an IndexBase of MinAllGID from the user map
  int Dir_NumGlobalElements = MaxAllGID - MinAllGID + 1;

  // Create a uniform linear map to contain the directory
  DirectoryMap_ = new Epetra_Map( Dir_NumGlobalElements, MinAllGID, Map_->Comm() );

  int Dir_NumMyElements = DirectoryMap_->NumMyElements(); // Get NumMyElements



  // Allocate Processor list and Local Index List.  Initialize to -1s.

  if (Dir_NumMyElements>0) {
    ProcList_ = new int[ Dir_NumMyElements ];
    LocalIndexList_ = new int[ Dir_NumMyElements ];
    
    // Initialize values to -1 in case the user global element list does
    // fill all IDs from MinAllGID to MaxAllGID (e.g., allows global indices to be 
    // all even integers.
    for (i=0; i<Dir_NumMyElements; i++) {
      ProcList_[i] = -1;
      LocalIndexList_[i] = -1;
    }
  }

  
  // Get list of processors owning the directory entries for the Map GIDs

  int NumProc = Map_->Comm_->NumProc();
  int MyPID = Map_->Comm_->MyPID();

  int Map_NumMyElements = Map_->NumMyElements();
  int * send_procs = 0;
  if (Map_NumMyElements>0) send_procs = new int[Map_NumMyElements];
  int * Map_MyGlobalElements = Map_->MyGlobalElements();

  assert(DirectoryMap_->RemoteIDList(Map_NumMyElements, Map_MyGlobalElements, 
				     send_procs, 0)==0); 



  bool det_flag = true;

  int num_recvs;
    
  Epetra_Distributor * Distor = Map_->Comm().CreateDistributor();

  EPETRA_CHK_ERR(Distor->CreateFromSends( Map_NumMyElements, send_procs, det_flag, num_recvs ));

  if (Map_NumMyElements>0) delete [] send_procs;

  int * export_elements = 0;
  int * import_elements = 0;

  if (Map_NumMyElements>0) {
    export_elements = new int[ 3 * Map_NumMyElements ];

    for( i = 0; i < Map_NumMyElements; i++ )
      {
	export_elements[3*i] = Map_MyGlobalElements[i];
	export_elements[3*i+1] = MyPID;
	export_elements[3*i+2] = i;
      }
  }

  if (num_recvs>0) import_elements = new int[ 3 * num_recvs ];

  EPETRA_CHK_ERR(Distor->Do(reinterpret_cast<char *> (export_elements), 
		  3 * sizeof( int ),
		  reinterpret_cast<char *> (import_elements) ));

  
  int curr_LID;
  for( i = 0; i < num_recvs; i++ )
  {
    curr_LID = DirectoryMap_->LID(import_elements[3*i]); // Convert incoming GID to Directory LID
    assert(curr_LID !=-1); // Internal error
    ProcList_[ curr_LID ] = import_elements[3*i+1];
    LocalIndexList_[ curr_LID ] = import_elements[3*i+2];
  }

  if (import_elements!=0) delete [] import_elements;
  if (export_elements!=0) delete [] export_elements;
  
  delete Distor;
  return(0);
}
//==============================================================================
// GetDirectoryEntries: Get non-local GID references ( procID and localID )
// 			Space should already be allocated for Procs and
//     			LocalEntries.
int Epetra_Directory::GetDirectoryEntries( const int NumEntries,
					  const int * GlobalEntries, int * Procs,
					  int * LocalEntries, int * EntrySizes ) {



  int ierr = 0;
  int j;
  int i;
  int MyPID = Map_->Comm().MyPID();
  int NumProc = Map_->Comm().NumProc();
  int n_over_p = Map_->NumGlobalElements() / NumProc;
  int remainder = Map_->NumGlobalElements() % NumProc;


  // Test for simple cases

  // Uniprocessor and local map cases

  if (!Map_->DistributedGlobal()) {
    for (i=0; i<NumEntries; i++) {
      if (LocalEntries!=0) LocalEntries[i] = Map_->LID(GlobalEntries[i]);

      // If GID is not valid, return -1 for both the proc and local entry info
      if (LocalEntries!=0)  if (LocalEntries[i]==-1) Procs[i] = -1; 
      else Procs[i] = MyPID;

    }
    if (EntrySizes!=0) {
      if (Map_->ConstantElementSize()) {
	int ElementSize = Map_->MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
      }
      else {
	int * ElementSizeList = Map_->ElementSizeList();
	for (i=0; i<NumEntries; i++) if (LocalEntries[i]>-1) EntrySizes[i] = ElementSizeList[LocalEntries[i]];
      }
    }
    EPETRA_CHK_ERR(ierr);
    return(0);
  }

  // Linear Map case
  if (Map_->LinearMap()) {
    
    int MinAllGID = Map_->MinAllGID(); // Get Min of all GID
    int MaxAllGID = Map_->MaxAllGID(); // Get Max of all GID
    for (i=0; i<NumEntries; i++) {
      int LID = -1; // Assume not found
      int Proc = -1;
      int GID = GlobalEntries[i];
      if (GID<MinAllGID) ierr = 1;
      else if (GID>MaxAllGID) ierr = 1;
      else {
	// Guess uniform distribution and start a little above it
	int Proc1 = EPETRA_MIN(GID/EPETRA_MAX(n_over_p,1) + 2, NumProc-1);
	bool found = false;
	while (Proc1 >= 0 && Proc1< NumProc) {
	  if (AllMinGIDs_[Proc1]<=GID) {
	    if (GID <AllMinGIDs_[Proc1+1]) {
	    found = true;
	    break;
	    }
	    else Proc1++;
	  }
	  else Proc1--;
	}
	if (found) {
	  Proc = Proc1;
	  LID = GID - AllMinGIDs_[Proc];
	}
      }
      Procs[i] = Proc;
      if (LocalEntries!=0) LocalEntries[i] = LID;
    }
    if (EntrySizes!=0) {
      if (Map_->ConstantElementSize()) {
	int ElementSize = Map_->MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
      }
      else {
	int * ElementSizeList = Map_->ElementSizeList(); // We know this exists
	
	
	Epetra_Distributor * Size_Distor = Map_->Comm().CreateDistributor();
	
	int Size_num_sends;
	int * Size_send_gids = 0;
	int * Size_send_procs = 0;

	
	EPETRA_CHK_ERR(Size_Distor->CreateFromRecvs( NumEntries, GlobalEntries, Procs, true,
						       Size_num_sends, Size_send_gids, Size_send_procs ));
	
	int * Size_exports = 0;
	int * Size_imports = 0;
	if (Size_num_sends>0) {
	  Size_exports = new int[ 2 * Size_num_sends ];
	  for( i = 0; i < Size_num_sends; i++ )
	    {
	      int Size_curr_GID = Size_send_gids[i];
	      int Size_curr_LID = Map_->LID(Size_curr_GID);
	      assert(Size_curr_LID!=-1); // Internal error 
	      Size_exports[2*i] = Size_curr_GID;
	      int Size_curr_size = ElementSizeList[Size_curr_LID];
	      Size_exports[2*i+1] = Size_curr_size;
	    }
	}
	
	if (NumEntries>0) Size_imports = new int[ 2 * NumEntries ];
	EPETRA_CHK_ERR(Size_Distor->Do( reinterpret_cast<char*> (Size_exports),
			 2 * sizeof( int ), reinterpret_cast<char*> (Size_imports)));
	
	for( i = 0; i < NumEntries; i++ )
	  {

	    // Need to change !!!!
	    //bool found = false;
	    int Size_curr_LID = Size_imports[2*i];
	    for( j = 0; j < NumEntries; j++ )
	      if( Size_curr_LID == GlobalEntries[j] )
		{
		  EntrySizes[j] = Size_imports[2*i+1];
		  // found = true;
		  break;
		}
	    //	if (!found) cout << "Internal error:  Epetra_Directory::GetDirectoryEntries: Global Index " << curr_LID
	    //	     << " not on processor " << MyPID << endl; abort();
	  }
	
	if( Size_send_gids != 0 ) delete [] Size_send_gids;
	if( Size_send_procs != 0 ) delete [] Size_send_procs;
	
	if( Size_imports != 0 ) delete [] Size_imports;
	if( Size_exports != 0 ) delete [] Size_exports;
	
	delete Size_Distor;
      }
    }
    EPETRA_CHK_ERR(ierr);
    return(0);
  }

  // General case (need to set up an actual directory structure)
  
  int * ElementSizeList = 0;
  int PacketSize = 2; // We will send at least the GID and PID.  Might also send LID and Size info
  bool DoSizes = false;
  if (EntrySizes!=0) {
    if (Map_->ConstantElementSize()) {
      int ElementSize = Map_->MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
    }
    else {
      ElementSizeList = Map_->ElementSizeList(); // We know this exists
      DoSizes = true;
      PacketSize++; // Sending Size info
    }
  }

  bool DoLIDs = (LocalEntries!=0); // Do LIDs?
  if (DoLIDs) PacketSize++; // Sending LIDs also

  
  Epetra_Distributor * Distor = DirectoryMap_->Comm().CreateDistributor();
  
  
  int * dir_procs = 0;
  if (NumEntries>0) dir_procs = new int[ NumEntries ];
  
  // Get directory locations for the requested list of entries
  DirectoryMap_->RemoteIDList(NumEntries, GlobalEntries, dir_procs, 0);
  int num_sends;
  int * send_gids = 0;
  int * send_procs = 0;
  
  EPETRA_CHK_ERR(Distor->CreateFromRecvs( NumEntries, GlobalEntries, dir_procs, true,
					   num_sends, send_gids, send_procs));

  if (NumEntries>0) delete [] dir_procs;


  int curr_LID;
  int * exports = 0;
  int * imports = 0;
  if (num_sends>0) {
    exports = new int[ PacketSize * num_sends ];
    int * ptr = exports;
    for( i = 0; i < num_sends; i++ )
      {
	int curr_GID = send_gids[i];
	*ptr++ = curr_GID;
	curr_LID = DirectoryMap_->LID(curr_GID);
	assert(curr_LID!=-1); // Internal error 
	*ptr++ = ProcList_[ curr_LID ];
	if (DoLIDs) *ptr++ = LocalIndexList_[curr_LID];
	if (DoSizes) *ptr++ = ElementSizeList[curr_LID];
      }
  }

  if (NumEntries>0) imports = new int[PacketSize*NumEntries];
  EPETRA_CHK_ERR(Distor->Do(reinterpret_cast<char*> (exports),
	     PacketSize * sizeof( int ), reinterpret_cast<char*> (imports)));

  int * ptr = imports;
  for( i = 0; i < NumEntries; i++ )
    {
      //bool found = false;
      curr_LID = *ptr++;
      for( j = 0; j < NumEntries; j++ )
	if( curr_LID == GlobalEntries[j] )
	  {
	    Procs[j] = *ptr++;
	    if (DoLIDs) LocalEntries[j] = *ptr++;
	    if (DoSizes) EntrySizes[j] = *ptr++;
	    // found = true;
	    break;
	  }
      //	if (!found) cout << "Internal error:  Epetra_Directory::GetDirectoryEntries: Global Index " << curr_LID
      //	     << " not on processor " << MyPID << endl; abort();
    }
  
  if( send_gids != 0 ) delete [] send_gids;
  if( send_procs != 0 ) delete [] send_procs;
  
  if( imports != 0 ) delete [] imports;
  if( exports != 0 ) delete [] exports;
  
  delete Distor;
  return(0);
}

//==============================================================================
void Epetra_Directory::Print( ostream & os) const {
  
  int MyPID;
  if( DirectoryMap_ != 0 )
  {
    MyPID = DirectoryMap_->Comm().MyPID();
    os << MyPID << " Epetra_Directory Object: "
      << DirectoryMap_->NumMyElements() << endl;
    for( int i = 0; i < DirectoryMap_->NumMyElements(); i++ )
      os << " " << i << " " << ProcList_[i] << " "
        << LocalIndexList_[i] << endl;
     os << endl;
  }
  else
  {
    cout << "Epetra_Directory not setup<<<<<<" << endl;
  }

  return;
}
