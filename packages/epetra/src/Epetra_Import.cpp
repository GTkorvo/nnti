
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


#include "Epetra_Import.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Distributor.h"
#include "Epetra_Comm.h"


//==============================================================================
// Epetra_Import constructor for a Epetra_BlockMap object
Epetra_Import::Epetra_Import( const Epetra_BlockMap &  TargetMap, const Epetra_BlockMap & SourceMap)
  : Epetra_Object("Epetra::Import"),
    TargetMap_(TargetMap),
    SourceMap_(SourceMap),
    NumSameIDs_(0),
    NumPermuteIDs_(0),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(0),
    RemoteLIDs_(0),
    NumExportIDs_(0),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(0),
    NumRecv_(0),
    Distor_(0)
{

  int i;
  
  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumRemoteIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be imported.
  
  int NumSourceIDs = SourceMap.NumMyElements();
  int NumTargetIDs = TargetMap.NumMyElements();
  
  int *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int[NumTargetIDs];
    TargetMap.MyGlobalElements(TargetGIDs);
  }
  
  int * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int[NumSourceIDs];
    SourceMap.MyGlobalElements(SourceGIDs);
  }
  
  int MinIDs = EPETRA_MIN(NumSourceIDs, NumTargetIDs);
  
  
  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;
  
  
  // Find count of Target IDs that are truly remote and those that are local but permuted

  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) 
    if (SourceMap.MyGID(TargetGIDs[i])) NumPermuteIDs_++; // Check if Target GID is a local Source GID
    else NumRemoteIDs_++; // If not, then it is remote
  
  
  
  // Define remote and permutation lists
  
  int * RemoteGIDs;
  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    RemoteGIDs = new int[NumRemoteIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }
  
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) {
    if (SourceMap.MyGID(TargetGIDs[i])) {
      PermuteToLIDs_[NumPermuteIDs_] = i;
      PermuteFromLIDs_[NumPermuteIDs_++] = SourceMap.LID(TargetGIDs[i]);
    }
    else {
      //NumRecv_ +=TargetMap.ElementSize(i); // Count total number of entries to receive
      NumRecv_ +=TargetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
      RemoteGIDs[NumRemoteIDs_] = TargetGIDs[i];
      RemoteLIDs_[NumRemoteIDs_++] = i;
    }
  }
  
  if ( NumRemoteIDs_>0 && !SourceMap.DistributedGlobal())
    throw ReportError("Error in Epetra_Export: Serial Export has remote IDs.", -1);
  
  // Test for distributed cases
  
  if (SourceMap.DistributedGlobal()) {
    
    int * RemotePIDs = 0;
    
    if (NumRemoteIDs_>0)  RemotePIDs = new int[NumRemoteIDs_];
    int ierr = SourceMap.RemoteIDList(NumRemoteIDs_, RemoteGIDs, RemotePIDs, 0); // Get remote PIDs
    if (ierr!=0) throw ReportError("Error in SourceMap.RemoteIDList call", ierr);
    Distor_ = SourceMap.Comm().CreateDistributor();
    
    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.
    
    bool Deterministic = true;
    ierr = Distor_->CreateFromRecvs( NumRemoteIDs_, RemoteGIDs,
					 RemotePIDs, Deterministic,
					 NumExportIDs_, ExportLIDs_, ExportPIDs_ );
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromRecvs()", ierr);
    
    // Use comm plan with Export GIDs (stored in ExportLIDs_) to
    // get proper ordering of GIDs for remote entries 
    // (that we will convert to LIDs when done).
    
    ierr = Distor_->Do( reinterpret_cast<char *> (ExportLIDs_), 
		sizeof( int ),
		reinterpret_cast<char *> (RemoteGIDs));
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.Do()", ierr);

        
    // Export IDs come in as GIDs, convert to LIDs
    for (i=0; i< NumExportIDs_; i++) {
      if (ExportPIDs_[i] < 0) throw ReportError("TargetMap requested a GID that is not in the SourceMap.", -1);
      
      ExportLIDs_[i] = SourceMap.LID(ExportLIDs_[i]);
      //NumSend_ += SourceMap.ElementSize(ExportLIDs_[i]); // Count total number of entries to send
      NumSend_ += SourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
    }
    
    // Remote IDs come in as GIDs, convert to LIDs in proper order

    // for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = TargetMap.LID(RemoteGIDs[i]); // Only works when target map has no repeated GIDs
    
    if (NumRemoteIDs_>0) {
      int * ReorderedRemoteLIDs = RemotePIDs; // Reuse some temp space
      for (i=0; i< NumRemoteIDs_; i++) {
	int CurrentGID = RemoteGIDs[i];
	bool Found = false;
	for (int j=0; j < NumRemoteIDs_; j++) {
	  if (RemoteLIDs_[j]!= -1) {
	    if (CurrentGID==TargetGIDs[RemoteLIDs_[j]]) {
	      ReorderedRemoteLIDs[i] = RemoteLIDs_[j];
	      RemoteLIDs_[j] = -1;
	      Found = true;
	      break;
	    }
	  }
	}
	if (!Found) throw ReportError("Internal error.  Cannot map incoming GID to Target Map", -2);
      }
      
      // Clean up and leave....
      
      delete [] RemoteLIDs_;
      delete [] RemoteGIDs;
      RemoteLIDs_ = ReorderedRemoteLIDs;
    }
  }

  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;
  
  return;
}

//==============================================================================
// Epetra_Import copy constructor 
Epetra_Import::Epetra_Import(const Epetra_Import & Importer)
  : Epetra_Object(Importer),
    TargetMap_(Importer.TargetMap_),
    SourceMap_(Importer.SourceMap_),
    NumSameIDs_(Importer.NumSameIDs_),
    NumPermuteIDs_(Importer.NumPermuteIDs_),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(Importer.NumRemoteIDs_),
    RemoteLIDs_(0),
    NumExportIDs_(Importer.NumExportIDs_),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(Importer.NumSend_),
    NumRecv_(Importer.NumRecv_),
    Distor_(0)
{
  int i;
  if (NumPermuteIDs_>0) {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (i=0; i< NumPermuteIDs_; i++) {
      PermuteToLIDs_[i] = Importer.PermuteToLIDs_[i];
      PermuteFromLIDs_[i] = Importer.PermuteFromLIDs_[i];
    }
  }

  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = Importer.RemoteLIDs_[i];
  }

  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (i=0; i< NumExportIDs_; i++) {
      ExportLIDs_[i] = Importer.ExportLIDs_[i];
      ExportPIDs_[i] = Importer.ExportPIDs_[i];
    }
  }

  if (Importer.Distor_!=0) Distor_ = Importer.Distor_->Clone();

}

//==============================================================================
// Epetra_Import destructor 
Epetra_Import::~Epetra_Import()
{
  if( Distor_ != 0 ) delete Distor_;
  if (RemoteLIDs_ != 0) delete [] RemoteLIDs_;
  if (PermuteToLIDs_ != 0) delete [] PermuteToLIDs_;
  if (PermuteFromLIDs_ != 0) delete [] PermuteFromLIDs_;

  if( ExportPIDs_ != 0 ) delete [] ExportPIDs_; // These were created by Distor_
  if( ExportLIDs_ != 0 ) delete [] ExportLIDs_;
}
//=============================================================================
void Epetra_Import::Print(ostream & os) const
{

  os << endl << endl << "Source Map:" << endl << endl;
  SourceMap_.Print(os);
  
  os << endl << endl << "Target Map:" << endl << endl;
  TargetMap_.Print(os);
  
  os << endl << endl << "Distributor:" << endl << endl;
  if (Distor_==0) os << "  Is empty...." << endl;
  else Distor_->Print(os);
  
  os << "Number of Same IDs = " << NumSameIDs_ << endl;
  
  os << "Epetra_Import Print Needs attention!!!!" << endl;
  return;
}

