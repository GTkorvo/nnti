
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


#include "Epetra_BlockMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Directory.h"


//==============================================================================
// Epetra_BlockMap constructor for a Epetra-defined uniform linear distribution of constant block size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    LID_(0),
    NumGlobalElements_(NumGlobalElements),
    MyGlobalElements_(0),
    FirstElementEntryList_(0),
    ElementSize_(ElementSize),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    MinMyElementSize_(ElementSize),
    MaxMyElementSize_(ElementSize),
    MinElementSize_(ElementSize),
    MaxElementSize_(ElementSize),
    ConstantElementSize_(true),
    LinearMap_(true)
{
  // Each processor gets roughly numGlobalEquations/p equations
  // This routine automatically defines a linear partitioning of a
  // map with numGlobalEquations across the processors
  // specified in the given Epetra_Comm

    if (NumGlobalElements_ < 0) 
      throw ReportError("NumGlobalElements = " + toString(NumGlobalElements_) + ".  Should be >= 0.", -1);
   if (ElementSize_ <= 0) 
     throw ReportError("ElementSize = " + toString(ElementSize_) + ".  Should be > 0.", -2);

  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();
  NumMyElements_ = NumGlobalElements_/NumProc;
  int remainder = NumGlobalElements_%NumProc;
  int start_index = MyPID*(NumMyElements_ + 1);

  if (MyPID<remainder) NumMyElements_++;
  else start_index -= (MyPID-remainder);

  NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
  NumMyEquations_ = NumMyElements_ * ElementSize_;

  MinMyElementSize_ = ElementSize_;
  MaxMyElementSize_ = ElementSize_;
  MinElementSize_ = ElementSize_;
  MaxElementSize_ = ElementSize_;

  MinAllGID_ = IndexBase_;
  MaxAllGID_ = MinAllGID_ + NumGlobalElements_ - 1;
  MinMyGID_ = start_index + IndexBase_;
  MaxMyGID_ = MinMyGID_ + NumMyElements_ - 1;
  MinLID_ = 0;
  MaxLID_ = MinLID_ + NumMyElements_ - 1;
  DistributedGlobal_ = IsDistributedGlobal(NumGlobalElements_, NumMyElements_);

  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}
//==============================================================================
// Epetra_BlockMap constructor for a user-defined linear distribution of constant block size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, 
			       int ElementSize, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    LID_(0),
    NumGlobalElements_(NumGlobalElements),
    NumMyElements_(NumMyElements),
    MyGlobalElements_(0),
    FirstElementEntryList_(0),
    ElementSize_(ElementSize),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    MinMyElementSize_(ElementSize),
    MaxMyElementSize_(ElementSize),
    MinElementSize_(ElementSize),
    MaxElementSize_(ElementSize),
    ConstantElementSize_(true),
    LinearMap_(true)
{
  // Each processor gets NumMyElements equations
  
  if (NumGlobalElements_ < -1) 
    throw ReportError("NumGlobalElements = " + toString(NumGlobalElements_) + ".  Should be >= -1.", -1);
  if (NumMyElements_ < 0) 
    throw ReportError("NumMyElements = " + toString(NumMyElements_) + ".  Should be >= 0.", -2);
  if (ElementSize_ <= 0) 
    throw ReportError("ElementSize = " + toString(ElementSize_) + ". Should be > 0.", -3);

  // Get processor information

  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();

  DistributedGlobal_ = IsDistributedGlobal(NumGlobalElements, NumMyElements);

   // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!DistributedGlobal_ || NumProc==1)
    {
      NumGlobalElements_ = NumMyElements_;
      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(NumGlobalElements) + 
			  ".  Should equal " + toString(NumGlobalElements_) + 
			  ", or be set to -1 to compute automatically", -4);

      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;
      
      MinAllGID_ = IndexBase_;
      MaxAllGID_ = MinAllGID_ + NumGlobalElements_ - 1;
      MinMyGID_ = IndexBase_;
      MaxMyGID_ = MinMyGID_ + NumMyElements_ - 1;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;
    }
  else if (NumProc > 1)
    {
      // Sum up all local element counts to get global count
     Comm_->SumAll(&NumMyElements_, &NumGlobalElements_, 1);

      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(NumGlobalElements) + 
			  ".  Should equal " + toString(NumGlobalElements_) + 
			  ", or be set to -1 to compute automatically", -4);

      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;

      MinAllGID_ = IndexBase_;
      MaxAllGID_ = MinAllGID_ + NumGlobalElements_ - 1;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;

      // Use the ScanSum function to compute a prefix sum of the number of equations
      Comm_->ScanSum(&NumMyElements_, &MaxMyGID_, 1);

      int start_index = MaxMyGID_ - NumMyElements_;
      MinMyGID_ = start_index + IndexBase_;
      MaxMyGID_ = MinMyGID_ + NumMyElements_ - 1;
    }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}
//==============================================================================
// Epetra_BlockMap constructor for a user-defined arbitrary distribution of constant block size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
			       int ElementSize, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    LID_(0),
    NumGlobalElements_(NumGlobalElements),
    NumMyElements_(NumMyElements),
    MyGlobalElements_(MyGlobalElements),
    FirstElementEntryList_(0),
    ElementSize_(ElementSize),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    MinMyElementSize_(ElementSize),
    MaxMyElementSize_(ElementSize),
    MinElementSize_(ElementSize),
    MaxElementSize_(ElementSize),
    ConstantElementSize_(true),
    LinearMap_(false)
{
  int i;
  // Each processor gets NumMyElements equations

  if (NumGlobalElements_ < -1) 
    throw ReportError("NumGlobalElements = " + toString(NumGlobalElements_) + ".  Should be >= -1.", -1);
  if (NumMyElements_ < 0) 
    throw ReportError("NumMyElements = " + toString(NumMyElements_) + ".  Should be >= 0.", -2);
  if (ElementSize_ <= 0) 
    throw ReportError("ElementSize = " + toString(ElementSize_) + ". Should be > 0.", -3);

  // Allocate storage for global index list information

  if (NumMyElements>0) MyGlobalElements_ = new int[NumMyElements];

  // Get processor information

  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();
  if (NumMyElements>0) {
    // Compute min/max GID on this processor
    MinMyGID_ = MyGlobalElements[0];
    MaxMyGID_ = MyGlobalElements[0];
    for (i = 0; i < NumMyElements; i++)
      {
	MyGlobalElements_[i] = MyGlobalElements[i];
	MinMyGID_ = EPETRA_MIN(MinMyGID_,MyGlobalElements[i]);
	MaxMyGID_ = EPETRA_MAX(MaxMyGID_,MyGlobalElements[i]);
      }
  }
  else {
    MinMyGID_ = IndexBase_;
    MaxMyGID_ = IndexBase_;
  }
    
  DistributedGlobal_ = IsDistributedGlobal(NumGlobalElements, NumMyElements);

    // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!DistributedGlobal_ || NumProc==1)
    {
      NumGlobalElements_ = NumMyElements_;
      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(NumGlobalElements) + 
			  ".  Should equal " + toString(NumGlobalElements_) + 
			  ", or be set to -1 to compute automatically", -4);
      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;
      
      MinAllGID_ = MinMyGID_;
      MaxAllGID_ = MaxMyGID_;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;
    }
  else if (NumProc > 1)
    {
      // Sum up all local element counts to get global count
      Comm_->SumAll(&NumMyElements_, &NumGlobalElements_, 1);
      // Check to see if user's value for NumGlobalElements is either -1
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(NumGlobalElements) + 
			  ".  Should equal " + toString(NumGlobalElements_) + 
			  ", or be set to -1 to compute automatically", -4);
      
      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;
      
      MinLID_ = 0;
      MaxLID_ = EPETRA_MAX(MinLID_ + NumMyElements_ - 1,MinLID_);
      
      // Use the Allreduce function to find min/max GID 
      int *tmp_send = new int[2];
      int *tmp_recv = new int[2];
      tmp_send[0] = - MinMyGID_; // Negative sign lets us do one reduction
      tmp_send[1] =   MaxMyGID_;
      Comm_->MaxAll(tmp_send, tmp_recv, 2);
      MinAllGID_ = - tmp_recv[0];
      MaxAllGID_ =   tmp_recv[1];
      delete [] tmp_send;
      delete [] tmp_recv;
      if (MinAllGID_ < IndexBase_)
	throw ReportError("Minimum global element index = " + toString(MinAllGID_) + " is less than index base = " + toString(IndexBase_) +".", -5);
    }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}

//==============================================================================
// Epetra_BlockMap constructor for a user-defined arbitrary distribution of variable block size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
			       int *ElementSizeList, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    LID_(0),
    NumGlobalElements_(NumGlobalElements),
    NumMyElements_(NumMyElements),
    FirstElementEntryList_(0),
    ElementSize_(0),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    ConstantElementSize_(false),
    LinearMap_(false)
{
  int i;
  // Each processor gets NumMyElements equations

  if (NumGlobalElements_ < -1) 
    throw ReportError("NumGlobalElements = " + toString(NumGlobalElements_) + ".  Should be >= -1.", -1);
  if (NumMyElements_ < 0) 
    throw ReportError("NumMyElements = " + toString(NumMyElements_) + ".  Should be >= 0.", -2);
  for (i=0; i<NumMyElements_; i++)
    if (ElementSizeList[i] <= 0) 
      throw ReportError("An entry in ElementSizeList = " + toString(ElementSizeList[i]) + ". Should be > 0.", -3);
  
  // Allocate storage for global index list and element size information

  if (NumMyElements>0) {
    MyGlobalElements_ = new int[NumMyElements];
    ElementSizeList_ = new int[NumMyElements];
  }
  // Get processor information

  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();
  
  if (NumMyElements>0) {
    // Compute min/max GID and element size, number of equations on this processor
    MinMyGID_ = MyGlobalElements[0];
    MaxMyGID_ = MyGlobalElements[0];
    MinMyElementSize_ = ElementSizeList[0];
    MaxMyElementSize_ = ElementSizeList[0];
    NumMyEquations_ = 0;
    for (i = 0; i < NumMyElements; i++)
      {
	MyGlobalElements_[i] = MyGlobalElements[i];
	ElementSizeList_[i] = ElementSizeList[i];
	MinMyGID_ = EPETRA_MIN(MinMyGID_,MyGlobalElements[i]);
	MaxMyGID_ = EPETRA_MAX(MaxMyGID_,MyGlobalElements[i]);
	MinMyElementSize_ = EPETRA_MIN(MinMyElementSize_,ElementSizeList[i]);
	MaxMyElementSize_ = EPETRA_MAX(MaxMyElementSize_,ElementSizeList[i]);
	NumMyEquations_ += ElementSizeList[i];
      }
  }
  else {
    MinMyGID_ = IndexBase_;
    MaxMyGID_ = IndexBase_;
    MinMyElementSize_ = 1;
    MaxMyElementSize_ = 1;
    NumMyEquations_ = 0;
  }

  DistributedGlobal_ = IsDistributedGlobal(NumGlobalElements, NumMyElements);  
  
  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!DistributedGlobal_ || NumProc==1)
    {
      NumGlobalElements_ = NumMyElements_;
      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(NumGlobalElements) + 
			  ".  Should equal " + toString(NumGlobalElements_) + 
			  ", or be set to -1 to compute automatically", -4);
      NumGlobalEquations_ = NumMyEquations_;
      
      MinAllGID_ = MinMyGID_;
      MaxAllGID_ = MaxMyGID_;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;
      MinElementSize_ = MinMyElementSize_;
      MaxElementSize_ = MaxMyElementSize_;
    }
  else if (NumProc > 1)
    {
      // Sum up all local element and equation counts to get global counts
      int *tmp_send = new int[4];
      int *tmp_recv = new int[4];
      tmp_send[0] = NumMyElements_;
      tmp_send[1] = NumMyEquations_;
      Comm_->SumAll(tmp_send, tmp_recv, 2);
      NumGlobalElements_ =  tmp_recv[0];
      NumGlobalEquations_ = tmp_recv[1];

      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(NumGlobalElements) + 
			  ".  Should equal " + toString(NumGlobalElements_) + 
			  ", or be set to -1 to compute automatically", -4);
      
      MinLID_ = 0;
      MaxLID_ = EPETRA_MAX(MinLID_ + NumMyElements_ - 1,MinLID_);
      
      // Use the MaxAll function to find min/max GID 
      tmp_send[0] = - MinMyGID_; // Negative signs lets us do one reduction
      tmp_send[1] =   MaxMyGID_;
      tmp_send[2] = - MinMyElementSize_;
      tmp_send[3] =   MaxMyElementSize_;

      Comm_->MaxAll(tmp_send, tmp_recv, 4);

      MinAllGID_ =      - tmp_recv[0];
      MaxAllGID_ =        tmp_recv[1];
      MinElementSize_ = - tmp_recv[2];
      MaxElementSize_ =   tmp_recv[3];

      delete [] tmp_send;
      delete [] tmp_recv;

      if (MinAllGID_ < IndexBase_)
	throw ReportError("Minimum global element index = " + toString(MinAllGID_) + " is less than index base = " + toString(IndexBase_) +".", -5);
    }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);

  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}
//==============================================================================
Epetra_BlockMap::Epetra_BlockMap(const Epetra_BlockMap& map)
  : Epetra_Object(map.Label()),
    LID_(0),
    NumGlobalElements_(map.NumGlobalElements_),
    NumMyElements_(map.NumMyElements_),
    MyGlobalElements_(0),
    FirstElementEntryList_(0),
    ElementSize_(map.ElementSize_),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(map.IndexBase_),
    Comm_(map.Comm_),
    Directory_(0),
    NumGlobalEquations_(map.NumGlobalEquations_),
    NumMyEquations_(map.NumMyEquations_),
    MinAllGID_(map.MinAllGID_),
    MaxAllGID_(map.MaxAllGID_),
    MinMyGID_(map.MinMyGID_),
    MaxMyGID_(map.MaxMyGID_),
    MinLID_(map.MinLID_),
    MaxLID_(map.MaxLID_),
    MinMyElementSize_(map.MinMyElementSize_),
    MaxMyElementSize_(map.MaxMyElementSize_),
    MinElementSize_(map.MinElementSize_),
    MaxElementSize_(map.MaxElementSize_),
    ConstantElementSize_(map.ConstantElementSize_),
    LinearMap_(map.LinearMap_),
    DistributedGlobal_(map.DistributedGlobal_)
{
  int i;

  if (map.MyGlobalElements_!=0)
    {
      MyGlobalElements_ = new int[NumMyElements_];
      
      for(i=0; i<NumMyElements_; i++)
	MyGlobalElements_[i] = map.MyGlobalElements_[i];
    }
  if (map.FirstElementEntryList_!=0)
    {
      FirstElementEntryList_ = new int[NumMyElements_+1];
      
      for(i=0; i<NumMyElements_+1; i++)
	FirstElementEntryList_[i] = map.FirstElementEntryList_[i];
    }
  if (map.ElementSizeList_!=0)
    {
      ElementSizeList_ = new int[NumMyElements_];
      
      for(i=0; i<NumMyElements_; i++)
	ElementSizeList_[i] = map.ElementSizeList_[i];
    }
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}

//==============================================================================
Epetra_BlockMap::~Epetra_BlockMap(void)  {

  if (MyGlobalElements_ != 0 && NumMyElements_>0) delete [] MyGlobalElements_;
  MyGlobalElements_ = 0;
  
  if (FirstElementEntryList_ != 0 && NumMyElements_>0) delete [] FirstElementEntryList_;
  FirstElementEntryList_ = 0;
  
  if (ElementSizeList_ != 0 && NumMyElements_>0) delete [] ElementSizeList_;
  ElementSizeList_ = 0;

  if (EquationToBlockList_ != 0 && NumMyEquations_>0) delete [] EquationToBlockList_;
  EquationToBlockList_ = 0;

  if (Directory_ !=0) delete Directory_;
  Directory_ = 0;

  if (LID_ !=0 && NumMyElements_>0) delete [] LID_;
  LID_ = 0;
}


//==============================================================================
bool Epetra_BlockMap::SameAs(const Epetra_BlockMap & Map) const
{

  if (this == &Map) return(true);
    
  if (MinAllGID_ != Map.MinAllGID() ||
      MaxAllGID_ != Map.MaxAllGID() ||
      NumGlobalElements_!=Map.NumGlobalElements() ||
      IndexBase_!=Map.IndexBase() ) return(false);
  
  if (ConstantElementSize_) {
    if (ElementSize_!=Map.ElementSize()) return(false);
    else return(true);
  }
  else {

    // If we get this far, we need to check local properties and then check across
    // all processors to see if local properties are all true

    int MySameMap = 1; // Assume not needed
    if (NumMyElements_!=Map.NumMyElements()) MySameMap = 0;
    
    if (MySameMap==1) 
      for (int i=0; i<NumMyElements_; i++) 
	if (GID(i) != Map.GID(i)) MySameMap = 0;

    // Now get min of MySameMap across all processors

    int GlobalSameMap = 0;
    assert(Comm().MinAll(&MySameMap, &GlobalSameMap, 1)==0);
    
    return(GlobalSameMap==1);
  }
  return(false);
}


//==============================================================================
int Epetra_BlockMap::MyGlobalElements(int * MyGlobalElements) const
{
  // If the global element list is not create, then do so.  This can only happen when
  // a linear distribution has been specified.  Thus we can easily construct the update
  // list in this case.

  int i;
  if (MyGlobalElements_==0)
      for (i = 0; i<NumMyElements_; i++)
	MyGlobalElements[i] = MinMyGID_ + i;
  else
    for (i = 0; i<NumMyElements_; i++)
      MyGlobalElements[i] = MyGlobalElements_[i];
  
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::MyGlobalElements() const {
  
  // If ElementSizeList not built, do so
  if (MyGlobalElements_==0 && NumMyElements_>0) {
    int * tmp = new int[NumMyElements_+1];
    MyGlobalElements(tmp);
    (int * &) MyGlobalElements_ = tmp;
  }
  return(MyGlobalElements_);
  
}
//==============================================================================
int Epetra_BlockMap::FirstElementEntryList(int * FirstElementEntryList) const
{
  // If the first element entry list is not create, then do so.  

  // Note: This array is of length NumMyElement+1

  int i;

  if (FirstElementEntryList_==0) {
    FirstElementEntryList[0] = 0; // First element of first entry is always zero
    
    if (ConstantElementSize_)
      for (i = 0; i<NumMyElements_; i++)
	FirstElementEntryList[i+1] = FirstElementEntryList[i] + ElementSize_;
    else
      for (i = 0; i<NumMyElements_; i++)
	FirstElementEntryList[i+1] = FirstElementEntryList[i] + ElementSizeList_[i];
  }
  else 
    for (i = 0; i<=NumMyElements_; i++)
      FirstElementEntryList[i] = FirstElementEntryList_[i];
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::FirstElementEntryList() const {

  // If ElementSizeList not built, do so
  if (FirstElementEntryList_==0 && NumMyElements_>0) {
    int * tmp = new int[NumMyElements_+1];
    FirstElementEntryList(tmp);
    (int * &) FirstElementEntryList_ = tmp;
 }
  return(FirstElementEntryList_);
  
}
//==============================================================================
int Epetra_BlockMap::ElementSizeList(int * ElementSizeList) const
{
  // If the element size list is not create, then do so.  This can only happen when
  // a constant element size has been specified.  Thus we can easily construct the block size
  // list in this case.

  int i;
  if (ElementSizeList_==0)
    for (i = 0; i<NumMyElements_; i++)
      ElementSizeList[i] = ElementSize_;
  else
    for (i = 0; i<NumMyElements_; i++)
      ElementSizeList[i] = ElementSizeList_[i];
  
  return(0);
}
//==============================================================================
int * Epetra_BlockMap::ElementSizeList() const {

  // If ElementSizeList not built, do so
  if (ElementSizeList_==0 && NumMyElements_>0) {
    int * tmp = new int[NumMyElements_];
    ElementSizeList(tmp);
    (int * &) ElementSizeList_ = tmp;
 }
  return(ElementSizeList_);
  
}
//==============================================================================
int Epetra_BlockMap::EquationToBlockList(int * EquationToBlockList) const
{
  // Build an array such that the local block ID is stored for each equation

  int i;
  if (EquationToBlockList_==0) {
    int * ptr = EquationToBlockList;
    for (i = 0; i<NumMyElements_; i++) {
      int Size = ElementSize(i);
      for (int j=0; j<Size; j++) *ptr++ = i;
    }
  }
  else
    for (i = 0; i<NumMyEquations_; i++)
      EquationToBlockList[i] = EquationToBlockList_[i];
  
  return(0);
}
//==============================================================================
int * Epetra_BlockMap::EquationToBlockList() const {

  // If EquationToBlockList not built, do so
  if (EquationToBlockList_==0 && NumMyEquations_>0) {
    int * tmp = new int[NumMyEquations_];
    EquationToBlockList(tmp);
    (int * &) EquationToBlockList_ = tmp;
 }
  return(EquationToBlockList_);
  
}
//==============================================================================
int Epetra_BlockMap::ElementSize(int LID) const {

  if (ConstantElementSize()) return(ElementSize_);
  else
    return(ElementSizeList_[LID]);
  
}
//==============================================================================
void Epetra_BlockMap::GlobalToLocalSetup() {

  int i;

  if (NumGlobalElements_==0) return; // Nothing to do

  else if (LinearMap() || (!DistributedGlobal()) || NumMyElements_==0) {
    if (Directory_ ==0) Directory_ = new Epetra_Directory(this); // Make directory
    return; // Nothing else to do
  }
  else {
    // Build LID_ vector to make look up of local index values fast
    
    int SpanGID = MaxMyGID_ - MinMyGID_ + 1;
    LID_ = new int[SpanGID];
    
    for (i=0; i<SpanGID; i++) LID_[i] = -1; // Fill all locations with -1
    
    for (i=0; i<NumMyElements_; i++) {
      int tmp = MyGlobalElements_[i]-MinMyGID_;
      assert(tmp>=0); assert(tmp <SpanGID);
      LID_[MyGlobalElements_[i]-MinMyGID_] = i; // Spread local indices
    }
    
    if (Directory_ ==0) Directory_ = new Epetra_Directory(this); // Make directory
  }
}

//==============================================================================
int Epetra_BlockMap::LID(int GID) const {

  if (GID<MinMyGID_ || GID > MaxMyGID_) return(-1); // Out of range
  else if (!DistributedGlobal()) return(GID-IndexBase_); // I own all indices
  else if (LinearMap()) return(GID-MinMyGID_); // Can compute with an offset
  else return(LID_[GID-MinMyGID_]); // Find it in LID array
}
//==============================================================================
int Epetra_BlockMap::GID(int LID) const {

  if (LID<MinLID_ || LID>MaxLID_) return(IndexBase_-1); // Out of range
  else if (!DistributedGlobal()) return(LID+IndexBase_); // I own all indices
  else if (LinearMap()) return(LID+MinMyGID_); // Can compute with an offset
  else return(MyGlobalElements_[LID]); // Find it in MyGlobalElements array
}
//==============================================================================
int Epetra_BlockMap::FindLocalBlockID(int EquationID, int & BlockID, int & BlockOffset) const {

  if (EquationID>=NumMyEquations_) return(-1); // Equation is out of range

  if (ConstantElementSize()) {
    BlockID = EquationID/MaxElementSize_;
    BlockOffset = EquationID%MaxElementSize_;
    return(0);
  }
  else {
    if (EquationToBlockList_==0) {
      int * tmp;
      if (EquationToBlockList_==0) tmp = new int[NumMyEquations_];
      EPETRA_CHK_ERR(EquationToBlockList(tmp));
      (int * &) EquationToBlockList_ = tmp;  // Allows assignment in a const method
    }
    if (FirstElementEntryList_==0) {
      int * tmp;
      if (FirstElementEntryList_==0) tmp = new int[NumMyElements_];
      EPETRA_CHK_ERR(FirstElementEntryList(tmp));
      (int * &) FirstElementEntryList_ = tmp;  // Allows assignment in a const method
    }

    BlockID = EquationToBlockList_[EquationID];
    BlockOffset = EquationID - FirstElementEntryList_[BlockID];
    return(0);
  }
}
//==============================================================================
int Epetra_BlockMap::RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList, int * SizeList) const {

  EPETRA_CHK_ERR(Directory_->GetDirectoryEntries(NumIDs, GIDList, PIDList, LIDList, SizeList));
  return(0);
}
//==============================================================================
bool Epetra_BlockMap::IsDistributedGlobal(int NumGlobalElements, int NumMyElements) const {

  bool DistributedGlobal = false; // Assume map is not global distributed
  if (Comm_->NumProc()>1) {
    int LocalReplicated = 0;
    int AllLocalReplicated;
    if (NumGlobalElements==NumMyElements) LocalReplicated=1;
    Comm_->MinAll(&LocalReplicated, &AllLocalReplicated, 1);

    // If any PE has LocalReplicated=0, then map is distributed global
    if (AllLocalReplicated!=1) DistributedGlobal = true;
  }
  return(DistributedGlobal);
}
//==============================================================================
void Epetra_BlockMap::Print(ostream & os) const
{
  int * MyGlobalElements1 = MyGlobalElements();
  int * FirstElementEntryList1 = 0;
  int * ElementSizeList1 = 0;
  if (!ConstantElementSize()) {
    FirstElementEntryList1 = FirstElementEntryList();
    ElementSizeList1 = ElementSizeList();
  }
  int MyPID = Comm().MyPID();
  int NumProc = Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);
      if (MyPID==0) {
	os <<  "\nNumber of Global Elements  = "; os << NumGlobalElements(); os << endl;
	os <<    "Number of Global Equations = "; os << NumGlobalEquations(); os << endl;
	os <<    "Maximum of all GIDs        = "; os << MaxAllGID(); os << endl;
	os <<    "Minimum of all GIDs        = "; os << MinAllGID(); os << endl;
	os <<    "Index Base                 = "; os << IndexBase(); os << endl;
	if (ConstantElementSize())
	  os <<  "Constant Element Size      = "; os << ElementSize(); os << endl;
      }
      os << endl;

	os <<    "Number of Local Elements   = "; os << NumMyElements(); os << endl;
	os <<    "Number of Local Equations  = "; os << NumMyEquations(); os << endl;
	os <<    "Maximum of my GIDs         = "; os << MaxMyGID(); os << endl;
	os <<    "Minimum of my GIDs         = "; os << MinMyGID(); os << endl;
      os << endl;

      os.width(14);
      os <<  "     MyPID"; os << "    ";
      os.width(14);
      os <<  "       Local Index "; os << " ";
      os.width(14);
      os <<  "      Global Index "; os << " ";
      if (!ConstantElementSize()) {
	os.width(14);
	os <<" FirstElementEntry "; os << " ";
	os.width(14);
	os <<"   ElementSize "; os << " ";
      }
      os << endl;
    
      for (int i=0; i < NumMyElements(); i++) {
	os.width(14);
	os <<  MyPID; os << "    ";
	os.width(14);
	os <<  i; os << "    ";
	os.width(14);
	os <<  MyGlobalElements1[i]; os << "    ";
	if (!ConstantElementSize()) {	  
	  os.width(14);
	  os << FirstElementEntryList1[i]; os << "    ";
	  os.width(14);
	  os << ElementSizeList1[i]; os << "    ";
	}
	os << endl;
      }
      
      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
  return;
}
