#ifndef _TPETRA_BLOCKELEMENTSPACE_HPP_
#define _TPETRA_BLOCKELEMENTSPACE_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include <Teuchos_RefCountPtr.hpp>

namespace Tpetra {

template<typename OrdinalType> class BlockElementSpaceData;

//! Tpetra::BlockElementSpace: A class for constructing and using template<OrdinalType> BlockElementSpaces.
/*! BlockElementSpace objects can have variable element sizes. (If variable element sizes are not needed, an ElementSpace object should probably be used instead.) Some BlockElementSpace methods throw exceptions, and should be enclosed in a try/catch block. All BlockElementSpace objects require an ElementSpace object, which requires a Comm object.  

BlockElementSpace error codes (positive for non-fatal, negative for fatal):
  <ol>
  <li> +1  Specified Point ID not found on this image.
  <li> +2  Specified Local ID not found on this image.
  <li> +3  elementSize requested in a variable-sized BlockElementSpace.
  <li> +4  Pointer passed to getFirstPointInElementList, getElementSizeList, or getPointToElementList does not have child allocated. (Null pointer).
  <li> -1  elementSize (or element in elementSizeList) <= 0.  Should be > 0.
  <li> -99 Internal BlockElementSpace error.  Contact developer.
  </ol>*/

template<typename OrdinalType> 
class BlockElementSpace : public Object {

public:

//@{ \name Constructor/Destructor Methods.

//! Tpetra::BlockElementSpace constructor with constant element size.
BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType elementSize);

//! Tpetra::BlockElementSpace constructor with arbitrary element sizes.
BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType* elementSizeList);

//! Tpetra::BlockElementSpace copy constructor.
BlockElementSpace(BlockElementSpace<OrdinalType> const& BlockElementSpace);

//! Tpetra::BlockElementSpace destructor.
~BlockElementSpace() {};

//@}


//@{ \name Local/Global ID Accessor Methods

//! Returns the image IDs, corresponding local index values, and element sizes for a given list of global indices.
/*! Theimage IDs, local index values, and element sizes are placed into arrays passed in by the user. The list of global indices used to create these is also passed in by the user. Exceptions might be thrown. */ 
void getRemoteIDList(OrdinalType numIDs, OrdinalType* GIDList, OrdinalType* imageIDList, 
										 OrdinalType* LIDList, OrdinalType* elementSizeList) const;

//! Returns the local ID of the element that contains the given local Point ID, and the offset of the point in that element.
/*! The local ID and offset are placed in OrdinalType variables passed in by reference by the user. */
void getLocalElementID(OrdinalType pointID, OrdinalType& elementID, OrdinalType& elementOffset) const;

//@}


//@{ \name Size & Dimension Accessor Methods

//! Returns the size of elements in the BlockElementSpace. Throws an exception of +2 if not all elements are the same size.
OrdinalType getElementSize() const;

//! Returns the size of the element whose local ID is passed in. Throws an exception of +1 if the local ID is not found on the calling image.
OrdinalType getElementSize(OrdinalType LID) const;

//! Returns the number of global points in the BlockElementSpace; equals the sum of all element sizes across all images.
OrdinalType getNumGlobalPoints() const {return(BlockElementSpaceData_->numGlobalPoints_);};

//! Returns the number of global points on this image; equals the sum of all element sizes on the calling image.
OrdinalType getNumMyPoints() const {return(BlockElementSpaceData_->numMyPoints_);};

//! Returns the minimum element size on the calling image.
OrdinalType getMinMyElementSize() const {return(BlockElementSpaceData_->minMySize_);};

//! Returns the maximum element size on the calling image.
OrdinalType getMaxMyElementSize() const {return(BlockElementSpaceData_->maxMySize_);};

//! Returns the minimum element size in the BlockElementSpace.
OrdinalType getMinElementSize() const {return(BlockElementSpaceData_->minGlobalSize_);};

//! Returns the maximum element size in the BlockElementSpace.
OrdinalType getMaxElementSize() const {return(BlockElementSpaceData_->maxGlobalSize_);};

//@}


//@{ \name Misc. Boolean Tests

//! Returns true if all elements have a constant size, returns false otherwise.
bool isConstantElementSize() const {return(BlockElementSpaceData_->constantSize_);};

//! Returns true if this BlockElementSpace is identical to the one passed in, returns false otherwise. Also implemented through the == and != operators.
bool isSameAs(BlockElementSpace<OrdinalType> const& BlockElementSpace) const;
bool operator==(BlockElementSpace<OrdinalType> const& BlockElementSpace) const {return(isSameAs(BlockElementSpace));};
bool operator!=(BlockElementSpace<OrdinalType> const& BlockElementSpace) const {return(!isSameAs(BlockElementSpace));};

//@}


//@{ \name Array Accessor Functions. Each of these methods is implemented twice, one that returns a pointer, and one that copies the array into one passed in by the user.

//! Returns a pointer to array of the sizes of all the elements that belong to the calling image.
OrdinalType const* getElementSizeList() const {return(BlockElementSpaceData_->elementSizeList_);};
void getElementSizeList(OrdinalType* elementSizeList) const;

//! Returns a pointer to the internal array of the mapping between the local elements, and the first local point number in each element.
OrdinalType const* getFirstPointInElementList() const;
void getFirstPointInElementList(OrdinalType* firstPointInElementList) const;

//! Returns a pointer to an array that lists the LID of the element that each point belongs to.
OrdinalType const* getPointToElementList() const;
void getPointToElementList(OrdinalType* pointToElementList) const;

//@}


//@{ Misc.

//! Prints the BlockElementSpace object to the output stream. (Used by the overloaded << operator inherited from Object).
void print(ostream& os) const;

//! Access function for ElementSpace object.
ElementSpace<OrdinalType> const& elementSpace() const {return(*BlockElementSpaceData_->ElementSpace_);};

//@}

private:
Teuchos::RefCountPtr< BlockElementSpaceData<OrdinalType> > BlockElementSpaceData_; // Teuchos smart pointer

}; // BlockElementSpace class

} // Tpetra namespace

#include "Tpetra_BlockElementSpaceData.hpp"


// begin Tpetra_BlockElementSpace.cpp
//=======================================================================

namespace Tpetra {

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType elementSize) 
  : Object("Tpetra::BlockElementSpace") 
	, BlockElementSpaceData_()
{
	// get ES data
  OrdinalType numMyElements = ElementSpace.getNumMyElements();
	OrdinalType numGlobalElements = ElementSpace.getNumGlobalElements();

	// initial throws
  if(elementSize <= 0)
    throw reportError("elementSize = " + toString(elementSize) + ".  Should be > 0.", -1);

	// initialize elementSizeList
	OrdinalType* elementSizeList = 0;
  if(numMyElements > 0) {
		elementSizeList = new OrdinalType[numMyElements];
		for(OrdinalType i = 0; i < numMyElements; i++)
			elementSizeList[i] = elementSize;
	}

	// set sizes
  OrdinalType minMySize = elementSize;
  OrdinalType maxMySize = elementSize;
  OrdinalType minGlobalSize = elementSize;
  OrdinalType maxGlobalSize = elementSize;

	// compute numGlobalPoints & numMyPoints
  OrdinalType numGlobalPoints = elementSize * numGlobalElements;
  OrdinalType numMyPoints = elementSize * numMyElements;

	// call BESData constructor
	BlockElementSpaceData_ = Teuchos::rcp(new BlockElementSpaceData<OrdinalType>(ElementSpace, true, elementSize, numMyPoints, 
																																							 numGlobalPoints, minMySize, maxMySize, 
																																							 minGlobalSize, maxGlobalSize, elementSizeList));
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType* elementSizeList) 
  : Object("Tpetra::BlockElementSpace")
	, BlockElementSpaceData_()
{
	// get ES data
  OrdinalType numMyElements = ElementSpace.getNumMyElements();

	// initial throws
  for(OrdinalType i = 0; i < numMyElements; i++)
    if(elementSizeList[i] <= 0)
      throw reportError("An element in elementSizeList = " + toString(elementSizeList[i]) + ".  Should be > 0.", -1);

	// initialize elementSizeList and compute minMySize, MaxMySize, & numMyPoints
	//   we copy elementSizeList into our own array because elementSizeList (the user's array)
	//   is not guaranteed to always hold the same values it does now.
	OrdinalType* myElementSizeList = 0;
	OrdinalType minMySize = 1;
	OrdinalType maxMySize = 1;
	OrdinalType numMyPoints = 0;
	if(numMyElements > 0) {
    myElementSizeList = new OrdinalType[numMyElements];
    minMySize = elementSizeList[0];
    maxMySize = elementSizeList[0];
    numMyPoints = 0;
    for(OrdinalType i = 0; i < numMyElements; i++) {
      myElementSizeList[i] = elementSizeList[i];
      minMySize = TPETRA_MIN(minMySize, elementSizeList[i]);
      maxMySize = TPETRA_MAX(maxMySize, elementSizeList[i]);
      numMyPoints += elementSizeList[i];
		}
	}

	// compute minGlobalSize & maxGlobalSize
	OrdinalType numGlobalPoints = numMyPoints;
	OrdinalType minGlobalSize = minMySize;
	OrdinalType maxGlobalSize = maxMySize;
  if(ElementSpace.isGlobal() == true) {
    ElementSpace.comm().sumAll(&numMyPoints, &numGlobalPoints, 1);
    ElementSpace.comm().minAll(&minMySize, &minGlobalSize, 1);
    ElementSpace.comm().maxAll(&maxMySize, &maxGlobalSize, 1);
  }

	// call BESData constructor
	BlockElementSpaceData_ = Teuchos::rcp(new BlockElementSpaceData<OrdinalType>(ElementSpace, false, 0, numMyPoints, 
																																							 numGlobalPoints, minMySize, maxMySize, 
																																							 minGlobalSize, maxGlobalSize, myElementSizeList));
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(BlockElementSpace<OrdinalType> const& BlockElementSpace) 
  : Object(BlockElementSpace.label())
	, BlockElementSpaceData_(BlockElementSpace.BlockElementSpaceData_) {}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getRemoteIDList(OrdinalType numIDs, OrdinalType* GIDList, OrdinalType* imageIDList, 
																										 OrdinalType* LIDList, OrdinalType* elementSizeList) const {
  elementSpace().getRemoteIDList(numIDs, GIDList, imageIDList, LIDList);

  if(isConstantElementSize())
    for(OrdinalType i = 0; i < numIDs; i++)
      elementSizeList[i] = getElementSize();
  else
    for(OrdinalType i = 0; i < numIDs; i++)
      elementSizeList[i] = getElementSize(LIDList[i]);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getLocalElementID(OrdinalType pointID, OrdinalType& elementID, OrdinalType& elementOffset) const {
  pointID -= elementSpace().getIndexBase(); // convert from indexBase-based to zero-based counting.
  if(pointID < 0 || pointID > getNumMyPoints())
    throw reportError("PointID " + toString(pointID) + " was not found on this processor.", 1);
  if(isConstantSize()) {
    elementID = pointID / getElementSize();
    elementOffset = pointID % getElementSize();
  }
  else {
    OrdinalType* tmp = getPointToElementList();
    elementID = tmp[pointID];
    tmp = getFirstPointInElementList();
    elementOffset = pointID - tmp[elementID];
    tmp = 0;
  }
}

//=======================================================================
template<typename OrdinalType>
OrdinalType BlockElementSpace<OrdinalType>::getElementSize() const {
  if(!isConstantElementSize()) 
    throw reportError("This BlockElementSpace does not have a constant element size.", 3);
  return(BlockElementSpaceData_->elementSize_);
}

//=======================================================================
template<typename OrdinalType>
OrdinalType BlockElementSpace<OrdinalType>::getElementSize(OrdinalType LID) const {
  if(elementSpace().isMyLID(LID) == false)
    throw reportError("Local ID " + toString(LID) + " was not found on this processor.", 2);
  if(isConstantSize())
    return(getElementSize());
  else {
    LID -= elementSpace().getIndexBase(); // convert to zero-based counting.
    return(BlockElementSpaceData_->elementSizeList_[LID]);
  }
}

//=======================================================================
template<typename OrdinalType>
bool BlockElementSpace<OrdinalType>::isSameAs(BlockElementSpace<OrdinalType> const& BlockElementSpace) const {
  if(this == &BlockElementSpace)
    return(true);
  // check to make sure ElementSpaces match.
  if(! elementSpace().isSameAs(BlockElementSpace.elementSpace()))
    return(false);
  if(isConstantElementSize() != BlockElementSpace.isConstantElementSize())
    return(false);

  if(isConstantElementSize()) {
    if(getElementSize() != BlockElementSpace.getElementSize())
      return(false);
    else
      return(true);
  }
  else {
    int mySameBES = 1;
    OrdinalType nME = elementSpace().getNumMyElements();
    for(OrdinalType i = 0; i < nME; i++)
      if(BlockElementSpaceData_->elementSizeList_[i] != BlockElementSpace.BlockElementSpaceData_->elementSizeList_[i])
				mySameBES = 0;
    
    int globalSameBES = 0;
    elementSpace().comm().minAll(&mySameBES, &globalSameBES, 1);
    return(globalSameBES == 1);
  }
}

//=======================================================================
// LID -> size of that element
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getElementSizeList(OrdinalType* elementSizeList) const {
	if(elementSizeList == 0) 
    throw reportError("This pointer does not have a child allocated.", 4);
  OrdinalType nME = elementSpace().getNumMyElements();
  for(OrdinalType i = 0; i < nME; i++)
    elementSizeList[i] = BlockElementSpaceData_->elementSizeList_[i];
}

//=======================================================================
// LID -> lowest PointID contained in that element
template<typename OrdinalType>
OrdinalType const* BlockElementSpace<OrdinalType>::getFirstPointInElementList() const {
	OrdinalType nME = elementSpace().getNumMyElements();
	if((BlockElementSpaceData_->firstPointList_ == 0) && (nME > 0)) {
		OrdinalType* tmpPtr = new OrdinalType[nME];
		getFirstPointInElementList(tmpPtr);
		BlockElementSpaceData_->firstPointList_ = tmpPtr;
	}
	return(BlockElementSpaceData_->firstPointList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getFirstPointInElementList(OrdinalType* firstPointInElementList) const {
	if(firstPointInElementList == 0) 
    throw reportError("This pointer does not have a child allocated.", 4);
  OrdinalType iB = elementSpace().getIndexBase();
	firstPointInElementList[0] = iB;
	OrdinalType nME = elementSpace().getNumMyElements();
	for(OrdinalType i = 1; i < nME; i++)
		firstPointInElementList[i] = firstPointInElementList[i-1] + BlockElementSpaceData_->elementSizeList_[i-1];
}

//=======================================================================
// pointID -> LID containing that point
template<typename OrdinalType>
OrdinalType const* BlockElementSpace<OrdinalType>::getPointToElementList() const {
	OrdinalType numPoints = getNumMyPoints();
	if((BlockElementSpaceData_->pointToElementList_ == 0) && (numPoints > 0)) {
		OrdinalType* tmpPtr = new OrdinalType[numPoints];
		//for(OrdinalType i = 0; i < numPoints; i++)
		//	tmpPtr[i] = i+1;
		getPointToElementList(tmpPtr);
		BlockElementSpaceData_->pointToElementList_ = tmpPtr;
	}
	return(BlockElementSpaceData_->pointToElementList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getPointToElementList(OrdinalType* pointToElementList) const {
	if(pointToElementList == 0) 
    throw reportError("This pointer does not have a child allocated.", 4);
	OrdinalType currPos = 0;
	OrdinalType nME = elementSpace().getNumMyElements();
	OrdinalType currLID = elementSpace().getIndexBase();
	OrdinalType currSize;
	for(OrdinalType i = 0; i < nME; i++) {
		currSize = BlockElementSpaceData_->elementSizeList_[i];
		for(OrdinalType j = 0; j < currSize; j++) {
			pointToElementList[currPos] = currLID;
			currPos++;
		}
		currLID++;
	}
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::print(ostream& os) const {
	OrdinalType const* elementSizeList1 = getElementSizeList();
  OrdinalType const* firstPointList1 = getFirstPointInElementList(); 
  OrdinalType const* pointToElementList1 = getPointToElementList(); 
 
  int myImageID = elementSpace().comm().getMyImageID();
  int numImages = elementSpace().comm().getNumImages();

  for(int imageCtr = 0; imageCtr < numImages; imageCtr++) {
    if(myImageID == imageCtr) {
      if(myImageID == 0) {
				os << "\nNumber of Global Points = "; os << getNumGlobalPoints(); os << endl;
				if(isConstantElementSize()) {
					os << "Constant Element Size   = "; os << getElementSize(); os << endl;
				}
				else {
					os << "Global Min Element Size = "; os << getMinElementSize(); os << endl;
					os << "Global Max Element Size = "; os << getMaxElementSize(); os << endl;
				}
      }
      os << endl;
			
      os <<     "Number of Local Points  = "; os << getNumMyPoints(); os << endl;
      if(!isConstantElementSize()) {
				os <<     "Min Local Element Size  = "; os << getMinMyElementSize(); os << endl;
				os <<     "Max Local Element Size  = "; os << getMaxMyElementSize(); os << endl;
      }
			
      os << "\nElementSizeList"; os << endl;
      OrdinalType nME = elementSpace().getNumMyElements();
      for(OrdinalType i = 0; i < nME; i++)
				os << elementSizeList1[i] << " ";
      os << endl;
      os << "\nFirstPointInElementList"; os << endl;
      for(OrdinalType i = 0; i < nME; i++)
				os << firstPointList1[i] << " ";
      os << endl;
      os << "\nPointToElementList"; os << endl;
      for(OrdinalType i = 0; i < BlockElementSpaceData_->numMyPoints_; i++)
				os << pointToElementList1[i] << " ";
			os << endl;
      os << flush;
    }
	}

  elementSpace().print(os);

  // global ops to let I/O complete done by elementSpace::print.
}

} // Tpetra namespace
//=======================================================================

//=======================================================================
// end Tpetra_BlockElementSpace.cpp

#endif // _TPETRA_BLOCKELEMENTSPACE_HPP_
