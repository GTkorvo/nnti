/*Paul
07-Oct-2002 ElementSpaceData move started
13-Oct-2002 Rewritten with massive constructor call
22-Oct-2002 Modified slightly - ESData constructor now takes Comm* argument
*/

#ifndef _TPETRA_ELEMENTSPACEDATA_H_
#define _TPETRA_ELEMENTSPACEDATA_H_

namespace Tpetra {

template<typename OrdinalType>
class ElementSpaceData : public Object {
	friend class ElementSpace<OrdinalType>;
 public:
	ElementSpaceData(const OrdinalType indexBase, 
									 const OrdinalType numGlobalElements,
									 const OrdinalType numMyElements,
									 const OrdinalType minAllGID,
									 const OrdinalType maxAllGID,
									 const OrdinalType minMyGID,
									 const OrdinalType maxMyGID,
									 const map<OrdinalType, OrdinalType> lgMap,
									 const map<OrdinalType, OrdinalType> glMap,
									 const bool contiguous,
									 const Platform<OrdinalType, OrdinalType>& Platform,
									 const Comm<OrdinalType, OrdinalType>* Comm)
		: Object("Tpetra::ElementSpaceData")
		, indexBase_(indexBase)
		, numGlobalElements_(numGlobalElements)
		, numMyElements_(numMyElements)
		, minMyGID_(minMyGID)
		, maxMyGID_(maxMyGID)
		, minAllGID_(minAllGID)
		, maxAllGID_(maxAllGID)
		, minLID_(indexBase)
		, maxLID_(indexBase + numMyElements)
		, lgMap_(lgMap)
		, glMap_(glMap)
		, contiguous_(contiguous)
		, global_(checkGlobalness())
		, myGlobalElements_(0)
		, Platform_(&Platform) 
		, Comm_(Comm) 
		, Directory_(0) 
		{};

	~ElementSpaceData() {
		if(Comm_ != 0) {
			delete Comm_;
			Comm_ = 0;
		}
		if(myGlobalElements_ != 0) {
			delete [] myGlobalElements_;
			myGlobalElements_ = 0;
		}
		if(Directory_ != 0) {
			delete Directory_;
			Directory_ = 0;
		}
	};

 protected:
	const Platform<OrdinalType, OrdinalType>* Platform_;
	const Comm<OrdinalType, OrdinalType>* Comm_;
	const OrdinalType numGlobalElements_;
	const OrdinalType numMyElements_;
	const OrdinalType indexBase_;
	const OrdinalType minLID_;
	const OrdinalType maxLID_;
	const OrdinalType minMyGID_;
	const OrdinalType maxMyGID_;
	const OrdinalType minAllGID_;
	const OrdinalType maxAllGID_;
	const bool contiguous_;
	const bool global_;
  map<OrdinalType, OrdinalType> lgMap_;
  const map<OrdinalType, OrdinalType> glMap_;
	OrdinalType* myGlobalElements_;
	Directory<OrdinalType>* Directory_;

 private:
	bool checkGlobalness() {
		bool global = false;
		if(Platform_->getNumImages() > 1) {
			int localRep = 0;
			int allLocalRep;
			if(numGlobalElements_ == numMyElements_)
				localRep = 1;
			Comm_->minAll(&localRep, &allLocalRep, 1);
			if(allLocalRep != 1)
				global = true;
		}
		return(global);
	}

}; // class ElementSpaceData

} // namespace Tpetra

#endif // _TPETRA_ELEMENTSPACEDATA_H_
