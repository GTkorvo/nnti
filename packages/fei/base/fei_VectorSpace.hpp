/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorSpace_hpp_
#define _fei_VectorSpace_hpp_

#include <fei_macros.hpp>
#include <fei_fwd.hpp>
#include <snl_fei_CommUtils.hpp>
#include <feiArray.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_Logger.hpp>
#include <fei_utils.hpp>
#include <fei_CommUtilsBase.hpp>
#include <fei_ctg_set.hpp>
#include <snl_fei_RaggedTable.hpp>

namespace fei {
  class FieldMask;
  class Lookup_Impl;
  class ParameterSet;
  class Pattern;
  class Record;
  class Record_Operator;
  class SharedIDs;

  /** Class containing the methods for defining a solution-space (a set of
      degrees-of-freedom) and mapping that space to a globally unique set of
      indices. The set of indices can then be used in defining vectors or in
      defining either the row-space or the column-space of matrices.

      Example: define a displacement field over a set of node-identifiers,
      and map that set to a set of equation-numbers.

      There are multiple ways to use an instance of this interface:

      For generating vectors:
      <ol>
      <li>Define fields and identifier-types
      <li>Initialize active fields over sets of identifiers
      <li>Obtain index offset and range information via one of the methods
      'getGlobalIndexOffsets()', 'getIndices_Owned()',
      'getIndices_SharedAndOwned()', etc., to use in constructing or
      initializing a vector object.
      </ol>

      For generating matrices:
      <ol>
      <li>Define fields and identifier-types
      <li>Construct an instance of fei::MatrixGraph (using this
      VectorSpace as a contructor argument or intialization argument) and
      proceed to initialize connectivities and other structural attributes on
      the MatrixGraph object.
      <li>Obtain matrix-graph information from the fei::MatrixGraph
      object to use in constructing or initializing a matrix object.
      </ol>
  */
  class VectorSpace : private fei::Logger {
  public:
    /** VectorSpace Factory interface */
    class Factory {
    public:
      /** Destructor */
      virtual ~Factory(){}

     /** Produce an instance of a VectorSpace. name may be NULL. */
     virtual fei::SharedPtr<VectorSpace> createVectorSpace(MPI_Comm,
							   const char* name);
    };

    /** Constructor.
	@param comm MPI_Communicator
	@param name String to be used in the name of a debug-log file,
	if any. This is an optional argument, defaults to NULL.
    */
    VectorSpace(MPI_Comm comm, const char* name = NULL);

    /** Destructor */
    virtual ~VectorSpace();

    //@{ \name Setup/initialization

    /** Set parameter values from a fei::ParameterSet object. */
    void setParameters(const fei::ParameterSet& paramset);

    /** Define fields that will occur in this solution space. <br>
	Example: a temperature field might be defined as fieldID 0, size 1.<br>
	Example: a velocity field might be defined as fieldID 5, size 3.<br>

	@param numFields Input. Length of the fieldIDs and fieldSizes lists.
	@param fieldIDs Input. List of user-supplied field-identifiers.
	Convention: Active solution-space fields should generally be denoted
	by non-negative field-identifiers, while "other" fields (such as
	geometric coordinates) should be denoted by negative field-identifiers.
	@param fieldSizes Input. List of user-specified field-sizes. A
	field-size is the number of scalar components that make up a field.
    */
    void defineFields(int numFields,
		      const int* fieldIDs,
		      const int* fieldSizes);

    /** Define identifier-types in this solution space. <br>
	For example, define node-identifiers to be type 0, edge-identifiers to
	be type 1, lagrange-multiplier identifiers to be type 2, etc.<br>
	identifier-types need not be zero-based or contiguous.

	@param numIDTypes Number of distinct identifier-types
	@param idTypes User-supplied list of identifier-types
	@return error-code 0 if successful
    */
    void defineIDTypes(int numIDTypes,
		       const int* idTypes);

    /** Add a set of identifiers to this solution-space. These solution-space
	entries consist of fields residing at identifiers.<br>
	Example: temperature field at a set of finite-element nodes.

	@param fieldID Input. The field-identifier to be added. Must be one of
	the fieldIDs defined previously via 'defineFields()'.
	@param numInstancesOfThisFieldPerID Input. It is possible to have
	multiple fields of the same fieldID at each mesh location. e.g., you
	could have 2 pressure fields at each edge.
	@param idType Input. The identifier-type over which the active field is
	being initialized. Must be one of the idTypes defined previously via
	'defineIDTypes()'.
	@param numIDs Input. Length of the IDs list.
	@param IDs Input List of identifiers over which 'fieldID' is active.
	@return error-code 0 if successful
     */
    int initSolutionEntries(int fieldID,
			    int numInstancesOfThisFieldPerID,
			    int idType,
			    int numIDs,
			    const int* IDs);

    /** Add a set of identifiers to the solution-space. These solution-space
	entries consist of identifiers that don't have associated fields.<br>
	Example: Lagrange-multiplier constraint identifiers.<br>
	This method may also be used for initializing a finite-element
	solution-space where the user knows that the entire problem contains only
	one scalar field (e.g., temperature) and so it is sufficient to define
	a solution space on identifiers without associating fields with those
	identifiers. (This will achieve a performance gain for the structure-
	definition, graph-generation and matrix/vector assembly.)

	@param idType Input. The identifier-type over which the solution-
	space is being defined. Must be one of the idTypes defined previously via
	'defineIDTypes()'.
	@param numIDs Input. Number of identifiers being added to the solution-
	space.
	@param IDs Input. List of length numIDs. Identifiers being added to the
	solution-space.
	@return error-code 0 if successful
    */
    int initSolutionEntries(int idType,
			    int numIDs,
			    const int* IDs);

    /** Identify a set of identifiers as being shared with other processors.
	The shared ids must be identified in a globally symmetric way. i.e., if
	the local processor identifies id x as being shared with processor p,
	then processor p MUST identify id x as being shared with the local
	processor.

	@param numShared Input. Length of the lists sharedIDs and 
	                        numSharingProcsPerID.
	@param idType Input. The identifier-type of the ids that are being 
	                        identified as shared.
        @param sharedIDs Input. List of shared identifiers.
	@param numSharingProcsPerID Input. List of length numShared, and the
	  i-th entry gives the number of processors being identified as sharing
	  the i-th sharedID.
	@param sharingProcs Input.
	Packed list of length sum(numSharingProcsPerID), containing the sharing
	processor ranks.
	@return error-code 0 if successful
    */
    int initSharedIDs(int numShared,
		      int idType,
		      const int* sharedIDs,
		      const int* numSharingProcsPerID,
		      const int* sharingProcs);

    /** Identify a set of identifiers as being shared with other processors.
	The shared ids must be identified in a globally symmetric way. i.e., if
	the local processor identifies id x as being shared with processor p,
	then processor p MUST identify id x as being shared with the local
	processor.

	@param numShared Input. Length of the lists sharedIDs and 
	                        numSharingProcsPerID.
	@param idType Input. The identifier-type of the ids that are being 
	                        identified as shared.
        @param sharedIDs Input. List of shared identifiers.
	@param numSharingProcsPerID Input. List of length numShared, and the
	  i-th entry gives the number of processors being identified as sharing
	  the i-th sharedID.
	@param sharingProcs Input.
	Table with 'numShared' rows, and each row is of length
	numSharingProcsPerID. This table contains the sharing processor ranks.
	@return error-code 0 if successful
    */
    int initSharedIDs(int numShared,
		      int idType,
		      const int* sharedIDs,
		      const int* numSharingProcsPerID,
		      const int* const* sharingProcs);

    /** Add the contents of another VectorSpace object to this one.
     */
    int addVectorSpace(fei::VectorSpace* inputSpace);

    /** Indicate that initialization is complete. This is a collective function,
	must be called on all processors. At this time ownership of shared IDs
	will be assigned, and the global index space calculated.

	@return error-code 0 if successful
    */
    int initComplete();
    //@}

    //@{ \name Attribute query methods

    /** Return the CommUtils object held by this vector-space. When built/run in
	serial mode, the MPI_Comm returned by CommUtils::getCommunicator() is
	#defined to be int.
    */
    fei::SharedPtr<snl_fei::CommUtils<int> > getCommUtils();

    /** Given a particular degree-of-freedom, request the corresponding global
	index. A particular degree-of-freedom is specified as a component of a
	particular field, residing at a particular location (ID).

	@param idType Input. Identifier-type of the location at which the 
	specified degree-of-freedom resides. Must be one of the identifier-types
	previously defined via a call to 'defineIDTypes()'.

	@param ID Input. Identifier for the location being specified, such as a
	node-identifier, etc.

	@param fieldID Input. Identifier for the field being specified.

	@param fieldOffset Input. In case there is more than one field with the
	specified fieldID residing at the specified ID, this provides an offset
	into those fields. If only one field with specified fieldID, then this
	parameter is 0.

	@param whichComponentOfField Input. Specifies a scalar component within
	the field. If the field only has 1 scalar component, then this parameter
	is 0.

	@param globalIndex Output. This is the global index of the specified
	degree-of-freedom. Not referenced if the specified degree-of-freedom is
	not found.

	@return error-code 0 if successful. If the specified degree-of-freedom
	is not found, -1 is returned.
     */
    int getGlobalIndex(int idType,
		       int ID,
		       int fieldID,
		       int fieldOffset,
		       int whichComponentOfField,
		       int& globalIndex);

    /** Given a particular degree-of-freedom, request the corresponding global
	index. A particular degree-of-freedom is specified as a component of a
	particular field, residing at a particular location (ID).

	@param idType Input. Identifier-type of the location at which the 
	specified degree-of-freedom resides. Must be one of the identifier-types
	previously defined via a call to 'defineIDTypes()'.

	@param ID Input. Identifier for the location being specified, such as a
	node-identifier, etc.

	@param fieldID Input. Identifier for the field being specified.

	@param globalIndex Output. This is the global index of the specified
	degree-of-freedom. Not referenced if the specified degree-of-freedom is
	not found.

	@return error-code 0 if successful. If the specified degree-of-freedom
	is not found, -1 is returned.
     */
    int getGlobalIndex(int idType,
		       int ID,
		       int fieldID,
		       int& globalIndex);

    /** Given a particular identifier, request the corresponding global
	block-index.

	@param idType Input. Identifier-type of the identifier being queried.

	@param ID Input. Identifier for which a block-index is being requested.

	@param globalBlkIndex Output. This is the global block-index of the
	specified identifier.

	@return error-code 0 if successful. If the specified degree-of-freedom
	is not found, -1 is returned.
     */
    int getGlobalBlkIndex(int idType,
			  int ID,
			  int& globalBlkIndex);

    /** Given a list of IDs, fill an output-list of the global-indices that
	correspond to the first instance of the specified field at each ID.

	@param numIDs Input. Length of the IDs.
	@param IDs Input. User-provided list of identifiers.
	@param idType Input. Type of the IDs for which indices are being
	requested.
	@param fieldID Input. Specified field
	@param globalIndices Output. User-allocated list which, on exit, will
	contain the requested indices. Note that the length of this list is
	assumed to be numIDs*getFieldSize(fieldID).

	@return error-code 0 if successful Note that for any IDs that are not
	found, or IDs which don't have the specified field, the corresponding
	global-index will be -1.
    */
    int getGlobalIndices(int numIDs,
			 const int* IDs,
			 int idType,
			 int fieldID,
			 int* globalIndices);

    /** Given a list of IDs, fill an output-list of the global-block-indices
	that correspond to each	ID.

	@param numIDs Input. Length of the IDs list and of the globalBlkIndices
	list.
	@param IDs Input. User-provided list of identifiers.
	@param idType Input. Type of the IDs for which block-indices are being
	requested.
	@param globalBlkIndices Output. User-allocated list which, on exit,
	will contain the requested indices. Note that the length of this list
	is assumed to be numIDs.

	@return error-code 0 if successful Note that for any IDs that are not
	found, the corresponding global-index will be -1.
    */
    int getGlobalBlkIndices(int numIDs,
			 const int* IDs,
			 int idType,
			 int* globalBlkIndices);

    /** Given a list of IDs, fill an output-list of the global-indices that
	correspond to the first instance of the specified field at each ID.
	Somewhat more general version of the getGlobalIndices() method above.

	@param numIDs Input. Length of the IDs list.
	@param IDs Input. User-provided list of identifiers.
	@param idTypes Input. List of length numIDs, specifying the types of
	the IDs for which indices are being requested.
	@param fieldIDs Input. List of length numIDs, specifying a field at
	each ID.
	@param globalIndices Output. User-allocated list which, on exit, will
	contain the requested indices. Note that the length of this list is
	assumed to be numIDs*getFieldSize(fieldID).

	@return error-code 0 if successful Note that for any IDs that are not
	found, or IDs which don't have the specified field, the corresponding
	global-index will be -1.
    */
    int getGlobalIndices(int numIDs,
			 const int* IDs,
			 const int* idTypes,
			 const int* fieldIDs,			 
			 int* globalIndices);

    /** Given a particular degree-of-freedom, request the corresponding global
	index. In this case, the degree-of-freedom is specified simply by an
	identifier and identifier-type, without specifying a field. This is
	intended to be used for requesting global indices for constraint-
	identifiers or other identifiers which don't have associated fields.
	If the specified identifier actually does have associated fields, then
	the output globalIndex will be the global-index corresponding to the
	first component of the first associated field.

	@param idType Input. Identifier-type of the location at which the 
	specified degree-of-freedom resides. Must be one of the identifier-types
	previously defined via a call to 'defineIDTypes()'.

	@param ID Input. Identifier for the location being specified, such as a
	node-identifier, etc.

	@param globalIndex Output. This is the global index of the specified
	degree-of-freedom. Not referenced if the specified degree-of-freedom is
	not found.

	@return error-code 0 if successful. If the specified degree-of-freedom
	is not found, -1 is returned.
     */
    int getGlobalIndex(int idType,
		       int ID,
		       int& globalIndex);

    /** Given a particular identifier, request the number of scalar degrees-of-
	freedom that are associated with that identifier.
    */
    int getNumDegreesOfFreedom(int idType,
			       int ID);

    /** Query the number of fields defined for this vector-space.
     */
    int getNumFields();

    /** Fill a std::vector with fieldIDs defined for this vector-space.

	@param fieldIDs On exit, contains fieldIDs.
     */
    void getFields(std::vector<int>& fieldIDs);

    /** Given a particular identifier, request the number of fields that are
	associated with that identifier.
    */
    int getNumFields(int idType,
		     int ID);

    /** Given a particular identifier, request the list of fields that are
	associated with that identifier.

	@param idType Identifier-type
	@param ID Specified identifier
	@param lenFieldIDs Input. Length of user-allocated 'fieldIDs' list.
	@param fieldIDs Input. User-allocated list, length must be at least 
	as large as the value produced by getNumFields() for this ID.
	@param numFields Output. Number of fields. If numFields > lenFieldIDs,
	then fieldIDs will contain the first 'lenFieldIDs' field identifiers.
    */
    void getFields(int idType, int ID, std::vector<int>& fieldIDs);

    /** Query for the number of identifier-types defined for this vector-space.
     */
    int getNumIDTypes();

    /** Query for the list of identifier-types defined for this vector-space.

	@param len Input, length of the user-allocated list 'idTypes'.
	@param idTypes Input/Output, user-allocated list, on exit contents will
	contain id-types that are defined for this vector-space.
	@param numIDTypes Output, number of id-types that are defined for this
	vector-space. If numIDTypes is less than user-provided 'len', then only
	'numIDTypes' positions in 'idTypes' are referenced. If numIDTypes is
	greater than user-provided len, then 'idTypes' is filled with the first
	'len' id-types that are defined for this vector-space.
     */
    int getIDTypes(int len, int* idTypes, int& numIDTypes);

    /** Request the number of partitions. (For MPI implementations, partitions
	is a synonym for processes.) The main purpose of this function is to
	give the user a way to calculate the length of the list that needs to
	be allocated before calling 'getGlobalIndexOffsets()'.

	@return number of partitions
    */
    int getNumPartitions();

    /** Request the global index offsets. Indices are zero-based.

	@param lenGlobalOffsets Input. This value gives the length of the
	user-allocated array globalOffsets. Should be numPartitions+1.
	@param globalOffsets Output. Caller-allocated array of length
	lenGlobalOffsets. On exit, contains global-offsets.<br>
	globalOffsets[i] is first global offset on processor i,
	for i in 0 .. numPartitions - 1<br>
	globalOffsets[i+1] - globalOffsets[i] is the number of indices on the
	i-th processor
	@return error-code 0 if successful
    */
    int getGlobalIndexOffsets(int lenGlobalOffsets,
			      int* globalOffsets);

    /** Request the global block-index offsets. Indices are zero-based.

	@param lenGlobalBlkOffsets Input. This value gives the length of the
	user-allocated array globalBlkOffsets. Should be numPartitions+1.
	@param globalBlkOffsets Output. Caller-allocated array of length
	lenGlobalBlkOffsets. On exit, contains global-block-offsets.<br>
	globalBlkOffsets[i] is first global block-offset on processor i,
	for i in 0 .. numPartitions - 1<br>
	globalBlkOffsets[i+1] - globalBlkOffsets[i] is the number of
	block-indices on the i-th processor
	@return error-code 0 if successful
    */
    int getGlobalBlkIndexOffsets(int lenGlobalBlkOffsets,
				 int* globalBlkOffsets);

    /** Given a global index in the point-equation space, return the
	owning processor. If the global index is not in the equation space,
	return -1.
    */
    int getOwnerProcPtIndex(int globalIndex);

    /** Given a global index in the block-equation space, return the
	owning processor. If the global index is not in the equation space,
	return -1.
    */
    int getOwnerProcBlkIndex(int globalIndex);

    /** Given an identifier (with identifier-type), return true if it resides
	on the local processor, false if not. This will return true if the 
	identifier is either owned or shared by the local processor.
    */
    bool isLocal(int idType, int ID);

    /** Given an identifier (with identifier-type), return true if it resides
	on the local processor and is locally owned, false if not.
    */
    bool isLocallyOwned(int idType, int ID);

    /** Request the field-size for a specified field-identifier. If the specified
	field-identifier is not found, std::runtime_error is thrown.
	@param fieldID Input. Specified field-identifier
    */
    unsigned getFieldSize(int fieldID);

    /** Query the number of locally owned-or-shared identifiers. */
    int getNumOwnedAndSharedIDs(int idType);

    /** Query the number of locally-owned identifiers. */
    int getNumOwnedIDs(int idType);

    /** Obtain a list of the local identifiers. Note that this includes
     identifiers that are locally shared but not owned. */
    int getOwnedAndSharedIDs(int idtype,
		    int lenList,
		    int* IDs,
		    int& numLocalIDs);

    /** Obtain a list of the locally owned identifiers.
     */
    int getOwnedIDs(int idtype,
			   int lenList,
			   int* IDs,
			   int& numLocalIDs);

    /** Query number of indices on local processor, including ones that are
	locally owned as well as shared-but-not-owned. Only available after
	initComplete has been called. (returns 0 before that)
    */
    int getNumIndices_SharedAndOwned() const;

    /** Obtain list of global indices on local processor, including ones that
        are locally owned as well as shared-but-not-owned. Only available
        after initComplete has been called.

	@param lenIndices Input. Length of user-allocated 'globalIndices' list.
	@param globalIndices User-allocated list. On output, will contain all
	indices owned or shared by local processor.
	@param numIndices Output. Number of indices. If 'numIndices' is
            different than 'lenIndices', then globalIndices will contain
	'min(lenIndices, numIndices)' of the local processor's indices.
    */
    int getIndices_SharedAndOwned(int lenIndices,
				  int* globalIndices,
				  int& numIndices) const;

    /** Query number of block indices on local processor, including ones that
        are locally owned as well as shared-but-not-owned. Only available after
	initComplete has been called.
    */
    int getNumBlkIndices_SharedAndOwned(int& numBlkIndices) const;

    /** Obtain list of global block indices on local processor, including ones
	that are locally owned as well as shared-but-not-owned. Only available
        after initComplete has been called.

	@param lenBlkIndices Input. Length of user-allocated 'globalBlkIndices'
	list.
	@param globalBlkIndices User-allocated list. On output, will contain all
	indices owned or shared by local processor.
	@param blkSizes User-allocated list. On output, will contain the number
         of scalars (point-indices) associated with each corresponding
         block-index.
	@param numBlkIndices Output. Number of indices. If 'numBlkIndices' is
	different than 'lenBlkIndices', then globalBlkIndices will contain
	'min(lenBlkIndices, numBlkIndices)' of the local processor's indices.
    */
    int getBlkIndices_SharedAndOwned(int lenBlkIndices,
				     int* globalBlkIndices,
				     int* blkSizes,
				     int& numBlkIndices);

    /** Query number of indices owned by local processor.
    */
    int getNumIndices_Owned() const;

    /** Obtain list of global indices owned by local processor. Only
	available after initComplete has been called.

	@param lenIndices Input. Length of user-allocated 'globalIndices' list.
	@param globalIndices User-allocated list. On output, will contain all
	indices owned by local processor.
	@param numIndices Output. Number of indices. If 'numIndices' is different
	than 'lenIndices', then globalIndices will contain
	'min(lenIndices, numIndices)' of the local processor's indices.
    */
    int getIndices_Owned(int lenIndices,
			 int* globalIndices,
			 int& numIndices) const;

    /** Query number of block indices owned by local processor.
    */
    int getNumBlkIndices_Owned() const;

    /** Obtain list of global block indices owned by local processor. Only
	available after	initComplete has been called.

	@param lenBlkIndices Input. Length of user-allocated 'globalBlkIndices'
	list.
	@param globalBlkIndices User-allocated list. On output, will contain all
	indices owned by local processor.
	@param blkSizes User-allocated list. On output, will contain the number of
	scalars (point-indices) associated with each corresponding block-index.
	@param numBlkIndices Output. Number of indices. If 'numBlkIndices' is
	different than 'lenBlkIndices', then globalBlkIndices will contain
	'min(lenBlkIndices, numBlkIndices)' of the local processor's indices.
    */
    int getBlkIndices_Owned(int lenBlkIndices,
			    int* globalBlkIndices,
			    int* blkSizes,
			    int& numBlkIndices);

    /** Query the number of shared identifiers of a given id-type. */
    int getNumSharedIDs(int idType, int& numShared);

    /** Query the number of eqn indices across all processors.
     */
    int getGlobalNumIndices() const;

    /** Query the number of block-eqn indices across all processors.
     */
    int getGlobalNumBlkIndices() const;

    /** Intended to be used by other snl_fei:: classes.
    */
    int getRecordCollection(int idType, snl_fei::RecordCollection*& records);

    /** Intended to be used only by other snl_fei:: classes.
    */
    std::vector<int>& getEqnNumbers();

    /** Intended to be used by other implementation classes.
    */
    snl_fei::PointBlockMap* getPointBlockMap();

    void getGlobalIndices(const fei::Pattern* pattern,
			  const fei::Record*const* records,
			  std::vector<int>& indices);

    void getGlobalBlkIndices(const fei::Pattern* pattern,
			     const fei::Record*const* records,
			     std::vector<int>& indices);

    void getGlobalIndices(int numRecords,
			  const fei::Record*const* records,
			  int fieldID,
			  int fieldSize,
			  int indicesAllocLen,
			  int* indices,
			  int& numIndices);

    void getGlobalIndices(int numRecords,
			  const fei::Record*const* records,
			  const int* numFieldsPerID,
			  const int* fieldIDs,
			  const int* fieldSizes,
			  int indicesAllocLen,
			  int* indices,
			  int& numIndices);

    void getGlobalBlkIndices(int numRecords,
			     const fei::Record*const* records,
			     int indicesAllocLen,
			     int* indices,
			     int& numIndices);

    int initSolutionEntries(int fieldID,
			    int numInstancesOfThisFieldPerID,
			    int idType,
			    int numIDs,
			    const int* IDs,
			    fei::Record** records);

    int initSolutionEntries(int idType,
			    int numIDs,
			    const int* IDs,
			    fei::Record** records);

    std::vector<fei::FieldMask*> fieldMasks_;

  private:
    friend class fei::Lookup_Impl;

  private:
    VectorSpace(const VectorSpace& src);
    VectorSpace& operator=(const VectorSpace& src);

    inline void check_version() { fei::utils::version(); }

    int setOwners_lowestSharing();

    int calculateGlobalIndices();

    int runRecords(fei::Record_Operator& record_op);

    int synchronizeSharedRecords();

    int setLocalEqnNumbers();

    int exchangeGlobalIndices();

    int exchangeFieldInfo(fei::comm_map* ownerPattern,
			  fei::comm_map* sharerPattern,
			  snl_fei::RecordCollection* recordCollection,
			  std::vector<fei::FieldMask*>& fieldMasks);

    int getSharedIDs_private(int idType, fei::SharedIDs*& shIDs);

    void setName(const char* name);

  private:
    fei::SharedPtr<snl_fei::CommUtils<int> > intCommUtils_;

    feiArray<int> idTypes_;
    std::map<int,unsigned> fieldDatabase_;
    int maxFieldSize_;
    feiArray<snl_fei::RecordCollection*> recordCollections_;

    feiArray<int> sharedIDTypes_;
    feiArray<fei::SharedIDs*> sharedIDTables_;
    feiArray<fei::comm_map*> ownerPatterns_;
    feiArray<fei::comm_map*> sharerPatterns_;

    bool sharedRecordsSynchronized_;

    snl_fei::PointBlockMap* ptBlkMap_;

    std::vector<int> globalOffsets_;
    std::vector<int> globalIDOffsets_;

    bool simpleProblem_;

    int firstLocalOffset_, lastLocalOffset_;

    std::vector<int> eqnNumbers_;

    bool newInitData_;

    std::string name_;
    std::string dbgprefix_;
    bool checkSharedIDs_;
  }; // class fei::VectorSpace

  inline fei::SharedPtr<snl_fei::CommUtils<int> > VectorSpace::getCommUtils()
    { return( intCommUtils_ ); }

  inline std::vector<int>& VectorSpace::getEqnNumbers()
    {
      return( eqnNumbers_ );
    }

  inline snl_fei::PointBlockMap* VectorSpace::getPointBlockMap()
    {
      return( ptBlkMap_ );
    }

} // namespace fei

#endif // _fei_VectorSpace_hpp_
