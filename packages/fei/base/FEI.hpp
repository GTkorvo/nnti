#ifndef _FEI_hpp_
#define _FEI_hpp_
/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
//=========================================================================
// public C++ interface definition for FEI, version 2.1, Jan. 25, 2001
//
// This interface header is the C++ expression of the Finite Element
// Interface to Linear Solvers.
//
//=========================================================================
// Work has begun, in specifying FEI 3.0. FEI 3.0 will consist of a number
// of separate interfaces, in contrast to FEI 2.x which is a single large
// interface that manages the entire assembly/solution process. Thus the
// FEI 3.0 interfaces will be spread among several files, reflecting the
// class names of the individual interfaces. This header, containing FEI 2.x,
// will continue to exist for users who don't wish to switch to the new
// specification.
//
// Below is a list of the functions in FEI 2.x, grouped roughly by
// usage categories, with the initialization and load groups (1 and 2)
// also reflecting common calling order.
//
//===================================================================
//
//  (0) Construction
//  ----------------
//
//  <Implementation-dependent>
//
//
//  (1) initialization
//  ------------------
//
//      parameters
//      setIDLists
//      initSolveStep
//
//      initFields
//      initElemBlock
//
//      initElem
//
//      initSharedNodes
//
//      initCRMult
//      initCRPen
//      initSlaveVariable
//
//      initCoefAccessPattern
//      initCoefAccess
//
//      initComplete
//
//
//  (2) load data
//  -------------
//
//      setLHSID
//      setRHSID
//
//      loadNodeBCs
//      loadElemBCs
//
//      sumInElem
//      sumInElemMatrix
//      sumInElemRHS
//      loadElemTransfer
//
//      loadCRMult
//      loadCRPen
//
//      sumIntoMatrix
//      putIntoMatrix
//      getFromMatrix
//      sumIntoRHS
//      putIntoRHS
//      getFromRHS
//
//      setMatScalars
//      setRHSScalars
//
//      loadComplete
//
//  (3) equation solution
//  ---------------------
//
//      residualNorm
//      solve
//
//
//  (4) return of solution data
//  ---------------------------
//
//      getBlockNodeSolution
//      getBlockFieldNodeSolution
//      getBlockElemSolution
//      getCRMultipliers
//
//  (5) initial guess 'put'ing
//  --------------------------
//
//      putBlockNodeSolution
//      putBlockFieldNodeSolution
//      putBlockElemSolution
//      putCRMultParam
//      putNodalFieldData
//
//  (6) resetting functions
//  -------------------------
//
//      resetSystem
//      resetMatrix
//      resetRHSVector
//      resetInitialGuess
//      deleteMultCRs
//
//  (7) query functions
//  -------------------
//
//      version
//      cumulative_cpu_times
//
//      iterations
//      getNumSolnParams
//      getNumElemBlocks
//      getNumBlockActNodes
//      getNumBlockActEqns
//      getNumNodesPerElement
//      getNumEqnsPerElement
//      getNumBlockElements
//      getNumBlockElemDOF
//      getNumCRMultipliers
//      getCRMultIDList
//      getBlockNodeIDList
//      getBlockElemIDList
//      getFieldSize
//      getEqnNumbers
//      getNodalFieldSolution
//      getNumLocalNodes
//      getLocalNodeIDList
//
 
#include <fei_macros.hpp>
//Include some #defines
#include <fei_defs.h>

/** public Finite Element Interface specification, version 2.1, in C++.
Note: all FEI functions return an int error code. A value of 0 indicates
that there was no error. Errors are usually indicated by -1, except where
noted in the documentation for a particular function.
*/
class FEI { 
 
  public: 

   /** Destructor. This will call through to the implementation's destructor.
   */
   virtual ~FEI() {} 

    /** Set parameters associated with solver choice, etc. This function may be
       called at any time after the FEI object is instantiated. This function
       may be called repeatedly with different parameters, which will accumulate
       those parameters into an internal 'master'-list of parameters.
       @param numParams Number of parameters being supplied.
       @param paramStrings List of 'numParams' strings. Each string usually
            contains a key-value pair, separated by a space.
    */
   virtual int parameters(int numParams, 
                          const char *const* paramStrings) = 0;

   /**Specify matrixIDs and rhsIDs to be used in cases where multiple matrices
      and/or rhs vectors are being assembled. This function does not need
      to be called if only one matrix and rhs are being assembled. Note: the
      values of the matrix and rhs identifiers must be non-negative. Important
      Note: If this function is called, it must be called BEFORE setSolveType
      is called. setSolveType must then be called with the parameter
      FEI_AGGREGATE_SUM (eigen-solves and product-solves aren't supported
      yet).
      @param numMatrices length of matrixIDs parameter
      @param matrixIDs list of user-defined identifiers for separate matrices
                       to be assembled
      @param numRHSs length of rhsIDs parameter
      @param rhsIDs list of user-defined identifiers for separate rhs vectors
                    to be assembled
   */

   virtual int setIDLists( int numMatrices,
                           const int* matrixIDs,
                           int numRHSs,
                           const int* rhsIDs ) = 0;


    /** Set the type of solve to be performed.
        This distinguishes between a 'standard' single solve of Ax=b,
        an eigen-solve (not yet supported), an 'aggregate-sum' solve
        (a linear-combination of several separate A's and b's), and an
        'aggregate-product' solve (not supported).
        @param solveType currently supported values for this are:
                          FEI_SINGLE_SOLVE, FEI_AGGREGATE_SUM
    */
   virtual int setSolveType( int solveType ) = 0;    

    /** Identify all the fields present in the analysis. A field may be a
        scalar such as temperature or pressure, or a 3-vector for velocity, etc.
        Non-solution fields may be denoted by a negative fieldID. This allows
        for situations where the application wants to pass data that the FEI
        doesn't need, (such as geometric coordinates) through to the
        underlying linear algebra library. This may be done via the various
        put*Solution functions.

        @param numFields Global number of fields in the entire problem, on all
               processors. (This is the length of the fieldSizes and fieldIDs
               lists.)
        @param fieldSizes Number of scalars contained in each field.
        @param fieldIDs User-supplied identifiers for each field.
    */
   virtual int initFields( int numFields,
                           const int *fieldSizes, 
                           const int *fieldIDs ) = 0;
                                
    /** Initialize the description of an element-block. This function informs
        the fei implementation of the defining characteristics for a block of
        elements. An element-block must be homogeneous -- all elements in the
        block must have the same number of nodes, same number of solution fields
        per node, etc.
        @param elemBlockID The user-defined identifier for this element-block.
        @param numElements The number of elements in this block.
        @param numNodesPerElement Length of the numFieldsPerNode list.
        @param numFieldsPerNode Lengths of the rows of the nodalFieldIDs table.
        @param nodalFieldIDs Table where row 'i' is the list of field ids for
               the ith node on every element in this element-block.
        @param numElemDofFieldsPerElement Length of the elemDOFFieldIDs list.
        @param elemDOFFieldIDs list of field identifiers for the element-
               centered degrees-of-freedom in the elements in this block.
        @param interleaveStrategy Indicates the ordering of solution-components
               in the element-wise (e.g., stiffness) contribution arrays. Valid
               values are FEI_NODE_MAJOR (all field-components for first node
               are followed by all field-components for second node, ...) or
               FEI_FIELD_MAJOR (first-field for all nodes, then second-field for
               all nodes, ...)
    */
   virtual int initElemBlock( GlobalID elemBlockID,
                              int numElements,
                              int numNodesPerElement,
                              const int *numFieldsPerNode,
                              const int *const *nodalFieldIDs,
                              int numElemDofFieldsPerElement,
                              const int* elemDOFFieldIDs,
                              int interleaveStrategy ) = 0;
 
    /** Initialize an element's connectivity. Provide a nodal connectivity list
        for inclusion in the sparse matrix structure being constructed.
        @param elemBlockID Which element-block this element belongs to.
        @param elemID A user-provided identifier for this element.
        @param elemConn List of nodeIDs connected to this element. Length of
               this list must be 'numNodesPerElement' provided to the function
               'initElemBlock' for this elemBlockID.
    */
   virtual int initElem( GlobalID elemBlockID,
                         GlobalID elemID, 
                         const GlobalID *elemConn ) = 0; 
 
    /** Identify a list of nodes that are shared by processors other than the
        local processor. This function must be called symmetrically on
        sharing processors -- if a node is identified to processor 0 as being
        shared with processor 1, then that node must also be identified to
        processor 1 as being shared with processor 0. This function may be
        called repeatedly, and shared nodes may appear in multiple
        calls to this function.
        @param numSharedNodes Length of the sharedNodeIDs and numProcsPerNode
               lists, and the number of rows in the sharingProcIDs table.
        @param sharedNodeIDs List of nodes that are shared by this processor and
               one or more other processors.
        @param numProcsPerNode Lengths of the rows of the sharingProcIDs table.
        @param sharingProcIDs List, for each sharedNode, of processors which
               share that node.
    */
   virtual int initSharedNodes( int numSharedNodes,
                                const GlobalID *sharedNodeIDs,  
                                const int *numProcsPerNode,
                                const int *const *sharingProcIDs ) = 0;
 
    /** Constraint relation initialization, Lagrange Multiplier formulation.
        @param numCRNodes Length of the CRNodeIDs and CRFieldIDs lists.
        @param CRNodeIDs Nodes involved in this constraint relation.
        @param CRFieldIDs List of the the field being constrained at each node.
        @param CRID Output. An identifier by which this constraint relation may
               be referred to later, when loading its weight data and recovering
               its Lagrange Multiplier after the solve.
    */
   virtual int initCRMult( int numCRNodes,
                           const GlobalID* CRNodeIDs,
                           const int *CRFieldIDs,
                           int& CRID ) = 0;
 
    /** Constraint relation initialization, Penalty function formulation .
        @param numCRNodes Length of the CRNodeIDs and CRFieldIDs lists.
        @param CRNodeIDs Nodes involved in this constraint relation.
        @param CRFieldIDs List of the the field being constrained at each node.
        @param CRID Output. An identifier by which this constraint relation may
               be referred to later, when loading its weight data and penalty
               value.
    */
   virtual int initCRPen( int numCRNodes,
                          const GlobalID* CRNodeIDs, 
                          const int *CRFieldIDs,
                          int& CRID ) = 0; 

   /** Advise the FEI that a nodal variable is slaved to a linear combination
       of other variables, plus a right-hand-side value (note that the rhsValue
       will often be zero). Since a field may contain more than one scalar
       component, the particular scalar equation that's being slaved must be
       specified by not only a nodeID and fieldID, but also an offset into the
       slave field.

       The general form of the dependency being specified is:
       seqn = sum ( weight_i * meqn_i ) + rhsValue
       where 'seqn' means slave-equation and 'meqn' means master equation.

       Example: to specify that a slave-equation is the average of two master-
       equations: seqn = 0.5*meqn_1 + 0.5*meqn_2 + 0.0
       (Where 0.0 is the rhsValue in this case.)

       The list of weights supplied will be assumed to be of length 
       sum(masterFieldSizes). i.e., the slave equation may be dependent on more
       than one component of the specified master field. In cases where a
       master field contains more than one scalar component, but only one of 
       those components is relevant to the dependency being specified, then 
       positions in the weights list corresponding to the non-relevant 
       field-components should contain zeros.

       This mechanism can also be used as an alternative way to specify
       essential boundary conditions, where the rhsValue will be the boundary
       condition value, with no master nodes or weights.

       Note: This is a new and experimental capability, and is not compatible
       with all other FEI capabilities. In particular, the following precaution
       should be taken:
       Don't identify both slave variables and constraint-relations for the
       same degree of freedom. They are mutually exclusive.

       @param slaveNodeID Node identifier of the node containing the slave eqn.
       @param slaveFieldID Field identifier corresponding to the slave eqn.
       @param offsetIntoSlaveField Denotes location, within the field, of the
       slave eqn.
       @param numMasterNodes Number of nodes containing variables on which the
       slave depends.
       @param masterNodeIDs Node identifiers of the master nodes.
       @param masterFieldIDs List, length numMasterNodes, of the field at each
       master-node which contains the scalar variable(s) on which the slave
       depends.
       @param weights List, length sum-of-master-field-sizes, containing the
       weighting coefficients described above.
       @param rhsValue 
   */
   virtual int initSlaveVariable(GlobalID slaveNodeID,
				 int slaveFieldID,
				 int offsetIntoSlaveField,
				 int numMasterNodes,
				 const GlobalID* masterNodeIDs,
				 const int* masterFieldIDs,
				 const double* weights,
				 double rhsValue) = 0;

   /** Initialization of a pattern by which assembled coefficients will be
       accessed. This is analogous to the 'initElemBlock' function, whereby the
       mask or description of elements in an element-block is provided.
       Locations in the matrix or vector are specified using pairs of either
       nodeIDs and fieldIDs or elementIDs and fieldIDs. For elementIDs, the
       coefficient-access refers to element-centered coefs, or elem-dof.
       An access pattern can be used to specify access to a rectangular region
       of coefficients, with the following limitation. The mask of fields-per-
       column-entity must be the same for each row in the pattern. The number of coefficients in each row of the accessed region will
       be the same.
       @param patternID User-supplied identifier, to be referred to later when
              calling 'initCoefAccess', 'getCoefMatrix', etc.
       @param numRowIDs Number of node IDs or element IDs for which
       corresponding rows will be accessed. (Length of the numFieldsPerRow
       list, and number of rows in the rowFieldIDs table.)
       @param numFieldsPerRow
       @param rowFieldIDs For each row-'entity' (nodeID or elemID), a list of 
       fieldIDs for which rows of coeficients will be accessed.
       @param numColIDsPerRow Number of IDs for which corresponding columns will
              be accessed in each row. (Length of the numFieldsPerCol list, and
	      number of rows in the colFieldIDs table.)
       @param numFieldsPerCol
       @param colFieldIDs For each column-'entity', a list of fieldIDs for which
              columns of the matrix will be accessed.
       @param interleaveStrategy Same valid values here as for the parameter of
              the same name in the function 'initElemBlock'.
   */
   virtual int initCoefAccessPattern( int patternID,
                                      int numRowIDs,
                                      const int* numFieldsPerRow,
                                      const int* const* rowFieldIDs,
                                      int numColIDsPerRow,
                                      const int* numFieldsPerCol,
                                      const int* const* colFieldIDs,
                                      int interleaveStrategy ) = 0;

   /** Initialization of coefficient access for specific rows/columns, using a
       user-defined pattern. This will
       cause positions to be inserted into the underlying sparse matrix 
       structure if they don't already exist.
       @param patternID Access pattern which must already have been defined by
              a call to the function 'initCoefAccessPattern'.
       @param rowIDTypes List whose entries can take the value FEI_NODE or
              FEI_ELEMENT. Entries indicate whether the corresponding entry in
              'rowIDs' is a nodeID or an elementID.
       @param rowIDs List of node or element identifiers for which
       corresponding rows are to be accessed. Length was defined when
       'initCoefAccessPattern' was  called.
       @param colIDTypes List whose entries can take the value FEI_NODE or
              FEI_ELEMENT. Entries indicate whether the corresponding entry in
              'colIDs' is a nodeID or an elementID. The length of this list is
	      assumed to be 'numRowIDs' X 'numColIDsPerRow' (these parameters
	      were supplied when the pattern was initialized via
	      'initCoefAccessPattern').
       @param colIDs List of identifiers for which corresponding columns are to
       be accessed. Length is the same as 'colIDTypes'.
   */ 
   virtual int initCoefAccess( int patternID,
			       const int* rowIDTypes,
			       const GlobalID* rowIDs,
			       const int* colIDTypes,
			       const GlobalID* colIDs ) = 0;

    /** Indicate that initialization phase is complete.
        This function will internally finish calculating the structure of the
        sparse matrix, and provide that structure to the underlying linear
        algebra library. At that point the linear algebra library may or may not
        go ahead and allocate the matrix and vectors for the linear system(s).
    */
   virtual int initComplete() = 0;

     
// Load Phase.......................................... 

   /** Set current matrix data 'context'. When assembling multiple matrices,
     this allows the application to specify which matrix should receive data
     provided by any subsequent data-loading function calls.
     @param matrixID One of the identifiers provided earlier via the function
            'setIDLists'.
    */
   virtual int setCurrentMatrix(int matrixID) = 0;

   /** Set current rhs data 'context'. When assembling multiple rhs vectors,
       this allows the application to specify which rhs should receive data
       provided by any subsequent data-loading function calls.
       @param rhsID One of the identifiers provided earlier via the function
            'setIDLists'.
   */
   virtual int setCurrentRHS(int rhsID) = 0;

    /** Set a value (usually zero) througout the linear system.
        @param s The value to be written into the linear system. 
    */
   virtual int resetSystem(double s=0.0) = 0;

    /** Set a value (usually zero) througout the matrix.
        @param s The value to be written into the matrix.
    */

   virtual int resetMatrix(double s=0.0) = 0;

    /** Set a value (usually zero) througout the rhs vector.
        @param s The value to be written into the rhs vector.
    */
   virtual int resetRHSVector(double s=0.0) = 0;

   /** Set a value (usually, if not always, 0.0) throughout the initial guess
       (solution) vector.
       @param s Input. Scalar value to use in filling the solution vector.
       @return error-code 0 if successful
   */
   virtual int resetInitialGuess(double s) = 0;

   /** Request that any existing Lagrange-Multiplier constraints be deleted.
       (Intended to be called in preparation for loading new/different
       constraints.)
   */
   virtual int deleteMultCRs() = 0;

    /** Load nodal boundary condition data. This allows the application to
       specify a
       boundary condition (dirichlet, neumann or mixed) on a list of nodes.
       The form of these boundary conditions is as follows for a given
       node/field upon which the condition is to be imposed (description
       courtesy of Kim Mish). If the primary
       field solution unknown is denoted by u, and the dual of the solution
       (e.g., force as opposed to displacement) is denoted by q, then a generic
       boundary condition can be specified by alpha*u + beta*q = gamma. A
       dirichlet boundary condition is given when alpha is non-zero but beta is
       zero. A neumann boundary condition is given when alpha is zero but beta
       is non-zero. A mixed boundary condition is given when alpha and beta are
       both non-zero.
       The boundary condition specified via this function applies to the same
       solution field on all nodes in the nodeIDs list.
       @param numNodes Length of the nodeIDs list.
       @param nodeIDs List of nodes upon which a boundary condition is to be
               imposed.
       @param fieldID The solution field that will receive the boundary
                   condition.
       @param alpha Table, with 'numNodes' number-of-rows, and the length of
               each row being the number of scalars that make up the solution
               field 'fieldID'.
       @param beta Table, same dimensions as alpha.
       @param gamma Table, same dimensions as alpha.
    */
   virtual int loadNodeBCs( int numNodes,
                            const GlobalID *nodeIDs,  
                            int fieldID,
                            const double *const *alpha,  
                            const double *const *beta,  
                            const double *const *gamma ) = 0;

    /** Load boundary condition data for element-dof.
        Similar to the function 'loadNodeBCs', as far as the definition of the
        coefficients alpha, beta and gamma.
        @param numElems Length of the elemIDs list.
        @param elemIDs List of elements for which a boundary condition is to be
               specified.
        @param fieldID The solution field for which to apply the boundary
              condition.
        @param alpha Table, as in 'loadNodeBCs', but with 'numElems' number-of-
              rows.
        @param beta Table, same dimensions as alpha.
        @param gamma Table, same dimensions as alpha.
    */
   virtual int loadElemBCs( int numElems,
                            const GlobalID *elemIDs,
                            int fieldID,
                            const double *const *alpha,  
                            const double *const *beta,  
                            const double *const *gamma ) = 0;

    /** Element-stiffness/load data loading. This function accumulates element
       stiffness and load data into the underlying matrix and rhs vector.
       @param elemBlockID Which element-block this element belongs to.
       @param elemID User-supplied identifier for this element.
       @param elemConn Connectivity list of nodes that are connected to this
             element.
       @param elemStiffness Table of element-stiffness data. Dimensions of this
             table defined by the sum of the sizes of the fields associated with
             this element. (This information supplied earlier via
             'initElemBlock'.)
       @param elemLoad Element-load vector.
       @param elemFormat Designates the way in which the 'elemStiffness' 
              stiffness-matrix data is laid out. Valid values for this parameter
              can be found in the file fei_defs.h.
   */
   virtual int sumInElem(GlobalID elemBlockID,
                            GlobalID elemID,  
                            const GlobalID* elemConn,
                            const double* const* elemStiffness,
                            const double* elemLoad,
                            int elemFormat) = 0;
     
    /** Element-stiffness data loading. This function is the same as 'sumInElem'
       but only accepts stiffness data, not the load data for the rhs.
       @param elemBlockID Which element-block this element belongs to.
       @param elemID User-supplied identifier for this element.
       @param elemConn Connectivity list of nodes that are connected to this
             element.
       @param elemStiffness Table of element-stiffness data. Dimensions of this
             table defined by the sum of the sizes of the fields associated with
             this element. (This information supplied earlier via
             'initElemBlock'.)
       @param elemFormat Designates the way in which the 'elemStiffness'
              stiffness-matrix data is laid out. Valid values for this parameter
              can be found in the file fei_defs.h.
   */
   virtual int sumInElemMatrix(GlobalID elemBlockID,
                               GlobalID elemID,  
                               const GlobalID* elemConn,
                               const double* const* elemStiffness,
                               int elemFormat) = 0;
     
    /** Element-load data loading. This function is the same as 'sumInElem',
       but only accepts the load for the rhs, not the stiffness matrix.
       @param elemBlockID Which element-block this element belongs to.
       @param elemID User-supplied identifier for this element.
       @param elemConn Connectivity list of nodes that are connected to this
             element.
       @param elemLoad Element-load vector.
   */
   virtual int sumInElemRHS(GlobalID elemBlockID,
                            GlobalID elemID,  
                            const GlobalID* elemConn,
                            const double* elemLoad) = 0;
     
    /** Load a transfer operator between a particular element two different 
       grids (not currently implemented).
       @param elemBlockID Which element-block this element is associated with.
       @param elemID User-supplied identifier for this element.
       @param coarseNodeList Nodal connectivity in the 'coarse' mesh.
       @param fineNodesPerCoarseElem Length of fineNodeList.
       @param fineNodeList Nodal connectivity in the 'fine' mesh.
       @param elemProlong The prolongation operator.
       @param elemRestrict The restriction operator.
   */
   virtual int loadElemTransfer(GlobalID elemBlockID,
                                GlobalID elemID,
                                const GlobalID* coarseNodeList,
                                int fineNodesPerCoarseElem,
                                const GlobalID* fineNodeList,
                                const double* const* elemProlong,
                                const double* const* elemRestrict) = 0;

    /** Load weight/value data for a Lagrange Multiplier constraint relation.
       @param CRMultID Identifier returned from an earlier call to 'initCRMult'.
       @param numCRNodes Length of CRNodeIDs and CRFieldIDs lists.
       @param CRNodeIDs List of nodes in this constraint relation.
       @param CRFieldIDs List of fields, one per node, to be constrained.
       @param CRWeights Weighting coefficients. This length of this list is the
       sum of the sizes associated with the fields identified in CRFieldIDs.
       @param CRValue The constraint's rhs value. Often (always?) zero.
    */
   virtual int loadCRMult(int CRMultID, int numCRNodes,
                          const GlobalID* CRNodeIDs,
                          const int* CRFieldIDs,
                          const double* CRWeights,
                          double CRValue) = 0;
 
    /** Load weight/value data for a Penalty constraint relation.
       @param CRPenID Identifier returned from an earlier call to 'initCRPen'.
       @param numCRNodes Length of CRNodeIDs and CRFieldIDs lists.
       @param CRNodeIDs List of nodes in this constraint relation.
       @param CRFieldIDs List of fields, one per node, to be constrained.
       @param CRWeights Weighting coefficients. This length of this list is the
       sum of the sizes associated with the fields identified in CRFieldIDs.
       @param CRValue The constraint's rhs value. Often (always?) zero.
       @param penValue The penalty value.
    */
   virtual int loadCRPen(int CRPenID, int numCRNodes,
                         const GlobalID* CRNodeIDs,
                         const int* CRFieldIDs,
                         const double* CRWeights,
                         double CRValue,
                         double penValue) = 0;

   /** Accumulate coefficients into the matrix. 
      @param patternID User-defined identifier used in an earlier call to the
            function 'initCoefAccessPattern'.
      @param rowIDTypes List whose entries can take the value FEI_NODE or
             FEI_ELEMENT. Entries indicate whether the corresponding entry in
             'rowIDs' is a nodeID or an elementID.
      @param rowIDs List of nodes or elements for which the corresponding
      rows in the matrix will be accessed.
       @param colIDTypes List whose entries can take the value FEI_NODE or
              FEI_ELEMENT. Entries indicate whether the corresponding entry in
              'rowIDs' is a nodeID or an elementID.
      @param colIDs List of nodes or elements for which the corresponding
      columns in the matrix will be accessed. This list is assumed to be of
      length 'numRowIDs' X 'numColIDsPerRow' (parameters passed when the
      access-pattern was initialized).
      @param matrixEntries Table of coefficient values. Number of rows given by
           the sum of the sizes of the fields accessed for rowIDs (data
           provided in the earlier call to 'initCoefAccessPattern'). Number of
           columns given by the sum of the sizes of the fields accessed for
           colIDs.
   */
   virtual int sumIntoMatrix(int patternID,
			     const int* rowIDTypes,
			     const GlobalID* rowIDs,
			     const int* colIDTypes,
			     const GlobalID* colIDs,
			     const double* const* matrixEntries) = 0;

   /** Get a copy of coefficients from the matrix. A special semantic applies
       to the getFromMatrix (and getFromRHS) function if any of the IDs in
       the argument 'rowIDs' correspond to shared nodes that are not locally
       owned. Equations for shared nodes reside on only 1 'owning' processor,
       with the owning processor being chosen by the FEI implementation. The
       calling code has no way to know which processor 'owns' a shared node.
       Thus if this call is made only on this processor, then the coefficients
       for equations (matrix rows) that correspond to remotely-owned shared 
       nodes will not be returned. Those positions in the 'matrixEntries'
       argument will not be referenced. If this function call is made on all
       sharing processors, then an attempt will be made to perform the
       communication to get those remotely owned coefficients from the owning
       processor. (This will probably be done using MPI_Iprobe calls to see if
       the owning/sharing processor is in this function and able to be
       communicated with. So not all processors need to make this call, only 
       all processors that share the node(s) for which matrix row entries are
       being requested.)
      @param patternID User-defined identifier used in an earlier call to the
            function 'initCoefAccessPattern'.
       @param rowIDTypes List whose entries can take the value FEI_NODE or
              FEI_ELEMENT. Entries indicate whether the corresponding entry in
              'rowIDs' is a nodeID or an elementID.
      @param rowIDs List of nodes or elements for which the corresponding rows
      in the matrix will be accessed.
       @param colIDTypes List whose entries can take the value FEI_NODE or
              FEI_ELEMENT. Entries indicate whether the corresponding entry in
              'rowIDs' is a nodeID or an elementID.
      @param colIDs List of nodes for which the corresponding columns in the
            matrix will be accessed. This list is assumed to be of
      length 'numRowIDs' X 'numColIDsPerRow' (parameters passed when the
      access-pattern was initialized).
      @param matrixEntries Table of coefficient values. Number of rows given by
           the sum of the sizes of the fields accessed for rowIDs (data
           provided in the earlier call to 'initCoefAccessPattern'). Number of
           columns given by the sum of the sizes of the fields accessed for
           colIDs.
   */
   virtual int getFromMatrix(int patternID,
			     const int* rowIDTypes,
			     const GlobalID* rowIDs,
			     const int* colIDTypes,
			     const GlobalID* colIDs,
			     double** matrixEntries) = 0;

   /** Put a copy of coefficients into the matrix.
      @param patternID User-defined identifier used in an earlier call to the
            function 'initCoefAccessPattern'.
      @param rowIDTypes Same meaning as 'sumIntoMatrix' argument of same name
      @param rowIDs List of nodes for which the corresponding rows in the
            matrix will be accessed.
      @param colIDTypes Same meaning as 'sumIntoMatrix' argument of same name
      @param colIDs List of nodes for which the corresponding columns in the
            matrix will be accessed. This list is assumed to be of
      length 'numRowIDs' X 'numColIDsPerRow' (parameters passed when the
      access-pattern was initialized).
      @param matrixEntries Table of coefficient values. Number of rows given by
           the sum of the sizes of the fields accessed for rowIDs (data
           provided in the earlier call to 'initCoefAccessPattern'). Number of
           columns given by the sum of the sizes of the fields accessed for
           colIDs.
   */
   virtual int putIntoMatrix(int patternID,
			     const int* rowIDTypes,
			     const GlobalID* rowIDs,
			     const int* colIDTypes,
			     const GlobalID* colIDs,
			     const double* const* matrixEntries) = 0;


   /** Accumulate coefficient data into the rhs vector.
      @param patternID User-defined identifier used in an earlier call to the
            function 'initCoefAccessPattern'.
      @param IDTypes Corresponds to the 'rowIDTypes' argument of the matrix-
      access functions
      @param IDs Corresponds to the 'rowIDs' argument of the functions that
            access matrix coefficients.
      @param coefficients List of coefficient values, length is the same as the
            number of rows in the 'matrixEntries' argument of the functions that
           access matrix coefficients.
   */
   virtual int sumIntoRHS(int patternID,
			  const int* IDTypes,
			  const GlobalID* IDs,
			  const double* coefficients) = 0;

                                  
   /** Get a copy of coefficient data from the rhs vector.
      @param patternID User-defined identifier used in an earlier call to the
            function 'initCoefAccessPattern'.
      @param IDTypes
      @param IDs
      @param coefficients List of coefficient values, length is the same as the
            number of rows in the 'matrixEntries' argument of the functions that
           access elem-dof matrix coefficients.
   */
   virtual int getFromRHS(int patternID,
			  const int* IDTypes,
			  const GlobalID* IDs,
			  double* coefficients) = 0;

                                  
   /** Put a copy of coefficient data into the rhs vector.
      @param patternID User-defined identifier used in an earlier call to the
            function 'initCoefAccessPattern'.
      @param IDTypes Corresponds to the 'rowIDTypes' argument of the matrix-
      access functions.
      @param IDs Corresponds to the 'rowIDs' argument of the functions that
            access matrix coefficients.
      @param coefficients List of coefficient values, length is the same as the
            number of rows in the 'matrixEntries' argument of the functions that
           access matrix coefficients.
   */
   virtual int putIntoRHS(int patternID,
			  const int* IDTypes,
                          const GlobalID* IDs,
                          const double* coefficients) = 0;

   /** Put a copy of coefficient data into the rhs vector. */
   virtual int putIntoRHS(int IDType,
			  int fieldID,
			  int numIDs,
			  const GlobalID* IDs,
			  const double* coefficients) = 0;

   /** Sum a copy of coefficient data into the rhs vector. */
   virtual int sumIntoRHS(int IDType,
			  int fieldID,
			  int numIDs,
			  const GlobalID* IDs,
			  const double* coefficients) = 0;

    /** Set scalars by which to multiply matrices, in cases where a linear-
     combination of several matrices is to be solved.
       @param numScalars Length of the IDs and scalars lists.
       @param IDs Matrix identifiers which must be a subset of those supplied
           via an earlier call to 'setIDLists'.
       @param scalars The coefficients by which to multiply the matrices.
    */
   virtual int setMatScalars( int numScalars,
                              const int* IDs,
                              const double* scalars ) = 0;

    /** Set scalars by which to multiply rhs vectors, in cases where a linear-
     combination of several rhs vectors is to be solved.
       @param numScalars Length of the IDs and scalars lists.
       @param IDs RHS-vector identifiers which must be a subset of those
           supplied via an earlier call to 'setIDLists'.
       @param scalars The coefficients by which to multiply the rhs vectors.
    */
   virtual int setRHSScalars( int numScalars,
                              const int* IDs,
                              const double* scalars) = 0;

   /** Indicate that all data loading is complete, and that the underlying
       linear system can now be "finalized". e.g., boundary conditions enforced,
       shared contributions exchanged among processors, etc.
   */
   virtual int loadComplete() = 0;

// Equation solution services..................................... 
 
   /** Request residual norms of the underlying linear-system, broken down
      by solution field. This function should obviously not be called until
      after the matrix and vector data has been fully assembled. Calculates
      the residual vector r = b - Ax, then takes a norm of r, separating out
      the components corresponding to the various fields in 'fieldIDs'. If
      the solution vector x contains zeros, then r == b.
      @param whichNorm Determines which norm is calculated. 1 -> 1-norm, 2 ->
          2-norm, 0 -> infinity-norm.
      @param numFields Length of the fieldIDs and norms lists.
      @param fieldIDs The fields for which residual norms are to be returned.
      @param norms A norm corresponding to each field in fieldIDs.
   */
   virtual int residualNorm(int whichNorm, int numFields,
                             int* fieldIDs, double* norms) = 0;

   /** Launch the solver to solve the assembled linear system.
      @param status Indicates the status of the solve. Varies depending on the
            specific underlying solver library.
      @return err is 0 if the solve was acheived according to any specified
           parameters (i.e., reached specified tolerance within any specified 
                 iteration-limit). Non-zero error return might not indicate a
           fatal error, but rather might indicate that the iteration-limit was
           reached before convergence, etc. A non-zero error return basically
          indicates that the status parameter should be checked and interpreted.
   */
   virtual int solve(int& status) = 0;

 
// Solution return services....................................... 

   /** Query number of iterations taken for last solve.
      @param itersTaken Iterations performed during any previous solve.
   */
   virtual int iterations(int& itersTaken) const = 0;
 
   /**Query the size of a field. This info is supplied to the FEI (initFields)
    by the application, but may not be readily available on the app side at
    all times. Thus, it would be nice if the FEI could answer this query.
   */
   virtual int getFieldSize(int fieldID, int& numScalars) = 0;

   /**Since the ultimate intent for matrix-access is to bypass the FEI and go
     straight to the underlying data objects, we need a translation
     function to map between the IDs that the FEI deals in, and equation
     numbers that linear algebra objects deal in.
     @param ID Identifier of either a node or an element.
     @param idType Can take either of the values FEI_NODE or FEI_ELEMENT.
     @param fieldID Identifies a particular field at this [node||element].
     @param numEqns Output. Number of equations associated with this
     node/field (or element/field) pair.
     @param eqnNumbers Caller-allocated array. On exit, this is filled with the
     above-described equation-numbers. They are global 0-based numbers.
   */
   virtual int getEqnNumbers(GlobalID ID,
			     int idType, 
			     int fieldID,
			     int& numEqns,
			     int* eqnNumbers) = 0;

   /**Get the solution data for a particular field, on an arbitrary set of
      nodes.
      @param fieldID Input. field identifier for which solution data is being
      requested.
      @param numNodes Input. Length of the nodeIDs list.
      @param nodeIDs Input. List specifying the nodes on which solution
      data is being requested.
      @param results Allocated by caller, but contents are output.
      Solution data for the i-th node/element starts in position i*fieldSize,
      where fieldSize is the number of scalar components that make up 'fieldID'.
      @return error-code 0 if successful
   */
   virtual int getNodalFieldSolution(int fieldID,
				     int numNodes,
				     const GlobalID* nodeIDs,
				     double* results) = 0;

   /**Get the number of nodes that are local to this processor (includes nodes
      that are shared by other processors).
      @param numNodes Output. Number of local nodes.
      @return error-code 0 if successful
   */
   virtual int getNumLocalNodes(int& numNodes) = 0;

   /**Get a list of the nodeIDs that are local to this processor (includes nodes
      that are shared by other processors).
      @param numNodes Output. Same as the value output by 'getNumLocalNodes'.
      @param nodeIDs Caller-allocated array, contents to be filled by this
      function.
      @param lenNodeIDs Input. Length of the caller-allocated nodeIDs array. If
      lenNodeIDs is less than numNodes, then only 'lenNodeIDs' nodeIDs are
      provided, of course. If lenNodeIDs is greater than numNodes, then only
      'numNodes' positions of the nodeIDs array are referenced.
      @return error-code 0 if successful
   */
   virtual int getLocalNodeIDList(int& numNodes,
				  GlobalID* nodeIDs,
				  int lenNodeIDs) = 0;

   /** Input data associated with the specified field on a specified list
       of nodes.
   */
   virtual int putNodalFieldData(int fieldID,
				 int numNodes,
				 const GlobalID* nodeIDs,
				 const double* data) = 0;

    /** Return nodal solutions for an element-block. The user supplies the list
        of nodes for which solutions are required. This list may be a subset of
        the nodes in the element-block.
       @param elemBlockID Element-block for which the nodal solutions are
             desired.
       @param numNodes Length of the nodeIDs list.
       @param nodeIDs User-supplied list of nodes in this element-block.
            The list of all nodes in an element-block may be obtained from the
            FEI via the function 'getBlockNodeIDList'.
       @param offsets List of offsets into the results list. The first solution
            value for node i is results[offsets[i]]. The number of solution
            values for node i is offsets[i+1]-offsets[i]. The offsets list is of
            length numNodes+1.
       @param results List of the nodal solution values. The results list should
              be allocated with length 'getNumBlockActEqns', if 'nodeIDs'
              contains all nodes in the element-block.
   */ 
   virtual int getBlockNodeSolution( GlobalID elemBlockID,  
                                     int numNodes,
                                     const GlobalID *nodeIDs, 
                                     int *offsets,  
                                     double *results ) = 0;

   /** Return nodal solutions for an arbitrary list of nodes. The user supplies
       the list of node-identifiers for which solutions are required. This list
       may be any subset of the nodes that reside on the local processor.

       @param numNodes Number of nodes
       @param nodeIDs Node-identifiers, list of length numNodes
       @param offsets List of length numNodes+1, allocated by caller. On exit,
       this list will contain offsets into the results array at which nodal
       solutions are located. The solution values for the nodeIDs[i] will begin
       at position offsets[i]. The last solution value for nodeIDs[i] will be
       located at position offsets[i+1]-1.
       @param results List, allocated by caller, of
       length sum-of-num-dof-per-node.
   */
   virtual int getNodalSolution(int numNodes,
				const GlobalID* nodeIDs,
				int* offsets,
				double* results) = 0;
 
    /** Return nodal solutions for one field.
       @param elemBlockID Element-block identifier.
       @param fieldID The field for which solution values are desired.
       @param numNodes Length of the nodeIDs list.
       @param nodeIDs User-supplied list of nodes for which solution values are
           desired.  The list of all nodes in an element-block may be obtained
           from the FEI via the function 'getBlockNodeIDList'.
       @param results List of solution values. The solution values for node i
           start at results[i*fieldSize] where fieldSize is the size of fieldID.
    */
   virtual int getBlockFieldNodeSolution(GlobalID elemBlockID,
                                          int fieldID,
                                          int numNodes, 
                                          const GlobalID *nodeIDs, 
                                          double *results) = 0;
         
    /** Return elem-dof solution params.
        @param elemBlockID Element-block identifier.
        @param numElems Length of the elemIDs list.
        @param elemIDs User-supplied list of elements for which solution
                    values are desired.
        @param numElemDOFPerElement Output. Number of elem-dof per element.
        @param results List of solution values. The solution values for element 
             i start at results[i*numElemDOFPerElement]
    */
   virtual int getBlockElemSolution( GlobalID elemBlockID,  
                                     int numElems,
                                     const GlobalID *elemIDs,
                                     int& numElemDOFPerElement, 
                                     double *results ) = 0;

   /** Number of constraint relation multipliers. The number of Lagrange
    Multiplier constraint relations in the problem.
     @param numMultCRs Output
   */
   virtual int getNumCRMultipliers(int& numMultCRs) = 0;

   /** List of identifiers for Lagrange Multipliers.
     @param numMultCRs Input. Value obtained from 'getNumCRMultipliers'.
     @param multIDs Output. User-allocated list, to be filled, on exit, with
        the identifiers for the first 'numMultCRs' multiplier constraints in
        the problem.
   */
   virtual int getCRMultIDList(int numMultCRs, int* multIDs) = 0;

   /** Return Lagrange Multiplier solutions.
      @param numCRs number of constraint relations for which multipliers are
        being requested.
      @param CRIDs Identifiers of constraint-relations for which multipliers
           are being requested.
      @param results List of Lagrange Multipliers, one per constraint relation.
   */
   virtual int getCRMultipliers(int numCRs, 
                                 const int *CRIDs,  
                                 double *results) = 0;

 
// Some 'put' methods paralleling the solution 'get' functions.
// the int sizing parameters are passed for error-checking purposes, so
// that the interface implementation can tell if the passed estimate
// vectors make sense -before- an attempt is made to utilize them as
// initial guesses by unpacking them into the solver's native solution
// vector format.

    /** Put nodal-based solution estimates for an element-block.
       @param elemBlockID Element-block identifier.
       @param numNodes Length of nodeIDs list.
       @param nodeIDs Those nodes for which solutions are being supplied.
       @param offsets List, length numNodes+1, of offsets into the estimates
           list.
       @param estimates List of solution estimates. The solution estimates for
            node i will be assumed to start in estimates[offsets[i]]. The
            length of the estimates list should be offsets[numNodes].
    */
   virtual int putBlockNodeSolution(GlobalID elemBlockID, 
                                     int numNodes, 
                                     const GlobalID *nodeIDs, 
                                     const int *offsets,
                                     const double *estimates) = 0;

    /** Put nodal-based solution estimates for one field for an element-block.
       @param elemBlockID Element-block identifier.
       @param fieldID The field for which estimates are being supplied.
       @param numNodes Length of the nodeIDs list.
       @param nodeIDs The nodes for which solution estimates are being supplied.
       @param estimates List of initial guesses. Should be of length numNodes *
           fieldSize, where fieldSize is the size of field 'fieldID'.
    */
   virtual int putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                          int fieldID, 
                                          int numNodes, 
                                          const GlobalID *nodeIDs, 
                                          const double *estimates) = 0;
         
    /** Put element-dof solution guesses for an element-block.
        @param elemBlockID Element-block identifier.
        @param numElems Length of the elemIDs list.
        @param elemIDs Those elements for which elem-dof initial guesses are being
             supplied.
        @param dofPerElem Number of degrees-of-freedom per element.
        @param estimates List of length numElems*dofPerElem. The estimates for element
             i should start in estimates[i*dofPerElem].
    */
   virtual int putBlockElemSolution( GlobalID elemBlockID,  
                                     int numElems, 
                                     const GlobalID *elemIDs, 
                                     int dofPerElem,
                                     const double *estimates) = 0;
  
    /** Put Lagrange Multiplier solution guesses.
       @param numMultCRs Length of the CRMultIDs and multEstimates lists.
       @param CRMultIDs Identifiers obtained from earlier calls to 'initCRMult'.
       @param multEstimates Initial guesses for the Lagrange Multipliers.
    */
   virtual int putCRMultipliers( int numMultCRs,
                                 const int* CRMultIDs,
                                 const double *multEstimates ) = 0;


// utility query functions..............

    /** Return list of nodes associated with an element-block.
        @param elemBlockID Element-block identifier.
        @param numNodes Allocated length of the nodeIDs list. If this is less
               than the value obtained via 'getNumBlockActNodes', then the
               first numNodes node id's will be placed in the nodeIDs list,
               ordered in ascending value.
        @param nodeIDs Output. All nodes associated with this element-block.
    */
   virtual int getBlockNodeIDList( GlobalID elemBlockID,
                                   int numNodes,
                                   GlobalID *nodeIDs ) = 0;
                               
    /** Return llist of all elements associated with an element-block.
        @param elemBlockID Element-block identifier.
        @param numElems Length of the elemIDs list.
        @param elemIDs Output. All elements associated with this element-block.
    */
   virtual int getBlockElemIDList(GlobalID elemBlockID, 
                                   int numElems,
                                   GlobalID *elemIDs ) = 0;
 
// miscellaneous self-explanatory query functions............ 

   /** Return a version string. This string is owned by the FEI implementation,
       the calling application should not delete/free it.
      This string will contain the FEI implementation's version number, and
      if possible, a build time/date.
      @param versionString Output reference to a char*. The C interface will
         have a char** here. This function is simply setting versionString to
         point to the internal version string.
   */
   virtual int version(const char*& versionString) = 0;

   /** Return cumulative cpu-time spent in each of 4 major FEI phases.
     @param initPhase Time in seconds, spent in the initialization phase.
     @param loadPhase Time in seconds, spent loading data, up until the solver
         was launched.
     @param solve Time in seconds, spent in the call to the underlying solver.
     @param solnReturn Time in seconds, spent in the solution-return functions.
   */
   virtual int cumulative_cpu_times(double& initPhase,
				    double& loadPhase,
				    double& solve,
				    double& solnReturn) = 0;
 
   /** Return the number of scalar degrees of freedom associated with a node.
      @param globalNodeID Globally unique node identifier
      @param numSolnParams Sum of the sizes of the solution fields associated
          with the node. This will be the union of the set of fields defined on
         this node over all element-blocks it is associated with.
   */
   virtual int getNumSolnParams( GlobalID globalNodeID,
                                 int& numSolnParams) const = 0;

    /** Return the number of element-blocks.
        @param numElemBlocks
   */
   virtual int getNumElemBlocks(int& numElemBlocks) const = 0;

    /**  Return the number of active nodes in an element-block.
       @param elemBlockID Input.
       @param numNodes Output.
   */
   virtual int getNumBlockActNodes( GlobalID elemBlockID,
                                    int& numNodes) const = 0;

    /**  Return the number of active equations in an element block.
       @param elemBlockID Input.
       @param numEqns Output. Includes both nodal equations and elem-dof
           equations.
   */
   virtual int getNumBlockActEqns( GlobalID elemBlockID,
                                   int& numEqns) const = 0;

    /**  Return the number of nodes associated with elements of an
        element-block.
        @param elemBlockID Input.
        @param nodesPerElem Output.
    */
   virtual int getNumNodesPerElement( GlobalID elemBlockID,
                                      int& nodesPerElem) const = 0;
    
    /** Return the number of equations at each element in an element-block.
       Includes elem-dof equations.
       @param elemBlockID Input.
       @param eqnsPerElem Output.
   */
   virtual int getNumEqnsPerElement( GlobalID elemBlockID,
                                     int& eqnsPerElem) const = 0;

    /**  Return the number of elements in an element-block.
       @param blockID Input.
       @param numElems Output.
    */
   virtual int getNumBlockElements( GlobalID blockID,
                                    int& numElems) const = 0;

    /**  Return the number of element-dof at elements in this element-block.
       @param blockID Input.
       @param DOFPerElem Output.
   */
   virtual int getNumBlockElemDOF( GlobalID blockID,
                                   int& DOFPerElem) const = 0;

};

#endif
