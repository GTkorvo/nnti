
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

#ifndef _EPETRA_VBRMATRIX_H_
#define _EPETRA_VBRMATRIX_H_

#include "Epetra_DistObject.h" 
#include "Epetra_CompObject.h" 
#include "Epetra_BLAS.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsGraph.h"
class Epetra_BlockMap;
class Epetra_Map;
class Epetra_Import;
class Epetra_Export;
class Epetra_Vector;
class Epetra_MultiVector;

//! Epetra_VbrMatrix: A class for constructing and using real-valued double-precision sparse compressed row matrices.

/*! The Epetra_VbrMatrix enable the piecewise construction and use of real-valued double-precision sparse matrices
    where matrix entries are intended for row access.

    At this time, the primary function provided by Epetra_VbrMatrix is matrix time vector and matrix 
    times multi-vector multiplication.  It is also possible to extract matrix rows from a constructed matrix.

<b>Constructing Epetra_VbrMatrix objects</b>

Constructing Epetra_VbrMatrix objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Epetra_VbrMatrix instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
</ol>

Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
create new entries.

<b> Counting Floating Point Operations </b>

Each Epetra_VbrMatrix object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
numbers.  The ResetFlops() function resets the floating point counter.

\warning A Epetra_BlockMap is required for the Epetra_VbrMatrix constructor.

*/    

class Epetra_VbrMatrix: public Epetra_DistObject, public Epetra_CompObject, public Epetra_BLAS, public virtual Epetra_RowMatrix {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Epetra_VbrMatrix constuctor with variable number of indices per row.
  /*! Creates a Epetra_VbrMatrix object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Epetra_BlockMap listing the block rows that this processor will contribute to.
    \param In
           NumBlockEntriesPerRow - An integer array of length NumRows
	   such that NumBlockEntriesPerRow[i] indicates the (approximate) number of Block entries in the ith row.
  */
  Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int *NumBlockEntriesPerRow);
  
  //! Epetra_VbrMatrix constuctor with fixed number of indices per row.
  /*! Creates a Epetra_VbrMatrix object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - An Epetra_BlockMap listing the block rows that this processor will contribute to.
    \param In
           NumBlockEntriesPerRow - An integer that indicates the (approximate) number of Block entries in the each Block row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumBlockEntriesPerRow);

  //! Epetra_VbrMatrix constuctor with variable number of indices per row.
  /*! Creates a Epetra_VbrMatrix object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Epetra_BlockMap listing the block rows that this processor will contribute to.
    \param In 
           ColMap - A Epetra_BlockMap.
    \param In
           NumBlockEntriesPerRow - An integer array of length NumRows
	   such that NumBlockEntriesPerRow[i] indicates the (approximate) number of Block entries in the ith row.
  */
  Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, int *NumBlockEntriesPerRow);
  
  //! Epetra_VbrMatrix constuctor with fixed number of indices per row.
  /*! Creates a Epetra_VbrMatrix object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Epetra_BlockMap listing the block rows that this processor will contribute to.
    \param In 
           ColMap - An Epetra_BlockMap listing the block columns that this processor will contribute to.
    \param In
           NumBlockEntriesPerRow - An integer that indicates the (approximate) number of Block entries in the each Block row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, int NumBlockEntriesPerRow);

  //! Construct a matrix using an existing Epetra_CrsGraph object.
  /*! Allows the nonzero structure from another matrix, or a structure that was
      constructed independently, to be used for this matrix.
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In
           Graph - A Epetra_CrsGraph object, extracted from another Epetra matrix object or constructed directly from
	   using the Epetra_CrsGraph constructors.
  */

  Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph & Graph);

  //! Copy constructor.
  Epetra_VbrMatrix(const Epetra_VbrMatrix & Matrix);

  //! Epetra_VbrMatrix Destructor
  virtual ~Epetra_VbrMatrix();
  //@}
  
  //@{ \name Insertion/Replace/SumInto methods.

  //! Initialize all values in graph of the matrix with constant value.
  /*!
    \param In
           ScalarConstant - Value to use.

    \return Integer error code, set to 0 if successful.
  */
    int PutScalar(double ScalarConstant);

  //! Initiate insertion of a list of elements in a given global row of the matrix, values are inserted via SubmitEntry().
  /*!
    \param In
           BlockRow - Block Row number (in global coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginInsertGlobalValues(int BlockRow, int NumBlockEntries, int * BlockIndices);

  //! Initiate insertion of a list of elements in a given local row of the matrix, values are inserted via SubmitEntry().
  /*!
    \param In
           BlockRow - Block Row number (in local coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginInsertMyValues(int BlockRow, int NumBlockEntries, int * BlockIndices);

  //! Initiate replacement of current values with this list of entries for a given global row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in global coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginReplaceGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

  //! Initiate replacement of current values with this list of entries for a given local row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in local coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginReplaceMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

  //! Initiate summing into current values with this list of entries for a given global row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in global coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginSumIntoGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

  //! Initiate summing into current values with this list of entries for a given local row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in local coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginSumIntoMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

    //! Submit a block entry to the indicated block row and column specified in the Begin routine.
    /* Submit a block entry that will recorded in the block row that was initiated by one of the
       Begin routines listed above.  Once a one of the following routines: BeginInsertGlobalValues(),
       BeginInsertMyValues(), BeginReplaceGlobalValues(), BeginReplaceMyValues(), BeginSumIntoGlobalValues(),
       BeginSumIntoMyValues(), you \e must call SubmitBlockEntry() NumBlockEntries times to register the values 
       corresponding to the block indices passed in to the Begin routine.  If the Epetra_VbrMatrix constuctor
       was called in Copy mode, the values will be copied.  However, no copying will be done until the EndSubmitEntries()
       function is call to complete submission of the current block row.  If the constructor was called in View mode, all
       block entries passed via SubmitBlockEntry() will not be copied, but a pointer will be set to point to the argument Values
       that was passed in by the user.

       For performance reasons, SubmitBlockEntry() does minimal processing of data.  Any processing that can be
       delayed is performed in EndSubmitEntries().

    \param In
           Values - The starting address of the values.
    \param In
           LDA - The stride between successive columns of Values.
    \param In
           NumRows - The number of rows passed in.
    \param In
           NumCols - The number of columns passed in.

    \return Integer error code, set to 0 if successful.
    */
    int SubmitBlockEntry(double *Values, int LDA, int NumRows, int NumCols);

    //! Completes processing of all data passed in for the current block row.
    /*! This function completes the processing of all block entries submitted via SubmitBlockEntry().  
        It also checks to make sure that SubmitBlockEntry was called the correct number of times as
	specified by the Begin routine that initiated the entry process.
    */

    int EndSubmitEntries();

    //! Signal that data entry is complete, perform transformations to local index space.
    /* This version of TransformToLocal assumes that the domain and range distributions are
       identical to the matrix row distributions.
    */
    int TransformToLocal();

    //! Signal that data entry is complete, perform transformations to local index space.
    /* This version of TransformToLocal requires the explicit specification of the domain
       and range distribution maps.  These maps are used for importing and exporting vector
       and multi-vector elements that are needed for distributed matrix computations.  For
       example, to compute y = Ax in parallel, we would specify the DomainMap as the distribution
       of the vector x and the RangeMap as the distribution of the vector y.
    \param In
           DomainMap - Map that describes the distribution of vector and multi-vectors in the
	               matrix domain.
    \param In
           RangeMap - Map that describes the distribution of vector and multi-vectors in the
	               matrix range.
    */
    int TransformToLocal(Epetra_BlockMap *DomainMap, Epetra_BlockMap *RangeMap);

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const {return(Graph_->Filled());};
    //@}

  //@{ \name Extraction methods.

    //! Copy the block indices into user-provided array, set pointers for rest of data for specified global block row.
    /*! 
      This function provides the lightest weight approach to accessing a global block row when the matrix may be
      be stored in local or global index space.  In other words, this function will always work because the block
      indices are returned in user-provided space.  All other array arguments are independent of whether or not
      indices are local or global.  Other than the BlockIndices array, all other array argument are returned as 
      pointers to internal data.

    \param In
           BlockRow - Global block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices array.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted global column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that pointed to by Values.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that is pointed to by Values.
    \param Out
	   Values - Pointer to list of pointers to block entries. Note that the actual values are not copied.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractGlobalBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
				      int & RowDim,  int & NumBlockEntries, 
				      int * BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;

    //! Copy the block indices into user-provided array, set pointers for rest of data for specified local block row.
    /*! 
      This function provides the lightest weight approach to accessing a local block row when the matrix may be
      be stored in local or global index space.  In other words, this function will always work because the block
      indices are returned in user-provided space.  All other array arguments are independent of whether or not
      indices are local or global.  Other than the BlockIndices array, all other array argument are returned as 
      pointers to internal data.

    \param In
           BlockRow - Local block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices array.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted local column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that pointed to by Values.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that is pointed to by Values.
    \param Out
	   Values - Pointer to list of pointers to block entries. Note that the actual values are not copied.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
				       int & RowDim, int & NumBlockEntries, 
				       int * BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;

    //! Initiates a copy of the specified global row in user-provided arrays.
    /*! 
    \param In
           BlockRow - Global block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices, ColDims, and LDAs arrays.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted global column indices for the corresponding block entries.
    \param Out
	   ColDim - List of column dimensions for each corresponding block entry that will be extracted.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractGlobalBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
				       int & RowDim,  int & NumBlockEntries, 
				       int * BlockIndices, int * ColDims) const;

    //! Initiates a copy of the specified local row in user-provided arrays.
    /*! 
    \param In
           BlockRow - Local block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices, ColDims, and LDAs arrays.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted local column indices for the corresponding block entries.
    \param Out
	   ColDim - List of column dimensions for each corresponding block entry that will be extracted.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractMyBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
				       int & RowDim, int & NumBlockEntries, 
				       int * BlockIndices, int * ColDims) const;

    //! Extract a copy of an entry from the block row specified by one of the BeginExtract routines.
    /*! Once BeginExtractGlobalBlockRowCopy() or BeginExtractMyBlockRowCopy() is called, you can extract
        the block entries of specified block row one-entry-at-a-time.  The entries will be extracted
	in an order corresponding to the BlockIndices list that was returned by the BeginExtract routine.

    \param In
           SizeOfValues - Amount of memory associated with Values.  This must be at least as big as
	                  LDA*NumCol, where NumCol is the column dimension of the block entry being copied
    \param InOut
           Values - Starting location where the block entry will be copied.  
    \param In
           LDA - Specifies the stride that will be used when copying columns into Values.
    \param In
           SumInto - If set to true, the block entry values will be summed into existing values.
    */

    int ExtractEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto) const;

    //! Initiates a view of the specified global row, only works if matrix indices are in global mode.
    /*! 
    \param In
           BlockRow - Global block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to global column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
				       int * & BlockIndices, int * & ColDims, int * & LDAs) const;

    //! Initiates a view of the specified local row, only works if matrix indices are in local mode.
    /*! 
    \param In
           BlockRow - Local block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to local column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
				       int * & BlockIndices, int * & ColDims, int * & LDAs) const;


    //! Returns a pointer to, and leading dimension o, the current block entry.
    /*! After a call to BeginExtractGlobal() or BlockRowViewBeginExtractMyBlockRowView(),
        ExtractEntryView() can be called up to NumBlockEntries times to get pointer and stride
	information for each block entry in the specified block row.
    \param InOut
           Values - A pointer that will be set to point to the starting address of the current block entry.
    */
    
    int ExtractEntryView(double * & Values) const;

    //! Initiates a view of the specified global row, only works if matrix indices are in global mode.
    /*! 
    \param In
           BlockRow - Global block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to global column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
    \param Out
	   Values - Pointer to an array of pointers to the block entries in the specified block row.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
				  int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;

    //! Initiates a view of the specified local row, only works if matrix indices are in local mode.
    /*! 
    \param In
           BlockRow - Local block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to local column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
    \param Out
	   Values - Pointer to an array of pointers to the block entries in the specified block row.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
			      int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;


    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
	   Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

    //! Initiates a copy of the block diagonal entries to user-provided arrays.
    /*! 
    \param In
	   MaxNumBlockDiagonalEntries - Length of user-provided RowColDims array.
    \param Out
	   NumBlockDiagonalEntries - Number of block diagonal entries that can actually be extracted.
    \param Out
	   RowColDim - List  of row and column dimension for corresponding block diagonal entries.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractBlockDiagonalCopy(int MaxNumBlockDiagonalEntries, 
				      int & NumBlockDiagonalEntries, int * RowColDims ) const;
    //! Extract a copy of a block diagonal entry from the matrix.
    /*! Once BeginExtractBlockDiagonalCopy() is called, you can extract
        the block diagonal entries one-entry-at-a-time.  The entries will be extracted
	in ascending order.

    \param In
           SizeOfValues - Amount of memory associated with Values.  This must be at least as big as
	                  LDA*NumCol, where NumCol is the column dimension of the block entry being copied
    \param InOut
           Values - Starting location where the block entry will be copied.  
    \param In
           LDA - Specifies the stride that will be used when copying columns into Values.
    \param In
           SumInto - If set to true, the block entry values will be summed into existing values.
    */

    int ExtractBlockDiagonalEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto) const;

    //! Initiates a view of the block diagonal entries.
    /*! 
    \param Out
	   NumBlockDiagonalEntries - Number of block diagonal entries that can be viewed.
    \param Out
	   RowColDim - Pointer to list  of row and column dimension for corresponding block diagonal entries.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractBlockDiagonalView(int & NumBlockDiagonalEntries, int * & RowColDims ) const;

    //! Extract a view of a block diagonal entry from the matrix.
    /*! Once BeginExtractBlockDiagonalView() is called, you can extract a view of
        the block diagonal entries one-entry-at-a-time.  The views will be extracted
	in ascending order.

    \param Out
           Values - Pointer to internal copy of block entry.  
    \param Out
           LDA - Column stride of Values.
    */

    int ExtractBlockDiagonalEntryView(double * & Values, int & LDA) const;
    //@}

  //@{ \name Computational methods.


    //! Returns the result of a Epetra_VbrMatrix multiplied by a Epetra_Vector x in y.
    /*! 
    \param In
	   TransA - If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   x - A Epetra_Vector to multiply by.
    \param Out
	   y - A Epetra_Vector containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const;

    //! Returns the result of a Epetra_VbrMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a solve using the Epetra_VbrMatrix on a Epetra_Vector x in y.
    /*! 
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   x -A Epetra_Vector to solve for.
    \param Out
	   y -A Epetra_Vector containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const;

    //! Returns the result of a Epetra_VbrMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


    //! Computes the sum of absolute values of the rows of the Epetra_VbrMatrix, results returned in x.
    /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvRowSums(Epetra_Vector& x) const;

    //! Scales the Epetra_VbrMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    int LeftScale(const Epetra_Vector& x);

    //! Computes the sum of absolute values of the columns of the Epetra_VbrMatrix, results returned in x.
    /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvColSums(Epetra_Vector& x) const ;

    //! Scales the Epetra_VbrMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    int RightScale(const Epetra_Vector& x);
  //@}

  //@{ \name Matrix Properties Query Methods.

    //! Sort column entries, row-by-row, in ascending order.
    int SortEntries();

    //! If SortEntries() has been called, this query returns true, otherwise it returns false.
    bool Sorted() const {return(Graph_->Sorted());};

    //! Add entries that have the same column index. Remove redundant entries from list.
    int MergeRedundantEntries();

    //! If MergeRedundantEntries() has been called, this query returns true, otherwise it returns false.
    bool NoRedundancies() const {return(Graph_->NoRedundancies());};
    //! Eliminates memory that is used for construction.  Make consecutive row index sections contiguous.
    int OptimizeStorage();

    //! If OptimizeStorage() has been called, this query returns true, otherwise it returns false.
    bool StorageOptimized() const {return(Graph_->StorageOptimized());};

    //! If matrix indices has not been transformed to local, this query returns true, otherwise it returns false.
    bool IndicesAreGlobal() const {return(Graph_->IndicesAreGlobal());};

    //! If matrix indices has been transformed to local, this query returns true, otherwise it returns false.
    bool IndicesAreLocal() const {return(Graph_->IndicesAreLocal());};

    //! If matrix indices are packed into single array (done in OptimizeStorage()) return true, otherwise false.
    bool IndicesAreContiguous() const {return(Graph_->IndicesAreContiguous());};

    //! If matrix is lower triangular, this query returns true, otherwise it returns false.
    bool LowerTriangular() const {return(Graph_->LowerTriangular());};

    //! If matrix is upper triangular, this query returns true, otherwise it returns false.
    bool UpperTriangular() const {return(Graph_->UpperTriangular());};

    //! If matrix is lower triangular, this query returns true, otherwise it returns false.
    bool NoDiagonal() const {return(Graph_->NoDiagonal());};

  //@}
  
  //@{ \name Atribute access functions

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
    */ 
    double NormInf() const;

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1 = \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].
    */ 
    double NormOne() const;


    //! Returns the maximum row dimension of all block entries on this processor.
    int MaxRowDim() const {return(Graph_->MaxRowDim());};

    //! Returns the maximum column dimension of all block entries on this processor.
    int MaxColDim() const {return(Graph_->MaxColDim());};

    //! Returns the maximum row dimension of all block entries across all processors.
    int GlobalMaxRowDim() const {return(Graph_->GlobalMaxRowDim());};

    //! Returns the maximum column dimension of all block entries across all processors.
    int GlobalMaxColDim() const {return(Graph_->GlobalMaxColDim());};

    //! Returns the number of matrix rows owned by the calling processor.
    int NumMyRows() const {return(Graph_->NumMyRows());};
    //! Returns the number of matrix columns owned by the calling processor.
    int NumMyCols() const {return(Graph_->NumMyCols());};

    //! Returns the number of nonzero entriesowned by the calling processor .
    int NumMyNonzeros() const {return(Graph_->NumMyNonzeros());};

    //! Returns the number of global matrix rows.
    int NumGlobalRows() const {return(Graph_->NumGlobalRows());};

    //! Returns the number of global matrix columns.
    int NumGlobalCols() const {return(Graph_->NumGlobalCols());};

    //! Returns the number of nonzero entries in the global matrix.
    int NumGlobalNonzeros() const {return(Graph_->NumGlobalNonzeros());};

    //! Returns the number of Block matrix rows owned by the calling processor.
    int NumMyBlockRows() const {return(Graph_->NumMyBlockRows());};

    //! Returns the number of Block matrix columns owned by the calling processor.
    int NumMyBlockCols() const {return(Graph_->NumMyBlockCols());};
    
    //! Returns the number of nonzero block entries in the calling processor's portion of the matrix.
    int NumMyBlockEntries() const {return(Graph_->NumMyEntries());};

    //! Returns the number of local nonzero block diagonal entries.
    int NumMyBlockDiagonals() const {return(Graph_->NumMyBlockDiagonals());};
    
    //! Returns the number of local nonzero diagonal entries.
    int NumMyDiagonals() const {return(Graph_->NumMyDiagonals());};
    
    //! Returns the number of global Block matrix rows.
    int NumGlobalBlockRows() const {return(Graph_->NumGlobalBlockRows());};
    
    //! Returns the number of global Block matrix columns.
    int NumGlobalBlockCols() const {return(Graph_->NumGlobalBlockCols());};
    
    //! Returns the number of nonzero block entries in the global matrix.
    int NumGlobalBlockEntries() const {return(Graph_->NumGlobalEntries());};
    
    //! Returns the number of global nonzero block diagonal entries.
    int NumGlobalBlockDiagonals() const {return(Graph_->NumGlobalBlockDiagonals());};
    
    //! Returns the number of global nonzero diagonal entries.
    int NumGlobalDiagonals() const {return(Graph_->NumGlobalDiagonals());};

    //! Returns the current number of nonzero Block entries in specified global row on this processor.
    int NumGlobalBlockEntries(int Row) const {return(Graph_->NumGlobalIndices(Row));};

    //! Returns the allocated number of nonzero Block entries in specified global row on this processor.
    int NumAllocatedGlobalBlockEntries(int Row) const{return(Graph_->NumAllocatedGlobalIndices(Row));};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumBlockEntries() const {return(Graph_->MaxNumIndices());};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int GlobalMaxNumBlockEntries() const {return(Graph_->GlobalMaxNumIndices());};

    //! Returns the current number of nonzero Block entries in specified local row on this processor.
    int NumMyBlockEntries(int Row) const {return(Graph_->NumMyIndices(Row));};

    //! Returns the allocated number of nonzero Block entries in specified local row on this processor.
    int NumAllocatedMyBlockEntries(int Row) const {return(Graph_->NumAllocatedMyIndices(Row));};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumNonzeros() const {return(Graph_->MaxNumNonzeros());};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int GlobalMaxNumNonzeros() const {return(Graph_->GlobalMaxNumNonzeros());};

    //! Returns the index base for row and column indices for this graph.
    int IndexBase() const {return(Graph_->IndexBase());};

    //! Returns a pointer to the Epetra_CrsGraph object associated with this matrix.
    const Epetra_CrsGraph & Graph() const {return(*Graph_);};

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    const Epetra_Import * Importer() const {return(Graph_->Importer());};

    //! Returns the Epetra_Export object that contains the export operations for distributed operations.
    const Epetra_Export * Exporter() const {return(Graph_->Exporter());};

    //! Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
    const Epetra_BlockMap & DomainMap() const {return(Graph_->DomainMap());};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    const Epetra_BlockMap & RangeMap() const  {return(Graph_->RangeMap());};

    //! Returns the RowMap object as an Epetra_BlockMap (the Epetra_Map base class) needed for implementing Epetra_RowMatrix.
    const Epetra_BlockMap & RowMap() const {return(Graph_->RowMap());};

    //! Returns the ColMap as an Epetra_BlockMap (the Epetra_Map base class) needed for implementing Epetra_RowMatrix.
    const Epetra_BlockMap & ColMap() const {return(Graph_->ColMap());};

    //! Fills a matrix with rows from a source matrix based on the specified importer.

    //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
    const Epetra_Comm & Comm() const {return(Graph_->Comm());};

  //@}
  
  //@{ \name Local/Global ID methods
    //! Returns the local row index for given global row index, returns -1 if no local row for this global row.
    int LRID( int GRID) const {return(Graph_->LRID(GRID));};

    //! Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
    int GRID( int LRID) const {return(Graph_->GRID(LRID));};

    //! Returns the local column index for given global column index, returns -1 if no local column for this global column.
    int LCID( int GCID) const {return(Graph_->LCID(GCID));};

    //! Returns the global column index for give local column index, returns IndexBase-1 if we don't have this local column.
    int GCID( int LCID) const {return(Graph_->GCID(LCID));};
 
    //! Returns true if the GRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGRID(int GRID) const {return(Graph_->MyGRID(GRID));};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLRID(int LRID) const {return(Graph_->MyLRID(LRID));};

    //! Returns true if the GCID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGCID(int GCID) const {return(Graph_->MyGCID(GCID));};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLCID(int LCID) const {return(Graph_->MyLCID(LCID));};

    //! Returns true of GID is owned by the calling processor, otherwise it returns false.
    bool MyGlobalBlockRow(int GID) const {return(Graph_->MyGlobalRow(GID));};
  //@}
  
  //@{ \name I/O Methods.

  //! Print method
  virtual void Print(ostream & os) const;
  //@}

  //@{ \name Additional methods required to support the Epetra_Operator interface.

    //! Returns a character string describing the operator
    char * Label() const {return(Epetra_Object::Label());};
    
    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! In this implementation, we use several existing attributes to determine how virtual
        method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	the Epetra_VbrMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	be accounted for when doing a triangular solve.

    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns true because this class can compute an Inf-norm.
    bool HasNormInf() const {return(true);};

    //! Returns the current UseTranspose setting.
		bool UseTranspose() const {return(UseTranspose_);};

    //! Returns the Epetra_Map object associated with the domain of this matrix operator.
    const Epetra_Map & OperatorDomainMap() const 
			{ if (!HavePointObjects_) GeneratePointObjects(); return(*OperatorDomainMap_); };

    //! Returns the Epetra_Map object associated with the range of this matrix operator.
    const Epetra_Map & OperatorRangeMap() const 
			{ if (!HavePointObjects_) GeneratePointObjects(); return(*OperatorRangeMap_); };

  //@}
  //@{ \name Additional methods required to implement RowMatrix interface.

    //! Returns a copy of the specified global row in user-provided arrays.
    /*! 
    \param In
           GlobalRow - Global row to extract.
    \param In
	   Length - Length of Values and Indices.
    \param Out
	   NumEntries - Number of nonzero entries extracted.
    \param Out
	   Values - Extracted values for this row.
    \param Out
	   Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractGlobalRowCopy(int GlobalRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Returns a copy of the specified local row in user-provided arrays.
    /*! 
    \param In
           MyRow - Local row to extract.
    \param In
	   Length - Length of Values and Indices.
    \param Out
	   NumEntries - Number of nonzero entries extracted.
    \param Out
	   Values - Extracted values for this row.
    \param Out
	   Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Return the current number of values stored for the specified local row.
    /*! 
    \param In
           MyRow - Local row.
    \param Out
	   NumEntries - Number of nonzero values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int NumMyRowEntries(int MyRow, int & NumEntries) const;

    //! Returns the EpetraMap object associated with the rows of this matrix.
    const Epetra_Map & RowMatrixRowMap() const 
			{ if (!HavePointObjects_) GeneratePointObjects(); return(*RowMatrixRowMap_); };

    //! Returns the Epetra_Map object associated with columns of this matrix.
    const Epetra_Map & RowMatrixColMap() const 
			{ if (!HavePointObjects_) GeneratePointObjects(); return(*RowMatrixColMap_); };

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    const Epetra_Import * RowMatrixImporter() const 
			{ if (!HavePointObjects_) GeneratePointObjects(); return(RowMatrixImporter_); };

  //@}

  //@{ \name Deprecated methods:  These methods still work, but will be removed in a future version.

    //! Use BlockColMap() instead. 
    const Epetra_BlockMap & BlockImportMap() const {return(Graph_->ImportMap());};
  //@}

 protected:
    bool Allocated() const {return(Allocated_);};
    int SetAllocated(bool Flag) {Allocated_ = Flag; return(0);};
    double *** Values() const {return(Values_);};

  // Internal utilities

  void InitializeDefaults();
  int Allocate();
  int BeginInsertValues(int BlockRow, int NumBlockEntries, 
			int * BlockIndices, bool IndicesAreLocal);
  int BeginReplaceValues(int BlockRow, int NumBlockEntries, 
			 int *BlockIndices, bool IndicesAreLocal);
  int BeginSumIntoValues(int BlockRow, int NumBlockEntries, 
			 int *BlockIndices, bool IndicesAreLocal);
  int SetupForSubmits(int BlockRow, int NumBlockEntries, int * BlockIndices, 
		      bool IndicesAreLocal, Epetra_CombineMode SubmitMode);
  int EndReplaceSumIntoValues();
  int EndInsertValues();

  int CopyMat(double * A, int LDA, int NumRows, int NumCols, 
	      double * B, int LDB, bool SumInto) const;
  int BeginExtractBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
			       int & RowDim, int & NumBlockEntries, 
			       int * BlockIndices, int * ColDims, 
			       bool IndicesAreLocal) const;
  int SetupForExtracts(int BlockRow, int & RowDim, int NumBlockEntries,
		       bool ExtractView, bool IndicesAreLocal) const;
  int ExtractBlockDimsCopy(int NumBlockEntries, int * ColDims) const;
  int ExtractBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
				  int & RowDim, int & NumBlockEntries, 
				  int * BlockIndices, int * & ColDims, 
				  int * & LDAs, double ** & Values, bool IndicesAreLocal) const;
  int BeginExtractBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
			       int * & BlockIndices, int * & ColDims, 
			       int * & LDAs, bool IndicesAreLocal) const;
  int ExtractBlockDimsView(int NumBlockEntries, int * & ColDims, int * & LDAs) const;
  int CopyMatDiag(double * A, int LDA, int NumRows, int NumCols, 
		  double * Diagonal) const;

  void BlockRowMultiply(bool TransA, int RowDim, int NumEntries, 
			int * BlockIndices, int RowOff,
			int * FirstPointInElementList, int * ElementSizeList,
			double Alpha, double ** As, int * LDAs, 
			double ** X, double Beta, double ** Y, int NumVectors) const;
  int InverseSums(bool DoRows, Epetra_Vector& x) const;
  int Scale(bool DoRows, const Epetra_Vector& x);
  void BlockRowNormInf(int RowDim, int NumEntries, 
		       int * ColDims, int * LDAs, double ** As, 
		       double * Y) const;
  void BlockRowNormOne(int RowDim, int NumEntries, int * BlockRowIndices,
		       int * ColDims, int * LDAs, double ** As, 
		       int * ColFirstPointInElementList, double * x) const;
  void SetStaticGraph(bool Flag) {StaticGraph_ = Flag;};

  int CheckSizes(const Epetra_DistObject& A){return(0);};
  int CopyAndPermute(const Epetra_DistObject & Source,
		     int NumSameIDs, 
		     int NumPermuteIDs, int * PermuteToLIDs,
		     int *PermuteFromLIDs);
  
  int PackAndPrepare(const Epetra_DistObject & Source,int NumExportIDs, int * ExportLIDs,
				      int Nsend, int Nrecv,
				      int & LenExports, char * & Exports, int & LenImports, 
				      char * & Imports, 
				      int & SizeOfPacket, Epetra_Distributor & Distor);
  
  int UnpackAndCombine(const Epetra_DistObject & Source, 
		       int NumImportIDs, int * ImportLIDs, 
		       char * Imports, int & SizeOfPacket, 
		       Epetra_Distributor & Distor, Epetra_CombineMode CombineMode);

  bool StaticGraph() const {return(StaticGraph_);};

	int GeneratePointObjects() const;
	int BlockMap2PointMap(const Epetra_BlockMap & BlockMap, Epetra_Map * & PointMap) const;
	int Epetra_VbrMatrix::UpdateOperatorXY(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  Epetra_CrsGraph * Graph_;
  bool Allocated_;
  bool StaticGraph_;
  bool UseTranspose_;
  
  int NumMyBlockRows_;

  int * NumEntriesPerRow_;
  int * NumAllocatedEntriesPerRow_;
  Epetra_DataAccess CV_;


  int * NumBlockEntriesPerRow_;
  int * NumAllocatedBlockEntriesPerRow_;
  int ** Indices_;
  int * ElementSizeList_;
  int * FirstPointInElementList_;

  double ***Values_;
  int ** ColDims_;
  int ** LDAs_;
  double **All_Values_;
  mutable double NormInf_;
  mutable double NormOne_;

  mutable Epetra_MultiVector * ImportVector_;
  mutable Epetra_MultiVector * ExportVector_;

  // State variables needed for constructing matrix entry-by-entry
  mutable int *TempRowDims_;
  mutable int *TempColDims_;
  mutable int *TempLDAs_;
  mutable double **TempValues_;
  mutable int LenTemps_;
  mutable int CurBlockRow_;
  mutable int CurNumBlockEntries_;
  mutable int * CurBlockIndices_;
  mutable int CurEntry_;
  mutable bool CurIndicesAreLocal_;
  mutable Epetra_CombineMode CurSubmitMode_;
  
  // State variables needed for extracting entries
  mutable int CurExtractBlockRow_;
  mutable int CurExtractEntry_; 
  mutable int CurExtractNumBlockEntries_;
  mutable bool CurExtractIndicesAreLocal_;
  mutable bool CurExtractView_;
  mutable int CurRowDim_;

  // State variable for extracting block diagonal entries
  mutable int CurBlockDiag_;

	// Maps and importer that support the Epetra_RowMatrix interface
	mutable Epetra_Map * RowMatrixRowMap_;
	mutable Epetra_Map * RowMatrixColMap_;
	mutable Epetra_Import * RowMatrixImporter_;

	// Maps that support the Epetra_Operator interface
	mutable Epetra_Map * OperatorDomainMap_;
	mutable Epetra_Map * OperatorRangeMap_;
	mutable Epetra_MultiVector * OperatorX_;
	mutable Epetra_MultiVector * OperatorY_;

	// bool to indicate if above four point maps and importer have already been created
	mutable bool HavePointObjects_;

};

#endif /* _EPETRA_VBRMATRIX_H_ */
