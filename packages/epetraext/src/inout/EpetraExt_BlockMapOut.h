//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
#include <Epetra_ConfigDefs.h>
class Epetra_BlockMap;
namespace EpetraExt {
 
  //! Writes an Epetra_BlockMap or Epetra_Map object to a Matrix Market format file
  /*! This function takes an Epetra_BlockMap or Epetra_Map object and writes it
      to the specified file.  The map can be distributed or serial.  The user can provide
      a strings containing the object name, a description, and specify that header information
      should or should not be printed to the file.

      Special information is encoded in the comment field of this matrix that allows for identical reproduction
      of the map, including distribution across processors and element size information.
      
      The first column of the output file will be the list of GIDs in the map.
      If the block map has non-uniform sizes, a second column will be generated containing the element sizes.

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contained any requested header information
		      followed by the matrix coefficients.  The file will contain a row for each entry.  All entries
		      for a column are listed before going to the next column.
      \param A (In) An Epetra_BlockMap Object containing the user matrix to be dumped to file.
      \param matrixName (In) A C-style string pointer to a name that will be stored in the comment field of the file.
                         This is not a required argument.  Note that it is possible to pass in the method A.Label().
      \param matrixDescription (In) A C-style string pointer to a matrix description that will be stored in the comment 
                                    field of the file.
      \param writeHeader (In) If true, the header will be written, otherwise only the matrix entries will be written.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int BlockMapToMatrixMarketFile( const char *filename, const Epetra_BlockMap & A, 
				   const char * matrixName=0,
				   const char *matrixDescription=0, 
				   bool writeHeader=true);

   

  //! Writes an Epetra_BlockMap or Epetra_Map object to a file handle.
  /*! This function takes an Epetra_BlockMap or Epetra_Map object and writes it
      to the specified file handle.

      \param handle (In) A C-style file handle, already opened.  On exit, the file associated with this handle will
                      have appended to it a row for each multivector row.
      \param A (In) An Epetra_BlockMap object containing the user object to be dumped to file.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int BlockMapToHandle(FILE * handle, const Epetra_BlockMap & A);

  // Internal function
  int writeBlockMap(FILE * handle, int length const int * v1, const int * v2, bool doSizes);

} // namespace EpetraExt
