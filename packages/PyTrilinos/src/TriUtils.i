// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %triutils_docstring
"
PyTrilinos.TriUtils is the python interface to the Trilinos utilities
package TriUtils:

    http://trilinos.sandia.gov/packages/triutils

The purpose of TriUtils is to provide some utilities typically needed
when testing Trilinos software.

TriUtils provides the following user-level classes:

    * CrsMatrixGallery  - Provide example CrsMatrix objects
    * VbrMatrixGallery  - Provide example VbrMatrix objects

and function:

    * ReadHB            - Obtain a problem from an HB file

For an examples of usage, please consult the following script in the
example subdirectory of the PyTrilinos package:

    * exIFPACK.py
"
%enddef

%module(package = "PyTrilinos",
	autodoc = "1",
	docstring = %triutils_docstring) TriUtils

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"

// Epetra wrapper helper includes
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"

// Trilinos utility includes
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Trilinos_Util_Version.h"
%}

// Auto-documentation feature
%feature("autodoc", "1");

// Include the TriUtils documentation
%include "TriUtils_dox.i"    // Doxygen-generated documentation

// General ignore directives
#pragma SWIG nowarn=503
%ignore *::operator<< ;

// Epetra interface includes
%import "Epetra.i"

///////////////////////////////////
// Trilinos_Util_Version support //
///////////////////////////////////
%rename (TriUtils_Version) Triutils_Version;
%include "Trilinos_Util_Version.h"
%pythoncode %{
__version__ = TriUtils_Version().split()[2]
%}

/////////////////////////////////////////
// Trilinos_Util_ReadHb2Epetra support //
/////////////////////////////////////////
%rename (ReadHB) Trilinos_Util_ReadHb2Epetra;
%include "Trilinos_Util_ReadHb2Epetra.cpp"

////////////////////////////////////////////
// Trilinos_Util_CrsMatrixGallery support //
////////////////////////////////////////////
%ignore Trilinos_Util::CrsMatrixGallery::operator<<(ostream&,
						    const Trilinos_Util::CrsMatrixGallery&);
%include "Trilinos_Util_CrsMatrixGallery.h"
