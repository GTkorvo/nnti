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

%define %komplex_docstring
"
PyTrilinos.Komplex is the python interface to Trilinos package Komplex:

    http://software.sandia.gov/trilinos/packages/komplex

The purpose of Komplex is to provide complex vectors and operators by
encapsulating two Epetra objects to represent real and imaginary
components.  Note that the C++ version of Komplex uses the prefix
'Komplex_' which has been stripped from the python version.

Komplex provides the following user-level classes:

    * MultiVector    - A class for complex-valued double-precison
                       multi-vectors stored in equivalent real formation
    * MultiVector    - A class for complex-valued double-precison
                       vectors stored in equivalent real formation
    * Ordering       - A class for manipulating the KForm of various Komplex
                       objects
    * Operator       - A class for using complex-valued double-precision
                       operators stored in equivalent real formation
    * RowMatrix      - A class for using complex-valued double-precision
                       row matrices stored in equivalent real formation
    * LinearProblem  - A class for explicitly eliminating matrix rows
                       and columns

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exKomplex.py
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %komplex_docstring) Komplex

%{
// System includes
//#include <iostream>
//#include <sstream>
//#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Epetra includes
// #ifdef HAVE_MPI
// #include "Epetra_MpiComm.h"
// #endif
// #include "Epetra_LocalMap.h"
// #include "Epetra_FEVector.h"
// #include "Epetra_Operator.h"
// #include "Epetra_InvOperator.h"
// #include "Epetra_RowMatrix.h"
// #include "Epetra_VbrMatrix.h"
// #include "Epetra_BasicRowMatrix.h"
// #include "Epetra_JadMatrix.h"
// #include "Epetra_JadOperator.h"
// #include "Epetra_FECrsMatrix.h"
// #include "Epetra_FEVbrMatrix.h"

// Epetra python includes
// #include "NumPyImporter.h"
// #include "Epetra_NumPyMultiVector.h"
// #include "Epetra_NumPyVector.h"
// #include "Epetra_NumPyFEVector.h"

// Teuchos Python utility code
//#include "Teuchos_PythonParameter.h"

// Komplex includes
#include "Komplex_Version.h"
#include "Komplex_MultiVector.h"
#include "Komplex_Vector.h"
#include "Komplex_Ordering.h"
#include "Komplex_Operator.h"
#include "Komplex_RowMatrix.h"
#include "Komplex_LinearProblem.h"

%}

// Auto-documentation feature
%feature("autodoc", "1");

// General ignores
%ignore *::operator=;
%ignore *::operator[];
%ignore *::operator();

// External Trilinos modules
%import "Epetra.i"

///////////////////////////////////
// Komplex configuration support //
///////////////////////////////////
%include "Komplex_config.h"

/////////////////////////////
// Komplex Version support //
/////////////////////////////
%include "Komplex_Version.h"
%pythoncode %{
__version__ = Komplex_Version().split()[2]
%}

/////////////////////////////////
// Komplex MultiVector support //
/////////////////////////////////
//%rename(MultiVector) Komplex_MultiVector;
%include "Komplex_MultiVector.h"

////////////////////////////
// Komplex Vector support //
////////////////////////////
 //%rename(Vector) Komplex_Vector;
%include "Komplex_Vector.h"

//////////////////////////////
// Komplex Ordering support //
//////////////////////////////
 //%rename(Ordering) Komplex_Ordering;
%include "Komplex_Ordering.h"

//////////////////////////////
// Komplex Operator support //
//////////////////////////////
 //%rename(Operator) Komplex_Operator;
%include "Komplex_Operator.h"

///////////////////////////////
// Komplex RowMatrix support //
///////////////////////////////
 //%rename(RowMatrix) Komplex_RowMatrix;
%include "Komplex_RowMatrix.h"

///////////////////////////////////
// Komplex LinearProblem support //
///////////////////////////////////
 //%rename(LinearProblem) Komplex_LinearProblem;
%include "Komplex_LinearProblem.h"
