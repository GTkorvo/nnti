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

%define %aztecoo_docstring
"
PyTrilinos.AztecOO is the python interface to the Trilinos iterative
linear solver package AztecOO:

    http://trilinos.sandia.gov/packages/aztecoo

AztecOO is the object-oriented interface to Aztec, Sandia's venerable
Krylov-space linear system solver package.  Note that the C++ version
of AztecOO uses the prefix 'AztecOO_' which has been stripped from the
python version.  AztecOO requires the Epetra module. The IFPACK and ML
modules can extend the preconditioning capabilities of AztecOO.

AztecOO has a single class:

    * AztecOO  - Object-oriented interface to Aztec package

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exAztecOO.py
    * exAztecOO_Operator.py
    * exAztecOO_RowMatrix.py
    * exAztecOO_BasicRowMatrix.py
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %aztecoo_docstring) AztecOO

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_FEVector.h"
#include "Epetra_InvOperator.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"

// Epetra python includes
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"

// AztecOO includes
#include "AztecOO.h"
#include "AztecOO_Version.h"

// Optional Teuchos support
#ifdef HAVE_AZTECOO_TEUCHOS
#include "Teuchos_PythonParameter.h"
#endif

%}

// SWIG does not support wrapping nested classes.  To suppress the
// swig warning that would otherwise result, we use the following:
#pragma SWIG nowarn=312

// Auto-documentation feature
%feature("autodoc", "1");

// Include AztecOO documentation
%include "AztecOO_dox.i"

// Include the NumPy typemaps
%include "numpy.i"

// External Trilinos interface imports
%import "Epetra.i"
#ifdef HAVE_AZTECOO_TEUCHOS
%import "Teuchos.i"
#endif

// Exception handling
%define %aztecoo_exception(className,methodName)
%exception className::methodName {
  try {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  } catch(int errCode) {
    PyErr_Format(PyExc_RuntimeError, "Error code = %d\nSee stderr for details", errCode);
    SWIG_fail;
  }
}
%enddef

// Macro for methods that return C arrays
%define %aztecoo_return_array(className,methodName,type,typeName,length)
%extend className {
  PyObject * methodName() const {
    intp dims[ ] = { (intp) length };
    return PyArray_SimpleNewFromData(1, dims, typeName, (void*)self->methodName());
  }
}
%ignore className::methodName() const;
%enddef

//////////////////////////////////////
// AztecOO enumerated types support //
//////////////////////////////////////
%include "az_aztec_defs.h"

/////////////////////////////
// AztecOO Version support //
/////////////////////////////
%include "AztecOO_Version.h"
%pythoncode %{
__version__ = AztecOO_Version().split()[2]
%}

/////////////////////
// AztecOO support //
/////////////////////
%aztecoo_exception(AztecOO,SetParameters)
%aztecoo_exception(AztecOO,Iterate      )
%ignore AztecOO::GetAllAztecStatus(double*);
%aztecoo_return_array(AztecOO, GetAllAztecOptions, int,    NPY_INT,    AZ_OPTIONS_SIZE)
%aztecoo_return_array(AztecOO, GetAllAztecParams,  double, NPY_DOUBLE, AZ_PARAMS_SIZE )
%aztecoo_return_array(AztecOO, GetAztecStatus,     double, NPY_DOUBLE, AZ_STATUS_SIZE )
%include "AztecOO.h"
