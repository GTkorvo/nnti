// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//            PyTrilinos.Epetra: Python Interface to Epetra
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%{
// Epetra includes
#include "Epetra_ConfigDefs.h"
#include "Epetra_Version.h"
#include "Epetra_CombineMode.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Flops.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Epetra_MapColoring.h"

// Local includes
#include "Epetra_DArray.h"
#include "Epetra_IArray.h"

// PyTrilinos includes
#include "FileStream.h"
#include "NumPyArray.h"
#include "NumPyWrapper.h"
#include "PyEpetra_Utils.h"  
%}

// Ignore directives
%ignore *::operator=;         // Not overrideable in python
%ignore *::operator[];        // Replaced with __setitem__ method
%ignore *::operator[] const;  // Replaced with __getitem__ method
%ignore *::print;             // Replaced with __str__ method
%ignore operator<<(ostream &, const Epetra_Object &);// From python, use __str__
%ignore NumPyArrayBase::getDataArray() const;
%ignore NumPyArrayBase::getArrayObject() const;
%ignore Epetra_Object::Print(ostream &) const;
%ignore Epetra_MapColoring::Epetra_MapColoring(const Epetra_BlockMap &, int*, const int);
%ignore Epetra_MapColoring::operator()(int) const;
%ignore Epetra_MapColoring::ListOfColors() const;
%ignore Epetra_MapColoring::ColorLIDList(int) const;
%ignore Epetra_MapColoring::ElementColors() const;
%ignore *::UpdateFlops(int) const;   // Use long int version
%ignore *::UpdateFlops(float) const; // Use double version
// These are Import/Export methods that get extended below
%ignore *::PermuteFromLIDs() const;
%ignore *::PermuteToLIDs() const;
%ignore *::RemoteLIDs() const;
%ignore *::ExportLIDs() const;
%ignore *::ExportPIDs() const;

// Rename directives
%rename(Version      ) Epetra_Version;
%rename(Object       ) Epetra_Object;
%rename(SrcDistObject) Epetra_SrcDistObject;
%rename(DistObject   ) Epetra_DistObject;
%rename(CompObject   ) Epetra_CompObject;
%rename(BLAS         ) Epetra_BLAS;
%rename(LAPACK       ) Epetra_LAPACK;
%rename(FLOPS        ) Epetra_Flops;
%rename(Time         ) Epetra_Time;
%rename(Util         ) Epetra_Util;
%rename(MapColoring  ) Epetra_MapColoring;
%rename(DArray       ) Epetra_DArray;
%rename(IArray       ) Epetra_IArray;

// Include directives
%include "Epetra_Version.h"
%include "Epetra_CombineMode.h"
%include "Epetra_DataAccess.h"
%include "Epetra_Object.h"
%include "Epetra_SrcDistObject.h"
%include "Epetra_DistObject.h"
%include "Epetra_CompObject.h"
%import  "Epetra_BLAS.h"       // These two classes are not included because I do not
%import  "Epetra_LAPACK.h"     // want to expose their functionality to python
%include "Epetra_Flops.h"
%include "Epetra_Time.h"
%include "Epetra_Util.h"
%include "Epetra_MapColoring.h"
%include "Epetra_DArray.h"
%include "Epetra_IArray.h"

// Extensions
%extend Epetra_Object {

  // Define the __str__() method, used by the python str() operator on any
  // object given to the python print command.
  string __str__() {
    stringstream os;
    self->Print(os);             // Put the output in os
    string s = os.str();         // Extract the string from os
    int last = s.length();       // Get the last index
    if (s.substr(last) == "\n")
      last-=1;                   // Ignore any trailing newline
    return s.substr(0,last);     // Return the string
  }

  // The Epetra_Object::Print(ostream) method is ignored and replaced by a
  // Print() method here that takes a python file as its argument.  If no
  // argument is given, then output is to standard out.
  PyObject * Print(PyObject*pf=NULL) const {
    if (pf == NULL) {
      self->Print(std::cout);
    } else {
      if (!PyFile_Check(pf)) {
	PyErr_SetString(PyExc_IOError, "Print() method expects file object");
	return NULL;
      } else {
	std::FILE*   f = PyFile_AsFile(pf);
	FileStream   buffer(f);
	std::ostream os(&buffer);
	self->Print(os);
      }
    }
    Py_INCREF(Py_None);
    return Py_None;
  }
}

%extend Epetra_MapColoring {

  Epetra_MapColoring(const Epetra_BlockMap & map,
		     PyObject * elementColors,
		     const int defaultColor=0) {
    Epetra_MapColoring * mapColoring;
    int * colors = 0;
    if (PyInt_Check(elementColors)) {
      int myDefaultColor = (int) PyInt_AsLong(elementColors);
      mapColoring = new Epetra_MapColoring(map,myDefaultColor);
    } else {
      PyArrayObject * colorArray = (PyArrayObject*) PyArray_ContiguousFromObject(elementColors,
										 'i',0,0);
      if (colorArray) colors = (int*)(colorArray->data);
      mapColoring = new Epetra_MapColoring(map,colors,defaultColor);
      Py_XDECREF(colorArray);
    }
    return mapColoring;
  }

  int __getitem__(int i) {
    return self->operator[](i);
  }

  void __setitem__(int i, int color) {
    self->operator[](i) = color;
  }

  PyObject * ListOfColors() {
    int      * list    = self->ListOfColors();
    int        dims[ ] = { self->NumColors() };
    int      * data;
    PyObject * retObj  = PyArray_FromDims(1,dims,'i');
    if (retObj == NULL) goto fail;
    data = (int*)(((PyArrayObject*)(retObj))->data);
    for (int i = 0; i<dims[0]; i++) data[i] = list[i];
    return retObj;
  fail:
    Py_XDECREF(retObj);
    return NULL;
  }

  PyObject * ColorLIDList(int color) {
    int      * list    = self->ColorLIDList(color);
    int        dims[ ] = { self->NumElementsWithColor(color) };
    int      * data;
    PyObject * retObj  = PyArray_FromDims(1,dims,'i');
    if (retObj == NULL) goto fail;
    data = (int*)(((PyArrayObject*)(retObj))->data);
    for (int i = 0; i<dims[0]; i++) data[i] = list[i];
    return retObj;
  fail:
    Py_XDECREF(retObj);
    return NULL;
  }

  PyObject * ElementColors() {
    int      * list    = self->ElementColors();
    int        dims[ ] = { self->Map().NumMyElements() };
    int      * data;
    PyObject * retObj  = PyArray_FromDims(1,dims,'i');
    if (retObj == NULL) goto fail;
    data = (int*)(((PyArrayObject*)(retObj))->data);
    for (int i = 0; i<dims[0]; i++) data[i] = list[i];
    return retObj;
  fail:
    Py_XDECREF(retObj);
    return NULL;
  }
}

// Python code.  Here we set the __version__ string
%pythoncode %{
__version__ = Version().split()[2]
%}
