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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef NUMPY_INCLUDE_H
#define NUMPY_INCLUDE_H

// This include file takes care of five of the six things necessary
// when including the numpy header file arrayobject.h.  First, the
// Python.h header file is included.  Second, the
// PY_ARRAY_UNIQUE_SYMBOL is defined, which will allow PyTrilinos to
// work with other extension modules that are compiled against NumPy.
// Third, the NPY_NO_DEPRECATED macro is set to NPY_1_7_API_VERSION to
// ensure that no deprecated NumPy code is used.  Fourth, the
// numpy/arrayobject.h header file is included.  Fifth and finally,
// the NPY_API_VERSION macro from arrayobject.h is checked, and if it
// is old enough, macros are defined so that PyTrilinos will compile
// with older versions of NumPy.

// The user is responsible for defining the macro NO_IMPORT_ARRAY in
// those source files that do not call the numpy routine
// import_array().

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL PyTrilinos_NumPy
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#if NPY_API_VERSION < 0x00000007
#define NPY_ANYORDER      PyArray_ANYORDER
#define NPY_DOUBLE        PyArray_DOUBLE
#define NPY_INT           PyArray_INT
#define NPY_ARRAY_FARRAY  NPY_FARRAY
#endif

#endif // NUMPY_INCLUDE_H
