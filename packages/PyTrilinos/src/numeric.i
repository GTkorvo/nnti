// -*- C -*-  (not really, but good for syntax highlighting)
%{
#ifndef SWIG_FILE_WITH_INIT
#  define NO_IMPORT_ARRAY
#endif
#include "stdio.h"
#include <Numeric/arrayobject.h>

/* The following code originally appeared in
 * enthought/kiva/agg/src/numeric.i, author unknown.  It was
 * translated from C++ to C by John Hunter.  I have modified it
 * slightly to fix some minor bugs.
 */

#define is_array(a)            ((a) && PyArray_Check((PyArrayObject *)a))
#define array_type(a)          (int)(((PyArrayObject *)a)->descr->type_num)
#define array_dimensions(a)    (((PyArrayObject *)a)->nd)
#define array_size(a,i)        (((PyArrayObject *)a)->dimensions[i])
#define array_is_contiguous(a) (PyArray_ISCONTIGUOUS(a))

char* pytype_string(PyObject* py_obj) {
  if (py_obj == NULL          ) return "C NULL value";
  if (PyCallable_Check(py_obj)) return "callable"    ;
  if (PyString_Check(  py_obj)) return "string"      ;
  if (PyInt_Check(     py_obj)) return "int"         ;
  if (PyFloat_Check(   py_obj)) return "float"       ;
  if (PyDict_Check(    py_obj)) return "dict"        ;
  if (PyList_Check(    py_obj)) return "list"        ;
  if (PyTuple_Check(   py_obj)) return "tuple"       ;
  if (PyFile_Check(    py_obj)) return "file"        ;
  if (PyModule_Check(  py_obj)) return "module"      ;
  if (PyInstance_Check(py_obj)) return "instance"    ;

  return "unkown type";
}

char* typecode_string(int typecode) {
  char* type_names[20] = {"char","unsigned byte","byte","short",
			  "unsigned short","int","unsigned int","long",
			  "float","double","complex float","complex double",
			  "object","ntype","unkown"};
  return type_names[typecode];
}

int type_match(int actual_type, int desired_type) {
  // Make sure input has correct numeric type. Allow character and byte to 
  // match also allow int and long to match.
  int match = 1;
  if ( actual_type != desired_type &&
       !(desired_type == PyArray_CHAR  && actual_type == PyArray_SBYTE) &&
       !(desired_type == PyArray_SBYTE && actual_type == PyArray_CHAR)  &&
       !(desired_type == PyArray_INT   && actual_type == PyArray_LONG)  &&
       !(desired_type == PyArray_LONG  && actual_type == PyArray_INT)) 
    match = 0;
  return match;
}

PyArrayObject* obj_to_array_no_conversion(PyObject* input, int typecode) {
  PyArrayObject* ary = NULL;
  if (is_array(input) && (typecode == PyArray_NOTYPE || array_type(input) == typecode)) {
        ary = (PyArrayObject*) input;
    }
    else if is_array(input)
    {
      char* desired_type = typecode_string(typecode);
      char* actual_type = typecode_string(array_type(input));
      PyErr_Format(PyExc_TypeError, 
		   "Array of type '%s' required.  Array of type '%s' given", 
		   desired_type, actual_type);
      ary = NULL;
    }
    else
    {
      char * desired_type = typecode_string(typecode);
      char * actual_type = pytype_string(input);
      PyErr_Format(PyExc_TypeError, 
		   "Array of type '%s' required.  A %s was given", 
		   desired_type, actual_type);
      ary = NULL;
    }
  return ary;
}

PyArrayObject* obj_to_array_allow_conversion(PyObject* input, int typecode,
                                             int* is_new_object)            {
  // Convert object to a Numeric array with the given typecode.
  //
  // Return: 
  //   On Success, return a valid PyArrayObject* with the correct type.
  //   On failure, return NULL.  A python error will have been set.
       
  PyArrayObject* ary = NULL;
  PyObject* py_obj;
  if (is_array(input) && (typecode == PyArray_NOTYPE ||type_match(array_type(input),typecode))) {
    ary = (PyArrayObject*) input;
    *is_new_object = 0;
  }
  else 
  {
    py_obj = PyArray_FromObject(input, typecode, 0, 0);
    // If NULL, PyArray_FromObject will have set python error value.
    ary = (PyArrayObject*) py_obj;
    *is_new_object = 1;
  }
  return ary;
}

PyArrayObject* make_contiguous(PyArrayObject* ary, int* is_new_object,
                               int min_dims, int max_dims)             {
  PyArrayObject* result;
  if (array_is_contiguous(ary)) {
    result = ary;
    *is_new_object = 0;
  }
  else
  {
    result = (PyArrayObject*) PyArray_ContiguousFromObject((PyObject*)ary, 
							   array_type(ary), 
							   min_dims,
							   max_dims);
    *is_new_object = 1;
  }
  return result;
}

PyArrayObject* obj_to_array_contiguous_allow_conversion(PyObject* input,
                                                        int typecode,
                                                        int* is_new_object) {
  int is_new1 = 0;
  int is_new2 = 0;
  PyArrayObject* ary2;
  PyArrayObject* ary1 = obj_to_array_allow_conversion(input, typecode, 
						      &is_new1);
  if (ary1) {
    ary2 = make_contiguous(ary1, &is_new2, 0, 0);
    if ( is_new1 && is_new2) {
      Py_DECREF(ary1);
    }
    ary1 = ary2;    
  }
  *is_new_object = is_new1 || is_new2;
  return ary1;
}

int require_contiguous(PyArrayObject* ary) {
  // Test whether a python object is contiguous.
  //
  // Return:
  //     1 if array is contiguous.
  //     Otherwise, return 0 and set python exception.
  int contiguous = 1;
  if (!array_is_contiguous(ary)) {
    PyErr_SetString(PyExc_TypeError, "Array must be contiguous.  A discontiguous array was given");
    contiguous = 0;
  }
  return contiguous;
}


int require_dimensions(PyArrayObject* ary, int exact_dimensions) {
  int success = 1;
  if (array_dimensions(ary) != exact_dimensions) {
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have %d dimensions.  Given array has %d dimensions", 
		 exact_dimensions, array_dimensions(ary));
    success = 0;
  }
  return success;
}

int require_dimensions_n(PyArrayObject* ary, int* exact_dimensions, int n) {
  int success = 0;
  int i;
  char dims_str[255] = "";
  char s[255];
  for (i = 0; i < n && !success; i++) {
    if (array_dimensions(ary) == exact_dimensions[i]) {
      success = 1;
    }
  }
  if (!success) {
    for (i = 0; i < n-1; i++) {
      sprintf(s, "%d, ", exact_dimensions[i]);                
      strcat(dims_str,s);
    }
    sprintf(s, " or %d", exact_dimensions[n-1]);            
    strcat(dims_str,s);
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have %s dimensions.  Given array has %d dimensions",
		 dims_str, array_dimensions(ary));
  }
  return success;
}    

int require_size(PyArrayObject* ary, int* size, int n) {
  int i;
  int success = 1;
  int len;
  char desired_dims[255] = "[";
  char s[255];
  char actual_dims[255] = "[";
  for(i=0; i < n;i++) {
    if (size[i] != -1 &&  size[i] != array_size(ary,i)) {
      success = 0;    
    }
  }
  if (!success) {
    for (i = 0; i < n; i++) {
      if (size[i] == -1) {
	sprintf(s, "*,");                
      }
      else
      {
	sprintf(s, "%d,", size[i]);                
      }    
      strcat(desired_dims,s);
    }
    len = strlen(desired_dims);
    desired_dims[len-1] = ']';
    for (i = 0; i < n; i++) {
      sprintf(s, "%d,", array_size(ary,i));                            
      strcat(actual_dims,s);
    }
    len = strlen(actual_dims);
    actual_dims[len-1] = ']';
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have shape of %s.  Given array has shape of %s",
		 desired_dims, actual_dims);
  }
  return success;
}
/* End John Hunter translation */

%}

/* TYPEMAP_IN macros
 *
 * This family of typemaps allows pure input C arguments of the form
 *
 *     (type* IN_ARRAY1, int DIM1)
 *     (type* IN_ARRAY2, int DIM1, int DIM2)
 *
 * where "type" is any type supported by the Numeric module, to be
 * called in python with an argument list of a single array (or any
 * python object that can be passed to the Numeric.array constructor
 * to produce an arrayof te specified shape).  This can be applied to
 * a existing functions using the %apply directive:
 *
 *     %apply (double* IN_ARRAY1, int DIM1) {double* series, int length}
 *     %apply (double* IN_ARRAY2, int DIM1, int DIM2) {double* mx, int rows, int cols}
 *     double sum(double* series, int length);
 *     double max(double* mx, int rows, int cols);
 *
 * or with
 *
 *     double sum(double* IN_ARRAY1, int DIM1);
 *     double max(double* IN_ARRAY2, int DIM1, int DIM2);
 */
%define TYPEMAP_IN1(type,typecode)
%typemap(in) (type* IN_ARRAY1, int DIM1)
             (PyArrayObject* array=NULL, int is_new_object) {
  array = obj_to_array_contiguous_allow_conversion($input, typecode, &is_new_object);
  int size[1] = {-1};
  if (!array || !require_dimensions(array,1) || !require_size(array,size,1)) SWIG_fail;
  $1 = (type*) array->data;
  $2 = array->dimensions[0];
}
%typemap(freearg) (type* IN_ARRAY1, int DIM1) {
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}
%enddef

%define TYPEMAP_IN2(type,typecode)
  %typemap(in) (type* IN_ARRAY2, int DIM1, int DIM2)
               (PyArrayObject* array=NULL, int is_new_object) {
  array = obj_to_array_contiguous_allow_conversion($input, typecode, &is_new_object);
  int size[2] = {-1,-1};
  if (!array || !require_dimensions(array,2) || !require_size(array,size,1)) SWIG_fail;
  $1 = (type*) array->data;
  $2 = array->dimensions[0];
  $3 = array->dimensions[1];
}
%typemap(freearg) (type* IN_ARRAY2, int DIM1, int DIM2) {
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}
%enddef

/* Define concrete examples of the TYPEMAP_IN macros */
TYPEMAP_IN1(char,          PyArray_CHAR  )
TYPEMAP_IN1(unsigned char, PyArray_UBYTE )
TYPEMAP_IN1(signed char,   PyArray_SBYTE )
TYPEMAP_IN1(short,         PyArray_SHORT )
TYPEMAP_IN1(int,           PyArray_INT   )
TYPEMAP_IN1(long,          PyArray_LONG  )
TYPEMAP_IN1(float,         PyArray_FLOAT )
TYPEMAP_IN1(double,        PyArray_DOUBLE)
TYPEMAP_IN1(PyObject,      PyArray_OBJECT)
TYPEMAP_IN2(char,          PyArray_CHAR  )
TYPEMAP_IN2(unsigned char, PyArray_UBYTE )
TYPEMAP_IN2(signed char,   PyArray_SBYTE )
TYPEMAP_IN2(short,         PyArray_SHORT )
TYPEMAP_IN2(int,           PyArray_INT   )
TYPEMAP_IN2(long,          PyArray_LONG  )
TYPEMAP_IN2(float,         PyArray_FLOAT )
TYPEMAP_IN2(double,        PyArray_DOUBLE)
TYPEMAP_IN2(PyObject,      PyArray_OBJECT)

#undef TYPEMAP_IN1
#undef TYPEMAP_IN2

/* TYPEMAP_INOUT macros
 *
 * This family of typemaps allows input/output C arguments of the form
 *
 *     (type* INOUT_ARRAY1, int DIM1)
 *     (type* INOUT_ARRAY2, int DIM1, int DIM2)
 *
 * where "type" is any type supported by the Numeric module, to be
 * called in python with an argument list of a single contiguous
 * Numeric array.  This can be applied to an existing function using
 * the %apply directive:
 *
 *     %apply (double* INOUT_ARRAY1, int DIM1) {double* series, int length}
 *     %apply (double* INOUT_ARRAY2, int DIM1, int DIM2) {double* mx, int rows, int cols}
 *     void negate(double* series, int length);
 *     void normalize(double* mx, int rows, int cols);
 *     
 *
 * or with
 *
 *     void sum(double* INOUT_ARRAY1, int DIM1);
 *     void sum(double* INOUT_ARRAY1, int DIM1, int DIM2);
 */
%define TYPEMAP_INOUT1(type,typecode)
%typemap(in) (type* INOUT_ARRAY1, int DIM1) (PyArrayObject* temp=NULL) {
  temp = obj_to_array_no_conversion($input,typecode);
  if (!temp) SWIG_fail;
  if (!require_contiguous(temp)) SWIG_fail;
  $1 = (type*) temp->data;
  $2 = 1;
  for (int i=0; i<temp->nd; ++i) $2 *= temp->dimensions[i];
}
%enddef

%define TYPEMAP_INOUT2(type,typecode)
%typemap(in) (type* INOUT_ARRAY1, int DIM1) (PyArrayObject* temp=NULL) {
  temp = obj_to_array_no_conversion($input,typecode);
  if (!temp) SWIG_fail;
  if (!require_contiguous(temp)) SWIG_fail;
  $1 = (type*) temp->data;
  $2 = 1;
  for (int i=0; i<temp->nd; ++i) $2 *= temp->dimensions[i];
}
%enddef

/* Define concrete examples of the TYPEMAP_OUT macro */
TYPEMAP_INOUT1(char,          PyArray_CHAR  )
TYPEMAP_INOUT1(unsigned char, PyArray_UBYTE )
TYPEMAP_INOUT1(signed char,   PyArray_SBYTE )
TYPEMAP_INOUT1(short,         PyArray_SHORT )
TYPEMAP_INOUT1(int,           PyArray_INT   )
TYPEMAP_INOUT1(long,          PyArray_LONG  )
TYPEMAP_INOUT1(float,         PyArray_FLOAT )
TYPEMAP_INOUT1(double,        PyArray_DOUBLE)
TYPEMAP_INOUT1(PyObject,      PyArray_OBJECT)
TYPEMAP_INOUT2(char,          PyArray_CHAR  )
TYPEMAP_INOUT2(unsigned char, PyArray_UBYTE )
TYPEMAP_INOUT2(signed char,   PyArray_SBYTE )
TYPEMAP_INOUT2(short,         PyArray_SHORT )
TYPEMAP_INOUT2(int,           PyArray_INT   )
TYPEMAP_INOUT2(long,          PyArray_LONG  )
TYPEMAP_INOUT2(float,         PyArray_FLOAT )
TYPEMAP_INOUT2(double,        PyArray_DOUBLE)
TYPEMAP_INOUT2(PyObject,      PyArray_OBJECT)

#undef TYPEMAP_INOUT1
#undef TYPEMAP_INOUT2
