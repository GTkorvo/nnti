// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef TPETRA_MPITRAITS_HPP
#define TPETRA_MPITRAITS_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_ScalarTraits.hpp>
#include <mpi.h>

/*! \file Tpetra_MpiTraits.hpp
    \brief Defines basic traits for MPI parameters
*/

namespace Tpetra {

  /*! \struct Tpetra::MpiTraits
			\brief This structure defins some basic traits for use with MPI.

      This is a traits file used by MpiComm. It provides traits
      specializations for the built-in types that MPI can handle, 
      and a potentially unsafe workaround for other types.
      
      The built-in types currently supported intrinsicly are:
      int, unsigned int, long int, unsigned long int, float,
      double, long double, and complex<T>.
  */
  

  template<typename T>
  struct UndefinedMpiTraits {
    //! This function should not compile if there is an attempt to instantiate!
    static inline T notDefined() {return(T::this_type_is_missing_a_specialization());};
  };

	template<class T>
	struct MpiTraits {
    //! Returns the MPI_Datatype that should be used for this type
    static inline MPI_Datatype datatype() {return(UndefinedMpiTraits<T>::notDefined());};
    static inline MPI_Op sumOp()           {return(UndefinedMpiTraits<T>::notDefined());};
    static inline MPI_Op maxOp()           {return(UndefinedMpiTraits<T>::notDefined());};
    static inline MPI_Op minOp()           {return(UndefinedMpiTraits<T>::notDefined());};

	};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
	
  //======================================================================
  // intrinsic types
  //======================================================================

	template<>
	struct MpiTraits<int> {
    static inline MPI_Datatype datatype() {return(MPI_INT);};
    static inline MPI_Op sumOp() {return(MPI_SUM);};
    static inline MPI_Op maxOp() {return(MPI_MAX);};
    static inline MPI_Op minOp() {return(MPI_MIN);};
	};

	template<>
	struct MpiTraits<unsigned int> {
    static inline MPI_Datatype datatype() {return(MPI_UNSIGNED);};
    static inline MPI_Op sumOp() {return(MPI_SUM);};
    static inline MPI_Op maxOp() {return(MPI_MAX);};
    static inline MPI_Op minOp() {return(MPI_MIN);};
  };

	template<>
	struct MpiTraits<long int> {
    static inline MPI_Datatype datatype() {return(MPI_LONG);};
    static inline MPI_Op sumOp() {return(MPI_SUM);};
    static inline MPI_Op maxOp() {return(MPI_MAX);};
    static inline MPI_Op minOp() {return(MPI_MIN);};
	};

	template<>
	struct MpiTraits<unsigned long int> {
    static inline MPI_Datatype datatype() {return(MPI_UNSIGNED_LONG);};
    static inline MPI_Op sumOp() {return(MPI_SUM);};
    static inline MPI_Op maxOp() {return(MPI_MAX);};
    static inline MPI_Op minOp() {return(MPI_MIN);};
	};

	template<>
	struct MpiTraits<float> {
    static inline MPI_Datatype datatype() {return(MPI_FLOAT);};
    static inline MPI_Op sumOp() {return(MPI_SUM);};
    static inline MPI_Op maxOp() {return(MPI_MAX);};
    static inline MPI_Op minOp() {return(MPI_MIN);};
	};

	template<>
	struct MpiTraits<double> {
    static inline MPI_Datatype datatype() {return(MPI_DOUBLE);};
    static inline MPI_Op sumOp() {return(MPI_SUM);};
    static inline MPI_Op maxOp() {return(MPI_MAX);};
    static inline MPI_Op minOp() {return(MPI_MIN);};
	};

	template<>
	struct MpiTraits<long double> {
    static inline MPI_Datatype datatype() {return(MPI_LONG_DOUBLE);};
    static inline MPI_Op sumOp() {return(MPI_SUM);};
    static inline MPI_Op maxOp() {return(MPI_MAX);};
    static inline MPI_Op minOp() {return(MPI_MIN);};
	};

  //======================================================================
  // custom specialization for complex<double>
  //======================================================================

  void complexSum(complex<double>* in, complex<double>* inout, int* len, MPI_Datatype* dptr) {
    for(int i = 0; i < *len; i++) {
      *inout += *in;
      in++;
      inout++;
    }
  }
  
  void complexMax(complex<double>* in, complex<double>* inout, int* len, MPI_Datatype* dptr) {
    for(int i = 0; i < *len; i++) {
      if(Teuchos::ScalarTraits< complex<double> >::magnitude(*in) > Teuchos::ScalarTraits< complex<double> >::magnitude(*inout))
        *inout = *in;
      in++;
      inout++;
    }
  }
  
  void complexMin(complex<double>* in, complex<double>* inout, int* len, MPI_Datatype* dptr) {
    for(int i = 0; i < *len; i++) {
      if(Teuchos::ScalarTraits< complex<double> >::magnitude(*in) < Teuchos::ScalarTraits< complex<double> >::magnitude(*inout))
        *inout = *in;
      in++;
      inout++;
    }
  }
  
  // we assume that complex is implemented like a struct, so that the real
  // and imaginary parts are contiguously stored, and there is no other
  // data stored with it. i.e. sizeof(complex<T> == (sizeof(T) * 2)
	template<>
	struct MpiTraits< complex<double> > {

    static MPI_Datatype datatype() {
      MPI_Datatype complex_type;
      MPI_Type_contiguous(2, MPI_DOUBLE, &complex_type);
      MPI_Type_commit(&complex_type);
      return(complex_type);
      };

    static MPI_Op sumOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)complexSum, true, &myOp);
      return(myOp);
    };

    static MPI_Op maxOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)complexMax, true, &myOp);
      return(myOp);
    };

    static MPI_Op minOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)complexMin, true, &myOp);
      return(myOp);
    };

	};

#endif // DOXYGEN_SHOULD_SKIP_THIS
   
} // namespace Tpetra

#endif // TPETRA_MPITRAITS_HPP
