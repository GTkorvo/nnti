#ifndef CEPETRA_SERIALDENSEMATRIX_CPP_HPP
#define CEPETRA_SERIALDENSEMATRIX_CPP_HPP

/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"

#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_SerialDenseMatrix from non-const table using CT_Epetra_SerialDenseMatrix_ID */
const RCP<Epetra_SerialDenseMatrix>
getSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id );

/*! get Epetra_SerialDenseMatrix from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_SerialDenseMatrix>
getSerialDenseMatrix( CTrilinos_Universal_ID_t id );

/*! get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CT_Epetra_SerialDenseMatrix_ID */
const RCP<const Epetra_SerialDenseMatrix>
getConstSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id );

/*! get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_SerialDenseMatrix>
getConstSerialDenseMatrix( CTrilinos_Universal_ID_t id );

/*! store Epetra_SerialDenseMatrix (owned) in non-const table */
CT_Epetra_SerialDenseMatrix_ID_t
storeNewSerialDenseMatrix( Epetra_SerialDenseMatrix *pobj );

/*! store Epetra_SerialDenseMatrix in non-const table */
CT_Epetra_SerialDenseMatrix_ID_t
storeSerialDenseMatrix( Epetra_SerialDenseMatrix *pobj );

/*! store const Epetra_SerialDenseMatrix in const table */
CT_Epetra_SerialDenseMatrix_ID_t
storeConstSerialDenseMatrix( const Epetra_SerialDenseMatrix *pobj );

/* remove Epetra_SerialDenseMatrix from table using CT_Epetra_SerialDenseMatrix_ID */
void
removeSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t *id );

/* remove Epetra_SerialDenseMatrix from table using CTrilinos_Universal_ID_t */
void
removeSerialDenseMatrix( CTrilinos_Universal_ID_t *aid );

/* purge Epetra_SerialDenseMatrix table */
void
purgeSerialDenseMatrix(  );

/* store Epetra_SerialDenseMatrix in non-const table */
CTrilinos_Universal_ID_t
aliasSerialDenseMatrix( const Teuchos::RCP< Epetra_SerialDenseMatrix > & robj );

/* store const Epetra_SerialDenseMatrix in const table */
CTrilinos_Universal_ID_t
aliasConstSerialDenseMatrix( const Teuchos::RCP< const Epetra_SerialDenseMatrix > & robj );

} // namespace CEpetra


#endif // CEPETRA_SERIALDENSEMATRIX_CPP_HPP


