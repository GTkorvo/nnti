#include "CEpetra_SerialComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_SERIALCOMM_CPP_HPP
#define CEPETRA_SERIALCOMM_CPP_HPP


class Epetra_SerialComm;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_SerialComm>
getSerialComm( CT_Epetra_SerialComm_ID_t id );

void
purgeSerialCommTable(  );


} // namespace CEpetra


#endif // CEPETRA_SERIALCOMM_CPP_HPP

