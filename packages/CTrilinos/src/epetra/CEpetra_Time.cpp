
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"

#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_Time.h"
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Time */
Table<Epetra_Time>& tableOfTimes()
{
    static Table<Epetra_Time>
        loc_tableOfTimes(CT_Epetra_Time_ID, "CT_Epetra_Time_ID", FALSE);
    return loc_tableOfTimes;
}

/* table to hold objects of type const Epetra_Time */
Table<const Epetra_Time>& tableOfConstTimes()
{
    static Table<const Epetra_Time>
        loc_tableOfConstTimes(CT_Epetra_Time_ID, "CT_Epetra_Time_ID", TRUE);
    return loc_tableOfConstTimes;
}


} // namespace


//
// Definitions from CEpetra_Time.h
//


extern "C" {


CT_Epetra_Time_ID_t Epetra_Time_Cast ( CTrilinos_Universal_ID_t id )
{
    CTrilinos_Universal_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstTimes(), id);
    } else {
        newid = CTrilinos::cast(tableOfTimes(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(newid);
}

CTrilinos_Universal_ID_t Epetra_Time_Abstract ( 
  CT_Epetra_Time_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id);
}

CT_Epetra_Time_ID_t Epetra_Time_Create ( CT_Epetra_Comm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tableOfTimes().store(new Epetra_Time(*CEpetra::getConstComm(
        CommID))));
}

CT_Epetra_Time_ID_t Epetra_Time_Duplicate ( 
  CT_Epetra_Time_ID_t TimeID )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tableOfTimes().store(new Epetra_Time(*CEpetra::getConstTime(
        TimeID))));
}

void Epetra_Time_Destroy ( CT_Epetra_Time_ID_t * selfID )
{
    CTrilinos_Universal_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Time_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstTimes().remove(&aid);
    } else {
        tableOfTimes().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Time_ID_t>(aid);
}

double Epetra_Time_WallTime ( CT_Epetra_Time_ID_t selfID )
{
    return CEpetra::getConstTime(selfID)->WallTime();
}

void Epetra_Time_ResetStartTime ( CT_Epetra_Time_ID_t selfID )
{
    CEpetra::getTime(selfID)->ResetStartTime();
}

double Epetra_Time_ElapsedTime ( CT_Epetra_Time_ID_t selfID )
{
    return CEpetra::getConstTime(selfID)->ElapsedTime();
}

void Epetra_Time_Assign ( 
  CT_Epetra_Time_ID_t selfID, CT_Epetra_Time_ID_t srcID )
{
    Epetra_Time& self = *( CEpetra::getTime(selfID) );

    self = *CEpetra::getConstTime(srcID);
}


} // extern "C"


//
// Definitions from CEpetra_Time_Cpp.hpp
//


/* get Epetra_Time from non-const table using CT_Epetra_Time_ID */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CT_Epetra_Time_ID_t id )
{
    CTrilinos_Universal_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id);
    return tableOfTimes().get(aid);
}

/* get Epetra_Time from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CTrilinos_Universal_ID_t id )
{
    return tableOfTimes().get(id);
}

/* get const Epetra_Time from either the const or non-const table
 * using CT_Epetra_Time_ID */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CT_Epetra_Time_ID_t id )
{
    CTrilinos_Universal_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id);
    if (id.is_const) {
        return tableOfConstTimes().get(aid);
    } else {
        return tableOfTimes().get(aid);
    }
}

/* get const Epetra_Time from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CTrilinos_Universal_ID_t id )
{
    if (id.is_const) {
        return tableOfConstTimes().get(id);
    } else {
        return tableOfTimes().get(id);
    }
}

/* store Epetra_Time in non-const table */
CT_Epetra_Time_ID_t
CEpetra::storeTime( Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
            tableOfTimes().storeShared(pobj));
}

/* store const Epetra_Time in const table */
CT_Epetra_Time_ID_t
CEpetra::storeConstTime( const Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
            tableOfConstTimes().storeShared(pobj));
}

/* dump contents of Epetra_Time and const Epetra_Time tables */
void
CEpetra::purgeTimeTables(  )
{
    tableOfTimes().purge();
    tableOfConstTimes().purge();
}



