/*
 * Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf
#include <sys/types.h>                  // for int64_t
#include "exodusII.h"                   // for exerrval, ex_err, etc
#include "exodusII_int.h"               // for EX_FATAL, EX_NOERR, etc
#include "netcdf.h"                     // for NC_NOERR, nc_get_vara_int, etc

/*
 * reads the id map
 */

int ex_get_partial_id_map ( int   exoid,
			    ex_entity_type map_type,
			    int64_t   start_entity_num,
			    int64_t   num_entities,
			    void_int*  map )
{
  int dimid, mapid, status;
  size_t i;
  size_t num_entries;
  size_t start[1], count[1];
  char errmsg[MAX_ERR_LENGTH];
  const char* dnumentries;
  const char* vmap;
  const char* tname;

  switch (map_type) {
  case EX_NODE_MAP:
    tname = "node";
    dnumentries = DIM_NUM_NODES;
    vmap = VAR_NODE_NUM_MAP;
    break;
  case EX_EDGE_MAP:
    tname = "edge";
    dnumentries = DIM_NUM_EDGE;
    vmap = VAR_EDGE_NUM_MAP;
    break;
  case EX_FACE_MAP:
    tname = "face";
    dnumentries = DIM_NUM_FACE;
    vmap = VAR_FACE_NUM_MAP;
    break;
  case EX_ELEM_MAP:
    tname = "element";
    dnumentries = DIM_NUM_ELEM;
    vmap = VAR_ELEM_NUM_MAP;
    break;
  default:
    exerrval = EX_BADPARAM;
    sprintf( errmsg,
	     "Error: Bad map type (%d) specified for file id %d",
	     map_type, exoid );
    ex_err( "ex_get_partial_id_map", errmsg, exerrval );
    return (EX_FATAL);
  }
  exerrval = 0; /* clear error code */

  /* See if any entries are stored in this file */
  if (nc_inq_dimid(exoid, dnumentries,&dimid) != NC_NOERR) {
    return (EX_NOERR);
  }

  if (nc_inq_varid (exoid, vmap, &mapid) != NC_NOERR) {
    if ((status = nc_inq_dimlen(exoid, dimid, &num_entries)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to get number of %ss in file id %d", tname, exoid);
      ex_err("ex_get_partial_id_map",errmsg,exerrval);
      return (EX_FATAL);
    }
    
    /* generate default map of 1..n, where n is num_entries */
    if (ex_int64_status(exoid) & EX_MAPS_INT64_API) {
      int64_t *lmap = (int64_t*)map;
      for (i=0; i < num_entities; i++) {
	lmap[i] = start_entity_num+i;
      }
    } else {
      int *lmap = (int*)map;
      for (i=0; i<num_entities; i++) {
	lmap[i] = start_entity_num+i;
      }
    }

    return (EX_NOERR);
  }

  start[0] = start_entity_num-1;
  count[0] = num_entities;

  /* read in the id map  */
  if (ex_int64_status(exoid) & EX_MAPS_INT64_API) {
    status = nc_get_vara_longlong(exoid, mapid, start, count, map);
  } else {
    status = nc_get_vara_int(exoid, mapid, start, count, map);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get %s id map in file id %d",
	    tname, exoid);
    ex_err("ex_get_partial_id_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  return(EX_NOERR);
}
