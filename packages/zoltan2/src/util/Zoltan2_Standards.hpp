// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_Standards.hpp
    \brief Gathering definitions used in software development.

     \todo Should we allow data types for part ID to be set as
         cmake configure options?  Part ID lists in the PartitioningSolution
         are of length "number of objects".  If part ID could be short
         or int, we save significant memory.  For now - typedef'd to int
         so it is easy to change.  It seems data type for proc should
         be int - since it is int in the rest of Trilinos.
*/

#ifndef _ZOLTAN2_STANDARDS_HPP_
#define _ZOLTAN2_STANDARDS_HPP_

#include <Zoltan2_Version.hpp>

//////////////////////////////////////////
// Generated by CMake
#include <Zoltan2_config.h>

//////////////////////////////////////////
// Omit time consuming actions?

#ifdef Z2_OMIT_ALL_OPTIONAL_ACTIONS
#define Z2_OMIT_ALL_STATUS_MESSAGES
#define Z2_OMIT_ALL_PROFILING
#define Z2_OMIT_ALL_ERROR_CHECKING
#endif

//////////////////////////////////////////
// Frequently used Trilinos symbols

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Zoltan2{

using Teuchos::ENull;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_implicit_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::Array;
using Teuchos::Tuple;
using Teuchos::ArrayRCP;
using Teuchos::arcp_const_cast;
using Teuchos::arcp_reinterpret_cast;
using Teuchos::arcp;
using Teuchos::ArrayView;
using Teuchos::av_const_cast;
#ifdef HAVE_ZOLTAN2_MPI
using Teuchos::MpiComm;
#endif
using Teuchos::Comm;
using Teuchos::SerialComm;
using Teuchos::DefaultComm;
using Teuchos::CommRequest;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;
using Teuchos::reduceAll;
using Teuchos::gatherAll;

typedef size_t global_size_t;

}

//////////////////////////////////////////////////////
// Our data types
//   Prepend API types with zoltan2_.
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// For debugging
//////////////////////////////////////////////////////

#define HELLO
//#define HELLO printf("HELLO from %s:%i\n", __FILE__, __LINE__); // Turn on for debug 

//////////////////////////////////////////////////////
// Internal macros and methods
//////////////////////////////////////////////////////

#include <Zoltan2_Exceptions.hpp>

//////////////////////////////////////////////////////
// Until Kokkos node types are supported, use default
//////////////////////////////////////////////////////

#include <Kokkos_DefaultNode.hpp>

namespace Zoltan2{
typedef KokkosClassic::DefaultNode::DefaultNodeType default_node_t;
}

#endif
