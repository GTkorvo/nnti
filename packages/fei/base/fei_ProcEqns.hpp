/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _ProcEqns_hpp_
#define _ProcEqns_hpp_


/** Internal implementation class.
A class for keeping equation numbers grouped according to processor. This
is useful when a set of equations is to be exchanged among processors. This
class holds a table of equation numbers, and optionally a companion table with
the lengths of those equations, organized such that each row of the table 
contains equations that are associated with a remote processor. (They are to be
sent to that proc, or they are to be received from that proc.) A list of those
processors is also maintained, of course.

Usage of this class is intended to consist of adding equations and associated
processors using the addEqn member, and then later retrieving the list of
processors and the associated table (list of lists) of equations, and their
lengths (if length data was supplied).
*/
 
class ProcEqns {
 public:
  /** Default constructor */
   ProcEqns();
   /** Destructor */
   virtual ~ProcEqns();

   /** Similar to a copy constructor, puts a copy of all data and state into
       the new instance.
   */
   ProcEqns* deepCopy();

   /** Return the number of processors for which equations are held. */
   size_t getNumProcs() {return(procs_.size());}

   /** Return a list of processors. */
   std::vector<int>& procsPtr() {return(procs_);}

   /** Return a list containing the number of equations corresponding to each
       processor. The length of this list should be 'getNumProcs()'.
   */
   std::vector<int>& eqnsPerProcPtr() {return(eqnsPerProc_);}

   /** Table of equation numbers. number-of-rows = 'getNumProcs()', length of
       row i is eqnsPerProcPtr()[i].
   */
   std::vector<std::vector<int>*>& procEqnNumbersPtr() {return(procEqnNumbers_);}

   /** Table containing the lengths of the equations in 'procEqnNumbersPtr()'.
       Returns NULL if no length data is present.
    */
   std::vector<std::vector<int>*>& procEqnLengthsPtr() {return(procEqnLengths_);}

   /** Add an equation-number and associated processor number to the internal
       data structures. Equations may be added for multiple different procs.
   */
   void addEqn(int eqnNumber, int proc);

   /** Add an equation-number/equation-length pair, and associated processor
       number, to the internal data structures. Equations may be added for
       multiple different procs.
   */
   void addEqn(int eqnNumber, int eqnLength, int proc);

   /** Replace the 'procEqnLengthsPtr()' table with these lengths. There is a
       large potential for user error here. The number of equation/length pairs
       being provided here should equal the total number of equations already
       identified to this object. The internal table of lengths will be 
       destroyed if it already exists, and replaced by this incoming data.
       @param eqnNumbers Equations for which lengths are being provided. These
       equations must already have been identified and associated with procs.
       @param eqnLengths Equation lengths.
       @param len The number of equations and lengths.
   */
   void setProcEqnLengths(int* eqnNumbers, int* eqnLengths, int len);

 private:
   void deleteMemory();
   void internalAddEqn(int eqnNumber, int eqnLength, int proc);

   std::vector<int> procs_;
   std::vector<int> eqnsPerProc_;
   std::vector<std::vector<int>* > procEqnNumbers_;
   std::vector<std::vector<int>* > procEqnLengths_;
};

#endif

