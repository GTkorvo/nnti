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

/*! \file Zoltan2_ColoringProblem.hpp
    \brief Defines the ColoringProblem class.
*/

#ifndef _ZOLTAN2_COLORINGPROBLEM_HPP_
#define _ZOLTAN2_COLORINGPROBLEM_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_ColoringAlgorithms.hpp>
#include <Zoltan2_ColoringSolution.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <string>

#include <bitset>

using Teuchos::rcp_dynamic_cast;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////

/*! \brief ColoringProblem sets up coloring problems for the user.
 *
 *  The ColoringProblem is the core of the Zoltan2 coloring API.
 *  Based on the the user's input and parameters, the ColoringProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the ColoringProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo - Should Problems and Solution have interfaces for returning
 *          views and for returning RCPs?  Or just one?  At a minimum, 
 *          we should have the word "View" in function names that return views.
 *
 *  \todo - Currently, only serial and shared-memory coloring are supported.
 */

template<typename Adapter>
class ColoringProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::zgid_t zgid_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

#ifdef HAVE_ZOLTAN2_MPI
   typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
#endif

  /*! \brief Destructor
   */
  virtual ~ColoringProblem() {};


#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor that takes an MPI communicator
   */
  ColoringProblem(Adapter *A, ParameterList *p, MPI_Comm comm) 
                      : Problem<Adapter>(A, p, comm) 
  {
    HELLO;
    createColoringProblem();
  };
#endif

  /*! \brief Constructor that uses a default communicator
   */
  ColoringProblem(Adapter *A, ParameterList *p) : Problem<Adapter>(A, p) 
  {
    HELLO;
    createColoringProblem();
  };

  //!  \brief Direct the problem to create a solution.
  //
  //    \param updateInputData   If true this indicates that either
  //          this is the first attempt at solution, or that we
  //          are computing a new solution and the input data has
  //          changed since the previous solution was computed.
  //          If false, this indicates that we are computing a
  //          new solution using the same input data was used for
  //          the previous solution, even though the parameters
  //          may have been changed.
  //
  //  For the sake of performance, we ask the caller to set \c updateInputData
  //  to false if he/she is computing a new solution using the same input data,
  //  but different problem parameters, than that which was used to compute
  //  the most recent solution.
  
  void solve(bool updateInputData=true); 

  //!  \brief Get the solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  ColoringSolution<Adapter> *getSolution() {
    // Get the raw ptr from the rcp
    return solution_.getRawPtr();
  };

private:
  void createColoringProblem();

  RCP<ColoringSolution<Adapter> > solution_;

  RCP<Comm<int> > problemComm_;
  RCP<const Comm<int> > problemCommConst_;

#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpiComm_;
#endif
};


////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void ColoringProblem<Adapter>::solve(bool newData)
{
  HELLO;

  size_t nVtx = this->baseModel_->getLocalNumObjects();

  try
  {
      this->solution_ = rcp(new ColoringSolution<Adapter>(nVtx));
  }
  Z2_FORWARD_EXCEPTIONS;

  // Determine which algorithm to use based on defaults and parameters.
  // Need some exception handling here, too.

  std::string method = this->params_->template get<std::string>("color_method", "SerialGreedy");

  try
  {
  // TODO: Ignore case
  if (method.compare("SerialGreedy") == 0)
  {
      AlgSerialGreedy<Adapter> alg(this->graphModel_, this->params_,
                                   this->env_, problemComm_);
      alg.color(this->solution_);
  }
#if 0 // TODO later
  else if (method.compare("speculative") == 0) // Gebremedhin-Manne
  {
      AlgGM<base_adapter_t> alg(this->graphModel_, problemComm_);
      alg.color(this->solution_, this->params_);
  }
#endif
  }
  Z2_FORWARD_EXCEPTIONS;

#ifdef HAVE_ZOLTAN2_MPI

  // The algorithm may have changed the communicator.  Change it back.
  // EGB: This seems excessive. Algorithms should never change the comm?!

  RCP<const mpiWrapper_t > wrappedComm = rcp(new mpiWrapper_t(mpiComm_));
  problemComm_ = rcp(new Teuchos::MpiComm<int>(wrappedComm));
  problemCommConst_ = rcp_const_cast<const Comm<int> > (problemComm_);

#endif

}

////////////////////////////////////////////////////////////////////////
//template <typename Adapter>
//void ColoringProblem<Adapter>::redistribute()
//{
//  HELLO;
//}

////////////////////////////////////////////////////////////////////////
//! createColoringProblem 
//  Method with common functionality for creating a ColoringProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter>
void ColoringProblem<Adapter>::createColoringProblem()
{
  HELLO;
  using Teuchos::ParameterList;

//  cout << __func__ << " input adapter type " 
//       << this->inputAdapter_->inputAdapterType() << " " 
//       << this->inputAdapter_->inputAdapterName() << endl;

  // Create a copy of the user's communicator.

  problemComm_ = this->comm_->duplicate();
  problemCommConst_ = rcp_const_cast<const Comm<int> > (problemComm_);


#ifdef HAVE_ZOLTAN2_MPI

  // TPLs may want an MPI communicator

  Comm<int> *c = problemComm_.getRawPtr();
  Teuchos::MpiComm<int> *mc = dynamic_cast<Teuchos::MpiComm<int> *>(c);
  if (mc){
    RCP<const mpiWrapper_t> wrappedComm = mc->getRawMpiComm();
    mpiComm_ = (*wrappedComm.getRawPtr())();
  }
  else{
    mpiComm_ = MPI_COMM_SELF;   // or would this be an error?
  }

#endif

  // Only graph model supported.
  // TODO: Allow hypergraph later?

  ModelType modelType = GraphModelType; 

  // Select Model based on parameters and InputAdapter type

  std::bitset<NUM_MODEL_FLAGS> graphFlags;
  std::bitset<NUM_MODEL_FLAGS> idFlags;

  switch (modelType) {

  case GraphModelType:
    graphFlags.set(SELF_EDGES_MUST_BE_REMOVED);
    graphFlags.set(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
    this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, problemCommConst_, graphFlags));

    this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->graphModel_);

    break;


  case IdentifierModelType:
  case HypergraphModelType:
  case CoordinateModelType:
    cout << __func__ << " Model type " << modelType << " not yet supported." 
         << endl;
    break;

  default:
    cout << __func__ << " Invalid model" << modelType << endl;
    break;
  }
}
} //namespace Zoltan2

#endif
