// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef MOOCHOPACK_MOOCHO_THYRA_SOLVER_HPP
#define MOOCHOPACK_MOOCHO_THYRA_SOLVER_HPP

#include "MoochoPack_MoochoSolver.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_DefaultFinalPointCaptureModelEvaluator.hpp"
#include "Thyra_MultiVectorFileIOBase.hpp"
#include "Thyra_ParameterDrivenMultiVectorInput.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief MOOCHO NLP Solver class for models represented through
 * <tt>Thyra::ModelEvaluator</tt>.
 *
 * The parameters this class accepts are shown below in different format:
 * <ul>
 * <li> \ref HumanReadableWithDocumentation "Human readable format (with documentation) for valid parameters accepted by this class"
 * <li> \ref HumanReadableWithoutDocumentation "Human readable format (without documentation) for valid parameters accepted by this class"
 * <li> \ref XmlFormat "XML format for valid parameters accepted by this class"
 * </ul>
 *
 * <b>\anchor HumanReadableWithDocumentation Human readable format (with documentation) for valid parameters accepted by this class</b>
 *
 * \verbinclude NLPThyraEpetraAdvDiffReactOpt.params.readabledoc.out
 *
 * <b>\anchor HumanReadableWithoutDocumentation Human readable format (without documentation) for valid parameters accepted by this class</b>
 *
 * \verbinclude NLPThyraEpetraAdvDiffReactOpt.params.readable.out
 *
 * <b>\anchor XmlFormat XML format for valid parameters accepted by this class</b>
 *
 * \verbinclude NLPThyraEpetraAdvDiffReactOpt.params.xml.out
 *
 * ToDo: Finish documetation!
 */
class MoochoThyraSolver
  : virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \name Public types */
  //@{

  /** \brief . */
  enum ESolveMode {
    SOLVE_MODE_FORWARD    ///< Do a forward solve for the states only.
    ,SOLVE_MODE_OPTIMIZE  ///< Solve an optimization problem for states and parameters.
  };

  //@}

  /** \name Constructors/initialization. */
  //@{

  /** \brief Construct with default settings.
   *
   * <b>Warning!</b> Do not change the defaults by passing then into this
   * constructor.  Instead, use the member functions to set them after
   * <tt>*this</tt> is constructed.  This will help to avoid problems with
   * updates to the ordering of the arguments.
   */
  MoochoThyraSolver(
    const std::string    &paramsXmlFileName                = ""
    ,const std::string   &extraParamsXmlString             = ""
    ,const std::string   &paramsUsedXmlOutFileName         = ""
    ,const std::string   &paramsXmlFileNameOption          = "moocho-thyra-params-file"
    ,const std::string   &extraParamsXmlStringOption       = "extra-moocho-thyra-params"
    ,const std::string   &paramsUsedXmlOutFileNameOption   = "moocho-thyra-params-used-file"
    );
  
  /** \brief . */
  ~MoochoThyraSolver();
  
  /** \brief The name an XML file that will be read to get XML parameters (if
   * not "").
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsXmlFileName);
    
  /** \brief An XML string that will be used to update the parameters (if not
   * "").
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,extraParamsXmlString);

  /** \brief The name of an XML file that will be written (if not "") for the
   * parameters actually used.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsUsedXmlOutFileName);

  /** \brief The name of the option that will be added the the commandline
   * processor that will set <tt>paramsXmlFileName()</tt> .
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsXmlFileNameOption);

  /** \brief The name of the option that will be added the the commandline
   * processor that will set <tt>extraParamsXmlString()</tt> .
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,extraParamsXmlStringOption);

  /** \brief The name of the option that will be added the the commandline
   * processor that will set <tt>paramsUsedXmlOutFileName()</tt> .
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsUsedXmlOutFileNameOption);

  /** \brief MultiVectorFileIOBase object used to read and write state vectors
   * to and from files.
   *
   * The default implementation used is
   * <tt>Thyra::DefaultSpmdMultiVectorFileIO</tt>.  See
   * <tt>Thyra::DefaultSpmdMultiVectorFileIO::DefaultSpmdMultiVectorFileIO()</tt>
   * for the default conditions of such an object.
   */
  STANDARD_NONCONST_COMPOSITION_MEMBERS( Thyra::MultiVectorFileIOBase<value_type>, stateVectorIO );

  /** \brief MultiVectorFileIOBase object used to read and write parameter vectors
   * to and from files.
   *
   * The default implementation used is
   * <tt>Thyra::DefaultSpmdMultiVectorFileIO</tt>.  See
   * <tt>Thyra::DefaultSpmdMultiVectorFileIO::DefaultSpmdMultiVectorFileIO()</tt>
   * for the default conditions.
   *
   * ToDo: You may need to change this to set different the
   * MultiVectorFileIOBase objects for each parameter subvector for each index
   * l=0...Np-1!
   */
  STANDARD_NONCONST_COMPOSITION_MEMBERS( Thyra::MultiVectorFileIOBase<value_type>, parameterVectorIO );
  
  /** \brief Sets up the commandline for reading in the parameter list for
   * this object (minus the options for MoochoSolver).
   *
   * The commandline arguments <tt>--moocho-thyra-params-file</tt> and
   * <tt>--moocho-thyra-extra-params</tt> will be added to <tt>*clp</tt> that
   * will allow the parameters to be read from a file or from the commandline
   * itself in XML format.
   */
  void setupCLP(
    Teuchos::CommandLineProcessor *clp
    );

  /** \brief Force the parameters to be read from a file or from the
   * commandline arguments.
   *
   * If no parameter XML file name or parameters XML string was specified on
   * the commandline, then this function will exist immediately and do
   * nothing.
   *
   * If an XML file name or XML parameter list string was specified on the
   * commandline.
   *
   * If <tt>this->getParameterList().get()==0</tt> before this function is
   * called then 
   * 
   * <b>Postconditions:</b><ul>
   * <li>If the parameters were read from a file or a commandline argument string
   *     then <tt>this->getParameterList()</tt> will be non-NULL and will contain
   *     the updated parameters.
   * 
   */
  void readParameters( std::ostream *out );

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Misc Access/Setup */
  //@

  /** \brief . */
  void setSolveMode( const ESolveMode solveMode );

  /** \brief . */
  ESolveMode getSolveMode() const;
  
  /** \brief . */
  MoochoSolver& getSolver();
  
  /** \brief . */
  const MoochoSolver& getSolver() const;

  //@}

  /** \name Model specification, setup, solve, and solution extraction. */
  //@{
  
  /** \brief . */
  void setModel(
    const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> > &model
    ,const int                                                     p_idx  = 0
    ,const int                                                     g_idx  = 0
    );
    
  /** \brief . */
  void readInitialGuess(
    std::ostream *out = NULL
    );
    
  /** \brief . */
  void setInitialGuess(
    const Teuchos::RefCountPtr<const Thyra::ModelEvaluatorBase::InArgs<value_type> > &initialGuess
    );
    
  /** \brief . */
  void setInitialGuess(
    const Thyra::ModelEvaluatorBase::InArgs<value_type> &initialGuess
    );
  
  /** \brief . */
  MoochoSolver::ESolutionStatus	solve();
    
  /** \brief Return the final point. */
  const Thyra::ModelEvaluatorBase::InArgs<value_type>& getFinalPoint() const;
  
  /** \brief Write the final solution to a file specified by the parameter
   * list option ???.
   */
  void writeFinalSolution(
    std::ostream *out = NULL
    ) const;

  /** \brief Write the parameters list for a this object to a file after the
   * parameters are read in order to show defaults and create a new list for
   * input the next time.
   *
   * If <tt>outputXmlFileName!=""</tt> then the parameter list with be written
   * to the file <tt>outputXmlFileName</tt> in XML format. If
   * <tt>outputXmlFileName==""</tt>, but
   * <tt>this->paramsUsedXmlOutFileNameOption()!=""</tt> then the parameter
   * list will be written to the file
   * <tt>this->paramsUsedXmlOutFileNameOption()</tt>.  If both
   * <tt>outputXmlFileName==""</tt> and
   * <tt>this->paramsUsedXmlOutFileNameOption()==""</tt> then no file is
   * written.
   */
  void writeParamsFile(
    const std::string &outputXmlFileName  = "" 
    ) const;

  //@}

public:

  // I am just making these public so that I can access them from
  // the anonymous namespace in the *.cpp file
  enum ENLPType { NLP_TYPE_FIRST_ORDER, NLP_TYPE_DIRECT };

private:

  typedef value_type Scalar;

  MoochoSolver                                               solver_;
  
  Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >   origModel_;
  int                                                        p_idx_;
  int                                                        g_idx_;

  mutable Thyra::ParameterDrivenMultiVectorInput<value_type> x_reader_;
  mutable Thyra::ParameterDrivenMultiVectorInput<value_type> p_reader_;

  Teuchos::RefCountPtr<Teuchos::ParameterList>               paramList_;
  
  Teuchos::RefCountPtr<Thyra::DefaultNominalBoundsOverrideModelEvaluator<value_type> >
  nominalModel_;
  
  Teuchos::RefCountPtr<Thyra::DefaultFinalPointCaptureModelEvaluator<value_type> >
  finalPointModel_;

  Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> > outerModel_;

  ESolveMode          solveMode_;
  ENLPType            nlpType_;
  bool                nonlinearlyElimiateStates_;
  bool                use_finite_diff_for_obj_;
  bool                use_finite_diff_for_con_;
  double              fwd_newton_tol_;
  int                 fwd_newton_max_iters_;
  bool                useInvObjFunc_;
  bool                showModelEvaluatorTrace_;
  std::string         stateSoluFileBase_;
  std::string         paramSoluFileBase_;

};

} // namespace MoochoPack

#endif	// MOOCHOPACK_MOOCHO_THYRA_SOLVER_HPP
