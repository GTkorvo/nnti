// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_Factory.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Parrot.hpp"
#include "Shards_CellTopology.hpp"

// *******************************************************************

void panzer::buildPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
				const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_ids_to_cell_topo,
				const Teuchos::RCP<Teuchos::ParameterList>& physics_blocks_plist,
				const int default_integration_order,
				const std::size_t workset_size,
				const panzer::EquationSetFactory & eqset_factory,
				const Teuchos::RCP<panzer::GlobalData>& global_data,
				const bool build_transient_support,
				std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::map;
  using std::string;
  
  TEUCHOS_ASSERT(nonnull(physics_blocks_plist));

  // Create a physics block for each element block
  map<string,string>::const_iterator itr;
  for (itr = block_ids_to_physics_ids.begin(); itr!=block_ids_to_physics_ids.end();++itr) {
    string element_block_id = itr->first;
    string physics_block_id = itr->second;
    
    map<string,RCP<const shards::CellTopology> >::const_iterator ct_itr =
      block_ids_to_cell_topo.find(element_block_id);
    TEUCHOS_TEST_FOR_EXCEPTION(ct_itr==block_ids_to_cell_topo.end(),
			       std::runtime_error,
			       "Falied to find CellTopology for element block id: \""
			       << element_block_id << "\"!");
    RCP<const shards::CellTopology> cellTopo = ct_itr->second; 
    
    const panzer::CellData volume_cell_data(workset_size,cellTopo);
    
    // find physics block parameter sublist
    TEUCHOS_TEST_FOR_EXCEPTION(physics_blocks_plist->isSublist(physics_block_id),
			       std::runtime_error,
			       "Failed to find physics id: \""
			       << physics_block_id 
			       << "\" requested by element block: \"" 
			       << element_block_id << "\"!");
    
    RCP<panzer::PhysicsBlock> pb = rcp(new panzer::PhysicsBlock(Teuchos::sublist(physics_blocks_plist,physics_block_id,true), 
								element_block_id,
								default_integration_order,
								volume_cell_data, 
								eqset_factory, 
								global_data,
								build_transient_support));
    physicsBlocks.push_back(pb);
  }
}

// *******************************************************************
panzer::PhysicsBlock::
PhysicsBlock(const Teuchos::RCP<Teuchos::ParameterList>& physics_block_plist,
             const std::string & element_block_id,
	     const int default_integration_order,
	     const panzer::CellData & cell_data,
	     const panzer::EquationSetFactory& factory,
	     const Teuchos::RCP<panzer::GlobalData>& global_data,
	     const bool build_transient_support) :
  m_element_block_id(element_block_id),
  m_default_integration_order(default_integration_order),
  m_cell_data(cell_data),
  m_input_parameters(physics_block_plist),
  m_build_transient_support(build_transient_support),
  m_global_data(global_data)
{
  TEUCHOS_ASSERT(nonnull(physics_block_plist));
  m_physics_id = physics_block_plist->name();

  initialize(m_input_parameters,
	     m_default_integration_order,
	     m_element_block_id,
	     m_cell_data,
	     factory, 
	     build_transient_support);
}

panzer::PhysicsBlock::
PhysicsBlock(const panzer::PhysicsBlock& pb,
	     const panzer::CellData & cell_data,
             const panzer::EquationSetFactory& factory) :
  m_physics_id(pb.m_physics_id),
  m_element_block_id(pb.m_element_block_id),
  m_default_integration_order(pb.m_default_integration_order),
  m_cell_data(cell_data),  // NOT copied from pb
  m_input_parameters(pb.m_input_parameters),
  m_build_transient_support(pb.m_build_transient_support),
  m_global_data(pb.m_global_data)
{
  initialize(m_input_parameters,
	     m_default_integration_order,
	     m_element_block_id,
	     m_cell_data,
	     factory,
	     m_build_transient_support);
}

void panzer::PhysicsBlock::initialize(const Teuchos::RCP<Teuchos::ParameterList>& input_parameters,
				      const int& default_integration_order,
                                      const std::string & element_block_id,
	                              const panzer::CellData & cell_data,
	                              const panzer::EquationSetFactory& factory,
				      const bool build_transient_support)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  
  TEUCHOS_TEST_FOR_EXCEPTION(input_parameters->numParams() < 1, std::runtime_error,
			     "The physics block \"" << input_parameters->name() 
			     << "\" required by element block \"" << element_block_id
			     << "\" does not have any equation sets associated with it."
			     << " Please add at least one equation set to this physics block!");

  m_equation_sets.clear();

  // Loop over equation sets
  typedef ParameterList::ConstIterator pl_iter;
  for (pl_iter eq = input_parameters->begin(); eq != input_parameters->end(); ++eq) {

    TEUCHOS_TEST_FOR_EXCEPTION( !(eq->second.isList()), std::logic_error,
				"All entries in the physics block \"" << m_physics_id 
				<< "\" must be an equation set sublist!" );

    RCP<ParameterList> eq_set_pl = Teuchos::sublist(input_parameters,eq->first,true);

    RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set
      = factory.buildEquationSet(eq_set_pl, default_integration_order, cell_data, m_global_data, build_transient_support);

    // add this equation set in
    m_equation_sets.push_back(eq_set);

    // figure out the DOF names from equation set
    const std::vector<std::string> & eqNames = eq_set->begin()->getDOFNames();
    for(std::size_t j=0;j<eqNames.size();j++)
       m_dof_names.push_back(eqNames[j]);

    // figure out the DOF names from equation set
    const std::vector<StrPureBasisPair> & sbNames = eq_set->begin()->getProvidedDOFs();
    for(std::size_t j=0;j<sbNames.size();j++)
       m_provided_dofs.push_back(sbNames[j]);

    // Get unique list of bases
    for(std::size_t j=0;j<sbNames.size();j++)
      m_bases[sbNames[j].second->name()] = sbNames[j].second;

  }

  // build up field library
  m_field_lib = Teuchos::rcp(new FieldLibrary);
  for(std::vector<StrPureBasisPair>::const_iterator itr=m_provided_dofs.begin();
      itr!=m_provided_dofs.end();++itr) 
     m_field_lib->addFieldAndBasis(itr->first,itr->second);

  // setup element blocks: for each evaluation type
  for(std::size_t eq_i=0;eq_i<m_equation_sets.size();eq_i++) {
     RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set = m_equation_sets[eq_i];
     for(panzer::EquationSet_TemplateManager<panzer::Traits>::iterator itr=eq_set->begin();
         itr!=eq_set->end();++itr) {
        itr->setElementBlockId(element_block_id);
        itr->setFieldLayoutLibrary(*m_field_lib); // this will build marriage between basis
                                                  // and integration rule.
     }
  }
 
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const Teuchos::ParameterList& user_data) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    std::vector<StrBasisPair> providedDOFs;

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {

      if(providedDOFs.size()==0) {
         Teuchos::RCP<IntegrationRule> intRule = eval_type->getIntegrationRule();
         for(std::size_t i=0;i<m_provided_dofs.size();i++) {
            Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
            providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
         }
      }

      eval_type->buildAndRegisterEquationSetEvaluators(fm, providedDOFs, user_data);
    }
  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					       const LinearObjFactory<panzer::Traits> & lof,
					       const Teuchos::ParameterList& user_data) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    std::vector<StrBasisPair> providedDOFs;

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {

      if(providedDOFs.size()==0) {
         Teuchos::RCP<IntegrationRule> intRule = eval_type->getIntegrationRule();
         for(std::size_t i=0;i<m_provided_dofs.size();i++) {
            Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
            providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
         }
      }

      eval_type->buildAndRegisterGatherAndOrientationEvaluators(fm, providedDOFs, lof, user_data);
    }
  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					     const Teuchos::ParameterList& user_data) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    std::vector<StrBasisPair> providedDOFs;

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {

      if(providedDOFs.size()==0) {
         Teuchos::RCP<IntegrationRule> intRule = eval_type->getIntegrationRule();
         for(std::size_t i=0;i<m_provided_dofs.size();i++) {
            Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
            providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
         }
      }

      eval_type->buildAndRegisterDOFProjectionsToIPEvaluators(fm, providedDOFs, user_data);
    }
  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				  const LinearObjFactory<panzer::Traits> & lof,
				  const Teuchos::ParameterList& user_data) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    std::vector<StrBasisPair> providedDOFs;

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {

      if(providedDOFs.size()==0) {
         Teuchos::RCP<IntegrationRule> intRule = eval_type->getIntegrationRule();
         for(std::size_t i=0;i<m_provided_dofs.size();i++) {
            Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
            providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
         }
      }

      eval_type->buildAndRegisterScatterEvaluators(fm, providedDOFs, lof, user_data);
    }
  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				       const Teuchos::ParameterList& models,
				       const Teuchos::ParameterList& user_data) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    std::vector<StrBasisPair> providedDOFs;

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {
      if(providedDOFs.size()==0) {
         Teuchos::RCP<IntegrationRule> intRule = eval_type->getIntegrationRule();
         for(std::size_t i=0;i<m_provided_dofs.size();i++) {
            Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
            providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
         }
      }

      eval_type->buildAndRegisterClosureModelEvaluators(fm, providedDOFs, factory, models, user_data);
    }
  }
}

// *******************************************************************

void panzer::PhysicsBlock::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				   const std::string& model_name,
				   const Teuchos::ParameterList& models,
				   const Teuchos::ParameterList& user_data) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    std::vector<StrBasisPair> providedDOFs;

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {
      if(providedDOFs.size()==0) {
         Teuchos::RCP<IntegrationRule> intRule = eval_type->getIntegrationRule();
         for(std::size_t i=0;i<m_provided_dofs.size();i++) {
            Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
            providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
         }
      }

      eval_type->buildAndRegisterClosureModelEvaluators(fm, providedDOFs, factory, model_name, models, user_data);
    }
  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
					   const std::string& model_name,
					   const Teuchos::ParameterList& models,
					   const panzer::LinearObjFactory<panzer::Traits> & lof,
					   const Teuchos::ParameterList& user_data) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Only use the <Residual> evaluation type, so pass through to type specific call
  this->buildAndRegisterInitialConditionEvaluatorsForType<panzer::Traits::Residual>(fm, factory, model_name, models, lof, user_data);

}


// *******************************************************************
const std::vector<std::string>& panzer::PhysicsBlock::getDOFNames() const
{
  return m_dof_names;
}

// *******************************************************************
const std::vector<panzer::StrPureBasisPair>& panzer::PhysicsBlock::getProvidedDOFs() const
{
  return m_provided_dofs;
}

// *******************************************************************
const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& 
panzer::PhysicsBlock::getBases() const
{
  return m_bases;
}

// *******************************************************************
const shards::CellTopology panzer::PhysicsBlock::getBaseCellTopology() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(m_bases.size() == 0, std::runtime_error,
		     "Cannot return a basis since none exist in this physics block.");
  return m_bases.begin()->second->getIntrepidBasis()->getBaseCellTopology();
}

// *******************************************************************
std::string panzer::PhysicsBlock::physicsBlockID() const
{
  return m_physics_id;
}

// *******************************************************************
std::string panzer::PhysicsBlock::elementBlockID() const
{
  return m_element_block_id;
}

// *******************************************************************
const panzer::CellData & panzer::PhysicsBlock::cellData() const
{
  return m_cell_data;
}

// *******************************************************************
Teuchos::RCP<panzer::PhysicsBlock> panzer::PhysicsBlock::copyWithCellData(const panzer::CellData & cell_data,
                                                                          const panzer::EquationSetFactory & factory) const
{
  return Teuchos::rcp(new panzer::PhysicsBlock(*this,cell_data,factory));
}

// *******************************************************************
Teuchos::RCP<panzer::PhysicsBlock> panzer::PhysicsBlock::copyWithCellData(const panzer::CellData & cell_data) const
{
  // use parrot to avoid needing an equation set factory
  EquationSet_FactoryParrot factory(m_equation_sets);
  return Teuchos::rcp(new panzer::PhysicsBlock(*this,cell_data,factory));
}

// *******************************************************************
Teuchos::RCP<panzer::GlobalData> panzer::PhysicsBlock::globalData() const
{
  return m_global_data;
}

// *******************************************************************
