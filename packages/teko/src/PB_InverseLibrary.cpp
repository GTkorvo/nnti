#include "PB_InverseLibrary.hpp"

#include <algorithm>

using Teuchos::RCP;
using Teuchos::rcp;

namespace PB {

InverseLibrary::InverseLibrary()
{
   // setup some valid Stratimikos parameters
   /////////////////////////////////////////////

   // set valid solve factory names
   stratValidSolver_.push_back("Belos"); 
   stratValidSolver_.push_back("Amesos"); 
   stratValidSolver_.push_back("AztecOO"); 

   // set valid preconditioner factory name
   stratValidPrecond_.push_back("ML"); 
   stratValidPrecond_.push_back("Ifpack"); 
}

//! add an unspecified inverse to the library
void InverseLibrary::addInverse(const std::string & label,const Teuchos::ParameterList & pl)
{
   // strip out the label
   const std::string type  = pl.get<std::string>("Type");

   // copy the parameter list so we can modify it
   Teuchos::ParameterList settingsList;
   settingsList.set(type,pl);
   settingsList.sublist(type).remove("Type");
   
   // is this a Stratimikos preconditioner or solver
   if(std::find(stratValidPrecond_.begin(),stratValidPrecond_.end(),type)!=stratValidPrecond_.end()) {
      // this is a Stratimikos preconditioner factory
      addStratPrecond(label,type,settingsList);
   }
   else if(std::find(stratValidSolver_.begin(),stratValidSolver_.end(),type)!=stratValidSolver_.end()) {
      // this is a Stratimikos preconditioner factory
      addStratSolver(label,type,settingsList);
   }
   else 
      TEUCHOS_ASSERT(false);
}

//! Add a Stratimikos solver with a label to the library
void InverseLibrary::addStratSolver(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl)
{
   // add some additional parameters onto the list
   RCP<Teuchos::ParameterList> stratList = rcp(new Teuchos::ParameterList());
   stratList->set("Linear Solver Type",type);
   stratList->set("Linear Solver Types",pl);

   stratSolver_[label] = stratList;
}

//! Add a Stratimikos preconditioner with a label to the library
void InverseLibrary::addStratPrecond(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl)
{
   // add some additional parameters onto the list
   RCP<Teuchos::ParameterList> stratList = rcp(new Teuchos::ParameterList());
   stratList->set("Preconditioner Type",type);
   stratList->set("Preconditioner Types",pl);

   stratPrecond_[label] = stratList;
}

//! Get the fully constructed parameter list for a particular label
Teuchos::RCP<const Teuchos::ParameterList> InverseLibrary::getParameterList(const std::string & label) const
{
   std::map<std::string,RCP<const Teuchos::ParameterList> >::const_iterator itr;

   // check preconditioners
   itr = stratPrecond_.find(label);
   if(itr!=stratPrecond_.end()) return itr->second;
 
   // check solvers
   itr = stratSolver_.find(label);
   if(itr!=stratSolver_.end()) return itr->second;

   return Teuchos::null;
}

//! Get the inverse factory associated with a particular label
Teuchos::RCP<const InverseFactory> InverseLibrary::getInverseFactory(const std::string & label) const
{
   std::map<std::string,RCP<const Teuchos::ParameterList> >::const_iterator itr;

   bool isStratSolver=false,isStratPrecond=false;

   // is this a Stratimikos solver?
   itr = stratPrecond_.find(label);
   isStratPrecond = itr!=stratPrecond_.end();

   // is this a Stratimikos preconditioner?
   if(not isStratPrecond) {
      itr = stratSolver_.find(label);
      isStratSolver = itr!=stratSolver_.end();
   }

   // must be a solver or preconditioner
   TEUCHOS_ASSERT(isStratSolver || isStratPrecond);
   
   RCP<const Teuchos::ParameterList> pl = itr->second;

   Stratimikos::DefaultLinearSolverBuilder strat;
   strat.setParameterList(rcp(new Teuchos::ParameterList(*pl)));

   // build inverse factory
   if(isStratPrecond) {
      // try to build a preconditioner factory
      std::string type = pl->get<std::string>("Preconditioner Type");
      RCP<Thyra::PreconditionerFactoryBase<double> > precFact = strat.createPreconditioningStrategy(type);

      // string must map to a preconditioner
      return rcp(new PreconditionerInverseFactory(precFact));
   }
   else if(isStratSolver) {
      // try to build a solver factory
      std::string type = pl->get<std::string>("Linear Solver Type");
      RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solveFact = strat.createLinearSolveStrategy(type);

      // if its around, build a InverseFactory
      return rcp(new SolveInverseFactory(solveFact));
   }
}

/** \brief Build an inverse library from a parameter list.
  * 
  * Build an inverse library from a parameter list. This will
  * contain all the labeled inverses specified.
  *
  * \param[in] pl Parameter list to build the library from
  *
  * \returns A pointer to the inverse library created.
  */
RCP<InverseLibrary> InverseLibrary::buildFromParameterList(const Teuchos::ParameterList & pl,bool useStratDefaults)
{
   // build from Stratimikos or allocate a new inverse library
   RCP<InverseLibrary> invLib;
   if(useStratDefaults)
      invLib = InverseLibrary::buildFromStratimikos();
   else
      invLib = rcp(new InverseLibrary());

   // to convert the void* like entry
   Teuchos::ParameterList * temp;

   // loop over all entries in parameter list
   Teuchos::ParameterList::ConstIterator itr;
   for(itr=pl.begin();itr!=pl.end();++itr) {
      // get current entry
      std::string label             = itr->first;
      Teuchos::ParameterList & list = itr->second.getValue(temp);
      
      // add to library
      invLib->addInverse(label,list);
   }
   
   return invLib;
}

/** \brief Build an inverse library from Stratimikos
  * 
  * Build an inverse library from Stratimkos. The labels
  * will just be the names in Stratimikos.
  *
  * \param[in] strat Stratimikos object to use
  *
  * \returns A pointer to the inverse library created.
  */
Teuchos::RCP<InverseLibrary> InverseLibrary::buildFromStratimikos(const Stratimikos::DefaultLinearSolverBuilder & strat)
{
   RCP<InverseLibrary> invLib = rcp(new InverseLibrary());

   // get default inveres in Stratimikos
   RCP<const Teuchos::ParameterList> pl = strat.getValidParameters();
   Teuchos::ParameterList lst(pl->sublist("Linear Solver Types"));
   Teuchos::ParameterList pft(pl->sublist("Preconditioner Types"));

   Teuchos::ParameterList::ConstIterator itr;
   Teuchos::ParameterList * temp;

   // loop over all entries in solver list
   for(itr=lst.begin();itr!=lst.end();++itr) {
      // get current entry
      std::string label             = itr->first;
      Teuchos::ParameterList & list = itr->second.getValue(temp);
      list.set("Type",label);
      
      // add to library
      invLib->addInverse(label,list);
   }

   // loop over all entries in preconditioner list
   for(itr=pft.begin();itr!=pft.end();++itr) {
      // get current entry
      std::string label             = itr->first;
      Teuchos::ParameterList & list = itr->second.getValue(temp);
      list.set("Type",label);
      
      // add to library
      invLib->addInverse(label,list);
   }

   return invLib;
}

} // end namespace PB
