#include "Thyra_LOWSFactoryBuilder.hpp"

#ifndef TRILINOS_6

#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_IfpackPreconditionerFactory.hpp"
#include "Thyra_MLPreconditionerFactory.hpp"


using namespace Thyra;
using namespace Teuchos;

RefCountPtr<LinearOpWithSolveFactoryBase<double> >
LOWSFactoryBuilder::createLOWSFactory(const ParameterList& params)
{
  /* check that we have a linear solver parameter list */
 //  TEST_FOR_EXCEPTION(params.name() != "Linear Solver",
//                      runtime_error,
//                      "Expected \"Linear Solver\" as name of parameter list input "
//                      "to createLOWSFactory()");

  
  RefCountPtr<LinearOpWithSolveFactoryBase<double> > rtn;  
  RefCountPtr<PreconditionerFactoryBase<double> > prec;  

  if (params.isSublist("Amesos"))
    {
      RefCountPtr<ParameterList> p = rcp(new ParameterList(params.sublist("Amesos")));
      rtn = rcp(new AmesosLinearOpWithSolveFactory());
      rtn->setParameterList(p);
    }
  else if (params.isSublist("Aztec"))
    {
      RefCountPtr<ParameterList> p = rcp(new ParameterList(params.sublist("Aztec")));
      rtn = rcp(new AztecOOLinearOpWithSolveFactory());
      rtn->setParameterList(p);
    }
  else if (params.isSublist("Belos"))
    {
      RefCountPtr<ParameterList> p = rcp(new ParameterList(params.sublist("Belos")));
      rtn = rcp(new BelosLinearOpWithSolveFactory<double>());
      rtn->setParameterList(p);
    }
  else
    {
      TEST_FOR_EXCEPTION(true, runtime_error, 
                         "solver parameter list did not contain one of [Aztec, Amesos, "
                         "Belos]");
    }

  if (params.isSublist("ML"))
    {
      RefCountPtr<ParameterList> p = rcp(new ParameterList(params.sublist("ML")));
      string probType = getParameter<string>(*p, "Problem Type");
      prec = rcp(new MLPreconditionerFactory(p));
    }
  else if (params.isSublist("Ifpack"))
    {
      RefCountPtr<ParameterList> p = rcp(new ParameterList(params.sublist("Ifpack")));
      prec = rcp(new IfpackPreconditionerFactory());
      prec->setParameterList(p);
    }

  TEST_FOR_EXCEPTION(prec.get() != 0 && !rtn->acceptsPreconditionerFactory(),
                     runtime_error,
                     "Huh? You have provided a preconditioner for a solver that cannot "
                     "accept a preconditioner!");

  if (prec.get() != 0 && rtn->acceptsPreconditionerFactory())
    {
      rtn->setPreconditionerFactory(prec, "precond");
    }
  

  return rtn;
  
}


#endif
