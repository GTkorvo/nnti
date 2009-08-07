// Use case for reusing and updating a preconditioner.
//
// We are simply addressing here a way that we would like to be able
// to reuse a preconditioner object. In this use case, we're assuming
// that Thyra will have something similar to the
// TSFPreconditionerFactory idea (whether a factory or a strategy...).

// In this example, we are using a combo of old TSF and Thyra
// notation.  The objects here such as ThyraLinearOp,
// ThyraVectorSpace, etc., are handles a la TSFExtended, not raw Thyra
// objects.  


// set up PreconditionerFactory as usual  
ThyraPreconditionerFactory pFac = new HotNewPreconditioner(...);

// Now set up a ReuseablePreconditionerFactory using
// pFac. ReuseablePreconditionerFactory extends PreconditionerFactory with a
// method that allows one to set the preconditioner in advance so that
// any call to createPrecond(anything) returns the already set
// preconditioner.  It takes a preconditioner factory in its ctor.

ThyraPreconditionerFactory reuseablePrecFactory = 
  new ReuseablePreconditionerFactory(pFac);

/* Now let's set the preconditioner  */

reuseablePrecFactory.updatePrecond(A); 
//This creates the precond and stores result

/* set up the solver  */
ThyraSolver mySolve = new GMRESSolver(reuseablePrecFactory, ...);

Ainv = A.backslashOperator(mySolve);

x = Ainv * b;

/* Update A  */

A = A + U;

/* re-solve without changing the preconditioner  */

xU = Ainv * b;

/* update the preconditioner if necessary  */

if (U.norm() > tol) reuseablePrecFactory.updatePrecond(A); 

/* now re-solve with the updated preconditioner */

newxU = Ainv * b;

