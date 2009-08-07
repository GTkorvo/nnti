// Use case for an "inverse" operator.
// We would like to be able to attach a solver to a linear operator,
// then have a clean way to using the "inverse" as a linear operator.
//
// In this example, we are using a combo of old TSF and Thyra
// notation.  The objects here such as ThyraLinearOp,
// ThyraVectorSpace, etc., are handles a la TSFExtended, not raw Thyra
// objects.  

// Something similar to this is what we'd like to have available in
// Thyra.
//
// As for naming this operation, to Paul and me (both heavy Matlab users),
// the names backslash and backslashOperator are completely understandable. 
// Apparently TSFExtended already uses some Matlab terminology
// (e.g., dotSlash). So we'll use "backslash" here.
// Unlike Matlab's backslash, though, this one can also do an iterative solve.


// Create some vector spaces.
ThyraVectorSpace rangeSpace = new ThyraVectorSpace(...);
ThyraVectorSpace domainSpace = new ThyraVectorSpace(...);

// Define a linear solver (with or without a preconditioning factory/strategy) 
ThyraLinearSolver mySolver = new ThyraGMRESSolver(myPrecStrategy, options...)

ThyraLinearOp A = new ThyraLinearOp(rangeSpace, domainSpace,...);
ThyraLinearOp B = new ThyraLinearOp(rangeSpace, domainSpace,...);

// Implicit multiplication of linear operators
ThyraLinearOp C = A * B;

// If we want to use backslash the way Matlab does, i.e., as an
// operation on a vector that returns another vector, we would use
x = A.backslash(b,mySolver);

// Or, we can create a linear operator whose apply() performs the
// linear solve (the back slash operation) with a specified solver.
// This type of linear operator is called an InverseOp.

ThyraLinearOp Ainv = A.backslashOperator(mySolver);

// Note that the method "backslashOperator" returns a LinearOp of type
// Inverse.

x = Ainv * b;

// A nice notational feature of this is that a solver can be attached
// to a nonsquare matrix, e.g., a least squares solver.  This is
// consistent with Matlab where backslash automatically reverts to a
// least squares solver in the case of a nonsquare matrix.

// What we would like in the end is to be able to use the resulting
// Ainv as a linear operator, e.g.,
ThyraLinearOp D = B * Ainv * B.transpose();

// Then I could attach another solver to this implicit product
ThyraLinearOp Dinv = D.backslashOperator(anotherSolver);

// The result should be that when I use Dinv,
ThyraVector x = Dinv*b;
// it should do a solve with "anotherSolver". If that solver is a
// Krylov method, at each iteration it will do an apply on D, which
// will in turn apply the transpose of B, then do solve with mySolver
// on A, then apply B.

// Similar things come up in my Schur complement solvers. Generally
// the point of the methods we're exploring is to avoid taking inverses
// of inverses like this, but it could come up, and we could certainly
// want to do something like this for testing and debugging.

// Also, having something that represents a solve (or an inverse
// operator) is needed for the block preconditioners such as:
ThyraVectorSpace blockVecSpace = new ThyraProductSpace(space1, space2);
ThyraLinearOp myPrecOp = new ThryaBlockLinearOp(blockVecSpace,blockVecSpace);
myPrecOp(1,1) = Ainv;
myPrecOp(1,2) = C.transpose();
myPrecOp(2,2) = Dinv;
// Again, transpose returns a TransposeOp whose apply is to apply the
// transpose. It does NOT explicitly form the transpose. Thus, if C
// changes, the LinearOp created by C.transpose() will continue to
// point to C and will thus give the transpose of the updated C.

// To explicitly form the transpose of a matrix, we recommend the use
// of a method like (C.transpose()).formExplicit(vecType)
