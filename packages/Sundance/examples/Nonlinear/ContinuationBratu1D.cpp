#include "Sundance.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"

/* 
 * Solve the Bratu problem in 1D using fixed-point iteration 
 */

int main(int argc, char** argv)
{
  try
  {
    int nx = 32;
    double convTol = 1.0e-8;
    int nSteps = 5;
    double lambdaMax = 0.5;
    Sundance::setOption("nx", nx, "Number of elements");
    Sundance::setOption("tol", convTol, "Convergence tolerance");
    Sundance::setOption("lambda-max", lambdaMax, 
      "final lambda (parameter in Bratu's equation)");
    Sundance::setOption("nSteps", nSteps, 
      "number of steps in lambda (continuation from lambda=0 to lambda=lambdaMax)");

    Sundance::init(&argc, &argv);

    Out::root() << "Bratu problem with continuation (lambda=[0, " << lambdaMax << "] in " 
                << nSteps << " steps)" << endl;
    Out::root() << "Newton's method with automated linearization" 
                << endl << endl;

    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
    Mesh mesh = mesher.getMesh();

    CellFilter interior = new MaximalCellFilter();
    CellFilter sides = new DimensionalCellFilter(mesh.spatialDim()-1);
    CellFilter left = sides.subset(new CoordinateValueCellPredicate(0, 0.0));
    CellFilter right = sides.subset(new CoordinateValueCellPredicate(0, 1.0));
    
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "w");
    Expr v = new TestFunction(basis, "v");

    Expr grad = gradient(1);

    Expr x = new CoordExpr(0);

    Expr lambda = new Sundance::Parameter(0.0);
    const double pi = 4.0*atan(1.0);
    Expr uExact = sin(pi*x);
    Expr R = pi*pi*uExact - lambda*exp(uExact);

    QuadratureFamily quad4 = new GaussianQuadrature(4);
    QuadratureFamily quad2 = new GaussianQuadrature(2);

    DiscreteSpace discSpace(mesh, basis, vecType);
    Expr uPrev = new DiscreteFunction(discSpace, 0.5);

    Expr eqn 
      = Integral(interior, (grad*v)*(grad*u) - v*lambda*exp(u) - v*R, quad4);

    Expr h = new CellDiameterExpr();
    Expr bc = EssentialBC(left+right, v*u/h, quad2); 

    NonlinearProblem prob(mesh, eqn, bc, v, u, uPrev, vecType);

    NonlinearSolver<double> solver 
      = NonlinearSolverBuilder::createSolver("playa-newton-amesos.xml");

    Expr soln = uPrev;
    for (int n=0; n<nSteps; n++)
    {
      double lambdaVal = n*lambdaMax/(nSteps-1.0);
      /* update the value of the parameter */
      lambda.setParameterValue(lambdaVal);
      Out::root() << "continuation step n=" << n
                  << " of " << nSteps << ", lambda="
                  << lambdaVal << endl;

      SolverState<double> state = prob.solve(solver);
    
      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
        std::runtime_error,
        "Nonlinear solve failed to converge: message=" << state.finalMsg());

      Expr soln = uPrev;
      FieldWriter writer = new DSVWriter("ContinuationBratu-" 
        + Teuchos::toString(n) + ".dat");
      writer.addMesh(mesh);
      writer.addField("soln", new ExprFieldWrapper(soln[0]));
      writer.write();
    }


    double L2Err = L2Norm(mesh, interior, soln-uExact, quad4);
    Out::root() << "L2 Norm of error: " << L2Err << endl;
    
    Sundance::passFailTest(L2Err, 1.5/((double) nx*nx));
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
}

