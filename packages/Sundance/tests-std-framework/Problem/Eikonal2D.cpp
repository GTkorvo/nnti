/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"
/** 
 * Solves the Eikonal equation in 2D
 * | \nabla u |^2 = 1 + \epsilon \nabla^2 u
 */

bool bdryPointTest(const Point& x) 
{
  static int nBdryPts = 0;

  double r2 = x[0]*x[0]+x[1]*x[1];
  bool rtn = ::fabs(r2 - 100.0) < 1.0e-6;
  return rtn;
}

extern "C" {double besi0_(double* x);}

class BesselI0Func : public UserDefFunctor
{
public: 
  /** */
  BesselI0Func() : UserDefFunctor("I0") {;}

  /** */
  virtual double eval0(const Array<double>& vars) const 
  {
    double* x = const_cast<double*>(&(vars[0]));
    return besi0_(x);
  }

};

Expr BesselI0(const Expr& x)
{
  return  new UserDefOp(x, rcp(new BesselI0Func()));
}

int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      int precision = 3;    // precision when printing vectors
      double epsilon = 1.0;   

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh object to fill in from a file. It will be of type BasisSimplicialMesh. */
      MeshType meshType = new BasicSimplicialMeshType();

      /* read the mesh from the file disk.1*/
      MeshSource meshReader = new TriangleMeshReader("../../../tests-std-framework/Problem/disk.1", meshType);

      Mesh mesh = meshReader.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain.  also define edges and boundaries.*/
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellPredicate bdryPointFunc = new PositionalCellPredicate(bdryPointTest);
      CellFilter bdry = edges.subset(bdryPointFunc);

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create a discrete space, and discretize the function 1.0 on it to represent rhs */
      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      Expr S = new DiscreteFunction(discSpace, 1.0, "S");

      /* Discretize the function 0.0 on the same discrete space as an initial guess*/
      Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");

      /* Create differential operators, gradient, and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = SundanceCore::List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);
     
      /* Define the weak form */
      Expr eqn = Integral(interior, (dx*u)*(dx*u)*v+(dy*u)*(dy*u)*v+epsilon*(grad*v)*(grad*u)-S*v, quad2);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bdry, v*u, quad4);

      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);

      /* Get the initial guess */      
      Vector<double> x0 = F.getInitialGuess();

      /* Create an Aztec solver for solving the linear subproblems*/
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azParams[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-13;

      LinearSolver<double> linSolver = new AztecSolver(azOptions,azParams);
     
      /* Now let's create a NOX solver */

      NOX::TSF::Group grp(x0, F, linSolver, precision);

      grp.verbosity() = VerbSilent;
      
      // Set up the status tests
      NOX::StatusTest::NormF statusTestA(grp, 1.0e-10);
      NOX::StatusTest::MaxIters statusTestB(20);
      NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

      // Create the list of solver parameters
      NOX::Parameter::List solverParameters;

      // Set the solver (this is the default)
      solverParameters.setParameter("Nonlinear Solver", "Line Search Based");

      // Create the line search parameters sublist
      NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");

      // Set the line search method
      lineSearchParameters.setParameter("Method","More'-Thuente");

      // Create the line printing parameters sublist
      NOX::Parameter::List& printingParameters = solverParameters.sublist("Printing");

      // Set the line search method
      printingParameters.setParameter("Output Precision",precision);

      // Create the solver
      NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = solver.solve();

      // Print the answer
      cout << "\n" << "-- Parameter List From Solver --" << "\n";
      solver.getParameterList().print(cout);

      // Get the answer
      grp = solver.getSolutionGroup();

     

    


      /* Check against exact solution */
      double R0 = 10.0;
      double z = R0/epsilon;
      double pi = 4.0*atan(1.0);
      double I0 = besi0_(&z);
      Expr r = sqrt(x*x + y*y);
      Expr exactSoln = epsilon * (log(I0) - log(BesselI0(r/epsilon)));
      Expr exactDisc = L2Projector(discSpace, exactSoln).project();

      //      Evaluator::classVerbosity() = VerbExtreme;
      
      /* this code writes the result to a file so we can visualize using paraview */
      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Eikonal2D");
      w.addMesh(mesh);
      w.addField("numerical soln", new ExprFieldWrapper(u0[0]));
      w.addField("exact soln", new ExprFieldWrapper(exactDisc));
      w.write();

      Expr errExpr = Integral(interior, 
                              pow(u0[0]-exactSoln, 2.0),
                              new GaussianQuadrature(4) );
      double errorSq = evaluateIntegral(mesh, errExpr)/pi/R0/R0;
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double tol = 1.0e-6;
      Sundance::passFailTest(errorSq, tol);
      
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
