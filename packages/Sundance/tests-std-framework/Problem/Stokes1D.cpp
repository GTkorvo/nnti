#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

using SundanceCore::List;
/** 
 * Solves the Stokes equation in 2D
 */

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1.0) < 1.0e-10;}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int nx = 60;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx*np,
                                                    meshType);

      



      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellPredicate leftPointFunc = new PositionalCellPredicate(leftPointTest);
      CellPredicate rightPointFunc = new PositionalCellPredicate(rightPointTest);
      CellFilter left = points.subset(leftPointFunc);
      CellFilter right = points.subset(rightPointFunc);

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr grad = List(dx);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      double beta = 0.02;
      Expr h = new CellDiameterExpr();
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          - p*(dx*vx)
                          - (h*h*beta)*(grad*q)*(grad*p) - q*(dx*ux),
                          quad2);
        
      /* Define the Dirichlet BC */
      Expr uInflow = cos(x);
      Expr bc =  EssentialBC(left, vx*(ux-uInflow) , quad4);


      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, q), 
                         List(ux, p), vecType);


     

      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_Neumann;
      azOptions[AZ_poly_ord] = 4;
      //      azOptions[AZ_subdomain_solve] = AZ_ilu;
      //      azOptions[AZ_graph_fill] = 1;
      azOptions[AZ_max_iter] = 5000;
      azParams[AZ_tol] = 1.0e-12;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new MatlabWriter("Stokes1d.dat");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("p", new ExprFieldWrapper(soln[1]));
      w.write();

      Expr uxErr = soln[0] - 1.0;
      Expr errExpr = Integral(interior, 
                              uxErr*uxErr,
                              new GaussianQuadrature(4));

      FunctionalEvaluator errInt(mesh, errExpr);

      double errorSq = errInt.evaluate();
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      
      double tol = 1.0e-10;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
