/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"


#include "TSFNOXSolver.H"

/** 
 * Solves Kepler's equation x = u(x) + e*sin(u(x))
 */



int main(int argc, void** argv)
{
 
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      const double pi = 4.0*atan(1.0);
      MeshSource mesher = new PartitionedLineMesher(0.0, pi, 200*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

      /* Create coordinate function */
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

     
      /* Define the weak form */
      Expr ecc = new Parameter(0.5, "e");
      Expr eqn = Integral(interior, v*(u - ecc*sin(u) - x), quad);
      Expr bc;

      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);
      F.verbosity() = VerbLow;
      NOX::TSF::Group::classVerbosity() = VerbLow;

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/nox.xml");
      ParameterList noxParams = reader.getParameters();

      cerr << "solver params = " << noxParams << endl;

      NOXSolver solver(noxParams, F);

      int numEcc = 10;
      double finalEcc = 0.95;
      double maxErr = 0.0;
      for (int r=1; r<=numEcc; r++)
        {
          double e = r*finalEcc/((double) numEcc);
          ecc.setParameterValue(e);
          cerr << "--------------------------------------------------------- " << endl;
          cerr << " solving at eccentricity = " << ecc << endl;
          cerr << "--------------------------------------------------------- " << endl;
          // Solve the nonlinear system
          NOX::StatusTest::StatusType status = solver.solve();

          int maxTerms = 400;
          Expr exactSoln = x;
          int nUsed = 1;
          for (int n=1; n<=maxTerms; n++)
            {
              double coeff = 2.0/((double) n) * jn(n, n*e);
              exactSoln = exactSoln + coeff*sin(n*x);
              nUsed++;
              if (fabs(coeff) < 1.0e-8) break;
            }
          cerr << "used " << nUsed << " terms in Kapteyn series" << endl;
          Expr exactSolnDisc = L2Projector(discSpace, exactSoln).project();
          Expr errExpr = Integral(interior, 
                                  pow(u0-exactSoln, 2),
                                  new GaussianQuadrature(8));
          double errorSq = evaluateIntegral(mesh, errExpr);
          cerr << "error norm = " << sqrt(errorSq) << endl << endl;
          if (sqrt(errorSq) > maxErr) maxErr = sqrt(errorSq);

          /* Write the field in matlab format */
          FieldWriter w = new MatlabWriter("kepler-e" + Teuchos::toString(e) + ".dat");
          w.addMesh(mesh);
          w.addField("eccentric anomaly", new ExprFieldWrapper(u0[0]));
          w.addField("exact", new ExprFieldWrapper(exactSolnDisc));
          w.write();
        }

      double tol = 1.0e-4;
      Sundance::passFailTest(maxErr, tol);
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
