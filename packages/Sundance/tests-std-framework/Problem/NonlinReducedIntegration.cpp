#include "Sundance.hpp"
#include "SundanceReducedQuadrature.hpp"
#include "SundanceProblemTesting.hpp"
using Sundance::List;


/* 
 * Solves the Bratu equation in 2D using reduced integration
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;}) 
  CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;}) 
  CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
  CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;}) 



  int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);
    int np = MPIComm::world().getNProc();

    int n = 4;
    bool increaseProbSize = true;
    if ( (np % 4)==0 ) increaseProbSize = false;

    Array<double> h;
    Array<double> errQuad;
    Array<double> errReduced;

    for (int i=0; i<4; i++)
    {
      n *= 2;
      int nx = n;
      int ny = n;

      VectorType<double> vecType = new EpetraVectorType();

      MeshType meshType = new BasicSimplicialMeshType();
      
      int npx = -1;
      int npy = -1;
      PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
      TEUCHOS_TEST_FOR_EXCEPT(npx < 1);
      TEUCHOS_TEST_FOR_EXCEPT(npy < 1);
      TEUCHOS_TEST_FOR_EXCEPT(npx * npy != np);
      if (increaseProbSize)
      {
        nx = nx*npx;
        ny = ny*npy;
      }
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx, 
        0.0,  1.0, ny, npy, meshType);
      Mesh mesh = mesher.getMesh();


      WatchFlag watchMe("watch eqn");
      watchMe.setParam("integration setup", 0);
      watchMe.setParam("integration", 0);
      watchMe.setParam("fill", 0);
      watchMe.setParam("evaluation", 0);
      watchMe.deactivate();

      WatchFlag watchBC("watch BCs");
      watchBC.setParam("integration setup", 0);
      watchBC.setParam("integration", 0);
      watchBC.setParam("fill", 0);
      watchBC.setParam("evaluation", 0);
      watchBC.deactivate();
    


      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
      
      BasisFamily basis = new Lagrange(1);
      Expr u = new UnknownFunction(basis, "u");
      Expr v = new TestFunction(basis, "v");

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      QuadratureFamily quad = new ReducedQuadrature();
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      /* Define the weak form */
      const double pi = 4.0*atan(1.0);

      Expr c = cos(pi*x);
      Expr s = sin(pi*x);
      Expr ch = cosh(y);
      Expr sh = sinh(y);
      Expr s2 = s*s; 
      Expr c2 = c*c;
      Expr sh2 = sh*sh;
      Expr ch2 = ch*ch;
      Expr pi2 = pi*pi;
      Expr uEx = s*ch;
      Expr eu = exp(uEx);
      Expr f = -(ch*eu*(-1 + pi2)*s) + ch2*(c2*eu*pi2 - s2) + eu*s2*sh2;

      Expr eqn = Integral(interior, exp(u)*(grad*u)*(grad*v)
        + v*f + v*u*u, quad, watchMe)
        + Integral(right, v*exp(u)*pi*cosh(y), quad,watchBC);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(left+top, v*(u-uEx), quad, watchBC);

      Expr eqn2 = Integral(interior, exp(u)*(grad*u)*(grad*v)
        + v*f + v*u*u, quad2, watchMe)
        + Integral(right, v*exp(u)*pi*cosh(y), quad2,watchBC);
      /* Define the Dirichlet BC */
      Expr bc2 = EssentialBC(left+top, v*(u-uEx), quad2, watchBC);


      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      Expr soln1 = new DiscreteFunction(discSpace, 0.0, "u0");
      Expr soln2 = new DiscreteFunction(discSpace, 0.0, "u0");
      L2Projector proj(discSpace, uEx);
      Expr uEx0 = proj.project();

      NonlinearProblem nlp(mesh, eqn, bc, v, u, soln1, vecType);
      NonlinearProblem nlp2(mesh, eqn2, bc2, v, u, soln2, vecType);
    
      ParameterXMLFileReader reader("nox-aztec.xml");
      ParameterList noxParams = reader.getParameters();
      NOXSolver solver(noxParams);
      nlp.solve(solver);
      nlp2.solve(solver);

      FieldWriter w = new VTKWriter("NonlinReduced-n" + Teuchos::toString(n));
      w.addMesh(mesh);
      w.addField("soln1", new ExprFieldWrapper(soln1[0]));
      w.addField("soln2", new ExprFieldWrapper(soln2[0]));
      w.addField("exact", new ExprFieldWrapper(uEx0[0]));
      w.write();

      Expr err1 = uEx - soln1;
      Expr errExpr1 = Integral(interior, 
        err1*err1,
        new GaussianQuadrature(4));

      Expr err2 = uEx - soln2;
      Expr errExpr2 = Integral(interior, 
        err2*err2,
        new GaussianQuadrature(4));

      Expr err12 = soln2 - soln1;
      Expr errExpr12 = Integral(interior, 
        err12*err12,
        new GaussianQuadrature(4));
      
      double error1 = ::sqrt(evaluateIntegral(mesh, errExpr1));
      double error2 = ::sqrt(evaluateIntegral(mesh, errExpr2));
      double error12 = ::sqrt(evaluateIntegral(mesh, errExpr12));

      Out::root() << "final result: " << n << " "  << error1 << " " << error2 << " " << error12
                  << endl;
        
      h.append(1.0/((double) n));
      errQuad.append(error2);
      errReduced.append(error1);
    }
    
    double pQuad = fitPower(h, errQuad);
    double pRed = fitPower(h, errReduced);
    Out::root() << "exponent (reduced integration) " << pRed << endl;
    Out::root() << "exponent (full integration) " << pQuad << endl;
    

    Sundance::passFailTest(::fabs(pRed-2.0), 0.1);

  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 

  return Sundance::testStatus();
}


