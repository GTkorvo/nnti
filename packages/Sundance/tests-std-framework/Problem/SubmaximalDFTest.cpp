#include "Sundance.hpp"


bool SubmaximalDF()
{
  int nx = 2;
  VectorType<double> vecType = new EpetraVectorType();

  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, 
    0.0, 1.0, nx, meshType);
  Mesh mesh = mesher.getMesh();

  CellFilter interior = new MaximalCellFilter();
  CellFilter pts = new DimensionalCellFilter(0);
  CellFilter middle = pts.subset(new CoordinateValueCellPredicate(0, 0.5));
    
  BasisFamily basis = new Lagrange(1);
  Expr u = new UnknownFunction(basis, "w");
  Expr v = new TestFunction(basis, "v");

  Expr grad = gradient(2);

  QuadratureFamily quad2 = new GaussianQuadrature(2);

  DiscreteSpace dsRes(mesh, basis, middle, vecType);
  Out::root() << "here's the map" << endl;
  dsRes.map()->print(Out::os());

  Expr uRes = new DiscreteFunction(dsRes, 1.0, "uRes");

  WatchFlag watch("watch");
  watch.setParam("discrete function evaluation", 5);
  watch.setParam("dof map setup", 5);
  watch.setParam("dof map access", 5);
  watch.setParam("eval mediator", 5);
  watch.setParam("integration", 1);
  watch.setParam("evaluation", 1);
  watch.deactivate();
  double f1 = L2Norm(mesh, middle, uRes, quad2, watch);

  Out::root() << "f1 = " << f1 << endl;

  bool test1 = std::fabs(f1 - std::sqrt(nx+1)) < 1.0e-10;
  Out::root() << "test1 pass=" << test1 << endl;

  Expr eqn = Integral(interior, (grad*u)*(grad*v), quad2);
  Expr bc = EssentialBC(middle, v*uRes*(u-uRes), quad2);

  LinearSolver<double> solver = LinearSolverBuilder::createSolver("amesos.xml");
    
  LinearProblem prob(mesh, eqn, bc, v, u, vecType);
  Expr u0 = prob.solve(solver);
  double f2 = L2Norm(mesh, interior, u0-1.0, quad2);
  Out::root() << "f2 = " << f2 << endl;

  bool test2 = std::fabs(f2) < 1.0e-10;
  Out::root() << "test2 pass=" << test2 << endl;

#ifdef CHECK_PROJECTION
  Expr y = new CoordExpr(1);
  Out::root() << "projecting by hand" << endl;

  watch.activate();
  Expr eqn_p = Integral(middle, v*(u-y), quad2, watch);
  Expr bc_p;

  LinearProblem prob_p(mesh, eqn_p, bc_p, v, u, vecType);
  Expr y_p = prob_p.solve(solver);



  Out::root() << "projecting by L2Proj obj" << endl;
  L2Projector proj(dsRes, y);
  Expr yRes = proj.project();

  double f3 = L2Norm(mesh, middle, yRes, quad2, watch);
  Out::root() << "f3 = " << f3 << endl;
  double g3 = 0.0;
  for (int i=0; i<=nx; i++) 
  {
    double z = i/(nx+1.0);
    g3 += z*z;
  }
  g3 = std::sqrt(g3);

  bool test3 = std::fabs(f3-g3) < 1.0e-10;
  Out::root() << "test3 pass=" << test3 << endl;
#endif

  return SundanceGlobal::checkTest(!test1 || !test2, 0.1);
}

