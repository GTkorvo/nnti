/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

using SundanceCore::List;
/** 
 * Solves the Poisson equation in 2D on the unit disk
 */

int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();



      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Get a mesh */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource meshReader = new ExodusNetCDFMeshReader("../../../tests-std-framework/Problem/disk.ncdf", meshType);
      Mesh mesh = meshReader.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter bdry = new BoundaryCellFilter();


      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      const double pi = 4.0*atan(1.0);
      /* Define the weak form */
      Expr eqn = Integral(interior, (grad*v)*(grad*u)  + v, quad2);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bdry, v*u, quad4);


      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/bicgstab.xml");
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);

      double R = 1.0;
      Expr exactSoln = 0.25*(x*x + y*y - R*R);

      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      Expr du = L2Projector(discSpace, exactSoln-soln).project();

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("PoissonOnDisk");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.addField("error", new ExprFieldWrapper(du));
      w.write();

      
      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2),
                              new GaussianQuadrature(4));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double tol = 1.0e-4;
      Sundance::passFailTest(sqrt(errorSq), tol);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
