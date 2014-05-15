//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ARISTOS_SQPALGO_H
#define ARISTOS_SQPALGO_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Aristos_Vector.hpp"
#include "Aristos_DataPool.hpp"
#include "Aristos_Objective.hpp"
#include "Aristos_Constraints.hpp"
#include "Aristos_HessVec.hpp"
#include "Aristos_LagMult.hpp"
#include "Aristos_FeasStep.hpp"

/** \class Aristos::SQPAlgo
    \brief Implements a trust-region (TR) sequential quadratic programming (SQP) algorithm for
           equality-constrained optimization.

           Aristos::SQPAlgo solves equality-constrained nonlinear optimization problems of the type:
           \f[ \min f(x) \f]
           \f[ \mbox{subject to }  c(x) = 0 \f]

           where \f$f : \mathcal{X} \rightarrow \mathbb{R}\f$, \f$c : \mathcal{X} \rightarrow \mathcal{Y}\f$,
           and \f$\mathcal{X}\f$ and \f$\mathcal{Y}\f$ are Hilbert spaces. We utilize a trust-region
           SQP approach in which the quadratic subproblems are solved using a composite step strategy.
           The TR quadratic subproblems are of the form
           \f[ \min  \frac{1}{2} \langle H_k s^x_k, s^x_k \rangle_{\mathcal{X}} +
                     \langle \nabla_x \mathcal{L} (x_k, \lambda_k), s^x_k \rangle_{\mathcal{X}} \f]
           subject to
           \f[ c_x(x_k)s^x_k + c(x_k) = 0 \f]
           \f[ \|s^x_k\|_{\mathcal{X}} \le \Delta_k \f]
           where \f$\mathcal{L} (x_k, \lambda_k)\f$ is the Lagrangian functional with multipliers
           \f$\lambda_k\f$, \f$H_k\f$ is a representation of the Hessian of the Lagrangian at \f$x_k\f$,
           and \f$\Delta_k\f$ is the trust-region radius at iteration \f$k\f$.
           The step \f$s_k\f$ consists of a \e quasi-normal step that moves the model toward feasibility and 
           a \e tangential step that moves the model toward optimality while staying in the null space of the
           linearized constraints (at least in absence of inexact linear system solves).
           Based on this approach, the SQP algorithm consists of the following steps:
           \li Lagrange multiplier estimation,
           \li computation of the quasi-normal step,
           \li computation of the tangential step,
           \li decision on whether to accept the computed total step and modify the trust-region radius.
           For convenience, as a "sanity check", we also provide a module that checks first-order objective
           and constraint function derivatives and the Hessian of the Lagrangian computation.
*/


namespace Aristos {

class SQPAlgo {
private:

  Teuchos::RefCountPtr<DataPool>    dat_;
  Teuchos::RefCountPtr<Objective>   obj_;
  Teuchos::RefCountPtr<Constraints> constr_;
  Teuchos::RefCountPtr<HessVec>     hessvec_;
  Teuchos::RefCountPtr<LagMult>     lagmult_;
  Teuchos::RefCountPtr<FeasStep>    feasstep_;

public:

  SQPAlgo(DataPool &dat, Objective &obj, Constraints &constr,
          HessVec &hessvec, LagMult &lagmult, FeasStep &feasstep);

  /** \brief Runs the trust-region SQP algorithm.

      \param x [in,out]    - Initial SQP iterate vector, final solution vector.
      \param c [out]       - The vector of constraint values.
      \param l [out]       - The vector of Lagrange multipliers.
      \param iter [out]    - Total number of SQP iterations.
      \param iflag [out]   - Return code: if 0, the algorithm has converged to desired tolerances
                             within the maximum number of iterations; if 1, the maximum number of
                             iterations has been exceeded; if 2, the trust-region radius has fallen
                             below its allowed minimum.
      \param parlist [in]  - Parameter list.

      \return <tt>true</tt> if the SQP run is successful.

      \par Detailed Description:

      The SQP algorithm consists of the following steps:
      \li Lagrange multiplier estimation,
      \li computation of the quasi-normal step,
      \li computation of the tangential step,
      \li decision on whether to accept the computed total step and modify the trust-region radius.

      Please read the documentation provided for the corresponding member functions.
  */
  bool run(Vector &x, Vector &c, Vector &l, int &iter, int &iflag, Teuchos::ParameterList &parlist);



  /** \brief Computes a Lagrange multiplier estimate.

      \param x [in]        - Current SQP iterate.
      \param g [in]        - Gradient of the objective.
      \param l [out]       - The vector of Lagrange multipliers.
      \param tol [in]      - Tolerance for inexact computations.

      \return <tt>true</tt> if runMultiplierEst is successful.

      \par Detailed Description:

      Calls a user-defined multiplier estimate function which is passed in through the data pool object.
      We might provide a default routine in the future.
      
      A common practice in SQP methods is to approximate the Lagrange multipliers by a least-squares
      estimate based on the KKT stationarity conditions (\ref{eq:eqstationarity}). Thus at each iterate $x_k$,
      we solve a problem related to the least-squares problem
      \f{eqnarray*}
        \min_{\lambda} & & \| \nabla_x f(x_k) + c_x(x_k)^* \lambda \|_{\mathcal{X}}.
      \f}
      
      
      INEXACT VERSION:
      
      It should be noted that there are many ways in which Lagrange multiplier estimation problems can
      be formulated, due to specific problem features and data representations provided by the user. The 
      convergence theory for SQP methods poses rather loose requirements on the Lagrange multipliers. In fact,
      our SQP algorithm merely requires that the sequence of Lagrange multipliers be uniformly bounded.

      This is guaranteed through our choice of the stopping criteria for the iterative linear system solver
      that computes the Lagrange multipliers.
  */
  bool runMultiplierEst(const Vector &x, const Vector &g, Vector &l, double tol);



  /** \brief Computes the quasi-normal step.

      \param x [in]        - Current SQP iterate.
      \param c [in]        - The vector of constraint values.
      \param s [out]       - The computed quasi-normal step vector.
      \param delta [in]    - Scaled trust-region radius.
      \param tol [in]      - Tolerance for inexact computations.

      \return <tt>true</tt> if runQuasiNormalStep is successful.

      \par Detailed Description:

      We ask that the quasi-normal step \f$n_k\f$ approximately solve a problem of the following type:
      
      \f{eqnarray*}
        \min         & & \|c_x(x_k) n + c(x_k)\|^2_{\mathcal{Y}}           \\
        \mbox{s.t.}  & & \|n\|_{\mathcal{X}} \le \zeta \Delta_k
      \f}

      The approximate solution is computed using the Powell dogleg method. The dogleg path is computed
      using the Cauchy point \f$n^{cp}_k = \alpha^{cp} c_x(x_k)^*c(x_k)\f$, where the step length
      \f$\alpha^{cp}\f$ is given by
      
      \f{equation*}
      \begin{cases}
        \frac{\| c_x(x_k)^*c(x_k) \|^2_{\mathcal{X}}}{\| c_x(x_k)c_x(x_k)^* c(x_k) \|^2_{\mathcal{Y}}} & \mbox{if }
        \frac{\| c_x(x_k)^*c(x_k) \|^3_{\mathcal{X}}}{\| c_x(x_k)c_x(x_k)^* c(x_k) \|^2_{\mathcal{Y}}} \le \zeta \Delta_k \\
        \frac{\zeta \Delta_k}{\| c_x(x_k)^* c(x_k) \|_{\mathcal{X}}} & \mbox{otherwise},
      \end{cases}
      \f}
      
      and the feasibility step \f$n^N_k\f$ given, for example, as the minimum norm minimizer of the
      2-norm of the linearized constraints (other choices are possible):

      \f{equation*}
        n^N_k = - c_x(x_k)^* \left( c_x(x_k)c(x_k)^* \right)^{-1} c(x_k).
      \f}

      The intersection of the dogleg path (origin - Cauchy point - feasibility point) with the trust-region
      constraint gives the quasi-normal step.

      INEXACT VERSION:

      We require from the quasi-normal step \f$n_k\f$ to give as much decrease as \f$n^{cp}_k\f$, and

      \f{equation*}
        \|c(x_k)\|^2_{\mathcal{Y}} - \|c_x(x_k)n_k+c(x_k)\|^2_{\mathcal{Y}} \ge
        \sigma_1 \left( \|c(x_k)\|^2_{\mathcal{Y}} - \|c_x(x_k)n^{cp}_k+c(x_k)\|^2_{\mathcal{Y}} \right),
      \f}

      for some \f$0 < \sigma_1 \le 1\f$ independent of \f$k\f$. Additionally, due to possible nonnormality of
      \f$n_k\f$ with respect to the tangent space of the constraints, the following condition must hold for
      \f$\kappa_1>0\f$ independent of \f$k\f$:

      \f{equation*}
        \|n_k\|_{\mathcal{X}} \le \kappa_1 \|c(x_k)\|_{\mathcal{Y}}.
      \f}

      All of the above is guaranteed through our choice of the stopping criteria for an iterative linear
      system solver that computes the feasibility step \f$n^N_k\f$.
  */
  bool runQuasiNormalStep(const Vector &x, const Vector &c, Vector &s, double delta, double tol);



  /** \brief Computes the tangential step.

      \param x [in]        - Current SQP iterate.
      \param g [in]        - Current gradient of the objective vector.
      \param v [in]        - Current quasi-normal step.
      \param l [in]        - Current Lagrange multiplier estimate.
      \param s [out]       - The computed tangential step vector.
      \param delta [in]    - Trust-region radius.
      \param cgitmax [in]  - Maximum number of conjugate gradient iterations (internal use only).
      \param cgtol [in]    - Stopping tolerance for the conjugate gradient method (internal use only).
      \param tol [in]      - Tolerance for inexact computations.
      \param cgiter [out]  - Required number of conjugate gradient iterations.
      \param iflag [out]   - Return flag:
                             0 - normal convergence
                             1 - negative curvature detected
                             2 - trust-region radius active
                             3 - maximum number of conjugate gradient iterations exceeded

      \return <tt>true</tt> upon successful termination (not the same as convergence!) of runTangentialStep.

      \par Detailed Description:
      
      The standard requirements on the tangential step are that it lie in the tangent space of the constraints
      and that it move the quadratic model of the Lagrangian toward optimality. Thus we ask that the
      tangential step \f$t_k\f$ approximately solve a problem of the type

      \f{eqnarray*}
      \min        & &  \frac{1}{2} \langle H_k (t+n_k), t+n_k \rangle_{\mathcal{X}} +
                       \langle \nabla_x \mathcal{L} (x_k, \lambda_k), t+n_k \rangle_{\mathcal{X}}     \\
      \mbox{s.t.} & &  c_x(x_k) t = 0                                          \\
                  & &  \|t+n_k\|_{\mathcal{X}} \le \Delta_k.
      \f}
       
      It is solved using the Steihaug variant of the conjugate gradient (CG) method. We use a projected CG
      method in the full space, rather than explicitly computing the null space of \f$c_x(x_k)\f$.
      This method relies on applications of the orthogonal projection operator
      \f$P_k:\mathcal{X} \rightarrow \mathcal{X}\f$, defined by
      \f{equation*}
        P_k = I - c_x(x_k)^* \left( c_x(x_k) c_x(x_k)^* \right) ^{-1} c_x(x_k),
      \f}
      where \f$I:\mathcal{X} \rightarrow \mathcal{X}\f$ is the identity operator.
      
      INEXACT VERSION:
      
      As in the case of the quasi-normal step computation, the tangential step only needs to satisfy a
      fraction of the Cauchy decrease condition. The model function is in this case related to a
      quadratic model of the Lagrangian functional (we leave out the details).
      
      Convergence in presence of inexactness is guaranteed through our choice of the stopping criteria for an 
      iterative linear system solver that applies the projection operator \f$P_k\f$.
      
  */
  bool runTangentialStep(const Vector &x, const Vector &g, const Vector &v, const Vector &l, Vector &s,
                         double delta, int cgitmax, double cgtol, double tol, int &cgiter, int &iflag);



  /** \brief Computes the tangential step using the inexact CG algorithm with full orthogonalization
             and inexactness control.

      \param x [in]          - Current SQP iterate.
      \param g [in]          - Current gradient of the objective vector.
      \param v [in]          - Current quasi-normal step.
      \param l [in]          - Current Lagrange multiplier estimate.
      \param s [out]         - The computed tangential step vector.
      \param delta [in]      - Trust-region radius.
      \param cgitmax [in]    - Maximum number of conjugate gradient iterations (internal use only).
      \param cgtol [in]      - Stopping tolerance for the conjugate gradient method (internal use only).
      \param fixedtol [in]   - Tolerance for inexact computations.
      \param istolfixed [in] - true, if fixed inexact solver tolerances are desired
                               false, if the SQP algorithm is managing inexactness
      \param fullortho [in]  - true, if full orthogonalization of search directions is desired
                               false, if standard CG 2-vector orthogonalization is sufficient
      \param orthocheck [in] - true, if orthoganality of projected residuals is to be verified
                               false, otherwise
      \param fcdcheck [in]   - true, if the fraction of Cauchy decrease condition is to be verified
                               false, otherwise
      \param cgiter [out]    - Returns the required number of conjugate gradient iterations.
      \param iflag [out]     - Return flag:
                               0 - normal convergence
                               1 - negative curvature detected
                               2 - trust-region radius active
                               3 - maximum number of conjugate gradient iterations exceeded

      \return <tt>true</tt> upon successful termination (not the same as convergence!) of runTangentialStep.

      \par Detailed Description:
      
      The standard requirements on the tangential step are that it lie in the tangent space of the constraints
      and that it move the quadratic model of the Lagrangian toward optimality. Thus we ask that the
      tangential step \f$t_k\f$ approximately solve a problem of the type

      \f{eqnarray*}
      \min        & &  \frac{1}{2} \langle H_k (t+n_k), t+n_k \rangle_{\mathcal{X}} +
                       \langle \nabla_x \mathcal{L} (x_k, \lambda_k), t+n_k \rangle_{\mathcal{X}}     \\
      \mbox{s.t.} & &  c_x(x_k) t = 0                                          \\
                  & &  \|t+n_k\|_{\mathcal{X}} \le \Delta_k.
      \f}
       
      It is solved using the Steihaug variant of the conjugate gradient (CG) method. We use a projected CG
      method in the full space, rather than explicitly computing the null space of \f$c_x(x_k)\f$.
      This method relies on applications of the orthogonal projection operator
      \f$P_k:\mathcal{X} \rightarrow \mathcal{X}\f$, defined by
      \f{equation*}
        P_k = I - c_x(x_k)^* \left( c_x(x_k) c_x(x_k)^* \right) ^{-1} c_x(x_k),
      \f}
      where \f$I:\mathcal{X} \rightarrow \mathcal{X}\f$ is the identity operator.
      
      The tangential step only needs to satisfy a
      fraction of the Cauchy decrease condition. The model function is in this case related to a
      quadratic model of the Lagrangian functional (we leave out the details).
      
      Convergence in presence of inexactness is guaranteed through our choice of the stopping criteria for an 
      iterative linear system solver that applies the projection operator \f$P_k\f$.
      
  */
  bool runTangentialStepInx(const Vector &x, const Vector &g, const Vector &v, const Vector &l, Vector &s,
                                     double delta, int cgitmax, double cgtol, double fixedtol,
                                     bool istolfixed, bool fullortho, bool orthocheck, bool fcdcheck,
                                     int &cgiter, int &iflag);
  


  /** \brief Checks whether the computed step is acceptable, and accordingly updates the trust-region
             radius and the merit function penalty parameter.

      \param x [in,out]     - SQP iterate.
      \param l [in]         - Current Lagrange multiplier estimate.
      \param f [in]         - Current value of the objective.
      \param g [in]         - Current gradient of the objective vector.
      \param c [in]         - The vector of constraint values.
      \param v [in]         - Current quasi-normal step.
      \param w [in]         - Current tangential step.
      \param s [out]        - The computed tangential step vector.
      \param delta [in,out] - Trust-region radius.
      \param rho [in,out]   - Penalty parameter.
      \param tol [in]       - Tolerance for inexact computations.
      \param iaccept [out]  - Integer return code:
                              0 - step rejected
                              1 - step accepted

      \return <tt>true</tt> if runTangentialStep terminates successfully.

      \par Detailed Description:

      Global convergence in trust-region SQP methods is usually ensured through the use of a merit
      function, which helps us determine whether a step is acceptable and whether the trust-region
      radius needs to be modified. We use the augmented Lagrangian merit function
  
      \f{equation*}
      \phi(x,\lambda;\rho) = f(x) + \langle \lambda, c(x) \rangle_{\mathcal{Y}} +
                             \rho \|c(x)\|^2_{\mathcal{Y}} =
                             \mathcal{L}(x,\lambda) + \rho \|c(x)\|^2_{\mathcal{Y}}.
      \f}
      
      Let \f$s_k^x=n_k+t_k\f$ be a trial step where \f$n_k\f$ is the computed quasi-normal step, and \f$t_k\f$
      is the computed tangential step, and let \f$\lambda_{k+1} = \lambda_{k} + \Delta\lambda_k\f$ be an updated
      Lagrange multiplier.
      To measure the improvement in the merit function \f$\phi\f$, we compare the \it actual \it reduction and
      the \it predicted \it reduction in moving from the current iterate \f$x_k\f$ to the trial iterate
      \f$x_k+s_k^x\f$. The actual reduction is defined by
      
      \f{equation*}
        ared(s_k^x; \rho_k) = \phi(x_k,\lambda_k;\rho_k) - \phi(x_k+s_k,\lambda_{k+1};\rho_k),
      \f}
      
      and the predicted reduction is given by

      \f{equation*}
        pred(s_k^x; \rho_k) = \phi(x_k,\lambda_k;\rho_k) - \widetilde{\phi}(x_k,\Delta\lambda_k;\rho_k),
      \f}

      where
      \f{equation*}
               \widetilde{\phi}(x_k,\Delta\lambda_k;\rho_k) =
               \frac{1}{2} \langle H_k s_k^x, s_k^x \rangle_{\mathcal{X}} +
                 \langle \nabla_x \mathcal{L}(x_k,\lambda_k), s_k^x \rangle_{\mathcal{X}} +
               \langle \Delta\lambda_k, c_x(x_k)s_k^x+c(x_k) \rangle_{\mathcal{Y}} +
               \rho_k \|c_x(x_k)s_k^x+c(x_k)\|^2_{\mathcal{Y}} .
      \f}

      The following decision procedure is then used (we leave out the details of the penalty parameter update):
      \f{enumerate}
        \addtocounter{enumi}{-1}
        \item Assume $0 < \eta_1 < \eta_2 < 1$, $0 < \gamma_1 < \gamma_2 < 1$.
        \item Acceptance Test.
              \begin{enumerate}
              \item Update penalty parameter $\rho_k$ (details omitted).
              \item Compute actual reduction $ared(s_k^x; \rho_k)$ and predicted reduction $pred(s_k^x; \rho_k)$,
                    and their ratio $\theta_k=\frac{ared(s_k^x; \rho_k)}{pred(s_k^x; \rho_k)}$.
              \item If $\theta_k \ge \eta_1$, set $x_{k+1} = x_k+s_k$,
                    otherwise set $x_{k+1} = x_k$ and reset $\lambda_{k+1} = \lambda_k$.
              \end{enumerate}
        \item Trust-Region Radius Update. Set
              \[
                \Delta_{k+1}  \in
                                  \begin{cases}
                                    [ \Delta_k, \infty] & \mbox{ if } \; \theta_k \ge \eta_2, \\
                                    [ \gamma_2 \Delta_k, \Delta_k ] & \mbox{ if } \;
                                                               \theta_k \in [\eta_1, \eta_2), \\
                                    [ \gamma_1 \Delta_k, \gamma_2 \Delta_k ] & \mbox{ if } \; \theta_k < \eta_1.
                                  \end{cases}
              \]
      \f}

  */
  bool runAcceptStep(Vector &x, const Vector &l, double f, const Vector &g, const Vector &c,
                     const Vector &v, const Vector &w, double &delta, double &rho, double tol, int &iaccept);



  /** \brief Checks whether the computed step is acceptable, and accordingly updates the trust-region
             radius and the merit function penalty parameter. Should be used whenever the SQP algorithm
             manages inexactness in linear system solves.

      \param x [in,out]      - SQP iterate.
      \param l [in]          - Current Lagrange multiplier estimate.
      \param f [in]          - Current value of the objective.
      \param g [in]          - Current gradient of the objective vector.
      \param c [in]          - The vector of constraint values.
      \param v [in]          - Current quasi-normal step.
      \param w [in]          - Current tangential step.
      \param s [out]         - The computed tangential step vector.
      \param delta [in,out]  - Trust-region radius.
      \param rho [in,out]    - Penalty parameter.
      \param addproj[in]     - if <tt>true</tt>, run an additional null-space projection
      \param projtol[in,out] - Tolerance for inexact null-space projections.
      \param tol [in]        - Tolerance for inexact derivative, etc computations.
      \param iaccept [out]   - Integer return code:
                               0 - step rejected
                               1 - step accepted

      \return <tt>true</tt> if runTangentialStepInx terminates successfully.

      \par Detailed Description:

      Global convergence in trust-region SQP methods is usually ensured through the use of a merit
      function, which helps us determine whether a step is acceptable and whether the trust-region
      radius needs to be modified. We use the augmented Lagrangian merit function
  
      \f{equation*}
      \phi(x,\lambda;\rho) = f(x) + \langle \lambda, c(x) \rangle_{\mathcal{Y}} +
                             \rho \|c(x)\|^2_{\mathcal{Y}} =
                             \mathcal{L}(x,\lambda) + \rho \|c(x)\|^2_{\mathcal{Y}}.
      \f}
      
      Let \f$s_k^x=n_k+t_k\f$ be a trial step where \f$n_k\f$ is the computed quasi-normal step, and \f$t_k\f$
      is the computed tangential step, and let \f$\lambda_{k+1} = \lambda_{k} + \Delta\lambda_k\f$ be an updated
      Lagrange multiplier.
      To measure the improvement in the merit function \f$\phi\f$, we compare the \it actual \it reduction and
      the \it predicted \it reduction in moving from the current iterate \f$x_k\f$ to the trial iterate
      \f$x_k+s_k^x\f$. The actual reduction is defined by
      
      \f{equation*}
        ared(s_k^x; \rho_k) = \phi(x_k,\lambda_k;\rho_k) - \phi(x_k+s_k,\lambda_{k+1};\rho_k),
      \f}
      
      and the predicted reduction is given by

      \f{equation*}
        pred(s_k^x; \rho_k) = \phi(x_k,\lambda_k;\rho_k) - \widetilde{\phi}(x_k,\Delta\lambda_k;\rho_k),
      \f}

      where
      \f{equation*}
               \widetilde{\phi}(x_k,\Delta\lambda_k;\rho_k) =
               \frac{1}{2} \langle H_k s_k^x, s_k^x \rangle_{\mathcal{X}} +
                 \langle \nabla_x \mathcal{L}(x_k,\lambda_k), s_k^x \rangle_{\mathcal{X}} +
               \langle \Delta\lambda_k, c_x(x_k)s_k^x+c(x_k) \rangle_{\mathcal{Y}} +
               \rho_k \|c_x(x_k)s_k^x+c(x_k)\|^2_{\mathcal{Y}} .
      \f}

      The following decision procedure is then used (we leave out the details of the penalty parameter update):
      \f{enumerate}
        \addtocounter{enumi}{-1}
        \item Assume $0 < \eta_1 < \eta_2 < 1$, $0 < \gamma_1 < \gamma_2 < 1$.
        \item Acceptance Test.
              \begin{enumerate}
              \item Update penalty parameter $\rho_k$ (details omitted).
              \item Compute actual reduction $ared(s_k^x; \rho_k)$ and predicted reduction $pred(s_k^x; \rho_k)$,
                    and their ratio $\theta_k=\frac{ared(s_k^x; \rho_k)}{pred(s_k^x; \rho_k)}$.
              \item If $\theta_k \ge \eta_1$, set $x_{k+1} = x_k+s_k$,
                    otherwise set $x_{k+1} = x_k$ and reset $\lambda_{k+1} = \lambda_k$.
              \end{enumerate}
        \item Trust-Region Radius Update. Set
              \[
                \Delta_{k+1}  \in
                                  \begin{cases}
                                    [ \Delta_k, \infty] & \mbox{ if } \; \theta_k \ge \eta_2, \\
                                    [ \gamma_2 \Delta_k, \Delta_k ] & \mbox{ if } \;
                                                               \theta_k \in [\eta_1, \eta_2), \\
                                    [ \gamma_1 \Delta_k, \gamma_2 \Delta_k ] & \mbox{ if } \; \theta_k < \eta_1.
                                  \end{cases}
              \]
      \f}

  */
  bool runAcceptStepInx(Vector &x, const Vector &l, double f, const Vector &g, const Vector &c,
                        const Vector &v, Vector &w, double &delta, double &rho, bool addproj,
                        double &projtol, double tol, int &iaccept);



  /** \brief Derivative checks..

      \param x [in]         - Iterate vector.
      \param l [in]         - Lagrange multiplier vector.
      \param dir [in]       - Direction of differentiation.

      \return <tt>true</tt> if runDerivativeCheck terminates successfully.

      \par Detailed Description:

      Checks validity of the computation of the following quantities:
      \li objective gradient
      \li constraint Jacobian
      \li Hessian of the Lagrangian

      It is best to run the checks once before the SQP run, with random input vectors.
  */
  bool runDerivativeCheck(const Vector &x, const Vector &l, const Vector &dir);
  
};

}

#endif
