#ifndef MLAPI_MULTILEVEL_H
#define MLAPI_MULTILEVEL_H

#include "ml_common.h"
#include "ml_agg_genP.h"
#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"
#include <vector>

namespace MLAPI {

/*!
\class MultiLevelSA

\brief Black-box multilevel smoothed aggregation preconditioner.

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.

*/

class MultiLevelSA : public BaseOperator, public CompObject, public TimeObject {

public:

  // @{ \name Constructors and destructors

  //! Constructs the hierarchy for given Operator and parameters.
  MultiLevelSA(const Operator FineMatrix, Teuchos::ParameterList& List,
               const bool ConstructNow = true) :
    IsComputed_(false)
  {
    FineMatrix_ = FineMatrix;
    List_ = List;
    if (ConstructNow) Compute();
  }

  //! Destructor.
  virtual ~MultiLevelSA()
  { }

  // @}
  // @{ \name Set and Get methods

  //! Returns a copy of the internally stored domain space.
  const Space GetOperatorDomainSpace() const 
  {
    return(FineMatrix_.GetDomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  const Space GetOperatorRangeSpace() const 
  {
    return(FineMatrix_.GetRangeSpace());
  }

  //! Returns a copy of the internally stored domain space.
  inline const Space GetDomainSpace() const 
  {
    return(FineMatrix_.GetDomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  inline const Space GetRangeSpace() const 
  {
    return(FineMatrix_.GetRangeSpace());
  }

  //! Returns a reference to the restriction operator of level \c i.
  inline const Operator& R(const int i) const
  {
    return(R_[i]);
  }

  //! Returns a reference to the operator of level \c i.
  inline const Operator& A(const int i) const
  {
    return(A_[i]);
  }

  //! Returns a reference to the prolongator operator of level \c i.
  inline const Operator& P(const int i) const
  {
    return(P_[i]);
  }

  //! Returns a reference to the inverse operator of level \c i.
  inline const InverseOperator& S(const int i) const
  {
    return(S_[i]);
  }

  //! Returns the actual number of levels
  inline int GetMaxLevels() const
  {
    return(MaxLevels_);
  }

  //! Returns \c true if the hierarchy has been successfully computed, \c false otherwise.
  inline bool IsComputed() const
  {
    return(IsComputed_);
  }

  // @}
  // @{ \name Mathematical methods

  //! Computes the hierarchy.
  void Compute() 
  {

    ResetTimer();
    StackPush();
    IsComputed_ = false;

    // get parameter from the input list
    int         MaxLevels     = List_.get("max levels", 10);
    double      Damping       = List_.get("aggregation: damping factor", 1.3333);
    string      EigenAnalysis = List_.get("eigen-analysis: type", "Anorm");
    int         MaxCoarseSize = List_.get("coarse: max size", 32);
    MultiVector EmptySpace;
    MultiVector ThisNS        = List_.get("aggregation: null space", EmptySpace);
    int         NumPDEEqns    = List_.get("PDE equations", 1);
    string      SmootherType  = List_.get("smoother: type", "symmetric Gauss-Seidel");
    string      CoarseType    = List_.get("coarse: type", "Amesos-KLU");
    
    // build up the default null space
    if (ThisNS.GetNumVectors() == 0) {
      ThisNS.Reshape(FineMatrix_.GetDomainSpace(),NumPDEEqns);
      ThisNS = 0.0;
      for (int i = 0 ; i < ThisNS.GetMyLength() ; ++i)
        for (int j = 0 ; j < NumPDEEqns ;++j)
          if (i % NumPDEEqns == j)
            ThisNS(i,j) = 1.0;
    }

    MultiVector NextNS;     // contains the next-level null space

    A_.resize(MaxLevels);
    R_.resize(MaxLevels);
    P_.resize(MaxLevels);
    S_.resize(MaxLevels);

    // work on increasing hierarchies only.
    A_[0] = FineMatrix_;

    double LambdaMax;
    Operator A;
    Operator C;
    Operator R;
    Operator P;
    Operator Ptent;
    Operator IminusA;
    InverseOperator S;

    int level;

    for (level = 0 ; level < MaxLevels - 1 ; ++level) {

      // only an alias
      A = A_[level];

      if (level)
        List_.set("PDE equations", ThisNS.GetNumVectors());

      if (GetPrintLevel()) {
      ML_print_line("-", 80);
        cout << "current working level   = " << level << endl;
        cout << "number of global rows   = " << A.GetNumGlobalRows() << endl;
        cout << "number of global nnz    = " << A.GetNumGlobalNonzeros() << endl;
        cout << "threshold               = " << List_.get("aggregation: threshold", 0.0) << endl;
        cout << "number of PDE equations = " << NumPDEEqns << endl;
        cout << "null space dimension    = " << ThisNS.GetNumVectors() << endl;
      }

      // load current level into database
      List_.set("workspace: current level", level);

      GetPtent(A, List_, ThisNS, Ptent, NextNS);
      ThisNS = NextNS;
      
      if (Damping) {

        if (EigenAnalysis == "Anorm")
          LambdaMax = MaxEigAnorm(A,true);
        else if (EigenAnalysis == "cg")
          LambdaMax = MaxEigCG(A,true);
        else if (EigenAnalysis == "power-method")
          LambdaMax = MaxEigPowerMethod(A,true);
        else
          ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

#if 0
        MultiVector Diag = GetDiagonal(A);
        Diag.Reciprocal();
        Diag.Scale(Damping / LambdaMax);
        Operator Dinv = GetDiagonal(Diag);
        Operator DinvA = Dinv * A;
        Operator I = GetIdentity(A.GetDomainSpace(),A.GetRangeSpace());
        Operator IminusA = I - DinvA;
#else
        IminusA = GetJacobiIterationOperator(A,Damping / LambdaMax);
#endif
        P = IminusA * Ptent;
      }
      else {
        P = Ptent;
        LambdaMax = -1.0;
      }

      if (GetPrintLevel()) {
        cout << "omega                   = " << Damping << endl;
        if (LambdaMax != -1.0) {
          cout << "lambda max              = " << LambdaMax << endl;
          cout << "damping factor          = " << Damping / LambdaMax << endl;
        }
        cout << "smoother type           = " << SmootherType << endl;
        cout << "relaxation sweeps       = " << List_.get("smoother: sweeps", 1) << endl;
        cout << "smoother damping        = " << List_.get("smoother: damping factor", 0.67) << endl;
      }

      R = GetTranspose(P);
      C = GetRAP(R,A,P);
      // build smoothers
      S.Reshape(A, SmootherType, List_);

      // put operators and inverse in hierarchy
      R_[level    ] = R;
      P_[level    ] = P;
      A_[level + 1] = C;
      S_[level    ] = S;

      // break if coarse matrix is below specified tolerance
      if (C.GetNumGlobalRows() <= MaxCoarseSize) {
        ++level;
        break;
      }
    }

    // set coarse solver
    S.Reshape(A_[level], CoarseType, List_);
    S_[level] = S;
    MaxLevels_ = level + 1;

    // set the label
    SetLabel("SA, L = " + GetString(MaxLevels_) +
             ", smoother = " + SmootherType);

    if (GetPrintLevel()) {
      ML_print_line("-", 80);
      cout << "final level             = " << level << endl;
      cout << "number of global rows   = " << A_[level].GetNumGlobalRows() << endl;
      cout << "number of global nnz    = " << A_[level].GetNumGlobalNonzeros() << endl;
      cout << "coarse solver           = " << CoarseType << endl;
      cout << "time spent in constr.   = " << GetTime() << " (s)" << endl;
      ML_print_line("-", 80);
    }

    IsComputed_ = true;
    StackPop();
    
    // FIXME: update flops!
    UpdateTime();

  }

  //! Applies the preconditioner to \c b_f, returns the result in \c x_f.
  int Apply(const MultiVector& b_f, MultiVector& x_f) const
  {
    ResetTimer();
    StackPush();

    if (IsComputed() == false)
      ML_THROW("Method Compute() must be called before Apply()", -1);
    SolveMultiLevelSA(b_f,x_f,0);
    UpdateTime();

    StackPop();
    return(0);
  }

  //! Recursively called core of the multi level preconditioner.
  int SolveMultiLevelSA(const MultiVector& b_f,MultiVector& x_f, int level) const 
  {
    if (level == MaxLevels_ - 1) {
      x_f = S(level) * b_f;
      return(0);
    }

    MultiVector r_f(P(level).GetRangeSpace());
    MultiVector r_c(P(level).GetDomainSpace());
    MultiVector z_c(P(level).GetDomainSpace());

    // reset flop counter
    S(level).SetFlops(0.0);
    A(level).SetFlops(0.0);
    R(level).SetFlops(0.0);
    P(level).SetFlops(0.0);
    
    // apply pre-smoother
    x_f = S(level) * b_f;
    // new residual
    r_f = b_f - A(level) * x_f;
    // restrict to coarse
    r_c = R(level) * r_f;
    // solve coarse problem
    SolveMultiLevelSA(r_c,z_c,level + 1);
    // prolongate back and add to solution
    x_f = x_f + P(level) * z_c;
    // apply post-smoother
    S(level).Apply(b_f,x_f); 

    UpdateFlops(2.0 * S(level).GetFlops());
    UpdateFlops(A(level).GetFlops());
    UpdateFlops(R(level).GetFlops());
    UpdateFlops(P(level).GetFlops());
    UpdateFlops(2.0 * x_f.GetGlobalLength());
    
    return(0);
  }

  // @}
  // @{ \name Miscellaneous methods

  //! Prints basic information about \c this preconditioner.
  std::ostream& Print(std::ostream& os, 
                      const bool verbose = true) const
  {
    if (GetMyPID() == 0) {
      os << endl;
      os << "*** MLAPI::MultiLevelSA, label = `" << GetLabel() << "'" << endl;
      os << endl;
      os << "Number of levels = " << GetMaxLevels() << endl;
      os << "Flop count       = " << GetFlops() << endl;
      os << "Cumulative time  = " << GetTime() << endl;
      if (GetTime() != 0.0)
        os << "MFlops rate      = " << 1.0e-6 * GetFlops() / GetTime() << endl;
      else
        os << "MFlops rate      = 0.0" << endl;
      os << endl;
    }
    return(os);
  }

  // @}

private:

  //! Maximum number of levels.
  int MaxLevels_;
  //! Fine-level matrix.
  Operator FineMatrix_;
  //! Contains the hierarchy of operators.
  vector<Operator> A_;
  //! Contains the hierarchy of restriction operators.
  vector<Operator> R_;
  //! Contains the hierarchy of prolongator operators.
  vector<Operator> P_;
  //! Contains the hierarchy of inverse operators.
  vector<InverseOperator> S_;
  //! Contains a copy of the input list.
  Teuchos::ParameterList List_;
  //! \c true if the hierarchy has been successfully computed, \c false otherwise.
  bool IsComputed_;

};

} // namespace MLAPI

#endif
