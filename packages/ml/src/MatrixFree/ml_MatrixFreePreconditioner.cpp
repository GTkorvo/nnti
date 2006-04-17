#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_IFPACK)
#include "ml_MatrixFreePreconditioner.h"
#include "ml_aggregate.h"
#include "ml_epetra_utils.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSVD.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "Ifpack_Chebyshev.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_MapColoring.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "ml_RowMatrix.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_lapack.h"

using namespace Teuchos;
using namespace EpetraExt;

const int ML_MFP_PRESMOOTHER_ONLY = 0;
const int ML_MFP_ADDITIVE = 1;
const int ML_MFP_HYBRID = 2;

const int ML_MFP_NONE = -1;
const int ML_MFP_JACOBI = 0;
const int ML_MFP_BLOCK_JACOBI = 1;
const int ML_MFP_CHEBY = 2;

// ============================================================================ 
ML_Epetra::MatrixFreePreconditioner::
MatrixFreePreconditioner(const Epetra_Operator& Operator,
                         const Epetra_CrsGraph& Graph,
                         Epetra_MultiVector& NullSpace,
                         const Epetra_Vector& PointDiagonal,
                         Teuchos::ParameterList& List) :
  verbose_(false),
  Comm_ML_(0),
  Comm_(Operator.Comm()),
  Label_("ML matrix-free preconditioner"),
  IsComputed_(false),
  PrecType_(ML_MFP_HYBRID),
  SmootherType_(ML_MFP_BLOCK_JACOBI),
  omega_(1.00),
  Operator_(Operator),
  InvPointDiagonal_(0),
  PreSmoother_(0),
  PostSmoother_(0),
  R_(0),
  C_(0),
  C_ML_(0),
  MLP_(0),
  NumPDEEqns_(0),
  NumMyBlockRows_(0),
  Time_(0)
{
  InvPointDiagonal_ = new Epetra_Vector(PointDiagonal);

  List_ = List;

  // ML communicator, here based on MPI_COMM_WORLD
  ML_Comm_Create(&Comm_ML_);

  Time_ = new Epetra_Time(Comm());

  ML_CHK_ERRV(Compute(Graph, NullSpace));
}

// ============================================================================ 
ML_Epetra::MatrixFreePreconditioner::
~MatrixFreePreconditioner()
{
  if (InvPointDiagonal_) delete InvPointDiagonal_;

  if (PreSmoother_) delete PreSmoother_;

  if (PostSmoother_) delete PostSmoother_;

  if (Time_) delete Time_;

  if (R_) delete R_;

  if (MLP_) delete MLP_;

  if (C_) delete C_;

  if (C_ML_) ML_Operator_Destroy(&C_ML_);

  ML_Comm_Destroy(&Comm_ML_);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  ResetStartTime();

  if (!X.Map().SameAs(R_->OperatorDomainMap())) ML_CHK_ERR(-1);
  if (X.NumVectors() != B.NumVectors()) ML_CHK_ERR(-1);

  Epetra_MultiVector B_c(R_->OperatorRangeMap(), B.NumVectors());
  Epetra_MultiVector X_c(R_->OperatorRangeMap(), B.NumVectors());

  if (PrecType_ == ML_MFP_PRESMOOTHER_ONLY)
  {
    ML_CHK_ERR(ApplyPreSmoother(X));
  }
  else if (PrecType_ == ML_MFP_ADDITIVE)
  {
    // ================================= //
    // ADDITIVE TWO-LEVEL PRECONDITIONER //
    // ================================= //

    Epetra_MultiVector tmp(B.Map(), B.NumVectors());

    ML_CHK_ERR(R_->Multiply(false, B, B_c));

    ML_CHK_ERR(MLP_->ApplyInverse(B_c, X_c));

    ML_CHK_ERR(R_->Multiply(true, X_c, tmp));

    // apply smoother with zero starting solution
    ML_CHK_ERR(ApplyPreSmoother(X)); 

    // sum up the two contributions
    ML_CHK_ERR(X.Update(1.0, tmp, 1.0));
  }
  else if (PrecType_ == ML_MFP_HYBRID)
  {
    // =============================== //
    // HYBRID TWO-LEVEL PRECONDITIONER //
    // =============================== //

    Epetra_MultiVector tmp(B.Map(), B.NumVectors());
    Epetra_MultiVector sol(B.Map(), B.NumVectors());
    sol = B;

    // apply pre-smoother
    ML_CHK_ERR(ApplyPreSmoother(sol));

    // new residual
    ML_CHK_ERR(Operator_.Apply(sol, tmp));
    ML_CHK_ERR(tmp.Update(1.0, B, -1.0));

    // restrict to coarse
    ML_CHK_ERR(R_->Multiply(false, tmp, B_c));

    X_c.PutScalar(0.0);
    // solve coarse problem
    ML_CHK_ERR(MLP_->ApplyInverse(B_c, X_c));

    // prolongate back
    ML_CHK_ERR(R_->Multiply(true, X_c, tmp));

    // add to solution, X now has the correction
    ML_CHK_ERR(sol.Update(1.0, tmp, 1.0));

    // apply post-smoother
    /////ML_CHK_ERR(PostSmoother_->ApplyInverse(B, sol));
    ML_CHK_ERR(ApplyPostSmoother(sol, B, tmp));
    X = sol;
  }
  else
    ML_CHK_ERR(-3); // type not recognized

  AddAndResetStartTime("ApplyInverse()");

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
SetUseTranspose(bool UseTranspose)
{
  ML_RETURN(-1);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  ML_RETURN(-1);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Coarsen(ML_Operator*A, ML_Aggregate** MLAggr, ML_Operator** P, 
        ML_Operator** R, ML_Operator** C, int NumPDEEqns, int NullSpaceDim,
        double* NullSpace)
{
  // Aggregate object, with settings
  ML_Aggregate_Create(MLAggr);
  
  string CoarsenType = List_.get("aggregation: type", "Uncoupled");
  double Threshold   = List_.get("aggregation: threshold", 0.0);

  ML_Aggregate_Set_MaxLevels(*MLAggr, 2);
  ML_Aggregate_Set_StartLevel(*MLAggr, 0);
  ML_Aggregate_Set_Threshold(*MLAggr, Threshold);
  (*MLAggr)->cur_level = 0;
  ML_Aggregate_Set_Reuse(*MLAggr);
  (*MLAggr)->keep_agg_information = 1;
  
  *P = ML_Operator_Create(Comm_ML());

  if (CoarsenType == "Uncoupled") 
    (*MLAggr)->coarsen_scheme = ML_AGGR_UNCOUPLED;
  else if (CoarsenType == "METIS")
    (*MLAggr)->coarsen_scheme = ML_AGGR_METIS;
  else 
  {
    ML_CHK_ERR(-1);
  }

  ML_Aggregate_Set_NullSpace(*MLAggr, NumPDEEqns, NullSpaceDim, NullSpace,
                             A->invec_leng);

  int NumAggregates = ML_Aggregate_Coarsen(*MLAggr, A, P, Comm_ML());
  if (NumAggregates == 0)
  {
    cerr << "Found 0 aggregates, perhaps the problem is too small." << endl;
    ML_CHK_ERR(-2);
  }

  *R = ML_Operator_Create(Comm_ML());

  ML_Operator_Transpose_byrow(*P, *R);

  *C = ML_Operator_Create(Comm_ML());

  // FIXME: try to build an Epetra_CrsMatrix directly to save memory
  // Note: I must create an Epetra_CrsMatrix object because I need the graph
  // to use EpetraExt coloring routines.
  ////ML_rap(*R, A, *P, *C, ML_EpetraCRS_MATRIX);
  ML_rap(*R, A, *P, *C, ML_MSR_MATRIX);

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Compute(const Epetra_CrsGraph& Graph, Epetra_MultiVector& NullSpace)
{
  Epetra_Time TotalTime(Comm());

  const int NullSpaceDim = NullSpace.NumVectors();
  // get parameters from the list
  string PrecType = List_.get("prec: type", "hybrid");
  string SmootherType = List_.get("smoother: type", "Jacobi");
  string ColoringType = List_.get("coloring: type", "JONES_PLASSMAN");
  int PolynomialDegree = List_.get("smoother: degree", 3);
  string DiagonalColoringType = List_.get("diagonal coloring: type", "JONES_PLASSMAN");
  int MaximumIterations = List_.get("eigen-analysis: max iters", 10);
  string EigenType_ = List_.get("eigen-analysis: type", "cg");
  int OutputLevel = List_.get("output", 10);
  omega_ = List_.get("smoother: damping", omega_);
  ML_Set_PrintLevel(OutputLevel);
  bool LowMemory = List_.get("low memory", true);

  verbose_ = (MyPID() == 0 && ML_Get_PrintLevel() > 5);

  // ================ //
  // check parameters //
  // ================ //

  if (PrecType == "presmoother only")
    PrecType_ = ML_MFP_PRESMOOTHER_ONLY;
  else if (PrecType == "hybrid")
    PrecType_ = ML_MFP_HYBRID;
  else if (PrecType == "additive")
    PrecType_ = ML_MFP_ADDITIVE;
  else
    ML_CHK_ERR(-3); // not recognized

  if (SmootherType == "none")
    SmootherType_ = ML_MFP_NONE;
  else if (SmootherType == "Jacobi")
    SmootherType_ = ML_MFP_JACOBI;
  else if (SmootherType == "block Jacobi")
    SmootherType_ = ML_MFP_BLOCK_JACOBI;
  else if (SmootherType == "Chebyshev")
    SmootherType_ = ML_MFP_CHEBY;
  else
    ML_CHK_ERR(-4); // not recognized

  // =============================== //
  // basic checkings and some output //
  // =============================== //
  
  int OperatorDomainPoints =  Operator_.OperatorDomainMap().NumGlobalPoints();
  int OperatorRangePoints =  Operator_.OperatorRangeMap().NumGlobalPoints();
  int GraphBlockRows = Graph.NumGlobalBlockRows();
  int GraphNnz = Graph.NumGlobalNonzeros();
  NumPDEEqns_ = OperatorRangePoints / GraphBlockRows;
  NumMyBlockRows_ = Graph.NumMyBlockRows();

  if (OperatorDomainPoints != OperatorRangePoints)
    ML_CHK_ERR(-1); // only square matrices

  if (OperatorRangePoints % NumPDEEqns_ != 0)
    ML_CHK_ERR(-2); // num PDEs seems not constant

  if (verbose_)
  {
    ML_print_line("=",78);
    cout << "*** " << endl;
    cout << "*** ML_Epetra::MatrixFreePreconditioner" << endl;
    cout << "***" << endl;
    cout << "The operator domain map has " << OperatorDomainPoints;
    cout << " points, the operator range map "  << OperatorRangePoints << endl;
    cout << "The graph has " << GraphBlockRows << " rows and " << GraphNnz << " nonzeros" << endl;
    cout << "Processors used in computation = " << Comm().NumProc() << endl;
    cout << "Number of PDE equations        = " << NumPDEEqns_ << endl;
    cout << "Null space dimension           = " << NullSpaceDim << endl;
    cout << "Preconditioner type            = " << PrecType << endl;
    cout << "Smoother type                  = " << SmootherType << endl;
    cout << "Coloring type                  = " << ColoringType << endl;
    cout << "Number of V-cycles for C       = " << List_.sublist("ML list").get("cycle applications", 1) << endl;
  }

  ResetStartTime();

  // ==================================== //
  // compute the inverse of the diagonal, //
  // control that no elements are zero.   //
  // ==================================== //
  
  for (int i = 0; i < InvPointDiagonal_->MyLength(); ++i)
    if ((*InvPointDiagonal_)[i] != 0.0)
      (*InvPointDiagonal_)[i] = 1.0 / (*InvPointDiagonal_)[i];

  // ========================================================= //
  // Setup the smoother. I need to extract the block diagonal  //
  // only if block Jacobi is used. For Chebyshev, I scale with //
  // the point diagonal only. In this latter case, I need to   //
  // compute lambda_max of the scaled operator.                //
  // ========================================================= //
  
  // probes for the block diagonal of the matrix.
  if (SmootherType_ == ML_MFP_JACOBI ||
      SmootherType_ == ML_MFP_NONE)
  {
    // do-nothing here
  }
  else if (SmootherType_ == ML_MFP_BLOCK_JACOBI)
  {
    if (verbose_);
      cout << "Diagonal coloring type         = " << DiagonalColoringType << endl;
    ML_CHK_ERR(GetBlockDiagonal(Graph, DiagonalColoringType));

    AddAndResetStartTime("block diagonal construction", true);
  }
  else if (SmootherType_ == ML_MFP_CHEBY)
  {
    double lambda_min = 0.0;
    double lambda_max = 0.0;
    Teuchos::ParameterList IFPACKList;

    if (EigenType_ == "power-method")
    {
      ML_CHK_ERR(Ifpack_Chebyshev::PowerMethod(Operator_, *InvPointDiagonal_,
                                               MaximumIterations, lambda_max));
    }
    else if(EigenType_ == "cg")
    {
      ML_CHK_ERR(Ifpack_Chebyshev::CG(Operator_, *InvPointDiagonal_,
                                      MaximumIterations, lambda_min, 
                                      lambda_max));
    }
    else
      ML_CHK_ERR(-1); // not recognized

    if (verbose_)
    {
      cout << "Using Chebyshev smoother of degree " << PolynomialDegree << endl;
      cout << "Estimating eigenvalues using " <<  EigenType_ << endl;
      cout << "lambda_min = " << lambda_min << ", ";
      cout << "lambda_max = " << lambda_max << endl;
    }

    IFPACKList.set("chebyshev: min eigenvalue", lambda_min);
    IFPACKList.set("chebyshev: max eigenvalue", lambda_max);
    // FIXME: this allocates a new vector inside
    IFPACKList.set("chebyshev: operator inv diagonal", InvPointDiagonal_);
    IFPACKList.set("chebyshev: degree", PolynomialDegree);

    PreSmoother_ = new Ifpack_Chebyshev((Epetra_Operator*)(&Operator_));
    if (PreSmoother_ == 0) ML_CHK_ERR(-1); // memory error?

    IFPACKList.set("chebyshev: zero starting solution", true);
    ML_CHK_ERR(PreSmoother_->SetParameters(IFPACKList));
    ML_CHK_ERR(PreSmoother_->Initialize());
    ML_CHK_ERR(PreSmoother_->Compute());

    PostSmoother_ = new Ifpack_Chebyshev((Epetra_Operator*)(&Operator_));
    if (PostSmoother_ == 0) ML_CHK_ERR(-1); // memory error?

    IFPACKList.set("chebyshev: zero starting solution", false);
    ML_CHK_ERR(PostSmoother_->SetParameters(IFPACKList));
    ML_CHK_ERR(PostSmoother_->Initialize());
    ML_CHK_ERR(PostSmoother_->Compute());
  }

  // ========================================================= //
  // building P and R for block graph. This is done by working //
  // on the Graph_ object. Support is provided for local       //
  // aggregation schemes only so that all is basically local.  //
  // Then, build the block graph coarse problem.               //
  // ========================================================= //
  
  // ML wrapper for Graph_
  ML_Operator* Graph_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraCrsGraph(const_cast<Epetra_CrsGraph*>(&Graph), Graph_ML);

  ML_Aggregate* BlockAggr_ML = 0;
  ML_Operator* BlockPtent_ML = 0, *BlockRtent_ML = 0,* CoarseGraph_ML = 0;

  if (verbose_) cout << endl;

  ML_CHK_ERR(Coarsen(Graph_ML, &BlockAggr_ML, &BlockPtent_ML, &BlockRtent_ML, 
                     &CoarseGraph_ML));

  if (verbose_) cout << endl;

  Epetra_CrsMatrix* GraphCoarse;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(CoarseGraph_ML, GraphCoarse));

  int NumAggregates = BlockPtent_ML->invec_leng;
  ML_Operator_Destroy(&BlockRtent_ML);
  ML_Operator_Destroy(&CoarseGraph_ML);

  AddAndResetStartTime("construction of block C, R, and P", true);
  if (verbose_) cout << endl;

  // ================================================== //
  // coloring of block graph:                           //
  // - color of block row `i' is given by `ColorMap[i]' //
  // - number of colors is ColorMap.NumColors().        //
  // ================================================== //
  
  ResetStartTime();

  CrsGraph_MapColoring* MapColoringTransform;
  
  if (ColoringType == "JONES_PLASSMAN")
    MapColoringTransform = new CrsGraph_MapColoring (CrsGraph_MapColoring::JONES_PLASSMAN,
                                                     0, false, 0);
  else if (ColoringType == "PSEUDO_PARALLEL")
    MapColoringTransform = new CrsGraph_MapColoring (CrsGraph_MapColoring::PSEUDO_PARALLEL,
                                                     0, false, 0);
  else if (ColoringType == "GREEDY")
    MapColoringTransform = new CrsGraph_MapColoring (CrsGraph_MapColoring::GREEDY,
                                                     0, false, 0);
  else if (ColoringType == "LUBY")
    MapColoringTransform = new CrsGraph_MapColoring (CrsGraph_MapColoring::LUBY,
                                                     0, false, 0);
  else 
    ML_CHK_ERR(-1);

  Epetra_MapColoring* ColorMap = &(*MapColoringTransform)(const_cast<Epetra_CrsGraph&>(GraphCoarse->Graph()));

  // move the information from ColorMap to vector Colors
  const int NumColors = ColorMap->MaxNumColors();
  RefCountPtr<Epetra_IntSerialDenseVector> Colors = rcp(new Epetra_IntSerialDenseVector(GraphCoarse->Graph().NumMyRows()));
  for (int i = 0; i < GraphCoarse->Graph().NumMyRows(); ++i)
    (*Colors)[i] = (*ColorMap)[i];

  delete MapColoringTransform;
  delete ColorMap; ColorMap = 0;
  delete GraphCoarse;

  AddAndResetStartTime("coarse graph coloring", true);
  if (verbose_) cout << endl;

  // get some other information about the aggregates, to be used
  // in the QR factorization of the null space. NodesOfAggregate
  // contains the local ID of block rows contained in each aggregate.

  // FIXME: make it faster
  vector< vector<int> > NodesOfAggregate(NumAggregates);

  for (int i = 0; i < Graph.NumMyBlockRows(); ++i)
  {
    int AID = BlockAggr_ML->aggr_info[0][i];
    NodesOfAggregate[AID].push_back(i);
  }

  int MaxAggrSize = 0;
  for (int i = 0; i < NumAggregates; ++i)
  {
    const int& MySize = NodesOfAggregate[i].size();
    if (MySize > MaxAggrSize) MaxAggrSize = MySize;
  }

  // collect aggregate information, and mark all nodes that are
  // connected with each aggregate. These nodes will have a possible
  // nonzero entry after the matrix-matrix product between the Operator_
  // and the tentative prolongator.

  vector<vector<int> > aggregates(NumAggregates);
  vector<int>::iterator iter;

  for (int i = 0; i < NumAggregates; ++i)
    aggregates[i].reserve(MaxAggrSize);

  for (int i = 0; i < Graph.NumMyBlockRows(); ++i)
  {
    int AID = BlockAggr_ML->aggr_info[0][i];

    int NumEntries;
    int* Indices;

    Graph.ExtractMyRowView(i, NumEntries, Indices);

    for (int k = 0; k < NumEntries; ++k)
    {
      // FIXME: use hash??
      const int& GCID = Graph.ColMap().GID(Indices[k]);

      iter = find(aggregates[AID].begin(), aggregates[AID].end(), GCID);
      if (iter == aggregates[AID].end())
        aggregates[AID].push_back(GCID);
    }
  }
  
  int* BlockNodeList = Graph.ColMap().MyGlobalElements();

  // finally get rid of the ML_Aggregate structure.
  ML_Aggregate_Destroy(&BlockAggr_ML);

  const Epetra_Map& FineMap = Operator_.OperatorDomainMap();
  Epetra_Map CoarseMap(-1, NumAggregates * NullSpaceDim, 0, Comm());
  Epetra_Map* BlockNodeListMap = new Epetra_Map(-1, Graph.ColMap().NumMyElements(),
                                                BlockNodeList, 0, Comm());

  vector<int> NodeList(Graph.ColMap().NumMyElements() * NumPDEEqns_);
  for (int i = 0; i < Graph.ColMap().NumMyElements(); ++i)
    for (int m = 0; m < NumPDEEqns_; ++m)
      NodeList[i * NumPDEEqns_ + m] = BlockNodeList[i] * NumPDEEqns_ + m;
  Epetra_Map* NodeListMap = new Epetra_Map(-1, NodeList.size(), &NodeList[0], 0, Comm());

  AddAndResetStartTime("data structures", true);

  // ====================== //
  // process the null space //
  // ====================== //

  // CHECKME
  Epetra_MultiVector NewNullSpace(CoarseMap, NullSpaceDim);
  NewNullSpace.PutScalar(0.0);

  if (NullSpaceDim == 1)
  {
    double* ns_ptr = NullSpace.Values();

    for (int AID = 0; AID < NumAggregates; ++AID)
    {
      double dtemp = 0.0;
      for (int j = 0; j < NodesOfAggregate[AID].size(); j++)
        for (int m = 0; m < NumPDEEqns_; ++m)
        {
          const int& pos = NodesOfAggregate[AID][j] * NumPDEEqns_ + m;
          dtemp += (ns_ptr[pos] * ns_ptr[pos]);
        }
      dtemp = sqrt(dtemp);

      NewNullSpace[0][AID] = dtemp;

      dtemp = 1.0 / dtemp;

      for (int j = 0; j < NodesOfAggregate[AID].size(); j++)
        for (int m = 0; m < NumPDEEqns_; ++m)
          ns_ptr[NodesOfAggregate[AID][j] * NumPDEEqns_ + m] *= dtemp;
    }
  }
  else
  {
    // FIXME
    vector<double> qr_ptr(MaxAggrSize * NumPDEEqns_ * MaxAggrSize * NumPDEEqns_);
    vector<double> tmp_ptr(MaxAggrSize * NumPDEEqns_ * NullSpaceDim);

    vector<double> work(NullSpaceDim);
    int info;

    for (int AID = 0; AID < NumAggregates; ++AID)
    {
      int MySize = NodesOfAggregate[AID].size();
      int MyFullSize = NodesOfAggregate[AID].size() * NumPDEEqns_;
      int lwork = NullSpaceDim;

      for (int k = 0; k < NullSpaceDim; ++k)
        for (int j = 0; j < MySize; ++j)
          for (int m = 0; m < NumPDEEqns_; ++m)
            qr_ptr[k * MyFullSize + j * NumPDEEqns_ + m] = 
              NullSpace[k][NodesOfAggregate[AID][j] * NumPDEEqns_ + m];

      DGEQRF_F77(&MyFullSize, (int*)&NullSpaceDim, &qr_ptr[0], 
                 &MyFullSize, &tmp_ptr[0], &work[0], &lwork, &info);

      ML_CHK_ERR(info);

      if (work[0] > lwork) work.resize((int) work[0]);

      // the upper triangle of qr_tmp is now R, so copy that into the 
      //  new nullspace

      for (int j = 0; j < NullSpaceDim; j++)
        for (int k = j; k < NullSpaceDim; k++)
          NewNullSpace[k][AID * NullSpaceDim + j] = qr_ptr[j + MyFullSize * k];
		 
      // to get this block of P, need to run qr_tmp through another LAPACK 
      // function:

      DORGQR_F77(&MyFullSize, (int*)&NullSpaceDim, (int*)&NullSpaceDim, 
                 &qr_ptr[0], &MyFullSize, &tmp_ptr[0], &work[0], &lwork, &info);
      ML_CHK_ERR(info); // dgeqtr returned a non-zero

      if (work[0] > lwork) work.resize((int) work[0]);

      // insert the Q block into the null space

      for (int k = 0; k < NullSpaceDim; ++k)
        for (int j = 0; j < MySize; ++j)
          for (int m = 0; m < NumPDEEqns_; ++m)
          {
            int LRID = NodesOfAggregate[AID][j] * NumPDEEqns_ + m;
            double& val = qr_ptr[k * MyFullSize + j * NumPDEEqns_ + m];
            NullSpace[k][LRID] = val;
          }
    }
  }

  AddAndResetStartTime("null space setup", true);

  if (verbose_)
    cout << "# colors on processor " << Comm().MyPID() << " = "
        << NumColors << endl;
  if (verbose_)
    cout << "Maximum # of colors = " << NumColors << endl;

  Epetra_FECrsMatrix* AP = 0;
  AP = new Epetra_FECrsMatrix(Copy, FineMap, NumColors * NullSpaceDim);
  // FIXME: the above declaration is a bit overestimated

  if (!LowMemory)
  {
    // ================================================= //
    // allocate one big chunk of memory, and use View    //             
    // to create Epetra_MultiVectors. Note that          //
    // NumColors * NullSpace can indeed be a quite large //
    // value. To reduce the memory consumption, both     //
    // ColoredAP and ExtColoredAP use the same memory    //
    // array.                                            //
    // ================================================= //
    
    Epetra_MultiVector* ColoredP = new Epetra_MultiVector(FineMap, NumColors * NullSpaceDim);
    vector<double> ColoredAP_ptr(NumColors * NullSpaceDim * NodeListMap->NumMyPoints());
    int ColoredP_LDA = FineMap.NumMyPoints();
    int ColoredAP_LDA = NodeListMap->NumMyPoints();

    ColoredP->PutScalar(0.0);

    for (int i = 0; i < BlockPtent_ML->outvec_leng; ++i)
    {
      int allocated = 1;
      int NumEntries;
      int Indices;
      double Values;
      int ierr = ML_Operator_Getrow(BlockPtent_ML, 1 ,&i, allocated,
                                    &Indices,&Values,&NumEntries);
      if (ierr < 0)
        ML_CHK_ERR(-1);

      assert (NumEntries == 1); // this is the block P
      const int& Color = (*Colors)[Indices] - 1;
      const int& MyLength = ColoredP->MyLength();
      for (int k = 0; k < NumPDEEqns_; ++k)
        for (int j = 0; j < NullSpaceDim; ++j)
          (*ColoredP)[(Color * NullSpaceDim + j)][i * NumPDEEqns_ + k] = 
            NullSpace[j][i * NumPDEEqns_ + k];
    }

    ML_Operator_Destroy(&BlockPtent_ML);

    Epetra_MultiVector ColoredAP(View, Operator_.OperatorRangeMap(), 
                                 &ColoredAP_ptr[0], ColoredAP_LDA, 
                                 NumColors * NullSpaceDim);
    // move ColoredAP into ColoredP. This should not be required.
    // but I prefer to skip strange games with View pointers
    Operator_.Apply(*ColoredP, ColoredAP);
    *ColoredP = ColoredAP;

    // FIXME: only if NumProc > 1
    Epetra_MultiVector ExtColoredAP(View, *NodeListMap, 
                                 &ColoredAP_ptr[0], ColoredAP_LDA, 
                                 NumColors * NullSpaceDim);
    Epetra_Import* Importer;
    Importer = new Epetra_Import(*NodeListMap, Operator_.OperatorRangeMap());
    ExtColoredAP.Import(*ColoredP, *Importer, Insert);
    delete Importer;
    delete ColoredP;

    AddAndResetStartTime("computation of AP", true); 

    // populate the actual AP operator, skip some controls to make it faster

    for (int i = 0; i < NumAggregates; ++i)
    {
      for (int j = 0; j < aggregates[i].size(); ++j)
      {
        int GRID = aggregates[i][j];
        int LRID = BlockNodeListMap->LID(GRID); // this is the block ID
        //assert (LRID != -1);
        int GCID = CoarseMap.GID(i * NullSpaceDim);
        //assert (GCID != -1); 
        int color = (*Colors)[i] - 1;
        for (int k = 0; k < NumPDEEqns_; ++k)
          for (int j = 0; j < NullSpaceDim; ++j)
          {
            double val = ExtColoredAP[color * NullSpaceDim + j][LRID * NumPDEEqns_ + k];
            if (val != 0.0)
            {
              int GRID2 = GRID * NumPDEEqns_ + k;
              int GCID2 = GCID + j;
              AP->InsertGlobalValues(1, &GRID2, 1, &GCID2, &val);
              //if (ierr < 0) ML_CHK_ERR(ierr);
            }
          }
      }
    }
  }
  else
  {
    // =============================================================== //
    // apply the operator one color at-a-time. This requires NumColors //
    // cycles over BlockPtent. However, the memory requirements are    //
    // drastically reduced. As for low-memory == false, both ColoredAP //
    // and ExtColoredAP point to the same memory location.             //
    // =============================================================== //
    
    if (verbose_)
      cout << "Using low-memory computation for AP" << endl;

    Epetra_MultiVector ColoredP(FineMap, NullSpaceDim);
    vector<double> ColoredAP_ptr(NullSpaceDim * NodeListMap->NumMyPoints());
    Epetra_MultiVector ColoredAP(View, Operator_.OperatorRangeMap(), 
                                 &ColoredAP_ptr[0], NodeListMap->NumMyPoints(), 
                                 NullSpaceDim);
    Epetra_MultiVector ExtColoredAP(View, *NodeListMap, 
                                 &ColoredAP_ptr[0], NodeListMap->NumMyPoints(), 
                                 NullSpaceDim);
    Epetra_Import Importer(*NodeListMap, Operator_.OperatorRangeMap());

    for (int ic = 0; ic < NumColors; ++ic)
    {
      ColoredP.PutScalar(0.0);

      for (int i = 0; i < BlockPtent_ML->outvec_leng; ++i)
      {
        int allocated = 1;
        int NumEntries;
        int Indices;
        double Values;
        int ierr = ML_Operator_Getrow(BlockPtent_ML, 1 ,&i, allocated,
                                      &Indices,&Values,&NumEntries);
        if (ierr < 0 ||  // something strange in getrow
            NumEntries != 1) // this is the block P
          ML_CHK_ERR(-1);

        const int& Color = (*Colors)[Indices] - 1;
        if (Color != ic)
          continue; // skip this color for this cycle

        const int& MyLength = ColoredP.MyLength();

        for (int k = 0; k < NumPDEEqns_; ++k)
          for (int j = 0; j < NullSpaceDim; ++j)
            ColoredP[j][i * NumPDEEqns_ + k] = 
              NullSpace[j][i * NumPDEEqns_ + k];
      }

      Operator_.Apply(ColoredP, ColoredAP);
      ColoredP = ColoredAP; // just to be safe

      ExtColoredAP.Import(ColoredP, Importer, Insert);

      // populate the actual AP operator, skip some controls to make it faster

      for (int i = 0; i < NumAggregates; ++i)
      {
        for (int j = 0; j < aggregates[i].size(); ++j)
        {
          int GRID = aggregates[i][j];
          int LRID = BlockNodeListMap->LID(GRID); // this is the block ID
          //assert (LRID != -1);
          int GCID = CoarseMap.GID(i * NullSpaceDim);
          //assert (GCID != -1); 
          int color = (*Colors)[i] - 1;
          if (color != ic)
            continue;

          for (int k = 0; k < NumPDEEqns_; ++k)
            for (int j = 0; j < NullSpaceDim; ++j)
            {
              double val = ExtColoredAP[j][LRID * NumPDEEqns_ + k];
              if (val != 0.0)
              {
                int GRID2 = GRID * NumPDEEqns_ + k;
                int GCID2 = GCID + j;
                AP->InsertGlobalValues(1, &GRID2, 1, &GCID2, &val);
                //if (ierr < 0) ML_CHK_ERR(ierr);
              }
            }
        }
      }
    }

    ML_Operator_Destroy(&BlockPtent_ML);
  }

  aggregates.resize(0);
  delete BlockNodeListMap;
  delete NodeListMap;

  Colors = Teuchos::null;

  AP->GlobalAssemble(false);
  AP->FillComplete(CoarseMap, FineMap);
  AP->OptimizeStorage();

  AddAndResetStartTime("computation of the final AP", true); 

  ML_Operator* AP_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraMatrix(AP, AP_ML);

  // ======== //
  // create R //
  // ======== //
  
  vector<int> REntries(NumAggregates * NullSpaceDim);
  for (int AID = 0; AID < NumAggregates; ++AID)
  {
    for (int m = 0; m < NullSpaceDim; ++m)
      REntries[AID * NullSpaceDim + m] = NodesOfAggregate[AID].size() * NumPDEEqns_;
  }

  R_ = new Epetra_CrsMatrix(Copy, CoarseMap, &REntries[0], true);
  REntries.resize(0);

  for (int AID = 0; AID < NumAggregates; ++AID)
  {
    const int& MySize = NodesOfAggregate[AID].size();

    // FIXME: make it faster
    for (int j = 0; j < MySize; ++j)
      for (int m = 0; m < NumPDEEqns_; ++m)
        for (int k = 0; k < NullSpaceDim; ++k)
        {
          int LCID = NodesOfAggregate[AID][j] * NumPDEEqns_ + m;
          int GCID = FineMap.GID(LCID);
          assert (GCID != -1);

          double& val = NullSpace[k][LCID];

          int GRID = CoarseMap.GID(AID * NullSpaceDim + k);
          assert (GRID != -1); // FIXME
          int ierr = R_->InsertGlobalValues(GRID, 1, &val, &GCID);
          if (ierr < 0)
            ML_CHK_ERR(-1);
        }
  }

  NodesOfAggregate.resize(0);

  R_->FillComplete(FineMap, CoarseMap);
  R_->OptimizeStorage();

  ML_Operator* R_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraMatrix(R_, R_ML);

  AddAndResetStartTime("computation of the R", true); 

  // ======== //
  // Create C //
  // ======== //

  C_ML_ = ML_Operator_Create(Comm_ML());
  ML_2matmult(R_ML, AP_ML, C_ML_, ML_MSR_MATRIX);

  ML_Operator_Destroy(&AP_ML);
  ML_Operator_Destroy(&R_ML);
  delete AP;

  C_ = new ML_Epetra::RowMatrix(C_ML_, &Comm(), false);
  assert (R_->OperatorRangeMap().SameAs(C_->OperatorDomainMap()));

  double SetupTime = TotalTime.ElapsedTime();
  TotalTime.ResetStartTime();

  AddAndResetStartTime("computation of the C", true); 

  if (verbose_)
  {
    cout << "Matrix-free preconditioner built. Now building solver for C..." << endl; 
  }

  Teuchos::ParameterList& sublist = List_.sublist("ML list");
  sublist.set("PDE equations", NullSpaceDim);
  sublist.set("null space: type", "pre-computed");
  sublist.set("null space: dimension", NewNullSpace.NumVectors());
  sublist.set("null space: vectors", NewNullSpace.Values());

  MLP_ = new MultiLevelPreconditioner(*C_, sublist, true);

  IsComputed_ = true;

  assert (MLP_ != 0);

  AddAndResetStartTime("computation of the preconditioner for C", true); 

  if (verbose_)
  {
    cout << endl;
    cout << "Total CPU time for construction (all included) = ";
    cout << TotalCPUTime() << endl;
    ML_print_line("=",78);
  }

  return(0);
}

// ============================================================================ 
double ML_Epetra::MatrixFreePreconditioner::
TotalCPUTime() const
{
  double TotalCPUTime = 0.0;
  map<string, double>::iterator iter2;

  for (iter2 = TimeTable.begin(); iter2 != TimeTable.end(); ++iter2)
  {
    TotalCPUTime += iter2->second;
  }

  return(TotalCPUTime);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
GetBlockDiagonal(const Epetra_CrsGraph& Graph, string DiagonalColoringType)
{
  int OperatorDomainPoints = Operator_.OperatorDomainMap().NumGlobalPoints();
  int OperatorRangePoints =  Operator_.OperatorRangeMap().NumGlobalPoints();
  int GraphBlockRows = Graph.NumGlobalBlockRows();
  int GraphNnz = Graph.NumGlobalNonzeros();

  CrsGraph_MapColoring MapColoringTransform(CrsGraph_MapColoring::JONES_PLASSMAN,
                                            0, true, 0);

  Epetra_MapColoring* ColorMap = &(MapColoringTransform(const_cast<Epetra_CrsGraph&>(Graph)));

  const int NumColors = ColorMap->MaxNumColors();

  Epetra_MultiVector X(Operator_.OperatorDomainMap(), NumPDEEqns_ * NumColors);
  X.PutScalar(0.0);

  for (int i = 0; i < Graph.NumMyBlockRows(); ++i)
  {
    int color = (*ColorMap)[i] - 1;
    for (int j = 0; j < NumPDEEqns_; ++j)
    {
      X[color * NumPDEEqns_ + j][i * NumPDEEqns_ + j] = 1.0;
    }
  }

  Epetra_MultiVector AX(Operator_.OperatorRangeMap(), NumPDEEqns_ * NumColors);

  Operator_.Apply(X, AX);

  InvBlockDiag_.resize(Operator_.OperatorRangeMap().NumMyElements() * NumPDEEqns_);
  
  // extract the diagonals

  Epetra_SerialDenseMatrix V(NumPDEEqns_, NumPDEEqns_);
  Epetra_SerialDenseSVD SVD;
  SVD.SetMatrix(V);

  for (int i = 0; i < Graph.NumMyBlockRows(); ++i)
  {
    int color = (*ColorMap)[i] - 1;
    int offset = i * NumPDEEqns_ * NumPDEEqns_;

    // extract the block
    for (int j = 0; j < NumPDEEqns_; ++j)
    {
      for (int k = 0; k < NumPDEEqns_; ++k)
      {
        V(j, k) = AX[color * NumPDEEqns_ + j][i * NumPDEEqns_ + k];
      }
    }

    // invert the block
    SVD.Invert();
    
    // set the inverted block
    for (int j = 0; j < NumPDEEqns_; ++j)
    {
      for (int k = 0; k < NumPDEEqns_; ++k)
      {
        InvBlockDiag_[offset + j * NumPDEEqns_ + k] = (*SVD.InvertedMatrix())(j, k);
      }
    }
  }

  delete ColorMap;

  /* some possible output for debugging
  Epetra_MultiVector XXX(Copy, Operator_.OperatorRangeMap(), &InvBlockDiag_[0],
                         Operator_.OperatorRangeMap().NumMyElements(), NumPDEEqns_);
  */
  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyInvBlockDiag(const double alpha, Epetra_MultiVector& X,
                  const double beta, const Epetra_MultiVector& B) const
{
  assert (X.NumVectors() == 1);
  int OperatorRangePoints =  Operator_.OperatorRangeMap().NumGlobalPoints();
  int NumPDEEqns2 = NumPDEEqns_ * NumPDEEqns_;

  char trans = 'N';
  int NumVectorsX = X.NumVectors();
  vector<double> tmp(NumPDEEqns_);

  size_t len = sizeof(double) * NumPDEEqns_;
  for (int i = 0; i < NumMyBlockRows_; ++i)
  {
    memcpy(&tmp[0], &(B[0][i * NumPDEEqns_]), len);

    int offset = i * NumPDEEqns2;

    DGEMM_F77(&trans, &trans, (int*)&NumPDEEqns_, &NumVectorsX, (int*)&NumPDEEqns_,
              (double*)&alpha, (double*)&InvBlockDiag_[offset], (int*)&NumPDEEqns_, 
              &tmp[0], (int*)&NumPDEEqns_, (double*)&beta, 
              (double*)&X[0][i * NumPDEEqns_], (int*)&NumPDEEqns_);
  }

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyJacobi(Epetra_MultiVector& X, const double omega) const
{
  ML_CHK_ERR(X.Multiply((double)omega, *InvPointDiagonal_, X, 0.0));

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyJacobi(Epetra_MultiVector& X, const Epetra_MultiVector& B,
            const double omega, Epetra_MultiVector& tmp) const
{
  Operator_.Apply(X, tmp);
  tmp.Update(1.0, B, -1.0);
  ML_CHK_ERR(X.Multiply((double)omega, *InvPointDiagonal_, tmp, 1.0));
  ///ML_CHK_ERR(X.Multiply('T', 'N', (double)omega, *InvPointDiagonal_, tmp, 1.0));

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyBlockJacobi(Epetra_MultiVector& X, const double omega) const
{
  ML_CHK_ERR(ApplyInvBlockDiag(omega, X, 0.0, X));

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyBlockJacobi(Epetra_MultiVector& X, const Epetra_MultiVector& B,
            const double omega, Epetra_MultiVector& tmp) const
{
  Operator_.Apply(X, tmp);
  tmp.Update(1.0, B, -1.0);
  ML_CHK_ERR(ApplyInvBlockDiag(omega, X, 1.0, tmp));
  ///ML_CHK_ERR(X.Multiply('T', 'N', omega, *InvBlockDiag_, tmp, 1.0));

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyPreSmoother(Epetra_MultiVector& X) const
{
  switch (SmootherType_) {
  case ML_MFP_NONE:
    break;
  case ML_MFP_JACOBI:
    ML_CHK_ERR(ApplyJacobi(X, omega_));
    break;
  case ML_MFP_BLOCK_JACOBI:
    ML_CHK_ERR(ApplyBlockJacobi(X, omega_));
    break;
  case ML_MFP_CHEBY:
    PreSmoother_->ApplyInverse(X, X);
    break;
  default:
    ML_CHK_ERR(-3); // not recognized
  }

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyPostSmoother(Epetra_MultiVector& X, const Epetra_MultiVector& Y,
                  Epetra_MultiVector& tmp) const
{
  switch (SmootherType_) {
  case ML_MFP_NONE:
    break;
  case ML_MFP_JACOBI:
    ML_CHK_ERR(ApplyJacobi(X, Y, omega_, tmp));
    break;
  case ML_MFP_BLOCK_JACOBI:
    ML_CHK_ERR(ApplyBlockJacobi(X, Y, omega_, tmp));
    break;
  case ML_MFP_CHEBY:
    ML_CHK_ERR(PostSmoother_->ApplyInverse(Y, X));
    break;
  default:
    ML_CHK_ERR(-3); // not recognized
  }

  return(0);
}

// ============================================================================ 
bool ML_Epetra::MatrixFreePreconditioner::
CheckSPD(const Epetra_Operator& A, const bool UseApply,
         const int NumChecks, 
         const int NumVectors) const
{
  bool res = true;
  vector<double> norm(NumVectors);

  if (!IsComputed())
    return(false);

  if (MyPID() == 0)
    cout << "Checking SPD property of the operator... " << endl;

  Epetra_MultiVector X(A.OperatorDomainMap(), NumVectors);
  Epetra_MultiVector AX(A.OperatorRangeMap(), NumVectors);

  try
  {
    for (int i = 0; i < NumChecks; ++i)
    {
      int ierr;

      if (X.Random()) res = false;
      if (UseApply)
        ierr = A.Apply(X, AX);
      else
        ierr = A.ApplyInverse(X, AX);

      if (ierr < 0)
        throw(-1);

      AX.Dot(X, &norm[0]);

      for (int v = 0; v < NumVectors; ++v)
      {
        cout << norm[v] << endl;
        if (norm[v] <= 0.0)
          throw(-2);
      }
    }
  }
  catch(...)
  {
    res = false;
  }

  if (MyPID() == 0)
  {
    if (res)
      cout << "Passed: all x * A * x are positive." << endl;
    else
      cout << "Failed: some  x * A * x are negative or zero!" << endl;
  }

  return(res);
}

#endif
