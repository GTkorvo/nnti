
#include <fei_impl_utils.hpp>
#include <fei_FillableMat.hpp>
#include <fei_CSVec.hpp>
#include <fei_CommUtils.hpp>

#include <fei_iostream.hpp>

#include <fei_unit_impl_utils.hpp>

#include <cmath>
#include <limits>

void test_pack_unpack_FillableMat()
{
  FEI_COUT << "testing fei::impl_utils::pack_FillableMat, unpack_FillableMat...";

  fei::FillableMat fm, fm2;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(0, 1, 0.1);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(1, 2, 1.2);
  fm.putCoef(2, 1, 2.1);
  fm.putCoef(2, 2, 2.2);

  std::vector<int> intdata;
  std::vector<double> doubledata;

  fei::impl_utils::pack_FillableMat(fm, intdata, doubledata); 

  fei::impl_utils::unpack_FillableMat(intdata, doubledata, fm2);

  if (fm.getNumRows() != fm2.getNumRows()) {
    throw std::runtime_error("pack/unpack FillableMat, wrong number of rows");
  }

  if (fei::count_nnz(fm) != fei::count_nnz(fm2)) {
    throw std::runtime_error("pack/unpack FillableMat, wrong number of nonzeros");
  }

  if (fm != fm2) {
    throw std::runtime_error("pack/unpack FillableMat test failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_separateBCEqns()
{
  FEI_COUT << "testing fei::impl_utils::separateBCEqns...";

  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(2, 2, 2.0);
  fm.putCoef(3, 3, 3.0);
  fm.putCoef(4, 4, 4.0);
  fm.putCoef(5, 5, 5.0);

  std::vector<int> bcEqns;
  std::vector<double> bcVals;

  fei::impl_utils::separate_BC_eqns(fm, bcEqns, bcVals);

  if (bcEqns.size() != 5 || bcEqns.size() != bcVals.size()) {
    throw std::runtime_error("separate_BC_eqns test 1 failed.");
  }

  const double eps = std::numeric_limits<double>::epsilon();

  if (bcEqns[1] != 2 || std::abs(bcVals[1] - 2.0) > eps) {
    throw std::runtime_error("separate_BC_eqns test 2 failed.");
  }

  fei::FillableMat fm2;

  fm2.putCoef(1, 1, 1.0);

  fei::impl_utils::separate_BC_eqns(fm2, bcEqns, bcVals);

  if (bcEqns.size() != 6 || bcEqns.size() != bcVals.size()) {
    throw std::runtime_error("separate_BC_eqns test 3 failed.");
  }

  if (bcEqns[2] != 2 || std::abs(bcVals[2] - 2.0) > eps) {
    throw std::runtime_error("separate_BC_eqns test 4 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_create_col_to_row_map()
{
  FEI_COUT << "testing fei::impl_utils::create_col_to_row_map...";

  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(0, 1, 0.1);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(2, 2, 2.2);

  std::multimap<int,int> crmap;

  fei::impl_utils::create_col_to_row_map(fm, crmap);

  if (crmap.size() != 4) {
    FEI_COUT << "ERROR, crmap.size()=="<<crmap.size()<<", expected 4."<<FEI_ENDL;
    throw std::runtime_error("create_col_to_row_map failed 1");
  }

  //next make sure that the col-to-row-map indicates that col 1 appears in 2 rows.
  typedef std::multimap<int,int>::iterator MM_Iter;
  std::pair<MM_Iter,MM_Iter> mm = crmap.equal_range(1);
  int num = 0;
  for(; mm.first!=mm.second; ++mm.first) ++num;

  if (num != 2) {
    FEI_COUT << "ERROR, size of equal_range(1)=="<<num<<", expected 2."<<FEI_ENDL;
    throw std::runtime_error("create_col_to_row_map failed 2");
  }

  FEI_COUT << "ok" << FEI_ENDL;
}

void test_remove_couplings()
{
  FEI_COUT << "testing fei::impl_utils::remove_couplings...";

  fei::FillableMat fm;

  fm.putCoef(2,  0, 0.5);
  fm.putCoef(2, 10, 0.5);
  fm.putCoef(8,  2, 0.5);
  fm.putCoef(8, 10, 0.5);

  int levels = fei::impl_utils::remove_couplings(fm);

  if (levels != 1) {
    throw std::runtime_error("remove_couplings test failed 1.");
  }

  //after remove_couplings, the matrix-row for 8 should have
  //2 column-indices, and they should be 0 and 10. Also, the
  //coefficients should be 0.25 and 0.75.
  fei::FillableVec* matrow = fm.getRow(8);

  if (matrow->size() != 2) {
    throw std::runtime_error("matrow 8 has wrong length");
  }

  fei::CSVec csrow(*matrow);
  std::vector<int>& indices = csrow.indices();
  std::vector<double>& coefs = csrow.coefs();
  if (indices[0] != 0 || indices[1] != 10 ||
      std::abs(coefs[0] -0.25) > 1.e-49 || std::abs(coefs[1] -0.75) > 1.e-49) {
    throw std::runtime_error("matrow 8 has wrong contents after remove_couplings");
  }

  levels = fei::impl_utils::remove_couplings(fm);
  if (levels > 0) {
    throw std::runtime_error("remove_couplings test 2 failed");
  }

  FEI_COUT << "ok" << FEI_ENDL;
}

void test_global_union_mat(MPI_Comm comm)
{
  FEI_COUT << "testing fei::impl_utils::global_union(,mat,)...";

  int numProcs = fei::numProcs(comm);
  int localProc = fei::localProc(comm);

  int numlocalrows = 5;
  int rowlen = 5;

  fei::FillableMat globalmat0;
  fei::FillableMat localmat;
  int row=0;
  for(int p=0; p<numProcs; ++p) {
    for(int i=0; i<numlocalrows; ++i) {
      for(int j=0; j<rowlen; ++j) {
        globalmat0.putCoef(row, j, 1.0);
        if (p == localProc) {
          localmat.putCoef(row, j, 1.0);
        }
      }
      ++row;
    }
  }

  fei::FillableMat globalmat;

  fei::impl_utils::global_union(comm, localmat, globalmat);

  if (globalmat != globalmat0) {
    throw std::runtime_error("globalUnion test (mat) failed");
  }

  FEI_COUT << "ok" << FEI_ENDL;
}

bool test_impl_utils::run(MPI_Comm comm)
{
  test_pack_unpack_FillableMat();

  test_separateBCEqns();

  test_create_col_to_row_map();

  test_remove_couplings();

  test_global_union_mat(comm);

  return true;
}

