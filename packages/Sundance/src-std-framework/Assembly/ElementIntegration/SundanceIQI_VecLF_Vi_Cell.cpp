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

#include "SundanceIQI_VecLF_Vi_Cell.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

IQI_VecLF_Vi_Cell::IQI_VecLF_Vi_Cell( int spatialDim ,
					const CellType & maxCellType ,
					const BasisFamily &testBasis ,
					int component ,
					const QuadratureFamily &quad ,
					const ParameterList& verbParams ):
  ElementIntegralLinearFormCell( spatialDim ,
				 maxCellType ,
				 testBasis ,
				 quad ,
				 verbParams )
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_VecLF_Vi_Cell is not implemented" );
}

void IQI_VecLF_Vi_Cell::evaluate( CellJacobianBatch& JTrans,
				   const double* const coeff,
				   RefCountPtr<Array<double> >& A) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_VecLF_Vi_Cell is not implemented" );
}

