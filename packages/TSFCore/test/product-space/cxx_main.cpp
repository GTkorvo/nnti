// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ///////////////////////////////
// cxx_main.cpp

#include "TSFCoreProductVectorSpace.hpp"
#include "TSFCoreSerialVectorSpace.hpp"
#include "TSFCoreTestingTools.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

namespace {

template <class Scalar>
Scalar rel_err( const Scalar &s1, const Scalar &s2 )
{
	return
		Teuchos::ScalarTraits<Scalar>::magnitude( s1 - s2 )
		/
		std::max( Teuchos::ScalarTraits<Scalar>::magnitude(s1), Teuchos::ScalarTraits<Scalar>::magnitude(s1) );
}

}

namespace TSFCore {

///
/** Main test driver
 */
template <class Scalar>
bool run_tests( const int n, const int numBlocks, const bool dumpAll, std::ostream* out )
{

	typedef Teuchos::ScalarTraits<Scalar> ST;

	if(out) *out << "\n*** Entering run_tests<"<<ST::name()<<">) ...\n";

	bool success = true, result;
	Scalar sresult1, sresult2;

	std::vector<Teuchos::RefCountPtr<const VectorSpace<Scalar> > >
		vecSpaces(numBlocks);
	for( int i = 0; i < numBlocks; ++i ) {
		vecSpaces[i] = Teuchos::rcp(new SerialVectorSpace<Scalar>(n));
	}

	if(out) *out
		<< "\nA) Performing basic tests on product vectors with serial contituent vectors ...\n";

	if(out) *out << "\nCreating a product space ps with numBlocks="<<numBlocks<<" and n="<<n<<"vector elements per block ...\n";

	ProductVectorSpace<Scalar> ps(numBlocks,&vecSpaces[0]);

	if(out) *out << "\nps.numBlocks()=";
	result = ps.numBlocks() == numBlocks;
	if(!result) success = false;
	if(out) *out
		<< ps.numBlocks() << " == numBlocks=" << numBlocks
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;

	if(out) *out << "\nps.dim()=";
	result = ps.dim() == n*numBlocks;
	if(!result) success = false;
	if(out) *out
		<< ps.dim() << " == n*numBlocks=" << n*numBlocks
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;

	if(out) *out << "\nps.isCompatible(ps)=";
	result = ps.isCompatible(ps);
	if(!result) success = false;
	if(out) *out
		<< result << " == true"
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;

	if(out) *out << "\nCreating product vectors; pv1, pv2 ...\n";
	Teuchos::RefCountPtr<Vector<Scalar> >
		pv1 = ps.createMember(),
		pv2 = ps.createMember();

	const Scalar
		one   = ST::one(),
		two   = Scalar(2)*one,
		three = Scalar(3)*one;

	if(out) *out << "\nassign(&*pv1,2.0) ...\n";
	assign( &*pv1, two );

	if(out) *out << "\nsum(pv1)=";
	sresult1 = sum(*pv1);
	sresult2 = two*Scalar(ps.dim());
	result = ( ST::magnitude( rel_err( sresult1, sresult2 ) )
						 < ST::magnitude( Scalar(10)*ST::eps() ) );
	if(!result) success = false;
	if(out) *out
		<< sresult1 << " == 2*ps.dim()=" << sresult2
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;
	
	if(out && dumpAll) *out
		<< "\npv1 =\n" << *pv1;
	
	if(out) *out << "\nassign(&*pv2,3.0) ...\n";
	assign( &*pv2, three );

	if(out) *out << "\nsum(pv2)=";
	sresult1 = sum(*pv2);
	sresult2 = three*Scalar(ps.dim());
	result = ( ST::magnitude( rel_err( sresult1, sresult2 ) )
						 < ST::magnitude( Scalar(10)*ST::eps() ) );
	if(!result) success = false;
	if(out) *out
		<< sresult1 << " == 3*ps.dim()=" << sresult2
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;
	
	if(out && dumpAll) *out
		<< "\npv2 =\n" << *pv2;

	if(out) *out << "\nps.scalarProd(*pv1,*pv2)=";
	const Scalar
		scalarProdTarget  = two*three*Scalar(ps.dim()),
		scalarProd        = ps.scalarProd(*pv1,*pv2);
	result = ( ST::magnitude( rel_err( scalarProd, scalarProdTarget ) )
						 < ST::magnitude( Scalar(10)*ST::eps() ) );
	if(!result) success = false;
	if(out) *out
		<< scalarProd << " == 2*3*ps.dim()=" << scalarProdTarget
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;

	if(out) *out
		<< "\nCreating a serial vector space ss with numBlocks*n=" << numBlocks*n << " vector elements ...\n";

	SerialVectorSpace<Scalar> ss(numBlocks*n);

	if(out) *out
		<< "\nB) Test the compatibility of serial vectors and product vectors with serial blocks."
		<< "\n   These tests demonstrate the principle of how all incore vectors are compatible ...\n";

	if(out) *out << "\nCreating serial vectors; sv1, sv2 ...\n";
	Teuchos::RefCountPtr<Vector<Scalar> >
		sv1 = ss.createMember(),
		sv2 = ss.createMember();

	if(out) *out << "\nassign(&sv1,*pv1) ...\n";
	assign( &*sv1, *pv1 );

	if(out) *out << "\nsum(sv1)=";
	sresult1 = sum(*sv1);
	sresult2 = two*Scalar(ps.dim());
	result = ( ST::magnitude( rel_err( sresult1, sresult2 ) )
						 < ST::magnitude( Scalar(10)*ST::eps() ) );
	if(!result) success = false;
	if(out) *out
		<< sresult1 << " == 2*ps.dim()=" << sresult2
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;
	
	if(out && dumpAll) *out
		<< "\nsv1 =\n" << *sv1;

	if(out) *out << "\nassign(&pv2,*sv1) ...\n";
	assign( &*pv2, *sv1 );

	if(out) *out << "\nsum(pv2)=";
	sresult1 = sum(*pv2);
	sresult2 = two*Scalar(ps.dim());
	result = ( ST::magnitude( rel_err( sresult1, sresult2 ) )
						 < ST::magnitude( Scalar(10)*ST::eps() ) );
	if(!result) success = false;
	if(out) *out
		<< sresult1 << " == 2*ps.dim()=" << sresult2
		<< " : " << ( result ? "passed" : "failed" ) << std::endl;
	
	if(out && dumpAll) *out
		<< "\npv2 =\n" << *pv2;

	// ToDo: Finish tests!

  if(out) *out
		<< "\n*** Leaving run_tests<"<<ST::name()<<">) ...\n";

	return success;

}

} // namespace TSFCore

int main( int argc, char* argv[] ) {

	using Teuchos::CommandLineProcessor;

	bool success = true;

	bool verbose = true;

	std::ostream &out = std::cout;

	try {

		//
		// Read options from commandline
		//

		int n         = 4;
		int numBlocks = 2;
		bool dumpAll  = false;

		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		clp.setOption( "n", &n, "Number of elements in each constituent vector." );
		clp.setOption( "num-blocks", &numBlocks, "blocks to create." );
		clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

		//
		// Run the tests
		//

		if( !TSFCore::run_tests<float>(n,numBlocks,dumpAll,verbose?&out:NULL) ) success = false;
		if( !TSFCore::run_tests<double>(n,numBlocks,dumpAll,verbose?&out:NULL) ) success = false;
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
		if( !TSFCore::run_tests<std::complex<float> >(n,numBlocks,dumpAll,verbose?&out:NULL) ) success = false;
		if( !TSFCore::run_tests<std::complex<double> >(n,numBlocks,dumpAll,verbose?&out:NULL) ) success = false;
#endif

	} // end try
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
		success = false;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception!\n";
		success = false;
	}

	if(verbose) {
		if(success)
			out << "\nAll of the tests seem to have run sucessfully!\n";
		else
			out << "\nOh no! at least one of the test failed!\n";	
	}
	
	return success ? 0 : 1;

}
