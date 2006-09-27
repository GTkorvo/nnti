#include "RBGen_Params.h"
#include "RBGen_ConfigDefs.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Utils.hpp"

Teuchos::RefCountPtr<Teuchos::ParameterList> RBGen::createParams( int argc, char* argv[] )
{
  // Create initial empty parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> params = Teuchos::rcp( new Teuchos::ParameterList() );
/*
  //
  // --------------------------------------------------------------
  //  This function scans the input arguments for relevant information
  //  used to specialize the reduced-basis method and its parameters.
  // --------------------------------------------------------------
  //
  //  params->set("Number of Input Arguments", argc );
  //  params->set("Input Argument Vector", argv );
  //
  // --------------------------------------------------------------
  // GET FILE FORMAT INFORMATION
  // --------------------------------------------------------------
  //
  // Absolute path where the input files (snapshot files, format file, steady state)
  //
  std::string in_path = "/home/hkthorn/TestCode/DataSets/Salsa/";
  params->set( "Data Input Path", in_path );
  //
  // Absolute path where the files containing the computed reduced basis should be placed
  //
  std::string out_path = "/home/hkthorn/TestCode/ROM/POD/test_output/";
  params->set( "Data Output Path", out_path );
  //
  // Input / Output file format type
  // Choices:  "Matrix Market", "Burkardt", or "netCDF"
  //
  params->set( "File IO Type", "netCDF" );
  // 
  // Name of nodal file (containing XY coordinates)
  // Note:  Only needed if the file IO type is "Burkardt"
  //
  params->set( "Data Format File", in_path + "xy.txt" );
  //
  // Name of output file 
  // 
  //params->set( "Reduced Basis Output File", out_path + "pod_basis.nc" );
  //
  // --------------------------------------------------------------
  // GENERATE SNAPSHOT FILENAMES
  // --------------------------------------------------------------
  //
  // Vector of snapshot filenames
  //
  // If the file format is "Burkardt", construct all the 500 input filenames.
  //
  if (Teuchos::getParameter<std::string>( *params, "File IO Type" ) == "Burkardt") {
    std::vector< std::string > filenames;
    int num_files = 500;
    for (int i=0; i<num_files; i++) {
      if (i < 9)
	filenames.push_back( in_path + "up00" + Teuchos::Utils::toString(i + 1) + ".txt" );
      else if (i < 99)
	filenames.push_back( in_path + "up0" + Teuchos::Utils::toString(i + 1) + ".txt" );
      else if (i < 999)
	filenames.push_back( in_path + "up" + Teuchos::Utils::toString(i + 1) + ".txt" );
      else
	cout << "There are more than 1000 files!" << endl;
    }
    //
    // Insert the vector of snapshot filenames into the parameter list.
    //
    params->set("Snapshot Files", filenames);
  }
  //
  // If the file format is "netCDF", then input each of the filenames that holds any snapshots
  //
  if (Teuchos::getParameter<std::string>( *params, "File IO Type" ) == "netCDF") {
    std::vector< std::string > filenames;
    //
    // List as many files as necessary that hold the snapshots
    //
    filenames.push_back( in_path + "snap-time1200-2400-stride8.nc" );
    //filenames.push_back( in_path + "inout-41x41-1_3.nc" );
    //
    // Insert the vector of snapshot filenames into the parameter list.
    //
    params->set("Snapshot Files", filenames);
  }
  //
  // --------------------------------------------------------------
  // GET PREPROCESSING INFORMATION
  // --------------------------------------------------------------
  //
  // Preprocessing method used on the input snapshot set
  // Choices:  "none" = no preprocessing 
  //           "ModifiedSS" = creates modified snapshot set using steady state file, scalings, and scaling_idx 
  //
  params->set( "Preprocessing", "none" );
  //
  // Name of steady state file to be used if the preprocessing method is "ModifiedSS"
  //
  params->set( "Steady State File", in_path + "snap-time_ave1200-7200.nc" );
  //
  // Scaling vector for subtracting steady state from solutions given in the snapshot files
  // if the preprocessing method is "ModifiedSS".
  // List as many scalings as needed for all snapshots
  //
  std::vector< double > scalings;
  //scalings.push_back( 1.0/3.0 );
  //scalings.push_back( 5.0/3.0 );
  scalings.push_back( 1.0 );
  //
  params->set("Snapshot Scaling", scalings);
  // 
  // Index vector for subtracting steady state from solutions given in the snapshot files
  // if the preprocessing method is "ModifiedSS".
  // This vector contains pairs indicating the beginning and ending of each scaling section.
  //
  // NOTE:  If the "Snapshot Scaling Indices" are not specified here, and you are using the 
  // netCDF file i/o handlers, then a scaling_idx vector will be created by the netCDF file i/o 
  // handler that is consistent with the snapshot size of each file.  For example, say you have
  // two files "snap1.nc" which has 100 snapshots and "snap2.nc" which has 500 snapshots, then
  // the netCDF file i/o handler will create "scaling_idx" which has the pairs [0,99] and [100,599].
  // This ALSO assumes that you have two entries in the scaling vector above, one for each file.
  //
  //std::vector< std::pair<int,int> > scaling_idx;
  //std::pair< int, int > first_pair( 0, 249 );
  //scaling_idx.push_back( first_pair );
  //std::pair< int, int > second_pair( 250, 499 );
  //scaling_idx.push_back( second_pair );
  //
  //params->set("Snapshot Scaling Indices", scaling_idx);
  //
  // --------------------------------------------------------------
  // GET REDUCED BASIS METHOD AND SIZE INFORMATION
  // --------------------------------------------------------------
  // 
  // Reduced basis method that RBGen should use
  // Choices:  "LAPACK POD" or "Anasazi POD"
  //
  params->set( "Reduced Basis Method", "Anasazi POD" );   // "LAPACK POD" or "Anasazi POD"
  //
  // Size of reduced basis (number of vectors that should be computed)
  //
  params->set( "Basis Size", 64 );  // Any size of basis you'd like
  //
*/
  return params;
}





