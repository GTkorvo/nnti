//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

/* ======================================================================== */
/*!
 \class ShellOptions

 \brief ShellOptions: a class to manage the input arguments and shell variables.

 With this class, it is easy to handle input line arguments and shell
 varibles. For instance, the user can write
 \verbatim
 $ ./a.out -nx 10 -tol 1e-6 -solver=cg
 \endverbatim
 and then easily retrive the value of \c nx, \c tol, and \c solver.
 
 A simple code using this class is as follows:
 \verbatim
 int main(int argc, char *argv[])
  {

   ShellOptions Args(argc,argv);
   int nx = Args.GetIntOption("-nx", 123);
   int ny = Args.GetIntOption("-ny", 145);
   double tol = Args.GetDoubleOption("-tol", 1e-12);
   string solver = Args.GetIntOption("-solver");

   cout << "nx = " << nx << endl;
   cout << "ny = " << ny << " (default value)" << endl;
   cout << "tol = " << tol << endl;
   cout << "solver = " << solver << endl;

   return 0;
   
 }
 \endverbatim

 Each line option can have a value or not. For options with a value,
 the user can specify this values as follows. Let \c -tolerance be the
 name of the option and \c 1e-12 its value. Both choices are valid:
 - \c -option \c 1e-12 (with one or more spaces)
 - \c -option=1e-12 (an `=' sign and no spaces)

 Options are indentified with one or more dashes (`-'). Each option
 cannot have more than one value.

 Note that the user can specify some values without giving them a name.
 This can be done as follows:
  \verbatim
  $ ./a.out value1 value2 value 3 -nx 10 -tol 1e-6 -solver=cg
  \endverbatim
 Here, \c valueX, (X=1,...,9) is stored in the database entry
 \c ARGV_X. 

 To use this class, the user has to build the database using the \c
 argc,argv input arguments. Then, to retrive the option value, the user
 has to use one of the following functions:
 - GetIntOption
 - GetDoubleOption
 - GetStringOption
 
 If option name is not found in the database, a value of 0, 0.0 or an
 empty string is returned. If needed, the user can also specify a
 default value to return when the option name is not found in the
 database. Method \c HaveOption can be used to query the database for
 an option.

 The user can modify the database as well, using
 - SetOption
 - AddOption

 Finally, the user can retrive the integer, double or string value of a
 shell environmental variable using:
 - GetIntShellVariable
 - GetDoubleShellVariable
 - GetStringShellVariable
 
 \date Albuquerque, 01.-Oct-03
 
 \author Marzio Sala, SNL 9214

*/
/* ------------------------------------------------------------------------ */


/* 
 * OptionDataBase : a class to manage the input arguments
 *
 * example:
 *
 * $ ./a.out -nx 10 -ny 20 -tol 1e-6
 *
 * a simple code using this class is as follows:

 int main(int argc, char *argv[]) {
   ShellOptions Args(argc,argv);
   int nx = Args.GetIntOption("-nx");
   int ny = Args.GetIntOption("-ny");
   double tol = Args.GetDoubleOption("-tol");

   std::cout << nx << " " << ny << " " << tol << std::endl;
  }
  * 
  * Marzio Sala, August 20, 2003
  */

#include "Trilinos_Util.h"

using namespace std;

// ===== ===== =========== //
// basic class definitions //
// ===== ===== =========== //

class ShellOptions {

public:

  // ============ //
  // constructors //
  // ------------ //
  
  ShellOptions( void );
  ShellOptions( int argc, char *argv[] );

  // ========== //
  // destructor //
  // ---------- //
  
   ~ShellOptions() {}

  // === //
  // get //
  // --- //
  
  int    GetIntOption( const string input) ;
  int    GetIntOption( const string input, const int def_value) ;

  double GetDoubleOption( const string input) ;
  double GetDoubleOption( const string input, const double def_value) ;

  string GetStringOption( const string input) ;
  string GetStringOption( const string input, const string def_value ) ;

  // =========== //
  // set options //
  // ----------- //

  bool SetOption(const string input, const string value);
  bool SetOption(const string input, const int value);

  // ============= //
  // query and add //
  // ------------- //

  bool   HaveOption(const string input);
  string GetProgramName(void);
  bool   AddOption( const string input, const string value );

  // =============== //
  // shell variables //
  // --------------- //
  
  int    GetIntShellVariable( const char *str );
  double GetDoubleShellVariable( const char *str );
  string GetCharShellVariable( const char *str );

  // ====== //
  // others //
  // ------ //
  
  void    ShowAll() const;
  void    ShowReallyAll() const;

private:
  
  // map containing the arguments
  map<string,string> OptionDatabase;
  
};

/* ======================================================================== */
/*!
 \brief Initialize the database.

*/
/* ------------------------------------------------------------------------ */

ShellOptions::ShellOptions( void )
{

  OptionDatabase["_PROGRAM_NAME_"] = "_NOT_SET_";
  OptionDatabase["_N_ARGS_"] = "0";

}

/* ======================================================================== */
/*!
 \brief Initialize the database using the options given at the shell line.

*/
/* ------------------------------------------------------------------------ */

ShellOptions::ShellOptions(int argc, char *argv[])
{

  OptionDatabase["_PROGRAM_NAME_"] = argv[0];
  OptionDatabase["_N_ARGS_"] = argc;

  // first, manage possible arguments without specifier
  // (e.g., a.out 12 -nx 10 -ny 20)
  // Those arguments are called _ARGV_1_, ... , _ARGV_N_
  // and the value of N is given by _N_UNNAMED_ARGS_

  int N_args = 0;
  char str[80];
  int i=1;
  
  for( i=1 ; i<argc ; ++i ) {
    if( *(argv[i]) == '-' ) break;
    N_args++;
    sprintf( str, "ARGV_%d", N_args );
    OptionDatabase[str] = argv[i];
  }

  OptionDatabase["_N_UNNAMED_ARGS_"] = N_args;
  
  // now only arguments with a dash (possibly followed by one
  // other specifier)
  
  for( ; i<argc ; ++i ) {
    // check if the option has a `=' inside.
    // If so, split the string into two substrings
    char * pos = strchr( argv[i], '='); 
    if( pos != NULL ) {
      *pos = '\0';
      OptionDatabase[argv[i]] = pos+1;
    } else if( i<argc-1 ) {
      if( *(argv[i+1]) != '-' ) {
	OptionDatabase[argv[i]] = argv[i+1];
	++i;
      } else {
	OptionDatabase[argv[i]] = "";
      }
    }
    
  }
}

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as an integer

 This method returns the integer value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns 0.

*/
/* ------------------------------------------------------------------------ */

int ShellOptions::GetIntOption( const string input )
{

  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atoi(OptionDatabase[input].c_str()) );
  }
  return 0;
  
} /* GetIntOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as an integer

 This method returns the integer value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns the specified default value.

*/
/* ------------------------------------------------------------------------ */

int ShellOptions::GetIntOption( const string input, const int def_value)
{

  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atoi(OptionDatabase[input].c_str()) );
  }
  
  return def_value;
   
} /* GetIntOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a double.

 This method returns the double value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns 0.0.

*/
/* ------------------------------------------------------------------------ */

double ShellOptions::GetDoubleOption( const string input)
{

  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atof(OptionDatabase[input].c_str()) );
  }
  
  return 0.0;

} /* GetDoubleOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a double.

 This method returns the double value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns the specified default value.

*/
/* ------------------------------------------------------------------------ */

double ShellOptions::GetDoubleOption( const string input, const double def_value)
{

  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atof(OptionDatabase[input].c_str()) );
  }
  
  return def_value;

} /* GetDoubleOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a C++ string.

 This method returns the string value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns an empty string ("").

*/
/* ------------------------------------------------------------------------ */

string ShellOptions::GetStringOption( const string input)
{

  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( OptionDatabase[input] );
  }
  
  return "";
  
} /* GetStringOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a C++ string.

 This method returns the string value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns the default value.

*/
/* ------------------------------------------------------------------------ */

string ShellOptions::GetStringOption( const string input, const string def_value)
{

  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( OptionDatabase[input] );
  }
  
  return def_value;
  
} /* GetStringOption */

/* ======================================================================== */
/*!
 \brief Check wheter an option is in the database or not

 This method checks whether option \c input is in the databse or not.
 It returns \c true if it is, \c false otherwise.

*/
/* ------------------------------------------------------------------------ */

bool ShellOptions::HaveOption( const string input)
{

  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return true;
  }
  return false;
  
} /* HaveOption */

/* ======================================================================== */
/*!
 \brief Show all the databse entries

*/
/* ------------------------------------------------------------------------ */

void ShellOptions::ShowAll() const 
{

  cout << "\nShellOptions :: \n";
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first.at(0) != '_' ) 
      cout << (*ci).first << " = " << (*ci).second << endl;
  }

} /* ShowAll */

/* ======================================================================== */
/*!
 \brief Show all the databse entries

*/
/* ------------------------------------------------------------------------ */

void ShellOptions::ShowReallyAll() const 
{

  cout << "\nShellOptions :: \n";
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    cout << (*ci).first << " = " << (*ci).second << endl;
  }

} /* ShowReallyAll */

/* ======================================================================== */
/*!
 \brief Add an entry to the databse

 This method add an entry to the databse. First, it checks that this
 entry does not exist. If it exists, the method returns \c
 false. Otherwise, it adds the entry and returns \c true.
 
*/
/* ------------------------------------------------------------------------ */

bool ShellOptions::AddOption( const string input, const string value )
{

  // check that "input" has not been already inserted
  if( this->HaveOption(input) == true )
    return false;
  
  OptionDatabase[input] = value;
  return true;

} /* AddOption */

/* ======================================================================== */
/*!
 \brief Modify the value of a database entry.

 This method modifies the value of a database entry. If the entry does
 not exist in the database, return \c false. Otherwise, returns \c true.
 
*/
/* ------------------------------------------------------------------------ */

bool ShellOptions::SetOption( const string input, const string value )
{

  // check that "input" has not been already inserted
  if( this->HaveOption(input) == true )
    return false;
  
  OptionDatabase[input] = value;
  return true;

} /* SetOption */

/* ======================================================================== */
/*!
 \brief Returns the name of the program as a C++ string.

*/
/* ------------------------------------------------------------------------ */

string  ShellOptions::GetProgramName( void )
{
  return OptionDatabase["_PROGRAM_NAME_"];
}

/* ======================================================================== */
/*!
  
 \brief Returns the value of the environmenta variable \c str as an integer.

 This methods returns the value of the environmenta variable \c
 str. If the variable does not exists, returns \c 0.

*/
/* ------------------------------------------------------------------------ */

int ShellOptions::GetIntShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0;
  
} /* GetIntShellVariable */

/* ======================================================================== */
/*!
  
 \brief Returns the value of the environmenta variable \c str as a double.
 
 This methods returns the value of the environmenta variable \c
 str. If the variable does not exists, returns \c 0.0.

*/
/* ------------------------------------------------------------------------ */

double ShellOptions::GetDoubleShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0.0;
  
} /* GetDoubleShellVariable */

/* ======================================================================== */
/*!
  
 \brief Returns the value of the environmenta variable \c str as a C++ string.

 This methods returns the value of the environmenta variable \c
 str. If the variable does not exists, returns \c "".

*/
/* ------------------------------------------------------------------------ */

string ShellOptions::GetCharShellVariable( const char *str ) 
{

  char * buffer;

  buffer = getenv( str );
  if( buffer == NULL )
    return( "" );

  return( buffer );
  
} /* GetCharShellVariable */






#define ONE

// MAIN driver
#ifdef ONE
 int main(int argc, char *argv[])
  {

   ShellOptions Args(argc,argv);
   
   int nx = Args.GetIntOption("-nx", 123);
   int ny = Args.GetIntOption("-ny", 145);
   double tol = Args.GetDoubleOption("-tol", 1e-12);
   string solver = Args.GetStringOption("-solver");

   cout << "nx = " << nx << endl;
   cout << "ny = " << ny << " (default value)" << endl;
   cout << "tol = " << tol << endl;
   cout << "solver = " << solver << endl;

   return 0;
   
 }
#else
int main(int argc, char *argv[])
{

   ShellOptions Args(argc,argv);
   int nx = Args.GetIntOption("-nx", 123);
   int ny = Args.GetIntOption("-ny", 145);
   double tol = Args.GetDoubleOption("-tol", 1e-12);

   cout << "nx = " << nx << " and ny = " << ny << endl;
   
   cout << "have arg -nx : " << Args.HaveOption("-nx") << std::endl;
   cout << "have arg -nz : " << Args.HaveOption("-nz") << std::endl;

   Args.SetOption("-residyal","12.3");
   Args.SetOption("-nz","12");

   cout << "after setting, have arg -nz : " << Args.HaveOption("-nz") << std::endl;

   cout << "Program name : " << Args.GetProgramName() << endl;

   Args.ShowAll();

   cout << "now getting the shell variable HOSTNAME\n";
   string hostname( Args.GetCharShellVariable("HOSTNAME") );
   cout << hostname << endl;
   
   cout << "now getting the shell variable HOST_NAME\n";
   string hostname2( Args.GetCharShellVariable("HOST_NAME") );
   cout << hostname2 << endl;

   Args.ShowAll();
   
   return 0;
   
}

#endif
