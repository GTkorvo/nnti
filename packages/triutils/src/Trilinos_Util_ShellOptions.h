#ifndef __TRILINOS_UTIL_SHELL_OPTIONS
#define __TRILINOS_UTIL_SHELL_OPTIONS

#include <iostream>
#include <string>
// if STL cannot be used, uncomment the line below
//#define TRILINOS_UTIL_MAP_WITH_STL
#ifdef TRILINOS_UTIL_MAP_WITH_STL
#include <map>
#endif

using namespace std;

class Trilinos_Util_Map {

public:


  
  Trilinos_Util_Map();

  ~Trilinos_Util_Map()
  {
#ifndef TRILINOS_UTIL_MAP_WITH_STL
    delete [] MapName_;
    delete [] MapValue_;
    Allocated_ = 0;
#endif
  }

  //@}

  //@{ \name Insertion/Removal methods.
  
  //! Gets the value of the specified option as an integer
  /*! This method returns the integer value assigned to option \c input.
    If \c input is not in the database, or it cannot be converted to an
    integer, returns 0.
  */
  virtual int    GetInt( const string input) ;

  //! Gets the value of the specified option as an integer. If not found, returns the specified default value.
  virtual int    GetInt( const string input, const int def_value) ;

  //! Gets the value of the specified option as a double  
  /*! This method returns the integer value assigned to option \c input.
    If \c input is not in the database, or it cannot be converted to an
    double, returns 0.
  */  
  virtual double GetDouble( const string input) ;
  
  //! Gets the value of the specified option as a double. If not found, returns the specified default value.
  virtual double GetDouble( const string input, const double def_value) ;

  //! Get the value of the specified option as a C++ string. 
  virtual string GetString( const string input) ;

 //! Gets the value of the specified option as a string. If not found, returns the specified default value.
  virtual string GetString( const string input, const string def_value ) ;

  //! Modify the value of a database entry.
  /*!  This method modifies the value of a database entry. If the entry
    does not exist in the database, return \c false. Otherwise, returns
    \c true.
  */
  virtual bool Set(const string input, const string value);
  virtual bool Set(const string input, const int value);

  //! Add an entry to the databse
  /*!  This method add an entry to the databse. First, it checks that
    this entry does not exist. If it exists, the method returns \c
    false. Otherwise, it adds the entry and returns \c true.
  */
  virtual bool   Add( const string input, const string value );
    
  inline bool SetLabel(string Label) {
    Label_ = Label;
  }
  
  inline string GettLabel(string Label) {
    return( Label_ );
  }

  //@}

  //@{ \name Query methods.

  /*! \brief Check wheter an option is in the database or not
    
  This method checks whether option \c input is in the databse or not.
  It returns \c true if it is, \c false otherwise.
  */
  virtual bool   Have(const string input);

  //@}

  //@{ \name Miscellaneous methods.

  //! Show all the databse entries
  virtual void    ShowAll() const;

  //! Show all the databse entries, including entries beginning with "_"
  virtual void    ShowReallyAll() const;

  virtual void Reset(void);
  
  //@}

  //@{ \name Friend functions .

  friend ostream & operator << (ostream & os,
				const Trilinos_Util_Map & S);

  //@}
  
private:

  string Label_;
  
#ifdef TRILINOS_UTIL_MAP_WITH_STL
  // map containing the arguments
  map<string,string> Map_;
#else
  // for non-STL, I use a very simple array of strings
  // for the option names and parameters. The dimension
  // of the arrays is hardwired
  string * MapName_;
  string * MapValue_;
  int NumEntries_;
  int Allocated_;
#endif
  
};

//! Trilinos_Util_ShellOptions: A class for managing the input arguments and variables.

/* ======================================================================== */
/*!
 \class Trilinos_Util_ShellOptions

 Using Trilinos_Util_ShellOptions, it is easy to handle input line arguments and shell
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
   int nx = Args.GetInt("-nx", 123);
   int ny = Args.GetInt("-ny", 145);
   double tol = Args.GetDouble("-tol", 1e-12);
   string solver = Args.GetInt("-solver");

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
 - GetInt
 - GetDouble
 - GetString
 
 If option name is not found in the database, a value of 0, 0.0 or an
 empty string is returned. If needed, the user can also specify a
 default value to return when the option name is not found in the
 database. Method \c HaveOption can be used to query the database for
 an option.

 The user can modify the database as well, using
 - Set
 - Add

 (GetInt, GetDouble, GetString, Set and Add are derived from the base class,
 Trilinos_Util_Map).
 
 Finally, the user can retrive the integer, double or string value of a
 shell environmental variable using:
 - GetIntShellVariable
 - GetDoubleShellVariable
 - GetStringShellVariable
 
 \date Albuquerque, 19-Jan-04
 
 \author Marzio Sala, SNL 9214

*/
// ================================================ ====== ==== ==== == =

class Trilinos_Util_ShellOptions : public Trilinos_Util_Map 
{

 public:

  //@{ \name Constructors/Destructor.
  //! Trilinos_Util_ShellOptions constructor using the options given at the shell line.
  Trilinos_Util_ShellOptions(int argc, char *argv[] );

  // ============= //
  // query and add //
  // ------------- //
  
  //@}

  //@{ \name Query methods.

  //!  Returns the name of the program as a C++ string.
  virtual string GetProgramName(void);

  // =============== //
  // shell variables //
  // --------------- //
  
  //!Returns the value of the environmenta variable \c str as an integer.
  /*! This methods returns the value of the environmenta variable \c
    str. If the variable does not exists, returns \c 0. */
  virtual int    GetIntShellVariable( const char *str );
  
  //!Returns the value of the environmenta variable \c str as an double.
  /*! This methods returns the value of the environmenta variable \c
    str. If the variable does not exists, returns \c 0.0. */
  virtual double GetDoubleShellVariable( const char *str );

  //!Returns the value of the environmenta variable \c str as a C++ string.
  /*! This methods returns the value of the environmenta variable \c
    str. If the variable does not exists, returns \c "". */
  virtual string GetStringShellVariable( const char *str );

  //@}
  
};

void   Trilinos_Util_ShellOptions_Set(int argc, char * argv[]);
int    Trilinos_Util_ShellOptions_GetInt( const string input) ;
int    Trilinos_Util_ShellOptions_GetInt( const string input, const int def_value) ;
double Trilinos_Util_ShellOptions_GetDouble( const string input) ;
double Trilinos_Util_ShellOptions_GetDouble( const string input, const double def_value) ;
string Trilinos_Util_ShellOptions_GetString( const string input) ;
string Trilinos_Util_ShellOptions_GetString( const string input, const string def_value ) ;
bool   Trilinos_Util_ShellOptions_Set(const string input, const string value);
bool   Trilinos_Util_ShellOptions_Set(const string input, const int value);
bool   Trilinos_Util_ShellOptions_Have(const string input);
bool   Trilinos_Util_ShellOptions_Add( const string input, const string value );
void   Trilinos_Util_ShellOptions_ShowAll();
void   Trilinos_Util_ShellOptions_ShowReallyAll();

// ================================================ ====== ==== ==== == =

class Trilinos_Util_FileOptions : public Trilinos_Util_Map 
{

 public:

  Trilinos_Util_FileOptions(const char FileName[] );
  ~Trilinos_Util_FileOptions();

  virtual string GetFileName(void) const;

  //  virtual bool SetCommentChar(char c);
  virtual void SetCommentChars(const string c);

  //  virtual bool SetSeparationChar(char c);
  virtual void SetSeparationChars(const string c);

  virtual int ReadFile();
  virtual int ReadFile(const char FileName[]);

private:
  bool FileHasBeenRead_;
  string FileName_;
  string CommentChars_;
  string SeparationChars_;
};
#endif
