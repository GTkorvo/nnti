// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_OBJECT_HPP_
#define _TEUCHOS_OBJECT_HPP_

/*! \file Teuchos_Object.hpp
    \brief The base Teuchos object.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DataAccess.hpp"

/*! \class Teuchos::Object
    \brief The base Teuchos class.

    The Object class provides capabilities common to all Teuchos objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
*/

namespace Teuchos
{

class Object
{
  public:
  //@{ \name Constructors/Destructor.
  //! Default Constructor.
  /*! Object is the primary base class in Teuchos.  All Teuchos class
      are derived from it, directly or indirectly.  This class is seldom
      used explictly.
  */
  Object(int tracebackModeIn = -1);

  //! Labeling Constructor.
  /*! Creates an Object with the given label.
  */
  Object(const char* label, int tracebackModeIn = -1);

  //! Copy Constructor.
  /*! Makes an exact copy of an existing Object instance.
  */
  Object(const Object& obj);

  //! Destructor.
  /*! Completely deletes an Object object.  
  */
  virtual ~Object();

  //@}
  
  //@{ \name Set methods.

  //! Define object label using a character string.
  /*! Defines the label used to describe \c this object.
  */
  virtual void setLabel(const char* label);

  //! Set the value of the Object error traceback report mode.
  /*! Sets the integer error traceback behavior.
      TracebackMode controls whether or not traceback information is printed when run time
      integer errors are detected:

      <= 0 - No information report

       = 1 - Fatal (negative) values are reported

      >= 2 - All values (except zero) reported.

      \note Default is set to -1 when object is constructed.
  */
  static void setTracebackMode(int tracebackModeValue);

  //@}

  //@{ \name Accessor methods.

  //! Access the object label.
  /*! Returns the string used to define \e this object.
  */
  virtual char* label() const;  

  //! Get the value of the Object error traceback report mode.
  static int getTracebackMode();

  //@}

  //@{ \name I/O method.

  //! Print method for placing the object in an output stream
  virtual void print(ostream& os) const;
  //@}

  //@{ \name Error reporting method.

  //!  Method for reporting errors with Teuchos objects.
  virtual int reportError(const string message, int errorCode) const 
  {
  // NOTE:  We are extracting a C-style string from Message because 
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"
	if ( (tracebackMode==1) && (errorCode < 0) )
	{  // Report fatal error
	   cerr << endl << "Error in Teuchos Object with label: " << label_ << endl 
		 << "Teuchos Error:  " << message.c_str() << "  Error Code:  " << errorCode << endl;
	   return(errorCode);
        }
	if ( (tracebackMode==2) && (errorCode != 0 ) ) 
	{
	   cerr << endl << "Error in Teuchos Object with label: " << label_ << endl 
		 << "Teuchos Error:  " << message.c_str() << "  Error Code:  " << errorCode << endl;
	   return(errorCode);
	}
	return(errorCode);
  }

  //@}

  static int tracebackMode;  

 protected:

 private:

  char* label_;

}; // class Object

/*! \relates Object
    Output stream operator for handling the printing of Object.
*/
inline ostream& operator<<(ostream& os, const Teuchos::Object& Obj)
{
  os << Obj.label() << endl;
  Obj.print(os);
 
  return os;
}

} // namespace Teuchos

// #include "Teuchos_Object.cpp"


#endif /* _TEUCHOS_OBJECT_HPP_ */
