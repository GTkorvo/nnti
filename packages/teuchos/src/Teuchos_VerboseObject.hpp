// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_VERBOSE_OBJECT_HPP
#define TEUCHOS_VERBOSE_OBJECT_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"

namespace Teuchos {

/** \brief Non-templated base class for objects that can print their
 * activities to a stream.
 *
 * Objects that derive from this interface print to a default class-owned
 * (i.e. static) output stream object (set using <tt>setDefaultOStream()</tt>)
 * or the output stream can be set on an object-by-object basis using
 * <tt>setOStream()</tt>.
 *
 * The output stream type is <tt>FancyOStream</tt> which allows for automated
 * indentation (using the <tt>OSTab</tt> class) and has other useful features.
 */
class VerboseObjectBase {
public:

  /** \name Public static member functions */
  //@{

  /** \brief Set the default output stream object.
   *
   * If this function is not called, then a default stream based on
   * <tt>std::cout</tt> is used.
   */
  static void setDefaultOStream( const RefCountPtr<FancyOStream> &defaultOStream );

  /** \brief Get the default output stream object. */
  static RefCountPtr<FancyOStream> getDefaultOStream();

  //@}

  /** \name Constructors/Initializers */
  //@{
  
  /** \brief . */
  virtual ~VerboseObjectBase() {}
  
  /** \brief Calls <tt>initializeVerboseObject()</tt>.
   */
  explicit
  VerboseObjectBase(
    const RefCountPtr<FancyOStream>   &oStream  = Teuchos::null
    );
  
  /** \brief Calls <tt>initializeVerboseObject()</tt>.
   */
  virtual void initializeVerboseObjectBase(
    const RefCountPtr<FancyOStream>   &oStream  = Teuchos::null
    );

  /** \brief Override the output stream for <tt>*this</tt> object */
  virtual VerboseObjectBase& setOStream(const RefCountPtr<FancyOStream> &oStream);

  /** \brief Set line prefix name for this object */
  virtual VerboseObjectBase& setLinePrefix(const std::string &linePrefix);

  //@}

  /** \name Query functions */
  //@{

  /** \brief Return the output stream to be used.
   *
   * If <tt>setOStream(
   */
  virtual RefCountPtr<FancyOStream> getOStream() const;

  /** \brief Get the line prefix for this object */
  virtual std::string getLinePrefix() const;

  //@}

  /** \name Utilities */
  //@{

  /** \brief Create a tab object which sets the number of tabs and optionally the line prefix.
   *
   * \param tabs  [in] The number of relative tabs to add (if <tt>tabs > 0</tt>) or remove (if <tt>tabs < 0</tt>).
   *              If <tt>tabs == OSTab::DISABLE_TABBING</tt> then tabbing will be turned off temporarily.
   *
   * \param linePrefix
   *              [in] Sets a line prefix that overrides <tt>this->getLinePrefix()</tt>.
   *
   * The side effects of these changes go away as soon as the returned
   * <tt>OSTab</tt> object is destroyed at the end of the block of code.
   *
   * Returns <tt>OSTab( this->getOStream(), tabs, linePrefix.length() ? linePrefix : this->getLinePrefix() )</tt>
   */
  virtual OSTab getOSTab(const int tabs = 1, const std::string &linePrefix = "") const;

  //@}
  
private:

  static RefCountPtr<FancyOStream>   defaultOStream_;

  RefCountPtr<FancyOStream>   thisOStream_;
  std::string                 thisLinePrefix_;

};

/** \brief Templated base class for objects that can print their activities to
 * a stream and have a verbosity level.
 *
 * Objects that derive from this interface print to a default class-owned
 * (i.e. static) output stream object (set using <tt>setDefaultOStream()</tt>)
 * or the output stream can be set on an object-by-object basis using
 * <tt>setOStream()</tt> .  In addition, each object, by default, has a
 * verbosity level that is shared by all objects (set using
 * <tt>setDefaultVerbosityLevel()</tt>) or can be set on an object-by-object
 * basis using <tt>setVerbLevel()</tt>.
 *
 * The output stream type is <tt>FancyOStream</tt> which allows for automated
 * indentation (using the <tt>OSTab</tt> class) and has other useful features.
 */
template<class ObjectType>
class VerboseObject : virtual public VerboseObjectBase {
public:

  /** \name Public static member functions */
  //@{

  /** \brief Set the default verbosity level.
   *
   * If not called, then the default verbosity level is <tt>VERB_DEFAULT</tt>
   */
  static void setDefaultVerbLevel( const EVerbosityLevel defaultVerbLevel);

  /** \brief Get the default verbosity level. */
  static EVerbosityLevel getDefaultVerbLevel();

  //@}

  /** \name Constructors/Initializers */
  //@{
  
  /** \brief Calls <tt>initializeVerboseObject()</tt>.
   */
  explicit
  VerboseObject(
    const EVerbosityLevel              verbLevel = VERB_DEFAULT  // Note, this must be the same as the default value for defaultVerbLevel_
    ,const RefCountPtr<FancyOStream>   &oStream  = Teuchos::null
    );
  
  /** \brief Calls <tt>initializeVerboseObject()</tt>.
   */
  virtual void initializeVerboseObject(
    const EVerbosityLevel              verbLevel = VERB_DEFAULT  // Note, this must be the same as the default value for defaultVerbLevel_
    ,const RefCountPtr<FancyOStream>   &oStream  = Teuchos::null
    );

  /** \brief Override the verbosity level for <tt>*this</tt> object */
  virtual VerboseObject& setVerbLevel(const EVerbosityLevel verbLevel);

  //@}

  /** \name Query functions */
  //@{

  /** \brief Get the verbosity level */
  virtual EVerbosityLevel getVerbLevel() const;

  //@}
  
private:

  static EVerbosityLevel  defaultVerbLevel_;
  EVerbosityLevel         thisVerbLevel_;

};

// //////////////////////////////////
// Template defintions

// Private static data members

template<class ObjectType>
EVerbosityLevel
VerboseObject<ObjectType>::defaultVerbLevel_ = VERB_DEFAULT;

// Public static member functions

template<class ObjectType>
void VerboseObject<ObjectType>::setDefaultVerbLevel( const EVerbosityLevel defaultVerbLevel)
{
  defaultVerbLevel_ = defaultVerbLevel;
}

template<class ObjectType>
EVerbosityLevel VerboseObject<ObjectType>::getDefaultVerbLevel()
{
  return defaultVerbLevel_;
}

// Constructors/Initializers

template<class ObjectType>
VerboseObject<ObjectType>::VerboseObject(
  const EVerbosityLevel              verbLevel
  ,const RefCountPtr<FancyOStream>   &oStream
  )
{
  this->initializeVerboseObject(verbLevel,oStream);
}

template<class ObjectType>
void VerboseObject<ObjectType>::initializeVerboseObject(
  const EVerbosityLevel              verbLevel
  ,const RefCountPtr<FancyOStream>   &oStream
  )
{
  thisVerbLevel_ = verbLevel;
  this->initializeVerboseObjectBase(oStream);
}

template<class ObjectType>
VerboseObject<ObjectType>&
VerboseObject<ObjectType>::setVerbLevel(const EVerbosityLevel verbLevel)
{
  thisVerbLevel_ = verbLevel;
  return *this;
}

// Query functions

template<class ObjectType>
EVerbosityLevel VerboseObject<ObjectType>::getVerbLevel() const
{
  if(thisVerbLevel_ == VERB_DEFAULT)
    return defaultVerbLevel_;
  return thisVerbLevel_;
}

// ////////////////////////////////////////
// Initialization classes to make sure
// that static data is initialized
// before it is used in other pre-main()
// code.  See Item #47 in Effective C++
// ////////////////////////////////////////

/*

class InitializeVerboseObjectBase {
public:
  InitializeVerboseObjectBase();
private:
  static unsigned short int count_;
};

static InitializeVerboseObjectBase initializeVerboseObjectBase;

*/

} // namespace Teuchos

#endif // TEUCHOS_VERBOSE_OBJECT_HPP
