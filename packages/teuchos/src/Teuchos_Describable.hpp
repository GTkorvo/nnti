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

#ifndef TEUCHOS_DESCRIBABLE_HPP
#define TEUCHOS_DESCRIBABLE_HPP

#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace Teuchos {

/** \brief Base class for all objects that can describe themselves and
 * their current state.
 * 
 * This base class is designed to be a minimally invasive approach for
 * allowing subclasses to optionally provide detailed debug-style information
 * about their current state.  This interface has just two virtual member
 * functions, <tt>describe(void)</tt> and <tt>description()</tt>, which both
 * have default implementations.  The shorter version of
 * <tt>describe(void)</tt> (which takes no arguments and returns an
 * <tt>std::string</tt> object) is meant for very short descriptions while the
 * longer version of <tt>description()</tt> takes and returns a
 * <tt>FancyOStream</tt> object and is designed for more detailed formated
 * output.
 *
 * Since both of these <tt>description()</tt> functions have reasonable
 * default implementations, when a subclass inherits from this base
 * class, no virtual functions need to be overridden to start with.
 * However, when debugging time comes, one or both of these functions
 * should be overridden to provide more useful information.
 *
 * ToDo: Include an example/testing function for a few different use
 * cases to demonstrate how to use this interface properly.
 *
 * \ingroup teuchos_outputting_grp
 */
class Describable {
public:

	/// Default value for <tt>verLevel</tt> in <tt>description()</tt>
	static const EVerbosityLevel   verbLevel_default;

	/** \brief . */
	virtual ~Describable() {}

	/** \name Public virtual member functions */
	//@{

	/** \brief Return a simple description (usually just one line) of this object.
	 *
	 * The default implementation just returns <tt>typeid(*this).name()</tt> but
	 * a subclass can modify this if needed.  Note that some compilers return a
	 * mangled name from <tt>std::type_info::name()</tt> (e.g. g++ version 3.4.x
	 * and before) that is hard for non-g++ developers to read.  Therefore, it
	 * is usually beneficial to override this function to build a more
	 * human-readable name for a subclass, especially if templating is used.
	 */
	virtual std::string description() const;

	/** \brief Print the object with some verbosity level to an
	 * <tt>FancyOStream</tt> object.
	 *
	 * \param  out   
	 *               [in] The <tt>FancyOStream</tt> object that output is sent to.
	 * \param  verbLevel
	 *               [in] Determines the level of verbosity for which the
	 *               the object will be printed.  If <tt>verbLevel==VERB_DEFAULT</tt>
	 *               (which is the default value), then the verbosity level will
	 *               be determined by the <tt>*this</tt> object (i.e. perhaps through the
	 *               <tt>ObjectWithVerbosity</tt> interface).  It is up to <tt>*this</tt>
	 *               how to interpret the level represented by <tt>verbLevel</tt>.
	 *               The default value is <tt>VERB_DEFAULT</tt>.
	 *
	 * In order for this function to work effectively for independently
	 * developed classes, a general consensus needs be reached as to
	 * what the various verbosity levels represented in
	 * <tt>verbLevel</tt> mean in relation to the amount of output
	 * produced.
	 *
	 * A default implementation of this function is provided that simply
	 * performs:
   \verbatim

   return out << this->description() << std::endl; \endverbatim
	 *
	 * A subclass should override this function to provide more
	 * interesting and more useful information about the object.
	 */
	virtual void describe(
		FancyOStream                &out
		,const EVerbosityLevel      verbLevel     = verbLevel_default
		) const;
  
};

// Describable stream manipulator state class
//
// This is not a class that a user needs to see and that is why it is not
// being given doxygen documentation!
struct DescribableStreamManipulatorState {
  const Describable          &describable;
	const EVerbosityLevel      verbLevel;
  DescribableStreamManipulatorState(
    const Describable          &_describable
    ,const EVerbosityLevel     _verbLevel=VERB_MEDIUM
    )
    :describable(_describable)
    ,verbLevel(_verbLevel)
    {}
};

/** \brief Describable output stream maniuplator.
 *
 * This simple function allows you to insert ouptut from
 * <tt>Describable::describe()</tt> right in the middle of a chain of
 * insertion operations.  For example, you can write:
 
 \code

  void someFunc( const Teuchos::Describable &obj )
  {
    ...
    std::cout << "The object is described as " << describe(obj,Teuchos::VERB_MEDIUM);
    ...
  }

 \endcode

 * \relates Describable
 */
inline DescribableStreamManipulatorState describe(
  const Describable          &describable
  ,const EVerbosityLevel     verbLevel     = Describable::verbLevel_default
  )
{
  return DescribableStreamManipulatorState(describable,verbLevel);
}

/** \brief Output stream operator for Describable manipulator.
 *
 * To call this function use something like:
 
 \code

  void someFunc( const Teuchos::Describable &obj )
  {
    ...
    std::cout << "The object is described as " << describe(obj,Teuchos::VERB_MEDIUM);
    ...
  }

 \endcode

 * Note: The input <tt>std::ostream</tt> is casted to a <tt>FancyOStream</tt>
 * object before calling <tt>Describable::describe()</tt> on the underlying
 * <tt>Describable</tt> object.  There is no way around this since this
 * function must be written in terms of <tt>std::ostream</tt> rather than
 * <tt>FancyOStream</tt> if one is to write compound output statements
 * involving primitive data types.
 *
 * \relates Describable
 */
inline
std::ostream& operator<<( std::ostream& os, const DescribableStreamManipulatorState& d )
{
  d.describable.describe(*getFancyOStream(Teuchos::rcp(&os,false)),d.verbLevel);
  return os;
}

//
// RAB: Note: The above function works with an std::ostream object even
// through Describable::describe(...) requires a FancyOStream object.  We must
// write the stream manipulator in terms of std::ostream, or compound right
// statements like:
//
//  void foo( FancyOStream &out, Describable &d, EVerbLevel verbLevel )
//    {
//      out << "\nThis is the describable object d:" << describe(d,verbLevel);
//    }
//
// will not work currectly.  The problem is that the first output
//
//   out << "\nThis is the describable object d:"
//
//  must return a reference to an std::ostream object.  This should mean that
//  the next statement, which is basically:
//
//    static_cast<std::ostream&>(out) << DescribableStreamManipulatorState
//
// should not even compile.  However, under gcc 3.4.3, the code did compile
// but did not call the above function.  Instead, it set up some type of
// infinite recursion that resulted in a segfault.
//

} // namespace Teuchos

#endif // TEUCHOS_DESCRIBABLE_HPP
