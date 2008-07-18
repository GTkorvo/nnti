/* @HEADER@ */
/* @HEADER@ */

#ifndef TSF_OBJECTWITHVERBOSITY_H
#define TSF_OBJECTWITHVERBOSITY_H

#include "TSFConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace TSFExtended
{
using Teuchos::RefCountPtr;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;

/** \enum Verbosity settings */
enum VerbositySetting {VerbSilent=0, VerbLow=1, VerbMedium=2, 
                       VerbHigh=3, VerbExtreme=4};

/** 
 * Defines traits for getting default verbosity parameters from a class
 */
template <class X>
class VerbosityTraits
{
public:
  static RefCountPtr<ParameterList> defaultVerbParams()
    {return X::defaultVerbParams();}
};

/**
 * ObjectWithVerbosity and the related verbosity() method
 * provide an interface for getting/setting
 * verbosity flags for classes or instances. 
 *
 * All objects start out with a verbosity setting of VerbSilent.
 * 
 * You can set verbosity for a single instance of a class, or for
 * the whole class. To set for an instance, use the verbosity()
 * member function, for example,
 * \code
 * Mesh mesh1 = reader1.getMesh();
 * Mesh mesh2 = reader2.getMesh();
 * Mesh mesh3 = reader3.getMesh();
 * mesh1.verbosity() = VerbHigh;
 * \endcode
 * which sets the verbosity of <tt>mesh1</tt> to VerbHigh and leaves
 * those of <tt>mesh2</tt> and <tt>mesh3</tt> unchanged.
 *
 * Alternatively, you can set a default verbosity for an entire
 * class, for example,
 * \code
 * Mesh mesh1 = reader1.getMesh();
 * Mesh mesh2 = reader2.getMesh();
 * Mesh mesh3 = reader3.getMesh();
 * mesh1.verbosity() = VerbHigh;
 * verbosity<Mesh>() = VerbMedium;
 * \endcode
 * which sets the default verbosity to VerbMedium. Since <tt>mesh1</tt>
 * has its own verbosity setting of VerbHigh, 
 * it will use it rather than the
 * default, but <tt>mesh2</tt> and <tt>mesh3</tt> will use VerbMedium.
 * 
 */
template <class X>
class ObjectWithVerbosity
{
public:
  /** \deprecated Construct, starting silent */
  ObjectWithVerbosity() : verbosity_(classVerbosity()), setLocally_(false) {;}

  /** \deprecated Read-only access to the verbosity */
  VerbositySetting verbosity() const 
    {
      if (setLocally_) return verbosity_;
      return classVerbosity();
    }

  /** \deprecated Writeable access to the verbosity setting */
  VerbositySetting& verbosity() 
    {
      setLocally_ = true; 
      return verbosity_;
    }

  /** \deprecated Writeable access to the default verbosity for the class */
  static VerbositySetting& classVerbosity() 
    {
      static VerbositySetting rtn = VerbSilent;
      return rtn;
    }


private:
  /** */
  VerbositySetting verbosity_;

  /** */
  bool setLocally_;
};

/** 
 * \relates ObjectWithVerbosity
 * Global method for setting verbosity of a class
 */
template <class X> VerbositySetting& verbosity() 
{
  return X::classVerbosity();
}

template <class X> 
class ParameterControlledObjectWithVerbosity 
  : public ObjectWithVerbosity<X>
{
public:
  /** \deprecated Construct, starting silent */
  ParameterControlledObjectWithVerbosity() : ObjectWithVerbosity<X>() {;}

  /** Construct with a parameter list controlling the verbosity settings */
  ParameterControlledObjectWithVerbosity(const std::string& objName, const ParameterList& p)
    : ObjectWithVerbosity<X>(),
    verbControlParams_() 
    {
      RefCountPtr<ParameterList> defaults = VerbosityTraits<X>::defaultVerbParams();
      TEST_FOR_EXCEPTION(defaults->name() != objName, std::runtime_error,
        "mismatched ParameterList names for verbosity control: expected "
        << defaults->name() << ", got " << objName);
      TEST_FOR_EXCEPTION(defaults->name() != p.name(), std::runtime_error,
        "mismatched ParameterList names for verbosity control: expected "
        << defaults->name() << ", got " << p.name());
      verbControlParams_ = mergeParams(*defaults, p);
    }

  /** */
  int verbLevel(const std::string& context) const 
    {
      const ParameterEntry* pe = verbControlParams_.getEntryPtr(context);
      TEST_FOR_EXCEPTION(pe==0, std::runtime_error,
        "parameter with name \"" << context << "\" not found in verbosity "
        "control parameter list " << verbControlParams_);
      TEST_FOR_EXCEPTION(pe->getAny().type() != typeid(int),
        std::runtime_error,
        "context parameter name \"" 
        << context << "\" does not have an integer value in verbosity "
        "control parameter list " << verbControlParams_);

      return Teuchos::any_cast<int>(pe->getAny());
    }

  /** */
  const ParameterList& verbSublist(const std::string& name) const 
    {
      TEST_FOR_EXCEPTION(!verbControlParams_.isSublist(name),
        std::runtime_error,
        "context parameter name \"" 
        << name << "\" is not a sublist in verbosity "
        "control parameter list " << verbControlParams_);

      return verbControlParams_.sublist(name);
    }

  /** */
  ParameterList mergeParams(const ParameterList& pDef, const ParameterList& pIn) const
    {
      ParameterList rtn = pDef;
      using namespace Teuchos;

      /* replace any defaults with overriden values */
      ParameterList::ConstIterator i;
      for (i=pDef.begin(); i!=pDef.end(); i++)
      {
        const ParameterEntry& eDef = pDef.entry(i);
        const std::string& name = pDef.name(i);
        const ParameterEntry* eIn = pIn.getEntryPtr(name);
        if (eIn != NULL)
        {
          if (eIn->isList() && eDef.isList())
          {
            ParameterList sub = mergeParams(
              getValue<ParameterList>(eDef),
              getValue<ParameterList>(*eIn));
          }
          else if (eIn->isType<int>() && eDef.isType<int>())
          {
            rtn.set(name, Teuchos::any_cast<int>(eIn->getAny()));
          }
          else
          {
            TEST_FOR_EXCEPTION(eIn->isList() && !eDef.isList(), 
              std::runtime_error, "mismatched parameters in mergeParams()");
            TEST_FOR_EXCEPTION(!eIn->isList() && eDef.isList(), 
              std::runtime_error, "mismatched parameters in mergeParams()");
            TEST_FOR_EXCEPT(1);
          }
        }
      }
      return rtn;
    }
       

private:
  ParameterList verbControlParams_;
};



}





#endif
