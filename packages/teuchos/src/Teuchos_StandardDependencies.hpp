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



#ifndef TEUCHOS_STANDARDDEPENDCIES_HPP_
#define TEUCHOS_STANDARDDEPENDCIES_HPP_

/*! \file Teuchos_StandardDependencies.hpp
    \brief A collection of standard dependencies.
*/

#include "Teuchos_Dependency.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_StandardFunctionObjects.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace Teuchos{


/**
 * \brief An abstract parent class for all visual dependencies.
 *
 * IMPORTANT NOTE:
 * If a parameter becomes hidden, it's validity will not necessarily 
 * be checked. This means that it
 * is indeed possible for a non-valid ParameterList to occur. Make 
 * sure that you program code takes
 * this into account.
 */
class VisualDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   */
  VisualDependency(
    RCP<const ParameterEntry> dependee, 
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependent parameters.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   */
  VisualDependency(
    RCP<const ParameterEntry> dependee,
    ParameterEntryList dependents,
    bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependees The dependees.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   */
  VisualDependency(
    ConstParameterEntryList dependees, 
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependees The dependees.
   * @param dependents The dependets.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   */
  VisualDependency(
    ConstParameterEntryList dependees,
    ParameterEntryList dependents,
    bool showIf=true);

  /**
   * \brief Desctructor
   *
   * Simply declaring the descrutor as virtual.
   */
  virtual ~VisualDependency(){}
  
  //@}

  //! @name Attribute/Query Methods 
  //@{

  /**
   * \brief Get the state of the dependee in order to evaluate the
   * dependency.
   *
   * @return The state of the dependee.
   */
  virtual bool getDependeeState() const = 0;
  
  /**
   * \brief Determines whether or not the dependent is currently visible.
   */
  bool isDependentVisible() const;

  inline
  bool getShowIf() const{
    return showIf_;
  }

  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief Whether or not the dependent is currently visible.
   */
  bool dependentVisible_;

  /**
   * \brief Whether or not to show the dependent if the dependee is set to the value.
   */
  bool showIf_;
  
  //@}

};

/**
 * \brief An abstract base class for all validator dependencies.
 */
class ValidatorDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a ValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   */
  ValidatorDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent);

  /**
   * \brief Constructs a ValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   */
  ValidatorDependency(
    RCP<const ParameterEntry> dependee, 
    ParameterEntryList dependents);

  /**
   * \brief Desctructor
   *
   * Simply declaring the descrutor as virtual.
   */
  virtual ~ValidatorDependency(){}
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  virtual void evaluate() = 0;
  
  //@}

};

/**
 * \brief A string visual depdencies says the following about the 
 * relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a particular value, 
 * the dependent may or may not be displayed to the user in a UI.
 * 
 * The dependee of a StringVisualDependency must be of type string and 
 * can't be an array. The dependent may be any type of
 * parameter or parameter list.
 */
class StringVisualDependency : public VisualDependency{

public:

  /** \name Public types */
  //@{

  /**
   * Convience typedef representing an array of strings.
   */
  typedef Array<std::string> ValueList; 
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * Constructs a StringVisualDependency.
   *
   * @param dependee The dependee paramter.
   * @parame dependent The dependent parameter.
   * @param value The value of the depndee that affects the visiblity 
   * of the dependent.
   * @param showIf When true, the depndent will be be shown 
   * if the dependee 
   * is set to the same value as specified by the value parameter.
   * If false, the dependent will be shown only when the dependee is 
   * set to a value other than the one specified by the value parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    std::string value,
    bool showIf=true);

  /**
   * Constructs a StringVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param values The values of the depndee that affect the 
   * visiblity of the dependent.
   * @param showIf When true, the depndent will be be shown if 
   * the dependee is set to one of the values specified by the 
   * values parameter.
   * If false, the dependent will be shown only when the dependee is set 
   * to a value other than the ones specified by the values parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    const ValueList& values,
    bool showIf=true);

  /**
   * Constructs a StringVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents
   * @param value The value of the depndee that affects the visiblity 
   * of the dependent.
   * @param showIf When true, the depndent will be be shown if 
   * the dependee is set to one of the values specified by the values 
   * parameter. If false, the dependent will be shown only when the 
   * dependee is set to a value other than the ones specified by the 
   * values parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents, 
    const std::string& value,
    bool showIf=true);

  /**
   * Constructs a StringVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents
   * @param values The values of the depndee that affect the 
   * visiblity of the dependent.
   * @param showIf When true, the depndent will be be 
   * shown if the dependee is 
   * set to one of the values specified by the values parameter.
   * If false, the dependent will be shown only when the dependee 
   * is set to a value other than the ones specified by 
   * the values parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents, 
    const ValueList& values,
    bool showIf=true);
  
  //@}

  /** \name Attribute/Query Functions */
  //@{

  inline
  const ValueList& getValues() const{
    return values_;
  }

  //@}

  /** \name Overridden from VisualDependency */
  //@{

  inline bool getDependeeState() const{
    return find(values_.begin(), values_.end(), 
      getFirstDependeeValue<std::string>()) != values_.end();
  }
  
  //@}
  
  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "stringVisualDependency";
  }
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * The value used to deteremine the visiblity of the dependent.
   */
  const ValueList values_;
  
  //@}
  
};


/** \brief Speicialized class for retrieving a dummy object of type
 * StringVisualDependency.
 *
 * \relates StringVisualDependency
 */
template<>
class DummyObjectGetter<StringVisualDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * StringVisualDependency.
  */
  static RCP<StringVisualDependency> getDummyObject(){
    static RCP<StringVisualDependency > dummyObject;
     if(dummyObject.is_null()){
      dummyObject = rcp(new StringVisualDependency(
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      "i'm a dummy"));
    }
    return dummyObject;
  }
  
  //@}
  
};


/**
 * \brief A bool visual dependency says the following about the 
 * relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee is true or false, the 
 * dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a BoolVisualDependency must be of type bool and can't 
 * be an array. The dependent may be any type of parameter
 * or parameter list.
 */
class BoolVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * Constructs a BoolVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   * If false, the dependent will be shown only when the dependee is false.
   */
  BoolVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * Constructs a BoolVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependent parameters.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   * If false, the dependent will be shown only when the dependee is false.
   */
  BoolVisualDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents, 
    bool showIf=true);
  
  //@}

  /** \name Overridden from VisualDependency */
  //@{

  /** \brief . */
  inline bool getDependeeState() const{
    return getFirstDependeeValue<bool>();
  }
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "boolVisualDependency";
  }
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{
  
  /** \brief . */
  void validateDep() const;
  
  //@}

};


/** \brief Speicialized class for retrieving a dummy object of type
 * BoolVisualDependency.
 *
 * \relates BoolVisualDependency
 */
template<>
class DummyObjectGetter<BoolVisualDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * BoolVisualDependency.
  */
  static RCP<BoolVisualDependency> getDummyObject(){
    static RCP<BoolVisualDependency > dummyObject;
    if(dummyObject.is_null()){
      dummyObject = rcp(new BoolVisualDependency(
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      DummyObjectGetter<ParameterEntry>::getDummyObject()));
    }
    return dummyObject;
  }
  
  //@}
  
};

/**
 * \brief A condition visual dependency says the following about the 
 * relationship between elements in a Parameter List:
 * Depending on whether or not the dependee(s) statisfy 
 * a particual condition, the dependent may or may not be displayed to 
 * the user in a UI.
 *
 * Condition Visual Dependencies are unique in that via 
 * the Condition class, they allow for multiple dependees.
 * The dependee(s) of a ConditionVisualDependency must be expressed as a 
 * Condition and are subject to the consquential constraints. The 
 * dependent may be any type of parameter or parameter list.
 */
class ConditionVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * Constructs a ConditionVisualDependency.
   *
   *
   * @param condition The condition that must be satisfied in 
   * order to display the dependent parameter.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if 
   * the condition is true. If false, the dependent will be shown 
   * only when the condition is false.
   */
  ConditionVisualDependency(
    RCP<const Condition> condition,
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * Constructs a ConditionVisualDependency.
   *
   * @param condition The condition that must be satisfied in 
   * order to display the dependent parameter.
   * @param dependents The dependent parameters.
   * @param showIf When true, the depndent will be be shown if 
   * the condition is true. If false, the dependent will be shown 
   * only when the condition is false.
   */
  ConditionVisualDependency(
    RCP<const Condition> condition, 
    Dependency::ParameterEntryList dependents,
    bool showIf=true);
  
  //@}

  /** \name Overridden from VisualDependency */
  //@{

  inline bool getDependeeState() const{
    return condition_->isConditionTrue();
  }
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "conditionVisualDependency";
  }
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const {}
  
  //@}


private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The Condition to determine whether or not the dependent is displayed.
   */
  RCP<const Condition> condition_;
  
  //@}

};


/** \brief Speicialized class for retrieving a dummy object of type
 * ConditionVisualDependency.
 *
 * \relates ConditionVisualDependency
 */
template<>
class DummyObjectGetter<ConditionVisualDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * ConditionVisualDependency.
  */
  static RCP<ConditionVisualDependency> getDummyObject(){
    static RCP<ConditionVisualDependency> dummyObject;
    if(dummyObject.is_null()){
      dummyObject = rcp(new ConditionVisualDependency(
      DummyObjectGetter<NotCondition>::getDummyObject(),
      DummyObjectGetter<ParameterEntry>::getDummyObject()));
    }
    return dummyObject;
  }
  
  //@}
  
};


/**
 * \brief A number visual dependency says the following about 
 * the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a certain value, 
 * the dependent may or may not be displayed to the user in a UI.
 *
 * The dependee of a NumberVisualDependency must 
 * be a number type and can't be an array. The dependent may be 
 * any type of parameter or parameter list.
 */
template <class T>
class NumberVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a NumberVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param func A function that takes the dependees value, does some 
   * calculations on it, and then returns a value. If this value is 
   * greater than 0, the dependent is show. If the value returned is
   * less than or equal to zero, the dependent is not shown. If no 
   * fuction is specified, the direct value of the dependee will be used 
   * to determine the dependents visibility in a similar fashion (postive
   * numbers causing the dependent to be displayed and 0 or 
   * negative numbers causing the dependent to be hidden).
   */
  NumberVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    RCP<SingleArguementFunctionObject<T,T> > func=null)
    :VisualDependency(dependee, dependent),
    func_(func)
  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param func A function that takes the dependees value, does some 
   * calculations on it, and then returns a value. If this value is 
   * greater than 0, the dependent is show. If the value returned is
   * less than or equal to zero, the dependent is not shown. If no 
   * fuction is specified, the direct value of the dependee will be used 
   * to determine the dependents visibility in a similar fashion (postive
   * numbers causing the dependent to be displayed and 0 or 
   * negative numbers causing the dependent to be hidden).
   */
  NumberVisualDependency(
    RCP<const ParameterEntry> dependee,
    ParameterEntryList dependents,
    RCP<SingleArguementFunctionObject<T,T> > func=null)
    :VisualDependency(dependee, dependents),
    func_(func)
  {
    validateDep();
  }

  //@}

  /** \name Overridden from VisualDependency */
  //@{

  inline bool getDependeeState() const{
    if(!func_.is_null()){
      func_->setParameterValue(getFirstDependeeValue<T>());
      return func_->runFunction() > ScalarTraits<T>::zero() ? true : false;
    }
    return getFirstDependeeValue<T>();
  }
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name() + " NumberVisualDependency";
  }
  
  //@}
  
  /** \name Getter Functions */
  //@{
  
  /** \brief Gets the function associated with this dependency. */
  RCP<SingleArguementFunctionObject<T,T> > getFunctionObject(){
    return func_;
  }

  /** \brief Const version of function getter. */
  RCP<const SingleArguementFunctionObject<T,T> > getFunctionObject() const{
    return func_.getConst();
  }


  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const;
  
  //@}
  
private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief the function used to determine the
   * visibility of the dependent.
   */
    RCP<SingleArguementFunctionObject<T,T> > func_;
  
  //@}
  //
};

template<class T>
void NumberVisualDependency<T>::validateDep() const{
  RCP<const ParameterEntry> dependee = getFirstDependee();
  TEST_FOR_EXCEPTION(
    !dependee->isType<int>()
    && !dependee->isType<short>()
    && !dependee->isType<double>()
    && !dependee->isType<float>(),
    InvalidDependencyException,
    "The dependee of a "
    "Number Visual Dependency must be of a supported number type!\n"
    "Type Encountered: " << dependee->getAny().typeName() << "\n");
}
  

/** \brief Speicialized class for retrieving a dummy object of type
 * NumberVisualDependency.
 *
 * \relates NumberVisualDependency
 */
template<class T>
class DummyObjectGetter<NumberVisualDependency<T> >{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * NumberVisualDependency.
  */
  static RCP<NumberVisualDependency<T> >
    getDummyObject()
  {
    static RCP<NumberVisualDependency<T> > dummyObject;
    if(dummyObject.is_null()){
      dummyObject = rcp(new NumberVisualDependency<T>(
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      DummyObjectGetter<ParameterEntry>::getDummyObject()));
    }
    return dummyObject;
  }
  
  //@}
  
};

/**
 * \brief A NumberArrayLengthDependency says the following about the 
 * relationship between two parameters:
 * The length of the dependent's array depends on the value 
 * of the dependee.
 *
 * A NumberArrayLengthDependency must have the following characteristics:
 *
 *   \li The dependee must be either of type int or short.
 *
 *   \li The dependent must be an array.
 *
 */
template<class DependeeType, class DependentType>
class NumberArrayLengthDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a NumberArrayLengthDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param func A function specifying how the arrays length 
   * should be calculated from the dependees value.
   */
  NumberArrayLengthDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    RCP<SingleArguementFunctionObject<Teuchos_Ordinal, DependeeType> > func=null);

  /**
   * \brief Constructs a NumberArrayLengthDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param func A function specifying how the arrays length 
   * should be calculated from the dependees value.
   */
  NumberArrayLengthDependency(
    RCP<const ParameterEntry> dependee,
    ParameterEntryList dependent,
    RCP<SingleArguementFunctionObject<Teuchos_Ordinal, DependeeType> > func=null);

  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "numberArrayLengthDependency<" +
      TypeNameTraits<DependeeType>::name() + ", " +
      TypeNameTraits<DependentType>::name() +">";
  }
  
  //@}

  /** \name Getters */
  //@{

  /** \brief Gets the function associated with this dependency. */
  RCP<SingleArguementFunctionObject<Teuchos_Ordinal, DependeeType> >
    getFunctionObject()
  {
    return func_;
  }

  /** \brief Const version of function getter. */
  RCP<const SingleArguementFunctionObject<Teuchos_Ordinal, DependeeType> > 
    getFunctionObject() const
  {
    return func_.getConst();
  }

  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}
  //
private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The function used to calculate the new value of the
   * arrays length.
   */
    RCP<SingleArguementFunctionObject<Teuchos_Ordinal, DependeeType> > func_;
  
  /**
   * \brief Modifies the length of an array.
   *
   * @param newLength The new length the array should be.
   * @param dependentValue The index of the dependent array that is going to be changed.
   */
  void modifyArrayLength(
    Teuchos_Ordinal newLength, RCP<ParameterEntry> dependentToModify);
  
  //@}
  
};


/** \brief Speicialized class for retrieving a dummy object of type
 * NumberArrayLengthDependency.
 *
 * \relates NumberArrayLengthDependency
 */
template<class DependeeType, class DependentType>
class DummyObjectGetter<NumberArrayLengthDependency<DependeeType, DependentType> >{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * NumberArrayLengthDependency.
  */
  static RCP<NumberArrayLengthDependency<DependeeType, DependentType> >
    getDummyObject()
  {
    static RCP<NumberArrayLengthDependency<DependeeType, DependentType> > 
      dummyObject;
    if(dummyObject.is_null()){
      dummyObject = rcp(
        new NumberArrayLengthDependency<DependeeType, DependentType>(
        DummyObjectGetter<ParameterEntry>::getDummyObject(),
        DummyObjectGetter<ParameterEntry>::getDummyObject()));
    }
    return dummyObject;
  }
  
  //@}
  
};

template<class DependeeType, class DependentType>
NumberArrayLengthDependency<DependeeType, DependentType>::NumberArrayLengthDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  RCP<SingleArguementFunctionObject<Teuchos_Ordinal, DependeeType> > func):
  Dependency(dependee, dependent),
  func_(func)
{
  validateDep();
}

template<class DependeeType, class DependentType>
NumberArrayLengthDependency<DependeeType, DependentType>::NumberArrayLengthDependency(
  RCP<const ParameterEntry> dependee,
  ParameterEntryList dependents,
  RCP<SingleArguementFunctionObject<Teuchos_Ordinal, DependeeType> > func):
  Dependency(dependee, dependents),
  func_(func)
{
  validateDep();
}

template <class DependeeType, class DependentType>
void 
NumberArrayLengthDependency<DependeeType, DependentType>::modifyArrayLength(
  Teuchos_Ordinal newLength, RCP<ParameterEntry> dependentToModify)
{
  const Array<DependentType> originalArray = 
    any_cast<Array<DependentType> >(dependentToModify->getAny()); 
  Array<DependentType> newArray(newLength);
  Teuchos_Ordinal i;
  for(i=OrdinalTraits<Teuchos_Ordinal>::zero(); i<originalArray.size() && i<newLength; ++i){
    newArray[i] = originalArray[i];
  }

  dependentToModify->setValue(newArray,
    false, dependentToModify->docString(), dependentToModify->validator());
}

template<class DependeeType, class DependentType>
void 
NumberArrayLengthDependency<DependeeType, DependentType>::evaluate(){
  Teuchos_Ordinal newLength;
  if(!func_.is_null()){
    func_->setParameterValue(getFirstDependeeValue<DependeeType>());
    newLength = func_->runFunction();
  }
  else{
    newLength = getFirstDependeeValue<DependeeType>();
  }

  TEST_FOR_EXCEPTION(newLength < OrdinalTraits<Teuchos_Ordinal>::zero(),
    Exceptions::InvalidParameterValue,
    "Ruh Roh Shaggy! Looks like a dependency tried to set the length "
    "of the Array(s) to a negative number. Silly. You can't have "
    "an Array with a negative length!\n\n" <<
    "Error:\n" <<
    "An attempt was made to set the length of an Array to a negative "
    "number by a NumberArrayLengthDependency\n");
  for(
    ParameterEntryList::iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    modifyArrayLength(newLength, *it);
  }
}

template<class DependeeType, class DependentType>
void 
NumberArrayLengthDependency<DependeeType, DependentType>::validateDep() 
  const
{
  TEST_FOR_EXCEPTION(
    typeid(DependeeType) != getFirstDependee()->getAny().type(),
    InvalidDependencyException,
    "Ay no! The dependee parameter types don't match." << std::endl <<
    "Dependee Template Type: " << TypeNameTraits<DependeeType>::name() <<
    std::endl <<
    "Dependee Parameter Type: " << getFirstDependee()->getAny().typeName()
    << std::endl << std::endl);

  for(
    ConstParameterEntryList::const_iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    TEST_FOR_EXCEPTION(
      typeid(Teuchos::Array<DependentType>) != (*it)->getAny().type(),
        InvalidDependencyException,
        "Ay no! The dependent parameter types don't match." << std::endl <<
        "Dependent Template Type: " << 
        TypeNameTraits<DependentType>::name() << std::endl <<
        "Dependent Parameter Type: " << 
        (*it)->getAny().typeName() << std::endl << std::endl);
  }
}

/**
 * \brief A StringValidatorDependency says the following about 
 * the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should 
 * use a particular validator from
 * a given set of validators.
 *
 *
 * A StringValidatorDependency must have the following characterisitics:
 * 
 *   \li The dependee must be of type string
 *
 *   \li The validators in the ValueToValidatorMap must all be the same
 *   type.
 *
 * If the dependee takes on a value not in the valuesAndValidators
 * map, then the default validator is assigned to the dependent.
 */
class StringValidatorDependency : public ValidatorDependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Conveniece typedef
   */
  typedef std::map<std::string, RCP<const ParameterEntryValidator> > 
    ValueToValidatorMap;

  /**
   * \brief Conveniece typedef
   */
  typedef std::pair<std::string, RCP<const ParameterEntryValidator> > 
    ValueToValidatorPair;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a StringValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param valuesAndValidators A map associating string values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should 
   * be applied to the dependent based on the dependees value.
   * @param defaultValidator If a value is entered in the 
   * dependee that is not in the valuesAndValidators map,
   * this is the validator that will be assigned to the dependent.
   */
  StringValidatorDependency(
    RCP<const ParameterEntry> dependee, 
    RCP<ParameterEntry> dependent,
    ValueToValidatorMap valuesAndValidators, 
    RCP<ParameterEntryValidator> defaultValidator=null);

  /**
   * \brief Constructs a StringValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param valuesAndValidators A map associating string values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied to 
   * the dependent based on the dependees value.
   * @param defaultValidator If a value is entered in the dependee 
   * that is not in the valuesAndValidators map,
   * this is the validator that will be assigned to the dependent.
   */
  StringValidatorDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents,
    ValueToValidatorMap valuesAndValidators, 
    RCP<ParameterEntryValidator> defaultValidator = null);

  //@}

  /** \name Getters */
  //@{

  /** \brief retrieve a const reference to the ValueToValidator map being 
   * used by this StringValidatorDependency */
  const ValueToValidatorMap& getValuesAndValidators(){
    return valuesAndValidators_;
  }
  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "stringValidatorDependency";
  }
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief A map associating particular dependee values with validators 
   * that could be placed on the dependent.
   */
  ValueToValidatorMap valuesAndValidators_;

  /**
   * \brief The default validator to be used if a request is made 
   * for a value that does not
   * appear in the valuesAndValidators map.
   */
  RCP<ParameterEntryValidator> defaultValidator_;
  
  //@}
  
};


/** \brief Speicialized class for retrieving a dummy object of type
 * StringValidatorDependency.
 *
 * \relates StringValidatorDependency
 */
template<>
class DummyObjectGetter<StringValidatorDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * StringValidatorDependency.
  */
  static RCP<StringValidatorDependency >
    getDummyObject()
  {
    static RCP<StringValidatorDependency> dummyObject;
    if(dummyObject.is_null()){
      dummyObject = rcp(new StringValidatorDependency(
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      StringValidatorDependency::ValueToValidatorMap()));
    }
    return dummyObject;
  }
  
  //@}
  
};

/**
 * \brief A BoolValidatorDependency says the following about the 
 * relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a 
 * particular validator from a given set of validators.
 *
 * A BoolValidatorDependency must have the following characterisitics:
 *
 *   \li The dependee must be of type bool
 *
 *   \li The false and true validators must be the same type.
 *
 */
class BoolValidatorDependency : public ValidatorDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a BoolValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param trueValidator The validator to be used on the dependent 
   * if the dependee is set to true.
   * @param falseValidator The validator to be used on the 
   * dependent if the dependee is set to false.
   */
  BoolValidatorDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    RCP<const ParameterEntryValidator> trueValidator,
    RCP<const ParameterEntryValidator> falseValidator);

  /**
   * \brief Constructs a BoolValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param trueValidator The validator to be used on the dependent 
   * if the dependee is set to true.
   * @param falseValidator The validator to be used on the dependent 
   * if the dependee is set to false.
   */
  BoolValidatorDependency(
    RCP<const ParameterEntry> dependee,
    Dependency::ParameterEntryList dependents,
    RCP<const ParameterEntryValidator> trueValidator,
    RCP<const ParameterEntryValidator> falseValidator);

  //@}

  /** \name Overridden from Dependency */
  //@{

  void evaluate();
  
  //@}

  /** \name Getters */
  //@{
    
  /** \brief Gets the true validator */
  RCP<const ParameterEntryValidator> getTrueValidator(){
    return trueValidator_;
  }

  /** \brief Gets the false validator */
  RCP<const ParameterEntryValidator> getFalseValidator(){
    return falseValidator_;
  }
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "boolValidatorDependency";
  }
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The validators to be used when the dependee is either
   * true or false.
   */
  RCP<const ParameterEntryValidator> trueValidator_, falseValidator_;
  
  //@}

};

/** \brief Speicialized class for retrieving a dummy object of type
 * BoolValidatorDependency.
 *
 * \relates BoolValidatorDependency
 */
template<>
class DummyObjectGetter<BoolValidatorDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * BoolValidatorDependency.
  */
  static RCP<BoolValidatorDependency >
    getDummyObject()
  {
    static RCP<BoolValidatorDependency > dummyObject;
    if(dummyObject.is_null()){
      dummyObject = rcp(new BoolValidatorDependency(
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      DummyObjectGetter<ParameterEntry>::getDummyObject(),
      null, null));
    }
    return dummyObject;
  }
  
  //@}
  
};

/**
 * \brief A RangeValidatorDependency says the following about the
 * relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should 
 * use a particular validator from a given set of validators.
 *
 * A RangeValidatorDependency achieves this by associating ranges of 
 * values with validators.
 * If the dependees value falls within the one of the ranges, 
 * the validator associated with the range is
 * used on the dependent. If the value doesn't fall within
 * any of the ranges, the dependent's validator is set to null.
 * All ranges are inclusive.
 *
 * A RangeValidatorDependency must have the following characterisitics:
 *
 *   \li The dependee type must be the same as the template type.
 *
 *   \li All the validators in the rangesAndValidators_ map must be
 *   the same type.
 */
template<class T>
class RangeValidatorDependency : public ValidatorDependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Convenience typedef
   */
  typedef std::pair<T,T> Range;

  /**
   * \brief Convenience typedef
   */
  typedef std::map<Range, RCP<const ParameterEntryValidator> > 
    RangeToValidatorMap;

  /**
   * \brief Convenience typedef
   */
  typedef std::pair<Range, RCP<const ParameterEntryValidator> > 
    RangeValidatorPair;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a RangeValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param rangesAndValidators A map associating ranges of values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied 
   * to the dependent based on the dependees value.
   */
  RangeValidatorDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    RangeToValidatorMap rangesAndValidators)
    :ValidatorDependency(dependee, dependent),
    rangesAndValidators_(rangesAndValidators)
  {
    validateDep();
  }

  /**
   * \brief Constructs a RangeValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param rangesAndValidators A map associating ranges of values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied 
   * to the dependent based on the dependees value.
   */
  RangeValidatorDependency(
    RCP<const ParameterEntry> dependee,
    Dependency::ParameterEntryList dependents,
    RangeToValidatorMap rangesAndValidators)
    :ValidatorDependency(dependee, dependents),
    rangesAndValidators_(rangesAndValidators)
  {
    validateDep();
  }

  //@}

  /** \name Getters */
  //@{

  /** \brief . */
  const RangeToValidatorMap& getRangeToValidatorMap() const{
    return rangesAndValidators_;
  }
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name() + " RangeValidatorDependency";
  }
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{
  
  /** \brief . */
  void validateDep() const;
  
  //@}

  
private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief A map associating ranges with validators.
   */
  RangeToValidatorMap rangesAndValidators_;
  
  void setDependentsToValidator(RCP<const ParameterEntryValidator> toSet);

  //@}

};

template<class T>
void RangeValidatorDependency<T>::evaluate(){
  typename RangeToValidatorMap::const_iterator it;
  T dependeeValue = getFirstDependeeValue<T>();
  for(
    it = rangesAndValidators_.begin(); 
    it != rangesAndValidators_.end(); 
    ++it)
  {
    T min = it->first.first;
    T max = it->first.second;
    if(dependeeValue >= min && dependeeValue <=max){
       setDependentsToValidator(it->second);
      return;
    }
  }
  setDependentsToValidator(null); 
}

template<class T>
void RangeValidatorDependency<T>::validateDep() const{
  RCP<const ParameterEntry> dependee = getFirstDependee();
  TEST_FOR_EXCEPTION(dependee->getAny().type() != typeid(T),
    InvalidDependencyException,
    "The dependee of a RangeValidatorDependency must be the same type as " <<
    "The RangeValidatorDependency template type!" << std::endl <<
    "Dependee Type: " << dependee->getAny().typeName() << std::endl <<
    "Templated Type: " << TypeNameTraits<T>::name() << std::endl << std::endl);

  typename RangeToValidatorMap::const_iterator it = 
    rangesAndValidators_.begin();
  RCP<const ParameterEntryValidator> firstValidator = it->second;
  for(; it!=rangesAndValidators_.end(); ++it){
    TEST_FOR_EXCEPTION( typeid(*firstValidator) != typeid(*(it->second)),
      InvalidDependencyException,
      "Ay no! All of the validators in a RangeValidatorDependency "
      "must have the same type.");
  }
}

template<class T>
void RangeValidatorDependency<T>::setDependentsToValidator(
  RCP<const ParameterEntryValidator> toSet)
{
  typename ParameterEntryList::const_iterator it;
  for(
    it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    (*it)->setValidator(toSet);
  }
}
/** \brief Speicialized class for retrieving a dummy object of type
 * RangeValidatorDependency.
 *
 * \relates RangeValidatorDependency
 */
template<class T>
class DummyObjectGetter<RangeValidatorDependency<T> >{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * RangeValidatorDependency.
  */
  static RCP<RangeValidatorDependency<T> >
    getDummyObject()
  {
    static RCP<RangeValidatorDependency<T> > dummyObject;
    if(dummyObject.is_null()){
      typename RangeValidatorDependency<T>::RangeToValidatorMap dummyMap;
      dummyObject = rcp(new RangeValidatorDependency<T>(
        DummyObjectGetter<ParameterEntry>::getDummyObject(),
        DummyObjectGetter<ParameterEntry>::getDummyObject(),
        dummyMap));
    }
    return dummyObject;
  }
  
  //@}
  
};


} //namespace Teuchos
#endif //TEUCHOS_STANDARDDEPENDCIES_HPP_
