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


#ifndef TEUCHOS_STANDARDCONDITION_HPP_
#define TEUCHOS_STANDARDCONDITION_HPP_

/*! \file Teuchos_StandardConditions.hpp
    \brief Standard Conditions to be used.
*/

#include "Teuchos_Condition.hpp"
#include "Teuchos_InvalidConditionException.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Teuchos{

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<BinaryLogicalCondition>{
public:
  static std::string name(){ return "BinaryLogicalCondition"; }
  static std::string concreteName(const BinaryLogicalCondition&){ return name(); }
};

/**
 * \brief An abstract parent class for all Binary Logic Conditions.
 *
 * Binary Logic Conditions return the result of performing some
 * Logical operation on a set of conditions. Note that although the
 * name implies the evaluation of two conditions, Binary Logic Conditions
 * can actually evaluate and arbiturary number of conditions.
 */
class BinaryLogicalCondition : public Condition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a BinaryLogicCondition
   *
   * \param conditions The conditions to be evaluated.
   */
  BinaryLogicalCondition(ConditionList& conditions);

  //@}

  /**
   * \brief Deconstructor for a BinaryLogicCondition
   */
  virtual ~BinaryLogicalCondition(){}

  //@}
  
  /** \name Modifier Functions */

  //@{

  /**
   * \brief Adds a Condition to the list of conditions that will
   * be evaluated by this Binary Logical Condition.
   *
   * \param toAdd The condition to be added to the list of
   * conditions this Binary Logic Condition will evaluate.
   */
  void addCondition(RCP<Condition> toAdd);

  //@}

  //! @name Attribute/Query Methods 
  //@{

  /**
   * \brief Applies a Binary Logic operator to two operands and returns the
   * result.
   *
   * \param op1 The first operand.
   * \param op2 The second operand.
   * \return The result of applying a binary logical operator to
   * the two operands.
   */
  virtual bool applyOperator(bool op1, bool op2) const = 0;

  /**
   * \brief Gets a list of all conditions that are a part of this 
   * BinaryLogicalCondition/
   */
  inline
  const ConditionList& getConditions() const{
    return conditions_;
  }

  //@}

  /** \name Overridden from Condition */
  //@{

  /** \brief . */
  virtual bool isConditionTrue() const;

  /** \brief . */
  bool containsAtLeasteOneParameter() const;

  /** \brief . */
  Dependency::ParameterParentMap getAllParameters() const;

  //@}

private:

  /** \name Private Members */
  //@{
  
  /*
   * \brief A list of conditions on which to perform some logic operation.
   */
  ConditionList conditions_;

  //@}

};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<OrCondition>{
public:
  static std::string name(){ return "OrCondition"; }
  static std::string concreteName(const OrCondition&){ return name(); }
};

/**
 * \brief A Binary Logic Condition that returns the result
 * or perfroming a logical OR on the conditions.
 */
class OrCondition : public BinaryLogicalCondition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs an Or Condition
   *
   * @param conditions The conditions to be evaluated.
   */
  OrCondition(ConditionList& conditions);

  /**
   * \brief Deconstructs an Or Condition.
   */
  virtual ~OrCondition(){}

  //@}

  /** \name Overridden from BinaryLogicalCondition */
  //@{

  /** \brief . */
  bool applyOperator(bool op1, bool op2) const;
  
  //@}

};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<AndCondition>{
public:
  static std::string name(){ return "AndCondition"; }
  static std::string concreteName(const AndCondition&){ return name(); }
};

/**
 * \brief A Binary Logic Condition that returns the result
 * or perfroming a logical AND on the conditions.
 */
class AndCondition : public BinaryLogicalCondition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs an And Condition
   *
   * @param conditions The conditions to be evaluated.
   */
  AndCondition(ConditionList& conditions);

  /**
   * \brief Deconstructs an And Condition.
   */
  virtual ~AndCondition(){}
  
  //@}

  /** \name Overridden from BinaryLogicalCondition */
  //@{

  /** \brief . */
  bool applyOperator(bool op1, bool op2) const;
  
  //@}

};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<EqualsCondition>{
public:
  static std::string name(){ return "EqualsCondition"; }
  static std::string concreteName(const EqualsCondition&){ return name(); }
};

/**
 * \brief A Binary Logic Condition that returns the result
 * or perfroming a logical EQUALS on the conditions.
 */
class EqualsCondition : public BinaryLogicalCondition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs an Equals Condition
   *
   * @param conditions The conditions to be evaluated.
   */
  EqualsCondition(ConditionList& conditions);

  /**
   * \brief Deconstructs an Equals Condition.
   */
  virtual ~EqualsCondition(){}
  
  //@}

  /** \name Overridden from Condition */
  //@{

  /** \brief . */
  bool applyOperator(bool op1, bool op2) const;
  
  //@}

};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<NotCondition>{
public:
  static std::string name(){ return "NotCondition"; }
  static std::string concreteName(const NotCondition&){ return name(); }
};

/**
 * \brief A Not condition returns the result of
 * performing a logical NOT on a given
 * condition.
 */
class NotCondition : public Condition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a Not Condition
   *
   * @param condition The condition to be evaluated.
   */
  NotCondition(RCP<Condition> condition);

  /**
   * \brief Deconstructs a Not Condition.
   */
  virtual ~NotCondition(){}
  
  //@}

  /** \name Attribute/Query Functions */
  //@{

  /** \brief Retrieve the child condition */
  RCP<const Condition> getChildCondition(){
    return childCondition_;
  }
  
  //@}

  /** \name Overridden from Condition */
  //@{

  /** \brief . */
  bool isConditionTrue() const;

  /** \brief . */
  bool containsAtLeasteOneParameter() const;

  /** \brief . */
  Dependency::ParameterParentMap getAllParameters() const;

  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * The condition on which to perfrom the logical NOT.
   */
  RCP<Condition> childCondition_;
  
  //@}


};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<ParameterCondition>{
public:
  static std::string name(){ return "ParameterCondition"; }
  static std::string concreteName(const ParameterCondition&){ return name(); }
};

/**
 * \brief An Abstract Base class for all ParameterConditions.
 *
 * A Parmaeter Condition examines the value of a given
 * parameter and returns a bool based on the condition of
 * that value.
 */
class ParameterCondition : public Condition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a Parameter Condition.
   *
   * @param parameterName The name of the parameter to be evaluated.
   * @param parentList The parent Parameter List of the parameter to be evaluated.
   * @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
   * results in a true or when the evaluation results in a false. When set to true, if the parameter
   * evaluates to true then the condition will evaluate to true. If set to false if the parameter
   * evaluates to false, then the condition will evaluate to true.
   */
  ParameterCondition(std::string parameterName, RCP<ParameterList> parentList, bool whenParamEqualsValue);

  virtual ~ParameterCondition(){}
  
  //@}

  //! @name Attribute/Query Methods 
  //@{

  /**
   * Evaluate the current condition of a paramtere and
   * return the result.
   *
   * @param The result of evaluating the current condition
   * of the parameter.
   */
  virtual bool evaluateParameter() const = 0;

  /** \brief Gets a const pointer to the Parameter being
   *  evaluated by this ParameterCondition
   */
  inline const ParameterEntry* getParameter() const{
    return parameter_;
  }

  /** \brief Gets the parameter name */
  inline
  const std::string& getParameterName() const{
    return parameterName_;
  }

  /** \brief Gets the parent parameter list */
  inline
  const RCP<const ParameterList> getParentList() const{
    return parentList_;
  }

  /** \brief Gets the WhenParamEqualsValue */
  inline
  bool getWhenParamEqualsValue() const{
    return whenParamEqualsValue_;
  }
  
  //@}

  /** \name Overridden from Condition */
  //@{

  bool isConditionTrue() const{
    if((whenParamEqualsValue_ && evaluateParameter()) || 
      (!whenParamEqualsValue_ && evaluateParameter()))
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  bool containsAtLeasteOneParameter() const{
    return true;
  }

  Dependency::ParameterParentMap getAllParameters() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * Name of parameter to be evaluated.
   */
  std::string parameterName_;

  /**
   * Parent List of the parameter to be evaluated.
   */
  RCP<ParameterList> parentList_;

  /**
   * Wether or not the condition should evaluate to true if the parameter evaluated to true.
   */
  bool whenParamEqualsValue_;

  /**
   * A pointer to the actual parameter to be evaluated.
   */
  ParameterEntry* parameter_;
  
  //@}

};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<StringCondition>{
public:
  static std::string name(){ return "StringCondition"; }
  static std::string concreteName(const StringCondition&){ return name(); }
};

/**
 * \brief A String Condition is a Parameter Condition that evaluates
 * whether or not a string parameter has taken on a particular
 * value or set of values.
 */
class StringCondition : public ParameterCondition{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Convience typedef representing an array of strings.
   */
  typedef Array<std::string> ValueList; 
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a String Condition.
   *
   * @param parameterName The name of the parameter to be evaluated.
   * @param parentList The parent Parameter List of the parameter to be evaluated.
   * @param value The value to compare the parameter's value against.
   * @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
   * results in a true or when the evaluation results in a false. When set to true, if the parameter
   * evaluates to true then the condition will evaluate to true. If set to false if the parameter
   * evaluates to false, then the condition will evaluate to true.
   */
  StringCondition(std::string parameterName, RCP<ParameterList> parentList, std::string value, bool whenParamEqualsValue=true);

  /**
   * \brief Constructs a String Condition.
   *
   * @param parameterName The name of the parameter to be evaluated.
   * @param parentList The parent Parameter List of the parameter to be evaluated.
   * @param values The values to compare the parameter's value against.
   * @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
   * results in a true or when the evaluation results in a false. When set to true, if the parameter
   * evaluates to true then the condition will evaluate to true. If seet to false if the parameter
   * evaluates to false, then the condition will evaluate to true.
   */
  StringCondition(std::string parameterName, RCP<ParameterList> parentList, ValueList values, bool whenParamEqualsValue=true);

  virtual ~StringCondition(){}
  
  //@}

  /** \name Overridden from ParameterCondition */
  //@{

  bool evaluateParameter() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * A list of values against which to evaluate the parameter's value.
   */
  ValueList values_;
  
  //@}
  
};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<NumberCondition>{
public:
  static std::string name(){ return "NumberCondition"; }
  static std::string concreteName(const NumberCondition&){ return name(); }
};

/**
 * \brief A Number Condition is a Parameter Condition that evaluates
 * whether or not a number parameter is greater 0. 
 *
 * If the parameter is
 * greater than 0 this is interperted as the condition being "true".
 * Otherwise the oncidiont is interperted as false.
 */
template<class T>
class NumberCondition : public ParameterCondition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a Number Condition.
   *
   * @param parameterName The name of the parameter to be evaluated.
   * @param parentList The parent Parameter List of the parameter to be evaluated.
   * @param func A function to run the value of the parameter through. If the function returns a value
   * greater than 0, this will be interperted as the condition being "true". If the 
   * function returns a value of 0 or less, this will be interperted as the condition being false.
   * @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
   * results in a true or when the evaluation results in a false. When set to true, if the parameter
   * evaluates to true then the condition will evaluate to true. If seet to false if the parameter
   * evaluates to false, then the condition will evaluate to true.
   */
  NumberCondition(std::string parameterName, RCP<ParameterList> parentList, T (*func)(T)=0, bool whenParamEqualsValue=true):
    ParameterCondition(parameterName, parentList, whenParamEqualsValue), func_(func)
  {
    checkForNumberType();
  }
    
  /**
   * Constructs a Number Condition.
   *
   * @param parameterName The name of the parameter to be evaluated.
   * @param parentList The parent Parameter List of the parameter to be evaluated.
   * @param func A function to run the value of the parameter through. If the function returns a value
   * greater than 0, this will be interperted as the parameter's current state being "true". If the 
   * function returns a value of 0 or less, this will be interperted as the parameter's current state
   * being "false".
   * @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
   * results in a true or when the evaluation results in a fals. When set to true, if the parameter
   * evaluates to true then the condition will evaluate to true. If seet to false if the parameter
   * evaluates to false, then the condition will evaluate to true.
   */
  NumberCondition(std::string parameterName, RCP<ParameterList> parentList, bool whenParamEqualsValue=true):
    ParameterCondition(parameterName, parentList, whenParamEqualsValue), func_(0)
  {
    checkForNumberType();
  }

  virtual ~NumberCondition(){}
  
  //@}

  /** \name Overridden from ParameterCondition */
  //@{

  bool evaluateParameter() const{
    return (runFunction(getValue<T>(*getParameter())) > 0);
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief . */
  T (*func_)(T);   

  /**
   * \brief Runs the function associated with this condition and
   * returns the result.
   *
   * @param argument The value upon which to run the function.
   */
  inline T runFunction(T argument) const{
    if(func_ !=0)
      return (*func_)(argument);
    else
      return argument;
  }  

  /** \brief Checks to make sure the given parameter is a 
   * number type
   */
  void checkForNumberType() const{
  const ParameterEntry* toCheck = getParameter();
    TEST_FOR_EXCEPTION(
      !toCheck->isType<int>() &&
      !toCheck->isType<short>() &&
      !toCheck->isType<double>() &&
      !toCheck->isType<float>(),
      InvalidConditionException,
      "The parameter of a Number Condition "
      "must be of a supported number type!" << std::endl <<
      "Actual Parameter type: " << getParameter()->getAny().typeName() );
  }
  
  //@}

};

template<>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<BoolCondition>{
public:
  static std::string name(){ return "BoolCondition"; }
  static std::string concreteName(const BoolCondition&){ return name(); }
};

/**
 * \brief A Bool Condition is a Parameter Condition that evaluates
 * whether or not a Boolean parameter is ture.
 * */
class BoolCondition : public ParameterCondition{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a Bool Condition.
   *
   * @param parameterName The name of the parameter to be evaluated.
   * @param parentList The parent Parameter List of the parameter to be evaluated.
   * @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
   * results in a true or when the evaluation results in a false. When set to true, if the parameter
   * evaluates to true then the condition will evaluate to true. If set to false if the parameter
   * evaluates to false, then the condition will evaluate to true.
   */
  BoolCondition(std::string parameterName, RCP<ParameterList> parentList, bool whenParamEqualsValue=true);

  virtual ~BoolCondition(){}
  
  //@}

  /** \name Overridden from ParameterCondition */
  //@{

  bool evaluateParameter() const;
  
  //@}

};


} //namespace Teuchos


#endif //TEUCHOS_STANDARDCONDITION_HPP_
