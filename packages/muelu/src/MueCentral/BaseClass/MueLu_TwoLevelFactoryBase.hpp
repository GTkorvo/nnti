#ifndef MUELU_TWOLEVELFACTORY_HPP
#define MUELU_TWOLEVELFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  /*!
    @class Base class for factories (e.g., R, P, and A_coarse).
    @brief Base class for factories that use two levels (fineLevel and coarseLevel).
  */

  class TwoLevelFactoryBase : public FactoryBase {

  public:

    //@{ Constructors/Destructors.

    //! Constructor.
    TwoLevelFactoryBase() {}

    //! Destructor.
    virtual ~TwoLevelFactoryBase() {}

    //@}

    //! Input
    //@{

    virtual void DeclareInput(Level &fineLevel, Level &coarseLevel) const = 0;

    //!
    virtual void CallDeclareInput(Level & requestedLevel) const {
      TEUCHOS_TEST_FOR_EXCEPTION(requestedLevel.GetPreviousLevel() == Teuchos::null, Exceptions::RuntimeError, "LevelID = " << requestedLevel.GetLevelID());
      DeclareInput(*requestedLevel.GetPreviousLevel(), requestedLevel);
    }

    //@}

    //@{

    //! @name Build methods.

    //! Build an object with this factory.
    virtual void Build(Level & fineLevel, Level & coarseLevel) const = 0;

    //!
    virtual void CallBuild(Level & requestedLevel) const {
      TEUCHOS_TEST_FOR_EXCEPTION(requestedLevel.GetPreviousLevel() == Teuchos::null, Exceptions::RuntimeError, "LevelID = " << requestedLevel.GetLevelID());
      Build(*requestedLevel.GetPreviousLevel(), requestedLevel);
    }

    //@}

  }; //class TwoLevelFactoryBase

} //namespace MueLu

#define MUELU_TWOLEVELFACTORY_SHORT
#endif //ifndef MUELU_TWOLEVELFACTORY_HPP
