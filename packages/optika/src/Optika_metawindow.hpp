// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#ifndef METAWINDOW_HPP_
#define METAWINDOW_HPP_
#include <QMainWindow>
#include <QDialog>
#include <QModelIndex>
#include "Optika_treeview.hpp"

/*! \file Optika_metawindow.hpp
    \brief The main window users interact
    with along with a small search widget used
    in the main window.
*/

class QAction;
class QMenu;
class QLabel;
class QPushButton;
class QLineEdit;
namespace Optika{

/**
 * \brief A small widget that searchs through a parameter list for a particular name
 * of either a parameter or another parameter list.
 */
class SearchWidget : public QDialog{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

	/**
	 * \brief Constructs a SearchWidget.
	 *
	 * @param treeModel The TreeModel being searched.
	 * @param treeView The TreeView being used to display the model.
	 * @param parent The parent widget.
	 */
	SearchWidget(TreeModel *treeModel, TreeView *treeView, QWidget *parent=0);

  //@}

private slots:

  /** \name Private Slots */
  //@{
  
	/**
	 * \brief Searches the for a parameter or parameter list containing the string enterd
	 * in the search terms box.
	 */
	void search();

	/**
	 * \brief Highlights the next result in the list of results that are set
	 * by the search function.
	 */
	void next();

	/**
	 * \brief Highlights the previous result in the list of results that are set
	 * by the search function.
	 */
	void previous();

  //@}

private:

  /** \name Private Functions */
  //@{
  
	/**
	 * \brief Removes any indicies in a QModelIndexList that correspond to a
	 * hidden item.
	 *
	 * @param items A list of indicies from which all hidden items
	 * will be removed
	 * @return A QModelIndexList identical to the one specified in the
	 * items parameter except that all indicies corresponding to hidden
	 * items have been removed.
	 */
	QModelIndexList removeHiddenItems(QModelIndexList& items);

  //@}

  /** \name Private Members */
  //@{
  
	/**
	 * \brief Widgets comprising a search widget
	 */
	QPushButton *searchButton, *closeButton, *nextButton, *previousButton;
	QLineEdit *searchTermsEdit;
	QLabel *matchesLabel;
	TreeModel *treeModel;
	TreeView *treeView;

	/**
	 * \brief The results of the search last performed.
	 */
	QList<QModelIndex> currentSearchResults;

	/**
	 * \brief An iterator over the results of the last search performed.
	 */
	QList<QModelIndex>::const_iterator currentSearchIterator;

  //@}

};

/**
 * \brief The Main Window that contains all other widgets in the Optika GUI.
 * For all undocumented functions please refer to the Qt API.
 */
class MetaWindow : public QMainWindow{
	Q_OBJECT
public:

  /** \name Constructors/Destructor */
  //@{

	/**
	 * \brief Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(RCP<ParameterList> validParameters, 
	QString fileName=QString());

	/**
	 * \brief Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param customFunc The function to run whenever the user clicks the action button.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(RCP<ParameterList> validParameters, 
	void (*customFunc)(RCP<const ParameterList>),
	QString fileName=QString());

	/**
	 * \brief Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(RCP<ParameterList> validParameters, 
	RCP<DependencySheet> dependencySheet,
	QString fileName=QString());

	/**
	 * \brief Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param customFunc The function to run whenever the user clicks the action button.
	 * @param fileName The name of a save file that may store previous values used by a user for the 
	 * Parameter List specified by validParameters.
	 */
	MetaWindow(RCP<ParameterList> validParameters, 
	RCP<DependencySheet> dependencySheet,
	void (*customFunc)(RCP<const ParameterList>),
	QString fileName=QString(),
  const std::string actionButtonText="submit");

	/**
	 * \brief Deconstructer for the metawindow
	 */
	~MetaWindow();

  //@}

  //! @name Getters and Setters
  //@{
  
	/**
	 * \brief Adds the information specified to the about dialog of the GUI.
	 *
	 * @param aboutInfo Information to be added to the about dialog of the GUI.
	 */
	 void setAboutInfo(QString aboutInfo);

	/**
	 * \brief Gets the information to be added to the about dialog of the GUI.
	 *
	 * @return the information to be added to the about dialog of the GUI.
	 */
	 QString getAboutInfo();

   /**
    * \brief Sets the action button text.
    *
    * @param text The text to put in the action button.
    */
   void setActionButtonText(QString newText);

   /**
    * \brief Gets the text being displayed int he action button.
    *
    * @return The text being displayed in the action button.
    */
   QString getActionButtonText();

protected:

  /** \name Overridden from QWidget */
  //@{

	/**
	 * \brief Handles any QCloseEvents for the metawindow.
	 *
	 * @param event The QCloseEvent that was issued.
	 */
	void closeEvent(QCloseEvent *event);

  //@}

private:
  /** \name Private Members */
  //@{
  
	/**
	 * \brief Widgets comprising the MetaWindow
	 */
	SearchWidget *searchWidget;

  /** \brief Various actions. */
	QAction *resetAct, *loadAct, *saveAct, *saveAsAct, *quitAct, *aboutAct, *searchAct;
  /** \brief Various menus. */
	QMenu *fileMenu, *recentMenu, *helpMenu;

  /** \brief The button the user pushes that either closes
   * the MetaWindow or runs the custom function.
   */
  QPushButton *actionButton;

	/**
	 * \brief Any additional about information that should be displayed in the about dialog.
	 */
	QString aboutInfo;

	/**
	 * \brief Load and save directory paths
	 */
	QString currentLoadDir, currentSaveDir;

	/**
	 * \brief A list of recently used documents.
	 */
	QStringList recentDocsList;

	/**
	 * \brief The TreeView being used in the metawindow.
	 */
	TreeView *view;

	/**
	 * \brief The TreeModel being used to display the inputs.
	 */
	TreeModel *model;

	/**
	 * \brief The deleages being used to modify any input values.
	 */
	Delegate *delegate;

  //@}

  /** \name Private Functions */
  //@{
  
	/**
	 * \brief The custom function to run when the user clicks the action button.
	 */
	void (*customFunc)(RCP<const ParameterList>);

  /**
   * \brief Returns the name used to store refernce the settings file.
   */
  static QString getSettingsFileName(){
    static QString settingsFileName("OptikaSettings.xml");
    return settingsFileName;
  }

	/**
	 * \brief Common initialization shared by all constructors.
	 *
	 * @param customFunc The function to run whenever the user clicks the action 
   * button.
   * @param actionButtonText Text to be placed in the action button.
	 */
	void initilization(
    void (*customFunc)(RCP<const ParameterList>)=0, 
    const std::string actionButtonText="submit");

	/**
	 * \brief Creates all the menus for the metawindow.
	 */
	void createMenus();

	/**
	 * \brief Creates all necessary actions used in the menut items.
	 */
	void createActions();

	/**
	 * \brief Loads previous parameter settings.
	 */
	void load();

	/**
	 * \brief Loads the last state of the MetaWindow (things like window size and screen position).
	 */
	void loadLastSettings();

	/**
	 * \brief Saves the state of the MetaWindow (things like window size and screen position).
	 */
	void saveSettings();

	/**
	 * \brief Currently under developement
	 */
	void addRecentDocument(QString recentDocument);

	/**
	 * \brief Currently under developement
	 */
	void updateRecentDocsMenu();

private slots:
  /** \name Private Slots */
  //@{
  
	/**
	 * \brief Resets the treemodel to its default state.
	 */
	void resetModel();

	/**
	 * \brief Saves the parameter list settings to a user specified file.
	 */
	bool saveFileAs();

	/**
	 * \brief Saves the current solver to the file the user has already specified.
	 */
	void saveFile();

	/**
	 * \brief Loads a solver the user was previously working on and had saved.
	 */
	void loadFile();

	/**
	 * \brief Asks the user whether or not they would like to currently save the file they are working on.
	 * Should be used when the user has modified the file and is about to perform an action that would cause those modifiation to be lost.
	 */
	bool saveCurrentUnsavedFile();

	/**
	 * \brief Loads a document from the set of recent documents
	 */
	void loadRecentDoc();

	/**
	 * \brief Shows information about the program.
	 */
	void showAbout();

	/**
	 * \brief Starts a search for a parituclar Parameter or ParameterList.
	 */
	void initiateSearch();
	
	/**
	 * \brief What should happen when the user clicks the action button.
	 */
	void doAction();

  //@}
};



}
#endif /* METAWINDOW_HPP_ */
