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
#ifndef OPTIKA_DELEGATE_HPP_
#define OPTIKA_DELEGATE_HPP_
#include <QItemDelegate>
#include <QModelIndex>
#include <QSize>
#include "Optika_ArrayWidget.hpp"


namespace Optika{

/**
 * The delegate used for the Optika package. For non-documented functions please refer to the Qt API.
 */
class Delegate : public QItemDelegate{
	Q_OBJECT
public:
	/**
	 * Constructs a Delegate.
	 * 
	 * @param parent The parent object of the Delegate.
	 */
	Delegate(QObject *parent = 0);

	QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
	virtual void setEditorData(QWidget *editor, const QModelIndex &index) const;
	void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;
	virtual void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const;

  /**
   * Gets the value the delegate uses to represent 
   * "true" when constructing comboboxes for
   * boolean parameters.
   *
   * @return The value the delegate uses to represent 
   * "true" when constructing comboboxes for
   * boolean parameters.
   */
  static const QString& getBoolEditorTrue(){
    static const QString boolEditorTrue("true");
    return boolEditorTrue;
  }

  /**
   * Gets the value the delegate uses to represent 
   * "false" when constructing comboboxes for
   * boolean parameters.
   *
   * @return The value the delegate uses to represent 
   * "false" when constructing comboboxes for
   * boolean parameters.
   */
  static const QString& getBoolEditorFalse(){
    static const QString boolEditorFalse("false");
    return boolEditorFalse;
  }


private:
  /**
   * Creates an ArrayWidget to edit the parameter at the given QModelIndex.
   *
   * @param index Index of the parameter to edit.
   * @param type The template type of the array.
   * @param parent The parent widget.
   * @return An array widget with which to edit the parameter.
   */
	QWidget* getArrayEditor(const QModelIndex& index, QString type, QWidget *parent) const;

  /**
   * Sets the data in an array widget to the current values found at index.
   *
   * @param editor ArrayWidget whose values are to be set.
   * @param type The template type of the array.
   * @param index The index of the parameter whose values are going to be
   * put in the ArrayWidget.
   */
  void setArrayWidgetData(QWidget* editor, QString type, const QModelIndex& index) const;

  /**
   * Gets an array from an ArrayWidget and puts it into a QVariant.
   *
   * @param editor The ArrayWidget which has the array we want to put into
   * a QVariant.
   * @param type the template type of the array.
   * @return A Qvariant containing an array that is equal to the array in
   * the given ArrayWidget.
   */
  QVariant extractValueFromArray(QWidget* editor, QString type) const;
};


}
#endif /* OPTIKA_DELEGATE_HPP_ */
