#include "side_bar.h"
#include <QFileDialog>
#include <QInputDialog>
#include <fstream>
#include <string>
#include <QtGui>
#include <QTreeView>
#include <iostream>
#include "configure.h"

namespace lq {

using namespace std;
//class MainWindow;

QString select_txt;
void side_bar::init()
{
  root = new QTreeWidgetItem(this->tree_widget,
                             QStringList(QString(MODEL)));
  this->tree_widget->setColumnCount(1);
  this->tree_widget->setHeaderLabels(QStringList()<<TREE_TITLE);
  this->std_spin_box->setRange(STD_ZERO, STD_MAX_RANGE);
  this->std_slider->setRange(STD_ZERO, STD_MAX_RANGE);
  this->std_slider->setSliderPosition(STD_DEFAULT_VALUE);
  this->color_legend_btn->setEnabled(false);
  this->delete_btn_action->setEnabled(false);
  this->std_slider->setEnabled(false);
  this->std_spin_box->setEnabled(false);
}


side_bar::side_bar(QWidget* parent):QDockWidget(parent)
{
  //Just call this to setup the GUI using QT UI file
  setupUi(this);
  create_action();
  init();
}

side_bar::~side_bar()
{

}

void side_bar::add_tree_leaf(const string &file_name)
{
  this->delete_btn_action->setEnabled(true);
  QTreeWidgetItem *leaf = new QTreeWidgetItem(root,
                                              QStringList(QString(file_name.c_str())));
  root->addChild(leaf);
  tree_widget->insertTopLevelItem(0, root);
}

//////////////////////////////////////////////////////////////////
//begin slot
void side_bar::select_item(QTreeWidgetItem *item, int column)
{
  //get deleted file name
  select_txt = item->text(column);
  
}

void side_bar::delete_btn_down()
{
  //send file name which will be deleted to the model view
  emit delete_model(tree_widget->currentItem()->text(0));
  delete tree_widget->currentItem();
}
void side_bar::adjust_color()
{
  emit send_adjust_color(this->std_spin_box->value());
}

void side_bar::show_color_legend()
{
  emit send_color_legend();
}

//end slot
//////////////////////////////////////////////////////////////
void side_bar::set_tree(const vector<string> &model_list)
{
 

  for(size_t i = 0; i < model_list.size(); ++i)
  {
    QTreeWidgetItem *leaf = new QTreeWidgetItem(root,
                                                QStringList(QString(model_list[i].c_str())));
    root->addChild(leaf);
  }
  QList<QTreeWidgetItem *> rootList;
  rootList << root;
  tree_widget->insertTopLevelItems(0, rootList);
}

void side_bar::create_action()
{
  connect(tree_widget, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
           this, SLOT(select_item(QTreeWidgetItem*, int)));

  connect(delete_btn_action, SIGNAL(clicked()), this,
          SLOT(delete_btn_down()));

  connect(std_spin_box, SIGNAL(valueChanged(int)),
          std_slider, SLOT(setValue(int)));

  connect(std_slider, SIGNAL(valueChanged(int)),
          std_spin_box, SLOT(setValue(int)));

  connect(std_spin_box, SIGNAL(valueChanged(int)),
          this, SLOT(adjust_color()));

  connect(color_legend_btn, SIGNAL(clicked()),
          this, SLOT(show_color_legend()));
}

}
