/********************************************************************************
** Form generated from reading UI file 'BorderWidgetQt.ui'
**
** Created: Sun Jan 6 14:38:51 2013
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BORDERWIDGETQT_H
#define UI_BORDERWIDGETQT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QWidget>
#include "QVTKWidget.h"

QT_BEGIN_NAMESPACE

class Ui_BorderWidgetQt
{
public:
    QAction *actionOpenFile;
    QAction *actionExit;
    QAction *actionPrint;
    QAction *actionHelp;
    QAction *actionSave;
    QAction *open_file_action;
    QAction *add_box_action;
    QAction *delete_one_action;
    QAction *delete_all_action;
    QAction *exit_action;
    QAction *get_tet_action;
    QAction *delete_model_action;
    QAction *rotate_x_action;
    QAction *rotate_y_action;
    QAction *rotate_z_action;
    QAction *save_action;
    QAction *import_action;
    QAction *add_plane_action;
    QAction *get_plane_tet;
    QAction *delete_all_plane_action;
    QAction *add_line_action;
    QAction *blue_action;
    QAction *white_action;
    QAction *red_action;
    QAction *surface_edge_action;
    QAction *export_cell_action;
    QAction *dual_screen_action;
    QAction *single_screen_action;
    QAction *link_camera_action;
    QAction *single_pick_action;
    QAction *multi_pick_action;
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    QVTKWidget *qvtk_widget_left;
    QVTKWidget *qvtk_widget_right;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menu_Edit;
    QMenu *menuRotation;
    QMenu *menu_Line;
    QMenu *menu_Color;
    QMenu *menu_View;
    QMenu *menu_Pick;
    QStatusBar *status_bar;

    void setupUi(QMainWindow *BorderWidgetQt)
    {
        if (BorderWidgetQt->objectName().isEmpty())
            BorderWidgetQt->setObjectName(QString::fromUtf8("BorderWidgetQt"));
        BorderWidgetQt->resize(822, 692);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/Icons/help.png"), QSize(), QIcon::Normal, QIcon::Off);
        BorderWidgetQt->setWindowIcon(icon);
        BorderWidgetQt->setIconSize(QSize(22, 22));
        actionOpenFile = new QAction(BorderWidgetQt);
        actionOpenFile->setObjectName(QString::fromUtf8("actionOpenFile"));
        actionOpenFile->setEnabled(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/Icons/fileopen.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionOpenFile->setIcon(icon1);
        actionExit = new QAction(BorderWidgetQt);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(""), QSize(), QIcon::Normal, QIcon::Off);
        actionExit->setIcon(icon2);
        actionPrint = new QAction(BorderWidgetQt);
        actionPrint->setObjectName(QString::fromUtf8("actionPrint"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/Icons/print.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPrint->setIcon(icon3);
        actionHelp = new QAction(BorderWidgetQt);
        actionHelp->setObjectName(QString::fromUtf8("actionHelp"));
        actionHelp->setIcon(icon);
        actionSave = new QAction(BorderWidgetQt);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/Icons/filesave.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSave->setIcon(icon4);
        open_file_action = new QAction(BorderWidgetQt);
        open_file_action->setObjectName(QString::fromUtf8("open_file_action"));
        add_box_action = new QAction(BorderWidgetQt);
        add_box_action->setObjectName(QString::fromUtf8("add_box_action"));
        delete_one_action = new QAction(BorderWidgetQt);
        delete_one_action->setObjectName(QString::fromUtf8("delete_one_action"));
        delete_all_action = new QAction(BorderWidgetQt);
        delete_all_action->setObjectName(QString::fromUtf8("delete_all_action"));
        exit_action = new QAction(BorderWidgetQt);
        exit_action->setObjectName(QString::fromUtf8("exit_action"));
        get_tet_action = new QAction(BorderWidgetQt);
        get_tet_action->setObjectName(QString::fromUtf8("get_tet_action"));
        delete_model_action = new QAction(BorderWidgetQt);
        delete_model_action->setObjectName(QString::fromUtf8("delete_model_action"));
        rotate_x_action = new QAction(BorderWidgetQt);
        rotate_x_action->setObjectName(QString::fromUtf8("rotate_x_action"));
        rotate_y_action = new QAction(BorderWidgetQt);
        rotate_y_action->setObjectName(QString::fromUtf8("rotate_y_action"));
        rotate_z_action = new QAction(BorderWidgetQt);
        rotate_z_action->setObjectName(QString::fromUtf8("rotate_z_action"));
        save_action = new QAction(BorderWidgetQt);
        save_action->setObjectName(QString::fromUtf8("save_action"));
        import_action = new QAction(BorderWidgetQt);
        import_action->setObjectName(QString::fromUtf8("import_action"));
        add_plane_action = new QAction(BorderWidgetQt);
        add_plane_action->setObjectName(QString::fromUtf8("add_plane_action"));
        get_plane_tet = new QAction(BorderWidgetQt);
        get_plane_tet->setObjectName(QString::fromUtf8("get_plane_tet"));
        delete_all_plane_action = new QAction(BorderWidgetQt);
        delete_all_plane_action->setObjectName(QString::fromUtf8("delete_all_plane_action"));
        add_line_action = new QAction(BorderWidgetQt);
        add_line_action->setObjectName(QString::fromUtf8("add_line_action"));
        blue_action = new QAction(BorderWidgetQt);
        blue_action->setObjectName(QString::fromUtf8("blue_action"));
        white_action = new QAction(BorderWidgetQt);
        white_action->setObjectName(QString::fromUtf8("white_action"));
        red_action = new QAction(BorderWidgetQt);
        red_action->setObjectName(QString::fromUtf8("red_action"));
        surface_edge_action = new QAction(BorderWidgetQt);
        surface_edge_action->setObjectName(QString::fromUtf8("surface_edge_action"));
        surface_edge_action->setEnabled(true);
        export_cell_action = new QAction(BorderWidgetQt);
        export_cell_action->setObjectName(QString::fromUtf8("export_cell_action"));
        dual_screen_action = new QAction(BorderWidgetQt);
        dual_screen_action->setObjectName(QString::fromUtf8("dual_screen_action"));
        single_screen_action = new QAction(BorderWidgetQt);
        single_screen_action->setObjectName(QString::fromUtf8("single_screen_action"));
        link_camera_action = new QAction(BorderWidgetQt);
        link_camera_action->setObjectName(QString::fromUtf8("link_camera_action"));
        single_pick_action = new QAction(BorderWidgetQt);
        single_pick_action->setObjectName(QString::fromUtf8("single_pick_action"));
        multi_pick_action = new QAction(BorderWidgetQt);
        multi_pick_action->setObjectName(QString::fromUtf8("multi_pick_action"));
        centralwidget = new QWidget(BorderWidgetQt);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout = new QGridLayout(centralwidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        qvtk_widget_left = new QVTKWidget(centralwidget);
        qvtk_widget_left->setObjectName(QString::fromUtf8("qvtk_widget_left"));

        gridLayout->addWidget(qvtk_widget_left, 0, 1, 1, 1);

        qvtk_widget_right = new QVTKWidget(centralwidget);
        qvtk_widget_right->setObjectName(QString::fromUtf8("qvtk_widget_right"));

        gridLayout->addWidget(qvtk_widget_right, 0, 2, 1, 1);

        BorderWidgetQt->setCentralWidget(centralwidget);
        menuBar = new QMenuBar(BorderWidgetQt);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 822, 25));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menu_Edit = new QMenu(menuBar);
        menu_Edit->setObjectName(QString::fromUtf8("menu_Edit"));
        menuRotation = new QMenu(menu_Edit);
        menuRotation->setObjectName(QString::fromUtf8("menuRotation"));
        menu_Line = new QMenu(menuBar);
        menu_Line->setObjectName(QString::fromUtf8("menu_Line"));
        menu_Color = new QMenu(menuBar);
        menu_Color->setObjectName(QString::fromUtf8("menu_Color"));
        menu_View = new QMenu(menuBar);
        menu_View->setObjectName(QString::fromUtf8("menu_View"));
        menu_Pick = new QMenu(menuBar);
        menu_Pick->setObjectName(QString::fromUtf8("menu_Pick"));
        BorderWidgetQt->setMenuBar(menuBar);
        status_bar = new QStatusBar(BorderWidgetQt);
        status_bar->setObjectName(QString::fromUtf8("status_bar"));
        BorderWidgetQt->setStatusBar(status_bar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menu_Edit->menuAction());
        menuBar->addAction(menu_Line->menuAction());
        menuBar->addAction(menu_Color->menuAction());
        menuBar->addAction(menu_View->menuAction());
        menuBar->addAction(menu_Pick->menuAction());
        menuFile->addAction(open_file_action);
        menuFile->addAction(save_action);
        menuFile->addAction(import_action);
        menuFile->addAction(export_cell_action);
        menuFile->addSeparator();
        menuFile->addAction(exit_action);
        menu_Edit->addAction(add_box_action);
        menu_Edit->addAction(delete_one_action);
        menu_Edit->addAction(delete_all_action);
        menu_Edit->addAction(get_tet_action);
        menu_Edit->addAction(delete_model_action);
        menu_Edit->addAction(menuRotation->menuAction());
        menu_Edit->addAction(add_plane_action);
        menu_Edit->addAction(get_plane_tet);
        menu_Edit->addAction(delete_all_plane_action);
        menuRotation->addAction(rotate_x_action);
        menuRotation->addAction(rotate_y_action);
        menuRotation->addAction(rotate_z_action);
        menu_Line->addAction(add_line_action);
        menu_Color->addAction(blue_action);
        menu_Color->addAction(white_action);
        menu_Color->addAction(red_action);
        menu_View->addAction(surface_edge_action);
        menu_View->addAction(dual_screen_action);
        menu_View->addAction(single_screen_action);
        menu_View->addSeparator();
        menu_View->addAction(link_camera_action);
        menu_Pick->addAction(single_pick_action);
        menu_Pick->addAction(multi_pick_action);

        retranslateUi(BorderWidgetQt);

        QMetaObject::connectSlotsByName(BorderWidgetQt);
    } // setupUi

    void retranslateUi(QMainWindow *BorderWidgetQt)
    {
        BorderWidgetQt->setWindowTitle(QApplication::translate("BorderWidgetQt", "SimpleView", 0, QApplication::UnicodeUTF8));
        actionOpenFile->setText(QApplication::translate("BorderWidgetQt", "Open File...", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("BorderWidgetQt", "Exit", 0, QApplication::UnicodeUTF8));
        actionPrint->setText(QApplication::translate("BorderWidgetQt", "Print", 0, QApplication::UnicodeUTF8));
        actionHelp->setText(QApplication::translate("BorderWidgetQt", "Help", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("BorderWidgetQt", "Save", 0, QApplication::UnicodeUTF8));
        open_file_action->setText(QApplication::translate("BorderWidgetQt", "&Open", 0, QApplication::UnicodeUTF8));
        open_file_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        add_box_action->setText(QApplication::translate("BorderWidgetQt", "&Add Box", 0, QApplication::UnicodeUTF8));
        add_box_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+A", 0, QApplication::UnicodeUTF8));
        delete_one_action->setText(QApplication::translate("BorderWidgetQt", "&Delete One Box", 0, QApplication::UnicodeUTF8));
        delete_one_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+D", 0, QApplication::UnicodeUTF8));
        delete_all_action->setText(QApplication::translate("BorderWidgetQt", "&Delete All Box", 0, QApplication::UnicodeUTF8));
        delete_all_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+X, A", 0, QApplication::UnicodeUTF8));
        exit_action->setText(QApplication::translate("BorderWidgetQt", "&Exit", 0, QApplication::UnicodeUTF8));
        exit_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+Q", 0, QApplication::UnicodeUTF8));
        get_tet_action->setText(QApplication::translate("BorderWidgetQt", "&Get Tet", 0, QApplication::UnicodeUTF8));
        get_tet_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+G", 0, QApplication::UnicodeUTF8));
        delete_model_action->setText(QApplication::translate("BorderWidgetQt", "&Delete model", 0, QApplication::UnicodeUTF8));
        delete_model_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+X, D", 0, QApplication::UnicodeUTF8));
        rotate_x_action->setText(QApplication::translate("BorderWidgetQt", "Rotate X", 0, QApplication::UnicodeUTF8));
        rotate_x_action->setShortcut(QApplication::translate("BorderWidgetQt", "X", 0, QApplication::UnicodeUTF8));
        rotate_y_action->setText(QApplication::translate("BorderWidgetQt", "Rotate Y", 0, QApplication::UnicodeUTF8));
        rotate_y_action->setShortcut(QApplication::translate("BorderWidgetQt", "Y", 0, QApplication::UnicodeUTF8));
        rotate_z_action->setText(QApplication::translate("BorderWidgetQt", "Rotate Z", 0, QApplication::UnicodeUTF8));
        rotate_z_action->setShortcut(QApplication::translate("BorderWidgetQt", "Z", 0, QApplication::UnicodeUTF8));
        save_action->setText(QApplication::translate("BorderWidgetQt", "&Save", 0, QApplication::UnicodeUTF8));
        save_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+S", 0, QApplication::UnicodeUTF8));
        import_action->setText(QApplication::translate("BorderWidgetQt", "&Import box", 0, QApplication::UnicodeUTF8));
        add_plane_action->setText(QApplication::translate("BorderWidgetQt", "&Add Plane", 0, QApplication::UnicodeUTF8));
        add_plane_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+P", 0, QApplication::UnicodeUTF8));
        get_plane_tet->setText(QApplication::translate("BorderWidgetQt", "&Get Plane Tet", 0, QApplication::UnicodeUTF8));
        delete_all_plane_action->setText(QApplication::translate("BorderWidgetQt", "&Delete All Plane", 0, QApplication::UnicodeUTF8));
        add_line_action->setText(QApplication::translate("BorderWidgetQt", "&Add Line", 0, QApplication::UnicodeUTF8));
        blue_action->setText(QApplication::translate("BorderWidgetQt", "Blue", 0, QApplication::UnicodeUTF8));
        blue_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+0", 0, QApplication::UnicodeUTF8));
        white_action->setText(QApplication::translate("BorderWidgetQt", "White", 0, QApplication::UnicodeUTF8));
        white_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+1", 0, QApplication::UnicodeUTF8));
        red_action->setText(QApplication::translate("BorderWidgetQt", "Red", 0, QApplication::UnicodeUTF8));
        red_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+2", 0, QApplication::UnicodeUTF8));
        surface_edge_action->setText(QApplication::translate("BorderWidgetQt", "&Surface with edge", 0, QApplication::UnicodeUTF8));
        surface_edge_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+Shift+S", 0, QApplication::UnicodeUTF8));
        export_cell_action->setText(QApplication::translate("BorderWidgetQt", "&Export cell", 0, QApplication::UnicodeUTF8));
        export_cell_action->setShortcut(QApplication::translate("BorderWidgetQt", "Ctrl+E", 0, QApplication::UnicodeUTF8));
        dual_screen_action->setText(QApplication::translate("BorderWidgetQt", "&Dual screen", 0, QApplication::UnicodeUTF8));
        single_screen_action->setText(QApplication::translate("BorderWidgetQt", "&Single screen", 0, QApplication::UnicodeUTF8));
        link_camera_action->setText(QApplication::translate("BorderWidgetQt", "&Link camera", 0, QApplication::UnicodeUTF8));
        single_pick_action->setText(QApplication::translate("BorderWidgetQt", "&Single pick", 0, QApplication::UnicodeUTF8));
        multi_pick_action->setText(QApplication::translate("BorderWidgetQt", "&Multi pick", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("BorderWidgetQt", "&File", 0, QApplication::UnicodeUTF8));
        menu_Edit->setTitle(QApplication::translate("BorderWidgetQt", "&Edit", 0, QApplication::UnicodeUTF8));
        menuRotation->setTitle(QApplication::translate("BorderWidgetQt", "Rotation", 0, QApplication::UnicodeUTF8));
        menu_Line->setTitle(QApplication::translate("BorderWidgetQt", "&Line", 0, QApplication::UnicodeUTF8));
        menu_Color->setTitle(QApplication::translate("BorderWidgetQt", "&Color", 0, QApplication::UnicodeUTF8));
        menu_View->setTitle(QApplication::translate("BorderWidgetQt", "&View", 0, QApplication::UnicodeUTF8));
        menu_Pick->setTitle(QApplication::translate("BorderWidgetQt", "&Pick", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class BorderWidgetQt: public Ui_BorderWidgetQt {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BORDERWIDGETQT_H
