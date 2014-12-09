/********************************************************************************
** Form generated from reading UI file 'side_bar.ui'
**
** Created: Mon Dec 17 21:13:50 2012
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SIDE_BAR_H
#define UI_SIDE_BAR_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDockWidget>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QScrollArea>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QTreeWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_side_bar
{
public:
    QWidget *dockWidgetContents;
    QGridLayout *gridLayout_3;
    QTreeWidget *tree_widget;
    QScrollArea *scrollArea;
    QWidget *scrollAreaWidgetContents;
    QGridLayout *gridLayout_2;
    QVBoxLayout *verticalLayout_2;
    QGridLayout *gridLayout;
    QLabel *label;
    QPushButton *apply_btn;
    QSpacerItem *verticalSpacer_2;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *verticalSpacer;
    QSpacerItem *verticalSpacer_4;
    QSpinBox *std_spin_box;
    QSpacerItem *horizontalSpacer_2;
    QSpacerItem *verticalSpacer_3;
    QSlider *std_slider;
    QPushButton *delete_btn_action;
    QVBoxLayout *verticalLayout;
    QSpacerItem *verticalSpacer_5;
    QPushButton *color_legend_btn;

    void setupUi(QDockWidget *side_bar)
    {
        if (side_bar->objectName().isEmpty())
            side_bar->setObjectName(QString::fromUtf8("side_bar"));
        side_bar->resize(250, 553);
        side_bar->setFloating(false);
        side_bar->setFeatures(QDockWidget::NoDockWidgetFeatures);
        side_bar->setAllowedAreas(Qt::LeftDockWidgetArea|Qt::RightDockWidgetArea);
        dockWidgetContents = new QWidget();
        dockWidgetContents->setObjectName(QString::fromUtf8("dockWidgetContents"));
        gridLayout_3 = new QGridLayout(dockWidgetContents);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        tree_widget = new QTreeWidget(dockWidgetContents);
        QTreeWidgetItem *__qtreewidgetitem = new QTreeWidgetItem();
        __qtreewidgetitem->setText(0, QString::fromUtf8("1"));
        tree_widget->setHeaderItem(__qtreewidgetitem);
        tree_widget->setObjectName(QString::fromUtf8("tree_widget"));

        gridLayout_3->addWidget(tree_widget, 0, 0, 1, 1);

        scrollArea = new QScrollArea(dockWidgetContents);
        scrollArea->setObjectName(QString::fromUtf8("scrollArea"));
        scrollArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 230, 250));
        gridLayout_2 = new QGridLayout(scrollAreaWidgetContents);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(scrollAreaWidgetContents);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 3, 0, 1, 1);

        apply_btn = new QPushButton(scrollAreaWidgetContents);
        apply_btn->setObjectName(QString::fromUtf8("apply_btn"));

        gridLayout->addWidget(apply_btn, 1, 4, 1, 2);

        verticalSpacer_2 = new QSpacerItem(20, 18, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer_2, 0, 0, 1, 1);

        horizontalSpacer = new QSpacerItem(28, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer, 1, 2, 1, 2);

        verticalSpacer = new QSpacerItem(20, 48, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer, 2, 0, 1, 1);

        verticalSpacer_4 = new QSpacerItem(20, 48, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer_4, 2, 5, 1, 1);

        std_spin_box = new QSpinBox(scrollAreaWidgetContents);
        std_spin_box->setObjectName(QString::fromUtf8("std_spin_box"));

        gridLayout->addWidget(std_spin_box, 3, 3, 1, 2);

        horizontalSpacer_2 = new QSpacerItem(28, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 3, 1, 1, 2);

        verticalSpacer_3 = new QSpacerItem(20, 18, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer_3, 0, 5, 1, 1);

        std_slider = new QSlider(scrollAreaWidgetContents);
        std_slider->setObjectName(QString::fromUtf8("std_slider"));
        std_slider->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(std_slider, 4, 0, 1, 6);

        delete_btn_action = new QPushButton(scrollAreaWidgetContents);
        delete_btn_action->setObjectName(QString::fromUtf8("delete_btn_action"));
        delete_btn_action->setMinimumSize(QSize(85, 27));
        delete_btn_action->setMaximumSize(QSize(85, 16777215));

        gridLayout->addWidget(delete_btn_action, 1, 0, 1, 2);


        verticalLayout_2->addLayout(gridLayout);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalSpacer_5 = new QSpacerItem(20, 28, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer_5);

        color_legend_btn = new QPushButton(scrollAreaWidgetContents);
        color_legend_btn->setObjectName(QString::fromUtf8("color_legend_btn"));

        verticalLayout->addWidget(color_legend_btn);


        verticalLayout_2->addLayout(verticalLayout);


        gridLayout_2->addLayout(verticalLayout_2, 0, 0, 1, 1);

        scrollArea->setWidget(scrollAreaWidgetContents);

        gridLayout_3->addWidget(scrollArea, 1, 0, 1, 1);

        side_bar->setWidget(dockWidgetContents);

        retranslateUi(side_bar);
        QObject::connect(std_spin_box, SIGNAL(valueChanged(int)), std_slider, SLOT(setValue(int)));

        QMetaObject::connectSlotsByName(side_bar);
    } // setupUi

    void retranslateUi(QDockWidget *side_bar)
    {
        side_bar->setWindowTitle(QApplication::translate("side_bar", "contral panel", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("side_bar", "Standard:", 0, QApplication::UnicodeUTF8));
        apply_btn->setText(QApplication::translate("side_bar", "Apply", 0, QApplication::UnicodeUTF8));
        delete_btn_action->setText(QApplication::translate("side_bar", "Delete", 0, QApplication::UnicodeUTF8));
        color_legend_btn->setText(QApplication::translate("side_bar", "show color legend", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class side_bar: public Ui_side_bar {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SIDE_BAR_H
