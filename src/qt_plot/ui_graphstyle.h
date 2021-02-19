/********************************************************************************
** Form generated from reading UI file 'graphstyle.ui'
**
** Created by: Qt User Interface Compiler version 5.12.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GRAPHSTYLE_H
#define UI_GRAPHSTYLE_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Graphstyle
{
public:
    QGridLayout *gridLayout_5;
    QDialogButtonBox *buttonBox;
    QTabWidget *tabWidget;
    QWidget *symbolTab;
    QGridLayout *gridLayout_3;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_2;
    QRadioButton *cross;
    QRadioButton *plus;
    QRadioButton *circle;
    QRadioButton *disc;
    QRadioButton *square;
    QRadioButton *diamond;
    QRadioButton *star;
    QRadioButton *triangle;
    QRadioButton *telescope;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_4;
    QRadioButton *dashdot;
    QRadioButton *dashdotdot;
    QRadioButton *dash;
    QRadioButton *dot;
    QRadioButton *solid;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_6;
    QRadioButton *stepright;
    QRadioButton *impulse;
    QRadioButton *line;
    QRadioButton *stepleft;
    QRadioButton *none;
    QGridLayout *gridLayout;
    QLabel *label;
    QDoubleSpinBox *markersize;
    QLabel *label_2;
    QDoubleSpinBox *linewidth;
    QLabel *label_3;
    QDoubleSpinBox *markerwidth;
    QCheckBox *filled;

    void setupUi(QDialog *Graphstyle)
    {
        if (Graphstyle->objectName().isEmpty())
            Graphstyle->setObjectName(QString::fromUtf8("Graphstyle"));
        Graphstyle->resize(645, 415);
        gridLayout_5 = new QGridLayout(Graphstyle);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        buttonBox = new QDialogButtonBox(Graphstyle);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_5->addWidget(buttonBox, 2, 0, 1, 1);

        tabWidget = new QTabWidget(Graphstyle);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        symbolTab = new QWidget();
        symbolTab->setObjectName(QString::fromUtf8("symbolTab"));
        gridLayout_3 = new QGridLayout(symbolTab);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        groupBox = new QGroupBox(symbolTab);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        gridLayout_2 = new QGridLayout(groupBox);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        cross = new QRadioButton(groupBox);
        cross->setObjectName(QString::fromUtf8("cross"));

        gridLayout_2->addWidget(cross, 0, 0, 1, 1);

        plus = new QRadioButton(groupBox);
        plus->setObjectName(QString::fromUtf8("plus"));

        gridLayout_2->addWidget(plus, 1, 0, 1, 1);

        circle = new QRadioButton(groupBox);
        circle->setObjectName(QString::fromUtf8("circle"));

        gridLayout_2->addWidget(circle, 2, 0, 1, 1);

        disc = new QRadioButton(groupBox);
        disc->setObjectName(QString::fromUtf8("disc"));

        gridLayout_2->addWidget(disc, 3, 0, 1, 1);

        square = new QRadioButton(groupBox);
        square->setObjectName(QString::fromUtf8("square"));

        gridLayout_2->addWidget(square, 4, 0, 1, 1);

        diamond = new QRadioButton(groupBox);
        diamond->setObjectName(QString::fromUtf8("diamond"));

        gridLayout_2->addWidget(diamond, 5, 0, 1, 1);

        star = new QRadioButton(groupBox);
        star->setObjectName(QString::fromUtf8("star"));

        gridLayout_2->addWidget(star, 6, 0, 1, 1);

        triangle = new QRadioButton(groupBox);
        triangle->setObjectName(QString::fromUtf8("triangle"));

        gridLayout_2->addWidget(triangle, 7, 0, 1, 1);

        telescope = new QRadioButton(groupBox);
        telescope->setObjectName(QString::fromUtf8("telescope"));

        gridLayout_2->addWidget(telescope, 8, 0, 1, 1);


        gridLayout_3->addWidget(groupBox, 0, 0, 2, 1);

        groupBox_2 = new QGroupBox(symbolTab);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        gridLayout_4 = new QGridLayout(groupBox_2);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        dashdot = new QRadioButton(groupBox_2);
        dashdot->setObjectName(QString::fromUtf8("dashdot"));

        gridLayout_4->addWidget(dashdot, 3, 0, 1, 1);

        dashdotdot = new QRadioButton(groupBox_2);
        dashdotdot->setObjectName(QString::fromUtf8("dashdotdot"));

        gridLayout_4->addWidget(dashdotdot, 4, 0, 1, 1);

        dash = new QRadioButton(groupBox_2);
        dash->setObjectName(QString::fromUtf8("dash"));

        gridLayout_4->addWidget(dash, 1, 0, 1, 1);

        dot = new QRadioButton(groupBox_2);
        dot->setObjectName(QString::fromUtf8("dot"));

        gridLayout_4->addWidget(dot, 2, 0, 1, 1);

        solid = new QRadioButton(groupBox_2);
        solid->setObjectName(QString::fromUtf8("solid"));

        gridLayout_4->addWidget(solid, 0, 0, 1, 1);


        gridLayout_3->addWidget(groupBox_2, 0, 1, 1, 1);

        groupBox_3 = new QGroupBox(symbolTab);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        gridLayout_6 = new QGridLayout(groupBox_3);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        stepright = new QRadioButton(groupBox_3);
        stepright->setObjectName(QString::fromUtf8("stepright"));

        gridLayout_6->addWidget(stepright, 3, 0, 1, 1);

        impulse = new QRadioButton(groupBox_3);
        impulse->setObjectName(QString::fromUtf8("impulse"));

        gridLayout_6->addWidget(impulse, 4, 0, 1, 1);

        line = new QRadioButton(groupBox_3);
        line->setObjectName(QString::fromUtf8("line"));

        gridLayout_6->addWidget(line, 1, 0, 1, 1);

        stepleft = new QRadioButton(groupBox_3);
        stepleft->setObjectName(QString::fromUtf8("stepleft"));

        gridLayout_6->addWidget(stepleft, 2, 0, 1, 1);

        none = new QRadioButton(groupBox_3);
        none->setObjectName(QString::fromUtf8("none"));

        gridLayout_6->addWidget(none, 0, 0, 1, 1);


        gridLayout_3->addWidget(groupBox_3, 0, 2, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(symbolTab);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        markersize = new QDoubleSpinBox(symbolTab);
        markersize->setObjectName(QString::fromUtf8("markersize"));

        gridLayout->addWidget(markersize, 0, 1, 1, 1);

        label_2 = new QLabel(symbolTab);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 0, 2, 1, 1);

        linewidth = new QDoubleSpinBox(symbolTab);
        linewidth->setObjectName(QString::fromUtf8("linewidth"));

        gridLayout->addWidget(linewidth, 0, 3, 1, 1);

        label_3 = new QLabel(symbolTab);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 1, 0, 1, 1);

        markerwidth = new QDoubleSpinBox(symbolTab);
        markerwidth->setObjectName(QString::fromUtf8("markerwidth"));

        gridLayout->addWidget(markerwidth, 1, 1, 1, 1);

        filled = new QCheckBox(symbolTab);
        filled->setObjectName(QString::fromUtf8("filled"));

        gridLayout->addWidget(filled, 1, 2, 1, 2);


        gridLayout_3->addLayout(gridLayout, 1, 1, 1, 2);

        tabWidget->addTab(symbolTab, QString());

        gridLayout_5->addWidget(tabWidget, 0, 0, 1, 1);


        retranslateUi(Graphstyle);
        QObject::connect(buttonBox, SIGNAL(accepted()), Graphstyle, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Graphstyle, SLOT(reject()));

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(Graphstyle);
    } // setupUi

    void retranslateUi(QDialog *Graphstyle)
    {
        Graphstyle->setWindowTitle(QApplication::translate("Graphstyle", "Dialog", nullptr));
        groupBox->setTitle(QApplication::translate("Graphstyle", "Symbol", nullptr));
        cross->setText(QApplication::translate("Graphstyle", "Cross", nullptr));
        plus->setText(QApplication::translate("Graphstyle", "Plus", nullptr));
        circle->setText(QApplication::translate("Graphstyle", "Circle", nullptr));
        disc->setText(QApplication::translate("Graphstyle", "Disc", nullptr));
        square->setText(QApplication::translate("Graphstyle", "Square", nullptr));
        diamond->setText(QApplication::translate("Graphstyle", "Diamond", nullptr));
        star->setText(QApplication::translate("Graphstyle", "Star", nullptr));
        triangle->setText(QApplication::translate("Graphstyle", "Triangle", nullptr));
        telescope->setText(QApplication::translate("Graphstyle", "Telescope", nullptr));
        groupBox_2->setTitle(QApplication::translate("Graphstyle", "Line Pen", nullptr));
        dashdot->setText(QApplication::translate("Graphstyle", "DashDot", nullptr));
        dashdotdot->setText(QApplication::translate("Graphstyle", "DashDotDot", nullptr));
        dash->setText(QApplication::translate("Graphstyle", "Dash", nullptr));
        dot->setText(QApplication::translate("Graphstyle", "Dot", nullptr));
        solid->setText(QApplication::translate("Graphstyle", "Solid", nullptr));
        groupBox_3->setTitle(QApplication::translate("Graphstyle", "Line Style", nullptr));
        stepright->setText(QApplication::translate("Graphstyle", "StepRight", nullptr));
        impulse->setText(QApplication::translate("Graphstyle", "Impulse", nullptr));
        line->setText(QApplication::translate("Graphstyle", "Line", nullptr));
        stepleft->setText(QApplication::translate("Graphstyle", "StepLeft", nullptr));
        none->setText(QApplication::translate("Graphstyle", "None", nullptr));
        label->setText(QApplication::translate("Graphstyle", "Marker Size:", nullptr));
        markersize->setPrefix(QString());
        label_2->setText(QApplication::translate("Graphstyle", "Line Width:", nullptr));
        linewidth->setPrefix(QString());
        label_3->setText(QApplication::translate("Graphstyle", "Marker Width:", nullptr));
        markerwidth->setPrefix(QString());
        filled->setText(QApplication::translate("Graphstyle", "Marker filled", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(symbolTab), QApplication::translate("Graphstyle", "Style", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Graphstyle: public Ui_Graphstyle {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GRAPHSTYLE_H
