/********************************************************************************
** Form generated from reading UI file 'graphstyle.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GRAPHSTYLE_H
#define UI_GRAPHSTYLE_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
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
            Graphstyle->setObjectName(QStringLiteral("Graphstyle"));
        Graphstyle->resize(645, 415);
        gridLayout_5 = new QGridLayout(Graphstyle);
        gridLayout_5->setObjectName(QStringLiteral("gridLayout_5"));
        buttonBox = new QDialogButtonBox(Graphstyle);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_5->addWidget(buttonBox, 2, 0, 1, 1);

        tabWidget = new QTabWidget(Graphstyle);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        symbolTab = new QWidget();
        symbolTab->setObjectName(QStringLiteral("symbolTab"));
        gridLayout_3 = new QGridLayout(symbolTab);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        groupBox = new QGroupBox(symbolTab);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout_2 = new QGridLayout(groupBox);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        cross = new QRadioButton(groupBox);
        cross->setObjectName(QStringLiteral("cross"));

        gridLayout_2->addWidget(cross, 0, 0, 1, 1);

        plus = new QRadioButton(groupBox);
        plus->setObjectName(QStringLiteral("plus"));

        gridLayout_2->addWidget(plus, 1, 0, 1, 1);

        circle = new QRadioButton(groupBox);
        circle->setObjectName(QStringLiteral("circle"));

        gridLayout_2->addWidget(circle, 2, 0, 1, 1);

        disc = new QRadioButton(groupBox);
        disc->setObjectName(QStringLiteral("disc"));

        gridLayout_2->addWidget(disc, 3, 0, 1, 1);

        square = new QRadioButton(groupBox);
        square->setObjectName(QStringLiteral("square"));

        gridLayout_2->addWidget(square, 4, 0, 1, 1);

        diamond = new QRadioButton(groupBox);
        diamond->setObjectName(QStringLiteral("diamond"));

        gridLayout_2->addWidget(diamond, 5, 0, 1, 1);

        star = new QRadioButton(groupBox);
        star->setObjectName(QStringLiteral("star"));

        gridLayout_2->addWidget(star, 6, 0, 1, 1);

        triangle = new QRadioButton(groupBox);
        triangle->setObjectName(QStringLiteral("triangle"));

        gridLayout_2->addWidget(triangle, 7, 0, 1, 1);

        telescope = new QRadioButton(groupBox);
        telescope->setObjectName(QStringLiteral("telescope"));

        gridLayout_2->addWidget(telescope, 8, 0, 1, 1);


        gridLayout_3->addWidget(groupBox, 0, 0, 2, 1);

        groupBox_2 = new QGroupBox(symbolTab);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        gridLayout_4 = new QGridLayout(groupBox_2);
        gridLayout_4->setObjectName(QStringLiteral("gridLayout_4"));
        dashdot = new QRadioButton(groupBox_2);
        dashdot->setObjectName(QStringLiteral("dashdot"));

        gridLayout_4->addWidget(dashdot, 3, 0, 1, 1);

        dashdotdot = new QRadioButton(groupBox_2);
        dashdotdot->setObjectName(QStringLiteral("dashdotdot"));

        gridLayout_4->addWidget(dashdotdot, 4, 0, 1, 1);

        dash = new QRadioButton(groupBox_2);
        dash->setObjectName(QStringLiteral("dash"));

        gridLayout_4->addWidget(dash, 1, 0, 1, 1);

        dot = new QRadioButton(groupBox_2);
        dot->setObjectName(QStringLiteral("dot"));

        gridLayout_4->addWidget(dot, 2, 0, 1, 1);

        solid = new QRadioButton(groupBox_2);
        solid->setObjectName(QStringLiteral("solid"));

        gridLayout_4->addWidget(solid, 0, 0, 1, 1);


        gridLayout_3->addWidget(groupBox_2, 0, 1, 1, 1);

        groupBox_3 = new QGroupBox(symbolTab);
        groupBox_3->setObjectName(QStringLiteral("groupBox_3"));
        gridLayout_6 = new QGridLayout(groupBox_3);
        gridLayout_6->setObjectName(QStringLiteral("gridLayout_6"));
        stepright = new QRadioButton(groupBox_3);
        stepright->setObjectName(QStringLiteral("stepright"));

        gridLayout_6->addWidget(stepright, 3, 0, 1, 1);

        impulse = new QRadioButton(groupBox_3);
        impulse->setObjectName(QStringLiteral("impulse"));

        gridLayout_6->addWidget(impulse, 4, 0, 1, 1);

        line = new QRadioButton(groupBox_3);
        line->setObjectName(QStringLiteral("line"));

        gridLayout_6->addWidget(line, 1, 0, 1, 1);

        stepleft = new QRadioButton(groupBox_3);
        stepleft->setObjectName(QStringLiteral("stepleft"));

        gridLayout_6->addWidget(stepleft, 2, 0, 1, 1);

        none = new QRadioButton(groupBox_3);
        none->setObjectName(QStringLiteral("none"));

        gridLayout_6->addWidget(none, 0, 0, 1, 1);


        gridLayout_3->addWidget(groupBox_3, 0, 2, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        label = new QLabel(symbolTab);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        markersize = new QDoubleSpinBox(symbolTab);
        markersize->setObjectName(QStringLiteral("markersize"));

        gridLayout->addWidget(markersize, 0, 1, 1, 1);

        label_2 = new QLabel(symbolTab);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 0, 2, 1, 1);

        linewidth = new QDoubleSpinBox(symbolTab);
        linewidth->setObjectName(QStringLiteral("linewidth"));

        gridLayout->addWidget(linewidth, 0, 3, 1, 1);

        label_3 = new QLabel(symbolTab);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout->addWidget(label_3, 1, 0, 1, 1);

        markerwidth = new QDoubleSpinBox(symbolTab);
        markerwidth->setObjectName(QStringLiteral("markerwidth"));

        gridLayout->addWidget(markerwidth, 1, 1, 1, 1);

        filled = new QCheckBox(symbolTab);
        filled->setObjectName(QStringLiteral("filled"));

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
        Graphstyle->setWindowTitle(QApplication::translate("Graphstyle", "Dialog", Q_NULLPTR));
        groupBox->setTitle(QApplication::translate("Graphstyle", "Symbol", Q_NULLPTR));
        cross->setText(QApplication::translate("Graphstyle", "Cross", Q_NULLPTR));
        plus->setText(QApplication::translate("Graphstyle", "Plus", Q_NULLPTR));
        circle->setText(QApplication::translate("Graphstyle", "Circle", Q_NULLPTR));
        disc->setText(QApplication::translate("Graphstyle", "Disc", Q_NULLPTR));
        square->setText(QApplication::translate("Graphstyle", "Square", Q_NULLPTR));
        diamond->setText(QApplication::translate("Graphstyle", "Diamond", Q_NULLPTR));
        star->setText(QApplication::translate("Graphstyle", "Star", Q_NULLPTR));
        triangle->setText(QApplication::translate("Graphstyle", "Triangle", Q_NULLPTR));
        telescope->setText(QApplication::translate("Graphstyle", "Telescope", Q_NULLPTR));
        groupBox_2->setTitle(QApplication::translate("Graphstyle", "Line Pen", Q_NULLPTR));
        dashdot->setText(QApplication::translate("Graphstyle", "DashDot", Q_NULLPTR));
        dashdotdot->setText(QApplication::translate("Graphstyle", "DashDotDot", Q_NULLPTR));
        dash->setText(QApplication::translate("Graphstyle", "Dash", Q_NULLPTR));
        dot->setText(QApplication::translate("Graphstyle", "Dot", Q_NULLPTR));
        solid->setText(QApplication::translate("Graphstyle", "Solid", Q_NULLPTR));
        groupBox_3->setTitle(QApplication::translate("Graphstyle", "Line Style", Q_NULLPTR));
        stepright->setText(QApplication::translate("Graphstyle", "StepRight", Q_NULLPTR));
        impulse->setText(QApplication::translate("Graphstyle", "Impulse", Q_NULLPTR));
        line->setText(QApplication::translate("Graphstyle", "Line", Q_NULLPTR));
        stepleft->setText(QApplication::translate("Graphstyle", "StepLeft", Q_NULLPTR));
        none->setText(QApplication::translate("Graphstyle", "None", Q_NULLPTR));
        label->setText(QApplication::translate("Graphstyle", "Marker Size:", Q_NULLPTR));
        markersize->setPrefix(QString());
        label_2->setText(QApplication::translate("Graphstyle", "Line Width:", Q_NULLPTR));
        linewidth->setPrefix(QString());
        label_3->setText(QApplication::translate("Graphstyle", "Marker Width:", Q_NULLPTR));
        markerwidth->setPrefix(QString());
        filled->setText(QApplication::translate("Graphstyle", "Marker filled", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(symbolTab), QApplication::translate("Graphstyle", "Style", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Graphstyle: public Ui_Graphstyle {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GRAPHSTYLE_H
