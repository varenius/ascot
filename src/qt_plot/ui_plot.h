/********************************************************************************
** Form generated from reading UI file 'plot.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PLOT_H
#define UI_PLOT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_Plot
{
public:
    QGridLayout *gridLayout;
    QCustomPlot *customPlot;

    void setupUi(QWidget *Plot)
    {
        if (Plot->objectName().isEmpty())
            Plot->setObjectName(QStringLiteral("Plot"));
        Plot->resize(484, 384);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(Plot->sizePolicy().hasHeightForWidth());
        Plot->setSizePolicy(sizePolicy);
        gridLayout = new QGridLayout(Plot);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        customPlot = new QCustomPlot(Plot);
        customPlot->setObjectName(QStringLiteral("customPlot"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(customPlot->sizePolicy().hasHeightForWidth());
        customPlot->setSizePolicy(sizePolicy1);

        gridLayout->addWidget(customPlot, 0, 0, 1, 1);


        retranslateUi(Plot);

        QMetaObject::connectSlotsByName(Plot);
    } // setupUi

    void retranslateUi(QWidget *Plot)
    {
        Plot->setWindowTitle(QApplication::translate("Plot", "Form", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Plot: public Ui_Plot {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PLOT_H
