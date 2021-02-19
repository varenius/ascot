/********************************************************************************
** Form generated from reading UI file 'axisrange.ui'
**
** Created by: Qt User Interface Compiler version 5.12.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_AXISRANGE_H
#define UI_AXISRANGE_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Range
{
public:
    QWidget *widget;
    QGridLayout *gridLayout;
    QLabel *label;
    QLineEdit *LE_min;
    QLineEdit *LE_max;
    QLabel *label_2;
    QLineEdit *LE_rot;
    QDialogButtonBox *buttonBox;
    QLabel *label_3;

    void setupUi(QDialog *Range)
    {
        if (Range->objectName().isEmpty())
            Range->setObjectName(QString::fromUtf8("Range"));
        Range->resize(361, 139);
        widget = new QWidget(Range);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(10, 29, 333, 91));
        gridLayout = new QGridLayout(widget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(widget);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        LE_min = new QLineEdit(widget);
        LE_min->setObjectName(QString::fromUtf8("LE_min"));

        gridLayout->addWidget(LE_min, 0, 1, 1, 1);

        LE_max = new QLineEdit(widget);
        LE_max->setObjectName(QString::fromUtf8("LE_max"));

        gridLayout->addWidget(LE_max, 0, 2, 1, 1);

        label_2 = new QLabel(widget);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        LE_rot = new QLineEdit(widget);
        LE_rot->setObjectName(QString::fromUtf8("LE_rot"));

        gridLayout->addWidget(LE_rot, 1, 1, 1, 1);

        buttonBox = new QDialogButtonBox(widget);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout->addWidget(buttonBox, 2, 1, 1, 2);

        label_3 = new QLabel(widget);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 1, 2, 1, 1);


        retranslateUi(Range);
        QObject::connect(buttonBox, SIGNAL(accepted()), Range, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Range, SLOT(reject()));

        QMetaObject::connectSlotsByName(Range);
    } // setupUi

    void retranslateUi(QDialog *Range)
    {
        Range->setWindowTitle(QApplication::translate("Range", "Dialog", nullptr));
        label->setText(QApplication::translate("Range", "Min / Max:", nullptr));
        label_2->setText(QApplication::translate("Range", "Rotation:", nullptr));
        label_3->setText(QApplication::translate("Range", "(-90\302\260 to 90\302\260)", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Range: public Ui_Range {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_AXISRANGE_H
