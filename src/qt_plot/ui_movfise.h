/********************************************************************************
** Form generated from reading UI file 'movfise.ui'
**
** Created by: Qt User Interface Compiler version 5.12.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MOVFISE_H
#define UI_MOVFISE_H

#include <QtCore/QLocale>
#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Movfise
{
public:
    QTextBrowser *textBrowser;
    QWidget *layoutWidget;
    QGridLayout *gridLayout;
    QLabel *label_2;
    QDoubleSpinBox *stepsize_spinbox;
    QLabel *label_3;
    QWidget *layoutWidget1;
    QGridLayout *gridLayout_3;
    QCheckBox *same_checkbox;
    QDialogButtonBox *buttonBox;
    QWidget *layoutWidget2;
    QGridLayout *gridLayout_2;
    QLabel *label;
    QDoubleSpinBox *size_spinbox;
    QPlainTextEdit *calculatorPlainTextEdit;
    QLabel *calculatorOutput;

    void setupUi(QDialog *Movfise)
    {
        if (Movfise->objectName().isEmpty())
            Movfise->setObjectName(QString::fromUtf8("Movfise"));
        Movfise->resize(500, 237);
        Movfise->setMaximumSize(QSize(500, 300));
        textBrowser = new QTextBrowser(Movfise);
        textBrowser->setObjectName(QString::fromUtf8("textBrowser"));
        textBrowser->setGeometry(QRect(270, 10, 221, 221));
        textBrowser->setStyleSheet(QString::fromUtf8("background-color: rgba(255, 255, 255, 0);"));
        layoutWidget = new QWidget(Movfise);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(4, 50, 251, 42));
        gridLayout = new QGridLayout(layoutWidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(layoutWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 0, 0, 1, 1);

        stepsize_spinbox = new QDoubleSpinBox(layoutWidget);
        stepsize_spinbox->setObjectName(QString::fromUtf8("stepsize_spinbox"));
        stepsize_spinbox->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom));
        stepsize_spinbox->setDecimals(5);
        stepsize_spinbox->setMaximum(999999999.000000000000000);
        stepsize_spinbox->setValue(10.000000000000000);

        gridLayout->addWidget(stepsize_spinbox, 0, 1, 1, 1);

        label_3 = new QLabel(layoutWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 1, 0, 1, 2);

        layoutWidget1 = new QWidget(Movfise);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(10, 170, 251, 52));
        gridLayout_3 = new QGridLayout(layoutWidget1);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setContentsMargins(0, 0, 0, 0);
        same_checkbox = new QCheckBox(layoutWidget1);
        same_checkbox->setObjectName(QString::fromUtf8("same_checkbox"));
        same_checkbox->setChecked(true);

        gridLayout_3->addWidget(same_checkbox, 0, 0, 1, 1);

        buttonBox = new QDialogButtonBox(layoutWidget1);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_3->addWidget(buttonBox, 1, 0, 1, 1);

        layoutWidget2 = new QWidget(Movfise);
        layoutWidget2->setObjectName(QString::fromUtf8("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(4, 12, 251, 24));
        gridLayout_2 = new QGridLayout(layoutWidget2);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(layoutWidget2);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout_2->addWidget(label, 0, 0, 1, 1);

        size_spinbox = new QDoubleSpinBox(layoutWidget2);
        size_spinbox->setObjectName(QString::fromUtf8("size_spinbox"));
        size_spinbox->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom));
        size_spinbox->setDecimals(5);
        size_spinbox->setMaximum(999999999.990000009536743);
        size_spinbox->setValue(10.000000000000000);

        gridLayout_2->addWidget(size_spinbox, 0, 1, 1, 1);

        calculatorPlainTextEdit = new QPlainTextEdit(Movfise);
        calculatorPlainTextEdit->setObjectName(QString::fromUtf8("calculatorPlainTextEdit"));
        calculatorPlainTextEdit->setGeometry(QRect(0, 110, 251, 31));
        calculatorOutput = new QLabel(Movfise);
        calculatorOutput->setObjectName(QString::fromUtf8("calculatorOutput"));
        calculatorOutput->setGeometry(QRect(10, 140, 241, 21));
        layoutWidget->raise();
        layoutWidget->raise();
        layoutWidget->raise();
        textBrowser->raise();
        calculatorPlainTextEdit->raise();
        calculatorOutput->raise();

        retranslateUi(Movfise);
        QObject::connect(buttonBox, SIGNAL(accepted()), Movfise, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Movfise, SLOT(reject()));

        QMetaObject::connectSlotsByName(Movfise);
    } // setupUi

    void retranslateUi(QDialog *Movfise)
    {
        Movfise->setWindowTitle(QApplication::translate("Movfise", "moving filter settings", nullptr));
        label_2->setText(QApplication::translate("Movfise", "Step Size:", nullptr));
        label_3->setText(QApplication::translate("Movfise", "Leave 0 for original data steps.", nullptr));
        same_checkbox->setText(QApplication::translate("Movfise", "Plot in same window", nullptr));
        label->setText(QApplication::translate("Movfise", "Window Size:", nullptr));
        calculatorPlainTextEdit->setPlainText(QApplication::translate("Movfise", "Calculator: e.g. type 60*60*24", nullptr));
        calculatorOutput->setText(QApplication::translate("Movfise", "= ", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Movfise: public Ui_Movfise {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MOVFISE_H
