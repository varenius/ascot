/********************************************************************************
** Form generated from reading UI file 'movfise.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MOVFISE_H
#define UI_MOVFISE_H

#include <QtCore/QLocale>
#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
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
            Movfise->setObjectName(QStringLiteral("Movfise"));
        Movfise->resize(500, 237);
        Movfise->setMaximumSize(QSize(500, 300));
        textBrowser = new QTextBrowser(Movfise);
        textBrowser->setObjectName(QStringLiteral("textBrowser"));
        textBrowser->setGeometry(QRect(270, 10, 221, 221));
        textBrowser->setStyleSheet(QStringLiteral("background-color: rgba(255, 255, 255, 0);"));
        layoutWidget = new QWidget(Movfise);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(4, 50, 251, 42));
        gridLayout = new QGridLayout(layoutWidget);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(layoutWidget);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 0, 0, 1, 1);

        stepsize_spinbox = new QDoubleSpinBox(layoutWidget);
        stepsize_spinbox->setObjectName(QStringLiteral("stepsize_spinbox"));
        stepsize_spinbox->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom));
        stepsize_spinbox->setDecimals(5);
        stepsize_spinbox->setMaximum(1e+9);
        stepsize_spinbox->setValue(10);

        gridLayout->addWidget(stepsize_spinbox, 0, 1, 1, 1);

        label_3 = new QLabel(layoutWidget);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout->addWidget(label_3, 1, 0, 1, 2);

        layoutWidget1 = new QWidget(Movfise);
        layoutWidget1->setObjectName(QStringLiteral("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(10, 170, 251, 52));
        gridLayout_3 = new QGridLayout(layoutWidget1);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        gridLayout_3->setContentsMargins(0, 0, 0, 0);
        same_checkbox = new QCheckBox(layoutWidget1);
        same_checkbox->setObjectName(QStringLiteral("same_checkbox"));
        same_checkbox->setChecked(true);

        gridLayout_3->addWidget(same_checkbox, 0, 0, 1, 1);

        buttonBox = new QDialogButtonBox(layoutWidget1);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_3->addWidget(buttonBox, 1, 0, 1, 1);

        layoutWidget2 = new QWidget(Movfise);
        layoutWidget2->setObjectName(QStringLiteral("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(4, 12, 251, 24));
        gridLayout_2 = new QGridLayout(layoutWidget2);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        gridLayout_2->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(layoutWidget2);
        label->setObjectName(QStringLiteral("label"));

        gridLayout_2->addWidget(label, 0, 0, 1, 1);

        size_spinbox = new QDoubleSpinBox(layoutWidget2);
        size_spinbox->setObjectName(QStringLiteral("size_spinbox"));
        size_spinbox->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom));
        size_spinbox->setDecimals(5);
        size_spinbox->setMaximum(1e+9);
        size_spinbox->setValue(10);

        gridLayout_2->addWidget(size_spinbox, 0, 1, 1, 1);

        calculatorPlainTextEdit = new QPlainTextEdit(Movfise);
        calculatorPlainTextEdit->setObjectName(QStringLiteral("calculatorPlainTextEdit"));
        calculatorPlainTextEdit->setGeometry(QRect(0, 110, 251, 31));
        calculatorOutput = new QLabel(Movfise);
        calculatorOutput->setObjectName(QStringLiteral("calculatorOutput"));
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
        Movfise->setWindowTitle(QApplication::translate("Movfise", "moving filter settings", Q_NULLPTR));
        label_2->setText(QApplication::translate("Movfise", "Step Size:", Q_NULLPTR));
        label_3->setText(QApplication::translate("Movfise", "Leave 0 for original data steps.", Q_NULLPTR));
        same_checkbox->setText(QApplication::translate("Movfise", "Plot in same window", Q_NULLPTR));
        label->setText(QApplication::translate("Movfise", "Window Size:", Q_NULLPTR));
        calculatorPlainTextEdit->setPlainText(QApplication::translate("Movfise", "Calculator: e.g. type 60*60*24", Q_NULLPTR));
        calculatorOutput->setText(QApplication::translate("Movfise", "= ", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Movfise: public Ui_Movfise {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MOVFISE_H
