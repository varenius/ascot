/********************************************************************************
** Form generated from reading UI file 'highlight.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_HIGHLIGHT_H
#define UI_HIGHLIGHT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpinBox>

QT_BEGIN_NAMESPACE

class Ui_Highlight
{
public:
    QGridLayout *gridLayout_2;
    QGroupBox *groupBox;
    QGridLayout *gridLayout;
    QRadioButton *type1_radiobutton;
    QLineEdit *identify_lineedit;
    QLabel *label_2;
    QSpinBox *row1_spinbox;
    QRadioButton *type2_radiobutton;
    QComboBox *relation_combobox;
    QDoubleSpinBox *value_spinbox;
    QLabel *label_7;
    QSpinBox *row2_spinbox;
    QLabel *label;
    QCheckBox *remove_checkbox;
    QCheckBox *allsubs_checkbox;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *Highlight)
    {
        if (Highlight->objectName().isEmpty())
            Highlight->setObjectName(QStringLiteral("Highlight"));
        Highlight->resize(603, 145);
        gridLayout_2 = new QGridLayout(Highlight);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        groupBox = new QGroupBox(Highlight);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout = new QGridLayout(groupBox);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        type1_radiobutton = new QRadioButton(groupBox);
        type1_radiobutton->setObjectName(QStringLiteral("type1_radiobutton"));
        type1_radiobutton->setChecked(true);

        gridLayout->addWidget(type1_radiobutton, 0, 0, 1, 2);

        identify_lineedit = new QLineEdit(groupBox);
        identify_lineedit->setObjectName(QStringLiteral("identify_lineedit"));

        gridLayout->addWidget(identify_lineedit, 0, 2, 1, 2);

        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 0, 4, 1, 2);

        row1_spinbox = new QSpinBox(groupBox);
        row1_spinbox->setObjectName(QStringLiteral("row1_spinbox"));

        gridLayout->addWidget(row1_spinbox, 0, 6, 1, 1);

        type2_radiobutton = new QRadioButton(groupBox);
        type2_radiobutton->setObjectName(QStringLiteral("type2_radiobutton"));
        type2_radiobutton->setCheckable(true);

        gridLayout->addWidget(type2_radiobutton, 1, 0, 1, 1);

        relation_combobox = new QComboBox(groupBox);
        relation_combobox->setObjectName(QStringLiteral("relation_combobox"));

        gridLayout->addWidget(relation_combobox, 1, 1, 1, 2);

        value_spinbox = new QDoubleSpinBox(groupBox);
        value_spinbox->setObjectName(QStringLiteral("value_spinbox"));
        value_spinbox->setDecimals(4);
        value_spinbox->setMaximum(1e+8);

        gridLayout->addWidget(value_spinbox, 1, 3, 1, 2);

        label_7 = new QLabel(groupBox);
        label_7->setObjectName(QStringLiteral("label_7"));

        gridLayout->addWidget(label_7, 1, 5, 1, 2);

        row2_spinbox = new QSpinBox(groupBox);
        row2_spinbox->setObjectName(QStringLiteral("row2_spinbox"));
        row2_spinbox->setValue(1);

        gridLayout->addWidget(row2_spinbox, 1, 7, 1, 1);

        label = new QLabel(groupBox);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 0, 7, 1, 1);


        gridLayout_2->addWidget(groupBox, 0, 0, 1, 3);

        remove_checkbox = new QCheckBox(Highlight);
        remove_checkbox->setObjectName(QStringLiteral("remove_checkbox"));
        remove_checkbox->setChecked(false);

        gridLayout_2->addWidget(remove_checkbox, 1, 0, 1, 1);

        allsubs_checkbox = new QCheckBox(Highlight);
        allsubs_checkbox->setObjectName(QStringLiteral("allsubs_checkbox"));

        gridLayout_2->addWidget(allsubs_checkbox, 1, 1, 1, 1);

        buttonBox = new QDialogButtonBox(Highlight);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_2->addWidget(buttonBox, 1, 2, 1, 1);


        retranslateUi(Highlight);
        QObject::connect(buttonBox, SIGNAL(accepted()), Highlight, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Highlight, SLOT(reject()));

        QMetaObject::connectSlotsByName(Highlight);
    } // setupUi

    void retranslateUi(QDialog *Highlight)
    {
        Highlight->setWindowTitle(QApplication::translate("Highlight", "Dialog", Q_NULLPTR));
        groupBox->setTitle(QApplication::translate("Highlight", "Highlight Points by Tooltips", Q_NULLPTR));
        type1_radiobutton->setText(QApplication::translate("Highlight", "Highlight data containing", Q_NULLPTR));
        label_2->setText(QApplication::translate("Highlight", "in Tooltip-Row", Q_NULLPTR));
        type2_radiobutton->setText(QApplication::translate("Highlight", "Highlight data", Q_NULLPTR));
        relation_combobox->clear();
        relation_combobox->insertItems(0, QStringList()
         << QApplication::translate("Highlight", "bigger than", Q_NULLPTR)
         << QApplication::translate("Highlight", "smaller than", Q_NULLPTR)
         << QApplication::translate("Highlight", "equal to", Q_NULLPTR)
        );
        label_7->setText(QApplication::translate("Highlight", "in Tooltip-Row", Q_NULLPTR));
        label->setText(QApplication::translate("Highlight", " (0 = all rows)", Q_NULLPTR));
        remove_checkbox->setText(QApplication::translate("Highlight", "remove data from old graph", Q_NULLPTR));
        allsubs_checkbox->setText(QApplication::translate("Highlight", "include all subplots", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Highlight: public Ui_Highlight {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_HIGHLIGHT_H
