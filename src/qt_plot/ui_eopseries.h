/********************************************************************************
** Form generated from reading UI file 'eopseries.ui'
**
** Created by: Qt User Interface Compiler version 5.12.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_EOPSERIES_H
#define UI_EOPSERIES_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDateEdit>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QTableWidget>

QT_BEGIN_NAMESPACE

class Ui_EopSeries
{
public:
    QGridLayout *gridLayout_3;
    QLabel *label;
    QTableWidget *eops_tablewidget;
    QGridLayout *gridLayout_2;
    QCheckBox *ut1_zonal_checkbox;
    QCheckBox *pmnut_checkbox;
    QComboBox *nutation_combobox;
    QLabel *label_2;
    QCheckBox *hfocean_checkbox;
    QCheckBox *ut_lib_checkbox;
    QComboBox *interpol_combobox;
    QLabel *label_3;
    QComboBox *pt_version_checkbox;
    QLabel *label_4;
    QGridLayout *gridLayout;
    QLabel *end_label;
    QDateEdit *end_dateedit;
    QDateEdit *start_dateedit;
    QLabel *start_label;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *EopSeries)
    {
        if (EopSeries->objectName().isEmpty())
            EopSeries->setObjectName(QString::fromUtf8("EopSeries"));
        EopSeries->resize(578, 389);
        gridLayout_3 = new QGridLayout(EopSeries);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        label = new QLabel(EopSeries);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout_3->addWidget(label, 0, 0, 1, 1);

        eops_tablewidget = new QTableWidget(EopSeries);
        eops_tablewidget->setObjectName(QString::fromUtf8("eops_tablewidget"));

        gridLayout_3->addWidget(eops_tablewidget, 1, 0, 1, 2);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        ut1_zonal_checkbox = new QCheckBox(EopSeries);
        ut1_zonal_checkbox->setObjectName(QString::fromUtf8("ut1_zonal_checkbox"));

        gridLayout_2->addWidget(ut1_zonal_checkbox, 0, 0, 1, 1);

        pmnut_checkbox = new QCheckBox(EopSeries);
        pmnut_checkbox->setObjectName(QString::fromUtf8("pmnut_checkbox"));

        gridLayout_2->addWidget(pmnut_checkbox, 0, 1, 1, 1);

        nutation_combobox = new QComboBox(EopSeries);
        nutation_combobox->addItem(QString());
        nutation_combobox->addItem(QString());
        nutation_combobox->addItem(QString());
        nutation_combobox->setObjectName(QString::fromUtf8("nutation_combobox"));

        gridLayout_2->addWidget(nutation_combobox, 0, 2, 1, 1);

        label_2 = new QLabel(EopSeries);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_2->addWidget(label_2, 0, 3, 1, 1);

        hfocean_checkbox = new QCheckBox(EopSeries);
        hfocean_checkbox->setObjectName(QString::fromUtf8("hfocean_checkbox"));

        gridLayout_2->addWidget(hfocean_checkbox, 1, 0, 1, 1);

        ut_lib_checkbox = new QCheckBox(EopSeries);
        ut_lib_checkbox->setObjectName(QString::fromUtf8("ut_lib_checkbox"));

        gridLayout_2->addWidget(ut_lib_checkbox, 1, 1, 1, 1);

        interpol_combobox = new QComboBox(EopSeries);
        interpol_combobox->addItem(QString());
        interpol_combobox->addItem(QString());
        interpol_combobox->addItem(QString());
        interpol_combobox->setObjectName(QString::fromUtf8("interpol_combobox"));

        gridLayout_2->addWidget(interpol_combobox, 1, 2, 1, 1);

        label_3 = new QLabel(EopSeries);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_2->addWidget(label_3, 1, 3, 1, 1);

        pt_version_checkbox = new QComboBox(EopSeries);
        pt_version_checkbox->addItem(QString());
        pt_version_checkbox->addItem(QString());
        pt_version_checkbox->addItem(QString());
        pt_version_checkbox->setObjectName(QString::fromUtf8("pt_version_checkbox"));

        gridLayout_2->addWidget(pt_version_checkbox, 2, 2, 1, 1);

        label_4 = new QLabel(EopSeries);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_2->addWidget(label_4, 2, 3, 1, 1);


        gridLayout_3->addLayout(gridLayout_2, 2, 0, 1, 2);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        end_label = new QLabel(EopSeries);
        end_label->setObjectName(QString::fromUtf8("end_label"));

        gridLayout->addWidget(end_label, 1, 0, 1, 1);

        end_dateedit = new QDateEdit(EopSeries);
        end_dateedit->setObjectName(QString::fromUtf8("end_dateedit"));
        end_dateedit->setDateTime(QDateTime(QDate(2015, 1, 1), QTime(0, 0, 0)));

        gridLayout->addWidget(end_dateedit, 1, 1, 1, 1);

        start_dateedit = new QDateEdit(EopSeries);
        start_dateedit->setObjectName(QString::fromUtf8("start_dateedit"));
        start_dateedit->setDateTime(QDateTime(QDate(2080, 1, 1), QTime(0, 0, 0)));
        start_dateedit->setMinimumDateTime(QDateTime(QDate(1752, 9, 16), QTime(0, 0, 0)));

        gridLayout->addWidget(start_dateedit, 0, 1, 1, 1);

        start_label = new QLabel(EopSeries);
        start_label->setObjectName(QString::fromUtf8("start_label"));

        gridLayout->addWidget(start_label, 0, 0, 1, 1);


        gridLayout_3->addLayout(gridLayout, 3, 0, 2, 1);

        buttonBox = new QDialogButtonBox(EopSeries);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_3->addWidget(buttonBox, 4, 1, 1, 1);


        retranslateUi(EopSeries);
        QObject::connect(buttonBox, SIGNAL(accepted()), EopSeries, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), EopSeries, SLOT(reject()));

        QMetaObject::connectSlotsByName(EopSeries);
    } // setupUi

    void retranslateUi(QDialog *EopSeries)
    {
        EopSeries->setWindowTitle(QApplication::translate("EopSeries", "Dialog", nullptr));
        label->setText(QApplication::translate("EopSeries", "Eop Series:", nullptr));
        ut1_zonal_checkbox->setText(QApplication::translate("EopSeries", "UT1 Zonal Tides", nullptr));
        pmnut_checkbox->setText(QApplication::translate("EopSeries", "PM Nutation", nullptr));
        nutation_combobox->setItemText(0, QApplication::translate("EopSeries", "MODFILE", nullptr));
        nutation_combobox->setItemText(1, QApplication::translate("EopSeries", "IAU2000/2006", nullptr));
        nutation_combobox->setItemText(2, QApplication::translate("EopSeries", "FCN", nullptr));

        label_2->setText(QApplication::translate("EopSeries", "Nutation Type", nullptr));
        hfocean_checkbox->setText(QApplication::translate("EopSeries", "HF Ocean", nullptr));
        ut_lib_checkbox->setText(QApplication::translate("EopSeries", "UT Libration", nullptr));
        interpol_combobox->setItemText(0, QApplication::translate("EopSeries", "linear", nullptr));
        interpol_combobox->setItemText(1, QApplication::translate("EopSeries", "cspline", nullptr));
        interpol_combobox->setItemText(2, QApplication::translate("EopSeries", "polynomial", nullptr));

        label_3->setText(QApplication::translate("EopSeries", "Interpolation Type", nullptr));
        pt_version_checkbox->setItemText(0, QApplication::translate("EopSeries", "2015", nullptr));
        pt_version_checkbox->setItemText(1, QApplication::translate("EopSeries", "2010", nullptr));
        pt_version_checkbox->setItemText(2, QApplication::translate("EopSeries", "2003", nullptr));

        label_4->setText(QApplication::translate("EopSeries", "Pole Tide Version", nullptr));
        end_label->setText(QApplication::translate("EopSeries", "End:", nullptr));
        end_dateedit->setDisplayFormat(QApplication::translate("EopSeries", "dd.MMM.yy", nullptr));
        start_dateedit->setDisplayFormat(QApplication::translate("EopSeries", "dd MMM yy", nullptr));
        start_label->setText(QApplication::translate("EopSeries", "Start:", nullptr));
    } // retranslateUi

};

namespace Ui {
    class EopSeries: public Ui_EopSeries {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EOPSERIES_H
