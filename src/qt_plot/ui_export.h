/********************************************************************************
** Form generated from reading UI file 'export.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_EXPORT_H
#define UI_EXPORT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>

QT_BEGIN_NAMESPACE

class Ui_Export
{
public:
    QGridLayout *gridLayout_3;
    QGroupBox *GB_export;
    QFormLayout *formLayout;
    QGridLayout *gridLayout;
    QRadioButton *RB_pdf;
    QRadioButton *RB_png;
    QRadioButton *RB_bmp;
    QRadioButton *RB_jpg;
    QGridLayout *gridLayout_2;
    QLabel *label;
    QLineEdit *LE_pathfile;
    QPushButton *PB_file;
    QLabel *label_2;
    QCheckBox *CB_cospen;
    QLabel *label_3;
    QSpinBox *SB_width;
    QSpinBox *SB_height;
    QLabel *label_6;
    QSpinBox *SB_scale;
    QSpacerItem *horizontalSpacer;
    QLabel *label_7;
    QSpinBox *SB_quality;
    QSpacerItem *horizontalSpacer_2;
    QLabel *label_5;
    QLineEdit *LE_pdftitle;
    QDialogButtonBox *BB_exit;

    void setupUi(QDialog *Export)
    {
        if (Export->objectName().isEmpty())
            Export->setObjectName(QStringLiteral("Export"));
        Export->resize(449, 255);
        gridLayout_3 = new QGridLayout(Export);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        GB_export = new QGroupBox(Export);
        GB_export->setObjectName(QStringLiteral("GB_export"));
        formLayout = new QFormLayout(GB_export);
        formLayout->setObjectName(QStringLiteral("formLayout"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        RB_pdf = new QRadioButton(GB_export);
        RB_pdf->setObjectName(QStringLiteral("RB_pdf"));
        RB_pdf->setChecked(true);

        gridLayout->addWidget(RB_pdf, 0, 0, 1, 1);

        RB_png = new QRadioButton(GB_export);
        RB_png->setObjectName(QStringLiteral("RB_png"));

        gridLayout->addWidget(RB_png, 1, 0, 1, 1);

        RB_bmp = new QRadioButton(GB_export);
        RB_bmp->setObjectName(QStringLiteral("RB_bmp"));

        gridLayout->addWidget(RB_bmp, 2, 0, 1, 1);

        RB_jpg = new QRadioButton(GB_export);
        RB_jpg->setObjectName(QStringLiteral("RB_jpg"));

        gridLayout->addWidget(RB_jpg, 3, 0, 1, 1);


        formLayout->setLayout(0, QFormLayout::LabelRole, gridLayout);


        gridLayout_3->addWidget(GB_export, 0, 0, 1, 1);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        label = new QLabel(Export);
        label->setObjectName(QStringLiteral("label"));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label, 0, 0, 1, 1);

        LE_pathfile = new QLineEdit(Export);
        LE_pathfile->setObjectName(QStringLiteral("LE_pathfile"));

        gridLayout_2->addWidget(LE_pathfile, 0, 1, 1, 2);

        PB_file = new QPushButton(Export);
        PB_file->setObjectName(QStringLiteral("PB_file"));

        gridLayout_2->addWidget(PB_file, 0, 3, 1, 1);

        label_2 = new QLabel(Export);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_2, 1, 0, 1, 1);

        CB_cospen = new QCheckBox(Export);
        CB_cospen->setObjectName(QStringLiteral("CB_cospen"));

        gridLayout_2->addWidget(CB_cospen, 1, 1, 1, 1);

        label_3 = new QLabel(Export);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_3, 2, 0, 1, 1);

        SB_width = new QSpinBox(Export);
        SB_width->setObjectName(QStringLiteral("SB_width"));
        SB_width->setMaximum(5000);
        SB_width->setSingleStep(50);

        gridLayout_2->addWidget(SB_width, 2, 1, 1, 1);

        SB_height = new QSpinBox(Export);
        SB_height->setObjectName(QStringLiteral("SB_height"));
        SB_height->setMaximum(5000);
        SB_height->setSingleStep(50);

        gridLayout_2->addWidget(SB_height, 2, 2, 1, 2);

        label_6 = new QLabel(Export);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setLayoutDirection(Qt::LeftToRight);
        label_6->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_6, 3, 0, 1, 1);

        SB_scale = new QSpinBox(Export);
        SB_scale->setObjectName(QStringLiteral("SB_scale"));
        SB_scale->setReadOnly(true);
        SB_scale->setMinimum(-1000);
        SB_scale->setMaximum(1000);
        SB_scale->setSingleStep(1);
        SB_scale->setValue(1);

        gridLayout_2->addWidget(SB_scale, 3, 1, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_2->addItem(horizontalSpacer, 3, 2, 1, 1);

        label_7 = new QLabel(Export);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_7, 4, 0, 1, 1);

        SB_quality = new QSpinBox(Export);
        SB_quality->setObjectName(QStringLiteral("SB_quality"));
        SB_quality->setReadOnly(true);
        SB_quality->setMinimum(-1);
        SB_quality->setMaximum(100);
        SB_quality->setSingleStep(1);
        SB_quality->setValue(-1);

        gridLayout_2->addWidget(SB_quality, 4, 1, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_2->addItem(horizontalSpacer_2, 4, 2, 1, 1);

        label_5 = new QLabel(Export);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_5, 5, 0, 1, 1);

        LE_pdftitle = new QLineEdit(Export);
        LE_pdftitle->setObjectName(QStringLiteral("LE_pdftitle"));

        gridLayout_2->addWidget(LE_pdftitle, 5, 1, 1, 3);


        gridLayout_3->addLayout(gridLayout_2, 0, 1, 1, 1);

        BB_exit = new QDialogButtonBox(Export);
        BB_exit->setObjectName(QStringLiteral("BB_exit"));
        BB_exit->setOrientation(Qt::Horizontal);
        BB_exit->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_3->addWidget(BB_exit, 1, 1, 1, 1);


        retranslateUi(Export);
        QObject::connect(BB_exit, SIGNAL(accepted()), Export, SLOT(accept()));
        QObject::connect(BB_exit, SIGNAL(rejected()), Export, SLOT(reject()));

        QMetaObject::connectSlotsByName(Export);
    } // setupUi

    void retranslateUi(QDialog *Export)
    {
        Export->setWindowTitle(QApplication::translate("Export", "Export Plot", Q_NULLPTR));
        GB_export->setTitle(QApplication::translate("Export", "Export", Q_NULLPTR));
        RB_pdf->setText(QApplication::translate("Export", "PDF", Q_NULLPTR));
        RB_png->setText(QApplication::translate("Export", "PNG", Q_NULLPTR));
        RB_bmp->setText(QApplication::translate("Export", "BMP", Q_NULLPTR));
        RB_jpg->setText(QApplication::translate("Export", "JPG", Q_NULLPTR));
        label->setText(QApplication::translate("Export", "Path + Filename:", Q_NULLPTR));
        PB_file->setText(QApplication::translate("Export", "File", Q_NULLPTR));
        label_2->setText(QApplication::translate("Export", "Cosmetic Pen:", Q_NULLPTR));
        CB_cospen->setText(QString());
        label_3->setText(QApplication::translate("Export", "Width / Height:", Q_NULLPTR));
        label_6->setText(QApplication::translate("Export", "Scale:", Q_NULLPTR));
        label_7->setText(QApplication::translate("Export", "Quality:", Q_NULLPTR));
        label_5->setText(QApplication::translate("Export", "PDF Title:", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Export: public Ui_Export {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EXPORT_H
