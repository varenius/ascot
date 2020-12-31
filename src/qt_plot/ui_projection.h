/********************************************************************************
** Form generated from reading UI file 'projection.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PROJECTION_H
#define UI_PROJECTION_H

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

class Ui_Projection
{
public:
    QGridLayout *gridLayout_5;
    QDialogButtonBox *buttonBox;
    QTabWidget *tabWidget;
    QWidget *symbolTab;
    QGridLayout *gridLayout_3;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_2;
    QRadioButton *mollweide;
    QRadioButton *naturalearth;
    QRadioButton *robinson;
    QRadioButton *mercator;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_4;
    QRadioButton *qrb_line;
    QRadioButton *qrb_disc;
    QRadioButton *qrb_flat;
    QRadioButton *qrb_spike;
    QRadioButton *qrb_none;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_6;
    QRadioButton *stepright;
    QRadioButton *impulse;
    QRadioButton *line;
    QRadioButton *stepleft;
    QRadioButton *none;
    QGridLayout *gridLayout;
    QLabel *label;
    QDoubleSpinBox *gridwidth;
    QLabel *label_3;
    QDoubleSpinBox *gridextension;
    QLabel *label_2;
    QDoubleSpinBox *lamspace;
    QLabel *label_4;
    QDoubleSpinBox *phispace;
    QCheckBox *textlabels;

    void setupUi(QDialog *Projection)
    {
        if (Projection->objectName().isEmpty())
            Projection->setObjectName(QStringLiteral("Projection"));
        Projection->resize(645, 415);
        gridLayout_5 = new QGridLayout(Projection);
        gridLayout_5->setObjectName(QStringLiteral("gridLayout_5"));
        buttonBox = new QDialogButtonBox(Projection);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        gridLayout_5->addWidget(buttonBox, 2, 0, 1, 1);

        tabWidget = new QTabWidget(Projection);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        symbolTab = new QWidget();
        symbolTab->setObjectName(QStringLiteral("symbolTab"));
        gridLayout_3 = new QGridLayout(symbolTab);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        groupBox = new QGroupBox(symbolTab);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout_2 = new QGridLayout(groupBox);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        mollweide = new QRadioButton(groupBox);
        mollweide->setObjectName(QStringLiteral("mollweide"));

        gridLayout_2->addWidget(mollweide, 0, 0, 1, 1);

        naturalearth = new QRadioButton(groupBox);
        naturalearth->setObjectName(QStringLiteral("naturalearth"));

        gridLayout_2->addWidget(naturalearth, 1, 0, 1, 1);

        robinson = new QRadioButton(groupBox);
        robinson->setObjectName(QStringLiteral("robinson"));
        robinson->setCheckable(false);

        gridLayout_2->addWidget(robinson, 2, 0, 1, 1);

        mercator = new QRadioButton(groupBox);
        mercator->setObjectName(QStringLiteral("mercator"));
        mercator->setCheckable(false);

        gridLayout_2->addWidget(mercator, 3, 0, 1, 1);


        gridLayout_3->addWidget(groupBox, 0, 0, 1, 1);

        groupBox_2 = new QGroupBox(symbolTab);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        gridLayout_4 = new QGridLayout(groupBox_2);
        gridLayout_4->setObjectName(QStringLiteral("gridLayout_4"));
        qrb_line = new QRadioButton(groupBox_2);
        qrb_line->setObjectName(QStringLiteral("qrb_line"));
        qrb_line->setCheckable(true);

        gridLayout_4->addWidget(qrb_line, 3, 0, 1, 1);

        qrb_disc = new QRadioButton(groupBox_2);
        qrb_disc->setObjectName(QStringLiteral("qrb_disc"));
        qrb_disc->setCheckable(true);

        gridLayout_4->addWidget(qrb_disc, 4, 0, 1, 1);

        qrb_flat = new QRadioButton(groupBox_2);
        qrb_flat->setObjectName(QStringLiteral("qrb_flat"));
        qrb_flat->setCheckable(true);

        gridLayout_4->addWidget(qrb_flat, 1, 0, 1, 1);

        qrb_spike = new QRadioButton(groupBox_2);
        qrb_spike->setObjectName(QStringLiteral("qrb_spike"));
        qrb_spike->setCheckable(true);

        gridLayout_4->addWidget(qrb_spike, 2, 0, 1, 1);

        qrb_none = new QRadioButton(groupBox_2);
        qrb_none->setObjectName(QStringLiteral("qrb_none"));
        qrb_none->setCheckable(true);

        gridLayout_4->addWidget(qrb_none, 0, 0, 1, 1);


        gridLayout_3->addWidget(groupBox_2, 0, 1, 1, 1);

        groupBox_3 = new QGroupBox(symbolTab);
        groupBox_3->setObjectName(QStringLiteral("groupBox_3"));
        gridLayout_6 = new QGridLayout(groupBox_3);
        gridLayout_6->setObjectName(QStringLiteral("gridLayout_6"));
        stepright = new QRadioButton(groupBox_3);
        stepright->setObjectName(QStringLiteral("stepright"));
        stepright->setCheckable(false);

        gridLayout_6->addWidget(stepright, 3, 0, 1, 1);

        impulse = new QRadioButton(groupBox_3);
        impulse->setObjectName(QStringLiteral("impulse"));
        impulse->setCheckable(false);

        gridLayout_6->addWidget(impulse, 4, 0, 1, 1);

        line = new QRadioButton(groupBox_3);
        line->setObjectName(QStringLiteral("line"));
        line->setCheckable(false);

        gridLayout_6->addWidget(line, 1, 0, 1, 1);

        stepleft = new QRadioButton(groupBox_3);
        stepleft->setObjectName(QStringLiteral("stepleft"));
        stepleft->setCheckable(false);

        gridLayout_6->addWidget(stepleft, 2, 0, 1, 1);

        none = new QRadioButton(groupBox_3);
        none->setObjectName(QStringLiteral("none"));
        none->setCheckable(false);

        gridLayout_6->addWidget(none, 0, 0, 1, 1);


        gridLayout_3->addWidget(groupBox_3, 0, 2, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        label = new QLabel(symbolTab);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        gridwidth = new QDoubleSpinBox(symbolTab);
        gridwidth->setObjectName(QStringLiteral("gridwidth"));

        gridLayout->addWidget(gridwidth, 0, 1, 1, 1);

        label_3 = new QLabel(symbolTab);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout->addWidget(label_3, 0, 2, 1, 1);

        gridextension = new QDoubleSpinBox(symbolTab);
        gridextension->setObjectName(QStringLiteral("gridextension"));

        gridLayout->addWidget(gridextension, 0, 3, 1, 1);

        label_2 = new QLabel(symbolTab);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 0, 4, 1, 1);

        lamspace = new QDoubleSpinBox(symbolTab);
        lamspace->setObjectName(QStringLiteral("lamspace"));

        gridLayout->addWidget(lamspace, 0, 5, 1, 1);

        label_4 = new QLabel(symbolTab);
        label_4->setObjectName(QStringLiteral("label_4"));

        gridLayout->addWidget(label_4, 1, 4, 1, 1);

        phispace = new QDoubleSpinBox(symbolTab);
        phispace->setObjectName(QStringLiteral("phispace"));

        gridLayout->addWidget(phispace, 1, 5, 1, 1);

        textlabels = new QCheckBox(symbolTab);
        textlabels->setObjectName(QStringLiteral("textlabels"));

        gridLayout->addWidget(textlabels, 2, 5, 1, 1);


        gridLayout_3->addLayout(gridLayout, 1, 0, 1, 3);

        tabWidget->addTab(symbolTab, QString());

        gridLayout_5->addWidget(tabWidget, 0, 0, 1, 1);


        retranslateUi(Projection);
        QObject::connect(buttonBox, SIGNAL(accepted()), Projection, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Projection, SLOT(reject()));

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(Projection);
    } // setupUi

    void retranslateUi(QDialog *Projection)
    {
        Projection->setWindowTitle(QApplication::translate("Projection", "Dialog", Q_NULLPTR));
        groupBox->setTitle(QApplication::translate("Projection", "Projection", Q_NULLPTR));
        mollweide->setText(QApplication::translate("Projection", "Mollweide", Q_NULLPTR));
        naturalearth->setText(QApplication::translate("Projection", "Natural Earth", Q_NULLPTR));
        robinson->setText(QApplication::translate("Projection", "Robinson", Q_NULLPTR));
        mercator->setText(QApplication::translate("Projection", "Mercator", Q_NULLPTR));
        groupBox_2->setTitle(QApplication::translate("Projection", "Arrow Ending", Q_NULLPTR));
        qrb_line->setText(QApplication::translate("Projection", "line", Q_NULLPTR));
        qrb_disc->setText(QApplication::translate("Projection", "disc", Q_NULLPTR));
        qrb_flat->setText(QApplication::translate("Projection", "flat", Q_NULLPTR));
        qrb_spike->setText(QApplication::translate("Projection", "spike", Q_NULLPTR));
        qrb_none->setText(QApplication::translate("Projection", "none", Q_NULLPTR));
        groupBox_3->setTitle(QApplication::translate("Projection", "Line Style", Q_NULLPTR));
        stepright->setText(QApplication::translate("Projection", "StepRight", Q_NULLPTR));
        impulse->setText(QApplication::translate("Projection", "Impulse", Q_NULLPTR));
        line->setText(QApplication::translate("Projection", "Line", Q_NULLPTR));
        stepleft->setText(QApplication::translate("Projection", "StepLeft", Q_NULLPTR));
        none->setText(QApplication::translate("Projection", "None", Q_NULLPTR));
        label->setText(QApplication::translate("Projection", "Grid Width:", Q_NULLPTR));
        gridwidth->setPrefix(QString());
        label_3->setText(QApplication::translate("Projection", "Grid Extension:", Q_NULLPTR));
        gridextension->setPrefix(QString());
        label_2->setText(QApplication::translate("Projection", "Lambda Space:", Q_NULLPTR));
        lamspace->setPrefix(QString());
        label_4->setText(QApplication::translate("Projection", "Phi Space:", Q_NULLPTR));
        textlabels->setText(QApplication::translate("Projection", "Textlabels", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(symbolTab), QApplication::translate("Projection", "Style", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Projection: public Ui_Projection {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PROJECTION_H
