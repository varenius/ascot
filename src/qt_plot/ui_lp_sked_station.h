/********************************************************************************
** Form generated from reading UI file 'lp_sked_station.ui'
**
** Created by: Qt User Interface Compiler version 5.12.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_LP_SKED_STATION_H
#define UI_LP_SKED_STATION_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLCDNumber>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_StationDialog
{
public:
    QVBoxLayout *verticalLayout;
    QGraphicsView *graphicsView;
    QGridLayout *gridLayout;
    QCheckBox *checkBox_6;
    QCheckBox *checkBox_3;
    QCheckBox *checkBox_7;
    QCheckBox *checkBox_5;
    QCheckBox *checkBox;
    QHBoxLayout *horizontalLayout_6;
    QCheckBox *checkBox_4;
    QCheckBox *checkBox_2;
    QGroupBox *groupBox_2;
    QHBoxLayout *horizontalLayout_2;
    QHBoxLayout *horizontalLayout_8;
    QPushButton *pushButton_2;
    QLabel *label;
    QSpinBox *spinBox_2;
    QLabel *label_2;
    QSpinBox *spinBox;
    QGroupBox *groupBox;
    QHBoxLayout *horizontalLayout_3;
    QHBoxLayout *horizontalLayout_7;
    QPushButton *pushButton_3;
    QLabel *label_3;
    QSpinBox *spinBox_3;
    QLabel *label_4;
    QSpinBox *spinBox_4;
    QCheckBox *checkBox_8;
    QGroupBox *groupBox_3;
    QHBoxLayout *horizontalLayout_5;
    QSlider *horizontalSlider;
    QLCDNumber *lcdNumber;
    QLabel *label_5;
    QHBoxLayout *horizontalLayout;
    QLineEdit *lineEdit;
    QPushButton *pushButton;

    void setupUi(QDialog *StationDialog)
    {
        if (StationDialog->objectName().isEmpty())
            StationDialog->setObjectName(QString::fromUtf8("StationDialog"));
        StationDialog->setEnabled(true);
        StationDialog->resize(549, 753);
        verticalLayout = new QVBoxLayout(StationDialog);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        graphicsView = new QGraphicsView(StationDialog);
        graphicsView->setObjectName(QString::fromUtf8("graphicsView"));
        graphicsView->setMouseTracking(true);

        verticalLayout->addWidget(graphicsView);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        checkBox_6 = new QCheckBox(StationDialog);
        checkBox_6->setObjectName(QString::fromUtf8("checkBox_6"));

        gridLayout->addWidget(checkBox_6, 3, 1, 1, 1);

        checkBox_3 = new QCheckBox(StationDialog);
        checkBox_3->setObjectName(QString::fromUtf8("checkBox_3"));
        checkBox_3->setChecked(true);

        gridLayout->addWidget(checkBox_3, 0, 0, 1, 1);

        checkBox_7 = new QCheckBox(StationDialog);
        checkBox_7->setObjectName(QString::fromUtf8("checkBox_7"));
        checkBox_7->setChecked(true);

        gridLayout->addWidget(checkBox_7, 4, 0, 1, 1);

        checkBox_5 = new QCheckBox(StationDialog);
        checkBox_5->setObjectName(QString::fromUtf8("checkBox_5"));

        gridLayout->addWidget(checkBox_5, 0, 1, 1, 1);

        checkBox = new QCheckBox(StationDialog);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));

        gridLayout->addWidget(checkBox, 4, 1, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        checkBox_4 = new QCheckBox(StationDialog);
        checkBox_4->setObjectName(QString::fromUtf8("checkBox_4"));
        checkBox_4->setChecked(true);

        horizontalLayout_6->addWidget(checkBox_4);

        checkBox_2 = new QCheckBox(StationDialog);
        checkBox_2->setObjectName(QString::fromUtf8("checkBox_2"));

        horizontalLayout_6->addWidget(checkBox_2);


        gridLayout->addLayout(horizontalLayout_6, 3, 0, 1, 1);


        verticalLayout->addLayout(gridLayout);

        groupBox_2 = new QGroupBox(StationDialog);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(groupBox_2->sizePolicy().hasHeightForWidth());
        groupBox_2->setSizePolicy(sizePolicy);
        groupBox_2->setAutoFillBackground(false);
        groupBox_2->setCheckable(true);
        horizontalLayout_2 = new QHBoxLayout(groupBox_2);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        pushButton_2 = new QPushButton(groupBox_2);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));

        horizontalLayout_8->addWidget(pushButton_2);

        label = new QLabel(groupBox_2);
        label->setObjectName(QString::fromUtf8("label"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy1);

        horizontalLayout_8->addWidget(label);

        spinBox_2 = new QSpinBox(groupBox_2);
        spinBox_2->setObjectName(QString::fromUtf8("spinBox_2"));
        spinBox_2->setMaximum(20);
        spinBox_2->setValue(4);

        horizontalLayout_8->addWidget(spinBox_2);

        label_2 = new QLabel(groupBox_2);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        sizePolicy1.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy1);

        horizontalLayout_8->addWidget(label_2);

        spinBox = new QSpinBox(groupBox_2);
        spinBox->setObjectName(QString::fromUtf8("spinBox"));
        spinBox->setValue(2);

        horizontalLayout_8->addWidget(spinBox);


        horizontalLayout_2->addLayout(horizontalLayout_8);


        verticalLayout->addWidget(groupBox_2);

        groupBox = new QGroupBox(StationDialog);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        sizePolicy.setHeightForWidth(groupBox->sizePolicy().hasHeightForWidth());
        groupBox->setSizePolicy(sizePolicy);
        groupBox->setCheckable(true);
        groupBox->setChecked(false);
        horizontalLayout_3 = new QHBoxLayout(groupBox);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        pushButton_3 = new QPushButton(groupBox);
        pushButton_3->setObjectName(QString::fromUtf8("pushButton_3"));

        horizontalLayout_7->addWidget(pushButton_3);

        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        sizePolicy1.setHeightForWidth(label_3->sizePolicy().hasHeightForWidth());
        label_3->setSizePolicy(sizePolicy1);

        horizontalLayout_7->addWidget(label_3);

        spinBox_3 = new QSpinBox(groupBox);
        spinBox_3->setObjectName(QString::fromUtf8("spinBox_3"));

        horizontalLayout_7->addWidget(spinBox_3);

        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        sizePolicy1.setHeightForWidth(label_4->sizePolicy().hasHeightForWidth());
        label_4->setSizePolicy(sizePolicy1);

        horizontalLayout_7->addWidget(label_4);

        spinBox_4 = new QSpinBox(groupBox);
        spinBox_4->setObjectName(QString::fromUtf8("spinBox_4"));
        spinBox_4->setValue(5);

        horizontalLayout_7->addWidget(spinBox_4);

        checkBox_8 = new QCheckBox(groupBox);
        checkBox_8->setObjectName(QString::fromUtf8("checkBox_8"));
        checkBox_8->setChecked(true);

        horizontalLayout_7->addWidget(checkBox_8);


        horizontalLayout_3->addLayout(horizontalLayout_7);


        verticalLayout->addWidget(groupBox);

        groupBox_3 = new QGroupBox(StationDialog);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        groupBox_3->setEnabled(true);
        groupBox_3->setCheckable(true);
        groupBox_3->setChecked(false);
        horizontalLayout_5 = new QHBoxLayout(groupBox_3);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        horizontalSlider = new QSlider(groupBox_3);
        horizontalSlider->setObjectName(QString::fromUtf8("horizontalSlider"));
        horizontalSlider->setEnabled(true);
        horizontalSlider->setPageStep(1);
        horizontalSlider->setOrientation(Qt::Horizontal);
        horizontalSlider->setTickPosition(QSlider::TicksBelow);

        horizontalLayout_5->addWidget(horizontalSlider);

        lcdNumber = new QLCDNumber(groupBox_3);
        lcdNumber->setObjectName(QString::fromUtf8("lcdNumber"));
        lcdNumber->setEnabled(true);
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        lcdNumber->setFont(font);
        lcdNumber->setAutoFillBackground(false);
        lcdNumber->setFrameShadow(QFrame::Raised);
        lcdNumber->setSmallDecimalPoint(false);
        lcdNumber->setDigitCount(3);
        lcdNumber->setSegmentStyle(QLCDNumber::Flat);

        horizontalLayout_5->addWidget(lcdNumber);

        label_5 = new QLabel(groupBox_3);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setEnabled(true);
        label_5->setMinimumSize(QSize(80, 0));
        label_5->setFrameShape(QFrame::Box);
        label_5->setAlignment(Qt::AlignCenter);

        horizontalLayout_5->addWidget(label_5);


        verticalLayout->addWidget(groupBox_3);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        lineEdit = new QLineEdit(StationDialog);
        lineEdit->setObjectName(QString::fromUtf8("lineEdit"));

        horizontalLayout->addWidget(lineEdit);

        pushButton = new QPushButton(StationDialog);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));

        horizontalLayout->addWidget(pushButton);


        verticalLayout->addLayout(horizontalLayout);


        retranslateUi(StationDialog);

        QMetaObject::connectSlotsByName(StationDialog);
    } // setupUi

    void retranslateUi(QDialog *StationDialog)
    {
        StationDialog->setWindowTitle(QApplication::translate("StationDialog", "Dialog", nullptr));
        checkBox_6->setText(QApplication::translate("StationDialog", "Show Sun Transit", nullptr));
        checkBox_3->setText(QApplication::translate("StationDialog", "Show Transits", nullptr));
        checkBox_7->setText(QApplication::translate("StationDialog", "Show Pole", nullptr));
        checkBox_5->setText(QApplication::translate("StationDialog", "Show Elevation Mask", nullptr));
        checkBox->setText(QApplication::translate("StationDialog", "observation order", nullptr));
        checkBox_4->setText(QApplication::translate("StationDialog", "Show Observations", nullptr));
        checkBox_2->setText(QApplication::translate("StationDialog", "Twin", nullptr));
        groupBox_2->setTitle(QApplication::translate("StationDialog", "Grid", nullptr));
        pushButton_2->setText(QApplication::translate("StationDialog", "Refresh Grid", nullptr));
        label->setText(QApplication::translate("StationDialog", "az:", nullptr));
        label_2->setText(QApplication::translate("StationDialog", "el:", nullptr));
        groupBox->setTitle(QApplication::translate("StationDialog", "Tree", nullptr));
        pushButton_3->setText(QApplication::translate("StationDialog", "Refresh Tree", nullptr));
        label_3->setText(QApplication::translate("StationDialog", "min:", nullptr));
        label_4->setText(QApplication::translate("StationDialog", "max:", nullptr));
        checkBox_8->setText(QApplication::translate("StationDialog", "obs", nullptr));
        groupBox_3->setTitle(QApplication::translate("StationDialog", "Interval", nullptr));
        label_5->setText(QApplication::translate("StationDialog", "-", nullptr));
        pushButton->setText(QApplication::translate("StationDialog", "save pdf", nullptr));
    } // retranslateUi

};

namespace Ui {
    class StationDialog: public Ui_StationDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_LP_SKED_STATION_H
