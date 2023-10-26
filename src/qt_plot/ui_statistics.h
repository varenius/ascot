/********************************************************************************
** Form generated from reading UI file 'statistics.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_STATISTICS_H
#define UI_STATISTICS_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableView>
#include <QtWidgets/QWidget>
#include "plot.h"

QT_BEGIN_NAMESPACE

class Ui_Statistics
{
public:
    QWidget *central_widget;
    QGridLayout *gridLayout_4;
    QGridLayout *gridLayout_3;
    Plot *sky_widget;
    Plot *crf_widget;
    Plot *trf_widget;
    QMenuBar *menubar;
    QStatusBar *statusbar;
    QDockWidget *stats_dockwidget;
    QWidget *stats_contentwidget;
    QGridLayout *gridLayout_2;
    QTableView *data_tableview;
    QGridLayout *gridLayout;
    QCheckBox *errorbar_checkbox;
    QCheckBox *histogram_checkbox;
    QCheckBox *crosscorr_checkbox;
    QLabel *yaxis_label;
    QSpinBox *histogram_spinbox;
    QCheckBox *absolute_checkbox;
    QLabel *xaxis_label;
    QCheckBox *outlier_checkbox;
    QCheckBox *boxplot_checkbox;
    QComboBox *yaxis_combobox;
    QComboBox *xaxis_combobox;
    QDialogButtonBox *plot_buttonbox;
    QPushButton *cbreak_pushbutton;
    QHBoxLayout *horizontalLayout;
    QPushButton *residUp_pushbotton;
    QPushButton *residDown_pushbotton;
    QPushButton *save_pushbutton;
    QDockWidget *plot_dockwidget;
    QWidget *plot_contentwidget;
    QGridLayout *gridLayout_5;
    QTabWidget *plot_tabwidget;
    QWidget *tab;
    QWidget *tab_2;
    QDockWidget *hist_dockwidget;
    QWidget *dockWidgetContents;

    void setupUi(QMainWindow *Statistics)
    {
        if (Statistics->objectName().isEmpty())
            Statistics->setObjectName(QString::fromUtf8("Statistics"));
        Statistics->resize(1235, 784);
        Statistics->setMinimumSize(QSize(0, 0));
        central_widget = new QWidget(Statistics);
        central_widget->setObjectName(QString::fromUtf8("central_widget"));
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(central_widget->sizePolicy().hasHeightForWidth());
        central_widget->setSizePolicy(sizePolicy);
        gridLayout_4 = new QGridLayout(central_widget);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        sky_widget = new Plot(central_widget);
        sky_widget->setObjectName(QString::fromUtf8("sky_widget"));
        sizePolicy.setHeightForWidth(sky_widget->sizePolicy().hasHeightForWidth());
        sky_widget->setSizePolicy(sizePolicy);
        sky_widget->setSizeIncrement(QSize(1, 1));

        gridLayout_3->addWidget(sky_widget, 0, 0, 1, 1);

        crf_widget = new Plot(central_widget);
        crf_widget->setObjectName(QString::fromUtf8("crf_widget"));
        sizePolicy.setHeightForWidth(crf_widget->sizePolicy().hasHeightForWidth());
        crf_widget->setSizePolicy(sizePolicy);

        gridLayout_3->addWidget(crf_widget, 0, 1, 1, 1);

        trf_widget = new Plot(central_widget);
        trf_widget->setObjectName(QString::fromUtf8("trf_widget"));
        sizePolicy.setHeightForWidth(trf_widget->sizePolicy().hasHeightForWidth());
        trf_widget->setSizePolicy(sizePolicy);

        gridLayout_3->addWidget(trf_widget, 0, 2, 1, 1);


        gridLayout_4->addLayout(gridLayout_3, 0, 0, 1, 1);

        Statistics->setCentralWidget(central_widget);
        menubar = new QMenuBar(Statistics);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1235, 24));
        Statistics->setMenuBar(menubar);
        statusbar = new QStatusBar(Statistics);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        Statistics->setStatusBar(statusbar);
        stats_dockwidget = new QDockWidget(Statistics);
        stats_dockwidget->setObjectName(QString::fromUtf8("stats_dockwidget"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(20);
        sizePolicy1.setVerticalStretch(20);
        sizePolicy1.setHeightForWidth(stats_dockwidget->sizePolicy().hasHeightForWidth());
        stats_dockwidget->setSizePolicy(sizePolicy1);
        stats_dockwidget->setMinimumSize(QSize(472, 246));
        stats_dockwidget->setMaximumSize(QSize(472, 734));
        stats_contentwidget = new QWidget();
        stats_contentwidget->setObjectName(QString::fromUtf8("stats_contentwidget"));
        gridLayout_2 = new QGridLayout(stats_contentwidget);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        data_tableview = new QTableView(stats_contentwidget);
        data_tableview->setObjectName(QString::fromUtf8("data_tableview"));

        gridLayout_2->addWidget(data_tableview, 0, 0, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        errorbar_checkbox = new QCheckBox(stats_contentwidget);
        errorbar_checkbox->setObjectName(QString::fromUtf8("errorbar_checkbox"));
        errorbar_checkbox->setChecked(true);

        gridLayout->addWidget(errorbar_checkbox, 2, 1, 1, 1);

        histogram_checkbox = new QCheckBox(stats_contentwidget);
        histogram_checkbox->setObjectName(QString::fromUtf8("histogram_checkbox"));

        gridLayout->addWidget(histogram_checkbox, 0, 2, 1, 1);

        crosscorr_checkbox = new QCheckBox(stats_contentwidget);
        crosscorr_checkbox->setObjectName(QString::fromUtf8("crosscorr_checkbox"));

        gridLayout->addWidget(crosscorr_checkbox, 2, 0, 1, 1);

        yaxis_label = new QLabel(stats_contentwidget);
        yaxis_label->setObjectName(QString::fromUtf8("yaxis_label"));

        gridLayout->addWidget(yaxis_label, 1, 0, 1, 1);

        histogram_spinbox = new QSpinBox(stats_contentwidget);
        histogram_spinbox->setObjectName(QString::fromUtf8("histogram_spinbox"));
        histogram_spinbox->setMinimumSize(QSize(71, 25));
        histogram_spinbox->setMaximumSize(QSize(71, 25));
        histogram_spinbox->setMinimum(0);
        histogram_spinbox->setMaximum(1000);
        histogram_spinbox->setValue(0);

        gridLayout->addWidget(histogram_spinbox, 0, 3, 1, 1);

        absolute_checkbox = new QCheckBox(stats_contentwidget);
        absolute_checkbox->setObjectName(QString::fromUtf8("absolute_checkbox"));

        gridLayout->addWidget(absolute_checkbox, 2, 3, 1, 1);

        xaxis_label = new QLabel(stats_contentwidget);
        xaxis_label->setObjectName(QString::fromUtf8("xaxis_label"));

        gridLayout->addWidget(xaxis_label, 0, 0, 1, 1);

        outlier_checkbox = new QCheckBox(stats_contentwidget);
        outlier_checkbox->setObjectName(QString::fromUtf8("outlier_checkbox"));
        outlier_checkbox->setChecked(true);

        gridLayout->addWidget(outlier_checkbox, 2, 2, 1, 1);

        boxplot_checkbox = new QCheckBox(stats_contentwidget);
        boxplot_checkbox->setObjectName(QString::fromUtf8("boxplot_checkbox"));

        gridLayout->addWidget(boxplot_checkbox, 1, 2, 1, 1);

        yaxis_combobox = new QComboBox(stats_contentwidget);
        yaxis_combobox->setObjectName(QString::fromUtf8("yaxis_combobox"));
        QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(yaxis_combobox->sizePolicy().hasHeightForWidth());
        yaxis_combobox->setSizePolicy(sizePolicy2);
        yaxis_combobox->setMinimumSize(QSize(111, 25));

        gridLayout->addWidget(yaxis_combobox, 1, 1, 1, 1);

        xaxis_combobox = new QComboBox(stats_contentwidget);
        xaxis_combobox->setObjectName(QString::fromUtf8("xaxis_combobox"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(xaxis_combobox->sizePolicy().hasHeightForWidth());
        xaxis_combobox->setSizePolicy(sizePolicy3);
        xaxis_combobox->setMinimumSize(QSize(111, 25));

        gridLayout->addWidget(xaxis_combobox, 0, 1, 1, 1);

        plot_buttonbox = new QDialogButtonBox(stats_contentwidget);
        plot_buttonbox->setObjectName(QString::fromUtf8("plot_buttonbox"));
        sizePolicy2.setHeightForWidth(plot_buttonbox->sizePolicy().hasHeightForWidth());
        plot_buttonbox->setSizePolicy(sizePolicy2);
        plot_buttonbox->setStandardButtons(QDialogButtonBox::Ok);

        gridLayout->addWidget(plot_buttonbox, 2, 5, 1, 1);

        cbreak_pushbutton = new QPushButton(stats_contentwidget);
        cbreak_pushbutton->setObjectName(QString::fromUtf8("cbreak_pushbutton"));
        cbreak_pushbutton->setEnabled(false);
        sizePolicy2.setHeightForWidth(cbreak_pushbutton->sizePolicy().hasHeightForWidth());
        cbreak_pushbutton->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(cbreak_pushbutton, 1, 3, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        residUp_pushbotton = new QPushButton(stats_contentwidget);
        residUp_pushbotton->setObjectName(QString::fromUtf8("residUp_pushbotton"));
        sizePolicy.setHeightForWidth(residUp_pushbotton->sizePolicy().hasHeightForWidth());
        residUp_pushbotton->setSizePolicy(sizePolicy);
        residUp_pushbotton->setMaximumSize(QSize(35, 16777215));

        horizontalLayout->addWidget(residUp_pushbotton);

        residDown_pushbotton = new QPushButton(stats_contentwidget);
        residDown_pushbotton->setObjectName(QString::fromUtf8("residDown_pushbotton"));
        sizePolicy.setHeightForWidth(residDown_pushbotton->sizePolicy().hasHeightForWidth());
        residDown_pushbotton->setSizePolicy(sizePolicy);
        residDown_pushbotton->setMaximumSize(QSize(35, 16777215));

        horizontalLayout->addWidget(residDown_pushbotton);


        gridLayout->addLayout(horizontalLayout, 1, 5, 1, 1);

        save_pushbutton = new QPushButton(stats_contentwidget);
        save_pushbutton->setObjectName(QString::fromUtf8("save_pushbutton"));
        sizePolicy2.setHeightForWidth(save_pushbutton->sizePolicy().hasHeightForWidth());
        save_pushbutton->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(save_pushbutton, 0, 5, 1, 1);


        gridLayout_2->addLayout(gridLayout, 1, 0, 1, 1);

        stats_dockwidget->setWidget(stats_contentwidget);
        Statistics->addDockWidget(Qt::TopDockWidgetArea, stats_dockwidget);
        plot_dockwidget = new QDockWidget(Statistics);
        plot_dockwidget->setObjectName(QString::fromUtf8("plot_dockwidget"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy4.setHorizontalStretch(70);
        sizePolicy4.setVerticalStretch(70);
        sizePolicy4.setHeightForWidth(plot_dockwidget->sizePolicy().hasHeightForWidth());
        plot_dockwidget->setSizePolicy(sizePolicy4);
        plot_dockwidget->setMinimumSize(QSize(157, 94));
        plot_dockwidget->setFocusPolicy(Qt::ClickFocus);
        plot_contentwidget = new QWidget();
        plot_contentwidget->setObjectName(QString::fromUtf8("plot_contentwidget"));
        plot_contentwidget->setMinimumSize(QSize(0, 0));
        gridLayout_5 = new QGridLayout(plot_contentwidget);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        plot_tabwidget = new QTabWidget(plot_contentwidget);
        plot_tabwidget->setObjectName(QString::fromUtf8("plot_tabwidget"));
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        plot_tabwidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        plot_tabwidget->addTab(tab_2, QString());

        gridLayout_5->addWidget(plot_tabwidget, 0, 0, 1, 1);

        plot_dockwidget->setWidget(plot_contentwidget);
        Statistics->addDockWidget(Qt::TopDockWidgetArea, plot_dockwidget);
        hist_dockwidget = new QDockWidget(Statistics);
        hist_dockwidget->setObjectName(QString::fromUtf8("hist_dockwidget"));
        QSizePolicy sizePolicy5(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(hist_dockwidget->sizePolicy().hasHeightForWidth());
        hist_dockwidget->setSizePolicy(sizePolicy5);
        hist_dockwidget->setMinimumSize(QSize(79, 44));
        dockWidgetContents = new QWidget();
        dockWidgetContents->setObjectName(QString::fromUtf8("dockWidgetContents"));
        dockWidgetContents->setMinimumSize(QSize(0, 0));
        hist_dockwidget->setWidget(dockWidgetContents);
        Statistics->addDockWidget(Qt::TopDockWidgetArea, hist_dockwidget);

        retranslateUi(Statistics);

        QMetaObject::connectSlotsByName(Statistics);
    } // setupUi

    void retranslateUi(QMainWindow *Statistics)
    {
        Statistics->setWindowTitle(QCoreApplication::translate("Statistics", "MainWindow", nullptr));
        errorbar_checkbox->setText(QCoreApplication::translate("Statistics", "Errorbar", nullptr));
        histogram_checkbox->setText(QCoreApplication::translate("Statistics", "Histogram", nullptr));
        crosscorr_checkbox->setText(QCoreApplication::translate("Statistics", "CrCo", nullptr));
        yaxis_label->setText(QCoreApplication::translate("Statistics", "Y-Axis:", nullptr));
        absolute_checkbox->setText(QCoreApplication::translate("Statistics", "absolute", nullptr));
        xaxis_label->setText(QCoreApplication::translate("Statistics", "X-Axis:", nullptr));
        outlier_checkbox->setText(QCoreApplication::translate("Statistics", "Outlier", nullptr));
        boxplot_checkbox->setText(QCoreApplication::translate("Statistics", "BoxPlot", nullptr));
        cbreak_pushbutton->setText(QCoreApplication::translate("Statistics", "CBreak", nullptr));
        residUp_pushbotton->setText(QCoreApplication::translate("Statistics", "^", nullptr));
        residDown_pushbotton->setText(QCoreApplication::translate("Statistics", "v", nullptr));
        save_pushbutton->setText(QCoreApplication::translate("Statistics", "save", nullptr));
        plot_tabwidget->setTabText(plot_tabwidget->indexOf(tab), QCoreApplication::translate("Statistics", "Tab 1", nullptr));
        plot_tabwidget->setTabText(plot_tabwidget->indexOf(tab_2), QCoreApplication::translate("Statistics", "Tab 2", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Statistics: public Ui_Statistics {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_STATISTICS_H
