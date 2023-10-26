/********************************************************************************
** Form generated from reading UI file 'analyzer.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ANALYZER_H
#define UI_ANALYZER_H

#include <QtCore/QDate>
#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDateEdit>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QTreeView>
#include <QtWidgets/QWidget>
#include "plot.h"

QT_BEGIN_NAMESPACE

class Ui_Analyzer
{
public:
    QWidget *centralwidget;
    QGridLayout *gridLayout_12;
    QPushButton *directory_pushbutton;
    QPushButton *clear_pushbutton;
    QPushButton *refresh_pushbutton;
    QCheckBox *extern_checkbox;
    QWidget *progressbar_widget;
    QTreeView *folder_treeview;
    QTextBrowser *info_textbrowser;
    QTableWidget *session_tablewidget;
    QTableWidget *param_tablewidget;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_3;
    QLabel *label_2;
    QDateEdit *start_dateedit;
    QDateEdit *end_dateedit;
    QPushButton *clear_param_pushbutton;
    QPushButton *load_pushbutton;
    QTabWidget *tabWidget;
    QWidget *tab;
    QGridLayout *gridLayout_9;
    Plot *widget;
    QWidget *tab_2;
    QTabWidget *configuration_tabwidget;
    QWidget *highligh_tab;
    QGridLayout *gridLayout_6;
    QGroupBox *highlight_groupbox;
    QGridLayout *gridLayout_5;
    QListWidget *highlight_listwidget;
    QListWidget *highlight_site_listwidget;
    QPushButton *highlight_pushbutton;
    QPushButton *highlight_clear_pushbutton;
    QCheckBox *highlight_remove_checkbox;
    QRadioButton *or_site_radiobutton;
    QRadioButton *and_site_radiobutton;
    QWidget *selection_tab;
    QGridLayout *gridLayout_7;
    QTableWidget *select_tablewidget;
    QSpacerItem *horizontalSpacer_3;
    QTabWidget *analysis_tabwidget;
    QWidget *rf_tab;
    QGridLayout *gridLayout_13;
    QGroupBox *analysis_groupbox;
    QGridLayout *gridLayout_2;
    QRadioButton *trf_radiobutton;
    QRadioButton *crf_radiobutton;
    QGroupBox *refframe_groupbox;
    QGridLayout *gridLayout;
    QRadioButton *natural_radiobutton;
    QRadioButton *mollweide_radiobutton;
    QRadioButton *mercator_radiobutton;
    QDoubleSpinBox *scale_arrows_spinbox;
    QLabel *label_5;
    QDoubleSpinBox *ref_arrow_spinbox;
    QLabel *ref_arrow_label;
    QPushButton *plot_ref_pushbutton;
    QCheckBox *residuals_checkbox;
    QGroupBox *transformation_groupbox;
    QGridLayout *gridLayout_4;
    QCheckBox *tx_checkbox;
    QCheckBox *rx_checkbox;
    QCheckBox *s_checkbox;
    QCheckBox *ty_checkbox;
    QCheckBox *ry_checkbox;
    QPushButton *transform_pushbutton;
    QCheckBox *tz_checkbox;
    QCheckBox *rz_checkbox;
    QCheckBox *def_only_checkbox;
    QCheckBox *except_sh_checkbox;
    QDoubleSpinBox *arrow_width_spinbox;
    QWidget *sta_tab;
    QGridLayout *gridLayout_8;
    QCheckBox *sta_timeseries_checkbox;
    QSpacerItem *horizontalSpacer;
    QDoubleSpinBox *maxest_spinbox;
    QLabel *label_6;
    QSpacerItem *horizontalSpacer_2;
    QCheckBox *sta_baselines_checkbox;
    QSpinBox *blreps_spinbox;
    QLabel *label_7;
    QCheckBox *sta_helmert_checkbox;
    QGridLayout *gridLayout_11;
    QCheckBox *tx_sta_checkbox;
    QCheckBox *ty_sta_checkbox;
    QCheckBox *tz_sta_checkbox;
    QCheckBox *rx_sta_checkbox;
    QCheckBox *ry_sta_checkbox;
    QCheckBox *rz_sta_checkbox;
    QCheckBox *s_sta_checkbox;
    QPushButton *plot_sta_pushbutton;
    QWidget *eop_tab;
    QListWidget *eops_listwidget;
    QPushButton *c04_pushbutton;
    QGroupBox *groupBox_5;
    QGridLayout *gridLayout_10;
    QComboBox *yaxis_combobox;
    QComboBox *xaxis_combobox;
    QLineEdit *timeformat_lineedit;
    QComboBox *interpolateComboBox;
    QLabel *label;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_8;
    QPushButton *plot_eop_pushbutton;
    QWidget *src_tab;
    QMenuBar *menubar;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *Analyzer)
    {
        if (Analyzer->objectName().isEmpty())
            Analyzer->setObjectName(QString::fromUtf8("Analyzer"));
        Analyzer->resize(1115, 851);
        centralwidget = new QWidget(Analyzer);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout_12 = new QGridLayout(centralwidget);
        gridLayout_12->setObjectName(QString::fromUtf8("gridLayout_12"));
        directory_pushbutton = new QPushButton(centralwidget);
        directory_pushbutton->setObjectName(QString::fromUtf8("directory_pushbutton"));
        directory_pushbutton->setMinimumSize(QSize(0, 0));
        directory_pushbutton->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_12->addWidget(directory_pushbutton, 0, 0, 1, 1);

        clear_pushbutton = new QPushButton(centralwidget);
        clear_pushbutton->setObjectName(QString::fromUtf8("clear_pushbutton"));
        clear_pushbutton->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_12->addWidget(clear_pushbutton, 0, 1, 1, 1);

        refresh_pushbutton = new QPushButton(centralwidget);
        refresh_pushbutton->setObjectName(QString::fromUtf8("refresh_pushbutton"));

        gridLayout_12->addWidget(refresh_pushbutton, 0, 2, 1, 1);

        extern_checkbox = new QCheckBox(centralwidget);
        extern_checkbox->setObjectName(QString::fromUtf8("extern_checkbox"));

        gridLayout_12->addWidget(extern_checkbox, 0, 3, 1, 1);

        progressbar_widget = new QWidget(centralwidget);
        progressbar_widget->setObjectName(QString::fromUtf8("progressbar_widget"));

        gridLayout_12->addWidget(progressbar_widget, 0, 4, 1, 2);

        folder_treeview = new QTreeView(centralwidget);
        folder_treeview->setObjectName(QString::fromUtf8("folder_treeview"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(folder_treeview->sizePolicy().hasHeightForWidth());
        folder_treeview->setSizePolicy(sizePolicy);
        folder_treeview->setMinimumSize(QSize(0, 0));
        folder_treeview->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_12->addWidget(folder_treeview, 1, 0, 1, 5);

        info_textbrowser = new QTextBrowser(centralwidget);
        info_textbrowser->setObjectName(QString::fromUtf8("info_textbrowser"));
        info_textbrowser->setMinimumSize(QSize(0, 0));
        info_textbrowser->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_12->addWidget(info_textbrowser, 1, 5, 1, 2);

        session_tablewidget = new QTableWidget(centralwidget);
        if (session_tablewidget->columnCount() < 7)
            session_tablewidget->setColumnCount(7);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        session_tablewidget->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        session_tablewidget->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        session_tablewidget->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        session_tablewidget->setHorizontalHeaderItem(3, __qtablewidgetitem3);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        session_tablewidget->setHorizontalHeaderItem(4, __qtablewidgetitem4);
        QTableWidgetItem *__qtablewidgetitem5 = new QTableWidgetItem();
        session_tablewidget->setHorizontalHeaderItem(5, __qtablewidgetitem5);
        QTableWidgetItem *__qtablewidgetitem6 = new QTableWidgetItem();
        session_tablewidget->setHorizontalHeaderItem(6, __qtablewidgetitem6);
        session_tablewidget->setObjectName(QString::fromUtf8("session_tablewidget"));
        sizePolicy.setHeightForWidth(session_tablewidget->sizePolicy().hasHeightForWidth());
        session_tablewidget->setSizePolicy(sizePolicy);
        session_tablewidget->setMinimumSize(QSize(0, 0));
        session_tablewidget->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_12->addWidget(session_tablewidget, 2, 0, 1, 5);

        param_tablewidget = new QTableWidget(centralwidget);
        param_tablewidget->setObjectName(QString::fromUtf8("param_tablewidget"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(param_tablewidget->sizePolicy().hasHeightForWidth());
        param_tablewidget->setSizePolicy(sizePolicy1);
        param_tablewidget->setMinimumSize(QSize(0, 0));
        param_tablewidget->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_12->addWidget(param_tablewidget, 2, 5, 1, 2);

        groupBox = new QGroupBox(centralwidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setMinimumSize(QSize(0, 0));
        groupBox->setMaximumSize(QSize(566, 16777215));
        gridLayout_3 = new QGridLayout(groupBox);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setMaximumSize(QSize(42, 16777215));
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_2, 0, 0, 1, 1);

        start_dateedit = new QDateEdit(groupBox);
        start_dateedit->setObjectName(QString::fromUtf8("start_dateedit"));
        start_dateedit->setMinimumSize(QSize(83, 0));
        start_dateedit->setMaximumSize(QSize(95, 16777215));
        start_dateedit->setDateTime(QDateTime(QDate(2079, 1, 1), QTime(0, 0, 0)));
        start_dateedit->setMaximumDateTime(QDateTime(QDate(7950, 12, 31), QTime(23, 59, 59)));
        start_dateedit->setMinimumDateTime(QDateTime(QDate(1779, 1, 1), QTime(0, 0, 0)));

        gridLayout_3->addWidget(start_dateedit, 0, 1, 1, 1);

        end_dateedit = new QDateEdit(groupBox);
        end_dateedit->setObjectName(QString::fromUtf8("end_dateedit"));
        end_dateedit->setMinimumSize(QSize(82, 0));
        end_dateedit->setMaximumSize(QSize(95, 16777215));
        end_dateedit->setDateTime(QDateTime(QDate(2020, 1, 1), QTime(0, 0, 0)));
        end_dateedit->setMaximumDate(QDate(7950, 12, 31));
        end_dateedit->setMinimumDate(QDate(1779, 1, 1));
        end_dateedit->setDisplayFormat(QString::fromUtf8("dd MMM yy"));

        gridLayout_3->addWidget(end_dateedit, 0, 2, 1, 1);

        clear_param_pushbutton = new QPushButton(groupBox);
        clear_param_pushbutton->setObjectName(QString::fromUtf8("clear_param_pushbutton"));
        clear_param_pushbutton->setMaximumSize(QSize(111, 16777215));

        gridLayout_3->addWidget(clear_param_pushbutton, 0, 4, 1, 1);

        load_pushbutton = new QPushButton(groupBox);
        load_pushbutton->setObjectName(QString::fromUtf8("load_pushbutton"));
        load_pushbutton->setMaximumSize(QSize(85, 16777215));

        gridLayout_3->addWidget(load_pushbutton, 0, 3, 1, 1);


        gridLayout_12->addWidget(groupBox, 3, 0, 1, 5);

        tabWidget = new QTabWidget(centralwidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setMinimumSize(QSize(0, 0));
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        gridLayout_9 = new QGridLayout(tab);
        gridLayout_9->setObjectName(QString::fromUtf8("gridLayout_9"));
        widget = new Plot(tab);
        widget->setObjectName(QString::fromUtf8("widget"));
        sizePolicy1.setHeightForWidth(widget->sizePolicy().hasHeightForWidth());
        widget->setSizePolicy(sizePolicy1);

        gridLayout_9->addWidget(widget, 0, 0, 1, 1);

        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        tabWidget->addTab(tab_2, QString());

        gridLayout_12->addWidget(tabWidget, 3, 5, 3, 2);

        configuration_tabwidget = new QTabWidget(centralwidget);
        configuration_tabwidget->setObjectName(QString::fromUtf8("configuration_tabwidget"));
        highligh_tab = new QWidget();
        highligh_tab->setObjectName(QString::fromUtf8("highligh_tab"));
        gridLayout_6 = new QGridLayout(highligh_tab);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        highlight_groupbox = new QGroupBox(highligh_tab);
        highlight_groupbox->setObjectName(QString::fromUtf8("highlight_groupbox"));
        sizePolicy.setHeightForWidth(highlight_groupbox->sizePolicy().hasHeightForWidth());
        highlight_groupbox->setSizePolicy(sizePolicy);
        highlight_groupbox->setMaximumSize(QSize(448, 16777215));
        gridLayout_5 = new QGridLayout(highlight_groupbox);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        highlight_listwidget = new QListWidget(highlight_groupbox);
        highlight_listwidget->setObjectName(QString::fromUtf8("highlight_listwidget"));
        sizePolicy.setHeightForWidth(highlight_listwidget->sizePolicy().hasHeightForWidth());
        highlight_listwidget->setSizePolicy(sizePolicy);

        gridLayout_5->addWidget(highlight_listwidget, 0, 0, 4, 1);

        highlight_site_listwidget = new QListWidget(highlight_groupbox);
        highlight_site_listwidget->setObjectName(QString::fromUtf8("highlight_site_listwidget"));
        sizePolicy.setHeightForWidth(highlight_site_listwidget->sizePolicy().hasHeightForWidth());
        highlight_site_listwidget->setSizePolicy(sizePolicy);

        gridLayout_5->addWidget(highlight_site_listwidget, 0, 1, 4, 1);

        highlight_pushbutton = new QPushButton(highlight_groupbox);
        highlight_pushbutton->setObjectName(QString::fromUtf8("highlight_pushbutton"));
        highlight_pushbutton->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_5->addWidget(highlight_pushbutton, 0, 2, 1, 2);

        highlight_clear_pushbutton = new QPushButton(highlight_groupbox);
        highlight_clear_pushbutton->setObjectName(QString::fromUtf8("highlight_clear_pushbutton"));

        gridLayout_5->addWidget(highlight_clear_pushbutton, 1, 2, 1, 2);

        highlight_remove_checkbox = new QCheckBox(highlight_groupbox);
        highlight_remove_checkbox->setObjectName(QString::fromUtf8("highlight_remove_checkbox"));

        gridLayout_5->addWidget(highlight_remove_checkbox, 2, 2, 1, 2);

        or_site_radiobutton = new QRadioButton(highlight_groupbox);
        or_site_radiobutton->setObjectName(QString::fromUtf8("or_site_radiobutton"));

        gridLayout_5->addWidget(or_site_radiobutton, 3, 2, 1, 1);

        and_site_radiobutton = new QRadioButton(highlight_groupbox);
        and_site_radiobutton->setObjectName(QString::fromUtf8("and_site_radiobutton"));
        and_site_radiobutton->setChecked(true);

        gridLayout_5->addWidget(and_site_radiobutton, 3, 3, 1, 1);


        gridLayout_6->addWidget(highlight_groupbox, 0, 0, 1, 1);

        configuration_tabwidget->addTab(highligh_tab, QString());
        selection_tab = new QWidget();
        selection_tab->setObjectName(QString::fromUtf8("selection_tab"));
        gridLayout_7 = new QGridLayout(selection_tab);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        select_tablewidget = new QTableWidget(selection_tab);
        select_tablewidget->setObjectName(QString::fromUtf8("select_tablewidget"));
        sizePolicy.setHeightForWidth(select_tablewidget->sizePolicy().hasHeightForWidth());
        select_tablewidget->setSizePolicy(sizePolicy);

        gridLayout_7->addWidget(select_tablewidget, 0, 0, 1, 1);

        configuration_tabwidget->addTab(selection_tab, QString());

        gridLayout_12->addWidget(configuration_tabwidget, 5, 0, 1, 5);

        horizontalSpacer_3 = new QSpacerItem(629, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_12->addItem(horizontalSpacer_3, 0, 6, 1, 1);

        analysis_tabwidget = new QTabWidget(centralwidget);
        analysis_tabwidget->setObjectName(QString::fromUtf8("analysis_tabwidget"));
        sizePolicy.setHeightForWidth(analysis_tabwidget->sizePolicy().hasHeightForWidth());
        analysis_tabwidget->setSizePolicy(sizePolicy);
        analysis_tabwidget->setMinimumSize(QSize(431, 245));
        analysis_tabwidget->setMaximumSize(QSize(16777215, 16777215));
        rf_tab = new QWidget();
        rf_tab->setObjectName(QString::fromUtf8("rf_tab"));
        gridLayout_13 = new QGridLayout(rf_tab);
        gridLayout_13->setObjectName(QString::fromUtf8("gridLayout_13"));
        analysis_groupbox = new QGroupBox(rf_tab);
        analysis_groupbox->setObjectName(QString::fromUtf8("analysis_groupbox"));
        QSizePolicy sizePolicy2(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(analysis_groupbox->sizePolicy().hasHeightForWidth());
        analysis_groupbox->setSizePolicy(sizePolicy2);
        analysis_groupbox->setMaximumSize(QSize(131, 51));
        gridLayout_2 = new QGridLayout(analysis_groupbox);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        trf_radiobutton = new QRadioButton(analysis_groupbox);
        trf_radiobutton->setObjectName(QString::fromUtf8("trf_radiobutton"));
        trf_radiobutton->setChecked(true);

        gridLayout_2->addWidget(trf_radiobutton, 0, 0, 1, 1);

        crf_radiobutton = new QRadioButton(analysis_groupbox);
        crf_radiobutton->setObjectName(QString::fromUtf8("crf_radiobutton"));
        crf_radiobutton->setChecked(false);

        gridLayout_2->addWidget(crf_radiobutton, 0, 1, 1, 1);


        gridLayout_13->addWidget(analysis_groupbox, 0, 0, 1, 1);

        refframe_groupbox = new QGroupBox(rf_tab);
        refframe_groupbox->setObjectName(QString::fromUtf8("refframe_groupbox"));
        refframe_groupbox->setMaximumSize(QSize(16777215, 16777215));
        gridLayout = new QGridLayout(refframe_groupbox);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        natural_radiobutton = new QRadioButton(refframe_groupbox);
        natural_radiobutton->setObjectName(QString::fromUtf8("natural_radiobutton"));
        QFont font;
        font.setPointSize(8);
        natural_radiobutton->setFont(font);
        natural_radiobutton->setChecked(true);

        gridLayout->addWidget(natural_radiobutton, 0, 0, 1, 1);

        mollweide_radiobutton = new QRadioButton(refframe_groupbox);
        mollweide_radiobutton->setObjectName(QString::fromUtf8("mollweide_radiobutton"));
        mollweide_radiobutton->setFont(font);

        gridLayout->addWidget(mollweide_radiobutton, 1, 0, 1, 3);

        mercator_radiobutton = new QRadioButton(refframe_groupbox);
        mercator_radiobutton->setObjectName(QString::fromUtf8("mercator_radiobutton"));
        mercator_radiobutton->setFont(font);

        gridLayout->addWidget(mercator_radiobutton, 2, 0, 1, 2);

        scale_arrows_spinbox = new QDoubleSpinBox(refframe_groupbox);
        scale_arrows_spinbox->setObjectName(QString::fromUtf8("scale_arrows_spinbox"));
        scale_arrows_spinbox->setMaximumSize(QSize(85, 25));
        scale_arrows_spinbox->setDecimals(1);
        scale_arrows_spinbox->setMinimum(0.100000000000000);
        scale_arrows_spinbox->setMaximum(100.000000000000000);
        scale_arrows_spinbox->setSingleStep(10.000000000000000);
        scale_arrows_spinbox->setValue(1.000000000000000);

        gridLayout->addWidget(scale_arrows_spinbox, 3, 0, 1, 1);

        label_5 = new QLabel(refframe_groupbox);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout->addWidget(label_5, 3, 1, 1, 3);

        ref_arrow_spinbox = new QDoubleSpinBox(refframe_groupbox);
        ref_arrow_spinbox->setObjectName(QString::fromUtf8("ref_arrow_spinbox"));
        ref_arrow_spinbox->setDecimals(1);
        ref_arrow_spinbox->setMinimum(0.100000000000000);
        ref_arrow_spinbox->setMaximum(10000.000000000000000);
        ref_arrow_spinbox->setSingleStep(100.000000000000000);
        ref_arrow_spinbox->setValue(1000.000000000000000);

        gridLayout->addWidget(ref_arrow_spinbox, 4, 0, 1, 2);

        ref_arrow_label = new QLabel(refframe_groupbox);
        ref_arrow_label->setObjectName(QString::fromUtf8("ref_arrow_label"));

        gridLayout->addWidget(ref_arrow_label, 4, 2, 1, 2);

        plot_ref_pushbutton = new QPushButton(refframe_groupbox);
        plot_ref_pushbutton->setObjectName(QString::fromUtf8("plot_ref_pushbutton"));

        gridLayout->addWidget(plot_ref_pushbutton, 5, 0, 1, 3);

        residuals_checkbox = new QCheckBox(refframe_groupbox);
        residuals_checkbox->setObjectName(QString::fromUtf8("residuals_checkbox"));

        gridLayout->addWidget(residuals_checkbox, 5, 3, 1, 1);


        gridLayout_13->addWidget(refframe_groupbox, 0, 1, 2, 1);

        transformation_groupbox = new QGroupBox(rf_tab);
        transformation_groupbox->setObjectName(QString::fromUtf8("transformation_groupbox"));
        transformation_groupbox->setEnabled(true);
        transformation_groupbox->setMaximumSize(QSize(217, 141));
        gridLayout_4 = new QGridLayout(transformation_groupbox);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        tx_checkbox = new QCheckBox(transformation_groupbox);
        tx_checkbox->setObjectName(QString::fromUtf8("tx_checkbox"));
        tx_checkbox->setChecked(true);

        gridLayout_4->addWidget(tx_checkbox, 0, 0, 1, 1);

        rx_checkbox = new QCheckBox(transformation_groupbox);
        rx_checkbox->setObjectName(QString::fromUtf8("rx_checkbox"));
        rx_checkbox->setChecked(true);

        gridLayout_4->addWidget(rx_checkbox, 0, 1, 1, 2);

        s_checkbox = new QCheckBox(transformation_groupbox);
        s_checkbox->setObjectName(QString::fromUtf8("s_checkbox"));
        s_checkbox->setChecked(true);

        gridLayout_4->addWidget(s_checkbox, 0, 3, 1, 1);

        ty_checkbox = new QCheckBox(transformation_groupbox);
        ty_checkbox->setObjectName(QString::fromUtf8("ty_checkbox"));
        ty_checkbox->setChecked(true);

        gridLayout_4->addWidget(ty_checkbox, 1, 0, 1, 1);

        ry_checkbox = new QCheckBox(transformation_groupbox);
        ry_checkbox->setObjectName(QString::fromUtf8("ry_checkbox"));
        ry_checkbox->setChecked(true);

        gridLayout_4->addWidget(ry_checkbox, 1, 1, 1, 2);

        transform_pushbutton = new QPushButton(transformation_groupbox);
        transform_pushbutton->setObjectName(QString::fromUtf8("transform_pushbutton"));

        gridLayout_4->addWidget(transform_pushbutton, 1, 3, 1, 1);

        tz_checkbox = new QCheckBox(transformation_groupbox);
        tz_checkbox->setObjectName(QString::fromUtf8("tz_checkbox"));
        tz_checkbox->setChecked(true);

        gridLayout_4->addWidget(tz_checkbox, 2, 0, 1, 1);

        rz_checkbox = new QCheckBox(transformation_groupbox);
        rz_checkbox->setObjectName(QString::fromUtf8("rz_checkbox"));
        rz_checkbox->setChecked(true);

        gridLayout_4->addWidget(rz_checkbox, 2, 1, 1, 2);

        def_only_checkbox = new QCheckBox(transformation_groupbox);
        def_only_checkbox->setObjectName(QString::fromUtf8("def_only_checkbox"));
        def_only_checkbox->setChecked(true);

        gridLayout_4->addWidget(def_only_checkbox, 3, 0, 1, 2);

        except_sh_checkbox = new QCheckBox(transformation_groupbox);
        except_sh_checkbox->setObjectName(QString::fromUtf8("except_sh_checkbox"));

        gridLayout_4->addWidget(except_sh_checkbox, 3, 2, 1, 2);

        arrow_width_spinbox = new QDoubleSpinBox(transformation_groupbox);
        arrow_width_spinbox->setObjectName(QString::fromUtf8("arrow_width_spinbox"));
        QSizePolicy sizePolicy3(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(arrow_width_spinbox->sizePolicy().hasHeightForWidth());
        arrow_width_spinbox->setSizePolicy(sizePolicy3);
        arrow_width_spinbox->setDecimals(1);
        arrow_width_spinbox->setMinimum(1.000000000000000);
        arrow_width_spinbox->setMaximum(20.000000000000000);

        gridLayout_4->addWidget(arrow_width_spinbox, 2, 3, 1, 1);


        gridLayout_13->addWidget(transformation_groupbox, 1, 0, 1, 1);

        analysis_tabwidget->addTab(rf_tab, QString());
        sta_tab = new QWidget();
        sta_tab->setObjectName(QString::fromUtf8("sta_tab"));
        gridLayout_8 = new QGridLayout(sta_tab);
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        sta_timeseries_checkbox = new QCheckBox(sta_tab);
        sta_timeseries_checkbox->setObjectName(QString::fromUtf8("sta_timeseries_checkbox"));
        sta_timeseries_checkbox->setChecked(true);

        gridLayout_8->addWidget(sta_timeseries_checkbox, 0, 0, 1, 3);

        horizontalSpacer = new QSpacerItem(29, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_8->addItem(horizontalSpacer, 1, 0, 1, 1);

        maxest_spinbox = new QDoubleSpinBox(sta_tab);
        maxest_spinbox->setObjectName(QString::fromUtf8("maxest_spinbox"));
        maxest_spinbox->setMaximumSize(QSize(85, 25));
        maxest_spinbox->setDecimals(0);
        maxest_spinbox->setMinimum(0.000000000000000);
        maxest_spinbox->setMaximum(1000000000.000000000000000);
        maxest_spinbox->setSingleStep(1.000000000000000);
        maxest_spinbox->setValue(50.000000000000000);

        gridLayout_8->addWidget(maxest_spinbox, 1, 1, 1, 2);

        label_6 = new QLabel(sta_tab);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout_8->addWidget(label_6, 1, 3, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(141, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_8->addItem(horizontalSpacer_2, 1, 4, 1, 1);

        sta_baselines_checkbox = new QCheckBox(sta_tab);
        sta_baselines_checkbox->setObjectName(QString::fromUtf8("sta_baselines_checkbox"));
        sta_baselines_checkbox->setCheckable(true);
        sta_baselines_checkbox->setChecked(false);

        gridLayout_8->addWidget(sta_baselines_checkbox, 2, 0, 1, 4);

        blreps_spinbox = new QSpinBox(sta_tab);
        blreps_spinbox->setObjectName(QString::fromUtf8("blreps_spinbox"));
        blreps_spinbox->setMinimum(30);
        blreps_spinbox->setMaximum(5000);

        gridLayout_8->addWidget(blreps_spinbox, 3, 1, 1, 1);

        label_7 = new QLabel(sta_tab);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout_8->addWidget(label_7, 3, 2, 1, 2);

        sta_helmert_checkbox = new QCheckBox(sta_tab);
        sta_helmert_checkbox->setObjectName(QString::fromUtf8("sta_helmert_checkbox"));
        sta_helmert_checkbox->setChecked(false);

        gridLayout_8->addWidget(sta_helmert_checkbox, 4, 0, 1, 4);

        gridLayout_11 = new QGridLayout();
        gridLayout_11->setObjectName(QString::fromUtf8("gridLayout_11"));
        tx_sta_checkbox = new QCheckBox(sta_tab);
        tx_sta_checkbox->setObjectName(QString::fromUtf8("tx_sta_checkbox"));
        tx_sta_checkbox->setChecked(true);

        gridLayout_11->addWidget(tx_sta_checkbox, 0, 0, 1, 1);

        ty_sta_checkbox = new QCheckBox(sta_tab);
        ty_sta_checkbox->setObjectName(QString::fromUtf8("ty_sta_checkbox"));
        ty_sta_checkbox->setChecked(true);

        gridLayout_11->addWidget(ty_sta_checkbox, 0, 1, 1, 1);

        tz_sta_checkbox = new QCheckBox(sta_tab);
        tz_sta_checkbox->setObjectName(QString::fromUtf8("tz_sta_checkbox"));
        tz_sta_checkbox->setChecked(true);

        gridLayout_11->addWidget(tz_sta_checkbox, 0, 2, 1, 1);

        rx_sta_checkbox = new QCheckBox(sta_tab);
        rx_sta_checkbox->setObjectName(QString::fromUtf8("rx_sta_checkbox"));
        rx_sta_checkbox->setChecked(false);

        gridLayout_11->addWidget(rx_sta_checkbox, 0, 3, 1, 1);

        ry_sta_checkbox = new QCheckBox(sta_tab);
        ry_sta_checkbox->setObjectName(QString::fromUtf8("ry_sta_checkbox"));
        ry_sta_checkbox->setChecked(false);

        gridLayout_11->addWidget(ry_sta_checkbox, 0, 4, 1, 1);

        rz_sta_checkbox = new QCheckBox(sta_tab);
        rz_sta_checkbox->setObjectName(QString::fromUtf8("rz_sta_checkbox"));
        rz_sta_checkbox->setChecked(false);

        gridLayout_11->addWidget(rz_sta_checkbox, 0, 5, 1, 1);

        s_sta_checkbox = new QCheckBox(sta_tab);
        s_sta_checkbox->setObjectName(QString::fromUtf8("s_sta_checkbox"));
        s_sta_checkbox->setChecked(true);

        gridLayout_11->addWidget(s_sta_checkbox, 0, 6, 1, 1);


        gridLayout_8->addLayout(gridLayout_11, 5, 1, 1, 4);

        plot_sta_pushbutton = new QPushButton(sta_tab);
        plot_sta_pushbutton->setObjectName(QString::fromUtf8("plot_sta_pushbutton"));

        gridLayout_8->addWidget(plot_sta_pushbutton, 6, 0, 1, 3);

        analysis_tabwidget->addTab(sta_tab, QString());
        eop_tab = new QWidget();
        eop_tab->setObjectName(QString::fromUtf8("eop_tab"));
        eops_listwidget = new QListWidget(eop_tab);
        eops_listwidget->setObjectName(QString::fromUtf8("eops_listwidget"));
        eops_listwidget->setGeometry(QRect(240, 50, 131, 121));
        sizePolicy.setHeightForWidth(eops_listwidget->sizePolicy().hasHeightForWidth());
        eops_listwidget->setSizePolicy(sizePolicy);
        c04_pushbutton = new QPushButton(eop_tab);
        c04_pushbutton->setObjectName(QString::fromUtf8("c04_pushbutton"));
        c04_pushbutton->setGeometry(QRect(240, 10, 111, 27));
        groupBox_5 = new QGroupBox(eop_tab);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        groupBox_5->setGeometry(QRect(9, 9, 211, 161));
        gridLayout_10 = new QGridLayout(groupBox_5);
        gridLayout_10->setObjectName(QString::fromUtf8("gridLayout_10"));
        yaxis_combobox = new QComboBox(groupBox_5);
        yaxis_combobox->addItem(QString());
        yaxis_combobox->addItem(QString());
        yaxis_combobox->addItem(QString());
        yaxis_combobox->setObjectName(QString::fromUtf8("yaxis_combobox"));

        gridLayout_10->addWidget(yaxis_combobox, 0, 1, 1, 1);

        xaxis_combobox = new QComboBox(groupBox_5);
        xaxis_combobox->addItem(QString());
        xaxis_combobox->addItem(QString());
        xaxis_combobox->addItem(QString());
        xaxis_combobox->setObjectName(QString::fromUtf8("xaxis_combobox"));

        gridLayout_10->addWidget(xaxis_combobox, 1, 1, 1, 1);

        timeformat_lineedit = new QLineEdit(groupBox_5);
        timeformat_lineedit->setObjectName(QString::fromUtf8("timeformat_lineedit"));

        gridLayout_10->addWidget(timeformat_lineedit, 3, 1, 1, 1);

        interpolateComboBox = new QComboBox(groupBox_5);
        interpolateComboBox->addItem(QString());
        interpolateComboBox->addItem(QString());
        interpolateComboBox->setObjectName(QString::fromUtf8("interpolateComboBox"));

        gridLayout_10->addWidget(interpolateComboBox, 2, 1, 1, 1);

        label = new QLabel(groupBox_5);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout_10->addWidget(label, 0, 0, 1, 1);

        label_3 = new QLabel(groupBox_5);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_10->addWidget(label_3, 1, 0, 1, 1);

        label_4 = new QLabel(groupBox_5);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_10->addWidget(label_4, 2, 0, 1, 1);

        label_8 = new QLabel(groupBox_5);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout_10->addWidget(label_8, 3, 0, 1, 1);

        plot_eop_pushbutton = new QPushButton(eop_tab);
        plot_eop_pushbutton->setObjectName(QString::fromUtf8("plot_eop_pushbutton"));
        plot_eop_pushbutton->setGeometry(QRect(10, 180, 141, 27));
        analysis_tabwidget->addTab(eop_tab, QString());
        src_tab = new QWidget();
        src_tab->setObjectName(QString::fromUtf8("src_tab"));
        analysis_tabwidget->addTab(src_tab, QString());

        gridLayout_12->addWidget(analysis_tabwidget, 4, 0, 1, 5);

        Analyzer->setCentralWidget(centralwidget);
        menubar = new QMenuBar(Analyzer);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1115, 21));
        Analyzer->setMenuBar(menubar);
        statusbar = new QStatusBar(Analyzer);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        Analyzer->setStatusBar(statusbar);

        retranslateUi(Analyzer);

        tabWidget->setCurrentIndex(1);
        configuration_tabwidget->setCurrentIndex(0);
        analysis_tabwidget->setCurrentIndex(2);


        QMetaObject::connectSlotsByName(Analyzer);
    } // setupUi

    void retranslateUi(QMainWindow *Analyzer)
    {
        Analyzer->setWindowTitle(QCoreApplication::translate("Analyzer", "MainWindow", nullptr));
        directory_pushbutton->setText(QCoreApplication::translate("Analyzer", "Directory", nullptr));
        clear_pushbutton->setText(QCoreApplication::translate("Analyzer", "Clear SNXs", nullptr));
        refresh_pushbutton->setText(QCoreApplication::translate("Analyzer", "Refresh Tree", nullptr));
        extern_checkbox->setText(QCoreApplication::translate("Analyzer", "ext Plot", nullptr));
        QTableWidgetItem *___qtablewidgetitem = session_tablewidget->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QCoreApplication::translate("Analyzer", "Name", nullptr));
        QTableWidgetItem *___qtablewidgetitem1 = session_tablewidget->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QCoreApplication::translate("Analyzer", "TRF", nullptr));
        QTableWidgetItem *___qtablewidgetitem2 = session_tablewidget->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QCoreApplication::translate("Analyzer", "CRF", nullptr));
        QTableWidgetItem *___qtablewidgetitem3 = session_tablewidget->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QCoreApplication::translate("Analyzer", "EOP", nullptr));
        QTableWidgetItem *___qtablewidgetitem4 = session_tablewidget->horizontalHeaderItem(4);
        ___qtablewidgetitem4->setText(QCoreApplication::translate("Analyzer", "Estimate", nullptr));
        QTableWidgetItem *___qtablewidgetitem5 = session_tablewidget->horizontalHeaderItem(5);
        ___qtablewidgetitem5->setText(QCoreApplication::translate("Analyzer", "Apriori", nullptr));
        QTableWidgetItem *___qtablewidgetitem6 = session_tablewidget->horizontalHeaderItem(6);
        ___qtablewidgetitem6->setText(QCoreApplication::translate("Analyzer", "SNX", nullptr));
        groupBox->setTitle(QCoreApplication::translate("Analyzer", "Initializer", nullptr));
        label_2->setText(QCoreApplication::translate("Analyzer", "Range:", nullptr));
        start_dateedit->setDisplayFormat(QCoreApplication::translate("Analyzer", "dd MMM yy", nullptr));
        clear_param_pushbutton->setText(QCoreApplication::translate("Analyzer", "Clear Param", nullptr));
        load_pushbutton->setText(QCoreApplication::translate("Analyzer", "Load SNX", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab), QCoreApplication::translate("Analyzer", "Tab 1", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QCoreApplication::translate("Analyzer", "Tab 2", nullptr));
        highlight_groupbox->setTitle(QCoreApplication::translate("Analyzer", "highlight series in plot", nullptr));
        highlight_pushbutton->setText(QCoreApplication::translate("Analyzer", "Highlight", nullptr));
        highlight_clear_pushbutton->setText(QCoreApplication::translate("Analyzer", "Clear all", nullptr));
        highlight_remove_checkbox->setText(QCoreApplication::translate("Analyzer", "remove old data", nullptr));
        or_site_radiobutton->setText(QCoreApplication::translate("Analyzer", "OR", nullptr));
        and_site_radiobutton->setText(QCoreApplication::translate("Analyzer", "AND ", nullptr));
        configuration_tabwidget->setTabText(configuration_tabwidget->indexOf(highligh_tab), QCoreApplication::translate("Analyzer", "Highlighting", nullptr));
        configuration_tabwidget->setTabText(configuration_tabwidget->indexOf(selection_tab), QCoreApplication::translate("Analyzer", "Selection", nullptr));
        analysis_groupbox->setTitle(QCoreApplication::translate("Analyzer", "Type", nullptr));
        trf_radiobutton->setText(QCoreApplication::translate("Analyzer", "TRF", nullptr));
        crf_radiobutton->setText(QCoreApplication::translate("Analyzer", "CRF", nullptr));
        refframe_groupbox->setTitle(QCoreApplication::translate("Analyzer", "Reference Frame Plots", nullptr));
        natural_radiobutton->setText(QCoreApplication::translate("Analyzer", "Natural", nullptr));
        mollweide_radiobutton->setText(QCoreApplication::translate("Analyzer", "Mollweide", nullptr));
        mercator_radiobutton->setText(QCoreApplication::translate("Analyzer", "Mercator", nullptr));
        label_5->setText(QCoreApplication::translate("Analyzer", "scale", nullptr));
        ref_arrow_label->setText(QCoreApplication::translate("Analyzer", "[microas]", nullptr));
        plot_ref_pushbutton->setText(QCoreApplication::translate("Analyzer", "Plot", nullptr));
        residuals_checkbox->setText(QCoreApplication::translate("Analyzer", "resids", nullptr));
        transformation_groupbox->setTitle(QCoreApplication::translate("Analyzer", "Transformation", nullptr));
        tx_checkbox->setText(QCoreApplication::translate("Analyzer", "tx", nullptr));
        rx_checkbox->setText(QCoreApplication::translate("Analyzer", "rx", nullptr));
        s_checkbox->setText(QCoreApplication::translate("Analyzer", "s", nullptr));
        ty_checkbox->setText(QCoreApplication::translate("Analyzer", "ty", nullptr));
        ry_checkbox->setText(QCoreApplication::translate("Analyzer", "ry", nullptr));
        transform_pushbutton->setText(QCoreApplication::translate("Analyzer", "Transform", nullptr));
        tz_checkbox->setText(QCoreApplication::translate("Analyzer", "tz", nullptr));
        rz_checkbox->setText(QCoreApplication::translate("Analyzer", "rz", nullptr));
        def_only_checkbox->setText(QCoreApplication::translate("Analyzer", "def only", nullptr));
        except_sh_checkbox->setText(QCoreApplication::translate("Analyzer", "except sh", nullptr));
        analysis_tabwidget->setTabText(analysis_tabwidget->indexOf(rf_tab), QCoreApplication::translate("Analyzer", "Reference Frames", nullptr));
        sta_timeseries_checkbox->setText(QCoreApplication::translate("Analyzer", "Time Series:", nullptr));
        label_6->setText(QCoreApplication::translate("Analyzer", "Max Estimation [mm]", nullptr));
        sta_baselines_checkbox->setText(QCoreApplication::translate("Analyzer", "Baseline Repeatabilities: ", nullptr));
        label_7->setText(QCoreApplication::translate("Analyzer", "Min Baseline No", nullptr));
        sta_helmert_checkbox->setText(QCoreApplication::translate("Analyzer", "Helmert Parameter:", nullptr));
        tx_sta_checkbox->setText(QCoreApplication::translate("Analyzer", "tx", nullptr));
        ty_sta_checkbox->setText(QCoreApplication::translate("Analyzer", "ty", nullptr));
        tz_sta_checkbox->setText(QCoreApplication::translate("Analyzer", "tz", nullptr));
        rx_sta_checkbox->setText(QCoreApplication::translate("Analyzer", "rx", nullptr));
        ry_sta_checkbox->setText(QCoreApplication::translate("Analyzer", "ry", nullptr));
        rz_sta_checkbox->setText(QCoreApplication::translate("Analyzer", "rz", nullptr));
        s_sta_checkbox->setText(QCoreApplication::translate("Analyzer", "s", nullptr));
        plot_sta_pushbutton->setText(QCoreApplication::translate("Analyzer", "Plot", nullptr));
        analysis_tabwidget->setTabText(analysis_tabwidget->indexOf(sta_tab), QCoreApplication::translate("Analyzer", "Station Analysis", nullptr));
        c04_pushbutton->setText(QCoreApplication::translate("Analyzer", "Load EOP Series", nullptr));
        groupBox_5->setTitle(QCoreApplication::translate("Analyzer", "plot setup", nullptr));
        yaxis_combobox->setItemText(0, QCoreApplication::translate("Analyzer", "total", nullptr));
        yaxis_combobox->setItemText(1, QCoreApplication::translate("Analyzer", "differences", nullptr));
        yaxis_combobox->setItemText(2, QCoreApplication::translate("Analyzer", "standard diviation", nullptr));

        xaxis_combobox->setItemText(0, QCoreApplication::translate("Analyzer", "time", nullptr));
        xaxis_combobox->setItemText(1, QCoreApplication::translate("Analyzer", "num of stations", nullptr));
        xaxis_combobox->setItemText(2, QCoreApplication::translate("Analyzer", "network volume", nullptr));

        timeformat_lineedit->setText(QCoreApplication::translate("Analyzer", "dd-MMM-yy", nullptr));
        interpolateComboBox->setItemText(0, QCoreApplication::translate("Analyzer", "subtrahend", nullptr));
        interpolateComboBox->setItemText(1, QCoreApplication::translate("Analyzer", "minuend", nullptr));

        label->setText(QCoreApplication::translate("Analyzer", "y-axis", nullptr));
        label_3->setText(QCoreApplication::translate("Analyzer", "x-axis", nullptr));
        label_4->setText(QCoreApplication::translate("Analyzer", "interpolate", nullptr));
        label_8->setText(QCoreApplication::translate("Analyzer", "time format", nullptr));
        plot_eop_pushbutton->setText(QCoreApplication::translate("Analyzer", "Plot", nullptr));
        analysis_tabwidget->setTabText(analysis_tabwidget->indexOf(eop_tab), QCoreApplication::translate("Analyzer", "EOP Analysis", nullptr));
        analysis_tabwidget->setTabText(analysis_tabwidget->indexOf(src_tab), QCoreApplication::translate("Analyzer", "Source Analysis", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Analyzer: public Ui_Analyzer {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ANALYZER_H
