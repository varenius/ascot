/********************************************************************************
** Form generated from reading UI file 'ascot.ui'
**
** Created by: Qt User Interface Compiler version 5.12.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ASCOT_H
#define UI_ASCOT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableView>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QToolBox>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Ascot
{
public:
    QAction *cfg_load_action;
    QAction *eph_load_action;
    QAction *warning_action;
    QAction *result_action;
    QAction *info_action;
    QAction *detail_action;
    QAction *all_action;
    QWidget *centralwidget;
    QGridLayout *gridLayout_7;
    QLabel *cfg_label;
    QLineEdit *cfg_lineedit;
    QPushButton *cfg_pushbutton;
    QSpacerItem *horizontalSpacer;
    QComboBox *loglevel_combobox;
    QTabWidget *main_tabwidget;
    QWidget *log_tab;
    QGridLayout *gridLayout;
    QPlainTextEdit *output_plaintextedit;
    QWidget *properties_tab;
    QGridLayout *gridLayout_3;
    QTableWidget *paraminfo_tablewidget;
    QTabWidget *refframe_tabwidget;
    QWidget *tab_4;
    QWidget *tab_5;
    QWidget *config_tab;
    QToolBox *control_toolbox;
    QWidget *definitions_page;
    QComboBox *comboBox;
    QWidget *page_3;
    QWidget *page_2;
    QWidget *tab;
    QGridLayout *gridLayout_4;
    QTableView *data_tableview;
    QTabWidget *plot_tabwidget;
    QWidget *tab_2;
    QWidget *tab_3;
    QGridLayout *gridLayout_2;
    QLabel *xaxis_label;
    QComboBox *xaxis_combobox;
    QCheckBox *histogram_checkbox;
    QSpinBox *histogram_spinbox;
    QLabel *yaxis_label;
    QComboBox *yaxis_combobox;
    QCheckBox *boxplot_checkbox;
    QCheckBox *errorbar_checkbox;
    QCheckBox *outlier_checkbox;
    QCheckBox *absolute_checkbox;
    QDialogButtonBox *plot_buttonbox;
    QTableWidget *session_tablewidget;
    QMenuBar *menubar;
    QMenu *main_menu;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *Ascot)
    {
        if (Ascot->objectName().isEmpty())
            Ascot->setObjectName(QString::fromUtf8("Ascot"));
        Ascot->resize(777, 443);
        cfg_load_action = new QAction(Ascot);
        cfg_load_action->setObjectName(QString::fromUtf8("cfg_load_action"));
        eph_load_action = new QAction(Ascot);
        eph_load_action->setObjectName(QString::fromUtf8("eph_load_action"));
        warning_action = new QAction(Ascot);
        warning_action->setObjectName(QString::fromUtf8("warning_action"));
        warning_action->setCheckable(false);
        result_action = new QAction(Ascot);
        result_action->setObjectName(QString::fromUtf8("result_action"));
        result_action->setCheckable(false);
        result_action->setChecked(false);
        info_action = new QAction(Ascot);
        info_action->setObjectName(QString::fromUtf8("info_action"));
        detail_action = new QAction(Ascot);
        detail_action->setObjectName(QString::fromUtf8("detail_action"));
        all_action = new QAction(Ascot);
        all_action->setObjectName(QString::fromUtf8("all_action"));
        centralwidget = new QWidget(Ascot);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout_7 = new QGridLayout(centralwidget);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        cfg_label = new QLabel(centralwidget);
        cfg_label->setObjectName(QString::fromUtf8("cfg_label"));

        gridLayout_7->addWidget(cfg_label, 0, 0, 1, 1);

        cfg_lineedit = new QLineEdit(centralwidget);
        cfg_lineedit->setObjectName(QString::fromUtf8("cfg_lineedit"));

        gridLayout_7->addWidget(cfg_lineedit, 0, 1, 1, 1);

        cfg_pushbutton = new QPushButton(centralwidget);
        cfg_pushbutton->setObjectName(QString::fromUtf8("cfg_pushbutton"));

        gridLayout_7->addWidget(cfg_pushbutton, 0, 2, 1, 1);

        horizontalSpacer = new QSpacerItem(428, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_7->addItem(horizontalSpacer, 0, 3, 1, 1);

        loglevel_combobox = new QComboBox(centralwidget);
        loglevel_combobox->addItem(QString());
        loglevel_combobox->addItem(QString());
        loglevel_combobox->addItem(QString());
        loglevel_combobox->addItem(QString());
        loglevel_combobox->addItem(QString());
        loglevel_combobox->setObjectName(QString::fromUtf8("loglevel_combobox"));

        gridLayout_7->addWidget(loglevel_combobox, 0, 4, 1, 1);

        main_tabwidget = new QTabWidget(centralwidget);
        main_tabwidget->setObjectName(QString::fromUtf8("main_tabwidget"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(main_tabwidget->sizePolicy().hasHeightForWidth());
        main_tabwidget->setSizePolicy(sizePolicy);
        log_tab = new QWidget();
        log_tab->setObjectName(QString::fromUtf8("log_tab"));
        gridLayout = new QGridLayout(log_tab);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        output_plaintextedit = new QPlainTextEdit(log_tab);
        output_plaintextedit->setObjectName(QString::fromUtf8("output_plaintextedit"));
        QFont font;
        font.setFamily(QString::fromUtf8("Arial"));
        output_plaintextedit->setFont(font);
        output_plaintextedit->setReadOnly(true);

        gridLayout->addWidget(output_plaintextedit, 0, 0, 1, 1);

        main_tabwidget->addTab(log_tab, QString());
        properties_tab = new QWidget();
        properties_tab->setObjectName(QString::fromUtf8("properties_tab"));
        gridLayout_3 = new QGridLayout(properties_tab);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        paraminfo_tablewidget = new QTableWidget(properties_tab);
        paraminfo_tablewidget->setObjectName(QString::fromUtf8("paraminfo_tablewidget"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(70);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(paraminfo_tablewidget->sizePolicy().hasHeightForWidth());
        paraminfo_tablewidget->setSizePolicy(sizePolicy1);

        gridLayout_3->addWidget(paraminfo_tablewidget, 0, 0, 1, 1);

        refframe_tabwidget = new QTabWidget(properties_tab);
        refframe_tabwidget->setObjectName(QString::fromUtf8("refframe_tabwidget"));
        sizePolicy1.setHeightForWidth(refframe_tabwidget->sizePolicy().hasHeightForWidth());
        refframe_tabwidget->setSizePolicy(sizePolicy1);
        tab_4 = new QWidget();
        tab_4->setObjectName(QString::fromUtf8("tab_4"));
        refframe_tabwidget->addTab(tab_4, QString());
        tab_5 = new QWidget();
        tab_5->setObjectName(QString::fromUtf8("tab_5"));
        refframe_tabwidget->addTab(tab_5, QString());

        gridLayout_3->addWidget(refframe_tabwidget, 0, 1, 1, 1);

        main_tabwidget->addTab(properties_tab, QString());
        config_tab = new QWidget();
        config_tab->setObjectName(QString::fromUtf8("config_tab"));
        control_toolbox = new QToolBox(config_tab);
        control_toolbox->setObjectName(QString::fromUtf8("control_toolbox"));
        control_toolbox->setGeometry(QRect(10, 10, 961, 441));
        definitions_page = new QWidget();
        definitions_page->setObjectName(QString::fromUtf8("definitions_page"));
        definitions_page->setGeometry(QRect(0, 0, 961, 354));
        comboBox = new QComboBox(definitions_page);
        comboBox->setObjectName(QString::fromUtf8("comboBox"));
        comboBox->setGeometry(QRect(0, 10, 161, 25));
        control_toolbox->addItem(definitions_page, QString::fromUtf8("Page 1"));
        page_3 = new QWidget();
        page_3->setObjectName(QString::fromUtf8("page_3"));
        page_3->setGeometry(QRect(0, 0, 100, 30));
        control_toolbox->addItem(page_3, QString::fromUtf8("Seite"));
        page_2 = new QWidget();
        page_2->setObjectName(QString::fromUtf8("page_2"));
        page_2->setGeometry(QRect(0, 0, 100, 30));
        control_toolbox->addItem(page_2, QString::fromUtf8("Page 2"));
        main_tabwidget->addTab(config_tab, QString());
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        gridLayout_4 = new QGridLayout(tab);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        data_tableview = new QTableView(tab);
        data_tableview->setObjectName(QString::fromUtf8("data_tableview"));
        sizePolicy.setHeightForWidth(data_tableview->sizePolicy().hasHeightForWidth());
        data_tableview->setSizePolicy(sizePolicy);

        gridLayout_4->addWidget(data_tableview, 0, 0, 1, 1);

        plot_tabwidget = new QTabWidget(tab);
        plot_tabwidget->setObjectName(QString::fromUtf8("plot_tabwidget"));
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        plot_tabwidget->addTab(tab_2, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        plot_tabwidget->addTab(tab_3, QString());

        gridLayout_4->addWidget(plot_tabwidget, 0, 1, 2, 1);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        xaxis_label = new QLabel(tab);
        xaxis_label->setObjectName(QString::fromUtf8("xaxis_label"));
        QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(xaxis_label->sizePolicy().hasHeightForWidth());
        xaxis_label->setSizePolicy(sizePolicy2);

        gridLayout_2->addWidget(xaxis_label, 0, 0, 1, 1);

        xaxis_combobox = new QComboBox(tab);
        xaxis_combobox->setObjectName(QString::fromUtf8("xaxis_combobox"));
        sizePolicy2.setHeightForWidth(xaxis_combobox->sizePolicy().hasHeightForWidth());
        xaxis_combobox->setSizePolicy(sizePolicy2);
        xaxis_combobox->setMinimumSize(QSize(111, 25));

        gridLayout_2->addWidget(xaxis_combobox, 0, 1, 1, 1);

        histogram_checkbox = new QCheckBox(tab);
        histogram_checkbox->setObjectName(QString::fromUtf8("histogram_checkbox"));
        sizePolicy2.setHeightForWidth(histogram_checkbox->sizePolicy().hasHeightForWidth());
        histogram_checkbox->setSizePolicy(sizePolicy2);

        gridLayout_2->addWidget(histogram_checkbox, 0, 2, 1, 1);

        histogram_spinbox = new QSpinBox(tab);
        histogram_spinbox->setObjectName(QString::fromUtf8("histogram_spinbox"));
        sizePolicy2.setHeightForWidth(histogram_spinbox->sizePolicy().hasHeightForWidth());
        histogram_spinbox->setSizePolicy(sizePolicy2);
        histogram_spinbox->setMinimumSize(QSize(20, 25));
        histogram_spinbox->setMaximumSize(QSize(71, 25));
        histogram_spinbox->setMinimum(0);
        histogram_spinbox->setMaximum(1000);
        histogram_spinbox->setValue(0);

        gridLayout_2->addWidget(histogram_spinbox, 0, 3, 1, 1);

        yaxis_label = new QLabel(tab);
        yaxis_label->setObjectName(QString::fromUtf8("yaxis_label"));
        sizePolicy2.setHeightForWidth(yaxis_label->sizePolicy().hasHeightForWidth());
        yaxis_label->setSizePolicy(sizePolicy2);

        gridLayout_2->addWidget(yaxis_label, 1, 0, 1, 1);

        yaxis_combobox = new QComboBox(tab);
        yaxis_combobox->setObjectName(QString::fromUtf8("yaxis_combobox"));
        sizePolicy2.setHeightForWidth(yaxis_combobox->sizePolicy().hasHeightForWidth());
        yaxis_combobox->setSizePolicy(sizePolicy2);
        yaxis_combobox->setMinimumSize(QSize(111, 25));

        gridLayout_2->addWidget(yaxis_combobox, 1, 1, 1, 1);

        boxplot_checkbox = new QCheckBox(tab);
        boxplot_checkbox->setObjectName(QString::fromUtf8("boxplot_checkbox"));
        sizePolicy2.setHeightForWidth(boxplot_checkbox->sizePolicy().hasHeightForWidth());
        boxplot_checkbox->setSizePolicy(sizePolicy2);

        gridLayout_2->addWidget(boxplot_checkbox, 1, 2, 1, 1);

        errorbar_checkbox = new QCheckBox(tab);
        errorbar_checkbox->setObjectName(QString::fromUtf8("errorbar_checkbox"));
        sizePolicy2.setHeightForWidth(errorbar_checkbox->sizePolicy().hasHeightForWidth());
        errorbar_checkbox->setSizePolicy(sizePolicy2);
        errorbar_checkbox->setChecked(true);

        gridLayout_2->addWidget(errorbar_checkbox, 2, 1, 1, 1);

        outlier_checkbox = new QCheckBox(tab);
        outlier_checkbox->setObjectName(QString::fromUtf8("outlier_checkbox"));
        sizePolicy2.setHeightForWidth(outlier_checkbox->sizePolicy().hasHeightForWidth());
        outlier_checkbox->setSizePolicy(sizePolicy2);
        outlier_checkbox->setChecked(true);

        gridLayout_2->addWidget(outlier_checkbox, 2, 2, 1, 1);

        absolute_checkbox = new QCheckBox(tab);
        absolute_checkbox->setObjectName(QString::fromUtf8("absolute_checkbox"));
        sizePolicy2.setHeightForWidth(absolute_checkbox->sizePolicy().hasHeightForWidth());
        absolute_checkbox->setSizePolicy(sizePolicy2);

        gridLayout_2->addWidget(absolute_checkbox, 2, 3, 1, 1);

        plot_buttonbox = new QDialogButtonBox(tab);
        plot_buttonbox->setObjectName(QString::fromUtf8("plot_buttonbox"));
        sizePolicy2.setHeightForWidth(plot_buttonbox->sizePolicy().hasHeightForWidth());
        plot_buttonbox->setSizePolicy(sizePolicy2);
        plot_buttonbox->setStandardButtons(QDialogButtonBox::Ok);

        gridLayout_2->addWidget(plot_buttonbox, 1, 3, 1, 1);


        gridLayout_4->addLayout(gridLayout_2, 1, 0, 1, 1);

        main_tabwidget->addTab(tab, QString());

        gridLayout_7->addWidget(main_tabwidget, 1, 0, 1, 5);

        session_tablewidget = new QTableWidget(centralwidget);
        if (session_tablewidget->columnCount() < 6)
            session_tablewidget->setColumnCount(6);
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
        session_tablewidget->setObjectName(QString::fromUtf8("session_tablewidget"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(session_tablewidget->sizePolicy().hasHeightForWidth());
        session_tablewidget->setSizePolicy(sizePolicy3);
        session_tablewidget->setMaximumSize(QSize(16777215, 65));

        gridLayout_7->addWidget(session_tablewidget, 2, 0, 1, 5);

        Ascot->setCentralWidget(centralwidget);
        menubar = new QMenuBar(Ascot);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 777, 23));
        main_menu = new QMenu(menubar);
        main_menu->setObjectName(QString::fromUtf8("main_menu"));
        Ascot->setMenuBar(menubar);
        statusbar = new QStatusBar(Ascot);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        Ascot->setStatusBar(statusbar);

        menubar->addAction(main_menu->menuAction());
        main_menu->addAction(cfg_load_action);

        retranslateUi(Ascot);

        main_tabwidget->setCurrentIndex(3);
        refframe_tabwidget->setCurrentIndex(0);
        control_toolbox->setCurrentIndex(0);
        plot_tabwidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(Ascot);
    } // setupUi

    void retranslateUi(QMainWindow *Ascot)
    {
        Ascot->setWindowTitle(QApplication::translate("Ascot", "MainWindow", nullptr));
        cfg_load_action->setText(QApplication::translate("Ascot", "Load *.cfg", nullptr));
        eph_load_action->setText(QApplication::translate("Ascot", "Load Ephemerides", nullptr));
        warning_action->setText(QApplication::translate("Ascot", "WARNING", nullptr));
        result_action->setText(QApplication::translate("Ascot", "RESULT", nullptr));
        info_action->setText(QApplication::translate("Ascot", "INFO", nullptr));
        detail_action->setText(QApplication::translate("Ascot", "DETAIL", nullptr));
        all_action->setText(QApplication::translate("Ascot", "ALL", nullptr));
        cfg_label->setText(QApplication::translate("Ascot", "selected *.cfg:", nullptr));
        cfg_pushbutton->setText(QApplication::translate("Ascot", "Load", nullptr));
        loglevel_combobox->setItemText(0, QApplication::translate("Ascot", "WARNING", "1"));
        loglevel_combobox->setItemText(1, QApplication::translate("Ascot", "RESULT", "2"));
        loglevel_combobox->setItemText(2, QApplication::translate("Ascot", "INFO", "3"));
        loglevel_combobox->setItemText(3, QApplication::translate("Ascot", "DETAIL", "4"));
        loglevel_combobox->setItemText(4, QApplication::translate("Ascot", "ALL", "5"));

        main_tabwidget->setTabText(main_tabwidget->indexOf(log_tab), QApplication::translate("Ascot", "Tab 1", nullptr));
        refframe_tabwidget->setTabText(refframe_tabwidget->indexOf(tab_4), QApplication::translate("Ascot", "Tab 1", nullptr));
        refframe_tabwidget->setTabText(refframe_tabwidget->indexOf(tab_5), QApplication::translate("Ascot", "Tab 2", nullptr));
        main_tabwidget->setTabText(main_tabwidget->indexOf(properties_tab), QApplication::translate("Ascot", "Seite", nullptr));
        control_toolbox->setItemText(control_toolbox->indexOf(definitions_page), QApplication::translate("Ascot", "Page 1", nullptr));
        control_toolbox->setItemText(control_toolbox->indexOf(page_3), QApplication::translate("Ascot", "Seite", nullptr));
        control_toolbox->setItemText(control_toolbox->indexOf(page_2), QApplication::translate("Ascot", "Page 2", nullptr));
        main_tabwidget->setTabText(main_tabwidget->indexOf(config_tab), QApplication::translate("Ascot", "Tab 2", nullptr));
        plot_tabwidget->setTabText(plot_tabwidget->indexOf(tab_2), QApplication::translate("Ascot", "Tab 1", nullptr));
        plot_tabwidget->setTabText(plot_tabwidget->indexOf(tab_3), QApplication::translate("Ascot", "Tab 2", nullptr));
        xaxis_label->setText(QApplication::translate("Ascot", "X-Axis:", nullptr));
        histogram_checkbox->setText(QApplication::translate("Ascot", "Hist", nullptr));
        yaxis_label->setText(QApplication::translate("Ascot", "Y-Axis:", nullptr));
        boxplot_checkbox->setText(QApplication::translate("Ascot", "BoxPlot", nullptr));
        errorbar_checkbox->setText(QApplication::translate("Ascot", "Errorbar", nullptr));
        outlier_checkbox->setText(QApplication::translate("Ascot", "Outlier", nullptr));
        absolute_checkbox->setText(QApplication::translate("Ascot", "abs", nullptr));
        main_tabwidget->setTabText(main_tabwidget->indexOf(tab), QApplication::translate("Ascot", "Seite", nullptr));
        QTableWidgetItem *___qtablewidgetitem = session_tablewidget->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("Ascot", "DBNAME", nullptr));
        QTableWidgetItem *___qtablewidgetitem1 = session_tablewidget->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QApplication::translate("Ascot", "ANALYSE", nullptr));
        QTableWidgetItem *___qtablewidgetitem2 = session_tablewidget->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QApplication::translate("Ascot", "RESIDUALS", nullptr));
        QTableWidgetItem *___qtablewidgetitem3 = session_tablewidget->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QApplication::translate("Ascot", "RESULTS", nullptr));
        QTableWidgetItem *___qtablewidgetitem4 = session_tablewidget->horizontalHeaderItem(4);
        ___qtablewidgetitem4->setText(QApplication::translate("Ascot", "SNX", nullptr));
        main_menu->setTitle(QApplication::translate("Ascot", "Main", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Ascot: public Ui_Ascot {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ASCOT_H
