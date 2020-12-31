/********************************************************************************
** Form generated from reading UI file 'statistics_info.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_STATISTICS_INFO_H
#define UI_STATISTICS_INFO_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QTextBrowser>

QT_BEGIN_NAMESPACE

class Ui_Statistics_info
{
public:
    QGridLayout *gridLayout;
    QTextBrowser *textBrowser;

    void setupUi(QDialog *Statistics_info)
    {
        if (Statistics_info->objectName().isEmpty())
            Statistics_info->setObjectName(QStringLiteral("Statistics_info"));
        Statistics_info->resize(236, 164);
        Statistics_info->setWindowOpacity(1);
        gridLayout = new QGridLayout(Statistics_info);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        textBrowser = new QTextBrowser(Statistics_info);
        textBrowser->setObjectName(QStringLiteral("textBrowser"));
        textBrowser->setFocusPolicy(Qt::StrongFocus);
        textBrowser->setStyleSheet(QLatin1String("background-color: rgba(255, 255, 255, 0);\n"
"border-color: rgba(255, 255, 255, 0);"));

        gridLayout->addWidget(textBrowser, 0, 0, 1, 1);


        retranslateUi(Statistics_info);

        QMetaObject::connectSlotsByName(Statistics_info);
    } // setupUi

    void retranslateUi(QDialog *Statistics_info)
    {
        Statistics_info->setWindowTitle(QApplication::translate("Statistics_info", "Dialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Statistics_info: public Ui_Statistics_info {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_STATISTICS_INFO_H
