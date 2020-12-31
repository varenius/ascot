#ifndef EXPORTDIALOG_H
#define	EXPORTDIALOG_H

#include <iostream>
#include <QDialog>
#include <QPalette>
#include <QFileDialog>

#include "ui_export.h"

class ExportDialog : public QDialog
{
    Q_OBJECT
    
public:

    ExportDialog(QWidget *parent = 0);
    
    Ui::Export ui;
    
private slots:
    
    void _updateInterface(bool checked);
    void _openQFileDialog();

};



#endif	/* EXPORTDIALOG_H */

