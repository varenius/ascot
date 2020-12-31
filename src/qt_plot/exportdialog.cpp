#include "exportdialog.h"

ExportDialog::ExportDialog(QWidget *parent) : QDialog(parent)
{
    ui.setupUi(this);
    connect(ui.RB_png, SIGNAL(toggled(bool)), this, SLOT(_updateInterface(bool))); 
    connect(ui.PB_file, SIGNAL(clicked()), this, SLOT(_openQFileDialog())); 
}

void ExportDialog::_openQFileDialog()
{
    QString fileName = QFileDialog::getSaveFileName(this,tr("Export path"), "/home/", tr("Image Files (*.png *.jpg *.bmp *.pdf)"));
    ui.LE_pathfile->setText(fileName);
}

void ExportDialog::_updateInterface(bool checked)
{

    if(ui.RB_pdf->isChecked())
    {
        ui.CB_cospen->setCheckable(true);
        ui.LE_pdftitle->setEnabled(true);
        ui.LE_pdftitle->setStyleSheet("QLineEdit{background: white;}");
        ui.SB_quality->setReadOnly(true);
        ui.SB_scale->setReadOnly(true);
    }
    else if(ui.RB_png->isChecked()) 
    {
        ui.CB_cospen->setCheckable(false);
        ui.LE_pdftitle->setDisabled(true);
        ui.LE_pdftitle->setStyleSheet("QLineEdit{background: lightGray;}");
        ui.SB_scale->setReadOnly(false);
        ui.SB_quality->setReadOnly(false);
    }
    else if(ui.RB_bmp->isChecked()) 
    {
        ui.CB_cospen->setCheckable(false);
        ui.LE_pdftitle->setDisabled(true);
        ui.LE_pdftitle->setStyleSheet("QLineEdit{background: lightGray;}");
        ui.SB_scale->setReadOnly(false);
        ui.SB_quality->setReadOnly(true);
    }
    else if(ui.RB_jpg->isChecked()) 
    {
        ui.CB_cospen->setCheckable(false);
        ui.LE_pdftitle->setDisabled(true);
        ui.LE_pdftitle->setStyleSheet("QLineEdit{background: lightGray;}");
        ui.SB_scale->setReadOnly(false);
        ui.SB_quality->setReadOnly(false);
    }
}


