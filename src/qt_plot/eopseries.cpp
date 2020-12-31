#include "eopseries.h"

EopSeries::EopSeries() {
    
    _ui.setupUi(this);
    
    _ui.eops_tablewidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    _ui.eops_tablewidget->resizeColumnsToContents();
    
    _ui.eops_tablewidget->setColumnCount(3);
    
//    connect(ui.RB_png, SIGNAL(toggled(bool)), this, SLOT(_updateInterface(bool))); 
//    connect(ui.PB_file, SIGNAL(clicked()), this, SLOT(_openQFileDialog())); 
    
}

EopSeries::~EopSeries() {
    
    
}
void EopSeries::set_start_end_epochs(ivg::Date start, ivg::Date end)
{
    _ui.start_dateedit->setDate(QDate(start.get_int_year(), start.get_int_month(), start.get_int_day()));
    _ui.end_dateedit->setDate(QDate(end.get_int_year(), end.get_int_month(), end.get_int_day()));
}

void EopSeries::set_eop_file(string name, string type, string path)
{
    _ui.eops_tablewidget->setRowCount(_ui.eops_tablewidget->rowCount()+1);

    QTableWidgetItem *name_item1 = new QTableWidgetItem(QString::fromStdString(name));
    _ui.eops_tablewidget->setItem(_ui.eops_tablewidget->rowCount()-1, 0, name_item1);
    
    QTableWidgetItem *name_item2 = new QTableWidgetItem(QString::fromStdString(type));
    _ui.eops_tablewidget->setItem(_ui.eops_tablewidget->rowCount()-1, 1, name_item2);
    
    QTableWidgetItem *name_item3 = new QTableWidgetItem(QString::fromStdString(path));
    _ui.eops_tablewidget->setItem(_ui.eops_tablewidget->rowCount()-1, 2, name_item3);
        
    // make a nice appearance
    _ui.eops_tablewidget->resizeColumnsToContents();
}

ivg::Eop_series EopSeries::get_eop_series(string &name)
{
    
    QModelIndexList list = _ui.eops_tablewidget->selectionModel()->selectedRows();
    
    // only if one row is selected
    if(list.size() == 1)
    {       
        int sY,sM,sD,eY,eM,eD;
        _ui.start_dateedit->date().getDate(&sY,&sM,&sD);
        _ui.end_dateedit->date().getDate(&eY,&eM,&eD);
        
        // adjust year to correct century
        if(sY >= 2070)
            sY -= 100;
        
        if(eY >= 2070)
            eY -= 100;
        
        // save name by call by reference
        name = _ui.eops_tablewidget->item(list.at(0).row(),0)->text().toStdString();
        
        //create EOP series
        ivg::Eop_series eops( _ui.eops_tablewidget->item(list.at(0).row(),2)->text().toStdString(),
                              _ui.eops_tablewidget->item(list.at(0).row(),1)->text().toStdString(),
                              ivg::Date(sY,sM,sD), ivg::Date(eY,eM,eD));

        eops.init( _ui.interpol_combobox->currentText().toStdString(), _ui.ut1_zonal_checkbox->isChecked(),
                   _ui.hfocean_checkbox->isChecked(), "iers2010",
                   _ui.ut_lib_checkbox->isChecked(), _ui.pmnut_checkbox->isChecked(),
                   _ui.nutation_combobox->currentText().toStdString(),
                   _ui.pt_version_checkbox->currentText().toInt());
        
        
        return eops;
    }
    else if(list.size() == 0)
    {
        return ivg::Eop_series();
    }
    else
        throw runtime_error( "ivg::Eop_series EopSeries::get_eop_series(): More than one EOP-file selected.");
            
}    
    
