#include "projection.h"
// ...........................................................................
Projection::Projection(QWidget *parent) : QDialog(parent)
// ...........................................................................
{
    // initialize dialog
    _ui.setupUi(this);
    
    // if user clicks OK, signal is emitted
    connect(_ui.buttonBox, SIGNAL(accepted()), this, SLOT(_setConfigurationFromUI()));
    
     // arrow end
    connect(_ui.qrb_disc, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.qrb_flat, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.qrb_line, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.qrb_none, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.qrb_spike, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    
}
// ...........................................................................
void Projection::setProjection(QCPProjection *pro)
// ...........................................................................
{
    _pro_ptr = pro;
    _updateInterface();
}
// ...........................................................................
void Projection::_setConfigurationFromUI()
// ...........................................................................
{
    if(_ui.mollweide->isChecked())
        _pro_ptr->type = protype::mollweide;
    else if(_ui.naturalearth->isChecked())
        _pro_ptr->type = protype::naturalearth;
    
    _pro_ptr->grid_width = _ui.gridwidth->value();
    _pro_ptr->grid_extension = _ui.gridextension->value();

    _pro_ptr->lam_space = _ui.lamspace->value();
    _pro_ptr->phi_space = _ui.phispace->value();
    
    _pro_ptr->plot_text = _ui.textlabels->isChecked();
    
}
// ...........................................................................
void Projection::_updateInterface()
// ...........................................................................
{
    // set interface based on the properties of the selected graph
    if(_pro_ptr->type == protype::mollweide)
        _ui.mollweide->setChecked(true);
    else if(_pro_ptr->type == protype::naturalearth)
        _ui.naturalearth->setChecked(true);
    
    _ui.gridextension->setValue(_pro_ptr->grid_extension);
    _ui.gridwidth->setValue(_pro_ptr->grid_width);
    
    _ui.lamspace->setValue(_pro_ptr->lam_space);
    _ui.phispace->setValue(_pro_ptr->phi_space);
    
    if(_pro_ptr->plot_text == true)
        _ui.textlabels->setChecked(true);
}

// ...........................................................................
void Projection::_setGraphStyle(bool checked)
// ...........................................................................
{
     if(_ui.qrb_disc->isChecked())
        _pro_ptr->arrow_end = QCPLineEnding::EndingStyle::esDisc;
    else if(_ui.qrb_flat->isChecked())
        _pro_ptr->arrow_end = QCPLineEnding::EndingStyle::esFlatArrow;
    else if(_ui.qrb_line->isChecked())
        _pro_ptr->arrow_end = QCPLineEnding::EndingStyle::esLineArrow;
    else if(_ui.qrb_none->isChecked())
        _pro_ptr->arrow_end = QCPLineEnding::EndingStyle::esNone;
    else if(_ui.qrb_spike->isChecked())
        _pro_ptr->arrow_end = QCPLineEnding::EndingStyle::esSpikeArrow;
}