#include <./qpen.h>

#include "graphdialog.h"

// ...........................................................................
GraphDialog::GraphDialog(QWidget *parent) : QDialog(parent)
// ...........................................................................
{
    // initialize dialog
    _ui.setupUi(this);

    // add colordialog as tab to the dialog
    QColorDialog *colorDialog = new QColorDialog(this);
    colorDialog->setWindowFlags(Qt::Widget);
    colorDialog->setOptions(QColorDialog::DontUseNativeDialog | QColorDialog::NoButtons );
    colorDialog->setObjectName(QString::fromUtf8("Color"));
    _ui.tabWidget->addTab(colorDialog, QString());
    _ui.tabWidget->setTabText(_ui.tabWidget->indexOf(colorDialog), QString::fromUtf8("Color"));
    
    // color changes
    connect(colorDialog, SIGNAL(currentColorChanged(QColor)), this, SLOT(_setGraphColor(QColor)));
    
    // symbols
    connect(_ui.cross, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.plus, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.circle, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.disc, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.square, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.diamond, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.star, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.triangle, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.telescope, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    
    // line pen
    connect(_ui.solid, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.dash, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.dot, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.dashdot, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.dashdotdot, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    
    // line style
    connect(_ui.none, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.line, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.stepleft, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.stepright, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    connect(_ui.impulse, SIGNAL(toggled(bool)), this, SLOT(_setGraphStyle(bool)));
    
}
// ...........................................................................
void GraphDialog::_setGraphColor(QColor color)
// ...........................................................................
{
    // set the same color for everything
    // could be changed in the future
    _line_pen.setColor(color);
    _marker_pen.setColor(color);
    _marker_brush.setColor(color);
}
// ...........................................................................
void GraphDialog::setGraph(QCPGraph * graph)
// ...........................................................................
{
    // selected graph for changes
    _graph = graph;
    
    // properties of the line
    _line_style = _graph->lineStyle();
    _line_pen = _graph->pen();
    _line_width = _line_pen.width();

    // properties of the marker
    _marker_shape = _graph->scatterStyle().shape();
    _marker_pen = _graph->scatterStyle().pen();
    _marker_brush = _graph->scatterStyle().brush();
    _marker_size = _graph->scatterStyle().size();
    _marker_width = _graph->scatterStyle().pen().width();
    
    _updateInterface();
}
// ...........................................................................
void GraphDialog::updateGraph( QCPGraph * graph )
// ...........................................................................
{
    // in case of updating a different graph than _graph
    if(graph != NULL)
        _graph = graph;
    
    // fill the brush if desired
    if(_ui.filled->isChecked())
    {
        _marker_brush.setColor(_marker_pen.color());
        _marker_brush.setStyle(Qt::SolidPattern);
    }
    else
        _marker_brush.setStyle(Qt::NoBrush);
    
    // always make penstyle of marker as solidline
    _marker_pen.setStyle(Qt::PenStyle::SolidLine);
    _marker_pen.setWidth((int)_ui.markerwidth->value());
    
    
    if(_ui.telescope->isChecked())
        _graph->setScatterStyle(QCPScatterStyle(QPixmap("/opt/bakkari/source/qt_plot/radio.gif")));
    else
        _graph->setScatterStyle(QCPScatterStyle(_marker_shape, _marker_pen, _marker_brush, _ui.markersize->value()));
    
    
    // set new properties of the line
    _line_pen.setWidth((int)_ui.linewidth->value());
    _graph->setPen(_line_pen);
    _graph->setLineStyle(_line_style);
  
}
// ...........................................................................
void GraphDialog::_setGraphStyle(bool checked)
// ...........................................................................
{

    if(_ui.cross->isChecked())
        _marker_shape = QCPScatterStyle::ScatterShape::ssCross;
    else if(_ui.plus->isChecked())
        _marker_shape = QCPScatterStyle::ScatterShape::ssPlus;
    else if(_ui.circle->isChecked())
        _marker_shape = QCPScatterStyle::ScatterShape::ssCircle;
    else if(_ui.disc->isChecked())
        _marker_shape = QCPScatterStyle::ScatterShape::ssDisc;
    else if(_ui.square->isChecked())
        _marker_shape = QCPScatterStyle::ScatterShape::ssSquare;
    else if(_ui.diamond->isChecked())
        _marker_shape = QCPScatterStyle::ScatterShape::ssDiamond;
    else if(_ui.star->isChecked())
        _marker_shape = QCPScatterStyle::ScatterShape::ssStar;
    else if(_ui.triangle->isChecked())
       _marker_shape = QCPScatterStyle::ScatterShape::ssTriangle;
    else if(_ui.telescope->isChecked())
       _marker_shape = QCPScatterStyle::ScatterShape::ssPixmap;
    
    if(_ui.solid->isChecked())
        _line_pen.setStyle(Qt::PenStyle::SolidLine);
    else if(_ui.dash->isChecked())
        _line_pen.setStyle(Qt::PenStyle::DashLine);
    else if(_ui.dot->isChecked())
        _line_pen.setStyle(Qt::PenStyle::DotLine);
    else if(_ui.dashdot->isChecked())
        _line_pen.setStyle(Qt::PenStyle::DashDotLine);
    else if(_ui.dashdotdot->isChecked())
        _line_pen.setStyle(Qt::PenStyle::DashDotDotLine);
    
    if(_ui.none->isChecked())
        _line_style = QCPGraph::lsNone;
    else if(_ui.line->isChecked())
        _line_style = QCPGraph::lsLine;
    else if(_ui.stepleft->isChecked())
        _line_style = QCPGraph::lsStepLeft;
    else if(_ui.stepright->isChecked())
        _line_style = QCPGraph::lsStepRight;
    else if(_ui.impulse->isChecked())
        _line_style = QCPGraph::lsImpulse;
    
}
// ...........................................................................
void GraphDialog::_updateInterface()
// ...........................................................................
{
    // set interface based on the properties of the selected graph
    
    if(_marker_shape == QCPScatterStyle::ScatterShape::ssCross)
        _ui.cross->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssPlus )
        _ui.plus->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssCircle )
        _ui.circle->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssDisc )
        _ui.disc->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssSquare )
        _ui.square->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssDiamond )
        _ui.diamond->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssStar )
        _ui.star->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssTriangle )
        _ui.triangle->setChecked(true);
    else if(_marker_shape == QCPScatterStyle::ScatterShape::ssPixmap )
        _ui.telescope->setChecked(true);
    
    if(_line_pen.style() == Qt::PenStyle::SolidLine)
        _ui.solid->setChecked(true);
    else if(_line_pen.style() == Qt::PenStyle::DashLine)
        _ui.dash->setChecked(true);
    else if(_line_pen.style() == Qt::PenStyle::DotLine)
        _ui.dot->setChecked(true);
    else if(_line_pen.style() == Qt::PenStyle::DashDotLine)
        _ui.dashdot->setChecked(true);
    else if(_line_pen.style() == Qt::PenStyle::DashDotDotLine)
        _ui.dashdotdot->setChecked(true);
    
    if(_line_style == QCPGraph::lsNone)
        _ui.none->setChecked(true);
    else if(_line_style == QCPGraph::lsLine)
        _ui.line->setChecked(true);
    else if(_line_style == QCPGraph::lsStepLeft)
        _ui.stepleft->setChecked(true);
    else if(_line_style == QCPGraph::lsStepRight)
        _ui.stepright->setChecked(true);
    else if(_line_style == QCPGraph::lsImpulse)
        _ui.impulse->setChecked(true);
    
    if(_marker_brush.style() == Qt::NoBrush)
        _ui.filled->setChecked(false);
    else
        _ui.filled->setChecked(true);
        
    // adjust values from spinboxes to current graph properties
    _ui.markersize->setValue(_marker_size);
    _ui.markerwidth->setValue(_marker_width);
    _ui.linewidth->setValue(_line_width);
        
}
