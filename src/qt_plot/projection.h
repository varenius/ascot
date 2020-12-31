/* 
 * File:   projection.h
 * Author: iddink
 *
 * Created on 8. September 2015, 09:12
 */

#ifndef PROJECTION_H
#define	PROJECTION_H

#include <iostream>
#include <QDialog>
#include <QColor>
#include <QPalette>
#include <QTabWidget>
#include <QFontDialog>
#include <QColorDialog>
#include <iostream>

#include "ui_projection.h"
#include "qcustomplot.h"

enum protype { mollweide, naturalearth, mercator, polarstereo, MAXPRO };

struct QCPProjection
{
    protype type;
    double lam_space;
    double phi_space; 
    bool plot_text;
    bool plot_coast;
    double coast_width;
    QColor coast_color;
    double grid_width;
    QColor grid_color;
    double grid_extension;
    QCPLineEnding::EndingStyle arrow_end;
    double arrow_width;
    double arrow_length;
    double arrow_scale;
    double arrow_value;
};

class Projection : public QDialog
{
    Q_OBJECT

public:
    
    Projection(QWidget *parent =0);
    
    void setProjection(QCPProjection *pro);
    
private slots:
    
    void _setConfigurationFromUI();
    
    void _setGraphStyle(bool checked);
    
private:
    
    void _updateInterface();
    
    Ui::Projection _ui;
    
    QCPProjection *_pro_ptr;
};

#endif	/* PROJECTION_H */

