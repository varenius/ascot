#ifndef ASCOT_H
#define	ASCOT_H

#include <iostream>
#include <string>
#include <vector>

#include <./qfiledialog.h>
#include <QtCore/qthread.h>
#include <QMainWindow>
#include <QFileDialog>
#include <QTableWidget>
#include <QPushButton>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QMessageBox>

#include <cstdlib>
#include <libconfig.h++>

#include "ui_ascot.h"

#include "plot.h"
#include "logger.h"
#include "ascotizer.h"
#include "session.h"
#include "logger.h"
#include "session_inout.h"

using namespace std;

struct Axislayout
{
    string name;
    int data_idx;
    string unit;
    double fak;
    QCPAxis::LabelType label_type;
    string number_format;
    string time_format;
    
};

class Ascot  : public QMainWindow
{
  Q_OBJECT
    
public:
    
    Ascot(libconfig::Config *cfg, string cfg_file = "");
    
private:
    
    Ui::Ascot *_ui;
    libconfig::Config* _cfg;
    void *_ephem;
    Ascotizer *_thread;
    // saves loaded sessions related to a fixed index from tablewidget
    map<int,ivg::Session> _loaded_sessions, _analysed_sessions;
    map< int, map< delaytype , vector<Residual> > > _latest_residuals;
    int _latests_idx;
    
    vector<Axislayout> _axis_layouts;
    
    QCPProjection _trf_projection,_crf_projection;
    
private slots:
        
    void _loadCfg();
    
    void _startLoadThread(const int row, const int column);
    void _startAnalyseThread();
    void _startResidualsThread();
    void _startExportThread();
    
    void _loadThreadFinished();
    void _analyseThreadFinished();
    void _residualsThreadFinished();
    void _exportThreadFinished();
    
    void _updateResidualTable(map< delaytype , vector<Residual> > resids);
    void _createResidualPlot();
    void _appendLogInGUI(const QString &text);
    
    int _getClickedTableRow(QObject *sender);
    void _init_axis_layout();
    
    void _setEliminationTable(ivg::Session &session);
    
    void _plotTrf(ivg::Session &session);
    void _plotCrf(ivg::Session &session);
    
    void closeTabWidgetTab(int index);    
};

#endif	/* ASCOT_H */

