/* 
 * File:   analyzer.h
 * Author: iddink
 *
 * Created on 14. Juli 2015, 10:57
 */

#ifndef ANALYZER_H
#define	ANALYZER_H

#include <string>
#include <iostream>
#include <sstream>
#include <QMainWindow>
#include <QDirModel>
#include <QFileInfo>
#include <QModelIndexList>
#include <QModelIndex>
#include <QStandardItem>
#include <QStringList>
#include <QDesktopServices>
#include <QProgressBar>
#include "ui_analyzer.h"


#include "station_analysis.h"
#include "masterfile.h"
#include "plot.h"
#include "eopseries.h"
#include "session_inout.h"
#include "transformation.h"
#include "sinex.h"

#include <cstdlib>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;

class Analyzer : public QMainWindow
{
  Q_OBJECT
  
  struct access{
      string full;
      string shorty;
  };
    
public:
    
    Analyzer(Setting *setup, QWidget *parent = 0);
    
    Analyzer(const Analyzer& orig);
    
    virtual ~Analyzer();
    
private:
    
    
    void plotSelectedTrf(Plot *plot, string ref_name, vector<string> other_sessions);
    void plotSelectedCrf(Plot *plot, string ref_name, vector<string> other_sessions);
    
    void updateTableWidget(string sess_name, string end_type);
    
    void setupParamTableView();
    
    QStringList _header_labels;
    
    Setting *_setup;
    
    string _latest_selected;
    vector<string> _latest_dir_files;
    
    int _latest_reference_row;
    
    ivg::Masterfile _masterfile;
    
    map<std::string, ivg::Session> _loaded_sessions;
    map<std::string, ivg::Sinex> _loaded_sinex;
    map<std::string, ivg::Eop_series> _loaded_eops;
    
    map<std::string, vector<ivg::Sinex> > _loaded_sinex_series;
    
    map<string, ivg::Matrix> _trf_residuals;
    map<string, ivg::Matrix> _crf_residuals;
    
    vector<ivg::Param> _params_in_tablewidget;
    
    int _highlight_color_cnt;
    
    Ui::Analyzer _ui;
    QDirModel *_model;
    
public slots:
    
    void objectTypeChanged();
    
    void clearLoadedSessions();
    void clearHighlightedItems();
    
    void updateDiffsGroupBox(const QModelIndex& index);
    
    
    void addObjectParameter();
    
    void updateSelectedReference(const int row, const int column);
    
    void selectedFileChanged(const QModelIndex &current, const QModelIndex &previous);
    
    void plotReferenceFrame();
    
    
    void loadEopSeries();
    
    void performStationAnalysis();
    
    void plotSelectedStations();
    void plotSelectedEops();
    void plotSelectedSources();
    
    void loadLatestSelectedSinex();
    
    void transformSelectedSessions();
    
    void closeTabWidgetTab(int index);
    void analysisTabChanged(int index);
    
    void highlightCurrentPlot();
    
private slots:
    
    bool stringCompare(const string & l, const string & r)                                                                                   
    {                                                                                                                                    
        return (l==r);                                                                                                                         
    }
    void _openDirectoryQFileDialog();
    void _refreshFolderTreeview();
    

};

#endif	/* ANALYZER_H */

