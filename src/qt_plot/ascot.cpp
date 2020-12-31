#include "ascot.h"

Ascot::Ascot(libconfig::Config *cfg, string cfg_file) {
    
    // logger definition in order to redirect logs to GUI
    g_verbose = (loglevel) 7;
    ivg::Logger::qt_output = new QLineEdit();
    
    // setting main configuration pointer, no configuration initialized up to now
    _cfg = cfg;
    
    // setting up GUI
    _ui = new Ui::Ascot;
    _ui->setupUi(this);
    
    // define default configuration for the projection
    _trf_projection = {protype::naturalearth, 45.0, 45.0, false, true, 1.0, Qt::black, 1.0, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 8.0, 8.0, 20.0};
    _crf_projection = {protype::mollweide, 60.0, 30.0, true, false, 1.0, Qt::black, 1.0, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 8.0, 8.0, 20.0};
        
    _init_axis_layout();
    
    // set correct names of tabs
    _ui->main_tabwidget->setTabText(0,"Logger");
    _ui->main_tabwidget->setTabText(1,"Properties");
    _ui->main_tabwidget->setTabText(2,"Configuration");
    _ui->main_tabwidget->setTabEnabled(2,false);
    _ui->main_tabwidget->setTabText(3,"Results");
    _ui->main_tabwidget->setCurrentIndex(0);
    
    // plot tabwidget for results
    _ui->plot_tabwidget->clear();
    _ui->plot_tabwidget->setTabsClosable(true);
    
    _ui->refframe_tabwidget->clear();
    
    // setting default GUI-browser-loglevel to DETAILS and horizontal scrolling
    _ui->loglevel_combobox->setCurrentIndex(3); // 3 = DETAILS
    _ui->session_tablewidget->setSelectionMode(QAbstractItemView::SelectionMode::NoSelection);   
    _ui->output_plaintextedit->setLineWrapMode(QPlainTextEdit::NoWrap);
    if(!cfg_file.empty())
        _ui->cfg_lineedit->setText(QString::fromStdString(cfg_file));
    
    // be sure to set monospace font for log-browser
    QTextDocument *doc = _ui->output_plaintextedit->document();
    QFont font = doc->defaultFont();
    font.setPointSize(8);
    font.setFamily("Monospace");
    doc->setDefaultFont(font);
    
    // load the selected cfg-file if action is triggered
    connect(_ui->cfg_load_action, SIGNAL(triggered()), this, SLOT(_loadCfg()));
    // if something is logged with log<XXX>, it has to appear in the GUI
    connect(ivg::Logger::qt_output, SIGNAL(textChanged(QString)), this, SLOT(_appendLogInGUI(QString)));
    // double clicking the dbname leads to start reading the vgosdb
    connect(_ui->session_tablewidget, SIGNAL(cellDoubleClicked(const int, const int)), this, SLOT(_startLoadThread(const int, const int)));
    // connect cfg load pushbutton with loading of the cfg
    connect(_ui->cfg_pushbutton, SIGNAL(clicked()), this, SLOT(_loadCfg()));
    // by clicking the close button, close the plot-tabs
    connect(_ui->plot_tabwidget, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTabWidgetTab(int)));
    
    // show start of ASCOT
    log<INFO>("===================================================================");
    log<INFO>("        > > > > > ivg::ASCOT (interactive modus) < < < < <");
    
    this->show();
}// ...........................................................................
void Ascot::_loadCfg()
// ...........................................................................
{    
    QString cfg_file;
    if(_ui->cfg_lineedit->text().size() == 0)
        cfg_file = QFileDialog::getOpenFileName(this,tr("Config File"));
    else
        cfg_file = _ui->cfg_lineedit->text();
    
    try
    {
       _cfg->readFile( cfg_file.toStdString().c_str() );
       _ui->cfg_lineedit->setText(cfg_file);
       
       // load setup containing all pre-definied configuration
       log<INFO>("*** ") % cfg_file.toStdString() % " loaded successfully.";
       
       Setting& setup= _cfg->lookup( "setup" );
       // load ephemeris only once to avoid multiple jpl_init_ephemeris
       char nams[400][6];
       double vals[400];
       string ephfile_name = setup["definitions"]["ephemerides"][(const char*)setup["ephemerides"]];
       _ephem = jpl_init_ephemeris(ephfile_name.c_str(), nams, vals);
       log<INFO>("*** ") % ephfile_name % " ephemerides loaded successfully.";
        
        // loop over databases from main-cfg to put them into tablewidget
        for( int i=0; i<setup[ "sessions" ].getLength(); ++i )
        {
            string dbname  = setup[ "sessions" ][ i ][ "dbname" ];
//            string version = setup[ "sessions" ][ i ][ "version" ]; // for NGS
            
            _ui->session_tablewidget->setRowCount(_ui->session_tablewidget->rowCount()+1);
             int current_row = _ui->session_tablewidget->rowCount()-1;

            QTableWidgetItem *name_item = new QTableWidgetItem(QString::fromStdString(dbname));
            name_item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsUserCheckable);
            _ui->session_tablewidget->setItem(current_row,0,name_item);
        }
       
    }
    catch( libconfig::ParseException & err )
    {
        log<WARNING>("!!! Error while reading ") % cfg_file.toStdString();
        log<WARNING>("!!! libconfig::") % err.what() % " in " % err.getFile() % " at line " % err.getLine();
    }
    
}
// ...........................................................................
void Ascot::_startLoadThread(const int row, const int column)
// ...........................................................................
{
  
    QTableWidgetItem *clicked_item = _ui->session_tablewidget->item(row,column);
    
    // only if the dbname in column 0 is double-clicked
    if( column == 0 && clicked_item->backgroundColor() != Qt::green )
    {
        // highlight item yellow during import
       
        clicked_item->setBackgroundColor(Qt::yellow);
        
        // dbname of current double clicked session
        string dbname = clicked_item->text().toStdString();
        
        // initialize sessions based on setup, dbname and already initialized ephemerides
        ivg::Session S( &(_cfg->lookup( "setup" )), dbname, &_ephem );

        // set, connect and start loading thread
	 
        _thread = new Ascotizer(threadtype::LOAD, row, S);
        connect(_thread, SIGNAL(finished()), this, SLOT(_loadThreadFinished()));
	
        _thread->start();
	
    }
    else if(column == 0 && clicked_item->backgroundColor() == Qt::green)
    {
        _plotTrf(_loaded_sessions[row]);
    }
}
// ...........................................................................
void Ascot::_startAnalyseThread()
// ...........................................................................
{ 
    int row = _getClickedTableRow( sender() );
    
    if(row != -1 && !_thread->isRunning())
    {
        // set, connect and start loading thread
        _thread = new Ascotizer(threadtype::INIT, row, _loaded_sessions[row]);
        connect(_thread, SIGNAL(finished()), this, SLOT(_analyseThreadFinished()));
        _thread->start();
    }
}
// ...........................................................................
void Ascot::_startResidualsThread()
// ...........................................................................
{ 
    int row = _getClickedTableRow( sender() );
    
    if(row != -1 && !_thread->isRunning())
    {
        // set, connect and start loading thread
        _thread = new Ascotizer(threadtype::RESIDUALS, row, _analysed_sessions[row]);
        connect(_thread, SIGNAL(finished()), this, SLOT(_residualsThreadFinished()));
        _thread->start();
    }
}
// ...........................................................................
void Ascot::_startExportThread()
// ...........................................................................
{ 
    int row = _getClickedTableRow( sender() );
    
    if(row != -1 && !_thread->isRunning())
    {
        // set, connect and start loading thread
        _thread = new Ascotizer(threadtype::EXPORT, row, _analysed_sessions[row]);
        connect(_thread, SIGNAL(finished()), this, SLOT(_exportThreadFinished()));
        _thread->start();
    }
}
// .........................................................................
void Ascot::_loadThreadFinished()
{
    // when loading is finished, set cell green and allow analysis
    int idx = _thread->get_session_idx();
    _ui->session_tablewidget->item(idx,0)->setBackgroundColor(Qt::green);
    
    _loaded_sessions[idx] = _thread->get_session();
    
    QPushButton* analyse = new QPushButton("Analyse");
    _ui->session_tablewidget->setCellWidget(idx, 1, analyse);
    connect(analyse, SIGNAL(clicked()),this,SLOT(_startAnalyseThread()));
    
    _plotTrf(_loaded_sessions[idx]);
    _plotCrf(_loaded_sessions[idx]);
    _setEliminationTable(_loaded_sessions[idx]);

}
// .........................................................................
void Ascot::_analyseThreadFinished()
{
    // when analyse is finished, add residuals pushbutton
    int idx = _thread->get_session_idx();
    _analysed_sessions[idx] = _thread->get_session();
    
    QPushButton* residuals = new QPushButton("Residuals");
    _ui->session_tablewidget->setCellWidget(idx, 2, residuals);
    connect(residuals, SIGNAL(clicked()),this,SLOT(_startResidualsThread()));
    
    QPushButton* export_snx = new QPushButton("Export");
    _ui->session_tablewidget->setCellWidget(idx, 3, export_snx);
    connect(export_snx, SIGNAL(clicked()),this,SLOT(_startExportThread()));
}// .........................................................................
void Ascot::_residualsThreadFinished()
{
    int idx = _thread->get_session_idx();
    if(_latest_residuals[idx].size() == 0)
        // after residuals were calculated, the results can be viewed
        _latest_residuals[idx] = _thread->get_session().get_residuals();   
    
    // go to the results-tab and show results
    _ui->main_tabwidget->setCurrentIndex(3);
    _updateResidualTable(_latest_residuals[idx]);
    
    _latests_idx = idx;
}
// .........................................................................
void Ascot::_exportThreadFinished()
{
    QMessageBox msgBox;
    msgBox.setText("The SINEX file has been exported.");
    msgBox.exec();
}
// ...........................................................................
void Ascot::_appendLogInGUI(const QString &text)
{
    // compare level with level set in GUI(combobox) and append it
    int level = text.mid(0,1).toInt();
    if(level <= _ui->loglevel_combobox->currentIndex())
        _ui->output_plaintextedit->appendPlainText(text.mid(1,text.size()-1));
}
// ........................................................................
void Ascot::_setEliminationTable(ivg::Session &session)
// ........................................................................
{  
    QTableWidget * W = _ui->paraminfo_tablewidget;
    
    QStringList header;
    header << "Name" << "Type" << "# obs-all" << "# obs-removed" << "# obs-left" << "% obs-removed";
    W->setRowCount(1);
    W->setColumnCount(6);
    W->setHorizontalHeaderLabels(header);
    W->verticalHeader()->setVisible(false);
    W->setEditTriggers(QAbstractItemView::NoEditTriggers);
    W->setSelectionBehavior(QAbstractItemView::SelectRows);
    W->setSelectionMode(QAbstractItemView::SingleSelection);
    W->setShowGrid(false);
    W->setStyleSheet("QTableView {selection-background-color: red;}");
    
    // use trf and obsstats to show the statistics of the session (same alphabetical order!!)
    ivg::Trf *trf = session.get_trf_ptr();
    ivg::Matrix ori = session.get_obsstats()["TRF_ORI"];
    ivg::Matrix rem = session.get_obsstats()["TRF_REM"];
    ivg::Matrix per = session.get_obsstats()["TRF_PER"];
    
    int row = 0;
    int idx1 = 0;
    for(vector<ivg::Analysis_station>::iterator sta1 = trf->begin(); sta1 != trf->end(); sta1++)
    {
        int idx2 = idx1;
        for(vector<ivg::Analysis_station>::iterator sta2 = sta1; sta2 != trf->end(); sta2++)
        {
            QString station1 = QString::fromStdString(sta1->get_name(ivg::staname::ivs_name));
            QString station2 = QString::fromStdString(sta2->get_name(ivg::staname::ivs_name));

            QString name,type;
            if(idx1 == idx2)
            {
                name = station1;
                type = "Station";
            }
            else
            {
                name = station1+"-"+station2;
                type = "Baseline";
            }
            
            W->setItem(row, 0, new QTableWidgetItem(name));
            W->setItem(row, 1, new QTableWidgetItem(type));
            W->setItem(row, 2, new QTableWidgetItem(QString::number(ori(idx1,idx2))));
            W->setItem(row, 3, new QTableWidgetItem(QString::number(rem(idx1,idx2))));
            W->setItem(row, 4, new QTableWidgetItem(QString::number(ori(idx1,idx2)-rem(idx1,idx2))));
            W->setItem(row, 5, new QTableWidgetItem(QString::number(per(idx1,idx2),'f',1)));

            // only plot baseline if observations are left
            if(type == "Baseline" && ori(idx1,idx2)-rem(idx1,idx2) > 0)
            {
                ivg::Matrix llh;
                llh.append_rows(sta1->calc_lat_lon_h().transpose());
                llh.append_rows(sta2->calc_lat_lon_h().transpose());
    
                Plot *trf_plot = (Plot *) _ui->refframe_tabwidget->currentWidget();
                trf_plot->projection(_trf_projection, llh(":",1), llh(":",0), 
                                   { Qt::red, 1.5, Qt::DotLine, QCPGraph::lsLine, QCPScatterStyle::ssDiamond, name, 0.0, 0.0 }, 
                                   "Terrestrial Reference Frame"); 
            }
            
            idx2++;
            row++;
            W->setRowCount(W->rowCount()+1);
        }
    idx1++;
    }
    
    W->resizeColumnsToContents();
}
// ........................................................................
void Ascot::_updateResidualTable(map< delaytype , vector<Residual> > resids)
// ........................................................................
{                
    // create a new model
    QStandardItemModel *model = new QStandardItemModel(0,4,this);
    
    // set header of columns from model
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("Object")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("#in")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QString("#out")));
    model->setHorizontalHeaderItem(3, new QStandardItem(QString("wrms")));
    
    for(auto const &it: resids[delaytype::group] ) 
    {   
        // generate rows for table view from model
        QList<QStandardItem*> new_row;
        new_row.append(new QStandardItem(QString::fromStdString(it.name)));
        new_row.append(new QStandardItem(QString::number(it.idx_in.size())));
        new_row.append(new QStandardItem(QString::number(it.idx_out.size())));
        new_row.append(new QStandardItem(QString::number(it.wrms)));
        
        for(auto *item: new_row)
            item->setFont(QFont("Arial", 7));
        
        model->appendRow(new_row);
    }
    
    // attach the model to the view
    _ui->data_tableview->setSelectionBehavior(QAbstractItemView::SelectRows);
    _ui->data_tableview->setModel(model);
    
    // fit generated table to view
    _ui->data_tableview->resizeColumnsToContents();
    _ui->data_tableview->horizontalHeader()->setStretchLastSection(true);
    
    _ui->xaxis_combobox->setCurrentIndex(0);
    _ui->yaxis_combobox->setCurrentIndex(2);
    
    connect(_ui->plot_buttonbox, SIGNAL(accepted()), this, SLOT(_createResidualPlot()));
}
// ........................................................................
void Ascot::_createResidualPlot()
// ........................................................................
{
    // create new plot
    Plot *dialog = new Plot();
    _ui->plot_tabwidget->addTab(dialog, QString::fromUtf8("R"));
    _ui->plot_tabwidget->setCurrentIndex(_ui->plot_tabwidget->indexOf(dialog));
    
    Plot *histogram;
    if(_ui->histogram_checkbox->isChecked())
    {
        histogram = new Plot();
        _ui->plot_tabwidget->addTab(histogram, QString::fromUtf8("H"));
    }
    
    Plot *boxplot;
    if(_ui->boxplot_checkbox->isChecked())
    {
        boxplot = new Plot();
        _ui->plot_tabwidget->addTab(boxplot, QString::fromUtf8("B"));
    }
    
    QCustomPlot * plot = dialog->get_plot();
    
    string xaxis_combobox = _ui->xaxis_combobox->currentText().toLatin1().data();
    string yaxis_combobox = _ui->yaxis_combobox->currentText().toLatin1().data();
        
    // find specific axis layout for selected x/y axis
    Axislayout xaxis,yaxis;
    for(auto const &it : _axis_layouts ) 
    {
        if(it.name == xaxis_combobox)
            xaxis = it;
        if(it.name == yaxis_combobox)
            yaxis = it;
    }
    
    // set labels for plot
    plot->xAxis->setLabel(QString::fromStdString(xaxis.name + xaxis.unit));
    plot->yAxis->setLabel(QString::fromStdString(yaxis.name + yaxis.unit));
    
    plot->xAxis->setTickLabelType(xaxis.label_type);
    plot->xAxis->setDateTimeFormat(QString::fromStdString(xaxis.time_format));
    
    plot->yAxis->setTickLabelType(yaxis.label_type);
    plot->yAxis->setDateTimeFormat(QString::fromStdString(yaxis.time_format));
            
    QModelIndexList indexes = _ui->data_tableview->selectionModel()->selectedRows();
    
    map< delaytype , vector<Residual> > residuals = _latest_residuals[_latests_idx];
    
    for (int i = 0; i < indexes.count(); ++i)
    {
        // get new color for each new graph
        QColor plot_color = QColor(ivg::color_values.at(i).c_str());
        
        // detect which objects have been selected (e.g. WETTZELL and 0059+581)
        QModelIndex index = indexes.at(i);
        string selection =  _ui->data_tableview->model()->data(index,Qt::DisplayRole).toString().toStdString();
        
        for(auto const &it : residuals[delaytype::group] ) 
        { 
            if(selection == it.name)
            {
                //if outliers selected, add them
                vector< vector<int> > idx_inout = { it.idx_in };
                if( it.idx_out.size() > 0 && _ui->outlier_checkbox->isChecked() )
                    idx_inout.push_back(it.idx_out);
                                    
                int loop=1;
                for(auto const &idx : idx_inout ) 
                {
                    ivg::Matrix x = it.data.get_sub( idx, {xaxis.data_idx} ) * xaxis.fak;
                    ivg::Matrix y = it.data.get_sub( idx, {yaxis.data_idx} ) * yaxis.fak;
                    ivg::Matrix sig = it.data.get_sub( idx, {5} )*1e12;
                    
                    if(_ui->absolute_checkbox->isChecked())
                        y = y.abs();

                    if(_ui->errorbar_checkbox->isChecked())
                    {
                        dialog->plot_yerr( x, y, sig, { plot_color, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, selection.c_str(), 1.0 , 4.0  });
                        
                        //if outliers, plot the errorbars red
                        if(loop == 2)
                            plot->graph()->setErrorPen(QPen(QColor(Qt::red)));
                    }
                    else
                        dialog->plot_data( x, y, { plot_color, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, selection.c_str(), 1.0 , 4.0  });
                        
                     if(_ui->histogram_checkbox->isChecked())
                    {
                        int bins;
                        if(_ui->histogram_spinbox->value() == 0)
                            bins = (int) 1+3.3 * log(y.length())/log(10);
                        else
                            bins = _ui->histogram_spinbox->value();
                        
                        double spacing = histogram->histogram(y,bins,  { plot_color, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, "", 1.0 , 4.0  });
                        double mean = y.meanD();
                        double sigma = y.stdD();
                        ivg::Matrix keys(y.min(),(y.max()-y.min())/1000.0,y.max(),1);
                        ivg::Matrix tmp =  keys-mean;
                        tmp = tmp / sigma;
                        tmp = tmp.pow(2);
                        tmp = tmp * -0.5;
                        tmp = tmp.exp();
                        tmp = tmp / (sigma * 2.506628274631 );
                        tmp = tmp * y.length() * spacing;
                        
//                        cerr << "Mean: " << mean << " Std: " << sigma << endl;
                        histogram->plot_data(keys,tmp,  { plot_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, selection.c_str(), 1.0 , 4.0  });
                        histogram->get_plot()->replot();
                        histogram->get_plot()->rescaleAxes();
                    }    
                        
                    if(_ui->boxplot_checkbox->isChecked())
                    {
                        boxplot->boxplot(y,{ plot_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, "", 1.0 , 4.0  });
                    }
                    loop++;
                }
            }  
        }   
    }
    
    plot->rescaleAxes();
    plot->replot();
}
// ...........................................................................
void Ascot::_plotTrf(ivg::Session &session)
// ...........................................................................
{
    // get plot for the trf and plot each site individually with tooltip    
    Plot *pro_plot = new Plot();
    _ui->refframe_tabwidget->addTab(pro_plot, QString::fromUtf8("TRF"));
       
    ivg::Matrix lat;
    ivg::Matrix lon;
    vector<string> tooltips;
    for(auto &sta: (*session.get_trf_ptr()))
    {
        ivg::Matrix latlonh = sta.calc_lat_lon_h();
        lat.append_rows(latlonh(1));
        lon.append_rows(latlonh(0));
        
        stringstream tt;
        tt << sta.get_name(ivg::staname::ivs_name) << "      \n";
        tt << "DOMES:  " << sta.get_name(ivg::staname::domes_no) << "\n";
        tt << "CDP:    " << sta.get_name(ivg::staname::cdp) << "\n";
        tt << "Height: " << latlonh(2); 
        tooltips.push_back(tt.str());
    }
    
    pro_plot->projection(_trf_projection, lat, lon, { Qt::red, 5.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDiamond, "Sites", 4.0, 6.0 }, "Terrestrial Reference Frame", tooltips);
}

// ...........................................................................
void Ascot::_plotCrf(ivg::Session &session)
// ...........................................................................
{
    Plot *pro_plot = new Plot();
    _ui->refframe_tabwidget->addTab(pro_plot, QString::fromUtf8("CRF"));
    
    struct crfdata{ivg::Matrix ra; ivg::Matrix dec; ivg::Matrix d_ra; ivg::Matrix d_dec; vector<string> tooltips;};
    map<string,crfdata> crf_info;

    ivg::Crf crf = (*session.get_crf_ptr());

    ivg::Matrix ra,dec,d_ra,d_dec;
    for(vector<ivg::Source>::iterator src = crf.begin(); src != crf.end(); ++src)
    {
        ra.append_rows(src->get_ra0() - M_PI);
        dec.append_rows(src->get_dec0());

        string src_type="N";
        if(src->is_defining())
            src_type = "D";
        else if(src->is_special_handling())
            src_type = "S";

        stringstream tt;
        double s,sec;
        int h,m,deg,min;
        tt << "IVS: " << src->get_name(ivg::srcname::ivs) << "\n";
        tt << "IERS: " << src->get_name(ivg::srcname::iers) << "\n";
        tt << "ICRF: " << src->get_name(ivg::srcname::icrf) << "\n";
        tt << "RA: " << (src->get_ra0()  - M_PI)*(180.0/M_PI) << " / DEC: " << src->get_dec0()*(180.0/M_PI)<< "\n";
        src->get_position(h,m,s,deg,min,sec);
        tt << "RA: " << h << "h" << m << "m" << setprecision(7) << fixed << s << "s / DEC: ";
        tt << deg << QString::fromUtf8("\u00B0").toStdString() << min << "'" << setprecision(7) << fixed << sec << "''" << endl;
        tt << "ICRF2 Defining: " << src->is_defining() << endl;
        tt << "ICRF2 Special: " << src->is_special_handling();

        crf_info[src_type].ra.append_rows(src->get_ra0() - M_PI);
        crf_info[src_type].dec.append_rows(src->get_dec0());
        crf_info[src_type].tooltips.push_back(tt.str());
    }   
    
    // go through different types of sources (definings, special handlings, other)
    for(auto &srces: crf_info)
    {
        QString legend_name = QString::fromStdString(srces.first);
        QFormat format =  { Qt::blue, 1.5, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssDiamond, legend_name, 2.0 , 3.0 };

        if(srces.first == "D")
            format =  { Qt::green, 1.5, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssDisc, legend_name, 1.0 , 6.0 };
        else if(srces.first == "S")
            format =  { Qt::red, 1.5, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssCross, legend_name, 2.0 , 4.0 };

        pro_plot->projection(_crf_projection, srces.second.ra, srces.second.dec, format ,"Celestial Reference Frame", srces.second.tooltips );
    }
    
    _ui->refframe_tabwidget->setCurrentIndex(1);
    _ui->refframe_tabwidget->setCurrentIndex(0);
}
// ........................................................................
int Ascot::_getClickedTableRow( QObject *sender )
// ...........................................................................
{ 
    for(int row=0; row<_ui->session_tablewidget->rowCount(); row++){
        for(int col=0; col<_ui->session_tablewidget->columnCount(); col++){
            if(sender == _ui->session_tablewidget->cellWidget(row,col))
                return row;
        }
    }
    return -1;
}
// ...........................................................................
void Ascot::_init_axis_layout()
// ...........................................................................
{
    _axis_layouts.push_back({"Time - hh:mm:ss",0,"",1.0,QCPAxis::ltDateTime,"gb","hh:mm:ss"});
    _axis_layouts.push_back({"Time - dd-MMM-yy",0,"",1.0,QCPAxis::ltDateTime,"gb","dd-MMM-yy"});
    _axis_layouts.push_back({"Residuals",1,"[ps]",1e12,QCPAxis::ltNumber,"gb"});
    _axis_layouts.push_back({"Sigma Correlation",2,"[ps]",1e12,QCPAxis::ltNumber,"gb"});
    _axis_layouts.push_back({"Sigma LSA",3,"[ps]",1e12,QCPAxis::ltNumber,"gb"});
    _axis_layouts.push_back({"Sigma Residuals",5,"[ps]",1e12,QCPAxis::ltNumber,"gb"});
    _axis_layouts.push_back({"Redudancy Values",6,"",1.0,QCPAxis::ltNumber,"gb"});
    _axis_layouts.push_back({"Azimuth STA1",7,"[degree]",ivg::rad2d,QCPAxis::ltNumber,"gb"}); // ivg::rad2d
    _axis_layouts.push_back({"Elevation STA1",8,"[degree]",ivg::rad2d,QCPAxis::ltNumber,"gb"}); // ivg::rad2d
    _axis_layouts.push_back({"Azimuth STA2",9,"[degree]",ivg::rad2d,QCPAxis::ltNumber,"gb"});  // ivg::rad2d
    _axis_layouts.push_back({"Elevation STA2",10,"[degree]",ivg::rad2d,QCPAxis::ltNumber,"gb"}); // ivg::rad2d
    _axis_layouts.push_back({"Raydistance",11,"[km]",1e-3,QCPAxis::ltNumber,"gb"});
    _axis_layouts.push_back({"SNR X-band",12,"[-]",1.0,QCPAxis::ltNumber,"gb"});
    _axis_layouts.push_back({"SNR S-band",13,"[-]",1.0,QCPAxis::ltNumber,"gb"});
    
    for(auto const &it : _axis_layouts )
    {
        _ui->xaxis_combobox->addItem(QString::fromStdString(it.name));
        _ui->yaxis_combobox->addItem(QString::fromStdString(it.name));
    }
}
// ...........................................................................
void Ascot::closeTabWidgetTab(int index){
    _ui->plot_tabwidget->removeTab(index);
}
// ...........................................................................
