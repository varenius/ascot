/* 
 * File:   analyzer.cpp
 * Author: iddink
 * 
 * Created on 14. Juli 2015, 10:57
 */

#include <qdirmodel.h>
#include "analyzer.h"

Analyzer::Analyzer(Setting *setup, QWidget *parent) : QMainWindow(parent)
{
    // in order to be able to read_snx we need the setup ( e.g. TRF initialization based on setup)
    _setup = setup;
    
    // initialize the masterfiles
    _masterfile = ivg::Masterfile((*_setup)["masterfiles"], ivg::mastertype::both);
    
    // setting up graphical user interface
    _ui.setupUi(this);

    // create and populate the model for the directories
    _model = new QDirModel(this);
        
    // enable modifying file system
    _model->setReadOnly(false);

    // setup sort criteria directory first, ignore case, and sort by name
    _model->setSorting(QDir::DirsFirst | QDir::IgnoreCase | QDir::Name);
    
    // allow horizontal scroll bar
    _ui.info_textbrowser->setLineWrapMode(QTextBrowser::NoWrap);
    
    QTextDocument *doc = _ui.info_textbrowser->document();
    QFont font = doc->defaultFont();
    font.setPointSize(9);
    font.setFamily("Monospace");
    doc->setDefaultFont(font);
    
    // no "latest reference" for CRF, TRF transformation up to now
    _latest_reference_row = -1;

    // QTreeView::setModel(QAbstractItemModel * model)
    _ui.folder_treeview->setModel(_model);
    _ui.folder_treeview->setColumnHidden(1,true);
    _ui.folder_treeview->setColumnHidden(2,true);
    
    // clean up tabwidget for plots
    _ui.tabWidget->clear();
    _ui.tabWidget->setTabsClosable(true);
    
    // set start value for the arrow scale
    _ui.scale_arrows_spinbox->setValue(1.0);
    _ui.ref_arrow_label->setText(QString::fromUtf8("[\u00B5as]"));
    
    // remove source-analysis tab because not implemented yet
    _ui.analysis_tabwidget->removeTab(3);

    
    // set starting point of foldertree
    QString path;
    if((bool)(*_setup).exists("foldertree") && (*_setup)["foldertree"] != "" )
        path = QString::fromStdString((*_setup)["foldertree"]);
    else
        path = QDir::homePath();
    
    QModelIndex index = _model->index(path);

    // set initial view of directory
    _ui.folder_treeview->expand(index);
    // scroll to the selected
    _ui.folder_treeview->scrollTo(index);
    // highlight the selected
    _ui.folder_treeview->setCurrentIndex(index);
    _ui.folder_treeview->setSelectionMode(QAbstractItemView::SingleSelection);
    // some further settings
    _ui.folder_treeview->resizeColumnToContents(0);
    _ui.param_tablewidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    _ui.param_tablewidget->setColumnCount(1);
    _header_labels << "Parameter";
    _ui.param_tablewidget->setHorizontalHeaderLabels(_header_labels);
    _ui.param_tablewidget->resizeColumnsToContents();
    
    // selection table for TRF-station-helmert transformation
    _ui.select_tablewidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    
    // widget holding all loaded sessions: no selection for the sessions!
    _ui.session_tablewidget->setSelectionMode(QAbstractItemView::SelectionMode::NoSelection);
    _ui.session_tablewidget->resizeColumnsToContents();
    
    // update reference session in case of double clicking (red/black)
    connect(_ui.session_tablewidget, SIGNAL(cellDoubleClicked(const int, const int)), this, SLOT(updateSelectedReference(const int, const int)));
    // open QFileDialog to change root view
    connect(_ui.directory_pushbutton, SIGNAL(clicked()), this, SLOT(_openDirectoryQFileDialog())); 
    // refresh folder tree for new results
    connect(_ui.refresh_pushbutton, SIGNAL(clicked()), this, SLOT(_refreshFolderTreeview())); 
    // if the selected file in the treeview changed, show the FILECOMMENT block from the sinex file
    connect(_ui.folder_treeview->selectionModel(), SIGNAL(currentRowChanged(const QModelIndex& , const QModelIndex&)), this, SLOT(selectedFileChanged(const QModelIndex& , const QModelIndex&)));
    // load the selected file if the LOAD button is clicked
    connect(_ui.load_pushbutton, SIGNAL(clicked()), this, SLOT(loadLatestSelectedSinex()));
    // clear all loaded sessions from data and listwidgets
    connect(_ui.clear_pushbutton, SIGNAL(clicked()), this, SLOT(clearLoadedSessions()));
    // clear tablewidget containing the double clicked parameter
    connect(_ui.clear_param_pushbutton, SIGNAL(clicked()), this, SLOT(objectTypeChanged()));
    // a loaded eop series can be set as reference eop series by double clicking on the entry in the listwidget
    connect(_ui.eops_listwidget, SIGNAL(doubleClicked(const QModelIndex&)), this, SLOT(updateDiffsGroupBox(const QModelIndex&)));
    // by clicking the plot reference frame button, plot hole content (TRF or CRF) of selected session
    connect(_ui.plot_ref_pushbutton, SIGNAL(clicked()), this, SLOT(plotReferenceFrame()));
    // by clicking the plot peries buttom, plot selected parameter from param_tablewidget
    connect(_ui.plot_eop_pushbutton, SIGNAL(clicked()), this, SLOT(plotSelectedEops()));
    connect(_ui.plot_sta_pushbutton, SIGNAL(clicked()), this, SLOT(plotSelectedStations()));
    connect(_ui.c04_pushbutton, SIGNAL(clicked()), this, SLOT(loadEopSeries()));
    connect(_ui.crf_radiobutton, SIGNAL(clicked()), this, SLOT(objectTypeChanged()));
    connect(_ui.trf_radiobutton, SIGNAL(clicked()), this, SLOT(objectTypeChanged()));
    
    // by clicking the transform button, all selected sessions (crf or trf), are compared
    connect(_ui.transform_pushbutton, SIGNAL(clicked()), this, SLOT(transformSelectedSessions()));
    // by clicking the close button, close the plot-tabs
    connect(_ui.tabWidget, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTabWidgetTab(int)));
    // update GUI in case of changing analysis-tabwidget
    connect(_ui.analysis_tabwidget, SIGNAL(currentChanged(int)), this, SLOT(analysisTabChanged(int)));
    
    // Highlight Series PushButtons
    connect(_ui.highlight_pushbutton, SIGNAL(clicked()), this, SLOT(highlightCurrentPlot()));
    connect(_ui.highlight_clear_pushbutton, SIGNAL(clicked()), this, SLOT(clearHighlightedItems()));
    
    // fill highlight-left-listwidget with groups from masterfile
    _highlight_color_cnt = 0;
    map<string, vector<string> > groups = _masterfile.get_groups();
    for(auto &grp: groups)
    {
        QListWidgetItem *item = new QListWidgetItem(QString::fromStdString(grp.first), _ui.highlight_listwidget);
        
        stringstream tt;
        for(auto &type: grp.second)
            tt << type << ",";
        
        item->setToolTip(QString::fromStdString(tt.str()));
    }
    
    _ui.highlight_listwidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    _ui.highlight_listwidget->setSelectionMode(QAbstractItemView::MultiSelection);
    
    
//    _ui.highlight_site_listwidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    _ui.highlight_site_listwidget->setSelectionMode(QAbstractItemView::MultiSelection);
    
    // set tab to reference-frame-tab
    _ui.analysis_tabwidget->setCurrentIndex(2);
    
    this->show();
    
    objectTypeChanged();
}
// ...........................................................................
Analyzer::Analyzer(const Analyzer& orig) {}
// ...........................................................................
Analyzer::~Analyzer(){}
// ...........................................................................
void Analyzer::_refreshFolderTreeview()
// ...........................................................................
{
    // DOES NOT WORK RIGHT NOW !!!!!!!
//    // refreshes the table treeview after new results have been computed
//    QFileInfo fileInfo = _model->fileInfo(_ui.folder_treeview->currentIndex());
//   
//    // create and populate the model for the directories
//    _model = new QDirModel(this);
//    // enable modifying file system
//    _model->setReadOnly(false);
//    // setup sort criteria directory first, ignore case, and sort by name
//    _model->setSorting(QDir::DirsFirst | QDir::IgnoreCase | QDir::Name);
//    
//    _model->index(fileInfo.absoluteFilePath());
////    QModelIndex index = _model->index(fileInfo.absoluteFilePath());
////    _ui.folder_treeview->clearFocus();
//    _ui.folder_treeview->setModel(_model);
//    
////    QModelIndex index = _model->index(fileInfo.absoluteFilePath());
//////            
//////    _ui.folder_treeview->
////
////    // set initial view of directory
////    _ui.folder_treeview->expand(index);
////    // scroll to the selected
////    _ui.folder_treeview->scrollTo(index);
////    // highlight the selected
////    _ui.folder_treeview->setCurrentIndex(index);
}
// ...........................................................................
void Analyzer::_openDirectoryQFileDialog()
{
    QString fileName = QFileDialog::getExistingDirectory(this,tr("Export path"));
    
    // set initial selection
    QModelIndex index = _model->index(fileName);

    // set initial view of directory
    _ui.folder_treeview->expand(index);
    // scroll to the selected
    _ui.folder_treeview->scrollTo(index);
    // highlight the selected
    _ui.folder_treeview->setCurrentIndex(index);
}
// ...........................................................................
void Analyzer::selectedFileChanged(const QModelIndex &current, const QModelIndex &previous)
// ...........................................................................
{
    // create stringstream containing information being showed in the textbrowser
    stringstream info_ss;
    QFileInfo fileInfo = _model->fileInfo(current);

    // save absolute file path to the latest selected sinex file
    _latest_selected = fileInfo.absoluteFilePath().toStdString();
    
    //clear old infos
    _ui.info_textbrowser->clear();
    
    // in case of selected folder
    if(fileInfo.isDir())
    {
        _latest_dir_files.clear();
        QDir dir(QString::fromStdString(_latest_selected));
        dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
        dir.setSorting(QDir::Name | QDir::Reversed);
        dir.setNameFilters(QStringList()<<"*.snx");

        info_ss << dir.count() << " SNX-files included (last modified): " << endl;
        info_ss << "------------------------------" << endl;
        QFileInfoList list = dir.entryInfoList();
        for (int i = 0; i < list.size(); ++i) {
            QFileInfo fileInfo = list.at(i);
            double mb = (double)(fileInfo.size()/1024.0/1024.0);
            info_ss << setfill(' ') << setw(22) << left << fileInfo.lastModified().toString("yyyy/MM/dd - hh:mm:ss").toStdString() << "| ";
            info_ss << setfill('0') << setw(5) << right << fixed << setprecision(2) << mb << "MByte | ";
            info_ss << left << fileInfo.fileName().toStdString() << "\n";
            _latest_dir_files.push_back(fileInfo.fileName().toStdString());
        }
        _ui.info_textbrowser->setText(QString::fromStdString(info_ss.str()));
    }
    // in case of selected file
    else
    {               
        // check for .err errorfiles and .snx SINEX-files
        string ending = _latest_selected.substr(_latest_selected.size()-4,4);
        
        string line;
        ifstream inStream(_latest_selected.c_str(), ios::in);
        
        if(ending == ".snx")
        {
            while (getline(inStream, line, '\n')){
                info_ss << line << endl;
                if(line.find("+SOLUTION/N")!=string::npos || line.find("+SOLUTION/MA")!=string::npos ){ 
                    _ui.info_textbrowser->setText(QString::fromStdString(info_ss.str()));
                    break;
                }        
            } 
        }
        else
        {
            while (getline(inStream, line, '\n'))
                info_ss << line << endl;
            
            _ui.info_textbrowser->setText(QString::fromStdString(info_ss.str()));
        }
    }
}
// ...........................................................................
void Analyzer::loadLatestSelectedSinex()
// ...........................................................................
{
   tictoc tim;
   tim.tic();
   
    //check if a dir is selected or not
    QFileInfo info(QString::fromStdString(_latest_selected));
    int snx_cnt = 1;
    int column, current_row;
    // save all ivs station names for highlight list widget
    vector<string> site_ivsname;
    if(info.isDir())
    {
        // in case of a dir, we need to load all sinex-files in this directory
        string name = _latest_selected.substr(_latest_selected.find_last_of("/")+1); // e.g. bkg2014a (equal to the version)
              
        // increase the name if name is already existing
        for(int i=0; i<=10; i++)
        {
            stringstream ss;
            ss << name << "-" << i;
            if(_loaded_sinex_series.find(ss.str()) == _loaded_sinex_series.end())
            {
                name = ss.str();
                break;
            } 
 	   
        }
        
        // get information about loading-range
        int sY,sM,sD,eY,eM,eD;
        _ui.start_dateedit->date().getDate(&sY,&sM,&sD);
        _ui.end_dateedit->date().getDate(&eY,&eM,&eD);
        
        // adjust year to correct century
        if(sY >= 2070)
            sY -= 100;
        
        if(eY >= 2070)
            eY -= 100;
        
        ivg::Date start(sY,sM,sD);
        ivg::Date end(eY,eM,eD);
        
//        ivg::Eop_series merged_series_apri, merged_series_esti;
        map<double, ivg::Eop_series> apri_eops,esti_eops; // fast alternative
        
        stringstream info_ss;
        
        QProgressBar *prog_bar = new QProgressBar(_ui.progressbar_widget);
        prog_bar->setMinimum(0);
        prog_bar->setMaximum(_latest_dir_files.size());
        prog_bar->setTextVisible(true);
        prog_bar->show();
                
        int sess_cnt=0;
        for(auto &sess: _latest_dir_files)
        {
            ivg::Date db_date(sess.substr(0,7),"YYMMMDD");
            if(db_date.get_double_mjd() >= start.get_double_mjd() && db_date.get_double_mjd() <= end.get_double_mjd())
            {
                // get all stations in this version-folder
                ivg::Sinex sinex = ivg::Sinex(_latest_selected+"/"+sess, false);
                vector<string> tmp_stations = sinex.get_trf(ivg::reftype::estimate).get_station_names(ivg::staname::ivs_name);
                site_ivsname.insert(site_ivsname.end(), tmp_stations.begin(), tmp_stations.end());
                
                //create based on the loaded sinex file series ONE merged eop series
                ivg::Eop_series tmp_apri = sinex.get_eop_series(ivg::reftype::apriori);
                ivg::Eop_series tmp_esti = sinex.get_eop_series(ivg::reftype::estimate);
                // check if the series is initialized or not
                if(tmp_apri.is_initialized())
                {
                    apri_eops[tmp_apri.get_data().operator()(0,0)] = tmp_apri;
                    esti_eops[tmp_esti.get_data().operator()(0,0)] = tmp_esti;
//                    merged_series_esti.merge(tmp_esti);
//                    merged_series_apri.merge(tmp_apri);
                }
                else
                    cerr << sess << " EOPs not successfully initialized" << endl;
                
                _loaded_sinex_series[name].push_back(sinex);
                snx_cnt++;
            }
            
            prog_bar->setValue(sess_cnt);
            sess_cnt++;
        }    
        
        // store both eop_series in vector
        vector<ivg::Eop_series> apri_vec, esti_vec;
        for(auto &A:apri_eops)
            apri_vec.push_back(A.second);
        
        for(auto &E: esti_eops)
            esti_vec.push_back(E.second);
        
        // use accumlating constructor for merging
        // (only works with eop_series containing one row)
        ivg::Eop_series merged_series_apri(apri_vec);
        ivg::Eop_series merged_series_esti(esti_vec);
        
        // save merged eop series for eop-analysis
        _loaded_eops[name+"_A"] = merged_series_apri;
        _loaded_eops[name+"_E"] = merged_series_esti;
        
        remove_duplicates(site_ivsname);
        
        _ui.session_tablewidget->setRowCount(_ui.session_tablewidget->rowCount()+1);
        current_row = _ui.session_tablewidget->rowCount()-1;

        QTableWidgetItem *name_item = new QTableWidgetItem(QString::fromStdString(name));
        name_item->setCheckState(Qt::CheckState::Unchecked);
        name_item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsUserCheckable);
        _ui.session_tablewidget->setItem(current_row,0,name_item);

        //create information about TRF,CRF,EOPs in loaded sinex file
        column = 1;
        vector<ivg::objtype> objects = {ivg::objtype::trf, ivg::objtype::crf, ivg::objtype::eop};
        for(auto &obj: objects)
        {
            QTableWidgetItem *item;
            if(obj == ivg::objtype::trf)
            {
                item = new QTableWidgetItem(QString::number(site_ivsname.size()));
                if(site_ivsname.size() > 0)
                    item->setBackgroundColor(Qt::green);
                else
                    item->setBackgroundColor(Qt::red);
            }
            else if(obj == ivg::objtype::eop)
            {
                item = new QTableWidgetItem(QString::number(merged_series_esti.length()));
                if(merged_series_esti.length() > 0)
                    item->setBackgroundColor(Qt::green);
                else
                    item->setBackgroundColor(Qt::red);
            }
            else if(obj == ivg::objtype::crf)
            {
                item = new QTableWidgetItem(QString("0"));
                item->setBackgroundColor(Qt::red);
            }
            
            _ui.session_tablewidget->setItem(current_row, column, item);
            column++;
        }
        
        // close progress bar when finished
        prog_bar->close();
    }
    else
    {

      string tmp= _latest_selected.substr(_latest_selected.find_last_of("/")+1);
     
      if (tmp.find("_psd") != std::string::npos || tmp.find("_PSD") != std::string::npos ||tmp.find("-psd") != std::string::npos || tmp.find("-PSD") != std::string::npos) {

	// PSD correction file
	 if(_latest_reference_row != -1)
	   {
	     string ref_name = _ui.session_tablewidget->item(_latest_reference_row,0)->text().toStdString();
	     _psd_file[ref_name]=_latest_selected;
	   }
	 else {
	   std::cout << "No reference selected for applying the PSD data to" << endl;
	 }
	
      } else {

        // if the LOAD button is clicked, load the selected sinex file
        string name = _latest_selected.substr(_latest_selected.find_last_of("/")+1,_latest_selected.find_first_of("_")-_latest_selected.find_last_of("/")-1);
        
        // increase the name if name is already existing
        for(int i=0; i<=10; i++)
        {
            stringstream ss;
            ss << name << "-" << i;
            if(_loaded_sinex.find(ss.str()) == _loaded_sinex.end())
            {
                _loaded_sinex[ss.str()] = ivg::Sinex(_latest_selected, true);
                name = ss.str();
                break;
            }
        }

        _ui.session_tablewidget->setRowCount(_ui.session_tablewidget->rowCount()+1);
        current_row = _ui.session_tablewidget->rowCount()-1;

        QTableWidgetItem *name_item = new QTableWidgetItem(QString::fromStdString(name));
        name_item->setCheckState(Qt::CheckState::Unchecked);
        name_item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsUserCheckable);
        _ui.session_tablewidget->setItem(current_row,0,name_item);


        //create information about TRF,CRF,EOPs in loaded sinex file
        column = 1;
        vector<ivg::objtype> objects = {ivg::objtype::trf, ivg::objtype::crf, ivg::objtype::eop};
        for(auto &obj: objects)
        {
            int cnt = _loaded_sinex[name].get_num_objects(obj);

            QTableWidgetItem *item = new QTableWidgetItem(QString::number(cnt));

            item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);

            if(cnt > 0 )
                item->setBackgroundColor(Qt::green);
            else
                item->setBackgroundColor(Qt::red);

            if(obj == ivg::objtype::trf)
            {
//                site_ivsname = _loaded_sinex[name].get_trf(ivg::reftype::estimate).get_station_names(ivg::staname::ivs_name);
//                string all_stations;
//                for(auto &sta: stations)
//                    all_stations = all_stations+" "+sta;
//                
//                item->setToolTip(QString::fromStdString(all_stations));
            }
            else if(obj == ivg::objtype::eop && cnt > 0)
            {
                _loaded_eops[name+"_A"] = _loaded_sinex[name].get_eop_series(ivg::reftype::apriori);
                _loaded_eops[name+"_E"] = _loaded_sinex[name].get_eop_series(ivg::reftype::estimate);
            }
            
            _ui.session_tablewidget->setItem(current_row, column, item);

            column++;
        }
      }
    }
    // add new signal-slots to each new pushbutton
    QPushButton* add_estimate = new QPushButton("Add");
    _ui.session_tablewidget->setCellWidget(current_row, column, add_estimate);
    connect(add_estimate, SIGNAL(clicked()),this,SLOT(addObjectParameter()));

    QPushButton* add_apriori = new QPushButton("Add");
    _ui.session_tablewidget->setCellWidget(current_row, column+1, add_apriori);
    connect(add_apriori, SIGNAL(clicked()),this,SLOT(addObjectParameter()));

    QTableWidgetItem *snx_item = new QTableWidgetItem(QString::number(snx_cnt));
    snx_item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    _ui.session_tablewidget->setItem(current_row, 6, snx_item);
    
    // make a nice appearance
    _ui.session_tablewidget->resizeColumnsToContents();
    
    // fill highlight-left-listwidget with site used in the sinex files
    
    // get all different types of names 
    vector< map<ivg::staname,string> > nscodes = ivg::parser::nscodes_parser( site_ivsname, ivg::staname::ivs_name, (*_setup)["nscodes"]);
    for(auto &site: nscodes)
    {
        QListWidgetItem *item = new QListWidgetItem(QString::fromStdString(site[ivg::staname::ivs_name]), _ui.highlight_site_listwidget);
        item->setToolTip(QString::fromStdString(site[ivg::staname::lettercode]));
    }
    
    objectTypeChanged();
}
// ...........................................................................
void Analyzer::updateDiffsGroupBox(const QModelIndex& index)
// ...........................................................................
{ 
    string eops_name = index.data().toString().toStdString();
    
    // put eop series to param_tablewidget
    updateTableWidget(eops_name, "");
}
// ...........................................................................
void Analyzer::addObjectParameter()
// ...........................................................................
{
    // an object is a TRF, a CRF or EOPs
    // add for pushed "add" button related parameter to param_tablewidget
    for(int row=0; row<_ui.session_tablewidget->rowCount(); row++){
        for(int col=0; col<_ui.session_tablewidget->columnCount(); col++){
 
            if(sender() == _ui.session_tablewidget->cellWidget(row,col))
            {
                string session =  _ui.session_tablewidget->item(row,0)->text().toStdString();

                string end_type;
                if(col == 4)
                   end_type = "_E";
                else if(col == 5)
                   end_type = "_A";

                // add the parameters right now to the param_tablewidget
                updateTableWidget(session, end_type);
            }
        }
    }    
}
// ...........................................................................
void Analyzer::updateSelectedReference(const int row, const int column)
// ...........................................................................
{      
    // only in case of double clicked session name (first column) in session_tablewidget
    // defines the reference for the transformation in TRF and CRF
    if( column == 0 )
    {
        QTableWidgetItem *clicked_item = _ui.session_tablewidget->item(row,column);
        
        if( row != _latest_reference_row && _latest_reference_row == -1)
        {
            clicked_item->setTextColor(Qt::red);
            _latest_reference_row = row;
        }
        else if( row != _latest_reference_row && _latest_reference_row != -1)
        {
            _ui.session_tablewidget->item(_latest_reference_row,0)->setTextColor(Qt::black);
            clicked_item->setTextColor(Qt::red);
            _latest_reference_row = row;
        }
        else if( row == _latest_reference_row)
        {
            clicked_item->setTextColor(Qt::black);
            _latest_reference_row = -1;
        }
        
        if(_latest_reference_row != -1)
        {
            // if a sinex file contains stations, show them in the selection table for helmert parameter defining stations
            _ui.select_tablewidget->setColumnCount(1);
            _ui.select_tablewidget->setRowCount(0);
            _ui.select_tablewidget->setHorizontalHeaderLabels(QStringList("Sites"));
            
            // ref_name is the current double clicked and red marked session
            string ref_name = _ui.session_tablewidget->item(_latest_reference_row,0)->text().toStdString();
            ivg::Trf ref_trf = _loaded_sinex[ref_name].get_trf(ivg::reftype::estimate);
            
            for(auto &sta: ref_trf)
            {
                _ui.select_tablewidget->setRowCount(_ui.select_tablewidget->rowCount()+1);
                int current_row = _ui.select_tablewidget->rowCount()-1;

                QTableWidgetItem *name_item = new QTableWidgetItem(QString::fromStdString(sta.get_name(ivg::staname::ivs_name)));
                name_item->setCheckState(Qt::CheckState::Checked);
                name_item->setFlags(Qt::ItemIsEnabled|Qt::ItemIsUserCheckable);
                _ui.select_tablewidget->setItem(current_row,0,name_item);
            }
            _ui.select_tablewidget->resizeColumnsToContents();
        }
    }
}
// ...........................................................................
void Analyzer::updateTableWidget(string sess_name, string end_type)
// ...........................................................................
{   
    // depending on the selected analysis (CRF,TRF,EOP) it need to be focused on different parameters
    // struct defining parametertype including orders
    struct selection{vector<ivg::paramtype> types; vector<int> orders;};
        
    // only if its not "Reference Frames"-Tab
    if(_ui.analysis_tabwidget->currentIndex() != 0)
    {
        selection selected;
        if(_ui.analysis_tabwidget->currentIndex() == 1)
            selected = { {ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz}, {0,1} };
        else if(_ui.analysis_tabwidget->currentIndex() == 2 )
            selected = { {ivg::paramtype::ut1, ivg::paramtype::xpo, ivg::paramtype::ypo, ivg::paramtype::nutx, ivg::paramtype::nuty}, {0,1} };
        else if(_ui.analysis_tabwidget->currentIndex() == 3 )
            selected = { {ivg::paramtype::ra,ivg::paramtype::dec}, {0,1} };

        // add new column for new session
        _ui.param_tablewidget->setColumnCount(_ui.param_tablewidget->columnCount()+1);
        int current_col = _ui.param_tablewidget->columnCount()-1;

        // create header labels base on session name and end_type (_E or _A)
        _header_labels << QString::fromStdString(sess_name+end_type);

        vector<ivg::Param> parameter = _loaded_sinex[sess_name].get_parameter_vetor();

        // in case of empty parameter vector, we need to take a look at the series-map
        // and loop over all sinex files within the series
        if(parameter.empty())
        {
            for(auto &sinex: _loaded_sinex_series[sess_name])
            {
                vector<ivg::Param> tmp = sinex.get_parameter_vetor();
                parameter.insert(parameter.end(), tmp.begin(), tmp.end());
            }
        }
        
        // if parameter vector is still empty, we need to focus on the loaded eop series
        if(parameter.empty())
        {              
            for(int i=0; i<=_loaded_eops[sess_name].length(); i++)
            {
                parameter.push_back(ivg::Param(ivg::paramtype::xpo, "EOP", ivg::Date(), i, 0));
                parameter.push_back(ivg::Param(ivg::paramtype::xpo, "EOP", ivg::Date(), i, 1));
                parameter.push_back(ivg::Param(ivg::paramtype::ypo, "EOP", ivg::Date(), i, 0));
                parameter.push_back(ivg::Param(ivg::paramtype::ypo, "EOP", ivg::Date(), i, 1));
                parameter.push_back(ivg::Param(ivg::paramtype::ut1, "EOP", ivg::Date(), i, 0));
                parameter.push_back(ivg::Param(ivg::paramtype::ut1, "EOP", ivg::Date(), i, 1));
                parameter.push_back(ivg::Param(ivg::paramtype::nutx, "EOP", ivg::Date(), i, 0));
                parameter.push_back(ivg::Param(ivg::paramtype::nuty, "EOP", ivg::Date(), i, 0));
            }
        }

        // go through the parameter list and show information in param_tablewidget
        for(auto &param_iter: parameter)
        {
            // check on which parameter should be focused on, e.g. sources with ra and dec, defined by the radiobuttons
            if(param_iter.is_type(selected.types, selected.orders))
            {
                    stringstream ss;
                    ss << param_iter.get_name() << "_" << param_iter.get_typename() << "_" << param_iter.get_order();
                    string compare_name = ss.str();

                    bool add_row = true;
                    for(int row=0; row<_ui.param_tablewidget->rowCount(); row++)
                    {
                        if(compare_name == _ui.param_tablewidget->item(row,0)->text().toStdString())
                        {
                            if(_ui.param_tablewidget->item(row,current_col))
                            {
                                int increment = _ui.param_tablewidget->item(row,current_col)->text().toInt() + 1;
                                _ui.param_tablewidget->setHorizontalHeaderLabels(_header_labels);
                                _ui.param_tablewidget->item(row,current_col)->setText(QString::number(increment));
                            }
                            else
                            {
                                QTableWidgetItem *num_item = new QTableWidgetItem(QString::number(1));
                                num_item->setBackgroundColor(Qt::green);
                                _ui.param_tablewidget->setItem(row, current_col, num_item);
                            }
                            add_row = false;
                        }
                    }

                    if(add_row)
                    {
                        // add new row for each new station, source or eop
                        _ui.param_tablewidget->setRowCount(_ui.param_tablewidget->rowCount()+1);

                        QTableWidgetItem *name_item = new QTableWidgetItem(QString::fromStdString(compare_name));
                        _ui.param_tablewidget->setItem(_ui.param_tablewidget->rowCount()-1, 0, name_item);
                        _params_in_tablewidget.push_back(param_iter);

                        QTableWidgetItem *num_item = new QTableWidgetItem(QString::number(1));
                        num_item->setBackgroundColor(Qt::green);
                        _ui.param_tablewidget->setItem(_ui.param_tablewidget->rowCount()-1, current_col, num_item);
                    }
            }
        }

        // make a nice appearance
        _ui.param_tablewidget->resizeColumnsToContents();
    }
}
// ...........................................................................
void Analyzer::plotReferenceFrame()
// ...........................................................................
{    
    // create new plot
    Plot *plot = new Plot();
    // in case of checked checkbox, the plots will be generated extern
    // to enable good axis scaling for plot-export
    if(!_ui.extern_checkbox->isChecked())
    {
        _ui.tabWidget->addTab(plot, QString::fromUtf8("RF"));
        _ui.tabWidget->setCurrentIndex(_ui.tabWidget->indexOf(plot));
    }
    else
        plot->show();
    
    // only if an item is selected, we can get a name
    string ref_name = "";
    if(_latest_reference_row != -1)
        ref_name = _ui.session_tablewidget->item(_latest_reference_row,0)->text().toStdString();
    
    // get selected sessions
    vector<string> other_sessions;
    for(int row=0; row<_ui.session_tablewidget->rowCount(); row++)
    {
        if(_ui.session_tablewidget->item(row,0)->checkState() == Qt::CheckState::Checked)
            other_sessions.push_back(_ui.session_tablewidget->item(row,0)->text().toStdString());
    }
   
   // only plot if at least one sessions is selected
   if(other_sessions.size()>0)
   {
        if(_ui.trf_radiobutton->isChecked())
            plotSelectedTrf(plot, ref_name, other_sessions);
        else if(_ui.crf_radiobutton->isChecked())
            plotSelectedCrf(plot, ref_name, other_sessions);
   }
   else
       throw runtime_error("void Analyzer::plotReferenceFrame(): No sessions selected with checkbox.");
}
// ...........................................................................
void Analyzer::plotSelectedTrf(Plot *plot, string ref_name, vector<string> other_sessions)
// ...........................................................................
{ 
    // loop over selected other_sessions and plot them
    int color_cnt=1;
    for(auto &other_name: other_sessions)
    {
         ivg::Trf other_trf = _loaded_sinex[other_name].get_trf(ivg::reftype::estimate);

         ivg::Matrix lat, d_lat;
         ivg::Matrix lon, d_lon;
         vector<string> tooltips;
         for(vector<ivg::Analysis_station>::iterator sta = other_trf.begin(); sta != other_trf.end(); ++sta)
         {
             ivg::Matrix latlonh = sta->calc_lat_lon_h();

            // if residuals are selected and a transformation has been done (= residuals not empty)
            if(_ui.residuals_checkbox->isChecked() && _trf_residuals[ref_name+"_"+other_name].rows() != 0)
            {
                ivg::Matrix residuals = _trf_residuals[ref_name+"_"+other_name];
                int sta_index = sta - other_trf.begin() ;

                for(int i=0; i<residuals.rows(); i++)
                {
                    if((int)residuals(i,1) == sta_index)
                    {
                        ivg::Matrix d_xyz = residuals.get_sub(i,2,i,4);
                        ivg::Analysis_station sta_moved = (*sta);
                        sta_moved.move_station(d_xyz.transpose());
                        ivg::Matrix latlonh_moved = sta_moved.calc_lat_lon_h();
                        ivg::Matrix d_latlonh = latlonh_moved - latlonh;
                        d_lat.append_rows(d_latlonh(0)); // 1.56789571971343e-09 rad = 8.9833807e-08 degree = 1cm 
                        d_lon.append_rows(d_latlonh(1)); // 1.56789571971343e-09 rad = 8.9833807e-08 degree = 1cm 
                    }
                }
            }

             lat.append_rows(latlonh(1));
             lon.append_rows(latlonh(0));

             stringstream tt;
             tt << "IVS:    " << sta->get_name(ivg::staname::ivs_name) << "\n";
             tt << "DOMES:  " << sta->get_name(ivg::staname::domes_no) << "\n";
             tt << "CDP:    " << sta->get_name(ivg::staname::cdp) << "\n";
             tt << "Height: " << latlonh(2); 
             tooltips.push_back(tt.str());
        }

        // add residuals as second column
        if(lat.rows() == d_lat.rows())
        {
           lat.append_cols(d_lat);
           lon.append_cols(d_lon);
        }

        protype pro_type = protype::mollweide;
        if(_ui.mollweide_radiobutton->isChecked())
            pro_type = protype::mollweide;
        else if(_ui.natural_radiobutton->isChecked())
            pro_type = protype::naturalearth;
        else if(_ui.mercator_radiobutton->isChecked())
            pro_type = protype::mercator;

        double scale_arrows = _ui.scale_arrows_spinbox->value();
        double ref_arrow_value = _ui.ref_arrow_spinbox->value();
        
        // used for grid, coast and residual arrows
        double line_width = _ui.arrow_width_spinbox->value(); 
        
        QCPProjection projection = {pro_type, 45.0, 45.0, false, true, line_width/1.8, Qt::black, line_width/1.2, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 5.0*(1+(line_width/2.8)), 5.0*(1+(line_width/2.8)), scale_arrows, ref_arrow_value};
        
        QColor marker_color = QColor(ivg::color_values.at(color_cnt).c_str());
        // plot last one with neutral grey (combi)
        if(color_cnt == other_sessions.size())
            marker_color = QColor("#ACACAC");
        plot->projection(projection, lat, lon, { marker_color, line_width, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, QString::fromStdString(other_name), 1.0, 10.0 }, "", tooltips);

        color_cnt += 1;         
    }    
}
// ...........................................................................
void Analyzer::plotSelectedCrf(Plot *plot, string ref_name, vector<string> other_sessions)
// ...........................................................................
{ 
        struct crfdata{ivg::Matrix ra; ivg::Matrix dec; ivg::Matrix d_ra; ivg::Matrix d_ra_cos; ivg::Matrix d_dec; vector<string> tooltips;};
        map< string, map<string,crfdata> > crfs;
        for(auto &name: other_sessions)
        { 
            ivg::Crf crf = _loaded_sinex[name].get_crf(ivg::reftype::estimate);
            
            ivg::Matrix ra,dec,d_ra, d_ra_cos,d_dec;
            for(vector<ivg::Source>::iterator src = crf.begin(); src != crf.end(); ++src)
            {
                ra.append_rows(src->get_ra0() - M_PI);
                dec.append_rows(src->get_dec0());
                
                string src_type="N";
                if(src->is_defining())
                    src_type = "D";
                else if(src->is_special_handling())
                    src_type = "S";
                
                if(_ui.residuals_checkbox->isChecked() && _crf_residuals[ref_name+"_"+name].rows() != 0)
                {
                    ivg::Matrix residuals = _crf_residuals[ref_name+"_"+name];
                    int src_index = src - crf.begin() ;
                    
                    for(int i=0; i<residuals.rows(); i++)
                    {
                        if((int)residuals(i,1) == src_index)
                        {                            
                            crfs[name][src_type].d_ra.append_rows(residuals(i,2));
                            crfs[name][src_type].d_ra_cos.append_rows(residuals(i,2)*cos(residuals(i,3)));
                            crfs[name][src_type].d_dec.append_rows(residuals(i,3));
                            
                            // RAW SYSTEM DIFFERENCES (columns 5 & 6) INSTEAD OF TRANSFORMED (columns 2 & 3)
//                            crfs[name][src_type].d_ra.append_rows(residuals(i,5));
//                            crfs[name][src_type].d_ra_cos.append_rows(residuals(i,5)*cos(residuals(i,6)));
//                            crfs[name][src_type].d_dec.append_rows(residuals(i,6));
                        }
                    }

                }

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


                crfs[name][src_type].ra.append_rows(src->get_ra0() - M_PI);
                crfs[name][src_type].dec.append_rows(src->get_dec0());
                crfs[name][src_type].tooltips.push_back(tt.str());
            }     
        }
        
//            if(_ui.local_sources_pushbutton->isChecked())
//                plotSelectedSources();

            // in case of plotting residuals we need an additional plot for source differences
            Plot *diff_plot = new Plot();
            Plot *dradec_plot = new Plot();
            int diff_subplot = diff_plot->add_rowplot("Y","X");
            if(_ui.residuals_checkbox->isChecked() && !_ui.extern_checkbox->isChecked())
            {
                _ui.tabWidget->addTab(diff_plot, QString::fromUtf8("DecDiff"));
                _ui.tabWidget->addTab(dradec_plot, QString::fromUtf8("\u0394\u03B4 vs. \u0394\u03B1"));
            }// plot as external plot for exporting in good porportion
            else if(_ui.residuals_checkbox->isChecked() && _ui.extern_checkbox->isChecked())
            {
                diff_plot->show();
                dradec_plot->show();
            }
        
            int color_cnt=1;
            int ac_cnt = 0;
            for(auto &name: other_sessions)
            {
               for(auto &srces: crfs[name])
               {
                   if(srces.second.d_ra.rows() == srces.second.ra.rows())
                   {
                       srces.second.ra.append_cols(srces.second.d_ra);
                       srces.second.dec.append_cols(srces.second.d_dec);
                   }
//                   srces.second.ra.show();
//                   srces.second.dec.show();
                                   
                // used for grid and residual arrows
                double line_width = _ui.arrow_width_spinbox->value(); 
                   
                QString legend_name = QString::fromStdString(name)+"_"+QString::fromStdString(srces.first);
                
                QColor marker_color = QColor(ivg::color_values.at(color_cnt).c_str());
                if(color_cnt == crfs.size())
                    marker_color = QColor("#ACACAC");
                
                QFormat format =  { marker_color, line_width, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssCross, legend_name, 3.0 , 10.0 };
                
                if(srces.first == "D")
                    format =  { marker_color, line_width, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssDisc, legend_name, 1.0 , 12.0 };
                else if(srces.first == "S")
                    format =  { marker_color, line_width, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssDiamond, legend_name, 3.0 , 10.0 };

//                if(ac_cnt == crfs.size()-1 && crfs.size() == 1)
//                {
//                    format =  { Qt::blue, line_width, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssDiamond, legend_name, 4.0 , 3.0 };
//                    
//                    if(srces.first == "D")
//                        format =  { Qt::green, line_width, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssDisc, legend_name, 1.0 , 6.0 };
//                    else if(srces.first == "S")
//                        format =  { Qt::red, line_width, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssCross, legend_name, 2.0 , 4.0 };
//                }
                
                
                protype pro_type = protype::mollweide;
                if(_ui.mollweide_radiobutton->isChecked())
                    pro_type = protype::mollweide;
                else if(_ui.natural_radiobutton->isChecked())
                    pro_type = protype::naturalearth;
                else if(_ui.mercator_radiobutton->isChecked())
                    pro_type = protype::mercator;
                
                
                double scale_arrows = _ui.scale_arrows_spinbox->value();
                double ref_arrow_value = _ui.ref_arrow_spinbox->value();
        
                QCPProjection projection = {pro_type, 60.0, 30.0, true, false, line_width/1.8, Qt::red, line_width/1.8, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 5.0*(1+(line_width/2.8)), 5.0*(1+(line_width/2.8)), scale_arrows, ref_arrow_value};
                plot->projection(projection, srces.second.ra, srces.second.dec, format ,"", srces.second.tooltips );
                
                // in case of plotting residuals, we also plot source differences
                if(_ui.residuals_checkbox->isChecked())
                {
                    // Generating one plot-figure containing two subplots
                    diff_plot->plot_data(srces.second.dec*ivg::rad2d,  srces.second.d_ra_cos*ivg::rad2mas*1000, format, 0 , srces.second.tooltips );
                    diff_plot->get_plot()->axisRect(0)->axis(QCPAxis::atLeft)->setLabel(QString::fromUtf8("\u0394\u03B1*cos(\u03B4) [\u00B5as]"));
                    diff_plot->get_plot()->axisRect(0)->axis(QCPAxis::atBottom)->setLabel(QString::fromUtf8("Declination [\u00B0]"));
                    diff_plot->get_plot()->axisRect(0)->axis(QCPAxis::atBottom)->setRange(-91.0,91.0);
                    diff_plot->plot_data(srces.second.dec*ivg::rad2d,  srces.second.d_dec*ivg::rad2mas*1000, format, diff_subplot, srces.second.tooltips);
                    diff_plot->get_plot()->axisRect(diff_subplot)->axis(QCPAxis::atLeft)->setLabel(QString::fromUtf8("\u0394\u03B4 [\u00B5as]"));
                    diff_plot->get_plot()->axisRect(diff_subplot)->axis(QCPAxis::atBottom)->setLabel(QString::fromUtf8("Declination [\u00B0]"));
                    diff_plot->get_plot()->axisRect(diff_subplot)->axis(QCPAxis::atBottom)->setRange(-91.0,91.0);
                    
                    // plotting ra and dec differences against each other
                    dradec_plot->plot_data( srces.second.d_ra_cos*ivg::rad2mas*1000,srces.second.d_dec*ivg::rad2mas*1000, format, 0 , srces.second.tooltips );
                    dradec_plot->get_plot()->xAxis->setLabel(QString::fromUtf8("\u0394\u03B1*cos(\u03B4) [\u00B5as]"));
                    dradec_plot->get_plot()->yAxis->setLabel(QString::fromUtf8("\u0394\u03B4 [\u00B5as]"));
                }
                
               }

            color_cnt += 1;
            
            ac_cnt++;
            }
            
            diff_plot->get_plot()->replot();
}
// ...........................................................................
void Analyzer::performStationAnalysis()
// ...........................................................................
{ 
    
//    std::vector<ivg::Trf> trf_list;
//    ivgat::Station_analysis sta_at(trf_list);
//    
//    
//     _loaded_sinex_series
    
//    sta_at.calc_station_positions( trfs_ref );
            
//    ivg::Matrix REN; 
//    ivg::Matrix epoch;
//    sta_at.get_station_position( sta, type, REN, epoch );
//
//    if( REN.rows() != 0 && REN.cols() != 0 )
//    {
//        REN *= 1e3;      
//        plot1.plot_data( epoch.transpose(), REN.transpose()(":",0), 
//                         {QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
//                          QCPGraph::lsLine, QCPScatterStyle::ssDisc, version.c_str(), 3.0, 3.0});
//        plot1.set_title( "Station positions" );
//        QCustomPlot *plt_ptr = plot1.get_plot();
//        plt_ptr->xAxis->setLabel("");
//        plt_ptr->yAxis->setLabel("Up est. [mm]");
//        plt_ptr->legend->setVisible(false);
//        QFont font = plt_ptr->xAxis->selectedLabelFont();
//        plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
//        plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) ); 
//        plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
//        plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                    
//
//        if(iter==sols.begin()) 
//            subplot1 = plot1.add_rowplot("East est. [mm]","");   
//
//        plot1.plot_data( epoch.transpose(), REN.transpose()(":",1),
//                         {QColor(ivg::color_values.at( p ).c_str()), 2.0, Qt::SolidLine, 
//                          QCPGraph::lsLine, QCPScatterStyle::ssDisc, version.c_str(), 3.0, 3.0}, subplot1);
//        plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
//        plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );  
//        plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
//        plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                    
//
//        if(iter==sols.begin())
//            subplot2 = plot1.add_rowplot("North est. [mm]","MJD");   
//
//        plot1.plot_data( epoch.transpose(), REN.transpose()(":",2),
//                         {QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
//                          QCPGraph::lsLine, QCPScatterStyle::ssDisc, version.c_str(), 3.0, 3.0}, subplot2);     
//
//        plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
//        plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );    
//        plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
//        plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                                       
//
//        plot1.get_plot()->replot();
    
}
// ...........................................................................
void Analyzer::plotSelectedEops()
// ...........................................................................
{     
    // get selected parameter, row by row, like "EOP_xpo_0" and "EOP_xpo_1"
    QModelIndexList list = _ui.param_tablewidget->selectionModel()->selectedRows();
    // also get selected eop-series; might be empty
    QModelIndexList eops = _ui.eops_listwidget->selectionModel()->selectedRows();
    
    struct grpdata{ivg::Matrix ra; ivg::Matrix dec; ivg::Matrix d_ra; ivg::Matrix d_dec; vector<string> tooltips;};
    
    // only plot if a row in param_tablewidget is selected
    if(list.size()>0)
    {       
        // for every selected parameter (=row)
        for(auto &index: list)
        {
            // sel_param: EOP_xpo_0
            string sel_param = index.data().toString().toStdString();
            ivg::Param param = _params_in_tablewidget.at(index.row());

            Plot *plot = new Plot();
            // plot as external plot for exporting in good porportion if checked
            if(!_ui.extern_checkbox->isChecked())
            {
                _ui.tabWidget->addTab(plot, QString::fromStdString(sel_param));
                _ui.tabWidget->setCurrentIndex(_ui.tabWidget->indexOf(plot));
            }
            else
                plot->show();
            
            // for every loaded session (=column) loaded in param_tablewidget
            int color_cnt = 1;

            // in case of relative plotting to first column, skip the first column (because will be zero)
            int start=1;
            
            int x_index =  _ui.xaxis_combobox->currentIndex();
            int y_index =  _ui.yaxis_combobox->currentIndex();
            if(y_index == 1)
                start=2;

            for(int i=start; i<_header_labels.size(); i++)
            {
                // selected session-unit (single file or folder, e.g. "cgs2014a_E")
                string sel_sess = _header_labels.at(i).toStdString();

                ivg::Matrix mjd, values, std;
                vector<string> tooltips;

                // get requested eop series
                ivg::Eop_series eop_series = _loaded_eops[sel_sess];
                // get raw data of requested parameter [mjd value std]
                ivg::Matrix raw_data = eop_series.get_time_series(param.get_typename(), param.get_order());
                
                // generate label depending on parameter (name, order, unit)
                stringstream label;
                label << param.get_typename();
                if(param.get_order() == 1)
                    label << "-rate ";

                label << " [" << param.get_unit() << "]" << endl;
                string label_str = label.str();

                // in case of relative plotting
                string title = "solution";
                
                // case difference
                if(y_index == 1)
                {
                    string relative_to =  _header_labels.at(1).toStdString();
                    
                    ivg::Eop_series relative_eop_series = _loaded_eops[relative_to];
                    // create title
                    title += " minus "+relative_to;
                    // compute difference of both eop-series
                    
                    int interpolate_index =  _ui.interpolateComboBox->currentIndex();
                    
                    ivg::Eop_series diff;
                                        
                    if(interpolate_index == 0){
                        // interpolate subtrahend;
                        diff = eop_series - relative_eop_series;
                        raw_data = diff.get_time_series(param.get_typename(), param.get_order());
                    } else {
                        // interpolate minuend;
                        diff = relative_eop_series - eop_series;
                        raw_data = diff.get_time_series(param.get_typename(), param.get_order());
                        raw_data.set_col(1,raw_data.get_col(1)*-1.0);
                    }

                    
                    // adjust unit from milli to micro (*1000)
                    if(param.get_typename() == "ut1" && param.get_order() == 0)
                    {
                        // raw_data from seconds to microseconds
                        raw_data.set_col(1,raw_data.get_col(1)*1e6);
                        raw_data.set_col(2,raw_data.get_col(2)*1e6);
                    }
                    else
                    {
                        raw_data.set_col(1,raw_data.get_col(1)*1000.0);
                        raw_data.set_col(2,raw_data.get_col(2)*1000.0);
                    }

                    replace_string_in_place(label_str, "m", "micro" );  
                }

                // the final data to be plotted
                mjd = raw_data.get_col(0);
                values = raw_data.get_col(1);   
                std = raw_data.get_col(2);

                // generate tooltips
                vector<double> x_values;
                string x_label="";
                for(int mjd_i=0; mjd_i<mjd.rows(); mjd_i++ )
                {
                    ivg::sessinfo info;
                    if(eop_series.get_origin_db(mjd(mjd_i)) == "UNKNOWN")
                        info = _masterfile.get_session_info(ivg::Date(mjd(mjd_i)));
                    else
                        info = _masterfile.get_session_info(eop_series.get_origin_db(mjd(mjd_i)));
                    
                    // if still unknown, try to check for intensives
                    if(info.tooltip == "UNKNOWN")
                        info = _masterfile.get_session_info(ivg::Date(mjd(mjd_i)),0.0,ivg::mastertype::intensive);
                    
                    //in case of num of stations ornetwork volume we need to store the information
                    if(x_index == 1)
                    {
                        x_values.push_back(info.trf.get_number_stations());
                        x_label = "number of stations [#]";
                    }
                    else if(x_index == 2)
                    {
                        x_values.push_back(info.volume);
                        x_label = "network volume [Mm^3]";
                    }
                    
                    tooltips.push_back(info.tooltip);
                }

                // generate plot
                QFormat format1 = { QColor(ivg::color_values.at(color_cnt).c_str()), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, sel_sess.c_str(), 1.0 , 6.0  };
                
                // in case of time series plotting
                if(x_index == 0)
                {
                    string time_format = _ui.timeformat_lineedit->text().toStdString();
                    if(y_index == 0 || y_index == 1)
                        plot->plot_mjd_series_std(mjd,values,std, format1, 0, time_format, tooltips);
                    else if(y_index == 2)
                        plot->plot_mjd_series(mjd,std, format1, 0, time_format, tooltips);
                }
                // in case of num of stations or network volume
                else if(x_index == 1 || x_index == 2)
                {
                    if(y_index == 2)
                    {
                        values = std;
                        label_str = "standard deviation";
                    }
                    
                    plot->plot_data(x_values,values,format1,0,tooltips);
                }
                
                plot->get_plot()->xAxis->setLabel(QString::fromStdString(x_label));                
                plot->get_plot()->yAxis->setLabel(QString::fromStdString(label_str));
                plot->set_title(title);

                color_cnt += 1;
            }

            // update view after every new insert plot
            plot->switchErrorbars();
            plot->get_plot()->rescaleAxes();
            plot->get_plot()->replot();      
        }
    }
    
}
// ...........................................................................
void Analyzer::plotSelectedStations()
// ...........................................................................
{            
    // get selected parameter, row by row, like "WETTZELL_stax_0"
    QModelIndexList list = _ui.param_tablewidget->selectionModel()->selectedRows();
    
    // only in case of atleast one selected row (=parameter)
    if(list.size() > 0)
    {
        // load each column as a vector of trfs
        vector< vector<ivg::Trf> > trfs_vec;
        for(int i=1; i<_header_labels.size(); i++)
        {
            // version v (e.g. itrf2014_E)
            string v = _header_labels.at(i).toStdString(); // itrf2014_E

            vector<ivg::Trf> trfs;     
            for(auto &sinex: _loaded_sinex_series[v.substr(0,v.length()-2)])
            {
                if(v.substr(v.length()-2,2) == "_A")
                    trfs.push_back(sinex.get_trf(ivg::reftype::apriori));
                else if(v.substr(v.length()-2,2) == "_E")
                    trfs.push_back(sinex.get_trf(ivg::reftype::estimate));
            }
            trfs_vec.push_back(trfs);
        }
        
        // use first column as vector of reference trfs
        ivgat::Station_analysis sta_at(trfs_vec.at(0));
        
        if(_ui.sta_timeseries_checkbox->isChecked())
        {
            // check stations
            double max_est = _ui.maxest_spinbox->value() / 1000.0;

            for(int i=1; i<trfs_vec.size(); i++)
            {
                //loop over all other trfs left in trf_vec
                Plot *plot = new Plot();
                stringstream tabtitle;
                tabtitle << _header_labels.at(0).toStdString() << " vs " << _header_labels.at(i+1).toStdString();
                _ui.tabWidget->addTab(plot, QString::fromStdString(tabtitle.str()));
                int tab_index = _ui.tabWidget->indexOf(plot);
                _ui.tabWidget->setCurrentIndex(tab_index);
                plot->get_plot()->xAxis->setLabel("");
                plot->get_plot()->yAxis->setLabel("Up est. [mm]");
                int subplot1 = plot->add_rowplot("East est. [mm]",""); 
                int subplot2 = plot->add_rowplot("North est. [mm]","");   

                plot->set_title( tabtitle.str() );
                sta_at.check_stations( trfs_vec.at(i), max_est );   
                sta_at.calc_station_positions( trfs_vec.at(i), max_est );

                // loop over all stations (selected row => list)
                int color_cnt=0;
                for(auto &index: list)
                {
                    string station = _params_in_tablewidget.at(index.row()).get_name();

                    ivg::Matrix REN; 
                    ivg::Matrix REN_std; 
                    ivg::Matrix epoch;
                    sta_at.get_station_position(station, "REN", REN, REN_std, epoch );

                    vector<string> tooltips;
                    // generate tooltip for each session
                    for(int mjd_i=0; mjd_i<epoch.cols(); mjd_i++ )
                        tooltips.push_back(_masterfile.get_session_info(ivg::Date(epoch(mjd_i))).tooltip);

                    if( REN.rows() != 0 && REN.cols() != 0 )
                    {
                        REN *= 1e3;     
                        QFormat format = {QColor(ivg::color_values.at( color_cnt ).c_str() ), 1.5, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, QString::fromStdString(station), 1.5, 3.0};
                        plot->plot_mjd_series( epoch.transpose(), REN.transpose()(":",0), format, 0, "dd-MMM-yy",tooltips);                
                        plot->plot_mjd_series( epoch.transpose(), REN.transpose()(":",1), format, 1, "dd-MMM-yy",tooltips);
                        plot->get_plot()->legend->itemWithPlottable(plot->get_plot()->graph())->setVisible(false);                
                        plot->plot_mjd_series( epoch.transpose(), REN.transpose()(":",2), format, 2, "dd-MMM-yy",tooltips);
                        plot->get_plot()->legend->itemWithPlottable(plot->get_plot()->graph())->setVisible(false);
                        // update view after every new insert plot
                        plot->get_plot()->rescaleAxes();
                        plot->get_plot()->replot();         
                    }

                    color_cnt++;
                }
            }
        }
        
        if(_ui.sta_baselines_checkbox->isChecked())
        {
            int min_bl_no = _ui.blreps_spinbox->value();
            
            // check if minimum amount of baselines reached
            if( sta_at.get_session_number() > min_bl_no )
            {
                std::map< std::string, ivg::Matrix > bl_reps;
                
                // calculate baseline lengths
                sta_at.calc_baseline_length();
                //sta_at.show_baselines();

                // calculate baseline lengths repeatabilities
                ivg::Matrix bl_rep;
                std::vector<std::string> bl_name;
                sta_at.calc_bl_rep( bl_name, bl_rep, min_bl_no );
                // sta_at.show_bl_repeatabilities();

                // remove selected baselines
//                sta_at.remove_baselines( rm_bls );

                // fit trend to baseline lengths repeatabilities
                ivg::Matrix fit_mat;
                sta_at.fit_bl_rep( fit_mat );

                bl_reps = sta_at.get_bl_rep();
                ivg::Matrix X( bl_reps.size(), 1, 0.0 );
                ivg::Matrix Y( bl_reps.size(), 1, 0.0 );
                
                int counter = 0;
                std::vector<std::string> tooltips;
                for (auto &it: bl_reps)
                {
                   tooltips.push_back( it.first );
                   X(counter,0) = it.second(0,0);
                   Y(counter,0) = it.second(0,1);
                   counter++;
                }
                
                int color_cnt = 0;
                for(int i=1; i<trfs_vec.size(); i++)
                {
                    //loop over all other trfs left in trf_vec
                    Plot *plot = new Plot();
                    stringstream tabtitle;
                    tabtitle << _header_labels.at(0).toStdString() << " vs " << _header_labels.at(i+1).toStdString();
                    _ui.tabWidget->addTab(plot, QString::fromStdString(tabtitle.str()));
                    int tab_index = _ui.tabWidget->indexOf(plot);
                    _ui.tabWidget->setCurrentIndex(tab_index);
                    plot->get_plot()->xAxis->setLabel("baseline length [km]");
                    plot->get_plot()->yAxis->setLabel("RMS [mm]");
    //                plot->set_title( "baseline length repeatabilities"+tabtitle.str() );

                    // PLOT baseline repeatabilities 
                    QFormat format_reps = {QColor(ivg::color_values.at( color_cnt ).c_str() ), 1.5, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, QString::fromStdString("blub"), 1.5, 3.0};
                    plot->plot_data( X*1e-3, Y*1e3,format_reps, 0, tooltips );
                    QFormat format_fit = {QColor(ivg::color_values.at( color_cnt ).c_str() ), 1.5, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "", 0.0, 0.0};
                    plot->plot_data( X*1e-3, fit_mat(":",1)*1e3,format_fit);
                                       
                    color_cnt++;
                }
            }
        }

        // if the user wants helmert parameter time series
        if(_ui.sta_helmert_checkbox->isChecked() )
        {
            // get the parameter the user wants to estimate
            vector<ivgat::t_param> trans_params;
            vector<string> tp_labels;
            vector<double> tp_fak;
            if(_ui.tx_sta_checkbox->isChecked())
            {
                trans_params.push_back(ivgat::t_param::tx);
                tp_labels.push_back("Tx [mm]");
                tp_fak.push_back(1000.0);
            }
            if(_ui.ty_sta_checkbox->isChecked())
            {
                trans_params.push_back(ivgat::t_param::ty);
                tp_labels.push_back("Ty [mm]");
                tp_fak.push_back(1000.0);
            }
            if(_ui.tz_sta_checkbox->isChecked())
            {
                trans_params.push_back(ivgat::t_param::tz);
                tp_labels.push_back("Tz [mm]");
                tp_fak.push_back(1000.0);
            }            
            if(_ui.rx_sta_checkbox->isChecked())
            {
                trans_params.push_back(ivgat::t_param::rx);
                tp_labels.push_back("Rx [microas]");
                tp_fak.push_back(ivg::rad2mas*1000);
            }
            if(_ui.ry_sta_checkbox->isChecked())
            {
                trans_params.push_back(ivgat::t_param::ry);
                tp_labels.push_back("Ry [microas]");
                tp_fak.push_back(ivg::rad2mas*1000);
            }
            if(_ui.rz_sta_checkbox->isChecked())
            {
                trans_params.push_back(ivgat::t_param::rz);
                tp_labels.push_back("Rz [microas]");
                tp_fak.push_back(ivg::rad2mas*1000);
            }
            if(_ui.s_sta_checkbox->isChecked())
            {
                trans_params.push_back(ivgat::t_param::s);
                tp_labels.push_back("S [-]");
                tp_fak.push_back(1.0);
            }

            if(trans_params.size() == 0)
                _ui.info_textbrowser->setText(QString("No transformation parameters selected. Checkbox them!"));
            else
            {
                // calcluate the parameter series once and save them
                vector<ivg::Matrix> raw_data;
                for(auto &tmp_trf: trfs_vec)
                    raw_data.push_back(sta_at.calc_helmert_params(tmp_trf, trans_params));

                for(int i=0; i<trans_params.size();i++)
                {
                    int color_cnt=0;
                    for(int j=1; j<trfs_vec.size(); j++)
                    {
                        Plot *plot;
                        int subplot=0;
                        if(j == 1)
                        {
                            if( i == 0 || i == 3 || i == 6 )
                            {
                                //loop over all other trfs left in trf_vec
                                plot = new Plot();
                                _ui.tabWidget->addTab(plot, QString::fromStdString("HelmertParams"));
                                int tab_index = _ui.tabWidget->indexOf(plot);
                                _ui.tabWidget->setCurrentIndex(tab_index);
                                plot->get_plot()->xAxis->setLabel(QString::fromStdString(""));
                                plot->get_plot()->yAxis->setLabel(QString::fromStdString(tp_labels.at(i)));
                            }
                            else
                                subplot = plot->add_rowplot(tp_labels.at(i)," ");
                        }
                        else
                            subplot = i % 3;

                        cerr << "Param: " << tp_labels.at(i) << " / Subplot: " << subplot << " / i: " << i << " / j: " << j << endl;

                        raw_data.at(j).show();

                        ivg::Matrix mjd = raw_data.at(j).get_col(0);
                        ivg::Matrix param_values = raw_data.at(j).get_col(i+1) * tp_fak.at(i);   

                        // generate tooltip for each session
                        vector<string> tooltips;
                        for(int mjd_i=0; mjd_i<mjd.rows(); mjd_i++ )
                            tooltips.push_back(_masterfile.get_session_info(ivg::Date(mjd(mjd_i))).tooltip);

                        mjd.show();
                        param_values.show();

                        QFormat format = {QColor(ivg::color_values.at( color_cnt ).c_str() ), 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, QString::fromStdString(tp_labels.at(i)), 3.0, 3.0};
                        plot->plot_mjd_series( mjd, param_values, format, subplot, "dd-MMM-yy",tooltips);
                        // update view after every new insert plot
                        plot->get_plot()->rescaleAxes();
                        plot->get_plot()->replot();    

                        color_cnt++;
                    }
                }
            }
        }
    }
}
// ...........................................................................
void Analyzer::plotSelectedSources()
// ...........................................................................
{ 
    QModelIndexList list = _ui.param_tablewidget->selectionModel()->selectedRows();
    
    // only plot if something is selected
    if(list.size()>0)
    {       
        // for every selected source
        for(auto &index: list)
        {
                
            string sel_src = index.data().toString().toStdString();
            
            QWidget *newWidget1 = new QWidget(); 
            Plot *plot = new Plot(newWidget1);
            _ui.tabWidget->addTab(plot, QString());
            int tab_index = _ui.tabWidget->indexOf(plot);
            _ui.tabWidget->setTabText(tab_index, QString::fromStdString(sel_src));
            _ui.tabWidget->setCurrentIndex(tab_index);
            
            // for every session loaded in param_tablewidget
            int color_cnt = 3;
            double mean_ra,mean_dec;
            for(int i=1; i<_header_labels.size(); i++)
            {
                string sel_sess = _header_labels.at(i).toStdString();
                ivg::Param_list *param_list = _loaded_sessions[sel_sess].get_param_list_ptr();
                
                vector<int> ra_idx = param_list->get_idx(ivg::paramtype::ra,sel_src);
                vector<int> dec_idx = param_list->get_idx(ivg::paramtype::dec,sel_src);
                                
                vector<double> ra_vec,dec_vec,mjd_vec;
                if(ra_idx.size() == dec_idx.size())
                {
                    for(int idx=0; idx<ra_idx.size(); idx++)
                    {
                        ra_vec.push_back(param_list->get_param(ra_idx.at(idx))->get_estimate() * ivg::rad2mas * 1000);
                        dec_vec.push_back(param_list->get_param(dec_idx.at(idx))->get_estimate() * ivg::rad2mas * 1000);
                        mjd_vec.push_back(param_list->get_param(dec_idx.at(idx))->get_epoch().get_double_mjd());
                    }
                    
                    ivg::Matrix ra(ra_vec);
                    ivg::Matrix dec(dec_vec);
                    ivg::Matrix mjd(mjd_vec);
                    
                    if(i == 1)
                    {
                        mean_ra = ra.meanD();
                        mean_dec = dec.meanD();
                    }
                    
                    ra = ra - mean_ra;
                    dec = dec - mean_dec;
                    
                    
                    cerr << "HIER "<<endl;
                    ra.show();
                    dec.show();
                    mjd.show();
                    
                    int subplot = 0;
                    QFormat format1 = { QColor(ivg::color_values.at(color_cnt).c_str()), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, sel_sess.c_str(), 1.0 , 4.0  };
                    plot->plot_mjd_series(mjd,ra, format1, subplot, "dd-MMM-yy");
                    plot->get_plot()->xAxis->setLabel(QString(""));
                    
                    plot->get_plot()->yAxis->setLabel(QString("Right Ascension [microas]"));
//                    
                    if(i == 1)
                        subplot = plot->add_rowplot("Declination [microas]","Time");
                    else
                        subplot = 1;
                    QFormat format2 = { QColor(ivg::color_values.at(color_cnt).c_str()), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, sel_sess.c_str(), 1.0 , 4.0  };
                    plot->plot_mjd_series(mjd,dec, format2, subplot, "dd-MMM-yy");
//                    
                    plot->set_title(sel_src);
                                        
                }
                
                color_cnt += 2;
            }
              
            plot->get_plot()->rescaleAxes();
            plot->get_plot()->replot();
        }
    }
}
// ...........................................................................
void Analyzer::transformSelectedSessions()
// ...........................................................................
{
    // get reference
    string ref_name = "";
    if(_latest_reference_row != -1)
    {
        ref_name = _ui.session_tablewidget->item(_latest_reference_row,0)->text().toStdString();
	std::ofstream transf_file("transf_params.out",std::ofstream::out);
        // get selected sessions
        vector<string> other_sessions;
        for(int row=0; row<_ui.session_tablewidget->rowCount(); row++)
        {
            if(_ui.session_tablewidget->item(row,0)->checkState() == Qt::CheckState::Checked)
                other_sessions.push_back(_ui.session_tablewidget->item(row,0)->text().toStdString());
        }

        vector<ivgat::t_param> trans_params;
        if(_ui.tx_checkbox->isChecked())
            trans_params.push_back(ivgat::t_param::tx);
        if(_ui.ty_checkbox->isChecked())
            trans_params.push_back(ivgat::t_param::ty);
        if(_ui.tz_checkbox->isChecked())
            trans_params.push_back(ivgat::t_param::tz);            
        if(_ui.rx_checkbox->isChecked())
            trans_params.push_back(ivgat::t_param::rx);            
        if(_ui.ry_checkbox->isChecked())
            trans_params.push_back(ivgat::t_param::ry);            
        if(_ui.rz_checkbox->isChecked())
            trans_params.push_back(ivgat::t_param::rz);
        if(_ui.s_checkbox->isChecked())
            trans_params.push_back(ivgat::t_param::s);

        if(trans_params.size() == 0)
            _ui.info_textbrowser->setText(QString("No transformation parameters selected. Checkbox them!"));
        else if( other_sessions.size() > 0 && !ref_name.empty() && trans_params.size()>0)
        {
            stringstream trans_info;
            for(auto &other_name: other_sessions)
            {
	      vector<ivg::Sinex> othersnxs;
	      othersnxs=_loaded_sinex_series[other_name];
	      if (othersnxs.size()==0)
		othersnxs.push_back(_loaded_sinex[other_name]);
	     
	      for (auto &other_snx: othersnxs )
		{
                if(_ui.trf_radiobutton->isChecked())
                {
		  
                    ivg::Trf ref_trf = _loaded_sinex[ref_name].get_trf(ivg::reftype::estimate);
		    if (_psd_file[ref_name].length()>0)
		      ivg::parser::psd_coefficients(&ref_trf,_psd_file[ref_name]);
                    //ivg::Trf other_trf = _loaded_sinex[other_name].get_trf(ivg::reftype::estimate);
		    ivg::Trf other_trf = other_snx.get_trf(ivg::reftype::estimate);

                    ivg::Matrix indexes = ref_trf.get_corresponding_stations(other_trf);
                    ivg::Matrix indexes_orig = indexes;
                                  
                    ivg::Matrix est_system1,est_system2;
                    ivg::Matrix all_system1,all_system2;
                    
                    for(int i=0; i<indexes.rows(); i++)
                    {
		      // all_system1.append_cols(ref_trf.get_station(indexes(i,0))->get_xyz(ivg::Date(2013,123.5)));
                      //  all_system2.append_cols(other_trf.get_station(indexes(i,1))->get_xyz(ivg::Date(2013,123.5)));
		      
		      ivg::Date tmpep=other_trf.get_station(indexes(i,1))->get_refepoch();
		      all_system1.append_cols(ref_trf.get_station(indexes(i,0))->calc_xyz(tmpep,{"PSD"}));
		      all_system2.append_cols(other_trf.get_station(indexes(i,1))->calc_xyz(tmpep));
                    }
                    
                    // now only use stations selected in table for transformation (default = all)
                    for(int row=0; row<_ui.select_tablewidget->rowCount(); row++)
                    {
                        if(_ui.select_tablewidget->item(row,0)->checkState() == Qt::CheckState::Unchecked)
                        {
                            ivg::Matrix ref_column = indexes(":",0);
                            indexes.rem_r(ref_column.find_idx(row));
                        }
                    }

                    // now only set systems used for computation of transformatio parameters
                    for(int i=0; i<indexes.rows(); i++)
                    {
		       ivg::Date tmpep=other_trf.get_station(indexes(i,1))->get_refepoch();
		       //est_system1.append_cols(ref_trf.get_station(indexes(i,0))->get_xyz(ivg::Date(2000,123.5)));
		       //est_system2.append_cols(other_trf.get_station(indexes(i,1))->get_xyz(ivg::Date(2000,123.5)));
		        est_system1.append_cols(ref_trf.get_station(indexes(i,0))->calc_xyz(tmpep,{"PSD"}));
                        est_system2.append_cols(other_trf.get_station(indexes(i,1))->calc_xyz(tmpep));
                    }
		    
                    est_system1 = est_system1.transpose();
                    est_system2 = est_system2.transpose(); 
		    
                    all_system1 = all_system1.transpose();
                    all_system2 = all_system2.transpose();
		   
                    ivgat::Transformation trans(ivgat::t_type::cart_3D);
                    
                    // set the two systems for the parameter estimation
                    trans.set_systems(est_system1,est_system2);

                    // define return matrices (call by reference) and the parameter which should be estimated
                    ivg::Matrix params,accu,corr;
                    // get results and a nice info block for output
                    string info_block = trans.estimate_parameter(trans_params, params, accu, corr);
                    
                    // now use ALL stations calculating residuals
                    trans.set_systems(all_system1,all_system2);
                    
                    // get the residuals after transformation on each other between both systems
                    ivg::Matrix resids = trans.get_residuals(trans_params,params);
                    cerr << "TRF residuals between system1 / system2 based on cartesian coordinates [m]" << endl;
                    resids.show(4);
		    
                    indexes_orig.append_cols(resids);
                    _trf_residuals[ref_name+"_"+other_name] = indexes_orig;
		    string curpath=other_snx.get_path();
		    string sess = curpath.substr(curpath.find_last_of("/")+1,curpath.find_first_of("_")-curpath.find_last_of("/")-1);
                    trans_info << "------------------------------------------------------------------" << endl;
                    trans_info << ref_name << " vs " << other_name << " (" << sess << ")" << endl;
                    trans_info << info_block;

                    _ui.info_textbrowser->setText(QString::fromStdString(trans_info.str()));
		    double fak=1,add=0;
		    int cnt=0;
		    transf_file << sess << " " << setw(14) << fixed << setprecision(4) << other_trf.get_station(indexes(0,1))->get_refepoch().get_double_mjd() << " " ;
		    for (auto &tp: trans_params)
		      {
			 if( tp <= 2 )
			   {
			     fak = 1000.0;
			     add=0;
			   }
			 else if( tp > 2 && tp <= 5)
			   {
			     fak = ivg::rad2mas*1000;
			     add=0;
			   }
			 else if( tp == 6 )
			   {
			     fak = 1.0e9;
			     add=-1;
			   }
			 transf_file << right << setw(14) << fixed << setprecision(6) << (params(cnt)+add)*fak << " ";
			 cnt++;
		      }
		    transf_file << endl;
                }
                else if(_ui.crf_radiobutton->isChecked())
                {
                    ivg::Crf ref_crf = _loaded_sinex[ref_name].get_crf(ivg::reftype::estimate);
                    //ivg::Crf other_crf = _loaded_sinex[other_name].get_crf(ivg::reftype::estimate);
                    ivg::Crf other_crf = other_snx.get_crf(ivg::reftype::estimate);
                    ivg::Matrix indexes = ref_crf.get_corresponding_sources(other_crf);
                    
                    ivg::Matrix est_system1,est_system2;
                    ivg::Matrix all_system1,all_system2;

                    for(int i=0; i<indexes.rows(); i++)
                    {
                        if( ( _ui.def_only_checkbox->isChecked() && ref_crf.get_source(indexes(i,0))->is_defining() ) ||
                            ( _ui.except_sh_checkbox->isChecked() && !ref_crf.get_source(indexes(i,0))->is_special_handling()))
                        {
                            est_system1.append_cols(ref_crf.get_source(indexes(i,0))->get_sph());
                            est_system2.append_cols(other_crf.get_source(indexes(i,1))->get_sph());
                        }

                        all_system1.append_cols(ref_crf.get_source(indexes(i,0))->get_sph());
                        all_system2.append_cols(other_crf.get_source(indexes(i,1))->get_sph());
                    }

                    if(!_ui.def_only_checkbox->isChecked() &&  !_ui.except_sh_checkbox->isChecked())
                    {
                        est_system1 = all_system1;
                        est_system2 = all_system2;
                    }

                    est_system1 = est_system1.transpose();
                    est_system2 = est_system2.transpose();

                    all_system1 = all_system1.transpose();
                    all_system2 = all_system2.transpose();
                    
                    ivg::Matrix raw_resids_sph = all_system2 - all_system1;
                    
                    ivgat::Transformation trans(ivgat::t_type::sph_3D);

                    // set the two systems for the parameter estimation
                    trans.set_systems(est_system1,est_system2);

                    // define return matrices (call by reference) and the parameter which should be estimated
                    ivg::Matrix params,accu,corr;
                    // get results and a nice info block for output
                    string info_block = trans.estimate_parameter(trans_params, params, accu, corr);

                    // get the residuals after transformation on each other between both systems
                    trans.set_systems(all_system1,all_system2);
                    ivg::Matrix all_trans_system1_cart = trans.transform_sph_to_cart(all_system2) + trans.get_residuals(trans_params,params);
                    ivg::Matrix all_trans_system1_sph = trans.transform_cart_to_sph(all_trans_system1_cart);

                    ivg::Matrix trans_resids_sph = all_trans_system1_sph - all_system1;
                    indexes.append_cols(trans_resids_sph); // transformed differences as columns 3 & 4
                    indexes.append_cols(raw_resids_sph); // raw differences as columns 5 & 6
                    
                    for(int i=0; i<indexes.rows(); i++)
                    {
                        if(indexes(i,2) < -6.28)
                            indexes(i,2) += 2*M_PI;
                    }
                    
                    indexes.show();
                    _crf_residuals[ref_name+"_"+other_name] = indexes;

                    trans_info << "------------------------------------------------------------------" << endl;
                    trans_info << ref_name << " vs " << other_name << endl;
                    trans_info << info_block;
		    
		      
                }
		}
            }

            _ui.info_textbrowser->setText(QString::fromStdString(trans_info.str()));
        }
        else
            _ui.info_textbrowser->setText(QString("Only reference session selected. No other sessions. Checkbox session names for transformation!"));
	transf_file.close();
    }
    else
        _ui.info_textbrowser->setText(QString("No reference session selected. Double click on a session name!"));
    
}
// ...........................................................................
void Analyzer::loadEopSeries()
// ...........................................................................
{  
    // load a different eop-series, e.g. finals or c04 in order to be able to compare them
    EopSeries *load_eops = new EopSeries();

    // in case of defined eop-series from config file, set them as default
    for(int i=0; i<(*_setup)["eop_files"].getLength(); i++)
        load_eops->set_eop_file((*_setup)["eop_files"][i][0],(*_setup)["eop_files"][i][1],(*_setup)["eop_files"][i][2]);

    // in case of selected session in session_tablewidget, use the epochs for eop-series initialization
    QModelIndexList list =_ui.session_tablewidget->selectionModel()->selectedIndexes();
    
    // get information about loading-range
    int sY,sM,sD,eY,eM,eD;
    _ui.start_dateedit->date().getDate(&sY,&sM,&sD);
    _ui.end_dateedit->date().getDate(&eY,&eM,&eD);

    // adjust year to correct century
    if(sY >= 2070)
        sY -= 100;

    if(eY >= 2070)
        eY -= 100;

    ivg::Date start_epoch(sY,sM,sD);
    ivg::Date end_epoch(eY,eM,eD);

    load_eops->set_start_end_epochs(start_epoch.add_days(-10.0), end_epoch.add_days(10.0));

    // only if the user pushed OK, apply changes
    if (load_eops->exec() == QDialog::Accepted) 
    {
        string eops_name;
        ivg::Eop_series eops = load_eops->get_eop_series(eops_name);
        _loaded_eops[eops_name] = eops;

        new QListWidgetItem(QString::fromStdString(eops_name), _ui.eops_listwidget);
    }
    
}
// ...........................................................................
void Analyzer::highlightCurrentPlot()
// ...........................................................................
{
    string identify;
    bool remove = _ui.highlight_remove_checkbox->isChecked();
    
    QModelIndexList grp_list =_ui.highlight_listwidget->selectionModel()->selectedIndexes();
    QList<QListWidgetItem*> site_list = _ui.highlight_site_listwidget->selectedItems();
    
    if(grp_list.size()>0 || site_list.size()>0)
    {       
        // get current plot widget
        Plot *current_plot = (Plot*) _ui.tabWidget->currentWidget();
        if(current_plot != NULL)
        {
            // several rows could be selected in the listwidget
            stringstream ss;
            for(int i=0; i< grp_list.size(); i++)
                ss << "{9}" << grp_list.at(i).data().toString().toStdString() << "\n";
            
            if(_ui.or_site_radiobutton->isChecked())
            {
                for(int i=0; i< site_list.size(); i++)
                    ss << "{6}" << site_list.at(i)->toolTip().toStdString() << "\n";
            }
            else if(_ui.and_site_radiobutton->isChecked() && site_list.size()>0)
            {
                ss << "{6}";
                for(int i=0; i< site_list.size(); i++)
                    ss << site_list.at(i)->toolTip().toStdString() << "$";
            }
            
            (*current_plot).highlight(ss.str(),{ QColor(ivg::color_values.at(_highlight_color_cnt).c_str()), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssCircle, QString::fromStdString(ss.str()), 2.0, 7.0 }, remove);
        }
    }
    
    _highlight_color_cnt++;
}
// ...........................................................................
void Analyzer::closeTabWidgetTab(int index){
    _ui.tabWidget->removeTab(index);
}
// ...........................................................................
void Analyzer::analysisTabChanged(int index){
    objectTypeChanged();
    
    // highlighting only supported for time series of stations and eops
    if(index == 1 || index == 2)
        _ui.highlight_groupbox->setEnabled(true);
    else
        _ui.highlight_groupbox->setEnabled(false);
}
// ...........................................................................
void Analyzer::objectTypeChanged()
// ...........................................................................
{
    // adjusts the appearence of the GUI 
    if(_ui.crf_radiobutton->isChecked() )
    {
        _ui.rx_checkbox->setChecked(true);
        _ui.ry_checkbox->setChecked(true);
        _ui.rz_checkbox->setChecked(true);
        _ui.tx_checkbox->setChecked(false);
        _ui.ty_checkbox->setChecked(false);
        _ui.tz_checkbox->setChecked(false);
        _ui.s_checkbox->setChecked(false);
        _ui.def_only_checkbox->setChecked(true);
        _ui.scale_arrows_spinbox->setValue(1.0);
        _ui.ref_arrow_spinbox->setValue(1000.0);
        _ui.scale_arrows_spinbox->setSingleStep(5.0);
        _ui.ref_arrow_spinbox->setSingleStep(100.0);
        _ui.ref_arrow_label->setText(QString::fromUtf8("[\u00B5as]"));
    }
    else if(_ui.trf_radiobutton->isChecked())
    {
        _ui.rx_checkbox->setChecked(true);
        _ui.ry_checkbox->setChecked(true);
        _ui.rz_checkbox->setChecked(true);
        _ui.tx_checkbox->setChecked(true);
        _ui.ty_checkbox->setChecked(true);
        _ui.tz_checkbox->setChecked(true);
        _ui.s_checkbox->setChecked(false);
        _ui.def_only_checkbox->setChecked(false);
        _ui.scale_arrows_spinbox->setValue(2.0);
        _ui.ref_arrow_spinbox->setValue(1.0);
        _ui.ref_arrow_spinbox->setSingleStep(1.0);
        _ui.scale_arrows_spinbox->setSingleStep(1.0);
        _ui.ref_arrow_label->setText(QString::fromUtf8("[cm]"));
    }
    
    
    for(int row=0; row<_ui.session_tablewidget->rowCount(); row++){
        
        bool button=true;
        if(_ui.analysis_tabwidget->currentIndex()== 0)
            button=false;
        
        if(_ui.analysis_tabwidget->currentIndex()== 2 && _ui.session_tablewidget->item(row,3)->backgroundColor()==Qt::red)
            button=false;
        
        if(_ui.analysis_tabwidget->currentIndex()== 1 && _ui.session_tablewidget->item(row,1)->backgroundColor()==Qt::red)
            button=false;
        
        if(_ui.analysis_tabwidget->currentIndex()== 3 && _ui.session_tablewidget->item(row,2)->backgroundColor()==Qt::red)
            button=false;
        
        for(int col=4; col<= 5; col++){
            QWidget *widget = _ui.session_tablewidget->cellWidget(row,col);
            if(widget != NULL)
                _ui.session_tablewidget->cellWidget(row,col)->setVisible(button);
        }
    }
    
    //clear whole table widget
    _ui.param_tablewidget->clear();
    _ui.param_tablewidget->setRowCount(0);
    _ui.param_tablewidget->setColumnCount(1);
    _ui.param_tablewidget->setHorizontalHeaderLabels(QStringList("Parameter"));
    _ui.param_tablewidget->resizeColumnsToContents();
    
    _params_in_tablewidget.clear();
    
    _header_labels.clear();
    _header_labels << "Parameter";
}
// ...........................................................................
void Analyzer::clearLoadedSessions()
// ...........................................................................
{
    //clear whole table widget
//    clearParameterTableWidget();
    objectTypeChanged();
    
//    _ui.refsession_listwidget->clear();
//    _ui.sessions_listwidget->clear();
    _ui.info_textbrowser->clear();
    _ui.session_tablewidget->clear();
    _ui.session_tablewidget->setRowCount(0);
    _ui.session_tablewidget->setColumnCount(6);
    QList<QString> labels = {"Name","TRF","CRF","EOP","Estimate","Apriori","SNX"};
    _ui.session_tablewidget->setHorizontalHeaderLabels(QStringList(labels));
    _ui.session_tablewidget->resizeColumnsToContents();
    _loaded_sinex.clear();
    _loaded_sinex_series.clear();
    _loaded_eops.clear();
    _trf_residuals.clear();
    _crf_residuals.clear();
}
// ...........................................................................
void Analyzer::clearHighlightedItems()
// ...........................................................................
{
 
    QList<QListWidgetItem*> site_list = _ui.highlight_site_listwidget->selectedItems();
    for(auto *item: site_list)
        item->setSelected(false);
    
    site_list = _ui.highlight_listwidget->selectedItems();
    for(auto *item: site_list)
        item->setSelected(false);
    
}
