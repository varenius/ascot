#include "statistics.h"
#include "station.h"
#include "vgosdb.h"

// ...........................................................................
Statistics::Statistics(QWidget *parent) : QMainWindow(parent), ui(new Ui::Statistics)
// ...........................................................................
{   
#if DEBUG_VLBI >=2
   cerr << "+++ Statistics::Statistics(QWidget *)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ui->setupUi(this);
    ui->data_tableview->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->plot_tabwidget->clear();
    
    _bl_selected = false;
    
    // define default configuration for the projection
    _trf_projection = {protype::naturalearth, 45.0, 45.0, false, true, 1.0, Qt::black, 1.0, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 8.0, 8.0, 20.0};
    _crf_projection = {protype::mollweide, 60.0, 30.0, true, false, 1.0, Qt::black, 1.0, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 8.0, 8.0, 20.0};
    _sky_projection = {protype::polarstereo, 90.0, 30.0, true, false, 1.0, Qt::black, 1.0, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 8.0, 8.0, 20.0};
    
    _init_axis_layout();
    
    for(auto const &it : _axis_layouts )
    {
        ui->xaxis_combobox->addItem(QString::fromStdString(it.name));
        ui->yaxis_combobox->addItem(QString::fromStdString(it.name));
    }
    
    //Set combobox defaults -> hh:mm:ss vs. residuals
    ui->xaxis_combobox->setCurrentIndex(0);
    ui->yaxis_combobox->setCurrentIndex(2);
    
    connect(ui->plot_buttonbox, SIGNAL(accepted()), this, SLOT(createFigure()));
    connect(ui->cbreak_pushbutton, SIGNAL(clicked()), this, SLOT(correctCBreak()));
    
    connect(ui->residDown_pushbotton, SIGNAL(clicked()), this, SLOT(shiftResidDown() ));
    connect(ui->residUp_pushbotton, SIGNAL(clicked()), this, SLOT(shiftResidUp() ));
    
    connect(ui->save_pushbutton, SIGNAL(clicked()), this, SLOT(saveResid() ));
    
    // initialize skyplot without plotting a point
    ui->sky_widget->projection(_sky_projection, ivg::Matrix(1,1,0.0)*(M_PI/180.0), ivg::Matrix(1,1,30.0)*(M_PI/180.0), 
                        { Qt::red, 5.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssNone, "NONE", 4.0, 6.0 }, "Skyplot");
                        
    // connect moving of mouse and clicking of baseline to some actions
    connect(ui->trf_widget->get_plot(), SIGNAL(plottableClick(QCPAbstractPlottable*, QMouseEvent*)), this, SLOT(showBaselineInfo(QCPAbstractPlottable*, QMouseEvent*))); 
    // connect signal from Plot-Class from TRF-Plot with plotting Sky-Plot of selected site
    // OUT OF ORDER DUE, BUGGY
    //connect(ui->trf_widget, SIGNAL(tooltipShowed(QString&)), this, SLOT(plotSiteSkyplot(QString&))); 
        
#if DEBUG_VLBI >=2
   cerr << "--- Statistics::Statistics(QWidget *)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
 Statistics::Statistics(ivg::Session *session_ptr, QWidget *parent): Statistics(parent)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Statistics::Statistics(ivg::Matrix, map, map, map, QWidget)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   _session = session_ptr;
   _residuals = _session->get_residuals();
   if(_session->getAmbigRes()){
      ui->errorbar_checkbox->setChecked(false);
   }
       
    // create a new model
    QStandardItemModel *model = new QStandardItemModel(0,4,this);
  
    // set header of columns from model
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("Object")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("#in")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QString("#out")));
    model->setHorizontalHeaderItem(3, new QStandardItem(QString("wrms")));
  
    for(auto const &it : _residuals[delaytype::group] ) 
    {   
        // generate rows for table view from model
        QList<QStandardItem*> new_row;
        new_row.append(new QStandardItem(QString::fromStdString(it.name)));
        new_row.append(new QStandardItem(QString::number(it.idx_in.size())));
        new_row.append(new QStandardItem(QString::number(it.idx_out.size())));
        new_row.append(new QStandardItem(QString::number(it.wrms)));
        model->appendRow(new_row);
    }
   
    // attach the model to the view
    ui->data_tableview->setModel(model);
   
    // fit generated table to view
    ui->data_tableview->resizeColumnsToContents();
    ui->data_tableview->horizontalHeader()->setStretchLastSection(true);

    // plot sites, sources, baselines
    _plot_sites();
  
    _plot_baselines();

    _plot_sources();
    
   //open widget when everything is done
    this->show();
    
#if DEBUG_VLBI >=2
   cerr << "--- Statistics::Statistics(ivg::Matrix, map, map, map, QWidget)" << " : " << tim.toc() << " s " << endl;
#endif   
}
// ...........................................................................
void Statistics::_plot_sites()
// ...........................................................................
{       
    ivg::Matrix lat;
    ivg::Matrix lon;
    vector<string> tooltips;
    for(auto &sta: (*_session->get_trf_ptr()))
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
    
    ui->trf_widget->projection(_trf_projection, lat, lon, { Qt::red, 5.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDiamond, "Sites", 4.0, 6.0 }, "Terrestrial Reference Frame", tooltips);
}
// ...........................................................................
void Statistics::_plot_sources()
// ...........................................................................
{
    struct crfdata{ivg::Matrix ra; ivg::Matrix dec; ivg::Matrix d_ra; ivg::Matrix d_dec; vector<string> tooltips;};
    map<string,crfdata> crf_info;

    ivg::Crf crf = (*_session->get_crf_ptr());

    ivg::Matrix ra,dec,d_ra,d_dec;
    for(vector<ivg::Source>::iterator src = crf.begin(); src != crf.end(); ++src)
    {
        // ignore trash-sources (only generated within scheduling)
        if(src->use_me())
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

        ui->crf_widget->projection(_crf_projection, srces.second.ra, srces.second.dec, format ,"Celestial Reference Frame", srces.second.tooltips );
    }
}
// ...........................................................................
void Statistics::_plot_baselines()
// ...........................................................................
{    
    // go through all baselines and split by name
    for(auto const &it : _residuals[delaytype::group])
    {
        if(it.type == residtype::baseline) 
        {
            // split by '-' delimter to ge station names
	  string sta1,sta2,tmp;
            stringstream sta1sta2(it.name);
            getline(sta1sta2, sta1, '-');
            getline(sta1sta2, sta2, '-');
	    // QAD solution  for VLBA stations, e.g. BR-VLBA. Better solution, including other stions with - in the name, is still needed
	    if (sta2 == "VLBA") {
	      sta1=sta1+"-"+sta2;
	      getline(sta1sta2, sta2, '-');
	    }
	    getline(sta1sta2, tmp, '-');
	    if (tmp == "VLBA") {
	      sta2=sta2+"-"+tmp;
	    }
            // get both site positions for plotting basline between them
            ivg::Matrix llh;
            ivg::Analysis_station *sta1_ptr,*sta2_ptr;
            _session->get_trf_ptr()->get_station(&sta1_ptr,sta1);
            _session->get_trf_ptr()->get_station(&sta2_ptr,sta2);
	    
            llh.append_rows(sta1_ptr->calc_lat_lon_h().transpose());
	    
            llh.append_rows(sta2_ptr->calc_lat_lon_h().transpose());
	    
            ui->trf_widget->projection(_trf_projection, llh(":",1), llh(":",0), { Qt::red, 1.5, Qt::DotLine, QCPGraph::lsLine, QCPScatterStyle::ssDiamond, QString::fromStdString(it.name), 0.0, 0.0 }, "Terrestrial Reference Frame");  
        }
    }   
}
// ...........................................................................
void Statistics::createFigure()
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Statistics::createFigure()" << endl; 
   tictoc tim;
   tim.tic();
#endif
      
    Plot *hist = new Plot(ui->hist_dockwidget);
    Plot *dialog = new Plot(ui->plot_dockwidget);
    
    QModelIndexList indexes = ui->data_tableview->selectionModel()->selectedRows();
    
    if(!ui->histogram_checkbox->isChecked())
        ui->hist_dockwidget->hide();
    else
        ui->hist_dockwidget->show();
    
    ui->plot_tabwidget->clear();
    
    int subplot = 0;
    if(ui->boxplot_checkbox->isChecked() && ui->histogram_checkbox->isChecked())
        subplot = hist->add_rowplot("boxplot","key");
    else if(ui->boxplot_checkbox->isChecked() && !ui->histogram_checkbox->isChecked())
        ui->hist_dockwidget->show();
    else if(ui->crosscorr_checkbox->isChecked() && indexes.count() == 2 
            && !ui->boxplot_checkbox->isChecked() && !ui->histogram_checkbox->isChecked())
        ui->hist_dockwidget->show();    
        
    ui->hist_dockwidget->setWidget(hist);
    
    // create single-band and multi-band plots
    // _residuals[delaytype::group] and _residuals[delaytype::single]
    for(int resid_i=0; resid_i<_residuals.size(); resid_i++)
    {
        if( resid_i == 0 )
        {
            Plot *mb = new Plot(ui->plot_dockwidget);
            ui->plot_tabwidget->addTab(mb, QString::fromUtf8("Multi-Band"));
            dialog = mb;
            _last_dialog = dialog;
        }
        else
        {
            Plot *sb = new Plot(ui->plot_dockwidget);
            ui->plot_tabwidget->addTab(sb, QString::fromUtf8("Single-Band"));
            dialog = sb;
        }
        
        QCustomPlot * plot = dialog->get_plot();

        connect(plot, SIGNAL(mousePress(QMouseEvent*)), this, SLOT(showPointToolTip(QMouseEvent*))); 

        string xaxis_combobox = ui->xaxis_combobox->currentText().toLatin1().data();
        string yaxis_combobox = ui->yaxis_combobox->currentText().toLatin1().data();

        // find specific axis layout for selected x/y axis
        AxisLayout xaxis,yaxis;
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

        // only if one station is selected plot skyplot
        if(indexes.count() == 1)
        {
            QString selection =  ui->data_tableview->model()->data(indexes.at(0),Qt::DisplayRole).toString();
            plotSiteSkyplot(selection);
        }
        if(ui->crosscorr_checkbox->isChecked() && resid_i == 0 && indexes.count() == 2)
        {          
            string bl1 =  ui->data_tableview->model()->data(indexes.at(0),Qt::DisplayRole).toString().toStdString(); 
            string bl2 =  ui->data_tableview->model()->data(indexes.at(1),Qt::DisplayRole).toString().toStdString();

            Plot plotMaxCorr;
            ivgat::Tsa tsa1;
            ivg::Matrix tout, mout;
            ivg::Matrix data1,data2;    
            std::string name1,name2;

            for( auto & it : _residuals[delaytype::group] )
            {
                if( it.type == residtype::baseline )
                {
                    if( it.name == bl1 )
                    {
                        data1 = it.data;
                        name1 = it.name;
                    }
                    else if( it.name == bl2 )
                    {
                        data2 = it.data;
                        name2 = it.name;
                    }
                }
            }

            tsa1.crosscorrNEQD( data1.get_col(0), data1.get_col(1), data2.get_col(0), data2.get_col(1), 100, tout, mout);
            
            hist->plot_data(tout, mout,{
                            QColor(ivg::color_values.at(7).c_str()), 1.0, Qt::SolidLine,
                            QCPGraph::lsLine, QCPScatterStyle::ssDot, ("cross-correlation"), 3.0, 12.0 }, 0);                     
        }

        for (int i = 0; i < indexes.count(); ++i)
        {
            // get new color for each new graph
            QColor plot_color = QColor(ivg::color_values.at(i).c_str());

            // detect which objects have been selected (e.g. WETTZELL and 0059+581)
            QModelIndex index = indexes.at(i);
            string selection =  ui->data_tableview->model()->data(index,Qt::DisplayRole).toString().toStdString();

            for(auto const &it : _residuals[(delaytype) resid_i] ) 
            { 
                if(selection == it.name)
                {                    
                    //if outliers selected, add them
                    vector< vector<int> > idx_inout = { it.idx_in };
                    if( it.idx_out.size() > 0 && ui->outlier_checkbox->isChecked() )
                        idx_inout.push_back(it.idx_out);

                    int loop=1;
                    
                    // idx_inout can contain two vector<int> (idx))
                    // the first contains the inliers
                    // the second the outliers
                    for(auto const &idx : idx_inout ) 
                    {
                        ivg::Matrix x = it.data.get_sub( idx, {xaxis.data_idx} ) * xaxis.fak;
                        ivg::Matrix y = it.data.get_sub( idx, {yaxis.data_idx} ) * yaxis.fak;
                        ivg::Matrix sig = it.data.get_sub( idx, {5} )*1e12;
                        
                 
                        // TEST OUTPUT to check indices
//                        std::cerr << "loop: " << loop << " x.size "<< x.size(1) << " " << x.size(2)  <<  std::endl;
//                        std::cout << "figure idx "; show_vector(idx);
//                        std::cout << "ncfile idx "; show_vector(it.idx_ncfile);

                        if(ui->absolute_checkbox->isChecked())
                            y = y.abs();

                        if(ui->errorbar_checkbox->isChecked())
                        {
                            dialog->plot_yerr( x, y, sig, { plot_color, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, selection.c_str(), 1.0 , 4.0  });

                            //if outliers, plot the errorbars red
                            if(loop == 2)
                                plot->graph()->setErrorPen(QPen(QColor(Qt::red)));
			    
                        }
                        else
                        {
                            // manual ambiguity resolution
                            // the indices of the observations in the nc file have to be connnected with the datatpoints in the plot
                            // After creating a graph a pointer to the graph is returned ( QCPGraph * graph_ptr )
                            // this pointer is the key for a map, that contains for each graph the original indices
                            // the map is a private attribute in the plot class, because when creating new garphs out of existing graphs the infoprmation is needed there
                            
                            // plot graph amd receive the pointer to the graph
                            QCPGraph * graph_ptr = dialog->plot_data( x, y, { plot_color, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, selection.c_str(), 1.0 , 4.0  });
                            
                            // insert into map (private attribute of class plot)
                            (*dialog->get_index_ptr())[graph_ptr] = it.idx_ncfile;
                        }

                        if(ui->histogram_checkbox->isChecked() && resid_i == 0)
                        {
                            int bins;
                            if(ui->histogram_spinbox->value() == 0)
                                bins = (int) 1+3.3 * log(y.length())/log(10);
                            else
                                bins = ui->histogram_spinbox->value();

                            double spacing = hist->histogram(y,bins,  { plot_color, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, "", 1.0 , 4.0  });
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

                            cerr << "Mean: " << mean << " Std: " << sigma << endl;
                            hist->plot_data(keys,tmp,  { plot_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, selection.c_str(), 1.0 , 4.0  });
                        }

                        if(ui->boxplot_checkbox->isChecked() && resid_i == 0)
                        {
                            hist->boxplot(y,{ plot_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, "", 1.0 , 4.0  }, subplot);
                        }
                        loop++;
                    }
                }  
            }   
        }

    //    hist->get_plot()->rescaleAxes();
        hist->get_plot()->replot();
        plot->rescaleAxes();
        plot->replot();
        ui->plot_tabwidget->setCurrentIndex(ui->plot_tabwidget->indexOf(_last_dialog));
    }
    
    ui->cbreak_pushbutton->setEnabled(true);
    
#if DEBUG_VLBI >=2
   cerr << "--- void Statistics::createFigure()" << " : " << tim.toc() << " s " << endl;
#endif  
}



// create two new graphs with STRG+CAPS than use the CBreak Button

// ...........................................................................
void Statistics::correctCBreak()
// ...........................................................................
{
    
    Plot *dialog = (Plot*) ui->plot_tabwidget->widget(0);
    QList<QCPGraph*> graph_list = dialog->get_plot()->selectedGraphs();

    if( graph_list.size() == 2 )
    {
        
        QDateTime cbTime;
        
        // loop over selected graphs
        std::vector<double> mean(2, 0.0);
        for( unsigned short i = 0; i < 2; ++i){ 
            QCPDataMap* data = graph_list.at(i)->data(); 

            // loop over all pairs in QCPDataMap data
            for(auto &point: (*data))
            {
                mean[i] += point.value;
            }
            mean[i] /= data->size();
            std::cout << graph_list.at(i)->name().toStdString() << " mean " << mean[i] << std::endl;
        }
        
        cbTime.setTime_t( (graph_list.at(0)->data()->end()-1).key() );
        //cbTime.setTimeSpec(Qt::UTC);
        double clockBreak = mean[0]-mean[1];
        
        std::cout << "  clockbreak time: " << cbTime.toUTC().toString("dd.MM.yyyy. HH:mm:ss.zzz t").toStdString() << "   offset: " << clockBreak << std::endl;
        
    
        ivg::Date d( cbTime.date().year(), cbTime.date().month(),cbTime.date().day(),
                     cbTime.time().hour(), cbTime.time().minute(), 
                     double(cbTime.time().second()) + double(cbTime.time().msec())  );
              
        
        _clock_break_station.push_back(graph_list.at(0)->name().toStdString());
        _clock_break_epoch_mjd.push_back(d.get_double_mjd());
        std::cout << "              mjd: " <<  std::setprecision(12) <<d.get_double_mjd() << endl;
    
        
    } else {
        std::cerr << "You have to select 2 graphs using STRG" << std::endl;
    }
    
    
} 


// ...........................................................................
void Statistics::plotSiteSkyplot(QString &tooltip)
// ...........................................................................
{  
    // use already initialized sky-plot for plotting of the individual residuals
    for( int i=0; i<=ui->sky_widget->get_plot()->graphCount(); ++i )
        ui->sky_widget->get_plot()->removeGraph(ui->sky_widget->get_plot()->graph(i));
    
    // get ivs_site name of selected station, e.g. GILCREEK
    string ivs_site = remove_spaces_end(tooltip.toStdString());
   
    // only proceed if it's not a baseline or empty
    if(ivs_site.find("-") == string::npos && !ivs_site.empty() && ivs_site != "ALL")
    {
        // to get sure
        ivs_site = ivs_site.substr(0,8);
        
        // take all residuals of the slected site and save azimuth and elevation
        vector<int> indexes;
        ivg::Matrix azi,ele;
        for(auto const &it : _residuals[delaytype::group] ) 
        { 
            if(it.type == residtype::station && ivs_site == it.name)
            {
                indexes = it.idx_origin;   
                for(int r=0; r<it.data.rows(); r++)
                {
                    // check whether selected staion was "sta1" or "sta2" to get correct azi/ele values
                    if(ivs_site == it.first_station.at(r))
                    {
                        azi.append_rows(it.data(r,7));
                        ele.append_rows(it.data(r,8));
                    }
                    else
                    {
                        azi.append_rows(it.data(r,9));
                        ele.append_rows(it.data(r,10));
                    }
                }
                
                // plot corresponding sources with azimuth and elevation from selected site
                QFormat format =  { Qt::red, 1.5, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssDisc, QString("Sources"), 1.0 , 4.0 };
                ui->sky_widget->projection(_sky_projection, azi , ele, format ,"Skyplot "+ivs_site, it.source_names);   
                break;
            }
        }
    }
    else if( ivs_site == "WETTZ13N-WETTZELL" && !ivs_site.empty() && ivs_site != "ALL" )
    {
        ivg::Matrix data;
        for( auto const & it : _residuals[delaytype::group] )
        {
            if( it.type == residtype::baseline && ivs_site == it.name )
                data = it.data;
        }      

        ivg::Matrix res = data.get_col(1).absD() * 1e12;
        int bins = 10;
        int dt = ( res.max()-res.min() ) / bins;

        ivg::Matrix idx(bins,1,0.0);
        idx(0,0) = res.min(); idx(bins-1,0) = res.max();
        for( int j=1; j<bins-1; ++j )
            idx(j,0) = idx(j-1,0) + dt;

        for( int i=1; i<idx.rows(); ++i )
        {
            vector<int> tmp = res.find_idx( gt, idx(i-1), le, idx(i) );
            stringstream ss;
            ss << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(5) << (int)idx(i-1)<<" [ps] < r < " << setw(5) << (int)idx(i) << " [ps]";

            if( tmp.size() > 0 )
            {                            
                QFormat format = { QColor( ivg::color_values.at( i ).c_str() ), 1.5, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, QString::fromStdString(ss.str()), 1.0 , 4.0 };
                ui->sky_widget->projection(_sky_projection, data.get_col(7).get_rows(tmp), data.get_col(8).get_rows(tmp), format, "Skyplot of post-fit residuals"); 
            }
            else
            {
                ivg::Matrix foo(1,1,0.0);
                QFormat format = { QColor( ivg::color_values.at( i ).c_str() ), 1.5, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssNone, QString::fromStdString(ss.str()), 1.0 , 1.0 };
                ui->sky_widget->projection(_sky_projection, foo, foo, format, "Skyplot of post-fit residuals"); 
            }
        }        
    }
    else
    {
        // skyplot without plotting a point in case of baselines and sources
        ui->sky_widget->projection(_sky_projection, ivg::Matrix(1,1,0.0)*(M_PI/180.0), ivg::Matrix(1,1,30.0)*(M_PI/180.0), 
                    { Qt::red, 5.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssNone, "NONE", 4.0, 6.0 }, "Skyplot");
    }
}
// ...........................................................................
void Statistics::showBaselineInfo(QCPAbstractPlottable *plottable, QMouseEvent *event)
// ...........................................................................
{    
    // get current crf plot
    Plot * crf_plot = ui->crf_widget;

    // only remove last graph if already one plotted before
    if(_bl_selected == true)
        crf_plot->get_plot()->removeGraph(crf_plot->get_plot()->graph());

    // get corresponding observation indexes to selected baseline
    vector<int> indexes;
    QCPGraph *graph= (QCPGraph*)plottable;
    for(auto const &it : _residuals[delaytype::group] ) 
    { 
        if(graph->name().toStdString() == it.name)
        {
            indexes = it.idx_origin;   
            break;
        }
    }

    // find out which sources are related to the same observations
    vector<string> sources;
    for(auto const &it : _residuals[delaytype::group] ) 
    { 
        if(it.type == residtype::source)
        {
            for(auto &idx: indexes)
            {
                vector<int> tmp = it.idx_origin;
                if(std::find(tmp.begin(),tmp.end(), idx) != tmp.end())
                    sources.push_back(it.name);
            }
        }
    }
    // delete sources which have been found several times
    remove_duplicates(sources);

    // create matrices containing positions of corresponding sources
    ivg::Matrix ra,dec;
    for(auto &src: sources)
    {
        ivg::Source *src_ptr;
        if(_session->get_crf_ptr()->get_source(&src_ptr,src))
        {
            ra.append_rows(src_ptr->get_ra0() - M_PI);
            dec.append_rows(src_ptr->get_dec0());
        }
    }

    // plot corresponding sources
    QFormat format =  { Qt::black, 1.5, Qt::SolidLine, QCPGraph::lsNone,  QCPScatterStyle::ssCircle, graph->name(), 2.0 , 6.0 };
    crf_plot->projection(_crf_projection, ra, dec, format ,"Celestial Reference Frame" );    
    
    _bl_selected = true;
}
// ...........................................................................
void Statistics::showPointToolTip(QMouseEvent *event)
// ...........................................................................
{
    QCustomPlot * plot = _last_dialog->get_plot();
    
    double x = plot->xAxis->pixelToCoord(event->pos().x());
    double y = plot->yAxis->pixelToCoord(event->pos().y());
    
    QList<QCPGraph*> selected_graphs = plot->selectedGraphs();
        
    
    // get all data plotted in the figure
    for(auto const &it : selected_graphs ) 
    { 
        if( it->selectTest(event->pos(), false) < plot->selectionTolerance() )
        {
            string graph_name = it->name().toLatin1().data();
            
                    
            for(auto const &resid : _residuals[delaytype::group] ) 
            { 
                if(resid.name == graph_name)
                {
                    QCPDataMap *dataMap = it->data();

                    double a = abs((dataMap->lowerBound(x)-1).key() - x);
                    double b = abs((dataMap->lowerBound(x)).key() - x);
                    
                    stringstream ss;
                    
                    double key;
                    if(b < a)
                        key = (double)(dataMap->lowerBound(x)).key();
                    else
                        key = (double)(dataMap->lowerBound(x)-1).key();
                    
                    double xfak,yfak;
                    for(auto const &layers : _axis_layouts ) 
                    {
                        if(QString::fromStdString(layers.name + layers.unit) == plot->xAxis->label())
                            xfak = layers.fak;
                        if(QString::fromStdString(layers.name + layers.unit) == plot->yAxis->label())
                            yfak = layers.fak;
                    }
                    
                    vector<int> indexes_vec  = resid.data.find_idx( key / xfak  );
                    
                    vector<int> rows_vec;
                    for(int i=0; i<indexes_vec.size(); i++)
                            rows_vec.push_back(indexes_vec.at(i) % resid.data.rows());
                        
                    
                    ivg::Matrix candidates = resid.data.get_rows(rows_vec);
                    ivg::Matrix diff = candidates - ivg::Matrix(candidates.rows(),candidates.cols(), (y/yfak));
                    ivg::Matrix abs = diff.absD();
                    
                    int row = abs.minIdx() % abs.rows();
                    
                    for(auto const &layers : _axis_layouts )
                    {
                        if(layers.name == "Time - hh:mm:ss")
                        {
                            ivg::Date date( 1970, candidates(row,layers.data_idx)/86400.0 + (7200.0/86400.0) );
                            std::streamsize p = ss.precision();
                            ss << "Time : " << date.get_date_time("DD MON YYYY / HH:MI:SS") << endl;
                            ss << "MJD : " << setprecision(10) << fixed << date.get_double_mjd() << endl;
                            ss << setprecision(p);
                        }
                        else if(layers.name != "Time - dd-MMM-yy")
                            ss << layers.name << layers.unit << ": " << candidates(row,layers.data_idx) * layers.fak << endl;
                    }
                    
                    QToolTip::showText(QCursor::pos(), QString::fromStdString(ss.str())); //QString("%1 , %2").arg((double)(dataMap->lowerBound(x)).key()).arg(y));
                }
            }  
        }
    }

}
// ...........................................................................
Statistics::~Statistics()
// ...........................................................................
{
  delete ui;
}
// ...........................................................................
void Statistics::_init_axis_layout()
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
    
}


// ...........................................................................
 void Statistics::_shiftResid(char direction){
// ...........................................................................
     
    Plot *dialog = (Plot*) ui->plot_tabwidget->widget(0);
    QList<QCPGraph*> graph_list = dialog->get_plot()->selectedGraphs();
   
    
    if( graph_list.size() == 1 )
    {
 
        QCPDataMap* data = graph_list.at(0)->data();
        

        double ambiguity_spacing = _session->get_max_ambiguity_spacing()*1e12;
        int sign = 0;
        

        if (direction == 'U'){
             std::cerr << "shifting up " << ambiguity_spacing << " ps" << std::endl;
             sign = 1;
        } else if (direction == 'D'){
             std::cerr << "shifting down " << ambiguity_spacing << " ps" << std::endl;
             ambiguity_spacing *= -1;
             sign = -1;
        }
        
        // original indices
        std::vector<int>* idx_orig = &(*dialog->get_index_ptr())[graph_list.at(0)];
//        show_vector((*idx_orig));
	
        QVector<double> x;
        QVector<double> y;
         
        //loop over all points in QCPDataMap data and add the ambiguity spacing
        int c = 0;
        for(QCPData &point: (*data))
        {
            
//            std::cout << idx_orig->at(c) << " "<< point.value << std::endl;
            x.push_back( point.key );
            y.push_back( point.value + ambiguity_spacing);
            
            ++c;
        }
        
        // refresh the graph
        graph_list.at(0)->setData(x,y);       

	dialog->get_plot()->replot();
        
        ambiguity_spacing *= 1e-12;
        
        // update num_ambig Matrix in session which connects the 'solve_ambiguities' program with the GUI 
        ivg::Matrix*  num_ambig_ptr = _session->get_num_ambig_ptr();
        //std::cerr << "size num_ambig Mat " << num_ambig_ptr->size(1) << " " << num_ambig_ptr->size(2) << std::endl;
        
        for( auto i: *idx_orig) {
            (*num_ambig_ptr)(i) += sign;
        }
        
        
        //update the residuals        
        // loop over all elemenmts in vector<Residual>
        for( Residual& r :    _residuals[delaytype::group] ){
           
            std::vector<int>::iterator low;
            
//            std::cout << r.name << std::endl;
            
            int d = 0;
            // loop over all original indices in the Residual struct and compare than with the original indices in the selected graph
            for( unsigned int i = 0; i < r.idx_ncfile.size(); ++i  ){
               
                low = std::lower_bound (  idx_orig->begin()+d,   idx_orig->end(), r.idx_ncfile.at(i) ); 

                if(low != idx_orig->end() && !(r.idx_ncfile.at(i)  < *low) ){
                    
                    d = distance (idx_orig->begin(), low); //position of the found element idx_orig vector of selected graph
                    
//                    std::cout << "index " << r.idx_ncfile.at(i) << " has been found at position " << d << "( " << idx_orig->at(d) << " )" << std::endl; 
//                    std::cout << r.data(i,1)*1e+12 << std::endl;
                   
                    r.data(i,1) += ambiguity_spacing;
                }
                
            }
        }
      
            
    } else {
        std::cerr << "selcet exact one graph in multi-band delay window" << std::endl;
    }

        
 }
 
 // ...........................................................................
 void Statistics::saveResid(){
// ...........................................................................
     //_session->get_residuals() = _residuals;
     
     if( _clock_break_station.size() > 0 ){
        ivg::Vgosdb vgos( _session->get_session_path() );
        vgos.create_ClockBreak_file(_clock_break_station, _clock_break_epoch_mjd, _session->get_name());
     }
     
     
 }


