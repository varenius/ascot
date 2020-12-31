#include "lp_sked_gui.h"

namespace lps{

    
GUI::GUI(ivg::Session* session, QObject *parent) : QObject(parent)
{
    this->session = session;
    this->twinMap = session->get_trf_ptr()->get_twins_map();
        
    lps::Seconds duration = (session->getEnd().get_double_mjd()-session->getStart().get_double_mjd())*3600*24;
    lps::Seconds temporal_grid_resolution = (int)(*session->get_setup())["SKED"]["temporal_resolution"] * 60;
    lps::Seconds temporal_grid_shift = (int)(*session->get_setup())["SKED"]["temporal_shift"] * 60;
    this->tg = TemporalGrid(duration, temporal_grid_resolution, temporal_grid_shift, session->getStart());
}

int GUI::run(int argc, char *argv[]){
    
    bool show = (bool)(*session->get_setup())["SKED"]["createPlots"];
    
    QApplication a(argc, argv);

    // initialize worker (Worker is inherited from QObject)
    lps::Worker worker(session);
    worker.moveToThread(&workerThread);
        
    connect(&workerThread, SIGNAL (started()), &worker, SLOT (process()));
    
    if (show){
        initGUI();
        
        connect(&worker,SIGNAL(foundSolution(std::vector<lps::StationActivity>)),
            this, SLOT( foundSolution(std::vector<lps::StationActivity>) ));
        connect(&worker,SIGNAL(selectedCell(int, int,lps::Path,int,bool)),
            this,SLOT(selectedCell(int, int, lps::Path,int,bool)));
        
        connect(&workerThread, SIGNAL (finished()), SLOT (print_skyplots()) );
        
    } else {
        connect(&workerThread, SIGNAL (finished()), &a, SLOT (quit()) );
    }
  
   
    // calls the function Worker::process(){
    workerThread.start();
        
    return a.exec();
}


 void GUI::initGUI(){
    bool autoClose = (bool)(*session->get_setup())["SKED"]["auto_close_windows"];
    double min_ele = (double)(*session->get_setup())["SKED"]["min_elevation"]; // in deg
    std::string outdir = (const char*)(*session->get_setup())[ "outdir" ];
    
    // create plot for each station  
    for(size_t sta_idx = 0; sta_idx <  session->get_trf_ptr()->get_number_stations(); ++sta_idx){

        ivg::Analysis_station* sta = session->get_trf_ptr()->get_station(sta_idx);

        stationViews.push_back(new StationDialog(session, &tg, nullptr) );
        StationDialog * sd = stationViews.back();

        if(!autoClose)
            sd->show();

        sd->addTransits(sta_idx, session->get_trasintsVisibleFromStation(), QPen(Qt::lightGray,2) );
        sd->addTransits(sta_idx, session->get_trasintsVisibleFrom2Stations(), QPen(QColor(95,95,107),2) );
        sd->addTransits(sta_idx, session->get_complete_transits(), QPen(QColor(135,206,235),2) );
        sd->addTransits(sta_idx, session->get_common_transits(), QPen(Qt::blue,2) );

        sd->setWindowTitle(QString(sta->get_name(ivg::staname::ivs_name).c_str()));
        sd->set_default_path( outdir + "/"+ sta->get_name(ivg::staname::ivs_name) + ".pdf" );
        ivg::Matrix lat_lon_h = sta->calc_lat_lon_h();
        double elePol =  lat_lon_h(0)*ivg::rad2d;
        double azPole = 0;
        if(elePol < 0){
            azPole = 180.0;
            elePol *=-1;
        }
        sd->drawPoint( azPole, elePol);

        // Elevation mask
        std::vector<lps::Position> elevationMask;
        for( int azi = 0; azi <= 360; ++azi ){

            double ele = max(min_ele, sta->get_ele_mask(azi*ivg::d2rad)*ivg::rad2d);
            elevationMask.push_back( Position(azi, ele, session->getStart())  );
        }
        sd->addElevationMask(elevationMask, QPen(Qt::black,2.5) );

        // Sun Transit
        std::vector<lps::Position> sunPath;
        ivg::Date ep = session->getStart();
        while(ep < session->getEnd()){
            // get current crf2trf and partials matrix for each scan
            ivg::Partials_t2c tmp { ivg::Matrix( 3,3,0.0 ) };
            ivg::Partials_t2c * deriv_ptr = &tmp;

            ivg::Matrix crf2trf = session->get_eops().form_crf2trf(ep, true, deriv_ptr );

            ivg::Matrix azel = sta->compute_az_el_sun(ep, session->get_ephem(), crf2trf)*ivg::rad2d;
            if(azel(1) > 0){
                sunPath.push_back( Position(azel(0), azel(1), ep) );
            }
            ep.add_secs(300);
        }
        sd->addSunTransit(sunPath,QPen(QColor(255,127,80,255),3));
    }
}

void GUI::showSession(){
    
    
    bool ea_grid  = (*session->get_setup())["SKED"].exists("EqualAreaGrid");
    
    std::vector<std::vector<unsigned>> ea_grid_setup;
    std::vector<unsigned> n_segments; 
    if (ea_grid){
        Setting &ea = (*session->get_setup())["SKED"]["EqualAreaGrid"];
        n_segments.resize(ea.getLength());
        ea_grid_setup.resize(ea.getLength());
    
        for(int i=0; i<ea.getLength(); i++){
            for(int j=0; j<ea[i].getLength(); j++){
                ea_grid_setup[i].push_back( (unsigned)ea[i][j] );
                n_segments[i] += (unsigned)ea[i][j];
            }
        }
    }
    

    
    // add Observations
    std::vector<std::vector < ObservationItem*>> observations;
    observations.resize(session->get_crf_ptr()->get_number_sources_inc_unused());

    for (ivg::Scan& scan : *session->get_scan_ptr()) {
        ivg::Source* const sou = scan.get_source();
        ivg::Date d = scan.get_epoch();

        lps::Seconds start = (d.get_double_mjd()-session->getStart().get_double_mjd())*3600*24;
  
        for (ivg::Analysis_station* station : scan.get_stations()){
            for( ivg::Analysis_station* sta: twinMap[station->get_idx()] ){
                StationDialog * sd = stationViews[sta->get_idx()];

                ivg::Matrix azel = sta->calc_az_el(d, *sou);
                azel*=ivg::rad2d;
                lps::Position pos ( azel(0), azel(1), d);

                ObservationItem * item = sd->addObservation(*sou, pos.azimuth(),pos.elevation(), start, tg.getTemporalIndices(start), sta->get_idx() == station->get_idx());
                observations[sou->get_idx()].push_back(item);
                if(sta->get_idx() != station->get_idx()){
                   item->setPen(QPen(Qt::magenta,2));
                }
            }
        }
    }
    
        
    double objective = 0;

    std::vector< std::vector< lps::Node<ivg::Obs*, lps::Wedge>  > >  roots;

    std::map<int, std::vector<double> > coverage = session->compute_sky_coverage(ea_grid_setup, tg, roots, objective);

    unsigned observatory_idx = 0;
    for(std::vector<ivg::Analysis_station*> observatory : session->get_trf_ptr()->get_station_twin_list() ){
        for(unsigned int t = 0; t < tg.get_number_of_intervals(); ++t){

            // coverage
            roots[observatory_idx][t].visitBoundAndBool([&]( int level, const lps::Wedge& rect, bool valid ){                
                for(ivg::Analysis_station* sta : observatory ){
                    StationDialog* sd = stationViews[sta->get_idx()];
                    QGraphicsPathItem * item = sd->addCell(rect.toPath(0.5),level,t ,valid);
                    item->setZValue(-1);
                    item->setBrush(QBrush(QColor(0,0,255, 20 )));
                }

            });
        }
        observatory_idx++;
    }

    for(StationDialog * sd : stationViews){
        sd->addTraverse( QPen(Qt::red,1) );
    }

    session->show_coverage(coverage);
    
     // save objective info to file
    std::ofstream objfile;
    std::stringstream ss;
    ss << (const char*)(*session->get_setup())[ "outdir" ] << "/"
    << session->get_name() << + "_objective_info.txt";
    std::string objfile_name = ss.str();

    objfile.open( objfile_name, ios::out);

    objfile << setprecision(12) << objective << std::endl;
    objfile << session->get_nobs_orig() << std::endl;
    objfile << session->coverage2string(coverage,10);
    objfile.close();

    
} 
 
void GUI::foundSolution(std::vector<lps::StationActivity> activity)
{   
    for(StationDialog * sd : stationViews){
        sd->clear();
    }

    std::vector<std::vector<ObservationItem*>> observations;
    observations.resize(session->get_crf_ptr()->get_number_sources_inc_unused());

    std::vector< std::set<unsigned >> added_intervals( session->get_trf_ptr()->get_number_stations() );

    for(const StationActivity & sa : activity){

        // Do not at dublicate observation. Due to data structure of StationActivity vector (each baseline is represented with two elements in this vector)
        // an observation (point on the map) would be added several times without this
        std::pair<std::set<unsigned>::iterator,bool> ret = added_intervals[sa.station_idx].insert(sa.index);
        if(ret.second){
    //        sa.print();
            ivg::Source* sou = session->get_crf_ptr()->get_source( sa.source_idx);
            ivg::Analysis_station* station = session->get_trf_ptr()->get_station( sa.station_idx );
            ivg::Date d = session->getStart();      
            d.add_secs(sa.begin);

            ivg::Matrix azel = station->calc_az_el(d, *sou);
            azel*=ivg::rad2d;
            lps::Position pos ( azel(0), azel(1), d);

            for( ivg::Analysis_station* sta: twinMap[sa.station_idx] ){
                StationDialog * sd = stationViews[sta->get_idx()];                
                ObservationItem * item = sd->addObservation(*sou, pos.azimuth(),pos.elevation(), sa.begin, tg.getTemporalIndices(sa.begin), sta->get_idx() == sa.station_idx);
                observations[sou->get_idx()].push_back(item);
                if(sta->get_idx() != sa.station_idx){
                    item->setPen(QPen(Qt::magenta,2));
                }
            }
        }
    }

    for(std::vector<ObservationItem*> & items : observations){
        for(ObservationItem * item : items){
            for(ObservationItem * copy : items){
                assert(item->quasar() == copy->quasar());
                item->copies.push_back(copy);
            }
        }
    }

    for(StationDialog * sd : stationViews){
        sd->addTraverse( QPen(Qt::red,1) );
    }
}

void GUI::selectedCell(int level, int temporal_idx, lps::Path rect, int sta_idx, bool valid)
{
    StationDialog * sd = stationViews[sta_idx];
    QGraphicsPathItem * item = sd->addCell(rect, level, temporal_idx, valid);
    item->setZValue(-1);
    item->setBrush(QBrush(QColor(0,0,255,20)));
}

void GUI::print_skyplots()
{
    
    bool autoClose = (bool)(*session->get_setup())["SKED"]["auto_close_windows"];
    
    // save pdf
     for(size_t sta_idx = 0; sta_idx <  session->get_trf_ptr()->get_number_stations(); ++sta_idx){

        StationDialog* sd = stationViews[sta_idx];

        std::stringstream ss;
        ss << (const char*)(*session->get_setup())[ "outdir" ] << "/"
           << session->get_name() << "_" << session->get_trf_ptr()->get_station(sta_idx)->get_name(ivg::staname::ivs_name) << ".pdf"; 
        std::string path = ss.str();

        sd->print_pdf( path );
        
        if(autoClose)
            sd->close();

    }
     
   

}

}
