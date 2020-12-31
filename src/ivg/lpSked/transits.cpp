
/* 
 * File:   transits.cpp
 * Author: corbin
 * 
 * Created on April 25, 2018, 9:52 AM
 */

#include "transits.h"
namespace lps{

Transits::Transits() {
}

Transits::Transits( ivg::Session *session_ptr, unsigned int dt ){
    
    log<DETAIL>("*** computing transits");
    
    _session_ptr = session_ptr;
    Setting* const setup = _session_ptr->get_setup();
        

    // read required quantities from controlfile
    double min_ele = (double)(*setup)["SKED"]["min_elevation"]*ivg::d2rad;
    if (min_ele < 0) {
        log<WARNING>("!!! min_elevation: ") % min_ele % " has to be larger than zero. Setting it to 5 deg";
        min_ele = 5*ivg::d2rad;
    }
    
    double min_dist_2_sun = (double)(*setup)["SKED"]["min_sun_dist"];
    
    if( dt == 0){
        _dt = (unsigned int) (*setup)["SKED"]["transit_interval_length"]; 
    } else {
        _dt = dt;
    }
    
    
    // Compute for each station the visible source trails
    
    // incremented  for each visible source/station combination transit
    unsigned int transit_counter = 0;
    
    // loop over all stations
    for( unsigned int sta_idx = 0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx ){
        ivg::Analysis_station* sta = _session_ptr->get_trf_ptr()->get_station(sta_idx);

        // loop over all sources
        for( unsigned int sou_idx = 0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx ){
            ivg::Source* sou = _session_ptr->get_crf_ptr()->get_source(sou_idx);
            
            // first position of a transit corresponding to the source/station combination
            unsigned int first = transit_counter;

            // maximal number of elements in path vector
            unsigned int n = ceil( (_session_ptr->getEnd().get_double_mjd() - _session_ptr->getStart().get_double_mjd())*86400 / _dt) +1;
                        
            ivg::Date current_epoch = _session_ptr->getStart();
            bool previousEpochIsValid = false;
            int start_i = -1;
            
            std::vector<lps::Position> positions;
            // loop over all time steps
            for( unsigned int i = 0; i < n; ++i){
                
                // calculate azimuth and elevation in radian               
                // TODO which calc_az_el function should be used?
                
                ivg::Matrix azel = sta->calc_az_el(current_epoch, *sou);
                
                
//                ivg::Matrix crf2trf = _session_ptr->_eops.form_crf2trf( current_epoch  );
//                ivg::Matrix k = sou->get_unit_vector_ssb();
//                ivg::Matrix azel = sta->calc_az_el(current_epoch, k, crf2trf); 
                                            
                // calculate arc between sun and source in degree
                double sun_dist = sou->calc_arclength_to_sun( current_epoch )*ivg::rad2d;
                
                // Only add to trail if all condotions are satisfied
                // (1) elevation is larger than global elevation limit
                // (2) source is locally visible (elevation mask from configfile)
                // (3) distance to sun is large enough
                if( azel(1) > min_ele  && sta->check_visibility( azel(0), azel(1) ) && sun_dist > min_dist_2_sun ){
                    
                    if(!previousEpochIsValid){
                        start_i = i;
                    }
                    
                    // az/el is saved in degree in Position
                    azel*=ivg::rad2d;
                    positions.push_back(lps::Position( azel(0), azel(1), current_epoch));
         
                    previousEpochIsValid = true;
                    

                } else {
                    if( previousEpochIsValid ){
                        // finish trail after invalid observation
                        
                        // a trail must have at least 2 positions
                        if( positions.size() >=2 ){
                            lps::Transit tra (transit_counter, sta_idx, sou_idx, positions.begin()->date(), positions.back().date(), start_i, i-1, positions, _session_ptr->getStart(), _session_ptr  );
                            container.push_back( tra );
                            transit_counter++;
                            positions.clear();
                        }

                    }
                    previousEpochIsValid = false;
                    
                }
                      
                current_epoch.add_secs(_dt);
            }
            // finish trail if last observation is valid
            if( previousEpochIsValid ){
                // a trail must have at least 2 positions
                if( positions.size() >=2 ){
                    lps::Transit tra (transit_counter, sta_idx, sou_idx, positions.begin()->date(), positions.back().date(), start_i, n-1, positions, _session_ptr->getStart(), _session_ptr  );
                    container.push_back( tra );
                    transit_counter++;
                    positions.clear();
                }
            }
            
            // only create entries in mapTransits_ if a transit for the source station combination exists
            if(transit_counter > first  ){
                mapTransits_[ std::make_pair(sou_idx, sta_idx ) ] = std::make_pair(first, transit_counter );
            }
        }
    }
    
    for(lps::Transit& tra : container){
        tra.approxVelocity(_dt);
    }    
    
//    for(int j = 0; j < container.size(); ++j){
//        std::cout << container[j].id() << " " << container[j].station() << " " << container[j].quasar() << " --------------------- "<< std::endl;
//        for(lps::Position& pos : container[j].path()){
//            std::cout << pos.azimuth() << " " << pos.elevation() << std::endl;
//        }
//        
//    }
    
  
}


Transits::TransitIterator Transits::transit(int source_idx, int station_idx) const{
    auto it = mapTransits_.find({source_idx, station_idx});
    if(it == mapTransits_.end()){
        return make_const_range(container,0,0);
    }
    return make_const_range(container,it->second.first,it->second.second);
}

bool Transits::hasTransits(int source_idx, int station_idx) const{
    auto it = mapTransits_.find({source_idx, station_idx});
    return it != mapTransits_.end();
}


Transits Transits::seenByAtLeastBaselines( unsigned n_bl ) const
{
    
    ivg::Schedule s(_session_ptr);
    
    unsigned n_sta = ceil( 0.5 +sqrt(0.25 + 2.0*n_bl) );
    log<DETAIL>("*** compute transits seen by at least ") % n_bl % " baselines ( " % n_sta % " stations )" ;
    
    Setting* const setup = _session_ptr->get_setup();
    double min_scan_duration = (double)(*setup)["SKED"]["min_scan"];
    double max_scan_duration = (double)(*setup)["SKED"]["max_scan"];
    
    double x_thres = (double)(*_session_ptr->get_setup())["SKED"]["snr_min_x"]-(double)(*_session_ptr->get_setup())["SKED"]["snr_margin_x"];
    double s_thres = (double)(*_session_ptr->get_setup())["SKED"]["snr_min_s"]-(double)(*_session_ptr->get_setup())["SKED"]["snr_margin_s"];
    
    bool twinObserveTogether = false;
    if(  (*setup)["SKED"].exists("allow_twin_same_source") ){
        twinObserveTogether = (bool)(*_session_ptr->get_setup())["SKED"]["allow_twin_same_source"];
    }
        
     // create new Transits object copy session pointer ane compute observable transits from existing transits
    lps::Transits observableTransits;
    observableTransits._session_ptr = this->_session_ptr;
    observableTransits._dt = this->_dt;
       
    // number of time stemps
    unsigned int n = ceil( (_session_ptr->getEnd().get_double_mjd() - _session_ptr->getStart().get_double_mjd())*86400 / _dt)+1;
    
    std::vector<std::pair<std::string, std::string>> exclude_bl (0);
    if(  (*setup)["SKED"].exists("exclude_bl") ){
        Setting &ebl = (*setup)["SKED"]["exclude_bl"];
        exclude_bl.resize(ebl.getLength());

        for(int i=0; i < ebl.getLength(); i++){  
            std::string first = ebl[i][0];
            std::string second = ebl[i][1];
            exclude_bl[i] = make_pair( first, second  );
        }
    }
    
    _session_ptr->get_trf_ptr()->create_station_indices();
     
    // loop over all sources
    for(size_t sou_idx =0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx){
        ivg::Source* sou = _session_ptr->get_crf_ptr()->get_source(sou_idx);
                
        // check if at least n_sta stations have a transit involving the same source
        unsigned staionsWithTransit = 0;
        for(size_t sta_idx=0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx){
            if( this->hasTransits(sou_idx, sta_idx) == true ){
                staionsWithTransit++;
            }
        }
        
        if(staionsWithTransit < n_sta)
            continue;
        
        // check for each station in each interval the visibility
        std::vector<std::vector<bool>> visibilityMatrix(_session_ptr->get_trf_ptr()->get_number_stations());
        for(size_t sta_idx=0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx){
            visibilityMatrix[sta_idx] = std::vector<bool>(n, false);
            
            TransitIterator it = this->transit(sou_idx, sta_idx);
            for(const Transit& trans : it){
                int s = trans.get_start_idx(); // or round(trans.begin() / _dt);
                int e = trans.get_end_idx();    // or round(trans.end() / _dt);

                for( int i = s; i <= e; ++i ){                        
                    visibilityMatrix[sta_idx][i] = true;
                }
            }
        }
         
        // check baseline visibility
        for(size_t sta_idx_1=0; sta_idx_1 < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx_1){
            
            std::vector<short> nVisibleBl(n, 0); // number of baselines including sta_idx_1 seeing the same source simultaneously
            for(size_t sta_idx_2 = 0; sta_idx_2 < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx_2){          
                if( !(twinObserveTogether == false && _session_ptr->get_trf_ptr()->areTwins(sta_idx_1, sta_idx_2) == true) && sta_idx_1 != sta_idx_2 ){                   
                    std::string sta_name_1 = _session_ptr->get_trf_ptr()->get_station( sta_idx_1 )->get_name(ivg::staname::ivs_name);
                    std::string sta_name_2 = _session_ptr->get_trf_ptr()->get_station( sta_idx_2 )->get_name(ivg::staname::ivs_name);
                    bool exclude = false;
                    for(std::pair<std::string, std::string>& p : exclude_bl ){
                        if( (sta_name_1 == p.first && sta_name_2 == p.second) || (sta_name_2 == p.first && sta_name_1 == p.second) ){
                            exclude = true;
                            break;
                        }
                    }
                    if(exclude)
                        continue;
                    
                    for(unsigned i=0; i < n; ++i){
                        
                        if(visibilityMatrix[sta_idx_1][i] && visibilityMatrix[sta_idx_2][i]){
                            // compute SNR
                            ivg::Date d = _session_ptr->getStart();
                            d.add_secs(i*_dt);

                            ivg::BandInfo xband, sband;
                            s.compute_band_info (*_session_ptr->get_crf_ptr()->get_source( sou_idx ), *_session_ptr->get_trf_ptr()->get_station(sta_idx_1),
                                                 *_session_ptr->get_trf_ptr()->get_station(sta_idx_2), d, xband, sband );
                             double snr_x(0.0), snr_s(0.0);
                            s.compute_snr(max_scan_duration, xband, sband, snr_x, snr_s);
                            nVisibleBl[i] += (int)(snr_x > x_thres && snr_s > s_thres);
                        }
                    }
                }
            }

             // find start and end of transits 
            bool newTransit = false;
            int s = 0;
            int e = 0;
            std::vector<std::pair<int,int>> start_end;
            for (int i = 0; i < n; ++i) {
                if (nVisibleBl[i] >= n_bl) {
                    if (newTransit == false) {
                        s = i;
                        newTransit = true;
                    }
                } else {
                    if (newTransit) {
                        e = i - 1;
                        newTransit = false;
                        start_end.push_back(std::pair<int, int>(s, e));
                    }
                }

            }
            if( newTransit ){
                    start_end.push_back( std::pair<int,int>(s,n-1));
            }

            // add new transit to container
            TransitIterator it = this->transit(sou_idx, sta_idx_1);
            int first = observableTransits.container.size();
            for(const std::pair<int,int>& se : start_end){
                                
                for(const Transit& trans : it){
                    if( trans.get_start_idx() <= se.first && trans.get_end_idx() >= se.second ){
                        
                        std::vector<lps::Position> path; 
                                
                        for(int i = se.first - trans.get_start_idx(); i <= se.second - trans.get_start_idx(); ++i){
                            path.push_back( trans.path()[i] );
                        }
                        

                        ivg::Date start = trans.path()[se.first - trans.get_start_idx()].date();
                        ivg::Date end = trans.path()[se.second - trans.get_start_idx()].date();
                        
//                        ivg::Date start2 = _session_ptr->getStart();
//                        start2.add_secs( se.first*_dt);
//                        ivg::Date end2 = _session_ptr->getStart();
//                        end2.add_secs( se.second*_dt);
//                        assert( abs(start.get_double_mjd() - start2.get_double_mjd()) < 1e-3 );
//                        assert( abs(end.get_double_mjd() - end2.get_double_mjd()) < 1e-3 );
                        
                        
                        if( (se.second-se.first)*_dt  >= min_scan_duration ){

                            observableTransits.container.push_back( Transit(trans.id(),trans.get_sta_idx(),trans.get_scr_idx(),
                                                                   start, end, se.first, se.second, path, _session_ptr->getStart(), _session_ptr) );
                        }
                                                       
                        
                        break;
                    }
                }
            }
            if(observableTransits.container.size() > first  ){
                observableTransits.mapTransits_[std::make_pair(sou_idx, sta_idx_1)] = std::make_pair(first, observableTransits.container.size());
            }
        }

    }
    return observableTransits;
}


Transits Transits::seenByAtLeastStations( unsigned n_sta ) const
{
    log<DETAIL>("*** compute transits seen by at least ") % n_sta % " stations";
    int counter=0;
    
    Setting* const setup = _session_ptr->get_setup();
    double min_scan_duration = (double)(*setup)["SKED"]["min_scan"];
    double max_scan_duration = (double)(*setup)["SKED"]["max_scan"];
        
     // create new Transits object copy session pointer ane compute observable transits from existing transits
    lps::Transits observableTransits;
    observableTransits._session_ptr = this->_session_ptr;
    observableTransits._dt = this->_dt;
       
    // number of time stemps
    unsigned int n = ceil( (_session_ptr->getEnd().get_double_mjd() - _session_ptr->getStart().get_double_mjd())*86400 / _dt)+1;     
     
    // loop over all sources
    for(size_t sou_idx =0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx){
        ivg::Source* sou = _session_ptr->get_crf_ptr()->get_source(sou_idx);
                
        // check if at least n_sta stations have a transit involving the same source
        unsigned staionsWithTransit = 0;
        for(size_t sta_idx=0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx){
            if( this->hasTransits(sou_idx, sta_idx) == true ){
                staionsWithTransit++;
            }
        }
        
        if(staionsWithTransit < n_sta)
            continue;
        
        // compute for each time stemp the number of stations seeing the current source
        
        std::vector<short> nVisibleSta(n, 0); // number of stations seeing the same source simultaneously
        
        _session_ptr->get_trf_ptr()->create_station_indices();
        // TODO allow_twin_same_source option
        for( std::vector<ivg::Analysis_station*>& stas : _session_ptr->get_trf_ptr()->get_station_twin_list() ){
            std::vector<bool> observatorySees(n, false); // true if at least one station at the obsrvatory sees the source
            for(ivg::Analysis_station* sta: stas){
                int sta_idx = sta->get_idx();
                TransitIterator it = this->transit(sou_idx, sta_idx);
                for(const Transit& trans : it){
                    int s = trans.get_start_idx(); // or round(trans.begin() / _dt);
                    int e = trans.get_end_idx();    // or round(trans.end() / _dt);

                    for( int i = s; i <= e; ++i ){                        
                        observatorySees[i] = true;
                    }
                }
            }
            for(int i = 0; i<n; ++i){
                nVisibleSta[i] += (int)observatorySees[i];
            }
        }
        

//        for(size_t sta_idx=0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx){
//            //ivg::Analysis_station* sta = _session_ptr->get_trf_ptr()->get_station(sta_idx);
//                        
//            TransitIterator it = this->transit(sou_idx, sta_idx);
//            for(const Transit& trans : it){
//                int s = trans.get_start_idx(); // or round(trans.begin() / _dt);
//                int e = trans.get_end_idx();    // or round(trans.end() / _dt);
//                                
//                for( int i = s; i <= e; ++i ){
//                    //std::cerr << i << " ";
//                    nVisibleSta[i]++;
//                }
//            }
//
//        }
                
        // find start and end of transits 
        bool newTransit = false;
        int s = 0;
        int e = 0;
        std::vector<std::pair<int,int>> start_end;
        for (int i = 0; i < n; ++i) {
            if (nVisibleSta[i] >= n_sta) {
                if (newTransit == false) {
                    s = i;
                    newTransit = true;
                }
            } else {
                if (newTransit) {
                    e = i - 1;
                    newTransit = false;
                    start_end.push_back(std::pair<int, int>(s, e));
                }
            }
            
        }
        if( newTransit ){
                start_end.push_back( std::pair<int,int>(s,n-1));
        }
        
        // add new transit to container
              
        
        for(size_t sta_idx=0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx){
            TransitIterator it = this->transit(sou_idx, sta_idx);
            int first = observableTransits.container.size();
            for(const std::pair<int,int>& se : start_end){
                                
                for(const Transit& trans : it){
                    if( trans.get_start_idx() <= se.first && trans.get_end_idx() >= se.second ){
                        
                        std::vector<lps::Position> path; 
                                
                        for(int i = se.first - trans.get_start_idx(); i <= se.second - trans.get_start_idx(); ++i){
                            path.push_back( trans.path()[i] );
                        }
                        

                        ivg::Date start = trans.path()[se.first - trans.get_start_idx()].date();
                        ivg::Date end = trans.path()[se.second - trans.get_start_idx()].date();
                        
//                        ivg::Date start2 = _session_ptr->getStart();
//                        start2.add_secs( se.first*_dt);
//                        ivg::Date end2 = _session_ptr->getStart();
//                        end2.add_secs( se.second*_dt);
//                        assert( abs(start.get_double_mjd() - start2.get_double_mjd()) < 1e-3 );
//                        assert( abs(end.get_double_mjd() - end2.get_double_mjd()) < 1e-3 );
                        
                        
                        if( (se.second-se.first)*_dt  >= min_scan_duration ){

                            observableTransits.container.push_back( Transit(trans.id(),trans.get_sta_idx(),trans.get_scr_idx(),
                                                                   start, end, se.first, se.second, path, _session_ptr->getStart(), _session_ptr) );
                        }
                                                       
                        
                        break;
                    }
                }
            }
            if(observableTransits.container.size() > first  ){
                observableTransits.mapTransits_[std::make_pair(sou_idx, sta_idx)] = std::make_pair(first, observableTransits.container.size());
            }
        } 
    }
    return observableTransits;
}

Transits Transits::getBaselineWiseTransits( ) const
{
    
    ivg::Schedule s(_session_ptr);
    
    log<DETAIL>("*** compute observable transits");
    int counter=0;  
    
    Setting *setup = _session_ptr->get_setup();
    double min_scan_duration = (double)(*setup)["SKED"]["min_scan"];
    double max_scan_duration = (double)(*setup)["SKED"]["max_scan"];
        
    // create new Transits object copy session pointer ane compute observable transits from existing transits
    lps::Transits observableTransits;
    observableTransits._session_ptr = this->_session_ptr;
    observableTransits._dt = this->_dt;
    
    // loop over all sources
    for(size_t sou_idx =0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx){
        ivg::Source* sou = _session_ptr->get_crf_ptr()->get_source(sou_idx);
        bool observable=false;
        
        // loop over all stations
        for(size_t sta_idx_1=0; sta_idx_1 < _session_ptr->get_trf_ptr()->get_number_stations() ;++sta_idx_1){
            ivg::Analysis_station* sta_1 = _session_ptr->get_trf_ptr()->get_station(sta_idx_1);
            
            size_t first =  observableTransits.container.size();
            
            // loop over all transits corresponding to current source/station combination
            TransitIterator it1 = this->transit(sou_idx, sta_idx_1);
            for(const Transit& trans1 : it1){
                                               
                // loop over all stations excluding the current station from superior loop
                for(size_t sta_idx_2=0; sta_idx_2 < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx_2){
                    if(sta_idx_1==sta_idx_2){continue;}
                    ivg::Analysis_station* sta_2 = _session_ptr->get_trf_ptr()->get_station(sta_idx_2);
                                                            
                    // loop over all transits corresponding to current source/station2 combination
                    TransitIterator it2 = this->transit(sou_idx, sta_idx_2);
                    for(const Transit& trans2 : it2){
                        
                        ivg::Date left, right;
                        // left is latest epoch of both transits
                        if( trans1.get_startEpoch() > trans2.get_startEpoch() ){
                            left = trans1.get_startEpoch();
                        } else {
                            left = trans2.get_startEpoch();
                        }
                        
                        // right is former epoch of both transits
                        if( trans1.get_endEpoch() < trans2.get_endEpoch() ){
                            right = trans1.get_endEpoch();
                        } else {
                            right = trans2.get_endEpoch();
                        }
                        
                        if( (right.get_double_mjd() - left.get_double_mjd())*DurationDay >= min_scan_duration ){
                            
                            observable = true;
                            
                            int start_idx = trans1.get_start_idx();
                             
                            std::vector<lps::Position> path;
                            path.reserve( (right.get_double_mjd() - left.get_double_mjd()) * DurationDay / _dt +1 );

                            
                            unsigned k = 0;
                            for(const lps::Position& pos: trans1.path() ){
                                // > operator is t1 -t2 > 1ms   -> t1==t2    
                                if( pos.date() >= left && pos.date() <= right){
                                    path.push_back(pos);
                                    
                                    ivg::Date d = pos.date();

                                    ivg::BandInfo xband, sband;
                                    s.compute_band_info (*_session_ptr->get_crf_ptr()->get_source( sou_idx ), *_session_ptr->get_trf_ptr()->get_station(sta_idx_1),
                                                         *_session_ptr->get_trf_ptr()->get_station(sta_idx_2), d, xband, sband );
                                    double duration = s.compute_obs_duration(*_session_ptr->get_crf_ptr()->get_source( sou_idx ), *_session_ptr->get_trf_ptr()->get_station(sta_idx_1),
                                                         *_session_ptr->get_trf_ptr()->get_station(sta_idx_2),d, xband, sband, false );  
                                    
                                    path.back().set_minObsDur(duration);
                                    
                                    double snrx(0.0), snrs(0.0);
                                    s.compute_snr( max_scan_duration, xband, sband, snrx, snrs );  
                                    
                                    path.back().set_SNR(snrx, snrs);
                                    
                                    
                                    k++;
                                }
                                
                                if(k == 0 ){
                                    start_idx++;
                                }
                               
                            }
                            int end_idx = start_idx + path.size() -1;
                             
                          
                            observableTransits.container.push_back(Transit(trans1.id(),trans1.get_sta_idx(),trans1.get_scr_idx(),
                                                               left, right, start_idx, end_idx, path,
                                                               _session_ptr->getStart(), _session_ptr, sta_idx_2));

                        }
                    }
                }
            }
            if(observableTransits.container.size() > first  ){
                observableTransits.mapTransits_[std::make_pair(sou_idx, sta_idx_1)] = std::make_pair(first, observableTransits.container.size());   
            }
        }
        if(observable){
            counter++;
        }
    }

    log<RESULT>("*** number of observable quasars ") % counter;
    return observableTransits;
}

bool Transits::get_border_pos(const int sou_idx, const int sta_idx, const lps::Seconds begin, 
        const lps::Seconds end, const lps::Seconds intLen, std::pair<lps::Position, lps::Position>& pos ) const{
    
    TransitIterator it = this->transit(sou_idx, sta_idx);
    for(const Transit& trans : it){
        unsigned int transit_start_idx = trans.get_start_idx(); 
        unsigned int transit_end_idx = trans.get_end_idx();
                
        bool t_starts_in_obsInt = transit_start_idx*_dt >= begin && transit_start_idx*_dt <= begin + intLen;
        bool t_ends_in_endInt = transit_end_idx*_dt <= end   &&   transit_end_idx*_dt >=   end - intLen;
        
        bool t_starts_before_obsInt = transit_start_idx*_dt <= begin;
        bool t_ends_after_endInt = transit_end_idx*_dt >= end;

        
        // case: transit is overlapping entire obs interval
        if( t_starts_before_obsInt && t_ends_after_endInt){
            pos.first  = trans.path().at( floor( begin/(double)_dt ) - transit_start_idx );
            pos.second  = trans.path().at( ceil( end/(double)_dt ) - transit_start_idx );
            return true;
        }
        // case: transit starts in obs interval and ends after obs end
        if( t_starts_in_obsInt && t_ends_after_endInt){
            pos.first  = trans.path().at( 0 );
            pos.second  = trans.path().at( ceil( end/(double)_dt ) - transit_start_idx );
            return true;
        }
        // case: transit starts befor obs interval and ends in the end interval
        if( t_starts_before_obsInt && t_ends_in_endInt ){
            pos.first  = trans.path().at( floor( begin/(double)_dt ) - transit_start_idx );
            pos.second  = trans.path().back();
            return true;
        }
        if( t_starts_before_obsInt && t_ends_in_endInt ){
            pos.first  = trans.path().at( 0 );
            pos.second  = trans.path().back();
            return true;
        }
    }
    
    return false;
}

bool Transits::sourceIsVisible( const int sou_idx, const int sta_idx, const lps::Seconds t ) const{
    TransitIterator it = this->transit(sou_idx, sta_idx);
    for(const Transit& trans : it){
        if( t >= trans.begin() && t <= trans.end() ){
            return true;
        }
    }
    return false;
}

void Transits::compute_possible_coverage( std::map<int, std::vector<double> >& max_possible_coverage, std::map<int, std::vector<double> >& max_possible_surface_area, int max_level ) const{
    
    int max_tree_level = 0;
    if (max_level == 0){
        max_tree_level = (int)(*_session_ptr->get_setup())["SKED"]["max_tree_level"];
    } else {
        max_tree_level = max_level;
    }
    Seconds intervalLength = (*_session_ptr->get_setup())["SKED"]["min_scan"];
    int n  = (_session_ptr->getEnd().get_double_mjd() - _session_ptr->getStart().get_double_mjd())*lps::DurationDay / intervalLength;
    
    max_possible_coverage.clear();
    max_possible_surface_area.clear();
    
    for(size_t sta_idx =0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx){
        std::vector<lps::DataPoint<bool>> points;
        max_possible_coverage[sta_idx] = std::vector<double>(max_tree_level, 0.0);
        max_possible_surface_area[sta_idx] = std::vector<double>(max_tree_level, 0.0);
        for(size_t sou_idx =0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx){
//            TransitIterator it = this->transit(sou_idx, sta_idx);
//            for(const Transit& trans : it){
//                for( const lps::Position& pos : trans.path()){
//                    lps::Point p (pos.azimuth(), pos.elevation());
//                    points.push_back( lps::DataPoint<bool>(p, true) );
//                }
//            }
            for(int i = 0; i < n; ++i){
                std::pair<lps::Position, lps::Position> pos;
                if( this->get_border_pos(sou_idx, sta_idx, i*intervalLength, (i+1)*intervalLength, intervalLength, pos ) ){
                    lps::Point p (pos.first.azimuth(), pos.first.elevation());
                    points.push_back(lps::DataPoint<bool>(p, true) );
                }
            }
        }
 
        // level zero is full circle, point is in the center
        lps::Wedge wedge (lps::Point(0,0),0,90,0,360);
        
        lps::Node<bool, lps::Wedge>  root = lps::Node<bool, lps::Wedge>(points, wedge, 0 , max_tree_level);
        
        
        root.visitBoundAndBool([&](int level,  const lps::Wedge& w, bool valid){         
            if( valid ){
                max_possible_coverage[sta_idx][level-1] += 1.0;
                max_possible_surface_area[sta_idx][level-1] += w.surfaceArea();
            } 
        });
    }
        
}

Transit::Transit(int id, int sta_idx, int sou_idx, const ivg::Date& startEpoch, const ivg::Date& endEpoch,
    int start_idx, int end_idx, const std::vector<Position>& path,ivg::Date referenceEpoch, ivg::Session* session){
    if(endEpoch < startEpoch){
        log<WARNING>("!!! end time of transit before start ");
        exit(EXIT_FAILURE);
    }
    
    this->id_ = id;
    this->sta_idx = sta_idx;
    this->scr_idx = sou_idx;
    this->_startEpoch = startEpoch;
    this->_endEpoch = endEpoch;
    this->_start_idx = start_idx;
    this->_end_idx = end_idx;
    this->path_ = path;
    this->_referenceEpoch = referenceEpoch;
    this->_source_name =  session->get_crf_ptr()->get_source(scr_idx)->get_name(ivg::srcname::ivs);
}

Transit::Transit(int id, int sta_idx, int sou_idx, const ivg::Date& startEpoch, const ivg::Date& endEpoch, int start_idx, int end_idx,
            const std::vector<Position>& path, ivg::Date referenceEpoch, ivg::Session* session, int visbibleFromStation) : 
            Transit(id, sta_idx, sou_idx, startEpoch, endEpoch, start_idx, end_idx, path, referenceEpoch, session){
    
    this->visbibleFromStation = visbibleFromStation;
}

std::vector<Position> Transit::get_path_part(ivg::Date begin, ivg::Date end) const{

    std::vector<Position> path;
    for( const lps::Position&  p : this->path_){
        if( p.date() >= begin && p.date() < end ){
            path.push_back(p);
        }
    }
    return path;
}

std::string Transit::get_sou_name() const{
    return _source_name;
}

void Transit::approxVelocity(double dt){
    
    unsigned n = path_.size();
    if(n >= 3){
        for(unsigned i = 1; i < n-1; ++i){
            path_[i].set_dazi(  azimuth_diff(path_[i+1].azimuth(), path_[i-1].azimuth())    / (2*dt)  );
            path_[i].set_dele(  (path_[i+1].elevation() - path_[i-1].elevation())  / (2*dt)  );
        }
    }
    
    path_[0].set_dazi(  azimuth_diff(path_[1].azimuth(), path_[0].azimuth())    / dt  );
    path_[0].set_dele(  (path_[1].elevation() - path_[0].elevation())  / dt  );
    
    path_[n-1].set_dazi(  azimuth_diff(path_[n-1].azimuth(), path_[n-2].azimuth())    / dt  );
    path_[n-1].set_dele(  (path_[n-1].elevation() - path_[n-2].elevation())  / dt  );
    
}

void Transits::print_transits() const{
    
    ivg::Trf* trf = _session_ptr->get_trf_ptr();
    ivg::Crf* crf = _session_ptr->get_crf_ptr();
   std::cout << std::string(145,'-') << std::endl; 
   std::cout << "number of transits: " << container.size() << std::endl; 
   std::cout << "station 1     station 2       source              time [s]       index       az[deg]     el[deg]   obs time[s]           snr x           snr s" << std::endl;
   std::cout << "                                               start-end      start-end    start-end   start-end     start-end       start-end       start-end" << std::endl;
   std::cout << std::string(145,'-') << std::endl; 
    
    for( int tr_ind = 0; tr_ind < container.size(); ++tr_ind ){
        const Transit& tr = container[tr_ind];
        
        std::string sta_name = trf->get_station( tr.get_sta_idx() )->get_name(ivg::staname::ivs_name);
        std::string src_name = crf->get_source(  tr.get_scr_idx()  )->get_name(ivg::srcname::ivs);
        
        std::string visbleFrom;
        
        if( tr.get_visbibleFromSta_idx() >= 0 ){
            visbleFrom = trf->get_station( tr.get_visbibleFromSta_idx() )->get_name(ivg::staname::ivs_name);
        } else {
            visbleFrom = "n/a";
        }
        
        std::cout << setw(8)  << sta_name << " (" << setw(2) << tr.get_sta_idx() << ") "
                  << setw(8)  << visbleFrom << " (" << setw(2) << tr.get_visbibleFromSta_idx() << ") "
                  << setw(10) << src_name << " (" << setw(4) << tr.get_scr_idx() << ") " 
                  << setw(5) << tr.begin() << " - " << setw(5) << tr.end() << "   "
                  << setw(4) << tr.get_start_idx() << " - " << setw(4) << tr.get_end_idx() << "   "
                  << setw(3) << round(tr.path()[0].azimuth()) << " - " <<  setw(3) << round(tr.path().back().azimuth()) << "   "  
                  << setw(3) << round(tr.path()[0].elevation()) << " - " <<  setw(3) << round(tr.path().back().elevation()) << "   "
                  << setw(5) << round(tr.path()[0].minObsDur()) << " - " <<  setw(5) << round(tr.path().back().minObsDur()) << "   "
                  << setw(5) << round(tr.path()[0]. snrx()) << " - " <<  setw(5) << round(tr.path().back().snrx()) << "   "
                  << setw(5) << round(tr.path()[0]. snrs()) << " - " <<  setw(5) << round(tr.path().back().snrs()) << "   "
                  <<   std::endl;
    }
}


void Transits::save_transits(std::string path) const {

    double sesStartmjd = _session_ptr->getStart().get_double_mjd();
    
    // loop over all stations
    for (size_t sta_idx_1 = 0; sta_idx_1 < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx_1) {
        ivg::Analysis_station* sta_1 = _session_ptr->get_trf_ptr()->get_station(sta_idx_1);

        std::ofstream ofs; // write only
        std::string output = path + '/' + _session_ptr->get_name() + '_' + sta_1->get_name(ivg::staname::ivs_name) + ".dat";
        ofs.open(output, ios::out);
        if (!ofs.good()) {
            log<WARNING>("!!! can not write transit to ") % output; 
        } else {
            log<INFO>("*** write transits to: ") % output; 
            ofs << "source  epoch(sec) az(deg) el(deg) Tmin(sec) snrx snrs vaz(deg/s) vel(deg/s)" ;
            // loop over all sources
            for (size_t sou_idx = 0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx) {
                ivg::Source* sou = _session_ptr->get_crf_ptr()->get_source(sou_idx);
                std::string souname = sou->get_name(ivg::srcname::ivs);

                // loop over all transits corresponding to current source/station combination
                TransitIterator it = this->transit(sou_idx, sta_idx_1);
                for (const Transit& trans : it) {
                    for(const lps::Position& pos: trans.path()){
                        ofs << setw(8) << souname << " " << setw(5) << setprecision(0) << std::fixed << round( (pos.date().get_double_mjd() - sesStartmjd)*86400) 
                            << " " << setw(13) << setprecision(9) << pos.azimuth() << " " << setw(13) << setprecision(9) << pos.elevation()
                            << " " << setw(6) << setprecision(0) << std::fixed <<  pos.minObsDur()
                            << " " << setw(7) << setprecision(3) << pos.snrx() << " " << setw(7) << setprecision(3) << pos.snrs()
                            << " " << setw(13) << setprecision(9) << pos.dazi() << " " << setw(13) << setprecision(9) << pos.dele() << std::endl;
                    }
                }
            }
        }
        ofs.close();
    }

}

void  Transits::remove_source( int source_idx ){
        
    for (size_t sta_idx = 0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations(); ++sta_idx) {
        
        auto it = mapTransits_.find({source_idx, sta_idx}); 
        
        // key exists in map mapTransits_
        if( it != mapTransits_.end() ){
            
            int first = it->second.first;
            int numDeleted = it->second.second - it->second.first;
            container.erase (container.begin()+it->second.first, container.begin()+it->second.second);
            
            // update map
            mapTransits_.erase(it);
            
            for(auto& a : mapTransits_ ){
                if(a.second.first > first){
                    a.second.first -= numDeleted;
                    a.second.second -= numDeleted;
                }
            }
        }
    }
    
}

Position::Position(double azimuth, double elevation, const ivg::Date& date, lps::Seconds minObsDur){
    this->azimuth_ = azimuth; // in degree
    this->elevation_ = elevation; // in degree
    this->date_ = date;
    this->minObsDur_ = minObsDur;
}


} //namespace

