#include "schedule.h"
#include "station.h"

namespace ivg
{

// ...........................................................................  
Schedule::Schedule( )
// ...........................................................................
{
    
}
// ...........................................................................   
Schedule::Schedule( ivg::Session *session_ptr ) : Schedule()
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ Schedule::Schedule( ivg::Session *session_ptr )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
_session_ptr = session_ptr;
Setting& setup = *_session_ptr->_setup;
    
// recorder name syntax: tracks-channels-fanout-bits
string rec_name = setup["SKED"]["rec_name"]; // e.g. "00-16-0-1"
_bit_sampling = std::stoi(rec_name.substr(rec_name.size()-1,1)); // last number of rec_name

if( setup["SKED"].exists("broadband") && setup["SKED"]["broadband"]["apply"] )
    _bandwidth = ( (double)setup["SKED"]["broadband"]["bandwidth"])*1.0e6; // in Hz
else
    _bandwidth = ( (double)setup["SKED"]["bandwidth"])*1.0e6; // in Hz

if( (int)(*_session_ptr->_setup)["SKED"]["approach"] == 1){

    Setting &ea = setup["SKED"]["EqualAreaGrid"];

    _ea_grid_setup.resize(ea.getLength());

    for(int i=0; i<ea.getLength(); i++){
        for(int j=0; j<ea[i].getLength(); j++){
            _ea_grid_setup[i].push_back( (unsigned)ea[i][j] );
        }
    }

    lps::Seconds duration = (_session_ptr->getEnd().get_double_mjd()-_session_ptr->getStart().get_double_mjd())*3600*24;
    lps::Seconds temporal_grid_resolution = (int)setup["SKED"]["temporal_resolution"] * 60;
    lps::Seconds temporal_grid_shift = (int)setup["SKED"]["temporal_shift"] * 60;
    _tg = lps::TemporalGrid(duration, temporal_grid_resolution, temporal_grid_shift, _session_ptr->getStart());
    
}

    
    
#ifdef DEBUG_VLBI
   cerr << "--- Schedule::Schedule( ivg::Session *session_ptr )" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
void Schedule::start_scheduling(){
// ...........................................................................
#ifdef DEBUG_VLBI
cerr << "+++ void Schedule::start_scheduling()" << endl; 
tictoc tim;
tim.tic();
#endif 
    
    if( (int)(*_session_ptr->_setup)["SKED"]["approach"] >= 4){ // use global approach
                
        _session_ptr->_schedule_ptr = this;
        
           
        // remove sources that are not visible
        lps::Transits coarseTransits(_session_ptr, 120);
        _session_ptr->_trasintsVisibleFromStation = coarseTransits;
        
//        lps::Transits infoTr = coarseTransits.getBaselineWiseTransits();
//        infoTr.print_transits();
//        infoTr.save_transits((std::string)(const char*)(*_session_ptr->_setup)[ "outdir" ] + "/");
        
        std::vector<unsigned> visble_from_stas = number_of_stations_source_is_visible(coarseTransits);
                        
        int erase_counter = 0;
        for(size_t sou_idx =0; sou_idx < visble_from_stas.size(); ++sou_idx){
            if (visble_from_stas[sou_idx] <2){
                _session_ptr->_crf.remove_source(sou_idx-erase_counter);
                erase_counter++;
            }
        }
        
        log<INFO>("*** removed ") % erase_counter % " of " % (visble_from_stas.size()) % " sources because they are not visible from at least two stations";
        
        // remove sources with bad SNR
        coarseTransits = lps::Transits(_session_ptr, 150);
        _session_ptr->_trasintsVisibleFrom2Stations = coarseTransits; 
        
        std::vector<bool> has_good_snr = quasar_has_at_least_one_bl_with_good_SNR(coarseTransits);

        erase_counter = 0;
        for(size_t sou_idx =0; sou_idx < has_good_snr.size(); ++sou_idx){
            if (has_good_snr[sou_idx] == false){
                _session_ptr->_crf.remove_source(sou_idx-erase_counter);
                erase_counter++;
            }
        }
        
        log<INFO>("*** removed ") % erase_counter % " of " % (has_good_snr.size()) % " sources because the minimum SNR constraints can not be achieved";
                

        _session_ptr->_complete_transits = lps::Transits (_session_ptr);
        _session_ptr->_common_transits = _session_ptr->_complete_transits.seenByAtLeastBaselines( 1 ); // getBaselineWiseTransits(); // 
        
//        _session_ptr->get_common_transits().getBaselineWiseTransits().print_transits();
                
        _session_ptr->create_crf_trf_inidices();
        
//        for(size_t sta_idx=0; sta_idx < _session_ptr->get_trf_ptr()->get_number_stations() ;++sta_idx){
//            for(size_t scr_idx=0; scr_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++scr_idx){
//                int counter = 0;
//                for(const lps::Transit& transit : _session_ptr->get_common_transits().transit(scr_idx, sta_idx)){
//                    std::cout << sta_idx << " " << scr_idx << " " << transit.get_visbibleFromSta_idx() << " " << transit.duration() << std::endl;
//                }
//            }
//        }

             
        
    } else { // use sequential approach
        
        // we need a pseudo-neutralpoint-scan as origin for the following scans
        ivg::Scan last_scan = create_neutralpoint_scan(_session_ptr->_start);

        int cnt = 0;
        while(true)
        {
            // get potential candidates for next scan based on the last_scan
            vector<ivg::Scan> potential_scans = get_next_potential_scans(last_scan);

            show_scans(potential_scans, "potential_scans1");

            // remove a potential scan if source within scan has been observed within the last XXX seconds
            remove_min_time_src(potential_scans, (double)(*_session_ptr->_setup)["SKED"]["min_time_src"]);
            if(potential_scans.size() == 0)
                throw runtime_error("ERROR: No potential scans left for schedulung.");

            // based on all left candidates, we need to determine the optimal scan, using different philosophies
            ivg::Scan next_scan;
            if(cnt < (int)(*_session_ptr->_setup)["SKED"]["initial_scans"] )
                next_scan = determine_optimal_scan((ivg::skedtype) (int)(*_session_ptr->_setup)["SKED"]["initial_approach"], potential_scans);
            else
                next_scan = determine_optimal_scan((ivg::skedtype) (int)(*_session_ptr->_setup)["SKED"]["approach"], potential_scans);

            // if next_scan is later than timerange for scheduling, stop scheduling
            double scan_end_mjd = next_scan.get_epoch().get_double_mjd() + next_scan.get_scheduled_duration()/86400.0;
            if(scan_end_mjd > _session_ptr->_end.get_double_mjd())
                break;

            // store selected scan in member-variable -> always everywhere available within schedule-class
            _scheduled_scans.push_back(next_scan);
            _session_ptr->_scans = _scheduled_scans; // update _scans within session
            _session_ptr->_nobs += next_scan.get_nobs(); // increase number of observations
            _session_ptr->_nobs_orig += next_scan.get_nobs(); // increase number of observations

            // great function to show a vector of scans
            show_scans(_scheduled_scans, "_schedulded_scans");

            last_scan = next_scan;        

            cnt++;
        }

        show_scans(_scheduled_scans, "FINAL SCHEDULE");

        // determine final end of session
        ivg::Date end = _scheduled_scans.back().get_epoch();
        end.add_secs(_scheduled_scans.back().get_scheduled_duration());
        _session_ptr->_end = end;

        _session_ptr->find_and_mark_unused_sources();

        // initialize _param_list within session for completeness
        _session_ptr->init_param_list();
        
        show_scans(_session_ptr->_scans, "FINAL SCHEDULE");
        
    }
#ifdef DEBUG_VLBI
cerr << "--- void Schedule::start_scheduling()" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................   
vector<ivg::Scan> Schedule::get_next_potential_scans(ivg::Scan last_scan)
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ vector<ivg::Scan> Schedule::get_next_potential_scans(ivg::Scan last_scan)" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
    vector<ivg::Scan> potential_scans; // this will be filled and finally returned
    
    ivg::Simulation simulator; // we need this for some privat simulation functions
    Setting *setup = _session_ptr->_setup;
    
    // we need this information for the scan meterology and simulation later on
    string ext_data_type = (*setup)["troposphere"]["external_meteo_data"][1];
    string ext_met_data = (*setup)["troposphere"]["external_meteo_data"][2];
    string gpt2filename = (*setup)["troposphere"]["gpt2_grid_file"];
    
    
    // create all possible baselines based on all stations within trf
    multimap<ivg::Analysis_station *, ivg::Analysis_station *> baselines;
    for(int i=0; i< _session_ptr->_trf.get_number_stations(); i++)
    {
        ivg::Analysis_station *sta_i = _session_ptr->_trf.get_station(i);
        for(int j=i+1; j< _session_ptr->_trf.get_number_stations(); j++)
        {
            ivg::Analysis_station *sta_j = _session_ptr->_trf.get_station(j);
            baselines.insert ( std::pair<ivg::Analysis_station *,ivg::Analysis_station *>(sta_i,sta_j) );
        }   
    }
    // get current crf2trf and partials matrix for each scan
    ivg::Partials_t2c tmp { ivg::Matrix( 3,3,0.0 ) };
    ivg::Partials_t2c * deriv_ptr = &tmp;
    ivg::Matrix crf2trf = _session_ptr->_eops.form_crf2trf( last_scan.get_epoch(), true, deriv_ptr );
    
    // to save some information about the reasons why a source is not usable
    map<string, int> reason_cnt;
    // go through all sources in crf (e.g. 301 good geodetic sources)
    // try to be close to ./jsked -f control/int2/2h/2014/k14103.cnf
    for(auto &src: _session_ptr->_crf)
    {
        // calculate possible start epoch for next scan
        // slewing time is calculated and added a few lines below within the loop over all sources
        ivg::Date epoch_new_scan = last_scan.get_epoch();
        epoch_new_scan.add_secs(last_scan.get_scheduled_duration()); // [sec]

        // for each candidate-source we need to calculate the potential slew time
        // Tslew = time required to slew to the next source
        string infostr;
        double slew_time = last_scan.calc_slew_time(&src, &_session_ptr->_eops, infostr);
        
        // in case of first scan we skip this part
        if(last_scan.get_source() != NULL)
        {
            // here we might need to consider some more const values like it's done in sked
            // Tproc = MIDTAPE + CHANGE + PREPASS + PARITY + spin + SETUP + HEAD
            // Tmax  = MAX (Tslew, Tproc + MAX(EARLY-CAL,0), MINIMUM)
            epoch_new_scan.add_secs(max(slew_time, (double)(*setup)["SKED"]["const_setup"]));
            
            // T(a2) = T(a1) + DUR + TAPE + IDLE + SOURCE + Tmax + CAL + TAPE
            epoch_new_scan.add_secs((double)(*setup)["SKED"]["const_idle"]);
            epoch_new_scan.add_secs((double)(*setup)["SKED"]["const_source"]);
            epoch_new_scan.add_secs((double)(*setup)["SKED"]["const_calib"]); // constant after slewing and before observation [sec]
        }
        
        double sun_dist = src.calc_arclength_to_sun( epoch_new_scan )*(180/M_PI);
        // only proceed if distance between quasar and sun is big enough (e.g. > 15.0 degree)
        if(sun_dist > (double)(*setup)["SKED"]["min_sun_dist"])
        {
            // based on each source, we create a scan, filling the scan with possible observations
            ivg::Scan scan(epoch_new_scan, &src, crf2trf.transpose(), deriv_ptr);
            
            // go through all baselines based on two Analysis_stations pointer
            // bl.first = Station 1
            // bl.second = Station 2
            int bl_cnt = 0;
            // in order to gain some scheduling-performance
            vector<string> ignore_station;
            for(auto &bl: baselines)
            {
                string sta1 = bl.first->get_name(ivs_name);
                string sta2 = bl.second->get_name(ivs_name);
                
                if(find( ignore_station.begin(), ignore_station.end(), sta1) != ignore_station.end() ||
                   find( ignore_station.begin(), ignore_station.end(), sta2) != ignore_station.end() )  
                {
                    // skip iteration in case of station with not useable elevation to source (e.g. negative)
                    continue;
                }

                string bl_name = sta1+"-"+sta2;
                
                // for each baseline we can create a observation
                int sta1_idx = scan.add_sta( bl.first );
                int sta2_idx = scan.add_sta( bl.second );
                
                ivg::Obs obs_new( _session_ptr, &scan, sta1_idx, sta2_idx );
                
                ivg::Matrix k = src.get_unit_vector_ssb();
                ivg::Matrix azel_sta1 = bl.first->calc_az_el(epoch_new_scan, k, crf2trf ); // [rad]
                ivg::Matrix azel_sta2 = bl.second->calc_az_el(epoch_new_scan, k, crf2trf ); // [rad]
                
                double min_ele = (double)(*setup)["SKED"]["min_elevation"]*(M_PI/180.0);
                
                // store stations with unusable elevation; not suitable for any baseline
                if(azel_sta1(1) < min_ele )
                    ignore_station.push_back(sta1);
                if(azel_sta2(1) < min_ele )
                    ignore_station.push_back(sta2);
                
                // here we need to check if quasar is aimable and above min. elevation mask from configfile
                if( azel_sta1(1) > min_ele && azel_sta2(1) > min_ele &&
                    bl.first->check_visibility(azel_sta1(0),azel_sta1(1)) && bl.second->check_visibility(azel_sta2(0),azel_sta2(1)) )
                {   
                      
                    ivg::BandInfo xband, sband;
                    compute_band_info (src, *bl.first, *bl.second, epoch_new_scan, xband, sband );
                    double duration = compute_obs_duration(src, *bl.first, *bl.second, epoch_new_scan, xband, sband );
                                 
                    double max_scan = (*setup)["SKED"]["max_scan"];
                    
                    
                    // only proceed if duration is less than max_scan from configfile
                    if(duration <= max_scan)
                    {
                        // stored information about duration in seconds within each scan
                        scan.set_schedulded_duration(duration);
                        
                        // check if potential source is still visible after final duration => end of scan
                        ivg::Date scan_end = epoch_new_scan;
                        scan_end.add_secs(duration);
                        ivg::Matrix crf2trf_end = _session_ptr->_eops.form_crf2trf( scan_end );
                        ivg::Matrix azel_sta1_end = bl.first->calc_az_el(scan_end, k, crf2trf_end ); // [rad]
                        ivg::Matrix azel_sta2_end = bl.second->calc_az_el(scan_end, k, crf2trf_end ); // [rad]

                        //  check if source is still visible at the end of the scheduled observation duration
                        if(azel_sta1_end(1) > min_ele && azel_sta2_end(1) > min_ele &&
                           bl.first->check_visibility(azel_sta1_end(0),azel_sta1_end(1)) && bl.second->check_visibility(azel_sta2_end(0),azel_sta2_end(1)))
                        {
                            
                            double snr_x, snr_s;
                            double sigma_final = compute_sigma_snr(duration, *bl.first, *bl.second, azel_sta1(1), azel_sta2(1), epoch_new_scan,
                                                                  xband, sband, snr_x, snr_s);
                            
                           
                            // set calculated SNR for this observation
                            obs_new.set_snr( snr_x, snr_s );
                            
                            //std::cout  << sta1 << "-" << sta2 << " " << snr_x << " " << snr_s << std::endl;

                            double new_wrap1,d_azi1,new_wrap2,d_azi2;
                            // check wrap-tolerance as well (e.g. sometimes too close to wrap edge)
                            bool wrap1_OK = bl.first->calc_wrap_and_rotangle( last_scan.get_cable_wrap(bl.first), azel_sta1_end(0), new_wrap1, d_azi1);
                            bool wrap2_OK = bl.second->calc_wrap_and_rotangle( last_scan.get_cable_wrap(bl.second), azel_sta2_end(0), new_wrap2, d_azi2);
                            
                            // check max and min degree of slewing (e.g. 160° and 15°)
                            if(d_azi1 <= (double)(*setup)["SKED"]["max_slew"]*(M_PI/180.0) && d_azi2 <= (double)(*setup)["SKED"]["max_slew"]*(M_PI/180.0)
                               && d_azi1 >= (double)(*setup)["SKED"]["min_slew"]*(M_PI/180.0) && d_azi2 >= (double)(*setup)["SKED"]["min_slew"]*(M_PI/180.0)     )
                            {
                                if(wrap1_OK && wrap2_OK)
                                {
                                    // set new wraps for this source in this scan for these two stations
                                    scan.set_cable_wrap(bl.first, new_wrap1);
                                    scan.set_cable_wrap(bl.second, new_wrap2);

                                    // save aiming position of telescopes at the end of the scan
                                    // -> later on we need this for skyplots
                                    scan.set_aiming(bl.first, azel_sta1_end);
                                    scan.set_aiming(bl.second, azel_sta2_end);

                                    // set fake delay and final sigma
                                    obs_new.set_delay( 0.0, sigma_final, 0.0 );
                                    // further features have to be set for calc_delay in case of impact factors
                                    obs_new.set_cable_cal(0.0,0.0);
                                    // adding meteorological data for simulation
                                    scan.add_scan_meteorology(sta1_idx, -999.000, -999.000, -99900.000, 0, ext_data_type, ext_met_data, gpt2filename );
                                    scan.add_scan_meteorology(sta2_idx, -999.000, -999.000, -99900.000, 0, ext_data_type, ext_met_data, gpt2filename );
                                    obs_new.set_ion_corr( 0.0 , 0.0, 0.0 , 0.0 );

                                    // command window output
    //                                cerr << src.get_name(ivg::srcname::ivs) << " " << setprecision(3) << sun_dist << " ";
    //                                cerr << " wrap1: " << new_wrap1*(180/M_PI) << " daz1: " << d_azi1*(180/M_PI) << " ";
    //                                cerr << " wrap2: " << new_wrap2*(180/M_PI) << " daz2: " << d_azi2*(180/M_PI) << " ";
    //                                cerr << " az: " << azel_sta1(0)*180.0/M_PI << "|" << azel_sta1_end(0)*180.0/M_PI <<  " el: " << azel_sta1(1)*180.0/M_PI << "|" << azel_sta1_end(1)*180.0/M_PI<< " ";
    //                                cerr << " az: " << azel_sta2(0)*180.0/M_PI << "|" << azel_sta2_end(0)*180.0/M_PI << " el: " << azel_sta2(1)*180.0/M_PI << "|" << azel_sta2_end(1)*180.0/M_PI<< " ";
    //                                cerr << " x_fl: " << x_flux << " s_fl: " << s_flux << " dur_x: " << dur_x << " dur_s: " << dur_s;
    //                                cerr << " snr_x: " << snr_x << " snr_s: " << snr_s << " duration: " << duration << " slew: " << slew_time << " sigma: " << setprecision(10) << sigma;
    //                                cerr << " f_sigma: " << sigma_final << endl;
//                                    if(bl_cnt == 0)
//                                        cerr << infostr;
                                    
                                    //store observation as candidate
                                    obs_new.set_use_flag( true );
                                    scan.add_obs( obs_new );
                                }
                                else
                                reason_cnt[" not suitable due too wrap-tolerance ("+bl_name+")"]++;
                            }
                            else
                            reason_cnt[" not suitable because max slew degree reached ("+bl_name+")"]++;
                        }
                        else
                            reason_cnt[" not visible AFTER scheduled observation duration ("+bl_name+")"]++;
                    }
                    else
                        reason_cnt[" longer scheduled than max_duration ("+bl_name+")"]++;
                }
                else
                    reason_cnt[" not visible ("+bl_name+")"]++;
                
            bl_cnt++;
            }
            
            // BETA STATUS! THIS NEED TO BE OPTIMIZED! -> e.g. SUBNET
            // scan is ONLY a candidate if there is atleast one observation
            // and ONLY if all telescopes can participate
            if(scan.get_nobs() > 0 &&  scan.get_nobs() == baselines.size())
                potential_scans.push_back(scan);
        }
        else
            reason_cnt[" too close to the sun"]++;
        
        reason_cnt[" all potential sources"]++;
    }
    
    // information output
    for(auto &rea: reason_cnt)
        cerr << "*** " << rea.second << rea.first << endl;
    
    cerr << "*** " << potential_scans.size() << " potential sources left" << endl;;
    
    return potential_scans;
    
#ifdef DEBUG_VLBI
   cerr << "--- vector<ivg::Scan> Schedule::get_next_potential_scans(ivg::Scan last_scan)" << " : " << tim.toc() << " s " << endl;
#endif
}
// ...........................................................................
ivg::Scan Schedule::determine_optimal_scan(ivg::skedtype type, vector<ivg::Scan> potential_scans)
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ ivg::Scan Schedule::determine_optimal_scan(ivg::skedtype, vector<ivg::Scan> )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
    // different approaches will determine different optimal scans based on a vector of potential scans
    ivg::Scan optimal_scan;
    if(type == ivg::skedtype::minsigma)
    {
        double min_scan_sigma = 10e10;
        for(auto &scan: potential_scans)
        {
            // only choose scan if scan mean group delay sigma is smaller than from all other scans
            double current_mean_sigma = scan.mean_group_delay_sigma();
            if( current_mean_sigma < min_scan_sigma)
            {
                min_scan_sigma = current_mean_sigma;
                optimal_scan = scan;
            }
        }
    }
    else if(type == ivg::skedtype::random)
    {
        int X = rand()%(potential_scans.size());
        return potential_scans.at(X);
    }
    else if(type == ivg::skedtype::coverage)
    {
        double max_scan_sky_coverage = 0.0;
        double dur_time = 1.0e10;
        double min_align_time = 1.0e10;
        
        std::vector< std::vector< lps::Node<ivg::Obs*, lps::Wedge>  > >  roots;
        double objective = 0.0;
        _session_ptr->compute_sky_coverage(_ea_grid_setup, _tg, roots, objective);
        
        ivg::Date last_epoch = _scheduled_scans.back().get_epoch();
        
        for(auto &scan: potential_scans)
        {
            vector<ivg::Scan> next_scheduled_scans = _scheduled_scans;
            next_scheduled_scans.push_back(scan);
            
            double align_time = scan.get_epoch().get_double_mjd() - last_epoch.get_double_mjd(); 
            
            // TODO everything is computed from scratch. better implementation with hemisphere data (still not the slowest part)
            std::map<int, std::vector<double> > coverage = _session_ptr->compute_sky_coverage(_ea_grid_setup, _tg, next_scheduled_scans, roots, objective);
            bool isBetter = false;
            
            if( objective > max_scan_sky_coverage)
            {
                isBetter = true;
            } else if (objective == max_scan_sky_coverage  && scan.get_scheduled_duration() < dur_time) {
                isBetter = true;
            } else if (objective == max_scan_sky_coverage  && scan.get_scheduled_duration() == dur_time && align_time < min_align_time) {
                isBetter = true;
            }
            
            if( isBetter ){
               max_scan_sky_coverage = objective;
               dur_time = scan.get_scheduled_duration();
               min_align_time = align_time;
               optimal_scan = scan; 
            }
        }
        
        
    }
    else if(type == ivg::skedtype::impact)
    {
        //suppress log output
        loglevel tmp_verbose = g_verbose;
        g_verbose = loglevel::NOTHING;
        
        // BETA STATUS  ||  BETA STATUS  || BETA STATUS
        // TODO too slow! unecessary work. with each scan a new line has to be added to jacobian matrix. 
        // But in each iteration entire matrix is computed!
        ivg::Date start = _scheduled_scans.at(0).get_epoch();
        ivg::Date end = _scheduled_scans.back().get_epoch();
        ivg::Date mean = ivg::Date(0.5*(start.get_double_mjd() + end.get_double_mjd()));
        
        double min_scan_impact_std = 10e10;
        for(auto &scan: potential_scans)
        {
            vector<ivg::Scan> next_scheduled_scans = _scheduled_scans;
            next_scheduled_scans.push_back(scan);
            _session_ptr->_scans = next_scheduled_scans; // update _scans within session
            _session_ptr->_nobs += scan.get_nobs(); // increase number of observations
            _session_ptr->_nobs_orig += scan.get_nobs(); // increase number of observations
            
            _session_ptr->_param_list = Param_list(_session_ptr->_trf, _session_ptr->_crf, _session_ptr->_eops, mean);
            _session_ptr->_param_list.set_start_end_epoch(start,end);
             
            _session_ptr->init_vgosdb_ngs_solution();
            _session_ptr->modify_parameterization();
            
            // _session_ptr->reduce_and_constrain(); // Do we need this?!
            ivg::Matrix impacts = _session_ptr->_lsa_solution.get_design0_ptr()->calc_impact_factors();
            
            // only choose scan if impactfactor of last insert observation is the greatest one
            double current_scan_impact = impacts(impacts.rows()-1);
            if( current_scan_impact < min_scan_impact_std)
            {
                min_scan_impact_std = current_scan_impact;
                optimal_scan = scan;
            }
            
            // we have to re-adjust the previously changed member variables
            _session_ptr->_scans = _scheduled_scans;
            _session_ptr->_nobs -= scan.get_nobs();
            _session_ptr->_nobs_orig -= scan.get_nobs();
        }
        // restore loglevel
        g_verbose = tmp_verbose;
    }
    
#ifdef DEBUG_VLBI
   cerr << "--- ivg::Scan Schedule::determine_optimal_scan(ivg::skedtype , vector<ivg::Scan> )" << " : " << tim.toc() << " s " << endl;
#endif
    return optimal_scan;
}
// ...........................................................................
void Schedule::remove_min_time_src(vector<ivg::Scan> &potential_scans, double min_time_src)
// ...........................................................................
{
    // compare potential_scans and all already _schedulded_scans
    // is there a source observed within the last min_time_src? if yes -> remove from potential_scans
    for(auto &scan_done: _scheduled_scans)
    {
        for (vector<ivg::Scan>::iterator scan_pot = potential_scans.begin(); scan_pot != potential_scans.end(); ++scan_pot)
        {
            // in case of identical source
            if(scan_pot->get_source()->get_name(ivs) == scan_done.get_source()->get_name(ivs) )
            {
                // check if source has been observed only a short time before (min_time_src)
                double timespan = abs(scan_pot->get_epoch().get_double_mjd() - scan_done.get_epoch().get_double_mjd());
                timespan *= 86400.0; // day to seconds
                if(timespan < min_time_src)
                {
                    //remove scan from potential list
                    potential_scans.erase(scan_pot);
                    scan_pot--;
                }
            }
        }
    }
}
// ...........................................................................
ivg::Scan Schedule::create_neutralpoint_scan(ivg::Date epoch)
// ...........................................................................
{
    // using default constructor for neutralpoint-scan
    // -> _source is NULL
    ivg::Scan neutralpoint;  
    neutralpoint.set_epoch(epoch);
    neutralpoint.set_source(NULL);
    
    // calculate and set cable wrap of neutralpoint
    for(int i=0; i< _session_ptr->_trf.get_number_stations(); i++)
    {
        ivg::Analysis_station *sta_i = _session_ptr->_trf.get_station(i);
        neutralpoint.set_cable_wrap(sta_i, (sta_i->set_antenna_info().azi_max + sta_i->set_antenna_info().azi_min) * 0.5);
    }
    
    return neutralpoint;
}
// ...........................................................................
void Schedule::get_skyplot_data(map<ivg::Analysis_station*, ivg::Matrix> &data, map<ivg::Analysis_station*, vector<string> > &tooltips)
// ...........................................................................
// CALLED WITHIN INDEP_MAIN.CPP FOR DEVELOPING USING:
// 
//    map<ivg::Analysis_station*, ivg::Matrix> azel_stations;
//    map<ivg::Analysis_station*, vector<string> > tooltips;
//    schedulator.get_skyplot_data(azel_stations,tooltips);
//
//    int cnt = 0;
//    Plot skyplot;
//    for(auto &station: azel_stations)
//    {
//        string name = station.first->get_name(ivg::staname::ivs_name);
//        ivg::Matrix azel = station.second;
//
//        QCPProjection projection2 = {protype::polarstereo, 30.0, 30.0, true, false, 1.0, Qt::red, 1.0, Qt::black, 0.0, QCPLineEnding::EndingStyle::esSpikeArrow, 8.0, 8.0, 20.0};
//        skyplot.projection(projection2, azel.get_col(0), azel.get_col(1), {QColor(ivg::color_values.at(cnt).c_str()), 2.5, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, QString::fromStdString(name), 4.0 , 2.0 }, "Skyplot", tooltips[station.first]);
//
//        cnt++;
//    }
//    a.exec();
{
    for(auto &scan: _scheduled_scans)
    { 
        map<ivg::Analysis_station*, ivg::Matrix> azel_stations = scan.get_aiming();
        for(auto &sta: azel_stations)
        {
            data[sta.first].append_rows(sta.second.transpose());
            tooltips[sta.first].push_back(scan.get_source()->get_name(ivg::srcname::ivs));
        }
    }
}
// ...........................................................................
void Schedule::show_scans(vector<ivg::Scan> &scans, string info)
// ...........................................................................
{
    ostringstream ss;
    ss << "+++ SCANS [" << scans.size() << "]: " << info << endl;
    ss << "    Source            Startepoch  Dur." << endl; 
    
    double mjd_pre = scans.at(0).get_epoch().get_double_mjd();
    for(auto &scan: scans)
    {
        if( _session_ptr->_trf.get_number_stations() > 3 & abs(scan.get_epoch().get_double_mjd() - mjd_pre) > 10/86400 ){
            ss << "......................................." << std::string(3*_session_ptr->_trf.get_number_stations(),'.') << std::endl;
        } 
        ss << setfill(' ') << setw(10) << scan.get_source()->get_name(ivs) << " ";        
        ss << setfill(' ') << setw(21) << scan.get_epoch().get_date_time("YYYY-MON-DD HH:MI:SS") << " ";
        ss << setfill(' ') << setw(5) << right << setprecision(1) << fixed << scan.get_scheduled_duration() << " ";
        for(int i=0; i< _session_ptr->_trf.get_number_stations(); i++)
        {
            if( scan.includes_station( _session_ptr->_trf.get_station(i) ) ){
                ss << setw(3) << _session_ptr->_trf.get_station(i)->get_name(ivg::staname::lettercode); 
            } else {
                ss << "   ";
            }
        }
        ss << " | ";
        for(int i = 0; i < scan.get_nobs(); ++i){
            ivg::Obs* obs = scan.get_obs_ptr(i);
            std::string sta1(""), sta2("");
            obs->get_station_names(sta1, sta2, ivg::staname::lettercode);
            ss << sta1 << sta2 << " "; 
        }
        ss << endl;
        
        mjd_pre = scan.get_epoch().get_double_mjd();
    }
    ss << "--- SCANS [" << scans.size() << "]: " << info << endl;
    cerr << ss.str();
}

// ...........................................................................
void Schedule::compute_band_info (ivg::Source& source, ivg::Analysis_station& sta1, Analysis_station& sta2,
                                  ivg::Date& epoch, ivg::BandInfo& xband, ivg::BandInfo& sband ) const{
// ...........................................................................
    ivg::Simulation simulator; 
    simulator.compute_band_info( source, sta1, sta2, epoch, _bit_sampling, xband, sband );
    
}

// ...........................................................................
double Schedule::compute_obs_duration(ivg::Source& source, ivg::Analysis_station& sta1, Analysis_station& sta2, ivg::Date& epoch, ivg::BandInfo& xband, ivg::BandInfo& sband, bool applyMinDur) const{
// ...........................................................................
    
    ivg::Simulation simulator; // we need this for some privat simulation functions
    Setting *setup = _session_ptr->_setup;
    
    bool broadband = (*setup)["SKED"].exists("broadband") && (*setup)["SKED"]["broadband"]["apply"];
    
    double duration = 0;
    
    if( broadband ){
        duration = simulator._calc_obs_duration(xband.sefd_sta2, xband.sefd_sta1, xband.flux, _bandwidth, (*setup)["SKED"]["broadband"]["channels"], _bit_sampling, (*setup)["SKED"]["snr_min_x"], (*setup)["SKED"]["const_sync"]);
    } else {
        double dur_x = simulator._calc_obs_duration(xband.sefd_sta2, xband.sefd_sta1, xband.flux, _bandwidth, xband.num_chan, _bit_sampling, (*setup)["SKED"]["snr_min_x"], (*setup)["SKED"]["const_sync"]);
        double dur_s = simulator._calc_obs_duration(sband.sefd_sta2, sband.sefd_sta1, sband.flux, _bandwidth, sband.num_chan, _bit_sampling, (*setup)["SKED"]["snr_min_s"], (*setup)["SKED"]["const_sync"]);

        duration = max( dur_x, dur_s );
        duration = ceil(duration);
    }
    
    if(applyMinDur){
        double min_scan = (*setup)["SKED"]["min_scan"];

        // get final observation-duration for current source
        duration = max( duration, min_scan ); 
    }
    
    return duration;   
    
}


// ...........................................................................
void Schedule::compute_snr(double duration, const ivg::BandInfo& xband, const ivg::BandInfo& sband, double& snr_x, double& snr_s) const{
// ...........................................................................
    Setting *setup = _session_ptr->_setup;
    ivg::Simulation simulator; // we need this for some privat simulation functions
    
    bool broadband = (*setup)["SKED"].exists("broadband") && (*setup)["SKED"]["broadband"]["apply"];
    if( broadband ){
        snr_x = simulator._snr_achieved(xband.sefd_sta1, xband.sefd_sta2, xband.flux, _bandwidth, (*setup)["SKED"]["broadband"]["channels"], _bit_sampling, duration, (*setup)["SKED"]["const_sync"]);
        snr_s = snr_x;
    } else {
        snr_x = simulator._snr_achieved(xband.sefd_sta1, xband.sefd_sta2, xband.flux, _bandwidth, xband.num_chan, _bit_sampling, duration, (*setup)["SKED"]["const_sync"]);
        snr_s = simulator._snr_achieved(sband.sefd_sta1, sband.sefd_sta2, sband.flux, _bandwidth, sband.num_chan, _bit_sampling, duration, (*setup)["SKED"]["const_sync"]);
    }
   
}


// ...........................................................................
double Schedule::compute_sigma_snr(double duration, ivg::Analysis_station& sta1, Analysis_station& sta2, double el_sta1,
        double el_sta2,  ivg::Date& epoch, ivg::BandInfo& xband, ivg::BandInfo& sband, double& snr_x, double& snr_s){
// ...........................................................................
    
    ivg::Simulation simulator; // we need this for some privat simulation functions
    Setting *setup = _session_ptr->_setup;
    
    
    bool broadband = (*setup)["SKED"].exists("broadband") && (*setup)["SKED"]["broadband"]["apply"];
    if( broadband ){
        snr_x = simulator._snr_achieved(xband.sefd_sta1, xband.sefd_sta2, xband.flux, _bandwidth, (*setup)["SKED"]["broadband"]["channels"], _bit_sampling, duration, (*setup)["SKED"]["const_sync"]);
        snr_s = snr_x;
    } else {
        snr_x = simulator._snr_achieved(xband.sefd_sta1, xband.sefd_sta2, xband.flux, _bandwidth, xband.num_chan, _bit_sampling, duration, (*setup)["SKED"]["const_sync"]);
        snr_s = simulator._snr_achieved(sband.sefd_sta1, sband.sefd_sta2, sband.flux, _bandwidth, sband.num_chan, _bit_sampling, duration, (*setup)["SKED"]["const_sync"]);
    }

    // additional noise (very dominant for sigma calculation)
    double noise = (double)(*setup)["SKED"]["add_ps"] * 1e-12; // [ps]
    double sigma = simulator._std_sked( noise, xband.mean_freq, sband.mean_freq, sta1.calc_bandwidth_rms(ivg::band::X), sta1.calc_bandwidth_rms(ivg::band::S), snr_x, snr_s );

    double mfh_sta1, mfw_sta1, mfh_sta2, mfw_sta2;
    ivg::Troposphere( &sta1, epoch ).calc_gmf(el_sta1, mfh_sta1, mfw_sta1);
    ivg::Troposphere( &sta2, epoch ).calc_gmf(el_sta2, mfh_sta2, mfw_sta2);

    double sigma_final = sqrt( pow(sigma,2.0) + pow(6.0e-12*mfw_sta1,2.0) + pow(6.0e-12*mfw_sta2,2.0) );
     
    return sigma_final;
    
}

 std::vector<unsigned> Schedule::number_of_stations_source_is_visible(const lps::Transits& transits) const{
    _session_ptr->get_trf_ptr()->create_station_indices();
    std::vector<unsigned> visble_from_stas(_session_ptr->get_crf_ptr()->get_number_sources_inc_unused(), 0);
    for(size_t sou_idx =0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx){

        // check if all observatories can observe the same source
        for( std::vector<ivg::Analysis_station*>& stas : _session_ptr->get_trf_ptr()->get_station_twin_list() ){
            for(ivg::Analysis_station* sta: stas){
                unsigned sta_idx = sta->get_idx();
                if( transits.hasTransits(sou_idx, sta_idx) ){
                    visble_from_stas[sou_idx]++;
                    break;
                }
            }
        }
    }
    return visble_from_stas;
 }
 
std::vector<bool> Schedule::quasar_has_at_least_one_bl_with_good_SNR(const lps::Transits& transits) const{
    std::vector<bool> has_good_snr(_session_ptr->get_crf_ptr()->get_number_sources_inc_unused(), false);

    double maxScanLength = (*_session_ptr->get_setup())["SKED"]["max_scan"];
    double x_thres = (double)(*_session_ptr->get_setup())["SKED"]["snr_min_x"]-(double)(*_session_ptr->get_setup())["SKED"]["snr_margin_x"];
    double s_thres = (double)(*_session_ptr->get_setup())["SKED"]["snr_min_s"]-(double)(*_session_ptr->get_setup())["SKED"]["snr_margin_s"];
        
    lps::Seconds duration = (_session_ptr->getEnd().get_double_mjd() - _session_ptr->getStart().get_double_mjd())*3600*24;
    for(size_t sou_idx =0; sou_idx < _session_ptr->get_crf_ptr()->get_number_sources_inc_unused(); ++sou_idx){
        for(size_t sta_idx_1=0; sta_idx_1 < _session_ptr->get_trf_ptr()->get_number_stations() ;++sta_idx_1){
            for(size_t sta_idx_2=sta_idx_1+1; sta_idx_2 < _session_ptr->get_trf_ptr()->get_number_stations() ;++sta_idx_2){
                for( double t = 0.0; t < duration; t+=30.0 ){
                    if( transits.sourceIsVisible(sou_idx, sta_idx_1, t) && transits.sourceIsVisible(sou_idx, sta_idx_2, t) ){
                        ivg::Date epoch = _session_ptr->getStart();
                        epoch.add_secs(t);
                        ivg::BandInfo xband, sband;
                        compute_band_info (*_session_ptr->get_crf_ptr()->get_source( sou_idx ), *_session_ptr->get_trf_ptr()->get_station(sta_idx_1),
                                              *_session_ptr->get_trf_ptr()->get_station(sta_idx_2), epoch, xband, sband );
                        double snr_x(0.0), snr_s(0.0);
                        compute_snr(maxScanLength, xband, sband, snr_x, snr_s);
                        
                        if( snr_x > x_thres && snr_s > s_thres){
                            has_good_snr[sou_idx] = true;
                            break;
                        }
                    }
                }
                if(has_good_snr[sou_idx])
                    break;
            }
            if(has_good_snr[sou_idx])
                break;
        }
    }
    return has_good_snr;
}



std::map<std::string, ivg::Matrix> Schedule::minElevation(ivg::Date epoch) const{
     
    size_t n_sta = _session_ptr->_trf.get_number_stations();
    _session_ptr->_trf.create_station_indices();
    
    std::vector<std::pair<std::string, std::string>> exclude_bl (0);
    if(  (*_session_ptr->_setup)["SKED"].exists("exclude_bl") ){
        exclude_bl = get_baselines( (*_session_ptr->_setup)["SKED"]["exclude_bl"]);
    }
    
    std::vector< std::vector<std::string> > networks;
    for(size_t sta_idx_1 = 0; sta_idx_1 < n_sta; ++sta_idx_1) {
        ivg::Analysis_station* sta_1 = _session_ptr->get_trf_ptr()->get_station(sta_idx_1);
        std::string sta_name_1 = sta_1->get_name( ivg::staname::ivs_name );
        for(size_t sta_idx_2 = sta_idx_1+1; sta_idx_2 < n_sta; ++sta_idx_2) {
            ivg::Analysis_station* sta_2 = _session_ptr->get_trf_ptr()->get_station(sta_idx_2);
            std::string sta_name_2 = sta_2->get_name( ivg::staname::ivs_name );
            bool exclude = includes_baseline(exclude_bl, sta_name_1, sta_name_2 );
            if(exclude == false && _session_ptr->get_trf_ptr()->areTwins(sta_idx_1, sta_idx_2) == false){
                std::vector<string> bl = { sta_name_1 , sta_name_2 };
                networks.push_back( bl );
            }
        }
    }
      
    return minElevation( epoch, networks);
}

std::map<std::string, ivg::Matrix> Schedule::minElevation(ivg::Date epoch, std::vector< std::vector<std::string> > networks) const{
    size_t n_sta = _session_ptr->_trf.get_number_stations();
        
    std::map<std::string, ivg::Matrix > minE;
    for(size_t sta_idx_1 = 0; sta_idx_1 < n_sta; ++sta_idx_1) {
        std::string sta_name_1 = _session_ptr->get_trf_ptr()->get_station( sta_idx_1 )->get_name(ivg::staname::ivs_name);
        minE[sta_name_1] = ivg::Matrix(360, 1, M_PI/2.0);
    }
        
    for( std::vector<std::string>& station_network :  networks ) {
        for( std::string& sta_name_1 : station_network){
            
            ivg::Analysis_station* sta_1;
            _session_ptr->get_trf_ptr()->get_station(&sta_1, sta_name_1);
            
            for(int az = 0; az < 360; az++ ){
                for(double el = 0.0; el < M_PI/2.0; el+=(.1*ivg::d2rad) ){
                    if( sta_1->check_visibility( az*ivg::d2rad, el ) ){
                        ivg::Matrix k_trf = sta_1->azel2k(az*ivg::d2rad, el, epoch);

                        bool seenFromOthers = true;
                        for( std::string& sta_name_2 : station_network){
                                                        
                            ivg::Analysis_station* sta_2;
                            _session_ptr->get_trf_ptr()->get_station(&sta_2, sta_name_2);
                            
                            if( sta_name_1 == sta_name_2)
                                continue;
                            
                            ivg::Matrix azel_2 = sta_2->calc_az_el( k_trf );
                            if( sta_2->check_visibility( azel_2(0), azel_2(1) ) == false ){
                                seenFromOthers = false;
                                break;
                            }
                        }
                        
                        if(seenFromOthers){
                            minE[sta_name_1](az) = min( minE[sta_name_1](az), el );
                            break;
                        }
                    }
                }                
            }
        }
    }
        
    return minE;
    
}

} //namespace