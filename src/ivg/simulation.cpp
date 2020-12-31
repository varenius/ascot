#include "simulation.h"
#include "matrix.h"

namespace ivg
{

// ===========================================================================
// 		constructors
// ===========================================================================

// ...........................................................................
Simulation::Simulation()
// ...........................................................................
{
}
// ...........................................................................
Simulation::Simulation( ivg::Session *session_ptr )
// ...........................................................................
{
   _session_ptr = session_ptr;
}

// ===========================================================================
// 		public methods
// ===========================================================================


void Simulation::compute_reference_tropo_triangular_mat(  ivg::Trf& trf, const Setting& setup ) const{
    
     //{1, 6, 12, 18, 24, 29, 35, 40, 44, 48, 51, 54, 56, 57, 58}

    std::string path = (const char*)(setup)["SIM"]["Reference"]["troposphere"]["VCMpath"];

    
    ivg::Matrix azel; //in rad
    std::vector<ivg::Date> epochs;
    
    ivg::HemisphereData hemisphere_data = init_hemisphere_data(setup["SIM"]["Reference"]["troposphere"]);
    hemisphere_data.set_name("grid");
    
    hemisphere_data.getGrid(azel,  epochs);

    hemisphere_data.save_grid(path);
    
    ivg::Turbulence turb;
    turb.set_param_form_cfg( setup, trf  );
    
    for( int sta_idx = 0; sta_idx < trf.get_number_stations(); ++sta_idx ) {

        std::string sta_name = trf.get_station(sta_idx)->get_name(ivg::staname::ivs_name);

        ivg::Matrix VCM = turb.calc_vcm_vienna_model( trf.get_station(sta_idx), azel, epochs, 200.0 );

          // scale: mm^2 to s^2
        VCM = VCM * ( 1e-6 / pow(ivg::c,2.0) );
        
        ivg::Matrix R = VCM;
        R.chol( R, 'L' );

        R.save_bin( path + "/" + sta_name + "_R.dat");

    }
}




// ...........................................................................
void Simulation::simulate_reference_troposhere_and_clock( ivg::Trf& trf, double zwd0, const std::vector<std::string>& databases, const Setting& setup)
// ...........................................................................
{
    std::string path = (const char*)(setup)["SIM"]["Reference"]["path"];
 
    std::string VCMpath = (const char*)(setup)["SIM"]["Reference"]["troposphere"]["VCMpath"];
    
    // Troposphere -------------------------------------------------------------
       
    ivg::HemisphereData hemisphere_data = init_hemisphere_data(setup["SIM"]["Reference"]["troposphere"]);
    unsigned npos(0), nepochs(0);
    hemisphere_data.getDimension(npos, nepochs );
    
    for(int sta_idx = 0; sta_idx < trf.get_number_stations(); ++sta_idx){
        
        ivg::Matrix R;
        std::string sta_name = trf.get_station(sta_idx)->get_name(ivg::staname::ivs_name);
        R.load_bin(VCMpath + sta_name + "_R.dat");

        for(const std::string& db : databases){

            std::cerr << std::string(10,'\b') << setw(10) << db;

            std::string dbpath = path + "/" + db;
            if (!directory_exists(dbpath)){
                mkdir(dbpath.c_str() , 0700);
            }

            ivg::Matrix w( R.size(1), 1, 0.0 );
            w.rand_norm( 0.0, 1.0 ) ;   

            // get initial equivalent zenith wet delays 
            // (3ps ~ 1mm)
            ivg::Matrix iEZWD( R.size(1), 1, zwd0 );    

            // simulate equivalent zenith wet delays
            ivg::Matrix EZWD = iEZWD + R * w ;

            EZWD.vec2Mat( nepochs );

            EZWD.save_bin(dbpath + "/" + sta_name + "_EZWD.dat");
        }
    }
    std::cout << std::endl;
    
    // clocks ------------------------------------------------------------------
    
    ivg::Matrix mjd = hemisphere_data.getIntervalBeginningsMJD();
     
    // list of all stations in solution
    vector<std::string> sta_names = trf.get_station_names( ivg::staname::ivs_name );
    std::sort( sta_names.begin(),sta_names.end() );  
    std::vector< std::string > grp_names;
     

    // get clock setting from cfg
    Setting &params = (setup)[ "SIM" ];
    Setting &groups = (setup)[ "groups" ];

    std::map<string, std::set<std::string> > twinmap = trf.get_twins_map_including_all();
       
    for( int j=0;j<params[ "clocks" ].getLength();++j ){
        Setting &setup = params[ "clocks" ][ j ];
        std::vector< std::string > grp_names;
        group2names( setup, groups, "stations",sta_names,grp_names );

        for( int i=0;i<grp_names.size();++i )
        { 
            double astd = (double)setup[ "allan_std_dev" ];  // Allan standard deviation 
            double tau = (double)setup[ "interval" ];        // time span for which ASD is valid
            
            
            // create clocks for all co-located station too
            // TODO config file entry 
            for(std::string sta: twinmap[ grp_names[ i ] ]){
//                 std::cerr << sta << std::endl;
                for(const std::string& db : databases){
                    ivg::Matrix clock = _sim_station_clock( astd, tau, mjd );

                    clock.save_bin(path + "/" + db + "/" + sta + "_clock.dat");

                }
            }
        }
    }
        
}

ivg::HemisphereData Simulation::init_hemisphere_data(const Setting& setup, ivg::Date start) const{
    Setting &numCelsCfg = (setup)["numCells"];
    std::vector<unsigned> numCells( numCelsCfg.getLength(), 0 );

    for(int i=0; i<numCelsCfg.getLength(); i++){
            numCells[i] = (unsigned)numCelsCfg[i];
    }

    double dt = (double)(setup)["dt"]; // in seconds
    double duration = (double)(setup)["duration"]; // in seconds

    return ivg::HemisphereData( start, duration, dt, numCells );
}

// ...........................................................................
void Simulation::simulate()
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ void Simulation::simulate()" << endl; 
   tictoc tim;
   tim.tic();
#endif
//   
//    ivg::HemisphereData test( ivg::Date(0.0), 60, 30, {1,4,8} );
//       
//    std::vector<double> d({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25});
//    ivg::Matrix D (d);
//    D.vec2Mat(2);
//    D.show();
//    test.set_data( D  );
//    
//    ivg::Matrix azel = test.get_azel();
//    
//    for(int r = 0; 0 < azel.rows(); ++r){
//        std::cerr <<azel(r,0)*ivg::rad2d << " " << azel(r,1)*ivg::rad2d << " t1 " <<   test.get_data(azel(r,0), azel(r,1), 0.0 ) << std::endl;
//        std::cerr <<azel(r,0)*ivg::rad2d << " " << azel(r,1)*ivg::rad2d << " t2 " <<   test.get_data(azel(r,0), azel(r,1), 30.0 ) << std::endl;
//    }
//    
//    exit(0);
//   
   Setting &params = (*_session_ptr->_setup)[ "SIM" ];
   Setting &groups = (*_session_ptr->_setup)[ "groups" ];
   
    bool referenceMode = (bool)(*_session_ptr->_setup)["SIM"].exists( "Reference" ) &&
                          (int)(*_session_ptr->_setup)["SIM"]["Reference"]["step"] == 2 &&
                         (bool)(*_session_ptr->_setup)["SIM"]["Reference"]["apply"];
    std::string ReferenceSimpath = ""; 
    ivg::HemisphereData hemisphere_data;
    ivg::Matrix refSimEpochs;
    if ( referenceMode ){
        
        log<INFO>("*** using pre computed reference values for clocks and zwd" ); 
        
        ReferenceSimpath = (const char*)(*_session_ptr->_setup)["SIM"]["Reference"]["path"];
        
        hemisphere_data = init_hemisphere_data((*_session_ptr->_setup)["SIM"]["Reference"]["troposphere"], _session_ptr->_start);
        
        refSimEpochs = hemisphere_data.getIntervalBeginningsMJD();
//        refSimEpochs.show(12);
    }
   
    int obs_counter = 0;
    ivg::Obs* obs_ptr = nullptr;
    
    ivg::Matrix whiteNoise( _session_ptr->_nobs, 1, 0.0 ), sigmaTau( _session_ptr->_nobs, 1, 0.0),
                tau( _session_ptr->_nobs, 1, 0.0 ), clockObs(  _session_ptr->_nobs, 1, 0.0), tropo( _session_ptr->_nobs, 1, 0.0 );
    
    std::map< std::string, std::vector<double> > ezwd_per_sta;
    std::map< std::string, std::vector<double> > mjd_per_sta;
   
    if ( !params.exists("load_saved_observations") || ! params["load_saved_observations"] ){

        log<INFO>("*** Simulate observations #") % _session_ptr->_nobs;

        // loop over scans and observations therein and create matrix of epochs
        // the last clock value is wrong => add another one 30s after last observation (Corinna will check this)
        ivg::Matrix mjd( _session_ptr->_nobs+1,1,0.0 );
        int num_bl = ( _session_ptr->_trf.get_number_stations()* (_session_ptr->_trf.get_number_stations() - 1 ) ) / 2;
        int num_of_params = _session_ptr->_trf.get_number_stations()*7  +_session_ptr->_crf.get_number_sources()*2 + 5 + num_bl;

        // transpose of the jacobian matrix and an iterator to first element of a column of A^T
        ivg::Matrix AT( num_of_params,1,0.0 );
        vector<double>::iterator at_iter = AT.begin();

        ivg::Matrix apr0( num_of_params, 1,0.0 );
        vector<double>::iterator apr0_iter = apr0.begin();       

        for( vector<ivg::Scan>::iterator scan_it = _session_ptr->_scans.begin(); scan_it != _session_ptr->_scans.end(); ++scan_it )
        { 
           for(int obs_i=0; obs_i < scan_it->get_nobs(); ++obs_i)
           { 
              mjd( obs_counter ) = scan_it->get_epoch().get_double_mjd();
              obs_ptr = scan_it->get_obs_ptr(obs_i);
              
              // group delay (observed) is set to zero --> o-c calculated in obs_ptr->calc_delay is -c 'minus computed'
              obs_ptr->set_delay( 0.0, 1.0e-11, 0.0 );
              obs_ptr->calc_delay( at_iter, apr0_iter );
              at_iter = AT.begin();
              obs_counter++;
           }
        }

        // the last clock value is wrong => add another one 30s after last observation (Corinna will check this)
        mjd( mjd.rows()-1 ) = mjd( mjd.rows()-2 ) + 30.0 / 86400.0 ;

        // list of all stations in solution
        vector<std::string> sta_names = _session_ptr->_trf.get_station_names( ivg::staname::ivs_name );
        std::sort( sta_names.begin(),sta_names.end() );  
        std::vector< std::string > grp_names;

        // create list of all possible baselines within solution
        vector<std::string> bl_names;  
        for( int i=0;i<sta_names.size();++i )
           for( int j=i+1;j<sta_names.size();++j )
               bl_names.push_back(sta_names.at( i )+"-"+sta_names.at( j ));

        // -------------------------------------------------
        // simulate CLOCKS for each station for _all_ epochs 
        std::map< std::string,ivg::Matrix > clocks;  

        for( int j=0;j<params[ "clocks" ].getLength();++j )
        {
           Setting &setup = params[ "clocks" ][ j ];
           std::vector< std::string > grp_names;
           group2names( setup, groups, "stations",sta_names,grp_names );

           for( int i=0;i<grp_names.size();++i )
           { 
                if( !referenceMode){
               
                    double astd = (double)setup[ "allan_std_dev" ];  // Allan standard deviation 
                    double tau = (double)setup[ "interval" ];        // time span for which ASD is valid
          //            cerr << "CLOCK " <<  grp_names.at( i ) << " " << astd << " " << tau << endl;
                    clocks[ grp_names.at( i ) ] = _sim_station_clock( astd, tau, mjd );
                } else {
                    ivg::Matrix simulatedClockValues;
                    
                    ivg::Matrix clock(mjd.rows(),1,0.0);
                   // cerr << "CLOCK " <<  grp_names.at( i ) << endl;
                    
                    simulatedClockValues.load_bin( ReferenceSimpath +"/" + _session_ptr->_name_in_masterfile + "/" + grp_names.at( i ) +  "_clock.dat" );
                    unsigned cidx = 0;
                    for(double m : mjd){
//                        std::cerr << std::setprecision(16) << m << std::endl;
                        vector<int> leq = refSimEpochs.find_idx(le, m+1e-6);
                        //show_vector(leq);
                        if(leq.size() < 1){
                            throw runtime_error( "void Simulation::simulate() no simulated clock values found" );
                        } 
                        
                        unsigned idx = leq.back();
                        clock(cidx) = simulatedClockValues(idx);
                        
                        cidx++;
                    }
                    clocks[ grp_names.at( i ) ] = clock;

                }
           }
        }

        // -------------------------------------------------
        // simulate WHITE NOISE for all baselines for _all_ epochs
        bool wnFromSNR = params["white_noise_from_SNR"];
        std::map< std::string, ivg::Matrix > bl_wn;   // maybe these should be initialized first!!!
        if(!wnFromSNR){
            for( int j=0;j<params[ "white_noise" ].getLength();++j )
            {
               Setting &setup = params[ "white_noise" ][ j ];
               group2names( setup, groups,"baselines",bl_names,grp_names );

               for( int i=0;i<grp_names.size();++i )
               { 
                   double sigma = (double)setup[ "std_dev" ];  // white noise std from config file
        //             cerr << "WN " <<  grp_names.at( i ) << " " << sigma << endl;
                   bl_wn[ grp_names.at( i ) ] = _sim_baseline_white_noise( sigma );
               }
            }
        }

        // -------------------------------------------------
        // simulate TROPO
        
        if( !referenceMode){
            std::map< std::string, turbulence_data > turb_sta = _session_ptr->_turbulence.set_param_form_cfg( *_session_ptr->get_setup(), *_session_ptr->get_trf_ptr()  );
            std::string model = params["troposphere"][ "model" ];    // turbulence model 
            tropo = _sim_tropo( model, turb_sta ); // ZWD from sessions covariance matrix is not ZWD! (Sebastians issue)
        } else {
            
            std::map<string, std::set<std::string> > twinmap = _session_ptr->_trf.get_twins_map_including_all();

            std::map<std::string, ivg::HemisphereData> referenceTroposhere;
            for(unsigned sta_idx = 0; sta_idx < _session_ptr->_trf.get_number_stations(); ++sta_idx){
                std::string sta_name = _session_ptr->_trf.get_station(sta_idx)->get_name(ivg::staname::ivs_name);
                referenceTroposhere[sta_name] = hemisphere_data;
                referenceTroposhere[sta_name].set_name(sta_name);
                
                std::string path = ReferenceSimpath +"/" + _session_ptr->_name_in_masterfile + "/" + sta_name +  "_EZWD.dat";
                
                bool referenceFileExists = true;
                if( !file_exists(path) ){
                    referenceFileExists = false;                    

                    for(std::string sta: twinmap[sta_name]){
                        path = ReferenceSimpath +"/" + _session_ptr->_name_in_masterfile + "/" + sta +  "_EZWD.dat";
                        if ( file_exists(path) ){
                            referenceFileExists = true;
                            break;
                        }
                    }
                }
                
                if(referenceFileExists){
                    referenceTroposhere[sta_name].load(path);
                } else {
                    throw runtime_error( " no reference troposphere found for station " + sta_name );
                }
                
            }
                
            
            obs_counter = 0;
            for( vector<ivg::Scan>::iterator scan_it = _session_ptr->_scans.begin(); scan_it != _session_ptr->_scans.end(); ++scan_it )
            { 
                for(int obs_i=0; obs_i < scan_it->get_nobs(); ++obs_i)
                { 
                    ivg::Date epoch = scan_it->get_epoch();
                    epoch.add_secs( scan_it->get_scheduled_duration()/2.0 );

                    obs_ptr = scan_it->get_obs_ptr(obs_i);
                    std::string sta_name1, sta_name2;
                    obs_ptr->get_station_names( sta_name1, sta_name2 );

                    ivg::Matrix azel1 = obs_ptr->get_az_el(1);
                    ivg::Matrix azel2 = obs_ptr->get_az_el(2);
                    
                    double ezwd1 = referenceTroposhere[sta_name1].get_data( azel1(0), azel1(1), epoch );
                    double ezwd2 = referenceTroposhere[sta_name2].get_data( azel2(0), azel2(1), epoch );
                    
                    ezwd_per_sta[sta_name1].push_back(ezwd1);
                    ezwd_per_sta[sta_name2].push_back(ezwd2);
                    
                    mjd_per_sta[sta_name1].push_back(epoch.get_double_mjd());
                    mjd_per_sta[sta_name2].push_back(epoch.get_double_mjd());
                                        
//                    double ezwd1 = referenceTroposhere[sta_name1].get_data( 0.0, M_PI/2.0, epoch );
//                    double ezwd2 = referenceTroposhere[sta_name2].get_data( 0.0, M_PI/2.0, epoch );
                    
                    double mfh, mfw1, mfw2;
                    scan_it->get_data(obs_ptr->get_scan_idx(1)).tropo.get_mapping_function(mfh,mfw1);
                    scan_it->get_data(obs_ptr->get_scan_idx(2)).tropo.get_mapping_function(mfh,mfw2);
                                        
                    tropo(obs_counter) = ezwd2*mfw2  - ezwd1*mfw1 ;

                    obs_counter++;
                }
            }              
           
        }
        // -------------------------------------------------
        // calc std and create wn vector

        log<INFO>("*** Calculating standard deviation and storing oc + std");
        
        if( wnFromSNR ){
            whiteNoise.rand_norm(0.0, 1.0);
        }
        
        // loop over scans and observations therein and set o-c values
        obs_counter = 0;
        std::string sta1, sta2;
        double sigma_ion;
        //double mfh, mfw1, mfw2;
        ivg::Analysis_station* sta1_ptr = nullptr;
        ivg::Analysis_station* sta2_ptr = nullptr;
        for( vector<ivg::Scan>::iterator scan_it = _session_ptr->_scans.begin(); scan_it != _session_ptr->_scans.end(); ++scan_it )
        {
           for(int obs_i=0; obs_i < scan_it->get_nobs(); ++obs_i)
           { 
              ivg::Date epo = scan_it->get_epoch();
              obs_ptr = scan_it->get_obs_ptr( obs_i );
              obs_ptr->get_station_names( sta1, sta2 );

              // get station pointer
              _session_ptr->_trf.get_station(&sta1_ptr,sta1);
              _session_ptr->_trf.get_station(&sta2_ptr,sta2);

              //scan_it->get_data(obs_ptr->get_scan_idx(1)).tropo.get_mapping_function(mfh,mfw1);
              //scan_it->get_data(obs_ptr->get_scan_idx(2)).tropo.get_mapping_function(mfh,mfw2);

              clockObs(obs_counter) = clocks[ sta2 ]( obs_counter ) - clocks[ sta1 ]( obs_counter );

               // calculate standard devation of simulated group delay
              double duration = scan_it->get_scheduled_duration();
              if(duration < 5.0)
                  throw runtime_error("!!! Observation-Duration of "+std::to_string(duration)+"secs cannot be correct.");

              _calc_standard_deviation(scan_it->get_source(),sta1_ptr, sta2_ptr, &epo, duration, sigmaTau(obs_counter), sigma_ion);

              if( wnFromSNR ){
                  whiteNoise(obs_counter) = whiteNoise(obs_counter)*sigmaTau(obs_counter) ;
              } else {
                //determine correct baseline-name
                string bl_name = sta1+"-"+sta2;
                if(bl_wn[ bl_name ].rows() == 0)
                    bl_name = sta2+"-"+sta1;
                whiteNoise(obs_counter) = bl_wn[ bl_name ]( obs_counter );
              }

              // negative computed delay
              tau(obs_counter) = obs_ptr->get_observed_minus_computed();

              obs_counter++;
           }
        }
        
        // save simulated observations
        std::string outdir = (const char*)(*_session_ptr->get_setup())[ "outdir" ];
        
        log<INFO>("*** save simulated observations to ") % outdir;
        
        whiteNoise.save_bin(outdir + "/simWhiteNoise.bin");
        sigmaTau.save_bin(outdir + "/simSigmaTau.bin");
        tau.save_bin(outdir + "/simTau.bin");
        clockObs.save_bin(outdir + "/simClocks.bin");
        tropo.save_bin(outdir + "/simTropo.bin");
//        
//        for( std::string& sta : _session_ptr->_trf.get_station_names(ivg::staname::ivs_name) ){
//            ivg::Matrix tmp ( mjd_per_sta[sta] );
//            ivg::Matrix tmp2 ( ezwd_per_sta[sta] );
//            tmp.append_cols(tmp2);
//            tmp.save_bin(outdir + "/sim" + sta + "ezwd.bin");
//        }
        
    } else {
        std::string outdir = (const char*)(*_session_ptr->get_setup())[ "outdir" ];
        
        log<INFO>("*** Load simulated observations from ") % outdir;
        
        whiteNoise.load_bin(outdir + "/simWhiteNoise.bin");
        sigmaTau.load_bin(outdir + "/simSigmaTau.bin");
        tau.load_bin(outdir + "/simTau.bin");
        clockObs.load_bin(outdir + "/simClocks.bin");
        tropo.load_bin(outdir + "/simTropo.bin");
    }
   
    // -------------------------------------------------
    // Add simulated values to session
    obs_counter = 0;
    double oc = 0.0;
    for( vector<ivg::Scan>::iterator scan_it = _session_ptr->_scans.begin(); scan_it != _session_ptr->_scans.end(); ++scan_it )
    {
       for(int obs_i=0; obs_i < scan_it->get_nobs(); ++obs_i)
       { 
          obs_ptr = scan_it->get_obs_ptr( obs_i );
                 
           //sum up stochastic components
          oc = clockObs( obs_counter ) + whiteNoise(obs_counter) + tropo(obs_counter);

          // only delay needed oc is recomputed anyway
          //obs_ptr->set_observed_minus_computed( oc );
           
//          oc = 0.0;
//          if( (bool)(*_session_ptr->get_setup())[ "SIM" ]["comp"]["zwd"] )
//              oc += tropo(obs_counter);
//          if( (bool)(*_session_ptr->get_setup())[ "SIM" ]["comp"]["clo"] )
//              oc +=  clockObs( obs_counter );
//          if( (bool)(*_session_ptr->get_setup())[ "SIM" ]["comp"]["wn"] )
//              oc += whiteNoise(obs_counter);
                      
          obs_ptr->set_delay( -tau(obs_counter) + oc, sigmaTau(obs_counter), 0.0 );
          obs_counter++;
       }
    }

    log<INFO>("*** Simulation finished");
    
#ifdef DEBUG_VLBI
   cerr << "--- void Simulation::simulate()" << " : " << tim.toc() << " s " << endl;
#endif
}

// ===========================================================================
// 		private methods
// ===========================================================================

// ...........................................................................
ivg::Matrix Simulation::_sim_station_clock( double allan_std, double tau, ivg::Matrix epochs )
/*  This method simulates station clock variations by a sum of power law
 *  time series:
 *     - white phase modulation
 *     - flicker phase modulation
 *     - white frequency modulation
 *     - flicker frequency modulation
 *     - random walk frequency modulation
 *
 *  Sigmas used for the simulation werde adopted to represent a 1e-14@50 min
 *  clock, however, the scaling is also valid for other noise levels.
 *  The process realizes clocks which are dominated by white and flicker PM
 *  within the first 30 sec. Until approximately 45 min. white FM dominates,
 *  finally flicker and random walk FM lead to a smoth trend change of the
 *  Allan STD plot.
 *
 *  References:
 *   Jeremy Kasdin,
 *   Discrete Simulation of Colored Noise and Stochastic Processes
 *   and 1/f^a Power Law Noise Generation,
 *   Proceedings of the IEEE, Volume 83, Number 5, 1995, pages 802-827.
 *
 *   Jacques Rutman
 *   Characterization of Phase and Frequency Instabilities in Precision
 *   Frequency Sources: Fifteen Years of Progress, In: PROCEEDINGS OF THE
 *   IEEE, VOL. 66, NO. 9, SEPTEMBER 1978
 * .........................................................................*/
{
#ifdef DEBUG_VLBI
   cerr << "+++ ivg::Matrix Simulation::_sim_station_clock( double , double , ivg::Matrix )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   ivg::Matrix var = _calc_noise_variances( pow( allan_std,2.0 ), tau*60.0, 1.0 );

   // SIMULATION of ...
   // ... white phase modulation
   ivg::Matrix x_0 = _f_alpha( epochs,var( 0 ),0.0 );

   // ... flicker phase modulation
   ivg::Matrix x_1 = _f_alpha( epochs,var( 1 ),1.0 );

   // ... white frequency modulation
   ivg::Matrix x_2 = _f_alpha( epochs,var( 2 ),2.0 );

   // ... flicker frequency modulation
   ivg::Matrix x_3 = _f_alpha( epochs,var( 3 ),3.0 );

   // ... random walk frequency modulation
   ivg::Matrix x_4 = _f_alpha( epochs,var( 4 ),4.0 );

#ifdef DEBUG_VLBI
   cerr << "--- ivg::Matrix Simulation::_sim_station_clock( double , double , ivg::Matrix )" << " : " << tim.toc() << " s " << endl;
#endif
   return( x_0 + x_1 + x_2 + x_3 + x_4 );
}

// ...........................................................................
ivg::Matrix Simulation::_sim_baseline_white_noise( double sigma )
// ...........................................................................
{
   ivg::Matrix noise( _session_ptr->_nobs, 1, 0.0 );
   noise.rand_norm( 0.0, sigma );

   return noise;
}
// ...........................................................................
ivg::Matrix Simulation::_sim_tropo(std::string model, std::map< std::string, turbulence_data > turb_sta)
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ ivg::Matrix Simulation::_sim_tropo(std::string , std::map< std::string, turbulence_data > )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
    // calculate variance covariance matrix
    _session_ptr->_turbulence.set_turb_params(turb_sta);
    //_session_ptr->_turbulence.show();
    ivg::Matrix VCM = _session_ptr->_model_turbulence(model,turb_sta);
    // VCM.save_bin( "/home/corbin/VCM.dat");
    // cholesky decomposition of the variance-covariance matrix and 	
    ivg::Matrix D = VCM;
    D.chol(D);

    // create vector of zero mean Gaussian random numbers with variance one	
    ivg::Matrix w(D.size(2),1,0.0);
    w.rand_norm(0.0,1.0);

    // get initial equivalent zenith wet delays   
    ivg::Matrix iEZWD(D.size(1),1,0.0);

    // simulate equivalent zenith wet delays
    ivg::Matrix EZWD = iEZWD+D* w;
    
#ifdef DEBUG_VLBI
   cerr << "--- ivg::Matrix Simulation::_sim_tropo(std::string , std::map< std::string, turbulence_data > )" << " : " << tim.toc() << " s " << endl;
#endif
   return EZWD;
}
// ...........................................................................
ivg::Matrix Simulation::_calc_noise_variances( double allan_var, double tau, double f_h )
/* determine variances of the power law processes (Rutman, 1978) and scale 
 * them by empirical factors
 *  
 *   Jacques Rutman
 *   Characterization of Phase and Frequency Instabilities in Precision
 *   Frequency Sources: Fifteen Years of Progress, In: PROCEEDINGS OF THE
 *   IEEE, VOL. 66, NO. 9, SEPTEMBER 1978
// .........................................................................*/
{
   ivg::Matrix var( 5,1,0.0 );

   // contribution of white PM to entire noise floor
   var( 0 ) = 4*pow( M_PI,2.0 )*pow( tau,2.0 )/3/f_h;

   // contribution of flicker PM to entire noise floor
   var( 1 ) = 4*pow( M_PI,2.0 )*pow( tau,2.0 )/3/log(8.88*f_h*tau);

   // contribution of white PM to entire noise floor
   var( 2 ) = 2*tau;

   // contribution of white PM to entire noise floor
   var( 3 ) = 1/log(4);

   // contribution of white PM to entire noise floor
   var( 4 ) = 3/(2*pow( M_PI,2.0 )*tau);

   var = var*allan_var;

   // empirical adjustment
   var( 0 ) = var( 0 )*1e-6;
   var( 1 ) = var( 1 )*1e-5;
   var( 2 ) = var( 2 )*1e-2;
   var( 3 ) = var( 3 )*1e-1;
   var( 4 ) = var( 4 )*5e-1;

   return var;
}
// ...........................................................................
ivg::Matrix Simulation::_f_alpha( ivg::Matrix epochs, double Q_d, double alpha )
// Generate power law noise f^alpha. Simulation is done on a 1 sec basis and
// then interpolated to given epochs.
//
//  Reference:
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes and 1/f^a
//    Power Law Noise Generation, Proceedings of the IEEE, Volume 83,
//    Number 5, 1995, pages 802-827.
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ Matrix vlbiSim::_f_alpha( Matrix, double, double )" << endl;
   tictoc tim;
   tim.tic();
#endif
   // initialize random number generator with zero mean and standard
   // deviation = sqrt( Q_d )
   Q_d = sqrt( Q_d ); /* find the deviation of the noise */

   double ha = alpha/2.0;

   // determine number of points for a 1 sec spacing and convert epochs
   // vector to seconds
   epochs = epochs-epochs.min();
   epochs = epochs * 86400.0;
   int n_pts = (int)epochs.max()+1;

   // generate the coefficients hk and create the sequence wk with white
   // noise
   ivg::Matrix hfa( 2*n_pts,1,0.0 );
   ivg::Matrix wfa( n_pts, 1, 0.0 );
   wfa.rand_norm( 0.0, Q_d );
   wfa.resize( 2*n_pts, 1, 0.0 );

   hfa( 0 ) = 1.0;
   for( int i=1; i<n_pts; ++i )
      hfa( i ) = hfa( i-1 )*( ha+(double)( i-1 ) )/(double)( i );

   // perform the discrete Fourier transform
   complex<double>* hfa_freq = new std::complex<double> [2*n_pts];
   fftw_plan p1 =  fftw_plan_dft_r2c_1d( hfa.length(), hfa.data_ptr(),
                                         reinterpret_cast<fftw_complex*>(hfa_freq),
                                         FFTW_ESTIMATE );
   fftw_execute(p1);
   complex<double>* wfa_freq = new std::complex<double> [2*n_pts];
   fftw_plan p2 =  fftw_plan_dft_r2c_1d( wfa.length(), wfa.data_ptr(),
                                         reinterpret_cast<fftw_complex*>(wfa_freq),
                                         FFTW_ESTIMATE );
   fftw_execute(p2);

   // multiply the two complex vectors
   for( int i=0; i<=n_pts; ++i )
      wfa_freq[i] = complex<double>
      (
         wfa_freq[i].real()*hfa_freq[i].real()-wfa_freq[i].imag()*hfa_freq[i].imag(),
         wfa_freq[i].real()*hfa_freq[i].imag()+wfa_freq[i].imag()*hfa_freq[i].real()
      );

   // Scale to match the conventions of the Numerical Recipes FFT code
   wfa_freq[0] /= 2.0;
   wfa_freq[n_pts-1] /= 2.0;

   // inverse Fourier transform the result
   ivg::Matrix x( 2*n_pts,1,0.0 );
   fftw_plan p3 =  fftw_plan_dft_c2r_1d( 2*n_pts,
                                         reinterpret_cast<fftw_complex*>(wfa_freq),
                                         x.data_ptr(), FFTW_ESTIMATE );
   fftw_execute( p3 );

   // create the final result
   x = x.get_sub( 0,0,n_pts-1,0 ) * ( 1.0/(double)n_pts );

   // free memory
   fftw_destroy_plan(p1);
   fftw_destroy_plan(p2);
   fftw_destroy_plan(p3);
   delete [] hfa_freq;
   delete [] wfa_freq;

   // interpolate to given epochs
   ivg::Matrix cl( epochs.rows(),1 );
   ivg::Matrix t0( 0.0,1.0,x.length(),1 );
   for( int i=0; i<epochs.length(); ++i )
      cl(i) = x.interpolate( t0, epochs(i), "linear" )(0);

#ifdef DEBUG_VLBI
   cerr << "--- Matrix vlbiSim::_f_alpha( Matrix, double, double )"
        << ": " << tim.toc() << " s " << endl;
#endif
   return cl;
}


// ...........................................................................
void Simulation::compute_band_info (ivg::Source& source, ivg::Analysis_station& sta1, Analysis_station& sta2,
                                    ivg::Date& epoch, double bit_sampling, ivg::BandInfo& xband, ivg::BandInfo& sband ) const{
// ...........................................................................
    
    // mean frequencies (based on freq.cat -> ivg::parser::freq_cat)
    // x-band: in general 8 channels used for mean-calculation
    // s-band: in general 5 channels used for mean-calculation
    xband.mean_freq  = sta1.calc_mean_frequency(ivg::band::X); //8595.49;
    sband.mean_freq  = sta1.calc_mean_frequency(ivg::band::S); //2290.99;
    double xwave = 1.0e6 * xband.mean_freq / ivg::c;
    double swave = 1.0e6 * sband.mean_freq / ivg::c;

    // calculate flux density in X and S-band (based on flux.cat -> ivg::parser::flux_cat)
    ivg::Matrix bl_vector = sta1.calc_xyz(epoch) - sta2.calc_xyz(epoch);
    xband.flux = source.calc_flux(ivg::band::X,xwave,bl_vector,&epoch);
    sband.flux = source.calc_flux(ivg::band::S,swave,bl_vector,&epoch);

    // Overall  antenna  performance  is  measured  by  System Equivalent Flux Density = SEFD 
    // map zenith SEFD to source direction (based on equip.cat -> ivg::parser::equip_cat)
    xband.sefd_sta1 = sta1.calc_slant_sefd(ivg::band::X, &source, &epoch);
    sband.sefd_sta1 = sta1.calc_slant_sefd(ivg::band::S, &source, &epoch);
    xband.sefd_sta2 = sta2.calc_slant_sefd(ivg::band::X, &source, &epoch);
    sband.sefd_sta2 = sta2.calc_slant_sefd(ivg::band::S, &source, &epoch); 
    
    xband.num_chan = bit_sampling * (sta1.get_nfreq(ivg::band::X) + 2); // additional 2 for sideband
    sband.num_chan = bit_sampling * (sta1.get_nfreq(ivg::band::S));
}


// ...........................................................................
void Simulation::_calc_standard_deviation(ivg::Source *src_ptr, ivg::Analysis_station *&sta1_ptr, 
                                          ivg::Analysis_station *&sta2_ptr, ivg::Date *epoch_ptr,
                                          double duration, double &sigma_snr_x, double &sigma_ion)
// ...........................................................................
{
    
    Setting& setup = *_session_ptr->_setup;
    // recorder name syntax: tracks-channels-fanout-bits
    string rec_name = setup["SKED"]["rec_name"]; // e.g. "00-16-0-1"
    int bit_sampling = std::stoi(rec_name.substr(rec_name.size()-1,1)); // last number of rec_name
    
    ivg::BandInfo xband, sband;
    compute_band_info (*src_ptr, *sta1_ptr, *sta2_ptr, *epoch_ptr, bit_sampling, xband, sband );
    
    double snrx(0.0), snrs(0.0);
    
    bool broadband = setup["SKED"].exists("broadband") && setup["SKED"]["broadband"]["apply"];
    if( broadband ){
        
        double bandwidth =  ((double) setup["SKED"]["broadband"]["bandwidth"] )*1.0e6 ; // in Hz
        snrx = _snr_achieved(xband.sefd_sta1, xband.sefd_sta2, xband.flux, bandwidth, setup["SKED"]["broadband"]["channels"], bit_sampling, duration, setup["SKED"]["const_sync"]);
        
        // TODO remove hard coded freq setup
        // vgos setp in MHz
        std::vector<double> tmp = {3448.4, 3416.4, 3352.4, 3288.4, 3192.4, 3064.4, 3032.4, 3000.4,
                                   5688.4, 5656.4, 5592.4, 5528.4, 5432.4, 5304.4, 5272.4, 5240.4,
                                   6808.4, 6776.4, 6712.4, 6648.4, 6552.4, 6424.4, 6392.4, 6360.4,
                                  10648.4, 10616.4, 10552.4, 10488.4, 10392.4, 10264.4, 10232.4, 10200.4};
        ivg::Matrix freq(tmp);
        double bbrms = sta1_ptr->calc_bandwidth_rms(freq);
        sigma_snr_x =  1.0e6/(2.0*M_PI*bbrms*snrx)*1e-12; // in ps
        sigma_ion = 0.0;
    } else {
        
        double bandwidth = ((double)setup["SKED"]["bandwidth"])*1.0e6; // in Hz
        
        snrx = _snr_achieved(xband.sefd_sta1, xband.sefd_sta2, xband.flux, bandwidth, xband.num_chan, bit_sampling, duration, setup["SKED"]["const_sync"]);
        snrs = _snr_achieved(sband.sefd_sta1, sband.sefd_sta2, sband.flux, bandwidth, sband.num_chan, bit_sampling, duration, setup["SKED"]["const_sync"]);
        // calculate standard deviation based on SNR
        double bwrmsx = sta1_ptr->calc_bandwidth_rms(ivg::band::X);
        sigma_snr_x = 1.0e6/(2.0*M_PI*bwrmsx*snrx)*1e-12;
        double bwrmss = sta1_ptr->calc_bandwidth_rms(ivg::band::S);
        double sigma_snr_s = 1.0e6/(2.0*M_PI*bwrmss*snrs)*1e-12;

        // calculate standard deviation of ionosperic correction (see, e.g.,  phd 
        // thesis Tesmer Chapter 2)
        double ffact = pow(sband.mean_freq,2.0)/(pow(xband.mean_freq,2.0)-pow(sband.mean_freq,2.0));
        sigma_ion = ffact*sqrt(pow(sigma_snr_x,2.0)+pow(sigma_snr_s,2.0));

    }

    
}
// ...........................................................................
double Simulation::_calc_obs_duration(const double & sefd1,const double & sefd2, 
                                      const double & flux,double bw, int num_channels, int bit_samp,
                                      const double & snr_min,const double & sync_time)
// ...........................................................................
{
    // empirical values based on a combination of bit- and correlator-efficiency
    // see snrac.f or SkedManual_v2016Dec09
    // using these values lead to no SNR-warning during sked check and seems to be correct!
    double nu;
    if(bit_samp == 1)
        nu = 0.62;
    else if(bit_samp == 2)
        nu = 0.58;
    
    double duration = pow((snr_min/(nu * flux)),2.0)*((sefd1*sefd2)/(2.0 *bw*(double)num_channels))+sync_time;
    
    // old version based on jSked
    // double duration = (sefd1*sefd2) /(2.0*bw*(double)num_channels)*pow(((1.75*snr_min)/flux),2.0)+0.5+sync_time;
    
    return duration;
}
// ...........................................................................
double Simulation::_snr_achieved(const double & sefd1, const double & sefd2,
                                 const double & flux,double bw, int num_channels, int bit_samp,
                                 const double & duration,const double & sync_time )
// calculate  actual SNRs achieved in the schedule
// !!! has to be compared with SKEDs snrac !!!
// ...........................................................................
{
    // see _calc_obs_duration
    double nu;
    if(bit_samp == 1)
        nu = 0.62;
    else if(bit_samp == 2)
        nu = 0.58;
    
   double snr = (flux*nu) * sqrt((2.0*bw*(double)num_channels*duration-sync_time)/(sefd1*sefd2));
   return snr;
}

// ...........................................................................
// sigma like sked
double Simulation::_std_sked( double noise, double mean_freq_x, double mean_freq_s, double bw_rms_x, double bw_rms_s, double snr_x, double snr_s  )
// ...........................................................................
{
   double ffact = pow(mean_freq_s,2.0) / ( pow(mean_freq_x,2.0) - pow(mean_freq_s,2.0) );
   
   double sigma_snr_x = 1.0e6 / ( 2.0*M_PI * bw_rms_x * snr_x )*1e-12;
   double sigma_snr_s = 1.0e6 / ( 2.0*M_PI * bw_rms_s * snr_s )*1e-12;
   
   double sigma_ion = ffact * sqrt( pow(sigma_snr_x,2.0) + pow(sigma_snr_s,2.0) );
   
   double sigma_snr = sqrt( pow(sigma_ion,2.0) + pow(sigma_snr_x,2.0) );
   
   double sigma = sqrt( pow(sigma_snr,2.0) + pow(noise,2.0) );
   
   return sigma;
}

} // namespace


