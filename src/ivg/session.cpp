#include "session.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Session::Session()
// ...........................................................................
{
}
// ...........................................................................
Session::Session( const ivg::Session &other )
// ...........................................................................
{
    (*this) = other;
}
// ...........................................................................
Session & Session::operator=( const Session &other )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Session & Session::operator=( const Session &other )" << endl;
   tictoc tim;
   tim.tic();
#endif

        _setup = other._setup;  // Pointer!!!
        _handling = other._handling; // Pointer!!!
        _name = other._name;      
        _name_in_masterfile = other._name_in_masterfile;
        _file_comment << other._file_comment.str();
        _type = other._type;

        _scans.resize( other._scans.size() );
        copy( other._scans.begin(), other._scans.end(), _scans.begin() );

        _trf = other._trf;
        _crf = other._crf;

        _eops = other._eops;
        _ephem = other._ephem; // Pointer!!!
        
        _sta_geophys_effects.resize( other._sta_geophys_effects.size() );
        copy( other._sta_geophys_effects.begin(), other._sta_geophys_effects.end(), _sta_geophys_effects.begin() );

        _start = other._start;
        _end = other._end;

        _nobs  = other._nobs;
        
        _nobs_orig = other._nobs_orig;
        
        _origin_obs_idxs.resize( other._origin_obs_idxs.size() );
        copy( other._origin_obs_idxs.begin(), other._origin_obs_idxs.end(), _origin_obs_idxs.begin() );

        // Least squares method object
        _lsa_solution = other._lsa_solution;
        _neq_solution = other._neq_solution;
        _icls_solution = other._icls_solution;
        _lsc_solution = other._lsc_solution;
        
        
        if(_type.find("snx") != string::npos)
            _solution = &_neq_solution; // Pointer!!!
        else if(_type == "ngs")
            _solution = &_lsa_solution; // Pointer!!!

        // parameter list
        _param_list = other._param_list;

        // turbulence 
        _turbulence = other._turbulence;;
        
        _disconts = other._disconts;
        
        _residuals = other._residuals;
        
        _obsstats = other._obsstats;
        
        _session_path = other._session_path;
        
        _nnr_nnt_set = other._nnr_nnt_set;

        _band_type = other._band_type;
        _ambigRes = other._ambigRes;
	_phaseDelay = other._phaseDelay;
	_SB_solution = other._SB_solution;
        _ambiguity_spacing = other._ambiguity_spacing;
        _num_ambig = other._num_ambig;
        _eff_freq = other._eff_freq;
        

#if DEBUG_VLBI >=2
   cerr << "--- Session & Session::operator=( const Session &other )" << " : " << tim.toc() << " s " << endl;
#endif 
   
    return *this;
}
// ...........................................................................
Session & Session::operator+=( Session &other )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Session & Session::operator+=( Session & )" << endl;
   tictoc tim;
   tim.tic();
#endif
        
    /*
     This is one of the most important function for stacking and combination!
     It's called after the session was initialized using the constructor Session(...) and init_snx_solution(...) was called (see global_main.cpp)
    */
   
    // operator only implemented in case of snx-sessions (not vgosDB or ngs)
    if( _type.find("snx") != string::npos )
    {        
            log<INFO>("*** Stacking SNX systems");

            vector<int> positions_in_this;

            int new_params_cnt=0;
            for(std::vector<ivg::Param>::iterator param_other = other._param_list.begin(); param_other != other._param_list.end(); ++param_other)
            {
                // get positions of parameters in reference NEQ (this-neq) and new system to be added (other-neq)
                int pos_other = param_other- other._param_list.begin();
                int pos_this = _param_list.get_index((*param_other));

                // if parameter is a xpo, xpor, ypo, ypor, ut1, lod, then we need to add these parameter from the other-neq to the this-neq
                // this is only the case, if it is a new session with a different date!
                // this is not the case, if it is a session with the same date but a different analysiscenter! (will be  considered later on)
                // bool new_erp = param_other->is_type({ivg::paramtype::xpo, ivg::paramtype::ypo, ivg::paramtype::ut1}, {0,1});
                bool new_erp = param_other->is_type({ivg::paramtype::xpo, ivg::paramtype::ypo, ivg::paramtype::ut1, ivg::paramtype::nutx, ivg::paramtype::nuty}, {0,1});
                
                // check if paramater is a source is a special handling source -> will be set up as local parameter (arc-parameter)
                bool sh_source = (find(special_handlings.begin(), special_handlings.end(), remove_spaces_end(param_other->get_name())) != special_handlings.end());
                // in case of ivg::ASCOT snx files it might be a clockbreak-parameter
                bool new_clockbreak = param_other->is_type({ivg::paramtype::cbr},{0});
                sh_source=false;                
                // check if special handling source really need to be set up new or if the identical source is already set up within a adjustable time range
                // ATTENTION: time interval is set manually! default is 45+/- days
                if(sh_source == true)
                {
                   vector<int> indexes = _param_list.get_idx(param_other->get_type(), param_other->get_name());
                   double other_doy = param_other->get_epoch().get_double_mjd();
                   
                   if(indexes.size()>0)
                   {
                       double this_doy = _param_list.get_param(indexes.back())->get_epoch().get_double_mjd();
                       
                       if(other_doy <= this_doy + 45.0) // use 45.0 instead
                       {
                           // source is already set up, so sh_source needs to be set to false
                           sh_source = false;
                           pos_this = indexes.back();
                       }
                   }
                }
                
                // check if a new eop parameter really need to be set up in the big system or if the identical eop is already set up
                // this might be the case if a session of a different analysis center is already imported
                if(new_erp == true)
                {
                    for(std::vector<ivg::Param>::iterator param_this = _param_list.begin(); param_this != _param_list.end(); ++param_this)
                    {
                        bool equal_name,equal_type,equal_order,equal_epoch;
                        param_other->compare_to((*param_this), equal_name, equal_type, equal_order, equal_epoch);    
                        if(equal_name && equal_type && equal_order && equal_epoch)
                        {
                            new_erp = false;
                            pos_this = param_this - _param_list.begin();
                            break;
                        }   
//                        cerr << "Name: " << equal_name << " | Type: " << equal_type << " | Order: " << equal_order << " | Epoch: " << equal_epoch << endl;
                    } 
                }
                
                // always set up as new parameter
                // there have been some analysis center contributions which included mistakenly zwds and gradients
                // these parameters always need to be set up new
                bool new_zwd_nrg_egr = param_other->is_type({zwd, ngr, egr},{0,1});

                // checking if cut off date for sources is reached, only for NON-DEFININGS 
                // this was done to stack the non-defining sources of all sessions until the cut off date of the ICRF2 ( 03.03.2009 = 54893.0mjd )
                // in order to compare the results with the official ones
//                if(param_other->is_type({ivg::paramtype::ra,ivg::paramtype::dec}, {0}))
//                {
//                    bool special = (find(special_handlings.begin(), special_handlings.end(), remove_spaces_end(param_other->get_name())) != special_handlings.end());
//                    bool defining = (find(icrf2_definings.begin(), icrf2_definings.end(), remove_spaces_end(param_other->get_name())) != icrf2_definings.end());
//                    if( !special && !defining )
//                    {
//                        if(pos_this != -1)
//                        {
//                            vector<int> posis = _param_list.get_idx(param_other->get_type(), param_other->get_name());
//                            if(param_other->get_epoch().get_double_mjd() > 54893.0 &&  posis.size() == 1)
//                            {
//                                pos_this = -1;
//                            }
//                            else if(param_other->get_epoch().get_double_mjd() > 54893.0 &&  posis.size() == 2)
//                            {
//                                pos_this = posis.at(1);
//                            }
//                        }
//                    }
//                }
                
                // MANUAL SWITCH IF SH SOURCE SHOULD BE SET UP LOCAL OR NOT
                // sh_source = false, then all special handling sources will be always set up as global!
//                sh_source = false;
                // detect parameter which need to be set up new in the big system
                if(new_erp == true || sh_source == true || new_clockbreak == true || new_zwd_nrg_egr == true || pos_this == -1  )
                {
                    stringstream ss;
                    if(new_erp)
                        ss << "*** New EOP " << param_other->get_typename() << " " << param_other->get_order() << " added at " << param_other->get_epoch().get_double_mjd();
                    else if(sh_source)
                        ss << "*** New Source " << param_other->get_typename() << " added at " << param_other->get_epoch().get_double_mjd();
                    else
                        ss << "*** New Param " << param_other->get_typename() << " " << param_other->get_order() << " added at " << param_other->get_epoch().get_double_mjd();

                    // insert the new parameter into the reference this-parameterlist 
                    _param_list.insert_param(_param_list.end(), (*param_other));

                    // keep trf and crf up-to-date because the method modify_parameterization() works on _crf and _trf
                    // in case of a station
                    if(param_other->is_type({ivg::paramtype::stax}, {0}))
                    {
                        ivg::Analysis_station *station_other;
                        if(other._trf.get_station(&station_other,param_other->get_name()))
                        {
                            ivg::Analysis_station* tmp;
			    if(!(_trf.get_station(&tmp,param_other->get_name())))
			      {
				ivg::Analysis_station new_station_this = (*station_other);
				_trf.push_back(new_station_this);
			     

				ss << ". Pushed " << new_station_this.get_name(ivg::staname::ivs_name) << " to _trf ";
			      }
                            
                        }
                        else
                            throw runtime_error( "Session & Session::operator+=( Session &other ): SHOULD NOT BE POSSIBLE (TRF)");
//                        {
//                            ivg::Matrix pos0(3,1,0.0);
//                            pos0(0) = param_other->get_apriori();
//                            pos0(1) = (param_other+1)->get_apriori();
//                            pos0(2) = (param_other+2)->get_apriori();
//                            
//                            vector<ivg::Date> none;
//                            map<ivg::staname, std::string> names = { {ivg::staname::ivs_name, param_other->get_name()} };
//                            Analysis_station unknown_station( pos0, ivg::Matrix(3,1,0.0), param_other->get_epoch(), none, names);
//                            _trf.push_back(unknown_station);
//                            
//                           log<WARNING>("*** CRITICAL: Unkown new Analysis_station initialized and added to _trf: ") % unknown_station.get_name(ivg::staname::ivs_name) % " at " % param_other->get_epoch().get_double_mjd();
//                           log<WARNING>("*** ") % unknown_station.get_name(ivg::staname::ivs_name) % " only contains pos0 and ivs_name!";
//                        }

                    }
                    // in case of a source
                    else if(param_other->is_type({ivg::paramtype::ra}, {0}))
                    {
                        ivg::Source *source_other;
                        if(other._crf.get_source(&source_other,param_other->get_name()))
                        {
                            ivg::Source new_source_this = (*source_other);
                            _crf.push_back(new_source_this);

                            ss << ". Pushed " << new_source_this.get_name(ivg::srcname::iers) << " to _crf [#" << _crf.get_number_sources() << "]";
                        }
                        else
                            throw runtime_error( "Session & Session::operator+=( Session &other ): SHOULD NOT BE POSSIBLE (CRF)");

                    }
                    // generate output
                    log<DETAIL>(ss.str());
                    
                    // save the position where this new detected parameter has been inserted
                    pos_this = _param_list.size()-1;

                    // save how many new columns/rows are needed for the neq-system
                    new_params_cnt++;
                }
                // in case of already existing parameter, we stack and don't insert!
                else
                {
                    // update the session-counter to know how often a parameter has been stacked
                    if(param_other->is_type({ivg::paramtype::ra}, {0}))
                    {
                        ivg::Source *source;
                        if(_crf.get_source(&source,param_other->get_name()))
                        {
                            // (right now only implemented in case of sources within a CRF)
                            source->increase_n_sessions(1);
			    source->increase_n_obs_in_sess(1);
                            log<DETAIL>("*** Old Source increased: ") % source->get_name(ivg::srcname::iers) % " at " % param_other->get_epoch().get_double_mjd();
                        }
                        else
                            throw runtime_error( "Session & Session::operator+=( Session &other ): Should never occur!");
                    }
                    
                    // this counter will be displayed in the results table (get_resultline())
                    _param_list.get_param(pos_this)->increment_stacked();
                    
                    //check if aprioris are equal
                    double apriori_other = param_other->get_apriori();
                    double apriori_this = _param_list.get_param(pos_this)->get_apriori();
                    
                    // the apriori values transformation is done in advance, so these errors shouldn't occur if the
                    // apriori transformation is activated (TRF2APR , EOP2APR, CRF2APR) in configfile, e.g. snx.cfg
                    if(param_other->is_type({ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz}, {0}))
                    {
                        double apriori_diff =  abs(apriori_this - apriori_other);
			
                        if(apriori_diff > 5.0e-5)
                            log<WARNING>("!!! Session & Session::operator+=( Session &other ): Difference of apriori coordinates bigger than 0.05mm ") % to_string(apriori_diff);
                        
                    }
                    else if (apriori_this != apriori_other)
                    {
                        double apriori_diff =  abs(apriori_this - apriori_other);
                        if(apriori_diff > 1.0e-8)
                            log<WARNING>("!!! Session & Session::operator+=( Session &other ): Difference of apriori EOPs apriori values bigger than 1.0e-8") % apriori_diff;
                    }
                    
                }
                
                positions_in_this.push_back(pos_this);
            }
            // up to now the book keeping is done. Now the NEQ (containing the numbers) must be adjusted!
            
            log<INFO>("*** ") % new_params_cnt % " Parameter added to _param_list";
            
            // prepare neq by resizing it if new parameters are found in other neq
            _neq_solution.enlarge_neq(new_params_cnt);

            // stack each element by element from other to this
            _neq_solution.stack_neq(other._neq_solution, positions_in_this);
            
            // end time of stacked need to be adjusted
            _end = other._end;            

#if DEBUG_VLBI >=2
   cerr << "--- Session & Session::operator+=( Session & )" << " : " << tim.toc() << " s " << endl;
#endif 
   
        return *this;
    }
    else
        throw runtime_error( "Session & Session::operator+=( Session &other ): Stacking databases with _type == ngs not implemented yet");
}

// ...........................................................................
Session::Session( Setting *setup, string name, void **ephem, int session_cnt )
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ Session::Session( Setting *setup, string name, void **ephem, int session_cnt)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // definitions from configfile
    _setup = setup;
    _name = name;
    _name_in_masterfile = name;
    if (session_cnt!=-1)
      _handling = &((*_setup)[ "sessions" ][ session_cnt ]);
    else
      {
	_handling= &((*_setup)[ "sessions" ][0]);
	if((bool)(_handling)->exists("handling"))
	   _handling->remove("handling");
      }
    // write database-name to filecomment
    _file_comment << "IVS_DB_NAME " << name << endl;
    
    for(int i=0; i<(*_setup)["stadisp"].getLength(); i++)
    {
        if((*_setup)["stadisp"][i][1])
            _sta_geophys_effects.push_back((*_setup)["stadisp"][i][0]);
    }
    
    // always add the "ECCENTRICITY" to the station displacements!
    // they need to be considered in every case!
    _sta_geophys_effects.push_back("ECCENTRICITY");

    // only load ephemerides once if pointer is NULL (default)
    if(ephem == NULL)
    {
        Setting &definitions = (*_setup)["definitions"];

        string ephfile_name = definitions["ephemerides"][(const char*)(
                                  *_setup)["ephemerides"]];
        char nams[400][6];
        double vals[400];

        _ephem = jpl_init_ephemeris(ephfile_name.c_str(), nams, vals);
    }
    else
        _ephem = (*ephem);

    // create EOP series
    ivg::Date date;
    if (_name.size()>9)
      date = ivg::Date( _name.substr(0,8), "YYYYMMDD" );
    else
      date = ivg::Date( _name.substr(0,7), "YYMMMDD" );
    
    if ( (*_setup).exists( "SIM" ) && (bool)(*_setup)[ "SIM" ]["apply"] && (*_setup)[ "SIM" ].exists("eop") ){
        init_session_eop( (*_setup)["SIM"]["eop"], date );
        
        if((*_setup)[ "SIM" ].exists("notut1")){
            ivg::Eop_series notut = init_eop((*_setup)["SIM"]["notut1"], date, 17.0);
            _eops.replace(notut, true, false, true );
        }
        
    } else {
        init_session_eop( (*_setup)["eop"], date );
    }
     
    _niter = 0;

    _ambigRes = (bool)(*_setup)["solve_ambig"];
    if (_setup->exists("phase_solution"))
      _phaseDelay = (bool)(*_setup)["phase_solution"];
    else
      _phaseDelay= false;
     if (_setup->exists("SB_solution"))
      _SB_solution = (bool)(*_setup)["SB_solution"];
    else
      _SB_solution= false;
    _nnr_nnt_set = false;
    
    // read band type from control file if possible. else set it to X band
    if (_setup->exists("band"))
    {
        std::string band =  (const char *)(*_setup)["band"];
        if( band.compare("X") == 0)
            _band_type = ivg::band::X;
        else if( band.compare("S") == 0)
            _band_type = ivg::band::S;
        else
        {
            log<WARNING>("!!! Unepected band-type. Using X band");
            _band_type = ivg::band::X;
        }
    }
    else
    {
        _band_type = ivg::band::X;
        log<INFO>("*** Band-type not specified in control file. Initialising using X-band. This setting could be changed during runtime by some programs.");
    }
    
#ifdef DEBUG_VLBI
   cerr << "--- Session::Session( Setting *setup, string name, void **ephem, int session_cnt )" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ===========================================================================
// public methods
// ===========================================================================
// ...........................................................................
std::string Session::get_name()
// ...........................................................................
{
    return _name;
}
// ...........................................................................
void Session::init_vgosdb_ngs_solution()
// ...........................................................................
{
#if DEBUG_VLBI >=1
   cerr << "+++ void Session::init_vgosdb_ngs_solution( )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
   log<INFO>("*** START-Session::init_vgosdb_ngs_solution with _name: ") % _name; 

    // number of initial parameters, i.e, 
    // * for each station: 
    //      - 3 station positions
    //      - 1 ZWD 
    //      - 1 clock offset 
    //      - 2 gradients 
    // * for each source
    //      - right ascension
    //      - declination
    // * EOPs
    //      - 2 polar motion values
    //      - UT1-TAI
    //      - 2 nutation components
    // * for each baseline:
    //      - 1 clock offset   
    int num_bl = ( _trf.get_number_stations()* (_trf.get_number_stations() - 1 ) ) / 2;
    int num_of_params = _trf.get_number_stations()*7 + _crf.get_number_sources()*2 + 5 + num_bl;
    
    ivg::Matrix* A_ptr   = nullptr;
    ivg::Matrix* oc_ptr  = nullptr;
    ivg::Matrix* epo_ptr = nullptr;
    ivg::Matrix* wgt_ptr = nullptr;            
    ivg::Matrix* aplo_ptr = nullptr; 
      
    // in case of ambiguity resolution group and single band delays are used
    // so there are twice as much observations
    if(_ambigRes)
    {
        _nobs *= 2;
    }
    
    // we need this emtpy-init due to impact-factor scheduling
    // -> has no influence on regular analysis
    _lsa_solution = ivg::Lsa();
    _lsa_solution.resize(_nobs,num_of_params,false);
    A_ptr = _lsa_solution.get_design_ptr();

    ivg::Matrix apr( num_of_params,_nobs,0.0 );
    vector<double>::iterator apr0_iter = apr.begin();
        
    // o-c vector, vector of variances as well as epochs of observations
    oc_ptr = _lsa_solution.get_obs_ptr();
    epo_ptr = _lsa_solution.get_tobs_ptr();
    wgt_ptr = _lsa_solution.get_wgt_ptr();
    aplo_ptr= _lsa_solution.get_aplo_ptr();
    // initialize index vector for origin length of observation
    vector<int> v(_nobs) ; // vector with _nobs ints
    iota (begin(v), end(v), 0); // fill from 0 to _nobs-1
    _origin_obs_idxs = v;
    
    // time for log files
    ivg::Date now;
   
    now.now();
    ivg::Matrix proc_times(1,2,0.0);
    
    // loop over scans and observations therein
    ivg::Obs* obs_ptr;
    int obs_counter = 0;
    bool apply_ion = (bool)(*_setup)["ionosphere"]["apply"];
    for(vector<ivg::Scan>::iterator scan_it = _scans.begin(); scan_it != _scans.end(); ++scan_it )
    { 
       for(int obs_i=0; obs_i < scan_it->get_nobs(); ++obs_i)
       { 
          ivg::Matrix ai( 1,num_of_params,0.0 ); // i-th row of jacobian matrix
	  
          if(obs_counter % 1000 == 0)
             log<INFO>("*** Delays calculated: ") % obs_counter % "/" % _nobs;

          // calculate delay and entries of the jacobian matrix
          obs_ptr = scan_it->get_obs_ptr(obs_i);
	  double tmp_aplo;
	  if (_phaseDelay)
	    obs_ptr->calc_delay( ai.begin(), apr0_iter, apply_ion,'p', &tmp_aplo); // i-th row of JACOBIAN MATRIX
	  else if (_SB_solution)
	    obs_ptr->calc_delay( ai.begin(), apr0_iter, apply_ion,'s', &tmp_aplo); 
	  else
	    obs_ptr->calc_delay( ai.begin(), apr0_iter, apply_ion,'g', &tmp_aplo); // i-th row of JACOBIAN MATRIX
	  
          A_ptr->set_sub(obs_counter,0,ai);
          
          ////////////////////////////////////////////////////////////////////
          if((bool)(*_setup).exists("vascc2015") && (bool)(*_setup)["vascc2015"] )
          {
            // stringstream for storing theoretical delay in VASCC2015 format
            stringstream tau_ss;
    
            // generate output for VASCC2015
            ivg::Date epoch = obs_ptr->get_epoch();
            string src, sta1, sta2;
            obs_ptr->get_source_name(src);
            obs_ptr->get_station_names(sta1,sta2);

            tau_ss << setfill('0') << setw(6) << right << obs_counter+1
            << " " << epoch.get_date_time("YYYY/MO/DD HH:MI:")
            << setfill('0') << setw(5) << right << fixed << setprecision(2) << epoch.get_double_sec()
            << " " << setfill(' ') << setw(8) << right << src
            << " " << setfill(' ') << setw(8) << right << sta1
            << " " << setfill(' ') << setw(8) << right << sta2;
            // ALL FURTHER INFORMATION FOR VASCC_FILE ADDED WITH "add_to_saveline" IN obs.cpp
            string vascc_file = string((const char*)(*_setup)["outdir"])+"/VASCC2015_"+now.get_date_time("YYYY_MO_DD_HH_MI_SS")+"_tau.txt";

            log<SAVE>(vascc_file, tau_ss.str());
          }

          // set iterator to first element of the next column of A^T
          apr0_iter += num_of_params;
          
          oc_ptr->operator()( obs_counter ) = obs_ptr->get_observed_minus_computed();
          aplo_ptr->operator()( obs_counter ) = tmp_aplo;
          //down weight group delay for ambiguity resolution
          if(_ambigRes)
              wgt_ptr->operator()( obs_counter ) = 1e-18;
          else if (_SB_solution)
	    wgt_ptr->operator()( obs_counter ) = 1.0/obs_ptr->get_obs_variance( apply_ion,false,true );
	  else
	    wgt_ptr->operator()( obs_counter ) = 1.0/obs_ptr->get_obs_variance( apply_ion,_phaseDelay );
          
          epo_ptr->operator()( obs_counter ) = scan_it->get_epoch().get_double_mjd();
          
          // calculate also single band delays
          if(_ambigRes)
          {
            obs_counter++;
	    if (_phaseDelay) {
	      
	      obs_ptr->calc_delay( ai.begin(), apr0_iter, apply_ion, 'g', &tmp_aplo );
	    }
	    else
	      {
	      
	      obs_ptr->calc_delay( ai.begin(), apr0_iter, apply_ion, 's', &tmp_aplo);
	      }
            A_ptr->set_sub(obs_counter,0,ai);
            apr0_iter += num_of_params;
            oc_ptr->operator()( obs_counter ) = obs_ptr->get_observed_minus_computed();
            wgt_ptr->operator()( obs_counter ) = 1.0/obs_ptr->get_obs_variance( apply_ion,false,!(_phaseDelay));
            epo_ptr->operator()( obs_counter ) = scan_it->get_epoch().get_double_mjd();
	    
          }
          obs_counter++;
	  (scan_it->get_source())->increase_n_obs_in_sess(1);
	 
       }
    }
    log<INFO>("*** Delays calculated: ") % obs_counter % "/" % _nobs;
    
    _aprioris = apr.transpose();
    
    // if there are clock breaks or CPWLF breaks for atmospheric parameters, include zero colomns in the apriori matrix
//    multimap<string,double> cl_br = _param_list.get_clock_breaks();
//    multimap<string,double> atm_br = _param_list.get_atmo_breaks();
    map< ivg::paramtype,multimap<string,double> > br = _param_list.get_breaks();
    if( !br[ivg::paramtype::cbr].empty() )
    {
        ivg::Matrix zeros( _aprioris.rows(),br[ivg::paramtype::cbr].size(),0.0 );
        _aprioris.append_cols(zeros);
    }
    if( !br[ivg::paramtype::atbr].empty() )
    {
        ivg::Matrix zeros( _aprioris.rows(),br[ivg::paramtype::atbr].size(),0.0 );
        _aprioris.append_cols(zeros);
    }
    
    _solution = &_lsa_solution;
    
    _eliminate_data();

    // eliminate observations in specified period
    if((bool)(*_handling).exists("handling") && (bool)((*_handling)["handling"]).exists("elim_obs_period"))
    {    
        _lsa_solution.reinit();
        _eliminate_obs_period();
    }
    if((bool)(*_setup).exists("elevation_angle_cutoff"))
    {    
        _lsa_solution.reinit();
        _apply_elevation_cutoff();
    }

    // check whether the session has a specific reference lock or not
    if((bool)(*_handling).exists("handling") && (bool)((*_handling)["handling"]).exists("ref_clock"))
    {   
        string new_ref_clock = (const char *)(*_handling)["handling"]["ref_clock"];
        
        if(!(*_setup).exists("RefClockStationList"))
            (*_setup).add("RefClockStationList",Setting::TypeString);

        (*_setup)["RefClockStationList"] = new_ref_clock;
        
        log<WARNING>("!!! Overwriting original vgosDB reference clock due to [handling][ref_clock] to: ") % new_ref_clock;
    }
    // check for session-specific baseline clocks
    if((bool)(*_handling).exists("handling") && (bool)((*_handling)["handling"]).exists("bl_clocks"))
    {   
              
      //if(!(*_setup).exists("BaselineClockList"))
      //       (*_setup).add("BaselineClockList",Setting::TypeList);
	for(int i=0;i<(*_handling)["handling"]["bl_clocks"].getLength();i++)
	  _param_list.add_bl_clock((const char*)(*_handling)["handling"]["bl_clocks"][i]);
	  //	  (*_setup)["BaselineClockList"].add(Setting::TypeString)=(const char*)(*_handling)["handling"]["bl_clocks"][i];
        
        
    }
    
    // check if stations with no apriori coordinates (=0.0) are still left
    string missing_aprioris;
    for(auto &param: _param_list)
        if(param.is_type({stax},{0,1}) &&  param.get_apriori() == 0.0)
            missing_aprioris += param.get_name()+" ";
          
    if(missing_aprioris.size()>0)
        throw runtime_error("void Session::init_vgosdb_ngs_solution(): Stations "+missing_aprioris+" without apriori values.");

    //includes _model_turbulence
    _create_weight_matrix(_lsa_solution.get_wgt_ptr());
    _lsa_solution.reinit();
               
    string type =  (*_setup)[ "solver" ];
    if( type == "LSC" )
    {
        ivg::Lsc lsc( &_lsa_solution );
        _lsc_solution = lsc;
        std::vector< ivg::Matrix > G = _create_collocation_assignment();
        _lsc_solution.set_G( G );
        _solution = &_lsc_solution;
    }
        
    // show information about network volume in megacubicmeter
    double V = _trf.calculate_network_volume()/1e18 ;
    log<INFO>("*** Scheduled networkvolume ") % V % "Mm³";
    
    log<INFO>("*** END-Session::init_vgosdb_ngs_solution with _name: ") % _name;
    
#if DEBUG_VLBI >=1
   cerr << "--- void Session::init_vgosdb_ngs_solution( )" << " : " << tim.toc() << " s " << endl;
#endif 
}
  
// ...........................................................................
void Session::init_snx_solution(string adjustment_options, ivg::Date t_0)
// ...........................................................................
{
#if DEBUG_VLBI >=1
   cerr << "+++ void Session::init_snx_solution" << endl; 
   tictoc tim;
   tim.tic();
#endif 
        /* Beside the Session+=Constructor this is one of the most important stacking / combination functions!
         Some adjustments are made manually inside this function! Sometimes you need to activate/deactivate comments!
         */
   
        log<INFO>("*** START-Session::init_snx_solution with _name: ") % _name % " and adjustments: " % adjustment_options;
        
        // clear FILE/COMMENTS from original snx-file
        _file_comment.str("");
   
        // write information from setup to FILE/COMMENT block in snx-file
        _file_comment << _create_configuration_textblock(_setup);
        
        _solution = &_neq_solution;
        // in case of debugging, this output might be interesting
//         _param_list.show("Session::init_snx_solution(): RAW DATA FROM SNX FILE BEFORE ANY FURTHER OPTERATIONS");        
       
        // get information from configfile, what to do with TRF,CRF,EOPs and parameter list concerning apriori transformation
        int option = ivg::parser::init_options(adjustment_options);
        
        // check for zero columns/rows in _neq_solution (_N and _n)
        // full session will be ignored and will be junked!
        vector<int> removed_indexes = _neq_solution.detect_erase_unparameterized();
        if(!removed_indexes.empty())
        {
            std::reverse(removed_indexes.begin(),removed_indexes.end());
            log<INFO>("*** Unparameterized source in _neq_solution. Erasing affected parameter from CRF");
            string infostr;
            for(auto &idx: removed_indexes)
            {
                ivg::Param *param = _param_list.get_param(idx);
                infostr += param->get_name()+" "+param->get_typename()+" "+std::to_string(param->get_order())+"|";
                log<WARNING>("!!! NEQ contains zeros at parameter ") % param->get_name() % " " % param->get_typename() % " " % param->get_order();
                
                // there was a time, when the whole session still was used but only the affected source was fixed
                // this is not a good idea anymore, because a single fixed source will define the datum
//                for(vector<ivg::Source>::iterator src = _crf.begin(); src != _crf.end(); ++src)
//                {
//                    if(src->get_name(ivg::srcname::iers) == _param_list.get_param(idx)->get_name())
//                    {
//                        int idx_plus = idx+1;
//                        log<DETAIL>("*** Erasing ") % src->get_name(ivg::srcname::iers) % " from CRF and _param_list at " % idx % "/" % idx_plus;
//                        _crf.remove_source(src);
//                        src--;
//                    }
//                }
//                _param_list.remove_param(idx);
            }
            
            show_vector(removed_indexes);
            throw runtime_error("WARNING: Critical session contains #"+std::to_string(removed_indexes.size())+" zeros within NEQ: "+infostr);
        }
                
        // if necessary transform from two cpwlf to offset+rate
        // right now only transformation using two cpwlfs is implemented, not three!
        // TODO: implement transformation of three cpwlfs to offset / rate
        _param_list.transform_cpwlf2offsetrate(_neq_solution, t_0);
        
        // get complete apriori vector before any adjustments
        ivg::Matrix old_aprioris = _param_list.extract_apriori_vector();
        
        //adjusting CRF,TRF,EOPs to parameter list, depends on options set in configfile
        // TRF2APR -> e.g. vtrf2014 to aprioris in param_list
        // CRF2APR -> e.g. icrf2 to aprioris in param_list
        // EOP2APR -> e.g. c04 to aprioris im param_list
        // APR2TRF / EST2TRF -> param_list aprioris/estimates to _trf
        // APR2CRF / EST2CRF -> param_list aprioris/estimates to _crf
        // APR2EOP / EST2EOP -> param_list aprioris/estimates to _eop_series
        _adjust_data_storage(option);
         
        // we only need to perform a aprioi transformation if external TRF,CRF,EOPs aprioris have been used (e.g. vtrf2014, icrf2, c04)
        if( option & TRF2APR || option & CRF2APR || option & EOP2APR )
        {
            log<INFO>("*** Running apriori transformation");
            
            ivg::Matrix new_aprioris = _param_list.extract_apriori_vector();  
	    
            //show old, new and difference of aprioris
            // just uncomment this block and you can see the differences between old and new aprioris
            //ivg::Matrix tmp = old_aprioris;
            //tmp.append_cols(new_aprioris);
            //ivg::Matrix diff = new_aprioris - old_aprioris;
            //tmp.append_cols(diff);
            //tmp.show(15);

            // now _neq_solution (_n) need to be adjusted (apriori transformation)
            _neq_solution.apriori_transformation(new_aprioris, old_aprioris);
           
            // in case of global modus, we use a specific reference epoch
            // and might want to introduce velocities for the stations
            string global_str = (const char*)(*_setup)["PARAMS"]["stations"][0]["stacking"][ "type" ];
	    if ((*_setup).exists("global"))
	      {
		string breakfile=(*_setup)["global"]["break_file"];
		if ( global_str == "global" && breakfile !="")
		  _param_list.set_breaks(breakfile);
	      }
            if( global_str == "global" && (bool)(*_setup)["PARAMS"]["stations"][0]["stacking"][ "insert_velocities" ])
            {
	      
                ivg::Date ref_epoch;
		
                ref_epoch.set_decimal_date((double)(*_setup)["PARAMS"]["stations"][0]["stacking"][ "ref_epoch" ]);
		 
		_param_list.insert_station_velocities(_neq_solution,ref_epoch, _trf);
		
            }
        }
        
        // calculate correct btPb 
        // this is mandatory to get correct standard deviations within the 
        // right now it's deactivated because fix_eops_sources is missing after some code changes
        // we need the information, where (which index) the eops and sources have, in order to fix them.
        // we can get this informaton within adjust_data_storage
        if(0)
        {
//            // -> generate a solution with fixed EOPs and SOURCES and NNR/NNT for STATIONS
//            log<INFO>("*** Calculating correct btPb");
//
//            // copy solution and parameterlist to avoid disorder in original data
//            ivg::Ls_neq _btPb_solution = _neq_solution;
//            ivg::Param_list _btPb_list = _param_list;
//
//            // fix EOPs and SOURCES and remove them from the list
//            _btPb_solution.fix_param( fix_eops_sources );
//            for( int i=fix_eops_sources.size()-1; i>=0; --i )
//                _btPb_list.remove_param(fix_eops_sources.at(i));
//
//            // add NNR/NNT constraint for all STATIONS
//            ivg::Matrix w( 6,1,1.0/pow( 1e-6, 2.0 ) );
//            vector<string> station_names = _trf.get_station_names( ivg::staname::ivs_name );
//            ivg::Matrix B = _btPb_list.generate_station_nnr_nnt_B(station_names, w);
//            _btPb_solution.update_constraints( B,w );
//
//            // solve the system to obtain parametervector _x
//            _btPb_solution.solve(ivg::solutiontype::neq_lu, false);
//
//            // calculate correct btPb based on: _simga0_post, _nobs, _nparam, _x, _N, _n 
//            // -> sets correct _btPb
//            _btPb_solution.calc_correct_btPb();
//            _neq_solution.set_btPb(_btPb_solution.get_btPb());
        }
        
        if(1) // sometimes we deactive it manually
        {
        // epoch transformation based on t_0 
        old_aprioris = _param_list.extract_apriori_vector();
	
        // getindexes of XPO, YPO, UT1 and XPOR, YPOR, LOD
        vector<ivg::paramtype> types = {xpo,ypo,ut1};
        vector<int> idxO_vec, idxR_vec;
        for(std::vector<ivg::Param>::iterator param = _param_list.begin(); param != _param_list.end(); ++param)
        {
            int idx = param-_param_list.begin();
            for(auto &type: types)
            {
                if(param->is_type({type},{0,1}))
                {
                    if(param->get_order() == 0)
                        idxO_vec.push_back(idx);
                    else if(param->get_order() == 1)
                        idxR_vec.push_back(idx);
                }
            }
        }
	
        if(idxO_vec.size() != idxR_vec.size())
            throw runtime_error("void Session::init_snx_solution(): EOPs and rates does not correspond for epoch transformation.");

        if(idxO_vec.size() != 3)
            throw runtime_error("void Session::init_snx_solution(): Unexpected length of EOP vector. No epoch transformation possible.");

        for(int i=0;i<idxO_vec.size(); i++)
        {
            int idxO = idxO_vec.at(i);
            int idxR = idxR_vec.at(i);

            double delta_t = _param_list.get_param(idxO)->get_epoch().get_double_mjd() - t_0.get_double_mjd();
            log<INFO>("*** Running epoch transformation over ") % delta_t % " days";
            // check if epoch transformation should be longer than half a day
            if(delta_t > 0.5)
                throw runtime_error("WARNING: Epoch transformation too extended");
            
            _neq_solution.epoch_transformation( idxO, idxR, delta_t );

            double new_apri_idxO =  _param_list.get_param(idxO)->get_apriori() - delta_t * _param_list.get_param(idxR)->get_apriori();
            _param_list.get_param(idxO)->set_apriori(new_apri_idxO);
            _param_list.get_param(idxO)->set_epoch(t_0);
            _param_list.get_param(idxR)->set_epoch(t_0);
        }

        ivg::Matrix new_aprioris = _param_list.extract_apriori_vector();  

         // show some output
//            ivg::Matrix tmp = old_aprioris;
//            tmp.append_cols(new_aprioris);
//            ivg::Matrix diff = new_aprioris - old_aprioris;
//            tmp.append_cols(diff);
//            tmp.show(15);            
        }
           
        // nutation and ut1 leap seconds correction are done without consideration of apriori transformation
        // setting nutation to zero (X,Y as well as LN,OB)
        for( std::vector<ivg::Param>::iterator param = _param_list.begin(); param != _param_list.end(); ++param )
        {
            if( option & NUT2ZER )
            {
                if(param->is_type({nutx, nuty, nutln, nutob},{0}))
                {
                    param->set_apriori(0.0);
                    log<DETAIL>("*** EOPs _param_list [nut (x,y,ln,ob)] setting to 0.0");
                }
            }
        }
        
        // adjust leap seconds in case of UT / UT1
        for( std::vector<ivg::Param>::iterator param = _param_list.begin(); param != _param_list.end(); ++param )
        {
            if(param->is_type({ut1},{0}))
            {
                double apri = param->get_apriori();
                double esti = param->get_estimate();
                if( apri > -1500.0/ivg::param_unit_fac.at(ut1) && apri < 1500.0/ivg::param_unit_fac.at(ut1) )
                {
                    double new_apri = apri - (param->get_epoch().get_leap_sec()*1000)/ivg::param_unit_fac.at(ut1);
                    double new_esti = esti - (param->get_epoch().get_leap_sec()*1000)/ivg::param_unit_fac.at(ut1);
                    param->set_apriori(new_apri);
                    param->set_estimate(new_esti);
                    log<DETAIL>("*** EOP _eops --> _param_list[") % param->get_typename() % " " % param->get_order() % "] Correcting for leap seconds: " % apri % " to "  % new_apri;
                    log<DETAIL>("*** EOP _eops --> _param_list[") % param->get_typename() % " " % param->get_order() % "] Correcting for leap seconds: " % esti % " to " % new_esti;
                }
            }
        }
        
        //only if session specific handling exists, use it
        if((bool)(*_handling).exists("handling"))
        {
            //getting new handling if EOPs of session need to be fixed for example
            vector<string> eops = {"pm","ut1","nut"};
            for(int i=0; i<eops.size(); i++){
                string new_handling = (*_handling)["handling"];
                (*_setup)["PARAMS"].lookup(eops.at(i))[0]["handling"] = new_handling;
            }
        }
        
        // setting nutation to same epoch (X,Y as well as LN,OB) in order to become combined
        vector<int> nut_idx = _param_list.get_indexes({nutx, nuty, nutln, nutob},"EOP");
        for(int i=0; i<nut_idx.size(); i++)
            if( nut_idx.at(i) != -1 )
                _param_list.get_param(nut_idx.at(i))->set_epoch(t_0);
              
        // detect if clockbreaks are existent (in case of ivg::ascot sinex files)
        for(std::vector<ivg::Param>::iterator param = _param_list.begin(); param != _param_list.end(); ++param)
            if(param->is_type({ivg::paramtype::cbr},{0}))
                param->set_reduce_flag(true);
            
//        _neq_solution.set_scales(scales);
//        _neq_solution.scale_system();
        
      log<INFO>("*** END-Session::init_snx_solution with _name: ") % _name;
      
#if DEBUG_VLBI >=1
   cerr << "--- void Session::init_snx_solution( )" << " : " << tim.toc() << " s " << endl;
#endif 
} 

// ...........................................................................
void Session::modify_parameterization()
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ void Session::modify_parameterization() " << endl;
   tictoc tim;
   tim.tic();
#endif
   
    // in case of processing vgosdb / ngs
    if( _type.find("snx") == string::npos )
    {
        ivg::Matrix obs_mjd = get_obs_epochs();
        _param_list.modify_parameterization(_name, (*_setup),_trf,_crf, (*_solution), _aprioris, obs_mjd );
    }
    else
    // in case of using snx and global_main
    {
        ivg::Matrix dummy(1,_param_list.size(),0.0); 
       _param_list.modify_parameterization(_name, (*_setup),_trf,_crf, (*_solution), dummy, dummy ); 
    }  
   
#ifdef DEBUG_VLBI
   cerr << "--- void Session::modify_parameterization() "
        << ": " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Session::reduce_and_constrain()
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ void Session::reduce_and_constrain() " << endl;
   tictoc tim;
   tim.tic();
#endif   
    tictoc tim;
    tim.tic();
    _param_list.create_constraint_conditions( (*_solution) );
    
    // find params w/o observations
    if((bool)(*_setup).exists("remove_undefined_params") &&  
       (bool)(*_setup)["remove_undefined_params"])
    {
      
        vector<int> id_zero = _lsa_solution.find_undefined_param();
	
        // and eliminate them
        if(id_zero.size()>0)
        {
	  
            for( int i=id_zero.size()-1; i>=0; --i )
            {
                log<WARNING>("!!! Fixing parameter [") 
                   % _param_list.get_param(id_zero.at( i ))->get_name() % " "
                   % _param_list.get_param(id_zero.at( i ))->get_typename() % " " 
                   % _param_list.get_param(id_zero.at( i ))->get_order() 
                   % "] --> no observations!!!!";
               _param_list.remove_param(id_zero.at( i ));
            }
	    
	    
            _lsa_solution.fix_param(id_zero);
        }
    }
    
    _param_list.reduce_params( (*_solution));
   
    // if session_type from configfile != snx
    if( _type.find("snx") == string::npos )
    {
        string type =  (*_setup)[ "solver" ];
        if( type == "LSC" )
        {
            _lsc_solution.get_lsa_ptr()->build_neq();
        }   
        else
            _lsa_solution.build_neq();
    }
    
         
#ifdef DEBUG_VLBI
   cerr << "--- void Session::reduce_and_constrain() "
        << ": " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
string Session::solve( bool use_only_indep_bl)
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ string Session::solve() " << endl;
   tictoc tim;
   tim.tic();
#endif
   
   // only allow nnr/nnt datum definition once per session (avoids overconstraining while strenghten weak sources in global solution)
   if(_nnr_nnt_set == false)
   {
    _param_list.create_nnr_nnt_equations((*_setup),_trf,_crf, (*_solution));
    _nnr_nnt_set = true;
   }
    if ((*_setup).exists("global"))
     {
       string constrfile=(*_setup)["global"]["velocity_constraints"];
       if ( constrfile !="")
	 _param_list.create_vel_constr(constrfile, (*_setup),  (*_solution) );
     }
    _param_list.create_common_clock_equations( (*_setup),(*_solution) );
            
    // OUTLIER HANDLING from config file:
    // set number iterations to zero, thereshold for detecting gross errors to 1ns 
    // and the quantile of statistical distribution to 3 sigma
    int iterations = 0;
    double thresh = 1.0e-9;
    double quantile = 3.0;

    bool pre_cond = true;
    if((bool)(*_setup).exists("pre-conditioning"))
       pre_cond = (bool)(*_setup)[ "pre-conditioning" ];
   
    // using a simple outlier test: threshold for gross errors and a quantile (inverse cdf) 
    // Only do it if it is not snx session_type
    if(  (bool)(*_setup)["outliers"]["detect"] && _type.find("snx") == string::npos )
    {
       std::string test_method = (*_setup)["outliers"]["data_snooping"];

       if( test_method == "simple" )
       {
          iterations = (int)get_list_element( (*_setup)["test_method"],(*_setup)["outliers"]["data_snooping"])[1]["iterations"];
          thresh = (double)get_list_element( (*_setup)["test_method"],(*_setup)["outliers"]["data_snooping"])[1]["threshold"];
          quantile = (double)get_list_element( (*_setup)["test_method"],(*_setup)["outliers"]["data_snooping"])[1]["quantile"];
       }
    }
   
    std::vector<int> at_idx;       
    string type =  (*_setup)[ "solver" ];       
    double tikh_lambda = 0.0;
    if((bool)(*_handling).exists("handling") && (bool)((*_handling)["handling"]).exists("tikhonov_lambda"))
        tikh_lambda = (double)(*_handling)["handling"]["tikhonov_lambda"];

    std::cout << "Type: " << _type << endl;
    // if session_type from configfile != snx
    if( _type.find("snx") == string::npos )
    {
        if(  (bool)(*_setup)["export_gmm"] )
           write_backend_files( (const char *)(*_setup)["outdir"] );
        
        string type =  (*_setup)[ "solver" ];
        if( type == "LSM" ){
            if(use_only_indep_bl)
            {
                // problems may occur if baseline dependend parameters like baseline clocks are estimated
                
                log<WARNING> ("!!! only using independent baselines to calculate parameters" );
                

                
                // copy current lsa object (including designmatix, observation vector, ...)
                ivg::Lsa lsa_tmp = _lsa_solution;

                // get all row indices of A that do not correspond to the reference clock
                std::vector<int> rows = get_star_formation_indices(string( (const char*)(*_setup)["RefClockStationList"] ));
                
                // delete those indices and solve the neq
                lsa_tmp.remove_observations( rows );
                lsa_tmp.reinit();
                lsa_tmp.solve_neq( iterations, thresh, quantile, pre_cond, tikh_lambda );
                
                // use parameters to compute residuals for all observation
                _lsa_solution.repalce_param(lsa_tmp.get_parameters(), lsa_tmp.get_vcm());
            }
            else{
	     
                _lsa_solution.solve_neq( iterations, thresh, quantile, pre_cond, tikh_lambda );
	
            }
        }
           
        else if( type == "ICLS" )
        {
            _lsa_solution.solve_neq( iterations, thresh, quantile, pre_cond );

            // check if there are negative ZWDs
            if( _exist_negative_zwd() )
            {                           
                // create inequality constraints and initialize ICLS object
                _create_inequality_constraints( at_idx );

                int MC = (int)(*_setup)["MC"];
                _icls_solution.estimate_quality( MC );
                _icls_solution.solve_with_active_set();

                _solution = &_icls_solution; 
            }

        }
        else if( type == "LSC" )
        {
            // solve deterministic normal equation system
            _lsc_solution.get_lsa_ptr()->solve_neq( iterations, thresh, quantile, pre_cond );
            
            // add stochastic parameters
//            _param_list.add_stoch_param();
                        
            // create stochastic variance-covariance matrix by means of a
            // suitable correlation function and solve the least squares
            // collocation method
            ivg::Matrix Qy = _create_correlation_fct();
            if( Qy.rows() != 0 && Qy.cols() != 0 )
            {
                _lsc_solution.calc_stoch_VCM( Qy );
                _lsc_solution.solve();
            }
        }
//        else
//           _lsa_solution.solve( iterations, thresh, quantile, pre_cond );
    }
    // if session_type = snx
    else
    {
        // in case of new_solve in config file, solve with inv()
        if((bool)(*_setup).exists("new_solve"))
            _neq_solution.new_solve();
        else
            _neq_solution.solve(ivg::solutiontype::neq_lu, false);
    }
    double perc_outliers = 0.0;
    int n = 0;
    // OUTLIER ELIMINATION (DATA SNOOPING): use BAARDA test if the variance factor is known; use POPE if not
    // Only do it if it is not snx session_type
    if(  (bool)(*_setup)["outliers"]["detect"]  && _type.find("snx") == string::npos )
    {
       std::string test_method = (*_setup)["outliers"]["data_snooping"];

       if( test_method == "baarda" || test_method == "pope" )
       {
            double alpha = (double)get_list_element( (*_setup)["test_method"],(*_setup)["outliers"]["data_snooping"])[1]["significance_level"];
          
            if( type == "LSM" )
            {
                _lsa_solution.data_snooping( alpha, test_method, n, perc_outliers, _exp_perc_out);
                
                if( (bool)(*_setup)["outliers"]["restoration"][0] && n > 0 )
                {
                    n = _lsa_solution.restore_outliers( (double)(*_setup)["outliers"]["restoration"][1] );
                    
                    if( n > 0)
                        _lsa_solution.data_snooping( alpha, test_method, n, perc_outliers, _exp_perc_out );                    
                }
            }
            else if( type == "LSC" )
            {
                _lsc_solution.data_snooping( alpha, test_method, n, perc_outliers, _exp_perc_out );
//                To-Do: outlier restauration for LSC!                
            }            
            else if( type == "ICLS" )
                throw runtime_error( "!!! outlier elimination not possible with an ICLS solution !!!" );
       }
    }

    ivg::Matrix x = _solution->get_parameters();  
    
    //save Sxx
    if( (bool)(*_setup)["export_gmm"] ){
        (_solution->get_vcm()/_solution->calc_posterior_vfac()).save_bin((const char *)(*_setup)["outdir"]+_name+"gmm"+"XY"+"_Qxx.dat");
    }
    
      
    double vfac, wrms, rms;
    if( type == "LSM" )
    {
        // set standard deviations in param_list
        ivg::Matrix Sxx = _solution->get_vcm();
        _param_list.set_estimates( x, Sxx.diag().sqrt() );
        
        // show results: estimates, VFAC, (W)RMS
        _param_list.show_estimates();
        
        vfac = sqrt(_solution->calc_posterior_vfac());
        wrms = _solution->calc_wrms();
        rms = _solution->calc_rms();
        _solution->rm_nnt_nnr_constraints();
        log<RESULT>("*** Database: ") % _name % " VFAC: " % vfac  % " WRMS: " % wrms % " RMS: " % rms;
    }
    else if( type == "ICLS")
    {

        if( _exist_negative_zwd() )
        {         
            // set highest probability density (HPD) intervals instead of
            // standard deviations in param_list
            ivg::Matrix Sxx = _icls_solution.get_confidence_limits().absD();
            // To-Do: confidence interval in param_list and sinex file
            _param_list.set_estimates( x, Sxx(":",1) );        
                       
            // scale confidece limits (highest probability density)
//            ivg::Matrix units( x.rows(),1, 0.0 );
//            int counter = 0;
//            for( std::vector<ivg::Param>::iterator it = _param_list.begin(); it != _param_list.end(); ++it )
//            {            
//                units( counter, 0 ) = ivg::param_unit_fac.at( it->get_type() );
//                counter++;
//            } 
//            Sxx.set_sub( 0,0, Sxx.get_col(0).mult_elem( units ) );
//            Sxx.set_sub( 0,1, Sxx.get_col(1).mult_elem( units ) );          
//  
//            string dir = (*_setup)[ "outdir" ];
//            std::string hpd = dir+_name+"_hpd";
//            Sxx.save_bin( hpd );
            
            // get Lagrange multipliers to search for the active constraints
            // which should be written into the sinex file (see below)
            ivg::Matrix k_idx = _icls_solution.get_lagrange_mult();             
            ivg::Matrix cnstr_idx2( x.rows(), 1, 0.0 );
            cnstr_idx2.setIdx( at_idx, k_idx );
            std::vector<int> cnstr_idx = cnstr_idx2.find_idx( gt, 0.0 );
          
            // show results: estimates, VFAC, (W)RMS
            _param_list.show_estimates();

            vfac = sqrt( _icls_solution.get_lsa_ptr()->calc_posterior_vfac() );
            wrms = _icls_solution.get_lsa_ptr()->calc_wrms();
            rms = _icls_solution.get_lsa_ptr()->calc_rms();
            _icls_solution.set_statistics( wrms, rms, vfac*vfac );

            log<RESULT>("*** Database: ") % _name % " VFAC: " % vfac % " WRMS: " % wrms % " RMS: " % rms;
            
            // if there is an active inequality constraint, write parameter
            // name and epoch to sinex file (FILE/COMMENT block)                     
            for( int k = 0; k < cnstr_idx.size(); ++k )
            { 
                _file_comment << " Inequality Constraint " << k+1 << ": " << _name << " " 
                              << cnstr_idx.at(k)+1 << " " << setprecision(15) 
                              << _param_list.get_param(cnstr_idx.at(k))->get_epoch().get_double_mjd() << " "
                              << _param_list.get_param(cnstr_idx.at(k))->get_name() << endl;                            
            }
        }
        else
        {
            ivg::Matrix Sxx = _solution->get_vcm();
            _param_list.set_estimates( x, Sxx.diag().sqrt() );

            // show results: estimates, VFAC, (W)RMS
            _param_list.show_estimates();
            
            vfac = sqrt( _solution->calc_posterior_vfac() );
            wrms = _solution->calc_wrms();
            rms = _solution->calc_rms();
            _icls_solution.set_statistics( wrms, rms, vfac*vfac );
            
            log<RESULT>("*** Database: ") % _name % " VFAC: " % vfac % " WRMS: " % wrms  % " RMS: " % rms;            
        }
    }    
    else if( type == "LSC" )
    {
        ivg::Matrix x = _lsc_solution.get_lsa_ptr()->get_parameters();
        ivg::Matrix Sxx = _lsc_solution.get_lsa_ptr()->get_vcm();
        ivg::Matrix y = _lsc_solution.get_stoch_params();
        ivg::Matrix Syy = _lsc_solution.get_stoch_vcm();
        
//        // add stochastic parameters to the (original, deterministic) parameter
//        // vector and corresponding variance-covariance matrix
//        if( y.rows() > 0 )
//        {
//            x.append_rows(y);
//            ivg::Matrix tmp( Sxx.rows()+Syy.rows(),Sxx.rows()+Syy.rows(),0.0);
//            tmp.set_sub(0,0,Sxx);
//            tmp.set_sub( Sxx.rows(),Sxx.rows(),Syy );
//            Sxx = tmp;
//        
//            // update number of params and variance-covariance matrix in
//            // ivg::Ls_solution for SINEX output files
//            _solution->set_nparam(x.rows());
//        }
//        _solution->set_vcm(Sxx);
        
        // set standard deviations in param_list        
        _param_list.set_estimates( x, Sxx.diag().sqrt() );
                
        // set preciction vector and corresponding standard deviation
        _param_list.set_estimates( y, Syy.diag().sqrt(), true );
        
        // show results: estimates, VFAC, (W)RMS
        _param_list.show_estimates();
        _param_list.show_estimates( true );
        
        vfac = sqrt( _lsc_solution.calc_posterior_vfac() );
        wrms = _lsc_solution.calc_wrms();
        rms = _lsc_solution.calc_rms();
        
        log<RESULT>("*** Database: ") % _name % " VFAC: " % vfac  % " WRMS: " % wrms % " RMS: " % rms;
    }    
    
    // create session info for logfile
    stringstream sess_info;
    if(std::isnan(vfac))
        sess_info << "VFAC:       NaN ";
    else
        sess_info << "VFAC: " << setw(6) << setfill(' ') << scientific << setprecision(3) << vfac << " ";
    
    sess_info << "WRMS: " << setw(6) << setfill(' ') << scientific << setprecision(3) << wrms << " ";
    sess_info << "OUT: " << setw(5) << setfill(' ') << fixed << setprecision(2) << perc_outliers << "%";
        
#ifdef DEBUG_VLBI
   cerr << "--- string Session::solve() "
        << ": " << tim.toc() << " s " << endl;
#endif
   
   return sess_info.str();
}
// ...........................................................................
vector<string> Session::check_stations_estimates()
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ vector<string> Session::check_stations_estimates() " << endl;
   tictoc tim;
   tim.tic();
#endif
   
   map<double, string > bad_stas;
   for(auto &sta: _trf)
   {
       // get corresponding X,Y,Z from station
       vector<int> idx = _param_list.get_indexes({stax, stay, staz}, sta.get_name(ivg::staname::ivs_name));
       
       if(idx.size() != 3 && idx.size() != 0)
            throw runtime_error( "vector<string> Session::check_stations_estimates(): Unexpected number of station coordinates for "+sta.get_name(ivg::staname::ivs_name));
       else if(idx.at(0) != -1 && idx.at(1) != -1 && idx.at(2) != -1)
       {
           ivg::Matrix xyz(3,1,0.0);
           xyz(0,0) = _param_list.get_param(idx.at(0))->get_estimate()*ivg::param_unit_fac.at(_param_list.get_param(idx.at(0))->get_type());
           xyz(1,0) = _param_list.get_param(idx.at(1))->get_estimate()*ivg::param_unit_fac.at(_param_list.get_param(idx.at(1))->get_type());
           xyz(2,0) = _param_list.get_param(idx.at(2))->get_estimate()*ivg::param_unit_fac.at(_param_list.get_param(idx.at(2))->get_type());
           
           // tranform from geocentric to topocentric
           ivg::Matrix ren = sta.form_topo2geo().transpose() * xyz;
           
           // check if estimates in up, eath, north fit the limitations or not
           if(ren(0,0) > (double)(*_setup)["station_threshold"]["vertical"] || 
              ren(1,0) > (double)(*_setup)["station_threshold"]["horizontal"] ||
              ren(2,0) > (double)(*_setup)["station_threshold"]["horizontal"] )
           {
               stringstream ss;
               ss << sta.get_name(ivg::staname::ivs_name) << "[";
               ss << scientific << setprecision(1) << ren(0,0) << ",";
               ss << scientific << setprecision(1) << ren(1,0) << ",";
               ss << scientific << setprecision(1) << ren(2,0) << "]";
               
               double distance = sqrt(pow(ren(0,0),2) + pow(ren(1,0),2) + pow(ren(2,0),2));
               bad_stas[distance] = ss.str();
           }
       }
   }
   
   // sort bad stations by increasing "estimate-vector-distance" order
   vector<string> bad_stations;
   for(auto &tmp: bad_stas)
       bad_stations.push_back(tmp.second);
   
#ifdef DEBUG_VLBI
   cerr << "--- vector<string> Session::check_stations_estimates() "
        << ": " << tim.toc() << " s " << endl;
#endif
   return(bad_stations);
}
// ...........................................................................
vector<double> Session::_adjust_data_storage(int opt)
// ...........................................................................
{
#ifdef DEBUG_VLBI
   cerr << "+++ void Session::_adjust_data_storage(uint32_t opt) " << endl;
   tictoc tim;
   tim.tic();
#endif
   
    /*
     This is again one of the most important stacking/combination functions!
     Depending on the selected options, an external apriori information is used (e.g. VTRF2014 and ICRF2) or the given aprioris from the snx files are used!
     A lot of book keeping is done here!
     */
   
        // vector containing scale-factors for all parameters
        // scaling was done to avoid numerical problems
        // there have been different approaches, using different scales for different parameter groups or even depending on the amount of stacking
        // at the moment only for sources scale factors are used
        vector<double> scales(_param_list.size(),1.0);

        // save indexes from all eops and sources for correct btPb calculation
        vector<int> fix_eops_sources;

        // to be able to verify if all aprioris/estimates have been adjusted
        int adjustments = 0;
        
        struct disconti
        {
            ivg::Date refepoch;
            ivg::Matrix pos0;
            ivg::Matrix vel0;
            vector<ivg::Date> epochs;
        };
        
        map<ivg::Analysis_station *, disconti> AC_discontis;
        
        // adjustments concerning the !TRF!
        if(opt & EST2TRF || opt & APR2TRF || opt & TRF2APR)
        {    
            // go through ALL parameter
            for(vector<ivg::Param>::iterator param_iter = _param_list.begin(); param_iter!= _param_list.end(); ++param_iter)
            {
                // only if the paramter is a X-coordinate
                if(param_iter->get_type() == ivg::paramtype::stax && param_iter->get_order() == 0)
                {
                    //only if X,Y,Z are in ascending paramter-order
                    if((param_iter+1)->get_type() == ivg::paramtype::stay && (param_iter+2)->get_type() == ivg::paramtype::staz )
                    {
                        // only if epochs of X,Y,Z are equal
                        if( (param_iter)->get_epoch() == (param_iter+1)->get_epoch() && (param_iter+1)->get_epoch() == (param_iter+2)->get_epoch() )
                        {
                            // then we have the three X,Y,Z parameters of a station
                            ivg::Analysis_station * sta_iter;
                            // now we need the corresponding analysis_station (sta_iter) to these three param_iter
                            if(!_trf.get_station(&sta_iter,param_iter->get_name(), ivg::staname::description))
                                throw runtime_error( "vector<double> Session::_adjust_data_storage(int opt): Corresponding Analysis_station NOT found in _trf: "+param_iter->get_name());
                            
                            // we need to check if velocities are existent in _param_list
                            // in general this shouldn't be the case
                            // it's only the case if we use ITRF2014.SNX as apriori TRF, because this function is then also used
                            if( param_iter->get_name() == (param_iter+3)->get_name() &&
                               (param_iter+3)->get_type() == ivg::paramtype::stax && (param_iter+3)->get_order() == 1 )
                            {
                                // adjust internal _trf to values from parameter-vector from snx file
                                if(opt & EST2TRF || opt & APR2TRF)
                                {
                                    ivg::Matrix pos0_tmp, vel0_tmp;
                                    if(opt & EST2TRF)
                                    {
                                        pos0_tmp = ivg::Matrix(vector<double>{param_iter->get_estimate(),(param_iter+1)->get_estimate(),(param_iter+2)->get_estimate()});
                                        vel0_tmp = ivg::Matrix(vector<double>{(param_iter+3)->get_estimate(),(param_iter+4)->get_estimate(),(param_iter+5)->get_estimate()});
                                    }
                                    else
                                    {
                                        pos0_tmp = ivg::Matrix(vector<double>{param_iter->get_apriori(),(param_iter+1)->get_apriori(),(param_iter+2)->get_apriori()});
                                        vel0_tmp = ivg::Matrix(vector<double>{(param_iter+3)->get_apriori(),(param_iter+4)->get_apriori(),(param_iter+5)->get_apriori()});
                                    }
                                    
                                    AC_discontis[sta_iter].refepoch =  param_iter->get_epoch();
                                    AC_discontis[sta_iter].pos0.append_cols(pos0_tmp);
                                    AC_discontis[sta_iter].vel0.append_cols(vel0_tmp);

                                    log<DETAIL>("*** ") % sta_iter->get_name(ivg::staname::ivs_name) % " _param_list [stax,stay,staz] and [velx,vely,velz]--> _trf";
                                    
                                }
                                else if(opt & TRF2APR )
                                    throw runtime_error( "vector<double> Session::_adjust_data_storage(int opt): Velocities and TRF2APR - NOT IMPLEMENTED YET");
                                
                              adjustments += 6;
                            }
                            // in case of NO velocity-parameter
                            else
                            {
                                // adjust internal _trf to values from parameter-vector from snx file
                                if(opt & EST2TRF || opt & APR2TRF)
                                {
                                    double snx_x, snx_y, snx_z;
                                    if(opt & EST2TRF)
                                    {
                                        snx_x = param_iter->get_estimate();
                                        snx_y = (param_iter+1)->get_estimate();
                                        snx_z = (param_iter+2)->get_estimate();
                                    }
                                    else if(opt & APR2TRF )
                                    {
                                        snx_x = param_iter->get_apriori();
                                        snx_y = (param_iter+1)->get_apriori();
                                        snx_z = (param_iter+2)->get_apriori();
                                    }

                                    ivg::Matrix snx_xyz(vector<double>({snx_x, snx_y, snx_z}));
                                    // adjusting xyz from internal ascot units to [m]
                                    snx_xyz *= ivg::param_unit_fac.at(param_iter->get_type());
                                    ivg::Date epoch = param_iter->get_epoch();

                                    ivg::Matrix ssc_xyz = sta_iter->calc_xyz(epoch);
                                    ivg::Matrix xyz0 = sta_iter->get_xyz0();
                                    ivg::Matrix delta = snx_xyz - ssc_xyz;

                                    // reproduce fitting matrix for pos0
                                    ivg::Matrix delta_rep(1,1,1.0);
                                    delta_rep.repmat(delta,1,xyz0.cols());

                                    // set corrected pos0 matrix
                                    ivg::Matrix new_xyz = delta_rep + xyz0;
                                    sta_iter->set_xyz0(new_xyz);

                                    log<DETAIL>("*** ") % sta_iter->get_name(ivg::staname::ivs_name) % " _param_list [stax,stay,staz] --> _trf";

                                }
                                // now the other way around!
                                // adjust parameter-vector from snx file to values from internal _trf
                                else if(opt & TRF2APR )
                                {
                                    ivg::Date epoch = param_iter->get_epoch();

                                    // in case of global modus, we use a specific ref_epoch to be transformed on
                                    // if we run a global solution, we want a specific reference epoch to which we refer
                                    string global_str = (const char*)(*_setup)["PARAMS"]["stations"][0]["stacking"][ "type" ];
                                    ivg::Matrix catalog_xyz;
                                    if( global_str == "global")
                                    {
                                        ivg::Date tmp_ref_epoch;
                                        tmp_ref_epoch.set_decimal_date((double)(*_setup)["PARAMS"]["stations"][0]["stacking"][ "ref_epoch" ]);
                                        catalog_xyz = sta_iter->calc_xyz(tmp_ref_epoch);
					catalog_xyz = sta_iter->calc_xyz(epoch);
                                    }
                                    else
                                        catalog_xyz = sta_iter->calc_xyz(epoch);
                                    
                                    // in case of itrf2014 and if PSD is set to true, add the PSD
                                    if( (bool)get_list_element((*_setup)["stadisp"], "PSD" )[1] )
                                    {
                                        ivg::Matrix d_tmp( 3,1,0.0 );
                                        d_tmp = sta_iter->calc_psd_displacement( epoch );
                                        catalog_xyz += d_tmp;
                                    }
                                    
                                    if(catalog_xyz(0) == 0.0 || catalog_xyz(1) == 0.0 || catalog_xyz(2) == 0.0)
                                        throw runtime_error( "ERROR: No xyz-position for station "+sta_iter->get_name(ivg::staname::ivs_name)+". No apriori transformation possible.");
                                        
                                    // adjusting xyz in [m] to internal ascot units
                                    catalog_xyz *= 1.0/ivg::param_unit_fac.at(param_iter->get_type());
                                    
                                    // this needs to reconsidered...
                                    double max_diff = 0.5; // [m] maximal difference in x,y,z for aprioris (of NOT global modus)
                                    double x_diff = abs(catalog_xyz(0)- param_iter->get_apriori()) *  ivg::param_unit_fac.at(param_iter->get_type());
                                    double y_diff = abs(catalog_xyz(1)- (param_iter+1)->get_apriori()) *  ivg::param_unit_fac.at(param_iter->get_type());
                                    double z_diff = abs(catalog_xyz(2)- (param_iter+2)->get_apriori()) *  ivg::param_unit_fac.at(param_iter->get_type());
                                    
                                    if( global_str != "global" && (x_diff > max_diff || y_diff > max_diff || z_diff > max_diff) )
                                    {
                                        stringstream ss;
                                        ss <<  "!!! Apriori differences quite big for ";
                                        ss << sta_iter->get_name(ivg::staname::ivs_name) << "[";
                                        ss << fixed << setprecision(1) << x_diff << "," << y_diff << "," << z_diff << " m]";
                                        log<WARNING>(ss.str());
                                    }
                                          
                                    param_iter->set_apriori(catalog_xyz(0));
                                    (param_iter+1)->set_apriori(catalog_xyz(1));
                                    (param_iter+2)->set_apriori(catalog_xyz(2));

                                    log<DETAIL>("*** ") % sta_iter->get_name(ivg::staname::ivs_name) % " _trf -->  _param_list [stax,stay,staz]";
                                }

                                adjustments += 3;
                            }
                        }
                        else
                            throw runtime_error( "vector<double> Session::_adjust_data_storage(int opt): Epochs of X,Y,Z of station "+param_iter->get_name()+" not equal.");
                    }
                    else
                        throw runtime_error( "vector<double> Session::_adjust_data_storage(int opt): X,Y,Z of station "+param_iter->get_name()+" not in correct order.");
                }
            }
        }
        
        // in case of snx files with VEL and option EST2TRF or APR2TRF, complete pos0 and vel0 matrices including epochs are saved in AC_discontis
        // and need to be set right now for each analysis stations. Alle these analysis stations are instances of _trf!
        if(!AC_discontis.empty())
            for(auto &AC: AC_discontis)
                AC.first->set_discontinuity(AC.second.refepoch, AC.second.pos0, AC.second.vel0, _disconts[AC.first->get_name(ivg::staname::cdp)]);
        
        
        // adjustments concerning the !CRF!
        if(opt & EST2CRF || opt & APR2CRF || opt & CRF2APR)
        {
            for(vector<ivg::Source>::iterator src_iter = _crf.begin(); src_iter!= _crf.end(); ++src_iter)
            {
                // get indexes of ra, dec position for each source included in _crf
                vector<int> position =  _param_list.get_indexes({ivg::paramtype::ra, ivg::paramtype::dec}, src_iter->get_name(ivg::srcname::iers) );
                
                // due to IERS name change in crf, source might not be found anymore in paramlist. checking for ivs name correspondence
                // due to some optimization this should be obsolete
//                if(position.at(0) == -1 && position.at(1) == -1 && position.size() == 2)
//                {
//                    position =  _param_list.get_indexes({ivg::paramtype::ra, ivg::paramtype::dec}, src_iter->get_name(ivg::srcname::ivs) );
//                    if(position.at(0) != -1 && position.at(1) != -1 && position.size() == 2)
//                    {
//                        _param_list.get_param(position.at(0))->set_name(src_iter->get_name(ivg::srcname::iers));
//                        _param_list.get_param(position.at(1))->set_name(src_iter->get_name(ivg::srcname::iers)); 
//                   }
//                }
                // for calculating correct btpb we need to know the positions of the sources (and later on also of the EOPs)
                fix_eops_sources.insert(fix_eops_sources.end(), position.begin(), position.end());

                // only if RA,DEC (2 parameters) from _crf have been found in _param_list -> proceed
                if(position.at(0) != -1 && position.at(1) != -1 && position.size() == 2)
                {
                    // only if epochs of RA, DEC are equal -> proceed
                    if( _param_list.get_param(position.at(0))->get_epoch() == _param_list.get_param(position.at(1))->get_epoch()  )
                    {     
                        // adjust internal _crf to values from parameter-vector from snx file
                        if(opt & EST2CRF || opt & APR2CRF)
                        {
                            double ra_rad, dec_rad;
                            if(opt & EST2CRF)
                            {
                                ra_rad = _param_list.get_param(position.at(0))->get_estimate();
                                dec_rad = _param_list.get_param(position.at(1))->get_estimate();
                            }
                            else if(opt & APR2CRF)
                            {
                                ra_rad = _param_list.get_param(position.at(0))->get_apriori();
                                dec_rad = _param_list.get_param(position.at(1))->get_apriori();
                            }

                            src_iter->set_ra0(ra_rad);
                            src_iter->set_dec0(dec_rad);

                            _param_list.get_param(position.at(0))->set_apriori(src_iter->get_ra0());
                            _param_list.get_param(position.at(1))->set_apriori(src_iter->get_dec0());

                            log<DETAIL>("*** ") % src_iter->get_name(ivg::srcname::iers) % " _param_list [ra,dec] --> _crf";
                        }
                        // adjust parameter-vector values from snx file to values from internal _crf
                        else if(opt & CRF2APR)
                        {
                            double old_ra0 = _param_list.get_param(position.at(0))->get_apriori();
                            double old_dec0 = _param_list.get_param(position.at(1))->get_apriori();
                            
                            double new_ra0 = src_iter->get_ra0();
                            double new_dec0 = src_iter->get_dec0();
                            // it might happen that no apriori positions exist in the selected apriori source catalog
                            if(new_ra0 == 0.0 || new_dec0 == 0.0)
                            {     
                                int h,m,deg,min;
                                double s,sec;
                                ivg::Source(ivg::srcname::ivs,src_iter->get_name(ivg::srcname::ivs),old_ra0,old_dec0).get_position(h,m,s,deg,min,sec);
                                
                                stringstream ss;
                                ss <<  "vector<double> Session::_adjust_data_storage(int opt): No ra/dec-position for ";
                                ss << src_iter->get_name(ivg::srcname::ivs) << "|";
                                ss << src_iter->get_name(ivg::srcname::iers);
                                ss << ". No apriori transformation possible. Could use: ";
                                ss << src_iter->get_name(ivg::srcname::icrf) << " " << src_iter->get_name(ivg::srcname::ivs) << " ";
                                ss << h << " " << m << " " << setprecision(8) << fixed << s << " ";
                                ss << deg << " " << min << " " << setprecision(7) << fixed << sec; 
                                throw runtime_error(ss.str());
                            }
                            
                            // check differences between OLD and NEW aprioris in right ascension and declination
                            // it's a bad sign if there are too big! name problems?!
                            if(abs(new_ra0 - old_ra0) > 1.5e-8 || abs(new_dec0 - old_dec0) > 1.5e-8)
                            {
                                
                                cerr << scientific << setprecision(10) << old_ra0 << " vs " << new_ra0 << endl;
                                cerr << scientific <<  setprecision(10) << old_dec0 << " vs " << new_dec0 << endl;
                                
                                stringstream ss;
                                ss <<  "vector<double> Session::_adjust_data_storage(int opt): apriori differences too big [ivs: ";
                                ss << src_iter->get_name(ivg::srcname::ivs) << "],[iers: ";
                                ss << src_iter->get_name(ivg::srcname::iers) << "], [icrf: ";
                                ss << src_iter->get_name(ivg::srcname::icrf) << "] ";
                                ss << "diff_ra = " << abs(new_ra0 - old_ra0)*ivg::rad2mas << " / diff_dec = " <<  abs(new_dec0 - old_dec0)*ivg::rad2mas;
                                ss << " [mas]";
                                throw runtime_error(ss.str());
                            }
                            
                            _param_list.get_param(position.at(0))->set_apriori(new_ra0);
                            _param_list.get_param(position.at(1))->set_apriori(new_dec0);

                            log<DETAIL>("*** ") % src_iter->get_name(ivg::srcname::iers) % " _crf --> _param_list [ra,dec]";
                        }

                        // scale RA and DEC of sources      
                        scales.at(position.at(0)) =  rad2mas;
                        scales.at(position.at(1)) =  rad2mas;

                        adjustments += 2; // because ra and dec
                    }
                    else
                        throw runtime_error( "vector<double> Session::_adjust_data_storage(int opt)): Epochs of RA.DEC of source "+src_iter->get_name(ivg::srcname::iers)+" not equal.");
                }
                //if the same source is parameterized in _param_list at more than one epoch
                else if(position.size() > 2)
                    throw runtime_error( "vector<double> Session::_adjust_data_storage(int opt): Source parameterized in snx file at more than one epoch.");
                else
                {                  
                    // find out which type of source for output
                    string src_type = "         default";
                    if(src_iter->is_special_handling())
                        src_type = "special handling";
                    else if(src_iter->is_defining())
                        src_type = "        defining";
                    
                    // if a source is already in the _crf because it is stored in the sinex SOURCE/ID block (used for initialization),
                    // it might happen that this source is NOT in the NEQ of the sinex file. For example BKG did this. These sources were fixed to there aprioris!
                    // In general, these sources, which are not in _param_list but in _crf, need to be deleted from _crf
                    int num_src = _crf.get_number_sources()-1;
                    log<WARNING>("!!! ") % src_type % " " % src_iter->get_name(ivg::srcname::iers) % " in _crf[#" % num_src % "] not parametrized. Seems to be already fixed in NEQ.";
                    _crf.remove_source(src_iter);
                    src_iter--;
                    // it's critical and should be reconsidered what to do with this kind of session!
//                    throw runtime_error("vector<double> Session::_adjust_data_storage(int opt): Session not useable due to prefixed sources.");
                }

            }
        }
            
        // now after CRF and TRF is considered, we focus on the EOPs
        if(opt & EST2EOP || opt & APR2EOP || opt & EOP2APR)
        {
            // get indexes of xpo, ypo, ut1, nutx, nuty. This does not include the rates!!!
            vector<int> erp_idx = _param_list.get_indexes({xpo, ypo, ut1, nutx, nuty},"EOP"); // ATTENTION: returns 5 even if the ERPs are set up as cpwlf!
            // only if ALL 5 EOPs exist, the session is suitable for combination
            if(erp_idx.at(0) != -1 && erp_idx.at(1) != -1 && erp_idx.at(2) != -1 && erp_idx.at(3) != -1 &&  erp_idx.at(4) != -1 && erp_idx.size() == 5 )
            {
                // check if epochs of xpo,ypo,ut1 are equal (in general the case)
                if(_param_list.get_param(erp_idx.at(0))->get_epoch() == _param_list.get_param(erp_idx.at(1))->get_epoch() 
                        && _param_list.get_param(erp_idx.at(1))->get_epoch() == _param_list.get_param(erp_idx.at(2))->get_epoch() )
                {
                    ivg::Date epoch = _param_list.get_param(erp_idx.at(0))->get_epoch();
                    ivg::Matrix erp_now, erp_rate;
                    if( opt & EOP2APR )
                    {
                        ivg::Matrix erp_before = _eops.calc_erp(epoch.add_days(-2.0), "linear", true);
                        erp_now = _eops.calc_erp(epoch, "linear", true);
                        ivg::Matrix erp_after = _eops.calc_erp(epoch.add_days(2.0), "linear", true);

                        // all rates in radiant
                        erp_rate = ( erp_after - erp_before ) / 4.0;
                    }

                    // now begin to loop over params
                    for( std::vector<ivg::Param>::iterator param = _param_list.begin(); param != _param_list.end(); ++param )
                    {
                        if( opt & EOP2APR )
                        {
                            if(param->is_type({xpo, ypo, ut1},{0,1}) || param->is_type({nutx,nuty},{0}))
                            {
                                // XPO
                                if(param->get_type() == ivg::paramtype::xpo && param->get_order() == 0)
                                    param->set_apriori(erp_now(0));
                                // XPOR
                                else if(param->get_type() == ivg::paramtype::xpo && param->get_order() == 1)
                                    param->set_apriori(erp_rate(0));
                                // YPO
                                else if(param->get_type() == ivg::paramtype::ypo && param->get_order() == 0)
                                    param->set_apriori(erp_now(1));
                                // YPOR
                                else if(param->get_type() == ivg::paramtype::ypo && param->get_order() == 1)
                                    param->set_apriori(erp_rate(1));
                                // UT1
                                else if(param->get_type() == ivg::paramtype::ut1 && param->get_order() == 0)
                                {
                                    // correction for leap seconds
                                    double apri = param->get_apriori();
                                    if( apri > -1500.0/ivg::param_unit_fac.at(ut1) && apri < 1500.0/ivg::param_unit_fac.at(ut1) )
                                    {
                                        double new_apri = erp_now(2)*ivg::param_unit_fac.at(ut1) + (param->get_epoch().get_leap_sec()*1000);
                                        param->set_apriori(new_apri / ivg::param_unit_fac.at(ut1));
                                    }
                                    else
                                        param->set_apriori(erp_now(2));
                                    
                                }
                                // LOD
                                else if(param->get_type() == ivg::paramtype::ut1 && param->get_order() == 1)
                                    param->set_apriori(erp_rate(2));

                                log<DETAIL>("*** EOP _eops --> _param_list[") % param->get_typename() % " " % param->get_order() % "]";

                                adjustments += 1;
                            }

                        }
                        else if(opt & EST2EOP || opt & APR2EOP)
                        {
                            if(param->is_type({xpo, ypo, ut1},{0, 1}) || param->is_type({nutx,nuty},{0}))
                            {
                                log<WARNING>("*** EOP _param_list[") % param->get_typename() % " " % param->get_order() % "] --> _eops NOT IMPLEMENTED YET!!!!";
                                adjustments += 1;
                            }
                        }  
                    }
                }
                else
                    throw runtime_error("vector<double> Session::_adjust_data_storage(int opt): Epochs of XPO, YPO, UT1 not equal - Unexpected");
            }
            // if only NUT_X and NUT_Y is parameterized
            else if(erp_idx.at(0) == -1 && erp_idx.at(1) == -1 && erp_idx.at(2) == -1 && erp_idx.at(3) != -1 && erp_idx.at(4) != -1 &&  erp_idx.size() == 5)
                throw runtime_error("ERROR: Only [nutx, nuty] existent");
            else if(erp_idx.at(0) == -1 && erp_idx.at(1) == -1 && erp_idx.at(2) == -1 && erp_idx.at(3) == -1 && erp_idx.at(4) == -1 && erp_idx.size() == 5)
                throw runtime_error("ERROR: No EOPs existent");
            else
                throw runtime_error("vector<double> Session::_adjust_data_storage(int opt): Unexpected EOP parameterization in _param_list");
        }
           
        // Setting nutation to zero right here, leads to a consideration within the apriori transformation
        // -> we don't want this
        // setting nutation to zero (X,Y as well as LN,OB)
//        vector<int> nut_idx = _param_list.get_indexes({nutx, nuty, nutln, nutob},"EOP");
//        for(int i=0; i<nut_idx.size(); i++)
//        {
//            if( nut_idx.at(i) != -1 )
//            {
//                if( opt & NUT2ZER )
//                {
//                    _param_list.get_param(nut_idx.at(i))->set_apriori(0.0);
//                    adjustments += 1;
//                    log<DETAIL>("*** EOPs _param_list [nut (x,y,ln,ob)] setting to 0.0");
//                }
//                else
//                    adjustments += 1;
//            }
//        } 

        // check for other parameters left in NEQ, like zwd, ngr, egr and REDUCE them
        for(auto &param: _param_list)
        {
            if(param.is_type({zwd, ngr, egr},{0,1}))
            {
                // problems reducing singular system
                param.set_reduce_flag(true);
                log<DETAIL>("*** ") % param.get_typename() % " of " % param.get_name() % " at " % param.get_epoch().get_double_mjd() % " --> REDUCE";
                adjustments += 1;
            }
        }
        
        if(adjustments == _param_list.size())
            log<INFO>("*** All parameter in _param_list adjusted [") % adjustments % "/" % _param_list.size() % "]";
        else
            throw runtime_error( "ERROR: Not all parameter adjusted in _param_list ["+std::to_string(adjustments)+"/"+std::to_string(_param_list.size())+"]" );
    
        // returning scales vector. Up to now only for sources implemented.
        return scales;
#ifdef DEBUG_VLBI
   cerr << "--- Session::_adjust_data_storage(uint32_t opt) "
        << ": " << tim.toc() << " s " << endl;
#endif
}
// ...........................................................................
void Session::write_backend_files( std::string dir )
// ...........................................................................
{
    // copy param list
    ivg::Param_list param_list = _param_list;
    std::vector<int> param_idx;

    // get original index vector of param list
    for( int i=0; i < param_list.size(); i++ )
        param_idx.push_back(i);

    int idx;
    int counter = 0;

    // loop over all params
    for( std::vector<ivg::Param>::iterator it = _param_list.begin();
            it != _param_list.end(); ++it )
    {
        idx = it - _param_list.begin();

/* FUER CORINNA  
          // search clock parameters
        if( it->get_type() == 3 && it->get_order() != 0 )
        {
            // re-sort clock parameters and corresponding index vector
            std::vector<int> indexes = param_list.get_idx( ivg::paramtype::clo,
                                       it->get_name() );

            int rem_idx = idx-counter;
            if( indexes.size() == 1 )
            {
                if( indexes.at(0)+it->get_order() < idx )
                   rem_idx = rem_idx+1;
                else
                    counter++;
                param_list.insert_param( indexes.at(0)+it->get_order(), *it );
                param_idx.insert( param_idx.begin()+indexes.at(0)+it->get_order(), idx );
            }
            else
            {
                if( indexes.at(1) < idx )
                   rem_idx = rem_idx+1;
                else
                    counter++;
                param_list.insert_param( indexes.at(1), *it );
                param_idx.insert( param_idx.begin()+indexes.at(1), idx );
            }

            // remove old clock parameters and refresh index vector
            param_list.remove_param( rem_idx );
            param_idx.erase( param_idx.begin()+rem_idx );
        }
*/
    }

    // write modified design matrix, o-c vector and weight matrix/vector
    _lsa_solution.write_backend_gmm( dir, _name, param_idx );


    // ----- observation dependent: write 'obs_info' file
    std::string obs_file = dir+_name+"obs_info"+"XY";
    ofstream out( obs_file.c_str() );

    std::string zhd_file = dir+_name+"apriorisAt"+"XY";
    ofstream zhd_out( zhd_file.c_str() );

    double jd;
    std::string sta1, sta2, src;
    ivg::Matrix k;
    int sta_idx1, sta_idx2;
    ivg::Analysis_station * sta_ptr1;
    ivg::Analysis_station * sta_ptr2;
    ivg::Matrix azel1, azel2, trf2crf;
    counter = 1;

    //stringstream out;
    for( int i=0; i <= _scans.size()-1; i++ )
    {
        for( int j=0; j <= _scans.at(i).get_nobs()-1; j++ )
        {

//            jd = _scans.at(i).get_obs_ptr(j)->get_epoch().get_jd();
            jd = _scans.at(i).get_obs_ptr(j)->get_epoch().get_double_mjd();

            _scans.at(i).get_obs_ptr(j)->get_station_names( sta1, sta2 );
            _scans.at(i).get_obs_ptr(j)->get_source_name( src );

            k = _scans.at(i).get_source()->get_unit_vector_ssb();

            sta_idx1 = _scans.at(i).get_obs_ptr(j)->get_scan_idx(1);
            sta_idx2 = _scans.at(i).get_obs_ptr(j)->get_scan_idx(2);

            sta_ptr1 = _scans.at(i).get_data( sta_idx1 ).sta_ptr;
            sta_ptr2 = _scans.at(i).get_data( sta_idx2 ).sta_ptr;

            trf2crf = _scans.at(i).get_trf2crf_matrix();

            azel1 = sta_ptr1->calc_az_el( _scans.at(i).get_obs_ptr(j)->get_epoch(), k, trf2crf.transpose() );
            azel2 = sta_ptr2->calc_az_el( _scans.at(i).get_obs_ptr(j)->get_epoch(), k, trf2crf.transpose() );

            // write output obs_file
            out << setiosflags(ios::right) << setiosflags(ios::fixed)
                << setw(5) << counter << "  "
                << setprecision(16) << setw(24) << jd << " ";
            out << left << setw(8) << src << " "
                << setw(8) << sta1 << " "
                << setw(15) << sta2 << " "
                << setprecision(16) << setw(24) << azel1(0) << "  "
                << setprecision(16) << setw(24) << azel2(0) << "  "
                << setprecision(16) << setw(24) << azel1(1) << "  "
                << setprecision(16) << setw(24) << azel2(1) << endl;


            // write output file 'apriorisAT' (ZHD apriori values)
            zhd_out << left << setw(8) << sta1 << " "
                    << fixed << setprecision(8) << setw(17) << jd << " "
                    << right << scientific << setprecision(6) << setw(13) << 0.0 << endl;
            zhd_out	<< left << setw(8) << sta2 << " "
                    << fixed << setprecision(8) << setw(17) << jd << " "
                    << right << scientific << setprecision(6) << setw(13) << 0.0 << endl;

            counter++;
        }
    }
    out.close();
    zhd_out.close();

    // ----- parameter dependent: write files for parameter, apriori values and reference epochs
    std::string param_file = dir+_name+"param"+"XY";
    ofstream param_out( param_file.c_str() );

    std::string apriori_file = dir+_name+"aprioris"+"XY";
    ofstream apriori_out( apriori_file.c_str() );

    std::string epoch_file = dir+_name+"ref_epochs"+"XY";
    ofstream epoch_out( epoch_file.c_str() );

    double apr;
    ivg::paramtype type;
    std::string ptype, name;
    int order;

    counter = 1;
    int cl_counter = 0;

    for( int i=0; i < param_list.size(); i++ )
    {
        jd = param_list.get_param(i)->get_epoch().get_jd();
        name = param_list.get_param(i)->get_name();
        type = param_list.get_param(i)->get_type();
        order = param_list.get_param(i)->get_order();
        std::string order_str = std::to_string(order);

        if( type == 0 )
            ptype = "  X";
        else if( type == 1 )
            ptype = "  Y";
        else if( type == 2 )
            ptype = "  Z";
        else if( type == 3 )
        {
            if( order != 0 )
                ptype = " CL";
            else
            {
                if( cl_counter == 0 )
                    ptype = " CL";
                else
                    ptype = " cl";
                cl_counter++;
            }
        }
        else if( type == 4 )
        {
            if( order != 0 )
                ptype = " AT";
            else
                ptype = " at";
        }
        else if( type == 5 )
        {
            ptype = " NG";
        }
        else if( type == 6 )
        {
            ptype = " EG";
        }
        else if( type == 7 )
        {
            ptype = " "+order_str;
            name = "X Wobble";
        }
        else if( type == 8 )
        {
            ptype = " "+order_str;
            name = "Y Wobble";
        }
        else if( type == 9 )
        {
            ptype = " "+order_str;
            name = "UT1-TAI";
        }
        else if( type == 10 )
        {
            ptype = "";
            name = "Nut longitude";
        }
        else if( type == 11 )
        {
            ptype = "";
            name = "Nut obliquity";
        }
        else if( type == 14 )
        {
            ptype = " BR";
        }

        // write output file 'param'
        param_out << setiosflags(ios::right) << setiosflags(ios::fixed)
                  << setw(5) << counter << left
                  << " " << setw(8) << name << ptype
                  << endl;

        // write output file 'aprioris'
        apr = param_list.get_param(i)->get_apriori();

        apriori_out << setiosflags(ios::right) << setiosflags(ios::fixed)
                    << setw(7) << counter << "  "
                    << scientific << setprecision(40) << setw(17) << apr << " " << endl;

        // write output file 'ref_epochs'
        epoch_out << setiosflags(ios::right) << setiosflags(ios::fixed)
                  << setw(7) << counter << "  "
                  << setprecision(8) << setw(17) << jd << " " << endl;

        counter++;
    }

    param_out.close();
    apriori_out.close();
    epoch_out.close();
}

// ...........................................................................
bool Session::operator==( const Session test ) const
// ...........................................................................
{
    bool out = false;
    if( _name == test._name )
        out = true;

    return out;
}

// ...........................................................................
void Session::show()
// ...........................................................................
{
    cerr << "------------------- Session.show() ------------------------" << endl;
    cerr << "DB: " << _name << endl;
    cerr << "#scans: " << _scans.size() << endl;
    cerr << "#obs:   " << _nobs  << endl;
    vector<string> station_names = _trf.get_station_names(ivg::staname::ivs_name);
    cerr << "TRF (#stations " << station_names.size() << ")" << endl;
//    _trf.show();
    show_vector(station_names);
    vector<string> source_names = _crf.get_source_names(ivg::srcname::ivs);
    cerr << "CRF (#sources " << source_names.size() << ")" << endl;
//    _crf.show();
    show_vector(source_names);
    cerr << "_param_list (#params " << _param_list.size() << ")" << endl;
    _param_list.show();
    cerr << "_neq_solution: " << endl;
    _neq_solution.show();
    
    cerr << "------------------- Session.show() ------------------------" << endl;
}


// ...........................................................................
string Session::coverage2string( std::map<int, std::vector<double> > coverage, unsigned precision, bool percent ){
// ...........................................................................
    stringstream ss;
    int max_tree_level = (int)(*_setup)["SKED"]["max_tree_level"];
    ss << "level:  ";
    
    int p = 5 + precision;
    if(precision == 0 && percent == true){
        p--;
    }
        
    for(unsigned level = 1; level <= max_tree_level ; level++){
        ss << setw(p+1) << level;
    }
    ss << std::endl;
    
    
    for(size_t i=0; i < _trf.get_number_stations(); ++i){
        ss << setw(9) << left << _trf.get_station(i)->get_name(ivg::staname::ivs_name);
        for(unsigned level = 1; level <= max_tree_level ; level++){
            if(percent){
                ss << right << setw(p) << setprecision(precision) << fixed << coverage[i][level-1]*100 << "%";
            } else {
                ss << right << setw(p) << setprecision(0) << fixed << coverage[i][level-1];
            }
        }
        ss << std::endl;
    }
    
    return ss.str();
}

// ...........................................................................
void Session::show_coverage( std::map<int, std::vector<double> > coverage ){
// ...........................................................................
    std::cout << coverage2string(coverage);
}

// ...........................................................................
ivg::Date Session::get_mid_epoch()
// ...........................................................................
{
    double mjd = (_start.get_double_mjd() + _end.get_double_mjd()) / 2.0;
    return ivg::Date(mjd);
}

// ...........................................................................
ivg::Matrix Session::get_obs_epochs( std::string station )
// ...........................................................................
{
    std::string sta1, sta2;  
    std::vector<int> obs_idx;
    
    ivg::Matrix obs_epochs( _nobs,1,0.0 );
    int c = 0;
    
    for( int i=0; i <= _scans.size()-1; i++ )
    {
        for( int j=0; j <= _scans.at(i).get_nobs()-1; j++ )
        {            
            obs_epochs(c,0) = _scans.at(i).get_obs_ptr(j)->get_epoch().get_double_mjd();
            c++;

            if( station != "ALL")
            {
                _scans.at(i).get_obs_ptr(j)->get_station_names( sta1, sta2 );            

                if( station == sta1 || station == sta2 )
                    obs_idx.push_back(i);
            }
        }
    }
    if( station != "ALL")
        obs_epochs = obs_epochs(obs_idx);
    
    return obs_epochs;
}

// ...........................................................................
void Session::_determine_nobs_ngs( std::string filename )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Session::_determine_nobs_ngs( std::string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream fin( filename.c_str() );
    if( !fin.is_open() )
        throw runtime_error( "int Session::_determine_nobs_ngs( std::string filename ): Failed to open file "+filename+"\n");

    // go to 5 spots before EOF, i.e., before identifier of final card entry
    // ( e.g., '09' for '325209' where the final result should be nobs = 3252)
    fin.seekg( -10,ios_base::end );
    string l;
    getline( fin, l, '\n' );
    int off;
    if( isdigit( l.back() ) )
       off = 2;
    else
       off = 3;
    fin.seekg( -off,ios_base::end );

    bool keepLooping = true;
    while( keepLooping )
    {
        char ch;
        // read card number
        fin.get( ch );
        if( ch != '1' )
            fin.seekg( -(80+off),ios_base::cur);
        else
        {
            fin.seekg(-80,ios_base::cur);
            string line;
            getline( fin, line, '\n' );

            int y = std::stoi( line.substr( 29,4 ) );
            int m = std::stoi( line.substr( 34,2 ) );
            int d = std::stoi( line.substr( 37,2 ) );
            int h = std::stoi( line.substr( 40,2 ) );
            int min = std::stoi( line.substr( 43,2 ) );
            double s = s2d( line.substr( 47,13 ) );

            ivg::Date epoch( y, m, d, h, min, s );
            _end = epoch;

            string num = line.substr( 72,6 );
            std::transform( num.begin(), num.end(), num.begin(), [](char ch) {
               return ch == ' ' ? '0' : ch;
            });
            _nobs = std::stoi( num );

            keepLooping = false;
        }
    }
    fin.close();
    
#if DEBUG_VLBI >=2
    cerr << "--- void Session::_determine_nobs_ngs( std::string )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Session::create_solution_info()
// ...........................................................................
{
    
       std::map< std::string,std::vector<int> > idx_sta, idx_src, idx_bl;
       ivg::Matrix data;
       //
       // matrix for obs data 
       // (mjd, residual, sigma_tau, sigma_oc, weight_factors, sigma_resid, partial redundancies, az1, el1, az2, el2, dist_rays)
       data.resize( _nobs,16,0.0 );
       vector<string> first_station,source_names;
       
       // get residuals
       ivg::Matrix resid;
       ivg::Matrix Qvv;
       log<INFO>("*** Creating Solution-Infos...");
       _lsa_solution.get_resid( resid,Qvv );

       // get partial redundancies 
       ivg::Matrix redundancies = _lsa_solution.get_partial_redundancies();
       // get weight matrix of LS adjustment
       ivg::Matrix W;
       ivg::Matrix wgt_facs;
       _lsa_solution.get_wgt_matrix( W,wgt_facs );
       if( W.cols() != 1 )
          W = W.diag();

       
       int maxiter = 1; // Just one pass through the for loop. Only group delay is saved in data
       
       if(_ambigRes == true) //check whether singleband is used
       {
           // one pass for group delay and another for singleband delay
           maxiter = delaytype::MAXDELAYTYPE;
           // _nobs is multiplied with 2 in void Session::init_vgosdb_ngs_solution()
           // in _ambigRes mode. In _ambigRes mode there is one data Matrix for group delay and one
           // for single band delay. So the size for data is _nobs/2.
           data.resize( _nobs/2,16,0.0 );    
       }
       
       // vector containing the position of each observation in the ncfile;
       vector<int> idx_obs_nc = this->get_obs_idxs_in_ncfile();
       
       // loop over delaytype. One pass if _ambigRes is false, two if not
       for(int cur_delay_type = delaytype::group; cur_delay_type < maxiter; cur_delay_type++)
       {   
            // loop over scans and observations therein
            int obs_counter = 0;
            double mjd;
            ivg::Source* src_ptr;
            ivg::Obs* obs_ptr;
            std::string sta1, sta2, src;
            ivg::Matrix azel1;
            ivg::Matrix azel2;
            for( vector<ivg::Scan>::iterator scan_it = _scans.begin(); scan_it != _scans.end(); ++scan_it )
            { 
         //      mjd = scan_it->get_epoch().get_double_mjd();
               mjd = scan_it->get_epoch().get_qcustomplot_date();
               for( int obs_i=0; obs_i < scan_it->get_nobs(); ++obs_i )
               {   
                  obs_ptr = scan_it->get_obs_ptr(obs_i);
                  
                  //the idx_ vectors should not contain the same entries twice!
                  //just push back once
                  if(cur_delay_type == 0)
                  {
                      obs_ptr->get_station_names( sta1, sta2 );
                      obs_ptr->get_source_name( src );
                      idx_sta[ sta1 ].push_back( obs_counter );
                      idx_sta[ sta2 ].push_back( obs_counter );
                      idx_src[ src  ].push_back( obs_counter );
                      idx_bl[ sta1+"-"+sta2 ].push_back( obs_counter );
                  }
                  
                  first_station.push_back(sta1);
                  source_names.push_back(src);
                  
                  azel1 = obs_ptr->get_az_el( 1 );
                  azel2 = obs_ptr->get_az_el( 2 );
                  
                  int idx = obs_counter;
                  if(_ambigRes) //ckeck whether singleband is used
                  {
                      // The entries in the residual vector are alternating group- and single band delays
                      // group delay: all even numbers -> obs_counter=0,2,4,... 2*nobs
                      // single band delay: all odd numbers -> obs_counter=1,3,5,...2*nobs
                      idx = obs_counter*2 + cur_delay_type;
                  }
                  
                  data( obs_counter,0 ) = mjd;
                  data( obs_counter,1 ) = resid( idx );
//                  data( obs_counter,1 ) = resid( idx )* sin( azel1( 1 ) );                  
                  data( obs_counter,2 ) = sqrt(obs_ptr->get_obs_variance( false,_phaseDelay,false ));
                  data( obs_counter,3 ) = sqrt( 1.0/W( idx,0 ) );
                  data( obs_counter,4 ) = wgt_facs( idx );
                  data( obs_counter,5 ) = sqrt( Qvv( idx,idx ) );
                  data( obs_counter,6 ) = redundancies( idx );
                  data( obs_counter,7 ) = azel1( 0 );
                  data( obs_counter,8 ) = azel1( 1 );
                  data( obs_counter,9 ) = azel2( 0 );
                  data( obs_counter,10 ) = azel2( 1 );
                  data( obs_counter,11 ) = calc_dist_between_rays( obs_ptr );
                  data( obs_counter,12 ) = obs_ptr->get_snr_bx();
                  data( obs_counter,13 ) = obs_ptr->get_snr_bs();
                  data( obs_counter,14 ) = scan_it-_scans.begin();
		  data( obs_counter,15 ) = obs_i;
		    
                  obs_counter++;
               }
            }

         // we got 4 different categories of residuals -> ALL, STATION, SOURCE, BASELINE 
         // preparing _residuals for easy plotting in interactive ascot
         map< residtype, map< string,vector<int> > > resid_categories;
         ivg::Matrix data_rows(0.0,1.0,data.rows()-1,1);       
         resid_categories[residtype::all] = {{"ALL", vector<int>(data_rows.begin(),data_rows.end()) }};
         resid_categories[residtype::station] = idx_sta;         
         resid_categories[residtype::source] = idx_src;
         resid_categories[residtype::baseline] = idx_bl;

         //generate cols vector based on number of cols from data matrix
         ivg::Matrix data_cols(0.0,1.0,data.cols()-1,1);
         vector<int> cols = vector<int>(data_cols.begin(),data_cols.end());
         
         for(auto &resid_group: resid_categories)
         {
             for(auto &it: resid_group.second ) 
             {            
                 Residual tmp;
                 tmp.name = it.first; // e.g. GILCREEK
                 tmp.type = resid_group.first; // e.g. residtype::station
                 tmp.data = data.get_sub( it.second, cols );
                 tmp.idx_in  = tmp.data( ":",4 ).find_idx( eq, 1 );
                 tmp.idx_out = tmp.data( ":",4 ).find_idx( ne, 1 );
                 tmp.idx_origin = it.second;
                 
                 vector<string> sta1_tmp, src_tmp;
                 for(auto &index: it.second)
                 {
                     sta1_tmp.push_back(first_station.at(index));
                     src_tmp.push_back(source_names.at(index));
                 }
                 
                 // save name of first sessions and source name in order to get relation of azi and ele (skyplot)
                 tmp.first_station = sta1_tmp;
                 tmp.source_names = src_tmp;
                 
                 // only calculate wrms if atleast two observations are not outliers
                 if(tmp.idx_in.size() > 1)
                 {
                     // get residuals where obs are not outliers
                     ivg::Matrix r_in = tmp.data.get_sub( tmp.idx_in, {1} )*1e12;
                     ivg::Matrix p_in = (tmp.data.get_sub( tmp.idx_in, {3} )*1e12)^(-2);

                     // calculate WRMS
                     tmp.wrms = sqrt( (r_in.transpose()*p_in.diag()*r_in)(0)/(p_in.sum_col())(0) );
                 }
                 else
                     tmp.wrms = 0.0;

                 vector<int> idx_ncfile(tmp.idx_origin.size());
                 for(unsigned int i = 0; i < tmp.idx_origin.size(); ++i){
                     idx_ncfile[i] = idx_obs_nc[tmp.idx_origin[i]];
                 }
                 
                 tmp.idx_ncfile = idx_ncfile;

                 _residuals[delaytype(cur_delay_type)].push_back(tmp);
             }
         }
        }
    log<INFO>("*** Solution-Infos created.");
}
// ...........................................................................
double Session::calc_dist_between_rays( ivg::Obs* obs_ptr )
// ...........................................................................
{

      ivg::Matrix az_el;
      ivg::Date epoch; 
      int sta_idx1, sta_idx2;
      ivg::Matrix sta1, sta2;
      ivg::Analysis_station* sta_ptr1, sta_ptr2;

      // get observation epoch
      epoch = obs_ptr->get_epoch();

      // get azimuth and elevation   
      az_el = obs_ptr->get_az_el(1);
      
      // get station coordinates and calculate baseline lengths
      sta_idx1 = obs_ptr->get_scan_idx(1);
      sta_idx2 = obs_ptr->get_scan_idx(2);
      
      sta1 = obs_ptr->get_scan()->get_sta_ptr( sta_idx1 )->calc_xyz(epoch);
      sta2 = obs_ptr->get_scan()->get_sta_ptr( sta_idx2 )->calc_xyz(epoch);

      // sta1 = _scans.at(i).get_sta_ptr( sta_idx1 )->calc_xyz(epoch);
      // sta2 = _scans.at(i).get_sta_ptr( sta_idx2 )->calc_xyz(epoch);
      ivg::Matrix bl = sta1 - sta2;
      
      // (1) calculate local source vector (homogeneous coordinates) using elevation and azimuth
      double kx = cos(az_el(1)) * cos(az_el(0));
      double ky = cos(az_el(1)) * sin(az_el(0));
      double kz = sin(az_el(1));
      
      ivg::Matrix k_homog( 1,4, 1.0 );
      k_homog( 0,0 ) = kx;
      k_homog( 0,1 ) = ky;
      k_homog( 0,2 ) = kz;
      
      // (2) transform to global system
      ivg::Matrix T1( 1,3,0.0 ); T1( 0,2 ) = 1.0 ;
      ivg::Matrix T2( 1,3,0.0 ); T2( 0,0 ) = 1.0 ;
      ivg::Matrix tmp = sta1;
      tmp(2,0) = 0.0;
   
      // calculate rotation angles
      double beta = acos( ( T1 * sta1 )(0) / (sta1.norm())(0) ) ;
      double delta = M_PI - acos( ( T2 * tmp )(0) / (tmp.norm())(0) );
      
      // set rotation matrices: R1 = rotation about y-axis; R2 = rotation about z-axis
      ivg::Matrix R1( 4,4,0.0 );
      R1( 0,0 ) = cos( beta );  R1( 0,1 ) = 0.0;  R1( 0,2 ) = sin( beta );  R1( 0,3 ) = 0.0 ;
      R1( 1,0 ) = 0.0; 	     R1( 1,1 ) = 1.0;  R1( 1,2 ) = 0.0; 	 R1( 1,3 ) = 0.0 ;
      R1( 2,0 ) = -sin( beta ); R1( 2,1 ) = 0.0;  R1( 2,2 ) = cos( beta );  R1( 2,3 ) = 0.0 ;
      R1( 3,0 ) = 0.0; 	     R1( 3,1 ) = 0.0;  R1( 3,2 ) = 0.0; 	 R1( 3,3 ) = 1.0 ;
      
      ivg::Matrix R2( 4,4,0.0 );
      R2( 0,0 ) = cos( delta );  R2( 0,1 ) = -sin( delta ); R2( 0,2 ) = 0.0;  R2( 0,3 ) = 0.0 ;
      R2( 1,0 ) = sin( delta );  R2( 1,1 ) = cos( delta );  R2( 1,2 ) = 0.0;  R2( 1,3 ) = 0.0 ;
      R2( 2,0 ) = 0.0;  	      R2( 2,1 ) = 0.0;  	 R2( 2,2 ) = 1.0;  R2( 2,3 ) = 0.0 ;
      R2( 3,0 ) = 0.0; 	      R2( 3,1 ) = 0.0;  	 R2( 3,2 ) = 0.0;  R2( 3,3 ) = 1.0 ;
      
      // transform local source vector to global (cartesian) system
      ivg::Matrix M = R2.transpose() * R1.transpose() ;
      ivg::Matrix ks = M * k_homog.transpose() ;
      
      // "re-transformation" from homogeneous coordinates
      ivg::Matrix k = ks.get_sub( 0, 0, 2, ks.size(2)-1 );
      
      // (3) calculate separation distance between two rays
      ivg::Matrix tau = bl.transpose() * k;
      ivg::Matrix ones( tau.size(2), 1, 1.0 );
      	
      double dist = ( ( (ones* (bl.transpose()* bl)) - tau.transpose().mult_elem( tau.transpose() ) ).sqrt() )(0);
    
      return dist; 
}


// ...........................................................................
void Session::_eliminate_data()
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Session::_eliminate_data()" << endl; 
   tictoc tim;
   tim.tic();
#endif

    ivg::Matrix outlier_indicator( _nobs,1,1.0 );
    string out_file = "/data/bakkari_outliers/"+_name+"_out.txt";
    if( file_exists( out_file ) &&  (bool)(*_setup)["outliers"]["load"][0] && _type == "ngs")
       outlier_indicator.load_ascii( out_file );

    
    // create list of stations and source in this solution
    vector<string>sta_names = _trf.get_station_names( ivg::staname::ivs_name );
    vector<string>src_names = _crf.get_source_names( ivg::srcname::iers );
    vector<string>src_names_ivs = _crf.get_source_names( ivg::srcname::ivs );
    
    vector<string>sta_elim,src_elim,bl_elim;

    for(auto &sta_elim_tmp: _trf)
    {
        ivg::Matrix pos0 = sta_elim_tmp.get_xyz0();

        if(pos0(0) == 0.0 && pos0(1) == 0.0 && pos0(2))
        {
            log<INFO>("*** Data elimination: Eliminate station ") % sta_elim_tmp.get_name(ivg::staname::ivs_name) % " because pos0 = [0.0, 0.0, 0.0]";
            sta_elim.push_back(sta_elim_tmp.get_name(ivg::staname::ivs_name));
        }
    }

    // check whether one of them are in the eliminate lists and create new lists 
    // optimized for this session
    // create also indices of corresponding parameters
    vector<int> rem_par;
    map<string, vector<ivg::paramtype> > assignment =
    {
       { "stations", {ivg::paramtype::stax, ivg::paramtype::stay, 
                      ivg::paramtype::staz, ivg::paramtype::clo,
                      ivg::paramtype::ngr, ivg::paramtype::egr,
                      ivg::paramtype::zwd}}, 
       { "sources",  {ivg::paramtype::ra,ivg::paramtype::dec}}
    };


    // find stations which should be eliminated
    // (1) from sta_elim block
    for( int i=0; i<(*_setup)["elim_sta"].getLength(); ++i )
    {
       string sta = (*_setup)["elim_sta"][i];
       vector<string>::iterator it = find( sta_names.begin(), sta_names.end(), sta );
       if( it != sta_names.end() )
       {
          sta_elim.push_back( *it );
          log<INFO>("*** Data elimination: Eliminate station ") % *it % " because of elim_sta in elim-block";
          vector<int> idx = _param_list.get_indexes( assignment[ "stations" ], sta );
          rem_par.insert( rem_par.end(), idx.begin(), idx.end() );
       }
    }
    // (2) from DB handling block
    if((bool)(*_handling).exists("handling") && (bool)((*_handling)["handling"]).exists("elim_sta"))
    {
       for( int i=0; i<(*_handling)["handling"]["elim_sta"].getLength(); ++i )
       {
          string sta = (*_handling)["handling"]["elim_sta"][i];
          // check whether station has been used in this session
          vector<string>::iterator it1 = find( sta_names.begin(), sta_names.end(), sta );
          // check whether station was already inlcuded in general sta_elim block
          vector<string>::iterator it2 = find( sta_elim.begin(), sta_elim.end(), sta );
          if( it1 != sta_names.end() && it2 == sta_elim.end() )
          {
             sta_elim.push_back( *it1 );
             log<INFO>("*** Data elimination: Eliminate station ") % *it1 % " because of elim_sta in handling";
             vector<int> idx = _param_list.get_indexes( assignment[ "stations" ], sta );
             rem_par.insert( rem_par.end(), idx.begin(), idx.end() );
          }
       }
    }
    for (int i=0;i<sta_elim.size();i++) {
      ivg::Analysis_station *tmp;
    
      vector<string>::iterator it1 = find( sta_names.begin(), sta_names.end(), sta_elim[i] );
      tmp=_trf.get_station(it1-sta_names.begin());
      tmp->set_num_obs(0);
      
      
    }

    for( int i=0; i<(*_setup)["elim_src"].getLength(); ++i )
    {
       string src = (*_setup)["elim_src"][i];
       for (int j=0;j<src_names_ivs.size();j++)
	 {
	   if (src == src_names_ivs.at(j))
	       src=src_names.at(j);
	 }
       vector<string>::iterator it = find( src_names.begin(), src_names.end(), src );
       vector<string>::iterator it2 = find( src_elim.begin(), src_elim.end(), src );
       if( it != src_names.end() && it2 == src_elim.end())
       {
          src_elim.push_back( *it );
          log<INFO>("*** Data elimination: Eliminate source ") % *it % " because of elim_src in elim-block";
          vector<int> idx = _param_list.get_indexes( assignment[ "sources" ], src );
          rem_par.insert( rem_par.end(), idx.begin(), idx.end() );
       }
    }
    // (2) from DB handling block
    if((bool)(*_handling).exists("handling") && (bool)((*_handling)["handling"]).exists("elim_src"))
    {
       for( int i=0; i<(*_handling)["handling"]["elim_src"].getLength(); ++i )
       {
          string src = (*_handling)["handling"]["elim_src"][i];
	   
	  for (int j=0;j<src_names_ivs.size();j++)
	    {
	      if (src == src_names_ivs.at(j))
		src=src_names.at(j);
	    }
          // check whether source has been used in this session
          vector<string>::iterator it1 = find( src_names.begin(), src_names.end(), src );
          // check whether sources was already inlcuded in general elim_src block
          vector<string>::iterator it2 = find( src_elim.begin(), src_elim.end(), src );
          if( it1 != src_names.end() && it2 == src_elim.end() )
          {
             src_elim.push_back( *it1 );
             log<INFO>("*** Data elimination: Eliminate source ") % *it1 % " because of elim_src in handling";
            vector<int> idx = _param_list.get_indexes( assignment[ "sources" ], src );
            rem_par.insert( rem_par.end(), idx.begin(), idx.end() );
          }
       }
    }
    
    // in case of keeping baselines (>0), everything except the "keeps" have to be removed
    if((*_setup)["keep_baseline"].getLength() > 0)
    {
        //generate baseline-vector
        vector<string> baselines;
        for( int i=0; i<(*_setup)["keep_baseline"].getLength(); ++i )
            baselines.push_back((*_setup)["keep_baseline"][i]);
        
        vector<int> keep_par, keep_obs;
        for( auto &bl: baselines)
        {
            for( auto &scan: _scans )
            { 
               for( int j=0; j<scan.get_nobs(); j++ )
               { 
                   string sta1,sta2,src;
                   scan.get_obs_ptr(j)->get_station_names( sta1, sta2 );
                   scan.get_obs_ptr(j)->get_source_name( src );
                   
                   // using keep_baseline with "WETTZELL-*" alls baselines to WETTZELL will be used
                   bool starform = bl.find("*") != std::string::npos;
                   if( sta1+"-"+sta2 == bl || sta2+"-"+sta1 == bl || (starform == true && (sta1+"-*" == bl || sta2+"-*" == bl)) )
                   {
                        vector<int> idx_src = _param_list.get_indexes( assignment[ "sources" ], src );
                        vector<int> idx_sta1 = _param_list.get_indexes( assignment[ "stations" ], sta1 );
                        vector<int> idx_sta2 = _param_list.get_indexes( assignment[ "stations" ], sta2 );
                        
                        keep_par.insert( keep_par.end(), idx_src.begin(), idx_src.end() );
                        keep_par.insert( keep_par.end(), idx_sta1.begin(), idx_sta1.end() );
                        keep_par.insert( keep_par.end(), idx_sta2.begin(), idx_sta2.end() );
                   }
                   else if(find(baselines.begin(), baselines.end(), sta1+"-"+sta2) == baselines.end())
                       bl_elim.push_back(sta1+"-"+sta2);
               }
            }
        }
        remove_duplicates(keep_par);
        remove_duplicates(bl_elim);
        for(int i=0; i<_param_list.size(); i++)
        {
            bool remove = (find(keep_par.begin(), keep_par.end(), i) == keep_par.end()); 
            if(remove && !_param_list.get_param(i)->is_type({xpo, ypo, ut1, nutx, nuty},{0, 1}))
                rem_par.push_back(i);
        }
    }
    
    // eliminate baseline because of eliminate block
    if((*_setup)["elim_baseline"].getLength() > 0)
    {
        //generate baseline-vector
        vector<string> baselines;
        for( int i=0; i<(*_setup)["elim_baseline"].getLength(); ++i )
        {
            bl_elim.push_back((*_setup)["elim_baseline"][i]);
            log<INFO>("*** Data elimination: Eliminate baseline ") % bl_elim.back() % " because of elim_baseline in elim-block";
        }
    }
    // (2) from DB handling block
    if((bool)(*_handling).exists("handling") && (bool)((*_handling)["handling"]).exists("elim_baseline"))
    {
       for( int i=0; i<(*_handling)["handling"]["elim_baseline"].getLength(); ++i )
       {
          string bl = (*_handling)["handling"]["elim_baseline"][i];
          // check whether baseline was already inlcuded from general sta_elim block
          vector<string>::iterator it1 = find( bl_elim.begin(), bl_elim.end(), bl );
          if( it1 == bl_elim.end() )
          {
             bl_elim.push_back( bl );
             log<INFO>("*** Data elimination: Eliminate baseline ") % bl % " because of elim_baseline in handling";
          }
       }
    }
    //remove_duplicates(rem_par);
    sort( rem_par.begin(), rem_par.end() );
    remove_duplicates(rem_par);
    for( int i=rem_par.size()-1; i>=0; --i )
      {
       _param_list.remove_param( rem_par.at( i ) );
      }
    // loop over scans and observations therein and remove observations
    // remember indexes to remove them later from Lsa-object
    map<string,int> reason_counter;
    vector<int> rem_obs;
    int cur_idx = 0;
    string sta1, sta2, src;
  
    for( int i=0; i <= _scans.size()-1; i++ )
    { 
       int nobs = _scans.at(i).get_nobs();
       cur_idx += nobs-1;  // index of last obs in current scan
       for( int j=_scans.at(i).get_nobs()-1; j>=0; --j )
       { 
          _scans.at(i).get_obs_ptr(j)->get_station_names( sta1, sta2 );
          //_scans.at(i).get_obs_ptr(j)->get_source_name( src );
	  src=_scans.at(i).get_source()->get_name( ivg::srcname::iers);
	  string bl12 = sta1+"-"+sta2;
          string bl21 = sta2+"-"+sta1;

          vector<string>::iterator it1 = find( sta_elim.begin(), sta_elim.end(), sta1 );
          vector<string>::iterator it2 = find( sta_elim.begin(), sta_elim.end(), sta2 );
          vector<string>::iterator it3 = find( src_elim.begin(), src_elim.end(), src );
          vector<string>::iterator it4 = find( bl_elim.begin(), bl_elim.end(), bl12 );
          vector<string>::iterator it5 = find( bl_elim.begin(), bl_elim.end(), bl21 );

	  
	  
          if(!_scans.at(i).get_obs_ptr(j)->get_use_flag())
             reason_counter["use_flag"]++;

          if( it1 != sta_elim.end() || it2 != sta_elim.end() || it3 != src_elim.end() || it4 != bl_elim.end() || it5 != bl_elim.end() ||
              outlier_indicator( cur_idx ) != 1.0 || !_scans.at(i).get_obs_ptr(j)->get_use_flag() )
          {
             // in case of ambiguity resolving we also need to eliminate single band observations
             if(_ambigRes)
             {
                // observation related to group_delay
                _origin_obs_idxs.at(2*cur_idx) = -1;
                rem_obs.push_back(2*cur_idx);
                // observation related to singleband_delay
                _origin_obs_idxs.at(2*cur_idx+1) = -1;
                rem_obs.push_back(2*cur_idx+1);
             }
             else
             {
                _origin_obs_idxs.at(cur_idx) = -1;
               rem_obs.push_back( cur_idx );
             }

             
             _scans.at(i).rem_obs(j);
          }

          // index of prior observation; after for loop it is the index of the last observation in the previous scan
          cur_idx--;  
       }
       cur_idx += nobs+1;    // index of 1st obs in next scan
    }
    sort( rem_obs.begin(), rem_obs.end() );
 //   rem_obs.erase( unique( rem_obs.begin(), rem_obs.end() ), rem_obs.end() );

    int new_nobs = _nobs - rem_obs.size();
    log<INFO>("*** Data elimination: ") % rem_obs.size() % " of " % _nobs % " observations removed. New #nobs " % new_nobs;
    log<INFO>("*** Data elimination: ") % reason_counter["use_flag"] % " removed due to use_flag == false";

    if( new_nobs == 0 )
        throw runtime_error( "ERROR: void Session::_eliminate_data(): No observations left (#obs = 0)" );        
    
    // re-calculate number of observations, stations and sources
    _nobs -= rem_obs.size();

    // get original positions of the observations after observations have been eliminated
    sort( _origin_obs_idxs.begin(), _origin_obs_idxs.end() );
    _origin_obs_idxs.erase( unique( _origin_obs_idxs.begin(), _origin_obs_idxs.end() ), _origin_obs_idxs.end() );
    if( _origin_obs_idxs.at( 0 ) == -1 )
         _origin_obs_idxs.erase( _origin_obs_idxs.begin() );

   
   // remove rows (observations) and columns (parameters) from Lsa-object
   _lsa_solution.remove_data( rem_obs,rem_par );
   _aprioris.rem_r( rem_obs );
   _aprioris.rem_c( rem_par );
      
   // we also (SHOULD) need to remove the station from _trf in case of elim_sta to be able to calculate correct network volume
   // ATTENTION! DOESN'T WORK WITH -r BECAUSE OF STATION POINTER IN SCAN
//   for(vector<ivg::Analysis_station>::iterator ac = _trf.begin(); ac != _trf.end(); ac++ )
//   {
//       if(find( sta_elim.begin(), sta_elim.end(), (*ac).get_name(ivg::staname::ivs_name)) != sta_elim.end())
//       {
//           log<INFO>("!!! Removing ") % (*ac).get_name(ivg::staname::ivs_name) % " from _trf due to elim_sta";
//           _trf.remove_station(ac);
//           ac--;
//       }
//   }
   
    for(std::string& sta: sta_elim){
       _trf.remove_from_station_twin_list(sta);
    }
    _trf.refresh_station_twin_list();
   
    // check if max_obs defined in configfile is reached or not
    if(new_nobs > (int)(*_setup)["max_obs"])
    {
        stringstream ss;
        ss << "void Session::_eliminate_data(): Maximum number of observation reached after elimination(" << new_nobs << ">" << (int)(*_setup)["max_obs"] << ")";
        throw runtime_error(ss.str());
    }
   
#if DEBUG_VLBI >=2
    cerr << "--- void Session::_eliminate_data()" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
void Session::_eliminate_obs_period()
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Session::_eliminate_obs_period( ivg::Date, ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif

    std::vector<ivg::Date> start;
    std::vector<ivg::Date> end;
   
    for( int i=0; i<(*_handling)["handling"]["elim_obs_period"].getLength(); ++i )
    {
        string start_time = (*_handling)["handling"]["elim_obs_period"][i]["start"];
        ivg::Date start_epo( std::stoi( start_time.substr( 0,4 ) ),std::stoi( start_time.substr( 5,2 ) ),
                             std::stoi( start_time.substr( 8,2 ) ),std::stoi( start_time.substr( 11,2 ) ),
                             std::stoi( start_time.substr( 14,2 ) ),s2d( start_time.substr( 17,2 ) ) ); 

        start.push_back( start_epo );

        string end_time = (*_handling)["handling"]["elim_obs_period"][i]["end"];
        ivg::Date end_epo( std::stoi( end_time.substr( 0,4 ) ),std::stoi( end_time.substr( 5,2 ) ),
                           std::stoi( end_time.substr( 8,2 ) ),std::stoi( end_time.substr( 11,2 ) ),
                           std::stoi( end_time.substr( 14,2 ) ),s2d( end_time.substr( 17,2 ) ) ); 

        end.push_back(end_epo);
    }

    // loop over scans and observations therein and remove observations
    // remember indexes to remove them later form Lsa-object
    vector<int> rem_obs;
    for( int k=0; k<start.size(); ++k )
    {
        int cur_idx = 0;
        for( int i=0; i <= _scans.size()-1; i++ )
        { 
            int nobs = _scans.at(i).get_nobs();
            cur_idx += nobs-1;  // index of last obs in current scan
            for( int j=_scans.at(i).get_nobs()-1; j>=0; --j )
            { 
                if( _scans.at(i).get_obs_ptr(j)->get_epoch().get_double_mjd() > start.at(k).get_double_mjd() 
                    && _scans.at(i).get_obs_ptr(j)->get_epoch().get_double_mjd() < end.at(k).get_double_mjd() )
                {
                    _origin_obs_idxs.at(cur_idx) = -1;
                    rem_obs.push_back( cur_idx );
                    _scans.at(i).rem_obs(j);
                }
                // index of prior observation; after for loop it is the index of the last observation in the previous scan
                cur_idx--;  
            }
            // index of 1st obs in next scan
            cur_idx += nobs+1;    
        }
    }
    
    sort( rem_obs.begin(), rem_obs.end() );

    int new_nobs = _nobs - rem_obs.size();
    log<INFO>("*** Data elimination: ") % rem_obs.size() % " of " % _nobs % " observations. New #nobs " % new_nobs; 
    for( int k=1; k<=start.size(); ++k )
        log<INFO>("*** interval (") % k % "): " % setprecision(20) % start.at(k-1).get_double_mjd() % " - " % end.at(k-1).get_double_mjd();

    // re-calculate number of observations, stations and sources
    _nobs -= rem_obs.size();

    // get original positions of the observations after observations have been eliminated
    sort( _origin_obs_idxs.begin(), _origin_obs_idxs.end() );
    _origin_obs_idxs.erase( unique( _origin_obs_idxs.begin(), _origin_obs_idxs.end() ), _origin_obs_idxs.end() );
    if( _origin_obs_idxs.at( 0 ) == -1 )
         _origin_obs_idxs.erase( _origin_obs_idxs.begin() );

    // remove rows (observations) and columns (parameters) from Lsa-object
    _lsa_solution.remove_observations( rem_obs );
    _aprioris.rem_r( rem_obs );
     
    // check if max_obs defined in configfile is reached or not
    if( new_nobs > (int)(*_setup)["max_obs"] )
    {
        stringstream ss;
        ss << "void Session::_eliminate_obs_period( ivg::Date start, ivg::Date end ): Maximum number of observation reached after elimination(" 
           << new_nobs << ">" << (int)(*_setup)["max_obs"] << ")";
        throw runtime_error(ss.str());
    }
   
#if DEBUG_VLBI >=2
    cerr << "--- void Session::_eliminate_obs_period( ivg::Date, ivg::Date )" 
         << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
void Session::_apply_elevation_cutoff()
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Session::_eliminate_obs_period( ivg::Date, ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif

    
    
    double elev_cutoff=((double)(*_setup)["elevation_angle_cutoff"])*M_PI/180.0;

    // loop over scans and observations therein and remove observations
    // remember indexes to remove them later form Lsa-object
    vector<int> rem_obs;
    
        int cur_idx = 0;
        for( int i=0; i <= _scans.size()-1; i++ )
        { 
            int nobs = _scans.at(i).get_nobs();
            cur_idx += nobs-1;  // index of last obs in current scan
            for( int j=_scans.at(i).get_nobs()-1; j>=0; --j )
            {
	      ivg::Matrix azel1=_scans.at(i).get_obs_ptr(j)->get_az_el(1);
	      ivg::Matrix azel2=_scans.at(i).get_obs_ptr(j)->get_az_el(2);
	      
	      if( (azel1(1)<elev_cutoff) || (azel2(1)<elev_cutoff))
                {
                    _origin_obs_idxs.at(cur_idx) = -1;
                    rem_obs.push_back( cur_idx );
                    _scans.at(i).rem_obs(j);
                }
                // index of prior observation; after for loop it is the index of the last observation in the previous scan
                cur_idx--;  
            }
            // index of 1st obs in next scan
            cur_idx += nobs+1;    
        }
    
    
    sort( rem_obs.begin(), rem_obs.end() );
    
    int new_nobs = _nobs - rem_obs.size();
    log<INFO>("*** Data elimination: ") % rem_obs.size() % " of " % _nobs % " observations. New #nobs " % new_nobs; 
   
    log<INFO>("*** Removed data below ") % (elev_cutoff*180/M_PI) % " degrees elevation angle";

    // re-calculate number of observations, stations and sources
    _nobs -= rem_obs.size();

    // get original positions of the observations after observations have been eliminated
    sort( _origin_obs_idxs.begin(), _origin_obs_idxs.end() );
    _origin_obs_idxs.erase( unique( _origin_obs_idxs.begin(), _origin_obs_idxs.end() ), _origin_obs_idxs.end() );
    if( _origin_obs_idxs.at( 0 ) == -1 )
         _origin_obs_idxs.erase( _origin_obs_idxs.begin() );

    // remove rows (observations) and columns (parameters) from Lsa-object
    _lsa_solution.remove_observations( rem_obs );
    _aprioris.rem_r( rem_obs );
     
    // check if max_obs defined in configfile is reached or not
    if( new_nobs > (int)(*_setup)["max_obs"] )
    {
        stringstream ss;
        ss << "void Session::_eliminate_obs_period( ivg::Date start, ivg::Date end ): Maximum number of observation reached after elimination(" 
           << new_nobs << ">" << (int)(*_setup)["max_obs"] << ")";
        throw runtime_error(ss.str());
    }
   
#if DEBUG_VLBI >=2
    cerr << "--- void Session:_apply_elevation_cutoff()" 
         << " : " << tim.toc() << " s " << endl;
#endif 
}


// ...........................................................................
string Session::_create_configuration_textblock(Setting *setup)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Session::_create_configuration_textblock(Setting *)" << endl; 
   tictoc tim;
   tim.tic();
#endif

   Setting &S = (*setup);
   
   stringstream ss;
   ss << " ----------------------- " << endl;
   ss << " Stacking Solution Setup " << endl;
   ss << " ----------------------- " << endl;
   ss << " Configuration: " << (const char *)get_list_element(S["datadirs"],S["session_type"])[3] << endl;
   ss << " ----------------------- " << endl;
   ss << " apriori TRF: " << (const char *)S["trf"] << endl;
   ss << " apriori CRF: " << (const char *)S["crf"] << endl;
   ss << " apriori EOP: " << (const char *)S["eop"]["erp_aprioris"] << endl;
   ss << " ----------------------- " << endl;
   ss << " Used Sessions: " << endl;
   
   for( int i_sess=0; i_sess < S[ "sessions" ].getLength(); ++i_sess )
    ss << " " << (const char *)S[ "sessions" ][ i_sess ][ "dbname" ] << " of " << (const char *)S[ "sessions" ][ i_sess ][ "version" ] << endl;
   
   ss << " ----------------------- " << endl;
   
   vector<string> station_names,source_names;
   if((bool)S["no_net_cnstr"]["stations"]["apply"])
   {
       int idx = S["no_net_cnstr"][ "stations" ][ "stations" ];
       
       if( idx == 0 )
       {
           ss << " NNR/NNT on ALL stations: " << endl;
           ss << "--------------------------" << endl;
           station_names = _trf.get_station_names(ivg::staname::ivs_name);
           std::copy(station_names.begin(), station_names.end(), std::ostream_iterator<string>(ss, "\n"));
       }
       else if (idx>0)
       {
            ss << " NNR/NNT on following stations: " << endl;
            ss << "--------------------------------" << endl << " ";
            for( int i=0; i<S["groups"]["stations"][ idx-1 ].getLength(); ++i )
            {
                string sta_name = S["groups"]["stations"][idx-1 ][ i ];
                ivg::Analysis_station *tmp;
                if(_trf.get_station(&tmp, sta_name, ivg::staname::lettercode))
                {
                    ss << (const char *)S["groups"]["stations"][idx-1 ][ i ] << "\n ";
                }
            }
       }
       else
       {
            ss << " NNR/NNT on all except the following stations: " << endl;
            ss << "--------------------------------" << endl << " ";
            for( int i=0; i<S["groups"]["stations"][ -idx-1 ].getLength(); ++i )
            {
                string sta_name = S["groups"]["stations"][-idx-1 ][ i ];
                ivg::Analysis_station *tmp;
                if(_trf.get_station(&tmp, sta_name, ivg::staname::lettercode))
                {
                    ss << (const char *)S["groups"]["stations"][-idx-1 ][ i ] << "\n ";
                }
            }
       }
       ss << endl;
       ss << "-----------------------------------------------" << endl;
   }
   
    if((bool)S["no_net_cnstr"]["sources"]["apply"])
    {
       int idx = S["no_net_cnstr"][ "sources" ][ "sources" ];
       
       if( idx == 0 )
       {
           ss << " NNR on ALL sources: " << endl;
           ss << "---------------------" << endl << " ";
           source_names = _crf.get_source_names(ivg::srcname::ivs);
           std::copy(source_names.begin(), source_names.end()-1, std::ostream_iterator<string>(ss, " "));
       }
       else
       {
            ss << " NNR on following sources: " << endl;
            ss << "---------------------------" << endl << " ";
            int src_cnt=0;
            for( int i=0; i<S["groups"]["sources"][ idx-1 ].getLength(); ++i )
            {
                string src_name = S["groups"]["sources"][idx-1 ][ i ];
                ivg::Source *tmp;
                if(_crf.get_source(&tmp, src_name))
                {
                    if(src_cnt == 5)
                    {
                        ss << endl << "| ";
                        src_cnt = 0;
                    }
                    src_cnt++;
                    
                    if(tmp->is_defining())
                        ss << "D" << (const char *)S["groups"]["sources"][idx-1 ][ i ] << " ";
                    else if(tmp->is_special_handling())
                        ss << "S" << (const char *)S["groups"]["sources"][idx-1 ][ i ] << " ";
                    else
                        ss << "N" << (const char *)S["groups"]["sources"][idx-1 ][ i ] << " ";
                }
            }
       }
       
       ss << endl;
    }
   
    return ss.str();
   
#if DEBUG_VLBI >=2
    cerr << "--- void Session::_create_configuration_textblock(Setting *)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
ivg::Matrix Session::_model_turbulence( std::string turbulence_model, std::map< std::string, turbulence_data > turb_sta )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ ivg::Matrix Session::_model_turbulence( std::string, std::map< std::string, turbulence_data > )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
   
    // TURBULENCE MODEL; get additional variance-covariance matrix due to atmospheric turbulence / refractive fluctuations
    ivg::Matrix Q_turb;
    ivg::Matrix coeffs(3,1,1.0);		
    ivg::Matrix C;
    double c;

    // set turbulence data per station and spectral coefficients (in case of SIGMA-C or MATERN model)
    _turbulence.set_turb_params( turb_sta );

    if( turbulence_model == "matern_model" || turbulence_model == "matern_sta_model" || turbulence_model == "sigma_c_model" )
    {
       for( int j=0;j<(*_setup)["STOCHASTIC_MODEL"][ "turbulence" ][ "spec_coeffs" ].getLength();++j )
       {
          c = (double)(*_setup)[ "STOCHASTIC_MODEL" ]["turbulence"][ "spec_coeffs" ][ j ];
          coeffs(j,0) = c;
       }  
    
       _turbulence.set_spectral_coefficients( coeffs );
    }

    // define turbulence model	
    if( turbulence_model == "matern_model" )
    {
       Q_turb = _turbulence.calc_matern_vcm_model( _nobs, _scans, &_trf, coeffs, C );
    }
    else if( turbulence_model == "matern_sta_model" )
    {
       Q_turb = _turbulence.calc_station_matern_vcm_model( _nobs, _scans, &_trf, coeffs, C );
    }
    else if( turbulence_model == "sigma_c_model" )
    {
       Q_turb = _turbulence.calc_sigma_c_model( &_trf, _scans, _nobs, coeffs, 200.0, C ) ;
    }
    else if( turbulence_model == "onsala_model" )
    {
       Q_turb = _turbulence.calc_vcm_onsala_model( &_trf, _scans, _nobs, 200.0, 2.0, C ); 
    }
    else if( turbulence_model == "vienna_model" )
    {
       Q_turb = _turbulence.calc_vcm_vienna_model( &_trf, _scans, _nobs, 200.0, 24.0, C ); 
    }    
    else
       throw runtime_error( "Session::_model_turbulence( string turbulence_model ): Selected turbulence model not supported!" );
    
//      cerr << "****** SESSION:::: Q aus MATERN MODEL:  " << endl;
//      Q_turb.get_sub(0,0,10,10).show();  
    
    // scale: mm^2 to s^2
    Q_turb = Q_turb * ( 1e-6 / pow(ivg::c,2.0) );

//      cerr << "****** SESSION:::: Q aus MATERN MODEL ! scaled ! :  " << endl;
//      Q_turb.get_sub(0,0,10,10).show();
      
#if DEBUG_VLBI >=2
    cerr << "--- ivg::Matrix Session::_model_turbulence( std::string, std::map< std::string, turbulence_data > )" << " : " << tim.toc() << " s " << endl;
#endif 

    return Q_turb;
}

// ...........................................................................
void Session::_create_weight_matrix( ivg::Matrix *wgt_ptr )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ ivg::Matrix Session::_create_weight_matrix( ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
   Setting &stoch_model = (*_setup)[ "STOCHASTIC_MODEL" ];
   ivg::Matrix sigma_squ = wgt_ptr->pow(-1.0);

   bool diagonal_vcm = true;
    
   // ADDITIONAL NOISE; Gipson model of station dependent noise
   if( (bool)(stoch_model)["additional_noise"][0]["apply"] == true )
   {
      // prepare VCM for off-diagonal terms
      sigma_squ.to_diag_matrix();

      vector<string>sta_names = _trf.get_station_names( ivg::staname::ivs_name );
      map< string, map< string,double > > add_variances;

      // initialize additional variances to zero for all stations
      for( int i=0;i<sta_names.size();++i )
      { 
         add_variances[ sta_names.at( i ) ][ "cl" ] = 0.0;
         add_variances[ sta_names.at( i ) ][ "at" ] = 0.0;
      }

      // loop over additional noise entries in cnf-file and modify variance-map
      vector<string> selected_sta;
      for( int j=0; j < stoch_model[ "additional_noise" ].getLength(); ++j )
      {
         if( (bool)(stoch_model)["additional_noise"][j]["apply"] )
         {
            // set flag wether off-diagonal terms should be included 
            if( (bool)(stoch_model)["additional_noise"][j]["correlated"] )
               diagonal_vcm = false;

            // group index 0  -> all stations, source or what ever
            //             >0 -> this group
            //             <0 -> all but this group
            int name_idx = stoch_model[ "additional_noise" ][ j ][ "stations" ];
	    selected_sta.clear();
            // aply settings for all stations
            if( name_idx == 0)
            {
               selected_sta = sta_names;
            }
            // select only those param-names within the requested group
            else if( name_idx > 0 )
            {
                for( int i=0;
                        i<(*_setup)["groups"][ "stations" ][ name_idx-1 ].getLength(); ++i )
                    selected_sta.push_back( (*_setup)["groups"][ "stations" ][ name_idx-1 ][ i ] );
            }
            // remove entries of given group from the entire list
            else if( name_idx < 0)
            {
                name_idx *= -1;
                selected_sta = sta_names;
                for( int i=0; i<(*_setup)["groups"][ "stations" ][ name_idx-1].getLength(); ++i )
                {
                    vector<string>::iterator iter = find( selected_sta.begin(), selected_sta.end(),
                                                          (const char *)(*_setup)["groups"][ "stations" ][ name_idx-1][ i ] );
                    if (iter != sta_names.end())
                        selected_sta.erase( iter );
                }
            }

            for( int i=0;i<selected_sta.size();++i )
            {
	      
               add_variances[ selected_sta.at( i ) ][ "cl" ] 
                  = pow( (double)(stoch_model)["additional_noise"][j]["constant_sigma"] * 1e-12,2.0 );
               add_variances[ selected_sta.at( i ) ][ "at" ]
                  = pow( (double)(stoch_model)["additional_noise"][j]["elevation_dependent_sigma"] * 1e-12,2.0 );
            }
         } // if( apply )
      } // loop over entries of additional_noise 


      // loop over observations and add variances as well as covariances within each scan
      int obs_idx = -1;
      string sta1, sta2, sta3, sta4;
      for( int i=0; i < _scans.size(); ++i )                  // loop over scans
      {        
         for( int j=0; j < _scans.at(i).get_nobs(); ++j )     // loop over obs in scan
         { 
            obs_idx++;

            // set variance
            ivg::Obs *obs_ptr = _scans.at(i).get_obs_ptr(j);
            obs_ptr->get_station_names( sta1,sta2 );
                  
            ivg::Analysis_station * stat1 = _scans.at(i).get_sta_ptr( obs_ptr->get_scan_idx(1) );
            ivg::Analysis_station * stat2 = _scans.at(i).get_sta_ptr( obs_ptr->get_scan_idx(2) );
            
            // calculate azimuth and elevation 
            ivg::Matrix AzEl1,AzEl2;
            // if source is a regular source
            if(_scans.at(i).get_source()->get_type() == ivg::srctype::source)
            {
                AzEl1 = stat1->calc_az_el( _scans.at(i).get_epoch(), 
                                       _scans.at(i).get_source()->get_unit_vector_ssb(),
                                       ( _scans.at(i).get_trf2crf_matrix() ).transpose() );
                AzEl2 = stat2->calc_az_el( _scans.at(i).get_epoch(), 
                                       _scans.at(i).get_source()->get_unit_vector_ssb(),
                                       ( _scans.at(i).get_trf2crf_matrix() ).transpose() );
            }
            else if(_scans.at(i).get_source()->get_type() == ivg::srctype::moon || _scans.at(i).get_source()->get_type() == ivg::srctype::satellite )
            // in case of satellite or moon observations
            {
                ivg::Matrix sat_crs = _scans.at(i).get_source()->get_vector_crs(_scans.at(i).get_epoch());
                AzEl1 = stat1->calc_az_el( _scans.at(i).get_epoch(), sat_crs, ( _scans.at(i).get_trf2crf_matrix() ).transpose() );
                AzEl2 = stat2->calc_az_el( _scans.at(i).get_epoch(), sat_crs, ( _scans.at(i).get_trf2crf_matrix() ).transpose() );
            }
                    
            double mf1 = 1.0/sin( AzEl1( 1 ) );
            double mf2 = 1.0/sin( AzEl2( 1 ) );

            sigma_squ( obs_idx,obs_idx ) += add_variances[ sta1 ][ "cl" ]+add_variances[ sta1 ][ "at" ]*pow( mf1,2.0 )
                                           +add_variances[ sta2 ][ "cl" ]+add_variances[ sta2 ][ "at" ]*pow( mf2,2.0 );

            if(diagonal_vcm == false)
            {
                // set covariances
                for( int k=0; k < _scans.at(i).get_nobs()-j-1; ++k )     // loop over remaining obs in scan
                { 
                   int idx = obs_idx+k+1;

                   ivg::Obs *obs2_ptr = _scans.at(i).get_obs_ptr(j+k+1);
                   obs2_ptr->get_station_names( sta3,sta4 );
                   double cov = 0.0;
                   if( sta3 == sta1 || sta4 == sta1 )
                      cov = add_variances[ sta1 ][ "cl" ]+add_variances[ sta1 ][ "at" ]*pow( mf1,2.0 );
                   if( sta3 == sta2 || sta4 == sta2 )
                      cov = add_variances[ sta2 ][ "cl" ]+add_variances[ sta2 ][ "at" ]*pow( mf2,2.0 );

                   sigma_squ( idx,obs_idx ) = cov;
                   sigma_squ( obs_idx,idx ) = cov;
                }
            }
         }
      }
   }

   ivg::Matrix Q;      // weight matrix to be returned
   ivg::Matrix Q_turb; // additional VCM due to atmospheric turbulence


   if( (bool)(stoch_model)["turbulence"]["apply"] && (const char *)(*_setup)["solver"] == "LSM" )
   {
      log<INFO>("*** correlation due to atmospheric turbulence (case 3) used ");
      double Cn, H, v, v_dir;
      std::string turbulence_model = (stoch_model)["turbulence"][ "model" ];  

      Setting &groups = (*_setup)[ "groups" ];
      ivg::Matrix vel(1,1,0.0);
      turbulence_data turb_params;
      std::map< std::string, turbulence_data > turb_sta;
   
      for( int j=0;j<stoch_model[ "turbulence" ]["turb_params"].getLength();++j )
      {
         Setting &setup_turb = stoch_model[ "turbulence" ]["turb_params"][ j ];

         std::vector< std::string > grp_names;
         vector<std::string> sta_names = _trf.get_station_names( ivg::staname::ivs_name );
         group2names( setup_turb,groups,"stations",sta_names,grp_names );

         for( int i=0;i<grp_names.size();++i )
         { 
            Cn = (double)setup_turb[ "Cn" ];       // structure constant 
            H = (double)setup_turb[ "H" ];         // effective tropospheric height
            v = (double)setup_turb[ "v" ];         // wind velocity
            vel.resize(2,1,v);
            v_dir = (double)setup_turb[ "v_dir" ]; // wind direction
            
            turb_params = { Cn, H, v, vel, v_dir, 0.0, 1.0, 1.0 }; 
            turb_sta[ grp_names.at( i ) ] = turb_params;
         }
      }

      Q_turb = _model_turbulence( turbulence_model, turb_sta );

      if( diagonal_vcm && sigma_squ.rows() == sigma_squ.cols() )
      {
         ivg::Matrix Qtmp = sigma_squ.diag();
         Q = Qtmp.diag();

         // Q = sigma_squ;
      }
      else
      {
         Q = sigma_squ.diag();
      }
      
//      cerr << "****** Q  " << endl;
//      Q.get_sub(0,0,10,10).show();
//      cerr << "****** TURB  " << endl;
//      Q_turb.get_sub(0,0,10,10).show();

//      double f = 2.0;
//      double theta = 0.8;
//      double gamma = 0.2;
//      Q = ( Q* gamma + Q_turb* theta ) * f;
      Q += Q_turb ;

//      cerr << "****** Q+TURB  " << endl;
//      Q.get_sub(0,0,10,10).show();

      ivg::Matrix I;
      I.eye( Q.rows() );
      I.solve_neq( Q*1e21 );
      *(wgt_ptr) = I*1e21;

//      cerr << "****** inv Q  " << endl;
//      Q.get_sub(0,0,10,10).show();
   }
   else if( diagonal_vcm && sigma_squ.rows() == sigma_squ.cols()  )
   {
      log<INFO>("*** no correlations (case 1) used ");
      // INITIAL weight matrix with measurement errors returned as a vector
      // Q = ( sigma_squ.diag() )^(-1.0);
      sigma_squ.from_diag_matrix(true);
      *(wgt_ptr) = sigma_squ;
   }
   else if( diagonal_vcm && sigma_squ.rows() != sigma_squ.cols()  )
   {
      log<INFO>("*** no correlations (case 2) used ");
      *(wgt_ptr) = ( sigma_squ )^(-1.0);
   }
   else
   {
      log<INFO>("*** correlated additional noise (case 4) used ");
      // full VCM matrix
      ivg::Matrix I;
      I.eye( sigma_squ.rows() );
      I.solve_neq( sigma_squ*1e21 );
      *(wgt_ptr) = I*1e21;   
   }


#if DEBUG_VLBI >=2
   cerr << "--- ivg::Matrix Session::_create_weight_matrix( ivg::Matrix & )" << " : " << tim.toc() << " s " << endl;
#endif 
}


// ...........................................................................
void Session::_create_inequality_constraints( std::vector<int> & at_idx )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Session::_create_inequality_constraints()" << endl; 
   tictoc tim;
   tim.tic();
#endif 
  
    // find all atmopshere parameters (ZWDs)
//    std::vector<int> at_idx, at_sta_idx;    
    std::vector<int> at_sta_idx;    
    std::vector<std::string> stas = _trf.get_station_names( ivg::staname::ivs_name );
    for( int i = 0; i<stas.size(); i++ )
    {  
        at_sta_idx = _param_list.get_idx( ivg::paramtype::zwd, stas.at(i) );
        at_idx.insert( at_idx.end(), at_sta_idx.begin(), at_sta_idx.end() );
    }   
    std::sort( at_idx.begin(), at_idx.end() );
      
    // assemble ineqaulity constraints   
    Matrix B ;
    B.eye( _param_list.size() );
    B *= -1.0 ;
    B = B( ":", at_idx );
    Matrix b( at_idx.size(),1, 0.0 );
  
    // scale constraints for numerical stability
    B *= 1e-9;
    b *= 1e-9;				
   
    Matrix x0( _param_list.size(), 1, 1.0 );
    ivg::Icls icls1( &_lsa_solution, B, b, x0 );
    
    _icls_solution = icls1;
    //_solution = &_icls_solution; 
    
#if DEBUG_VLBI >=2
   cerr << "--- void Session::_create_inequality_constraints()" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
bool Session::_exist_negative_zwd()
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ bool Session::_exist_negative_zwd()" << endl; 
   tictoc tim;
   tim.tic();
#endif 
  
    // find all atmopshere parameters (ZWDs)
    std::vector<int> at_idx, at_sta_idx;    
    std::vector<std::string> stas = _trf.get_station_names( ivg::staname::ivs_name );
    for( int i = 0; i<stas.size(); i++ )
    {  
        at_sta_idx = _param_list.get_idx( ivg::paramtype::zwd, stas.at(i) );
        at_idx.insert( at_idx.end(), at_sta_idx.begin(), at_sta_idx.end() );
    }   
    std::sort( at_idx.begin(), at_idx.end() );
      
    ivg::Matrix x = _lsa_solution.get_parameters();   
    ivg::Matrix x_atm = x( at_idx );
    if( ( x_atm.find_idx( le, 0.0 ) ).empty() )
        return false;
    else
        return true;
    
#if DEBUG_VLBI >=2
   cerr << "--- bool Session::_exist_negative_zwd()" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
std::vector< ivg::Matrix > Session::_create_collocation_assignment()
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ std::vector< ivg::Matrix > Session::_create_collocation_assignment()" << endl; 
   tictoc tim;
   tim.tic();
#endif 
  
//    std::map< std::string, ivg::Matrix > G;
    std::vector<ivg::Matrix> G;
    std::map<std::string,ivg::Matrix> sta_epochs;
    ivg::Matrix Gi( _nobs,_scans.size(),0.0 );
    ivg::Obs* obs_ptr;
    std::string sta1, sta2;
    int i,k;
                   
    for( vector<ivg::Analysis_station>::iterator sta_it = _trf.begin(); sta_it != _trf.end(); ++sta_it )
    {
        i = 0;
        k = 0;
        std::vector<int> empty_scans;
        ivg::Matrix scan_epo( _scans.size(),1,0.0 );
            
        ivg::Matrix Gi( _nobs,_scans.size(),0.0 );

        for(vector<ivg::Scan>::iterator scan_it = _scans.begin(); scan_it != _scans.end(); ++scan_it )
        { 
            scan_epo(k,0) = scan_it->get_epoch().get_double_mjd();
            
            for(int obs_i=0; obs_i < scan_it->get_nobs(); ++obs_i)
            {    
                obs_ptr = scan_it->get_obs_ptr(obs_i);
                obs_ptr->get_station_names( sta1, sta2 );
              
                if( sta_it->get_name( ivg::staname::ivs_name ) == sta1 || 
                    sta_it->get_name( ivg::staname::ivs_name ) == sta2 )
                    Gi( i,k ) = 1.0;

                ++i;
            }     
            if( scan_it->get_nobs() == 0 )
                empty_scans.push_back(k);
            
            ++k;
        }
                
//        Gi.rem_c(empty_scans);
        std::vector<int> idx = Gi.sum_col().find_idx( eq, 0.0 );
        Gi.rem_c(idx);            
    
        scan_epo.rem_r(idx);
        sta_epochs[ sta_it->get_name( ivg::staname::ivs_name) ] = scan_epo;
        _param_list.set_stoch_param_epoch( sta_epochs );
        
//        G[ sta_it->get_name( ivg::staname::ivs_name) ] = Gi; 
        G.push_back(Gi);
        
    }
    
#if DEBUG_VLBI >=2
   cerr << "--- std::vector< ivg::Matrix > Session::_create_collocation_assignment()" 
        << " : " << tim.toc() << " s " << endl;
#endif 
   
   return G;
}

// ...........................................................................
ivg::Matrix Session::_create_correlation_fct()
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ ivg::Matrix _create_correlation_fct();" << endl; 
   tictoc tim;
   tim.tic();
#endif 
  
    bool turbulence = (bool)(*_setup)["STOCHASTIC_MODEL"]["turbulence"]["apply"];
    std::string solver = (const char *)(*_setup)["solver"];
   
    ivg::Matrix t = *( _lsc_solution.get_lsa_ptr()->get_tobs_ptr() );
    t = ( t - t(0) )* 24.0;
    ivg::Matrix ti = t;
    ivg::Matrix c( t.rows(),1,0.0 );
    ivg::Matrix Q;
    ivg::Matrix Q2;
    ivg::Matrix Qtmp;
    ivg::Matrix Qcc(0,0,0.0);     
    ivg::Matrix Qt(0,0,0.0);      
    double phi, alpha, beta;
    double var = 1.0;
    int cnt = 0;
    
    for( int i=0; i<ivg::paramtype::MAXPARAM; ++i )
    {   
        if( solver == "LSC" && (ivg::paramtype) i == ivg::paramtype::zwd && turbulence )
        {
            double Cn, H, v, v_dir;
            Setting &groups = (*_setup)[ "groups" ];
            ivg::Matrix vel(1,1,0.0);
            turbulence_data turb_params;
            std::map< std::string, turbulence_data > turb_sta;

            for( int j=0;j<(*_setup)["STOCHASTIC_MODEL"][ "turbulence" ]["turb_params"].getLength();++j )
            {
               Setting &setup_turb = (*_setup)["STOCHASTIC_MODEL"][ "turbulence" ]["turb_params"][ j ];

               std::vector< std::string > grp_names;
               vector<std::string> sta_names = _trf.get_station_names( ivg::staname::ivs_name );
               group2names( setup_turb,groups,"stations",sta_names,grp_names );

               for( int i=0;i<grp_names.size();++i )
               { 
                  Cn = (double)setup_turb[ "Cn" ];       // structure constant 
                  H = (double)setup_turb[ "H" ];         // effective tropospheric height
                  v = (double)setup_turb[ "v" ];         // wind velocity
                  vel.resize(2,1,v);
                  v_dir = (double)setup_turb[ "v_dir" ]; // wind direction

                  turb_params = { Cn, H, v, vel, v_dir, 0.0, 1.0, 1.0 }; 
                  turb_sta[ grp_names.at( i ) ] = turb_params;
               }
            }

            _turbulence.set_turb_params( turb_sta );

            ivg::Matrix coeffs(3,1,1.0);
            for( int j=0;j<(*_setup)["STOCHASTIC_MODEL"][ "turbulence" ][ "spec_coeffs" ].getLength();++j )
            {
                coeffs(j,0) = (double)(*_setup)[ "STOCHASTIC_MODEL" ]["turbulence"][ "spec_coeffs" ][ j ];
            }  
            _turbulence.set_spectral_coefficients( coeffs );

            ivg::Matrix cc;
            Qtmp = _turbulence.calc_matern_covariance_function( _nobs, _scans, &_trf, coeffs, cc );            
        }     
        
        for( vector<ivg::Analysis_station>::iterator sta_it = _trf.begin(); sta_it != _trf.end(); ++sta_it )
        {                
            if( _param_list.exist_stoch_param( (ivg::paramtype) i , sta_it->get_name(ivg::staname::ivs_name) ) )
            {   
                ivg::Matrix G = _lsc_solution.get_G().at( sta_it-_trf.begin() );
                std::vector<int> idx = G.transpose().sum_col().find_idx( eq, 0.0 );
                G.rem_r(idx);
                
                if( solver == "LSC" && turbulence )
                {
                    Q2 = Qtmp;
                    Q2.rem_r(idx);
                    Q2.rem_c(idx);

                    ivg::Matrix tmp(cnt+G.cols(),cnt+G.cols(),0.0);
                    tmp.set_sub( 0,0,Qt );
                    tmp.set_sub( cnt,cnt,G.transpose()*Q2*G);
                    Qt = tmp;
                }
                
                ti = t;
                ti.rem_r(idx);
                Q.eye( ti.rows() );
                
                if( (ivg::paramtype) i == ivg::paramtype::zwd )
                {
                    alpha = 6.24;
                    beta = 6.48;
                    phi = 0.82;
                }
                else if( (ivg::paramtype) i == ivg::paramtype::clo )
                {
                    alpha = 2.64;
                    beta = 8.64;
                    phi = 0.33;          
                }
                
                ivg::Matrix dt(ti.rows(),1,0.0);
                for( int j=0; j<ti.rows()-1; ++j )
                {                   
                    dt = ti.get_sub(j+1,0,ti.rows()-1,0);
                    dt = (dt - ti.get_sub(0,0,dt.rows()-1,0) ) / 24.0;                    
                    c.resize( dt.rows(),1,0.0 );
                    
                    for( int k=0; k<dt.rows(); ++k )
                        Q(k+j+1,k) = var / cos( phi ) * exp( dt(k,0) * alpha * -1.0 ) 
                                   * cos( dt(k,0) *beta + phi );  
                }
                
                // fill upper right part of the covariance matrix
                ivg::Matrix D = Q.diag();
                Q = Q + Q.transpose() - D.diag();
                
                if( (bool)(*_setup)["lsc_scaling"] )
                {
                    std::vector<double> foo = G.sum_col().get_data_vec();
                    for( int k=0; k<foo.size(); ++k )
                    {
                        if(foo.at(k)>1)
                            G.set_col(k, G.get_col(k)* 1/foo.at(k) );
                    }
                    
//                    ivg::Matrix foo = G.sum_col();
//                    std::vector<int> idxx = foo.find_idx( gt, 1.0 );
//                    ivg::Matrix E = G.transpose() * Q * G;
//                    E.set_sub( idxx, idxx, E.get_sub( idxx, idxx ).sqrt() );                    
                    
                }
                
                // create VCM for stochastic parameters of a certain group (e.g. ZWDs)
                ivg::Matrix tmp(cnt+G.cols(),cnt+G.cols(),0.0);
                tmp.set_sub( 0,0,Qcc );
                tmp.set_sub( cnt,cnt,G.transpose()*Q*G);
                Qcc = tmp;
                cnt = cnt+G.cols();
                
            }
        }
        
        if( solver == "LSC" && (ivg::paramtype) i == ivg::paramtype::zwd && turbulence )
        {
            double scale = pow( 1e-2 / ivg::c, 2.0 );
            Qcc *= scale;            
             
            double scale2 = pow( 1e-3 / ivg::c, 2.0 );
            Qt *= scale2;
            
            Qcc = Qcc + Qt;            
        }
        else if( solver == "LSC" && (ivg::paramtype) i == ivg::paramtype::zwd && !turbulence )
        {
            double scale = pow( 1e-2 / ivg::c, 2.0 );
            Qcc *= scale;            
        }
    }

#if DEBUG_VLBI >=2
   cerr << "--- ivg::Matrix _create_correlation_fct();" 
        << " : " << tim.toc() << " s " << endl;
#endif 
   
   return Qcc;
}

// ...........................................................................
vector<int> Session::get_obs_idxs_in_ncfile(){
// ...........................................................................
    //vector containing the index of the observations as they are in the nc file
    vector<int> idxs_in_ncfile;
    if(_ambigRes)
    {
        idxs_in_ncfile.resize(_nobs/2); 
        // get only every second element and divide is with two
        // Reason: alternating group and singleband delays in obs. we want the index
        // before adding single band delays 
        for(unsigned i = 0; i < _nobs/2; ++i )
        {
            idxs_in_ncfile[i] = (*get_origin_obs_idxs())[2*i]/2; 
        }
    } else {
        idxs_in_ncfile = *get_origin_obs_idxs();
    }
    
    return idxs_in_ncfile;
    
}

// ...........................................................................
vector<int> Session::get_star_formation_indices(std::string station) {
// ...........................................................................

            
    vector<int> idx;
    
    int cur_idx = 0;
    string sta1, sta2;
    for( ivg::Scan& scan : _scans )
    { 
     
       for( int j=0; j < scan.get_nobs(); ++j )
       { 
          scan.get_obs_ptr(j)->get_station_names( sta1, sta2 );
          
          if( sta1 != station && sta2 != station  ){

                // in case of ambiguity resolving we also need to eliminate single band observations
                if(_ambigRes)
                {
                   // observation related to group_delay
                   idx.push_back(2*cur_idx);
                   // observation related to singleband_delay
                   idx.push_back(2*cur_idx+1);
                } else {
                    idx.push_back(cur_idx);
                }

             
          } 

          cur_idx ++;
       }
       
    }
  
    return idx;
       
}

// ...........................................................................
void Session::find_and_mark_unused_sources(  ){
// ...........................................................................
    for(ivg::Source& src : _crf.get_sources() ){

        bool found = false;
        for(ivg::Scan& scan: _scans){
            if(&src == scan.get_source()){
                found = true;
                break;
            }
        }

        if(found == false){
            src.use_me() = false;
//            std::cerr << "unused " << src.get_name(ivg::srcname::ivs) << std::endl;
        }
    }
    
}

// ...........................................................................
void Session::cable_wrap_time_series(std::map<std::string, ivg::Matrix>& wrap,
                                     std::map<std::string, ivg::Matrix>& time){
// ...........................................................................
    
    std::map<std::string, double> previous_az;
    
    for( int sta_idx = 0; sta_idx < _trf.get_number_stations(); ++ sta_idx){
        ivg::Analysis_station* sta = _trf.get_station(sta_idx);
        double neutral_point = (sta->get_antenna_info().azi_max + sta->get_antenna_info().azi_min) * 0.5;

        wrap[ sta->get_name( ivg::staname::ivs_name ) ] = ivg::Matrix(1,1, neutral_point );
        time[ sta->get_name( ivg::staname::ivs_name ) ] = ivg::Matrix(1,1, _start.get_double_mjd()-0.0004 );
               
    }
    

    for( ivg::Scan& scan : _scans )
    {
       ivg::Date obsStart = scan.get_epoch();
       ivg::Date obsEnd = scan.get_epoch();
       obsEnd.add_secs(scan.get_scheduled_duration());
       
       for( ivg::Analysis_station* sta: scan.get_stations()){
           string sta_name = sta->get_name( ivg::staname::ivs_name );
           
           ivg::Matrix azel_start = sta->calc_az_el(obsStart, *scan.get_source() );
           ivg::Matrix azel_end = sta->calc_az_el(obsEnd, *scan.get_source() );
           
           // slew from last to current source
  
           std::string zone = scan.get_cable_wrapzone(sta);

            double new_wrap, d_azi;
            sta->calc_wrap_and_rotangle( wrap[sta_name].back() ,azel_start(0), new_wrap, d_azi );
            double slew  = new_wrap-wrap[sta_name].back();
            
           if( sta->determine_wrap_zone(new_wrap).compare(zone) != 0  ){
                log<WARNING>("!!! wrap not in zone ") % zone % " at " % sta_name % " " %obsStart.get_date_time("HH:MI:SS");
              slew = -sign(slew)*(2*M_PI-abs(slew));
           }
           
           wrap[sta_name].append_rows( wrap[sta_name].back() + slew );
           time[sta_name].append_rows( obsStart.get_double_mjd() );
           
           // slew after observation
           wrap[sta_name].append_rows( wrap[sta_name].back() + (azel_end(0) - azel_start(0)) );
           time[sta_name].append_rows( obsEnd.get_double_mjd() );
           
       }
       
    }
    
    for( int sta_idx = 0; sta_idx < _trf.get_number_stations(); ++ sta_idx){
        ivg::Analysis_station* sta = _trf.get_station(sta_idx);
                    
        wrap[ sta->get_name( ivg::staname::ivs_name ) ]*=rad2d;

    }
}

// ...........................................................................
void Session::intSked2Latex(ivg::Masterfile masterfile, std::map<int, std::vector<double> >& ca,
                                   std::map<int, std::vector<double> >& cr, std::map<int, std::vector<double> >& caa,
                                   std::map<int, std::vector<double> >& car, double objfun){
// ...........................................................................  
    
    
    
    Setting* setup = this->get_setup();
    
    double minScanLength = (*setup)["SKED"]["min_scan"];
    double maxScanLength = (*setup)["SKED"]["max_scan"];
    double min_time_src = (double)(*setup)["SKED"]["min_time_src"];
    
    double snr_x = (double)(*setup)["SKED"]["snr_min_x"];
    double snr_s = (double)(*setup)["SKED"]["snr_min_s"];
    
    double min_ele = (double)(*setup)["SKED"]["min_elevation"];
    double min_dist_2_sun = (double)(*setup)["SKED"]["min_sun_dist"];
    
    double max_time_without_observation = (double)(*setup)["SKED"]["max_time_without_observation"];
    
    int max_tree_level = (int)(*setup)["SKED"]["max_tree_level"];
    bool objective_surface =  (bool)(*setup)["SKED"]["objective_surface"];
    bool objective_relative =  (bool)(*setup)["SKED"]["objective_relative"];
    double objective_nobs_weight =  (double)(*setup)["SKED"]["objective_nobs_weight"];
    
    int tree_level = ca[0].size();
    
    sessinfo info = masterfile.get_session_info( this->_name_in_masterfile );
    ivg::Date start = info.date;
    ivg::Date end = start;
    end.add_secs(info.duration * 3600.0);
    
    std::string sta1 = _trf.get_station(0)->get_name(ivg::staname::ivs_name);
    std::string sta2 = _trf.get_station(1)->get_name(ivg::staname::ivs_name);
    
    std::ofstream ofs; // write only
    std::string outdir = (const char*)(*setup)[ "outdir" ];
    
    double gap = 0.0;
    double runtime = 0.0;
    // parse  GUROBI options
    for( int param_idx = 0; param_idx < (int)(*setup)["SKED"]["gurobi_param"].getLength(); ++param_idx){
        std::string param = (const char *)(*setup)["SKED"]["gurobi_param"][param_idx][0];
        std::string value = (const char *)(*setup)["SKED"]["gurobi_param"][param_idx][1];
        if( param.compare("MIPGap") == 0){
            gap = stod(value);
        } else if (  param.compare("TimeLimit") == 0 ){
            runtime = stod(value);
        }

    }
    
    
    ofs.open(outdir+"/" + this->_name+"_latexSum.tex", ios::out);
    if (!ofs.good()) {
        std::cerr << "can not write latex sked summary file" << std::endl;
    }else{
    
        // setup block
        ofs << "\\begin{table}" << std::endl
            << "\t\\small" << std::endl
            << "\t\\begin{tabular}{llll}" << std::endl
            << "\t\t\\textbf{setup} \\\\ \\midrule" << std::endl
            << "\t\tscan length & " << (int)minScanLength << " - " << (int)maxScanLength << "s & max level & " << max_tree_level << "\\\\" << std::endl
            << "\t\t\\ac{snr}  X  / S & " << (int)snr_x <<"/" << (int)snr_s  <<"& observation weight & " << objective_nobs_weight << "\\\\" << std::endl
            << "\t\tmin distance to the Sun &" <<  (int)min_dist_2_sun << "$ ^{\\circ} $ & maximized & ";

        int case_obj = 0;
        if( objective_relative == false && objective_surface == false ){
            ofs << "absolute coverage \\\\" << std::endl;
            case_obj = 1;
        } else if( objective_relative == true && objective_surface == false ){
            ofs << "relative coverage \\\\" << std::endl;
            case_obj = 2;
        } else if( objective_relative == false && objective_surface == true ){
            ofs << "absolute area coverage \\\\" << std::endl;
            case_obj = 3;
        } else if( objective_relative == true && objective_surface == true ){
            ofs << "relative area coverage \\\\" << std::endl;
            case_obj = 4;
        }

        ofs << "\t\tmin duration  same source  &" << min_time_src/60 << " m & max runtime & " << setprecision(1) << std::fixed << runtime/3600.0 << " h \\\\" << std::endl
            << "\t\tglobal elevation mask &"<< (int)min_ele << "$ ^{\\circ} $ & min gap & " << gap << " \\%\\\\" << std::endl
            << "\t\ttime period (\\ac{ut}) &" << start.get_date_time("HH:MI:SS") << " - " << end.get_date_time("HH:MI:SS") << "& max duration without obs & " << std::fixed << max_time_without_observation << "s \\\\" << std::endl
            << "\t\t\\\\" << std::endl;

        // source block
        ofs << "\t\t\\textbf{source catalog} \\\\ \\midrule" << std::endl
            << "\t\t\\multicolumn{4}{p{\\textwidth}}{" << std::endl
            << "\t\t";
        for(ivg::Source& sou: _crf.get_sources()){
            ofs << sou.get_name(ivg::srcname::ivs) << " ";
        }
        ofs << std::endl << "\t\t}\\\\"
            << std::endl << "\t\t\\\\" << std::endl;

        // sky plot block
        ofs << "\t\t\\textbf{sky plots} \\\\ \\midrule" << std::endl
            << "\t\t\\multicolumn{2}{c}{" << sta1 << "} & \\multicolumn{2}{c}{" << sta2 << "}\\\\" << std::endl
            << "\t\t\\multicolumn{2}{c}{\\includegraphics[height=.35\\textwidth]{figures/experiments/xxx}} & \\multicolumn{2}{c}{\\includegraphics[height=.35\\textwidth]{figures/experiments/xxx}}\\\\ \\\\"  << std::endl;

        // sky coverage block
        ofs << "\t\t\\textbf{sky coverage} \\\\ \\midrule" << std::endl;

        for(int sta = 0; sta < 2; ++sta){
            ofs << "\t\t\\multicolumn{2}{c}{" << std::endl
                << "\t\t\\begin{footnotesize}" << std::endl
                << "\t\t\\addtolength{\\tabcolsep}{-.5em}" << std::endl
                << "\t\t\\begin{tabular}{"<< std::string(tree_level+1-sta, 'r') << "}" << std::endl
                << "\t\t\t\\multicolumn{" << tree_level+1-sta << "}{l}{" << _trf.get_station(sta)->get_name(ivg::staname::ivs_name) << "}\\\\";
            for(int l = 1; l <= tree_level; ++l){
                if( !(l ==1 && sta==1) )
                     ofs << "&";
                if (l <= max_tree_level){
                    ofs << " \\textbf{" << l << "}";
                }
                else{
                    ofs << " " << l;
                } 
            }
            ofs << "\\\\" << std::endl
                << "\t\t\t";
            if(sta==0){
                if (case_obj == 2){
                    ofs << "\\textbf{rel} &";
                } else {
                    ofs << "rel &";
                }
            }
            for(int l = 0; l < tree_level; ++l){
                if(case_obj == 2 && l < max_tree_level ){
                    ofs << "\\textbf{"<< std::fixed << std::setprecision(1) << setw(5) <<  cr[sta][l]*100 << "\\%}";
                } else{
                    ofs << std::fixed << std::setprecision(1) << setw(5) <<  cr[sta][l]*100 << "\\%";
                }
                if(l < tree_level-1 )
                        ofs << "&";
            }
            ofs << "\\\\" << std::endl
                << "\t\t\t";
            if(sta==0){
                if (case_obj == 1){
                    ofs << "\\textbf{abs} &";
                } else {
                    ofs << "abs &";
                }
            }
            for(int l = 0; l < tree_level; ++l){
                if(case_obj == 1 && l < max_tree_level ){
                    ofs << "\\textbf{"<< std::fixed << std::setprecision(1) << setw(5) <<  ca[sta][l]*100 << "\\%}";
                } else{
                    ofs << std::fixed << std::setprecision(1) << setw(5) <<  ca[sta][l]*100 << "\\%";
                }
                    if(l < tree_level-1 )
                        ofs << "&";
            }
            ofs << "\\\\" << std::endl
                << "\t\t\t";
            if(sta==0){
                if (case_obj == 4){
                    ofs << "\\textbf{rel area} &";
                } else {
                    ofs << "rel area &";
                }
            }
            for(int l = 0; l < tree_level; ++l){
                if(case_obj == 4 && l < max_tree_level ){
                    ofs << "\\textbf{"<< std::fixed << std::setprecision(1) << setw(5) <<  car[sta][l]*100 << "\\%}";
                } else{
                    ofs << std::fixed << std::setprecision(1) << setw(5) <<  car[sta][l]*100 << "\\%";
                }
                    if(l < tree_level-1 )
                        ofs << "&";
            }
            ofs << "\\\\" << std::endl
                << "\t\t\t";
            if(sta==0){
                if (case_obj == 3){
                    ofs << "\\textbf{abs area} &";
                } else {
                    ofs << "abs area &";
                }
            }
            for(int l = 0; l < tree_level; ++l){
                if(case_obj == 3 && l < max_tree_level ){
                    ofs << "\\textbf{"<< std::fixed << std::setprecision(1) << setw(5) <<  caa[sta][l]*100 << "\\%}";
                } else{
                    ofs << std::fixed << std::setprecision(1) << setw(5) <<  caa[sta][l]*100 << "\\%";
                }
                    if(l < tree_level-1 )
                        ofs << "&";
            }
            ofs << "\\\\" << std::endl
                << "\t\t\\end{tabular}" << std::endl
                << "\t\t\\addtolength{\\tabcolsep}{+.5em}" << std::endl
                << "\t\t\\end{footnotesize}" << std::endl
                << "\t\t} ";
            if(sta==0){
                ofs << "&" << std::endl;
            } else {
                ofs << "\\\\" << std::endl;
            }
        }
        ofs << "\t\t\\\\" << std::endl;

         double se = sqrt(_lsa_solution.calc_posterior_vfac());  

        // result block
        ofs << "\t\t\\textbf{results} \\\\ \\midrule " << std::endl
            << "\t\tobjective $ \\Phi $ & "<< objfun << "& runtime & s \\\\" << std::endl
            << "\t\t$ \\sigma_{UT1} $ &" << std::setprecision(0) <<
                _param_list.get_param(_param_list.get_index(ivg::paramtype::ut1,"EOP"))->get_standard_deviation()*ivg::param_unit_fac.at(ivg::paramtype::ut1)*1e3/se <<" micros"
            << "& reached gap & \\% \\\\" << std::endl
            << "number of observations & " << this->_nobs << "\\\\" << std::endl
            << "\t\t\\\\" << std::endl;
        // schedule block
        ofs << "\t\t\\textbf{schedule} \\\\ \\midrule" << std::endl
            << "\t\t\\multicolumn{4}{c}{" << std::endl
            << "\t\t\\begin{tabular}{llr@{\\hskip 1em}|@{\\hskip 1em}llr@{\\hskip 1em}|@{\\hskip 1em}llr}" << std::endl
            << "\t\t\ttime & source & dur & time & source & dur & time & source & dur\\\\" << std::endl;
 
        int nr = ceil(this->_nobs/3.0);
        int count = 0;
        for(int r = 0; r < nr; ++r){
            ofs << "\t\t\t";
            for( int c = 0; c <3; ++c){
                int idx = nr*c +r;
                if(_nobs%3==1 && c ==2){
                    idx--;
                }
                if(idx < _nobs && count < _nobs){
                    count++;
                    ivg::Scan& scan = _scans[idx];
                    ofs << scan.get_epoch().get_date_time("HH:MI:SS") << " & "<< scan.get_source()->get_name(ivg::srcname::ivs) << " & " << (int)scan.get_scheduled_duration();
                    if(c<2){
                        ofs << " & ";
                    } else {
                        ofs << " \\\\" << std::endl;
                    }
                }
            }
        }
        ofs << std::endl << "\t\t\\end{tabular}" << std::endl
            << "\t\t}" << std::endl
            << "\t\\end{tabular}" << std::endl
            << "\t\\caption{Schedule of session " << this->_name << "}" << std::endl
            << "\\end{table}" << std::endl;
    
    }
}

// ...........................................................................
std::string Session::get_best_refstation( std::string preferredStation) {
// ...........................................................................

    vector<string> station_names = _trf.get_station_names(ivg::staname::ivs_name);

    ivg::Matrix has_obs_to_all(station_names.size(), 1, 1.0); // true if a station has  observations to all other station
    ivg::Matrix stat = _obsstats["TRF_LEF"];

    // find stations that have not observations to all other stations
    for (unsigned int c = 0; c < stat.cols(); ++c) {
        for (unsigned int r = 0; r < stat.rows(); ++r) {
            
            // only the upper triangle is used
            int cur_stat = stat(r, c);
            if (r > c) {
                cur_stat = stat(c, r);
            }

            if (cur_stat == 0) {
                has_obs_to_all(c) = 0.0;
                break;
            }
            
        }
    }
    
    
    if ( has_obs_to_all.sum_col()(0) == 0){
         log<WARNING>("!!! There is no station with observations to all other stations.");
         return "";
    }
    
    // if there is a preferred station check whether it has observations to all other stations
    if( !preferredStation.empty() ){
        // find index in _obsstats Matrix
         int pos = std::find(station_names.begin(), station_names.end(), preferredStation) - station_names.begin();
         if(pos >= station_names.size()) {
             log<WARNING>("!!! The station ") % preferredStation % " was not found in observation statistics matrix";
        } else {
            if (has_obs_to_all(pos) == true){
                log<INFO>("*** the preferred station: ") % preferredStation % " has observations to all other stations";
                return preferredStation;
            } else {
                log<INFO>("*** the preferred station: ") % preferredStation % " has not observations to all other stations. It is not used";
            }
        }
    }

    // if no preferred station is declared use the one with the most observations
    // get number of observations for each station
    ivg::Matrix tot_obs = stat.diag();

    // the station with the most observations that has observations to all other stations is the best one;
    int best = has_obs_to_all.mult_elem( tot_obs ).maxIdx();
    
    log<INFO>("*** the station with the most observations that has observations to all other stations is: ") % station_names[best];
    
    return station_names[best];
}

  void Session::calc_transits(  unsigned int dt ){
            
    //  transits at each station -----------------------------------------------
    _trasintsVisibleFromStation = lps::Transits (this, dt);

    // transits visible from at least by 2 stations ----------------------------
    _trasintsVisibleFrom2Stations = _trasintsVisibleFromStation;

    ivg::Schedule s(this);

    std::vector<unsigned> visble_from_stas = s.number_of_stations_source_is_visible(_trasintsVisibleFrom2Stations);
    for( unsigned scr_idx = 0; scr_idx < _crf.get_number_sources_inc_unused(); ++scr_idx){
        if(visble_from_stas[scr_idx] < 2)
            _trasintsVisibleFrom2Stations.remove_source(scr_idx);
    }
    
    //_trasintsVisibleFrom2Stations.getBaselineWiseTransits().print_transits();
    
    // transits with at least one baseline with good SNR at one point in time --
    _complete_transits = _trasintsVisibleFrom2Stations;
    std::vector<bool> has_good_snr = s.quasar_has_at_least_one_bl_with_good_SNR(_complete_transits);
    for( unsigned scr_idx = 0; scr_idx < _crf.get_number_sources_inc_unused(); ++scr_idx){
        if(!has_good_snr[scr_idx])
            _complete_transits.remove_source(scr_idx);
    }
    
    // final transits. good snr any time ---------------------------------------
    _common_transits = _complete_transits.seenByAtLeastBaselines( 1 ); // 
};


 void Session::init_session_eop( const Setting& eopblock, ivg::Date date ){
                
        _eops = init_eop( eopblock,  date );

    };

ivg::Eop_series Session::init_eop( const Setting& eopblock, ivg::Date date, double days ){
    
    std::string eopfile = (const char *)get_list_element((*_setup)["eop_files"], eopblock["erp_aprioris"])[2];

    ivg::Eop_series eops( eopfile ,
                  (const char *)get_list_element((*_setup)["eop_files"], eopblock["erp_aprioris"])[1],
                  date.add_days(-days), date.add_days(days) );



    eops.init( eopblock["interpolation_type"], eopblock["ut1_zonal_tides"],
               eopblock["hf_ocean"],
 	       eopblock["hf_ocean_model"],
               eopblock["ut_libration"],eopblock["pm_nutation"],
               eopblock["nutation_type"],
               get_list_element((*_setup)["stadisp"], "POLE TIDE")[2] );
    
    return eops;
}
    
void Session::init_param_list(  ){
    ivg::Date mean = ivg::Date( 0.5*( _start.get_double_mjd() +_end.get_double_mjd()) );
    _param_list = Param_list(_trf, _crf, _eops, mean);
    _param_list.set_start_end_epoch(_start, _end); 
}

void Session::change_eop_series( const Setting& eopblock, ivg::Date date ) {
            
    init_session_eop( eopblock, date );
    init_param_list();

    for(ivg::Scan& scan : _scans) {
        ivg::Partials_t2c tmp { ivg::Matrix( 3,3,0.0 ) };
        ivg::Partials_t2c * deriv_ptr = &tmp;
        ivg::Matrix crf2trf = _eops.form_crf2trf( scan.get_epoch(), true, deriv_ptr );

        scan.set_trf2crf_matrix( crf2trf.transpose(), deriv_ptr );
    }

}

std::map<int, std::vector<double> > Session::compute_sky_coverage(const std::vector<std::vector<unsigned> >& ea_grid_setup,
                                    const lps::TemporalGrid& tg,
                                    std::vector< std::vector< lps::Node<ivg::Obs*, lps::Wedge>  > >&  roots,
                                    double& objective){
    return compute_sky_coverage(ea_grid_setup, tg, this->_scans, roots, objective);
}  


std::map<int, std::vector<double> > Session::compute_sky_coverage(const std::vector<std::vector<unsigned> >& ea_grid_setup,
                                    const lps::TemporalGrid& tg,
                                    std::vector<ivg::Scan>& scans,
                                    std::vector< std::vector< lps::Node<ivg::Obs*, lps::Wedge>  > >&  roots,
                                    double& objective){
    objective = 0;
    std::map<int, std::vector<double> > coverage;
    
    std::vector<unsigned> n_segments( ea_grid_setup.size() ); 
    for(int i=0; i < ea_grid_setup.size(); i++){
        for(int j=0; j<ea_grid_setup[i].size(); j++){
            n_segments[i] += ea_grid_setup[i][j];
        }
    }
    
    // Tree

    std::vector< std::vector<ivg::Analysis_station*> > twinlist = _trf.get_station_twin_list();
    roots.resize( twinlist.size() );
    
    unsigned observatory_idx = 0;
    for(std::vector<ivg::Analysis_station*> observatory : twinlist ){
        
        std::vector< std::vector<lps::DataPoint<ivg::Obs*>> > points;
        points.resize(tg.get_number_of_intervals());
        

        for(ivg::Analysis_station* sta : observatory ){
            coverage[sta->get_idx()].resize(ea_grid_setup.size());
            for (ivg::Scan& scan : scans) {
                for(int obs_idx = 0; obs_idx < scan.get_nobs(); ++obs_idx){
                     ivg::Obs* const obs = scan.get_obs_ptr(obs_idx);
                     string sta1 = "";
                     string sta2 = "";
                     obs->get_station_names(sta1, sta2);

                    if( sta1.compare( sta->get_name(ivg::staname::ivs_name) ) == 0 || sta2.compare( sta->get_name(ivg::staname::ivs_name) ) == 0  ){
                        ivg::Source* const src = scan.get_source();
                        ivg::Date epo = scan.get_epoch();

                        ivg::Matrix azel = sta->calc_az_el(epo, *src);
                        azel*=ivg::rad2d;
                        lps::Position pos ( azel(0), azel(1), epo);

                        lps::Point p (pos.azimuth(),pos.elevation());

                        lps::Seconds start = (epo.get_double_mjd() - _start.get_double_mjd())*3600*24;
                        for( int t : tg.getTemporalIndices(start)){
                            points[t].push_back( lps::DataPoint<ivg::Obs*>(p, obs) );
                        }
                    }
                }
            }

        }

        // build tree
        lps::Wedge wedge (lps::Point(0,0),0,90,0,360);
        roots[observatory_idx].resize(tg.get_number_of_intervals());
        for(unsigned int t = 0; t < tg.get_number_of_intervals(); ++t){

            roots[observatory_idx][t] = lps::Node<ivg::Obs*, lps::Wedge>(points[t], wedge, ea_grid_setup);

            // coverage
            roots[observatory_idx][t].visitBoundAndBool([&]( int level, const lps::Wedge& rect, bool valid ){                
                for(ivg::Analysis_station* sta : observatory ){
                       if( valid ){
                        coverage[sta->get_idx()][level-1] += 1.0/n_segments[level-1];
                        objective += 1.0/n_segments[level-1];
                    }
                }
                
            });
        }
        
        observatory_idx++;
    }

    return coverage;
}

} // # namespace ivg
