#include "obs.h"
#include "session.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Obs::Obs()
// ...........................................................................
{
    _station1 = "";
    _station2 = "";
    _source = "";

    _group_delay      = -999.99;
    _delay_rate       = -999.99;
    _geoc_delay       = -999.99;
    _sigma_delay      = -999.99;
    _sigma_delay_rate = -999.99;
    
    _sb_delay = -999.99;
    _sigma_sb_delay = -999.99;

    _delta_ion       = -999.99;
    _sigma_delta_ion = -999.99;
    _delta_ion_rate  = -999.99;
    _sigma_delta_ion_rate = -999.99;
    
    _oc = 0.0;
    
    _use_me = false;
}


// ...........................................................................
Obs::Obs( ivg::Session * session_ptr, ivg::Scan * scan_ptr,
          int sta1_scan_idx, int sta2_scan_idx )
// ...........................................................................
{
    _session = session_ptr;
    _scan = scan_ptr;
    _sta1_scan_idx = sta1_scan_idx;
    _sta2_scan_idx = sta2_scan_idx;

    _epoch = _scan->_epoch;

    _group_delay      = -999.99;
    _delay_rate       = -999.99;
    _geoc_delay       = -999.99;
    _sigma_delay      = -999.99;
    _sigma_delay_rate = -999.99;
    
    _sb_delay = -999.99;
    _sigma_sb_delay = -999.99;

    _delta_ion       = -999.99;
    _sigma_delta_ion = -999.99;
    _delta_ion_rate  = -999.99;
    _sigma_delta_ion_rate = -999.99;
    
    _use_me = false;
}


// ===========================================================================
// public methods
// ===========================================================================

// ...........................................................................
void Obs::show()
// ...........................................................................
{
    cout << "-----------------OBS-----------------------" << endl;
    cout << setprecision(16) << _epoch.get_double_mjd() << " " << _station1 << " "
         << _station2 << " " << _source << " " << _group_delay << endl;
    cout << setprecision(16)
         << _group_delay << endl;
    _scan->show();
    cout << "-----------------OBS-----------------------" << endl;
}


// ...........................................................................
void Obs::set_delay( double d, double sigmaD, double gcd, double r, double sigmaR )
// ...........................................................................
{
  if (d>0.5)
    _group_delay = d-1;
  else if (d<-0.5)
    _group_delay=d+1;
  else
    _group_delay=d;

  
    _sigma_delay = sigmaD;

    _geoc_delay = gcd;
    
    _delay_rate = r;
    _sigma_delay_rate = sigmaR;
}
// ...........................................................................
void Obs::set_sb_delay( double sbd, double sigmaSbd )
// ...........................................................................
{
    _sb_delay = sbd;
    _sigma_sb_delay = sigmaSbd;
}
// ...........................................................................
void Obs::set_phase_delay( double Phd, double sigmaPhd )
// ...........................................................................
{
    _phase_delay = Phd;
    _sigma_phase_delay = sigmaPhd;
}  
// ...........................................................................
void Obs::set_snr( double snr_bx, double snr_bs )
// ...........................................................................
{
    _snr_bx = snr_bx;
    _snr_bs = snr_bs;
}
// ...........................................................................
void Obs::set_ion_corr( double id, double sigmaId, double ir,
                        double sigmaIr )
// ...........................................................................
{
    _delta_ion = id;
    _sigma_delta_ion = sigmaId;
    _delta_ion_rate = ir;
    _sigma_delta_ion_rate = sigmaIr;
}

// ...........................................................................
void Obs::set_feed_rotation( double feed1, double feed2 )
// ...........................................................................
{
    _scan->_data.at( _sta1_scan_idx ).feed_rotation = feed1;
    _scan->_data.at( _sta2_scan_idx ).feed_rotation = feed2;
}
// ...........................................................................
void Obs::set_cable_cal( double cc1, double cc2 )
// ...........................................................................
{
    _scan->_data.at( _sta1_scan_idx ).cable_cal = cc1;
    _scan->_data.at( _sta2_scan_idx ).cable_cal = cc2;
}


// ...........................................................................
void Obs::set_scan( ivg::Scan* scan )
// ...........................................................................
{
    _scan = scan;
}

// ...........................................................................
void Obs::set_epoch( ivg::Date epoch )
// ...........................................................................
{
    _epoch = epoch;
}

// ...........................................................................
int Obs::get_scan_idx( int idx )
// ...........................................................................
{
    if( idx == 1 )
    {
        return _sta1_scan_idx;
    }
    else if( idx == 2 )
    {
        return _sta2_scan_idx;
    }
    else
    {
        stringstream errormessage;
        errormessage << "int Obs::_get_scan_idx( int idx ): "
                     << "in valid station index: " << idx << "!= (1,2)! Exiting";
        throw runtime_error( errormessage.str() );
    }
}


// ...........................................................................
ivg::Matrix Obs::get_az_el( int idx )
// ...........................................................................
{
    if( idx == 1 )
    {
        return _AzEl1;
    }
    else if( idx == 2 )
    {
        return _AzEl2;
    }
    else
    {
        stringstream errormessage;
        errormessage << "int Obs::get_az_el( int idx ): "
                     << "in valid station index: " << idx << "!= (1,2)! Exiting";
        throw runtime_error( errormessage.str() );
    }
}


// ...........................................................................
ivg::Matrix Obs::get_mapping_functions( int idx )
// ...........................................................................
{
    if( idx == 1 )
    {
        return _mfs1;
    }
    else if( idx == 2 )
    {
        return _mfs2;
    }
    else
    {
        stringstream errormessage;
        errormessage << "int Obs::get_mapping_functions( int idx ): "
                     << "in valid station index: " << idx << "!= (1,2)! Exiting";
        throw runtime_error( errormessage.str() );
    }
}

// ...........................................................................
bool Obs::operator==( const Obs test ) const
// ...........................................................................
{
    bool out = false;
    if( _station1 == test._station1 &&
            _station2 == test._station2
      )
        out = true;
    else if( _station2 == test._station1 &&
             _station1 == test._station2
           )
        out = true;

    return out;
}
// ...........................................................................
void Obs::calc_delay( vector<double>::iterator design_iter, vector<double>::iterator apriori_iter, bool ion, bool singleband, double * par_aplo )
// ...........................................................................
{
  if (singleband)
    calc_delay(design_iter,apriori_iter,ion, 's', par_aplo );
  else
    calc_delay(design_iter,apriori_iter,ion, 'g', par_aplo );
}

// ...........................................................................
void Obs::calc_delay( vector<double>::iterator design_iter, vector<double>::iterator apriori_iter, bool ion, char delaytype, double * par_aplo )
// ...........................................................................
{
   
    tictoc time;
    time.tic();
    
    //   log<DETAIL>("*** Calculating delay for ") % _scan->_source->get_name( ivg::srcname::ivs ) % " at " % _epoch.get_double_mjd();

    // following variables (call-by-reference) are calculated within calc_total_delay.
    ivg::Matrix v_earth_ssb(3, 1);
    ivg::Matrix vel1( 3,1,0.0 );
    ivg::Matrix vel2( 3,1,0.0 );
    ivg::Matrix b_gcrs,b_trs;
    ivg::Matrix k;
    
    double tau = 0.0;
    // if source is a regular radio-source use standard consensus model
    
    if(_scan->_source->get_type() == ivg::srctype::source)
        tau = calc_total_delay(ivg::delaymodel::consensus, v_earth_ssb, vel1, vel2, b_gcrs, b_trs, k);
    // if NEARFIELD-block exists and source is moon or satellite
    else if((*(_session->_setup)).exists("NEARFIELD") && _scan->_source->get_type() != ivg::srctype::source)
    {
        k.resize(9, 1);
        
        if( _scan->_source->get_type() == ivg::srctype::moon )
            tau = calc_total_delay(get_delaymodel((*(_session->_setup))["NEARFIELD"]["moon_model"]), v_earth_ssb, vel1, vel2, b_gcrs, b_trs, k);
        else if( _scan->_source->get_type() == ivg::srctype::satellite )
            tau = calc_total_delay(get_delaymodel((*(_session->_setup))["NEARFIELD"]["sat_model"]), v_earth_ssb, vel1, vel2, b_gcrs, b_trs, k); 
        
        tau *= -1.0;
    }
    else
       throw runtime_error("void Obs::calc_delay( ... ): Unexpected model configuration.");

   
    ivg::Matrix trf2crf = _scan->_trf2crf;
    ivg::Analysis_station * sta1 = _scan->_data.at( _sta1_scan_idx ).sta_ptr;
    ivg::Analysis_station * sta2 = _scan->_data.at( _sta2_scan_idx ).sta_ptr;

    ivg::Matrix dx1_ntaplo = sta1->calc_nontidal_aplo( _epoch, "cspline" );
    ivg::Matrix dx2_ntaplo = sta2->calc_nontidal_aplo( _epoch, "cspline" );
    
    
    ivg::Matrix k1(3, 1);
    ivg::Matrix k2(3, 1);
    
    if ((*(_session->_setup)).exists("NEARFIELD") && _scan->_source->get_type() != ivg::srctype::source)
    {
        // Source vectors k1 and k2 for near-field observations
        ivg::Matrix ktemp(3, 1);
        
        for (int i = 0; i < 3; i++)
            ktemp(i) = k(i);
        
        k1 = ktemp +( v_earth_ssb+vel1 )/ivg::c - ktemp*(ktemp.transpose()*( v_earth_ssb+vel1 )/ivg::c);

        for (int i = 0; i < 3; i++)
            ktemp(i) = k(i+3);
        
        k2 = ktemp +( v_earth_ssb+vel2 )/ivg::c - ktemp*(ktemp.transpose()*( v_earth_ssb+vel2 )/ivg::c);
        
        for (int i = 0; i < 3; i++)
            ktemp(i) = k(i+6);
        
        k.resize(3, 1);
        for (int i = 0; i < 3; i++)
            k(i) = ktemp(i);
    } 
    else
    {    
    // (7) calculate the abberated source vector for use in the calculation of the tropospheric propagation delay. (Eq 11.15)
    // (7.1) annual abberation
        k1 = k+( v_earth_ssb+vel1 )/ivg::c - k*(k.transpose()*( v_earth_ssb+vel1 )/ivg::c);
        k2 = k+( v_earth_ssb+vel2 )/ivg::c - k*(k.transpose()*( v_earth_ssb+vel2 )/ivg::c);
    }
        
    // (7.2) correct for tropospheric bending and transform to TRF calculate azimuth and elevation for both stations
    ivg::Matrix azel1 = sta1->calc_az_el( _epoch, k1, trf2crf.transpose() );
    ivg::Matrix azel2 = sta2->calc_az_el( _epoch, k2, trf2crf.transpose() );
    _AzEl1 = azel1;
    _AzEl2 = azel2;
    
//    log<DETAIL>("*** [Step 7of12] Az and El of 1&2: ") % _AzEl1(0) % "/" % _AzEl1(1) % " | " % _AzEl2(0) % "/" % _AzEl2(1);

    ivg::Matrix k1_trf = trf2crf.transpose()*k1;
    ivg::Matrix k2_trf = trf2crf.transpose()*k2;   
    ivg::bendingmode bend_type=ivg::bendingmode::INSITU;
    if((*(_session->_setup))["troposphere"].exists("bending_mode"))
    {
        std::string bend = (*(_session->_setup))["troposphere"]["bending_mode"];
        if(bend== "none")
            bend_type = ivg::bendingmode::NO;
        else if(bend=="insitu")
            bend_type = ivg::bendingmode::INSITU;
        else if(bend=="approx")
            bend_type = ivg::bendingmode::APPROX;
        else if(bend=="height")
            bend_type = ivg::bendingmode::HEIGHT;
        else
            bend_type = ivg::bendingmode::INSITU;
    }
    if(bend_type != NO)
    {
        ivg::Matrix llh = sta1->calc_lat_lon_h();
        k1_trf = _scan->_data.at(_sta1_scan_idx).tropo.build_apparent_source_vector
            (k1_trf,azel1,llh(2),bend_type);
        llh = sta2->calc_lat_lon_h();
        k2_trf = _scan->_data.at(_sta2_scan_idx).tropo.build_apparent_source_vector
            (k2_trf,azel2,llh(2),bend_type);
    }    
    
    // (8) add the geometric part of the tropospheric propagation delay to the vacuum delay. (Eq 11.11)
    Setting &definitions = ( *(_session->_setup) )["definitions"];
    std::string mf_type = ( *(_session->_setup) )["troposphere"]["mapping_function"];
    std::string interpolation_type = ( *(_session->_setup) )["troposphere"]["interpolation_type"];

    // determine slant hydrostatic tropospheric delay (wet?????)
    double zhd1 = _scan->_data.at( _sta1_scan_idx ).tropo.get_zhd();
    double zhd2 = _scan->_data.at( _sta2_scan_idx ).tropo.get_zhd();

    double zwd1 = _scan->_data.at( _sta1_scan_idx ).tropo.get_zwd();
    double zwd2 = _scan->_data.at( _sta2_scan_idx ).tropo.get_zwd();
    
    // calculate mapping function coefficients
    double mf_hydr1, mf_wet1, mf_hydr2, mf_wet2;
    
    if( mf_type == "gpt2" )
    {
        _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf1_gpt2( azel1(1), mf_hydr1, mf_wet1 );
        _scan->_data.at( _sta2_scan_idx ).tropo.calc_vmf1_gpt2( azel2(1), mf_hydr2, mf_wet2 );
    }
    else if( mf_type == "gpt3" )
    {
        _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf3_gpt3( azel1(1), mf_hydr1, mf_wet1 );
        _scan->_data.at( _sta2_scan_idx ).tropo.calc_vmf3_gpt3( azel2(1), mf_hydr2, mf_wet2 );
    }
    else if( mf_type == "vmf1" )
    {
        // if not enough VMF1 data is available for a station use gpt2
        if(sta1->get_data_status("VMF1") == "W")
            _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf1_gpt2( azel1(1), mf_hydr1, mf_wet1 );
        else
            _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf1( azel1(1), interpolation_type, mf_hydr1, mf_wet1 );
        
        if(sta2->get_data_status("VMF1") == "W")
            _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf1_gpt2( azel2(1), mf_hydr2, mf_wet2 );
        else
            _scan->_data.at( _sta2_scan_idx ).tropo.calc_vmf1( azel2(1), interpolation_type, mf_hydr2, mf_wet2 );
    }
    else if( mf_type == "vmf3" )
    {
        // if not enough VMF1 data is available for a station use gpt2
        if(sta1->get_data_status("VMF3") != "X")
            _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf3_gpt3( azel1(1), mf_hydr1, mf_wet1 );
        else
	  {
	    try {
	      _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf3( azel1(1), interpolation_type, mf_hydr1, mf_wet1 );
	    }
	    catch (int e) {
	       _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf3_gpt3( azel1(1), mf_hydr1, mf_wet1 );
	    }
	  }
        
        if(sta2->get_data_status("VMF3") != "X")
            _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf3_gpt3( azel2(1), mf_hydr2, mf_wet2 );
        else
	  {
	    try {
            _scan->_data.at( _sta2_scan_idx ).tropo.calc_vmf3( azel2(1), interpolation_type, mf_hydr2, mf_wet2 );
	    }
	    catch (int e) {
	      _scan->_data.at( _sta1_scan_idx ).tropo.calc_vmf3_gpt3( azel2(1), mf_hydr2, mf_wet2 );
	    }
	  }
    }
    else if( mf_type == "gmf" )
    {
        _scan->_data.at( _sta1_scan_idx ).tropo.calc_gmf( azel1(1), mf_hydr1, mf_wet1 );
        _scan->_data.at( _sta2_scan_idx ).tropo.calc_gmf( azel2(1), mf_hydr2, mf_wet2 );
    }
    else if( mf_type == "raytrace" )
    {
      double tmp_mfh,tmp_mfw;
      _scan->_data.at( _sta1_scan_idx ).tropo.get_mapping_function(tmp_mfh,tmp_mfw );
      mf_hydr1=tmp_mfw;
      mf_wet1=tmp_mfw;
      _scan->_data.at( _sta2_scan_idx ).tropo.get_mapping_function(tmp_mfh,tmp_mfw );
      mf_hydr2=tmp_mfw;
      mf_wet2=tmp_mfw;
    }
    else
        throw runtime_error("void Obs::calc_delay( vector<double>::iterator design_iter ): Unknown mapping function chosen.");
    
    _mfs1.resize(2,1,0.0);
    _mfs2.resize(2,1,0.0);

    _mfs1(0) = mf_hydr1;
    _mfs1(1) = mf_wet1;
    _mfs2(0) = mf_hydr2;
    _mfs2(1) = mf_wet2;  

    // determine the asymmetric delay caused by gradients as well as east and north gradients
    double dgr1, ngr1, egr1, dgr2, ngr2, egr2;
    
    dgr1 = _scan->_data.at( _sta1_scan_idx ).tropo.calc_gradient_delay( azel1(1), azel1(0) );
    dgr2 = _scan->_data.at( _sta2_scan_idx ).tropo.calc_gradient_delay( azel2(1), azel2(0) );
   
//    _scan->_data.at( _sta1_scan_idx ).tropo.calc_gradients( azel1(1), azel1(0), dgr1, ngr1, egr1 );
//    _scan->_data.at( _sta2_scan_idx ).tropo.calc_gradients( azel1(1), azel1(0), dgr2, ngr2, egr2 );
    
    // total tropospheric delay
    double tau_atm1 = ( mf_hydr1 * zhd1 + mf_wet1 * zwd1 + dgr1 )/ivg::c;
    double tau_atm2 = ( mf_hydr2 * zhd2 + mf_wet2 * zwd2 + dgr2 )/ivg::c;

//    double tau_atm1 = ( mf_hydr1 * zhd1 + dgr1 )/ivg::c;
//    double tau_atm2 = ( mf_hydr2 * zhd2 + dgr2 )/ivg::c;    
    
    tau += tau_atm1*(k.transpose()*(vel2-vel1))(0)/ivg::c;
    
    // (9) compute total delay by adding the best estimate of the tropospheric propagation delay (Eq 11.12)
    tau += ( tau_atm2-tau_atm1 );
   
    // log<DETAIL>("*** [Step 9of12] Total delay: ") % tau;

    // (10) correct time delay from step 9 for "post-model" changes in the baseline. (Eq 11.13)
    // ???????????????????????????
    //double dtau =
    // ???????????????????????????

    // ===================
    // FURTHER CORRECTIONS
    // ===================
    
    // axis offsets
    double ao1 = _scan->calc_axis_offset_delay( _sta1_scan_idx, k1_trf );
    double ao2 = _scan->calc_axis_offset_delay( _sta2_scan_idx, k2_trf );
    tau += ( ao2-ao1 );
   
    // thermal expansion
    double thermal_1 = 0.0;
    double thermal_2 = 0.0;
    if((*(_session->_setup)).exists("thermal_expansion")&&(bool)(*(_session->_setup))["thermal_expansion"])
    {
        thermal_1 = _scan->calc_thermal_expansion_delay( _sta1_scan_idx, k1_trf );
        thermal_2 = _scan->calc_thermal_expansion_delay( _sta2_scan_idx, k2_trf );
    }
    
    tau += ( thermal_1-thermal_2 );
   
    double gravdef_1 = 0.0;
    double gravdef_2 = 0.0;
    if((*(_session->_setup)).exists("gravitational_deformation")&&(bool)(*(_session->_setup))["gravitational_deformation"])
    {
      gravdef_1 = _scan->calc_gravdef_delay( _sta1_scan_idx, k1_trf );
      gravdef_2 = _scan->calc_gravdef_delay( _sta2_scan_idx, k2_trf  );
    }
    tau += ( gravdef_2-gravdef_1 );
   
//    // TEMPORARY LUCIA SATELLITE COMPARISON OUTPUT
//    stringstream tmp_out2;
//    tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << -1.0 * tau << endl;
//    ivg::Logger::add_to_saveline = tmp_out2.str();
    
    // cable calibration
    double cc = _scan->_data.at( _sta2_scan_idx ).cable_cal - _scan->_data.at( _sta1_scan_idx ).cable_cal;
    
    // distinguish between singleband delay and groupdelay
   
    if(delaytype=='s')
    {
        double sbd = _sb_delay;
        sbd += cc;
        _oc = sbd-tau;

    }
    else if (delaytype=='p')
      {
	double phd = _phase_delay; 
	double feedc = _scan->_data.at( _sta2_scan_idx ).feed_rotation - _scan->_data.at( _sta1_scan_idx ).feed_rotation;
        phd += cc+feedc;
        _oc = phd-tau;
	
      }
    else
    {
        double gd = _group_delay;

        // ionosphere correction
        if( ion )
           gd -= _delta_ion;

//        stringstream tmp_out2;
//        tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << _AzEl1(1) ;
//        tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << _AzEl2(1) ;
//        tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << _delta_ion << endl;
//        ivg::Logger::add_to_saveline = tmp_out2.str();          
        
        gd += cc;
        _oc = gd-tau;
	

	//if (_use_me)
	// std::cout << (_epoch.get_double_mjd()-57204)*24 << " " <<_oc << " " << tau << " " << gd << " " << cc <<" "<< _delta_ion << " " << sta1->get_name(ivg::staname::ivs_name) << " " << sta2->get_name(ivg::staname::ivs_name) << " " << _use_me << endl;
    }
    //   std::cout << _scan->get_epoch().get_double_doy() << sta1->get_name(ivg::staname::ivs_name) <<" " << sta2->get_name(ivg::staname::ivs_name) <<" " <<_scan->_source->get_name(ivg::srcname::iers) <<" "<< _oc*1e9 << " " << _group_delay*1e9 << " " << cc*1e9 << " " << _delta_ion*1e9 << " " << tau*1e9 << endl;
    
//    // TEMOPORARY OUTPUT: GEOCENTER DELAY TEST
//    stringstream tmp_out2;
//    tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) 
//             << _group_delay << ", " << _sigma_delay << ", " << _geoc_delay << ", " << _sigma_delay 
//             << ", " << _delay_rate << ", " << _sigma_delay_rate << endl;
//    ivg::Logger::add_to_saveline = tmp_out2.str();      
    
//    // TEMPORARY VASCC2015 OUTPUT (1 COLUMN TAU)
//    stringstream tmp_out2;
//    tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << -1.0 * _oc << endl;
//    ivg::Logger::add_to_saveline = tmp_out2.str();    
//    log<DETAIL>("*** [Step 11bof12] O minus C: ") % _oc % "(gd: " % gd % " / tau: " % tau % " )";    

    // ===================
    // build design matrix
    // ===================
    double dt = _scan->get_epoch().get_double_mjd() - _session->_start.get_double_mjd();

    // intermediate values according to MODEST Handbook 1994 p.46 ff
    ivg::Matrix beta = v_earth_ssb*(1.0/ivg::c);
    ivg::Matrix b2 = ( vel2+v_earth_ssb )*(1.0/ivg::c);
    double gam = 1/sqrt(1.0-(beta.transpose()*beta)(0));
    double rho = 1.0+( k.transpose()*b2 )(0);
    ivg::Matrix dij;
    dij.eye(3);

    ivg::Matrix m1 = ( k*( gam/rho*( 1.0-( beta.transpose()*b2 )(0) ) ) );
    ivg::Matrix m2 = beta*gam;
    ivg::Matrix psi = ( m1+m2 )*-1.0;
    ivg::Matrix E = dij+(beta*(gam-1.0)/(beta.transpose()*beta)(
                             0)-b2*gam)*beta.transpose();
    ivg::Matrix K = E*psi;
    ivg::Matrix B = K.transpose()*trf2crf;
    ivg::Matrix T = (dij-k*b2.transpose()*1.0/rho);
    double y = (-gam*( 1.0-( b2.transpose()*beta )(0) ) );
    ivg::Matrix t = ( E*b_gcrs*(1.0/rho) );
    ivg::Matrix M = T*y*t;

    // insert values in design matrix
    // design_iter is an iterator to the first element of the row of A, thus,
    // the offset of any parameter to this iterator is its index within the
    // parameter list

    // station 1
    // =========
    // coordinates  [m]
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::stax,
            _station1)) = -B(0)*1e-2;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::stay,
            _station1)) = -B(1)*1e-2;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::staz,
            _station1)) = -B(2)*1e-2;
    // clock [0.1 s]
    //*(design_iter+_session->_param_list.get_index(ivg::paramtype::clo, _station1)) = -1.0 * 1e-5;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::clo,
            _station1)) = -1.0e-1;
    // ZWD [0.01 s]
    //*(design_iter+_session->_param_list.get_index(ivg::paramtype::zwd, _station1)) = -mf_hydr1 / ivg::c;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::zwd,
            _station1)) = -mf_wet1*1e-1;
    // gradients [0.01 s]
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::ngr,
            _station1)) = -1.0/( tan( azel1(1) )*sin( azel1(1) )+0.0032 )
                          *cos(azel1(0) )*1e-2;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::egr,
            _station1)) = -1.0/( tan( azel1(1) )*sin( azel1(1) )+0.0032 )
                          *sin(azel1(0) )*1e-2;

    // station 2
    // =========
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::stax,
            _station2)) = B(0)*1e-2;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::stay,
            _station2)) = B(1)*1e-2;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::staz,
            _station2)) = B(2)*1e-2;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::clo,
            _station2)) = 1.0e-1;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::zwd,
            _station2)) =  mf_wet2*1e-1;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::ngr,
            _station2)) = 1.0/( tan( azel2(1) )*sin( azel2(1) )+0.0032 )
                          *cos( azel2(0) )*1e-2;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::egr,
            _station2)) = 1.0/( tan( azel2(1) )*sin( azel2(1) )+0.0032 )
                          *sin( azel2(0) )*1e-2;

    if (par_aplo)
      {
	*par_aplo=((dx2_ntaplo(0)-dx1_ntaplo(0))*B(0)+(dx2_ntaplo(1)-dx1_ntaplo(1))*B(1)+(dx2_ntaplo(2)-dx1_ntaplo(2))*B(2))/ivg::c;

      }

    // EOP
    // ===
    // partial derivative wrt PM [rad]
    //*(design_iter+_session->_param_list.get_index(ivg::paramtype::xpo, "EOP"))  = ( ( K.transpose()*( _scan->_partials_t2c.pmx*b_trs ) )(0)/ivg::c ) * ivg::mas2rad;
    //*(design_iter+_session->_param_list.get_index(ivg::paramtype::ypo, "EOP")) = ( ( K.transpose()*( _scan->_partials_t2c.pmy*b_trs ) )(0)/ivg::c ) * ivg::mas2rad;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::xpo,"EOP"))  
       = ( ( K.transpose()*( _scan->_partials_t2c.pmx*b_trs ) )(0)/ivg::c );

    *(design_iter+_session->_param_list.get_index(ivg::paramtype::ypo,"EOP")) 
       = ( ( K.transpose()*( _scan->_partials_t2c.pmy*b_trs ) )(0)/ivg::c );

    // partial derivative wrt UT1 [rad]
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::ut1, "EOP")) 
       = ( ( K.transpose()*( _scan->_partials_t2c.ut1*b_trs ) )(0)/ivg::c );
    //*(design_iter+_session->_param_list.get_index(ivg::paramtype::ut1,
    //        "EOP")) = ( ( K.transpose()*( _scan->_partials_t2c.ut1*b_trs ) )(
    //                        0)/ivg::c*M_PI/43200.0 );

    // nutation  X/Y [rad]]
    //*(design_iter+_session->_param_list.get_index(ivg::paramtype::nutx, "EOP")) = ( ( K.transpose()*( _scan->_partials_t2c.pnx*b_trs ) )(0)/ivg::c ) * ivg::mas2rad;
    //*(design_iter+_session->_param_list.get_index(ivg::paramtype::nuty, "EOP")) = ( ( K.transpose()*( _scan->_partials_t2c.pny*b_trs ) )(0)/ivg::c ) * ivg::mas2rad;
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::nutx,"EOP")) 
       = ( ( K.transpose()*( _scan->_partials_t2c.pnx*b_trs ) )(0)/ivg::c );

    *(design_iter+_session->_param_list.get_index(ivg::paramtype::nuty,"EOP")) 
       = ( ( K.transpose()*( _scan->_partials_t2c.pny*b_trs ) )(0)/ivg::c );

    // sources [rad]
    ivg::Matrix deriv = _scan->_source->get_unit_vector_ssb_partials_ra();
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::ra,
            _scan->_source->get_name(ivg::srcname::iers)))= ( (deriv.transpose() * M)(0) ) / ivg::c;
    deriv = _scan->_source->get_unit_vector_ssb_partials_dec();
    *(design_iter+_session->_param_list.get_index(ivg::paramtype::dec,
            _scan->_source->get_name(ivg::srcname::iers)))= ( (deriv.transpose() * M)(0) ) / ivg::c;
   
    // baseline-dependent parameters
    // =============================
    // baseline-dependent clock offsets [s]
    // in case of wrong baselinename, turn the name and try to find it again! (e.g. Sv-Yj to Yj-Sv)
    int bl_pos = _session->_param_list.get_index(ivg::paramtype::blcl,sta1->get_name(lettercode)+"-"+sta2->get_name(lettercode));
    if( bl_pos == -1 )
        bl_pos = _session->_param_list.get_index(ivg::paramtype::blcl,sta2->get_name(lettercode)+"-"+sta1->get_name(lettercode));
        
    *(design_iter+bl_pos) = 1.0e-1;
    
    // =============================================================================
    // build apriori vector (only AT/clocks/gradients/EOP, no sources or stations)
    // =============================================================================
    // station 1
    // =========
    // clock [s]
//    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::clo,_station1)) = _scan->_data.at( _sta1_scan_idx ).clo0;
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::clo,_station1)) = 0.0;
    
    // ZWD [s]
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::zwd,_station1)) 
            = (zhd1+zwd1) / ivg::param_unit_fac.at( ivg::paramtype::zwd );
    
    // gradients [mm]
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::ngr,_station1)) 
            = _scan->_data.at( _sta1_scan_idx ).tropo.get_north_gradient()  / ivg::param_unit_fac.at( ivg::paramtype::ngr );
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::egr,_station1)) 
            = _scan->_data.at( _sta1_scan_idx ).tropo.get_east_gradient() / ivg::param_unit_fac.at( ivg::paramtype::egr );

    // station coordinates
    ivg::Matrix sta1xyz =sta1->calc_xyz(_epoch,{"PSD","SEASONALS"});
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::stax,_station1))=sta1xyz(0)/ ivg::param_unit_fac.at( ivg::paramtype::stax );
      *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::stay,_station1))=sta1xyz(1)/ ivg::param_unit_fac.at( ivg::paramtype::stay );
      *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::staz,_station1))=sta1xyz(2)/ ivg::param_unit_fac.at( ivg::paramtype::staz );
     // station 2
    // =========
//    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::clo,_station2)) = _scan->_data.at( _sta2_scan_idx ).clo0;
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::clo,_station2)) = 0.0;
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::zwd,_station2)) 
            = (zhd2+zwd2) / ivg::param_unit_fac.at( ivg::paramtype::zwd );
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::ngr,_station2)) 
            = _scan->_data.at( _sta2_scan_idx ).tropo.get_north_gradient() / ivg::param_unit_fac.at( ivg::paramtype::ngr );
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::egr,_station2)) 
            = _scan->_data.at( _sta2_scan_idx ).tropo.get_east_gradient() / ivg::param_unit_fac.at( ivg::paramtype::egr ); 
    // station coordinates
    ivg::Matrix sta2xyz =sta2->calc_xyz(_epoch,{"PSD","SEASONALS"});
     *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::stax,_station2))=sta2xyz(0)/ ivg::param_unit_fac.at( ivg::paramtype::stax );
     *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::stay,_station2))=sta2xyz(1)/ ivg::param_unit_fac.at( ivg::paramtype::stay );
     *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::staz,_station2))=sta2xyz(2)/ ivg::param_unit_fac.at( ivg::paramtype::staz );
    // EOP apriori setting
    //*(apriori_iter+_session->_param_list.get_index(ivg::paramtype::xpo,"EOP"))
    //        = _session->_param_list.get_param(_session->_param_list.get_index(ivg::paramtype::xpo,"EOP"))->get_apriori();
    //*(apriori_iter+_session->_param_list.get_index(ivg::paramtype::ypo,"EOP"))
    //        = _session->_param_list.get_param(_session->_param_list.get_index(ivg::paramtype::ypo,"EOP"))->get_apriori();
    //*(apriori_iter+_session->_param_list.get_index(ivg::paramtype::ut1,"EOP"))
    //        = _session->_param_list.get_param(_session->_param_list.get_index(ivg::paramtype::ut1,"EOP"))->get_apriori();

    ivg::Eop_series eoptmp=_session->get_eops();
    ivg::Matrix tmperp=eoptmp.calc_erp(_scan->get_epoch());
    ivg::Matrix tmpnut=eoptmp.calc_nut(_scan->get_epoch());
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::xpo,"EOP"))=tmperp(0);
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::ypo,"EOP"))=tmperp(1);
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::ut1,"EOP"))=tmperp(2);
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::nutx,"EOP"))=tmpnut(0);
    *(apriori_iter+_session->_param_list.get_index(ivg::paramtype::nuty,"EOP"))=tmpnut(1);
    
}

// ...........................................................................
double Obs::calc_total_delay(ivg::delaymodel model, ivg::Matrix &v_earth_ssb, ivg::Matrix &vel1, ivg::Matrix &vel2, ivg::Matrix &b_gcrs, ivg::Matrix &b_trs, ivg::Matrix &k)
// ...........................................................................
{
    tictoc time;
    time.tic();
    const double c = ivg::c;
    const double c2 = pow(ivg::c, 2.0);
    const double c3 = pow(ivg::c, 3.0);
    const double LC = 1.48082686741e-8;
        
    double tau = 0.0; // FJ: Initialize return value.
    
    // init some variables and do some calculations needed for all delay-models
    ivg::Matrix trf2crf = _scan->_trf2crf;

    // JPL ephemeris
    double r[6];
    int err_code = -1;
    
    std::map<string, int> jpl_args = {{"Mercury", 1},{"Venus", 2},{"Earth", 3}, {"Mars", 4},{"Jupiter", 5},
                                      {"Saturn", 6},{"Uranus", 7},{"Neptune", 8},{"Moon", 10},{"Sun", 11},{"SSB", 12}};

    double au = jpl_get_double( _session->_ephem, JPL_EPHEM_AU_IN_KM )*1e3;
    
    // position and velocity of the Earth w.r.t. SSB
    ivg::Matrix x_earth_ssb(3, 1);
    err_code = jpl_pleph( _session->_ephem, _epoch.get_jd_tdb(), jpl_args["Earth"], jpl_args["SSB"], r, 1 );
    if( err_code != 0)
        throw runtime_error("double Obs::calc_total_delay( ... ): jpl_pleph error code " + std::to_string(err_code));
    for (int i = 0; i < 3; ++i)
    {
        // save in matrix and convert from au and au/d to m and m/s
        x_earth_ssb(i) = r[i] * au;
        v_earth_ssb(i) = r[i+3] * au/86400.0;
    }
    // (1) estimate the barycentric station vectors (Eq 11.6)

    // station 1 w.r.t. SSB
    ivg::Analysis_station * sta1 = _scan->_data.at( _sta1_scan_idx ).sta_ptr ;
   
    ivg::Matrix x1_trs  = sta1->calc_xyz( _epoch, _session->_sta_geophys_effects,
					  _session->_ephem, &(_session->_eops), trf2crf );
    ivg::Matrix x1_gcrs = trf2crf*x1_trs;
   
    ivg::Matrix x1_ssb = x_earth_ssb + x1_gcrs;
   
    vel1( 0 ) = -iers::omega*x1_trs(1);
    vel1( 1 ) =  iers::omega*x1_trs(0);  // velocity in TRS
    vel1 = trf2crf*vel1;                 // velocity in GCRS
   
    // station 2 w.r.t. SSB
    ivg::Analysis_station * sta2 = _scan->_data.at( _sta2_scan_idx ).sta_ptr;
    ivg::Matrix x2_trs  = sta2->calc_xyz( _epoch, _session->_sta_geophys_effects,
                                                 _session->_ephem, &(_session->_eops), trf2crf  );
        ivg::Matrix x2_gcrs = trf2crf*x2_trs;
    ivg::Matrix x2_ssb = x_earth_ssb + x2_gcrs;
    
    vel2( 0 ) = -iers::omega*x2_trs(1);
    vel2( 1 ) =  iers::omega*x2_trs(0);  // velocity in TRS
    vel2 = trf2crf*vel2;                 // velocity in GCRS
    
    _station1 = sta1->get_name( ivg::staname::ivs_name );
    _station2 = sta2->get_name( ivg::staname::ivs_name );
        
    // baseline in GSRC
    b_gcrs = x2_gcrs-x1_gcrs;

    // baseline in TRS
    b_trs = x2_trs-x1_trs;
   
    if( model == ivg::delaymodel::consensus)
    {
            // =============================================================================
            // compute delay following IERS Conv. 2010, ch 11, item (d)
            //
            // Assuming that the reference time is the UTC arrival time of the VLBI signal
            // at receiver 1, and that it is transformed to the appropriate timescale to be
            // used to compute each element of the geometric model, the following steps are
            // recommended to compute the VLBI time delay.
            // =============================================================================

        
            // unit source vector barycenter-source

      //   Setting *setup=_session->_setup;
      // Setting &ga=(*setup)["galactic_abberation"];
      //if(*ga["apply"])
      //     {
      if ((bool)(*_session->_setup)["galactic_abberation"]["apply"])
	{
	  Setting * ga;
	  ga=&((*_session->_setup)["galactic_abberation"]);
	  k = _scan->_source->get_unit_vector_ssb_ga(_scan->get_epoch().get_double_mjd(),ga);
	}
	    //     }
	    else
            k = _scan->_source->get_unit_vector_ssb();
	   
           
           // (2) estimate the vectors from the Sun, the Moon, and each planet except the Earth to receiver 1 & 2. (Eq 11.3-11.5)
           // (3) estimate the differential gravitational delay for each of those bodies.(Eq 11.1)
           // (5) sum to find the total differential gravitational delay. (eq 11.7)

           double grav_delay = 0.0;
           ivg::Matrix x_body_ssb( 3,1 );
           ivg::Matrix v_body_ssb( 3,1 );
	   
           for ( std::map<string, int>::iterator iter = jpl_args.begin(); iter != jpl_args.end(); ++iter )
           {
               if( iter->first != "Earth" && iter->first != "SSB" )
               {
                   // step (2)
                   // ========
                   err_code = jpl_pleph( _session->_ephem, _epoch.get_jd_tdb(), iter->second, jpl_args["SSB"], r, 1 );
                   if( err_code != 0)
                       throw runtime_error("double Obs::calc_total_delay( ... ): jpl_pleph error code " + std::to_string(err_code));
                   for (int i = 0; i < 3; ++i)
                   {
                       x_body_ssb(i) = r[i] * au;
                       v_body_ssb(i) = r[i+3] * au/86400.0;
                   }

       //            // Eq 11.3
       //            double tim = _epoch.get_double_mjd()-( ( k.transpose()*(x_body_ssb-x1_ssb) )/ivg::c )(0);
       //            ivg::Date ref_time( min( tim,_epoch.get_double_mjd() ) );
       //
       //            // determine body position for new epoch
       //            err_code = jpl_pleph( ephem, ref_time.get_jd_tt(), iter->second, jpl_cent_ssb, r, 1 );
       //            for (int i = 0; i < 3; ++i)
       //            {
       //               x_body_ssb(i) = r[i] * au;
       //               v_body_ssb(i) = r[i+3] * au/86400.0;
       //            }

                   // time of closest approach (from VieVS grav_delay.m)
                   double dt = ( (x_body_ssb-x1_ssb).norm() )(0)/ivg::c;
                   ivg::Matrix xbody = x_body_ssb-v_body_ssb*dt;
                   dt    = ( (xbody-x1_ssb).norm() )(0)/ivg::c;
                   xbody = xbody-v_body_ssb*dt;
                   dt    = ( (xbody-x1_ssb).norm() )(0)/ivg::c;
                   xbody = xbody-v_body_ssb*dt;
                   x_body_ssb = xbody;

                   // vector from planet to receiver 1 (Eq 11.4)
                   ivg::Matrix r1_body = x1_ssb-x_body_ssb;

                   ivg::Matrix r2_body = x2_ssb-x_body_ssb - (v_earth_ssb*(1.0/ivg::c)*
                                         (k.transpose()*b_gcrs)(0));

                   // step (3)
                   double dtau = 2.0*iers::jpl_gm[iter->first]/pow( ivg::c,3 )*
                                 log( (r1_body.norm()+k.transpose()*r1_body)(0)/
                                      (r2_body.norm()+k.transpose()*r2_body)(0) );

                   // step (5)
                   grav_delay += dtau;

                   // higher order relativistic time delay (Eq 11.14)
                   if( iter->first == "Sun ")
                   {
                       // double dgrav = ...
                       // grav_delay += dgrav;
                   }
               }
	        
           }
	   
           // (4) find the differential gravitational delay due to the Earth. (Eq 11.2)
           double dtau = 2.0*iers::GM/pow( ivg::c,3 )*
                         log( ( x1_gcrs.norm()+k.transpose()*x1_gcrs )(0)/
                              ( x2_gcrs.norm()+k.transpose()*x2_gcrs )(0) );

           // (5) sum to find the total differential gravitational delay. (eq 11.7)
           grav_delay += dtau;

           // log<DETAIL>("*** [Step 5of12] Total differential gravitational delay: ") % grav_delay;

           // (6) compute the vacuum delay. (Eq 11.9)
           // determine gravitational delay at the geocenter
           ivg::Matrix r_geo_sun( 3,1,0.0 );
           err_code = jpl_pleph( _session->_ephem, _epoch.get_jd_tdb(), jpl_args["Sun"], jpl_args["Earth"], r, 1 );
           if( err_code != 0)
               throw runtime_error("double Obs::calc_total_delay( ... ): jpl_pleph error code " + std::to_string(err_code));
           
           for (int i = 0; i < 3; ++i)
               r_geo_sun(i) = r[i]*au;
           double U = iers::jpl_gm["Sun"]/(r_geo_sun.norm())(0);

           double term1 = (k.transpose()*b_gcrs)(0)/ivg::c*
                          ( 1.0-( 1+iers::ppn_gamma)*U/pow( ivg::c,2 )
                            -(v_earth_ssb.transpose()*v_earth_ssb)(0)/(2*pow( ivg::c,2 ))
                            -(v_earth_ssb.transpose()*vel2)(0)/pow( ivg::c,2 ));
           double term2 = ( (v_earth_ssb.transpose()*b_gcrs)(0) )/pow( ivg::c,2 )*
                          (1+(k.transpose()*v_earth_ssb)(0)/( 2.0*ivg::c ) );
           tau = ( grav_delay-term1-term2 )/
                        ( ( 1.0+(k.transpose()*( v_earth_ssb+vel2) )(0)/ivg::c ) );
	    
           // TEMPORARY LUCIA SATELLITE COMPARISON OUTPUT
           // stringstream tmp_out2;
           // tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << -1.0 * tau << endl;
           // tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << ( ao2-ao1 );
           // tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << ( ao1 );
           // tmp_out2 << " " << setfill(' ') << setw(21)<< right << scientific << setprecision(14) << ( ao2 ) ;
           // ivg::Logger::add_to_saveline = tmp_out2.str();
    
    }
    else if( model == ivg::delaymodel::duev )
    {
#if 0
        /**********************************************************************
         * Implementation of Duev et al. (2012)
         */
        
        const int MAXIT = 10;
        const double precision = 1e-6;
        const double dt = 1e-7;
        
        ivg::Date T0;
        ivg::Date T2;
        double UE;
        ivg::Matrix rdot2gc(3, 1);
        int n;
        
        // Solve light-time equation for T0
        n = 0;
        do {
            double f;
            ivg::Matrix R01(3, 1);
            ivg::Matrix R0T0(3, 1);

            ++n;
            
            // Source position at T0 (SSB)
            if (_scan->_source->get_type() == moon)
            {
                R0T0 = _scan->_source->calc_rover_position(T0);
            } 
            else if (_scan->_source->get_type() == satellite)
            {
                ivg::Matrix X0_T0_trf(3, 1);
                
                jpl_pleph(_session->_ephem, T0.get_jd_tdb(), jpl_args["Earth"], jpl_args["SSB"], r, 0);
                for (int i = 0; i < 3; ++i)
                {
                    x_earth_ssb_T0(i) = r[i]*au;     // m
                }
                
                T0_corr = T0;
                T0_corr.add_secs(T0.get_leap_sec()*leap + Toffset + (C0 + dC0*T0.get_frac_day()));
                X0_T0_trf = _scan->_source->get_vector_trs(ivg::Date(T0_corr));
                X0_T0 = trf2crf*X0_T0_trf + x_earth_ssb_T0;
                
                T0_corr = T0;
                T0_corr.add_secs(dt*86400.);
                jpl_pleph(_session->_ephem, T0_corr.get_jd_tdb(), jpl_args["Earth"], jpl_args["SSB"], r, 0);
                for (int i = 0; i < 3; ++i)
                {
                    x_earth_ssb_T0(i) = r[i]*au;     // m
                }
                
                T0_corr = T0;
                T0_corr.add_secs(T0.get_leap_sec()*leap + Toffset + (C0 + dC0*T0.get_frac_day()) + dt*86400.);
                X0_T0_trf = _scan->_source->get_vector_trs(ivg::Date(T0_corr));
                X0_T0dt = trf2crf*X0_T0_trf + x_earth_ssb_T0;
            }

            
            R01 = x1_ssb - R0T0;
            
            
            
        } while (fabs(f) < precision && n <= MAXITER)
        
        tau = (
                (T2 - T1)/(1.0 - LC)*(
                    1.0 - 1.0/c2*(pow(v_earth_ssb.norm()(0), 2.0)/2.0 + UE)
                ) - v_earth_ssb.transpose()*b_gcrs/c2
              )*1.0/(
                    1.0 + v_earth_ssb.transpose()*rdot2gc/c2
              );
        
        /* End of Duev et al. (2012)
         */
#endif
    }
    else if( model == ivg::delaymodel::sekido )
    {
        /*********************************************************************
         * Implementation of Sekido & Fukushima (2006)
         * Author: Frederic Jaron
         * Last modified: 08/04/2016
         * Comments:
         *   - The sign of tau is opposite to the one of Eq. (20) in the paper.
         *   - Tested with the observations of GPS satellites (August 2015)
         *     and compared to the results of Lucia Plank.
         */
        
        const double gamma = 1.0;
        const double precision = 1e-7;
        const double dt = 1./(24.*15.*60.); // dt for numerical differentiation of f (in days)
        const int Niter = 10;   // Maximum number of iterations for Newton
        
        double WE = 0.0;             // Gravitational potential 
        ivg::Matrix K(3, 1);          // Pseudo source vector (TT)
        ivg::Matrix R2_unit(3, 1);    // Unit vector R2 (TBD)
        ivg::Matrix V2(3, 1);         // Coordinate velocity of station 2 (TDB)
        double DeltaTg21 = 0.0;       // Gravitational effect
        double H = 0.0;               // Correction term of second order   
        ivg::Matrix XP(3, 1);         // Position of gravitating body (GCRF)
        ivg::Matrix VP(3, 1);         // Velocity of gravitating body (GCRF)
        ivg::Matrix X0_T0(3, 1);    // Source position at T0 (SSB)
        ivg::Matrix X0_T0dt(3, 1, 0.0);    // Source velocity at T0 (SSB)
        ivg::Matrix R1E_T1(3, 1);   // Vector from geocenter to station 1 (SSB)
        ivg::Matrix R2E_T1(3, 1);   // Vector from geocenter to station 2 (SSB)
        double f = 0.0;
        double fdt = 0.0;
        double df = 0.0;
        double DeltaT0 = 0.0;
        double DeltaTg01 = 0.0;   // Gravitationnal delay
        double dt_close = 0.0;            // Time of closest approad wrt T0
        ivg::Matrix x_earth_ssb_T0(3, 1);
        ivg::Matrix R1_T1;
        ivg::Matrix R2_T1;
        ivg::Date T0(_epoch); 
        ivg::Date T0_corr; 
        int n = 0;
        ivg::Matrix k1(3, 1);
        ivg::Matrix k2(3, 1);

        // Gravitational potential at the geocenter
        WE = 0.0;
        for (std::map<string, int>::iterator iter = jpl_args.begin(); 
                iter != jpl_args.end(); ++iter )
        {
            if (iter->first != "Earth" && iter->first != "SSB")
            {
                // Get GCRS planet vector at T1
                jpl_pleph(_session->_ephem, _epoch.get_jd_tdb(), iter->second, 
                            jpl_args["Earth"], r, 1);
                for (int i = 0; i < 3; ++i)
                {
                    XP(i) = r[i]*au;     // m
                }
                WE += (iers::jpl_gm[iter->first])/XP.norm()(0);
            }
        }
        
        // Coordinate velocity of station 2 (TDB)
        V2 = vel2*(1.0 - (1.0 + gamma)*WE/c2 
                - pow(v_earth_ssb.norm()(0), 2.0)/(2.0*c2)
                - (v_earth_ssb.transpose()*vel2)(0)/c2) + v_earth_ssb*(1.0 - 
                1.0/(2*c2)*(v_earth_ssb.transpose()*vel2)(0));
        
        // Solving the light time equation...
        R1E_T1 = x1_gcrs*(1.0 - WE/c2 - LC) 
                - v_earth_ssb*(v_earth_ssb.transpose()*x1_gcrs)(0)/(2.0*c2);
        R2E_T1 = x2_gcrs*(1.0 - WE/c2 - LC) 
                - v_earth_ssb*(v_earth_ssb.transpose()*x2_gcrs)(0)/(2.0*c2);
        do {
            ++n;
            // Source position at T0 (SSB)
            if (_scan->_source->get_type() == moon)
            {
                X0_T0 = _scan->_source->calc_rover_position(T0);
            } 
            else if (_scan->_source->get_type() == satellite)
            {
                // Get satellite position at T0
                jpl_pleph(_session->_ephem, T0.get_jd_tdb(), 
                            jpl_args["Earth"], jpl_args["SSB"], r, 0);
                for (int i = 0; i < 3; ++i)
                {
                    x_earth_ssb_T0(i) = r[i]*au;     // m
                }
                T0_corr = T0;
                X0_T0 = _scan->_source->get_vector_crs(T0_corr) 
                        + x_earth_ssb_T0;     
                
                // Get satellite position at T0 + dt
                T0_corr.add_secs(dt*86400.);
                jpl_pleph(_session->_ephem, T0_corr.get_jd_tdb(), 
                            jpl_args["Earth"], jpl_args["SSB"], r, 0);
                for (int i = 0; i < 3; ++i)
                {
                    x_earth_ssb_T0(i) = r[i]*au;     // m
                }
                X0_T0dt = _scan->_source->get_vector_crs(T0_corr) 
                            + x_earth_ssb_T0;
            }
            
            // Gravitational delay
            ivg::Matrix j(3, 1);    // Unit vector from source to receiver
            DeltaTg01 = 0.0;
            j = (x1_ssb - X0_T0);
            j = j/j.norm()(0);                
            // Loop over all gravitating bodies except for the Earth
            for (std::map<string, int>::iterator iter = jpl_args.begin(); 
                 iter != jpl_args.end(); ++iter )
            {
                if (iter->first != "Earth" && iter->first != "SSB")
                {
                    // 1. Get position and velocity of gravitating body
                    jpl_pleph(_session->_ephem, T0.get_jd_tdb(), iter->second, 
                                jpl_args["SSB"], r, 1);
                    for (int i = 0; i < 3; ++i)
                    {
                        XP(i) = r[i]*au;
                        VP(i) = r[i + 3]*au/86400.0;
                    }                                  
                
                    // 2. Compute time of closest approach (i.e., time offset)
                    dt_close = ((j*c - VP).transpose()*(X0_T0 - XP))(0)/
                                pow((j*c - VP).norm()(0), 2.0);
                    // dt_close must be within boundaries
                    if (dt_close < 0.0)
                    {
                        dt_close = 0.0;
                    }
                    else if (dt_close > (x1_ssb - X0_T0).norm()(0)/c)
                    {
                        dt_close = (x1_ssb - X0_T0).norm()(0)/c;
                    }
                
                    // 3. Get position of grav. body at tclose
                    // get ssb planet vector at T0 + dt
                    jpl_pleph(_session->_ephem, T0.get_jd_tdb() 
                                + dt_close/86400.0 + 2400000.5,iter->second, 
                                jpl_args["SSB"], r, 1);
                    for (int i = 0; i < 3; ++i)
                    {
                        XP(i) = r[i]*au;
                    }                                  
                
                    // 4. Compute gravitational delay
                    double R1J = (XP - x1_ssb).norm()(0);
                    double R0J = (XP - x_earth_ssb).norm()(0);
                    double R01 = R1E_T1.norm()(0);
                
                    DeltaTg01 += 2.0*iers::jpl_gm[iter->first]/c3*log((R1J + 
                                    R0J + R01)/(R1J + R0J - R01));
            }
        }
            
            ivg::Matrix V = X0_T0 - x_earth_ssb - R1E_T1;
            ivg::Matrix Vdt = X0_T0dt - x_earth_ssb - R1E_T1;

            f = (_epoch.get_jd_tdb() - T0.get_jd_tdb())*86400. 
                    - V.norm()(0)/c - DeltaTg01;
            fdt = (_epoch.get_jd_tdb() - (T0.get_jd_tdb() + dt))*86400. 
                    - Vdt.norm()(0)/c - DeltaTg01;
            df = (fdt - f)/(dt*86400.);
            DeltaT0 = -f/df;
            if (fabs(f) > precision  || fabs(DeltaT0) > precision)
                T0.add_secs(DeltaT0);
        } while (fabs(f) > precision && fabs(DeltaT0) > precision 
                    && n < Niter);
        // ...light time equation solved.
                          
        R1_T1 = X0_T0 - x_earth_ssb - R1E_T1;
        R2_T1 = X0_T0 - x_earth_ssb - R2E_T1;        
        K = (R1_T1 + R2_T1)/(R1_T1.norm()(0) + R2_T1.norm()(0));
        R2_unit = R2_T1/R2_T1.norm()(0); // Unit vector R2 (TBD)
        
        // Calculate cross product V2 x R2_unit. There has to be a more
        // elegant way...
        ivg::Matrix V2xR2_unit(3, 1);
        V2xR2_unit(0) = V2(1)*R2_unit(2)-V2(2)*R2_unit(1);
        V2xR2_unit(1) = V2(2)*R2_unit(0)-V2(0)*R2_unit(2);
        V2xR2_unit(2) = V2(0)*R2_unit(1)-V2(1)*R2_unit(0);
        
        // Second order correction term
        H = 1.0/c2*pow(V2xR2_unit.norm()(0), 2.0)*(K.transpose()*b_gcrs)(0)/
                (2.0*R2_T1.norm()(0));
        
        // Gravitational delay
        for (std::map<string, int>::iterator iter = jpl_args.begin(); 
             iter != jpl_args.end(); ++iter )
        {
            if (iter->first != "Earth" && iter->first != "SSB")
            {
            // 1. Get position and velocity of gravitating body
                jpl_pleph(_session->_ephem, T0.get_jd_tdb(),
                          iter->second, jpl_args["SSB"], r, 1);
                for (int i = 0; i < 3; ++i)
                {
                    XP(i) = r[i]*au;
                }                                  
            double R2J = (XP - x2_ssb).norm()(0);
            double R0J = (XP - X0_T0).norm()(0);
            double R20J = (XP - X0_T0 - x2_ssb).norm()(0);
            double R1J = (XP - x1_ssb).norm()(0);
            double R10J = (XP - X0_T0 - x1_ssb).norm()(0);          
                
            DeltaTg21 += (1.0 + gamma)*iers::jpl_gm[iter->first]/c3*log((R2J 
                            + R0J + R20J)/(R2J + R0J - R20J)*(R1J + R0J 
                            + R10J)/(R1J + R0J - R10J));
         }
        }

        // NOTE: The sign of tau is different from (20) in Sekido        
        tau = -(
                    -(
                        1.0 - 2.0*WE/c2
                        - (
                            pow(v_earth_ssb.norm()(0), 2.0) 
                            + 2.0*(v_earth_ssb.transpose()*vel2)(0) 
                          )/(2.0*c2)
                    )*(K.transpose()*b_gcrs)(0)/c
                    - (v_earth_ssb.transpose()*b_gcrs)(0)/c2*(
                            1.0
                            + (R2_unit.transpose()*V2)(0)/c
                            - (K.transpose()*(v_earth_ssb 
                                + vel2*2.0))(0)/(2.0*c)
                    )
                    + DeltaTg21
                )/(
                    (1.0 + (R2_unit.transpose()*V2)(0)/c)*(1.0 + H)
                );
                
        // Prepare source vectors k1 and k2 for atmospheric corrections
        k1 = (X0_T0 - x1_ssb)/(X0_T0 - x1_ssb).norm()(0);
        k2 = (X0_T0 - x2_ssb)/(X0_T0 - x2_ssb).norm()(0);
        k.resize(9, 1);
        for (int i = 0; i < 3; i++)
        {
            k(i) = k1(i);
            k(i+3) = k2(i);
            k(i+6) = K(i);
        }
        /*
         * End of Sekido & Fukushima (2006)         
         *********************************************************************/
    }
    else if( model == ivg::delaymodel::simple )
        throw runtime_error("double Obs::calc_total_delay( ... ): ivg::delaymodel::simple not yet implemented!");
    
    
    return tau;
}

/*
// ...........................................................................
void calc_delay_finite_distance()
// ...........................................................................
// calculate delay according to Sekido and Fukushima (2006)
{
   // arrival time at station 1 T1 in TDB frame
   double T1 = _epoch.get_jd_tdb();

   
   // position and velocity of the Earth XE in TDB-frame at arrival time at station 1 T1
   double r[6];
   int err_code = -1;
   std::map<string, int> jpl_args =
   {
       {"Mercury", 1},{"Venus", 2},{"Earth", 3}, {"Mars", 4},{"Jupiter", 5},
       {"Saturn", 6},{"Uranus", 7},{"Neptune", 8},{"Moon", 10},{"Sun", 11},
       {"SSB", 12}
   };
   double au = jpl_get_double( _session->_ephem, JPL_EPHEM_AU_IN_KM )*1e3;

   ivg::Matrix::XE_T1;
   ivg::Matrix::VE;
   err_code = jpl_pleph( _session->_ephem, _epoch.get_jd_tdb(),
                         jpl_args["Earth"], jpl_args["SSB"], r, 1 );
   for (int i = 0; i < 3; ++i)
   {
       // save in matrix and convert from au and au/d to m and m/day
       XE_T1(i) = r[i] * au;
       VE = r[i+3] * au/86400.0;
   }
   
   // geocentric station positions in TT frame at T1
   ivg::Analysis_station * sta1 = _scan->_data.at( _sta1_scan_idx ).sta_ptr ;
   ivg::Matrix x1_t1  = sta1->calc_xyz( _epoch, _session->_sta_geophys_effects,
                                        _session->_ephem, &(_session->_eops) );
   x1_t1 = trf2crf*x1_t1;
   ivg::Analysis_station * sta2 = _scan->_data.at( _sta2_scan_idx ).sta_ptr ;
   ivg::Matrix x2_t1  = sta2->calc_xyz( _epoch, _session->_sta_geophys_effects,
                                        _session->_ephem, &(_session->_eops) );
   x2_t1 = trf2crf*x2_t1;
   
   // geocentric station positions in TDB frame at T1 (eq 15)
   ivg::Matrix R1E_T1 = ( 1.0-ivg::wE/pow( ivg::c,2.0 ) )*x1_t1
                       -( VE*x1_t1/( 2.0*pow( ivg::c,2.0 ) ) )*VE;
   ivg::Matrix R2E_T1 = ( 1.0-ivg::wE/pow( ivg::c,2.0 ) )*x2_t1
                       -( VE*x2_t1/( 2.0*pow( ivg::c,2.0 ) ) )*VE;
   
   // solve light time equation by numerical iteration (Newton-Raphson method)
   ivg::Matrix R1_T1;
   ivg::Matrix R2_T1;

   // Moon's velocity
   v_moon( 3,1,0.0 );
   err_code = jpl_pleph( _session->_ephem, _epoch.get_jd_tdb(),
                         jpl_args["Moon"], jpl_args["SSB"], r, 1 );
   for (int i = 3; i < 6; ++i)
       v_moon = r[i]*au/86400.0;

   // departure time T0 is unknown, use arrival time at station 1 as a first 
   // approximation
   ivg::Date send_epoch = _epoch;
   // position of source X0 in TDB frame
   ivg::Matrix X0_T0;

   double dt_tdb = 1e6;
   while( dt_tdb > 1e-12 )
   {
      // position of source X0 in TDB frame at T0 in current iteration
      X0_T0 = _scan->_source->calc_rover_position( _epoch ) ;

      // TDB frame vector between source and stations (eq 13)
      R1_T1 = X0_T0-XE_T1-R1_T1;
      R2_T1 = X0_T0-XE_T1-R2_T1;

      // gravitational effect on ray path (eq 17)
      double dTg 
      ivg::Matrix x_body_ssb( 3,1 );
      ivg::Matrix v_body_ssb( 3,1 );
      for ( std::map<string, int>::iterator iter = jpl_args.begin();
              iter != jpl_args.end(); ++iter )
      {
          if( iter->first != "Earth" && iter->first != "SSB" )
          {
              // step (2)
              // ========
              err_code = jpl_pleph( _session->_ephem, _epoch.get_jd_tdb(), iter->second,
                                    jpl_args["SSB"], r, 1 );
              for (int i = 0; i < 3; ++i)
              {
                  x_body_ssb(i) = r[i] * au;
                  v_body_ssb(i) = r[i+3] * au/86400.0;
              }

              // time of closest approach (from VieVS grav_delay.m)
              double dt = ( (x_body_ssb-x1_ssb).norm() )(0)/ivg::c;
              ivg::Matrix xbody = x_body_ssb-v_body_ssb*dt;
              dt    = ( (xbody-x1_ssb).norm() )(0)/ivg::c;
              xbody = xbody-v_body_ssb*dt;
              dt    = ( (xbody-x1_ssb).norm() )(0)/ivg::c;
              xbody = xbody-v_body_ssb*dt;
              x_body_ssb = xbody;
  
              // vector from planet to receiver 1 (Eq 11.4)
              ivg::Matrix r1_body = X1_T1-x_body_ssb;
              ivg::Matrix r0_body = X0_T0-x_body_ssb;
              ivg::Matrix r01 = r0_body-r1_body;
  
              // step (3)
              // ========
              double dtau = 2.0*iers::jpl_gm[iter->first]/pow( ivg::c,2 )*
                            log( (r1_body.norm()+r0_body.norm()+r01.norm())(0)/
                                 ( r1_body.norm()+r0_body.norm()-r01.norm())(0) );
  
              // step (5)
              // ========
              dTg += dtau;
          }
      }


      // light time equation (eq 16) and its derivative wrt T0
      double f  = T1-((X0_T0-XE_T1-R1E_T1).norm())(0)/ivg::c-dTg;
      double df = 1.0-(v_moon.norm())(0)/ivg::c-dTg;
      dt_tdb = -f/df;
      cerr << "iterate light time euation: " << dt_tdb << endl;

      // update send epoch (scale difference between TDB and TT is ignored)
      send_epoch.add_secs( dt_tdb*86400.0 );
   }
      
   // pseudo source vector (eq 12)
   ivg::Matrix k = ( R1_T1+R2_T2 )/( R1_T1.norm()(0)+R2_T1.norm(0) );
   


   // delay in TT frame (eq 20)
   double tau = 
}
*/

// ...........................................................................
double Obs::get_obs_variance( bool ion, bool phase,bool sb)
// ...........................................................................
{
  double variance;
  if (phase)
    variance  = pow( _sigma_phase_delay,2.0 );
  else if (sb)
    variance  = pow( _sigma_sb_delay,2.0 );
  else
    variance  = pow( _sigma_delay,2.0 );
  if( ion )
        variance += pow( _sigma_delta_ion,2.0 );

  return variance;
}


// ...........................................................................
double Obs::get_obs_std( bool ion )
// ...........................................................................
{
    return sqrt( get_obs_variance( ion ) );
}

// ...........................................................................
void Obs::get_station_names( string & sta1, string & sta2, ivg::staname staname ) const
// ...........................................................................
{
    sta1 = _scan->_data.at( _sta1_scan_idx ).sta_ptr->get_name( staname );
    sta2 = _scan->_data.at( _sta2_scan_idx ).sta_ptr->get_name( staname );
}

// ...........................................................................
void Obs::get_source_name( string & src ) const
// ...........................................................................
{
    src = _scan->_source->get_name( ivg::srcname::ivs );
}

// ...........................................................................
ivg::delaymodel Obs::get_delaymodel( string model )
// ...........................................................................
{
    if( model == "duev_2012" )
        return ivg::delaymodel::duev;
    else if( model == "sekido_2006" )
        return ivg::delaymodel::sekido;
    else if( model == "sekido_simple_2006" )
        return ivg::delaymodel::simple;
    else if( model == "consensus" )
        return ivg::delaymodel::consensus;
    else
        throw runtime_error("ivg::delaymodel Obs::get_delaymodel( string model ): Unknown delay-model selected.");
    
}

}
