#include "scan.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Scan::Scan()
// ...........................................................................
{
    _source = NULL;
    _epoch = ivg::Date(ivg::fake_mjd);
    
    _duration = 0.0;
}

// ...........................................................................
Scan::Scan( ivg::Date epoch, ivg::Source* src, const ivg::Matrix t2c, Partials_t2c * deriv_ptr  )  : Scan()
// ...........................................................................
{
    _source = src;
    _epoch  = epoch;

    _trf2crf = t2c;
    _partials_t2c = *(deriv_ptr);
    
}

// ...........................................................................
Scan::Scan( const ivg::Scan &other )
// ...........................................................................
{
    _epoch = other._epoch;
    _source = other._source;

    // observation list
    _observations.resize( other._observations.size() );
    copy( other._observations.begin(), other._observations.end(),
          _observations.begin() );

    // list of station pointer
    _stations.resize( other._stations.size() );
    copy( other._stations.begin(), other._stations.end(), _stations.begin() );

    // troposphere list
    _tropo.resize( other._tropo.size() );
    copy( other._tropo.begin(), other._tropo.end(), _tropo.begin() );

    // list of cable calibration data
    _cable_cal.resize( other._cable_cal.size() );
    copy( other._cable_cal.begin(), other._cable_cal.end(), _cable_cal.begin());

    // rotation matrix from TRF to CRF
    _trf2crf = other._trf2crf;
    _partials_t2c = other._partials_t2c;

   // data struct
   _data = other._data;
   
   for(vector<ivg::Obs>::iterator it=_observations.begin(); it<_observations.end(); ++it)
       it->set_scan(this);
   
    // stores information about cablewrap for each station
   _cable_wrap = other._cable_wrap;
   
   _azel_aiming = other._azel_aiming;
   
   _wrap_zone = other._wrap_zone;
   
   // schedulded duration
   _duration = other._duration;
}
// ...........................................................................
Scan & Scan::operator=( const Scan &other )
// ...........................................................................
{
    _epoch = other._epoch;
    _source = other._source;

    // observation list
    _observations.resize( other._observations.size() );
    copy( other._observations.begin(), other._observations.end(),
          _observations.begin() );

    // list of station pointer
    _stations.resize( other._stations.size() );
    copy( other._stations.begin(), other._stations.end(), _stations.begin() );

    // troposphere list
    _tropo.resize( other._tropo.size() );
    copy( other._tropo.begin(), other._tropo.end(), _tropo.begin() );

    // list of cable calibration data
    _cable_cal.resize( other._cable_cal.size() );
    copy( other._cable_cal.begin(), other._cable_cal.end(), _cable_cal.begin());

    // rotation matrix from TRF to CRF
    _trf2crf = other._trf2crf;
    _partials_t2c = other._partials_t2c;

   // data struct
   _data = other._data;
   
   for(vector<ivg::Obs>::iterator it=_observations.begin(); it<_observations.end(); ++it)
       it->set_scan(this);
   
    // stores information about cablewrap for each station within scan
   _cable_wrap = other._cable_wrap;
   
   // stores information about telescope aiming for each station within scan
   _azel_aiming = other._azel_aiming;
   
   _wrap_zone = other._wrap_zone;
   
   // schedulded duration
   _duration = other._duration;

    return *this;
}

// ===========================================================================
// public methods
// ===========================================================================

// ...........................................................................
void Scan::show()
// ...........................................................................
{
    cout << "-----------------SCAN-----------------------" << endl;
    cout << "Epoch (mjd):  " << setprecision(16) << _epoch.get_double_mjd() << endl;
    cout << "Epoch (text): " << _epoch.get_date_time("YYYY-MON-DD HH:MI:SS") << endl;
    if(_source != NULL)
        cout << "Source: " << _source->get_name( ivg::srcname::ivs ) << endl;
    cout << "Stations: ";
    for( int i=0; i<_stations.size(); ++i )
        cout << _stations.at( i )->get_name( ivg::staname::ivs_name ) << " " ;
    cout << endl;
    cout << "#Obs:   " << _observations.size() << endl;
    cout << "Data-Size " << _data.size() << endl;
    cout << "trf2crf:";
    _trf2crf.show();
    cout << "Scheduled duration: " << _duration << endl;
    cout << "Cable-Wraps:" << endl;
    for(auto &wrap: _cable_wrap)
        cout << wrap.first->get_name(ivg::staname::ivs_name) << " " << wrap.second*(180.0/M_PI)<< endl;
        
    cout << "-----------------SCAN-----------------------" << endl;
}

// ...........................................................................
ivg::Date Scan::get_epoch()
// ...........................................................................
{
    return _epoch;
}

// ...........................................................................
ivg::Source* Scan::get_source()
// ...........................................................................
{
    return _source;
}

// ...........................................................................
void Scan::set_epoch( ivg::Date epo )
// ...........................................................................
{
    _epoch = epo;
}

// ...........................................................................
void Scan::set_source( ivg::Source* src )
// ...........................................................................
{
    _source = src;
}


// ...........................................................................
void Scan::set_trf2crf_matrix( ivg::Matrix t2c, Partials_t2c * deriv_ptr )
// ...........................................................................
{
    _trf2crf = t2c;
    _partials_t2c = *(deriv_ptr);
}


// ...........................................................................
bool Scan::includes_station(ivg::Analysis_station *station )
// ...........................................................................
{
    if( find(_stations.begin(), _stations.end(), station) != _stations.end() )
        return true;
    else
        return false;
}


// ...........................................................................
ivg::Obs* Scan::add_obs( Obs new_obs )
// ...........................................................................
{
    _observations.push_back( new_obs );
    return &_observations.back();
}

// ...........................................................................
int Scan::add_sta( ivg::Analysis_station* new_sta )
// ...........................................................................
{
    if(!includes_station(new_sta))
        _stations.push_back( new_sta );

    bool sta_inside = false;
    int ret_sta;
    for( int i=0; i<_data.size(); ++i )
    {
        if( _data.at( i ).sta_ptr == new_sta )
        {
            sta_inside = true;
            ret_sta = i;
            break;
        }
    }

    if( !sta_inside )
    {
        Sta_epo sta;
        sta.sta_ptr = new_sta;

        _data.push_back( sta );
        ret_sta = _data.size()-1;
    }

    return ret_sta;
}

// ...........................................................................
void Scan::add_tropo( ivg::Troposphere new_tropo )
// ...........................................................................
{
    _tropo.push_back( new_tropo );
}

// ...........................................................................
void Scan::add_scan_meteorology( int idx, double t, double p, double h, 
                                 int code, std::string data_type, 
				 std::string met_data, std::string grdname, std::string grdname3,
				 std::string grad_type, std::string grad_name)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ void Scan::add_scan_meteorology( ... ) "<<endl;
    tictoc tim;
    tim.tic();
#endif  
    
    // in case of default value -999.9 (set in Troposphere constructor following NGS card files)
    if( _data.at( idx ).tropo.get_temperature() == -999.9 )
    {
        ivg::Troposphere tro( _data.at(idx).sta_ptr, _epoch );
        tro.set_gpt2_grdfilename(grdname);
        tro.set_gpt3_grdfilename(grdname3);
	
        // use meteorological data (air pressure) and the modified Saastamoinen
        // to calculate the zenith hydrostatic delay
        if( data_type == "meteo")
        {
            // use in-situ data (from vgosDB or NGS card file)
            if( met_data == "insitu" )
            {      
                // use GPT2 if no meteorological data is available (t = -999.0; p = -999.0, h = -99900.00)
                // (set humidity to default value; not given in GPT2)
                //if( t == -999.0 && p == -999.0 && h == -99900.00 )
	        if( t <= -99.0 || p <= 0.0 || h <= -10 )
	        {
                   double pres, temp, dtemp, e, ah, aw; 
                   tro.use_gpt2( pres, temp, dtemp, e, ah, aw );
                   tro.set_meteorology( temp, pres, 99.9, 0 );
                }
                else
                   tro.set_meteorology( t, p, h, code );
            }

            // use ECMWF model data (from VMF1-files)
            else if( met_data == "ecmwf")
            {
                ivg::Matrix Met = _data.at(idx).sta_ptr->interpolate_ext_met_data( "hydrostatic", _epoch, "cspline" );   
                tro.set_meteorology( Met(7), Met(6), Met(8), 0 );   
            }

            // use GPT2 blind model data (set humidity to default value; not given in GPT2)
            else if( met_data == "gpt2" )
            {      
                double pres, temp, dtemp, e, ah, aw; 
                tro.use_gpt2( pres, temp, dtemp, e, ah, aw );
                tro.set_meteorology( temp, pres, 99.9, 0 );
            }      

            // use in-situ data (from vgosDB or NGS card file) and ECWMF data
            // if either no in-situ meteorological data is available or
            // the difference between in-situ and model data are greater than 10cm
            else if( met_data == "insitu_ecmwf_met" )
            {                    
                ivg::Matrix Met = _data.at(idx).sta_ptr->interpolate_ext_met_data( "hydrostatic", _epoch, "cspline" );          
                
                if( ( t == -999.0 && p == -999.0 && h == -99900.00 ) 
                    || ( abs(t-Met(7)) > 0.1 || abs(p-Met(6)) > 0.1  ) )
                {    
                    tro.set_meteorology( Met(7), Met(6), Met(8), 0 );  
                }
                else
                   tro.set_meteorology( t, p, h, code );
            }        

            _data.at(idx).tropo = tro;

            // get zenith hydrostatic delay
            _data.at(idx).tropo.calc_zhd_davis();  
            _data.at(idx).tropo.set_zwd(0.0); 
        }
        // directly use external zenith hydrostatic delays...
        else if( data_type == "zhd")
        {
            // ... from raytracing
            if( met_data == "raytracing" )
            {
                std::map< ivg::Date,std::map<std::string,ivg::Matrix> > rtmap = 
                        _data.at(idx).sta_ptr->get_raytracing_data();
                
                ivg::Matrix rt_delay = rtmap[ _epoch ][ _source->get_name(ivg::srcname::ivs) ];
                tro.set_wet_mapping_factor( rt_delay(0,0) );
		//double mfh,mfw;
		//tro.get_mapping_function( mfh, mfw );
                // tro.set_zwd( rt_delay(0,2) ); // a priori ZWDs???
                tro.set_zwd(0.0);
                //tro.set_zhd( rt_delay(0,1) );
		
		tro.set_zhd(rt_delay(0,3)/rt_delay(0,0));
                _data.at(idx).tropo = tro;                
            }
            // ... from numerical weather models (e.g., the ECMWF)
            else
            {
                // set zenith hydrostatic delay
                ivg::Matrix Met = _data.at(idx).sta_ptr->get_zhd_data();
                int pos  = (Met.get_col(0) - _epoch.get_double_mjd()).abs().minIdx();
                tro.set_zhd( Met(pos,3) );

                // set meteorological data from insitu measurements (only for 
                // refractivity bending effect used for thermal expansion)
                if( t == -999.0 && p == -999.0 && h == -99900.00 )
                {
                   double pres, temp, dtemp, e, ah, aw; 
                   tro.use_gpt2( pres, temp, dtemp, e, ah, aw );
                   tro.set_meteorology( temp, pres, 99.9, 0 );
                }
                else
                   tro.set_meteorology( t, p, h, code );            

                _data.at(idx).tropo = tro;   
                _data.at(idx).tropo.set_zwd(0.0);
            }
        }

        // calculate azimuth and elevation angle        
        ivg::Matrix AzEl;
        // in case of a regular source
        if(_source->get_type() == ivg::srctype::source)
        {
            ivg::Matrix k = _source->get_unit_vector_ssb();
            AzEl = _data.at( idx ).sta_ptr->calc_az_el( _epoch, k, _trf2crf.transpose() );
        }
        else if(_source->get_type() == ivg::srctype::moon || _source->get_type() == ivg::srctype::satellite )
        // in case of a satellite or moon
        {
            AzEl = _data.at( idx ).sta_ptr->calc_az_el( _epoch, _source->get_vector_crs(_epoch), _trf2crf.transpose() );
        }
        
        // calculate asymmetric delay caused by gradients
        double gdel, grN, grE;
        _data.at( idx ).tropo.calc_gradients( AzEl(1), AzEl(0), gdel, grN, grE, grad_type, grad_name);
    }
    
#if DEBUG_VLBI >=2
    cerr<<"--- void Scan::add_scan_meteorology( ... )" <<" : "<<tim.toc()<<" s "<<endl;
#endif 
}

//// ...........................................................................
//void Scan::add_scan_meteorology( int idx, double t, double p, double h,
//                                 int code )
//// ...........................................................................
//{
//    // in case of default value -999.9
//    if( _data.at( idx ).tropo.get_temperature() == -999.9 )
//    {
//        ivg::Troposphere tro( _data.at(idx).sta_ptr, _epoch );
//
//        // use GPT2 if no meteorological data is available (set humidity to default value; not given in GPT2)
//        if( p == -999.0 )
//        {
//           double pres, temp, dtemp, e, ah, aw; 
//           tro.use_gpt2( pres, temp, dtemp, e, ah, aw );
//           tro.set_meteorology( temp, pres, 99.9, 0 );
//        }
//        else
//        {
//           tro.set_meteorology( t, p, h, code );
//        }
//
//        _data.at(idx).tropo = tro;
//
//        // get zenith hydrostatic delay
//        double zhd = _data.at(idx).tropo.calc_zhd_davis();
//
//        // calculate azimuth and elevation angle
////	    ivg::Matrix AzEl( 2,1,0.0 );
////	    AzEl = _data.at( idx ).sta_ptr->calc_az_el( _epoch, *(_source) );
//        ivg::Matrix k = _source->get_unit_vector_ssb();
//        ivg::Matrix AzEl = _data.at( idx ).sta_ptr->calc_az_el( _epoch, k,
//                           _trf2crf.transpose() );
//
//        // calculate vienna mapping functions 1 using GPT2 coefficients
////	    double mf_hydr, mf_wet;
////	    _data.at( idx ).tropo.calc_vmf1_gpt2( AzEl(1), mf_hydr, mf_wet );
//
//        // calculate asymmetric delay caused by gradients
//        double gdel, grN, grE;
//        _data.at( idx ).tropo.calc_gradients( AzEl(1), AzEl(0), gdel, grN, grE );
//    }
//}


// ...........................................................................
void Scan::add_cable_cal( double cc )
// ...........................................................................
{
    _cable_cal.push_back( cc );
}

// ...........................................................................
double Scan::calc_axis_offset_delay_vievs( int idx, ivg::Matrix abb_src_vec )
// compared code to VIEVS axis_stat.m to detect difference in VASCC2015
// ...........................................................................
{    
    // get latitude, longitude and height of current station    
    ivg::Matrix llh( 3,1,0.0 );
    llh = _data.at( idx ).sta_ptr->calc_lat_lon_h( _epoch );
    double phi =  llh(0);
    
    // get source vector
    ivg::Matrix source_unit_vec( 3,1,0.0 );
    source_unit_vec = abb_src_vec; 
    
    ivg::Matrix azel = _data.at( idx ).sta_ptr->calc_az_el( _epoch,source_unit_vec, _trf2crf.transpose() );
    double zd = M_PI_2 - azel(1);  
    double corz = 3.13e-4/tan( azel(1) );
    if(corz > M_PI_2)
        corz*=-1.0;    
    
    double CAZ = cos(azel(0));
    double SAZ = sin(azel(0));
    double SZ  = sin (zd-corz);
    double CZ  = cos (zd-corz);
    cerr << "VieVS "<< (zd-corz) << " " << zd << " "  << (-1.0*corz) << endl;
    
    double offs = _data.at(idx ).sta_ptr->get_antenna_info().axis_offset_length;
    
    double axkt0;
    if( (_data.at( idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_AZEL" ) == 0 )
    {
        axkt0 = - (offs * SZ)/ivg::c;
    }
    else if( (_data.at( idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_EQUA" ) == 0 )
    {
        double csx = SZ*CAZ*cos(phi) + CZ*sin(phi);
        double sn  = acos(csx);    
        double snx = sin(sn);
        axkt0 = -offs * snx/c;
    }
    else if( (_data.at(idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_XYNO" ) == 0 )
    {
        axkt0 = -(offs * sqrt(1 - pow(SZ*CAZ,2)))/ivg::c; 
    }
    else if( (_data.at( idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_XYEA" ) == 0 )
    {
        axkt0 = -(offs * sqrt(1 - pow(SZ*SAZ,2)))/ivg::c; 
    }
    else if( (_data.at(idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_RICH" ) == 0 )
    {
        double lat_rich = 39.06 * ivg::d2rad;
        double lon_rich = 0.12  * ivg::d2rad;

//        fixed_axis_vec(0,0) = cos( lat_rich )* cos( lon_rich ) ;
//        fixed_axis_vec(1,0) = cos( lat_rich )* sin( lon_rich ) ;
//        fixed_axis_vec(2,0) = sin( lat_rich );
        return 0.0;
    }
    else if( (_data.at(idx ).sta_ptr->get_antenna_info().mounting_type).empty() )
    {
//        log<WARNING>("!!! Mounting Type not available of site ") % _data.at( idx ).sta_ptr->get_name(ivg::staname::ivs_name) % ". Setting axis offset delay to 0.0";
        return 0.0;
    }
    else
    {
        stringstream errormessage;
        errormessage << "void calc_axis_offset_delay( int idx ): "
                     << "invalid antenna  mounting_type: " << _data.at(
                         idx ).sta_ptr->get_antenna_info().mounting_type <<
                     "!= (MO_EQUA, MO_AZEL, MO_XYNO)! Exiting";
        throw runtime_error( errormessage.str() );
    }
    return axkt0;
}

// ...........................................................................
double Scan::calc_axis_offset_delay( int idx, ivg::Matrix apparent_src_vec_trf )
// ...........................................................................
{
    // get latitude, longitude and height of current station
    ivg::Matrix llh( 3,1,0.0 );
    llh = _data.at( idx ).sta_ptr->calc_lat_lon_h( _epoch );

    // select antenna mounting type and set unit source vector
    ivg::Matrix fixed_axis_vec( 3,1,0.0);

    //    Mounting type: MO_AZEL (azimuthal), MO_EQUA (equatorial),
    //                   MO_XYNO (XY north), MO_XYEA (XY east),
    //                   MO_RICH (misplaced equatorial RICHMOND)
    double axis_off_del;
    if( (_data.at( idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_EQUA" ) == 0 )
    {
        if( llh(0) > 0.0 )
        {
            fixed_axis_vec(0,0) = 0;
            fixed_axis_vec(1,0) = 0;
            fixed_axis_vec(2,0) = 1;
        }
        else
        {
            fixed_axis_vec(0,0) = 0;
            fixed_axis_vec(1,0) = 0;
            fixed_axis_vec(2,0) = -1;
        }
    }
    else if( (_data.at( idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_AZEL" ) == 0 )
    {
        fixed_axis_vec(0,0) = cos( llh(0) )* cos( llh(1) ) ;
        fixed_axis_vec(1,0) = cos( llh(0) )* sin( llh(1) ) ;
        fixed_axis_vec(2,0) = sin( llh(0) );
    }
    else if( (_data.at(idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_XYNO" ) == 0 )
    {
        fixed_axis_vec(0,0) = -sin( llh(0) )* cos( llh(1) ) ;
        fixed_axis_vec(1,0) = -sin( llh(0) )* sin( llh(1) ) ;
        fixed_axis_vec(2,0) =  cos( llh(0) );
    }
    else if( (_data.at( idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_XYEA" ) == 0 )
    {
        fixed_axis_vec(0,0) = -sin( llh(1) ) ;
        fixed_axis_vec(1,0) =  cos( llh(1) ) ;
        fixed_axis_vec(2,0) =  0;
    }
    else if( (_data.at(idx ).sta_ptr->get_antenna_info().mounting_type).compare( "MO_RICH" ) == 0 )
    {
        double lat_rich = 39.06 * ivg::d2rad;
        double lon_rich = -0.12  * ivg::d2rad;
	 
        fixed_axis_vec(0,0) = cos( lat_rich )* cos( lon_rich ) ;
        fixed_axis_vec(1,0) = cos( lat_rich )* sin( lon_rich ) ;
        fixed_axis_vec(2,0) = sin( lat_rich );
		ivg::Matrix azel=_data.at(idx ).sta_ptr->calc_az_el( apparent_src_vec_trf );
	double axis_offset_delay = -1.0/ivg::c*_data.at(idx).sta_ptr->get_antenna_info().axis_offset_length
	  *sqrt(1-pow(sin(azel(1))*fixed_axis_vec(2,0)+cos(azel(1))*(cos(azel(0))*fixed_axis_vec(0,0)+sin(azel(0))*fixed_axis_vec(1,0)),2));
	return axis_offset_delay;
    }
    else if( (_data.at(idx ).sta_ptr->get_antenna_info().mounting_type).empty() )
    {
//        log<WARNING>("!!! Mounting Type not available of site ") % _data.at( idx ).sta_ptr->get_name(ivg::staname::ivs_name) % ". Setting axis offset delay to 0.0";
        return 0.0;
    }
    else
    {
        stringstream errormessage;
        errormessage << "void calc_axis_offset_delay( int idx ): "
                     << "invalid antenna  mounting_type: " << _data.at(
                         idx ).sta_ptr->get_antenna_info().mounting_type <<
                     "!= (MO_EQUA, MO_AZEL, MO_XYNO)! Exiting";
        throw runtime_error( errormessage.str() );
    }
    ivg::Matrix azel=_data.at(idx ).sta_ptr->calc_az_el( apparent_src_vec_trf );
   
    double axis_offset_delay = -1.0/ivg::c
                               *_data.at(idx).sta_ptr->get_antenna_info().axis_offset_length
                               *sqrt(1.0-pow((apparent_src_vec_trf.transpose()*fixed_axis_vec)(0),2));
    return axis_offset_delay;
}
double Scan::calc_gravdef_delay( int idx,ivg::Matrix apparent_src_vec_trf )
{
  ivg::Matrix azel = _data.at( idx ).sta_ptr->calc_az_el
            (_epoch, _trf2crf*apparent_src_vec_trf,_trf2crf.transpose() );
  
  double gd_delay=_data.at( idx ).sta_ptr->calc_gravdef(azel(1));
  return gd_delay*1e-12;
  
}
// ...........................................................................
double Scan::calc_thermal_expansion_delay( int idx, ivg::Matrix apparent_src_vec_trf )
// reference: Nothnagel (2009), Conventions on thermal expansion modelling of
//            radio telescopes for geodetic and astrometric VLBI, J Geod (2009),
//            83:787–792, DOI 10.1007/s00190-008-0284-z
// ...........................................................................
{
    double defo = 0.0;

    // get latitude, longitude and height of current station
    ivg::Matrix llh = _data.at( idx ).sta_ptr->calc_lat_lon_h( _epoch );

    // get azimuth and elevation
    //ivg::Matrix AzEl = _data.at( idx ).sta_ptr->calc_az_el( _epoch, *(_source) );
    //ivg::Matrix k = _source->get_unit_vector_ssb();
    ivg::Matrix azel = _data.at( idx ).sta_ptr->calc_az_el
            (_epoch, _trf2crf*apparent_src_vec_trf,_trf2crf.transpose() );
   
    // approximate bending angle (already included in apparent source vector)
//    double zencor = _data.at(idx).tropo.calc_bending_angle(azel(1),llh(2),ivg::bendingmode::INSITU);
    double zencor = 0.0;            // as zencor is now zero, we might want to change intermediate values from zd to elev!?
    double z = M_PI/2.0-azel(1);

    // get declination
    //double dec = _source->get_dec0();
    ivg::Matrix k =  apparent_src_vec_trf;
    k = k / (k.norm())(0);
    double dec = atan2( k(2),sqrt( pow( k(0),2.0 )+pow( k(1),2.0 ) ) );

    // calculate some intermediate values
    double sz = sin( z-zencor );        // sin(zenith)
    double cz = cos( z-zencor );        // cos(zenith)
    double sa = sin( azel(0) );         // sin(azimuth)
    double ca = cos( azel(0) );         // cos(azimuth)
    double cd = cos( dec );             // cos(declination)
    double sd = sin( dec );             // sin(declination)


    // extract coefiicients and length for calculation of thermal expansion
    double temp0    = _data.at( idx ).sta_ptr->get_antenna_info().ref_temp;
    double gamma_f  = _data.at( idx ).sta_ptr->get_antenna_info().found_coeff;
    double hf       = _data.at( idx ).sta_ptr->get_antenna_info().found_height;
    double gamma_a  = _data.at(
                          idx ).sta_ptr->get_antenna_info().fixed_axis_coeff;
    double hp       = _data.at(
                          idx ).sta_ptr->get_antenna_info().fixed_axis_length;
    double hv       = _data.at( idx ).sta_ptr->get_antenna_info().vertex_height;
    double hs       = _data.at( idx ).sta_ptr->get_antenna_info().focus_height;
    double axis_off = _data.at(
                          idx ).sta_ptr->get_antenna_info().axis_offset_length;

    // antenna focus factor depending on focus type (Otoshi and Young 1982))
    double fa;
    if( _data.at( idx ).sta_ptr->get_antenna_info().focus_type == "FO_PRIM" )
    {
        fa = .9;
    }
    else if( _data.at( idx ).sta_ptr->get_antenna_info().focus_type == "FO_SECN" )
    {
        fa = 1.8;
    }
    else if( _data.at( idx ).sta_ptr->get_antenna_info().focus_type.empty() )
    {
//        log<WARNING>("!!! Antenna Focus not available of site ") % _data.at( idx ).sta_ptr->get_name(ivg::staname::ivs_name) % ". Setting thermal expansion delay to 0.0";
        return 0.0;
    }
    else
    {
        stringstream errormessage;
        errormessage << "double Scan::calc_thermal_expansion_delay( int idx, "
                     << "ivg::Matrix abb_src_vec ): invalid antenna  focus: "
                     << _data.at( idx ).sta_ptr->get_antenna_info().focus_type
                     << "!= (FO_PRIM, FO_SECN)! Exiting";
        throw runtime_error( errormessage.str() );
    }

    // get current temperature
    double temp = _data.at( idx ).tropo.get_temperature();


    // actual calculation (eq 2-6 of of Nothnagel (2009) assuming that time lags
    // [t-delta_t] are zero)

    // (1) determine contribution of axis offset which differs for the different
    //     mounts
    double ao_contribution = 0.0;
    // Az-El mount
    if( _data.at( idx ).sta_ptr->get_antenna_info().mounting_type == "MO_AZEL" )
        ao_contribution = axis_off*sz;
    // polar mount
    else if( _data.at( idx ).sta_ptr->get_antenna_info().mounting_type ==
             "MO_EQUA" &&
             axis_off == 0.0 )
        ao_contribution = axis_off*cd;
    // polar mount with axis displaced
    else if( _data.at( idx ).sta_ptr->get_antenna_info().mounting_type ==
             "MO_EQUA" &&
             axis_off != 0.0 )
      ao_contribution = axis_off*cd;
	  //                    *sqrt( 1.0-( cz*sin(llh(0))
          //                             +sz*cos(llh(0))*( ca*cos(llh(1))
          //                                     -sz*sin(llh(1)) ) ) );
    //  X/Y mounts, primary axis north–south
    else if( _data.at( idx ).sta_ptr->get_antenna_info().mounting_type ==
             "MO_XYNO" )
        ao_contribution = axis_off*sqrt( 1-pow( sz*ca,2.0 ) );
    //  X/Y mounts, primary axis east–west
    else if( _data.at( idx ).sta_ptr->get_antenna_info().mounting_type ==
             "MO_XYEA" )
        ao_contribution = axis_off*sqrt( 1-pow( cz*sa,2.0 ) );
    // special case: Richmond, i.e., polar mount at Richmont which was designed
    // for  another location
    else if( (_data.at( idx ).sta_ptr->get_antenna_info().mounting_type) ==
             "MO_RICH" )
    {
        double lat_rich = 39.06 * ivg::d2rad;
        double lon_rich = 0.12  * ivg::d2rad;

        ao_contribution = axis_off*sqrt( 1.0-pow( cz*sin(lat_rich)
						  +sz*cos(lat_rich)*( ca*cos(lon_rich)-sa*sin(lon_rich) ) ,2) );
    }
    else
    {
        stringstream errormessage;
        errormessage << "double Scan::calc_thermal_expansion_delay( int idx ): "
                     << "invalid antenna  mounting_type: "
                     << _data.at( idx ).sta_ptr->get_antenna_info().mounting_type
                     << "!= (MO_EQUA, MO_AZEL, MO_XYNO, MO_XYEA, MO_RICH)! Exiting";
        throw runtime_error( errormessage.str() );
    }

    // calculate thermal expansion
    defo = ( gamma_f*( temp-temp0 )*hf*cz
             +gamma_a*( temp-temp0 )*( hp*cz+ao_contribution+hv-fa*hs ) )/ivg::c;
   
    return defo;
}

// ...........................................................................
double Scan::calc_slew_time(ivg::Source* src, ivg::Eop_series* eops, string &infostr)
// ...........................................................................
{
    //calculates slewtimes for all telescopes within the scan
    double slew_time = 0.0;
    
    // we don't need this epoch
    ivg::Date fake_epoch(ivg::fake_mjd);
    
    // in case of initial scan, there is no source defined and 
    // no slew time needs to be calculated
    if(_source != NULL)
    {
        // just for the generation of some info-output for sked-comparisons (e.g. slew times, duration, etc. pp)
        stringstream ss;
        
        // direction of next source
        ivg::Matrix k_new = src->get_unit_vector_ssb();
        ivg::Matrix k_old = _source->get_unit_vector_ssb();
        
        // up to now, there is no big difference between the EOPs at the beginning
        // of the scan and the end of the scan. Thus the same EOPs can be used for determination of aiming position
        // at the end of the scan. Either way the slew time is rounded to the next integer (ceil).
        // -> iteration can be defined as 1.
        // maybe iteration becomes important for scans with long durations
        
        int iteration = 1;
        for( int iter=0; iter<iteration; iter++)
        {
            // we also need to take into account the slewing time for elevation.
            // after determination of slewing time, new _trf2crf needs to be calculated in order to get proper EOPs
            ivg::Matrix crf2trf;
            if(iter == 0)
                crf2trf = _trf2crf.transpose();
            else
            {
                ivg::Date updated_epoch;
                updated_epoch = _epoch;
                updated_epoch.add_secs(_duration);
                updated_epoch.add_secs(slew_time);
                crf2trf = eops->form_crf2trf(updated_epoch);   
            }
            
            for( int i=0; i<_stations.size(); ++i )// loop over TSUKUBA and WETTZELL
            {   
                // azimut and elevation of station i to source k with EOPs from previous/this scan (_trf2crf)
                ivg::Matrix azel_new = _stations.at(i)->calc_az_el(fake_epoch, k_new, crf2trf ); // [rad]
                ivg::Matrix azel_old = _stations.at(i)->calc_az_el(fake_epoch, k_old, crf2trf ); // [rad]

                double old_wrap = _cable_wrap[_stations.at(i)]; // in case of WETTZELL e.g. 610°
                double new_wrap,d_azi;
                _stations.at(i)->calc_wrap_and_rotangle( old_wrap, azel_new(0), new_wrap, d_azi);

                // e.g. new_wrap 420° , d_azi = 190°
                double d_ele = abs(azel_new(1)-azel_old(1));
                // here we need to calculate how long it takes to turn from old_wrap to new_wrap  --> for azi and ele
                ivg::Antenna antenna_i = _stations.at(i)->set_antenna_info();

                // [sec] = [rad] / [rad/sec] + [sec]
                double slewt_azi = abs(d_azi)/antenna_i.azi_rate + antenna_i.azi_const;
                double slewt_ele = abs(d_ele)/antenna_i.ele_rate + antenna_i.ele_const;

                double station_slew_time = max(slewt_azi, slewt_ele);           

                // information stringstream
                ss << _stations.at(i)->get_name(ivs_name) << " " << src->get_name(ivs) << " OldW: " << old_wrap*(180.0/M_PI) << " NewW: " << new_wrap*(180.0/M_PI) << " dAzi: " << d_azi*(180.0/M_PI) << " dEle: " << d_ele*(180.0/M_PI) << " S: " << station_slew_time << " ||| ";

                // store slower telescope
                if(station_slew_time > slew_time)
                    slew_time = station_slew_time;
            }
        }  
        ss << endl;
        // return infostr by reference
        infostr = ss.str();
    }

    // round to next integer
    return ceil(slew_time);
}
// ...........................................................................
double Scan::mean_group_delay_sigma()
// ...........................................................................
{
    double sum = 0.0;
    for(vector<ivg::Obs>::iterator it = _observations.begin(); it<_observations.end(); ++it)
        sum += it->get_sigma_group_delay();
    
    return (sum/(double)get_nobs());
}

/*
bool Scan::operator==( const Scan test ) const
{
   bool out = false;

   if( _source == test._source &&
      abs( _epoch - test._epoch ) < 60.0/86400.0 )
      out = true;

   return out;
}


int Scan::find_obs( const Obs test ) const
{
   int out = -1;

   vector<Obs*>::const_iterator it = _observations.begin();
   for( it; it < _observations.end(); it++ )
   {
      if( *(*it) == test )
      {
         out = it - _observations.begin();
         it = _observations.end();
      }
   }

   return out;
}
*/


}
