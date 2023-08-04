#include "source.h"

namespace ivg
{
// ...........................................................................
Source::Source()
// ...........................................................................
{
    _vcs=false;
    _defining=false;
    _found=false;
    _special_handling=false;
    
    _ra0 = 0.0;
    _dec0 = 0.0;
    
    _tle = tle_t();
    _sp3 = ivg::Matrix();
    // default-src type is a regular source
    _type = ivg::srctype::source;
    
    _ephem = NULL;
    
    _ra_series.push_back(_ra0);
    _dec_series.push_back(_dec0);
    _epoch_series.push_back(ivg::Date(1970,1.0));
    
    _n_sessions = 1;
    _n_obs = 0;
    _n_obs_in_sess = 0;
    
    _first = ivg::Date(ivg::fake_mjd);
    _mean = ivg::Date(ivg::fake_mjd);
    _last = ivg::Date(ivg::fake_mjd);
    
    _use_me = true;
}

// ...........................................................................
Source::Source( ivg::srcname type, string value ) : Source()
    // ...........................................................................
{
    // in case of iers name, check if source is a special handling source
    // we need this information for crf_analyzer!
    if(type == ivg::srcname::iers)
        _special_handling = (find(special_handlings.begin(), special_handlings.end(), remove_spaces_end(value)) != special_handlings.end());
    
    if(type == ivg::srcname::iers)
        _defining = (find(icrf2_definings.begin(), icrf2_definings.end(), remove_spaces_end(value)) != icrf2_definings.end());
    
    if(_defining == true && _special_handling == true)
        throw runtime_error("Source::Source( ivg::srcname type, string value ): Not possible. Source cannot be _defining=true AND _special_handling=true");
    
    set_name(type, value);
}
// ...........................................................................
Source::Source( ivg::srcname type, string value, double ra0, double dec0) : Source( type, value )
// ...........................................................................
{
    _ra0 = ra0;
    _dec0 = dec0;
    
    _ra_series.at(0) = _ra0;
    _dec_series.at(0) = _dec0;
}
// ...........................................................................
Source::Source(ivg::srcname type, string value, double ra0, double dec0, double sig_ra0, double sig_dec0, double corr) : Source( type, value, ra0, dec0 )
    // ...........................................................................
{
        _sig_ra0 = sig_ra0 * cos(_dec0);
        _sig_dec0 = sig_dec0;
        _corr = corr;
}
// ...........................................................................
Source::Source(ivg::srcname type, string value, int h, int m, double s, bool negative,
               int deg, int min, double sec, double sig_ra_mas, double sig_dec_mas, double corr) : Source( type, value, h, m, s, negative, deg, min, sec )
    // ...........................................................................
{
        _sig_ra0 = (sig_ra_mas*ivg::mas2rad) * cos(_dec0);
        _sig_dec0 = sig_dec_mas * ivg::mas2rad;
        _corr = corr;
}
// ...........................................................................
Source::Source(ivg::srcname type, string value, int h, int m, double s, bool negative,
               int deg, int min, double sec) : Source( type, value )
    // ...........................................................................
{
    set_ra0(h, m, s);
    set_dec0(negative, deg, min, sec);
}

// ...........................................................................
void Source::set_ra0(int h, int m, double s)
// ...........................................................................
{
    _ra0 = ((double)h + ((double)m / 60) + (s / 3600)) * ivg::h2rad;
    
    _ra_series.at(0) = _ra0;
}

// ...........................................................................
void Source::set_dec0(bool negative, int deg, int min, double sec)
// ...........................................................................
{
    if( negative )
        _dec0 = (((double)deg) + ((double)min / 60) + (sec / 3600)) * (-1) * ivg::d2rad;
    else
        _dec0 = (((double)deg) + ((double)min / 60) + (sec / 3600)) * ivg::d2rad;
    
    
    _dec_series.at(0) = _dec0;
}
// ...........................................................................
void Source::set_bck_ra0(int h, int m, double s)
// ...........................................................................
{
    _bck_ra0 = ((double)h + ((double)m / 60) + (s / 3600)) * ivg::h2rad;
    
    _ra_series.at(0) = _ra0;
}

// ...........................................................................
void Source::set_bck_dec0(bool negative, int deg, int min, double sec)
// ...........................................................................
{
    if( negative )
        _bck_dec0 = (((double)deg) + ((double)min / 60) + (sec / 3600)) * (-1) * ivg::d2rad;
    else
        _bck_dec0 = (((double)deg) + ((double)min / 60) + (sec / 3600)) * ivg::d2rad;
    
    
    _dec_series.at(0) = _dec0;
}

// ...........................................................................
void Source::add_local_position( double local_ra, double local_dec, ivg::Date epoch)
// ...........................................................................
{
    _ra_series.push_back(local_ra);
    _dec_series.push_back(local_dec);
    _epoch_series.push_back(epoch);
}
// ...........................................................................
void Source::set_name(ivg::srcname type, string value)
// ...........................................................................
{
    if(type == ivg::srcname::iers)
        _special_handling = (find(special_handlings.begin(), special_handlings.end(), remove_spaces_end(value)) != special_handlings.end());
    
    if(type == ivg::srcname::iers)
        _defining = (find(icrf2_definings.begin(), icrf2_definings.end(), remove_spaces_end(value)) != icrf2_definings.end());
    
    if(_defining == true && _special_handling == true)
        throw runtime_error("Source::Source( ivg::srcname type, string value ): Not possible. Source cannot be _defining=true AND _special_handling=true");
    
    value = remove_spaces_end(value);
    if(value.size()>0)
    {
        _names[type] = value;
    }
}

// ...........................................................................
string Source::get_name( srcname type )
// ...........................................................................
{
    if(type == MAXSRC)
    {
        for ( int i = 0; i<MAXSRC; i++ )
        {
            if (!get_name((ivg::srcname) i).empty())
                return(_names[(ivg::srcname) i]);
        }
    }
    else
    {
        return(_names[type]);
    }
}

// ...........................................................................
ivg::Matrix Source::get_vector_crs(ivg::Date epoch)
// ...........................................................................
{    
    // CALL "calc_rover_position()" IN CASE OF MOON
    if( _type == ivg::srctype::moon )
        throw runtime_error("ivg::Matrix Source::get_vector_crs(ivg::Date epoch): Not supported for ivg::srctype::moon");
    else if( _type == ivg::srctype::satellite )
    {
        // check whether TLE-elements or SP3-IGS-final-orbits are loaded
        // only ONE type of information should be loaded!
        if( _tle.epoch == 0.0 && _sp3.rows() == 0 )
            throw runtime_error("ivg::Matrix Source::get_vector_crs(ivg::Date epoch): No orbital information available for satellite " + _names[ivg::srcname::ivs]);
        else if( _tle.epoch != 0.0 && _sp3.rows() > 0)
            throw runtime_error("ivg::Matrix Source::get_vector_crs(ivg::Date epoch): TLE-elements and sp3-orbits loaded. Only one choice valid.");
        else if( _tle.epoch != 0.0 && _sp3.rows() == 0 )
        {
            // do we need a time offset here?!
            double mjd = epoch.get_double_mjd();
            double t_since = ( mjd+2400000.5 -_tle.epoch )*1440.0;
            
            // check whether it is a deep space or near Earth satellite
            double pos[3]; // satellite position vector
            double sat_params[N_SAT_PARAMS];

            bool is_deep = select_ephemeris( &_tle );
            if( is_deep )
            {
                SDP4_init( sat_params, &_tle );    // switch to an SDx 
                SDP4( t_since, &_tle, sat_params, pos, NULL );
            }
            else
            {
                SGP4_init( sat_params, &_tle );    // switch to an SGx
                SGP4( t_since, &_tle, sat_params, pos, NULL );
            }
//             cerr << "POS: " << _names[ivg::srcname::ivs] << "(" << is_deep << ") : " << _tle.norad_number << ".no " << fixed << setprecision(20) << pos[0] << "|" << pos[1] << "|" << pos[2] << " at " << _tle.epoch << endl;

            ivg::Matrix tmp(3,1,0.0);
            tmp(0) = pos[0]*1000.0;
            tmp(1) = pos[1]*1000.0;
            tmp(2) = pos[2]*1000.0;
            return tmp;
        }
        // in case of sp3_crf data
        else if(_tle.epoch == 0.0 && _sp3.rows() > 0 )
        {
            // correct for difference between GPS and UTC
            epoch.add_secs(epoch.get_leap_sec() - 19.0);
    
            double add = 0.5; // fraction of day +/-
            
            // don't use the whole sp3 information but only a window of e.g. 0.25 days +/- 
            int start_idx  = (_sp3.get_col(0) - epoch.add_days(-add).get_double_mjd()).abs().minIdx();
            int end_idx  = (_sp3.get_col(0) - epoch.add_days(add).get_double_mjd()).abs().minIdx();
                        
            // get time and XYZ position matrix
            ivg::Matrix sub_sp3_txyz = _sp3.get_sub(start_idx,0,end_idx,3);
            
            return(sub_sp3_txyz.interpolate( sub_sp3_txyz(":",0), epoch.get_double_mjd(), "neville" ).get_sub(0,1,0,3).transpose()*1000.0);
        }
    }
    else
        throw runtime_error("ivg::Matrix Source::get_vector_crs(ivg::Date epoch): Function only valid for srctype::satellite");

}
// ...........................................................................
ivg::Matrix Source::get_unit_vector_ssb()
// ...........................................................................
{
    if(_type != ivg::srctype::source)
        throw runtime_error("ivg::Matrix Source::get_unit_vector_ssb(ivg::Date epoch): Function only valid for srctype::source");
    
    ivg::Matrix k(3, 1, 0.0);
    k(0) = cos(_dec0) * cos(_ra0);
    k(1) = cos(_dec0) * sin(_ra0);
    k(2) = sin(_dec0);

    // unit source vector barycenter-source
    k = k / (k.norm())(0);

    return(k);
}
// ...........................................................................
ivg::Matrix Source::get_unit_vector_ssb_ga(double mjd, Setting * ga)
// ...........................................................................
{
    if(_type != ivg::srctype::source)
        throw runtime_error("ivg::Matrix Source::get_unit_vector_ssb(ivg::Date epoch): Function only valid for srctype::source");
    
    ivg::Matrix k(3, 1, 0.0);
    k(0) = cos(_dec0) * cos(_ra0);
    k(1) = cos(_dec0) * sin(_ra0);
    k(2) = sin(_dec0);
    ivg::Matrix A(3,1,0.0), dk(3,1,0.0);
    
    double ra_gc,de_gc;
    if (((*ga).exists("galacitic_center_ra"))&&(!((*ga).exists("galactic_center_ra"))))
      {
	ra_gc=(double) (*ga)["galacitic_center_ra"];
	de_gc=(double) (*ga)["galacitic_center_de"];
      } else
      {
	ra_gc=(double) (*ga)["galactic_center_ra"];
	de_gc=(double) (*ga)["galactic_center_de"];
	
      }
    double refep=(double) (*ga)["ref_epoch"];
    double acc_ssb=((double) (*ga)["ssb_acceleration"])*(mjd-refep)*M_PI*1e-6/(3600*180*365.25);

    
    A(0)=acc_ssb*cos(M_PI*ra_gc/180)*cos(M_PI*de_gc/180);
    A(1)=acc_ssb*sin(M_PI*ra_gc/180)*cos(M_PI*de_gc/180);
    A(2)=acc_ssb*sin(M_PI*de_gc/180);

    dk(0)=A(0)*(k(1)*k(1)+k(2)*k(2))-A(1)*k(0)*k(1)-A(2)*k(0)*k(2);
    dk(1)=A(1)*(k(0)*k(0)+k(2)*k(2))-A(0)*k(0)*k(1)-A(2)*k(1)*k(2);
    dk(2)=A(2)*(k(1)*k(1)+k(0)*k(0))-A(0)*k(0)*k(2)-A(1)*k(1)*k(2);

    k=k+dk;
    // unit source vector barycenter-source
    k = k / (k.norm())(0);

    return(k);
}

// ...........................................................................
ivg::Matrix Source::get_unit_vector_ssb_partials_ra()
// ...........................................................................
{
    ivg::Matrix k(3, 1, 0.0);
    // partials of source unit vector wrt right ascension
    k(0) = -cos(_dec0) * sin(_ra0);
    k(1) =  cos(_dec0) * cos(_ra0);

    return(k);
}

// ...........................................................................
ivg::Matrix Source::get_unit_vector_ssb_partials_dec()
// ...........................................................................
{
    ivg::Matrix k(3, 1, 0.0);
    // partials of source unit vector wrt declination
    k(0) = -sin(_dec0) * cos(_ra0);
    k(1) = -sin(_dec0) * sin(_ra0);
    k(2) =  cos(_dec0);

    return(k);
}

// ...........................................................................
bool Source::check_name(string name)
// ...........................................................................
{
    for ( int i = 0; i<MAXSRC; i++ )
    {
        if (get_name((ivg::srcname) i) == name)
            return(true);
    }
    return(false);
}

// ...........................................................................
double Source::calc_greenwich_hour_angle( ivg::Date gst )
// ...........................................................................
{
    double gha = gst.get_greenwich_sidereal_time() - _ra0;
    if( gha > (2.0*M_PI) )
        gha -= (2.0*M_PI);
    else if (gha < 0.0 )
        gha += (2.0*M_PI);
    return gha;
}

// ...........................................................................
ivg::Matrix Source::calc_rover_position(ivg::Date epoch)
// ...........................................................................
{
    if(_ephem == NULL)
        throw runtime_error("ivg::Matrix Source::calc_rover_position(ivg::Date epoch): No _ephem loaded. No calculation of rover-positions possible");
    
    //Position rover on moon
    // double lambda = -19.5088*ivg::d2rad;
    // double phi = 44.1197*ivg::d2rad;
    // double r = 1734721.5;
    // 2015-07-08 update based on http://lroc.sese.asu.edu/posts/637 
    // 44.1214°N, 340.4884°E, -2640 meters elevation with a reference radius of 1737.4 km.
    double lambda = -19.5116*ivg::d2rad;
    double phi = 44.1214*ivg::d2rad;
    double r = 1734760.0;
    
    //Position Rover ME-System (MeanEarth)
    ivg::Matrix pos_ME(vector<double>{r*cos(lambda)*cos(phi), r*sin(lambda)*cos(phi), r*sin(phi)});

    // Rotations angles from ME to PrincipalSystem only valid for DE403
    // Different rotations should be used for different JPL DE's !!!!!!
    // DE403: P = RZ(63.8986˝) RY(79.0768˝) RX(0.1462˝) M
    // DE430: P = RZ(67.573") RY(78.580") RX(0.285") M 
    ivg::Matrix rot_x, rot_y, rot_z;
    rot_x.rot3D_x(0.1462*ivg::as2rad);
    rot_y.rot3D_y(79.0768*ivg::as2rad);
    rot_z.rot3D_z(63.8986*ivg::as2rad);
    
    //generated rotation matrix
    ivg::Matrix ME2PA = rot_z*rot_y*rot_x; 
    
    //rotation of vector from ME to PA system
    ivg::Matrix pos_PA = ME2PA * pos_ME;

    //get libration from JPL ephems
    double rrd[6];
    int err_code = jpl_pleph( _ephem, epoch.get_jd_tt(), 15, 0, rrd, 0 );
    
    ivg::Matrix phi_l, theta_l, psi_l;;
    phi_l.rot3D_z(rrd[0]);
    theta_l.rot3D_x(rrd[1]);
    psi_l.rot3D_z(rrd[2]);
    
    //generate rotation matrix for libration
    ivg::Matrix libration = psi_l*theta_l*phi_l;
    libration.inv();
    
    //transformation rover position into Inertial Moon System (IM)
    ivg::Matrix pos_IM = libration * pos_PA;;
        
    // position of the Moon(10) w.r.t. to the SSB(12)
    double au = jpl_get_double( _ephem, JPL_EPHEM_AU_IN_KM )*1e3;
    err_code = jpl_pleph( _ephem, epoch.get_jd_tt(), 10, 12, rrd, 0 );
        
    // save in ivg::matrix and convert from au to m
    ivg::Matrix pos_MO(3, 1);
    for (int i = 0; i < 3; ++i)
        pos_MO(i) = rrd[i] * au;;
        
    // position rover in SSB
    ivg::Matrix pos_SSB   = pos_MO + pos_IM; 
    
    return pos_SSB;
}
// ...........................................................................
double Source::calc_arclength_to_sun( ivg::Date epoch )
// ...........................................................................
{
    if(_ephem == NULL)
        throw runtime_error("double Source::calc_arclength_to_sun(): No _ephem loaded. No calculation of sun-positions possible for source "+_names[ivg::srcname::ivs]);
    
    // Example how to call the method:
    // _scan->_source->calc_arclength_to_sun( _epoch ); 
    
    // get unit source vector
    ivg::Matrix k_source = get_unit_vector_ssb();
    k_source = k_source / (k_source.norm())(0);
    
    // get unit sun vector w.r.t. ssb
    double au = jpl_get_double( _ephem, JPL_EPHEM_AU_IN_KM )*1e3;
    
    double r[6];
    ivg::Matrix k_sun(3, 1);
    int err_code = jpl_pleph( _ephem, epoch.get_jd_tt(), 11, 3, r, 0 );
    for (int i = 0; i < 3; ++i)
        k_sun(i) = r[i] * au; // save in matrix and convert from au to m
    
    k_sun = k_sun / (k_sun.norm())(0);
    
    // transformation from cartesian to spherical coordinates (both geocentric; radius = 1)
    ivg::Matrix sun = k_sun;
    ivg::Matrix src = k_source;
    
    double theta_src = atan2(src(1),src(0));
    double phi_src = atan2(src(2),(sqrt((src(0)*src(0))+(src(1)*src(1)))));
    
    double theta_sun = atan2(sun(1),sun(0));
    double phi_sun = atan2(sun(2),(sqrt((sun(0)*sun(0))+(sun(1)*sun(1)))));
    
    // orthodrome between quasar and sun on unit ball w.r.t. geocenter
    double arc = acos(sin(phi_sun)*sin(phi_src)+cos(phi_sun)*cos(phi_src)*cos((theta_sun-theta_src)));
    
    return arc;
}


// ...........................................................................
double Source::calc_arclength_to_source( ivg::Source& other )
// ...........................................................................
{
    ivg::Matrix k_source_1 = this->get_unit_vector_ssb();
    k_source_1 = k_source_1 / (k_source_1.norm())(0);
    
    ivg::Matrix k_source_2 = other.get_unit_vector_ssb();
    k_source_2 = k_source_2 / (k_source_2.norm())(0);
    
    // transformation from cartesian to spherical coordinates (both geocentric; radius = 1)
    double theta_1 = atan2(k_source_1(1),k_source_1(0));
    double phi_1 = atan2(k_source_1(2),(sqrt((k_source_1(0)*k_source_1(0))+(k_source_1(1)*k_source_1(1)))));
    
    double theta_2 = atan2(k_source_2(1),k_source_2(0));
    double phi_2 = atan2(k_source_2(2),(sqrt((k_source_2(0)*k_source_2(0))+(k_source_2(1)*k_source_2(1)))));
    
    // orthodrome between quasara on unit ball w.r.t. geocenter
    return acos(sin(phi_1)*sin(phi_2)+cos(phi_1)*cos(phi_2)*cos((theta_2-theta_1)));
}
// ...........................................................................
void Source::get_position(int &h, int &m, double &s, int &deg, int &min, double &sec)
// ...........................................................................
{
    
    double hours = _ra0 * rad2h;
    double hours_floor = floor(hours);
    double min_in_hours = hours - hours_floor; 
    double minutes = min_in_hours * 60;
    double minutes_floor = floor(minutes);
    double seconds = minutes - minutes_floor;
    
    h = (int) hours_floor;
    m = (int) minutes_floor;
    s = seconds*60;
    
    double degrees = _dec0 * rad2d;
    
    bool negative = false;
    if(degrees<0)
    {
        negative = true;
        degrees = abs(degrees);
    }
    
    double degrees_floor = floor(degrees);
    double degreemin_in_hours = degrees - degrees_floor;
    double degreeminutes = degreemin_in_hours * 60;
    double degreeminutes_floor = floor(degreeminutes);
    double secs = degreeminutes- degreeminutes_floor;
    
    if(negative)
        deg = (int) -degrees_floor;
    else
        deg = (int) degrees_floor;
    
    min = (int) degreeminutes_floor;
    sec = secs*60;
    
}
// ...........................................................................
void Source::show()
// ...........................................................................
{
    cerr << "------------------------------------------------" << endl;
    cerr << " Names:" << endl;
    vector<string> cout_names = {"iers", "icrf", "ivs", "ngs","jpl"};
    for(map<ivg::srcname, string>::const_iterator iter = _names.begin(),end = _names.end(); iter != end; ++iter )
        cerr << " " << cout_names.at(iter->first) << "->" << iter->second << endl;

    cerr << setprecision(15) << " _ra0/_dec0 [rad]: " << _ra0 << " / " << _dec0 <<endl;
    cerr << "_ra_series: ";
    show_vector(_ra_series);
    cerr << "_dec_series: ";
    show_vector(_dec_series);
    cerr << "_epoch_series: ";
    for(auto &epoch: _epoch_series)
        cerr << epoch.get_date_time("YY:DOY:SSSSS") << "|";
    cerr << endl;
    cerr << " First | Mean | Last: " << _first.get_date_time("YYYY/MO/DD") << " | " << _mean.get_date_time("YYYY/MO/DD") << " | " << _last.get_date_time("YYYY/MO/DD") << endl;
    cerr << " #sessions | #obs: " << _n_sessions << " | " << _n_obs << endl;
    cerr << " VCS-only [true/false]: " << _vcs << endl;
    cerr << " Defining [true/false]: " << _defining << endl;
    cerr << " SpecialHandling [true/false]: " << _special_handling << endl;
    cerr << " SourceType (source=0, moon=1, satellite=2): " << _type << endl;
    cerr << " Unit Source Vector Barycenter:";
    get_unit_vector_ssb().show();
    
    vector<ivg::band> bands = {ivg::band::X, ivg::band::S};
    for(auto &info: _band_flux_info)
    {
        if(info.first == ivg::band::X)
             cerr << "Flux Info X-Band with type " << info.second.type << " ";
        else if(info.first == ivg::band::S)
            cerr << "Flux Info S-Band with type " << info.second.type << " ";
        
        if(info.second.type == "M")
        {
            cerr << "Flux: " << info.second.flux;
            cerr << " | MajAx: " << info.second.major_axis;
            cerr << " | Ratio: " << info.second.ratio;
            cerr << " | PA: " << info.second.pa;
            cerr << " | Off1/Off2: " << info.second.off1 << "/" << info.second.off2 << endl;
        }
        else if(info.second.type == "B")
            info.second.baseline_flux.show();
    }
    cerr << "------------------------------------------------" << endl;
}
// ...........................................................................
void Source::add_band_flux_info(ivg::band band, ivg::flux_info info)
// ...........................................................................
{
    // check if band is already set up
    if(_band_flux_info[band].type == "M" || _band_flux_info[band].type == "S")
    {
        log<WARNING>("!!! More than one fluxinfo for ") % _names[ivg::srcname::iers] % ". Overwriting old fluxinfo for this band and type " % _band_flux_info[band].type;
        _band_flux_info[band] = info;
    }
    else
        _band_flux_info[band] = info;
}
// ...........................................................................
double Source::calc_flux(ivg::band band, double wave, ivg::Matrix bl,
                              ivg::Date *epoch_ptr)
// ...........................................................................
{
    // might be necessary to switch from runtime_error to log<WARNING>
    if(_band_flux_info[band].type.empty() || _band_flux_info.empty() )
        throw runtime_error("double Source::calc_flux(): No _band_flux_info set. No flux-calculation possible for source "+_names[ivg::srcname::ivs]);
    
    double obs_flux = 0.0;
    
    // calculate baseline in uv-plane (has to be improved!!!)
    double u = sqrt(pow(bl(0),2)+pow(bl(1),2))*wave;
    double v = bl(2)*cos(_dec0)*wave;
    double gha = calc_greenwich_hour_angle(*epoch_ptr);
    double un = bl(0)*sin(gha)+bl(1)*cos(gha);
    double vn = bl(2)+sin(_dec0)*(-bl(0)*cos(gha)+bl(1)*sin(gha));
    double pbasen = sqrt(pow(un,2.0)+pow(vn,2.0))/1000.0;
    un *= wave;
    vn *= wave;
    
    // in case of type B, the corresponding row related to the baseline-length in baseline_flux has to be found
    if (_band_flux_info[band].type=="B")
    {
        vector<int> idx = _band_flux_info[band].baseline_flux.find_idx(gt,pbasen);
        if(idx.size() == 0)
            obs_flux = _band_flux_info[band].baseline_flux(_band_flux_info[band].baseline_flux.rows()-1,1);
        else
            obs_flux = _band_flux_info[band].baseline_flux(idx.at(0),1);
    }
    else if (_band_flux_info[band].type=="M")
    {
        double flmax = _band_flux_info[band].flux;
        double size  = _band_flux_info[band].major_axis*(M_PI/(3600.0*180.0*1000.0));
        double ratio = _band_flux_info[band].ratio;
        double pa    = _band_flux_info[band].pa*M_PI/180.0;
        
        double arg1 = pow((vn*cos(pa)+un*sin(pa)),2.0);
        double arg2 = pow((ratio*(un*cos(pa)-vn*sin(pa))),2.0);

        double arg = -((pow(M_PI,2.0))/(4.0*0.6931471))*(arg1+arg2)*size*size;
        double fl = 0.0;
        if (arg<20.0)
            fl = flmax*exp(arg);
        if (fl>0.001)
            obs_flux += fl;
    }
    return obs_flux;
}

} // # namespace ivg
