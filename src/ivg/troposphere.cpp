#include "troposphere.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Troposphere::Troposphere()
// ...........................................................................
{
    _epoch.set_decimal_date( 1970.0 );

    _temperature = -999.9;
    _pressure    = -999.9;
    _humidity    = -999.9;
    
    _zwd = -999.9;
    _grE = -999.9;
    _grN = -999.9;
    
    _gpt2_grd_file = "/data/bakkari/gpt2_5.grd";
}

// ...........................................................................
Troposphere::Troposphere( ivg::Analysis_station * station, ivg::Date epoch )
// ...........................................................................
{
    _sta = station;
    _epoch = epoch;

    _temperature = -999.9;
    _pressure    = -999.9;
    _humidity    = -999.9;
    
    _zwd = -999.9;
    _grE = -999.9;
    _grN = -999.9;
    
    _gpt2_grd_file = "/data/bakkari/gpt2_5.grd";
}
// ...........................................................................
Troposphere::~Troposphere()
// ...........................................................................
{

}


// ===========================================================================
// public methods
// ===========================================================================

// getter and setter
// ...........................................................................
double Troposphere::get_temperature() const
// ...........................................................................
{
    return _temperature;
}
// ...........................................................................
double Troposphere::get_pressure() const
// ...........................................................................
{
    return _pressure;
}
// ...........................................................................
double Troposphere::get_humidity() const
// ...........................................................................
{
    return _humidity;
}
// ...........................................................................
double Troposphere::get_zwd() const
// ...........................................................................
{
    return _zwd;
}
// ...........................................................................
double Troposphere::get_zhd() const
// ...........................................................................
{
    return _zhd;
}
// ...........................................................................
double Troposphere::get_ztd() const
// ...........................................................................
{
    return _ztd;
}
// .....................................................................
void Troposphere::get_mapping_function( double &mfh, double &mfw )
// .....................................................................
{
    mfh = _mfh;
    mfw = _mfw;
}
// ...........................................................................
double Troposphere::get_gradient_delay( double &grN, double &grE )
// ...........................................................................
{
    grN = _grN;
    grE = _grE;

    return _grad_del;
}
// ...........................................................................
void Troposphere::set_temperature( double t )
// ...........................................................................
{
    _temperature = t ;
}
// ...........................................................................
void Troposphere::set_pressure( double p )
// ...........................................................................
{
    _pressure = p ;
}
// ...........................................................................
void Troposphere::set_humidity( double h )
// ...........................................................................
{
    _humidity = h ;
}
// ...........................................................................
void Troposphere::set_meteorology( double t, double p, double h, int code )
// ...........................................................................
{
    _temperature = t;
    _pressure = p;
    _humidity = h ;
    _humCode = code;
}
// ...........................................................................
void Troposphere::set_station( ivg::Analysis_station * station )
// ...........................................................................
{
    _sta = station;
}
// ...........................................................................
void Troposphere::set_epoch( ivg::Date epo )
// ...........................................................................
{
    _epoch = epo;
}
// ...........................................................................
void Troposphere::add_zwd( double zwd )
// ...........................................................................
{
    if( _zwd == -999.9 )
        _zwd = zwd;
    else
        _zwd += zwd; 
}
// ...........................................................................
void Troposphere::add_egr( double egr )
// ...........................................................................
{
    if( _grE == -999.9 )
        _grE = egr;
    else
        _grE += zwd; 
}
// ...........................................................................
void Troposphere::add_ngr( double ngr )
// ...........................................................................
{
    if( _grN == -999.9 )
        _grN = ngr;
    else
        _grN += ngr; 
}
// ...........................................................................
void Troposphere::show_met_data()
// ...........................................................................
{
    std::cout << "\n Station: " << _sta->get_name( ivg::staname::ivs_name )
              << " at epoch (mjd): " << _epoch.get_double_mjd()
              << "\n temperature [°C]: " << _temperature
              << "\n pressure [hPa]: "   << _pressure
              << "\n humidity [%]: "     << _humidity
              << endl;
}
// ...........................................................................
double Troposphere::calc_zhd_davis( )
// ...........................................................................
{
    // get station latitude, longitude and height
    ivg::Matrix llh;
    llh = _sta->calc_lat_lon_h( "TIDE FREE" );

    // calc ZHD according to mod. Saastamoinen model, see Davis (1985)
//    _zhd = 0.0022768 * _pressure / ( 1 - 0.00266 * cos( 2*llh(0) ) - 0.00028 *
//                                     ( llh(2) / 1000.0 ) ) ;
    _zhd = calc_zhd_davis(_pressure, llh(0), llh(2));

    return _zhd;
}
// ...........................................................................
double Troposphere::calc_zhd_davis(double pressure, double phi, double h )
// ...........................................................................
{
    // calc ZHD according to mod. Saastamoinen model, see Davis (1985)
    return 0.0022768 * pressure / ( 1 - 0.00266 * cos( 2*phi ) - 0.00028 *
                                     ( h / 1000.0 ) ) ;
}
// ...........................................................................
void Troposphere::calc_ztd()
// ...........................................................................
{
    _ztd = _zhd + _zwd;
}
// ...........................................................................
void Troposphere::calc_vmf1( double el, std::string interpolation_type,
                             double &vmf1h, double &vmf1w )
// ...........................................................................
{

    // zenith distance
    double v = M_PI/2.0 - el;

    // get station latitude, longitude and height
    ivg::Matrix llh;
    llh = _sta->calc_lat_lon_h();

    // get mjd
    double mjd = _epoch.get_double_mjd();

    ivg::Matrix vmf_data = _sta->interpolate_ext_met_data( "mapping", _epoch, interpolation_type );
    double ah = vmf_data(0,1);
    double aw = vmf_data(0,2);

    // calculate vienna mapping function 1 coefficients
    iers::vmf1_( &ah, &aw, &mjd, &llh(0), &v, &vmf1h, &vmf1w );
    
    _mfh = vmf1h;
    _mfw = vmf1w;
}
// ...........................................................................
void Troposphere::calc_vmf3( double el, std::string interpolation_type,
                             double &vmf3h, double &vmf3w )
// ...........................................................................
{

    // zenith distance
    double v = M_PI/2.0 - el;

    // get station latitude, longitude and height
    ivg::Matrix llh;
    llh = _sta->calc_lat_lon_h();

    // get mjd
    double mjd = _epoch.get_double_mjd();

    ivg::Matrix vmf_data = _sta->interpolate_ext_met_data( "mapping3", _epoch, interpolation_type );
    double ah = vmf_data(0,1);
    double aw = vmf_data(0,2);

    // calculate vienna mapping function 1 coefficients
    iers::vmf3_( &ah, &aw, &mjd, &llh(0), &llh(1),&v, &vmf3h, &vmf3w );
    
    _mfh = vmf3h;
    _mfw = vmf3w;
}

// ...........................................................................
void Troposphere::calc_gmf( double el, double &gmfh, double &gmfw )
// ...........................................................................
{

    double v = M_PI/2.0 - el;
    double mjd = _epoch.get_double_mjd();

    // get station latitude, longitude and height
    ivg::Matrix llh;
    llh = _sta->calc_lat_lon_h();

    // calculate global mapping function coefficients
    iers::gmf_( &mjd, &llh(0), &llh(1), &llh(2), &v, &gmfh, &gmfw );

    _mfh = gmfh;
    _mfw = gmfw;
}
// ...........................................................................
void Troposphere::calc_vmf1_gpt2( double el, double &vmf1h, double &vmf1w )
// ...........................................................................
{
    // zenith distance
    double v = M_PI/2.0 - el;

    // get station latitude, longitude and height
    ivg::Matrix llh = _sta->calc_lat_lon_h( "TIDE FREE" );
    double mjd = _epoch.get_double_mjd();
    double ah,aw;
    _sta->get_gpt2_ah_aw(&ah,&aw);
    // use gpt2 to get mapping function coefficients
    if ((ah==0)&&(aw==0))
      {
	int it = 0;
	int nstat = 1;
	double pres, temp, dtemp, e, undu;

	char * gpt2_grd_file = new char[_gpt2_grd_file.size() + 1];
	strcpy(gpt2_grd_file, _gpt2_grd_file.c_str());
	iers::gpt2_( &mjd, &llh(0), &llh(1), &llh(2), &nstat, &it, gpt2_grd_file, &pres, &temp,
		     &dtemp, &e, &ah, &aw, &undu, strlen(gpt2_grd_file) );
	_sta->set_gpt2_ah_aw(ah, aw);
      }
    // calculate vienna mapping functions 1
    iers::vmf1_ht_( &ah, &aw, &mjd, &llh(0), &llh(2), &v, &vmf1h, &vmf1w );

    _mfh = vmf1h;
    _mfw = vmf1w;
}
// ...........................................................................
void Troposphere::calc_vmf3_gpt3( double el, double &vmf3h, double &vmf3w )
// ...........................................................................
{
    // zenith distance
    double v = M_PI/2.0 - el;

    // get station latitude, longitude and height
    ivg::Matrix llh = _sta->calc_lat_lon_h( "TIDE FREE" );
    double mjd = _epoch.get_double_mjd();
    double ah,aw;
    _sta->get_gpt3_ah_aw(&ah,&aw);
    // use gpt2 to get mapping function coefficients
    if ((ah==0)&&(aw==0))
      {
	int it = 0;
	int nstat = 1;
	double pres, temp, dtemp, e, undu,Tm,la,Gn_h,Ge_h,Gn_w,Ge_w;

	char * gpt3_grd_file = new char[_gpt3_grd_file.size() + 1];
	strcpy(gpt3_grd_file, _gpt3_grd_file.c_str());

	iers::gpt3_1_( &mjd, &llh(0), &llh(1), &llh(2), &nstat, &it, gpt3_grd_file, &pres, &temp,
		       &dtemp, &Tm, &e, &ah, &aw, &la,&undu, &Gn_h,&Ge_h, &Gn_w, &Ge_w, strlen(gpt3_grd_file) );
	_sta->set_gpt3_ah_aw(ah, aw);
      }
    // calculate vienna mapping functions 1
    iers::vmf3_ht_( &ah, &aw, &mjd, &llh(0), &llh(1), &llh(2), &v, &vmf3h, &vmf3w );

    _mfh = vmf3h;
    _mfw = vmf3w;
}


// ...........................................................................
void Troposphere::calc_gradient_mf( double el )
// ...........................................................................
{
    // calculate gradient mapping function using the gradient model by
    // Chen and Herring (1997)

    double C = 0.0032;

    _mg = 1.0 / ( tan(el)* sin(el) + C );
}
// ...........................................................................
void Troposphere::calc_gradients( double el, double az,
                                  double &delay, double &grN, double &grE,
				  std::string grad_type, std::string grad_file)
// ...........................................................................
{
    // get station latitude, longitude and height
    ivg::Matrix llh;
    llh = _sta->calc_lat_lon_h();

    // calculate north and east gradients
    if ( grad_type == "none" ) {
      grN=0;
      grE=0;
      delay=0;
    } else if ( grad_type == "apg" ){
      
      iers::apg_( &llh(0), &llh(1), &az, &el, &delay, &grN, &grE );
    } 

    // set delay due to asymmetric gradients
    _grad_del = delay;
    _grN = grN*1e-3;
    _grE = grE*1e-3;
}
// ...........................................................................
double Troposphere::calc_gradient_delay( double el, double az )
// ...........................................................................
{
      _grad_del = ( 1.0 / ( sin(el) * tan(el) + 0.0031 ) 
                  * ( _grN * cos(az) + _grE * sin(az) ) ) ;              
      return _grad_del;            
}
// ...........................................................................
void Troposphere::calc_tropo_delay( std::string mf, double el, double az )
// ...........................................................................
{
    // seclect and calculate hydrostatic and wet mapping function
    if( mf.compare( "vmf1" ) == 0 )
    {
        double vmf1h, vmf1w;
        calc_vmf1_gpt2( el, vmf1h, vmf1w );
    }
    else if( mf.compare( "gmf" ) == 0  )
    {
        double gmfh, gmfw;
        calc_gmf( el, gmfh, gmfw );
    }
    else
    {
        stringstream errormessage;
        errormessage << "void calc_tropo_delay( std::string mf ): "
                     << "No valid mapping function selected!" << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    // calculate gradiant mapping function; according to Chen and Herring (1997)
    calc_gradient_mf( el );

    // calculate gradients
    double delay, grN, grE;
    calc_gradients( el, az, delay, grN, grE );

    // calculate tropospheric delay (in slant direction)
    _std = _mfh * _zhd + _mfw * _zwd + _mg * ( _grN * cos(az) + _grE * sin(az) );
}


// blind models
// .....................................................................
void Troposphere::use_gpt2( double &pres, double &temp, double &dtemp,
                            double &e, double &ah, double &aw )
// .....................................................................
{
    double undu;
    int nstat = 1;
    int it = 0;

    // get station latitude, longitude and height
    ivg::Matrix llh;
    llh = _sta->calc_lat_lon_h( "TIDE FREE" );

    // get mjd
    double mjd = _epoch.get_double_mjd();

    //cout << "DIE SELBE?: " << mjd << endl;
    //tictoc timeCount;
    //timeCount.tic();
    
    char * gpt2_grd_file = new char[_gpt2_grd_file.size() + 1];
    strcpy(gpt2_grd_file, _gpt2_grd_file.c_str());
    iers::gpt2_( &mjd, &llh(0), &llh(1), &llh(2), &nstat, &it, gpt2_grd_file, &pres, &temp,
                 &dtemp, &e, &ah, &aw, &undu, strlen(gpt2_grd_file) );

    // timeCount.cerr_toc("ZEITFRESSER gpt2_");
    //if( fabs( mjd - 51923.0000 )< 1e-6 )
    //   cerr << "GPT " << setprecision(21) <<llh(0) << " " << llh(1) << " " << llh(2) << " " << ah << endl;
    // ah-diff von 7e-8 gesucht
}
// .....................................................................
// Troposphere tie corrections
void Troposphere::calc_tropospheric_ties( double H, double H0, double phi0,
        double p0, double t0, double e0 )
// .....................................................................
{

    // calculate pressure of co-located site
    double p = p0 * pow( ( 1.0 - ivg::gamma * (H-H0) / t0 ),
                         ( ivg::g / (ivg::gamma * ivg::R)) ) ;

    // correction term for ZHD and ZWD
    double dZHD = ( 0.0022768* ( p - p0 ) / ( 1.0 - 0.00266 * cos(
                        2* phi0 ) - 0.28 * H0 * 1.0e-6 ) ) * 1000.0;
    double dZWD = ( -2.789* e0 / pow( t0,
                                      2.0 ) * ( 5383 / t0 - 0.7803 ) * gamma * ( H - H0 ) ) * 1000.0 ;

    cout << "dzhd " << dZHD << ", dzwd " << dZWD << endl;

    // correction term for ZTD
    double dZTD = dZHD + dZWD;
}
// .....................................................................
// Troposphere tie corrections
double Troposphere::calc_tropospheric_ties_zwd( double H, double H0,
        double p0, double t0, double e0 )
// .....................................................................
{

    // calculate pressure of co-located site
    double p = p0 * pow( ( 1.0 - ivg::gamma * (H-H0) / t0 ),
                         ( ivg::g / (ivg::gamma * ivg::R)) ) ;

    // correction term for ZWD
    double dZWD = ( -2.789* e0 / pow( t0,
                                      2.0 ) * ( 5383 / t0 - 0.7803 ) * gamma * ( H - H0 ) ) * 1000.0 ;


    // correction term for ZWD
    return dZWD;
}
// .....................................................................
double Troposphere::calc_bending_angle(double elevation,double sitheight, bendingmode type)
// calculate tropospheric bending angle 
// (according to CALC11 $MK5_ROOT/progs/calc11/caxom.f)
// .....................................................................
{
    if(type==bendingmode::APPROX)
    {
        return 3.13e-4/tan(elevation);
    }
    
    double temperature, pressure, humid;
    
    if(type==bendingmode::INSITU)
    {
        temperature = _temperature+273.15;
        pressure = _pressure/1.33322368;
        humid = _humidity;
    }
    else if(type==bendingmode::HEIGHT)
    {
        // temperature in Kelvin
        temperature = 293.15-6.5e-3*sitheight;
        // pressure in mm of Mercury (Hg)
        pressure = 760.0*pow(1.0-6.5e-3*sitheight/293.15,5.26);
        // rel. humidity 
        humid = 0.5;
    }
    
    // calculate zenith distance in degree
    double zd = 90.0-elevation*ivg::rad2d;
    
    // reference values
    double a1 = 0.40816;
    double a2 = 112.30;
    double b1 = 0.12820;
    double b2 = 142.88;
    double c1 = 0.80000;
    double c2 = 99.344;
    double p1 = 760.0;
    double t1 = 273.0;
    
    double z1 =  91.870;   
    vector<double> w = { 22000.0,17.149,4684.1,38.450 };
    
    // evaluate statement functions
    double d3 = _delta_fun(z1,zd,c1,c2,zd);
    d3 += 1.0;
    double fpx = _delta_fun(p1,pressure,a1,a2,zd);
    double fp = (pressure/p1)*(1.0-fpx/d3);
    double ftx = _delta_fun(t1,temperature,b1,b2,zd);
    double ft = (t1/temperature)*(1.0-ftx/d3);
    double fw = 1.0+(w.at(0)*humid*exp((w.at(1)*temperature-w.at(2))/(temperature-w.at(3)))/(temperature*pressure));
    
    // calculate optical refraction
    vector<double> e = {46.625, 45.375, 4.1572, 1.4468,0.25391,2.2716,-1.3465, 
                        -4.3877, 3.1484,4.5201,-1.8982,0.89000};
    double u = (zd-e.at(0))/e.at(1);
    double x = e.at(10);
    for(int i=1;i<=8;++i)
        x = e.at(10-i)+u*x;
    
    // COMBINE FACTORS AND FINISH OPTICAL FACTOR
    double zencor = ft*fp*fw*(exp(x/d3)-e.at(11));    
    if(zencor > M_PI_2) 
        zencor*=-1.0;
    
    return -1.0*zencor*ivg::as2rad;
}

// statement function for calculation of tropospheric bending
double Troposphere::_delta_fun(double ad1, double ad2, double bd1,double bd2,double zd2) 
{
     return ((ad2-ad1)*exp(bd1*(zd2-bd2)));
}

// .....................................................................
ivg::Matrix Troposphere::build_apparent_source_vector(const ivg::Matrix &src,
                                                      const ivg::Matrix &azel,
                                                      double sitheight,
                                                      bendingmode type)
// .....................................................................
{
    ivg::Matrix app_src = src;
    
    // approximate bending angle
    double  rho_trop = calc_bending_angle(azel(1),sitheight,type);
    
    // construct source vector corrected for trop. bending
    double az = azel(0);
    double zd = M_PI_2-azel(1)-rho_trop;  
    
    // topocentric unit vector
    app_src(0) = cos(zd);
    app_src(1) = sin(zd)*sin(az);
    app_src(2) = sin(zd)*cos(az);
    app_src = app_src/(app_src.norm())(0);
    
    // convert to geocentric system
    return _sta->form_topo2geo()*app_src;
}
// ...........................................................................
void Troposphere::set_gpt2_grdfilename(string fname)
// ...........................................................................
{
    _gpt2_grd_file = fname;
}
// ...........................................................................
void Troposphere::set_gpt3_grdfilename(string fname)
// ...........................................................................
{
    _gpt3_grd_file = fname;
}

} // # namespace ivg















