# include "eop_series.h"

namespace ivg
{
// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Eop_series::Eop_series()
// ...........................................................................
{
    _intrp_mode   = "linear";
    _rg_zont      = true;
    _hf_ocean     = true;
    _hf_ocean_model= "iers2010";
    _ut_libration = true;
    _pm_nutation  = true;
    _nut_type     = "IAU2000/2006";
    _pt_version   = 2015;
}

// .............................................................................
Eop_series::Eop_series(vector<ivg::Eop_series> series) : Eop_series()
// .............................................................................
{
    // calculate required length in advance
    int length=0;
    for(auto &eops: series)
        length += eops.length();
    
    // only if each eop_series contains only one row, this constructor works properly
    if(length != series.size())
        throw logic_error("Eop_series::Eop_series(vector<ivg::Eop_series> series): Unexpected length of accumulated rows in series");
    
    // resize matrices with correct size
    _mjd.resize(length,1);
    _erp.resize(length,3);
    _erp_rates.resize(length,3);
    _nut.resize(length,2);
    _std_erp.resize(length,3);
    _std_erp_rates.resize(length,3);
    _std_nut.resize(length,2);
    _dut_zonal.resize(length,1);
    
    // fill matrices
    int i = 0;
    double last_mjd = 0.0;
    for(auto &eops: series)
    {
        if(eops._mjd(0,0) > last_mjd)
        {
            _mjd.set_sub(i,0,eops._mjd);
            _erp.set_sub(i,0,eops._erp);
            _erp_rates.set_sub(i,0,eops._erp_rates);
            _nut.set_sub(i,0,eops._nut);
            _std_erp.set_sub(i,0,eops._std_erp);
            _std_erp_rates.set_sub(i,0,eops._std_erp_rates);
            _std_nut.set_sub(i,0,eops._std_nut);     
            
            string tmp = eops._type.substr(eops._type.find_last_of("/")+1,9);
            replace_string_in_place(tmp, "_", " " );
            _origin[_mjd(i,0)] = remove_spaces_end(tmp);
            
            last_mjd = _mjd(i,0);
            i++;
        }
        else
            throw logic_error("Eop_series::Eop_series(vector<ivg::Eop_series> series): MJDs of series not in ascending order");
    }
}
// ...........................................................................
Eop_series::Eop_series( const std::string &type, ivg::Matrix data ) : Eop_series()
// ...........................................................................
{
    if(data.cols() != 17)
        throw logic_error("Eop_series::Eop_series( const std::string &type, ivg::Matrix data ): Unexpected number of data-cols, 17 expected");
    
        _type = type;
        _mjd = data.get_col(0);
        _erp = data.get_cols(1,3);
        _erp_rates = data.get_cols(4,6);
        _nut = data.get_cols(7,8);
        _std_erp = data.get_cols(9,11);
        _std_erp_rates = data.get_cols(12,14);
        _std_nut = data.get_cols(15,16);

        _calc_dUT1_zonal_tides();
}
// ...........................................................................
Eop_series::Eop_series( const std::string &file, const std::string &type,
                        const ivg::Date &start, const ivg::Date &end ) : Eop_series()
// ...........................................................................
{
    // c04, finals, etc...
    _type = type;
    
    log<INFO>("*** Initializing ivg::Eop_series with ") % file % " between " % start.get_date_time("DD/MON/YYYY") % " - "  % end.get_date_time("DD/MON/YYYY");
            
    if(start == ivg::Date( ivg::fake_mjd ) && end == ivg::Date( ivg::fake_mjd ))
        read_eop( file, type, ivg::Date(1950,1), ivg::Date(2050,1) );
    else
        read_eop( file, type, start, end );
}

// ...........................................................................
Eop_series::~Eop_series()
// ...........................................................................
{

}
// ===========================================================================
// 		operators
// ===========================================================================
Eop_series Eop_series::operator-( const Eop_series &other ) const
{
    ivg::Eop_series out = *this;
    
    // common start and end epoch 
    double mjd;
    _mjd(0) < other._mjd(0) ? mjd=other._mjd(0) : mjd=_mjd(0);
    ivg::Date start(mjd);
    _mjd(_mjd.rows()-1) > other._mjd(other._mjd.rows()-1) ? mjd=other._mjd(other._mjd.rows()-1) : mjd=_mjd(_mjd.rows()-1);
    ivg::Date end(mjd);
    
    int i=0;
    ivg::Matrix d_erp;
    ivg::Matrix d_nut;
    ivg::Matrix sigma_erp;
    ivg::Matrix sigma_nut;
    std::vector<int> idx;
    while(i<_mjd.rows()&&_mjd(i)<=other._mjd(other._mjd.rows()-1))
    {
        if(_mjd(i)>=other._mjd(0))
        {
            ivg::Date epoch(_mjd(i));
            // ADDED July 2016: "linear", false, in order to be able to subtract UT1 Eop_series within sinex-analyzer (for armin INTs)
            // UT1-UTC has jumps -> problem for not linear interpolation. TODO add leap seconds and remove after interpolation again
            d_erp.append_rows(_erp.get_sub(i,0,i,2)-other.calc_erp(epoch,"linear",false).transpose());
            d_nut.append_rows(_nut.get_sub(i,0,i,1)-other.calc_nut(epoch).transpose());
            // standard deviation: not correct but better than nothing
            ivg::Matrix std = other._std_erp.interpolate(_mjd,epoch.get_double_mjd(),other._intrp_mode);
            
            // HOW DO WE NEED TO CALCULATE THE SIGMA DIFFERENCES?
//            sigma_erp.append_rows(((_std_erp.get_sub(i,0,i,2)^2.0)+(std^2.0))^0.5);
            sigma_erp.append_rows(((std^2.0))^0.5);
            std = other._std_nut.interpolate(_mjd,epoch.get_double_mjd(),other._intrp_mode);
//            sigma_nut.append_rows(((_std_nut.get_sub(i,0,i,1)^2.0)+(std^2.0))^0.5);
            sigma_nut.append_rows(((std^2.0))^0.5);
            
            idx.push_back(i);
        }
        ++i;
    }
    
    // if there is no overlapping
    if(i == 0)
        throw runtime_error("ERROR: No overlapping for Eop_series subtraction.");
    
    out._erp = d_erp;
    out._nut = d_nut;
    out._std_erp = sigma_erp;
    out._std_nut = sigma_nut;
    
    out._erp_rates = ivg::Matrix(d_erp.rows(),3,0.0);
    out._std_erp_rates = ivg::Matrix(d_erp.rows(),3,0.0);
    
    out._mjd = _mjd.get_sub(idx,{0});
    
    return out;
}
// ===========================================================================
// public methods
// ===========================================================================

void Eop_series::init( const string &intrp_mode, const bool &rg_zont,
                       const bool &hf_ocean, const string &hf_ocean_model,
                       const bool &ut_libration, const bool &pm_nutation, const string nut_type,
                       int pt_version )
{
    if( intrp_mode == "linear" ||intrp_mode == "cspline" || intrp_mode == "polynomial" )
    {
        _intrp_mode = intrp_mode;
    }
    else
    {
        stringstream errormessage;
        errormessage << "void Eop_series::init: "
                     << "ERROR: unknown interpolation type: " << intrp_mode
                     << "! Exiting";
        throw logic_error( errormessage.str() );
    }

    _rg_zont = rg_zont;
    _hf_ocean = hf_ocean;
    _hf_ocean_model = hf_ocean_model;
    _ut_libration = ut_libration;
    _pm_nutation = pm_nutation;

    if( nut_type == "IAU2000/2006" || nut_type == "FCN" || nut_type == "MODFILE")
    {
        _nut_type = nut_type;
    }
    else
    {
        stringstream errormessage;
        errormessage << "void Eop_series::init: "
                     << "ERROR: unknown nutation type: " << nut_type
                     << "! Exiting";
        throw logic_error( errormessage.str() );
    }

    _pt_version = pt_version;
    
    // check if last mjd of eop series is not older than pole tide version
    double year = ivg::Date(_mjd(_mjd.rows()-1)).get_decimal_date();
    if(_pt_version == 2003 && year > 2003.0 )
        log<WARNING>( "!!! Eop_series::init(...): Requesting an epoch out of the range 1975.0-2003.0 in version 2003. Extrapolation performed." );
    else if( _pt_version == 2010 && year > 2010.0 )
        log<WARNING>( "!!! Eop_series::init(...): Requesting an epoch out of the range 1975.0-2010.0 in version 2010. Extrapolation performed." );
    else if( _pt_version == 2015 && year > 2016.2 )
        log<WARNING>( "!!! Eop_series::init(...): Requesting an epoch out of the range 1970.0-2016.2 in version 2015. Extrapolation performed." );
    
}
// ...........................................................................
ivg::Matrix Eop_series::calc_mean_pole( const ivg::Date &epoch ) const
// ...........................................................................
{
/*    // load model coefficients from IERS Conventions 2010
    ivg::Matrix M;
    M.load_ascii( "/data/bakkari/iers2010_ptide_ascii.ivg" );

    // calculate mean pole with model coefficients from IERS Conventions 2010
    ivg::Matrix mean_pole_model( 4,2,0.0 );

    if( epoch.get_int_year() < 2010 )
        mean_pole_model = M.get_sub( 0,0,3,1) ;
    else
        mean_pole_model = M.get_sub( 0,2,3,3) ;

    // calculate polynomial at current epoch
    ivg::Matrix dt( 1,4,1.0 );
    dt(1) = epoch.get_decimal_date()-2000.0;
    dt(1) = epoch.get_tmp_date()-2000.0;

    for( int i=2; i<4; ++i )
        dt(i) = pow( dt(1),i );

    // calculate mean pole
    ivg::Matrix mean_pole;
    mean_pole = dt * mean_pole_model;

    // convert mas -> rad
    mean_pole *= ivg::mas2rad;
    
    ivg::Matrix mean_pole0 = mean_pole;
*/
    // -TA- (2015.07.07) make use of IERS subroutine which is available since
    //                   the 2015.06.19 Conventions update; differences to the
    //                   code above are approx 1^-10 rad
    
    // check whether versions is valid
    if( !( _pt_version == 2015 || _pt_version == 2010 || _pt_version == 2003  || _pt_version == 2018) )
    {
        stringstream errormessage;
        errormessage << "ivg::Matrix Eop_series::calc_mean_pole( const ivg::Date "
                     << "const int ): "
                     << "wrong version: " << _pt_version 
                     << ". Valid versions: 2003, 2010, 2015! Exiting";
        throw runtime_error( errormessage.str() );
    }


    // calculate mean pole
    ivg::Matrix mean_pole( 4,2,0.0 );
    int err;
    double epo = epoch.get_decimal_date();
    int version = _pt_version;
    iers::iers_cmp_2015_( &version, &epo, mean_pole.data_ptr(), (mean_pole.data_ptr()+1), &err );

    // convert as -> rad
    mean_pole *= ivg::as2rad;

    // give warning depending on error code
    if( err != 0 )
    {
       if( err == -1 )
       {
           stringstream errormessage;
           errormessage << "ivg::Matrix Eop_series::calc_mean_pole( const ivg::Date "
                        << "const int ): "
                        << "epoch for mean pole calculation < 1970: " <<  epoch.get_decimal_date()
                        << "! Exiting";
           throw runtime_error( errormessage.str() );
       }
       else
       {
          std::vector<std::string> errors 
             { "",
               "Requesting an epoch out of the range 1975.0-2003.0 in version 2003. Extrapolation performed.",  
               "Requesting an epoch out of the range 1975.0-2010.0 in version 2010. Extrapolation performed.",  
               "Requesting an epoch out of the range 1970.0-2016.2 in version 2015. Extrapolation performed."
             };
             // error warning moved to eop series init() in order to get warning once
             // log<WARNING>( "!!! calc_mean_pole: " ) % errors.at( err );
       }
    }

    return mean_pole;
}

// ...........................................................................
ivg::Matrix Eop_series::calc_erp( const ivg::Date &epoch ) 
// ...........................................................................
{
  if (!(epoch.get_double_mjd() == _last_erp_epoch)){
    
    _last_erp=calc_erp( epoch, _intrp_mode, _rg_zont );
    _last_erp_epoch=epoch.get_double_mjd();
  }
  return _last_erp;
}

// ...........................................................................
ivg::Matrix Eop_series::calc_erp( const ivg::Date &epoch,
                                  const std::string &type,
                                  const bool &rg_zont )  const
// ...........................................................................
{
    if( _mjd.rows() < 1 )
    {
        stringstream errormessage;
        errormessage << "ivg::Matrix Eop_series::calc_erp( ivg::Date "
                     << "epoch, std::string type ): "
                     << "no EOP file has been read, call read_eop first! Exiting";
        throw runtime_error( errormessage.str() );
    }

    // reduce effect of zonal tides
    ivg::Matrix erp0 = _erp;
    if( rg_zont )
        erp0.set_sub( 0,2,( erp0( ":",2)-_dut_zonal ) );

    if(epoch.get_double_mjd() >= _mjd(_mjd.rows()-1)+2.0 )
        throw runtime_error("ivg::Matrix Eop_series::calc_erp(...): Not enough data available for interpolation. Update EOP file.");




    // interpolate
    ivg::Matrix erp = (erp0.interpolate( _mjd, epoch.get_double_mjd(),
                                         type )).transpose();
    //erp(2)-=(epoch.get_leap_sec()-ivg::Date(_mjd(0)).get_leap_sec())*ivg::s2rad;
    // restore effect of zonal tides
    if( rg_zont )
    {
        vector<double> dut( 3,0.0 );
        double t = (epoch.get_double_mjd()-51544.0)/36525.0;
        iers::rg_zont2_( &t, &dut[0], &dut[1], &dut[2] );
        erp(2) += dut.at( 0 )*ivg::s2rad;
    }
    
    return erp;
    
}

// ...........................................................................
ivg::Matrix Eop_series::calc_nut( const ivg::Date &epoch ) const
// ...........................................................................
{
    return calc_nut( epoch, _intrp_mode, _nut_type );
}

// ...........................................................................
ivg::Matrix Eop_series::calc_nut( const ivg::Date &epoch,
                                  const std::string &int_type,
                                  const std::string &nut_type  ) const
// ...........................................................................
{

    ivg::Matrix nut(2,1,0.0);
    if (nut_type == "MODFILE")
        nut = (_nut.interpolate(_mjd, epoch.get_double_mjd(), int_type)).transpose();
    else if (nut_type == "FCN")
        nut = calc_free_core_nutation_iers(epoch);

    return nut;
}

// ....................................................................................................
ivg::Matrix Eop_series::calc_wobble_params( const ivg::Date &epoch,
        const ivg::Matrix &mean_pole, const ivg::Matrix &erp ) const
// ....................................................................................................
{
    // wobble parameters
    double m1 = erp(0) - mean_pole(0);
    double m2 = -( erp(1) - mean_pole(1) );

    ivg::Matrix w(2,1,0.0);
    w(0,0) = m1;
    w(1,0) = m2;

    return w;
}
// ....................................................................................................
void Eop_series::read_eop( const std::string &file, const std::string &type,
                           const ivg::Date &start, const ivg::Date &end )
// ....................................................................................................
{
    ivg::Matrix data;
    if( type == "C04" )
        data = ivg::parser::c04(file, start, end);
    else if( type == "finals")
        data = ivg::parser::finals(file, start, end);
    else if( type == "cs_erp") // calc_solve_erp
        data = ivg::parser::cs_erp(file, start, end);
    else if( type == "eops") // eops format
        data = ivg::parser::eops(file, start, end);
    else if( type == "igs") // igs erp format
        data = ivg::parser::igs_erp(file, start, end);
    else
        throw logic_error( "void Eop_series::read_eop( std::string file ): ERROR: unknown type: "+type+"! Exiting" );

    //If the selected file doesnt contain data, abort!
    if(data.rows()==0)
        throw logic_error( "void Eop_series::read_eop( std::string file ): ERROR: No EOP data existent! Maybe file not up to date! ("+file+")" );

    _mjd = data.get_col(0);
    _erp = data.get_cols(1,3);
    _nut = data.get_cols(4,5);
    _std_erp = data.get_cols(6,8);
    _std_nut = data.get_cols(9,10);
    
    _calc_dUT1_zonal_tides();

}
// ....................................................................................................
ivg::Matrix Eop_series::calc_subdaily_erp_iers( const ivg::Date &epoch ) const
// ....................................................................................................
{
  return calc_subdaily_erp_iers( epoch, _hf_ocean, _hf_ocean_model, _ut_libration,
                                   _pm_nutation );
}

// ....................................................................................................
ivg::Matrix Eop_series::calc_subdaily_erp_iers( const ivg::Date &epoch,
        const bool &ocean,
	const std::string &ocean_model,				
        const bool &ut_libration,
        const bool &pm_nutation ) const
// ....................................................................................................
{
//    double mjd = epoch.get_double_mjd();
  double mjd = epoch.get_mjd_tt();

    // effect from ocean tidal model (Ch. 8)
    ivg::Matrix erp_ocean( 3,1,0.0 );
    if( ocean )
    {
        if( ocean_model == "iers2010" )
	  iers::ortho_eop_( &mjd, erp_ocean.data_ptr() );
	else
	  {
	    double dt;
	    dt=epoch.get_leap_sec()+iers::d_tai_tt_sec;
	    iers::calc_hfeop_(&mjd,(char *) ocean_model.c_str(),&dt,erp_ocean.data_ptr() );
	  }
    }
    // tri-axiality of the Earth
    // (1) UT1
    ivg::Matrix ut_lib( 3,1,0.0);
    if( ut_libration )
    {
        vector<double> ut( 2,0.0 );
        iers::utlibr_( &mjd, &ut[0], &ut[1] );
        ut_lib( 2 ) = ut.at( 0 );
    }

    // (2) subdaily nutations in PM
    ivg::Matrix pm_nut( 3,1,0.0);
    if( pm_nutation )
    {
        iers::pmsdnut2_( &mjd, pm_nut.data_ptr() );
    }

    ivg::Matrix sd_erp = erp_ocean+ut_lib+pm_nut;
    sd_erp( 0 ) *= ivg::mas2rad*1e-3;
    sd_erp( 1 ) *= ivg::mas2rad*1e-3;
    sd_erp( 2 ) *= ivg::s2rad*1e-6;

    return( sd_erp );
}

// ....................................................................................................
ivg::Matrix Eop_series::calc_free_core_nutation_iers( const ivg::Date &epoch )
const
// ....................................................................................................
{
    double mjd = epoch.get_double_mjd();

    double x, y, std_x, std_y;
    iers::fcnnut_( &mjd, &x, &y, &std_x, &std_y );

    ivg::Matrix fcn( 2,1,0.0 );
    fcn(0) = x*ivg::mas2rad;
    fcn(1) = y*ivg::mas2rad;

    return fcn;
}

// ....................................................................................................
ivg::Matrix Eop_series::form_crf2trf( const ivg::Date &epoch,
                                      const bool derivative_flag, Partials_t2c * deriv_ptr ) const
{
  return form_crf2trf( epoch, _nut_type, _intrp_mode, _rg_zont, _hf_ocean,_hf_ocean_model,
                         _ut_libration, _pm_nutation, derivative_flag, deriv_ptr );
}
// ....................................................................................................

// ....................................................................................................
ivg::Matrix Eop_series::form_crf2trf( const ivg::Date &epoch,
                                      const string &nut_type,
                                      const std::string &interp_type, const bool &rg_zont, const bool &ortho_eop,
				      const std::string &hf_ocean_model, const bool &ut1libr,
                                      const bool &pmsdnut2, const bool derivative_flag,
                                      Partials_t2c * deriv_ptr ) const
// ....................................................................................................
{
    ivg::Matrix c2t( 3,3,0.0 );
    double rc2t[3][3];
    double mjd = epoch.get_double_mjd();

    // determine polar motion
    ivg::Matrix erp = calc_erp( epoch, interp_type, rg_zont );
    
    ivg::Matrix d_erp = calc_subdaily_erp_iers( epoch, ortho_eop, hf_ocean_model,ut1libr,
                        pmsdnut2 );
    erp += d_erp;
    // determine nutation corrections to IAU2000/2006
    ivg::Matrix dnut( 2,1,0.0 );
    if(nut_type == "MODFILE"||nut_type == "FCN")
    {
        dnut = calc_nut( epoch, interp_type, nut_type );
        c2t = form_crf2trf( epoch, erp, dnut, derivative_flag, deriv_ptr );
    }
    else if(nut_type == "IAU2000/2006")
        c2t = form_crf2trf( epoch, erp, derivative_flag, deriv_ptr );
    else
    {
        stringstream errormessage;
        errormessage <<
                     "ivg::Matrix Eop_series::form_crf2trf( ivg::Date epoch ): "
                     << "ERROR: unknown nutation type: " << nut_type
                     << "! Exiting";
        throw logic_error( errormessage.str() );
    }
    
    return c2t;
}

// ....................................................................................................
ivg::Matrix Eop_series::form_crf2trf( const ivg::Date &epoch,
                                      const ivg::Matrix &erp, const bool derivative_flag,
                                      Partials_t2c * deriv_ptr ) const
// ....................................................................................................
{
    ivg::Matrix dnut( 2,1,0.0 );
    ivg::Matrix c2t = form_crf2trf( epoch, erp, dnut, derivative_flag,
                                    deriv_ptr );
    return c2t;
}

// ....................................................................................................
ivg::Matrix Eop_series::form_crf2trf( const ivg::Date &epoch,
                                      const ivg::Matrix &pm_ut1, const ivg::Matrix &dnut,
                                      const bool derivative_flag, Partials_t2c * deriv_ptr ) const
// ....................................................................................................
{
    
    tictoc time;
    time.tic();
    
    ivg::Matrix c2t( 3,3,0.0 );
    double mjd = epoch.get_double_mjd();
    double jd_utc_0h = epoch.get_int_mjd()+2400000.5;
    double frac_utc = epoch.get_frac_day();

    double frac_tt = epoch.get_jd_tt()-jd_utc_0h;

    ivg::Matrix erp = pm_ut1;
    erp(2) *= ivg::rad2s;
    erp(2) += epoch.get_leap_sec();
    
    // nutation corrections should be applied => stepwise procedure

    // calculate precession/nutation and CIO locator
    double pn_x, pn_y, s_cio;
    iauXy06( jd_utc_0h, frac_tt, &pn_x, &pn_y );
    s_cio = iauS06a( jd_utc_0h, frac_tt );
    // apply nutation corrections
    double x_angle = asin( pn_x );
    double y_angle = asin( pn_y );
    x_angle += dnut(0);
    y_angle += dnut(1);
    pn_x = sin( x_angle );
    pn_y = sin( y_angle );

    // rotate to intermediate frame
    double rc2i[3][3];
    iauC2ixys( pn_x, pn_y, s_cio, rc2i );
    ivg::Matrix c2i = _2darr2mat( rc2i );

    // Earth rotation angle
    double era = iauEra00( jd_utc_0h, frac_utc+erp(2)/86400.0 );
    ivg::Matrix spin( 3,3,0.0 );
    spin.rot3D_z( era );

    // TIO locator
    double s_tio  = iauSp00( 2451545.0, epoch.get_jd_tt()-2451545.0 );

    // polar motion
    double rpom[3][3];
    iauPom00( erp(0), erp(1), s_tio, rpom );
    ivg::Matrix pom = _2darr2mat( rpom );

    // celestial to terrestrial rotation matrix
    c2t = pom*spin*c2i;

    // build matrices with partial derivatives wrt UT1, PM and nutation (XY))
    // several parameters have to be re-calculated as we use the SOFA routines above
    // and now we need internals of theese functions
    if( derivative_flag )
    {
        // UT1
        // (1) inner partial derivative of ERA wrt UT1 (see era00.c from SOFA)
        double dera_dut = 0.00273781191135448+1.0;
        // (2) derivative of the spin matrix
        ivg::Matrix dspin_dut = _3d_deriv_rotation_matrix( 2,era )*dera_dut;
        // (3) derivative of the C2T matrix wrt UT1
        ivg::Matrix dc2t_dut = pom*dspin_dut*c2i;
        deriv_ptr->ut1 = dc2t_dut.transpose();

        // polar motion
        ivg::Matrix rs( 3,3,0.0 );
        rs.rot3D_z( -s_tio );
        ivg::Matrix rx( 3,3,0.0 );
        rx.rot3D_y( erp(0) );
        ivg::Matrix ry( 3,3,0.0 );
        ry.rot3D_x( erp(1) );

        // test what we have re-constructed
        //(rs*rx*ry).transpose().show(21);
        //(ry.transpose()*rx.transpose()*rs.transpose()).show(21);
        //pom.show(21);

        deriv_ptr->pmx = ( ( rs*_3d_deriv_rotation_matrix( 1,erp(0) )*ry ).transpose()*spin*c2i ).transpose();
        deriv_ptr->pmy = ( ( rs*rx*_3d_deriv_rotation_matrix( 0,erp(1) ) ).transpose()*spin*c2i ).transpose();

        // nutation
        // Obtain the spherical angles E and d. (see SOFA c2ixys.c)
        double r2 = x_angle*x_angle + y_angle*y_angle;
        double e = (r2 != 0.0) ? atan2(y_angle, x_angle) : 0.0;
        double d = atan(sqrt(r2 / (1.0 - r2)));

        ivg::Matrix r_1( 3,3,0.0 );
        ivg::Matrix r_2( 3,3,0.0 );
        ivg::Matrix r_3( 3,3,0.0 );
        ivg::Matrix r_4( 3,3,0.0 );
        //r_1.rot3D_z( -(e+s_cio) );  replaced by two rotations
        r_1.rot3D_z( -e );
        r_4.rot3D_z( -s_cio );
        r_2.rot3D_y( d );
        r_3.rot3D_z( e );
        ivg::Matrix pn = ( r_4*r_1*r_2*r_3 );
        // test: below 1e-12, max. deviation for c2i(0,2)
        //(pn-c2i).show(21);

        // build partials
        // (1) spherical angles and cio_locator wrt X and Y
        double tmp  = sqrt( 1.0-( pow( x_angle,2.0 )+pow( y_angle,2.0 ) ) );
        double dedX = -y_angle/( pow( x_angle,2.0 )+pow( y_angle,2.0 ) );
        double dsdX = -y_angle/2.0;
        double dddX =  x_angle/( tmp*sqrt(pow( x_angle,2.0 )+pow( y_angle,2.0 ) ) );
        double dedY =  x_angle/( pow( x_angle,2.0 )+pow( y_angle,2.0 ) );
        double dsdY = -x_angle/2.0;
        double dddY =  y_angle/( tmp*sqrt( pow( x_angle,2.0 )+pow( y_angle,2.0 ) ) );

        // (2) C2I matrix wrt X and Y
        deriv_ptr->pnx = ( pom*spin*
                           ( _3d_deriv_rotation_matrix( 2,-s_cio )*r_1*r_2*r_3*dsdX*(-1.0)
                             +r_4*_3d_deriv_rotation_matrix( 2,-e )*r_2*r_3*dedX*(-1.0)
                             +r_4*r_1*_3d_deriv_rotation_matrix( 1,d )*r_3*dddX
                             +r_4*r_1*r_2*_3d_deriv_rotation_matrix( 2,e )*dedX ) ).transpose();
        deriv_ptr->pny = ( pom*spin*
                           ( _3d_deriv_rotation_matrix( 2,-s_cio )*r_1*r_2*r_3*dsdY*(-1.0)
                             +r_4*_3d_deriv_rotation_matrix( 2,-e )*r_2*r_3*dedY*(-1.0)
                             +r_4*r_1*_3d_deriv_rotation_matrix( 1,d )*r_3*dddY
                             +r_4*r_1*r_2*_3d_deriv_rotation_matrix( 2,e )*dedY ) ).transpose();
    }

    return c2t;
}

// .............................................................................
void ivg::Eop_series::merge(const ivg::Eop_series &other)
// .............................................................................
{
    // merge two Eop_series. Three opssible cases: (1) new series is entrirely
    // before, or (2) after the old series. In case (3) the epochs are mixed =>
    // data has to be sorted
    
    // in case of merging with uninitialized eop series
    if(_mjd.rows() == 0)
        (*this) = other;
    else
    {
        // determine range of new EOPs
        double min_epoch = other._mjd(0);
        double max_epoch = other._mjd(other._mjd.rows()-1);

        // initialize new Eop_series once again with setting of current EOP-series
        ivg::Eop_series new_eops = other;
        new_eops.init(_intrp_mode,_rg_zont,_hf_ocean,_hf_ocean_model,_ut_libration,_pm_nutation, 
                      _nut_type,_pt_version);
        // case (1)
        if(max_epoch<_mjd(0))
        {        
            new_eops._mjd.append_rows(_mjd);
            new_eops._erp.append_rows(_erp);
            new_eops._erp_rates.append_rows(_erp_rates);  
            new_eops._nut.append_rows(_nut);
            new_eops._std_erp.append_rows(_std_erp);
            new_eops._std_erp_rates.append_rows(_std_erp_rates);
            new_eops._std_nut.append_rows(_std_nut);
            new_eops._dut_zonal.append_rows(_dut_zonal);

            *this = new_eops;
        }
        // case (2)
        else if(min_epoch>_mjd(_mjd.rows()-1))
        {
            _mjd.append_rows(new_eops._mjd);
            _erp.append_rows(new_eops._erp);
            _erp_rates.append_rows(new_eops._erp_rates);
            _nut.append_rows(new_eops._nut);
            _std_erp.append_rows(new_eops._std_erp);
            _std_erp_rates.append_rows(new_eops._std_erp_rates);
            _std_nut.append_rows(new_eops._std_nut);
            _dut_zonal.append_rows(new_eops._dut_zonal);
        }
        // case (3)
        else
        {
            // get data from both series in matrices
            ivg::Matrix data0 = get_data();
            if(_rg_zont)
                data0.append_cols(_dut_zonal);
            ivg::Matrix data1 = new_eops.get_data();
            if(_rg_zont)
                data1.append_cols(new_eops._dut_zonal);

            // concatenate the data and sort by MJD
            data0.append_rows(data1);
            data0.sort_cols(0);

            // check whether there are identical epochs
            _mjd = data0(":",0);
            vector<int> idx = _mjd.diff().find_idx(0.0);
            if(idx.size()>0) // identical MJDs => calculate weighted means
            {            
                std::vector<int> idx1(data0.cols()-1);
                std::iota(std::begin(idx1),std::end(idx1),0);
                ivg::Matrix d1 = data0.get_sub(idx,idx1);
                std::transform(std::begin(idx),std::end(idx),std::begin(idx),[](int x){return ++x;});
                ivg::Matrix d2 = data0.get_sub(idx,idx1);
                for(int i=1;i<9;++i)
                {
                    for(int j=0;j<d1.rows();++j)
                    {
                        double w1 = pow(d1(j,8+i),-2.0);
                        double w2 = pow(d2(j,8+i),-2.0);
                        // set the new values but keep in mind that we have inremented 
                        // the indexes
                        data0(idx.at(j)-1,i) = (d1(j,i)*w1+d2(j,i)*w2)/(w1+w2);
                        data0(idx.at(j)-1,i+8) = sqrt( w1*pow(d1(j,i),2.0)
                                                       +w2*pow(d2(j,i),2.0)
                                                       -1.0/(w1+w2)*pow( w1*d1(j,i)
                                                                        +w2*d2(j,i),2));
                    }
                }
                data0.rem_r(idx);
            }

            // replace the EOPs of this
            int n = data0.rows()-1;
            _mjd = data0(":",0);
            _erp = data0.get_sub(0,1,n,3);
            _erp_rates = data0.get_sub(0,4,n,6);
            _nut = data0.get_sub(0,7,n,8);
            _std_erp = data0.get_sub(0,9,n,11);
            _std_erp_rates = data0.get_sub(0,12,n,14);
            _std_nut = data0.get_sub(0,15,n,16);
            if(_rg_zont)
                _dut_zonal = data0.get_sub(0,17,n,17);
        }
    }
}


// .............................................................................
void ivg::Eop_series::replace(ivg::Eop_series &other, bool wob, bool ut1, bool nut){
// .............................................................................
    
    
    std::stringstream ss;
    
    ss << "EOP series: replacing ";
    if(wob)
        ss << "polar motion ";
    if(ut1)
        ss << "UT1-UTC ";
    if(nut)
        ss << "nutation ";

    log<WARNING>("!!! ") % ss.str();
    
    for(unsigned i = 0; i < this->_mjd.length(); ++i){
        ivg::Date epoch( this->_mjd(i) );
                
        ivg::Matrix erp = other.calc_erp(epoch); // X-pole, Y-pole, UT1-UTC
        ivg::Matrix std = other._std_erp.interpolate(_mjd, epoch.get_double_mjd(), "nearest_neighbor" );
        
        if(wob){
            this->_erp(i,0) = erp(0);
            this->_erp(i,1) = erp(1);
            
            this->_std_erp(i,0) = std(0);
            this->_std_erp(i,1) = std(1);

        }
        if(ut1){
            this->_erp(i,2) = erp(2);
            
            this->_std_erp(i,2) = std(2);
        }
        
        if(nut){
            ivg::Matrix dnut = other.calc_nut(epoch);
            ivg::Matrix std = other._std_nut.interpolate(_mjd, epoch.get_double_mjd(), "nearest_neighbor" );
            
            this->_nut(i,0) = dnut(0);
            this->_nut(i,1) = dnut(1);
            
            this->_std_nut(i,0) = std(0);
            this->_std_nut(i,1) = std(1);
        }
    }
    
    if(wob){
        this->_pt_version = other._pt_version;
        this->_pm_nutation = other._pm_nutation;
    }
    if(ut1){
        this->_ut_libration = other._ut_libration;
        this->_hf_ocean = other._hf_ocean;
	this->_hf_ocean_model = other._hf_ocean_model;
        this->_rg_zont = other._rg_zont;
        _calc_dUT1_zonal_tides();
    }
    
    if(nut){
        this->_nut_type = other._nut_type;
    }
 
}

// ....................................................................................................
template <size_t rows, size_t cols>
ivg::Matrix  Eop_series::_2darr2mat( const double (&array)[rows][cols] ) const
// ....................................................................................................
{
    ivg::Matrix mat( rows,cols,0.0 );
    for( int col=0; col<cols; ++col )
    {
        for( int row=0; row<rows; ++row )
        {
            mat( row,col ) = array[row][col];
        }
    }

    return mat;
}

// ....................................................................................................
void Eop_series::_calc_dUT1_zonal_tides()
// ....................................................................................................
{
    // calculate impact of zonal tides on UT1 (independent of _rg_zont)
    _dut_zonal.resize( _erp.rows(),1,0.0 );
    vector<double> dut( 3,0.0 );
    double t;

    for( int i=0; i<_erp.rows(); ++i )
    {
        t = (_mjd(i)-51544.0)/36525.0;
        iers::rg_zont2_( &t, &dut[0], &dut[1], &dut[2] );

        // scale to radian
        _dut_zonal(i) = dut[0]*ivg::s2rad;
    }
}

// ....................................................................................................
ivg::Matrix Eop_series::_3d_deriv_rotation_matrix( const int axis,
        const double phi ) const
// ....................................................................................................
{
    ivg::Matrix r( 3,3,0.0 );
    if( axis == 0 )
    {
        r(1,1) = -sin( phi );
        r(1,2) =  cos( phi );
        r(2,1) = -cos( phi );
        r(2,2) = -sin( phi );
    }
    else if( axis == 1 )
    {
        r(0,0) = -sin( phi );
        r(0,2) = -cos( phi );
        r(2,0) =  cos( phi );
        r(2,2) = -sin( phi );
    }
    else if( axis == 2 )
    {
        r(0,0) = -sin( phi );
        r(0,1) =  cos( phi );
        r(1,0) = -cos( phi );
        r(1,1) = -sin( phi );
    }
    else
    {
        stringstream errormessage;
        errormessage
                << "ivg::Matrix Eop_series::_3d_deriv_rotation_matrix( int axis, "
                << "double phi ): ERROR: wrong axis index: " << axis
                << " != (x=0, y=2, z=2 )! Exiting";
        throw logic_error( errormessage.str() );
    }

    return r;
}
// ...........................................................................
void Eop_series::show(string out)
// ...........................................................................
{
    cerr << "++++++++++++ Eop_series.show(" << out << ") +++++++++++++++++" << endl;
    cerr << "Initialized with type: " << _type << endl;
    cerr << "MJD / XPO / YPO / UT1-UTC / NUTX / NUTY / XPOR / YPOR / UT1R" << endl;
    ivg::Matrix show(_mjd);
    show.append_cols(_erp);
    show.append_cols(_nut);
    if(_erp_rates.cols() == _mjd.cols() && _erp_rates.rows() == _mjd.rows())
        show.append_cols(_erp_rates);
//    show.append_cols(_std_erp);
//    show.append_cols(_std_nut);
    
    show.show(15);
    cout << "------------ Eop_series.show(" << out << ") -----------------" << endl;
}
// ...........................................................................
ivg::Matrix Eop_series::get_time_series( string eop, int order )
// ...........................................................................
{
   ivg::Matrix out = _mjd;
   
    int column = 0;
    if(eop == "xpo" && order == 0)
    {
        out.append_cols( _erp.get_col(0)*ivg::rad2mas );
        out.append_cols( _std_erp.get_col(0)*ivg::rad2mas );
    }
    else if(eop == "xpo" && order == 1)
    {
        out.append_cols( _erp_rates.get_col(0)*ivg::rad2mas );
        out.append_cols( _std_erp_rates.get_col(0)*ivg::rad2mas );
    }
    else if(eop == "ypo" && order == 0)
    {
        out.append_cols( _erp.get_col(1)*ivg::rad2mas );
        out.append_cols( _std_erp.get_col(1)*ivg::rad2mas );
    }
    else if(eop == "ypo" && order == 1)
    {
        out.append_cols( _erp_rates.get_col(1)*ivg::rad2mas );
        out.append_cols( _std_erp_rates.get_col(1)*ivg::rad2mas );
    }
    else if(eop == "ut1" && order == 0)
    {
        out.append_cols( _erp.get_col(2)*ivg::rad2s );
        out.append_cols( _std_erp.get_col(2)*ivg::rad2s );
    }
    else if(eop == "ut1" && order == 1)
    {
        out.append_cols( _erp_rates.get_col(2)*ivg::rad2s*1000 );
        out.append_cols( _std_erp_rates.get_col(2)*ivg::rad2s*1000 );
    }
    else if(eop == "nutx" && order == 0)
    {
       out.append_cols( _nut.get_col(0)*ivg::rad2mas );
       out.append_cols( _std_nut.get_col(0)*ivg::rad2mas );
    }
    else if(eop == "nuty" && order == 0)
    {
       out.append_cols( _nut.get_col(1)*ivg::rad2mas );
       out.append_cols( _std_nut.get_col(1)*ivg::rad2mas );
    }
    else
       throw runtime_error("ivg::Matrix Eop_series::get_time_series( string eop, int order ): Unknown eop type or order requested: "+eop);
   
   return out;
}
// ...........................................................................
ivg::Matrix Eop_series::get_data()
// ...........................................................................
{
   ivg::Matrix out(_mjd.rows(),17); 
   out.set_sub(0,0,_mjd);
   out.set_sub(0,1,_erp);
   out.set_sub(0,4,_erp_rates);
   out.set_sub(0,7,_nut);
   out.set_sub(0,9,_std_erp);
   out.set_sub(0,12,_std_erp_rates);
   out.set_sub(0,15,_std_nut);
   
   return out;
}
// ...........................................................................
string Eop_series::get_origin_db(double mjd)
{
    if(_origin[mjd].empty())
        return("UNKNOWN");
    else
        return(_origin[mjd]); 
}
// ...........................................................................
bool Eop_series::is_initialized()
// ...........................................................................
{
    if(_mjd.rows() == 0)
        return false;
    else
        return true;
}

} // namepsace ivg

