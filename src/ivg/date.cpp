#include"date.h"

namespace ivg
{
// static vector to realize correspondence of month names and their integer
// numbers
std::vector<std::string> Date::_month_names = { "JAN", "FEB", "MAR", "APR", "MAY",
                                                "JUN", "JUL", "AUG", "SEP", "OCT",
                                                "NOV", "DEC"
                                              };

// ===========================================================================
// constructors
// ===========================================================================
// ...........................................................................
Date::Date()
// w/o arguments, i.e., 1.1.0000
// ...........................................................................
{
    _year  = 0;
    _month = 1;
    _day   = 1;
    _hour  = 0;
    _min   = 0;
    _sec   = 0.0;
    _doy   = 1;
    _mjd   = -678941;
}

// ...........................................................................
Date::Date( const int y, const int m, const int d, const int h, const int min,
            const double sec )
// Date( 2000, 1, 1, 12, 5, 10.5 )
// ...........................................................................
{
    _year  = y;
    _month = m;
    _day   = d;
    _hour  = h;
    _min   = min;
    _sec   = sec;

    _calc_mjd();
    _calc_doy();
    _init_leap_sec();
}

// ...........................................................................
Date::Date( const int y, const double doy )
// w/ year and day of year where day of yer includes fraction of day
// ...........................................................................
{
    _year = y;

    Date tmp( _year, 1, 1 );
    _mjd = tmp.get_int_mjd() + doy - 1.0;

    _doy = _mjd-tmp.get_int_mjd()+1.0;

    _calc_date();
    _init_leap_sec();
}

// ...........................................................................
Date::Date( const double mjd )
// calculate date and time from decimal MJD
// ...........................................................................
{
    _mjd = mjd;

    _calc_date();
    _calc_doy();
    _init_leap_sec();
}

// ...........................................................................
Date::Date( const std::string yy, const std::string mmm,
            const std::string dd )
// create date object from strings: date( "00", "JAN", "12" )
// ...........................................................................
{
    if( atoi(yy.c_str()) < 70 )
    {
        _year=atoi(yy.c_str())+2000;
    }
    else
    {
        _year=atoi(yy.c_str())+1900;
    }

    std::vector<std::string>::iterator it;
    it = find( _month_names.begin(), _month_names.end(), mmm );
    if( it == _month_names.end() )
    {
        stringstream errormessage;
        errormessage <<
                     "Date::Date( const std::string yy, const std::string mmm, const std::string dd ): "
                     << "ERROR: wrong month name " << mmm << " Exiting";
        throw logic_error( errormessage.str() );
    }
    _month = it-_month_names.begin()+1;

    _day=atoi(dd.c_str());

    _hour  = 0;
    _min   = 0;
    _sec   = 0;

    _calc_mjd();
    _calc_doy();
    _init_leap_sec();
}



// ...........................................................................
Date::Date( const std::string yymmmdd, const std::string type )
// ...........................................................................
{
    if( type == "YYMMMDD" )
    {
        Date tmp(yymmmdd.substr(0,2),yymmmdd.substr(2,3),yymmmdd.substr(5,2));
        *this = tmp;

    }
    else if( type == "SINEX" )
    {
        int y   = atoi( yymmmdd.substr( 0,2 ).c_str() );
        int doy = atoi( yymmmdd.substr( 3,3 ).c_str() );
        int sec = atoi( yymmmdd.substr( 7,5 ).c_str() );

        ( y < 70 ) ? y += 2000 : y += 1900;

        Date tmp( y, (double)doy + ( (double)sec ) / 86400.0 );

        *this=tmp;
    }
    else
    {
        stringstream errormessage;
        errormessage << "Date::Date( string yymmmdd, string type ): "
                     << "ERROR: unknown type " << type << " Exiting";
        throw logic_error( errormessage.str() );
    }
}


// ===========================================================================
// Operators
// ===========================================================================
bool Date::operator==( const Date &other ) const
{
    if( _year == other._year && _doy == other._doy && _month == other._month &&
            _hour == other._hour && _min == other._min && abs( _sec-other._sec)<1e-3 )
    {
        return true;
    }
    else
    {
        //std::cerr << _sec-other._sec << std::endl;
        return false;
    }
}

bool Date::operator<( const Date &other ) const
{
    if( _mjd-other._mjd<-1e-3 )
    {
        return true;
    }
    else
    {
        //std::cerr << _sec-other._sec << std::endl;
        return false;
    }
}

bool Date::operator>( const Date &other ) const
{
    if( _mjd-other._mjd>1e-3 )
    {
        return true;
    }
    else
    {
        //std::cerr << _sec-other._sec << std::endl;
        return false;
    }
}

bool Date::operator<=( const Date &other ) const
{
    if( _mjd <= other._mjd )
    {
        return true;
    }
    else
    {
        //std::cerr << _sec-other._sec << std::endl;
        return false;
    }
}

bool Date::operator>=( const Date &other ) const
{
    if( _mjd >= other._mjd )
    {
        return true;
    }
    else
    {
        //std::cerr << _sec-other._sec << std::endl;
        return false;
    }
}

// ===========================================================================
// public methods
// ===========================================================================

// ...........................................................................
void Date::set_string_date( const std::string yymmmdd )
// ...........................................................................
{
    Date tmp(yymmmdd.substr(0,2),yymmmdd.substr(2,3),yymmmdd.substr(5,2));
    *this = tmp;
}

// ...........................................................................
void Date::set_date( const std::string yy, const std::string mmm,
                     const std::string dd)
// ...........................................................................
{
    Date temp(yy,mmm,dd);
    *this = temp;
}


// ...........................................................................
bool Date::is_leap_year() const
// ...........................................................................
{
    Date tmp( _year, 12, 31 );
    if( tmp.get_int_doy() == 365 )
        return false;
    else
        return true;
}


double Date::get_tmp_date() const
{
    double ifac = 365.0;
    if( is_leap_year() )
        ifac++;

    double t = _year + (_doy + _hour/23.93447 + _min/1440.0 + _sec/86400.0 ) / ((
                   ifac) + 0.2422);

    return t;
}

// ...........................................................................
double Date::get_decimal_date() const
// ...........................................................................
{
    Date tmp( _year, 12, 31 );
    double yDec = _year + ( _doy + _day_frac - 1.0 ) / ((double)
                  tmp.get_int_doy());


    return yDec;
}

// ...........................................................................
void Date::set_decimal_date( const double d )
// ...........................................................................
{
    _year = (int) d;

    Date tmp( _year, 12, 31 );
    double doy = ( d - _year ) * tmp.get_int_doy()+1.0;

    Date newDate( _year, doy );
    *this = newDate;
}

// ...........................................................................
void Date::show() const
// ...........................................................................
{
    cerr <<  get_date_time("DD/MON/YYYY") << endl;
    cout << _year << "-" << _month << "-" << _day << " "
         << _hour << ":" << _min << ":" << _sec << endl;
    cout << _mjd << ", " << _doy << ", " << setprecision(16) << get_decimal_date()
         << endl;
    cout << _leap_sec << endl;
}

// ...........................................................................
void Date::now()
// ...........................................................................
{
    time_t curDate;
    tm *nun;

    curDate = time(0);
    nun = localtime( &curDate );

    Date tmp( nun->tm_year+1900, nun->tm_mon+1, nun->tm_mday, nun->tm_hour,
              nun->tm_min, nun->tm_sec );

    *this = tmp;
}


// ...........................................................................
double Date::get_jd_tt() const
// ...........................................................................
{
    return ( _mjd+_tt_corr_frac_days()+2400000.5 );
}

// ...........................................................................
double Date::get_jd_tdb() const
// ...........................................................................
{
    return ( _mjd+_tdb_corr_frac_days()+2400000.5 );
}

// ...........................................................................
double Date::get_mjd_tt() const
// ...........................................................................
{
    return ( _mjd+_tt_corr_frac_days() );
}

// ...........................................................................
double Date::get_mjd_tdb() const
// ...........................................................................
{
    return ( _mjd+_tdb_corr_frac_days() );
}

// ...........................................................................
double Date::get_greenwich_sidereal_time0() const
// ...........................................................................
{
    double jd_tdb = get_jd_tdb();

    double T = ( int(jd_tdb) - 2451545.0 ) / 36525.0;
    double hsga = 0.2790572733 + 100.002139037801 * T
                  + 1.077592593e-6 * pow( T, 2 )
                  - 7.1759e-11 * pow( T, 3 );

    hsga = hsga - (int)hsga;

    if( hsga <= 0.0 )
        hsga += 1.0;

    double gst0 = hsga * 2.0 * M_PI;
    return gst0;
}

// ...........................................................................
double Date::get_greenwich_sidereal_time() const
// ...........................................................................
{
    double jd_tdb = get_jd_tdb();
    double gst0 = get_greenwich_sidereal_time0();

    double gst = gst0 + ( jd_tdb - int(jd_tdb) - 0.5 )
                 * 2.0 * M_PI * 1.002737909265;
    return gst;

}

// ...........................................................................
double Date::_tt_corr_frac_days() const
// ...........................................................................
{
    return ( (iers::d_tai_tt_sec+_leap_sec)/86400.0 );
}


// ...........................................................................
double Date::_tdb_corr_frac_days() const
// ...........................................................................
// reference: Seidelmann and Fukushima 1992, 'Why new time scales', table 1
{
    double d_utc_tt = _tt_corr_frac_days();

    // elapsed time from J2000.0 in Julian centuries
    double t = ( _mjd-51544.5+d_utc_tt)/36525.0;

    double p = 0.0016568*sin( (35999.37*t+357.5 )*ivg::d2rad )+( 0.0000224*sin( (
                   32964.5*t+246.0 )*ivg::d2rad )
               +( 0.0000138*sin( ( 71998.7*t+355.0 )*ivg::d2rad )+( 0.0000048*sin( (
                           3034.9*t+25.0 )*ivg::d2rad )
                       +( 0.0000047*sin( ( 34777.3*t+230.0 )*ivg::d2rad ) ) ) ) ); // [sec]

    double tdb = d_utc_tt+p/86400.0;

    // modest / VieVS
    double g = ( 357.528+35999.05*t )*ivg::d2rad;
    tdb= d_utc_tt+( 0.001658*sin( g+0.0167*sin(g) ) )/86400.0;

    return tdb;
}



std::string Date::get_date_time( std::string format ) const
// return date and time string in arbitrary format. The following placeholders
// will be filled
// YYYY     - 4 digit year
// YY       - 2 digit year w/ leading zero (YY > 70 for pre 2k dates)
// DOY      - day of year
// MON      - month in uppercase 3 character representation
// MO       - month
// DD       - day of week in 2 character format w/ leading zero
// HH       - hours in 2 character
// MI       - minute in 2 character format w/ leading zero
// SS       - integer second in 2 character format w/ leading zero
// SS.S     - double seconds with given number of decimal places
// SSSSS    - int seconds of day
// SSSSS.SS - double seconds of day with given number of decimal places
{
    char datech [100];

    //
    int y = _year;
    ( y >= 2000 ) ? y -= 2000 : y -= 1900;

    // replace 4 digit year first
    sprintf( datech, "%04d", _year );
    replace_string_in_place( format, "YYYY", std::string( datech ) );

    // replace 2 digit year if YY is remaining
    sprintf( datech, "%02d", y );
    replace_string_in_place( format, "YY", std::string( datech ) );

    sprintf( datech, "%03d", _doy );
    replace_string_in_place( format, "DOY", std::string( datech ) );
    replace_string_in_place( format, "MON", _month_names.at( _month-1 ) );
    sprintf( datech, "%02d", _month );
    replace_string_in_place( format, "MO", std::string( datech ) );
    sprintf( datech, "%02d", _day );
    replace_string_in_place( format, "DD", std::string( datech ) );
    sprintf( datech, "%02d", _hour );
    replace_string_in_place( format, "HH", std::string( datech ) );
    sprintf( datech, "%02d", _min );
    replace_string_in_place( format, "MI", std::string( datech ) );

    // seconds

    // (1.1) seconds of day

    // check whether decimal point is given
    std::size_t found = format.find( "SSSSS.S" );
    if( found!=std::string::npos )
    {
        // then replace fraction of second
        std::size_t idx1 = found;
        std::size_t idx2 = idx1+6;
        while( format.at( idx2 ) == 'S' )
            idx2++;
        string fmt = "%0"+std::to_string(idx2-idx1)+"."+std::to_string(
                         idx2-idx1-6)+"f";

        double sec = (double)_day_frac * 86400.0;
        sprintf( datech, fmt.c_str(), sec );
        std::string datestring = std::string( datech );
        format.replace( idx1, datestring.length(), datestring );
    }
    else
    {
        int sec = _day_frac * 86400.0;
        sprintf( datech, "%05d", sec );
        replace_string_in_place( format, "SSSSS", std::string( datech ) );
    }

    // (1.1) seconds within minute
    found = format.find( "SS.S" );
    if( found!=std::string::npos )
    {
        // then replace fraction of second
        std::size_t idx1 = found;
        std::size_t idx2 = idx1+3;
        while( format.at( idx2 ) == 'S' )
            idx2++;
        string fmt = "%0"+std::to_string(idx2-idx1)+"."+std::to_string(
                         idx2-idx1-3)+"f";
        sprintf( datech, fmt.c_str(), _sec );
        std::string datestring = std::string( datech );
        format.replace( idx1, datestring.length(), datestring );
    }
    else
    {
        sprintf( datech, "%02d", (int)_sec );
        replace_string_in_place( format, "SS", std::string( datech ) );
    }

    return format;
}
// ...........................................................................
double Date::get_qcustomplot_date() const
// ...........................................................................
{
    Date ref_date( 1970, 1.0 );
    double sec_diff = ( _mjd-ref_date.get_double_mjd() )*86400.0;

    return sec_diff;
}
// ...........................................................................
void Date::add_secs( const double secs )
// ...........................................................................
{
    Date newDate( _mjd + (secs / 86400.0) );
    *this = newDate;
}
// ...........................................................................
Date Date::add_days( const double days )
// ...........................................................................
{
    return(ivg::Date( _mjd + days ));
}
// ===========================================================================
// private methods
// ===========================================================================

// ...........................................................................
void Date::_calc_mjd( )
// ...........................................................................
{
    double frac, gyr;
    int iy0, im0;
    int ia, ib;
    int jd;

    /* decimal day fraction */
    frac = ((double)_hour/ 24.0)
           + ((double)_min / 1440.0)
           + (_sec / 86400.0);
    _day_frac = frac;
    /* convert date to format YYYY.MMDDdd  */
    gyr = (double) _year
          + (0.01 * (double)_month)
          + (0.0001 * (double)_day)
          + (0.0001 * frac) + 1.0e-9;
    /* conversion factors */
    if ( _month <= 2 )
    {
        iy0 = _year - 1L;
        im0 = _month + 12;
    }
    else
    {
        iy0 = _year;
        im0 = _month;
    }
    ia = iy0 / 100L;
    ib = 2L - ia + (ia >> 2);
    /* calculate julian date   */
    if ( _year <= 0L )
        jd = (int) ((365.25 * (double) iy0) - 0.75)
             + (int) (30.6001 * (im0 + 1L) )
             + (int) _day + 1720994L;
    else
        jd = (int) (365.25 * (double) iy0)
             + (int) (30.6001 * (double) (im0 + 1L))
             + (int) _day + 1720994L;
    if ( gyr >= 1582.1015 ) /* on or after 15 October 1582   */
        jd += ib;

    _mjd = jd + .5 - 2400000.5 + frac;
    double j_date = (double) jd + frac + 0.5;
    jd = (int) (j_date + 0.5);
    _weekday = (jd + 1L) % 7L;
    //_mjd = j_date - 2400000.5;

}

// ...........................................................................
void Date::_calc_date( )
// ...........................................................................
{
    double frac;
    int jd;
    int ka;
    int kb;
    int kc;
    int kd;
    int ke;
    int ialp;
    double d_hour;
    double d_minute;

    double j_date = _mjd + 2400000.5;
    jd = (int) (j_date + 0.5); /* integer julian date */
    frac = j_date + 0.5 - (double) jd + 1.0e-10; /* day fraction */
    _day_frac = frac;
    ka = (int) jd;
    if ( jd >= 2299161L )
    {
        ialp = ( (double) jd - 1867216.25 ) / ( 36524.25 );
        ka = jd + 1L + ialp - ( ialp >> 2 );
    }
    kb = ka + 1524L;
    kc =  ( (double) kb - 122.1 ) / 365.25;
    kd = (double) kc * 365.25;
    ke = (double) ( kb - kd ) / 30.6001;
    _day = kb - kd - ((int) ( (double) ke * 30.6001 ));
    if ( ke > 13L )
        _month = ke - 13L;
    else
        _month = ke - 1L;
    if ( (_month == 2) && (_day > 28) )
        _day = 29;
    if ( (_month == 2) && (_day == 29) && (ke == 3L) )
        _year = kc - 4716L;
    else if ( _month > 2 )
        _year = kc - 4716L;
    else
        _year = kc - 4715L;
    _hour = d_hour = frac * 24.0; /* hour */
    _min = d_minute =
               ( d_hour - (double) _hour ) * 60.0; /* minute */
    _sec =
        ( d_minute - (double) _min ) * 60.0;/* second */
    _weekday = (jd + 1L) % 7L; /* day of week */
}

// ...........................................................................
void Date::_calc_doy()
// ...........................................................................
{
    if ( _year == ((_year >> 2) << 2) )
        _doy =
            ( ( 275 * _month ) / 9)
            - ((_month + 9) / 12)
            + _day - 30;
    else
        _doy =
            ( ( 275 * _month ) / 9)
            - (((_month + 9) / 12) << 1)
            + _day - 30;

}
// ....................................................................................................
void Date::_init_leap_sec()
// ....................................................................................................
{
    vector<double> ls = ivg::leap_secs;
    ivg::Matrix leap_sec( ls.begin(),ls.end(),ls.size()/2,2 );

    std::vector<int> idx = leap_sec( ":",0 ).find_idx( le, _mjd );
    if( idx.size() > 0 )
    {
        int leap_idx = idx.at( idx.size()-1);
        _leap_sec = leap_sec( leap_idx,1 );
    }
}

}
