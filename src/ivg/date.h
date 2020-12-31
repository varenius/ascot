#ifndef IVGDATE_H
#define IVGDATE_H

#include<vector>
#include<iostream>
#include<iterator>
#include<limits>
#include<string>
#include<fstream>
#include<sstream>
#include<stdexcept>
#include<sys/stat.h>
#include<numeric>
#include<cmath>
#include<time.h>
#include<auxfunc.h>
#include<iers_wrapper.h>
#include<ivg_const.h>


/**
*
* @brief Handling of date and time in various formats
* @author bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

using namespace std;

namespace ivg
{
// ===========================================================================
class Date
// ===========================================================================
{
    public:
        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
         *  \b Description: \n
         *        Default constructor
         *  \param [in] no input parameters needed
         *  \return An instance of the class 'Date'
         */
        Date( );

        /**
        *  \b Description: \n
        *        Constructor using one parameter: Modified Julian Date
        *  \param [in] [int] mjd
        *  \return An instance of the class 'Date'
        */
        Date( const double mjd );

        /**
        *  \b Description: \n
        *        Constructor using two parameters: 4 digit year and day of year (including fraction of day)
        *  \param [in] [int] year
        *              [double] doy
        *  \return An instance of the class 'Date'
        */
        Date( const int yyyy, const double doy );

        /**
        *  \b Description: \n
        *        Constructor using six parameters: year, month, day of month, hours, minutes and seconds
        *  \param [in] [int] year
        *              [int] month
        *              [int] day
        *              [int] h
        *              [int] min
        *              [double] sec
        *  \return An instance of the class 'Date'
        */
        Date( const int y, const int m, const int d, const int h = 0,
              const int min = 0, const double sec = 0.0 );

        /**
        *  \b Description: \n
        *        Constructor using three parameters: 2 digit year (assuming dates > 70 are before 2k), 3 char month string
        *        and day of month
        *  \param [in] [string] year
        *              [string] mmm
        *              [string] dd
        *  \return An instance of the class 'Date'
        */
        Date( const std::string yy, const std::string mmm, const std::string dd );

        /**
        *  \b Description: \n
        *        Constructor using two parameters: date string, and type of string (yymmmdd, SINEX)
        *  \param [in] [string] date
        *              [string] type
        *  \return An instance of the class 'Date'
        */
        Date( const std::string yymmmdd, const std::string type );


        // ==============================================
        // ================ Operators: ==================
        // ==============================================
        /**
         *  \b Description: \n
         *        equals operator '==' using another date object, two objects are considered being equal by 1 ms.
         *  \param [in] [Date] &other: other date
         */
        bool operator==( const Date &other ) const;

        /**
         *  \b Description: \n
         *        less than operator '<' using another date object (two objects are considered being equal by 1 ms).
         *  \param [in] [Date] &other: other date
         */
        bool operator<( const Date &other ) const;
        
        /**
        *  \b Description: \n
        *        greater than operator '>' using another date object (two objects are considered being equal by 1 ms).
        *  \param [in] [Date] &other: other date
        */
        bool operator>( const Date &other ) const;
        
        bool operator>=( const Date &other ) const;
        bool operator<=( const Date &other ) const;
        
        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        /**
        *  \b Description: \n
        *        Method to modify date according to decimal year
        *  \param [in] no input parameters needed
        */
        void set_decimal_date( const double d );

        /**
        *  \b Description: \n
        *        Method to modify date according to 2 character year, 3 character month and day (e.g. "00", "JAN", "01")
        *        assuming that dates > 70 are before 2k
        *  \param [in] no input parameters needed
        */
        void set_date( const std::string yy, const std::string mmm,
                       const std::string dd);

        /**
        *  \b Description: \n
        *        Method to modify date according to one single string which is composed of
        *        2 character year, 3 character month and day (e.g. "00JAN01")
        *        assuming that dates > 70 are before 2k
        *  \param [in] no input parameters needed
        */
        void set_string_date( const std::string yymmmdd );

        /**
        *  \b Description: \n
        *        Method to modify date according to current epoch (local hardware time)
        *  \param [in] no input parameters needed
        */
        void now();

        /**
        *  \b Description: \n
        *        Method to get day of year
        *  \param [in] no input parameters needed
        *  \return [int] day of year
        */
        int get_int_doy() const
        {
            return _doy;
        }

        /**
        *  \b Description: \n
        *        Method to get day of year including fractions of day
        *  \param [in] no input parameters needed
        *  \return [double] day of year
        */
        double get_double_doy() const
        {
            return ( _doy + _day_frac );
        }

        /**
        *  \b Description: \n
        *        Method to get integer year
        *  \param [in] no input parameters needed
        *  \return [int] year
        */
        int get_int_year() const
        {
            return _year;
        }

        /**
        *  \b Description: \n
        *        Method to get integer month
        *  \param [in] no input parameters needed
        *  \return [int] year
        */
        int get_int_month() const
        {
            return _month;
        }

        /**
        *  \b Description: \n
        *        Method to get integer day
        *  \param [in] no input parameters needed
        *  \return [int] year
        */
        int get_int_day() const
        {
            return _day;
        }

        /**
        *  \b Description: \n
        *        Method to get integer hour
        *  \param [in] no input parameters needed
        *  \return [int] year
        */
        int get_int_hour() const
        {
            return _hour;
        }

        /**
        *  \b Description: \n
        *        Method to get integer minute
        *  \param [in] no input parameters needed
        *  \return [int] year
        */
        int get_int_min() const
        {
            return _min;
        }

        /**
        *  \b Description: \n
        *        Method to get double second
        *  \param [in] no input parameters needed
        *  \return [double] sec
        */
        double get_double_sec() const
        {
            return _sec;
        }

        /**
        *  \b Description: \n
        *        Method to get double fraction of day [UTC]
        *  \param [in] no input parameters needed
        *  \return [double] fractions of day
        */
        double get_frac_day() const
        {
            return _day_frac;
        }

        /**
        *  \b Description: \n
        *        Method to get modified julian date
        *  \param [in] no input parameters needed
        *  \return [int] mjd
        */
        int get_int_mjd() const
        {
            return (int)_mjd;
        }

        /**
        *  \b Description: \n
        *        Method to get modified julian date including fractions of day
        *  \param [in] no input parameters needed
        *  \return [double] mjd
        */
        double get_double_mjd() const
        {
            return ( _mjd );
        }

        /**
        *  \b Description: \n
        *        Method to get julian date including fractions of day
        *  \param [in] no input parameters needed
        *  \return [double] jd
        */
        double get_jd() const
        {
            return ( _mjd + 2400000.5 );
        }

        /**
        *  \b Description: \n
        *        Method to get julian date including fractions of day in TT time scale
        *  \param [in] no input parameters needed
        *  \return [double] jd
        */
        double get_jd_tt() const;
        double get_mjd_tt() const;

        /**
        *  \b Description: \n
        *        Method to get julian date including fractions of day in TDB time scale
        *  \param [in] no input parameters needed
        *  \return [double] jd
        */
        double get_jd_tdb() const;
        double get_mjd_tdb() const;

        double get_greenwich_sidereal_time0() const;

        double get_greenwich_sidereal_time() const;


        /**
        *  \b Description: \n
        *        Method to get number of leap seconds for this epoch
        *  \param [in] no input parameters needed
        *  \return [double] TAI-UTC
        */
        double get_leap_sec() const
        {
            return ( _leap_sec );
        }

        /**
        *  \b Description: \n
        *        Method to get day of week
        *  \param [in] no input parameters needed
        *  \return [int] dayofweek
        */
        int get_weekday() const
        {
            return _weekday;
        }

        /**
        *  \b Description: \n
        *        Method to get decimal year including fraction of year
        *  \param [in] no input parameters needed
        *  \return [double] year
        */
        double get_decimal_date() const;

        /**
        *  \b Description: \n
        *        Method to get date as used in VieVS within ptide calculations
        *  \param [in] no input parameters needed
        *  \return [double] date
        */
        double get_tmp_date() const;

        /**
        *  \b Description: \n
        *        Method to get date/time string in arbitrary format. The following placeholders
        *        will be filled
        *        YYYY     - 4 digit year
        *        YY       - 2 digit year w/ leading zero (YY > 70 for pre 2k dates)
        *        DOY      - day of year
        *        MON      - month in uppercase 3 character representation
        *        MO       - month
        *        DD       - day of week in 2 character format w/ leading zero
        *        HH       - hours in 2 character
        *        MI       - minute in 2 character format w/ leading zero
        *        SS       - integer second in 2 character format w/ leading zero
        *        SS.S     - double seconds with given number of decimal places
        *        SSSSS    - int seconds of day
        *        SSSSS.SS - double seconds of day with given number of decimal places
        *  \param [in] format string
        *  \return [std::string] modified format string
        */
        string get_date_time( std::string format ) const;

        /**
        *  \b Description: \n
        *        Method to determine the date for qcustomplot which is related to 1970
        *  \param [in] no input parameters needed
        *  \return [double] milliseconds that have passed since 1970-01-01T00:00:00.000
        */
        double get_qcustomplot_date() const;

        /**
        *  \b Description: \n
        *        Method adds secs to existing date-instance.
        *  \param [in] seconds, eg. 86400.00 for adding one day
        *  \return no return
        */
        void add_secs( const double secs );

        /**
        *  \b Description: \n
        *        Method adds days to existing date-instance.
        *  \param [in] days, eg. 14 for adding 14 days
        *  \return no return
        */
        Date add_days( const double days );

        /**
        *  \b Description: \n
        *        Method to determine if the date belongs to a leap year, i.e., 366 days per year if true
        *  \param [in] no input parameters needed
        *  \return [bool] year_is_leap_year
        */
        bool is_leap_year() const;

        /**
        *  \b Description: \n
        *        Method to display all components of the object in STDOUT
        *  \param [in] no input parameters needed
        */
        void show() const;


    private:

        // ==============================================
        // ======== private MEMBER-Functions: ===========
        // ==============================================


        /**
        *  \b Description: \n
        *        Method to calculate date components from modified julian day
        *  \param [in] no input parameters needed
        */
        void _calc_date( );

        /**
        *  \b Description: \n
        *        Method to calculate modified julian day from year, month, day, hour, minute and second
        *  \param [in] no input parameters needed
        */
        void _calc_mjd( );


        /**
        *  \b Description: \n
        *        Method to calculate day of year year, month and day
        *  \param [in] no input parameters needed
        */
        void _calc_doy( );

        void _init_leap_sec();

        /**
        *  \b Description: \n
        *        Method to calculate difference between UTC and TDB timescale in fraction of days
        *  \param [in] no input parameters needed
        */
        double _tdb_corr_frac_days() const;

        /**
        *  \b Description: \n
        *        Method to calculate difference between UTC and TT timescale in fraction of days
        *  \param [in] no input parameters needed
        */
        double _tt_corr_frac_days() const;

        // ==============================================
        // ======== class variables / attributes: =======
        // ==============================================

        // year, mont, day of month, hour, minute and second
        int _year;
        int _month;
        int _day;
        int _hour;
        int _min;
        double _sec;
        double _mjd;

        // redundant information
        double _day_frac;         // hour minute and second as fraction of day
        int _doy;                // day of year
        int _weekday;            // day of week

        // vector containing 3 character descriptions of month names
        static std::vector<std::string> _month_names;

        double _leap_sec;
};
}
#endif  // IVGDATE_H:
