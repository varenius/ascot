#ifndef OBS_H
#define OBS_H

#include "logger.h"
#include "date.h"
#include "troposphere.h"
#include "matrix.h"
#include "trf.h"
#include "eop_series.h"
#include "auxfunc.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "param.h"
#include <cstdlib>
#include <libconfig.h++>

/**
* @brief Session class
* @author SH - bakkari developer team
* @date 2015-03-24
* @version 0.1
*
* @todo calc_delay: (1) revise time of closest approach,
*                   (2) higher order relativistic time delay for the Sun
*                   (3) "post-model" changes in the baseline
*                   (4) new units for parameters
*/

using namespace std;
using namespace libconfig;

namespace ivg
{

enum delaymodel{ consensus, duev, sekido, simple };

class Scan; // forward declaration
class Session;

// ===========================================================================
class Obs
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
        *  \return An instance of the class 'Obs'
        */
        Obs( );

        /**
        *  s\b Description: \n
        *        Constructor using five parameters:
        *        two station names, source name, epoch, mjd
        *  \param [in] [std::string] station name 1
        *              [std::string] station name 2
        *              [std::string] source name
        *              [ivg::Date] epoch
        *              [double] mjd
        *  \return An instance of the class 'Obs'
        */
        Obs( string s1, string s2, string src, ivg::Date epochIn,
             double mjdIn = -999.99 );

        /**
        *  s\b Description: \n
        *        Constructor using four parameters:
        *        session and scan pointer as well as two station indexes
        *  \param [in] [ivg::Session *] session pointer
        *              [ivg::Scan *] scan pointer
        *              [int] index for station 1 in scan
        *              [int] index for station 2 in scan
        *              [double] mjd
        *  \return An instance of the class 'Obs'
        */
        Obs( ivg::Session * session_ptr, ivg::Scan * scan_ptr,
             int sta1_scan_idx, int sta2_scan_idx );


        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Method to show observation information
        *  \param [in] no input parameters needed
        */
        void show();

        /**
        *  \b Description: \n
        *        Method to set observed delay (ns) and delay rate (ps/sec)
        *  \param [in] [double] observed (baseline-based) group delay (ns)
        *              [double] formal error of observed delay (ns)
        *              [double] geocenter group delay (ns); formal error is the same as for the baseline-based delay
        *              [double] delay rate (ps/sec)
        *              [double] formal error of delay rate (ps/sec)
        */
        void set_delay( double d, double sigmaD, double gcd, double r = -999.99, double sigmaR = -999.99 );
        
        /**
        *  \b Description: \n
        *        Method to set observed singleband delay (ns)
        *  \param [in] [double] observed single band delay (ns)
        *              [double] formal error of observed single band delay (ns)
        */
        void set_sb_delay(double sbd, double sigmaSbd);

        /**
        *  \b Description: \n
        *        Method to set observed phase delay (ns)
        *  \param [in] [double] observed phase delay (ns)
        *              [double] formal error of observed phase delay (ns)
        */
        void set_phase_delay(double phd, double sigmaPhd);
        
        /**
        *  \b Description: \n
        *        Method to set SNR for X and S band
        */
        void set_snr( double snr_bx, double snr_bs );

        /**
        *  \b Description: \n
        *        Method to set ionosphere correction
        *  \param [in] [double] delay ionosphere correction (ns)
        *              [double] formal error of delay ionosphere correction (ns)
        *              [double] delay rate ionosphere correction (ps/sec)
        *              [double] formal error of delay rate ionosphere correction (ps/sec)
        */
        void set_ion_corr( double id, double sigmaId, double ir,
                           double sigmaIr );

        /**
        *  \b Description: \n
        *        Method to set se feed rotation correction
        *  \param [in] [double] feed rotation correction for station 1
        *              [double] feed rotation correction for station 2
        */
        void set_feed_rotation( double feed1, double feed2 );
        /**
        *  \b Description: \n
        *        Method to set cable calibration correction
        *  \param [in] [double] cable calibration correction for station 1
        *              [double] cable calibration correction for station 2
        */
        void set_cable_cal( double cc1, double cc2 );

        /**
        *  \b Description: \n
        *        Method to set scan
        *  \param [in] [ivg::Scan *] scan
        */
        void set_scan( ivg::Scan* scan );
        
        /**
        *  \b Description: \n
        *        Method to set epoch of observation
        *        (Overwrites epoch which was originally taken from _scan) 
        *  \param [in] [ivg::Date] epoch of obs 
        */
        void set_epoch( ivg::Date epoch );

        /**
        *  \b Description: \n
        *        Method to set scan
        *  \param [in] [ivg::Scan *] scan
        */
        ivg::Scan* get_scan( ){ return _scan; };

        /**
        *  \b Description: \n
        *        Method to get scan index for station 1 or 2
        *  \param [in] [int] station 1 or station 2
        *  \return [int] scan index for station 1 or 2
        */
        int get_scan_idx( int sta_idx );

        /**
        *  \b Description: \n
        *        Method to get 'observed-minus-computed'
        *  \param [in] no input parameters needed
        *  \return [double] 'observed-minus-computed'
        */
        double get_observed_minus_computed()
        {
            return _oc;
        };

        /**
        *  \b Description: \n
        *        Method to set 'observed-minus-computed'
        *  \param [in][double] (simulated) observed minus computed value
        */
        void set_observed_minus_computed( double oc )
        {
            _oc = oc;
        };

        /**
        *  \b Description: \n
        *        Method to get the variance of the observation
        *  \param [in] [bool] ion: use delay ionosphere correction formal error?
        *  \return [double] variance of the observation
        */
        double get_obs_variance( bool ion=true, bool phase=false , bool sb=false);

        /**
        *  \b Description: \n
        *        Method to get the standard deviation of the observation
        *  \param [in] [bool] ion: use delay ionosphere correction formal error?
        *  \return [double] standard deviation of the observation
        */
        double get_obs_std( bool ion=true );

        /**
        *  \b Description: \n
        *        Method to get names of station 1 and 2
        *  \param [in] no input parameters needed
        *  \param [out] [std::string] name of station 1
        *               [std::string] name of station 2
        */
        void get_station_names( string & sta1, string&  sta2, ivg::staname staname = ivg::staname::ivs_name ) const;

        /**
        *  \b Description: \n
        *        Method to get source name
        *  \param [in] no input parameters needed
        *  \param [out] [std::string] source name
        */
        void get_source_name( string & src ) const;

        /**
        *  \b Description: \n
        *        Method to get azimuth and elevation
        *  \param [in] [int] select station 1 or station 2
        *  \return [ivg::Matrix] azimuth and elevation
        */
        ivg::Matrix get_az_el( int idx );

        /**
        *  \b Description: \n
        *        Method to get the mapping functions (hydrostatic and wet part) 
        *  \param [in] [int] select station 1 or station 2
        *  \return [ivg::Matrix] mapping functions (hydrostatic and wet part) 
        */
        ivg::Matrix get_mapping_functions( int idx );
     
        /**
        *  \b Description: \n
        *        Method to get epoch
        *  \param [in] no input parameters needed
        *  \return [ivg::Date] epoch
        */
        ivg::Date get_epoch()
        {
            return _epoch;
        };


        /**
        *  \b Description: \n
        *        Method to calculate delay and fill the designmatrix A
        *  \param [in]
        */
        void calc_delay( vector<double>::iterator design_iter, 
                         vector<double>::iterator apriori_iter, bool ion,
                         bool singleband,double *par_aplo =NULL);

	  
        /**
        *  \b Description: \n
        *        Method to calculate delay and fill the designmatrix A
        *  \param [in]
        */
        void calc_delay( vector<double>::iterator design_iter, 
                         vector<double>::iterator apriori_iter, bool ion = true,
                         char delaytype = 'g',double *par_aplo =NULL);
        /**
        *  \b Description: \n
        *        Method to calculate the total delay related to the selected delaymodel
        *  \param [in] delaymodel and some call-by-reference variables need in calc_delay(...)
        */
        double calc_total_delay(ivg::delaymodel model, ivg::Matrix &v_earth_ssb, ivg::Matrix &vel1, ivg::Matrix &vel2, ivg::Matrix &b_gcrs, ivg::Matrix &b_trs, ivg::Matrix &k);

        bool operator==( const Obs test ) const;

     
        /**
        *  \b Description: \n
        *        Method to get use_me
        */
        bool get_use_flag(){ return _use_me; };
        

        /**
        *  \b Description: \n
        *        Method to set use_me
        */
        void set_use_flag( bool flag ){ _use_me = flag; };

        
        double get_snr_bx() { return _snr_bx; };
        double get_snr_bs() { return _snr_bs; };
        double get_group_delay(){ return _group_delay; };
        double get_sigma_group_delay(){ return _sigma_delay; };
        double get_phase_delay(){ return _phase_delay; };
        double get_sigma_phase_delay(){ return _sigma_phase_delay; };
        /**
        *  \b Description: \n
        *        Method to translate between model defined as string in 
        *        the config-file and the enum ivg::delaymodel.
        */
        ivg::delaymodel get_delaymodel( string model );

    private:

        // class attributes

        // scan indexes
        int _sta1_scan_idx;
        int _sta2_scan_idx;

        // observation: two stations, one radio source, one epoch
        string _station1;
        string _station2;
        string _source;
        ivg::Date _epoch;
        double _mjd;
        
        // SNR
        double _snr_bx, _snr_bs;

        // (baseline-based) group delay and formal error
        double _group_delay;
        double _sigma_delay;
        
        // geocenter group delay (formal error is the same as for the baseline-based group delay)
        double _geoc_delay;
        
        // single band delay and formal error
        double _sb_delay;
        double _sigma_sb_delay;

        // phase delay and formal error
        double _phase_delay;
        double _sigma_phase_delay;
  
        // delay rate and formal error
        double _delay_rate;
        double _sigma_delay_rate;

        // delay ionososphere correction and formal error
        double _delta_ion;
        double _sigma_delta_ion;

        // delay rate ionosphere correation and formal error
        double _delta_ion_rate;
        double _sigma_delta_ion_rate;

        // pointer to scan and session object
        ivg::Scan* _scan;
        ivg::Session* _session;

        // observed minus computed
        double _oc;

        // azimuth and elevation
        ivg::Matrix _AzEl1;
        ivg::Matrix _AzEl2;

        // mapping function
        ivg::Matrix _mfs1;
        ivg::Matrix _mfs2;
        
        bool _use_me;
};



}

#endif // OBS_H
