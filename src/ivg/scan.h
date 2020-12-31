#ifndef SCAN_H
#define SCAN_H

#include "logger.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "obs.h"
#include "date.h"
#include "station.h"
#include "analysis_station.h"
#include "auxfunc.h"
#include "source.h"
#include "troposphere.h"
#include "matrix.h"
#include "ivg_const.h"
#include "structs.h"
#include "parser.h"
#include "eop_series.h"


/**
*
* @brief Session class
* @author SH - bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

using namespace std;


namespace ivg
{

class Session; // forward declaration



// ===========================================================================
class Scan
// ===========================================================================
{
        friend class Obs;

    public:

        struct Sta_epo
        {
            ivg::Analysis_station* sta_ptr;
            ivg::Troposphere tropo;
            double cable_cal;
            double feed_rotation;
            double clo0;
        };

        // ==============================================
        // =============== Constructors: ================
        // ==============================================

        /**
        *  \b Description: \n
        *        Default constructor
        *  \param [in] no input parameters needed
        *  \return An instance of the class 'Scan'
        */
        Scan();

        /**
        *  s\b Description: \n
        *        Constructor using two parameters:
        *        two station names, source name, epoch, mjd
        *  \param [in] [ivg::Date] epoch
        *              [std::string] source name
        *              [ivg::Matrix] TRF2CRF rotation matrix
         *          
        *  \return An instance of the class 'Scan'
        */
        Scan( ivg::Date epoch, ivg::Source* src, const ivg::Matrix t2c, Partials_t2c * deriv_ptr );

        /**
        *  \b Description: \n
        *        Copy constructor
        *  \param [in] [ivg::Scan] other observation
        *  \return An instance of the class 'Scan'
        */
        Scan( const ivg::Scan &other );
        
        /**
         *  \b Description: \n
         *        assignment operator '=' using another scan
         *  \param [in] [Session] &other: other scan
         */
        Scan & operator=( const Scan &other );

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
        *        Method to get epoch
        *  \param [in] no input parameters needed
        *  \return [ivg::Date] epoch
        */
        ivg::Date get_epoch();

        /**
        *  \b Description: \n
        *        Method to get source name
        *  \param [in] no input parameters needed
        *  \return [std::string] source name
        */
        ivg::Source* get_source();

        /**
        *  \b Description: \n
        *        Method to set epoch
        *  \param [in] [ivg::Date] epoch
        */
        void set_epoch( ivg::Date epo );

        /**
        *  \b Description: \n
        *        Method to set source name
        *  \param [in] [std::string] source name
        */
        void set_source( ivg::Source* src );

        /**
        *  \b Description: \n
        *        Method to set rotation matrix from TRF to CRF
        *  \param [in] [ivg::Matrix] trf2crf matrix
        */
        void set_trf2crf_matrix( ivg::Matrix t2c, Partials_t2c * deriv_ptr = NULL );

        /**
        *  \b Description: \n
        *        Method to check if station is already in station list
        *  \param [in] [ivg::Analysis_station*] pointer to station
        *  \return bool \n return true if station is already in station list
        *                  return false if station is not in station list
        */
        bool includes_station(ivg::Analysis_station *station );

        /**
        *  \b Description: \n
        *        Method to add new observation to observation list
        *  \param [in] [ivg::Obs] new observation
        *  \return ivg::Obs* returns pointer to added observation within scan
        */
        ivg::Obs* add_obs( ivg::Obs new_obs );

        /**
        *  \b Description: \n
        *        Method to remove an observation from the observation list
        *  \param [in] [int] Index
        */
        void rem_obs( int idx ){ _observations.erase( _observations.begin()+idx ); };


        /**
        *  \b Description: \n
        *        Method to add new analysis station to station list
        *  \param [in] [ivg::Analysis_station*] new analysis station
        */
        int add_sta( ivg::Analysis_station * new_sta );

        /**
        *  \b Description: \n
        *        Method to add new troposphere to troposphere list
        *  \param [in] [ivg::Troposphere] new troposphere
        */
        void add_tropo( ivg::Troposphere new_tropo );

        /**
        *  \b Description: \n
        *        Method to add either meteorological data (data_type = meteo) 
        *        or external zenith hydrostatic delays (data_type = zhd). 
        *        In case of the meteorological data,
        *        the the modified Saastamoinen model (Davis et al., 1985) is 
        *        used to calculate the zenith hydrostatic delays. In case of
        *        directly using external data, zenith hydrostatic delays from
        *        raytracing (e.g., from RADIATE, TU Vienna) or numerical weather
        *        models (e.g., ECMWF) are used.
        *  \param [in] [double] atmospheric temperature at site (deg C)
        *              [double] atmospheric pressure at site (hPa)
        *              [double] atmospheric humidity at site (hPa)
        *              [int] Humidity parameter definition code:
        * 			0 = humidity parameter is relative humidity (%)
        * 			1 = humidity parameter is dew point (deg. C)
        * 			2 = humidity parameter is wet bulb temperature (deg. C)
        *              [std::string] data_type: use meteorological data or external ZHDs
        *              [std::string] met_data: selecte data source (e.g., ECMWF, insitu, GPT2, ...)
        *              [std::string] grdname: absolute path to gpt2-grid file needed in GPT2.F
        *              [std::string] grdname3: absolute path to gpt3-grid file needed in GPT3_1.F90
        */
        void add_scan_meteorology( int idx, double t, double p,
                                   double h = -999.99, int code = 0,
                                   std::string data_type = "meteo",
                                   std::string met_data = "database",
                                   std::string grdname = "",
				   std::string grdname3 = "",
				   std::string grad_type = "apg",
				   std::string grad_file = "");

        /**
        *  \b Description: \n
        *        Method to set the hydrostatic delay in zenith direction
        *  \param [in] [int] scan index
        *              [double] zenith hydrostatic delay [m]
        */
        void set_scan_zhd( int idx, double zhd ){ _data.at( idx ).tropo.set_zhd( zhd ); };      

        /**
        *  \b Description: \n
        *        Method to set the wet delay in zenith direction
        *  \param [in] [int] scan index
        *              [double] zenith wet delay [m]
        */
        void set_scan_zwd( int idx, double zwd ){ _data.at( idx ).tropo.add_zwd( zwd ); };      
        
        /**
        *  \b Description: \n
        *        Method to set the east component of azimuthal gradients
        *  \param [in] [int] scan index
        *              [double] east component of azimuthal gradients
        */
        void set_scan_egr( int idx, double egr ){ _data.at( idx ).tropo.add_egr( egr ); };      

        /**
        *  \b Description: \n
        *        Method to set the north component of azimuthal gradients
        *  \param [in] [int] scan index
        *              [double] north component of azimuthal gradients
        */
        void set_scan_ngr( int idx, double ngr ){ _data.at( idx ).tropo.add_ngr( ngr ); };

        /**
        *  \b Description: \n
        *        Method to set the clock offset parameter per scan and station
        *  \param [in] [int] scan index
        *              [double] clock offset parameter
        */
        void set_scan_clo( int idx, double cl0 ){ _data.at( idx ).clo0 = cl0; };
        
        /**
        *  \b Description: \n
        *        Method to add cable calibration correltion for both sites
        *  \param [in] [double] cable calibration correltion
        */
        void add_cable_cal( double cc );

        /**
        *  \b Description: \n
        *        Method to get the pointer to the observation at index 'idx'
        *        in the observation vector
        *  \param [in] no input parameter needed
        *  \return [ivg::Obs*] observation pointer to observation at index 'idx'
        */
        ivg::Obs* get_obs_ptr( int idx )
        {
            return &_observations.at( idx );
        };

        /**
        *  \b Description: \n
        *        Method to get the data (ivg::Analysis_station* sta_ptr;
        *                                ivg::Troposphere tropo;
        *                                double cable_cal)
        *        of the 'Sta_epo' struct at index 'idx'
        *  \param [in] no input parameter needed
        *  \return [ivg::Obs*] observation pointer to observation at index 'idx'
        */
        Sta_epo get_data( int idx )
        {
            return _data.at( idx );
        };

        std::vector<Sta_epo> get_data( )
        {
            return _data;
        };      
        
        /**
        *  \b Description: \n
        *        Method to get the pointer to the station at index 'idx'
        *        in the analysis station vector
        *  \param [in] no input parameter needed
        *  \return [ivg::Analysis_station*] analysis station pointer to the station at index 'idx'
        */
        ivg::Analysis_station* get_sta_ptr( int idx )
        {
            return _data.at( idx ).sta_ptr;
        };

        /**
        *  \b Description: \n
        *        Method to get the rotation matrix from TRF to CRF
        *  \param [in] no input parameter needed
        *  \return [ivg::Matrix] rotation matrix from TRF to CRF
        */
        ivg::Matrix get_trf2crf_matrix()
        {
            return _trf2crf;
        };

        /**
        *  \b Description: \n
        *        Method to get the number of observations
        *  \param [in] no input parameter needed
        *  \return [int] number of observations
        */
        int get_nobs()
        {
            return _observations.size();
        };
        
        /**
        *  \b Description: \n
        *        Method to return all contributing stations within scan.
        *  \param [in] no input parameter needed
        *  \return [vector<ivg::Analysis_station*>] vector of station pointers
        */
        vector<ivg::Analysis_station*> get_stations(){ return _stations; };

        /**
        *  \b Description: \n
        *        Method to calculate the additional delay due to axis offsets
        *        including correction for tropospheric bending
        *  \param [in] [int] scan idx of station 1 or 2
        *              [ivg::Matrix] apparent source vector (including abberation and troposherical bending) in TRF
        *  \return [double] delay due to axis offset
        */
        double calc_axis_offset_delay( int idx, ivg::Matrix apparent_src_vec_trf );
        double calc_axis_offset_delay_vievs( int idx, ivg::Matrix abb_src_vec);

        /**
        *  \b Description: \n
        *        Method to calculate the additional delay due to
        *        gravitational deformation of the antenna
        *  \param [in] [int] scan idx of station 1 or 2
        *              [ivg::Matrix] apparent source vector (including abberation and troposherical bending) in TRF
        *  \return [double] delay due to gravitational deforamtion
        */
        double calc_gravdef_delay( int idx, ivg::Matrix apparent_src_vec_trf);
        /**
        *  \b Description: \n
        *        Method to calculate the additional delay due to
        *        thermal expansion of the antenna
        *  \param [in] [int] scan idx of station 1 or 2
        *              [ivg::Matrix] apparent source vector (including abberation and troposherical bending) in TRF
        *  \return [double] delay due to thermal expansion
        */
        double calc_thermal_expansion_delay( int idx, ivg::Matrix apparent_src_vec_trf );

        /**
        *  \b Description: \n
        *        Method to calculate the mean (group delay) sigma of all observations within a scan. This needs to be done for scheduling
        *        and is based on jSked procedure. Maybe better criteria necessary?!
        *  \param [in]  none
        * 
        *  \return [double] mean sigma of all observations (using _sigma_delay from obs.h)
        */
        double mean_group_delay_sigma();
        
        /**
        *  \b Description: \n
        *        Method to calculate the greates slew time of any of the participating telescopes to a new source.
        *  \param [in] [ivg::Source*] next source for that the slew time should be calculated
        *         [in] [ivg::Eop_series*] eop_series for proper crf2trf
        *         [in/out] [string] information string containing slew times and wraps for each station 
        * 
        *  \return [double] limiting slewtime
        */
        double calc_slew_time(ivg::Source* src, ivg::Eop_series* eops, string &infostr);
        
        /**
        *  \b Description: \n
        *        Method to set the cable wrap for a specific station within the scan
        *  \param [in] [ivg::Analysis_station*] corresponding station
        *         [in] [double] wrap to be stored
        */
        void set_cable_wrap(ivg::Analysis_station* ac, double wrap){ _cable_wrap[ac]=wrap; };
        
        /**
        *  \b Description: \n
        *        Method to get the cable wrap for a specific station
        *  \param [in] [ivg::Analysis_station*] station of desire
        * 
        *  \return [double] cable wrap for a specific station
        */
        double get_cable_wrap(ivg::Analysis_station* ac){ return _cable_wrap[ac]; };
        
        void set_cable_wrapzone(ivg::Analysis_station* ac, string wrap){ _wrap_zone[ac]=wrap; };
        string get_cable_wrapzone(ivg::Analysis_station* ac) { return _wrap_zone.at(ac); };
        
        // stores the last telescope aiming position (az/el) at the end of the scan
        void set_aiming(ivg::Analysis_station* ac, ivg::Matrix aiming){ _azel_aiming[ac]=aiming; };
        map<ivg::Analysis_station*, ivg::Matrix> get_aiming(){ return _azel_aiming; };
        
        /**
        *  \b Description: \n
        *        Method to set the scheduled duration within a scan.
        *        Only stores the duration if the new duration is bigger than the old one. 
        *  \param [in] [double] duration in seconds
        */
        void set_schedulded_duration(double duration)
        { 
            if(duration > _duration)
                _duration = duration; 
        };
        
        /**
        *  \b Description: \n
        *        Method to get the longest scheduled duration within a scan of all observations.
        * 
        *  \return [double] scheduled duration
        */
        double get_scheduled_duration(){return _duration;};
        
        // to-do: still in use?

//	bool operator==( const Scan test ) const;
//	int find_obs( const Obs test ) const;



    private:
        // class attributes
        ivg::Date _epoch;
        ivg::Source* _source;

        // observation list
        std::vector<ivg::Obs> _observations;

        // list of station pointer
        std::vector<ivg::Analysis_station*> _stations;
        
        // troposphere list
        std::vector<ivg::Troposphere> _tropo;

        // list of cable calibration data
        std::vector<double> _cable_cal;

        vector<Sta_epo> _data;

        // rotation matrix from TRF to CRF
        ivg::Matrix _trf2crf;
        Partials_t2c _partials_t2c;

        // stores information about cablewrap for each station
        map<ivg::Analysis_station*, double> _cable_wrap;
        
        // stores information about cablewrap zone for each station
        map<ivg::Analysis_station*, string> _wrap_zone;
        
        // aiming position (az/el) of telescopes at the end of the scan
        // used for skyplots
        map<ivg::Analysis_station*, ivg::Matrix> _azel_aiming;
        
        // stores duration schedulded within ivg::ASCOT schedule part
        double _duration;

//  Session* _session;

//  friend function
//	friend void Session::add_scan( Scan* in );

};

}

#endif // SCAN_H
