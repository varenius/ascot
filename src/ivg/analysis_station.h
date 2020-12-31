#ifndef ANALYSISSTATIONS_H
#define ANALYSISSTATIONS_H

#include "station.h"
#include "matrix.h"
#include "date.h"
#include "eop_series.h"
#include "source.h"
#include "iers_wrapper.h"
#include <iterator>
#include "jpleph.h"
#include "ivg_const.h"
#include "tictoc.h"
#include "auxfunc.h"
#include "math.h"
#include "structs.h"
#include "logger.h"

extern "C"
{
#include <sofa.h>
}


/**
*
* @brief Analysis_station - derived class from Stations class:
*                           Station class for VLBI-analysis
* @author SH - bakkari developer team
* @date 2015-03-24
* @version 0.1
*/


namespace ivg
{

class Eop_series;  // forward declaration

// ===========================================================================
class Analysis_station : public Station
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
        *  \return An instance of the class 'Analysis_station'
        */
        Analysis_station();

        /**
        *  s\b Description: \n
        *        Constructor using four parameters:
        *        position and velocity matrix, reference epoch and names
        *  \param [in] [ivg::Matrix] position vector
        *               [ivg::Matrix] velocity vector
        *               [ivg::Matrix] reference epoch
        *               [std::map<ivg::staname, std::string>] station names
        *               (possible types: ivs_name, domes_no, cdp, lettercode, description)
        *  \return An instance of the class 'Analysis_station'
        */
        Analysis_station( ivg::Matrix pos0, ivg::Matrix vel0, ivg::Date refepoch,
                          std::vector<ivg::Date> discontinuity,
                          std::map<ivg::staname, std::string> names );
        
        Analysis_station( ivg::Matrix pos0, ivg::Matrix vel0,
                          ivg::Matrix pos0_std, ivg::Matrix vel0_std, 
                          ivg::Date refepoch,
                          std::vector<ivg::Date> discontinuity,
                          std::map<ivg::staname, std::string> names );

        /**
        *  \b Description: \n
        *        Default deconstructor
        *  \param [in] no input parameters needed
        */
        ~Analysis_station();


        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        // getter and setter
        /**
        *  \b Description: \n
        *        Method to set hydrology loading
        *  \param [in] [ivg::Matrix] nx4 (Epoch, UP, EASTWEST, NORTSOUTH)
        */
        void set_hydlo( ivg::Matrix hydlo );
        /**
        *  \b Description: \n
        *        Method to set non tidal atmospheric pressure loading
        *  \param [in] [ivg::Matrix] 4xn (Epoch, X, Y, Z) x (#datasets)
        */
        void set_nontidal_aplo( ivg::Matrix aplo );

        /**
        *  \b Description: \n
        *        Method to set tidal atmospheric pressure loading coeffs
        *  \param [in] [map<string,ivg::Wave>]
        */
        void set_tidal_aplo( map<string,ivg::Wave> waves );

        /**
        *  \b Description: \n
        *        Method to set ocean pole tide loading coefficients
        *  \param [in] [ivg::Matrix] ocean pole tide loading coefficients
        */
        void set_optl_coeff( ivg::Matrix optlc );

        /**
        *  \b Description: \n
        *        Method to set ocean loading coefficients
        *  \param [in] [std::vector<float>] ocean loading coefficients
        */
        void set_ocean_loading_coeff( std::vector<float> olc );

        /**
        *  \b Description: \n
        *        ITRF2014 specific CATREF psd-model coefficients
        *  \param [in] [std::vector<double>] ocean loading coefficients
        */
        void set_psd_coeff(std::vector<Psd> data){ _psd = data; };
        
        /**
        *  \b Description: \n
        *        Method to set eccentricity values
        *  \param [in] [ivg::Matrix] eccentricity values
        *              [std::vector<ivg::Date>] reference epoch
        */
        void set_eccentricity( ivg::Matrix ecc, std::vector<ivg::Date> refepoch );

        /**
        *  \b Description: \n
        *        Method to set gravitational deformation coefficients
        *  \param [in] [ivg::Matrix] gravdef coefficients
        *          
        */
        void set_gravdef( ivg::Matrix gd );

	/**
        *  \b Description: \n
        *        Method
        */
        void set_antenna_info(ivg::Antenna antenna_info);
        
        /**
        *  \b Description: \n
        *        Method return a reference in order to be able to set variables 
        *        of antenna_info
        */
        ivg::Antenna & set_antenna_info();
        
        /**
        *  \b Description: \n
        *        Method return a reference in order to be able to set/get variables 
        *        within Channels-struct, in different parser functions (loif,tracks,freq,hdpos)
        */
        ivg::Channels & set_channel_setup();

        /**
        *  \b Description: \n
        *        Method to set VMF1 data
        *  \param [in] [ivg::Matrix] vmf1 data
        */
        void set_vmf1_data( ivg::Matrix vmf1_data );

        /**
        *  \b Description: \n
        *        Method to set VMF3 data
        *  \param [in] [ivg::Matrix] vmf3 data
        */
        void set_vmf3_data( ivg::Matrix vmf3_data );
  
        /**
        *  \b Description: \n
        *        Method to set external hydrological meteorological data
        *  \param [in] [ivg::Matrix] external hydrological meteorological data
        */
        void set_zhd_data( ivg::Matrix zhd_data );

        /**
        *  \b Description: \n
        *        Method to set external raytraced delays
        *        derived from RADIATE raytracing algorithm (TU Vienna)
        *  \param [in] [std::map< ivg::Date,std::map<std::string,ivg::Matrix> >] external raytraced delays
        */
        void set_raytracing_data( std::map< ivg::Date,std::map<std::string,ivg::Matrix> > rt_data );

        /**
        *  \b Description: \n
        *        Method to get external raytraced delays
        *        derived from RADIATE raytracing algorithm (TU Vienna)
        *  \return [std::map< ivg::Date,std::map<std::string,ivg::Matrix> >] external raytraced delays
        */
        std::map< ivg::Date,std::map<std::string,ivg::Matrix> > get_raytracing_data(){ return _raytracing_data; };
        
        /**
        *  \b Description: \n
        *        Method to get external hydrological meteorological data
        *  \return [ivg::Matrix] external hydrological meteorological data
        */
        ivg::Matrix get_zhd_data( ){ return _zhd_data; };        
        
        /**
        *  \b Description: \n
        *        Method to get the antenna information, based on antenna.cat
        */
        ivg::Antenna get_antenna_info();
        
        /**
        *  \b Description: \n
        *        Method to get the equipment information, based on equip.cat
        */
        ivg::Equip get_equip_info();

        /**
        *  \b Description: \n
        *        Method to get the eccentricity, based on ECCDAT.ecc
        */

        ivg::Matrix get_eccentricity( ivg::Date epoch );
        
        /**
        *  \b Description: \n
        *        Method to add eccentricity values
        *  \param [in] [ivg::Matrix] eccentricity values
        *              [std::vector<ivg::Date>] reference epoch
        */
        void add_eccentricity( ivg::Matrix ecc, ivg::Date epoch, string epoch_type );

        /**
        *  \b Description: \n
        *        Method to add complete equipment information including SEFD values for different bands
        *  \param [in] [ivg::Equip] equip, containing SEFD[band]
        */
        void add_equip_info( ivg::Equip equip);

        /**
        *  \b Description: \n
        *        Method to add frequency sequence for a corresponding band
        *  \param [in] [ivg::band] band type, like X or S
        *              [ivg::Matrix] column-matrix containing channel frequencies
        */
        void add_band_frequency_sequence( ivg::band band, ivg::Matrix freq_seq );
        
        /**
        *  \b Description: \n
        *        Method to add mask information. Matrix contains azimuth and elevation in two columns
        *  \param [in] [ivg::band] ivg::Matrix column-matrix [nx2]
        */
        void add_mask_info( ivg::Matrix mask );
        
        
        // output functions
        /**
        *  \b Description: \n
        *        Method to show station
        *  \param [in] no input parameters needed
        */
        void show();


        // position and baseline functions
        /**
        *  \b Description: \n
        *        Method to add a discontinuity to an existing station instance
        *  \param [in] pos, vel, ref_epoch, dicontinuity
        *  \return void
        */
        void add_discontinuity( const ivg::Matrix pos, const ivg::Matrix vel,
                                const ivg::Date ref_epoch, const ivg::Date discontinuity,
                                ivg::Matrix pos_std = ivg::Matrix(3,1,0.0),
                                ivg::Matrix vel_std = ivg::Matrix(3,1,0.0) );
        
        /**
        *  \b Description: \n
        *        Method to set a complete set of discontinuities to an existing station instance
        *  \param [in] refepoch, pos0, vel0, vector of dates (discontinuties)
        *  \return void
        */
        void set_discontinuity( const ivg::Date refepoch, const ivg::Matrix pos0,
        const ivg::Matrix vel0, const vector<ivg::Date> discons);

        /**
        *  \b Description: \n
        *        Method to detect weather there are discontinuities or not
        *  \param [in] no input parameters needed
        *  \return [bool] true = discontinuities; false = no discontinuities
        */
        bool has_discontinuity();
        
        /**
        *  \b Description: \n
        *        Method to select correct XYZ vector in case of discontinuties
        *  \param [in] [ivg::Date] epoch of interest
        *  \return [Matrix] [3 x 1] position vector
        */
        ivg::Matrix get_xyz( ivg::Date epoch );
        
        /**
        *  \b Description: \n
        *        Method to select correct XYZ velocity vector in case of discontinuties
        *  \param [in] [ivg::Date] epoch of interest
        *  \return [Matrix] [3 x 1] velocitiy vector
        */
        ivg::Matrix get_vel( ivg::Date epoch );

        /**
        *  \b Description: \n
        *        Method to calculate position vector
        *  \param [in] [ivg::Date] epoch of interest
        *  \return [Matrix] [3 x 1] position vector
        */
        ivg::Matrix calc_xyz( ivg::Date epoch );

        /**
        *  \b Description: \n
        *        Method to calculate position vector
        *  \param [in] [ivg::Date] epoch of interest
        *  \return [ivg::Matrix] [3x1] position vector
        *          [ivg::Matrix] [3x1] standard deviations
        */        
        ivg::Matrix calc_xyz( ivg::Date epoch, ivg::Matrix & std );

        /**
        *  \b Description: \n
        *        Method to calculate position vector including displacements
        *  \param [in] [ivg::Date] epoch of interest
        *              [std::vector<string>] strings for defining displacements that are applied: ("SOLID EARTH TIDES", "POLE TIDE", "OCEAN LOADING", "OCEAN POLE TIDE LOADING"       )
        *              [ivg::Eop_series *] pointer to initialized Eop_series instance
        *              [void pointer] Ephemeris
        *  \return [Matrix] [3 x 1] position vector
        */
        ivg::Matrix calc_xyz( ivg::Date epoch, std::vector<string> disp, void * ephem = NULL,
                              ivg::Eop_series * eops = NULL );

         /**
        *  \b Description: \n
        *        Method to calculate position vector including displacements
        *  \param [in] [ivg::Date] epoch of interest
        *              [std::vector<string>] strings for defining displacements that are applied: ("SOLID EARTH TIDES", "POLE TIDE", "OCEAN LOADING", "OCEAN POLE TIDE LOADING"       )
        *              [ivg::Eop_series *] pointer to initialized Eop_series instance
        *              [void pointer] Ephemeris
        *              [ivg:Matrix] trf2crf
        *  \return [Matrix] [3 x 1] position vector
        */
        ivg::Matrix calc_xyz( ivg::Date epoch, std::vector<string> disp, void * ephem,
                              ivg::Eop_series * eops , ivg::Matrix trf2crf );

        /**
        *  \b Description: \n
        *        Method d => delta
        */
        ivg::Matrix calc_dxyz( ivg::Date epoch, std::vector<string> disp,
                               void * ephem, ivg::Eop_series * eops = NULL);

        /**
        *  \b Description: \n
        *        Method d => delta
        */
        ivg::Matrix calc_dxyz( ivg::Date epoch, std::vector<string> disp,
                               void * ephem, ivg::Eop_series * eops, ivg::Matrix trf2crf );
        /**
        *  \b Description: \n
        *        Method d => delta
        */
        ivg::Matrix calc_dren( ivg::Date epoch, std::vector<string> disp,
                               void * ephem, ivg::Eop_series * eops = NULL );

        /**
        *  \b Description: \n
        *        Method to adjust the zenith SEFDs for the source elevation.
        *        If no SEFD parameters are given for this station, the zenith SEFD is used.
        *  \param [in] [ivg::band] X-/S- or whatever-band
        *              [ivg::Source *] observed source
        *              [ivg::Date *] epoch of interest
        *  \return [double] slant SEFD
        */
        double calc_slant_sefd(ivg::band band, ivg::Source * src_ptr,ivg::Date * epo_ptr);

        
        ivg::Matrix compute_az_el_sun( ivg::Date epoch, void* ephem, ivg::Matrix c2t );

        
        /**
        *  \b Description: \n
        *        Method to move this-station for a delta_xyz
        *  \param [in] [ivg::Matrix] d_xyz
        */
        void move_station( ivg::Matrix d_xyz );
        
        
        /**
        *  \b Description: \n
        *        Method to calculate the baseline length between two
        *        stations (this station and another station)
        *  \param [in] [ivg::Analysis_station] other: other station
        *               [ivg::Matrix] epoch: epoch
        *  \return [double] baseline length
        */
        double calc_baseline_length( ivg::Analysis_station other, ivg::Date epoch, double & std ) ;

        /**
        *  \b Description: \n
        *        Method to calc latitude (lat), longitude (lon) of reference epoch
        *        and height (h) of a station
        *  \param [in] no input parameters needed
        *  \return [Matrix] [n x 3] matrix with lat-lon-h
        */
        Matrix calc_lat_lon_h( string ref_ell = "GRS80" );

        /**
        *  \b Description: \n
        *        Method to calc latitude (lat), longitude (lon) of current epoch
        *        and height (h) of a station
        *  \param [in] [ivg::Date] epoch
        *  \return [Matrix] [n x 3] matrix with lat-lon-h
        */
        Matrix calc_lat_lon_h( const ivg::Date epoch,  string ref_ell = "GRS80" );


        /**
        *  \b Description: \n
        *        Method to get elevation angle (el) and azimuth (az)
        *  \param [in] [std::string] station name: ivs_name
        *  \return [Matrix] [n x 2] matrix with el-az
        */
        ivg::Matrix calc_az_el( ivg::Date epoch, ivg::Source src );
        

        ivg::Matrix calc_az_el( ivg::Date epoch, ivg::Matrix src_vec,
                                ivg::Matrix crf2trf );
        
        ivg::Matrix calc_az_el( const ivg::Matrix& src_vec_trf );
        
   
        // rotate vector k given in az el (radian) into trf (station coordiantes xyz)
        ivg::Matrix azel2k( double az, double el, ivg::Matrix xyz );
        
        // rotate vector k given in az el (radian) into trf
        ivg::Matrix azel2k( double az, double el, ivg::Date epoch );
        
        /**
        *  \b Description: \n
        *        Method to check if a point defined by azimut and elevation is visibile
        *        at the site. Based on the elevation-mask and telescope limitations.
        *  \param [in] [double] azimut and elevation of the point to be checked in rad
        *  \return [bool] true/false 
        */
        bool check_visibility( double azi, double ele );
        /**
        *  \b Description: \n
        *        Method to calualte delay due to gravitational deformation.
        *  \param [in] [double] elevation angle in rad
        *  \return [double] graviational deformation delay in ps
        */
        double calc_gravdef( double el);
        /**
        *  \b Description: \n
        *        Method to get the elevation mask of a site at a given azimut.
        *  \param [in] [double] azimut in rad
        *  \return [double] elevation angle in rad
        */
        double get_ele_mask( double azi );
        
        /**
        *  \b Description: \n
        *        Method calculates wrap (e.g. between 270 and 810) and azimut angle to turn
        *        from an old position "old_wrap" to a new azimut. Gives feedback if wrap might leads
        *        to issues concerning lower/upper wrap edge and 180° turns 
        *  \param [in] [double] old wrap, new_azi
        *  \param [out] [double] new wrap and azimut rotation angle
        *   \return [bool] suitability info concerning wrap-tolerance, TRUE or FALSE
        */
        bool calc_wrap_and_rotangle( double old_wrap, double new_azi, double &new_wrap, double &d_azi);
        
        /**
        *  \b Description: \n
        *        Method to determine wrap zone (-,C,W) based on a wrap.
        *  \param [in] [double] wrap
        *  \param [out] [string] wrap zone
        */
        string determine_wrap_zone(double wrap);
        
        /**
        *  \b Description: \n
        *        Method to calc displacements due to pole tide
        *  \param [in] [ivg::Date] epoch
         *             [ivg::Eop_series *] pointer to initialized Eop_series instance
        *  \return [ivg::Matrix] displacements due to pole tide
        */
        ivg::Matrix calc_pole_tide( const ivg::Date epoch,
                                    ivg::Eop_series * eops = NULL );

        /**
        *  \b Description: \n
        *        Method to calc displacements due to ocean pole tide loading
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] displacements due to ocean pole tide loading
        */
        ivg::Matrix calc_optl( const ivg::Date epoch, ivg::Matrix op_coeff_rne,
                               ivg::Eop_series * eops = NULL ) ;

        /**
        *  \b Description: \n
        *        Method to calc displacements due to ocean pole tide loading
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] displacements due to ocean pole tide loading
        */
        ivg::Matrix calc_optl( const ivg::Date epoch,
                               ivg::Eop_series * eops = NULL ) ;

        /**
        *  \b Description: \n
        *        Method to calc displacements due to solid Earth tides
        *  \param [in] [ivg::Date] epoch
        *              [void pointer] Ephemeris
        *  \return [ivg::Matrix] displacements due to solid Earth tides
        */
        ivg::Matrix calc_solid_earth_tides( ivg::Date epoch, void *ephem,
                                            ivg::Matrix c2t );


        /**
        *  \b Description: \n
        *        Method to calc displacements due to hydrology
        *  \param [in] [ivg::Date] epoch
         *             [string interpolation_type] e.g. "cspline" oer "linear"
        *  \return [ivg::Matrix] displacements due to ocean pole tide loading
        */
        ivg::Matrix calc_hydlo( ivg::Date epoch, const string interpolation_type   );


        /**
        *  \b Description: \n
        *        Method to calc displacements due to ocean loading
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] displacements due to ocean loading
        */
        ivg::Matrix calc_ocean_loading( ivg::Date epoch );

        /**
        *  \b Description: \n
        *        Method to calc displacements due to non tidal aplos
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] displacement
        */
        ivg::Matrix calc_nontidal_aplo( ivg::Date epoch,
                                        const string interpolation_type  );

        /**
        *  \b Description: \n
        *        Method to calc displacements due to tidal aplos
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] displacements
        */
        ivg::Matrix calc_tidal_aplo( ivg::Date epoch );

        ivg::Matrix interpolate_ext_met_data( std::string type, ivg::Date epoch,
                                              std::string interpolation_type );

        std::string get_data_status(string key){return _data_status[key];};

                /**
        *  \b Description: \n
        * Method to calc displacements due parametric post-seismic 
        * deformation model (CATREF psd) introduced with ITRF2014
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] displacements (X,Y,Z) [m]
        */
        ivg::Matrix calc_psd_displacement( ivg::Date epo );
        
        /**
        *  \b Description: \n
        *        Method to build rotation matrix from tpoocentric (up, east north)
        *        to geocentric system (XYZ).
        *  \param [in] no input parameters needed
        *  \return [ivg::Matrix] rotation matrix topo2geo (3x1)
        */
        ivg::Matrix form_topo2geo();
        
        ivg::Matrix calc_pos_diff( ivg::Analysis_station other );
        
        double calc_bandwidth_rms(ivg::band band);
        double calc_bandwidth_rms(ivg::Matrix freq);


        double calc_mean_frequency(ivg::band band){return _freq_sequences[band].meanD();};
        
        int get_nfreq(ivg::band band){return _freq_sequences[band].length();};
        
        int get_idx() const {return _idx;}; 
        void set_idx( int idx) { _idx = idx;};
        /**
        *  \b Description: \n
        *        Returns the ah and aw values from GPT2. 
        *        
        *  \param [in] none
        *  \return doube ah, double aw
        */ 
        void get_gpt2_ah_aw(double *ah, double *aw);
        /**
        *  \b Description: \n
        *        Sets the ah and aw values from GPT2. 
        *        
        *  \param [in] double ah,aw
        *  \return none
        */ 
        void set_gpt2_ah_aw(double ah, double aw);

        /**
        *  \b Description: \n
        *        Returns the ah and aw values from GPT3. 
        *        
        *  \param [in] none
        *  \return doube ah, double aw
        */ 
        void get_gpt3_ah_aw(double *ah, double *aw);
        /**
        *  \b Description: \n
        *        Sets the ah and aw values from GPT3. 
        *        
        *  \param [in] double ah,aw
        *  \return none
        */ 
        void set_gpt3_ah_aw(double ah, double aw);

        int get_num_obs() {return _num_obs;};
        ivg::Date get_first_epoch() {return _first_epoch;};
        ivg::Date get_last_epoch() {return _last_epoch;};
        void inc_num_obs() {_num_obs++;};
        void dec_num_obs() {_num_obs--;};
        void set_num_obs(int n) {_num_obs=n;};
        void set_first_epoch(ivg::Date ep){_first_epoch=ep;};
        void set_last_epoch(ivg::Date ep){_last_epoch=ep;};
  
    private:
        /**
        *  \b Description: \n
        *        convert radial, east north displacements to XYZ
        *  \param [in] [ivg::Matrix] dUEN (3x1)
        *  \return [ivg::Matrix] dXYZ (3x1)
        */
        ivg::Matrix _ren2xyz( ivg::Matrix uen , string type = "ren");


        /**
        *  \b Description: \n
        *        convert radial, east north displacements to XYZ
        *  \param [in] [ivg::Matrix] dUEN (3x1)
        *  \return [ivg::Matrix] dXYZ (3x1)
        */
        ivg::Matrix _xyz2ren( ivg::Matrix uen );


        /**
        *  \b Description: \n
        *        Method
        */
        ivg::Matrix _calc_lat_lon_h( ivg::Matrix xyz,
                                     const std::string ref_ell = "GRS80" );

        // ==============================================
        // ==============================================

        // Hydrology Loading
        ivg::Matrix _hydlo; // (nx4)

        // Ocean Pole Tide Loading Coefficients
        ivg::Matrix _optl_coeff;

        // Atmospheric Pressure Loading (non tidal)
        ivg::Matrix _aplo_xyz; //(4xn)

        // Atmospheric Pressure Loading (tidal)
        map<string,ivg::Wave> _aplo_tidal;

        // Ocean Loading Coefficients
        vector<float> _ol_coeff;

        // eccentricity values
        ivg::Matrix _ecc;
        vector<ivg::Date> _ecc_epochs;

        // gravitational deformations
        ivg::Matrix _grav_deform;
        //antenna infos
        ivg::Antenna _antenna_info;
        
        //station loif,rx,rec,tracks,freq information for scheduling
        ivg::Channels _channel_setup;

        // VMF1 data (VMF1 mapping function coeff., hydrostatic delay from ECMWF,...)
        ivg::Matrix _vmf1_data;

        // VMF3 data (VMF3 mapping function coeff., hydrostatic delay from ECMWF,...)
        ivg::Matrix _vmf3_data;

        // external hydrological (ZHD) meteorological data
        ivg::Matrix _zhd_data;
        
        std::map< ivg::Date,std::map<std::string,ivg::Matrix> > _raytracing_data;
        
        // Saves latest crf2trf transformation matrix for solid earth tides (performance!)
        map<ivg::Date, ivg::Matrix> _last_crf2trf;
        
        // saves information about set data, e.g. if HYDLO is set or ECC etc.
        map<string,string> _data_status;
        
        // saves information about post-seismic deformation (for Iers::parametric_)
        std::vector<Psd> _psd;
        
        // save sefd information for each band from equip.cat
        ivg::Equip _equip;
        map< ivg::band, vector<double> > _sefd;
        
        // save frequency sequence for each band information from freq.cat
        map< ivg::band, ivg::Matrix > _freq_sequences;
        
        // save mask elevation and azimut information from mask.cat
        // ivg::Matrix used instead of a map for easy interpolation
        ivg::Matrix _mask; // stores values in rad
        
        // index in TRF vector. Indices have to be created first with member function in TRF
        int _idx = -1;

        //  GPT2 a-vlues
        double _gpt2_ah=0, _gpt2_aw=0;
        double _gpt3_ah=0, _gpt3_aw=0;

        int _num_obs=0;
        ivg::Date _first_epoch;
        ivg::Date _last_epoch;
};


} // # namespace ivg

#endif // ANALYSISSTATIONS_H
