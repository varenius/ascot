#ifndef TURBULENCE_H	// Include-Wachter
#define TURBULENCE_H

#include <iostream>
#include <string>
#include <iomanip>
#include <cmath> //-conf_hyperg
#include <iterator>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
#include <fstream>
#include <climits>
#include <fstream>
#include <random>
#include "matrix.h"
#include "obs.h"
#include "scan.h"
#include "analysis_station.h"
#include "ivg_const.h"
#include "date.h"
#include <ctime>
#include <time.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"    
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/bessel.hpp> 
#include <boost/math/special_functions/factorials.hpp>

//#include <omp.h>


using namespace std;

/**
*  \brief Class to investigate the turbulence behaviour and simulate turbulent parameters. \n
*         References: 
*         Treuhaft/Lanyi model and extensions: Treuhaft and Lanyi (1987); Nilsson and Haas (2010)
*         SIGMA-C model: Schoen and Brunner (2008a); Schoen and Brunner (2008b)
*         Matern Covariance model: Kermarrec and Schoen (2014)
*
*  \author    Sebastian Halsig
*  \version   v0.2
*  \date      March 2015
*/ 

/* 	======================================================================
			   class 'Turbulence' by S. Halsig 
				    HEADER-file
	======================================================================
*/

namespace ivg
{

/*
// struct for turbulence data
struct turbulence_data
{
   double Cn;     // refractive index structure constant [m^(-1/3)] 
   double h;      // effective height of the troposphere [m]
   ivg::Matrix vel;    // wind vector [m/s] --> [m/h]
   double zwd0;   // initial zenith wet delay [mm]
   double dh_seg; // time segments over which observations are to be correlated [h]
   double dh;     // height increment for numeric integration [m]
};
*/

// struct for observation data
struct observation_data
{
   double mjd;      // modified julian date [-]
   string source;   // source ID
   string sta1;     // name of first station 
   string sta2;     // name of second station
   double el1;      // elevation angle of first station [rad]
   double el2;      // elevation angle of second station [rad]
   double az1;      // azimuth of first station [rad]
   double az2;      // azimuth of second station [rad]
};



    class Turbulence
    {
	public:
	
	// ==============================================
	// =============== Konstruktoren: ===============
	// ==============================================

	/**
	 *  \b Description: \n
	 *        Default constructor
	 *  \param [in] no input parameters needed
	 *  \return An instance of the class 'Turbulence'
	 */
	Turbulence( );

        /**
         *  \b Description: \n
         *        Constructor using the full information 
         *        (station dependent and global turbulence parameters)
         *  \param [in] [std::map< std::string, turbulence_data >] 
         *                                   turb_sta: station dependent turbulence data \n
         *              [ivg::Matrix] spectral_coeffs: spectral coefficients \n
         *                                             (used for MATERN and SIGMA-C model)
         *  \return An instance of the class 'Turbulence'
         */
        Turbulence( std::map< std::string, turbulence_data > turb_sta, ivg::Matrix spectral_coeffs );

	/**
	 *  \b Description: \n
	 *        Constructor using two files and the station name
	 *  \param [in] [string] fnTurbData: filename for turbulence data \n
	 *               [string] fnAzEl: filename for azimuth and elevation angle \n
	 *               [string] name of an arbitrary station for turbulence investigations
	 *  \return An instance of the class 'Turbulence'
	 */
        Turbulence( std::vector<ivg::Analysis_station> stations, 
                    double Cn, double H, double v, ivg::Matrix vel, double v_dir );

	
	// ==============================================
	// =============== Destruktoren: ================
	// ==============================================

	/**
	 *  \b Description: \n
	 *        Default deconstructor
	 *  \param [in] no input parameters needed
	 */
	~Turbulence( );


	// ==============================================
	// ========== MEMBER-Funktionen: ================
	// ==============================================
	/**
	 *  \b Description: \n
	 *        This method generates a formatted output on console
	 *  \param [in] [int] nk: value for the decimal precision
	*/
	void show( int nk = 15 );
	 
	/**
	 *  \b Description: \n
	 *        This method sets the saturation length scale 
	 *  \param [in] [double] saturation length scale
	*/	
	void set_saturation_length( double L );

	/**
	 *  \b Description: \n
	 *        This method sets the outer scale length
	 *  \param [in] [double] outer scale length
	*/	
	void set_outer_scale_length( double L0 );

	/**
	 *  \b Description: \n
	 *        This method sets the station dependent turbulence parameter
         *        (i.e., one turbulence struct per station including the
         *        following parameters: structure constant, effective tropospheric
         *        height, wind vector and direction)
	 *  \param [in] [std::map< std::string, turbulence_data >] turb_sta: station dependent turbulence parameter
	*/	
        void set_turb_params( std::map< std::string, turbulence_data > turb_sta );

	/**
	 *  \b Description: \n
	 *        This method sets the spectral coefficients (a,b,c), which are used
         *        to consider anisotrpy (used in MATERN and SIGMA-C turbulence model)
	 *  \param [in] [ivg::Matrix] spectral coefficients
	*/	
        void set_spectral_coefficients( ivg::Matrix coeffs );

	/**
	 *  \b Description: \n
	 *        This method simulates the Equivalent Zenith Wet Delays (EZWD) \n
	 *        using initial Equivalent Zenith Wet Delays and a variance-covariance matrix
	 *  \param [in] [ivg::Matrix] VCM: variance covariance matrix
         *              [double] zwd0: apriori value 
	 *  \return [ivg::Matrix] simulated Equivalent Zenith Wet Delays (EZWD)
	*/
	ivg::Matrix simulate_ezwd( ivg::Matrix & VCM, double zwd0 );
			
	/**
	 *  \b Description: \n
	 *        This method calculates the Wet Delays in slant direction \n
	 *        using a mapping function
	 *  \param [in] [ivg::Matrix] simulated Equivalent Zenith Wet Delays (EZWD)
	 *  \return [ivg::Matrix] mapped Equivalent Zenith Wet Delays (EZWD)
	 *
	 *  \todo insert more complex and more accurate mapping function 
	*/	
	ivg::Matrix map_ezwd ( ivg::Matrix & EZWD ) ;	


	// ****************************************************************************
	// *********** (1) M A T E R N   C O V A R I A N C E   F A M I L Y ************
	// ****************************************************************************	
	/**
	 *  \b Description: \n
	 *        This method calculates a variance-covariance matrix due to refractivity 
         *        fluctuations using the Matern covariance family (see Kermarrec and
         *        Schoen, 2014 and Halsig et al., 2016).
	 *  \param [in]  [int] nobs: number of observations in the session \n
         *               [std::vector<ivg::Scan>] scans: vector of all scans
         *               [ivg::Trf*] trf: trf pointer
         *               [std::map< std::string, turbulence_data >] turb_sta: station dependent turbulence parameters
         *               [ivg::Matrix] spectral coefficients (a,b,c) to consider anisotropy 
	 *  \param [out] [ivg::Matrix] correlation matrix
	 *  \return [ivg::Matrix] Variance covariance matrix due to refractive fluctuations 	 
	*/
	ivg::Matrix calc_matern_vcm_model( int nobs, std::vector<ivg::Scan> scans, ivg::Trf * trf,
                                      ivg::Matrix spectral_coeffs, ivg::Matrix & C );

	/**
	 *  \b Description: \n
	 *        This method calculates a variance-covariance matrix due to refractivity 
         *        fluctuations using the Matern covariance family (see Kermarrec and
         *        Schoen, 2014 and Halsig et al., 2016). 
         *        This method works station-dependent, i.e. each station is taken into account
         *        individually, before a NxN session-wise variance-covariance 
         *        matrix is filled based on these station-dependent variance-covariance matrices.
	 *  \param [in]  [int] nobs: number of observations in the session \n
         *               [std::vector<ivg::Scan>] scans: vector of all scans
         *               [ivg::Trf*] trf: trf pointer
         *               [std::map< std::string, turbulence_data >] turb_sta: station dependent turbulence parameters
         *               [ivg::Matrix] spectral coefficients (a,b,c) to consider anisotropy 
	 *  \param [out] [ivg::Matrix] correlation matrix
	 *  \return [ivg::Matrix] Variance covariance matrix due to refractive fluctuations 	 
	*/        
	ivg::Matrix calc_station_matern_vcm_model( int nobs, std::vector<ivg::Scan> scans, ivg::Trf * trf,
                                                   ivg::Matrix spectral_coeffs, ivg::Matrix & C );        
        
	/**
	 *  \b Description: \n
	 *        This method calculates the separation distance between to two rays (for all observations)
	 *  \param [in] [std::string] sta1     - station 1 \n
	 *              [std::string] sta2     - station 2 \n
	 *              [ivg::Matrix] el       - elevation angles \n
	 *              [ivg::Matrix] az       - azimuth angles \n
	 *              [ivg::Matrix] xyz_sta1 - position of station 1 \n
	 *              [ivg::Matrix] xyz_sta2 - position of station 2 \n
	 *  \return [ivg::Matrix] separation distance between two rays
	*/		
	ivg::Matrix get_dist_between_rays( std::string sta1, std::string sta2, 
                                           ivg::Matrix el, ivg::Matrix az, 
                                           ivg::Matrix xyz_sta1, ivg::Matrix xyz_sta2 );

	/**
	 *  \b Description: \n
	 *        This method calculates the separation distance between to two rays 
         *        (for one observation) at height h (height h1 and h2 for 
         *        station 1 and 2, respectively).
         *        [this function is used in 'ivg::Matrix calc_matern_vcm_model( ... )']
         * 
	 *  \param [in] [double] el        - elevation angle \n
	 *              [double] az        - azimuth angle \n
	 *              [ivg::Matrix] sta1 - position of station 1 \n
	 *              [ivg::Matrix] sta2 - position of station 2 \n
	 *              [double] h1        - height h over station 1 \n
	 *              [double] h2        - height h over station 2 \n
	 *  \return [ivg::Matrix] separation distance between two rays
	*/		
        double get_dist_between_rays( double el, double az, ivg::Matrix sta1, ivg::Matrix sta2,
                                      double h1, double h2 );

        // test: matern covariance function
        ivg::Matrix calc_matern_covariance_function( int nobs, std::vector<ivg::Scan> scans, ivg::Trf *trf, 
                                                     ivg::Matrix spectral_coeffs, ivg::Matrix & C );

	// ****************************************************************************
	// ********** (2) B R U N N E R / S C H OE N   O R I G.   M O D E L ***********
	// ****************************************************************************
	/**
	 *  \b Description: \n
	 *        This method calculates a variance-covariance matrix due to refractive fluctuations
         *        according to Schoen and Brunner (2008a) and Schoen and Brunner (2008b). 
	 *  \param [in]  [double] dheight: height increment for numeric integration \n
	 *               [double] dh_Seg: time segments over which observations are to be correlated [h]
	 *  \param [out] [ivg::Matrix] reference: correlation matrix
	 *  \return [ivg::Matrix] Variance covariance matrix due to refractive fluctuations 	 
	*/
	ivg::Matrix calc_sigma_c_model( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, ivg::Matrix spectral_coeffs, double dh, ivg::Matrix & C ) ;	
	
        ivg::Matrix calc_sigma_c_model( ivg::Analysis_station* sta , ivg::Matrix azel, std::vector<ivg::Date> epochs, ivg::Matrix spectral_coeffs, double dh );
        
	/**
	 *  \b Description: \n
	 *        This method calculates the hypergeometric function
	 *  \param [in] [ivg::Matrix] a, b: n x 1 matrix (2x) \n
	 *              [double] z: value \n
	 *  \return [double] result of the hypergeometric function
	*/	
	double hypergeom( ivg::Matrix a, ivg::Matrix b, double z );
	
	
	// ****************************************************************************
	// ******************* (3) N I E L S O N  et al   M O D E L *******************
	// ****************************************************************************	
	/**
	 *  \b Description: \n
	 *        This method calculates a variance-covariance matrix due to
         *        atmospheric turbulence according to Nilsson and Haas (2010).
	 *  \param [in] [ivg::Matrix] structure function \n
	 *               [double] dheight: height increment for numeric integration \n
	 *               [double] dh_Seg: time segments over which observations are to be correlated [h]
	 *  \return [matrix] simulated Equivalent Zenith Wet Delays (EZWD)
	*/
	ivg::Matrix calc_vcm_onsala_model( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, double dheight, double dh_Seg, ivg::Matrix & C ) ;
        
        ivg::Matrix calc_vcm_onsala_model( ivg::Analysis_station* sta , ivg::Matrix azel, std::vector<ivg::Date> epochs, double dheight) ;

	/**
	 *  \b Description: \n
	 *        This method calculates a variance-covariance matrix due to
         *        atmospheric turbulence;
         *        modified version of the Nilsson and Haas (2010) model by J. Boehm  
         *        due to optimization purposes. 
	 *  \param [in] [ivg::Matrix] structure function \n
	 *               [double] dheight: height increment for numeric integration \n
	 *               [double] dh_Seg: time segments over which observations are to be correlated [h]
	 *  \return [matrix] simulated Equivalent Zenith Wet Delays (EZWD)
	*/        
        ivg::Matrix calc_vcm_vienna_model( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, double dheight, double dh_Seg, ivg::Matrix & C ) ;
        
        ivg::Matrix calc_vcm_vienna_model ( ivg::Analysis_station* sta , ivg::Matrix azel, std::vector<ivg::Date> epochs, double dh );
        
        ivg::Matrix calc_vcm_vienna_model ( const std::vector<double>& el, const std::vector<double>& az, const std::vector<double>& mjd,
                                                double dh, double dh_seg, double Cn, double h, ivg::Matrix vel );
        
        
	/**
	 *  \b Description: \n
	 *        This method calculates the structure function depending on
	 *        the structure constant and saturation legth scale
	 *  \param [in] no input parameters needed
	 *  \return [matrix] structure function
	*/
	ivg::Matrix calc_structure_function( int nobs ); 
        
        std::map< std::string, turbulence_data > set_param_form_cfg( const Setting& setup, ivg::Trf& trf );

	
	private:

	// ==============================================
	// ================= Methoden: ==================
	// ==============================================	

	/**
	 *  \b Description: \n
	 *        This method calculates the covariances between two observations i and j (given by its azimuth and elevation)
	 *  \param [in]  [double]      tau:               time difference between the two observations \n
	 *               [ivg::Matrix] az_el_i:           matrix contraining azimuth and elevation for observation i
	 *               [ivg::Matrix] az_el_j:           matrix contraining azimuth and elevation for observation j
         *               [std::string] staA/staB:         two stations of a baseline
         *               [ivg::Matrix] spectral_coeffs:   matrix containing the 3D spectral coefficients (a,b,c)
         *               [double]      dist_between_rays: distance between the two rays of observation i and j
	 *  \return [double] covariance between observation i and j
	*/
        double _calc_matern_covariance( double tau, ivg::Matrix az_el_i, ivg::Matrix az_el_j, 
                                        std::string staA, std::string staB, 
                                        ivg::Matrix spectral_coeffs, double dist_between_rays );  

	/**
	 *  \b Description: \n
	 *        This method calculates the covariances between two observations 
         *        i and j (given by its azimuth and elevation) for a certain station
	 *  \param [in]  [double]      tau:               time difference between the two observations \n
	 *               [ivg::Matrix] az_el_i:           matrix contraining azimuth and elevation for observation i
	 *               [ivg::Matrix] az_el_j:           matrix contraining azimuth and elevation for observation j
         *               [std::string] station:           station
         *               [ivg::Matrix] spectral_coeffs:   matrix containing the 3D spectral coefficients (a,b,c)
         *               [double]      dist_between_rays: distance between the two rays of observation i and j
	 *  \return [double] covariance between observation i and j for a certain station
	*/   
        double _calc_station_matern_covariance( double tau, ivg::Matrix az_el_i, ivg::Matrix az_el_j, 
                                                std::string station, 
                                                ivg::Matrix spectral_coeffs, double dist_between_rays );     
   
   
	/**
	 *  \b Description: \n
	 *        This method calculates the covariances between two observations i and j (given by its azimuth and elevation)
	 *  \param [in]  [double]      tau:               time difference between the two observations \n
	 *               [ivg::Matrix] az_el_i:           matrix contraining azimuth and elevation for observation i
	 *               [ivg::Matrix] az_el_j:           matrix contraining azimuth and elevation for observation j
         *               [ivg::Matrix] spectral_coeffs:   matrix containing the 3D spectral coefficients (a,b,c)
	 *  \return [double] covariance between observation i and j
	*/
        double _calc_sigma_c_covariance( double tau_A, double tau_B, ivg::Matrix az_el_i, ivg::Matrix az_el_j, double dh );
        double _calc_sigma_c_variance( ivg::Matrix az_el_i, ivg::Matrix spectral_coeffs );

        double _calc_onsala_covariance( double tau_i0, double tau_j0, double tau_ij, ivg::Matrix az_el_i, ivg::Matrix az_el_j, ivg::Matrix az_el_0, double dh );
        
        ivg::Matrix _azel2r (const ivg::Matrix& azel) const{ 
            ivg::Matrix r (3, 1, 1.0);
            double tanEl = tan( azel(1) );
            r(0) = cos( azel(0) ) / tanEl ;
            r(1) = sin( azel(0) ) / tanEl ;
            return r;
        };
        
        double _calc_onsala_structure_function (const ivg::Matrix& ri, const ivg::Matrix& rj, const ivg::Matrix& vd, double L) const{ 
            
            ivg::Matrix r = (ri - rj) + vd;

            double sp = ( r.transpose() * r )(0,0);
            double norm =  pow( sp, 1.0/3.0 );
            
            return norm / ( 1.0 + ( norm / pow( L, 2.0/3.0 ) ) );
        };
        

    // ==============================================
    // ================= Attribute: =================
    // ==============================================	

    // turbulence data ( effective tropo. height, structure constant, wind velocity,... )
    turbulence_data _turb_data;

    // STATION DEPENDENT TURBULENCE DATA (station IVS name is used as key): 
    // (a) structure constant, 
    // (b) effective tropospheric height, 
    // (c) wind velocity, 
    // (d) wind direction 
    std::map< std::string, turbulence_data  > _turb_sta;

    // GLOBAL TURBULENCE DATA: 
    // (a) saturation length scale
    double _L ;			

    // (b) outer scale length
    double _L0 ;
    
    // (c) 3D spectral coefficients (a,b,c) to consider anisotrophy
    ivg::Matrix _spectral_coeffs;
	
};

} // # namespace ivg

#endif // TURBULENCE_H










