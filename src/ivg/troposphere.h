#ifndef TROPOSPHERE_H
#define TROPOSPHERE_H

#include "matrix.h"
#include "date.h"
#include "analysis_station.h"
#include "iers_wrapper.h"
#include "ivg_const.h"
#include <iterator>
#include "math.h"
#include "tictoc.h"
#include <sstream>
#include <string>

/**
*
* @brief Troposphere class
* @author SH - bakkari developer team
* @date 2015-03-24
* @version 0.1
*/


namespace ivg
{
    enum bendingmode{NO,APPROX,INSITU,HEIGHT};

// ===========================================================================
class Troposphere
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
        *  \return An instance of the class 'Troposphere'
        */
        Troposphere();

        /**
        *  \b Description: \n
        *        Constructor using two input parameters (station pointer, epoch)
        *  \param [in] [ivg::Analysis_station*] station pointer
        *              [ivg::Date] epoch
        *  \return An instance of the class 'Troposphere'
        */
        Troposphere( ivg::Analysis_station * station, ivg::Date epoch );

        /**
        *  \b Description: \n
        *        Default deconstructor
        *  \param [in] no input parameters needed
        */
        ~Troposphere();


        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        // getter and setter
        /**
        *  \b Description: \n
        *        Method to show observation information
        *  \param [in] no input parameters needed
        */
        void show_met_data();

        /**
        *  \b Description: \n
        *        Method to get temperature [deg C]
        *  \param [in] no input parameters needed
        *  \return [double] temperature [deg C]
        */
        double get_temperature() const;

        /**
        *  \b Description: \n
        *        Method to get pressure [hPa]
        *  \param [in] no input parameters needed
        *  \return [double] pressure [hPa]
        */
        double get_pressure() const;

        /**
        *  \b Description: \n
        *        Method to get humidity [%]
        *  \param [in] no input parameters needed
        *  \return [double] humidity [%]
        */
        double get_humidity() const;

        /**
        *  \b Description: \n
        *        Method to get Zenith Wet Delay (ZWD) [m]
        *  \param [in] no input parameters needed
        *  \return [double] Zenith Wet Delay (ZWD) [m]
        */
        double get_zwd() const;

        /**
        *  \b Description: \n
        *        Method to get Zenith Hydrostatic Delay (ZHD) [m]
        *  \param [in] no input parameters needed
        *  \return [double] Zenith Hydrostatic Delay (ZHD) [m]
        */
        double get_zhd() const;

        /**
        *  \b Description: \n
        *        Method to get Zenith Total Delay (ZTD) [m]
        *  \param [in] no input parameters needed
        *  \return [double] Zenith Total Delay (ZTD) [m]
        */
        double get_ztd() const;
        
         /**
        *  \b Description: \n
        *        Method to get Slant Total Delay (STD) [m]
        *  \param [in] no input parameters needed
        *  \return [double] Slant Total Delay (STD) [m]
        */
        double get_std() const {return _std;};

        /**
        *  \b Description: \n
        *        Method to get mapping function coeffizients
        *        for the hydrostatic and wet part
        *  \param [in] no input parameters needed
        *  \param [out] [double] hydrostatic mapping function coefficient
        *  \param [out] [double] wet mapping function coefficient
        */
        void get_mapping_function( double &mfh, double &mfw );

        /**
        *  \b Description: \n
        *        Method to get delay due to asymmetric gradients
        *  \param [out] [double] north gradient
        *  \param [out] [double] east gradient
        */
        double get_gradient_delay( double &grN, double &grE );

        /**
        *  \b Description: \n
        *        Method to get the North-South gradient 
        *  \return [double] North-South gradient
        */
        double get_north_gradient() const { return _grN; };        
        
        /**
        *  \b Description: \n
        *        Method to get the East-West gradient 
        *  \return [double] East-West gradient
        */
        double get_east_gradient() const { return _grE; };        
        
        /**
        *  \b Description: \n
        *        Method to set temperature [deg C]
        *  \param [in] [double] temperature [deg C]
        */
        void set_temperature( double t );

        /**
        *  \b Description: \n
        *        Method to set pressure [hPa]
        *  \param [in] [double] pressure [hPa]
        */
        void set_pressure( double p );

        /**
        *  \b Description: \n
        *        Method to set humidity [%]
        *  \param [in] [double] humidity [%]
        */
        void set_humidity( double h );

        /**
        *  \b Description: \n
        *        Method to set meteorological data:
        *         - temperature [deg C]
        *         - pressure [hPa]
        *         - humidity [%]
        *  \param [in] [double] temperature [deg C]
        *              [double] pressure [hPa]
        *              [double] humidity [%]
        */
        void set_meteorology( double t, double p, double h, int code );

        /**
        *  \b Description: \n
        *        Method to set station pointer
        *  \param [in] [ivg::Analysis_station*] analysis_station pointer
        */
        void set_station( ivg::Analysis_station * station );

        /**
        *  \b Description: \n
        *        Method to set epoch
        *  \param [in] [ivg::Date] epoch
        */
        void set_epoch( ivg::Date epo );

        /**
        *  \b Description: \n
        *        Method to calc zenith hydrostatic delay (ZHD) [m]
        *        according to the modified Saastamoinen model
        *        (see Davis, 1985)
        *  \param [in] no input parameters needed
        */
        void set_zhd( double zhd ){ _zhd = zhd; };        

        /**
        *  \b Description: \n
        *        Method to set zenith wet delay [m]
        *  \param [in] [double] zenith wet delay [m]
        */
        void set_zwd( double zwd ){ _zwd = zwd; };    
        void add_zwd( double zwd );    

        /**
        *  \b Description: \n
        *        Method to set the east component of the azimuthal gradients
        *  \param [in] [double] east component of the azimuthal gradients
        */
        void set_egr( double egr ){ _grE = egr; };  
        void add_egr( double egr );  

        /**
        *  \b Description: \n
        *        Method to set the north component of the azimuthal gradients
        *  \param [in] [double] east component of the azimuthal gradients
        */
        void set_ngr( double ngr ){ _grN = ngr; };
        void add_ngr( double ngr );
        
        /**
        *  \b Description: \n
        *        Method to set the wet mapping factor
        *  \param [in] [double] wet mapping factor
        */
        void set_wet_mapping_factor( double mfw ){ _mfw = mfw; };         
        
        // calcuate tropospheric delays
        /**
        *  \b Description: \n
        *        Method to calc zenith hydrostatic delay (ZHD) [s]
        *        according to the modified Saastamoinen model
        *        (see Davis, 1985)
        *  \param [in] no input parameters needed
        */
        double calc_zhd_davis( );
        static double calc_zhd_davis(double pressure, double phi, double h );

        /**
        *  \b Description: \n
        *        Method to calc zenith total delay (ZTD) [m]
        *  \param [in] no input parameters needed
        */
        void calc_ztd( );

        /**
        *  \b Description: \n
        *        Method to determine the asymmetric delay d in meters caused
        *        by gradients. The north and east gradients should be used
        *        with the gradient model by Chen and Herring (1997).
        *  \param [in] [std::string] mf - mapping function type
        *                 "vmf1" for Vienna Mapping Functions 1
        *                 "gmf"  for Global Mapping Functions
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [in] [double] az - azimuth from north [rad]
        */
        void calc_tropo_delay( std::string mf, double el, double az );


        // mapping functions
        /**
        *  \b Description: \n
        *        Methode to determine the Vienna Mapping Function 1
        *        (VMF1, site dependent version). The coefficients can be obtained
        *        from the website http://ggosatm.hg.tuwien.ac.at/DELAY/SITE/
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [out] [double] hydrostatic vienna mapping function coefficient
        *  \param [out] [double] wet mapping vienna function coefficient
        */
        void calc_vmf1( double el, std::string interpolation_type, double &vmf1h,
                        double &vmf1w );

        // mapping functions
        /**
        *  \b Description: \n
        *        Methode to determine the Vienna Mapping Function 3
        *        (VMF3, site dependent version). The coefficients can be obtained
        *        from the website http://vmf.geo.tuwien.ac.at/
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [out] [double] hydrostatic vienna mapping function coefficient
        *  \param [out] [double] wet mapping vienna function coefficient
        */
        void calc_vmf3( double el, std::string interpolation_type, double &vmf3h,
                        double &vmf3w );

        /**
        *  \b Description: \n
        *        Method to determine the Global Mapping Functions GMF (Boehm et al. 2006)
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [out] [double] hydrostatic global mapping function coefficient
        *  \param [out] [double] wet mapping global function coefficient
        */
        void calc_gmf( double el, double &gmfh, double &gmfw );

        /**
        *  \b Description: \n
        *        Methode to determine the Vienna Mapping Function 1
        *        (VMF1, site dependent version). The coefficients are
        *        obtained using the GPT2 model
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [out] [double] hydrostatic vienna mapping function coefficient
        *  \param [out] [double] wet mapping vienna function coefficient
        */
        void calc_vmf1_gpt2( double el, double &vmf1h, double &vmf1w );


        /**
        *  \b Description: \n
        *        Methode to determine the Vienna Mapping Function 3
        *        (VMF1, site dependent version). The coefficients are
        *        obtained using the GPT3 model
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [out] [double] hydrostatic vienna mapping function coefficient
        *  \param [out] [double] wet mapping vienna function coefficient
        */
        void calc_vmf3_gpt3( double el, double &vmf3h, double &vmf3w );

        // azimutal gradients
        /**
        *  \b Description: \n
        *        Method to determine the gradient mapping function with the
        *        gradient model by Chen and Herring (1997)
        *  \param [in] [double] el - elevation angle [rad]
        */
        void calc_gradient_mf( double el );

        /**
        *  \b Description: \n
        *        Method to determine the asymmetric delay d in meters caused
        *        by gradients. The north and east gradients should be used
        *        with the gradient model by Chen and Herring (1997).
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [in] [double] az - azimuth from north [rad]
        */
        void calc_gradients( double el, double az,
                             double &delay, double &_grN, double &_grE,
			     std::string grad_type="apg", std::string grad_name="");

        /**
        *  \b Description: \n
        *        Method to determine the asymmetric delay d [m] caused
        *        by gradients using the gradient model by Chen and Herring (1997).
        *  \param [in] [double] el - elevation angle [rad]
        *  \param [in] [double] az - azimuth from north [rad]
        *  \param [in] [double] grN - North-South gradient
        *  \param [in] [double] grE - East-West gradient
        *
        *  \return [double] gradient delay [m]  
        */
        double calc_gradient_delay( double el, double az );
        
        
        // blind models
        /**
        *  \b Description: \n
        *        Method to determine pressure, temperature, temperature lapse rate, water
        *        vapour pressure, hydrostatic and wet mapping function coefficients ah and aw.
        *  \param [out] [double *] press - Pressure given in hPa
        *  \param [out] [double *] temp  - Temperature in degrees Celsius
        *  \param [out] [double *] dtemp - Temperature lapse rate in degrees per km
        *  \param [out] [double *] e     - Water vapour pressure in hPa
        *  \param [out] [double *] ah    - hydrostatic mapping function coefficient at zero height (VMF1)
        *  \param [out] [double *] aw    - wet mapping function coefficient (VMF1)
        */
        void use_gpt2( double &pres, double &temp, double &dtemp,
                       double &e, double &ah, double &aw );


        // correction terms (troposphere ties,...)
        /**
        *  \b Description: \n
        *        Method to calculate troposphere ties [mm] (cf. Teke et al. 2011, p8)
        *  \param [in] [double] station height
        *              [double] reference station height
        *              [double] latitude of reference station
        *              [double] pressure of reference station [hPa]
        *              [double] temperature of reference station [deg C]
        *              [double] humidity of reference station [%]
        */
        static void calc_tropospheric_ties( double H, double H0, double phi0,
                                     double p0, double t0, double e0 );

        static double calc_tropospheric_ties_zwd( double H, double H0,
                                     double p0, double t0, double e0 );
        
        /**
        *  \b Description: \n
        *        Method to calculate tropospheric bending angle (according to CALC11 $MK5_ROOT/progs/calc11/caxom.f).
        *  \param [in] [double] elevation [rad]
        *              [double] station height [m]
        *              [bendingmode] APPROX|INSITU|HEIGHT 
        *  \param [out] [double *] bending angle [rad]
        */
        double calc_bending_angle(double elevation,double sitheight,bendingmode type);

        /**
        *  \b Description: \n
        *        Method to incorporate tropospheric bending angle into source vector.
        *  \param [in] [double] source vector in TRF
        *              [Matrix] Azimuth/Elevation [rad]
        *              [bendingmode] APPROX|INSITU|HEIGHT 
        *  \param [out] [double *] corrected source vector in TRF [rad]
        */
        ivg::Matrix build_apparent_source_vector(const ivg::Matrix &src,const ivg::Matrix &azel,
                                          double sitheight,bendingmode type);
        /**
        *  \b Description: \n
        *        Method to absolute path to gpt2-grid file needed in GPT2.F
        *  \param [in] [string] grd-file 
        */
        void set_gpt2_grdfilename(string fname);
        /**
        *  \b Description: \n
        *        Method to absolute path to gpt3-grid file needed in GPT3_1.F90
        *  \param [in] [string] grd-file 
        */
        void set_gpt3_grdfilename(string fname);

    private:
        // statement function for calculation of tropospheric bending angle
        double _delta_fun(double ad1,double ad2,double bd1,double bd2,double zd2);

        // ==============================================
        // ==============================================

        // station and epoch
        ivg::Analysis_station * _sta;
        ivg::Date _epoch;

        // meteorology
        double _temperature;
        double _pressure;
        double _humidity;
        int _humCode;

        // tropospheric delay [s]
        double _zwd;
        double _zhd;
        double _ztd;
        double _std;

        // mapping function
        double _mfh;
        double _mfw;

        // azimutal gradients (north and east component, gradient mapping function)
        double _grN;
        double _grE;
        double _grad_del;
        double _mg;

        // turbulence parameter
        double _Cn;
        double _H;
        ivg::Matrix v;
        double _L0;
        double _a;
        double _b;
        double _c;
        double _lambda;

        string _gpt2_grd_file;
        string _gpt3_grd_file;

};


} // # namespace ivg

#endif // TROPOSPHERE_H
