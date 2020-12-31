#ifndef EOP_SERIES_H
#define EOP_SERIES_H

#include "logger.h"
#include "matrix.h"
#include "date.h"
#include <map>
#include <iostream>
#include <string>
#include "iers_wrapper.h"
#include "ivg_const.h"
#include "structs.h"
#include "parser.h"

extern "C"
{
#include <sofa.h>
}

/**
*
* @brief class Eop_series: EOPs are stored in radian
* @author bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

namespace ivg
{


// ===========================================================================
class Eop_series
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
        *  \return An instance of the class 'Eop_series'
        */
        Eop_series();
        
        
        /**
        *  \b Description: \n
        *    Constructor accumulating series of Eop_serieses. 
         *   IMPORTANT:
         *   MUST be already in ascending order
         *   NO overlapping of time
        *  \return An instance of the class 'Eop_series'
        */
        Eop_series(vector<ivg::Eop_series> series);
                
        /**
        *  \b Description: \n
        *        Constructor using type and matrix containing all EOPs.
        *  \param [in] [std::string] type
        *              [ivg::Matrix] ivg::Matrix with mjd, erp, erp_rates, nut, std_erp, std_erp_rates, std_nut
         *                           The data-matrix must have 17 cols!
        *  \return An instance of the class 'Eop_series'
        */
        Eop_series( const std::string &type, ivg::Matrix data );

        /**
        *  \b Description: \n
        *        Constructor using filename of EOP series and its type.
        *  \param [in] [std::string] file
        *              [std::string] type
        *              [ivg::Date] first epoch
        *              [ivg::Date] last epoch
        *  \return An instance of the class 'Eop_series'
        *  \todo Import of USNO finals and SOLVE modfiles
        */
        Eop_series( const std::string &file, const std::string &type,
                    const ivg::Date &start = ivg::Date( ivg::fake_mjd ),
                    const ivg::Date &end = ivg::Date( ivg::fake_mjd ) );

        /**
        *  s\b Description: \n
        *        Constructor using two parameters:
        *        start and end time tag
        *  \param [in] [ivg::Date] start
        *               [ivg::Date] end
        *  \return An instance of the class 'Eop_series'
        */
        // Eop_series( ivg::Date start, ivg::Date end );

        /**
        *  \b Description: \n
        *        Default deconstructor
        *  \param [in] no input parameters needed
        */
        ~Eop_series();


        // ==============================================
        // =========== OPERATORS: =======================
        // ==============================================
        Eop_series operator-( const Eop_series &other ) const;

        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        /**
        *  \b Description: \n
        *        Method to calculate mean pole with model coefficients
        *        from IERS Conventions 2010
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] [3 x 1] mean pole [rad]
        */
        ivg::Matrix calc_mean_pole( const ivg::Date &epoch ) const;

        /**
        *  \b Description: \n
        *        Method to calculate ERP (from C04 time series)
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] [3 x 1] ERPs [rad]
        */
        ivg::Matrix calc_erp( const ivg::Date &epoch );

        /**
        *  \b Description: \n
        *        Method to calculate ERP (from C04 time series)
        *  \param [in] [ivg::Date] epoch
        *              [std::string] interpolation type (default: linear)
        *              [bool] apply effect of zonal tides on UT1 before interpolation and restore it afterwards
        *  \return [ivg::Matrix] [3 x 1] ERPs [rad]
        */
        ivg::Matrix calc_erp( const ivg::Date &epoch, const std::string &type,
                              const bool &rg_zont ) const;

        /**
        *  \b Description: \n
        *        Method to calculate nutation corrections to IAU2000/2006 (from C04 time series)
        *  \param [in] [ivg::Date] epoch
        *  \return [ivg::Matrix] [2 x 1] nutation(X,Y) [rad]
        */
        ivg::Matrix calc_nut( const ivg::Date &epoch ) const;

        /**
        *  \b Description: \n
        *        Method to calculate nutation corrections to IAU2000/2006 (from C04 time series)
        *  \param [in] [ivg::Date] epoch
        *               [std::string] interpolation type
        *  \return [ivg::Matrix] [2 x 1] nutation(X,Y) [rad]
        */
        ivg::Matrix calc_nut( const ivg::Date &epoch, const std::string &int_type, const std::string &nut_type ) const;

        /**
        *  \b Description: \n
        *        Method to calculate wobble parameters
        *  \param [in] [ivg::Date] epoch
        *               [ivg::Matrix] [3 x 1] mean pole [rad]
        *               [ivg::Matrix] [3 x 1] ERPs [rad]
        *  \return [Matrix] [2 x 1] matrix with wobble parameters [rad]
        */
        ivg::Matrix calc_wobble_params( const ivg::Date &epoch,
                                        const ivg::Matrix &mean_pole,
                                        const ivg::Matrix &erp ) const;

        /**
        *  \b Description: \n
        *        Method to read EOP series from file. Supported types are: ("C04")
        *  \param [in] [std::string] c04_filename
        *              [ivg::Date] first epoch
        *              [ivg::Date] last epoch
        */
        void read_eop( const std::string &file, const std::string &type,
                       const ivg::Date &start, const ivg::Date &end );

        /**
        *  \b Description: \n
        *        Method to calculate subdaily Earth rotation parameters following models from IERS2010 Conventions. Effects of Earths tri-axiality can be turned off.
        *  \param [in] [ivg::Date] epoch
        *              [bool] ocean
        *              [bool] ut_libration
        *              [bool] pm_nutation
        *  \return [Matrix] [3 x 1] matrix with ERPs [rad]
        */
        ivg::Matrix calc_subdaily_erp_iers( const ivg::Date &epoch, const bool &ocean,const std::string &ocean_model,	
                                            const bool &ut_libration, const bool &pm_nutation ) const;
        /**
        *  \b Description: \n
        *        Method to calculate subdaily Earth rotation parameters following models from IERS2010 Conventions. Effects of Earths tri-axiality can be turned off.
        *  \param [in] [ivg::Date] epoch
        *  \return [Matrix] [3 x 1] matrix with ERPs [rad]
        */
        ivg::Matrix calc_subdaily_erp_iers( const ivg::Date &epoch ) const;

        /**
        *  \b Description: \n
        *        Method to calculate free core nutation following models from IERS2010 Conventions.
        *  \param [in] [ivg::Date] epoch
        *  \return [Matrix] [2 x 1] matrix with XY-nutation contributions [rad]
        */
        ivg::Matrix calc_free_core_nutation_iers( const ivg::Date &epoch ) const;

        /**
        *  \b Description: \n
        *        Method to build rotation matrix from CRF to TRF.
        *  \param [in] [ivg::Date] epoch
        *              [bool] nutation type
        *              [std::string] interpolation type
        *              [bool] apply zonal tide correction when interpolating UT1
        *              [bool] ocean
        *              [bool] ut_libration
        *              [bool] pm_nutation
        *  \return [Matrix] [3 x 3] CRF to TRF rotation matrix
        *  \todo test various nutation corrections when parameter estimation is sucessfull
        */
        ivg::Matrix form_crf2trf( const ivg::Date &epoch, const string &nut_type,
                                  const string &interp_type, const bool &rg_zont, const bool &ortho_eop, const std::string &hf_ocean_model,
                                  const bool &ut1libr, const bool &pmsdnut2, const bool derivative_flag=false,
                                  Partials_t2c * foo_ptr = NULL ) const;

        /**
        *  \b Description: \n
        *        Method to build rotation matrix from CRF to TRF based on EOPs.
        *  \param [in] [ivg::Date] epoch
        *  \return [Matrix] [3 x 3] CRF to TRF rotation matrix
        */
        ivg::Matrix form_crf2trf( const ivg::Date &epoch,
                                  const bool derivative_flag=false, Partials_t2c * foo_ptr = NULL ) const;

        /**
        *  \b Description: \n
        *        Method to build rotation matrix from CRF to TRF based on EOPs.
        *  \param [in] [ivg::Date] epoch
        *              [ivg::Matrix] ERP (X-/Y-pole UT1-TAI in rad)
         *             [ivg::Matrix] nuztation corretions to IAU2000/2006
        *  \return [Matrix] [3 x 3] CRF to TRF rotation matrix
        */
        ivg::Matrix form_crf2trf( const ivg::Date &epoch, const ivg::Matrix &pm_ut1,
                                  const ivg::Matrix &nut, const bool derivative_flag=false,
                                  Partials_t2c * foo_ptr = NULL ) const;

        /**
        *  \b Description: \n
        *        Method to build rotation matrix from CRF to TRF based on EOPs.
        *  \param [in] [ivg::Date] epoch
        *              [ivg::Matrix] ERP (X-/Y-pole, UT1-TAI in rad)
        *  \return [Matrix] [3 x 3] CRF to TRF rotation matrix
        */
        ivg::Matrix form_crf2trf( const ivg::Date &epoch, const ivg::Matrix &erp,
                                  bool derivative_flag=false, Partials_t2c * foo_ptr = NULL ) const;

        /**
        *  \b Description: \n
        *        Method to initialize Eop_series which has been created by default constructor.
        *  \param [in] [std::string] interpolation type
        *              [bool] apply zonal tide correction when interpolating UT1
        *              [bool] ocean
        *              [bool] ut_libration
        *              [bool] pm_nutation
        *              [bool] nutation type
        *              [int] pole tide version (2003, 2010 or 2015, 2015 = default)
        *  \return [Matrix] [3 x 3] CRF to TRF rotation matrix
        */
        void init( const string &intrp_mode, const bool &rg_zont,
                   const bool &hf_ocean,
		   const string &hf_ocean_model,
                   const bool &ut_libration, const bool &pm_nutation, const string nut_type, 
                   int pt_version=2015 );
        
        /**
        *  \b Description: \n
        *    Merge the data of two ivg::Eop_series. Settings like
        *    interpolation mode or handling of subdaily variations will be kept from 
        *    this.
        * \param [in] [ivg::Eop_series] new data to be added
        */
        void merge(const ivg::Eop_series &other);
        
        
        // experimental
        // better solution; mjd matrix in Eop_series for each EOP
        void replace(ivg::Eop_series &other, bool wob, bool ut1, bool nut);
        
        /**
        *  \b Description: \n
        *  Shows eop series data in form of a matrix.
        */
        void show(string out="");
        
        /**
        *  \b Description: \n
        *        Method to get the time series of a specific eop series.
        *        The time series contains [MJD, VALUE, STD] in sinex common units (mas, ms, masD, msD).
        *        Requesting time series with a single eop-string: "xpo", "ypo", "ut1" and order 0 or 1
        *  \param [in] [string] "xpo", "ypo" or "ut1"
        *  \param [in] [int] 0 or 1
        *  \return [ivg::Matrix] time series of requested eop [nx3]
        */
        ivg::Matrix get_time_series( string eop, int order );
        
        /**
        *  \b Description: \n
        *        Get all data in form of a matrix.Method to get the time series of a specific eop series.
        *        The time series contains [MJD, VALUE, STD] in sinex common units (mas, ms, masD, msD).
        *        Requesting time series with a single eop-string: "xpo", "ypo", "ut1" and order 0 or 1
        *  \return [ivg::Matrix] time series of all data with rows: mjd, erp, erprates, nut, erp_std, erprates_std, nut_std
        */
        ivg::Matrix get_data();
        
        /**
        *  \b Description: \n
        *        Checks whether eop series is initiailized or not. Based on _mjd matrix.
        *  \return [bool] in case of existing data (_mjd.rows > 0 ) == true, else false
        */
        bool is_initialized();
        
        /**
        *  \b Description: \n
        *        In case of accumlated Eop_series (special constructor), the database origin can requested
         *       where the mjd-corresponding EOPs come from; from which dbname.
        *  \return [string] corresponding dbname to mjd
        */
        string get_origin_db(double mjd);
        
        /**
        *  \b Description: \n
        *        Length of eop series.
        *  \return [int] number of rows from _mjd matrix
        */
        int length(){ return _mjd.rows(); };
        
    private:
        // ==============================================
        // ======== methods: ============================
        // ==============================================

        /**
        *  \b Description: \n
        *        Method to calculate the effect of zonal tides on UT1 for each entry of _mjd using IERS routine RG_ZONT.
        *        Values are stored in _dut_zonal
        */
        void _calc_dUT1_zonal_tides();

        template <size_t rows, size_t cols>
        ivg::Matrix _2darr2mat( const double (&array)[rows][cols] ) const;

        /**
        *  \b Description: \n
        *        Method to derive 3x3 matrix with partial derivatives of the rotation matrix wrt rotation angle.
        *  \param [in] [int] rotation axis (0=X, 1=Y, 2=Z)
        *              [double] rotation angle phi
        *  \return [Matrix] [3 x 3] drivative matrix dR/dphi
        */
        ivg::Matrix _3d_deriv_rotation_matrix( const int axis,
                                               const double phi ) const;


        // ==============================================
        // ======== class variables / attributes: =======
        // ==============================================

        std::string _type;

        // EOP data with daily resolution
        ivg::Matrix _mjd;
        ivg::Matrix _erp;           // X-pole, Y-pole, UT1-UTC
        ivg::Matrix _erp_rates;     
        ivg::Matrix _nut;           // differences to IAU2000/2006 in X and Y
        ivg::Matrix _std_erp;       // standard deviations
        ivg::Matrix _std_erp_rates;
        ivg::Matrix _std_nut;    
        
        map<double, string> _origin; // saves from which session (dbname) the eop-values come from

        string _intrp_mode;   // interpolation mode according to interpolation types in Matrix class
        bool _rg_zont;        // reduce effect of zonal tides befor interpolation of UT1
        bool _hf_ocean;       // account for oceanic high frequency variations (diurnal and semi-diurnal)
        string _hf_ocean_model;
        bool _ut_libration;   //             libration in UT1
        bool _pm_nutation;    //                       in PM
        string _nut_type;     // how to handle nutation (see comments above)

        ivg::Matrix _dut_zonal; // effect of zonal tides for each epoch

        int _pt_version;     // version for pole tide loading

        double _last_erp_epoch;
        ivg::Matrix _last_erp; 
};

} // # namespace ivg

#endif // EOP_SERIES_H
