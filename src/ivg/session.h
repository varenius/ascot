#ifndef SESSION_H
#define SESSION_H

#include "logger.h"
#include "obs.h"
#include "scan.h"
#include "trf.h"
#include "crf.h"
#include "date.h"
#include "auxfunc.h"
#include "structs.h"
#include "ls_solution.h"
#include "lsa.h"
#include "ls_neq.h"
#include "icls.h"
#include "troposphere.h"
#include "turbulence.h"
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "tictoc.h"
#include "param_list.h"
#include "lsc.h"
#include "lpSked/transits.h"
#include "masterfile.h"

#include <cstdlib>
#include <libconfig.h++>

/**
*
* @brief Session class
* @author SH - bakkari developer team
* @date 2015-03-24
* @version 0.1
*
*/
using namespace libconfig;
using namespace std;

// each session got resdiuals after analysis
enum residtype{ all, station, baseline, source };

// obsdelaytype (group-delay or single-band-delay)
// MAXDELAYTYPE is used as stop condition in for loops
enum delaytype{ group, single, MAXDELAYTYPE };

struct Residual
{
    string name; // e.g. GILCREEK or GILCREEK-WETTZELL
    double wrms;
    int nobs;
    vector<int> idx_in;     // index of inliers in data
    vector<int> idx_out;    // index of outliers in data
    vector<int> idx_origin; // index in observation vector
    vector<int> idx_ncfile;  // index in ncfile;
    ivg::Matrix data;
    residtype type; // e.g. residtype::station or residtype::baseline
    vector<string> first_station, source_names;
};   

namespace ivg
{  
    
// ===========================================================================
 class Session
// ===========================================================================
{
        friend class Obs;
        friend class Session_inout;
        friend class Simulation;
        friend class Schedule;

    public:

        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Default constructor
        *  \param [in] no input parameters needed
        *  \return An instance of the class 'Session'
        */
        Session( );
        
        /**
        *  \b Description: \n
        *        Copy constructor
        *  \param [in] [ivg::Scan] other observation
        *  \return An instance of the class 'Scan'
        */
        Session( const ivg::Session &other );
        
        /**
         *  \b Description: \n
         *        assignment operator '=' using another session
         *  \param [in] [Session] &other: other session
         */
        Session & operator=( const Session &other );
        
         /**
         *  \b Description: \n
         *        operator for adding a session + assignment to this-session '+='
         *  \param [in] [Session] S: other session
         */
        Session & operator+=( Session &other );

        /**
        *  \b Description: \n
        *        ONE and ONLY conctructor for session. Includes the initialization of the Eop_series.
        *        All other member variables (e.g. trf, crf) need to be initiailzed with session_inout.load() 
        *  \param [in] [std::string] session name
        *  \return An instance of the class 'Session'
        */
        Session( Setting *setup, string name, void **ephem = NULL, int session_cnt = 0);
        
        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Method to initialize session from vgosdb or ngs
        */
        void init_vgosdb_ngs_solution();
        
        /**
        *  \b Description: \n
        *        Method to initialize session from sinex. We need to define what to to with 
         *       param_list, TRF, CRF, EOPs
        *  \param [in]  adjustment_options string 
        *  \return
        */
        void init_snx_solution(string adjustment_options, ivg::Date t_0);
        
        void init_session_eop( const Setting& eopblock, ivg::Date date );
        
        ivg::Eop_series init_eop( const Setting& eopblock, ivg::Date date, double days = 15.0 );
       
        void init_param_list(  );
        
        void change_eop_series( const Setting& eopblock, ivg::Date date ) ;
        
        /**
        *  \b Description: \n
        *        Method we don't want
        */
        void reduce_and_constrain();
        
        /**
        *  \b Description: \n
        *        Method we don't want
        */
        void modify_parameterization();
        
        /**
        *  \b Description: \n
        *        Method to get session name
        *  \param [in] no input parameters needed
        *  \return [std::string] session name
        */
        string get_name();
        
        /**
        *  \b Description: \n
        *        Method to get the setup-part of the config-file
        *  \param [in] no input parameters needed
        *  \return [Setting *] libconfig Setting Pointer
        */
        Setting *get_setup(){return _setup;};
        
        /**
        *  \b Description: \n
        *        Method to get the mid epoch of the session
        *  \return [ivg::Date] mid epoch of session
        */
        ivg::Date get_mid_epoch();
        
        /**
        *  \b Description: \n
        *        Method to get the observation epochs for specific stations or
        *        the epochs of all observations
        *  \param [in] [std::string] station (if "ALL" selected; return 
        *                            the epochs of all observations)
        *  \return [ivg::Matrix] observation epochs 
        */        
        ivg::Matrix get_obs_epochs( std::string station = "ALL" );
        

        /**
        *  \b Description: \n
        *        Method to get TRF pointer
        *  \param [in] no input parameters needed
        *  \return [ivg::Trf*] TRF pointer
        */
        ivg::Trf *get_trf_ptr(){ return &_trf; };

        /**
        *  \b Description: \n
        *        Method to get CRF pointer
        *  \param [in] no input parameters needed
        *  \return [ivg::Crf*] CRF pointer
        */
        ivg::Crf *get_crf_ptr(){ return &_crf; };
        /**
        *  \b Description: \n
        *        Method to get parameter list
        *  \param [in] no input parameters needed
        *  \return [ivg::Param_list] parameter list
        */
        ivg::Param_list *get_param_list_ptr(){ return &_param_list; };
        
        ivg::Ls_solution *get_solution_ptr(){ return _solution; };
        
        std::vector<ivg::Scan> *get_scan_ptr(){ return &_scans; };
        
        std::string get_session_path() const { return _session_path; };
        
        ivg::Turbulence*  getTurbulencePtr(){ return &_turbulence;};
         
        /**
        *  \b Description: \n
        *        Method to show observation information
        *  \param [in] no input parameters needed
        */
        void show();

        /**
        *  \b Description: \n
        *        Method to check wheather a test session is equal to this session
        *  \param [in] [ivg::Session] test session
        *  \return [bool] true: same session
        *                 false: other session
        */
        bool operator==( const Session test ) const;

        /**
        *  \b Description: \n
        *        Method to manage
        *           ... the modification of the parametrization
        *               according to the config file,
        *           ... to create the constraints and NNR/NNT equations,
        *           ... and to build and solve the normal equation system
        *  \param [in] [bool] if true and 'LSA' is used only observations corresponding to
        *           the reference clock are used to calculate parameters. However, 
        *           residuals are computed for all observations
        *  \param [out] returning session info string for logfile 
        */
        string solve( bool use_only_indep_bl = false);
        
        /**
        *  \b Description: \n
        *        Method to check if the estimates in Up,East,North of each station don't
        *        exceed the threshold in configfile (vertical and horizontal)
        *  \param [out] string vector containing affected stations 
        */
        vector<string> check_stations_estimates();

        /**
        *  \b Description: \n
        *        Method to write output files for the Calc/Solve C++ backend:
        *           obs_info:   observation file (incl. epochs, stations, sources, elevation and azimuth)
        *           param:      list of parameters (incl. index number)
        *           apriori:    a priori values for the parameters
        *           ref_epoch:  link between epochs and indexes
        *           A, wgt, oc: design matrix, weight vector or matrix and o-c vector (binary)
        *  \param [in] [std::string] direction, where the output files should be stored
        */
        void write_backend_files( std::string dir );

        double calc_dist_between_rays( ivg::Obs* obs_ptr );

        void create_solution_info();

        stringstream *get_file_comment(){return &_file_comment; };
        
        map< string, ivg::Matrix > get_obsstats() const { return _obsstats; };
                
        map< delaytype , vector<Residual> > get_residuals() const {return _residuals; };
        void set_residuals(  map< delaytype , vector<Residual> > residuals ){ _residuals = residuals; };
        
        int get_nobs_orig(){ return _nobs_orig; };
               
        vector<int> *get_origin_obs_idxs(){ return &_origin_obs_idxs; };
        
        vector<int> get_obs_idxs_in_ncfile();
        
        void setAmbigRes(bool _ambigRes){this->_ambigRes = _ambigRes;};
        
        bool getAmbigRes(){return _ambigRes ;};
        
        double get_max_ambiguity_spacing() const{ return _ambiguity_spacing.max(); };
        ivg::Matrix get_ambiguity_spacing() const { return _ambiguity_spacing; };
        
        ivg::Matrix get_eff_freq() const { return _eff_freq; };
                
        void set_band_type(ivg::band type){_band_type = type;};
        
        ivg::Matrix* get_num_ambig_ptr() {return &_num_ambig;};
        
        /**
        *  \b Description: \n
        *        returns the row indices of the design matrix that do not 
        *        correspond to the given station
        *        the station 
        *  \param [in] [std::string] station
        * \return [vector<int>] indices of observations not corresponding to a
        *          certain station
        */
        vector<int> get_star_formation_indices(std::string station) ;    
        
        const lps::Transits& get_trasintsVisibleFromStation() const{ return _trasintsVisibleFromStation;};
        const lps::Transits& get_trasintsVisibleFrom2Stations() const{ return _trasintsVisibleFrom2Stations;};
        const lps::Transits& get_complete_transits() const{ return _complete_transits;}; // all visible transits
        const lps::Transits& get_common_transits() const{ return _common_transits; }; // transits visbile from all stations
        ivg::Date getEnd() const { return _end;};
        ivg::Date getStart() const {return _start;};
        
        void set_end(ivg::Date end) {_end = end;};
        void set_start(ivg::Date start) {_start = start;};
        
        ivg::Eop_series& get_eops() { return _eops;};
        ivg::Schedule* get_schedule_ptr() const {return _schedule_ptr;};
            
        void obsCounterPlus(unsigned int n){_nobs+=n; _nobs_orig+=n;};
        
         // in case dt = 0 dt is read from configfile
        void calc_transits(  unsigned int dt = 0);
        
        void create_crf_trf_inidices(){
            _crf.create_source_indices();
            _trf.create_station_indices();
        }
        
        string coverage2string( std::map<int, std::vector<double> > coverage, unsigned precision = 1 ,  bool percent=true);
        void show_coverage( std::map<int, std::vector<double> > coverage );
        
        void find_and_mark_unused_sources();
        
        void cable_wrap_time_series( std::map<std::string, ivg::Matrix>& wrap,
                                     std::map<std::string, ivg::Matrix>& time);
        
        void intSked2Latex(ivg::Masterfile masterfile, std::map<int, std::vector<double> >& ca,
                                   std::map<int, std::vector<double> >& cr, std::map<int, std::vector<double> >& caa,
                                   std::map<int, std::vector<double> >& car, double objfun);
        
        const void* const get_ephem() const {return _ephem;};
        void*  get_ephem() {return _ephem;};
                
        /**
        *  \b Description: \n
        *        returns the ivs name of the station that has observations to all other stations.
        *        If several stations have observations to all other stations the one with the most observations is choosen.
        *        If there is no station with observations to all other stations an empty string is returned
        *        If preferredStation is declared this station is used if it has observations to all other stations, 
        *        otherwise the one with most observations is used
        *  \param [in] [std::string] preferredStation 
        * \return [std::string] ivs name of best refernce station
        */
        std::string get_best_refstation( std::string preferredStation = "");
        
        std::map<int, std::vector<double> > compute_sky_coverage(const std::vector<std::vector<unsigned> >& ea_grid_setup,
                                    const lps::TemporalGrid& tg,
                                    std::vector<ivg::Scan>& scans,
                                    std::vector< std::vector< lps::Node<ivg::Obs*, lps::Wedge>  > >&  roots,
                                    double& objective);
        
        std::map<int, std::vector<double> > compute_sky_coverage(const std::vector<std::vector<unsigned> >& ea_grid_setup,
                                    const lps::TemporalGrid& tg,
                                    std::vector< std::vector< lps::Node<ivg::Obs*, lps::Wedge>  > >&  roots,
                                    double& objective);
        
        
    private:

        // ==============================================
        // ============ private methods: ================
        // ==============================================

        /**
        *  \b Description: \n
        *        Method to determine the number of observations from NGS card file
        *  \param [in] [std::string] filename of NGS card file
        */
        void _determine_nobs_ngs( std::string filename );

        /**
        *  \b Description: \n
        *        Method to eliminate stations or sources as well as all corresponding observations.
        *        The selection of the stations or sources is done in the configuration file.
        *  \param [in] no input parameters needed
        */
        void _eliminate_data();

        /**
        *  \b Description: \n
        *        Method to eliminate observations in a defined time period. If the 
        *        time interval is larger than a specific threshold, e.g., 1 hour,
        *        the piece-wise linear parameters in the corresponding have to be
        *        deleted as well. 
        *  \param [in] no input parameters needed
        */
        void _eliminate_obs_period();        
        
        /**
        *  \b Description: \n
        *        Method to handle the modeling of correlations due to atmospheric 
        *        refraction fluctuations (turbulences). The resulting variance-covariance
        *        matrix is used as additional information in the VLBI stochastic model.
        *  \param [in] [std::string] turbulence model
        *                              (0) Matern covariance model (Kermarrec and Schoen, 2014)
	     *                              (1) SIGMA-C model (Schoen and Brunner, 2008)
	     *                              (2) Turbulence model by Treuhaft/Lanyi (1987) and Nilsson et al. (2010)
        *              [double]       Cn - structure constant [m**-2/3] (scaling for turbulence)
        *              [double]       H - effective tropospheric height (integration height)
        *              [ivg::Matrix]  v - wind velocity [m/s**2]
        *              [double]       v_dir - wind direction [gon]; 200.0 == horizontal; 0.0 == zenith
        * \return [ivg::Matrix] variance covariance matrix due to atmospheric turbulences 
        */
        ivg::Matrix _model_turbulence( std::string turbulence_model, std::map< std::string, turbulence_data > turb_sta );

        void _create_weight_matrix( ivg::Matrix * wgt_ptr );

        /**
        *  \b Description: \n
        *        Method to adjust _param_list and TRF, CRF, EOPs to each other. Depends on selected options in configfile.
        *        The selection of the stations or sources is done in the configuration file.
        *  \param [in] int option need to be computed with int ivg::parser::init_options(string option_str)
        */
        vector<double> _adjust_data_storage(int opt);
        
        /**
        *  \b Description: \n
        *        Method to create a string-text-block which contains information about the used configuration setup.
         *       e.g. which sources for NNR, which stations NNR/NNT and so on!
         * 
        *  \param [in] Setting setup, which contains all information about the configuration
         * \return [string] string which can be streamed to cout or whatever
        */
        string _create_configuration_textblock(Setting *setup);
        
        
        void _create_inequality_constraints( std::vector<int> & at_idx );
        bool _exist_negative_zwd();
     
        /**
        *  \b Description: \n
        *        Method to create a matrix for each station including the 
        *        least squares collocation assignmenet.
        *  \param [in] no input parameters needed
        *  \return [std::map< std::string, ivg::Matrix >] G
        */        
        std::vector< ivg::Matrix > _create_collocation_assignment();

        /**
        *  \b Description: \n
        *        Method to create a correlation matrix for stochastic parameters
        *        (e.g., ZWDs or clocks) for filter estimation or least squares collocation
        *  \param [in] no input parameters needed
        *  \return [ivg::Matrix] correlation matrix (for stochastic parameter)
        */        
        ivg::Matrix _create_correlation_fct();
                
        // ==============================================
        // ==============================================

        // settings from config file
        Setting *_setup;
        
        // special handling for specific session (e.g. fix eops)
        Setting *_handling;

        // database name
        std::string _name;
        std::string _name_in_masterfile;
        
        // stringstream containing infos for SNX Block +FILE/COMMENT
        std::stringstream _file_comment;
        
        // session_type ( ngs/ snx / vgosdb )
        std::string _type;

        // vector of scans
        std::vector<ivg::Scan> _scans;

        // TRF and CRF
        ivg::Trf _trf;
        ivg::Crf _crf;

        ivg::Eop_series _eops;

        void *_ephem;

        // vector of strings listing the geophysical effects to be used for the correction of the station position
        std::vector<std::string> _sta_geophys_effects;

        //ivg::Date: start and end epoch of the session
        ivg::Date _start;
        ivg::Date _end;

        // number of observations
        int _nobs;
        
        // number of observations before Elimination of invalid observations
        int _nobs_orig;
        
        // percentage of outliers to be expected, based on delay flag
        double _exp_perc_out;
        
        // vector saving original positions of observations
        // important for correct outlier reading/writing handling
        vector<int> _origin_obs_idxs;

        // Least squares method object
        ivg::Lsa _lsa_solution;
        ivg::Ls_neq _neq_solution;
        ivg::Icls _icls_solution;
        ivg::Lsc _lsc_solution;
        
        ivg::Ls_solution *_solution;
        
        // number of LSQ iterations
        int _niter;

        // apriori values (per observation)
        ivg::Matrix _aprioris;
        
        // parameter list
        ivg::Param_list _param_list;

        // turbulence 
        ivg::Turbulence _turbulence;

        // we need this information if we initialize a trf
        // based on a snx-file. This file must contain STAs and VELs
        map< string,vector<ivg::Date> > _disconts;
        
        // after analysis each station, baseline, source got corresponding residuals for all observations
        // including information like time,az,ele,snr... (used for plotting). These residuals exist for
        // group delay and singleband delay.
        map< delaytype , vector<Residual> > _residuals;
        
        // general statistical information between stations/baselines/sources and corresponding observations
        map< string, ivg::Matrix > _obsstats;
        
        // e.g. /data/vgosDB/2002/02OCT25-A/ in case of _type = vgosdb
        string _session_path;
        
         // default = false; to check whether nnr/nnt is already used or not (for global solution only)
        bool _nnr_nnt_set;
        
        // X or S band
        ivg::band _band_type;
        
        /*
         *  Attributes only needed for ambiguity resolution
         */
        
        //if  true ambiguities are resolved, else group delay is used
        bool _ambigRes;
        bool _phaseDelay;
        bool _SB_solution;
        ivg::Matrix _ambiguity_spacing;
        ivg::Matrix _num_ambig;
        ivg::Matrix _eff_freq;
        
        /*
         *  attributes only needed for lp_sked
         */
        
        lps::Transits _trasintsVisibleFromStation; // transits befor erasing unseen
        lps::Transits _trasintsVisibleFrom2Stations; // transits befor erasing unseen
        lps::Transits _complete_transits; // all visible transits
        lps::Transits _common_transits; // transits seen by at least tow stations
        ivg::Schedule* _schedule_ptr;
};

}
#endif // SESSION_H
