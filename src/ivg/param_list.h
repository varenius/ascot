#ifndef PARAM_LIST_H
#define	PARAM_LIST_H

/**
*
* @brief class Param_list - base class
* @author AI - bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

#include "logger.h"
#include "trf.h"
#include "crf.h"
#include "ls_solution.h"
#include "lsa.h"
#include "lsc.h"
#include "param.h"

#include <cstdlib>
#include <libconfig.h++>

namespace ivg
{

enum param_def{ aprioris, estimates, totals };     
    
class Param_list
{

    public:

        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Default constructor
        *  \param [in] no input parameters needed
        *  \return An instance of the class 'Param_list'
        */
        Param_list();

        /**
        *  \b Description: \n
        *        copy constructor
        *  \param [in] [ivg::Param_list] other 'Param_list' object
        *  \return An instance of the class 'Param_list'
        */
        Param_list(const Param_list& orig);
        
        /**
        *  \b Description: \n
        *        initialize Param_list object by existing param vector
        *  \param [in] [vector<ivg::Param>] vector of Param
        *  \return An instance of the class 'Param_list'
        */
        Param_list(vector<ivg::Param> params);

        /**
        *  \b Description: \n
        *        initialize Param_list object by existing param vector
        *  \param [in] [vector<ivg::Param>] vector of (deterministic) parameters
        *              [vector<ivg::Param>] vector of stochastic prediction parameters  
        *  \return An instance of the class 'Param_list'
        */
        Param_list(vector<ivg::Param> params,vector<ivg::Param> stoch_params);        
        
        /**
        *  \b Description: \n
        *        constructor using four input parameter: TRF, CRF, EOPs and epoch
        *  \param [in] [ivg::Trf] TRF
        *              [ivg::Crf] CRF
         *             [ivg::Eop_series] EOPs
        *              [ivg::Date] epoch
        *  \return An instance of the class 'Param_list'
        */
        Param_list(ivg::Trf &trf, ivg::Crf &crf, ivg::Eop_series &eops, ivg::Date epoch);

        /**
        *  \b Description: \n
        *        Default deconstructor
        *  \param [in] no input parameters needed
        */
        virtual ~Param_list();


        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Method to modifiy the parameter list according to the config file,
        *        e.g, basic parameters are fixed or reduced or constraint equations are created.
        *  \param [in] [std::string] database name
        *              [Setting]  setup from config file
        *              [ivg::Trf] TRF
        *              [ivg::Crf] CRF
        *              [ivg::Lsa] least squares solver
        */
        void modify_parameterization( std::string dbname, Setting &setup, ivg::Trf &trf, 
                                      ivg::Crf &crf, ivg::Ls_solution &solver, 
                                      ivg::Matrix apriori = ivg::Matrix(),
                                      ivg::Matrix obs_epochs = ivg::Matrix() );

        /**
        *  \b Description: \n
        *        Method to create No-Net-Rotation (NNR) and No-Net-Translation (NNT)
        *        equations to manage for the rank defiecency. The constraint matrix
        *        for the NNR- and NNT-equations is build and the overall
        *        constraint matrix (in ivg::Lsa) is updated.
        *  \param [in] [Setting]  setup from config file
        *              [ivg::Trf] TRF
        *              [ivg::Crf] CRF
        *              [ivg::Lsa] least squares solver
        */
        void create_nnr_nnt_equations( Setting &setup, ivg::Trf &trf, ivg::Crf &crf,
                                       ivg::Ls_solution &solver );
        
        /**
        *  \b Description: \n
        *        Method to create No-Net-Rotation (NNR) and No-Net-Translation (NNT)
        *        constraint matrix B for stations
        *  \param [in] [vector<string>] stations names als ivs_name
        *              [ivg::Matrix] weight matrix w 
        *  \return [out][ivg::Matrix] matrix B [6 x #nparams] 
        */
        ivg::Matrix generate_station_nnr_nnt_B(vector<string> &stations, ivg::Matrix &w);

        /**
        *  \b Description: \n
        *        Method to create (equality) equations to estimate common
        *        clocks (e.g., equal estimates for linear and quadratic term) for
        *        two stations. The constraint matrix for the equations is build 
        *        and the overall constraint matrix (in ivg::Lsa) is updated.
        *  \param [in] [Setting]  setup from config file
        *              [ivg::Ls_solution] least squares solver
        */        
        void create_common_clock_equations( Setting &setup, ivg::Ls_solution &solver );
        
        void merge_zwd_param( ivg::Trf &trf, ivg::Lsa &solver  );
        void merge_grad_param( ivg::Trf &trf, ivg::Lsa &solver  );
        /**
        *  \b Description: \n
        *        Method to create offset and rate constraints. A constraint is build
        *        and the overall constraint matrix (in ivg::Lsa) is updated.
        *  \param [in] [ivg::Ls_solution] least squares solution
        */
        void create_constraint_conditions( ivg::Ls_solution &solver );
        
        /**
        *  \b Description: \n
        *        Method to transform the cpwlf in neq to offset + rate
        *  \param [in] [ivg::Ls_solution] least squares solution
        */
        void transform_cpwlf2offsetrate( ivg::Ls_neq &neq_solution, ivg::Date &t_0 );
        
        /**
        *  \b Description: \n
        *        Method to insert station velocities in param_list
        *  \param [in] [ivg::Ls_solution] least squares solution
        *  \param [in] [ivg::Date] reference epoch for velocity
        *  \param [in] [ivg::Trf] for apriori velocity 
        */
        void insert_station_velocities( ivg::Ls_neq &neq_solution, ivg::Date &ref_epoch, ivg::Trf &trf  );

        /**
        *  \b Description: \n
        *        Method to set start and end epoch
        *  \param [in] [ivg::Date] start epoch
        *              [ivg::Date] end epoch
        */
        void set_start_end_epoch( ivg::Date &start_epoch, ivg::Date &end_epoch );

        /**
        *  \b Description: \n
        *        Method to set the stochastic parameter epochs per station.
        *  \param [in] [std::map<std::string,ivg::Matrix>] stochastic parameter epochs (per station)
        */        
        void set_stoch_param_epoch( std::map<std::string,ivg::Matrix> epochs ){ _stoch_params_epochs = epochs; };
        
        /**
        *  \b Description: \n
        *        Method to set the estimates and the corresponding standard deviation
        *  \param [in] [ivg::Matrix] estimates
        *              [ivg::Matrix] standard deviation of the estimates
        */
        void set_estimates( ivg::Matrix x, ivg::Matrix std, bool prediction = false );
        
        /**
        *  \b Description: \n
        *        Method to set the CPWLF breaks for atmospheric parameters if existent. 
        *        Contains the station and corresponding JD.
        *  \param [in] [multimap<string,double>] clock_breaks
        */
        void set_breaks( map< ivg::paramtype,multimap<string,double> > br){ _breaks = br; };

        /**
        *  \b Description: \n
        *        Method to get the CPWLF breaks for atmospheric parameters if existent. 
        *        Contains the station and corresponding JD.
        *  \param [out] [multimap<string,double>] clock_breaks
        */
        map< ivg::paramtype,multimap<string,double> > get_breaks(){ return _breaks; };      
        
        /**
        *  \b Description: \n
        *        Method to check if at least one parameter of a special paramtype
        *        is existent in the param list vector or not
        *  \param [in] [ivg::paramtype] parameter type
        *  \return [bool] is or is not existent
        */
        bool does_include( ivg::paramtype type );
        
        /**
        *  \b Description: \n
        *        Method to get the index of a parameter selected by a parameter type and name
        *  \param [in] [ivg::paramtype] parameter type
        *              [std::string]    parameter name
        *  \return [int] index of a selected parameter in the parameter list
        */
        int get_index( ivg::paramtype type, string name );
        
        /**
        *  \b Description: \n
        *        Method to get the index of a parameter
        *  \param [in] [ivg::Param] parameter searched for
        *  \return [int] index of the searched parameter in the parameter list
        */
        int get_index( ivg::Param parameter );

        /**
        *  \b Description: \n
        *        Method to get the indices of a list of parameters selected by a list of
        *        (different) parameter types and a name. The size of the indices vector
        *        is equal to the size of the vector including the parameter types.
        *  \param [in] [std::vector<ivg::paramtype>] list of parameter types
        *              [std::string]                 parameter name
        *  \return [std::vector<int>] vector with indices of selected parameters in the parameter list
        */
        vector<int> get_indexes( vector<ivg::paramtype> types, string name );

        /**
        *  \b Description: \n
        *        Method to get the indices of all parameters in the parameter list of a
        *        selected parameter type ('ivg::paramtype type').
        *  \param [in] [std::vector<ivg::paramtype>] list of parameter types
        *              [std::string]                 parameter name
        *  \return [std::vector<int>] vector with indices of all parameters of a selected parameter type
        */
        std::vector<int> get_idx(ivg::paramtype type, std::string name );
        
        /**
        *  \b Description: \n
        *        Method to get a parameter pointer to a given index.
        *  \param [in] [idx] index of the required parameter
        *  \return [ivg::Param] parameter pointer at the given index
        */
        ivg::Param *get_param( int idx );

        /**
        *  \b Description: \n
        *        Method to get station dependent parameter values (e.g., all clock
        *        or atmospheric parameters). Choose between parameter estimates,
        *        apriori values or the total parameter estimates
        *  \param [in] [ivg::paramtype] parameter type (e.g., ivg::paramtype::zwd)
        *              [ivg::param_def] parameter definition (e.g., ivg::param_def::apriori)
        *              [bool] sinex: specify if parameter list was filled based on a SINEX file (default: false)  
        *  \return [map< string, ivg::Matrix >] station dependent parameters
        */
        map< string, ivg::Matrix > get_station_dependent_param( ivg::paramtype type, param_def pdef = param_def::estimates,
                                                                bool prediction = false, bool sinex = false );
        
        /**
        *  \b Description: \n
        *        Method to reduce parameter, whichs has already been flagged
        */
        void reduce_params( ivg::Ls_solution &solver );
        
        /**
        *  \b Description: \n
        *        Method to get the current apriori values from all
        *        parameters in the param_list as a ivg::Matrix (nx1)
        *  \param [in] none
        *  \return [ivg::Matrix] vector including all current apriori values
        */
        ivg::Matrix extract_apriori_vector();
        
        /**
        *  \b Description: \n
        *        Method to get the current estimate values from all
        *        parameters in the param_list as a ivg::Matrix (nx1)
        *  \param [in] none
        *  \return [ivg::Matrix] vector including all current estimate values
        */
        ivg::Matrix extract_estimate_vector();        
        
        /**
        *  \b Description: \n
        *        Method to get datablock containing all data from one specific
        *        parameter, e.g. timeseries of xpo (order=0)
        *        
        *  \param [in] [ivg::paramtype] selected parameter
        *              [int order] order of selected paramater
        *  \return [ivg::Matrix] matrix containing epoch, estimate, aprioi, std
        */
        ivg::Matrix get_param_data(ivg::paramtype type, int order);

        /**
        *  \b Description: \n
        *        Method to get datablock containing all data from one specific
        *        parameter type and name
        *        
        *  \param [in] [ivg::paramtype] selected parameter
        *              [std::string name] paramater name
        *  \return [ivg::Matrix] matrix containing epoch, estimate, aprioi, std
        */
        ivg::Matrix get_param_cpwlf_data( ivg::paramtype type, std::string name );
        
        ivg::Matrix get_param_poly_data( ivg::paramtype type, std::string name );

        vector<ivg::Param> get_stoch_param( ){ return _stoch_params; };
        
        /**
        *  \b Description: \n
        *        Method to check whether there are continuous piece-wise linear 
        *        parameters for a certain parameter type (e.g., clocks or 
        *        ZWDs) and a certain station
        *  \param [in] no input parameters needed
        *  \return [bool] true if there are CPWLF parameters; false if not
        */                
        bool exist_cpwlf_param( ivg::paramtype type, std::string name );

        /**
        *  \b Description: \n
        *        Method to check whether there are polynomial parameters 
        *        for a certain parameter type (e.g., clocks or 
        *        ZWDs) and a certain station
        *  \param [in] no input parameters needed
        *  \return [bool] true if there are polynomial parameters; false if not
        */        
        bool exist_polynomial_param( ivg::paramtype type, std::string name );
        
        /**
        *  \b Description: \n
        *        Method to check whether there are parameters to be estimated 
        *        stochastically for a certain parameter type (e.g., clocks or 
        *        ZWDs) and a certain station
        *  \param [in] no input parameters needed
        *  \return [bool] true if there are stochastic parameters; false if not
        */        
        bool exist_stoch_param( ivg::paramtype type, std::string name );
        
        /**
        *  \b Description: \n
        *        Method to get the size of the parameter list
        *  \param [in] no input parameters needed
        *  \return [int] size of the parameter list
        */
        int size(){return _params.size();};

        /**
        *  \b Description: \n
        *        Method to show the parameter list
        *  \param [in] no input parameters needed
        */
        void show(string out="");

        /**
        *  \b Description: \n
        *        Method to show the results of the estimation process. This includes:
        *        (a) parameter index
        *        (b) parameter name
        *        (c) parameter type
        *        (d) polynomial degree
        *        (e) epoch (MJD)
        *        (f) a priori value
        *        (g) estimate
        *        (h) standard deviation
        *  \param [in] [bool] prediction: show (deterministic) parameters (prediction = false)
        *                                 or (stochastic) prediction (prediction = true)  
        */
        void show_estimates( bool prediction = false );
        
        /**
        *  \b Description: \n
        *        Method to insert a new parameter at a given index.
        *  \param [in] [idx] index for the new parameter
        *              [ivg::Param] new parameter
        */
        void insert_param( int idx, ivg::Param param );

        /**
        *  \b Description: \n
        *        Method to insert a new parameter at a given location (using an iterator).
        *  \param [in] [std::vector<ivg::Param>::iterator] location for the new parameter
        *              [ivg::Param] new parameter
        */
        void insert_param( std::vector<ivg::Param>::iterator it, ivg::Param param );


        /**
        *  \b Description: \n
        *        Method to insert a new baseline clock
        *  \param [in] [std::string] Name of baseline
        */
        void add_bl_clock( std::string bl_clk )
        {
	  _bl_clocks.push_back(bl_clk);
        }
  
        /**
        *  \b Description: \n
        *        Method to remove the parameter at a given index.
        *  \param [in] [idx] index of the parameter to be removed
        */
        void remove_param( int idx );
        
        /**
        *  \b Description: \n
        *        Method to remove a parameter at a given location (using an iterator).
        *  \param [in] [std::vector<ivg::Param>::iterator] location of the parameter to be removed
        */
        void remove_param( std::vector<ivg::Param>::iterator it );

        /**
        *  \b Description: \n
        *        Method to get the pointer (vector<ivg::Param>::iterator) to the first
        *        element of the parameter list
        *  \param [in] [std::vector<ivg::Param>::iterator] pointer to the first element of the parameter list
        */
        vector<ivg::Param>::iterator begin()
        {
            return _params.begin();
        };

        /**
        *  \b Description: \n
        *        Method to get the pointer (vector<ivg::Param>::iterator) to the last
        *        element of the parameter list
        *  \param [in] [std::vector<ivg::Param>::iterator] pointer to the last element of the parameter list
        */
        vector<ivg::Param>::iterator end()
        {
            return _params.end();
        };

        /**
        *  \b Description: \n
        *        Method to get the pointer (vector<ivg::Param>::iterator) to the first
        *        element of the parameter list
        *  \param [in] [std::vector<ivg::Param>::iterator] pointer to the first element of the parameter list
        */
        vector<ivg::Param>::iterator stoch_params_begin()
        {
            return _stoch_params.begin();
        };

        /**
        *  \b Description: \n
        *        Method to get the pointer (vector<ivg::Param>::iterator) to the last
        *        element of the parameter list
        *  \param [in] [std::vector<ivg::Param>::iterator] pointer to the last element of the parameter list
        */
        vector<ivg::Param>::iterator stoch_params_end()
        {
            return _stoch_params.end();
        };        
        void create_vel_constr( string constr_file, Setting &setup, ivg::Ls_solution &solver );
        void set_breaks( string breakfile );
 
       
    private:

        // ==============================================
        // ============ private methods: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Method to get the indices of station groups seleceted in the config file (e.g., to fix the clock for
        *        one station or to use different parametrizations for different groups of stations)
        *  \param [in] [Settings] group settings from config file
        *              [Settings] parameter settings from config file
        *              [std::vector<std::string>] parameter names
        *              [std::string] -
        *              [int] group indices (0 = all)
        *              [std::vector<ivg::paramtype>] parameter types
        *  \return [ivg::Param] parameter at the given index
        */
        vector<int> _get_group_idx( Setting & group, Setting & params,
                                    vector<string> names, string selected_group_name, int name_idx,
                                    vector<ivg::paramtype> type );

        
        void _insert_breaks( ivg::Lsa * solver, ivg::paramtype type );

        // ==============================================
        // ====== class variables / attributes:==========
        // ==============================================

        // list (std::vector) of parameters
        vector<ivg::Param> _params;
        vector<ivg::Param> _stoch_params;
        std::map<std::string,ivg::Matrix> _stoch_params_epochs;

        vector<std::string> _bl_clocks;
        // start and end epoch of the parameter list
        ivg::Date _start_epoch;
        ivg::Date _end_epoch;
       
        // map to save clock breaks (ivg::paramtype::cbr) and CPWLF breaks for 
        // atmospheric parameters (ivg::paramtype::atbr) with station and epoch
        // multiple breaks for a single station are possible        
        map< ivg::paramtype,multimap<string,double> > _breaks;        
};


}
#endif	/* PARAM_LIST_H */

