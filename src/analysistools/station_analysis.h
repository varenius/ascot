#ifndef STATION_ANALYIS_H
#define STATION_ANALYIS_H

#include <vector>
#include <iterator>
#include <sstream>
#include <cstdlib>
#include <map>
#include <utility>
#include "tictoc.h"
#include "matrix.h"
#include "trf.h"
#include "fit.h"
#include "transformation.h"
#include <libconfig.h++>

/**
* @brief Station_analysis class
* @author SH - bakkari developer team
* @date 2015-04-22
* @version 0.1
*/

//using namespace ivg;
using namespace libconfig;
using namespace std;

namespace ivgat
{
    
struct station_epochs
{
    ivg::Matrix ren;
    ivg::Matrix xyz;
    ivg::Matrix xyz_std;
    ivg::Date epoch; 
    int _trf_idx;
};    
    

// ===========================================================================
class Station_analysis
// ===========================================================================
{
   public:

   // ==============================================
   // =============== Constructors: ================
   // ==============================================
   Station_analysis( );
   /**
   *  \b Description: \n
   *        Default constructor
   *  \param [in] no input parameters needed
   *  \return An instance of the class 'Station_analysis'
   */
   Station_analysis( std::vector<ivg::Trf> trf_list );


   // ==============================================
   // =========== MEMBER-functions: ================
   // ==============================================

   /**
   *  \b Description: \n
   *        This method returns the number of sessions for a version/solution
   *  \param [in] no input parameters needed
   */   
   int get_session_number(){ return _trf_list.size(); };

   /**
   *  \b Description: \n
   *        This method returns a vector of sessions in the TRF
   *  \param [in] no input parameters needed
   *  \return [std::vector<std::string> vector of sessions in TRF
   */   
   std::vector<std::string> get_sessions();
   
   /**
   *  \b Description: \n
   *        Method to get the baseline length repeatabilities
   *  \param [in] no input parameters needed
   */
   std::map<std::string, ivg::Matrix> get_bl_rep(){ return _bl_rep; };

   /**
   *  \b Description: \n
   *        This method checks if the estimate values are within a certain
   *        threshold; if not, exclude session from analyis tools
   *  \param [in] [vector<ivg::Trf> trf0: list of TRFs with a-priori values
   *              (needed for calculation of the estimate value from sinex file)
   *              [double] threshold     
   */   
   void check_stations( vector<ivg::Trf> trf0, double th = 0.05 );
   
   /**
   *  \b Description: \n
   *        This method calculates time series of station positions for all stations
   *  \param [in] [vector<ivg::Trf> trf0: list of TRFs with a-priori values
   *              (needed for calculation of the estimate value from sinex file)
   */      
   void calc_station_positions( vector<ivg::Trf> trf0, double th = 0.05  );

   /**
   *  \b Description: \n
   *        This method gets the time series of station positions for a selected station
   *  \param [in] [std::string sta: station name
   *              [ivg::Matrix] Matrix with station positions (call by reference)
   *              [ivg::Matrix] Matrix with standard devitiations of station positions (call by reference)
   *              [ivg::Matrix] Matrix with epochs (call by reference)   
   */         
   void get_station_position( std::string sta, std::string type, 
                              ivg::Matrix & sta_ts, ivg::Matrix & sta_std_ts, 
                              ivg::Matrix & epo_ts );
   /**
   *  \b Description: \n
   *        This method gets the level of uncertainties in terms of the
   *        standard deviations of the station positions
   *  \param [in] [std::string] sta: station name
   *              [std::string] type: XYZ or REN  
   *              [ivg::Matrix] Matrix with standard devitiations of station positions (call by reference)
   */           
  void get_level_of_uncertainties( std::string sta, std::string type,
                                   ivg::Matrix & xyz_std, ivg::Matrix  & epo );
   
   /**
   *  \b Description: \n
   *        This method calculates the baseline lengths for all observations in all sessions (=TRFs)
   *  \param [in] no input parameters needed
   */
   void calc_baseline_length();

   /**
   *  \b Description: \n
   *        This method calculates the baseline length repeatabilites (all sessions = TRFs)
   *  \param [in] no input parameters needed
   */
   void calc_bl_rep( std::vector<std::string> & bl_rep_name, ivg::Matrix & bl_rep, 
                     int n_bl_length, std::string type = "RMS" );

   /**
   *  \b Description: \n
   *        This method gives information about all baselines (baseline name, mjd, baseline length)
   *  \param [in] no input parameters needed
   */
   void show_baselines();

   /**
   *  \b Description: \n
   *        Method to show the baseline lengths repeatabilities (baseline name, mean baseline length, baseline length)
   *  \param [in] no input parameters needed
   */
   void show_bl_repeatabilities();

   /**
   *  \b Description: \n
   *        Method to exclude selected baselines from the solution (i.e., these baselines
    *       are not used for the calculation of the baseline lengths repeatabilities)
   *  \param [in] [std::vector<std::string>] vector with baselines to exclude
   */   
   void remove_baselines( std::vector<std::string> rm_bls ); 
   
   /**
   *  \b Description: \n
   *        Method to fit a quadratic polynomial to the baseline lengths repeatabilities
   *  \param [in] no input parameters needed
   *  \return [map<std::string,ivg::Matrix>] baseline lengths and fitted baseline lengths repeatabilities 
   *                                         for all baselines
   */
   std::map<std::string, ivg::Matrix> fit_bl_rep( ivg::Matrix & bl );
   
   /**
   *  \b Description: \n
   *        Method to calculate the difference in RMS between two different 
   *        solution (i.e., baseline lengths repeatabilities of two solutions)
   *  \param [in] [std::map<std::string,ivg::Matrix> other solution for comparison with this-solution
   *              [std::string] version (name) of solution 1
   *              [std::string] version (name) of solution 2
   */   
   void calc_RMS_differences( std::map< std::string, ivg::Matrix > other, 
                              std::string sol1 = "sol1", std::string sol2 = "sol2",
                              double amount = 1e-3 );

   /**
   *  \b Description: \n
   *        Method to calculate the Helmert parameters between two TRF realizations
   *        (used in different solutions (versions)  
   *  \param [in] [std::vector<ivg::Trf> vector of TRFs (one TRF per session) 
   *              of the second realization
   *              [std::vector<t_param>] list of Helmert parameters to be calculated
   *  \param [out][ivg::Matrix] matrix [mjd x #params] containing the estimated parameter 
   */      
   ivg::Matrix calc_helmert_params( vector<ivg::Trf> other, vector<t_param> tp );
   

   private:

   // ==============================================
   // ==============================================

   // baseline lengths and repeatabilities
   std::vector<ivg::Trf> _trf_list;
   std::vector<bool> _use_trf;
   std::map< std::string, std::vector<station_epochs> > _stations;
   std::vector< std::map<std::string, ivg::Matrix> > _baselines;
   std::vector< std::string > _all_bls;
   std::map<std::string, ivg::Matrix> _bl_rep;
   std::map< std::string, int > _counter;

};

} // # namespace ivgat

#endif // STATION_ANALYIS_H
