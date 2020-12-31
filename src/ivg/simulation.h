#ifndef SIMULATION_H
#define SIMULATION_H

#include "trf.h"
#include "crf.h"
#include "session.h"
#include "scan.h"
#include "date.h"
#include "auxfunc.h"
#include "structs.h"
#include "troposphere.h"
#include "hemisphereData.h"
#include "turbulence.h"
#include "schedule.h"
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "tictoc.h"
#include "logger.h"
#include <complex>
#include <fftw3.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <cstdlib>
#include <libconfig.h++>

/**
*
* @brief Simulation class
* @author TA - bakkari developer team
* @date 2015-07-02
* @version 0.1
*
*/
using namespace libconfig;

namespace ivg
{

    class HemisphereData;

// used for encapsulating baseline ( and some epoch) dependend info
struct BandInfo{
    double mean_freq;
    double flux;
    double sefd_sta1;
    double sefd_sta2;
    int num_chan;
};
    
// ===========================================================================
class Simulation
// ===========================================================================
{
    friend class Session;
    friend class Schedule;

    public:

        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Default constructor
        *  \param [in] no input parameters needed
        *  \return An instance of the class 'Simulation'
        */
        Simulation( );
        
        /**
        *  \b Description: \n
        *  \param [in] [ivg::Session] Session for which the simulation 
        *                             should be performed ()
        *  \return An instance of the class 'Simulation'
        */
        Simulation( ivg::Session *session_ptr );

        /**
        *  \b Description: \n
        *        Method to simulate observations (based on entries of config 
        *        file, which is a member variable of the session pointer
        *  \param [in] no input parameters needed
        */
        void simulate();

        void compute_reference_tropo_triangular_mat(  ivg::Trf& trf, const Setting& setup ) const;
        void simulate_reference_troposhere_and_clock( ivg::Trf& trf, double zwd0, const std::vector<std::string>& databases, const Setting& setup );
        ivg::HemisphereData init_hemisphere_data(const Setting& setup,  ivg::Date start = ivg::Date( 0.0 )) const;
        
    private:
//        void _group2names( Setting &params, std::string selected_group_name,
//                           std::vector< std::string > &names, 
//                           std::vector< std::string > &param_lst );
        ivg::Matrix _sim_station_clock( double allan_std, double tau, ivg::Matrix epochs );
        ivg::Matrix _sim_baseline_white_noise( double sigma );
        ivg::Matrix _sim_tropo( std::string model, std::map< std::string, turbulence_data > turb_sta );

        ivg::Matrix _calc_noise_variances( double aVar, double tau, double f_h );
        ivg::Matrix _f_alpha( ivg::Matrix epochs, double Q_d, double alpha );
        
        
        void compute_band_info (ivg::Source& source, ivg::Analysis_station& sta1, Analysis_station& sta2,
                                ivg::Date& epoch, double bit_sampling, ivg::BandInfo& xband, ivg::BandInfo& sband ) const;

        void _calc_standard_deviation(ivg::Source *src_ptr, ivg::Analysis_station *&sta1_ptr, 
                                      ivg::Analysis_station *&sta2_ptr, ivg::Date *epoch_ptr,
                                      double duration, double &sigma_snr_x, double &sigma_ion);
        double _calc_obs_duration(const double & sefd1,const double & sefd2, 
                                  const double & flux,double bw, int num_channels, int bit_samp,
                                  const double & snr_min,const double & sync_time=0.0);
        double _snr_achieved(const double & sefd1, const double & sefd2,
                             const double & flux,double bw, int num_channels, int bit_samp,
                             const double & duration,const double & sync_time=0.0 );

        double _std_sked( double noise, double mean_freq_x, double mean_freq_s, double bw_rms_x, double bw_rms_s, double snr_x, double snr_s  );
                
        ivg::Session* _session_ptr = nullptr;
        
}; // class
}; // namespace

#endif  // SIMULATION_H
