#ifndef SCHEDULE_H
#define	SCHEDULE_H

#include <vector>
#include <iostream>
#include "matrix.h"
#include "simulation.h"
#include "logger.h"
#include "session.h"
#include "analysis_station.h"
#include "crf.h"
#include "scan.h"
#include "obs.h"
#include "tictoc.h"

#include <cstdlib>
#include <libconfig.h++>
#include "lpSked/transits.h"

using namespace libconfig;
/**
*
* @brief class Schedule - base class
* @author AI
* @date 2016-10-20
*/

namespace ivg
{
    enum skedtype{ impact, coverage, minsigma, random };
    
    class Simulation;  // forward declaration
    struct BandInfo;    
    
class Schedule 
{
    
    
    friend class Simulation;
    
    public:

       Schedule();
        
       Schedule( ivg::Session *session_ptr );
       
       void start_scheduling();
        
       void remove_min_time_src(vector<ivg::Scan> &potential_scans, double min_time_src);
        
       ivg::Scan create_neutralpoint_scan(ivg::Date epoch);
        
       vector<ivg::Scan> get_next_potential_scans(ivg::Scan last_scan);
        
       ivg::Scan determine_optimal_scan(ivg::skedtype type, vector<ivg::Scan> potential_scans);
        
       void get_skyplot_data(map<ivg::Analysis_station*,ivg::Matrix> &data, map<ivg::Analysis_station*,vector<string> > &tooltips);
        
       void show_scans(vector<ivg::Scan> &scans, string info = "");
             

       //  ivg::BandInfo& xband, ivg::BandInfo& sband are output variables
       void compute_band_info(ivg::Source& source, ivg::Analysis_station& sta1, Analysis_station& sta2,
                              ivg::Date& epoch, ivg::BandInfo& xband, ivg::BandInfo& sband ) const;
       
       // computes the minimal duration to reach a specific snr
       // if duration is shorter than min duration in config file the latter is used
        double compute_obs_duration(ivg::Source& source, ivg::Analysis_station& sta1, Analysis_station& sta2,
                                    ivg::Date& epoch, ivg::BandInfo& xband, ivg::BandInfo& sband, bool applyMinDur = true) const;
     
        // double& snr_x, double& snr_s are output variables
        double compute_sigma_snr(double duration, ivg::Analysis_station& sta1, Analysis_station& sta2, double el_sta1,
                                double el_sta2,  ivg::Date& epoch, ivg::BandInfo& xband, ivg::BandInfo& sband, double& snr_x, double& snr_s);
        
        void compute_snr(double duration, const ivg::BandInfo& xband, const ivg::BandInfo& sband, double& snr_x, double& snr_s) const;
        
              
        std::vector<unsigned> number_of_stations_source_is_visible(const lps::Transits& transits) const;
        std::vector<bool> quasar_has_at_least_one_bl_with_good_SNR(const lps::Transits& transits) const;
        
        std::map<std::string, ivg::Matrix> minElevation( ivg::Date epoch ) const;
        std::map<std::string, ivg::Matrix> minElevation( ivg::Date epoch, std::vector< std::vector<std::string> > networks) const;
    private:
    
        ivg::Session* _session_ptr;
        
        // final schedulded scans are stored in here
        std::vector<ivg::Scan> _scheduled_scans;
        
        int _bit_sampling;
        double _bandwidth; // in Hz
        
        lps::TemporalGrid _tg;
        std::vector<std::vector<unsigned>> _ea_grid_setup;
};


} // # namespace ivg
#endif

