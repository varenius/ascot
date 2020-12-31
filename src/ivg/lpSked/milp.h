#ifndef MILP_H
#define MILP_H

#include <gurobi_c++.h>
#include <vector>
#include <functional>

#include "transits.h"
#include "grid.h"
#include "logger.h"
#include "date.h"
namespace lps{
    
    class Solver; //forward declaration

struct StationActivity{

    int station_idx;
    int source_idx;
    Seconds begin; // seconds since start of session
    Seconds end;

    int index;  //index of interval
    ivg::Date start_epoch;
    ivg::Date end_epoch;
    
    double wrap;
    lps::Seconds obs_duration;
    lps::Seconds align_duration;
    lps::Seconds slew_time;
    
    void print() const{
        std::cerr << "sta: " << station_idx << "  sou:" << source_idx << "  int: " 
                  << index  << "  align: " 
                  << (int)align_duration << "s  slew:"  << (int)slew_time << "s  begin:"
                  << begin << "s  end: " << end << std::endl;
    };
};


// =============================================================================
/**
 * Models the activity of a station.
 */

class Activity{
public:
    
    Activity( Solver* const solver, int index, int numberOfSources, int numberOfIntervals,
                    Seconds lb, Seconds ub, Seconds slewDuration);
    
    Seconds rangeStart() const{return rangeStart_;}
    Seconds rangeEnd()   const {return rangeEnd_;}
        
    int index() const {return index_;}
    
    // Gurobi Getter
    GRBVar & selected() {return grb_selected_;}
    const GRBVar & selected() const {return grb_selected_;}
    
    GRBVar & source_flags(int i) {
        assert(0<=i);
        assert(i<numberOfSources_);
        return grb_source_flags_[i];
    }
    const GRBVar & source_flags(int i) const{
        assert(0<=i);
        assert(i<numberOfSources_);
        return grb_source_flags_[i];
    }
    
    int getActiveSource() const {
        for(int i=0; i < numberOfSources_; ++i){
            if(grb_source_flags_[i].get(GRB_DoubleAttr_X) > 0.001){
                return i;
            }
        }
        return -1;
    }
        
    // These gurobi vars are not existent in the base model
    virtual GRBVar & start_time() {accessError("start_time");};
    virtual const GRBVar & start_time() const{accessError("start_time");};
        
    virtual GRBVar& end_idx() {accessError("end_idx");};
    virtual const GRBVar& end_idx() const{accessError("end_idx");};
    
    virtual GRBVar& start_idx() {accessError("start_idx");};
    virtual const GRBVar& start_idx() const{accessError("start_idx");};
        
    virtual GRBVar & obs_duration() {accessError("obs_duration");};
    virtual const  GRBVar & obs_duration()const {accessError("obs_duration");};
    
    virtual GRBVar & align_duration() {accessError("align_duration");};
    virtual const  GRBVar & align_duration()const {accessError("align_duration");};
       
    virtual GRBVar& wrap() {accessError("wrap");};
    virtual const GRBVar& wrap() const {accessError("wrap");};
    
    virtual GRBVar& slew_align() {accessError("slew_align");};
    virtual const GRBVar& slew_align() const {accessError("slew_align");};
    
    virtual GRBVar& slew_observe() {accessError("slew_observe");};
    virtual const GRBVar& slew_observe() const {accessError("slew_observe");};
    
    virtual GRBVar& slew_dir() {accessError("slew_dir");};
    virtual const GRBVar& slew_dir() const {accessError("slew_dir");};
    
    virtual GRBVar& slew_ele_align() {accessError("slew_ele_align");};
    virtual const GRBVar& slew_ele_align() const {accessError("slew_ele_align");};
    
    virtual GRBVar& slew_duration() {accessError("slew_duration");};
    virtual const GRBVar& slew_duration() const {accessError("slew_duration");} ;

    virtual const GRBVar & end_flags(int i) const{ accessError("end_flags");  };
    virtual GRBVar & end_flags(int i){  accessError("end_flags");  };
    
    // Solution Getter
    virtual double sol_begin( bool tmp = false ) const {return rangeStart_+slewDuration_;};
    virtual double sol_end( bool tmp = false) const  {return rangeEnd_;};
    
    virtual double sol_start_time() const { return rangeStart_+slewDuration_;};
    virtual double sol_obs_duration( bool tmp = false) const { return rangeEnd_ - rangeStart_ - slewDuration_;};
    virtual double sol_wrap() const { return 0.0;};
    virtual double sol_slew_duration() const { return slewDuration_;};
    virtual double sol_align_duration() const { return slewDuration_;};
    
    virtual std::string times2String() const;
    virtual std::string wrap2String() const {return "";};
    
protected:
    // Attributes
    int    index_ = -1;            // index in observations  vector in StationEnvironment Class. correspondes to interval_idx
    int    numberOfSources_ = -1;  // Number of  all sources
    int    numberOfIntervals_ = -1; 
        
    // lower and upper bound for start of observation. reference is start of session 
    Seconds rangeStart_ = -1.0; 
    Seconds rangeEnd_ = -1.0;
    
    Seconds slewDuration_ = 0.0;
    
    // base gurobi vars
    GRBVar*  grb_source_flags_; // array of binary variables each corresponding to one source
    GRBVar   grb_selected_; // 1 if used. else 0 
    
    Solver* solver_;
    
private:
    void accessError(std::string var) const { std::cerr << var << " variable is not accessible in base model!"; exit(EXIT_FAILURE);};
};

class ExtendedActivity : public Activity{

public:
    
    // Constructor 
    ExtendedActivity( Solver* const solver,
                    int index, int numberOfSources, int numberOfIntervals, Seconds lb, Seconds ub,
                    Seconds minLength, Seconds maxLength, Seconds minSlewDuration,
                    std::pair<double, double> wrap_limits );
     
    // Getter
    GRBVar & start_time() override {return grb_obs_start_time_;}
    const GRBVar & start_time() const override {return grb_obs_start_time_;}
    
    GRBVar& end_idx() override {return grb_end_idx_;}
    const GRBVar& end_idx() const override {return grb_end_idx_;}
    
    GRBVar& start_idx() override {return grb_start_idx_;}
    const GRBVar& start_idx() const override {return grb_start_idx_;}
    
    GRBVar & obs_duration() override {return grb_obs_duration_;}
    const  GRBVar & obs_duration()const override {return grb_obs_duration_;}
    
    GRBVar & align_duration() override {return grb_align_duration_;}
    const  GRBVar & align_duration()const override {return grb_align_duration_;}
    
    GRBVar& wrap() override {return grb_wrap_;}
    const GRBVar& wrap() const override {return grb_wrap_;}
    
    GRBVar& slew_align() override {return grb_slew_az_align_;}
    const GRBVar& slew_align() const override {return grb_slew_az_align_;}
    
    GRBVar& slew_observe() override {return grb_slew_az_observe_;}
    const GRBVar& slew_observe() const override {return grb_slew_az_observe_;}
    
    GRBVar& slew_dir() override {return grb_slew_az_dir_;}
    const GRBVar& slew_dir() const override {return grb_slew_az_dir_;}
    
    GRBVar& slew_ele_align() override {return grb_slew_ele_align_;}
    const GRBVar& slew_ele_align() const override {return grb_slew_ele_align_;}
    
    GRBVar& slew_duration() override {return grb_slew_duration_;}
    const GRBVar& slew_duration() const override {return grb_slew_duration_;}

    const GRBVar & end_flags(int i) const override{
        assert(0<=i);
        assert(i<numberOfIntervals_);
        return end_flags_[i];
    }
    GRBVar & end_flags(int i) override{
        assert(0<=i);
        assert(i<numberOfIntervals_);
        return end_flags_[i];
    }

    
    double sol_begin( bool tmp = false ) const override;
    double sol_end( bool tmp = false ) const override;
    double sol_obs_duration( bool tmp = false) const override;
    
    double sol_start_time()  const override{ return grb_obs_start_time_.get(GRB_DoubleAttr_X);};
    double sol_wrap() const override{ return grb_wrap_.get(GRB_DoubleAttr_X);};
    double sol_slew_duration() const override{ return grb_slew_duration_.get(GRB_DoubleAttr_X);};
    double sol_align_duration() const override{ return grb_align_duration_.get(GRB_DoubleAttr_X);};

    std::string times2String() const override;
    std::string wrap2String() const override;
    
private:
    
    GRBVar   grb_obs_start_time_;  // integer variable: start time of activity. t=0 is session start. unit in seconds
    GRBVar   grb_obs_duration_;   // integer variable: duration of observation in seconds
    GRBVar   grb_align_duration_;   // integer variable: total duration required for alignment  (constant + slew part) in seconds
    
    GRBVar grb_start_idx_; // index of alignment start (in observations  vector in StationEnvironment Class)
    GRBVar grb_end_idx_; // index of observation end ( in observations  vector in StationEnvironment Class)

    // Slewing and wrap....................................
    
	GRBVar*  end_flags_; // array of binary variables

    // movment in azimuth
    GRBVar grb_wrap_; // wrap at start of observation (degree)
    GRBVar grb_slew_az_align_; // slew during alignment.  wrap at start of alignment : o(i).wrap - o(i).slew_align (degree)
    GRBVar grb_slew_az_observe_; // slew during observation. wrap at end of observation : o(i).wrap + o(i).slew_observe (degree)
    GRBVar grb_slew_az_dir_; // rotate clockwise or counterclockwise
    
    // movment in elevation
    GRBVar grb_slew_ele_align_;
    
    GRBVar grb_slew_duration_;   // integer variable: duration required to slew telescopes (without constant part )  in seconds
    
};


// =============================================================================
using Activities = std::vector< std::shared_ptr<Activity> >;

// =============================================================================

struct Observation{
    unsigned int station1; 
    unsigned int station2; 
    unsigned int interval_idx; // index in Activities vector in StationEnvironment 
    
    GRBVar  selected; // binary variable 
};


// =============================================================================
/**
 *
 */
struct StationEnvironment{
    Activities activities;
    unsigned int sta_idx;
    
    // constructor for base model
    StationEnvironment(int sta_idx, Seconds internalBegin, Seconds internalEnd, Seconds totalLength, Seconds slewLength,
                                        Solver* const solver, int numberOfQuasars);
    
    // constructor for extended model
    StationEnvironment(int sta_idx, Seconds internalBegin, Seconds internalEnd, Seconds minScanLength, Seconds maxScanLength,
                       Seconds minSlewDuration, Solver* const solver, int numberOfQuasars,
                        std::pair<double, double> wrap_limits);
    
    


};

// =============================================================================
/**
 *   An observatory is an assembly of antennas
 *   The sky above each antenna is nearly identical
 *  
 */
struct Observatory{
    // each element is a root of a tree modeling the partitioning of the sky. For temporal coverage there are several of these roots for each interval
    std::vector< lps::Node<GRBVar,lps::Wedge> >  roots;
    std::vector<GRBVar> cells;
    std::vector<StationEnvironment> stations;
    
    // vector with entries for each tree level(except level 0. First element corresponds to level 1)
    // each entry is the theoretical maximal number of selected cells at a certain level.
    std::vector<double>  max_possible_selected_cells;
    std::vector<double>  max_possible_surface_area;

};

// =============================================================================
class Solver : public GRBCallback
{
public:
    
    
    // constructor
    Solver(ivg::Session * session, const GRBEnv& env);

    // member functions
    void createVariables();
    void createConstraints();
    void createObjective();
    bool optimize(std::function<void(const std::vector<lps::StationActivity> &)> callback,
                  std::function<void (int, int, const lps::Wedge &, int, bool)> treeCallback);

    std::vector<lps::StationActivity> collectResults();
    
    std::vector<lps::StationActivity> collectTemporaryResults();
      
    void printActivityTable() const;
    void milp2session();
    
    // maximal number of intervals that an observation in interval index can reach
    int maxIntsAffectedByObservation(int index )
          { return min( (int)ceil( maxScanLength/intervalLength )+1, numberOfIntervals-index); };
      
    // distance to current interval
    bool objectIsWithinRange( int indexCurrentObject, int index){ return abs(indexCurrentObject-index) <= maxIntsWithoutObs+2; };
//  bool objectIsWithinRange( int indexCurrentObject, int index){ return true; };   
 
         
    GRBModel& getGRBModel() {return grbModel;}
    
    int getMaxIntsWithoutObs() const {return maxIntsWithoutObs;};
    
    double get_solution( const GRBVar& var, bool tmp);
    
private:
    
   double M;
   
    // two stations have to observe simultaneously the same source
    void addGeneralObservationConstrains( bool baseModel = false);
    
   
    /*
     * only visible sources can be observed
     * only one Quasar can be observed at the same time
     * specified has to be past before same source is observed again
     */    
    void addGeneralStationConstraints(bool baseModel = false);
    
    void addSubNetAndSourceRatioConstraints();
    
    void addNumSourceConstraints();
        
    void add_SNR_durationConstraints( bool baseModel = false);
    
    void addNumberOfObservationConstraints();
    
    void addMinTimeWithoutObsConstraints();
    
    void addTwinTelescopeConstraints();
    
    void addSlewConstraints();
    
    // start at neutral point
    void addNeutralPointConstraints();
   
    
   // Attributes
   GRBModel grbModel;
   ivg::Session* const session;
   Setting* setup;
   
   lps::TemporalGrid tg;
   std::vector<std::vector<unsigned>> ea_grid_setup;
   std::vector<unsigned> n_segments;   
   
    // t=0 is start of session
   lps::Seconds internalBegin;
   lps::Seconds internalEnd;
   
   lps::Seconds intervalLength;
   lps::Seconds maxScanLength;
   
   lps::Seconds max_time_without_observation;
   
   bool modelAlignment;

   std::vector<lps::Observatory> observatories;
   // the stations are stored in observatories. This vector contains references to this stations and enables easy looping
   std::vector< lps::StationEnvironment* > stations;
   
   // contains for each baseline and each observation
   //   -a scalar binary variable :true if both stations can observe the same source
   //   -a vector of binary variables: true if quasar is observed (only one quasar can be observed during one scan)
   std::vector<lps::Observation> observations;
   
   std::vector<GRBVar*> subNetFlag;
   
   GRBVar* sourceObservedAtleastOnce;
   
   int numberOfIntervals;
   int n_src;
   int n_sta;
   
   int maxIntsWithoutObs;
      
   // member functions
   std::function<void(const std::vector<lps::StationActivity> &)> callbackFunction;
   std::function<void (int level, int temporal_idx, const lps::Wedge &, int station, bool valid)> treeCallback;

   void createTrees();
   lps::StationActivity createTempActivity(const lps::Activity &a, int station);
   lps::StationActivity createActivity(const lps::Activity &a, int station);
   void checkResults(std::vector<lps::StationActivity> &activity);
   
   void printTimeStats(  const std::vector<lps::StationActivity>& activities ) const;
         
protected:
    void callback () override;
};


} //namespace

#endif // MILP_H
