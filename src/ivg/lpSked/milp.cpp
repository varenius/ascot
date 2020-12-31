#include "milp.h"
#include <gurobi_c++.h>
#include <set>
#include <functional>
#include "logger.h"

namespace lps{

// ...........................................................................
StationEnvironment::StationEnvironment(int sta_idx, Seconds internalBegin, Seconds internalEnd, Seconds minScanLength, Seconds maxScanLength,
                                       Seconds minSlewDuration, Solver* const solver,
                                       int numberOfQuasars, std::pair<double, double> wrap_limits)
// ...........................................................................
{
    int numberOfIntervals = ceil((internalEnd - internalBegin) / minScanLength);
    activities.reserve(numberOfIntervals);
    this->sta_idx = sta_idx; 
    
    for(int i=0; i < numberOfIntervals; ++i){
        
        activities.push_back( std::make_shared<ExtendedActivity>( ExtendedActivity( solver, i, numberOfQuasars, numberOfIntervals,
                              internalBegin +   i   * minScanLength,      
                              internalBegin + (i+1) * minScanLength,
                              minScanLength, maxScanLength,
                              minSlewDuration, wrap_limits) ));
        
    }
}

// ...........................................................................
StationEnvironment::StationEnvironment(int sta_idx, Seconds internalBegin, Seconds internalEnd, Seconds totalLength, Seconds slewLength,
                                        Solver* const solver, int numberOfQuasars)
// ...........................................................................
{
    int numberOfIntervals = ceil((internalEnd - internalBegin) / totalLength);
    activities.reserve(numberOfIntervals);
    this->sta_idx = sta_idx; 
    
    for(int i=0; i < numberOfIntervals; ++i){
        
        activities.push_back( std::make_shared<Activity>( Activity( solver, i, numberOfQuasars, numberOfIntervals,
                              internalBegin +   i   * totalLength,      
                              internalBegin + (i+1) * totalLength,
                              slewLength))
                            );
        
    }
}

// ...........................................................................
Activity::Activity( Solver* const solver, int index, int numberOfSources, int numberOfIntervals,
                    Seconds lb, Seconds ub, Seconds slewDuration) : 
// ...........................................................................
    solver_(solver),
    index_(index),
    numberOfSources_(numberOfSources),
    numberOfIntervals_(numberOfIntervals),
    rangeStart_(lb),
    rangeEnd_ (ub),
    slewDuration_(slewDuration)

{
    grb_source_flags_   =  solver->getGRBModel().addVars(numberOfSources, GRB_BINARY); // array of binary variables each corresponding to one source
    grb_selected_ =  solver->getGRBModel().addVar(0.0,1.0,0.0,GRB_BINARY); 

}
// ...........................................................................
ExtendedActivity::ExtendedActivity( Solver* const solver,
                    int index, int numberOfSources, int numberOfIntervals, Seconds lb, Seconds ub,
                    Seconds minLength, Seconds maxLength, Seconds minSlewDuration,
                    std::pair<double, double> wrap_limits ) :
// ...........................................................................
    Activity(solver, index, numberOfSources, numberOfIntervals, lb, ub, minSlewDuration) {
    
    grb_obs_duration_   =  solver->getGRBModel().addVar(minLength, maxLength, 0.0,GRB_CONTINUOUS);
    
    
    // lb is 0.0 and not minSlewDuration. limited by slew_duration_ implicitly
    // if lb is set to minSlewDuration model is infeasible if minSlewDuration is larger than interval length/min scan length
    grb_align_duration_ =  solver->getGRBModel().addVar(0.0, 1800.0, 0.0,GRB_CONTINUOUS); 
    grb_obs_start_time_ =  solver->getGRBModel().addVar(lb, ub, 0.0, GRB_CONTINUOUS);
    
    grb_start_idx_ =  solver->getGRBModel().addVar(max(0, index - solver->getMaxIntsWithoutObs()-1 ), index, 0.0,GRB_INTEGER);
    grb_end_idx_ =  solver->getGRBModel().addVar(index,  index+solver->maxIntsAffectedByObservation(index) , 0.0,GRB_INTEGER);

	end_flags_ =  solver->getGRBModel().addVars( solver->maxIntsAffectedByObservation(index), GRB_BINARY);
               
    grb_wrap_ =  solver->getGRBModel().addVar(wrap_limits.first, wrap_limits.second, 0.0, GRB_CONTINUOUS);
    grb_slew_az_align_ =  solver->getGRBModel().addVar(-360.0, 360.0, 0.0, GRB_CONTINUOUS);
    grb_slew_az_observe_ =  solver->getGRBModel().addVar(-10.0, 10.0, 0.0, GRB_CONTINUOUS);
    grb_slew_az_dir_  =  solver->getGRBModel().addVar(0.0,1.0,0.0,GRB_BINARY);
    grb_slew_ele_align_ =  solver->getGRBModel().addVar(-90.0, 90.0, 0.0, GRB_CONTINUOUS);
    
    grb_slew_duration_ =  solver->getGRBModel().addVar(minSlewDuration, 1800.0, 0.0, GRB_CONTINUOUS);
    
}

double ExtendedActivity::sol_begin( bool tmp ) const { return Activity::solver_->get_solution(grb_obs_start_time_, tmp); }
double ExtendedActivity::sol_end( bool tmp ) const { return Activity::solver_->get_solution(grb_obs_start_time_, tmp) + Activity::solver_->get_solution(grb_obs_duration_, tmp); }
double ExtendedActivity::sol_obs_duration( bool tmp) const { return Activity::solver_->get_solution(grb_obs_duration_, tmp);}

// ...........................................................................
double Solver::get_solution( const GRBVar& var, bool tmp){
// ...........................................................................
    if(tmp){
        return getSolution(var);
    }
    return var.get(GRB_DoubleAttr_X);
}

// ...........................................................................
std::string Activity::times2String() const{
// ...........................................................................
    stringstream ss;
    
    
    std::string c_norm = ivg::Logger::get_color("white");
    std::string c_high;
  
    if( this->grb_selected_.get(GRB_DoubleAttr_X) > 0.1 ){
        c_high = ivg::Logger::get_color("boldgreen");
    } else {
        c_high = c_norm;
    }
    
    ss << c_high
       << setw(2) << (unsigned) abs(this->grb_selected_.get(GRB_DoubleAttr_X))
       << c_norm << " " << c_high
       << setw(5) << this->getActiveSource()
       << c_norm 
       << " ";


    return ss.str();
}

// ...........................................................................
std::string ExtendedActivity::times2String() const{
// ...........................................................................
    stringstream ss;
    
    
    std::string c_norm = ivg::Logger::get_color("white");
    std::string c_high;
  
    if( this->grb_selected_.get(GRB_DoubleAttr_X) > 0.1 ){
        c_high = ivg::Logger::get_color("boldgreen");
    } else {
        c_high = c_norm;
    }
    
    ss << c_high
       << setw(2) << (unsigned) abs(this->grb_selected_.get(GRB_DoubleAttr_X))
       << c_norm << " : " << c_high
       << setw(5) << round( (this->grb_obs_start_time_.get(GRB_DoubleAttr_X) - this->grb_align_duration_.get(GRB_DoubleAttr_X)) ) << " "
       << setw(5) << round( this->grb_obs_start_time_.get(GRB_DoubleAttr_X) )<< " "
       << setw(5) << round( (this->grb_obs_start_time_.get(GRB_DoubleAttr_X) + this->grb_obs_duration_.get(GRB_DoubleAttr_X) ))
       << c_norm << " : " << c_high
       << setw(4) << abs(this->grb_start_idx_.get(GRB_DoubleAttr_X)) << " "
       << setw(4) << this->index_ << " "
       << setw(4) << abs(this->grb_end_idx_.get(GRB_DoubleAttr_X))
       << c_norm << " : " << c_high
       << setw(4) << round( this->grb_align_duration_.get(GRB_DoubleAttr_X) ) << "/"
       << setw(4) << left <<  round( this->grb_slew_duration_.get(GRB_DoubleAttr_X) ) << right << " "
       << setw(4) << round( this->grb_obs_duration_.get(GRB_DoubleAttr_X) ) << " "
       << c_norm << " : " << c_high
       << setw(5) << this->getActiveSource()
       << c_norm;


    return ss.str();
}

// ...........................................................................
std::string ExtendedActivity::wrap2String() const{
// ...........................................................................
    stringstream ss;
    
    
    std::string c_norm = ivg::Logger::get_color("white");
    std::string c_high;
  
    if( this->grb_selected_.get(GRB_DoubleAttr_X) > 0.1 ){
        c_high = ivg::Logger::get_color("boldblue");
    } else {
        c_high = c_norm;
    }
    
    ss << c_high
       << setw(2) << (unsigned) abs(this->grb_selected_.get(GRB_DoubleAttr_X))
       << c_norm << " : " << c_high
       << setw(7) << setprecision(1) << std::fixed << (this->grb_wrap_.get(GRB_DoubleAttr_X) - this->grb_slew_az_align_.get(GRB_DoubleAttr_X)) << " "
       << setw(7) << setprecision(1) << std::fixed << this->grb_wrap_.get(GRB_DoubleAttr_X) << " "
       << setw(7) << setprecision(1) << std::fixed << (this->grb_wrap_.get(GRB_DoubleAttr_X) + this->grb_slew_az_observe_.get(GRB_DoubleAttr_X))
       << c_norm << " : " << c_high
       << setw(7) << setprecision(2) << std::fixed << this->grb_slew_az_align_.get(GRB_DoubleAttr_X) << " "
       << setw(8) << setprecision(3) << std::fixed << this->grb_slew_az_observe_.get(GRB_DoubleAttr_X) << " "
       << c_norm << " : " << c_high
       << setw(6) << setprecision(1) << std::fixed << this->grb_slew_ele_align_.get(GRB_DoubleAttr_X) << " "
       << c_norm << " : " << c_high
       << setw(5) << this->getActiveSource()
       << c_norm;

    return ss.str();
}

// ...........................................................................
Solver::Solver(ivg::Session * ses, const GRBEnv& env): session(ses), grbModel(env)
// ...........................................................................
{
    
   // t=0 is start of session
   internalBegin = 0;
   internalEnd   = (session->getEnd().get_double_mjd() - session->getStart().get_double_mjd())*lps::DurationDay;
   
    n_src = session->get_crf_ptr()->get_number_sources_inc_unused();
    n_sta = session->get_trf_ptr()->get_number_stations(); 
    
    setup = session->get_setup();
    
    lps::Seconds temporal_grid_resolution = (int)(*session->get_setup())["SKED"]["temporal_resolution"] * 60;
    lps::Seconds temporal_grid_shift = (int)(*session->get_setup())["SKED"]["temporal_shift"] * 60;
    this->tg = lps::TemporalGrid(internalEnd-internalBegin, temporal_grid_resolution, temporal_grid_shift, session->getStart());
    
    //  M > session dur (86400 is good value)
    // However, for shorter sessions smaller value leads to smaller Matrix range. Better numeric.
    if( (int) (*session->get_setup())["SKED"]["approach"] == 4 )
        M = max(internalEnd-internalBegin + 100, 10000.0); // used for 'constraint only if constraints
    else 
        M = 10000.0;
    
    
    bool ea_grid  = (*session->get_setup())["SKED"].exists("EqualAreaGrid");
    

    if (ea_grid){
        Setting &ea = (*session->get_setup())["SKED"]["EqualAreaGrid"];
        n_segments.resize(ea.getLength());
        ea_grid_setup.resize(ea.getLength());
    
        for(int i=0; i<ea.getLength(); i++){
            for(int j=0; j<ea[i].getLength(); j++){
                ea_grid_setup[i].push_back( (unsigned)ea[i][j] );
                n_segments[i] += (unsigned)ea[i][j];
            }
        }
    }
    
        
}

// ...........................................................................
void Solver::createVariables()
// ...........................................................................
{
           
    Seconds minScanLength = (*setup)["SKED"]["min_scan"];
    maxScanLength = (*setup)["SKED"]["max_scan"];
    
    max_time_without_observation = (double)(*setup)["SKED"]["max_time_without_observation"];    

    lps::Seconds min_slew_time =  (double)(*setup)["SKED"]["const_setup"];
    
    modelAlignment =  (bool)(*setup)["SKED"]["model_alignment"];
    
    //only for 5
    Seconds alignmentTime = 30.0;
                        
    log<INFO> ("*** CREATE VARIABLES ");
    log<DETAIL> ("*** number of stations: ") % n_sta;
    log<DETAIL> ("*** number of sources: ") % n_src;
    
    if( (int)(*setup)["SKED"]["approach"] == 4){
        intervalLength = minScanLength;
    } else if( (int)(*setup)["SKED"]["approach"] == 5){
        intervalLength = minScanLength + alignmentTime;
    }
    
    numberOfIntervals = ceil((internalEnd-internalBegin)/intervalLength);
    
    
    maxIntsWithoutObs = ceil (max_time_without_observation/intervalLength);
    std::cout << "maximal number of intervals without observation: " << maxIntsWithoutObs << " interval length " << intervalLength 
              << " maximal duration without observation " << max_time_without_observation << std::endl;
        
    log<DETAIL>( "*** begin: ") % internalBegin % " [sec] end: " % internalEnd % " [sec]";
    log<DETAIL>( "*** minimum scan length: ") % minScanLength % " [sec]";
    log<DETAIL>( "*** maximum scan length: ") % maxScanLength % " [sec]";
    log<DETAIL>( "*** number of intervals: ") % numberOfIntervals;
    
    // create variables
    
    observatories.reserve(session->get_trf_ptr()->get_station_twin_list().size());
    
    for(std::vector<ivg::Analysis_station*>& stations_at_observatory : session->get_trf_ptr()->get_station_twin_list()){
        observatories.push_back(lps::Observatory());
        observatories.back().stations.reserve(stations_at_observatory.size());
        for(ivg::Analysis_station* sta: stations_at_observatory){
            int sta_idx = sta->get_idx();
            log<DETAIL> ("*** create environment for station ") % session->get_trf_ptr()->get_station(sta_idx)->get_name(ivg::staname::ivs_name);

            ivg::Antenna antenna = session->get_trf_ptr()->get_station(sta_idx)->get_antenna_info();
            
            if(!(bool)(*setup)["SKED"]["model_alignment"] && (bool)(*setup)["SKED"]["assume_worst_slew"]){
                min_slew_time = 2*M_PI/antenna.azi_rate;
                log<DETAIL> ("*** 360 slew time: ") % min_slew_time % "s";
            }
            
            std::pair<double, double> wrap_limits (antenna.azi_min*ivg::rad2d+1, antenna.azi_max*ivg::rad2d-1 ); // degree
            log<DETAIL>( "*** wrap limits (1 deg buffer): ") % wrap_limits.first % " " % wrap_limits.second;
            log<DETAIL>( "*** slew rate: ") %  (antenna.azi_rate*ivg::rad2d) % " [deg/sec]"; 
            
            // General case
            if( (int)(*setup)["SKED"]["approach"] == 4){ 
                observatories.back().stations.push_back( StationEnvironment(sta_idx, internalBegin, internalEnd,
                                                  minScanLength, maxScanLength, min_slew_time, this, n_src,
                                                  wrap_limits));
                
            } else if( (int)(*setup)["SKED"]["approach"] == 5){
                observatories.back().stations.push_back( StationEnvironment(sta_idx, internalBegin, internalEnd,
                                                  intervalLength, alignmentTime, this, n_src ));
            }
        }
    }
    
    // stations stored here have to be in same order as in TRF class
    stations.resize(n_sta);
    
    for(lps::Observatory& observatory : observatories){
        for(lps::StationEnvironment& se: observatory.stations){
            stations[se.sta_idx] = &se;
        }
    }
    
    double ratioUsedSources = 0.0;
    if( (*setup)["SKED"].exists("RatioUsedSources") ){
        ratioUsedSources = (double)(*setup)["SKED"]["RatioUsedSources"]; 
    }

    if( (int)(*setup)["SKED"]["MaxSubNets"] > 0 || ratioUsedSources > 0.0){
        for(int k=0; k < numberOfIntervals; ++k){
            subNetFlag.push_back( grbModel.addVars(n_src, GRB_BINARY) );
        }
    }
    
    
    sourceObservedAtleastOnce = grbModel.addVars(n_src, GRB_BINARY);
        
    std::vector<std::pair<std::string, std::string>> exclude_bl (0);
    if(  (*setup)["SKED"].exists("exclude_bl") ){
        exclude_bl = get_baselines( (*setup)["SKED"]["exclude_bl"]);
    }
    
    std::vector< std::vector<std::pair<std::string, std::string>>> force_bl (0);
    
    if(  (*setup)["SKED"].exists("force_bl") ){
        for(int i=0; i < (*setup)["SKED"]["force_bl"].getLength(); i++){
            force_bl.push_back( get_baselines( (*setup)["SKED"]["force_bl"][i])  );
        }  
    }
    
    bool twinObserveTogether = false;
    if(  (*setup)["SKED"].exists("allow_twin_same_source") ){
        twinObserveTogether = (bool)(*setup)["SKED"]["allow_twin_same_source"];
    }
    
    // loop over all possible baselines (station_i - station_j)
    for(size_t sta_idx_1 = 0; sta_idx_1 < n_sta; ++sta_idx_1) {
        for(size_t sta_idx_2 = sta_idx_1 + 1; sta_idx_2 < n_sta; ++sta_idx_2){
            
            std::string sta_name_1 = session->get_trf_ptr()->get_station( sta_idx_1 )->get_name(ivg::staname::ivs_name);
            std::string sta_name_2 = session->get_trf_ptr()->get_station( sta_idx_2 )->get_name(ivg::staname::ivs_name);
            bool exclude = includes_baseline(exclude_bl, sta_name_1, sta_name_2 );
            if(exclude)
                continue;
            
            int nfbl =  force_bl.size();
            std::vector<bool> force_bl_at( nfbl );
            for(int i=0; i < nfbl; i++){
                force_bl_at[i] = includes_baseline(  force_bl[i], sta_name_1, sta_name_2 );
            }
                        
            // Do not add observations between TWIN STATIONS
            if( !(twinObserveTogether == false && session->get_trf_ptr()->areTwins(sta_idx_1, sta_idx_2) == true) ){
                log<DETAIL>( "*** adding baseline: ") % sta_name_1 % " " % sta_name_2;

                // loop over all observations
                for(int k=0; k < numberOfIntervals; ++k){
                    Observation obs;
                    obs.station1 = sta_idx_1;
                    obs.station2 = sta_idx_2;
                    obs.interval_idx = k;
                    
                    bool isForced(false), isForcedAtK(false);
                    for(int i=0; i < nfbl; i++){
                        if (force_bl_at[i] ){
                            isForced = true;
                            if ( (k+i) % nfbl  == 0 ){
                                isForcedAtK = true;
                                break;
                            }
                        }
                    }
                    
                    if(isForced){
                        if (isForcedAtK){
                            obs.selected = grbModel.addVar(1.0, 1.0, 0.0, GRB_BINARY);
                        } else {
                            obs.selected = grbModel.addVar(0.0, 0.0, 0.0, GRB_BINARY);
                        }
                    
                    } else {
                        obs.selected = grbModel.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
       
                    observations.push_back(obs);
                }
            }
        }
    }

    createTrees();

    grbModel.update();

}

// ...........................................................................
void Solver::createTrees(){
// ...........................................................................
     
    bool ea_grid  = (*session->get_setup())["SKED"].exists("EqualAreaGrid");  

    int max_tree_level = (int)(*session->get_setup())["SKED"]["max_tree_level"];
    
    for(lps::Observatory& observatory : observatories){

        observatory.max_possible_selected_cells = std::vector<double>(max_tree_level, 0.0);
        observatory.max_possible_surface_area = std::vector<double>(max_tree_level, 0.0);
        
        unsigned num_temp_ints = tg.get_number_of_intervals();
        observatory.roots.reserve(num_temp_ints);
        
        int u = 0;
        for(int t = 0; t < num_temp_ints ; ++t){
            std::vector<lps::DataPoint<GRBVar>> points;
            for(lps::StationEnvironment& se: observatory.stations){

                std::pair<unsigned, unsigned> activity_range =  tg.getIndexOtherSampling(t, intervalLength); 
//                std::cout << activity_range.first << " " << activity_range.second << std::endl;

                
                for(unsigned activity_idx = activity_range.first; activity_idx <= activity_range.second && activity_idx < se.activities.size(); ++ activity_idx ){
                    Activity& o = *se.activities[activity_idx];
                    for(size_t src_idx=0; src_idx <n_src; ++src_idx){
                        ivg::Source& src = *session->get_crf_ptr()->get_source(src_idx);
                        std::pair<lps::Position, lps::Position> pos;
                        if( session->get_common_transits().get_border_pos(src_idx, se.sta_idx, o.index()*intervalLength, (o.index()+1)*intervalLength, intervalLength, pos ) ){
                            lps::Point p (pos.first.azimuth(), pos.first.elevation());
                            points.push_back(lps::DataPoint<GRBVar>(p, o.source_flags(src_idx) ) );
                        }
                    }
                }

            }
            
            //level zero is full circle
            lps::Wedge wedge (lps::Point(0,0),0,90,0,360);

            if(ea_grid){
                observatory.roots.push_back( lps::Node<GRBVar,lps::Wedge>(points, wedge, ea_grid_setup) );
            } else{
                observatory.roots.push_back( lps::Node<GRBVar,lps::Wedge>(points, wedge, 0 ,max_tree_level) );
            }

            int u = 0;
            observatory.roots.back().visitBoundAndBool([&](int level, const lps::Wedge& w ,bool valid){  
                if( valid ){
                    observatory.cells.push_back( grbModel.addVar(0,1,0,GRB_BINARY) );
                    observatory.max_possible_selected_cells[level-1] += 1.0;
                    observatory.max_possible_surface_area[level-1] += w.surfaceArea();
                    ++u;
                    //w.print();
                } 
            });
        }
        
    } 
   
}

// ...........................................................................
void Solver::createObjective(){
// ...........................................................................
    
    log<INFO>( "*** CREATE OBJECTIVE "); 
    
    bool objective_surface =  (bool)(*setup)["SKED"]["objective_surface"];
    bool objective_relative =  (bool)(*setup)["SKED"]["objective_relative"];
    double objective_nobs_weight =  (double)(*setup)["SKED"]["objective_nobs_weight"];
    
    GRBLinExpr objectiveFunction;
    for(lps::Observatory& observatory : observatories){
        int u = 0;
        for(lps::Node<GRBVar,lps::Wedge>& root: observatory.roots){

            root.visitAll([&](const std::vector<lps::DataPoint<GRBVar>> & points, int level, const lps::Wedge& w , bool valid ){

                if(valid){

                    // only add gain if cell has observation
                    GRBLinExpr expr;
                    for(const lps::DataPoint<GRBVar>& point : points){
                        expr += point.object();
                    }
                    grbModel.addConstr( observatory.cells[u]  <=  expr );

                    if( objective_surface ){
                        if (objective_relative){
                            objectiveFunction += observatory.cells[u] * w.surfaceArea() / observatory.max_possible_surface_area[level-1];
                        } else {
                            objectiveFunction += observatory.cells[u] * w.surfaceArea() / (2*M_PI);
                        }
                    } else {
                        if (objective_relative){
                            objectiveFunction += observatory.cells[u] / observatory.max_possible_selected_cells[level-1];
                        } else {
                            objectiveFunction += observatory.cells[u] / n_segments[level-1];
                        }
                    }
                    u++;
                }

            });
        }

    }
    
    /*
     *  maximize number of observations, too
     *  important for schedules with more than 2 stations
     *  example: 3 Stations -> 3 baselines (observations)
     *     -> 2 obs are enough to create a point on all 3 sky plots
     *     however, maybe a third observation is possible and missed 
     * 
     * also important for sessions with only 2 stations
     *  no need to create possible observation if it is in same rect
     */
    if( objective_nobs_weight > 0.0){
        for(Observation& obs : observations){
            objectiveFunction += objective_nobs_weight * obs.selected;
        }
    }
        
    grbModel.setObjective(objectiveFunction,GRB_MAXIMIZE);
}

// ...........................................................................
void Solver::createConstraints()
// ...........................................................................
{
    log<INFO>( "*** CREATE CONSTRAINTS ");      
    

    // General case
    if( (int)(*setup)["SKED"]["approach"] == 4){ 
 
        addNumberOfObservationConstraints();                              
        addMinTimeWithoutObsConstraints();
        
         // ["SKED"]["const_setup"] is not added because it is already  included in variable limits
        lps::Seconds constAlignTime =  (double)(*setup)["SKED"]["const_idle"]
                                 + (double)(*setup)["SKED"]["const_source"]
                                 + (double)(*setup)["SKED"]["const_calib"]; // constant after slewing and before observation [sec] 
        
        for(size_t sta_idx = 0; sta_idx < stations.size(); ++sta_idx){
            StationEnvironment & se = *stations[sta_idx];  
            ivg::Antenna antenna = session->get_trf_ptr()->get_station(sta_idx)->get_antenna_info();

            for(int i=0; i < numberOfIntervals; ++i){  

                GRBLinExpr sumSelected = 0;

                Activity& o1 = *se.activities[i];
                
                grbModel.addConstr( o1.start_time() - o1.align_duration() >= 0.0 );

                /*******************************************************************
                 * end_idx and start_idx 
                 * 
                 */
                grbModel.addConstr(   o1.end_idx()      * intervalLength + 0.1 <= o1.start_time() + o1.obs_duration() );
                grbModel.addConstr( ( o1.end_idx() + 1 ) * intervalLength>= o1.start_time() + o1.obs_duration() );

                grbModel.addConstr(   o1.start_idx()      * intervalLength <= o1.start_time() - o1.align_duration() );
                grbModel.addConstr( ( o1.start_idx() + 1 ) * intervalLength -0.1 >= o1.start_time() - o1.align_duration() );

             	/*******************************************************************
                 * end_flag
                 * 
                 * i + j*o(i).end_flag(j)  = o(i).end_idx
                 * 
                 */
                GRBLinExpr sumEndFlag = 0;
                for(int j=0; j <  maxIntsAffectedByObservation(i); ++j){
                    grbModel.addConstr( o1.end_idx() - (i+j)  - ( 1-o1.end_flags(j) )*M <= M-M*o1.selected() );
                    grbModel.addConstr( (i+j) - o1.end_idx()  - ( 1-o1.end_flags(j) )*M <= M-M*o1.selected() );
                    sumEndFlag += o1.end_flags(j);
                }
                grbModel.addConstr( 1-sumEndFlag  <=  M-M*o1.selected() );



                /*******************************************************************
                 * obs between end_idx and current index can not be used
                 * o1.end_idx() - (i+j) <= M*(1-se.observations[i+j].selected())
                 * 
                 * use o1.end_idx()+1 if end_idx obs should not be used too.
                 * 
                 */

                for(int j=1; j < maxIntsAffectedByObservation(i); ++j){
                    grbModel.addConstr( o1.end_idx() - (i+j) - M*(1-se.activities[i+j]->selected()) <= 2*M-2*M*o1.selected() );
                }

                /*******************************************************************
                 * obs between start and current index can not be used
                 * j - o1.start_idx() <= M*(1-se.observations[j].selected())
                 * 
                 * use o1.start_idx()-1 if start_idx obs should not be used too. 
                 */

                for(int k=i-1; k >= 0 && objectIsWithinRange(i,k); --k){
                    grbModel.addConstr(  k - o1.start_idx() - M*(1-se.activities[k]->selected()) <= 2*M-2*M*o1.selected() );
                }

                /*******************************************************************
                 * end of an observation == start of next alignment
                 * 
                 * If o1.selected & o2.selected & o2. is the next selected observation after o1
                 *      -> o1.start_time() + o1.obs_duration() == o2.start_time() - o2.align_duration()
                 * 
                 * 
                 * 
                 */
                sumSelected = 0;
                for(int k=i+1; k < numberOfIntervals && objectIsWithinRange(i,k); ++k){
                    Activity& o2 = *se.activities[k];
                    // this expression is zero if and only if o1 and o2 are selected and o2 is the first selected observation after o1. Else this expression is M or even larger
                    //GRBLinExpr firstSelectedAfterCurrentSelected = M*( 1-o1.selected() ) + M*( 1-o2.selected() )  + M*sumSelected;

                    grbModel.addConstr( (o1.start_time() + o1.obs_duration()  ) - (o2.start_time() - o2.align_duration()) <=  M*( 1-o1.selected() ) + M*( 1-o2.selected() )  );
                   
                     
                }


                /*******************************************************************
                 * align_duration duration
                 *  slew_duration + constant parts
                 */
                // simplification 
                lps::Seconds antenna_const = max(antenna.azi_const, antenna.ele_const);

                grbModel.addConstr( o1.align_duration() - o1.slew_duration() - constAlignTime -antenna_const  >= -M+M*o1.selected() );


            }

        }

        addGeneralStationConstraints();
        addGeneralObservationConstrains();
        add_SNR_durationConstraints();
        addTwinTelescopeConstraints(); 

        if(modelAlignment){
           addNeutralPointConstraints();
           addSlewConstraints();
        }  
   
    } else if ( (int)(*setup)["SKED"]["approach"] == 5){ // vgos easy model
            

        addNumberOfObservationConstraints();
        addGeneralStationConstraints(true);
        addGeneralObservationConstrains(true);
        add_SNR_durationConstraints(true);
        addTwinTelescopeConstraints(); 
        addSubNetAndSourceRatioConstraints();
        addNumSourceConstraints();
    }
   
    
}

// ...........................................................................
void Solver::addGeneralStationConstraints( bool baseModel){
// ...........................................................................
    
    log<DETAIL>( "*** add general station constraints");
    
    lps::Seconds min_time_src = (double)(*setup)["SKED"]["min_time_src"];
    
    for(size_t sta_idx = 0; sta_idx < stations.size(); ++sta_idx){
        StationEnvironment & se = *stations[sta_idx];  
                
        for(int i=0; i < numberOfIntervals; ++i){  
            
            Activity& o1 = *se.activities[i];
            
            /*******************************************************************
             * Observations must start when source is visible
             * Observations must end before source is not visible any longer
             * 
             * if j-th quasar is observed add two constraints:
             *     (1) begin of the transit has to be before the observation is started:
             *          transit.begin <  o1.start_time
             *     (2) end of transit has to be after end of observation
             *          transit.end >  o1.start_time + o1.duration
             * else
             *      constraints are always satisfied because M is large
             * 
             * Four possibilities
             * Transit
             *  (1) is present in complete interval -> usable
             *  (2) starts in interval              -> maybe usable
             *  (3) ends in the interval            -> not usable
             *  (4) starts and ends in the interval -> not usable
             */

            Seconds lowerBound = i*intervalLength;
            Seconds upperBound = (i+1)*intervalLength;
            
            if(baseModel)
                lowerBound += o1.sol_slew_duration();
            
            // loop over all sources
            for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
                int counter = 0;
                for(const lps::Transit& transit : session->get_common_transits().transit(scr_idx, sta_idx)){

                    Seconds tbegin = transit.begin();
                    Seconds tend = transit.end();

                    if(baseModel){
                        if( tbegin-1E-4 <= lowerBound && tend+1E-4 >= upperBound){
                            assert(counter==0);
                            ++counter;
                        }
                    } else {
                        if( ( tbegin <= lowerBound && tend >= upperBound) ||
                            ( tbegin >  lowerBound && tbegin < upperBound && tend >= upperBound ) ){
                            assert(counter==0);
                                grbModel.addConstr( tbegin  -o1.start_time()  <=  M - M * o1.source_flags(scr_idx));
                                grbModel.addConstr( o1.start_time() + o1.obs_duration() - tend  <=  M - M * o1.source_flags(scr_idx));
                            ++counter;
                        }        
                    }               
                }

                // Do not use source for the current observation
                //
                // source_flags <= 0 is equal to source_flags == 0 because source_flags is binary [0,1]
                //
                if(counter == 0){
                    grbModel.addConstr( o1.source_flags(scr_idx) <= 0);
                } 
            }

            /*******************************************************************
             * only one quasar can be selected at the same time
             * if activity is selected exactly one entry in source array is one
             * Otherwise all entries are zero
             *  sum(source_flags) <= 1 if selected
             *  sum(source_flags) <= 0 if not selected
             */
            GRBLinExpr exclusiveObservation;
            for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
                exclusiveObservation += o1.source_flags(scr_idx);
            }
            grbModel.addConstr( exclusiveObservation == o1.selected() );
            
            
            /******************************************************************* 
             * time between observations of the same source has to be larger than specified duration
             */
            for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
                for(int j = i-1; j>=0; --j){ 

                    // condition is alway satisfied if certain number of observations has past
                    if((i - j)*intervalLength > min_time_src )
                        break;
                    
                    if(baseModel){
                        grbModel.addConstr( se.activities[j]->sol_begin() - o1.sol_begin()   + min_time_src <=  
                                        2*M-M*o1.source_flags(scr_idx)  - M*se.activities[j]->source_flags(scr_idx));
                    } else {
                        grbModel.addConstr( se.activities[j]->start_time() - o1.start_time()   + min_time_src <=  
                                        2*M-M*o1.source_flags(scr_idx)  - M*se.activities[j]->source_flags(scr_idx));
                    }
                }
            }
        }
    }
}
    
// ...........................................................................
void Solver::addGeneralObservationConstrains( bool baseModel ){
// ...........................................................................
    
    log<DETAIL>( "*** add general observation constraints");
    
    for(Observation& obs : observations){
        Activity& o1 = *(*stations[obs.station1]).activities[obs.interval_idx];
        Activity& o2 = *(*stations[obs.station2]).activities[obs.interval_idx];
        
        /***********************************************************************
         * if obs.selected is true then  o1 and o2 have to be selected too
         * 
         * obs.slected == 1 --> o1.selected = 1 && o1.selected = 2
         */
        grbModel.addConstr( o1.selected() -1 >= M*(obs.selected-1) );
        grbModel.addConstr( o2.selected() -1 >= M*(obs.selected-1) );
                   
        for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
            
            /***********************************************************************
             * o1.source_flags(scr_idx) && o2.source_flags(scr_idx) --> obs.selected = 1
             */
            grbModel.addConstr(obs.selected - 1 <= 2*M-M*o1.source_flags(scr_idx)-M*o2.source_flags(scr_idx));
            grbModel.addConstr(1 - obs.selected <= 2*M-M*o1.source_flags(scr_idx)-M*o2.source_flags(scr_idx));
            
            /***********************************************************************
             * if obs.selected is true then both stations have to observe the same source
             * 
             * obs.selected --> o1.source_flags(scr_idx) == o2.source_flags(scr_idx)
             * 
             */
            grbModel.addConstr( o1.source_flags(scr_idx) - o2.source_flags(scr_idx) <= M*(1-obs.selected) );
            grbModel.addConstr( o2.source_flags(scr_idx) - o1.source_flags(scr_idx) <= M*(1-obs.selected) );
        }
        
         /***********************************************************************
         * observation of one source starts simultaniously at both stations
         *  start1 <= start2
         *  start1 >= start2   ---> start1 == start2
         */
        if(!baseModel){
            grbModel.addConstr(o1.start_time()-o2.start_time()<=M-M*obs.selected);
            grbModel.addConstr(o2.start_time()-o1.start_time()<=M-M*obs.selected);
        }
       
    }
    
    /***************************************************************************
     * if observation is selected at least one baseline including the station has to be used at the same scan
     * 
     * For each observation
     *  o.selected() <= sum obs.selected
     * 
     * if selected == 0 --> expr >= 0
     * if selected == 1 --> expr >= 1     
     * -> if observation is selected at least one baseline including the station has to be used at the same scan
     * 
     */
    for(size_t sta_idx=0; sta_idx < n_sta; ++sta_idx) {
        StationEnvironment & se = *stations[sta_idx];
        for(size_t k=0; k < numberOfIntervals; ++k){
            GRBLinExpr expr;
            Activity& o = *se.activities[k];
            for(Observation& obs : observations) {
                if(obs.interval_idx == k && (obs.station1==sta_idx || obs.station2==sta_idx)){
                    expr += obs.selected;
                }
            }
            grbModel.addConstr(o.selected()<=expr);
        }
    }
}

// ...........................................................................
void Solver::addSubNetAndSourceRatioConstraints(){
// ...........................................................................
    
    int maxSubNets = (int)(*setup)["SKED"]["MaxSubNets"];
    
    double ratioUsedSources = 0.0;
    if( (*setup)["SKED"].exists("RatioUsedSources") ){
        ratioUsedSources = (double)(*setup)["SKED"]["RatioUsedSources"]; 
    }
    
    if(maxSubNets > 0 || ratioUsedSources > 0){
        for(int i=0; i < numberOfIntervals; ++i){
            for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
                GRBLinExpr sum = 0;
                for(size_t sta_idx = 0; sta_idx < stations.size(); ++sta_idx){
                    grbModel.addConstr( subNetFlag[i][scr_idx] >= (*stations[sta_idx]).activities[i]->source_flags(scr_idx) );
                    sum += (*stations[sta_idx]).activities[i]->source_flags(scr_idx);
                }
                grbModel.addConstr( subNetFlag[i][scr_idx] <= sum );
            }
        }
    }
    
    if(maxSubNets > 0){
        log<DETAIL>( "*** add subnet constraints");
        for(int i=0; i < numberOfIntervals; ++i){
            GRBLinExpr sum = 0;
            for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
                sum += subNetFlag[i][scr_idx];
            }
            grbModel.addConstr( sum <= maxSubNets );
        }
    }

    if(ratioUsedSources > 0){
        log<DETAIL>( "*** add percentage of used sources constraints");
        GRBLinExpr numObservedSources = 0;
        for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
            GRBLinExpr sum = 0;
            for(int i=0; i < numberOfIntervals; ++i){
                sum += subNetFlag[i][scr_idx];
                grbModel.addConstr( sourceObservedAtleastOnce[scr_idx] >= subNetFlag[i][scr_idx] );
            }
            grbModel.addConstr( sourceObservedAtleastOnce[scr_idx] <= sum );
            numObservedSources += sourceObservedAtleastOnce[scr_idx];
        }
        grbModel.addConstr( numObservedSources/n_src >= ratioUsedSources );
    }
    
}

// ...........................................................................
void Solver::addNumSourceConstraints(){
// ...........................................................................
    
    int maxNumSource(0), minNumSource(0);
    
    if( (*setup)["SKED"].exists("MaxNumSource") ){
        maxNumSource = (int)(*setup)["SKED"]["MaxNumSource"];
    }
    if( (*setup)["SKED"].exists("MinNumSource") ){
        minNumSource = (int)(*setup)["SKED"]["MinNumSource"];
    }
    
    if( minNumSource > 0 || maxNumSource > 0 ){
        
        log<DETAIL>( "*** add min/max number of observations per source");
        
        for(std::vector<ivg::Analysis_station*>& stas : session->get_trf_ptr()->get_station_twin_list()){
            for(size_t scr_idx=0; scr_idx < n_src; ++scr_idx){
                GRBLinExpr sum;
                for(ivg::Analysis_station* sta : stas){
                    int sta_idx = sta->get_idx();
                    StationEnvironment & se = *stations[sta_idx];  
                    for(int i=0; i < numberOfIntervals; ++i){  
                        sum += se.activities[i]->source_flags(scr_idx);
                    }

                }
                if( minNumSource > 0 ){
                    grbModel.addConstr( sum >= minNumSource);
                }
                if( maxNumSource > 0 ){
                    grbModel.addConstr( sum <= maxNumSource);
                }
            }
        }
    }
}
// ...........................................................................
void Solver::add_SNR_durationConstraints( bool baseModel ){
// ...........................................................................
    
    log<DETAIL>( "*** add SNR / observation duration constraints");
    
     /***************************************************************************
     * observation duration for all baselines
     */ 
    ivg::Schedule* const s = session->get_schedule_ptr();
    
    double x_thres = (double)(*setup)["SKED"]["snr_min_x"]-(double)(*setup)["SKED"]["snr_margin_x"];
    double s_thres = (double)(*setup)["SKED"]["snr_min_s"]-(double)(*setup)["SKED"]["snr_margin_s"];
    
    // compute observation duration required to reach specified snr
     for(Observation& obs : observations){
        Activity& o1 = *(*stations[obs.station1]).activities[obs.interval_idx];
        Activity& o2 = *(*stations[obs.station2]).activities[obs.interval_idx];

        for(size_t sou_idx = 0; sou_idx < n_src; ++sou_idx){
            
            ivg::Date lowerBoundEpoch = session->getStart();
            lowerBoundEpoch.add_secs( obs.interval_idx*intervalLength );
            ivg::BandInfo xband, sband;
            s->compute_band_info (*session->get_crf_ptr()->get_source( sou_idx ), *session->get_trf_ptr()->get_station(obs.station1),
                                  *session->get_trf_ptr()->get_station(obs.station2),lowerBoundEpoch, xband, sband );
            double duration = s->compute_obs_duration(*session->get_crf_ptr()->get_source( sou_idx ), *session->get_trf_ptr()->get_station(obs.station1),
                                  *session->get_trf_ptr()->get_station(obs.station2),lowerBoundEpoch, xband, sband );                

            
            // if duration is longer than maximal time compute snr
            if ( duration > maxScanLength){
                // compute achieved snr
                double snr_x(0.0), snr_s(0.0);
                s->compute_snr(maxScanLength,xband,sband, snr_x, snr_s);
                               
                if( snr_x >  x_thres && snr_s > s_thres ){
                    if(!baseModel){
                        grbModel.addConstr( o1.obs_duration() - maxScanLength >= -2*M + M*obs.selected + M*o1.source_flags( sou_idx ) );
			//grbModel.addConstr( maxScanLength - o1.obs_duration() >= -2*M + M*obs.selected + M*o1.source_flags( sou_idx ) );
                        grbModel.addConstr( o2.obs_duration() - maxScanLength >= -2*M + M*obs.selected + M*o2.source_flags( sou_idx ) );
			//grbModel.addConstr( maxScanLength - o2.obs_duration() >= -2*M + M*obs.selected + M*o2.source_flags( sou_idx ) );
                    }
                } else {
                    // Do not use source if computed duration is larger than maxScanLength and snr is lower than SNR_min - margin
                    grbModel.addConstr( o1.source_flags( sou_idx ) == 0 );
                    grbModel.addConstr( o2.source_flags( sou_idx ) == 0 );
                }
            } else {
                if(!baseModel){
                    grbModel.addConstr( o1.obs_duration() - duration >= -2*M + M*obs.selected + M*o1.source_flags( sou_idx ) );
		    //grbModel.addConstr( duration - o1.obs_duration() >= -2*M + M*obs.selected + M*o1.source_flags( sou_idx ) );
                    grbModel.addConstr( o2.obs_duration() - duration >= -2*M + M*obs.selected + M*o2.source_flags( sou_idx ) );
		    //grbModel.addConstr( duration - o2.obs_duration() >= -2*M + M*obs.selected + M*o2.source_flags( sou_idx ) );
                }
            }
        }         
    }
    
    if(!baseModel){
        /***************************************************************************
         * scan duration has to be equal at all stations 
         */ 
        for(size_t sta_idx_1 = 0; sta_idx_1 < n_sta; ++sta_idx_1) {
            StationEnvironment& se1 = *stations[sta_idx_1];
            for(size_t sta_idx_2 = sta_idx_1+1; sta_idx_2 < n_sta; ++sta_idx_2){
                StationEnvironment& se2 = *stations[sta_idx_2];  
                for(int k=0; k < numberOfIntervals; ++k){
                    grbModel.addConstr( se1.activities[k]->obs_duration() == se2.activities[k]->obs_duration() );
                }  
            }
        }
    }
}

// ...........................................................................
void Solver::addNumberOfObservationConstraints(){
// ...........................................................................
    

    if( (*setup)["SKED"].exists("minNumObs") && (int)(*setup)["SKED"]["minNumObs"] > 0 ){
        log<DETAIL>( "*** add min number of observation constraints");
        GRBLinExpr numObs;
        for(Observation& obs : observations){
            numObs += obs.selected;
        }
        grbModel.addConstr(numObs >= (int)(*setup)["SKED"]["minNumObs"] );
    }
    if( (*setup)["SKED"].exists("maxNumObs") && (int)(*setup)["SKED"]["maxNumObs"] > 0 ){
        log<DETAIL>( "*** add max number of observation constraints");
        GRBLinExpr numObs;
        for(Observation& obs : observations){
            numObs += obs.selected;
        }
        grbModel.addConstr(numObs <= (int)(*setup)["SKED"]["maxNumObs"] );
    }
}

// ...........................................................................
void Solver::addMinTimeWithoutObsConstraints(){
// ...........................................................................
    
    log<DETAIL>( "*** add minimal time without observation constraints");
    
    for(size_t sta_idx = 0; sta_idx < stations.size(); ++sta_idx){
        StationEnvironment & se = *stations[sta_idx];  
        

        for(int i=0; i < numberOfIntervals; ++i){  
            
            GRBLinExpr sumSelected = 0;
            
//            /*******************************************************************
//             * time between two suceeding selected obs starts has to be smaller than given value
//             * used to limit the number of constraints 
//             * within a window of size max_time_without_observation at least on observation has to take place
//             *  
//             */
//            if( i <= numberOfIntervals - maxIntsWithoutObs ){
//                GRBLinExpr sumObsWindowSelected = 0;
//                for(int j=0; j < maxIntsWithoutObs ; ++j){
//                    sumObsWindowSelected +=  se.activities[i+j].selected();
//                }
//                grbModel.addConstr( sumObsWindowSelected >= 1 );
//            } 
//                     

            Activity& o1 = *se.activities[i];
 
            
            /*******************************************************************
             * time between two suceeding selected obs starts has to be smaller than given value
             * used to limit the number of constraints added for modelling the slew of alignment
             * 
             *  
             */
            sumSelected = 0;
                for(int k=i+1; k < numberOfIntervals; ++k){
                    Activity& o2 = *se.activities[k];
                    
                    //  this expression is zero if and only if o1 and o2 are selected and o2 is the first selected observation after o1. Else this expression is M or even larger
                    GRBLinExpr firstSelectedAfterCurrentSelected = M*( 1-o1.selected() ) + M*( 1-o2.selected() )  + M*sumSelected;
                    
                    grbModel.addConstr( o2.start_time()-o1.start_time()-max_time_without_observation  <=  firstSelectedAfterCurrentSelected  );
                    
                    sumSelected += o2.selected();
                    
            }
          
            // at least one observation has to be within the time period from 0-max_time_without_observation
            if( i == 0){
                GRBLinExpr sumObsWindowSelected = 0;
                for(int j=0; j < maxIntsWithoutObs ; ++j){
                    sumObsWindowSelected +=  se.activities[j]->selected();
                }
                grbModel.addConstr( sumObsWindowSelected >= 1 );
            }
        }
    }
}
    
// ...........................................................................
void Solver::addTwinTelescopeConstraints(){
// ...........................................................................
    
    
     /***************************************************************************
     * TWIN TELESCOPES
     */
    double twin_arc_dist =  (double)(*setup)["SKED"]["twin_arc_dist"]; // maximal spherical distance between two sources observed at the same time by two twins in degree
    
    if(twin_arc_dist > 0.0 ){
        log<DETAIL>( "*** add twin telescope constraints");
        for(lps::Observatory& observatory : observatories){
            unsigned n_sta_obs = observatory.stations.size();
            if( n_sta_obs == 2 ){
                unsigned sta_idx_1 = observatory.stations[0].sta_idx;
                unsigned sta_idx_2 = observatory.stations[1].sta_idx;
                for(size_t sou_idx_1 = 0; sou_idx_1 < n_src; ++sou_idx_1){
                    
                    for(size_t k=0; k < numberOfIntervals; ++k){
                        Activity& o1 = *(*stations[sta_idx_1]).activities[k];
                        Activity& o2 = *(*stations[sta_idx_2]).activities[k];
                        grbModel.addConstr( o2.source_flags(sou_idx_1) <= 1-o1.source_flags(sou_idx_1) );
                    }
                    
                    for(size_t sou_idx_2 = sou_idx_1+1; sou_idx_2 < n_src; ++sou_idx_2){
                        ivg::Source* sou_1 = session->get_crf_ptr()->get_source(sou_idx_1);
                        ivg::Source* sou_2 = session->get_crf_ptr()->get_source(sou_idx_2);
                        double arc = sou_1->calc_arclength_to_source(*sou_2)*ivg::rad2d;
                        if( arc < twin_arc_dist){
                            for(size_t k=0; k < numberOfIntervals; ++k){
                                Activity& o1 = *(*stations[sta_idx_1]).activities[k];
                                Activity& o2 = *(*stations[sta_idx_2]).activities[k];
                                if( session->get_common_transits().sourceIsVisible(sou_idx_1, sta_idx_1, k*intervalLength) && session->get_common_transits().sourceIsVisible(sou_idx_2, sta_idx_2, k*intervalLength)){
                                   /*
                                    *  o1.sou1 == 1 --> o2.sou2 = 0  \ 
                                    *                                  --->  o2.sou2 <= 1 - o1.sou1  <=>  o1.sou1 <= 1 - o1.sou2    
                                    *  o2.sou2 == 1 --> o1.sou1 = 0  /
                                    * 
                                    *
                                    *  o1.sou2 == 1 --> o2.sou1 = 0  \ 
                                    *                                  --->  o2.sou1 <= 1 - o1.sou2  <=>  o1.sou2 <= 1 - o1.sou1    
                                    *  o2.sou1 == 1 --> o2.sou1 = 0  /
                                    */
                                   grbModel.addConstr( o2.source_flags(sou_idx_2) <= 1-o1.source_flags(sou_idx_1) );
                                   grbModel.addConstr( o2.source_flags(sou_idx_1) <= 1-o1.source_flags(sou_idx_2) );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ...........................................................................
void Solver::addSlewConstraints(){
// ...........................................................................  
    
    log<DETAIL>( "*** add telescope slewing and cable wrap constraints");
    
    for(size_t sta_idx = 0; sta_idx < stations.size(); ++sta_idx){
        StationEnvironment & se = *stations[sta_idx];  
        ivg::Antenna antenna = session->get_trf_ptr()->get_station(sta_idx)->get_antenna_info();
        
        for(int i=0; i < numberOfIntervals; ++i){  
            
            Activity& o1 = *se.activities[i];
                  
            /*******************************************************************
             * wrap has always to be within limits
             * 
             */
            grbModel.addConstr( o1.wrap() - o1.slew_align()  - o1.wrap().get(GRB_DoubleAttr_LB) >= -M*(1-o1.selected()) );
            grbModel.addConstr( o1.wrap() - o1.slew_align()  - o1.wrap().get(GRB_DoubleAttr_UB) <= M*(1-o1.selected()) );
            grbModel.addConstr( o1.wrap() + o1.slew_observe()  - o1.wrap().get(GRB_DoubleAttr_UB) <= M*(1-o1.selected()) );
            grbModel.addConstr( o1.wrap() + o1.slew_observe()  - o1.wrap().get(GRB_DoubleAttr_LB) >= -M*(1-o1.selected())  );
            
            
            /*******************************************************************
             * wrap at end of an observation ==  wrap at start of next alignment
             */
            GRBLinExpr sumSelected = 0;
            for(int k=i+1; k < numberOfIntervals && objectIsWithinRange(i,k); ++k){
                Activity& o2 = *se.activities[k];
                // this expression is zero if and only if o1 and o2 are selected and o2 is the first selected observation after o1. Else this expression is M or even larger
                GRBLinExpr firstSelectedAfterCurrentSelected = M*( 1-o1.selected() ) + M*( 1-o2.selected() )  + M*sumSelected;

                grbModel.addConstr( (o1.wrap() + o1.slew_observe()) - (o2.wrap() - o2.slew_align()   )  <=  firstSelectedAfterCurrentSelected  );
                grbModel.addConstr( (o2.wrap() - o2.slew_align()  ) - (o1.wrap() + o1.slew_observe() )  <=  firstSelectedAfterCurrentSelected  );
                
                sumSelected += o2.selected();
                
            }
            
            

            /*******************************************************************
             * slew in azimuth during observation. Required to keep track of cable wrap
             * 
             */
                for(int j=0; j < maxIntsAffectedByObservation(i); ++j){
                    for(int scr_idx = 0; scr_idx < n_src; ++scr_idx){

                        double slew = 0;
                        //find position closest to start of observation but not after the start
                        std::pair<lps::Position, lps::Position> border;
                        if( session->get_common_transits().get_border_pos(scr_idx, sta_idx, i*intervalLength, (i+j+1)*intervalLength, intervalLength, border ) ){


                            double alpha1 = border.second.azimuth() - border.first.azimuth();

                            double alpha2 = -sign(alpha1)*(360-abs(alpha1));

                            // it is assumed the shorter angle is the correct one. Only not satisfied in pathological cases
                            if( abs(alpha1) <= abs(alpha2) ){
                                slew = alpha1;
                            } else {
                                slew = alpha2;
                            }
                                                     
                            grbModel.addConstr( o1.slew_observe()-slew <= 3*M - M*o1.selected() - M*o1.source_flags(scr_idx) - M*o1.end_flags(j) );
                            grbModel.addConstr( slew-o1.slew_observe() <= 3*M - M*o1.selected() - M*o1.source_flags(scr_idx) - M*o1.end_flags(j) );
                            

                        }                    
                    }
                }

/*
            for(int scr_idx = 0; scr_idx < n_src; ++scr_idx){

                //find position closest to start of observation but not after the start
                std::pair<lps::Position, lps::Position> border;
                if( session->get_common_transits().get_border_pos(scr_idx, sta_idx, i*intervalLength, (i+1)*intervalLength, intervalLength, border ) ){

                    double vaz = border.first.dazi();
                    // numeric
                    if(abs(vaz) < 1e-3){
                        if(vaz > 0) 
                            vaz = 1e-3;
                        else
                            vaz = -1e-3;
                    }
                    
                    grbModel.addConstr(  o1.slew_observe() - ( vaz*o1.obs_duration() ) <= 2*M - M*o1.selected() - M*o1.source_flags(scr_idx) );
                    grbModel.addConstr( ( vaz*o1.obs_duration() ) - o1.slew_observe()  <= 2*M - M*o1.selected() - M*o1.source_flags(scr_idx) );
                }
            }
*/
			 

             /*******************************************************************
             * aligment slew
             * 
             */
            sumSelected = 0;
            for(int j=1; i+j < numberOfIntervals && objectIsWithinRange(i,i+j); ++j ){                  
                Activity& o2 = *se.activities[i+j];

                // this expression is zero if and only if o1 and o2 are selected and o2 is the first selected observation after o1. Else this expression is M or even larger
                GRBLinExpr firstSelectedAfterCurrentSelected = M*( 1-o1.selected() ) + M*( 1-o2.selected() )  + M*sumSelected;

                for(size_t scr_idx_1 = 0; scr_idx_1 < n_src; ++scr_idx_1){


                    std::pair<lps::Position, lps::Position> azel_1;
                    if( session->get_common_transits().get_border_pos(scr_idx_1, sta_idx, i*intervalLength, (i+j+1)*intervalLength, intervalLength,azel_1 ) ){
                        for(size_t scr_idx_2 = 0; scr_idx_2 < n_src; ++scr_idx_2){

                            std::pair<lps::Position, lps::Position> azel_2;
                            if( session->get_common_transits().get_border_pos(scr_idx_2, sta_idx, i*intervalLength, (i+j+1)*intervalLength, intervalLength, azel_2 ) ){
                                double alpha1 = azel_2.second.azimuth() - azel_1.first.azimuth();
                                double alpha2 = -sign(alpha1)*(360-abs(alpha1));
                                // to interpret slew_dir as counter clockwise and clockwise, the values of alpha1 and alpha2 must be switch if necessary, such that alpha1 is positiv and alpha2 negativ
                                // not done in this implementation
                                
//                                std::cout << " i:" << i << " i+j:" << i+j << " s1:" << scr_idx_1 << " s2:" << scr_idx_2 << " a1:" << alpha1 << "a2:" << alpha2 << std::endl;


                                /* Constraint slew in AZIMUTH
                                 * 
                                 * slew_observe is already determined with constrains (fixed)
                                 * 
                                 * alpha is azimuth from start of the first observation  to the start of the next observation
                                 * 
                                 * simplified alignment is constrained with
                                 *   o2.slew_align = alpha - o1.slew_observe
                                 * 
                                 * Rotation direction has to be considered
                                 */
                                grbModel.addConstr( (o1.slew_observe() + o2.slew_align()) - o2.slew_dir()*alpha1 - (1-o2.slew_dir())*alpha2 <= 
                                                    firstSelectedAfterCurrentSelected + M*(1-o1.source_flags(scr_idx_1)) + M*(1-o2.source_flags(scr_idx_2))  );
                                grbModel.addConstr( o2.slew_dir()*alpha1 + (1-o2.slew_dir())*alpha2 - (o1.slew_observe() + o2.slew_align()) <= 
                                                    firstSelectedAfterCurrentSelected + M*(1-o1.source_flags(scr_idx_1)) + M*(1-o2.source_flags(scr_idx_2))  );

                                double epsilon = azel_2.second.elevation() - azel_1.first.elevation();

                                /* Constraint slew in ELEVATION
                                 * 
                                 * change in elevation during observation is neglected here. additional angle should be small
                                 * 
                                 */
                                grbModel.addConstr( o2.slew_ele_align() - epsilon <= 
                                                    firstSelectedAfterCurrentSelected + M*(1-o1.source_flags(scr_idx_1)) + M*(1-o2.source_flags(scr_idx_2)) );
                                grbModel.addConstr( epsilon - o2.slew_ele_align() <= 
                                                    firstSelectedAfterCurrentSelected + M*(1-o1.source_flags(scr_idx_1)) + M*(1-o2.source_flags(scr_idx_2)) );


                            }
                        }
                    }                    
                }

                sumSelected += o2.selected();

            }


            /*******************************************************************
             * slew_duration duration
             *  slew_duration >= |  slew_align/slew_rate  |
             * 
             *         y >= | x |
             *            y >=  x
             *            y >= -x
             * 
             * 
             */
            // azimuth
            grbModel.addConstr( o1.slew_duration() - (o1.slew_align()/(antenna.azi_rate*ivg::rad2d) ) >= -M+M*o1.selected() );
            grbModel.addConstr( o1.slew_duration() + (o1.slew_align()/(antenna.azi_rate*ivg::rad2d) ) >= -M+M*o1.selected() );
            // elevation
            grbModel.addConstr( o1.slew_duration() - (o1.slew_ele_align()/(antenna.ele_rate*ivg::rad2d) ) >= -M+M*o1.selected() );
            grbModel.addConstr( o1.slew_duration() + (o1.slew_ele_align()/(antenna.ele_rate*ivg::rad2d) ) >= -M+M*o1.selected() );

        }
    }
}

// ...........................................................................
void Solver::addNeutralPointConstraints(){
// ........................................................................... 
    
    log<DETAIL>( "*** add start at neutral point constraints");

    for(size_t sta_idx = 0; sta_idx < stations.size(); ++sta_idx){
        StationEnvironment & se = *stations[sta_idx];

        double neutralPoint = (se.activities[0]->wrap().get(GRB_DoubleAttr_LB) + se.activities[0]->wrap().get(GRB_DoubleAttr_UB))/2;

        int intpart = 0;
        modf (neutralPoint/360.0 , &intpart);
        double azi_np = neutralPoint - 360*intpart; //azimuth of neutral point (neutral point in range 0...360 )

        GRBLinExpr sumSelected = 0;
        for(int i=0; i < numberOfIntervals && objectIsWithinRange(0,i); ++i){
            Activity& o1 = *se.activities[i];

            grbModel.addConstr( ( o1.wrap() - o1.slew_align() )  - neutralPoint  <= M*(1-o1.selected()) + M*sumSelected );
            grbModel.addConstr( neutralPoint - ( o1.wrap() - o1.slew_align() )   <= M*(1-o1.selected()) + M*sumSelected );

            for(int scr_idx = 0; scr_idx < n_src; ++scr_idx){
                std::pair<lps::Position, lps::Position> border;
                if( session->get_common_transits().get_border_pos(scr_idx, sta_idx, i*intervalLength, (i+1)*intervalLength, intervalLength, border ) ){
                    double slew_az = border.first.azimuth() - azi_np;
                    grbModel.addConstr( o1.slew_align()  - slew_az  <= M*(1-o1.selected()) + M*sumSelected + M*(1-o1.source_flags(scr_idx)) );
                    grbModel.addConstr( slew_az - o1.slew_align()   <= M*(1-o1.selected()) + M*sumSelected + M*(1-o1.source_flags(scr_idx)) );

                    double slew_el = border.first.elevation() - 90.0;
                    grbModel.addConstr( o1.slew_ele_align()  - slew_el  <= M*(1-o1.selected()) + M*sumSelected + M*(1-o1.source_flags(scr_idx)) );
                    grbModel.addConstr( slew_el - o1.slew_ele_align()   <= M*(1-o1.selected()) + M*sumSelected + M*(1-o1.source_flags(scr_idx)) );
                }

            }

            sumSelected += o1.selected();
        }
    }
}

// ...........................................................................
bool Solver::optimize(std::function<void (const std::vector<lps::StationActivity> &)> callback,
                      std::function<void (int, int, const lps::Wedge &, int, bool)> treeCallback)
// ...........................................................................
{

    callbackFunction = callback;
    this->treeCallback = treeCallback;
    grbModel.setCallback(this);


    grbModel.optimize();
    
    int optimstatus = grbModel.get(GRB_IntAttr_Status);
    
    bool success = false;
    
    switch (optimstatus){
        case GRB_OPTIMAL:{
            success = true;
            break;
        }
        case GRB_INF_OR_UNBD:{
            log<WARNING>( "!!! Model is infeasible or unbounded");
            break;
        }
        case GRB_INFEASIBLE:{
            log<WARNING>( "!!! Model is infeasible");
            break;
        }
        case GRB_UNBOUNDED:{
            log<WARNING>( "!!! Model is unbounded");
            break;
        }
        case GRB_TIME_LIMIT:{
            log<WARNING>( "!!! Optimization was stopped because time limit has been reached");
            if(  grbModel.get( GRB_IntAttr_SolCount) > 0 ){
                success = true;
            }
            break;
        }
        case GRB_NUMERIC:{
            log<WARNING>( "!!! Optimization was terminated due to unrecoverable numerical difficulties");
            break;
        }
        default:{
            log<WARNING>( "!!! Optimization was stopped with status = ") % optimstatus;
            break;
        }
    }
    

    if (success){
        callback(collectResults());
        
        int max_tree_level = (int) (*session->get_setup())["SKED"]["max_tree_level"];
        std::map<int, std::vector<double> > coverage_absolut;
        std::map<int, std::vector<double> > coverage_relativ;
        
        std::map<int, std::vector<double> > coverage_absolut_surface;
        std::map<int, std::vector<double> > coverage_relativ_surface;
        
        for(lps::Observatory& observatory : observatories){
            
            for(lps::StationEnvironment& se : observatory.stations){
                int i = se.sta_idx;
                coverage_absolut[i] = std::vector<double>(max_tree_level, 0.0);
                coverage_relativ[i] = std::vector<double>(max_tree_level, 0.0);

                coverage_absolut_surface[i] = std::vector<double>(max_tree_level, 0.0);
                coverage_relativ_surface[i] = std::vector<double>(max_tree_level, 0.0);
            }
            
            int t = 0;
            int u = 0;
            for(lps::Node<GRBVar,lps::Wedge>& root: observatory.roots){
                root.visitBoundAndBool([&](int level, const lps::Wedge & rect, bool valid) {
                    
                    if (valid) {    
                        
                        bool cell_contains_observation = observatory.cells[u].get(GRB_DoubleAttr_X) > 0.5;
                        
                        for(lps::StationEnvironment& se : observatory.stations){
                            int i = se.sta_idx;
                            treeCallback(level, t, rect, i, cell_contains_observation);
                            if(cell_contains_observation){
                                coverage_absolut[i][level - 1] += 1.0 / n_segments[level-1];
                                coverage_relativ[i][level - 1] += 1.0 / observatory.max_possible_selected_cells[level - 1];

                                coverage_absolut_surface[i][level - 1] += rect.surfaceArea() / (2*M_PI);
                                coverage_relativ_surface[i][level - 1] += rect.surfaceArea() / observatory.max_possible_surface_area[level - 1];
                            }
                        }
                        
                        

                        u++;
                    }
                });
                ++t;
            }
        }

        bool objective_surface =  (bool)(*session->get_setup())["SKED"]["objective_surface"];
        bool objective_relative =  (bool)(*session->get_setup())["SKED"]["objective_relative"];
        
        if( objective_relative == false && objective_surface == false ){
            std::cout << ivg::Logger::get_color("boldblue");
        }
        std::cout << "absolute coverage --------------------------------------" << std::endl;
        session->show_coverage(coverage_absolut);
        std::cout << ivg::Logger::get_color("white");
        
        if( objective_relative == true && objective_surface == false ){
            std::cout << ivg::Logger::get_color("boldblue");
        }
        std::cout << "relativ coverage ---------------------------------------" << std::endl;
        session->show_coverage(coverage_relativ);
        std::cout << ivg::Logger::get_color("white");
        
        if( objective_relative == false && objective_surface == true ){
            std::cout << ivg::Logger::get_color("boldblue");
        }
        std::cout << "absolute surface area coverage -------------------------" << std::endl;
        session->show_coverage(coverage_absolut_surface);
        std::cout << ivg::Logger::get_color("white");
        
        if( objective_relative == true && objective_surface == true ){
            std::cout << ivg::Logger::get_color("boldblue");
        }
        std::cout << "relativ surface area coverage --------------------------" << std::endl;
        session->show_coverage(coverage_relativ_surface);
        std::cout << ivg::Logger::get_color("white");
        std::cout << "--------------------------------------------------------" << std::endl;
    }
    return success;

}

// ...........................................................................
lps::StationActivity Solver::createActivity(const lps::Activity & a, int station){
// ...........................................................................
    lps::StationActivity sa;
    
    sa.begin = a.sol_begin();
    sa.end   = a.sol_end();
    sa.source_idx = a.getActiveSource();
    sa.station_idx = station;
   
    sa.index = a.index();
    
    ivg::Date s,e;
    s =  session->getStart();
    e =  session->getStart();
    s.add_secs(a.sol_begin());
    e.add_secs(a.sol_end());
    sa.start_epoch = s;
    sa.end_epoch = e;
    
    if(modelAlignment)
        sa.wrap = a.sol_wrap();
    else{
        sa.wrap = 0.0;
    }
    sa.slew_time = a.sol_slew_duration();
    sa.align_duration = a.sol_align_duration();
    sa.obs_duration = a.sol_obs_duration();
    
    return sa;
}

// ...........................................................................
std::vector<lps::StationActivity> Solver::collectResults(){
// ...........................................................................
    std::vector<lps::StationActivity> results(0);

    for(Observation& obs : observations){
        if(obs.selected.get(GRB_DoubleAttr_X) > 0.0001){
        StationEnvironment &se1 = *stations[obs.station1];
        StationEnvironment &se2 = *stations[obs.station2];
            results.push_back(createActivity(*se1.activities[obs.interval_idx], obs.station1));
            results.push_back(createActivity(*se2.activities[obs.interval_idx], obs.station2));
        }
    }

    return results;
}

// ...........................................................................
lps::StationActivity Solver::createTempActivity(const lps::Activity & a, int station){
// ...........................................................................
    lps::StationActivity sa;

    sa.begin = a.sol_begin(true);
    sa.end   = a.sol_end(true);
    sa.source_idx = -1;
    for(size_t i=0; i < n_src; ++i){
        if(getSolution(a.source_flags(i)) > 0.001){
            sa.source_idx = i;
            break;
        }
    }
    sa.station_idx = station;
    
    sa.obs_duration = a.sol_obs_duration(true);
            
    sa.index = a.index();
    return sa;
}

// ...........................................................................
std::vector<lps::StationActivity> Solver::collectTemporaryResults(){
// ...........................................................................
    std::vector<lps::StationActivity> results;

    for(Observation & obs : observations){
        if(getSolution(obs.selected) > 0.0001){
        StationEnvironment &se1 = *stations[obs.station1];
        StationEnvironment &se2 = *stations[obs.station2];
        
            results.push_back(createTempActivity(*se1.activities[obs.interval_idx],obs.station1));
            results.push_back(createTempActivity(*se2.activities[obs.interval_idx],obs.station2));
        }
    }

    return results;
}

// ...........................................................................
void assert_msg(bool value,const std::string & msg){
// ...........................................................................
    if(!value){
        std::cout << msg << std::endl;
        assert(false);
    }
}

// ...........................................................................
void Solver::checkResults(std::vector<lps::StationActivity> & activity){
// ...........................................................................
    log<INFO> ("*** Check Results ");
    
    using Observations = std::set<lps::StationActivity,std::function<bool (const lps::StationActivity &a1, const lps::StationActivity &a2)>>;
    std::vector<Observations> observations;
    std::vector<std::vector<lps::StationActivity>> stationActivity;
    stationActivity.resize(2);
    for(size_t i=0; i < stations.size(); ++i){
        observations.push_back(Observations([](const lps::StationActivity & a1, const lps::StationActivity &a2){
            return a1.source_idx < a2.source_idx;
        }));

    }
    for(lps::StationActivity & a : activity){
      
            Observations & obs = observations[a.station_idx];
            assert_msg(obs.find(a)==obs.end(),"source is observed several times by the same stations");
            obs.insert(a);
        
        stationActivity[a.station_idx].push_back(a);
    }
    for(std::vector<lps::StationActivity> & sa : stationActivity){

        std::sort(sa.begin(),sa.end(),[](const lps::StationActivity& a1, const lps::StationActivity& a2){
           return a1.begin < a2.begin;
        });
        for(size_t i=1; i<sa.size(); ++i){
            assert_msg(sa[i-1].end <= sa[i].begin,"Found overlap");
        }
    }


}

// ...........................................................................
void Solver::printActivityTable() const{
// ...........................................................................
    
// time table *********************************************************************************************** 
   
    if( (int)(*setup)["SKED"]["approach"] == 4){
        // header
        std::cout << "number of intervals:" << numberOfIntervals << std::endl;

        std::cout << "station column info:" << std::endl
                  << "( 1) selected" << std::endl
                  << "( 2) alignment start time" << std::endl
                  << "( 3) observation start time" << std::endl
                  << "( 4) end time" << std::endl
                  << "( 5) alignment start interval" << std::endl
                  << "( 6) interval" << std::endl
                  << "( 7) observation end interval" << std::endl
                  << "( 8) alignment duration" << std::endl
                  << "( 9) slew duration" << std::endl
                  << "(10) observation duration" << std::endl
                  << "(11) source index" << std::endl;

        std::cout << "               ";
        for(size_t sta_idx = 0; sta_idx < session->get_trf_ptr()->get_number_stations(); ++sta_idx) {
            std::cout << " | " << left << setw(65) << session->get_trf_ptr()->get_station(sta_idx)->get_name(ivg::staname::ivs_name);
        }
        std::cout << std::endl;

        for(int k=0; k < numberOfIntervals; ++k){
            std::cout << setw(4) << k;

                std::cout << std::setw(5) << (k)*(int)this->intervalLength << " " << std::setw(5) << (k+1)*(int)this->intervalLength;

                for(size_t sta_idx = 0; sta_idx < session->get_trf_ptr()->get_number_stations(); ++sta_idx) {
                    std::cout  << " | " << (*stations[sta_idx]).activities[k]->times2String();
                }
                std::cout << std::endl;

        }
        std::cout << "-------------------------------------------------------------" << std::endl;

    // wrap table ***********************************************************************************************
        if(modelAlignment){
            std::cout << "station column info:" << std::endl
                      << "(1) selected" << std::endl
                      << "(2) wrap at start of alignment" << std::endl
                      << "(3) wrap at start of observation" << std::endl
                      << "(4) wrap at end of observation" << std::endl
                      << "(5) slew alignment" << std::endl
                      << "(6) slew observation " << std::endl
                      << "(7) slew elevation " << std::endl
                      << "(8) source index" << std::endl;

            std::cout << "               ";
            for(size_t sta_idx = 0; sta_idx < session->get_trf_ptr()->get_number_stations(); ++sta_idx) {
                std::cout << " | " << left << setw(66) << session->get_trf_ptr()->get_station(sta_idx)->get_name(ivg::staname::ivs_name);
            }
            std::cout << std::endl;

            for(int k=0; k < numberOfIntervals; ++k){
                std::cout << setw(4) << k;

                    std::cout << std::setw(5) << (k)*(int)this->intervalLength << " " << std::setw(5) << (k+1)*(int)this->intervalLength;

                    for(size_t sta_idx = 0; sta_idx < session->get_trf_ptr()->get_number_stations(); ++sta_idx) {
                        std::cout  << " | " << (*stations[sta_idx]).activities[k]->wrap2String();
                    }
                    std::cout << std::endl;

            }
            std::cout << "-------------------------------------------------------------" << std::endl;
        }
    } else {
        
        std::cout << "               ";
        for(size_t sta_idx = 0; sta_idx < session->get_trf_ptr()->get_number_stations(); ++sta_idx) {
            std::cout << " | " << left << setw(9) << session->get_trf_ptr()->get_station(sta_idx)->get_name(ivg::staname::ivs_name);
        }
        std::cout << std::endl;
        
        for(int k=0; k < numberOfIntervals; ++k){
            std::cout << setw(4) << k;

                std::cout << std::setw(5) << (k)*(int)this->intervalLength << " " << std::setw(5) << (k+1)*(int)this->intervalLength;

                for(size_t sta_idx = 0; sta_idx < session->get_trf_ptr()->get_number_stations(); ++sta_idx) {
                    std::cout  << " | " << (*stations[sta_idx]).activities[k]->times2String();
                }
                std::cout << std::endl;
        }
    }
    
//    for(int s = 0; s < n_src; ++s){
//        std::cout << sourceObservedAtleastOnce[s].get(GRB_DoubleAttr_X) << std::endl;
//    }
//    
//    for(int k=0; k < numberOfIntervals; ++k){
//        for(int src_idx = 0; src_idx < n_src; ++src_idx){
//            std::cout << subNetFlag[k][src_idx].get(GRB_DoubleAttr_X) << " ";
//        }
//        std::cout << std::endl;
//    }
    
    
}


void Solver::printTimeStats( const std::vector<lps::StationActivity>& activities) const{
    
    lps::Seconds sessdur = intervalLength*numberOfIntervals;
    
    lps::Seconds min_slew_time =  (lps::Seconds)(*session->get_setup())["SKED"]["const_setup"];
    
    std::vector<lps::Seconds> obs_time (n_sta, 0.0);
    std::vector<lps::Seconds> slew_time (n_sta, 0.0);
    std::vector<lps::Seconds> idle_time (n_sta, sessdur);
    std::vector<lps::Seconds> const_time (n_sta, 0.0);
    
    // Const part
    std::vector<lps::Seconds> const_part (n_sta, 0.0);
    lps::Seconds constAlignTime =  (lps::Seconds)(*session->get_setup())["SKED"]["const_idle"]
                                 + (lps::Seconds)(*session->get_setup())["SKED"]["const_source"]
                                 + (lps::Seconds)(*session->get_setup())["SKED"]["const_calib"];
    
    for(unsigned int i =  0; i < n_sta; ++i){
        ivg::Antenna antenna = session->get_trf_ptr()->get_station(i)->get_antenna_info();
        const_part[i] = ( max(antenna.azi_const, antenna.ele_const) + constAlignTime);
    }
    
    for(const lps::StationActivity& sa : activities ){
        lps::Seconds obs_dur = (sa.end-sa.begin);
        obs_time[sa.station_idx] += obs_dur;
        slew_time[sa.station_idx] += sa.slew_time;
        const_time[sa.station_idx] += const_part[sa.station_idx];
                
        idle_time[sa.station_idx] -= (obs_dur + max( sa.slew_time + const_part[sa.station_idx], min_slew_time) );
    }
    
    std::cout << "------ activity statistics -------" << std::endl;
    std::cout << "Station  obs           slew         const       idle         sum" << std::endl;
    for(unsigned int i =  0; i < n_sta; ++i){

        std::cout <<        setw(8) << right << (session->get_trf_ptr()->get_station(i)->get_name(ivg::staname::ivs_name)) << " "
                  <<        setw(5) << right << fixed << setprecision(0) << obs_time[i] << " "
                  << "(" << setw(3) << right << fixed << setprecision(0) << obs_time[i]/sessdur*100 << "%) "
                  <<        setw(5) << right << fixed << setprecision(0) << slew_time[i] << " "
                  << "(" << setw(3) << right << fixed << setprecision(0) << slew_time[i]/sessdur*100 << "%) "
                  <<        setw(5) << right << fixed << setprecision(0) << const_time[i] << " "
                  << "(" << setw(3) << right << fixed << setprecision(0) << const_time[i]/sessdur*100 << "%) "
                  <<        setw(5) << right << fixed << setprecision(0) << idle_time[i] << " "
                  << "(" << setw(3) << right << fixed << setprecision(0) << idle_time[i]/sessdur*100 << "%) "
                  <<        setw(5) << right << fixed << setprecision(0) << obs_time[i]+idle_time[i]+slew_time[i]+const_time[i] << std::endl;
    }
    
}


// ...........................................................................
void Solver::milp2session(){
// ...........................................................................
    
    Setting* const setup = session->get_setup();
    
    // we need this information for the scan meterology and simulation later on
    string ext_data_type = (*setup)["troposphere"]["external_meteo_data"][1];
    string ext_met_data = (*setup)["troposphere"]["external_meteo_data"][2];
    string gpt2filename = (*setup)["troposphere"]["gpt2_grid_file"];
    
    std::vector<lps::StationActivity> res = collectResults();
        
    printTimeStats(res);
    
    log<DETAIL> ("*** extracting milp results");
    
    // key: interval index, value: corresponding entries in std::vector<lps::StationActivity> res
    std::map<unsigned, std::vector<unsigned>> obsMap;
    
    for(unsigned int i=0; i < res.size(); i+=2){
//        std::cout << res[i].index << " " << res[i].begin - offset << " " << res[i].station_idx << "-" << res[i+1].station_idx << " " << res[i].source_idx << std::endl;
        obsMap[ res[i].index ].push_back(i);
    }
    
    std::vector<ivg::Scan>* const scans = session->get_scan_ptr();
    scans->clear();
    
    for( unsigned int interval_idx = 0; interval_idx < numberOfIntervals; ++interval_idx){
        if ( obsMap.count(interval_idx) == 1){
            
            // find all scans
            // key: source index, value: corresponding entries in std::vector<lps::StationActivity> res
            std::map<unsigned, std::vector<unsigned>> scanMap;
            for( unsigned int j = 0; j < obsMap[interval_idx].size(); ++j){
                scanMap[  res[obsMap[interval_idx][j]].source_idx ].push_back( obsMap[interval_idx][j] );
            }
            
            // kv first is source idx
            // kv second is vector with idx of observation with same sources starting at same interval in  std::vector<lps::StationActivity> res
            for( const auto& kv : scanMap ){

                ivg::Source* const scr = session->get_crf_ptr()->get_source(kv.first);
                
                lps::StationActivity sa = res[ kv.second[0] ];
                ivg::Date epoch = sa.start_epoch;
                
                ivg::Partials_t2c tmp { ivg::Matrix( 3,3,0.0 ) };
                ivg::Partials_t2c* deriv_ptr = &tmp;
                ivg::Matrix crf2trf = session->get_eops().form_crf2trf( epoch, true, deriv_ptr );

                ivg::Scan scan(epoch, scr, crf2trf.transpose(), deriv_ptr);

                Seconds duration = sa.obs_duration;

                scan.set_schedulded_duration(round(duration));
                                
                // add observations to scan
                for( const unsigned& obs_idx: kv.second ){
                    
                    // add scan attributes

                    ivg::Analysis_station* const sta1 = session->get_trf_ptr()->get_station( res[obs_idx].station_idx );
                    ivg::Analysis_station* const sta2 = session->get_trf_ptr()->get_station( res[obs_idx+1].station_idx );

                    int sta1_idx = scan.add_sta( sta1 );
                    int sta2_idx = scan.add_sta( sta2 );
                    
                    

                    // used to check results. can be removed in final version
                    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                    
//                    std::cout << "adding observation (" << obs_idx << ") " << sta1->get_name(ivg::staname::ivs_name) 
//                                                                   << " - " << sta2->get_name(ivg::staname::ivs_name) << std::endl;
                    
                    bool isValid = false;
                    for(const lps::Transit& transit : session->get_common_transits().transit(kv.first, sta1->get_idx())){
                        if(     sa.start_epoch.get_double_mjd()+1E-5 > transit.get_startEpoch().get_double_mjd()
                            &&  sa.end_epoch.get_double_mjd()-1E-5 < transit.get_endEpoch().get_double_mjd()){
                            isValid = true;
                        } 
                    }
                    if( !isValid){
                        log<WARNING> ("!!! source: ") % scr->get_name(ivg::srcname::ivs) %"("% kv.first % ") not visible from station: "
                                                      % sta1->get_name(ivg::staname::ivs_name) % "(" %sta1->get_idx() %") at " % res[obs_idx].index; 
                    }
                    isValid = false;
                    for(const lps::Transit& transit : session->get_common_transits().transit(kv.first, sta2->get_idx())){
                        if(     sa.start_epoch.get_double_mjd()+1E-5 > transit.get_startEpoch().get_double_mjd()
                            &&  sa.end_epoch.get_double_mjd()-1E-5 < transit.get_endEpoch().get_double_mjd()){
                            isValid = true;
                        } 
                    }
                    if( !isValid){
                        log<WARNING> ("!!! source: ") % scr->get_name(ivg::srcname::ivs) %"("% kv.first % ") not visible from station: "
                                                      % sta2->get_name(ivg::staname::ivs_name) % "(" %sta2->get_idx() %") at " % res[obs_idx].index; 
                    }
                    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                    scan.add_scan_meteorology(sta1_idx, -999.000, -999.000, -99900.000, 0, ext_data_type, ext_met_data, gpt2filename );
                    scan.add_scan_meteorology(sta2_idx, -999.000, -999.000, -99900.000, 0, ext_data_type, ext_met_data, gpt2filename );

                    // set new wraps for this source in this scan for these two stations
                    scan.set_cable_wrap(sta1, res[obs_idx].wrap*M_PI/180.0 );
                    scan.set_cable_wrap(sta2, res[obs_idx+1].wrap*M_PI/180.0);
    //
    //                // save aiming position of telescopes at the end of the scan
    //                // -> later on we need this for skyplots
    //                scan.set_aiming(bl.first, azel_sta1_end);
    //                scan.set_aiming(bl.second, azel_sta2_end);

                    // add obs attributes

                    ivg::Obs obs_new( session, &scan, sta1_idx, sta2_idx );

                    ivg::BandInfo xband, sband;
                    session->get_schedule_ptr()->compute_band_info (*scr, *sta1, *sta2, epoch, xband, sband );

                    ivg::Matrix azel_sta1 = sta1->calc_az_el(epoch, *scr);
                    ivg::Matrix azel_sta2 = sta2->calc_az_el(epoch, *scr);

                    double snr_x, snr_s;
                    double sigma_final = session->get_schedule_ptr()->compute_sigma_snr(duration, *sta1, *sta2, azel_sta1(1), azel_sta2(1), epoch,
                                                                      xband, sband, snr_x, snr_s);

                    // set calculated SNR for this observation
                    obs_new.set_snr( snr_x, snr_s );


                    // set fake delay and final sigma
                    obs_new.set_delay( 0.0, sigma_final, 0.0 );

                    // further features have to be set for calc_delay in case of impact factors
                    obs_new.set_cable_cal(0.0,0.0);
                    obs_new.set_ion_corr( 0.0 , 0.0, 0.0 , 0.0 );

                    obs_new.set_use_flag( true );
                    scan.add_obs( obs_new );

                    session->obsCounterPlus(1);
                    
                }
                scans->push_back(scan);
                
            }  
        }
    }
     
    session->get_schedule_ptr()->show_scans(*scans, "FINAL SCHEDULE");
    
    // determine final end of session
    ivg::Date end = scans->back().get_epoch();
    ivg::Date start = scans->begin()->get_epoch();
    end.add_secs(scans->back().get_scheduled_duration());
    session->set_end(end);
    session->set_start(start);

    // set all unused sources within _crf to use_me = false
    session->find_and_mark_unused_sources();

    // initialize _param_list within session for completeness
    session->init_param_list();
    
}

// ...........................................................................
void Solver::callback()
// ...........................................................................
{
    try {
        if (where == GRB_CB_MIPSOL){
            
            int max_tree_level = (int)(*session->get_setup())["SKED"]["max_tree_level"];
            std::map<int, std::vector<double> > coverage_absolut;
            std::map<int, std::vector<double> > coverage_relativ;
            
            std::map<int, std::vector<double> > coverage_absolut_surface;
            std::map<int, std::vector<double> > coverage_relativ_surface;
            
            std::vector<lps::StationActivity> activity = collectTemporaryResults();
            callbackFunction(activity);

            for(lps::Observatory& observatory : observatories){

                for(lps::StationEnvironment& se : observatory.stations){
                    int sta_idx = se.sta_idx;
                    coverage_absolut[sta_idx] = std::vector<double>(max_tree_level, 0.0);
                    coverage_relativ[sta_idx] = std::vector<double>(max_tree_level, 0.0);
                    coverage_absolut_surface[sta_idx] = std::vector<double>(max_tree_level, 0.0);
                    coverage_relativ_surface[sta_idx] = std::vector<double>(max_tree_level, 0.0);
                }
                
                int t=0; 
                int u=0;
                for(lps::Node<GRBVar,lps::Wedge>& root: observatory.roots){
                    root.visitBoundAndBool([&](int level, const lps::Wedge & rect, bool valid){
                        if(valid){  
                            if(getSolution(observatory.cells[u]) > 0.5){
                                
                                for(lps::StationEnvironment& se : observatory.stations){
                                    int sta_idx = se.sta_idx;

                                    treeCallback(level, t, rect,sta_idx, valid);
                                    coverage_absolut[sta_idx][level-1] += 1.0/n_segments[level-1];
                                    coverage_relativ[sta_idx][level - 1] += 1.0 / observatory.max_possible_selected_cells[level - 1];

                                    coverage_absolut_surface[sta_idx][level - 1] += rect.surfaceArea() / (2*M_PI);
                                    coverage_relativ_surface[sta_idx][level - 1] += rect.surfaceArea() / observatory.max_possible_surface_area[level - 1];
                                }
                            }
                            u++;
                        }
                    });
                    ++t;
                }
            }
            
            bool objective_surface =  (bool)(*session->get_setup())["SKED"]["objective_surface"];
            bool objective_relative =  (bool)(*session->get_setup())["SKED"]["objective_relative"];

            if( objective_relative == false && objective_surface == false ){
                std::cout << ivg::Logger::get_color("boldblue");
            }
            std::cout << "absolute coverage --------------------------------------" << std::endl;
            session->show_coverage(coverage_absolut);
            std::cout << ivg::Logger::get_color("white");

            if( objective_relative == true && objective_surface == false ){
                std::cout << ivg::Logger::get_color("boldblue");
            }
            std::cout << "relativ coverage ---------------------------------------" << std::endl;
            session->show_coverage(coverage_relativ);
            std::cout << ivg::Logger::get_color("white");

            if( objective_relative == false && objective_surface == true ){
                std::cout << ivg::Logger::get_color("boldblue");
            }
            std::cout << "absolute surface area coverage -------------------------" << std::endl;
            session->show_coverage(coverage_absolut_surface);
            std::cout << ivg::Logger::get_color("white");

            if( objective_relative == true && objective_surface == true ){
                std::cout << ivg::Logger::get_color("boldblue");
            }
            std::cout << "relativ surface area coverage --------------------------" << std::endl;
            session->show_coverage(coverage_relativ_surface);
            std::cout << ivg::Logger::get_color("white");
            std::cout << "--------------------------------------------------------" << std::endl;
           
        }
    }catch(...){
        std::cerr << "ERROR in callback" << std::endl;
    }

}


} //namespace
