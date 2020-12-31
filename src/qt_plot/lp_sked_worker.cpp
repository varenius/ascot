#include "lp_sked_worker.h"

namespace lps{

Worker::Worker(QObject *parent) : QObject(parent)
{

}



void Worker::process(){

    try{
        
        GRBEnv env;     
        
        // parse and set GUROBI options
        for( int param_idx = 0; param_idx < (int)(*session->get_setup())["SKED"]["gurobi_param"].getLength(); ++param_idx){
            std::string param = (const char *)(*session->get_setup())["SKED"]["gurobi_param"][param_idx][0];
            std::string value = (const char *)(*session->get_setup())["SKED"]["gurobi_param"][param_idx][1];
            env.set(param, value);
        }
        
        // build the mixed integer LP  model
        lps::Solver solver(session, env);
        solver.createVariables();
        solver.createObjective();
        solver.createConstraints();
     
        std::string outdir = (const char*)(*session->get_setup())[ "outdir" ];
        
        if( (*session->get_setup())["SKED"].exists("grb_export_model") && (bool)(*session->get_setup())["SKED"]["grb_export_model"]  ){
            solver.getGRBModel().write( outdir + "/" + session->get_name() + ".rew" );
            solver.getGRBModel().write( outdir + "/" + session->get_name() + ".mps" );
        }
   
        if( (*session->get_setup())["SKED"].exists("grbtune") && (bool)(*session->get_setup())["SKED"]["grbtune"]  ){
            // Set the TuneResults parameter to 1
            solver.getGRBModel().set(GRB_IntParam_TuneResults, 1);
            // Tune the model
            solver.getGRBModel().tune();
            // Get the number of tuning results
            int resultcount =  solver.getGRBModel().get(GRB_IntAttr_TuneResultCount);
            if (resultcount > 0){
                // Load the tuned  parameters into the model's environment
                solver.getGRBModel().getTuneResult(0);
                // save tuned parameters
                solver.getGRBModel().write(outdir + "/" + session->get_name() + "_tune.prm");
            }
        } else if ( (*session->get_setup())["SKED"].exists("tune_file") ) {
            std::string tune_file = (const char*)(*session->get_setup())["SKED"]["tune_file"];
            if(!tune_file.empty()){
                solver.getGRBModel().read(tune_file);
            }
        }
        
        if( (*session->get_setup())["SKED"].exists("import_start_solution") ) {
            std::string start_solution = (const char*)(*session->get_setup())["SKED"]["import_start_solution"];
            if(!start_solution.empty()){
                solver.getGRBModel().read(start_solution);
            }
        }
        
        bool success = solver.optimize([&]( const std::vector<lps::StationActivity>& activity ){
            
            std::map<unsigned, std::vector<unsigned>> obsMap;
            for(unsigned int i=0; i < activity.size(); i+=2){
                obsMap[ activity[i].index ].push_back(i);
            }
            
            int numScans = 0;
            
            for( const auto& obs: obsMap ){
                std::vector<int> scr;
                for( unsigned int j = 0; j < obs.second.size(); ++j){
                    scr.push_back( activity[obs.second[j]].source_idx );
                }
                std::sort(scr.begin(), scr.end());
                std::vector<int>::iterator it = std::unique(scr.begin(), scr.end());
                numScans += std::distance(scr.begin(), it);
            }
            
            std::cout << ivg::Logger::get_color("boldblue") << "FOUND SOLUTION WITH OBSERVATIONS:" << activity.size()/2 << " SCANS:"
                      << numScans << " INTERVALS:" << obsMap.size() << ivg::Logger::get_color("white") << std::endl;
            emit foundSolution(activity);           
        },
        [&](int level, int temporal_idx, const lps::Wedge & rect, int station, bool valid){
            emit selectedCell(level, temporal_idx, rect.toPath(),station, valid);
        });
        
	if ( success ){
	        solver.printActivityTable();
        	solver.milp2session();
	        std::cerr << "worker finished" << std::endl;
      
	        if( (*session->get_setup())["SKED"].exists("grb_export_solution") && (bool)(*session->get_setup())["SKED"]["grb_export_solution"]  ){
        	    solver.getGRBModel().write( outdir + "/" + session->get_name() + ".sol" );
        	}	
	}
      
    }catch(GRBException & e){
        std::cout << e.getMessage() << std::endl;
        return;
    }

    QThread::currentThread()->terminate();

}

} //namespace
