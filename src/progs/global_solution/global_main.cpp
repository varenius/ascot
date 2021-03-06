
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <iterator>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <climits>
#include "matrix.h"
#include "date.h"
#include "lapack_wrapper.h"
#include "tictoc.h"
#include "session.h"
#include "scan.h"
#include "ivg_const.h"
#include "ls_neq.h"
#include "logger.h"
#include "session_inout.h"
#include "masterfile.h"
#include <tclap/CmdLine.h>

#include <curl/curl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <cstdlib>
#include <libconfig.h++>
#include <QApplication>
//#include "plot.h"

using namespace libconfig;

loglevel g_verbose;
void * g_ephis;

struct handling
{
    string type;// g for geodetic
    string pm; //  keep_none
    string ut1; // keep_none
    string nut; // keep_fix
};

// ...........................................................................
 map<string, handling> read_handling_file(string path, string default_eop_handling)
// ...........................................................................
{
     map<string, handling> sesshand;
    
    ifstream inStream;
    inStream.open(path.c_str(), ios::in);
    if( !inStream.is_open() )
        throw runtime_error( "map<string, handling> session_handling = read_handling_file(string path): Failed to open file: " + path );
    else
    {
        string line;
        while( getline(inStream, line))
        {
            vector<string> tokens = get_tokens(line);
            
            string sess = tokens.at(0);
            sesshand[sess].type = tokens.at(1);
            sesshand[sess].pm = default_eop_handling;
            sesshand[sess].ut1 = default_eop_handling;
            sesshand[sess].nut = default_eop_handling;
            
            string config = tokens.at(2); // e.g. YYYYNN
            if(config.substr(0,2) == "YY")
                sesshand[sess].pm = "keep_fix";
            if(config.substr(2,2) == "YY")
                sesshand[sess].ut1 = "keep_fix";
            if(config.substr(4,2) == "YY")
                sesshand[sess].nut = "keep_fix";
        }
    }
     
    return sesshand;
}

int main( int argc, char *argv[] )
{ 
    // .........................................................................
    // handle commandline arguments
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
    
    // name of config file
    TCLAP::ValueArg<std::string> m1Arg( "c","config","Name of the control file",true,"","string");
    cmd.add( m1Arg );
    
    // verbose level
    TCLAP::ValueArg<int> m2Arg( "v","verbose","Verbose Level (0=NOTHING, 1=INFO, 2=DETAIL, 3=RESULT, 4=WARNING, 5=ALL)",false,0,"int");
    cmd.add( m2Arg );
    
    // plot and export residuals
    TCLAP::SwitchArg  resid_switch( "r","residuals","plot and save residuals", false );
    cmd.add( resid_switch );
    
    // parse commandline arguments
    cmd.parse( argc, argv );
    string controlfile = m1Arg.getValue();
    g_verbose = (loglevel) m2Arg.getValue();
    bool resid = resid_switch.getValue();
    
    // .........................................................................
    // initialize timer    
    tictoc timeCount;
    timeCount.tic();
    
    std::cout << " \n =================================================================== " << std::endl;
    std::cout << "     > > > > > BAKKARI GLOBAL SNX SOLUTION < < < < <     " << std::endl << std::endl;

    // loading controlfile
    Config cfg;
    try
    {
       cfg.readFile( controlfile.c_str() );
    }
    catch( libconfig::ParseException & err )
    {
        cerr << "libconfig::" << err.what() << " in " << err.getFile() << " at line " << err.getLine() << endl;
        exit( -1 );
    }
    
    Setting& setup= cfg.lookup( "setup" );
    
    // directory to store output data
    string outdir = setup[ "outdir" ];
    
    // get session_type  
    string session_type = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[1];
    
    // initialize ephemerides only once to save time 
    // and to avoid pointer problems in the assignment-operator of Session!!!
    string ephfile_name = setup["definitions"]["ephemerides"][(const char*)setup["ephemerides"]];
    char nams[400][6];double vals[400];
    void *ephem = jpl_init_ephemeris(ephfile_name.c_str(), nams, vals);
    
    ivg::Masterfile masterfile(setup["definitions"]["masterfiles"], ivg::mastertype::both);
    ivg::Session_inout nonsense("snx",masterfile); // NGS / SNX
    
    // X/KA catalog import - BETA-STATUS
//    ivg::Session XKA; 
//    ivg::Session_inout catanizer("cata");
//    ivg::Matrix Qxx = catanizer.init_xka_catalog(&XKA, &setup, "/data/bakkari/xka.cata");
    
    string default_eop_handling =  "keep_none";
    map<string, handling> session_handling = read_handling_file("/opt/ascot/apriori_files/global.hlf", default_eop_handling);
    try
    {
        cerr << "*** Generating " << setup[ "arc_files" ].getLength() << " different global solutions" << endl;
        
        // create folder for results and errorfiles
        ivg::Date now;
        now.now();
        string folder_name = setup[ "output_folder" ]; 
        folder_name = folder_name + "_" + now.get_date_time("DOY-HH-MI-SS");

        if(chdir(outdir.c_str()) == -1)
          cerr << "!!! ERROR in chdir: " << strerror(errno) << endl; 
        else if(mkdir(folder_name.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
          cerr << "!!! ERROR in mkdir: " << strerror(errno) << endl; 
        else
            cerr << "*** Folder created: " << folder_name << endl;
                
        // while we got new arc files
        int sess_start = 0;
        int sess_end = 0;
        for( int arc=0; arc < setup[ "arc_files" ].getLength(); arc++ )
        {
            string arc_name = setup[ "arc_files" ][ arc ][0];
            int number_sess = setup[ "arc_files" ][ arc ][1];
            
            sess_end += number_sess;
                    
            cerr << "::: Running global solution based on " << arc_name << " arc-file and #" << number_sess << " sessions" <<endl;
            
            ivg::Session_inout sessionizer(session_type,masterfile);
    
            // declaration of session-vector
            vector<ivg::Session> sessions;

            // declare global Session
            ivg::Session G;
            
            // Write errorfile if exceptions fly
            string errorfile = outdir + "/" + folder_name + "/" + arc_name + ".err";
            ofstream outstream ( errorfile.c_str(), ios::out ); 

            // to be able to assign first session as global session
            int sess_cnt = 1;
            // loop over databases
            for( int i_sess=sess_start; i_sess < sess_end; ++i_sess )
            {
                bool session_ok = true;

                string versions = setup[ "sessions" ][ i_sess ][ "version" ];
                string dbname  = setup[ "sessions" ][ i_sess ][ "dbname" ];

                string version;
                istringstream iss(versions);
                int i_version = 0;
                ivg::Date mid_S;
                // loop over analysis centers
                while (std::getline(iss, version, ','))
                {
                    cerr << "--- Using " << dbname << " and " << version << "(" << i_sess-sess_start+1 << "/" << sess_end-sess_start << ")" << endl;

                    // initialize _setup / _name / _date / _ephem / _eops
                    ivg::Session S( &setup, dbname, &ephem, i_sess );

                    try
                    {
                        // load information from each database into session
                        // sets following member variables of the session object S
                        // in case of session_type = snx : _type _start, _end, _trf, _crf, _param_list  
                        sessionizer.load( &S, &setup, dbname, version);   

                        // epoch for EOP transformation based on first AC, used also for all others
                        if(i_version == 0)
                            mid_S = S.get_mid_epoch();
//                            mid_S = S.get_param_list_ptr()->get_param(0)->get_epoch();
			mid_S=ivg::Date(roundf(mid_S.get_double_mjd()*4.0)/4.0);
                        log<INFO>("*** mid epoch of current session ") % mid_S.get_date_time("YY:DOY:SSSSS");

                        ivg::Session T;
                        T = S;
                        
                        // in case of e.g. 93APR20X_ , rename it to 93APR20X
                        if(dbname.substr(8,1) == "_")
                            dbname = dbname.substr(0,8);
                        
                        // e.g. YYYYNN lead to fixing XPO,YPO,
                        if( !session_handling[dbname].type.empty() )
                        {
                            (*T.get_setup())[ "PARAMS" ]["pm"][0]["handling"] = session_handling[dbname].pm;
                            (*T.get_setup())[ "PARAMS" ]["ut1"][0]["handling"] = session_handling[dbname].ut1;
                            (*T.get_setup())[ "PARAMS" ]["nut"][0]["handling"] = session_handling[dbname].nut;
                        }
                        else
                        {
                            (*T.get_setup())[ "PARAMS" ]["pm"][0]["handling"] = default_eop_handling;
                            (*T.get_setup())[ "PARAMS" ]["ut1"][0]["handling"] = default_eop_handling;
                            (*T.get_setup())[ "PARAMS" ]["nut"][0]["handling"] = default_eop_handling;
                        }
//                            throw runtime_error("ERROR: major handling information not existent");
                        
//                        (*T.get_setup())[ "no_net_cnstr" ]["velocities"]["apply"] = false;
//                        (*T.get_setup())["PARAMS"]["stations"][0]["stacking"]["type"] = "independent";
//                        (*T.get_setup())["PARAMS"]["stations"][0]["stacking"]["insert_velocities"]  = false;
//                        
			 (*T.get_setup())["PARAMS"]["stations"][0]["handling"] = "keep_none";
			 (*T.get_setup())["PARAMS"]["stations"][1]["handling"] = "keep_none";
//                        
//                        (*T.get_setup())["PARAMS"]["sources"][0]["sources"] = 1;
//                        (*T.get_setup())["PARAMS"]["sources"][0]["handling"] = "keep_fix";
//                        
//                        (*T.get_setup())["PARAMS"]["sources"][1]["sources"] = 1;
//                        (*T.get_setup())["PARAMS"]["sources"][1]["handling"] = "keep_fix";
                        
                        // set up solution
                        T.init_snx_solution(get_list_element(setup["datadirs"],setup["session_type"])[3], mid_S);
			
                        T.modify_parameterization();
			
			  T.reduce_and_constrain();
			
			// T.solve();
                       
                        // check if station estimates are too big, create output for errorfile
			// vector<string> bad_stations = T.check_stations_estimates();
                        //if(bad_stations.size() > 1)
			  // {
			  //    double percent = (T.get_trf_ptr()->get_number_stations() / bad_stations.size() ) * 100.0;
			    //    stringstream ss;
                            //ss << "WARNING: STA_CHECK " << bad_stations.size() << "/" << T.get_trf_ptr()->get_number_stations() << " stations exceeding threshold: ";
			    // std::copy(bad_stations.begin(), bad_stations.end(),std::ostream_iterator<std::string>(ss," "));
			    // ss << "[UP,EAST,NORTH](m)";
			    //   throw runtime_error(ss.str());
			    //        }
                        
                        
                        // write result files
//                        string resultfile = outdir + "/" + folder_name + "/" + dbname+"_"+version + ".res";
//                        sessionizer.write_results(&T, resultfile );
//                        string snxfile = outdir + "/" + folder_name + "/" + dbname + ".snx";
//                        sessionizer.write_snx(&T, snxfile );
                        
//                        (*S.get_setup())[ "no_net_cnstr" ]["velocities"]["apply"] = false;
//                        (*S.get_setup())["PARAMS"]["stations"][0]["stacking"]["type"] = "global";
//                        (*S.get_setup())["PARAMS"]["stations"][0]["stacking"]["insert_velocities"]  = false;
//                        
                        //(*S.get_setup())["PARAMS"]["stations"][0]["handling"] = "keep_fix";
                        //(*S.get_setup())["PARAMS"]["stations"][1]["handling"] = "keep_fix";
//                        
//                        (*S.get_setup())["PARAMS"]["sources"][0]["sources"] = -6;
//                        (*S.get_setup())["PARAMS"]["sources"][0]["handling"] = "keep_none";
//                        
//                        (*S.get_setup())["PARAMS"]["sources"][1]["sources"] = 6;
//                        (*S.get_setup())["PARAMS"]["sources"][1]["handling"] = "keep_reduce";
//                        
                        S.init_snx_solution(get_list_element(setup["datadirs"],setup["session_type"])[3], mid_S);
		
			S.modify_parameterization();
		
			S.reduce_and_constrain();
		
                        
                    }
                    catch(std::runtime_error& e)
                    {
                        cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                        session_ok = false;
                        
                        stringstream sstmp;
                        vector<string> station_names = S.get_trf_ptr()->get_station_names(ivg::staname::lettercode);
                        copy ( station_names.begin () , station_names.end () , ostream_iterator <string>( sstmp , "" ) );

                        outstream << setw(9) << setfill(' ') << left << dbname << ": " << setw(42) << setfill(' ') << sstmp.str() << " " << e.what() << endl;; 
                        cerr << setw(9) << setfill(' ') << left << dbname << ": " << setw(42) << setfill(' ') << sstmp.str() << " " << e.what() << endl;;

                        //S.show();
                        cerr << "\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                    }


                    if(session_ok)
                    {
                        cerr << "*** Session Nr. " << sess_cnt << " imported" << endl;

                        if(sess_cnt == 1)
                            G = S;
                        else
                            G += S;

                        sess_cnt++;
			
                    }

                i_version++;
                }

            }

            outstream << endl;
            outstream.close();
            
//            (*G.get_setup())["PARAMS"]["stations"][0]["handling"] = "keep_none";
//            (*G.get_setup())["PARAMS"]["stations"][1]["handling"] = "keep_none";
//            
//            (*G.get_setup())[ "no_net_cnstr" ]["stations"]["stations"] = 1;
//            (*G.get_setup())[ "no_net_cnstr" ]["stations"]["apply"] = true;
//            (*G.get_setup())[ "no_net_cnstr" ]["velocities"]["apply"] = false;
//
//            (*G.get_setup())["PARAMS"]["sources"][0]["sources"] = 1;
//            (*G.get_setup())["PARAMS"]["sources"][0]["handling"] = "keep_fix";
//            (*G.get_setup())["PARAMS"]["stations"][0]["stacking"]["type"] = "global";
//            (*G.get_setup())["PARAMS"]["stations"][0]["stacking"]["insert_velocities"]  = true;

            sess_start = sess_end;
            
            int max_iter=0;
            ivg::Param *weak_param = NULL;
            ivg::Param *old_weak_param = NULL;
            double cnstr_sigma = 1e-7; // 1e-7 equal to 20mas
	
	    while(max_iter <= 100)
            {
                try
                {
                    cerr << "///////////////////////////////////////////////////////" << endl;
                    cerr << "::: Estimating global solution" << endl;
                    
                    if(weak_param != NULL)
                    {
                        cerr << "::: Trying to strengthen weak parameter " << weak_param->get_name();
                        if(weak_param->is_type({ivg::paramtype::ra, ivg::paramtype::dec},{0,1}))
                        {
                            vector<int> weakis = G.get_param_list_ptr()->get_indexes({ivg::paramtype::ra, ivg::paramtype::dec}, weak_param->get_name() );
                            //show_vector(weakis);
                            //G.get_param_list_ptr()->get_param(weakis.at(0))->show(); G.get_param_list_ptr()->get_param(weakis.at(1))->show();
                            
                            
                            if(weak_param == old_weak_param)
                            {
                                cnstr_sigma -= (5e-8/10.0);
                                cerr << "!!! Identical weak source. Need to increase constraint to " << cnstr_sigma << endl;
                            }
                            else
                                cnstr_sigma = 5e-8;
                                
                            if(cnstr_sigma<0.0)
                                cnstr_sigma = 5e-012;

                            cerr << " using " << cnstr_sigma << "rad sigma";
                            
                            G.get_param_list_ptr()->get_param(weakis.at(0))->set_offset_cnstr_sigma(cnstr_sigma);
                            G.get_param_list_ptr()->get_param(weakis.at(1))->set_offset_cnstr_sigma(cnstr_sigma);
                        }
                        else if(weak_param->is_type({ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz},{0,1}))
                        { 
                            cerr << endl;
                            string sta_name = weak_param->get_name(); // e.g. HART15M
                            // in case of weak station, constraint all of it
                            ivg::Param_list *list = G.get_param_list_ptr();
                            for(vector<ivg::Param>::iterator param_iter = list->begin(); param_iter!= list->end(); ++param_iter)
                            {
                                if(param_iter->get_name() == sta_name)
                                {
                                    param_iter->show();
                                    if(param_iter->get_order() == 0)
                                        param_iter->set_offset_cnstr_sigma(1.6678e-8);
                                    else
                                        param_iter->set_offset_cnstr_sigma(1.6678e-9);
                                }
                            }
                        }
                        
                        G.get_param_list_ptr()->create_constraint_conditions((*G.get_solution_ptr()));
                    }
		    //std::cout <<"solve"<< endl;
                    // modify and perform least squares solution
                    G.solve();
                    // write result of global solution in SNX
                    string snxfile = outdir + "/" + folder_name + "/" + arc_name + ".snx";
		   
                    sessionizer.write_snx(&G, snxfile );
		    
                    string resultfile = outdir + "/" + folder_name + "/" + arc_name + ".res";
                    sessionizer.write_results(&G, resultfile );
                    cerr << "///////////////////////////////////////////////////////" << endl;
                    
                    // if solving was successfull, we can break the while
                    break;

                }
                catch(std::exception& e)
                {
                    cerr << "!!! ERROR: " << e.what() << endl;

                    // show parameter in case of LAPACK-error
                    std::string tmp = e.what();
                    if(tmp.find("LAPACK")!=std::string::npos)
                    {
                        if( max_iter == 10 )
                            G.get_param_list_ptr()->show();
                        
                        std::stringstream tokenizer(tmp);
                        std::string a, b, c;
                        int pos;
                        tokenizer >> a >> b >> c >> pos;
                        cerr << "!!! " << setfill('0') << setw(5) << right << pos-1 << " " << G.get_param_list_ptr()->get_param(pos-1)->get_resultline(true) << endl;
                        old_weak_param = weak_param;
                        weak_param = G.get_param_list_ptr()->get_param(pos-1);
                    }
                    else
                        break;
                }
                
                max_iter++;
            }
        }
        
    }
    catch( libconfig::SettingException &err )
    {
        cerr << "libconfig::SettingException: " << err.what() << " at " << err.getPath() << endl;
        exit( -1 );
    }
    catch(std::exception& e)
    {
        cerr << "std::exception: " << e.what();
    }
    
    std::cout << "Time required: " << fixed << timeCount.toc() << " seconds." << std::endl;

    return ( 0 );
}
	
	
	
	
	


//            sessionizer.write_time_series( &S,  "/opt/bakkari/output/time_series/" + now.get_date_time("DOY_HH_MI_SS") +, {ivg::paramtype::ra, ivg::paramtype::dec} );   
            
//            ivg::Param_list * param_list_ptr = S.get_param_list_ptr();
//            ivg::Ls_solution * solution_ptr = S.get_solution_ptr();
//            string src = "0119+041"; // 0059+581 , 0923+392
//            
//            int ra_idx = param_list_ptr->get_index(ivg::paramtype::ra, src );
//            int dec_idx = param_list_ptr->get_index(ivg::paramtype::dec, src );
//            
//            if(ra_idx != -1 && dec_idx != -1)
//            {
//                double ra = param_list_ptr->get_param( ra_idx )->get_estimate();
//                double dec = param_list_ptr->get_param( dec_idx )->get_estimate();
//                double epoch = param_list_ptr->get_param( dec_idx )->get_epoch().get_double_mjd();
//                
//                
//              ostringstream ss;
//              ss << "SOLVED: " << setw(9) << setfill(' ') << left << dbname << " "
//              << src << " " << setw(17) << setfill(' ')<< scientific << setprecision(10)  << right << epoch
//              << " " << setw(17) << setfill(' ')<< scientific << setprecision(10)  << right << ra
//              << " " << setw(17) << setfill(' ') << scientific << setprecision(10) << right << dec << endl;  
//              
//              
//               log<SAVE>(" ") % ss.str();
//               
//               outstream << ss.str();
//                
//               ss.clear();
//            }
//            else
//                log<SAVE>(" ") % "\n";
//            
//            
//               
//            // write results in ivg specific format
//            sessionizer.write_results(&S, dir+"/"+dbname+"_result.txt" );
