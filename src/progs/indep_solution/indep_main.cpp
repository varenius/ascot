
/*****************************************************************************
 * independent VLBI solution                                                 *
 * *
 * *
 * *
 * 2014-12-16 - SH                                                           *
 ****************************************************************************/

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
#include "simulation.h"
#include <tclap/CmdLine.h>
#include "masterfile.h"
#include "wrapper.h"

#include <cstdlib>
#include <libconfig.h++>
#include <QApplication>
#include <QColor>
#include "statistics.h"
#include "plot.h"
#include "ascot.h"
#include "tsa.h"
#include "lp_sked_gui.h"
#include "lp_sked_station.h"
#include "lp_sked_worker.h"
#include "transits.h"

#include "milp.cpp"
#include "grid.h"
#include "geometry.h"
#include "definitions.h"

#include "masterfile.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace libconfig;

loglevel g_verbose;
void * g_ephis;

int main( int argc, char *argv[] )
{
    
    
    #ifdef USE_MPI 
    // Initialisiere MPI
    MPI_Init( &argc, &argv ); 
    
    char procName [MPI_MAX_PROCESSOR_NAME];
    
    int mpi_size, mpi_rank, mpi_len;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);//get number processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); //rank/ID of current process
    MPI_Get_processor_name(procName, &mpi_len); //hostname process runs
    
    if(mpi_size < 2)
        std::exit(EXIT_FAILURE);
    
    #endif
   
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
    
    // dbname
    TCLAP::ValueArg<std::string> dbArg("e", "explicitDB", "e.g. 93NOV05XU", false, "", "string");
    cmd.add(dbArg);

    // using interactive ascot wit -I
    TCLAP::SwitchArg  interactive_switch( "I","interactive","interactive ascot", false );
    cmd.add( interactive_switch );

    // parse commandline arguments
    cmd.parse( argc, argv );
    string controlfile = m1Arg.getValue();
    g_verbose = (loglevel) m2Arg.getValue();
    bool resid = resid_switch.getValue();

    //////////////////////////////////////////
    // start the super new cool interactive ascot
    if( interactive_switch.getValue() )
    {
        QApplication a(argc, argv);
        // cfg instance necessary here because copy constructor not supported
        libconfig::Config cfg;
        Ascot main(&cfg, controlfile);
        return a.exec();
    }
    //////////////////////////////////////////
    // .........................................................................
    
    
    #ifdef USE_MPI
    tictoc timeCount;
    if(mpi_rank==0){
        std::cout << " \n =================================================================== " << std::endl;
        std::cout << "     > > > > > ivg::ASCOT (independent solution) " << ivg::ASCOT_version << " < < < < <     " << std::endl << std::endl;
        // initialize timer
        timeCount.tic();
    }
    #else
    // initialize timer
    tictoc timeCount;
    timeCount.tic();
    std::cout << " \n =================================================================== " << std::endl;
    std::cout << "     > > > > > ivg::ASCOT (independent solution) " << ivg::ASCOT_version << " < < < < <     " << std::endl << std::endl;
    #endif


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
   
    // load ephemeris only once to avoid multiple jpl_init_ephemeris
    char nams[400][6];
    double vals[400];
   
    string ephfile_name = setup["definitions"]["ephemerides"][(const char*)setup["ephemerides"]];
     
    void *ephem = jpl_init_ephemeris(ephfile_name.c_str(), nams, vals);

    // initialization of session_inout with defined session_type
     
    string session_type = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[1];
     
    // initializing masterfiles for more information in logfiles
    ivg::Date now;
    now.now();
    ivg::Masterfile masterfile(setup["definitions"]["masterfiles"], ivg::mastertype::both,1979,s2d(now.get_date_time("YYYY")),setup["definitions"]["nscodes"],(const char *)get_list_element((setup)["refframes"],(setup)["trf"])[2],(const char *)get_list_element((setup)["refframes"],(setup)["trf"])[1]);
     
    ivg::Session_inout sessionizer(session_type,masterfile); // NGS / SNX
     
    string outdir = setup[ "outdir" ];
   
    qRegisterMetaType<std::vector<lps::StationActivity>>("std::vector<lps::StationActivity>");
    qRegisterMetaType<lps::Rect>("lps::Rect");
    qRegisterMetaType<lps::Path>("lps::Path");
    
    std::vector<std::string> database;
    std::vector<std::string> versions;
    if (dbArg.isSet()) {
    // -e option : sessions are parsed from argv
        database.push_back(dbArg.getValue());
        versions.push_back("-1");
    } else {
        database.reserve( setup[ "sessions" ].getLength() );
        versions.reserve( setup[ "sessions" ].getLength() );
        for( int i=0; i<setup[ "sessions" ].getLength(); ++i )
        {
            database.push_back( setup[ "sessions" ][ i ][ "dbname" ] );
            versions.push_back( setup[ "sessions" ][ i ][ "version" ] );
        }  
    }        

    
    #ifdef USE_MPI
    stringstream outstream;
    ofstream master_stream;
  
    
    
    string dir = setup[ "logdir" ];
    
    if(mpi_rank==0){
        // write 'bad-session' file
        // name of the file based on the time of execution
        string dir = setup[ "logdir" ];
        ivg::Date now;
        now.now();
        std::string ascot_log = dir+"/"+"ascot_"+now.get_date_time("DOY_HH_MI_SS")+".log";
        master_stream.open( ascot_log.c_str() );

    }
     

    #else
    // write 'bad-session' file
    // name of the file based on the time of execution
    string dir = setup[ "logdir" ];
    //ivg::Date now;
    //now.now();
    std::string ascot_log = dir+"/"+"ascot_"+now.get_date_time("DOY_HH_MI_SS")+".log";
    ofstream outstream( ascot_log.c_str() );
    #endif

  
    if( setup.exists( "SIM" ) && (bool)setup[ "SIM" ]["apply"] 
        && setup[ "SIM" ].exists( "Reference" )
        && (bool)setup[ "SIM" ]["Reference"]["apply"] ){
        
              
        if( (int)(setup)["SIM"]["Reference"]["step"] < 2 ){
            
            ivg::Simulation RefSimulator;
            
            vector<string> stations;
            for( int i=0; i<setup["SKED"]["stations"].getLength(); ++i )
                stations.push_back(setup["SKED"]["stations"][i]);
            std::sort( stations.begin(), stations.end() );  
            
            ivg::Date dummyDate(2018,1);
            ivg::Trf trf( setup, stations, ivg::staname::ivs_name, true, dummyDate.add_days(-10.0), dummyDate.add_days(10.0) );
            
            if( (int)(setup)["SIM"]["Reference"]["step"] == 0 ){
            
                RefSimulator.compute_reference_tropo_triangular_mat( trf,  setup );
                
            } else if ( (int)(setup)["SIM"]["Reference"]["step"] == 1 ){
                double zwd0 = (double)(setup)["SIM"]["Reference"]["troposphere"]["initalEZWD"];
                RefSimulator.simulate_reference_troposhere_and_clock(trf, zwd0, database, setup);
               
            }
            
            exit(0);
            
        }
    }
   
  
    #ifdef USE_MPI

    std::vector<int> last_v, first_v;
        
    unsigned int numberSessions = setup[ "sessions" ].getLength();
        
    unsigned int inc = floor(numberSessions/mpi_size);
    unsigned int rest = numberSessions - inc*mpi_size;
        
    inc;
        
    unsigned int first = 0;
    unsigned int last = inc-1;
    for (unsigned int i = 0 ; i < mpi_size; ++i){
        if (i < rest){
            last++;
        } 
            
        last_v.push_back(last);
        first_v.push_back(first);
            
        first = last+1;
        last = first+inc-1;
    }  
    
    
    for( int i = first_v.at(mpi_rank); i <= last_v.at(mpi_rank); ++i )
    {
        std::cerr << "rank: " << mpi_rank << "/" << mpi_size-1 << " @" << procName << " solving the " << i << "-th session" << std::endl; 
    
    
    #else
    // .........................................................................
    // loop over databases
    for( int i=0; i < database.size(); ++i )
    {
    #endif
        string dbname  = database[ i ];
        string version = versions [ i ];
        
        log<NOTHING>("********* processing session ") % dbname % " " % std::string('*',9);
       
        // write dbname in logfile for each session
        outstream << setw(9) << setfill(' ') << left << dbname << "|";
      
        // get and write some more information concerning the session (e.g. duration, networkvolume, name)
        ivg::sessinfo info = masterfile.get_session_info(dbname);
        outstream << setw(10) << setfill(' ') << left << info.name << "|";
        outstream << setw(2) << setfill(' ') << right << setprecision(0) << fixed << info.duration << "h|";
        outstream << setw(4) << setfill(' ') << right << setprecision(0) << fixed << info.volume << "|";
      
        // initialize _setup / _name / _date / _ephem / _eops
        string sess_info = "";

        stringstream sstmp;
	int sesnum=-1;
	if (dbArg.isSet()) {
	  // -e option : sessions are parsed from argv
	  for( int j=0; j<setup[ "sessions" ].getLength(); ++j )
	    {
	      
	      string tmpstr=setup[ "sessions" ][ j ][ "dbname" ];
	      
	      if (dbname==tmpstr){
		sesnum=j;
		break;
	      }
	    }
		
	} else {
	  sesnum=i;  
        
	}        
	
	ivg::Session S( &setup, dbname, &ephem, sesnum );
        
        try
        {
            // construct wrapper
	    
	    int year = masterfile.get_session_info(dbname).date.get_int_year();
            if(year == 0){
	      if (dbname.size()<=9){
                year = stoi(dbname.substr(0,2));
                if (year<79)
                    year += 2000;
                else
                    year += 1900;
	      } else {
		year=stoi(dbname.substr(0,4));
	      }
	      
            }
	    
            string vgos_dir = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[2];
            vgos_dir += "/" + std::to_string(year) + "/" + dbname + "/";
            
            ivg::Wrapper wrapper;
            if ( (bool) setup ["use_wrapper"]){
                wrapper = ivg::Wrapper(vgos_dir,dbname,(const char *)setup["vgosdb_editing"]);
                sessionizer.setWrapper_ptr(&wrapper);
            }
            // load information from each database into session
            // setup, S, "02OCT16XA", "5" or "opa2014a"
            // sets following member variables of the session object S
            // in case of session_type = ngs
            // _start, _end, _nobs, _trf, _crf, _scans, _param_list
            // in case of session_type = snx
            // _start, _end, _trf, _crf, _param_list
            // in case of session_type = skd
            // _start, _end, _trf, _crf
	    
            std::string dir = sessionizer.load( &S, &setup, dbname, version);
	    
	    bool only_create_sky_plot = session_type == "skd" && ( !setup.exists( "SKED" ) || !(bool)setup["SKED"][ "apply" ]) && ( !setup.exists( "SIM" ) || !(bool)setup[ "SIM" ]["apply"]);
	   
	    
	    if( !only_create_sky_plot )
	    {
	     
            ivg::Schedule schedulator( &S );
            
            
            if( setup.exists( "SKED" ) && (bool)(setup)[ "SKED" ]["use_year_ses_structure"]
                && ((bool)setup[ "SKED" ]["apply"] || (setup.exists( "SIM" ) && (bool)setup[ "SIM" ]["apply"]) )){
                    std::string skedout = outdir+"/"+S.getStart().get_date_time("YYYY")+"/";

                    if( !directory_exists(skedout) ){
                        mkdir(skedout.c_str() , 0700);
                    }
                    skedout += S.get_name() + "/";
                    if( !directory_exists(skedout) ){
                        mkdir(skedout.c_str() , 0700);
                    }
                    setup[ "outdir" ].operator =( skedout );
            }
            
            // in case of scheduling, the simulation is also performed
	    
	    if( setup.exists( "SKED" ) && (bool)setup[ "SKED" ]["apply"] )
            {        
                
//                std::vector< std::vector<std::string> > networks = { {"WETTZ13N","KOKEE12M","ONSA13NE"},
//                                                                     {"ISHIOKA", "ONSA13SW","WETTZ13S"},
//                                                                     {"WETTZ13N","ISHIOKA","ONSA13NE"},
//                                                                     {"KOKEE12M", "ONSA13SW","WETTZ13S"}
//                                                                   };
//                
//                std::map<std::string, ivg::Matrix> emap = schedulator.minElevation( S.getStart(), networks );
//
//                for(auto& a : emap){
//                    if( a.second.numel() > 0 )
//                        a.second.save_bin("/home/corbin/Documents/visibleSkyArea/" + a.first +"_mine.bin");
//                }
//                
//                exit(0);
                                
                schedulator.start_scheduling();
                
                if( (int)(setup)["SKED"]["approach"] >= 4){
                    
                    lps::GUI gui(&S);
                    gui.run(argc, argv);
                }

                string skd_file = setup[ "outdir" ];
                sessionizer.write_skd( &S, skd_file + S.get_name() + ".skd");
                
                if( (setup.exists( "SIM" ) && !(bool)setup[ "SIM" ]["apply"]) ){
                    continue;
                }
            }
            
            
            if( (setup.exists( "SIM" ) && (bool)setup[ "SIM" ]["apply"] && (bool)setup["SIM"]["load"][0] == false ))
            {
               
                ivg::Simulation simulator( &S );
                // compute
                // delay
                // o-c
                // sigma
                simulator.simulate();

                if( setup[ "SIM" ].exists("eop") ){
                    // in case of SIM the session constructor uses different settings for eops.
                    // load new eop time series  for analysis 
                    // update aprioris
                    // update crf2trf Matrix in all scans
                    S.change_eop_series( setup["eop"], S.getStart() );
                }
            }
                       
            // calculate theoretical delay, and setup least squares solution
	    
	    S.init_vgosdb_ngs_solution();
	    
            // write stations of session in logfile (lettercode)
            vector<string> station_names = S.get_trf_ptr()->get_station_names(ivg::staname::lettercode);
            outstream << setw(2) << setfill(' ') << right << station_names.size() << "|";
            copy ( station_names.begin () , station_names.end () , ostream_iterator <string>( sstmp , "" ) );
            outstream << setw(42) << setfill(' ') << left << sstmp.str() << "|";
	   
            // only proceed if we are not in VASCC2015 mode
            if(!(bool)setup.exists("vascc2015") || ((bool)setup.exists("vascc2015") && (bool)setup["vascc2015"] == false))
            {
                // modify the parameterization
	       
	        S.modify_parameterization();
		 
                // reduce parameter and constrain them
                S.reduce_and_constrain();
		
                // export datum free normal equations in SINEX format; parameters
                // others than stations or sources might be constrainted
                std::string snx_type = (const char*)setup["export_snx"]["type"];
		//   if((bool)setup["export_snx"]["save"] && snx_type == "NEQ")
                //{
                //    std::string snx_dir = setup["export_snx"]["dir"];
                //    if( snx_dir == "" )
                //        snx_dir = dir;
                //    std::string name = setup["export_snx"]["version"];
		    
                //    if( name == "" )
                //        sessionizer.write_snx(&S, snx_dir+"/"+dbname+".snx" );
                //    else
                //        sessionizer.write_snx(&S, snx_dir+"/"+dbname+"_"+name+".snx" );
                //}
		
                // perform least squares solution
                // contains create_nnr_nnt_equations
                sess_info = S.solve();
		
                // if solving is successfull, we write the session information to the logfile
                outstream << sess_info;
		
                // check if station estimates are too big, create output for logfile
                vector<string> bad_stations = S.check_stations_estimates();
                if(bad_stations.size() > 0)
                {
                    stringstream ss;
                    ss << " WARNING: STA_CHECK " << bad_stations.size() << " stations exceeding threshold: ";
                    std::copy(bad_stations.begin(), bad_stations.end(),std::ostream_iterator<std::string>(ss," "));
                    ss << "[UP,EAST,NORTH](m)";
                    // create output for logfile and display
                    cerr << "!!!" << ss.str() << endl;
                    outstream << ss.str();
                }
		
                // path to save outputs (snx&residuals)
                string dir = setup[ "outdir" ];
		
                // writing residuals for extern analysis/plotting         
                if( setup.exists("export_resid") && (bool)setup["export_resid"]["save"] )
                {
                    std::string resid_dir = setup["export_resid"]["dir"];
                    if( resid_dir == "" )
                        resid_dir = dir;

                    std::string resid_format = setup["export_resid"]["format"];
                    sessionizer.write_residuals(&S, resid_dir+"/"+dbname+"_resid.txt",resid_format );
                }                
		
                // writing outlier information files if selected in configfile
                if( (bool)setup["outliers"]["save"][0] ) // && !( (bool)setup["outliers"]["load"] ) )
                    sessionizer.write_outliers(&S, "/data/bakkari_outliers/"+dbname+"_out.txt");

                if( (bool)setup["SIM"]["save"][0] )
                    sessionizer.write_groupdelay( &S );
		 
                // write results in ivg specific format
                sessionizer.write_results(&S, dir+"/"+dbname+"_result.txt", true );
		

		// write each session (solution) results in SNX
		if((bool)setup["export_snx"]["save"] && snx_type == "NEQ")
                {
                    std::string snx_dir = setup["export_snx"]["dir"];
                    if( snx_dir == "" )
                        snx_dir = dir;
                    std::string name = setup["export_snx"]["version"];
		    
                    if( name == "" )
		      sessionizer.write_snx(&S, snx_dir+"/"+dbname+".snx" ,true);
                    else
		      sessionizer.write_snx(&S, snx_dir+"/"+dbname+"_"+name+".snx", true );
                }
		
		if((bool)setup["export_snx"]["save"] && snx_type == "COVA")
                {
                    std::string snx_dir = setup["export_snx"]["dir"];
                    if( snx_dir == "" )
                        snx_dir = dir;

                    std::string name = setup["export_snx"]["version"];
		    
                    if( name == "" )
                        sessionizer.write_snx(&S, snx_dir+"/"+dbname+".snx" );
                    else
                        sessionizer.write_snx(&S, snx_dir+"/"+dbname+"_"+name+".snx" );
                }
		
                // write troposphere results for each session (solution) SNX TRO file
                if((bool)setup["export_snx_tro"]["save"])
                {
                    std::string tro_dir = setup["export_snx"]["dir"];
                    if( tro_dir == "" )
                        tro_dir = dir;

                    std::string name = setup["export_snx_tro"]["version"];

                    if( name == "" )
                        sessionizer.write_snx_tro( &S, tro_dir+"/"+dbname+".tro" );
                    else
                        sessionizer.write_snx_tro( &S, tro_dir+"/"+dbname+"_"+name+".tro" );
                }                
		
                // show residuals analysis application in case of "-r"-flag
                if( resid )
                {                
                    // creating residuals plus corresponding information (e.g. azi/ele of stations, stds, etc.)
                    
		    S.create_solution_info();
                   
                    // generate mainwindow for residual plotting
                    QApplication a(argc, argv);
		   
                    Statistics stats_figure(&S);
		    
                    a.exec();
                }
		if( setup.exists("export_eop") && (bool)setup["export_eop"]["save"] )
		  {
		    std::string name = setup["export_eop"]["filename"];
		    sessionizer.write_eop_file(&S, dir+"/"+name,info.code);
		  }
		}
	    }
                // create skyplots and compute coverage vvvvvvvvvvvvvvvvvvvvvvvvvvvv
	    if( only_create_sky_plot || (
			(   setup.exists( "SKED" ) 
                        && (    ( (int)setup[ "SKED" ]["approach"] < 4 && (bool)setup[ "SKED" ]["apply"] ) 
                                || !(bool)setup[ "SKED" ]["apply"]  )
                    )
                    && 
                    ( setup.exists( "SIM" ) && (bool)setup[ "SIM" ]["apply"] ) && (bool)setup[ "SKED" ]["createPlots"] )){
                    
                    S.calc_transits(  0  );
                    S.create_crf_trf_inidices();
                    S.get_complete_transits().print_transits();
                    S.get_common_transits().print_transits();
                    //S.get_complete_transits().getBaselineWiseTransits().save_transits((const char*)setup[ "outdir" ] );
                    
                    lps::GUI gui(&S);
                    QApplication a(argc, argv);
                    gui.initGUI();
                    gui.showSession();

                    // save pdf
                    for(size_t sta_idx = 0; sta_idx <  S.get_trf_ptr()->get_number_stations(); ++sta_idx){

                       StationDialog* sd = gui.stationViews[sta_idx];
//                           sd->show_tree();

                       std::stringstream ss;
                       ss << (const char*)setup[ "outdir" ] << "/"
                          << S.get_name() << "_" << S.get_trf_ptr()->get_station(sta_idx)->get_name(ivg::staname::ivs_name) << "_sim.pdf"; 
                       std::string path = ss.str();
                       sd->print_pdf( path );

                   }

                    if ( (bool)setup[ "SKED" ]["auto_close_windows"]){
                        a.closeAllWindows();
                        a.exit(0);
                    } else {
                        a.exec();
                    }

                }
                // create skyplots and compute coverage ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            
            
        }
        catch( libconfig::SettingException &err )
        {
            cerr << "libconfig::SettingException: " << err.what() << " at " << err.getPath() << endl;
            exit( -1 );
        }
        catch(std::exception& e)
        {
            // create bad session file containing the dbname, the affected stations and the occured error
            cerr << "std::exception: " << e.what() << " in " << dbname << endl;
            
            // if error occured before init_vgosdb_ngs_solution() finished, substitute stationlist
            if(sstmp.str().empty())
                outstream << "  |                                          |";
            
            // if error occured before solve() finished, substitute
            if(sess_info.empty())
                outstream << "VFAC: 0.000e+00 WRMS: 0.000e+00 OUT:  0.00%";
            
            outstream << " ERROR: " << e.what();

            // show parameter in case of LAPACK-error
            std::string tmp = e.what();
            if(tmp.find("LAPACK")!=std::string::npos)
            {
                std::stringstream tokenizer(tmp);
                std::string a, b, c;
                int pos;
                tokenizer >> a >> b >> c >> pos;
                ivg::Param *p = S.get_param_list_ptr()->get_param(pos-1);
                p->show();
                outstream << " at " << p->get_name() << " " << p->get_typename() << " order " << p->get_order() << " " << p->get_epoch().get_double_mjd() << "mjd";
            }
        }
        
        outstream << endl;
        
         #ifdef USE_MPI
        // receive outstreams from other processes 
        if(mpi_rank == 0){
            
            master_stream << outstream.str();
            
            MPI_Status st;
            for(unsigned int rank = 1; rank < mpi_size; ++rank){
                
                if( i < inc  || rank  < rest)
                {
                    int numstr;
                    MPI_Recv(&numstr, 1, MPI_INT, rank, 12, MPI_COMM_WORLD,&st);
                    char cstr[numstr];
                    MPI_Recv(cstr, numstr, MPI_CHAR, rank, 12, MPI_COMM_WORLD,&st);    

                    master_stream << cstr;
                
                }
            }

        // send outstream to master
        } else {
            std::string outstring = outstream.str();
            const char* cstr = outstring.c_str();
            int numstr =  outstring.length();

            MPI_Send (&numstr, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
            MPI_Send (cstr, numstr, MPI_CHAR, 0, 12, MPI_COMM_WORLD);
        
        }
        
        outstream.str(std::string());
        #endif
    }
  

    
    #ifdef USE_MPI
    if(mpi_rank == 0){
        master_stream << endl;
        master_stream.close();
        std::cout << " =================================================================== " << std::endl;
        std::cout << "Time required: " << fixed << timeCount.toc() << " seconds." << std::endl;
    }
    MPI_Finalize(); // end MPI
    #else
    std::cout << " =================================================================== " << std::endl;
    std::cout << "Time required: " << fixed << timeCount.toc() << " seconds." << std::endl;
    outstream << endl;
    outstream.close();
    #endif
    
    return ( 0 );
}




	
