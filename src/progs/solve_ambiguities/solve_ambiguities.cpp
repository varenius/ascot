
/* 
 * File:   solve_ambiguities.cpp
 * Author: corbin
 *
 * Created on 3. Juni 2016, 15:45
 */


#include "logger.h"
#include "ivg_const.h"
#include "db_download.h"
#include "session.h"
#include "session_inout.h"
#include "tictoc.h"
#include "matrix.h"
#include "lsa.h"
#include "vgosdb.h"
#include "date.h"
#include "statistics.h"
#include "network.h"


#include <boost/algorithm/string/predicate.hpp>
#include <libconfig.h++>
#include <QApplication>
#include <string>
#include <tclap/CmdLine.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>



loglevel g_verbose;

int main(int argc, char** argv) {

    /*
     * 
     * 14APR06XK TsWz nice ambigs visible in S and X band
     * 15JAN05XK NyTsWz  ambigs in Sband for baseline NyWz + ambigs between baselines
     * 15JAN05XU KkWz Xband has some ambigs
     * 15JAN06XU KkWz no ambigs
     * 15JAN12XK NyTsWz good session to show why using the smallest resids is not always a good idea!!!
     * 15JAN19XK NyShTsWz -> 3 triangles with ambiguities
     * 15JUN29XK S band has different ambiguity spacings
     */
    
    // get the directory this program is called from. Needed for the default value of the configfile
    std::string call_dir = "";
    std::string::size_type lastSlash = std::string(argv[0]).find_last_of("/");
    if (lastSlash !=  string::npos){
        call_dir = std::string(argv[0]).substr(0,lastSlash); 
    }
    

    try {
        // handle commandline arguments-----------------------------------------
        // ---------------------------------------------------------------------
        
        //Construktor needs: message, delimeter, version
        TCLAP::CmdLine cmd("auto_intensive", ' ', "alpha 0.1");

        //flag 	- The one character flag that identifies this argument on the command line.
        //name 	- A one word name for the argument. Can be used as a long flag on the command line.
        //desc 	- A description of what the argument is for or does.
        //req 	- Whether the argument is required on the command line.
        //value 	- The default value assigned to this argument if it is not present on the command line.
        //typeDesc 	- A short, human readable description of the type that this object expects. This is used in the generation of the USAGE statement. The goal is to be helpful to the end user of the program.
        //v 	- An optional visitor. You probably should not use this unless you have a very good reason. 

        // verbose level
        TCLAP::ValueArg<int> verboseArg("v", "verbose", "Verbose Level (0=NOTHING, 1=INFO, 2=DETAIL, 3=RESULT, 4=WARNING, 5=ALL)", false, 4, "int");
        cmd.add(verboseArg);

        // dbname
        TCLAP::ValueArg<std::string> dbArg("e", "explicitDB", "e.g. 93NOV05XU", false, "", "string");
        cmd.add(dbArg);

        // controlfile
        TCLAP::ValueArg<std::string> cntArg("c", "controlfile", "/home/ascot/cnt", false, call_dir + std::string("/../src/progs/solve_ambiguities/solve_ambiguities.cfg"), "string");
        cmd.add(cntArg);

        // arcfile
        TCLAP::SwitchArg write_arc_switch("a", "arcfile", "create an arc file with all sessions with successfull ionosphere corretion", false);
        cmd.add(write_arc_switch);

        TCLAP::SwitchArg load_switch("d", "download", "download the databases", false);
        cmd.add(load_switch);

        TCLAP::SwitchArg force_switch("f", "force", "force download the databases even if it is already on local directory", false);
        cmd.add(force_switch);

        TCLAP::SwitchArg plot_switch("p", "plot", "open a window to show the residuals", false);
        cmd.add(plot_switch);
        
        TCLAP::SwitchArg use_arc_switch("u", "usearc", "uses the databases spefified in the controlfile", false);
        cmd.add(use_arc_switch);

        // Parse the argv array.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg. 
        g_verbose = (loglevel) verboseArg.getValue();

        std::string controlfile = cntArg.getValue();

        bool load = load_switch.getValue();

        bool force = force_switch.getValue();

        bool plot = plot_switch.getValue();

        bool write_arc = write_arc_switch.getValue();
        
        bool use_arc = use_arc_switch.getValue();

        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------

        std::cout << " \n =================================================================== " << std::endl;
        std::cout << "     > > > > > ivg::ASCOT (auto_intensive) < < < < <     " << std::endl << std::endl;

        // initialize timer
        tictoc timeCount;
        timeCount.tic();

        // loading control file (only clocks estimated, group and single band)
        libconfig::Config cfg;
        try {
            cfg.readFile(controlfile.c_str());
        } catch (libconfig::ParseException & err) {
            std::cerr << "libconfig::" << err.what() << " in " << err.getFile() << " at line " << err.getLine() << std::endl;
            exit(-1);
        }
        libconfig::Setting& setup = cfg.lookup("setup");


        try {
            //TODO check wether local paths are existent

            // get paths from controlfile ---------------------------------------------------------
            std::string local_log_path = (const char*) setup["ltndir"]["log"];
            std::string local_vgos_path = (const char *) get_list_element(setup["datadirs"], setup["session_type"])[2];
            std::string remote_log_path = (const char*) setup["remotedir"]["log"];
            std::string remote_vgos_path = (const char*) setup["remotedir"]["vgos"];

            std::string ltn_data = (const char*) setup["ltndir"]["data"];

            std::string masterfile_path = (const char*) setup["masterfiles"];

            int start_year = setup["period"]["start"];
            int end_year = setup["period"]["end"];

	    bool PhaseSolution;
	    if (setup.exists("phase_solution"))
	      PhaseSolution= setup["phase_solution"];
	    else
	      PhaseSolution=false;
	    
            // initzialise masterfiles ------------------------------------------------------------
            ivg::Masterfile master(masterfile_path, ivg::mastertype::both, start_year, end_year);

            // fill database vector ---------------------------------------------------------------
            
            // For all sessions in the database vector the ambiguities are resolved
            std::vector<std::string> database;
            
            if (dbArg.isSet()) {
            // -e option : sessions are parsed from argv
                database.push_back(dbArg.getValue());
            } else if( use_arc ){
            // -u option : sessions specfified in contolfile are used
                database.reserve( setup[ "sessions" ].getLength() );
                for( int i=0; i<setup[ "sessions" ].getLength(); ++i )
                {
                    database.push_back( setup[ "sessions" ][ i ][ "dbname" ] );
                }          
            } else {
            // no optioon : all sessions in masterfile between start and end year are used
                for (ivg::sessinfo &info : *master.get_sessions()) {
                    database.push_back(info.dbname);
                }

            }

            // Download databases -----------------------------------------------------------------
            // if -d option is set missing databases are loaded from the remote server
            // if -f option is set  already existing databases a loaded too
            if (load) {
                // initialize timer for download time
                tictoc timeCount2;

                ivg::Db_download d(local_log_path, local_vgos_path, remote_log_path, remote_vgos_path, &master, ltn_data);

                unsigned n = database.size();
                for (unsigned i = 0; i < n; ++i) {
                    timeCount2.tic();
                    std::cout << "> > > > > > > > > > " << database[i] << "(" << i << "/" << n << ") < < < < < < < < < <" << std::endl;
                    d.download_vgosDB(database[i], force);
                    std::cout << "> > > > > > > > > > " << timeCount2.toc() << " secondes < < < < < < < < < <" << std::endl;
                }
            }

            // process databases ------------------------------------------------------------------


            // the rescue files contains the database names that are allready processed. If a unforseen errror occures 
            // not all sessions have to be calculated again.
            std::string rescue_file = controlfile + ".rescue";
            
            // if -a option is passed all sucessfully processed session ( ionosheric correction could be calculated )
            // are saved in an arcfile
            std::string arc_file = (const char*) setup["logdir"];
            arc_file += "/auto_int_iono_ok.arc";

            std::ofstream rf; //rescue file, write only
            std::fstream af; // arcfile, read and write

            const char *rescueFile = rescue_file.c_str();
            const char *arcFile = arc_file.c_str();

            bool rescue = file_exists(rescue_file);
            rf.open(rescueFile, ios::app);

            if (!rf.good()) {
                std::cerr << "can not write rescue file" << std::endl;
            } else {
                // if a rescue file is found all databses that are in the rescue file are removed from the database vector and therfore not processed again.
                if (rescue) {
                    std::cout << "A rescue file has been found. Do you want to use the file <y> or to start a new solution <n>?" << std::endl;
                    std::cout << "You can skip the last processed session and proceed (use rescue file) with <s> ...";
                    std::string input;
                    std::cin >> input;
                    if (input.find("y") != std::string::npos || input.find("s") != std::string::npos) {
                        std::ifstream g;
                        g.open(rescue_file);
                        if (!g.good()) {
                            std::cerr << "resucue file could not be opened!!!" << std::endl;
                        } else {
                            // parse the rescue file
                            std::string line;
                            while (std::getline(g, line)) {
                                std::istringstream iss(line);
                                std::string db;
                                iss >> db;
                                vector<string>::iterator it = std::find(database.begin(), database.end(), db);
                                if (it != database.end()) {
                                    database.erase(it);
                                }
                            }
                            // if skip is selected the first entry in the databse vector is removed and written into the rescue file
                            if (input.find("s") != std::string::npos) {
                                std::cout << "skipping session " << database.front() << std::endl;
                                rf << database.front() << std::endl;
                                database.erase(database.begin());
                            }
                        }
                        g.close();
                    }
                }

                // each entry in the arcfile ecept the last have to end with an comma
                // if the last entry in an existing file has not a comma the bool is set true;
                bool add_comma = false;

                // if -a option is used the arcfile outstream will be opened
                if (write_arc) {
                    //if arcfile is already existing get the last char 
                    if (file_exists(arc_file)) {
                        af.open(arcFile, ios::in); // open af for reading
                        if (af.good()) {
                            //get last line in arcfile
                            std::string lastLine, line;
                            while (getline(af, line)) {
                                lastLine = line;
                            }

                            // if last line ends not with a comma and is not a comment add a comma
                            if (!boost::algorithm::ends_with(lastLine, ",") && !boost::algorithm::starts_with(lastLine, "#"))
                                add_comma = true;

                        } else {
                            std::cerr << "can not read arc file" << std::endl;
                        }
                        af.close();
                    }

                    // open af for writing
                    af.open(arcFile, ios::app);
                    std::cout << "writing arc file to" << arcFile << std::endl;
                }
                if (af.good() == false && write_arc == true) {
                    std::cerr << "can not write arc file" << std::endl;
                } else {
                    // write date to arc file ----
                    //Creation time
                    ivg::Date d;
                    d.now();
                    if (add_comma)
                        af << "," << std::endl;

                    af << "# run: " << d.get_date_time("YYYY/MO/DD HH:MI:SS");

                    // load ephemeris only once to avoid multiple jpl_init_ephemeris
                    char nams[400][6];
                    double vals[400];
                    string ephfile_name = setup["definitions"]["ephemerides"][(const char*) setup["ephemerides"]];
                    void *ephem = jpl_init_ephemeris(ephfile_name.c_str(), nams, vals);

                    // initialization of session_inout with defined session_type
                    string session_type = (const char *) get_list_element(setup["datadirs"], setup["session_type"])[1];

                    ivg::Session_inout sessionizer = ivg::Session_inout(session_type);

                    // write 'bad-session' file
                    // name of the file based on the time of execution
                    std::string dir = setup[ "logdir" ];
                    ivg::Date now;
                    now.now();
                    std::string ascot_log = dir + "/" + "ascot_" + now.get_date_time("DOY_HH_MI_SS") + ".log";
                    ofstream outstream(ascot_log.c_str());
                    
                    int algorithem = (int)setup["resolve_type"];
                    
                    std::string new_editing = (const char*) setup["new_editing"];

                    // for loop over dbs where they are finally processed =============================================
                    // ================================================================================================
                    for (string &db : database) {
                        
                        // start session  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                        if (g_verbose >= 2) {
                            std::cout << "______________________________________ _________ ______________________________________" << std::endl;
                            std::cout << "______________________________________ " << db << " ______________________________________" << std::endl;
                        }
                        
                        // construct wrapper
                        int year = master.get_session_info(db).date.get_int_year();
                        if(year == 0){
                            
			  if (db.size()<=9){
			    year = stoi(db.substr(0,2));
			    if (year<79)
			      year += 2000;
			    else
			      year += 1900;
			  } else {
			    year=stoi(db.substr(0,4));
			  }
                        }
                        string vgos_dir = (const char *) get_list_element(setup["datadirs"], setup["session_type"])[2];
                        vgos_dir += "/" + std::to_string(year) + "/" + db + "/";
                        ivg::Wrapper wrapper(vgos_dir, db, (const char *) setup["vgosdb_editing"]);
                        ivg::Vgosdb vgos(vgos_dir, &wrapper);                    
                        
                        // calculate residuals for X and S band -------------------------------------------------             
                        std::map<ivg::band, Network> networks;
                        std::array<bool, 2> band_ok({true, true});
                        
                        for (int band = ivg::band::X; band < ivg::band::MAXBANDTYPE; band++) {

                            std::string version = "-1";

                            // write dbname in logfile for each session
                            outstream << setw(9) << setfill(' ') << left << db << ": ";


                            // initialize _setup / _name / _date / _ephem / _eops
                            string sess_info = "";
                            stringstream sstmp;

                            ivg::Session S (&setup, db, &ephem, 0);
                            S.set_band_type((ivg::band)band);
                            S.setAmbigRes(true);
                            
                            try {

                                sessionizer.setWrapper_ptr(&wrapper);

                                // load information from each database into session
                                // _setup, S, "02OCT16XA", "5" or "opa2014a"
                                // sets following member variables of the session object S
                                // in case of session_type = ngs
                                // _start, _end, _nobs, _trf, _crf, _scans, _param_list
                                // in case of session_type = snx
                                // _start, _end, _trf, _crf, _param_list
                                sessionizer.load(&S, &setup, db, version);


                                // write stations of session in logfile (lettercode)
                                vector<string> station_names = S.get_trf_ptr()->get_station_names(ivg::staname::lettercode);
                                copy(station_names.begin(), station_names.end(), ostream_iterator <string>(sstmp, ""));
                                outstream << setw(42) << setfill(' ') << sstmp.str();
                                
                                /*
                                 * find the best reference clock. Because a star-form adjustment is performed
                                 * the reference station needs observations to any other involved station.
                                 */

                                 std::string best_ref = S.get_best_refstation(  setup["RefClockStationList"] );
                                 setup["RefClockStationList"] = best_ref;
				 std::cout << "Ref: " << best_ref << endl;
				 
                                /*
                                 * calculate theoretical delay, and _setup least squares solution
                                 - object of Class ivg::Lsa() is created
                                 - Jacobian Matrix is calculated
                                 - observed minus computed is calculated 
                                 - data is eliminated ( rows in A, and corresponding observations are deleted)
                                 - weight matrix is calculated
                                 */
                                S.init_vgosdb_ngs_solution();

                                /*
                                 * modify the parameterization
                                 - clock breaks are inserted (new cols in A, some entries are set to 0)
                                 - parameters are fixed (cols in A and corresponding entries in paramvector are deleted)
                                 - reduce_flag is set
                                 */
                                S.modify_parameterization();

                                /* 
                                 * reduce parameter and constrain them
                                 - A is divided in two parts: A1 remaines, A2 is reduced
                                 - fully de-correlated jacobian matrix and observations are calculated
                                 - system of normal equations is created
                                 */
                                S.reduce_and_constrain();
                                
                                /* 
                                 * perform least squares solution
                                 * contains create_nnr_nnt_equations
                                 - A might change size if outliers are detected
                                 - post fit residual are calculated
                                 */
                                //S.solve(station_names.size() > 2);
				S.solve();
                                // if solving is successfull, we write the sesseion information to the logfile
                                outstream << " " << sess_info;

                                // check if station estimates are too big, create output for logfile
                                vector<string> bad_stations = S.check_stations_estimates();
                                if (bad_stations.size() > 0) {
                                    stringstream ss;
                                    ss << " WARNING: STA_CHECK " << bad_stations.size() << " stations exceeding threshold: ";
                                    std::copy(bad_stations.begin(), bad_stations.end(), std::ostream_iterator<std::string>(ss, " "));
                                    ss << "[UP,EAST,NORTH](m)";
                                    // create output for logfile and display
                                    cerr << "!!!" << ss.str() << endl;
                                    outstream << ss.str();
                                }

                                S.create_solution_info();

                                if( algorithem == Network::Algorithem::manual){                    
                                    // generate mainwindow for residual plotting
                                    QApplication a(argc, argv);
                                    Statistics stats_figure( &S );
                                    a.exec();                             
                                }  
                                
                                // solve ambiguities  -----------------------------------------------------------------
                                networks[(ivg::band)band] = Network( &S );
                                
                                networks[(ivg::band)band].compute_baselinewise_integer_ambiguities((Network::Algorithem)algorithem, &S);
                                networks[(ivg::band)band].apply_closure_condition();
                                networks[(ivg::band)band].print_Network();
                                networks[(ivg::band)band].fill_all_delay_and_integer_matrix();


                                // show residuals analysis application in case of "-p"-flag
                                if (plot) {
                                    networks[(ivg::band)band].update_resids_in_session(&S);
                                     
                                    // generate mainwindow for residual plotting
                                    QApplication a(argc, argv);
                                    Statistics stats_figure( &S );
                                    a.exec();
                                }
                                

                                std::vector<short>ia = networks[(ivg::band)band].get_all_integerAmbiguities();
                                std::vector<double> gd = networks[(ivg::band)band].get_all_group_delays().get_data_vec();

                                vgos.create_NumGroupAmbig_file(ia, db, (ivg::band)band, new_editing, PhaseSolution);
                                vgos.create_GroupDelayFull_file(gd, db, (ivg::band)band, new_editing, PhaseSolution);
                                

                            } catch (libconfig::SettingException &err) {
                                cerr << "libconfig::SettingException: " << err.what() << " at " << err.getPath() << endl;
                                exit(-1);
                                band_ok[band] = false;
                            } catch (std::exception& e) {
                                // create bad session file containing the dbname, the affected stations and the occured error
                                cerr << "std::exception: " << e.what() << " in " << db << endl;

                                // if error occured before init_vgosdb_ngs_solution() finished, substitute stationlist
                                if (sstmp.str().empty())
                                    outstream << "                                          ";

                                // if error occured before solve() finished, substitute
                                if (sess_info.empty())
                                    outstream << " VFAC: 0.000e+00 WRMS: 0.000e+00 OUT:  0.00%";

                                outstream << " ERROR: " << e.what();

                                // show parameter in case of LAPACK-error
                                std::string tmp = e.what();
                                if (tmp.find("LAPACK") != std::string::npos) {
                                    std::stringstream tokenizer(tmp);
                                    std::string a, b, c;
                                    int pos;
                                    tokenizer >> a >> b >> c >> pos;
                                    ivg::Param *p = S.get_param_list_ptr()->get_param(pos - 1);
                                    p->show();
                                    outstream << " at " << p->get_name() << " " << p->get_typename() << " order " << p->get_order() << " " << p->get_epoch().get_double_mjd() << "mjd";
                                }
                                
                                band_ok[band] = false;
                                
                            }

                            outstream << endl;

                        }

                        

                        // ionospheric correction  -------------------------------------------------
                        if( band_ok[0] &&  band_ok[1] ){
                            ivg::Matrix delta_tau_x, delta_tau_x_sigma;
                            std::vector<short> error_flag;
			    std::cout << "iono corr" <<endl;
                            Network::get_ionospheric_correction(networks[ivg::band::X], networks[ivg::band::S], delta_tau_x, delta_tau_x_sigma, error_flag);
			    std::cout << "iono corr done" <<endl;
                            vgos.create_IonoGroup_file(delta_tau_x, delta_tau_x_sigma, error_flag, db, ivg::band::X, new_editing, PhaseSolution);
                        }
                        
                        wrapper.write_wrapper(new_editing);
                          
                         // write session to arc file if ionsphere is corrected
                        if( band_ok[0] &&  band_ok[1]){
                            if(write_arc){
                                af << "," << std::endl <<  "{ dbname = \"" << db << "\", version = \"-1\" }";
                            }
                        }else{
                             log<WARNING> ("!!! can not correct ionosphere for session ") % db % " X-band: " %  band_ok[0] % " S-band: " %  band_ok[1];
                        }

                        rf << db << std::endl;
                        
                        
                        // session finished ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    }
                }
            }
            if (write_arc) {
                af.close();
            }

            rf.close();
            remove(rescue_file.c_str());

        } catch (libconfig::SettingException &err) {
            cerr << "libconfig::SettingException: " << err.what() << " at " << err.getPath() << endl;
            exit(-1);
        } catch (std::exception& e) {
            // create bad session file containing the dbname, the affected stations and the occured error
            cerr << "std::exception: " << e.what() << endl;
        }

        std::cout << " =================================================================== " << std::endl;
        std::cout << "Time required: " << fixed << timeCount.toc() << " seconds." << std::endl;


    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

    return 0;
}

