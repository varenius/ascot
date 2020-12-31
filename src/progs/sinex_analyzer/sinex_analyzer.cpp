#include <cstdlib>
#include <QApplication>
#include <libconfig.h++>
#include <tclap/CmdLine.h>
#include "analyzer.h"
#include "session_inout.h"
#include "transformation.h"
#include "masterfile.h"
#include "station_analysis.h"

using namespace libconfig;
using namespace ivg;
using namespace ivgat;
using namespace std;

loglevel g_verbose;

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
    
    // parse commandline arguments
    cmd.parse( argc, argv );
    string controlfile = m1Arg.getValue();
    g_verbose = (loglevel) m2Arg.getValue();

    // loading controlfile
    Config cfg;
    cfg.readFile( controlfile.c_str() );
    Setting& setup= cfg.lookup( "setup" );
    
    // nonsense initializations because otherwise compilation error
    ivg::Masterfile supernonsense;
    ivg::Session_inout nonsense("snx",supernonsense);
    ivgat::Station_analysis meganonsense;
    
    // start analyzer GUI
    QApplication a(argc, argv);
    Analyzer main(&setup);

    return a.exec();
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
