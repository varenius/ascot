/***************************************************************************** 
 * analysis tools for clock parameters                                       *
 *                                                                           *
 *                                                                           *
 * 2016-05-23 - SH                                                           *
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
#include "lapack_wrapper.h"
#include "tictoc.h"
#include "session.h"
#include "session_inout.h"
#include "date.h"
#include "trf.h"
#include "station_analysis.h"
#include "lsa.h"
#include "fit.h"
#include <tclap/CmdLine.h>
#include "logger.h"
#include <cstdlib>
#include <libconfig.h++>
#include <QApplication>
#include "plot.h"
#include "sinex.h"
#include "auxfunc.h"

using namespace libconfig;

loglevel g_verbose;
void * g_ephis;

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
    
   // .........................................................................
   // initialize timer    
   tictoc timeCount;
   timeCount.tic();
    
   std::cout << " \n =================================================================== " 
             << std::endl;
   std::cout << "     > > > > > ivg::ASCOT (analysis tools) < < < < <     " 
             << std::endl;
   std::cout << "         A N A L Y S I S   C L O C K S " 
             << std::endl << std::endl;
   
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
   
   std::vector<std::string> stations;
   for( int i=0; i<(setup)[ "clock_time_series" ]["station"].getLength(); ++i )
   {
      std::string sta = (const char *)setup[ "clock_time_series" ]["station"][i]; 
      stations.push_back( sta );
   }
   
   try
   {       
    // create map containing different solutions with all databases
    std::map< std::string, std::vector<std::string> > sols;    
    for( int i=0; i<setup[ "sessions" ].getLength(); ++i )
    {  
        std::string versions = setup[ "sessions" ][ i ][ "version" ];
        std::string version;
        istringstream iss(versions);
        std::string dbname  = setup[ "sessions" ][ i ][ "dbname" ];

        while (std::getline(iss, version, ','))                    
            sols[ version ].push_back( dbname );
    }
     
    const std::vector<QCPScatterStyle::ScatterShape> symbols = {
            QCPScatterStyle::ssTriangle, QCPScatterStyle::ssCircle,
            QCPScatterStyle::ssSquare, QCPScatterStyle::ssStar,        
            QCPScatterStyle::ssDiamond, QCPScatterStyle::ssCross, 
            QCPScatterStyle::ssTriangle, QCPScatterStyle::ssCircle,
            QCPScatterStyle::ssSquare, QCPScatterStyle::ssStar,        
            QCPScatterStyle::ssDiamond, QCPScatterStyle::ssCross
    };
    
    const std::vector<string> color_values = {
            "#000000", "#00FF00", "#ab5e1a", "#000000", "#8a38be", "#FF34FF", "#ffa200", "#FF0000",
            "#A30059", "#FFDBE5", "#0000A6", "#997D87", "#B79762", "#004D43", "#8FB0FF", "#63FFAC",
            "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
            "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
            "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
            "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
            "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",      
    };    
    
    typedef std::map<std::string, std::vector<std::string> >::iterator mapIter;
    std::string version, dbname;
    QApplication fig( argc, argv );   
    int p = 0;   
    map< std::string, ivg::Matrix > clocks;
    std::string type; 
    int max_order;
    
    // .........................................................................
    // loop over solutions
    for (mapIter iter = sols.begin(); iter != sols.end(); iter++)
    {
        log<RESULT>("*** Solution: '") % iter->first % "' including " % iter->second.size() % " dbs.";
        
        int n = iter->second.size();
        std::string data_path;
        
        for( int i=0; i < n; ++i )
        {    
            version = iter->first;
            dbname = iter->second.at(i); 

            // load information from each database into session
            data_path = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[2];
            std::string path = data_path+version+"/"+dbname+".snx";
            
            // get statistics (variance factor, WRMS, RMS)
            ivg::Sinex snx( path, false );
            if( file_exists( path ) )
            {
                if( (bool)setup[ "clock_time_series" ]["apply"] )
                    clocks = snx.get_clocks( type, max_order );
            }
               
            for( int k=0; k<stations.size(); ++k )
            {      
                Plot plot0; 
                ivg::Matrix clo;
                if( type == "cpwlf" )
                {
                    clo.resize(clocks[stations.at(k)].rows()-2,2,0.0);
                    ivg::Matrix dt = ( clocks[stations.at(k)](":",0) - clocks[stations.at(k)](0,0) ); 

                    for( int i=2; i< clocks[stations.at(k)].rows(); ++i )
                    {
                        clo(i-2,0) = clocks[stations.at(k)](i,0);
                        clo(i-2,1) = clocks[stations.at(k)](i,1) 
                                   + clocks[stations.at(k)](0,1)*dt(i,0) 
                                   + clocks[stations.at(k)](1,1)* pow( dt(i,0), 2.0 );
                    }
                    log<RESULT>("*** clock parameters for station: ") % stations.at(k) 
                              % "; parametrization type: " % type % ".";                
                }
                else
                {             
                    double start = snx.get_start_epoch().get_double_mjd();
                    double end = snx.get_end_epoch().get_double_mjd();
                    double inc = (end-start) / 10.0;
                    ivg::Matrix T(start, inc, end, 1);
                    ivg::Matrix dT = T(":",0) - T(0,0);
                    clo.resize(dT.rows(),2,0.0);

                    for( int i=0; i<dT.rows(); ++i )
                    {
                        clo(i,0) = dT(i,0);     
                        for( int j=0; j<=max_order; ++j )
                            clo(i,1) += clocks[stations.at(k)](j,1)* pow( dT(i,0), j );
                    }            
                    log<RESULT>("*** clock parameters for station: ") % stations.at(k) 
                              % "; parametrization type: " % type 
                              % " of order " % max_order % ".";                
                }
                clo.set_sub(0,1,clo(":",1)*1e9);
                
                ivg::Lsa nonsense; // for compilation only
                ivgat::Fit fitter;
                fitter.polyfit( clo(":",0), clo(":",1), 1 );
                ivg::Matrix x = fitter.get_param();               
                                                
                log<DETAIL>("*** max.: ") % clo(":",1).max() 
                          % ", min.: " % clo(":",1).min() 
                          % ", diff.: " % (clo(":",1).max()-clo(":",1).min()) % ".";

                ivg::Date epoch( 2014, 7, 8, 19, 0, 0.0 );
                double y = x(1)* epoch.get_double_mjd() + x(0);
                double y2 = x(1) / 24.0;
                log<DETAIL>("*** coefficients (linear fit) at ") 
                         % epoch.get_date_time("YYYY-MO-DD,HH:MI:SS") 
                         % ": (1) " % y % " ns, (2) " % y2 % " ns/h";
                
                // PLOT data
                if( (bool)setup[ "clock_time_series" ]["apply"] )
                {      
                    plot0.plot_mjd_series( clo(":",0),clo(":",1),
                                     { QColor(color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                                       QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, "hh:mm:ss");
                    plot0.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at( p ).c_str() )));
                    plot0.get_plot()->graph()->setErrorBarSize(20.0);            
                    plot0.set_title( ("clocks for "+stations.at(k)).c_str() );
                    QCustomPlot *plt_ptr0 = plot0.get_plot();
                    plt_ptr0->xAxis->setLabel("");
                    plt_ptr0->yAxis->setLabel( "[ns]" );
                    QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
                    plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font0.family(), 20) );
                    plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font0.family(), 20) ); 
                    plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font0.family(), 20) );
                    plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font0.family(), 20) );
                    plt_ptr0->replot();     
                    p++;
                }
            fig.exec();
            }
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
   

   
   std::cerr << " =================================================================== \n" 
             << "Time required: " << fixed << timeCount.toc() << " seconds." << std::endl;
	    
// ============== ENDE 

   return ( 0 );
}
	
	
	
	
	
