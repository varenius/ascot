/***************************************************************************** 
 * independent VLBI solution                                                 *
 *                                                                           *
 *                                                                           *
 *                                                                           *
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
   std::cout << "         A N A L Y S I S   A T M O S P H E R E " 
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
   std::string type = (const char *)setup[ "tropo_time_series" ]["type"];
   std::string station = (const char *)setup[ "tropo_time_series" ]["station"];

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
            "#000000", "#00FF00", "#ab5e1a", "#FF0000", "#8a38be", "#FF34FF", "#ffa200", "#886F4C",
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
    Plot plot0;    
    int p = 0;   
    map< std::string, ivg::Matrix > trop_delay;
    map< std::string, ivg::Matrix > trop_pred;
    ivg::Session_inout Sio;
    ivg::Session S;    
    
    // .........................................................................
    // loop over solutions
    for (mapIter iter = sols.begin(); iter != sols.end(); iter++)
    {
        log<RESULT>("*** Solution: '") % iter->first % "' including " % iter->second.size() % " dbs.";
        
        int n = iter->second.size();
        std::string data_path,session_type;
        
        for( int i=0; i < n; ++i )
        {    
            version = iter->first;
            dbname = iter->second.at(i); 

            // load information from each database into session
            session_type = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[1];
            data_path = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[2];
            
            std::string infile;
            if( session_type == "snx")
                infile = data_path+"/"+dbname+"_ivg2015x.snx";
            else if( session_type == "resultfile")
                infile = data_path+version+"/"+dbname+"_result.txt";
             
            // read data from SINEX file 
            if( file_exists( infile ) )
            {
                if( session_type == "snx")
                {
                    ivg::Sinex snx( infile, false );

                    if( (bool)setup[ "tropo_time_series" ]["apply"] )
                    {
                        if( type == "zwd" )
                            trop_delay = snx.get_tropo_delays( ivg::reftype::delta );
                        else if( type == "zhd" )
                            trop_delay = snx.get_tropo_delays( ivg::reftype::apriori ); 
                        else if( type == "ztd" )
                            trop_delay = snx.get_tropo_delays( ivg::reftype::estimate ); 
                    }
                }
                // read data from internal ivg::ASCOT result file
                else if( session_type == "resultfile")
                {
                    Sio.read_results(&S,infile,true);
                    
                    if( type == "zwd" )
                        trop_delay = S.get_param_list_ptr()->get_station_dependent_param( ivg::paramtype::zwd, ivg::param_def::estimates );
                    if( type == "zhd" )
                        trop_delay = S.get_param_list_ptr()->get_station_dependent_param( ivg::paramtype::zwd, ivg::param_def::aprioris );
                    if( type == "ztd" )
                        trop_delay = S.get_param_list_ptr()->get_station_dependent_param( ivg::paramtype::zwd, ivg::param_def::totals );
                    
                    // check if there are stochastic prediction parameters in the ivg::ASCOT result file   
                    if( S.get_param_list_ptr()->exist_stoch_param( ivg::paramtype::zwd, station ) )
                    {
                        if( type == "zwd" )
                            trop_pred = S.get_param_list_ptr()->get_station_dependent_param( ivg::paramtype::zwd, ivg::param_def::estimates, true );
                        
                        S.get_param_list_ptr()->show_estimates(true);
                    }
                }
            }
        }
       
        // PLOT data
        if( (bool)setup[ "tropo_time_series" ]["apply"] )
        {
            // only in case of Inequality Constrained Least Squares adjustment
            // since standard deviations are replaced by HPD intervals
            if( (bool)setup[ "icls" ]["apply"] )
            {
                if( p == 0)
                {
                    ivg::Matrix zeros( trop_delay[station].rows(),1,0.0 );
                    plot0.plot_mjd_series( trop_delay[station](":",0), zeros, 
                                     { QColor( "008941" ), 10.0, Qt::DashDotDotLine, 
                                       QCPGraph::lsLine, symbols.at( p ), version.c_str(), 0.0, 0.0});        
                }            

                if( iter->first == "lsm" || iter->first == "lsm_insitu" || iter->first == "lsm_combi_zhd" )
                {
                    plot0.plot_mjd_series_std( trop_delay[station](":",0),trop_delay[station](":",1)*1e3,trop_delay[station](":",2)*1e3,
                                     { QColor(color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                                       QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, "hh:mm:ss");
                }
                else
                {
//                    string hpd_file = (const char *)setup[ "icls" ]["hpd_file"];
//                    ivg::Matrix hpd, hpd_at;
//                    hpd.load_bin( data_path+version+"/"+dbname+"_hpd" );
//                    hpd_at = hpd.get_sub(310,0,334,1);
//                    
//                    plot0.plot_mjd_series_hpd( trop_delay[station](":",0),trop_delay[station](":",1)*1e3,
//                                               hpd_at(":",0)*1e3, hpd_at(":",1)*1e3,
//                                     { QColor(color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
//                                       QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, "hh:mm:ss"); 

                    plot0.plot_mjd_series( trop_delay[station](":",0),trop_delay[station](":",1)*1e3,
                                           { QColor(color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                                             QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, "hh:mm:ss"); 
                }
            }
            // check if there is a stochastic signal (e.g., when using a Least Squares Collocation 
            // adjustement or a filter techniqe), which has to be added to 
            // the deterministic representation of the tropospheric delay
            else if( S.get_param_list_ptr()->exist_stoch_param( ivg::paramtype::zwd, station ) )
            {
                ivg::Matrix ones( trop_pred[station].rows(),1,1.0 );
                ivg::Matrix at_delays = ones*trop_delay[station].get_col(1) + trop_pred[station].get_col(1);
                ivg::Matrix at_std = ones*trop_delay[station].get_col(2) + trop_pred[station].get_col(2);

                plot0.plot_mjd_series_std( trop_pred[station](":",0),at_delays*1e3,at_std*1e3,
                                 { QColor(color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                                   QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, "hh:mm:ss");                    

            }
            else
            {
                plot0.plot_mjd_series_std( trop_delay[station](":",0),trop_delay[station](":",1)*1e3,trop_delay[station](":",2)*1e3,
                                 { QColor(color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                                   QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, "hh:mm:ss");
            }
                
            
            plot0.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at( p ).c_str() )));
            plot0.get_plot()->graph()->setErrorBarSize(20.0);            
            plot0.set_title( "Tropospheric delays" );
            QCustomPlot *plt_ptr0 = plot0.get_plot();
            plt_ptr0->xAxis->setLabel("MJD");
            plt_ptr0->yAxis->setLabel( ( type+" [mm]").c_str() );
            QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
            plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font0.family(), 14) );
            plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font0.family(), 14) ); 
            plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font0.family(), 14) );
            plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font0.family(), 14) );
            plt_ptr0->replot();     

            p++;
        }
    }   
    fig.exec();

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
	
	
	
	
	
