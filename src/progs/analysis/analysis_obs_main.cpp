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
   std::cout << "         A N A L Y S I S   O B S E R V A T I O N S " 
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
    
    typedef std::map<std::string, std::vector<std::string> >::iterator mapIter;
    std::string version, dbname;
    double stax;
    QApplication fig( argc, argv );
    Plot plot0, plot1, plot2;
    int p = 0;   
    int c = 0;
    
    // .........................................................................
    // loop over solutions
    for (mapIter iter = sols.begin(); iter != sols.end(); iter++)
    {

        log<RESULT>("*** Solution: '") % iter->first % "' including " % iter->second.size() % " dbs.";
        
        int n = iter->second.size();
        ivg::Matrix vf(n,1,0.0);
        ivg::Matrix wrms(n,1,0.0);
        ivg::Matrix rms(n,1,0.0);
        ivg::Matrix mjd(n,1,0.0);
 
        for( int i=0; i < n; ++i )
        {    
            version = iter->first;
            dbname = iter->second.at(i); 

            // load information from each database into session
            string data_path = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[2];
            std::string path = data_path+version+"/"+dbname+".snx";
            
            // get statistics (variance factor, WRMS, RMS)
            if( file_exists( path ) )
            {
                ivg::Sinex snx( path, false );              
                vf(i,0) = snx.get_statistics("VF"); 
                wrms(i,0) = snx.get_statistics("WRMS"); 
                rms(i,0) = snx.get_statistics("RMS"); 
                mjd(i,0) = snx.get_trf( ivg::reftype::estimate ).get_reference_epoch().get_double_mjd();
            }
        }
       
        // remove null elements (if sinex file does not exist)
        std::vector<int> null_idx = wrms.find_idx( lt, 1e-15 );
        if( null_idx.begin() != null_idx.end() )
        {
            vf.rem_r( null_idx );
            wrms.rem_r( null_idx );
            rms.rem_r( null_idx );
            mjd.rem_r( null_idx );            
        }
        
        // rescale (W)RMS
        wrms *= 1e12;
        rms *= 1e12;                       
       
        // calculate mean and median value and standard deviation for all statistical numbers
        double mean_vf = vf.meanD();
        double mean_wrms = wrms.meanD();
        double mean_rms = rms.meanD();

        double std_vf = vf.stdD();
        double std_wrms = wrms.stdD();
        double std_rms = rms.stdD();
        
        double median_vf = vf.median();
        double median_wrms = wrms.median();
        double median_rms = rms.median();     

        log<RESULT>("*** Solution: '") % version % "' >>> mean variance factor: " % mean_vf  
                                       % "   |    mean WRMS: " % mean_wrms % "   |    mean RMS: " % mean_rms ;
        log<RESULT>("*** Solution: '") % version % "' >>> median variance factor: " % median_vf  
                                       % "   |    median WRMS: " % median_wrms % "   |    median RMS: " % median_rms;        
        log<RESULT>("*** Solution: '") % version % "' >>> std. variance factor: " % std_vf  
                                       % "   |    std. WRMS: " % std_wrms % "   |    std. RMS: " % std_rms;
        
//        cerr << "variance factor - Solution: '" << version << " " << mean_vf << " | " << std_vf << endl;
//        cerr << "WRMS - Solution: '" << version << " " << mean_wrms << " | " << std_wrms << endl;        
//        cerr << "RMS - Solution: '" << version << " " << mean_rms << " | " << std_rms << endl; 
                
        // PLOT data
        if( c == 0 )
        {
            ivg::Matrix ones( vf.rows(),1,1.0 );
            plot0.plot_mjd_series( mjd, ones, 
                             { QColor( "008941" ), 5.0, Qt::SolidLine, 
                               QCPGraph::lsLine, symbols.at( p ), version.c_str(), 0.0, 0.0}, 0, "yyyy-MM-dd");        
        }
        
        plot0.plot_mjd_series( mjd,vf,
                         { QColor(ivg::color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                           QCPGraph::lsLine, symbols.at( p ), version.c_str(), 4.0, 3.0}, 0, "yyyy-MM-dd");
        plot0.set_title( "Chi-Squared Variance factor" );
        QCustomPlot *plt_ptr0 = plot0.get_plot();
        plt_ptr0->yAxis->setLabel("variance factor [-]");
        //plt_ptr0->legend->setVisible(false);
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font0.family(), 14) );
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font0.family(), 14) ); 
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font0.family(), 14) );
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font0.family(), 14) );
        plt_ptr0->replot();     

//        if( c == sols.size()-1 )
//        {
//            ivg::Matrix ones( vf.rows(),1,1.0 );
//            plot0.plot_mjd_series( mjd, ones, 
//                             { QColor( "008941" ), 5.0, Qt::SolidLine, 
//                               QCPGraph::lsLine, symbols.at( p ), version.c_str(), 0.0, 0.0}, 0, "yyyy-MM-dd");        
//        }        
              
        plot1.plot_mjd_series( mjd, wrms, 
                         { QColor(ivg::color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                           QCPGraph::lsLine, symbols.at( p ), version.c_str(), 4.0, 3.0}, 0, "yyyy-MM-dd");
        plot1.set_title( "Weighted Root Mean Square error" );
        QCustomPlot *plt_ptr1 = plot1.get_plot();
        plt_ptr1->yAxis->setLabel("WRMS [ps]");
        //plt_ptr1->legend->setVisible(false);
        QFont font1 = plt_ptr1->xAxis->selectedLabelFont();
        plt_ptr1->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font1.family(), 14) );
        plt_ptr1->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font1.family(), 14) ); 
        plt_ptr1->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font1.family(), 14) );
        plt_ptr1->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font1.family(), 14) );
        plt_ptr1->replot();        
        
        plot2.plot_mjd_series( mjd, rms, 
                         { QColor(ivg::color_values.at( p ).c_str() ), 1.0, Qt::SolidLine, 
                           QCPGraph::lsLine, symbols.at( p ), version.c_str(), 4.0, 3.0}, 0, "yyyy-MM-dd");
        plot2.set_title( "Root Mean Square error" );
        QCustomPlot *plt_ptr2 = plot2.get_plot();
        plt_ptr2->yAxis->setLabel("RMS [ps]");
        //plt_ptr2->legend->setVisible(false);
        QFont font2 = plt_ptr2->xAxis->selectedLabelFont();
        plt_ptr2->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font2.family(), 14) );
        plt_ptr2->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font2.family(), 14) ); 
        plt_ptr2->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font2.family(), 14) );
        plt_ptr2->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font2.family(), 14) );
        plt_ptr2->replot();          
       
        p++;
        c++;
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
	
	
	
	
	
