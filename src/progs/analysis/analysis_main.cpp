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
   std::cout << "         S T A T I O N     A N A L Y S I S   "            
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

   // initialization of session_inout with defined session_type
   string session_type = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[1];
   ivg::Session_inout sessionizer(session_type); // NGS / SNX

   // get analysis tools setup parameters
   std::string sta = (const char *)setup[ "station" ]["plot_sta_pos"];
   std::string type = (const char *)setup[ "station" ]["type"];

   std::string sta_std = (const char *)setup[ "std_niveau" ]["plot_sta_pos"];
   std::string type_std = (const char *)setup[ "std_niveau" ]["type"];   
   
   double max_est = (double)setup[ "bl_rep" ]["max_est"];
   int min_bl_no = (int)setup[ "bl_rep" ]["min_bl_no"];
   std::string rms_type = (const char *)setup[ "bl_rep" ]["type"];
   
   std::vector<std::string> rm_bls;
   for( int j = 0; j < setup[ "bl_rep" ]["exclude"].getLength(); j++ )
       rm_bls.push_back( (const char *)setup[ "bl_rep" ]["exclude"][j] );
  
    const std::vector<QCPScatterStyle::ScatterShape> symbols = {
            QCPScatterStyle::ssDisc, QCPScatterStyle::ssTriangle,
            QCPScatterStyle::ssCross, QCPScatterStyle::ssDiamond, 
            QCPScatterStyle::ssSquare,QCPScatterStyle::ssStar,        
            QCPScatterStyle::ssDisc, QCPScatterStyle::ssTriangle,
            QCPScatterStyle::ssCross, QCPScatterStyle::ssDiamond, 
            QCPScatterStyle::ssSquare,QCPScatterStyle::ssStar,
    };   
   

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
     
    
    typedef std::map<std::string, std::vector<std::string> >::iterator mapIter;
    std::map< std::string, vector<ivg::Trf> > dtrf;
    std::map< std::string, ivg::Matrix > dbl_rep;
    
    QApplication fig( argc, argv );
    Plot plot0, plot1, plot2;
    int subplot1, subplot2;
    int p = 0;
    int c1 = 0;
    std::string version, ver, dbname;
    ivg::Matrix lvl_ref;
    ivg::Matrix Y_ref; 

    // .........................................................................
    // loop over solutions
    for (mapIter iter = sols.begin(); iter != sols.end(); iter++)
    {
        vector<ivg::Trf> trfs_est, trfs_apr;
        std::map< std::string, ivg::Matrix > bl_reps;

        log<RESULT>("*** Solution: '") % iter->first % "' including " % iter->second.size() % " dbs.";
        
        int n = iter->second.size();

        for( int i=0; i < n; ++i )
        {    
            ivg::Trf trf_est; 
            ivg::Trf trf_apr;

            version = iter->first;
            dbname = iter->second.at(i); 

            // load information from each database into session
            string data_path = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[2];
            std::string path = data_path+version+"/"+dbname+".snx";

            if( file_exists( path ) )
            {
                ivg::Sinex snx( path, false );
                trf_est = snx.get_trf( ivg::reftype::estimate );
                trf_apr = snx.get_trf( ivg::reftype::apriori );

                // fill vector of TRF and Date pointer
                trfs_est.push_back( trf_est );
                trfs_apr.push_back( trf_apr );   
            }
        }    

        // initialize STATION ANALYSIS TOOLS
        ivgat::Station_analysis sta_at( trfs_est );

        // check stations
        int no_sta = sta_at.get_session_number();
        sta_at.check_stations( trfs_apr, max_est );
        log<NOTHING>("*** Solution: '") % iter->first % ": #sessions: " % no_sta 
               % ", #sessions with max. estimate value of " % max_est % "[m]: " % sta_at.get_session_number();
          
        std::vector<std::string> tmp = sta_at.get_sessions();
        
//        string dir = setup[ "outdir" ];
//        std::string sessions = dir+"/"+iter->first+"_sessions.txt";
//        ofstream newfile( sessions.c_str() );        
//        for (string &s : tmp)
//        {
//           newfile << s << endl;
//        }

  
        // handle baseline lengths repeatabilities
        if( (bool)setup[ "bl_rep" ]["apply"] )
        {
            // calculate baseline lengths
            sta_at.calc_baseline_length();
            //sta_at.show_baselines();

            // calculate baseline lengths repeatabilities
            ivg::Matrix bl_rep;
            std::vector<std::string> bl_name;
            sta_at.calc_bl_rep( bl_name, bl_rep, min_bl_no, rms_type );
            // sta_at.show_bl_repeatabilities();

            // remove selected baselines
            sta_at.remove_baselines( rm_bls );


            typedef std::map<std::string, ivg::Matrix>::const_iterator mapIter; 
            bl_reps = sta_at.get_bl_rep();
            ivg::Matrix X( bl_reps.size(), 1, 0.0 );
            ivg::Matrix Y( bl_reps.size(), 1, 0.0 );
            std::vector<std::string> bls;
            int counter = 0;

            for (mapIter it = bl_reps.begin(); it != bl_reps.end(); it++)
            {
               bls.push_back( it->first );
               X(counter,0) = it->second(0,0);
               Y(counter,0) = it->second(0,1);
               counter++;
            }
            
            if( (bool)setup[ "bl_rep" ]["ref"] )
            {   
                if( p == 0 )
                    Y_ref = Y;
                
                Y -= Y_ref;   
                             
                // PLOT baseline repeatabilities                 
                plot0.plot_data( X*1e-3, Y*1e3, 
                                 { QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
                                   QCPGraph::lsImpulse, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, bls );
                plot0.set_title( "differences in baseline length repeatabilities" );
                QCustomPlot *plt_ptr = plot0.get_plot();
                plt_ptr->xAxis->setLabel("baseline length [km]");
                plt_ptr->yAxis->setLabel( ("differences in "+rms_type+" [mm]").c_str() );
                //plt_ptr->legend->setVisible(false);
                QFont font = plt_ptr->xAxis->selectedLabelFont();
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) ); 
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->replot();                                  
            }
            else
            {               
                // fit trend to baseline lengths repeatabilities
                std::map< std::string, ivg::Matrix > fit;
                ivg::Matrix fit_mat;
                fit = sta_at.fit_bl_rep( fit_mat );                
                
                // PLOT baseline repeatabilities                 
                plot0.plot_data( X*1e-3, Y*1e3, 
                                 { QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
                                   QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 12.0}, 0, bls );
                plot0.plot_data( X*1e-3, fit_mat(":",1)*1e3, 
                                 { QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
                                   QCPGraph::lsLine, QCPScatterStyle::ssDisc, "", 0.0, 0.0 });
                plot0.set_title( "baseline length repeatabilities" );
                QCustomPlot *plt_ptr = plot0.get_plot();
                plt_ptr->xAxis->setLabel("baseline length [km]");
                plt_ptr->yAxis->setLabel( (rms_type+" [mm]").c_str() );
                //plt_ptr->legend->setVisible(false);
                QFont font = plt_ptr->xAxis->selectedLabelFont();
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) ); 
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->replot();                          
            }
  

            // compare two versions/solutions and calculate the percentage
            // of baselines which have been improved, degradated by at least 1mm
            // or which remain unchanged
            if( (bool)setup[ "cmp_sol" ]["apply"] && 
                ( version == (const char *)setup[ "cmp_sol" ]["versions"][0] || 
                  version == (const char *)setup[ "cmp_sol" ]["versions"][1] ) )
            {                     
                double amount = (double)setup[ "cmp_sol" ]["amount"];
                
                if( c1 == 0 )
                {
                    dbl_rep = bl_reps;   
                    ver = version;
                }
                else if( c1 == 1 )
                {
                    sta_at.calc_RMS_differences( dbl_rep, ver, version, amount );
                }
                c1++;
            }             
        } // bl_rep
            
        
        // handle station time series      
        if( (bool)setup[ "station" ]["apply"] )
        {
            sta_at.calc_station_positions( trfs_apr, max_est );
            
            ivg::Matrix pos; 
            ivg::Matrix pos_std; 
            ivg::Matrix epoch;
            sta_at.get_station_position( sta, type, pos, pos_std, epoch );

            if( pos.rows() != 0 && pos.cols() != 0 )
            {
                pos *= 1e3;      
                plot1.plot_data( epoch.transpose(), pos.transpose()(":",0), 
                                 {QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
                                  QCPGraph::lsLine, QCPScatterStyle::ssDisc, version.c_str(), 3.0, 3.0});
                plot1.set_title( "Station positions" );
                QCustomPlot *plt_ptr = plot1.get_plot();
                plt_ptr->yAxis->setLabel( "UP est. [mm]" );
                plt_ptr->xAxis->setLabel( "" );

                plt_ptr->legend->setVisible(false);
                QFont font = plt_ptr->xAxis->selectedLabelFont();
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) ); 
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                    

                if(iter==sols.begin()) 
                    subplot1 = plot1.add_rowplot("East est. [mm]","");   

                plot1.plot_data( epoch.transpose(), pos.transpose()(":",1),
                                 {QColor(ivg::color_values.at( p ).c_str()), 2.0, Qt::SolidLine, 
                                  QCPGraph::lsLine, QCPScatterStyle::ssDisc, version.c_str(), 3.0, 3.0}, subplot1);
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );  
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                    

                if(iter==sols.begin())
                    subplot2 = plot1.add_rowplot("North est. [mm]","MJD");   

                plot1.plot_data( epoch.transpose(), pos.transpose()(":",2),
                                 {QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
                                  QCPGraph::lsLine, QCPScatterStyle::ssDisc, version.c_str(), 3.0, 3.0}, subplot2);     

                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );    
                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                                       

                plot1.get_plot()->replot();
            }
            else
                throw runtime_error( "Selected IVS session is not part of the database. \n" );
        }

        // handle niveau of uncertainties (in terms of standard deviations of XYZ station positions)
        if( (bool)setup[ "std_niveau" ]["apply"] )
        {
            ivg::Matrix lvl_std, epo; 
            sta_at.get_level_of_uncertainties( sta_std, type_std, lvl_std, epo );

            if( (bool)setup[ "std_niveau" ]["lvl_ref"] )
            {   
                if( p == 0 )
                    lvl_ref = lvl_std;
                
                lvl_std -= lvl_ref;   
            }
            
            if( lvl_std.rows() != 0 && lvl_std.cols() != 0 )
            {            
                lvl_std *= 1e3;      
                plot2.plot_mjd_series( epo, lvl_std.transpose()(":",0), 
                                 {QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
                                  QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 8.0}, 0, "yyyy-MM-dd");
                plot2.set_title( "level of uncertainties" );
                QCustomPlot *plt_ptr = plot2.get_plot();
                plt_ptr->yAxis->setLabel( "std.dev. X [mm]" );
                plt_ptr->xAxis->setLabel( "" );

                //plt_ptr->legend->setVisible(false);
                QFont font = plt_ptr->xAxis->selectedLabelFont();
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) ); 
                plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                    

                if(iter==sols.begin()) 
                    subplot1 = plot2.add_rowplot("std.dev. Y [mm]","");   

                plot2.plot_mjd_series( epo, lvl_std.transpose()(":",1),
                                       {QColor(ivg::color_values.at( p ).c_str()), 2.0, Qt::SolidLine, 
                                        QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 8.0}, subplot1, "yyyy-MM-dd");
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );  
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                    

                if(iter==sols.begin())
                    subplot2 = plot2.add_rowplot("std.dev. Z [mm]","");   

                plot2.plot_mjd_series( epo, lvl_std.transpose()(":",2),
                                       {QColor(ivg::color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, 
                                        QCPGraph::lsNone, QCPScatterStyle::ssDisc, version.c_str(), 1.0, 8.0}, subplot2, "yyyy-MM-dd");     

                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );    
                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
                plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );                                       

                plot2.get_plot()->replot();
            }
            else
                throw runtime_error( "Selected IVS session is not part of the database. \n" );
        }      
         
        p++;          
        
        // save trf vector for each version
        dtrf[ version ] = trfs_est;             
    }   
    fig.exec();

    // handle TRF comparison
    if( (bool)setup[ "cmp_trf" ]["apply"] )
    {       
        // calculate Helmert parameter between two trf realizations, defined as
        // reference and comparison trf in config file 
        if( dtrf.find( (const char *)setup[ "cmp_trf" ]["ref"] ) != dtrf.end() && 
            dtrf.find( (const char *)setup[ "cmp_trf" ]["cmp"] ) != dtrf.end() )
        {
            std::vector<ivg::Trf> ref_trf = dtrf.find( (const char *)setup[ "cmp_trf" ]["ref"] )->second;
            std::vector<ivg::Trf> cmp_trf = dtrf.find( (const char *)setup[ "cmp_trf" ]["cmp"] )->second;

            ivgat::Station_analysis ato1( ref_trf );              
            vector<ivgat::t_param> trans_params = { ivgat::t_param::tx, ivgat::t_param::rx, 
                                                    ivgat::t_param::ry, ivgat::t_param::rz, ivgat::t_param::s };

            ato1.calc_helmert_params( cmp_trf, trans_params );
        }
        else
            throw runtime_error( "One or more TRF realization(s) selected for the calculation "
                    "of the Helmert parameters are not selected as versions in the "
                    "session defintion. \n" );    
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
	
	
	
	
	
