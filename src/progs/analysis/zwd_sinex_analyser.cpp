/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   zwd_sinex_analyser.cpp
 * Author: corbin
 *
 * Created on February 14, 2019, 9:34 AM
 */

#include <cstdlib>
#include <iostream>
#include <string>

#include "lapack_wrapper.h"
#include "tictoc.h"
#include "session.h"
#include "session_inout.h"
#include "date.h"
#include "trf.h"
#include "station_analysis.h"
#include "lsa.h"

#include <tclap/CmdLine.h>
#include "logger.h"

#include <libconfig.h++>
#include <QApplication>
#include "plot.h"
#include "sinex.h"
#include "auxfunc.h"

#include <gsl/gsl_bspline.h>

ivg::Matrix get_spline_jacobian( const ivg::Matrix& mjd, double start, double end, int n_basisfuns, int order){
    
    int nbreak = n_basisfuns - order + 2;

    gsl_bspline_workspace* bw;
    gsl_vector *B;

    B = gsl_vector_alloc(order);

    bw = gsl_bspline_alloc(order, nbreak);
    gsl_bspline_knots_uniform(start, end, bw);
    
    std::vector<double> mjd_valid;
    for( double m : mjd ){
        if(m >= start && m <= end )
            mjd_valid.push_back(m);
    }

    ivg::Matrix A(mjd_valid.size(), n_basisfuns, 0.0);

    for( int r = 0; r < mjd_valid.size(); ++r){
        size_t istart, iend;
        gsl_bspline_eval_nonzero(mjd_valid[r], B, &istart, &iend, bw);

        for(int b = 0; b < order; ++b){
            int c = istart+b;
            A(r, c) = gsl_vector_get(B,b);
        }
    }
    
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    return A;
}

ivg::Matrix compute_clock(const ivg::Matrix& polyC, const ivg::Matrix& splineC, const ivg::Matrix& T, double start, double end, double refepoch){
    
    ivg::Matrix res(T.rows(), 1, 0.0);
    
    if( polyC.rows() > 0 ){

        ivg::Matrix dT = T - refepoch;

        ivg::Matrix A(dT.rows(), polyC.rows(), 1.0);
        for( int power = 1; power < polyC.rows(); ++power ){
            A.set_col(power, dT.pow((double)power) );
        }

        res += A*polyC;
    }

    if( splineC.rows() > 0 ){

        ivg::Matrix B = get_spline_jacobian(  T, start, end, splineC.rows(), 2);

        res += B*splineC;

    }
    
    return res;
}

void plot_hist(  ivg::Matrix& data, Plot& plot, std::string save_path, std::string histname,
                 double scale_1, double scale_2, std::string ylabel, std::string yunit_1, std::string yunit_2, std::string title, int cid){
 

    QFormat format = { Qt::black, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, "", 1.0 , 4.0  };


    std::string color_hex = ivg::get_color_value(cid);
    format.color = color_hex.c_str();
    format.legend = histname.c_str();

    plot.get_plot()->legend->setVisible(false);

//    data.save_bin( save_path + histname + ".dat");

    ivg::Matrix scaled_data = data*scale_1;

    plot.get_plot()->xAxis->setLabel( (ylabel + "(" + yunit_1 + ")").c_str());
    plot.get_plot()->yAxis->setLabel("n");
    double spacing = plot.histogram(scaled_data, sqrt(scaled_data.numel()), format,0,1);

    double mean = scaled_data.meanD();
    double sigma = scaled_data.stdD();
    ivg::Matrix keys(scaled_data.min(),(scaled_data.max()-scaled_data.min())/1000.0,scaled_data.max(),1);
    ivg::Matrix tmp =  keys-mean;
    tmp = tmp / sigma;
    tmp = tmp.pow(2);
    tmp = tmp * -0.5;
    tmp = tmp.exp();
    tmp = tmp / (sigma * 2.506628274631 );
    tmp = tmp * scaled_data.length() * spacing;

    std::ostringstream leg;
    leg << " Mean: " << fixed << setw(6) << setprecision(2) << mean << yunit_1 << " Std: " << fixed << setw(6) << setprecision(1) <<  sigma << yunit_1;
    std::cout << title << " " << leg.str() << std::endl;
    plot.plot_data(keys,tmp,  { color_hex.c_str(), 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, leg.str().c_str(), 1.0 , 4.0  });

    plot.set_title(title );
    plot.setWindowTitle(title.c_str());

    plot.get_plot()->rescaleAxes();

    // percent axis
    unsigned int n = scaled_data.length();
   
    double upperY = plot.get_plot()->yAxis->range().upper;

    plot.get_plot()->yAxis2->setVisible(true);
    plot.get_plot()->yAxis2->setRange(0, upperY);
//    plot.get_plot()->yAxis2->setLabel(" % ");

    // secondary axis
    double lowerX = plot.get_plot()->xAxis->range().lower;
    double upperX = plot.get_plot()->xAxis->range().upper;

    plot.get_plot()->xAxis2->setVisible(true);
    plot.get_plot()->xAxis2->setRange(lowerX*scale_2, upperX*scale_2);
    plot.get_plot()->xAxis2->setLabel( (ylabel + "(" + yunit_2 + ")").c_str() );

    plot.get_plot()->replot();


}


void plot_hist(  std::map< std::string, ivg::Matrix>& data, std::vector<Plot>& plots, std::string save_path,
                 double scale_1, double scale_2, std::string ylabel, std::string yunit_1, std::string yunit_2, std::string title){

    int c = 0;
    for(auto& a : data){
        plot_hist(  a.second, plots[c], save_path, a.first, scale_1,  scale_2,  ylabel,  yunit_1, yunit_2, title, c);
        ++c;
    }
      
}

void plot_hist(  std::map< std::string, ivg::Matrix>& data, std::map<std::string,Plot>& plots, std::string save_path,
                 double scale_1, double scale_2, std::string ylabel, std::string yunit_1, std::string yunit_2,
                 std::string title, int solution_id, std::string solution_name){

    for(auto& a : data){
        plot_hist(  a.second, plots[a.first], save_path, solution_name, scale_1,  scale_2,  ylabel,  yunit_1, yunit_2, title+a.first, solution_id+1);
    }
      
}


bool get_azel(std::string path, std::string station,  std::vector<double>& mjd, std::vector<double>& az, std::vector<double>& el){

    ifstream inStream_init(path.c_str(), ios::in);
    if( !inStream_init.is_open() ){
        return false;
    }else{
        string line; 
        while (getline(inStream_init,line,'\n')){
            vector<string> tokens = get_tokens(line);
            if( tokens.at(3) == station ){
                az.push_back( std::stod(tokens.at(5)) );
                el.push_back( std::stod(tokens.at(7)) );
                mjd.push_back(  std::stod(tokens.at(1)) );
            } else if( tokens.at(4) == station ){
                az.push_back( std::stod(tokens.at(6)) );
                el.push_back( std::stod(tokens.at(8)) );
                mjd.push_back(  std::stod(tokens.at(1)) );
            }
        
        }
    }
    
    return true;
    
}

using namespace libconfig;

loglevel g_verbose;
void * g_ephis;

/*
 * 
 */
int main(int argc, char** argv) {
    
    
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
   std::string controlfile = m1Arg.getValue();
   g_verbose = (loglevel) m2Arg.getValue();
   
 
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

    std::string save_path = setup["outdir"];
    std::string ReferenceSimpath = setup["Reference"]["path"];
    
    bool analyse_zwd = setup["analyse_zwd_per_session"];
    bool analyse_clo = setup["analyse_clo_per_session"];
    bool zwd_hist = setup["analyse_zwd_hist"];
    bool clo_hist = setup["analyse_clo_hist"];
    
    ivg::Simulation sim;
    
    // unsigned: solution_id    string: station    Matrix:differences
    std::map< unsigned, std::map< std::string, ivg::Matrix> > zwd_deltas;
    std::map< unsigned, std::map< std::string, ivg::Matrix> > clo_deltas;
    
    try
    {     
       
        ivg::Masterfile masterfile(setup["definitions"]["masterfiles"], ivg::mastertype::both);
        map< std::string, std::vector<unsigned> > snxfiles;
   
        std::vector<string> sinex_path;
        std::vector<string> solution_name;
        std::vector<string> obsInfo_path;
           
        for(unsigned i = 0; i < setup["sinex"].getLength(); ++i ){
            std::vector<std::string> ls = list_local_dir( (const char*)setup["sinex"][i][0]  );
            sinex_path.push_back( (const char*)setup["sinex"][i][0] );
            solution_name.push_back( (const char*)setup["sinex"][i][1] );
            if(setup["sinex"][i].getLength() > 2){
                obsInfo_path.push_back( (const char*)setup["sinex"][i][2] );
            } else {
                obsInfo_path.push_back("");
            }
            for( std::string& s : ls ){
                if( ends_with(s, ".snx") == false  )
                    continue;
                snxfiles[ s.substr(0,9) ].push_back(i);

            }
        }
        
        // loop over dbs
        QApplication fig( argc, argv );
        for(auto& session : snxfiles ){

            std::string dbname = session.first;
            
           
            std::map<std::string, Plot> clo_plots;
            std::map<std::string, Plot> zwd_plots;
            
            // init sinex files
            std::map<unsigned, ivg::Sinex> snxs;
            for( unsigned solution_id : session.second ){
                std::string infile = sinex_path[solution_id] + "/" + dbname + "_ivg2015x.snx";
                snxs[solution_id] = ivg::Sinex( infile, false );
            }

            // get ref epoch
            // ASSUMPTION: All solutions the same start epoch
            ivg::Date epoch = snxs.begin()->second.get_start_epoch();
            
            // get unique list of all stations in all solution corresponding to the same database
            std::vector<std::string> stations; 
            for( auto snx : snxs){
                std::vector<std::string> sta_names = snx.second.get_station_names();
                stations.insert(stations.end(), sta_names.begin(), sta_names.end());
            }            
            sort( stations.begin(), stations.end() );
            stations.erase( unique( stations.begin(), stations.end() ), stations.end() );     

            // get true clock values
            std::map<std::string, ivg::Matrix > true_clo;
            for(std::string& sta_name : stations){
                    ivg::Matrix tmp;
                    tmp.load_bin( ReferenceSimpath +"/" + dbname + "/" + sta_name +  "_clock.dat" );
                    true_clo[sta_name] = tmp;
            }
            
            // get true troposphere 
            std::map<std::string, ivg::HemisphereData> hemisphere_data;
            for(std::string& sta_name : stations){
                hemisphere_data[sta_name] = sim.init_hemisphere_data(setup["Reference"]["troposphere"], epoch);
//                std::cout <<  dbname << " " <<  epoch.get_date_time("YYYY-MO-DD,HH:MI:SS") << std::endl;

                ivg::Trf trf( setup, stations, ivg::staname::ivs_name, false);
                std::map<string, std::set<std::string> > twinmap = trf.get_twins_map_including_all();

                
                hemisphere_data[sta_name].set_name(sta_name);

                std::string path = ReferenceSimpath +"/" + dbname + "/" + sta_name +  "_EZWD.dat";

                bool referenceFileExists = true;
                if( !file_exists(path) ){
                    referenceFileExists = false;                    

                    for(std::string sta: twinmap[sta_name]){
                        path = ReferenceSimpath +"/" + dbname + "/" + sta +  "_EZWD.dat";
                        if ( file_exists(path) ){
                            referenceFileExists = true;
                            break;
                        }
                    }
                }

                if(referenceFileExists){
                    hemisphere_data[sta_name].load(path);
                } else {
                    throw runtime_error( " no reference troposphere found for station " + sta_name );
                }

            }
            
            //get reference clock station
            // ASSUMPTION all -solutions have the same refclock!
            //                -only one refclock
            std::set<std::string> ref_clock_stations = snxs.begin()->second.get_ref_clock_station( );
            
            if(ref_clock_stations.size() > 1){
                 log<WARNING> ("!!! MORE THAN ONE REFSTATION ");
            } else if ( ref_clock_stations.size() == 0 ){
                log<WARNING> ("!!! NO REFSTATION ");
            }
            std::string ref_clock_station = *ref_clock_stations.begin();
            log<INFO> ("*** ref station:" ) % ref_clock_station;

            // plot true data ---------------------------------------------
            if(analyse_clo){
                for(std::string& sta_name : stations){
                    ivg::Matrix mjd = hemisphere_data[sta_name].getIntervalBeginningsMJD();
                    
                    clo_plots[sta_name].get_plot()->xAxis->setLabel("time");
                    clo_plots[sta_name].get_plot()->yAxis->setLabel("clocks (ps)");
                    
                    clo_plots[sta_name].get_plot()->legend->setVisible(false);

                    QString title = "clocks " + QString::fromStdString(sta_name);
                    clo_plots[sta_name].setWindowTitle( title );

                    clo_plots[sta_name].set_title( dbname + " CLOCKS " + sta_name );

                    clo_plots[sta_name].plot_mjd_series( mjd, true_clo[sta_name]*1e+12, { Qt::gray, 1.0, Qt::DotLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, "true clock", 1.0 , 4.0  }, 0, std::string("HH:mm:ss"));
                    clo_plots[sta_name].plot_mjd_series( mjd, (true_clo[sta_name]-true_clo[ref_clock_station])*1e+12, { Qt::black, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, "true clock w.r.t. ref station", 1.0 , 4.0  }, 0, std::string("HH:mm:ss"));
    
                }
            }
            
            if(analyse_zwd){
                for(std::string& sta_name : stations){
                    
                    zwd_plots[sta_name].get_plot()->yAxis->setTickLabelFont(QFont("",14));
                    zwd_plots[sta_name].get_plot()->yAxis2->setTickLabelFont(QFont("",14));
                    zwd_plots[sta_name].get_plot()->yAxis->setLabelFont(QFont("",16));
                    zwd_plots[sta_name].get_plot()->yAxis2->setLabelFont(QFont("",16));
                    
                    
                    ivg::Matrix mjd = hemisphere_data[sta_name].getIntervalBeginningsMJD();
                    ivg::Matrix true_zwd_series = hemisphere_data[sta_name].get_data(0.0, M_PI/2.0) * ivg::c;
                    
                    // plot ezwd of all cells
                    for( unsigned cell_id = 1; cell_id < hemisphere_data[sta_name].get_number_of_cells(); ++cell_id ){
                        zwd_plots[sta_name].plot_mjd_series( mjd, hemisphere_data[sta_name].get_data(cell_id)*1000* ivg::c, 
                                                  { Qt::gray, 1.0, Qt::DotLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, "", 1.0 , 4.0  },
                                                  0, std::string("HH:mm"));
                    }
                    zwd_plots[sta_name].get_plot()->legend->clearItems();

                    zwd_plots[sta_name].get_plot()->xAxis->setLabel("");
                    zwd_plots[sta_name].get_plot()->yAxis->setLabel("zwd (mm)");
                    
                    zwd_plots[sta_name].get_plot()->legend->setVisible(false);
                    
                    QString title = "ZWD " + QString::fromStdString(sta_name);
                    zwd_plots[sta_name].setWindowTitle(title);

                    //zwd_plots[sta_name].set_title( dbname + " ZWD " + sta_name );

                    zwd_plots[sta_name].plot_mjd_series( mjd, true_zwd_series.transpose()*1000, 
                                                  { Qt::black, 2.0, Qt::DashLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, "true zwd", 1.0 , 4.0  },
                                                  0, std::string("HH:mm"));
                }
            }
            
            // get estimated data ----------------------------------------------
            
            for( unsigned solution_id : session.second){

                ivg::Sinex& snx = snxs[solution_id];
                
                if(analyse_clo || clo_hist){

                    // get clocks ---------------------------------------------------------
                    std::string clotype;
                    int clo_max_order;
                    std::map< std::string, ivg::Matrix > clo = snx.get_clocks( clotype, clo_max_order );

                    unsigned c = 0;

                    for(const auto& a : clo){

                        std::string sta_name = a.first;

                        double start = snx.get_start_epoch().get_double_mjd();
                        double end = snx.get_end_epoch().get_double_mjd();
                        
                        ivg::Matrix mjd = hemisphere_data[sta_name].getIntervalBeginningsMJD();

                        // get polynomial coefficient and cpwlf coefficient
                        ivg::Matrix polyC(0,0); // in sec / day
                        ivg::Matrix splineCmjd(0,0); 
                        ivg::Matrix splineC(0,0);

                        if( clotype == "polynomial" ){
                            polyC = a.second.get_sub( 0, 1, a.second.rows()-1, 1 );
                        } else if ( clotype == "cpwlf" && clo_max_order > 0 ) {
                            polyC = ivg::Matrix(clo_max_order+1, 1, 0.0);
                            for( int order = 1; order <= clo_max_order; ++order ){
                                polyC(order) = a.second(order-1,1);
                            }
                            splineC = a.second.get_sub( clo_max_order, 1, a.second.rows()-1, 1 );
                            splineCmjd = a.second.get_sub( clo_max_order, 0, a.second.rows()-1, 0);
                            
                            start = *splineCmjd.begin();
                            end = splineCmjd.back();
                            
                        } else {
                            splineC = a.second.get_sub( 0, 1, a.second.rows()-1, 1 );
                            splineCmjd = a.second.get_sub( clo_max_order, 0, a.second.rows()-1, 0);
                            
                            start = *splineCmjd.begin();
                            end = splineCmjd.back();
                        }
                        
                        ivg::Matrix T(start, mjd(1)-mjd(0), end, 1); // times for evaluation of polynom and cplf

    //                    polyC.show();
    //                    splineC.show();
                        
                        ivg::Matrix res = compute_clock(polyC, splineC, T, start, end, a.second(0,0) );
                        
                        ivg::Matrix clockAtKnot;
                        if( clotype == "cpwlf" ){
                            clockAtKnot = compute_clock(polyC, splineC, splineCmjd, start, end, a.second(0,0) );                            
                        }
                        ivg::Matrix diff = true_clo[sta_name].get_sub( 0, 0, res.rows()-1, 0 ) - true_clo[ref_clock_station].get_sub( 0, 0, res.rows()-1, 0 ) - res;

                        if( clo_deltas[solution_id].count(sta_name) ){
                            clo_deltas[solution_id][sta_name].append_rows( diff );
                        } else {
                            clo_deltas[solution_id][sta_name] = diff;
                        }

                        // PLOT CLO --------------------------------------------------------
                        if(analyse_clo){      
                            
                            QString label = (const char*)setup["sinex"][ solution_id ][1];
                            QColor color(ivg::get_color_value(solution_id+1).c_str()); 
                            clo_plots[sta_name].plot_mjd_series( T, res*1e+12, { color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, label , 1.0 , 4.0  }, 0, std::string("HH:mm:ss"));

                            if( clotype == "cpwlf" ){
                                clo_plots[sta_name].plot_mjd_series( splineCmjd, clockAtKnot*1e+12, { color, 2.0, Qt::NoPen, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "" , 1.0 , 4.0  }, 0, std::string("HH:mm:ss"));
                                clo_plots[sta_name].get_plot()->legend->removeItem( clo_plots[sta_name].get_plot()->legend->rowCount()-1 );
                            }
                            
                            
                            clo_plots[sta_name].get_plot()->rescaleAxes();

                            double lower = clo_plots[sta_name].get_plot()->yAxis->range().lower;
                            double upper = clo_plots[sta_name].get_plot()->yAxis->range().upper;

                            clo_plots[sta_name].get_plot()->yAxis2->setVisible(true);
                            clo_plots[sta_name].get_plot()->yAxis2->setRange(lower*ivg::c*1e-9, upper*ivg::c*1e-9);
                            clo_plots[sta_name].get_plot()->yAxis2->setLabel("clocks (mm)");

                            clo_plots[sta_name].get_plot()->replot();
                        }

                        c++;
                    }

                }

                if(analyse_zwd || zwd_hist){

                    //key: station value: [mjd, zwd, std]
                    std::map< std::string, ivg::Matrix > zwd;

                    zwd = snx.get_tropo_delays( ivg::reftype::delta );


                    // PLOT ZWD --------------------------------------------------------
                    // loop over each station
                    unsigned c = 0;
                    for(auto& a : zwd){
                        std::string sta_name = a.first;
  
                        ivg::Matrix mjd = hemisphere_data[sta_name].getIntervalBeginningsMJD();
                        double end = snx.get_end_epoch().get_double_mjd();
                        
                        ivg::Matrix true_zwd(a.second.rows(), 1);
                        
                        for( unsigned i = 0; i < a.second.rows(); ++i ){
                            true_zwd(i) = hemisphere_data[sta_name].get_data(0.0, M_PI/2.0, ivg::Date( a.second(i,0) ) ) * ivg::c ;
                        }
                        
                        ivg::Matrix diff;
                        
                        // compute difference
                        std::string obsinfoPath = obsInfo_path[solution_id];

                        if(obsinfoPath.size() > 0){

                            ivg::sessinfo si = masterfile.get_session_info(dbname);

                            std::string code = si.code;
                            std::transform(code.begin(), code.end(), code.begin(), ::tolower);

                            obsinfoPath += code + "/sim" + sta_name + "ezwd.bin";

                            if( file_exists(obsinfoPath) ){

                                //std::cout << sta_name << std::endl;

                                ivg::Matrix observed_true_zwd;
                                observed_true_zwd.load_bin(obsinfoPath);
                               
                                 // plot actual observed zwd
                                 QColor color(ivg::get_color_value(solution_id+1).c_str()); 
                                if (analyse_zwd){
                                    zwd_plots[sta_name].plot_mjd_series( observed_true_zwd.get_col(0), observed_true_zwd.get_col(1)*ivg::c*1000, 
                                                              { color, 1.0, Qt::NoPen, QCPGraph::lsLine, QCPScatterStyle::ssDiamond, "", 2.0 , 2.0  },
                                                              0, std::string("HH:mm"));
                                    zwd_plots[sta_name].get_plot()->legend->removeItem( zwd_plots[sta_name].get_plot()->legend->rowCount()-1 );
                                }
                                
                                
                                ivg::Matrix SplineCoeff = zwd[sta_name].get_col(1);
                                ivg::Matrix SplineCoeffmjd = zwd[sta_name].get_col(0);


                                ivg::Matrix B = get_spline_jacobian(  observed_true_zwd.get_col(0), *SplineCoeffmjd.begin(), SplineCoeffmjd.back(), SplineCoeff.rows(), 2);

                                ivg::Matrix estimated_zwd = B*SplineCoeff;


                                diff = observed_true_zwd.get_col(1)*ivg::c  - estimated_zwd;
                                



                            } else {
                                std::cerr << "file not existent " <<  obsinfoPath << std::endl;
                                 diff = true_zwd - zwd[sta_name].get_col(1);
                            }
                        } else {
                            diff = true_zwd - zwd[sta_name].get_col(1);

                           
                        }
                        
                        if( zwd_deltas[solution_id].count(sta_name) ){
                            zwd_deltas[solution_id][sta_name].append_rows( diff );
                        } else {
                            zwd_deltas[solution_id][sta_name] = diff;
                        }
                        
                        if (analyse_zwd){

                            QString label = (const char*)setup["sinex"][ solution_id ][1];
                            QColor color(ivg::get_color_value(solution_id+1).c_str()); 
                                 
//                            zwd_plots[sta_name].plot_mjd_series ( a.second.get_col(0), a.second.get_col(1)*1000,
//                                                         {  color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc,label, 3.0 , 3.0  },
//                                                          0, std::string("HH:mm"));     
                            
                             
                            zwd_plots[sta_name].plot_mjd_series_std ( a.second.get_col(0), a.second.get_col(1)*1000, a.second.get_col(2)*1000,
                                                         {  color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc,label, 3.0 , 3.0  },
                                                          0, std::string("HH:mm"));     
                                                         
                            zwd_plots[sta_name].get_plot()->graph()->setErrorPen(QPen(color, 1.5));
                                                        

                            // change appereance of plot
                                         
                            zwd_plots[sta_name].get_plot()->rescaleAxes();
                                                        
                            double lower = zwd_plots[sta_name].get_plot()->yAxis->range().lower;
                            double upper = zwd_plots[sta_name].get_plot()->yAxis->range().upper;

                            zwd_plots[sta_name].get_plot()->yAxis2->setVisible(true);
                            zwd_plots[sta_name].get_plot()->yAxis2->setRange(lower/ivg::c*1e+9, upper/ivg::c*1e+9);
                            zwd_plots[sta_name].get_plot()->yAxis2->setLabel("zwd (ps)");
                                                        
                            zwd_plots[sta_name].get_plot()->xAxis->setRangeLower( zwd_plots[sta_name].get_plot()->xAxis->range().lower - 60.0 );
                                 
                            zwd_plots[sta_name].get_plot()->replot();
                        }
                        c++;
                    }

                }
            
            }
            
            if (analyse_zwd || analyse_clo){
                fig.exec();
            }
        }
     
        // station name, plot
        std::map<std::string,Plot> zwd_hists, clo_hists;
        
        std::set<string> stations_with_zwd;
        for(auto& a : zwd_deltas){
            for(auto& s :  a.second){
                stations_with_zwd.insert(s.first);
            }
        }
        for(const std::string& s: stations_with_zwd ){
            zwd_hists[s];
        }
        
        
        std::set<string> stations_with_clo;
        for(auto& a : clo_deltas){
            for(auto& s :  a.second){
                stations_with_clo.insert(s.first);
            }
        }
        for(const std::string& s : stations_with_clo ){
            clo_hists[s];
        }
        
        if(zwd_hist){ 
            for(unsigned solution_id = 0; solution_id < solution_name.size(); ++solution_id){
                std::cerr << "histogram zwd" << solution_name[solution_id] << std::endl;
                plot_hist( zwd_deltas[solution_id],  zwd_hists, save_path + "/zwd_delta_", 
                           1.0e12/ivg::c, 1e-9*ivg::c, "true - estimated ", "ps", "mm", "ZWD ", solution_id, solution_name[solution_id] );
            }
        }
        
        if(clo_hist){ 
            for(unsigned solution_id = 0; solution_id < solution_name.size(); ++solution_id){
                std::cerr << "histogram clo" << solution_name[solution_id] << std::endl;
                plot_hist( clo_deltas[solution_id],  clo_hists, save_path + "/clo_delta_", 
                           1.0e12, 1e-9*ivg::c, "true - estimated ", "ps", "mm", "CLOCKS ", solution_id, solution_name[solution_id] );
            }
        }
        
                
//        std::map <unsigned, std::vector<Plot> > zwd_hists2, clo_hists2;
//        if(zwd_hist){ 
//            for(unsigned solution_id = 0; solution_id < solution_name.size(); ++solution_id){
//                zwd_hists2[solution_id].resize( zwd_deltas[solution_id].size() );
//                plot_hist( zwd_deltas[solution_id],  zwd_hists2[solution_id], save_path + "/zwd_delta_", 
//                           1.0e12/ivg::c, 1e-9*ivg::c, "true - estimated ", "ps", "mm", "ZWD " + solution_name[solution_id]);
//            }
//        }
//        
//        if(clo_hist){
//            for(unsigned solution_id = 0; solution_id < solution_name.size(); ++solution_id){
//                clo_hists2[solution_id].resize( clo_deltas[solution_id].size() );
//                plot_hist( clo_deltas[solution_id],  clo_hists2[solution_id], save_path + "/clo_delta_",
//                           1.0e12, 1e-9*ivg::c, "true - estimated ", "ps", "mm", "clock " + solution_name[solution_id]);
//            }
//        }
        
        if( zwd_hist || clo_hist ){
            
            fig.exec();
        }
        
    } 
   catch( libconfig::SettingException &err )   {
      cerr << "libconfig::SettingException: " << err.what() << " at " << err.getPath() << endl;
      exit( -1 );
   } catch(std::exception& e) {
      cerr << "std::exception: " << e.what();
   }
    
    return 0;
}

