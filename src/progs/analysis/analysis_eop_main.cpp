/* 
 * File:   analysis_eop_main.cpp
 * Author: TA
 *
 * Created on 12. Oktober 2015, 09:40
 */

#include <cstdlib>
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

using namespace libconfig;

loglevel g_verbose;
void * g_ephis;

using namespace std;
/*
 * 
 */
int main(int argc,char** argv)
{
    // .........................................................................
    // handle commandline arguments
    TCLAP::CmdLine cmd("Command description message",' ',"0.9");

    // name of config file
    TCLAP::ValueArg<std::string> m1Arg("c","config","Name of the control file",true,"","string");
    cmd.add(m1Arg);

    // verbose level
    TCLAP::ValueArg<int> m2Arg("v","verbose","Verbose Level (0=NOTHING, 1=INFO, 2=DETAIL, 3=RESULT, 4=WARNING, 5=ALL)",false,0,"int");
    cmd.add(m2Arg);

    // parse commandline arguments
    cmd.parse(argc,argv);
    string controlfile = m1Arg.getValue();
    g_verbose = (loglevel) m2Arg.getValue();

    // .........................................................................
    // initialize timer    
    tictoc timeCount;
    timeCount.tic();

    std::cout<<" \n =================================================================== "
             <<std::endl;
    std::cout<<"     > > > > > ivg::ASCOT (anaylsis tools) < < < < <     "
             <<std::endl<<std::endl;


    // loading controlfile
    Config cfg;

    try
    {
        cfg.readFile(controlfile.c_str());
    }
    catch (libconfig::ParseException & err)
    {
        cerr<<"libconfig::"<<err.what()<<" in "<<err.getFile()<<" at line "
            <<err.getLine()<<endl;
        exit(-1);
    }

    vector<std::string> color_values = {
        "#000000", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    };
    
    try
    {
        Setting& setup = cfg.lookup("setup");

        // initialization of session_inout with defined session_type
        string session_type = (const char *) get_list_element(setup["datadirs"],setup["session_type"])[1];
        ivg::Session_inout sessionizer(session_type); // NGS / SNX

        //string session_type_ref = (const char *) get_list_element(setup["datadirs"],setup["session_type_ref"])[1];   ivg::Session_inout sessionizer_ref(session_type_ref); // NGS / SNX   

        //double th = (double)*setup["th_sta"]; 

        // create map containing different solutions with all databases
        std::map< std::string,std::vector<std::string> > sols;
        for (int i = 0; i<setup[ "sessions" ].getLength(); ++i)
        {
            std::string versions = setup[ "sessions" ][ i ][ "version" ];
            std::string version;
            istringstream iss(versions);
            std::string dbname = setup[ "sessions" ][ i ][ "dbname" ];

            while (std::getline(iss,version,','))
                sols[ version ].push_back(dbname);
        }

        typedef std::map<std::string,std::vector<std::string> >::iterator mapIter;
        QApplication fig(argc,argv);
        Plot plot0;        
        int subplot1,subplot2,subplot;
        int p = 0;

        // .........................................................................
        // loop over solutions
        for (mapIter iter = sols.begin(); iter!=sols.end(); iter++)
        {
            //ivg::Param_list * params_ref;

            cerr<<"START SOLUTION: "<<iter->first<<" including "
                <<iter->second.size()<<" dbs."<<endl;
            int n = iter->second.size();
       
            ivg::Matrix eops(n,6);
            //ivg::Matrix eops0(n,6);

            vector<int> outliers;
            for (int i = 0; i<n; ++i)
            {
                //cerr<<"********************************************"<<endl;

                std::string version = iter->first;
                std::string dbname = iter->second.at(i);

                cerr<<"*** Using "<<dbname<<" and "<<version<<endl;
                // load information from each database into session
                string data_path = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[2];
                std::string path = data_path+version+"/"+dbname+".snx";

                ivg::Sinex snx( path, false );

                ivg::Param_list params(snx.get_parameter_vetor());
                int idx = params.get_index(ivg::paramtype::xpo,"EOP");

                ivg::Param *para = params.get_param(idx);
                ivg::Date epo = para->get_epoch();
                double val = para->get_estimate();

                                
                eops(i,0)  = epo.get_decimal_date();
                eops(i,1)  = val*ivg::rad2mas;

                para = params.get_param(params.get_index(ivg::paramtype::ypo,"EOP"));
                eops(i,2) = para->get_estimate()*ivg::rad2mas;
                para = params.get_param(params.get_index(ivg::paramtype::ut1,"EOP"));
                eops(i,3) = para->get_estimate()*ivg::rad2s*1e-3;
                para = params.get_param(params.get_index(ivg::paramtype::nutx,"EOP"));
                eops(i,4) = para->get_estimate()*ivg::rad2mas;
                para = params.get_param(params.get_index(ivg::paramtype::nuty,"EOP"));
                eops(i,5) = para->get_estimate()*ivg::rad2mas;
                
                for(int j=1;j<4;++j)
                {
                    if(abs(eops(i,j))>1)
                    {
                        outliers.push_back(i);
                        break;
                    }
                }
            }
            eops.rem_r(outliers);   
            
            //Plot Data  
            plot0.plot_data(eops(":",0),eops(":",1),
                            {QColor(color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "X", 3.0, 3.0});
            plot0.set_title( "EOP w.r.t. IERS 08C04" );
            QCustomPlot *plt_ptr = plot0.get_plot();
            plt_ptr->xAxis->setLabel("");
            plt_ptr->yAxis->setLabel("X-pole [mas]");
            plt_ptr->legend->setVisible(false);
            QFont font = plt_ptr->xAxis->selectedLabelFont();
            plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) ); 
            plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setRange(-1.0,1.0);
            
            if(iter==sols.begin())
                subplot1 = plot0.add_rowplot("Y-pole [mas]","");    
            plot0.plot_data(eops(":",0),eops(":",2), 
                             {QColor(color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "X", 3.0, 3.0},
                             subplot1);
            plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );  
            plt_ptr->axisRect(subplot1)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(subplot1)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setRange(-1.0,1.0);
            if(iter==sols.begin())
                subplot2 = plot0.add_rowplot("UT1-TAI [ms]","MJD");    
            plot0.plot_data(eops(":",0),eops(":",3), 
                             {QColor(color_values.at( p ).c_str() ), 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "X", 3.0, 3.0},
                             subplot2); 
            plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setLabelFont( QFont (font.family(), 14) );    
            plt_ptr->axisRect(subplot2)->axis(QCPAxis::atLeft)->setTickLabelFont( QFont (font.family(), 14) );
            plt_ptr->axisRect(subplot2)->axis(QCPAxis::atBottom)->setTickLabelFont( QFont (font.family(), 14) ); 
            plt_ptr->axisRect(0)->axis(QCPAxis::atLeft)->setRange(-1.0,1.0);
            
            plot0.get_plot()->replot();
            p++;
            
            cerr << "#sessions " << eops.rows() << endl;
            for(int j=1;j<6;++j)
            {
                ivg::Matrix v = eops(":",j);
                double rms = sqrt((v.transpose()*v)(0)/(double)v.rows());
                cerr << j << " " << rms <<  endl;
            }
        }
        fig.exec();

    }
    catch (libconfig::SettingException &err)
    {
        cerr<<"libconfig::SettingException: "<<err.what()<<" at "<<err.getPath()<<endl;
        exit(-1);
    }
    catch (std::exception& e)
    {
        cerr<<"std::exception: "<<e.what();
    }

    std::cerr<<" =================================================================== \n"
            <<"Time required: "<<fixed<<timeCount.toc()<<" seconds."<<std::endl;

    // ============== ENDE 


    return 0;
}

