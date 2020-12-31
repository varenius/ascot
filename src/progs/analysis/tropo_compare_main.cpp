/***************************************************************************** 
 * independent VLBI solution                                                 *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 * 2016-07-01 - TS                                                           *
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
#include "tsa.h"
#include <regex>

using namespace libconfig;

loglevel g_verbose;
void * g_ephis;

bool contains(const Setting &setup, string str1, string str2);
bool contains(const Setting &setup, string str1);
int find(const Setting &setup, string str1, string str2, string cmp);
std::string ReplaceString(std::string subject, const std::string& search, const std::string& replace);

int main(int argc, char *argv[]) {

    // .........................................................................
    // handle commandline arguments
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");

    // name of config file
    TCLAP::ValueArg<std::string> m1Arg("c", "config", "Name of the control file", true, "", "string");
    cmd.add(m1Arg);

    // verbose level
    TCLAP::ValueArg<int> m2Arg("v", "verbose", "Verbose Level (0=NOTHING, 1=INFO, 2=DETAIL, 3=RESULT, 4=WARNING, 5=ALL)", false, 0, "int");
    cmd.add(m2Arg);

    // parse commandline arguments
    cmd.parse(argc, argv);
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

    try {
        cfg.readFile(controlfile.c_str());
    } catch (libconfig::ParseException & err) {
        cerr << "libconfig::" << err.what() << " in " << err.getFile() << " at line " << err.getLine() << endl;
        exit(-1);
    }

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


    QApplication fig(argc, argv);
    Plot plot0;
    Plot plotDiff;
    Plot plotDiffBias;
    Plot plotBiasReduced;
    Plot plotAllan;
    Plot plotAllan2;
    Plot plotMaxCorr;
    Plot plotPowerSp;
    map< std::string, ivg::Matrix > compare;
    int p = 0;

    std::string type_var = "zwd";

    std::string stuff = "compare";
    Setting& setup = cfg.lookup("setup");
    std::string typeParam = (const char *) setup[ "tropo_time_series" ]["type"];
    std::string stationName = (const char *) setup[ "tropo_time_series" ]["station"];

    std::string time_start = (const char *) setup[ "tropo_time_series" ]["time_start"];
    std::string time_end = (const char *) setup[ "tropo_time_series" ]["time_end"];
    ivg::Date d_start(time_start, "YYMMMDD");
    ivg::Date d_end(time_end, "YYMMMDD");


    for (int i = 0; i < setup.lookup( stuff ).getLength(); i++) {
        std::string name = setup.lookup( stuff )[ i ][ "name" ];
        
        std::string path = setup.lookup( stuff )[ i ][ "path" ];

        path= ReplaceString(path, "$station$",stationName);

        
        std::string type = setup.lookup( stuff )[ i ][ "type" ];
        std::string apriori_type = setup.lookup( stuff )[ i ][ "apriori_type" ];
        double unit_factor;
        stringstream stream(setup.lookup( stuff )[ i ][ "unit_factor" ]);
        stream >> unit_factor;
        string local_ties_from = setup.lookup( stuff )[ i ][ "local_ties_from" ];
        string local_ties_to = setup.lookup( stuff )[ i ][ "local_ties_to" ];

//        double d_zwd = 0;
                ivg::Matrix d_zwd(0,0,0);
        ivg::Matrix data;
        
        cerr << endl << endl << endl<<"starte Einlesen für "<< name<<" "<< typeParam<<" " << type<<endl;
        
        if (type == "cfg") {

        } else if (type == "bin") {

            data.load_bin(path);
//                        d_zwd = ivg::Matrix(data.size(1),1,0);
            
        } else if (type == "snx_tro" || type == "resultfile"  || type == "snx") {
            
            std::string ending = setup.lookup( stuff )[ i ][ "ending" ];
            std::string statName = setup.lookup( stuff )[ i ][ "statName" ];
            
            std::string gpspath = path + "/";
            std::string pathrest = "(.+)("  + ending + ")";
            cerr<<"Location: "<<gpspath<<", regexp: "<<pathrest<<endl;
            
            ivg::Matrix wholeSolution;
            map< std::string, ivg::Matrix > trop_delay;
            std::map<std::string, ivg::Matrix >::iterator mapSessionsIter;
            mapSessionsIter = trop_delay.begin();
            
            int k = 0;

            DIR *dir;
            struct dirent *ent;
            const char* gpspathchr = gpspath.c_str();
            if ((dir = opendir(gpspathchr)) != NULL) {
                while ((ent = readdir(dir)) != NULL) {

                    std::regex e(pathrest);

                    if (std::regex_match(ent->d_name, e)) {
                        string pathi = gpspath + ent->d_name;
                        if (file_exists(pathi)) {
//                            cerr<<"pathi:"<<pathi<<endl;
                            std::string station = stationName;
                            boost::to_upper(station);
                                
                            if (type == "snx_tro") {


                                string test;
                                if (typeParam == "zwd") {
                                    test = "TROWET";
                                } else if (typeParam == "zhd") {
                                    test = "TROTOT";
                                }


                                ivg::Sinex snx = ivg::Sinex(pathi, test);
//                                cerr << "SINEX files einlesen fertig:" << endl;

                                if (typeParam == "zwd") {
                                    trop_delay = snx.get_tropo_delays(ivg::reftype::estimate);
                                } else if (typeParam == "zhd") {
                                    trop_delay = snx.get_tropo_delays(ivg::reftype::apriori);
                                }


                            } else if (type == "snx") {

                                ivg::Sinex snx = ivg::Sinex(pathi, false);
//                                cerr << "SINEX_TRO files einlesen fertig:" << endl;

                                if (typeParam == "zwd") {
                                    trop_delay = snx.get_tropo_delays(ivg::reftype::delta);
                                } else if (typeParam == "zhd") {
                                    trop_delay = snx.get_tropo_delays(ivg::reftype::estimate);
                                } else if (typeParam == "ztd") {
                                    trop_delay = snx.get_tropo_delays(ivg::reftype::apriori);
                                }
                                
                            } 
//                            else if (type == "resultfile") {
//
//                                ivg::Session_inout sin;
//                                ivg::Session session;
//                                sin.read_results(&session, pathi); // wenn _order -1, dann stochastischer Zuschlag. Dieser auch immer pro beobachtung.
//
//                                // Wenn stochastische Zuschläge vorhanden, immer einlesen und darstellen (Zuschlag plus nearest neighbor offset)
//                                // Wenn nicht vorhanden, nur Offsets darstellen.
//                                //
//                                // Sessions können ausfallen. Oder auch ganze Session oder Teile enthalten falsche Luftdruckwerte wegen Sensorausfall.
//                                // Möglichkeit 1: Zeitlicher Schwellwert (0.5 mjd) für nearest neighbor einführen bzw. Sessionzugehörigkeit prüfen,
//                                // ob Offset nicht zu falscher Session gehört.
//                                // Möglichkeit 2: Überprüfung der Funktionswerte. Wenn Sessionübergänge großen Sprung aufweisen. Robuste Funktion reinlegen, Ausreißerdetektion.
//                                // 
//
//                                if (session.get_param_list_ptr()->exist_stoch_param(ivg::paramtype::zwd, station)) {
//
//                                    if (typeParam == "zwd") {
//                                        trop_delay = ivg::Sinex::get_stoch_tropo_delays(ivg::reftype::estimate, session.get_param_list_ptr()->begin(), session.get_param_list_ptr()->end());
//                                    } else if (typeParam == "zhd") {
//                                        trop_delay = ivg::Sinex::get_stoch_tropo_delays(ivg::reftype::apriori, session.get_param_list_ptr()->begin(), session.get_param_list_ptr()->end());
//                                    }
//                                } else {
//                                    if (typeParam == "zwd") {
//                                        trop_delay = ivg::Sinex::get_tropo_delays(ivg::reftype::estimate, session.get_param_list_ptr()->begin(), session.get_param_list_ptr()->end());
//                                    } else if (typeParam == "zhd") {
//                                        trop_delay = ivg::Sinex::get_tropo_delays(ivg::reftype::apriori, session.get_param_list_ptr()->begin(), session.get_param_list_ptr()->end());
//                                    }
//                                }
//                            }
                            
//                            cerr<<"Eingelesene Daten "<<trop_delay.size()<<endl;
//                            for(map<std::string, ivg::Matrix >::const_iterator it = trop_delay.begin();it != trop_delay.end(); it++)
//                            {
//                                std::cerr << ">"<<it->first << "< " << it->second.size(1) << " "<< it->second.size(2)<<" "  <<
//                                       " "<<   ivg::Date((it->second)(":",0).min()).get_date_time("DD/MON/YYYY")<<" " <<ivg::Date(it->second.max()).get_date_time("DD/MON/YYYY")<< "\n";
//                            }
//                            cerr << endl;

                            if ((trop_delay[station]).size(1) != 0 && (trop_delay[station]).size(2) != 0) {
                                if (k == 0) {
                                    wholeSolution = ivg::Matrix(trop_delay[station]);
                                } else
                                    wholeSolution.append_rows(trop_delay[station]);
                                k++;
                            }else{
                                if (k == 0) {
                                    wholeSolution = ivg::Matrix(trop_delay[""]);
                                } else
                                    wholeSolution.append_rows(trop_delay[""]);
                                k++;
                            }
                        }
                    }

                    //   

                }
                closedir(dir);
            } else {
                /* could not open directory */
                perror("direrror");
                return EXIT_FAILURE;
            }
                            
            data = wholeSolution;
        }
   
        cerr<<endl;
        cerr<<"DATEN EINGELESEN: Größe Daten: "<<data.size(1)<<"x"<<data.size(2)<<endl;
        cerr<<"Von: "<<ivg::Date(data(":",0).min()).get_date_time("DD/MON/YYYY")<<" bis " <<ivg::Date(data(":",0).max()).get_date_time("DD/MON/YYYY")<<endl;
        cerr<<endl;
        
        
        // Zeitbereich ausschneiden
        cerr<<"ZEITBEREICH AUSSCHNEIDEN: "<<endl;
        d_start.get_date_time("DD/MON/YYYY");
        d_end.get_date_time("DD/MON/YYYY");
        
        vector<int> idx;
        ivg::Matrix t_bereich = data(":", 0).find_elem(le, d_end.get_double_mjd() + 1, ge, d_start.get_double_mjd(), idx);
        std::vector<int> ix(data.size(2)-1);
        std::iota(ix.begin(), ix.end(), 1);
        ivg::Matrix(ix).disp();
        t_bereich.append_cols((data(":", ix)).get_rows(idx));
        data = t_bereich;
        cerr<<"Größe Daten: "<<data.size(1)<<"x"<<data.size(2)<<endl;
        cerr<<endl;
        
        
        // Sortieren
        data.sort_cols(0);
        

        d_zwd = ivg::Matrix(data.size(1), 1, 0);
        int loc = -1;
        if ((bool)setup[ "tropo_time_series" ]["applyTies"]) {
            if (setup.exists(local_ties_to) && setup.exists(local_ties_from)) {


                double type = (double) setup.lookup( local_ties_from )["lat"];
                double local_ties_from_lat, local_ties_from_lon, local_ties_from_h;
                local_ties_from_lat = (double) setup.lookup( local_ties_from )["lat"];
                local_ties_from_lon = (double) setup.lookup( local_ties_from )["lon"];
                local_ties_from_h = (double) setup.lookup( local_ties_from )["h_ell"];


                double local_ties_to_lat, local_ties_to_lon, local_ties_to_h;
                local_ties_to_lat = (double) setup.lookup( local_ties_to )["lat"];
                local_ties_to_lon = (double) setup.lookup( local_ties_to )["lon"];
                local_ties_to_h = (double) setup.lookup( local_ties_to )["h_ell"];


                std::string meteo_source = setup.lookup( stuff )[ i ][ "meteo_source" ];
                loc = find(setup, stuff, "name", meteo_source);
                if (loc != -1) {
                    std::string meteo_path = setup.lookup( stuff )[ loc ][ "path" ];
                    ivg::Matrix meteo;
                    meteo_path = ReplaceString(meteo_path, "$station$", stationName);
                    meteo.load_bin(meteo_path);

                    int mitte = floor(data.size(1) / 2);

                    // TODO: efficient Interpolation for long time series (out of memory for >10years)

                    //        double p0=(meteo(":",3).interpolate(meteo(":",0),data(mitte, 0),"linear"))(0,0);

                    //      int startR = floor(meteo.size(1)*19/20);
                    //      int endR = meteo.size(1)-1;
                    //      cerr<<startR<<" "<<endR<<" "<<meteo(startR,0)<<" "<<meteo(endR,0)<<" "<<data(mitte, 0)<<endl;
                    //        double p0=(meteo.get_sub(startR, 3, endR, 3 ) .interpolate(meteo.get_sub(startR, 0, endR, 0 ),data(mitte, 0),"linear"))(0,0);
                    //        double t0=(meteo.get_sub(startR, 2, endR, 2 ) .interpolate(meteo.get_sub(startR, 0, endR, 0 ),data(mitte, 0),"linear"))(0,0);
                    //        cerr<<p0<<" "<<t0<<" "<<endl;



                    //         vector<int> idx_meteo;
                    //         ivg::Matrix meteo_bereich = meteo(":",0).find_elem(lt,data(data.size(1)-1, 0),ge,data(0, 0),idx_meteo );
                    //         int ixA = idx_meteo(0,0);
                    //         (data(":",0)+(ixA-meteo(ixA,0))).round();


                    // TODO: Noch wird ein konstanter Wert für den gesamten Bereich abgezogen
                    //                    double p0 = (meteo(":", 3).interpolate(meteo(":", 0), data(mitte, 0), "linear"))(0, 0);
                    //                    double t0 = (meteo(":", 2).interpolate(meteo(":", 0), data(mitte, 0), "linear"))(0, 0);
                    //                    double e0 = 4.76;
                    //                    d_zwd = ivg::Troposphere::calc_tropospheric_ties_zwd(local_ties_from_h, local_ties_to_h, local_ties_to_lat, p0, t0, e0);



                    for (int ii = 0; ii < data.size(1); ii++) {
                        double p0 = (meteo(":", 3).interpolate(meteo(":", 0), data(ii, 0), "linear"))(0, 0);
                        double t0 = (meteo(":", 2).interpolate(meteo(":", 0), data(ii, 0), "linear"))(0, 0);
                        double e0 = 4.76;


                        if (apriori_type == "ZTD") {
                            double davis = ivg::Troposphere::calc_zhd_davis(p0, local_ties_to_lat * M_PI / 180.0, local_ties_to_h);
                            d_zwd(ii, 0) = ivg::Troposphere::calc_tropospheric_ties_zwd(local_ties_from_h, local_ties_to_h, p0, t0, e0) / 1000.0 - davis;
                        } else {
                            d_zwd(ii, 0) = ivg::Troposphere::calc_tropospheric_ties_zwd(local_ties_from_h, local_ties_to_h, p0, t0, e0) / 1000.0;
                        }
                        //        
                        //        data(ii,1) -= ivg::Troposphere::calc_tropospheric_ties_zwd( local_ties_from_h, local_ties_to_h, local_ties_to_lat,p0, t0,  e0);
                        //            data(ii,1) -= d_zwd;
                        //                
                        ////        ivg::Matrix delta = ivg::Troposphere.calc_tropospheric_ties( local_ties_from_h, local_ties_to_h, local_ties_to_lat,p0, t0,  e0,"ZWD");
                    }
                }
            } else { // still not reduced for d_zhd_davis

                // TODO: use flag in cfg to perform ties or not. In case of not still need to reduce totals or sth. else.

                //            std::string meteo_source = setup.lookup( stuff )[ i ][ "meteo_source" ];
                //            loc = find(setup, stuff, "name", meteo_source);
                //            cerr << "aa" << endl;
                //            if (loc != -1) {
                //                std::string meteo_path = setup.lookup( stuff )[ loc ][ "path" ];
                //                ivg::Matrix meteo;
                //                meteo_path = ReplaceString(meteo_path, "$station$", stationName);
                //                meteo.load_bin(meteo_path);
                //
                //                if (apriori_type == "ZTD") {
                //                    cerr << "still not reduced for d_zhd_davis -> REDUCE" << endl;
                //                    for (int ii = 0; ii < data.size(1); ii++) {
                //                        double p0 = (meteo(":", 3).interpolate(meteo(":", 0), data(ii, 0), "linear"))(0, 0);
                //                        double t0 = (meteo(":", 2).interpolate(meteo(":", 0), data(ii, 0), "linear"))(0, 0);
                //                        double e0 = 4.76;
                //                        //                            std::cerr<<ii<<std::endl;
                //
                //
                //                        double davis = ivg::Troposphere::calc_zhd_davis(p0, local_ties_to_lat * M_PI / 180.0, local_ties_to_h);
                //                        data(ii, 1) = data(ii, 1) * unit_factor - davis;
                //                    }
                //                }
                //            }
            }
        }
        if (data.size(1) > 0) {
            for (int j = 0; j < data.size(1); j++) {
                data(j, 1) = data(j, 1) * unit_factor + d_zwd(j, 0);
            }
        }
        cerr << "EINLESEN BEENDET" << endl;

        



        compare[name] = data;

        if (contains(setup, "time_series_analysis", "plot")) {

            std::cerr << "Starte plotten für " <<name.c_str()<< std::endl;

            // plot in mm
            plot0.plot_mjd_series(data(":", 0), data(":", 1)*1e3, {
                QColor(color_values.at(p).c_str()), 1.0, Qt::SolidLine,
                QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
            }, 0, "hh:mm:ss");
            p++;

            if (contains(setup, "time_series_analysis", "diff_biases")) {


                ivg::Matrix data_t_biases, data_m_biases;
                ivgat::Tsa t1;


                // nur 1 Tag!!!
                t1.intervalMeanValues(data(":", 0), data(":", 1), 1, data_t_biases, data_m_biases, ivgat::stepMode::daily);
//
//                plot0.plot_mjd_series(data_t_biases, data_m_biases * 1e3,{
//                    QColor(color_values.at(p).c_str()), 1.0, Qt::DotLine,
//                    QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
//                }, 0);
//                p++;



                int k = 0;
                ivg::Matrix time_bias_reduced(data.size(1), 1, 0);
                ivg::Matrix data_bias_reduced(data.size(1), 1, 0);
                for (int ii = 0; ii < data.size(1); ii++) {
                    if (data(ii, 0) >= data_t_biases(0) && data(ii, 0) <= data_t_biases(data_t_biases.length() - 1)) {
                        data_bias_reduced(k, 0) = data(ii, 1)-(data_m_biases.interpolate(data_t_biases, data(ii, 0), "linear"))(0);
                        time_bias_reduced(k, 0) = data(ii, 0);
                        k++;
                    }
                }
                data_bias_reduced = data_bias_reduced.get_sub(0, 0, k - 1, 0);
                time_bias_reduced = time_bias_reduced.get_sub(0, 0, k - 1, 0);

                plotBiasReduced.plot_mjd_series(time_bias_reduced, data_bias_reduced * 1e3,{
                    QColor(color_values.at(p).c_str()), 1.0, Qt::SolidLine,
                    QCPGraph::lsLine, QCPScatterStyle::ssDot, (name + " bias_red").c_str(), 1.0, 12.0
                }, 0);
                p++;
            }
        }

    }

    if (contains(setup, "time_series_analysis", "plot")) {
        plot0.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plot0.get_plot()->graph()->setErrorBarSize(20.0);
        plot0.set_title("Tropospheric delays");
        QCustomPlot *plt_ptr0 = plot0.get_plot();
        plt_ptr0->xAxis->setLabel("MJD");
        plt_ptr0->yAxis->setLabel((type_var + " [mm]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->replot();
    }
    
    if (contains(setup, "time_series_analysis", "diff_biases")) {
        plotBiasReduced.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plotBiasReduced.get_plot()->graph()->setErrorBarSize(20.0);
        plotBiasReduced.set_title("Tropospheric delays bias reduced");
        QCustomPlot *plt_ptr0 = plotBiasReduced.get_plot();
        plt_ptr0->xAxis->setLabel("MJD");
        plt_ptr0->yAxis->setLabel((type_var + " [mm]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->replot();
    }

    if (contains(setup, "time_series_analysis", "allan")) {
        int p=0;
        for (int i = 0; i < setup.lookup( stuff ).getLength(); i++) {
            std::string name = setup.lookup( stuff )[ i ][ "name" ];
            std::cerr << "Starte allan plot für " <<name.c_str() << std::endl;

            ivg::Matrix data = compare[name];

            ivg::Matrix freq      = data(":",1).diff();
            ivg::Matrix time_diff = data(":",0).diff();
            freq.div_elem( time_diff );
            ivg::Matrix t2 = data(":",0).get_sub( 0,0,data(":",0).length()-2,0 );
      
            ivg::Matrix  avar;
            ivg::Matrix tau;
            
            ivgat::Tsa ts1;

//            t2.save_bin("t2.bin");
//            freq.save_bin("freq.bin");
//            ts1.calcTwoSampleAllanVarianceFrequency(t2,freq,tau, avar);
            ts1.calcTwoSampleAllanVarianceFrequency(data(":",0),data(":",1),tau, avar); // so ist es richtig

            
            // Variance to standard deviation
            avar = avar.sqrt();
            
            plotAllan.plot_data(tau, avar,{
                QColor(color_values.at(p ).c_str()), 1.0, Qt::SolidLine,
                QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
            }, 0);
            p++;
            
            
//   
//            ivg::Matrix  avarPhase;
//            ivg::Matrix tauPhase;
////            data(":",0).save_bin("t.bin");
////            data(":",1).save_bin("obs.bin");
//            ts1.calcTwoSampleAllanVariancePhase(data(":",0),data(":",1),tauPhase, avarPhase);
//            
//            
//            // Already standard deviation
//            // Variance to standard deviation
//            avarPhase = avarPhase.sqrt();
//            
//            plotAllan.plot_data(tauPhase, avarPhase,{
//                QColor(color_values.at(p ).c_str()), 1.0, Qt::DotLine,
//                QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
//            }, 0);
//            p++;
            
        }

        
        plotAllan.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plotAllan.get_plot()->graph()->setErrorBarSize(20.0);
        plotAllan.set_title("Allan standard deviation");
        QCustomPlot *plt_ptr0 = plotAllan.get_plot();
        plt_ptr0->xAxis->setLabel("tau");
        plt_ptr0->yAxis->setLabel((""+type_var + " Allan standard deviation [mm]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        
        plt_ptr0->yAxis->setScaleType(QCPAxis::stLogarithmic);
        plt_ptr0->yAxis->setScaleLogBase(10);
        plt_ptr0->xAxis->setScaleType(QCPAxis::stLogarithmic);
        plt_ptr0->xAxis->setScaleLogBase(10);
        plt_ptr0->replot();
    }
    
    
    
    if (contains(setup, "time_series_analysis", "allan_of_diff")) {
        p = 0;
        ivg::Matrix data_ref = compare[setup.lookup( stuff )[ 0 ][ "name" ]];
        string ref_name = (setup.lookup( stuff )[ 0 ][ "name" ]);

        for (int i = 1; i < setup.lookup( stuff ).getLength(); i++) {
            std::string name = setup.lookup( stuff )[ i ][ "name" ];
            std::cerr << "Starte diff plot für " <<name.c_str()<<" w.r.t " <<ref_name <<std::endl;

            ivg::Matrix data = compare[name];
            ivg::Matrix interpForDiff(0, 2, 0);

            if (data.rows() != 0 && data.cols() != 0) {
                //            int k = 0;
                for (int j = 0; j < data_ref.size(1); j++) {

                    if (data_ref(j, 0) > data(":", 0).min() && data_ref(j, 0) < data(":", 0).max()) {
                        ivg::Matrix interp = data(":", 1).interpolate(data(":", 0), data_ref(j, 0), "linear");
                        ivg::Matrix append(1, 2, 0.0);
                        append(0, 0) = data_ref(j, 0);
                        append(0, 1) = interp(0, 0) - data_ref(j, 1);
                        interpForDiff.append_rows(append);
                    }


                }

                

            ivg::Matrix freq      = interpForDiff(":",1).diff();
            ivg::Matrix time_diff = interpForDiff(":",0).diff();
            freq.div_elem( time_diff );
            ivg::Matrix t2 = interpForDiff(":",0).get_sub( 0,0,interpForDiff(":",0).length()-2,0 );
      
            ivg::Matrix  avar;
            ivg::Matrix tau;
            
            ivgat::Tsa ts1;

//            t2.save_bin("t2.bin");
//            freq.save_bin("freq.bin");
            ts1.calcTwoSampleAllanVarianceFrequency(t2,freq,tau, avar);

            
            // Variance to standard deviation
            avar = avar.sqrt();
            
            plotAllan2.plot_data(tau, avar,{
                QColor(color_values.at(p ).c_str()), 1.0, Qt::SolidLine,
                QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
            }, 0);
            p++;
            
   
            ivg::Matrix  avarPhase;
            ivg::Matrix tauPhase;

//            data(":",0).save_bin("t.bin");
//            data(":",1).save_bin("obs.bin");
            ts1.calcTwoSampleAllanVariancePhase(interpForDiff(":",0),interpForDiff(":",1),tauPhase, avarPhase);
            
            
            // Already standard deviation
            // Variance to standard deviation
            avarPhase = avarPhase.sqrt();
            
            plotAllan2.plot_data(tauPhase, avarPhase,{
                QColor(color_values.at(p ).c_str()), 1.0, Qt::DotLine,
                QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
            }, 0);
            p++;
            
//            plotDiff.plot_mjd_series(interpForDiff(":", 0), interpForDiff(":", 1)*1e3, {
//                    QColor(color_values.at(p + 1).c_str()), 1.0, Qt::SolidLine,
//                    QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
//                }, 0);
//            }

            p++;
        }
        }

        
        plotAllan2.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plotAllan2.get_plot()->graph()->setErrorBarSize(20.0);
        plotAllan2.set_title("Allan standard deviation");
        QCustomPlot *plt_ptr0 = plotAllan2.get_plot();
        plt_ptr0->xAxis->setLabel("tau");
        plt_ptr0->yAxis->setLabel((""+type_var + " Allan standard deviation [mm]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        
        plt_ptr0->yAxis->setScaleType(QCPAxis::stLogarithmic);
        plt_ptr0->yAxis->setScaleLogBase(10);
        plt_ptr0->xAxis->setScaleType(QCPAxis::stLogarithmic);
        plt_ptr0->xAxis->setScaleLogBase(10);
        plt_ptr0->replot();
    }



    if (contains(setup, "time_series_analysis", "diff")) {
        p = 0;
        ivg::Matrix data_ref = compare[setup.lookup( stuff )[ 0 ][ "name" ]];
        string ref_name = (setup.lookup( stuff )[ 0 ][ "name" ]);

        for (int i = 1; i < setup.lookup( stuff ).getLength(); i++) {
            std::string name = setup.lookup( stuff )[ i ][ "name" ];
            std::cerr << "Starte diff plot für " <<name.c_str()<<" w.r.t " <<ref_name <<std::endl;

            ivg::Matrix data = compare[name];
            ivg::Matrix interpForDiff(0, 2, 0);

            if (data.rows() != 0 && data.cols() != 0) {
                //            int k = 0;
                
                if(data.rows() == data_ref.rows() && ((data(":",0)-data_ref(":",0)).sum_col())(0)==0.0){
                        
                    cerr<<"x_data identical"<<endl;
                    interpForDiff = data(":",0);
                    interpForDiff.append_cols(data(":",1)-data_ref(":",1));
                            
                } else {

                    for (int j = 0; j < data_ref.size(1); j++) {

                        if (data_ref(j, 0) > data(":", 0).min() && data_ref(j, 0) < data(":", 0).max()) {
                            
                            
                            std::vector<int> found ;
                            found = data(":", 0).find_idx( eq, data_ref(j, 0)) ;
                            if(found.size()>0){
                                ivg::Matrix append(1, 2, 0.0);
                                append(0, 0) = data_ref(j, 0);
                                append(0, 1) = data(found.at(0), 1) - data_ref(j, 1);
                                interpForDiff.append_rows(append);
                            } else {

                                ivg::Matrix interp = data(":", 1).interpolate(data(":", 0), data_ref(j, 0), "linear");
                                ivg::Matrix append(1, 2, 0.0);
                                append(0, 0) = data_ref(j, 0);
                                append(0, 1) = interp(0, 0) - data_ref(j, 1);
                                interpForDiff.append_rows(append);
                            }
                        }


                    }
                }

                plotDiff.plot_mjd_series(interpForDiff(":", 0), interpForDiff(":", 1)*1e3, {
                    QColor(color_values.at(p + 1).c_str()), 1.0, Qt::SolidLine,
                    QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
                }, 0);
            }

            p++;
        }

        plotDiff.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plotDiff.get_plot()->graph()->setErrorBarSize(20.0);
        plotDiff.set_title("Tropospheric delays differences");
        QCustomPlot *plt_ptr0 = plotDiff.get_plot();
        plt_ptr0->xAxis->setLabel("MJD");
        plt_ptr0->yAxis->setLabel(("Delta " + type_var + " [mm]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->replot();
    }
    
    
    if (contains(setup, "time_series_analysis", "diff_biases")) {
        p = 0;
        ivg::Matrix data_ref = compare[setup.lookup( stuff )[ 0 ][ "name" ]];
        
//        ivg::Matrix data_ref_t_biases,data_ref_m_biases;
//        ivgat::Tsa t1;
//        t1.intervalMeanValues(data_ref(":",0),data_ref(":",1), 1,data_ref_t_biases,data_ref_m_biases, ivgat::stepMode::even_mjd);

        for (int i = 1; i < setup.lookup( stuff ).getLength(); i++) {
            std::string name = setup.lookup( stuff )[ i ][ "name" ];
            std::cerr << "Starte diff_biases plot für " <<name.c_str() << std::endl;

            ivg::Matrix data = compare[name];
            
            
            
            ivg::Matrix interpForDiff(0, 2, 0);

            if (data.rows() != 0 && data.cols() != 0) {
                
                ivgat::Tsa t1;
                ivg::Matrix data_t_biases,data_m_biases;
                t1.intervalMeanBiases(data_ref(":",0),data_ref(":",1),data(":",0),data(":",1), 1,data_t_biases,data_m_biases,ivgat::stepMode::daily);

                
                //            int k = 0;
//                for (int j = 0; j < data_ref.size(1); j++) {
//
//                    
//                    if (data_ref(j, 0) > data(":", 0).min() && data_ref(j, 0) < data(":", 0).max()) {
//                        ivg::Matrix interp = data(":", 1).interpolate(data(":", 0), data_ref(j, 0), "linear");
//                        ivg::Matrix append(1, 2, 0.0);
//                        append(0, 0) = data_ref(j, 0);
//                        append(0, 1) = interp(0, 0) - data_ref(j, 1);
//                        interpForDiff.append_rows(append);
//                    }
//
//
//                }

                plotDiffBias.plot_mjd_series(data_t_biases, data_m_biases*1e3, {
                    QColor(color_values.at(p + 1).c_str()), 1.0, Qt::SolidLine,
                    QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
                }, 0);
                
                
                if (contains(setup, "time_series_analysis", "plot")) {
                    int k=0;
                    ivg::Matrix time_bias_reduced(data.size(1),1,0);
                    ivg::Matrix data_bias_reduced(data.size(1),1,0);
                    for(int i=0;i<data.size(1);i++){
                        if(data(i,0)>=data_t_biases(0) && data(i,0)<=data_t_biases(data_t_biases.length()-1)){
                            data_bias_reduced(k,0) = data(i,1)-(data_m_biases.interpolate(data_t_biases, data(i,0), "linear"))(0);
                            time_bias_reduced(k,0) = data(i,0);
                            k++;
                        }
                    }
                    data_bias_reduced = data_bias_reduced.get_sub(0,0,k-1,0);
                    time_bias_reduced = time_bias_reduced.get_sub(0,0,k-1,0);
                    
                    plot0.plot_mjd_series(time_bias_reduced, data_bias_reduced*1e3, {
                        QColor(color_values.at(p).c_str()), 1.0, Qt::SolidLine,
                        QCPGraph::lsLine, QCPScatterStyle::ssDot, (name+" bias_red").c_str(), 1.0, 12.0
                    }, 0);
//                    p++;
            
                }
            }

            p++;
        }
        
        

        plotDiffBias.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plotDiffBias.get_plot()->graph()->setErrorBarSize(20.0);
        string ref_name = (setup.lookup( stuff )[ 0 ][ "name" ]);
        plotDiffBias.set_title("Tropospheric delays biases w.r.t. " +ref_name);
        QCustomPlot *plt_ptr0 = plotDiff.get_plot();
        plt_ptr0->xAxis->setLabel("MJD");
        plt_ptr0->yAxis->setLabel(("Delta " + type_var + " [mm]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->replot();
    }


    if (contains(setup, "time_series_analysis", "max_corr")) {
        
        // TODO: crosscorrelation wird für nichtsymmetrische Abstände berechnet. 
        // Achse geht also von -x1 bis x2.
        
        p = 0;
        ivg::Matrix data_ref = compare[setup.lookup( stuff )[ 0 ][ "name" ]];
        double min_ref = data_ref(0, 0);
        double max_ref = data_ref(data_ref.size(1) - 1, 0);
        double max_corr = (max_ref - min_ref); ///10;
//        cerr << "min_ref" << min_ref << " " << max_ref << " " << max_corr << " " << endl;

        for (int i = 1; i < setup.lookup( stuff ).getLength(); i++) {
//            std::cerr << "Starte diff plot" << std::endl;
            std::string name = setup.lookup( stuff )[ i ][ "name" ];

            ivg::Matrix data = compare[name];
            //            int k = 0;

            ivgat::Tsa tsa1;
            ivg::Matrix tout;
            ivg::Matrix mout;
//            std::cerr << "Start crosscorr" << std::endl;
            int numVal = 1000;
            int start = 30000;



            //            tsa1.crosscorrNEQD2(data_ref(":",0), data_ref(":",1),data(":",0), data(":",1),100,tout ,  mout);
            //                tsa1.crosscorrNEQD2(data_ref.get_sub(start, 0, start + numVal, 0), data_ref.get_sub(start, 1, start + numVal, 1), data.get_sub(start, 0, start + numVal, 0), data.get_sub(start, 1, start + numVal, 1), 100, tout, mout);
            vector<int> idx;
//            cerr << "start_ber" << endl;
            ivg::Matrix t_bereich = data(":", 0).find_elem(le, max_ref + max_corr, ge, min_ref - max_corr, idx);
//            cerr << "ber" << t_bereich.size(1) << " gr " << t_bereich(0) << " " << data_ref(0, 0) << " idx " << idx.at(1) << endl;
            int rows = data.rows();
            //                std::for_each(idx.begin(), idx.end(), [rows](int &d) { d+=rows;});// nächste Spalte
//            cerr << "idx" << idx.at(1) << " " << idx.size() << " " << (data(":", 0))(idx.at(0)) << endl;
//            cerr << "t1" << data_ref(0, 0) << " " << data_ref(data_ref.size(1) - 1, 0) << " t2 " << t_bereich(0, 0) << " " << t_bereich(t_bereich.size(1) - 1, 0) << endl;
            tsa1.crosscorrNEQD2(data_ref(":", 0), data_ref(":", 1), t_bereich, (data(":", 1))(idx), 11, tout, mout);
//            std::cerr << "End crosscorr" << std::endl;
            int maxIndex;
            double maxValue = mout.max(maxIndex);
            plotMaxCorr.plot_data(tout, mout,{
                QColor(color_values.at(p + 1).c_str()), 1.0, Qt::SolidLine,
                QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
            }, 0);
            p++;
        }

        plotMaxCorr.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plotMaxCorr.get_plot()->graph()->setErrorBarSize(20.0);
        plotMaxCorr.set_title("Tropospheric delays crosscorrelation");
        QCustomPlot *plt_ptr0 = plotMaxCorr.get_plot();
        plt_ptr0->xAxis->setLabel("delta time");
        plt_ptr0->yAxis->setLabel(("CrossCorr " + type_var + " [-]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->replot();
    }

    if (contains(setup, "time_series_analysis", "powerspektrum")) {
        p = 0;
        ivg::Matrix data_ref = compare[setup.lookup( stuff )[ 0 ][ "name" ]];

        //        for (int i = 1; i < setup.lookup( stuff ).getLength(); i++) {
        //            std::cerr << "Starte diff plot" << std::endl;
        //            std::string name = setup.lookup( stuff )[ i ][ "name" ];
        //            std::cerr << "line " << i << std::endl;
        //
        //            ivg::Matrix data = compare[name];
        ////            int k = 0;
        //            
        //            ivgat::Tsa tsa1;
        //            ivg::Matrix tout;
        //            ivg::Matrix mout;
        //            std::cerr<<"Start powerspektrum"<<std::endl;
        //            int numVal = 5000;
        ////            tsa1.crosscorrNEQD2(data_ref(":",0), data_ref(":",1),data(":",0), data(":",1),100,tout ,  mout);
        //            tsa1.crosscorrNEQD2(data_ref.get_sub(0,0,numVal,0), data_ref.get_sub(0,1,numVal,1),data.get_sub(0,0,numVal,0), data.get_sub(0,1,numVal,1),100,tout ,  mout);
        //            std::cerr<<"End powerspektrum"<<std::endl;
        //            int maxIndex;
        //            double maxValue = mout.max(maxIndex);
        //            plotPowerSp.plot_data(tout, mout, {
        //                QColor(color_values.at(p).c_str()), 1.0, Qt::SolidLine,
        //                QCPGraph::lsLine, QCPScatterStyle::ssDot, name.c_str(), 1.0, 12.0
        //            }, 0);
        //            p++;
        //        }

        plotPowerSp.get_plot()->graph()->setErrorPen(QPen(QColor(color_values.at(p).c_str())));
        plotPowerSp.get_plot()->graph()->setErrorBarSize(20.0);
        plotPowerSp.set_title("Tropospheric delays powerspektrum");
        QCustomPlot *plt_ptr0 = plotPowerSp.get_plot();
        plt_ptr0->xAxis->setLabel("MJD");
        plt_ptr0->yAxis->setLabel(("CrossCorr " + type_var + " [-]").c_str());
        QFont font0 = plt_ptr0->xAxis->selectedLabelFont();
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atLeft)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->axisRect(0)->axis(QCPAxis::atBottom)->setTickLabelFont(QFont(font0.family(), 14));
        plt_ptr0->replot();
    }


    fig.exec();

    

    return ( 0);
}

bool contains(const Setting &setup, string str1, string str2) {
    bool b = false;
    for (int i = 0; i < setup.lookup( str1 ).getLength(); i++) {
        std::string name = setup.lookup( str1 )[ i ];
        if (name == str2) {
            b = true;
            break;
        }
    }
    return b;
}

int find(const Setting &setup, string str1, string str2, string cmp) {
    int b = -1;
    for (int i = 0; i < setup.lookup( str1 ).getLength(); i++) {
        std::string name = setup.lookup( str1 )[ i ].lookup(str2);
        if (name == cmp) {
            b = i;
            break;
        }
    }
    return b;
}

bool contains(const Setting &setup, string str1) {
    bool b = false;

    for (int i = 0; i < setup.getLength(); i++) {

        Setting &params = setup[ i ];
        std::string name = params.getName();
//        std::cerr << name << std::endl;
        if (name == str1) {
            b = true;
            break;
        }
    }
    return b;
}


std::string ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
    return subject;
}
