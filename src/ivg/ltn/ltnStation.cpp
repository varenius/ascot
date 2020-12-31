/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Station.cpp
 * Author: armin
 * 
 * Created on 9. Januar 2016, 13:51
 */

#include "ltnStation.h"
#include "date.h"
#include "vgosdb.h"
#include "auxfunc.h"
#include "logger.h"

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <map>

#include "boost/format.hpp"
#include <boost/regex.hpp>
#include <netcdfcpp.h>

namespace ltn{

//consturctor -----------------------------------------------------------------
LtnStation::LtnStation(std::string logfile, std::map<std::string, std::string> &nscodes,
                 std::map<std::string, std::vector<std::string>> &regExp,
                 std::map<std::string, std::vector<double>> &pressCor,
                 std::map<std::string, int> &sig) {
    
    this->logfile = logfile;
    
    // 2-letter code
    int pos = logfile.find(".log");
    this->name2 = logfile.substr(pos-2,2);
//    std::cout << "2-letter code: " << this->name2;
    
    // 8-letter code
    if( nscodes.count(this->name2) == 1 ){
        this->name8 = nscodes[this->name2];
//        std::cout << "\t8-letter code: " << this->name8 << endl;
    }
    else{
        log<WARNING>( "!!! ") % this->name2 % " IS NOT in the file ./Data/nscodes ! ";

    }
    // Regular Expression
    if( regExp.count(this->name2) == 1 ){
        std::string d, w, c;
        d = regExp[this->name2][0];
        w = regExp[this->name2][1];
        c = regExp[this->name2][2];
        
        this->wx = boost::regex(d+w);
        this->cb = boost::regex(d+c);
        
        this->readable = true;

    }
    else{
        log<WARNING> ("!!! ") % this->name2 % " IS NOT in the file ./Data/regex ! ";
        this->readable = false;

    }
    
    //Pressure Correction
    if( pressCor.count(this->name2) == 1 ){
        std::vector<double> d = pressCor[this->name2];
        this->gamma = d[0];
        this->g = d[1];
        this->h = d[2];
        this->h0 = d[3];
        this->pressCor = true;
        log<INFO> ( "*** \tusing Pressure Correction for: " ) % this->name8;
        //std::cout << this->gamma << " " << this->g << " " << this->h << " " << this->h0 << endl;
    }
    else{
        this->pressCor = false;
    }

    // Cable Sign
    if( sig.count(this->name2) == 1 ){

        this->sign = sig[this->name2];    
        
    }
    else{
        log<WARNING> ("!!! ") % this->name2 % " IS NOT in the file ./Data/sign ! assuming +";
        this->sign = 1; 

    }    
    

}


// Getter Setter --------------------------------------------------------------

 std::string LtnStation::getName8() const {
     return name8;
 }

 std::string LtnStation::getName2() const {
     return name2;
 }
 
  bool LtnStation::isPressCor() const {
     return pressCor;
 }
  
  bool LtnStation::isReadable() const {
     return readable;
 }
  

 // Functions -----------------------------------------------------------------

// function extracting the information from the logfile
void LtnStation::readLogFile(std::string path){
    try{
        path = path + this->logfile;

        //convert string to char.
        const char *file = path.c_str(); 

        std::ifstream f;
        f.open(file);
        if(f.good()){
            log<INFO> ("*** opened ") % path % "successfully";
            std::string line;


            // Class that saves match results/groups. 
            boost::smatch m;

            while(std::getline(f,line)){

                //std::cout << line << std::endl;             
                if(boost::regex_search(line,m,this->wx)){
                    //std::cout << "FOUND wx" << std::endl;

                    ivg::Date time = YDHMS2date(m);

                    //double temp = boost::lexical_cast<double>(m[6]);

                    double temp = std::stod(m[6]);

                    double pressure = std::stod(m[7]);
                    double humidity  = std::stod(m[8]);

                    //std::cout << m[8] << std::endl;
                    
                    ltn::AtmosphereParam para(time,temp,pressure,humidity);
                    if(this->pressCor){
                        para.pressCor(this->gamma, this->g, this->h, this->h0);
                    }

                    this->atp.push_back(para);
                    //std::cout << para.GetTime().getDoubleMjd() << " " << para.GetTemperature() << " " << para.GetPressure() << " " << para.GetHumidity() << std::endl;
                    continue;

                }
                else if(boost::regex_search(line,m,this->cb)){
                    //std::cout << "FOUND cb" << std::endl;
                    ivg::Date time = YDHMS2date(m); 
                    double val = this->sign * std::stod(m[6]);

                    ltn::CableCor corr(time,val);
                    this->cbc.push_back(corr);

                    //for debug
                    //std::cout << corr.GetDate().getDoubleMjd() << " " << corr.GetCorr() << std::endl;
                }
            }   
        }
        else{
            log<WARNING> ("!!! opening ") % path % " failed";
            f.clear();
        }

        f.close();
    }
    catch(std::exception& e)
    {
             std::cerr << "std::exception: " << e.what() << "  during reading " << this->logfile << std::endl;
    }
}


 ivg::Date LtnStation::YDHMS2date(boost::smatch m) const{
     //std::cout << "YDHMS2date" << std::endl;
     int year      = std::stoi(m[1]);
     int doyInt    = std::stoi(m[2]);
     int hour      = std::stoi(m[3]);
     int min       = std::stoi(m[4]);
     double sec    = std::stod(m[5]);
     double fracDoy = (sec + 60*min + 3600*hour)/86400;
     double doy = doyInt + fracDoy;
     
     //std::cout << "Y:" << year << "  doy:"  << doyInt << "  h:" << hour << "  m:" << min << "  s:" << sec << "   fravDoy:>" << fracDoy << std::endl;
     
    ivg::Date date(year,doy);


    return date;
 }
 
 
void LtnStation::writeNc(std::string path, std::string sesName) const{
    //std::cout << "write Nc file" << std::endl;
    
    string fullpath = path + "/" + name8;
    if( !directory_exists(fullpath) ) // does the directory containing the vgosdb exist
    {
         log<WARNING>( "!!! directory does not exist ") % fullpath;
    }
    else
    {
        ivg::Vgosdb vgosdb(path);

        fullpath += "/TimeUTC.nc";

        log<INFO>( "try reading ") % fullpath;
        if( !is_file_exist(fullpath.c_str()) ) // is the TimeUTC file aviable ?
        { 
            log<WARNING> ("!!! file does not exist ") % fullpath;
        }
        else
        {
            //read time
            std::vector< std::vector<double> >time = vgosdb.get_vector_2d_data<double>( name8, "TimeUTC", "YMDHM" );
            std::vector<double> seconds = vgosdb.get_vector<double>( name8, "TimeUTC", "Second" );

            //create instances of class date
            std::vector<ivg::Date> dates;
            for(int i = 0; i < time.size(); ++i)
            {
                ivg::Date tmp(time[i][0], time[i][1], time[i][2], time[i][3], time[i][4], seconds[i] );
                dates.push_back(tmp );
            }

            if(cbc.size() > 0)
            {
                std::vector<double> cb;

                //loop over every time stemp (every observation)
                int s = 0; // start var for second loop. not all entries must looked up again
                for(int i = 0; i < dates.size(); i++)
                { 
                    //loop over every found cable cal correction found in the log file 
                    //and find the onethat is closest to the time tag
                    double pre = std::numeric_limits<double>::infinity();
                    bool found_min = false;
                    for(int j = s; j < cbc.size(); ++j)
                    {
                        double dif =  cbc[j].GetDate().get_double_mjd() - dates[i].get_double_mjd();
                        //std::cout << "i,j" << i << "," << j <<" pre " << pre << "   dif " << dif << std::endl;
                        if(abs(dif) > pre)
                        {
                            //std::cout << "  is minimum " << dif << std::endl;
                            cb.push_back(cbc[j-1].GetCorr()*1e-9);
                            s = j-1;
                            found_min = true;
                            break;
                        }
                        pre = abs(dif);
                    }
                    if(found_min == false) //in case the last element is the closest
                    {
                        cb.push_back(cbc.back().GetCorr()*1e-9);
                        //std::cout << "  last is minimum " << pre << std::endl;
                    }

                }

                vgosdb.create_cal_file(cb, name8, sesName);
            }
            else
            {
                log<WARNING>("!!! no cable information found for station ") % name8;
            }


            if(atp.size() > 0)
            {
                std::vector<double> T;
                std::vector<double> P;
                std::vector<double> h;

                //loop over every time stemp (every observation)
                int s = 0; // start var for second loop. not all entries must looked up again
                for(int i = 0; i < dates.size(); i++)
                { 
                    //loop over every found cable cal correction found in the log file 
                    //and find the onethat is closest to the time tag
                    double pre = std::numeric_limits<double>::infinity();
                    bool found_min = false;
                    for(int j = s; j < atp.size(); ++j)
                    {
                        double dif =  atp[j].GetTime().get_double_mjd() - dates[i].get_double_mjd();
                        //std::cout << "i,j" << i << "," << j <<" pre " << pre << "   dif " << dif << std::endl;
                        if(abs(dif) > pre)
                        {
                            //std::cout << "  is minimum " << dif << std::endl;
                            T.push_back(atp[j-1].GetTemperature());
                            P.push_back(atp[j-1].GetPressure());
                            h.push_back(atp[j-1].GetHumidity()/100);
                            s = j-1;
                            found_min = true;
                            break;
                        }
                        pre = abs(dif);
                    }
                    if(found_min == false) //in case the last element is the closest
                    {
                        T.push_back(atp.back().GetTemperature());
                        P.push_back(atp.back().GetPressure());
                        h.push_back(atp.back().GetHumidity()/100);
                        //std::cout << "  last is minimum " << pre << std::endl;
                    }

                }


                vgosdb.create_met_file(T,P,h, name8, sesName);
                }
            else
            {
                log<WARNING>("!!! no atmosphere information found for station ") % name8;
            }
        }
    }
 

}

}

