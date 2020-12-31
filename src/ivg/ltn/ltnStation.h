/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Station.h
 * Author: armin
 *
 * Created on 9. Januar 2016, 13:51
 */

#ifndef LTNSTATION_H
#define LTNSTATION_H

#include "atmosphereParam.h"
#include "cableCor.h"
#include "date.h"

#include <iostream>
//#include <regex>
#include <boost/regex.hpp>
#include <fstream>
#include <vector>
#include <string>
#include <map>

namespace ltn{

class LtnStation{
public:
    
    /*
     *  logfile : name of the logfile (without path)
     *  nscodes : Map key: 2-letter code    value: 8-letter code
     *  regExp  : Map key: 2-letter code    value: Vector with date wx cb
     * pressCor : Map key: 2-letter code    value: Vector with gamma g h h h0
     */
    LtnStation(std::string logfile, std::map<std::string, std::string> &nscodes,
            std::map<std::string, std::vector<std::string>> &regExp,
            std::map<std::string, std::vector<double>> &pressCor,
            std::map<std::string, int> &sign);
    
    /*
     *  path : Path of the folder containing the logfiles (e.g. /data/logs/s15234/ )
     */
    void readLogFile(std::string path);
    
    /*   Creates .nc files used by ascot
     *   path : Output Path without filename
     *   sesName : Name of the Session (Attribute in Session)
     *   hasCB cal by Reference : true 
     */
    void writeNc(std::string path, std::string sesName) const;
    
    
    
    
    // Getter Setter
    std::string getName8() const;
    std::string getName2() const;
    bool isPressCor() const;
    bool isReadable() const;

private:
    std::string logfile; //name des Logfiles (ohne Pfad))
    std::string name2; //2 stellige nscode
    std::string name8; //8 stellige nscode
    
    boost::regex wx; //regular expression for atmosphere
    boost::regex cb; //regular expression for cable cal
    
    std::vector<ltn::AtmosphereParam> atp; //vector containing all wx entries
    std::vector<ltn::CableCor> cbc; //vector containing all cb entries
    
    
    //Parameter for Pressure Correction
    bool readable;
    bool pressCor;
    double gamma; //average temperature lapse rate
    double g; //gravity
    double h; //Height of the barometer
    double h0; //h0: height of the VLBI Reference point
    
    int sign;
    
    // help functions
    ivg::Date YDHMS2date(boost::smatch m) const;
    
};

}

#endif /* LTNSTATION_H */

