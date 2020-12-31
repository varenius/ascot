/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Session.h
 * Author: armin
 *
 * Created on 9. Januar 2016, 13:50
 */

#ifndef LTNSESSION_H
#define LTNSESSION_H

#include "ltnStation.h"
#include "logger.h"

#include <string>
#include <vector>
#include <map>

namespace ltn{

class LtnSession{

public:
    LtnSession(std::string name, std::string path, std::string outputPath,
            std::map<std::string, std::string> &nscodes,
            std::map<std::string, std::vector<std::string>> &regExp,
            std::map<std::string, std::vector<double>> &pressCor,
            std::map<std::string, int> &sign);
   // virtual ~Session();
    
private:
    std::string name;
    std::string path;
    std::string outputPath;
    
    std::vector<std::string> logfiles;
    std::vector<ltn::LtnStation> stations;
    
    
    void processAll();

    
    //findes all .log files in the given path and saves them in the vector logfiles
    void findLogFiles();
    
    

};

}

#endif /* LTNSESSION_H */

