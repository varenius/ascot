/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   log2nc.h
 * Author: corbin
 *
 * Created on 11. Juni 2016, 12:11
 */

#ifndef LOG2NC_H
#define LOG2NC_H

#include <string>
#include <vector>
#include <map>

namespace ltn{

class Log2nc {
public:
    
    Log2nc(){};
    
    Log2nc(std::string ltn_data_path);
    
    void createNc(std::string code, std::string logpath, std::string outpath);
    

private:
    // Maps containing data read from ./Data/* files
    std::map<std::string, std::string> nscodes; 
    std::map<std::string, std::vector<std::string>> regExp;
    std::map<std::string, std::vector<double>> pressCor;
    std::map<std::string, int> sign; // cableCal sign 
    
    std::string ltn_data_path;
    
     //read File containing 2-letter and 8-letter codes
    void readNScode();
    
    //read File containing 2-letter code and Regular Expressions
    void readRegExp();
    
    //read file containing Pressure correction parameter
    void readPressCorFile();
    
    //read file containing Cable sign parameter
    void readCableSignFile();
    

};

}

#endif /* LOG2NC_H */

