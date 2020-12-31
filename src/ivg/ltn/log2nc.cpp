/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   log2nc.cpp
 * Author: corbin
 * 
 * Created on 11. Juni 2016, 12:11
 */

#include "log2nc.h"
#include "ltnSession.h"

namespace ltn{


    
Log2nc::Log2nc(std::string ltn_data_path) {
    this->ltn_data_path = ltn_data_path;
    readNScode();
    readRegExp();
    readPressCorFile();
    readCableSignFile();
}

void Log2nc::createNc(std::string code, std::string logpath, std::string outpath){
    ltn::LtnSession ses(code,logpath,outpath,this->nscodes,this->regExp,this->pressCor,this->sign);
}

void Log2nc::readNScode(){
   
    std::ifstream f;
    std::string path = this->ltn_data_path + "/nscodes";
    f.open(path.c_str());
    
    if(f.good()){
        log<INFO>( "*** " ) % path % " opend successfully";
        std::string line;
        
        
        while(std::getline(f,line)){
            
            std::istringstream iss(line);
            std::string ns1,ns2;
            iss >>  ns1;
            iss >>  ns2;

            this->nscodes[ns1] = ns2;
            }
    }
    else{
        log<WARNING>( "!!! ") % path  % " could not be opend ";
        f.clear();
    }
    
    f.close();
}

void Log2nc::readRegExp(){
    std::string path = this->ltn_data_path + "/regex";
    std::ifstream f;
    f.open(path.c_str());
    
    if(f.good()){
        log<INFO>( "*** " ) % path % " opend successfully";
        std::string line;
        
        
        while(std::getline(f,line)){
            if(line.find_first_of("#") != 0){
                std::istringstream iss(line);
                std::string sta,date,wx,cb;
                iss >> sta;
                iss >> date;
                iss >> wx;
                iss >> cb;
                
                //std::cout << sta << " " << date << " " << wx << " " << cb << std::endl;

                std::vector<std::string> vec;

                vec.push_back(date);
                vec.push_back(wx);
                vec.push_back(cb);


                this->regExp[sta] = vec;
            }
        }
    }
    else{
        log<WARNING>( "!!! ") % path  % " could not be opend ";
        f.clear();
    }
    
    f.close();
}

void Log2nc::readPressCorFile(){
    
    std::string path = this->ltn_data_path + "/presscor";
    std::ifstream f;
    f.open(path.c_str());
    
    if(f.good()){
        log<INFO>( "*** " ) % path % " opend successfully";
        std::string line;
        
        
        while(std::getline(f,line)){
            if(line.find_first_of("#") != 0){
                std::istringstream iss(line);
                
                std::vector<double> vec;
                
                std::string sta,s;
                iss >> sta;
                for(int i = 0; i < 4; ++i){
                    iss >> s;
                    vec.push_back(std::stod(s));
                }
                
                this->pressCor[sta] = vec;

            }
        }
    }
    else{
        log<WARNING>( "!!! ") % path  % " could not be opend ";
        f.clear();
    }
    
    f.close();
}

void Log2nc::readCableSignFile(){
    
    std::string path = this->ltn_data_path + "/sign";
    std::ifstream f;
    f.open(path.c_str());
    
    if(f.good()){
        log<INFO>( "*** " ) % path % " opend successfully";
        std::string line;
        
        
        while(std::getline(f,line)){
            
            std::istringstream iss(line);
            std::string ns, sig_str;
            int sig;
            iss >>  ns;
            iss >>  sig_str;
            if(sig_str.compare("-") == 0){
                sig = -1;
            }
            else if(sig_str.compare("+") == 0){
                sig = 1;
            }
            else{
                std::cout << sig_str << " is not a valid sign. assuming +" << std::endl;
                sig = 1; 
            }

            this->sign[ns] = sig;
            }
    }
    else{
        log<WARNING>( "!!! ") % path  % " could not be opend ";
        f.clear();
    }
    
    f.close();
}


}



