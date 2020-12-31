/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Session.cpp
 * Author: armin
 * 
 * Created on 9. Januar 2016, 13:50
 */


#include "ltnSession.h"

#include "ltnStation.h"
#include "auxfunc.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

namespace ltn{
   
//Construktor
LtnSession::LtnSession(std::string name, std::string path, std::string outputPath,
                 std::map<std::string, std::string> &nscodes,
                 std::map<std::string, std::vector<std::string>> &regExp,
                 std::map<std::string, std::vector<double>> &pressCor,
                 std::map<std::string, int> &sign){
        this->name = name;
        this->path = path;
        this->outputPath = outputPath;
        findLogFiles();

        
        // Initialisierung der Stationen zu denen log Files gefunden wurden
        // Die entsprechen ns codes und Regular Expressions zum einlesen werden
        // vom Konstruktor der Klasse Station zugewiesen
        for(int i=0; i < this->logfiles.size(); i++){
            ltn::LtnStation temp(this->logfiles[i],nscodes,regExp,pressCor,sign);
            this->stations.push_back(temp);
         }
        
        processAll();  
    }

void LtnSession::processAll(){    

         for(int i=0; i < this->stations.size(); i++){
            // Einlesen der LogFiles und schreiben der .wx .cb Dateien
            if(this->stations[i].isReadable()){
                this->stations[i].readLogFile(this->path);
                this->stations[i].writeNc(this->outputPath,this->name);
            }
         }
         
    
}


void LtnSession::findLogFiles(){
    struct stat sb;
    
    const char * c = this->path.c_str();
    //ckeck wether dir is an existing dir
    if (stat(c, &sb) == 0 && S_ISDIR(sb.st_mode)){
        struct dirent **namelist;
        int i,n;

            n = scandir(c, &namelist, 0, alphasort); //get all entries in the directory

            // check wether there are any entries
            if (n < 0)
                log<WARNING>( "!!! scandir: Empty Directory" );
            else {
                // save all entries including ".log" and the session name
                if( g_verbose >= 4){
                    std::cout << "the following log files (.log) were found in " 
                              << this->path << " for the Session " << this->name 
                              << " :" <<std::endl;
                }
                for (i = 0; i < n; i++) {
                    std::string s = namelist[i]->d_name;
                    if( ends_with(s,".log") && s.find(this->name) != std::string::npos){
                            this->logfiles.push_back(s);
                            log<INFO>( s );
                        }
                    free(namelist[i]);
                    }
                }
           free(namelist);
    }
     else{
        log<WARNING>("!!! ") % this->path % " is not a valid path";
    }
}



          
   
}

