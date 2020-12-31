
/* 
 * File:   db_download.cpp
 * Author: corbin
 * 
 * Created on 3. Juni 2016, 15:21
 */

#include "db_download.h"
#include "logger.h"
#include "masterfile.h"
#include "auxfunc.h"

#include "ltn/log2nc.h"

#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include <curl/curl.h>
#include <string.h>
#include <vector>

#include <sstream>  //for std::istringstream
#include <iterator> //for std::istream_iterator

#include <dirent.h> //for scandir mkdir
#include <algorithm> //for find and transform
#include <sys/stat.h> //for stat

#include <boost/algorithm/string/predicate.hpp>

#include <spawn.h>

namespace ivg
{

//Constructor/..................................................................  
Db_download::Db_download()
// .............................................................................
{
    this->_curl = curl_easy_init(); //create handle
}
// .............................................................................
Db_download::Db_download(std::string local_log_path, std::string local_vgos_path, std::string remote_log_path,
                         std::string remote_vgos_path, ivg::Masterfile * master, std::string ltn_data): Db_download()
// .............................................................................
{
    this->_local_log_path = local_log_path;
    this->_local_vgos_path = local_vgos_path;
    this->_remote_log_path = remote_log_path;
    this->_remote_vgos_path = remote_vgos_path;
        
    _masterfile = master;
    
    // initialize log2nc
    _log2nc = ltn::Log2nc(ltn_data);
    
}

//CALLBACK functions for curl
// .............................................................................
size_t Db_download::write_to_file(void *ptr, size_t size, size_t nmemb, FILE *stream)
// .............................................................................
{
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

// .............................................................................
size_t Db_download::write_to_string(void *ptr, size_t size, size_t count, void *stream)
// .............................................................................
{
  ((std::string*)stream)->append((char*)ptr, 0, size*count);
  return size*count;
}

// .............................................................................
std::vector<std::string> Db_download::list_remote_dir(std::string dir, bool * ls_ok)
// .............................................................................
{
    CURLcode res;
    if(_curl) {
        log<DETAIL> ("*** ls ") % dir;
         
        char errbuf[CURL_ERROR_SIZE];
        // provide a buffer to store errors in 
        curl_easy_setopt(_curl, CURLOPT_ERRORBUFFER, errbuf);
        // set the error buffer as empty before performing a request 
        errbuf[0] = 0;
        // ask libcurl to show us the verbose output (0,1L) 
        curl_easy_setopt(_curl, CURLOPT_VERBOSE, 0);
        
        curl_easy_setopt(_curl, CURLOPT_URL, dir.c_str()); //set URL
        curl_easy_setopt(_curl, CURLOPT_WRITEFUNCTION, &Db_download::write_to_string);
        std::string response;
        curl_easy_setopt(_curl, CURLOPT_WRITEDATA, &response);
        curl_easy_setopt(_curl, CURLOPT_DIRLISTONLY, true);
        
        res = curl_easy_perform(_curl);
        
        //write string to string vec
        std::istringstream ss(response);
        std::istream_iterator<std::string> begin(ss), end;
        std::vector<std::string> ls(begin, end); 
        
        *ls_ok  = curl_OK(res,errbuf);
        if(*ls_ok){
            log<DETAIL> ("*** remote folder content: ") % dir;
            
          

            //putting all the tokens in the vector
             if( g_verbose >= 4){
                for( std::string &str : ls ){ 
                    std::cout << str << "\t";
                }
                std::cout << std::endl;
            }
 
        }
        else
        {
            log<WARNING> ("!!! Error during listing remote directory: ") % dir;
        }
        curl_easy_reset(_curl);
        return ls;

    }
}

// .............................................................................
std::vector<std::string> Db_download::not_in_local_dir(std::string local, std::string remote)
// .............................................................................
{
    bool ls_ok;
    std::vector<std::string> rem_list = list_remote_dir(remote, &ls_ok);
    std::vector<std::string> loc_list = list_local_dir(local);
    std::vector<std::string> not_in_local;
    
    if(ls_ok)
    {
        for( std::string &str : rem_list ){ 
                if (std::find(loc_list.begin(), loc_list.end(), str) == loc_list.end()){
                    not_in_local.push_back(str);
                }
            }
    }
    return not_in_local;   
}

// .............................................................................
void Db_download::update_log(int start_year, int end_year)
// .............................................................................
{
    log<DETAIL> ("*** updating log files");
    for(int i = start_year; i <= end_year;++i){
        //ckeck weather the folder "year" exits and creats it if not
        make_local_dir(_local_log_path + std::to_string(i));
        
        // get folders that does not exist on local sytem
        std::vector<std::string> not_in_local = not_in_local_dir(_local_log_path + std::to_string(i) + "/",_remote_log_path + std::to_string(i)+ "/");
        for( std::string &str : not_in_local ){
            
            download_log(i, str);
                       
//            std::string newdir = _local_log_path + std::to_string(i) + "/" + str + "/";     
//            make_local_dir(newdir);
//            add_remote_directory_to_local(newdir,_remote_log_path + std::to_string(i)+ "/" + str + "/", ".log");
        }
    }
}

// .............................................................................
void Db_download::update_vgosDB(int start_year, int end_year, bool force_load)
// .............................................................................
{
    log<DETAIL>("*** updating vgosDB");
    for(int i = start_year; i <= end_year;++i){
        std::string dir = _local_vgos_path + std::to_string(i) + "/";
        //ckeck weather the folder "year" exits and creats it if not
        make_local_dir(dir);
        
        std::vector<std::string> sessions;
        if(force_load)
        {
            bool ls_ok;
            sessions = list_remote_dir(_remote_vgos_path + std::to_string(i)+ "/", &ls_ok);
            if(!ls_ok)
            {
                log<WARNING>("!!! Error during listing remote content ");
                break;
            }
        }
        else
        {
            // get folders that does not exist on local sytem
            sessions = not_in_local_dir(dir ,_remote_vgos_path + std::to_string(i)+ "/");
        }
            
        for( std::string &str : sessions )
        {
            
            std::size_t ending = str.find(".tar.gz");
            str.erase(ending, str.length());
            std::cout << str << std::endl;
            
            //download_vgosDB(str,false);
            add_remote_archiv_to_local(dir,_remote_vgos_path + std::to_string(i) + "/", str + ".tar.gz");
        }
    }
    
}

// .............................................................................
void Db_download::download_vgosDB(std::string db, bool force_load, bool load_log_too)
// .............................................................................
{
        
    // Load vgosDB
    int year = (*_masterfile).get_session_info(db).date.get_int_year();
    log<INFO> ("*** Downloading vgosDB ") % db;
    std::string dir = _local_vgos_path + std::to_string(year) + "/";
    make_local_dir(dir);
    std::string dir2 = dir + db + "/";
    if(directory_exists(dir2) == false || force_load == true) // database is already in local dir
    {
        //make_local_dir(dir);
        //if(add_remote_to_local_rec(dir,_remote_vgos_path + std::to_string(year) + "/" + db + "/"))
        if(add_remote_archiv_to_local(dir,_remote_vgos_path + std::to_string(year) + "/", db + ".tar.gz"))
        {
            _new_dbs.push_back(db);
            
            if(load_log_too)
            {
            // Load corrsponding log file
            string code = (*_masterfile).get_session_info(db).code;
            std::transform(code.begin(), code.end(), code.begin(), ::tolower); //capital letter to lower letter
            download_log(year,code);

            // creat Met and Cal-Cable ncfiles
            _log2nc.createNc(code,_local_log_path + "/" + std::to_string(year) + "/" + code + "/",dir2);
            }
        }
        else
        {
            remove( (dir +  db + ".tar.gz").c_str());
            log<WARNING> ("!!! ") % db % " could not be downloaded. Maybe it is not in the remote directory";
        }
    }
    else
    {
        log<INFO> ("*** ") % dir % "already in local dirctory. Is not loaded again";
    }
}


// .............................................................................
void Db_download::download_log(int year ,std::string code)
// .............................................................................
{
    log<INFO> ("*** downloading logfiles for session ") % code;
    std::string dir = _local_log_path + std::to_string(year) + "/";
    make_local_dir(dir);
    dir += code + "/";
    make_local_dir(dir);
    if (! add_remote_directory_to_local(dir,_remote_log_path + std::to_string(year) + "/" + code + "/",".log"))
    {
        log<WARNING> ("!!! Error at session ") % code ;
    }

}

// .............................................................................
void Db_download::make_local_dir(std::string dir)
// .............................................................................
{
    // does dir already exists ?
    if (directory_exists(dir)){
        log<DETAIL> ("*** ") % dir % " already exists";
    }
    else{
        log<DETAIL> ("*** creating new directory ") % dir;
        mkdir(dir.c_str() , 0700);
    }
    
}

// .............................................................................
bool Db_download::add_remote_directory_to_local(std::string local, std::string remote, std::string filter)
// .............................................................................
{
    bool load_ok = true;
    bool ls_ok;
    std::vector<std::string> rem = list_remote_dir(remote, &ls_ok);
    
    //only perform if listing remote directory was successfull
    if(ls_ok){
        CURLcode res;
        if(_curl) {
            log<DETAIL> ("*** copy *") % filter % " from " % remote % " to " % local;

            char errbuf[CURL_ERROR_SIZE];
            // provide a buffer to store errors in 
            curl_easy_setopt(_curl, CURLOPT_ERRORBUFFER, errbuf);

            // ask libcurl to show us the verbose output (0,1L) 
            curl_easy_setopt(_curl, CURLOPT_VERBOSE, 0);

            curl_easy_setopt(_curl, CURLOPT_WRITEFUNCTION, &Db_download::write_to_file);
            curl_easy_setopt(_curl, CURLOPT_DIRLISTONLY, false);

            // write to file
            for( std::string &str : rem ){

                // set the error buffer as empty before performing a request 
                errbuf[0] = 0;

                // only loade files that ends withe the string filter (e.g. to only load '.log' files)
                if( ends_with(str,".log") ){

                    //craete file
                    FILE *fp;
                    std::string file = local+str;
                    fp = fopen(file.c_str(),"wb");
                        std::string url = remote + str;
                        curl_easy_setopt(_curl, CURLOPT_URL, url.c_str()); //set URL
                        curl_easy_setopt(_curl, CURLOPT_WRITEDATA, fp);
                        res = curl_easy_perform(_curl); // load file

                        if(curl_OK(res,errbuf)){
                            log<DETAIL> ("*** ") %  str % " copied successfullly"; 
                        }
                        else
                        {
                            load_ok = false;
                        }
                    fclose(fp);
                }
            }
            curl_easy_reset(_curl);
        }
    }
    else
    {
        load_ok = false;
    }
    
    return load_ok;
 }

// .............................................................................
bool Db_download::add_remote_archiv_to_local(std::string local, std::string remote, std::string archiv_name)
// .............................................................................
{
    bool load_ok = true;
    
    CURLcode res;
        if(_curl) {
            log<DETAIL> ("*** copy ") % archiv_name %  " from " % remote % " to " % local;

            char errbuf[CURL_ERROR_SIZE];
            // provide a buffer to store errors in 
            curl_easy_setopt(_curl, CURLOPT_ERRORBUFFER, errbuf);

            // ask libcurl to show us the verbose output (0,1L) 
            curl_easy_setopt(_curl, CURLOPT_VERBOSE, 0);

            curl_easy_setopt(_curl, CURLOPT_WRITEFUNCTION, &Db_download::write_to_file);
            curl_easy_setopt(_curl, CURLOPT_DIRLISTONLY, false);
            
            // set the error buffer as empty before performing a request 
            errbuf[0] = 0;
            
            //craete file
            FILE *fp;
            std::string file = local + archiv_name;
            fp = fopen(file.c_str(),"wb");
            std::string url = remote + archiv_name;
            curl_easy_setopt(_curl, CURLOPT_URL, url.c_str()); //set URL
            curl_easy_setopt(_curl, CURLOPT_WRITEDATA, fp);
            res = curl_easy_perform(_curl); // load file

            if(curl_OK(res,errbuf)){
                log<DETAIL> ("*** ") %  archiv_name % " copied successfullly";
                chmod_urw_grw_or(file.c_str());
                
                
                // With posix spawn
                // Does not wait until extraction is ready before deleting
                
                    pid_t processID;
                    char **environ;
                    int status;
                    
                    char *argV[] = {const_cast<char*>("tar"), const_cast<char*>("-xzf"), const_cast<char*>(file.c_str()),
                                    const_cast<char*>("-C"), const_cast<char*>(local.c_str()), NULL};
                    status = posix_spawn(&processID,"/bin/tar",NULL,NULL,argV,environ);
                    
                    if( status != 0)
                    {
                        load_ok = false;
                        log<WARNING> ("!!! ") %  archiv_name % " could not be extracted with tar";
                    } 
            }
            else
            {
                load_ok = false;
            }

            fclose(fp);
            
        }
        else
        {
             load_ok = false;
        }

    return load_ok;
}
// .............................................................................
bool Db_download::add_remote_to_local_rec(std::string local, std::string remote)
// .............................................................................
{
    bool ls_ok;
    bool load_ok = true;
    std::vector<std::string> rem = list_remote_dir(remote, &ls_ok); // list remote directory
    
    //only perform if listing remote directory was successfull
    if(ls_ok){
        char errbuf[CURL_ERROR_SIZE];
        // provide a buffer to store errors in 
        curl_easy_setopt(_curl, CURLOPT_ERRORBUFFER, errbuf);
        // ask libcurl to show us the verbose output (0,1L) 
        curl_easy_setopt(_curl, CURLOPT_VERBOSE, 0);

        curl_easy_setopt(_curl, CURLOPT_WRITEFUNCTION, &Db_download::write_to_file);
        curl_easy_setopt(_curl, CURLOPT_DIRLISTONLY, false);

        //loop over every element in rem. If its a file load it if its a directory call this function recursiv
        for( std::string &str : rem ){
            //is a file
            if(str.find(".") != std::string::npos && str.compare(".") != 0 && str.compare("..") != 0){           
                CURLcode res;
                if(_curl) {
                    log<DETAIL> ("*** copy ") % str % " from " % remote % " to " % local;

                    // set the error buffer as empty before performing a request 
                    errbuf[0] = 0;

                    FILE *fp;
                    std::string file = local+str;
                    fp = fopen(file.c_str(),"wb");
                        std::string url = remote + str;
                        curl_easy_setopt(_curl, CURLOPT_URL, url.c_str()); //set URL
                        curl_easy_setopt(_curl, CURLOPT_WRITEDATA, fp);
                        res = curl_easy_perform(_curl); // load file

                        if(curl_OK(res,errbuf)){
                            log<DETAIL> ("*** copy success");
                        }
                        else{
                            load_ok = false;
                        }
                    fclose(fp);
                }
            }
            //is a dir
            else if(str.compare(".") != 0 && str.compare("..") != 0){
                std::string newdir = local + str + "/";  
                make_local_dir(newdir);
                
                // recursive function call
                if(! add_remote_to_local_rec(newdir, remote + str + "/") ){
                    load_ok = false; //never change the bool to true. Its just true at the initialisation
                }
            }
        }
        curl_easy_reset(_curl);
    }
    else
    {
        load_ok = false;
    }
    
    return load_ok;
}
// .............................................................................
bool Db_download::curl_OK(CURLcode &res, char errbuf[CURL_ERROR_SIZE])
// .............................................................................
{
    // if the request did not complete correctly, show the error
    // information. if no detailed error information was written to errbuf
    // show the more generic information from curl_easy_strerror instead.
    if(res != CURLE_OK && g_verbose >= 1) {
        size_t len = strlen(errbuf);
        fprintf(stderr, "\nlibcurl: (%d) ", res);
        if(len)
            fprintf(stderr, "%s%s", errbuf,
                ((errbuf[len - 1] != '\n') ? "\n" : ""));
        else{
            fprintf(stderr, "%s\n", curl_easy_strerror(res));
        }   
        return false;
    }
    else{
        return true;
    }
}
// .............................................................................
void Db_download::clean_up(){
// .............................................................................
log<DETAIL> ("*** cleaning up vgos dir. Removing archives");

// get list of year folders
std::vector<std::string> folders = list_local_dir(_local_vgos_path);
    // loop over all folders
    for( std::string &year : folders ){
        if(year.compare(".") != 0 && year.compare("..") != 0){
            std::string dir = _local_vgos_path + "/" + year + "/";
            
            //get folder content and loop over the files
            std::vector<std::string> files = list_local_dir(dir);
            for( std::string &file : files ){
                //delte file if it is ending with .tar.gz
                if( boost::algorithm::ends_with(file, ".tar.gz") ){
                    log<DETAIL> ("*** Deleting ") % dir % file;
                    std::remove((dir + file).c_str());
                }
            }
            
        }
     }
}


// .............................................................................
std::vector<std::string> Db_download::getNew_dbs() const {
// .............................................................................
    return _new_dbs;
}
    
}
