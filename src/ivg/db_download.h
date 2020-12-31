
/* 
 * File:   db_download.h
 * Author: corbin
 *
 * Created on 3. Juni 2016, 15:21
 */

#ifndef DB_DOWNLOAD_H
#define DB_DOWNLOAD_H

#include <curl/curl.h>
#include <string>
#include <vector>

#include "masterfile.h"
#include "ltn/log2nc.h"


/**
*
* @brief Downloading vgosDB and log files from ftp server
* @author AC
* @date 2016-06-09
*/

namespace ivg
{

class Db_download {
public:
    //constructor---------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    /**
    *  \b Description: \n
    *        Default constructor
    *  \param [in] no input parameters needed
    *  \return An instance of the class 'Db_download'
    */
    Db_download();
    
    /**
    *  \b Description: \n
    *        using six parameters: remote and local path to log files and vgosDB directory and masterfile
    *  \param [in] [std::string] local_log_path
    *              [std::string] local_vgos_path
    *              [std::string] remote_log_path
    *              [std::string] remote_vgos_path
    *              [ivg::Masterfile] masterfile
    *              [std::string] path to log2nc config files 
    *  \return An instance of the class 'Db_download'
    */
    Db_download(std::string local_log_path, std::string local_vgos_path, std::string remote_log_path, std::string remote_vgos_path, ivg::Masterfile * master, std::string ltn_data);
    
    //destructor----------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    /**
    *  \b Description: \n
    *        Default destructor
    *  \param [in] no input parameters needed
    *  \return nothing
    */
    ~Db_download() {
        curl_easy_cleanup(_curl);
        clean_up();
    }
    
    //public functions----------------------------------------------------------
    //--------------------------------------------------------------------------
    
    /**
    *  \b Description: \n
    *        Function saving remote(ftp) filelist in a vector
    *  \param [in] [std::string] path e.g. "ftp://ivs.bkg.bund.de/"
    *              [bool pointer] true if there occured no error
    *  \return [vector<std::string>] all entries in remote location
    */
    std::vector<std::string> list_remote_dir(std::string dir, bool * ls_ok);
    
   /**
    *  \b Description: \n
    *        Function downloading all logfiles between the start and end year.
    *        If there is allready an folder with the session name it will 
    *        not be updated 
    *  \param [in] [int] start_year
    *              [int] end_year
    *  \return nothing
    */    
    void update_log(int start_year, int end_year);
    
    
    /**
    *  \b Description: \n
    *       Function downloading  vgos databases in the years between start and
    *       end year. 
    *        This function does not use the masterfile entries!
    *  \param [in] [int] start_year
    *              [int] end_year
    *              [bool] force_load
    *  \return nothing
    */    
    void update_vgosDB(int start_year, int end_year, bool force_load = false);
    
    /**
    *  \b Description: \n
    *        Function downloading a  certain vgos database. The corresponding
    *        log file is also loaded if load_log_too is true.
    *        This function uses the masterfile entries. 
    *  \param [in] [string] db
    *              [bool] force_load
    *              [bool] load_log_too
    *  \return nothing
    */    
    void download_vgosDB(std::string db, bool force_load = false,  bool load_log_too = true);
    
     /**
    *  \b Description: \n
    *        Getter for vector new dbs
    *  \param [in] nothing
    *  \return [vector<string>] all loaded dbs
    */    
    std::vector<std::string> getNew_dbs() const;
    
    // removes archives
    void clean_up();
    
    
    
private:
    //private attributes---------------------------------------------------------
    //--------------------------------------------------------------------------
    std::string _local_log_path;
    std::string _local_vgos_path;
    std::string _remote_log_path;
    std::string _remote_vgos_path;
    
    CURL *_curl; //curl handele
        
    ivg::Masterfile * _masterfile;
    ltn::Log2nc _log2nc;
    
    std::vector<std::string> _new_dbs;
    
    //private functions---------------------------------------------------------
    //--------------------------------------------------------------------------
    
    /**
    *  \b Description: \n
    *        Function comparing a remote and a local folder. All files that are
    *        in the remote folder but not in the local folder are returned
    *  \param [in] [std::string] path to local folder
    *  \param      [std::string] path to remote folder
    *  \return [vector<std::string>] all entries that are in remote location but not in the local
    */
    std::vector<std::string> not_in_local_dir(std::string local, std::string remote);
    
    /**
    *  \b Description: \n
    *        Function that creates a new folder in a local directory (if it does
    *        not exist). Note that all folders except the last one have to exist.    
    *  \param [in] [std::string] dirrectory
    *  \return nothing
    */
    void make_local_dir(std::string dir);
    
    /**
    *  \b Description: \n
    *        Function downloading log files beloning to a certain database
    *  \param [in] [int] year
    *              [string] db
    *  \return nothing
    */    
    void download_log(int year ,std::string code);
    
    
    /**
    *  \b Description: \n
    *        copys all files in a remote directory to a local directory. It will
    *        not download recursive folders only files containing the string 
    *        filter will be downloaded. Make sure there are no folders in the
    *        remote directory
    *  \param [in] [std::string] path to local directory
    *              [std::string] path to remote directory
    *              [std::string] only files ending with this string are loaded
    *  \return [bool] true if no error occured
    */
    bool add_remote_directory_to_local(std::string local, std::string remote,std::string filter);
    
    /**
    *  \b Description: \n
    *        copys an entire folder with all subfolders and all files from a 
    *        remote directory to a local directory. Missing folders in the local
    *        directory are created automatically. Note that all entries 
    *        containing a '.' are treted as file. everything else is treated as 
    *        a folder.
    *  \param [in] [std::string] path to local directory
    *              [std::string] path to remote directory
    *  \return [bool] true if no error occured
    */
    bool add_remote_to_local_rec(std::string local, std::string remote);
    
    
    /**
    *  \b Description: \n
    *        copys one file and extracts it (only tar.gz)
    *  \param [in] [std::string] path to local directory
    *              [std::string] path to remote directory
    *              [std::string] name of archive 
    *  \return [bool] true if no error occured
    */
    bool add_remote_archiv_to_local(std::string local, std::string remote, std::string archiv_name);
    
    
    // Function to print Curl error message if necessary
    bool curl_OK(CURLcode &res, char errbuf[CURL_ERROR_SIZE]);

    //CALLBACK functions for curl
    static size_t write_to_file(void *ptr, size_t size, size_t nmemb, FILE *stream);
    static size_t write_to_string(void *ptr, size_t size, size_t count, void *stream);
    
};

}

#endif /* DB_DOWNLOAD_H */

