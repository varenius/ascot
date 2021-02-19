/***************************************************************************** 
 * download VMF1 files and convert them to one single GDK matrix for each    *
 * site                                                                      *
 *                                                                           *
 * 2015-02-25 - TA                                                           *
 ****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setfill, std::setw
#include <curl/curl.h>
#include <sys/stat.h>
#include "date.h"
#include "matrix.h"
#include "tictoc.h"
#include "logger.h"

#include <tclap/CmdLine.h>

loglevel g_verbose;

using namespace std;

// ...........................................................................
size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream)
// ...........................................................................
{
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

// ...........................................................................
void download_file(string http_file_path, string save_path)
// ...........................................................................
{
    CURL *curl;
    FILE *fp;
    CURLcode res;
    
    cerr << "+++ Start downloading " << http_file_path << " to " << save_path << endl;

    curl = curl_easy_init();
    if (curl)
    {
        fp = fopen(save_path.c_str(),"wb");
        curl_easy_setopt(curl, CURLOPT_URL, http_file_path.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        fclose(fp);
    }

    //if file is ".tar.bz2" untar it and delete tared file
    if (save_path.find(".tar.bz2") != std::string::npos)
    {
        string extract_folder = save_path.substr(0,save_path.find_last_of("/\\")+1)+"/bds/";
        
        // check if bds-folder already exists
        if(mkdir(extract_folder.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
            cerr << "*** bds folder " << extract_folder << " already exists " << endl;
        else
            cerr << "*** bds folder created " << extract_folder << endl;
        
        string call = "tar jxf "+save_path+" -C "+extract_folder+" ; rm "+save_path;
        system(call.c_str());
        
//        chmod_ax(extract_folder);
//        chmod_urw_grw_or(extract_folder);
        chmod_urw_grw_or(extract_folder+"*");
    }
    //if file is ".zip"unzip it and delete zip file
    else if(save_path.find(".zip") != std::string::npos)
    {
        string call = "unzip -qq -o "+save_path+" -d "+save_path.substr(0,
                      save_path.find_last_of("/\\")+1)+" ; rm "+save_path;
        system(call.c_str());
    }
    else
        chmod_urw_grw_or(save_path);
    
    cerr << "--- Finished downloading " << http_file_path << " to " << save_path << endl;
}

// ...........................................................................
void download_save_vmf1(string http_path, string save_path,
                        string output_path, int start_year=1979)
// ...........................................................................
{
    CURL *curl;
    FILE *fp;
    CURLcode res;

    ivg::Date now;
    now.now();
    
    // check if VMF-folder already exists
    if(mkdir(save_path.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) { 
        cerr << "*** VMF folder " << save_path << " already exists " << endl;
    }
    else
        cerr << "*** VMF folder created " << save_path << endl;

    map<string,ivg::Matrix> additional;

    for(int year=start_year; year<=now.get_int_year(); year++) //now.get_int_year()
    {
        ivg::Date tmp(year,1.0);
        tmp.is_leap_year();

        int final_day;
        if(tmp.is_leap_year())
            final_day = 366;
        else
            final_day = 365;

        if(year == now.get_int_year())
            final_day = now.get_int_doy()-2; // Ensure data has been processed, is released in evening
        
        string year_folder = save_path + "/" + std::to_string(year);
        
        // check if year-folder already exists
        if(mkdir(year_folder.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) { 
            cerr << "*** Folder " << year_folder << " already exists " << endl;
        }
        else
            cerr << "*** Created new folder " << year_folder << endl;

        for(int doy=1; doy<=final_day; doy++)
        {
            stringstream from, to, ending;
            ending << year << "/" << year << setfill('0') << setw(3) << right << doy <<
                   ".vmf1_r";
            from << http_path << ending.str();
            to << save_path << ending.str();

            struct stat sb;
            if (!(stat(to.str().c_str(), &sb) == 0 && S_ISREG(sb.st_mode)))
            {
                cout << "*** Downloading and saving data from " << from.str() << endl;

                curl = curl_easy_init();
                if (curl)
                {
                    fp = fopen(to.str().c_str(),"wb");
                    curl_easy_setopt(curl, CURLOPT_URL,from.str().c_str());
                    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
                    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
                    res = curl_easy_perform(curl);
                    curl_easy_cleanup(curl);
                    fclose(fp);

                    string line;
                    ifstream inStream(to.str().c_str(),ios::in);
                    if (!inStream.is_open())
                        throw runtime_error("void download_save_vmf1(): Failed to open file "+to.str());
                    else
                    {
                        while(getline(inStream,line,'\n'))
                        {
                            string ivs_name = remove_spaces_end(line.substr(0,8));

                            ivg::Matrix data( 1, 10 );
                            istringstream data_line( line.substr( 10, 87 ) );
                            data_line >> data( 0 ) >> data( 1 ) >> data( 2 ) >> data( 3 ) >> data( 4 ) >> data( 5 );
                            data_line >> data( 6 ) >> data( 7 ) >> data( 8 ) >> data( 9 );

                            string matrix_path = output_path + ivs_name + ".bin";

                            if(additional[matrix_path].size(1) == 0)
                            {
                                additional[matrix_path] = data;
                            }
                            else
                            {
                                //MJD of new dataline has to be larger than last dataline
                                if(data(0,0) > additional[matrix_path](additional[matrix_path].rows()-1,0))
                                {
                                    if(additional[matrix_path].get_col(0).find_idx(data(0,0)).size() > 0)
                                    {
                                        cout << to.str() << ": Error: MJD of new dataline already existent in dataset" << endl;
                                    }
                                    else
                                    {
                                        additional[matrix_path].append_rows(data);
                                    }
                                }
                                else
                                {
                                    cout << to.str() << ": Error: MJD of new dataline younger than existing data" << endl;
                                    cout << matrix_path << setprecision(10) << " " << data(0,0) << " < " << additional[matrix_path](additional[matrix_path].rows()-1,0) << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // check if year-folder already exists
    if(mkdir(output_path.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) { 
        cerr << "*** Folder " << output_path << " already exists " << endl;
    }
    else
        cerr << "*** Created new folder " << output_path << endl;
    
    for(map<string, ivg::Matrix>::const_iterator iter = additional.begin(),
            end = additional.end(); iter != end; ++iter )
    {
        ivg::Matrix old_data;
        
        ifstream inStream ( iter->first.c_str(), ios::in | ios::binary);
        if( inStream.is_open() )
            old_data.load_bin(iter->first);
        else
            old_data = ivg::Matrix(1,10,0.0);
        
        old_data.append_rows((iter->second));
        cout << "Save new VMF1 data to " << iter->first << endl;
        old_data.save_bin(iter->first);
    }

}

// ...........................................................................
void download_save_vmf3(string http_path, string save_path,
                        string output_path, int start_year=2008)
// ...........................................................................
{
    CURL *curl;
    FILE *fp;
    CURLcode res;

    ivg::Date now;
    now.now();
    
    // check if VMF-folder already exists
    if(mkdir(save_path.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) { 
        cerr << "*** VMF folder " << save_path << " already exists " << endl;
    }
    else
        cerr << "*** VMF folder created " << save_path << endl;

    map<string,ivg::Matrix> additional;

    for(int year=start_year; year<=now.get_int_year(); year++) //now.get_int_year()
    {
        ivg::Date tmp(year,1.0);
        tmp.is_leap_year();

        int final_day;
        if(tmp.is_leap_year())
            final_day = 366;
        else
            final_day = 365;

        if(year == now.get_int_year())
            final_day = now.get_int_doy()-1;
        
        string year_folder = save_path + "/" + std::to_string(year);
        
        // check if year-folder already exists
        if(mkdir(year_folder.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) { 
            cerr << "*** Folder " << year_folder << " already exists " << endl;
        }
        else
            cerr << "*** Created new folder " << year_folder << endl;

        for(int doy=1; doy<=final_day; doy++)
        {
            stringstream from, to, ending;
            ending << year << "/" << year << setfill('0') << setw(3) << right << doy <<
                   ".vmf3_r";
            from << http_path << ending.str();
            to << save_path << ending.str();

            struct stat sb;
            if (!(stat(to.str().c_str(), &sb) == 0 && S_ISREG(sb.st_mode)))
            {
                cout << "*** Downloading and saving data from " << from.str() << endl;

                curl = curl_easy_init();
                if (curl)
                {
                    fp = fopen(to.str().c_str(),"wb");
                    curl_easy_setopt(curl, CURLOPT_URL,from.str().c_str());
                    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
                    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
                    res = curl_easy_perform(curl);
                    curl_easy_cleanup(curl);
                    fclose(fp);

                    string line;
                    ifstream inStream(to.str().c_str(),ios::in);
                    if (!inStream.is_open())
                        throw runtime_error("void download_save_vmf3(): Failed to open file "+to.str());
                    else
                    {
                        while(getline(inStream,line,'\n'))
                        {
                            string ivs_name = remove_spaces_end(line.substr(0,8));

                            ivg::Matrix data( 1, 10 );
                            istringstream data_line( line.substr( 10, 87 ) );
                            data_line >> data( 0 ) >> data( 1 ) >> data( 2 ) >> data( 3 ) >> data( 4 ) >> data( 5 );
                            data_line  >> data( 7 ) >> data( 8 );
			    data(6)=0;data(9)=0;
                            string matrix_path = output_path + ivs_name + ".bin";

                            if(additional[matrix_path].size(1) == 0)
                            {
                                additional[matrix_path] = data;
                            }
                            else
                            {
                                //MJD of new dataline has to be larger than last dataline
                                if(data(0,0) > additional[matrix_path](additional[matrix_path].rows()-1,0))
                                {
                                    if(additional[matrix_path].get_col(0).find_idx(data(0,0)).size() > 0)
                                    {
                                        cout << to.str() << ": Error: MJD of new dataline already existent in dataset" << endl;
                                    }
                                    else
                                    {
                                        additional[matrix_path].append_rows(data);
                                    }
                                }
                                else
                                {
                                    cout << to.str() << ": Error: MJD of new dataline younger than existing data" << endl;
                                    cout << matrix_path << setprecision(10) << " " << data(0,0) << " < " << additional[matrix_path](additional[matrix_path].rows()-1,0) << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // check if year-folder already exists
    if(mkdir(output_path.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) { 
        cerr << "*** Folder " << output_path << " already exists " << endl;
    }
    else
        cerr << "*** Created new folder " << output_path << endl;
    
    for(map<string, ivg::Matrix>::const_iterator iter = additional.begin(),
            end = additional.end(); iter != end; ++iter )
    {
        ivg::Matrix old_data;
        
        ifstream inStream ( iter->first.c_str(), ios::in | ios::binary);
        if( inStream.is_open() )
            old_data.load_bin(iter->first);
        else
            old_data = ivg::Matrix(1,10,0.0);
        
        old_data.append_rows((iter->second));
        cout << "Save new VMF3 data to " << iter->first << endl;
        old_data.save_bin(iter->first);
    }

}


int main( int argc, char *argv[] )
{ 
    // .........................................................................
    // handle commandline arguments
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");

    // name of config file
    TCLAP::ValueArg<std::string> m1Arg( "d","directory","Directory path for apriori files",true,"","string");
    cmd.add( m1Arg );

    // parse commandline arguments
    cmd.parse( argc, argv );
    string directory = m1Arg.getValue();
    
    // ALWAYS CHECK IF PERMISSIONS ARE CORRECT OF STATION MATRICES! READ AND WRITE!
    //check and get new VMF files from server and update VMF station matrices files (single ascii file each day)
    //VMF1
    download_save_vmf1("https://vmf.geo.tuwien.ac.at/trop_products/VLBI/VMF1/VMF1_OP/daily/", directory+"/VMF/", directory+"/VMF/station_matrices/");
    //VMF3, two locations to get data from both before and after 2008
    download_save_vmf3("https://vmf.geo.tuwien.ac.at/trop_products/VLBI/VMF3/VMF3_EI/daily/", directory+"/VMF3/", directory+"/VMF3/station_matrices/",1980);
    download_save_vmf3("https://vmf.geo.tuwien.ac.at/trop_products/VLBI/VMF3/VMF3_OP/daily/", directory+"/VMF3/", directory+"/VMF3/station_matrices/");
    //get latest non tidal atmospheric pressure loading files (single tar-file, containing station-wise bindisp files) 
    //creates "bds" folder
    // NO LONGER USED?
    // DON'T UPDATE AUTOMATICALLY DUE TO MANUAL CHANGES AND COPIES OF BDS FILES
    // download_file("http://lacerta.gsfc.nasa.gov/aplo/aplo_bds.tar.bz2", directory+"/aplo_bds.tar.bz2");
    
    //get latest ocars file (single ascii file)
    // NO LONGER USED?
    // DON'T UPDATE OCARS AUTOMATICALLY DUE TO FORMAT CHANGES
    // download_file("http://www.gao.spb.ru/english/as/ac_vlbi/ocars.txt", directory+"/ocars.txt");
    
    //get latest antenna-info file from vlbi.geod.uni-bonn.de (single ascii file)
    // DON'T UPDATE: DONE BY CRONJOB
    // download_file("http://vlbi.geod.uni-bonn.de/Analysis/Thermal/antenna-info.txt", directory+"/antenna_info.txt");
    
    //get latest ns-codes file from CDDIS (single ascii file)
    // DON'T UPDATE: DONE BY CRONJOB
    // download_file("ftp://ivs.bkg.bund.de/pub/vlbi/ivscontrol/ns-codes.txt", directory+"/ns-codes.txt");
    
    //get latest last.erp file from GSFC (single ascii file)
    // NOT USED; USING USNOFINALS INSTEAD
    // download_file("http://gemini.gsfc.nasa.gov/500/oper/solve_save_files/last.erp", directory+"/last.erp");

    //get latest c04 eop series file from IERS (single ascii file)
    // NO LONGER USED?
    // download_file("http://datacenter.iers.org/eop/-/somos/5Rgv/latest/214", directory+"/eopc04_IAU2000.txt");
    
    //get latest c04 opa eop series file from IERS (single ascii file)
    // DON'T UPDATE: DONE BY CRONJOB
    // download_file("ftp://hpiers.obspm.fr/iers/series/opa/eopc04_IAU2000", directory+"/eopc04_IAU2000");
    
    //get latest c04 eop series file from IERS (single ascii file)
    // DON'T UPDATE: DONE BY CRONJOB
    // download_file("http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now", directory+"/eopc04_IAU2000.62-now");

    //get latest finals Bulletin A (single ascii file)
    // DON'T UPDATE: DONE BY CRONJOB
    // download_file("ftp://ftp.iers.org/products/eop/rapid/standard/finals2000A.all", directory+"/finals2000A.all");
    
    //get latest hydrology data of total earth (single zip-file, containing station-wise ascii files)
    //creates "hydlo" folder
    // WE DON'T USE HYDROLOGICAL LOADING
    //download_file("http://lacerta.gsfc.nasa.gov/hydlo/loadingfiles/cmte_series.zip", directory+"/cmte_series.zip");
    
    //get latest hydrology data of solid earth (single zip-file, containing station-wise ascii files)
    //creates "hydlo" folder
    // WE DON'T USE HYDROLOGICAL LOADING
    // download_file("http://lacerta.gsfc.nasa.gov/hydlo/loadingfiles/cmse_series.zip", directory+"/cmse_series.zip");

    return 0;
}
