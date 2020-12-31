#ifndef VGOSDB_H
#define	VGOSDB_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "ivg_const.h"
#include "parser.h"
#include "matrix.h"
#include "auxfunc.h"
#include "logger.h"
#include "netcdfcpp.h"
#include "dirent.h"
#include "date.h"
#include "wrapper.h"

//TODO: The new netcdf version is not compatible with the here used version. But the new version is in the repositories --> upgrade code to new version

using namespace std;

namespace ivg
{
    
    
class Vgosdb {
    
#define VGOSDB_HEAD  "Head.nc"
#define VGOSDB_OBS   "Observables"
    
public:
    
    Vgosdb();
    
    Vgosdb(const string directory);
    Vgosdb(const string directory, ivg::Wrapper* wrapper_ptr);
       
    vector<string> get_string( string folder, string nc, string var );
    
    ivg::Matrix get_matrix( string folder, string nc, string var );
        
    bool does_file_exist( string folder, string nc );
    
    bool does_variable_exist( string folder, string nc, string var);
    
    void copy_file( string from_folder, string from_nc, string to_folder, string to_nc);
    
    NcFile* create_file( string folder, string nc, NcFile::FileMode mode = NcFile::Replace );
    
    
    /**
    *  \b Description: \n
    *        Function creating a Cal-Cable.nc file in the vgos Database (in _directory)
    *  \param [in] [vector<double>] cable correction values (s)
    *              [string] 8-letter Station name
    *              [string] session name
    *  \return nothing
    */
    void create_cal_file( vector<double> &cbc, string station, string sesion );
    
    
     /**
    *  \b Description: \n
    *        Function creating a Met.nc file in the vgos Database (in _directory)
    *  \param [in] [vector<double>] Temperatur values (C)
    *              [vector<double>] Preassure values (hPa)
    *              [vector<double>] relative humidity (%) 
    *              [string] 8-letter Station name
    *              [string] session name
    *  \return nothing
    */
    void create_met_file( vector<double> &T,vector<double> &P,vector<double> &h, string station, string session);
    
    void create_NumGroupAmbig_file( vector<short> &ambig, string session, ivg::band type, string editing = "", bool phase = false);
    
    void create_GroupDelayFull_file( vector<double> &gd_full, string session, ivg::band type, string editing = "", bool phase = false );

    void create_EffFreq_file( std::vector<double> &eff_freq_iono_g, std::vector<double> &eff_freq_iono_ph,std::vector<double> &eff_freq_iono_rt,
			      string session, ivg::band type,  std::string editing  );

    void create_Edit_file( std::vector<int> &flag_g, std::vector<int> &flag_ph,std::vector<int> &flag_rt,
			   string session,  std::string editing  );
  
    void create_IonoGroup_file( ivg::Matrix &delta_tau_x, ivg::Matrix &delta_tau_x_sigma, std::vector<short> error_flag,
                                string session, ivg::band type = ivg::band::X,  std::string editing = "", bool phase = false );
    
    void create_ClockBreak_file( std::vector<string> stations, std::vector<double> epochs_mjd, string session );
    
    template <typename T> T get_scalar(string folder, string nc, string var, int i=-1, int j=-1) {
   
         // create path to file which need to be open
        string path = _directory + "/" + folder + "/" + nc + ".nc";   
        
        log<DETAIL>("*** read ") % var % " in " % path;
                
        // open the file for read access
        NcFile file(path.c_str(), NcFile::ReadOnly);
        
        if(!file.is_valid())
            throw runtime_error( "template <typename T> T get_scalar(string, string, string, int i=-1, int j=-1): Could not open "+path);
        
        // check if desired variable is existent in opened file
        NcVar *tmp_var = NULL;
        for(int l=0; l<file.num_vars(); l++)
        {
            if(string(file.get_var(l)->name()) == var)
            {
                tmp_var = file.get_var(var.c_str());
                break;
            }
        }
        
        if(tmp_var == NULL)
            throw runtime_error( "template <typename T> T get_scalar(string, string, string, int i=-1, int j=-1): Variable "+var+" not existent in "+path);
       
        
        if(i == -1 && j == -1)
        {
            if(tmp_var->num_dims() == 1 && tmp_var->get_dim(0)->size() == 1)
            {
                T value;
                tmp_var->get(&value,1,0);
                return value;
            }
        }
        else if(i >= 0 && j == -1)
        {
            if(tmp_var->num_dims() == 1 && tmp_var->get_dim(0)->size() >= 1)
            {
                int n = tmp_var->get_dim(0)->size(); // e.g. 2131 delays

                T values[n];                
                tmp_var->get(&values[0],n);               
                
                return values[i];
            }    
        }
        else if(i >= 0 && j >=0)
        {
            if(tmp_var->num_dims() == 2 && tmp_var->get_dim(0)->size() >= 1 && tmp_var->get_dim(1)->size() >= 1)
            {
                
                int n1 = tmp_var->get_dim(0)->size(); // e.g. 7 stations
                int n2 = tmp_var->get_dim(1)->size(); // e.g. 8 chars for station-name

                T values[n1][n2];                
                tmp_var->get(&values[0][0],n1,n2);

                return values[i][j];
            }
        }  
    }  
    // ...........................................................................
    template <typename T> vector<T> get_vector( string folder, string nc, string var ) {
        
         // create path to file which need to be open
        string path = _directory + "/" + folder + "/" + nc + ".nc"; 
        
        log<DETAIL>("*** read ") % var % " in " % path;

        // open the file for read access
        NcFile file(path.c_str(), NcFile::ReadOnly);
        
        if(!file.is_valid())
            throw runtime_error( "template <typename T> vector<T> get_vector( string, string, string ): Could not open "+path);
        
        // check if desired variable is existent in opened file
        NcVar *tmp_var = NULL;
        for(int i=0; i<file.num_vars(); i++)
        {
            if(string(file.get_var(i)->name()) == var)
            {
                tmp_var = file.get_var(var.c_str());
                break;
            }
        }
        
        if(tmp_var == NULL)
            throw runtime_error( "template <typename T> vector<T> get_vector( string folder, string nc, string var ): Variable "+var+" not existent in "+path);
                
        
        if(tmp_var->num_dims() == 1 && tmp_var->get_dim(0)->size() > 1)
        {
            int n = tmp_var->get_dim(0)->size(); // e.g. 2131 delays

            T values[n];                
            tmp_var->get(&values[0],n);

            vector<T> output (values, values + sizeof(values) / sizeof(values[0]) );
            
            return output;
        }
        // if attribute REPEAT is existent, the vector will contain n-times the same value
        else if(tmp_var->num_dims() == 1  && tmp_var->get_dim(0)->size() == 1)
        {   
            int repeats = 1;
            // if attribute REPEAT is not existent, repeats keeps = 1
            for(int i=0; i<tmp_var->num_atts(); i++)
                if(string(tmp_var->get_att(i)->name()) == "REPEAT")
                    repeats = tmp_var->get_att(i)->as_int(0);

            T values_tmp[1];   
            tmp_var->get(&values_tmp[0],1);
            
            vector<T> out(repeats,values_tmp[0]);
            
            return out;
        }
    }
    
    
    // ...........................................................................
    template <typename T> vector<vector<T>> get_vector_2d_data( string folder, string nc, string var ) {
        
        using vecvec = std::vector< std::vector<T> >;
            
         // create path to file which need to be open
        string path = _directory + "/" + folder + "/" + nc + ".nc";  
        
        log<DETAIL>("*** read ") % var % " in " % path;

        // open the file for read access
        NcFile file(path.c_str(), NcFile::ReadOnly);
        
        if(!file.is_valid())
            throw runtime_error( "template <typename T> vector<T> get_vector_2d_data( string, string, string ): Could not open "+path);
        
        // check if desired variable is existent in opened file
        NcVar *tmp_var = NULL;
        for(int i=0; i<file.num_vars(); i++)
        {
            if(string(file.get_var(i)->name()) == var)
            {
                tmp_var = file.get_var(var.c_str());
                break;
            }
        }
        
        if(tmp_var == NULL)
            throw runtime_error( "template <typename T> vector<T> get_vector_2d_data( string folder, string nc, string var ): Variable "+var+" not existent in "+path);
                
        

        

            int n = tmp_var->get_dim(0)->size(); // e.g. 2131 delays
            int m = tmp_var->get_dim(1)->size();
            
            T values[n][m];  
            
            //Retrieve the variable var
            tmp_var->get(&values[0][0],n,m);

            vecvec output;
            for(int i = 0; i < n; ++i)
            {
                vector<T> col;
                for(int j = 0; j < m; ++j){
                    col.push_back(values[i][j]);
                }
                output.push_back(col);
            }
            

            return output;


    }
    //
   template <typename T> vector<vector<vector<T> > > get_vector_3d_data( string folder, string nc, string var ) {
        
     using vecvecvec = std::vector<std::vector< std::vector<T> > >;
        using vecvec = std::vector< std::vector<T> >;    
         // create path to file which need to be open
        string path = _directory + "/" + folder + "/" + nc + ".nc";  
        
        log<DETAIL>("*** read ") % var % " in " % path;

        // open the file for read access
        NcFile file(path.c_str(), NcFile::ReadOnly);
        
        if(!file.is_valid())
            throw runtime_error( "template  <typename T> <typename T> vector<T> get_vector_3d_data( string, string, string ): Could not open "+path);
        
        // check if desired variable is existent in opened file
        NcVar *tmp_var = NULL;
        for(int i=0; i<file.num_vars(); i++)
        {
            if(string(file.get_var(i)->name()) == var)
            {
                tmp_var = file.get_var(var.c_str());
                break;
            }
        }
        
        if(tmp_var == NULL)
            throw runtime_error( "template  <typename T> <typename T> vector<T> get_vector_3d_data( string folder, string nc, string var ): Variable "+var+" not existent in "+path);
                
        

        

            int n = tmp_var->get_dim(0)->size(); // e.g. 2131 delays
            int m = tmp_var->get_dim(1)->size();
            int p = tmp_var->get_dim(2)->size();
            T values[n][m][p];  
            
            //Retrieve the variable var
            tmp_var->get(&values[0][0][0],n,m,p);

            vecvecvec output;
            for(int i = 0; i < n; ++i)
            {
	        vecvec col;
                for(int j = 0; j < m; ++j){
		  vector<T> col2;
		  for(int k = 0; k < p; ++k){
                    col2.push_back(values[i][j][k]);
		  }
		  col.push_back(col2);
                }
                output.push_back(col);
            }
            

            return output;


    }
    // ...........................................................................
    template <typename T> void set_vector( vector<T> vec, string folder, string nc, string var) {
        
         // create path to file which need to be open
        string path = _directory + "/" + folder + "/" + nc + ".nc";     

        // open the file for write access
        NcFile file(path.c_str(), NcFile::Write);
        
        if(!file.is_valid())
            throw runtime_error( "template <typename T> void set_vector( vector<T>, string , string, string ): Could not open "+path);
        
        NcVar *data = file.get_var(var.c_str());
        
        // in case of REPEAT in DelayFlag, create new vector
        if(data->num_dims() == 1  && data->get_dim(0)->size() == 1)
        {   
            data->rename("RepeatDelayFlag");
            NcDim *repeatDim = file.add_dim("REPEATS", data->get_att("REPEAT")->as_int(0));
            NcVar *newDelayFlag = file.add_var("DelayFlag", ncShort, repeatDim);
            data = newDelayFlag;
            log<WARNING>("!!! Replacing REPEAT variable for") % var;
        }
        
        data->put(&vec[0], vec.size());
        
    }
    
    
    
    // ...........................................................................
private:
    
    
    string _create_file_name(string nc, string band, bool extension =true , string editing = "", string version = "");
    
    /**
    *  \b Description: \n
    *        Function adding the attributes "LCODE","CreateTime","Definition",
    *        "Units"  to a nc variable
    *  \param [in] [NcVar *] pointer to the nc Variable
    *              [string] lcode
    *              [string] time
    *              [string] definition
    *              [string] unit
    *  \return nothing
    */    
    void add_attribute(NcVar *var, string lcode, string time, string def, string unit);
    
    /**
    *  \b Description: \n
    *        Function adding attributes to a Ncvar. 
    *  \param [in] [NcVar *] pointer to the nc Variable
    *              [vector<std::string>] attributes names
    *              [vector<std::string>] attributes values
    *  \return nothing
    */
    void add_attribute(NcVar *var, std::vector<std::string> &name,  std::vector<std::string> &val );
    
    /**
    *  \b Description: \n
    *        Function adding a char Variable (without attributes) to a nc file. 
    *  \param [in] [NcFile] NcFile by reference
    *              [string] name of nc variable
    *              [string] name of nc Dimension
    *              [string] content of nc variable
    *  \return nothing
    */
    void add_char_var(NcFile &dataFile, string ncVarName, string ncDimName, string content);
    
        /**
    *  \b Description: \n
    *        Function adding a char Variable with attributes to a nc file. 
    *  \param [in] [NcFile] NcFile by reference
    *              [string] name of nc variable
    *              [string] name of nc Dimension
    *              [string] content of nc variable
    *              [vector<std::string>] attributes names
    *              [vector<std::string>] attributes values
    *  \return nothing
    */
    void add_char_var(NcFile &dataFile, string ncVarName, string ncDimName, string content, std::vector<std::string> &name,  std::vector<std::string> &val);
    
      
    // saves all creates NcFiles
    vector<NcFile> _new_files;
    
    // the current working vgosdb
    string _directory; 
    
    ivg::Wrapper* _wrapper_ptr;

};


}

#endif	/* VGOSDB_H */

