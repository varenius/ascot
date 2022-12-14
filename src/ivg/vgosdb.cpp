#include "vgosdb.h"
#include "date.h"
#include "ivg_const.h"
#include "logger.h"
#include "auxfunc.h"
#include "wrapper.h"



namespace ivg
{
// ...........................................................................
Vgosdb::Vgosdb() {}
// ...........................................................................

// ...........................................................................
Vgosdb::Vgosdb(const string directory)
// ...........................................................................
{
    // vgosdb directory chosen
    _directory = directory;
    
    _wrapper_ptr = nullptr;
    
     // check if directory looking for is in place
    struct stat sb;
    if (stat(directory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)){
        
        string headpath = directory + "/" + VGOSDB_HEAD;
        if (stat(headpath.c_str(), &sb) != 0 && S_ISREG(sb.st_mode))
            throw runtime_error( "Vgosdb::Vgosdb(): Head.nc in vgosDB not existent: "+headpath );
    }
    else
        throw runtime_error( "Vgosdb::Vgosdb(): vgosDB not existent: "+directory );
}


// ...........................................................................
Vgosdb::Vgosdb(const string directory, ivg::Wrapper* wrapper_ptr) 
// ...........................................................................
{
    
    // vgosdb directory chosen
    _directory = directory;
    
    _wrapper_ptr = wrapper_ptr;
    
     // check if directory looking for is in place
    struct stat sb;
    if (stat(directory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)){
        
        string headpath = directory + "/" + _wrapper_ptr->get_file(ivg::wrapper_entries::Head,ivg::band::X);
        if (stat(headpath.c_str(), &sb) != 0 && S_ISREG(sb.st_mode))
            throw runtime_error( "Vgosdb::Vgosdb(): Head.nc in vgosDB not existent: "+headpath );
    }
    else
        throw runtime_error( "Vgosdb::Vgosdb(): vgosDB not existent: "+directory );
    
    
}

// ...........................................................................
vector<string> Vgosdb::get_string( string folder, string nc, string var )
// ...........................................................................
{
    // create path to file which need to be open
    string path = _directory + "/" + folder + "/" + nc + ".nc";  
    
    log<DETAIL>("*** read ") % var % " in " % path;
    
    // open the file for read access
    NcFile file(path.c_str(), NcFile::ReadOnly);
    
    if(!file.is_valid())
        throw runtime_error( "vector<string> Vgosdb::get_string( string, string, string): Could not open "+path);
    
    NcVar *tmp_var=file.get_var(var.c_str());
        
    // depending if its a single world or several worlds, we need to read the variable differently
    // 1. single word like "2013/06/11-16:51:23.0", returned in a vector
    if(tmp_var->num_dims() == 1 && tmp_var->get_dim(0)->size() >= 1)
    {
        int n = tmp_var->get_dim(0)->size(); // e.g. 5 chars for session-name
        
        char values[n];
        tmp_var->get(&values[0],n);
        
        vector<string> tmp = {values};
        
        return tmp;
    }
    // 2. several words like "Here","I","am", returned in a vector
    else if(tmp_var->num_dims() == 2 && tmp_var->get_dim(0)->size() >= 1 && tmp_var->get_dim(1)->size() >= 1)
    {
        
        int n1 = tmp_var->get_dim(0)->size(); // e.g. 7 stations
        int n2 = tmp_var->get_dim(1)->size(); // e.g. 8 chars for station-name
	
        char values[n1][n2];                
        tmp_var->get(&values[0][0],n1,n2);
        
        vector<string> out;
        for(int i=0; i<n1; ++i)
        {
            char tmp_name[n2];
            strncpy( tmp_name, values[i], n2 );
            string tmp_string(tmp_name);
	    
	    tmp_string.resize(n2);
            out.push_back(remove_spaces_end(tmp_string));
        }
        
        return out;
    }
    
}
// ...........................................................................
ivg::Matrix Vgosdb::get_matrix( string folder, string nc, string var )
// ...........................................................................
{
    // create path to file which need to be open
    string path = _directory + "/" + folder + "/" + nc + ".nc";
    
    log<DETAIL>("*** read ") % var % " in " % path;
    
    // open the file for read access
    NcFile file(path.c_str(), NcFile::ReadOnly);
    
    if(!file.is_valid())
        throw runtime_error( "ivg::Matrix Vgosdb::get_matrix( string, string, string ): Could not open "+path);
    
    NcVar *tmp_var=file.get_var(var.c_str());
    
    //check if there is the REPEAT attribute
    int repeat=0;
    for(int i=0; i<tmp_var->num_atts(); i++)
    {
        stringstream ss;
        ss << tmp_var->get_att(i)->name();
        if(ss.str() == "REPEAT")
            repeat = tmp_var->get_att("REPEAT")->as_int(0);
    }
      
    // e.g. getting Scan2Station matrix [3 1 4 5; 2 1 4 2; 1 2 5 2]
    if(tmp_var->num_dims() == 2 && tmp_var->get_dim(0)->size() >= 1 && tmp_var->get_dim(1)->size() >= 1)
    {
    
        int n1 = tmp_var->get_dim(0)->size();
        int n2 = tmp_var->get_dim(1)->size();

        double values[n1][n2];                
        tmp_var->get(&values[0][0],n1,n2);

        ivg::Matrix out(n1,n2,0);
        for(int i=0; i<n1; ++i)
        {
            for(int j=0; j<n2; j++)
                out(i,j) = values[i][j];
        }

        // in case of attribute REPEAT
        if(repeat>0)
        {
            ivg::Matrix repmat = out;
            for(int i=1; i<repeat; i++)
                out.append_rows(repmat);
        }
        
        return out;
    }
    // if attribute REPEAT is existent, the matrix will contain n-times the same value
    else if(tmp_var->num_dims() == 1 )
    {
        int repeats = tmp_var->get_att("REPEAT")->as_int(0);
        
        int n1 = tmp_var->get_dim(0)->size();

        double values[n1];                
        tmp_var->get(&values[0],n1);

        ivg::Matrix out(repeats,n1,0);
        for(int i=0; i<repeats; ++i)
        {
            for(int j=0; j<n1; j++)
                out(i,j) = values[j];
        }
        
        return out;
    }
    else
        throw runtime_error( "ivg::Matrix Vgosdb::get_matrix( string , string , string ): Dimension of variable ("+var+") unexpected in "+path );
    
}

// ...........................................................................
bool Vgosdb::does_file_exist( string folder, string nc )
// ...........................................................................
{
    // create path to file asked for
    string path = _directory + "/" + folder + "/" + nc + ".nc";     
    
    // open the file for read access
    NcFile file(path.c_str(), NcFile::ReadOnly);
    
    // check if file could be opened
    if(file.is_valid())
        return true;
    else
        return false;
}
// ...........................................................................
bool Vgosdb::does_variable_exist( string folder, string nc, string var)
// ...........................................................................
{
    // create path to file asked for
    string path = _directory + "/" + folder + "/" + nc + ".nc";     
    
    // open the file for read access
    NcFile file(path.c_str(), NcFile::ReadOnly);
    
    // check if file could be opened
    if(file.is_valid())
    {
        // check if desired variable is existent in opened file
        NcVar *tmp_var = NULL;
        for(int i=0; i<file.num_vars(); i++)
        {
            // if we find a variable with the name, the variable exists
            if(string(file.get_var(i)->name()) == var)
                return true;
        }
        
        // if we get here, the variable does not exist
        return false;
            
    }
    // if the file does not exist, the variable does not exist
    else
        return false;
}
// ...........................................................................
void Vgosdb::copy_file( string from_folder, string from_nc, string to_folder, string to_nc)
// ...........................................................................
{
    string from_path = _directory + "/" + from_folder + "/" + from_nc + ".nc";  
    string to_path = _directory + "/" + to_folder + "/" + to_nc + ".nc";  
        
    ifstream in (from_path); // open original file
    ofstream out(to_path); // open target file
 
    if(!in) 
        throw runtime_error( "Vgosdb::copy_file(  string from_folder, string from_nc, string to_folder, string to_nc ): Can't open source file for copying: "+from_path );
    
    if(!out)
        throw runtime_error( "Vgosdb::copy_file(  string from_folder, string from_nc, string to_folder, string to_nc ): Can't open target file for copying: "+to_path );
 
    // copy file by by buffering
    out << in.rdbuf();
    out.close();
    in.close();
    
    // change permissions
    chmod_urw_grw_or( to_path );
}
// ...........................................................................
NcFile* Vgosdb::create_file( string folder, string nc , NcFile::FileMode mode)
{
    string path = _directory + "/" + folder + "/" + nc + ".nc";  
    
    _new_files.push_back(NcFile(path.c_str(), mode));
    
    return(&(_new_files.at(_new_files.size()-1)));
}


// ...........................................................................
void Vgosdb::create_cal_file( vector<double> &cbc, string station, string session )
// ...........................................................................
{
   string path = _directory + "/" + station + "/Cal-Cable.nc";   
   
    
   // Create the file. The Replace parameter tells netCDF to overwrite
   // this file, if it already exists.
   NcFile dataFile(path.c_str(), NcFile::Replace);
   
   // check whether  netCDF file creation or open constructor succeeded.
   if (!dataFile.is_valid())
       log<WARNING>("!!! Couldn't open file within create_cal_file. No Cal-Cable.nc will be written. ") % path;
   else
   {
   
        // When we create netCDF dimensions, we get back a pointer to an
        // NcDim for each one.
        NcDim* nDim = dataFile.add_dim("NumScans", cbc.size());

        // Define zhe netCDF variable. The type of the variable in this case
        // is ncDouble 
        NcVar *cal_data = dataFile.add_var("Cal-Cable", ncDouble, nDim);
        
        //Creation time
        ivg::Date d;
        d.now();
        string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");
        
        add_attribute(cal_data, "CABL DEL",now, "Cable calibration data", "second");       

       
        // Write the  data to the file. 
        cal_data->put(&cbc[0], cbc.size());
               
        // add char vars to file
        add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
        add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
        add_char_var(dataFile, "DataOrigin", "DataOriginLen", "log");
        add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
        add_char_var(dataFile, "Session", "SessionLen", session);
        add_char_var(dataFile, "Station", "StationLen", station);
        add_char_var(dataFile, "Stub", "StubLen", "Cal-Cable");
        add_char_var(dataFile, "Subroutine", "SubroutineLen", "create_cal_file AC");
        add_char_var(dataFile, "TimeTag", "TimeTagLen", "StationScan");
        add_char_var(dataFile, "TimeTagFile", "TimeTagFileLen", "TimeUTC.nc");
        add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");
        
        vector<string>temp_name={"CreateTime"};
        vector<string>temp_val={now};
        add_char_var(dataFile, "Cal-CableHistory", "Char00060", "???",temp_name,temp_val);
       
        // The file will be automatically close when the NcFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.
	   if (_wrapper_ptr){
	     if(_wrapper_ptr->file_exists(ivg::wrapper_entries::CalCable,station)){
                _wrapper_ptr->set_created_flag(ivg::wrapper_entries::CalCable,station);
                _wrapper_ptr->set_file(ivg::wrapper_entries::CalCable,station,"Cal-Cable.nc");
	      }
	      else
		{
		  _wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::CalCable,station,"Cal-Cable.nc");
		  log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::CalCable) % " found in wrapper";
		}
	    }
	log<INFO>("*** success writing Cal-Cable.nc file for ") % station;
   }
   
   
    
}
// ...........................................................................
void Vgosdb::create_met_file( vector<double> &T,vector<double> &P,vector<double> &h, string station, string session )
// ...........................................................................
{
    
   if( T.size() == P.size() && T.size() == h.size())
   {
        string path = _directory + "/" + station + "/Met.nc";   


        // Create the file. The Replace parameter tells netCDF to overwrite
        // this file, if it already exists.
        NcFile dataFile(path.c_str(), NcFile::Replace);

        // check whether  netCDF file creation or open constructor succeeded.
        if (!dataFile.is_valid())
        {
           log<WARNING>("!!! Couldn't open file within create_met_file. No Met.nc will be written. ") % path;
        }
        else
        {

             // When we create netCDF dimensions, we get back a pointer to an
             // NcDim for each one.
             NcDim* nDim = dataFile.add_dim("NumScans", T.size());

             // Define zhe netCDF variable. The type of the variable in this case
             // is ncDouble 
             NcVar *pres_data = dataFile.add_var("AtmPres", ncDouble, nDim);
             NcVar *temp_data = dataFile.add_var("TempC", ncDouble, nDim);
             NcVar *relh_data = dataFile.add_var("RelHum", ncDouble, nDim);
             
             
            //Creation time
            ivg::Date d;
            d.now();
            string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");

            add_attribute(pres_data, "ATM PRES",now, "Pressure in hPa at site", "hPa");
            add_attribute(temp_data, "TEMP C  ",now, "Temp in C at local WX station", "Celsius");
            add_attribute(relh_data, "REL.HUM.",now, "Rel.Hum. at local WX st (50%=.5)", "%");

            pres_data->put(&P[0], P.size());
            temp_data->put(&T[0], P.size());
            relh_data->put(&h[0], P.size());
            
            // add char vars to file
            add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
            add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
            add_char_var(dataFile, "DataOrigin", "DataOriginLen", "log");
            add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
            add_char_var(dataFile, "Session", "SessionLen", session);
            add_char_var(dataFile, "Station", "StationLen", station);
            add_char_var(dataFile, "Stub", "StubLen", "Met");
            add_char_var(dataFile, "Subroutine", "SubroutineLen", "create_met_file AC");
            add_char_var(dataFile, "TimeTag", "TimeTagLen", "StationScan");
            add_char_var(dataFile, "TimeTagFile", "TimeTagFileLen", "TimeUTC.nc");
            add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");



             // The file will be automatically close when the NcFile object goes
             // out of scope. This frees up any internal netCDF resources
             // associated with the file, and flushes any buffers.
	    if (_wrapper_ptr){
	      if(_wrapper_ptr->file_exists(ivg::wrapper_entries::Met,station)){
                _wrapper_ptr->set_created_flag(ivg::wrapper_entries::Met,station);
                _wrapper_ptr->set_file(ivg::wrapper_entries::Met,station,"Met.nc");
	      }
	      else
		{
		  _wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::Met,station,"Met.nc");
		  log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::Met) % " found in wrapper";
		}
	    }


	    log<INFO> ("*** success writing Met.nc file for ") % station;
        }
    }
    else
    {
       log<WARNING> ("!!! Dimensions of AtmPres, TempC and RelHum not equal. Writing of Met.nc not possible");
    }
     
     
}

// ...........................................................................
void Vgosdb::create_NumGroupAmbig_file( vector<short> &ambig, string session, ivg::band type, string editing, bool phase )
// ...........................................................................
{   
    std::string band = ivg::band_to_string( type );
    
    if (!editing.empty()){
        editing = "_" + editing;
    }
    
    std::string file = "NumGroupAmbig" + editing + "_b"+ band + ".nc";
    if (phase)
      file = "NumPhaseAmbig" + editing + "_b"+ band + ".nc";
    std::string path = _directory + "/ObsEdit/" + file;
    
   // Create the file. The Replace parameter tells netCDF to overwrite
   // this file, if it already exists.
   NcFile dataFile(path.c_str(), NcFile::Replace);
   
   // check whether  netCDF file creation or open constructor succeeded.
   if (!dataFile.is_valid())
   {
     if (phase)
       log<WARNING>("!!! Couldn't open file within create_NumPhaseAmbig_file. No NumPhaseAmbig****.nc will be written. ") % path;
     else
       log<WARNING>("!!! Couldn't open file within create_NumGroupAmbig_file. No NumGroupAmbig****.nc will be written. ") % path;
   }
   else
   {
   
        // When we create netCDF dimensions, we get back a pointer to an
        // NcDim for each one.
        NcDim* nDim = dataFile.add_dim("NumObs", ambig.size());

        // Define zhe netCDF variable. The type of the variable in this case
        // is ncInt
	std::string VarName;
	if (phase)
	  VarName="NumPhaseAmbig";
	else
	  VarName="NumGroupAmbig";
        NcVar *ambig_data = dataFile.add_var(VarName.c_str(), ncShort, nDim);
        
        //Creation time
        ivg::Date d;
        d.now();
        string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");
        
        //add Attributes to the NcVar
        std::vector<std::string> name = {"LCODE","CreateTime","Band","Definition"};
        if (phase)
	  VarName="Number of phase delay ambiguities";
	else
	  VarName="Number of group delay ambiguities";
        std::vector<std::string> val = {"# AMBIG ",now,band,VarName};
        add_attribute(ambig_data,name,val);       

        // Write the  data to the file. 
        ambig_data->put(&ambig[0], ambig.size());
               
        // add char vars to file
        add_char_var(dataFile, "Band", "BandLen", band);
        add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
        add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
        add_char_var(dataFile, "DataOrigin", "DataOriginLen", " ");
        add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
        add_char_var(dataFile, "Session", "SessionLen", session);
	if (phase)
	  VarName="NumPhaseAmbig";
	else
	  VarName="NumGroupAmbig";
        add_char_var(dataFile, "Stub", "StubLen", VarName);
        add_char_var(dataFile, "Subroutine", "SubroutineLen", "resolve_ambiguity AC");
        add_char_var(dataFile, "TimeTag", "TimeTagLen", "Obs");
        add_char_var(dataFile, "TimeTagFile", "TimeTagFileLen", "TimeUTC.nc");
        add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");
       
        // The file will be automatically close when the NcFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.
        if (phase) {
	  if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::NumPhaseAmbig, type)){
	      _wrapper_ptr->set_created_flag(ivg::wrapper_entries::NumPhaseAmbig, type);
	      _wrapper_ptr->set_file(ivg::wrapper_entries::NumPhaseAmbig, type, file);
            }
            else
	      {
		_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::NumPhaseAmbig,type,file);
		log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::NumPhaseAmbig) % " found in wrapper";
	      }
	  }
	}
        else {
	  if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::NumGroupAmbig, type)){
	      _wrapper_ptr->set_created_flag(ivg::wrapper_entries::NumGroupAmbig, type);
	      _wrapper_ptr->set_file(ivg::wrapper_entries::NumGroupAmbig, type, file);
            }
            else
	      {
		_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::NumGroupAmbig,type,file);
		log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::NumGroupAmbig) % " found in wrapper";
	      }
	  }
	}
        
        log<INFO> ("*** success writing ") % path % " for session " % session ;
   }
}

// ...........................................................................
void Vgosdb::create_GroupDelayFull_file( vector<double> &gd_full, string session, ivg::band type, string editing, bool phase )
// ...........................................................................
{   
    std::string band = ivg::band_to_string( type );
    
    if (!editing.empty()){
        editing = "_" + editing;
    }
    
    std::string file = "GroupDelayFull" + editing +  "_b"+ band + ".nc";
    if (phase)
      file = "PhaseDelayFull" + editing +  "_b"+ band + ".nc";
    std::string path = _directory + "/ObsEdit/" + file;
    
   // Create the file. The Replace parameter tells netCDF to overwrite
   // this file, if it already exists.
   NcFile dataFile(path.c_str(), NcFile::Replace);
   
   // check whether  netCDF file creation or open constructor succeeded.
   if (!dataFile.is_valid())
   {
     if (phase)
       log<WARNING>("!!! Couldn't open file within create_PhaseDelayFull_file. No PhaseDelayFull****.nc will be written. ") % path;
     else
       log<WARNING>("!!! Couldn't open file within create_GroupDelayFull_file. No GroupDelayFull****.nc will be written. ") % path;
   }
   else
   {
        // When we create netCDF dimensions, we get back a pointer to an
        // NcDim for each one.
        NcDim* nDim = dataFile.add_dim("NumObs", gd_full.size());

        // Define zhe netCDF variable. The type of the variable in this case
        // is ncInt
	string tmp;
	if (phase)
	  tmp="PhaseDelayFull";
	else
	  tmp="GroupDelayFull";
        NcVar *gd_full_data = dataFile.add_var(tmp.c_str(), ncDouble, nDim);
        
        //Creation time
        ivg::Date d;
        d.now();
        string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");
        
        //add Attributes to the NcVar
        std::vector<std::string> name = {"CreateTime","Band","Definition","Units"};
	if (phase) {
	  tmp="Phase Delay";
	  //	  for (int j=0;j<gd_full.size();j++)
	  //  gd_full[j]*=1e6;
	}
	else
	  tmp="Group Delay";
        std::vector<std::string> val = {now,band,tmp,"seconds"};
        add_attribute(gd_full_data,name,val);       
	for (int j=0;j<gd_full.size();j++){
	  if (isnan(gd_full[j]))
	    gd_full[j]=0;
	}
        // Write the  data to the file. 
        gd_full_data->put(&gd_full[0], gd_full.size());
               
        // add char vars to file
        add_char_var(dataFile, "Band", "BandLen", band);
        add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
        add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
        add_char_var(dataFile, "DataOrigin", "DataOriginLen", " ");
        add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
        add_char_var(dataFile, "Session", "SessionLen", session);
	if (phase)
	  tmp="PhaseDelay";
	else
	  tmp="GroupDelayFull";
        add_char_var(dataFile, "Stub", "StubLen", tmp);
        add_char_var(dataFile, "Subroutine", "SubroutineLen", "resolve_ambiguity AC");
        add_char_var(dataFile, "TimeTag", "TimeTagLen", "Obs");
        add_char_var(dataFile, "TimeTagFile", "TimeTagFileLen", "TimeUTC.nc");
        add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");
       
        // The file will be automatically close when the NcFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.
        if (phase){
	  if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::PhaseDelayFull,type)){
	      _wrapper_ptr->set_created_flag(ivg::wrapper_entries::PhaseDelayFull,type);
                _wrapper_ptr->set_file(ivg::wrapper_entries::PhaseDelayFull,type,file);
            }
            else
	      {
		_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::PhaseDelayFull,type,file);
		log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::PhaseDelayFull) % " found in wrapper";
		
	      }
	  }
	}
	else {
	  if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::GroupDelayFull,type)){
	      _wrapper_ptr->set_created_flag(ivg::wrapper_entries::GroupDelayFull,type);
                _wrapper_ptr->set_file(ivg::wrapper_entries::GroupDelayFull,type,file);
            }
            else
	      {
		_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::GroupDelayFull,type,file);
		log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::GroupDelayFull) % " found in wrapper";
		
	      }
	  }
        }
        log<INFO> ("*** success writing ") % path % " for session " % session ;
   }
}


void Vgosdb::create_IonoGroup_file( ivg::Matrix &delta_tau_x, ivg::Matrix &delta_tau_x_sigma, std::vector<short> error_flag, 
                                    string session, ivg::band type,  std::string editing  , bool phase){
    
    std::string band = ivg::band_to_string( type );
    
    if (!editing.empty()){
        editing = "_" + editing;
    }
    
    std::string file = "Cal-SlantPathIonoGroup" + editing + "_b" + band + ".nc";
    if (phase)
      file = "Cal-SlantPathIonoPhase" + editing + "_b" + band + ".nc";
    std::string path = _directory + "/ObsDerived/" + file;
    
   // Create the file. The Replace parameter tells netCDF to overwrite
   // this file, if it already exists.
   NcFile dataFile(path.c_str(), NcFile::Replace);
   
   // check whether  netCDF file creation or open constructor succeeded.
   if (!dataFile.is_valid())
   {
       log<WARNING> ("Couldn't open file! ") % path;
   }
   else
   {
       
        const int NX = delta_tau_x.rows();
        const int NY = delta_tau_x.cols();
        // Create netCDF dimensions
        NcDim *xDim = dataFile.add_dim("NumObs", NX);
        NcDim *yDim = dataFile.add_dim("TimeDim2", NY);

        //---- delta tau and rate ----
	string tmp;
	if (phase)
	  tmp="Cal-SlantPathIonoPhase";
	else
	  tmp="Cal-SlantPathIonoGroup";
        NcVar *values = dataFile.add_var(tmp.c_str(), ncDouble, xDim, yDim);
      
        //Creation time
        ivg::Date d;
        d.now();
        string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");

        //add Attributes to the NcVar
        std::vector<std::string> name = {"LCODE","CreateTime","Band","Definition","Units"};
        std::vector<std::string> val = {"ION CORR", now, band, "Ion correction. Add to theo. sec", "second"};
        add_attribute(values,name,val);  
      
        values->put(&delta_tau_x.transpose().get_data_vec()[0], NX, NY);
        
        //---- sigma delta tau and sigma rate ----
        if (phase)
	  tmp="Cal-SlantPathIonoPhaseSigma";
	else
	  tmp="Cal-SlantPathIonoGroupSigma";
        NcVar *sigmas = dataFile.add_var(tmp.c_str(), ncDouble, xDim, yDim);
        
        //add Attributes to the NcVar
        name = {"LCODE","CreateTime","Band","Definition","Units"};
        val = {"IONRMS  ", now, band, "Ion correction to sigma. sec", "second"};
        add_attribute(sigmas,name,val);  
      
        sigmas->put(&delta_tau_x_sigma.transpose().get_data_vec()[0], NX, NY);
        
        //--- Flag
	if (phase)
	  tmp="Cal-SlantPathIonoPhaseDataFlag";
	else
	  tmp="Cal-SlantPathIonoGroupDataFlag";
        NcVar *flags = dataFile.add_var(tmp.c_str(), ncDouble, xDim);
        
        //add Attributes to the NcVar
        name = {"CreateTime","Band","Definition"};
        val = { now, band, "0=OK, -1=Missing,  -2=bad"};
        add_attribute(flags,name,val);  
      
        flags->put(&error_flag[0], NX, NY);
        
        //---- char vars ----
        add_char_var(dataFile, "Band", "BandLen", band);
        add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
        add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
        add_char_var(dataFile, "DataOrigin", "DataOriginLen", " ");
        add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
        add_char_var(dataFile, "Session", "SessionLen", session);
	if (phase)
	  tmp="Cal-SlantPathIonoPhase";
	else
	  tmp="Cal-SlantPathIonoGroup";
        add_char_var(dataFile, "Stub", "StubLen", tmp);
        add_char_var(dataFile, "Subroutine", "SubroutineLen", "calc_iono_cor AC");
        add_char_var(dataFile, "TimeTag", "TimeTagLen", "Obs");
        add_char_var(dataFile, "TimeTagFile", "TimeTagFileLen", "TimeUTC.nc");
        add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");
        
        // The file will be automatically close when the NcFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.

	if (phase){
	  if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::CalSlantPathIonoPhase,type)){
	      _wrapper_ptr->set_created_flag(ivg::wrapper_entries::CalSlantPathIonoPhase,type);
	      _wrapper_ptr->set_file(ivg::wrapper_entries::CalSlantPathIonoPhase,type,file);
            }
            else
	      {
		_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::CalSlantPathIonoPhase,type,file);
		log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::CalSlantPathIonoPhase) % " found in wrapper";
	      }
	  }
	}
	else{
	  
	  if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::CalSlantPathIonoGroup,type)){
	      _wrapper_ptr->set_created_flag(ivg::wrapper_entries::CalSlantPathIonoGroup,type);
	      _wrapper_ptr->set_file(ivg::wrapper_entries::CalSlantPathIonoGroup,type,file);
            }
            else
	      {
		_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::CalSlantPathIonoGroup,type,file);
		log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::CalSlantPathIonoGroup) % " found in wrapper";
	      }
	  }
        }
        log<INFO> ("*** success writing ") % path;
   }
    
}

void Vgosdb::create_EffFreq_file( std::vector<double> &eff_freq_iono_g, std::vector<double> &eff_freq_iono_ph,std::vector<double> &eff_freq_iono_rt,
                                    string session, ivg::band type,  std::string editing  ){
    
    std::string band = ivg::band_to_string( type );
    
    if (!editing.empty()){
        editing = "_" + editing;
    }
    
    std::string file = "EffFreq" + editing + "_b" + band + ".nc";
    std::string path = _directory + "/ObsDerived/" + file;
    
   // Create the file. The Replace parameter tells netCDF to overwrite
   // this file, if it already exists.
   NcFile dataFile(path.c_str(), NcFile::Replace);
   
   // check whether  netCDF file creation or open constructor succeeded.
   if (!dataFile.is_valid())
   {
       log<WARNING> ("Couldn't open file! ") % path;
   }
   else
   {
        // When we create netCDF dimensions, we get back a pointer to an
        // NcDim for each one.
        NcDim* nDim = dataFile.add_dim("NumObs", eff_freq_iono_g.size());

        // Define zhe netCDF variable. The type of the variable in this case
        // is ncInt
        NcVar *efffreq_g_data = dataFile.add_var("FreqGroupIono", ncDouble, nDim);
        
        //Creation time
        ivg::Date d;
        d.now();
        string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");
        
        //add Attributes to the NcVar
        std::vector<std::string> name = {"LCODE","CreateTime","Band","Definition","Units"};
        std::vector<std::string> val = {"GRIONFRQ",now,band,"Effective Group Delay Ionospheric Frequency","MHz"};
        add_attribute(efffreq_g_data,name,val);       

	NcVar *efffreq_p_data = dataFile.add_var("FreqPhaseIono", ncDouble, nDim);
	std::vector<std::string> namep = {"LCODE","CreateTime","Band","Definition","Units"};
        std::vector<std::string> valp = {"PHIONFRQ",now,band,"Effective Phase Delay Ionospheric Frequency","MHz"};
	add_attribute(efffreq_p_data,namep,valp);

	NcVar *efffreq_r_data = dataFile.add_var("FreqRateIono", ncDouble, nDim);
	std::vector<std::string> namer = {"CreateTime","Band","Definition","Units"};
        std::vector<std::string> valr = {now,band,"Effective Group Rate Ionospheric Frequency","MHz"};
	add_attribute(efffreq_r_data,namer,valr);

	
        // Write the  data to the file. 
        efffreq_g_data->put(&eff_freq_iono_g[0], eff_freq_iono_g.size());
        efffreq_p_data->put(&eff_freq_iono_ph[0], eff_freq_iono_ph.size());
	efffreq_r_data->put(&eff_freq_iono_rt[0], eff_freq_iono_rt.size());
        // add char vars to file
        add_char_var(dataFile, "Band", "BandLen", band);
        add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
        add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
        add_char_var(dataFile, "DataOrigin", "DataOriginLen", " ");
        add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
        add_char_var(dataFile, "Session", "SessionLen", session);
        add_char_var(dataFile, "Stub", "StubLen", "GroupDelayFull");
        add_char_var(dataFile, "Subroutine", "SubroutineLen", "LogReader");
        add_char_var(dataFile, "TimeTag", "TimeTagLen", "Obs");
        add_char_var(dataFile, "TimeTagFile", "TimeTagFileLen", "TimeUTC.nc");
        add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");
       
        // The file will be automatically close when the NcFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.
        
        if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::EffFreq,type)){
                _wrapper_ptr->set_created_flag(ivg::wrapper_entries::EffFreq,type);
                _wrapper_ptr->set_file(ivg::wrapper_entries::EffFreq,type,file);
            }
            else{
	      	_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::EffFreq,type,file);
                 log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::EffFreq) % " found in wrapper";
	    }
        }
        
        log<INFO> ("*** success writing ") % path % " for session " % session ;
   }
       
}

void Vgosdb::create_Edit_file( std::vector<int> &flag_g, std::vector<int> &flag_ph,std::vector<int> &flag_rt,
                                    string session,  std::string editing  ){
    
    
    
    if (!editing.empty()){
        editing = "_" + editing;
    }
    
    std::string file = "Edit" + editing  + ".nc";
    std::string path = _directory + "/ObsEdit/" + file;
    
   // Create the file. The Replace parameter tells netCDF to overwrite
   // this file, if it already exists.
   NcFile dataFile(path.c_str(), NcFile::Replace);
   
   // check whether  netCDF file creation or open constructor succeeded.
   if (!dataFile.is_valid())
   {
       log<WARNING> ("Couldn't open file! ") % path;
   }
   else
   {
        // When we create netCDF dimensions, we get back a pointer to an
        // NcDim for each one.
        NcDim* nDim = dataFile.add_dim("NumObs", flag_g.size());

        // Define zhe netCDF variable. The type of the variable in this case
        // is ncInt
        NcVar *flag_r_data = dataFile.add_var("RateFlag", ncDouble, nDim);
        
        //Creation time
        ivg::Date d;
        d.now();
        string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");
        
        //add Attributes to the NcVar
        std::vector<std::string> name = {"LCODE","CreateTime","Definition"};
        std::vector<std::string> val = {"RATUFLAG",now,"Delay rate unweight flag........"};
        add_attribute(flag_r_data,name,val);       

	NcVar *flag_p_data = dataFile.add_var("PhaseFlag", ncDouble, nDim);
	std::vector<std::string> namep = {"LCODE","CreateTime","Definition"};
        std::vector<std::string> valp = {"PHSUFLAG",now,"Phase Delay unweight flag"};
	add_attribute(flag_p_data,namep,valp);

	NcVar *flag_g_data = dataFile.add_var("DelayFlag", ncDouble, nDim);
	std::vector<std::string> namer = {"LCODE","CreateTime","Definition"};
        std::vector<std::string> valr = {"DELUFLAG",now,"Delay unweight flag"};
	add_attribute(flag_g_data,namer,valr);

	
        // Write the  data to the file. 
        flag_g_data->put(&flag_g[0], flag_g.size());
        flag_p_data->put(&flag_ph[0], flag_ph.size());
	flag_r_data->put(&flag_rt[0], flag_rt.size());
        // add char vars to file
        add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
        add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
        add_char_var(dataFile, "DataOrigin", "DataOriginLen", " ");
        add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
        add_char_var(dataFile, "Session", "SessionLen", session);
        add_char_var(dataFile, "Stub", "StubLen", "GroupDelayFull");
        add_char_var(dataFile, "Subroutine", "SubroutineLen", "LogReader");
        add_char_var(dataFile, "TimeTag", "TimeTagLen", "Obs");
        add_char_var(dataFile, "TimeTagFile", "TimeTagFileLen", "TimeUTC.nc");
        add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");
       
        // The file will be automatically close when the NcFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.
        
        if (_wrapper_ptr){
            if(_wrapper_ptr->file_exists(ivg::wrapper_entries::Edit,ivg::band::X)){
	      _wrapper_ptr->set_created_flag(ivg::wrapper_entries::Edit,ivg::band::X);
	      _wrapper_ptr->set_file(ivg::wrapper_entries::Edit,ivg::band::X,file);
            }
            else{
	      	_wrapper_ptr->add_wrapper_entry(ivg::wrapper_entries::Edit,ivg::band::X,file);
                 log<WARNING> ("!!! no entry for ")% _wrapper_ptr->wrapper_entries_to_string(ivg::wrapper_entries::Edit) % " found in wrapper";
	    }
        }
        
        log<INFO> ("*** success writing ") % path % " for session " % session ;
   }
       
}

void Vgosdb::create_ClockBreak_file( std::vector<string> stations, std::vector<double> epochs_mjd, string session ){
    
    if(this->does_file_exist("Session","ClockBreak"))
    {
        vector<string> stations_old = this->get_string("Session","ClockBreak","ClockBreakStationList");
        vector<double> epochs_mjd_old = this->get_vector<double>("Session","ClockBreak","ClockBreakEpoch");
        stations.insert( stations.end(), stations_old.begin(), stations_old.end() );
        epochs_mjd.insert( epochs_mjd.end(), epochs_mjd_old.begin(), epochs_mjd_old.end() ); 
    }
     
    
    
    
    std::string file = "ClockBreak.nc";
    std::string path = _directory + "/Session/" + file;
    
   // Create the file. The Replace parameter tells netCDF to overwrite
   // this file, if it already exists.
   NcFile dataFile(path.c_str(), NcFile::Replace);
   
   // check whether  netCDF file creation or open constructor succeeded.
   if (!dataFile.is_valid())
   {
       log<WARNING> ("Couldn't open file! ") % path;
   }
   else
   {
       
        const int M = stations.size();
        const int N = epochs_mjd.size();
        
        if( M != N ){
            log<WARNING> ("create_ClockBreak: number of stations and epochs do not match ") % path;
        }
        
        // convert vector of strings to array of chars
        char charField[N*8];
        for(int i = 0; i < N; ++i){
            const char* tmp = stations[i].c_str();
            for(int j = 0; j <  8; ++j){
                if( j < stations[i].size())
                    charField[i*8+j] = tmp[j];
                else
                    charField[i*8+j] = ' ';
            }
        }

        
        // Create netCDF dimensions
        NcDim* nDim = dataFile.add_dim("NumBreaks", N);
        NcDim* oDim = dataFile.add_dim("len", 1);
        NcDim* Dim8 = dataFile.add_dim("stalen", 8);
         

        //----  stations ----
        NcVar *sta = dataFile.add_var("ClockBreakStationList", ncChar, nDim, Dim8);
      
        //Creation time
        ivg::Date d;
        d.now();
        string now = d.get_date_time("YYYY/MO/DD HH:MI:SS");

        //add Attributes to the NcVar
        std::vector<std::string> name = {"LCODE","CreateTime","Definition"};
        std::vector<std::string> val = {"BRK_SNAM", now, "Batchmode clock break stations"};
        add_attribute(sta,name,val);  
      
        sta->put(&charField[0], N, 8);
        //---- epoch ----
        NcVar *epochs = dataFile.add_var("ClockBreakEpoch", ncDouble, nDim);
        
        //add Attributes to the NcVar
        name = {"LCODE","CreateTime","Definition"};
        val = {"BRK_SNAM", now, "Batchmode clock break epochs"};
        add_attribute(epochs,name,val);  

        epochs->put(&epochs_mjd[0], N);
        
        //---- BRK_NUMB ----
        NcVar *num = dataFile.add_var("BRK_NUMB", ncShort, oDim);
        
        //add Attributes to the NcVar
        name = {"LCODE","CreateTime","Definition"};
        val = {"BRK_NUMB", now, "Number of batchmode clock breaks"};
        add_attribute(num,name,val);  
      
        num->put(&N, 1);
        
       
        
        //---- char vars ----
        add_char_var(dataFile, "CreatedBy", "CreatedBy", "AC");
        add_char_var(dataFile, "CreatedTime", "CreatedTimeLen", now);
        add_char_var(dataFile, "DataOrigin", "DataOriginLen", " ");
        add_char_var(dataFile, "Program", "ProgramLen", "ivg::ASCOT");
        add_char_var(dataFile, "Session", "SessionLen", session);
        add_char_var(dataFile, "Stub", "StubLen", "ClockBreak");
        add_char_var(dataFile, "Subroutine", "SubroutineLen", "calc_iono_cor AC");
        add_char_var(dataFile, "vgosDB_Version", "vgosDB_VersionLen", "???");
        
        // The file will be automatically close when the NcFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.
       
        
        log<INFO> ("*** success writing ") % path;
   }
    
}


// ...........................................................................
void Vgosdb::add_attribute(NcVar *var, string lcode, string time, string def, string unit){
// ...........................................................................               

        std::vector<string> name = {"LCODE","CreateTime","Definition","Units"};
        std::vector<string> val = {lcode,time,def,unit};
        
        add_attribute(var, name,val );

}

// ...........................................................................
void Vgosdb::add_attribute(NcVar *var, std::vector<std::string> &name,  std::vector<std::string> &val ){
// ...........................................................................

    
    //ckeck whether the vectors have the same length
    if (name.size() != val.size())
        throw runtime_error( "void Vgosdb::add_attribute(): Vectors must have the same length.");
    else
    {
        // add attributes to NcVar
        for( int i = 0; i < name.size(); ++i){ 
            var->add_att(name[i].c_str(), val[i].c_str());
        }
    }

}
// ...........................................................................
void Vgosdb::add_char_var(NcFile &dataFile, string ncVarName, string ncDimName, string content){
 // ........................................................................... 
    
    //creating empty vector, so that no attributes are added in the called function
    std::vector<std::string> empty;
    // add the variable withou attributes
    add_char_var(dataFile, ncVarName, ncDimName, content, empty, empty);
}


// ...........................................................................
void Vgosdb::add_char_var(NcFile &dataFile, string ncVarName, string ncDimName, string content, std::vector<std::string> &att_name,  std::vector<std::string> &att_val ){
 // ...........................................................................   
        //TODO: check wether Nc var already exists!
    
        //convert the string content to a char variable
        const char *cstr = content.c_str();
         
        //define NcDimension
        NcDim* dim = dataFile.add_dim(ncDimName.c_str(), content.length());
        
        //create Nc variable
        NcVar *var = dataFile.add_var(ncVarName.c_str(), ncChar, dim);
        
        //adding the attributes
        if(att_name.size() > 0 && att_val.size() > 0){
            add_attribute(var, att_name,  att_val );
        }
        
        //Add var to file
        var->put(&cstr[0], content.length());
}




// ...........................................................................
string Vgosdb::_create_file_name(string nc, string band, bool extension, string editing, string version ){
// ...........................................................................
    string res = nc;
    
    // add editing if available e.g. "iIVS"
    if(!editing.empty())
    {
        res = res + "_" + editing;   
    }
    
    // add band type if availabl e.g. "bX"
    if(!band.empty())
    {
        res = res + "_b" + band;
    }
    
    // add version if available e.g. "V003"
    if(!version.empty())
    {
        res = res + "_" + version;   
    }
    
    // file extension .nc
    if(extension){
        res += ".nc";
    }
    
    return res;
}



}

