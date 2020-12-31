/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   wrapper.cpp
 * Author: corbin
 * 
 * Created on 25. April 2017, 09:16
 */

#include "wrapper.h"
#include "logger.h"
#include "auxfunc.h"
#include "parser.h"

#include <boost/regex.hpp>
#include "boost/format.hpp"
#include "ivg_const.h"
#include <iostream>


namespace ivg
{
// ...........................................................................
Wrapper::Wrapper() {
// ...........................................................................    
    _wrapper_found = false;
}

// ...........................................................................
Wrapper::Wrapper(std::string directory, std::string dbName, std::string editing)
// ...........................................................................
{
    _wrapper_found = false;
    _directory = directory;
    _dbName = dbName;
    _editing = editing;
    if(_get_latest_version()){
        _read_wrapper(_directory + _wrapper_name,  dbName, editing);
    }
    
}

// ...........................................................................
bool Wrapper::file_exists(ivg::wrapper_entries entry, ivg::band band){
// ........................................................................... 
    if(_wrapper_found && _association.count(entry) == 1)
    {
        if( _association[entry].count(band) == 1)
        {
            return true;
        }
    }
    
    log<WARNING> ("!!! no entry found for \"") % wrapper_entries_to_string(entry) % "\" (" % ivg::band_to_string(band) % "-band) in wrapper";
    return false;
}

// ...........................................................................
bool Wrapper::file_exists(ivg::wrapper_entries entry, std::string sta){
// ........................................................................... 
    if(_wrapper_found && _association_sta.count(entry) == 1)
    {
        if( _association_sta[entry].count(sta) == 1)
        {
            return true;
        }
    }
    
    log<WARNING> ("!!! no entry found for \"") % wrapper_entries_to_string(entry) % "\" " % sta % " in wrapper";
    return false;
}  

// ...........................................................................
bool Wrapper::_get_latest_version(short version)
// ...........................................................................
{
    std::vector< std::string > ls = list_local_dir(_directory);
    std::string latest = "";
    int max_ver = -1;
    int ver = -1;
    
    std::string ver_reg;
    if (version == 0)
    { // aribitrary version
        ver_reg = "V(\\d{3})";
    }
    else
    { // specific version
        std::stringstream ss;
        ss << boost::format("%3u") % version;
        ver_reg = "V(" + ss.str() + ")";
    }
        
    std::string end_reg = "kall\\.wrp";
    
    boost::regex wrapper_naming;
    
    boost::smatch m;
  
    bool wrapper_found = false;
    
    try{
    
        // two cases from specific to general
        for (int i = 0; i < 2 ; i++)
        {
            switch(i){
                case 0:
                    wrapper_naming = boost::regex(_dbName + "_" + ver_reg + "(?:_" + _editing + "_|_)" + end_reg);
                    break;
                case 1:
                    // arbitrary editing
                    wrapper_naming = boost::regex(_dbName + "_" + ver_reg + "(?:_.*_|_)" + end_reg);
                    log<WARNING> ("!!! No wrapper found for: ") % _dbName % " with editing: " % _editing % " and ending kall.wrp";
                    break;
            }

            //loop over all files in the directory
            for(std::string str : ls)
            {

                if(boost::regex_match(str,m,wrapper_naming))
                {
                    //get the number of the version (3 digits after V));
                    ver = stoi(m[1]);
                    // check if version number is higher than the former highest
                    if(ver > max_ver){
                       max_ver = ver;
                       latest = str;
                    }

                    wrapper_found = true;
                }
            }

            // end for loop if a file was found.
            if(wrapper_found)
                break;

        }
    }
    catch(boost::regex_error re){
        std::cerr << "boost::regex_error: " << re.what();
    }
    catch(std::exception& e){
        std::cerr << "std::exception: " << e.what();
    }

    if (wrapper_found)
    {
        log<DETAIL> ("*** Using latest wrapper file ") % latest;
        _version = max_ver;
        _wrapper_name = latest;
    }   
    else
        log<WARNING> ("!!! No wrapper found ") % version;
    
    _wrapper_found = wrapper_found;
    return wrapper_found;
        
}

// ...........................................................................
void Wrapper::_read_wrapper(std::string wrp_path, std::string dbName, std::string editing)
// ...........................................................................
{
        ifstream inStream;
        std::string line;
        std::string current_station; 
        try{
            
            unsigned short row = 0;
            while( ivg::parser::get_line(wrp_path,inStream, line))
            {   
                ++row;
                enum blocks {Session, Observation, Station,Scan};
                int current_block;
		
                // exclamation mark is comment symbol
                if(line.substr(0,1) != "!")
                {
                    // Detrmine the current block
                    if(line.find("Begin Session") !=std::string::npos){
                        current_block = blocks::Session;
                        continue;
                    }
                    else if(line.find("Begin Observation") !=std::string::npos){
                        current_block = blocks::Observation;
                        continue;
                    }
		     else if(line.find("Begin Station") !=std::string::npos){
                        current_block = blocks::Station;
			current_station = remove_spaces(line.substr(13,line.length()-13));
                        continue;
                    }
		     else if(line.find("Begin Scan") !=std::string::npos){
                        current_block = blocks::Scan;
                        continue;
                    }
                    // remove the ending (.nc)
                    std::size_t found = line.find(".nc");
                    if( found !=std::string::npos ){
                        line.erase(static_cast<int>(found),line.length());
                    }
                    // serach in current block for the files
                    switch(current_block)
                    {
                    case Session:
                        if( line.find("Head")  !=std::string::npos ){
                            _association[ivg::wrapper_entries::Head][ivg::band::X] = _create_wrapper_entry(line, row );
                            _association[ivg::wrapper_entries::Head][ivg::band::S] = _create_wrapper_entry(line, row );
                        }
			if( line.find("StationCrossRef")  !=std::string::npos ){
                            _association[ivg::wrapper_entries::StationCrossRef][ivg::band::X] = _create_wrapper_entry(line, row );
                            _association[ivg::wrapper_entries::StationCrossRef][ivg::band::S] = _create_wrapper_entry(line, row );
                        }
                        if( line.find("ClockBreak")  !=std::string::npos ){
                            _association[ivg::wrapper_entries::ClockBreak][ivg::band::X] = _create_wrapper_entry(line, row );
                            _association[ivg::wrapper_entries::ClockBreak][ivg::band::S] = _create_wrapper_entry(line, row );
                        }
			if( line.find("Source")  !=std::string::npos ){
                            _association[ivg::wrapper_entries::Sources][ivg::band::X] = _create_wrapper_entry(line, row );
                            _association[ivg::wrapper_entries::Sources][ivg::band::S] = _create_wrapper_entry(line, row );
                        }   
                        break;
                    case Observation:
                        for ( int i = 0; i <  ivg::band::MAXBANDTYPE; i++ )
                        {
			  if ( line.find("Default_Dir")==std::string::npos) {
                            if( line.find("Cal-SlantPathIonoGroup")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band)i ))!=std::string::npos )
                            {
                                _association[ivg::wrapper_entries::CalSlantPathIonoGroup][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if( line.find("Cal-SlantPathIonoPhase")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band)i ))!=std::string::npos )
                            {
                                _association[ivg::wrapper_entries::CalSlantPathIonoPhase][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
                            else if( line.find("CorrInfo")  !=std::string::npos && line.find( "_b" + ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::CorrInfo][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
                            else if( line.find("GroupDelayFull")  !=std::string::npos && line.find("_b" + ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {    
                                _association[ivg::wrapper_entries::GroupDelayFull][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if( line.find("GroupDelay")  !=std::string::npos && line.find("_b" + ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {    
                                _association[ivg::wrapper_entries::GroupDelay][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if( line.find("SBDelay")  !=std::string::npos && line.find("_b" + ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {
			     
                                _association[ivg::wrapper_entries::SBDelay][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if( line.find("SNR")  !=std::string::npos && line.find("_b" + ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {    
                                _association[ivg::wrapper_entries::SNR][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
                            else if(line.find("NumGroupAmbig")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {    
                                _association[ivg::wrapper_entries::NumGroupAmbig][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if(line.find("NumPhaseAmbig")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {    
                                _association[ivg::wrapper_entries::NumPhaseAmbig][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
                            else if(line.find("AmbigSize")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::AmbigSize][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
                            else if(line.find("GroupRate")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::GroupRate][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
                            else if(line.find("EffFreq")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::EffFreq][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if(line.find("QualityCode")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::QualityCode][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if(line.find("RefFreq")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::RefFreq][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if(line.find("PhaseDelayFull")  !=std::string::npos && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::PhaseDelayFull][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    else if(line.find("Phase_")  ==0 && line.find("_b"+ ivg::band_to_string((ivg::band) i) )!=std::string::npos )
                            {   
                                _association[ivg::wrapper_entries::Phase][(ivg::band) i] = _create_wrapper_entry(line, row, (ivg::band) i );
                            }
			    
                            else if(line.find("Edit")  !=std::string::npos ){
                                _association[ivg::wrapper_entries::Edit][ivg::band::X] = _create_wrapper_entry(line, row);
                                _association[ivg::wrapper_entries::Edit][ivg::band::S] = _create_wrapper_entry(line, row);
                            }
			    else if(line.find("TimeUTC")  !=std::string::npos ){
			      _association[ivg::wrapper_entries::TimeUTC][ivg::band::X] = _create_wrapper_entry(line, row);
			      _association[ivg::wrapper_entries::TimeUTC][ivg::band::S] = _create_wrapper_entry(line, row);
			    }
			    else if(line.find("ObsCrossRef")  !=std::string::npos ){
			      _association[ivg::wrapper_entries::ObsCrossRef][ivg::band::X] = _create_wrapper_entry(line, row);
			      _association[ivg::wrapper_entries::ObsCrossRef][ivg::band::S] = _create_wrapper_entry(line, row);
			    }
			  }
                            
                        }
                        break;
		    case Station:
		      {
			if( line.find("Cal-Cable")  !=std::string::npos)
                            {
                                _association_sta[ivg::wrapper_entries::CalCable][current_station] = _create_wrapper_entry(line, row );
                            }
			if( line.find("Met")  !=std::string::npos)
                            {
                                _association_sta[ivg::wrapper_entries::Met][current_station] = _create_wrapper_entry(line, row );
                            }
			if( line.find("TimeUTC")  !=std::string::npos)
                            {
                                _association_sta[ivg::wrapper_entries::TimeUTC][current_station] = _create_wrapper_entry(line, row );
                            }
			if( line.find("FeedRotation")  !=std::string::npos)
                            {
                                _association_sta[ivg::wrapper_entries::FeedRotation][current_station] = _create_wrapper_entry(line, row );
                            }
		      }
		    case Scan:
                        if( line.find("ScanName")  !=std::string::npos ){
                            _association[ivg::wrapper_entries::ScanName][ivg::band::X] = _create_wrapper_entry(line, row );
                            _association[ivg::wrapper_entries::ScanName][ivg::band::S] = _create_wrapper_entry(line, row );
                        }
                    }
                }


            }

            inStream.close();
        }
        catch(std::exception& e)
        {
            cerr << "std::exception: " << e.what()  << endl;
        }
        
        
        if (g_verbose > 3){
            std::cout << "*** The following file associations have been found:" << std::endl;
            for ( int j = 0; j <  ivg::band::MAXBANDTYPE; j++ ){
                std::cout << ivg::band_to_string((ivg::band)j) << "-band:"<< std::endl;
                for ( int i = 0; i <  ivg::wrapper_entries::MAX_WRAPPER_ENTRIES; i++ ){
                    if( _association[(wrapper_entries)i].count((ivg::band)j) == 1)
                        std::cout << "   " << wrapper_entries_to_string((wrapper_entries)i) << ": " << get_file((wrapper_entries)i,(ivg::band)j) << std::endl;
                }
            }
        }
          
        
}

// ...........................................................................
wrapper_entry Wrapper::_create_wrapper_entry(std::string ncfile, unsigned short row, bool new_version_created,bool newfile )
// ...........................................................................
{
   wrapper_entry  entry;
   entry.ncfile = ncfile;
   entry.row = row;
   entry.new_version_created = new_version_created;
   entry.new_file=newfile;
   return entry;
}

// ...........................................................................
std::string Wrapper::wrapper_entries_to_string(ivg::wrapper_entries entry )
// ...........................................................................
{
        switch(entry)
        {
            case ivg:: wrapper_entries::Head:
                return "Head";
	    case ivg:: wrapper_entries::ObsCrossRef:
	        return "ObsCrossRef";
	    case ivg:: wrapper_entries::StationCrossRef:
	        return "StationCrossRef";
	    case ivg:: wrapper_entries::Sources:
	        return "Sources";		
            case ivg:: wrapper_entries::CalSlantPathIonoGroup:
                return "Cal-SlantPathIonoGroup";
            case ivg:: wrapper_entries::CorrInfo:
                return "CorrInfo";
            case ivg:: wrapper_entries::GroupDelayFull:
                return "GroupDelayFull";
	    case ivg:: wrapper_entries::GroupDelay:
                return "GroupDelay";
	    case ivg:: wrapper_entries::SBDelay:
                return "SBDelay";
	    case ivg:: wrapper_entries::SNR:
                return "SNR";
            case ivg:: wrapper_entries::NumGroupAmbig:
                return "NumGroupAmbig";
	    case ivg:: wrapper_entries::AmbigSize:
                return "AmbigSize";
	    case ivg:: wrapper_entries::QualityCode:
                return "QualityCode";	
            case ivg:: wrapper_entries::GroupRate:
                return "GroupRate"; 
            case ivg:: wrapper_entries::Edit:
                return "Edit";
            case ivg:: wrapper_entries::ScanName:
                return "ScanName";
	    case ivg:: wrapper_entries::TimeUTC:
                return "TimeUTC";	
            case ivg:: wrapper_entries::EffFreq:
                return "EffFreq";
	    case ivg:: wrapper_entries::CalCable:
                return "CalCable";
	    case ivg:: wrapper_entries::ClockBreak:
                return "ClockBreak";
	    case ivg:: wrapper_entries::Met:
                return "Met"; 	
            default:
                return "n/a";
        }
}

void Wrapper::write_wrapper( std::string editing, unsigned short version)
{
    
    ifstream inStream;
    std::string line, output;
    
    try{
        
        // key: row in wrapper file to be replaced with string stored in value
        map<unsigned short, string> replace;
	
	vector<string> add_obs;
	std::map<string,vector<string>> add_sta;
        for( auto &a: _association)
        {
            for ( int i = 0; i <  ivg::band::MAXBANDTYPE; i++ )
            {
	      if( _association[a.first].count((ivg::band)i) == 1){
                if (a.second[(ivg::band)i].new_version_created){
		  std::size_t found = a.second[(ivg::band)i].ncfile.find(".nc");
		  if( found !=std::string::npos ){
		    replace[a.second[(ivg::band)i].row] = a.second[(ivg::band)i].ncfile;
		  } else {
		    replace[a.second[(ivg::band)i].row] = a.second[(ivg::band)i].ncfile+".nc";
		  }
		  
		}
	      }
	      if( _association[a.first].count((ivg::band)i) == 1)
                if (a.second[(ivg::band)i].new_file)
		  {
		    add_obs.push_back(a.second[(ivg::band)i].ncfile);
		  }
            }
            
        }
	for( auto &a: _association_sta)
        {
	  for( auto &b: a.second) {
            if (b.second.new_version_created){
	      std::size_t found = b.second.ncfile.find(".nc");
	      if( found !=std::string::npos ){
		replace[b.second.row] =  b.second.ncfile;
	      } else {
		replace[b.second.row] =b.second.ncfile+".nc";
	      }
	    }
	    if (b.second.new_file) {
	      add_sta[b.first].push_back(b.second.ncfile);
	      
	    }
	  }
        }
        
        unsigned short row = 0;
	bool new_obs_written=false;
        // read existing wrapper file and replace line if

        while( ivg::parser::get_line(_directory + _wrapper_name,inStream, line))
	  {  
            ++row;
	    
	    if (line.find("End Station") !=std::string::npos){
	      
	      string station = remove_spaces(line.substr(11,line.length()-11));
	      
	      for (std::map<string, vector<string>>::iterator it=add_sta.begin();it!=add_sta.end();it++)
	      	{
		 
	        if ((*it).first == station) {
		  
		  for (int i=0;i<add_sta[station].size();i++)
	      	      output+=add_sta[station][i]+"\n";
	      	  }
	      	}
	      
	    }
	    
	    if ((line.find("End Observation") !=std::string::npos)&&(!new_obs_written)){
	      char curdir=' ';
	      for (int i=0;i<add_obs.size();i++) {
		if( ((add_obs[i].find("Cal-SlantPathIonoGroup")!=std::string::npos)||(add_obs[i].find("EffFreq")!=std::string::npos))&&(curdir!='D')) {
		  output+="!\nDefault_Dir ObsDerived\n";
		  curdir='D';
		} else if (((add_obs[i].find("Edit")!=std::string::npos)||(add_obs[i].find("NumGroupAmig")!=std::string::npos)||(add_obs[i].find("GroupDelayFull")!=std::string::npos))&&(curdir!='E')){
		  output+="!\nDefault_Dir ObsEdit\n";
		  curdir='E';
		}
		output+=add_obs[i]+"\n";
		new_obs_written=true;
	      }
	      
	    }
            if(replace.count(row) == 1)
	      {
                output += replace[row] + "\n";
	      }
            else
	      {
                output += line+"\n";
           
		
	      }
	  }

        inStream.close();

        std::ofstream f;


        char ver [4];
	if (version>999) version=_version;
        sprintf(ver, "%03d",version );
        string wrapperFile = _dbName + "_V" + ver;
        if(!editing.empty())
        {
            wrapperFile = wrapperFile + "_" + editing;   
        }
        wrapperFile += "_kall.wrp";
        
        log<DETAIL> ("*** writing wrapper file") % _directory % wrapperFile;

        f.open(_directory + wrapperFile, ios::trunc);

        f << output;

        f.close();

    }
    catch(std::exception& e)
    {
        cerr << "std::exception: " << e.what()  << endl;
    }

}



}

