/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   wrapper.h
 * Author: corbin
 *
 * Created on 25. April 2017, 09:16
 */

#ifndef WRAPPER_H
#define WRAPPER_H

#include "ivg_const.h"

#include "logger.h"
#include "ivg_const.h"


namespace ivg
{
    
enum wrapper_entries{
    Head,
    ObsCrossRef,
    StationCrossRef,
    Sources,
    CalSlantPathIonoGroup,
    CalSlantPathIonoPhase,
    CorrInfo,
    GroupDelayFull,
    GroupDelay,
    SBDelay,
    SNR,
    NumGroupAmbig,
    NumPhaseAmbig,
    AmbigSize,
    QualityCode,
    GroupRate,
    Phase,
    PhaseDelayFull,
    RefFreq,
    EffFreq,
    ScanName,
    Edit,
    TimeUTC,
    CalCable,
    ClockBreak,
    Met,
    FeedRotation,
    MAX_WRAPPER_ENTRIES
};
            

struct wrapper_entry
{
    std::string ncfile; //file name of corresponding nc file
    unsigned short row; // row in wrapper file
    bool new_version_created; // true if another version of this file has been created
    bool new_file; // true a new nc-file has been added
}; 
    

class Wrapper {
public:
    Wrapper();
    Wrapper(string directory, string dbName, string editing);
    
    bool file_exists(ivg::wrapper_entries entry, ivg::band band);
    bool file_exists(ivg::wrapper_entries entry, std::string sta);
    std::string get_file(ivg::wrapper_entries entry, ivg::band band){return _association[entry][band].ncfile;};
    std::string get_file(ivg::wrapper_entries entry, std::string sta){return _association_sta[entry][sta].ncfile;};
    void set_file(ivg::wrapper_entries entry, ivg::band band, std::string file){ _association[entry][band].ncfile = file;};
  void set_file(ivg::wrapper_entries entry,  std::string station, std::string file){ _association_sta[entry][station].ncfile = file;};
    void set_created_flag(ivg::wrapper_entries entry, ivg::band band, bool b = true){ _association[entry][band].new_version_created = b;};
    void set_created_flag(ivg::wrapper_entries entry, std::string station, bool b = true){ _association_sta[entry][station].new_version_created = b;};
    
    bool isWrapper_found(){return _wrapper_found;};
    std::string get_dbName(){return _dbName;};
    
    std::string wrapper_entries_to_string(ivg::wrapper_entries entry );
    
    void write_wrapper(std::string editing, unsigned short version=9999);
    void add_wrapper_entry(ivg::wrapper_entries entry,ivg::band band,std::string ncfile)
    {
      _association[entry][band]=_create_wrapper_entry(ncfile,0,false,true);
    };
    void add_wrapper_entry(ivg::wrapper_entries entry,std::string sta,std::string ncfile)
    {
      _association_sta[entry][sta]=_create_wrapper_entry(ncfile,0,false,true);
    };
private:
    bool _get_latest_version(short version = 0);
    void _read_wrapper(std::string wrp_path, std::string dbName, std::string editing);
    
    wrapper_entry _create_wrapper_entry(std::string ncfile, unsigned short row, bool new_version_created = false , bool newfile=false );
    
    unsigned short _version;
    std::string _directory;   
    std::string _wrapper_name;
    bool _wrapper_found;
    
    std::string _dbName;
    std::string _editing;
    
    /*
     * first enum is wrapper_entries, second enum is ivg::band
     * wrapper_entry that are same for X and S band are stored in both;
     */
    std::map<ivg::wrapper_entries, std::map<ivg::band, wrapper_entry> > _association;
    std::map<ivg::wrapper_entries, std::map<std::string, wrapper_entry> > _association_sta;
};

}

#endif /* WRAPPER_H */

