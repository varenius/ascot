/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   transits.h
 * Author: corbin
 *
 * Created on April 25, 2018, 9:52 AM
 */

#include "definitions.h"
#include "date.h"
#include "logger.h"
#include <map>

#include <vector>
#include <iostream>
#include "matrix.h"
#include "simulation.h"
#include "schedule.h"
#include "logger.h"
#include "session.h"
#include "analysis_station.h"
#include "crf.h"
#include "scan.h"
#include "obs.h"
#include "tictoc.h"
#include "auxfunc.h"
#include "geometry.h"
#include "grid.h"

#include <cstdlib>
#include <libconfig.h++>

#ifndef TRANSITS_H
#define TRANSITS_H

namespace ivg{
class Schedule ; //forward declaration
}

namespace lps{
        
/**
 * @brief The Postition class describes the position of a quasar with respect to
 * a given station and time.
 *
 */
class Position
{
public:
    Position(){};
    Position(double azimuth, double elevation, const ivg::Date & date, lps::Seconds minObsDur = 0.0);
    double azimuth() const {return azimuth_;} // in degree
    double elevation() const {return elevation_;} // in degree
    const  ivg::Date & date() const {return date_;}
    lps::Seconds minObsDur() const {return minObsDur_;}
    double dazi() const {return dazi_;}
    double dele() const {return dele_;}
    float snrx() const {return SNRxs_.first;}
    float snrs() const {return SNRxs_.second;}
    void set_minObsDur(lps::Seconds s)  { minObsDur_ = s;}
    void set_dazi(lps::Seconds dazi)  { dazi_ = dazi;}
    void set_dele(lps::Seconds dele)  { dele_ = dele;}
    
    void set_SNR(float x, float s)  { SNRxs_ = make_pair(x,s);}
    
    std::string toString() const {return "az=" +std::to_string(azimuth())+" el="+std::to_string(elevation());}
        
private:
    double azimuth_=-1; // in degree
    double elevation_=-1; // in degree
    double dazi_ = 0.0; // velocity in azimuth in degree/sec
    double dele_ = 0.0; // velocity in elevation in degree/sec
    ivg::Date   date_;
    
    lps::Seconds minObsDur_; // in sec
    std::pair<float,float> SNRxs_ = make_pair(0.0, 0.0);

};

/**
 * @brief The Transit class respresent the movement of a quasar with respect to
 * a given station and its view on the stars. Hence, each transit has a unique
 * quasar and a unique station. The transit is described by a path on the star map.
 */
class Transit
{       
public:
    
    // constructor
    Transit();

    Transit(int id, int sta_idx, int sou_idx, const ivg::Date& startEpoch, const ivg::Date& endEpoch,
            int start_idx, int end_idx, const std::vector<Position>& path, ivg::Date referenceEpoch,
            ivg::Session* session);
    
    Transit(int id, int sta_idx, int sou_idx, const ivg::Date& startEpoch, const ivg::Date& endEpoch, 
            int start_idx, int end_idx, const std::vector<Position>& path, ivg::Date referenceEpoch,
            ivg::Session* session, int visbibleFromStation);

    
    std::vector<Position> get_path_part(ivg::Date begin, ivg::Date end) const;
    
    // getter
    int id()      const  {return id_;}
    int get_sta_idx() const  {return sta_idx;}
    int get_scr_idx()  const  {return scr_idx;}
    int get_visbibleFromSta_idx()  const  {return visbibleFromStation;}
    const std::vector<Position>&  path() const {return path_;}
    std::vector<Position>&  path() {return path_;}

    // Reference is the start of the session
    const Seconds begin() const {return ( _startEpoch.get_double_mjd() - _referenceEpoch.get_double_mjd())*DurationDay; } 
    const Seconds end() const {return ( _endEpoch.get_double_mjd() - _referenceEpoch.get_double_mjd())*DurationDay; }
    
    Seconds duration() const {return end() - begin();}
    
    const ivg::Date get_startEpoch() const {return _startEpoch; };
    const ivg::Date get_endEpoch() const { return _endEpoch; };
    
    int get_start_idx() const{ return _start_idx;};
    int get_end_idx() const { return _end_idx;};
    
    // setter
    void set_endEpoch ( ivg::Date date ){ _endEpoch = date;}
    
    void set_referenceEpoch( ivg::Date date ){ _referenceEpoch = date;}
    
    std::string get_sou_name() const;
    
    void approxVelocity(double dt);

private:
    int  id_ =-1;
    int  sta_idx=-1; // transit is seen from station with this index
    int  scr_idx=-1;
    
    int visbibleFromStation = -1; // transit can be seen from this station too  from start until end;
    
    ivg::Date _startEpoch;
    ivg::Date _endEpoch;
    
   
    int _start_idx; // _start_idx * _dt = _startEpoch - session.start
    int _end_idx; // _end_idx * _dt = _endEpoch - session.start

    std::vector<Position> path_; // the path describes the transit of a quasar in az el
    
    ivg::Date _referenceEpoch; // Start of session
    
    std::string _source_name;
};


class Transits {
public:
    Transits();
    
    // in case dt = 0 dt is read from configfile
    Transits( ivg::Session *session_ptr, unsigned int dt = 0 );
    
    // Iterators for accessing std::vector<lps::Transit> container
   typedef range<typename std::vector<Transit>::const_iterator> TransitIterator;

   TransitIterator transit(int source_idx, int station_idx)const;

   // Getter
   const std::vector<Transit> & get() const {return container;}
   

   bool hasTransits(int source_idx, int station_idx) const;

   // returns all Transits that are visible from 2 stations
   Transits getBaselineWiseTransits( ) const;  
     
   Transits seenByAtLeastBaselines( unsigned n_bl) const; 
   Transits seenByAtLeastStations( unsigned n_sta) const;

   void print_transits() const;
   
   bool get_border_pos(const int sou_idx, const int sta_idx, const lps::Seconds begin, 
        const lps::Seconds end, const lps::Seconds intLen, std::pair<lps::Position, lps::Position>& pos) const;
   
   bool sourceIsVisible( const int sou_idx, const int sta_idx, const lps::Seconds t ) const;
   
   // if int max_level = 0 the value from the config file is used
   void compute_possible_coverage(std::map<int, std::vector<double> >& max_possible_coverage, std::map<int, std::vector<double> >& max_possible_surface_area, int max_level = 0 ) const;
   
   void save_transits(std::string path) const;
   
   void remove_source( int source_idx );
   
private:
    // maps a pair of stations and quasars on the corresponding transits
    // key(index station, index source)
    // value is first entry for sou/sta comibination in transits.container and last one+1 
    std::map<std::pair<int,int>, std::pair<int,int>> mapTransits_;
    
    // vector containg all transits
    std::vector<lps::Transit> container;
    ivg::Session* _session_ptr;
    
    unsigned int _dt; // time between two positions in Transits vector 
    

};



}

#endif /* TRANSITS_H */

