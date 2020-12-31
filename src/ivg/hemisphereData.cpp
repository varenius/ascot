/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TroposphereData.cpp
 * Author: corbin
 * 
 * Created on January 18, 2019, 12:31 PM
 */

#include <vector>

#include "hemisphereData.h"
#include "troposphere.h"
#include "session.h"

namespace ivg{

HemisphereData::HemisphereData(){
    //std::vector<unsigned> k = {1, 10, 14, 20, 26, 30, 35, 38, 41, 42, 43}; // 300
    std::vector<unsigned> k = {1, 6, 12, 18, 24, 29, 35, 40, 44, 48, 51, 54, 56, 57, 58}; // 533
    
    ivg::Date d;
    d.now();
    
    HemisphereData(d, 3600, 3600, k);
}
    
HemisphereData::HemisphereData(ivg::Date startEpoch, lps::Seconds duration, lps::Seconds intervalLength, std::vector<unsigned> k ) {
    
    _cells = lps::createEqualAreaPartition( k );
    _tg = lps::TemporalGrid( duration, intervalLength, intervalLength, startEpoch );
    _data.resize( _cells.size(), _tg.get_number_of_intervals(), 0.0 );
   
   _idxOfFirstCell.resize(k.size());
   unsigned ksum = 0;
   for(int ring = 0; ring < k.size(); ++ring){
       _idxOfFirstCell[ring] = std::make_pair( _cells[  ksum  ].e1(),  ksum);
       ksum += k[ring];
   }
   
}

HemisphereData::HemisphereData(ivg::Date startEpoch, lps::Seconds duration, lps::Seconds intervalLength, std::vector<unsigned> k,
                                std::string path ) : HemisphereData( startEpoch, duration, intervalLength, k ) {

    this->load(path);
  
}

int HemisphereData::get_cell_idx( double az, double el ) const{
    lps::Point p(az*ivg::rad2d, el*ivg::rad2d);
    // find idx of first cell in the ring that contains the cell
    unsigned start_idx = 0;
    for( const std::pair<double, unsigned>& a : _idxOfFirstCell ){
        if( el*ivg::rad2d > a.first ){
            start_idx = a.second;
            break;
        }
    }
    
    // find cell containing point p
    for( unsigned cell_idx = start_idx ; cell_idx < _cells.size(); ++cell_idx){
        if( _cells[cell_idx].contains( p ) ){
                return (int)cell_idx;
        }
    }
    
    return -1;
}

double HemisphereData::get_data( double az, double el, unsigned tIdx ) const{
    int cell_idx = get_cell_idx(az, el);
    
    if( cell_idx < 0){
        throw runtime_error( "TroposphereData::get_data something went wrong");
    }
    
    return _data( cell_idx, (int)tIdx);

}

ivg::Matrix HemisphereData::get_data(double az, double el) const{
    int cell_idx = get_cell_idx(az, el);
    
    if( cell_idx < 0){
        throw runtime_error( "TroposphereData::get_data something went wrong");
    }
    
    return get_data( cell_idx );

}

ivg::Matrix HemisphereData::get_data( unsigned cell_idx ) const{
    return _data.get_sub(cell_idx, 0, cell_idx, _data.cols()-1);
}

double HemisphereData::get_data(double az, double el, lps::Seconds t) const{
    
    std::vector<unsigned> tIdxRet = _tg.getTemporalIndices(t);
    unsigned tIdx = tIdxRet[0];
    
    return get_data(  az,  el,  tIdx );
}

double HemisphereData::get_data( double az, double el, ivg::Date date) const{
    
    std::vector<unsigned> tIdxRet = _tg.getTemporalIndices(date);

    unsigned tIdx = tIdxRet[0];
    return get_data(  az,  el,  tIdx );
}

 ivg::Matrix HemisphereData::get_azel() const{
    ivg::Matrix azel(_cells.size(), 2);
    
    for( int i = 0; i < _cells.size(); ++i ){
        azel(i,0) = _cells[i].am();
        azel(i,1) = _cells[i].em();
    }
    
    azel*=ivg::d2rad;
    
    return azel;
 }
 
 void HemisphereData::getDimension( unsigned& npos, unsigned& nepochs ) const{
    npos = _cells.size();
    nepochs = _tg.get_number_of_intervals();
 }


void HemisphereData::save( const std::string path ) const{
    
    _data.save_bin( path + "/" + _name + "_troposphere.dat"  );
    save_grid(path);
}

 void HemisphereData::save_grid( const std::string path ) const{
    get_azel().save_bin( path + "/" + _name + "_troposphere_azel.dat"  );
    
  
    ivg::Matrix epoMat = getIntervalBeginningsMJD( );
    
    epoMat.save_bin( path + "/" + _name + "_troposphere_epoch.dat" );
}

void HemisphereData::load( const std::string path ){
    log<DETAIL> ("*** loading hemisphere data  from ") % path;
    _data.load_bin( path );

}

void HemisphereData::set_data( const ivg::Matrix& data){
    if( data.rows() != _data.rows() || data.cols() != _data.cols() ){
        std::cerr << "void HemisphereData::set_data. Dimension missmatch";
    }

    _data = data;
}

}