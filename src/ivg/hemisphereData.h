/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TroposphereData.h
 * Author: corbin
 *
 * Created on January 18, 2019, 12:31 PM
 */

#ifndef HEMISPHEREDATA_H
#define HEMISPHEREDATA_H

#include "lpSked/grid.h"
#include "session.h"
#include "analysis_station.h"
#include "date.h"

namespace ivg{

class HemisphereData {
public:
    HemisphereData();
    HemisphereData( ivg::Date startEpoch, lps::Seconds duration, lps::Seconds intervalLength, std::vector<unsigned> k );
    HemisphereData( ivg::Date startEpoch, lps::Seconds duration, lps::Seconds intervalLength, std::vector<unsigned> k, std::string path );
    
    void set_data(const ivg::Matrix& data  );
    
    void set_name(const std::string name){ _name = name;};
    
    int get_cell_idx( double az, double el ) const;
    
    // az and el in radian
    double get_data( double az, double el, unsigned tIdx ) const;
    double get_data( double az, double el, lps::Seconds t ) const;
    double get_data( double az, double el, ivg::Date date) const;
    
    ivg::Matrix get_data( double az, double el) const;
    ivg::Matrix get_data( unsigned cell_idx ) const;
    
    unsigned get_number_of_cells( ) const { return _cells.size(); };
    
    ivg::Matrix get_azel() const;
        
    void save( const std::string path ) const;
    
    void save_grid( const std::string path ) const;
    
    void load( const std::string path );
    
    void getGrid( ivg::Matrix& azel, std::vector<ivg::Date>& epochs) const{
        azel = get_azel();
        epochs = _tg.get_intervalBeginningsDate();
    }
    
    ivg::Matrix getIntervalBeginningsMJD( ) const{
        
        std::vector<ivg::Date> epochs = _tg.get_intervalBeginningsDate();
        ivg::Matrix epoMat(epochs.size(), 1, 0.0);
        for(unsigned int i = 0; i < epochs.size(); ++i){
            epoMat(i) = epochs[i].get_double_mjd();
        }
        
        return epoMat;
    }
    
    void getDimension( unsigned& npos, unsigned& nepochs ) const;
private:
    
    //spacial grid
    std::vector<lps::Wedge> _cells;
    
    //temporal grid
    lps::TemporalGrid _tg;
    
    //each col corresponds to one time interval
    ivg::Matrix _data;
    
    std::string _name;
    
    // info to speed up lookup of delay
    //first: elevation, second: idx of first possible cell
    std::vector< std::pair<double, unsigned> > _idxOfFirstCell;
    

};
}
#endif /* HEMISPHEREDATA_H */

