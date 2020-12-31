
/* 
 * File:   triangle.h
 * Author: corbin
 *
 * Created on 12. September 2017, 12:28
 */

#ifndef TRIANGLE_H
#define TRIANGLE_H


#include <array>
#include <iomanip>

#include "baseline.h"
#include "logger.h"

/*
 *           st(0)
 *  bl(2) /        \ bl(0)
 *       /          \
 *      v            v
 *    st(2) <-----  st(1)
 *           bl(1)
 * 
 */

class Triangle {
public:
   Triangle( std::array< Baseline* ,3> baselines, std::array<std::string, 3> stations);
   
   Triangle(const Triangle& orig);
   
   ~Triangle();
           
   void print_Triangle() const;
   
   void close_loop(std::string ref_sta);
      
       
    #if DEBUG
    static unsigned int count;
    unsigned int ID;
    #endif
    
private:
    
    double _calculateTriangleClosure( ) const;
    
    double _calculateTriangleClosure(  bool& satisfied ) const;
    
    std::array< Baseline* ,3> _baselines;
    std::array<std::string,3> _stations;
    
    double _threshold;
    unsigned int _n_obs;

};

#endif /* TRIANGLE_H */

