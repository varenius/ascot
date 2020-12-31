
/* 
 * File:   triangle.cpp
 * Author: corbin
 * 
 * Created on 12. September 2017, 12:28
 */

#include "triangle.h"

#if DEBUG
unsigned int Triangle::count = 0;
#endif

// ----------------------------------------------------------------------------
Triangle::Triangle( std::array< Baseline* ,3> baselines, std::array<std::string, 3> stations){
// ----------------------------------------------------------------------------
    #if DEBUG
    ID = count;
    std::cerr<<"+++ Triangle(  Baseline* bl1,  Baseline* bl2,  Baseline* bl3, std::array<std::string, 3> stations) ID: " << ID << std::endl;
    count++;
    #endif

    _stations = stations;
    _baselines = baselines;
    
    _n_obs = 0.0;
    for(unsigned int c = 0; c <3; c++){
        _n_obs += _baselines[c]->residuals.numel();
    }
    
    _threshold = _baselines[1]->ambiguity_spacing.max()/2;
}

// ----------------------------------------------------------------------------
Triangle::Triangle(const Triangle& orig){
// ----------------------------------------------------------------------------
    #if DEBUG
    this->ID = count;
    count++;
    std::cerr<<"+++ Triangle(const Triangle& orig) ID: " << ID << std::endl;
    #endif

    _baselines = std::array<Baseline*,3>(orig._baselines);
    _stations = orig._stations;
    
    _threshold = orig._threshold;
    _n_obs = orig._n_obs;
    
}

// ----------------------------------------------------------------------------
Triangle::~Triangle(){
// ----------------------------------------------------------------------------
    #if DEBUG
    std::cerr<<"--- ~Triangle() ID: " << ID << std::endl;
    #endif
}

// ----------------------------------------------------------------------------
void Triangle::close_loop(std::string ref_sta){
// ----------------------------------------------------------------------------
    #if DEBUG
    std::cerr <<"... Triangle::close_loop(std::string ref_sta)" << ID << std::endl;
    #endif
    
    bool satisfied;
    

    int ref = std::find (_stations.begin(), _stations.end(), ref_sta) - _stations.begin();
    switch(ref){
        case 0:{
	  _threshold = _baselines[1]->ambiguity_spacing.max()/2;
	}
        case 1:{
	  _threshold = _baselines[2]->ambiguity_spacing.max()/2;
	}
	case 2:{
	  _threshold = _baselines[0]->ambiguity_spacing.max()/2;
	}
    }
    double closureCondition = _calculateTriangleClosure( satisfied );
    stringstream lc, ts;
    lc << std::fixed  << std::setprecision(1) << closureCondition*1E+9 << " ns";
    ts << std::fixed  << std::setprecision(1) << _threshold*1E+9 << " ns";
    
    if( satisfied ){
        log<INFO> ("*** loop closure  | ") %  lc.str() % " | is smaller than " % ts.str() % " nothing to do";  
    } else {
        log<INFO> ("*** loop closure | ") % lc.str() % " | is lager than " % ts.str();  

        //int ref = std::find (_stations.begin(), _stations.end(), ref_sta) - _stations.begin();
        
        double shift = 0.0;
                
        switch(ref){
            case 0:{
                // 
                shift = -round( closureCondition / _baselines[1]->ambiguity_spacing.max() );
                
                _baselines[1]->shift += round(shift);
                log<INFO> ("*** 0 ") % _stations[1] % " - " % _stations[2] % " is shifted : " % round(shift) ;
                break;
            }
            case 1:{
                shift = +round( closureCondition / _baselines[2]->ambiguity_spacing.max() );
                
                _baselines[2]->shift += round(shift);
                log<INFO> ("*** 1 ") % _stations[0] % " - " % _stations[2] % " is shifted : " % round(shift) ;                
                break;
            }
            case 2:{
                shift = -round( closureCondition / _baselines[0]->ambiguity_spacing.max() );
                
                _baselines[0]->shift += round(shift);
                log<INFO> ("*** 2 ") % _stations[0] % " - " % _stations[1] % " is shifted : " % round(shift) ;
                break;
            }
            default:{
                log<WARNING> ("!!! the Reference Station \"") % ref_sta % "\" is not in this triangle.";
                break;
            }
        }
        
        
    }

}

// ----------------------------------------------------------------------------
double Triangle::_calculateTriangleClosure( ) const{
// ----------------------------------------------------------------------------
    
    return   _baselines[0]->residuals.meanD()+
             _baselines[1]->residuals.meanD()-
             _baselines[2]->residuals.meanD();
    
}   

// ----------------------------------------------------------------------------
double Triangle::_calculateTriangleClosure( bool& satisfied) const{
// ----------------------------------------------------------------------------
    
    double condition =  _calculateTriangleClosure( );
    
    satisfied = (abs (condition) < _threshold);
    
    return condition;

}   

// ----------------------------------------------------------------------------
void Triangle::print_Triangle() const{
// ----------------------------------------------------------------------------
    
    std::cout << "  " << _stations[0] << " " << _stations[1] << " "<< _stations[2] << std::endl; 
    std::cout << "    " << "observations: " << _n_obs << std::endl;
    
    bool satisfied;
    double closureCondition = _calculateTriangleClosure( satisfied );
    std::string color;
    if(satisfied)
        color = ivg::Logger::get_color("green");
    else
        color = ivg::Logger::get_color("boldred");
    
    stringstream ss;
    ss << std::fixed  << std::setprecision(1) << closureCondition*1E+9 << " ns";
    std::cout << color << "    " << "loop closure: " << ss.str() << ivg::Logger::get_color("white") << std::endl; 
   
    if(g_verbose > 4){
        
        std::cout << std::endl
                  << "          " << setw(8)  <<_stations[0] << std::endl
                  << "           /    \\" << std::endl
                  << "       " << std::setprecision(0) << setw(4) << _baselines[2]->residuals.meanD()*1E+9 << "      " << _baselines[0]->residuals.meanD()*1E+9 << std::endl
                  << "         /        \\"  << std::endl
                  << "        v          v" << std::endl
                  << "         <- "<< setw(4) <<_baselines[1]->residuals.meanD()*1E+9 << " ---" << std::endl
                  << "   " << setw(8) << _stations[2]  << "        " <<  _stations[1] << std::endl  << std::endl;
        
    }
       
}

