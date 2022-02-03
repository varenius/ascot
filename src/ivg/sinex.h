#include "logger.h"
#include "parser.h"
#include "matrix.h"
#include "ivg_const.h"
#include "date.h"
#include "param.h"
#include "trf.h"
#include "crf.h"
#include "eop_series.h"
#include "ls_neq.h"
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#ifndef SINEX_H
#define	SINEX_H

namespace ivg
{
    
enum reftype{ apriori, estimate, delta }; 
enum objtype{ trf, crf, eop };
    
struct statistics
{
    string name;
    string shorty;
    double value;
};
   
class Trf;
class Crf;
class Eop_series;

class Sinex {
    
public:
    
    /**
    *  \b Description: \n
    *        Default constructor
    *  \param [in] no input parameters needed
    *  \return An instance of the class 'Session'
    */
    Sinex( );

    /**
    *  \b Description: \n
    *        Copy constructor
    *  \param [in] [ivg::Scan] other observation
    *  \return An instance of the class 'Scan'
    */
    Sinex( const ivg::Sinex &other );
    /**
    *  \b Description: \n
    *        assignment operator '=' using another sinex
    *  \param [in] [Sinex] &other: other sinex
    */
    Sinex & operator=( const Sinex &other );
    
    /**
    *  \b Description: \n
    *       Constructor based on full path to the sinex file.
    *       If desired, reading of N and n_side can be skipped by setting 
    *       init_neq to false: performance boost!  
    *  \param [in] path to sinex file
    *  \return An instance of the class 'Sinex'
    */
    Sinex(string path, bool init_neq = true);
    
    /**
    *  \b Description: \n
    *       Constructor based on full path to the sinex_tro file.
    *       CAUTION: This function reads in a snx_tro file, but is a constructor for sinex.
    *       Function will not work on snx Files AND other constructors will not work on snx_tro.
    *       Intance can only call 
    *  \param [in] path to sinex file
    *              string s: random just to give unique constructor
    *  \return An instance of the class 'Sinex'
    */
    Sinex(string path, string tropo_snx);
    
    /**
    *  \b Description: \n
    *       Function to get am ivg::Trf based on the sinex file.
    *       You can choose apriori or estimate block by reftype  
    *  \param [in] reftype::apriori or reftype::estimate
    *  \return An instance of the class 'Trf'
    */
    ivg::Trf get_trf(reftype type);
    
    /**
    *  \b Description: \n
    *       Function to get am ivg::Trf based on the sinex file.
    *       You can choose apriori or estimate block by reftype  
    *  \param [in] reftype::apriori or reftype::estimate
    *  \return An instance of the class 'Trf'
    */
    ivg::Crf get_crf(reftype type);
    
    /**
    *  \b Description: \n
    *       Function to get am ivg::Trf based on the sinex file.
    *       You can choose apriori or estimate block by reftype  
    *  \param [in] reftype::apriori or reftype::estimate
    *  \return An instance of the class 'Trf'
    */
    ivg::Eop_series get_eop_series(reftype type);
    
    
    /**
    *  \b Description: \n
    *       Function to get start / end epoch
    *  \return An instance of the class ivg::Date
    */
    ivg::Date get_start_epoch(){ return _start; };
    ivg::Date get_end_epoch(){ return _end; };
    
    vector<ivg::Param> get_parameter_vetor(){ return _parameter;};
    
    ivg::Ls_neq get_neq(){return ivg::Ls_neq(_N,_n_side,_N.rows());};
    int get_num_objects(objtype type);
    
    double get_statistics( std::string shorty );
    
    // returns ivs names
    std::vector<std::string> get_station_names() const { return _station_names; } ;

    /**
    *  \b Description: \n
    *       Method to get the tropospheric delays for  
    *       each station based on the sinex file. 
    *  \param [in] [ivg::reftype] apriori (i.e., ZHD), estimate (i.e., ZTD) or estimate-total (i.e., ZWD)
    *  \return [map<string, ivg::Matrix>] tropospheric delays for each station
    */    
    map< string, ivg::Matrix > get_tropo_delays( reftype type );
    
    /**
    *  \b Description: \n
    *       Method to get the clock parameter (total values) for  
    *       each station based on the sinex file.
    *  \param [in] none
    *  \return [map<string, ivg::Matrix>] clock parameters for each station
    */    
    map< string, ivg::Matrix > get_clocks( std::string & type, int & max_order );
    
    
    std::set<std::string> get_ref_clock_station( );

    string get_path() {return _path;};
private:
    
    // funtion to replace D to E within a exponential value
    string _DtoE(string str);
    
    // function to convert std::string to double (std::stod gives only int if 
    // string is in scientific notation)
    double _stod(string str);
    
    // full path to snx file
    string _path;
    
    // stores epochs from header-line
    ivg::Date _start,_end;
    
    // stores all comments in FILE/COMMENT block
    std::stringstream _file_comment;
    
    // stores source information from SOURCE/ID block
    map<string, map<ivg::srcname,string> > _src_assignment;
    
    // stores station information from SITE/ID block
    map< string, map<ivg::staname,string> > _sta_assignment;
    
    std::vector<string> _station_names;
    
    // stores parameter estimates, aprioris, deviation, from ESTIMATE and APRIORI block
    vector<ivg::Param> _parameter;
    
    // stores information about discontinuties from EPOCH block
    map< string, vector<ivg::Date> > _discontinuities;
    
    // stores information about statistics from STATISTICS block
    map<string,statistics> _stats;
    
    // stores information about ne NEQ-system
    ivg::Matrix _N, _n_side;

};

 
}

#endif	/* SINEX_H */

