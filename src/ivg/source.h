#ifndef SOURCE_H
#define	SOURCE_H

#include <map>

#include "logger.h"
#include "date.h"
#include "ivg_const.h"
#include "auxfunc.h"
#include "jpleph.h"
#include "norad.h"


/**
*
* @brief class Crf - base class
* @author AI - bakkari developer team
* @date 2015-03-24
* @version 0.1
* @todo DOXIGEN!!!
*/

namespace ivg
{

//Source-Name-Type
enum srcname { iers, icrf, ivs, ngs, jpl, MAXSRC };

//defines the type of source
enum srctype { source, moon, satellite };

//* Source Band Type 0.0  Flux   Baseline Flux   Baseline
//* Name          B       (Jy)    limit    (Jy)   limit
//   3C84   X     B  0.0  10.00   600.0    0.00   13000.0
//   3C84   S     B  0.0  10.00  1000.0    0.00   13000.0
// 0202+149 X     B  0.0  0.94   900.0     0.94    1530.0  0.82  2600.0  0.63  4420.0  0.52  7520.0  0.44 10400.0  0.44 12800.0

//* Source Band Type Flux   MajAx  Ratio  PA   Off1   Off2
//* Name        M    (Jy)   (mas)              (mas)  (mas)
//0003-066  X   M    3.82   0.80   0.40   10.   0.0   0.0

struct flux_info{
    string type; // B or M
    double flux; 
    double major_axis;
    double ratio;
    double pa;
    double off1,off2;
    ivg::Matrix baseline_flux;
    string flux_line; // stores flux-infos in a single line for easy skd-file writing
};

class Source
{

    public:

        Source();
        Source( ivg::srcname type, string value );
        Source( ivg::srcname type, string value, double ra0, double dec0);
        Source( ivg::srcname type, string value, int h, int m, double s, bool negative, int deg, int min, double sec );
        Source( ivg::srcname type, string value, double ra0, double dec0, double sig_ra0, double sig_dec0, double corr);
        Source( ivg::srcname type, string value, int h, int m, double s, bool negative, int deg, int min, double sec, double sig_ra_mas, double sig_dec_mas, double corr );
        
        void set_tle( tle_t tle)
        {
            _tle=tle;
        }
        void set_sp3( ivg::Matrix sp3)
        {
            _sp3=sp3;
        }
        void set_type( ivg::srctype type )
        {
            _type=type;
        }
        void set_vcs( bool value )
        {
            _vcs=value;
        }
        void set_defining( bool value )
        {
            _defining=value;
        }

        void set_ra0( double ra )
        {
            _ra0 = ra;
            _ra_series.at(0) = _ra0;
        }
        void set_dec0( double dec )
        {
            _dec0 = dec;
            _dec_series.at(0) = _dec0;
        }
  
        void set_bck_ra0(double ra) {_bck_ra0=ra;}

        void set_bck_dec0(double dec) {_bck_dec0=dec;}

        void set_ra_dec_bck()
        {
	  _ra0=_bck_ra0;
	  _dec0=_bck_dec0;
	}
	  
        void set_found( bool value )
        {
            _found=value;
        }
        void set_refepoch( ivg::Date epoch )
        {
            _epoch_series.at(0) = epoch;
        }
        void set_jpl_ephem( void **ephem )
        {
            _ephem = (*ephem);
        }
        
        void set_sigma_ra0dec0(double sig_ra0, double sig_dec0)
        {
            _sig_ra0 = sig_ra0;
            _sig_dec0 = sig_dec0;
        }
        
        void set_additional_information(ivg::Date first, ivg::Date mean, ivg::Date last, int n_sessions, int n_obs )
        {
            _first = first;
            _mean = mean;
            _last = last;
            _n_sessions = n_sessions;
            _n_obs = n_obs;
        }

        tle_t get_tle()
        {
            return(_tle);
        }
        
        ivg::Matrix get_sp3()
        {
            return(_sp3);
        }
        
        ivg::srctype get_type()
        {
            return(_type);
        }
        
        double get_ra0()
        {
            return(_ra0);
        }
        double get_dec0()
        {
            return(_dec0);
        }
  
        double get_sigma_ra0()
        {
            return(_sig_ra0);
        }
        double get_sigma_dec0()
        {
            return(_sig_dec0);
        }        
        double get_corr()
        {
            return(_corr);
        }
        ivg::Matrix get_sph()
        {
            ivg::Matrix sph(3,1,1);
            sph(0) = _ra0;
            sph(1) = _dec0;
            
            return sph;
        }
        bool get_found( )
        {
            return(_found);
        }
        bool is_special_handling( )
        {
            return(_special_handling);
        }
        bool is_defining( )
        {
            return(_defining);
        }
        int get_n_sessions( )
        {
            return(_n_sessions);
        }

        int get_n_obs( )
	{
            return(_n_obs);
        }
        int get_n_obs_in_sess( )
	{
            return(_n_obs_in_sess);
        }
        void increase_n_sessions(int plus)
        {
            _n_sessions += plus;
        }
        void increase_n_obs_in_sess(int plus)
        {
            _n_obs_in_sess += plus;
        }
        void set_ra0( int h, int m, double s );
        void set_dec0( bool negative, int deg, int min, double sec );
        void set_bck_ra0( int h, int m, double s );
        void set_bck_dec0( bool negative, int deg, int min, double sec );
        
        void add_local_position( double local_ra, double local_dec, ivg::Date epoch);

        void set_name( ivg::srcname type, string value );
        string get_name( ivg::srcname type = ivg::srcname::MAXSRC);

        // for satellites and moon
        ivg::Matrix get_vector_crs(ivg::Date epoch);
        
        // for regular source definied with ra/dec in SSB
        ivg::Matrix get_unit_vector_ssb();
        ivg::Matrix get_unit_vector_ssb_ga(double mjd, Setting * ga);
        ivg::Matrix get_unit_vector_ssb_partials_ra();
        ivg::Matrix get_unit_vector_ssb_partials_dec();

        bool check_name(string name);

        double calc_greenwich_hour_angle( ivg::Date gst );
        
        ivg::Matrix calc_rover_position( ivg::Date epoch );
        
        double calc_arclength_to_sun( ivg::Date epoch );
        
        double calc_arclength_to_source( ivg::Source& other );
        
        void get_position(int &h, int &m, double &s, int &deg, int &min, double &sec);
        
        void add_band_flux_info(ivg::band band, ivg::flux_info info);
        
        flux_info get_band_flux_info(ivg::band band){ return(_band_flux_info[band]); };

        double calc_flux(ivg::band band, double wave, ivg::Matrix bl, ivg::Date *epoch_ptr);

        bool& use_me(){ return _use_me; };
        
        void show();
        
        int get_idx() const {return _idx;}; 
        void set_idx( int idx) { _idx = idx;}; 

    private:

        //different names of the source: ivs, iers, icrf, ngs
        map<ivg::srcname,string> _names;
        
        // is it a moon, a satellite or just a regular source
        ivg::srctype _type;
        // in case of moon or satellite, the orbit information need to be stored
        // as TLE-elements
        tle_t _tle;
        // or in case of GPS-satellites as IGS-final-orbits
        ivg::Matrix _sp3;

        // ephems to caluclate sun-position and moon orbits
        void *_ephem;
        
        // reference position
        double _ra0;
        double _dec0;
        double _bck_ra0;
        double _bck_dec0;

        // standard deviation of the reference position
        double _sig_ra0;
        double _sig_dec0;
        double _corr;
        
        // in case of local source, there will be several positions with epochs
        vector<double> _ra_series;
        vector<double> _dec_series;
        vector<ivg::Date> _epoch_series;

        string _obj_type; //Galaxy, AGN, Quasar, etc
        bool _defining, _vcs, _found;
        bool _special_handling; // special handling

        ivg::Date _first, _mean, _last;
        int _n_sessions, _n_obs, _n_obs_in_sess;
        
        map<ivg::band, flux_info> _band_flux_info;

        bool _use_me;
        
        // index in CRF vector. Indices have to be created first with member function in CRF
        int _idx = -1;
};

} // # namespace ivg

#endif	/* SOURCE_H */

