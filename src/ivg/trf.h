#ifndef TRF_H
#define	TRF_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

#include <vector>
#include <algorithm>
#include "matrix.h"

#include "logger.h"
#include "station.h"
#include "date.h"
#include "analysis_station.h"
#include "parser.h"
#include "auxfunc.h"
#include "sinex.h"

#include <cstdlib>
#include <libconfig.h++>

using namespace libconfig;

// necessary definitions for volume calculations
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
struct Plane_from_facet {
  CGAL::Polyhedron_3<K>::Plane_3 operator()(CGAL::Polyhedron_3<K>::Facet& f) {
      CGAL::Polyhedron_3<K>::Halfedge_handle h = f.halfedge();
      return CGAL::Polyhedron_3<K>::Plane_3( h->vertex()->point(),
                                    h->next()->vertex()->point(),
                                    h->opposite()->vertex()->point());
  }
};

/**
*
* @brief class Trf - base class
* @author AI - bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

namespace ivg
{

class Analysis_station;  // forward declaration

// ===========================================================================
class Trf
{
// ===========================================================================
    public:

        /**
        *  \b Description: \n
        *        Constructor just to define"
        *  @param[in] name
        *  \return An instance of the class 'Trf'
        *
        *  \b Example: \n  ivg::Trf itrf2008;
        */
        Trf( );
        
        /**
        *  \b Description: \n
        *        Constructor based on name, reference epoch and stations vector
         *       These are all trf-membervariables need to be set.
        */
        Trf( string name, ivg::Date ref_epoch, vector<ivg::Analysis_station> &stations );

        /**
        *  \b Description: \n
        *        Constructor using the TRF name and a vector with station names. TRF only contains chosen stations.
        *        The names of the stations are specified by the type, e.g. "ivs_name". Also selectable if displacements defined in stadisp
        *       should be loaded and added to the analysis stations or not. For NTAPL and VMF we need a time range start and end.
        *  @param[in] name
        *  @param[in] vector<string> station_names
        *  @param[in] type
        *  @param[in] decision: loading stadisps or not
        *  \return An instance of the class 'Trf'
        *
        *  \b Example: \n      (Pseudocode) \n
        *                      std::vector<std::string> station_names = { "WETTZELL", "WESTFORD", "GILCREEK", "NOTO" }; \n
        *                      ivg::Trf itrf2008("ITRF2008", station_names, "ivs_name");
        */
        Trf(  Setting &setup, const vector<string> station_names,
              const ivg::staname type, bool init_disps = true, ivg::Date start = ivg::Date(1970,1.0),
              ivg::Date end = ivg::Date(2070,1.0));
        
        
        void push_back(ivg::Analysis_station &new_station)
        {
            _stations.push_back(new_station);
        }
        

        /**
        *  \b Description: \n
        *        Function searches if a explicit station is present in the Trf.
        *        Any identKey ("ivs_name","domes_no", "cdp", etc.) can be used.
        *
        *   @param[in,out] station \n
        *   @param[in] identKey \n
        *   @param[in] identValu \n
        *
        *  \return bool \n if station was found: return is true and station pointer points to the found station\n
        *                  if station was not found: return is false and station pointer is not changed
        *
        *  \b Example: \n  ivg::Analysis_station * station; \n
        *               bool found = itrf2008.get_station(&station,"ivs_name","TSUKUB32");
        */
        bool get_station(ivg::Analysis_station **station, string name, ivg::staname maximal = MAXSTA);
        bool get_station(ivg::Analysis_station **station, string name, ivg::staname minimal, ivg::staname maximal){
            return(_get_station(_stations, station, name, minimal, maximal));
        };
        
        /**
        *  \b Description: \n
        *        Function generates a string vector containing the station names specified with the type.
        *        Any type ("ivs_name","domes_no", "cdp", etc.) can be used.
        *
        *   @param[in] type as string \n
        *
        *  \return vector<string> containing all station names of the stations included in the specific Trf \n
        *
        */
        vector<string>get_station_names(ivg::staname type);
        
        /**
        *  \b Description: \n
        *        Get indexes of corresponding stations in other-trf.
        *
        *   @param[in] ivg::Trf, another trf \n
        *
        *  \return ivg::Matrix, containing the indexes of corresponding stations in two columns for both trfs\n
        *
        */
        ivg::Matrix get_corresponding_stations(ivg::Trf &other);

        /**
        *  \b Description: \n
        *        Method to get the number of stations included in TRF
        *
        *  \param [in] no input parameters needed \n
        *
        * \return [int] number of stations
        */
        int get_number_stations() const
        {
            return _stations.size();
        }

        /**
        *  \b Description: \n
        *        Method to get the reference epoch
        *
        *  \param [in] no input parameters needed \n
        *
        * \return [ivg::Date] reference epoch
        */
        ivg::Date get_reference_epoch() const
        {
            return _ref_epoch;
        }

        /**
        *  \b Description: \n
        *        Method to get the TRF name
        *
        *  \param [in] no input parameters needed \n
        *
        * \return [std::string] TRF name
        */
        string get_name() const
        {
            return _name;
        }
        
        /**
        *  \b Description: \n
        *        Method to set the TRF name
        *
        *  \param [in] string name \n
        */
        void set_name(string name){_name = name;};

        /**
        *  \b Description: \n
        *        Method to set the reference epoch
        *
        *  \param [in] [ivg::Date] reference epoch \n
        */
        void set_reference_epoch( ivg::Date ref_epoch ){ _ref_epoch = ref_epoch; };
        
        /**
        *  \b Description: \n
        *        Method to get a station by a known index
        *
        *  \param [in] index of station \n
        *
        * \return [ivg::Analysis_station *] Analysis_station pointer
        */
        ivg::Analysis_station * get_station(int index)
        { 
            return &_stations.at(index); 
        };
        
        /**
        *  \b Description: \n
        *        Method to show(logger used) the data loaded within the specific analysis center
        */
        void log_data_info_table();
        
        /**
        *  \b Description: \n
        *        Method to initialize the station displacements within a given time range
        */
        void init_displacements(Setting &setup, ivg::Date start, ivg::Date end);
        
        /**
        *  \b Description: \n
        *        Method to initialize some sked catalogs
        */
        void init_sked_catalogs(Setting &setup);
        
        /**
        *  \b Description: \n
        *        Method to calculate the network volume based on convex hull
        *  \param [out] network volume in cubicmeter m³
        */
        double calculate_network_volume();
        
        /**
        *  \b Description: \n
        *        Method to keep a list of stations in trf, based on any staname type.
        *        All other stations will be erased from the trf
        * 
        *  \param [vector<string>] vector of station names with type of ivg::staname \n
        *
        */
        void keep_stations(vector<string> stations, ivg::staname type);
        
        /**
        *  \b Description: \n
        *        Method to remove an Analysis_station from TRF based on an index.
        * 
        *  \param [int idx] index of station in _stations which should be removed \n
        *
        */
        void remove_station( int idx );
        
        /**
        *  \b Description: \n
        *        Method to remove an Analysis_station from TRF based on an iterator
        * 
        *  \param [vector<ivg::Analysis_station>::iterator remove] iterator of station in _stations which should be removed \n
        *
        */
        void remove_station( vector<ivg::Analysis_station>::iterator remove);

        /**
        *  \b Description: \n
        *        Method to show the stations included in TRF
        *  \param [in] verbose - true or false
        */
        void show(bool verbose=false);

        vector<ivg::Analysis_station>::iterator begin()
        {
            return _stations.begin();
        };

        vector<ivg::Analysis_station>::iterator end()
        {
            return _stations.end();
        };
        
        void create_station_indices();
        
        std::vector< std::vector<ivg::Analysis_station*> > get_station_twin_list(){return _twin_stations;};
        
        void refresh_station_twin_list();
        
        void remove_from_station_twin_list(std::string station);
        
        bool areTwins(unsigned sta1, unsigned sta2);
        
        std::vector<ivg::Analysis_station*> get_twins(ivg::Analysis_station* sta);
        
        std::map<unsigned, std::vector<ivg::Analysis_station*> > get_twins_map();
        
        std::map<string, std::vector<ivg::Analysis_station*> > get_twins_map(ivg::staname type) ;
        
        // get_twins_map_including_all returns also stations that are not in the instance of the trf
        std::map<string, std::set<std::string> > get_twins_map_including_all() ;
       
    private:

        bool _get_station(vector<ivg::Analysis_station> &stations,ivg::Analysis_station **station, string name, ivg::staname minimal = MINSTA, ivg::staname maximal = MAXSTA);

        // _stations vector containing all station-specific information
        vector<ivg::Analysis_station> _stations;

        // _name like ITRF2008
        string _name;

        // _ref_epoch based on e.g. 01.01.2005
        ivg::Date _ref_epoch;
        
        std::vector< std::set<std::string> > _twin_station_names;
        
        std::vector< std::vector<ivg::Analysis_station*> > _twin_stations;

};

} // # namespace ivg

#endif	/* TRF_H */

