#ifndef CRF_H
#define	CRF_H

#include <vector>
#include <iostream>
#include "matrix.h"
#include "date.h"
#include "parser.h"
#include "source.h"
#include "logger.h"

#include <cstdlib>
#include <libconfig.h++>

using namespace libconfig;

/**
*
* @brief class Crf - base class
* @author AI - bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

namespace ivg
{
        
class Crf
{

    public:

        /**
        *  \b Description: \n
        *  Default constructor
        *  \param [in] no input parameters needed
        *  \return An instance of the class 'Crf'
        */
        Crf();
        
        /**
        *  \b Description: \n
        *        constructor using crf-name and a vector of sources. This defines all crf-member-variables
        *  \param [in] [string name] crf-name, e.g. ICRF2 or 04JAN05XA or whatever
        *              [vector<ivg::Source>] vector of sources
        *  \return An instance of the class 'Crf'
        */
        Crf( string name, vector<ivg::Source> sources );

        /**
        *  \b Description: \n
        *        constructor using two input parameter: config file setup and source vector
        *  \param [in] [Setting] setup from config file
        *              [std::vector<ivg::Source>] vector of sources
        *              [void **ephem] to be able to calculate transformations 
        *  \return An instance of the class 'Crf'
        */
        Crf( Setting &setup, vector<ivg::Source> sources, void **ephem = NULL, std::string skdfile = "");
        
        /**
        *  \b Description: \n
        *        Method to push back a new ivg::Source at the end of the crf
        *  \param [in] [ivg::Source]
        */
        void push_back(ivg::Source &new_source)
        {
            _sources.push_back(new_source);
        }
        
        /**
        *  \b Description: \n
        *        Method to remove a ivg::Source at a given location (using an iterator).
        *  \param [in] [std::vector<ivg::Param>::iterator] location of the parameter to be removed
        */
        void remove_source( vector<ivg::Source>::iterator it );
                
        /**
        *  \b Description: \n
        *        Method to remove a source from CRF based on an index.
        * 
        *  \param [int idx] index of source in _sources which should be removed \n
        *
        */
        void remove_source( int idx );

        /**
        *  \b Description: \n
        *        Method to get a vector with the source names of a selected type
        *  \param [in] [ivg::Source] source
        *              [std::string] source name
        *  \param [out] [ivg::Source] source
        * \return [bool] true if source is found in CRF
        */
        bool get_source(ivg::Source **source, string name);
        
        /**
        *  \b Description: \n
        *        Method to get a source by a known index
        *
        *  \param [in] index of source \n
        *
        * \return [ivg::Source *] Source pointer
        */
        ivg::Source * get_source(int index)
        { 
            return &_sources.at(index); 
        };
        
        vector<ivg::Source>& get_sources() { return _sources; };

        /**
        *  \b Description: \n
        *        Get indexes of corresponding sources in other-trf.
        *
        *   @param[in] ivg::Crf, another crf \n
        *
        *  \return ivg::Matrix, containing the indexes of corresponding sources in two columns for both crfs\n
        *
        */
        ivg::Matrix get_corresponding_sources(ivg::Crf &other);
        
        /**
        *  \b Description: \n
        *        Method to get a vector with the source names of a selected type
        *  \param [in] [ivg::srcname] source type
        * \return [int] vector with the source names of a selected type
        */
        vector<string>get_source_names(ivg::srcname type);

        /**
        *  \b Description: \n
        *        Method to get the number of sources included in CRF.
        *        Ignoring use_me == false sources. 
        *  \param [in] no input parameters needed
        * \return [int] number of sources
        */
        int get_number_sources();
        
        int get_number_sources_inc_unused() const {return _sources.size();};
        
        /**
        *  \b Description: \n
        *        Method to show the sources in CRF
        *  \param [in] no input parameters needed
        */
        void show(bool full_info=false);

        vector<ivg::Source>::iterator begin()
        {
            return _sources.begin();
        };

        vector<ivg::Source>::iterator end()
        {
            return _sources.end();
        };
        
        // adds each source in _sources an attribute corresponding to its position in the vector
        void create_source_indices();
        
    private:

        // CRF name
        string _name;
        
        // vector of sources in the CRF
        vector<ivg::Source> _sources;

        void _ivssrc_parser(vector<ivg::Source> &sources, const string path,
			    bool set_pos);
        void _ocars_parser(vector<ivg::Source> &sources, const string path,
                           bool set_pos);
        void _icrf2_parser(vector<ivg::Source> &sources, const string path_sources);

        void _icrf3_parser(vector<ivg::Source> &sources, const string path_sources);

        bool _get_source(vector<ivg::Source> &sources, ivg::Source **source,
                         string name);

};

} // # namespace ivg

#endif	/* CRF_H */

