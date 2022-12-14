#ifndef MASTERFILE_H
#define	MASTERFILE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

#include <iostream>
#include <fstream>
#include "auxfunc.h"
#include "date.h"
#include "trf.h"
#include "analysis_station.h"

using namespace std;

namespace ivg {
    
    enum mastertype {
        intensive,
        regular,
        both
    };
    
    // struct storing all session realted information from master file
    struct sessinfo{
        string name; // IVS-R1104
        string code; // R1104
        string skdname; // r1104.skd on ftp://ivs.bkg.bund.de/pub/vlbi/ivsdata/aux/
        ivg::Date date; // 2004-JAN-05 17:00:00
        int duration; // 24h
        string stations; // BdFtHbHtKeKkNySvWwWzYgYj
        vector<string> stationnames; // string vector of lettercodes 
        ivg::Trf trf; // based on stations
        double volume;
        string sked; // NASA
        string corr; // BONN
        string dbc; // XA
        string dbname; // 04JAN05XA
        vector<string> groups; // IVS, R, R1
        string tooltip; // complete information from above seperated with \n
        mastertype type;
    };

    class Masterfile {

    public:
        

        /**
        *  \b Description: \n
        *        Default constructor
        */
        Masterfile();
        
        /**
        *  \b Description: \n
        *        Constructor initializing masterfile instance based on given directory.
        *        The directory must contain all masterXX.txt files between requestet time range.
        *        Timerange can be defined with start- and end-year. Directory has to contain ns-codes.txt file.
        *  \param [in] [string] directory to masterfiles
        *  \param [in] [int] start year, default 1992
        *  \param [in] [int] end year, default 2022
        */
        Masterfile(string directory,ivg::mastertype type ,int start_year = 1992, int end_year = 2022, string nscodespath="", string trfpath="", string trftype="");
        
        /**
        *  \b Description: \n
        *        Method to get all information of a single session based on DBNAME of session.
        *        The information is of struct type sessinfo. Further details see sessinfo struct defined in this headerfile. 
        *        Returns a default unknown sessinfo if session is not found.
        *  \param [in] [string] session name based on database name, e.g. 04JAN05XA, length of 9
        *  \return [ivg::sessinfo] sessinfo-struct containing information based on masterfile
        */
        ivg::sessinfo get_session_info(string dbname);
        
        /**
        *  \b Description: \n
        *        Method to get all information of a single session based on MIDEPOCH of session.
        *        The information is of struct type sessinfo. Further details see sessinfo struct defined in this headerfile. 
        *        Returns a default unknown sessinfo if session is not found.
        *  \param [in] [ivg::Date] mid-epoch of session
        *  \param [in] [double] threshold to match mid-epoch with sessions-epoch from database
        *  \return [ivg::sessinfo] sessinfo-struct containing information based on masterfile
        */
        ivg::sessinfo get_session_info(ivg::Date midsess, double hour_threshold = 4.0, mastertype type=mastertype::both);
        
        /**
        *  \b Description: \n
        *        Method to check whether session is of a specific group.
        *        Based on DBNAME (e.g. 04JAN05XA) of session. 
        *  \param [in] [string] session name based on database name, e.g. 04JAN05XA, length of 9
        *  \param [in] [string] group name, e.g. IVS-R, IVS-R1, IVS, CORE
        *  \return [bool] is of group or not
        */
        bool is_group(string dbname, string grpname);
        
        /**
        *  \b Description: \n
        *        Method to check whether session is of a specific group.
        *        Based on MIDEPOCH of session. 
        *  \param [in] [string] session name based on database name, e.g. 04JAN05XA, length of 9
        *  \param [in] [string] group name, e.g. IVS-R, IVS-R1, IVS, CORE
        *  \return [bool] is of group or not
        */
        bool is_group(ivg::Date midsess, string grpname);
        
        /**
        *  \b Description: \n
        *        Method to get the manually predefined groups of sessions
        *  \return [map] containing all groups and corresponding session types
        */
        map<string, vector<string> > get_groups(){ return _groups; };
        
        /**
        *  \b Description: \n
        *        Shows information about ALL initialized sessions.
        */
        void show();
        
        /**
        *  \b Description: \n
        *        Method to get the  sessions
        *  \return [vector*] containing  session
        */
        vector<sessinfo>* get_sessions(){ return &_sessions; };
        

    private:
        
        void parse_file(string path, int year, vector<ivg::Analysis_station> &all_stations, ivg::mastertype type);

	void parse_file_v2(string path, int year, vector<ivg::Analysis_station> &all_stations, ivg::mastertype type);
	
        // vector storing all session related information
        vector<sessinfo> _sessions;
        // map containing all manually defined group relationships
        map<string, vector<string> > _groups;
        // default sessinfo if a specific requested session is not found
        ivg::sessinfo _unknown;

    };

} // end ivg namespace

#endif	/* MASTERFILE_H */

