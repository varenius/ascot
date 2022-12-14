#include "masterfile.h"

namespace ivg {
   
// ...........................................................................
Masterfile::Masterfile()
// ...........................................................................
{
    // not 100% correct and complete right now
    map<string, vector<string> >  tmp ={
        {"APT",{"APT"}},
        {"APSG",{"APSG"}},
        {"AOV",{"AOV"}},
        {"AUSTRAL",{"AUS"}},
        {"AUST",{"A1"}},
        {"AUS-ALL",{"AUG","AUA"}},
        {"AUS-GEO",{"AUG"}},
        {"AUS-AST",{"AUA"}},
        {"CDP",{"CDP","CD","C9"}},
        {"IVS",{"OHI","CRD","CRF","R1","R4","T2","E3","RD"}},
        {"CORE",{"COH","CA","CB","C10","C3","CC","CNA"}},
        {"IVS-OHIG",{"OHI"}},
        {"CORE-OHIG",{"COH"}},
        {"CORE-A",{"CA"}},
        {"CORE-B",{"CB"}},
        {"CORE-1",{"C10"}},
        {"CORE-3",{"C3"}},
        {"CONT",{"C02","C05","C06","C08","C11","C14"}},
        {"IVS-CRD",{"CRD"}},
//        {"CONT96",{"CRD"}},
        {"CRF",{"CRF"}},
        {"CRMS",{"CRMS"}},
        {"CRL",{"CRL"}},
        {"EUROPE",{"EUR"}},
        {"IVS-R",{"R1","R4"}},
        {"IVS-R1",{"R1"}},
        {"IVS-R4",{"R4"}},
        {"IVS-T2",{"T2"}},
        {"IVS-E3",{"E3"}},
        {"IRIS",{"IS"}},
        {"JADE",{"JD"}},
        {"NAVEX",{"NAX"}},
        {"NEOS",{"NA","NB"}},
        {"NEOS-A",{"NA"}},
        {"NEOS-B",{"NB"}},
        {"SURVEY",{"SUR"}},
        {"WEST",{"WPC"}},
        {"VLBA",{"RV","RD","BK","BP","BB","BL"}},
        {"INT1",{"I"}},
        {"INT2",{"Q"}},
        {"UNDEFINED",{"XXX"}}
    };
    _groups = tmp;
    
    // define unknown 
    ivg::sessinfo unknown;
    unknown.name = "UNKNOWN";
    unknown.tooltip = "UNKNOWN";
    unknown.duration = 0;
    _unknown = unknown;
}  
// ...........................................................................
Masterfile::Masterfile(string directory, ivg::mastertype type, int start_year, int end_year,string nscodespath, string trfpath, string trftype) : Masterfile()
// ...........................................................................
{
    vector<string> empty;
    if (nscodespath=="")
      nscodespath=directory+"/ns-codes.txt";
    if (trfpath=="")
      {
	trfpath=directory+"/TRF.SSC";
	trftype="SSC";
	
      }
    vector< map<ivg::staname,string> > nscodes = ivg::parser::nscodes_parser( empty , ivg::staname::MAXSTA, nscodespath);
    // TRF.SSC only used for calculating network volume
     vector<ivg::Analysis_station> all_stations;
    if ( trftype == "SSC" )
      all_stations = ivg::parser::ssc_parser(trfpath, nscodes);
    else if (trftype == "SNX" )
      {
	ivg::Sinex snx(trfpath, false);
        ivg::Trf trf_esti = snx.get_trf(reftype::estimate);
	//trf_esti.keep_stations(nscodes, ivg::staname::ivs_name);
	
        for( int i=0;i<trf_esti.get_number_stations();i++ )
	  {
	    ivg::Analysis_station tmp= *trf_esti.get_station(i);
	    
	    for(auto &tmp_map: nscodes)
	      if(tmp_map[ivg::staname::domes_no] == tmp.get_name(ivg::staname::domes_no))
		tmp.set_name(ivg::staname::lettercode,tmp_map[ivg::staname::lettercode]);
			     
	    all_stations.push_back(tmp);
	    
	  }
	
      }
   
    // loop over all requested years of masterfiles
    for(int year=start_year; year<= end_year; year++)
    {
        // adjust year in order to find correct filename
      stringstream year_ss,year_ssss;
        if(year >= 1970 && year < 2000)
            year_ss << year-1900;
        else if(year >= 2000 && year < 2010)
            year_ss << "0" << year-2000;
        else if(year >= 2010)
            year_ss << year-2000;
        year_ssss << year;
        string path;
        switch(type){
            case intensive:
	       path = directory + "/master"+year_ssss.str()+"-int.txt";
	       if (file_exists(path))
		 parse_file_v2(path, year, all_stations, ivg::mastertype::intensive);
	       else {
                path = directory + "/master"+year_ss.str()+"-int.txt";
		if (file_exists(path))
		  parse_file(path, year, all_stations, ivg::mastertype::intensive);
	       }
	       break;
	    case regular:
	      path = directory + "/master"+year_ssss.str()+".txt";
	      if (file_exists(path))
		parse_file_v2(path, year, all_stations, ivg::mastertype::regular);
	      else {
		path = directory + "/master"+year_ss.str()+".txt";
		if (file_exists(path))
		  parse_file(path, year, all_stations, ivg::mastertype::regular);
	      }
	      break;
	    case both:
	      path = directory + "/master"+year_ssss.str()+".txt";
	      if (file_exists(path))
		parse_file_v2(path, year, all_stations, ivg::mastertype::regular);
	      else {
		path = directory + "/master"+year_ss.str()+".txt";
		if (file_exists(path))
		  parse_file(path, year, all_stations, ivg::mastertype::regular);
              }
	      path = directory + "/master"+year_ssss.str()+"-int.txt";
	      if (file_exists(path))
		parse_file_v2(path, year, all_stations, ivg::mastertype::intensive);
	      else {
		path = directory + "/master"+year_ss.str()+"-int.txt";
		if (file_exists(path))
		  parse_file(path, year, all_stations, ivg::mastertype::intensive);
	      }
	      break;
        }
 
    }    
}   




// ...........................................................................
void Masterfile::parse_file(string path, int year, vector<ivg::Analysis_station> &all_stations, ivg::mastertype type){
// ...........................................................................
ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "Masterfile::Masterfile(string directory, int start_year, int end_year): Failed to open file: "+path );
    }else{

        string line;
        while (getline(inStream,line,'\n'))
        {
            // skip all lines except the line holding the information
            if(line.substr(0,1) == "|")
            {
                string token;
                stringstream line_stream(line);
                vector<string> tokens;
                while (getline(line_stream, token, '|'))
                    tokens.push_back(token);

                // erase first trash token
                tokens.erase(tokens.begin());
		// For some new mastefiles, there is an empty token at the end
		// Quickfix:
                if(tokens.size() == 16)
		{
                    // erase last trash token
                    tokens.erase(tokens.end());
		}
                // only in case of correct number of fields extract the information
                if(tokens.size() == 15)
                {
                    sessinfo tmp;
                    tmp.name = remove_spaces(tokens.at(0));
                    tmp.code = remove_spaces(tokens.at(1));
                    
                    // create skd-file-name with lower character
                    string lower = tmp.code.substr(0,1); 
                    *lower.begin() = tolower(*lower.begin());
                    tmp.skdname = lower+tmp.code.substr(1);

                    double hour=0.0;
                    double min=0.0;
                    // sometimes in old masterfiles the exact epoch is not given (=empty)
                    if(tokens.at(4) != "     ")
                    {
                        hour = s2d(tokens.at(4).substr(0,2)) / 24.0;
                        min = s2d(tokens.at(4).substr(3,2)) / 60.0 / 24.0;
                    }
                    double doy = s2d(tokens.at(3)) +  hour + min;
                    tmp.date = ivg::Date(year,doy);
                    // also the duration might be empty
                    if(!(tokens.at(5) == "  " || tokens.at(5) == "   "))
                        tmp.duration = stoi(tokens.at(5));

                    tmp.stations = tokens.at(6);

                    // generate trf to be able to calculate network volume
                    string obs_stations = remove_spaces_end(tmp.stations.substr(0,tmp.stations.find_last_of("-")));

                    // in case of Va, replace it with all 10 individual stations
                    size_t posi = obs_stations.find("Va");
                    if (posi!=string::npos)
                        obs_stations.replace(posi,2,"BrFdGtHnKpLaMkNlOvSc");
                    
                    vector<string> lc_stations;
                    for(int i=0; i< obs_stations.size(); i+=2)
                        lc_stations.push_back(obs_stations.substr(i,2));        
                    
                    // stations in a vector using letter codes
                    tmp.stationnames = lc_stations;

                    vector<ivg::Analysis_station> sess_stations;
                    for(auto &sta: all_stations)
                        if(find(lc_stations.begin(), lc_stations.end(), sta.get_name(ivg::staname::lettercode)) != lc_stations.end())
                            sess_stations.push_back(sta);

                    tmp.sked = remove_spaces(tokens.at(7));
                    tmp.corr = remove_spaces(tokens.at(8));
                    tmp.dbc = remove_spaces(tokens.at(11));
                    tmp.dbname = tmp.date.get_date_time("YY")+remove_spaces(tokens.at(2))+tmp.dbc;

                    // initialize trf based on station-vector
                    tmp.trf = ivg::Trf(tmp.dbname, tmp.date, sess_stations);
                    if (! type == 0) //donnot calc volume for intensive sessions
                    {
                        tmp.volume = tmp.trf.calculate_network_volume()/1e18;
                    }

                    // find out which groups the current session is related to
                    vector<string>::iterator it;
                    for(auto &grp: _groups)
                    {
                        // the session code can be 2 or 3 letters
                        bool check_2l = find( grp.second.begin(), grp.second.end(), tmp.code.substr(0,2) ) != grp.second.end(); 
                        bool check_3l = find( grp.second.begin(), grp.second.end(), tmp.code.substr(0,3) ) != grp.second.end(); 
                        bool check_int = type==ivg::mastertype::intensive && find( grp.second.begin(), grp.second.end(), tmp.code.substr(0,1) ) != grp.second.end(); 
                        if( (check_2l || check_3l) || check_int )
                            tmp.groups.push_back(grp.first);
                    }

                    // if no group has been found
                    if(tmp.groups.empty())
                        tmp.groups.push_back("UNDEFINED");

                    //generate tooltip-string based on all selected information
                    stringstream tt;
                    tt << setfill(' ') << setw(8) << left << "DB:" << tmp.dbname << "\n";
                    tt << setfill(' ') << setw(8) << left << "Name:" << tmp.name << "\n";
                    tt << setfill(' ') << setw(8) << left << "Code:" << tmp.code << "\n";
                    tt << setfill(' ') << setw(8) << left << "Date:" << tmp.date.get_date_time("YYYY-MON-DD HH:MI:SS") << "\n";
                    tt << setfill(' ') << setw(8) << left << "Plan:" << tmp.stations << "\n";
                    tt << setfill(' ') << setw(8) << left << "Obs:" << tmp.stations.substr(0,tmp.stations.find_last_of("-")) << "\n";
                    tt << setfill(' ') << setw(8) << left << "Volume:" << fixed << setprecision(1) << tmp.volume << "\n";
                    tt << setfill(' ') << setw(8) << left << "Sked:" << tmp.sked << "\n";
//                        tt << setfill(' ') << setw(8) << left << "Corr:" << tmp.corr << "\n";
                    tt << setfill(' ') << setw(8) << left << "Groups:";
                    for(auto &grp: tmp.groups)
                        tt << grp << ",";

                    tmp.tooltip = tt.str();
                    tmp.type = type;
                    _sessions.push_back(tmp);
                }
                else
		{
	            // Debug print tokens in case it fails
		    std::string a = std::accumulate(tokens.begin(), tokens.end(), std::string(""));
                    throw runtime_error("Masterfile::Masterfile(string directory, int start_year, int end_year): Unexpected field length of " + std::to_string(tokens.size()) + " with tokens " + a + " in masterfile: "+path);
		}
            }
        }
    }
}

// ...........................................................................
void Masterfile::parse_file_v2(string path, int year, vector<ivg::Analysis_station> &all_stations, ivg::mastertype type){
// ...........................................................................
ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "Masterfile::Masterfile(string directory, int start_year, int end_year): Failed to open file: "+path );
    }else{

        string line;
	std::vector<std::string> month_names = { "JAN", "FEB", "MAR", "APR", "MAY",
                                                "JUN", "JUL", "AUG", "SEP", "OCT",
                                                "NOV", "DEC"
                                              };
        while (getline(inStream,line,'\n'))
        {
            // skip all lines except the line holding the information
            if(line.substr(0,1) == "|")
            {
                string token;
                stringstream line_stream(line);
		
                vector<string> tokens;
                while (getline(line_stream, token, '|'))
                    tokens.push_back(token);

                // erase first trash token
                tokens.erase(tokens.begin());
                
                // only in case of correct number of fields extract the information
                if(tokens.size() >= 13)
                {
                    sessinfo tmp;
                    tmp.name = remove_spaces(tokens.at(1))+"-"+remove_spaces(tokens.at(2));
                    tmp.code = remove_spaces(tokens.at(2));
		    
                    // create skd-file-name with lower character
                    string lower = tmp.code.substr(0,1); 
                    *lower.begin() = tolower(*lower.begin());
                    tmp.skdname = lower+tmp.code.substr(1);

                    double hour=0.0;
                    double min=0.0;
                    // sometimes in old masterfiles the exact epoch is not given (=empty)
                    if(tokens.at(4) != "     ")
                    {
                        hour = s2d(tokens.at(4).substr(0,2)) / 24.0;
                        min = s2d(tokens.at(4).substr(3,2)) / 60.0 / 24.0;
                    }
                    double doy = s2d(tokens.at(3)) +  hour + min;
                    tmp.date = ivg::Date(year,doy);
                    // also the duration might be empty
                    if(!(tokens.at(5) == "  " || tokens.at(5) == "     "))
                        tmp.duration = s2d(tokens.at(5).substr(0,2))+s2d(tokens.at(5).substr(3,2)) / 60.0;

                    tmp.stations = tokens.at(6);

                    // generate trf to be able to calculate network volume
                    string obs_stations = remove_spaces_end(tmp.stations.substr(0,tmp.stations.find_last_of("-")));

                    // in case of Va, replace it with all 10 individual stations
                    size_t posi = obs_stations.find("Va");
                    if (posi!=string::npos)
                        obs_stations.replace(posi,2,"BrFdGtHnKpLaMkNlOvSc");
                    
                    vector<string> lc_stations;
                    for(int i=0; i< obs_stations.size(); i+=2)
                        lc_stations.push_back(obs_stations.substr(i,2));        
                    
                    // stations in a vector using letter codes
                    tmp.stationnames = lc_stations;

                    vector<ivg::Analysis_station> sess_stations;
                    for(auto &sta: all_stations)
                        if(find(lc_stations.begin(), lc_stations.end(), sta.get_name(ivg::staname::lettercode)) != lc_stations.end())
                            sess_stations.push_back(sta);

                    tmp.sked = remove_spaces(tokens.at(7));
                    tmp.corr = remove_spaces(tokens.at(8));
                    tmp.dbc = remove_spaces(tokens.at(10));
                    //tmp.dbname = tmp.date.get_date_time("YY")+remove_spaces(tokens.at(2))+tmp.dbc;
                    tmp.dbname = tokens.at(1).substr(2,2)+month_names.at(s2d(tokens.at(1).substr(4,2))-1)+ tokens.at(1).substr(6,2)+tmp.dbc;
                    // initialize trf based on station-vector
                    tmp.trf = ivg::Trf(tmp.dbname, tmp.date, sess_stations);
                    if (! type == 0) //donnot calc volume for intensive sessions
                    {
                        tmp.volume = tmp.trf.calculate_network_volume()/1e18;
                    }

                    // find out which groups the current session is related to
                    vector<string>::iterator it;
                    for(auto &grp: _groups)
                    {
                        // the session code can be 2 or 3 letters
                        bool check_2l = find( grp.second.begin(), grp.second.end(), tmp.code.substr(0,2) ) != grp.second.end(); 
                        bool check_3l = find( grp.second.begin(), grp.second.end(), tmp.code.substr(0,3) ) != grp.second.end(); 
                        bool check_int = type==ivg::mastertype::intensive && find( grp.second.begin(), grp.second.end(), tmp.code.substr(0,1) ) != grp.second.end(); 
                        if( (check_2l || check_3l) || check_int )
                            tmp.groups.push_back(grp.first);
                    }

                    // if no group has been found
                    if(tmp.groups.empty())
                        tmp.groups.push_back("UNDEFINED");

                    //generate tooltip-string based on all selected information
                    stringstream tt;
                    tt << setfill(' ') << setw(8) << left << "DB:" << tmp.dbname << "\n";
                    tt << setfill(' ') << setw(8) << left << "Name:" << tmp.name << "\n";
                    tt << setfill(' ') << setw(8) << left << "Code:" << tmp.code << "\n";
                    tt << setfill(' ') << setw(8) << left << "Date:" << tmp.date.get_date_time("YYYY-MON-DD HH:MI:SS") << "\n";
                    tt << setfill(' ') << setw(8) << left << "Plan:" << tmp.stations << "\n";
                    tt << setfill(' ') << setw(8) << left << "Obs:" << tmp.stations.substr(0,tmp.stations.find_last_of("-")) << "\n";
                    tt << setfill(' ') << setw(8) << left << "Volume:" << fixed << setprecision(1) << tmp.volume << "\n";
                    tt << setfill(' ') << setw(8) << left << "Sked:" << tmp.sked << "\n";
//                        tt << setfill(' ') << setw(8) << left << "Corr:" << tmp.corr << "\n";
                    tt << setfill(' ') << setw(8) << left << "Groups:";
                    for(auto &grp: tmp.groups)
                        tt << grp << ",";

                    tmp.tooltip = tt.str();
                    tmp.type = type;
                    _sessions.push_back(tmp);
		    
                }
                else {
		  
                    throw runtime_error("Masterfile::Masterfile(string directory, int start_year, int end_year): Unexpected field length of in masterfile: "+path);
		}
            }
        }
    }
}
  
// ...........................................................................
ivg::sessinfo Masterfile::get_session_info(ivg::Date midsess, double hour_threshold, mastertype type)
// ...........................................................................
{
    double add = hour_threshold / 24.0;
    // get complete session info based on db_name (04JAN05XA)
    for(auto &sess: _sessions)
    {
        // only if type fits desired type
        if(type == mastertype::both || sess.type == type)
        {
            if(hour_threshold != 0.0)
            {
                double half_duration = ((double)sess.duration / 2.0)/24.0;
                if(midsess.get_double_mjd() > sess.date.get_double_mjd()+half_duration-add && midsess.get_double_mjd() < sess.date.get_double_mjd()+half_duration+add)
                    return sess;
            }
            else
            {
                // in case of 0.0, just check if it's the idential day, without considering fraction of day
                if(midsess.get_int_mjd() == sess.date.get_int_mjd())
                    return sess;
            }
        }
    }
    
    // if no session found, return unknown one    
    return _unknown;
}
// ...........................................................................
ivg::sessinfo Masterfile::get_session_info(string dbname)
// ...........................................................................
{
    //replace '-' in databasename with 'X'
    size_t posi = dbname.substr(0,8).find("-");
    if (posi!=string::npos)
        dbname.replace(posi,1,"X");
    
    // get complete session info based on db_name (04JAN05XA)
    for(auto &sess: _sessions)
    {
        if((sess.dbname == dbname)||(sess.name == dbname))
            return sess;
    }
    
    // if no session found, return unknown one
    log<WARNING> ("!!! No session found within masterfiles for dbname: ") % dbname % " are start and end year specified correct?";
    return _unknown;
}
// ...........................................................................
bool Masterfile::is_group(ivg::Date midsess, string grpname)
// ...........................................................................
{
    // get the groups which the database is related to
    vector<string> groups = get_session_info(midsess).groups;
    
    // check if database is in requested group
    bool found = find( groups.begin(), groups.end(), grpname ) != groups.end();
    
    return found;
}
// ...........................................................................
bool Masterfile::is_group(string dbname, string grpname)
// ...........................................................................
{
    // get the groups which the database is related to
    vector<string> groups = get_session_info(dbname).groups;
    
    // check if database is in requested group
    bool found = find( groups.begin(), groups.end(), grpname ) != groups.end();
    
    return found;
}
// ...........................................................................
void Masterfile::show()
// ...........................................................................
{
    for(auto &tmp_session: _sessions)
    {
        cerr << "-------------------" << endl;
        cerr << tmp_session.tooltip << endl;
    }
}

} // end ivg namespace


