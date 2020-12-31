#include "trf.h"


namespace ivg
{

Trf::Trf( )
{

}

// ...........................................................................
Trf::Trf( string name, ivg::Date ref_epoch, vector<ivg::Analysis_station> &stations ) 
{
    _name = name;
    _ref_epoch = ref_epoch;
    _stations = stations;
}
// ...........................................................................
//Trf::Trf( vector<ivg::Analysis_station> &stations)
//{
//    _stations = stations;
//}
// ...........................................................................
Trf::Trf( Setting &setup, vector<string> station_names,
          const ivg::staname type, bool init_disps, ivg::Date start, ivg::Date end)
// ...........................................................................
{
#ifdef DEBUG_REFFRAME
   cerr << "+++ Trf::Trf( Setting , const vector<string> , const ivg::staname , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    string file_type = (const char *)get_list_element((setup)["refframes"],(setup)["trf"])[1];
    string file_path = (const char *)get_list_element((setup)["refframes"],(setup)["trf"])[2];
    
    if(init_disps)
        log<INFO>("*** Initializing ivg::Trf with ") % file_type % " file " % file_path % " including station displacements";
    else
        log<INFO>("*** Initializing ivg::Trf with ") % file_type % " file " % file_path % " WITHOUT station displacements";
    
    if (start > end)
        throw runtime_error( "Trf::Trf(...) Contructor: start-date later than end-date" );
   
    Setting &definitions = setup["definitions"];
    
    // save original order of station_names
    vector<string> orig_station_names = station_names;
            
    // replace unknown stations with known stations
    map<string,string> restore;
    Setting &identical_stations = definitions["identical_stations"];
    for(int i=0; i<identical_stations.getLength(); i++)
    {
        string unknown = (const char *)identical_stations[i][0];
        string known = (const char *)identical_stations[i][1];
        // check if unknown station is in current sessions, e.g. "MEDILIFT" to "MEDICINA" (15SEP07-L)
        vector<string>::iterator found_iter = find(station_names.begin(), station_names.end(), remove_spaces_end(unknown));
        if(found_iter != station_names.end())
        {
            log<INFO>("*** Duplicating station ") % unknown % " based on " % known % " due to definition of identical_stations";
            restore[unknown] = known;
            (*found_iter) = known;
        }
    }
    
    // in case of duplicates after "identical_stations", e.g. two times "MEDICINA", save indexes for erasing
    map<int,int> duplicates; // duplicates[unknown_idx] = known_idx
    for(int i=0; i<station_names.size(); i++)
        for(int j=i+1; j<station_names.size(); j++)
            if( station_names.at(i) ==  station_names.at(j))
                duplicates[j] = i;
        
    // erase duplicate temporary for displacement initilizations
    for (std::map<int,int>::reverse_iterator iter = duplicates.rbegin(); iter != duplicates.rend(); ++iter)
        station_names.erase(station_names.begin() + (*iter).first);
    
        
    Setting &twins = definitions["twin_stations"];
    _twin_station_names.resize(twins.getLength());
    
    for(int i=0; i<twins.getLength(); i++)
    {
        for(int j=0; j<twins[i].getLength(); j++){
            _twin_station_names[i].insert( (const char *)twins[i][j] );
        }
    }
        
    vector< map<ivg::staname,string> > nscodes = ivg::parser::nscodes_parser( station_names, type, definitions["nscodes"]);
    
    if( file_type == "SSC" )
        _stations = ivg::parser::ssc_parser(file_path, nscodes); // sets _ref_epoch
    else if( file_type == "SNX")
    {
        // to get a trf based on a snx file
        ivg::Sinex snx(file_path, false);
        ivg::Trf trf_esti = snx.get_trf(reftype::estimate);
        // erase stations we do not want
        trf_esti.keep_stations(station_names, type);
        // set this as snx-based trf
        (*this) = trf_esti;
    }
    
    _name = (const char *)setup["trf"];

    // if a station could not be initialized
    if(_stations.size() != station_names.size())
    {   
        for(auto &station : station_names)
        {
            ivg::Analysis_station * pointer_station;
	    
            // only compare station up to staname::description-level. Comparing staname::corres would be critical!!
            if(!_get_station(_stations, &pointer_station, station, ivg::staname::MINSTA, ivg::staname::description))
            {
                int cnt=0;
                for(auto &tmp_map: nscodes)
                {
                    if(tmp_map[type] == station)
                    {                        
                        log<DETAIL>("!!! CRITICAL: Not found in SSC-File. _pos0/_velo0 set to zero. Only ns-codes-names. Station initialized with ns-codes-names: ") % station;
       
                        vector<ivg::Date> none;
                        ivg::Date refepoch( ivg::fake_mjd );
                        Analysis_station unknown_station( ivg::Matrix(3,1,0.0), ivg::Matrix(3,1,0.0), refepoch, none, tmp_map);
                        _stations.push_back(unknown_station);
                        
                        break;
                    }
                    
                    cnt++;
                    
                    if(cnt == nscodes.size())
                    {
                        log<DETAIL>("!!! CRITICAL: Not found in SSC-File AND in ns-codes. _pos0/_velo0 set to zero. No ns-codes-names. Station initialized with single name: ") % station;
                                                
                        vector<ivg::Date> none;
                        ivg::Date refepoch( ivg::fake_mjd );
                        map<ivg::staname,string> tmp_names;
                        tmp_names[ivg::staname::ivs_name] = station;
                        Analysis_station unknown_station( ivg::Matrix(3,1,0.0), ivg::Matrix(3,1,0.0), refepoch, none, tmp_names);
                        _stations.push_back(unknown_station);
                        
                        break;
                    }
                }
            } 
        }
        
        if(_stations.size() != station_names.size()) {
	    std::cout << _stations.size() << " " << station_names.size() << endl;
	    for(auto &station : station_names) std::cout << station << endl;
	    for(auto &station : _stations) std::cout << station.get_name(ivg::staname::ivs_name) << endl;
            throw runtime_error( "Trf::Trf():Constructor: Not all stations could be initialized correctly!" );
	}
    }
    
    // add all activated displacements to the analysis_stations
    if(init_disps)
        init_displacements(setup, start, end);
    else
    {
        // in case of itrf2014 and if PSD is set to true, parse the psd_coefficients
        if( _name == "itrf2014" && (bool)get_list_element(setup["stadisp"], "PSD" )[1] )
        {
            log<WARNING>("!!! USING ITRF2014P - Correction of post-seismic deformation - SPECIAL itrf2014 CASE ");
            ivg::parser::psd_coefficients(this,(const char *)get_list_element(setup["stadisp"], "PSD" )[2]);
        }
    }
    
    // only add further information like equip.cat and antenna.cat if simulation is activated
    if( ( setup.exists( "SIM" ) && (bool)setup[ "SIM" ]["apply"]) || (setup.exists( "SKED" ) && (bool)setup[ "SKED" ]["apply"]) 
        || ((setup).exists( "SKED" ) && !(bool)setup[ "SKED" ]["apply"] &&  (bool)(setup)[ "SKED" ]["createPlots"] && ( !(setup).exists( "SIM" ) || !(bool)(setup)[ "SIM" ]["apply"]))
      )
    {
        // not in refernce sim step 1 and 2
        if(  ( setup[ "SIM" ].exists( "Reference" ) && setup[ "SIM" ]["Reference"]["apply"] && (int)(setup)["SIM"]["Reference"]["step"] < 2) == false  ){
        
            log<INFO>("*** Initializing ivg::Trf with sked catalogs.");
            init_sked_catalogs(setup);
        }
    }
        
    // restoration of erased station
    for (std::map<string,string>::iterator iter = restore.begin(); iter != restore.end(); ++iter)
    {
        // get existing known station
        ivg::Analysis_station *old_station;
        _get_station(_stations, &old_station, (*iter).second);
        
        ivg::Analysis_station new_station = (*old_station);
        new_station.set_name(ivg::staname::ivs_name, (*iter).first);
        _stations.push_back(new_station);
    }
    
    // store stations in same order as order from original station_names
    vector<ivg::Analysis_station> tmp_stations = _stations;
    vector<string> old_order_names = get_station_names(type);
    for(int old_idx=0; old_idx< old_order_names.size(); old_idx++ )
    {
        int new_idx = find(orig_station_names.begin(),orig_station_names.end(),old_order_names.at(old_idx)) - orig_station_names.begin();
        _stations.at(new_idx) = tmp_stations.at(old_idx);
    }
    
    // Finally fill twin list
    refresh_station_twin_list();
    
//    for( std::vector<ivg::Analysis_station*>& a : _twin_stations ){
//        for(ivg::Analysis_station* sta : a){
//            std::cout << sta->get_name(ivg::staname::ivs_name) << " ";
//        }
//        std::cout << std::endl;
//    }
    
#ifdef DEBUG_REFFRAME
   cerr << "--- Trf::Trf( Setting , const vector<string> , const ivg::staname , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Trf::init_sked_catalogs(Setting &setup)
// ...........................................................................
{
#ifdef DEBUG_REFFRAME
   cerr << "+++ void Trf::init_sked_catalogs(Setting &setup)" << endl; 
   tictoc tim;
   tim.tic();
#endif 
    
   Setting &definitions = setup["definitions"];
   
   // only sets the sked_antenna_name (ivg::staname::ant_name)
   ivg::parser::stations_cat(this, definitions["sked"]["stations"]);
   // sets further antenna info based on ivg::staname::ant_name
   // ivg::Antenna already set from antenna_info.txt but completed from antenna.cat (e.g. elevation/azimut limits and turn-rates)
   // seperated by "ANT_I" and "ANT_C"
   ivg::parser::antenna_cat(this, definitions["sked"]["antenna"]);
   // sets _sefd[band], known as "SEFD"
   ivg::parser::equip_cat(this, definitions["sked"]["equip"]);
   // sets _mask, known as "MASK", includes elevation-azimut relation
   ivg::parser::mask_cat(this, definitions["sked"]["mask"]);
   // sets _freq_sequences[band], known as "FREQ", similar for each analysiscenter in trf
   // returns rxname like "SX_WIDE" to get relation from rxcat to loifcat (e.g. SX_WIDE for WETTZELL -> CDP_WIDE) ;
   string rx_name = ivg::parser::freq_cat(this, definitions["sked"]["freq"], setup["SKED"]["freq_name"]);
   // we need all further catalogs in order to be able to WRITE a regular skd-file
   ivg::parser::loif_cat(this, definitions["sked"]["loif"], definitions["sked"]["rx"], rx_name);
   ivg::parser::tracks_cat(this, definitions["sked"]["tracks"], definitions["sked"]["rec"], setup["SKED"]["rec_name"]);
   ivg::parser::hdpos_cat(this, definitions["sked"]["hdpos"]);
   
#ifdef DEBUG_REFFRAME
   cerr << "--- void Trf::init_sked_catalogs(Setting &setup)" << " : " << tim.toc() << " s " << endl;
#endif    
}
// ...........................................................................
void Trf::init_displacements(Setting &setup, ivg::Date start, ivg::Date end)
// ...........................................................................
{
#ifdef DEBUG_REFFRAME
   cerr << "+++ void Trf::init_displacements(Setting &setup, ivg::Date start, ivg::Date end)" << endl; 
   tictoc tim;
   tim.tic();
#endif 
                  
    Setting &definitions = setup["definitions"];
    
    //Antenna Information
    // -> based on ivs_name
    ivg::parser::antenna_info(this, definitions["stadisp"]["antenna"]);
    //Eccentricity Values
    // -> based on cdp
    ivg::parser::ecc(this, definitions["stadisp"]["ecc"]);
    //Gravititation deformation coefficients
    ivg::parser::gravdef(this, definitions["stadisp"]["grav_deform"], start);
    //Ocean Pole Tide Loading Coefficients
    // -> based on ivs_name
    if((bool)get_list_element(setup["stadisp"], "OCEAN POLE TIDE LOADING" )[1])
        ivg::parser::optl(this, definitions["stadisp"]["optl"]);
    //Ocean Loading Coefficients
    // -> based on ivs_name
    if((bool)get_list_element(setup["stadisp"], "OCEAN LOADING" )[1])
        ivg::parser::blq(this, definitions["stadisp"]["ol"][(const char *)get_list_element(setup["stadisp"], "OCEAN LOADING")[2]]);
    //Non-Tidal Athmospheric Pressure Loading
    // -> based on ivs_name or corres (corres not stored anymore!! only used once here) (from vlbi_to_vsgd.inp file)
    if((bool)get_list_element(setup["stadisp"], "NON TIDAL APLO" )[1])
        ivg::parser::bindisp(this, definitions["stadisp"]["ntapl"], ivg::parser::correspondence(definitions["corres"]), start, end);
    //Tidal Athmospheric Pressure Loading
    // -> based on difference less than 3km in case of hps
    if((bool)get_list_element(setup["stadisp"], "TIDAL APLO" )[1])
    {
        string tapl_path = definitions["stadisp"]["tapl"];
        if(tapl_path.substr(tapl_path.size()-3) == "dat" )
            ivg::parser::dat(this, definitions["stadisp"]["tapl"]);
        else
            ivg::parser::hps(this, definitions["stadisp"]["tapl"]);
    }
    // Mapping functions (VMF)
    // -> based on ivs_name
    if((bool)setup.exists("troposphere"))
    {
        if( string((const char*)setup["troposphere"]["mapping_function"]) == "vmf1")
            ivg::parser::external_met_data(this, definitions["troposphere"]["mapping_function"]["vmf1"], start, end, ivg::extdata::MAPPING);
	if( string((const char*)setup["troposphere"]["mapping_function"]) == "vmf3")
            ivg::parser::external_met_data(this, definitions["troposphere"]["mapping_function"]["vmf3"], start, end, ivg::extdata::MAPPING3);
    }
    // External Meteorology Data (e.g., ECMWF/VMF, MERRA, GPS,...)
    // -> based on ivs_name   
    if((bool)setup.exists("troposphere"))
    { 
        // use of external meteorological data or zenith hydrostatic delays (ZHDs);
        // in case of GPT2 the meteorological data or the ZHDs are directly calculated
        // in case of raytraced delays, the external data are loaded immediatly 
        // after the TRF initialization (since there is a need of the session name)
        if( (bool)setup["troposphere"]["external_meteo_data"][0]  && 
            string((const char*)setup["troposphere"]["external_meteo_data"][2]) != "gpt2" && 
            string((const char*)setup["troposphere"]["external_meteo_data"][2]) != "raytracing" )
        {
            ivg::parser::external_met_data(this, 
                                           //definitions["troposphere"]["external_meteo_data"][string((const char*)setup["troposphere"]["external_meteo_data"][2])], 
                                           definitions["troposphere"]["external_meteo_data"].lookup(string((const char*)setup["troposphere"]["external_meteo_data"][2])), 
                                           start, end, ivg::extdata::HYDROSTATIC);        
        }
    }
    //Hydrology Loading
    // -> based on ivs_name / only read files if HYDROLOGY LOADING set true in config file
    if((bool)get_list_element(setup["stadisp"], "HYDROLOGY LOADING" )[1])
        ivg::parser::hydlo(this,definitions["stadisp"]["hydlo"][(const char *)get_list_element(setup["stadisp"], "HYDROLOGY LOADING")[3]]);
    
    // in case of itrf2014 and if PSD is set to true, parse the psd_coefficients
    if( _name == "itrf2014" && (bool)get_list_element(setup["stadisp"], "PSD" )[1] )
    {
        log<WARNING>("!!! USING ITRF2014P - Correction of post-seismic deformation - SPECIAL itrf2014 CASE ");
        ivg::parser::psd_coefficients(this,(const char *)get_list_element(setup["stadisp"], "PSD" )[2]);
    }
    
#ifdef DEBUG_REFFRAME
   cerr << "--- void Trf::init_displacements(Setting &setup, ivg::Date start, ivg::Date end)" << " : " << tim.toc() << " s " << endl;
#endif    
}
// ...........................................................................
void Trf::log_data_info_table()
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ void Trf::log_data_info_table()" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // generate output-table
    vector<string> data_types = {"XYZ","OPTL","OLC","VMF1","NTAPL","TAPL","ECC","ANT_I","HYDLO","ANT_C","SEFD","FREQ","MASK","CHA_I"};
    stringstream ssds; // stringstream_data_status
    ssds << "*** TRF " << _name << " Data Information Table (W: Warning, E: Error, X: Set, Empty: Not Set)" << endl;
    ssds << "***          ";
    for(auto &type: data_types)
        ssds << setw(5) << type << "|";
    ssds << endl;
    
    for(auto &sta_tmp: _stations)
    {
        ssds << "*** " << setw(8) << setfill(' ') << sta_tmp.get_name(staname::ivs_name) << "|";
        for(auto &type: data_types)
            ssds << setw(5) << sta_tmp.get_data_status(type) << "|";
        
        ssds << endl;
    }
    ssds << "***          ";
    for(auto &type: data_types)
        ssds << setw(5) << type << "|";
    
    log<INFO>(ssds.str());
    
#if DEBUG_REFFRAME >=3
    cerr << "--- void Trf::log_data_info_table()" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Trf::keep_stations(vector<string> stations, ivg::staname type)
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ void Trf::keep_stations(vector<string> stations, ivg::staname type)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    for(vector<ivg::Analysis_station>::iterator sta = _stations.begin(); sta != _stations.end(); ++sta)
    {
        if(find( stations.begin(), stations.end(), sta->get_name(type)) == stations.end())
        {
            _stations.erase(sta);
            sta--;
        }
    }
   
#if DEBUG_REFFRAME >=3
    cerr << "--- void Trf::keep_stations(vector<string> stations, ivg::staname type)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Trf::remove_station( int idx )
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ void Trf::remove_station( int idx )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    _stations.erase( _stations.begin()+idx );
    
#if DEBUG_REFFRAME >=3
    cerr << "--- void Trf::remove_station( int idx )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Trf::remove_station( vector<ivg::Analysis_station>::iterator remove)
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ void Trf::remove_station( int idx )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    _stations.erase( remove );
    
#if DEBUG_REFFRAME >=3
    cerr << "--- void Trf::remove_station( int idx )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
bool Trf::get_station(ivg::Analysis_station **station, string name, ivg::staname maximal)
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ bool Trf::get_station(vector<ivg::Analysis_station> , ivg::Analysis_station , string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    return(_get_station(_stations, station, name, ivg::staname::MINSTA, maximal));
    
#if DEBUG_REFFRAME >=3
    cerr << "--- bool Trf::get_station(vector<ivg::Analysis_station> , ivg::Analysis_station , string )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
bool Trf::_get_station(vector<ivg::Analysis_station> &stations,
                       ivg::Analysis_station **station, string name, ivg::staname minimal, ivg::staname maximal)
// ...........................................................................
{

    if(name != "")
    {
        vector<ivg::Analysis_station>::iterator iter;
        for(iter=stations.begin(); iter < stations.end(); iter++)
        {
            for ( int i = minimal; i<maximal; i++ )
            {
                if ((*iter).get_name((ivg::staname) i) == name)
                {
                    *station = &(*iter);
                    return(true);
                }
            }
        }
    }
   
    return(false);
}

// ...........................................................................
vector<string>Trf::get_station_names(ivg::staname type) 
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ vector<string>Trf::get_station_names(ivg::staname)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    vector<string> names;
    vector<ivg::Analysis_station>::iterator iter;
    for(iter=_stations.begin(); iter < _stations.end(); iter++)
        names.push_back((*iter).get_name(type));

    
#if DEBUG_REFFRAME >=3
    cerr << "--- vector<string>Trf::get_station_names(ivg::staname)" << " : " << tim.toc() << " s " << endl;
#endif 
    return(names);
}
// ...........................................................................
ivg::Matrix Trf::get_corresponding_stations(ivg::Trf &other)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ ivg::Matrix Trf::get_corresponding_stations(ivg::Trf &other)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ivg::Matrix indexes;
    
    // definition which types of names should be used for the comparison
    vector<ivg::staname> compare = {ivg::staname::cdp, ivg::staname::domes_no, ivg::staname::ivs_name};
    
    // we need to detect if a source is equal to another source
    int this_index=0;
    for(auto &ac_this: (*this))
    {
        ivg::Matrix row(1,2,0.0);
        int other_index=0;
        for(auto &ac_other: other)
        {
            int cmp_cnt=0;
            for ( auto &cmp_name: compare )
            {
                if (ac_this.get_name(cmp_name) == ac_other.get_name(cmp_name))
                    cmp_cnt++;
            }      
        
            if(cmp_cnt == compare.size())
            {
                row(0,0) = this_index;
                row(0,1) = other_index;
                indexes.append_rows(row);
                log<DETAIL>("*** Station ") % ac_this.get_name(ivg::staname::ivs_name) % " successfully matched.";
                break;
            }
            else if(cmp_cnt > 0 && cmp_cnt < compare.size())
                throw runtime_error( "ivg::Matrix Trf::get_corresponding_stations(ivg::Trf &other): Not all stanames correspond!" );
            else
            {
                row(0,0) = 0.0;
                row(0,1) = 0.0;
            }
            
        other_index++;
        }
        
    this_index++; 
    }
    return indexes;
   
#if DEBUG_REFFRAME >=2
    cerr << "--- ivg::Matrix Trf::get_corresponding_stations(ivg::Trf &other)" << " : " << tim.toc() << " s " << endl;
#endif  
}
// ...........................................................................
double Trf::calculate_network_volume()
// ...........................................................................
{
    double V = 0.0;
    
    if(_stations.size() > 3)
    {
        // store coordinates of all stations in correct datatype
        std::vector<K::Point_3> points;
        double mean_x=0.0,mean_y=0.0,mean_z=0.0;
        for(auto &sta: _stations)
        {
            // in order to have slightly different network volumes for 
            // correct tooltip assignment we use the _ref_epoch
            ivg::Matrix xyz = sta.calc_xyz(_ref_epoch);// sta.get_xyz0();
            mean_x += xyz(0,0);
            mean_y += xyz(1,0);
            mean_z += xyz(2,0);
            points.push_back(K::Point_3(xyz(0,0),xyz(1,0),xyz(2,0)));
        }

        // calculate center of all stations as centerpeak of the pyramids
        mean_x /= _stations.size();
        mean_y /= _stations.size();
        mean_z /= _stations.size();

        // define polyhedron to hold convex hull
        CGAL::Polyhedron_3<K>  poly;

        // compute convex hull of non-collinear points
        CGAL::convex_hull_3(points.begin(), points.end(), poly); 

        // assign a plane equation to each polyhedron facet using functor Plane_from_facet
        std::transform( poly.facets_begin(), poly.facets_end(), poly.planes_begin(),Plane_from_facet() );

        // calculate volume for each pyramid, using the base area A and the height d to the center
        for(CGAL::Polyhedron_3<K>::Facet_iterator it=poly.facets_begin(); it != poly.facets_end(); it++)
        { 
            double A = K::Compute_area_3()( (*it).halfedge()->vertex()->point(), (*it).halfedge()->next()->vertex()->point(), (*it).halfedge()->opposite()->vertex()->point() );
            double d = std::sqrt(CGAL::squared_distance((*it).plane(),K::Point_3(mean_x, mean_y, mean_z)));
            V += d * A / 3.0;
        } 
    }
    return(V);
}
// ...........................................................................
void Trf::show(bool verbose)
// ...........................................................................
{
    cout << "------------------------------------------------" << endl;
    if(verbose)
    {
        for(int i=0; i<_stations.size(); i++)
            _stations.at(i).show();
    }

    vector<string> station_names = get_station_names(ivg::staname::ivs_name);
    cout << " Following " << station_names.size() << " stations are included based on " << _name << " with reference epoch at " << _ref_epoch.get_date_time("YYYY/MON/DD") << endl << " ";
    copy(station_names.begin(), station_names.end(),
         ostream_iterator<string>(cout, " | "));
    cout << endl;
    cout << "------------------------------------------------" << endl;
}

// ...........................................................................
void Trf::create_station_indices()
// ...........................................................................
{
    for(int i = 0; i < _stations.size(); ++i){
        _stations[i].set_idx(i);
    }
}

// ...........................................................................
void Trf::refresh_station_twin_list(){
// ...........................................................................
    _twin_stations.clear();
    std::vector< std::vector<ivg::Analysis_station*> > twin;
    twin.resize(_twin_station_names.size());
    
    // loop over each station in TRF
    for( ivg::Analysis_station & sta : _stations ){
        
        // check if station is in twin station
        // if so, add station to temporary vector
        bool isTwin = false;
        for( unsigned i = 0; i < _twin_station_names.size(); ++i ){
            std::string sta_name = sta.get_name(ivg::staname::ivs_name);
            if(_twin_station_names[i].find( sta_name ) != _twin_station_names[i].end()){
                isTwin = true;
                twin[i].push_back(&sta);
                break;
            }
        }
        
        // if station is not in  twin list add it directly
        if(!isTwin){
            std::vector<ivg::Analysis_station*>  tmp = {&sta};
            _twin_stations.push_back( tmp );
        }

    }
    
    // Now add found twin stations to the result
    for( std::vector<ivg::Analysis_station*>& t : twin ){
        if(t.size() > 0){
            _twin_stations.push_back(t);
        }
    }

}

// ...........................................................................
void Trf::remove_from_station_twin_list(std::string station){
// ...........................................................................
    int del(-1), i(0);
    for(std::set<std::string>& twins: _twin_station_names){
        std::set<std::string>::iterator it = twins.find( station );
        if( it != twins.end()){
            twins.erase(it);
            if(twins.size() < 2){
                del = i;
                break;
            }
        }
        ++i;
    }
    
    if(del>=0){
        _twin_station_names.erase(_twin_station_names.begin()+del);
    }
}

// ...........................................................................
bool Trf::areTwins(unsigned sta1, unsigned sta2){
// ...........................................................................
    bool isTwin = false;
    std::string sta_name1 = _stations.at(sta1).get_name(ivg::staname::ivs_name);
    std::string sta_name2 = _stations.at(sta2).get_name(ivg::staname::ivs_name);
    
    for( unsigned i = 0; i < _twin_station_names.size(); ++i ){
        if(_twin_station_names[i].find( sta_name1 ) != _twin_station_names[i].end() &&
           _twin_station_names[i].find( sta_name2 ) != _twin_station_names[i].end()    ){
            return true;
        }
    }
    
    return false;
}

// ...........................................................................
std::vector<ivg::Analysis_station*> Trf::get_twins(ivg::Analysis_station* sta){
// ...........................................................................
    for(std::vector<ivg::Analysis_station*>& stas: this->_twin_stations){
        if( std::find(stas.begin(), stas.end(), sta) != stas.end()){
            return stas;
        }

    }  
    
    std::vector<ivg::Analysis_station*> tmp;
    return tmp;
}

// ...........................................................................
std::map<unsigned, std::vector<ivg::Analysis_station*> > Trf::get_twins_map(){
// ...........................................................................
    std::map<unsigned, std::vector<ivg::Analysis_station*> > twin_map;
    for( std::vector<ivg::Analysis_station*>& stas: this->_twin_stations ){
        for(ivg::Analysis_station* sta: stas){
            twin_map[sta->get_idx()] = stas;
        }
    } 
    
    return twin_map;
}

// ...........................................................................
std::map<std::string, std::vector<ivg::Analysis_station*> > Trf::get_twins_map(ivg::staname type) {
// ...........................................................................
    std::map<std::string, std::vector<ivg::Analysis_station*> > twin_map;
    for( std::vector<ivg::Analysis_station*>& stas: this->_twin_stations ){
        for(ivg::Analysis_station* sta: stas){
            std::string name = sta->get_name(type);
            twin_map[name] = stas;
        }
    } 
    
    return twin_map;
}

std::map<string, std::set<std::string> > Trf::get_twins_map_including_all() {
    std::map<std::string, std::set<std::string> > twin_map;
    for( const std::set<std::string>& twins: _twin_station_names){
        for( const std::string& sta : twins){
            twin_map[sta] = twins;
        }
    }
    
    std::map<std::string, std::vector<ivg::Analysis_station*> > twinsInInstance = get_twins_map(ivg::staname::ivs_name);
    for(auto a : twinsInInstance){
        if( !twin_map.count( a.first ) ){
            std::set<std::string> s;
            for( ivg::Analysis_station* sta : a.second )
                s.insert(sta->get_name(ivg::staname::ivs_name));
            twin_map[a.first] = s;
        }
    }
    
    
    return twin_map;
}

}
