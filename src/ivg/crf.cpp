#include <string>

#include "crf.h"
#include "scan.h"

namespace ivg
{
// ..........................................................................
Crf::Crf()
// ..........................................................................
{
}
// ...........................................................................
Crf::Crf( string name, vector<ivg::Source> sources ) 
{
    _name = name;
    _sources= sources;
}
// ..........................................................................
Crf::Crf( Setting &setup, vector<ivg::Source> sources, void **ephem, std::string skdfile )
// ..........................................................................
{
#ifdef DEBUG_REFFRAME
   cerr << "+++ Crf::Crf( Setting &, vector<ivg::Source> )" << endl; 
   tictoc tim;
   tim.tic();
#endif

    _name = (const char *)setup["crf"];
    _sources = sources;
   
    if( skdfile.empty() )
        log<INFO>("*** Initializing ivg::Crf with ") % sources.size() % " sources based on " % _name;
    else
        log<INFO>("*** Initializing ivg::Crf with ") % sources.size() % " sources based on " % skdfile;
    
    if(_name == "ocars")
    {
        _ocars_parser(_sources, (const char *)get_list_element(setup["refframes"],_name)[2], true);
    }
    else if(_name == "ivssrc")
    {
      
        _ivssrc_parser(_sources, (const char *)get_list_element(setup["refframes"],_name)[2], true);
    }
    else if(_name == "icrf2")
    {
        _ivssrc_parser(_sources, (const char *)get_list_element(setup["refframes"],"ivssrc")[2], false);
        _icrf2_parser(_sources, (const char *)get_list_element(setup["refframes"],_name)[3]);
        _icrf2_parser(_sources, (const char *)get_list_element(setup["refframes"],_name)[2]);
    }
    else if(_name == "icrf3")
    {
        _ivssrc_parser(_sources, (const char *)get_list_element(setup["refframes"],"ivssrc")[2], false);
        _icrf3_parser(_sources, (const char *)get_list_element(setup["refframes"],_name)[2]);
    }
    else if(_name == "ngs" || _name == "snx" || _name == "vgosdb")
    {
        _ivssrc_parser(_sources, (const char *)get_list_element(setup["refframes"],"ivssrc")[2],false);
    }
    ivg::Source* tmpsou;
    
    bool sats_init_once = false;
    for(vector<ivg::Source>::iterator iter=_sources.begin(); iter < _sources.end(); iter++)
    {
        // in case of initialized ephemeris, set the pointer for each source
        // necessary for calculation of sun-position for sun-quasar-arcdistance
        if(ephem != NULL)
            (*iter).set_jpl_ephem(ephem);
	
        // if source has not been found yet
        if(!(*iter).get_found())
        {
            string src_notfound = (*iter).get_name();
            // if a source couldn't be initialized it might be satellite or moon
            if(src_notfound.substr(0,2) == "PG")
            {
                if(sats_init_once == false)
                { 
                    // e.g. "sp3"
                    string sat_ephem_type = get_list_element(setup["NEARFIELD"]["ephemerides"],(const char *)setup["NEARFIELD"]["sat_ephem"])[1];
                    // e.g. "/home/ascot/apriori_files/sat_catalogs/15AUG24X.sp3"
                    string sat_ephem_path = get_list_element(setup["NEARFIELD"]["ephemerides"],(const char *)setup["NEARFIELD"]["sat_ephem"])[2];
                    if(sat_ephem_type == "tle")
                        // initialize satellite orbits based on TLEs
                        ivg::parser::tle(this, sat_ephem_path);
                    else if(sat_ephem_type == "sp3")
                        // initialize satellite orbits based on IGS finals
                        // (in case of session 15AUG24X, manually merged igs18590.sp3.Z  igs18591.sp3.Z  igs18592.sp3.Z 
                        // from ftp://cddis.gsfc.nasa.gov/gps/products/1859/ to 15AUG24X.sp3) 
                        ivg::parser::sp3(this, sat_ephem_path);
                    else
                        throw runtime_error( "Crf::Crf( Setting &setup, vector<ivg::Source> sources ): Unknown sat_ephem_type: "+sat_ephem_type );
                    
                    // load ephems only once in a session
                    sats_init_once = true;
                }
                
                (*iter).set_type(ivg::srctype::satellite);
                (*iter).set_name(ivg::srcname::ivs, src_notfound);
                (*iter).set_name(ivg::srcname::iers, src_notfound);
                log<INFO>("*** CRF Initialization: Source ") % (*iter).get_name() % " initialized as satellite.";
            }
            // need to check for correct moonname 
            // up to now, a specific defintion of the moonname is not defined. Might be epoch-wise?
            else if(src_notfound== "MOONNAME ?!?!?")
            {
                (*iter).set_type(ivg::srctype::moon);
                
                // NO FURTHER IMPLEMENTATION YET
                // IT'S ONLY POSSIBLE TO CALCULATE ROVER POSITION BASED ON LOADED JPL-EPHEMERIDES
                // NO TLE-SUPPORT YET
                // Rover Position [x,y,z]' w.r.t. SSB at epoch
                // e.g. call ivg::Matrix pos_SSB = (*iter).calc_rover_position( ivg::Date(2014, 04, 11, 18, 3, 31) );
                
                log<INFO>("*** CRF Initialization: Source ") % (*iter).get_name() % " initialized as moon.";
            }
            else
            {
                if(_name == "ngs")
                {
                    cerr << "!!! CRF Initialization: Source " << (*iter).get_name()
                    << " not found in database. Positions taken from NGS File" << endl;
                    
                     (*iter).set_name(ivg::srcname::ivs, (*iter).get_name(ivg::srcname::ngs));
                     (*iter).set_name(ivg::srcname::iers, (*iter).get_name(ivg::srcname::ngs));
                }
                else
                {
		  if (iter->get_name(ivg::srcname::iers).empty())
		    {
		      cerr << "!!! CRF Initialization: Source " << (*iter).get_name()
                         << " not found in selected apriori file or in IVS source table (" << _name << "). Setting apriori value to zero" <<  endl;
		      iter->set_ra0(0.0);
		      iter->set_dec0(0.0);
		      iter->set_name(ivg::srcname::iers,iter->get_name(ivg::srcname::ivs));
		    } else
		    {
                    cerr << "!!! CRF Initialization: Source " << (*iter).get_name()
                         << " not found in selected apriori file (" << _name << "). Setting apriori value to IVS source table" <<  endl;
		    // iter->set_ra0(1.0);
                    //iter->set_dec0(1.0);
		    iter->set_ra_dec_bck();
                    iter->set_found(true);
//                    _sources.erase(iter);
//                    iter--;
		    }
                }
            }
        }
        // if source was already found in the apriori files it must be a source
        else
                (*iter).set_type(ivg::srctype::source);
	
    }
    
    if( (setup.exists( "SIM" ) && (bool)setup[ "SIM" ]["apply"]) || (setup.exists( "SKED" ) && (bool)setup[ "SKED" ]["apply"])
        || ( setup.exists( "SKED" ) &&  (bool)setup[ "SKED" ]["createPlots"] && ( !setup.exists( "SIM" ) || !(bool)setup[ "SIM" ]["apply"]))  )
    {
        bool init_from_sked = (setup)["SKED"].exists( "init_from_sked" ) &&  (bool)(setup)["SKED"]["init_from_sked"]["apply"] &&  (bool)(setup)["SKED"]["init_from_sked"]["flux"] ;
        if(init_from_sked){
            log<INFO>("*** Initializing ivg::Crf with flux catalog from sked file.");
                ivg::parser::flux_cat(this, skdfile, true);   
        } else {
            log<INFO>("*** Initializing ivg::Crf with flux catalog.");
            ivg::parser::flux_cat(this, setup["definitions"]["sked"]["flux"], false);   
        }
    }
    
    if((bool)setup.exists("vascc2015") && (bool)setup["vascc2015"] == true)
    {
        // 		RA [rad]              DE[rad] (from VASCC2015.pdf)
        // 0039+568    +1.84674138320125E-01 +9.97342153613740E-01
        // 0230-790    +6.52676534369252E-01 -1.37524964611452E+00
        ivg::Source *src;
        if(get_source(&src, "0039+568"))
        {
            log<WARNING>("!!! SPECIAL VASCC2015 CASE: Setting specific VASCC2015 0039+568 source position in crf");
            src->set_ra0(+1.84674138320125E-01);
            src->set_dec0(+9.97342153613740E-01);
        }
        if(get_source(&src, "0230-790"))
        {
            log<WARNING>("!!! SPECIAL VASCC2015 CASE: Setting specific VASCC2015 0230-790 source position in crf");
            src->set_ra0(+6.52676534369252E-01);
            src->set_dec0(-1.37524964611452E+00);
        }
    }
    
#ifdef DEBUG_REFFRAME
   cerr << "--- Crf::Crf( Setting &, vector<ivg::Source> )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Crf::remove_source( std::vector<ivg::Source>::iterator it )
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ void Crf::remove_source( std::vector<ivg::Source>::iterator )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    _sources.erase( it );
    
#if DEBUG_REFFRAME >=3
   cerr << "--- void Crf::remove_source( std::vector<ivg::Source>::iterator )" << " : " << tim.toc() << " s " << endl;
   tim.tic();
#endif
}
// ...........................................................................
void Crf::remove_source( int idx )
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ void Crf::remove_source( int idx )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    _sources.erase( _sources.begin()+idx );
    
#if DEBUG_REFFRAME >=3
    cerr << "--- void Crf::remove_source( int idx )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
bool Crf::get_source(ivg::Source **source, string name)
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ bool Crf::get_source(ivg::Source ** , string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    return(_get_source(_sources, source, name));
    
#if DEBUG_REFFRAME >=3
    cerr << "--- bool Crf::get_source(ivg::Source ** , string )" << " : " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
bool Crf::_get_source(vector<ivg::Source> &sources, ivg::Source **source, string name)
// ...........................................................................
{
    if(name != "")
    {
        vector<ivg::Source>::iterator iter;
        for(iter=sources.begin(); iter < sources.end(); iter++)
        {
            for ( int i = 0; i<MAXSRC; i++ )
            {
                if ((*iter).get_name((ivg::srcname) i) == name)
                {
                    *source = &(*iter);
                    return(true);
                }
            }
        }
    }
    return(false);
}
// ...........................................................................
ivg::Matrix Crf::get_corresponding_sources(ivg::Crf &other)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ ivg::Matrix Crf::get_corresponding_sources(ivg::Crf &other)" << endl; 
   tictoc tim;
   tim.tic();
#endif
      
    ivg::Matrix indexes;
    
    // definition which types of names should be used for the comparison
    vector<ivg::srcname> compare = {ivg::srcname::iers, ivg::srcname::ivs};
    
    // we need to detect if a source is equal to another source
    int this_index=0;
    for(auto &src_this: (*this))
    {
        ivg::Matrix row(1,2,0.0);
        int other_index=0;
        for(auto &src_other: other)
        {
            int cmp_cnt=0;
            for ( auto &cmp_name: compare )
            {
                if (src_this.get_name(cmp_name) == src_other.get_name(cmp_name))
                    cmp_cnt++;
            }      
        
            // only if a corresponding name is found
            if(cmp_cnt > 0 && cmp_cnt <= compare.size())
            {
                row(0,0) = this_index;
                row(0,1) = other_index;
                indexes.append_rows(row);
                log<DETAIL>("*** Source ") % src_this.get_name(ivg::srcname::iers) % " successfully matched (" % cmp_cnt % "). [" 
                            % indexes.rows() % "/" % (*this)._sources.size() % "]";
                break;
            }
            else
            {
                row(0,0) = -1.0;
                row(0,1) = -1.0;
            }
            
        other_index++;
        }
        
    this_index++; 
    }
    return indexes;
   
#if DEBUG_REFFRAME >=2
    cerr << "--- ivg::Matrix ivg::Matrix Crf::get_corresponding_sources(ivg::Crf &other)" << " : " << tim.toc() << " s " << endl;
#endif  
}
// ...........................................................................
vector<string> Crf::get_source_names(ivg::srcname type)
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ vector<string> Crf::get_source_names(ivg::srcname)" << endl; 
   tictoc tim;
   tim.tic();
#endif
    vector<string> names;
    vector<ivg::Source>::iterator iter;
    for(iter=_sources.begin(); iter < _sources.end(); iter++)
        if(iter->use_me() == true)
        names.push_back((*iter).get_name(type));

#if DEBUG_REFFRAME >=3
    cerr << "--- vector<string> Crf::get_source_names(ivg::srcname)" << " : " << tim.toc() << " s " << endl;
#endif
    return(names);
}
// ...........................................................................
int Crf::get_number_sources()
// ...........................................................................
{
    int cnt=0;
    for(auto &src: _sources)
        if(src.use_me() == true)
            cnt++;

    return cnt;
}


// ...........................................................................
void Crf::create_source_indices()
// ...........................................................................
{
    for(int i = 0; i < _sources.size(); ++i){
        _sources[i].set_idx(i);
    }
}

// ...........................................................................
void Crf::_ocars_parser(vector<ivg::Source> &sources, const string path,
                        bool set_pos)
// ..........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void Crf::_ocars_parser(vector<ivg::Source> & , const string , bool )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    vector<ivg::Source> crf;

    ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "void Crf::_ocars_parser(vector<ivg::Source> &sources, const string path,bool set_pos): Failed to open file: "+path );
    }else{

        string line;
        while (getline(inStream,line,'\n'))
        {
            if(line.substr(0,1) != "#")
            {
                string ivs = remove_spaces_end(line.substr(0,8));
                string iers = remove_spaces_end(line.substr(9,8));

                ivg::Source * tmp_source;
                if(_get_source(sources, &tmp_source, iers) || _get_source(sources, &tmp_source, ivs) )
                {
                    tmp_source->set_found(true);

                    // only set iers name if it isn't already defined
                    if(tmp_source->get_name(ivg::srcname::iers).empty())
                        tmp_source->set_name(ivg::srcname::iers, iers);
                    else
                    {
                        // check if already defined name is equal to iers name within ocars catalog
                        if(tmp_source->get_name(ivg::srcname::iers) != iers )
                        {
                            log<WARNING>("!!! ivg::Crf init: IERS name ") % tmp_source->get_name(ivg::srcname::iers) % " in SINEX != " % iers % " in ocars.";
//                            tmp_source->set_name(ivg::srcname::iers, iers);
                        }
                    }

                    if(ivs.size() == 0)
                        tmp_source->set_name(ivg::srcname::ivs, tmp_source->get_name(ivg::srcname::iers));
                    else
                    {
                        if(tmp_source->get_name(ivg::srcname::ivs).empty())
                            tmp_source->set_name(ivg::srcname::ivs, ivs);
                    }


                    //position in right ascension and declination
                    if(set_pos)
                    {
                        tmp_source->set_ra0(atoi(line.substr(19,2).c_str()), atoi(line.substr(22,2).c_str()), s2d(line.substr(25,7)));

                        bool negative = false;
                        if(line.substr(34,1) == "-")
                            negative = true;

                        tmp_source->set_dec0(negative, atoi(line.substr(35,2).c_str()), atoi(line.substr(38, 2).c_str()), s2d(line.substr(41,6)));
                    }

                    // if ICRF name is available in ocars
                    if(line.substr(71,4) == "ICRF")
                    {
                        string icrf = remove_spaces_end(line.substr(76,16));
                        if(tmp_source->get_name(ivg::srcname::icrf).empty())
                            tmp_source->set_name(ivg::srcname::icrf, icrf);
                        else
                        {
                            // check if already defined name is equal to icrf name within ocars
                            if(tmp_source->get_name(ivg::srcname::icrf) != icrf )
                            {
                                log<WARNING>("!!! ivg::Crf init: ICRF name ") % tmp_source->get_name(ivg::srcname::iers) % " in SINEX != " % iers % " in ocars. Using ocars.";
                                tmp_source->set_name(ivg::srcname::icrf, icrf);
                            }
                        }
                    }

                    if(line.size() > 102)
                    {
                        if(line.substr(94,8) == "VCS-only")
                            tmp_source->set_vcs(true);
                    }
                }
            }
        }
    }
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void Crf::_ocars_parser(vector<ivg::Source> & , const string , bool )" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
void Crf::_ivssrc_parser(vector<ivg::Source> &sources, const string path,
                        bool set_pos)
// ..........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void Crf::_ivssrc_parser(vector<ivg::Source> & , const string , bool )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    vector<ivg::Source> crf;
    
    ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "void Crf::_ivssrc_parser(vector<ivg::Source> &sources, const string path,bool set_pos): Failed to open file: "+path );
    }else{

        string line;
        while (getline(inStream,line,'\n'))
        {
            if(line.substr(0,1) != "#")
            {
                string ivs = remove_spaces_end(line.substr(0,8));
                string iers = remove_spaces_end(line.substr(40,8));
		if ((iers=="")||(iers=="-"))
		  iers=ivs;
                ivg::Source * tmp_source;
                if(_get_source(sources, &tmp_source, iers) || _get_source(sources, &tmp_source, ivs) )
                {
                    

                    // only set iers name if it isn't already defined
                    if(tmp_source->get_name(ivg::srcname::iers).empty())
		      {
			tmp_source->set_name(ivg::srcname::iers, iers);
		      }
                    else
                    {
                        // check if already defined name is equal to iers name within ocars catalog
                        if(tmp_source->get_name(ivg::srcname::iers) != iers )
                        {
                            log<WARNING>("!!! ivg::Crf init: IERS name ") % tmp_source->get_name(ivg::srcname::iers) % " in SINEX != " % iers % " in IVS_SrcNamesTable.";
//                            tmp_source->set_name(ivg::srcname::iers, iers);
                        }
                    }

		    if(tmp_source->get_name(ivg::srcname::ivs).empty())
		      tmp_source->set_name(ivg::srcname::ivs, ivs);
 


                    //position in right ascension and declination
                    if(set_pos)
                    {
                        tmp_source->set_found(true);
		        tmp_source->set_ra0(atoi(line.substr(64,2).c_str()), atoi(line.substr(67,2).c_str()), s2d(line.substr(70,9)));
			
                        bool negative = false;
                        if(line.substr(81,1) == "-")
                            negative = true;

                        tmp_source->set_dec0(negative, atoi(line.substr(82,2).c_str()), atoi(line.substr(85, 2).c_str()), s2d(line.substr(88,8)));
			tmp_source->set_bck_ra0(atoi(line.substr(64,2).c_str()), atoi(line.substr(67,2).c_str()), s2d(line.substr(70,9)));
			tmp_source->set_bck_dec0(negative, atoi(line.substr(82,2).c_str()), atoi(line.substr(85, 2).c_str()), s2d(line.substr(88,8)));
                    } else
		    {
                        tmp_source->set_bck_ra0(atoi(line.substr(64,2).c_str()), atoi(line.substr(67,2).c_str()), s2d(line.substr(70,9)));
			
                        bool negative = false;
                        if(line.substr(81,1) == "-")
                            negative = true;

                        tmp_source->set_bck_dec0(negative, atoi(line.substr(82,2).c_str()), atoi(line.substr(85, 2).c_str()), s2d(line.substr(88,8)));
		    }

                        string icrf = remove_spaces_end(line.substr(10,16));
                        if(tmp_source->get_name(ivg::srcname::icrf).empty())
                            tmp_source->set_name(ivg::srcname::icrf, icrf);
                        else
                        {
                            // check if already defined name is equal to icrf name within ocars
                            if(tmp_source->get_name(ivg::srcname::icrf) != icrf )
                            {
                                log<WARNING>("!!! ivg::Crf init: ICRF name ") % tmp_source->get_name(ivg::srcname::iers) % " in SINEX != " % iers % " in IVS_SrcNamesTable. Using IVS_SrcNamesTable.";
                                tmp_source->set_name(ivg::srcname::icrf, icrf);
                            }
                        }
                    

                    
                }
            }
        }
    }
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void Crf::_ivssrc_parser(vector<ivg::Source> & , const string , bool )" << " : " << tim.toc() << " s " << endl;
#endif 
}


// ...........................................................................
void Crf::_icrf2_parser(vector<ivg::Source> &sources,
                        const string path_sources)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void Crf::_icrf2_parser(vector<ivg::Source> & , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;

    int found = 0;
    while (ivg::parser::get_line(path_sources,inStream, line))
    {
            // names in apriori file 
            string icrf = remove_spaces_end(line.substr(5,16));
            string iers = remove_spaces_end(line.substr(23,8));

            ivg::Source * tmp_source_icrf, * tmp_source_iers;
            bool icrf_found = _get_source(sources, &tmp_source_icrf, icrf);
            bool iers_found = _get_source(sources, &tmp_source_iers, iers);

            if( (icrf_found||iers_found)&( tmp_source_icrf == tmp_source_iers) || (icrf_found == true && iers_found == false) )
            {
                found++;

                string is_vcsfile = line.substr(33,1);

                int offset=0;
                if(is_vcsfile == "D" || is_vcsfile == " ")
                {
                    offset = 3;
                    if(is_vcsfile == "D")
                        tmp_source_icrf->set_defining(true);
                }

                //position in right ascension and declination
                tmp_source_icrf->set_ra0(atoi(line.substr(33+offset,2).c_str()),atoi(line.substr(36+offset,2).c_str()), s2d(line.substr(39+offset,11)));

                bool negative = false;
                if(line.substr(52+offset,1) == "-")
                    negative = true;
		
                tmp_source_icrf->set_dec0(negative, atoi(line.substr(53+offset,2).c_str()),atoi(line.substr(56+offset,2).c_str()), s2d(line.substr(59+offset,10)));
                tmp_source_icrf->set_found(true);
               
                // set formal errors of source position
                tmp_source_icrf->set_sigma_ra0dec0(s2d(line.substr(71+offset,10)), s2d(line.substr(82+offset,9)));
                tmp_source_icrf->set_additional_information(ivg::Date(s2d(line.substr(109+offset,7))),ivg::Date(s2d(line.substr(101+offset,7)))
                                                            ,ivg::Date(s2d(line.substr(117+offset,7))), atoi(line.substr(126+offset,5).c_str()),
                                                            atoi(line.substr(132+offset,6).c_str()));
            }
        }
   
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void Crf::_icrf2_parser(vector<ivg::Source> & , const string )" << " : " << tim.toc() << " s " << endl;
#endif 
}

void Crf::_icrf3_parser(vector<ivg::Source> &sources,
                        const string path_sources)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void Crf::_icrf2_parser(vector<ivg::Source> & , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;

    int found = 0;
    while (ivg::parser::get_line(path_sources,inStream, line))
    {
            // names in apriori file 
            string icrf = remove_spaces_end(line.substr(5,16));
            string iers = remove_spaces_end(line.substr(25,8));

            ivg::Source * tmp_source_icrf, * tmp_source_iers;
            bool icrf_found = _get_source(sources, &tmp_source_icrf, icrf);
            bool iers_found = _get_source(sources, &tmp_source_iers, iers);

            if( (icrf_found||iers_found)&(( tmp_source_icrf == tmp_source_iers) || (icrf_found == true && iers_found == false) ))
            {
                found++;
                string is_defining = line.substr(35,1);

                if(is_defining == "D")
                        tmp_source_icrf->set_defining(true);

                //position in right ascension and declination
                tmp_source_icrf->set_ra0(atoi(line.substr(40,2).c_str()),atoi(line.substr(43,2).c_str()), s2d(line.substr(46,11)));

                bool negative = false;
                if(line.substr(61,1) == "-")
                    negative = true;

                tmp_source_icrf->set_dec0(negative, atoi(line.substr(62,2).c_str()),atoi(line.substr(65,2).c_str()), s2d(line.substr(68,10)));
                tmp_source_icrf->set_found(true);
		
                // set formal errors of source position
                tmp_source_icrf->set_sigma_ra0dec0(s2d(line.substr(83,10)), s2d(line.substr(98,9)));
                tmp_source_icrf->set_additional_information(ivg::Date(s2d(line.substr(127,7))),ivg::Date(s2d(line.substr(118,7)))
                                                            ,ivg::Date(s2d(line.substr(136,7))), atoi(line.substr(144,5).c_str()),
                                                            atoi(line.substr(150,6).c_str()));
		
	    }
        }
    
    inStream.close();
   
#if DEBUG_REFFRAME >=2
   cerr << "--- void Crf::_icrf2_parser(vector<ivg::Source> & , const string )" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
void Crf::show(bool full_info)
// ...........................................................................
{
    cout << "-----------------Crf::show------------------" << endl;
    cout << "#sources: " << _sources.size() << " including use_me == false ones" << endl;
    if(full_info)
    {
        for(int i=0; i<_sources.size(); i++)
            _sources.at(i).show();
    }
    vector<string> source_names = get_source_names(ivg::srcname::iers);
    cout << " IERS: Following sources are included based on " << _name << " positions (" << source_names.size() << ")" 
         << endl << " ";
    copy(source_names.begin(), source_names.end(), ostream_iterator<string>(cout,
            " | "));
    cout << endl;
    cout << "------------------------------------------------" << endl;
    source_names = get_source_names(ivg::srcname::ngs);
    cout << " NGS: Following sources are included based on " << _name << " positions (" << source_names.size() << ")" 
         << endl << " ";
    copy(source_names.begin(), source_names.end(), ostream_iterator<string>(cout,
            " | "));
    cout << endl;
    cout << "------------------------------------------------" << endl;
    source_names = get_source_names(ivg::srcname::ivs);
    cout << " IVS: Following sources are included based on " << _name << " positions (" << source_names.size() << ")" 
         << endl << " ";
    copy(source_names.begin(), source_names.end(), ostream_iterator<string>(cout,
            " | "));
    cout << endl;
        cout << "------------------------------------------------" << endl;
    source_names = get_source_names(ivg::srcname::icrf);
    cout << " ICRF: Following sources are included based on " << _name << " positions (" << source_names.size() << ")" 
         << endl << " ";
    copy(source_names.begin(), source_names.end(), ostream_iterator<string>(cout,
            " | "));
    cout << endl;
    cout << "-----------------Crf::show------------------" << endl;
}


} // # namespace ivg
