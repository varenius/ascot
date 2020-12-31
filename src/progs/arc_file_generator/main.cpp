 
#include "date.h"
#include "auxfunc.h"
#include <tclap/CmdLine.h>
#include <curl/curl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <cstdlib>
#include <libconfig.h++>

#include "trf.h"
#include "logger.h"
#include "masterfile.h"
#include "parser.h"

loglevel g_verbose;

int main( int argc, char *argv[] )
{ 
    // .........................................................................
    // handle commandline arguments
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
    
    // name of config file
    TCLAP::ValueArg<std::string> m1Arg( "c","config","Name of the control file",true,"","string");
    cmd.add( m1Arg );
    
    // parse commandline arguments
    cmd.parse( argc, argv );
    
    string control_file = m1Arg.getValue();
    
    std::cout << " \n ============================================ " << std::endl;
    std::cout << "     > > > > > ARC FILE GENERATOR  < < < < <     " << std::endl << std::endl;

    // loading controlfile
    Config cfg;
    cfg.readFile( control_file.c_str() );
    
    // get range between
    ivg::Date start_epoch(cfg.lookup( "start" ), "YYMMMDD");
    ivg::Date end_epoch(cfg.lookup( "end" ), "YYMMMDD");
    
    // get path for output arcfile
    string output_arcfile = cfg.lookup( "output_folder" );
    
    // get type for which the arc file should be generated (snx/vgosdb/ngs)
    string arcfile_type = cfg.lookup( "arcfile_type" );
    
    // get the directory where the masterfiles are stored
    string masterfile_dir = cfg.lookup("masterfiles");
    ivg::Masterfile masterfile(masterfile_dir, ivg::mastertype::both, start_epoch.get_int_year(), end_epoch.get_int_year());
    
    // get range of network volume
    double min_volume = cfg.lookup( "min_volume" );
    double max_volume = cfg.lookup( "max_volume" );
    
    // get minimal and maximal number of stations
    int min_stations = cfg.lookup( "min_stations" );
    int max_stations = cfg.lookup( "max_stations" );
    
    // get groups which should be excluded
    vector<string> excl_groups;
    for( int i=0; i<cfg.lookup("exclude_groups").getLength(); ++i )
        excl_groups.push_back(cfg.lookup("exclude_groups")[i]);
    
    // get database code which should be execluded
    vector<string> excl_db_codes;
    for( int i=0; i<cfg.lookup("exclude_db_codes").getLength(); ++i )
        excl_db_codes.push_back(cfg.lookup("exclude_db_codes")[i]);
    
    // get sessions which should be excluded
    vector<string> excl_sessions;
    for( int i=0; i<cfg.lookup("exclude_sessions").getLength(); ++i )
        excl_sessions.push_back(cfg.lookup("exclude_sessions")[i]);
       
    // get the file containing sessions which should be ignored for the arcfile
    string exclude_file = cfg.lookup("exclude_file");
    if(exclude_file.size() > 0)
    {
        ifstream inStream;
        string line;
        while( ivg::parser::get_line(exclude_file, inStream, line) )
        {
            replace_string_in_place( line , "-", "X" );
            excl_sessions.push_back(remove_spaces_end( line ));
        }
    }
    Setting& versions= cfg.lookup( "versions" );
    
    // in case of several versions, we need to create a directory
    string new_folder;
    new_folder = output_arcfile.substr(0,output_arcfile.find_last_of("/\\"));

    if(mkdir(new_folder.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) { 
        cerr << "*** Folder already exists " << new_folder << endl;
    }
    else
        cerr << "*** Created new folder " << new_folder << endl;
    
    // check whether scaling modus is activated
    bool scaling_modus = false;
    struct block{
        ivg::Date start;
        ivg::Date end;
        int cnt;
    };
    vector<block> block_vector; 
    if(cfg.exists("scaling_modus"))
    {
        // -> days-blocks will be created
        scaling_modus = cfg.lookup( "scaling_modus" );

        if(scaling_modus)
        {
            int block_interval = cfg.lookup( "block_interval" );
            ivg::Date block_start = start_epoch;
            do
            {
                block_vector.push_back({block_start,block_start.add_days(block_interval-(1/86400.0))});
                block_start.add_secs(block_interval*86400.0);
            }
            while(block_start.add_days(block_interval) < end_epoch );

            // finally add last possible block (might be smaller than defined)
            block_vector.push_back({block_start, end_epoch});

            cerr << "*** " << block_vector.size() << " blocks created based on interval of " << block_interval << " days" << endl;
        }
    }
    
    vector<string> arc_files_text, include_text;
    vector< vector<string> >ref_db_names;
    for(int i=0; i<versions.getLength(); i++)
    {
        // get version as first elemente
        string version = (const char *)versions[i][0];
        
        // in case of some add_version
        string add_version = "";
        if(versions[i].getLength() == 2)
            add_version = (const char *)versions[i][1];
        // in case of version == -1 (vgsodb use the full origin outfile)
        //in case of several versions, we need to define the arc file name
        if(add_version == "" && version != "-1")
            output_arcfile = new_folder + "/" + version + ".arc";
        else if(add_version == "" && version == "-1")
            output_arcfile = output_arcfile+"/arcfile.arc";
        else
            output_arcfile = new_folder + "/combi.arc";
        vector<string> paths;
        string data_path = get_list_element(cfg.lookup("datadirs"),arcfile_type)[1];
        if(arcfile_type == "snx")
        {
            data_path = data_path+"/"+version+"/";
            paths.push_back(data_path);
        }
        else if(arcfile_type == "ngs" || arcfile_type == "vgosdb")
        {
            for(int year = start_epoch.get_int_year(); year<= end_epoch.get_int_year(); year++)
                paths.push_back(data_path+to_string(year));
        }

        ////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////
        map<string,string> sources;
        if(arcfile_type == "snx")
        {
            bool use_super_file = cfg.lookup("use_super_file");
            if(use_super_file)
            {
                // import super_file
                string super_file = cfg.lookup("super_file");

                Setting &handling = cfg.lookup("super_file_handling");

                ofstream outstream_hlf ( "/home/iddink/handling_file.hlf", ios::out ); 
                
                ifstream inStream;
                inStream.open(super_file.c_str(), ios::in);
                if( !inStream.is_open() )
                    throw runtime_error( "Failed to open file: " + super_file );
                else
                {
                    string line;
                    while( getline(inStream, line))
                    {
                        string type,sess,conf;
                        conf = "NNNNNN";
                        if(line.find("SUPPRESS_XYULPE") != string::npos )
                        {
                            int pos = line.find("SUPPRESS_XYULPE");
                            conf = line.substr(pos+16,6);
                        }
                        else if(line.find("IN_EOP_CONSTRAINT") != string::npos)
                        {
                            conf = "CCCCCC";
                        }
                                
                        if(line.substr(0,1) == "*" && line.substr(0,2) != "* ")
                        {
                            type = line.substr(1,1);
                            sess = remove_spaces_end(line.substr( 5,9 ));
                        }
                        else if(line.substr(0,1) == "*" && line.substr(0,2) == "* ")
                        {
                            type = "g";
                            sess = remove_spaces_end(line.substr( 4,9 ));
                        }   
                        else
                        {
                            type = "g";
                            sess = remove_spaces_end(line.substr( 3,9 ));
                        }
                          
                        if(line.substr(0,2) != "**" && line.substr(0,6) != "*Earth" )
                        {
                            outstream_hlf << left << setw(9) << setfill(' ') << sess << " ";
                            outstream_hlf << type << " " << conf << endl;
                        }

                        // in case of gsfc arc file
                        if(line.substr(0,1) != "*" )
                        {
                            string session = remove_spaces_end(line.substr( 3,9 ));
                            sources[session] = "keep_none";
                            for(int i=0; i<handling.getLength(); i++)
                            {
                                if(line.find((const char *)handling[i][0]) != string::npos )
                                    sources[session] = (const char *)handling[i][2];
                            }
        //                    if(line.find("SUPPRESS_XYULPE") != string::npos || line.find("IN_EOP_CONSTRAINT") != string::npos)
        //                        sources[session] = "B";
                        }
                    }
                }
                
                outstream_hlf.close();
            }
        }
        ////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////
        // only use a specific defined XXX_Edit.nc file
        string defined_edit_nc = cfg.lookup("defined_nc");
        vector<string> db_names;
        for(int i=0; i<paths.size(); i++)
        {
            cerr << "*** Looking up directory " << paths.at(i);

            dirent *cd;
            DIR *mypath = opendir(paths.at(i).c_str());
            int dir_cnt=0;
            int use_cnt=0;
            while((cd = readdir(mypath))) {
                                
                string name = cd->d_name;

                // not . or .. for folders
                if(name.size() > 5 && name.substr(0,1) != ".")
                {
                    dir_cnt++;
                    name = name.substr(0,9);
                    // used for editfile_criteria
                    string tmp_name = name;
                    
                    // replacing "-" with "X" for vgosDB-names
                    replace_string_in_place( name , "-", "X" );
                    
                    //Abgleich mit DBC-Codes die im FlagFile definiert sind. Zur Auswahl der Sessions.
                    bool not_present, use_session;
                    if((bool)cfg.lookup( "invert_db_exlusion" ))
                    {
                        not_present = (find(excl_db_codes.begin(), excl_db_codes.end(), name.substr(7,string::npos)) != excl_db_codes.end());
                        use_session = (find(excl_sessions.begin(), excl_sessions.end(), name) == excl_sessions.end());
                    }
                    else
                    {
                        not_present = (find(excl_db_codes.begin(), excl_db_codes.end(), name.substr(7,string::npos)) == excl_db_codes.end());
                        use_session = (find(excl_sessions.begin(), excl_sessions.end(), name) == excl_sessions.end());
                    }

                    // check for existing edited *.nc files
                    bool editfile_criteria = true;
                    if(!defined_edit_nc.empty())
                    {
                        string editfile_path = paths.at(i)+"/"+tmp_name+"/ObsEdit/"+defined_edit_nc;
                        if((bool)cfg.lookup( "invert_definied_nc" ))
                            editfile_criteria = !file_exists(editfile_path);
                        else
                            editfile_criteria = file_exists(editfile_path);
                    }

                    ivg::Date db_date(name.substr(0,7),"YYMMMDD");
                    if(not_present && use_session && editfile_criteria ){
                        if(db_date.get_double_mjd() >= start_epoch.get_double_mjd() && db_date.get_double_mjd() <= end_epoch.get_double_mjd())
                        {
                            // check if session is in requested group
                            bool is_group = true;
                            if((bool)cfg.lookup( "invert_groups_exclusion" ))
                            {
                                int false_cnt=0;
                                for(auto &grp: excl_groups)
                                    if( masterfile.is_group(name,grp) == false )
                                        false_cnt++;
                                
                                if(false_cnt == excl_groups.size())
                                    is_group = false;
                            }
                            else
                            {
                                for(auto &grp: excl_groups)
                                {
                                    if( masterfile.is_group(name,grp) == true )
                                    {
                                        is_group = false;
                                        break;
                                    }
                                }
                            }
                            
                            // check if network volume fits and minimal/maximal number of stations
                            ivg::sessinfo info = masterfile.get_session_info(name);  
                            if(is_group && info.volume >= min_volume && info.volume <= max_volume 
                                && info.trf.get_number_stations() >= min_stations && info.trf.get_number_stations() <= max_stations)
                            {
                                db_names.push_back(name);
                                use_cnt++;
                                
                                // check how many sessions are in each generated block
                                if(scaling_modus)
                                {
                                    for(auto &blocki: block_vector)
                                    {
                                        if(db_date.get_double_mjd() >= blocki.start.get_double_mjd() && db_date.get_double_mjd() <= blocki.end.get_double_mjd())
                                        {
                                            blocki.cnt++;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            cerr << ". Result: " << use_cnt << " of " << dir_cnt << " sessions within this folder fit criterions." << endl;
        }
        
        cerr << "*** Finished Looking up ALL directories. " << db_names.size() << " databases found." << endl;
        
        // none handling for every session except already defined by super_file
        if((const char *)cfg.lookup( "default_handling" ) != "")
        {
            for(auto &session: db_names)
                if(sources[session].empty())
                    sources[session] = (const char *)cfg.lookup( "default_handling" );
        }
        ////////////////////////////////////////////////////////////
        // check if sessions need specific handling
        Setting &groups_handling = cfg.lookup( "groups_handling" );
        for(int i=0; i<groups_handling.getLength(); i++)
        {
            string group = (const char *)groups_handling[i][0];
            for(auto &session: db_names)
            {
                if(masterfile.is_group(session,group))
                    sources[session] = (const char *)groups_handling[i][1];
            }
        }  
        ////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////
        // check if sessions need specific handling
        if(arcfile_type != "snx")
        {
            Setting &session_handling = cfg.lookup( "session_handling" );
            for(int i=0; i<session_handling.getLength(); i++)
            {
                string sess = (const char *)session_handling[i][0];
                for(auto &session: db_names)
                {
                    if(sess == session)
                        sources[session] = (const char *)session_handling[i][1];
                }
            }  
        }
        ////////////////////////////////////////////////////////////
        // SORTING
        multimap<double,string> sorted_dbs;
        for(auto &db: db_names)
            sorted_dbs.insert(make_pair(ivg::Date(db.substr(0,7),"YYMMMDD").get_double_mjd(),db));
        
        vector<string> new_dbs;
        for(auto &db: sorted_dbs)
            new_dbs.push_back(db.second);
        
        db_names = new_dbs;
        
//        ref_db_names.push_back(db_names);
//        if(ref_db_names.size() == 2)
//        {
//            vector<string> A =ref_db_names.at(0);
//            vector<string> B =ref_db_names.at(1);
//            
//            for(auto &db_tmp: A)
//            {
//                
//                
//            }
//            
//            for(auto &tmp_db_name: ref_db_names)
//            {
//                (find(special_handlings.begin(), special_handlings.end(), remove_spaces_end(param_other->get_name())) != special_handlings.end());
//                
//                std::vector<string> common;
//                std::set_intersection(ref_db_names.at(0).begin(), ref_db_names.at(0).end(), tmp_db_name.begin(), tmp_db_name.end(), std::back_inserter(common));
//                show_vector(common);
//                cerr << common.size() << endl;
//            }
//            
//        }
        
        // WRITING DATABASES TO ARC FILE
        if(db_names.size() > 0)
        {
            cerr << "*** Writing " << db_names.size() << " selected databased to arcfile " << output_arcfile << endl;
            ofstream outstream ( output_arcfile.c_str(), ios::out ); 
            for(int i=0; i<db_names.size()-1; i++)
            {
                if( sources[db_names.at(i)] != "")
                    outstream << "{ dbname = \"" << db_names.at(i) << "\", version = \"" << version+add_version << "\", handling={" << sources[db_names.at(i)] << "}}," << endl;
                else
                    outstream << "{ dbname = \"" << db_names.at(i) << "\", version = \"" << version+add_version << "\" }," << endl;
            }

                if( sources[db_names.back()] != "")
                    outstream << "{ dbname = \"" << db_names.back() << "\", version = \"" << version+add_version << "\", handling={" << sources[db_names.back()] << "}}" << endl;
                else
                    outstream << "{ dbname = \"" << db_names.back() << "\", version = \"" << version+add_version << "\" }" << endl;
            outstream.close();
                                    
            stringstream ss;
            ss << db_names.size();
                        
            if(add_version != "")
                arc_files_text.push_back("(\"combi\","+ss.str()+")");
            else
                arc_files_text.push_back("(\""+version+"\","+ss.str()+")");
            
            include_text.push_back("@include \""+output_arcfile+"\"");
        }
        else
            throw runtime_error("Writing of arcfile not possible. No db_names which match the criteria.");
    
    }
    
    // only in case of snx-arc-files
    if(arcfile_type == "snx")
    {
        // in case of scaling modus, the arc_files variable needs to be written differently
        if(scaling_modus)
        {
            arc_files_text.clear();
            
            int cnt=1;
            for(auto &blocki: block_vector)
            {
                if(blocki.cnt>0)
                {
                    stringstream ss;
                    ss << "(\"" << "BLOCK" << cnt << "\"," << blocki.cnt << ")";
                    cerr << ss.str() << endl;
                    arc_files_text.push_back(ss.str());
                    cnt++;
                }
            }
        }
        
        // generating main.arc file containing the different arc_files for the different global solutions
        int N = arc_files_text.size();

        string main_arc_file = output_arcfile.substr(0,output_arcfile.find_last_of("/\\")+1)+"/main.arc";
        cerr << "*** Writing main.arc file" << main_arc_file << endl;
        ofstream outstream ( main_arc_file.c_str(), ios::out ); 
        outstream << "arc_files = (";
        for(int i=0; i<N-1; i++)
            outstream << arc_files_text.at(i) << ",";

        outstream << arc_files_text.at(N-1) << ");" << endl;

        outstream << endl;
        outstream << "output_folder = \"" << new_folder.substr(new_folder.find_last_of("/\\")+1) << "\";" << endl;
        outstream << endl;

        N = include_text.size();
        outstream << "sessions = (" << endl;
        for(int i=0; i<N-1; i++)
            outstream << include_text.at(i) << "," << endl;

        outstream << include_text.at(N-1) << endl;
        outstream << ");";

        outstream.close();
    }
    
    cerr << "*** Finished" << endl;
    
}

// CODE FROM SESSION TO USE SOURCE STATS FOR SOURCE-SPECIFIC SCALING

                        //only if session specific handling exists, use it
//                        if((bool)(*_setup).exists("source_stats"))
//                        {
//                            double max_num = 0.0;
//                            double src_num = 0.0;
//                            for(int s=0; s<(*_setup)["source_stats"].getLength(); s++){
//                                
//                                string tmp_src = (*_setup)["source_stats"][s]["source"];
//                                double tmp_num = (double)(*_setup)["source_stats"][s]["num"];
//                                
//                                if(src_iter->get_name(ivg::srcname::ivs) == tmp_src)
//                                    src_num = tmp_num;
//                                
//                                if(tmp_num> max_num)
//                                    max_num = tmp_num;
//                            }
//                            cerr << src_iter->get_name(ivg::srcname::ivs) << " : " << src_num << " : " << max_num <<  endl;
                        
//                            if(src_num == 0.0){
//                                log<WARNING>("!!! Source not found") % src_iter->get_name(ivg::srcname::iers) % " " % src_iter->get_name(ivg::srcname::ivs);
//                                src_num = max_num / 2.0;
//                            }