
#include <cstdlib>
#include <vector>

#include <tclap/CmdLine.h>
#include <libconfig.h++>

#include "logger.h"
#include "db_download.h"
#include "masterfile.h"
#include "ltn/log2nc.h"
#include "date.h"

loglevel g_verbose;

int main(int argc, char** argv) {
    
    /*
     *  Parse args -------------------------------------------------------------
     */
    
    //Construktor needs: message, delimeter, version
    TCLAP::CmdLine cmd("get_db", ' ', "0.1");
    
    //flag 	- The one character flag that identifies this argument on the command line.
    //name 	- A one word name for the argument. Can be used as a long flag on the command line.
    //desc 	- A description of what the argument is for or does.
    //req 	- Whether the argument is required on the command line.
    //value 	- The default value assigned to this argument if it is not present on the command line.
    //typeDesc 	- A short, human readable description of the type that this object expects. This is used in the generation of the USAGE statement. The goal is to be helpful to the end user of the program.
    //v 	- An optional visitor. You probably should not use this unless you have a very good reason. 

    // verbose level
    TCLAP::ValueArg<int> verboseArg( "v","verbose","Verbose Level (0=NOTHING, 1=INFO, 2=DETAIL, 3=RESULT, 4=WARNING, 5=ALL)",false,4,"int");
    cmd.add( verboseArg );
  
    // controlfile
    TCLAP::ValueArg<std::string> cntArg( "c","controlfile","/home/ascot/cnt",false,"../src/progs/get_db/get_db.cfg","string");
    cmd.add(cntArg);  
    
    // force downloading db
    TCLAP::SwitchArg  force_switch( "f","force","force download the databases even if it is already on local directory", false );
    cmd.add( force_switch );
    
    // dbname
    TCLAP::ValueArg<std::string> dbArg( "e","explicitDB","e.g. 93NOV05XU",false,"","string");       
    cmd.add(dbArg);

    // Parse the argv array.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg. 
    g_verbose = (loglevel) verboseArg.getValue();

    std::string controlfile = cntArg.getValue();
    
    bool force = force_switch.getValue();
    
    std::cout << " \n =================================================================== " << std::endl;
    std::cout << "     > > > > > > > > > > ivg::ASCOT (get_db) < < < < < < < < < <    " << std::endl << std::endl; 
    
    // read control file and initialize variables
    // loading control file (only clocks estimated, group and single band)
    libconfig::Config cfg;
    
    try
    {
        cfg.readFile( controlfile.c_str() );
    }
    catch( libconfig::ParseException & err )
    {
        std::cerr << "libconfig::" << err.what() << " in " << err.getFile() << " at line " << err.getLine() << std::endl;
        exit( -1 );
    }
    libconfig::Setting& setup= cfg.lookup( "setup" );
    
    // get controls from controlfile
    std::string local_log_path = (const char*)setup["ltndir"]["log"];
    std::string local_vgos_path = (const char *) setup["local_vgos"];
    std::string remote_log_path = (const char*)setup["remotedir"]["log"];
    std::string remote_vgos_path = (const char*)setup["remotedir"]["vgos"];

    std::string ltn_data = (const char*)setup["ltndir"]["data"];   

    std::string masterfile_path = (const char*)setup["masterfiles"];
    
    std::string reference = (const char*)setup["reference"];  

    int start_year = setup["period"]["start"];
    int end_year = setup["period"]["end"];
    
    bool load_log = (bool)setup["load_log"];
    
    std::string mastertype_str = (const char*)setup["mastertype"];
    ivg::mastertype mtype;
    if( mastertype_str.compare("both") == 0)
        mtype = ivg::mastertype::both;
    else if( mastertype_str.compare("regular") == 0)
        mtype = ivg::mastertype::regular;
    else if( mastertype_str.compare("intensive") == 0)
        mtype = ivg::mastertype::intensive;
    else
    {
        log<WARNING>("!!! Invalid mastertype. Using ivg::mastertype::both");
        mtype = ivg::mastertype::both;
    }
    
    //initzialise masterfiles
    ivg::Masterfile master(masterfile_path,mtype, start_year, end_year);   
    
    ivg::Db_download d(local_log_path,local_vgos_path,remote_log_path,remote_vgos_path,&master,ltn_data);

            // first case: Masterfile is reference
            if(reference.compare("masterfile") == 0)
            {  
                
                if (master.get_sessions()->empty()){
                    log<WARNING>("!!! Masterfile error. Make sure the path in controlfile is right and the files are existing");
                }
                else
                {
                    // fill session vector
                    std::vector<std::string> database;
                    if(dbArg.isSet())
                    {
                        database.push_back(dbArg.getValue());
                    }
                    else{
                        for(ivg::sessinfo &info : *master.get_sessions()){

                            ivg::Date d;
                            d.now();

                            // Only add session from masterfile to list if it is not in the future
                            double diff = info.date.get_double_mjd() - d.get_double_mjd();
                            if( diff < 0 ){
                                database.push_back(info.dbname);
                            }
                            else{
                                log<DETAIL>("*** session ") % info.dbname % " is " % to_string((int)diff) %" days in the future and not added to download list";
                            }

                        }

                    }

                    // download all sessions in session vector
                    unsigned n = database.size();
                    for(unsigned i = 0; i < n; ++i)
                    {
                        std::cout << "\n> > > > > > > > > > " << database[i] << "(" << i+1 << "/" << n+1 << ") < < < < < < < < < <" << std::endl; 
                        d.download_vgosDB(database[i],force,load_log);
                    }
                }
            }
            // second case: FTP Server is reference
            else if(reference.compare("remote") == 0)
            {
                if(dbArg.isSet())
                {
                    log<WARNING>("!!! option -e works only if refernce is set to masterfile in the config file");
                }
                else
                {
                    d.update_vgosDB(start_year,end_year);
                }
                
            }
            else
            {
                log<WARNING>("!!! invalid reference for  databases. exit");
                exit( -1 );
            }

    return 0;
}

