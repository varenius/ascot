#include "session_inout.h"
#include "wrapper.h"
#include <algorithm>

namespace ivg
{
    
Session_inout::Session_inout()
{
    
}
    
// ...........................................................................  
Session_inout::Session_inout(const string type, ivg::Masterfile masterfile)
{

    vector<string> types = {"ngs","vgosdb","snx","utas","skd","cata"};

    if ((find(types.begin(),types.end(),type)==types.end()))
        throw runtime_error("Session_inout::Session_inout( const string ): type "+type+" unknown. No initialization possible.");
    else
    {
        _masterfile = masterfile;
        _type = type;
        _wrapper_ptr = nullptr;
    }
    
}
// ...........................................................................   
std::string Session_inout::load(ivg::Session *session_ptr,Setting *setup,const string name,const string version)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::load( ivg::Session *, Setting *, const string , const string )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // setting the type of session for the session (little bit redundant)
    session_ptr->_type = _type;
   
    //build dir based on selected dir    
    string dir = (const char *) get_list_element((*setup)["datadirs"],(*setup)["session_type"])[2];
    
    //read ngs or snx or ...
    if (_type=="ngs")
    {

        // create 4 digit year from DB name, in case of ngs / vgosdb needed
        int year = stoi(name.substr(0,2));
        if (year<79)
            year += 2000;
        else
            year += 1900;

        // construct file name
        stringstream path;
        path<<dir<<"/"<<year<<"/"<<name
                <<"_N"<<setfill('0')<<setw(3)<<right<<version;

        _session_path = path.str();

        log<INFO>("*** Loading ngs database: ")%path.str();

        // read data from ngs into session object
        // sets following member variables of the session object
        // _start, _end, _nobs, _trf, _crf, _scans, _param_list
        _read_ngs(session_ptr,setup,path.str());

    }
    else if (_type=="vgosdb")
    {
        // create 4 digit year from DB name, in case of ngs / vgosdb needed
        int year;
	if (name.size()>9)
	  year=stoi(name.substr(0,4));
	else
	  {
	    year=stoi(name.substr(0,2));
	    if (year<79)
	      year += 2000;
	    else
	      year += 1900;
	  }
        // copy name to be able to change it
        string dbname = name;

        // replace X in databasename with -
        // not needed anymore, because vgosDB_iIVS uses e.g. XX, not -X
//        size_t posi = dbname.substr(0,8).find("X");
//        if (posi!=string::npos)
//            dbname.replace(posi,1,"-");
	
        // construct file name
        stringstream folderpath;
        folderpath<<dir<<"/"<<year<<"/"<<dbname<<"/";

        _session_path = folderpath.str();
        session_ptr->_session_path = _session_path;
	
        // read data from vgosdb into session object
        // sets following member variables of the session object
        _read_vgosdb(session_ptr,setup,folderpath.str());
//        _read_utas_apriori(session_ptr,setup,"/home/iddink/lucia.dat");
        
    }
    // in case of snx
    else if (_type=="snx")
    {
        stringstream path;
        //if the version is empty, use the name as full path
        if (version.empty())
            path << name;
        else
        {
            // construct file name      
            path << dir << "/" << version << "/" << name;

            if (name.size()==8)
                path<<"__"<<version<<".snx";
            else
                path<<"_"<<version<<".snx";

        }

        _session_path = path.str();
        // redundant storage in session and session-inout
        session_ptr->_session_path = path.str();

        // read data from snx into object
        // sets following member variables of the session object
        // _start, _end, _trf, _crf, _param_list
        _read_snx(session_ptr,setup,path.str());
    }
    // in case of gsfc specific catalog file
    else if (_type=="sou")
    {
        log<INFO>("*** Loading sou database: ") % name;
        _read_sou(session_ptr,setup,name);
    }
    else if (_type=="utas")
    {
        log<INFO>("*** Loading utas database: ") % name;
        _read_utas_apriori(session_ptr,setup,name);
    }
    else if (_type=="skd")
    {
        
        sessinfo si = _masterfile.get_session_info( session_ptr->_name );
        ivg::Date date = ivg::Date(name.substr(0,7),"YYMMMDD");
        std::string skdname;
        if( si.name == "UNKNOWN" ){
            log<WARNING>("!!! ") % session_ptr->_name % " not in masterfile. Create folder name bsed on date and name_prefix in controlfile";
            stringstream ss;
            ss << std::setw(3) << std::setfill('0') << date.get_int_doy();
            skdname =  (const char*)(*setup)[ "SKED" ]["name_prefix"] + name.substr(0,2) + ss.str();
        } else {
            //string skdname = _masterfile.get_session_info( date, 0.0, ivg::mastertype::intensive).skdname;
            skdname = _masterfile.get_session_info( session_ptr->_name ).skdname;
        }
        
        if( (bool)(*setup)[ "SKED" ]["use_year_ses_structure"] ){
            dir+= "/"+date.get_date_time("YYYY")+"/"+skdname+"/";
        } else {
            dir+= "/"+date.get_date_time("YYYY")+"/";
        }
               
        std::string skdfile = dir + "/"+skdname+".skd";
        session_ptr->_session_path = skdfile;
        
        // session will be created based on config setting in SKED-block
        if((*setup).exists( "SKED" ) && (bool)(*setup)[ "SKED" ]["apply"])
            _init_skd_session(session_ptr, setup, skdfile);
        else if((*setup).exists( "SKED" ) && !(bool)(*setup)[ "SKED" ]["apply"] && 
               ((*setup).exists( "SIM" ) && (bool)(*setup)[ "SIM" ]["apply"] && (bool)(*setup)["SIM"]["load"][0] == false ) 
               || ((*setup).exists( "SKED" ) &&  (bool)(*setup)[ "SKED" ]["createPlots"] && ( !(*setup).exists( "SIM" ) || !(bool)(*setup)[ "SIM" ]["apply"])) )
        {
            // define some default-name
            session_ptr->_name = skdname;
            _read_skd(session_ptr, setup, skdfile);
        }
        else
            throw runtime_error("void Session_inout::load(...): Unexpected SIM = false, SKED = false but type = skd setup.");
    }
    
    return dir;

#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::load( ivg::Session *, Setting *, const string , const string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif   
}
// ...........................................................................
void Session_inout::_read_sou(ivg::Session *session_ptr,Setting *setup,const string path)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::_read_sou(ivg::Session *, Setting *, const string )"<<endl;
    tictoc tim;
    tim.tic();
#endif   

//    ifstream inStream(path.c_str(),ios::in);
//    if (!inStream.is_open())
//        throw runtime_error("void _read_sou(ivg::Session *, Setting *, const string ): Failed to open file "+path);
//    else
//    {
//
//        vector<ivg::Source> sources;
//
//        string line;
//        while (getline(inStream,line,'\n'))
//        {
//
//            if (line.substr(0,7)=="SOU_GCO")
//            {
//
//                string ivs_name = remove_spaces_end(line.substr(10,8));
//
//                // reading right ascenion
//                int ra_h = stoi(line.substr(24,2));
//                int ra_m = stoi(line.substr(27,2));
//                double ra_s = s2d(line.substr(30,11));
//
//                //reading declination
//                int dec_g = stoi(line.substr(61,3));
//                int dec_m = stoi(line.substr(65,2));
//                double dec_s = s2d(line.substr(68,10));
//
//                double std_ra_mas = s2d(line.substr(47,8));
//                double std_dec_mas = s2d(line.substr(84,10));
//                double corr = s2d(line.substr(98,6));
//
//                ivg::Source tmp_src(ivg::srcname::ivs,ivs_name,ra_h,ra_m,ra_s,dec_g,dec_m,dec_s,std_ra_mas,std_dec_mas,corr);
//                sources.push_back(tmp_src);
//            }
//        }
//        // create CRF based on the creates source vector
//        session_ptr->_crf = ivg::Crf(*setup,sources);
//    }

#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::_read_sou(ivg::Session *, Setting *, const string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif    
}
// ...........................................................................
ivg::Matrix Session_inout::init_xka_catalog(ivg::Session *session_ptr, Setting *setup,  const string path)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::_read_cata(ivg::Session *session_ptr,Setting *setup,const string path)"<<endl;
    tictoc tim;
    tim.tic();
#endif   

    log<INFO>("*** START-Loading cata ") % path;
    
    ifstream inStream(path.c_str(),ios::in);
    if (!inStream.is_open())
        throw runtime_error("void _read_cata(ivg::Session *session_ptr,Setting *setup,const string path): Failed to open file "+path);
    else
    {
        vector<ivg::Source> sources;

        string line;
        // get first line containing number of sources
        getline(inStream,line,'\n');
        int num_src = stoi(remove_spaces(line.substr(0,9)));
        log<INFO>("*** Loading ") % num_src % " sources from X/KA catalog";
        // skip next trash-line containing "xka.srf"
        getline(inStream,line,'\n');
                     
        // 1. Section 1: IERS catalog with one source per line
        for(int i=1; i<=num_src; i++)
        {
            // get each line containing JPL-source-name, 
            getline(inStream,line,'\n');
            
            string jplname = remove_spaces_end(line.substr(0,13));
            
            // replacing all D to E because of different exponential appearence
            std::replace( line.begin(), line.end(), 'D', 'E');
            
            double ra,dec,ra_std,dec_std,corr,mean,low,high;
            int sessions,delays;
            istringstream data_line( line.substr(14) );
            data_line >> ra >> dec >> ra_std >> dec_std >> corr >> mean >> low >> high >> sessions >> delays;

            sources.push_back(ivg::Source(ivg::srcname::jpl,jplname,ra,dec,ra_std,dec_std,corr));
        }
        
        // next line after the first block contains number of parameters
        getline(inStream,line,'\n');
        if( std::stoi(remove_spaces(line)) != num_src*2 )
            throw runtime_error( "void Session_inout::_read_cata(...): Number of parameter and sources does not match.");

        // 2. covariance header and parameter list
        // containing all information already set in block 1. Just to check if block 1 and 2 corresponde to eachother
        // columns 1 to 3 = [ra, ra_std, mean_epoch]
        // columns 1 to 3 = [dec, dec_std, mean_epoch]
        for(int i=1; i<=num_src; i++)
        {
            // get line containing right ascension
            getline(inStream,line,'\n');

            // replacing all D to E because of different exponential appearence
            std::replace( line.begin(), line.end(), 'D', 'E');

            double ra,ra_std,ra_mean;
            istringstream ra_line( line.substr(26) );
            ra_line >> ra >> ra_std >> ra_mean;       

            // get line containing declination
            getline(inStream,line,'\n');

            // replacing all D to E because of different exponential appearence
            std::replace( line.begin(), line.end(), 'D', 'E');

            double dec,dec_std,dec_mean;
            istringstream dec_line( line.substr(26) );
            dec_line >> dec >> dec_std >> dec_mean;   

            if(ra != sources.at(i-1).get_ra0() || dec != sources.at(i-1).get_dec0() || ra_mean != dec_mean )
                throw runtime_error( "void Session_inout::_read_cata(...): Information from cata block 1 and 2 does not match.");
        }
                
        // 3. covariance information stored as inverse covariance matrix (i.e. weight matrix)
        // Number of elements ( = nsrc * (nsrc+1) /2 )
        // inverse_covariance values (radians**2) --duplicate parts of symmetric matrix not stored

        int num_para = num_src * 2; // = 1262 parameter
        int num_ele = num_para * num_para; // 796953*2
        
        // first all values are stored in one long vector
        vector<double> values;
        values.resize(num_ele);
        
        int i_row=-1;
        int i_col=-1;
        int maxposi = 0;
        while(getline(inStream,line,'\n'))
        {
            // replacing all D to E because of different exponential appearence
            std::replace( line.begin(), line.end(), 'D', 'E');
            double val;
            stringstream tokenizer(line);
            while(tokenizer >> val)
            {
                if(i_col == i_row)
                {
                    i_col++;
                    i_row=0;
                }
                else
                    i_row++;
                
                // save element to correct position (upper and lower symmetric matrix)
                values.at(i_row + i_col * num_para) = val;           
                values.at(i_col + i_row * num_para) = val;  
            }
        }
        
        // create Qxx based on values-vector
        ivg::Matrix Qxx( values.begin(), values.end(), num_para, num_para);
        
        // first: check if correlations from block 1 and block 3 corresponde to eachother
        int j = 0;
        for(int i=0; i<num_para; i+=2)
        {
            double test_A = Qxx(i,i+1) / ( sqrt(Qxx(i,i)) * sqrt(Qxx(i+1,i+1)));
            double test_B = sources.at(j).get_corr();
            
            if(abs(test_A-test_B)>10e-16)
                throw runtime_error( "void Session_inout::_read_cata(...): Correlation of block 1 and 3 not identical of source "+sources.at(j).get_name(ivg::srcname::jpl));
                
            j++;
        }
        
        // second: check of dec0 standard deviation from block 1 correspondes to it from block 3
        j=0;
        for(int i=0; i<num_para; i+=2)
        {
            double test_A_dec0 = sqrt(Qxx(i+1,i+1));
            double test_B_dec0 = sources.at(j).get_sigma_dec0();
            
            if(abs(test_A_dec0-test_B_dec0)>10e-16)
                throw runtime_error( "void Session_inout::_read_cata(...): dec0 standard deviation of block 1 and 3 not identical of source "+sources.at(j).get_name(ivg::srcname::jpl));

            j=j+1;
        }
        
        log<INFO>("*** END-Loading cata ") % path;
    
        return Qxx;
    }

#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::_read_cata(ivg::Session *session_ptr,Setting *setup,const string path)"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif    
}
// ...........................................................................
void Session_inout::_read_vgosdb(ivg::Session *session_ptr, Setting *setup, const string directory)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Session_inout::_read_vgosdb(ivg::Session *session_ptr, Setting *setup, const string directory)" << endl; 
   tictoc tim;
   tim.tic();
#endif

    string band_str = ivg::band_to_string(session_ptr->_band_type);

    log<NOTHING>("*** START-Loading vgosDB (") % band_str % "-band) : " % directory;
   
    string editing = (const char *)(*setup)["vgosdb_editing"];
    
    bool use_wrapper = (bool)(*session_ptr->_setup)["use_wrapper"];
    
    ivg::Vgosdb vgosdb;
   
    if(use_wrapper){
        if(_wrapper_ptr){
            if(session_ptr->_name.compare(_wrapper_ptr->get_dbName()) == 0 ){
                 vgosdb = ivg::Vgosdb(directory,_wrapper_ptr);
            }
            else{
                log<WARNING>("!!! session in wrapper and vgosdb are different: ") % session_ptr->_name % " vs. "% _wrapper_ptr->get_dbName();
                exit(0);
            }
        }
        else{
            log<WARNING>("!!! use wrapper is set in controlfile, but there is no valid pointer to wrapper");
        }
    }
    else
    {
       vgosdb = ivg::Vgosdb(directory);
       use_wrapper = false;
    }
   
    std::string head_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::Head,session_ptr->_band_type) )
        head_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::Head,session_ptr->_band_type);  
    else
        head_filename = "Head";
  
    // get effective_duration (only if wrapper is used because the names are highly variable)
    // not used right now; should be set for each observation for correct simulation
    vector<double> geocenter_delay;
    vector<double> effective_duration;
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::CorrInfo,session_ptr->_band_type) )
    {
        try
        {
            effective_duration = vgosdb.get_vector<double>("Observables",_wrapper_ptr->get_file(ivg::wrapper_entries::CorrInfo,session_ptr->_band_type),"EffectiveDuration");
            geocenter_delay = vgosdb.get_vector<double>("Observables",_wrapper_ptr->get_file(ivg::wrapper_entries::CorrInfo,session_ptr->_band_type),"GeocMBD");
        }
        catch(std::exception& e)
        {
            log<WARNING>("!!! ") % (string)e.what() % ". Setting it to 0.0";
        }
    }
    
    // get information about start and end epoch
    ivg::Matrix dates = vgosdb.get_matrix("",head_filename,"iUTCInterval");

    session_ptr->_start =  ivg::Date(dates(0,0),dates(0,1),dates(0,2),dates(0,3),dates(0,4));
    session_ptr->_end =  ivg::Date(dates(1,0),dates(1,1),dates(1,2),dates(1,3),dates(1,4));

    // set number of observations from Head.nc
    int nobs = vgosdb.get_scalar<int>("",head_filename,"NumObs");
    session_ptr->_nobs = nobs;

    // save number of observations before _nobs is changed due to invalid obs
    session_ptr->_nobs_orig = nobs;

    // read reference frequency
    vector<double> ref_freq;
    std::string reffreq_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::RefFreq,session_ptr->_band_type) )
        reffreq_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::RefFreq,session_ptr->_band_type);
    else
        reffreq_filename = "RefFreq_b"+band_str;
   
    ref_freq = vgosdb.get_vector<double>("Observables",reffreq_filename,"RefFreq");
    if ((ref_freq.size()>0)&&(ref_freq.size()<nobs))
      for (int j=ref_freq.size();j<nobs;j++)
	ref_freq.push_back(ref_freq.at(0));

    
    if ( session_ptr->_ambigRes ){
        // get ambiguity_spacing
      if ( session_ptr->_phaseDelay ){
	
	ivg::Matrix numambig(nobs,1,0.0);
	
	if (ref_freq.size()==1) {
	  
	  for (int j=0;j<nobs;j++){
	    
	    numambig(j,0)=(1/(ref_freq.at(0)*1.0e6));
	  }
	} else {
	  
	  for (int j=0;j<nobs;j++){
	    if (ref_freq.at(j)>0)
	      numambig(j,0)=(1/(ref_freq.at(j)*1.0e6));
	    else
	      numambig(j,0)=(1/(ref_freq.at(0)*1.0e6));
	  }
	}
	session_ptr->_ambiguity_spacing = numambig;
   
      } else {
       
        if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::AmbigSize, session_ptr->_band_type) ){
            session_ptr->_ambiguity_spacing =  ivg::Matrix(vgosdb.get_vector<double>("Observables", _wrapper_ptr->get_file(ivg::wrapper_entries::AmbigSize,session_ptr->_band_type) ,"AmbigSize"));
        } else {
            session_ptr->_ambiguity_spacing =  ivg::Matrix(vgosdb.get_vector<double>("Observables","AmbigSize_b"+band_str,"AmbigSize"));
        }

        
      }
      session_ptr->_num_ambig.resize( nobs, 1, 0.0 );
        /*
         * TODO wtf is wrong here? vgosdb.get_matrix says "NetCDF: Attribute not found" and exits and
         * get_vector<double> with the same parameters works fine 
         * reason: REPEAT attribute is used in get_Matrix if DIM==1 although it is not defined
         */
      string EffFreqType;
      if (session_ptr->_phaseDelay )
	EffFreqType="FreqPhaseIono";
      else
	 EffFreqType="FreqGroupIono";
       
        if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::EffFreq, session_ptr->_band_type) )
            session_ptr->_eff_freq = ivg::Matrix(vgosdb.get_vector<double>("ObsDerived",_wrapper_ptr->get_file(ivg::wrapper_entries::EffFreq,session_ptr->_band_type),EffFreqType));
        else if(vgosdb.does_file_exist("ObsDerived","EffFreq_b"+band_str))
            session_ptr->_eff_freq = ivg::Matrix(vgosdb.get_vector<double>("ObsDerived","EffFreq_b"+band_str,EffFreqType));
	else
	  {
	    vector<vector<vector<double> > > tst;
	    
	    tst=vgosdb.get_vector_3d_data<double>("Observables","ChannelInfo_b"+band_str,"ChanAmpPhase");
	    //std::cout << tst.size() << " " << tst[0].size() << " " << tst[0][0].size() <<endl;;
	    ivg::Matrix chanAmp = ivg::Matrix(tst.size(),tst[0].size(),0.0);
	    ivg::Matrix chanPhas = ivg::Matrix(tst.size(),tst[0].size(),0.0);
	    //std::cout << chanPhas.rows() << "x"  << chanPhas.cols() <<endl;
	    for (int i=0;i<tst.size();i++) {
	      for (int j=0;j<tst[i].size();j++) {
		chanAmp(i,j)=tst[i][j][0];
		chanPhas(i,j)=tst[i][j][1];
	      }
	      
	    }
	    
	    int bitsampl;
	    if (vgosdb.does_variable_exist("Observables","ChannelInfo_b"+band_str,"BITSAMPL"))
		bitsampl = vgosdb.get_scalar<int>("Observables","ChannelInfo_b"+band_str,"BITSAMPL");
	    else
	      bitsampl = 1;
	    int sampl_rate = vgosdb.get_scalar<int>("Observables","ChannelInfo_b"+band_str,"SampleRate");
	    ivg::Matrix chanfreq = vgosdb.get_matrix("Observables","ChannelInfo_b"+band_str,"ChannelFreq");
	    vector<vector<vector<int> > > numAp = vgosdb.get_vector_3d_data<int>("Observables","ChannelInfo_b"+band_str,"NumAp");

	    
	    double bw=double(sampl_rate)/(2*bitsampl)*1e-6;
	    session_ptr->_eff_freq = ivg::Matrix(1,tst.size(),0.0);
	    for (int i=0;i<tst.size();i++)
	      {
		double f0;
		if (ref_freq.size()==1)
		  f0=ref_freq.at(0);
		else
		  f0=ref_freq.at(i);
		double sumr=0,sumrf=0,sumrf2=0,sumrdf=0,sumrf2p=0,sumrfp=0,sumrfdfp=0;
		for (int j=0;j<tst[i].size();j++)
		  {
		    double tmpfreq=chanfreq(i,j),tmpr=0;
		    if (numAp[i][j][0]>0)
		      {
			tmpfreq-=bw/2;
			tmpr+=tst[i][j][0];
		      }
		    if (numAp[i][j][1]>0)
		      {
			tmpfreq+=bw/2;
			tmpr+=tst[i][j][0];
		      }
		    if (chanfreq(i,j)!=0)
		      {
			sumr+=tmpr;
			sumrf+=tmpr*tmpfreq;
			sumrf2+=tmpr*tmpfreq*tmpfreq;
			sumrdf+=tmpr/tmpfreq;
			 sumrf2p+=tmpr*(tmpfreq-f0)*(tmpfreq-f0);
			 sumrfp+=tmpr*(tmpfreq-f0);
			 sumrfdfp+=tmpr*(tmpfreq-f0)/tmpfreq;
		      }
		    
		  }
		
		if (sumr>0) {
		  if (session_ptr->_phaseDelay )
		    session_ptr->_eff_freq(0,i)=sqrt(f0*(sumrf2p*sumr-sumrfp*sumrfp)/(sumrf2p*sumrdf-sumrfp*sumrfdfp));
		  else
		     session_ptr->_eff_freq(0,i)=sqrt(-(sumr*sumrf2-sumrf*sumrf)/(sumr*sumr-sumrf*sumrdf));

		}
		else
		  session_ptr->_eff_freq(0,i)=0;
	        
	      }
	  }
	

    }

    // get TRF stations and CRF sources names
    ivg::Date date = session_ptr->_start;
    vector<string> station_names = vgosdb.get_string("",head_filename,"StationList");
    vector<string> source_names = vgosdb.get_string("",head_filename,"SourceList");
     
   // eliminating stations based on "elim_sta" in config file
//   for( int i=0; i<(*setup)["elim_sta"].getLength(); ++i )
//   {
//      string sta = (*setup)["elim_sta"][i];
//      vector<string>::iterator it = find( station_names.begin(), station_names.end(), sta );
//      if( it != station_names.end() )
//      {
//         log<INFO>("*** Eliminating Station ") % *it % " from vgosDB.";
//         station_names.erase(it);
//      }
//   }
   
    for(auto &st: station_names)
         replace_string_in_place( st , " ", "_" );
   
    // in case of SH the StationList and SourceList in Head.nc is not correct
    if(directory.find("VASCC2015_SH_VGOSDB")!=string::npos)
    {
        log<WARNING>("!!! EXTREME SPECIAL VASCC2015_SH_VGOSDB CASE: Replacing wrong station and source names manually.");
        station_names =  vector<string>({"HART15M","HOBART26","KATH12M","WARK12M"});
        source_names = vector<string>({"0230-790"});
    }
      
   // works with all diferent kind of trfs (vtrf2014, itrf2014, itrf2008)
    session_ptr->_trf =  ivg::Trf( (*setup), station_names, ivg::staname::ivs_name, true, date.add_days(-10.0), date.add_days(10.0) );
      
    // check if _trf corresponde to number of stations in original vgosdb
    if(session_ptr->_trf.get_number_stations() != station_names.size())
        throw runtime_error( "void Session_inout::_read_vgosdb(ivg::Session *session_ptr, Setting *setup, const string directory): _trf not correctly initialized. Wrong number of stations.");

    // the use of external meteorological data or zenith hydrostatic delays (ZHDs)
    // is usually considered in init_displacement(...) in ivg::Trf;
    // however, in case of raytraced delays, the external data have to be loaded 
    // immediatly after the TRF initialization, since there is a need of the session name    
    if( (bool)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][0] &&
        string((const char*)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][2]) == "raytracing" )           
    { 
        Setting &definitions = (*session_ptr->_setup)["definitions"];
        std::string path = definitions["troposphere"]["external_meteo_data"]["raytracing"]; 
        std::string rt_file = path+"/"+session_ptr->_name+".trp";
        ivg::parser::raytraced_delays( &session_ptr->_trf, rt_file );        
    }
    
    // shows the main feautres of the loaded trf
    session_ptr->_trf.log_data_info_table();

    // if source aprioris exist, they can be used with flag crf = "vgosdb"; the original apriori values from vgosDB will be used
    vector<ivg::Source> sources;
    std::string source_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::Sources,session_ptr->_band_type) )
        source_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::Sources,session_ptr->_band_type);
    else
        source_filename = "Source";
    
    if(vgosdb.does_file_exist("Apriori",source_filename))
    {
         vector<string> apri_source_names = vgosdb.get_string("Apriori", source_filename,"AprioriSourceList");
         ivg::Matrix apri_source_matrix = vgosdb.get_matrix("Apriori", source_filename,"AprioriSource2000RaDec");

         //create CRF
         for(int s=0; s<apri_source_names.size(); s++)
             sources.push_back(ivg::Source(ivg::srcname::ivs, apri_source_names.at(s), apri_source_matrix(s,0), apri_source_matrix(s,1)));
    }
    else
    {
         for(int s=0; s<source_names.size(); s++)
             sources.push_back(ivg::Source(ivg::srcname::ivs, source_names.at(s) ));
    }
    
    session_ptr->_crf = ivg::Crf( *setup, sources, &session_ptr->_ephem );
   
    // in case of satellite observations, we need to transform the sp3 orbits from trf to crf
    for(auto &src: session_ptr->_crf)
    {
        if(src.get_type() == ivg::srctype::satellite)
        {
            // get original sp3 xyz coordinates in itrf
            ivg::Matrix sp3 = src.get_sp3();
            for(int r=0; r<sp3.rows(); r++)
            {
                ivg::Date time(sp3(r,0)); // GPS-Time
                time.add_secs(-(time.get_leap_sec() - 19.0)); // correct for offset between GPS-Time and UTC
                // get rotation matrix, valid at specific time
                ivg::Matrix trf2crf = session_ptr->_eops.form_crf2trf( time, false ).transpose();
                ivg::Matrix xyz_trans = trf2crf * sp3.get_sub(r,1,r,3).transpose();
                sp3.set_sub(r,1,xyz_trans.transpose());
            }
            // save new transformed sp3
            src.set_sp3(sp3);
        }
    }  
   
    //Initializing parameter list (Stations, Sources, EOPs) including epoch and aprioris
    session_ptr->_param_list = Param_list(session_ptr->_trf, session_ptr->_crf, session_ptr->_eops, ivg::Date( 0.5*( session_ptr->_start.get_double_mjd()+session_ptr->_end.get_double_mjd() ) ));
   
   
    map< ivg::paramtype,multimap<string,double> > breaks;
    // check if a clock break need to be added manually, not stored in original vgosDB
    if((bool)(*(session_ptr->_handling)).exists("handling") && (bool)((*(session_ptr->_handling))["handling"]).exists("add_cbr"))
    {   
        for( int i=0; i<(*(session_ptr->_handling))["handling"]["add_cbr"].getLength(); ++i )
        {
            string station = (*(session_ptr->_handling))["handling"]["add_cbr"][i][0];
            string date_str = (*(session_ptr->_handling))["handling"]["add_cbr"][i][1];

            if(date_str.size() > 18) // in this case, it should be something like 2004/11/16-04:03:00
            {
                int y = std::stoi( date_str.substr( 0,4 ) );
                int m = std::stoi( date_str.substr( 5,2 ) );
                int d = std::stoi( date_str.substr( 8,2 ) );
                int h = std::stoi( date_str.substr( 11,2 ) );
                int min = std::stoi( date_str.substr( 14,2 ) );
                double sec = s2d(date_str.substr( 17,2 ) );

                breaks[ivg::paramtype::cbr].insert(std::pair<string,double>(station.substr(0,8), ivg::Date(y,m,d,h,min,sec).get_jd()) );
            }
            else // in this case, it should be MJD
                breaks[ivg::paramtype::cbr].insert(std::pair<string,double>(station.substr(0,8), ivg::Date(s2d(date_str)).get_jd()) );
            
            log<INFO>("*** Insert additional clock break due to defintion in [handling][add_cbr]");
        }
    }
   
    // check if the session contains clock breaks in vgosDB
    string cb_filename;
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::ClockBreak,ivg::band::X))
        cb_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::ClockBreak,ivg::band::X);
    else
        cb_filename = "ClockBreak";   
    bool use_cb_file=true;
    if(setup->exists("UseClockbreakFile"))
      use_cb_file=(*setup)["UseClockbreakFile"];

    if((vgosdb.does_file_exist("Session",cb_filename))&&use_cb_file)
    {

        int ncbr = vgosdb.get_scalar<short>("Session", cb_filename,"BRK_NUMB");
        log<WARNING>("*** Clock break(s) existent: #ncbr ") % ncbr;
	
        vector<string> break_stations = vgosdb.get_string("Session", cb_filename,"ClockBreakStationList");
        vector<double> break_epochs = vgosdb.get_vector<double>("Session", cb_filename,"ClockBreakEpoch");

        // lenght of epochs and of stations must be the same
        if(break_stations.size() != break_epochs.size() && break_stations.size() != ncbr && break_stations.size() != 1 )
             throw runtime_error( "void Session_inout::_read_vgosdb(): Inconsistency within ClockBreak.nc" );

        for(int i=0; i<ncbr; i++) {
	  string station;
	  ivg::Date startep,endep;
	  if (break_stations.size() != 1) 
	    station=break_stations.at(i).substr(0,8);
	  else
	    station=break_stations.at(0).substr(0,8);
	  startep=session_ptr->getStart();
	  endep=session_ptr->getEnd();
	  // If break epoch in JD, convert to MJD
	  if (break_epochs.at(i)>2400000.5)
	    break_epochs.at(i)-=2400000.5; 
	  if (startep.get_double_mjd()<break_epochs.at(i) && endep.get_double_mjd()>break_epochs.at(i))
	    breaks[ivg::paramtype::cbr].insert(std::pair<string,double>(remove_spaces_end(station), break_epochs.at(i)) );
	  else
	    std::cout <<"clock break out of range " <<endl;
	}
    }
    
    // check if a CPWLF break for atmospheric parameters is needed, e.g., due to too large observation gaps
    if((bool)(*(session_ptr->_handling)).exists("handling") && (bool)((*(session_ptr->_handling))["handling"]).exists("add_atbr"))
    {   
        for( int i=0; i<(*(session_ptr->_handling))["handling"]["add_atbr"].getLength(); ++i )
        {
            string station = (*(session_ptr->_handling))["handling"]["add_atbr"][i][0];
            string date_str = (*(session_ptr->_handling))["handling"]["add_atbr"][i][1];

            int y = std::stoi( date_str.substr( 0,4 ) );
            int m = std::stoi( date_str.substr( 5,2 ) );
            int d = std::stoi( date_str.substr( 8,2 ) );
            int h = std::stoi( date_str.substr( 11,2 ) );
            int min = std::stoi( date_str.substr( 14,2 ) );
            double sec = s2d(date_str.substr( 17,2 ) );
       
            breaks[ivg::paramtype::atbr].insert(std::pair<string,double>(station.substr(0,8), ivg::Date(y,m,d,h,min,sec).get_jd()) );
            
            log<INFO>("*** Insert CPWLF breaks for atmospheric parameters due to defintion in [handling][add_atbr]");
        }
    } 
    // in case of clock breaks, set them to the paramerter list
    if( breaks[ivg::paramtype::cbr].size() > 0 || breaks[ivg::paramtype::atbr].size() > 0 )
    {
       session_ptr->_param_list.set_breaks( breaks );   
       log<INFO>("*** Setting clock breaks to be parameterized [#") % breaks[ivg::paramtype::cbr].size() % "] and CPWLF breaks for atmospheric parameters [#" % breaks[ivg::paramtype::atbr].size() % "]";
    }

    // Check for Baseline clock offset that should be set up
    bool use_blc_file=false;
    if(setup->exists("UseBLClockFile"))
      use_blc_file=(*setup)["UseBLClockFile"];
    string blc_filename;
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::BLClock,ivg::band::X))
        blc_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::BLClock,ivg::band::X);
    else
        blc_filename = "BaselineClockSetup";   
    if((vgosdb.does_file_exist("Solve",blc_filename))&&use_blc_file)
    {
      
      vector<vector<vector<char>>> blc_char= vgosdb.get_vector_3d_data<char>("Solve", blc_filename,"BaselineClock");
      for (int i=0;i<blc_char.size();i++)
	{
	  std::string blsta1(blc_char[i][0].begin(),blc_char[i][0].end());
	  std::string blsta2(blc_char[i][1].begin(),blc_char[i][1].end());
	  
	  session_ptr->_param_list.add_bl_clock(remove_spaces_end(blsta1)+"-"+remove_spaces_end(blsta2));
	}
      
      //vector<vector<string>> blc_stations = vgosdb.get_vector_2d_data<string>("Solve", blc_filename,"BaselineClock");
      //	for ( int i=0;i<blc_stations.size();i++)
      //	  std::cout << (blc_stations.at(i))[0] << " -- " <<(blc_stations.at(i))[1]  << endl;
    }
    
    // get information which station-clock should be used as reference
    if(vgosdb.does_file_exist("Solve","ClockSetup") && vgosdb.does_variable_exist("Solve","ClockSetup","RefClockStationList"))
    {
       vector<string> refclock = vgosdb.get_string("Solve","ClockSetup","RefClockStationList"); 

       if(refclock.size() >= 1)
       {
           // in order to get correct name for e.g. "NRAO 3" -> "NRAO_3"
           string refclo = remove_spaces_end( refclock.at(0).substr(0,8) );
           replace_string_in_place(refclo, " ", "_" );

           log<INFO>("*** Reference clock in vgosDB (not necessarily used, can be overwritten by configfile): ") % refclo;

         if(!setup->exists("RefClockStationList"))
             setup->add("RefClockStationList",Setting::TypeString);

         (*setup)["RefClockStationList"] = refclo;

         if(refclock.size() > 1)
           log<WARNING>("!!!  More than one reference clock in vgosDB. Using first one.") % refclo;  
       }       
    }
    else
        log<WARNING>("!!! No information about ClockSetup (RefClockStationList) in vgosdb. Automatic setting of RefClock not possible. Must be done manually in configfile.");


    // // // // // // // // // // // // // // // // // // // // // // // // // //

    // (2844x1) (e.g. 1 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 ...)
    std::string obsxref_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::ObsCrossRef,session_ptr->_band_type) )
        obsxref_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::ObsCrossRef,session_ptr->_band_type);
    else
        obsxref_filename = "ObsCrossRef";
    
    vector<int> obs2scan = vgosdb.get_vector<int>("CrossReference",obsxref_filename,"Obs2Scan");

    // (2844x2)
    
    ivg::Matrix obs2baseline = vgosdb.get_matrix("CrossReference",obsxref_filename,"Obs2Baseline");

    // (660x7)(#nscans x #nstations)
    std::string staxref_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::StationCrossRef,session_ptr->_band_type) )
        staxref_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::StationCrossRef,session_ptr->_band_type);
    else
        staxref_filename = "StationCrossRef";
    ivg::Matrix scan2station = vgosdb.get_matrix("CrossReference",staxref_filename,"Scan2Station");

    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // now get all information for initializing scans
    // (2844x1) (e.g. 2013/06/11-16:51:23.0 1614+051 )
    vector<string> scan_names;
    std::string scanname_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::ScanName,session_ptr->_band_type) )
        scanname_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::ScanName,session_ptr->_band_type);
    else
        scanname_filename = "ScanName";
    
    if(vgosdb.does_file_exist("Scan",scanname_filename))
        scan_names = vgosdb.get_string("Scan",scanname_filename,"ScanNameFull");
    else
    {
        log<WARNING>("!!! SPECIAL VASCC2015 CASE: ScanName does not exist. Generating ScanNames based on TimeUTC and SourceList.");
        // we need to create a scan_names vector following this format:
        // 2013/06/11-16:51:23.0 1614+051
        ivg::Matrix time = vgosdb.get_matrix("Scan","TimeUTC","YMDHM");
        vector<int> seconds_vec = vgosdb.get_vector<int>("Scan","TimeUTC","Second");
        time.append_cols(ivg::Matrix(seconds_vec));

        for(int i=0; i<time.rows(); i++)
        {
            stringstream ss;
            ss << time(i,0) << "/";
            ss << setfill('0')<< setw(2) << time(i,1) << "/";
            ss << setfill('0')<< setw(2) << time(i,2) << "-";
            ss << setfill('0')<< setw(2) << time(i,3) << ":";
            ss << setfill('0')<< setw(2) << time(i,4) << ":";
            ss << setfill('0')<< setw(2) << time(i,5) << ".0";
            ss << " " << source_names.at(0);
            scan_names.push_back(ss.str());
        }       
        // in this special VASCC2015 case we need to add 1 because the station index start with zero
        obs2baseline = obs2baseline+1;
    }
   
    // in order to be able to use different epochs for each delay
    std::string utc_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::TimeUTC,session_ptr->_band_type) )
        utc_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::TimeUTC,session_ptr->_band_type);
    else
        utc_filename = "TimeUTC";
    
    ivg::Matrix time = vgosdb.get_matrix("Observables",utc_filename,"YMDHM");
    vector<double> seconds_vec = vgosdb.get_vector<double>("Observables",utc_filename,"Second");
    time.append_cols(ivg::Matrix(seconds_vec));

    vector<ivg::Date> obs_epochs;

    for(int i=0; i<time.rows(); i++){
        if( time(i,0) < 1900 ){
            if ( time(i,0) < 79 )
                time(i,0) += 2000;
            else
                time(i,0) += 1900;
        }
	//       obs_epochs.push_back(ivg::Date(time(i,0),time(i,1),time(i,2),time(i,3),time(i,4),time(i,5)));
	obs_epochs.push_back(ivg::Date(time(i,0),time(i,1),time(i,2),time(i,3),time(i,4),seconds_vec[i]));
    }
    double leap1=obs_epochs.at(1).get_leap_sec();
    if (obs_epochs.at(obs_epochs.size()-1).get_leap_sec()>leap1){
       for (int i=0;i<obs_epochs.size();i++)
		obs_epochs.at(i).add_secs(leap1-obs_epochs.at(i).get_leap_sec());
    }
    // create Obs-Baseline-Station-Statistic
    ivg::Matrix stats(scan2station.cols(),scan2station.cols(),0.0);
    for(int obs_cnt=0; obs_cnt<obs2baseline.rows(); obs_cnt++)
    {
         int idx_1 = (int)obs2baseline(obs_cnt,0) - 1; // minus 1 because vector begins with 1
         int idx_2 = (int)obs2baseline(obs_cnt,1) - 1; // minus 1 because vector begins with 1

         stats(idx_1,idx_2) += 1.0;
	 stats(idx_2,idx_1) += 1.0;
         stats(idx_1,idx_1) += 1.0;
         stats(idx_2,idx_2) += 1.0;
    }
    
// save original statistics to be able to compare to statistics with flagged observations
    ivg::Matrix orig_stats = stats;

    // need to store station dependent data (e.g. met, cablecal) in matrix to gain performance
    map< string, std::map< std::string, std::vector<double> > > aux_data;
    for(int i=0; i<station_names.size(); i++)
    {
        // check whether cable calibration and meteorological data are available for
        // a certain station; if not, set default values. If there are data
        // available, check whether cable calibration data should be excluded
        // manually (specified in config file)
        string cablecal_filename;
        if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::CalCable,station_names.at(i)))
           cablecal_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::CalCable,station_names.at(i));
        else
           cablecal_filename = "Cal-Cable";
        if(vgosdb.does_file_exist(station_names.at(i),cablecal_filename) && (*session_ptr->_setup).exists("modify_cable_cal") == false )
        {                      
            aux_data[station_names.at(i)]["Cal-Cable"] = vgosdb.get_vector<double>(station_names.at(i),cablecal_filename,"Cal-Cable");
        }
        else if(vgosdb.does_file_exist(station_names.at(i),cablecal_filename) && (*session_ptr->_setup).exists("modify_cable_cal"))
        {
             aux_data[station_names.at(i)]["Cal-Cable"] = vgosdb.get_vector<double>(station_names.at(i),cablecal_filename,"Cal-Cable");
             
             for(int j=0; j<(*session_ptr->_setup)["modify_cable_cal"].getLength(); j++)
             {
                 string sta = (*session_ptr->_setup)["modify_cable_cal"][j][0];

                 if( station_names.at(i) == sta )
                 {
                     aux_data[station_names.at(i)]["Cal-Cable"] = vector<double>(nobs, 0.0);
                     log<WARNING>("*** Cable cal data for station ") % station_names.at(i) % " manually set to 0.0";                    
                 }
             }            
        }
        else
        {
            aux_data[station_names.at(i)]["Cal-Cable"] = vector<double>(nobs, 0.0);
            log<WARNING>("*** No cable cal data available for station ") % station_names.at(i) % ". Setting cable cal to 0.0";
        }
	string feed_filename;
        if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::FeedRotation,station_names.at(i)))
           feed_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::FeedRotation,station_names.at(i));
        else
           feed_filename = "FeedRotation";
        if(vgosdb.does_file_exist(station_names.at(i),feed_filename))
        {
             aux_data[station_names.at(i)]["FeedRotation"] = vgosdb.get_vector<double>(station_names.at(i),feed_filename,"FeedRotation");
	}
	else {
	  log<WARNING>("*** No feed rotation data available for station ") % station_names.at(i) % ". Setting feed rotation to 0.0";
	   aux_data[station_names.at(i)]["FeedRotation"] = vector<double>(nobs, 0.0);
	}
        // only take meteorological data from Met.nc file if the file is existent
	string met_filename;
        if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::Met,station_names.at(i)))
           met_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::Met,station_names.at(i));
        else
           met_filename = "Met";
        if(vgosdb.does_file_exist(station_names.at(i),met_filename))
        {
             aux_data[station_names.at(i)]["TempC"] = vgosdb.get_vector<double>(station_names.at(i),met_filename,"TempC");
             aux_data[station_names.at(i)]["AtmPres"] = vgosdb.get_vector<double>(station_names.at(i),met_filename,"AtmPres");
             aux_data[station_names.at(i)]["RelHum"] = vgosdb.get_vector<double>(station_names.at(i),met_filename,"RelHum");
        }
        // if there are no meteorological data available, it is possible (at 
        // least for twin telescopes) to copy the data from another station       
        else if( ((*session_ptr->_setup)["troposphere"]).exists("copy_wx_met_data") 
                    && (bool)(*session_ptr->_setup)["troposphere"]["copy_wx_met_data"][0] )
        {
             std::string from_station = string((const char*)(*session_ptr->_setup)["troposphere"]["copy_wx_met_data"][2]);
             std::string to_station = string((const char*)(*session_ptr->_setup)["troposphere"]["copy_wx_met_data"][1]);            

             log<WARNING>("*** No meteorological data for station ") % station_names.at(i) % ". Copy data from station " % from_station ;              if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::Met,from_station))
                 met_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::Met,from_station);
             else
                 met_filename = "Met";
             aux_data[to_station]["TempC"] = vgosdb.get_vector<double>(from_station,met_filename,"TempC");
             aux_data[to_station]["AtmPres"] = vgosdb.get_vector<double>(from_station,met_filename,"AtmPres");
             aux_data[to_station]["RelHum"] = vgosdb.get_vector<double>(from_station,met_filename,"RelHum");              
        }
        // otherwise set values to "trash"-values like it it done in the ngs card 6
        else
        {
            log<WARNING>("*** No meteorological data for station ") % station_names.at(i) % ". Setting [TempC = -999.000, AtmPres = -999.000, RelHum = -99900.000]";
            aux_data[station_names.at(i)]["TempC"] = vector<double>(nobs, -999.000);
            aux_data[station_names.at(i)]["AtmPres"] = vector<double>(nobs, -999.000);
            aux_data[station_names.at(i)]["RelHum"] = vector<double>(nobs, -99900.000);
        }
    }

    // edit file to be loaded, e.g. Edit_IVG or Edit_iGSFC
    string edit_load_file = (*session_ptr->_setup)["outliers"]["load"][1];
    
    
    //read grouprate
    std::string groupRate_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::GroupRate,session_ptr->_band_type) )
        groupRate_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::GroupRate,session_ptr->_band_type);
    else
        groupRate_filename = "GroupRate_b" + band_str;
    
    std::vector<double> group_rate = vgosdb.get_vector<double>("Observables", groupRate_filename, "GroupRate"); // unit second
    std::vector<double> group_rate_sig = vgosdb.get_vector<double>("Observables", groupRate_filename, "GroupRateSig"); // unit second
    

    // (2844x1)
    vector<double> delay;
    // load simulated group delays
    // creating correct name for GroupDelayFull, e.g. GroupDelayFull_iIVS_bX.nc
    
    //read group delay
    string groupdelayfull;
    
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::GroupDelayFull,session_ptr->_band_type) )
        groupdelayfull = _wrapper_ptr->get_file(ivg::wrapper_entries::GroupDelayFull,session_ptr->_band_type);
    else 
    {// if no wrapper has been found, use_wrapper is false, or the wrapper has no the entrie for GroupDelayFull

         if(editing.empty())
             groupdelayfull = "GroupDelayFull_b" + band_str;
         else
             groupdelayfull = "GroupDelayFull_"+editing+"_b" + band_str;

         if(use_wrapper)
             log<WARNING>("!!! GroupDelayFull not found in wrapper. Using ") % groupdelayfull % ".nc instead";
    }
    std::string gd_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::GroupDelay,session_ptr->_band_type) )
        gd_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::GroupDelay,session_ptr->_band_type);
    else
        gd_filename = "GroupDelay_b"+band_str;
    
    if( (bool)(*session_ptr->_setup).exists("SIM") && (bool)(*session_ptr->_setup)["SIM"]["apply"] 
         && (bool)(*session_ptr->_setup)["SIM"]["load"][0] )
    {
        string edit_sim_load_file = (*session_ptr->_setup)["SIM"]["load"][1];
        delay = vgosdb.get_vector<double>("ObsEdit",edit_sim_load_file,"GroupDelayFull"); 
    }
     // standard case (second condition:  do not read group delay full for ambiguity resolution)
    else if(vgosdb.does_file_exist("ObsEdit",groupdelayfull) &&  (session_ptr->_ambigRes == false || session_ptr->_phaseDelay == true))
        delay = vgosdb.get_vector<double>("ObsEdit",groupdelayfull,"GroupDelayFull"); // unit second
    else
    {    
        if(session_ptr->_ambigRes)
            log<WARNING>("!!! Ambiguity resolution mode: Only GroupDelay (NOT GroupDelayFull) is used.");
        else
            log<WARNING>("!!! Ambiguities have not been resolved. Only GroupDelay (NOT GroupDelayFull) is used.");
        delay = vgosdb.get_vector<double>("Observables", gd_filename, "GroupDelay"); // unit second
    }
    
    vector<double> delay_sigma = vgosdb.get_vector<double>("Observables",gd_filename, "GroupDelaySig"); // unit second

    // read single band delay only if file exists
    vector<double> sb_delay, sb_delay_sigma;
     std::string sb_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::SBDelay,session_ptr->_band_type) )
        sb_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::SBDelay,session_ptr->_band_type);
    else
        sb_filename = "SBDelay"+band_str;
    if(vgosdb.does_file_exist("Observables",sb_filename))
    {
         sb_delay = vgosdb.get_vector<double>("Observables",sb_filename, "SBDelay"); // unit second
	 try
	   {
	     sb_delay_sigma = vgosdb.get_vector<double>("Observables",sb_filename, "SBDelaySig");// unit second
	   }
	 catch (exception& e) // in case there are no sigmas availible, set them to some value
	   {
	     sb_delay_sigma = vector<double>(sb_delay.size(),1.0e-6);
	   }
    }
    // read phases only if file exists
    vector<double> phase, phase_sigma;
     std::string phase_filename;
     
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::Phase,session_ptr->_band_type) )
        phase_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::Phase,session_ptr->_band_type);
    else
        phase_filename = "Phase_b"+band_str;
    
    if(vgosdb.does_file_exist("Observables",phase_filename))
    {
      
         phase = vgosdb.get_vector<double>("Observables",phase_filename, "Phase"); // unit second
	 try
	   {
	     phase_sigma = vgosdb.get_vector<double>("Observables",phase_filename, "PhaseSig");// unit second
	   }
	 catch (exception& e) // in case there are no sigmas availible, set them to some value
	   {
	     phase_sigma = vector<double>(phase.size(),2.0);
	   }
    } else {
      log<WARNING>("!!! No Phases available (/Observables/" +phase_filename+ ".nc not existent).");
    }
    
    // read phase delays only if file exists
    vector<double> phase_delay, phase_delay_sigma;
    std::string phase_delay_filename;
    
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::PhaseDelayFull,session_ptr->_band_type)  )
    {
        phase_delay_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::PhaseDelayFull,session_ptr->_band_type);
    }
    else
    {
        phase_delay_filename = "PhaseDelayFull_b"+band_str;
    }
    
    if(vgosdb.does_file_exist("ObsEdit",phase_delay_filename)&& !session_ptr->_ambigRes)
    {
      
	phase_delay = vgosdb.get_vector<double>("ObsEdit",phase_delay_filename, "PhaseDelayFull"); // unit second
	//phase_delay = vgosdb.get_vector<double>("ObsEdit",phase_delay_filename, "PhaseDelayFull"); // unit microsecond
	//for (int j=0;j<phase_delay.size();j++)
	  //phase_delay[j]*=1e-6;
	 try
	   {
	     phase_delay_sigma = vgosdb.get_vector<double>("ObsEdit",phase_delay_filename, "PhaseDelayFullSig");// unit second
	   }
	 catch (exception& e) // in case there are no sigmas availible, set them to some value
	   {
	     if (ref_freq.size()==1)
	       for (int j=0;j<phase_sigma.size();j++)
		 phase_delay_sigma.push_back(phase_sigma.at(j)/(2*M_PI*ref_freq.at(0)*1.0e6));
	     else {
	       for (int j=0;j<ref_freq.size();j++)
		 phase_delay_sigma.push_back(phase_sigma.at(j)/(2*M_PI*ref_freq.at(j)*1.0e6));
	     }
	   }
    } else if (vgosdb.does_file_exist("Observables",phase_filename)){ // Calulate phase delay from phases
      log<WARNING>("!!! No PhaseDelays found (/ObsEdit/" +phase_delay_filename+ ".nc not existent). Calculating from Phase.");
	std::string phase_numamb_filename;
	//	if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::NumPhaseAmbig,session_ptr->_band_type)  )
   	// {
      	// 	phase_numamb_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::NumPhaseAmbig,session_ptr->_band_type);
	//	 }
    //else
   	// {
	//   phase_numamb_filename = "NumPhaseAmbig_b"+band_str;
	//}
	//if(vgosdb.does_file_exist("ObsEdit",phase_numamb_filename)&& !session_ptr->_ambigRes)
	//{
	//	vector<short> ph_amb;
	//       ph_amb=vgosdb.get_vector<short>("ObsEdit",phase_numamb_filename, "NumPhaseAmbig"); // unit microsecond
       	      
	//	       if (ref_freq.size()==1) {
	//	for (int j=0;j<phase.size();j++) {
	    
	//	  phase_delay.push_back(phase.at(j)/(2*M_PI*ref_freq.at(0)*1.0e6));
	//	  phase_delay_sigma.push_back(phase_sigma.at(j)/(2*M_PI*ref_freq.at(0)*1.0e6));
	//	  phase_delay[j]+=round((delay[j]-phase_delay[j])*ref_freq.at(j)*1.0e6)/(ref_freq.at(j)*1.0e6);
		  //	  phase_delay[j]+=((double)ph_amb.at(j))/(ref_freq.at(j)*1.0e6);
	//	}
      	//	} else {


	//	    for (int j=0;j<ref_freq.size();j++) {
	//	
	//  		phase_delay.push_back(phase.at(j)/(2*M_PI*ref_freq.at(j)*1.0e6));
	//  		phase_delay_sigma.push_back(phase_sigma.at(j)/(2*M_PI*ref_freq.at(j)*1.0e6));
	//		phase_delay[j]+=round((delay[j]-phase_delay[j])*ref_freq.at(j)*1.0e6)/(ref_freq.at(j)*1.0e6);
			//	phase_delay[j]+=((double)ph_amb.at(j))/(ref_freq.at(j)*1.0e6);
	//	    }

	//     }	       
	//}
	//else
	//{
		try {
      		if (ref_freq.size()==1) {
		for (int j=0;j<phase.size();j++) {
		  phase_delay.push_back(phase.at(j)/(2*M_PI*ref_freq.at(0)*1.0e6));
		  phase_delay_sigma.push_back(phase_sigma.at(j)/(2*M_PI*ref_freq.at(0)*1.0e6));
		  phase_delay[j]+=round((delay[j]-phase_delay[j])*ref_freq.at(j)*1.0e6)/(ref_freq.at(j)*1.0e6);
		}
      		} else {


	  	    for (int j=0;j<ref_freq.size();j++) {
	  		phase_delay.push_back(phase.at(j)/(2*M_PI*ref_freq.at(j)*1.0e6));
	  		phase_delay_sigma.push_back(phase_sigma.at(j)/(2*M_PI*ref_freq.at(j)*1.0e6));
	  		phase_delay[j]+=round((delay[j]-phase_delay[j])*ref_freq.at(j)*1.0e6)/(ref_freq.at(j)*1.0e6);

		}
      

      		}
      		} catch (exception& e) {
			log<WARNING>("!!! Error in calculating phase delays from phases");
      		}
		//}

      } else {
           log<WARNING>("!!! No PhasesDelays found (/ObsEdit/" +phase_delay_filename+ ".nc not existent and no raw phases found (/Observables/" +phase_filename+ ".nc not existent). No phase delays availible!" );
    }
	 
    // information used to decide weather observation should be used or not
    vector<char> qcode;
    std::string qcode_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::QualityCode,session_ptr->_band_type) )
        qcode_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::QualityCode,session_ptr->_band_type);
    else
        qcode_filename = "QualityCode"+band_str;
    if(vgosdb.does_file_exist("Observables",qcode_filename))
         qcode = vgosdb.get_vector<char>("Observables",qcode_filename,"QualityCode");
     else
     {
        log<WARNING>("!!! No QualityCode information available (/Observables/QualityCode_b"+band_str+".nc not existent). Setting all observations to QualityCode = 8");
        qcode = vector<char>(nobs,'8');  
     }

    string CalSlantPathIonoGroup;
    
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::CalSlantPathIonoGroup,session_ptr->_band_type) )
        CalSlantPathIonoGroup = _wrapper_ptr->get_file(ivg::wrapper_entries::CalSlantPathIonoGroup,session_ptr->_band_type);
    else
    {
        CalSlantPathIonoGroup = "Cal-SlantPathIonoGroup_b" + band_str;
        if(use_wrapper)
             log<WARNING>("!!! Cal-SlantPathIonoGroup*.nc not found in wrapper. Using ") % CalSlantPathIonoGroup % ".nc instead";
    }


    vector<short> errFlag;

    if(vgosdb.does_variable_exist("ObsDerived", CalSlantPathIonoGroup, "Cal-SlantPathIonoGroupDataFlag"))
        errFlag = vgosdb.get_vector<short>("ObsDerived", CalSlantPathIonoGroup, "Cal-SlantPathIonoGroupDataFlag");
    else
    {
         log<WARNING>("!!! No ErrorFlag information available (variable /ObsDerived/" + CalSlantPathIonoGroup + ".nc/Cal-SlantPathIonoGroupDataFlag not existent). Setting all observations to ErrorFlag = 0");
         errFlag = vector<short>(nobs,0);  
    }
   
    ivg::Matrix ionCorr_matrix, ionCorrSigma_matrix;
    if(vgosdb.does_variable_exist("ObsDerived", CalSlantPathIonoGroup, "Cal-SlantPathIonoGroup"))
    {
        ionCorr_matrix = vgosdb.get_matrix("ObsDerived", CalSlantPathIonoGroup, "Cal-SlantPathIonoGroup");
        ionCorrSigma_matrix = vgosdb.get_matrix("ObsDerived", CalSlantPathIonoGroup, "Cal-SlantPathIonoGroupSigma");
	if (!vgosdb.does_variable_exist("ObsDerived", CalSlantPathIonoGroup, "Cal-SlantPathIonoGroupDataFlag"))
	  {
	    for (int j=0;j< ionCorrSigma_matrix.rows();j++)
	      {
		if (ionCorrSigma_matrix(j,0)==0) errFlag[j]=9;
	      }
	  }
    }
    else
    {
        log<WARNING>("!!! No IonoCorr information available (variable /ObsDerived/" + CalSlantPathIonoGroup + ".nc/Cal-SlantPathIonoGroup not existent). Setting all IonoCorr & IonoCorrSigma = 0");
        ionCorr_matrix = ivg::Matrix(nobs,2,0.0); 
        ionCorrSigma_matrix = ivg::Matrix(nobs,2,0.0); 
    }
    
    if( (bool)((*session_ptr->_setup)["ionosphere"]).exists("exclude_bls") && (*session_ptr->_setup)["ionosphere"]["apply"] )
    {
        for(int j=0; j<(*session_ptr->_setup)["ionosphere"]["exclude_bls"].getLength(); j++)
        {
            string bl = (*session_ptr->_setup)["ionosphere"]["exclude_bls"][j];
            log<WARNING>("*** Ionospheric correction disabled for baseline ") % bl % ".";
        }
    }
     std::string snr_filename;
        
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::SNR,ivg::band::X) )
        snr_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::SNR,ivg::band::X);
    else
        snr_filename = "SNR_bX";
    // try to get SNR information for each observation (X-BAND)
    vector<double> snr_bx(nobs,0.0);
    if(vgosdb.does_file_exist("Observables",snr_filename))
        snr_bx = vgosdb.get_vector<double>("Observables",snr_filename,"SNR"); // unit second
    else
        log<WARNING>("*** No SNR information available (/Observables/SNR_bX.nc not existent).");

     // try to get SNR information for each observation (S-BAND)
    vector<double> snr_bs(nobs,0.0);
    if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::SNR,ivg::band::S) )
        snr_filename = _wrapper_ptr->get_file(ivg::wrapper_entries::SNR,ivg::band::S);
    else
        snr_filename = "SNR_bS";
    if(vgosdb.does_file_exist("Observables",snr_filename))
      snr_bs = vgosdb.get_vector<double>("Observables",snr_filename,"SNR"); // unit second
    else
        log<WARNING>("*** No SNR information available (/Observables/SNR_bS.nc not existent).");

    // in case of presence of an edited Edit_IVG, use it!
    vector<short> delay_flag;
    if(vgosdb.does_file_exist("ObsEdit",edit_load_file) && (bool)(*session_ptr->_setup)["outliers"]["load"][0] )
    {
        log<WARNING>("*** Using ") % edit_load_file % " for DelayFlag";
        delay_flag = vgosdb.get_vector<short>("ObsEdit",edit_load_file,"DelayFlag");
    }
    else if((bool)(*session_ptr->_setup)["outliers"]["force_load"])
    {
        throw runtime_error("void Session_inout::_read_vgosdb(): setup[outliers][force_load] = true. Skipping session because "+edit_load_file+" not existent.");
    }
    else
    {

        string editname;
        
        if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::Edit, session_ptr->_band_type) )
            editname = _wrapper_ptr->get_file(ivg::wrapper_entries::Edit,session_ptr->_band_type);
        else
        {
             if(editing.empty())
                 editname = "Edit";
             else
                 editname = "Edit_"+editing;

             if(use_wrapper)
                 log<WARNING>("!!! Edit*.nc not found in wrapper. Using ") % editname % ".nc instead";
        }


        if(!(bool)(*session_ptr->_setup)["outliers"]["load"][0] && !vgosdb.does_file_exist("ObsEdit",edit_load_file))
            log<WARNING>("*** Using ") % editname % ".nc for DelayFlag because " % edit_load_file % " doesn't exist AND outliers_load == false";
        else if(!(bool)(*session_ptr->_setup)["outliers"]["load"][0])
            log<WARNING>("*** Using ") % editname % ".nc for DelayFlag because outliers load == false";
        else
            log<WARNING>("*** Using ") % editname % ".nc for DelayFlag because " % edit_load_file % " doesn't exist";
	if(vgosdb.does_file_exist("ObsEdit",editname))
	  delay_flag = vgosdb.get_vector<short>("ObsEdit",editname,"DelayFlag");
	else
	  delay_flag.assign(delay.size(),0);
    }

    if((bool)(*(session_ptr->_handling)).exists("handling") && (bool)((*(session_ptr->_handling))["handling"]).exists("ignore_cable_cal"))
        for( int i=0; i<(*(session_ptr->_handling))["handling"]["ignore_cable_cal"].getLength(); ++i )
            log<INFO>("*** Cable cal for station: ") % (const char*)(*(session_ptr->_handling))["handling"]["ignore_cable_cal"][i] % " set to zero.";
   
    log<RESULT>("*** Importing ") % scan_names.size() % " scans containing " % nobs % " nobs";

    map<string,int> info_cnter;
    
    if(obs_epochs.size() != delay.size())
        throw runtime_error("ERROR: Unequal length of epochs and delays ("+std::to_string(obs_epochs.size())+" vs "+std::to_string(delay.size())+")");
                
    int cnt=0;
    for(int scan_idx=0; scan_idx<scan_names.size(); scan_idx++){

        if(scan_idx % 100 == 0)
             log<INFO>("*** Scans imported: ") % scan_idx % "/" % scan_names.size();

        // scan_info looks like: 2013/06/11-16:51:23.0 1614+051 in case of regular vgosDB
        // scan_info looks like: 2013/06/11-16:51:23.04245 1614+051 in case of moon-specific-vgosDB
        string scan_info = scan_names.at(scan_idx);

        int y = std::stoi( scan_info.substr( 0,4 ) );
        int m = std::stoi( scan_info.substr( 5,2 ) );
        int d = std::stoi( scan_info.substr( 8,2 ) );
        int h = std::stoi( scan_info.substr( 11,2 ) );
        int min = std::stoi( scan_info.substr( 14,2 ) );

        // due to different amount of sec-fractions, it need to be tokenized by space
        vector<string> secs_and_src = get_tokens(scan_info.substr(17));
        double sec = s2d(secs_and_src.at(0));
        string srcname = secs_and_src.at(1);

        // epoch of scan
        ivg::Date epoch(y,m,d,h,min,sec);
	epoch=obs_epochs.at(cnt);
        // save epoch of first scan/observation as start-epoch
        if(scan_idx == 0)
            session_ptr->_start = epoch;

        // get source observerd in scan
        ivg::Source* source;
        session_ptr->_crf.get_source( &source, remove_spaces_end( srcname ) );

        //set crf2trf matrix for each scan
        ivg::Partials_t2c foo { ivg::Matrix( 3,3,0.0 ) };
        ivg::Partials_t2c * deriv_ptr = &foo;
        ivg::Matrix crf2trf = session_ptr->_eops.form_crf2trf( epoch, true, deriv_ptr );

        ivg::Scan scan(epoch, source, crf2trf.transpose(), deriv_ptr);
        // ATTENTION: duration only used in case of simulation based on vgsoDB
        // TODO: should be taken from .../Observables/CorrInfo-***_bX.nc
        scan.set_schedulded_duration(0.0);

        // now we need to handle every observation for current scan
        while(obs2scan.at(cnt) == (scan_idx+1) )
        {           
            // station index numbers within the vgosdb, not within the scan!!!
            double sta1_num = obs2baseline(cnt,0) - 1; // minus 1 because vector begins with 0
            double sta2_num = obs2baseline(cnt,1) - 1; // minus 1 because vector begins with 0

            string sta1 = station_names.at(sta1_num);
            string sta2 = station_names.at(sta2_num);

            int data_idx1 = scan2station(scan_idx,sta1_num)-1;
            int data_idx2 = scan2station(scan_idx,sta2_num)-1;

 //           cerr << scan_info << ": CNT: " << cnt << " SCAN_IDX: " << scan_idx << " | STA1: " << sta1_num << "/" << sta1 << "/" << data_idx1 << " | STA2: " << sta2_num << "/" << sta2 << "/" << data_idx2 << endl;

            // get Analysis_station pointer
            ivg::Analysis_station * station1;
            session_ptr->_trf.get_station( &station1, sta1 );

            ivg::Analysis_station * station2;
            session_ptr->_trf.get_station( &station2, sta2 );

            int sta1_idx = scan.add_sta( station1 );
            int sta2_idx = scan.add_sta( station2 );

            ivg::Obs obs_new( session_ptr, &scan, sta1_idx, sta2_idx );

            // overwriting epoch which was originally taken from &scan
            // this allows different epochs for each delay within a scan
            obs_new.set_epoch( obs_epochs.at(cnt) );


	    
            if( geocenter_delay.empty() )
                geocenter_delay = std::vector<double>(delay.size(),0.0);
            
            if(!group_rate.empty())
                obs_new.set_delay( delay.at(cnt), delay_sigma.at(cnt), geocenter_delay.at(cnt), group_rate.at(cnt), group_rate_sig.at(cnt) );
            else
                obs_new.set_delay( delay.at(cnt), delay_sigma.at(cnt), geocenter_delay.at(cnt) );
            
            // only set single-band delay if the information is available
            if(!sb_delay.empty())
                 obs_new.set_sb_delay( sb_delay.at(cnt), sb_delay_sigma.at(cnt) );
            obs_new.set_snr( snr_bx.at(cnt), snr_bs.at(cnt) );
            if(!phase_delay.empty())
	      obs_new.set_phase_delay( phase_delay.at(cnt), phase_delay_sigma.at(cnt) );

 //           cerr << "Delay_size: " << delay.size() << " " << " Obs_epoch_size: " << obs_epochs.size() << " ";
 //           cerr << "Delay: " << delay.at(cnt) << " SNX_BX: " << snr_bx.at(cnt) << endl;
 //           cerr << "IDX1: " << data_idx1 << " IDX2: " << data_idx2 << " LANG: " <<  aux_data[sta1]["Cal-Cable"].size() << " LANG2: " << aux_data[sta2]["Cal-Cable"].size()<< endl;
	    
	    double feedconv,feedrot1,feedrot2;
	    try {
	    if (ref_freq.size()==1)
	      feedconv=1/(2*M_PI*ref_freq.at(0)*1.0e6);
	    else
	      feedconv=1/(2*M_PI*ref_freq.at(cnt)*1.0e6);
	    
	    feedrot1 = aux_data[sta1]["FeedRotation"].at(data_idx1)*feedconv; // unit second
            feedrot2 = aux_data[sta2]["FeedRotation"].at(data_idx2)*feedconv; // unit second
	    } catch (exception& e) {
	      feedrot1=0;
	      feedrot2=0;
	    }
	    obs_new.set_feed_rotation( feedrot1, feedrot2 ); 
	    
            if(data_idx1 >= aux_data[sta1]["Cal-Cable"].size() || data_idx2 >= aux_data[sta2]["Cal-Cable"].size())
                throw runtime_error("void Session_inout::_read_vgosdb(....): Something wrong with the length of the cable-cal vector. Wrong wrapper version?");
            
            double cable_cal1 = aux_data[sta1]["Cal-Cable"].at(data_idx1); // unit second
            double cable_cal2 = aux_data[sta2]["Cal-Cable"].at(data_idx2); // unit second

            if((bool)(*(session_ptr->_handling)).exists("handling") && (bool)((*(session_ptr->_handling))["handling"]).exists("ignore_cable_cal"))
            {
                 for( int i=0; i<(*(session_ptr->_handling))["handling"]["ignore_cable_cal"].getLength(); ++i )
                 {
                     string station = (*(session_ptr->_handling))["handling"]["ignore_cable_cal"][i];                 
                     if( station == sta1 )
                         cable_cal1 = 0.0;
                     else if( station == sta2 )
                         cable_cal2 = 0.0;
                 }                                   
            }

            double temp1 = aux_data[sta1]["TempC"].at(data_idx1); // unit celsius
            double press1 = aux_data[sta1]["AtmPres"].at(data_idx1); // unit hPa
            double humidity1 = aux_data[sta1]["RelHum"].at(data_idx1); // unit %

            double temp2 = aux_data[sta2]["TempC"].at(data_idx2); // unit celsius
            double press2 = aux_data[sta2]["AtmPres"].at(data_idx2); // unit hPa
            double humidity2 = aux_data[sta2]["RelHum"].at(data_idx2); // unit %

            obs_new.set_cable_cal( cable_cal1, cable_cal2 );    

 //           cerr << sta1 << " (" << sta1_idx << "): " << cable_cal1 << "|" << temp1 << "|" << press1 << "|" << humidity1 << endl;
 //           cerr << sta2 << " (" << sta2_idx << "): " << cable_cal2 << "|" << temp2 << "|" << press2 << "|" << humidity2 << endl;

            std::string ext_data_type = string((const char*)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][1]);
            std::string ext_met_data = string((const char*)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][2]);
            std::string gpt2filename = string((const char*)(*session_ptr->_setup)["troposphere"]["gpt2_grid_file"]);
	    std::string gpt3filename = string((const char*)(*session_ptr->_setup)["troposphere"]["gpt3_grid_file"]);
	    std::string grad_type = string((const char*)(*session_ptr->_setup)["troposphere"]["gradient"]);
            if( !(bool)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][0] )           
            {
                ext_data_type = "meteo";
                ext_met_data = "insitu";
            }
            scan.add_scan_meteorology( sta1_idx, temp1, press1, humidity1, 0, ext_data_type, ext_met_data, gpt2filename, gpt3filename, grad_type );
            scan.add_scan_meteorology( sta2_idx, temp2, press2, humidity2, 0, ext_data_type, ext_met_data, gpt2filename, gpt3filename,grad_type);


            int qcode_num = qcode.at(cnt) - '0'; // e.g. 8

            if( qcode_num >= (int)(*session_ptr->_setup)["quality_code"])
                obs_new.set_use_flag( true );
            else
                info_cnter["qCode"]++;

            // setting ionospheric correction depending on ErrorFlag 0=OK , -1 = Missing, -2 = bad
            if( errFlag.at(cnt) == 0 )
            {
                // check whether the ionospheric correction should be disabled
                // for individual sessions
                if( (bool)((*session_ptr->_setup)["ionosphere"]).exists("exclude_bls") 
                     && (*session_ptr->_setup)["ionosphere"]["exclude_bls"].getLength()>0 )
                {
                    std::vector<std::string> bl;
                    for(int j=0; j<(*session_ptr->_setup)["ionosphere"]["exclude_bls"].getLength(); j++)
                        bl.push_back( (*session_ptr->_setup)["ionosphere"]["exclude_bls"][j] );
                    
                    std::vector<std::string>::iterator it = find( bl.begin(), bl.end(), sta1+"-"+sta2 );
                    
                    if( it != bl.end() )
                       obs_new.set_ion_corr( 0.0, 0.0, 0.0, 0.0 );
                    else               
                       obs_new.set_ion_corr( ionCorr_matrix(cnt,0), ionCorrSigma_matrix(cnt,0), ionCorr_matrix(cnt,1), ionCorrSigma_matrix(cnt,1) );                   
                }
                else
                {
                    obs_new.set_ion_corr( ionCorr_matrix(cnt,0), ionCorrSigma_matrix(cnt,0), ionCorr_matrix(cnt,1), ionCorrSigma_matrix(cnt,1) );
                }
            }
            else
            {
                obs_new.set_ion_corr( 666.66 , 1e-12, 0.0 , 0.0 );
                if((bool)(*session_ptr->_setup)["ionosphere"]["apply"])
                    obs_new.set_use_flag( false ); // if ionosphere correction is not applied use_flag is true
                info_cnter["Iono"]++;
            }

            if( (bool)(*session_ptr->_setup)["use_obs_flags"] && obs_new.get_use_flag() && delay_flag.at(cnt) != 0) 
                obs_new.set_use_flag( false );

            //  saving of outliers only possible if all observations are used
            //if( (bool)(*session_ptr->_setup)["outliers"]["save"][0] )
            //    obs_new.set_use_flag( true );

            // increment counter
            if(delay_flag.at(cnt) != 0) 
                info_cnter["Delay"]++;

            // check weather observation got bad iono error flag AND bad delay flag
            if(errFlag.at(cnt) != 0 && delay_flag.at(cnt) != 0)
                info_cnter["BothIonoDelay"]++;

            scan.add_obs( obs_new );

            // clear statistic from not useable (flagged) observations
            if(obs_new.get_use_flag() == false)
            {
                stats(sta1_num,sta2_num) -= 1.0;
		stats(sta2_num,sta1_num) -= 1.0;
                stats(sta1_num,sta1_num) -= 1.0;
                stats(sta2_num,sta2_num) -= 1.0;
            }

	    if (obs_new.get_use_flag() ==true)
	      {
		station1->inc_num_obs();
		station2->inc_num_obs();
		if (station1->get_num_obs()==1)
		  station1->set_first_epoch(obs_epochs.at(cnt));
		if (station2->get_num_obs()==1)
		  station2->set_first_epoch(obs_epochs.at(cnt));
		station1->set_last_epoch(obs_epochs.at(cnt));
		station2->set_last_epoch(obs_epochs.at(cnt));
	      }

            cnt++;  

            //stop importing if #nobs is reached
            if(cnt == nobs)
                break;
        }
        
        //add scan to scan vector
        session_ptr->_scans.push_back( scan );

         if(epoch.get_decimal_date() > session_ptr->_end.get_decimal_date())
         {
            log<INFO>("*** Observation ") % cnt % "/" % nobs % " later than session-end from Edit***.nc. Setting new _end to " % epoch.get_date_time("DD/MON/YYYY HH:MI:SS");
            session_ptr->_end = epoch;
         }
    }
   
    log<RESULT>("*** Observations with invalid qcode: ") % info_cnter["qCode"] % " of " % nobs;
    log<RESULT>("*** Observations with invalid delay flag: ") % info_cnter["Delay"] % " of " % nobs;
    log<RESULT>("*** Observations with invalid ionospheric correction: ") % info_cnter["Iono"] % " of " % nobs;
    log<RESULT>("*** Observations with invalid ionospheric correction AND invalid delay flag: ") % info_cnter["BothIonoDelay"] % " of " % nobs;
    
    // calculate percentage of outliers to be expected
    session_ptr->_exp_perc_out = ((double) info_cnter["Delay"] / (double)nobs )*100.0;
    
    if( (bool)(*session_ptr->_setup)["outliers"]["save"][0] )
        log<WARNING>("*** ATTENTION: Overwriting quality code and use flag. Using ANY observation because of outliers{ save = (true, ...) }");
    
    /////////////////////////////////////////////////////////////////////////////////
    // generating information table containing information about flagged observations
    ivg::Matrix perc_stats = (stats.div_elem(orig_stats)+(-1.0))*(-100.0);    
    stringstream rop; // Flagged Observations Proportion
    rop << "*** Percentage of flagged observations and remaining number of observations for each site and baseline" << endl;
    rop << "***             ";
    for(int sta_idx=0; sta_idx<station_names.size(); sta_idx++)
        rop << "|" << setw(9) << setfill(' ') << right << station_names.at(sta_idx);
    rop << "|" << endl;
    // loop over each station and each baseline and get the percentage of flagged obs
    for(int sta_idx=0; sta_idx<station_names.size(); sta_idx++)
    {
        rop << "*** " << setw(8) << setfill(' ') << station_names.at(sta_idx) << "| % |";
        for(int col=0; col<perc_stats.cols(); col++)
        {
            if(col<sta_idx)
                rop << "         |";
            else if(orig_stats(sta_idx,col) == 0.0 || stats(sta_idx,col) == 0.0 )
                rop << "   ---   |";
            else if(orig_stats(sta_idx,col) ==  stats(sta_idx,col) )
                rop << "   0.0 % |";
            else
                rop << setw(6) << setfill(' ') << setprecision(1) << fixed << perc_stats(sta_idx,col) << " % |";
        }
        rop << endl << "*** " << setw(8) << setfill(' ') << station_names.at(sta_idx) << "| # |";
        
        for(int col=0; col<perc_stats.cols(); col++)
        {
            if(col<sta_idx)
                rop << "         |";
            else
                rop << setw(6) << setfill(' ') << setprecision(0) << fixed << stats(sta_idx,col) << " # |";
        }
        
        if(sta_idx < station_names.size()-1)
            rop << endl;
    }
    log<RESULT>(rop.str());
    
    // checking for stations without too few observatios defined in config-file (elim_min_sta)
    int elim_min_sta = (*session_ptr->_setup)["elim_min_sta"];
    vector<int> zero_obs = stats.diag().find_idx(le,elim_min_sta);
    if(zero_obs.size()>0)
    {
        for(int i=0; i<zero_obs.size(); i++)
        {
            string station = station_names.at(zero_obs.at(i));
            log<WARNING>("!!! Station ") % station % " less than " % elim_min_sta % " observations. Automatically setting station temporary to [handling][elim_sta].";
            
            if((bool)(*session_ptr->_handling).exists("handling") && (bool)((*session_ptr->_handling)["handling"]).exists("elim_sta"))
                (*session_ptr->_handling)["handling"]["elim_sta"].add(Setting::TypeString) = station;
            else
            { 
                if(!(bool)(*session_ptr->_handling).exists("handling"))
                    Setting &handling = (*session_ptr->_handling).add("handling",Setting::TypeGroup);
                
                if((bool)(*session_ptr->_handling).exists("handling") && !(bool)((*session_ptr->_handling)["handling"]).exists("elim_sta"))
                {
                    Setting &elim_sta = (*session_ptr->_handling)["handling"].add("elim_sta",Setting::TypeList);
                    (*session_ptr->_handling)["handling"]["elim_sta"].add(Setting::TypeString) = station;
                }
            }
        }
    }
    if ((bool)(*session_ptr->_setup).exists("elim_min_baseline"))

	{
	  int elim_min_bl = (*session_ptr->_setup)["elim_min_baseline"];
	  
	  for (int i=0;i<station_names.size();i++)
	    {
	      for (int j=i+1;j<station_names.size();j++)
		{
		  if (stats(i,j)<=elim_min_bl)
		    {
		      string bline= station_names.at(i)+"-"+station_names.at(j);
		       log<WARNING>("!!! Baseline ") % bline % " less than " % elim_min_bl % " observations. Automatically setting baseline temporary to [handling][elim_baseline].";
		      if((bool)(*session_ptr->_handling).exists("handling") && (bool)((*session_ptr->_handling)["handling"]).exists("elim_baseline"))
			(*session_ptr->_handling)["handling"]["elim_baseline"].add(Setting::TypeString) = bline;
		      else
			{ 
			  if(!(bool)(*session_ptr->_handling).exists("handling"))
			    Setting &handling = (*session_ptr->_handling).add("handling",Setting::TypeGroup);
			  
			  if((bool)(*session_ptr->_handling).exists("handling") && !(bool)((*session_ptr->_handling)["handling"]).exists("elim_baseline"))
			    {
			      Setting &elim_sta = (*session_ptr->_handling)["handling"].add("elim_baseline",Setting::TypeList);
			      (*session_ptr->_handling)["handling"]["elim_baseline"].add(Setting::TypeString) = bline;
			    }
			}  
		    }
		}
	    }
	}
	
    
    /////////////////////////////////////////////////////////////////////////////////
    // save statistics in session
    session_ptr->_obsstats["TRF_REM"] = orig_stats-stats; // stats: removed/flagged observations
    session_ptr->_obsstats["TRF_LEF"] = stats; // stats: observations left after flagging/removing
    session_ptr->_obsstats["TRF_ORI"] = orig_stats; // original stats
    session_ptr->_obsstats["TRF_PER"] = perc_stats; // percentage of flagged/removed
    /////////////////////////////////////////////////////////////////////////////////
    
    // correct epoch ?!?!? or from first and last scan???
    session_ptr->_param_list.set_start_end_epoch(session_ptr->_start,session_ptr->_end);
   
    log<NOTHING>("*** END-Loading vgosDB: ") % directory;
   
#if DEBUG_VLBI >=2
   cerr << "--- Session_inout::_read_vgosdb(ivg::Session *session_ptr, Setting *setup, const string directory)" 
        << " : " << tim.toc() << " s " << endl; 
#endif
}

// ...........................................................................
void Session_inout::write_residuals(ivg::Session *session_ptr,string outfile,std::string resid_format)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ void Session_inout::write_residuals(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    ivg::Matrix resid = session_ptr->_lsa_solution.get_resid()*1e12;

    ofstream out(outfile.c_str());

    if (!out.is_open())
        throw runtime_error("void Session_inout::write_residuals(ivg::Session *, string ): Failed to open file for writing: "+outfile);

    double mjd;
    std::string sta1,sta2,src;
    ivg::Matrix AzEl1, AzEl2;
    int counter = 0;
    for (int i = 0; i<=session_ptr->_scans.size()-1; i++)
    {
        for (int j = 0; j<=session_ptr->_scans.at(i).get_nobs()-1; j++)
        {
            mjd = session_ptr->_scans.at(i).get_obs_ptr(j)->get_epoch().get_double_mjd();

            session_ptr->_scans.at(i).get_obs_ptr(j)->get_station_names(sta1,sta2);
            session_ptr->_scans.at(i).get_obs_ptr(j)->get_source_name(src);

            // write output
            if( resid_format == "ASCOT" )
            {
                // standard ASCOT residual file 
                // ( ID - MJD - source name - name station1 - name station2 - residual )
                out<<setiosflags(ios::right)<<setiosflags(ios::fixed)
                        <<setw(5)<<counter<<"  "
                        <<setprecision(16)<<setw(24)<<mjd<<" ";
                out<<left<<setw(8)<<src<<" "
                        <<setw(8)<<sta1<<" "
                        <<setw(15)<<sta2<<" "
                        <<setprecision(16)<<setw(24)<<resid(counter)<<endl;
            }
            else if( resid_format == "ASCOT_EXT" )
            {
                // extended ASCOT residual file 
                // ( ID - MJD - source name - name station1 - name station2 -residual -
                //   azimuth station1 - elevation station1 - azimuth station2 - elevation station2 )
                ivg::Obs *obs_ptr = session_ptr->_scans.at(i).get_obs_ptr(j);
                ivg::Analysis_station * sta_ptr1 = session_ptr->_scans.at(i).get_sta_ptr( obs_ptr->get_scan_idx(1) );
                ivg::Analysis_station * sta_ptr2 = session_ptr->_scans.at(i).get_sta_ptr( obs_ptr->get_scan_idx(2) );                

                AzEl1 = sta_ptr1->calc_az_el( session_ptr->_scans.at(i).get_obs_ptr(j)->get_epoch(), 
                                              session_ptr->_scans.at(i).get_source()->get_unit_vector_ssb(), 
                                              session_ptr->_scans.at(i).get_trf2crf_matrix().transpose() );
                AzEl2 = sta_ptr2->calc_az_el( session_ptr->_scans.at(i).get_obs_ptr(j)->get_epoch(), 
                                              session_ptr->_scans.at(i).get_source()->get_unit_vector_ssb(), 
                                              session_ptr->_scans.at(i).get_trf2crf_matrix().transpose() );                
                
                out <<  setiosflags(ios::right) << setiosflags(ios::fixed)
                            << setw(5) << counter << "  "
                            << setprecision(16) << setw(24) << mjd << " ";
                out << left << setw(8) << src << " "
                            << setw(8) << sta1 <<" "
                            << setw(15) <<sta2 <<" "
                            << setprecision(16)<<setw(24)<<resid(counter)<<" "
                            << setprecision(16) << setw(24) << AzEl1(0) << "  "
                            << setprecision(16) << setw(24) << AzEl2(0) << "  "
                            << setprecision(16) << setw(24) << AzEl1(1) << "  "
                            << setprecision(16) << setw(24) << AzEl2(1)
                            << endl;                
            }    
            else if( resid_format == "SUB_AMB" )
            {
                ivg::Date d = session_ptr->_scans.at(i).get_obs_ptr(j)->get_epoch();
                out<<"         "
                   <<setiosflags(ios::right)<<setiosflags(ios::fixed)
                   <<left<<setw(8)<<sta1<<" "<<setw(8)<<sta2<<" "<<setw(8)<<src
                   <<" "<<right<<setw(4)<<d.get_int_year()
                   <<" "<<setw(2)<<d.get_int_month()
                   <<" "<<setw(2)<<d.get_int_day()
                   <<" "<<setw(2)<<d.get_int_hour()
                   <<" "<<setw(2)<<d.get_int_min()
                   <<" "<<setw(4)<<setprecision(1)<<d.get_double_sec()
                   <<" "<<setprecision(0)<<setw(15)
                   <<session_ptr->_scans.at(i).get_obs_ptr(j)->get_group_delay()*1e12
                   <<" "<<setprecision(0)<<setw(8)<<resid(counter)
                   <<endl;
            }                     
            counter++;
        }
    }
    out.close();


#if DEBUG_VLBI >=2
    cerr<<"--- void Session_inout::write_residuals(ivg::Session *, string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Session_inout::write_eop_file(ivg::Session *session_ptr, string outfile, string sess_name)
// ...........................................................................
{
  int idxp=std::find(ivg::paramtype_str.begin(),ivg::paramtype_str.end(),"xpo")-ivg::paramtype_str.begin();
  ivg::Matrix data,sigdata,colnum;
  Setting *setup=session_ptr->get_setup();
  ivg::Matrix vcm= session_ptr->_solution->get_vcm();
  vector<ivg::Date> epochs;
  ivg::Date sess_st=session_ptr->_start;
  ivg::Date sess_en=session_ptr->_end;
  if (sess_name=="DUMMY") 
    sess_name=session_ptr->_name_in_masterfile;
  vector<double> parfactor={0,0,0,0,0};
  int parid=0;
  int num_old=0;

  vector<string> techniques;
  if ((*setup)["export_eop"].exists("technique"))
    {
      stringstream ss((*setup)["export_eop"]["technique"]);
      string tmp;
      while (getline(ss,tmp,'+'))
	techniques.push_back(tmp);
    }
  
  vector<string> oldline;
  vector<ivg::Date> olddate;
  map<string,bool> oldeoptypes;
  map<string,string> oldeopstrings;
  if ((bool) (*setup)["export_eop"]["append"])
    {
      ifstream instream(outfile.c_str(),ios::in);
      if (instream.is_open())
	{
	  bool datafield=false, headerfield=false;
	  string str;
       
	  while (getline(instream,str))
	    {
	      if (str.substr(0,5)=="+DATA")
		datafield=true;
	      else if (str.substr(0,5)=="-DATA")
		datafield=false;
	      else if (str.substr(0,7)=="+HEADER")
		headerfield=true;
	      else if (str.substr(0,7)=="-HEADER")
		headerfield=false;
	      else if (datafield && str.substr(0,1)!="#" && str.substr(0,1)!="*" && str.substr(0,1)!="!")
		{
		  stringstream ss(str);
		  double tmpmjd;
		  string cursess,tmp;
		  ss >> tmpmjd >> tmp >> tmp >> tmp>> tmp>> tmp>> tmp >> tmp >> tmp>> tmp>> tmp>> tmp >> tmp >> tmp>> tmp>> tmp>> tmp >> cursess;
		  if (cursess!=sess_name)
		    {
		      num_old++;
		      olddate.push_back(ivg::Date(tmpmjd));
		      oldline.push_back(str);
		    }
		  
		}
	      else if (headerfield && str.substr(0,1)!="#" && str.substr(0,1)!="*" && str.substr(0,1)!="!")
		{
		  stringstream ss(str);
		  string tmp1,tmp2;
		  ss >> tmp1 >> tmp2;
		  if (tmp1=="DATA_START")
		    {
		      ivg::Date tmpd(stoi(tmp2.substr(0,4)),stoi(tmp2.substr(5,2)),stoi(tmp2.substr(8,2)),stoi(tmp2.substr(11,2)),stoi(tmp2.substr(14,2)),stod(tmp2.substr(17,5)));
		      if (tmpd<sess_st)
			sess_st=tmpd;
		      
		    }
		  if (tmp1=="DATA_END")
		    {
		      ivg::Date tmpd(stoi(tmp2.substr(0,4)),stoi(tmp2.substr(5,2)),stoi(tmp2.substr(8,2)),stoi(tmp2.substr(11,2)),stoi(tmp2.substr(14,2)),stod(tmp2.substr(17,5)));
		      if (tmpd>sess_en)
			sess_en=tmpd;
		      
		    }
		  if (tmp1=="EOP_ESTIMATED")
		    {
		      
		      oldeoptypes[tmp2]=true;
		      oldeopstrings[tmp2]=str;
		      
		    }
		  if (tmp1=="TECHNIQUE")
		    {
		     
		      stringstream ss(tmp2);
		      string tmp;
		      while (getline(ss,tmp,'+'))
			{
			  if (std::find(techniques.begin(),techniques.end(),tmp) == techniques.end())
			    techniques.push_back(tmp);
			}
		    }
		}
	    }
	    instream.close();

	    
	}
     
    }
  
  for (auto &para : session_ptr->_param_list)
    {
      
      if ((para.get_typename()=="xpo")||(para.get_typename()=="ypo")||(para.get_typename()=="ut1")||(para.get_typename()=="nutx")||(para.get_typename()=="nuty"))
	{
	  ivg::Date curt=para.get_epoch();
	  vector<ivg::Date>::iterator did=std::find(epochs.begin(), epochs.end(), curt);
	  int id=did-epochs.begin();
	  if (did==epochs.end())
	    {
	      
	      epochs.push_back(curt);
	      data.append_rows(ivg::Matrix(1,10,-999999));//vector<double>({-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999}));
	      sigdata.append_rows(ivg::Matrix(1,10,-999999));//v(vector<double>({-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999}));
	      colnum.append_rows(ivg::Matrix(1,10,-999999));//v(vector<double>({-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999}));
	      id=epochs.size()-1;
	    }
	 
	  
	  int parmnum=para.get_type()-idxp;
	  									    
	  int order=para.get_order();
	  
	  if (order<=1) {
	    data(id,parmnum+order*5)=(para.get_apriori()+para.get_estimate())*ivg::param_unit_fac.at(para.get_type());
	    sigdata(id,parmnum+order*5)=para.get_standard_deviation()*ivg::param_unit_fac.at(para.get_type());
	    colnum(id,parmnum+order*5)=parid;
	  }
	  if (parmnum<=2)
	    parfactor[parmnum]=ivg::param_unit_fac.at(para.get_type())/1000;
	  else
	    parfactor[parmnum]=ivg::param_unit_fac.at(para.get_type());	
	}
      parid++;
    }
  log<INFO>("*** Writing eop to ")%outfile;
  ofstream outstream(outfile.c_str(),ios::out);
  if (!outstream.is_open())
        throw runtime_error("void Session_inout::write_eop_file(ivg::Session *,string): Failed to open file for writing: "+outfile);
  ivg::Date curdate;
  curdate.now();
  vector<ivg::Date> sortepochs=epochs;
  std::sort(sortepochs.begin(),sortepochs.end());
  
  string agency="OSO";
  if  ((*setup)["export_eop"].exists("agency"))
    agency=(const char*)(*setup)["export_eop"]["agency"];
  outstream <<"%=IVS-EOP 3.0 "<< agency << " " << curdate.get_date_time("YYYY-MO-DDTHH:MI:SS") << " " << agency << " " << sess_st.get_date_time("YYYY-MO-DDTHH:MI:SS") << " " << sess_en.get_date_time("YYYY-MO-DDTHH:MI:SS") << " TAI R" << endl;
  outstream<< "+HEADER"<< endl;
  outstream<< setw(20) << left << "GENERATION_TIME" << curdate.get_date_time("YYYY-MO-DDTHH:MI:SS")<< endl;
  outstream<< setw(20) << left << "DATA_START" << sess_st.get_date_time("YYYY-MO-DDTHH:MI:SS")<< endl;
  outstream<< setw(20) << left << "DATA_END" << sess_en.get_date_time("YYYY-MO-DDTHH:MI:SS")<< endl;													 
  if ((*setup)["export_eop"].exists("description"))
    outstream<< setw(20) << left << "DESCRIPTION" << (const char*)((*setup)["export_eop"]["description"])<<endl;
  else
    outstream<< setw(20) << left << "DESCRIPTION" << "EOP estiamted using ASCOT"<<endl;
  if ((*setup)["export_eop"].exists("analysis_center"))
    outstream<< setw(20) << left << "ANALYSIS_CENTER" << (const char*)(*setup)["export_eop"]["analysis_center"]<<endl;
  else
    outstream<< setw(20) << left << "ANALYSIS_CENTER" << "An analysis center using ASCOT"<<endl;
  if ((*setup)["export_eop"].exists("contact"))
    outstream<< setw(20) << left << "CONTACT" << (const char*)(*setup)["export_eop"]["contact"]<<endl;
  else
    outstream<< setw(20) << left << "CONTACT" << "ASCOT user <name@ac.xx>"<<endl;
  if ((*setup)["export_eop"].exists("software"))
    outstream<< setw(20) << left << "SOFTWARE" << (const char*)(*setup)["export_eop"]["software"]<<endl;													 
  else
    outstream<< setw(20) << left << "SOFTWARE" << "ASCOT"<<endl;	
    if (techniques.size()>0)
      {
	outstream<< setw(20) << left << "TECHNIQUE" << techniques.at(0);
	for (int j=1;j<techniques.size();j++)
	  outstream<< "+" <<techniques.at(j);
	outstream<<endl;
      }
  else
    outstream<< setw(20) << left << "TECHNIQUE" << "VLBI"<<endl;
  if ((*setup)["export_eop"].exists("nutation_type"))
    outstream<< setw(20) << left << "NUTATION_TYPE" << (const char*)(*setup)["export_eop"]["nutation_type"]<<endl;													 
  else
    outstream<< setw(20) << left << "NUTATION_TYPE" << "CIO-BASED"<<endl;
   if ((*setup)["export_eop"].exists("rotation_type"))
    outstream<< setw(20) << left << "ROTATION_TYPE" << (const char*)(*setup)["export_eop"]["rottion_type"]<<endl;													 
  else
    outstream<< setw(20) << left << "ROTATION_TYPE" << "UT1-TAI_LOD"<<endl;									  
  if ((*setup)["export_eop"].exists("TRF"))
    outstream<< setw(20) << left << "TRF_APRIORI" << (const char*)(*setup)["export_eop"]["TRF"]<<endl;													 
  else
    outstream<< setw(20) << left << "TRF_APRIORI" << (const char*)(*setup)["trf"]<<endl;	
  if ((*setup)["export_eop"].exists("CRF"))
    outstream<< setw(20) << left << "CRF_APRIORI" << (const char*)(*setup)["export_eop"]["CRF"]<<endl;													 
  else
    outstream<< setw(20) << left << "CRF_APRIORI" << (const char*)(*setup)["crf"]<<endl;	
  if ((*setup)["export_eop"].exists("eop_subdaily"))
    outstream<< setw(20) << left << "EOP_SUB-DAILY_MODEL" << (const char*)(*setup)["export_eop"]["eop_subdaily"]<<endl;													 
  else
    outstream<< setw(20) << left << "EOP_SUB-DAILY_MODEL" << (const char*)(*setup)["eop"]["hf_ocean_model"]<<endl;	
  if ((*setup)["export_eop"].exists("eop_apriori"))
    outstream<< setw(20) << left << "EOP_APRIORI" << (const char*)(*setup)["export_eop"]["eop_apriori"]<<endl;													 
  else
    outstream<< setw(20) << left << "EOP_APRIORI" << (const char*)(*setup)["eop"]["erp_aprioris"]<<endl;

  vector<string> eop_comp={"XPOL","YPOL","DUT1","DX","DY"};
   vector<string> eop_der={"XPOL_DER_1","YPOL_DER_1","LOD","DX_DER_1","DY_DER_1"};
  vector<string> eop_unit={"as","as","s","mas","mas"};
  vector<string> eop_conf={"pm","pm","ut1","nut","nut"};
  for (int i=0;i<5;i++)
   {
     string handling=(*setup)["PARAMS"].lookup(eop_conf.at(i))[0]["handling"];
     if (handling=="none" )
       {
	   if ((bool)(*setup)["PARAMS"].lookup(eop_conf[i])[0]["cpwlf"]["insert"])
	   {
	     outstream<< setw(20) << left << "EOP_ESTIMATED" << setw(11) << eop_comp.at(i)+"_BSP_1";
	     if (((double)(*setup)["PARAMS"].lookup(eop_conf.at(i))[0]["cpwlf"]["rate_cnstr"])==0)
	       outstream<< setw(7) << right << "NONE" << " " << eop_unit.at(i) <<endl;
	     else
		 outstream<< fixed <<setw(7) << right << setprecision(4)<<((double)(*setup)["PARAMS"].lookup(eop_conf[i])[0]["cpwlf"]["rate_cnstr"])*parfactor.at(i) << " " << eop_unit.at(i) <<endl;
	   }
	 else
	   {
	       outstream<< setw(20) << left << "EOP_ESTIMATED" << setw(11) << eop_comp.at(i);
	       if (((double)(*setup)["PARAMS"].lookup(eop_conf[i])[0]["polynom"]["cnstr"][0])==0)
	         outstream<< setw(7) << right << "NONE" << " " << eop_unit.at(i) <<endl;
	       else
		   outstream<< fixed <<setw(7) << right << setprecision(4)<< ((double)(*setup)["PARAMS"].lookup(eop_conf[i])[0]["polynom"]["cnstr"][0])*parfactor.at(i) << " " << eop_unit.at(i) <<endl;
	   }
	   if (((int)(*setup)["PARAMS"].lookup(eop_conf[i])[0]["polynom"]["order"])>=1)
	     {
	       outstream<< setw(20) << left << "EOP_ESTIMATED" << setw(11) << eop_der.at(i);
		 if (((double)(*setup)["PARAMS"].lookup(eop_conf[i])[0]["polynom"]["cnstr"][1]==0))
	         outstream<< setw(7) << right << "NONE" << " " << eop_unit.at(i) <<endl;
	       else
		   outstream<< fixed <<setw(7) << right << setprecision(4)<< ((double)((*setup)["PARAMS"].lookup(eop_conf[i])[0]["polynom"]["cnstr"][1]))*parfactor.at(i) << " " << eop_unit.at(i) << "/day" <<endl;								      
	     }
	 else if (oldeoptypes[eop_der.at(i)])
	   outstream<< oldeopstrings[eop_der.at(i)]<<endl;
       }
     else if (oldeoptypes[eop_comp.at(i)]){
	   outstream<< oldeopstrings[eop_comp.at(i)]<<endl;
	    if (oldeoptypes[eop_der.at(i)])
	      outstream<< oldeopstrings[eop_der.at(i)]<<endl;
     }
   }
      
  outstream<< setw(20) << left <<"NUMBER_OF_ENTRIES"  << num_old+sortepochs.size()<<endl;  
  outstream<< "-HEADER"<< endl;
  
  outstream<< "+DATA"<< endl;
  outstream<< "# All fields are in free format separated by blanks."<<endl;
  outstream <<"# If a parameter was not estimated a filler NA is placed."<<endl;
  outstream<<"#"<<setw(11) << right<< "1" <<" "<< setw(13) << "2" << " "<< setw(13) << "3" << " "<< setw(13) << "4" << " "<< setw(8) << "5" << " "<< setw(8) << "6" << " "<< setw(13) << "7" << " "<< setw(13) << "8" << " "<< setw(13) << "9" << " "<< setw(8) << "10" << " "<< setw(8) << "11" << " "<< setw(8) << "12" << " "<< setw(7) << "13" << " "<< setw(7) << "14" << " "<< setw(7) << "15" << " "<< setw(7) << "16" << " "<< setw(7) << "17" << " "<< setw(10) << "18" << " "<< setw(6) << "19" << " "<< setw(13) << "20" << " "<< setw(13) << "21" << " "<< setw(13) << "22" << " "<< setw(8) << "23" << " "<< setw(8) << "24" << " "<< setw(13) << "25" << " "<< setw(13) << "26" << " "<< setw(13) << "27" << " "<< setw(8) << "28" << " "<< setw(8) << "29" << setw(8)<< "30"<<endl;
  outstream <<"#"<< setw(11) << right<< "epoch" <<" "<< setw(13) << "xPol" << " "<< setw(13) << "yPol" << " "<< setw(13) << "dUT1" << " "<< setw(8) << "dX" << " "<< setw(8) << "dY" << " "<< setw(13) << "sig_xp" << " "<< setw(13) << "sig_yp" << " "<< setw(13) << "sig_UT1" << " "<< setw(8) << "sig_dX" << " "<< setw(8) << "sig_dY" << " "<< setw(8) << "wRMS" << " "<< setw(7) << "cc_xpyp" << " "<< setw(7) << "cc_xput" << " "<< setw(7) << "cc_yput" << " "<< setw(7) << "cc_dXdY" << " "<< setw(7) << "nObs" << " "<< setw(10) << "sessID" << " "<< setw(6) << "span" << " "<< setw(13) << "xPolR" << " "<< setw(13) << "yPolR" << " "<< setw(13) << "LOD" << " "<< setw(8) << "dXR" << " "<< setw(8) << "dYR" << " "<< setw(13) << "sig_xpR" << " "<< setw(13) << "sig_ypR" << " "<< setw(13) << "sig_LOD" << " "<< setw(8) << "sig_dXR" << " "<< setw(8) << "sig_dYR" << setw(8)<< "network"<<endl;
  
   outstream <<"#"<< setw(11) << right<< "[MJD]" <<" "<< setw(13) << "[as]" << " "<< setw(13) << "[as]" << " "<< setw(13) << "[s]" << " "<< setw(8) << "[mas]" << " "<< setw(8) << "[mas]" << " "<< setw(13) << "[as]" << " "<< setw(13) << "[as]" << " "<< setw(13) << "[s]" << " "<< setw(8) << "[mas]" << " "<< setw(8) << "[mas]" << " "<< setw(8) << "[ps]" << " "<< setw(7) << "[-]" << " "<< setw(7) << "[-]" << " "<< setw(7) << "[-]" << " "<< setw(7) << "[-]" << " "<< setw(7) << "[-]" << " "<< setw(10) << "[-]" << " "<< setw(6) << "[h]" << " "<< setw(13) << "[as/d]" << " "<< setw(13) << "[as/d]" << " "<< setw(13) << "[s]" << " "<< setw(8) << "[mas/d]" << " "<< setw(8) << "[mas/d]" << " "<< setw(13) << "[as/d]" << " "<< setw(13) << "[as/d]" << " "<< setw(13) << "[s]" << " "<< setw(8) << "[mas/d]" << " "<< setw(8) << "[mas/d]" << setw(8)<< "[-]"<<endl;
 
  double duration=(sess_en.get_mjd_tt()-sess_st.get_mjd_tt())*24;
  vector<string> allstats=session_ptr->get_trf_ptr()->get_station_names(ivg::staname::lettercode);
  vector<string> stats;
  for (int i=0;i<allstats.size();i++)
    {
      ivg::Analysis_station *sta=session_ptr->get_trf_ptr()->get_station(i);
      if (sta->get_num_obs()>0)
	stats.push_back(allstats.at(i));
    }
  int curoldid=0;
  for (int i=0;i<sortepochs.size();i++) {
    while (curoldid<num_old && olddate.at(curoldid)<=sortepochs.at(i))
      {
	outstream << oldline.at(curoldid) <<endl;
	curoldid++;
      }
     vector<ivg::Date>::iterator  did=std::find(epochs.begin(), epochs.end(), sortepochs.at(i));
    int id=did-epochs.begin();
    outstream <<setfill(' ') << fixed << setw(12) << right << setprecision(6)  << sortepochs.at(i).get_mjd_tt()<<" ";
														      
    for (int j=0;j<2;j++)
      {
	if (colnum(id,j)>=0 )
	  outstream<<setfill(' ') << fixed << setw(13) << right << setprecision(8) << data(id,j)/1000 << " ";
	else
	  outstream<<setfill(' ') << fixed << setw(13) << right<< "NA"<< " ";
      }
      if (colnum(id,2)>=0 )	    
	outstream <<setfill(' ') << fixed << setw(13) << right << setprecision(8)<< -data(id,2)/1000 << " ";
      else
	  outstream<<setfill(' ') << fixed << setw(13) << right<< "NA"<< " ";
	
     for (int j=3;j<5;j++)
      {
	if (colnum(id,j)>=0 )
	  outstream<<setfill(' ') << fixed << setw(8) << right << setprecision(4) << data(id,j) << " ";
	 else
	  outstream<<setfill(' ') << fixed << setw(8) << right<< "NA"<< " ";
      }
    for (int j=0;j<3;j++)
      {
	if (colnum(id,j)>=0 )
	  outstream<<setfill(' ') << fixed << setw(13) << right << setprecision(8) << sigdata(id,j)/1000 << " ";
	else
	  outstream<<setfill(' ') << fixed << setw(13) << right<< "NA"<< " ";
      }
    for (int j=3;j<5;j++)
      {
	if (colnum(id,j)>=0 )
	outstream<<setfill(' ') << fixed << setw(8) << right << setprecision(4) << sigdata(id,j) << " ";
	else
	  outstream<<setfill(' ') << fixed << setw(8) << right<< "NA"<< " ";
      }
    outstream <<setfill(' ') << fixed << setw(8) << right << setprecision(2) << session_ptr->_solution->calc_wrms()*1e12<<" ";
  
 
   
    if (colnum(id,0)>=0 && colnum(id,1)>=0)
      outstream<< setfill(' ') << fixed << setw(7) << right << setprecision(4) <<vcm(colnum(id,0),colnum(id,1))/sqrt(vcm(colnum(id,0),colnum(id,0))*vcm(colnum(id,1),colnum(id,1)))<< " ";
    else
      outstream<< setfill(' ') <<fixed <<  setw(7) << right << "NA"<< " ";
  
    if (colnum(id,0)>=0 && colnum(id,2)>=0)
      
      outstream<< setfill(' ') << fixed << setw(7) << right << setprecision(4)<<vcm(colnum(id,0),colnum(id,2))/sqrt(vcm(colnum(id,0),colnum(id,0))*vcm(colnum(id,2),colnum(id,2)))<< " ";
    else
      outstream<< setfill(' ') << fixed << setw(7) << right << "NA"<< " ";
    if (colnum(id,1)>=0 && colnum(id,2)>=0)
      outstream<< setfill(' ') << fixed << setw(7) << right << setprecision(4)<< vcm(colnum(id,1),colnum(id,2))/sqrt(vcm(colnum(id,1),colnum(id,1))*vcm(colnum(id,2),colnum(id,2)))<< " ";
    else
      outstream<< setfill(' ') << fixed << setw(7) << right << "NA"<< " ";
    if (colnum(id,3)>=0 && colnum(id,4)>=0)
      outstream<< setfill(' ') << fixed << setw(7) << right << setprecision(4) << vcm(colnum(id,3),colnum(id,4))/sqrt(vcm(colnum(id,3),colnum(id,3))*vcm(colnum(id,4),colnum(id,4)))<< " ";
    else
      outstream<< setfill(' ') << fixed << setw(7) << right << "NA"<< " ";
    outstream<< setfill(' ') << fixed << setw(7) << right << setprecision(0) << session_ptr->_solution->get_nobs()<<" ";
    outstream<< setfill(' ') << fixed << setw(10) << right <<sess_name<<" ";
    outstream<< setfill(' ') << fixed << setw(6)<< right << setprecision(2) << duration<<" ";
    for (int j=0;j<2;j++)
      {
	if (colnum(id,j+5)>=0 )
	outstream<<setfill(' ') << fixed << setw(13) << right << setprecision(8)  << data(id,j+5)/1000 << " ";
	else
	  outstream<<setfill(' ') << setw(13) << right<< "NA"<< " ";
      }
    if (colnum(id,7)>=0 )
     outstream <<setfill(' ') << fixed << setw(13) << right << setprecision(8)  << -data(id,7)/1000 << " ";
    else
	  outstream<<setfill(' ') << setw(13) << right<< "NA"<< " ";
	
     for (int j=3;j<5;j++)
      {
	if (colnum(id,j+5)>=0 )
	outstream<<setfill(' ') << fixed << setw(8) << right << setprecision(4) << data(id,j+5) << " ";
	else
	  outstream<<setfill(' ') << setw(8) << right<< "NA"<< " ";
      }
    for (int j=0;j<3;j++)
      {
	if (colnum(id,j+5)>=0 )
	outstream<<setfill(' ') << fixed << setw(13) << right << setprecision(8) << sigdata(id,j+5)/1000 << " ";
	else
	  outstream<<setfill(' ') << setw(13) << right<< "NA"<< " ";
      }
    for (int j=3;j<5;j++)
      {
	if (colnum(id,j+5)>=0 )
	outstream<<setfill(' ') << fixed << setw(8) << right << setprecision(4)  << sigdata(id,j+5) << " ";
	else
	  outstream<<setfill(' ') << setw(8) << right<< "NA"<< " ";
      }
    for (int j=0;j<stats.size();j++){
      
      outstream<<setw(2) << right<<stats.at(j);
      if (j!=stats.size()-1)
	outstream<<"-";
    }
    outstream<<endl;
  }
  while (curoldid<num_old)
    {
      outstream << oldline.at(curoldid) <<endl;
       curoldid++;
    }	       
  outstream<< "-DATA"<< endl;
  outstream<< "%IVS-EOP 3.0 END"<< endl;		  

}
// ...........................................................................
void Session_inout::write_results(ivg::Session *session_ptr, string outfile, bool write_prediction )
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ void Session_inout::write_results(ivg::Session *, string, bool )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    // save results    
    log<INFO>("*** Writing results to ")%outfile;
    ofstream outstream(outfile.c_str(),ios::out);

    if (!outstream.is_open())
        throw runtime_error("void Session_inout::write_results(ivg::Session *, string, bool ): Failed to open file for writing: "+outfile);

    outstream<<session_ptr->_name<<endl;
    outstream<<"---------"<<endl;
    outstream<<" IDX   #  STATION  TYPE O         MJD                  TOTAL                APRIORI               ESTIMATE     STD-DEVIATION   UNIT"<<endl;
    int i = 0;
    for (auto &para : session_ptr->_param_list)
    {
        outstream<<setfill('0')<<setw(4)<<right<<i<<" "<<para.get_resultline(true)<<endl;
        i++;
    }
    if( write_prediction )
    {
        for (auto &para : session_ptr->_param_list.get_stoch_param())
        {
            outstream<<setfill('0')<<setw(4)<<right<<i<<" "<<para.get_resultline(true)<<endl;
            i++;
        }
    }    
    outstream<<endl;


#if DEBUG_VLBI >=2
    cerr<<"--- void Session_inout::write_results(ivg::Session *, string, bool )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Session_inout::read_results(ivg::Session *session_ptr,string infile, bool apriori)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ void Session_inout::read_results(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    // read intern ivg::ASCOT result file and initialize the ivg::Param_list    
    log<INFO>("*** Reading results from ")%infile;
    
    ifstream inStream(infile.c_str(), ios::in);
    if( !inStream.is_open() )
    {
        throw runtime_error( "void read_results(ivg::Session *, Setting *, const string, bool ): Failed to open file: "+infile );
    }
    else
    {
        std::vector<ivg::Param> params, stoch_params;
        
        string line;
        getline(inStream,line,'\n');    //02OCT16XA
        getline(inStream,line,'\n');    //---------
        while (getline(inStream,line,'\n'))
        {
            if(remove_spaces_end(line).length()>0)
            {
                ivg::Param param( line.substr(5,line.length()-4), apriori );
                
                if( param.get_order() != -1 )
                    params.push_back( param );
                else
                    stoch_params.push_back( param );
            }
        }   
        
        session_ptr->_param_list = ivg::Param_list( params, stoch_params );        
    }

#if DEBUG_VLBI >=2
    cerr<<"--- void Session_inout::read_results(ivg::Session *, string, bool )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Session_inout::write_outliers(ivg::Session *session_ptr,string outfile)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ void Session_inout::write_outliers(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    if (_type=="ngs")
    {
        ivg::Matrix W;
        ivg::Matrix wgt_facs;
        session_ptr->_lsa_solution.get_wgt_matrix(W,wgt_facs);
        wgt_facs.save_ascii(outfile);
    }
    else if (_type=="vgosdb")
    {
        string editing = (*session_ptr->_setup)["vgosdb_editing"];
        bool use_wrapper = (bool)(*session_ptr->_setup)["use_wrapper"];
 
        // initialize vgosdb      
        ivg::Vgosdb vgosdb;
        if(use_wrapper){
            if(_wrapper_ptr){
                if(session_ptr->_name.compare(_wrapper_ptr->get_dbName()) == 0 ){
                     vgosdb = ivg::Vgosdb(_session_path,_wrapper_ptr);
                 }
                else{
                    log<WARNING>("!!! session in wrapper and vgosdb are different: ") % session_ptr->_name % " vs. "% _wrapper_ptr->get_dbName();
                    exit(0);
                }
            }
            else{
                log<WARNING>("!!! use wrapper is set in controlfile, but there is no valid pointer to wrapper");
            }
        }
        else
        {
           vgosdb = ivg::Vgosdb(_session_path);
           use_wrapper = false;
        }

        std::string old_edit;
        if( use_wrapper  && _wrapper_ptr->file_exists(ivg::wrapper_entries::Edit,session_ptr->_band_type) )
            old_edit = _wrapper_ptr->get_file(ivg::wrapper_entries::Edit,session_ptr->_band_type);
        else if(editing.empty())
            old_edit = "Edit";        
        else
            old_edit = "Edit_"+editing;
	
        // copy a existing file to a new one
        string edit_save_file = (*session_ptr->_setup)["outliers"]["save"][1];
        if (!vgosdb.does_file_exist("ObsEdit",edit_save_file))
            vgosdb.copy_file("ObsEdit",old_edit,"ObsEdit",edit_save_file);

        // get information about the outliers
        ivg::Matrix W;
        ivg::Matrix wgt_facs;
        session_ptr->_lsa_solution.get_wgt_matrix(W,wgt_facs);

        vector<short> delay_flag;
        string edit_load_file = (*session_ptr->_setup)["outliers"]["load"][1];
        if ((bool)(*session_ptr->_setup)["outliers"]["load"][0])
            delay_flag = vgosdb.get_vector<short>("ObsEdit",edit_load_file,"DelayFlag");
        else
            delay_flag = vgosdb.get_vector<short>("ObsEdit",old_edit,"DelayFlag");

        vector<short> delay_flag_new( delay_flag.size(), 55 );
        int cnt_new_outliers = 0;
        int cnt_bad_delayflag = 0;
        for (int i = 0; i<wgt_facs.rows(); i++)
        {
            if (wgt_facs(i,0)!=1.0&&delay_flag[session_ptr->_origin_obs_idxs.at(i)]==0)
            {
                delay_flag[session_ptr->_origin_obs_idxs.at(i)] = 55;
                cnt_new_outliers++;
            }
            else
            {
                delay_flag_new[session_ptr->_origin_obs_idxs.at(i)] 
                        = delay_flag[session_ptr->_origin_obs_idxs.at(i)];
                
                if( delay_flag_new[session_ptr->_origin_obs_idxs.at(i)] != 0 )
                    cnt_bad_delayflag++;
            }
        }
       
        log<RESULT>("*** #")%cnt_new_outliers%" new outliers and " 
                            %cnt_bad_delayflag%" new bad delay flags set in "
                            %edit_save_file%" with DelayFlag set to 1";
        
        // change a single variable (DelayFlag vector) within the new generated Edit_xxx.nc file
        vgosdb.set_vector<short>(delay_flag_new,"ObsEdit",edit_save_file,"DelayFlag");
        
        // log<WARNING>("*** Inequal length of wgt_facs and delay_flag. No outlier writing possible. MUST set: quality_code = 0 and use_obs_flags = false");
    }


#if DEBUG_VLBI >=2
    cerr<<"--- void Session_inout::write_outliers(ivg::Session *, string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Session_inout::write_groupdelay( ivg::Session *session_ptr )
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ void Session_inout::write_groupdelay(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    // initialize vgosdb based on the whole directory
    ivg::Vgosdb vgosdb(_session_path);
    
//     string editing = (*session_ptr->_setup)["vgosdb_editing"];
//     ivg::Vgosdb vgosdb(_session_path, session_ptr->_name, editing, ivg::band_to_string(session_ptr->_band_type)); 

    // copy a existing file to a new one
    string edit_save_file = (*session_ptr->_setup)["SIM"]["save"][1];
    if (!vgosdb.does_file_exist("ObsEdit",edit_save_file))
        vgosdb.copy_file("ObsEdit","GroupDelayFull_b"+ivg::band_to_string(session_ptr->_band_type),"ObsEdit",edit_save_file);

    // get information about the outliers
    ivg::Matrix W;
    ivg::Matrix wgt_facs;
    session_ptr->_lsa_solution.get_wgt_matrix(W,wgt_facs);

    std::vector<double> sim_group_del( session_ptr->_origin_obs_idxs.back()+1 ,0.0);
//        cerr << "HHHHH" << sim_group_del.size() << ", size of orig idx vec: " << session_ptr->_origin_obs_idxs.size() << endl;

    int c = 0;
    for( int i=0; i <= session_ptr->_scans.size()-1; i++ )
    {
        for( int j=0; j <= session_ptr->_scans.at(i).get_nobs()-1; j++ )
        {
            sim_group_del.at( session_ptr->_origin_obs_idxs.at( c ) ) = session_ptr->_scans.at(i).get_obs_ptr(j)->get_group_delay();
	    if (isnan(sim_group_del.at( session_ptr->_origin_obs_idxs.at( c ) )))
		sim_group_del.at( session_ptr->_origin_obs_idxs.at( c ) )=0.0;
            c++;
        }
    }

//    cerr << " orn vec " << endl;
//    show_vector( session_ptr->_origin_obs_idxs );      
//    cerr << " *******************   " << endl;
//    show_vector( sim_group_del );      

    // change a single variable (DelayFlag vector) within the new generated Edit_xxx.nc file
    vgosdb.set_vector<double>( sim_group_del, "ObsEdit", edit_save_file, "GroupDelayFull" );
        

#if DEBUG_VLBI >=2
    cerr<<"--- void Session_inout::write_groupdelay(ivg::Session *, string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Session_inout::write_obs_met_info( ivg::Session *session_ptr, string outfile )
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ void Session_inout::write_obs_met_info(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    struct obs_met_data{
        std::string session;
        ivg::Date epoch;
        ivg::Matrix met_data;
        double zhd;
    };    
    
    obs_met_data data;
    double mjd, dec_date;
    ivg::Date epo;
    std::string sta1, sta2;
    int sta_idx1, sta_idx2;
    ivg::Analysis_station * sta_ptr1;
    ivg::Analysis_station * sta_ptr2;
    std::map< std::string, std::vector<obs_met_data> > sta_map;
   
    for( int i=0; i <= session_ptr->_scans.size()-1; i++ )
    {
        for( int j=0; j <= session_ptr->_scans.at(i).get_nobs()-1; j++ )
        {
            epo = session_ptr->_scans.at(i).get_obs_ptr(j)->get_epoch();

            session_ptr->_scans.at(i).get_obs_ptr(j)->get_station_names( sta1, sta2 );

            sta_idx1 = session_ptr->_scans.at(i).get_obs_ptr(j)->get_scan_idx(1);
            sta_idx2 = session_ptr->_scans.at(i).get_obs_ptr(j)->get_scan_idx(2);

            sta_ptr1 = session_ptr->_scans.at(i).get_data( sta_idx1 ).sta_ptr;
            sta_ptr2 = session_ptr->_scans.at(i).get_data( sta_idx2 ).sta_ptr;
           
            data = { session_ptr->_name, 
                     epo, 
                     ivg::Matrix( vector<double>{ session_ptr->_scans.at(i).get_data( sta_idx1 ).tropo.get_temperature(),
                                  session_ptr->_scans.at(i).get_data( sta_idx1 ).tropo.get_pressure(),
                                  session_ptr->_scans.at(i).get_data( sta_idx1 ).tropo.get_humidity() }),
                     session_ptr->_scans.at(i).get_data( sta_idx1 ).tropo.get_zhd() };

            sta_map[ sta1 ].push_back( data );
  
            data = { session_ptr->_name, 
                     epo, 
                     ivg::Matrix( vector<double>{ session_ptr->_scans.at(i).get_data( sta_idx2 ).tropo.get_temperature(),
                                  session_ptr->_scans.at(i).get_data( sta_idx2 ).tropo.get_pressure(),
                                  session_ptr->_scans.at(i).get_data( sta_idx2 ).tropo.get_humidity() }),
                     session_ptr->_scans.at(i).get_data( sta_idx2 ).tropo.get_zhd() };

            sta_map[ sta2 ].push_back( data );
        }
    }
    
    ofstream out;
    for( std::map< std::string, std::vector<obs_met_data> >::iterator it = sta_map.begin(); it != sta_map.end(); ++it )
    {
        out.open( (outfile+"obs_info_"+it->first+".txt").c_str(), std::ofstream::out | std::ofstream::app );
        if (!out.is_open())
            throw runtime_error( "void Session_inout::write_obs_info(ivg::Session *, string ): Failed to open file for writing: "+outfile );
        
        for( int k = 0; k < it->second.size(); ++k)
        {
            // write output obs_file
            out << setiosflags(ios::right) << setiosflags(ios::fixed)
               << "$" << left << setw(9) << it->second.at(k).session << " "
               << setw(1) << 1 << " "                        
               << setprecision(6) << setw(11) << it->second.at(k).epoch.get_double_mjd() << " "
               << setprecision(6) << setw(11) << it->second.at(k).epoch.get_decimal_date() << " "
               << left << setw(8) << it->first << " "
               << setprecision(6) << setw(10) << it->second.at(k).met_data(0) << " " 
               << setprecision(6) << setw(10) << it->second.at(k).met_data(1) << " " 
               << setprecision(6) << setw(8) << it->second.at(k).met_data(2) << " " 
               << setprecision(14) << setw(16) << it->second.at(k).zhd << " "  
               << endl;
        }
        out.close();       
    }
  
#if DEBUG_VLBI >=2
    cerr<<"--- void Session_inout::write_obs_met_info(ivg::Session *, string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Session_inout::write_sou(ivg::Session *session_ptr,string outfile)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::write_sou(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    log<INFO>("*** Writing "+outfile+".sou /.lso / .spha");

    //assignment for easy use
    ivg::Trf * trf_ptr = &(session_ptr->_trf);
    ivg::Crf * crf_ptr = &(session_ptr->_crf);
    ivg::Param_list * param_list_ptr = &(session_ptr->_param_list);
    ivg::Ls_solution * sol_ptr = session_ptr->_solution;
    ivg::Date start = session_ptr->_start;
    ivg::Date end = session_ptr->_end;

    string sou_outfile = outfile+".sou";
    string lso_outfile = outfile+".lso";
    string sh_outfile = outfile+".spha";
    ofstream outstream_sou(sou_outfile.c_str(),ios::out);
    ofstream outstream_lso(lso_outfile.c_str(),ios::out);

    map<int,vector<string>> sh_data;

    for (std::vector<ivg::Param>::iterator param = (*param_list_ptr).begin(); param!=(*param_list_ptr).end(); ++param)
    {


        if (param->is_type({ivg::paramtype::ra},
        {
            0
        }))
    {
        if (!(param+1)->is_type({ivg::paramtype::dec},
        {
            0
        }))
        log<WARNING>("*** CRITICAL: declination of source not following right ascension. Not expected.")%param->get_name()%" | "%(param+1)->get_name();


        double ra = param->get_apriori()+param->get_estimate();
        double dec = (param+1)->get_apriori() + (param+1)->get_estimate();

        double std_ra = param->get_standard_deviation();
        double std_dec = (param+1)->get_standard_deviation();

        // at first the source name is the IERS name
        string src_name = param->get_name();
        string src_name_iers = src_name;

        ivg::Source src_tmp(ivg::srcname::iers,src_name,ra,dec);

        int h,m,deg,min;
        double s,sec;
        src_tmp.get_position(h,m,s,deg,min,sec);

        string sign;
        if (deg>0)
            sign = "+";
        else
            sign = "-";

        deg = abs(deg);

        // if possible, take the IVS name
        ivg::Source * src;
        if (crf_ptr->get_source(&src,src_name))
            src_name = src->get_name(ivg::srcname::ivs);

        int n_sessions = src->get_n_sessions();


        if (!src_tmp.is_special_handling())
        {
            outstream_sou<<"SOU_GCO:"
                    <<"  "<<setfill(' ')<<setw(8)<<left<<src_name
                    <<"  R:  "<<setfill('0')<<setw(2)<<right<<h
                    <<"_"<<setfill('0')<<setw(2)<<right<<m
                    <<"_"<<setfill('0')<<setw(11)<<fixed<<setprecision(8)<<right<<s
                    <<" -+ "<<setfill(' ')<<setw(10)<<right<<setprecision(4)<<std_ra*rad2mas
                    <<"  D:  "<<sign<<setfill('0')<<setw(2)<<right<<deg
                    <<"_"<<setfill('0')<<setw(2)<<right<<min
                    <<"_"<<setfill('0')<<setw(10)<<fixed<<setprecision(7)<<right<<sec
                    <<" -+ "<<setfill(' ')<<setw(10)<<right<<setprecision(4)<<std_dec*rad2mas
                    <<"  C:         Obs_used:         Obs_tot:         Ses_used: "<<setw(5)<<right<<n_sessions
                    <<" Ses_tot:       Date_beg: 2002.12.11 Date_end: 2002.12.11"<<endl;
        }
        else
        {
            outstream_lso<<"SOU_LCO:"
                    <<"  "<<setfill(' ')<<setw(8)<<left<<src_name
                    <<"  $00JAN00XA    00   "
                    <<"EPOCH:  "<<fixed<<setprecision(5)<<right<<param->get_epoch().get_decimal_date()
                    <<"  R:  "<<setfill('0')<<setw(2)<<right<<h
                    <<"_"<<setfill('0')<<setw(2)<<right<<m
                    <<"_"<<setfill('0')<<setw(11)<<fixed<<setprecision(8)<<right<<s
                    <<" -+ "<<setfill(' ')<<setw(10)<<right<<setprecision(4)<<std_ra*rad2mas
                    <<"  D:  "<<sign<<setfill('0')<<setw(2)<<right<<deg
                    <<"_"<<setfill('0')<<setw(2)<<right<<min
                    <<"_"<<setfill('0')<<setw(10)<<fixed<<setprecision(7)<<right<<sec
                    <<" -+ "<<setfill(' ')<<setw(10)<<right<<setprecision(4)<<std_dec*rad2mas
                    <<" Corr: -0.5206 Obs:   00 / "<<setw(4)<<right<<n_sessions<<endl;


            vector<string>::const_iterator sh_iter = find(special_handlings.begin(),special_handlings.end(),remove_spaces_end(src_name_iers));
            int pos = sh_iter-special_handlings.begin();

            stringstream part;
            part<<fixed<<setprecision(8)<<right<<param->get_epoch().get_decimal_date()<<setprecision(15)<<" "<<ra<<" "<<std_ra<<" "<<dec<<" "<<std_dec;
            sh_data[pos].push_back(part.str());

        }
    }
    }
    outstream_sou.close();
    outstream_lso.close();


    ofstream outstream_sh(sh_outfile.c_str(),ios::out);
    outstream_sh<<sh_data.size()<<endl;

    for (auto &sh : sh_data)
    {
        outstream_sh<<"# "<<special_handlings.at(sh.first)<<" "<<sh.second.size()<<endl;

        for (auto &line : sh.second)
            outstream_sh<<line<<endl;
    }

    outstream_sh.close();

#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::write_sou(ivg::Session *, string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Session_inout::_read_snx(ivg::Session *session_ptr,Setting *setup,const string path)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::_read_snx(ivg::Session *, Setting *, const string )"<<endl;
    tictoc tim;
    tim.tic();
#endif  
    log<INFO>("*** START-Loading SNX ") % path;
   
    // paramaeter vector for initializing Param_list
    vector<ivg::Param> parameter;
    // we need to scale the NEQ in case of different units (e.g. for EOPs)
    vector<double> scales;

    // data matrices
    ivg::Matrix N(1,1,0.0);
    ivg::Matrix n_side(1,1,0.0);
    vector<double> mjd_parameter;
    
    // to be able to check if N and n corresponde to param_list
    struct dims {
        int N_rows=0;
        int n_rows=0;
    } neq_dims;
    
    // extracting date from database name
    //    ivg::Date date = ivg::Date( path.substr(path.find_last_of("/")+1,9), "YYMMMDD" );

    // intern save structure and relation between apriori and estimate block
    map<string,string> src_assignment;
    map< string, map<ivg::staname,string> > sta_assignment;
    map<int,int> apriori_estimate_assignment;

    // internal storage structure for statistic block
    struct statistics
    {
        string name;
        string shorty;
        double value;
    };
    map<string,statistics> stats;
    
    // need to know if parameter vector initialized or not
    // (because of unknown order of APRIORI and ESTIMATE block)
    bool parameter_vec_init=false;

    // vector containing station names (ivs_names)
    // we need this information for ZZ-files
    vector<string> stations, stations_domes_no;
    vector<ivg::Source> sources;
    
    // if we have discontinuities parameterized in the snx file, we need to store the epochs
    map< string,vector<ivg::Date> > discontinuities;
    
    ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "void _read_snx(ivg::Session *, Setting *, const string ): Failed to open file: "+path );
    }else{

        string line;
        while (getline(inStream,line,'\n'))
        {
//            cerr << line << endl;
            if (line.substr(0,1)=="%"&&line.substr(1,6)!="ENDSNX")
            {

                int y1 = (int) s2d(line.substr(32,2));
                int y2 = (int) s2d(line.substr(45,2));
                int doy1 = (int) s2d(line.substr(35,3));
                int doy2 = (int) s2d(line.substr(48,3));
                int sec1 = (int) s2d(line.substr(39,5));
                int sec2 = (int) s2d(line.substr(52,5));

                (y1<70)?y1 += 2000:y1 += 1900;
                (y2<70)?y2 += 2000:y2 += 1900;

                ivg::Date start(y1,(double) doy1+((double) sec1)/86400.0);
                ivg::Date end(y2,(double) doy2+((double) sec2)/86400.0);

                // setting start and epoch of data from sinex file
                session_ptr->_start = start;
                session_ptr->_end = end;

                statistics NoPH; // Number of Parameters Header
                NoPH.name = "NUMBER OF PARAMETERS HEADER";
                NoPH.shorty = "NoPH";
                NoPH.value = s2d(line.substr(60,5));
                stats[NoPH.shorty] = NoPH;

            }
            else if (line.find("+FILE/COMMENT")!=string::npos)
            {
                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                    session_ptr->_file_comment<<line<<endl;
            }
            else if (line.find("+SOLUTION/STATISTICS")!=string::npos)
            {
                statistics tmp;

                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {

                    if (line.find("NUMBER OF DEGREES OF FREEDOM")!=string::npos)
                    { //NoDoF
                        tmp.name = "NUMBER OF DEGREES OF FREEDOM";
                        tmp.shorty = "NoDoF";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("SQUARE SUM OF O-C")!=string::npos)
                    { //LTPL
                        tmp.name = "SQUARE SUM OF O-C";
                        tmp.shorty = "LTPL";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("NUMBER OF OBSERVATIONS")!=string::npos)
                    { //NoO
                        tmp.name = "NUMBER OF OBSERVATIONS";
                        tmp.shorty = "NoO";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;

                        tmp.name = "NUMBER OF OBSERVATIONS+CONSTS";
                        tmp.shorty = "NoO+NoC";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("WRMS OF POSTFIT RESIDUALS")!=string::npos)
                    { //WoPR
                        tmp.name = "WRMS OF POSTFIT RESIDUALS";
                        tmp.shorty = "WoPR";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("NUMBER OF UNKNOWNS")!=string::npos)
                    { //NoU
                        tmp.name = "NUMBER OF UNKNOWNS";
                        tmp.shorty = "NoU";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("VARIANCE FACTOR")!=string::npos)
                    { //VF
                        tmp.name = "VARIANCE FACTOR";
                        tmp.shorty = "VF";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("SQUARE SUM OF RESIDUALS")!=string::npos)
                    { //VTPV
                        tmp.name = "SQUARE SUM OF RESIDUALS";
                        tmp.shorty = "VTPV";
                        tmp.value = s2d(_DtoE(line.substr(32,line.size()-32)));
                        stats[tmp.shorty] = tmp;
                    }

                }

              }else if(line.find("+SITE/ID")!=string::npos){  
                  
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    
                    if (line.substr(0,1) != "*") {
                        
                        string code = remove_spaces_end(line.substr(1,4));
                        string domes_no = remove_spaces_end(line.substr(9,9));
                        string ivs_name = remove_spaces_end(line.substr(21,8));
                        replace_string_in_place(ivs_name, " ", "_" );
                        
                        // in case of ISHIOKA the domes_no might be missing in sinex file
                        // we need to define it manually! (bad style here)
                        if(ivs_name == "ISHIOKA" )
                            domes_no = "21791S001";

                        map<ivg::staname,string> names;
                        names[ivg::staname::cdp] = code;
                        names[ivg::staname::domes_no] = domes_no;
                        names[ivg::staname::ivs_name] = ivs_name;
                        
                        sta_assignment[code] = names;  
                        sta_assignment[ivs_name] = names; 
                        
                        // avoid duplicated stations in trf
                        if(find(stations_domes_no.begin(), stations_domes_no.end(), domes_no) == stations_domes_no.end())
                            stations_domes_no.push_back(domes_no);
                    }
                }
                
                // SITE/ID block have to be in front of ESTIMATE/APRIORI/EPOCHS block in SINEX file
                session_ptr->_trf =  ivg::Trf( *setup, stations_domes_no, ivg::staname::domes_no, false );
                
            }else if(line.find("+SOLUTION/EPOCHS")!=string::npos){ 
                
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") {
                        
                        string code = remove_spaces_end(line.substr(1,4));
                        int soln =  stoi(remove_spaces_end(line.substr(9,4)));
                        
                        ivg::Date epoch_start( remove_spaces_end(line.substr(16,12)), "SINEX" );
                        ivg::Date epoch_end( remove_spaces_end(line.substr(29,12)), "SINEX" );
                        ivg::Date epoch_mean( remove_spaces_end(line.substr(42,12)), "SINEX" );
                        
                        discontinuities[code].push_back(epoch_start);
                    }
                }
                
            }else if(line.find("+SOURCE/ID")!=string::npos){ 
                
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") {
                        
                        //assignment between code and iers name of source
                        string iers_name = remove_spaces_end(line.substr(6,8));
                        string ivs_name;
                        
                        if( path.find("vie") != string::npos )
                        {
                            ivs_name = iers_name;
                            iers_name = remove_spaces_end(line.substr(15,8));
                        }
                        // in case of gfz-snx files we don't have IAU names
                        else if( path.find("gfz") != string::npos )
                            ivs_name = remove_spaces_end(line.substr(32,8));  
                        // in case of dgf files, the ivs name is located somewhere else
                        else if( path.find("dgf") != string::npos )
                              ivs_name = remove_spaces_end(line.substr(34));  
                        // all other cases
                        else 
                            ivs_name = remove_spaces_end(line.substr(43,8));
                        
                        src_assignment[line.substr(1,4)] = iers_name;

                        ivg::Source tmp_src(ivg::srcname::iers,iers_name);
                        tmp_src.set_name(ivg::srcname::ivs,ivs_name);
                        
                        // in case icrf name available (alsways begins with "J")
                        if(line.substr(15,1) == "J")
                        {
                            string icrf_name = remove_spaces_end(line.substr(15,16));
                            tmp_src.set_name(ivg::srcname::icrf, icrf_name);
                        }

                        sources.push_back(tmp_src);
                    }
                }                
            }
            else if(line.find("+SOLUTION/APRIORI")!=string::npos || 
                    line.find("+SOLUTION/ESTIMATE")!=string::npos)
            {
                // save which block we are parsing right now
                string solution;
                if(line.find("+SOLUTION/APRIORI")!=string::npos)
                    solution = "APRIORI";
                else
                    solution = "ESTIMATE";
                
                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {

                    if (line.substr(0,1)!="*")
                    {

                        // get information from file
                        int index = stoi(line.substr(1,5));
                        string type = remove_spaces_end(line.substr(7,6));
                        string code = line.substr(14,4).c_str();
                        ivg::Date epoch(line.substr(27,12),"SINEX");
                        mjd_parameter.push_back(epoch.get_double_mjd());
                        
                        // unit used in sinex file
                        string unit = remove_spaces_end(line.substr(40,4));
                        
                        double apri_or_esti, const_or_std;
                        
                        stringstream sstr1; 
                        sstr1 << remove_spaces_end(_DtoE(line.substr(47,21))); 
                        sstr1 >> apri_or_esti;
                        
                        stringstream sstr2;
                        sstr2 << remove_spaces_end(_DtoE(line.substr(69,11))); 
                        sstr2 >> const_or_std;
                        
                        // check if parameter vector already has been initialized or not
                        // depends on if the first block is ESTIMATE or APRIORI
                        if(!parameter_vec_init)
                        {
                            vector<string>::const_iterator pos_iter = find( paramtype_snx.begin(), paramtype_snx.end(), type );
                            int pos = pos_iter - paramtype_snx.begin();
                            
                            // we use a struct to save parameter specific information just temporary
                            struct param_init{
                                ivg::paramtype type;
                                string name;
                                int order;
                                double scale;
                            } p_tmp;
                                                       
                            // CRITICAL TO USE HARD CODED NUMBERS HERE !!!!!!!!
                            //   0       1       2       3       4         5         6        7      8      9      10        11       12      13        14       15        16
                           //{ "STAX", "STAY", "STAZ", "CLO", "TROTOT", "TGNTOT", "TGETOT", "XPO", "YPO", "UT1", "NUT_X", "NUT_Y", "RS_RA", "RS_DE", "CL_BR", "NUT_LN", "NUT_OB" } );
                            // stax, stay, staz
                            if(pos >= 0 && pos <=2 )
                                p_tmp = {(ivg::paramtype)pos, sta_assignment[code][ivg::staname::domes_no], 0, 1.0/ivg::param_unit_fac.at(pos) };
                            // trotot, tgntot, tgetot
                            else if(pos >=4 && pos <= 6)
                                p_tmp = {(ivg::paramtype)pos, sta_assignment[code][ivg::staname::domes_no], 0, 1.0/ivg::param_unit_fac.at(pos) };
                            // xpo/ypo
                            else if(pos >=7 && pos <=8)
                                p_tmp = {(ivg::paramtype)pos, "EOP", 0, 1.0/ivg::param_unit_fac.at(pos) };
                            // ut1
                            else if(pos == 9)
                                p_tmp = {(ivg::paramtype)pos, "EOP", 0, 1.0/ivg::param_unit_fac.at(pos)};
                            else if(pos >=10 && pos <=11)
                                p_tmp = {(ivg::paramtype)pos, "EOP", 0, 1.0/ivg::param_unit_fac.at(pos) };
                            else if(pos >= 12 && pos <=13 )
                                p_tmp = {(ivg::paramtype)pos, src_assignment[code], 0, 1.0/ivg::param_unit_fac.at(pos) };
                            else if(pos >= 14)
                            {
                                if(type == "VELX" || type == "VELY" || type == "VELZ" )
                                {
                                    
                                    if(type == "VELX")
                                        p_tmp = {ivg::paramtype::stax, sta_assignment[code][ivg::staname::ivs_name], 1, 1.0/ivg::param_unit_fac.at(0) };
                                    else if(type == "VELY")
                                        p_tmp = {ivg::paramtype::stay, sta_assignment[code][ivg::staname::ivs_name], 1, 1.0/ivg::param_unit_fac.at(1) };
                                    else if(type == "VELZ")
                                        p_tmp = {ivg::paramtype::staz, sta_assignment[code][ivg::staname::ivs_name], 1, 1.0/ivg::param_unit_fac.at(2) }; 
                                   
                                }
                                else if(type == "XPOR")
                                   p_tmp = {ivg::paramtype::xpo, "EOP", 1, 1.0/ivg::param_unit_fac.at(7) };
                                else if(type == "YPOR")
                                   p_tmp = {ivg::paramtype::ypo, "EOP", 1, 1.0/ivg::param_unit_fac.at(8) };
                                else if(type == "UT")
                                   p_tmp = {ivg::paramtype::ut1, "EOP", 0, 1.0/ivg::param_unit_fac.at(9) };
                                else if(type == "LOD")
                                   p_tmp = {ivg::paramtype::ut1, "EOP", 1, -1.0*1.0/ivg::param_unit_fac.at(9) };
                                else if(type == "NUT_LN")
                                     p_tmp = {ivg::paramtype::nutln, "EOP", 0, 1.0/ivg::param_unit_fac.at(15)};
                                else if(type == "NUT_OB")
                                      p_tmp = {ivg::paramtype::nutob, "EOP", 0, 1.0/ivg::param_unit_fac.at(16) };
                                else if(type.find("CLO") != string::npos)
                                    p_tmp = {ivg::paramtype::clo, "----", stoi(type.substr(3)), 1.0/ivg::param_unit_fac.at(3) };
                                else if(type.find("CL_BR") != string::npos)
                                    p_tmp = {ivg::paramtype::cbr, "----", 0, 1.0/ivg::param_unit_fac.at(14) };
                                else
                                    throw runtime_error("void _read_snx(ivg::Session *, Setting *, const string ): Unknown paramtype "+type+" in "+path);
                            }
                            
                            
                            // dgf uses always the same epoch (1.1.2000) for a source in every available snx file
                            // in order to be able to set up sources as local parameter, we need to assign a daily epoch
                            // furthermore we need to scale from mas to rad
                            // UPDATE OCT 2016: THIS SHOULD BE OBSOLETE RIGHT NOW
//                            if(path.find("dgf") != string::npos && (p_tmp.type == ivg::paramtype::ra || p_tmp.type == ivg::paramtype::dec) )
//                            {
//                                double mid_session_mjd = (session_ptr->_start.get_double_mjd() + session_ptr->_end.get_double_mjd() )/2.0;
//                                epoch = ivg::Date(mid_session_mjd);
//                                
//                                p_tmp.scale = 1.0/ivg::rad2mas;
//                            }
                            
                            // scale to internal ascot units
                            scales.push_back(p_tmp.scale);
                            apri_or_esti *= p_tmp.scale;
                            
                            if(solution == "APRIORI")
                                parameter.push_back(ivg::Param(p_tmp.type, p_tmp.name, epoch, apri_or_esti, p_tmp.order));
                            else if(solution == "ESTIMATE")
                            {
                                ivg::Param tmp_param = ivg::Param(p_tmp.type, p_tmp.name, epoch, 0.0, p_tmp.order);
                                tmp_param.set_estimate(apri_or_esti);
                                tmp_param.set_standard_deviation(const_or_std);
                                parameter.push_back(tmp_param);
                            }

                            int vec_position = parameter.size()-1;
                            apriori_estimate_assignment[index] = vec_position;
                            
                        }
                        // if parameter vector already has been initialized
                        else
                        {
                            int vec_position = apriori_estimate_assignment[index];
                            
                            // scale to internal ascot units
                            apri_or_esti *= scales.at(vec_position);
                            
                            if(solution == "APRIORI")
                                parameter.at(vec_position).set_apriori(apri_or_esti);
                            else if(solution == "ESTIMATE")
                            {
                                parameter.at(vec_position).set_estimate(apri_or_esti);
                                parameter.at(vec_position).set_standard_deviation(const_or_std);
                            }
                        }
                    }
                }
                
            statistics tmp; // Number of Parameters Apriori Block
            tmp.name = "NUMBER OF PARAMETERS APRIBLOCK";
            tmp.shorty = "NoPAB";
            tmp.value = parameter.size()/1.0;
            stats[tmp.shorty] = tmp;

            parameter_vec_init = true;
            }
                else if(line.substr(0,line.size()) == "+SOLUTION/CONSTRAINT_INFO"){
                
                double nCons=0;
                
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") 
                       nCons++;
                }

                statistics tmp;
                tmp.name = "NUMBER OF CONSTRAINTS";
                tmp.shorty = "NoC";
                tmp.value = nCons;
                stats[tmp.shorty] = tmp;

                stats["NoO+NoC"].value = nCons+stats["NoO"].value;

            }
            else if (line.find("+SOLUTION/NORMAL_EQUATION_MATRIX")!=string::npos||line.find("+SOLUTION/DECOMPOSED_NORMAL_MATRIX")!=string::npos)
            {

                N.resize(parameter.size(),parameter.size(),0.0);

                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {
                    if (line.substr(0,1)!="*" && line.find("*")==string::npos)
                    {
                        int n = stoi(line.substr(1,5));
                        int m = stoi(line.substr(7,5));
                        int lineSize = line.size();
                        
                        if(parameter.size() > 100)
                            if(n == m && n % 250 == 0)
                                log<INFO>("*** Importing NEQ: Parameter ") % n % " of " % parameter.size();
                                
                        if(n > neq_dims.N_rows)
                             neq_dims.N_rows = n;
                        
                        vector<string> tokens = get_tokens(line.substr(13));
                        int n_tokens = tokens.size();
                        
                        if(n_tokens == 1)
                        {
                            N(n-1,m-1) = s2d(_DtoE(tokens.at(0)));
                            N(m-1,n-1) = N(n-1,m-1);
                        }
                        else if(n_tokens == 2)
                        {
                            N(n-1,m-1) = s2d(_DtoE(tokens.at(0)));
                            N(n-1,m) = s2d(_DtoE(tokens.at(1)));
                            N(m-1,n-1) = N(n-1,m-1);
                            N(m,n-1) = N(n-1,m);
                        }
                        else if(n_tokens == 3)
                        {

                            N(n-1,m-1) = s2d(_DtoE(tokens.at(0)));
                            N(n-1,m) = s2d(_DtoE(tokens.at(1)));
                            N(n-1,m+1) = s2d(_DtoE(tokens.at(2)));

                            N(m-1,n-1) = N(n-1,m-1);
                            N(m,n-1) = N(n-1,m);
                            N(m+1,n-1) = N(n-1,m+1);
                        }
                        else
                        {
                            throw runtime_error("void _read_snx(ivg::Session *, Setting *, const string ): Unexpected line length of matrix block in "+path);
                        }
                    }
                }

            }
            else if (line.find("+SOLUTION/NORMAL_EQUATION_VECTOR")!=string::npos||line.find("+SOLUTION/DECOMPOSED_NORMAL_VECTOR")!=string::npos)
            {
                n_side.resize(parameter.size(),1,0.0);

                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*"){
                        int n = stoi(line.substr(1,5));

                        if(n > neq_dims.n_rows)
                            neq_dims.n_rows = n;
                        
                        // in some cases there is a different type of VECTOR-Block
                        if(line.size() > 50)
                            n_side(n-1,0) = s2d(_DtoE(line.substr(47,23)));
                        else
                            n_side(n-1,0) = s2d(_DtoE(line.substr(7,23)));
                    }
                }
            }
        }
    }
    // close stream after reading the snx-file
    inStream.close();
    
    // now setting all the left member variables from session object "session_ptr"
    // _param_list (_start_epoch, _end_epoch)
    // in case of: _neq_solution (_N, _n, _nparam, _tx, _btPb, _nobs)
    // in case of: _lsa_solution (_Sxx)
    // right now TRF, CRF, EOPs will be based on default aprioris like SSC, ICRF2, C04 aprioris.
    // not until init_snx_solution() the TRF, CRF and EOPs will be adjusted to the param_list
    // TRF based on sites from SITE/ID block -> will contain position and velocities from SSC file
    session_ptr->_crf = ivg::Crf( *setup, sources );
//    session_ptr->_trf =  ivg::Trf( *setup, stations, ivg::staname::ivs_name, false );
    //in case of unusual ivs_name in SITE/ID block (e.g. GFZ sinex files), we need to replace the
    //parameter name right now based on the initialized trf
    for(auto &para: parameter)
    {
        // only in case of a station (x,y,z)
        if(para.is_type({stax,stay,staz},{0,1}))
        {
            string domes_no = para.get_name(); // sta_assignment[para.get_name()][ivg::staname::domes_no];
            ivg::Analysis_station *station;
            session_ptr->_trf.get_station(&station,domes_no);
            para.set_name(station->get_name(ivg::staname::ivs_name));
        }
    }
    
    // now we can initialize the param_list containing correct station-ivs-names
    session_ptr->_param_list =  ivg::Param_list( parameter );
    // in case of several STA/VEL blocks, we have _disconts which will be used in init_snx_solution()
    session_ptr->_disconts = discontinuities;
    //session_ptr->_eops; // already initialized in session constructor 
    // setting start and end epoch of the parameter. information used from first header line   
    session_ptr->_param_list.set_start_end_epoch(session_ptr->_start, session_ptr->_end);
    
    if(neq_dims.N_rows > 0 && (neq_dims.N_rows != parameter.size() || neq_dims.n_rows != parameter.size()))
    {
        stringstream ss;
        ss << "void _read_snx(ivg::Session *, Setting *, const string ): Dimensions of the NEQ does not match with parameter block: N_rows=";
        ss << neq_dims.N_rows << " vs. parameter=" <<  parameter.size() << " in " << path; 
        throw runtime_error( ss.str() );
    }
    
    stringstream stats_warnings;
    // initiialize NEQ object
    if (stats["NoU"].value!=0.0)
    {
        int nparam = (int) stats["NoU"].value;
        session_ptr->_neq_solution = ivg::Ls_neq(N,n_side,nparam,ivg::Matrix(mjd_parameter));
        
        // generate multiplicative inverse
        for(auto &scale: scales)
            scale = 1.0/scale;
                
        // scale neq-solution to internal ascot units
        session_ptr-> _neq_solution.set_scales(scales);
        session_ptr->_neq_solution.scale_system();
        session_ptr->_neq_solution.clear_scales();
    }
    else
        stats_warnings << "|No NUMBER OF UNKNOWNS|";

    if (stats["LTPL"].value!=0.0)
        session_ptr->_neq_solution.set_btPb(stats["LTPL"].value);
    else
        stats_warnings << "|No WEIGHTED SQUARE SUM OF O-C|";

    if (stats["VTPV"].value!=0.0)
        session_ptr->_neq_solution.set_rtPr(stats["VTPV"].value);
    else
        stats_warnings << "|No SQUARE SUM OF RESIDUALS (VTPV)|";

    if (stats["VF"].value!=0.0)
        session_ptr->_neq_solution.set_sigma0_post(stats["VF"].value);
    else
        stats_warnings << "|No VARIANCE FACTOR|";

    if (stats["NoO+NoC"].value!=0.0)
        session_ptr->_neq_solution.set_nobs((int) stats["NoO+NoC"].value);
    else
        stats_warnings << "|No NUMBER OF UNKNOWNS +  SOLUTION/CONSTRAINT_INFO|";
        
    if(!stats_warnings.str().empty())
        log<WARNING>("!!! "+stats_warnings.str()+" in "+path );
    
    
    log<INFO>("*** END-Loading SNX ") % path;
    
#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::_read_snx(ivg::Session *, Setting *, const string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
}
// ...........................................................................
void Session_inout::write_snx(ivg::Session *session_ptr,string outfile, bool incl_results)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::write_snx(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    log<INFO>("*** Writing "+outfile);
 
    string snx_type = (const char*)(*session_ptr->_setup)["export_snx"]["type"];
    Setting &sinex_info = (*session_ptr->_setup)["export_snx"]["sinex_info"];
    //assignment for easy use
    ivg::Trf * trf_ptr = &(session_ptr->_trf);
    ivg::Crf * crf_ptr = &(session_ptr->_crf);
    ivg::Param_list * params = &(session_ptr->_param_list);
    ivg::Ls_solution * sol_ptr = session_ptr->_solution;
    ivg::Icls * icls_ptr = &(session_ptr->_icls_solution);
    ivg::Date start = session_ptr->_start;
    ivg::Date end = session_ptr->_end;

    std::string seperator = "*\n*-------------------------------------------------------------------------------\n*";


    ofstream outstream(outfile.c_str(),ios::out);

    if (!outstream.is_open())
        throw runtime_error("void Session_inout::write_snx(ivg::Session *, string ): Failed to open file for writing: "+outfile);

    stringstream types;
    bool siteID = false;
    bool sourceID = false;

    std::vector<ivg::Param>::iterator iter;

    // -----------------------
    // write sinex header line
    // -----------------------
    if (params->does_include(ivg::paramtype::stax))
    {
        types<<"S ";
        siteID = true;
	
    }
    if (params->does_include(ivg::paramtype::xpo))
        types<<"E ";
    if (params->does_include(ivg::paramtype::zwd))
    {
        types<<"T ";
        siteID = true;
	
    }
    if (params->does_include(ivg::paramtype::ra))
    {
        types<<"C ";
        sourceID = true;
    }

    std::string agency = sinex_info["agency"];
    ivg::Date now;
    now.now();

    //   outstream << "----|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0---|---8| \n" << endl;
    outstream<<"%=SNX 2.02 "<<agency<<" "<<now.get_date_time("YY:DOY:SSSSS")<<" "<<agency
            <<" "<<start.get_date_time("YY:DOY:SSSSS")<<" "<<end.get_date_time("YY:DOY:SSSSS")
            <<" R "<<setfill('0')<<setw(5)<<params->size()<<" "<<2<<" "<<types.str().c_str()<<endl;

    // -----------------------
    // write FILE/REFERENCE block
    // -----------------------
    struct utsname myuts;
    stringstream sys;
    if (!uname(&myuts))
    {
        sys<<myuts.sysname<<" "<<myuts.release<<" "<<myuts.machine
                <<" ("<<myuts.nodename<<")";
    }
    std::string filetype;
    if (session_ptr->_type=="vgosdb")
      filetype="vgosDB file";
    else if (session_ptr->_type=="ngs")
      filetype="NGS file";
    else if (session_ptr->_type=="snx")
      filetype="SINEX file";
    else
      filetype="data";
    outstream<<seperator<<endl;
    outstream<<"+FILE/REFERENCE"<<endl;
    outstream<<"*info_type__________ info_______________________________________________________"<<endl;
    std::string tmpstr=sinex_info["reference"]["description"];
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
              <<" "<<setfill(' ')<<setw(18)<<"DESCRIPTION"
	     << " "<<setfill(' ')<<setw(60)<< tmpstr << endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"OUTPUT"
            <<" "<<setfill(' ')<<setw(60)<<"VLBI solution"<<endl;
    tmpstr=(sinex_info["reference"]["contact"]).c_str();
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"CONTACT"
	     <<" "<<setfill(' ')<<setw(60)<< tmpstr <<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"SOFTWARE"
            <<" "<<setfill(' ')<<setw(60)<<"ASCOT"<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"HARDWARE"
            <<" "<<setfill(' ')<<setw(60)<<sys.str().c_str()<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"INPUT"
            <<" "<<setfill(' ')<< "VLBI session"
	     << " " << session_ptr->_name 
	     << ", "<< filetype<<endl;
    outstream<<"-FILE/REFERENCE"<<endl;


    // ------------------------
    // write FILE/COMMENT block
    // ------------------------   
    outstream<<seperator<<endl;
    outstream<<"+FILE/COMMENT"<<endl;
    
    // write all comments from membervariable into FILE/COMMENT-block
    if(!session_ptr->_file_comment.str().empty())
        outstream<<session_ptr->_file_comment.rdbuf();
    
    if ((snx_type == "COVA")||(incl_results))
    {
        string type =  (*session_ptr->_setup)[ "solver" ];
        if( type == "LSM" || type == "LSC" )
        {   
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(30)<<left<<"WRMS "
		     <<setfill(' ')<<setw(20)<<setprecision(14)<<right<<scientific<<sol_ptr->calc_wrms() <<endl;
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(30)<<left<<"RMS "
                    <<setfill(' ')<<setw(20)<<setprecision(14)<<right<<scientific<<sol_ptr->calc_rms()<<endl;
            
            if( type == "LSC" )
                outstream << setiosflags(ios::left) << setiosflags(ios::fixed)
                          << " " <<setfill(' ')<<setw(30)<<left<< "Least Squares Collocation: "
                          << setfill(' ')<<setw(24)<<right << "# deterministic params " <<scientific<<sol_ptr->get_nparams()
                          << setfill(' ')<<setw(22)<<right << "# stochastic params " <<scientific<<sol_ptr->get_nobs() << endl;
        }
        else if( type == "ICLS" )
        {
            double wrms, rms, vfac;
            icls_ptr->get_statistics( wrms, rms, vfac );
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(30)<<left<<"WRMS "
                    <<setfill(' ')<<setw(20)<<setprecision(14)<<right<<scientific<< wrms <<endl;
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(30)<<left<<"RMS "
                    <<setfill(' ')<<setw(20)<<setprecision(14)<<right<<scientific<< rms<<endl;            
        }

        // only in case of short baselines in WHISP sessions
        if((*session_ptr->_setup).exists("whisp") &&  (bool)(*session_ptr->_setup)["whisp"] )
        {  
            ivg::Matrix dxyz(3,1,0.0);
            dxyz(0) = params->get_param( params->get_index( ivg::paramtype::stax, "WETTZ13N" ) )->get_estimate();
            dxyz(1) = params->get_param( params->get_index( ivg::paramtype::stay, "WETTZ13N" ) )->get_estimate();
            dxyz(2) = params->get_param( params->get_index( ivg::paramtype::staz, "WETTZ13N" ) )->get_estimate();
            
            ivg::Analysis_station * wn;  
            ivg::Analysis_station * wz; 
            double std;

            bool found1 = trf_ptr->get_station(&wz,"WETTZELL");
            bool found2 = trf_ptr->get_station(&wn,"WETTZ13N");

            ivg::Matrix d_xyz = wz->calc_xyz( start ) - wn->calc_xyz( start ) - dxyz;
            ivg::Matrix d_ren = ( wz->form_topo2geo() ).transpose() * d_xyz;      
            double bl = sqrt( pow( d_xyz(0),2 ) + pow( d_xyz(1),2 ) + pow( d_xyz(2),2 ) );
            
            
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(30)<<left<<"WHISP SBAS (XYZ) "
                    <<setfill(' ')<<setw(20)<<setprecision(14)<<right << d_xyz(0) << ", " << d_xyz(1) << ", " << d_xyz(2) << endl;
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(30)<<left<<"WHISP SBAS (REN) "
                    <<setfill(' ')<<setw(20)<<setprecision(14)<<right << d_ren(0) << ", " << d_ren(1) << ", " << d_ren(2) << endl;
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(30)<<left<<"WHISP SBAS length "
                    <<setfill(' ')<<setw(20)<<setprecision(14)<<right << bl << endl;        
        } 
    }    
    outstream<<endl;
    outstream<<"-FILE/COMMENT"<<endl;

    // ------------------------
    // write SITE/ID block
    // ------------------------
    outstream<<seperator<<endl;
    outstream<<"+SITE/ID"<<endl;
    outstream<<"*CODE PT DOMES____ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_"<<endl;

    ivg::Matrix llh(3,1,0.0);

    map<string,string> sta_src_assignment;

    std::vector<Analysis_station>::iterator it_trf;
    for (it_trf = trf_ptr->begin(); it_trf<trf_ptr->end(); it_trf++)
    {
      if (it_trf->get_num_obs()>0) {
        //don't create redundant entries
        if (sta_src_assignment[it_trf->get_name(ivg::staname::ivs_name)].empty())
        {
            sta_src_assignment[it_trf->get_name(ivg::staname::ivs_name)] = it_trf->get_name(ivg::staname::cdp);

            llh = it_trf->calc_lat_lon_h();
            int lon_grad,lon_min,lat_grad,lat_min;
            double lon_sec,lat_sec;
            _dec_lat_2_grad_min_sec(llh(0)*180/M_PI,lat_grad,lat_min,lat_sec);
            _dec_lon_2_grad_min_sec(llh(1)*180/M_PI,lon_grad,lon_min,lon_sec);

            outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                    <<" "<<setfill(' ')<<setw(4)<<it_trf->get_name(ivg::staname::cdp)<<"  A"
                    <<" "<<setfill(' ')<<setw(9)<<it_trf->get_name(ivg::staname::domes_no)<<" R"
                    <<" "<<setfill(' ')<<setw(8)<<left<<it_trf->get_name(ivg::staname::ivs_name)
                    <<" "<<setfill(' ')<<setw(13)<<(it_trf->get_name(ivg::staname::description)).substr(0,13)
                    <<" "<<right<<setfill(' ')<<setw(3)<<fixed<<lon_grad
                    <<" "<<setfill(' ')<<setw(2)<<fixed<<lon_min
                    <<" "<<right<<setw(4)<<fixed<<setprecision(1)<<right<<lon_sec
                    <<" "<<right<<setfill(' ')<<setw(3)<<fixed<<lat_grad
                    <<" "<<setfill(' ')<<setw(2)<<fixed<<lat_min
                    <<" "<<right<<setw(4)<<fixed<<setprecision(1)<<right<<lat_sec
                    <<" "<<setfill(' ')<<setw(7)<<setprecision(1)<<llh(2)<<endl;
        }
      }
    }
    outstream<<"-SITE/ID"<<endl;

    // -------------------------
    // write SITE/RECEIVER block
    // -------------------------
    outstream<<seperator<<endl;
    outstream<<"+SITE/RECEIVER"<<endl;
    outstream<<"*Code PT SBIN T Data_Start__ Data_End____ Receiver Type_______ S/N__ Firmware_ID"<<endl;

    for (it_trf = trf_ptr->begin(); it_trf<trf_ptr->end(); it_trf++)
    {
       if (it_trf->get_num_obs()>0) {
        outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                <<" "<<setfill(' ')<<setw(4)<<it_trf->get_name(ivg::staname::cdp)<<"  A"
                <<" "<<setfill(' ')<<setw(4)<<1<<" R"
                <<" "<<start.get_date_time("YY:DOY:SSSSS")<<" "<<end.get_date_time("YY:DOY:SSSSS")
                <<" "<<setfill(' ')<<setw(20)<<"----VLBI Station----"
                <<" "<<setfill(' ')<<setw(5)<<"--NM-"
                <<" "<<setfill(' ')<<setw(11)<<"-----NA----"<<endl;
       }
    }
    outstream<<"-SITE/RECEIVER"<<endl;

    // -------------------------
    // write SITE/ANTENNA block
    // -------------------------
    outstream<<seperator<<endl;
    outstream<<"+SITE/ANTENNA"<<endl;
    outstream<<"*Code PT SBIN T Data_Start__ Data_End____ Receiver Type_______ S/N__"<<endl;

    for (it_trf = trf_ptr->begin(); it_trf<trf_ptr->end(); it_trf++)
    {
      if (it_trf->get_num_obs()>0) {
        outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                <<" "<<setfill(' ')<<setw(4)<<it_trf->get_name(ivg::staname::cdp)<<"  A"
                <<" "<<setfill(' ')<<setw(4)<<1<<" R"
                <<" "<<start.get_date_time("YY:DOY:SSSSS")<<" "<<end.get_date_time("YY:DOY:SSSSS")
                <<" "<<setfill(' ')<<setw(20)<<"----VLBI Station----"
                <<" "<<setfill(' ')<<setw(5)<<"--NM-"<<endl;
      }
    }
    outstream<<"-SITE/ANTENNA"<<endl;
    // -------------------------
    // write SITE/ECCENTRICITY block
    // -------------------------
    outstream<<seperator<<endl;
    outstream<<"+SITE/ECCENTRICITY"<<endl;
    outstream << "*Code PT SOLN T Data_Start__ Data_End____ typ Apr --> Benchmark (m)_______" <<endl;
    for (it_trf = trf_ptr->begin(); it_trf<trf_ptr->end(); it_trf++)
    {
      if (it_trf->get_num_obs()>0) {
        ivg::Matrix ecc = it_trf->get_eccentricity(start);
        outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                <<" "<<setfill(' ')<<setw(4)<<it_trf->get_name(ivg::staname::cdp)<<"  A"
                <<" "<<setfill(' ')<<setw(4)<<1<<" R"
                <<" "<<start.get_date_time("YY:DOY:SSSSS")<<" "<<end.get_date_time("YY:DOY:SSSSS")
		 <<" "<<"XYZ " << setfill(' ')<<setw(8)<<fixed  <<setprecision(4)<<right<< ecc(0)
	         << " " <<setfill(' ')<<setw(8)<<right<< ecc(1)
		 << " " <<setfill(' ')<<setw(8)<<right<< ecc(2)
		 <<endl;
      }
    }
    outstream<<"-SITE/ECCENTRICITY"<<endl;
    
    // ------------------------
    // write SOURCE/ID block
    // ------------------------
    outstream<<seperator<<endl;
    outstream<<"+SOURCE/ID"<<endl;
    outstream<<"*Code IERS nam ICRF designator  IAU name   IVS name"<<endl;

    std::vector<ivg::Source>::iterator it_crf;
    int counter = 1;

    for (it_crf = crf_ptr->begin(); it_crf<crf_ptr->end(); it_crf++)
    {
        // only write source if use_me == true
        if(it_crf->use_me())
        {
            //don't create redundant entries
            if (sta_src_assignment[it_crf->get_name(ivg::srcname::iers)].empty())
            {
                stringstream ss;
                ss<<setfill('0')<<setw(4)<<setiosflags(ios::right)<<counter;
                sta_src_assignment[it_crf->get_name(ivg::srcname::iers)] = ss.str();

                outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                        <<" "<<setfill('0')<<setw(4)<<setiosflags(ios::right)<<counter
                        <<" "<<setfill(' ')<<setw(8)<<left<<it_crf->get_name(ivg::srcname::iers)
                        <<" "<<setfill(' ')<<setw(16)<<left<<it_crf->get_name(ivg::srcname::icrf)
                        <<" "<<setfill(' ')<<setw(10)<<left<<"----------"
                        <<" "<<setfill(' ')<<setw(8)<<left<<it_crf->get_name(ivg::srcname::ivs)
                        <<endl;

                counter++;
            }
        }
    }
    outstream<<"-SOURCE/ID"<<endl;

    
    // ------------------------------
    // write NUTATION/DATA block
    // ------------------------------
    outstream<<seperator<<endl;
    outstream<<"+NUTATION/DATA"<<endl;
    outstream<< " IAU2000a"<<endl;
    outstream<<"-NUTATION/DATA"<<endl;
    // ------------------------------
    // write PRECESSION/DATA block
    // ------------------------------
    outstream<<seperator<<endl;
    outstream<<"+PRECESSION/DATA"<<endl;
    outstream<< " IAU2006"<<endl;
    outstream<<"-PRECESSION/DATA"<<endl;
    // ------------------------------
    // write SOLUTION/STATISTICS block
    // ------------------------------
    outstream<<seperator<<endl;
    outstream<<"+SOLUTION/STATISTICS"<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(30)<<left<<"NUMBER OF OBSERVATIONS "
	     <<" "<<setfill(' ')<<setw(20)<<right<<sol_ptr->get_nobs()<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(30)<<left<<"NUMBER OF UNKNOWNS "
	     <<" "<<setfill(' ')<<setw(20)<<right<<sol_ptr->get_nparams()<<endl;    
    if ((snx_type == "COVA")||(incl_results))
    {
        outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                <<" "<<setfill(' ')<<setw(30)<<left<<"NUMBER OF DEGREES OF FREEDOM "
		 <<" "<<setfill(' ')<<setw(20)<<right<<sol_ptr->get_degrees_of_freedom()<<endl;
        
        string type =  (*session_ptr->_setup)[ "solver" ];
        if( type == "LSM" || type == "LSC" )
        {   
            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                     <<" "<<setfill(' ')<<setw(30)<<left<<"VARIANCE FACTOR "
                     <<" "<<setfill(' ')<<setw(20)<<setprecision(14)<<right
                     <<scientific<<sol_ptr->calc_posterior_vfac()<<endl;        
        }
        else if( type == "ICLS" )
        {
            double wrms, rms, vfac;
            icls_ptr->get_statistics( wrms, rms, vfac );

            outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
                     <<" "<<setfill(' ')<<setw(30)<<left<<"VARIANCE FACTOR "
                     <<" "<<setfill(' ')<<setw(20)<<setprecision(14)<<right
                     <<vfac<<endl;
        }
    }
    outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
            <<" "<<setfill(' ')<<setw(30)<<left<<"WEIGHTED SQUARE SUM OF O-C "
	     <<" "<<setfill(' ')<<setw(20)<<scientific<<setprecision(14)<<right<<sol_ptr->get_btPb()<<endl;
 
/*    outstream<<setiosflags(ios::left)<<setiosflags(ios::scientific)
            <<" "<<setfill(' ')<<setw(30)<<left<<"SQUARE SUM OF O-C "
            <<setfill(' ')<<setw(20)<<setprecision(14)<<right<<sol_ptr->calc_rms()<<endl;*/
    outstream<<"-SOLUTION/STATISTICS"<<endl;

    // -------------------------
    // write SOLUTION/EPOCHS block
    // -------------------------
    outstream<<seperator<<endl;
    outstream<<"+SOLUTION/EPOCHS"<<endl;
    outstream << "*Code PT SOLN T Data_start__ Data_end____ Mean_epoch__" <<endl;
    for (it_trf = trf_ptr->begin(); it_trf<trf_ptr->end(); it_trf++)
    {
       if (it_trf->get_num_obs()>0) {
      ivg::Date start_st=it_trf->get_first_epoch();
      ivg::Date end_st=it_trf->get_last_epoch();
      ivg::Date mean_ep = ivg::Date((start_st.get_double_mjd()+end_st.get_double_mjd())/2);
        outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                <<" "<<setfill(' ')<<setw(4)<<it_trf->get_name(ivg::staname::cdp)<<"  A"
                <<" "<<setfill(' ')<<setw(4)<<1<<" R"
                <<" "<<start_st.get_date_time("YY:DOY:SSSSS")<<" "<<end_st.get_date_time("YY:DOY:SSSSS")
		 <<" "<<mean_ep.get_date_time("YY:DOY:SSSSS")<<endl;
       }
    }
    outstream<<"-SOLUTION/EPOCHS"<<endl;
    // ------------------------------
    // write SOLUTION/APRIORI or SOLUTION/ESTIMATE block
    // ------------------------------
    // build scales for each param from internal ivg::ASCOT representation to SINEX
    ivg::Matrix scales(params->size(),1,0.0);
    for(std::vector<ivg::Param>::iterator it = params->begin();it!=params->end();++it)
    {
        scales(it-params->begin()) = ivg::param_unit_fac.at(it->get_type());
        if(it->get_type()==ivg::paramtype::ut1&&it->get_order()==1)
            scales(it-params->begin()) *= -1.0;
    }

    int loops = 1;
    if ((snx_type=="COVA")||(incl_results)) // ESTIMATE block only if solution is exported
        loops = 2;
    ivg::Matrix val(params->size(),1,0.0);
    for (int i = 0; i<loops; ++i)
    {
        outstream<<seperator<<endl;
        if (i==0)
        {
            outstream<<"+SOLUTION/APRIORI"<<endl;
            outstream<<"*Index Type__ CODE PT SBIN Ref_epoch___ Unit S Apriori_value________ Constraint_"<<endl;
        }
        else
        {
            outstream<<"+SOLUTION/ESTIMATE"<<endl;
            outstream<<"*Index Type__ CODE PT SBIN Ref_epoch___ Unit S Total_value__________ Formal_erro"<<endl;
        }

        int counter = 1;
        std::string ptype,unit,cdp;
        int degree;

        // loop over all params
        for (std::vector<ivg::Param>::iterator it = params->begin(); it!=params->end(); ++it)
        {
            ptype = paramtype_snx.at((int) it->get_type());
            unit = paramtype_unit.at((int) it->get_type());
            if (i==0)
                val(counter-1) = scales(counter-1)*it->get_apriori();
            if (i==1)
                val(counter-1) += scales(counter-1)*it->get_estimate();

            if (it->get_type()==ivg::paramtype::xpo
                    ||it->get_type()==ivg::paramtype::ypo
                    ||it->get_type()==ivg::paramtype::ut1
                    ||it->get_type()==ivg::paramtype::nutx
                    ||it->get_type()==ivg::paramtype::nuty)
                cdp = "----";
            else
                cdp = sta_src_assignment[it->get_name()];
                
            if (it->get_type()==ivg::paramtype::clo)
            {
                degree = it->get_order();
                std::stringstream tmp;
                tmp<<degree;
                ptype = ptype+tmp.str();
            }
            else if (it->get_type()==ivg::paramtype::xpo&&it->get_order()==1
                    ||it->get_type()==ivg::paramtype::ypo&&it->get_order()==1)
            {
                ptype = ptype+"R";
                unit = "masD";
            }
            else if (it->get_type()==ivg::paramtype::ut1&&it->get_order()==1)
            {
                ptype = "LOD";
            }

            outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                    <<" "<<setfill(' ')<<setw(5)<<counter
                    <<" "<<setfill(' ')<<setw(6)<<left<<ptype
                    <<" "<<setfill(' ')<<setw(4)<<cdp<<"  A"
                    <<" "<<setfill(' ')<<setw(4)<<right<<it->get_stacked()
                    <<" "<<it->get_epoch().get_date_time("YY:DOY:SSSSS")
                    <<" "<<setfill(' ')<<setw(4)<<unit<<" 2"
                    <<setiosflags(ios::right)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)
                    <<val(counter-1)
                    <<setiosflags(ios::right)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(11)<<right<<scientific<<setprecision(5);
            if (i==0)
                outstream<<0.0<<endl;
            if (i==1)
                outstream<<abs(scales(counter-1)*it->get_standard_deviation())<<endl;
            counter++;
        }
        if (i==0)
            outstream<<"-SOLUTION/APRIORI"<<endl;
        if (i==1)
            outstream<<"-SOLUTION/ESTIMATE"<<endl;
    }

    if(snx_type == "COVA")
    {
        // -------------------------------------------
        // write SOLUTION/MATRIX_ESTIMATE L COVA block
        // -------------------------------------------
        ivg::Matrix S = scales.diag();
        ivg::Matrix vcm = S*sol_ptr->get_vcm()*S;

        outstream<<seperator<<endl;
        outstream<<"+SOLUTION/MATRIX_ESTIMATE L COVA"<<endl;
        _write_snx_matrix(outstream,vcm);
        outstream<<"-SOLUTION/MATRIX_ESTIMATE L COVA"<<endl;
    }
    else if(snx_type == "NEQ")
    {
        // -------------------------------------------
        // write  SOLUTION/NORMAL_EQUATION_VECTOR
        // -------------------------------------------
      ivg::Matrix N,n,aplo;
      sol_ptr->get_neq(N,n,aplo,true);
      ivg::Matrix S(N.rows(),N.rows(),0.0);
      for (int i=0;i<N.rows();i++)
	S(i,i)=1/scales(i);
      
      N=S*N*S;
     
      n=S*n;
      aplo=S*aplo;
      outstream<<seperator<<endl;
      outstream<<"+SOLUTION/NORMAL_EQUATION_VECTOR"<<endl;
      //_write_snx_vector(outstream,n);
      outstream<<"*Index Type__ CODE PT SBIN Ref_epoch___ Unit S Total_value"<<endl;
      int counter = 1;
      std::string ptype,unit,cdp;
      int degree;
      
      // loop over all params
      for (std::vector<ivg::Param>::iterator it = params->begin(); it!=params->end(); ++it)
        {
            ptype = paramtype_snx.at((int) it->get_type());
            unit = paramtype_unit.at((int) it->get_type());
            val(counter-1)=n(counter-1);

            if (it->get_type()==ivg::paramtype::xpo
                    ||it->get_type()==ivg::paramtype::ypo
                    ||it->get_type()==ivg::paramtype::ut1
                    ||it->get_type()==ivg::paramtype::nutx
                    ||it->get_type()==ivg::paramtype::nuty)
                cdp = "----";
            else
                cdp = sta_src_assignment[it->get_name()];
                
            if (it->get_type()==ivg::paramtype::clo)
            {
                degree = it->get_order();
                std::stringstream tmp;
                tmp<<degree;
                ptype = ptype+tmp.str();
            }
            else if (it->get_type()==ivg::paramtype::xpo&&it->get_order()==1
                    ||it->get_type()==ivg::paramtype::ypo&&it->get_order()==1)
            {
                ptype = ptype+"R";
                unit = "masD";
            }
            else if (it->get_type()==ivg::paramtype::ut1&&it->get_order()==1)
            {
                ptype = "LOD";
            }

            outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                    <<" "<<setfill(' ')<<setw(5)<<counter
                    <<" "<<setfill(' ')<<setw(6)<<left<<ptype
                    <<" "<<setfill(' ')<<setw(4)<<cdp<<"  A"
                    <<" "<<setfill(' ')<<setw(4)<<right<<it->get_stacked()
                    <<" "<<it->get_epoch().get_date_time("YY:DOY:SSSSS")
                    <<" "<<setfill(' ')<<setw(4)<<unit<<" 2"
                    <<setiosflags(ios::right)<<setiosflags(ios::scientific)
                    <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)
                    <<val(counter-1)
                    <<setiosflags(ios::right)<<setiosflags(ios::scientific)
		     << endl;
            counter++;
        }
        outstream<<"-SOLUTION/NORMAL_EQUATION_VECTOR"<<endl;

        // -------------------------------------------
        // write SOLUTION/NORMAL_EQUATION_MATRIX L block
        // -------------------------------------------
        //ivg::Matrix N(10,10,1.0);
        outstream<<seperator<<endl;
        outstream<<"+SOLUTION/NORMAL_EQUATION_MATRIX L"<<endl;
	outstream << "*Row__ Col__ Norm_Equ_Matrix_Value Norm_Equ_Matrix_Valu2 Norm_Equ_Matrix_Valu3"<<endl;
        _write_snx_matrix(outstream,N);
        outstream<<"-SOLUTION/NORMAL_EQUATION_MATRIX L"<<endl;
	outstream<<seperator<<endl;
	outstream<<"+SOLUTION/NORMAL_CALIBRATION LOADING_EFFECT"<<endl;
	outstream<<"*Index Value________________"<<endl;
	_write_snx_vector(outstream,aplo*-1);
	outstream<<"-SOLUTION/NORMAL_CALIBRATION"<<endl;
	
    }
    // ------------------------------
    // write end line
    // ------------------------------
    outstream<<seperator<<endl;
    outstream<<"%ENDSNX"<<endl;

    //   outstream.close();

#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::write_snx(ivg::Session *, string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Session_inout::write_snx_tro(ivg::Session *session_ptr,string outfile)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::write_snx(ivg::Session *, string )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    log<INFO>("*** Writing "+outfile);
    
    string snx_type = (const char*)(*session_ptr->_setup)["export_snx"]["type"];

    //assignment for easy use
    ivg::Trf * trf_ptr = &(session_ptr->_trf);
    ivg::Crf * crf_ptr = &(session_ptr->_crf);
    ivg::Param_list * params = &(session_ptr->_param_list);
    ivg::Ls_solution * sol_ptr = session_ptr->_solution;
    ivg::Icls * icls_ptr = &(session_ptr->_icls_solution);
    ivg::Date start = session_ptr->_start;
    ivg::Date end = session_ptr->_end;

    std::string seperator = "*\n*-------------------------------------------------------------------------------\n*";


    ofstream outstream(outfile.c_str(),ios::out);

    if (!outstream.is_open())
        throw runtime_error("void Session_inout::write_snx_tro(ivg::Session *, string ): Failed to open file for writing: "+outfile);

    stringstream types;
    bool siteID = false;
    bool sourceID = false;

    std::vector<ivg::Param>::iterator iter;

    // -----------------------
    // write sinex header line
    // -----------------------
    if (params->does_include(ivg::paramtype::stax))
    {
        types<<"S ";
        siteID = true;
    }
    if (params->does_include(ivg::paramtype::xpo))
        types<<"E ";
    if (params->does_include(ivg::paramtype::zwd))
    {
        types<<"T ";
        siteID = true;
    }
    if (params->does_include(ivg::paramtype::ra))
    {
        types<<"C ";
        sourceID = true;
    }

    std::string agency = "IVG";
    ivg::Date now;
    now.now();

    //   outstream << "----|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0---|---8| \n" << endl;
    outstream<<"%=TRO 1.00 "<<agency<<" "<<now.get_date_time("YY:DOY:SSSSS")<<" "<<agency
            <<" "<<start.get_date_time("YY:DOY:SSSSS")<<" "<<end.get_date_time("YY:DOY:SSSSS")
            <<" R "<<setfill('0')<<setw(5)<<params->size()<<" "<<0<<" "<<types.str().c_str()<<endl;

    // -----------------------
    // write FILE/REFERENCE block
    // -----------------------
    struct utsname myuts;
    stringstream sys;
    if (!uname(&myuts))
    {
        sys<<myuts.sysname<<" "<<myuts.release<<" "<<myuts.machine
                <<" ("<<myuts.nodename<<")";
    }

    outstream<<seperator<<endl;
    outstream<<"+FILE/REFERENCE"<<endl;
    outstream<<"*info_type__________ info_______________________________________________________"<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"DESCRIPTION"
            <<" "<<setfill(' ')<<setw(60)<<"IGG Bonn VLBI Group (IVG)"<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"OUTPUT"
            <<" "<<setfill(' ')<<setw(60)<<"VLBI analysis"<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"CONTACT"
            <<" "<<setfill(' ')<<setw(60)<<"bakkari@ivg.uni-bonn.de"<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"SOFTWARE"
            <<" "<<setfill(' ')<<setw(60)<<"ivg::ASCOT"<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"HARDWARE"
            <<" "<<setfill(' ')<<setw(60)<<sys.str().c_str()<<endl;
    outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
            <<" "<<setfill(' ')<<setw(18)<<"INPUT"
            <<" "<<setfill(' ')<<setw(23)<<"ivg::ASCOT VLBI session"
            <<" "<<setfill(' ')<<setw(9)<<endl;
    outstream<<"-FILE/REFERENCE"<<endl;


    // ------------------------
    // write SITE/ID block
    // ------------------------
    outstream<<seperator<<endl;
    outstream<<"+SITE/ID"<<endl;
    outstream<<"*CODE PT DOMES____ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_"<<endl;

    ivg::Matrix llh(3,1,0.0);

    map<string,string> sta_src_assignment;

    std::vector<Analysis_station>::iterator it_trf;
    for (it_trf = trf_ptr->begin(); it_trf<trf_ptr->end(); it_trf++)
    {
        //don't create redundant entries
        if (sta_src_assignment[it_trf->get_name(ivg::staname::ivs_name)].empty())
        {
            sta_src_assignment[it_trf->get_name(ivg::staname::ivs_name)] = it_trf->get_name(ivg::staname::cdp);

            llh = it_trf->calc_lat_lon_h();
            int lon_grad,lon_min,lat_grad,lat_min;
            double lon_sec,lat_sec;
            _dec_lat_2_grad_min_sec(llh(0),lat_grad,lat_min,lat_sec);
            _dec_lon_2_grad_min_sec(llh(1),lon_grad,lon_min,lon_sec);

            outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                    <<" "<<setfill(' ')<<setw(4)<<it_trf->get_name(ivg::staname::cdp)<<"  A"
                    <<" "<<setfill(' ')<<setw(9)<<it_trf->get_name(ivg::staname::domes_no)<<" R"
                    <<" "<<setfill(' ')<<setw(8)<<left<<it_trf->get_name(ivg::staname::ivs_name)
                    <<" "<<setfill(' ')<<setw(13)<<(it_trf->get_name(ivg::staname::description)).substr(0,13)
                    <<" "<<right<<setfill(' ')<<setw(3)<<lon_grad
                    <<" "<<setfill(' ')<<setw(2)<<lon_min
                    <<" "<<right<<setw(3)<<fixed<<setprecision(1)<<right<<lon_sec
                    <<" "<<right<<setfill(' ')<<setw(3)<<lat_grad
                    <<" "<<setfill(' ')<<setw(2)<<lat_min
                    <<" "<<right<<setw(3)<<fixed<<setprecision(1)<<right<<lat_sec
                    <<" "<<setfill(' ')<<setw(7)<<setprecision(1)<<llh(2)<<endl;
        }
    }
    outstream<<"-SITE/ID"<<endl;


    // ------------------------------
    // write TROP/SOLUTION block
    // ------------------------------
    double scale = ivg::param_unit_fac.at( ivg::paramtype::zwd );
    
    outstream << seperator << endl;
    outstream << "+TROP/SOLUTION" << endl;
    outstream << "*SITE ___EPOCH____ TROTOT STDDEV TROWET STDDEV PRESS_ _TEMP _HUMI" << endl;
    int counter = 1;
    std::string site;

    // loop over all params
    for (std::vector<ivg::Param>::iterator it = params->begin(); it!=params->end(); ++it)
    {
        if( it->get_type() == ivg::paramtype::zwd )
        {
            site = sta_src_assignment[ it->get_name() ];

            outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                    << " " << setfill(' ') << setw(4) << site
                    << " " << it->get_epoch().get_date_time("YY:DOY:SSSSS")
                    //<< setiosflags(ios::right) << setiosflags(ios::scientific)
                    << " " << setfill(' ') << setw(6) << right << setprecision(1)
                    << (it->get_apriori()+it->get_estimate())* scale * 1e3
                    << " " << setfill(' ') << setw(6) << right << setprecision(1)
                    << -0.0
                    << " " << setfill(' ') << setw(6) << right << setprecision(1)
                    << it->get_estimate() * scale * 1e3
                    << " " << setfill(' ') << setw(6) << right << setprecision(2)
                    << it->get_standard_deviation() * scale * 1e3
                    << endl;
            
            counter++;       
        }
    }
    
    outstream<<"-TROP/SOLUTION"<<endl;
    
    // ------------------------------
    // write end line
    // ------------------------------
    outstream<<seperator<<endl;
    outstream<<"%ENDSNX"<<endl;

    //   outstream.close();

#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::write_snx_tro(ivg::Session *, string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Session_inout::_write_snx_vector(ofstream &outstream,ivg::Matrix matrix)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::_write_snx_vector( std::string , ivg::Matrix )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    for (int i = 0; i<matrix.rows(); i++)
    {
        outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                 <<" "<<setfill(' ')<<setw(5)<<i+1
                 <<" "<<setfill(' ')<<setw(21)<<right<<scientific
                 <<setprecision(14)<<matrix(i,0)<<endl;
    }
#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::_write_snx_vector( std::string , ivg::Matrix  )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Session_inout::_write_snx_matrix(ofstream &outstream,ivg::Matrix matrix)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Session_inout::_write_snx_matrix( std::string , ivg::Matrix )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    for (int i = 0; i<matrix.rows(); i++)
    {
        for (int j = 0; j<=i; j = j+3)
        {
            if (j+2<=i)
            {
                outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                        <<" "<<setfill(' ')<<setw(5)<<i+1
                        <<" "<<setfill(' ')<<setw(5)<<j+1
                        <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)<<matrix(i,j)
                        <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)<<matrix(i,j+1)
                        <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)<<matrix(i,j+2)
                        <<endl;
            }
            else if (j+1<=i)
            {
                outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                        <<" "<<setfill(' ')<<setw(5)<<i+1
                        <<" "<<setfill(' ')<<setw(5)<<j+1
                        <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)<<matrix(i,j)
                        <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)<<matrix(i,j+1)
                        <<endl;
            }
            else
            {
                outstream<<setiosflags(ios::left)<<setiosflags(ios::fixed)
                        <<" "<<setfill(' ')<<setw(5)<<i+1
                        <<" "<<setfill(' ')<<setw(5)<<j+1
                        <<" "<<setfill(' ')<<setw(21)<<right<<scientific<<setprecision(14)<<matrix(i,j)
                        <<endl;
            }
        }
    }

#if DEBUG_VLBI >=2
    cerr<<"--- Session_inout::_write_snx_matrix( std::string , ivg::Matrix  )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Session_inout::_dec_lat_2_grad_min_sec(double coord,int &latGrad,int &latMin,double &latSec)
// ...........................................................................
{
    latGrad = int(coord);
    double minutes;

    if (coord<0)
        minutes = -(coord-latGrad)*60;
    else
        minutes = (coord-latGrad)*60;

    latMin = int(minutes);
    latSec = (minutes-latMin)*60;

    return;
}
// ...........................................................................
void Session_inout::_dec_lon_2_grad_min_sec(double coord,int & lonGrad,int & lonMin,double & lonSec)
// ...........................................................................
{
    if (coord<0)
        coord = coord+360.0;

    lonGrad = int(coord);
    double minutes = (coord-lonGrad)*60;
    lonMin = int(minutes);
    lonSec = (minutes-lonMin)*60;

    return;
}
// ...........................................................................
void Session_inout::_read_ngs(ivg::Session *session_ptr,Setting *setup,const string path)
// ...........................................................................
{
#if DEBUG_VLBI >=1
    cerr<<"+++ Session_inout::_read_ngs(const string path)"<<endl;
    tictoc tim;
    tim.tic();
#endif
   
    if(!file_exists(path))
        throw runtime_error( "void Session_inout::_read_ngs(ivg::Session *session_ptr, Setting *setup, const string path): Failed to open file "+path+"\n");
   
    // read card 1 of last observation and set class variables for #observations
    // _nobs as well as epoch of last observation _end
    session_ptr->_determine_nobs_ngs(path);
    log<INFO>("*** #nobs: ")%session_ptr->_nobs;

    ifstream inStream(path.c_str());
    if (!inStream.is_open())
        throw runtime_error("void _read_ngs( const string path ): Failed to open file "+path);

    // get first line and extract DB name and Version
    string line;
    getline(inStream,line,'\n');

    // read second line and forget about it
    getline(inStream,line,'\n');

    // read stations
    vector<string> stations;
    bool end = false;
    while (!end)
    {
        getline(inStream,line,'\n');
        if (line.substr(0,4)=="$END")
            end = true;
        else
        {
            stringstream tmp;
            string sta,staT,mount;
            double x,y,z,axOff;

            sta = remove_spaces_end(line.substr(0,8));
            tmp<<line;
            tmp>>staT>>x>>y>>z>>mount >> axOff;

            stations.push_back(sta);
        }
    }

    //extracting date from database name
    ivg::Date date = ivg::Date(path.substr(path.find_last_of("/")+1,9),"YYMMMDD");

    // create TRF
    for (auto &st : stations)
        replace_string_in_place(st," ","_");

    session_ptr->_trf = ivg::Trf(*setup,stations,ivg::staname::ivs_name,true,date.add_days(-10.0),date.add_days(10.0));

    // read sources
    vector<ivg::Source> sources;
    end = false;
    while (!end)
    {
        getline(inStream,line,'\n');
        if (line.substr(0,4)=="$END")
            end = true;
        else
        {
            string src;
            src = remove_spaces_end(line.substr(0,8));

            int raH,raMin,deDeg,deMin;
            double raSec,deSec;

            raH = (int) s2d(line.substr(10,2));
            raMin = (int) s2d(line.substr(13,2));
            raSec = s2d(line.substr(16,12));
            deDeg = (int) s2d(line.substr(30,2));
            deMin = (int) s2d(line.substr(33,2));
            deSec = s2d(line.substr(36,12));
            
            // in case of negative declination, we need to save this information
            // the problem is "- 0" because negative zero doesnt exist
            bool negative = false;
            if(line.substr(29,1) == "-")
                negative = true;

            sources.push_back(ivg::Source(ivg::srcname::ngs,src,raH,raMin,raSec,negative,deDeg,deMin,deSec));
        }
    }

    // create CRF based on the creates source vector
    session_ptr->_crf = ivg::Crf(*setup,sources);

    // skip auxilliary parameters
    end = false;
    while (!end)
    {
        getline(inStream,line,'\n');
        if (line.substr(0,4)=="$END")
            end = true;
    }

    //    // vector of o-c, variances of observations and design matrix (transposed)
    //    ivg::Matrix oc( session_ptr->_nobs,1,0.0 );
    //    ivg::Matrix sigma_squ( session_ptr->_nobs,1,0.0 );
    //    ivg::Matrix obs_mjd( session_ptr->_nobs,1,0.0 );
    //    ivg::Matrix AT;   // wil be re-sized when first observation is initialized

    // data cards
    int obsCounter = 0;
    int ionoCounter = 0;
    int scanCounter = 0;
    string tmpsrc;
    ivg::Date tmpepo;
    bool is_sta1,is_sta2;

    ivg::Obs obs;

    ivg::Obs* curr_obs_ptr;
    vector<double>::iterator at_iter;
    std::string ext_met_data;
    std::string ext_data_type;
    std::string gpt2filename,gpt3filename,grad_type;

    map< string, int > station_indexes;
    vector<string> all_stations = session_ptr->_trf.get_station_names(ivg::staname::ivs_name);
    for(int i=0; i<all_stations.size(); i++)
        station_indexes[all_stations.at(i)]=i;
    
   ivg::Matrix stats(stations.size(),stations.size(),0.0);
    
    while (getline(inStream,line,'\n'))
    {
        int seqNo = (int) s2d(line.substr(70,8));
        int cardNo = (int) s2d(line.substr(78,2));

        switch (cardNo)
        {
        case 1:
        {

            //            if(obsCounter % 1000 == 0){
            //                log<INFO>("*** Delays calculated: ") % obsCounter % "/" % session_ptr->_nobs;
            //            }
            // calculate delay for prior observation to scan (not all cards
            // might be contained in the NGS file, but, card 1
            // is always present)
            if (session_ptr->_scans.size()>0)
            {
                int num_of_params = session_ptr->_trf.get_number_stations()*7+session_ptr->_crf.get_number_sources()*2+5;

                if (session_ptr->_scans.size()==1&&obsCounter==1)
                {
                    session_ptr->_start = session_ptr->_scans.at(0).get_epoch();

                    //Initializing parameter list (Stations, Sources, EOPs) including epoch and aprioris
                    session_ptr->_param_list = Param_list(session_ptr->_trf,session_ptr->_crf, session_ptr->_eops,ivg::Date(0.5*(session_ptr->_start.get_double_mjd()+session_ptr->_end.get_double_mjd())));

                    //                    AT.resize( num_of_params,session_ptr->_nobs,0.0 );
                    //                    at_iter = AT.begin();
                }

                //                log<DETAIL>("*** Observation No.: ") % obsCounter;
                //                curr_obs_ptr->calc_delay( at_iter );
                //                
                //                // save some values in local variables
                //                oc( obsCounter-1 ) = curr_obs_ptr->get_observed_minus_computed();
                //                sigma_squ( obsCounter-1 ) = curr_obs_ptr->get_obs_variance( true );
                //                obs_mjd( obsCounter-1 ) = session_ptr->_scans.back().get_epoch().get_double_mjd();
                //                at_iter += num_of_params;
            }

            curr_obs_ptr = NULL;

            string sta1,sta2,src;

            // Card # 1:
            // Col 1-8: Eight character site name for site 1
            // Col 11-18: Eight character site name for site 2
            // Col 21-28: Eight character source name for radio source
            // Col 30-33: Year of observation (e.g. 1979)
            // Col 35-36: Month
            // Col 38-39: Day
            // Col 41-42: Hour
            // Col 44-45: Minute
            // Col 47-60: Seconds
            // Col 61-70: Run identification code (if desired)
            // Col 71-78: Sequence number
            // Col 79-80: 01

            sta1 = remove_spaces_end(line.substr(0,8));
            sta2 = remove_spaces_end(line.substr(10,8));
            src = remove_spaces_end(line.substr(20,8));

            int y = std::stoi(line.substr(29,4));
            int m = std::stoi(line.substr(34,2));
            int d = std::stoi(line.substr(37,2));
            int h = std::stoi(line.substr(40,2));
            int min = std::stoi(line.substr(43,2));
            double s = s2d(line.substr(47,13));

            ivg::Date epoch(y,m,d,h,min,s);

            // get Analysis_station pointer
            ivg::Analysis_station * station1;
            session_ptr->_trf.get_station(&station1,sta1);

            ivg::Analysis_station * station2;
            session_ptr->_trf.get_station(&station2,sta2);

            ivg::Scan scan;

            bool old_src = false;
            ivg::Date tmpepo(ivg::fake_mjd);
            if (session_ptr->_scans.size()>0)
            {
                old_src = session_ptr->_scans.back().get_source()->check_name(src);
                tmpepo = session_ptr->_scans.back().get_epoch();
            }

            //if( _scans.size() == 0 || ( !old_src && !(epoch == tmpepo) ) )
            if (session_ptr->_scans.size()==0|| !old_src|| !(epoch==tmpepo))
            {
                // this is a new scan
                scan.set_epoch(epoch);
                ivg::Source* source;
                session_ptr->_crf.get_source(&source,src);
                scan.set_source(source);

                ivg::Partials_t2c foo{ivg::Matrix(3,3,0.0)};
                ivg::Partials_t2c * deriv_ptr = &foo;
                ivg::Matrix crf2trf = session_ptr->_eops.form_crf2trf(epoch,true,deriv_ptr);
                scan.set_trf2crf_matrix(crf2trf.transpose(),deriv_ptr);
                session_ptr->_scans.push_back(scan);

            }

            int sta1_idx = session_ptr->_scans.back().add_sta(station1);
            int sta2_idx = session_ptr->_scans.back().add_sta(station2);
            
            // create Obs-Baseline-Station-Statistic
            stats(station_indexes[sta1],station_indexes[sta2]) += 1.0;
            stats(station_indexes[sta2],station_indexes[sta1]) += 1.0;
            stats(station_indexes[sta1],station_indexes[sta1]) += 1.0;
            stats(station_indexes[sta2],station_indexes[sta2]) += 1.0;

            ivg::Obs obs_new(session_ptr,&session_ptr->_scans.back(),sta1_idx,sta2_idx);
            session_ptr->_scans.back().add_obs(obs_new);
            //..............

            obsCounter++;

            // determine pointer to current observation
            int idx = session_ptr->_scans.back().get_nobs()-1;
            curr_obs_ptr = session_ptr->_scans.back().get_obs_ptr(idx);
            break;
        }
        case 2:
        {
            // Card # 2:
            // Col 1-20: Observed delay (ns)
            // Col 21-30: Formal error for the observed delay (ns)
            // Col 31-50: Observed delay rate (ps/sec)
            // Col 51-60: Formal error for the observed delay rate (ps/sec)
            // Col 61-62: Data quality flag (blank or 0 indicates good data)
            // Col 64-65: Delay type (blank if same as in Aux. Par.)
            // Col 67-68: Delay rate type (blank if same as in Aux. Par.)
            // Col 71-78: Sequence number
            // Col 79-80: 02

            double delay,sigmaDelay,delayRate,sigmaRate;
            string qFlag;

            delay = s2d(line.substr(0,19))*1e-9;
            sigmaDelay = s2d(line.substr(20,19))*1e-9;
            delayRate = s2d(line.substr(30,19))*1e-9;
            sigmaRate = s2d(line.substr(50,19))*1e-9;

            qFlag = line.substr(61,1);
            int qcode;
            if (!(bool)(*session_ptr->_setup)["use_obs_flags"])
                curr_obs_ptr->set_use_flag(true);
            else if ((bool)(*session_ptr->_setup)["use_obs_flags"] && (qFlag==" "||qFlag=="0"))
                curr_obs_ptr->set_use_flag(true);
            else
                curr_obs_ptr->set_use_flag(false);

            curr_obs_ptr->set_delay(delay,sigmaDelay,delayRate,sigmaRate);

            break;
        }
        case 3:
            // Card # 3:
            // Col 1-10: Correlation coefficient (0-1)
            // Col 11-20: Formal error for correlation coefficient
            // Col 21-30: Fringe amplitude (J)
            // Col 31-40: Formal error for fringe amplitude (J)
            // Col 41-60: Total fringe phase (radians)
            // Col 61-70: Formal error for total fringe phase (radians)
            // Col 71-78: Sequence number
            // Col 79-80: 03

            double corrCoeff,sigmaCorrCoeff,fringeAmp,sigmaFringeAmp,fringePhase,
                    sigmaFringePhase;

            corrCoeff = s2d(line.substr(0,10));
            sigmaCorrCoeff = s2d(line.substr(10,10));
            fringeAmp = s2d(line.substr(20,10));
            sigmaFringeAmp = s2d(line.substr(30,10));
            fringePhase = s2d(line.substr(40,10));
            sigmaFringePhase = s2d(line.substr(60,10));

            break;
        case 4:
            // Card # 4:
            // Col 1-10: System temperature at site 1 (K)
            // Col 11-15: Formal error for system temperature at site 1 (K)
            // Col 16-25: System temperature at site 2 (K)
            // Col 26-30: Formal error for system temperature at site 2 (K)
            // Col 31-40: Antenna temperature at site 1 (K)
            // Col 41-45: Formal error for antenna temperature at site 1 (K)
            // Col 46-55: Antenna temperature at site 2 (K)
            // Col 56-60: Formal error for antenna temperature at site 2 (K)
            // Col 71-78: Sequence number
            // Col 79-80: 04

            double antT1,sysT1,antT2,sysT2;

            antT1 = s2d(line.substr(30,10));
            sysT1 = s2d(line.substr(0,10));
            antT2 = s2d(line.substr(45,10));
            sysT2 = s2d(line.substr(15,10));

            break;
        case 5:
            // Card # 5:
            // Col 1-10: Cable calibration correction (one-way) for site 1 (ns)
            // Col 11-20: Cable calibration correction (one-way) for site 2 (ns)
            // Col 21-30: Water vapor radiometer parameter at site 1 (ns)
            // Col 31-40: Formal error for water vapor radiometer at site 1 (ns)
            // Col 41-50: Water vapor radiometer parameter at site 2 (ns)
            // Col 51-60: Formal error for water vapor radiometer at site 2 (ns)
            // Col 62: Water vapor radiometer parameter definition code for site 1 as follows:
            // 		0 = parameter is zenith path delay,
            // 		1 = parameter is path delay along line-of-sight
            // Col 64: Water vapor radiometer parameter definition code for site 2 (see above)
            // Col 71-78: Sequence number
            // Col 79-80: 05

            double cable_cal1,cable_cal2;
            cable_cal1 = s2d(line.substr(0,10)) * 1e-9;
            cable_cal2 = s2d(line.substr(10,10)) * 1e-9;

            curr_obs_ptr->set_cable_cal(cable_cal1,cable_cal2);

            break;
        case 6:
            // Card # 6:
            // Col 1-10: Ambient atmospheric temperature at site 1 (deg. C)
            // Col 11-20: Ambient atmospheric temperature at site 2 (deg. C)
            // Col 21-30: Ambient atmospheric barometric pressure at site 1 (mb)
            // Col 31-40: Ambient atmospheric barometric pressure at site 2 (mb)
            // Col 41-50: Ambient atmospheric humidity at site 1
            // Col 51-60: Ambient atmospheric humidity at site 2
            // Col 62: Humidity parameter definition code for site 1 as follows:
            // 		0 = humidity parameter is relative humidity (%)
            // 		1 = humidity parameter is dew point (deg. C)
            // 		2 = humidity parameter is wet bulb temperature (deg. C)
            // Col 64: Humidity parameter definition code for site 2 (see above)
            // Col 71-78: Sequence number
            // Col 79-80: 06
            double temp1,temp2,press1,press2,humidity1,humidity2;
            int humidityC;

            temp1 = s2d(line.substr(0,10));
            temp2 = s2d(line.substr(10,10));
            press1 = s2d(line.substr(20,10));
            press2 = s2d(line.substr(30,10));

            humidity1 = s2d(line.substr(40,10));
            humidity2 = s2d(line.substr(50,10));
            humidityC = (int)s2d(line.substr(61,1));

	    // ext_met_data = string((const char*)(*session_ptr->_setup)["troposphere"]["meteo_data"]);
	    ext_data_type = string((const char*)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][1]);
            ext_met_data = string((const char*)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][2]);
            gpt2filename = string((const char*)(*session_ptr->_setup)["troposphere"]["gpt2_grid_file"]);
	    gpt3filename = string((const char*)(*session_ptr->_setup)["troposphere"]["gpt3_grid_file"]);
	    grad_type = string((const char*)(*session_ptr->_setup)["troposphere"]["gradient"]);
            if( !(bool)(*session_ptr->_setup)["troposphere"]["external_meteo_data"][0] )           
            {
                ext_data_type = "meteo";
                ext_met_data = "insitu";
            }
	    session_ptr->_scans.back().add_scan_meteorology(curr_obs_ptr->get_scan_idx(1),temp1,press1,humidity1,humidityC,ext_data_type, ext_met_data, gpt2filename, gpt3filename,grad_type );
            session_ptr->_scans.back().add_scan_meteorology(curr_obs_ptr->get_scan_idx(2),temp2,press2,humidity2,humidityC,ext_data_type,ext_met_data, gpt2filename, gpt3filename,grad_type);
	    
            //scan.add_scan_meteorology( sta1_idx, temp1, press1, humidity1, 0, ext_data_type, ext_met_data, gpt2filename );
            //scan.add_scan_meteorology( sta2_idx, temp2, press2, humidity2, 0, ext_data_type, ext_met_data, gpt2filename );
            break;
        case 7:
            // Card # 7:
            // Col 1-10: Time difference between the reference epoch (card 1)
            //           and the start of the observation (e.g. -60.) (seconds)
            // Col 11-20: Duration of the observation (seconds)
            // Col 21-30: A priori UTC offset at site 1 (if any) (seconds)
            // Col 31-50: Observation frequency (MHz) (blank if same as in Aux. Par.)
            // Col 51-60: Group delay ambiguity (ns) (blank if same as in Aux. Par.)
            // Col 71-78: Sequence number
            // Col 79-80: 07

            //cerr << "NGS-import: card 7 is defined!!!" << endl;

            break;
        case 8:
            // Card # 8:
            // Col 1-20: Delay ionosphere correction (ns)
            // Col 21-30: Delay ionosphere correction formal error (ns)
            // Col 31-50: Delay rate ionosphere correction (ps/s)
            // Col 51-60: Delay rate ionosphere correction formal error (ps/s)
            // Col 62-63: Ionosphere error flag (0=ionosphere correction OK)
            // Col 71-78: Sequence number
            // Col 79-80: 08

            double ionCorr,ionCorrSigma,ionCorrRate,ionCorrRateSigma;
            int errFlag;

            ionCorr = s2d(line.substr(0,20))*1e-9;
            ionCorrSigma = s2d(line.substr(20,10))*1e-9;
            ionCorrRate = s2d(line.substr(30,20))*1e-9;
            ionCorrRateSigma = s2d(line.substr(50,10))*1e-9;

            errFlag = (int)s2d(line.substr(62,1));

            // sometimes errFlag is good but card_8 only contains zeros!
            // this need to be checked by ionCorr != 0.0
            if (errFlag==0&&((ionCorr!=0.0)||(!(bool)(*session_ptr->_setup)["ionosphere"]["apply"])))
                curr_obs_ptr->set_ion_corr(ionCorr,ionCorrSigma,ionCorrRate,ionCorrRateSigma);
            else
            {
                ionoCounter++;
                curr_obs_ptr->set_ion_corr(666.66,1e-12,0.0,0.0);
                if ((bool)(*session_ptr->_setup)["use_obs_flags"])
                    curr_obs_ptr->set_use_flag(false);
                else
                    curr_obs_ptr->set_use_flag(true);
            }

            break;
        case 9:

            // Card # 9:
            // Col 1-70: Comment text
            // Col 71-78: Sequence number
            // Col 79-80: 09 */

            // Seems to be something different like:
            // delay, sigma delay, delayRate, sigmaDelayRate, however, sigmas are
            // different from card 1!!!???
            break;
        }
    }

    log<RESULT>("*** Observations with invalid ionospheric correction: ")%ionoCounter%" of "%obsCounter;
    log<INFO>("*** #scans: ")%session_ptr->_scans.size();

    if( session_ptr->_nobs != obsCounter )
        throw runtime_error( "void Session_inout::_read_ngs(): Inconsistent number of observations (session_ptr->_nobs != obsCounter)" );

    session_ptr->_end = session_ptr->_scans.back().get_epoch();

    session_ptr->_param_list.set_start_end_epoch(session_ptr->_start,session_ptr->_end);
 
    // check if the session contains clock breaks based on br_info file
    string dbname = path.substr(path.find_last_of("/")+1,9);
    string brdir = (*setup)["brdir"];
    std::string br_file = brdir+"/"+dbname+"_br_info";
    
    // if the fiel exists, create the clockbreak-map which will be used later on in param_list.cpp
    map< ivg::paramtype,multimap<string,double> > breaks;
    if(file_exists(br_file))
    {
        ifstream inStream1;
        string line;

        // if the file exists, take the information and insert it in the multimap
        //map< ivg::paramtype,multimap<string,double> > breaks;
        while( ivg::parser::get_line(br_file, inStream1, line))
        {
            std::string station;
            double jd;

            istringstream info_line(line);
            info_line >> station >> jd;

            breaks[ivg::paramtype::cbr].insert(std::pair<string,double>(station, jd) );
        }

        inStream1.close();
         
        
        log<INFO>("*** Setting ngs clock break information [#") % breaks[ivg::paramtype::cbr].size() % "] based on " % br_file;
    }
    
    inStream.close();
    if((bool)(*(session_ptr->_handling)).exists("handling") && (bool)((*(session_ptr->_handling))["handling"]).exists("add_cbr"))
    {   
        for( int i=0; i<(*(session_ptr->_handling))["handling"]["add_cbr"].getLength(); ++i )
        {
            string station = (*(session_ptr->_handling))["handling"]["add_cbr"][i][0];
            string date_str = (*(session_ptr->_handling))["handling"]["add_cbr"][i][1];

            if(date_str.size() > 18) // in this case, it should be something like 2004/11/16-04:03:00
            {
                int y = std::stoi( date_str.substr( 0,4 ) );
                int m = std::stoi( date_str.substr( 5,2 ) );
                int d = std::stoi( date_str.substr( 8,2 ) );
                int h = std::stoi( date_str.substr( 11,2 ) );
                int min = std::stoi( date_str.substr( 14,2 ) );
                double sec = s2d(date_str.substr( 17,2 ) );

                breaks[ivg::paramtype::cbr].insert(std::pair<string,double>(station.substr(0,8), ivg::Date(y,m,d,h,min,sec).get_jd()) );
            }
            else // in this case, it should be MJD
                breaks[ivg::paramtype::cbr].insert(std::pair<string,double>(station.substr(0,8), ivg::Date(s2d(date_str)).get_jd()) );
            
            log<INFO>("*** Insert additional clock break due to defintion in [handling][add_cbr]");
        }
    }
    session_ptr->_param_list.set_breaks(breaks);
    ivg::Matrix orig_stats = stats;
    ivg::Matrix perc_stats = ivg::Matrix(stations.size(),stations.size(),0.0);
    session_ptr->_obsstats["TRF_REM"] = orig_stats-stats; // stats after flagging/removing
    session_ptr->_obsstats["TRF_ORI"] = orig_stats; // original stats
    session_ptr->_obsstats["TRF_PER"] = perc_stats; // percentage of flagged/removed
    
#if DEBUG_VLBI >=1
    cerr<<"--- void Session::read_ngs( string )"<<" : "<<tim.toc()<<" s "<<endl;
#endif   
}
// ...........................................................................
string Session_inout::_DtoE(string str)
// ...........................................................................
{
    size_t posi = str.find("D");
    if (posi!=string::npos)
    {
        str.replace(posi,1,"e");
    }

    return (str);
}
// ...........................................................................
void Session_inout::_init_skd_session(ivg::Session *session_ptr, Setting *setup, const string dir)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Session_inout::_init_skd_session(ivg::Session *session_ptr, Setting *setup, const string directory)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // start and and of schedule session is based on session name and time from SKED block
    ivg::Date date( session_ptr->_name.substr(0,7), "YYMMMDD" );
   
    // WARNING: Right now this will probably only work in case of sessions within ONE day
    string start_time = (*setup)["SKED"]["start_time"]; // e.g. "07:30:00"
    string end_time = (*setup)["SKED"]["end_time"]; // e.g. "09:30:00"
    
    double start_frac_day = ( s2d( start_time.substr(0,2) )*3600.0 + s2d(start_time.substr(3,2))*60.0 + s2d(start_time.substr(6,2)) )/86400.0;
    double end_frac_day   = ( s2d(   end_time.substr(0,2) )*3600.0 + s2d(  end_time.substr(3,2))*60.0 + s2d(  end_time.substr(6,2)) )/86400.0;
    
    if(end_frac_day < start_frac_day){
        end_frac_day += 1.0;
    }
    
    //sessinfo info = _masterfile.get_session_info(date, 0.0, mastertype::intensive);
    sessinfo info = _masterfile.get_session_info( session_ptr->_name );
    if( info.skdname.empty())
    {
        log<WARNING>("!!! Scheduling a day without planned observations. Using start/end epoch and name_prefix from SKED-block.");
        (*setup)["SKED"]["use_intensive_masterfile_infos"]["time"] = false;
        (*setup)["SKED"]["use_intensive_masterfile_infos"]["name"] = false;
        (*setup)["SKED"]["use_intensive_masterfile_infos"]["stations"] = false;
    }  


    if( (bool)(*setup)["SKED"]["use_intensive_masterfile_infos"]["time"] ){
        ivg::Date start = info.date;
        ivg::Date end = start;
        end.add_secs(info.duration * 3600.0);

        session_ptr->_start = start;
        session_ptr->_end = end;
    } else {
        session_ptr->_start = date.add_days(start_frac_day);
        session_ptr->_end = date.add_days(end_frac_day);
    }
    
    std::cerr << session_ptr->_start.get_date_time("YYYY-MO-DD,HH:MI:SS") << "   " << session_ptr->_end.get_date_time("YYYY-MO-DD,HH:MI:SS") << std::endl;
    
    if( (bool)(*setup)["SKED"]["use_intensive_masterfile_infos"]["name"] ){
        session_ptr->_name = info.skdname;
    } else {
        // define some default-name, e.g. k16013
        string prefix = (*setup)["SKED"]["name_prefix"];
        session_ptr->_name = prefix+session_ptr->_start.get_date_time("YYDOY");
        
    }
    
    
    if( (bool)(*setup)["SKED"]["use_intensive_masterfile_infos"]["stations"] ){
        session_ptr->_trf =  ivg::Trf( (*setup), info.stationnames, ivg::staname::lettercode, true, date.add_days(-10.0), date.add_days(10.0) );
    } else {
        // get stations to be scheduled
        vector<string> stations;
        for( int i=0; i<(*setup)["SKED"]["stations"].getLength(); ++i )
            stations.push_back((*setup)["SKED"]["stations"][i]);
        
        // date for TRF initialization is not that important in case of scheduling
        session_ptr->_trf =  ivg::Trf( (*setup), stations, ivg::staname::ivs_name, true, date.add_days(-10.0), date.add_days(10.0) );
    }
    
    
    bool init_from_sked = (bool)(*setup)["SKED"]["init_from_sked"]["apply"];
    const Setting* ifsSetup = &(*setup)["SKED"]["init_from_sked"];
    
    // sources to get used


    vector<ivg::Source> sources;
    std::string sources_file = dir;
    if(init_from_sked){
        string line;
        ifstream inStream;
        while (ivg::parser::get_line(dir, inStream, line)) {
            if ( boost::algorithm::starts_with(line, "$PARAM")) {
                while (ivg::parser::get_line(dir, inStream, line) && !boost::algorithm::starts_with(line, "$")) {
                    if (boost::algorithm::starts_with(line, "MAXSCAN")) {
                        vector<string> tokens = get_tokens(line);
                        if( (bool)(*ifsSetup)["max_scan"] )
                            (*setup)[ "SKED" ]["max_scan"].operator = ( stod(tokens.at(1))  );
                        if( (bool)(*ifsSetup)["min_scan"] )
                            (*setup)[ "SKED" ]["min_scan"].operator = ( stod(tokens.at(3))  );
                    } else if (boost::algorithm::starts_with(line, "SNR")) {
                        vector<string> tokens = get_tokens(line);
                        if( line.find("MARGIN") == std::string::npos && line.find("-") != std::string::npos){
                            // WARNING only working for sessions with one baseline
                            if( (bool)(*ifsSetup)["snr_min_x"] )
                                (*setup)[ "SKED" ]["snr_min_x"] .operator = ( stod(tokens.at(3))  );
                            if( (bool)(*ifsSetup)["snr_min_s"] )
                                (*setup)[ "SKED" ]["snr_min_s"] .operator = ( stod(tokens.at(6))  );
                        } else if (boost::algorithm::starts_with(line, "SNR MARGIN")) {
                            if( (bool)(*ifsSetup)["snr_margin_x"] )
                                (*setup)[ "SKED" ]["snr_margin_x"] .operator = ( stod(tokens.at(3))  );
                            if( (bool)(*ifsSetup)["snr_margin_s"] )
                                (*setup)[ "SKED" ]["snr_margin_s"] .operator = ( stod(tokens.at(6))  );
                        }                        
                    } else if (  (bool)(*ifsSetup)["min_elevation"] && boost::algorithm::starts_with(line, "ELEVATION")) {
                        vector<string> tokens = get_tokens(line);
                        (*setup)[ "SKED" ]["min_elevation"] .operator = ( stod(tokens.at(2))  );
                    }
                }
            } else if ( boost::algorithm::starts_with(line, "$MAJOR")) {
                   while (ivg::parser::get_line(dir, inStream, line) && !boost::algorithm::starts_with(line, "$")) {
                    if (  (bool)(*ifsSetup)["min_sun_dist"]  && boost::algorithm::starts_with(line, "MinSunDist")) {
                        vector<string> tokens = get_tokens(line);
                        (*setup)[ "SKED" ]["min_sun_dist"].operator = ( stod(tokens.at(1))  );
                    } else if ( (bool)(*ifsSetup)["min_time_src"]  && boost::algorithm::starts_with(line, "MinBetween")) {
                        vector<string> tokens = get_tokens(line);
                        (*setup)[ "SKED" ]["min_time_src"] .operator = ( stod(tokens.at(1))*60  );
                    }
                }
            } else if ( (bool)(*ifsSetup)["sources"] && boost::algorithm::starts_with(line, "$SOURCES")) {
                while (ivg::parser::get_line(dir, inStream, line) && !boost::algorithm::starts_with(line, "$")) {
                    vector<string> tokens = get_tokens(line);
                    if (tokens.size() > 0 && tokens.at(0) != "*" && line.substr(0, 1) != "*") {
                        // in case of a common name, use the common name, e.g. CTA26 instead of 0336-019
                        if (tokens.at(1) == "$")
                            sources.push_back(ivg::Source(ivg::srcname::ivs, tokens.at(0)));
                        else
                            sources.push_back(ivg::Source(ivg::srcname::ivs, tokens.at(1)));
                    }
                }
            }
        }
    } 
    if( init_from_sked == false || (bool)(*ifsSetup)["sources"] == false) {
        sources_file = (const char*)(*setup)["SKED"]["sources"];
        string line;
        ifstream inStream;
        while( ivg::parser::get_line(sources_file, inStream, line) )
        {
            vector<string> tokens = get_tokens( line );
            if(tokens.size() > 0 && tokens.at(0) != "*" && line.substr(0,1) != "*")
            {
                // in case of a common name, use the common name, e.g. CTA26 instead of 0336-019
                if(tokens.at(1) == "$")
                    sources.push_back(ivg::Source(ivg::srcname::ivs, tokens.at(0) ));
                else
                    sources.push_back(ivg::Source(ivg::srcname::ivs, tokens.at(1) ));
            }
        }
    }
    
    std::map<bool,std::string> col { {true, ivg::Logger::get_color("blue")}, {false, ivg::Logger::get_color("white")}};
    std::string white = ivg::Logger::get_color("white");
    
    if( g_verbose >= loglevel::DETAIL){
        std::cout << "*** schedule settings" << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["max_scan"] ]
                            << "max scan duration: " << (double)(*setup)[ "SKED" ]["max_scan"] << " s" << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["min_scan"] ]
                            << "min scan duration: " << (double)(*setup)[ "SKED" ]["min_scan"] << " s" << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["snr_min_x"] ]
                            << "min SNR X        : " << (double)(*setup)[ "SKED" ]["snr_min_x"] << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["snr_min_s"] ]
                            << "min SNR S        : " << (double)(*setup)[ "SKED" ]["snr_min_s"] << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["snr_margin_x"] ]                    
                            << "SNR margin X     : " << (double)(*setup)[ "SKED" ]["snr_margin_x"] << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["snr_margin_s"] ]                    
                            << "SNR margin S     : " << (double)(*setup)[ "SKED" ]["snr_margin_s"] << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["min_sun_dist"] ]                    
                            << "min sun distance : " << (double)(*setup)[ "SKED" ]["min_sun_dist"] << " deg" << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["min_time_src"] ]                    
                            << "min time source  : " << (double)(*setup)[ "SKED" ]["min_time_src"] << " s" << white << std::endl;
        std::cout << "    " << col[ (bool)(*ifsSetup)["min_elevation"] ]                    
                            << "min elevation    : " << (double)(*setup)[ "SKED" ]["min_elevation"] << " deg" << white << std::endl;
    }
    
    session_ptr->_trf.log_data_info_table();
    
    // initializing CRF
    // create crf based on all sources in vector (based on list from e.g. source.cat.geodetic)
    session_ptr->_crf = ivg::Crf( (*setup), sources, &session_ptr->_ephem, sources_file);
    
    // at the beginning there are no observations. They need to be scheduled.
    session_ptr->_nobs = 0;
    session_ptr->_nobs_orig = 0;
    
#if DEBUG_VLBI >=2 
   cerr << "--- Session_inout::_init_skd_session(ivg::Session *session_ptr, Setting *setup, const string directory)" 
        << " : " << tim.toc() << " s " << endl; 
#endif 
}
// ...........................................................................
void Session_inout::_read_skd(ivg::Session *session_ptr, Setting *setup, const string path)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Session_inout::_read_skd(ivg::Session *session_ptr, Setting *setup, const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
    // we need this information for the scan meterology and simulation later on
    string ext_data_type = (*setup)["troposphere"]["external_meteo_data"][1];
    string ext_met_data = (*setup)["troposphere"]["external_meteo_data"][2];
    string gpt2filename = (*setup)["troposphere"]["gpt2_grid_file"];
   
    // we need to initialize following variables in order to be able to simulate
    // _start, _end, _scans containing _observations, _nobs, _nobs_orig, _trf, _crf, _param_list
    
    
    const Setting* ifsSetup = &(*setup)["SKED"]["init_from_sked"];
    
    // in a first iteration we need to open the file and store ifnromation about crf and trf
    // this is because STATIONS-block could be at the end of the skd-file but is needed before the SKED-block
    map<string, ivg::Analysis_station* > sta_assignment;
    ifstream inStream_init(path.c_str(), ios::in);
    if( !inStream_init.is_open() ){
        throw runtime_error( "void Session_inout::_read_skd(ivg::Session *session_ptr, Setting *setup, const string path): Failed to open file: "+path );
    }else{
        log<INFO>("*** START-Loading SKD-file ") % path;
        string line; 
        while (getline(inStream_init,line,'\n'))
        {
            if ( boost::algorithm::starts_with(line, "$PARAM")) {
                while (getline(inStream_init,line,'\n') && !boost::algorithm::starts_with(line, "$")){
                    if (boost::algorithm::starts_with(line, "MAXSCAN")) {
                        vector<string> tokens = get_tokens(line);
                        if( (bool)(*ifsSetup)["max_scan"] )
                            (*setup)[ "SKED" ]["max_scan"].operator = ( stod(tokens.at(1))  );
                        if( (bool)(*ifsSetup)["min_scan"] )
                            (*setup)[ "SKED" ]["min_scan"].operator = ( stod(tokens.at(3))  );
                    } else if (boost::algorithm::starts_with(line, "SNR")) {
                        vector<string> tokens = get_tokens(line);
                        if( line.find("MARGIN") == std::string::npos && line.find("-") != std::string::npos){
                            // WARNING only working for sessions with one baseline
                            if( (bool)(*ifsSetup)["snr_min_x"] )
                                (*setup)[ "SKED" ]["snr_min_x"] .operator = ( stod(tokens.at(3))  );
                            if( (bool)(*ifsSetup)["snr_min_s"] )
                                (*setup)[ "SKED" ]["snr_min_s"] .operator = ( stod(tokens.at(6))  );
                        } else if (boost::algorithm::starts_with(line, "SNR MARGIN")) {
                            if( (bool)(*ifsSetup)["snr_margin_x"] )
                                (*setup)[ "SKED" ]["snr_margin_x"] .operator = ( stod(tokens.at(3))  );
                            if( (bool)(*ifsSetup)["snr_margin_s"] )
                                (*setup)[ "SKED" ]["snr_margin_s"] .operator = ( stod(tokens.at(6))  );
                        }                        
                    } else if (  (bool)(*ifsSetup)["min_elevation"] && boost::algorithm::starts_with(line, "ELEVATION")) {
                        vector<string> tokens = get_tokens(line);
                        (*setup)[ "SKED" ]["min_elevation"] .operator = ( stod(tokens.at(2))  );
                    } else if (line.find("START ")!=string::npos) {
                        // just to have the day toinitialize trf correctly
                        int pos = line.find("START")+6;
                        ivg::Date start(stoi(line.substr(pos,4)),stoi(line.substr(pos+4,3)));
                        start.add_secs(s2d(line.substr(pos+7,2))*3600 + s2d(line.substr(pos+9,2))*60 + s2d(line.substr(pos+11,2)));

                        session_ptr->_start = start;
                    }
                }
            } else if ( boost::algorithm::starts_with(line, "$MAJOR")) {
                while (getline(inStream_init,line,'\n') && !boost::algorithm::starts_with(line, "$")){
                    if (  (bool)(*ifsSetup)["min_sun_dist"]  && boost::algorithm::starts_with(line, "MinSunDist")) {
                        vector<string> tokens = get_tokens(line);
                        (*setup)[ "SKED" ]["min_sun_dist"].operator = ( stod(tokens.at(1))  );
                    } else if ( (bool)(*ifsSetup)["min_time_src"]  && boost::algorithm::starts_with(line, "MinBetween")) {
                        vector<string> tokens = get_tokens(line);
                        (*setup)[ "SKED" ]["min_time_src"] .operator = ( stod(tokens.at(1))*60  );
                    }
                }
            } else if(line.find("$STATIONS")!=string::npos) {
                vector<string> stations, sked_shortcut;
                while (getline(inStream_init,line,'\n') && line.substr(0,1) == "A")
                {
                    vector<string> tokens = get_tokens(line);
                    sked_shortcut.push_back(tokens.at(1));
                    stations.push_back(tokens.at(2));
                }
                // date for TRF initialization is not that important in case of scheduling
                session_ptr->_trf =  ivg::Trf( (*setup), stations, ivg::staname::ivs_name, true, session_ptr->_start.add_days(-10.0), session_ptr->_start.add_days(10.0) );
                session_ptr->_trf.log_data_info_table();

                // store relation in assignment-map, e.g. T = TSUKUB32
                for(int i=0; i<stations.size(); i++)
                {
                    ivg::Analysis_station *sta_ptr;
                    session_ptr->_trf.get_station(&sta_ptr,stations.at(i));
                    sta_assignment[sked_shortcut.at(i)] = sta_ptr;
                }

                if(session_ptr->_trf.get_number_stations() != 2)
                    log<WARNING>("void Session_inout::_read_skd(): Reading of skd-file with more than two stations is experimental.  Individual observation duration for each station is not read.");
            }
            else if(line.find("$SOURCES")!=string::npos)
            {
                vector<ivg::Source> sources;
                streampos oldpos;
                while (getline(inStream_init,line,'\n') && !boost::algorithm::starts_with(line, "$")){
                    vector<string> tokens = get_tokens( line );
                    if(tokens.size() > 0 && tokens.at(0) != "*" && line.substr(0,1) != "*")
                    {
                        // in case of a common name, use the common name, e.g. CTA26 instead of 0336-019
                        if(tokens.at(1) == "$")
                            sources.push_back(ivg::Source(ivg::srcname::ivs, tokens.at(0) ));
                        else
                            sources.push_back(ivg::Source(ivg::srcname::ivs, tokens.at(1) ));
                    }
                    oldpos = inStream_init.tellg();
                }
                inStream_init.seekg (oldpos);
                session_ptr->_crf = ivg::Crf( (*setup), sources, &session_ptr->_ephem, path);
            }
        }
    }
        
    ivg::Schedule schedule(session_ptr);
    
    int nobs = 0;
    ifstream inStream(path.c_str(), ios::in);
    string line;
    while (getline(inStream,line,'\n'))
    {
        if(line.find("$SKED")!=string::npos)
        {
             while (getline(inStream,line,'\n') && line.substr(0,1) != "$"){
                vector<string> tokens = get_tokens(line);
                string epoch = tokens.at(4);
                ivg::Date scan_epoch(session_ptr->_start.get_int_year(), stoi(epoch.substr(2,3)));
                scan_epoch.add_secs(s2d(epoch.substr(5,2))*3600 + s2d(epoch.substr(7,2))*60 + s2d(epoch.substr(9,2)));

                // get current crf2trf and partials matrix for each scan
                ivg::Partials_t2c tmp { ivg::Matrix( 3,3,0.0 ) };
                ivg::Partials_t2c * deriv_ptr = &tmp;
                ivg::Matrix crf2trf = session_ptr->_eops.form_crf2trf( scan_epoch, true, deriv_ptr );

                ivg::Source *src;
                session_ptr->_crf.get_source(&src, tokens.at(0));

                ivg::Scan scan(scan_epoch, src, crf2trf.transpose(), deriv_ptr);
                
                ivg::Matrix k = src->get_unit_vector_ssb();

                string stacodes = tokens.at(9); // e.g. T-VW
                vector<int> sta_indexes;
                for(int i=0; i<stacodes.size(); i+=2){
                    sta_indexes.push_back(scan.add_sta( sta_assignment[stacodes.substr(i,1)] ));
                    ivg::Analysis_station* sta = scan.get_sta_ptr( sta_indexes.back() );
                    scan.set_cable_wrapzone(sta, stacodes.substr(i+1, 1) );
                }           
                int n_sta = sta_indexes.size();

                // WARNING! this is only correct in case of two stations within a scan
                //double scan_duration = s2d(tokens.at(9+n_sta+2));
                // individual observation duration for each station is not read
                
                double scan_duration = s2d(tokens.at(5));
                
                for(int sta_idx_1 = 0; sta_idx_1 < n_sta; sta_idx_1++)
                {
                    for(int sta_idx_2 = sta_idx_1 + 1; sta_idx_2 < n_sta; sta_idx_2++)
                    {
                        // adding meteorological data for simulation
                        scan.add_scan_meteorology(sta_indexes.at(sta_idx_1), -999.000, -999.000, -99900.000, 0, ext_data_type, ext_met_data, gpt2filename );
                        scan.add_scan_meteorology(sta_indexes.at(sta_idx_2), -999.000, -999.000, -99900.000, 0, ext_data_type, ext_met_data, gpt2filename );

                        ivg::Obs obs_new( session_ptr, &scan, sta_indexes.at(sta_idx_1), sta_indexes.at(sta_idx_2));
                        
                        ivg::Analysis_station* sta1 = scan.get_sta_ptr(sta_indexes.at(sta_idx_1));
                        ivg::Analysis_station* sta2 = scan.get_sta_ptr(sta_indexes.at(sta_idx_2));
                        
                        ivg::Matrix azel_sta1 = sta1->calc_az_el(scan_epoch, k, crf2trf ); // [rad]
                        ivg::Matrix azel_sta2 = sta2->calc_az_el(scan_epoch, k, crf2trf ); // [rad]
                        
                        ivg::BandInfo xband, sband;
                        schedule.compute_band_info (*src, *sta1, *sta2, scan_epoch, xband, sband );
             
                        double snr_x, snr_s;
                        double sigma_final = schedule.compute_sigma_snr(scan_duration, *sta1, *sta2, azel_sta1(1), azel_sta2(1), scan_epoch,
                                                                      xband, sband, snr_x, snr_s);
                        
                        obs_new.set_snr( snr_x, snr_s );
                        
                        // set fake delay and final sigma
                        obs_new.set_delay( 0.0, sigma_final, 0.0 );
                        // further features have to be set for calc_delay in case of impact factors
                        obs_new.set_cable_cal(0.0,0.0);
                        
                        
                        obs_new.set_ion_corr( 0.0 , 0.0, 0.0 , 0.0 );

                        //store observation as candidate
                        obs_new.set_use_flag( true );
                        scan.add_obs( obs_new );
                                                
                        nobs++;
                    }
                }
                // store scan within session
                scan.set_schedulded_duration(scan_duration);
                session_ptr->_scans.push_back(scan);
            }
        }
    }
    
    if( session_ptr->_scans.size() == 0 ){
        log<WARNING>("void Session_inout::_read_skd(): No scans in sked file specified.");
        throw runtime_error("void Session_inout::_read_skd():: Failed to read sked file: "+path);
    }
    
    // save information to member variables
    session_ptr->_nobs = nobs;
    session_ptr->_nobs_orig = nobs;
    
    // use scans to determine correct start and end epochs
    // START and END within skd-file cannot be used!!! e.g. see i15197.skd
    session_ptr->_start = session_ptr->_scans.front().get_epoch();
    session_ptr->_end = session_ptr->_scans.back().get_epoch();
    session_ptr->_end.add_secs( session_ptr->_scans.back().get_scheduled_duration() );
    
    session_ptr->find_and_mark_unused_sources();
    
    // initialize _param_list within session for completeness
    session_ptr->init_param_list();
    
#if DEBUG_VLBI >=2 
   cerr << "--- Session_inout::_read_skd(ivg::Session *session_ptr, Setting *setup, const string path)" 
        << " : " << tim.toc() << " s " << endl; 
#endif   
}
// ...........................................................................
void Session_inout::write_skd(ivg::Session *session_ptr, string outfile)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Session_inout::write_skd(...)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
   log<INFO>("*** Writing "+outfile);

    // assignment for easy use
    Setting *setup = session_ptr->_setup;
    ivg::Trf * trf_ptr = &(session_ptr->_trf);
    ivg::Crf * crf_ptr = &(session_ptr->_crf);
    vector<ivg::Scan> scans = session_ptr->_scans;
    ivg::Param_list * param_ptr = &(session_ptr->_param_list);
    ivg::Ls_solution * sol_ptr = session_ptr->_solution;
    ivg::Date start = session_ptr->_start;
    ivg::Date end = session_ptr->_end;

    // initialize outstream for writing
    ofstream outstream(outfile.c_str(),ios::out);

    if (!outstream.is_open())
        throw runtime_error("void Session_inout::write_skd(ivg::Session *session_ptr, string outfile)): Failed to open file for writing: "+outfile);
   
    // ======= $PARAM BLOCK
    outstream << "$EXPER " << session_ptr->get_name() <<  endl;
    outstream << "$PARAM" << endl;
    
    string description = (*setup)["SKED"]["description"];
    outstream << "DESCRIPTION " << description << " " << session_ptr->get_name() << endl;   
    outstream << "SCHEDULING_SOFTWARE ivg::ASCOT" << endl;
    outstream << "SKED_CREATE_DATE ";   
    ivg::Date timestamp;
    timestamp.now();    
    outstream << timestamp.get_date_time("YYYY/MO/DD HH:MI:SS") << endl;
    
    string freq_lc = trf_ptr->get_station(0)->set_channel_setup().freq_lc; // get e.g. "SX" or "8F"
    outstream << "SCHEDULER" << setw(8) << setfill(' ') << right << (const char *)(*setup)["SKED"]["scheduler"] << " ";
    outstream << "CORRELATOR " << setw(3) << setfill(' ') << right << (const char *)(*setup)["SKED"]["correlator"] << " ";
    outstream << "START " << start.get_date_time(("YYYYDOYHHMISS ")) << "END   ";
    outstream << end.get_date_time("YYYYDOYHHMISS") << endl;
    
    outstream << "CALIBRATION " << right << setfill(' ') << setw(5);
    outstream << (double)(*setup)["SKED"]["const_calib"] << " CORSYNCH ";
    outstream << right << setfill(' ') << setw(8) << (double)(*setup)["SKED"]["const_sync"] << " DURATION ";
    outstream << right << setfill(' ') << setw(8) << "196" << endl;
    outstream << "EARLY " << right << setfill(' ') << setw(11) << "0" << " IDLE ";
    outstream << right << setfill(' ') << setw(12) << (double)(*setup)["SKED"]["const_idle"] << " LOOKAHEAD ";
    outstream << right << setfill(' ') << setw(7) << "20" << endl;
    outstream << "MAXSCAN " << right << setfill(' ') << setw(9) << (double)(*setup)["SKED"]["max_scan"];
    outstream << " MINSCAN " << right << setfill(' ') << setw(9) << (double)(*setup)["SKED"]["min_scan"];
    outstream << " MINIMUM " << right << setfill(' ') << setw(9) << "0" << endl; 
    outstream << "MIDTP " << right << setfill(' ') << setw(11) << "0";
    outstream << " MODULAR " << right << setfill(' ') << setw(9) << "1";
    outstream << " MODSCAN " << right << setfill(' ') << setw(9) << "1";
    outstream << " PARITY " << right << setfill(' ') << setw(10) << "100" << endl;
    outstream << "SETUP " << right << setfill(' ') << setw(11) << (double)(*setup)["SKED"]["const_setup"];
    outstream << " SOURCE " << right << setfill(' ') << setw(10) << (double)(*setup)["SKED"]["const_source"];
    outstream << " WIDTH " << right << setfill(' ') << setw(11) << "79" << endl;
    outstream << "CONFIRM " << right << setfill(' ') << setw(9) << "Y";
    outstream << " VSCAN " << right << setfill(' ') << setw(11) << "Y" << endl; 
    outstream << "DEBUG " << right << setfill(' ') << setw(11) << "N";
    outstream << " KEEP_LOG " << right << setfill(' ') << setw(8) << "N";
    outstream << " VERBOSE " << right << setfill(' ') << setw(9) << "N" << endl; 
    outstream << "PRFLAG " << right << setfill(' ') << setw(10) << "YYNN";
    outstream << " SNR " << right << setfill(' ') << setw(13) << "AUTO" << endl; 
    outstream << "FREQUENCY   " << freq_lc << " PREOB      PREOB  MIDOB     MIDOB  POSTOB     POSTOB" << endl;
    outstream << "ELEVATION _" << right << setfill(' ') << setw(6) << setprecision(1);
    outstream << fixed << (double)(*setup)["SKED"]["min_elevation"] << endl;
    outstream << "TAPE_MOTION _" << "  START&STOP" << endl;
    outstream << "TAPE_TYPE " ;
    for(auto &sta: (*trf_ptr))
            outstream << sta.get_name(ivg::staname::lettercode) << " " << sta.get_equip_info().recorder << " ";
    outstream << endl;
    outstream << "TAPE_ALLOCATION _" << "  SCHEDULED" << endl;
    outstream << "SNR ";
    
    vector<string> lcs = trf_ptr->get_station_names(ivg::staname::lettercode);
    for(int i=0; i<lcs.size(); i++)
    {
        for(int j=i+1; j<lcs.size(); j++)
        {
             outstream << lcs.at(i) << "-" << lcs.at(j) << " X" << setw(5) << setfill(' ') << right << setprecision(0) << fixed << (double)(*setup)["SKED"]["snr_min_x"] << " ";
             outstream << lcs.at(i) << "-" << lcs.at(j) << " S" << setw(5) << setfill(' ') << right << setprecision(0) << fixed << (double)(*setup)["SKED"]["snr_min_s"] << " ";
        }
    }
    outstream << endl;
    /*
        for(auto &sta: (*trf_ptr))
    {
            outstream << sta.get_name(ivg::staname::lettercode) << "-" << " X ";              // X abgreifen
            outstream << right << setfill(' ') << setw(13) << (double)(*setup)["SKED"]["snr_min_x"];
            outstream << sta.get_name(ivg::staname::lettercode) << "-" << " S ";              // S abgreifen
            outstream << right << setfill(' ') << setw(13) << (double)(*setup)["SKED"]["snr_min_s"];
    }
    outstream << endl;
   */
    // determine number of different sources
    map<ivg::Source*, bool> tmp_used;
    int num_src=0;
    for(auto &scan: scans)
    {
        if(tmp_used[scan.get_source()] == false)
        {
            num_src++;
            tmp_used[scan.get_source()] = true;
        }
    }

    outstream << "SNR MARGIN X" << setw(4) << setfill(' ') << right << setprecision(0) << fixed << (double)(*setup)["SKED"]["snr_margin_x"]
              << " MARGIN S"     << setw(4) << setfill(' ') << right << setprecision(0) << fixed << (double)(*setup)["SKED"]["snr_margin_s"] << endl;
    
    outstream << "SCAN  ";
    map<ivg::Source*, bool> already_used;
    int cnt = 1;
    for(auto &scan: scans)
    {
        if(already_used[scan.get_source()] == false)
        {
            outstream << right << setfill(' ') << setw(4) << cnt << " 196";
            if( cnt%8 == 0 && cnt != num_src)
                outstream << endl << "SCAN  ";
            cnt++;
            
            already_used[scan.get_source()] = true;
        }
    }
    
    outstream << endl;
   // ======= $OP BLOCK
    outstream << "$OP" << endl;
    for(int i=0; i<2; i++) // WHY TWO TIMES THE SAME?!
    {
        outstream << "XP F YP F DUT F PSI F EPS F" << endl; // WHAT IS THIS?!
        for(auto &sta: (*trf_ptr))
            outstream << sta.get_name(ivg::staname::lettercode) << " AOFF F ARAT F COFF F CRT1 F CRT2 F X F Y F Z F" << endl; // AND THIS?!
        
        for( int j = 0; j < session_ptr->_crf.get_number_sources_inc_unused(); ++j ){
            outstream << right << setfill(' ') << setw(3) << j+1 << " F  ";
                if( (j+1)%8 == 0 && (j+1) != session_ptr->_crf.get_number_sources_inc_unused() )
                    outstream << endl;
        }
        outstream << endl;
    }
    // ======= $MAJOR BLOCK 
    outstream << "$MAJOR" << endl;
    outstream << "Add_ps        " << right << setfill(' ') << setw(7) << fixed << setprecision(1) << (double)(*setup)["SKED"]["add_ps"] << endl;
    outstream << "MinSunDist    " << right << setfill(' ') << setw(7) << fixed << setprecision(1) << (double)(*setup)["SKED"]["min_sun_dist"] << endl;
    outstream << "MinBetween    " << right << setfill(' ') << setw(7) << fixed << setprecision(1) << (double)(*setup)["SKED"]["min_time_src"] / 60.0 << endl; 
    outstream << "MaxAngle      " << right << setfill(' ') << setw(7) << fixed << setprecision(1) << (double)(*setup)["SKED"]["max_slew"] << endl; 
    outstream << "MinAngle      " << right << setfill(' ') << setw(7) << fixed << setprecision(1) << (double)(*setup)["SKED"]["min_slew"] << endl; 
    
    // ======= $MINOR BLOCK
    outstream << "$MINOR" << endl;
    
    // ======= OTHER
    outstream << "$ASTROMETRIC" << endl;
    outstream << "$SRCWT" << endl;
    outstream << "$STATWT" << endl;
    outstream << "$BROADBAND" << endl;
    outstream << "$DOWNTIME" << endl;
    
   // ======= $STATONS BLOCK 
    // A P T H blocks needed
    outstream << "$STATIONS" << endl;
    vector<string> block_strings = {"A","P","T","H"};
    for(auto &block: block_strings)
    {
        for(auto &sta: (*trf_ptr))
        {
            if(block == "A")
            {
                outstream << block << "  ";
                outstream << sta.get_antenna_info().id << " ";
                outstream << setw(8) << setfill(' ') << left << sta.get_name(ivg::staname::ivs_name) << " AZEL   ";
                outstream << sta.get_antenna_info().A_line << endl;
            }
            else if(block == "P")
            {
                outstream << block << "  ";
                outstream << setw(2) << setfill(' ') << left << sta.get_name(ivg::staname::lettercode) << " ";
                outstream << setw(8) << setfill(' ') << left << sta.get_name(ivg::staname::ivs_name) << " ";
                
                ivg::Matrix xyz = sta.calc_xyz(ivg::Date(1997,1.5));
                outstream  << setw(15) << setfill(' ') << right << setprecision(5) << fixed << xyz(0) << " ";
                outstream  << setw(14) << setfill(' ') << right << setprecision(5) << fixed << xyz(1) << " ";
                outstream  << setw(14) << setfill(' ') << right << setprecision(5) << fixed << xyz(2) << "   ";
                outstream << "00000000" << " ";
                
                ivg::Matrix llh = sta.calc_lat_lon_h(ivg::Date(1997,1.5));
                outstream  << setw(7) << setfill(' ') << right << setprecision(2) << fixed << 360.0-llh(1)*(ivg::rad2d) << " ";
                outstream  << setw(7) << setfill(' ') << right << setprecision(2) << fixed << llh(0)*(ivg::rad2d) << " " << (const char *)get_list_element((*setup)["refframes"],(*setup)["trf"])[0] << endl;
                
            }
            else if(block == "T")
            {
                outstream << block << "  ";
                outstream << sta.get_equip_info().line;
                outstream << endl;
            }
            else if(block == "H" && !sta.get_antenna_info().H_line.empty())
            {
                outstream << block << "  ";
                outstream << setw(2) << setfill(' ') << left << sta.get_name(ivg::staname::lettercode) << " ";
                outstream << sta.get_antenna_info().H_line << endl;
            }
        }
    }
    
   // ======= $SKED BLOCK
    outstream << "$SKED" << endl;
    for(auto &scan: scans)
    {
        ivg::Date epoch = scan.get_epoch();
        outstream << left << setfill(' ') << setw(8) << scan.get_source()->get_name(ivg::srcname::ivs) << "  10 " << freq_lc << " PREOB  ";
        outstream << epoch.get_date_time("YY");
        outstream << right << setfill('0') << setw(3) << epoch.get_int_doy();
        outstream << epoch.get_date_time("HHMISS") << " ";
        outstream << right << setfill(' ') << setw(7) << setprecision(0) << fixed << scan.get_scheduled_duration() << " MIDOB       0 POSTOB ";
        
        string dummy;
        stringstream durations;
        for(auto *station: scan.get_stations())
        {
            outstream << station->get_antenna_info().id << station->determine_wrap_zone(scan.get_cable_wrap(station));
            dummy += " 1F000000";
            // IMPORTANT: Here we need (different) durations for each station.
            durations << " " << right << setfill(' ') << setw(5) << scan.get_scheduled_duration(); 
        }
        
        outstream << dummy << " YYNN" << durations.str();
        outstream << endl;
    }
        
   // ======= $SOURCES BLOCK
    outstream << "$SOURCES" << endl;
    for( ivg::Source& sou : session_ptr->_crf.get_sources()  ){
         outstream << " " << left << setfill(' ') << setw(8) << sou.get_name(ivg::srcname::iers) << " ";
            if(sou.get_name(ivg::srcname::ivs)  == sou.get_name(ivg::srcname::iers))
                outstream << left << setfill(' ') << setw(8) << "$" << " ";
            else
                outstream << left << setfill(' ') << setw(8) << sou.get_name(ivg::srcname::ivs) << " ";
            
            int h,m,deg,min;
            double s,sec;
            sou.get_position(h,m,s,deg,min,sec);
            outstream << right << setfill('0') << setw(2) << h << " ";
            outstream << right << setfill('0') << setw(2) << m << " ";
            outstream << right << setfill('0') << setw(9) << fixed << setprecision(6) << s << "     ";
            if(deg>=0)
                outstream << "+";
            else
                outstream << "-";
            
            outstream << right << setfill('0') << setw(2) << abs(deg) << " ";
            outstream << right << setfill('0') << setw(2) << min << " ";
            outstream << right << setfill('0') << setw(9) << fixed << setprecision(6) << sec << " ";
            
            outstream << "2000.0 0.0  ICRF2" << endl; // TO DO: WHAT SHOULD WE DO HERE?
    }

    
   // ======= $FLUX BLOCK
    outstream << "$FLUX" << endl;
    for( ivg::Source& sou : session_ptr->_crf.get_sources()  ){
        outstream << sou.get_band_flux_info(ivg::band::X).flux_line << endl;
        outstream << sou.get_band_flux_info(ivg::band::S).flux_line << endl;
    }
    
   // ======= $CODES BLOCK
    outstream << "$CODES" << endl;
    // cha_bbc.first or cha_bbc.second ?!?!?!
    for(auto &sta: (*trf_ptr))
    {
        outstream << "F " << setw(8) << setfill(' ') << left << sta.set_channel_setup().freq_sq_name << " ";
        outstream << setw(8) << setfill(' ') << left << sta.set_channel_setup().freq_lc << " ";
        outstream << sta.get_name(ivg::staname::ivs_name) << endl;
        
        for(auto &cha_bbc: sta.set_channel_setup().cha_bbc)
        {
            ivg::Frequency tmp = sta.set_channel_setup().freq_info[cha_bbc.first];
            outstream << "C " << sta.set_channel_setup().freq_lc << " " << tmp.band << "  " << setprecision(2) << fixed << tmp.frequency << "  " << setprecision(1) << fixed << tmp.phase_cal_freq << " ";
            outstream << right << setw(4) << setfill(' ')  << cha_bbc.first << " ";
            outstream << sta.set_channel_setup().rec_format << tmp.fanout_fac << "    8.00 " << tmp.track_info; 
            outstream << endl;
        }
    } 
    outstream << "R " << freq_lc << "   " << "16.000" << endl;
    outstream << "B " << freq_lc << endl;
    // cha_bbc.first or cha_bbc.second ?!?!?!
    for(auto &sta: (*trf_ptr))
    {   
        for(auto &cha_bbc: sta.set_channel_setup().cha_bbc)
        {
            ivg::Frequency tmp = sta.set_channel_setup().freq_info[cha_bbc.first];
            outstream << "L " << sta.get_antenna_info().id << " " << sta.set_channel_setup().freq_lc << " ";
            outstream << tmp.band << " ";
            outstream << right << setfill(' ') << setw(3) << tmp.if_channel << " ";
            outstream << right << setfill(' ') << setw(10) << fixed << tmp.lo_frequency << " ";
            outstream << right << setfill(' ') << setw(2) << cha_bbc.first << " " << tmp.sideband; 
            outstream << endl;
        }
    } 
    
   // ======= $HEAD BLOCK
    outstream << "$HEAD" << endl;
    for(auto &sta: (*trf_ptr))
        outstream << sta.set_channel_setup().hdpos_lines << endl; // originally from hdpos.cat
    
#if DEBUG_VLBI >=2 
   cerr << "--- Session_inout::write_skd(...)" 
        << " : " << tim.toc() << " s " << endl; 
#endif   
}

void Session_inout::setWrapper_ptr(ivg::Wrapper* _wrapper_ptr) {
    this->_wrapper_ptr = _wrapper_ptr;
}
// ...........................................................................
void Session_inout::_read_utas_apriori(ivg::Session *session_ptr,Setting *setup,const string path)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Session_inout::_read_utas_apriori(ivg::Session *session_ptr, Setting *setup, const string directory)" << endl; 
   tictoc tim;
   tim.tic();
#endif
    // read Lucias file to get station names, source names and epochs    
    vector<string> station_names, source_names;
    std::string line;
    std::ifstream in_stream;
    int y,m,d,h,min;
    double sec,tau;
    std::string sta1,sta2,src,type;
    vector<ivg::Date> dates;
    while(ivg::parser::get_line(path,in_stream, line))
    {
        std::stringstream tokenizer(line);
        tokenizer>>y>>m>>d>>h>>min>>sec>>sta1>>sta2>>src>>type>>tau;
        
        if(find(station_names.begin(),station_names.end(),sta1)==station_names.end())
            station_names.push_back(sta1);
        if(find(station_names.begin(),station_names.end(),sta2)==station_names.end())
            station_names.push_back(sta2);
        if(find(source_names.begin(),source_names.end(),src)==source_names.end())
            source_names.push_back(src);
        dates.push_back(ivg::Date(y,m,d,h,min,sec));
    }
    
   session_ptr->_start =  dates.at(0);
   session_ptr->_end =  *(dates.end()-1);
   
   // set number of observations 
   session_ptr->_nobs = dates.size();
   
   // create TRF
   ivg::Date date = session_ptr->_start;   
   for(auto &st: station_names)
        replace_string_in_place( st , " ", "_" );
   
   session_ptr->_trf =  ivg::Trf( (*setup), station_names, ivg::staname::ivs_name, true, date.add_days(-10.0), date.add_days(10.0) );
   
    // check if _trf correspondes to number of stations in original vgosdb
    if(session_ptr->_trf.get_number_stations() != station_names.size())
        throw runtime_error( "void Session_inout::_read_vgosdb(ivg::Session *session_ptr, Setting *setup, const string directory): _trf not correctly initialized. Wrong number of stations.");
    
   // shows the main feautres of the loaded trf
   session_ptr->_trf.log_data_info_table();
   
   //create CRF
   vector<ivg::Source> sources;
   for(auto &src: source_names)
       sources.push_back(ivg::Source(ivg::srcname::ivs, src));

   session_ptr->_crf = ivg::Crf( *setup, sources );
   
   //Initializing parameter list (Stations, Sources, EOPs) including epoch and aprioris
   session_ptr->_param_list = Param_list(session_ptr->_trf, session_ptr->_crf, session_ptr->_eops, ivg::Date( 0.5*( session_ptr->_start.get_double_mjd()+session_ptr->_end.get_double_mjd() ) ));
   
    // read Lucias file ones again to create scans and observations
    int obs_counter = 0;
    in_stream.close();
    while(ivg::parser::get_line(path,in_stream, line))
    {
        std::stringstream tokenizer(line);
        tokenizer>>y>>m>>d>>h>>min>>sec>>sta1>>sta2>>src>>type>>tau;
        
        // get source observerd in scan
        ivg::Source* source;
        session_ptr->_crf.get_source(&source,src);

        //set crf2trf matrix for each scan
        ivg::Partials_t2c foo{ivg::Matrix(3,3,0.0)};
        ivg::Partials_t2c * deriv_ptr = &foo;
        ivg::Matrix crf2trf = session_ptr->_eops.form_crf2trf(dates.at(obs_counter),true,deriv_ptr);

        ivg::Scan scan(dates.at(obs_counter),source,crf2trf.transpose(),deriv_ptr);

        // get Analysis_station pointer
        ivg::Analysis_station * station1;
        session_ptr->_trf.get_station(&station1,sta1);

        ivg::Analysis_station * station2;
        session_ptr->_trf.get_station(&station2,sta2);

        int sta1_idx = scan.add_sta(station1);
        int sta2_idx = scan.add_sta(station2);

        ivg::Obs obs_new(session_ptr,&scan,sta1_idx,sta2_idx);

        obs_new.set_delay(0.0,0.0,0.0);
        obs_new.set_cable_cal(0.0,0.0);
        scan.add_scan_meteorology(sta1_idx,-999.99,-999.99,-999.99,0);
        scan.add_scan_meteorology(sta2_idx,-999.99,-999.99,-999.99,0);
        obs_new.set_ion_corr( 666.66 , 1e-12, 0.0 , 0.0 );
        obs_new.set_use_flag( false );

        scan.add_obs( obs_new );       

        //add scan to scan vector
        session_ptr->_scans.push_back( scan );

        obs_counter++;
    }   
#if DEBUG_VLBI >=2 
   cerr << "--- Session_inout::_read_utas_apriori(ivg::Session *session_ptr, Setting *setup, const string directory)" 
        << " : " << tim.toc() << " s " << endl; 
#endif
}

} /* ivg namespace */
