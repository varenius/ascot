#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf
#include <vector>
#include <boost/algorithm/string.hpp>
#include "date.h"
#include <libconfig.h++>
#include <QApplication>
#include <QColor>
#include "matrix.h"
#include "date.h"
#include "lapack_wrapper.h"
#include "tictoc.h"
#include "session.h"
#include "scan.h"
#include "ivg_const.h"
#include "ls_neq.h"
#include "logger.h"
#include "session_inout.h"
#include "simulation.h"
#include <tclap/CmdLine.h>
#include "masterfile.h"
#include "wrapper.h"

#include <cstdlib>
#include <libconfig.h++>
#include "statistics.h"
#include "plot.h"
#include "ascot.h"
#include "tsa.h"
#include "lp_sked_gui.h"
#include "lp_sked_station.h"
#include "lp_sked_worker.h"
#include "transits.h"

#include "milp.cpp"
#include "grid.h"
#include "geometry.h"
#include "definitions.h"
#include "vgosdb.h"
#include "wrapper.h"
#include "masterfile.h"
#include <curl/curl.h>

using namespace libconfig;

using namespace std;
loglevel g_verbose;
void * g_ephis;

bool get_log_date(std::string str,ivg::Date *tim)
{
  std::stringstream ss(str);
  int yy,dd,hh,mi;
  double sec;
  if (!(ss >> yy) ) return false;
  ss.ignore();
  if (!(ss >> dd) ) return false;
  ss.ignore();
  if (!(ss >> hh) ) return false;
  ss.ignore();
  if (!(ss >> mi) ) return false;
  ss.ignore();
  if (!(ss >> sec) ) return false;
  
  *tim=ivg::Date(yy,1,dd,hh,mi,sec);;
  return true;
}

void read_logfile(std::string file, vector<ivg::Date> *cab_date,vector<double> *cable,vector<ivg::Date> *met_date, vector<double> *pres,vector<double> *temp,vector<double> *hum)
{
  std::ifstream infile(file);
  std::string line;
  std::vector<double> cablelong;
  int sign=0;
  std::vector<ivg::Date> long_date;
  while (std::getline(infile, line))
    {
      std::size_t id;
      boost::to_lower(line);
      id=line.find("/wx/");
      if (id !=std::string::npos)
	{
	  std::vector<double> tmp;
	  std::stringstream ss(line.substr(id+4));
	  ivg::Date tim;
	  for (double i;ss>>i;)
	    {
	      tmp.push_back(i);
	      if (ss.peek()==',')
		ss.ignore();
	    }
	  if ((tmp.size()>=3)&&( get_log_date(line,&tim))) {
	    met_date->push_back(tim);
	    pres->push_back(tmp[1]);
	    temp->push_back(tmp[0]);
	    hum->push_back(tmp[2]);
	  }
	    
	}
      id=line.find("/cable/");
      if (id !=std::string::npos)
	{
	  std::stringstream ss(line.substr(id+7));
	  double i;
	  ivg::Date tim;
	  if ((ss>>i)&&( get_log_date(line,&tim))){
	    cab_date->push_back(tim);
	    cable->push_back(i);
	  }
	}
     
      id=line.find("/cablelong/");
      if (id !=std::string::npos)
	{
	  std::stringstream ss(line.substr(id+11));
	  double i;
	  ivg::Date tim;
	  if ((ss>>i)&&( get_log_date(line,&tim))){
	    long_date.push_back(tim);
	    cablelong.push_back(i);
	  }
	}
    }
  for (int i=0;i<cablelong.size();i++)
    {
      for (int j=0;j<cable->size()-1;j++)
	{
	  if ((cab_date->at(j)<=long_date.at(i))&&(cab_date->at(j+1)>long_date.at(i)))
	    {
	      if (cable->at(j)>cablelong.at(i))
		sign-=1;
	      else
		sign+=1;
	    }
	}
    }
  if (sign<0)
    {
      for (int j=0;j<cable->size();j++)
	cable->at(j)=-cable->at(j);
    }
	     
    
}

void calc_eff_iono_freqs(std::vector<double> *fig,std::vector<double> *fip,std::vector<double> *fir, ivg::Vgosdb *vgosdb, std::string band_str)
{
  vector<vector<vector<double> > > tst;
  std::cout << "read stuff" << endl;
  tst=vgosdb->get_vector_3d_data<double>("Observables","ChannelInfo_b"+band_str,"ChanAmpPhase");
 
  ivg::Matrix chanAmp = ivg::Matrix(tst.size(),tst[0].size(),0.0);
  ivg::Matrix chanPhas = ivg::Matrix(tst.size(),tst[0].size(),0.0);
  
  for (int i=0;i<tst.size();i++) {
    for (int j=0;j<tst[i].size();j++) {
      chanAmp(i,j)=tst[i][j][0];
      chanPhas(i,j)=tst[i][j][1];
    }
	      
  }
 
  int bitsampl;
  if (vgosdb->does_variable_exist("Observables","ChannelInfo_b"+band_str,"BITSAMPL"))
    bitsampl = vgosdb->get_scalar<int>("Observables","ChannelInfo_b"+band_str,"BITSAMPL");
  else
    bitsampl = 1;
  int sampl_rate = vgosdb->get_scalar<int>("Observables","ChannelInfo_b"+band_str,"SampleRate");
  ivg::Matrix chanfreq = vgosdb->get_matrix("Observables","ChannelInfo_b"+band_str,"ChannelFreq");
  vector<vector<vector<int> > > numAp = vgosdb->get_vector_3d_data<int>("Observables","ChannelInfo_b"+band_str,"NumAp");
  vector<double> RefFreq = vgosdb->get_vector<double>("Observables","RefFreq_b"+band_str,"RefFreq");
  double bw=double(sampl_rate)/(2*bitsampl)*1e-6;
  std::cout << "loop" << endl;
  for (int i=0;i<tst.size();i++)
    {
     
      double sumr=0,sumrf=0,sumrf2=0,sumrdf=0,sumrf2p=0,sumrfp=0,sumrfdfp=0;
      double f0;
      if (RefFreq.size()==1)
	f0=RefFreq[0];
      else
	f0=RefFreq[i];
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
      
      if (sumr>0){
	fig->push_back(sqrt(-(sumr*sumrf2-sumrf*sumrf)/(sumr*sumr-sumrf*sumrdf)));
	fip->push_back(sqrt(f0*(sumrf2p*sumr-sumrfp*sumrfp)/(sumrf2p*sumrdf-sumrfp*sumrfdfp)));
	fir->push_back(sqrt(sumrf2/sumr));
      } else {
	fig->push_back(0);
	fip->push_back(0);
	fir->push_back(0);
      }
    }
  
}

// ...........................................................................
size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream)
// ...........................................................................
{
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

int main( int argc, char *argv[])
{
      TCLAP::CmdLine cmd("logreader", ' ', "0.1");
    
    // verbose level
    TCLAP::ValueArg<int> verboseArg( "v","verbose","Verbose Level (0=NOTHING, 1=INFO, 2=DETAIL, 3=RESULT, 4=WARNING, 5=ALL)",false,4,"int");
    cmd.add( verboseArg );
  
    // controlfile
    TCLAP::ValueArg<std::string> cntArg( "c","controlfile","/home/ascot/cnt",false,"../src/progs/logreader/logreader.cfg","string");
    cmd.add(cntArg);  
        
    // dbname
    TCLAP::ValueArg<std::string> dbArg( "e","explicitDB","e.g. 93NOV05XU",false,"","string");       
    cmd.add(dbArg);

    // parts
    TCLAP::ValueArg<std::string> partArg( "s","specify_files","e.g. emcf (edit, met, cable, iono freq)",false,"emcf","string");       
    cmd.add(partArg);

    // Parse the argv array.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg. 
    g_verbose = (loglevel) verboseArg.getValue();
    std::string parts = partArg.getValue();
    
    std::string controlfile = cntArg.getValue();
    libconfig::Config cfg;
    try
    {
       cfg.readFile( controlfile.c_str() );
    }
    catch( libconfig::ParseException & err )
    {
        cerr << "libconfig::" << err.what() << " in " << err.getFile() << " at line " << err.getLine() << endl;
        exit( -1 );
    }
    libconfig::Setting& setup= cfg.lookup( "setup" );
    std::vector<std::string> database;
    
    if (dbArg.isSet()) {
      // -e option : sessions are parsed from argv
      database.push_back(dbArg.getValue());
        
    } else {
      database.reserve( setup[ "sessions" ].getLength() );
      
      for( int i=0; i<setup[ "sessions" ].getLength(); ++i )
        {
            database.push_back( setup[ "sessions" ][ i ][ "dbname" ] );
	    
        }  
    }    
    for( int i=0; i <database.size(); ++i ){
      std::string db=database.at(i);
      std::string logfile;
      ivg::Vgosdb vgosdb;

      int year;
      if (db.size()<=9){
	year = stoi(db.substr(0,2));
	if (year<79)
	  year += 2000;
	else
	  year += 1900;
      } else {
	year=stoi(db.substr(0,4));
      }
 
      stringstream folderpath;
      string session_dir = (const char *)get_list_element(setup["datadirs"],setup["session_type"])[2];
      folderpath<<session_dir<<year<<"/"<<db<<"/";
      ivg::Wrapper wrapper(folderpath.str(), db, "iIVS");
      
      vgosdb = ivg::Vgosdb(folderpath.str(),&wrapper);
      string editpath=folderpath.str()+"/ObsEdit";
      mkdir(editpath.c_str(),0777);
      string derivedpath=folderpath.str()+"/ObsDerived";
      mkdir(derivedpath.c_str(),0777);
      std::string head_filename = wrapper.get_file(ivg::wrapper_entries::Head,ivg::band::X); 
      vector<string> station_names = vgosdb.get_string("",head_filename,"StationList");
  
      for(auto &st: station_names)
	replace_string_in_place( st , " ", "_" );
      string expname="";
      string headfile=folderpath.str()+"/"+head_filename+".nc";
      NcFile file(headfile.c_str(), NcFile::ReadOnly);
      NcVar *tmp_var;
  
      tmp_var=file.get_var("Session");
      
      int n = tmp_var->get_dim(0)->size();
      std::cout << n <<endl;
      char values[n];
      int nobs = vgosdb.get_scalar<int>("",head_filename,"NumObs");
      tmp_var->get(&values[0],n);
      for (int i=0;i<n;i++)
	expname+=values[i];
      std::cout << expname << endl;
      
      vector< map<ivg::staname,string> > nscodes = ivg::parser::nscodes_parser( station_names, ivg::staname::ivs_name,"apriori_files/ns-codes.txt");
      boost::to_lower(expname);
      std::string nscode=nscodes.at(1)[ivg::staname::lettercode];
      boost::to_lower(nscode);
      if ((parts.find("c")!=std::string::npos)||(parts.find("m")!=std::string::npos)) {
	CURLcode res;
	FILE *fp;
	CURL *curl;
	for (int i=0;i<nscodes.size();i++)
	  {
	    std::cout << "Met and cable data " << nscodes.at(i)[ivg::staname::lettercode] <<endl;
	    vector<double> cable,pres,temp,hum;
	    vector<ivg::Date> cab_date,met_date;
	    std::string nscode=nscodes.at(i)[ivg::staname::lettercode];
	    boost::to_lower(nscode);
	    logfile = "logfiles/"+expname+nscode+".log";
	    std::cout << " log: "<< logfile << endl;
	    
	    std::cout << "download log-file" << endl;
	    curl = curl_easy_init();
	    if (curl)
	      {
	      fp = fopen(logfile.c_str(),"wb");
	      stringstream cddis;
	      cddis<<"ftp://ivsopar.obspm.fr/vlbi/ivsdata/aux/"<<year<<"/"<<expname<<"/"<<expname << nscode << ".log";
		  curl_easy_setopt(curl, CURLOPT_URL,cddis.str().c_str());
		  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
		  curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
		  res = curl_easy_perform(curl);
		  curl_easy_cleanup(curl);
		  fclose(fp);
	      }
	    std::cout << "Read log" << endl;
	    read_logfile(logfile,&cab_date,&cable,&met_date,&pres,&temp,&hum);
	    
	    ivg::Matrix time = vgosdb.get_matrix(station_names.at(i),"TimeUTC","YMDHM");
	    vector<double> seconds_vec = vgosdb.get_vector<double>(station_names.at(i),"TimeUTC","Second");
	    time.append_cols(ivg::Matrix(seconds_vec));
	    
	    vector<double> cable2,pres2,temp2,hum2;
	    int cur_cab_id=0,cur_met_id=0; 
	    std::cout << "Interpolate" << endl;
	    for(int j=0; j<time.rows(); j++){
	      if( time(j,0) < 1900 ){
		if ( time(j,0) < 79 )
		  time(j,0) += 2000;
		else
		  time(j,0) += 1900;
	      }
	      ivg::Date obs_epoch(ivg::Date(time(j,0),time(j,1),time(j,2),time(j,3),time(j,4),time(j,5)));
	      if (cable.size()>0) {
		while ((cur_cab_id<cable.size()-1)&&(cab_date.at(cur_cab_id)<obs_epoch))
		  cur_cab_id++;
		if (cur_cab_id==0)
		  cable2.push_back(2.5e-6*cable.at(0));
		else if (cur_cab_id==cable.size())
		  cable2.push_back(2.5e-6*cable.at(cur_cab_id-1));
		else if (cable.size()!=0)
		  {
		    cable2.push_back(2.5e-6*(cable.at(cur_cab_id-1)*(cab_date.at(cur_cab_id).get_double_mjd()-obs_epoch.get_double_mjd())+cable.at(cur_cab_id)*(-cab_date.at(cur_cab_id-1).get_double_mjd()+obs_epoch.get_double_mjd()))/(cab_date.at(cur_cab_id).get_double_mjd()-cab_date.at(cur_cab_id-1).get_double_mjd()));
		  }
	      }
	      if (pres.size()>0) {
		while ((cur_met_id<pres.size()-1)&&(met_date.at(cur_met_id)<obs_epoch))
		  cur_met_id++;
		if (cur_met_id==0){
		  pres2.push_back(pres.at(0));
		  temp2.push_back(temp.at(0));
		  hum2.push_back(hum.at(0));
		}
		else if (cur_met_id==pres.size()){
		  pres2.push_back(pres.at(cur_met_id-1));
		  temp2.push_back(temp.at(cur_met_id-1));
		  hum2.push_back(hum.at(cur_met_id-1));
		}
		else
		  {
		    pres2.push_back((pres.at(cur_met_id-1)*(met_date.at(cur_met_id).get_double_mjd()-obs_epoch.get_double_mjd())+pres.at(cur_met_id)*(-met_date.at(cur_met_id-1).get_double_mjd()+obs_epoch.get_double_mjd()))/(met_date.at(cur_met_id).get_double_mjd()-met_date.at(cur_met_id-1).get_double_mjd()));
		    temp2.push_back((temp.at(cur_met_id-1)*(met_date.at(cur_met_id).get_double_mjd()-obs_epoch.get_double_mjd())+temp.at(cur_met_id)*(-met_date.at(cur_met_id-1).get_double_mjd()+obs_epoch.get_double_mjd()))/(met_date.at(cur_met_id).get_double_mjd()-met_date.at(cur_met_id-1).get_double_mjd()));
		    hum2.push_back((hum.at(cur_met_id-1)*(met_date.at(cur_met_id).get_double_mjd()-obs_epoch.get_double_mjd())+hum.at(cur_met_id)*(-met_date.at(cur_met_id-1).get_double_mjd()+obs_epoch.get_double_mjd()))/(met_date.at(cur_met_id).get_double_mjd()-met_date.at(cur_met_id-1).get_double_mjd()));
		  }
	      }
	      
	    }
	    if (parts.find("m")!=std::string::npos) {
	      std::cout << "Write Met" << endl;
	      if (pres2.size()>0)
		vgosdb.create_met_file(temp2,pres2,hum2, station_names.at(i), expname);
	    }
	     if (parts.find("c")!=std::string::npos) {
	       std::cout << "Write Cable" << endl;
	       if (cable2.size()>0)
		 vgosdb.create_cal_file(cable2, station_names.at(i), expname);
	     }
	  }
      }
      if (parts.find("e")!=std::string::npos) {
	std::cout << "reading Q-codes" <<endl;
	vector<char> qcodex,qcodes;
	std::string qcode_filename;
	if( wrapper.file_exists(ivg::wrapper_entries::QualityCode,ivg::band::X) )
	  qcode_filename = wrapper.get_file(ivg::wrapper_entries::QualityCode,ivg::band::X);
	else
	  qcode_filename = "QualityCode_bX";
	if(vgosdb.does_file_exist("Observables",qcode_filename))
	  qcodex = vgosdb.get_vector<char>("Observables",qcode_filename,"QualityCode");
	else
	  qcodex = vector<char>(nobs,'9');
	if( wrapper.file_exists(ivg::wrapper_entries::QualityCode,ivg::band::S) )
	  qcode_filename = wrapper.get_file(ivg::wrapper_entries::QualityCode,ivg::band::S);
	else
	  qcode_filename = "QualityCode_bS";
	if(vgosdb.does_file_exist("Observables",qcode_filename))
	  qcodes = vgosdb.get_vector<char>("Observables",qcode_filename,"QualityCode");
	else
	  qcodes = vector<char>(nobs,'9');
	std::cout << "calculating flags" <<endl;
      
	std::vector<int> flag;
	for (int i=0;i<nobs;i++){
	  
	  if ((qcodex[i]>='5')&&(qcodes[i]>='5')&&(qcodex[i]<='9')&&(qcodes[i]<='9'))
	    flag.push_back(0);
	  else
	    flag.push_back(5);
	  
	}
	std::cout << "Writing edit file" <<endl;
	vgosdb.create_Edit_file( flag, flag, flag,expname,"" );
      }
      
      if (parts.find("f")!=std::string::npos) {
	vector<double> fig,fip,fir;
	std::cout << "ion freqs X" << endl;
	calc_eff_iono_freqs(&fig,&fip,&fir, &vgosdb,  ivg::band_to_string(ivg::band::X));
	std::cout <<"write file" <<endl;
	vgosdb.create_EffFreq_file(fig,fip,fir,expname,ivg::band::X,"");
	vector<double> figs,fips,firs;
      
	std::cout << "ion freqs S" << endl;
	calc_eff_iono_freqs(&figs,&fips,&firs, &vgosdb,  ivg::band_to_string(ivg::band::S));
	std::cout <<"write file" <<endl;
	vgosdb.create_EffFreq_file(figs,fips,firs,expname,ivg::band::S,"");
	std::cout << fig.size() <<" " << figs.size()<< endl;
      }
      wrapper.write_wrapper((const char *)setup["vgosdb_editing"],3 );
    }
}
