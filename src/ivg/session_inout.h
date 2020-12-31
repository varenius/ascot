#ifndef SESSION_INOUT_H
#define	SESSION_INOUT_H


#include "session.h"
#include "trf.h"
#include "crf.h"
#include "source.h"
#include "scan.h"
#include "analysis_station.h"
#include "sys/utsname.h"
#include "auxfunc.h"
#include "param.h"
#include "vgosdb.h"
#include "tictoc.h"
#include "masterfile.h"
#include "wrapper.h"

#include <cstdlib>
#include <libconfig.h++>


namespace ivg
{
  
class Session_inout {
    
public:
    
    Session_inout();
    
    Session_inout(const string type, ivg::Masterfile masterfile = ivg::Masterfile());
       
    std::string load(ivg::Session *session_ptr, Setting *setup,  const string name, const string version = "");
    
    ivg::Matrix init_xka_catalog(ivg::Session *session_ptr, Setting *setup,  const string path);
  void write_snx(ivg::Session *session_ptr, string outfile) {
    write_snx(session_ptr, outfile, false);
  };
  void write_snx(ivg::Session *session_ptr, string outfile, bool incl_results);
    
    void write_snx_tro(ivg::Session *session_ptr, string outfile );
    
    /**
    *  \b Description: \n
    *        Method to export the results of the estimation process. This includes:
    *        (a) parameter index
    *        (b) parameter name
    *        (c) parameter type
    *        (d) polynomial degree
    *        (e) epoch (MJD)
    *        (f) a priori value
    *        (g) estimate
    *        (h) standard deviation
    *  \param [in] no input parameters needed
    */
    
    void write_results(ivg::Session *session_ptr, string outfile, bool write_prediction = false );
    void read_results(ivg::Session *session_ptr,string infile, bool apriori = true);
    void write_residuals(ivg::Session *session_ptr, string outfile,
                         std::string resid_format = "ASCOT");
    void write_outliers(ivg::Session *session_ptr, string outfile );
    void write_groupdelay( ivg::Session *session_ptr );
    void write_obs_met_info( ivg::Session *session_ptr, string outfile );
    void write_sou(ivg::Session *session_ptr, string outfile );
    void write_skd(ivg::Session *session_ptr, string outfile);
    void setWrapper_ptr(ivg::Wrapper* _wrapper_ptr);
    
private:

    void _write_snx_matrix( ofstream &outstream, ivg::Matrix matrix );
    void _write_snx_vector( ofstream &outstream, ivg::Matrix matrix );
    void _dec_lat_2_grad_min_sec(double coord, int &latGrad, int &latMin, double &latSec);
    void _dec_lon_2_grad_min_sec(double coord, int & lonGrad, int & lonMin, double & lonSec);
    
    
    void _read_ngs(ivg::Session *session_ptr, Setting *setup, const string path);
    void _read_vgosdb(ivg::Session *session_ptr, Setting *setup, const string directory);
    void _read_snx(ivg::Session *session_ptr, Setting *setup, const string path);
    void _read_sou(ivg::Session *session_ptr, Setting *setup, const string path);
    void _read_utas_apriori(ivg::Session *session_ptr,Setting *setup,const string path);
    void _init_skd_session(ivg::Session *session_ptr,Setting *setup,const string dir);
    void _read_skd(ivg::Session *session_ptr, Setting *setup, const string path);
        
    string _DtoE(string str);
    
    ivg::Masterfile _masterfile;
    ivg::Wrapper* _wrapper_ptr;
    string _type;
    string _session_path; // e.g. /data/vgosDB/2002/02OCT25-A/ in case of _type = vgosdb
    
};

  
} /* ivg namespace */

#endif	/* SESSION_INOUT_H */

