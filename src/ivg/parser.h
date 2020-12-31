#ifndef IVGPARSER_H
#define IVGPARSER_H

//#include "station.h"
//#include "analysis_station.h"
#include "logger.h"
#include "matrix.h"
#include "date.h"
#include "auxfunc.h"
#include "trf.h"
#include "structs.h"
#include "source.h"
#include "crf.h"
#include "sinex.h"
#include "analysis_station.h"
#include "norad.h" //used for TLE elements for moon and satellites


// ---------------------------------------------------------------------------
namespace ivg
{

class Trf; // forward declaration
class Analysis_station; // forward declaration

enum mode_i {
  NONE = 0,
  EST2TRF = 1<<0,
  EST2CRF = 1<<1,
  EST2EOP = 1<<2,
  APR2TRF = 1<<3,
  APR2CRF = 1<<4,
  APR2EOP = 1<<5,
  TRF2APR = 1<<6,
  CRF2APR = 1<<7,
  EOP2APR = 1<<8,
  NUT2ZER = 1<<9,
};

 enum extdata { MAPPING, HYDROSTATIC, MAPPING3 };

namespace parser
{    
    //Used in trf.h
    void optl(ivg::Trf * trf_ptr, const string path);
    void blq(ivg::Trf * trf_ptr, const string path);
    void ecc(ivg::Trf * trf_ptr, const string path);
    void gravdef(ivg::Trf * trf_ptr, const string path, ivg::Date ep);
    void bindisp(ivg::Trf * trf_ptr, const string folder_path, map<string,string> correspondence, ivg::Date start, ivg::Date end);
    void hps(ivg::Trf * trf_ptr, const string path);
    void dat(ivg::Trf * trf_ptr, const string path);
    void antenna_info(ivg::Trf * trf_ptr, const string path);
    void external_met_data( ivg::Trf * trf_ptr, const string path, ivg::Date start,ivg::Date end, ivg::extdata type );
    void hydlo(ivg::Trf * trf_ptr, const string folder_path);
    void psd_coefficients(ivg::Trf * trf_ptr, const string path);
    void raytraced_delays( ivg::Trf * trf_ptr, const string path );
    
    //used in eop_series.h
    ivg::Matrix c04(const string path, ivg::Date start, ivg::Date end);
    ivg::Matrix cs_erp(const string path, ivg::Date start, ivg::Date end);
    ivg::Matrix eops(const string path, ivg::Date start, ivg::Date end);
    ivg::Matrix finals(const string path, ivg::Date start, ivg::Date end);
    ivg::Matrix igs_erp(const string path, ivg::Date start, ivg::Date end);
    
    // methods for sked catalogs
    void stations_cat(ivg::Trf * trf_ptr,  const string path);
    void antenna_cat(ivg::Trf * trf_ptr,  const string path);
    void equip_cat(ivg::Trf * trf_ptr,  const string path);
    void flux_cat(ivg::Crf * crf_ptr,  const string path, bool isSkedFile = false);
    void flux_cat_line(ivg::Crf * crf_ptr, std::string line);
    string freq_cat(ivg::Trf * trf_ptr,  const string path, string freqname);
    void loif_cat(ivg::Trf * trf_ptr,  const string path, const string rx_path, string rx_name);
    void tracks_cat(ivg::Trf * trf_ptr,  const string path, const string rec_path, string rec_mode);
    void hdpos_cat(ivg::Trf * trf_ptr,  const string path);
    void mask_cat(ivg::Trf * trf_ptr,  const string path);
    
    // used in crf.cpp for initializing tle elements (for satellites and moon)
    void tle(ivg::Crf * crf_ptr,  const string path);
    // used in crf.cpp for initializing sp3-sat orbits from IGS
    void sp3(ivg::Crf * crf_ptr,  const string path);
    
    // used in trf.cpp and masterfile.cpp
    vector< map<ivg::staname,string> >  nscodes_parser(const vector<string>station_names, const ivg::staname type, const string path);
    vector<ivg::Analysis_station> ssc_parser(const string path, const vector< map<ivg::staname,string> > &nscodes);
    map<string,string> correspondence(const string path);

    istream& get_line(const string &path, ifstream  &inStream, string &line);
    
    int init_options(string option_str);
}
}
// ...........................................................................

#endif
