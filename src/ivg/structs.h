/*
 * File:   structs.h
 * Author: iddink
 *
 * Created on 24. Februar 2015, 11:30
 */

#ifndef STRUCTS_H
#define	STRUCTS_H

#include "matrix.h"
#include "ivg_const.h"

namespace ivg
{

struct Antenna
{
    string focus_type;                 // FO_PRIM: primary focus,
    // FO_SECN: secondary focus
    string mounting_type;
    string radome;
    string measure_type;
    double ref_temp;                   // reference temperature [°C]
    double sin_temp;
    double cos_temp;
    double ref_press;
    double ant_diameter;
    double found_height;               // h_f [m]
    double found_depth;                // h_d [m]
    double found_coeff;                // gamma_f; gf [1/K]
    double fixed_axis_length;          // h_p [m]]
    double fixed_axis_coeff;           // gamma_a; ga [1/K]
    double axis_offset_length;
    double axis_offset_coeff;
    double vertex_height;              // h_v [m]
    double struct_coeff;
    double focus_height;               // h_s [m]
    double focus_coeff;
    
    // new antenna information set from antenna.cat
    string id;
    string axis_type;
    string num_id;
    double offset; // m
    double diameter; // m
    double azi_rate; // rad/sec
    double azi_const; // sec
    double azi_min; // rad
    double azi_max; // rad
    double ele_rate; // rad/sec
    double ele_const; // sec
    double ele_min; // rad
    double ele_max; // rad
    
    // we store some sked-file-writing information in this struct
    string A_line; 
    string P_line;
    string H_line;
};

struct Equip
{
    string id;
    string dat_name;
    string head;
    string tape_length;
    map< ivg::band, vector<double> > sefd;
    string rack;
    string recorder;
    string line; // almost complete line for easy skd-file writing
};

struct Frequency
{
    // from freq.cat
    string band; // X or S
    string pol; // R or L
    double frequency; // 8210.99
    string sideband; // U or L
    string channel; // CH1
    double phase_cal_freq; // 10000.0
    
    // from loif.cat, relation by bbc-ID, e.g. 1 to 14
    string if_channel; // 1N
    double lo_frequency; // 8080
    string lo_sideband; // U
    
    // from tracks.cat, relation by channel-ID, CH1, not bbc-ID!
    string track_info; // 1(1,15)
    string fanout_fac; // 1:2
};

struct Channels
{
    string freq_sq_name; // e.g. GEOSX8N, frequency sequence name, as found in freq.cat.
    string freq_lc; // frequency sequence name shot version, e.g. 8F
    string rec_format; // from rec.cat, e.g. Mk34 for WETTZELL
    string rec_hdpos_format; // from rec.cat, e.g. MK3V-A for WETTZELL
    string loif_name; // from rx.cat, set during parse of loif_cat, e.g. CDP_WIDE for WETTZELL
    string hdpos_lines; // from hdpos.cat, including "\n"-delimiter (in general 3 lines)
    
    map<int, int> cha_bbc; // cross reference table between channel(cha)-ID and bbc-ID
    map<int, ivg::Frequency> freq_info; // bbc-ID <-> Frequency-struct  
};

struct Wave
{
    string name;
    double phase, freq, accel;
    double up_cos, east_cos, north_cos;
    double up_sin, east_sin, north_sin;
};


struct Partials_t2c
{
    Matrix ut1;
    Matrix pmx;
    Matrix pmy;
    Matrix pnx;
    Matrix pny;
};

// struct for turbulence data
struct turbulence_data
{
   double Cn;       // refractive index structure constant [m^(-1/3)] 
   double h;        // effective height of the troposphere [m]
   double v;        // wind velocity [m/s] (for Matern model)
   ivg::Matrix vel; // wind vector [m/s] --> [m/h] (for Onsala and SIGMA-C model)
   double v_dir;    // wind direction [rad]; pi == horizontal; 0.0 == zenith (for Onsala and SIGMA-C model)
   double zwd0;     // initial zenith wet delay [mm]
   double dh_seg;   // time segments over which observations are to be correlated [h]
   double dh;       // height increment for numeric integration [m]
};

struct Psd
{
    double mjd;
    int e_mode, n_mode, u_mode;
    double e_a1, e_a2, e_t1, e_t2;
    double n_a1, n_a2, n_t1, n_t2;
    double u_a1, u_a2, u_t1, u_t2;
};

} // namespace ivg


#endif	/* STRUCTS_H */

