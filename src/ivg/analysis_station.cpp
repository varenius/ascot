#include <vector>
#include "analysis_station.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Analysis_station::Analysis_station()
// ...........................................................................
{

}

// ...........................................................................
Analysis_station::Analysis_station( ivg::Matrix pos0, ivg::Matrix vel0,
                                    ivg::Date refepoch,
                                    std::vector<ivg::Date> discontinuity,
                                    std::map<ivg::staname, std::string> names )
// ...........................................................................
{
    _pos0 = pos0 ;
    _vel0 = vel0 ;
    
    // checking if sensible aprioris are set
    if(_pos0(0) != 0.0 && _pos0(1) != 0.0 && _pos0(2) != 0.0)
        _data_status["XYZ"] = "X";
    
    _pos0_std = ivg::Matrix(_pos0.rows(),_pos0.cols(),0.0) ;
    _vel0_std = ivg::Matrix(_vel0.rows(),_vel0.cols(),0.0) ;
    _refepoch = refepoch ;
    _discontinuity = discontinuity;
    _names = names;
    _last_crf2trf[ivg::Date(1970,1.0)]=ivg::Matrix(3,3,0.0);
}

// ...........................................................................
Analysis_station::Analysis_station( ivg::Matrix pos0, ivg::Matrix vel0,
                                    ivg::Matrix pos0_std, ivg::Matrix vel0_std, 
                                    ivg::Date refepoch,
                                    std::vector<ivg::Date> discontinuity,
                                    std::map<ivg::staname, std::string> names ):Analysis_station( pos0, vel0, refepoch, discontinuity, names )
// ...........................................................................
{
    _pos0_std = pos0_std;
    _vel0_std = vel0_std;
}


// ...........................................................................
Analysis_station::~Analysis_station() {}
// ...........................................................................



// ===========================================================================
// public methods
// ===========================================================================

// getter and setter
void Analysis_station::set_hydlo( ivg::Matrix hydlo )
{
    _hydlo = hydlo;
    _data_status["HYDLO"] = "X";
}

// ...........................................................................
void Analysis_station::set_nontidal_aplo( ivg::Matrix aplo )
// ...........................................................................
{
    if( aplo.rows() == 4 && aplo.cols() > 40)
    {
        _aplo_xyz = aplo;
        _data_status["NTAPL"] = "X";
    }
    // if there is to less data, the interpolation doesn't work
    else
    {
        log<WARNING>("!!! Not enough NTAPL data for ") % _names[staname::ivs_name] % ". No data set.";
        _data_status["NTAPL"] = "E";
    }
}

// ...........................................................................
void Analysis_station::set_tidal_aplo( map<string,ivg::Wave> waves )
// ...........................................................................
{
    if(!waves.empty())
    {
        _aplo_tidal = waves;
        _data_status["TAPL"] = "X";
    }
    else
    {
        log<WARNING>("!!! Dimensions of _aplo_tidal map not correct for ") % _names[staname::ivs_name];
        _data_status["TAPL"] = "E";
    }
}

// ...........................................................................
void Analysis_station::set_optl_coeff( ivg::Matrix optlc )
// ...........................................................................
{
    if( optlc.rows() == 3 && optlc.cols() == 2 )
    {
        _optl_coeff = optlc;
        _data_status["OPTL"] = "X";
    }
    else
    {
        log<WARNING>("!!! Dimensions of _optl_coeff matrix not correct for ") % _names[staname::ivs_name];
        _data_status["OPTL"] = "E";
    }
}

// ...........................................................................
void Analysis_station::set_ocean_loading_coeff( std::vector<float> olc )
// ...........................................................................
{
    if( olc.size() == 66 )
    {
        _ol_coeff = olc;
        _data_status["OLC"] = "X";
    }
    else
    {
        log<WARNING>("!!! Dimensions of _ol_coeff vector not correct for ") % _names[staname::ivs_name];
        _data_status["OLC"] = "E";
    }
}

// ...........................................................................
  void Analysis_station::set_ocean_loading_coeff_harpos( ivg::Matrix olc )
// ...........................................................................
{
    
        _ol_harpos = olc;
	
        _data_status["OLC"] = "X";
 
}
  
// .........................................................................................
void Analysis_station::set_eccentricity( ivg::Matrix ecc,
        std::vector<ivg::Date> refepoch )
// .........................................................................................
{
    if( ecc.rows() == 3 && ecc.cols() == refepoch.size() )
        _ecc = ecc;
}

// .........................................................................................
void Analysis_station::set_antenna_info(ivg::Antenna antenna_info)
// .........................................................................................
{
    _antenna_info = antenna_info;
    _data_status["ANT_I"] = "X";
}

// .........................................................................................
ivg::Antenna & Analysis_station::set_antenna_info()
// .........................................................................................
{
    // called from ivg::parser::antenna_cat to add further information to _antenna_info
    _data_status["ANT_C"] = "X";
    return(_antenna_info);
}

// .........................................................................................
ivg::Channels & Analysis_station::set_channel_setup()
// .........................................................................................
{
    _data_status["CHA_I"] = "X";
    return(_channel_setup);
}

// .........................................................................................
void Analysis_station::set_vmf1_data( ivg::Matrix vmf1_data )
// .........................................................................................
{
    _vmf1_data = vmf1_data;
    _data_status["VMF1"] = "X";
    
    // if not enough rows, it seems to be wrong and not suitable for interpolation
    if(_vmf1_data.rows() <= 15)
        _data_status["VMF1"] = "W";
}
// .........................................................................................
void Analysis_station::set_vmf3_data( ivg::Matrix vmf3_data )
// .........................................................................................
{
    _vmf3_data = vmf3_data;
    _data_status["VMF3"] = "X";
    
    // if not enough rows, it seems to be wrong and not suitable for interpolation
    if(_vmf3_data.rows() <= 15)
        _data_status["VMF3"] = "W";
}


// .........................................................................................
void Analysis_station::set_zhd_data( ivg::Matrix zhd_data )
// .........................................................................................
{
    _zhd_data = zhd_data;
//    _data_status["ZHD"] = "X";
//    
//    // if only one row, it seems to be wrong
//    if(_zhd_data.rows() == 1)
//        _data_status["ZHD"] = "W";
}

// .........................................................................................
void Analysis_station::set_raytracing_data( std::map< ivg::Date,std::map<std::string,ivg::Matrix> > rt_data )
// .........................................................................................
{
    _raytracing_data = rt_data;
}

// .........................................................................................
ivg::Antenna Analysis_station::get_antenna_info()
// .........................................................................................
{
    return(_antenna_info);
}

// .........................................................................................
ivg::Equip Analysis_station::get_equip_info()
// .........................................................................................
{
    return(_equip);
}
// .........................................................................................
ivg::Matrix Analysis_station::get_eccentricity( ivg::Date epoch )
// .........................................................................................
{    
    // if not valid eccentricity is available, return [0.0 0.0 0.0]
    if(_ecc_epochs.size() == 0 || epoch < _ecc_epochs.at(0))
        return ivg::Matrix(3,1,0.0);
    
    if(_ecc.cols() == 1)
    {
        return(_ecc);
    }
    else
    {
        for(int i=0; i<_ecc.cols()-1; i++)
        {

            if(epoch > _ecc_epochs.at(i) && epoch < _ecc_epochs.at(i+1) )
            {
                return(_ecc.get_col(i));
            }
        }
        return(_ecc.get_col(_ecc.cols()-1));
    }
}

// .........................................................................................
void Analysis_station::add_eccentricity( ivg::Matrix ecc, ivg::Date epoch,
        string epoch_type )
// .........................................................................................
{
    //If ecc not xyz but north, east, up -> transformation!
    if(epoch_type == "NEU")
        ecc = _ren2xyz(ecc,"ner");

    if( ecc.rows() == 3 && ecc.cols() == 1 )
    {
        if( _ecc.rows() == 0 && _ecc.cols() == 0 )
        {
            _ecc = ecc;
            _ecc_epochs.push_back( epoch );
        }
        else
        {
            _ecc.append_cols( ecc );
            _ecc_epochs.push_back( epoch );
        }
        
        _data_status["ECC"] = "X";
    }
    else
    {
        log<WARNING>("!!! Dimensions of _ecc matrix not correct for ") % _names[staname::ivs_name];
        _data_status["ECC"] = "E";
    }
    
    
}
void Analysis_station::set_gravdef( ivg::Matrix gd )
// .........................................................................................
{
    if( gd.cols() == 2 )
    {
      _grav_deform=gd;
      _data_status["GRAVD"] = "X";
    }
    else
    {
      log<WARNING>("!!! Dimensions of _grav_deform matrix not correct for ") % _names[staname::ivs_name];
       _data_status["GRAVD"] = "W";
    }
    
    
}
  
// ...........................................................................
void Analysis_station::add_equip_info( ivg::Equip equip )
// ...........................................................................
{
    // stores all information from equip.cat for scheduling
    _equip = equip;
    
    // check if band is already set up
    if(_sefd[ivg::band::X].size() > 0 || _sefd[ivg::band::S].size() > 0 )
    {
        log<WARNING>("!!! More than one sefd info for ") % get_name(ivg::staname::ivs_name) % ". Ignoring new information.";
        _data_status["SEFD"] = "W";
    }
    else
    {
        _sefd = _equip.sefd;
        _data_status["SEFD"] = "X";
    }
    
}
// ...........................................................................
void Analysis_station::add_band_frequency_sequence( ivg::band band, ivg::Matrix freq_seq )
// ...........................................................................
{
    // check if band is already set up
    if(_freq_sequences[band].rows() > 0 )
    {
        log<WARNING>("!!! More than one frequency sequence for ") % get_name(ivg::staname::ivs_name) % ". Overwriting old frequency sequence for this band";
        _data_status["FREQ"] = "W";
    }
    else
    {
        _freq_sequences[band] = freq_seq;
        _data_status["FREQ"] = "X";
    }
}
// ...........................................................................
void Analysis_station::add_mask_info( ivg::Matrix mask )
// ...........................................................................
{
    // check if band is already set up
    if(_mask.rows() > 0 )
    {
        log<WARNING>("!!! More than one azimut-elevation-mask for ") % get_name(ivg::staname::ivs_name) % ". Overwriting old mask.";
        _data_status["MASK"] = "W";
    }
    else
    {
        _mask = mask;
        _data_status["MASK"] = "X";
    }
}
// ...........................................................................
void Analysis_station::show()
// ...........................................................................
{
    cout << "------------------------------------------------" << endl;
    show_names();
    double d = _refepoch.get_double_mjd();
    cout << "----------------------" << endl;
    cout << " reference epoch (mjd): " << d << endl;
    cout << "----------------------" << endl;
    cout << " position matrix (xyz)" << endl;
    _pos0.show(12);
    cout << "----------------------" << endl;
    cout << " velocity matrix (xyz)" << endl;
    _vel0.show();
    cout << "----------------------" << endl;
    cout << " discontinuities (date)" << endl;
    for(auto &disc: _discontinuity)
        cout << disc.get_date_time("YYYY:DOY:SSSSS") << " | ";
    cout << endl;
    cout << "----------------------" << endl;
    cout << " ocean pole tide loading coefficients (optl)" << endl;
    _optl_coeff.show();
    cout << "----------------------" << endl;
    cout << " eccentricity values (ecc)" << endl;
    _ecc.show();
    cout << "----------------------" << endl;
    cout << " ocean loading coefficients (blq)" << endl;
    copy(_ol_coeff.begin(), _ol_coeff.end(), ostream_iterator<float>(cout, " "));
    cout << endl;
    cout << "----------------------" << endl;
    cout << " atmospheric pressure loading non tidal (bindisp)" << endl;
    cout << " rows: " << _aplo_xyz.rows() << " / cols: " << _aplo_xyz.cols() <<
         endl;
    cout << "----------------------" << endl;
    cout << " atmospheric pressure loading coeffs tidal (hps)" << endl;
    map < string, ivg::Wave >::iterator iter;
    for (iter = _aplo_tidal.begin(); iter!= _aplo_tidal.end(); iter++)
    {
        cout << scientific << iter->first <<
             " phase frequency acceleration up_cos east_cos north_cos up_sin east_sin north_sin"
             << endl;
        cout << "   " << iter->second.phase << " " <<  iter->second.freq << " " <<
             iter->second.accel;
        cout << " " << iter->second.up_cos << " " <<  iter->second.east_cos << " " <<
             iter->second.north_cos;
        cout << " " <<  iter->second.up_sin << " " <<  iter->second.east_sin << " " <<
             iter->second.north_sin << endl;
    }
    cout << "----------------------" << endl;
    cout << " antenna information" << endl;
    cout << " mounting type: " << _antenna_info.mounting_type <<
         " / Axis Offset Length: " << _antenna_info.axis_offset_length << endl;
    cout << "----------------------" << endl;
    cout << " vmf1 data matrix " << endl;
    cout << " rows: " << _vmf1_data.rows() << " / cols: " << _vmf1_data.cols() <<endl;
    cout << " vmf3 data matrix " << endl;
    cout << " rows: " << _vmf3_data.rows() << " / cols: " << _vmf3_data.cols() <<endl;
    cout << "------------------------------------------------" << endl;
    cout << " hydlo data matrix " << endl;
    cout << " rows: " << _hydlo.rows() << " / cols: " << _hydlo.cols() <<endl;
    cout << "------------------------------------------------" << endl;
    vector<ivg::band> bands = {ivg::band::X, ivg::band::S};
    for(auto &info: _sefd)
    {
        if(info.first == ivg::band::X)
             cerr << "SEFD Info for X-Band: ";
        else if(info.first == ivg::band::S)
            cerr << "SEFD Info for S-Band: ";
        
        show_vector(info.second);
    }
    cout << "------------------------------------------------" << endl;
    for(auto &info: _freq_sequences)
    {
        if(info.first == ivg::band::X)
             cerr << "Frequency Sequence for X-Band";
        else if(info.first == ivg::band::S)
            cerr << "Frequency Sequence for S-Band";
        
        info.second.show();
    }
    cout << "------------------------------------------------" << endl;
    cerr << "Mask Information \n ele | azi" << endl;
    _mask.show(5);
    cout << "------------------------------------------------" << endl;
    

}
// ...........................................................................
void Analysis_station::set_discontinuity( const ivg::Date refepoch, const ivg::Matrix pos0,
        const ivg::Matrix vel0, const vector<ivg::Date> discons)
// ...........................................................................
{
        _refepoch = refepoch;
        _pos0 = pos0;
        _vel0 = vel0;
        _discontinuity = discons;
}
// ...........................................................................
void Analysis_station::add_discontinuity( const ivg::Matrix pos,
        const ivg::Matrix vel,
        const ivg::Date refepoch, const ivg::Date discontinuity,
        ivg::Matrix pos_std, ivg::Matrix vel_std )
// ...........................................................................
{
    if( _refepoch.get_double_mjd() == refepoch.get_double_mjd() )
    {
        _pos0.append_cols(pos);
        _vel0.append_cols(vel);
        _pos0_std.append_cols(pos_std);
        _vel0_std.append_cols(vel_std);        
        _discontinuity.push_back( discontinuity );
    }
    else
    {
        stringstream errormessage;
        errormessage <<
                     "void Analysis_station::add_discontinuity(...): reference epochs do not match";
        throw runtime_error( errormessage.str() );
    }

}

// ...........................................................................
bool Analysis_station::has_discontinuity()
// ...........................................................................
{
    if(	_discontinuity.size() > 1 )
        return true;
    else
        return false;
}
// ...........................................................................
ivg::Matrix Analysis_station::get_xyz( ivg::Date epoch )
// ...........................................................................
{
    ivg::Matrix outPos;
    double mjd = epoch.get_double_mjd();
    double ref = _refepoch.get_double_mjd();
    
    int epo = _discontinuity.size()-1;

    if ( has_discontinuity() )
    {
        for( int i = 0; i < _discontinuity.size() -1 ; i++ )
        {
            if( mjd >= _discontinuity[i].get_double_mjd() && mjd < _discontinuity[i+1].get_double_mjd())
                epo = i;
        }
        outPos = _pos0(":",epo);
    }
    else
        outPos = _pos0;
    
    
    return outPos;
}
// ...........................................................................
ivg::Matrix Analysis_station::get_vel( ivg::Date epoch )
// ...........................................................................
{
    ivg::Matrix outPos;
    double mjd = epoch.get_double_mjd();
    double ref = _refepoch.get_double_mjd();
    
    int epo = _discontinuity.size()-1;

    if ( has_discontinuity() )
    {
        for( int i = 0; i < _discontinuity.size() -1 ; i++ )
        {
            if( mjd >= _discontinuity[i].get_double_mjd() && mjd < _discontinuity[i+1].get_double_mjd())
                epo = i;
        }
        outPos = _vel0(":",epo);
    }
    else
        outPos = _vel0;
    
    
    return outPos;
}
// ...........................................................................
ivg::Matrix Analysis_station::calc_xyz( ivg::Date epoch,
                                        std::vector<string> disp, void * ephem, ivg::Eop_series * eops)
// ...........................................................................
{
    ivg::Matrix pos0_vel0 = calc_xyz( epoch );
    ivg::Matrix d_pos = calc_dxyz(epoch, disp, ephem, eops);

    return pos0_vel0 + d_pos;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_xyz( ivg::Date epoch,
                                        std::vector<string> disp, void * ephem, ivg::Eop_series * eops, ivg::Matrix trf2crf)
// ...........................................................................
{
    ivg::Matrix pos0_vel0 = calc_xyz( epoch );
    ivg::Matrix d_pos = calc_dxyz(epoch, disp, ephem, eops,trf2crf);

    return pos0_vel0 + d_pos;
}
// ...........................................................................
ivg::Matrix Analysis_station::calc_xyz( ivg::Date epoch )
// ...........................................................................
{
    ivg::Matrix outPos;
    double mjd = epoch.get_double_mjd();
    double ref = _refepoch.get_double_mjd();
    
    int epo = _discontinuity.size()-1;

    if ( has_discontinuity() )
    {
        if( mjd < _discontinuity[0].get_double_mjd() )
            epo = 0;
        else
        {
            for( int i = 0; i < _discontinuity.size() -1 ; i++ )
            {
                if( mjd >= _discontinuity[i].get_double_mjd() && 
                    mjd < _discontinuity[i+1].get_double_mjd())
                    epo = i;
            }
        }
        
        // calculate position vector (xyz) to a given epoch
        outPos = _pos0(":",epo) + _vel0(":",epo) / 365.25 * ( mjd - ref );
    }
    else
    {
        // calculate position vector (xyz) to a given epoch
        outPos = _pos0 + _vel0 / 365.25 * ( mjd - ref );
    }

    return outPos;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_xyz( ivg::Date epoch, ivg::Matrix & std )
// ...........................................................................
{
    ivg::Matrix outPos;
    double mjd = epoch.get_double_mjd();
    double ref = _refepoch.get_double_mjd();
    
    int epo = _discontinuity.size()-1;

    if ( has_discontinuity() )
    {
        if( mjd < _discontinuity[0].get_double_mjd() )
            epo = 0;
        else
        {
            for( int i = 0; i < _discontinuity.size() -1 ; i++ )
            {
                if( mjd >= _discontinuity[i].get_double_mjd() && 
                    mjd < _discontinuity[i+1].get_double_mjd())
                    epo = i;
            }
        }
        
        // calculate position vector (xyz) to a given epoch
        outPos = _pos0(":",epo) + _vel0(":",epo) / 365.25 * ( mjd - ref );
        
        ivg::Matrix F(3,6,0.0);
        F(0,0) = 1.0;
        F(1,1) = 1.0;
        F(2,2) = 1.0;
        F(0,3) = ( mjd - ref ) / 365.25;
        F(1,4) = ( mjd - ref ) / 365.25;
        F(2,5) = ( mjd - ref ) / 365.25;
        
        ivg::Matrix Cov = _pos0_std(":",epo);
        Cov.append_rows( _vel0_std(":",epo) );
        
        std = (( F * Cov.diag() * F.transpose() ).sqrt());       
    }
    else
    {
        // calculate position vector (xyz) to a given epoch
        outPos = _pos0 + _vel0 / 365.25 * ( mjd - ref );
        
        ivg::Matrix F(3,6,0.0);
        F(0,0) = 1.0;
        F(1,1) = 1.0;
        F(2,2) = 1.0;
        F(0,3) = ( mjd - ref ) / 365.25;
        F(1,4) = ( mjd - ref ) / 365.25;
        F(2,5) = ( mjd - ref ) / 365.25;
        
        ivg::Matrix Cov = _pos0_std;
        Cov.append_rows(_vel0_std);
        
        std = (( F * Cov.diag() * F.transpose() ).sqrt());
    }
    
    std.from_diag_matrix();

    return outPos;
}


// ...........................................................................
ivg::Matrix Analysis_station::calc_dxyz( ivg::Date epoch,
        std::vector<string> disp, void * ephem, ivg::Eop_series * eops)
// ...........................................................................
{
#ifdef DEBUG_REFFRAME
   cerr << "+++ ivg::Matrix Analysis_station::calc_dxyz(...)" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
    ivg::Matrix d_tmp( 3,1,0.0 );
    ivg::Matrix d_pos( 3,1,0.0 );

    for( int i=0; i<disp.size(); ++i)
    {
        if( disp.at( i ) == "ECCENTRICITY" )
        {
            d_tmp = get_eccentricity( epoch );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "SOLID EARTH TIDES" )
        {
            //Only calculate new c2t if the epoch is different (performance boost!)
            ivg::Matrix c2t;
            if(epoch == _last_crf2trf.begin()->first ){
                c2t = _last_crf2trf.begin()->second;
            }else{
                c2t = eops->form_crf2trf( epoch );
                _last_crf2trf.erase(_last_crf2trf.begin());
                _last_crf2trf[epoch] = c2t;
            }
            
//            ivg::Matrix c2t = eops->form_crf2trf( epoch );
            d_tmp = calc_solid_earth_tides( epoch, ephem, c2t );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "POLE TIDE" )
        {
            d_tmp = calc_pole_tide( epoch, eops );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "OCEAN LOADING" )
        {
            d_tmp = calc_ocean_loading( epoch );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "OCEAN POLE TIDE LOADING" )
        {
            d_tmp = calc_optl( epoch, eops );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "NON TIDAL APLO" )
        {
            d_tmp = calc_nontidal_aplo( epoch, "cspline" );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "TIDAL APLO" )
        {
            d_tmp = calc_tidal_aplo( epoch );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "HYDROLOGY LOADING" )
        {
            d_tmp = calc_hydlo( epoch, "cspline"  );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "PSD" )
        {
            d_tmp = calc_psd_displacement( epoch );
            d_pos += d_tmp;
        }
	else if( disp.at( i ) == "SEASONALS" )
        {
	    d_tmp = calc_seas_displacement( epoch );
	    d_pos += d_tmp;
        }
        else
        {
            stringstream errormessage;
            errormessage << "ivg::Matrix Analysis_station::calc_xyz("
                         << "ivg::Date epoch, std::vector<string> disp, "
                         << "void * ephem, ivg::Eop_series * eops): "
                         << " displacement type" << disp.at( i )
                         <<  " not supported. Exiting";
            throw logic_error( errormessage.str() );
        }
    }
    
#ifdef DEBUG_REFFRAME
   cerr << "--- ivg::Matrix Analysis_station::calc_dxyz(...)" << " : " << tim.toc() << " s " << endl;
#endif 
    return d_pos;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_dxyz( ivg::Date epoch,
					 std::vector<string> disp, void * ephem, ivg::Eop_series * eops, ivg::Matrix trf2crf)
// ...........................................................................
{
#ifdef DEBUG_REFFRAME
   cerr << "+++ ivg::Matrix Analysis_station::calc_dxyz(...)" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
    ivg::Matrix d_tmp( 3,1,0.0 );
    ivg::Matrix d_pos( 3,1,0.0 );

    for( int i=0; i<disp.size(); ++i)
    {
        if( disp.at( i ) == "ECCENTRICITY" )
        {
            d_tmp = get_eccentricity( epoch );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "SOLID EARTH TIDES" )
        {
            //Only calculate new c2t if the epoch is different (performance boost!)
                        
            ivg::Matrix c2t = eops->form_crf2trf( epoch );
	    //	     d_tmp = calc_solid_earth_tides( epoch, ephem, trf2crf.transpose() );
	    d_tmp = calc_solid_earth_tides( epoch, ephem,c2t);
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "POLE TIDE" )
        {
            d_tmp = calc_pole_tide( epoch, eops );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "OCEAN LOADING" )
        {
            d_tmp = calc_ocean_loading( epoch );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "OCEAN POLE TIDE LOADING" )
        {
            d_tmp = calc_optl( epoch, eops );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "NON TIDAL APLO" )
        {
            d_tmp = calc_nontidal_aplo( epoch, "cspline" );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "TIDAL APLO" )
        {
            d_tmp = calc_tidal_aplo( epoch );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "HYDROLOGY LOADING" )
        {
            d_tmp = calc_hydlo( epoch, "cspline"  );
            d_pos += d_tmp;
        }
        else if( disp.at( i ) == "PSD" )
        {
            d_tmp = calc_psd_displacement( epoch );
            d_pos += d_tmp;
        }
	else if( disp.at( i ) == "SEASONALS" )
        {
	    d_tmp = calc_seas_displacement( epoch );
	    d_pos += d_tmp;
        }
        else
        {
            stringstream errormessage;
            errormessage << "ivg::Matrix Analysis_station::calc_xyz("
                         << "ivg::Date epoch, std::vector<string> disp, "
                         << "void * ephem, ivg::Eop_series * eops): "
                         << " displacement type" << disp.at( i )
                         <<  " not supported. Exiting";
            throw logic_error( errormessage.str() );
        }
    }
    
#ifdef DEBUG_REFFRAME
   cerr << "--- ivg::Matrix Analysis_station::calc_dxyz(...)" << " : " << tim.toc() << " s " << endl;
#endif 
    return d_pos;
}
// ...........................................................................
void Analysis_station::move_station( ivg::Matrix d_xyz )
// ...........................................................................
{
    
    if(d_xyz.rows() == 3 && d_xyz.cols() == 1)
    {
        ivg::Matrix delta;
        delta.repmat(d_xyz,1,_pos0.cols());

        _pos0 += delta;
    }
    else
        throw runtime_error("void Analysis_station::move_station( ivg::Matrix d_xyz ): Wrong matrix dimensions. [3x1] expected." );
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_dren( ivg::Date epoch,
        std::vector<string> disp, void * ephem, ivg::Eop_series * eops)
// ...........................................................................
{
    return _xyz2ren(calc_dxyz( epoch, disp, ephem, eops ));
}

// ..............................................................................................
double Analysis_station::calc_baseline_length( ivg::Analysis_station other,
        ivg::Date epoch, double & std )
// ..............................................................................................
{
    // calculate position matrices of both stations for the corresponding epoch
//    ivg::Matrix pos_std, other_std;
//    ivg::Matrix pos = calc_xyz( epoch, pos_std );
//    ivg::Matrix otherPos = other.calc_xyz( epoch, other_std );

    ivg::Matrix pos = calc_xyz( epoch );
    ivg::Matrix otherPos = other.calc_xyz( epoch );    
    
    // calculate baseline length
    ivg::Matrix d = pos - otherPos ;
    double bl = sqrt( pow( d(0),2 ) + pow( d(1),2 ) + pow( d(2),2 ) );

    ivg::Matrix F(1,6);
    F(0,0) = d(0) / bl;
    F(0,1) = -d(0) / bl;
    F(0,2) = d(1) / bl; 
    F(0,3) = -d(1) / bl;
    F(0,4) = d(2) / bl;
    F(0,5) = -d(2) / bl;    
    
    ivg::Matrix Cov = _pos0_std;
    Cov.append_rows(other._pos0_std);
    
//    ivg::Matrix Cov = pos_std;
//    Cov.append_rows( other_std );    
    
    std = (( F * Cov.diag() * F.transpose() ).sqrt())(0);
        
    return bl;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_lat_lon_h( const string ref_ell )
// ...........................................................................
{ 
  _llh0=_calc_lat_lon_h( _pos0, ref_ell );
  return _llh0;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_lat_lon_h( const ivg::Date epoch,
        const string ref_ell )
// ...........................................................................
{
    return _calc_lat_lon_h( calc_xyz( epoch ), ref_ell );
}

// ...........................................................................
ivg::Matrix Analysis_station::_calc_lat_lon_h( ivg::Matrix xyz,
        const std::string ref_ell )
// ...........................................................................
{
    
    
    double a, f;
    if( ref_ell == "GRS80" )
    {
        a = 6378137.0;
        f = 0.00335281068118;
    }
    else if( ref_ell == "TIDE FREE" )
    {
        a = iers::a_E;
        f = iers::f_E;
    }
    else
    {
        stringstream errormessage;
        errormessage << "ivg::Matrix Analysis_station::_calc_lat_lon_h( "
                     << "ivg::Matrix xyz, string ref_ell )"
                     << " invalid type of reference ellipsoid" << ref_ell
                     <<  ". Exiting!";
        throw logic_error( errormessage.str() );
    }

    ivg::Matrix lat_lon_h( 3,1,0.0 );
        
    iers::gconv2_( &a, &f, &xyz(0), &xyz(1), &xyz(2), &lat_lon_h(0),
                   &lat_lon_h(1), &lat_lon_h(2) );
    
    return lat_lon_h;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_az_el( ivg::Date epoch, ivg::Source src )
// ...........................................................................
{
    // get declination and greenwich hour angle
    double dec = src.get_dec0();
    double gha = src.calc_greenwich_hour_angle( epoch );

    // get latitude, longitude and height
    ivg::Matrix llh = calc_lat_lon_h( epoch );

    // calculate local hour angle
    double local_hour_angle = gha + llh(1);
    if( local_hour_angle < 0.0 )
        local_hour_angle += (2.0*M_PI);

    local_hour_angle = fmod( local_hour_angle, (2.0*M_PI) );

    double sine = sin( llh(0) ) * sin( dec )
                  + cos( llh(0) ) * cos( dec ) * cos( local_hour_angle );

    double elevation = asin( sine );
    double xb = sin( local_hour_angle );
    double yb = sin( llh(0) ) * cos( local_hour_angle )
                - cos( llh(0) ) * tan( dec );

    double az = atan2( xb, yb );
    double azimuth = az + M_PI;

    ivg::Matrix AzEl( 1, 2 );
    AzEl(0) = azimuth;
    AzEl(1) = elevation;

    return AzEl;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_az_el( ivg::Date epoch,
        ivg::Matrix src_vec,
        ivg::Matrix crf2trf )
// ...........................................................................
{
    ivg::Matrix azel( 2,1,0.0 );

    // transform source unit vector from CRS 2 TRS
    ivg::Matrix k = crf2trf*src_vec;
    k = k / (k.norm())(0);

    return calc_az_el( k );
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_az_el( const  ivg::Matrix& src_vec_trf )
// ...........................................................................
{
    ivg::Matrix azel( 2,1,0.0 );

    // rotate to topocentric coordinates (REN)
    ivg::Matrix k = _xyz2ren( src_vec_trf );

    // calculate azimuth
    azel(0) = atan2( k(1),k(2) );
    if( azel(0) < 0.0 )
        azel(0) += 2*M_PI;

    // calculate elevation
    azel(1) = asin( k(0) );  //zd = acos(lq(3));
    return azel;
}

// ...........................................................................
ivg::Matrix Analysis_station::azel2k( double az, double el, ivg::Matrix xyz )
// ...........................................................................
{
    

    // (1) calculate local source vector using elevation and azimuth
    // in North East Up component
   ivg::Matrix k_local(3, 1, 0.0);
   k_local(0) = cos(el) * cos(az);
   k_local(1) = cos(el) * sin(az);
   k_local(2) = sin(el);
   
//   // first rotation around y axis in local system (align z axis of local and terrestial system) 
//   // corresponds to latitude
//   double beta = acos( xyz(2) / (xyz.norm()(0)) );
//
//   // second rotation around z axis in local system (align all axis) 
//   // corresponds to longitude
//   double delta = atan2( xyz(1), xyz(0) );
//   
////   ivg::Matrix L2R( 3, 3, 0.0 );
////   L2R( 0,0 ) = 1.0;
////   L2R( 1,1 ) = -1.0;
////   L2R( 2,2 ) = 1.0;  
//   
//   
//    // set rotation matrices: R1 = rotation about y-axis; R2 = rotation about z-axis
//   ivg::Matrix R1( 3,3,0.0 );
//   R1( 0,0 ) = cos( beta );  R1( 0,1 ) = 0.0;  R1( 0,2 ) = sin( beta );  
//   R1( 1,0 ) = 0.0; 	     R1( 1,1 ) = 1.0;  R1( 1,2 ) = 0.0; 	 
//   R1( 2,0 ) = -sin( beta ); R1( 2,1 ) = 0.0;  R1( 2,2 ) = cos( beta );  
//
//
//   ivg::Matrix R2( 3,3,0.0 );
//   R2( 0,0 ) = cos( delta );  R2( 0,1 ) = -sin( delta ); R2( 0,2 ) = 0.0;  
//   R2( 1,0 ) = sin( delta );  R2( 1,1 ) = cos( delta );  R2( 1,2 ) = 0.0;  
//   R2( 2,0 ) = 0.0;  	      R2( 2,1 ) = 0.0;  	 R2( 2,2 ) = 1.0;  
//   
//   // transform local source vector to global (cartesian) system
//   ivg::Matrix M = R2.transpose() * R1.transpose();
//   
//   ivg::Matrix k = M*k_local;
//   
//   //TODO check correctness
//   k(0) = -k(0);
//   
//   std::cout << "beta: " << beta*ivg::rad2d << "   delta:" << delta*ivg::rad2d << std::endl;
   
   
   return  this->_ren2xyz( k_local, "ner");
   
   //return k;
}

// ...........................................................................
ivg::Matrix Analysis_station::azel2k( double az, double el, ivg::Date epoch )
// ...........................................................................
{
   ivg::Matrix xyz = calc_xyz( epoch );
   return azel2k( az, el, xyz );
}


// ...........................................................................
bool Analysis_station::check_visibility( double azi, double ele )
// ...........................................................................
{
    bool answer = false;
    // check if horizon allows aiming to position [azi,ele]
    if( ele >= get_ele_mask( azi ) )
    {
        // if point is generaly visibile at the horizon, check if telescope is able to aim due to limitations
        if( ele >= _antenna_info.ele_min && ele <= _antenna_info.ele_max )
            answer = true;
    }
    
    return answer;
}

double Analysis_station::calc_gravdef( double el)
// .........................................................................................
{
  if (_grav_deform.rows()==0)
    return 0;
  el=el*180/M_PI;
  for (int i=0;i<_grav_deform.rows();i++)
    {
      if (_grav_deform(i,0)>=el)
	{
	  if (i==0)
	    {
	      return _grav_deform(i,1);
	    }
	  else
	    {
	      double val;
	      val=(_grav_deform(i,1)*(el-_grav_deform(i-1,0))-_grav_deform(i-1,1)*(el-_grav_deform(i,0)))/(_grav_deform(i,0)-_grav_deform(i-1,0));
	      return val;
	    }
	  
	}
    }
  
  log<WARNING>("!!! Could not calulate gravitational deformations for station ") % _names[staname::ivs_name];
  return 0;
    
}

// ...........................................................................
double Analysis_station::get_ele_mask( double azi )
// ...........................................................................
{
    // if no elevation mask is set, return as always visibile (ele = 0.0)
    if(_mask.rows() == 0)
        return 0.0;
    else
    {
        // interpolate mask at a given azimut
        std::vector<double> azis = _mask.get_col(0).get_vec();
        
        
        std::vector<double>::iterator it = std::upper_bound( azis.begin(), azis.end(), azi );
        int idx  = std::distance(azis.begin(), it)-1;
        
        
        return _mask(idx,1);
    }
}
// ...........................................................................
bool Analysis_station::calc_wrap_and_rotangle( double old_wrap, double new_azi, double &new_wrap, double &d_azi)
// ...........................................................................
{
    // tolerance for wrap-edge and 180° slew
    // -> returning true/false for critical or not cirtical turn
    double tolerance = 2.5*ivg::d2rad; // 2.5° degree tolerance
    
    // azimut angle to turn from old aiming position to new one
    d_azi = 1000.0*(M_PI/180); // big fake angle
    new_wrap = 0.0;
    bool answer = true;
    for(int i=0; i<3; i++)
    {
        // calculates possible wraps for a given azimut
        // e.g. az = 100° -> wraps using i from 0 to 2 => [100°,460°,820°]
        double tmp_wrap = new_azi + i*(2*M_PI);
        // check if wrap is between telescope-limits => if yes, its a candidate
        if(tmp_wrap >= (_antenna_info.azi_min+tolerance) && tmp_wrap <= (_antenna_info.azi_max-tolerance))
        {
            // find out which wrap is the one with the shortest rotation angle to turn
            double diff = abs(tmp_wrap - old_wrap);
            if( diff < d_azi )
            {
                d_azi = diff;
                new_wrap = tmp_wrap;
            }
        }
    }
             
    // if turn-angle is too close to 180°, see tolerance
    if(abs(d_azi-M_PI) < tolerance )
        answer = false;
    
    if(abs(new_wrap-_antenna_info.azi_min) < tolerance || abs(new_wrap-_antenna_info.azi_max) < tolerance)
        answer = false;
    
    return answer;
}   
// ...........................................................................
string Analysis_station::determine_wrap_zone(double wrap)
// ...........................................................................
{
    // detect wrap zone (W,-,C)
    double overlap = ( _antenna_info.azi_max -_antenna_info.azi_min ) - (2*M_PI);
    
    string symbol = "-";
    if(wrap <= (_antenna_info.azi_min + overlap))
        symbol = "W";
    else if( wrap >= (_antenna_info.azi_max - overlap))
        symbol = "C";
    
    return symbol;
}
// ...........................................................................
ivg::Matrix Analysis_station::calc_pole_tide( const ivg::Date epoch,
        ivg::Eop_series * eops )
// ...........................................................................
{
    
    
    ivg::Eop_series eop_series;
    if(eops == NULL)
    {
        eop_series = ivg::Eop_series( "/data/bakkari/vievs_C04_08_1962_now.txt",
                                      "C04" );;
        eops = &eop_series;
    }

    // calculate mean pole
    ivg::Matrix mean_pole = eops->calc_mean_pole( epoch );
    
    // calculate ERPs
    ivg::Matrix eop = eops->calc_erp( epoch );
     
     //ivg::Matrix d_erp = eops->calc_subdaily_erp_iers( epoch );
    
     //eop += d_erp;

    // calculate wobble parameters
    ivg::Matrix w = eops->calc_wobble_params( epoch, mean_pole, eop );
     
    w *= ivg::rad2as;

    // calculate latitude, longitude and height
    ivg::Matrix llh = calc_lat_lon_h( epoch );
   
    // calculate displacements due to pole tide loading
    ivg::Matrix disp_clr( 3,1,0.0 );
    double colat = M_PI/2.0-llh(0);
    double lambda = llh(1);

    disp_clr(2) = -33.0*sin( 2.0*colat )* ( w(0)*cos( lambda ) + w(1)*sin(
            lambda ) );  // radial
    disp_clr(0) =  -9.0*cos( 2.0*colat )* ( w(0)*cos( lambda ) + w(1)*sin(
            lambda ) );  // co-latitude = -north
    disp_clr(1) =   9.0*cos( colat )    * ( w(0)*sin( lambda ) - w(1)*cos(
            lambda ) );  // longitude = east
    disp_clr *= 1e-3;
    
    // rotate displacements to cartesian coordinates
    ivg::Matrix R( 3,3,0.0 );
    R(0) = cos(colat)*cos(lambda);
    R(1) = -sin(lambda);
    R(2) = sin(colat)*cos(lambda);
    R(3) = cos(colat)*sin(lambda);
    R(4) = cos(lambda);
    R(5) = sin(colat)*sin(lambda);
    R(6) = -sin(colat);
    R(7) = 0.0;
    R(8) = cos(colat);

    ivg::Matrix disp_xyz = R.transpose()*disp_clr;
   
    return disp_xyz ;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_optl( const ivg::Date epoch,
        ivg::Matrix op_coeff_une, ivg::Eop_series * eops )
// ...........................................................................
{
    set_optl_coeff( op_coeff_une );

    ivg::Matrix disp = calc_optl( epoch, eops );

    return disp;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_optl( const ivg::Date epoch,
        ivg::Eop_series * eops )
// ...........................................................................
{
    if (_data_status["OPTL"] != "X")
    {
        log<WARNING>("!!! Ocean Pole Tide Loading of site ") % get_name(ivg::staname::ivs_name) % " not available. Using [0.0 / 0.0 / 0.0]";
        return ivg::Matrix(3,1,0.0);
    }
    
    ivg::Eop_series eop_series;
    if(eops == NULL)
    {
        eop_series = ivg::Eop_series( "/data/bakkari/vievs_C04_08_1962_now.txt",
                                      "C04" );;
        eops = &eop_series;
    }

    // calculate mean pole
    ivg::Matrix mean_pole = eops->calc_mean_pole( epoch );

    // calculate ERPs
    ivg::Matrix eop = eops->calc_erp( epoch );

    // calculate wobble parameters
    ivg::Matrix w = eops->calc_wobble_params( epoch, mean_pole, eop );
    w *= ivg::rad2as;

    // calculate latitude, longitude and height
    ivg::Matrix llh = calc_lat_lon_h( epoch );


    // calculate displacements due to ocean pole tide loading
    ivg::Matrix fac( 1,2,0.0 );
    fac(0) = w(0)*iers::optl_gamma2_re + w(1)*iers::optl_gamma2_im;
    fac(1) = w(1)*iers::optl_gamma2_re - w(0)*iers::optl_gamma2_im;

    fac *= iers::optl_k() * (M_PI/3600.0/180.0);

    ivg::Matrix disp_rne = ( _optl_coeff( ":",0 )*fac(0) + _optl_coeff( ":",
                             1 )*fac(1) );

    // change east and north component, then transform UNE->XYZ
    ivg::Matrix disp_xyz = _ren2xyz( disp_rne , "rne" );

    return disp_xyz;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_solid_earth_tides( ivg::Date epoch,
        void *ephem, ivg::Matrix c2t )
// ...........................................................................
{
    // check void pointer
    if( !ephem )
        cerr << "Ephemeris file not loaded." << endl;

    double r[6];

    double au = jpl_get_double( ephem, JPL_EPHEM_AU_IN_KM )*1e3;

    // get position of moon
    ivg::Matrix x_moon( 3,1 );
    int err_code = -1;
    err_code = jpl_pleph( ephem, epoch.get_jd_tdb(), 10, 3, r, 1 );

    for( int i=0; i<3; ++i )
        x_moon( i ) = r[i]*au;
    x_moon = c2t*x_moon;

    // get position of sun
    ivg::Matrix x_sun( 3,1 );
    err_code = jpl_pleph( ephem, epoch.get_jd_tdb(), 11, 3, r, 1 );

    for( int i=0; i<3; ++i )
        x_sun( i ) = r[i]*au;
    x_sun = c2t*x_sun;

    // calculate displacements due to solid earth tides
    ivg::Matrix disp( 3,1,0.0 );
    int year = epoch.get_int_year();
    int month = epoch.get_int_month();
    int day = epoch.get_int_day();
    double hr = epoch.get_frac_day()*24.0;

    ivg::Matrix pos = calc_xyz( epoch );
    iers::dehanttideinel_( pos.data_ptr(), &year, &month, &day, &hr,
                           x_sun.data_ptr(),
                           x_moon.data_ptr(), disp.data_ptr() );

    return disp;
}

ivg::Matrix Analysis_station::compute_az_el_sun( ivg::Date epoch,  void* ephem, ivg::Matrix c2t ){
    
    // check void pointer
    if( !ephem )
        cerr << "Ephemeris file not loaded." << endl;
    
     // get position of sun
    ivg::Matrix x_sun( 3,1 );
    
    double r[6];
    double au = jpl_get_double( ephem, JPL_EPHEM_AU_IN_KM )*1e3;
    
    int err_code = -1;
    err_code = jpl_pleph( ephem, epoch.get_jd_tdb(), 11, 3, r, 1 );

    for( int i=0; i<3; ++i )
        x_sun( i ) = r[i]*au;
    x_sun = c2t*x_sun;
    
    return this->calc_az_el(epoch, x_sun, c2t);
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_hydlo( ivg::Date epoch,
        const string interpolation_type )
// ...........................................................................
{
    if(_hydlo.rows() == 0)
    {
        stringstream errormessage;
        errormessage << "Analysis_station::calc_hydlo( ivg::Date,const string ): _hydlo not set! Exiting!" << endl;
        throw runtime_error( errormessage.str() );
    }
    
    ivg::Matrix d_ren = _hydlo.interpolate( _hydlo(":",0) ,
                                            epoch.get_double_mjd(), interpolation_type ).get_sub(0,1,0,3).transpose();
    return _ren2xyz( d_ren );;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_ocean_loading( ivg::Date epoch )
// ...........................................................................
{
    if (_data_status["OLC"] != "X")
    {
        log<WARNING>("!!! Ocean Loading Coefficients of site ") % get_name(ivg::staname::ivs_name) % " not available. Using [0.0 / 0.0 / 0.0]";
        return ivg::Matrix(3,1,0.0);
    }
    ivg::Matrix disp_xyz;
      if (_ol_harpos.rows()>0)  //Harpos format
     {
    	double dt=(epoch.get_mjd_tt()-51544.5)*86400;
    	vector<double> disp( 3,0.0 );
    	for (int i=0;i<_ol_harpos.rows();i++)
    	  {
    	    double arg=_ol_harpos(i,0)+_ol_harpos(i,1)*dt+0.5*_ol_harpos(i,2)*dt*dt;
	    
    	    disp[0]+=_ol_harpos(i,3)*cos(arg)+_ol_harpos(i,6)*sin(arg);
    	    disp[1]+=_ol_harpos(i,4)*cos(arg)+_ol_harpos(i,7)*sin(arg);
    	    disp[2]+=_ol_harpos(i,5)*cos(arg)+_ol_harpos(i,8)*sin(arg);
        
    	  }
    	ivg::Matrix d_ren( disp );
	//	std::cout << disp[0]*1000 << " " <<disp[1]*1000 << " " <<disp[2]*1000 <<endl;
    	disp_xyz = _ren2xyz( d_ren );
	
     }
    else
    {
    // divide ocean loading coefficients to amplitude and phase vectors
    vector<float> amp;
    amp.resize(33) ;
    vector<float> pha;
    pha.resize(33);

    copy( _ol_coeff.begin(), _ol_coeff.begin()+32, amp.begin());
    copy( _ol_coeff.begin()+33, _ol_coeff.end(), pha.begin());

    // calculate displacements due to ocean loading
    vector<double> disp( 3,0.0 );

    vector<int> epo = { epoch.get_int_year(), epoch.get_int_doy(), epoch.get_int_hour(), epoch.get_int_min(), (int) epoch.get_double_sec() };

    iers::hardisp_( &epo[0], &amp[0], &pha[0], &disp[0] );

    // transform std::vector to ivg::Matrix
    ivg::Matrix d_ren( disp );

    disp_xyz = _ren2xyz( d_ren );
    //std::cout << disp[0]*1000 << " " <<disp[1]*1000 << " " <<disp[2]*1000 <<endl;
    }
    return disp_xyz;
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_nontidal_aplo( ivg::Date epoch,
        const string interpolation_type)
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ ivg::Matrix Analysis_station::calc_nontidal_aplo( ivg::Date ,const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ivg::Matrix tmp = _aplo_xyz.transpose();
    ivg::Matrix out;
    if(tmp.size(1) == 0 && tmp.size(2) == 0)
    {
//        log<WARNING>("!!! Non Tidal Atmospheric Pressure Loading data of site " ) % get_name(ivg::staname::ivs_name) % " not available";
        out = ivg::Matrix(1,4,0.0);
    }
    else
        out = tmp.interpolate( tmp(":",0) , epoch.get_double_mjd(), interpolation_type );
       
#if DEBUG_REFFRAME >=3
   cerr << "--- ivg::Matrix Analysis_station::calc_nontidal_aplo( ivg::Date ,const string )" << " : " << tim.toc() << " s " << endl;
#endif
   return out.get_sub(0,1,0,3).transpose();
}

// ...........................................................................
ivg::Matrix Analysis_station::calc_tidal_aplo( ivg::Date epoch )
// ...........................................................................
{
#if DEBUG_REFFRAME >=3
   cerr << "+++ ivg::Matrix Analysis_station::calc_tidal_aplo( ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ivg::Matrix disp_ren(3,1,0.0);
    
    if(_aplo_tidal.begin() == _aplo_tidal.end())
    {
//        log<WARNING>("!!! Tidal Atmospheric Pressure Loading data of site ") % get_name(ivg::staname::ivs_name) % " not available";
    }
    else
    {

        ivg::Date t0(2000,1,1,12);
        std::map < string, ivg::Wave >::iterator iter;
        for (iter = _aplo_tidal.begin(); iter!= _aplo_tidal.end(); iter++)
        {
            ivg::Wave S = iter->second;

            double dt = (epoch.get_jd_tt() - t0.get_jd_tt() )*86400.0;
            double term = S.phase + S.freq * dt + (1/2)*S.accel * pow(dt,2);

            double cos_term = cos(term);
            double sin_term = sin(term);

            double d_radi = S.up_cos * cos_term + S.up_sin * sin_term ;
            double d_east = S.east_cos * cos_term + S.east_sin * sin_term ;
            double d_north = S.north_cos * cos_term + S.north_sin * sin_term ;

            disp_ren(0,0) = disp_ren(0,0) + d_radi;
            disp_ren(1,0) = disp_ren(1,0) + d_east;
            disp_ren(2,0) = disp_ren(2,0) + d_north;
        }
        
    }

//    cout << "disp_U: " << disp_ren(0,0) << " / disp_E: " << disp_ren(1,0) << " / disp_N: " << disp_ren(0,0) << endl;
#if DEBUG_REFFRAME >=3
   cerr << "--- ivg::Matrix Analysis_station::calc_tidal_aplo( ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    return _ren2xyz(disp_ren);
}

// ...........................................................................
ivg::Matrix Analysis_station::_ren2xyz( ivg::Matrix ren, string type)
// ...........................................................................
{
    //Expecting columnvector 3x1
    if( ren.rows() != 3 || ren.cols() != 1)
    {
        stringstream errormessage;
        errormessage <<
                     "ivg::Matrix Analysis_station::_ren2xyz( ivg::Matrix ren, string type): wrong vector dimension. should be 3x1";
        throw runtime_error( errormessage.str() );
    }


    if(type == "rne")
    {
        double tmp = ren(2);
        ren(2) = ren(1);
        ren(1) = tmp;
    }
    else if(type == "ner")
    {
        double tmp = ren(2);
        ren(2) = ren(0);
        ren(0) = tmp;
    }
    
    ivg::Matrix rot = form_topo2geo();

    return rot*ren;
}

// ...........................................................................
ivg::Matrix Analysis_station::_xyz2ren( ivg::Matrix xyz )
// ...........................................................................
{
    ivg::Matrix rot = form_topo2geo();
    return rot.transpose()*xyz;
}

// ...........................................................................
ivg::Matrix Analysis_station::form_topo2geo()
// ...........................................................................
{
    ivg::Matrix lat_lon_h = calc_lat_lon_h();

    ivg::Matrix rot(3,3,0.0);

    double cp = cos(lat_lon_h(0));
    double sp = sin(lat_lon_h(0));
    double cl = cos(lat_lon_h(1));
    double sl = sin(lat_lon_h(1));

    rot(0) = cp*cl;
    rot(1) = cp*sl;
    rot(2) = sp;
    rot(3) = -sl;
    rot(4) = cl;
    rot(5) = 0.0;
    rot(6) = -sp*cl;
    rot(7) = -sp*sl;
    rot(8) = cp;

    return rot;
}

// ...........................................................................
ivg::Matrix Analysis_station::interpolate_ext_met_data( std::string type, 
                                                        ivg::Date epoch, 
                                                        std::string interpolation_type )
// ...........................................................................
{
    ivg::Matrix out;
    
    if( type == "mapping" )
    {   
        if(_vmf1_data.rows() <= 10 )
        {
            stringstream ss;
            ss << "ivg::Matrix Analysis_station::interpolate_ext_met_data( ... ): Not enough _vmf1_data for station ";
            ss << _names[staname::ivs_name] << " set. Only " << _vmf1_data.rows() << " datarows in _vmf1_data.";
            throw runtime_error( ss.str() );
        }
        else if( _vmf1_data.get_col(0).max() <= epoch.get_double_mjd() )
            throw runtime_error( "ivg::Matrix Analysis_station::interpolate_ext_met_data( std::string, ivg::Date, std::string ): Requested epoch out of range from _vmf1_data for station "+_names[staname::ivs_name]);
        else
           out = _vmf1_data.interpolate( _vmf1_data(":",0) , epoch.get_double_mjd(), interpolation_type );
    }
    if( type == "mapping3" )
    {   
        if(_vmf3_data.rows() <= 10 )
        {
            stringstream ss;
            ss << "ivg::Matrix Analysis_station::interpolate_ext_met_data( ... ): Not enough _vmf3_data for station ";
            ss << _names[staname::ivs_name] << " set. Only " << _vmf3_data.rows() << " datarows in _vmf3_data.";
            throw runtime_error( ss.str() );
        }
        else if( _vmf3_data.get_col(0).max() <= epoch.get_double_mjd() )
            throw runtime_error( "ivg::Matrix Analysis_station::interpolate_ext_met_data( std::string, ivg::Date, std::string ): Requested epoch out of range from _vmf3_data for station "+_names[staname::ivs_name]);
        else
           out = _vmf3_data.interpolate( _vmf3_data(":",0) , epoch.get_double_mjd(), interpolation_type );
    }
    else if( type == "hydrostatic" )
    {           
        if(_zhd_data.rows() < 2 )
            throw runtime_error( "ivg::Matrix Analysis_station::interpolate_ext_met_data( std::string, ivg::Date, std::string ): No _zhd_data for station "+_names[staname::ivs_name]+" set.");
        else if( _zhd_data.get_col(0).max() <= epoch.get_double_mjd() )
            throw runtime_error( "ivg::Matrix Analysis_station::interpolate_ext_met_data( std::string, ivg::Date, std::string ): Requested epoch out of range from _zhd_data for station "+_names[staname::ivs_name]);
        else
           out = _zhd_data.interpolate( _zhd_data(":",0) , epoch.get_double_mjd(), interpolation_type );
    }
    return out;
}

// .............................................................................
ivg::Matrix Analysis_station::calc_seas_displacement( ivg::Date epo )
// .............................................................................
{
  
  ivg::Matrix disp( 3,1,0.0 );
  for (int i=0;i<_seasonals.rows();i++)
    {
      double mjd = epo.get_double_mjd();
      double period = _seasonals(i,0);
      double omegat = 2*M_PI*(mjd-51544)/period;
      
      disp(0,0)+=_seasonals(i,1)*cos(omegat)+_seasonals(i,3)*sin(omegat);
      disp(1,0)+=_seasonals(i,5)*cos(omegat)+_seasonals(i,7)*sin(omegat);
      disp(2,0)+=_seasonals(i,9)*cos(omegat)+_seasonals(i,11)*sin(omegat);
      
    }
  return disp/1000;
}

// .............................................................................
ivg::Matrix Analysis_station::calc_psd_displacement( ivg::Date epo )
// .............................................................................
{
    ivg::Matrix disp( 3,1,0.0 );
    double dpos;
    if(_psd.begin() != _psd.end())
    {
        double mjd = epo.get_double_mjd();
        
        // find seismic event prior to input date
        int idx = -1;
        double dt = 1e9;
        for(int i=0;i<_psd.size();++i)
        {
            double dt_i = (mjd-_psd.at(i).mjd)/365.25;
            // event prior to input epoch and nearer in time than priviously 
            // detected event
            if(dt_i>0) // && dt_i<dt)
            {
                dt = dt_i;
                idx = i;
		// East
		iers::parametric_(&_psd.at(idx).e_mode,&dt,&_psd.at(idx).e_a1,
				  &_psd.at(idx).e_t1,&_psd.at(idx).e_a2,
				  &_psd.at(idx).e_t2,&dpos);
		disp(1) += dpos;
		// North
		iers::parametric_(&_psd.at(idx).n_mode,&dt,&_psd.at(idx).n_a1,
				  &_psd.at(idx).n_t1,&_psd.at(idx).n_a2,
				  &_psd.at(idx).n_t2,&dpos);
		disp(2) += dpos;
		// Up
		iers::parametric_(&_psd.at(idx).u_mode,&dt,&_psd.at(idx).u_a1,
				  &_psd.at(idx).u_t1,&_psd.at(idx).u_a2,
				  &_psd.at(idx).u_t2,&dpos);
		disp(0) += dpos;

		
            }
        }
        // calculate displacements if there is an event prior to the input epoch
        if(idx != -1)
        {            
            // East
            //iers::parametric_(&_psd.at(idx).e_mode,&dt,&_psd.at(idx).e_a1,
			      // &_psd.at(idx).e_t1,&_psd.at(idx).e_a2,
			      //                  &_psd.at(idx).e_t2,&dpos);
	    // disp(1) = dpos;
            // North
            //iers::parametric_(&_psd.at(idx).n_mode,&dt,&_psd.at(idx).n_a1,
			      // &_psd.at(idx).n_t1,&_psd.at(idx).n_a2,
			      //                  &_psd.at(idx).n_t2,&dpos);
            //disp(2) = dpos;
            // Up
            //iers::parametric_(&_psd.at(idx).u_mode,&dt,&_psd.at(idx).u_a1,
			      //                  &_psd.at(idx).u_t1,&_psd.at(idx).u_a2,
			      //                  &_psd.at(idx).u_t2,&dpos);
	    //    disp(0) = dpos;

            disp *= 1e-3;  // mm -> m
            
            disp = _ren2xyz(disp);
        }
    }
    return disp;
}


// ...........................................................................
ivg::Matrix Analysis_station::calc_pos_diff( ivg::Analysis_station other )
// ...........................................................................
{
    //cerr << "*********** ATOOLS: XYZ0" << endl;
    //get_xyz0().show();
    //other.get_xyz0().show();
        
    ivg::Matrix diff = get_xyz0() - other.get_xyz0();
    //cerr << "******** Analysis_station: DIFF XYZ" << endl;
    //diff.show();
    return _xyz2ren( diff );
}

// ...........................................................................
double ivg::Analysis_station::calc_slant_sefd(ivg::band band, ivg::Source * src_ptr, 
                                              ivg::Date * epo_ptr)
// adjusts the zenith SEFD for the source elevation
// ...........................................................................
{
    double sefdel;

    if( _sefd[band].size()==0 )
        throw runtime_error("double ivg::Analysis_station::calc_slant_sefd(): No _sefd set. No sefd-calculation possible for station "+_names[ivg::staname::ivs_name]);
    
    // use SEFD in zenith direction if no parameters given (first vector entry 
    // is zenith SEFD, the following ones are SEFD parameters)
    if (_sefd[band].size()==1) 
        sefdel = _sefd[band].at(0);                                       
    else // calculate elevation dependent SEFD
    {
        // calc elevation to source at given epoch
        ivg::Matrix azel = calc_az_el(*epo_ptr,*src_ptr);
    
        double xel = 1.0/pow(sin(azel(1)),_sefd[band].at(1));
        double fac = 0.0;
        for (int i = 1; i<_sefd[band].size()-1; ++i)
            fac = fac+_sefd[band].at(i+1)*pow(xel,i-1);
        sefdel = _sefd[band].at(0)*fac;
    }      
   return sefdel;
}

// ...........................................................................
double Analysis_station::calc_bandwidth_rms(ivg::Matrix freq)
// ...........................................................................
{  
   double sum = 0.0;
   double frsqsum = 0.0;
   for( int i=0; i<freq.length(); i++ )
   {
      sum = sum + freq(i);
      frsqsum = frsqsum + freq(i)*freq(i);
   }
   
   double rms = sqrt( (frsqsum / freq.length())
                       - pow((sum / freq.length() ),2.0 ) );

   return rms;
}
// ...........................................................................
double Analysis_station::calc_bandwidth_rms(ivg::band band)
// ...........................................................................
{
   return calc_bandwidth_rms(_freq_sequences[band]);
}

// ...........................................................................
void Analysis_station::get_gpt2_ah_aw(double *ah, double *aw)
// ...........................................................................
{

      *ah=_gpt2_ah;
      *aw=_gpt2_aw;
}

// ...........................................................................
void Analysis_station::set_gpt2_ah_aw(double ah, double aw)
// ...........................................................................
{

      _gpt2_ah=ah;
      _gpt2_aw=aw;
}

// ...........................................................................
void Analysis_station::get_gpt3_ah_aw(double *ah, double *aw)
// ...........................................................................
{

      *ah=_gpt2_ah;
      *aw=_gpt2_aw;
}

// ...........................................................................
void Analysis_station::set_gpt3_ah_aw(double ah, double aw)
// ...........................................................................
{

      _gpt2_ah=ah;
      _gpt2_aw=aw;
}

} // # namespace ivg











