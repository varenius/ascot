#include "turbulence.h"

/* 	======================================================================
			   class 'Turbulence' by S. Halsig 
			   	     CODE-file
	======================================================================
*/

namespace ivg
{

// ==============================================
// =============== Konstruktoren: ===============
// ==============================================

// Default-Konstruktor
// .....................................................................
Turbulence::Turbulence( )
// .....................................................................
{
    
    // structure constant and tropospheric height
    double cn = 1.0 * 1e-7;
    double h = 2000.0;

    // ivg::Matrix for wind velocity (conversion from [m/s] to [m/h])
    ivg::Matrix vel( 1,2 ); 
    vel(0) = 8.00 * 3600;   // [m/h] north component of wind vector
    vel(1) = 8.00 * 3600;   // [m/h] east component of wind vector
    double v_dir = M_PI; 

    double v = 8.0;         // default wind velocity [m/s]
    double zwd0 = 0.0;
    double dh_seg = 1.0;
    double dh = 1.0;

    // set struct with turbulence data
    _turb_data = { cn, h, v, vel, v_dir, zwd0, dh_seg, dh };

    // set global turbulence parameters: 
    // (a) spectral coefficients
    ivg::Matrix spectral_coeffs(3,1,1.0);
    _spectral_coeffs = spectral_coeffs; 

    // (b) set outer scale lenth and saturation length scale
    set_outer_scale_length( ivg::L0 );
    set_saturation_length( ivg::L );
    
}

// .....................................................................
Turbulence::Turbulence( std::map< std::string, turbulence_data > turb_sta, ivg::Matrix spectral_coeffs )
// .....................................................................
{
   
    // set station dependent turbulence parameters
    _turb_sta = turb_sta;

    // set global turbulence parameters: 
    // (a) spectral coefficients
    _spectral_coeffs = spectral_coeffs; 

     // (b) set outer scale lenth and saturation length scale
     set_outer_scale_length( ivg::L0 );
     set_saturation_length( ivg::L );

}

// .....................................................................
Turbulence::Turbulence( std::vector<ivg::Analysis_station> stations, 
                        double Cn, double H, double v, 
                        ivg::Matrix vel, double v_dir )
// .....................................................................
{

    double zwd0 = 0.0;
    double dh_seg = 1.0;
    double dh = 1.0;

    _turb_data = { Cn, H, v, vel, v_dir, zwd0, dh_seg, dh }; 

    std::vector<Analysis_station>::iterator iter;
    for( iter = stations.begin(); iter < stations.end(); iter++ )
    {
       _turb_sta[ iter->get_name( ivg::staname::ivs_name ) ] = _turb_data;   
    }

    // set global turbulence parameters: 
    // (a) spectral coefficients
    ivg::Matrix spectral_coeffs(3,1,1.0);
    _spectral_coeffs = spectral_coeffs; 
       	
    // (b) set outer scale lenth and saturation length scale
    set_outer_scale_length( ivg::L0 );
    set_saturation_length( ivg::L ); 	
}


// ==============================================
// =============== Destruktoren: ================
// ==============================================

Turbulence::~Turbulence( )
{
}


// ======================================================
// ================= Memberfunktionen: ==================
// ======================================================

// .....................................................................
void Turbulence::show( int nk )
// .....................................................................
{
#if DEBUG_TROP >=3
   cerr << "+++ ivg::Matrix Turbulence::show( int )" << endl; 
   tictoc tim;
   tim.tic();
#endif  

   cerr << "--------------- Turbulence.show() --------------- " << endl;
   cerr << "TURBULENCE PARAMETER: \n" << endl;

   for( map<std::string,turbulence_data>::iterator it = _turb_sta.begin(); it != _turb_sta.end(); ++it )
   {
      cerr << "Station " << it->first << ": " << endl;
      cerr << "   Cn: " << it->second.Cn << endl;
      cerr << "   H: " << it->second.h << endl;
      cerr << "   v: " << it->second.v << endl;
      cerr << "   vel: " << endl;
      it->second.vel.show();
      cerr << "  \n v_dir: " << it->second.v_dir << " \n " << endl;
   }
   
   cerr << "--------------- Turbulence.show() --------------- " << endl;

#if DEBUG_TROP >=3
   cerr << "+++ ivg::Matrix Turbulence::show( int )"  
        << " : " << tim.toc() << " s " << endl; 
#endif  
}

// .....................................................................
void Turbulence::set_saturation_length( double L )
// .....................................................................
{
   _L = L;
}

// .....................................................................
void Turbulence::set_outer_scale_length( double L0 )
// .....................................................................
{
   _L0 = L0;
}

// .....................................................................
void Turbulence::set_turb_params( std::map< std::string, turbulence_data > turb_sta )
// .....................................................................
{
   _turb_sta = turb_sta;
}

// .....................................................................
void Turbulence::set_spectral_coefficients( ivg::Matrix coeffs )
// .....................................................................
{
   _spectral_coeffs = coeffs;
}


// simulation of Equivalent Zenith Wet Delay (EZWD) and mapping to Zenith Wet Delay (ZWD)
// .....................................................................
ivg::Matrix Turbulence::simulate_ezwd( ivg::Matrix & VCM, double zwd0 )
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ ivg::Matrix Turbulence::simulate_ezwd( ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif  

    // cholesky decomposition of the variance-covariance matrix and 	
    ivg::Matrix D = VCM ;
    D.chol( D, 'L' ) ;
    
    // create vector of zero mean Gaussian random numbers with variance one	
    ivg::Matrix w( D.size(2), 1, 0.0 );
    w.rand_norm( 0.0, 1.0 ) ;   

    // get initial equivalent zenith wet delays   
    ivg::Matrix iEZWD( D.size(1), 1, zwd0 );    

    // simulate equivalent zenith wet delays
    ivg::Matrix EZWD = iEZWD + D * w ;

#if DEBUG_TROP >=2
   cerr << "--- ivg::Matrix Turbulence::simulate_ezwd( ivg::Matrix & )"
        << " : " << tim.toc() << " s " << endl; 
#endif  

    return ( EZWD ) ;	
}

std::map< std::string, turbulence_data > Turbulence::set_param_form_cfg(const Setting& setup, ivg::Trf& trf){
    
    Setting &params = setup[ "SIM" ];
    Setting &groups = setup[ "groups" ];
    
    vector<std::string> sta_names = trf.get_station_names( ivg::staname::ivs_name );
    std::sort( sta_names.begin(), sta_names.end() );  
    
    turbulence_data turb_params;
    std::map< std::string, turbulence_data > turb_sta;

    for( int j=0;j<params[ "troposphere" ]["turb_params"].getLength();++j )
    {
       Setting &turb_params_block = params[ "troposphere" ]["turb_params"][ j ];
       
       std::vector< std::string > grp_names;
       group2names( turb_params_block, groups, "stations", sta_names, grp_names );
       
       for( int i=0;i<grp_names.size();++i )
       { 
           
          double Cn = (double)turb_params_block[ "Cn" ];      // structure constant 
          double H = (double)turb_params_block[ "H" ];        // effective tropospheric height
          double v = (double)turb_params_block[ "v" ];        // wind velocity
          double v_dir = (double)turb_params_block[ "v_dir" ]; // wind direction in elevation
          
          double v_az = (double)turb_params_block[ "vaz" ]*ivg::d2rad; // azimuthal wind direction in degree
          
          ivg::Matrix vel(2,1,0);
          vel(0) = cos(v_az);
          vel(1) = sin(v_az);
          
          vel*= v;
          
          turb_params = { Cn, H, v, vel, v_dir, 0.0, 1.0, 1.0 }; 
          turb_sta[ grp_names.at( i ) ] = turb_params;
       }
    }    

    _turb_sta = turb_sta;
    return turb_sta;
}


/*
// .....................................................................
ivg::Matrix Turbulence::map_ezwd ( ivg::Matrix & EZWD ) 
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ ivg::Matrix Turbulence::map_ezwd ( ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
	
	ivg::Matrix m ( _data.size(1), 1, 0.0 );
	
	// calculate mapping function
	for ( int i = 0; i < _data.size(1); i++ )
	{
		m( i ) = 1.0 / sin( _data( i,1 ) );
	}
	
	// map equivalent zenith wet delay (ZWD) to slant wet delay (SWD)
	ivg::Matrix ones( _data.size( 1 ), 1, 1.0 );
	ivg::Matrix SWD = ones / m * EZWD ;

#if DEBUG_TROP >=2
   cerr << "--- ivg::Matrix Turbulence::map_ezwd ( ivg::Matrix & )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  
	
	return SWD ;	
}
*/



// ****************************************************************************
// *********** (1) M A T E R N   C O V A R I A N C E   F A M I L Y ************
// ****************************************************************************
// .....................................................................
ivg::Matrix Turbulence::calc_matern_vcm_model( int nobs, std::vector<ivg::Scan> scans, ivg::Trf *trf, 
                                               ivg::Matrix spectral_coeffs, ivg::Matrix & C )
// .....................................................................
{
#if DEBUG_TROP >=1
   cerr << "+++ ivg::Matrix Turbulence::calc_matern_vcm_model( int, std::vector<ivg::Scan>, ivg::Trf*, "  
        << "std::map< std::string, turbulence_data >, ivg::Matrix, ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif  

    ivg::Matrix VCM( nobs,nobs,0.0 );

    double a = spectral_coeffs(0);
    double c = spectral_coeffs(2);

    // wavenumber of the outer scale length [1/m]
    double k0 = 2.0 * M_PI / _L0;

    // calculate elecromagnetic wavenumber
    double lambda = ivg::c / ivg::freq_x ;
    double elec_wn = 2.0 * M_PI / lambda ;

    double tau, Cn2;

    ivg::Obs* obs_ptr1;
    ivg::Obs* obs_ptr2;

    int counter1 = 0;
    int counter2 = 0;
    int J = 0;
    int I = 0;

    ivg::Matrix az_el_i1;
    ivg::Matrix az_el_i2;
    ivg::Matrix az_el_j1;
    ivg::Matrix az_el_j2;

    std::string sta1, sta2, sta3, sta4;

    // calculate a priori station positions
    std::map< std::string, ivg::Matrix > sta;

    std::vector<Analysis_station>::iterator iter;
    for( iter = trf->begin(); iter < trf->end(); iter++ )
    {
        sta[ iter->get_name( ivg::staname::ivs_name ) ] = iter->calc_xyz( trf->get_reference_epoch() );   
    }

    // fill upper right part of the variance covariance matrix
    for( int k = 0; k <= scans.size()-1; k++ )                    // loop over scans k
    {
        for( int i = 0; i <= scans.at(k).get_nobs()-1; i++ )       // loop over obs i
        {
            counter2 = 0;
            I = i+counter1;

            // get pointer to current observation i
            obs_ptr1 = scans.at(k).get_obs_ptr(i);

            // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
            az_el_i1 = obs_ptr1->get_az_el(1);
            az_el_i2 = obs_ptr1->get_az_el(2);
            obs_ptr1->get_station_names( sta1,sta2 );

            for( int l = 0; l <= scans.size()-1; l++ )              // loop over scans l
            {   
                // calculate time differences [hours]
                tau = abs( scans.at(k).get_epoch().get_double_mjd() - scans.at(l).get_epoch().get_double_mjd() ) * 24.0;
                if( tau == 0.0 )
                   tau = 1e-16;

                for( int j = 0; j <= scans.at(l).get_nobs()-1; j++ ) // loop over obs j
                {
                    J = j+counter2;

                    if( VCM(J,I) == 0.0 ) // only fill upper right triangle to reduce computational cost
                    {                   
                        // get pointer to current observation j
                        obs_ptr2 = scans.at(l).get_obs_ptr(j);

                        // get azimuth and elevation as well as station name of station 1 and 2 for observation j  
                        az_el_j1 = obs_ptr2->get_az_el(1);
                        az_el_j2 = obs_ptr2->get_az_el(2);
                        obs_ptr2->get_station_names( sta3,sta4 );

                        // VARIANCES
                        if ( I == J )   
                        {
//                            VCM(I,J) = 0.782e2 * ( pow(elec_wn,2.0) * _turb_data.h * Cn2 * c * pow(k0,-5.0/3.0) ) / ( pow( sin( az_el_i1(1) ), 2.0) )
//                                     + 0.782e2 * ( pow(elec_wn,2.0) * _turb_data.h * Cn2 * c * pow(k0,-5.0/3.0) ) / ( pow( sin( az_el_i2(1) ), 2.0) );

                            VCM(I,J) = 0.782e2 * ( pow(elec_wn,2.0) * _turb_sta[ sta1 ].h * pow(_turb_sta[sta1].Cn,2.0) * c * pow(k0,-5.0/3.0) ) / ( pow( sin( az_el_i1(1) ), 2.0) )
                                     + 0.782e2 * ( pow(elec_wn,2.0) * _turb_sta[ sta2 ].h * pow(_turb_sta[sta2].Cn,2.0) * c * pow(k0,-5.0/3.0) ) / ( pow( sin( az_el_i2(1) ), 2.0) );
                        }
                        // COVARIANCES
                        else   
                        {
                            // (*) check which observations are correlated (in time); 
                            //     spatial correlation is considered by the separation distance (= distance between two rays)  
                            // (1) calcualte separation distance between two rays
                            // (2) calculate covariance between observations i and j 
                            if( sta1 == sta3 && sta2 == sta4 )
                            {
                               double d1 = get_dist_between_rays( az_el_i1(1), az_el_i1(0), sta[ sta1 ], sta[ sta3 ], 
                                                                  _turb_sta[ sta1 ].h, _turb_sta[ sta3 ].h );
                               double d2 = get_dist_between_rays( az_el_i2(1), az_el_i2(0), sta[ sta2 ], sta[ sta4 ], 
                                                                  _turb_sta[ sta2 ].h, _turb_sta[ sta4 ].h );

                               VCM(I,J) = _calc_matern_covariance( tau, az_el_i1, az_el_j1, sta1, sta3, spectral_coeffs, d1 )
                                        + _calc_matern_covariance( tau, az_el_i2, az_el_j2, sta2, sta4, spectral_coeffs, d2 );
                            }
                            else if( sta1 == sta4 && sta2 == sta3 )
                            {
                               double d1 = get_dist_between_rays( az_el_i1(1), az_el_i1(0), sta[ sta1 ], sta[ sta4 ], 
                                                                  _turb_sta[ sta1 ].h, _turb_sta[ sta4 ].h );
                               double d2 = get_dist_between_rays( az_el_i2(1), az_el_i2(0), sta[ sta2 ], sta[ sta3 ], 
                                                                  _turb_sta[ sta2 ].h, _turb_sta[ sta3 ].h );

                               VCM(I,J) = _calc_matern_covariance( tau, az_el_i1, az_el_j2, sta1, sta4, spectral_coeffs, d1 )
                                        + _calc_matern_covariance( tau, az_el_i2, az_el_j1, sta2, sta3, spectral_coeffs, d2 );
                            }
                            else if( sta1 == sta3 )
                            {
                               double d = get_dist_between_rays( az_el_i1(1), az_el_i1(0), sta[ sta1 ], sta[ sta3 ], 
                                                                  _turb_sta[ sta1 ].h, _turb_sta[ sta3 ].h );

                               VCM(I,J) = _calc_matern_covariance( tau, az_el_i1, az_el_j1, sta1, sta3, spectral_coeffs, d );
                            } 
                            else if( sta1 == sta4 )
                            {
                               double d = get_dist_between_rays( az_el_i1(1), az_el_i1(0), sta[ sta1 ], sta[ sta4 ], 
                                                                  _turb_sta[ sta1 ].h, _turb_sta[ sta4 ].h );

                               VCM(I,J) = _calc_matern_covariance( tau, az_el_i1, az_el_j2, sta1, sta4, spectral_coeffs, d );
                            }
                            else if( sta2 == sta3 )
                            {
                               double d = get_dist_between_rays( az_el_i2(1), az_el_i2(0), sta[ sta2 ], sta[ sta3 ], 
                                                                  _turb_sta[ sta2 ].h, _turb_sta[ sta3 ].h );

                               VCM(I,J) = _calc_matern_covariance( tau, az_el_i2, az_el_j1, sta2, sta3, spectral_coeffs, d );
                            }
                            else if( sta2 == sta4 )
                            {
                               double d = get_dist_between_rays( az_el_i2(1), az_el_i2(0), sta[ sta2 ], sta[ sta4 ], 
                                                                  _turb_sta[ sta2 ].h, _turb_sta[ sta4 ].h );

                               VCM(I,J) = _calc_matern_covariance( tau, az_el_i2, az_el_j2, sta2, sta4, spectral_coeffs, d );
                            }
                            else
                               VCM(I,J) = 0.0;
                        }                
                    }               
                    J++;

                } // #j
                counter2 += scans.at(l).get_nobs();
            } // #l
        } // #i

        counter1 += scans.at(k).get_nobs();
     } // #k

    // fill lower left part of the covariance matrix
    ivg::Matrix D = VCM.diag();
    VCM = VCM + VCM.transpose() - D.diag();   

    // calculate correlation matrix
    ivg::Matrix ones( nobs, 1, 1.0 );
    ivg::Matrix tmp1 = ( ones.div_elem( ( VCM.diag() ).sqrt() ) ).diag();
    C = tmp1 * VCM * tmp1 ;

#if DEBUG_TROP >=1
   cerr << "--- ivg::Matrix Turbulence::calc_matern_vcm_model( int, std::vector<ivg::Scan>, ivg::Trf*, "  
        << "std::map< std::string, turbulence_data >, ivg::Matrix, ivg::Matrix & )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  
          
      return VCM;
}

// .....................................................................
double Turbulence::_calc_matern_covariance( double tau, ivg::Matrix az_el_i, ivg::Matrix az_el_j, 
                                            std::string staA, std::string staB, 
                                            ivg::Matrix spectral_coeffs, double dist_between_rays )
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ double Turbulence::_calc_matern_covariance( double, ivg::Matrix, ivg::Matrix, ivg::Matrix, double )" << endl; 
   tictoc tim;
   tim.tic();
#endif  

   double covar, besselZ, Kv;

   // calculate wavenumber of the outer scale length [1/m]
   double k0 = 2.0 * M_PI / _L0;
      
   // calculate elecromagnetic wavenumber
   double lambda = ivg::c / ivg::freq_x ;
   double elec_wn = 2 * M_PI / lambda ;

   double besselX = 5.0/6.0;

   double a = spectral_coeffs(0);
   double c = spectral_coeffs(2);

   // get structure constant
   //double Cn2 = pow( _turb_data.Cn, 2.0 );
              	      
   if( dist_between_rays < 6000.0 )
   {	
      // calculate modified bessel function of second kind
      // besselZ = k0 * _turb_data.v * tau / a ;
      besselZ = k0 * _turb_sta[ staA ].v * tau / a ; 
      Kv = boost::math::cyl_bessel_k( besselX, besselZ );
   
      covar = 0.7772 * ( pow(elec_wn,2.0) * _turb_sta[ staA ].h * _turb_sta[ staA ].Cn * _turb_sta[ staB ].Cn * c * pow(k0,-5.0/3.0) ) 
              / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) ) * pow( k0*_turb_sta[ staA ].v*tau/a, 5.0/6.0 ) * Kv;		
      //covar = 0.7772 * ( pow(elec_wn,2.0) * _turb_data.h * Cn2 * c * pow(k0,-5.0/3.0) ) / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) ) * pow( k0*_turb_data.v*tau/a, 5.0/6.0 ) * Kv;		
   }
   else if( dist_between_rays > 6000.0 && dist_between_rays <= 10000.0 )
   {
      // calculate modified structure constant according to Kermarrec and Schoen (2014), p. 1067
      //Cn2 *= 1e-1 ;		    
       
      // calculate modified bessel function of second kind
      besselZ = 2.0 * M_PI * _turb_sta[ staA ].v * tau / ( a* dist_between_rays ) ;
      Kv = boost::math::cyl_bessel_k( besselX, besselZ );

      covar = 0.7772 * ( pow(elec_wn,2.0) * _turb_sta[ staA ].h * _turb_sta[ staA ].Cn * _turb_sta[ staB ].Cn * c * pow(2.0*M_PI/dist_between_rays,-5.0/3.0) ) 
                     / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) ) * pow( 2.0*M_PI*_turb_sta[ staA ].v*tau/ (a*dist_between_rays), 5.0/6.0 ) * Kv ;					    
      //covar = 0.7772 * ( pow(elec_wn,2.0) * _turb_data.h * Cn2 * c * pow(2.0*M_PI/dist_between_rays,-5.0/3.0) ) / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) ) 
      //               * pow( 2.0*M_PI*_turb_data.v*tau/ (a*dist_between_rays), 5.0/6.0 ) * Kv ;					    
   }
   else
   {
      covar = 0.0;
   }

#if DEBUG_TROP >=2
   cerr << "--- double Turbulence::_calc_matern_covariance( double, ivg::Matrix, ivg::Matrix, ivg::Matrix, double )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  

   return covar;
}


// .....................................................................
ivg::Matrix Turbulence::calc_station_matern_vcm_model( int nobs, std::vector<ivg::Scan> scans, ivg::Trf *trf, 
                                                       ivg::Matrix spectral_coeffs, ivg::Matrix & C )
// .....................................................................
{
#if DEBUG_TROP >=1
   cerr << "+++ ivg::Matrix Turbulence::calc_station_matern_vcm_model( int, std::vector<ivg::Scan>, "  
        << "ivg::Trf*, ivg::Matrix, ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
    
    ivg::Matrix VC( nobs,nobs,0.0 );
    ivg::Matrix VCM( nobs,nobs,0.0 );

    double a = spectral_coeffs(0);
    double c = spectral_coeffs(2);

    // wavenumber of the outer scale length [1/m]
    double k0 = 2.0 * M_PI / _L0;

    // calculate elecromagnetic wavenumber
    double lambda = ivg::c / ivg::freq_x ;
    double elec_wn = 2.0 * M_PI / lambda ;

    double tau;
    ivg::Obs* obs_ptr1;
    ivg::Obs* obs_ptr2;
    std::string sta1, sta2, sta3, sta4;

    // calculate a priori station positions
    std::map< std::string, ivg::Matrix > sta_vcm;

    std::vector<Analysis_station>::iterator iter;
    for( iter = trf->begin(); iter < trf->end(); iter++ )
    {
        
        ivg::Matrix tmp(nobs,nobs,0.0);
        VCM = tmp;

        int counter1 = 0;
        int counter2 = 0;
        int J = 0;
        int I = 0;
        
        ivg::Matrix az_el_i;
        ivg::Matrix az_el_j;        
          
        std::string sta_name = iter->get_name( ivg::staname::ivs_name );
        ivg::Matrix sta_pos = iter->calc_xyz( trf->get_reference_epoch() );
        
        // fill upper right part of the variance covariance matrix
        for( int k = 0; k < scans.size(); ++k )                    // loop over scans k
        {
            for( int i = 0; i < scans.at(k).get_nobs(); ++i )       // loop over obs i
            {        
                counter2 = 0;
                I = i+counter1;

                // get pointer to current observation i
                obs_ptr1 = scans.at(k).get_obs_ptr(i);
                obs_ptr1->get_station_names( sta1,sta2 );

                // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
                if( sta_name == sta1 )
                    az_el_i = obs_ptr1->get_az_el(1);
                else if( sta_name == sta2 )
                    az_el_i = obs_ptr1->get_az_el(2);

                for( int l = 0; l < scans.size(); ++l )              // loop over scans l
                {   
                    // calculate time differences [hours]
//                    tau = abs( scans.at(k).get_epoch().get_double_mjd() - scans.at(l).get_epoch().get_double_mjd() ) * 24.0;
                    tau = abs( scans.at(l).get_epoch().get_double_mjd() - scans.at(k).get_epoch().get_double_mjd() ) * 86400.0;
                    if( tau == 0.0 )
                       tau = 1e-16;
                    
//                    cerr << ", diff: " << abs(scans.at(l).get_epoch().get_double_mjd()-scans.at(k).get_epoch().get_double_mjd()) *24.0
//                         << ", diff: " << abs(scans.at(l).get_epoch().get_double_mjd()-scans.at(k).get_epoch().get_double_mjd()) * 86400.0 
//                         << ", tau: " << tau << endl;

                    for( int j = 0; j < scans.at(l).get_nobs(); ++j ) // loop over obs j
                    {                        
                        J = j+counter2;
 
                        if( VCM(J,I) == 0.0 ) // only fill upper right triangle to reduce computational cost
                        {
                            // get pointer to current observation j
                            obs_ptr2 = scans.at(l).get_obs_ptr(j);
                            obs_ptr2->get_station_names( sta3,sta4 );                        

                            // only fill VCM if the a station joins the current observation
                            if( sta_name == sta1 || sta_name == sta2 )
                            {                   
                                // VARIANCES
                                if ( I == J )   
                                {
                                    VCM(I,J) = 0.782e3 * ( pow(elec_wn,2.0) * _turb_sta[ sta_name ].h * pow(_turb_sta[sta_name].Cn,2.0) * c * pow(k0,-5.0/3.0) ) / ( pow( sin( az_el_i(1) ), 2.0) );
//                                    VCM(I,J) = 0.782 * ( pow(elec_wn,2.0) * _turb_sta[ sta_name ].h * pow(_turb_sta[sta_name].Cn,2.0) * c * pow(k0,-5.0/3.0) ) / ( pow( sin( az_el_i(1) ), 2.0) );
                                }
                                // COVARIANCES
                                else   
                                {
                                    // (*) check which observations are correlated (in time); 
                                    //     spatial correlation is considered by the separation distance (= distance between two rays)  
                                    // (1) calcualte separation distance between two rays
                                    // (2) calculate covariance between observations i and j 
                                    ivg::Matrix sta_pos2;
                                    if( sta_name == sta3 )
                                    {
                                        az_el_j = obs_ptr2->get_az_el(1);
                                        sta_pos2 = scans.at(l).get_sta_ptr( obs_ptr2->get_scan_idx(2) )->calc_xyz( trf->get_reference_epoch() );
                                        
                                        double d = get_dist_between_rays( az_el_i(1), az_el_i(0), sta_pos, sta_pos2, 
                                                                          _turb_sta[ sta1 ].h, _turb_sta[ sta2 ].h ); 
                                        VCM(I,J) = _calc_station_matern_covariance( tau, az_el_i, az_el_j, sta_name, spectral_coeffs, d );                        
                                    }
                                    else if( sta_name == sta4 )
                                    {
                                        az_el_j = obs_ptr2->get_az_el(2);
                                        sta_pos2 = scans.at(l).get_sta_ptr( obs_ptr2->get_scan_idx(1) )->calc_xyz( trf->get_reference_epoch() );

                                        double d = get_dist_between_rays( az_el_i(1), az_el_i(0), sta_pos, sta_pos2, 
                                                                          _turb_sta[ sta1 ].h, _turb_sta[ sta2 ].h ); 
                                        VCM(I,J) = _calc_station_matern_covariance( tau, az_el_i, az_el_j, sta_name, spectral_coeffs, d ); 
                                    }                      
                                }
                            } 
                        }   
                        
                        J++;
                        
                    } // #j
                    counter2 += scans.at(l).get_nobs();
                } // #l
            } // #i
            counter1 += scans.at(k).get_nobs();
        } // #k

        // fill lower left part of the covariance matrix
        ivg::Matrix D = VCM.diag();
        VCM = VCM + VCM.transpose() - D.diag();   

        //VCM.save_bin( "/home/halsig/ascot/output/turbMatern_"+sta_name);
        
        VC += VCM;
        sta_vcm[ sta_name ] = VCM;         
    } // stations      

    // calculate correlation matrix
    ivg::Matrix ones( nobs, 1, 1.0 );
    ivg::Matrix tmp1 = ( ones.div_elem( ( VC.diag() ).sqrt() ) ).diag();
    C = tmp1 * VC * tmp1 ;

#if DEBUG_TROP >=1
   cerr << "--- ivg::Matrix Turbulence::calc_station_matern_vcm_model( int, std::vector<ivg::Scan>,"  
        << " ivg::Trf*, ivg::Matrix, ivg::Matrix & )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  
          
    return VC;
}

// .....................................................................
double Turbulence::_calc_station_matern_covariance( double tau, ivg::Matrix az_el_i, ivg::Matrix az_el_j, 
                                  std::string station, ivg::Matrix spectral_coeffs, double dist_between_rays )
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ double Turbulence::_calc_station_matern_covariance( double, ivg::Matrix, ivg::Matrix, "
           "std::string, ivg::Matrix, double )" << endl; 
   tictoc tim;
   tim.tic();
#endif  

    double covar, besselZ, Kv;

    // calculate wavenumber of the outer scale length [1/m]
    double k0 = 2.0 * M_PI / _L0;

    // calculate elecromagnetic wavenumber
    double lambda = ivg::c / ivg::freq_x ;
    double elec_wn = 2 * M_PI / lambda ;

    double besselX = 5.0/6.0;

    double a = spectral_coeffs(0);
    double c = spectral_coeffs(2);

//    besselZ = k0 * _turb_sta[ station ].v * tau / a ; 
//    Kv = boost::math::cyl_bessel_k( besselX, besselZ );
//    cerr << " teil1: " << pow( k0*_turb_sta[ station ].v*tau/a, 5.0/6.0 ) * Kv << endl;
    
    if( dist_between_rays < 6000.0 )
    {	
        // calculate modified bessel function of second kind
        besselZ = k0 * _turb_sta[ station ].v * tau / a ; 
        Kv = boost::math::cyl_bessel_k( besselX, besselZ );

        covar = 0.7772e3 * ( pow(elec_wn,2.0) * _turb_sta[ station ].h * pow(_turb_sta[ station ].Cn,2.0) * c * pow(k0,-5.0/3.0) ) 
                / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) ) * pow( k0*_turb_sta[ station ].v*tau/a, 5.0/6.0 ) * Kv;		

//        cerr << "TURB::: COV::: STA: " << station << ", " << dist_between_rays*1e-3 << " [km], tau: " << tau << ", " ;        
//        cerr << " FALL < 6km " << covar << endl;
        
        
    }
    else if( dist_between_rays > 6000.0 && dist_between_rays <= 10000.0 )
    {
        // calculate modified structure constant according to Kermarrec and Schoen (2014), p. 1067
        double Cn_mod = _turb_sta[ station ].Cn * 1e-1 ;

        // calculate modified bessel function of second kind
        besselZ = 2.0 * M_PI * _turb_sta[ station ].v * tau / ( a* dist_between_rays ) ;
        Kv = boost::math::cyl_bessel_k( besselX, besselZ );

        covar = 0.7772e3 * ( pow(elec_wn,2.0) * _turb_sta[ station ].h * pow(_turb_sta[ station ].Cn,2.0) * 1e-1 * c * pow(2.0*M_PI/dist_between_rays,-5.0/3.0) ) 
                       / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) ) * pow( 2.0*M_PI*_turb_sta[ station ].v*tau/ (a*dist_between_rays), 5.0/6.0 ) * Kv ;        
        
//        covar = 0.7772 * ( pow(elec_wn,2.0) * _turb_sta[ station ].h * pow(Cn_mod,2.0) * c * pow(2.0*M_PI/dist_between_rays,-5.0/3.0) ) 
//                       / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) ) * pow( 2.0*M_PI*_turb_sta[ station ].v*tau/ (a*dist_between_rays), 5.0/6.0 ) * Kv ;
        
        cerr << "TURB::: COV::: STA: " << station << ", " << dist_between_rays*1e-3 << " [km], tau: " << tau << ", " ;
        cerr << " FALL 6-10km " << covar << endl;
                
    }
    else
    {
        covar = 0.0;
        
//        cerr << " FALL >10km " << covar << endl;
//        cerr << endl;
    }

#if DEBUG_TROP >=2
   cerr << "--- double Turbulence::_calc_station_matern_covariance( double, ivg::Matrix, ivg::Matrix,"
           "std::string, ivg::Matrix, double )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  

   return covar;
}

// .....................................................................
double Turbulence::get_dist_between_rays( double el, double az, ivg::Matrix sta1,
                                          ivg::Matrix sta2, double h1, double h2 )
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ double Turbulence::get_dist_between_rays( double, double, ivg::Matrix, ivg::Matrix )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
   
   // transform station coordinates to virtual station coordinates at height h
   ivg::Matrix dsta1 = sta1 / ((sta1.transpose()*sta1).sqrt())(0);
   ivg::Matrix dsta2 = sta2 / ((sta2.transpose()*sta2).sqrt())(0);
   
   sta1 = sta1 + dsta1* h1;
   sta2 = sta2 + dsta2* h2;
   
   // calculate baseline length between station 1 and 2
   ivg::Matrix bl = sta1 - sta2;
   
   // (1) calculate local source vector (homogeneous coordinates) using elevation and azimuth
   double kx = cos(el) * cos(az);
   double ky = cos(el) * sin(az);
   double kz = sin(el);
   
   ivg::Matrix k_homog( 1,4, 1.0 );
   k_homog( 0,0 ) = kx;
   k_homog( 0,1 ) = ky;
   k_homog( 0,2 ) = kz;
   
   // (2) transform to global system
   ivg::Matrix T1( 1,3,0.0 ); T1( 0,2 ) = 1.0 ;
   ivg::Matrix T2( 1,3,0.0 ); T2( 0,0 ) = 1.0 ;
   ivg::Matrix tmp = sta1;
   tmp(2,0) = 0.0;

   // calculate rotation angles
   double beta = acos( ( T1 * sta1 )(0) / (sta1.norm())(0) ) ;
   double delta = M_PI - acos( ( T2 * tmp )(0) / (tmp.norm())(0) );
   
   // set rotation matrices: R1 = rotation about y-axis; R2 = rotation about z-axis
   ivg::Matrix R1( 4,4,0.0 );
   R1( 0,0 ) = cos( beta );  R1( 0,1 ) = 0.0;  R1( 0,2 ) = sin( beta );  R1( 0,3 ) = 0.0 ;
   R1( 1,0 ) = 0.0; 	     R1( 1,1 ) = 1.0;  R1( 1,2 ) = 0.0; 	 R1( 1,3 ) = 0.0 ;
   R1( 2,0 ) = -sin( beta ); R1( 2,1 ) = 0.0;  R1( 2,2 ) = cos( beta );  R1( 2,3 ) = 0.0 ;
   R1( 3,0 ) = 0.0; 	     R1( 3,1 ) = 0.0;  R1( 3,2 ) = 0.0; 	 R1( 3,3 ) = 1.0 ;
   
   ivg::Matrix R2( 4,4,0.0 );
   R2( 0,0 ) = cos( delta );  R2( 0,1 ) = -sin( delta ); R2( 0,2 ) = 0.0;  R2( 0,3 ) = 0.0 ;
   R2( 1,0 ) = sin( delta );  R2( 1,1 ) = cos( delta );  R2( 1,2 ) = 0.0;  R2( 1,3 ) = 0.0 ;
   R2( 2,0 ) = 0.0;  	      R2( 2,1 ) = 0.0;  	 R2( 2,2 ) = 1.0;  R2( 2,3 ) = 0.0 ;
   R2( 3,0 ) = 0.0; 	      R2( 3,1 ) = 0.0;  	 R2( 3,2 ) = 0.0;  R2( 3,3 ) = 1.0 ;
   
   // transform local source vector to global (cartesian) system
   ivg::Matrix M = R2.transpose() * R1.transpose() ;
   ivg::Matrix ks = M * k_homog.transpose() ;
   
   // "re-transformation" from homogeneous coordinates
   ivg::Matrix k = ks.get_sub( 0, 0, 2, ks.size(2)-1 );
   
   // (3) calculate separation distance between two rays
   ivg::Matrix tau = bl.transpose() * k;
   ivg::Matrix ones( tau.size(2), 1, 1.0 );
   	
   double dist = ( ( (ones* (bl.transpose()* bl)) - tau.transpose().mult_elem( tau.transpose() ) ).sqrt() )(0);
   	
#if DEBUG_TROP >=2
   cerr << "--- double Turbulence::get_dist_between_rays( double, double, ivg::Matrix, ivg::Matrix )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  

   return dist;
} 
  
  
/*
// .....................................................................
ivg::Matrix Turbulence::get_dist_between_rays( ivg::Matrix el, ivg::Matrix az, ivg::Matrix xyz_sta1, ivg::Matrix xyz_sta2 )
// .....................................................................
{ 
         #ifdef VERBOSE
    	cout << "INFO:: Methode 'ivg::Matrix Turbulence::get_dist_between_rays( ivg::Matrix el, ivg::Matrix az, ivg::Matrix xyz_sta1, ivg::Matrix xyz_sta2 )' aufgerufen"<< endl;
         cout << endl;
         #endif
         
         // calculate baseline length between station 1 and 2
         ivg::Matrix bl = xyz_sta1 - xyz_sta2;
         
         // (1) calculate local source vector (homogeneous coordinates) using elevation and azimuth
         ivg::Matrix kx = el.cos().mult_elem( az.cos() );
         ivg::Matrix ky = el.cos().mult_elem( az.sin() );
         ivg::Matrix kz = el.sin();
         
         ivg::Matrix k_homog = kx;
         k_homog.append_cols( ky );
         k_homog.append_cols( kz );
         ivg::Matrix ones( k_homog.size(1), 1, 1.0);
         k_homog.append_cols( ones );
         
         // (2) transform to global system
         ivg::Matrix T1( 1,3,0.0 ); T1( 0,2 ) = 1.0 ;
         ivg::Matrix T2( 1,3,0.0 ); T2( 0,0 ) = 1.0 ;
	ivg::Matrix tmp = xyz_sta1;
	tmp(2,0) = 0.0;
	
	// calculate rotation angles
	double beta = acos( ( T1 * xyz_sta1 )(0) / (xyz_sta1.norm())(0) ) ;
	double delta = M_PI - acos( ( T2 * tmp )(0) / (tmp.norm())(0) );
	
	// set rotation matrices: R1 = rotation about y-axis; R2 = rotation about z-axis
	ivg::Matrix R1( 4,4,0.0 );
	R1( 0,0 ) = cos( beta );  R1( 0,1 ) = 0.0;  R1( 0,2 ) = sin( beta );  R1( 0,3 ) = 0.0 ;
	R1( 1,0 ) = 0.0; 	  R1( 1,1 ) = 1.0;  R1( 1,2 ) = 0.0; 	      R1( 1,3 ) = 0.0 ;
	R1( 2,0 ) = -sin( beta ); R1( 2,1 ) = 0.0;  R1( 2,2 ) = cos( beta );  R1( 2,3 ) = 0.0 ;
	R1( 3,0 ) = 0.0; 	  R1( 3,1 ) = 0.0;  R1( 3,2 ) = 0.0; 	      R1( 3,3 ) = 1.0 ;
	
	ivg::Matrix R2( 4,4,0.0 );
	R2( 0,0 ) = cos( delta );  R2( 0,1 ) = -sin( delta ); R2( 0,2 ) = 0.0;  R2( 0,3 ) = 0.0 ;
	R2( 1,0 ) = sin( delta );  R2( 1,1 ) = cos( delta );  R2( 1,2 ) = 0.0;  R2( 1,3 ) = 0.0 ;
	R2( 2,0 ) = 0.0;  	   R2( 2,1 ) = 0.0;  	      R2( 2,2 ) = 1.0;  R2( 2,3 ) = 0.0 ;
	R2( 3,0 ) = 0.0; 	   R2( 3,1 ) = 0.0;  	      R2( 3,2 ) = 0.0;  R2( 3,3 ) = 1.0 ;
	
	// transform local source vector to global (cartesian) system
	ivg::Matrix invR1 = R1;
	invR1.inv();
	
	ivg::Matrix invR2 = R2;
	invR2.inv();
	
	ivg::Matrix M = invR2 * invR1 ;
		
	ivg::Matrix ks = M * k_homog.transpose() ;
	
	// "re-transformation" from homogeneous coordinates
	ivg::Matrix k = ks.get_sub( 0, 0, 2, ks.size(2)-1 );

	// (3) calculate separation distance between two rays
	ivg::Matrix tau = bl.transpose() * k;
	ones.resize( tau.size(2), 1, 1.0 );
		
	ivg::Matrix dist = ( (ones* (bl.transpose()* bl)) - tau.transpose().mult_elem( tau.transpose() ) ).sqrt();
	  	
	return dist;
}
*/


// .....................................................................
ivg::Matrix Turbulence::calc_matern_covariance_function( int nobs, std::vector<ivg::Scan> scans, ivg::Trf *trf, 
                                                         ivg::Matrix spectral_coeffs, ivg::Matrix & C )
// .....................................................................
{
#if DEBUG_TROP >=1
   cerr << "+++ ivg::Matrix Turbulence::calc_matern_covariance_function( int, std::vector<ivg::Scan>, "  
        << "ivg::Trf*, ivg::Matrix, ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
    
    ivg::Matrix VC( nobs,nobs,0.0 );
    ivg::Matrix VCM( nobs,nobs,0.0 );

    double a = spectral_coeffs(0);
    double c = spectral_coeffs(2);

    // wavenumber of the outer scale length [1/m]
    double k0 = 2.0 * M_PI / _L0;

    // calculate elecromagnetic wavenumber
    double lambda = ivg::c / ivg::freq_x ;
    double elec_wn = 2.0 * M_PI / lambda ;

    double tau;
    ivg::Obs* obs_ptr1;
    ivg::Obs* obs_ptr2;
    std::string sta1, sta2, sta3, sta4;

    // calculate a priori station positions
    std::map< std::string, ivg::Matrix > sta_vcm;

    std::vector<Analysis_station>::iterator iter;
    for( iter = trf->begin(); iter < trf->end(); iter++ )
    {
        
        ivg::Matrix tmp(nobs,nobs,0.0);
        VCM = tmp;

        int counter1 = 0;
        int counter2 = 0;
        int J = 0;
        int I = 0;
        
        ivg::Matrix az_el_i;
        ivg::Matrix az_el_j;        
          
        std::string sta_name = iter->get_name( ivg::staname::ivs_name );
        ivg::Matrix sta_pos = iter->calc_xyz( trf->get_reference_epoch() );
        
        // fill upper right part of the variance covariance matrix
        for( int k = 0; k < scans.size(); ++k )                    // loop over scans k
        {
            for( int i = 0; i < scans.at(k).get_nobs(); ++i )       // loop over obs i
            {        
                counter2 = 0;
                I = i+counter1;

                // get pointer to current observation i
                obs_ptr1 = scans.at(k).get_obs_ptr(i);
                obs_ptr1->get_station_names( sta1,sta2 );

                // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
                if( sta_name == sta1 )
                    az_el_i = obs_ptr1->get_az_el(1);
                else if( sta_name == sta2 )
                    az_el_i = obs_ptr1->get_az_el(2);

                for( int l = 0; l < scans.size(); ++l )              // loop over scans l
                {   
                    // calculate time differences [hours]
                    tau = abs( scans.at(l).get_epoch().get_double_mjd() - scans.at(k).get_epoch().get_double_mjd() ) * 86400.0;
                    if( tau == 0.0 )
                       tau = 1e-16;

                    for( int j = 0; j < scans.at(l).get_nobs(); ++j ) // loop over obs j
                    {                        
                        J = j+counter2;
 
                        if( VCM(J,I) == 0.0 ) // only fill upper right triangle to reduce computational cost
                        {
                            // get pointer to current observation j
                            obs_ptr2 = scans.at(l).get_obs_ptr(j);
                            obs_ptr2->get_station_names( sta3,sta4 );                        

                            // only fill VCM if the a station joins the current observation
                            if( sta_name == sta1 || sta_name == sta2 )
                            {                   
                                // VARIANCES
                                if ( I == J )   
                                {
                                    VCM(I,J) = 0.782e3 * ( pow(elec_wn,2.0) * _turb_sta[ sta_name ].h * pow(_turb_sta[sta_name].Cn,2.0) 
                                             * c * pow(k0,-5.0/3.0) ); // / ( pow( sin( az_el_i(1) ), 2.0) );
                                }
                                // COVARIANCES
                                else   
                                {
                                    // (*) check which observations are correlated (in time); 
                                    //     spatial correlation is considered by the separation distance (= distance between two rays)  
                                    // (1) calcualte separation distance between two rays
                                    // (2) calculate covariance between observations i and j 
                                    if( sta_name == sta3 || sta_name == sta4 )
                                    {
                                        // calculate modified bessel function of second kind
                                        double besselX = 5.0/6.0;
                                        double besselZ = k0 * _turb_sta[ sta_name ].v * tau / a ; 
                                        double Kv = boost::math::cyl_bessel_k( besselX, besselZ );

                                        VCM(I,J) = 0.7772e3 * ( pow(elec_wn,2.0) * _turb_sta[ sta_name ].h * pow(_turb_sta[ sta_name ].Cn,2.0) 
                                                 * c * pow(k0,-5.0/3.0) ) * pow( k0*_turb_sta[ sta_name ].v*tau/a, 5.0/6.0 ) * Kv;
                                    }                     
                                }
                            } 
                        }   
                        
                        J++;
                        
                    } // #j
                    counter2 += scans.at(l).get_nobs();
                } // #l
            } // #i
            counter1 += scans.at(k).get_nobs();
        } // #k

        // fill lower left part of the covariance matrix
        ivg::Matrix D = VCM.diag();
        VCM = VCM + VCM.transpose() - D.diag();   

        VCM.save_bin( "/home/halsig/ascot/output/turbMatern_"+sta_name);
        
        VC += VCM;
        sta_vcm[ sta_name ] = VCM;         
    } // stations      

    // calculate correlation matrix
    ivg::Matrix ones( nobs, 1, 1.0 );
    ivg::Matrix tmp1 = ( ones.div_elem( ( VC.diag() ).sqrt() ) ).diag();
    C = tmp1 * VC * tmp1 ;

#if DEBUG_TROP >=1
   cerr << "--- ivg::Matrix Turbulence::calc_matern_covariance_function( int, std::vector<ivg::Scan>,"  
        << " ivg::Trf*, ivg::Matrix, ivg::Matrix & )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  
          
    return VC;
}


// ****************************************************************************
// ********** (2) B R U N N E R / S C H OE N   O R I G.   M O D E L ***********
// ****************************************************************************
//// .....................................................................
//ivg::Matrix Turbulence::calc_sigma_c_model ( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, ivg::Matrix spectral_coeffs, double dh, ivg::Matrix & CORR ) 
//// .....................................................................
//{
//#if DEBUG_TROP >=1
//   cerr << "+++ ivg::Matrix Turbulence::calc_sigma_c_model ( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, ivg::Matrix & )" << endl; 
//   tictoc tim;
//   tim.tic();
//#endif  
//	
//	// initialization
//	ivg::Matrix VCM( nobs, nobs, 0.0 );
//	ivg::Matrix ri(3,1,0.0);
//	ivg::Matrix rj(3,1,0.0);
//	double rN, rE, rU, rX, rY, rZ;
//	double pk0, p ;
//	double factor, h1, h2, h3, z, f, sum ;
//	double tau_A, tau_B;
//	double Kv;
//	double besselX, besselZ;
//	
//	// wind field
//	double alpha_V = 230.0 / 180.0 * M_PI;
//	double eps_V = 0.0 ;
//	
//	double alpha_AB = 0.0 ;
//	double eps_AB = 0.0 ;
//	double rho_AB = 0.0 ;
//	
//	double t_start = 0.0 ;
//	
//	// ellipsoidal parameters of the anisotropic spectum [-]
//	double a = spectral_coeffs(0);	
//	double b = spectral_coeffs(1);	
//	double c = spectral_coeffs(2);	
//	
//        // rotation matrix Rv
//	ivg::Matrix Rv( 3, 3, 0.0 );
//	Rv(0,0) =  cos(alpha_V)*cos(eps_V);
//	Rv(0,1) = -sin(alpha_V);   
//	Rv(0,2) = -cos(alpha_V)*sin(eps_V) ;
//	Rv(1,0) =  sin(alpha_V)*cos(eps_V);
//	Rv(1,1) =  cos(alpha_V);
//	Rv(1,2) = -sin(alpha_V)*sin(eps_V) ;
//	Rv(2,0) =  sin(eps_V);
//	Rv(2,1) =  0.0 ;
//	Rv(2,2) =  cos(eps_V) ;	
//	Rv = Rv.transpose();
//	
//	ivg::Matrix v_vec(3,1,0.0);
//
//        int counter1 = 0;
//        int counter2 = 0;
//        int J = 0;
//        int I = 0;
//	
//        ivg::Matrix az_el_i1;
//        ivg::Matrix az_el_i2;
//        ivg::Matrix az_el_j1;
//        ivg::Matrix az_el_j2;
//
//        ivg::Obs* obs_ptr1;
//        ivg::Obs* obs_ptr2;
//
//        std::string sta1, sta2, sta3, sta4;
//
//	// constants
//	double k0 = 2.0 * M_PI / _L0 ;
//	
//	// integration constants
//	double dz1, dz2; 
//
//	ivg::Matrix height( 0.0, dh, _turb_data.h, 1); 
//	height( height.size(1)-1, 0 ) = _turb_data.h ;
//	
//	double Cn_coeff = pow( _turb_data.Cn, 2.0 ) ;
//	Cn_coeff *= 1.0e6 ;
//	
//
//        for( int k = 0; k <= scans.size()-1; k++ )                    // loop over scans k
//        {
//           tau_A = ( scans.at(k).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0  - t_start ;
//
//           for( int i = 0; i <= scans.at(k).get_nobs()-1; i++ )       // loop over obs i
//           {
//              counter2 = 0;
//              I = i+counter1;
//              J = I; 
//            
//              // get pointer to current observation i
//              obs_ptr1 = scans.at(k).get_obs_ptr(i);
//                     
//               // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
//               az_el_i1 = obs_ptr1->get_az_el(1);
//               az_el_i2 = obs_ptr1->get_az_el(2);
//               obs_ptr1->get_station_names( sta1,sta2 );
//
//	       // components of the topocentric vector used for integration coordinates
//               ri(0,0) = cos( az_el_i1(0) ) * cos( az_el_i1(1) );
//	       ri(1,0) = sin( az_el_i1(0) ) * cos( az_el_i1(1) );
//	       ri(2,0) = sin( az_el_i1(1) );
//
//
//              for( int l = k; l <= scans.size()-1; l++ )              // loop over scans l
//              {   
//                 tau_B = ( scans.at(l).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0  - t_start ;
//
//                 for( int j = i; j <= scans.at(l).get_nobs()-1; j++ ) // loop over obs j
//                 {		
//
//                    // get pointer to current observation j
//                    obs_ptr2 = scans.at(l).get_obs_ptr(j);
//                 
//                    // get azimuth and elevation as well as station name of station 1 and 2 for observation j  
//                    az_el_j1 = obs_ptr2->get_az_el(1);
//                    az_el_j2 = obs_ptr2->get_az_el(2);
//                    obs_ptr2->get_station_names( sta3,sta4 );
//
//                        // VARIANCES
//			if ( I == J )	
//			{				
//				// rotation of the topocentric vector to the wind speed system
//				rX =  cos( alpha_V ) * cos( eps_V )*ri(0,0)  +  sin( alpha_V ) * cos(eps_V) * ri(1,0)     + sin( eps_V ) * ri(2,0);
//				rY = -sin( alpha_V ) * ri(0,0)               +  cos( alpha_V ) * ri(1,0);
//				rZ = -cos( alpha_V ) * sin( eps_V )*ri(0,0)  -  sin( alpha_V ) * sin( eps_V ) * ri(1,0)   + cos( eps_V ) * ri(2,0);
//				
//				// integration distance 
//				pk0 = sqrt( pow( rX, 2.0 ) / pow( a, 2.0 ) + pow( rY, 2.0 ) / pow( b, 2.0 ) + pow( rZ, 2.0 ) / pow( c, 2.0 ) ) * k0 ;
//				
//				// factor of integral
//				factor = 4.0* pow( M_PI, 2.0 ) / sin( az_el_i1(1) ) * 0.033 * 3/5 * pow( k0, (-2.0/3.0) ) / ( pk0 * sqrt(M_PI) * std::tgamma(5.0/6.0)) * Cn_coeff;	// [m]
//				
//				// helping values
//				h1 = 2.0*M_PI / (sqrt(3.0) * std::tgamma(2.0/3.0));   		// c ~ 2.6789 [-]
//				h2 = 9 * pow( 2.0, (1.0/3.0) ) * std::tgamma(2.0/3.0) / 10;   	// d  ~1.5355 [-]
//				h3 = pk0 / sin( az_el_i1(1) );               			//            [1/m] 
//				z =  h3 * _turb_data.h;                          		//            [-] 
//				
//				// calculate hypergeometric function (F23, F12)
//				ivg::Matrix input_a(1,2); input_a(0,0) = 0.5; input_a(0,1) = 1.0;
//				ivg::Matrix input_b(1,3); input_b(0,0) = 2.0/3.0; input_b(0,1) = 3.0/2.0; input_b(0,2) = 2.0 ;
//				double input_z = 1.0 / 4.0 * pow( z, 2.0 );
//				double F23 = hypergeom( input_a, input_b, input_z ) ;
//				
//				ivg::Matrix input_a2(1,1); input_a2(0,0) = 5.0 / 6.0;
//				ivg::Matrix input_b2(1,2); input_b2(0,0) = 11.0 / 6.0; input_b2(0,1) = 7.0 / 3.0; 
//				double input_z2 = 1.0 / 4.0 * pow( z, 2.0 );
//				double F12 = hypergeom( input_a2, input_b2, input_z2 ) ;
//			
//				// closed solution of the integral				
//				VCM(I,J) = 1e2 * factor / 8.0 * pow( z, 2.0 ) / h3* ( 4.0 * h1 * F23 - 3.0 * h2 * pow( z, (2.0/3.0) ) * F12 ); 
//				
//			}
//                        // COVARIANCES
//			else 	
//			{
//                           // (*) check which observations are correlated (in time); 
//                           //     spatial correlation is considered by the separation distance (= distance between two rays)  
//                           // (1) calcualte separation distance between two rays
//                           // (2) calculate covariance between observations i and j 
//                           if( sta1 == sta3 && sta2 == sta4 )
//                              VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i1, az_el_j1, dh )
//                                       + _calc_sigma_c_covariance( tau_A, tau_B, az_el_i2, az_el_j2, dh );
//      
//                           else if( sta1 == sta4 && sta2 == sta3 )
//                              VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i1, az_el_j2, dh )
//                                       + _calc_sigma_c_covariance( tau_A, tau_B, az_el_i2, az_el_j1, dh );
//      
//                           else if( sta1 == sta3 )
//                              VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i1, az_el_j1, dh );
//      
//                           else if( sta1 == sta4 )
//                              VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i1, az_el_j2, dh );
//      
//                           else if( sta2 == sta3 )
//                              VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i2, az_el_j1, dh );
//      
//                           else if( sta2 == sta4 )
//                              VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i2, az_el_j2, dh );
//      
//                           else
//                     	        VCM(I,J) = 0.0;
//			}
//                        
//                        J++;
//			
//	         } // #j
//               
//                 if(J <= scans.at(l).get_nobs())
//                    counter2 += scans.at(l).get_nobs();
//
//	      } // #l
//           } // #i
//
//           counter1 += scans.at(k).get_nobs();
//        } // #k
//
//	// fill lower left part of the covariance matrix
//	ivg::Matrix D = VCM.diag();
//	VCM = VCM + VCM.transpose() - D.diag(); 
//	
//	// calculate correlation matrix
//	ivg::Matrix ones( nobs, 1, 1.0 );
//	ivg::Matrix tmp1 = ( ones.div_elem( ( VCM.diag() ).sqrt() ) ).diag();
//	CORR = tmp1 * VCM * tmp1 ;
//
//#if DEBUG_TROP >=1
//   cerr << "--- ivg::Matrix Turbulence::calc_sigma_c_model ( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, ivg::Matrix & )" 
//        << " : " << tim.toc() << " s " << endl; 
//#endif  
//	
//	return VCM ;
//}


// .....................................................................
ivg::Matrix Turbulence::calc_sigma_c_model ( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, ivg::Matrix spectral_coeffs, double dh, ivg::Matrix & CORR ) 
// .....................................................................
{
#if DEBUG_TROP >=1
   cerr << "+++ ivg::Matrix Turbulence::calc_sigma_c_model ( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
	
    // initialization
    ivg::Matrix VCM( nobs, nobs, 0.0 );
    ivg::Matrix VC( nobs, nobs, 0.0 );
    std::map< std::string, ivg::Matrix > sta_vcm;   

    double rX, rY, rZ;
    double pk0, p ;
    double factor, h1, h2, h3, z;
    double tau_A, tau_B;

    ivg::Obs* obs_ptr1;
    ivg::Obs* obs_ptr2;

    std::string sta1, sta2, sta3, sta4;
	
    
    std::vector<Analysis_station>::iterator iter;
    for( iter = trf->begin(); iter < trf->end(); iter++ )
    {
        
        _turb_data = _turb_sta[iter->get_name(ivg::staname::ivs_name)];
        
        ivg::Matrix tmp(nobs,nobs,0.0);
        VCM = tmp;

        int counter1 = 0;
        int counter2 = 0;
        int J = 0;
        int I = 0;
        
        ivg::Matrix az_el_i;
        ivg::Matrix az_el_j;        
          
        std::string sta_name = iter->get_name( ivg::staname::ivs_name );
       
        // fill upper right part of the variance covariance matrix
        for( int k = 0; k < scans.size(); ++k )                    // loop over scans k
        {
            tau_A = ( scans.at(k).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0  ;
            
            
            for( int i = 0; i < scans.at(k).get_nobs(); ++i )       // loop over obs i
            {        
                counter2 = 0;
                I = i+counter1;
      
                // get pointer to current observation i
                obs_ptr1 = scans.at(k).get_obs_ptr(i);
                obs_ptr1->get_station_names( sta1,sta2 );              
                
                // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
                if( sta_name == sta1 )
                    az_el_i = obs_ptr1->get_az_el(1);
                else if( sta_name == sta2 )
                    az_el_i = obs_ptr1->get_az_el(2);
                
                
                for( int l = 0; l < scans.size(); ++l )              // loop over scans l
                {   
                    tau_B = ( scans.at(l).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0 ;
                    
                    for( int j = 0; j < scans.at(l).get_nobs(); ++j ) // loop over obs j
                    {                        
                        J = j+counter2;
 
                        if( VCM(J,I) == 0.0 ) // only fill upper right triangle to reduce computational cost
                        {
                          // get pointer to current observation j
                          obs_ptr2 = scans.at(l).get_obs_ptr(j);
                          obs_ptr2->get_station_names( sta3,sta4 );                        

                            // only fill VCM if the a station joins the current observation
                            if( sta_name == sta1 || sta_name == sta2 )
                            {       
                                // VARIANCES
                                if ( I == J )	
                                {		
                                   VCM(I,J) = _calc_sigma_c_variance( az_el_i, spectral_coeffs );

                                }
                                // COVARIANCES
                                else 	
                                {
                                    // (*) check which observations are correlated (in time); 
                                    //     spatial correlation is considered by the separation distance (= distance between two rays)  
                                    // (1) calcualte separation distance between two rays
                                    // (2) calculate covariance between observations i and j                                  
                                    ivg::Matrix sta_pos2;
                                    if( sta_name == sta3 )
                                    {
                                      az_el_j = obs_ptr2->get_az_el(1);
                                      VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i, az_el_j, dh );
                                    }
                                    else if( sta_name == sta4 )
                                    {
                                      az_el_j = obs_ptr2->get_az_el(2);
                                      VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, az_el_i, az_el_j, dh );
                                    }                     
                                }
                            }
                        }
                                             
                        J++;
			
                    } // #j
               
                    counter2 += scans.at(l).get_nobs();
                    
	      } // #l
           } // #i

           counter1 += scans.at(k).get_nobs();
        } // #k

	// fill lower left part of the covariance matrix
	ivg::Matrix D = VCM.diag();
	VCM = VCM + VCM.transpose() - D.diag(); 
	
//        VCM.save_bin( "/home/halsig/ascot/output/turbVCM_"+sta_name);
        
        VC += VCM;
        sta_vcm[ sta_name ] = VCM;  
        
    } // stations
        
    // calculate correlation matrix
    ivg::Matrix ones( nobs, 1, 1.0 );
    ivg::Matrix tmp1 = ( ones.div_elem( ( VCM.diag() ).sqrt() ) ).diag();
    CORR = tmp1 * VCM * tmp1 ;

#if DEBUG_TROP >=1
   cerr << "--- ivg::Matrix Turbulence::calc_sigma_c_model ( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, ivg::Matrix & )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  
	
    return VC;
}

// .....................................................................
ivg::Matrix Turbulence::calc_sigma_c_model( ivg::Analysis_station* sta , ivg::Matrix azel, std::vector<ivg::Date> epochs, ivg::Matrix spectral_coeffs, double dh ){
// ..................................................................... 
    
    
    _turb_data = _turb_sta[sta->get_name(ivg::staname::ivs_name)];
        
    int npos = azel.rows();
    int nepochs = epochs.size();
    int nobs = npos * nepochs;

    // initialization
    Matrix VCM( nobs, nobs, 0.0 );

    double tau_A(0.0), tau_B(0.0);    
    
    /*
    std::cerr << "calc vcm" << std::string(9,' ');
    int I = -1;
    for( int epoch_idx_1 = 0; epoch_idx_1 < epochs.size(); ++epoch_idx_1 ){
        tau_A =  ( epochs[epoch_idx_1].get_double_mjd() - epochs[0].get_double_mjd() ) * 24.0;
        
        for( int pos_idx_1 = 0; pos_idx_1 < azel.rows(); ++pos_idx_1 ){
            I++;

            ivg::Matrix azel_1 = azel.get_sub (pos_idx_1, 0, pos_idx_1, 1);
            
            std::cerr << std::string(9,'\b') << setw(4) << I+1 << "/" << setw(4) << nobs;
            
            for( int epoch_idx_2 = 0; epoch_idx_2 < epochs.size(); ++epoch_idx_2 ){
                tau_B =  ( epochs[epoch_idx_2].get_double_mjd() - epochs[0].get_double_mjd() ) * 24.0;
                
                #pragma omp parallel for num_threads(3)
                for( int pos_idx_2 = 0; pos_idx_2 < azel.rows(); ++pos_idx_2 ){

                    int J = epoch_idx_2*npos + pos_idx_2;
                    if(J>=I ){
                        ivg::Matrix azel_2 =   azel.get_sub (pos_idx_2, 0, pos_idx_2, 1);
                        
                        if(J==I){                            
                            VCM(I,J) = _calc_sigma_c_variance( azel_1, spectral_coeffs );
                        } else {
                            VCM(I,J) = _calc_sigma_c_covariance( tau_A, tau_B, azel_1, azel_2, dh );
                            VCM(J,I) = VCM(I, J);
                        }

                    }
                }
            }
        }
    }
    std::cout << std::endl;
    */
    
    std::cerr << "calc vcm" << std::string(9,' ');
    std::vector<ivg::Matrix> blocks( nepochs, ivg::Matrix(npos, npos, 0.0) );
    
    for( int epoch_idx = 0; epoch_idx < nepochs; ++epoch_idx ){
        ivg::Matrix& block = blocks[epoch_idx];
        tau_A =  ( epochs[epoch_idx].get_double_mjd() - epochs[0].get_double_mjd() ) * 24.0;
        tau_B = 0.0;
//        #pragma omp parallel for num_threads(3)
        for( int pos_idx_1 = 0; pos_idx_1 < azel.rows(); ++pos_idx_1 ){
            ivg::Matrix azel_1 = azel.get_sub (pos_idx_1, 0, pos_idx_1, 1);
            for( int pos_idx_2 = 0; pos_idx_2 < azel.rows(); ++pos_idx_2 ){
 
                ivg::Matrix azel_2 =   azel.get_sub (pos_idx_2, 0, pos_idx_2, 1);

                if(epoch_idx == 0 ) {
                    if( pos_idx_1==pos_idx_2  ){
                        block(pos_idx_1, pos_idx_2) = _calc_sigma_c_variance( azel_1, spectral_coeffs );
                    } else if(pos_idx_2 > pos_idx_1) {
                        block(pos_idx_1, pos_idx_2) = _calc_sigma_c_covariance( tau_A, tau_B, azel_1, azel_2, dh );
                        block(pos_idx_2, pos_idx_1) = block(pos_idx_1, pos_idx_2);
                    }
                } else {
                    block(pos_idx_1, pos_idx_2) = _calc_sigma_c_covariance( tau_A, tau_B, azel_1, azel_2, dh );
                }
            }
        }
        std::cerr << std::string(9,'\b') << setw(4) << epoch_idx+1 << "/" << setw(4) << nepochs;
    }
    std::cout << std::endl;
    
    VCM.blockToepliz(blocks,'S');
    VCM.save_bin( "/home/corbin/VCM.dat");
    return VCM;

}

// .....................................................................
double Turbulence::_calc_sigma_c_covariance( double tau_A, double tau_B, ivg::Matrix az_el_i, ivg::Matrix az_el_j, double dh )
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ double Turbulence::_calc_sigma_c_covariance( double, ivg::Matrix, ivg::Matrix, ivg::Matrix, double )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
   
    ivg::Matrix ri(3,1,0.0);
    ivg::Matrix rj(3,1,0.0);
    double factor, p, sum ;
    double besselX, besselZ, Kv;

    // wind field
    double alpha_V = 230.0 / 180.0 * M_PI;
    double eps_V = 0.0 ;

    // rotation matrix Rv
    ivg::Matrix Rv( 3, 3, 0.0 );
    Rv(0,0) =  cos(alpha_V)*cos(eps_V);
    Rv(0,1) = -sin(alpha_V);   
    Rv(0,2) = -cos(alpha_V)*sin(eps_V) ;
    Rv(1,0) =  sin(alpha_V)*cos(eps_V);
    Rv(1,1) =  cos(alpha_V);
    Rv(1,2) = -sin(alpha_V)*sin(eps_V) ;
    Rv(2,0) =  sin(eps_V);
    Rv(2,1) =  0.0 ;
    Rv(2,2) =  cos(eps_V) ;	
    Rv = Rv.transpose();

    ivg::Matrix v_vec(3,1,0.0);

    // constants
    double k0 = 2.0 * M_PI / _L0 ;

    // integration constants
    double dz1, dz2; 

    ivg::Matrix height( 0.0, dh, _turb_data.h, 1); 
    height( height.size(1)-1, 0 ) = _turb_data.h ;

    double Cn_coeff = pow( _turb_data.Cn, 2.0 ) ;
    Cn_coeff *= 1.0e6 ;


    // components of the topocentric vector used for integration coordinates
    ri(0,0) = cos( az_el_i(0) ) * cos( az_el_i(1) );
    ri(1,0) = sin( az_el_i(0) ) * cos( az_el_i(1) );
    ri(2,0) = sin( az_el_i(1) );

    rj(0,0) = cos( az_el_j(0) ) * cos( az_el_j(1) );
    rj(1,0) = sin( az_el_j(0) ) * cos( az_el_j(1) );
    rj(2,0) = sin( az_el_j(1) );

    v_vec(0,0) = sin(alpha_V) * cos(eps_V);
    v_vec(1,0) = cos(alpha_V) * cos(eps_V);
    v_vec(2,0) = sin(eps_V);		

    v_vec *= _turb_data.v ;
    ivg::Matrix vt = v_vec * abs(tau_B - tau_A);

    // factor of integral
    factor = Cn_coeff * 4.0* M_PI * sqrt(M_PI) /  std::tgamma(5.0/6.0) * 0.033 * 3/5 ;

    // >>> numerical integration
    sum = 0.0 ;

    for ( int z1 = 0; z1 < height.size(1); z1++ )
    {
        for ( int z2 = 0; z2 < height.size(1); z2++ )
        {
            dz1 = height(z1,0) / sin( az_el_i(1) );
            dz2 = height(z2,0) / sin( az_el_j(1) );

            // calculate integration distance
            ivg::Matrix d_vec = rj * dz2  - vt - ri* dz1 ;
            ivg::Matrix d = ( d_vec.transpose() * Rv * Rv.transpose() * d_vec ).sqrt();							
            p = d(0,0);	

            // calculate modified bessel function of second kind
            besselX = -1.0/3.0;
            besselZ = p * k0 ;

            if( besselZ == 0.0 )
               besselZ = 1e-16;

            Kv = boost::math::cyl_bessel_k( besselX, besselZ );

            // solution of integral			
            sum += factor * pow( ( p / (2.0*k0) ), (1.0/3.0) ) * Kv / ( sin( az_el_i(1) ) * sin( az_el_j(1) ) );
        }
    }			
    double covar = sum ;

#if DEBUG_TROP >=2
   cerr << "--- double Turbulence::_calc_sigma_c_covariance( double, ivg::Matrix, ivg::Matrix, ivg::Matrix, double )"    
        << " : " << tim.toc() << " s " << endl; 
#endif 

   return covar;
}

double Turbulence::_calc_sigma_c_variance( ivg::Matrix az_el_i, ivg::Matrix spectral_coeffs ){
        
    ivg::Matrix ri(3,1,0.0);

    double rX, rY, rZ;
    double pk0 ;
    double factor, h1, h2, h3, z;

    // wind field
    double alpha_V = 230.0 / 180.0 * M_PI;
    double eps_V = 0.0 ;

    // ellipsoidal parameters of the anisotropic spectum [-]
    double a = spectral_coeffs(0);	
    double b = spectral_coeffs(1);	
    double c = spectral_coeffs(2);	

    // constants
    double k0 = 2.0 * M_PI / _L0 ;

    double Cn_coeff = pow( _turb_data.Cn, 2.0 ) ;
    Cn_coeff *= 1.0e6 ;
    
    // components of the topocentric vector used for integration coordinates
    ri(0,0) = cos( az_el_i(0) ) * cos( az_el_i(1) );
    ri(1,0) = sin( az_el_i(0) ) * cos( az_el_i(1) );
    ri(2,0) = sin( az_el_i(1) );                                                                

    // rotation of the topocentric vector to the wind speed system
    rX =  cos( alpha_V ) * cos( eps_V )*ri(0,0)  +  sin( alpha_V ) * cos(eps_V) * ri(1,0)     + sin( eps_V ) * ri(2,0);
    rY = -sin( alpha_V ) * ri(0,0)               +  cos( alpha_V ) * ri(1,0);
    rZ = -cos( alpha_V ) * sin( eps_V )*ri(0,0)  -  sin( alpha_V ) * sin( eps_V ) * ri(1,0)   + cos( eps_V ) * ri(2,0);

    // integration distance 
    pk0 = sqrt( pow( rX, 2.0 ) / pow( a, 2.0 ) + pow( rY, 2.0 ) / pow( b, 2.0 ) + pow( rZ, 2.0 ) / pow( c, 2.0 ) ) * k0 ;

    // factor of integral
    factor = 4.0* pow( M_PI, 2.0 ) / sin( az_el_i(1) ) * 0.033 * 3/5 * pow( k0, (-2.0/3.0) ) / ( pk0 * sqrt(M_PI) * std::tgamma(5.0/6.0)) * Cn_coeff;	// [m]

    // helping values
    h1 = 2.0*M_PI / (sqrt(3.0) * std::tgamma(2.0/3.0));   		// c ~ 2.6789 [-]
    h2 = 9 * pow( 2.0, (1.0/3.0) ) * std::tgamma(2.0/3.0) / 10;   	// d  ~1.5355 [-]
    h3 = pk0 / sin( az_el_i(1) );               			//            [1/m] 
    z =  h3 * _turb_data.h;                          		//            [-] 

    // calculate hypergeometric function (F23, F12)
    ivg::Matrix input_a(1,2); input_a(0,0) = 0.5; input_a(0,1) = 1.0;
    ivg::Matrix input_b(1,3); input_b(0,0) = 2.0/3.0; input_b(0,1) = 3.0/2.0; input_b(0,2) = 2.0 ;
    double input_z = 1.0 / 4.0 * pow( z, 2.0 );
    double F23 = hypergeom( input_a, input_b, input_z ) ;

    ivg::Matrix input_a2(1,1); input_a2(0,0) = 5.0 / 6.0;
    ivg::Matrix input_b2(1,2); input_b2(0,0) = 11.0 / 6.0; input_b2(0,1) = 7.0 / 3.0; 
    double input_z2 = 1.0 / 4.0 * pow( z, 2.0 );
    double F12 = hypergeom( input_a2, input_b2, input_z2 ) ;

    // closed solution of the integral				
    return factor / 8.0 * pow( z, 2.0 ) / h3* ( 4.0 * h1 * F23 - 3.0 * h2 * pow( z, (2.0/3.0) ) * F12 ); 
}


// calculate hypergeometric function
// .....................................................................
double Turbulence::hypergeom( ivg::Matrix a, ivg::Matrix b, double z )
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ double Turbulence::hypergeom( ivg::Matrix, ivg::Matrix, double )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
	
    double sum = 0.0 ;
    ivg::Matrix a2 = a.gamma_fct();
    ivg::Matrix b2 = b.gamma_fct();
    double factor = ( ( b2.transpose() ).prod_col() / ( a2.transpose() ).prod_col() )(0,0) ;

    ivg::Matrix ak, bk, K;
    ivg::Matrix ones1( a.size(1), a.size(2), 1.0 );
    ivg::Matrix ones2( b.size(1), b.size(2), 1.0 );

    for( int k = 0; k <= 20; k++ )
    {
        ak = (a+ones1*(double)k).gamma_fct();
        bk = (b+ones2*(double)k).gamma_fct();	

        sum += ( ( ak.transpose() ).prod_col() / ( bk.transpose() ).prod_col() )(0,0) * pow( z, (double)k ) / boost::math::factorial<double>(k) ;
    }

    double F = factor * sum ;

#if DEBUG_TROP >=2
   cerr << "--- double Turbulence::hypergeom( ivg::Matrix, ivg::Matrix, double )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  
	
    return F ;
}



// ****************************************************************************
// ******************* (3) N I E L S O N  et al   M O D E L *******************
// ****************************************************************************
//// .....................................................................
//ivg::Matrix Turbulence::calc_vcm_onsala_model ( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, double dh, double dh_seg, ivg::Matrix & CORR ) 
//// .....................................................................
//{
//#if DEBUG_TROP >=1
//   cerr << "+++ ivg::Matrix Turbulence::calc_vcm_onsala_model( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, double, ivg::Matrix & )" << endl; 
//   tictoc tim;
//   tim.tic();
//#endif  
//
//    // initialization
//    Matrix VCM( nobs, nobs, 0.0 );
//    CORR.resize( nobs, nobs, 0.0 );
//
//    double x, y, z, dz1, dz2;
//    Matrix ri_z1( 3,1,0.0 );
//    Matrix r0_z2( 3,1,0.0 );
//    Matrix rj_z1( 3,1,0.0 );
//    Matrix rj_z2( 3,1,0.0 );
//    Matrix r0_z1( 3,1,0.0 );
//    Matrix v_dti0( 3,1,0.0 );
//    Matrix v_dtj0( 3,1,0.0 );
//    Matrix v_dtij( 3,1,0.0 );
//
//    double tau_A, tau_B, tau;
//
//    double term1, term2, term3, term4, aux1, aux2, aux3, aux4 ;
//
//    // integration constants
//    Matrix height( 0.0, dh, _turb_data.h, 1); 
//    height( height.size(1)-1, 0 ) = _turb_data.h ;
//
//    double Cn_coeff = 0.5 * pow( _turb_data.Cn, 2.0 ) * 1.0e6 * pow( dh, 2.0 );
//    // double Cn_coeff = pow( _turb_data.Cn, 2.0 ) * 1.0e6 ;
//
//    int counter1 = 0;
//    int counter2 = 0;
//    int J = 0;
//    int I = 0;
//
//    std::string sta1, sta2, sta3, sta4;
//
//    ivg::Matrix az_el_i1;
//    ivg::Matrix az_el_i2;
//    ivg::Matrix az_el_j1;
//    ivg::Matrix az_el_j2;
//    
//    int ref = 0;
//    while( scans.at(ref).get_nobs() == 0 )
//        ref++;
//    
//    ivg::Matrix az_el_01 = scans.at(ref).get_obs_ptr(0)->get_az_el(1);
//    ivg::Matrix az_el_02 = scans.at(ref).get_obs_ptr(0)->get_az_el(2);
//
//    ivg::Obs* obs_ptr1;
//    ivg::Obs* obs_ptr2;
//    ivg::Matrix v(3,1,1.0);
//    v.set_sub( 0,0, _turb_data.vel.transpose() );
//   
//    // start calculation
//    // fill upper right part of the variance covariance matrix
//    for( int k = 0; k <= scans.size()-1; k++ )                    // loop over scans k
//    {
//       tau_A = ( scans.at(k).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0  ;
//
//       for( int i = 0; i <= scans.at(k).get_nobs()-1; i++ )       // loop over obs i
//       {
//          counter2 = 0;
//          I = i+counter1;
////              J = I; 
//
//          // get pointer to current observation i
//          obs_ptr1 = scans.at(k).get_obs_ptr(i);
//
//          // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
//          az_el_i1 = obs_ptr1->get_az_el(1);
//          az_el_i2 = obs_ptr1->get_az_el(2);
//          obs_ptr1->get_station_names( sta1,sta2 );
//
//          for( int l = 0; l <= scans.size()-1; l++ )              // loop over scans l
//          {   
//            tau_B = ( scans.at(l).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0;
//
//            // calculate time differences [hours]
//            tau = abs( scans.at(k).get_epoch().get_double_mjd() - scans.at(l).get_epoch().get_double_mjd() ) * 24.0;
//
//            for( int j = 0; j <= scans.at(l).get_nobs()-1; j++ ) // loop over obs j
//            { 
//                J = j+counter2;
//
//                if( VCM(J,I) == 0.0 ) // only fill upper right triangle to reduce computational cost
//                {
//                    // get pointer to current observation j
//                    obs_ptr2 = scans.at(l).get_obs_ptr(j);
//
//                    // get azimuth and elevation as well as station name of station 1 and 2 for observation j  
//                    az_el_j1 = obs_ptr2->get_az_el(1);
//                    az_el_j2 = obs_ptr2->get_az_el(2);
//                    obs_ptr2->get_station_names( sta3,sta4 );
//
//                    // VARIANCES
//                    if ( I == J )	
//                    {				
//                       VCM(I,J) = 1e2 * _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i1, az_el_j1, az_el_01, dh )
//                                + _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i2, az_el_j2, az_el_02, dh );
//                    }
//                    // COVARIANCES
//                    else 	
//                    {
//                       // (*) check which observations are correlated (in time); 
//                       //     spatial correlation is considered by the separation distance (= distance between two rays)  
//                       // (1) calcualte separation distance between two rays
//                       // (2) calculate covariance between observations i and j 
//                       if( sta1 == sta3 && sta2 == sta4 )
//                          VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i1, az_el_j1, az_el_01, dh )
//                                   + _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i2, az_el_j2, az_el_02, dh );
//
//                       else if( sta1 == sta4 && sta2 == sta3 )
//                          VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i1, az_el_j2, az_el_02, dh )
//                                   + _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i2, az_el_j1, az_el_02, dh );
//
//                       else if( sta1 == sta3 )
//                          VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i1, az_el_j1, az_el_01, dh );
//
//                       else if( sta1 == sta4 )
//                          VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i1, az_el_j2, az_el_02, dh );
//
//                       else if( sta2 == sta3 )
//                          VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i2, az_el_j1, az_el_01, dh );
//
//                       else if( sta2 == sta4 )
//                          VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i2, az_el_j2, az_el_02, dh );
//
//                       else
//                            VCM(I,J) = 0.0;
//                    }
//                }
//                J++;
//
//            } // #j
//            counter2 += scans.at(l).get_nobs();
//
//          } // #l
//       } // #i
//
//       counter1 += scans.at(k).get_nobs();
//    } // #k
//
//    // fill lower left part of the covariance matrix
//    ivg::Matrix D = VCM.diag();
//    VCM = VCM + VCM.transpose() - D.diag(); 
//
//    // calculate correlation matrix
//    ivg::Matrix ones( nobs, 1, 1.0 );
//    ivg::Matrix tmp1 = ( ones.div_elem( ( VCM.diag() ).sqrt() ) ).diag();
//    CORR = tmp1 * VCM * tmp1 ;
//
//#if DEBUG_TROP >=1
//    cerr << "--- ivg::Matrix Turbulence::calc_vcm_onsala_model( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, double, ivg::Matrix & )"
//         << " : " << tim.toc() << " s " << endl; 
//#endif  
//	
//    return VCM;
//}
//
//
//// .....................................................................
//double Turbulence::_calc_onsala_covariance( double tau_A0, double tau_B0, double tau_AB, ivg::Matrix az_el_i, ivg::Matrix az_el_j, ivg::Matrix az_el_0, double dh )
//// .....................................................................
//{
//#if DEBUG_TROP >=2
//   cerr << "+++ double Turbulence::_calc_onsala_covariance( double, double, ivg::Matrix, ivg::Matrix, double ) )" << endl; 
//   tictoc tim;
//   tim.tic();
//#endif  
//
//	double x, y, z, dz1, dz2;
//	Matrix ri_z1( 3,1,0.0 );
//	Matrix r0_z2( 3,1,0.0 );
//	Matrix rj_z1( 3,1,0.0 );
//	Matrix rj_z2( 3,1,0.0 );
//	Matrix r0_z1( 3,1,0.0 );
//	Matrix v_dti0( 3,1,0.0 );
//	Matrix v_dtj0( 3,1,0.0 );
//	Matrix v_dtij( 3,1,0.0 );
//	
//	double covar = 0.0;
//	double term1, term2, term3, term4, aux1, aux2, aux3, aux4 ;
//	
//	// integration constants
//	Matrix height( 0.0, dh, _turb_data.h, 1); 
//	height( height.size(1)-1, 0 ) = _turb_data.h ;
//	
//	double Cn_coeff = 0.5 * pow( _turb_data.Cn, 2.0 ) * 1.0e6 * pow( dh, 2.0 );
//	// double Cn_coeff = pow( _turb_data.Cn, 2.0 ) * 1.0e6 ;
//
//        ivg::Matrix v(3,1,1.0);
//        v.set_sub( 0,0, _turb_data.vel.transpose() );
//
//	// start numerical integration		
//	for ( int z2 = 0; z2 < height.size(1); z2++ )
//	{
//		for ( int z1 = 0; z1 < height.size(1); z1++ )
//		{			
//			dz1 = height(z1,0) ; // / sin( _data( j,1 ) );
//			dz2 = height(z2,0) ; // / sin( _data( j,1 ) );
//
//			// calculate directions
//			// term ri(z)
//			x = dz1 * sin( az_el_i(0) ) / tan( az_el_i(1) );
//			y = dz1 * cos( az_el_i(0) ) / tan( az_el_i(1) );
//			ri_z1(0,0) = x ; ri_z1(1,0) = y ; ri_z1(2,0) = dz1 ;
//			
//			// term r0(z')
//			x = dz2 * sin( az_el_0(0) ) / tan( az_el_0(1) );
//			y = dz2 * cos( az_el_0(0) ) / tan( az_el_0(1) );
//			r0_z2(0,0) = x ; r0_z2(1,0) = y ; r0_z2(2,0) = dz2 ;
//			
//			// term rj(z)
//			x = dz1 * sin( az_el_j(0) ) / tan( az_el_j(0) );
//			y = dz1 * cos( az_el_j(0) ) / tan( az_el_j(0) );
//			rj_z1(0,0) = x ; rj_z1(1,0) = y ; rj_z1(2,0) = dz1 ;
//			
//			// term rj(z')
//			x = dz2 * sin( az_el_j(0) ) / tan( az_el_j(1) );
//			y = dz2 * cos( az_el_j(0) ) / tan( az_el_j(1) );
//			rj_z2(0,0) = x ; rj_z2(1,0) = y ; rj_z2(2,0) = dz2 ;
//
//			// term r0(z)
//			x = dz1 * sin( az_el_0(0) ) / tan( az_el_0(1) );
//			y = dz1 * cos( az_el_0(0) ) / tan( az_el_0(1) );
//			r0_z1(0,0) = x ; r0_z1(1,0) = y ; r0_z1(2,0) = dz1 ;							
// 
//			// considering Taylor's frozen flow hypothesis
//			v_dti0 = v * tau_A0 ;
//			v_dtj0 = v * tau_B0 ;
//			v_dtij = v * tau_AB ;
//
//			// calculate the different terms of the structure function
//			aux1 = ( ( ( ri_z1 - r0_z2 + v_dti0 ).transpose() * ( ri_z1 - r0_z2 + v_dti0 ) ).sqrt() )(0,0);
//			aux1 = pow( aux1, 2.0/3.0 );
//			term1 =  Cn_coeff * ( aux1 / ( 1.0 + aux1 / pow( _L, 2.0/3.0 ) ) );
//
//			aux2 = ( ( ( rj_z1 - r0_z2 + v_dtj0 ).transpose() * ( rj_z1 - r0_z2 + v_dtj0 ) ).sqrt() )(0,0);
//			aux2 = pow( aux2, 2.0/3.0 );
//			term2 =  Cn_coeff * ( aux2 / ( 1.0 + aux2 / pow( _L, 2.0/3.0 ) ) );
//			
//			aux3 = ( ( ( ri_z1 - rj_z2 + v_dtij ).transpose() * ( ri_z1 - rj_z2 + v_dtij ) ).sqrt() )(0,0);
//			aux3 = pow( aux3, 2.0/3.0 );
//			term3 =  Cn_coeff * ( aux3 / ( 1.0 + aux3 / pow( _L, 2.0/3.0 ) ) );
//			
//			aux4 = ( ( ( r0_z1 - r0_z2 ).transpose() * ( r0_z1 - r0_z2 ) ).sqrt() )(0,0);
//			aux4 = pow( aux4, 2.0/3.0 );
//			term4 =  Cn_coeff * ( aux4 / ( 1.0 + aux4 / pow( _L, 2.0/3.0 ) ) );										
//
//
//			covar += term1 + term2 - term3 - term4 ;
//		}
//	}
//
//
//#if DEBUG_TROP >=2
//   cerr << "+++ double Turbulence::_calc_onsala_covariance( double, double, ivg::Matrix, ivg::Matrix, double ) )"
//        << " : " << tim.toc() << " s " << endl; 
//#endif  
//
//   return covar;
//}



// .....................................................................
ivg::Matrix Turbulence::calc_vcm_onsala_model ( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, double dh, double dh_seg, ivg::Matrix & CORR ) 
// .....................................................................
{
#if DEBUG_TROP >=1
   cerr << "+++ ivg::Matrix Turbulence::calc_vcm_onsala_model( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, double, ivg::Matrix & )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
   
    // initialization
    Matrix VC( nobs, nobs, 0.0 );
    Matrix VCM( nobs, nobs, 0.0 );
    CORR.resize( nobs, nobs, 0.0 );

    double tau_A, tau_B, tau;


    ivg::Matrix az_el_i0(2,1,0.0);
    az_el_i0(0,0) = 0.0; az_el_i0(1,0) = 0.5 *M_PI;
    ivg::Matrix az_el_j0(2,1,0.0);
    az_el_j0(0,0) = 0.0; az_el_j0(1,0) = 0.5 *M_PI;
    
    ivg::Obs* obs_ptr1;
    ivg::Obs* obs_ptr2;
    ivg::Matrix v(3,1,1.0);

    std::string sta1, sta2, sta3, sta4; 
    
    std::map< std::string, ivg::Matrix > sta_vcm;

    std::vector<Analysis_station>::iterator iter;
    for( iter = trf->begin(); iter < trf->end(); iter++ )
    {
        
        std::string sta_name = iter->get_name( ivg::staname::ivs_name );
        
        _turb_data = _turb_sta[sta_name];
        
        v.set_sub( 0, 0, _turb_data.vel );
                
        // integration constants
        Matrix height( 0.0, dh, _turb_data.h, 1); 
        height( height.size(1)-1, 0 ) = _turb_data.h ;
        
        ivg::Matrix tmp(nobs,nobs,0.0);
        VCM = tmp;

        int counter1 = 0;
        int counter2 = 0;
        int J = 0;
        int I = 0;
        
        ivg::Matrix az_el_i;
        ivg::Matrix az_el_j;          
        
        // start calculation
        // fill upper right part of the variance covariance matrix
        for( int k = 0; k <= scans.size()-1; k++ )                    // loop over scans k
        {
           tau_A = ( scans.at(k).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0  ;

           for( int i = 0; i <= scans.at(k).get_nobs()-1; i++ )       // loop over obs i
           {
                counter2 = 0;
                I = i+counter1;

                // get pointer to current observation i
                obs_ptr1 = scans.at(k).get_obs_ptr(i);
                obs_ptr1->get_station_names( sta1,sta2 );
                
                // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
                if( sta_name == sta1 )
                    az_el_i = obs_ptr1->get_az_el(1);
                else if( sta_name == sta2 )
                    az_el_i = obs_ptr1->get_az_el(2);
                
                
              for( int l = 0; l <= scans.size()-1; l++ )              // loop over scans l
              {   
                tau_B = ( scans.at(l).get_epoch().get_double_mjd() - scans.at(0).get_epoch().get_double_mjd() ) * 24.0;

                // calculate time differences [hours]
                tau = abs( scans.at(k).get_epoch().get_double_mjd() - scans.at(l).get_epoch().get_double_mjd() ) * 24.0;

                for( int j = 0; j <= scans.at(l).get_nobs()-1; j++ ) // loop over obs j
                { 
                    J = j+counter2;

                    if( VCM(J,I) == 0.0 ) // only fill upper right triangle to reduce computational cost
                    {
                        // get pointer to current observation j
                        obs_ptr2 = scans.at(l).get_obs_ptr(j);
                        obs_ptr2->get_station_names( sta3,sta4 );
                        
                        // only fill VCM if the a station joins the current observation
                        if( sta_name == sta1 || sta_name == sta2 )
                        {                        
                            // VARIANCES
                            if ( I == J )	
                            {				
                               VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i, az_el_i, az_el_i0, dh );
                            }
                            // COVARIANCES
                            else 	
                            {
                                // (*) check which observations are correlated (in time); 
                                //     spatial correlation is considered by the separation distance (= distance between two rays)  
                                // (1) calcualte separation distance between two rays
                                // (2) calculate covariance between observations i and j 
                                if( sta_name == sta3 )
                                {
                                   az_el_j = obs_ptr2->get_az_el(1);
                                   VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i, az_el_j, az_el_i0, dh );
                                }
                                else if( sta_name == sta4 )
                                {
                                   az_el_j = obs_ptr2->get_az_el(2);
                                   VCM(I,J) = _calc_onsala_covariance( tau_A, tau_B, tau, az_el_i, az_el_j, az_el_j0, dh );
                                }
                            }
                        }
                    }
                    J++;

                } // #j
                counter2 += scans.at(l).get_nobs();

              } // #l
           } // #i

           counter1 += scans.at(k).get_nobs();
        } // #k

        // fill lower left part of the covariance matrix
        ivg::Matrix D = VCM.diag();
        VCM = VCM + VCM.transpose() - D.diag(); 

//        VCM.save_bin( "/home/halsig/ascot/output/turbOnsala_"+sta_name);
        
        VC += VCM;
        sta_vcm[ sta_name ] = VCM;         
    } // stations    
    
    // calculate correlation matrix
    ivg::Matrix ones( nobs, 1, 1.0 );
    ivg::Matrix tmp1 = ( ones.div_elem( ( VC.diag() ).sqrt() ) ).diag();
    CORR = tmp1 * VC * tmp1 ;

#if DEBUG_TROP >=1
    cerr << "--- ivg::Matrix Turbulence::calc_vcm_onsala_model( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, double, ivg::Matrix & )"
         << " : " << tim.toc() << " s " << endl; 
#endif  
	
    return VC;
}

// .....................................................................
ivg::Matrix Turbulence::calc_vcm_onsala_model( ivg::Analysis_station* sta , ivg::Matrix azel, std::vector<ivg::Date> epochs, double dh ){
// .....................................................................     
    
    int npos = azel.rows();
    int nepochs = epochs.size();
    int nobs = npos * nepochs;
    
    log<INFO>("*** Simulate troposphere #") % nobs;

    // initialization
    Matrix VCM( nobs, nobs, 0.0 );

    double tau_i0(0.0),  tau_j0(0.0),  tau_ij(0.0);
    
    _turb_data = _turb_sta[sta->get_name(ivg::staname::ivs_name) ];

    ivg::Matrix az_el_0(2,1,0.0);
    az_el_0(0,0) = 0.0; az_el_0(1,0) = 0.5 *M_PI;
    
    
    std::cout << " h: " << _turb_data.h << "    Cn: " << _turb_data.Cn << "   L: " << _L << std::endl;
     
//    std::cerr << "calc vcm" << std::string(9,' ');
    unsigned I = 0;
    for( unsigned epoch_idx_1 = 0; epoch_idx_1 < nepochs; ++epoch_idx_1 ){
        tau_i0 =  ( epochs[epoch_idx_1].get_double_mjd() - epochs[0].get_double_mjd() ) * 24.0;
        
        for( int pos_idx_1 = 0; pos_idx_1 < azel.rows(); ++pos_idx_1 ){

            ivg::Matrix azel_1 = azel.get_sub (pos_idx_1, 0, pos_idx_1, 1);
            
//            std::cerr << std::string(9,'\b') << setw(4) << I+1 << "/" << setw(4) << nobs;
            
            for( unsigned epoch_idx_2 = 0; epoch_idx_2 < nepochs; ++epoch_idx_2 ){
                tau_j0 =  ( epochs[epoch_idx_2].get_double_mjd() - epochs[0].get_double_mjd() ) * 24.0;
                tau_ij =  ( epochs[epoch_idx_1].get_double_mjd() - epochs[epoch_idx_2].get_double_mjd() ) * 24.0;
                
//                #pragma omp parallel for num_threads(2)
                for( int pos_idx_2 = 0; pos_idx_2 < azel.rows(); ++pos_idx_2 ){
                    int J = epoch_idx_2*npos + pos_idx_2;
                                        
                    if(J>=I ){
                        ivg::Matrix azel_2 =   azel.get_sub (pos_idx_2, 0, pos_idx_2, 1);
                                
//                        std::cout << string(30,'-') << std::endl;
//                        std::cout << "i " << I << "    j " << J  << "    t1 "  <<  epoch_idx_1
//                                  << "   t2 " << epoch_idx_2  << "    p1 " << pos_idx_1 << "   p2 " << pos_idx_2
//                                  << "   taui0 " <<  tau_i0 << "   tauj0 "  << tau_j0 << "   tauij " << tau_ij << std::endl;
                        
                        VCM(I,J) = _calc_onsala_covariance( tau_i0, tau_j0, tau_ij, azel_1, azel_2, az_el_0, dh );
                        if(I!=J) {
                            VCM(J,I) = VCM(I, J);
                        }

                    }
                }
            }
            
            I++;
        }
    }
     
    
    // calculate correlation matrix
    ivg::Matrix ones( nobs, 1, 1.0 );
    ivg::Matrix tmp1 = ( ones.div_elem( ( VCM.diag() ).sqrt() ) ).diag();
    ivg::Matrix CORR = tmp1 * VCM * tmp1 ;
    
    
    
    std::cout << std::endl;
    VCM.save_bin( "/home/corbin/VCM.dat");
    CORR.save_bin( "/home/corbin/COR.dat");
    return VCM;
}

// .....................................................................
double Turbulence::_calc_onsala_covariance( double tau_i0, double tau_j0, double tau_ij, ivg::Matrix az_el_i, ivg::Matrix az_el_j, ivg::Matrix az_el_0, double dh )
// .....................................................................
{
    
    // tropospheric parameters
  
    const double h = _turb_data.h; // height in m
    const double Cn = _turb_data.Cn; // refractive index structure constant
    
    ivg::Matrix v(3, 1, 0.0); // wind vector in m/s
    v(0) = _turb_data.v;
        
    
    // directions with z = 1
    const ivg::Matrix r_i_z = _azel2r( az_el_i );
    const ivg::Matrix r_0_z = _azel2r( az_el_0 );
    const ivg::Matrix r_j_z = _azel2r( az_el_j );
    
//    r_i_z.show();
//    r_0_z.show();
//    r_j_z.show();    
    
    // heights in m
    dh = 200.0; // increment for numerical integration
    
    Matrix height( dh, dh, h, 1); 
    height( height.size(1)-1, 0 ) = h ;

    
    double Cn_coeff = 0.5 * pow( Cn, 2.0 ) * 1.0e6;
    
    // declarations
    ivg::Matrix r_i_z1(3, 1);
    
    ivg::Matrix r_j_z1(3, 1);
    ivg::Matrix r_j_z2(3, 1);
    
    ivg::Matrix r_0_z1(3, 1);
    ivg::Matrix r_0_z2(3, 1);
    
    const ivg::Matrix rv_i_j = v * tau_ij*3600.0;
    const ivg::Matrix rv_i_0 = v * tau_i0*3600.0;
    const ivg::Matrix rv_j_0 = v * tau_j0*3600.0;
    const ivg::Matrix rv_0_0(3, 1, 0.0);

    // numerical integration
    double sum = 0.0;
    for( const double z2 : *height.data_vec_ptr() ){
        
        r_j_z2 = r_j_z*z2;
        r_0_z2 = r_0_z*z2;
        
        for( const double z1 : *height.data_vec_ptr() ){

            r_i_z1 = r_i_z*z1;
            r_j_z1 = r_j_z*z1;
            r_0_z1 = r_0_z*z1;
            
            double term1 =  _calc_onsala_structure_function( r_i_z1, r_0_z2, rv_i_0, _L );
            double term2 =  _calc_onsala_structure_function( r_j_z1, r_0_z2, rv_j_0, _L );
            double term3 =  _calc_onsala_structure_function( r_i_z1, r_j_z2, rv_i_j, _L );
            double term4 =  _calc_onsala_structure_function( r_0_z1, r_0_z2, rv_0_0, _L );

            sum += ( term1 + term2 - term3 - term4 );
            
//            std::cout << string(20,'-') << std::endl;
//            std::cout << "z1: " << z1 << "   z2: " << z2 << std::endl;            
//            std::cout << "r_i_z1"<< std::endl;
//            r_i_z1.show();
//            std::cout << "r_j_z1"<< std::endl;
//            r_j_z1.show();
//            std::cout << "r_j_z2"<< std::endl;
//            r_j_z2.show();
//            std::cout << "r_0_z1"<< std::endl;
//            r_0_z1.show();
//            std::cout << "r_0_z2"<< std::endl;
//            r_0_z2.show();
//            std::cout << "term1 " << term1 << "   term2 " << term2  << "   term3 " << term3  << "   term4 " << term4 << "   su m" << sum << std::endl;
            
        }
    }
         
    // multiply with surface area of 'the top of pillars' (dh*dh) and with Cn term;
    sum *= pow( dh, 2.0 ) * Cn_coeff;
    
    return sum;
    
}



// ****************************************************************************
// ******************* (3b) Boehm version of onsala model  ********************
// ****************************************************************************
// .....................................................................
ivg::Matrix Turbulence::calc_vcm_vienna_model ( ivg::Trf * trf, std::vector<ivg::Scan> scans, int nobs, double dh, double dh_seg, ivg::Matrix & CORR ) 
// .....................................................................
{
#if DEBUG_TROP >=1
    cerr << "+++ ivg::Matrix Turbulence::calc_vcm_vienna_model( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, double, ivg::Matrix & )" << endl; 
    tictoc tim;
    tim.tic();
#endif  

    // initialization
    ivg::Matrix VC( nobs, nobs, 0.0 ); 
    ivg::Matrix VCM( nobs, nobs, 0.0 );

   
  
    ivg::Obs* obs_ptr1;
    std::string sta1, sta2; 
    
    std::map< std::string, ivg::Matrix > sta_vcm;

    std::vector<Analysis_station>::iterator iter;
    for( iter = trf->begin(); iter < trf->end(); iter++ )
    {
        ivg::Matrix tmp(nobs,nobs,0.0);
        VCM = tmp;

        int counter1 = 0;
        int counter2 = 0;
        ivg::Matrix az_el_i;
        std::string sta_name = iter->get_name( ivg::staname::ivs_name );
        

        std::vector<double> el;
        std::vector<double> az;
        std::vector<double> mjd;
        
       
        for( int k = 0; k <= scans.size()-1; k++ )                    // loop over scans k
        {
            for( int i = 0; i <= scans.at(k).get_nobs()-1; i++ )       // loop over obs i
            {
                counter2 = 0;

                // get pointer to current observation i
                obs_ptr1 = scans.at(k).get_obs_ptr(i);
                obs_ptr1->get_station_names( sta1,sta2 );
                
                // get azimuth and elevation as well as station name of station 1 and 2 for observation i  
                if( sta_name == sta1 || sta_name == sta2 )
                {
                    if( sta_name == sta1 )
                        az_el_i = obs_ptr1->get_az_el(1);
                    else if( sta_name == sta2 )
                        az_el_i = obs_ptr1->get_az_el(2);

                    az.push_back(az_el_i(0));
                    el.push_back(az_el_i(1));
                    mjd.push_back(obs_ptr1->get_epoch().get_double_mjd());
                }
            } // #i
           counter1 += scans.at(k).get_nobs();
        } // #k  
        
        double Cn = _turb_sta[iter->get_name(ivg::staname::ivs_name)].Cn;
        double h = _turb_sta[iter->get_name(ivg::staname::ivs_name)].h;
        ivg::Matrix vel = _turb_sta[iter->get_name(ivg::staname::ivs_name)].vel;
        
        log<INFO>("*** calc vcm for station " ) % iter->get_name(ivg::staname::ivs_name) % " with:   Cn " % Cn % "   h " % h % "   vx " % vel(0) % "   vy " % vel(1); 
        
        // start calculation
        // fill upper right part of the variance covariance matrix
        VCM = calc_vcm_vienna_model ( el, az, mjd, dh, dh_seg,  Cn,  h,  vel );
        

	// fill lower left part of the covariance matrix
	ivg::Matrix D = VCM.diag();
	VCM = VCM + VCM.transpose() - D.diag();         
        
//        VCM.save_bin( "/home/halsig/ascot/output/turbVienna_"+sta_name);
        
        VC += VCM;
        sta_vcm[ sta_name ] = VCM;        
    } // stations    
    
    // calculate correlation matrix
    ivg::Matrix ones( nobs, 1, 1.0 );
    ivg::Matrix tmp1 = ( ones.div_elem( ( VC.diag() ).sqrt() ) ).diag();
    CORR = tmp1 * VC * tmp1 ;
    
#if DEBUG_TROP >=1
   cerr << "--- ivg::Matrix Turbulence::calc_vcm_vienna_model( ivg::Trf *, std::vector<ivg::Scan>, int, ivg::Matrix, double, double, ivg::Matrix & )"
        << " : " << tim.toc() << " s " << endl; 
#endif  
	
    return VCM ;
}

// .....................................................................
ivg::Matrix Turbulence::calc_vcm_vienna_model ( ivg::Analysis_station* sta , ivg::Matrix azel, std::vector<ivg::Date> epochs, double dh )
// .....................................................................
{
    std::string sta_name = sta->get_name( ivg::staname::ivs_name );
    
    int npos = azel.rows();
    int nepochs = epochs.size();

    
    double Cn = _turb_sta[sta_name].Cn;
    double h = _turb_sta[sta_name].h;
    ivg::Matrix vel = _turb_sta[sta_name].vel;
    
    std::vector<double> el;
    std::vector<double> az;
    std::vector<double> mjd;
    
    for(unsigned ec = 0; ec < nepochs; ec++){
        double current_mjd =epochs[ec].get_double_mjd();
        for(unsigned pc = 0; pc < npos; pc++){
            mjd.push_back(current_mjd);
            az.push_back( azel(pc,0) );
            el.push_back( azel(pc,1) );
        }
    }
  
    log<INFO>("*** calc vcm for station " ) % sta_name % " with:   Cn " % Cn % "   h " % h % "   vx " % vel(0) % "   vy " % vel(1); 
    ivg::Matrix VCM = calc_vcm_vienna_model ( el, az, mjd, dh, 24.0,  Cn,  h,  vel );
    
    
    // fill lower left part of the covariance matrix
    ivg::Matrix D = VCM.diag();
    VCM = VCM + VCM.transpose() - D.diag();    
    
//    // calculate correlation matrix
//    ivg::Matrix ones( npos*nepochs, 1, 1.0 );
//    ivg::Matrix tmp1 = ( ones.div_elem( ( VCM.diag() ).sqrt() ) ).diag();
//    ivg::Matrix CORR = tmp1 * VCM * tmp1 ;
//    
//    VCM.save_bin( "/home/corbin/VCM.dat");
//    CORR.save_bin( "/home/corbin/COR.dat");
    
    return VCM;
  
}

// .....................................................................
ivg::Matrix Turbulence::calc_vcm_vienna_model ( const std::vector<double>& el, const std::vector<double>& az, const std::vector<double>& mjd,
                                                double dh, double dh_seg, double Cn, double h, ivg::Matrix vel ){
// .....................................................................
    
    
    unsigned nobs = el.size();
    ivg::Matrix VCM( nobs, nobs, 0.0 );
    
    // create reference direction r0
    // >>> discard first observation ( see Wresnik (2009); Boehm et al. (2007) )
    double az1 = 0.0;
    double el1 = M_PI / 2.0;
    ivg::Matrix R0( 1, 3, 0.0 );
    R0( 0 ) = cos( az1 ) / tan( el1 );
    R0( 1 ) = sin( az1 ) / tan( el1 );
    R0( 2 ) = 1;
    
     // coefficients and settings
    double nhi = floor( h / dh );		// number of height increment

    double Cn_coeff = pow( Cn, 2.0 ) / 2.0 ;	// squared structure constant, multiplied by factor 1/PI
    Cn_coeff *= 1.0e6 * pow( dh, 2 ) ;

    int n = static_cast<int> ( pow( nhi + 1, 2.0 ) ) ;	// type casting: conversion from double to int    

    ivg::Matrix zs( 1, n, 0.0 );
    ivg::Matrix z ( 1, n, 0.0 );
    ivg::Matrix v ( 3, n, 0.0 );  

    v.set_sub( 0,0, vel );
    int k1 = -1;
    // define differentials dz and dzs (for numerical integration)        
    for( int i=0; i <= nhi; i++ )
    {         
        for( int r=0; r <= nhi; r++ )
        {
            k1++;
            zs(k1) = dh * double( i );		// type casting: int to double 
            z (k1) = dh * double( r );		// type casting: int to double

            v.set_sub( 0, k1, vel );
        }
    }
        
    ivg::Matrix t(mjd.size(),1,0.0);
    ivg::Matrix Ri(mjd.size(),3,0.0);
    ivg::Matrix ti( mjd.size(),1,0.0);

    for( int j=0; j< mjd.size(); ++j)
    {
        // t(j,0) = (mjd.at(j) - mjd.front()) * 24.0;
        t(j,0) = (mjd.at(j) - floor(mjd.at(j))) * 24.0;
        ti(j,0) = floor( t(j,0) * ( 1.0 / dh_seg ) ) + 1.0;

        // get directions and time tag from data-matrix
        Ri(j,0) = cos( az.at(j) ) / tan( el.at(j) );
        Ri(j,1) = sin( az.at(j) ) / tan( el.at(j) );
        Ri(j,2) = 1.0;             
    }
    
    // get the last segment which is used to define how many segments are needed
    int isegments = ti( mjd.size() - 1 );

    ivg::Matrix tn( 1, isegments );
    for( int i=0; i < isegments; i++ )
    {
        vector<int> k = ti.find_idx( double( i+1 ) );
        tn( i ) = k.size();
    }               

    // --- calculate different terms of variance-covariance matrix (term1x ... term4x) ---
    // + + + + + + + + + + + + +
    // + + + +  term 4 + + + + +
    // + + + + + + + + + + + + +

    // compute position vectors
    ivg::Matrix R0z  = R0.transpose() * z;
    ivg::Matrix R0zs = R0.transpose() * zs;	

    // compute separations
    ivg::Matrix tmp4 = ( R0z - R0zs )^2.0; 

    // get term4 (only reference positions R0; no consideration of wind vector)
    tmp4 = ( tmp4.sum_col() )^0.5;
    tmp4 =  tmp4^( 2.0/3.0 );
    ivg::Matrix term4 = tmp4 / ( ( tmp4 * 1.0 / ( pow( _L, 2.0/3.0 ) ) ) + 1.0 ) ; 

    // if only one segment
    if( isegments == 1 )
    {
        int num3 = tn( 0 ); 
        VCM( 0,0 ) = 0.0;

        // how many epochs are before the segment i1
        for( int i=0; i < num3; i++ )
        {
            // + + + + + + + + + + + + +
            // + + + +  term 1 + + + + +
            // + + + + + + + + + + + + +

            // compute position vectors
            ivg::Matrix RiSub = Ri.get_sub( i, 0, i, 2);
            ivg::Matrix Riz = RiSub.transpose() * z;

            // compute time differences in hours
            double dti0 = t(i);

            // compute separations
            ivg::Matrix dd1 = Riz - R0zs + v* dti0;
            
            // get term1 (difference between position Ri and reference positions R0)	
            dd1 = dd1^2;
            ivg::Matrix tmp1  = ( dd1.sum_col() )^(1.0/3.0);
            ivg::Matrix term1 = tmp1.div_elem( tmp1/ pow( _L, 2.0/3.0 ) + 1.0 );

            for( int r=i; r<num3; r++ )
            {				
                // + + + + + + + + + + + + +
                // + + + terms 2 and 3 + + +
                // + + + + + + + + + + + + +

                // compute time differences in hours
                double dtj0 = t( r );
                double dtij = t( r ) - t( i );

                // compute position vectors
                ivg::Matrix RjSub = Ri.get_sub( r, 0, r, 2);
                ivg::Matrix Rjz   = RjSub.transpose() * z;
                ivg::Matrix Rjzs  = RjSub.transpose() * zs;

                // compute separations
                ivg::Matrix dd2 = Rjz - R0zs + v *dtj0;
                ivg::Matrix dd3 = Riz - Rjzs - v *dtij;
                
                // get terms 2 (difference between position Rj and reference positions R0)
                // and 3 (difference between position Ri and positions Rj)
                dd2 = dd2^2;
                dd3 = dd3^2;
                ivg::Matrix tmp2 = ( dd2.sum_col() )^(1.0/3.0);
                ivg::Matrix tmp3 = ( dd3.sum_col() )^(1.0/3.0);

                ivg::Matrix term2 = tmp2.div_elem( tmp2/( pow( _L, 2.0/3.0 ) ) + 1.0 );
                ivg::Matrix term3 = tmp3.div_elem( tmp3/( pow( _L, 2.0/3.0 ) ) + 1.0 );

                // join all terms	                
                ivg::Matrix out = term1 + term2 - term3 - term4;
                out = out.transpose();
                ivg::Matrix result = out.sum_col();

                VCM( i,r ) = result( 0 ) * Cn_coeff;
            }
        }	
    }
    else
    {
        // first two segments
        int i1 = 0;
        int num1 = tn( i1 );
        int num2 = tn( i1 + 1 );
        int num3 = num1 + num2;

        VCM.resize( num3, num3, 0.0 );

        // how many epochs are before the segment i1
        int k = 0;

        for( int i=0; i < num3; i++ )
        {
            // + + + + + + + + + + + + +
            // + + + +  term 1 + + + + +
            // + + + + + + + + + + + + +

            // compute position vectors
            ivg::Matrix RiSub = Ri.get_sub( i+k, 0, i+k, 2);
            ivg::Matrix Riz = RiSub.transpose() * z;

            // compute time differences in hours
            double dti0 = t(i+k);

            // compute separations
            ivg::Matrix dd1 = Riz - R0zs + v* dti0;

            // get term1 (difference between position Ri and reference positions R0)	
            dd1 = dd1^2;
            ivg::Matrix tmp1  = ( dd1.sum_col() )^(1.0/3.0);
            ivg::Matrix term1 = tmp1.div_elem( tmp1/ pow( _L, 2.0/3.0 ) + 1.0 );

            for( int r=i; r<num3; r++ )
            {				
                // + + + + + + + + + + + + +
                // + + + terms 2 and 3 + + +
                // + + + + + + + + + + + + +

                // compute time differences in hours
                double dtj0 = t( r+k );
                double dtij = t( r+k ) - t( i+k );

                // compute position vectors
                ivg::Matrix RjSub = Ri.get_sub( r+k, 0, r+k, 2);
                ivg::Matrix Rjz   = RjSub.transpose() * z;
                ivg::Matrix Rjzs  = RjSub.transpose() * zs;

                // compute separations
                ivg::Matrix dd2 = Rjz - R0zs + v *dtj0;
                ivg::Matrix dd3 = Riz - Rjzs - v *dtij;

                // get terms 2 (difference between position Rj and reference positions R0)
                // and 3 (difference between position Ri and positions Rj)
                dd2 = dd2^2;
                dd3 = dd3^2;
                ivg::Matrix tmp2 = ( dd2.sum_col() )^(1.0/3.0);
                ivg::Matrix tmp3 = ( dd3.sum_col() )^(1.0/3.0);

                ivg::Matrix term2 = tmp2.div_elem( tmp2/( pow( _L, 2.0/3.0 ) ) + 1.0 );
                ivg::Matrix term3 = tmp3.div_elem( tmp3/( pow( _L, 2.0/3.0 ) ) + 1.0 );

                // join all terms	
                ivg::Matrix out = term1 + term2 - term3 - term4;
                out = out.transpose();
                ivg::Matrix result = out.sum_col();

                VCM( i,r ) = result( 0 ) * Cn_coeff;		
            }
        }	

        ivg::Matrix VCM11;
        if( VCM.length() != 0 )
        {
            VCM11 = VCM.get_sub( num1, num1, num3-1, num3-1 );
        }
        else
        { 
            VCM11 = VCM;
        }

        // loop over two adjacent two-hour blocks
        for( int i1 = 1; i1 < isegments-1; i1++ )
        {
            int num1 = tn( i1 );
            int num2 = tn( i1 + 1 );
            int num3 = num1 + num2;

            VCM.resize( num3, num3, 0.0 );

            // how many epochs are before the segment i1?
            int k = 0;
            for( int k1 = 0; k1 < i1; k1++ )
                    k += tn( k1 );
            VCM.set_sub( 0, 0, VCM11 );

            for( int i=0; i < num1; i++ )
            {
                // + + + + + + + + + + + + +
                // + + + +  term 1 + + + + +
                // + + + + + + + + + + + + +

                // compute position vectors
                ivg::Matrix RiSub = Ri.get_sub( i+k, 0, i+k, 2);
                ivg::Matrix Riz = RiSub.transpose() * z;

                // compute time differences in hours
                double dti0 = t(i+k);

                // compute separations
                ivg::Matrix dd1 = Riz - R0zs + v* dti0;

                // get term1 (difference between position Ri and reference positions R0)	
                dd1 = dd1^2;
                ivg::Matrix tmp1  = ( dd1.sum_col() )^(1.0/3.0);
                ivg::Matrix term1 = tmp1.div_elem( tmp1/ pow( _L, 2.0/3.0 ) + 1.0 );

                for( int r = num1; r < num3; r++ )
                {				
                    // + + + + + + + + + + + + +
                    // + + + terms 2 and 3 + + +
                    // + + + + + + + + + + + + +

                    // compute time differences in hours
                    double dtj0 = t( r+k );
                    double dtij = t( r+k ) - t( i+k );

                    // compute position vectors
                    ivg::Matrix RjSub = Ri.get_sub( r+k, 0, r+k, 2);
                    ivg::Matrix Rjz   = RjSub.transpose() * z;
                    ivg::Matrix Rjzs  = RjSub.transpose() * zs;

                    // compute separations
                    ivg::Matrix dd2 = Rjz - R0zs + v *dtj0;
                    ivg::Matrix dd3 = Riz - Rjzs - v *dtij;

                    // get terms 2 (difference between position Rj and reference positions R0)
                    // and 3 (difference between position Ri and positions Rj)
                    dd2 = dd2^2;
                    dd3 = dd3^2;
                    ivg::Matrix tmp2 = ( dd2.sum_col() )^(1.0/3.0);
                    ivg::Matrix tmp3 = ( dd3.sum_col() )^(1.0/3.0);

                    ivg::Matrix term2 = tmp2.div_elem( tmp2/( pow( _L, 2.0/3.0 ) ) + 1.0 );
                    ivg::Matrix term3 = tmp3.div_elem( tmp3/( pow( _L, 2.0/3.0 ) ) + 1.0 );

                    // join all terms	
                    ivg::Matrix out = term1 + term2 - term3 - term4;
                    out = out.transpose();
                    ivg::Matrix result = out.sum_col();

                    VCM( i,r ) = result( 0 ) * Cn_coeff;				
                }
            }
            for( int i = num1; i < num3; i++ )
            {
                // + + + + + + + + + + + + +
                // + + + +  term 1 + + + + +
                // + + + + + + + + + + + + +

                // compute position vectors
                ivg::Matrix RiSub = Ri.get_sub( i+k, 0, i+k, 2);
                ivg::Matrix Riz = RiSub.transpose() * z;

                // compute time differences in hours
                double dti0 = t(i+k);

                // compute separations
                ivg::Matrix dd1 = Riz - R0zs + v* dti0;

                // get term1 (difference between position Ri and reference positions R0)	
                dd1 = dd1^2;
                ivg::Matrix tmp1  = ( dd1.sum_col() )^(1.0/3.0);
                ivg::Matrix term1 = tmp1.div_elem( tmp1/ pow( _L, 2.0/3.0 ) + 1.0 );

                for( int r = num1; r < num3; r++ )
                {				
                    // + + + + + + + + + + + + +
                    // + + + terms 2 and 3 + + +
                    // + + + + + + + + + + + + +

                    // compute time differences in hours
                    double dtj0 = t( r+k );
                    double dtij = t( r+k ) - t( i+k );

                    // compute position vectors
                    ivg::Matrix RjSub = Ri.get_sub( r+k, 0, r+k, 2);
                    ivg::Matrix Rjz   = RjSub.transpose() * z;
                    ivg::Matrix Rjzs  = RjSub.transpose() * zs;

                    // compute separations
                    ivg::Matrix dd2 = Rjz - R0zs + v *dtj0;
                    ivg::Matrix dd3 = Riz - Rjzs - v *dtij;

                    // get terms 2 (difference between position Rj and reference positions R0)
                    // and 3 (difference between position Ri and positions Rj)
                    dd2 = dd2^2;
                    dd3 = dd3^2;
                    ivg::Matrix tmp2 = ( dd2.sum_col() )^(1.0/3.0);
                    ivg::Matrix tmp3 = ( dd3.sum_col() )^(1.0/3.0);

                    ivg::Matrix term2 = tmp2.div_elem( tmp2/( pow( _L, 2.0/3.0 ) ) + 1.0 );
                    ivg::Matrix term3 = tmp3.div_elem( tmp3/( pow( _L, 2.0/3.0 ) ) + 1.0 );

                    // join all terms	
                    ivg::Matrix out = term1 + term2 - term3 - term4;
                    out = out.transpose();
                    ivg::Matrix result = out.sum_col();

                    VCM( i,r ) = result( 0 ) * Cn_coeff;				
                }
            }

            if( VCM.length() != 0 )
            {
                VCM11 = VCM.get_sub( num1, num1, num3-1, num3-1 );
            }
            else
            { 
                VCM11 = VCM;
            }
        }
    }
    
    return VCM;
}
/*
// calculate structure function
// .....................................................................
ivg::Matrix Turbulence::calc_structure_function( int nobs ) 
// .....................................................................
{
#if DEBUG_TROP >=2
   cerr << "+++ ivg::Matrix Turbulence::calc_structure_function() )" << endl; 
   tictoc tim;
   tim.tic();
#endif  
	
	ivg::Matrix D( nobs, nobs, 0.0 ) ;
	ivg::Matrix tmp, Ri, Rj;
	ivg::Matrix ones( 1,1,1.0 );
	double dt;
	
	for( int i = 0; i < nobs; i++ ) 
	{
		for( int j = 0; j < nobs; j++ )
		{
			dt = _data( i,1 ) - _data( j,1 );
				
			// get directions from data-matrix
			Ri = ( _data.get_sub( i, 3, i, 4 ) ).transpose() ;
			Rj = ( _data.get_sub( j, 3, j, 4 ) ).transpose() ;
			
			// create auxiliary matrix
			tmp = ( ( Ri - Rj + _turb_data.vel.transpose()* dt ).norm() )^( 2.0/3.0 ) ;							
	
			// calculate structure function			
			D( i,j ) = ( tmp / ( ones + ( tmp / pow( _L, 2/3 ) ) ) * pow( _turb_data.Cn, 2 ) )( 0 );			
		}
	}

#if DEBUG_TROP >=2
   cerr << "--- ivg::Matrix Turbulence::calc_structure_function() )" 
        << " : " << tim.toc() << " s " << endl; 
#endif  
	
	return D ;
}

*/


} // # namespace ivg









