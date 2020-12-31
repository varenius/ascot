/*****************************************************************************
* wrapper-file for using IERS fortran routines from C++
*
* -TA- 2014-12-03
*****************************************************************************/
#ifndef IERSWRAP_H
#define IERSWRAP_H

#include "cmath"
#include "map"

/**
*
* @brief wrapper functions to access routines provided with IERS 2010 Conventions.
*        Original source code can be derived at 
*        http://maia.usno.navy.mil/conv2010/software.html
*
*/ 

namespace iers
{
   extern"C" 
   {
      /**
      *  \b Description: \n
      *        Subroutine for transformation from cartesian (x, y, z)to geographic
      *        coordinates which uses GRS 1980 ellipsoid parameters by default.
      *        Thus, a and f need not be defined.
      *  \param [in] [double *] Equatorial Radius of the Earth [m]
      *              [double *] Flattening form factor 
      *              [double *] Rectangular X coordinate [m]
      *              [double *] Rectangular Y coordinate [m]
      *              [double *] Rectangular Z coordinate [m]
      *  \param [out] [double *] Latitude coordinate on the ellipsoid [rad]
      *               [double *] Longitude coordinate on the ellipsoid [rad]
      *               [double *] Height coordinate on the ellipsoid [m]
      */
      void gconv2_( double *a, double *f, double *x, double *y, double *z, 
      	            double *phi, double *lambda, double *h );
      
      /**
      *  \b Description: \n
      *        subroutine evaluates the model of polar motion for
      *        a nonrigid Earth due to tidal gravitation. This polar motion
      *        is equivalent to the so-called "subdiurnal nutation".
      *        The model is a sum of a first order polynomial and 25 trigonometric 
      *        terms (15 long periodic and 10 quasi diurnal) with coefficients 
      *        given  in Table 5.1a of the IERS Conventions (2010).
      *  \param [in] [double *] modified julian date
      *  \param [out] [double *] PM variations as array w/ 2 elements (X and Y)
      */
      void pmsdnut2_( double *rmjd, double *pm );
      
      /**
      *  \b Description: \n
      *        This subroutine evaluates the model of subdiurnal libration
      *        in the axial component of rotation, expressed by UT1 and LOD.
      *        This effect is due to the influence of tidal gravitation on the
      *        departures of the Earth's mass distribution from the rotational
      *        symmetry, expressed by the non-zonal components of geopotential.
      *        The amplitudes have been computed for an elastic Earth with liquid
      *        core. The adopted truncation level is 0.033 microseconds in UT1
      *        corresponding to the angular displacement of 0.5 microarcseconds
      *        or to 0.015 mm at the planet surface. With this truncation level
      *        the model contains 11 semidiurnal terms. The coefficients of
      *        the model are given in Table 5.1b of the IERS Conventions (2010).
      *  \param [in] [double *] modified julian date
      *  \param [out] [double *] dUT1
      *  \param [out] [double *] dLOD
      */
      void utlibr_( double *rmjd, double *dut, double *dlod );

      /**
      *  \b Description: \n
      *        This subroutine computes the effects of the free core nutation.  
      *        Please note that the table is updated each year.
      *        The parameter N needs to be incremented for each additional
      *        year in the table. The updated table is maintained at the website
      *        http://syrte.obspm.fr/~lambert/fcn/
      *  \param [in] [double *] modified julian date (strictly TDB but TT makes no difference)
      *  \param [out] [double *] X  - CIP offset x component, in microas 
      *  \param [out] [double *] Y  - CIP offset y component, in microas 
      *  \param [out] [double *] dX - Uncertainty of x component, in microas
      *  \param [out] [double *] dY - Uncertainty of y component, in microas
      */
      void fcnnut_( double *mjd, double *x, double *y, double *dx, double *dy );

      /**
      *  \b Description: \n
      *        This subroutine computes the station tidal displacement
      *        caused by lunar and solar gravitational attraction. 
      *        The computations are calculated by the following steps:
      *
      *        Step 1): General degree 2 and degree 3 corrections + CALL ST1IDIU 
      *                 + CALL ST1ISEM + CALL ST1L1.
      *  
      *        Step 2): CALL STEP2DIU + CALL STEP2LON
      *
      *        It has been decided that the Step 3 non-correction for permanent tide
      *        would not be applied in order to avoid a jump in the reference frame.
      *        This Step 3 must be added in order to get the non-tidal station position
      *        and to conform with the IAG Resolution.
      *  \param [in] [double *]  xSta - 3x1 Geocentric position of the IGS station [m]
      *  \param [in] [double *]  xSun - 3x1 Geocentric position of the Sun [m]
      *  \param [in] [double *]  xMon - 3x1 Geocentric position of the Moon [m]
      *  \param [in] [int *]  yr - Year [UTC]
      *  \param [in] [int *]  month - Month [UTC]
      *  \param [in] [int *]  day - Day of Month [UTC]
      *  \param [in] [double *]  hr - Hour in the day [UTC]
      *  \param [out] [double *] dXtide - 3x1 Displacement vector 
      */
      void dehanttideinel_( double *xSta, int *yr, int *month, int *day, double *hr,
                            double *xSun, double *xMon, double *dXtide );
      
      /**
      *  \b Description: \n
      *        This subroutine computes the tidal displacements, using an expanded set
      *        of tidal constituents, whose amplitudes and phases are found by
      *        spline interpolation of the tidal admittance.  A total of 342
      *        constituent tides are included, which gives a precision of about 0.1%.
      *  \param [in] [int *]  1x5 date (yyyy, doy, hh, mm, ss) [UTC]
      *  \param [in] [double *] 11x3 amplitudes in (transposed) BLQ format
      *  \param [in] [double *] 11x3 phases in (transposed) BLQ format
      *  \param [out] [double *] 3x1 Displacement vector 
      */
      void hardisp_( int *it, float *amp, float *pha, double *disp );
      
      /**
      *  \b Description: \n
      *        compute the diurnal and semidiurnal variations in Earth Orientation 
      *        Parameters (x,y, UT1) from ocean tides.
      *  \param [in] [double *] modified julian date
      *  \param [out] [double *] 3x1 - delta_x, delta_y [microarcseconds] and
      *                          delta_UT1 [microseconds ]
      */
      void ortho_eop_( double *mjd, double *dEop );

      /**
      *  \b Description: \n
      *        compute the diurnal and semidiurnal variations in Earth Orientation 
      *        Parameters (x,y, UT1) from ocean tides, using various models.
      *  \param [in] [double *] modified julian date
      *  \param [in] [char *] file with eop tidal coefficients
      *  \param [out] [double *] 3x1 - delta_x, delta_y [microarcseconds] and
      *                          delta_UT1 [microseconds ]
      */
      void calc_hfeop_( double *time, char * lhfeop_file, double *delta_t, double *eop );

      /**
      *  \b Description: \n
      *        This subroutine evaluates the effects of zonal Earth tides on the
      *        rotation of the Earth.  The model used is a combination of Yoder
      *        et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
      *        body tide, and Kantha et al. (1998) ocean tide models 
      *        as recommended by the IERS Conventions (2010).  Refer to
      *        Chapter 8 pp. 105 - 106.
      *  \param [in] [double *]  TT, Julian centuries since J2000 (strictly TDB 
      *                          but TT makes no difference)
      *  \param [out] [double *] dUT - Effect on UT1 [s]
      *  \param [out] [double *] dLod - Effect on excess length of day (LOD) [s/day]
      *  \param [out] [double *] dOmega - Effect on rotational speed [rad/s]
      */

     
      void rg_zont2_( double *t, double *dUt, double *dLod, double *dOmega );
      
      /**
      *  \b Description: \n
      *        This subroutine determines the Global Mapping Functions GMF (Boehm et al. 2006)
      *  \param [in] [double *] mjd - Modified Julian Date
      *  \param [in] [double *] lat - Latitude given in radians (North Latitude)
      *  \param [in] [double *] lon - Longitude given in radians (East Longitude)
      *  \param [in] [double *] hgt - Height in meters (mean sea level)
      *  \param [in] [double *] zd  - Zenith distance in radians                  
      *  \param [out] [double *] gmfh - Hydrostatic mapping function
      *  \param [out] [double *] gmfw - Wet mapping function
      */
      void gmf_( double *mjd, double *lat, double *lon, double *hgt, double *zd,
                 double *gmfh, double *gmfw );
      
      /**
      *  \b Description: \n
      *        This subroutine determines Global Pressure and Temperature (Boehm et al. 2007)
      *        based on Spherical Harmonics up to degree and order 9.
      *  \param [in] [double *] mjd - Modified Julian Date
      *  \param [in] [double *] lat - Latitude given in radians (North Latitude)
      *  \param [in] [double *] lon - Longitude given in radians (East Longitude)
      *  \param [in] [double *] hgt - Height in meters (mean sea level)
      *  \param [out] [double *] p - Pressure given in hPa
      *  \param [out] [double *] t - Temperature in degrees Celsius
      *  \param [out] [double *] undu - Geoid undulation in meters
      */
      void gpt_( double *mjd, double *lat, double *lon, double *hgt, double *p,
                 double *t, double *undu );
      
      /**
      *  \b Description: \n
      *        This subroutin determines pressure, temperature, temperature lapse rate, water
      *        vapour pressure, hydrostatic and wet mapping function coefficients ah and aw, 
      *        and geoid undulation for specific sites near the Earth surface. It is 
      *        based on a 5 x 5 degree external grid file ('gpt2_5.grd') 
      *  \param [in] [double *] mjd - Modified Julian Date
      *  \param [in] [double *] lat - Latitude given in radians (North Latitude)
      *  \param [in] [double *] lon - Longitude given in radians (East Longitude)
      *  \param [in] [double *] hgt - Height in meters (mean sea level)
      *  \param [in] [int *] nstat - Number of stations in DLAT, DLON, and HELL
      *  \param [in] [int *] it - case 1: no time variation but static quantities 
      *                           case 0: with time variation (annual and semiannual)
      *  \param [in] [char *] grd - absolute path to gpt2_5.grd file
      *  \param [out] [double *] p - Pressure given in hPa
      *  \param [out] [double *] t - Temperature in degrees Celsius
      *  \param [out] [double *] dt - Temperature lapse rate in degrees per km
      *  \param [out] [double *] e - Water vapour pressure in hPa
      *  \param [out] [double *] ah - hydrostatic mapping function coefficient at zero height (VMF1)
      *  \param [out] [double *] aw - wet mapping function coefficient (VMF1)
      *  \param [out] [double *] undu - Geoid undulation in meters
      *  \param [in] [int ] size of grd file
      */
      void gpt2_( double *mjd, double *lat, double *lon, double *hgt, int *nstat, 
                  int *it,char * grd, double *p, double *t, double *dt, double *e, double * ah,
                  double *aw, double *undu, int sizeof_grd );
       /**
      *  \b Description: \n
      *        This subroutin determines pressure, temperature, temperature lapse rate, water
      *        vapour pressure, hydrostatic and wet mapping function coefficients ah and aw, 
      *        and geoid undulation for specific sites near the Earth surface. It is 
      *        based on a 5 x 5 degree external grid file ('gpt2_5.grd') 
      *  \param [in] [double *] mjd - Modified Julian Date
      *  \param [in] [double *] lat - Latitude given in radians (North Latitude)
      *  \param [in] [double *] lon - Longitude given in radians (East Longitude)
      *  \param [in] [double *] hgt - Height in meters (mean sea level)
      *  \param [in] [int *] nstat - Number of stations in DLAT, DLON, and HELL
      *  \param [in] [int *] it - case 1: no time variation but static quantities 
      *                           case 0: with time variation (annual and semiannual)
      *  \param [in] [char *] grd - absolute path to gpt2_5.grd file
      *  \param [out] [double *] p - Pressure given in hPa
      *  \param [out] [double *] t - Temperature in degrees Celsius
      *  \param [out] [double *] dt - Temperature lapse rate in degrees per km
      *  \param [out] [double *] Tm - Weighted mean atm temperature 
      *  \param [out] [double *] e - Water vapour pressure in hPa
      *  \param [out] [double *] ah - hydrostatic mapping function coefficient at zero height (VMF1)
      *  \param [out] [double *] aw - wet mapping function coefficient (VMF1)
      *  \param [out] [double *] la - water vapour decrease factor
      *  \param [out] [double *] undu - Geoid undulation in meters
      *  \param [out] [double *] Gn_h -  hydrostatic north gradient in m
      *  \param [out] [double *] Ge_h -  hydrostatic east gradient in m
      *  \param [out] [double *] Gn_w -  wet north gradient in m
      *  \param [out] [double *] Ge_w -  wet east gradient in m
      *  \param [in] [int ] size of grd file
      */
      void gpt3_1_( double *mjd, double *lat, double *lon, double *hgt, int *nstat, 
		    int *it,char * grd, double *p, double *t, double *dt, double *Tm,
		    double *e, double * ah,
		    double *aw, double *la, double *undu,
		    double *Gn_h,  double *Ge_h, double *Gn_w, double *Ge_w,int sizeof_grd );     
      /**
      *  \b Description: \n
      *        This subroutin determines the Vienna Mapping Function 1 
      *        (VMF1, site dependent version). The coefficients can be obtained 
      *        from the website http://ggosatm.hg.tuwien.ac.at/DELAY/SITE/
      *  \param [in] [double *] ah - Hydrostatic coefficient a
      *  \param [in] [double *] aw - Wet coefficient a
      *  \param [in] [double *] mjd -Modified Julian Date
      *  \param [in] [double *] lat -Latitude given in radians (North Latitude)
      *  \param [in] [double *] zd - Zenith distance in radians
      *  \param [out] [double *] vmf1h - Hydrostatic mapping function
      *  \param [out] [double *] vmf1w - Wet mapping function
      */
      void vmf1_( double *ah, double *aw, double *mjd, double *lat,
                    double *zd, double *vmf1h, double *vmf1w );
      /**
      *  \b Description: \n
      *        This subroutin determines the Vienna Mapping Function 3
      *        (VMF3, site dependent version). The coefficients can be obtained 
      *        from the website http://vmf.geo.tuwien.ac.at/
      *  \param [in] [double *] ah - Hydrostatic coefficient a
      *  \param [in] [double *] aw - Wet coefficient a
      *  \param [in] [double *] mjd -Modified Julian Date
      *  \param [in] [double *] lat -Latitude given in radians (North Latitude)
      *  \param [in] [double *] lon -Longitude given in radians (East Longitude)
      *  \param [in] [double *] zd - Zenith distance in radians
      *  \param [out] [double *] vmf1h - Hydrostatic mapping function
      *  \param [out] [double *] vmf1w - Wet mapping function
      */
	  void vmf3_( double *ah, double *aw, double *mjd, double *lat, double *lon,
                    double *zd, double *vmf1h, double *vmf1w );
      /**
      *  \b Description: \n
      *        Determines the Vienna Mapping Function 1(VMF1, grid version) including
      *        height correction (VMF1, grid version). The coefficients can be obtained 
      *        from the website http://ggosatm.hg.tuwien.ac.at/DELAY/GRID
      *  \param [in] [double *] ah - Hydrostatic coefficient a
      *  \param [in] [double *] aw - Wet coefficient a
      *  \param [in] [double *] mjd -Modified Julian Date
      *  \param [in] [double *] lat -Latitude given in radians (North Latitude)
      *  \param [in] [double *] hgt - Ellipsoidal height given in meters
      *  \param [in] [double *] zd - Zenith distance in radians
      *  \param [out] [double *] vmf1h - Hydrostatic mapping function
      *  \param [out] [double *] vmf1w - Wet mapping function
      */
      void vmf1_ht_( double *ah, double *aw, double *mjd, double *lat, 
                     double *hgt, double *zd, double *vmf1h, double *vmf1w );
      /**
      *  \b Description: \n
      *        Determines the Vienna Mapping Function 3(VMF3, grid version) including
      *        height correction (VMF3, grid version). The coefficients can be obtained 
      *        from the website http://ggosatm.hg.tuwien.ac.at/DELAY/GRID
      *  \param [in] [double *] ah - Hydrostatic coefficient a
      *  \param [in] [double *] aw - Wet coefficient a
      *  \param [in] [double *] mjd -Modified Julian Date
      *  \param [in] [double *] lat -Latitude given in radians (North Latitude)
      *  \param [in] [double *] lon -Longitude given in radians (East Longitude)
      *  \param [in] [double *] hgt - Ellipsoidal height given in meters
      *  \param [in] [double *] zd - Zenith distance in radians
      *  \param [out] [double *] vmf1h - Hydrostatic mapping function
      *  \param [out] [double *] vmf1w - Wet mapping function
      */
      void vmf3_ht_( double *ah, double *aw, double *mjd, double *lat, double *lon, 
                     double *hgt, double *zd, double *vmf1h, double *vmf1w );
      /**
      *  \b Description: \n
      *        Determines the asymmetric delay d in meters caused by gradients. 
      *        The north and east gradients are also provided they should be used
      *        with the gradient model by Chen and Herring (1997). 
      *  \param [in] [double *] lat - Latitude given in radians (North Latitude)
      *  \param [in] [double *] lon - Longitude given in radians (East Longitude)
      *  \param [in] [double *] as - Azimuth from north in radians
      *  \param [in] [double *] al - Elevation angle in radians
      *  \param [out] [double *] delay - Delay in meters
      *  \param [out] [double *] grN - North gradient in mm
      *  \param [out] [double *] grE - East gradient in mm
      */
      void apg_( double *lat, double *lon, double * az, double *el, double *delay,
                 double *grN, double *grE );
   
      /**
      *  \b Description: \n
      *        This function computes TCB-TCG at the geocenter using an approximation 
      *        of the time ephemeris TE405. It approximates TE405 (including the trend) 
      *        with errors of 0.453 ns (RMS) and 2.248 ns (maximum) for the period 
      *        1600-2200. Although the Julian date is, formally, barycentric dynamical 
      *        time (TDB), the terrestrial dynamical time (TT) can be used with no practical
      *        effect on the accuracy of the prediction. 
      *  \param [in] [double *] Julian date in days 
      *  \param [return] [double] TCB-TCG [s]
      */
      double hf2002_iers_( double *jd );
   
      /**
      *  \b Description: \n
      *  This routine is part of the International Earth Rotation and
      *  Reference Systems Service (IERS) Conventions software collection.
      *
      * This subroutine provides the angular coordinates of the IERS Conventional Mean Pole (CMP)
      * to be used in the analysis of space geodesy data after 1970.
      * Starting with the version CMP(2015), the coordinates are
      * based on the table of values from ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab
      * See IERS Conventions Section 7.1.4 at http://tai.bipm.org/iers/convupdt/convupdt.html for details.
      * The subroutine also provides previous versions of the CMP in the IERS Conventions (2003) and (2010)
      *
      *  \param [in] [int *] Year of the conventional model.  Limited to integer values of 2003, 2010, 2015
      *  \param [in] [double *] Date for which the angular coordinates of the Conventional Mean Pole are 
      *                         desired. Units are decimal years and fraction, e.g.1942.65 
      *  \param [out] [double *] x - Angular coordinate x of conventional mean pole [as]
      *  \param [out] [double *] y - Angular coordinate y of conventional mean pole [as]
      *  \param [out] [int *] error - flag indicating possible error in requested epoch or version.  
      *                               Requesting an invalid version or an epoch before 1970 returns error code -1 and mean pole coordinates (0,0)
      *                               Requesting an epoch out of the range 1975.0-2003.0 in version 2003 returns error code 1   
      *                               Requesting an epoch out of the range 1975.0-2010.0 in version 2010 returns error code 2  
      *                               Requesting an epoch out of the range 1970.0-2016.2 in version 2015 returns error code 3 
      *                               Requesting a valid version and an epoch in the recommended range returns error code 0
      */
      void iers_cmp_2015_( int *version, double *epoch, double *x, double *y, int *error );
        
   /**
    *  \b Description: ITRF2014 - CATREF psd-model\n 
    * Compute the post-seismic deformation/correction using parametric models
    * #    Model
    * 0    PWL (Piece-Wise Linear Function)
    * 1    Logarithmic Function
    * 2    Exponential Function
    * 3    Logarithmic + Exponential
    * 4    Two Exponential Functions
    *
    * \param [in] [int *] modn: model #
    * \param [in] [double *] dtq : time difference (t - t_Earthquake) in decimal year (but see note below)
    * \param [in] [double *] a1: amplitude 1 of the parametric model, if modn = 1 or 2 (or 3 or 4, if a2 & t2 are supplied)
    * \param [in] [double *] t1: relaxation time 1, if modn = 1 or 2 (or 3 or 4, if a2 & t2 are supplied)
    * \param [in] [double *] a2: amplitude 2 of the parametric model, if modn = 3 or 4
    * \param [in] [double *] t2: relaxtaion time 2, if modn = 3 or 4
    * \param [out][double *] defo: post-seismic correction
    */
   void parametric_(int *modn,double *dtq,double *a1,double *t1,double *a2,
                    double *t2,double *defo);
   }
  
   // constants from chater 1 IERS Conventions 2010
   const double GM    = 3.986005e14; 
   const double G     = 6.67428e-11;
   const double a_E   = 6378136.6;
   const double f_E   = 1.0/298.25642;
   const double rho_w = 1025.0;
   const double g_E   = 9.7803278;
   const double omega = 7.292115e-5;
  
   const double optl_gamma2_re = 0.6870;
   const double optl_gamma2_im = 0.0036;
   
   // PPN parameter gamma = 1 (p 163 of IERS 2010 Conv.)
   const double ppn_gamma = 1.0;
   
   const double d_tai_tt_sec = 32.184;
   
   // IERS Conventions 2010: eq (7.29)
   static const double optl_h_p() 
   { 
	   return ( sqrt( 8*M_PI/15.0 )*pow( omega,2 )*pow( a_E, 4)/GM ); 
   };
   
   static const double optl_k() 
   { 
	   return ( 4.0*M_PI*G*a_E*rho_w*optl_h_p()/( 3.0*g_E ) ); 
   };
   
    // mass parameters from DE421 (IERS Conv. 2010 Tab 3.1)): GM [m^3/s^2]]
    static std::map<std::string, double> jpl_gm
       { {"Mercury", 22032.090000*1e9},
         {"Venus", 324858.592000*1e9},
         {"Earth", 398600.436233*1e9},
         {"Mars", 42828.375214*1e9},
         {"Jupiter", 126712764.800000*1e9},
         {"Saturn", 37940585.200000*1e9},
         {"Uranus", 5794548.600000*1e9},
         {"Neptune", 6836535.000000*1e9},
         {"Moon", 4902.800076*1e9},
         {"Sun", 1.32712442076e20}
       };
   
}

#endif // IERSWRAP_H
