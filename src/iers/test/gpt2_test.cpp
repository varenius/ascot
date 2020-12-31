/***************************************************************************** 
 * USE OF IERS SUBROUTINE GPT: determines pressure, temperature, temperature *
 *                             lapse rate, water vapour pressure, hydrostatic*
 *                             and wet mapping function coefficients ah and  *
 *                             aw, and geoid undulation for specific sites   *
 *                             near the Earth surface.                       *
 *                             It is based on a 5 x 5 degree external grid   *
 *                             file ('gpt2_5.grd')                           *
 *                                                                           *
 * Test case:                                                                *
 *                                                                           *
 *       Example 1 (Vienna, 2 August 2012, with time variation):             *
 *                                                                           *
 *      dmjd = 56141.d0                                                      *
 *      dlat(1) = 48.20d0*pi/180.d0                                          *
 *      dlon(1) = 16.37d0*pi/180.d0                                          *
 *      hell(1) = 156.d0                                                     *
 *      nstat = 1                                                            *
 *      it = 0                                                               *
 *                                                                           *
 *      output:                                                              *
 *      p = 1002.56 hPa                                                      *
 *      T = 22.12 deg Celsius                                                *
 *      dT = -6.53 deg / km                                                  *
 *      e = 15.63 hPa                                                        *
 *      ah = 0.0012647                                                       *
 *      aw = 0.0005726                                                       *
 *      undu = 44.06 m                                                       *
 *                                                                           *
 *      Example 2 (Vienna, 2 August 2012, without time variation,            *
 *      i.e. constant values):                                               *
 *                                                                           *
 *      dmjd = 56141.d0                                                      *
 *      dlat(1) = 48.20d0*pi/180.d0                                          *
 *      dlon(1) = 16.37d0*pi/180.d0                                          *
 *      hell(1) = 156.d0                                                     *
 *      nstat = 1                                                            *
 *      it = 1                                                               *
 *                                                                           *
 *      output:                                                              *
 *      p = 1003.49 hPa                                                      *
 *      T = 11.95 deg Celsius                                                *
 *      dT = -5.47 deg / km                                                  *
 *      e = 9.58 hPa                                                         *
 *      ah = 0.0012395                                                       *
 *      aw = 0.0005560                                                       *
 *      undu = 44.06 m                                                       *
 *                                                                           *
 *                                                                           *
 * 2014-12-08 - TA                                                           *
 ****************************************************************************/
#include<iostream> 
#include<vector> 
#include<math.h> 
#include"iers_wrapper.h"


using namespace std;
using namespace iers;
			  
int main( const int argc, const char **argv )
{
   double mjd  = 56141.0;
   double lat  = 48.20*M_PI/180.0;
   double lon  = 16.37*M_PI/180.0;
   double hell = 156.0;            
   int nstat = 1;                   
   int it = 0;                      

   double p0 = 1002.56;   
   double T0 = 22.12;    
   double dT0 = -6.53;   
   double e0 =  15.63;    
   double ah0 = 0.0012647;
   double aw0 = 0.0005726;
   double undu0 =44.06; 
  
   double pres, temp, dtemp, e, ah, aw, undu;
 
   gpt2_( &mjd, &lat, &lon, &hell, &nstat, &it, &pres, &temp, &dtemp, &e, 
          &ah, &aw, &undu );

   cout << "differences to test case 1" << endl;
   cout << "pressure     = " << pres-p0 << endl;
   cout << "temperature  = " << temp-T0 << endl;
   cout << "dtemperature = " << dtemp-dT0 << endl;
   cout << "undulation   = " << undu-undu0 << endl;
   cout << "mf hydr      = " << ah-ah0 << endl;
   cout << "mf wet       = " << aw-aw0 << endl;

   it = 1;
   p0 = 1003.49;
   T0 = 11.95; 
   dT0 = -5.47;
   e0 = 9.58;
   ah0 = 0.0012395;
   aw0 = 0.0005560;
   undu0 = 44.06;
  
   gpt2_( &mjd, &lat, &lon, &hell, &nstat, &it, &pres, &temp, &dtemp, &e, 
          &ah, &aw, &undu );

   cout << "differences to test case 2" << endl;
   cout << "pressure     = " << pres-p0 << endl;
   cout << "temperature  = " << temp-T0 << endl;
   cout << "dtemperature = " << dtemp-dT0 << endl;
   cout << "undulation   = " << undu-undu0 << endl;
   cout << "mf hydr      = " << ah-ah0 << endl;
   cout << "mf wet       = " << aw-aw0 << endl;

   return 0;
}
