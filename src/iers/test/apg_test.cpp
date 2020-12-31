/***************************************************************************** 
 * USE OF IERS SUBROUTINE GPT: determines Global Pressure and Temperature    *
 *                             (Boehm et al. 2007) based on Spherical        *
 *                             Harmonics up to degree and order 9            *
 *                                                                           *
 * Test case:                                                                *
 *                                                                           *
 *     Kashima 11 Station information retrieved at:                          *
 *     ftp://ivscc.gsfc.nasa.gov/pub/config/ns/kashim11.config.txt           *
 *                                                                           *
 *     given input: DLAT = 0.6274877539940092D0 radians (KASHIMA 11, Japan)  *
 *                  DLON = 2.454994088489240D0 radians                       *
 *                  AZ   = 0.2617993877991494D0 radians                      *
 *                  EL   = 0.8726646259971648D0 radians                      *
 *                                                                           *
 *     expected output: D   = -0.9677190006296187757D-4 meters               *
 *                      GRN = -0.1042668498001996791D0 mm                    *
 *                      GRE = 0.4662515377110782594D-1 mm                    *
 *                                                                           *
 *                                                                           *
 * 2014-12-08 - TA                                                           *
 ****************************************************************************/
#include<iostream> 
#include<vector> 
#include"iers_wrapper.h"


using namespace std;
using namespace iers;
			  
int main( const int argc, const char **argv )
{
  double lat = 0.6274877539940092;
  double lon = 2.454994088489240;
  double az  = 0.2617993877991494;
  double el  = 0.8726646259971648;

  double delay0 = -0.9677190006296187757e-4;
  double nGr0 = -0.1042668498001996791;
  double eGr0 = 0.4662515377110782594e-1; 
  double dtau, nGr, eGr;
  
  apg_( &lat, &lon, &az, &el, &dtau, &nGr, &eGr );

  cout << "difference to test case" << endl;
  cout << "delay          = " << dtau-delay0 << endl;
  cout << "north gradient = " << nGr-nGr0 << endl;
  cout << "east gradient  = " << eGr-eGr0 << endl;
}
