/***************************************************************************** 
 * USE OF IERS SUBROUTINE GPT: determines Global Pressure and Temperature    *
 *                             (Boehm et al. 2007) based on Spherical        *
 *                             Harmonics up to degree and order 9            *
 *                                                                           *
 * Test case:                                                                *
 *                                                                           *
 *       given input: DMJD = 55055D0                                         *
 *                    DLAT = 0.6708665767D0 radians (NRAO, Green Bank, WV)   *
 *                    DLON = -1.393397187D0 radians                          *
 *                    DHGT = 812.546 meters                                  *
 *       expected output: PRES = 918.0710638757363995D0 hPa                  *
 *                        TEMP = 19.31914181012882992D0 degrees Celsius      *
 *                        UNDU = -42.19185643717770517D0 meters              *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *                                                                           *
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
  double mjd = 55055.0;
  double lat = 0.6708665767;
  double lon = -1.393397187;
  double h   = 812.546;

  double pres0 = 918.0710638757363995;
  double temp0 = 19.31914181012882992; 
  double undu0 = -42.19185643717770517; 
  double pres, temp, undu;
  
  gpt_( &mjd, &lat, &lon, &h, &pres, &temp, &undu );

  cout << "difference to test case" << endl;
  cout << "pressure    = " << pres-pres0 << endl;
  cout << "temperature = " << temp-temp0 << endl;
  cout << "undulation  = " << undu-undu0 << endl;
}
