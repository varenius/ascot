/***************************************************************************** 
 * USE OF IERS SUBROUTINE GMF: determines the Global Mapping Functions GMF   *
 *                             (Boehm et al. 2006)                           *
 *                                                                           *
 * Test case:                                                                *
 *                                                                           *
 *    given input: DMJD = 55055D0                                            *
 *                 DLAT = 0.6708665767D0 radians (NRAO, Green Bank, WV)      *
 *                 DLON = -1.393397187D0 radians                             *
 *                 DHGT = 844.715D0 meters                                   *
 *                 ZD   = 1.278564131D0 radians                              *
 *                                                                           *
 *    expected output: GMFH = 3.425245519339138678D0                         *
 *                     GMFW = 3.449589116182419257D0                         *
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
  double h   = 844.715;
  double zd  = 1.278564131;

  double mfh0 = 3.425245519339138678;
  double mfw0 = 3.449589116182419257; 
  double mfh, mfw;
  
  gmf_( &mjd, &lat, &lon, &h, &zd, &mfh, &mfw );

  cout << "difference to test case" << endl;
  cout << "hydrostatic = " << mfh-mfh0 << endl;
  cout << "wet         = " << mfw-mfw0 << endl;
}
