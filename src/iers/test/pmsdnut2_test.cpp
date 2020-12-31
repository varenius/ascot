/***************************************************************************** 
 * USE OF IERS SUBROUTINE PMSDNUT2: calculates so-called "subdiurnal         *
 *                                  nutation" in polar motion                *
 *                                                                           *
 * Test case:                                                                *
 *     given input: rmjd = 54335D0 ( August 23, 2007 )                       *
 *                                                                           *
 *     expected output: (dx) pm(1)  = 24.83144238273364834D0 microarcseconds *
 *                      (dy) pm(2) = -14.09240692041837661D0 microarcseconds *
 *                                                                           *
 *                                                                           *
 * 2014-12-03 - TA                                                           *
 ****************************************************************************/
#include<iostream> 
#include<vector> 
#include"iers_wrapper.h"


using namespace std;
using namespace iers;
			  
int main( const int argc, const char **argv )
{
  double mjd = 54335.0;
  vector<double> pm( 2,0.0 );
  
  double x0 =  24.8314423827336483;
  double y0 = -14.0924069204183766;
  
  pmsdnut2_( &mjd, &pm[0] );

  cout << "difference to test case" << endl;
  cout << "dx-pole = " << pm.at(0)-x0 << endl;
  cout << "dy-pole = " << pm.at(1)-y0 << endl;
}
