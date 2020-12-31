/***************************************************************************** 
 * USE OF IERS SUBROUTINE PMSDNUT2: calculates so-called "subdiurnal         *
 *                                  nutation" in polar motion                *
 *                                                                           *
 * Test case:                                                                *
 *     given input:  rmjd_a = 44239.1 ( January 1, 1980 2:24.00 )            *
 *                   rmjd_b = 55227.4 ( January 31, 2010 9:35.59 )           *
 *                                                                           *
 *     expected output: dUT1_a =   2.441143834386761746D0 mus;               *
 *                      dLOD_a = -14.78971247349449492D0 mus / day           *
 *                      dUT1_b = - 2.655705844335680244D0 mus;               *
 *                      dLOD_b =  27.39445826599846967D0 mus / day           *
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
  vector<double> ut( 2,0.0 );

  // test case 1
  double mjd = 44239.1;
  
  double ut0  = 2.441143834386761746;
  double lod0 = -14.7897124734944949;
  
  utlibr_( &mjd, &ut[0], &ut[1] );

  cout << "difference to test case 1" << endl;
  cout << "dUT1 = " << ut.at(0)-ut0 << endl;
  cout << "dLOD = " << ut.at(1)-lod0 << endl;


  // test case 2
  mjd = 55227.4;
  
  ut0  = -2.65570584433568024;
  lod0 = 27.3944582659984696;
  
  utlibr_( &mjd, &ut[0], &ut[1] );

  cout << endl << "difference to test case 2" << endl;
  cout << "dUT1 = " << ut.at(0)-ut0 << endl;
  cout << "dLOD = " << ut.at(1)-lod0 << endl;
}
