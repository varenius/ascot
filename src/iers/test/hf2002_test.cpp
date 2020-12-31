/***************************************************************************** 
 * USE OF IERS SUBROUTINE HF2002: computes TCB-TCG at the geocenter using an *
 *                                approximation of the time ephemeris TE405  *
 *                                                                           *
 * Test case:                                                                *
 *                                                                           *
 *  output should be identical to IERS program XHF2002_IERS                  *
 *                                                                           *
 *      expected output:                                                     *
 *   Year      JD(TT)      TCB-TCG / s                                       *
 *   1600 2305445.0   -176.1773555610856D0                                   *
 *   1700 2341970.0   -129.4460764716997D0                                   *
 *   1800 2378495.0    -82.7147442912425D0                                   *
 *   1900 2415020.0    -35.9834220078961D0                                   *
 *   2000 2451545.0     10.7478546790489D0                                   *
 *   2100 2488070.0     57.4792210328145D0                                   *
 *   2200 2524595.0    104.2104805049757D0                                   *
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
  vector<double> jd = { 2305445.0, 2341970.0, 2378495.0, 2415020.0, 2451545.0,
                        2488070.0, 2524595.0 };

  vector<double> expect = { -176.1773555610856, -129.4460764716997,
                             -82.7147442912425,  -35.9834220078961,
                              10.7478546790489,   57.4792210328145,
                             104.2104805049757 };

  cout << "difference to test case" << endl;
  for( int i=0; i<jd.size(); ++i )
  {
     double tcb_tcg = iers::hf2002_iers_( &jd[i] );

     cout << "dTCB-TCG = " << tcb_tcg-expect.at(i) << endl;
  }
}
