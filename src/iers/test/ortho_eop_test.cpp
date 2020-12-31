/***************************************************************************** 
 * USE OF IERS SUBROUTINE ORTHO_EOP: computes the diurnal and semidiurnal    *
 *                                   variations in Earth Orientation         *
 *                                   Parameters (x,y, UT1) from ocean tides  *
 *                                                                           *
 * Test case:                                                                *
 *                                                                           *
 *     given input: MJD = 47100D0                                            *
 *     expected output: delta_x = -162.8386373279636530D0 microarcseconds    *
 *                      delta_y = 117.7907525842668974D0 microarcseconds     *
 *                      delta_UT1 = -23.39092370609808214D0 microseconds     *
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
  double mjd = 47100.0;
  vector<double> eop( 3,0.0 );
  
  double x0 = -162.8386373279636530;
  double y0 =  117.7907525842668974; 
  double u0 =  -23.39092370609808214;
  
  ortho_eop_( &mjd, &eop[0] );

  cout << "difference to test case" << endl;
  cout << "dx-pole = " << eop.at(0)-x0 << endl;
  cout << "dy-pole = " << eop.at(1)-y0 << endl;
  cout << "dUT1    = " << eop.at(2)-u0 << endl;
}
