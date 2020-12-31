/***************************************************************************** 
 * USE OF IERS SUBROUTINE FCNNUT: computes the effects of the free core      *
 *                                nutation                                   *
 *                                                                           *
 * Test case:                                                                *
 *     given input: MJD = 54790D0   Modified Julian Date, TDB                *
 *                                                                           *
 *     expected output:  X = -176.8012290066270680D0 microarcseconds         *
 *                       Y = -93.51855308903756736D0 microarcseconds         *
 *                       dX = 3.745573770491803067D0 microarcseconds         *
 *                       dY = 3.745573770491803067D0 microarcseconds         *
 *                                                                           *
 * data table has to be updated each year following                          *
 *       http://syrte.obspm.fr/~lambert/fcn/                                 *
 * => test case is only fullfilled with original data                        *
 *                                                                           *
 * 2014-12-03 - TA                                                           *
 ****************************************************************************/
#include<iostream> 
#include"iers_wrapper.h"


using namespace std;
using namespace iers;
			  
int main( const int argc, const char **argv )
{
  double dy, dx;
  
  double mjd = 54790.0;
  double X = -176.801229006627068;
  double Y = -93.5185530890375673;
  
  double dx0 = 3.74557377049180306;
  double dy0 = 3.74557377049180306;
  
  fcnnut_( &mjd, &X, &Y, &dx, &dy );

  cout << "difference to test case" << endl;
  cout << "dX = " << dx-dx0 << endl;
  cout << "dY = " << dy-dy0 << endl;
}
