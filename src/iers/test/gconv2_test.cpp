/***************************************************************************** 
 * USE OF IERS SUBROUTINE GCONV2: transforms from geocentric rectangular to  *
 *                                geodetic coordinates.                      *
 *                                                                           *
 * Test case:                                                                *
 *     given input: x = 4075579.496D0 meters  Wettzell (TIGO) station        *
 *                  y =  931853.192D0 meters                                 *
 *                  z = 4801569.002D0 meters                                 *
 *                                                                           *
 *     expected output: phi    =   0.857728298603D0 radians                  *
 *                      lambda =   0.224779294628D0 radians                  *
 *                      h      = 665.9207D0 meters                           *
 *                                                                           *
 *                                                                           *
 * 2014-12-03 - TA                                                           *
 ****************************************************************************/
#include<iostream> 
#include"iers_wrapper.h"


using namespace std;
using namespace iers;
			  
int main( const int argc, const char **argv )
{
  double phi, lambda, hell;
  double x = 4075579.496;
  double y =  931853.192;
  double z = 4801569.002;
  double a, e;
  
  gconv2_( &a,&e,&x,&y,&z, &phi, &lambda, &hell );
  
  double phi0    =   0.85772829860;
  double lambda0 =   0.22477929462;
  double h0      = 665.9207;
  
  cout << "difference to test case" << endl;
  cout << "dphi    = " << phi-phi0 << endl;
  cout << "dlambda = " << lambda-lambda0 << endl;
  cout << "dh      = " << hell-h0 << endl;
}
