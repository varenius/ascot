/***************************************************************************** 
 * USE OF IERS SUBROUTINE RG_ZONT2: evaluates the effects of zonal Earth     *
 *                                  tides on the  rotation of the Earth.     *
 *                                                                           *
 * Test case:                                                                *
 *                                                                           *
 *    given input: T = .07995893223819302 Julian centuries since J2000       *
 *                 (MJD = 54465)                                             *
 *    expected output: DUT    =  7.983287678576557467E-002 seconds           *
 *                     DLOD   =  5.035331113978199288E-005 seconds / day     *
 *                     DOMEGA = -4.249711616463017E-014 radians / second     *
 *                                                                           *
 *                                                                           *
 * 2014-12-08 - TA                                                           *
 ****************************************************************************/
#include<iostream> 
#include<vector> 
#include"iers_wrapper.h"


using namespace std;
			  
int main( const int argc, const char **argv )
{
  double t = .07995893223819302;
  vector<double> dut( 3,0.0 );
  
  double dut0  =  7.983287678576557467e-002;
  double dlod0 =  5.035331113978199288e-005; 
  double dom0  = -4.249711616463017e-014;
  
  iers::rg_zont2_( &t, &dut[0], &dut[1], &dut[2] );

  cout << "difference to test case" << endl;
  cout << "dUT     = " << dut.at(0)-dut0 << endl;
  cout << "dLOD    = " << dut.at(1)-dlod0 << endl;
  cout << "dOmega  = " << dut.at(2)-dom0 << endl;
}
