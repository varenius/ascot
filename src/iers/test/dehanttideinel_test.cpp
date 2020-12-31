/***************************************************************************** 
 * USE OF IERS SUBROUTINE DEHANTTIDEINEL: computes the station tidal         *
 *                                        displacement caused by lunar and   *
 *                                        solar gravitational attraction     *
 *                                                                           *
 *  Test case:                                                               *
 *     given input: XSTA(1) = 4075578.385D0 meters                           *
 *                  XSTA(2) =  931852.890D0 meters                           *
 *                  XSTA(3) = 4801570.154D0 meters                           *
 *                  XSUN(1) = 137859926952.015D0 meters                      *
 *                  XSUN(2) = 54228127881.4350D0 meters                      *
 *                  XSUN(3) = 23509422341.6960D0 meters                      *
 *                  XMON(1) = -179996231.920342D0 meters                     *
 *                  XMON(2) = -312468450.131567D0 meters                     *
 *                  XMON(3) = -169288918.592160D0 meters                     *
 *                  YR      = 2009                                           *
 *                  MONTH   = 4                                              *
 *                  DAY     = 13                                             *
 *                  FHR     = 0.00D0 hour                                    *
 *                                                                           *
 *     expected output:  DXTIDE(1) = 0.7700420357108125891D-01 meters        *
 *                       DXTIDE(2) = 0.6304056321824967613D-01 meters        *
 *                       DXTIDE(3) = 0.5516568152597246810D-01 meters        *
 *                                                                           *
 *  Test case:                                                               *
 *     given input: XSTA(1) =  1112189.660D0 meters                          * 
 *                  XSTA(2) = -4842955.026D0 meters                          *
 *                  XSTA(3) =  3985352.284D0 meters                          *
 *                  XSUN(1) = -54537460436.2357D0 meters                     *
 *                  XSUN(2) =  130244288385.279D0 meters                     *
 *                  XSUN(3) =  56463429031.5996D0 meters                     *
 *                  XMON(1) =  300396716.912D0 meters                        *
 *                  XMON(2) =  243238281.451D0 meters                        *
 *                  XMON(3) =  120548075.939D0 meters                        *
 *                  YR      = 2012                                           *
 *                  MONTH   = 7                                              *
 *                  DAY     = 13                                             *
 *                  FHR     = 0.00D0 hour                                    *
 *                                                                           *
 *     expected output:  DXTIDE(1) = -0.2036831479592075833D-01 meters       *
 *                       DXTIDE(2) =  0.5658254776225972449D-01 meters       *
 *                       DXTIDE(3) = -0.7597679676871742227D-01 meters       *
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
   vector<double> disp( 3,0.0 );
   vector<double> xSta( 3,0.0 );
   vector<double> xSun( 3,0.0 );
   vector<double> xMon( 3,0.0 );

   // test case 1
   xSta[0] = 4075578.385;
   xSta[1] =  931852.890;
   xSta[2] = 4801570.154; 
   xSun[0] = 137859926952.015;
   xSun[1] = 54228127881.4350;
   xSun[2] = 23509422341.6960;
   xMon[0] = -179996231.920342;
   xMon[1] = -312468450.131567;
   xMon[2] = -169288918.592160;
   int year    = 2009;
   int month   = 4;
   int day     = 13;
   double hr   = 0.00;
              
   double dX0  = 0.7700420357108125891e-1;
   double dY0  = 0.6304056321824967613e-1;
   double dZ0  = 0.5516568152597246810e-1;

   dehanttideinel_( &xSta[0], &year, &month, &day, &hr, &xSun[0], &xMon[0],
                    &disp[0] );
  
   cout << "difference to test case 1" << endl;
   cout << "dX = " << disp[0]-dX0 << endl;
   cout << "dY = " << disp[1]-dY0 << endl;
   cout << "dZ = " << disp[2]-dZ0 << endl;


   // test case 2
   xSta[0]  =  1112189.660;
   xSta[1]  = -4842955.026;
   xSta[2]  =  3985352.284;   
   xSun[0]  = -54537460436.2357;
   xSun[1]  =  130244288385.279;
   xSun[2]  =  56463429031.5996;
   xMon[0]  =  300396716.912;
   xMon[1]  =  243238281.451;
   xMon[2]  =  120548075.939;
   year  = 2012;
   month = 7;
   day   = 13;
   hr    = 0.00;
                   
   dX0   = -0.2036831479592075833e-1;
   dY0   =  0.5658254776225972449e-1;
   dZ0   = -0.7597679676871742227e-1;

   dehanttideinel_( &xSta[0], &year, &month, &day, &hr, &xSun[0], &xMon[0],
                    &disp[0] );
  
   cout << "difference to test case 2" << endl;
   cout << "dX = " << disp[0]-dX0 << endl;
   cout << "dY = " << disp[1]-dY0 << endl;
   cout << "dZ = " << disp[2]-dZ0 << endl;
}
