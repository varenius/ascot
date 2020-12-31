/***************************************************************************** 
 * USE OF IERS SUBROUTINE IERS_CMP_2015: computes coordinate x of 
 *                                       conventional mean pole              *
 *                                                                           *
 * Test cases:                                                               *
 *                                                                           *
 *   1)  INPUT:    2015   2015.35                                            *
 *       OUTPUT:     0.12523135017430187       0.34623285979490870    0      * 
 *                                                                           *
 *   2)  INPUT:    2003   2002.1D0                                           *
 *       OUTPUT:     5.57429999999999246E-002  0.36529499999999965    0      *
 *                                                                           *
 *   3)  INPUT:    2010   2002.1D0                                           *
 *       OUTPUT:     6.06820922651758257E-002  0.34962260201001832    0      *
 *                                                                           *
 *   4)  INPUT:    2015   2002.1D0                                           *
 *       OUTPUT:     6.07152994763746900E-002  0.34770609736572355    0      *
 *                                                                           *
 *   5)  INPUT:    2003   2003.1D0                                           *
 *       OUTPUT:     5.65729999999999220E-002  0.36924499999999960    1      *
 *                                                                           *
 *                                                                           *
 * 2015-07-07 - TA                                                           *
 ****************************************************************************/
#include<iostream> 
#include<vector> 
#include"iers_wrapper.h"


using namespace std;
using namespace iers;
			  
int main( const int argc, const char **argv )
{
   vector<int> versions( 5 );
   versions.at(0) = 2015;
   versions.at(1) = 2003;
   versions.at(2) = 2010;
   versions.at(3) = 2015;
   versions.at(4) = 2003;

   vector<double> epochs( versions.size() );
   epochs.at(0) = 2015.35;
   epochs.at(1) = 2002.10;
   epochs.at(2) = 2002.10;
   epochs.at(3) = 2002.10;
   epochs.at(4) = 2003.10;
  
   vector<double> x0( versions.size() );
   x0.at(0) = 0.12523135017430187;   
   x0.at(1) = 5.57429999999999246e-2;
   x0.at(2) = 6.06820922651758257e-2;
   x0.at(3) = 6.07152994763746900e-2;
   x0.at(4) = 5.65729999999999220e-2;
   vector<double> y0( versions.size() );
   y0.at(0) = 0.34623285979490870;
   y0.at(1) = 0.36529499999999965;
   y0.at(2) = 0.34962260201001832;
   y0.at(3) = 0.34770609736572355;
   y0.at(4) = 0.36924499999999960;
   vector<int> er( versions.size() );
   er.at(0) = 0;
   er.at(1) = 0;
   er.at(2) = 0;
   er.at(3) = 0;
   er.at(4) = 1;
  

   double x,y;
   int err;
   for( int i=0; i<versions.size();++i )
   {
      iers_cmp_2015_( &versions.at(i), &epochs.at(i), &x, &y, &err );
      cout << endl << "difference to test case " << i << endl;
      cout << "x-pole   = " << x-x0.at(i) << endl;
      cout << "y-pole   = " << y-y0.at(i) << endl;
      cout << "err code = " << err << "/" << er.at(i) << endl << endl;
   }
}
