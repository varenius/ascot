#include "../iers_wrapper.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
using namespace std; 

//extern"C" 
//  {
//     void calc_hfeop_( double *time, char * lhfeop_file, double *delta_t, double *eop );
//}

int main( int argc, char *argv[] )
{
  double mjd = 58000;
  double erp[3];
  
  string file="iers2010_xyu.txt";
  char* cstr=new char[file.length()+1];
  std::strcpy (cstr, file.c_str());
  //char* cstr=(char *)file.c_str();
  char* tst=(char *)"iers2010_xyu.txt";
  std::cout << cstr << endl;
  double delta_t=32.184+37;
   iers::ortho_eop_( &mjd, erp );
  std::cout << "ORTHO_EOP: Xpol: " << erp[0] << " Ypol: " << erp[1] << " dUT1: " << erp[2] << std::endl;
  iers::calc_hfeop_(&mjd,cstr,&delta_t,erp);
  std::cout << file << ": Xpol: " << erp[0] << " Ypol: " << erp[1] << " dUT1: " << erp[2] << std::endl;
  //  file=(char*)"desai_model_jgrb51665-sup-0002-ds01.txt";
  iers::calc_hfeop_(&mjd,(char *) file.c_str(),&delta_t,erp);
  std::cout << file << ": Xpol: " << erp[0] << " Ypol: " << erp[1] << " dUT1: " << erp[2] << std::endl;
}
