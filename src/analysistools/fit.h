#ifndef FIT_H
#define FIT_H

#include <vector>
#include <iterator>
#include <sstream>
#include <math.h> 
#include <cstdlib>
#include "tictoc.h"
#include "matrix.h"
#include "lsa.h"
#include <libconfig.h++>

/**
* @brief Fit class
* @author SH - bakkari developer team
* @date 2015-04-22
* @version 0.1
*/

//using namespace ivg;
using namespace libconfig;
using namespace std;

namespace ivgat
{

// ===========================================================================
class Fit
// ===========================================================================
{
   public:

   // ==============================================
   // =============== Constructors: ================
   // ==============================================
   /**
   *  \b Description: \n
   *        Default constructor
   *  \param [in] no input parameters needed
   *  \return An instance of the class 'Fit'
   */
   Fit( );



   // ==============================================
   // =========== MEMBER-functions: ================
   // ==============================================

   /**
   *  \b Description: \n
   *        Method gets residuals
   * \param [in] no input parameters needed
   * \return [ivg::Matrix] residuals
   */
   ivg::Matrix get_resid(){ return _r; };

   /**
   *  \b Description: \n
   *        Method gets parameters (coefficients)
   * \param [in] no input parameters needed
   * \return [ivg::Matrix] parameter (coefficients)
   */
   ivg::Matrix get_param(){ return _x; };   
   
   /**
   *  \b Description: \n
   *        This method calculates the coefficients of a polynomial of a given degree that 
   *        fits the data in a least squares sense.
   *  \param [in] [ivg::Matrix] Matrix of mesh points
   *              [ivg::Matrix] Matrix of observations
   *              [int] polynomial degree
   * \return [ivg::Matrix] curve fit functional values
   */
   ivg::Matrix polyfit( ivg::Matrix t, ivg::Matrix y, int degree );

   
   /**
   *  \b Description: \n
   *        This method calculates the coefficients of an exponential function 
   *        that fits the data in a least squares sense.
   *  \param [in] [ivg::Matrix] Matrix of mesh points
   *              [ivg::Matrix] Matrix of observations
   * \return [ivg::Matrix] curve fit  functional values
   */   
   ivg::Matrix expfit( ivg::Matrix t, ivg::Matrix y );

   ivg::Matrix pwpoly( ivg::Matrix x, ivg::Matrix b, int n_int, int deg, int deriv);
   ivg::Matrix pwpoly( ivg::Matrix x, ivg::Matrix b, double dx_int, int deg, int deriv );
   ivg::Matrix pwpoly( ivg::Matrix x, ivg::Matrix b, ivg::Matrix xint, int deg, int deriv );

   ivg::Matrix cspline( ivg::Matrix x, ivg::Matrix y, int m );

   private:

   // ==============================================
   // ======== class variables / attributes: =======
   // ==============================================

   // residuals
   ivg::Matrix _r;
   ivg::Matrix _x;
};

} // # namespace ivgat

#endif // FIT_H
