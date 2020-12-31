#include "fit.h"

namespace ivgat
{


// ===========================================================================
//              constructors and destructor
// ===========================================================================

// ...........................................................................
Fit::Fit( )
// ...........................................................................
{
}


// ...........................................................................
ivg::Matrix Fit::polyfit( ivg::Matrix t, ivg::Matrix b, int degree )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 1
   cerr << "+++ Fit::polyfit( ivg::Matrix t, ivg::Matrix b, int degree )" << endl; 
   tictoc tim;
   tim.tic();
#endif

   ivg::Matrix A( t.rows(),degree+1,1.0 );
   ivg::Matrix t_coeff( degree+1,1,t(0) );
   for( int deg=1; deg<=degree; ++deg )
      A.set_sub( 0,deg,t^deg );

   ivg::Matrix w( b.length(), 1, 1.0 );
   ivg::Lsa sol( A, b, w, t );
   
   sol.solve_neq( 0 );
   _x = sol.get_parameters();

   // save residuals
   _r = sol.get_resid();

#if DEBUG_ANALYSIS >= 1
   cerr << "--- Fit::polyfit( ivg::Matrix t, ivg::Matrix b, int degree )" 
        << " : " << tim.toc() << " s " << endl; 
#endif

   return ( b - _r );
}


// ...........................................................................
ivg::Matrix Fit::expfit( ivg::Matrix t, ivg::Matrix b )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 1
   cerr << "+++ Fit::expfit( ivg::Matrix t, ivg::Matrix b )" << endl; 
   tictoc tim;
   tim.tic();
#endif

   ivg::Matrix A( t.rows(),2,1.0 );
   A.set_sub( 0,1,t );
  
   ivg::Matrix w( b.length(), 1, 1.0 );
     
   // re-cast the general exponential form y = C exp( Ax ) to a linear problem
   // of type Y = DX + B by using ln(y) = ln( C exp( Ax ) ) = ln(C) + Ax;
   // solve linear problem and transform parameter: C = exp( coeff1 ); D = coeff2;
   ivg::Matrix b2 = b.log();
   ivg::Lsa sol( A, b2, w, t );
   
   sol.solve_neq( 0 );
   _x = sol.get_parameters();
   
   double C = exp( _x(0) );
   double D = _x(1);
   
   ivg::Matrix f = (t*D).exp() * C;
  

#if DEBUG_ANALYSIS >= 1
   cerr << "--- Fit::expfit( ivg::Matrix t, ivg::Matrix b )" 
        << " : " << tim.toc() << " s " << endl; 
#endif

   return ( f );
}


// ...........................................................................
ivg::Matrix Fit::cspline( ivg::Matrix x, ivg::Matrix y, int m )
// ...........................................................................
{
   // reduce mean value to increase numerical stability
   double ave = y.meanD();
   y = y-ave;

   // transform x-values in range of B-Splines
   double len = (double)(m-1)/(double)y.length();
   x = x-x.min();
   x = x*1.0/x.max()*(x.length()-1)*len+len/2.0;

   // Jacobianmatrix
   ivg::Matrix A( x.length(),m,0.0 );

   for( int i=0;i<x.length();++i )
   {
      for( int j=0;j<m;++j )
      {
         // calculate contribution of observation to B-Spline (displaced by j)
         if( x(i) <= j-2 )
            A(i,j) = 0.0;
         else if( x(i) <= j-1 )
            A(i,j) =  1.0/6.0*pow( x(i)-j+2.0,3.0 );    
         else if( x(i) <= j )
            A(i,j) = 1.0/6.0*pow(x(i)-j+2.0, 3.0 )-4.0/6.0*pow( x(i)-j+1,3.0 );
         else if( x(i) <= j + 1 )
            A(i,j) = 1.0/6.0*pow(-x(i)+j+2.0, 3.0 )-4.0/6.0*pow( -x(i)+j+1,3.0 );
         else if ( x(i) <= j + 2 )
            A(i,j) =  1.0/6.0*pow( -x(i)+j+2.0,3.0 );    
         else
            A(i,j) = 0.0;
      }
   }

   ivg::Matrix w( x.length(),1,1.0 );
   ivg::Lsa sol( A,y,w );
   sol.solve_neq( 0 );

   // save residuals
   _r = sol.get_resid();
   _x = sol.get_parameters();
   
   return( y-_r+ave );
}

                                                                                  

// ...........................................................................
ivg::Matrix Fit::pwpoly( ivg::Matrix x, ivg::Matrix b, ivg::Matrix xint, 
                         int deg, int deriv )
// ...........................................................................
{
    // reduce mean value to increase numerical stability
    double ave = b.meanD();
    b = b-ave;

    // vector of parameter epochs, i.e. mid of intervals
    ivg::Matrix xparam( xint.rows()-1,1,0.0 );

    // find oberservations in intervals
    std::map< int,vector<int> > int_idx;
    for( int i=0; i<xint.length()-1; ++i )
       int_idx[i] = x.find_idx( ge, xint(i), le, xint(i+1) );

    // loop over intervals and construct jacobian matrix for LS adjustment
    ivg::Matrix A( x.length(),(deg+1)*(xint.length()-1),0.0 );
    for( int i=0; i<xint.length()-1; ++i )
    {
       xparam( i ) = ( xint(i+1)+xint(i) )/2.0;
       ivg::Matrix ti = x.get_sub( int_idx[i], {0} )-xparam(i);
       for( int d=0; d<=deg; ++d )
          A.set_sub( int_idx[i],{i*(deg+1)+d},ti^d );
    }

    ivg::Matrix wgt( b.rows(),1,1.0 );
    ivg::Lsa solution( A,b,wgt,x );

    // build pseudo observations to force equal derivatives at interval borders
    ivg::Matrix F;
    for( int j=0; j<deriv; ++j )
    {
       ivg::Matrix E( xint.length()-2, A.cols(), 0.0 );
       for( int i=0; i<xint.length()-2; ++i )
       {
          double dt1 = xparam(i)-xint(i);
          double dt2 = xparam(i+1)-xint(i+1);
          for( int d=j; d<=deg; ++d )
          {
             double v1 = pow(  dt1,(d-j) );
             double v2 = pow( -dt2,(d-j) );
             for( int k=0; k<j; ++k )
             {
                v1 *= (double)(d-k);
                v2 *= (double)(d-k);
             }
             E( i,i*(deg+1)+d ) = v1;
             E( i,(i+1)*(deg+1)+d ) = -v2;
          }
       }
       if( F.rows() == 0 )
          F = E;
       else
          F.append_rows( E );
    }


    // set constraints for solution 
    if( F.rows() != 0 )
    {
       ivg::Matrix w( F.rows(),1,1e6 );
       solution.update_constraints( F,w );
    }

    // and solve
    solution.solve_neq( 0 );

    _r = solution.get_resid();

    return( b-_r+ave );
}

// ...........................................................................
ivg::Matrix Fit::pwpoly( ivg::Matrix x, ivg::Matrix b, double dx_int, int deg, 
                         int deriv )
// ...........................................................................
{                                                                                  
    // interval length => epochs of interval borders
    ivg::Matrix int_borders( x.min(), dx_int, x.max()+dx_int-dx_int/100.0, 1 );

    ivg::Matrix y = pwpoly( x,b,int_borders,deg,deriv );
    return y;
}                                                                                  
                                                                                   
// ...........................................................................
ivg::Matrix Fit::pwpoly( ivg::Matrix x, ivg::Matrix b, int n_int, int deg, 
                         int deriv)
// ...........................................................................
{                                                                                  
    // #intervals => interval length
    double dt_int = ( x.max()-x.min() )/(double)n_int;

    ivg::Matrix y = pwpoly( x,b,dt_int,deg,deriv );
    return y;
}                                                                                  

} // # namespace ivgat
