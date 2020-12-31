#include <math.h>
#include "norad.h"
#include "norad_in.h"

/*------------------------------------------------------------------*/

/* FMOD2P */
double FMod2p( const double x)
{
   double rval = fmod( x, twopi);

   if( rval < 0.)
      rval += twopi;
   return( rval);
} /* fmod2p */

#define EPHEM_TYPE_DEFAULT    '0'
#define EPHEM_TYPE_SGP        '1'
#define EPHEM_TYPE_SGP4       '2'
#define EPHEM_TYPE_SDP4       '3'
#define EPHEM_TYPE_SGP8       '4'
#define EPHEM_TYPE_SDP8       '5'
#define EPHEM_TYPE_HIGH       'h'

/*------------------------------------------------------------------*/

/* Selects the type of ephemeris to be used (SGP*-SDP*) */
int DLL_FUNC select_ephemeris( const tle_t *tle)
{
   double ao, xnodp, delo, a1, del1, r1, temp;
   int rval;

   if( tle->ephemeris_type == EPHEM_TYPE_HIGH)
      rval = 1;         /* force high-orbit state vector model */
   else if( tle->xno <= 0. || tle->eo > 1. || tle->eo < 0.)
      rval = -1;                 /* error in input data */
   else if( tle->ephemeris_type == EPHEM_TYPE_SGP4
         || tle->ephemeris_type == EPHEM_TYPE_SGP8)
      rval = 0;         /* specifically marked non-deep */
   else if( tle->ephemeris_type == EPHEM_TYPE_SDP4
         || tle->ephemeris_type == EPHEM_TYPE_SDP8)
      rval = 1;         /* specifically marked deep */
   else
      {
      /* Period > 225 minutes is deep space */
      a1 = pow( xke / tle->xno, two_thirds);
      r1 = cos(tle->xincl);
      temp = ck2 * 1.5 * (r1*r1*3.0-1.0) * pow( 1.0-tle->eo*tle->eo, -1.5);
      del1 = temp/(a1*a1);
      ao = a1 * (1.0 - del1 * (1./3. + del1 * (del1 * 1.654320987654321+1.0)));
      delo = temp/(ao*ao);
      xnodp = tle->xno / (delo + 1.0);

      /* Select a deep-space/near-earth ephemeris */
      /* If the object makes less than 6.4 revolutions around the earth... */
      if (twopi / (xnodp * minutes_per_day) >= (1. / 6.4))
         rval = 1;      /* yes,  it should be a deep-space (SDPx) ephemeris */
      else
         rval = 0;      /* no,  you can go with an SGPx ephemeris */
      }
   return( rval);
} /* End of select_ephemeris() */

/*------------------------------------------------------------------*/

long DLL_FUNC sxpx_library_version( void)
{
   return( 0x100);
}
