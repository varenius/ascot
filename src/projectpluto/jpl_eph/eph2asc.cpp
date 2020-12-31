/* eph2asc.cpp: convert binary JPL ephemerides back to ASCII

Copyright (C) 2011, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

   This program is possibly useless.  It is the opposite of the quite
useful 'asc2eph':  given a _binary_ JPL ephemeris file,  it can recreate
the ASCII ephemeris files from which the binary version was made.  Usage
is described in error_exit( ) below.

   The ASCII files will sometimes differ from the original input in the
last decimal place.  */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "jpleph.h"

#define jpl_get_cache( ephem) (*(double **)((char *)ephem + 224 + 48))

/* The format of doubles in JPL ASCII ephemeris files is a little strange. */
/* The usual 'E' exponent is replaced with 'D'.  A leading zero is used,   */
/* so that 3.14159E+000 would become 0.314159D+01.  Two-digit exponents    */
/* are used instead of three-digit ones.  Put it all together,  and you    */
/* have the need for the following:                                        */

static void put_in_buff( char *buff, const double val)
{
   sprintf( buff, "%25.15lE", val);
   sprintf( buff + 20, "D%+03d", atoi( buff + 21) + 1);
   buff[1] = buff[2];   /* grab minus sign if any */
   buff[4] = buff[3];
   buff[2] = '0';       /* insert leading zero */
   buff[3] = '.';
}

static void error_exit( const int exit_code)
{
   printf( "'eph2asc' requires as command line arguments the name of a\n");
   printf( "JPL DE binary ephemeris and a starting and ending year.  For\n");
   printf( "example:\n\neph2asc unix.405 1950 1975\n\n");
   printf( "would result in an ASCII ephemeris covering the years 1950 to 1975\n");
   printf( "being written to the standard output.\n");
   exit( exit_code);
}

/* The conversion works as follows:

   The (binary) ephemeris is opened,  and the command-line years are converted
to JD.  We figure out what the "real" starting JD is,  given the ephemeris
step size.  Then we step through,  one block at a time (i.e.,  jd += ephem_step)
and compute the positions of Jupiter and Saturn for a time at the center of
the current block.  The only reason we do so is that it results in the block
of coefficients being loaded up in a cache,  accessible via jpl_get_cache( ).
So we can just write them out in ASCII form.

   Eventually,  'jd' reaches the ending year specified on the command line,
and we're done.  */

int main( const int argc, const char **argv)
{
   void *eph = jpl_init_ephemeris( argv[1], NULL, NULL);
   int rval = -1;

   if( argc < 4)
      error_exit( -1);
   eph = jpl_init_ephemeris( argv[1], NULL, NULL);
   if( !eph)
      {
      printf( "Ephemeris '%s' not opened; error %d\n",
                       argv[1], jpl_init_error_code( ));
      error_exit( -2);
      }
   else
      {
      const double ephem_step = jpl_get_double( eph, JPL_EPHEM_STEP);
      const double start_jd = jpl_get_double( eph, JPL_EPHEM_START_JD);
      const double J2000 = 2451545.0;
      const int n_coeffs = jpl_get_long( eph, JPL_EPHEM_KERNEL_NCOEFF);
      const int n_to_write = ((n_coeffs + 2) / 3) * 3;
      double jd = (atof( argv[2]) - 2000.) * 365.25 + J2000, vect[6];
      const double jd_end = (atof( argv[3]) - 2000.) * 365.25 + J2000;
      double *coeffs = jpl_get_cache( eph);
      int zval, i, rec_num = 1;

      jd = floor( (jd - start_jd) / ephem_step) * ephem_step + start_jd;
      while( jd <= jd_end)
         {
         zval = jpl_pleph( eph, jd + ephem_step / 2., 5, 6, vect, 0);
         if( !zval)
            {
            printf( "%6d%6d\n", rec_num, n_coeffs);
            rec_num++;
            for( i = 0; i < n_to_write; i++)
               {
               char buff[40];

               put_in_buff( buff, (i < n_coeffs) ? coeffs[i] : 0.);
               printf( "  %s%s", buff, (i % 3 == 2) ? "\n" : "");
               }
            }
         jd += ephem_step;
         }
      jpl_close_ephemeris( eph);
      rval = 0;
      }
   return( rval);
}
