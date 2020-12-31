/* testeph.cpp: verify a JPL ephemeris

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

   2013 March:  (BJG) You can now specify the name of the test data
file on the command line,  and modify the test tolerance by a desired
scale factor.  (I found that for DE430,  the tolerance could be tightened
fifty-fold without difficulty;  i.e.,  adding the '-s.02' command line
option did not result in warnings.)

   2013 April:  (BJG) Modified to handle more ephemeris constants.   */

/*****************************************************************************
*        *****    jpl planetary and lunar ephemerides    *****     C ver.1.2 *
******************************************************************************
* program testeph                                                            *
*                                                                            *
*                                                                            *
* Testeph tests the jpl ephemeris reading and interpolating routine using    *
* examples computed from the original ephemeris.                             *
*                                                                            *
* Testeph contains the reading and interpolating subroutines that are of     *
* eventual interest to the user.  Once testeph is working correctly, the     *
* user can extract those subroutines and the installation process is         *
* complete.                                                                  *
*                                                                            *
* You must allow acces to "testpo.xxx" to testeph.                           *
* "testpo.xxx" is the specially formatted text file that contains the test   *
* cases for the ephmeris, dexxx.                                             *
*                                                                            *
* After the initial identifying text which is concluded by an "EOT" in       *
* columns 1-3, the test file contains the following quantities:              *
*                                                                            *
*     JPL ephemeris number                                                   *
*     calendar date                                                          *
*     julian ephemeris date                                                  *
*     target number (1-mercury, ...,3-earth, ,,,9-pluto, 10-moon, 11-sun,    *
*                    12-solar system barycenter, 13-earth-moon barycenter    *
*                    14-nutations, 15-librations)                            *
*     center number (same codes as target number)                            *
*     coordinate number (1-x, 2-y, ... 6-zdot)                               *
*     coordinate  [au, au/day].                                              *
*                                                                            *
* For each test case input, testeph                                          *
*                                                                            *
*     - computes the corresponding state from data contained                 *
*       in dexxx,                                                            *
*                                                                            *
*     - compares the two sets,                                               *
*                                                                            *
*     - writes an error message if the difference between                    *
*       any of the state components is greater than 10**(-13).               *
*                                                                            *
*     - writes state and difference information for every npt'th             *
*       test case processed.                                                 *
*                                                                            *
*                                                                            *
*  This program was written in standard fortran-77 and it was manually       *
*  translated to C language by Piotr A. Dybczynski (dybol@phys.amu.edu.pl).  *
*                                                                            *
*  This is version 1.2 of this C translation, use jplbin.h version 1.2       *
*
******************************************************************************
*                 Last modified: July 23, 1997 by PAD                        *
******************************************************************************
16 Mar 2001:  Revised by Bill J. Gray.  You can now use binary
ephemerides with either byte order ('big-endian' or 'small-endian');
the code checks to see if the data is in the "wrong" order for the
current platform,  and swaps bytes on-the-fly if needed.  (Yes,  this
can result in a slowdown... sometimes as much as 1%.  The function is
so mathematically intensive that the byte-swapping is the least of our
troubles.)  You can also use DE-200, 403, 404, 405,  or 406 and most,
if not all,  later ephemerides without recompiling (the constan( )
function now determines which ephemeris is in use and its byte order);
and you can set the TESTFILE and EPHFILE on the command line.

Also,  I did some minor optimization of the interp( ) (Chebyshev
interpolation) function,  resulting in a bit of a speedup.

   2013 Apr 6:  DE-430 has 572 constants.  We can only store the first
400 in the binary ephemeris;  some logic was needed to ensure that
constants after the first 400 were ignored.
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

/**** include variable and type definitions, specyfic for this C version */

#include "jpleph.h"

void error_exit( const int err_code)
{
   printf( "usage: testeph <ephemeris file> <options>\n\n");
   printf( "'testeph' requires the name of the JPL ephemeris file as a command\n");
   printf( "line argument.  It then looks for a 'testpo' (test positions) file\n");
   printf( "in the same folder with the same extension,  and checks the positions\n");
   printf( "against those computed from the ephemeris.\n\nOptions:\n");
   printf( "   -a       Report all errors,  not just the first incidence\n");
   printf( "   -fN      Report each Nth result.  Default is N=100.\n");
   printf( "   -tFILE   Get test input data from FILE instead of the default 'testpo'\n");
   exit( err_code);
}

/* The following maximum accepted errors were determined by actually looking at
what the maximum errors were,  for each ephemeris,  for values and rates,
for each axis.  In evaluating the "errors",  I accounted for machine precision
and roundoff in the test files.  (For example,  the 'testpo' files for DE-200,
DE-405,  and the earlier 'testpo.406' gave data rounded to 13 decimal places,
so you get errors up to 5e-14 from that source alone.  Other 'testpo' files
that I've seen go to 20 places.  See the comments about the 'roundoff_error'
variable below.)

   Accounting for machine precision error is trickier.  Somewhat arbitrarily,
I decided that after doing all the Chebyshev math,  it would be reasonable to
assume the values are good,  at best,  to one part in 1e+14.  So in testing,
the error abs( computed - value_from_file) must be less than the
'max_accepted_error' plus the roundoff error plus 1e-14 times the computed
value.  If it isn't,  an error message is shown.

   After doing this for all ephemerides I have from DE-102 to DE-432t,  I
found that the maximum error in a coordinate was 6.5e-15 AU,  or a little
under a millimeter.  The maximum error in a velocity was 1.2e-17 AU/day,  or
about two microns per day.  The nutations were all within zero error,  as were
the TT-TDB values and the lunar rotation (libration) angles.  However,  some
of the lunar rotation _rates_ had errors up to 8.9e-20 radians/day,  or about
2e-14 arcseconds/day.  (Note that by 'error',  I mean 'difference from the
testpo input file and values computed by this code'.  The actual difference
between these values and what the celestial objects are doing is certainly far
greater!)

   As you can see,  I set the 'max_accepted_error' to be a bit more than the
maximum error encountered in the various ephemerides.   */

static double max_accepted_error( const int ntarg, const int ncoord,
                              const int de_num)
{
   double rval = 0.;

   if( ntarg <= 13 && ncoord <= 3)
      rval = 1.e-14;    /* planet posn; 1e-14 AU = 1.5 mm */
   else if( ntarg <= 13 && ncoord >= 3)
      rval = 2e-17;     /* planet velocity; 2e-17 AU/day = 3 microns/day */
   else if( ntarg == 15 && ncoord > 3)       /* lunar libration rate: */
      rval = 1e-19;                    /* 1e-19 rad/day = 2e-14 arcsec/day */
   else
      rval = 0.;              /* everything else should be within limits */
   return( rval);
}

   /* At least at present,  there's no provision for storing more than 1018 */
   /* constants in a binary JPL ephemeris.  See 'asc2eph.cpp' for details.  */

#define JPL_MAX_N_CONSTANTS 1018

/***** THERE IS NO NEED TO MODIFY THE REST OF THIS SOURCE (I think) ******/

int main( const int argc, const char **argv)
{
  char nams[JPL_MAX_N_CONSTANTS][6], buff[102];
  double vals[JPL_MAX_N_CONSTANTS];
  double max_err_found[10][6];
  int i, j, line,  n_failures[10];
  int n_constants, n_columns;
  int output_frequency = 100;
  const char *ephfile_name = argv[1];
  double start_jd, end_jd;
  double roundoff_error = 0.;
  bool report_all_errors = false;
  bool pause_on_errors = true;
  FILE *testfile;
  clock_t timer;
  void *ephem;
  const char *test_file_name = NULL;
#if defined( __GNUC__) && !defined( __MINGW32__)
   const char path_separator = '/';
#else
   const char path_separator = '\\';
#endif
   const char *error_messages[8] = {
       "Result outside acceptable error tolerance",
       "Outside date range",
       "Read error",
       "No nutations in ephemeris",
       "No librations in ephemeris",
       "Invalid index",
       "Seek error",
       "No TT-TDB data in ephemeris" };
   const int n_errors = sizeof( error_messages) / sizeof( error_messages[0]);

/***** Write a fingerprint to the screen. ***********************************/

  setvbuf( stdout, NULL, _IONBF, 0);
  puts("\n JPL test-ephemeris program (v.1.2)\n"
       " C version translated from the original JPL FORTRAN code.\n");

   if( argc < 2)
      error_exit( -1);
   for( i = 0; i < 10; i++)
      for( j = 0; j < 6; j++)
         max_err_found[i][j] = 0.;

   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'a': case 'A':
               report_all_errors = true;
               break;
            case 'f': case 'F':
               output_frequency = atoi( argv[i] + 2);
               break;
            case 'n': case 'N':
               pause_on_errors = false;
               break;
            case 't': case 'T':
               test_file_name = argv[i] + 2;
               break;
            default:
               printf( "Unrecognized option '%s'\n", argv[i]);
               error_exit( -2);
               break;
            }

/****** Print the ephemeris constants. **************************************/

  ephem = jpl_init_ephemeris( ephfile_name, nams, vals);
  if( !ephem)
     {
     printf( "Ephemeris file '%s' not loaded\n", ephfile_name);
     printf( "Error code: %d\n", jpl_init_error_code( ));
     error_exit( -2);
     }
   else
      printf( "Ephemeris initialized\n");

   n_constants = (int)jpl_get_long( ephem, JPL_EPHEM_N_CONSTANTS);
   printf( "%d constants\n", n_constants);
   if( n_constants > JPL_MAX_N_CONSTANTS)
      n_constants = JPL_MAX_N_CONSTANTS;
   start_jd = jpl_get_double( ephem, JPL_EPHEM_START_JD),
   end_jd =   jpl_get_double( ephem, JPL_EPHEM_END_JD),
   printf("%.9f  %.9f  %.9f\n", start_jd, end_jd,
                           jpl_get_double( ephem, JPL_EPHEM_STEP));
   n_columns = (n_constants + 1) / 2;
   for( i = 0; i < n_columns; i++)
      {
      printf("%.6s  %24.16E",nams[i],vals[i]);
      if( i + n_columns < n_constants)
         printf("   %.6s  %24.16E",nams[i + n_columns],vals[i + n_columns]);
      printf( "\n");
      }

   printf( "emrat = %.15lf      AU = %.5lf\n",
                           jpl_get_double( ephem, JPL_EPHEM_EARTH_MOON_RATIO),
                           jpl_get_double( ephem, JPL_EPHEM_AU_IN_KM));

/****** Skip the test points file header comments.  *************************/
  if( test_file_name)
     strcpy( buff, test_file_name);
  else
     {
     const char *extension;

     strcpy( buff, ephfile_name);
     for( i = strlen( buff); i && buff[i - 1] != path_separator; i--)
        ;
     strcpy( buff + i, "testpo");
     extension = strchr( ephfile_name + i, '.');
     if( extension)
        strcat( buff + i, extension);
     else
        sprintf( buff + strlen( buff), ".%3ld",
            jpl_get_long( ephem, JPL_EPHEM_EPHEMERIS_VERSION));
     }

  testfile = fopen( buff, "r");
  if( !testfile)
     {
     printf( "Test data file '%s' not found\n", buff);
     error_exit( -3);
     }

   while( fgets( buff, 100, testfile) && memcmp( buff, "EOT", 3))
      ;

   puts(" LINE  JED    t# c# x#  --- JPL value ---   "
          "--- user value --   -- difference --");

   line=0;
   timer = clock( );
   for( i = 0; i < 10; i++)
      n_failures[i] = 0;
   while( fgets( buff, 100, testfile) != NULL)
   {
     int err_code, ntarg, nctr, ncoord;
     double del, et;
     double r[6], xi, xi_computed;
     bool fatal_error = false, report_this_error = false;

/*****  Read a value from the test case; Skip if not within the time-range
        of the present version of the ephemeris.                            */
     line++;
     if( sscanf( buff + 15," %lf %d %d %d %lf", &et, &ntarg, &nctr, &ncoord,
                       &xi) != 5)
        {
        printf( "Failure to parse line %d:\n%s\n", line, buff);
        fatal_error = true;
        }
     else
        {
        err_code = jpl_pleph(ephem, et, ntarg, nctr, r, 1);
        if( err_code > 0 || err_code < -n_errors)
           {
           printf( "Internal error:  unknown error code %d\n", err_code);
           fatal_error = true;
           }
        else if( !roundoff_error)
           {
           char *tptr = strchr( buff + 30, '.');

           assert( tptr);
           tptr++;
           roundoff_error = .5;
           while( *tptr >= '0' && *tptr <= '9')
              {
              tptr++;
              roundoff_error /= 10.;
              }
           }
        }
     if( fatal_error)
        {
        printf( "\nThis error shouldn't happen,  ever.  It indicates a bug\n");
        printf( "that should be fixed.\n");
        printf( "Please contact pluto@projectpluto.com and report this.\n");
        return( -1);
        }
     xi_computed = r[ncoord - 1];
     del = fabs( xi_computed - xi);
     if( err_code)
        {
        n_failures[-err_code]++;
        if( n_failures[-err_code] == 1 || report_all_errors)
           report_this_error = true;
        }
     else
        {
        const double tolerance = max_accepted_error( ntarg, ncoord,
                   jpl_get_long( ephem, JPL_EPHEM_EPHEMERIS_VERSION))
                               + fabs( xi) * 1e-14 + roundoff_error;
        const unsigned idx = (ntarg <= 13 ? 0 : ntarg - 13);

        if( max_err_found[idx][ncoord - 1] < del - tolerance)
           max_err_found[idx][ncoord - 1] = del - tolerance;
        if( del > tolerance)
           {
           n_failures[0]++;
           if( n_failures[-err_code] == 1 || report_all_errors)
              {
              report_this_error = true;
              printf( "*****  warning : next difference >= tolerance *****\n");
              printf( "%s", buff);
              }
           }
        }

     if( (!err_code && !(line % output_frequency)) || report_this_error)
        printf("%4d %10.1f %2d %2d %2d %25.20f %25.20f %22.20e\n",
               line,et,ntarg,nctr,ncoord,xi, xi_computed, xi_computed - xi);
     if( report_this_error)
        {
        printf( "Error message: '%s'\n", error_messages[-err_code]);
        if( err_code == JPL_EPH_OUTSIDE_RANGE)
           {
           const double J2000 = 2451545.;

           printf( "WARNING:  The test file tests items outside the range\n");
           printf( "of this ephemeris!\n");
           printf( "The input ephemeris file covers years from %.1lf to %.1lf.\n",
                 2000. + (start_jd - J2000) / 365.25,
                 2000. + (end_jd - J2000) / 365.25);
           printf( "The test is for the year %.1lf\n",
                 2000. + (et - J2000) / 365.25);
           printf( "It's common for DE files to be built that cover a subset of the\n");
           printf( "range of the original ephemerides,  so this may not be an error.\n");
           }
        if( !report_all_errors)
           {
           printf( "  Further errors of this type won't be shown,  but you'll get a count\n");
           printf( "  of how many are found. (Use the '-a' switch to report all errors.)\n");
           }
        if( pause_on_errors)
           {
           printf( "  Hit any key:\n");
           getchar( );
           }
        }
   }
   for( i = 0; i < 10; i++)
      for( j = 0; j < 6; j++)
         if( max_err_found[i][j])
            printf( "%d %d %.8le %c\n", i, j + 1, max_err_found[i][j],
                       (max_err_found[i][j] > 0.) ? '*' : ' ');

   printf( "%d lines read and tested in %.3lf seconds\n", line,
           (double)( clock( ) - timer) / (double)CLOCKS_PER_SEC);
   for( i = 0; i < n_errors; i++)
      if( n_failures[i])
         printf( "%d lines failed with error code %d ('%s').\n",
                          n_failures[i], -i, error_messages[i]);
   fclose( testfile);
   jpl_close_ephemeris( ephem);
   return( 0);
}
