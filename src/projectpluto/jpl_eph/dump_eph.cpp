/* dump_eph.cpp: dumps header data from a JPL ephemeris file

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
02110-1301, USA.    */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jpleph.h"

#define AU_IN_KM 1.49597870691e+8

/* A little test program to dump the header data from a JPL ephemeris file.
   Also,  if run in the form

dump_eph ephem_filename start_jd step n_steps

   then ephemerides will be shown for eleven objects (nine planets,  the
moon,  and the solar system barycenter) for 'n_steps' dates,  starting on
'start_jd',  for every 'step' days.   */

int main( const int argc, const char **argv)
{
   void *p;

   if( argc < 2)
      {
      printf( "dump_eph requires the name of a JPL ephemeris file as a\n");
      printf( "command-line argument.  Add another argument,  and the\n");
      printf( "ephemeris constants will be dumped.\n");
      return( -1);
      }

   p = jpl_init_ephemeris( argv[1], NULL, NULL);

   if( !p)
      {
      printf( "JPL data not loaded from '%s'\n", argv[1]);
      printf( "Error code: %d\n", jpl_init_error_code( ));
      return( -1);
      }
   else
      {
      const double start_jd = jpl_get_double( p, JPL_EPHEM_START_JD);
      const double end_jd = jpl_get_double( p, JPL_EPHEM_END_JD);
      const double j2000 = 2451545.;
      const int n_constants = (int)jpl_get_long( p, JPL_EPHEM_N_CONSTANTS);
      int i, j;

      printf( "Ephemeris runs from JD %.3lf to %.3lf (years %.3lf to %.3lf)\n",
                  start_jd, end_jd,
                  2000. + (start_jd - j2000) / 365.25,
                  2000. + (end_jd - j2000) / 365.25);
      printf( "Stepsize is %lf days\n", jpl_get_double( p, JPL_EPHEM_STEP));
      printf( "1 AU = %lf km\n", jpl_get_double( p, JPL_EPHEM_AU_IN_KM));
      printf( "Earth/Moon = %lf\n", jpl_get_double( p, JPL_EPHEM_EARTH_MOON_RATIO));
      printf( "Ephemeris version DE%ld\n", jpl_get_long( p, JPL_EPHEM_EPHEMERIS_VERSION));
      printf( "Kernel size: %ld\n", jpl_get_long( p, JPL_EPHEM_KERNEL_SIZE));
      printf( "Record size: %ld\n", jpl_get_long( p, JPL_EPHEM_KERNEL_RECORD_SIZE));
      printf( "N coeffs: %ld\n", jpl_get_long( p, JPL_EPHEM_KERNEL_NCOEFF));
      printf( "Byte swap: %ld\n", jpl_get_long( p, JPL_EPHEM_KERNEL_SWAP_BYTES));
      printf( "   0    1    2    3    4    5    6    7    8    9   10   11   12   13\n");
      printf( "  Mer  Ven  EMB  Mar  Jup  Sat  Ura  Nep  Plu  Moo  Sun  Nut  Lib  TT-TDB\n");
      for( j = 0; j < 3; j++)
         {
         for( i = 0; i < 15; i++)
            printf( "%5ld", jpl_get_long( p,
                        JPL_EPHEM_IPT_ARRAY + (i * 3 + j) * sizeof( int32_t)));
         printf( "\n");
         }
      printf( "%d constants\n", n_constants);
      if( argc > 2)        /* dump constants,  too */
         {

         for( i = 0; i < n_constants; i++)
            {           /* show constants in two columns: */
            const int idx = i / 2 + (i & 1) * (n_constants + 1) / 2;
            char constant_name[7];
            const double ephem_constant = jpl_get_constant( idx, p, constant_name);

            printf( "%.6s  %24.16E%s", constant_name, ephem_constant,
                        (i & 1) ? "\n" : "   ");
            }
         printf( "\n");
         }
      if( argc > 4)
         {
         double jd = atof( argv[2]);
         const double step = atof( argv[3]);
         int n_steps = atoi( argv[4]);

         while( n_steps--)
            {
            int i;
            double state_vect[6];

            printf( "JD %lf\n", jd);
            for( i = 0; i < 11; i++)
               {
               const char *format_str;
               static const char *object_names[] = {
                     "SSBar", "Mercu", "Venus", "EMB  ", "Mars ",
                     "Jupit", "Satur", "Uranu", "Neptu", "Pluto",
                     "Moon " };

               if( i == 10)       /* the moon */
                  {
                  jpl_pleph( p, jd, 3, 10, state_vect, 1);
                  state_vect[0] *= AU_IN_KM;  /* cvt AU to kilometers */
                  state_vect[1] *= AU_IN_KM;
                  state_vect[2] *= AU_IN_KM;
                  format_str = "%22.12lf %22.12lf %22.12lf\n";
                  }
               else
                  {
                  jpl_pleph( p, jd, (i ? i : 11), 12, state_vect, 1);
                  format_str = "%22.18lf %22.18lf %22.18lf\n";
                  }
               printf( "%s ", object_names[i]);
               printf( format_str,
                          state_vect[0], state_vect[1], state_vect[2]);
               }
            for( i = 0; i < 3; i++)
               {
               static const char *text[3] = {
                  "Sun posn", "Sun vel ", "Sun acc " };
               double *tptr = jpl_get_pvsun( p) + i * 3;

               printf( "%s ", text[i]);
               printf( "%22.18lf %22.18lf %22.18lf\n",
                        tptr[0], tptr[1], tptr[2]);
               }
            jd += step;
            }
         }
      jpl_close_ephemeris( p);
      }
   return( 0);
}
