/* sub_eph.cpp: extract a subsection of a JPL ephemeris

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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "jpleph.h"
#include "jpl_int.h"
#include "watdefs.h"
#include "date.h"

/* Code to extract a subsection of a JPL DE file.  Suppose you've got a
large DE file covering,  say,  the years -3000 to +3000 (span of the
full DE-406 file),  and would like to create a smaller one that just
covers the years 1900 to 2100.  This code would accomplish that task. */

/* The JPL DE file headers include the coverage span in 'calendar'
(year-month-day) form.  This 'adjust_epoch' function puts the date
into a buffer in the usual DE format. */

static void adjust_epoch( char *buff, const long jd)
{
   char tbuff[40];

   sprintf( tbuff, "%9ld.5 ", jd);
   full_ctime( tbuff + 12, (long)( jd + 1.001),
             FULL_CTIME_YMD | FULL_CTIME_DATE_ONLY);
   tbuff[18] -= 32;           /* capitalize month names */
   tbuff[19] -= 32;
   memcpy( buff, tbuff, strlen( tbuff));
}

int DLL_FUNC make_sub_ephem( void *ephem, const char *sub_filename,
                              double start_jd, double end_jd)
{
   struct jpl_eph_data *eph = (struct jpl_eph_data *)ephem;
   FILE *ofile = fopen( sub_filename, "wb");
   int rval = -1;

   if( start_jd < eph->ephem_start)
      start_jd = eph->ephem_start;
   if( end_jd > eph->ephem_end)
      end_jd = eph->ephem_end;

   if( !eph->ifile)
      rval = -99;
   else if( ofile)
      {
      char *buff = (char *)eph->cache;
      long i;
      long start_block = (long)floor(
                          (start_jd - eph->ephem_start) / eph->ephem_step);
      long end_block = (long)ceil(
                          (end_jd - eph->ephem_start) / eph->ephem_step);

      start_jd = eph->ephem_start + (double)start_block * eph->ephem_step;
      end_jd   = eph->ephem_start + (double)end_block   * eph->ephem_step;
      printf( "Start JD: %.4lf (block %ld)\n", start_jd, start_block);
      printf( "End   JD: %.4lf (block %ld)\n", end_jd, end_block);
      printf( "recsize %ld\n", (long)eph->recsize);

      rval = 0;

      if( fseek( eph->ifile, 0L, SEEK_SET))
         rval = -2;
                  /* read the very first file record... */
      else if( fread( buff, 1, eph->recsize, eph->ifile) != (size_t)eph->recsize)
         rval = -3;
                  /* ...adjust the starting and ending dates... */
      *(double *)( buff + 2652) = start_jd;
      *(double *)( buff + 2660) = end_jd;
                  /* ...swap bytes in those two doubles if we have to... */
      if( eph->swap_bytes)
         {
         char tbuff[16];

         for( i = 0; i < 16; i++)
            tbuff[i] = buff[2652 + (i ^ 7)];
         for( i = 0; i < 16; i++)
            buff[2652 + i] = tbuff[i];
         }
      adjust_epoch( buff + 101, (long)start_jd);
      adjust_epoch( buff + 185, (long)end_jd);
                  /* and write out the newly-revised first record: */
      if( !rval)
         {
         if( fwrite( buff, 1, eph->recsize, ofile) != (size_t)eph->recsize)
            rval = -4;
                  /* Copy over second block,  containing constants data: */
         else if( fread( buff, 1, eph->recsize, eph->ifile) != (size_t)eph->recsize)
            rval = -5;
         else if( fwrite( buff, 1, eph->recsize, ofile) != (size_t)eph->recsize)
            rval = -6;
                  /* Copy over the actual ephemeris records: */
         else if( fseek( eph->ifile, (start_block + 2L) * eph->recsize, SEEK_SET))
            rval = -7;
         }
      for( i = start_block; !rval && i < end_block; i++)
         {
         if( fread( buff, 1, eph->recsize, eph->ifile) != (size_t)eph->recsize)
            rval = -8;
         else if( fwrite( buff, 1, eph->recsize, ofile) != (size_t)eph->recsize)
            rval = -9;
         }
      fclose( ofile);
      eph->curr_cache_loc = -1;         /* mark the cache as 'dirty' */
      }
   return( rval);
}

#ifdef TEST_MAIN

#define J2000 2451545.

int main( const int argc, const char **argv)
{
   int rval = 0;

   if( argc < 5)
      {
      printf( "'sub_eph' takes a large binary DE ephemeris (such as the\n");
      printf( "files provided on the Willmann-Bell CD-ROM) and creates a\n");
      printf( "'sub-ephemeris' covering a specified range of years.  For\n");
      printf( "example,  to extract a small ephemeris 'sm_eph.406' from\n");
      printf( "the file 'unix.406' covering the years 1900 to 2100:\n\n");
      printf( "sub_eph unix.406 sm_eph.406 1900 2100\n");
      rval = -2;
      }
   else
      {
      void *ephem = jpl_init_ephemeris( argv[1], NULL, NULL);
      const double year1 = atof( argv[3]);
      const double year2 = atof( argv[4]);

      if( !ephem)
         printf( "Ephemeris '%s' not loaded\n", argv[1]);
      else
         {
         int err_val = make_sub_ephem( ephem, argv[2],
                                      J2000 + (year1 - 2000.) * 365.25,
                                      J2000 + (year2 - 2000.) * 365.25);

         if( err_val)
            printf( "Error %d occurred\n", err_val);
         }
      }
   return( rval);
}
#endif
