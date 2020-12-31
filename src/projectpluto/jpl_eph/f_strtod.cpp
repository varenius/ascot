/* f_strtod.cpp: "fast" version of strtod( ).

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

#include <stdint.h>

#ifndef INT64_MAX
   #ifdef _MSC_VER
      # define INT64_MAX      (9223372036854775807i64)
   #else
      # define INT64_MAX      (9223372036854775807LL)
   #endif
#endif

    /* I've found that most of the time in asc2eph,  the program   */
    /* to convert ASCII JPL ephemerides to binary,  was consumed   */
    /* in parsing floating-point values from the ASCII data.  That */
    /* part was done using sscanf().  strtod() and atof() proved   */
    /* faster,  but still quite slow.  The following was written   */
    /* with speed in mind;  as a result,  it's doubtless bigger    */
    /* and more complicated than your average library strtod().    */
    /* But it speeds asc2eph by two or three times,  depending on  */
    /* platform and compiler.                                      */
    /*    Also,  note that while it's been well-tested with the    */
    /* asc2eph program,  and somewhat tested with the TEST_CODE at */
    /* the bottom of this file,  it's entirely possible that odd   */
    /* cases may remain.  The TEST_CODE simplifies comparison with */
    /* results from the library strtod() and atof() functions.     */
    /*    2014 May 15:  streamlined fast_strtold()'s logic.        */

long double fast_strtold( const char *iptr, char **endptr)
{
   int n_decimal_places = 0;
   const int64_t max_rval = INT64_MAX / 10 - 1;
   const char *tptr;
   int64_t rval = 0;
   int exponent = 0;
   bool exponent_is_negative = false;
   bool is_negative = false;
   long double d_rval;

   while( *iptr == ' ')
      iptr++;
   if( *iptr == '+')
      iptr++;
   else if( *iptr == '-')
      {
      iptr++;
      is_negative = true;
      }
   while( rval < max_rval && *iptr >= '0' && *iptr <= '9')
      rval = rval * 10 + (int64_t)( *iptr++ - '0');
   if( rval >= max_rval)        /* meaningless digits provided; */
      {
      tptr = iptr;
      while( *iptr >= '0' && *iptr <= '9')    /* skip over them */
         iptr++;
      n_decimal_places = tptr - iptr;
      }
   if( *iptr == '.')
      {
      tptr = ++iptr;
      while( rval < max_rval && *iptr >= '0' && *iptr <= '9')
         rval = rval * 10 + (int64_t)( *iptr++ - '0');
      n_decimal_places += iptr - tptr;
      if( rval >= max_rval)        /* meaningless digits provided; */
         while( *iptr >= '0' && *iptr <= '9')    /* skip over them */
            iptr++;
      }

   if( *iptr == 'e' || *iptr == 'E')
      {
      iptr++;
      if( *iptr == '+')
         iptr++;
      else if( *iptr == '-')
         {
         exponent_is_negative = true;
         iptr++;
         }
      while( *iptr >= '0' && *iptr <= '9')
         exponent = exponent * 10 + (int)( *iptr++ - '0');
      }
   exponent += (exponent_is_negative ? n_decimal_places : -n_decimal_places);
   if( exponent < 0)
      {
      exponent = -exponent;
      exponent_is_negative = true;
      }
   d_rval = (long double)( is_negative ? -rval : rval);
   if( exponent >= 16)
      {
      int mask;
      long double ten_pow = (exponent_is_negative ? 1.0e-16 : 1.e+16);

      for( mask = 16; mask <= exponent; mask <<= 1, ten_pow *= ten_pow)
         if( exponent & mask)
            d_rval *= ten_pow;
      exponent %= 16;
      }
   if( exponent)
      {
      static const long double multipliers[16] = { 1., 10., 100., 1000.,
               10000., 100000., 1e+6, 1e+7, 1e+8, 1e+9, 1e+10, 1e+11,
               1e+12, 1e+13, 1e+14, 1e+15 };

      if( exponent_is_negative)
         d_rval /= multipliers[exponent];
      else
         d_rval *= multipliers[exponent];
      }
   if( endptr)
      *endptr = (char *)iptr;
   return( d_rval);
}

double fast_strtod( const char *iptr, char **endptr)
{
   return( (double)fast_strtold( iptr, endptr));
}

#ifdef TEST_CODE

#include <stdio.h>
#include <stdlib.h>

int main( const int argc, const char **argv)
{
   char *endptr, *endptr1;
   const double tval = fast_strtod( argv[1], &endptr);
   const double tval1 =     strtod( argv[1], &endptr1);
   const double tval2 = atof( argv[1]);

   printf( "new:  %.37lg %zd\n", tval, endptr - argv[1]);
   printf( "old:  %.37lg %zd\n", tval1, endptr1 - argv[1]);
   printf( "atof: %.37lg\n", tval2);
   if( tval != tval1)
      printf( "FAST_STRTOD DISAGREES WITH STRTOD!!!!\n");
   if( endptr != endptr1)
      printf( "FAST_STRTOD READ A DIFFERENT NUMBER OF BYTES THAN STRTOD!!!!\n");
   if( tval != tval2)
      printf( "FAST_STRTOD DISAGREES WITH ATOF!!!!\n");
   if( tval1 != tval2)
      printf( "STRTOD DISAGREES WITH ATOF!!!!\n");
   return( 0);
}
#endif
