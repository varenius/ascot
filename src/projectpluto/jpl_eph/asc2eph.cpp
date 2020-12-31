/* asc2eph.cpp: convert ASCII JPL ephemerides to binary

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

/***************************************************************************
***************             ASC2EPH.C v.1.2                *****************
****************************************************************************
** This program is a C version of an original JPL fortran code, slightly  **
** changed (interactive input added) for user convenience.                **
**                                                                        **
** You will be prompted for start and final epochs (in JED) of the        **
** ephemeris you want to obtain in binary form.                           **
** This epochs will be rounded to the nearest integer multiply of 32 days **
** counting from 2440400.5, outside of the interval given.                **
****************************************************************************
**  Written: May 28, 1997 by PAD   **  Last modified: June 23,1997 by PAD **
****************************************************************************
**  PAD: dr. Piotr A. Dybczynski,          e-mail: dybol@phys.amu.edu.pl  **
**   Astronomical Observatory of the A.Mickiewicz Univ., Poznan, Poland   **
****************************************************************************

   Since heavily modified by Bill J. Gray.  The current version can be used
for all DE versions without recompiling.  The program figures out which
DE files are available,  and there are command line options to specify
the path where the ASCII files are to be found;  the output file name;
and the range of years to be extracted,  as well as a 'verbose mode'.  I
found that feof( ) apparently didn't work as expected in MinGW;  a bit
of rearrangement to produce better error detection fixed this.

   I also revised it to make it more 'future proof'.  As more DE ephemerides
come out,  I think this will adjust automatically to changes in kernel size,
number of coefficients,  date ranges,  and so on.  */

/* 2011 Oct 4:  prompted by an inquiry from Tatjana Jaksic,  I revised
the code to work on 64-bit Linux.  This involved providing a #define for
_MAX_PATH,  changing the path separator from '\' to '/',  using int32_t
in some structures (instead of just relying on integers being 32 bits by
default),  and cleaning up a _lot_ of g++ warnings.  Also,  there was a
lowercase/uppercase filename issue that wasn't an issue in Windows.  In
the process of doing all this,  I added better error checking,  a
progress bar (dots across screen until done),  and fixed a performance
problem that caused slowdowns when asking for a subset of the entire
ephemeris.  */

/* 2011 Oct 7:  further testing revealed that almost all the processing
time was spent in sscanf() extracting double-precision floats from ASCII
data.  strtod() sped things up,  but I had to write a new fast_strtod()
function to get real speed (see f_strtod.cpp for details). */

/* 2013 Apr 6:  DE-430 has 572 constants.  Revising the code to handle
the actual constants wasn't so bad;  there is actually room in the record
for about 1018 constants (for DE-4xx).  The names of those constants get
a little trickier.  The binary header structure allows for 400 names.
However,  there's zero padding after that header structure allowing
for an additional 881 names (for DE-4xx).  The code now does that sort
of storage... eventually,  we'll hit the 1018 limit,  work around it,
and then hit the 1281 limit.  */

/* 2013 Aug 10:  DE-431 uses five-digit years,  because it runs from
-13200 to +17191,  so I revised 'get_ascii_de_file' to use the revised
file names.  The 1050 group is a bit different as well,  with extra
numbers at the end of all three lines,  presumably to accommodate some
future features;  that required some small changes.  To be done:  handling
the extra range (code currently works only from -4000 to +6000) and
probably some negative JD issues. */

/* 2013 Aug 17:  Comparison of the binary DE-431 ephems created with
this software to those on the JPL site revealed that in the first 84*3
bytes,  the lines are space-padded,  not zero-padded.  This software now
does the same thing. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#ifdef __linux
   #include <limits.h>
   #define _MAX_PATH   PATH_MAX
#endif

static void nxtgrp(char *header, FILE *ifile);
static void errprt( const int i, const char * msg);

#define NRECL 4
#define JPL_MAX_N_CONSTANTS             400

#pragma pack(1)

struct header_record {
         char ttl[3][84];
         char cnam[JPL_MAX_N_CONSTANTS][6];
         double ss[3];
         int32_t ncon;
         double au;
         double emrat;        /* Earth-Moon mass ratio */
         int32_t ipt[12][3];  /* Pointers to # coeffs for objects */
         int32_t numde;
         int32_t lpt[3];      /* Pointers to # coeffs for lunar librations */
       };

#pragma pack()

/* ASCII DE files have the form 'ascp2150.423',  for example,  for the
DE423 file for the year +2150;  or 'ascm0100.405',  for example,  for the
DE405 file for the year -100.  Except that DE431,  DE430t,  and I'll assume
subsequent ephemerides,  use five-digit years:  e.g., 'ascp02000.431'.
So check 'em both and see which one works.   */

static FILE *get_ascii_de_file(
            const char *path_to_ascii_files, const int year, const char *de_num)
{
   unsigned i;
   FILE *rval = NULL;

   for( i = 0; !rval && i < 2; i++)
      {
      char buff[_MAX_PATH];
      const char *format_string = (i ? "%sasc%1c%05d.%s":
                                       "%sasc%1c%04d.%s");

      sprintf( buff, format_string, path_to_ascii_files,
                 (year < 0 ? 'm' : 'p'), abs( year), de_num);
      rval = fopen( buff, "rb");
      }
   return( rval);
}

/* In determining which ASCII DE files are available,  the following logic
is used:  we start looking for files at the year -4000,  admittedly outside
the current range of any DE ephemeris (except DE-431,  and possibly DE-408,
and there may be future exceptions).  We keep looking at five-year intervals
by default.  When we find our first file,  we know what *year_start is.
When we find our second,  we know the step size (by subtracting the year
for the previously found file).

   To handle ephemerides that extend before -4000,  we then start counting
backward.  This does mean that binary ephems cannot contain only data for
dates before -4000,  but I don't consider that to be a really major
limitation.  Further revision may be needed to handle ephems that run into
the distant future.  */

static int determine_year_range( int *year_start, int *year_end,
             int *year_step, const char *path_to_ascii_files, const char *de_num)
{
   int year;
   const int max_year = (atoi( de_num) == 431 ? 19000 : 6000);
   unsigned n_found = 0;
   FILE *ifile;

   *year_step = 5;
   *year_start = *year_end = 0;   /* this just suppresses a compiler  */
                     /* warning about using these variables uninitialized. */
                     /* gcc isn't quite bright enough to realize that the  */
                     /* warning is a spurious one.                         */
   for( year = -4000; year < max_year; year += *year_step)
      {
      ifile = get_ascii_de_file( path_to_ascii_files, year, de_num);
      if( ifile)
         {
         fclose( ifile);
         if( !n_found)
            *year_start = year;
         if( n_found == 1)
            *year_step = year - *year_start;
         n_found++;
         *year_end = year;
         }
#ifdef DEBUGGING_STATEMENTS
      if( n_found)
         printf( "%d found; year %d, ifile %p, %d to %d (step %d)\n",
               n_found, year, ifile, *year_start, *year_end, *year_step);
#endif
      }
                  /* OK,  now start looking backward: */
   if( n_found)
      while( (ifile = get_ascii_de_file( path_to_ascii_files,
                           *year_start - *year_step, de_num)) != NULL)
         {
         *year_start -= *year_step;
         fclose( ifile);
         }
   return( (n_found ? 0 : -1));
}


static void error_msg( void)
{
   printf( "ASC2EPH takes as command-line arguments the name of the path\n");
   printf( "to the JPL ASCII files (optionally;  default is the current path)\n");
   printf( "followed by any of the following options:\n\n");
   printf( "   -d(number)     Specifies the DE number (200,  405,  406,  etc.)\n");
   printf( "   -o(filename)   Specifies the output filename (default is jpleph.xxx)\n");
   printf( "   -r(JD1,JD2)    (or year1, year2) Specifies output range.\n");
   printf( "   -v             Verbose mode\n");
   printf( "   -v2            Extra-verbose mode\n");
}

static void update_progress_bar( const double fraction_done)
{
   static int dots_shown = 0;
   int dots_to_show = (int)( fraction_done * 79.);

   while( dots_shown < 79 && dots_shown < dots_to_show)
      {
      printf( ".");
      dots_shown++;
      }
}

  /* Both the header files and the actual data files contain lines    */
  /* with three floating-point values,  stored using D to indicate     */
  /* the exponent (instead of the E expected in C/C++).  This code      */
  /* replaces Ds with Es,  then extracts the three values.               */
  /*    Previously,  they were extracted with sscanf:                     */
  /*                                                                      */
  /* return( sscanf( iline, "%lE %lE %lE", ovals, ovals + 1, ovals + 2)); */
  /*                                                                     */
  /*   Oddly,  this made things run at about half the speed as the      */
  /* strtod() version,  which in turn is two or three times slower     */
  /* than the fast_strtod() version!  (See f_strtod() for details.)   */
  /* This was true for various compilers under Windows and Linux,    */
  /* and on 32-bit and 64-bit systems.                              */

double fast_strtod( const char *iptr, char **endptr);

static int get_three_doubles( char *iline, double *ovals)
{
   int i;

   for( i = 0; i < 102 && iline[i]; i++)
      if( iline[i] == 'D')
         iline[i] = 'E';

#ifdef HORRIBLY_SLOW
   return( sscanf( iline, "%lE %lE %lE", ovals, ovals + 1, ovals + 2));
#else
   for( i = 0; iline && i < 3; i++)
      {
      char *endptr;
//    const double tval = strtod( iline, &endptr);
      const double tval = fast_strtod( iline, &endptr);

      if( endptr != iline)  /* yes,  a double was read */
         {
         ovals[i] = tval;
         iline = (char *)endptr;
         }
      else                  /* no new double was read; */
         iline = NULL;      /* break out of loop       */
      }
   return( i);
#endif
}

#define J2000 2451545.0
#define YEAR_TO_JD( year)   (J2000 + ((year)-2000.) * 365.25)
#define JD_TO_YEAR( jd)     (((jd)-J2000) / 365.25 + 2000.)

int main( const int argc, const char **argv)
{
    char header[14];
    char buff[102];
    char path_to_ascii_files[_MAX_PATH];
    char output_filename[_MAX_PATH];
    double jd1 = -99999999., jd2 = 99999999., db2z, *db;
    unsigned i, j, ksize;
    const char *de_num = "405";
    int year;
    unsigned n, nrout, ncoeff, nrw, out, last;
    unsigned n_padding_zeroes;
    int year_start, year_end, year_step;
    int verbose = 0;
    FILE *ifile, *ofile;
    struct header_record rec1;
    double *cval;
    char *zero_padding_buffer;
    char *cnames;
    int32_t rpt[3];      /* Pointers to # coeffs for lunar Euler angle rates */
    int32_t tpt[3];      /* Pointers to # coeffs for TT-TDB */

    memset( &rec1, 0, sizeof( rec1));
    if( argc < 2)
      {
      error_msg( );
      return( -1);
      }
    setvbuf( stdout, NULL, _IONBF, 0);
    *output_filename = *path_to_ascii_files = '\0';
    for( i = 1; i < (unsigned)argc; i++)
       if( argv[i][0] == '-')
          switch( argv[i][1])
             {
             case 'd':
                de_num = argv[i] + 2;
                break;
             case 'r':
                sscanf( argv[i] + 2, "%lf,%lf", &jd1, &jd2);
                if( jd1 < 30000.)         /* must be a year;  cvt to a JD */
                   jd1 = YEAR_TO_JD( jd1);
                if( jd2 < 30000.)         /* must be a year;  cvt to a JD */
                   jd2 = YEAR_TO_JD( jd2);
                break;
             case 'o':
                strcpy( output_filename, argv[i] + 2);
                j = strlen( argv[i]);
                while( j && argv[i][j - 1] != '.')
                  j--;
                if( j)
                  de_num = argv[i] + j;
                break;
             case 'v':
                verbose = atoi( argv[i] + 2);
                if( !verbose)
                   verbose = 1;
                break;
             default:
                error_msg( );
                return( -2);
             }
        else
           strcpy( path_to_ascii_files, argv[i]);

/****************************************************************************/
/* You can add path before filename (if necessary) in the following line:   */
   if( *path_to_ascii_files)
#if defined( __GNUC__) && !defined( __MINGW32__)
      strcat( path_to_ascii_files, "/");
#else
      strcat( path_to_ascii_files, "\\");
#endif
   sprintf( buff, "%sheader.%s", path_to_ascii_files, de_num);
/****************************************************************************/
   ifile = fopen( buff, "rb");
   if( !ifile)
      {
      printf("Cannot open header file: %s, aborted.\n",buff);
      return( -4);
      }

/*  write a fingerprint to the screen.  */

   if( verbose)
      puts("\n JPL ASCII-TO-DIRECT-I/O PROGRAM.\n C-version, translated from original fortran code from JPL\n");

/*  read the size and number of main ephemeris records. */

   if( fscanf( ifile,"KSIZE= %d NCOEFF= %d", &ksize, &ncoeff) != 2)
      errprt( 1000, "KSIZE/NCOEFF fail\n");
   if( verbose)
      printf( "KSIZE = %6d  NCOEFF= %6d", ksize, ncoeff);
   db = (double *)malloc( (ncoeff + 3) * sizeof( double));
   zero_padding_buffer = (char *)calloc( ksize, NRECL);
   if( !db || !zero_padding_buffer)
      errprt( 1000, "Out of memory");

/*  now for the alphameric heading records (group 1010) */

   nxtgrp( header, ifile);
   if( strcmp(header,"GROUP   1010"))
      errprt(1010,"NOT HEADER");

   for( i = 0; i < 3; i++)
      {
      if( !fgets( buff, sizeof( buff), ifile))
         errprt( 1010, "fgets fail\n");
      for( j = 0; j < sizeof( rec1.ttl[i]) && buff[j] >= ' '; j++)
         rec1.ttl[i][j] = buff[j];
      while( j < sizeof( rec1.ttl[i]))
         rec1.ttl[i][j++] = ' ';
      if( verbose)
         puts( buff);
      }

/*  read start, end and record span  (group 1030) */

   nxtgrp( header, ifile);
   if( strcmp( header, "GROUP   1030"))
      errprt( 1030, "NOT HEADER");

   if( fscanf( ifile, " %lf %lf %lf", &rec1.ss[0], &rec1.ss[1], &rec1.ss[2]) != 3)
      errprt( 1030, "sscanf fail\n");
   /* There is an error in DE403 header file: */
   if( atoi( de_num) == 403)
      rec1.ss[0] = 2305424.5;

/* read number of constants and names of constants (group 1040/4). */

   nxtgrp( header, ifile);
   if( strcmp( header, "GROUP   1040"))
      errprt( 1040, "NOT HEADER");

   if( fscanf( ifile, " %d", &n) != 1)
      errprt( 1040, "fscanf fail\n");
   if( verbose > 1)
      printf( "Reading %d constants\n", n);
   cnames = (char *)calloc( n, 6);
   for( i = 0; i < n; i++)
      if( fscanf( ifile, " %6c", buff) != 1)
         errprt( 1040, "fscanf fail (2)\n");
      else
         {
         memcpy( cnames + i * 6, buff, 6);
         if( i < JPL_MAX_N_CONSTANTS)
            memcpy( rec1.cnam[i], buff, 6);
         }
   if( verbose > 1)
      printf( "%d constants read\n", n);
   rec1.ncon = n;
            /* Most software written before DE-430 has a hard-coded    */
            /* limit of 400 ephemeris constants.  Some software may    */
            /* object to this;  if so,  uncomment the following line   */
            /* to 'ncon' at 400.                                       */
// rec1.ncon = (n > JPL_MAX_N_CONSTANTS ? JPL_MAX_N_CONSTANTS : n);

/*  read number of values and values (group 1041/4)  */
   nxtgrp( header, ifile);
   if( strcmp( header, "GROUP   1041"))
      errprt( 1041, "NOT HEADER");

   if( fscanf( ifile, " %d", &n) != 1 || !fgets( buff, 100, ifile))
      errprt( 1041, "File error");
   i = 0;
   cval = (double *)calloc( n, sizeof( double));
   while( i < n)
      {
      double temp[3];
      unsigned n_found;

      if( !fgets( buff, 100, ifile))
         errprt( 1041, "fgets error");
      n_found = get_three_doubles( buff, temp);
      if( n_found > n - i)
         n_found = n - i;
      memcpy( cval + i, temp, n_found * sizeof( double));
      i += n_found;
      }

   for( i = 0; i < n && i < JPL_MAX_N_CONSTANTS; i++)
      {
      for( j = 0; j < 6; j++)
         buff[j]=rec1.cnam[i][j];
      buff[6]='\0';
      if( !strcmp( buff, "AU    "))
         rec1.au    = cval[i];
      if( !strcmp( buff, "EMRAT "))
         rec1.emrat = cval[i];
      if( !strcmp( buff, "DENUM "))
         rec1.numde = (long)( cval[i] + .5);
      }

   if( verbose > 1)
      {
      const unsigned n_lines = n / 2 + (n & 1);

      for( i = 0; i < n_lines; i++)
         {
         printf("%.6s  %24.16E  ", cnames + i * 6, cval[i]);
         if( i + n_lines < n)
            printf("    %.6s  %24.16E\n",
                   cnames + (i + n_lines) * 6, cval[i + n_lines]);
         else
            printf( "\n");
         }
      }

/*  read pointers needed by interp (group 1050)  */
   nxtgrp( header, ifile);
   if( strcmp( header, "GROUP   1050"))
      errprt( 1050, "NOT HEADER");

   memset( buff, 0, sizeof( buff));
   for( i = 0; i < 3; i++)
      if( !fgets( buff, 100, ifile))
         errprt( 1051, "fgets error");
      else
         {
         for( j = 0; j < 12; j++)
            rec1.ipt[j][i] = (int32_t)atoi( buff + j * 6);
         rec1.lpt[i] = (int32_t)atoi( buff + 72);
         rpt[i] = (int32_t)atoi( buff + 78);
         tpt[i] = (int32_t)atoi( buff + 84);
         }

   if( verbose)
      for( i = 0; i < 3; i++)
         {
         for( j = 0; j < 12; j++)
            printf( " %5d", (int)rec1.ipt[j][i]);
         printf( " %5d %5d %5d\n", (int)rec1.lpt[i],
                        (int)rpt[i], (int)tpt[i]);
         }

            /* If we reset the range of the output ephemeris,  we'll need to */
            /* reset those values in the header,  too: */
   if( jd1 > rec1.ss[0])
      {
      jd1 -= 2440400.5;
      jd1 = rec1.ss[2] * floor( jd1 / rec1.ss[2]) + 2440400.5;
      rec1.ss[0] = jd1;
      }
   else
      jd1 = rec1.ss[0];

   if(jd2 < rec1.ss[1])
      {
      jd2 -= 2440400.5;
      jd2 = rec1.ss[2] * ceil( jd2 / rec1.ss[2]) + 2440400.5;
      rec1.ss[1] = jd2;
      }
   else
      jd2 = rec1.ss[1];

/*   open direct-access output file (defaults to 'jpleph.xxx') */

   if( !*output_filename)
      sprintf( output_filename, "jpleph.%s", de_num);
/***************************************************************************/

   ofile=fopen( output_filename, "wb");
   if( !ofile)
      {
      printf("Cannot create binary output file: %s,aborted.\n",
                                                output_filename);
      fclose( ifile);
      return( -3);
      }


/*  write header records onto output file.   */

   out = fwrite( &rec1, sizeof( rec1), 1,ofile);
   if( out != 1)
      errprt( 1, "ST RECORD NOT WRITTEN BECAUSE OF ERROR\n");
   n_padding_zeroes = ksize * NRECL - sizeof( rec1);
   assert( ksize * NRECL >= sizeof( rec1));
            /* The 'rec1' struct gives us room for 400 constant names.     */
            /* Extras can be stored in the area that would normally be     */
            /* zero padding.  I've included an assert in case we go too    */
            /* far... shouldn't happen until we have 1281 constants or so. */
            /* Could happen,  though,  so this assert should stay :        */
   if( n > JPL_MAX_N_CONSTANTS)
      {
      const int write_size = (n - JPL_MAX_N_CONSTANTS) * 6;

      out = fwrite( cnames + JPL_MAX_N_CONSTANTS * 6, write_size, 1, ofile);
      if( out != 1)
         errprt( 1, "ND RECORD: extra constant names write error\n");
      n_padding_zeroes -= write_size;
      }
   out  = fwrite( rpt, 1, sizeof( rpt), ofile);
   out += fwrite( tpt, 1, sizeof( tpt), ofile);
   if( out != sizeof( rpt) + sizeof( tpt))
      errprt( 1, "ND RECORD: write error on rpt/tpt\n");
   n_padding_zeroes -= sizeof( rpt) + sizeof( tpt);
   out = fwrite( zero_padding_buffer, n_padding_zeroes, 1, ofile);
   if( out != 1)
      errprt( 1, "ND RECORD: padding write error\n");

            /* There _should_ be room for about a thousand constants. If */
            /* they wouldn't all fit into ksize * NRECL bytes,  the      */
            /* following lines will chop n down to size.  This won't     */
            /* happen until we get to 1018 constants (for DE-4xx).  I've */
            /* added an assert for this case.  It's important to know    */
            /* about it (rather than just silently truncate constants    */
            /* after the 1018th) because other code is not set up to     */
            /* handle it.                                                */
   assert( n * sizeof( double) <= ksize * NRECL);
   if( n * sizeof( double) > ksize * NRECL)
      n = (ksize * NRECL) / sizeof( double);
   assert( ksize * NRECL > n * sizeof( double));
   n_padding_zeroes = ksize * NRECL - n * sizeof( double);
   out = fwrite( cval, n * sizeof( double), 1, ofile);
   if( out != 1)
      errprt( 2, "ND RECORD: cval write error\n");
   out = fwrite( zero_padding_buffer, n_padding_zeroes, 1, ofile);
   if( out != 1)
      errprt( 2, "ND RECORD: padding write error\n");
// printf( "%u after rec1; %u after cval\n",
//                ksize * NRECL - sizeof( rec1),
//                ksize * NRECL - n * sizeof( double));

/*  read and write the ephemeris data records (group 1070). */
   nxtgrp( header, ifile);
   fclose( ifile);
   if( strcmp( header, "GROUP   1070"))
      errprt( 1070, "NOT HEADER");

   nrout  = 0;
   out    = 0;
   db2z   = 0; /* formal only */

   if( determine_year_range( &year_start, &year_end, &year_step,
                                      path_to_ascii_files, de_num))
      {
      printf( "Unable to locate files in path '%s' for DE-%s\n",
                              path_to_ascii_files, de_num);
      exit( -1);
      }
   if( verbose)
      printf( "Years %d to %d found;  step size %d\n",
               year_start, year_end, year_step);
            /* If we're making a subset,  we may be able to skip a   */
            /* lot of the ASCII files.  Three years of "slop" are    */
            /* allowed for,  doubtless more than is really needed.   */
   year = (int)JD_TO_YEAR( jd1);
   while( year_start + year_step < year - 3)
      year_start += year_step;

   year = (int)JD_TO_YEAR( jd2);
   if( year_end > year + 3)
      year_end = year + 3;

/* main loop */
   for( year = year_start; year <= year_end; year += year_step)
      {
      ifile = get_ascii_de_file( path_to_ascii_files, year, de_num);
      if( verbose)
         printf( "File for year %d %s\n", year, (ifile ? "found" : "missing"));

      while( ifile && fgets( buff, 100, ifile))
         {
         if( sscanf( buff, " %d %d ",&nrw,&ncoeff) != 2)
             {
             printf( "Failed to read record number and number of coeffs!\n");
             printf( "Input line: '%s'\n", buff);
             printf( "Prev db: %lf %lf %lf\n", db[0], db[1], db[2]);
             errprt( nrw, "BAD END OF FILE");
             }
         if( 2 * ncoeff != ksize )
            errprt( ncoeff, "2*NCOEFF NOT EQUAL TO KSIZE");
/* there is an integer multiply of 3 coefficients in source file */
         last = ncoeff / 3 + 1;
         for( j = 0; j < last; j++)
            {
            if( !fgets( buff, 100, ifile))
               {
               printf( "Failed to read line %d!\n", j);
               errprt( nrw, "BAD END OF FILE");
               }
            if( !j || db[1] > jd1)
               get_three_doubles( buff, db + 3 * j);
#ifdef SHOW_WRONG_DATES
            if( !j)
               if( db[0] != floor( db[0]) + .5 || db[1] != floor( db[1]) + .5)
                  printf( "Bad dates:\n%s", buff);
#endif
            }

/*  The files initially posted for DE-431 had a problem in which the
    start/end time were sometimes misprinted;  JD 961967.5 could be shown
    as 9.61967500000001D+06 or 9.61967499999999D+06.  This happened for
    occasional dates between JD -1e+6 and JD +1e+6.  You can avoid the
    problem entirely by using the 2013 Aug 15 or later versions of the
    ASCII files,  in which the error is fixed.  Or you can just rely on
    these two lines to fix the problem :        */

    db[0] = floor( db[0]) + .5;
    db[1] = floor( db[1]) + .5;

/*  skip this data block if the end of the interval is less
    than the specified start time or if the it does not begin
    where the previous block ended.  */

         if( db[0] >= jd2)
            {
            year = year_end + 1;
            break;
            }
                                   /* all necessary blocks were processed */
         if( db[1] <= jd1)         /* wait for the first useful block     */
            continue;
         if( db[1] == db2z)        /* some blocks appears twice on the    */
            continue;              /* source files boundary !             */
         if( nrout && db[0] != db2z)
            {
                 /*  beginning of current interval is past the end
                     of the previous one.  */
            printf( "Jump from JD %lf to %lf\n", db[0], db2z);
            errprt(nrw,"RECORDS DO NOT OVERLAP OR ABUT");
            }

         db2z = db[1];
         if( !nrout)
            rec1.ss[0] = db[0];
         nrout++;
         out = fwrite( db, sizeof( double), ncoeff, ofile);
         if( out != ncoeff)
            errprt(nrout,"TH RECORD NOT WRITTEN BECAUSE OF ERROR");

/*  update the user as to our progress every 50th block.  */
         if( verbose && nrout % 50 == 1 )
            {
            if( db[0] >= jd1 )
               printf( "%d EPHEMERIS RECORDS WRITTEN. LAST JED= %12.2lf (%.2lf)",
                              nrout, db[1], JD_TO_YEAR( db[1]));
            else
               puts( "SEARCHING FOR THE FIRST REQUESTED RECORD...");
            printf( verbose > 1 ? "\n" : "\r");
            }
         if( !verbose)
            {
            const double fraction_done = (db[0] - jd1) / (jd2 - jd1);

            update_progress_bar( fraction_done);
            }
         }
      if( ifile)
         fclose( ifile);
      }
   if( verbose)
      printf( "\nLast rec: %d EPHEMERIS RECORDS WRITTEN. JED = %12.2f (%.2lf)\n",
           nrout, rec1.ss[1], JD_TO_YEAR( rec1.ss[1]));
   else        /* in non-verbose mode,  just end the progress bar: */
      printf( "\n");

   rec1.ss[1] = db2z;

 /* overwrite first record with epochs adjusted */
   fflush( ofile);
   fseek( ofile, 0L, SEEK_SET);
   fwrite( &rec1, sizeof( rec1), 1, ofile);
   free( db);

/*   we're through.  wrap it up.  */
   free( cval);
   free( cnames);

   fclose( ofile);
   printf( "%.2lf seconds elapsed\n",
                  (double)clock( ) / (double)CLOCKS_PER_SEC);
   return( 0);
}
/***************************************************************************/

static void errprt( const int i, const char * msg)
{
   printf( "\nERROR #%8d  %-50s", i, msg);
   exit( 1);
}

/****************************************************************************/
static void nxtgrp(char * header, FILE *ifile)
{
   char buff[102];

   while( fgets( buff, 100, ifile))
      if( buff[0] != ' ' && buff[0] != '\n')
         {
         buff[12]='\0';
         strcpy( header, buff);
         if( !fgets( buff, 100, ifile))
            *buff = '\0';
         return;
         }
}
/***************************  THE  END  ************************************/

