/* merge_de.cpp: merge two or more JPL binary ephemeris files

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
#include <string.h>
#include "jpleph.h"

struct jpl_file
   {
   char header[84 * 3];
   const char *filename;
   long jd_start, jd_end;
   };

/* Code to combine an arbitrary number of JPL binary ephemerides.  Opens
first file to get the kernel and step size;  then it opens all the files
to be combined,  to check that they are all of the same DE number and to
get their time spans.  These are then sorted,  and checked to ensure that
they meet or overlap.  Finally,  the output file is created. */

#define MERGE_ERR_INITIALIZING                  -1
#define MERGE_ERR_OPENING                       -2
#define MERGE_ERR_READING_HEADER                -3
#define MERGE_ERR_NOT_A_JPL_EPHEM               -4
#define MERGE_ERR_MISMATCHED_VERSIONS           -5
#define MERGE_ERR_GAP_IN_DATA                   -6
#define MERGE_ERR_OPENING_2                     -7
#define MERGE_ERR_OUT_OF_MEMORY                 -8
#define MERGE_ERR_OPENING_3                     -9
#define MERGE_ERR_DATA_READ                     -10
#define MERGE_ERR_DATA_WRITE                    -11
#define MERGE_ERR_DATA_READ_2                   -12
#define MERGE_ERR_DATA_WRITE_2                  -13
#define MERGE_ERR_DATA_READ_3                   -14
#define MERGE_ERR_DATA_WRITE_3                  -15
#define MERGE_ERR_DATA_WRITE_HEADER             -16

int merge_jpl_files( const char *output_filename, const int n_input_files,
                     const char **input_filenames)
{
   void *jpl_eph = jpl_init_ephemeris( input_filenames[0], NULL, NULL);
   int i, j, de_number, kernel_size, kernel_days, rval = 0;
   FILE *ifile, *ofile;
   struct jpl_file *idata = (struct jpl_file *)calloc( n_input_files,
                                             sizeof( struct jpl_file));
   char *buff = NULL;

   if( !jpl_eph || !idata)
      return( MERGE_ERR_INITIALIZING);
   de_number = (int)jpl_get_long( jpl_eph, JPL_EPHEM_EPHEMERIS_VERSION);
   kernel_days = (int)( jpl_get_double( jpl_eph, JPL_EPHEM_STEP) + .5);
   kernel_size = jpl_get_long( jpl_eph, JPL_EPHEM_KERNEL_SIZE);
   jpl_close_ephemeris( jpl_eph);

   for( i = 0; !rval && i < n_input_files; i++)
      {
      ifile = fopen( input_filenames[i], "rb");
      if( !ifile)
         rval = MERGE_ERR_OPENING;
      else
         {
         if( fread( idata[i].header, 84, 3, ifile) != 3)
            rval = MERGE_ERR_READING_HEADER;
         fclose( ifile);
         if( !rval)
            {
            if( memcmp( idata[i].header, "JPL Planetary Ephemeris", 23))
               rval = MERGE_ERR_NOT_A_JPL_EPHEM;
            else if( de_number != atoi( idata[i].header + 26))
               {
               printf( "You can't merge files from different DE ephemerides\n");
               printf( "'%s' is from DE-%d\n", input_filenames[i],
                                                atoi( idata[i].header + 26));
               printf( "'%s' is from DE-%d\n", input_filenames[0], de_number);
               rval = MERGE_ERR_MISMATCHED_VERSIONS;
               }
            }
         idata[i].jd_start = atol( idata[i].header + 102);
         idata[i].jd_end   = atol( idata[i].header + 186);
         idata[i].filename = input_filenames[i];
         }
      }

   if( rval)
      return( rval);
                   /* OK, now sort by date: */
   for( i = 0; i < n_input_files; i++)
      for( j = 0; j < i; j++)
         if( idata[i].jd_start < idata[j].jd_start)
            {
            struct jpl_file temp = idata[i];

            idata[i] = idata[j];
            idata[j] = temp;
            }
   printf( "Merging:\nJD Start  JD End     Filename\n");
   for( i = 0; i < n_input_files; i++)
      printf( "%8ld.5 %8ld.5 %s\n", idata[i].jd_start, idata[i].jd_end,
                                                       idata[i].filename);
                  /* Check to be sure the files overlap: */
   for( i = 0; i < n_input_files - 1; i++)
      if( idata[i + 1].jd_start > idata[i].jd_end)
         {
         printf( "ERROR:  there is a gap between '%s' and '%s'.\n",
                        idata[i].filename, idata[i + 1].filename);
         rval = MERGE_ERR_GAP_IN_DATA;
         }

   if( !rval)
      {
      ofile = fopen( output_filename, "wb");
      if( !ofile)
         {
         printf( "Couldn't open output file '%s'\n", output_filename);
         rval = MERGE_ERR_OPENING_2;
         }
      }
   if( !rval)
      {
      buff = (char *)malloc( kernel_size);
      if( !buff)
         rval = MERGE_ERR_OUT_OF_MEMORY;
      }
   for( i = 0; !rval && i < n_input_files; i++)
      {
      int n_blocks;

      ifile = fopen( idata[i].filename, "rb");
      if( !ifile)
         rval = MERGE_ERR_OPENING_3;
      if( !i)
         for( j = 0; !rval && j < 2; j++)       /* copy out two header blocks */
            {
            if( fread( buff, kernel_size, 1, ifile) != 1)
               rval = MERGE_ERR_DATA_READ;
            else if( fwrite( buff, kernel_size, 1, ofile) != 1)
               rval = MERGE_ERR_DATA_WRITE;
            }
         else
            fseek( ifile, 2L * kernel_size, SEEK_SET);
      if( i == n_input_files - 1)
         n_blocks = (idata[i].jd_end - idata[i].jd_start) / kernel_days;
      else
         n_blocks = (idata[i + 1].jd_start - idata[i].jd_start) / kernel_days;
      while( !rval && n_blocks--)
         {
         if( fread( buff, kernel_size, 1, ifile) != 1)
            rval = MERGE_ERR_DATA_READ_2;
         else if( fwrite( buff, kernel_size, 1, ofile) != 1)
            rval = MERGE_ERR_DATA_WRITE_2;
         }
      if( !rval && i == n_input_files - 1)      /* patch up ending date */
         {
         fseek( ofile, 2L * 84L, SEEK_SET);
         if( fwrite( idata[i].header + 2 * 84, 84, 1, ofile) != 1)
            rval = MERGE_ERR_DATA_WRITE_HEADER;
         fseek( ifile, 2660L, SEEK_SET);
         if( !rval)
            {
            if( fread( buff, 8, 1, ifile) != 1)
               rval = MERGE_ERR_DATA_READ_3;
            fseek( ofile, 2660L, SEEK_SET);
            if( !rval && fwrite( buff, 8, 1, ofile) != 1)
               rval = MERGE_ERR_DATA_WRITE_3;
            }
         }
      fclose( ifile);
      }
   // free memory and close filehandles
   if( buff)
      free( buff);
   free( idata);
   fclose( ofile);
   return( rval);
}

int main( int argc, const char **argv)
{
   int err_code = -99;

   if( argc < 4)
      {
      printf( "'merge_de' takes as command-line arguments the name of the\n");
      printf( "output JPL DE file to be created,  followed by the names of\n");
      printf( "the input files.  For example,  to merge DE1990.406, DE2000.406,\n");
      printf( "and DE2010.406 into the single file 30YEARS.406:\n\n");
      printf( "merge_de 30years.406 de1990.406 de2000.406 de2010.406\n");
      }
   else
      {
      err_code = merge_jpl_files( argv[1], argc - 2, argv + 2);
      if( err_code)
         printf( "Error code %d\n", err_code);
      }
   return( err_code);
}
