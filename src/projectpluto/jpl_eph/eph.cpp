#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"
#include "jpleph.h"

#define PI 3.141592653589793238462643383279
#define J2000_OBLIQUITY  (23.4392911 * PI / 180)

/* 12 Dec 1998:  Following an e-mail from Luc Desamore,  I 'offset' the */
/* time used in the ELP function by: */

/*  -0.000091 (n + 26)(year-1955)^2 seconds, n being here -23.8946.  */

/* This shift (in seconds) corresponds to a difference in Delta-T between */
/* that used for VSOP,  etc.,  and that used for ELP-82.                  */

static void show_state_vector( const double *state_vect)
{
   double r2 = state_vect[0] * state_vect[0] +
               state_vect[1] * state_vect[1] + state_vect[2] * state_vect[2];

   printf( "%.11lf %.11lf %.11lf  %.11lf\n", state_vect[0], state_vect[1],
            state_vect[2], sqrt( r2));
   printf( "%.11lf %.11lf %.11lf\n", state_vect[3], state_vect[4],
            state_vect[5]);
}

void main( int argc, char **argv)
{
   double state_vect[6], state_vect2[6], t0 = atof( argv[2]);
   int planet_no = atoi( argv[1]), i;
   void *p;
   char *filename;

   printf( "Year approx %.3lf\n", 2000. + (t0 - 2451545.) / 365.25);
   if( planet_no != 10)
      {
      FILE *ifile = fopen( "c:\\find_orb\\ps_1996.dat", "rb");

      if( !ifile)
         {
         printf( "ps_1996.dat not loaded\n");
         exit( -1);
         }
      p = load_ps1996_series( ifile, t0, planet_no);
      if( !p)
         {
         printf( "Series not loaded\n");
         exit( -1);
         }
      printf( "Series loaded\n");
      fclose( ifile);
      get_ps1996_position( t0, p, state_vect, 1);
      free( p);
      }
   else           /* lunar case */
      {
      FILE *ifile = fopen( "elp82big.dat", "rb");

      if( !ifile)
         {
         printf( "elp82big.dat not loaded\n");
         exit( -1);
         }
      compute_elp_xyz( ifile, (t0 - 2451545.0) / 36525., 0., state_vect);
      rotate_vector( state_vect, J2000_OBLIQUITY, 0);
      fclose( ifile);
      }

   show_state_vector( state_vect);
   rotate_vector( state_vect, -J2000_OBLIQUITY, 0);
   printf( "In ecliptic coords:\n");
   show_state_vector( state_vect);
   rotate_vector( state_vect,  J2000_OBLIQUITY, 0);

   filename = "..\\unxp2000.403";
   if( argc > 3)
      switch( atoi( argv[3]))
         {
         case 400:
            filename = "d:\\guide_b\\jpl_eph\\sub_de.406";
            break;
         case 406:
            filename = "h:\\unix.406";
            break;
         case 200:
            filename = "h:\\unix.200";
            break;
         default:
            filename = argv[3];
            break;
         }
   p = jpl_init_ephemeris( filename, NULL, NULL);
   if( !p)
      {
      printf( "JPL data not loaded\n");
      exit( -1);
      }
   if( planet_no != 10)
      jpl_pleph( p, t0, (planet_no == 3) ? 13 : planet_no, 11, state_vect2, 1);
   else
      {
//    const double t_minus_1955 = (t0 - 2435108.5) / 365.25;
//    const double t_elp = t0 - .000091 * (26 - 23.8946)
//                      * t_minus_1955 * t_minus_1955 / 86400.;

//    jpl_pleph( p, t_elp, 10, 3, state_vect2, 1);
      jpl_pleph( p, t0, 10, 3, state_vect2, 1);
      for( i = 0; i < 6; i++)
         state_vect2[i] *= AU_IN_KM;
      }
   jpl_close_ephemeris( p);

   printf( "DE state vector:\n");
   show_state_vector( state_vect2);
   for( i = 0; i < 6; i++)
      state_vect2[i] -= state_vect[i];
   printf( "Difference vector:\n");
   show_state_vector( state_vect2);
}
