#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int sat_id_main( const int argc, const char **argv);
void avoid_runaway_process( const int max_time_to_run);   /* cgi_func.c */
int get_multipart_form_data( const char *boundary, char *field,
                char *buff, char *filename, const size_t max_len);

int main( const int unused_argc, const char **unused_argv)
{
   const char *argv[20];
   const size_t max_buff_size = 40000;       /* room for 500 obs */
   char *buff = (char *)malloc( max_buff_size);
   char boundary[100], field[30];
   const char *temp_obs_filename = "sat_obs.txt";
   double search_radius = 2.;     /* look 2 degrees for matches */
   int argc = 3;
   FILE *lock_file = fopen( "lock.txt", "w");
   size_t bytes_written = 0;
   extern char **environ;

#ifndef _WIN32
   avoid_runaway_process( 15);
#endif         /* _WIN32 */
   printf( "Content-type: text/html\n\n");
   printf( "<html> <body> <pre>\n");
   if( !lock_file)
      {
      printf( "<p> Server is busy.  Try again in a minute or two. </p>");
      printf( "<p> Your orbit is very important to us! </p>");
      return( 0);
      }
   fprintf( lock_file, "We're in\n");
   for( size_t i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
   if( !fgets( boundary, sizeof( boundary), stdin))
      {
      printf( "<b> No info read from stdin</b>");
      printf( "This isn't supposed to happen.\n");
      return( 0);
      }
   while( get_multipart_form_data( boundary, field, buff, NULL, max_buff_size) >= 0)
      {
      if( !strcmp( field, "TextArea") || !strcmp( field, "upfile"))
         {
         if( strlen( buff) > 70)
            {
            FILE *ofile = fopen( temp_obs_filename,
                               (bytes_written ? "ab" : "wb"));

            bytes_written += fwrite( buff, 1, strlen( buff), ofile);
            fclose( ofile);
            }
         }
      if( !strcmp( field, "radius"))
         search_radius = atof( buff);
      }
   argv[0] = "sat_id";
   argv[1] = temp_obs_filename;
   argv[2] = "-ttle_list.txt";
   sprintf( field, "-r%.2f", search_radius);
   argv[3] = field;
   argv[4] = NULL;
   sat_id_main( argc, argv);
   free( buff);
   printf( "</pre> </body> </html>");
   return( 0);
}
