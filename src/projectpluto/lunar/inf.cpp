#include <stdio.h>
#include <stdlib.h>

int main( const int argc, const char **argv)
{
   double z1 = atof( argv[1]);
   double z2 = atof( argv[2]);

   printf( "%lf is %s %lf\n", z1,
         (z1 > z2 ? "greater than" : (z1 == z2 ? "equal to" : "less than")),
         z2);
   return( 0);
}
