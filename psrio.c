/****************************************************************************
 * PSRGEOM, psrio.c
 * ----------------
 * Sam McSweeney, 2018
 * sammy.mcsweeney@gmail.com
 *
 * This source file houses functions that relate to outputting quantities
 * relating to pulsar magnetic fields and pulsar geometry to file streams.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "psrgeom.h"

void print_psrg_header( FILE *f, int argc, char *argv[] )
{
    char comment = '#';
    fprintf( f, "%c Made with PSRGEOM v%s\n", comment, PSRGEOM_VERSION );
    fprintf( f, "%c Command line:", comment );
    int n;
    for (n = 0; n < argc; n++)
    {
        if (strstr( argv[n], " " ) == NULL)
            fprintf( f, " %s", argv[n] );
        else
            fprintf( f, " \"%s\"", argv[n] );
    }
    fprintf( f, "\n%c\n", comment );
}


void parse_range( char *str, double *start, double *stop, int *nsteps )
/* Parses a string that is assumed to be in one of the following formats:
 *
 *   start
 *   start:stop
 *   start:stop:nsteps
 *
 * where start and stop are floating point numbers and nsteps is an integer.
 * If the string fails to parse, an error is thrown and execution stops.
 */
{
    int nitems = sscanf( str, "%lf:%lf:%d", start, stop, nsteps );
    if (nitems == 3)
        ; // do nothing, all is well
    else if (nitems == 2)
    {
        *nsteps = 2;
    }
    else if (nitems == 1)
    {
        *stop = *start;
        *nsteps  = 1;
    }
    else
    {
        fprintf( stderr, "error: parse_range: failed to parse \"%s\" as "
                         "\"start:stop:nsteps\"\n", str );
        exit(EXIT_FAILURE);
    }

    if (*nsteps < 1)
    {
        fprintf( stderr, "error: parse_range: nsteps must be > 0\n" );
        exit(EXIT_FAILURE);
    }
}


void parse_direction( char *str, point *direction )
/* Parses a string that is assumed to be in the following format:
 *
 *   colatitude,longitude
 *
 * where both the colatitude and the longitude are in degrees.
 * If the string fails to parse, an error is thrown and execution stops.
 * Otherwise, the colatitude and longitude are converted into a unit
 * vector (cast as a "point") in 3D space.
 */
{
    double col_deg, long_deg;
    int nitems = sscanf( str, "%lf,%lf", &col_deg, &long_deg );
    if (nitems == 2)
    {
        fprintf( stderr, "error: parse_direction: failed to parse \"%s\" as "
                         "\"colatitude,longitude\"\n", str );
        exit(EXIT_FAILURE);
    }

    // Turn the angles into psr_angles
    psr_angle c, l; // (c)olatitude, (l)ongitude
    set_psr_angle_deg( &c, col_deg );
    set_psr_angle_deg( &l, long_deg );

    // Convert into a point struct
    set_point_sph( direction, 1.0, &c, &l, POINT_SET_ALL );
}
