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
}

