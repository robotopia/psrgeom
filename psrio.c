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
#include "psrgeom.h"

void print_header( FILE *f, int argc, char *argv[] )
{
    char comment = '#';
    fprintf( f, "%c Made with PSRGEOM %s\n", comment, PSRGEOM_VERSION );
    fprintf( f, "%c Command line:", comment );
    int n;
    for (n = 0; n < argc; n++)
        fprintf( f, " %s", argv[n] );
    fprintf( f, "\n%c\n", comment );
}
