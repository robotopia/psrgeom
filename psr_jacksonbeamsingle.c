/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program plots (i.e. outputs numbers for) the single-particle beam as a
 * function of theta and phi, for a given gamma.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "psrgeom.h"

struct opts
{
    double  gamma;       // The Lorentz factor
    char   *outfile;     // name of output file (NULL means stdout)
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

int main( int argc, char *argv[] )
{
    // Set up struct for command line options and set default values
    struct opts o;
    o.gamma     = NAN;
    o.outfile   = NULL;

    parse_cmd_line( argc, argv, &o );

    // Set up output file
    FILE *f;
    if (o.outfile == NULL)
        f = stdout;
    else
    {
        f = fopen( o.outfile, "w" );
        if (f == NULL)
        {
            fprintf( stderr, "error: could not open file %s\n", o.outfile );
            exit(EXIT_FAILURE);
        }
    }

    // Write the file header
    print_psrg_header( f, argc, argv );

    // Set up variables for the answers
    double power;

    // Loop over theta angles between -4/γ and 4/γ
    int xi, yi;
    double x, y;
    int X = 500, Y = 500;
    double vdot = 1.0;
    double pixel_x_size = 1000.0/o.gamma / X;
    double pixel_y_size = 1000.0/o.gamma / Y;
    psr_angle th, ph;
    for (xi = 0; xi < X; xi++)
    {
        x = (double)(xi - X/2)/(double)X * pixel_x_size;

        for (yi = 0; yi < Y; yi++)
        {
            y = (double)(yi - Y/2)/(double)Y * pixel_y_size;
            set_psr_angle_rad( &th, sqrt(x*x + y*y) );
            set_psr_angle_rad( &ph, atan2(y, x) );
            power = single_particle_power_perp( o.gamma, &th, &ph, vdot );

            fprintf( f, "%.15e %.15e %.15e\n", x*RAD2DEG, y*RAD2DEG, power );
        }
        fprintf( f, "\n" );
    }

    free( o.outfile );

    if (o.outfile != NULL)
        fclose( f );

    return 0;
}



void usage()
{
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -g  gamma    The Lorentz factor\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "g:ho:r:")) != -1)
    {
        switch (c)
        {
            case 'g':
                o->gamma = atof(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case '?':
                fprintf( stderr, "error: unknown option character '-%c'\n",
                         optopt );
                exit(EXIT_FAILURE);
                break;
            default:
                fprintf( stderr, "error: couldn't parse command line\n" );
                exit(EXIT_FAILURE);
        }
    }

    // Check that all the arguments are valid
    if (isnan(o->gamma))
    {
        fprintf( stderr, "error: -g option required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "#  x_deg  y_deg  power\n" );
}


