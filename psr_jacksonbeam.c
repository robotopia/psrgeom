/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program plots (i.e. outputs numbers for) the Jackson beam as a
 * function of theta, for a given gamma and frequency.
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
    double  freq;        // The frequency (in Hz)
    double  gamma;       // The Lorentz factor
    double  rho;         // The curvature (in m)
    char   *outfile;     // name of output file (NULL means stdout)
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

int main( int argc, char *argv[] )
{
    // Seed the random number generator
    srand( time( NULL ) );

    // Set up struct for command line options and set default values
    struct opts o;
    o.freq      = NAN;
    o.gamma     = NAN;
    o.rho       = NAN;
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
    double Ipos, Ineg;
    double I, L, V;

    // Loop over theta angles between -1/γ and 1/γ
    int n;
    int N = 1000;
    psr_angle th;
    for (n = 0; n < N; n++)
    {
        set_psr_angle_deg( &th, ((double)n/(double)N - 0.5) * 3.0 );
        particle_beam_intensity( o.freq, o.gamma, &th, o.rho, &Ipos, &Ineg );

        I = Ipos + Ineg;
        V = Ipos - Ineg;
        L = sqrt(4.0*Ipos*Ineg);

        fprintf( f, "%.15e %.15e %.15e %.15e %.15e %.15e\n",
                th.deg, Ipos, Ineg, I, L, V );
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
    printf( "  -f  freq     The emission frequency, in MHz.\n" );
    printf( "  -g  gamma    The Lorentz factor\n" );
    printf( "  -r  rho      The radius of curvature (in km)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "f:g:ho:r:")) != -1)
    {
        switch (c)
        {
            case 'f':
                o->freq = atof(optarg)*1e6; // convert to Hz
                break;
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
            case 'r':
                o->rho = atof(optarg)*1e3; // convert to m
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
    if (isnan(o->freq) || isnan(o->gamma) || isnan(o->rho))
    {
        fprintf( stderr, "error: -f, -g, and -r options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
/* The final output includes:
 *   1&2) the polar coordinates of the footpoint of the magnetic field line
 *   3&4) the polar coordinates of the photon direction (with retardation)
 *   5) the polarisation angle
 *   6) the phase retardation due to the emission height
 *   7) the critical frequency (in MHz)
 *   8) the emission height (in km)
 */
{
    // Print out a line to file handle f
    fprintf( f, "#  th_deg  Ipos  Ineg  I  L  V\n" );
}


