/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program prints out a series of points (in the magnetic frame)
 * representing the line of sight as the pulsar rotates.
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
    double  al_deg;      // alpha angle in deg
    double  ze_deg;      // zeta angle in deg
    char   *outfile;     // name of output file (NULL means stdout)
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

int main( int argc, char *argv[] )
{
    // Seed the random number generator
    srand( time( NULL ) );

    // Generic counter:
    int i;

    // Set up struct for command line options and set default values
    struct opts o;
    o.al_deg    = NAN;
    o.ze_deg    = NAN;
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

    // Set up pulsar
    pulsar psr;

    psr_angle al;
    psr_angle ze;
    set_psr_angle_deg( &al, o.al_deg );
    set_psr_angle_deg( &ze, o.ze_deg );

    // The following values don't affect the output of this program
    double P = 1.0;
    double r = 1e4;

    set_pulsar( &psr, NULL, NULL, P, r, &al, &ze );

    // Write the file header
    print_psrg_header( f, argc, argv );

    // Write the column headers
    print_col_headers( f );

    point LoS, LoS_mag;
    psr_angle LoS_ph;

    for (i = 0; i < 360; i++)
    {
        // Set up angle
        set_psr_angle_deg( &LoS_ph, (double)i );

        // Create point in the observer frame, on a unit sphere
        set_point_sph( &LoS, 1.0, &(psr.ze), &LoS_ph, POINT_SET_ALL );

        // Convert it into the magnetic frame
        obs_to_mag_frame( &LoS, &psr, NULL, &LoS_mag );

        // Print out the result
        fprintf( f, "%.15e %.15e %.15e\n",
                LoS_ph.deg,
                LoS_mag.th.deg,
                LoS_mag.ph.deg );
    }

    // Clean up

    free( o.outfile );

    if (o.outfile != NULL)
        fclose( f );

    return 0;
}

void usage()
{
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation axis and the "
                           "magnetic axis in degrees (required)\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:ho:z:")) != -1)
    {
        switch (c)
        {
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'z':
                o->ze_deg = atof(optarg);
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
    if (isnan(o->al_deg) || isnan(o->ze_deg))
    {
        fprintf( stderr, "error: -a and -z options required"
                         "\n" );
        usage();
        exit(EXIT_FAILURE);
    }

}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# phase_deg  th_deg  ph_deg\n" );
}


