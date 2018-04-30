/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program attempts to simulate the spread of polarisation angles
 * observed from a pulsar with given angles α and ζ, and period P.
 * Here, the idea is that emission is allowed to originate on any open field
 * line. No attempt is made to sample the active field lines realistically.
 * To keep things simple, we start at the magnetic pole and work our way
 * outwards in concentric circles on the polar cap until we reach a "radius"
 * where no open field lines are found.
 *
 * For each open field line, the line is followed outward in incremental steps
 * until the velocity field is found to make an angle ζ with the rotation
 * axis, which means that there is some rotation phase at which curvature
 * radiation from that point becomes visible. It is possible that there are
 * multiple emission heights that fit the ζ criterion, so even after a point
 * has been successfully found, the algorithm continues to search for other
 * points further along the field line.
 *
 * For each point that satisfies the ζ criterion, the acceleration vector is
 * evaluated and the polarisation angle found. Before printing out the
 * polarisation, the rotation phase is adjusted for retardation effects.
 *
 * The final output includes:
 *   1) the polar coordinates of the footpoint of the magnetic field line
 *   2) the Cartesian coordinates of the emission point
 *   3) the retarded (observed) emission phase
 *   4) the observed polarisation angle
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "psrgeom.h"

struct opts
{
    double  al_deg;      // alpha angle in deg
    double  ze_deg;      // zeta angle in deg
    double  P_sec;       // period, in sec
    int     rL_norm;     // bool: normalise to light cylinder radius?
    char   *outfile;     // name of output file (NULL means stdout)
    double  foot_ds_deg; // step size in "distance" from mag. pole in deg
    double  foot_dp_deg; // step size in "φ" around polar cap in deg
    double  tmult;       // step size along field lines
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

int main( int argc, char *argv[] )
{
    // Set up struct for command line options and set default values
    struct opts o;
    o.al_deg      = NAN;
    o.P_sec       = NAN;
    o.ze_deg      = NAN;
    o.rL_norm     = 0;
    o.outfile     = NULL;
    o.foot_ds_deg = NAN;
    o.foot_dp_deg = NAN;
    o.tmult       = NAN;

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

    psr_angle *ra  = NULL;
    psr_angle *dec = NULL;

    psr_angle *al = create_psr_angle_deg( o.al_deg );
    psr_angle *ze = create_psr_angle_deg( o.ze_deg );

    double P = o.P_sec;
    double r = 1e4; /* This will be used later as we move outwards from the
                       pulsar surface */

    set_pulsar( &psr, ra, dec, P, r, al, ze );

    // Set scaling factor in case user requests result in units of rL
    double xscale = (o.rL_norm ? 1.0/psr.rL : 1.0);

    // Write the file and column headers
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

/**** BOOKMARK ****/


/**** End calculation ****/

    // Clean up
    destroy_psr_angle( ra  );
    destroy_psr_angle( dec );
    destroy_psr_angle( al  );
    destroy_psr_angle( ze  );

    free( o.outfile );

    if (o.outfile != NULL)
        fclose( f );

    return 0;
}

void usage()
{
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -p  step     Step size for distance from magnetic pole, in "
                           "degrees (required)\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -s  step     Step size for moving around the polar cap, in "
                           "degrees (required)\n" );
    printf( "  -t  step     Step size for moving along magnetic field lines, "
                           "as a fraction of the light cylinder radius "
                           "(default: 0.01)\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:hLo:p:P:s:t:z:")) != -1)
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
            case 'L':
                o->rL_norm = 1;
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'p':
                o->foot_dp_deg = atof(optarg);
                break;
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 's':
                o->foot_ds_deg = atof(optarg);
                break;
            case 't':
                o->tmult = atof(optarg);
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
    if (isnan(o->al_deg) || isnan(o->P_sec) || isnan(o->ze_deg) ||
        isnan(o->foot_dp_deg) || isnan(o->foot_ds_deg) || isnan(o->tmult))
    {
        fprintf( stderr, "error: -a, -p, -P, -s, -t, and -z options required"
                         "\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
/* The final output includes:
 *   1) the polar coordinates of the footpoint of the magnetic field line
 *   2) the Cartesian coordinates of the emission point
 *   3) the retarded (observed) emission phase
 *   4) the observed polarisation angle
 */
{
    // Print out a line to file handle f
    fprintf( f, "#  s_deg  p_deg  x  y  z  phase  pol\n" );
}


