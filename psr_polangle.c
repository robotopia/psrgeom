#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "psrgeom.h"

struct opts
{
    double  al_deg;    // alpha angle in deg
    double  ze_deg;    // zeta angle in deg
    double  P_sec;     // period, in sec
    int     rL_norm;   // bool: normalise to light cylinder radius?
    int     direction; // either DIR_OUTWARD or DIR_INWARD
    char   *outfile;   // name of output file (NULL means stdout)
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

int main( int argc, char *argv[] )
{
    // Set up struct for command line options and set default values
    struct opts o;
    o.al_deg    = NAN;
    o.P_sec     = NAN;
    o.ze_deg    = NAN;
    o.rL_norm   = 0;
    o.direction = DIR_OUTWARD;
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

    psr_angle *ra  = NULL;
    psr_angle *dec = NULL;

    psr_angle *al = create_psr_angle_deg( o.al_deg );
    psr_angle *ze = create_psr_angle_deg( o.ze_deg );

    double P = o.P_sec;
    double r = 1e4; /* This doesn't actually make a difference to any of the
                       outputs of this program */

    set_pulsar( &psr, ra, dec, P, r, al, ze );

    // Write the file and column headers
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Calculate the polarisation angle for each phase
    point emit_pt, prev_pt;
    psr_angle ph, psi;

    // Initially, set the prev_pt to an approximate initial guess
    set_psr_angle_deg( &ph, 0.0 );
    find_approx_emission_point( &psr, &ph, o.direction, &prev_pt );

    // Make sure the initial "previous" point is at least 2 pulsar radii above
    // the pulsar's surface
    double min_height = 3.0*psr.r;
    if (prev_pt.r < min_height)
    {
        psr_angle za; // "zero angle"
        set_psr_angle_rad( &za, 0.0 );
        set_point_sph( &prev_pt, min_height,
                                 &psr.al,
                                 &za,
                                 POINT_SET_ALL );

        if (o.direction == DIR_INWARD)
        {
            set_point_xyz( &prev_pt, -prev_pt.x[0],
                                     -prev_pt.x[1],
                                     -prev_pt.x[2],
                                     POINT_SET_ALL );
        }
    }

    // Loop over the phases and calculate the new pos angle, using the
    // previous emission point as a seed for each new emission point
    double ph_deg;

    double ph_start = (o.direction == DIR_OUTWARD ? 0.0 : -180.0);

    int N = 360; // The number of points to sample
    int i;
    for (i = 0; i < N+1; i++)
    {
        ph_deg = i*360.0/(double)N + ph_start;
        set_psr_angle_deg( &ph, ph_deg );
        if (calc_pol_angle( &psr, &ph, o.direction, &prev_pt, &emit_pt, &psi ))
        {
            fprintf( f, "%.15e %.15e %.15e %.15e %.15e\n",
                        emit_pt.x[0], emit_pt.x[1], emit_pt.x[2],
                        (ph.deg > 180.0 ? ph.deg - 360.0 : ph.deg), psi.deg );
            if (i == 180 && o.direction == DIR_OUTWARD)
            {
                fprintf( f, "\n%.15e %.15e %.15e %.15e %.15e\n",
                        emit_pt.x[0], emit_pt.x[1], emit_pt.x[2],
                        ph.deg - 360.0, psi.deg );
            }
            copy_point( &emit_pt, &prev_pt );
        }
        else
        {
            fprintf( f, "\n" );
        }
    }

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
    printf( "usage: psr_polangle [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -i           Calculate the emission point for particles that "
                           "flow in the opposite direction to the magnetic "
                           "field\n"
            "               (default is to assume particles are flowing in "
                           "the same direction as the magnetic field)\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:hiLo:P:r:z:")) != -1)
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
            case 'i':
                o->direction = DIR_INWARD;
                break;
            case 'L':
                o->rL_norm = 1;
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'P':
                o->P_sec = atof(optarg);
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
    if (isnan(o->al_deg) || isnan(o->P_sec) || isnan(o->ze_deg))
    {
        fprintf( stderr, "error: -a, -P, and -z options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# x y z rotation_phase(φ)  pol_angle(Ψ)\n" );
}


