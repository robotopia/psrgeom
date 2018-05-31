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
    double  ph_deg;    // phase angle in deg
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
    o.ph_deg    = NAN;
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

    // Set up rotation phase
    psr_angle *ph = create_psr_angle_deg( o.ph_deg );

    // Set up point for answer, and for initial seed guess
    point emit_pt, init_pt;

    // Initially, set the init_pt to an approximate initial guess
    find_approx_emission_point( &psr, ph, o.direction, &init_pt );

    // Make sure the initial "previous" point is at least 2 pulsar radii above
    // the pulsar's surface
    double min_height = 3.0*psr.r;
    if (init_pt.r < min_height)
    {
        psr_angle za; // "zero angle"
        set_psr_angle_rad( &za, 0.0 );
        set_point_sph( &init_pt, min_height,
                                 &psr.al,
                                 &za,
                                 POINT_SET_ALL );
        if (o.direction == DIR_INWARD)
        {
            set_point_xyz( &init_pt, -init_pt.x[0],
                                     -init_pt.x[1],
                                     -init_pt.x[2],
                                     POINT_SET_ALL );
        }
    }

    // Write the file and column headers
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Calculate answer
    int status;
    status = find_emission_point_elevator( &psr, ph, o.direction,
                                           &init_pt, &emit_pt, f );

    if (status == EMIT_PT_TOO_HIGH)
    {
        fprintf( stderr, "# Stopped at light cylinder\n" );
        exit(EXIT_SUCCESS);
    }
    else if (status == EMIT_PT_TOO_LOW)
    {
        fprintf( stderr, "# Stopped at pulsar surface\n" );
        exit(EXIT_SUCCESS);
    }

    // Print out the features of the found point.
    point V1, V2;
    point LoS;

    calc_fields( &emit_pt, &psr, SPEED_OF_LIGHT, NULL, &V1, &V2, NULL, NULL, NULL );
    line_of_sight( &psr, ph, &LoS );

    fprintf( f, "# V   = ( %.15e, %.15e, %.15e )\n", V1.x[0], V1.x[1], V1.x[2] );
    fprintf( f, "# LoS = ( %.15e, %.15e, %.15e )\n", LoS.x[0], LoS.x[1], LoS.x[2] );

    point far_pt;
    farpoint( &emit_pt, &psr, NULL, 0, 1.1, &far_pt );

    fprintf( f, "# extreme point has ρ/r_L = %.15f\n", sqrt( far_pt.rhosq / psr.rL2 ) );

    // Clean up
    destroy_psr_angle( ra  );
    destroy_psr_angle( dec );
    destroy_psr_angle( al  );
    destroy_psr_angle( ze  );
    destroy_psr_angle( ph  );

    free( o.outfile );

    if (o.outfile != NULL)
        fclose( f );

    return 0;
}

void usage()
{
    printf( "usage: psr_emit [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -p  phase    The rotation phase of the observed emission, in "
                           "degrees (required)\n" );
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
    while ((c = getopt( argc, argv, "a:hiLo:p:P:r:z:")) != -1)
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
            case 'p':
                o->ph_deg = atof(optarg);
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
    if (isnan(o->al_deg) || isnan(o->P_sec) ||
        isnan(o->ze_deg) || isnan(o->ph_deg))
    {
        fprintf( stderr, "error: -a, -p, -P, and -z options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# x  y  z  psr_cost_lofl  psr_cost_los\n" );
}


