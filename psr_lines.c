#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "psrgeom.h"
#include <time.h>

struct opts
{
    double  al_deg;    // alpha angle in deg
    double  ze_deg;    // zeta angle in deg
    double  P_sec;     // period, in sec
    double  tmult;     // RK4 step size, as a fraction of lt cyl radius
    int     rL_norm;   // bool: normalise output to light cylinder radius?
    double  rho_lim;   // limit to within rho_lim * light cylinder radius?
    char   *outfile;   // name of output file (NULL means stdout)
    double  s_start;   // starting value of s
    double  s_stop;    // stopping value of s
    int     s_nstep;   // number of s steps
    double  p_start;   // starting value of p
    double  p_stop;    // stopping value of p
    int     p_nstep;   // number of p steps
    int     dipole;    // bool: use dipole
    int     mag_axis;  // bool: use magnetic axis as reference axis
    double  axis_col;  // colatitude of reference axis (deg)
    double  axis_long; // longitude of reference axis (deg)
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

int main( int argc, char *argv[] )
{
    // Seed random number generator
    srand( time( NULL ) );

    // Set up struct for command line options and set default values
    struct opts o;
    o.al_deg    = NAN;
    o.P_sec     = NAN;
    o.ze_deg    = NAN;
    o.tmult     = 0.01;
    o.rL_norm   = 0;
    o.rho_lim   = 1.2;
    o.outfile   = NULL;
    o.s_start   = NAN;
    o.s_stop    = NAN;
    o.s_nstep   = 0;
    o.p_start   = NAN;
    o.p_stop    = NAN;
    o.p_nstep   = 0;
    o.dipole    = 0;
    o.mag_axis  = 0;
    o.axis_col  = 0.0;
    o.axis_long = 0.0;

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

    double P = o.P_sec;  // Pulsar spin period
    double r = 1e4;      // Pulsar radius

    set_pulsar( &psr, ra, dec, P, r, al, ze );

    // Set up points
    point emit_pt;
    point foot_pt; // The end point of either footpoint() or farpoint()

    // Print header and column header to file
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Start at φ = 0° at work around
    psr_angle ph;
    set_psr_angle_deg( &ph, 0.0 );

    // Initially, set the init_pt to an approximate initial guess
    find_approx_emission_point( &psr, &ph, DIR_OUTWARD, &emit_pt );

    // Make sure the initial "previous" point is at least 2 pulsar radii above
    // the pulsar's surface
    double min_height = 3.0*psr.r;
    if (emit_pt.r < min_height)
    {
        psr_angle za; // "zero angle"
        set_psr_angle_rad( &za, 0.0 );
        set_point_sph( &emit_pt, min_height,
                                 &psr.al,
                                 &za,
                                 POINT_SET_ALL );
    }

    int N = 360;
    int i;
    for (i = 0; i < N; i++)
    {
        // Set the rotation phase
        set_psr_angle_deg( &ph, i*360.0/N );

        // Get the emission pt (guaranteed to be on a last open field line)
        int status;
        status = find_emission_point_elevator( &psr, &ph, DIR_OUTWARD,
                                               &emit_pt, &emit_pt, NULL );
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

        // Trace the line back and forwards to the footpoints
        footpoint( &emit_pt, &psr, o.tmult, DIR_INWARD, f, o.rL_norm,
                   o.rho_lim, &foot_pt );
        fprintf( f, "\n" );
        footpoint( &emit_pt, &psr, o.tmult, DIR_OUTWARD, f, o.rL_norm,
                   o.rho_lim, &foot_pt );
        fprintf( f, "\n\n" );
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
    printf( "usage: psr_lines [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -p  p[:P[:n]]   The azimuth relative to the reference axis, "
                           "in degrees. The range is from p to P with n "
                           "steps.\n"
            "                 p      ==> p:p:1\n"
            "                 p:P    ==> p:P:2\n" );
    printf( "  -s  s[:S[:n]]   The angular distance from the reference axis, "
                           "in degrees. The range is from s to S with n "
                           "steps.\n"
            "                 s      ==> s:s:1\n"
            "                 s:S    ==> s:S:2\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -d           Use the dipole model instead of the full Deutsch "
                           "solution\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -m           Set the reference axis to the magnetic axis "
                           "(default is reference axis = rotation axis)\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -r  rho_lim  The maximum distance allowed for lines. "
                           "(default = 1.2)\n" );
    printf( "  -t  tmult    The initial size of the RK4 steps, as a fraction "
                           "of the light cylinder radius (default: 0.01)\n" );
    printf( "  -X  col:long The reference axis from which s and p are "
                           "calculated (see -p and -s), in degrees (default ="
                           " 0:0)\n" );

}



void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:dhLmo:p:P:r:s:t:X:z:")) != -1)
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
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 'r':
                o->rho_lim = atof(optarg);
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
    fprintf( f, "# x y z\n" );
}


