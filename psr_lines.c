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
    point   ref_axis;  // the reference axis for s and p
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
    set_point_xyz( &(o.ref_axis), 0.0, 0.0, 1.0, POINT_SET_ALL );

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

    // Check if -m was supplied, and set o.ref_axis accordingly
    if (o.mag_axis)
    {
        psr_angle zero;
        set_psr_angle_deg( &zero, 0.0 );
        set_point_sph( &(o.ref_axis), 1.0, al, &zero, POINT_SET_ALL );
    }

    // Set up points
    point foot_pt;     // The footpoint relative to the ref_axis
    point obs_foot_pt; // The footpoint in the observer frame
    point far_pt;      // The furthest point

    // Print header and column header to file
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Draw the requested lines
    int linetype;  // either CLOSED_LINE or OPEN_LINE
    int s_idx, p_idx;
    double s_deg, p_deg;
    double ds, dp;
    psr_angle s, p;

    for (s_idx = 0; s_idx < o.s_nstep; s_idx++)
    {
        // Convert s_idx to an angle
        ds = (o.s_nstep == 1 ?
                0.0 :
                (o.s_stop - o.s_start)/(o.s_nstep - 1.0));
        s_deg = o.s_start + s_idx*ds;
        set_psr_angle_deg( &s, s_deg );

        for (p_idx = 0; p_idx < o.p_nstep; p_idx++)
        {
            // Convert p_idx to an angle
            dp = (o.p_nstep == 1 ?
                    0.0 :
                    (o.p_stop - o.p_start)/(o.p_nstep - 1.0));
            p_deg = o.p_start + p_idx*dp;
            set_psr_angle_deg( &p, p_deg );

            // From s and p, get a footpoint (rel. to reference axis)
            set_point_sph( &foot_pt, psr.r, &s, &p, POINT_SET_ALL );

            // Rotate it into the observer's frame
            rotate_about_axis( &foot_pt, &obs_foot_pt, &(o.ref_axis.th),
                    'y', POINT_SET_ALL );
            rotate_about_axis( &obs_foot_pt, &obs_foot_pt, &(o.ref_axis.ph),
                    'z', POINT_SET_ALL );

            // Follow the mag field line to the extreme
            linetype = get_fieldline_type( &obs_foot_pt, &psr, o.tmult,
                    o.rL_norm, f, &far_pt );
            if (linetype == CLOSED_LINE)
            {
                footpoint( &far_pt, &psr, o.tmult, DIR_OUTWARD, f, o.rL_norm,
                        1.0, NULL );
            }

            // Insert a blank line into the output
            fprintf( f, "\n\n" );
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
                           "solution (not yet implemented)\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -m           Set the reference axis to the magnetic axis. "
                           "This option overrides the -X option "
                           "(default is reference axis = rotation axis)\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -r  rho_lim  The maximum distance allowed for lines. "
                           "(default = 1.2)\n" );
    printf( "  -t  tmult    The initial size of the RK4 steps, as a fraction "
                           "of the light cylinder radius (default: 0.01)\n" );
    printf( "  -X  col,long The reference axis from which s and p are "
                           "calculated (see -p and -s), in degrees. If -m is "
                           "also given, this option is ignored (default = "
                           "0,0)\n" );

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
            case 'd':
                o->dipole = 1;
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'L':
                o->rL_norm = 1;
                break;
            case 'm':
                o->mag_axis = 1;
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'p':
                parse_range( optarg, &(o->p_start),
                                     &(o->p_stop),
                                     &(o->p_nstep) );
                break;
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 'r':
                o->rho_lim = atof(optarg);
                break;
            case 's':
                parse_range( optarg, &(o->s_start),
                                     &(o->s_stop),
                                     &(o->s_nstep) );
                break;
            case 't':
                o->tmult = atof(optarg);
                break;
            case 'X':
                parse_direction( optarg, &(o->ref_axis) );
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

    if (isnan(o->s_start) || isnan(o->p_start))
    {
        fprintf( stderr, "error: -p and -s options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# x y z tstep\n" );
}


