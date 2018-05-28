/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program attempts to simulate the profile as observed from a pulsar
 * with given angles α and ζ, period P, and frequency range f1 to f2.
 *
 * For each sampled field line, the line is followed outward in incremental
 * steps. At each step, the average beam pattern of emitting particles is
 * calculated, with the gamma factors of the particles drawn from a pre-
 * determined distribution.
 *
 * The final output is a two-column file:
 *   1) the rotation phase in degrees
 *   2) the profile power (in arbitrary units)
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "psrgeom.h"

struct opts
{
    double  al_deg;      // alpha angle in deg
    double  ze_deg;      // zeta angle in deg
    double  P_sec;       // period, in sec
    char   *outfile;     // name of output file (NULL means stdout)
    double  s_start;     // starting value of s
    double  s_stop;      // stopping value of s
    double  p_start;     // starting value of p
    double  p_stop;      // stopping value of p
    double  f_start;     // starting value of freq (MHz)
    double  f_stop;      // stopping value of freq (MHz)
    int     open_only;   // only consider open field lines
    int     num_lines;   // sample this many lines
    int     nsparks;     // number of sparks in carousel
    int     dipole;      // use dipole field?
    int     nbins;       // number of profile phase bins
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
    o.P_sec     = NAN;
    o.ze_deg    = NAN;
    o.outfile   = NULL;
    o.s_start   = NAN;
    o.s_stop    = NAN;
    o.p_start   = 0.0;
    o.p_stop    = 360.0;
    o.f_start   = NAN;
    o.f_stop    = NAN;
    o.open_only = 0;
    o.num_lines = 10000;
    o.nsparks   = 0;
    o.dipole    = 0; // use Deutsch field by default
    o.nbins     = 1024;

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

    if (o.dipole)
        psr.field_type = DIPOLE;

    // Write the file header
    print_psrg_header( f, argc, argv );

    // Some needed variables
    int linetype;   // either CLOSED_LINE or OPEN_LINE
    point foot_pt, foot_pt_mag;
    point init_pt;
    double profile[o.nbins];
    int bin_count[o.nbins];
    int centre_bin = o.nbins/2;

    // Some default values
    int    rL_norm = 0;
    double tmult   = 0.01;

    // Reset profile to zero
    for (i = 0; i < o.nbins; i++)
    {
        profile[i]   = 0.0;
        bin_count[i] = 0;
    }

    // Write the column headers
    print_col_headers( f );

//#pragma omp parallel for
    for (i = 0; i < o.num_lines; i++)
    {
        // Obtain a random point on the pulsar surface
        if (o.nsparks == 0)
        {
            random_direction_bounded( &foot_pt_mag, o.s_start*DEG2RAD,
                    o.s_stop*DEG2RAD, o.p_start*DEG2RAD, o.p_stop*DEG2RAD );
        }
        else /* a carousel of sparks! */
        {
            random_direction_spark( &foot_pt_mag,
                    (o.s_start + o.s_stop )/2.0 * DEG2RAD,
                    (o.s_stop  - o.s_start)/2.0 * DEG2RAD,
                    o.nsparks );
        }
        scale_point( &foot_pt_mag, psr.r, &foot_pt_mag );

        // Convert the foot_pt into observer coordinates
        mag_to_obs_frame( &foot_pt_mag, &psr, NULL, &foot_pt );

        // If requested, check that we're on an open field line
        if (o.open_only)
        {
            linetype = get_fieldline_type( &foot_pt, &psr, tmult, rL_norm,
                    NULL, NULL );
            if (linetype == CLOSED_LINE)
            {
                continue;
            }
        }

        // Now climb up the field line, emitting as we go
        // Start 1 metre above the surface
        Bstep( &foot_pt, &psr, 1.0, DIR_OUTWARD, &init_pt );
        set_point_xyz( &init_pt, init_pt.x[0],
                init_pt.x[1],
                init_pt.x[2],
                POINT_SET_ALL );

        fieldline_to_profile( &psr, &init_pt, o.f_start*1.0e6, o.f_stop*1.0e6,
                o.nbins, centre_bin, profile, bin_count );
    }

    // Print out the profile
    double phase_deg;
    double bin_width = 360.0 / (double)o.nbins;
    for (i = 0; i < o.nbins; i++)
    {
        // Normalise the profile
        if (bin_count[i] > 0)
            profile[i] /= (double)bin_count[i];

        // Convert bin number to phase
        phase_deg = (double)(i - centre_bin) * bin_width;

        fprintf( f, "%.15e %.15e\n", phase_deg, profile[i] );
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
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -f  f1:f2    The emission frequency, in MHz. "
                           "The range is from f1 to f2.\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -s  s1:s2    The angular distance from the magnetic axis, "
                           "in degrees. The range is from s1 to s2.\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -b  nbins    The number of bins in the output profile\n" );
    printf( "  -d           Use a dipole field instead of the default "
                           "Deutsch field\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -n  nlines   Sample nlines magnetic field lines "
                           "(default: 10000)\n" );
    printf( "  -N  nsparks  The number of sparks in the carousel. If nsparks "
                           "= 0 (default), the footpoints are sampled "
                           "uniformly in the range given by -s. Otherwise, "
                           "the s-range is used to define the spark size.\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -O           Only consider open field lines (default: off)\n" );
    printf( "  -p  p1:p2    The azimuth relative to the magnetic axis, "
                           "in degrees. The range is from p1 to p2. Ensure "
                           "p1 < p2 [default = 0:360]\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:b:df:hn:N:o:Op:P:s:S:z:")) != -1)
    {
        switch (c)
        {
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'b':
                o->nbins = atoi(optarg);
                break;
            case 'd':
                o->dipole = 1;
                break;
            case 'f':
                parse_range( optarg, &(o->f_start),
                                     &(o->f_stop),
                                     NULL );
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'n':
                o->num_lines = atoi(optarg);
                break;
            case 'N':
                o->nsparks = atoi(optarg);
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'O':
                o->open_only = 1;
                break;
            case 'p':
                parse_range( optarg, &(o->p_start),
                                     &(o->p_stop),
                                     NULL );
                break;
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 's':
                parse_range( optarg, &(o->s_start),
                                     &(o->s_stop),
                                     NULL );
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
        fprintf( stderr, "error: -a, -P and -z options required"
                         "\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (isnan(o->s_start) || isnan(o->f_start))
    {
        fprintf( stderr, "error: -f and -s options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (o->nsparks < 0)
    {
        fprintf( stderr, "error: -N (=%d) must be >= 0\n", o->nsparks );
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
/* The final output includes:
 *   1) the rotation phase in degrees
 *   2) the profile power (in arbitrary units)
 */
{
    // Print out a line to file handle f
    fprintf( f, "#  phase_deg  power\n" );
}


