/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program attempts to simulate the emerging beam pattern as
 * observed from a pulsar with given angles α and ζ, period P, within a
 * frequency range [f₁:f₂].
 *
 * Here, the idea is that emission is allowed to originate on any open field
 * line, and at any height. The photon that would be emitted at said field
 * line and height is collected in an array where each pixel (in Cartesian
 * coordinates) represents a particular direction on the (pulsar's) sky.
 *
 * Because of the dynamic nature of the spark carousel, the output is not
 * just a single array, but several. Each array represents the instantaneous
 * beam sampled at a particular time.
 *
 * The final output includes:
 *   1) the x-coordinate of the beam pixel
 *   2) the y-coordinate of the beam pixel
 *   3,4,...) the beam at different times
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
    double  P_sec;       // period, in sec
    char   *outfile;     // name of output file (NULL means stdout)
    double  s_deg;       // the spark angular size
    double  S_deg;       // the carousel angular radius
    double  f_lo;        // minimum value of freq (MHz)
    double  f_hi;        // maximum value of freq (MHz)
    int     open_only;   // only consider open field lines
    int     num_lines;   // sample this many lines
    int     nsparks;     // number of sparks in carousel
    int     dipole;      // use dipole field?
    double  P4_sec;      // the rotation time of the carousel
    int     csl_type;    // the type of spark profile (TOPHAT or GAUSSIAN)
    int     nframes;     // Create this many time frames
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
    o.outfile   = NULL;
    o.s_deg     = NAN;
    o.S_deg     = NAN;
    o.f_lo      = NAN;
    o.f_hi      = NAN;
    o.open_only = 0;
    o.num_lines = 10000;
    o.nsparks   = 0;
    o.dipole    = 0; // use Deutsch field by default
    o.P4_sec    = NAN;
    o.csl_type  = GAUSSIAN;
    o.nframes   = 10;

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

    // Set up pulsar's carousel
    psr_angle s, S;
    set_psr_angle_deg( &s, o.s_deg );
    set_psr_angle_deg( &S, o.S_deg );
    set_pulsar_carousel( &psr, o.nsparks, &s, &S, o.csl_type, o.P4_sec );

    // Calculate the time associated with each frame. Instead of making the
    // frames span one whole carousel rotation, make them span the time
    // for one spark to cross the distance to the next spark.
    double tstep = psr.csl.P4 / (double)psr.csl.n / (double)o.nframes;

    // Write the file header
    print_psrg_header( f, argc, argv );

    // Some needed variables
    int linetype;   // either CLOSED_LINE or OPEN_LINE
    point foot_pt;  // a randomly chosen foot_point
    double height;  // a randomly chosen height

    // Write the column headers
    print_col_headers( f, o.nframes, tstep );

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
            int rL_norm = 0;
            linetype = get_fieldline_type( &foot_pt, &psr, rL_norm,
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

        climb_and_emit( &psr, &init_pt, o.gamma, o.f_start*1.0e6,
                o.f_stop*1.0e6, f ); /* (1.0e6: convert frequencies to Hz) */
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

    o.scl_type  = GAUSSIAN;
void usage()
{
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -f  f1:f2    The observing frequency, in MHz. "
                           "The range is from f1 to f2.\n" );
    printf( "  -P  period   The rotation period of the pulsar (in sec)\n" );
    printf( "  -s  s_deg    The angular size of the sparks (in deg)\n" );
    printf( "  -S  S_deg    The carousel's angular radius (in deg)\n" );
    printf( "  -4  P4       The carousel's rotation period (in sec)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -c type      The spark profile type, either GAUSSIAN (default) "
                           "or TOPHAT\n" );
    printf( "  -d           Use a dipole field instead of the default "
                           "Deutsch field\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -n  photons  Sample this many photons (default: 10000)\n" );
    printf( "  -N  nsparks  The number of sparks in the carousel. If nsparks "
                           "= 0 (default), the carousel is a \"solid\" "
                           "annulus.\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -O           Only consider open field lines (default: off)\n" );
    printf( "  -t  nframes  Create this many time frames (default = 10)\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:c:df:hn:N:o:OP:s:S:t:4:")) != -1)
    {
        switch (c)
        {
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'c':
                if (strcmp( optarg, "GAUSSIAN" ) == 0)
                    o->csl_type = GAUSSIAN;
                else if (strcmp( optarg, "TOPHAT" ) == 0)
                    o->csl_type = TOPHAT;
                else
                {
                    fprintf( stderr, "error: -c argument must be either "
                            "GAUSSIAN or TOPHAT\n" );
                    exit(EXIT_FAILURE);
                }
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
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 's':
                o->s_deg = atof(optarg);
                break;
            case 'S':
                o->S_deg = atof(optarg);
                break;
            case 't':
                o->nframes = atoi(optarg);
                break;
            case '4':
                o->P4_sec = atof(optarg);
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
    if (isnan(o->al_deg) || isnan(o->P_sec) || isnan(o->P4_sec))
    {
        fprintf( stderr, "error: -a, -P and -4 options required"
                         "\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (isnan(o->s_deg) || isnan(o->S_deg) || isnan(o->f_start))
    {
        fprintf( stderr, "error: -f, -s and -S options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (o->nsparks < 0)
    {
        fprintf( stderr, "error: -N (=%d) must be >= 0\n", o->nsparks );
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f, int nframes, double tstep )
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
    fprintf( f, "#  x_deg  y_deg" );
    int i;
    for (i = 0; i < nframes; i++)
        fprintf( f, "  t=%fs", (double)i*tstep );
    fprintf( f, "\n" );
}


