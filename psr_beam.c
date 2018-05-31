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
    int     nphotons;    // sample this many photons
    int     nsparks;     // number of sparks in carousel
    int     dipole;      // use dipole field?
    double  P4_sec;      // the rotation time of the carousel
    int     csl_type;    // the type of spark profile (TOPHAT or GAUSSIAN)
    int     nframes;     // Create this many time frames
    double  g_idx;       // The power law index the gamma distribution
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f, int nframes, double tstep );

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
    o.nphotons  = 10000;
    o.nsparks   = 0;
    o.dipole    = 0; // use Deutsch field by default
    o.P4_sec    = NAN;
    o.csl_type  = GAUSSIAN;
    o.nframes   = 10;
    o.g_idx     = -6.2;

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
    set_psr_angle_deg( &ze, 0.0 ); // This isn't used here

    double P = o.P_sec;
    double r = 1e4; /* This will be used later as we move outwards from the
                       pulsar surface */

    set_pulsar( &psr, NULL, NULL, P, r, &al, &ze );

    if (o.dipole)
        psr.field_type = DIPOLE;

    // Set up pulsar's carousel
    psr_angle s, S;
    set_psr_angle_deg( &s, o.s_deg );
    set_psr_angle_deg( &S, o.S_deg );
    set_pulsar_carousel( &psr, o.nsparks, &s, &S, o.csl_type, o.P4_sec );

    // Calculate the time associated with each frame. For convenience,
    // the frames span one pulsar rotation period.
    double tstep = psr.P / (double)o.nframes;
    double t;
    int tbin;

    // Write the file header
    print_psrg_header( f, argc, argv );

    // Some needed variables
    int linetype;                // either CLOSED_LINE or OPEN_LINE
    point foot_pt, foot_pt_mag;  // a randomly chosen foot_point
    double height;               // a randomly chosen height
    double max_height;
    double freq;                 // a randomly chosen frequency (Hz)
    photon pn;
    double weight;
    int pulse_number;
    int rL_norm = 0;

    // Set up the output array
    double xmax = 90.0; // deg
    double ymax = 90.0; // deg
    double xmin = -xmax; // deg
    double ymin = -ymax; // deg
    int X = 200; // size of array in x-direction
    int Y = 200; // size of array in y-direction
    double xbinwidth = (xmax - xmin)/(double)X;
    double ybinwidth = (ymax - ymin)/(double)Y;

    double ***beam;
    int x, y;
    beam = (double ***)malloc( o.nframes * sizeof(double **) );
    for (i = 0; i < o.nframes; i++)
    {
        beam[i] = (double **)malloc( X * sizeof(double *) );
        for (x = 0; x < X; x++)
        {
            beam[i][x] = (double *)calloc( Y, sizeof(double) );
        }
    }

    // Write the column headers
    print_col_headers( f, o.nframes, tstep );

    for (i = 0; i < o.nphotons; i++)
    {
        // Obtain a random point on the pulsar surface, within the region
        // occupied by the carousel. For Gaussian sparks, choose a 3σ cutoff
        double inner_limit = (psr.csl.type == GAUSSIAN ?
                              o.S_deg - 3.0*o.s_deg :
                              o.S_deg - o.s_deg) * DEG2RAD;

        double outer_limit = (psr.csl.type == GAUSSIAN ?
                              o.S_deg + 3.0*o.s_deg :
                              o.S_deg + o.s_deg) * DEG2RAD;

        if (inner_limit < 0.0)  inner_limit = 0.0;

        random_direction_bounded( &foot_pt_mag, inner_limit, outer_limit,
                0.0, 2.0*PI );

        scale_point( &foot_pt_mag, psr.r, &foot_pt_mag );

        // Convert the foot_pt into observer coordinates
        mag_to_obs_frame( &foot_pt_mag, &psr, NULL, &foot_pt );

        // Check what kind of field line we're on, and how far it is to the
        // extreme point.
        linetype = get_fieldline_type( &foot_pt, &psr, rL_norm, NULL,
                &max_height, NULL );
        if (o.open_only && linetype == CLOSED_LINE)
        {
            continue;
        }

        // Choose a height anywhere up to the max_height for this field line
        height = RAND(max_height);

        // Choose a random frequency
        freq = (RAND(o.f_hi - o.f_lo) + o.f_lo) * 1.0e6;

        // Now climb up the field line to the specified height and emit a
        // photon, pn.
        climb_and_emit_photon( &psr, &foot_pt, height, freq, &pn );

        // Weight the photon (for now) only by the max_height of the field
        // line and the spark profiles
        pulse_number = 0;
        weight = weight_photon_total( &pn, &psr, &foot_pt, height,
                pulse_number, o.g_idx, WEIGHT_TOTAL );
        //weight *= max_height;

        // Calculate which time bin this photon will fall in
        t = psr.P * (pulse_number + pn.phase.deg/360.0);
        tbin = (int)floor( t / tstep );
        while (tbin < 0)           tbin += o.nframes;
        while (tbin >= o.nframes)  tbin -= o.nframes;

        // Bin it up!
        bin_photon_beam( &pn, &psr, weight, xmin, ymin, xbinwidth, ybinwidth,
                X, Y, beam[tbin] );
    }

    // Print out the results
    for (x = 0; x < X; x++)
    {
        for (y = 0; y < Y; y++)
        {
            fprintf( f, "%.15e %.15e", x*xbinwidth+xmin, y*ybinwidth+ymin );
            for (i = 0; i < o.nframes; i++)
                fprintf( f, " %.15e", beam[i][x][y] );
            fprintf( f, "\n" );
        }
        fprintf( f, "\n" );
    }

    // Clean up
    free( o.outfile );

    for (i = 0; i < o.nframes; i++)
    {
        for (x = 0; x < X; x++)
        {
            free( beam[i][x] );
        }
        free( beam[i] );
    }
    free( beam );

    if (o.outfile != NULL) // @suppress("Symbol is not resolved")
        fclose( f );

    return 0;
}

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
    printf( "  -g idx       The power law index for the gamma distribution "
                           "(default = -6.2)\n" );
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
    while ((c = getopt( argc, argv, "a:c:df:g:hn:N:o:OP:s:S:t:4:")) != -1)
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
                break;
            case 'd':
                o->dipole = 1;
                break;
            case 'f':
                parse_range( optarg, &(o->f_lo),
                                     &(o->f_hi),
                                     NULL );
                break;
            case 'g':
                o->g_idx = atof(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'n':
                o->nphotons = atoi(optarg);
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

    if (isnan(o->s_deg) || isnan(o->S_deg) || isnan(o->f_lo))
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


