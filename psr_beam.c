/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program attempts to simulate the emerging beam pattern as
 * observed from a pulsar with given angles α and ζ, period P, and Lorentz
 * factor γ.
 *
 * Here, the idea is that emission is allowed to originate on any open field
 * line. No attempt is made to sample the active field lines realistically.
 * To keep things simple, we start at the magnetic pole and work our way
 * outwards in concentric circles on the polar cap until we reach a "radius"
 * where no open field lines are found.
 *
 * For each open field line, the line is followed outward in incremental
 * steps. At each step, the frequency and direction of an emitted photon is
 * calculated.
 *
 * The final output includes:
 *   1&2) the polar coordinates of the footpoint of the magnetic field line
 *   3&4) the polar coordinates of the photon direction (with retardation)
 *   5) the polarisation angle
 *   6) the phase retardation due to the emission height
 *   7) the critical frequency (in MHz)
 *   8) the emission height (in km)
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
    double  P_sec;       // period, in sec
    double  gamma;       // The Lorentz factor
    char   *outfile;     // name of output file (NULL means stdout)
    double  tmult;       // step size along field lines
    double  s_start;     // starting value of s
    double  s_stop;      // stopping value of s
    double  p_start;     // starting value of p
    double  p_stop;      // stopping value of p
    double  f_start;     // starting value of freq (MHz)
    double  f_stop;      // stopping value of freq (MHz)
    int     open_only;   // only consider open field lines
    int     num_lines;   // sample this many lines
    int     print_los;   // print out the line of sight
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
    o.gamma     = NAN;
    o.outfile   = NULL;
    o.tmult     = 0.01;
    o.s_start   = NAN;
    o.s_stop    = NAN;
    o.p_start   = 0.0;
    o.p_stop    = 360.0;
    o.f_start   = NAN;
    o.f_stop    = NAN;
    o.open_only = 0;
    o.num_lines = 10000;
    o.print_los = 0;

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

    // Write the file header
    print_psrg_header( f, argc, argv );

    // Some needed variables
    int linetype;   // either CLOSED_LINE or OPEN_LINE
    point foot_pt, foot_pt_mag;
    point init_pt;

    // If the line of sight is requested, print it out
    if (o.print_los)
    {
        point LoS, LoS_mag;
        psr_angle LoS_ph;

        // Print out column headers
        fprintf( f, "# phase_deg  th_deg  ph_deg\n" );

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

        // Print out a couple of blank lines to separate it from what follows
        fprintf( f, "\n\n" );
    }

    // Write the column headers
    print_col_headers( f );

    for (i = 0; i < o.num_lines; i++)
    {
        // Obtain a random point on the pulsar surface
        random_direction_bounded( &foot_pt_mag, o.s_start*DEG2RAD,
                o.s_stop*DEG2RAD, o.p_start*DEG2RAD, o.p_stop*DEG2RAD );
        scale_point( &foot_pt_mag, psr.r, &foot_pt_mag );

        // Convert the foot_pt into observer coordinates
        mag_to_obs_frame( &foot_pt_mag, &psr, NULL, &foot_pt );

        // If requested, check that we're on an open field line
        if (o.open_only)
        {
            int rL_norm = 0;
            linetype = get_fieldline_type( &foot_pt, &psr, o.tmult, rL_norm,
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

        climb_and_emit( &psr, &init_pt, o.tmult, o.gamma, o.f_start*1.0e6,
                o.f_stop*1.0e6, f ); /* (1.0e6: convert frequencies to Hz) */
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
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -f  f1:f2    The emission frequency, in MHz. "
                           "The range is from f1 to f2.\n" );
    printf( "  -g  gamma    The Lorentz factor\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -s  s1:s2    The angular distance from the magnetic axis, "
                           "in degrees. The range is from s1 to s2.\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -l           Print out the line of sight first\n" );
    printf( "  -n  nlines   Sample nlines magnetic field lines "
                           "(default: 10000)\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -O           Only consider open field lines (default: off)\n" );
    printf( "  -p  p1:p2    The azimuth relative to the magnetic axis, "
                           "in degrees. The range is from p1 to p2. Ensure "
                           "p1 < p2 [default = 0:360]\n" );
    printf( "  -t  step     Step size for moving along magnetic field lines, "
                           "as a fraction of the light cylinder radius "
                           "(default: 0.01)\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:f:g:hln:o:Op:P:s:S:t:z:")) != -1)
    {
        switch (c)
        {
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'f':
                parse_range( optarg, &(o->f_start),
                                     &(o->f_stop),
                                     NULL );
                break;
            case 'g':
                o->gamma = atof(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'l':
                o->print_los = 1;
                break;
            case 'n':
                o->num_lines = atoi(optarg);
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
        isnan(o->gamma))
    {
        fprintf( stderr, "error: -a, -g, -P and -z options required"
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
    fprintf( f, "#  s_deg  p_deg  th_deg  ph_deg  polangle_deg  retard_deg  "
                "freqcrit_MHz  height_km\n" );
}


