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
    double  freq_MHz;  // The observing frequency
    double  tstep;     // RK4 step size, as a fraction of lt cyl radius
    int     rL_norm;   // bool: normalise output to light cylinder radius?
    char   *outfile;   // name of output file (NULL means stdout)
    double  s_deg;     // the spark angular size
    double  S_deg;     // the carousel angular radius
    int     nsparks;   // number of sparks in carousel
    double  P4_sec;    // the rotation time of the carousel
    int     csl_type;  // the type of spark profile (TOPHAT or GAUSSIAN)
    int     dipole;    // bool: use dipole
    int     nlines;    // number of lines to sample
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
    o.freq_MHz  = NAN;
    o.tstep     = MAX_BSTEP;
    o.rL_norm   = 0;
    o.outfile   = NULL;
    o.s_deg     = NAN;
    o.S_deg     = NAN;
    o.P4_sec    = NAN;
    o.csl_type  = GAUSSIAN;
    o.dipole    = 0;
    o.nlines    = 1000;
    o.nsparks   = 0;

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

    double P = o.P_sec;  // Pulsar spin period
    double r = 1e4;      // Pulsar radius

    set_pulsar( &psr, NULL, NULL, P, r, &al, &ze );

    if (o.dipole)
        psr.field_type = DIPOLE;

    // Set up pulsar's carousel
    psr_angle s, S;
    set_psr_angle_deg( &s, o.s_deg );
    set_psr_angle_deg( &S, o.S_deg );
    set_pulsar_carousel( &psr, o.nsparks, &s, &S, o.csl_type, o.P4_sec );
    double t;

    // Convert the observing frequency to Hz
    double freq = o.freq_MHz * 1.0e6;

    // Scale the tstep size according to light cylinder radius
    o.tstep *= psr.rL;
    double dist; // Keep track of the total distance travelled
    double step; // Each step size

    // Set up needed structs
    point emit_pt;     // The footpoint relative to the ref_axis
    photon pn;         // An emitted photon

    // Print header and column header to file
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Write out the requested lines, and their associated values
    int i;
    for (i = 0; i < o.nlines; i++)
    {
        // Choose a foot point from among the sparks
        t = 0.0; // Just use time = 0 exclusively for now -- i.e.
                 // ignore carousel rotation
        random_spark_footpoint( &emit_pt, &psr, t );

        // Reset the distance counter
        dist = 0.0;

        // Climb the line at regular intervals
        do
        {
            // Calculate the properties of a photon emitted at this point
            emit_pulsar_photon( &psr, &emit_pt, freq, &pn );

            // Print out the next line
            fprintf( f, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
                    emit_pt.r, dist, sqrt(emit_pt.rhosq), pn.curvature,
                    pn.retarded_LoS.th.deg, pn.gamma, pn.power );

            // Climb a bit more
            if (emit_pt.r < 2.0*MAX_BSTEP*psr.rL)
                step = 0.5*emit_pt.r;
            else
                step = o.tstep;

            Bstep( &emit_pt, &psr, step, DIR_OUTWARD, &emit_pt );
            set_point_xyz( &emit_pt, emit_pt.x[0], emit_pt.x[1], emit_pt.x[2],
                    POINT_SET_ALL );
            dist += step;
        }
        while ((emit_pt.rhosq < psr.rL2) && (emit_pt.r > psr.r));

        fprintf( f, "\n" );
    }

    // Clean up
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
    printf( "  -f  freq     The observing frequency, in MHz\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -s  s_deg    The angular size of the sparks (in deg)\n" );
    printf( "  -S  S_deg    The carousel's angular radius (in deg)\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "  -4  P4       The carousel's rotation period (in sec)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -c type      The spark profile type, either GAUSSIAN (default) "
                           "or TOPHAT\n" );
    printf( "  -d           Use a dipole field instead of the default "
                           "Deutsch field\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -n  lines    Sample this many field lines (default: 10000)\n" );
    printf( "  -N  nsparks  The number of sparks in the carousel. If nsparks "
                           "= 0 (default), the carousel is a \"solid\" "
                           "annulus.\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -t  tstep    The size of the RK4 steps, as a fraction of the "
                           "light cylinder radius (default: 0.005)\n" );

}



void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:c:df:hLn:N:o:P:s:S:t:z:4:")) != -1)
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
                o->freq_MHz = atof(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'L':
                o->rL_norm = 1;
                break;
            case 'n':
                o->nlines = atoi(optarg);
                break;
            case 'N':
                o->nsparks = atoi(optarg);
                break;
            case 'o':
                o->outfile = strdup(optarg);
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
                o->tstep = atof(optarg);
                break;
            case 'z':
                o->ze_deg = atof(optarg);
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
    if (isnan(o->al_deg) || isnan(o->P_sec) || isnan(o->ze_deg) ||
        isnan(o->P4_sec))
    {
        fprintf( stderr, "error: -a, -P, -z and -4 options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (isnan(o->s_deg) || isnan(o->S_deg) || isnan(o->freq_MHz))
    {
        fprintf( stderr, "error: -f, -S and -s options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# radial_height  line_height  perp_height  curvature  "
            "beam_opening_angle_deg  gamma  radiated_power\n" );
}


