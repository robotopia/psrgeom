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
    double  tstep;     // RK4 step size, in seconds
    double  ttotal;    // the total time
    int     rL_norm;   // bool: normalise to light cylinder radius?
    int     direction; // Direction of particles along magnetic field lines
    double  rho_lim;   // limit to within rho_lim * light cylinder radius?
    char   *outfile;   // name of output file (NULL means stdout)
    double  X[3];
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
    o.tstep     = NAN;
    o.ttotal    = NAN;
    o.rL_norm   = 0;
    o.direction = DIR_OUTWARD;
    o.rho_lim   = 1.2;
    o.outfile   = NULL;
    o.X[0]      = NAN;
    o.X[1]      = NAN;
    o.X[2]      = NAN;

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

    if (isnan(o.tstep))
        o.tstep = 0.01*psr.rL/SPEED_OF_LIGHT;

    if (isnan(o.ttotal))
        o.ttotal = psr.P;

    // Print header and column header to file
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Set up points
    point X;
    double xscale = (o.rL_norm ? 1.0/psr.rL : 1.0);
    set_point_xyz( &X, o.X[0], o.X[1], o.X[2], POINT_SET_ALL );
    fprintf( f, "%.15e %.15e %.15e %.15e\n",
                0.0, xscale*o.X[0], xscale*o.X[1], xscale*o.X[2] );

    // Start at t = 0 and step along the velocity field
    double t;

    for (t = o.tstep; t < o.ttotal; t += o.tstep)
    {
        Vstep( &X, &psr, o.tstep*SPEED_OF_LIGHT, o.direction, &X );
        fprintf( f, "%.15e %.15e %.15e %.15e\n",
                    t, xscale*X.x[0], xscale*X.x[1], xscale*X.x[2] );
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
    printf( "  -X  x,y,z    Starting point, in meters (required)\n" );
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
    printf( "  -r  rho_lim  The maximum distance allowed for lines. "
                           "(default = 1.2)\n" );
    printf( "  -t  tstep    The initial size of the RK4 steps, in seconds "
                           "of light travel time (default: light travel time "
                           "of 1%% of light cylinder radius)\n" );
    printf( "  -T  ttotal   The total time to integrate for, in seconds "
                           "(default: 1x spin period)\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:hiLo:P:r:t:T:X:z:")) != -1)
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
            case 'r':
                o->rho_lim = atof(optarg);
                break;
            case 't':
                o->tstep = atof(optarg);
                break;
            case 'T':
                o->ttotal = atof(optarg);
                break;
            case 'X':
                sscanf( optarg, "%lf,%lf,%lf", &o->X[0], &o->X[1], &o->X[2] );
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

    if (isnan(o->X[0]) || isnan(o->X[1]) || isnan(o->X[2]))
    {
        fprintf( stderr, "error: -X option required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# t_(sec) x y z\n" );
}


