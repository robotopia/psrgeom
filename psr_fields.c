#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "psrgeom.h"

// Define arbitrary numbers for grid types
#define  GT_CART  1
#define  GT_CYL2  2
#define  GT_CYL3  3

struct opts
{
    double  al_deg;    // alpha angle in deg
    double  ze_deg;    // zeta angle in deg
    double  P_sec;     // period, in sec
    char   *format;    // output format string
    int     grid_type; // grid type (cart, cyl2, or cyl3)
    int     rL_norm;   // bool: normalise to light cylinder radius?
    int     npoints;   // number of points on one side of grid
    char   *outfile;   // name of output file (NULL means stdout)
    double  rho_max;   // largest rho to consider for gridpoints
    double  vsize;     // ratio of vector size to grid cell size
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );

int main( int argc, char *argv[] )
{
    // Set up struct for command line options and set default values
    struct opts o;
    o.al_deg    = NAN;
    o.P_sec     = NAN;
    o.ze_deg    = NAN;
    o.format    = NULL;
    o.grid_type = GT_CART;
    o.rL_norm   = 0;
    o.npoints   = 0; /* Set to 0 now, because the true default depends on
                        value of grid_type */
    o.outfile   = NULL;
    o.rho_max   = 1.0;
    o.vsize     = 0.5;

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
            fprintf( stderr, "error: could open file %s\n", o.outfile );
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

    // Set up grid point
    point X;

    // Set up particle speed
    double v = SPEED_OF_LIGHT;

    // Get the number of velocity solutions
    int nsols;

    // Figure our grid spacing
    double dx     = 2.0 * o.rho_max * psr.RL / (double)(npoints - 1);
    double dtheta = 2.0 * o.rho_max * PI     / (double)(npoints);

    point B;
    point V1, V2;
    point A1, A2;

    print_psrg_header( f, argc, argv );
    fprintf( f, "#x y z Bx By Bz V1x V1y V1z V2x V2y V2z A1x A1y A1z A2x A2y A2z\n" );

    int i, j, k;
    for (i = -1; i <= 1; i += 0.2)
    for (j = -1; j <= 1; j += 0.2)
    for (k = -1; k <= 1; k += 0.2)
    {
        set_point_xyz( &X, (double)i*psr.rL,
                           (double)j*psr.rL,
                           (double)k*psr.rL,
                           POINT_SET_ALL );

        if (X.rhosq > psr.rL2) continue;

        calc_fields( &X, &psr, v, &B, &V1, &V2, &A1, &A2, &nsols );

        if (nsols > 0)
        {
            fprintf( f, "%lf %lf %lf  ", i, j, k );
            fprintf( f, "%.15e %.15e %.15e  ",  B.x[0],  B.x[1],  B.x[2] );
            fprintf( f, "%.15e %.15e %.15e  ", V1.x[0], V1.x[1], V1.x[2] );
            fprintf( f, "%.15e %.15e %.15e  ", V2.x[0], V2.x[1], V2.x[2] );
            fprintf( f, "%.15e %.15e %.15e  ", A1.x[0], A1.x[1], A1.x[2] );
            fprintf( f, "%.15e %.15e %.15e\n", A2.x[0], A2.x[1], A2.x[2] );
        }
    }

    // Clean up
    destroy_psr_angle( ra  );
    destroy_psr_angle( dec );
    destroy_psr_angle( al  );
    destroy_psr_angle( ze  );

    free( o.format );
    free( o.outfile );

    if (o.outfile != NULL)
        fclose( f );

    return 0;
}

void usage()
{
    printf( "usage: psr_fields [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -f  format   The output format string. For complete list of "
                           "format specifiers, see the man page. (default: "
                           "\"Px Py Pz Bx By Bz\")\n" );
    printf( "  -g  gt       The grid type; \"gt\" can be \"cart\", \"cyl2\", "
                           "or \"cyl3\" (default: \"cart\")\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -N  npoints  The number of points along a side of the grid.\n"
            "               With grid type \"cart\" or \"cyl3\", this "
                           "effectively sets the grid spacing size to be "
                           "2*rho*rL/(npoints - 1).\n"
            "               With grid type \"cyl2\", the grid spacing is set "
                           "to be 2*pi*rho/npoints.\n"
            "               (default: 11 for cart/cyl3, 24 for cyl2)\n");
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -r  rho      The largest rho to consider, as a fraction of "
                           "the light cylinder radius, i.e. rho = sqrt("
                           "x^2 + y^2)/rL (default: 1.0)\n" );
    printf( "  -v  vsize    Normalise vector outputs so that their length "
                           "is \"vsize\" times the distance between grid "
                           "points (default: 0.5)\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:f:g:hLN:o:P:r:v:z:")) != -1)
    {
        switch (c)
        {
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'f':
                o->format = strdup(optarg);
                break;
            case 'g':
                if (!strcmp( optarg, "cart" ))
                    o->grid_type = GT_CART;
                else if (!strcmp( optarg, "cyl2" ))
                    o->grid_type = GT_CYL2;
                else if (!strcmp( optarg, "cyl3" ))
                    o->grid_type = GT_CYL3;
                else
                {
                    fprintf( stderr, "error: unrecognised grid type \"%s\"",
                                     optarg );
                    exit(EXIT_FAILURE);
                }
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'L':
                o->rL_norm = 1;
                break;
            case 'N':
                o->npoints = atoi(optarg);
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 'r':
                o->rho_max = atof(optarg);
                break;
            case 'v':
                o->vsize = atof(optarg);
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

    if (o->format == NULL)
        o->format = strdup("Px Py Pz Bx By Bz");

    if (o->npoints == 0)
    {
        if (o->grid_type == GT_CART || o->grid_type == GT_CYL3)
            o->npoints = 11;
        else
            o->npoints = 24;
    }
}
