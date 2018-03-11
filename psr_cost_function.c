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
    double  tmult;     // RK4 step size, as a fraction of lt cyl radius
    int     rL_norm;   // bool: normalise to light cylinder radius?
    int     direction; // either DIR_OUTWARD or DIR_INWARD
    char   *outfile;   // name of output file (NULL means stdout)
    int     n[3];      // The number of samples in X,Y,Z
    double  lim[3][2]; // The lower and upper limits of X
    char    XYZ[3];    // The loop order
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
    o.tmult     = 0.01;
    o.rL_norm   = 0;
    o.direction = DIR_OUTWARD;
    o.outfile   = NULL;
    int i, j;
    for (i = 0; i < 3; i++)
    {
        o.n[i]   = -1;
        o.XYZ[i] = '\0';
    }
    for (i = 0; i < 3; i++)
    for (j = 0; j < 2; j++)
        o.lim[i][j] = NAN;

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

    // Set up point for sampling
    point X;

    // Set up the spacing between pixels
    double dX[3];
    for (i = 0; i < 3; i++)
        dX[i] = (o.lim[i][1] - o.lim[i][0]) / (double)o.n[i];

    // Write the file and column headers
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Iterate through the three dimensions, covering all pixels
    double a;        // A temporary placeholder for either x, y, or z
    double x, y, z;  // The location values
    double c1, c2;   // The cost functions
    int xi[3];       // The pixel indices

    x = y = z = NAN;

    for (xi[0] = 0; xi[0] < o.n[0]; xi[0]++)
    {
        a = xi[0]*dX[0] + o.lim[0][0];
        if      (o.XYZ[0] == 'X')  x = a; 
        else if (o.XYZ[0] == 'Y')  y = a;
        else if (o.XYZ[0] == 'Z')  z = a;

        for (xi[1] = 0; xi[1] < o.n[1]; xi[1]++)
        {
            a = xi[1]*dX[1] + o.lim[1][0];
            if      (o.XYZ[1] == 'X')  x = a; 
            else if (o.XYZ[1] == 'Y')  y = a;
            else if (o.XYZ[1] == 'Z')  z = a;

            for (xi[2] = 0; xi[2] < o.n[2]; xi[2]++)
            {
                a = xi[2]*dX[2] + o.lim[2][0];
                if      (o.XYZ[2] == 'X')  x = a; 
                else if (o.XYZ[2] == 'Y')  y = a;
                else if (o.XYZ[2] == 'Z')  z = a;

                // Set up point
                set_point_xyz( &X, x, y, z, POINT_SET_ALL );

                // Evaluate cost functions
                c1 = psr_cost_lofl( &X, &psr );
                c2 = psr_cost_los( &X, &psr, ph, o.direction );

                // Write out the result
                fprintf( f, "%.15e %.15e %.15e %.15e %.15e\n",
                        x, y, z, c1, c2 );
            }
            fprintf( f, "\n" );
        }
        fprintf( f, "\n" );
    }

    // Clean up
    destroy_psr_angle( ra  );
    destroy_psr_angle( dec );
    destroy_psr_angle( al  );
    destroy_psr_angle( ze  );
    destroy_psr_angle( ph  );

    free( o.outfile );

    if (o.outfile != NULL)
        fclose( f );

    return EXIT_SUCCESS;
}

void usage()
{
    printf( "usage: psr_cost_function [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -p  phase    The rotation phase of the observed emission, in "
                           "degrees (required)\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "  -X  x1,x2,n  Calculate the cost function in the range "
                           "[x1,x2] with n samples\n" );
    printf( "  -Y  y1,y2,n  Calculate the cost function in the range "
                           "[y1,y2] with n samples\n" );
    printf( "  -Z  z1,z2,n  Calculate the cost function in the range "
                           "[z1,z2] with n samples\n" );
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
    printf( "  -t  tmult    The initial size of the RK4 steps, as a fraction "
                           "of the light cylinder radius (default: 0.01)\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    int l = 0; // For determining the XYZ loop order
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:hiLo:p:P:r:t:X:Y:z:Z:")) != -1)
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
            case 't':
                o->tmult = atof(optarg);
                break;
            case 'X':
            case 'Y':
            case 'Z':
                sscanf( optarg, "%lf,%lf,%d",
                        &(o->lim[l][0]), &(o->lim[l][1]), &(o->n[l]) );
                o->XYZ[l++] = c;
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

    int i, j;
    for (i = 0; i < 3; i++)
        if (o->n[i] <= 0)
        {
            fprintf( stderr, "error: -n must be present and its arguments "
                             "must be positive\n" );
            exit(EXIT_FAILURE);
        }

    for (i = 0; i < 3; i++)
    for (j = 0; j < 2; j++)
    {
        if (isnan( o->lim[i][j] ))
        {
            fprintf( stderr, "error: failed to read/parse required -X,-Y,-Z "
                             "arguments\n" );
            exit(EXIT_FAILURE);
        }
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# x  y  z  psr_cost_lofl  psr_cost_los\n" );
}


