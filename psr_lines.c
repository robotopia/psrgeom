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
    double  tmult;     // RK4 step size, as a fraction of lt cyl radius
    char   *outfile;   // name of output file (NULL means stdout)
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
    o.tmult     = 0.01;
    o.outfile   = NULL;

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

    // Set up grid point
    point X;

    print_psrg_header( f, argc, argv );

    // Set the point at this gridpoint
    set_point_xyz( &X, 0.4*psr.rL,
                       0.4*psr.rL,
                       0.4*psr.rL,
                       POINT_SET_ALL );

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
                           "\"Xx Xy Xz Bx By Bz\")\n" );
    printf( "  -g  gt       The grid type; \"gt\" can be \"cart\", \"cyl2\", "
                           "or \"cyl3\" (default: \"cart\")\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -l           Limit to points within light cylinder\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -N  npoints  The number of points along a side of the grid.\n"
            "               With grid type \"cart\", this effectively sets "
                           "the grid spacing size to be 2*rho*rL/(npoints - "
                           "1).\n"
            "               With grid type \"cyl2\" or \"cyl3\", the grid "
                           "spacing is set to be 2*pi*rho/npoints.\n"
            "               (default: 11 for cart, 24 for cyl2/cyl3)\n");
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
    while ((c = getopt( argc, argv, "a:f:g:hlLN:o:P:r:v:z:")) != -1)
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
            case 'l':
                o->rL_lim = 1;
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
        o->format = strdup("Xx Xy Xz Bx By Bz");

    if (o->npoints == 0)
    {
        if (o->grid_type == GT_CART)
            o->npoints = 11;
        else
            o->npoints = 24;
    }
}


void parse_format( char *format, struct tokens *tok )
{
    char c;
    int char_num = 0;
    tok->n = 0;
    do
    {
        c = *format;
        if (char_num == 0) // Parse the first character
        {
            switch (c)
            {
                case 'A':
                    tok->calcA = 1;
                    __attribute__ ((fallthrough));
                case 'V':
                    tok->calcV = 1;
                    __attribute__ ((fallthrough));
                case 'B':
                    tok->calcB = 1;
                    __attribute__ ((fallthrough));
                case 'X':
                    tok->token[tok->n][char_num] = c;
                    char_num++;
                    break;
                default:
                    // char_num stays at 0
                    break;
            }
            format++;
        }
        else if (char_num == 1) // Parse the 2nd character
        {
            switch (c)
            {
                case 'x':
                case 'y':
                case 'z':
                    tok->token[tok->n][char_num] = c;
                    if (tok->token[tok->n][0] == 'X' ||
                        tok->token[tok->n][0] == 'B')
                    {
                        char_num = 0;
                        tok->n++;
                    }
                    else
                        char_num++;
                    format++;
                    break;
                case 'r':
                case 'h':
                case 't':
                case 'p':
                    if (tok->token[tok->n][0] == 'X')
                    {
                        tok->token[tok->n][char_num] = c;
                        tok->n++;
                        format++;
                    }
                    char_num = 0;
                    break;
                case 'd': // only allowed after 'B'
                    if (tok->token[tok->n][0] == 'B')
                    {
                        tok->token[tok->n][char_num] = c;
                        char_num++;
                        format++;
                    }
                    else
                        char_num = 0;
                    break;
                default:
                    char_num = 0;
                    break;
            }
        }
        else if (char_num == 2) // Parse the 3rd character
        {
            switch (c)
            {
                case 'P':
                case 'N':
                    if (tok->token[tok->n][0] == 'V')
                    {
                        tok->calcV = 1;
                        tok->calcB = 1;
                        tok->token[tok->n][char_num] = c;
                        tok->n++;
                        format++;
                    }
                    else if (tok->token[tok->n][0] == 'A')
                    {
                        tok->calcA = 1;
                        tok->calcV = 1;
                        tok->calcB = 1;
                        tok->token[tok->n][char_num] = c;
                        tok->n++;
                        format++;
                    }
                    break;
                case 'R':
                    if (tok->token[tok->n][0] == 'B' &&
                        tok->token[tok->n][1] == 'd')
                    {
                        tok->token[tok->n][char_num] = c;
                        tok->n++;
                        tok->calcB = 1;
                        format++;
                    }
            }
            char_num = 0;
        }
    } while (c != '\0' && tok->n < MAX_NTOKENS);
}


void print_col_headers( FILE *f, char *format )
{
    // Print out a line to file handle f
    fprintf( f, "# %s\n", format );
}


void print_token_value( FILE *f, struct tokens *tok, int n, point *X,
                        point *B, point *V1, point *V2, point *A1, point *A2,
                        double BdR, double xscale, double vscale )
{
    // Print X
    if (!strncmp( tok->token[n], "Xx", 2 ))
        fprintf( f, "%.15e", xscale*X->x[0] );
    else if (!strncmp( tok->token[n], "Xy", 2 ))
        fprintf( f, "%.15e", xscale*X->x[1] );
    else if (!strncmp( tok->token[n], "Xz", 2 ))
        fprintf( f, "%.15e", xscale*X->x[2] );
    else if (!strncmp( tok->token[n], "Xr", 2 ))
        fprintf( f, "%.15e", xscale*X->r );
    else if (!strncmp( tok->token[n], "Xh", 2 ))
        fprintf( f, "%.15e", xscale*sqrt(X->rhosq) );
    else if (!strncmp( tok->token[n], "Xt", 2 ))
        fprintf( f, "%.15e", X->th.rad );
    else if (!strncmp( tok->token[n], "Xp", 2 ))
        fprintf( f, "%.15e", X->ph.rad );

    // Print B's
    else if (!strncmp( tok->token[n], "Bx", 2 ))
        fprintf( f, "%.15e", vscale*B->x[0] );
    else if (!strncmp( tok->token[n], "By", 2 ))
        fprintf( f, "%.15e", vscale*B->x[1] );
    else if (!strncmp( tok->token[n], "Bz", 2 ))
        fprintf( f, "%.15e", vscale*B->x[2] );

    // Print V's
    else if (!strncmp( tok->token[n], "VxP", 3 ))
        fprintf( f, "%.15e", vscale*V1->x[0] );
    else if (!strncmp( tok->token[n], "VyP", 3 ))
        fprintf( f, "%.15e", vscale*V1->x[1] );
    else if (!strncmp( tok->token[n], "VzP", 3 ))
        fprintf( f, "%.15e", vscale*V1->x[2] );
    else if (!strncmp( tok->token[n], "VxN", 3 ))
        fprintf( f, "%.15e", vscale*V2->x[0] );
    else if (!strncmp( tok->token[n], "VyN", 3 ))
        fprintf( f, "%.15e", vscale*V2->x[1] );
    else if (!strncmp( tok->token[n], "VzN", 3 ))
        fprintf( f, "%.15e", vscale*V2->x[2] );

    // Print A's
    else if (!strncmp( tok->token[n], "AxP", 3 ))
        fprintf( f, "%.15e", vscale*A1->x[0] );
    else if (!strncmp( tok->token[n], "AyP", 3 ))
        fprintf( f, "%.15e", vscale*A1->x[1] );
    else if (!strncmp( tok->token[n], "AzP", 3 ))
        fprintf( f, "%.15e", vscale*A1->x[2] );
    else if (!strncmp( tok->token[n], "AxN", 3 ))
        fprintf( f, "%.15e", vscale*A2->x[0] );
    else if (!strncmp( tok->token[n], "AyN", 3 ))
        fprintf( f, "%.15e", vscale*A2->x[1] );
    else if (!strncmp( tok->token[n], "AzN", 3 ))
        fprintf( f, "%.15e", vscale*A2->x[2] );

    // Print B dot Rxy
    else if (!strncmp( tok->token[n], "BdR", 3 ))
        fprintf( f, "%.15e", BdR );

}



void print_all_tokens( FILE *f, struct tokens *tok, point *X, pulsar *psr,
                       double xscale, double vscale )
{
    point B, V1, V2, A1, A2;
    int nsols;
    double BdR = NAN;
    double v = SPEED_OF_LIGHT;

    if (tok->calcA)
        calc_fields( X, psr, v, &B, &V1, &V2, &A1, &A2, &nsols );
    else if (tok->calcV)
        calc_fields( X, psr, v, &B, &V1, &V2, NULL, NULL, &nsols );
    else if (tok->calcB)
        calc_fields( X, psr, v, &B, NULL, NULL, NULL, NULL, &nsols );

    if (tok->calcA || tok->calcV || tok->calcB)
    {
        BdR = B.x[0] * X->ph.cos +
              B.x[1] * X->ph.sin;
    }

    int n;
    for (n = 0; n < tok->n; n++)
    {
        print_token_value( f, tok, n, X, &B, &V1, &V2, &A1, &A2,
                BdR, xscale, vscale );
        fprintf( f, " " );
    }
    fprintf( f, "\n" );
}
