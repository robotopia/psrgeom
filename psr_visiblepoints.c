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
    double  P_sec;     // period, in sec
    double  dph_deg;   // step size for rotation phase, in deg
    double  dr;        // step size for radius, in fraction of rL
    int     rL_norm;   // bool: normalise to light cylinder radius?
    int     direction; // either DIR_OUTWARD or DIR_INWARD
    char   *outfile;   // name of output file (NULL means stdout)
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
    o.dph_deg   = 1.0;
    o.dr        = 0.01;
    o.rL_norm   = 0;
    o.direction = DIR_OUTWARD;
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
    double r = 1e4; /* This will be used later as we move outwards from the
                       pulsar surface */

    set_pulsar( &psr, ra, dec, P, r, al, ze );

    // Set scaling factor in case user requests result in units of rL
    double xscale = (o.rL_norm ? 1.0/psr.rL : 1.0);

    // Write the file and column headers
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Calculate the polarisation angle for each phase
    point emit_pt, dip_emit_pt, init_pt, LoS, LoS_mag;
    psr_angle ph, th;

    double ph_deg;
    for (ph_deg = 0.0; ph_deg < 360.0; ph_deg += o.dph_deg)
    {
        // Set the rotation phase
        set_psr_angle_deg( &ph, ph_deg );

        // Find the line of sight
        line_of_sight( &psr, &ph, &LoS );

        // Get the line of sight in the magnetic frame
        obs_to_mag_frame( &LoS, &psr, &ph, &LoS_mag );

        // Get the (dipole) position angle corresponding to that line of sight
        beamangle_to_posangle( &(LoS_mag.th), &th );

        // Construct the emission point for the dipole case
        // (in the magnetic frame)
        set_point_sph( &dip_emit_pt, psr.r,
                                     &th,
                                     &(LoS_mag.ph),
                                     POINT_SET_ALL );

        // And rotate it back to the observer frame
        mag_to_obs_frame( &dip_emit_pt, &psr, &ph, &init_pt );

        // This will now be used as the initial point for finding the visible
        // point at radius = stellar radius
        while (init_pt.r < psr.rL)
        {
            // Find the "true" visible point at this radius
            find_LoS_at_r( &init_pt, &psr, &ph,
                           o.direction, &emit_pt, NULL );

            // Print it out (x,y,z,φ)!
            fprintf( f, "%.15e %.15e %.15e %.15e\n",
                        emit_pt.x[0] * xscale,
                        emit_pt.x[1] * xscale,
                        emit_pt.x[2] * xscale,
                        ph.deg );

            // Move up to the next radius
            set_point_sph( &dip_emit_pt, emit_pt.r + o.dr,
                                         &emit_pt.th,
                                         &emit_pt.ph,
                                         POINT_SET_ALL );
        }

        // Put in a blank line to separate different phases
        fprintf( f, "\n" );
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
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
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
    printf( "  -p  step     Step size for rotation phase, in degrees "
                           "(default: 1.0)\n" );
    printf( "  -r  step     Step size for radius, as a fraction of the light"
                           "cylinder radius (default: 0.01)\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:hiLo:p:P:r:z:")) != -1)
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
                o->dph_deg = atof(optarg);
                break;
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 'r':
                o->dr = atof(optarg);
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
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# x y z φ(deg)\n" );
}


