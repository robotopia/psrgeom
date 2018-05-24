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
    double  P_sec;     // period, in sec
    char   *outfile;   // name of output file (NULL means stdout)
    double  p_start;   // starting value of p
    double  p_stop;    // stopping value of p
    int     p_nstep;   // number of p steps
    int     dipole;    // bool: use dipole
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
    o.outfile   = NULL;
    o.p_start   = NAN;
    o.p_stop    = NAN;
    o.p_nstep   = 0;
    o.dipole    = 0;

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
    psr_angle *ze = create_psr_angle_deg( 0.0 ); // This isn't used

    double P = o.P_sec;  // Pulsar spin period
    double r = 1e4;      // Pulsar radius

    set_pulsar( &psr, ra, dec, P, r, al, ze );

    if (o.dipole)
        psr.field_type = DIPOLE;

    point mag_axis;
    psr_angle zero;
    set_psr_angle_deg( &zero, 0.0 );
    set_point_sph( &mag_axis, 1.0, al, &zero, POINT_SET_ALL );

    // Set up points
    point foot_pt;     // The footpoint relative to the ref_axis
    point obs_foot_pt; // The footpoint in the observer frame
    point far_pt;      // The furthest point

    // Print header and column header to file
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Draw the requested lines
    double tmult = 0.001;
    int linetype;  // either CLOSED_LINE or OPEN_LINE
    int p_idx;
    double s_deg, p_deg;
    double s_step;
    double new_s_step;
    double dp;
    psr_angle s, p;

    for (p_idx = 0; p_idx < o.p_nstep; p_idx++)
    {
        // Convert p_idx to an angle
        dp = (o.p_nstep == 1 ?
                0.0 :
                (o.p_stop - o.p_start)/(o.p_nstep - 1.0));
        p_deg = o.p_start + p_idx*dp;
        set_psr_angle_deg( &p, p_deg );

        // Start in the middle (i.e. at the magnetic axis)
        // (we assume the field line whose foot point is at the magnetic axis
        // is an open field line)
        s_deg = 0.0;
        s_step = 1.0;
        new_s_step = s_step / 2.0;

        while (s_deg + s_step > s_deg)
        {
            // Convert s_idx to an angle
            set_psr_angle_deg( &s, s_deg + s_step );

            // From s and p, get a footpoint (rel. to reference axis)
            set_point_sph( &foot_pt, psr.r, &s, &p, POINT_SET_ALL );

            // Rotate it into the observer's frame
            rotate_about_axis( &foot_pt, &obs_foot_pt, &(mag_axis.th),
                    'y', POINT_SET_ALL );
            rotate_about_axis( &obs_foot_pt, &obs_foot_pt, &(mag_axis.ph),
                    'z', POINT_SET_ALL );

            // Find out if the field line is open or closed
            linetype = get_fieldline_type( &obs_foot_pt, &psr, tmult,
                    0, NULL, &far_pt );
            if (linetype == CLOSED_LINE)
            {
                // If closed, halve the trial step size and check again
                s_step = new_s_step;
                new_s_step /= 2.0;
            }
            else
            {
                // If open, take an actual step "forward" and try again
                s_deg += s_step;
            }
        }

        // Print out the (s, p) pair
        fprintf( f, "%.15e %.15e\n", p_deg, s_deg );
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
    printf( "usage: psr_lofl [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -p  p[:P[:n]]   The azimuth relative to the reference axis, "
                           "in degrees. The range is from p to P with n "
                           "steps.\n"
            "                 p      ==> p:p:1\n"
            "                 p:P    ==> p:P:2\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -d           Use a dipole field instead of the default "
                           "Deutsch field\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );

}



void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "a:dho:p:P:")) != -1)
    {
        switch (c)
        {
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'd':
                o->dipole = 1;
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'p':
                parse_range( optarg, &(o->p_start),
                                     &(o->p_stop),
                                     &(o->p_nstep) );
                break;
            case 'P':
                o->P_sec = atof(optarg);
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
    if (isnan(o->al_deg) || isnan(o->P_sec))
    {
        fprintf( stderr, "error: -a and -P options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (isnan(o->p_start))
    {
        fprintf( stderr, "error: -p option required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
{
    // Print out a line to file handle f
    fprintf( f, "# p_deg  s_deg\n" );
}


