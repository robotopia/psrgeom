/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program attempts to locate the polar caps for alpha angles
 * 0 ≤ α ≤ 90°. For each α, it outputs the "distance" (in deg) between
 * the centroid of the polar cap and the magnetic axis, and the "radius"
 * of the polar cap (in deg).
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "psrgeom.h"

void usage()
{
    printf( "PSRGEOM v%s\n", PSRGEOM_VERSION );
    printf( "psr_findcaps -h\n" );
    printf( "    Display this help and exit\n" );
    printf( "psr_findcaps PERIOD\n" );
    printf( "    Find the polar caps for a pulsar with period PERIOD (in "
                "seconds) for alpha angles 0° < α < 90°\n" );
    printf( "\n" );
    printf( "The output (stdout) is three columns of numbers [a,s,da], "
                "where:\n" );
    printf( "    a  = alpha angle\n" );
    printf( "    s  = colatitude with respect to magnetic axis\n" );
    printf( "    da = angular width of polar cap\n" );
    printf( "[a,s,da] are all given in degrees\n" );
}

int main( int argc, char* argv[] )
{
    // See if -h is passed as first argument on the command line
    if (argc <= 1)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    if (strcmp( argv[1], "-h" ) == 0)
    {
        usage();
        exit(EXIT_SUCCESS);
    }

    // Parse the first argument as the period (in seconds)
    // Set step size (along field lines) to 1% of radial distance from origin
    double P = atof( argv[1] );

    double tmult = 0.001;
    int rL_norm = 0;

    // Set up the pulsar struct
    pulsar psr;
    double alpha_deg;
    double r = 1e4;

    psr_angle *ra  = NULL;
    psr_angle *dec = NULL;

    psr_angle al, ze;
    set_psr_angle_deg( &ze, 0.0 ); // This can double up for zero
                                   // as well as zeta

    // Some angles, etc, that we'll need for manipulating start points
    psr_angle start_angle, step_angle;
    point     start_pt, zenith, this_pt, next_pt;
    set_point_sph( &zenith, r+1.0, &ze, &ze, POINT_SET_ALL );
    double s_step;
    double s_low, s_high;
    int linetype;

    double start_s = 0.0; // The first guess of the polar cap location

    // Output a header
    printf( "# alpha_deg  polarcap_colatitude_deg  polarcap_radius_deg\n" );

    // Loop over different alpha angles
    for (alpha_deg = 0.1/P; alpha_deg <= 90.0; alpha_deg += 0.1)
    {
        // For each α, start at the polar cap location of the previous α
        // and scan along the great circle in the xz-plane in both directions
        // for the last open field lines
        set_psr_angle_deg( &al, alpha_deg );
        set_pulsar( &psr, ra, dec, P, r, &al, &ze );

        // Check that the first guess is actually on an open field line
        set_psr_angle_deg( &start_angle, start_s );
        rotate_about_axis( &zenith, &start_pt, &start_angle, 'y', POINT_SET_ALL );
        if (get_fieldline_type( &start_pt, &psr, tmult, rL_norm, NULL, NULL ) == CLOSED_LINE)
        {
            fprintf( stderr, "error: main: initial guess was a closed field line: "
                             "α = %.1f°; s = %.2f°\n", alpha_deg, start_s );
            exit(EXIT_FAILURE);
        }

        // Forward direction, go up 1 degree at a time, until a closed field
        // line is found, then halve step size
        s_step = 1.0;
        copy_point( &start_pt, &this_pt );
        while(s_step >= 1e-6) // I'm happy with this precision
        {
            set_psr_angle_deg( &step_angle, s_step );
            rotate_about_axis( &this_pt, &next_pt, &step_angle, 'y', POINT_SET_ALL );
            linetype = get_fieldline_type( &next_pt, &psr, tmult, rL_norm, NULL, NULL );
            if (linetype == OPEN_LINE)
            {
                copy_point( &next_pt, &this_pt );
            }
            else
            {
                s_step /= 2.0;
            }
        }

        // Extract the answer for the "forward" direction
        if (fabs(this_pt.ph.deg) <= 1e-3)
        {
            s_high = this_pt.th.deg;
        }
        else if (fabs(this_pt.ph.deg - 180.0) <= 1e-3)
        {
            s_high = -this_pt.th.deg;
        }
        else
        {
            fprintf( stderr, "error: main: point somehow got off the xz-plane\n" );
            exit(EXIT_FAILURE);
        }

        // Now do the same thing for the "backward" direction
        s_step = -1.0;
        copy_point( &start_pt, &this_pt );
        while(s_step <= -1e-6) // I'm happy with this precision
        {
            set_psr_angle_deg( &step_angle, s_step );
            rotate_about_axis( &this_pt, &next_pt, &step_angle, 'y', POINT_SET_ALL );
            linetype = get_fieldline_type( &next_pt, &psr, tmult, rL_norm, NULL, NULL );
            if (linetype == OPEN_LINE)
            {
                copy_point( &next_pt, &this_pt );
            }
            else
            {
                s_step /= 2.0;
            }
        }

        // Extract the answer for the "forward" direction
        if (fabs(this_pt.ph.deg) <= 1e-3)
        {
            s_low = this_pt.th.deg;
        }
        else if (fabs(this_pt.ph.deg - 180.0) <= 1e-3)
        {
            s_low = -this_pt.th.deg;
        }
        else
        {
            fprintf( stderr, "error: main: point somehow got off the xz-plane\n" );
            exit(EXIT_FAILURE);
        }

        // Print out the results for this α
        start_s = (s_low + s_high) / 2.0; // = new start point
        printf( "%12f %12f %12f\n", alpha_deg, start_s - alpha_deg, (s_high - s_low)/2.0 );
    }
}
