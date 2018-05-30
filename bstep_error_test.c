#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "psrgeom.h"

int choose_nsteps()
/* Choose a random number of steps drawn from the values
 * {2, 4, 8, 16, ... 1024}
 */
{
    return 0x10 << (rand() % 12);
}

int main()
{
    // Seed the random number generator
    srand( time( NULL ) );

    int s;           // A counter for steps
    int n, N = 1000; // <-- number of points in magnetosphere to test
    double rad_curv; // The estimated radius of curvature at each point

    double min_step = 100.0; // The "baseline" step size, with which to
                             // compare errors
    double step; // The actual step size taken
    int nsteps; // The number of 'min_step's to take for testing
    double error;

    // Set up a pulsar (with the default deutsch field)
    pulsar psr;

    psr_angle *ra  = NULL;
    psr_angle *dec = NULL;

    psr_angle al, ze;
    set_psr_angle_deg( &al, 0.0 );
    set_psr_angle_deg( &ze, 0.0 );

    double P = 1.0;
    double r = 1e4; /* This doesn't actually make a difference to any of the
                       outputs of this program */

    set_pulsar( &psr, ra, dec, P, r, &al, &ze );

    // Set up points
    point init_pt; // The starting point
    point end_pt1; // The point reached by taking small steps
    point end_pt2; // The point reached by taking one large step
    point diff;    // end_pt1 - end_pt2

    // Test this many points
    for (n = 0; n < N; n++)
    {
        // Generate a random point within 1/2 the light cylinder radius
        random_point_in_lightcyl( &init_pt, &psr, 0.5, 1.0 );

        // Estimate the rad. of curv. at that point (as for dipole field)
        rad_curv = 1.0/calc_dipole_curvature( &init_pt );

        // Choose step size to test
        nsteps = choose_nsteps();
        step   = (double)nsteps * min_step;

        // Get a baseline measurement by taking several small steps
        copy_point( &init_pt, &end_pt1 );
        for (s = 0; s < nsteps; s++)
            Bstep( &end_pt1, &psr, min_step, DIR_OUTWARD, &end_pt1 );

        // Compare this to the point reached by taking one large step
        Bstep( &init_pt, &psr, step, DIR_OUTWARD, &end_pt2 );

        // Calculate the distance between the two point (= the error)
        set_point_xyz( &diff, end_pt1.x[0] - end_pt2.x[0],
                              end_pt1.x[1] - end_pt2.x[1],
                              end_pt2.x[2] - end_pt2.x[2],
                              POINT_SET_R );
        error = diff.r;

        // Print out a line with the curvature, the step size, and the error
        printf( "%.15e %.f %.15e\n", rad_curv, step, error );
    }
}
