/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the functions for calculating the magnetic
 *   field.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psrgeom.h"

void quad_solve_real( double a, double b, double c, double *x1, double *x2,
                      int *nsols )
/* Solve a quadratic equation ax^2 + bx + c = 0. If there are no real
 * solutions, neither x1 nor x2 are set. If there is one real solution,
 * it is put into x1.
 *
 * Inputs:
 *   double a, b, c  = quadratic coefficients
 * Outputs:
 *   double *x1, *x2 = real solutions
 *   int *nsols      = the number of real solutions (can be 0, 1, 2)
 */
{
    double b2  = b*b;
    double d   = 0.5/a;
    double det = sqrt(b2 - 4.0*a*c);

    if (det < 0.0)
        *nsols = 0;
    else if (det == 0.0)
    {
        *nsols = 1;
        *x1 = (b2 + det) * d;
    }
    else // (det > 0.0)
    {
        *nsols = 2;
        *x1 = (b2 + det) * d;
        *x2 = (b2 - det) * d;
    }
}

void Vfield( point *x, pulsar *psr, double v, point *V1, point *V2,
             int *nsols )
/* This returns the normalised vector V (cast as a "point"), i.e. the velocity
 * of a particle instantaneously at point x (observer frame) when the magnetic
 * axis is instantaneously in the xz-plane. It is understood that the
 * azimuthal component of the velocity is fixed by the rotation of the pulsar
 * and that the total velocity of the particle is v.
 *
 * This function assumes that both the Cartesian and spherical coordinates of
 * x are set.
 *
 * Inputs:
 *   point  *x      = the location of the particle in observer coordinates
 *   pulsar *psr    = the pulsar struct
 *   double  v      = the speed of the particle in the observer frame
 *
 * Outputs:
 *   int    *nsols  = the number of solutions found (0, 1, or 2)
 *   point  *V1     = the normalised velocity vector of the particle (sol #1)
 *   point  *V2     = the normalised velocity vector of the particle (sol #1)
 */
{
    // The vectors we will need (cast as points)
    point B;    // The magnetic field
    point az_V; // The azimuthal velocity

    // First, get the magnetic field direction at x
    Bfield( x, psr, &B );

    // Now, calculate the azimuth velocity at x due to the pulsar's rotation
    double rho = x->r * x->th.sin; // The perp dist from the z-axis
    double az_v = rho * psr->Om.rad;
    psr_angle ph_V;
    set_psr_angle_deg( &ph_V, x->ph.deg + 90.0 );
    set_point_cyl( &az_V, 1.0, &ph_V, 0.0, POINT_SET_ALL ); // A unit vector

    /* The total velocity is a weighted sum of the azimuthal velocity (whose
       speed we know) and the velocity along the magnetic field line. Although
       we don't know the magnitude of the latter, we can work it out
       geometrically (SAS triangle). As a quadratic equation, there are
       potentially two solutions, both of which we return, if present.
    */
    double a, b, c;  // The quadratic coefficients
    double vB1, vB2; // The solutions

    // (These are derived from squaring V = V_az + V_B)
    a = 1.0;
    b = az_v * (az_V.x[0] * B.x[0] +
                az_V.x[1] * B.x[1] +
                az_V.x[2] * B.x[2]);
    c = az_v*az_v - v*v;

    quad_solve_real( a, b, c, &vB1, &vB2, nsols );

    // Calculate the final velocity
    int i; // Generic counter
    if (*nsols > 0) // i.e. there is at least one solution
    {
        set_point_xyz( V1, az_v*az_V.x[0] + vB1*B.x[0],
                           az_v*az_V.x[1] + vB1*B.x[1],
                           az_v*az_V.x[2] + vB1*B.x[2], POINT_SET_ALL );
        // Normalise
        for (i = 0; i < NDEP; i++)
            V1->x[i] /= V1->r;
    }
    if (*nsols > 1) // i.e. there are two solutions
    {
        set_point_xyz( V2, az_v*az_V.x[0] + vB2*B.x[0],
                           az_v*az_V.x[1] + vB2*B.x[1],
                           az_v*az_V.x[2] + vB2*B.x[2], POINT_SET_ALL );
        // Normalise
        for (i = 0; i < NDEP; i++)
            V2->x[i] /= V2->r;
    }
}



