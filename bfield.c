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

void Bfield( point *x, pulsar *psr, point *B )
/* This returns the normalised vector B (cast as a "point"), i.e. the magnetic
 * field direction at the point *x around pulsar *psr, revolving about the
 * z-axis.
 */
{
    double xx    = x->x[0] * x->x[0];
    double yy    = x->x[1] * x->x[1];
    double zz    = x->x[2] * x->x[2];

    double xy    = x->x[0] * x->x[1];
    double xz    = x->x[0] * x->x[2];
    double yz    = x->x[1] * x->x[2];

    double rr    = xx + yy + zz;
    double r     = sqrt( rr );
    double r3    = r*rr;
    double r4    = r*r3;
    double r5    = r*r4;

    // Taking into account finite light travel time...
    // Set up so that pole is where we think it is
    psr_angle phase;
    set_psr_angle_rad( &phase, (r - psr->r)/psr->rL );

    // Paul Arendt's equations
    double a  =  1.0 / r5;
    double b  =  psr->al.sin * phase.cos / r3 / psr->rL2;
    double c  =  psr->al.sin * ( phase.cos/r5 + phase.sin/r4/psr->rL);
    double d  = -psr->al.sin * phase.sin / r3 / psr->rL2;
    double e  =  psr->al.sin * (-phase.sin/r5 + phase.cos/r4/psr->rL);

    set_point_xyz( B,
                   3.0*a*xz + b*(yy + zz) + c*(3.0*xx - rr) - d*xy + 3.0*e*xy,
                   3.0*a*yz - b*xy + 3.0*c*xy + d*(xx + zz) + e*(3.0*yy - rr),
                   a*(3.0*zz - rr) - b*xz + 3.0*c*xz - d*yz + 3.0*e*yz,
                   POINT_SET_XYZ | POINT_SET_R );

    // Normalise
    int i;
    for (i = 0; i < NDEP; i++)
        B->x[i] /= B->r;
}



double Bdotrxy( point *x, pulsar *psr )
/* Returns the dot product of B, the normalised magnetic field, and r_{xy},
 * the vector pointing cylindrically outward, at point *x.
 * This quantity is useful because the last open magnetic field lines have
 * this quantity equal to zero at the light cylinder.
 */
{
    point B;
    Bfield( x, psr, &B );
    double retval = B.x[0] * x->ph.cos +
                    B.x[1] * x->ph.sin;
    return retval;
}



void footpoint( point *start_pt, pulsar *psr, double tmult, int direction,
                FILE *write_xyz, point *foot_pt )
/* This uses Runge-Kutta (RK4) to follow the magnetic field line outward from
 * an initial starting point "start_pt"
 *
 * Inputs:
 *   start_pt  : the initial starting point
 *   psr       : the pulsar (includes geometric information)
 *   tmult     : the rate at which to progress in RK4 (fraction of psr radius)
 *   direction : progress inwards or outwards along field lines
 *               {DIR_INWARD, DIR_OUTWARD}
 *   write_xyz : file handle where to write x,y,z values
 *               (if NULL, no output is written)
 * Outputs:
 *   foot_pt  : the final point reached during RK algorithm
 */
{
    // Error checking: foot_pt must be a valid pointer to a point struct
    if (!foot_pt)
    {
        fprintf( stderr, "error: footpoint: foot_pt cannot be NULL\n" );
        exit(EXIT_FAILURE);
    }

    point x, old_x;

    double tstep;

    copy_point( start_pt, &x );
    tstep = tmult * x.r;

    // Trace this line outwards with a 4 stage Runge-Kutta

    int temp_drctn = direction; // In case we've gone too far
    double precision = 1.0e-14;

    while (1) // Indefinite loop until surface is reached
    {
        // Keep track of the previous point
        copy_point( &x, &old_x );

        // Write out the current xyz position, if requested,
        // but only if we're still moving "forward"
        if (write_xyz && (temp_drctn == direction))
            fprintf( write_xyz, "%.14e %.14e %.14e\n", x.x[0], x.x[1], x.x[2] );

        // Take a single RK4 step along the magnetic field
        Bstep( &old_x, psr, tstep, temp_drctn, &x );

        // Recalculate the various distances to the new point
        set_point_xyz( &x, x.x[0], x.x[1], x.x[2],
                POINT_SET_SPH | POINT_SET_RHOSQ );

        /* Figure out whether to stop or not */

        // Error checking: this algorithm should have stopped long before
        // tstep reaches underflow
        if ((tstep/2.0 <= 0.0) || (tstep/2.0 >= tstep))
        {
            fprintf( stderr, "error: Bline: tstep underflow\n" );
            exit(EXIT_FAILURE);
        }

        // Just check to see if we've happened to land exactly on the surface
        if (x.r == psr->r)
            break;

        // Otherwise, only do anything special if we've crossed the surface
        if ((x.r - psr->r) / (old_x.r - psr->r) < 0.0) // then x and old_x are
                                                       // on opposite sides of
                                                       // the surface
        {
            if ((fabs(x.r - old_x.r) <= precision*psr->r))
                break;

            if (temp_drctn == DIR_OUTWARD)
                temp_drctn = DIR_INWARD;
            else if (temp_drctn == DIR_INWARD)
                temp_drctn = DIR_OUTWARD;
            else
            {
                fprintf( stderr, "error: footpoint: unknown direction\n" );
                exit(EXIT_FAILURE);
            }
            tstep /= 2.0;
        }

        // Adjust tstep proportionally to how far away from the pulsar we are
        tstep *= x.r / old_x.r;
    }

    // Make the final point available to the caller
    copy_point( &x, foot_pt );
}


void Bstep( point *x1, pulsar *psr, double tstep, int direction, point *x2 )
/* This follows a magnetic field from a given starting point according to one
 * step of the RK4 algorithm.
 *
 * Inputs:
 *   x1       : the starting point
 *   psr      : the pulsar (includes geometric information)
 *   tstep    : the size of the RK4 step to take
 *   direction: whether to follow the line DIR_INWARD or DIR_OUTWARD
 *
 * Outputs:
 *   x2       : the ending point
 */
{
    point slop1, slop2, slop3, slope;
    point xp1, xp2;

    double sgn = (direction ? 1.0 : -1.0);
    int i; // Generic loop counter

    // First stage
    Bfield( x1, psr, &slop1 );
    for (i = 0; i < NDEP; i++)
        xp1.x[i] = x1->x[i] + 0.5*sgn*slop1.x[i]*tstep;

    // Second stage
    Bfield( &xp1, psr, &slop2 );
    for (i = 0; i < NDEP; i++)
        xp2.x[i] = x1->x[i] + 0.5*sgn*slop2.x[i]*tstep;

    // Third stage
    Bfield( &xp2, psr, &slop3 );
    for (i = 0; i < NDEP; i++)
        xp1.x[i] = x1->x[i] + sgn*slop3.x[i]*tstep;

    // Last stage
    Bfield( &xp1, psr, &slope );
    for (i = 0; i < NDEP; i++)
    {
        x2->x[i] = x1->x[i] +
            tstep*sgn*(    slope.x[i] +     slop1.x[i] +
                       2.0*slop2.x[i] + 2.0*slop3.x[i]) / 6.0;
    }
}
