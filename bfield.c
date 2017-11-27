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
    angle phase;
    set_angle_rad( &phase, (r - psr->r)/psr->rL );

    // Paul Arendt's equations
    double a  = psr->al.cos / r5;
    double b  = psr->al.sin * phase.cos / r3 / psr->rL2;
    double c  = psr->al.sin * (phase.cos/r5 + phase.sin/r4/psr->rL);
    double d  = psr->al.sin * phase.sin / r3 / psr->rL2;
    double e  = d*psr->rL2/rr - b*psr->rL/r;

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
    return (B.x[0] * x->ph.cos +
            B.x[1] * x->ph.sin);
}

