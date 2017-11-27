/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the struct for geometric points defined in
 *   psrgeom.h.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "psrgeom.h"

void set_point_xyz( point *p, double x, double y, double z, int flags )
/* Populate a point structure given x,y,z coordinates
 */
{
    // Set cartesian coordinates
    if (flags & POINT_SET_X)
        p->x[0] = x;
    if (flags & POINT_SET_Y)
        p->x[1] = y;
    if (flags & POINT_SET_Z)
        p->x[2] = z;

    // Calculate cylindrical distance rho squared
    if (flags & POINT_SET_RHOSQ)
        p->rhosq = x*x + y*y;

    // Calculate spherical coordinates
    if ((flags & POINT_SET_R) || (flags & POINT_SET_TH))
        p->r = sqrt( x*x + y*y + z*z );
    if (flags & POINT_SET_TH)
        set_angle_cos( &p->th, z/(p->r) );
    if (flags & POINT_SET_PH)
        set_angle_rad( &p->ph, atan2(y,x) );
}

void set_point_sph( point *p, double r, angle *th, angle *ph, int flags )
/* Populate a point structure given spherical coordinates
 */
{
    // Set spherical coordinates
    if (flags & POINT_SET_R)
        p->r = r;
    if (flags & POINT_SET_TH)
        copy_angle( th, &p->th );
    if (flags & POINT_SET_PH)
        copy_angle( ph, &p->ph );

    // Calculate cartesian coordinates
    if ((flags & POINT_SET_X) || (flags & POINT_SET_RHOSQ))
        p->x[0] = r * th->sin * ph->cos;
    if ((flags & POINT_SET_Y) || (flags & POINT_SET_RHOSQ))
        p->x[1] = r * th->sin * ph->sin;
    if (flags & POINT_SET_Z)
        p->x[2] = r * th->cos;

    // Calculate cylindrical distance rho squared
    if (flags & POINT_SET_RHOSQ)
        p->rhosq = p->x[0]*p->x[0] + p->x[1]*p->x[1];
}

void copy_point( point *src, point *dest )
/* Copy a point from one struct to another
 */
{
    dest->x[0]  = src->x[0];
    dest->x[1]  = src->x[1];
    dest->x[2]  = src->x[2];
    dest->r     = src->r;
    copy_angle( &src->th, &dest->th );
    copy_angle( &src->ph, &dest->ph );
    dest->rhosq = src->rhosq;
}


