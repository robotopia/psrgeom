/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the struct for mathematical angles defined in
 *   psrgeom.h. It's main motivation is to avoid needlessly repeating calls to
 *   the trig functions sin and cos. It achieves this by calculating these
 *   quantities (and a few others, for convenience) just once upon creation,
 *   so that the results exist in memory for the lifetime of the struct
 *   instance.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "psrgeom.h"

angle *create_angle()
/* Allocate memory for a new instance of the angle struct, but do not populate
 * it with any value
 */
{
    return (angle *)malloc( sizeof(angle) );
}

angle *create_angle_rad( double rad )
/* Allocate memory for a new instance of the angle struct, and populate its
 * members from the supplied angle in radians
 */
{
    angle *ang = (angle *)malloc( sizeof(angle) );
    set_angle_rad( ang, rad );

    return ang;
}

angle *create_angle_deg( double deg )
/* Allocate memory for a new instance of the angle struct, and populate its
 * members from the supplied angle in degrees
 */
{
    angle *ang = (angle *)malloc( sizeof(angle) );
    set_angle_deg( ang, deg );

    return ang;
}

void destroy_angle( angle *ang )
/* Free memory assigned to an angle struct */
{
    free( ang );
}

void copy_angle( angle *src, angle *dest )
/* Copy across the contents from angle src to angle dest */
{
    dest->rad = src->rad;
    dest->deg = src->deg;
    dest->sin = src->sin;
    dest->cos = src->cos;
}

void set_angle_rad( angle *ang, double rad )
/* Fill the angle struct with calculated trigonometric values derived from the
 * supplied angle in radians
 */
{
    ang->rad   = rad;
    ang->deg   = rad*RAD2DEG;
    ang->sin   = sin(ang->rad);
    ang->cos   = cos(ang->rad);
}

void set_angle_deg( angle *ang, double deg )
/* Fill the angle struct with calculated trigonometric values derived from the
 * supplied angle in degrees
 */
{
    ang->rad   = deg*DEG2RAD;
    ang->deg   = deg;
    ang->sin   = sin(ang->rad);
    ang->cos   = cos(ang->rad);
}

void set_angle_sin( angle *ang, double Sin )
/* Fill the angle struct with calculated trigonometric values derived from the
 * supplied sin of the angle
 */
{
    ang->rad   = asin( Sin );
    ang->deg   = ang->rad*RAD2DEG;
    ang->sin   = Sin;
    ang->cos   = cos(ang->rad);
}

void set_angle_cos( angle *ang, double Cos )
/* Fill the angle struct with calculated trigonometric values derived from the
 * supplied cos of the angle
 */
{
    ang->rad   = acos( Cos );
    ang->deg   = ang->rad*RAD2DEG;
    ang->sin   = sin(ang->rad);
    ang->cos   = Cos;
}

void rotatex( angle *th, angle *ph, angle *th_out, angle *ph_out, angle *rot )
/* Rotate a direction defined by (th,ph) about the x-axis by amount "rot"
 */
{
    // Make a local copy of theta so as not to break if output pointers
    // equal input pointers
    angle t;
    copy_angle( th, &t );

    set_angle_cos( th_out, t.sin*ph->sin*rot->sin + t.cos*rot->cos );
    set_angle_rad( ph_out, atan2( t.sin*ph->sin*rot->cos - t.cos*rot->sin,
                                  t.sin*ph->cos ) );
}
