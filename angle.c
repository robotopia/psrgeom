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
#include <stdio.h>
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
    // Error checking: NULLs disallowed
    if (!src || !dest)
    {
        fprintf( stderr, "error: copy_angle: src and dest cannot be NULL\n" );
        exit(EXIT_FAILURE);
    }

    // Copy values
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


void rotate_about_axis( point *in, point *out, angle *rot, char axis,
                        int flags )
/* Rotate a point about the specified "axis" by amount "rot".
 * Assumes only that the Cartesian points of *in are defined.
 */
{
    double newx, newy, newz;
    switch (axis)
    {
        case 'x':
        case 'X':
            newx = in->x[0];
            newy = in->x[1]*rot->cos - in->x[2]*rot->sin;
            newz = in->x[1]*rot->sin + in->x[2]*rot->cos;
            break;
        case 'y':
        case 'Y':
            newx =  in->x[0]*rot->cos + in->x[2]*rot->sin;
            newy =  in->x[1];
            newz = -in->x[0]*rot->sin + in->x[2]*rot->cos;
            break;
        case 'z':
        case 'Z':
            newx = in->x[0]*rot->cos - in->x[1]*rot->sin;
            newy = in->x[0]*rot->sin + in->x[1]*rot->cos;
            newz = in->x[2];
            break;
        default:
            fprintf( stderr, "error: rotate_about_axis: unrecognised axis "
                             "'%c'\n", axis );
            exit(EXIT_FAILURE);
            break;
    }

    set_point_xyz( out, newx, newy, newz, flags );
}


void min_phase_diff( angle *a1, angle *a2, angle *diff )
/* Calculate the minimum absolute phase difference between two angles.
 *
 * Inputs:
 *   angle *a1    : any arbitrary angle
 *   angle *a2    : any arbitrary angle
 *
 * Outputs:
 *   angle *diff  : the angle between a1 and a2
 */
{
    double phase_diff_deg = fabs( a1->deg - a2->deg );
    phase_diff_deg = fmod( phase_diff_deg, 360.0 ); //  0   < x < 360
    phase_diff_deg = (phase_diff_deg > 180.0 ?
                      phase_diff_deg - 360.0 :
                      phase_diff_deg);              // -180 < x < 180
    phase_diff_deg = fabs( phase_diff_deg );        //  0   < x < 180
    set_angle_deg( diff, phase_diff_deg );
}
