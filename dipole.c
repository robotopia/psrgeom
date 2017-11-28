/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the functions for calculating quantities
 *   associated with a dipole magnetic field.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "psrgeom.h"

void obs_to_mag_frame( point *xo, pulsar *psr, angle *ph, point *xm )
/* Given a point *x, get its coordinates in the magnetic frame.
 * This assumes *x has the r field defined (i.e. x.r).
 *
 * Inputs:
 *   point  *xo  : The point in the observer's frame
 *   pulsar *psr : The pulsar (needs only alpha "al" defined)
 *   angle  *ph  : The rotation phase (if NULL, assumed to be zero)
 *
 * Outputs:
 *   point  *xm  : The point *x in the magnetic frame coordinate system
 */
{
    angle iph, ial; // inverse ph and inverse al

    // Un-rotate the phase
    if (ph)
    {
        set_angle_deg( &iph, 360.0-ph->deg );
        rotate_about_axis( xo, xm, &iph, 'z', POINT_SET_ALL );
    }
    else
    {
        copy_point( xo, xm );
    }

    // Un-rotate the inclination alpha
    set_angle_deg( &ial, 360.0-psr->al.deg );
    rotate_about_axis( xm, xm, &ial, 'y', POINT_SET_ALL );
}


void mag_to_obs_frame( point *xm, pulsar *psr, angle *ph, point *xo )
/* Given a point *xm, get its coordinates in the observer's frame.
 * This assumes *x has the r field defined (i.e. x.r).
 *
 * Inputs:
 *   point  *xm  : The point in the observer's frame
 *   pulsar *psr : The pulsar (needs only alpha "al" defined)
 *   angle  *ph  : The rotation phase (if NULL, assumed to be zero)
 *
 * Outputs:
 *   point  *xo  : The point *x in the magnetic frame coordinate system
 */
{
    // Rotate the inclination alpha
    rotate_about_axis( xo, xm, &psr->al, 'y', POINT_SET_ALL );

    // Rotate the phase
    if (ph)
    {
        rotate_about_axis( xm, xm, ph, 'z', POINT_SET_ALL );
    }
}


double calc_dipole_R( point *xm )
/* Calculates the maximum extent, R, of a dipole magnetic field line that
 * passes through the given point, *xm, assumed to be in magnetic coordinates.
 * It is assumed that the given point's "r" value is defined.
 */
{
    return xm->r / (xm->th.sin * xm->th.sin);
}
