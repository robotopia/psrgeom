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

void obs_to_mag_frame( point *xo, pulsar *psr, psr_angle *ph, point *xm )
/* Given a point *x, get its coordinates in the magnetic frame.
 * This assumes *x has the r field defined (i.e. x.r).
 *
 * Inputs:
 *   point  *xo  : The point in the observer's frame
 *   pulsar *psr : The pulsar (needs only alpha "al" defined)
 *   psr_angle  *ph  : The rotation phase (if NULL, assumed to be zero)
 *
 * Outputs:
 *   point  *xm  : The point *x in the magnetic frame coordinate system
 */
{
    psr_angle iph, ial; // inverse ph and inverse al

    // Un-rotate the phase
    if (ph)
    {
        set_psr_angle_deg( &iph, 360.0-ph->deg );
        rotate_about_axis( xo, xm, &iph, 'z', POINT_SET_ALL );
    }
    else
    {
        copy_point( xo, xm );
    }

    // Un-rotate the inclination alpha
    set_psr_angle_deg( &ial, 360.0-psr->al.deg );
    rotate_about_axis( xm, xm, &ial, 'y', POINT_SET_ALL );
}


void mag_to_obs_frame( point *xm, pulsar *psr, psr_angle *ph, point *xo )
/* Given a point *xm, get its coordinates in the observer's frame.
 * This assumes *x has the r field defined (i.e. x.r).
 *
 * Inputs:
 *   point  *xm  : The point in the observer's frame
 *   pulsar *psr : The pulsar (needs only alpha "al" defined)
 *   psr_angle  *ph  : The rotation phase (if NULL, assumed to be zero)
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

void dipole_footpoint( pulsar *psr, double R, psr_angle *si, point *foot_pt )
/* Calculates the footpoint of a magnetic field line assuming a dipole
 * geometry.
 *
 * Inputs:
 *   pulsar *psr     : The pulsar
 *   double  R       : The maximum extent of the dipolar field line
 *   psr_angle  *si      : The magnetic azimuth of the dipolar field line
 *
 * Outputs:
 *   point  *foot_pt : The footpoint of the magnetic field line
 */
{
    // Convert magnetic azimuth to "phase"
    psr_angle ph;
    set_psr_angle_deg( &ph, si->deg + 180.0 );

    // Calculate the magnetic colatitude of the footpoint
    psr_angle th;
    set_psr_angle_sin( &th, sqrt( psr->r / R ) );

    // Get the point in magnetic coordinates
    set_point_sph( foot_pt, psr->r, &th, &ph, POINT_SET_ALL );

    // Convert from magnetic coordinates to the observer's frame
    mag_to_obs_frame( foot_pt, psr, NULL, foot_pt );
}
