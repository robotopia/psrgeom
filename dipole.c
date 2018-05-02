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
        set_psr_angle_deg( &iph, -ph->deg * psr->spin );
        rotate_about_axis( xo, xm, &iph, 'z', POINT_SET_ALL );
    }
    else
    {
        copy_point( xo, xm );
    }

    // Un-rotate the inclination alpha
    set_psr_angle_deg( &ial, -psr->al.deg );
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
    rotate_about_axis( xm, xo, &psr->al, 'y', POINT_SET_ALL );

    // Rotate the phase
    if (ph)
    {
        // Adjust for pulsar's spin direction
        psr_angle iph;
        set_psr_angle_deg( &iph, ph->deg * psr->spin );
        rotate_about_axis( xo, xo, &iph, 'z', POINT_SET_ALL );
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


void beamangle_to_posangle( psr_angle *ba, psr_angle *pa )
/* This function converts a beam opening angle, Γ, to a position angle, θ,
 * as described by Eq (4) of Gangadhara & Gupta (2001).
 *
 *   Z
 *   |
 *   |
 *   |            .  o     .
 *   X      .                  .
 *   |   .                       .
 *   | .                          .
 *   |.                           .
 *   Y-------------------------------------> x
 *
 * In the figure, "o" is a point on the dipole field line, denoted with dots
 * ("."), "X" is the point on the z-axis where the tangent line of "o"
 * intersects it, Y is at the origin, and Z is a point higher up on the
 * z-axis. Then we have the definitions:
 *
 *   psr_angle *ba = Γ = ∠oXZ  (input)
 *   psr_angle *pa = θ = ∠oYZ  (output)
 *
 * This function is the inverse of posangle_to_beamangle().
 */
{
    double b = 1.5 * ba->cos / ba->sin;           /* = 3/(2 tan Γ)   */
    double theta;                                 /* = θ             */
    if (ba->deg == 0)
        theta = 0.0;
    else if (ba->deg < 180.0)                     /* if (Γ < π)      */
        theta = atan( -b + sqrt( 2.0 + b*b ) );
    else
        theta = atan( -b - sqrt( 2.0 + b*b ) ) + PI;

    set_psr_angle_rad( pa, theta );
}

void posangle_to_beamangle( psr_angle *pa, psr_angle *ba )
/* This function converts a position angle, θ, to a beam opening angle, Γ.
 *
 *   Z
 *   |
 *   |
 *   |            .  o     .
 *   X      .                  .
 *   |   .                       .
 *   | .                          .
 *   |.                           .
 *   Y-------------------------------------> x
 *
 * In the figure, "o" is a point on the dipole field line, denoted with dots
 * ("."), "X" is the point on the z-axis where the tangent line of "o"
 * intersects it, Y is at the origin, and Z is a point higher up on the
 * z-axis. Then we have the definitions:
 *
 *   psr_angle *pa = θ = ∠oYZ  (input)
 *   psr_angle *ba = Γ = ∠oXZ  (output)
 *
 * This function is the inverse of beamangle_to_posangle().
 */
{
    double gamma = atan2( 3.0*pa->cos*pa->sin,
                          3.0*pa->cos*pa->cos - 1 );

    set_psr_angle_rad( ba, gamma );
}


/* BIBLIOGRAPHY:
 *
 * Gangadhara & Gupta (2001), ApJ, 555, 31-39
 *
 */
