/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the struct for a pulsar, which includes
 *   information about its period, its location, its size, and its viewing
 *   geometry.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "psrgeom.h"

double light_cylinder( double P )
/* Calculates the light cylinder radius.
 *
 * Inputs:
 *    double P : the spin period of the pulsar (sec)
 *
 * Outputs:
 *    double [return] : the light cylinder radius (m)
 */
{
    return SPEED_OF_LIGHT * P / (2.0*PI);
}


void set_pulsar( pulsar *psr, psr_angle *ra, psr_angle *dec, double P, double r,
        psr_angle *al, psr_angle *ze )
{
    // Only set non-NULL angle arguments
    if (ra)
        copy_psr_angle( ra,  &psr->ra  );
    if (dec)
        copy_psr_angle( dec, &psr->dec );

    // Set the pulsar period and the light cylinder radius
    set_pulsar_period( psr, P );

    // Set the pulsar radius
    psr->r = r;

    if (al)
        copy_psr_angle( al, &psr->al );
    if (ze)
        copy_psr_angle( ze, &psr->ze );
}

pulsar *create_pulsar( psr_angle *ra, psr_angle *dec, double P, double r,
        psr_angle *al, psr_angle *ze )
{
    // Dynamically allocate memory for a pulsar struct
    pulsar *psr = (pulsar *)malloc( sizeof(pulsar) );

    // Populate the values
    set_pulsar( psr, ra, dec, P, r, al, ze );

    // Return the pointer to the new pulsar struct
    return psr;
}

void destroy_pulsar( pulsar *psr )
{
    free( psr );
}

void set_pulsar_period( pulsar *psr, double P )
{
    psr->P   = P;
    psr->rL  = light_cylinder( P );
    psr->rL2 = psr->rL * psr->rL;

    set_psr_angle_rad( &psr->Om, 2.0*PI/P );
}


void line_of_sight( pulsar *psr, psr_angle *phase, point *LoS )
/* This function finds the line-of-sight unit vector (LoS) for a given pulsar
 * geometry and rotation phase angle. The LoS is given in rotating frame
 * coordinates, where the z-axis is the rotation axis, and where the magnetic
 * axis is fixed to lie in the xz-plane.
 *
 * The point LoS is set with the POINT_SET_ALL flag, and is guaranteed to have
 * length 1, both in the sense of sqrt(x^2 + y^2 + z^2) = 1, and r = 1.
 *
 * Inputs:
 *   pulsar    *psr    : the pulsar geometry
 *   psr_angle *phase  : the rotation phase
 * Outputs:
 *   point     *LoS    : the line of sight unit vector
 */
{
    // A positive rotation of the pulsar equals an apparent reverse rotation
    // of the line of sight in the rotating frame
    psr_angle rev_phase;
    set_psr_angle_deg( &rev_phase, 360.0 - phase->deg );

    // Set the spherical coordinates for the line of sight
    double     r  = 1.0; // A unit vector
    psr_angle *th = &(psr->ze);
    psr_angle *ph = &rev_phase;

    set_point_sph( LoS, r, th, ph, POINT_SET_ALL );
}

void pol_zero( pulsar *psr, psr_angle *phase, point *pz )
/* This function finds the vector that is defined to be the 0° reference angle
 * for the PA swing. It is equivalent to the line of sight rotated 90° towards
 * the pulsar's north pole.
 *
 * The point pz is set with the POINT_SET_ALL flag, and is guaranteed to have
 * length 1, both in the sense of sqrt(x^2 + y^2 + z^2) = 1, and r = 1.
 *
 * Inputs:
 *   pulsar    *psr    : the pulsar geometry
 *   psr_angle *phase  : the rotation phase
 * Outputs:
 *   point     *pz     : the pol-zero reference vector
 */
{
    // A positive rotation of the pulsar equals an apparent reverse rotation
    // of the line of sight in the rotating frame
    psr_angle rev_phase;
    set_psr_angle_deg( &rev_phase, 360.0 - phase->deg );

    psr_angle LoS90;
    set_psr_angle_deg( &LoS90, psr->ze.deg - 90.0 );

    // Set the spherical coordinates for the line of sight
    double     r  = 1.0; // A unit vector
    psr_angle *th = &LoS90;
    psr_angle *ph = &rev_phase;

    set_point_sph( pz, r, th, ph, POINT_SET_ALL );
}


void calc_retardation( point *X, pulsar *psr, point *LoS, point *retarded_LoS )
/* This function calculates the retardation due to the difference in flight
 * time of photons starting from different locations in the magnetosphere.
 * In particular, it returns a pseudo-line-of-sight, i.e. the line of sight
 * that has been rotated around the pulsar's rotation axis by an amount
 * corresponding to the pulsar's rotation rate and the time difference between
 * a photon starting at X and a photon starting atsome reference point (here,
 * we choose the origin).
 *
 * Thus, an emission point closer to the observer than the reference point
 * would be seen at an earlier time. Since the LoS (apparently) rotates in the
 * opposite sense to the pulsar rotation (i.e. in the rotating frame), we can
 * simulate an earlier time of arrival with an artificial rotation of the line
 * of sight in the *same* sense as the pulsar.
 *
 * This function assumes that X's Cartesian coordinates are set, that its
 * radius, r, is set, and that LoS is of length unity.
 *
 * Inputs:
 *   point  *X   : the point at which the retardation is to be calculated
 *   pulsar *psr : the pulsar struct
 *   point  *LoS : the (non-retarded) line of sight
 * Outputs:
 *   point *retarded_LoS : the retarded line of sight
 */
{
    // Calculate the prejection of the position vector of X onto the LoS
    // (we assume LoS has length 1)
    double proj = X->r * (X->x[0] * LoS->x[0] +
                          X->x[1] * LoS->x[1] +
                          X->x[2] * LoS->x[2]);

    // Convert the projected length into a photon travel time
    double t = proj / SPEED_OF_LIGHT;

    // Convert the photon travel time into a pulsar rotation phase
    psr_angle dph;
    set_psr_angle_rad( &dph, t*psr->Om.rad );

    // Rotate the line of sight by that amount (around the z-axis = pulsar's
    // rotation axis)
    rotate_about_axis( LoS, retarded_LoS, &dph, 'z', POINT_SET_ALL );
}

