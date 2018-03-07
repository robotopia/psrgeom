/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2018
 *
 * Description:
 *   This source file implements the functions required to find the emission
 *   locations for a given rotation phase. It uses a Nelder-Mead algorithm to
 *   converge on the desired point. The cost function incorporates both how
 *   close a point's field line is to the last open field line, and how
 *   close the emission direction is to the line of sight.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psrgeom.h"

double psr_cost_lofl( point *X, pulsar *psr )
/* This calculates the "cost" associated with how close X's magnetic field
 * line is to a last open field line. Because the field lines behave
 * erratically outside the light cylinder, it is impractical to always find
 * the line's most distant point (i.e. perp. dist. from the z-axis) and
 * use that in the calculation. Therefore, we break the problem up into two
 * cases:
 *   1) the field line extends beyond the light cylinder,
 *   2) the field line is closed (within the light cylinder)
 * In the first case, we find the point at which the field line penetrates the
 * light cylinder and use the "B dot r_{xy}" quantity as a proxy for the cost.
 * In the second case, we use fractional distance of the most extreme point to
 * the light cylinder as a proxy for the cost.
 *
 * Inputs:
 *   point *X    : The point whose cost we want to evaluate
 *   pulsar *psr : The pulsar whose magnetic field we want to use for the
 *                 calculation
 * Return value:
 *   cost        : The cost of the point X. It is guaranteed to be a number
 *                 between 0 and 1, where 0 means that X is on a last open
 *                 field line.
 */
{
    // Set up the return variable for the cost
    double cost;

    // Set up the parameters for the farpoint() function
    double tmult     = 0.01;  // Hard code the step size
    FILE  *no_output = NULL;  // Don't write out the intermediate points
    int    rL_norm   = 0;     // (Irrelevant, since not writing anything out)
    double rL_lim    = 1.0;   // i.e. Stop at the light cylinder
    point  far_pt;            // Where to store the final point

    // Find the most extreme point (or where it crosses the light cylinder)
    int stop_type = farpoint( X, psr, tmult, no_output, rL_norm, rL_lim,
                              &far_pt );

    // Consider the two cases separately:
    //   1) the field line extends beyond the light cylinder
    //   2) the extreme point is within the light cylinder

    if (stop_type == STOP_EXCEED)
    {
        // In this case, far_pt is the point where the magnetic field line
        // penetrates the light cylinder
        double Bdr = Bdotrxy( &far_pt, psr );

        // At the moment, -1 <= Bdr <= 1, so we can get a valid cost by
        // squaring it.
        cost = Bdr * Bdr;
    }
    else /* (stop_type == STOP_FOUND) */
    {
        // In this case, far_pt is the point on the field line with the
        // largest perpendicular distance to the rotation (z) axis
        double rhosq = far_pt.rhosq;

        // The fractional distance (squared) will always fall neatly into the
        // range: 0 <= (rho/rL)^2 <= 1, so we just need to invert it so that "0"
        // means the last open field line
        cost = 1.0 - (rhosq / psr->rL2);
    }

    // Return the calculated cost
    return cost;
}


double psr_cost_los( point *X, pulsar *psr, psr_angle *phase, int direction )
/* This calculates the "cost" associated with how close emission from point X
 * is from the line of sight, given a rotation phase. Each point has
 * associated with it two possible velocity vectors (which determine the
 * direction of the beamed emission), depending on whether the particle
 * is flowing "outward" or "inward" (i.e. in the same or in the opposite
 * direction as the magnetic field. Which vector field is used is set by the
 * "direction" parameter.
 *
 * Inputs:
 *   point *X         : The point whose cost we want to evaluate
 *   pulsar *psr      : The pulsar whose magnetic field we want to use for the
 *                      calculation
 *   psr_angle *phase : The rotation phase of the pulsar
 *   int direction    : Must be either DIR_OUTWARD or DIR_INWARD
 * Return value:
 *   cost        : The cost of the point X. It is guaranteed to be a number
 *                 between 0 and 1, where 0 means that X is on a last open
 *                 field line.
 */
{
    // First, check that the direction is either DIR_OUTWARD or DIR_INWARD
    if (direction != DIR_OUTWARD && direction != DIR_INWARD)
    {
        fprintf( stderr, "error: psr_cost_los: direction must be either "
                         "DIR_OUTWARD or DIR_INWARD\n" );
        exit(EXIT_FAILURE);
    }

    // Set up the return variable for the cost
    double cost;

    // Set up the points for the velocity vector field
    // (V1 is for DIR_OUTWARD, V2 is for DIR_INWARD)
    point V1, V2;
    point *V = (direction == DIR_OUTWARD ? &V1 : &V2);

    // Use the speed of light as the particle velocity, keeping it in line
    // with the assumption of a highly relativistic plasma
    double v = SPEED_OF_LIGHT;

    // Calculate the vector field at X
    int nsols;
    calc_fields( X, psr, v, NULL, &V1, &V2, NULL, NULL, &nsols );

    // Check to make sure we actually have two solutions. If not, exit with an
    // error
    if (nsols != 2)
    {
        fprintf( stderr, "error: psr_cost_los: Could not find two solutions "
                         "for the velocity field at point X (%.2e,%.2e,%.2e) "
                         "[nsols = %d]\n", X->x[0], X->x[1], X->x[2], nsols );
        exit(EXIT_FAILURE);
    }

    // Calculate the line of sight for the given rotation phase (remember that
    // the coordinate system is in the rotating frame, so that setting the
    // phase corresponds to rotating the line of sight about the rotation
    // axis). NB: The line of sight (apparently) rotates about the z-axis in
    // the opposite direction to the pulsar rotation.
    point LoS;

    psr_angle rev_phase;
    set_psr_angle_deg( &rev_phase, 360.0 - phase->deg );

    // Set the spherical coordinates for the line of sight
    double     r  = 1.0; // A unit vector
    psr_angle *th = &(psr->ze);
    psr_angle *ph = &rev_phase;

    set_point_sph( &LoS, r, th, ph, POINT_SET_ALL );

    // Calculate the normalised dot product of the line of sight with the
    // velocity vector
    double LdV = norm_dot( &LoS, V );

    // At the moment, -1 <= LdV <= 1, so squaring it will give us something in
    // the correct range. But we need 0 to represent when the two vectors are
    // parallel, so after squaring, invert.
    cost = 1.0 - LdV*LdV;

    // Return the calculated cost
    return cost;
}
