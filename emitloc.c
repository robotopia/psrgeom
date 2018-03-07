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
