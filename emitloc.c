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
#include <nmead.h>
#include <time.h>
#include "psrgeom.h"

// The following struct is for conveniently passing information to the
// Nelder-Mead optimisation function
struct nm_cost_args
{
    pulsar    *psr;
    psr_angle *phase;
    int        direction;
};

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
        point B;
        calc_fields( X, psr, 0.0, &B, NULL, NULL, NULL, NULL, NULL );
        double Bdr = B.x[0] * X->ph.cos +
                     B.x[1] * X->ph.sin;
        double Bz  = B.x[2];

        // At the moment, -1 <= Bdr <= 1, so we can get a valid cost by
        // taking the absolute value. But we also want to penalise points
        // on the light cylinder whose B direction is "up"
        cost = (Bz <= 0 ? fabs(Bdr) : 2.0 - fabs(Bdr));
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
        fprintf( stderr, "warning: psr_cost_los: Could not find two solutions "
                         "for the velocity field at point X (%.2e,%.2e,%.2e) "
                         "[nsols = %d]\n", X->x[0], X->x[1], X->x[2], nsols );
        return NAN;
        //exit(EXIT_FAILURE);
    }

    // Calculate the line of sight for the given rotation phase (remember that
    // the coordinate system is in the rotating frame).

    point LoS;
    line_of_sight( psr, phase, &LoS );

    // Calculate the normalised dot product of the line of sight with the
    // velocity vector
    double LdV = LoS.x[0] * V->x[0] +
                 LoS.x[1] * V->x[1] +
                 LoS.x[2] * V->x[2];

    // At the moment, -1 <= LdV <= 1, and we only want those with an answer
    // of 1 (which means, we have to map 1 --> 0, with all the worse solutions
    // being > 0).
    cost = 1.0 - LdV;

    // Return the calculated cost
    return cost;
}

double psr_cost_total_nmead( int n, const double *xyz, void *varg )
/* A cost function for use with the Nelder-Mead algorithm, as implemented in
 * the nmead library. This cost function combines the costs associated with
 * the last open field lines (psr_cost_lofl) and the line of sight
 * (psr_cost_los).
 *
 * Arguments:
 *   n    = number of parameters to be fitted = 3 (i.e., x,y,z)
 *   xyz  = an array of three doubles, representing a point's position in
 *          Cartesian coordinates
 *   arg  = a pointer to a custom struct that contains:
 *          (1) a pointer to a pulsar struct
 *          (2) a pointer to a psr_angle struct (the phase)
 *          (3) a direction (either DIR_OUTWARD or DIR_INWARD)
 */
{
    // Force n to equal 3
    if (n != 3)
    {
        fprintf( stderr, "error: psr_cost_total_nmead: n (=%d) must "
                         "equal 3\n", n );
        exit(EXIT_FAILURE);
    }

    // Convert the parameters to a point
    point X;
    set_point_xyz( &X, xyz[0], xyz[1], xyz[2], POINT_SET_ALL );

    // "Deconstruct" the data argument "arg"
    struct nm_cost_args *arg       = (struct nm_cost_args *)varg;
    pulsar              *psr       = arg->psr;
    psr_angle           *phase     = arg->phase;
    int                  direction = arg->direction;

    // Get the two associated costs
    double cost_lofl = psr_cost_lofl( &X, psr );
    double cost_los  = psr_cost_los( &X, psr, phase, direction );

    // Combine them in the simplest way possible, and return the result
    return cost_lofl + cost_los;
}

void find_emission_point( pulsar *psr, psr_angle *phase, int direction, 
                          point *emit_pt )
/* This function finds the point which satisfies the dual contraints of
 *   1) sitting on a last open field line, and
 *   2) producing emission that is beamed along the line of sight.
 * It uses Nelder-Mead optimisation, as implemented in the nmead library.
 *
 * The initial point to feed in should be as close to the correct point as
 * possible, to minimise the time spent converging on the answer. Here,
 * however, we just randomise the initial guess, which can be improved upon
 * later if need be.
 *
 * Inputs:
 *   pulsar *psr      : a pointer to a pulsar struct
 *   psr_angle *phase : the rotation phase of interest
 *   int direction    : either DIR_OUTWARD or DIR_INWARD (for particles
 *                      flowing along or against the magnetic field,
 *                      respectively)
 * Outputs:
 *   point *emit_pt   : the emission point that the Nelder-Mead algorithm
 *                      converged upon
 */
{
    // Assume that the random number generator has been already seeded

    // Set up the arguments for the call to the Nelder-Mead function
    int n = 3;     // Three parameters to fit (x,y,z)

    double p0[n];  // The initial guess
    int i;         // for looping over 0..(n-1)
    do
    {
        for (i = 0; i < n; i++)
            p0[i] = RANDU;
    } while ((p0[0]*p0[0] + p0[1]*p0[1]) >= 1.0); /* Make sure point is within
                                                     the light cylinder */
    // Scale up to physical units
    for (i = 0; i < n; i++)
        p0[i] *= psr->rL;

    // Use default Nelder-Mead algorithm parameters
    nm_optimset optimset;

    optimset.tolx     = NM_TOL_X;
    optimset.tolf     = NM_TOL_F;
    optimset.max_iter = NM_MAX_ITER;
    optimset.max_eval = NM_MAX_EVAL;

    // Set up a place to store the result
    nm_point solution;
    double solvals[n];
    solution.x = solvals;

    // Set up the data args that get passed to the cost function
    struct nm_cost_args arg;
    arg.psr       = psr;
    arg.phase     = phase;
    arg.direction = direction;

    // Call the Nelder-Mead function and find the point
    nelder_mead( p0, n, optimset, &solution, psr_cost_total_nmead, &arg );

    // Pass the solution out of this function
    set_point_xyz( emit_pt, solution.x[0],
                            solution.x[1],
                            solution.x[2],
                            POINT_SET_ALL );
}
