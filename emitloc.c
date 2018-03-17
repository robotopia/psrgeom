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
#include <newuoa.h>
#include <time.h>
#include "psrgeom.h"

// The following struct is for conveniently passing information to the
// Nelder-Mead optimisation function
struct cost_args
{
    pulsar    *psr;
    psr_angle *phase;
    int        direction;
    FILE      *f;
    double     r;  // <-- Only used by some cost functions
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
    double rL_lim    = 1.1;   // i.e. Stop at the light cylinder
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
        calc_fields( &far_pt, psr, 0.0, &B, NULL, NULL, NULL, NULL, NULL );
        double Bdr = B.x[0] * far_pt.ph.cos +
                     B.x[1] * far_pt.ph.sin;

        // At the moment, -1 <= Bdr <= 1, so we can get a valid cost by
        // taking the absolute value.
        cost = fabs(Bdr) + rL_lim - 1.0;
        //cost = 1.0; // Temporary, while debugging cost function
    }
    else /* (stop_type == STOP_FOUND) */
    {
        // In this case, far_pt is the point on the field line with the
        // largest perpendicular distance to the rotation (z) axis

        // The fractional distance (squared) will always fall neatly into the
        // range: 0 <= (rho/rL)^2 <= rL_lim^2, so we just need to invert it so
        // that "0" means the last open field line
        cost = fabs(1.0 - sqrt( far_pt.rhosq / psr->rL2 ));
    }

    // Return the calculated cost
    return cost;
}


double psr_cost_los( point *X, pulsar *psr, psr_angle *phase, int direction,
                     int retardation )
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
 *   int retardation  : Whether to adjust for retardation effects or not
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

    // If retardation effects are requested, offset the line of sight by the
    // appropriate amount
    if (retardation)
        calc_retardation( X, psr, &LoS, &LoS );

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

double psr_cost_los_at_r_optim( long int n, const double *thph, void *varg )
/* A cost function for use with an optimisation algorithm. This cost function
 * retrieves the cost associated with the line of sight (psr_cost_los).
 * Moreover, it assumes that the search is taking place on a constant sphere
 * centred on the origin.
 *
 * Arguments:
 *   n    = number of parameters to be fitted = 3 (i.e., x,y,z)
 *   thph = an array of two doubles, representing a point's position in
 *          spherical coordinates (r, in varg, is constant)
 *   varg = a pointer to a custom struct that contains:
 *          (1) a pointer to a pulsar struct
 *          (2) a pointer to a psr_angle struct (the phase)
 *          (3) a direction (either DIR_OUTWARD or DIR_INWARD)
 *          (4) a file handle to print out intermediate steps (can be NULL)
 *          (5) The radius, r
 */
{
    // Force n to equal 2
    if (n != 2)
    {
        fprintf( stderr, "error: psr_cost_los_at_r_optim: n (=%ld) must "
                         "equal 2\n", n );
        exit(EXIT_FAILURE);
    }

    // "Deconstruct" the data argument "arg"
    struct cost_args *arg       = (struct cost_args *)varg;
    pulsar           *psr       = arg->psr;
    psr_angle        *phase     = arg->phase;
    int               direction = arg->direction;
    FILE             *f         = arg->f;
    double            r         = arg->r;

    // Convert the parameters to a point
    point X;
    psr_angle th, ph;
    set_psr_angle_rad( &th, thph[0] );
    set_psr_angle_rad( &ph, thph[1] );
    set_point_sph( &X, r, &th, &ph, POINT_SET_ALL );

    // Get the "LoS" cost
    int retardation = 1; // Turn on retardation effects
    double cost_los = psr_cost_los( &X, psr, phase, direction, retardation );

    if (f != NULL)
        fprintf( f, "%.15e %.15e %.15e %.15e\n",
                    X.x[0], X.x[1], X.x[2], cost_los );

    // Return the result
    return cost_los;
}


double psr_cost_total( long int n, const double *xyz, void *varg )
/* A cost function for use with an optimisation algorithm. This cost function
 * combines the costs associated with the last open field lines
 * (psr_cost_lofl) and the line of sight (psr_cost_los).
 *
 * Arguments:
 *   n    = number of parameters to be fitted = 3 (i.e., x,y,z)
 *   xyz  = an array of three doubles, representing a point's position in
 *          Cartesian coordinates
 *   varg = a pointer to a custom struct that contains:
 *          (1) a pointer to a pulsar struct
 *          (2) a pointer to a psr_angle struct (the phase)
 *          (3) a direction (either DIR_OUTWARD or DIR_INWARD)
 *          (4) a file handle to print out intermediate steps (can be NULL)
 */
{
    // Force n to equal 3
    if (n != 3)
    {
        fprintf( stderr, "error: psr_cost_total: n (=%ld) must "
                         "equal 3\n", n );
        exit(EXIT_FAILURE);
    }

    // Convert the parameters to a point
    point X;
    set_point_xyz( &X, xyz[0], xyz[1], xyz[2], POINT_SET_ALL );

    // "Deconstruct" the data argument "arg"
    struct cost_args *arg       = (struct cost_args *)varg;
    pulsar              *psr       = arg->psr;
    psr_angle           *phase     = arg->phase;
    int                  direction = arg->direction;
    FILE                *f         = arg->f;

    // Get the two associated costs
    int retardation = 1; // Turn on retardation effects
    double cost_lofl = psr_cost_lofl( &X, psr );
    double cost_los  = psr_cost_los( &X, psr, phase, direction, retardation );

    if (f != NULL)
        fprintf( f, "%.15e %.15e %.15e %.15e %.15e\n",
                    xyz[0], xyz[1], xyz[2], cost_lofl, cost_los );

    // Combine the two costs, and return the result
    //double coeff = 1.0;
    return 1e10 * (cost_lofl*cost_lofl*cost_lofl + cost_los);

    /* The coefficient above seems to be a magic number; without it (or when
     * it has a different value), the Nelder-Mead algorithm function fails to
     * converge at the correct point. For the geometries:
     *
     *   α = 45°, ζ = 40°, φ = -45°, P = 1.0 sec,
     *   α = 30°, ζ = 40°, φ =   0°, P = 0.5 sec,
     *
     * a value around 0.04 or 0.05 performs best. However, for
     *
     *   α = 10°, ζ = 20°, φ = 110°, P = 0.01 sec,
     *
     * a value of 0.02 works, while 0.04 doesn't.
     */
}


void find_approx_emission_point( pulsar *psr, psr_angle *phase, int direction,
                                 point *emit_pt )
/* This function finds the emission point around a "simplified" pulsar, which
 * means a pulsar that
 *   1) is not rotating,
 *   2) has a perfectly dipolar magnetic field,
 *   3) emits along field lines whose maximum extent is the light cylinder
 *      radius (NB this is not generally the same as the last open field
 *      lines).
 * The main purpose of this function is to seed the input of the
 * find_emission_point() function, which solves for the full Deutsch (1955)
 * magnetic field, for a rotating pulsar.
 *
 * Inputs:
 *   pulsar *psr      : a pointer to a pulsar struct
 *   psr_angle *phase : the rotation phase of interest
 *   int direction    : whether the emitting particles are flowing along
 *                      (DIR_OUTWARD) or against (DIR_INWARD) the magnetic
 *                      field
 * Outputs:
 *   point *emit_pt   : the found emission point
 */
{
    // Make sure that direction is either DIR_OUTWARD or DIR_INWARD
    if (direction != DIR_OUTWARD && direction != DIR_INWARD)
    {
        fprintf( stderr, "error: find_approx_emission_point: "
                         "direction must be either DIR_OUTWARD or "
                         "DIR_INWARD\n" );
        exit(EXIT_FAILURE);
    }

    // Get the line of sight
    point LoS;
    line_of_sight( psr, phase, &LoS );

    // Get the magnetic axis (as a unit vector aka point)
    point mag;
    psr_angle angle_zero = ANGLE_ZERO;
    set_point_sph( &mag, 1.0, &(psr->al), &angle_zero, POINT_SET_ALL );

    // Get Γ, the angle that the line of sight makes with the magnetic axis
    psr_angle gamma;
    set_psr_angle_cos( &gamma, LoS.x[0] * mag.x[0] +
                               LoS.x[1] * mag.x[1] +
                               LoS.x[2] * mag.x[2] );

    // Convert the opening beam angle, Γ, to a position angle, θ
    psr_angle theta;
    beamangle_to_posangle( &gamma, &theta );

    // Calculate the emission point's distance from the origin, r
    double r = psr->rL * theta.sin * theta.sin;

    // Calculate the emission point's "phase" (in the magnetic frame, σ) by
    // finding the "phase" of the line of sight (because the magnetic axis,
    // line of sight, and emission point are co-planar).

    point LoS_mag;        /* The LoS in magnetic frame coordinates */
    obs_to_mag_frame( &LoS, psr, NULL, &LoS_mag );
    psr_angle *sigma = &(LoS_mag.ph);

    // Now I have the spherical coordinates (r,θ,φ) in the magnetic frame,
    // so set the emission point in the mag. frame and then convert back
    // to the observer frame.
    point emit_pt_mag;
    set_point_sph( &emit_pt_mag, r, &theta, sigma, POINT_SET_ALL );

    mag_to_obs_frame( &emit_pt_mag, psr, NULL, emit_pt );

    // In the case of DIR_INWARD, the emission point will be exactly opposite
    // the one just calculated (i.e. reflected in the origin)
    if (direction == DIR_INWARD)
    {
        set_point_xyz( emit_pt, -emit_pt->x[0],
                                -emit_pt->x[1],
                                -emit_pt->x[2],
                                POINT_SET_ALL );
    }

}


void find_emission_point_nmead( pulsar *psr, psr_angle *phase, int direction, 
                                point *emit_pt, FILE *f )
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

    point init_guess;
    find_approx_emission_point( psr, phase, direction, &init_guess );

    for (i = 0; i < n; i++)
        p0[i] = init_guess.x[i];

    // Use default Nelder-Mead algorithm parameters
    nm_optimset optimset;

    optimset.tolx     = NM_TOL_X * 1e-16;
    optimset.tolf     = NM_TOL_F * 1e-16;
    optimset.max_iter = NM_MAX_ITER;
    optimset.max_eval = NM_MAX_EVAL;

    // Set up a place to store the result
    nm_point solution;
    double solvals[n];
    solution.x = solvals;

    // Set up the data args that get passed to the cost function
    struct cost_args arg;
    arg.psr       = psr;
    arg.phase     = phase;
    arg.direction = direction;
    arg.f         = f;

    // Call the Nelder-Mead function and find the point
    nelder_mead( p0, n, optimset, &solution, psr_cost_total, &arg );

    // Pass the solution out of this function
    set_point_xyz( emit_pt, solution.x[0],
                            solution.x[1],
                            solution.x[2],
                            POINT_SET_ALL );
}


void find_emission_point_newuoa( pulsar *psr, psr_angle *phase, int direction, 
                                 point *emit_pt, FILE *f )
/* This function finds the point which satisfies the dual contraints of
 *   1) sitting on a last open field line, and
 *   2) producing emission that is beamed along the line of sight.
 * It uses Powell's NEWUOA optimisation, as implemented in the newuoa library.
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
    // Set up the arguments for the call to the NEWUOA function
    int n   = 3;     // Three parameters to fit (x,y,z)
    int npt = n+2;   // The number of interpolation conditions
                     // Must be in range [n+2,(n+1)(n+2)/2]
    newuoa_objfun *objfun = psr_cost_total; // The function to be evaluated

    // Set up the initial guess
    double x[n]; // <-- This is also where the answer will end up
    int i;

    point init_guess;
    find_approx_emission_point( psr, phase, direction, &init_guess );

    for (i = 0; i < n; i++)
        x[i] = init_guess.x[i];

    // Create memory space for working calculations
    int workspace_size = (npt+13)*(npt+n)+3*n*(n+3)/2; // Mandated by the docs
    double w[workspace_size];

    // NEWUOA expects two variables to be set: RHOBEG and RHOEND:
    // The documentation says:
    //
    // "RHOBEG and RHOEND must be set to the initial and final values of a
    // trust region radius, so both must be positive with RHOEND<=RHOBEG.
    // Typically RHOBEG should be about one tenth of the greatest expected
    // change to a // variable, and RHOEND should indicate the accuracy that
    // is required in the final values of the variables."
    //
    // I'll set RHOBEG to 10% of the radial distance of the initial point,
    // and RHOEND to 1 metre.
    double rhobeg = 0.1 * init_guess.r;
    double rhoend = 1e-12;

    // I've got my own printf'ing happening, so turn NEWUOA's internal
    // reporting off
    int iprint = 0;

    // Set the maximum number of calls to the cost function
    int maxfun = 1000;

    // Set up the data args that get passed to the cost function
    struct cost_args data;
    data.psr       = psr;
    data.phase     = phase;
    data.direction = direction;
    data.f         = f;

    // Call the NEWUOA function and find the minimum point
    int status = newuoa( n, npt, objfun, &data, x, rhobeg, rhoend,
                         iprint, maxfun, w );

    if (status != NEWUOA_SUCCESS)
    {
        fprintf( f, "warning: find_emission_point_newuoa: "
                    "%s\n", newuoa_reason( status ) );
    }

    // Pass the solution out of this function
    set_point_xyz( emit_pt, x[0], x[1], x[2], POINT_SET_ALL );
}


/* BIBLOGRAPHY
 *
 * Deutsch, A. (1955), Annales d'Astrophysique, 18, 1-10
 */


void psr_cost_deriv( point *X, pulsar *psr, psr_angle *phase, int direction,
                     double dx, point *grad )
/* This function finds the grad (∇) of the cost function (psr_cost_total)
 * numerically at the specified point X. The caller chooses how closely to
 * sample the nearby points (dx).
 */
{
    struct cost_args arg;
    arg.psr       = psr;
    arg.phase     = phase;
    arg.direction = direction;
    arg.f         = NULL;

    int n = 3;
    double xyz[n];
    double xyz_dx[n];
    double xyz_dy[n];
    double xyz_dz[n];

    int i;
    for (i = 0; i < n; i++)
        xyz[i] = xyz_dx[i] = xyz_dy[i] = xyz_dz[i] = X->x[i];

    xyz_dx[0] += dx;
    xyz_dy[1] += dx;
    xyz_dz[2] += dx;

    double c0, cx, cy, cz;

    // Get the function at the points X, X+dx, X+dy, X+dz
    c0 = psr_cost_total( n, xyz   , &arg );
    cx = psr_cost_total( n, xyz_dx, &arg );
    cy = psr_cost_total( n, xyz_dy, &arg );
    cz = psr_cost_total( n, xyz_dz, &arg );

    set_point_xyz( grad, (cx-c0)/dx, (cy-c0)/dx, (cz-c0)/dx, POINT_SET_ALL );
}


void find_LoS_at_r( point *init_pt, pulsar *psr, psr_angle *phase,
                    int direction, point *end_pt, FILE *f )
/* This function finds the point nearest* the initial point which satisfies
 * the criterion that the velocity field there is parallel to the line of
 * sight.
 *
 * Inputs:
 *   point *init_pt   : a point to start at
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
    // Set up the arguments for the call to the NEWUOA function
    int n   = 2;     // Two parameters to fit (θ,φ)
    int npt = n+2;   // The number of interpolation conditions
                     // Must be in range [n+2,(n+1)(n+2)/2]
    newuoa_objfun *objfun = psr_cost_los_at_r_optim;

    // Set up the initial guess
    double x[n]; // <-- This is also where the answer will end up

    x[0] = init_pt->th.rad;
    x[1] = init_pt->ph.rad;

    // Create memory space for working calculations
    int workspace_size = (npt+13)*(npt+n)+3*n*(n+3)/2; // Mandated by the docs
    double w[workspace_size];

    // NEWUOA expects two variables to be set: RHOBEG and RHOEND:
    // The documentation says:
    //
    // "RHOBEG and RHOEND must be set to the initial and final values of a
    // trust region radius, so both must be positive with RHOEND<=RHOBEG.
    // Typically RHOBEG should be about one tenth of the greatest expected
    // change to a // variable, and RHOEND should indicate the accuracy that
    // is required in the final values of the variables."
    //
    // Since I'm hoping that the initial guess point is reasonably close to
    // the final point on the sphere, and since my parameters are given in
    // radians, I'll set my RHOBEG to 0.1 (rad) and RHOEND to 1e-16 (rad).
    double rhobeg = 0.1;
    double rhoend = 1e-16;

    // I've got my own printf'ing happening, so turn NEWUOA's internal
    // reporting off
    int iprint = 0;

    // Set the maximum number of calls to the cost function
    int maxfun = 1000;

    // Set up the data args that get passed to the cost function
    struct cost_args data;
    data.psr       = psr;
    data.phase     = phase;
    data.direction = direction;
    data.f         = f;
    data.r         = init_pt->r;

    // Call the NEWUOA function and find the minimum point
    int status = newuoa( n, npt, objfun, &data, x, rhobeg, rhoend,
                         iprint, maxfun, w );

    if (status != NEWUOA_SUCCESS)
    {
        fprintf( f, "warning: find_LoS_at_r: "
                    "%s\n", newuoa_reason( status ) );
    }

    // Pass the solution out of this function
    psr_angle th, ph;
    set_psr_angle_rad( &th, x[0] );
    set_psr_angle_rad( &ph, x[1] );
    set_point_sph( end_pt, init_pt->r, &th, &ph, POINT_SET_ALL );
}


int get_fieldline_type( point *X, pulsar *psr, double tmult, FILE *f )
{
    int    rL_norm = 0;
    double rL_lim  = 1.0;  // Stop at the light cylinder

    int stop_type;
    point far_pt;

    stop_type = farpoint( X, psr, tmult, f, rL_norm, rL_lim, &far_pt );

    if (stop_type == STOP_FOUND)
        return CLOSED_LINE;
    else if (stop_type == STOP_EXCEED)
        return OPEN_LINE;
    else
    {
        fprintf( stderr, "error: get_fieldline_type: "
                         "unrecognised stop type (=%d)\n", stop_type );
        exit(EXIT_FAILURE);
    }

}


int find_emission_point_elevator( pulsar *psr, psr_angle *phase,
        int direction, point *init_pt, point *emit_pt, FILE *f )
/* This function attempts to find the emission point by using the "elevator"
 * method. I call it the "elevator" method, because it searches for solutions
 * up and down the set of points that satisfy the "line of sight" criterion,
 * which set forms a (non-straight) line stretching from the pulsar outwards.
 * Maybe "space elevator" would have been more apt, but that would've made the
 * function name tediously long.
 *
 * The "elevator" method starts at the supplied initial point (e.g. from the
 * find_approx_emission_point() function). It then uses NEWUOA to search for
 * a point at the same radius r that satisfies the "line of sight" criterion,
 * namely, that the velocity field is parallel to the line of sight. Then,
 * having found that point, it evaluates it for the "last open field line"
 * criterion.
 *
 * The same process is carried out at radius r/2. If the "r/2" field line is
 * of a different type as the "r" field line, then we know that the last open
 * field line is somewhere between them, and we can repeat this whole process
 * in a kind of elaborate "bisection" way, until the correct radius is found.
 *
 * If, however, the "r/2" field line is of the same type as the "r" line, then
 * we try (rL+r)/2 (i.e. halfway to the light cylinder). If this field line is
 * also of the same type, then we try r/4, and then (3*rL+r)/2, etc, until we
 * find a radius at which the corresponding field line is of the opposite
 * type. Then we proceed as explained in the previous paragraph.
 *
 * This function returns 0 if it thinks it has successfully found the correct
 * point, and -1 if it couldn't find anything within the bounds [rp,rL-rp],
 * where rp is the pulsar radius, and rL is the light cylinder radius.
 */
{
    // Make a few temporary points
    point half_down_pt, half_up_pt;

    // Find the point at this radius which satisfies the LoS criterion
    point rlo_pt, rhi_pt, temp_pt;
    find_LoS_at_r( init_pt, psr, phase, direction, &rlo_pt, NULL );

    // This first point is, initially, both the high and low point,
    // since we don't know which side of it our solution is.
    copy_point( &rlo_pt, &rhi_pt );

    // And find its type:
    int rlo_type  = -1; // Just initialise to a dummy value
    int rhi_type  = -1; //  "
    int temp_type = -1; //  "

    double tmult = 0.01;
    rlo_type = get_fieldline_type( &rlo_pt, psr, tmult, NULL );

    // Iterate both inwards and outwards until an opposite-type field line
    // is found
    while (1)
    {
        /* Go down the elevator, a half radius at a time */

        // Declare the search a failure, if we dip below the pulsar radius,
        // or (equivalently) come within a pulsar radius of the light cylinder
        if (rlo_pt.r/2 < psr->r)
            return -1;

        // Save off the previous lowest point to "temp_pt"
        copy_point( &rlo_pt, &temp_pt );
        temp_type = rlo_type;

        // Get the point at half the radius from the previous lowest point
        set_point_sph( &half_down_pt, temp_pt.r/2.0, &(temp_pt.th),
                       &(temp_pt.ph), POINT_SET_ALL );

        // Find the LoS criterion point at this new, lower, radius
        find_LoS_at_r( &half_down_pt, psr, phase, direction, &rlo_pt, NULL );

        // And the field line type there
        rlo_type = get_fieldline_type( &rlo_pt, psr, tmult, NULL );

        // If they differ, then we have a match!
        if (rlo_type != temp_type)
        {
            // Set the high point to the temp point, and quit the while loop
            copy_point( &temp_pt, &rhi_pt );
            break;
        }

        /* Otherwise, repeat the above, but for successively larger radii */

        // Save off the previous highest point to "temp_pt"
        copy_point( &rhi_pt, &temp_pt );
        temp_type = rhi_type;

        // Get the point at half the radius from the previous lowest point
        set_point_sph( &half_up_pt, (temp_pt.r + psr->rL)/2.0, &(temp_pt.th),
                       &(temp_pt.ph), POINT_SET_ALL );

        // Find the LoS criterion point at this new, lower, radius
        find_LoS_at_r( &half_up_pt, psr, phase, direction, &rhi_pt, NULL );

        // And the field line type there
        rhi_type = get_fieldline_type( &rhi_pt, psr, tmult, NULL );

        // If they differ, then we have a match!
        if (rhi_type != temp_type)
        {
            // Set the low point to the temp point, and quit the while loop
            copy_point( &temp_pt, &rlo_pt );
            break;
        }
    }

    /* At this point, we should have an upper and lower bound for radii,
       between which we expect the solution to be found. We now do an
       elaborate kind of bisection iteration, to zone in on the correct point.
    */

    double mid_r;
    point mid_pt;
    while (1)
    {
        // We quit this loop if we've reached machine precision for the radius
        mid_r = (rlo_pt.r + rhi_pt.r) / 2.0;
        if ((rlo_pt.r >= mid_r) || (mid_r >= rhi_pt.r))
        {
            // Declare the lower point the winner and exit!
            copy_point( &rlo_pt, emit_pt );
            break;
        }

        // Set the temp point to half-way between the upper and lower points
        spherical_midpoint( &rlo_pt, &rhi_pt, &mid_pt, POINT_SET_ALL );

        // As before, find the LoS criterion point at this radius
        find_LoS_at_r( &mid_pt, psr, phase, direction, &temp_pt, NULL );

        // And get it's field line type
        temp_type = get_fieldline_type( &temp_pt, psr, tmult, NULL );

        // Have the new point replace whichever of the two straddling points
        // it has the same type as
        if (temp_type == rlo_type)
            copy_point( &temp_pt, &rlo_pt );
        else /* if (temp_type == rhi_type) */
            copy_point( &temp_pt, &rhi_pt );

        // Print out the temp point, if requested
        if (f != NULL)
        {
            int retardation = 1; // Turn on retardation effects
            double cost_lofl = psr_cost_lofl( &temp_pt, psr );
            double cost_los  = psr_cost_los( &temp_pt, psr, phase, direction,
                                             retardation );
            fprintf( f, "%.15e %.15e %.15e %.15e %.15e\n",
                        temp_pt.x[0], temp_pt.x[1], temp_pt.x[2],
                        cost_lofl, cost_los );
        }
    }

    return 0; // = successful
}
