/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2018
 *
 * Description:
 *   This source file implements the functions required to find the emission
 *   height that corresponds to a given (main or inter-) pulse width.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <newuoa.h>
#include "psrgeom.h"

int fitwidth( pulsar *psr, int direction, double width_rad,
              psr_angle *ph1, psr_angle *ph2, FILE *f )
/* This function finds the pair of rotation phases that are width_rad radians
 * apart, and that have the same emission height. It does this by first
 * finding the emission height for the fiducial point (or anti-fiducial point
 * in the direction = DIR_INWARD case), and using that point as the initial
 * guess for the two points at ± width_rad/2 on either side. These two points
 * are then slid up and down in rotation phase until the heights are the same.
 *
 * Initially, it is not known whether the solution is at a lower or higher
 * phase than the first guess. We try lower phases if the lower of the two
 * phases has a lower height than the upper phase:
 *
 *   |
 * h | .                                            .
 * e |    .                                   .
 * i |       .                           .
 * g |            .                   .*
 * h |                *      .           φ2
 * t |              φ1
 *   |
 *   +-------------------------------------------------->
 *                      rotation phase, φ
 *
 * Inputs:
 *   pulsar *psr      : a pointer to a pulsar struct
 *   int direction    : either DIR_OUTWARD or DIR_INWARD (for particles
 *                      flowing along or against the magnetic field,
 *                      respectively)
 *   double width_rad : the phase difference between the two phases
 *   FILE *f          : where to write intermediate steps
 * Outputs:
 *   psr_angle *ph1   : the lower rotation phase
 *   psr_angle *ph2   : the upper rotation phase
 * Returns:
 *   The return value is one of the following:
 *     EMIT_PT_TOO_LOW, EMIT_PT_FOUND, EMIT_PT_TOO_HIGH
 */
{
    // Error check direction
    if (direction != DIR_OUTWARD && direction != DIR_INWARD)
    {
        fprintf( stderr, "error: fitwidth: direction must be either "
                         "DIR_OUTWARD or DIR_INWARD\n" );
        exit(EXIT_FAILURE);
    }

    // Set up initial guess emission point
    point     approx_init_pt;
    point     init_pt;
    psr_angle init_ph;

    set_psr_angle_rad( &init_ph, 0.0 );
    find_approx_emission_point( psr, &init_ph, direction, &approx_init_pt );

    // Find the true emission point at this initial phase
    if (direction == DIR_INWARD)
        set_psr_angle_deg( &init_ph, 180.0 );
    int status;
    status = find_emission_point_elevator(
                 psr, &init_ph, direction, &approx_init_pt, &init_pt, NULL );

    // Assert that a suitable point was found
    if (status != EMIT_PT_FOUND)   return status;

    // Keep track of the midpoint between the two phases
    double mid_ph = init_ph.rad;

    // Set ph1 and ph2 to straddle the init_pt with the specified width
    set_psr_angle_rad( ph1, mid_ph - 0.5*width_rad );
    set_psr_angle_rad( ph2, mid_ph + 0.5*width_rad );

    // Evaluate the height at these two points to determine the first shift
    // direction
    point ph1_pt;
    point ph2_pt;

    // ph1:
    status = find_emission_point_elevator(
                 psr, ph1, direction, &init_pt, &ph1_pt, NULL );
    if (status != EMIT_PT_FOUND)   return status;

    // ph2:
    status = find_emission_point_elevator(
                 psr, ph2, direction, &init_pt, &ph2_pt, NULL );
    if (status != EMIT_PT_FOUND)   return status;

    // Check the extremely unlikely possibility that we got it first go
    if (ph1_pt.r == ph2_pt.r)      return EMIT_PT_FOUND;

    // Otherwise, get the shift direction
    int init_shift_dir = (ph1_pt.r < ph2_pt.r ? -1 : 1);

    // Print out our progress so far, if requested
    if (f != NULL)
    {
        // First, print a column header line
        fprintf( f, "# φ1+φ2/2  h(φ1)  h(φ2)\n" );

        // Then print the result
        fprintf( f, "%.15e %.15e %.15e\n", mid_ph, ph1_pt.r, ph2_pt.r );
    }

    // Loop over different phases (always in one direction) until we've gone
    // too far.
    double ph_step = 0.1 * width_rad;
    int shift_dir  = init_shift_dir;
    while (shift_dir == init_shift_dir)
    {
        // Step to the next "mid" phase
        mid_ph += ph_step * init_shift_dir;

        // Evaluate the heights at the two (new) phases
        set_psr_angle_rad( ph1, mid_ph - 0.5*width_rad );
        set_psr_angle_rad( ph2, mid_ph + 0.5*width_rad );

        // ph1:
        status = find_emission_point_elevator(
                     psr, ph1, direction, &ph1_pt, &ph1_pt, NULL );
        if (status != EMIT_PT_FOUND)   return status;

        // ph2:
        status = find_emission_point_elevator(
                     psr, ph2, direction, &ph2_pt, &ph2_pt, NULL );
        if (status != EMIT_PT_FOUND)   return status;

        // Re-evaluate the shift direction
        shift_dir = (ph1_pt.r < ph2_pt.r ? -1 : 1);

        // Print out progress so far
        if (f != NULL)
            fprintf( f, "%.15e %.15e %.15e\n", mid_ph, ph1_pt.r, ph2_pt.r );
    }

    // Knowing that the solution lies between the last two sets of phases that
    // we tried, we can use bisection to zero in on the correct answer
    double lo_mid_ph, hi_mid_ph;
    if (init_shift_dir == -1)
    {
        lo_mid_ph = mid_ph;
        hi_mid_ph = mid_ph + ph_step;
    }
    else
    {
        lo_mid_ph = mid_ph - ph_step;
        hi_mid_ph = mid_ph;
    }

    while (1)
    {
        // Evaluate the phase in the middle
        mid_ph += (lo_mid_ph + hi_mid_ph) / 2.0;

        // Check if machine precision has been reached
        if (lo_mid_ph >= mid_ph || mid_ph >= hi_mid_ph)
            break;

        // Evaluate the heights at the two (new) phases
        set_psr_angle_rad( ph1, mid_ph - 0.5*width_rad );
        set_psr_angle_rad( ph2, mid_ph + 0.5*width_rad );

        // ph1:
        status = find_emission_point_elevator(
                     psr, ph1, direction, &ph1_pt, &ph1_pt, NULL );
        if (status != EMIT_PT_FOUND)   return status;

        // ph2:
        status = find_emission_point_elevator(
                     psr, ph2, direction, &ph2_pt, &ph2_pt, NULL );
        if (status != EMIT_PT_FOUND)   return status;

        // Print out progress so far
        if (f != NULL)
            fprintf( f, "%.15e %.15e %.15e\n", mid_ph, ph1_pt.r, ph2_pt.r );

        // Bisect!
        if      (ph1_pt.r == ph2_pt.r)     break;
        else if (ph1_pt.r <  ph2_pt.r)     hi_mid_ph = mid_ph;
        else /* (ph1_pt.r >  ph2_pt.r) */  lo_mid_ph = mid_ph;
    }

    return EMIT_PT_FOUND;
}
