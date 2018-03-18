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

void fitwidth( pulsar *psr, int direction, double width_rad,
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
 */
{
    // Set up initial guess emission point
    point     approx_init_pt;
    psr_angle init_ph;
    set_psr_angle_rad( &init_ph, 0.0 );
    find_approx_emission_point( psr, &init_ph, direction, &approx_init_pt );

    // Find the true emission point at this inital phase
    int status;
    status = find_emission_point_elevator( psr, init_ph
}
