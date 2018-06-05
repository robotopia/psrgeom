/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2018
 *
 * Description:
 *   This source file implements the functions for treating distributions of
 *   the Lorentz factors of outflowing particles, aka "gamma distributions".
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psrgeom.h"


double avg_power_single( gamma_distr *gd, double vdot )
/* This function calculates the average power that is output by a single
 * relativistic particle being instantaneously accelerated perpendicular to
 * its velocity:
 *
 *   P̅ = ∫ P(t′,γ) dN/dγ dγ,
 *
 * where P(t′,γ) is given by Jackson Eq. (14.46) and dN/dγ is given by the
 * supplied gamma distribution GD.
 *
 * This function does not include the scalar multiplicative factor in Jackson
 * Eq. (14.46).
 *
 * For power laws, the integration becomes:
 *
 *   P̅ = ∫ P(t′,γ) dN/dγ dγ
 *     = |v̇|² ∫ γ⁴γᵅ dγ
 *     = |v̇|² (γ₂⁵⁺ᵅ - γ₁⁵⁺ᵅ) / (5+α)
 *
 * For normal distributions, it becomes (integrating from -∞ to ∞):
 *
 *   P̅ = ∫ P(t′,γ) dN/dγ dγ
 *     = |v̇|² ∫ γ⁴ Norm(γ,μ,σ) dγ
 *     = |v̇|² (μ⁴ + 6μ²σ² + 3σ⁴)
 */
{
    double avgP = 0.0;
    double v2 = vdot*vdot;
    if (gd->type == POWER_LAW)
    {
        double idx = 5.0 + gd->idx;
        avgP = v2*(pow(gd->g_max, idx) - pow(gd->g_min, idx)) / idx;
    }
    else if (gd->type == NORMAL)
    {
        double m2 = gd->mean * gd->mean;
        double s2 = gd->std  * gd->std;
        avgP = v2*(m2*m2 + 6.0*m2*s2 + 3.0*s2*s2);
    }
    else
    {
        fprintf( stderr, "error: avg_power_single: unrecognised gamma "
                "distribution type\n" );
        exit(EXIT_FAILURE);
    }

    return avgP;
}
