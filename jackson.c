/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2018
 *
 * Description:
 *   This source file implements functions found within Jackson's "Classical
 *   Electrodynamics", 3rd ed.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psrgeom.h"


void particle_beam_intensity( double freq, double gamma, psr_angle *theta,
        double rho, double *Ipos, double *Ineg )
/* From Jackson, 3rd ed. pg 706. Problem 14.25:
 *
 * For a relativistic particle moving in a path with instantaneous radius of
 * curvature ρ, the frequency-angle spectra of radiations with positive and
 * negative helicity are
 *
 *    d²I±    e²  / ωρ \² / 1      \² |                   θ                  |²
 *   ----- = ---- |----|  |--- + θ²|  | K_{2/3}(ξ) ± ------------ K_{1/3}(ξ) |
 *   dω dΩ   6π²c \ c  /  \ γ²     /  |              √(1/γ² + θ²)            |
 *
 * where ξ is defined on pg 678 in Eq. (14.76):
 *
 *       ωρ /  1      \ 3/2
 *   ξ = -- | --- + θ²|
 *       3c \  y²     /
 *
 * For now, this function is implemented without the constant out the front,
 * since the important thing here is the shape of the pulsar, not its absolute
 * flux density.
 *
 * The first term and the second term in the brackets are treated separately,
 * as they calculate emission parallel to the synchrotron plane and emission
 * perpendicular to the emission plane respectively. Moreover, the answers are
 * given in terms of the electric field E, not the intensity I ∝ E².
 *
 * Inputs:
 *     double freq     : the radiation frequency, in Hz
 *     double gamma    : the particle's Lorentz factor
 *     point *V        : the (instantaneous) velocity of the particle
 *     point *A        : the (instantaneous) acceleration of the particle
 *     point *LoS      : the direction of the line of sight
 *     double rho      : the radius of curvature of the particle's trajectory,
 *                       in m
 * Outputs:
 *     double *Ipos    : the intensity with positive helicity (a.u.)
 *     double *Ineg    : the intensity with negative helicity (a.u.)
 */
{
    double gt      = 1.0/(gamma*gamma) + theta->rad*theta->rad;
    double sqrt_gt = sqrt(gt);    // i.e. √(gt)
    double gt3_2   = gt*sqrt_gt;  // i.e. (gt)^(3/2)

    double w    = 2.0*PI*freq;             // = ω, the frequency in radians
    double wp_c = w*rho / SPEED_OF_LIGHT;  // = ωρ/c
    double xi   = wp_c * gt3_2 / 3.0;      // = ξ

    double K13, K23;
    bessik(xi, 1.0/3.0, NULL, &K13, NULL, NULL);
    bessik(xi, 2.0/3.0, NULL, &K23, NULL, NULL);

    double A     = wp_c * wp_c * gt * gt;
    double B     = K13 * theta->rad / sqrt_gt;
    double Cpos  = K23 + B;
    double Cneg  = K23 - B;

    *Ipos = A*Cpos*Cpos;
    *Ineg = A*Cneg*Cneg;
}

/* The following scale factor is probably wrong since I used SI units and not
 * cgs units which are used in Jackson */
#define  SCALE  6.81381303799731e-48  /* (e²)/(4πc) */

double single_particle_power( point *n, point *beta, point *beta_dot )
/* From Jackson, 3rd ed. pg 669, Eq (14.38):
 *
 * In the relativistic limit (γ ≫ 1), the angular distribution [of a particle
 * accelerated instantaneously in a direction perpendicular to the velocity]
 * can be written approximately
 *
 *    dP(t′)    e²  |n × {(n - β) × dβ/dt}|²
 *    —————— ≃ ———  ————————————————————————
 *      dΩ     4πc         (1 - n∙β)⁵
 *
 * Only the Cartesian coordinates (i.e. x[0], etc) of the points are used.
 */
{
    double n_minus_B[] = { n->x[0] - beta->x[0],
                           n->x[1] - beta->x[1],
                           n->x[2] - beta->x[2] };
    
    double cross1[] = { n_minus_B[1] * beta_dot->x[2] -
                        n_minus_B[2] * beta_dot->x[1],
                        n_minus_B[2] * beta_dot->x[0] -
                        n_minus_B[0] * beta_dot->x[2],
                        n_minus_B[0] * beta_dot->x[1] -
                        n_minus_B[1] * beta_dot->x[0] };
    
    double cross2[] = { n->x[1] * cross1[2] - n->x[2] * cross1[1],
                        n->x[2] * cross1[0] - n->x[0] * cross1[2],
                        n->x[0] * cross1[1] - n->x[1] * cross1[0] };

    double num = cross2[0]*cross2[0] +
                 cross2[1]*cross2[1] +
                 cross2[2]*cross2[2];

    double n_dot_B = n->x[0] * beta->x[0] +
                     n->x[1] * beta->x[1] +
                     n->x[2] * beta->x[2];

    double den1 = 1.0 / (1.0 - n_dot_B);

    double den = den1*den1*den1*den1*den1;

    return SCALE * num * den;
}

#undef SCALE

/* The following scale factor is probably wrong since I used SI units and not
 * cgs units which are used in Jackson */
#define  SCALE  6.06511156693326e-64  /* (2e²)/(πc³) */

double single_particle_power_perp( double gamma, psr_angle *th, psr_angle *ph,
        double vdot )
/* From Jackson, 3rd ed. pg 671, Eq (14.45):
 *
 * In the relativistic limit (γ ≫ 1), the angular distribution [of a particle
 * accelerated instantaneously in a direction perpendicular to the velocity]
 * can be written approximately
 *
 *    dP(t′)   2 e²        |v̇|²     [     4γ²θ²cos²ϕ  ]
 *    —————— ≃ — —— γ⁶ ———————————  [ 1 - ——————————— ]
 *      dΩ     π c³    (1 + γ²θ²)³  [     (1 + γ²θ²)² ]
 */
{
    double g2    = gamma*gamma;
    double g6    = g2*g2*g2;
    double g2th2 = g2 * th->rad * th->rad;
    double den   = 1.0/(1.0 + g2th2);
    double den2  = den*den;
    double den3  = den*den2;

    return SCALE*g6*vdot*vdot*den3*(1.0 - 4.0*g2th2*ph->cos*ph->cos*den2);
}

#undef SCALE


/* The following scale factor is probably wrong since I used SI units and not
 * cgs units which are used in Jackson */
#define  SCALE  6.35136998062667e-64  /* (2e²)/(3c³) */

double single_particle_power_total( double gamma, double vdot )
/* From Jackson, 3rd ed. pg 671, Eq (14.46):
 *
 * In the relativistic limit (γ ≫ 1), the angular distribution of a particle
 * accelerated instantaneously in a direction perpendicular to the velocity
 * is given in Eq (14.45). Integrating over all angles gives the total power
 *
 *             2 e²|v̇|²
 *    dP(t′) ≃ — ————— γ⁴
 *             3   c³
 */
{
    return SCALE*vdot*vdot*gamma*gamma*gamma*gamma;
}

#undef SCALE


double calc_crit_freq( double gamma, double curvature )
/* Implementation of the formula (see Jackson, 3rd edition, Eqn 14.81):
 *
 *          3γ³cκ
 *    f_c = -----
 *           4π
 *
 * where γ = "gamma", the Lorentz factor,
 *       κ = "curvature" (in inverse metres),
 *       c = the speed of light (m/s)
 */
{
    return 0.75*gamma*gamma*gamma*SPEED_OF_LIGHT*curvature/PI;
}


double calc_crit_gamma( double crit_freq, double curvature )
/* Implementation of the formula (see Jackson, 3rd edition, Eqn 14.81):
 *
 *        |  4πf  |(1/3)
 *    γ = |-------|
 *        |  3cκ  |
 *
 * where γ = "gamma", the Lorentz factor,
 *       f = the critical frequency
 *       κ = "curvature" (in inverse metres),
 *       c = the speed of light (m/s)
 */
{
    return cbrt( 4.0*PI*crit_freq / (3.0*SPEED_OF_LIGHT*curvature) );
}

