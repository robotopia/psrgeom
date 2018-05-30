/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2018
 *
 * Description:
 *   This source file implements functions relevant to the photon struct.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "psrgeom.h"

void emit_pulsar_photon( pulsar *psr, point *pt, double freq, photon *pn )
/* Sets the properties of PN suitable for a photon emitted with frequency FREQ
 * at point PT around pulsar PSR
 */
{
    // Set the photon's frequency
    pn->freq = freq;

    // Set the photon's source location
    copy_point( pt, &pn->source );

    // Calculate the B, V, and A fields at the source location
    calc_fields( pt, psr, SPEED_OF_LIGHT,
            &pn->B, &pn->V, NULL, &pn->A, NULL, NULL );

    // Calculate the particle's trajectory curvature
    pn->curvature = calc_curvature( &pn->V, &pn->A );

    // Calculate the particles Lorentz factor, assuming that the photon's
    // frequency equals the critical frequency
    pn->gamma = calc_crit_gamma( pn->freq, pn->curvature );

    // Calculate the observed phase (incl. retardation)
    point retarded_LoS;
    psr_angle dph;
    calc_retardation( pt, psr, &pn->V, &dph, &retarded_LoS );

    if (psr->spin == SPIN_POS)
        reverse_psr_angle( &(retarded_LoS.ph), &pn->phase );
    else
        copy_psr_angle( &(retarded_LoS.ph), &pn->phase );

    // Calculate the polarisation angle
    accel_to_pol_angle( psr, &pn->A, &pn->phase, &pn->psi );
}


double weight_photon_by_particle_density( photon *pn )
/* This function produces a "weight" that is proportional to the relative
 * density of particles at the point of emission.
 *
 * The rationale is that, when summing up the contributions from multiple
 * particles, some locations should get downweighted if the density of
 * emitting particles is less there.
 *
 * This should only be applied in the case that the point of emission was
 * drawn from a uniform distribution over some region the magnetosphere. If
 * the locations were drawn from some "2D" sampling of the lines themselves,
 * the functino weight_photon_by line_density() should be used instead.
 */
{
    return pn->B.r;
}


double weight_photon_by_line_density( point *init_pt, pulsar *psr )
/* This function produces a "weight" that is proportional to the relative
 * density of field lines at the point of emission.
 *
 * This should only be applied in the case that the field line was drawn from
 * a uniform distribution over some 2D surface (which the field line
 * intersects). If the emission point was drawn from a uniform sampling over
 * some 3D region of the magnetosphere, the function
 * weight_photon_by_particle_density() should be used instead.
 */
{
    point B;
    calc_fields( init_pt, psr, SPEED_OF_LIGHT,
            &B, NULL, NULL, NULL, NULL, NULL );

    return B.r;
}


double weight_photon_by_power( photon *pn )
/* This function produces a "weight" that is proportional to the total
 * instantaneous power emitted by the particle. This is a function of the
 * particle's Lorentz factor, which is assumed to be that which corresponds to
 * a critical frequency equal to the photon's set frequency.
 */
{
    return single_particle_power_total( pn->gamma, pn->A.r );
}


double weight_photon_by_gamma_distr( photon *pn, double index )
/* This function produces a "weight" that is proportional to the number
 * density of Lorentz factors evaluated at the photon's gamma value. The
 * power law INDEX is supplied by the caller. In this function (unlike
 * neg_power_law_distr()), the index is not assumed to be negative, i.e.
 *
 *   dN/dγ ∝ γᵅ
 */
{
    return pow( pn->gamma, index );
}


double weight_photon_by_spark( point *foot_pt, photon *pn, pulsar *psr,
        double height, int pulse_number )
/* This function produces a "weight" that is proportional to the PDF of the
 * spark configuration at the time a spark would have led to the photon being
 * emitted (and subsequently observed).
 */
{
    // Calculated the fraction number of periods (incl. whole number of
    // rotations) that must have elapsed in order for this photon to be
    // observed. This includes the time it has taken the spark information
    // to climb the field line to the given height.
    double t = (double)pulse_number - (psr->spin)*(pn->V.ph.deg)/360.0;
    t += height*SPEED_OF_LIGHT;

    return spark_profile( psr, t, foot_pt );
}

