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
 */
{
    return pn->B.r;
}

double weight_photon_by_power( photon *pn )
/* This function produces a "weight" that is proportional to the total
 * instantaneous power emitted by the particle. This is a function of the
 * particle's Lorentz factor, which is assumed to be that which corresponds to
 * a critical frequency equal to the photon's set frequency.
 */
{
    double vdot  = pn->A.r;
    double gamma = calc_crit_gamma( pn->freq, pn->curvature );

    return single_particle_power_total( gamma, vdot );
}

//double weight_photon_by_footpoint( photon *pn, pulsar *psr )
//{
//}
