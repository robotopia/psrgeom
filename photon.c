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
    set_point_xyz( &pn->V, pn->V.x[0], pn->V.x[1], pn->V.x[2],
            POINT_SET_TH | POINT_SET_PH );

    // Calculate the particle's trajectory curvature
    pn->curvature = calc_curvature( &pn->V, &pn->A );

    // Calculate the particles Lorentz factor, assuming that the photon's
    // frequency equals the critical frequency
    pn->gamma = calc_crit_gamma( pn->freq, pn->curvature );

    // Calculate the observed phase (incl. retardation)
    psr_angle dph;
    calc_retardation( pt, psr, &pn->V, &dph, &(pn->retarded_LoS) );

    if (psr->spin == SPIN_POS)
        reverse_psr_angle( &(pn->retarded_LoS.ph), &pn->phase );
    else
        copy_psr_angle( &(pn->retarded_LoS.ph), &pn->phase );

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


double weight_photon_total( photon *pn, pulsar *psr, point *foot_pt,
        double height, int pulse_number, double index, int weight_flags )
/* This function combines the weights from all the other weight functions.
 * To choose which weights get applied, set the WEIGHT_FLAG, which can take
 * the following (logical) values:
 *
 *   WEIGHT_GAMMA | WEIGHT_POWER | WEIGHT_SPARK | WEIGHT_LDENS
 *
 * All of them combined can be abbreviated to WEIGHT_TOTAL
 */
{
    double w = 1.0;

    if (weight_flags & WEIGHT_GAMMA)
        w *= weight_photon_by_gamma_distr( pn, index );

    if (weight_flags & WEIGHT_POWER)
        w *= weight_photon_by_power( pn ) * 1.0e60; // (to avoid underflow)

    if (weight_flags & WEIGHT_SPARK)
        w *= weight_photon_by_spark( foot_pt, pn, psr, height, pulse_number );

    if (weight_flags & WEIGHT_LDENS)
        w *= weight_photon_by_line_density( foot_pt, psr );

    return w;
}


void climb_and_emit_photon( pulsar *psr, point *foot_pt, double height,
        double freq, photon *pn )
/* This function starts from a given footpoint, climbs the field line to a
 * given height, and calculates the properties of the photon emitted there
 *
 * Inputs:
 *   pulsar *psr     : pulsar struct
 *   point  *foot_pt : the footpoint of the magnetic field line to be climbed
 *   double  height  : the distance to climb up the field line
 *   double  freq    : the frequency of the emitted photon
 * Outputs:
 *   photon *pn      : the photon struct
 */
{
    // Climb the field line up to the speficied emission point
    point emit_pt;
    B_large_step( foot_pt, psr, height, DIR_OUTWARD, &emit_pt );

    // Calculate the properties of the emitted point
    emit_pulsar_photon( psr, &emit_pt, freq, pn );
}


void bin_photon_freq( photon *pn, double weight, double fmin,
        double fbinwidth, int nfbins, double *farray )
/* This function bins a photon into the appropriate frequency bin in FARRAY,
 * weighting it according to WEIGHT.
 *
 * This function assumes that FARRAY points to already allocated memory, and
 * that it has NFBINS elements. Any photon falling outside of the array bounds
 * is ignored.
 */
{
    int fbin = (int)floor( (pn->freq - fmin) / fbinwidth );
    if ((0 <= fbin) && (fbin < nfbins))
        farray[fbin] += weight;
}


void bin_photon_beam( photon *pn, pulsar *psr, double weight, double xmin,
        double ymin, double xbinwidth, double ybinwidth, int nxbins,
        int nybins, double **xyarray )
/* This function bins a photon into the appropriate beam bin in XYARRAY,
 * weighting it according to WEIGHT.
 *
 * This function assumes that XYARRAY points to already allocated memory, and
 * that it has [NXBINS][NYBINS] elements. Any photon falling outside of the
 * array bounds is ignored.
 *
 * When calculating the x and y coordinates, angular units of degrees are
 * assumed.
 *
 * The x and y values are calculated from the polar coordinates of the
 * photon's observed emission direction (taking into account retardation)
 * in magnetic coordinates.
 */
{
    // Convert the beam direction to magnetic coordinates
    point mag;
    obs_to_mag_frame( &(pn->retarded_LoS), psr, NULL, &mag );

    // Calculate the x and y coordinates of the beam point
    double x = mag.th.deg * mag.ph.cos;
    double y = mag.th.deg * mag.ph.sin;

    int xbin = (int)floor( (x - xmin) / xbinwidth );
    int ybin = (int)floor( (y - ymin) / ybinwidth );

    // Add the weight to the array
    if ((0 <= xbin) && (xbin < nxbins) &&
        (0 <= ybin) && (ybin < nybins))
    {
        xyarray[xbin][ybin] += weight;
    }
}


void bin_photon_pulsestack( photon *pn, double weight, int pulse_number,
        int pulsemin, int npulses, double phasemin, double phasebinwidth,
        int nphasebins, double **pulsestack )
/* This function bins a photon into the appropriate pulsestack bin in
 * PULSESTACK, weighting it according to WEIGHT.
 *
 * This function assumes that PULSESTACK points to already allocated memory,
 * and that it has [NPULSES][NPHASEBINS] elements. Any photon falling outside
 * of the array bounds is ignored.
 *
 * The rotation phase is assumed to be in degrees.
 */
{
    int phasebin = (int)floor( (pn->phase.deg - phasemin) / phasebinwidth );
    int pulsebin = pulse_number - pulsemin;

    // Add the weight to the array
    if ((0 <= phasebin) && (phasebin < nphasebins) &&
        (0 <= pulsebin) && (pulsebin < npulses))
    {
        pulsestack[pulsebin][phasebin] += weight;
    }
}


void bin_photon_profile( photon *pn, double weight, double phasemin,
        double phasebinwidth, int nphasebins, double *profile )
/* This function bins a photon into the appropriate profile bin in PROFILE,
 * weighting it according to WEIGHT.
 *
 * This function assumes that PROFILE points to already allocated memory,
 * and that it has NPHASEBINS elements. Any photon falling outside of the
 * array bounds is ignored.
 *
 * The rotation phase is assumed to be in degrees.
 */
{
    int phasebin = (int)floor( (pn->phase.deg - phasemin) / phasebinwidth );

    // Add the weight to the array
    if ((0 <= phasebin) && (phasebin < nphasebins))
        profile[phasebin] += weight;
}
