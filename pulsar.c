/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the struct for a pulsar, which includes
 *   information about its period, its location, its size, and its viewing
 *   geometry.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "psrgeom.h"

double light_cylinder( double P )
/* Calculates the light cylinder radius.
 *
 * Inputs:
 *    double P : the spin period of the pulsar (sec)
 *
 * Outputs:
 *    double [return] : the light cylinder radius (m)
 */
{
    return SPEED_OF_LIGHT * P / (2.0*PI);
}


void set_pulsar( pulsar *psr, psr_angle *ra, psr_angle *dec, double P,
        double r, psr_angle *al, psr_angle *ze )
{
    // Only set non-NULL angle arguments
    if (ra)
        copy_psr_angle( ra,  &psr->ra  );
    if (dec)
        copy_psr_angle( dec, &psr->dec );

    // Set the pulsar period and the light cylinder radius
    set_pulsar_period( psr, P );

    // Set the pulsar radius
    psr->r = r;

    if (al)
        copy_psr_angle( al, &psr->al );
    if (ze)
        copy_psr_angle( ze, &psr->ze );

    // Set to Deutsch field by default
    psr->field_type = DEUTSCH;
}

pulsar *create_pulsar( psr_angle *ra, psr_angle *dec, double P, double r,
        psr_angle *al, psr_angle *ze )
{
    // Dynamically allocate memory for a pulsar struct
    pulsar *psr = (pulsar *)malloc( sizeof(pulsar) );

    // Populate the values
    set_pulsar( psr, ra, dec, P, r, al, ze );

    // Return the pointer to the new pulsar struct
    return psr;
}

void destroy_pulsar( pulsar *psr )
{
    free( psr );
}

void set_pulsar_period( pulsar *psr, double P )
/* Sets the pulsar period, the light cylinder radius, and the spin direction.
 * The spin direction is set according to the sign of the period. If a
 * positive period is given, then the spin direction is positive in the sense
 * that the pulsar spins counterclockwise when viewed from "above" (i.e.
 * positive z).
 */
{
    psr->P    = fabs(P);

    psr->spin = (P >= 0.0 ? SPIN_POS : SPIN_NEG);

    psr->rL   = light_cylinder( psr->P );
    psr->rL2  = psr->rL * psr->rL;

    set_psr_angle_rad( &psr->Om, 2.0*PI/psr->P );
}


void line_of_sight( pulsar *psr, psr_angle *phase, point *LoS )
/* This function finds the line-of-sight unit vector (LoS) for a given pulsar
 * geometry and rotation phase angle. The LoS is given in rotating frame
 * coordinates, where the z-axis is the rotation axis, and where the magnetic
 * axis is fixed to lie in the xz-plane.
 *
 * The point LoS is set with the POINT_SET_ALL flag, and is guaranteed to have
 * length 1, both in the sense of sqrt(x^2 + y^2 + z^2) = 1, and r = 1.
 *
 * Inputs:
 *   pulsar    *psr    : the pulsar geometry
 *   psr_angle *phase  : the rotation phase
 * Outputs:
 *   point     *LoS    : the line of sight unit vector
 */
{
    // A positive rotation of the pulsar equals an apparent reverse rotation
    // of the line of sight in the rotating frame
    double apparent_phase_deg =
        (psr->spin == SPIN_POS ? -phase->deg : phase->deg);

    psr_angle apparent_phase;
    set_psr_angle_deg( &apparent_phase, apparent_phase_deg );

    // Set the spherical coordinates for the line of sight
    double     r  = 1.0; // A unit vector
    psr_angle *th = &(psr->ze);
    psr_angle *ph = &apparent_phase;

    set_point_sph( LoS, r, th, ph, POINT_SET_ALL );
}

void pol_zero( pulsar *psr, psr_angle *phase, point *pz )
/* This function finds the vector that is defined to be the 0° reference angle
 * for the PA swing. It is equivalent to the line of sight rotated 90° towards
 * the pulsar's north pole.
 *
 * The point pz is set with the POINT_SET_ALL flag, and is guaranteed to have
 * length 1, both in the sense of sqrt(x^2 + y^2 + z^2) = 1, and r = 1.
 *
 * Inputs:
 *   pulsar    *psr    : the pulsar geometry
 *   psr_angle *phase  : the rotation phase
 * Outputs:
 *   point     *pz     : the pol-zero reference vector
 */
{
    // A positive rotation of the pulsar equals an apparent reverse rotation
    // of the line of sight in the rotating frame
    double apparent_phase_deg =
        (psr->spin == SPIN_POS ? -phase->deg : phase->deg);

    psr_angle apparent_phase;
    set_psr_angle_deg( &apparent_phase, apparent_phase_deg );

    psr_angle LoS90;
    set_psr_angle_deg( &LoS90, psr->ze.deg - 90.0 );

    // Set the spherical coordinates for the line of sight
    double     r  = 1.0; // A unit vector
    psr_angle *th = &LoS90;
    psr_angle *ph = &apparent_phase;

    set_point_sph( pz, r, th, ph, POINT_SET_ALL );
}


void calc_retardation( point *X, pulsar *psr, point *LoS,
        psr_angle *dph, point *retarded_LoS )
/* This function calculates the retardation due to the difference in flight
 * time of photons starting from different locations in the magnetosphere.
 * In particular, it returns a pseudo-line-of-sight, i.e. the line of sight
 * that has been rotated around the pulsar's rotation axis by an amount
 * corresponding to the pulsar's rotation rate and the time difference between
 * a photon starting at X and a photon starting at some reference point (here,
 * we choose the origin).
 *
 * Thus, an emission point closer to the observer than the reference point
 * would be seen at an earlier time. Since the LoS (apparently) rotates in the
 * opposite sense to the pulsar rotation (i.e. in the rotating frame), we can
 * simulate an earlier time of arrival with an artificial rotation of the line
 * of sight in the *same* sense as the pulsar.
 *
 * This function assumes that X's Cartesian coordinates are set, that its
 * radius, r, is set, and that LoS is of length unity.
 *
 * Inputs:
 *   point  *X   : the point at which the retardation is to be calculated
 *   pulsar *psr : the pulsar struct
 *   point  *LoS : the (non-retarded) line of sight
 * Outputs:
 *   point *retarded_angle : the retardation angle
 *   point *retarded_LoS   : the retarded line of sight (can be NULL)
 */
{
    // Calculate the prejection of the position vector of X onto the LoS
    // (we assume LoS has length 1)
    double proj = (X->x[0] * LoS->x[0] +
                   X->x[1] * LoS->x[1] +
                   X->x[2] * LoS->x[2]);

    // Convert the projected length into a photon travel time
    double t = proj / SPEED_OF_LIGHT;

    // Convert the photon travel time into a pulsar rotation phase
    set_psr_angle_rad( dph, t * psr->Om.rad * psr->spin );

    // Rotate the line of sight by that amount (around the z-axis = pulsar's
    // rotation axis), if requested
    if (retarded_LoS != NULL)
        rotate_about_axis( LoS, retarded_LoS, dph, 'z', POINT_SET_ALL );
}


double neg_power_law_distr( double lo, double hi, double index )
/* This function generates a random number between lo and hi that follows a
 * (negative) power law distribution with the supplied index.
 *
 * dN/dx ∝ x⁻ᵅ
 *
 * This function does not seed the random number generator. The caller is
 * responsible for doing that.
 *
 * Inputs:
 *   double lo     : the lower bound for the generated random number
 *   double hi     : the upper bound for the generated random number
 *   double index  : the power law index for the distribution
 *
 * Returns:
 *   double        : the generated random number
 */
{
    // Generate a (uniform) random number between 0 and 1
    double C = RAND(1.0);

    // Convert C to a correctly distributed random number via the cumulative
    // distribution
    double exp  = 1.0 - index;
    double expi = 1.0 / exp;
    double lo_a = pow( lo, exp );
    double hi_a = pow( hi, exp );
    return pow( C*(hi_a - lo_a) + lo_a, expi );
}


void random_point_in_lightcyl( point *rand_pt, pulsar *psr, double frac,
        double z_max_frac )
/* Generate a random point (RAND_PT) within fraction FRAC of the light
 * cylinder of pulsar PSR. The returned point is guaranteed to be outside of
 * the pulsar itself. The absolute value of the point's z coordinate is
 * guaranteed to be not more than Z_MAX_FRAC × FRAC × the light cylinder
 * radius.
 */
{
    double rho   = frac * psr->rL;
    double max_z = z_max_frac * rho;

    rand_pt->r = 0.0;
    while (rand_pt->r < psr->r)
        random_point_in_cyl( rand_pt, rho, max_z );
}


void set_pulsar_carousel( pulsar *psr, int n, psr_angle *s, psr_angle *S,
        int type, double P4 )
/*
 * Parameters:
 *   pulsar    *psr   : The pulsar whose carousel parameters are to be set
 *   int        n     : Number of sparks in the carousel (0 = annulus)
 *   psr_angle  s     : The angular radius of a spark
 *   psr_angle  S     : The angular radius of the whole carousel
 *   int        type  : One of {TOPHAT, GAUSSIAN}
 *   double     P4    : The rotation rate of the carousel (in sec)
 */
{
    psr->csl.n    = n;
    psr->csl.type = type;
    psr->csl.P4   = P4;
    copy_psr_angle( s, &(psr->csl.s) );
    copy_psr_angle( S, &(psr->csl.S) );
}


double spark_profile( pulsar *psr, double t, point *foot_pt )
/* This function calculates the spark profile at the given FOOT_PT on the
 * surface of pulsar PSR at time T.
 *
 * Here, a positive P4 is interpreted to mean a "right-handed" rotation
 * about the magnetic axis, i.e. counterclockwise as viewed from above the
 * surface.
 */
{
    // Convert the time into a rotation of the carousel
    psr_angle csl_rot;
    set_psr_angle_deg( &csl_rot, 360.0*t/(psr->csl.P4) );

    double h = 0.0;     // The "height" of the profile at the given footpoint
    psr_angle spark_ph; // The magnetic azimuth of a spark
    point spark_pt;     // The location of the spark centre on the surface
    psr_angle sep;      // The angular separation (in azimuth) between sparks
    double dist;        // The angular separation between the foot point and
                        // the spark (in rad)
    double dist_norm;   // Same as 'dist', but normalised to the spark size

    if (psr->csl.n == 0)
    {
        dist      = fabs( psr->csl.S.rad - foot_pt->th.rad );
        dist_norm = dist / psr->csl.s.rad;
        
        // Add this spark's contribution
        switch (psr->csl.type)
        {
            case TOPHAT:
                if (dist_norm <= 1.0)
                    h += 1.0;
                break;
            case GAUSSIAN:
                h += exp( -dist_norm*dist_norm );
                break;
            default:
                fprintf( stderr, "error: spark_profile: unrecognised "
                        "carousel type\n" );
                exit(EXIT_FAILURE);
                break;
        }
    }
    else if (psr->csl.n > 0)
    {
        set_psr_angle_deg( &sep, 360.0 / psr->csl.n );

        // For each spark...
        int n;
        for (n = 0; n < psr->csl.n; n++)
        {
            // Calculate how far away the foot_point is from the spark centre
            set_psr_angle_deg( &spark_ph, (double)n*sep.deg + csl_rot.deg );
            set_point_sph( &spark_pt, psr->r,
                                      &(psr->csl.S),
                                      &spark_ph,
                                      POINT_SET_ALL );
            dist      = acos( norm_dot( &spark_pt, foot_pt ) );
            dist_norm = dist / psr->csl.s.rad;

            // Add this spark's contribution
            switch (psr->csl.type)
            {
                case TOPHAT:
                    if (dist_norm <= 1.0)
                        h += 1.0;
                    break;
                case GAUSSIAN:
                    h += exp( -dist_norm*dist_norm );
                    break;
                default:
                    fprintf( stderr, "error: spark_profile: unrecognised "
                                     "carousel type\n" );
                    exit(EXIT_FAILURE);
                    break;
            }
        }
    }
    else
    {
        fprintf( stderr, "error: spark_profile: carousel cannot have a "
                "negative number of sparks\n" );
        exit(EXIT_FAILURE);
    }

    return h;
}
