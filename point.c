/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the struct for geometric points defined in
 *   psrgeom.h.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "psrgeom.h"

void set_point_xyz( point *p, double x, double y, double z, int flags )
/* Populate a point structure given x,y,z coordinates
 */
{
    // Set cartesian coordinates
    if (flags & POINT_SET_X)
        p->x[0] = x;
    if (flags & POINT_SET_Y)
        p->x[1] = y;
    if (flags & POINT_SET_Z)
        p->x[2] = z;

    // Calculate cylindrical distance rho squared
    if (flags & POINT_SET_RHOSQ)
        p->rhosq = x*x + y*y;

    // Calculate spherical coordinates
    double r = 0.0;
    if ((flags & POINT_SET_R) || (flags & POINT_SET_TH))
        r = sqrt( x*x + y*y + z*z );
    if (flags & POINT_SET_R)
        p->r = r;
    if (flags & POINT_SET_TH)
        set_psr_angle_cos( &p->th, z/r );
    if (flags & POINT_SET_PH)
        set_psr_angle_rad( &p->ph, atan2(y,x) );
}

void set_point_sph( point *p, double r, psr_angle *th, psr_angle *ph, int flags )
/* Populate a point structure given spherical coordinates
 */
{
    // Set spherical coordinates
    if (flags & POINT_SET_R)
        p->r = r;
    if (flags & POINT_SET_TH)
        copy_psr_angle( th, &p->th );
    if (flags & POINT_SET_PH)
        copy_psr_angle( ph, &p->ph );

    // Calculate cartesian coordinates
    if ((flags & POINT_SET_X) || (flags & POINT_SET_RHOSQ))
        p->x[0] = r * th->sin * ph->cos;
    if ((flags & POINT_SET_Y) || (flags & POINT_SET_RHOSQ))
        p->x[1] = r * th->sin * ph->sin;
    if (flags & POINT_SET_Z)
        p->x[2] = r * th->cos;

    // Calculate cylindrical distance rho squared
    if (flags & POINT_SET_RHOSQ)
        p->rhosq = p->x[0]*p->x[0] + p->x[1]*p->x[1];
}

void set_point_cyl( point *p, double rh, psr_angle *ph, double z, int flags )
/* Populate a point structure given cylindrical coordinates
 */
{
    // Calculate cartesian coordinates
    if (flags & POINT_SET_X)
        p->x[0] = rh * ph->cos;
    if (flags & POINT_SET_Y)
        p->x[1] = rh * ph->sin;
    if (flags & POINT_SET_Z)
        p->x[2] = z;

    // Calculate spherical coordinates
    if (flags & POINT_SET_R)
        p->r = sqrt( rh*rh + z*z );
    if (flags & POINT_SET_TH)
        set_psr_angle_rad( &p->th, atan2( rh, z ) );
    if (flags & POINT_SET_PH)
        copy_psr_angle( ph, &p->ph );

    // Calculate cylindrical distance rho squared
    if (flags & POINT_SET_RHOSQ)
        p->rhosq = rh*rh;
}

void copy_point( point *src, point *dest )
/* Copy a point from one struct to another
 */
{
    dest->x[0]  = src->x[0];
    dest->x[1]  = src->x[1];
    dest->x[2]  = src->x[2];
    dest->r     = src->r;
    copy_psr_angle( &src->th, &dest->th );
    copy_psr_angle( &src->ph, &dest->ph );
    dest->rhosq = src->rhosq;
}


double norm_dot( point *p1, point *p2 )
/* Calculate the normlised dot product of p1 and p2.
 * That is,
 *           p1     p2
 *          ---- . ----
 *          |p1|   |p2|
 * This function assumes that the Cartesian coordinates and the radial
 * coordinate are set.
 */
{
    return ((p1->x[0] * p2->x[0] +
             p1->x[1] * p2->x[1] +
             p1->x[2] * p2->x[2]) / (p1->r * p2->r));
}


void spherical_midpoint( point *p1, point *p2, point *mid_pt, int flags )
/* This function calculates the "spherical midpoint", which is here defined as
 * the point whose radius is midway between the two given points, and whose
 * latitude and longitude are such that if the two given points had the same
 * radius, the midpoint would be on the great circle connecting the two
 * points, and would be equidistant from them.
 *
 * This implementation assumes that both the Cartesian coordinates of points
 * p1 and p2 and their radius, r, have been set.
 */
{
    // Once p1 and p2 have been scaled to the same radius, their vector sum
    // will have the correct latitude and longitute
    point sum_pt;
    set_point_xyz( &sum_pt, p1->x[0]/p1->r + p2->x[0]/p2->r,
                            p1->x[1]/p1->r + p2->x[1]/p2->r,
                            p1->x[2]/p1->r + p2->x[2]/p2->r,
                            POINT_SET_TH | POINT_SET_PH );
    set_point_sph( mid_pt, (p1->r + p2->r) / 2.0,
                           &sum_pt.th,
                           &sum_pt.ph,
                           flags );
}


void random_direction( point *rand_pt )
/* This function generates a uniformly distributed random point on the unit
 * sphere, which is set in rand_pt.
 * This function assumes that srand() has already been called.
 */
{
    double th_rad, ph_rad;
    psr_angle th, ph;

    th_rad = RANDTH;
    ph_rad = RAND(2.0*PI);

    set_psr_angle_rad( &th, th_rad );
    set_psr_angle_rad( &ph, ph_rad );

    set_point_sph( rand_pt, 1.0, &th, &ph, POINT_SET_ALL );
}


void random_direction_bounded( point *rand_pt, double lo_th_rad,
        double hi_th_rad, double lo_ph_rad, double hi_ph_rad )
/* This function generates a uniformly distributed random point on the swath
 * of the unit sphere in the range θ1 < θ < θ2, which is set in rand_pt.
 * (θ1 = "lo_th_rad" and θ2 = "hi_th_rad")
 * This function assumes that srand() has already been called.
 */
{
    double th_rad, ph_rad;
    psr_angle th, ph;

    th_rad = RANDTHAB(lo_th_rad, hi_th_rad);
    ph_rad = RAND(hi_ph_rad - lo_ph_rad) + lo_ph_rad;

    set_psr_angle_rad( &th, th_rad );
    set_psr_angle_rad( &ph, ph_rad );

    set_point_sph( rand_pt, 1.0, &th, &ph, POINT_SET_ALL );
}


void random_spark_footpoint( point *foot_pt_obs, point *foot_pt_mag,
        pulsar *psr, double t )
/* This function produces a random footpoint FOOT_PT within one of the sparks
 * in the carousel of pulsar PSR, at time T.
 * The footpoint is returned in the observer frame (not the magnetic frame).
 */
{
    point pt;
    psr_angle th, ph;

    int    n = psr->csl.n;      // Number of sparks
    double S = psr->csl.S.rad;  // Angular radius of carousel
    double s = psr->csl.s.rad;  // Angular size of sparks
    int type = psr->csl.type;   // Either TOPHAT or GAUSSIAN

    if (type != TOPHAT && type != GAUSSIAN)
    {
        fprintf( stderr, "error: random_spark_footpoint: unrecognised "
                "carousel type\n" );
        exit(EXIT_FAILURE);
    }

    // How the point is chosen depends on whether the carousel is a ring of
    // sparks (n > 0) or a "solid" annulus (n = 0)
    if (n == 0)  /* The annulus case */
    {
        // If the annulus has a "tophat" profile, try to ensure that the
        // resulting footpoint is drawn from a uniform distribution across the
        // allowed "swath"
        if (type == TOPHAT)
        {
            random_direction_bounded( &pt, S-s, S+s, 0.0, 2.0*PI );
            scale_point( &pt, psr->r, &pt );
        }
        // If the annulus has a "gaussian" profile, just draw θ directly from
        // a gaussian distribution and φ from a uniform distribution.
        else /* (type == GAUSSIAN) */
        {
            set_psr_angle_rad( &th, randn( S, s ) );
            set_psr_angle_rad( &ph, RAND(2.0*PI) );
            set_point_sph( &pt, psr->r, &th, &ph, POINT_SET_ALL );
        }
    }
    else if (n > 0) /* The discrete spark case */
    {
        // Generate a random point on a spark centred at the zenith, as
        // if there were only a single spark located there.
        if (type == TOPHAT)
        {
            random_direction_bounded( &pt, 0.0, s, 0.0, 2.0*PI );
            scale_point( &pt, psr->r, &pt );
        }
        else /* (type == GAUSSIAN) */
        {
            double x = randn( 0.0, s );
            double y = randn( 0.0, s );

            set_psr_angle_rad( &th, hypot( x, y ) );
            set_psr_angle_rad( &ph, atan2( y, x ) );

            set_point_sph( &pt, psr->r, &th, &ph, POINT_SET_ALL );
        }

        // Rotate the point down to the correct annulus and around to a random
        // spark
        set_psr_angle_rad( &th, S );
        set_psr_angle_deg( &ph, 360.0*((double)(rand() % n)/(double)n +
                    t/psr->csl.P4) );

        rotate_about_axis( &pt, &pt, &th, 'y', POINT_SET_ALL );
        rotate_about_axis( &pt, &pt, &ph, 'z', POINT_SET_ALL );
    }

    // Convert the footpoint to the observer frame
    if (foot_pt_mag != NULL)
        copy_point( &pt, foot_pt_mag );
    if (foot_pt_obs != NULL)
        mag_to_obs_frame( &pt, psr, NULL, foot_pt_obs );
}


void random_point_in_cyl( point *rand_pt, double max_rho, double max_z )
/* Generate a random point (RAND_PT) within a cylinder of dimensions:
 *   ρ = MAX_RHO,  z = ± MAX_Z
 */
{
    double rho     = max_rho * sqrt(RAND(1.0));
    double z       = max_z * RANDU;
    double ph_rad  = RAND(2.0*PI);
    psr_angle ph;
    set_psr_angle_rad( &ph, ph_rad );

    set_point_cyl( rand_pt, rho, &ph, z, POINT_SET_ALL );
}


void scale_point( point *in, double scale, point *out )
/* Scale a point by a given amount. The in and out pointers can point to the
 * same point. ;-)
 */
{
    int i;
    for (i = 0; i < 3; i++)
        out->x[i] = in->x[i] * scale;
    out->r = in->r * scale;
    out->rhosq = in->rhosq * scale;

    if (in != out)
    {
        copy_psr_angle( &(in->th), &(out->th) );
        copy_psr_angle( &(in->ph), &(out->ph) );
    }
}


void transform_new_xz( point *v, point *new_z, point *new_x, point *new_v )
/* Calculate the vector v in the (new) coordinate system in which new_z is the
 * new z-axis and new_x is the new x-axis. The new y-axis is assumed to be
 * new_x × new_z.
 *
 * new_x and new_z are assumed to be perpendicular to each other.
 *
 * Points are assumed to have spherical coordinates set.
 *
 * The resulting vector is written to new_v.
 */
{
    point nxp, nxpp;
    point vp, vpp;
    psr_angle tmp;

    // Euler rotation #1:
    // Rotate about z-axis to get (virtual) new_z in xz-plane
    reverse_psr_angle( &(new_z->ph), &tmp );

    rotate_about_axis( v,     &vp,  &tmp, 'z', POINT_SET_ALL );
    rotate_about_axis( new_x, &nxp, &tmp, 'z', POINT_SET_ALL );

    // Euler rotation #2:
    // Now rotate about the y-axis to get (virtual) new_z parallel to z-axis
    reverse_psr_angle( &(new_z->th), &tmp );

    rotate_about_axis( &vp,   &vpp,  &tmp, 'y', POINT_SET_ALL );
    rotate_about_axis( &nxp,  &nxpp, &tmp, 'y', POINT_SET_ALL );

    // Euler rotation #3:
    // Finally, rotate about the z-axis again to get new_x parallel to x-axis
    reverse_psr_angle( &(nxpp.ph), &tmp );

    rotate_about_axis( &vpp, new_v, &tmp, 'z', POINT_SET_ALL );
}

double randn( double mu, double sigma )
/* A function generate a number drawn from a normal distribution.
 * Obtained from https://phoxis.org/2013/05/04/generating-random-numbers-
 * from-normal-distribution-in-c/.
 *
 * This function assumes that the random number generator has already been
 * seeded.
 */
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double) X2);
    }

    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);

    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (double) X1);
}
