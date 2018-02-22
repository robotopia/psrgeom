#include <math.h>
#include "psrgeom.h"

void move_around_cyl( point *start_pt, psr_angle *xi, double dist,
                      point *end_pt )
/* From a starting point, step around the cylinder a fixed distance at a
 * specified angle. This function assumes that the starting point has
 * defined spherical coordinates.
 *
 * Inputs:
 *   start_pt : the starting point
 *   xi       : the angle the trajectory makes with the xy-plane
 *   dist     : the distance to travel
 *
 * Outputs:
 *   end_pt   : the final point
 */
{
    // "Unwrap" the cylinder
    double rh = start_pt->r * start_pt->th.sin; // = r sin θ
    double ph = start_pt->ph.rad * rh;

    double znew = start_pt->x[2] + dist*xi->sin;
    psr_angle phnew;
    set_psr_angle_rad( &phnew, (ph + dist*xi->cos)/rh );

    set_point_cyl( end_pt, rh, &phnew, znew, POINT_SET_ALL );
}

