/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2018
 *
 * Description:
 *   This source file implements the functions required to find the emission
 *   locations for a given rotation phase. It uses the NEWUOA algorithm to
 *   converge on the desired point. The cost function incorporates both how
 *   close a point's field line is to the last open field line, and how
 *   close the emission direction is to the line of sight.
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
    pn->freq = freq;
    copy_point( pt, &pn->source );
    calc_fields( pt, psr, SPEED_OF_LIGHT, &pn->B, &pn->V, NULL, &pn->A, NULL,
            NULL );
}
