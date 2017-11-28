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

void set_pulsar( pulsar *psr, angle *ra, angle *dec, double P, double r,
        angle *al, angle *ze )
{
    // Only set non-NULL angle arguments
    if (ra)
        copy_angle( ra,  &psr->ra  );
    if (dec)
        copy_angle( dec, &psr->dec );

    // Set the pulsar period and the light cylinder radius
    set_pulsar_period( psr, P );

    // Set the pulsar radius
    psr->r = r;

    if (al)
        copy_angle( al, &psr->al );
    if (ze)
        copy_angle( ze, &psr->ze );
}

pulsar *create_pulsar( angle *ra, angle *dec, double P, double r,
        angle *al, angle *ze )
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
{
    psr->P   = P;
    psr->rL  = light_cylinder( P );
    psr->rL2 = psr->rL * psr->rL;
}
