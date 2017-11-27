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
