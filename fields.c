/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the functions for calculating the magnetic
 *   field, the velocity field, and the acceleration field, as well as a
 *   handful of related quantities.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psrgeom.h"


void calc_fields( point *X, pulsar *psr, double v,
                  point *B1,
                  point *V1, point *V2,
                  point *A1, point *A2,
                  int *nsols )
/* This is a wrapper function for:
 *     calc_dipole_fields(), and
 *     calc_deutsch_fields().
 * It merely tests what field type has been assigned to the pulsar (psr)
 * and calls the appropriate function.
 */
{
    switch (psr->field_type)
    {
        case DIPOLE:
            calc_dipole_fields( X, psr, v, B1, V1, V2, A1, A2, nsols );
            break;
        case DEUTSCH:
            calc_deutsch_fields( X, psr, v, B1, V1, V2, A1, A2, nsols );
            break;
        default:
            fprintf( stderr, "error: calc_fields: unrecognised field type\n" );
            break;
    }
}

void calc_deutsch_fields( point *X, pulsar *psr, double v,
                  point *B1,
                  point *V1, point *V2,
                  point *A1, point *A2,
                  int *nsols )
/* This returns the normalised vectors corresponding to the magnetic field, B;
 * the velocity field, V; and the acceleration field, A. These are the fields
 * corresponding to a particle instantaneously at point X (observer frame)
 * when the magnetic axis is instantaneously in the xz-plane.
 *
 * There are two possible solutions for the velocity and acceleration, denoted
 * V1, V2, etc. For the velocity, it is possible for there to be exactly one
 * solution, in which case V1 will be set to it; however, in ths case, the
 * acceleration blows up, so the acceleration can only have 0 or 2 solutions.
 *
 * If any of B1, V1, V2, A1, or A2 are set to NULL, they are not recorded.
 * Moreover, this function tries to limit the calculation actually performed
 * to only what is needed to calculate the requested fields.
 *
 * If the pulsar is spinning in the negative direction, the strategy is to
 * flip the y-coordinate so that all the calculations are done in the "mirror"
 * universe where the pulsar is spinning in the positive direction, and then
 * flip back the y-coordinates of the derived B, V, and A vectors at the end
 * of the calculations.
 *
 * This function assumes that both the Cartesian and spherical coordinates of
 * X are set.
 *
 * Inputs:
 *   point  *X      = the location of the particle in observer coordinates
 *   pulsar *psr    = the pulsar struct
 *   double  v      = the speed of the particle in the observer frame
 *
 * Outputs:
 *   int    *nsols  = the number of solutions found (0, 1, or 2)
 *   point  *B1     = the normalised mag. field   of the particle
 *   point  *V1     = the normalised velocity     of the particle (sol #1)
 *   point  *V2     = the normalised velocity     of the particle (sol #2)
 *   point  *A1     = the normalised acceleration of the particle (sol #1)
 *   point  *A2     = the normalised acceleration of the particle (sol #2)
 */
{
    // First, figure out how much of the following we need to calculate
    int calcB = 0;
    int calcV = 0;
    int calcA = 0;

    if (B1)         calcB = 1;
    if (V1 || V2)   calcV = 1;
    if (A1 || A2)   calcA = 1;

    if (!calcB && !calcV && !calcA)
        return;

    // Generic loop counters
    int i, j, k;

    /**********************************************************************
     * Everything needs the calculation of B, so do this bit in any case. *
     **********************************************************************/

    // Paul Arendt's equations
    double x     = X->x[0];
    double y     = X->x[1] * psr->spin;
    double z     = X->x[2];

    double xx    = x*x;
    double yy    = y*y;
    double zz    = z*z;

    double xy    = x*y;
    double xz    = x*z;
    double yz    = y*z;

    double rr    = xx + yy + zz;
    double r     = sqrt( rr );
    double r3    = r*rr;
    double r4    = r*r3;
    double r5    = r*r4;    double r5i   = 1.0/r5;
    double r6    = r*r5;
    double r7    = r*r6;

    double rL    = psr->rL;
    double rL2   = psr->rL2;
    double rL3   = rL*rL2;

    double Om    = psr->Om.rad;

    // Taking into account finite light travel time...
    // Set up so that pole is where we think it is
    psr_angle phase;
    set_psr_angle_rad( &phase, (r - psr->r)/psr->rL );
    double phS = phase.sin;
    double phC = phase.cos;

    double a[] =
    {
        psr->al.cos * r5i,
        psr->al.sin * phC/(r3*rL2),
        psr->al.sin * ( phC*r5i + phS/(r4*rL)),
       -psr->al.sin * phS/(r3*rL2),
        psr->al.sin * (-phS*r5i + phC/(r4*rL))
    };

    // The following are the x,y,z terms that get multiplied to the a[] terms
    // to produce the magnetic field (see below)
    double c[5][3] =
    {
        {      3.0*xz,      3.0*yz, 3.0*zz - rr },
        {     rr - xx,         -xy,         -xz },
        { 3.0*xx - rr,      3.0*xy,      3.0*xz },
        {         -xy,     rr - yy,         -yz },
        {      3.0*xy, 3.0*yy - rr,      3.0*yz }
    };

    // Now we construct the magnetic field
    //   B[0] = a[0]*c[0][0] + a[1]*c[1][0] + ...
    double B[3];
    for (i = 0; i < 3; i++)
    {
        B[i] = 0.0;
        for (j = 0; j < 5; j++)
            B[i] += a[j]*c[j][i];
    }

    // Normalise B (--> Bn)
    double Blen = sqrt( B[0]*B[0] +
                        B[1]*B[1] +
                        B[2]*B[2] );
    double Bn[3];
    for (i = 0; i < 3; i++)
        Bn[i] = B[i] / Blen;

    // If requested, record B result:
    if (calcB)
    {
        for (i = 0; i < 3; i++)
            B1->x[i] = Bn[i];

        B1->x[1] *= psr->spin;
        B1->r = Blen;
    }

    // Return, if there's nothing more to calculate
    if (!calcV && !calcA)
        return;

    /***********************************************************
     * Both V and A need the calculation of V, so do this now. * 
     ***********************************************************/

    /*
     *     V_B = -Ω(φ∙B̂) ± √[Ω²(φ∙B̂)² - (ρ²Ω² - v²)]
     *
     *     V = Ωφ + (V_B)B̂
     */

    double pdBn  = -y*Bn[0] + x*Bn[1];   // ph dot Bn
    double Om2   = Om*Om;
    double pdBn2 = pdBn * pdBn;          // (ph dot Bn)^2
    double rho2  = x*x + y*y;
    double det   = Om2*pdBn2 - (rho2*Om2 - v*v);
    //      ^-- This is the bit under the sqrt sign

    // See how many V solutions we get
    int nsol; // Use this as a place holder for nsols for now
    if (det < 0.0)
    {
        nsol = 0;
        return;
    }
    else if (det == 1.0)
    {
        nsol = 1;
        calcA = 0; /* If there's only one V solution, then the A solution
                      corresponds to A = inf, which we don't care about */
    }
    else // (det >  0.0)
    {
        nsol = 2;
    }

    // If the caller wants to know nsols, they can request it
    if (nsols)
        *nsols = nsol;

    double sqrt_det = sqrt(det);
    double VBpos    = -Om*pdBn + sqrt_det;
    double VBneg    = -Om*pdBn - sqrt_det;

    // Calulate the velocity fields
    double Vpos[3] = { -y*Om + VBpos*Bn[0], x*Om + VBpos*Bn[1], VBpos*Bn[2] };
    double Vneg[3] = { -y*Om + VBneg*Bn[0], x*Om + VBneg*Bn[1], VBneg*Bn[2] };

    // The length of this velocity should be v, but calculate it explicitly
    // just to make sure we're doing it correctly.
    double Vpos_len = sqrt( Vpos[0] * Vpos[0] +
                            Vpos[1] * Vpos[1] +
                            Vpos[2] * Vpos[2] );
    double Vneg_len = sqrt( Vneg[0] * Vneg[0] +
                            Vneg[1] * Vneg[1] +
                            Vneg[2] * Vneg[2] );

    // If requested, record V results:
    for (i = 0; i < 3; i++)
    {
        if (V1)
        {
            V1->x[i] = Vpos[i] / Vpos_len;
            if (i == 1) V1->x[i] *= psr->spin;
        }
        if (V2)
        {
            V2->x[i] = Vneg[i] / Vneg_len;
            if (i == 1) V2->x[i] *= psr->spin;
        }
    }
    if (V1)  V1->r = Vpos_len;
    if (V2)  V2->r = Vneg_len;

    // Return, if there's nothing more to calculate
    if (!calcA)
        return;

    /********************************************
     * Now, all that remains is to calculate A. * 
     ********************************************/

    // Temporary values for calculating the derivates of a[] above.
    double r7rL2i = 1.0 / (r7 * rL2);
    double r5rL3i = 1.0 / (r5 * rL3);
    double a_temp[] =
    {
        -5.0*psr->al.cos / r7,
        -r5rL3i * psr->al.sin * (r*phS + 3.0*rL*phC),
        -r7rL2i * psr->al.sin * (5.0*r*rL*phS + (5.0*rL2 - rr)*phC),
        -r5rL3i * psr->al.sin * (r*phC - 3.0*rL*phS),
        -r7rL2i * psr->al.sin * (5.0*r*rL*phC - (5.0*rL2 - rr)*phS)
    };

    // The full derivatives of (a[]) with respect to (x,y,z)
    double da[5][3]; // da[4][2] means d(a[4])/dz
    for (i = 0; i < 5; i++)
    for (j = 0; j < 3; j++)
    {
        da[i][j] = a_temp[i]*X->x[j];
        if (j == 1) // The y-coordinate that may need to be flipped
            da[i][j] *= psr->spin;
    }

    // The derivatives of c[][] above, with respect to (x,y,z):
    // e.g. dc[4][2][1] means d(c[4][2])/dy
    double dc[5][3][3] =
    {
        {
            {  3.0*z,    0.0, 3.0*x },
            {    0.0,  3.0*z, 3.0*y },
            { -2.0*x, -2.0*y, 4.0*z }
        },
        {
            { 0.0, 2.0*y, 2.0*z },
            {  -y,    -x,   0.0 },
            {  -z,   0.0,    -x }
        },
        {
            { 4.0*x, -2.0*y, -2.0*z },
            { 3.0*y,  3.0*x,    0.0 },
            { 3.0*z,    0.0,  3.0*x }
        },
        {
            {    -y,  -x,   0.0 },
            { 2.0*x, 0.0, 2.0*z },
            {   0.0,  -z,    -y }
        },
        {
            {  3.0*y, 3.0*x,    0.0 },
            { -2.0*x, 4.0*y, -2.0*z },
            {    0.0, 3.0*z,  3.0*y }
        }
    };

    // And now the partial derivatives of the magnetic field with respect to
    // (x,y,z)
    double dB[3][4]; /* dB[0][1] means d(B_x)/dy; dB[][3] means d/dt, which is
                                                  calculated later. */
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
    {
        dB[i][j] = 0.0;
        for (k = 0; k < 5; k++)
            dB[i][j] += (da[k][j] * c[k][i]) + (a[k] * dc[k][i][j]);
    }

    // There's one more partial derivative we need: dB/dt
    for (i = 0; i < 3; i++)
        dB[i][3] = Om * (-y*dB[i][0] + x*dB[i][1]);

    // Now calculate how Bn changes with respect to (x,y,z,t)
    double dBn[3][4];
    double Bn_dot_dB;
    for (i = 0; i < 4; i++)
    {
        Bn_dot_dB = 0.0;
        for (j = 0; j < 3; j++)
            Bn_dot_dB += Bn[j] * dB[j][i];

        for (j = 0; j < 3; j++)
            dBn[j][i] = (dB[j][i] - Bn_dot_dB * Bn[j]) / Blen;
    }

    double dpdBn[4] =
    {
        -y*dBn[0][0] + x*dBn[1][0] + Bn[1],  // d(ph.Bn)/dx
        -y*dBn[0][1] + x*dBn[1][1] - Bn[0],  // d(ph.Bn)/dy
        -y*dBn[0][2] + x*dBn[1][2],          // d(ph.Bn)/dz
        -y*dBn[0][3] + x*dBn[1][3]           // d(ph.Bn)/dt
    };

    double chi = Om2 / sqrt_det;

    double dVBpos[4] =
    {
        -Om*dpdBn[0] + chi*(pdBn*dpdBn[0] - x), // d(VB)/dx
        -Om*dpdBn[1] + chi*(pdBn*dpdBn[1] - y), // d(VB)/dy
        -Om*dpdBn[2] + chi*(pdBn*dpdBn[2]),         // d(VB)/dz
        -Om*dpdBn[3] + chi*(pdBn*dpdBn[3])          // d(VB)/dt
    };

    double dVBneg[4] =
    {
        -Om*dpdBn[0] - chi*(pdBn*dpdBn[0] - 2.0*x), // d(VB)/dx
        -Om*dpdBn[1] - chi*(pdBn*dpdBn[1] - 2.0*y), // d(VB)/dy
        -Om*dpdBn[2] - chi*(pdBn*dpdBn[2]),         // d(VB)/dz
        -Om*dpdBn[3] - chi*(pdBn*dpdBn[3])          // d(VB)/dt
    };

    double dVpos[3][4], dVneg[3][4];
    for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
    {
        dVpos[i][j] = VBpos*dBn[i][j] + dVBpos[j]*Bn[i];
        dVneg[i][j] = VBneg*dBn[i][j] + dVBneg[j]*Bn[i];
    }
    dVpos[1][0] += Om;
    dVneg[1][0] += Om;
    dVpos[0][1] -= Om;
    dVneg[0][1] -= Om;

    // Calculate the acceleration fields
    double Apos[3], Aneg[3];
    for (i = 0; i < 3; i++)
    {
        Apos[i] = Vpos[0]*dVpos[i][0] +
                  Vpos[1]*dVpos[i][1] +
                  Vpos[2]*dVpos[i][2] +
                          dVpos[i][3];
        Aneg[i] = Vneg[0]*dVneg[i][0] +
                  Vneg[1]*dVneg[i][1] +
                  Vneg[2]*dVneg[i][2] +
                          dVneg[i][3];
    };

    // Calculate the magnitudes of the acceleration vectors.
    double Apos_len = sqrt( Apos[0] * Apos[0] +
                            Apos[1] * Apos[1] +
                            Apos[2] * Apos[2] );
    double Aneg_len = sqrt( Aneg[0] * Aneg[0] +
                            Aneg[1] * Aneg[1] +
                            Aneg[2] * Aneg[2] );

    // If requested, record A results:
    for (i = 0; i < 3; i++)
    {
        if (A1)
        {
            A1->x[i] = Apos[i] / Apos_len;
            if (i == 1) A1->x[i] *= psr->spin;
        }
        if (A2)
        {
            A2->x[i] = Aneg[i] / Aneg_len;
            if (i == 1) A2->x[i] *= psr->spin;
        }
    }

    if (A1)  A1->r = Apos_len;
    if (A2)  A2->r = Aneg_len;
}



void calc_dipole_fields( point *X, pulsar *psr, double v,
                  point *B1,
                  point *V1, point *V2,
                  point *A1, point *A2,
                  int *nsols )
/* This function is identical to calc_fields(), only instead of returning the
 * Deutsch solution, it returns a simple dipole field. For an explanation of
 * this function's inputs and outputs, refer to calc_fields().
 */
{
    // First, figure out how much of the following we need to calculate
    int calcB = 0;
    int calcV = 0;
    int calcA = 0;

    if (B1)         calcB = 1;
    if (V1 || V2)   calcV = 1;
    if (A1 || A2)   calcA = 1;

    if (!calcB && !calcV && !calcA)
        return;

    // Generic loop counters
    int i, j;

    /**********************************************************************
     * Everything needs the calculation of B, so do this bit in any case. *
     **********************************************************************/

    double x     = X->x[0];
    double y     = X->x[1] * psr->spin;
    double z     = X->x[2];

    double xx    = x*x;
    double yy    = y*y;
    double zz    = z*z;

    double xy    = x*y;
    double xz    = x*z;
    double yz    = y*z;

    double rr    = xx + yy + zz;
    double r     = sqrt( rr );
    double r3    = r*rr;
    double r4    = r*r3;
    double r5    = r*r4;
    double r6    = r*r5;
    double r7    = r*r6;

    double r5i   = 1.0/r5;
    double r3i   = r5i*rr;

    double Om    = psr->Om.rad;

    double B[3];
    double temp_sin = 3.0*psr->al.sin*r5i;
    double temp_cos = 3.0*psr->al.cos*r5i;
    B[0] = xz*temp_cos - psr->al.sin*r3i + xx*temp_sin;
    B[1] = yz*temp_cos                   + xy*temp_sin;
    B[2] = zz*temp_cos - psr->al.cos*r3i + xz*temp_sin;

    // Normalise B (--> Bn)
    double Blen = sqrt( B[0]*B[0] +
                        B[1]*B[1] +
                        B[2]*B[2] );
    double Bn[3];
    for (i = 0; i < 3; i++)
        Bn[i] = B[i] / Blen;

    // If requested, record B result:
    if (calcB)
    {
        for (i = 0; i < 3; i++)
            B1->x[i] = Bn[i];

        B1->x[1] *= psr->spin;
        B1->r = Blen;
    }

    // Return, if there's nothing more to calculate
    if (!calcV && !calcA)
        return;

    /***********************************************************
     * Both V and A need the calculation of V, so do this now. * 
     ***********************************************************/

    /*
     *     V_B = -Ω(φ∙B̂) ± √[Ω²(φ∙B̂)² - (ρ²Ω² - v²)]
     *
     *     V = Ωφ + (V_B)B̂
     */

    double pdBn  = -y*Bn[0] + x*Bn[1];   // ph dot Bn
    double Om2   = Om*Om;
    double pdBn2 = pdBn * pdBn;          // (ph dot Bn)^2
    double rho2  = x*x + y*y;
    double det   = Om2*pdBn2 - (rho2*Om2 - v*v);
    //      ^-- This is the bit under the sqrt sign

    // See how many V solutions we get
    int nsol; // Use this as a place holder for nsols for now
    if (det < 0.0)
    {
        nsol = 0;
        return;
    }
    else if (det == 1.0)
    {
        nsol = 1;
        calcA = 0; /* If there's only one V solution, then the A solution
                      corresponds to A = inf, which we don't care about */
    }
    else // (det >  0.0)
    {
        nsol = 2;
    }

    // If the caller wants to know nsols, they can request it
    if (nsols)
        *nsols = nsol;

    double sqrt_det = sqrt(det);
    double VBpos    = -Om*pdBn + sqrt_det;
    double VBneg    = -Om*pdBn - sqrt_det;

    // Calulate the velocity fields
    double Vpos[3] = { -y*Om + VBpos*Bn[0], x*Om + VBpos*Bn[1], VBpos*Bn[2] };
    double Vneg[3] = { -y*Om + VBneg*Bn[0], x*Om + VBneg*Bn[1], VBneg*Bn[2] };

    // The length of this velocity should be v, but calculate it explicitly
    // just to make sure we're doing it correctly.
    double Vpos_len = sqrt( Vpos[0] * Vpos[0] +
                            Vpos[1] * Vpos[1] +
                            Vpos[2] * Vpos[2] );
    double Vneg_len = sqrt( Vneg[0] * Vneg[0] +
                            Vneg[1] * Vneg[1] +
                            Vneg[2] * Vneg[2] );

    // If requested, record V results:
    for (i = 0; i < 3; i++)
    {
        if (V1)
        {
            V1->x[i] = Vpos[i] / Vpos_len;
            if (i == 1) V1->x[i] *= psr->spin;
        }
        if (V2)
        {
            V2->x[i] = Vneg[i] / Vneg_len;
            if (i == 1) V2->x[i] *= psr->spin;
        }
    }
    if (V1)  V1->r = Vpos_len;
    if (V2)  V2->r = Vneg_len;

    // Return, if there's nothing more to calculate
    if (!calcA)
        return;

    /********************************************
     * Now, all that remains is to calculate A. * 
     ********************************************/

    // Temporary values for calculating the derivates of a[] above.
    double r7i = 1.0 / r7;

    // And now the partial derivatives of the magnetic field with respect to
    // (x,y,z)
    double dB[3][4]; /* dB[0][1] means d(B_x)/dy; dB[][3] means d/dt, which is
                                                  calculated later. */
    dB[0][0] = 3.0*r7i*(
                 z*psr->al.cos*(    rr - 5.0*xx) +
                 x*psr->al.sin*(3.0*rr - 5.0*xx)
               );
    dB[0][1] = 3.0*r7i*y*(
                 psr->al.cos*(   - 5.0*xz) +
                 psr->al.sin*(rr - 5.0*xx) 
               );
    dB[0][2] = 3.0*r7i*(
                 x*psr->al.cos*(rr - 5.0*zz) +
                 z*psr->al.sin*(rr - 5.0*xx)
               );
    dB[1][0] = dB[0][1];
    dB[1][1] = 3.0*r7i*(
                 z*psr->al.cos*(rr - 5.0*yy) +
                 x*psr->al.sin*(rr - 5.0*yy)
               );
    dB[1][2] = 3.0*r7i*y*(
                 psr->al.cos*(rr - 5.0*zz) +
                 psr->al.sin*(   - 5.0*xz) 
               );
    dB[2][0] = dB[0][2];
    dB[2][1] = dB[1][2];
    dB[2][2] = 3.0*r7i*(
                 z*psr->al.cos*(3.0*rr - 5.0*zz) +
                 x*psr->al.sin*(    rr - 5.0*zz)
               );

    // There's one more partial derivative we need: dB/dt
    for (i = 0; i < 3; i++)
        dB[i][3] = Om * (-y*dB[i][0] + x*dB[i][1]);

    // Now calculate how Bn changes with respect to (x,y,z,t)
    double dBn[3][4];
    double Bn_dot_dB;
    for (i = 0; i < 4; i++)
    {
        Bn_dot_dB = 0.0;
        for (j = 0; j < 3; j++)
            Bn_dot_dB += Bn[j] * dB[j][i];

        for (j = 0; j < 3; j++)
            dBn[j][i] = (dB[j][i] - Bn_dot_dB * Bn[j]) / Blen;
    }

    double dpdBn[4] =
    {
        -y*dBn[0][0] + x*dBn[1][0] + Bn[1],  // d(ph.Bn)/dx
        -y*dBn[0][1] + x*dBn[1][1] - Bn[0],  // d(ph.Bn)/dy
        -y*dBn[0][2] + x*dBn[1][2],          // d(ph.Bn)/dz
        -y*dBn[0][3] + x*dBn[1][3]           // d(ph.Bn)/dt
    };

    double chi = Om2 / sqrt_det;

    double dVBpos[4] =
    {
        -Om*dpdBn[0] + chi*(pdBn*dpdBn[0] - x), // d(VB)/dx
        -Om*dpdBn[1] + chi*(pdBn*dpdBn[1] - y), // d(VB)/dy
        -Om*dpdBn[2] + chi*(pdBn*dpdBn[2]),         // d(VB)/dz
        -Om*dpdBn[3] + chi*(pdBn*dpdBn[3])          // d(VB)/dt
    };

    double dVBneg[4] =
    {
        -Om*dpdBn[0] - chi*(pdBn*dpdBn[0] - 2.0*x), // d(VB)/dx
        -Om*dpdBn[1] - chi*(pdBn*dpdBn[1] - 2.0*y), // d(VB)/dy
        -Om*dpdBn[2] - chi*(pdBn*dpdBn[2]),         // d(VB)/dz
        -Om*dpdBn[3] - chi*(pdBn*dpdBn[3])          // d(VB)/dt
    };

    double dVpos[3][4], dVneg[3][4];
    for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
    {
        dVpos[i][j] = VBpos*dBn[i][j] + dVBpos[j]*Bn[i];
        dVneg[i][j] = VBneg*dBn[i][j] + dVBneg[j]*Bn[i];
    }
    dVpos[1][0] += Om;
    dVneg[1][0] += Om;
    dVpos[0][1] -= Om;
    dVneg[0][1] -= Om;

    // Calculate the acceleration fields
    double Apos[3], Aneg[3];
    for (i = 0; i < 3; i++)
    {
        Apos[i] = Vpos[0]*dVpos[i][0] +
                  Vpos[1]*dVpos[i][1] +
                  Vpos[2]*dVpos[i][2] +
                          dVpos[i][3];
        Aneg[i] = Vneg[0]*dVneg[i][0] +
                  Vneg[1]*dVneg[i][1] +
                  Vneg[2]*dVneg[i][2] +
                          dVneg[i][3];
    };

    // Calculate the magnitudes of the acceleration vectors.
    double Apos_len = sqrt( Apos[0] * Apos[0] +
                            Apos[1] * Apos[1] +
                            Apos[2] * Apos[2] );
    double Aneg_len = sqrt( Aneg[0] * Aneg[0] +
                            Aneg[1] * Aneg[1] +
                            Aneg[2] * Aneg[2] );

    // If requested, record A results:
    for (i = 0; i < 3; i++)
    {
        if (A1)
        {
            A1->x[i] = Apos[i] / Apos_len;
            if (i == 1) A1->x[i] *= psr->spin;
        }
        if (A2)
        {
            A2->x[i] = Aneg[i] / Aneg_len;
            if (i == 1) A2->x[i] *= psr->spin;
        }
    }

    if (A1)  A1->r = Apos_len;
    if (A2)  A2->r = Aneg_len;
}


double Bdotrxy( point *x, pulsar *psr )
/* Returns the dot product of B, the normalised magnetic field, and r_{xy},
 * the vector pointing cylindrically outward, at point *x.
 * This quantity is useful because the last open magnetic field lines have
 * this quantity equal to zero at the light cylinder.
 */
{
    point B;
    calc_fields( x, psr, 0.0, &B, NULL, NULL, NULL, NULL, NULL );
    double retval = B.x[0] * x->ph.cos +
                    B.x[1] * x->ph.sin;
    return retval;
}



int footpoint( point *start_pt, pulsar *psr, double tmult, int direction,
               FILE *write_xyz, int rL_norm, double rL_lim, point *foot_pt )
/* This uses Runge-Kutta (RK4) to follow the magnetic field line outward from
 * an initial starting point "start_pt"
 *
 * Inputs:
 *   start_pt  : the initial starting point
 *   psr       : the pulsar (includes geometric information)
 *   tmult     : the rate at which to progress in RK4 (fraction of psr radius)
 *   direction : progress inwards or outwards along field lines
 *               {DIR_INWARD, DIR_OUTWARD}
 *   write_xyz : file handle where to write x,y,z values
 *               (if NULL, no output is written)
 *   rL_norm   : (bool) normalise written-out points to lt cyl radius
 *   rL_lim    : the rho limit, normalised to the light cylinder radius.
 *               If the points go beyond this limit, stop.
 * Outputs:
 *   foot_pt   : the final point reached during RK algorithm
 * Return values:
 *   STOP_FOUND  : The point was successfully found
 *   STOP_EXCEED : The outer limit (rL_lim) was reached
 */
{
    int retval = STOP_FOUND;

    point x, old_x;

    double tstep;

    copy_point( start_pt, &x );
    tstep = tmult * x.r;

    // Trace this line outwards with a 4 stage Runge-Kutta

    int temp_drctn = direction; // In case we've gone too far
    double precision = 1.0e-14;
    double xscale = (rL_norm ? 1.0/psr->rL : 1.0);
    double rho_lim = rL_lim * psr->rL;
    double rhosq_lim = rho_lim * rho_lim;

    while (1) // Indefinite loop until surface is reached
    {
        // Keep track of the previous point
        copy_point( &x, &old_x );

        // If the point has strayed too far afield, stop
        if (x.rhosq > rhosq_lim)
        {
            retval = STOP_EXCEED;
            break;
        }

        // Write out the current xyz position, if requested,
        // but only if we're still moving "forward"
        if (write_xyz && (temp_drctn == direction))
            fprintf( write_xyz, "%.15e %.15e %.15e %.15e\n",
                     xscale*x.x[0], xscale*x.x[1], xscale*x.x[2], tstep );

        // Take a single RK4 step along the magnetic field
        Bstep( &old_x, psr, tstep, temp_drctn, &x );

        // Recalculate the various distances to the new point
        set_point_xyz( &x, x.x[0], x.x[1], x.x[2],
                POINT_SET_SPH | POINT_SET_RHOSQ );

        /* Figure out whether to stop or not */

        // Error checking: this algorithm should have stopped long before
        // tstep reaches underflow
        if ((tstep/2.0 <= 0.0) || (tstep/2.0 >= tstep))
        {
            fprintf( stderr, "error: footpoint: tstep underflow\n" );
            exit(EXIT_FAILURE);
        }

        // Just check to see if we've happened to land exactly on the surface
        if (x.r == psr->r)
            break;

        // Otherwise, only do anything special if we've crossed the surface
        if ((x.r - psr->r) / (old_x.r - psr->r) < 0.0) // then x and old_x are
                                                       // on opposite sides of
                                                       // the surface
        {
            if ((fabs(x.r - old_x.r) <= precision*psr->r))
                break;

            if (temp_drctn == DIR_OUTWARD)
                temp_drctn = DIR_INWARD;
            else if (temp_drctn == DIR_INWARD)
                temp_drctn = DIR_OUTWARD;
            else
            {
                fprintf( stderr, "error: footpoint: unknown direction (%d)\n",
                                 direction );
                exit(EXIT_FAILURE);
            }
            tstep /= 2.0;
        }

        // Adjust tstep proportionally to how far away from the pulsar we are
        tstep *= x.r / old_x.r;
    }

    // Make the final point available to the caller
    if (foot_pt != NULL)
        copy_point( &x, foot_pt );

    return retval;
}


int cmp_extreme( point *x, pulsar *psr, double precision )
/* This function attempts to compare the point x with the location of the
 * extreme point of the magnetic field line of pulsar psr that passes through
 * x. If you would have to follow the magnetic field line in the same
 * direction as the magnetic field to get to the extreme from x, then 1 is
 * returned. If you would have to follow the magnetic field line in the
 * opposite direction, then -1 is returned. If you are at the extreme exactly
 * (i.e. within the specified precision) then 0 is returned.
 *
 * The extreme points are identified as those at which the magnetic field is
 * perpendicular to the "cylindrically outward pointing" vector (i.e. where
 * the quantity Bdotrxy() vanishes) and which is pointing generally downward
 * (i.e. B.z is negative).
 *
 * Inputs:
 *   x         : The point in question, whose relationship to the extreme
 *               point is sought.
 *   psr       : The pulsar (includes geometric information).
 *   precision : The precision with which we consider x to equal the extreme
 *               point.
 * Return values:
 *   DIR_OUTWARD : The extreme point is "in front" of x, in the sense that you
 *                 would have to progress along the field line in the same
 *                 direction as the magnetic field to reach it.
 *   DIR_INWARD  : The extreme point is "behind" x, in the sense that you
 *                 would have to progress along the field line in the opposite
 *                 direction as the magnetic field to reach it.
 *   DIR_STOP    : x is an extreme point (i.e. B dot r_{xy} equals 0.0 to
 *                 double precision.
 */
{
    // Get the magnetic field
    point B;
    calc_fields( x, psr, 0.0, &B, NULL, NULL, NULL, NULL, NULL );

    // Calling Bdotrxy() would re-calculate the magnetic field unnecessarily,
    // so duplicate the code here.
    double Bdr = B.x[0] * x->ph.cos +
                 B.x[1] * x->ph.sin;

    /* Figure out which segment of the magnetic field line x is on.
      
       If 0 < alpha < 90 (deg):
          Seg 1 = from origin to the point with maximum z value
          Seg 2 = from max z point to extreme point
          Seg 3 = from extreme point to point with minimum z value
          Seg 4 = from min z point to origin

       If 90 < alpha < 180 (deg):
          Seg 1 = from origin to the point with minimum z value
          Seg 2 = from min z point to extreme point
          Seg 3 = from extreme point to point with maximum z value
          Seg 4 = from max z point to origin

       We'll call the extreme point itself Seg 0.
    */

    int seg;
    if (psr->al.deg <= 90.0)
    {
        if (B.x[2] > 0.0)
        {
            if (x->x[2] > 0.0)          seg = 1;
            else                        seg = 4;
        }
        else
        {
            if (fabs(Bdr) <= precision) seg = 0;
            else if (Bdr > 0.0)         seg = 2;
            else /* (Bdr < 0.0) */      seg = 3;
        }
    }
    else // al > 90.0 deg
    {
        if (B.x[2] < 0.0)
        {
            if (x->x[2] < 0.0)  seg = 1;
            else                seg = 4;
        }
        else
        {
            if (fabs(Bdr) <= precision) seg = 0;
            else if (Bdr > 0.0)         seg = 2;
            else /* (Bdr < 0.0) */      seg = 3;
        }
    }

    // Now sum up
    if (seg == 0)
        return DIR_STOP; // We hit the extreme point!

    if (seg == 1 || seg == 2)
        return DIR_OUTWARD; // The extreme point is "in front" of x

    // seg must be either 3 or 4
    return DIR_INWARD; // The extreme point is "behind" x
}


int farpoint( point *start_pt, pulsar *psr, double tmult,
              FILE *write_xyz, int rL_norm, double rL_lim, point *far_pt )
/* This uses Runge-Kutta (RK4) to follow the magnetic field line outward from
 * an initial starting point "start_pt". It seeks the point (x,y,z) such that
 * x^2+y^2 is maximised.
 *
 * Inputs:
 *   start_pt  : the initial starting point
 *   psr       : the pulsar (includes geometric information)
 *   tmult     : the rate at which to progress in RK4 (fraction of psr radius)
 *   write_xyz : file handle where to write x,y,z values
 *               (if NULL, no output is written)
 *   rL_norm   : (bool) normalise written-out points to lt cyl radius
 *   rL_lim    : the rho limit, normalised to the light cylinder radius.
 *               If the points go beyond this limit, stop.
 * Outputs:
 *   far_pt    : the final point reached during RK algorithm
 * Return values:
 *   STOP_FOUND  : The point was successfully found
 *   STOP_EXCEED : The outer limit (rL_lim) was reached
 */
{
    int retval = STOP_FOUND;
    point x, old_x;

    double tstep;

    copy_point( start_pt, &x );
    tstep = tmult * x.r;

    // Trace this line along the magnetic field line with a 4 stage Runge-
    // Kutta in the specified direction.

    double precision = 1.0e-14;
    int direction, prev_direction, init_direction = DIR_STOP;
    double xscale = (rL_norm ? 1.0/psr->rL : 1.0);
    double rho_lim = rL_lim * psr->rL;
    double rhosq_lim = rho_lim * rho_lim;

    while (1) // Indefinite loop until surface is reached
    {
        // On the first time through the while loop, set up directions
        if (init_direction == DIR_STOP)
        {
            direction      = cmp_extreme( &x, psr, precision );
            init_direction = direction;
        }

        // If the point has strayed too far afield, stop
        if (x.rhosq > rhosq_lim)
        {
            retval = STOP_EXCEED;
            break;
        }

        // Check to see if we've landed on the extreme.
        if (direction == DIR_STOP)
            break;

        // Keep track of the previous point
        copy_point( &x, &old_x );
        prev_direction = direction;

        // Take a single RK4 step along the magnetic field
        Bstep( &old_x, psr, tstep, direction, &x );

        // When tstep gets too small, the step will be below machine precision
        // so old_x will equal x. If this happens assume we're at the extreme
        if (old_x.x[0] == x.x[0] &&
            old_x.x[1] == x.x[1] &&
            old_x.x[2] == x.x[2])
        {
            break;
        }

        // Recalculate the various distances to the new point
        set_point_xyz( &x, x.x[0], x.x[1], x.x[2],
                POINT_SET_SPH | POINT_SET_RHOSQ );

        // Find where the extreme point is in relation to the new x.
        direction = cmp_extreme( &x, psr, precision );

        // Write out the current xyz position, if requested,
        // but only if we're still moving in the original direction
        if (write_xyz && (direction == init_direction))
            fprintf( write_xyz, "%.15e %.15e %.15e %.15e\n",
                     xscale*x.x[0], xscale*x.x[1], xscale*x.x[2], tstep );

        // Error checking: this algorithm should have stopped long before
        // tstep reaches underflow
        if ((tstep/2.0 <= 0.0) || (tstep/2.0 >= tstep))
        {
            fprintf( stderr, "error: farpoint: tstep underflow\n" );
            exit(EXIT_FAILURE);
        }

        // If we've changed direction, then halve the step size
        if (direction != prev_direction)
            tstep /= 2.0;

        // Adjust tstep proportionally to how far away from the pulsar we are
        tstep *= x.r / old_x.r;
    }

    // Make the final point available to the caller
    if (far_pt != NULL)
        copy_point( &x, far_pt );

    return retval;
}


void Bstep( point *x1, pulsar *psr, double tstep, int direction, point *x2 )
/* This follows a magnetic field from a given starting point according to one
 * step of the RK4 algorithm.
 *
 * Inputs:
 *   x1       : the starting point
 *   psr      : the pulsar (includes geometric information)
 *   tstep    : the size of the RK4 step to take
 *   direction: whether to follow the line DIR_INWARD or DIR_OUTWARD
 *
 * Outputs:
 *   x2       : the ending point (only Cartesian coordinates are set)
 */
{
    point slop1, slop2, slop3, slope;
    point xp1, xp2;

    if (direction == DIR_STOP)
    {
        fprintf( stderr, "warning: Bstep: direction set to DIR_STOP\n" );
        copy_point( x1, x2 );
        return;
    }

    double sgn;
    if (direction == DIR_OUTWARD)
        sgn = 1.0;
    else if (direction == DIR_INWARD)
        sgn = -1.0;
    else
    {
        fprintf( stderr, "error: Bstep: unrecognised direction (%d)\n",
                         direction );
        exit(EXIT_FAILURE);
    }

    int i; // Generic loop counter

    // First stage
    calc_fields( x1, psr, 0.0, &slop1, NULL, NULL, NULL, NULL, NULL );
    for (i = 0; i < NDEP; i++)
        xp1.x[i] = x1->x[i] + 0.5*sgn*slop1.x[i]*tstep;

    // Second stage
    calc_fields( &xp1, psr, 0.0, &slop2, NULL, NULL, NULL, NULL, NULL );
    for (i = 0; i < NDEP; i++)
        xp2.x[i] = x1->x[i] + 0.5*sgn*slop2.x[i]*tstep;

    // Third stage
    calc_fields( &xp2, psr, 0.0, &slop3, NULL, NULL, NULL, NULL, NULL );
    for (i = 0; i < NDEP; i++)
        xp1.x[i] = x1->x[i] + sgn*slop3.x[i]*tstep;

    // Last stage
    calc_fields( &xp1, psr, 0.0, &slope, NULL, NULL, NULL, NULL, NULL );
    for (i = 0; i < NDEP; i++)
    {
        x2->x[i] = x1->x[i] +
            tstep*sgn*(    slope.x[i] +     slop1.x[i] +
                       2.0*slop2.x[i] + 2.0*slop3.x[i]) / 6.0;
    }
}


void traj_step( point *x1, double t, pulsar *psr, double tstep, int direction,
            point *x2, int rL_norm, FILE *f )
/* This follows a particle's trajectory from a given starting point to time t.
 * This necessarily cannot be a naive following of the V field, because the V
 * field takes a point to a neighbouring point and the neighbouring point
 * must be interpreted in the rotating frame.
 *
 * What we're really doing, then, is taking steps along the magnetic field,
 * but only after de-rotating the given point back to where it was after time
 * t. Unlike Bstep(), where we interpret tstep as a distance, here we
 * interpret tstep as a time.
 *
 * Inputs:
 *   x1       : the starting point
 *   t        : the time since rotation started
 *   psr      : the pulsar (includes geometric information)
 *   tstep    : the time step for RK4, assuming travelling at light speed
 *   direction: whether to follow the line DIR_INWARD or DIR_OUTWARD
 *   rL_norm  : whether to normalise output coords to light cyl units
 *   f        : where to write out the intermediate steps
 *
 * Outputs:
 *   x2       : the ending point
 */
{
    // De-rotate the given point x1 back time t
    psr_angle back_rot;
    set_psr_angle_rad( &back_rot, -t*psr->Om.rad );
    rotate_about_axis( x1, x1, &back_rot, 'z', POINT_SET_ALL );

    // Evaluate the B and V fields there
    double v = SPEED_OF_LIGHT;
    point B, V1, V2, A1, A2;
    point *V = (direction == DIR_OUTWARD ? &V1 : &V2);
    point *A = (direction == DIR_OUTWARD ? &A1 : &A2);
    calc_fields( x1, psr, v, &B, &V1, &V2, &A1, &A2, NULL );

    // Convert the supplied time step into a size step along B,
    // which is V dotted with B
    double V_B = v * (V->x[0] * B.x[0] +
                      V->x[1] * B.x[1] +
                      V->x[2] * B.x[2]); /* This has units of velocity */
    double dstep = tstep * V_B;

    // dstep is a distance, so now we can just call Bstep() to step along the
    // magnetic field by the appropriate amount.

    Bstep( x1, psr, dstep, direction, x2 );

    // Finally, re-rotate the point back to the correct rotation
    psr_angle rerotate;
    set_psr_angle_rad( &rerotate, (t+tstep)*psr->Om.rad );
    rotate_about_axis( x2, x2, &rerotate, 'z', POINT_SET_ALL );

    // If requested, write out the trajectory vectors (X, V, A)
    if (f != NULL)
    {
        // Remember that we have to re-rotate V and A as well!
        rotate_about_axis( V, V, &rerotate, 'z', POINT_SET_XYZ );
        rotate_about_axis( A, A, &rerotate, 'z', POINT_SET_XYZ );

        // Do we normalise point coordinates to light cylinder units?
        double s = (rL_norm ? 1.0/psr->rL : 1.0);

        fprintf( f, "%.15e  ", t );
        fprintf( f, "%.15e %.15e %.15e  ", s*x2->x[0], s*x2->x[1], s*x2->x[2]);
        fprintf( f, "%.15e %.15e %.15e  ", V->x[0],  V->x[1],  V->x[2]);
        fprintf( f, "%.15e %.15e %.15e\n", A->x[0],  A->x[1],  A->x[2]);
    }
}


int calc_pol_angle( pulsar *psr, psr_angle *phase, int direction,
                    point *init_pt, point *emit_pt, psr_angle *psi )
/* This function calculates the polarisation angle, Ψ, for a pulsar at a
 * given phase, φ.
 *
 * Inputs:
 *   pulsar *psr      : the pulsar geometry
 *   psr_angle *phase : the rotation phase
 *   int direction    : whether the particles are moving along or against the
 *                      magnetic field. Can be DIR_OUTWARD (along) or
 *                      DIR_INWARD (against).
 *   point *init_pt   : an initial guess of where the emission point it
 * Outputs:
 *   point *emit_pt   : the point of emission
 *   psr_angle *psi   : the polsarisation angle. Guaranteed to be in the range
 *                      0° ≤ Ψ < 180°
 */
{
    // First order of business: find the emission point corresponding to this
    // geometry and rotation phase
    int status;
    status = find_emission_point_elevator( psr, phase, direction,
                                           init_pt, emit_pt, NULL );

    if (status == EMIT_PT_TOO_HIGH || status == EMIT_PT_TOO_LOW)
        return 0; // = failed to find point

    // Second, find the acceleration vectors at that point
    point A1, A2;
    point *A = (direction == DIR_OUTWARD ? &A1 : &A2);
    int nsols;
    calc_fields( emit_pt, psr, SPEED_OF_LIGHT, NULL,
                 NULL, NULL, &A1, &A2, &nsols );

    if (nsols <= 0)
    {
        fprintf( stderr, "error: calc_pol_angle: did not find a solution\n" );
        return 0; // = fail
    }

    // Third, convert the acceleration to a polarisation angle
    accel_to_pol_angle( psr, A, phase, psi );

    return 1; // = success
}


void accel_to_pol_angle( pulsar *psr, point *A, psr_angle *phase,
        psr_angle *psi )
/* This function calculates the polarisation angle, Ψ, for a pulsar at a
 * given phase, φ.
 *
 * Inputs:
 *   pulsar *psr      : the pulsar geometry
 *   point *A         : the acceleration vector
 *   psr_angle *phase : the rotation phase
 * Outputs:
 *   psr_angle *psi   : the polsarisation angle. Guaranteed to be in the range
 *                      0° ≤ Ψ < 180°
 */
{
    point u; // (u)nrotated
    psr_angle uz; // For (u)nrotating the (z)eta angle
    reverse_psr_angle( &psr->ze, &uz );
    psr_angle ph;
    if (psr->spin == SPIN_POS)
        copy_psr_angle( phase, &ph );
    else
        reverse_psr_angle( phase, &ph );

    rotate_about_axis(  A, &u, &ph, 'z', POINT_SET_ALL );
    rotate_about_axis( &u, &u, &uz, 'y', POINT_SET_ALL );

    // Now we have a vector whose x and y coords give us the angle we need
    set_psr_angle_rad( psi, -atan( u.x[1] / u.x[0] ) );
}


double calc_curvature( point *V, point *A )
/* Calculate the curvature, given the velocity and acceleration vectors:
 *
 *       |v⃗×a⃗|         |a⃗|
 *   κ = ----- = |v̂×â| ---
 *        |v⃗|³         |v⃗|²
 *
 * V and A are expected to have the normalised components stored in V.x[0],
 * etc, and the lengths to be in V.r, etc.
 */
{
    // Calculate v̂×â (i.e. the cross product of the normalised v⃗ and a⃗)
    point VxA;
    set_point_xyz( &VxA, 
                   V->x[1] * A->x[2]  -  A->x[1] * V->x[2],
                   V->x[2] * A->x[0]  -  A->x[2] * V->x[0],
                   V->x[0] * A->x[1]  -  A->x[0] * V->x[1],
                   POINT_SET_R );

    // Now multiply the length by |a⃗|/|v⃗|² to get the final curvature
    return VxA.r * A->r / (V->r * V->r);
}


