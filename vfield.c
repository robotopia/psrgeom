/*****************************************************************************
 * A source file for the library "psrgeom", which includes structs and
 * functions for treating pulsar geometry.
 *
 * Author: Sam McSweeney
 * Date  : 2017
 *
 * Description:
 *   This source file implements the functions for calculating the magnetic
 *   field.
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
    double y     = X->x[1];
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

    double Om    = psr->Om.rad;

    // Taking into account finite light travel time...
    // Set up so that pole is where we think it is
    psr_angle phase;
    set_psr_angle_rad( &phase, (r - psr->r)/psr->rL );

    double a[] =
    {
        1.0 / r5,
        psr->al.sin * phase.cos / r3 / psr->rL2,
        psr->al.sin * ( phase.cos/r5 + phase.sin/r4/psr->rL),
        -psr->al.sin * phase.sin / r3 / psr->rL2,
        psr->al.sin * (-phase.sin/r5 + phase.cos/r4/psr->rL)
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

        B1->r = Blen;
    }

    // Return, if there's nothing more to calculate
    if (!calcV && !calcA)
        return;

    /***********************************************************
     * Both V and A need the calculation of V, so do this now. * 
     ***********************************************************/

    double pdBn  = -y*Bn[0] + x*Bn[1];   // ph dot Bn
    double Om2   = Om*Om;
    double pdBn2 = pdBn * pdBn;          // (ph dot Bn)^2
    double rho   = sqrt( x*x + y*y );
    double rho2  = rho*rho;
    double det   = Om2*pdBn2 - (rho2*Om2 - v*v);
    //      ^-- This is the bit under the sqrt sign

    // See how many V solutions we get
    if (det < 0.0)
    {
        *nsols = 0;
        return;
    }
    else if (det == 1)
    {
        *nsols = 1;
        calcA = 0; /* If there's only one V solution, then the A solution
                      corresponds to A = inf, which we don't care about */
    }
    else // (det >  0.0)
    {
        *nsols = 2;
    }

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
            V1->r    = Vpos_len;
        }
        if (V2)
        {
            V2->x[i] = Vneg[i] / Vneg_len;
            V2->r    = Vneg_len;
        }
    }

    // Return, if there's nothing more to calculate
    if (!calcA)
        return;

    /********************************************
     * Now, all that remains is to calculate A. * 
     ********************************************/

    // Temporary values for calculating the derivates of a[] above.
    double a_temp[] =
    {
        -5.0 / r7,
        -3.0 * psr->al.sin * phase.cos / r5 / psr->rL2,
        -psr->al.sin / r6 * ( 5.0*phase.cos/r + 4.0*phase.sin/psr->rL ),
        3.0 * psr->al.sin * phase.sin / r5 / psr->rL2,
        -psr->al.sin / r6 * ( 4.0*phase.cos/psr->rL - 5.0*phase.sin/r )
    };

    // The full derivatives of (a[]) with respect to (x,y,z)
    double da[5][3]; // da[4][2] means d(a[4])/dz
    for (i = 0; i < 5; i++)
    for (j = 0; j < 3; j++)
        da[i][j] = a_temp[i]*X->x[j];

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
            dB[i][j] = (da[k][j] * c[k][i]) + (a[k] * dc[k][i][j]);
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
        -Om*dpdBn[0] + chi*(pdBn*dpdBn[0] - 2.0*x), // d(VB)/dx
        -Om*dpdBn[1] + chi*(pdBn*dpdBn[1] - 2.0*y), // d(VB)/dy
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
            A1->r    = Apos_len;
        }
        if (A2)
        {
            A2->x[i] = Aneg[i] / Aneg_len;
            A2->r    = Aneg_len;
        }
    }
}


