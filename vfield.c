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

void quad_solve_real( double a, double b, double c, double *x1, double *x2,
                      int *nsols )
/* Solve a quadratic equation ax^2 + bx + c = 0. If there are no real
 * solutions, neither x1 nor x2 are set. If there is one real solution,
 * it is put into x1.
 *
 * Inputs:
 *   double a, b, c  = quadratic coefficients
 * Outputs:
 *   double *x1, *x2 = real solutions
 *   int *nsols      = the number of real solutions (can be 0, 1, 2)
 */
{
    double b2      = b*b;
    double d       = 0.5/a;
    double det     = b2 - 4.0*a*c;

    if (det < 0.0)
        *nsols = 0;
    else if (det == 0.0)
    {
        *nsols = 1;
        *x1 = b2 * d;
    }
    else // (det > 0.0)
    {
        *nsols = 2;
        double sqr_det = sqrt(det);
        *x1 = (b2 + sqr_det) * d;
        *x2 = (b2 - sqr_det) * d;
    }
}

void Vfield( point *x, pulsar *psr, double v, point *V1, point *V2,
             int *nsols )
/* This returns the normalised vector V (cast as a "point"), i.e. the velocity
 * of a particle instantaneously at point x (observer frame) when the magnetic
 * axis is instantaneously in the xz-plane. It is understood that the
 * azimuthal component of the velocity is fixed by the rotation of the pulsar
 * and that the total velocity of the particle is v.
 *
 * This function assumes that both the Cartesian and spherical coordinates of
 * x are set.
 *
 * Inputs:
 *   point  *x      = the location of the particle in observer coordinates
 *   pulsar *psr    = the pulsar struct
 *   double  v      = the speed of the particle in the observer frame
 *
 * Outputs:
 *   int    *nsols  = the number of solutions found (0, 1, or 2)
 *   point  *V1     = the normalised velocity vector of the particle (sol #1)
 *   point  *V2     = the normalised velocity vector of the particle (sol #2)
 */
{
    // The vectors we will need (cast as points)
    point B;    // The magnetic field
    point az_V; // The azimuthal velocity

    // First, get the magnetic field direction at x
    Bfield( x, psr, &B );

    // Now, calculate the azimuth velocity at x due to the pulsar's rotation
    double rho = x->r * x->th.sin; // The perp dist from the z-axis
    double az_v = rho * psr->Om.rad;
    psr_angle ph_V;
    set_psr_angle_deg( &ph_V, x->ph.deg + 90.0 );
    set_point_cyl( &az_V, 1.0, &ph_V, 0.0, POINT_SET_ALL ); // A unit vector

    /* The total velocity is a weighted sum of the azimuthal velocity (whose
       speed we know) and the velocity along the magnetic field line. Although
       we don't know the magnitude of the latter, we can work it out
       geometrically (SAS triangle). As a quadratic equation, there are
       potentially two solutions, both of which we return, if present.
    */
    double a, b, c;  // The quadratic coefficients
    double vB1, vB2; // The solutions

    // (These are derived from squaring V = V_az + V_B)
    a = 1.0;
    b = az_v * (az_V.x[0] * B.x[0] +
                az_V.x[1] * B.x[1] +
                az_V.x[2] * B.x[2]);
    c = az_v*az_v - v*v;

    quad_solve_real( a, b, c, &vB1, &vB2, nsols );

    // Calculate the final velocity
    int i; // Generic counter
    if (*nsols > 0) // i.e. there is at least one solution
    {
        set_point_xyz( V1, az_v*az_V.x[0] + vB1*B.x[0],
                           az_v*az_V.x[1] + vB1*B.x[1],
                           az_v*az_V.x[2] + vB1*B.x[2], POINT_SET_ALL );
        // Normalise
        for (i = 0; i < NDEP; i++)
            V1->x[i] /= V1->r;
    }
    if (*nsols > 1) // i.e. there are two solutions
    {
        set_point_xyz( V2, az_v*az_V.x[0] + vB2*B.x[0],
                           az_v*az_V.x[1] + vB2*B.x[1],
                           az_v*az_V.x[2] + vB2*B.x[2], POINT_SET_ALL );
        // Normalise
        for (i = 0; i < NDEP; i++)
            V2->x[i] /= V2->r;
    }
}

void Afield( point *X, pulsar *psr, double v, point *A1, point *A2,
             int *nsols )
/* This returns the normalised vector A (cast as a "point"), i.e. the accel-
 * eration of a particle instantaneously at point x (observer frame) when the
 * magnetic axis is instantaneously in the xz-plane.
 *
 * This function assumes that both the Cartesian and spherical coordinates of
 * x are set.
 *
 * Inputs:
 *   point  *X      = the location of the particle in observer coordinates
 *   pulsar *psr    = the pulsar struct
 *   double  v      = the speed of the particle in the observer frame
 *
 * Outputs:
 *   int    *nsols  = the number of solutions found (0, 1, or 2)
 *   point  *A1     = the normalised acceleration of the particle (sol #1)
 *   point  *A2     = the normalised acceleration of the particle (sol #2)
 */
{
    // Generic loop counters
    int i, j, k;

    // This initial bit is very similar to Bfield()
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

    // Paul Arendt's equations
    double a[] =
    {
        1.0 / r5,
        psr->al.sin * phase.cos / r3 / psr->rL2,
        psr->al.sin * ( phase.cos/r5 + phase.sin/r4/psr->rL),
        psr->al.sin * phase.sin / r3 / psr->rL2,
        psr->al.sin * (-phase.sin/r5 + phase.cos/r4/psr->rL)
    };

    // and their derivatives with respect to (x,y,z)
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

    // ... and their derivatives with respect to (x,y,z):
    // dc[4][2][1] means d(c[4][2])/dy
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

    // Now we construct the magnetic field
    //   B[0] = a[0]*c[0][0] + a[1]*c[1][0] + ...
    double B[3];
    for (i = 0; i < 3; i++)
    {
        B[i] = 0.0;
        for (j = 0; j < 5; j++)
            B[i] += a[j]*c[j][i];
    }

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

    // Normalise B (--> Bn)
    double Blen = sqrt( B[0]*B[0] +
                        B[1]*B[1] +
                        B[2]*B[2] );
    double Bn[3];
    for (i = 0; i < 3; i++)
        Bn[i] = B[i] / Blen;

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

    // Now prepare for the calculation of the change of velocity with respect
    // to (x,y,z,t)
    double pdBn  = -y*Bn[0] + x*Bn[1]; // ph dot Bn
    double Om2   = Om*Om;
    double pdBn2 = pdBn * pdBn;
    double rho   = sqrt( x*x + y*y );
    double rho2  = rho*rho;

    double det = Om2 * pdBn2 - 4.0*(rho2*Om2 - v*v);
    // Make sure that we get two valid solutions
    if (det <= 0.0)
    {
        *nsols = 0;
        return;
    }
    else // (det >  0.0)
    {
        *nsols = 2;
    }

    double sqrt_det = sqrt(det);
    double chi = Om2 / (2.0 * sqrt_det);
    double VBpos = 0.5*(-Om*pdBn + sqrt_det);
    double VBneg = 0.5*(-Om*pdBn - sqrt_det);

    double dpdBn[4] =
    {
        -y*dBn[0][0] + x*dBn[1][0] + Bn[1],  // d(ph.Bn)/dx
        -y*dBn[0][1] + x*dBn[1][1] - Bn[0],  // d(ph.Bn)/dy
        -y*dBn[0][2] + x*dBn[1][2],          // d(ph.Bn)/dz
        -y*dBn[0][3] + x*dBn[1][3]           // d(ph.Bn)/dt
    };

    double dVBpos[4] =
    {
        -0.5*Om*dpdBn[0] + chi*(pdBn*dpdBn[0] - 4.0*x), // d(VB)/dx
        -0.5*Om*dpdBn[1] + chi*(pdBn*dpdBn[1] - 4.0*y), // d(VB)/dy
        -0.5*Om*dpdBn[2] + chi*(pdBn*dpdBn[2]),         // d(VB)/dz
        -0.5*Om*dpdBn[3] + chi*(pdBn*dpdBn[3])          // d(VB)/dt
    };

    double dVBneg[4] =
    {
        -0.5*Om*dpdBn[0] - chi*(pdBn*dpdBn[0] - 4.0*x), // d(VB)/dx
        -0.5*Om*dpdBn[1] - chi*(pdBn*dpdBn[1] - 4.0*y), // d(VB)/dy
        -0.5*Om*dpdBn[2] - chi*(pdBn*dpdBn[2]),         // d(VB)/dz
        -0.5*Om*dpdBn[3] - chi*(pdBn*dpdBn[3])          // d(VB)/dt
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

    // Calulate the velocity fields
    double Vpos[3] = { -y*Om + VBpos*Bn[0], x*Om + VBpos*Bn[1], VBpos*Bn[2] };
    double Vneg[3] = { -y*Om + VBneg*Bn[0], x*Om + VBneg*Bn[1], VBneg*Bn[2] };

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

    // The flags to set for vectors
    int F = POINT_SET_XYZ | POINT_SET_R;

    set_point_xyz( A1, Apos[0], Apos[1], Apos[2], F );
    set_point_xyz( A2, Aneg[0], Aneg[1], Aneg[2], F );
}



