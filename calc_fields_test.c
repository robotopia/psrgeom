#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "psrgeom.h"

#define  NTESTS  1000

// define a random double between 0 and x
#define  RANDUNS(x)  ((x)*(double)rand()/(double)RAND_MAX)

// define a random double between -x and x
#define  RANDS(x)  ((2.0*(double)rand()/(double)RAND_MAX - 1.0)*(x))

int main()
{
    srand( time( NULL ) );

    // Declare needed variables for testing
    pulsar psr;
    psr_angle al;
    psr_angle ze;
    double P;
    double r;
    double al_deg;
    double ze_deg;
    double rL;
    point X;
    point B, V1, V2, A1, A2;
    int nsols;
    double rho; // = sqrt(x^2+y^2)
    point Vaz_norm;
    point V_sub_Vaz;
    double v = SPEED_OF_LIGHT;
    double SdB;
    double spin_speed;
    double V_dot_A;

    // do the whole test NTESTS times
    int n;
    int npassed = 0;
    for (n = 0; n < NTESTS; n++)
    {

        // choose a random pulsar geometry
        P      = RANDUNS(1.0) + 0.001; // 0.001 <= P < 1.001
        r      = 1e4;
        al_deg = RANDUNS(180.0);
        ze_deg = RANDUNS(180.0);

        // choose a set pulsar geometry
        /*
        P      = 0.5;
        r      = 1e4;
        al_deg = 30.0;
        ze_deg = 40.0; // <-- Not used in field calculation
        */

        set_psr_angle_deg( &al, al_deg );
        set_psr_angle_deg( &ze, ze_deg );
        set_pulsar( &psr, NULL, NULL, P, r, &al, &ze );

        // choose a random place (within the light cylinder)
        rL = psr.rL;
        X.r = 2.0*rL; // point is, at the moment X.r > rL
        while (X.r >= rL)
        {
            set_point_xyz( &X, RANDS(rL), RANDS(rL), RANDS(rL), POINT_SET_ALL );
        }

        // Get the B and V fields at point X
        calc_fields( &X, &psr, v, &B, &V1, &V2, &A1, &A2, &nsols, NULL );

        if (nsols != 2)
            continue;

        // Test whether V - Vaz is parallel to B
        rho = sqrt(X.rhosq);
        spin_speed = rho*2.0*PI/P;
        set_point_xyz( &Vaz_norm, -X.x[1]/rho, X.x[0]/rho, 0.0, POINT_SET_ALL );
        set_point_xyz( &V_sub_Vaz, v*V1.x[0] - spin_speed*Vaz_norm.x[0],
                                   v*V1.x[1] - spin_speed*Vaz_norm.x[1],
                                   v*V1.x[2] - spin_speed*Vaz_norm.x[2],
                                   POINT_SET_ALL );

        SdB = (V_sub_Vaz.x[0]/V_sub_Vaz.r) * B.x[0] +
              (V_sub_Vaz.x[1]/V_sub_Vaz.r) * B.x[1] +
              (V_sub_Vaz.x[2]/V_sub_Vaz.r) * B.x[2];

        // Test whether V (which has constant magnitude) is perpendicular to A
        V_dot_A = V1.x[0] * A1.x[0] +
                  V1.x[1] * A1.x[1] +
                  V1.x[2] * A1.x[2];

        if ((fabs(SdB - 1.0) < 1e-10) && (fabs(V_dot_A) < 1e-10))
            npassed++;
        else
            printf( "#fail: (cV - Ωφ) . B = %.15e,  V.A = %.15e\n", SdB, V_dot_A );

    }

    printf( "calc_fields() tests:\n" );
    printf( "  total:      %d\n", NTESTS );
    printf( "  successful: %d\n", npassed );
}
