#include <stdlib.h>
#include <stdio.h>
#include "psrgeom.h"
#include <time.h>

// Random number between -x and x
#define  RANDU(x)  ((2.0*(double)rand()/(double)RAND_MAX-1.0)*(x))

int main()
{
    // Set up pulsar
    pulsar psr;

    psr_angle *ra  = NULL;
    psr_angle *dec = NULL;

    psr_angle *al = create_psr_angle_deg( 10.0 );
    psr_angle *ze = create_psr_angle_deg( 15.0 );

    double P = 0.5;
    double r = 1e4;

    set_pulsar( &psr, ra, dec, P, r, al, ze );
    fprintf( stderr, "rL = %.15e\n", psr.rL );

    // Set up point in magnetosphere
    srand( time( NULL ) );
    point X;

    // Set up particle speed
    double v = SPEED_OF_LIGHT;

    // Get the number of velocity solutions
    int nsols;

    // Test both Bfield() and calc_fields() to make sure they return the same
    // thing
    point B;
    point V1, V2;

    /*
    do
    {
        set_point_xyz( &X,
                       RANDU( psr.rL ),
                       RANDU( psr.rL ),
                       RANDU( psr.rL ),
                       POINT_SET_ALL );
    } while (X.rhosq > psr.rL2);
    */

    printf( "# Made with PSRGEOM %s\n", PSRGEOM_VERSION );
    printf( "#x y z Bx By Bz V1x V1y V1z V2x V2y V2z\n" );

    double i, j, k;
    for (i = -1; i <= 1; i += 0.2)
    for (j = -1; j <= 1; j += 0.2)
    for (k = -1; k <= 1; k += 0.2)
    {
        set_point_xyz( &X, (double)i*psr.rL,
                           (double)j*psr.rL,
                           (double)k*psr.rL,
                           POINT_SET_ALL );

        if (X.rhosq > psr.rL2) continue;

        //Bfield( &X, &psr, &B_Bfield );
        calc_fields( &X, &psr, v, &B, &V1, &V2, NULL, NULL, &nsols );

        if (nsols > 0)
        {
            printf( "%lf %lf %lf  ", i, j, k );
            printf( "%.15e %.15e %.15e  ",  B.x[0],  B.x[1],  B.x[2] );
            printf( "%.15e %.15e %.15e  ", V1.x[0], V1.x[1], V1.x[2] );
            printf( "%.15e %.15e %.15e\n", V2.x[0], V2.x[1], V2.x[2] );
        }
    }

    // Clean up
    destroy_psr_angle( ra  );
    destroy_psr_angle( dec );
    destroy_psr_angle( al  );
    destroy_psr_angle( ze  );

    return 0;
}
