#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psrgeom.h"

void usage()
{
    printf( "spark_test [nsparks] [S_deg] [s_deg] [P4_sec] [time]\n" );
    printf( "    nsparks = number of sparks (0 = annulus)\n" );
    printf( "    S_deg   = angle between magnetic axis and carousel ring\n" );
    printf( "    s_deg   = angular size of sparks\n" );
    printf( "    P4_sec  = carousel rotation time in sec\n" );
    printf( "    time    = time in seconds\n" );
}

int main( int argc, char *argv[] )
{
    // Set up a pulsar
    pulsar psr;

    psr_angle *ra  = NULL;
    psr_angle *dec = NULL;

    psr_angle al, ze;
    set_psr_angle_deg( &al, 0.0 );
    set_psr_angle_deg( &ze, 0.0 );

    double P = 1.0;
    double r = 1e4; /* This doesn't actually make a difference to any of the
                       outputs of this program */

    set_pulsar( &psr, ra, dec, P, r, &al, &ze );

    // Set up its carousel
    if (argc < 6)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    int n = atoi(argv[1]);
    psr_angle S, s;
    set_psr_angle_deg( &S, atof(argv[2]) );
    set_psr_angle_deg( &s, atof(argv[3]) );
    double P4 = atof(argv[4]);

    set_pulsar_carousel( &psr, n, &s, &S, GAUSSIAN, P4 );

    double t = atof(argv[5]); // The time!
    double h;                 // The height of the profile

    // Sample the sparks!
    point foot_pt;

    psr_angle th, ph;
    double th_deg_max = S.deg + s.deg*3.0;

    int x_idx, y_idx;
    double x, y;
    int X = 100, Y = 100;

    double x_max = th_deg_max;
    double y_max = th_deg_max;
    double x_min = -x_max;
    double y_min = -y_max;

    double x_size = (x_max - x_min) / (double)X;
    double y_size = (y_max - y_min) / (double)Y;

    for (x_idx = 0; x_idx < X; x_idx++)
    {
        x = x_idx*x_size + x_min;
        for (y_idx = 0; y_idx < Y; y_idx++)
        {
            y = y_idx*y_size + y_min;

            set_psr_angle_deg( &th, sqrt(x*x + y*y) );
            set_psr_angle_rad( &ph, atan2( y, x ) );

            set_point_sph( &foot_pt, psr.r, &th, &ph, POINT_SET_ALL );
            h = spark_profile( &psr, t, &foot_pt );

            printf( "%.15e %.15e %.15e\n", x, y, h );
        }

        printf( "\n" );
    }
}
