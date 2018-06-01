#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "psrgeom.h"

void usage()
{
    printf( "spark_random_test [nsparks] [S_deg] [s_deg] [P4_sec] [time] [type] [npoints]\n" );
    printf( "    nsparks = number of sparks (0 = annulus)\n" );
    printf( "    S_deg   = angle between magnetic axis and carousel ring\n" );
    printf( "    s_deg   = angular size of sparks\n" );
    printf( "    P4_sec  = carousel rotation time in sec\n" );
    printf( "    time    = time in seconds\n" );
    printf( "    type    = profile type, either TOPHAT or GAUSSIAN\n" );
    printf( "    npoints = number of points to randomly sample from spark pattern\n" );
}

int main( int argc, char *argv[] )
{
    // Seed random number generator
    srand( time( NULL ) );

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
    if (argc < 8)
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
    if (strcmp( argv[6], "TOPHAT" ) == 0)
        psr.csl.type = TOPHAT;
    else if (strcmp( argv[6], "GAUSSIAN" ) == 0)
        psr.csl.type = GAUSSIAN;
    else
    {
        fprintf( stderr, "error: unknown carousel type %s\n", argv[6] );
        usage();
        exit(EXIT_FAILURE);
    }

    double t = atof(argv[5]); // The time!
    int N = atoi(argv[7]);

    // Sample the sparks!
    point foot_pt;
    double x, y;

    int i;
    for (i = 0; i < N; i++)
    {
        random_spark_footpoint( &foot_pt, &psr, t );
        x = foot_pt.th.deg * foot_pt.ph.cos;
        y = foot_pt.th.deg * foot_pt.ph.sin;
        printf( "%.15e %.15e\n", x, y );
    }

    return 0;
}
