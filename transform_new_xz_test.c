#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "psrgeom.h"

int main()
{
    srand( time( NULL ) );

    point new_x, new_y, new_z, v, new_v, tmp;

    random_direction( &new_z );
    random_direction( &tmp );
    random_direction( &v );

    set_point_xyz( &new_x, new_z.x[1] * tmp.x[2] - new_z.x[2] * tmp.x[1],
                           new_z.x[2] * tmp.x[0] - new_z.x[0] * tmp.x[2],
                           new_z.x[0] * tmp.x[1] - new_z.x[1] * tmp.x[0],
                           POINT_SET_XYZ | POINT_SET_R );
    set_point_xyz( &new_x, new_x.x[0]/new_x.r,
                           new_x.x[1]/new_x.r,
                           new_x.x[2]/new_x.r,
                           POINT_SET_ALL );

    set_point_xyz( &new_y, new_z.x[1] * new_x.x[2] - new_z.x[2] * new_x.x[1],
                           new_z.x[2] * new_x.x[0] - new_z.x[0] * new_x.x[2],
                           new_z.x[0] * new_x.x[1] - new_z.x[1] * new_x.x[0],
                           POINT_SET_ALL );

    transform_new_xz( &v, &new_z, &new_x, &new_v );

    // output answer in a octave-checkable way
    printf( "v = [%.15e %.15e %.15e]';\n", v.x[0], v.x[1], v.x[2] );
    printf( "new_v = [%.15e %.15e %.15e]';\n", new_v.x[0], new_v.x[1], new_v.x[2] );
    printf( "new_x = [%.15e %.15e %.15e]';\n", new_x.x[0], new_x.x[1], new_x.x[2] );
    printf( "new_y = [%.15e %.15e %.15e]';\n", new_y.x[0], new_y.x[1], new_y.x[2] );
    printf( "new_z = [%.15e %.15e %.15e]';\n", new_z.x[0], new_z.x[1], new_z.x[2] );
    printf( "T = inv([new_x, new_y, new_z]);\n" );
    printf( "T*v - new_v\n" );
}
