#ifndef PSRGEOM_H
#define PSRGEOM_H

#include <math.h>

#define _USE_MATH_DEFINES

#define  PI              M_PI
#define  SPEED_OF_LIGHT  299792458.0

#define  DEG2RAD         PI/180.0
#define  RAD2DEG         180.0/PI

#define  NDEP            3     // Number of dimensions

#define POINT_SET_NONE   0x0   /* 00000000 */
#define POINT_SET_X      0x1   /* 00000001 */
#define POINT_SET_Y      0x2   /* 00000010 */
#define POINT_SET_Z      0x4   /* 00000100 */
#define POINT_SET_R      0x8   /* 00001000 */
#define POINT_SET_TH     0x10  /* 00010000 */
#define POINT_SET_PH     0x20  /* 00100000 */
#define POINT_SET_RHOSQ  0x40  /* 01000000 */
#define POINT_SET_XYZ    0x7   /* 00000111 */
#define POINT_SET_SPH    0x38  /* 00111000 */
#define POINT_SET_ALL    0x7F  /* 01111111 */



/**** Structures ****/

typedef struct angle_t
{
    double deg; // The angle in degrees
    double rad; // The angle in radians
    double cos; // The cosine of the angle
    double sin; // The sine of the angle
} angle;

typedef struct point_t {
    double  x[3];   // x,y,z coordinates
    double  r;      // distance from origin = sqrt(x^2 + y^2 + z^2)
    angle   th;     // The colatitude
    angle   ph;     // The longitude
    double  rhosq;  // distance from z axis squared = x^2 + y^2
} point;

typedef struct pulsar_t {
    angle    ra;    // Right Ascension
    angle    dec;   // Declination
    double   P;     // Rotation period
    //double   Pdot;  // First time derivative of rotation period
    double   r;     // Stellar radius
    double   rL;    // Light cylinder radius
    angle    al;    // Angle between the rotation and magnetic axes
    angle    ze;    // Angle between the rotation axis and the line of sight
} pulsar;


/**** Angle functions ****/

angle *create_angle();
angle *create_angle_rad( double rad );
angle *create_angle_deg( double deg );

void copy_angle( angle *src, angle *dest );
void destroy_angle( angle *ang );

void set_angle_rad( angle *ang, double rad );
void set_angle_deg( angle *ang, double deg );
void set_angle_sin( angle *ang, double Sin );
void set_angle_cos( angle *ang, double Cos );

void rotatex( angle *th, angle *ph, angle *th_out, angle *ph_out, angle *rot );



/**** Point functions ****/

void set_point_xyz( point *p, double x, double y, double z, int flags );
void set_point_sph( point *p, double r, angle *th, angle *ph, int flags );
void copy_point( point *src, point *dest );



/**** Pulsar functions ****/
void set_pulsar( pulsar *psr, angle *ra, angle *dec, double P, double r,
        angle *al, angle *ze );
pulsar *create_pulsar( angle *ra, angle *dec, double P, double r,
        angle *al, angle *ze );
void destroy_pulsar( pulsar *psr );

void set_pulsar_period( pulsar *psr, double P );

/**** Other functions ****/

double light_cylinder( double P );

#endif
