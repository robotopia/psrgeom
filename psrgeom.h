#ifndef PSRGEOM_H
#define PSRGEOM_H

#define PSRGEOM_VERSION "1.1.7"

#include <stdio.h>
#include <math.h>

#define _USE_MATH_DEFINES

#define  PI              M_PI
#define  SPEED_OF_LIGHT  299792458.0

#define  DEG2RAD         PI/180.0
#define  RAD2DEG         180.0/PI

#define  NDEP            3     // Number of dimensions

// The ways to set the values in a point
#define  POINT_SET_NONE   0x0   /* 00000000 */
#define  POINT_SET_X      0x1   /* 00000001 */
#define  POINT_SET_Y      0x2   /* 00000010 */
#define  POINT_SET_Z      0x4   /* 00000100 */
#define  POINT_SET_R      0x8   /* 00001000 */
#define  POINT_SET_TH     0x10  /* 00010000 */
#define  POINT_SET_PH     0x20  /* 00100000 */
#define  POINT_SET_RHOSQ  0x40  /* 01000000 */
#define  POINT_SET_XYZ    0x7   /* 00000111 */
#define  POINT_SET_SPH    0x38  /* 00111000 */
#define  POINT_SET_ALL    0x7F  /* 01111111 */

// Directions for traversing field lines
#define  DIR_INWARD      -1
#define  DIR_OUTWARD     1
#define  DIR_STOP        0

// The possible return values of footpoint() and farpoint()
#define  STOP_FOUND   0
#define  STOP_EXCEED  1

// Field line types
#define  OPEN_LINE    0
#define  CLOSED_LINE  1

// Generate random numbers
#define RAND(x)   ((x)*(double)rand()/(double)RAND_MAX)   /* 0 < rand < x */
#define RANDU     (RAND(1.0) * 2.0 - 1.0)                 /* -1 < rand < 1 */

/**** Structures ****/

typedef struct psr_angle_t
{
    double deg; // The angle in degrees
    double rad; // The angle in radians
    double cos; // The cosine of the angle
    double sin; // The sine of the angle
} psr_angle;

#define  ANGLE_ZERO  (psr_angle){0.0,0.0,1.0,0.0}
#define  ANGLE_RIGHT (psr_angle){90.0,PI/2.0,0.0,1.0}
#define  ANGLE_HALF  (psr_angle){180.0,PI,-1.0,0.0}
#define  ANGLE_FULL  (psr_angle){360.0,2.0*PI,1.0,0.0}

typedef struct point_t {
    double    x[3];   // x,y,z coordinates
    double    r;      // distance from origin = sqrt(x^2 + y^2 + z^2)
    psr_angle th;     // The colatitude
    psr_angle ph;     // The longitude
    double    rhosq;  // distance from z axis squared = x^2 + y^2
} point;

typedef struct pulsar_t {
    psr_angle    ra;    // Right Ascension
    psr_angle    dec;   // Declination
    double   P;     // Rotation period
    psr_angle    Om;    // Rotation frequency = 2π/P
    //double   Pdot;  // First time derivative of rotation period
    double   r;     // Stellar radius
    double   rL;    // Light cylinder radius
    double   rL2;   // = rL^2 (because it comes up quite often)
    psr_angle    al;    // Angle between the rotation and magnetic axes
    psr_angle    ze;    // Angle between the rotation axis and the line of sight
} pulsar;



/**** Angle functions ****/

psr_angle *create_psr_angle();
psr_angle *create_psr_angle_rad( double rad );
psr_angle *create_psr_angle_deg( double deg );

void copy_psr_angle( psr_angle *src, psr_angle *dest );
void destroy_psr_angle( psr_angle *ang );

void set_psr_angle_rad( psr_angle *ang, double rad );
void set_psr_angle_deg( psr_angle *ang, double deg );
void set_psr_angle_sin( psr_angle *ang, double Sin );
void set_psr_angle_cos( psr_angle *ang, double Cos );

void rotate_about_axis( point *in, point *out, psr_angle *rot, char axis,
                        int flags );

void min_phase_diff( psr_angle *a1, psr_angle *a2, psr_angle *diff );



/**** Point functions ****/

void set_point_xyz( point *p, double x, double y, double z, int flags );
void set_point_sph( point *p, double r, psr_angle *th, psr_angle *ph, int flags );
void set_point_cyl( point *p, double rh, psr_angle *ph, double z, int flags );
void copy_point( point *src, point *dest );

double norm_dot( point *p1, point *p2 );



/**** Pulsar functions ****/

void set_pulsar( pulsar *psr, psr_angle *ra, psr_angle *dec, double P, double r,
        psr_angle *al, psr_angle *ze );
pulsar *create_pulsar( psr_angle *ra, psr_angle *dec, double P, double r,
        psr_angle *al, psr_angle *ze );
void destroy_pulsar( pulsar *psr );

void set_pulsar_period( pulsar *psr, double P );
double light_cylinder( double P );
void line_of_sight( pulsar *psr, psr_angle *phase, point *LoS );



/**** Magnetic field functions ****/

void calc_fields( point *X, pulsar *psr, double v,
                  point *B1,
                  point *V1, point *V2,
                  point *A1, point *A2,
                  int *nsols );

void Bstep( point *x1, pulsar *psr, double tstep, int direction, point *x2 );
double Bdotrxy( point *x, pulsar *psr );
int cmp_extreme( point *x, pulsar *psr, double precision );
int footpoint( point *start_pt, pulsar *psr, double tmult, int direction,
               FILE *write_xyz, int rL_norm, double rL_lim, point *foot_pt );
int farpoint( point *start_pt, pulsar *psr, double tmult,
              FILE *write_xyz, int rL_norm, double rL_lim, point *far_pt );


/**** Dipole field functions ****/

void obs_to_mag_frame( point *xo, pulsar *psr, psr_angle *ph, point *xm );
void mag_to_obs_frame( point *xm, pulsar *psr, psr_angle *ph, point *xo );

double calc_dipole_R( point *xm );
void dipole_footpoint( pulsar *psr, double R, psr_angle *si, point *foot_pt );

void beamangle_to_posangle( psr_angle *ba, psr_angle *pa );
void posangle_to_beamangle( psr_angle *pa, psr_angle *ba );


/**** Other functions ****/

void move_around_cyl( point *start_pt, psr_angle *xi, double zdist,
                      point *end_pt );

/**** PSRGEOM I/O ****/
void print_psrg_header( FILE *f, int argc, char *argv[] );

/**** Finding the emission point ****/

double psr_cost_lofl( point *X, pulsar *psr );
double psr_cost_los( point *X, pulsar *psr, psr_angle *phase, int direction );
int get_fieldline_type( point *X, pulsar *psr, double tmult );

void find_approx_emission_point( pulsar *psr, psr_angle *phase,
                                 point *emit_pt );
void find_emission_point_nmead( pulsar *psr, psr_angle *phase, int direction,
                                point *emit_pt, FILE *f );
void find_emission_point_newuoa( pulsar *psr, psr_angle *phase, int direction,
                                 point *emit_pt, FILE *f );
int find_emission_point_elevator( pulsar *psr, psr_angle *phase,
        int direction, point *emit_pt, FILE *f );

void psr_cost_deriv( point *X, pulsar *psr, psr_angle *phase, int direction,
                     double dx, point *grad );
void find_LoS_at_r( point *init_pt, pulsar *psr, psr_angle *phase,
                    int direction, point *end_pt, FILE *f );

#endif
