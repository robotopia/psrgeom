#ifndef PSRGEOM_H
#define PSRGEOM_H

#define PSRGEOM_VERSION "1.4a.5"

#include <stdio.h>
#include <math.h>

#define _USE_MATH_DEFINES

#define  PI              M_PI
#define  SPEED_OF_LIGHT  299792458.0

#define  DEG2RAD         PI/180.0
#define  RAD2DEG         180.0/PI

#define  NDEP            3     /* Number of dimensions */

#define  MAX_BSTEP       0.005 /* fraction of light cylinder radius */

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

// The possible return values of find_emission_point_elevator()
#define  EMIT_PT_TOO_LOW  -1
#define  EMIT_PT_FOUND     0
#define  EMIT_PT_TOO_HIGH  1

// Field types
#define DIPOLE   0
#define DEUTSCH  1

// Field line types
#define  OPEN_LINE    0
#define  CLOSED_LINE  1

// Different weighting options (used in weight_photon_total())
#define WEIGHT_GAMMA  0x1
#define WEIGHT_POWER  0x2
#define WEIGHT_SPARK  0x4
#define WEIGHT_LDENS  0x8
#define WEIGHT_TOTAL  0xF

// Generate random numbers
#define RAND(x)   ((x)*(double)rand()/(double)RAND_MAX)   /* 0 < rand < x */
#define RANDU     (RAND(1.0) * 2.0 - 1.0)                 /* -1 < rand < 1 */
#define RANDTH    (acos(1.0-2.0*(RAND(1.0))))
#define RANDTHAB(a_rad,b_rad)  (acos(cos(a_rad)+(cos(b_rad)-cos(a_rad))*(RAND(1.0))))

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

#define SPIN_POS   1
#define SPIN_NEG  -1

typedef struct point_t
{
    double    x[3];   // x,y,z coordinates
    double    r;      // distance from origin = sqrt(x^2 + y^2 + z^2)
    psr_angle th;     // The colatitude
    psr_angle ph;     // The longitude
    double    rhosq;  // distance from z axis squared = x^2 + y^2
} point;

// Carousel types
#define  TOPHAT    0   /* Uniform distribution, sharp cutoff */
#define  GAUSSIAN  1   /* Sparks have 2D Gaussian profiles */

typedef struct carousel_t
{
    int        n;    // Number of sparks in the carousel (0 = annulus)
    psr_angle  s;    // The angular radius of a spark
    psr_angle  S;    // The angular radius of the whole carousel
    int        type; // {TOPHAT, GAUSSIAN}
    double     P4;   // The rotation rate of the carousel (in sec)
} carousel;


typedef struct pulsar_t
{
    psr_angle    ra;    // Right Ascension
    psr_angle    dec;   // Declination
    double       P;     // Rotation period
    psr_angle    Om;    // Rotation frequency = 2π/P
    //double       Pdot;  // First time derivative of rotation period
    double       r;     // Stellar radius
    double       rL;    // Light cylinder radius
    double       rL2;   // = rL^2 (because it comes up quite often)
    psr_angle    pcr;   // the polar cap radius (as an angle)
    psr_angle    al;    // Angle between the rotation and magnetic axes
    psr_angle    ze;    // Angle between the rotation axis and the LoS
    int          spin;  // Spin direction (SPIN_POS or SPIN_NEG)
    int          field_type; // Can be DIPOLE or DEUTSCH
    carousel     csl;   // The carousel parameters
} pulsar;

typedef struct photon_t
{
    double    freq;         // The photon's frequency, in Hz
    point     source;       // The location of the particle when it emitted
    point     B, V, A;      // The B-, V-, and A-fields at the source location
    point     retarded_LoS; // Same as V, but adjusted for retardation
    double    curvature;    // The curvature of the particle's trajectory
    double    gamma;        // The Lorentz factor of the particle
    psr_angle phase;        // The phase at which the particle is observed
    psr_angle psi;          // The photon's polarisation angle (kind of...)
    double    power;        // The photon's power (Jackson Eq (14.46))
} photon;


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

void reverse_psr_angle( psr_angle *in, psr_angle *out );
void rotate_about_axis( point *in, point *out, psr_angle *rot, char axis,
                        int flags );

void min_phase_diff( psr_angle *a1, psr_angle *a2, psr_angle *diff );



/**** Point functions ****/

void set_point_xyz( point *p, double x, double y, double z, int flags );
void set_point_sph( point *p, double r, psr_angle *th, psr_angle *ph, int flags );
void set_point_cyl( point *p, double rh, psr_angle *ph, double z, int flags );
void copy_point( point *src, point *dest );

double norm_dot( point *p1, point *p2 );

void spherical_midpoint( point *p1, point *p2, point *mid_pt, int flags );

void scale_point( point *in, double scale, point *out );

void transform_new_xz( point *v, point *new_z, point *new_x, point *new_v );

double randn( double mu, double sigma );

/**** Pulsar functions ****/

void set_pulsar( pulsar *psr, psr_angle *ra, psr_angle *dec, double P, double r,
        psr_angle *al, psr_angle *ze );
pulsar *create_pulsar( psr_angle *ra, psr_angle *dec, double P, double r,
        psr_angle *al, psr_angle *ze );
void destroy_pulsar( pulsar *psr );

void set_pulsar_period( pulsar *psr, double P );
double light_cylinder( double P );
void line_of_sight( pulsar *psr, psr_angle *phase, point *LoS );
void pol_zero( pulsar *psr, psr_angle *phase, point *pz );

void calc_retardation( point *X, pulsar *psr, point *LoS,
        psr_angle *dph, point *retarded_LoS );
void set_polar_cap_radius( pulsar *psr );
void s_to_deg( pulsar *psr, double s, psr_angle *s_angle );

void set_pulsar_carousel( pulsar *psr, int n, psr_angle *s, psr_angle *S,
        int type, double P4 );

void random_direction( point *rand_pt );
void random_direction_bounded( point *rand_pt, double lo_th_rad,
        double hi_th_rad, double lo_ph_rad, double hi_ph_rad );
void random_spark_footpoint( point *foot_pt_obs, point *foot_pt_mag,
                pulsar *psr, double t );
void random_point_in_cyl( point *rand_pt, double max_rho, double max_z );
void random_point_in_lightcyl( point *rand_pt, pulsar *psr, double frac,
        double z_max_frac );

double neg_power_law_distr( double lo, double hi, double index );
double spark_profile( pulsar *psr, double t, point *foot_pt );

/**** Magnetic field functions ****/

void calc_fields( point *X, pulsar *psr, double v, point *B1, point *V1,
        point *V2, point *A1, point *A2, int *nsols, double *dfact );
void calc_dipole_fields( point *X, pulsar *psr, double v, point *B1,
        point *V1, point *V2, point *A1, point *A2, int *nsols,
        double *dfact );
void calc_deutsch_fields( point *X, pulsar *psr, double v, point *B1,
        point *V1, point *V2, point *A1, point *A2, int *nsols,
        double *dfact );

void Bstep( point *x1, pulsar *psr, double tstep, int direction, point *x2,
        double *ifdist );
void B_large_step( point *x1, pulsar *psr, double step, int direction,
        point *x2 );
void traj_step( point *x1, double t, pulsar *psr, double tstep, int direction,
                point *x2, int rL_norm, FILE *f );
double Bdotrxy( point *x, pulsar *psr );
int cmp_extreme( point *x, pulsar *psr, double precision );
int footpoint( point *start_pt, pulsar *psr, int direction, FILE *write_xyz,
        int rL_norm, double rL_lim, point *foot_pt );
int farpoint( point *start_pt, pulsar *psr, FILE *write_xyz, int rL_norm,
        double rL_lim, double *dist, point *far_pt );
int calc_pol_angle( pulsar *psr, psr_angle *phase, int direction,
                    point *init_pt, point *emit_pt, psr_angle *psi );
void accel_to_pol_angle( pulsar *psr, point *A, psr_angle *phase,
        psr_angle *psi );
double calc_curvature( point *V, point *A );


/**** Jackson's Classical Electrodynamics ****/

double calc_crit_freq( double gamma, double curvature );
double calc_crit_gamma( double crit_freq, double curvature );
void particle_beam_intensity( double freq, double gamma, psr_angle *theta,
        double rho, double *Ipos, double *Ineg );
double single_particle_power( point *n, point *beta, point *beta_dot );
double single_particle_power_perp( double gamma, psr_angle *th, psr_angle *ph,
        double vdot );
double single_particle_power_total( double gamma, double vdot );


/**** Dipole field functions ****/

void obs_to_mag_frame( point *xo, pulsar *psr, psr_angle *ph, point *xm );
void mag_to_obs_frame( point *xm, pulsar *psr, psr_angle *ph, point *xo );

double calc_dipole_R( point *xm );
double calc_dipole_curvature( point *xm );
void dipole_footpoint( pulsar *psr, double R, psr_angle *si, point *foot_pt );

void beamangle_to_posangle( psr_angle *ba, psr_angle *pa );
void posangle_to_beamangle( psr_angle *pa, psr_angle *ba );


/**** Other functions ****/

void move_around_cyl( point *start_pt, psr_angle *xi, double zdist,
                      point *end_pt );

/**** PSRGEOM I/O ****/
void print_psrg_header( FILE *f, int argc, char *argv[] );
void parse_range( char *str, double *start, double *stop, int *nsteps );
void parse_direction( char *str, point *direction );

/**** Finding the emission point ****/

double psr_cost_lofl( point *X, pulsar *psr );
double psr_cost_los( point *X, pulsar *psr, psr_angle *phase, int direction,
                     int retardation );
int get_fieldline_type( point *X, pulsar *psr, int rL_norm, FILE *f,
        double *dist, point *far_pt );

void find_approx_emission_point( pulsar *psr, psr_angle *phase, int direction,
                                 point *emit_pt );
void find_emission_point_nmead( pulsar *psr, psr_angle *phase, int direction,
                                point *emit_pt, FILE *f );
void find_emission_point_newuoa( pulsar *psr, psr_angle *phase, int direction,
                                 point *emit_pt, FILE *f );
int find_emission_point_elevator( pulsar *psr, psr_angle *phase,
        int direction, point *init_pt, point *emit_pt, FILE *f );

int find_next_line_emission_point( pulsar *psr, point *init_pt, int direction,
        point *emit_pt, double *dist, FILE *f );

void psr_cost_deriv( point *X, pulsar *psr, psr_angle *phase, int direction,
                     double dx, point *grad );
void find_LoS_at_r( point *init_pt, pulsar *psr, psr_angle *phase,
                    int direction, point *end_pt, FILE *f );

void climb_and_emit( pulsar *psr, point *init_pt, double gamma,
        double freq_lo, double freq_hi, FILE *f );

void fieldline_to_profile( pulsar *psr, point *init_pt, double freq_lo,
        double freq_hi, int nbins, int centre_bin, double *profile,
        int *bin_count );

/**** Photon functions ****/

void emit_pulsar_photon( pulsar *psr, point *pt, double freq, photon *pn );
double weight_photon_by_particle_density( photon *pn );
double weight_photon_by_line_density( point *init_pt, pulsar *psr );
double weight_photon_by_power( photon *pn );
double weight_photon_by_gamma_distr( photon *pn, double index );
double weight_photon_by_spark( point *foot_pt, photon *pn, pulsar *psr,
        double height, int pulse_number );
double weight_photon_total( photon *pn, pulsar *psr, point *foot_pt,
        double height, int pulse_number, double index, int weight_flags );

void climb_and_emit_photon( pulsar *psr, point *foot_pt, double height,
        double freq, photon *pn );

void bin_photon_freq( photon *pn, double weight, double fmin,
        double fbinwidth, int nfbins, double *farray );
void bin_photon_beam( photon *pn, pulsar *psr, double weight, double xmin,
        double ymin, double xbinwidth, double ybinwidth, int nxbins,
        int nybins, double **xyarray );
void bin_photon_pulsestack( photon *pn, double weight, int pulse_number,
        int pulsemin, int npulses, double phasemin, double phasebinwidth,
        int nphasebins, double **pulsestack );
void bin_photon_profile( photon *pn, double weight, double phasemin,
        double phasebinwidth, int nphasebins, double *profile );

/**** Finding the pulse width ****/

double fitwidth( pulsar *psr, int direction, double width_rad,
              psr_angle *ph1, psr_angle *ph2,
              point *ph1_pt, point *ph2_pt, FILE *f );

/**** Functions from Numerical Recipes for C ****/

void nrerror(char *error_text);
double chebev(double a, double b, double c[], int m, double x);
void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);

#endif
