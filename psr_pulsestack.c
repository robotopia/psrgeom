/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program simulates pulsestacks from 1D carousels. The user supplies
 * a basic carousel model including a carousel radius, s; a carousel rotation
 * rate, P₄; and a spark angular size, σ. The user also supplies the basic
 * pulsar model (P₁, rₚ, α, ζ/β).
 *
 * The footpoints are then sampled with some preset granularity. For each
 * footpoint, the magnetic field line is climbed, and the particle's
 * trajectory distance to the observed emission point is calculated. The
 * geometry of the emission point is solved to get the observed rotation phase
 * associated with the emission point. The trajectory distance is then used to
 * calculate the carousel rotation phase (modulo the pulsar rotation), and
 * thence the intensity. Finally, retardation effects are included.
 *
 * The final set of sampled intensities will not be distributed evenly in
 * rotation phase, so the final step is to interpolate the values at the
 * desired resolution.
 *
 * The default output of this program is a text file designed to imitate the
 * format of * the output of PSRCHIVE's "pdv" program. The columns are
 *   1) pulse number
 *   2) frequency bin (always 0)
 *   3) phase bin
 *   4) stokes I
 *   5) stokes Q
 *   6) stokes U
 *   7) stokes V (always 0)
 *   8) polarisation angle
 *   9) polarisation angle error (always 0)
 *
 * An alternative output is similar, but without doing the final
 * interpolation. In this case, the format of the output is the same, but the
 * third column becomes a real phase (0 < φ < 1) instead of an integer bin
 * number.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "psrgeom.h"

struct opts
{
    double  al_deg;      // alpha angle in deg
    double  ze_deg;      // zeta angle in deg
    double  be_deg;      // beta angle in deg
    double  P_sec;       // period, in sec
    double  rp;          // the pulsar radius
    char   *outfile;     // name of output file (NULL means stdout)
    double  s;           // normalised carousel radius
    int     npoints;     // number of sampled footpoints
    int     nphases;     // number of sampled phases
    int     npulses;     // number of pulses to output
    int     dipole;      // use a dipole field instead of Deutsch
    int     firstonly;   // for each field line, stop at the first solution
    int     no_interp;   // do not perform the final interpolation
    double  P4;          // the carousel rotation period
    int     nsparks;     // the number of sparks
    double  sigma;       // the spark angular size (measured from mag. pole)
    double  tstep;       // the increment travel distance along field lines
    char   *profile;     // file containing profile information
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

void set_default_options( struct opts *o );
FILE *open_output_file( struct opts *o );
void setup_pulsar( struct opts *o, pulsar *psr );

double *read_profile( char *filename, int *n );
void interp( double *x, double *y, int n,
        double *newx, double *newy, int newn );

void resample_array( double *in, int nin, double *out, int nout );

int main( int argc, char *argv[] )
{
    // Set up struct for command line options and set default values
    struct opts o;
    set_default_options( &o );

    // Parse the command line options
    parse_cmd_line( argc, argv, &o );

    // Set up output file
    FILE *f = open_output_file( &o );

    // Set up pulsar
    pulsar psr;
    setup_pulsar( &o, &psr );

    // If a profile is supplied, read it in
    double *profile_orig = NULL; // The "raw" profile as read in
    double profile[o.nphases];   // The resampled profile
    if (o.profile)
    {
        // Read in profile from file
        int profile_size;
        profile_orig = read_profile( o.profile, &profile_size );

        // Resample the profile for the desired time resolution
        resample_array( profile_orig, profile_size, profile, o.nphases );
    }

    // Write the file and column headers
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // Loop over the polar cap
    int linetype;   // either CLOSED_LINE or OPEN_LINE
    point foot_pt, foot_pt_mag;
    point init_pt;
    int find_emitpt_result;

    // Keep track of dist travelled along field lines
    double dist[o.npoints], dist_tmp;

    int p_idx;
    double p_deg;
    psr_angle p; // The azimuth angle around the polar cap

    // Loop around the polar cap in azimuth and find the emission points
    point emit_pts[o.npoints];

    for (p_idx = 0; p_idx < o.npoints; p_idx++)
    {
        // Convert p_idx to an angle
        p_deg = p_idx * 360.0 / o.npoints - 180.0; // Go from -180° to 180°
        set_psr_angle_deg( &p, p_deg );

        // Convert (s,p) into a point in the magnetic frame
        set_point_sph( &foot_pt_mag, psr.r, &psr.csl.S, &p, POINT_SET_ALL );

        // Convert the foot_pt into observer coordinates
        mag_to_obs_frame( &foot_pt_mag, &psr, NULL, &foot_pt );

        // Now check that we're on an open field line
        linetype = get_fieldline_type( &foot_pt, &psr, 0, NULL,
                NULL, NULL );
        if (linetype == CLOSED_LINE)
        {
            fprintf( stderr, "error: field line at magnetic azimuth p = %.1f "
                             "is closed. Try a smaller value of s.\n", p_deg );
            exit(EXIT_FAILURE);
        }

        // Now find the emission points along this line!
        // Start 1 metre above the surface
        Bstep( &foot_pt, &psr, 1.0, DIR_OUTWARD, &init_pt );
        set_point_xyz( &init_pt, init_pt.x[0],
                                 init_pt.x[1],
                                 init_pt.x[2],
                                 POINT_SET_ALL );
        dist[p_idx] = 0.0; // Start distance tracker

        // Climb up the field line to find the next emit_pt
        find_emitpt_result = find_next_line_emission_point( &psr,
                &init_pt, DIR_OUTWARD, &emit_pts[p_idx], &dist_tmp,
                NULL );
        dist[p_idx] += dist_tmp;

        // If no point was found, put x = y = z = 0
        if (find_emitpt_result != EMIT_PT_FOUND)
        {
            set_point_xyz( &emit_pts[p_idx], 0.0, 0.0, 0.0, POINT_SET_ALL );
        }
    }

    // Calculate the (retarded) phase at which the emission would be seen.
    // First, set the V to the velocity vector. While we're at it, get the
    // acceleration vector.

    point B[o.npoints]; // The magnetic field at the emission point
    point V[o.npoints]; // The velocity field at the emission point
    point A[o.npoints]; // The accelrtn field at the emission point
    point retarded_LoS[o.npoints];

    psr_angle dph[o.npoints]; // The retardation angle
    psr_angle psi[o.npoints]; // The polarisation angle

    psr_angle phase[o.npoints];     // The emission phase
    psr_angle ret_phase[o.npoints]; // The retarded (observed) phase

    psr_angle spark_phase[o.npoints]; // The phase when the particles left
                                      // the surface

    for (p_idx = 0; p_idx < o.npoints; p_idx++)
    {
        calc_fields( &emit_pts[p_idx], &psr, SPEED_OF_LIGHT, &B[p_idx],
                &V[p_idx], NULL, &A[p_idx], NULL, NULL );

        calc_retardation( &emit_pts[p_idx], &psr, &V[p_idx], &dph[p_idx],
                &retarded_LoS[p_idx] );

        // Now, the observed phase is the negative of the azimuthal
        // angle of the retarded line of sight
        set_point_xyz( &V[p_idx], V[p_idx].x[0], V[p_idx].x[1], V[p_idx].x[2],
                POINT_SET_PH );
        if (psr.spin == SPIN_POS)
            set_psr_angle_deg( &phase[p_idx], -V[p_idx].ph.deg );
        else
            copy_psr_angle( &(V[p_idx].ph), &phase[p_idx] );
        set_psr_angle_deg( &ret_phase[p_idx], -(retarded_LoS[p_idx].ph.deg) );

        // Calculate the observed polarisation angle at emit_pt
        accel_to_pol_angle( &psr, &A[p_idx], &phase[p_idx], &psi[p_idx] );

        // For each emission point, calculate the rotation phase at which the
        // particle must have left the surface in order to arrive at the emission
        // point at the correct phase for observing.
        set_psr_angle_rad( &spark_phase[p_idx],
                phase[p_idx].rad + dist[p_idx]/psr.rL );
    }

    // Calculate Stokes I from the 1D carousel model
    int pulse, spark;
    double t, x;

    double **In = (double **)malloc( o.npulses * sizeof(double *) );
    for (pulse = 0; pulse < o.npulses; pulse++)
        In[pulse] = (double *)malloc( o.npoints * sizeof(double) );

    for (pulse = 0; pulse < o.npulses; pulse++)
    {
        for (p_idx = 0; p_idx < o.npoints; p_idx++)
        {
            t = pulse*psr.P + spark_phase[p_idx].deg/360.0;

            p_deg = p_idx * 360.0 / (o.npoints-1) - 180.0; // From -180 to 180
            set_psr_angle_deg( &p, p_deg );

            In[pulse][p_idx] = 0.0;
            for (spark = 0; spark < o.nsparks; spark++)
            {
                x = p.rad -
                    2*PI*((double)spark/(double)o.nsparks + t/psr.csl.P4);
                while (x < -PI) x += 2.0*PI;
                while (x >= PI) x -= 2.0*PI;
                In[pulse][p_idx] += exp(-0.5*x*x/
                                        (psr.csl.s.rad*psr.csl.s.rad));
            }

        }
    }

    // Interpolate the phases
    double **stokesI = (double **)malloc( o.npulses * sizeof(double *) );
    for (pulse = 0; pulse < o.npulses; pulse++)
        stokesI[pulse] = (double *)malloc( o.nphases * sizeof(double) );

    double ph[o.npoints];
    for (p_idx = 0; p_idx < o.npoints; p_idx++)
    {
        ph[p_idx] = ret_phase[p_idx].deg;
    }

    double newph[o.nphases];
    double dp = 360.0 / (o.nphases-1);
    for (p_idx = 0; p_idx < o.nphases; p_idx++)
    {
        newph[p_idx] = p_idx*dp - 180.0;
    }

    for (pulse = 0; pulse < o.npulses; pulse++)
    {
        // Interpolate!
        interp( ph, In[pulse], o.npoints,
                newph, stokesI[pulse], o.nphases );

        // Print out the result
        for (p_idx = 0; p_idx < o.nphases; p_idx++)
        {
            // If requested, modulate with the supplied profile
            if (o.profile)
            {
                stokesI[pulse][p_idx] *= profile[p_idx];
            }

            // Output results
            fprintf( f, "%d %d %d %e\n",
                    pulse, 0, p_idx, stokesI[pulse][p_idx] );
        }
        fprintf( f, "\n" );
    }

    // Clean up

    free( o.outfile );
    free( o.profile );
    free( profile_orig );
    for (pulse = 0; pulse < o.npulses; pulse++)
    {
        free( In[pulse] );
        free( stokesI[pulse] );
    }
    free( In );
    free( stokesI );

    if (o.outfile != NULL)
        fclose( f );

    return 0;
}

void set_default_options( struct opts *o )
{
    // Set up struct for command line options and set default values
    o->al_deg    = NAN;
    o->P_sec     = NAN;
    o->ze_deg    = NAN;
    o->be_deg    = NAN;
    o->rp        = 1e4;
    o->outfile   = NULL;
    o->s         = NAN;
    o->npoints   = 361;
    o->nphases   = 1024;
    o->npulses   = 100;
    o->dipole    = 0;
    o->firstonly = 1;
    o->no_interp = 0;
    o->P4        = NAN;
    o->nsparks   = 0;
    o->sigma     = NAN;
    o->tstep     = 0.01;
    o->profile   = NULL;
}


FILE *open_output_file( struct opts *o )
{
    FILE *f;
    if (o->outfile == NULL)
        f = stdout;
    else
    {
        f = fopen( o->outfile, "w" );
        if (f == NULL)
        {
            fprintf( stderr, "error: could not open file %s\n", o->outfile );
            exit(EXIT_FAILURE);
        }
    }

    return f;
}


void setup_pulsar( struct opts *o, pulsar *psr )
{
    // Set basic pulsar properties
    psr_angle al, ze;
    set_psr_angle_deg( &al, o->al_deg );
    set_psr_angle_deg( &ze, o->ze_deg );

    //               RA    Dec
    set_pulsar( psr, NULL, NULL, o->P_sec, o->rp, &al, &ze );

    if (o->dipole)
        psr->field_type = DIPOLE;

    // Set up carousel
    psr_angle spark_size;
    psr_angle csl_radius;

    set_psr_angle_deg( &spark_size, o->sigma );
    s_to_deg( psr, o->s, &csl_radius );

    set_pulsar_carousel( psr, o->nsparks, &spark_size, &csl_radius, GAUSSIAN,
            o->P4 );
}


double *read_profile( char *filename, int *n )
/* Reads the (ascii) numbers found in the file named FILENAME, and stores
 * them in a newly allocated array of doubles, which must be freed by the
 * caller.
 */
{
    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error: read_profile: could not open file %s\n",
                filename );
        exit(EXIT_FAILURE);
    }

    // Find out how many numbers there are in the file
    int i = 0;
    double d;
    while ((fscanf( f, "%lf", &d) != EOF))  i++;
    *n = i;

    rewind(f);

    // Allocate memory for reading the whole file contents
    double *profile = (double *)malloc( i * sizeof(double) );

    // Go the the file again, and read everything in
    i = 0;
    while ((fscanf( f, "%lf", &d) != EOF))  profile[i++] = d;

    return profile;
}


void interp( double *x, double *y, int n,
        double *newx, double *newy, int newn )
/* Linear interpolation of the function defined by the arrays x and y,
 * evaluated at newx. The result is saved out to newy.
 *
 * It is assumed that x and y have at least size n, and that newx and newy
 * have at least size newn.
 *
 * The array x does not have to be ordered.
 */
{
    int i, newi;
    int i0, in; // indexes for lowest x, highest x
    int il, ir; // indexes for x's that straddle newx value

    // Get lowest and highest x values
    i0 = in = 0;
    for (i = 1; i < n; i++)
    {
        if (x[i] < x[i0])  i0 = i;
        if (x[i] > x[in])  in = i;
    }

    // Go through newx's and find their straddling points in x,
    // and then linearly interpolate
    for (newi = 0; newi < newn; newi++)
    {
        // Make sure we're within the boundaries
        if (newx[newi] < x[i0] || newx[newi] > x[in])
        {
            newy[newi] = NAN;
            continue;
        }

        // Now go through the x's to find the straddling points
        il = i0;
        ir = in;
        for (i = 0; i < n; i++)
        {
            if (x[i] <= newx[newi] && x[il] < x[i])  il = i;
            if (x[i] >= newx[newi] && x[ir] > x[i])  ir = i;
        }

        // If we find an identical x point, no interpolation needed!
        if (x[il] == newx[newi])
        {
            newy[newi] = y[il];
        }
        else
        {
            newy[newi] = (y[ir] - y[il]) * (newx[newi] - x[il]) /
                         (x[ir] - x[il]) + y[il];
        }
    }
}

void resample_array( double *in, int nin, double *out, int nout )
/* Resample the input array IN of size NIN to the output array OUT of size
 * NOUT. This function assumes that sufficient memory for the arrays has
 * already been allocated.
 * 
 * The resampling is done in Fourier space, and this algorithm is a straight-
 * forward one that uses the FFTW3 library.
 *
 * It is assumed that IN and OUT point to different arrays.
 */
{
    // Make sure supplied array sizes are positive definite
    if (nin <= 0 || nout <= 0)
    {
        fprintf( stderr, "error: resample_array: array sizes must be >= 1\n" );
        exit(EXIT_FAILURE);
    }

    int i; // Generic counter

    // If the sizes are the same, just copy the data across
    if (nin == nout)
    {
        for (i = 0; i < nin; i++)
        {
            out[i] = in[i];
        }

        return;
    }

    // In the general case, FFTs will be needed
    fftw_complex *F; // The Fourier representation of the input array

    fftw_plan pf, pb; // (f)orward and (b)ackward plans

    int max_size = (nin > nout ? nin : nout);

    F = (fftw_complex*)fftw_malloc( max_size * sizeof(fftw_complex) );

    pf = fftw_plan_dft_r2c_1d( nin,  in,  F, FFTW_ESTIMATE );
    pb = fftw_plan_dft_c2r_1d( nout, F, out, FFTW_ESTIMATE );

    // Switch to frequency domain
    fftw_execute( pf );

    if (nout > (nin/2+1)) // i.e. need to upsample --> zero-pad
    {
        for (i = nin/2+1; i < nout; i++)
            F[i] = 0.0 + 0.0*I;
    }

    // Now change back to time domain
    fftw_execute( pb );

    // Clean up
    fftw_destroy_plan( pf );
    fftw_destroy_plan( pb );
    fftw_free( F );
}


void usage()
{
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -4  period   The rotation period of the carousel, in "
                           "seconds\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees\n" );
    printf( "  -b  beta     The \"impact\" angle between the magnetic axis "
                           "and the line of sight in degrees (either -z or "
                           "-b is required)\n" );
    printf( "  -N  sparks   The number of sparks in the carousel. Must be "
                           "> 0\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds\n" );
    printf( "  -s  radius   The carousel radius, normalised to the (aligned) "
                           "polar cap radius\n" );
    printf( "  -S  sigma    The angular size of the spark gaussian "
                           "profiles, in degrees\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (either -z or -b is "
                           "required, but -z trumps -b)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -1           Only report one (i.e. the first) solution for "
                           "each field line (default: off)\n" );
    printf( "  -d           Use a dipole field instead of the default "
                           "Deutsch field\n" );
    printf( "  -e  file     File containing profile which will be used to "
                           "modulate the pulsestack.\n"
            "               The file must contain a single list of whitespace-"
                           "separated numbers indicating total intensity.\n"
            "               It is assumed that the numbers span all 360 deg "
                           "of rotation phase (from -180 to 180).\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -i           Don't interpolate in pulse phase\n" );
    printf( "  -l  pulses   Number of pulses to output (default: 100)\n" );
    printf( "  -n  points   The number of footpoints to sample "
                           "(default: 361)\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
    printf( "  -p  phases   The number of (interpolated) phase points "
                           "(default: 1024). Has no effect if -i is "
                           "supplied\n" );
    printf( "  -r  radius   The pulsar radius, in km (default: 10)\n" );
    printf( "  -t  step     Step size for moving along magnetic field lines, "
                           "as a fraction of the light cylinder radius "
                           "(default: 0.01)\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "14:a:b:de:hil:n:N:o:p:P:r:s:S:t:z:")) != -1)
    {
        switch (c)
        {
            case '1':
                o->firstonly = 1;
                break;
            case '4':
                o->P4 = atof(optarg);
                break;
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'b':
                o->be_deg = atof(optarg);
                break;
            case 'd':
                o->dipole = 1;
                break;
            case 'e':
                o->profile = strdup(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'i':
                o->no_interp = 1;
                break;
            case 'l':
                o->npulses = atoi(optarg);
                break;
            case 'n':
                o->npoints = atoi(optarg);
                break;
            case 'N':
                o->nsparks = atoi(optarg);
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'p':
                o->nphases = atoi(optarg);
                break;
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 'r':
                o->rp = atof(optarg)*1e3; // km -> m
                break;
            case 's':
                o->s = atof(optarg);
                break;
            case 'S':
                o->sigma = atof(optarg);
                break;
            case 't':
                o->tstep = atof(optarg);
                break;
            case 'z':
                o->ze_deg = atof(optarg);
                break;
            case '?':
                fprintf( stderr, "error: unknown option character '-%c'\n",
                         optopt );
                exit(EXIT_FAILURE);
                break;
            default:
                fprintf( stderr, "error: couldn't parse command line\n" );
                exit(EXIT_FAILURE);
        }
    }

    // Check that all the arguments are valid
    if ( isnan(o->al_deg) || isnan(o->P_sec) ||
        (isnan(o->ze_deg) && isnan(o->be_deg)))
    {
        fprintf( stderr, "error: -a, -P, and -z/-b options required"
                         "\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (isnan(o->s) || isnan(o->P4) || isnan(o->sigma))
    {
        fprintf( stderr, "error: -s, -S, and -4 options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (o->nsparks <= 0)
    {
        fprintf( stderr, "error: number of sparks (-N) must be > 0\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (o->npulses <= 0)
    {
        fprintf( stderr, "error: number of pulses (-l) must be > 0\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (o->nphases <= 0 && !o->no_interp)
    {
        fprintf( stderr, "error: number of phase bins (-p) must be > 0\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (isnan(o->ze_deg))
        o->ze_deg = o->al_deg + o->be_deg;
}


void print_col_headers( FILE *f )
/* The final output includes:
 *   1) pulse number
 *   2) frequency bin (always 0)
 *   3) phase bin / phase
 *   4) stokes I
 *   5) stokes Q
 *   6) stokes U
 *   7) stokes V (always 0)
 *   8) polarisation angle
 *   9) polarisation angle error (always 0)
 */
{
    // Print out a line to file handle f
    //fprintf( f, "#  pulse_no  freq_bin  phase  I  Q  U  V  PA  PA_err\n" );
    fprintf( f, "#  pulse_no  freq_bin  phase  I\n" );
}


