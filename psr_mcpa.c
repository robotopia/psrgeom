/*****************************************************************************
 * PSRGEOM
 * Sam McSweeney, 2018
 *
 * This program attempts to simulate the spread of polarisation angles
 * observed from a pulsar with given angles α and ζ, and period P.
 * Here, the idea is that emission is allowed to originate on any open field
 * line. No attempt is made to sample the active field lines realistically.
 * To keep things simple, we start at the magnetic pole and work our way
 * outwards in concentric circles on the polar cap until we reach a "radius"
 * where no open field lines are found.
 *
 * For each open field line, the line is followed outward in incremental steps
 * until the velocity field is found to make an angle ζ with the rotation
 * axis, which means that there is some rotation phase at which curvature
 * radiation from that point becomes visible. It is possible that there are
 * multiple emission heights that fit the ζ criterion, so even after a point
 * has been successfully found, the algorithm continues to search for other
 * points further along the field line.
 *
 * For each point that satisfies the ζ criterion, the acceleration vector is
 * evaluated and the polarisation angle found. Before printing out the
 * polarisation, the rotation phase is adjusted for retardation effects.
 *
 * The final output includes:
 *   1) the polar coordinates of the footpoint of the magnetic field line
 *   2) the Cartesian coordinates of the emission point
 *   3) the retarded (observed) emission phase
 *   4) the observed polarisation angle
 *   5) the curvature of the particle's trajectory at the emission point
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "psrgeom.h"

struct opts
{
    double  al_deg;      // alpha angle in deg
    double  ze_deg;      // zeta angle in deg
    double  P_sec;       // period, in sec
    int     rL_norm;     // bool: normalise to light cylinder radius?
    char   *outfile;     // name of output file (NULL means stdout)
    double  s_start;     // starting value of s
    double  s_stop;      // stopping value of s
    int     s_nstep;     // number of s steps
    double  p_start;     // starting value of p
    double  p_stop;      // stopping value of p
    int     p_nstep;     // number of p steps
    int     dipole;      // use a dipole field instead of Deutsch
    int     firstonly;   // for each field line, stop at the first solution
};

void usage();
void parse_cmd_line( int argc, char *argv[], struct opts *o );
void print_col_headers( FILE *f );

int main( int argc, char *argv[] )
{
    // Set up struct for command line options and set default values
    struct opts o;
    o.al_deg    = NAN;
    o.P_sec     = NAN;
    o.ze_deg    = NAN;
    o.rL_norm   = 0;
    o.outfile   = NULL;
    o.s_start   = NAN;
    o.s_stop    = NAN;
    o.s_nstep   = 0;
    o.p_start   = NAN;
    o.p_stop    = NAN;
    o.p_nstep   = 0;
    o.dipole    = 0;
    o.firstonly = 0;

    parse_cmd_line( argc, argv, &o );

    // Set up output file
    FILE *f;
    if (o.outfile == NULL)
        f = stdout;
    else
    {
        f = fopen( o.outfile, "w" );
        if (f == NULL)
        {
            fprintf( stderr, "error: could not open file %s\n", o.outfile );
            exit(EXIT_FAILURE);
        }
    }

    // Set up pulsar
    pulsar psr;

    psr_angle *ra  = NULL;
    psr_angle *dec = NULL;

    psr_angle *al = create_psr_angle_deg( o.al_deg );
    psr_angle *ze = create_psr_angle_deg( o.ze_deg );

    double P = o.P_sec;
    double r = 1e4; /* This will be used later as we move outwards from the
                       pulsar surface */

    set_pulsar( &psr, ra, dec, P, r, al, ze );

    if (o.dipole)
        psr.field_type = DIPOLE;

    // Set scaling factor in case user requests result in units of rL
    double xscale = (o.rL_norm ? 1.0/psr.rL : 1.0);

    // Write the file and column headers
    print_psrg_header( f, argc, argv );
    print_col_headers( f );

    // For calculating the (retarded) line of sight
    point B; // The magnetic field at the emission point
    point V, retarded_LoS;
    point A; // The acceleration vector at the emission point
    psr_angle dph; // The retardation angle
    psr_angle psi; // The polarisation angle

    // Loop over the polar cap
    int open_found; // Set to 1 if an open line is found for a given s
    int linetype;   // either CLOSED_LINE or OPEN_LINE
    point foot_pt, foot_pt_mag;
    point init_pt, emit_pt;
    psr_angle phase;     // The emission phase
    psr_angle ret_phase; // The retarded (observed) phase
    double kappa; // The curvature
    int find_emitpt_result;
    double dist, dist_tmp; // Keep track of distance travelled along field lines

    int s_idx, p_idx;
    double s_deg, p_deg;
    double ds, dp;
    psr_angle s; // The "polar cap distance" to the magnetic pole
    psr_angle p; // The azimuth angle around the polar cap

    for (s_idx = 0; s_idx < o.s_nstep; s_idx++)
    {
        // Convert s_idx to an angle
        ds = (o.s_nstep == 1 ?
                0.0 :
                (o.s_stop - o.s_start)/(o.s_nstep - 1.0));
        s_deg = o.s_start + s_idx*ds;
        set_psr_angle_deg( &s, s_deg );

        // Reset open_found to NOT found
        open_found = 0;

        // Loop around the polar cap in azimuth
        for (p_idx = 0; p_idx < o.p_nstep; p_idx++)
        {
            // Convert p_idx to an angle
            dp = (o.p_nstep == 1 ?
                    0.0 :
                    (o.p_stop - o.p_start)/(o.p_nstep - 1.0));
            p_deg = o.p_start + p_idx*dp;
            set_psr_angle_deg( &p, p_deg );

            // Convert (s,p) into a point in the magnetic frame
            set_point_sph( &foot_pt_mag, psr.r, &s, &p, POINT_SET_ALL );

            // Convert the foot_pt into observer coordinates
            mag_to_obs_frame( &foot_pt_mag, &psr, NULL, &foot_pt );

            // Now check that we're on an open field line
            linetype = get_fieldline_type( &foot_pt, &psr, o.rL_norm, NULL,
                    NULL, NULL );
            if (linetype == CLOSED_LINE)
            {
                continue;
            }
            open_found = 1;

            // Now find the emission points along this line!
            // Start 1 metre above the surface
            Bstep( &foot_pt, &psr, 1.0, DIR_OUTWARD, &init_pt );
            set_point_xyz( &init_pt, init_pt.x[0],
                                     init_pt.x[1],
                                     init_pt.x[2],
                                     POINT_SET_ALL );
            dist = 0.0; // Start distance tracker

            while (1)
            {
                // Climb up the field line to find the next emit_pt
                find_emitpt_result = find_next_line_emission_point( &psr,
                        &init_pt, DIR_OUTWARD, &emit_pt, &dist_tmp,
                        NULL );
                dist += dist_tmp;

                // If no point was found, exit the loop
                if (find_emitpt_result != EMIT_PT_FOUND)
                {
                    break;
                }

                // Calculate the (retarded) phase at which the emission would
                // be seen. First, set the V to the velocity vector. While
                // we're at it, get the acceleration vector and the curvature.
                calc_fields( &emit_pt, &psr, SPEED_OF_LIGHT, &B, &V,
                             NULL, &A, NULL, NULL );
                calc_retardation( &emit_pt, &psr, &V, &dph, &retarded_LoS );
                kappa = calc_curvature( &V, &A );

                // Now, the observed phase is the negative of the azimuthal
                // angle of the retarded line of sight
                set_point_xyz( &V, V.x[0], V.x[1], V.x[2], POINT_SET_PH );
                if (psr.spin == SPIN_POS)
                    set_psr_angle_deg( &phase, -V.ph.deg );
                else
                    copy_psr_angle( &(V.ph), &phase );
                set_psr_angle_deg( &ret_phase, -(retarded_LoS.ph.deg) );

                // Calculate the observed polarisation angle at emit_pt
                accel_to_pol_angle( &psr, &A, &phase, &psi );

                // Print out results!
                /*
                fprintf( f, "%.15e %.15e %.15e %.15e "
                            "%.15e %.15e %.15e %.15e\n",
                            s_deg, p_deg,
                            emit_pt.x[0] * xscale,
                            emit_pt.x[1] * xscale,
                            emit_pt.x[2] * xscale,
                            ret_phase.deg, psi.deg, kappa );
                */
                fprintf( f, "%.15e %.15e %.15e %.15e "
                            "%.15e %.15e %.15e %.15e "
                            "%.15e %.15e %.15e %.15e "
                            "%.15e %.15e %.15e %.15e "
                            "%.15e %.15e %.15e %.15e "
                            "%.15e\n",
                            s_deg, p_deg,
                            emit_pt.x[0] * xscale,
                            emit_pt.x[1] * xscale,
                            emit_pt.x[2] * xscale,
                            phase.deg, psi.deg, kappa,
                            B.x[0], B.x[1], B.x[2], B.r,
                            V.x[0], V.x[1], V.x[2], V.r,
                            A.x[0], A.x[1], A.x[2], A.r,
                            dist );

                // Set the emission point to the new initial point, go another
                // 1 km along, and then try to find the next emit_pt
                // It seems that anything much short than a 1 km jump
                // results in the next point being found in the same area. The
                // size of the volume that converges appears to be larger than
                // expected. This is probably a bug, but for now I'll just
                // move along 1 km before trying again.
                if (o.firstonly)
                    break;
                else
                {
                    copy_point( &emit_pt, &init_pt );
                    Bstep( &init_pt, &psr, 1000.0, DIR_OUTWARD, &init_pt );
                    set_point_xyz( &init_pt, init_pt.x[0],
                                             init_pt.x[1],
                                             init_pt.x[2],
                                             POINT_SET_ALL );
                }
            }

            // If s = 0°, then no need to do any more values of p
            if (s_deg == 0.0)
                break;
        }

        // Check to see if any open lines were found. If not, exit the loop
        if (open_found == 0)
            break;
    }


    // Clean up
    destroy_psr_angle( ra  );
    destroy_psr_angle( dec );
    destroy_psr_angle( al  );
    destroy_psr_angle( ze  );

    free( o.outfile );

    if (o.outfile != NULL)
        fclose( f );

    return 0;
}

void usage()
{
    printf( "usage: psr_visiblepoints [OPTIONS]\n\n" );
    printf( "REQUIRED OPTIONS:\n" );
    printf( "  -a  alpha    The angle between the rotation and magetic axes "
                           "in degrees (required)\n" );
    printf( "  -P  period   The rotation period of the pulsar, in seconds "
                           "(required)\n" );
    printf( "  -p  p[:P[:n]]   The azimuth relative to the magnetic axis, "
                           "in degrees. The range is from p to P with n "
                           "steps.\n"
            "                 p      ==> p:p:1\n"
            "                 p:P    ==> p:P:2\n" );
    printf( "  -s  s[:S[:n]]   The angular distance from the magnetic axis, "
                           "in degrees. The range is from s to S with n "
                           "steps.\n"
            "                 s      ==> s:s:1\n"
            "                 s:S    ==> s:S:2\n" );
    printf( "  -t  step     Step size for moving along magnetic field lines, "
                           "as a fraction of the light cylinder radius "
                           "(default: 0.01)\n" );
    printf( "  -z  zeta     The angle between the rotation axis and the line "
                           "of sight in degrees (required)\n" );
    printf( "\nOTHER OPTIONS:\n" );
    printf( "  -1           Only report one (i.e. the first) solution for "
                           "each field line (default: off)\n" );
    printf( "  -d           Use a dipole field instead of the default "
                           "Deutsch field\n" );
    printf( "  -h           Display this help and exit\n" );
    printf( "  -L           Normalise distances to light cylinder radius\n" );
    printf( "  -o  outfile  The name of the output file to write to. If not "
                           "set, output will be written to stdout.\n" );
}


void parse_cmd_line( int argc, char *argv[], struct opts *o )
{
    // Collect the command line arguments
    int c;
    while ((c = getopt( argc, argv, "1a:dhLo:p:P:s:S:z:")) != -1)
    {
        switch (c)
        {
            case '1':
                o->firstonly = 1;
                break;
            case 'a':
                o->al_deg = atof(optarg);
                break;
            case 'd':
                o->dipole = 1;
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'L':
                o->rL_norm = 1;
                break;
            case 'o':
                o->outfile = strdup(optarg);
                break;
            case 'p':
                parse_range( optarg, &(o->p_start),
                                     &(o->p_stop),
                                     &(o->p_nstep) );
                break;
            case 'P':
                o->P_sec = atof(optarg);
                break;
            case 's':
                parse_range( optarg, &(o->s_start),
                                     &(o->s_stop),
                                     &(o->s_nstep) );
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
    if (isnan(o->al_deg) || isnan(o->P_sec) || isnan(o->ze_deg))
    {
        fprintf( stderr, "error: -a, -P, and -z options required"
                         "\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    if (isnan(o->s_start) || isnan(o->p_start))
    {
        fprintf( stderr, "error: -p and -s options required\n" );
        usage();
        exit(EXIT_FAILURE);
    }
}


void print_col_headers( FILE *f )
/* The final output includes:
 *   1) the polar coordinates of the footpoint of the magnetic field line
 *   2) the Cartesian coordinates of the emission point
 *   3) the retarded (observed) emission phase
 *   4) the observed polarisation angle
 *   5) the curvature of the particle's trajectory at the emission point
 */
{
    // Print out a line to file handle f
    fprintf( f, "#  s_deg  p_deg  x  y  z  phase_deg  polangle_deg  curvature  "
                "Bx  By  Bz  B  "
                "Vx  Vy  Vz  V  "
                "Ax  Ay  Az  A  "
                "dist\n" );
}


