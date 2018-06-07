#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include "psrgeom.h"

#define MAX_NPOINTS 1000000
#define NBINS 1024

static int creation_rate; // particles per time step
static int npoints;
static photon pns[MAX_NPOINTS];
static double profile[NBINS];

static pulsar psr;
static gamma_distr gd;

static double freq_lo;
static double freq_hi;

static char status_str[256];

// Mouse-related variables
static int mouse_old_x;
static int mouse_old_y;
static int button_down;

// Determine which feature is nearer
static int nearest_feature;
static int selected_feature;

static int highlight;

static int repopulate; // bool message, whether to call init_population()

static double P4_scale; // Degrees (of circular P4 arrow) per second

static double t;
static double tstep;

static double max_power;

static int playforward;


enum
{
    NO_FEATURE,
    ALPHA_LINE,
    ZETA_LINE,
    CSL_CIRCLE,
    SPARK_CIRCLE,
    P4_ARROW,
    GAMMA_MEAN,
    GAMMA_STD,
    GAMMA_LO,
    GAMMA_HI,
    PERIOD_SLIDER
};

typedef struct scene_t
{
    point camera;
    point lookat;
    point up;
    int dim;
    double near_clip;
    double far_clip;
} scene;

typedef struct view_t
{
    int x, y, w, h;
    double FoV;
    double aspect_ratio;
    int scene_num;
    double left, right, bottom, top;
    int border; // boolean: border or not?
} view;

typedef struct graph_t
{
    double lv, rv, bv, tv; // The graph's boundary in view coords
    double xmin, xmax;
    double ymin, ymax;
    int logx, logy;
    double xtics, ytics;
    double xmtics, ymtics;
    char xlabel[256], ylabel[256];
} graph;

#define NVIEWS   11
#define NSCENES  8
static view  views[NVIEWS];
static scene scenes[NSCENES];
static graph gamma_graph;
static graph period_graph;

static int active_view;

enum
{
    SCENE_NONE       = -1,
    SCENE_ANGLES     = 0,
    SCENE_FOOTPTS    = 1,
    SCENE_BEAM       = 2,
    SCENE_PROFILE    = 3,
    SCENE_GAMMA      = 4,
    SCENE_FIELDLINES = 5,
    SCENE_STATUS     = 6,
    SCENE_TELESCOPE  = 7
};

enum
{
    VIEW_NONE         = -1,
    VIEW_THUMBNAIL_1  = 0,
    VIEW_THUMBNAIL_2  = 1,
    VIEW_THUMBNAIL_3  = 2,
    VIEW_THUMBNAIL_4  = 3,
    VIEW_THUMBNAIL_5  = 4,
    VIEW_THUMBNAIL_6  = 5,
    VIEW_THUMBNAIL_7  = 6,
    VIEW_THUMBNAIL_8  = 7,
    VIEW_LEFT_MAIN    = 8,
    VIEW_RIGHT_MAIN   = 9,
    VIEW_STATUS       = 10
};

enum
{
    MOUSE_LEFT_BUTTON = 0,
    MOUSE_MIDDLE_BUTTON = 1,
    MOUSE_RIGHT_BUTTON = 2,
    MOUSE_SCROLL_UP = 3,
    MOUSE_SCROLL_DOWN = 4
};

static int W, H;

void print_str( char *s, double x, double y, void *font )
{
    if (s == NULL)
    {
        fprintf( stderr, "error: print_str: char *s is NULL\n" );
        exit(EXIT_FAILURE);
    }

    glRasterPos2d( x, y );
    while (*s)
    {
        glutBitmapCharacter(font, *s);
        s++;
    }
}


void clear_profile()
{
    int n;
    for (n = 0; n < NBINS; n++)
        profile[n] = 0.0;
}


int which_view( int x, int y )
{
    view *vw;
    int v;
    for (v = 0; v < NVIEWS; v++)
    {
        vw = &views[v];
        if ((  x >= vw->x) && (  x < vw->x + vw->w) &&
            (H-y >= vw->y) && (H-y < vw->y + vw->h))
        {
            return v;
        }
    }

    return SCENE_NONE;
}


void respawn_particle( photon *pn )
{
    random_spark_footpoint( &(pn->source), NULL, &psr, t );
}


void update_powers()
{
    int n;
    for (n = 0; n < npoints; n++)
    {
        pns[n].power = binary_power_single( &gd, freq_lo, freq_hi,
                pns[n].curvature );
    }
}


void add_to_profile_bins()
/* Iterate over all existing points and add their powers (assumed to be
 * already calculated) to the appropriate profile bin.
 */
{
    int n;
    double retard_time;
    int bin;

    // Calculate the phase and the line of sight
    psr_angle phase;
    point LoS;

    psr_angle impact_angle;

    set_psr_angle_deg( &phase, 360.0*t/psr.P );
    line_of_sight( &psr, &phase, &LoS );

    for (n = 0; n < npoints; n++)
    {
        // Check to see if the particle is pointing in the right direction
        set_psr_angle_cos( &impact_angle,
                           pns[n].V.x[0]*LoS.x[0] +
                           pns[n].V.x[1]*LoS.x[1] +
                           pns[n].V.x[2]*LoS.x[2] );

        // Only worry about particles where θ <= 1/γ
        if (impact_angle.rad > 1.0/gd.mean)
            continue;

        // Calculate the retardation time
        retard_time = (pns[n].source.x[0]*LoS.x[0] +
                       pns[n].source.x[1]*LoS.x[1] +
                       pns[n].source.x[2]*LoS.x[2]) / SPEED_OF_LIGHT;

        // Calculate the profile bin (bin 0 is the fiducial point)
        bin = (int)round( (t - retard_time)/psr.P );
        while (bin <  0    )  bin += NBINS;
        while (bin >= NBINS)  bin -= NBINS;

        // Add the power to the profile bin
        profile[bin] += pns[n].power;
    }
}


int advance_particles_once()
/* Move each existing particle along a bit.
 * Returns 1 is at least one particle was killed and respawned,
 * 0 otherwise.
 */
{
    int respawned = 0;
    int n;
    double V_dot_r;
    for (n = 0; n < npoints; n++)
    {
        // Calculate the magnetic and velocity vectors at the
        // current point
        calc_fields( &pns[n].source, &psr, SPEED_OF_LIGHT,
                &pns[n].B, &pns[n].V, NULL, NULL, NULL, NULL );
        set_point_xyz( &pns[n].V,
                       pns[n].V.x[0], pns[n].V.x[1], pns[n].V.x[2],
                       POINT_SET_ALL );

        // Take a step forward in time
        traj_step( &(pns[n].source), &psr, tstep, DIR_OUTWARD,
                &(pns[n].source), &pns[n].B, &pns[n].V );

        set_point_xyz( &(pns[n].source),
                       pns[n].source.x[0],
                       pns[n].source.x[1],
                       pns[n].source.x[2],
                       POINT_SET_ALL );

        // See if we have to kill any particles off and respawn them
        // Are any about to crash into the pulsar?
        // Are any about to escape the light cylinder?
        // For the former, need to check that the velocity is pointing
        // inward as well
        V_dot_r = pns[n].V.x[0] * pns[n].source.x[0] +
                  pns[n].V.x[1] * pns[n].source.x[1] +
                  pns[n].V.x[2] * pns[n].source.x[2];
        if (((pns[n].source.r <= SPEED_OF_LIGHT*tstep) &&
             (V_dot_r < 0.0)) ||
            (pns[n].source.rhosq >= psr.rL2))
        {
            respawn_particle( &pns[n] );
            respawned = 1;
        }

        // Calculate photon properties
        emit_avg_pulsar_photon( &psr, &pns[n].source, freq_lo, freq_hi,
                                &gd, &pns[n] );

        // Keep track of the particle with the most power, for normalisation
        // purposes
        if (pns[n].power > max_power)
            max_power = pns[n].power;
    }

    // Actually change the time!
    t += tstep;

    return respawned;
}

void init_particles()
/* Populate the pns array until one of the created particles is killed off
 */
{
    npoints = 0;
    int n;
    t = 0.0;
    int respawned;

    while (1)
    {
        // Advance existing particles one time step, and check to see if any
        // were respawned. If so, declare that we have enough particles, and
        // stop the loop.
        respawned = advance_particles_once();
        if (respawned)   break;

        // Make sure our array is still big enough to accommodate a new batch
        if (npoints + creation_rate > MAX_NPOINTS)
        {
            fprintf( stderr, "warning: init_particles: exceeded MAX_NPOINTS. "
                    "Try lowering the creation rate.\n" );
            break;
        }

        // Create a new batch (size is determined by creation_rate)
        for (n = npoints; n < npoints + creation_rate; n++)
            respawn_particle( &pns[n] );
        npoints += creation_rate;
    }
}

void draw_2D_circle( double radius, double xc, double yc )
{
    glBegin(GL_LINE_LOOP);

    int i;
    for (i = 0; i < 360; i++)
    {
        double rad = i*DEG2RAD;
        glVertex2d( xc + cos(rad)*radius, yc + sin(rad)*radius );
    }

    glEnd();
}

void draw_2D_arc( double radius, double xc, double yc, double start_deg,
        double end_deg )
/* Draws a circle arc, counterclockwise from START_DEG to END_DEG */
{
    while (end_deg < start_deg)  end_deg += 360.0;

    // Calculate the biggest interval smaller than a degree
    double L = end_deg - start_deg;
    int nintervals = (int)ceil(L);
    double dl = L / (double)nintervals;

    glBegin(GL_LINE_STRIP);

    int i;
    for (i = 0; i < nintervals; i++)
    {
        double rad = ((double)i * dl + start_deg)*DEG2RAD;
        glVertex2d( xc + cos(rad)*radius, yc + sin(rad)*radius );
    }

    glEnd();
}

void draw_3D_circle( double radius, double xc, double yc, double z )
/* Draws a circle parallel to the z-plane */
{
   glBegin(GL_LINE_LOOP);

   int i;
   for (i = 0; i < 360; i++)
   {
      double rad = i*DEG2RAD;
      glVertex3d( xc + cos(rad)*radius, yc + sin(rad)*radius, z );
   }

   glEnd();
}

void reshape_views()
{
    int status_bar_height = 20;

    // Thumbnails
    int ntn = 8;      // Number of thumnails
    int tn_s = W/ntn; // Size of thumbnails

    int tn, v;
    for (tn = 0; tn < ntn; tn++)
    {
        v = tn + VIEW_THUMBNAIL_1; // This assumes the view numbers of the
                                   // thumbnails are contiguous
        views[v].x = tn*tn_s;
        views[v].y = H-tn_s;
        views[v].h = tn_s;
        views[v].w = tn_s;
    }

    views[VIEW_LEFT_MAIN].x = 0;
    views[VIEW_LEFT_MAIN].y = status_bar_height;
    views[VIEW_LEFT_MAIN].h = H - tn_s - status_bar_height;
    views[VIEW_LEFT_MAIN].w = W/2;

    views[VIEW_RIGHT_MAIN].x = W/2;
    views[VIEW_RIGHT_MAIN].y = status_bar_height;
    views[VIEW_RIGHT_MAIN].h = H - tn_s - status_bar_height;
    views[VIEW_RIGHT_MAIN].w = W/2;

    views[VIEW_STATUS].x = 0;
    views[VIEW_STATUS].y = 0;
    views[VIEW_STATUS].h = status_bar_height;
    views[VIEW_STATUS].w = W;

}

void init_scenes()
{
    /* Set up initial cameras */
    // (Some useful angles for initialisation)
    psr_angle z_angle; // "zero" angle
    psr_angle r_angle; // right angle
    psr_angle n_angle; // negative right angle
    set_psr_angle_deg( &z_angle,   0.0 );
    set_psr_angle_deg( &r_angle,  90.0 );
    set_psr_angle_deg( &n_angle, -90.0 );

    // SCENE_FIELDLINES
    scene *scn = &scenes[SCENE_FIELDLINES];
    set_point_xyz( &scn->camera, 0.0, -2.0*psr.rL, 0.0, POINT_SET_ALL );
    set_point_xyz( &scn->lookat, 0.0, 0.0,         0.0, POINT_SET_ALL );
    set_point_xyz( &scn->up,     0.0, 0.0,         1.0, POINT_SET_ALL );

}

void init_views()
{
    // Assign scenes to views
    views[VIEW_THUMBNAIL_1 ].scene_num = SCENE_FOOTPTS;
    views[VIEW_THUMBNAIL_2 ].scene_num = SCENE_FIELDLINES;
    views[VIEW_THUMBNAIL_3 ].scene_num = SCENE_ANGLES;
    views[VIEW_THUMBNAIL_4 ].scene_num = SCENE_GAMMA;
    views[VIEW_THUMBNAIL_5 ].scene_num = SCENE_TELESCOPE;
    views[VIEW_THUMBNAIL_6 ].scene_num = SCENE_BEAM;
    views[VIEW_THUMBNAIL_7 ].scene_num = SCENE_PROFILE;
    views[VIEW_THUMBNAIL_8 ].scene_num = SCENE_NONE;
    views[VIEW_LEFT_MAIN   ].scene_num = SCENE_FOOTPTS;
    views[VIEW_RIGHT_MAIN  ].scene_num = SCENE_BEAM;
    views[VIEW_STATUS      ].scene_num = SCENE_STATUS;

    // Turn on borders
    views[VIEW_THUMBNAIL_1 ].border = 1;
    views[VIEW_THUMBNAIL_2 ].border = 1;
    views[VIEW_THUMBNAIL_3 ].border = 1;
    views[VIEW_THUMBNAIL_4 ].border = 1;
    views[VIEW_THUMBNAIL_5 ].border = 1;
    views[VIEW_THUMBNAIL_6 ].border = 1;
    views[VIEW_THUMBNAIL_7 ].border = 1;
    views[VIEW_THUMBNAIL_8 ].border = 1;
    views[VIEW_LEFT_MAIN   ].border = 0;
    views[VIEW_RIGHT_MAIN  ].border = 0;
    views[VIEW_STATUS      ].border = 1;

    reshape_views();
}

void set_scene_properties()
{
    scenes[SCENE_FOOTPTS   ].dim = 2;
    scenes[SCENE_ANGLES    ].dim = 0;
    scenes[SCENE_GAMMA     ].dim = 0;
    scenes[SCENE_PROFILE   ].dim = 2;
    scenes[SCENE_BEAM      ].dim = 2;
    scenes[SCENE_STATUS    ].dim = 0;
    scenes[SCENE_FIELDLINES].dim = 3;

    scenes[SCENE_FIELDLINES].near_clip = 0.0;
    scenes[SCENE_FIELDLINES].far_clip  = 2.0*psr.rL;

}

void set_view_properties()
{
    view *vw;

    int v;
    for (v = 0; v < NVIEWS; v++)
    {
        vw = &views[v];
        switch (vw->scene_num)
        {
            case SCENE_FOOTPTS:
                vw->bottom = -(psr.csl.S.deg + psr.csl.s.deg)*1.4;
                vw->top    = -vw->bottom;
                vw->left   =  vw->bottom*(double)vw->w/(double)vw->h;
                vw->right  = -vw->left;
                break;
            case SCENE_FIELDLINES:
                vw->FoV = 60.0;
                vw->aspect_ratio = (GLfloat) vw->w/(GLfloat) vw->h;
                break;
            case SCENE_ANGLES:
                vw->left   = -0.1;
                vw->right  = 1.5;
                vw->bottom = 0.0;
                vw->top    = (vw->right - vw->left)*(double)vw->h/(double)vw->w;
                break;
            case SCENE_GAMMA:
                vw->left   = 0.0;
                vw->right  = 1.0;
                vw->bottom = 0.0;
                vw->top    = 1.0;
                break;
            case SCENE_BEAM:
                vw->left   = -100.0;
                vw->right  =  100.0;
                vw->bottom = -100.0;
                vw->top    =  100.0;
                break;
            case SCENE_STATUS:
                vw->left   = 0.0;
                vw->right  = 1.0;
                vw->bottom = 0.0;
                vw->top    = 1.0;
                break;
        }
    }
}

void apply_2D_camera( int view_num )
{
    view *vw = &views[view_num];

    glMatrixMode (GL_PROJECTION);

    glLoadIdentity ();
    gluOrtho2D( vw->left, vw->right, vw->bottom, vw->top );
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
}

void apply_3D_camera( int view_num )
{
    view *vw = &views[view_num];
    scene *scn = &scenes[vw->scene_num];
    point *cam = &scn->camera;

    glMatrixMode (GL_PROJECTION);

    glLoadIdentity ();
    gluPerspective( vw->FoV, vw->aspect_ratio, scn->near_clip, scn->far_clip );
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
    gluLookAt (cam->x[0], cam->x[1], cam->x[2],
               0.0, 0.0, 0.0,
               0.0, 0.0, 1.0);
}

void screen2world( int x, int y, double *xw, double *yw,
        psr_angle *th, psr_angle *ph, int view_num )
/* Convert from screen coordinates to world coordinates. The polar coordinates
 * th and ph assume that the x and y coordinates are in units of degrees, and
 * that the reference axis for rotation is the positive y axis, with positive
 * angles on the clockwise side of it.
 */
{
    if (view_num == VIEW_NONE)
    {
        fprintf( stderr, "error: screen2world: invalid view number\n" );
        exit(EXIT_FAILURE);
    }

    double x_tmp, y_tmp;

    view *vw = &views[view_num];

    x_tmp = (double)(  x - vw->x)/(double)vw->w*(vw->right - vw->left)
        + vw->left;
    y_tmp = (double)(H-y - vw->y)/(double)vw->h*(vw->top - vw->bottom)
        + vw->bottom;

    if (xw)
        *xw = x_tmp;
    if (yw)
        *yw = y_tmp;
    if (ph)
        set_psr_angle_rad( ph, atan2( x_tmp, y_tmp ) );
    if (th)
        set_psr_angle_deg( th, hypot( x_tmp, y_tmp ) );
}

void screen2csl( int x, int y, psr_angle *S, psr_angle *s, int *n,
        int view_num )
/* Convert from screen coordinates to "carousel" polar coordinates, which
 * means the angle S from the magnetic axis (assumed to be at the origin of
 * the world coordinates); the angle s from the centre of the nearest spark;
 * and the spark number n.
 * It is assumed that at time t=0, spark number 0 is on the negative y-axis.
 */
{
    int n_tmp; // Doing it this way allows the caller to pass in n = NULL
    double xw, yw;
    psr_angle th, ph;
    screen2world( x, y, &xw, &yw, &th, &ph, view_num );

    // Taking carousel rotation and the 90° shift into account
    double P4 = psr.csl.P4;
    double ph_t = 360.0*t/P4; // The amount of rotation since t=0, in deg
    psr_angle ph0;
    set_psr_angle_deg( &ph0, (180.0 - ph.deg) - ph_t );

    // Calculate which spark we're "in"
    int N = psr.csl.n;
    if (N == 0)
        n_tmp = 0;
    else
    {
        n_tmp = (int)floor((double)N * ph0.deg / 360.0 + 0.5);
        while (n_tmp <  0)  n_tmp += N;
        while (n_tmp >= N)  n_tmp -= N;
    }

    if (S)
        copy_psr_angle( &th, S );
    if (n)
        *n = n_tmp;
    if (s)
    {
        if (N == 0)
        {
            set_psr_angle_deg( s, fabs( psr.csl.S.deg - th.deg ) );
        }
        else
        {
            // Calculate how far away from the nearest spark we are
            double xs, ys; // The coordinates of the spark centre
            double dx, dy;

            psr_angle ph_s;
            set_psr_angle_deg( &ph_s, (double)n_tmp*360.0/(double)N +
                    ph_t - 90.0 );

            xs = psr.csl.S.deg * ph_s.cos;
            ys = psr.csl.S.deg * ph_s.sin;

            dx = xw - xs;
            dy = yw - ys;

            set_psr_angle_deg( s, hypot( dx, dy ) );
        }
    }
}


void screen2graph( int x, int y, double *xg, double *yg, graph *gr,
        int view_num )
/* Convert from screen coordinates to graph-world coordinates */
{
    double xv, yv; // "View" coordinates
    screen2world( x, y, &xv, &yv, NULL, NULL, view_num );

    // Convert to "plot" coords
    double xp, yp;
    xp = (xv - gr->lv) / (gr->rv - gr->lv);
    yp = (yv - gr->bv) / (gr->tv - gr->bv);

    // And finally to graph-world coords
    if (!gr->logx)
        *xg = xp*(gr->xmax - gr->xmin) + gr->xmin;
    else
        *xg = gr->xmin*pow( gr->xmax / gr->xmin, xp );

    if (!gr->logy)
        *yg = yp*(gr->ymax - gr->ymin) + gr->ymin;
    else
        *yg = gr->ymin*pow( gr->ymax / gr->ymin, yp );
}


void adjust_tstep()
{
    tstep = 0.005*psr.P/(2.0*PI);
}


void init(void) 
{
    // Set a white background
    glClearColor (1.0, 1.0, 1.0, 0.0);
    glShadeModel (GL_FLAT);

    // Set up a default pulsar
    psr_angle ze;
    psr_angle al;
    set_psr_angle_deg( &ze, 5.0 );
    set_psr_angle_deg( &al, 10.0 );
    double P = 1.0;
    double rp = 1.0e4;
    set_pulsar( &psr, NULL, NULL, P, rp, &al, &ze );

    // Set up default carousel
    int nsparks = 7;
    psr_angle s, S;
    set_psr_angle_deg( &S, 1.0 );
    set_psr_angle_deg( &s, 0.1 );
    double P4_sec = -10.0;
    set_pulsar_carousel( &psr, nsparks, &s, &S, GAUSSIAN, P4_sec );
    P4_scale = 4.0;

    // Set up the default gamma distribution
    gd.type  = NORMAL;
    gd.mean  = 500.0;
    gd.std   = 10.0;
    gd.idx   = -6.2;
    gd.g_min = 450.0;
    gd.g_max = 550.0;

    // Set up telescope properties
    freq_lo = 300.0e6;
    freq_hi = 350.0e6;

    // No status message to start with
    strcpy( status_str, "" );

    // Set up the default gamma distribution graph
    gamma_graph.lv     = 0.1;
    gamma_graph.rv     = 0.95;
    gamma_graph.bv     = 0.6;
    gamma_graph.tv     = 0.8;
    gamma_graph.xmin   = 400.0;
    gamma_graph.xmax   = 600.0;
    gamma_graph.ymin   = 0.0;
    gamma_graph.ymax   = 1.0;
    gamma_graph.logx   = 0;
    gamma_graph.logy   = 0;
    gamma_graph.xtics  = 50.0;
    gamma_graph.ytics  = 0.0;
    gamma_graph.xmtics = 10.0;
    gamma_graph.ymtics = 0.0;;
    strcpy( gamma_graph.xlabel, "γ" );
    strcpy( gamma_graph.ylabel, "" );

    // Set up the default period graph
    period_graph.lv     = 0.1;
    period_graph.rv     = 1.4;
    period_graph.bv     = 0.1;
    period_graph.tv     = 0.12;
    period_graph.xmin   = 0.001;
    period_graph.xmax   = 100.0;
    period_graph.ymin   = -1.0;
    period_graph.ymax   = 1.0;
    period_graph.logx   = 1;
    period_graph.logy   = 0;
    period_graph.xtics  = 10.0;
    period_graph.ytics  = 0.0;
    period_graph.xmtics = 10.0;
    period_graph.ymtics = 0.0;;
    strcpy( period_graph.xlabel, "P (s)" );
    strcpy( period_graph.ylabel, "" );

    // Set the rest of the scene properties based on the cameras
    init_views();
    init_scenes();
    set_scene_properties();
    set_view_properties();
    repopulate = 0;

    // Set the (initial) creation rate of particles
    adjust_tstep();
    creation_rate = 100; // per time step
    playforward = 0;

    // Clear the profile and pulsestack
    clear_profile();

    max_power = DBL_MIN; // A dummy, non-zero value to get things started

    // Unset nearest_feature
    nearest_feature = NO_FEATURE;
    selected_feature = NO_FEATURE;
    highlight = VIEW_NONE;

    // Unset active_view
    active_view = SCENE_NONE;
}


void draw_graph_frame( graph *gr, int draw_text )
{
    glPushMatrix();

    // Draw the axes (graph coords assumed)
    glColor3f( 0.0, 0.0, 0.0 );
    glBegin( GL_LINES );

    glVertex2d( 0.0, 0.0 );
    glVertex2d( 0.0, 1.0 );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.0, 0.0 );

    glEnd();

    // Now go to the graph's world coordinates
    glScaled( 1.0/(gr->xmax - gr->xmin), 1.0/(gr->ymax - gr->ymin), 1.0 );
    glTranslated( -gr->xmin, -gr->ymin, 0.0 );

    glBegin( GL_LINES );
    char ticlabel[32];

    // Calculate the tics
    double xtic_min, ytic_min;
    double xtic_max, ytic_max;
    double xtic, ytic;
    if (gr->xtics > 0.0)
    {
        xtic_min =  ceil( gr->xmin / gr->xtics ) * gr->xtics;
        xtic_max = floor( gr->xmax / gr->xtics ) * gr->xtics;
        for (xtic = xtic_min; xtic <= xtic_max; xtic += gr->xtics)
        {
            glVertex2d( xtic, gr->ymin );
            glVertex2d( xtic, 0.05*(gr->ymax - gr->ymin) + gr->ymin );
        }
    }
    if (gr->ytics > 0.0)
    {
        ytic_min =  ceil( gr->ymin / gr->ytics ) * gr->ytics;
        ytic_max = floor( gr->ymax / gr->ytics ) * gr->ytics;
        for (ytic = ytic_min; ytic <= ytic_max; ytic += gr->ytics)
        {
            glVertex2d( 0.0,  ytic );
            glVertex2d( 0.05, ytic );
        }
    }

    glEnd();

    if (draw_text)
    {
        for (xtic = xtic_min; xtic <= xtic_max; xtic += gr->xtics)
        {
            sprintf( ticlabel, "%d", (int)xtic );
            print_str( ticlabel, xtic-2.0*strlen(ticlabel), -0.20, GLUT_BITMAP_HELVETICA_12 );
        }
    }

    glPopMatrix();
}


void draw_border( int view_num, double linewidth )
{
    view *vw = &views[view_num];
    //if (highlight == view_num)
        glColor3f( 0.0, 0.0, 0.0 );
    //else
    //    glColor3f( 1.0, 0.5, 0.0 );
    glLineWidth( linewidth );
    glBegin( GL_LINE_LOOP );
        glVertex2d( vw->left , vw->bottom );
        glVertex2d( vw->right, vw->bottom );
        glVertex2d( vw->right, vw->top    );
        glVertex2d( vw->left , vw->top    );
    glEnd();
    glLineWidth( 1.0 );
}


void display_status( int view_num )
{
    apply_2D_camera( view_num );

    // Draw a border around the perimeter of the view
    view *vw = &views[view_num];
    if (vw->border)
    {
        double borderwidth = 2.0;
        glPushMatrix();
        draw_border( view_num, borderwidth );
        glPopMatrix();
    }

    print_str( status_str, 0.01, 0.2, GLUT_BITMAP_HELVETICA_18 );

    char time_str[64];
    if (playforward)
        sprintf( time_str, "t = %.3f    Push space to pause", t );
    else
        sprintf( time_str, "t = %.3f    Push space to play", t );
    print_str( time_str, 0.51, 0.2, GLUT_BITMAP_HELVETICA_18 );
}

void display_footpts( int view_num )
{
    apply_2D_camera( view_num );

    // Draw a border around the perimeter of the view
    view *vw = &views[view_num];
    if (vw->border)
    {
        double borderwidth = 2.0;
        glPushMatrix();
        draw_border( view_num, borderwidth );
        glPopMatrix();
    }

    // Draw lines of colatitude
    glPushMatrix();
    glColor3f( 0.65, 0.65, 0.65 );
    double dr = 0.5;
    int r;
    for (r = 1; r <= 5; r++)
        draw_2D_circle( r*dr, 0.0, 0.0 );

    // Draw lines of longitude
    int nl = 12;
    double max_r = 10.0;
    int l;
    glBegin( GL_LINES );
    for( l = 0; l < nl; l++)
    {
        glVertex2d( 0.0, 0.0 );
        glVertex2d( max_r*cos(l*2.0*PI/nl),
                    max_r*sin(l*2.0*PI/nl) );
    }
    glEnd();
    glPopMatrix();

    // Draw circles representing sparks
    if (nearest_feature == CSL_CIRCLE)
        glColor3f( 1.0, 0.0, 0.0 );
    else
        glColor3f( 0.0, 0.0, 1.0 );

    glPushMatrix();
    draw_2D_circle( psr.csl.S.deg, 0.0, 0.0 );
    glPopMatrix();

    if (nearest_feature == SPARK_CIRCLE)
        glColor3f( 1.0, 0.0, 0.0 );
    else
        glColor3f( 0.0, 0.0, 1.0 );

    if (psr.csl.n == 0)
    {
        glPushMatrix();
        draw_2D_circle( psr.csl.S.deg + psr.csl.s.deg, 0.0, 0.0 );
        draw_2D_circle( psr.csl.S.deg - psr.csl.s.deg, 0.0, 0.0 );
        glPopMatrix();
    }
    else
    {
        int n;
        for (n = 0; n < psr.csl.n; n++)
        {
            glPushMatrix();
            glRotated( 360.0*t/psr.csl.P4, 0.0, 0.0, 1.0 );
            glRotated( (double)n * 360.0 / (double)psr.csl.n, 0.0, 0.0, 1.0 );
            glTranslated( 0.0, -psr.csl.S.deg, 0.0 );
            draw_2D_circle( psr.csl.s.deg, 0.0, 0.0 );
            glPopMatrix();
        }
    }

    // Draw P4 arrow
    glPushMatrix();
    glColor3f( 0.0, 0.75, 0.0 );
    double R = (vw->top - vw->bottom)*13.0/28.0;
    double P4_L = psr.csl.P4*P4_scale;
    if (psr.csl.P4 > 0.0)
        draw_2D_arc( R, 0.0, 0.0, -90.0, P4_L - 90.0 );
    else
    {
        draw_2D_arc( R, 0.0, 0.0, P4_L - 90.0, -90.0 );
    }
    glPopMatrix();

    // The P4 arrow head
    glPushMatrix();
    glRotated( P4_L, 0.0, 0.0, 1.0 );
    glTranslated( 0.0, -R, 0.0 );
    if (psr.csl.P4 < 0.0)
        glRotated( 180.0, 0.0, 0.0, 1.0 );
    glBegin( GL_LINES );
        glVertex2d( 0.0, 0.0 );
        glVertex2d( -0.1*R, -0.05*R );
        glVertex2d( 0.0, 0.0 );
        glVertex2d( -0.1*R, 0.05*R );
    glEnd();
    glPointSize( 3.0 );
    if (nearest_feature == P4_ARROW)
        glColor3f( 1.0, 0.0, 0.0 );
    else
        glColor3f( 0.0, 0.65, 0.0 );
    glBegin( GL_POINTS );
        glVertex2d( 0.0, 0.0 );
    glEnd();
    glPopMatrix();

    // The controls for changing the number of sparks
    glPushMatrix();
    glColor3f( 0.0, 0.65, 0.0 );
    if (view_num == VIEW_LEFT_MAIN)
        print_str( "-/+", 0.0, 0.0, GLUT_BITMAP_TIMES_ROMAN_24 );

    glPopMatrix();
}


void display_fieldlines( int view_num )
{
    apply_3D_camera( view_num );

    // Draw a border around the perimeter of the view
    view *vw = &views[view_num];
    if (vw->border)
    {
        double borderwidth = 2.0;
        glPushMatrix();
        draw_border( view_num, borderwidth );
        glPopMatrix();
    }


    // Draw the (tiny) pulsar in the middle
    glPushMatrix();
    glRotated( 360.0*t/psr.P, 0.0, 0.0, 0.1 );
    glColor3f( 0.65, 0.65, 0.65 );
    glutSolidSphere(psr.r, 100, 100);

    // Draw particles
    glColor3f(0.0, 0.65, 0.0);
    int n;
    glPointSize( 2.0 );
    glBegin( GL_POINTS );
    for (n = 0; n < npoints; n++)
        glVertex3dv( pns[n].source.x );
    glEnd();
    glPopMatrix();

    // Draw the light cylinder
    glPushMatrix();
    glColor3f( 0.65, 0.65, 0.65 );
    double z;
    int zi;
    for (zi = -10; zi <= 10; zi++)
    {
        z = 0.5*zi*psr.rL;
        draw_3D_circle( psr.rL, 0.0, 0.0, z );
    }
    int l;
    int nl = 12;
    double minz = -5.0*psr.rL;
    double maxz =  5.0*psr.rL;
    glBegin( GL_LINES );
    for( l = 0; l < nl; l++)
    {
        glVertex3d( psr.rL*cos(l*2.0*PI/nl),
                    psr.rL*sin(l*2.0*PI/nl),
                    minz );
        glVertex3d( psr.rL*cos(l*2.0*PI/nl),
                    psr.rL*sin(l*2.0*PI/nl),
                    maxz );
    }
    glEnd();
    glPopMatrix();
}


void display_gamma( int view_num )
{
    apply_2D_camera( view_num );

    // Draw a border around the perimeter of the view
    view *vw = &views[view_num];
    if (vw->border)
    {
        double borderwidth = 2.0;
        glPushMatrix();
        draw_border( view_num, borderwidth );
        glPopMatrix();
    }

    glPushMatrix();

    // Translate to "graph coords"
    glTranslated( gamma_graph.lv, gamma_graph.bv, 0.0 );
    glScaled( gamma_graph.rv - gamma_graph.lv,
              gamma_graph.tv - gamma_graph.bv,
              1.0 );

    int draw_text = 0;
    if (view_num == VIEW_LEFT_MAIN)
        draw_text = 1;
    draw_graph_frame( &gamma_graph, draw_text );

    // Now go to the graph's world coordinates
    glScaled( 1.0/(gamma_graph.xmax - gamma_graph.xmin),
              1.0/(gamma_graph.ymax - gamma_graph.ymin),
              1.0 );
    glTranslated( -gamma_graph.xmin, -gamma_graph.ymin, 0.0 );

    // Draw vertical lines representing the distribution's parameters
    glBegin( GL_LINES );
    if (gd.type == NORMAL)
    {
        if (nearest_feature == GAMMA_MEAN)
            glColor3d( 1.0, 0.0, 0.0 );
        else
            glColor3d( 0.8, 0.8, 1.0 );
        if (gamma_graph.xmin <= gd.mean && gd.mean <= gamma_graph.xmax)
        {
            glVertex2d( gd.mean, 0.0 );
            glVertex2d( gd.mean, 1.0 );
        }

        if (nearest_feature == GAMMA_STD)
            glColor3d( 1.0, 0.0, 0.0 );
        else
            glColor3d( 0.8, 0.8, 1.0 );
        double h = exp(-0.5);
        if (gamma_graph.xmin <= gd.mean + gd.std &&
            gd.mean + gd.std <= gamma_graph.xmax)
        {
            glVertex2d( gd.mean + gd.std, 0.0 );
            glVertex2d( gd.mean + gd.std, h   );
        }
        if (gamma_graph.xmin <= gd.mean - gd.std &&
            gd.mean - gd.std <= gamma_graph.xmax)
        {
            glVertex2d( gd.mean - gd.std, 0.0 );
            glVertex2d( gd.mean - gd.std, h   );
        }
    }
    else if (gd.type == POWER_LAW)
    {
        if (nearest_feature == GAMMA_MEAN)
            glColor3d( 1.0, 0.0, 0.0 );
        else
            glColor3d( 0.8, 0.8, 1.0 );
        glVertex2d( gd.g_min, 0.0 );
        glVertex2d( gd.g_max, 1.0 );

        if (nearest_feature == GAMMA_STD)
            glColor3d( 1.0, 0.0, 0.0 );
        else
            glColor3d( 0.8, 0.8, 1.0 );
        glVertex2d( gd.g_max, 0.0 );
        glVertex2d( gd.g_max, 1.0 );
    }
    glEnd();

    // Plot the gamma curve
    int nsamples = 1000;
    int n;
    double x, y;
    glColor3f( 0.0, 0.0, 1.0 );
    glBegin( GL_LINE_STRIP );
    for (n = 0; n < nsamples; n++)
    {
        x = (double)n / (double)nsamples * (gamma_graph.xmax - gamma_graph.xmin) + gamma_graph.xmin;
        if (gd.type == NORMAL)
            y = exp(-(x-gd.mean)*(x-gd.mean)/(2.0*gd.std*gd.std));
        else if (gd.type == POWER_LAW)
        {
            if (x < gd.g_min || x > gd.g_max)
                y = 0.0;
            else
                y = pow( x, gd.idx );
        }
        glVertex2d( x, y );
    }
    glEnd();

    glPopMatrix();
}


void display_angles( int view_num )
{
    apply_2D_camera( view_num );

    // Draw a border around the perimeter of the view
    view *vw = &views[view_num];
    if (vw->border)
    {
        double borderwidth = 2.0;
        glPushMatrix();
        draw_border( view_num, borderwidth );
        glPopMatrix();
    }

    glLineWidth( 1.0 );

    // Draw a circle at the origin to represent the pulsar
    glColor3f( 0.0, 0.0, 0.0 );
    draw_2D_circle( 0.5, 0.0, 0.0 );

    // Draw lines for the various angles
    glColor3f( 0.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 0.0, 1.5 );
    glEnd();

    if (nearest_feature == ALPHA_LINE)
        glColor3f( 1.0, 0.0, 0.0 );
    else
        glColor3f( 0.0, 0.0, 1.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.5*sin(psr.al.rad), 1.5*cos(psr.al.rad) );
    glEnd();

    if (nearest_feature == ZETA_LINE)
        glColor3f( 1.0, 0.0, 0.0 );
    else
        glColor3f( 0.0, 0.65, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.5*sin(psr.ze.rad), 1.5*cos(psr.ze.rad) );
    glEnd();

    glPopMatrix();

    // Draw the period graph
    glPushMatrix();
    glTranslated( period_graph.lv, period_graph.bv, 0.0 );
    glScaled( period_graph.rv - period_graph.lv,
              period_graph.tv - period_graph.bv,
              1.0 );
    glColor3d( 0.3, 0.3, 0.3 );
    glLineWidth( 2.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.0, 0.0 );
    int n, N = 6;
    for (n = 0; n < N; n++)
    {
        glVertex2d( (double)n/((double)N-1), -1.0 );
        glVertex2d( (double)n/((double)N-1),  1.0 );
    }
    glEnd();
    if (view_num == VIEW_LEFT_MAIN)
    {
        print_str( "0.001",  0.0, -3.0, GLUT_BITMAP_TIMES_ROMAN_10 );
        print_str( "0.01",   0.2, -3.0, GLUT_BITMAP_TIMES_ROMAN_10 );
        print_str( "0.1",    0.4, -3.0, GLUT_BITMAP_TIMES_ROMAN_10 );
        print_str( "1",      0.6, -3.0, GLUT_BITMAP_TIMES_ROMAN_10 );
        print_str( "10",     0.8, -3.0, GLUT_BITMAP_TIMES_ROMAN_10 );
        print_str( "100",    1.0, -3.0, GLUT_BITMAP_TIMES_ROMAN_10 );
        print_str( "Period", 0.0,  3.0, GLUT_BITMAP_TIMES_ROMAN_10 );
    }
    glPointSize( 7.0 );
    if (nearest_feature == PERIOD_SLIDER)
        glColor3f( 1.0, 0.0, 0.0 );
    else
        glColor3f( 0.0, 0.65, 0.0 );
    glBegin( GL_POINTS );
    glVertex2d( log( psr.P             / period_graph.xmin ) /
                log( period_graph.xmax / period_graph.xmin ), 0.0 );
    glEnd();
    glPopMatrix();

    glPushMatrix();

}


void display_beam( int view_num )
{
    apply_2D_camera( view_num );

    glLineWidth( 1.0 );

    // Draw a border around the perimeter of the view
    view *vw = &views[view_num];
    if (vw->border)
    {
        double borderwidth = 2.0;
        glPushMatrix();
        draw_border( view_num, borderwidth );
        glPopMatrix();
    }

    // Draw concentric circles representing lines of latitude on a polar plot
    glPushMatrix();
    glColor3d( 0.65, 0.65, 0.65 );
    int n;
    double r;
    for (n = 1; n <= 9; n++)
    {
        r = (double)n*10.0;
        draw_2D_circle( r, 0.0, 0.0 );
    }
    glPopMatrix();

    // Draw the points, whose darkness is weighted by the power
    double c; // color
    point mag; // The photon point in magnetic frame coordinates
    double x, y; // The (polar) position in Cartesian coordinates
    glPushMatrix();
    glPointSize( 1.0 );
    glEnable( GL_BLEND );
    glBegin( GL_POINTS );
    for (n = 0; n < npoints; n++)
    {
        // Set the color
        c = 1.0 - pns[n].power / max_power;
        glColor3d( 1.0, c, 0.0 );

        // Convert the point to the magnetic frame
        obs_to_mag_frame( &pns[n].retarded_LoS, &psr, NULL, &mag );

        // Draw it on the screen!
        x =  mag.th.deg*mag.ph.sin;
        y = -mag.th.deg*mag.ph.cos;
        glVertex2d( x, y );
    }
    glEnd();
    glDisable( GL_BLEND );
    glPopMatrix();

    // Draw the line-of-sight line
    point LoS, LoS_mag;
    psr_angle a;
    glPushMatrix();
    glColor3d( 0.0, 0.0, 0.0 );
    glBegin( GL_LINE_LOOP );
    for (n = 0; n < 360; n++)
    {
        set_psr_angle_deg( &a, (double)n );
        set_point_sph( &LoS, 1.0, &psr.ze, &a, POINT_SET_ALL );
        obs_to_mag_frame( &LoS, &psr, NULL, &LoS_mag );
        x =  LoS_mag.th.deg*LoS_mag.ph.sin;
        y = -LoS_mag.th.deg*LoS_mag.ph.cos;
        glVertex2d( x, y );
    }
    glEnd();

    // Draw the line-of-sight dot
    psr_angle phase;
    set_psr_angle_deg( &phase, 360.0*t/psr.P );
    line_of_sight( &psr, &phase, &LoS );
    obs_to_mag_frame( &LoS, &psr, NULL, &LoS_mag );
    x =  LoS_mag.th.deg*LoS_mag.ph.sin;
    y = -LoS_mag.th.deg*LoS_mag.ph.cos;
    glColor3d( 0.0, 0.65, 0.0 );
    glPointSize( 5.0 );
    glBegin( GL_POINTS );
    glVertex2d( x, y );
    glEnd();

    glPopMatrix();
}


void display_view(int view_num)
{
    view *vw = &views[view_num];

    glViewport ((GLsizei)vw->x, (GLsizei)vw->y,
                (GLsizei)vw->w, (GLsizei)vw->h); 

    switch (vw->scene_num)
    {
        case SCENE_FOOTPTS:
            display_footpts( view_num );
            break;
        case SCENE_FIELDLINES:
            display_fieldlines( view_num );
            break;
        case SCENE_ANGLES:
            display_angles( view_num );
            break;
        case SCENE_GAMMA:
            display_gamma( view_num );
            break;
        case SCENE_BEAM:
            display_beam( view_num );
            break;
        case SCENE_STATUS:
            display_status( view_num );
            break;
    }
}


void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    int v;
    for (v = 0; v < NVIEWS; v++)
    {
        display_view(v);
    }

    glutSwapBuffers();

    if (repopulate)
    {
        init_particles();
        repopulate = 0;
        strcpy( status_str, "" );
        glutPostRedisplay();
    }
}


void reposition_3D_camera( int view_num, double dr, double dth_deg,
        double dph_deg )
{
    view *vw = &views[view_num];
    scene *scn = &scenes[vw->scene_num];
    point *cam = &scn->camera;

    // Calculate the new camera position
    double new_th_deg = cam->th.deg + dth_deg;
    double new_ph_deg = cam->ph.deg + dph_deg;

    if (new_th_deg <    0.0)  new_th_deg =    0.0;
    if (new_th_deg >= 180.0)  new_th_deg =  180.0;
    if (new_ph_deg <    0.0)  new_ph_deg += 360.0;
    if (new_ph_deg >= 360.0)  new_ph_deg -= 360.0;

    psr_angle new_th, new_ph;
    set_psr_angle_deg( &new_ph, new_ph_deg );
    set_psr_angle_deg( &new_th, new_th_deg );

    set_point_sph( cam, cam->r + dr, &new_th, &new_ph, POINT_SET_ALL );
    set_view_properties();
}


void reshape(int w, int h)
{
    W = w;
    H = h;

    reshape_views();
    set_view_properties();

    apply_3D_camera( SCENE_FIELDLINES );
    apply_3D_camera( SCENE_ANGLES );

    glutPostRedisplay();
}

void mouseclick( int button, int state, int x, int y)
{
    active_view = which_view( x, y );
    if (active_view == SCENE_NONE) return;

    view  *vw  = &views[active_view];
    scene *scn = &scenes[vw->scene_num];
    point *cam = &scn->camera;

    double xw, yw;

    mouse_old_x = x;
    mouse_old_y = y;

    double dr, r, mid;

    if (active_view >= VIEW_THUMBNAIL_1 &&
        active_view <= VIEW_THUMBNAIL_5)
    {
        if (button == MOUSE_LEFT_BUTTON &&
            state  == GLUT_DOWN)
        {
            views[VIEW_LEFT_MAIN].scene_num = views[active_view].scene_num;
            set_view_properties();
        }
    }
    else if (active_view >= VIEW_THUMBNAIL_6 &&
             active_view <= VIEW_THUMBNAIL_8)
    {
        if (button == MOUSE_LEFT_BUTTON &&
            state  == GLUT_DOWN)
        {
            views[VIEW_RIGHT_MAIN].scene_num = views[active_view].scene_num;
            set_view_properties();
        }
    }
    else
    {
        switch (button)
        {
            case MOUSE_LEFT_BUTTON:
                if (state == GLUT_DOWN)
                {
                    button_down = MOUSE_LEFT_BUTTON;
                    selected_feature = nearest_feature;
                }
                if (state == GLUT_UP)
                {
                    if (selected_feature == GAMMA_MEAN ||
                             selected_feature == GAMMA_STD  ||
                             selected_feature == GAMMA_LO   ||
                             selected_feature == GAMMA_HI)
                    {
                        update_powers();
                    }
                    selected_feature = NO_FEATURE;
                }
                break;
            case MOUSE_SCROLL_UP:
                if (scn->dim == 2)
                {
                    vw->left   *= 5.0/6.0;
                    vw->right  *= 5.0/6.0;
                    vw->bottom *= 5.0/6.0;
                    vw->top    *= 5.0/6.0;
                }
                else if (scn->dim == 3)
                {
                    dr = -(cam->r - psr.r)/4.0;
                    reposition_3D_camera( active_view, dr, 0.0, 0.0 );
                }
                else if (vw->scene_num == SCENE_GAMMA)
                {
                    screen2graph( x, y, &xw, &yw, &gamma_graph, active_view );
                    // Must be in graph region
                    if ((gamma_graph.xmin <= xw) && (xw <= gamma_graph.xmax) &&
                        (gamma_graph.ymin <= yw) && (yw <= gamma_graph.ymax))
                    {
                        mid = (gamma_graph.xmax + gamma_graph.xmin)/2.0;
                        r   = (gamma_graph.xmax - gamma_graph.xmin)/2.0;
                        gamma_graph.xmin = mid - r*5.0/6.0;
                        gamma_graph.xmax = mid + r*5.0/6.0;
                    }
                }
                break;
            case MOUSE_SCROLL_DOWN:
                if (scn->dim == 2)
                {
                    vw->left   *= 6.0/5.0;
                    vw->right  *= 6.0/5.0;
                    vw->bottom *= 6.0/5.0;
                    vw->top    *= 6.0/5.0;
                }
                else if (scn->dim == 3)
                {
                    dr = (cam->r - psr.r)/3.0;
                    reposition_3D_camera( active_view, dr, 0.0, 0.0 );
                }
                else if (vw->scene_num == SCENE_GAMMA)
                {
                    screen2graph( x, y, &xw, &yw, &gamma_graph, active_view );
                    // Must be in graph region
                    if ((gamma_graph.xmin <= xw) && (xw <= gamma_graph.xmax) &&
                        (gamma_graph.ymin <= yw) && (yw <= gamma_graph.ymax))
                    {
                        mid = (gamma_graph.xmax + gamma_graph.xmin)/2.0;
                        r   = (gamma_graph.xmax - gamma_graph.xmin)/2.0;
                        gamma_graph.xmin = mid - r*6.0/5.0;
                        gamma_graph.xmax = mid + r*6.0/5.0;
                        if (gamma_graph.xmin < 0.0)  gamma_graph.xmin = 0.0;
                    }
                }
                break;
        }
    }

    glutPostRedisplay();
}


void mousemove( int x, int y )
/* This happens when the mouse is moved while a button is being pressed */
{
    if (active_view == SCENE_NONE) return;

    if (active_view == VIEW_THUMBNAIL_1 ||
        active_view == VIEW_THUMBNAIL_2 ||
        active_view == VIEW_THUMBNAIL_3 ||
        active_view == VIEW_THUMBNAIL_4 ||
        active_view == VIEW_THUMBNAIL_5 ||
        active_view == VIEW_THUMBNAIL_6 ||
        active_view == VIEW_THUMBNAIL_7 ||
        active_view == VIEW_THUMBNAIL_8)
    {
        return;
    }

    double xw, yw, xw_old, yw_old;
    double shift;
    psr_angle th, ph;
    psr_angle S, s;

    double angle_x, angle_y;

    screen2world( x, y, &xw, &yw, &th, &ph, active_view );

    switch (views[active_view].scene_num)
    {
        case SCENE_FIELDLINES:
            angle_x = -(x-mouse_old_x)*PI/180.0*7.6;
            angle_y = -(y-mouse_old_y)*PI/180.0*7.6;

            reposition_3D_camera( active_view, 0.0, angle_y, angle_x );

            mouse_old_x = x;
            mouse_old_y = y;
            break;
        case SCENE_ANGLES:
            if (ph.deg <   0.0)  set_psr_angle_deg( &ph,  0.0 );
            if (ph.deg >= 90.0)  set_psr_angle_deg( &ph, 90.0 );

            // See if the period control is nearest
            if (selected_feature == PERIOD_SLIDER)
            {
                screen2graph( x, y, &xw, &yw, &period_graph, active_view );
                set_pulsar_period( &psr, xw );
                adjust_tstep();
                sprintf( status_str, "period = %.3f sec", psr.P );
            }
            if (selected_feature == ALPHA_LINE)
            {
                copy_psr_angle( &ph, &psr.al );
                sprintf( status_str, "alpha = %.3f deg", psr.al.deg );
            }
            else if (selected_feature == ZETA_LINE)
            {
                copy_psr_angle( &ph, &psr.ze );
                sprintf( status_str, "zeta = %.3f deg", psr.ze.deg );
            }
            break;
        case SCENE_FOOTPTS:
            screen2csl( x, y, &S, &s, NULL, active_view );
            if (selected_feature == CSL_CIRCLE)
            {
                copy_psr_angle( &S, &psr.csl.S );
                sprintf( status_str, "Carousel size = %.3f deg", psr.csl.S.deg );
            }
            else if (selected_feature == SPARK_CIRCLE)
            {
                copy_psr_angle( &s, &psr.csl.s );
                sprintf( status_str, "Spark size = %.3f deg", psr.csl.s.deg );
            }
            else if (selected_feature == P4_ARROW)
            {
                if (ph.deg > 0.0)
                    psr.csl.P4 = ( 180.0 - ph.deg) / P4_scale;
                else
                    psr.csl.P4 = (-180.0 - ph.deg) / P4_scale;
                sprintf( status_str, "P4 = %.2f sec", fabs(psr.csl.P4) );
            }
            break;
        case SCENE_GAMMA:
            screen2graph( x, y, &xw, &yw, &gamma_graph, active_view );
            // Shift mean
            if (selected_feature == GAMMA_MEAN)
            {
                sprintf( status_str, "Mean gamma = %.2f", gd.mean );
                gd.mean = xw;
            }
            // Shift std
            else if (selected_feature == GAMMA_STD)
            {
                sprintf( status_str, "Std gamma = %.2f", gd.std );
                gd.std = fabs( xw - gd.mean );
            }
            // For panning, mouse must be in graph region
            else if ((gamma_graph.xmin <= xw) && (xw <= gamma_graph.xmax) &&
                     (gamma_graph.ymin <= yw) && (yw <= gamma_graph.ymax))
            {
                // Get old mouse pos (in graph-world coords)
                screen2graph( mouse_old_x, mouse_old_y, &xw_old, &yw_old,
                        &gamma_graph, active_view );
                shift = xw - xw_old;
                if (gamma_graph.xmin - shift >= 0.0)
                {
                    gamma_graph.xmin -= shift;
                    gamma_graph.xmax -= shift;
                }
                mouse_old_x = x;
                mouse_old_y = y;
            }
            break;
    }
    glutPostRedisplay();
}


void mousepassivemove( int x, int y )
/* This happens when the mouse is moved while NO button is being pressed */
{
    int view_num = which_view( x, y );
    if (view_num == VIEW_NONE)
    {
        nearest_feature = NO_FEATURE;
        glutPostRedisplay();
        return;
    }

    // Some possibly needed variables
    double xw, yw, w;
    double dal, dze, dS, ds, dP4;
    double P4x, P4y, P4r;
    psr_angle th, ph;
    psr_angle S, s;
    int n;

    // Behaviour is different depending on whether the mouse is over a
    // thumbnail, or over one of the larger panes

    if (view_num >= VIEW_THUMBNAIL_1 &&
        view_num <= VIEW_THUMBNAIL_8)
    {
        highlight = view_num;
    }
    else /* mouse if not over a thumbnail */
    {
        int scene_num = views[view_num].scene_num;
        screen2world( x, y, &xw, &yw, &th, &ph, view_num );
        switch (scene_num)
        {
            case SCENE_ANGLES:
                // See if the period control is nearest
                screen2graph( x, y, &xw, &yw, &period_graph, view_num );
                if (psr.P*4.0/5.0 <= xw && xw <= psr.P*5.0/4.0 &&
                    period_graph.ymin <= yw && yw <= period_graph.ymax)
                {
                    nearest_feature = PERIOD_SLIDER;
                    sprintf( status_str, "period = %.3f sec", psr.P );
                }
                else // check for the other possibilties, alpha and zeta
                {
                    dal = fabs( ph.rad - psr.al.rad );
                    dze = fabs( ph.rad - psr.ze.rad );
                    nearest_feature = (dal <= dze ?  ALPHA_LINE : ZETA_LINE);
                    if (nearest_feature == ALPHA_LINE && dal < 5.0*PI/180.0)
                    {
                        sprintf( status_str, "alpha = %.3f deg", psr.al.deg );
                    }
                    else if (nearest_feature == ZETA_LINE && dze < 5.0*PI/180.0)
                    {
                        sprintf( status_str, "zeta = %.3f deg", psr.ze.deg );
                    }
                    else
                    {
                        nearest_feature = NO_FEATURE;
                        strcpy( status_str, "" );
                    }
                }
                glutPostRedisplay();
                break;
            case SCENE_FOOTPTS:
                screen2csl( x, y, &S, &s, &n, view_num );

                dS = fabs( psr.csl.S.deg - S.deg );
                ds = fabs( psr.csl.s.deg - s.deg );

                P4r = (psr.csl.S.deg + psr.csl.s.deg)*1.3;
                P4x = P4r *  sin(psr.csl.P4*P4_scale*DEG2RAD);
                P4y = P4r * -cos(psr.csl.P4*P4_scale*DEG2RAD);
                dP4 = hypot( P4x - xw, P4y - yw );

                if ((dS <= ds) && (dS <= dP4) && (dS/psr.csl.S.deg < 0.05))
                {
                    nearest_feature = CSL_CIRCLE;
                    sprintf( status_str, "Carousel size = %.3f deg", psr.csl.S.deg );
                }
                else if ((ds <= dP4) && (ds/psr.csl.S.deg < 0.05))
                {
                    nearest_feature = SPARK_CIRCLE;
                    sprintf( status_str, "Spark size = %.3f deg", psr.csl.s.deg );
                }
                else if (dP4/psr.csl.S.deg < 0.05)
                {
                    nearest_feature = P4_ARROW;
                    sprintf( status_str, "P4 = %.2f sec", fabs(psr.csl.P4) );
                }
                else
                {
                    nearest_feature = NO_FEATURE;
                    strcpy( status_str, "" );
                }

                glutPostRedisplay();
                break;
            case SCENE_GAMMA:
                screen2graph( x, y, &xw, &yw, &gamma_graph, view_num );
                // Must be in graph region
                if ((gamma_graph.xmin <= xw) && (xw <= gamma_graph.xmax) &&
                    (gamma_graph.ymin <= yw) && (yw <= gamma_graph.ymax))
                {
                    w = (gamma_graph.xmax - gamma_graph.xmin);
                    if (gd.type == NORMAL)
                    {
                        if (fabs(xw - gd.mean) / w <= 0.02)
                        {
                            nearest_feature = GAMMA_MEAN;
                            sprintf( status_str, "Mean gamma = %.2f", gd.mean );
                        }
                        else if (fabs(xw - (gd.mean + gd.std)) / w <= 0.02 ||
                                 fabs(xw - (gd.mean - gd.std)) / w <= 0.02)
                        {
                            nearest_feature = GAMMA_STD;
                            sprintf( status_str, "Std gamma = %.2f", gd.std );
                        }
                        else
                        {
                            nearest_feature = NO_FEATURE;
                            strcpy( status_str, "" );
                        }
                    }
                }
                else
                    nearest_feature = NO_FEATURE;

                glutPostRedisplay();
                break;
            default:
                strcpy( status_str, "" );
        }
    }
}


void play()
{
    advance_particles_once();
    glutPostRedisplay();
}


void keyboard( unsigned char key, int x, int y )
{
    int view_num = which_view( x, y );
    if (view_num == VIEW_NONE)
    {
        nearest_feature = NO_FEATURE;
        glutPostRedisplay();
        return;
    }

    view *vw = &views[view_num];

    switch (key)
    {
        case 'i':
            strcpy( status_str, "Initialising particle population..." );
            playforward = 0;
            repopulate = 1;
            glutPostRedisplay();
            break;
        case 'q':
            glutDestroyWindow( glutGetWindow() );
            break;
        case '+': // Change number of sparks
            psr.csl.n++;
            glutPostRedisplay();
            break;
        case '-': // Change number of sparks
            if (psr.csl.n >= 1)
                psr.csl.n--;
            glutPostRedisplay();
            break;
        case ' ': // Push "play/pause"
            playforward = !playforward;
            if (!playforward)
                glutIdleFunc( NULL );
            else
                glutIdleFunc( play );
            break;
        case 'c': // Clear profile or pulsestack
            if (vw->scene_num == SCENE_PROFILE)
                clear_profile();
            break;
        case 'h': // Display help
            printf( "Keyboard commands (mouse cursor must be within program "
                    "window):\n" );
            printf( "  c        Clear the profile\n" );
            printf( "  h        Display this help\n" );
            printf( "  i        Initialise particle population\n" );
            printf( "  q        Quit\n" );
            printf( "  +/-      Change number of sparks\n" );
            printf( "  [space]  Play/pause time\n" );
            break;
    }
}

int main(int argc, char** argv)
{
    srand( time( NULL ) );
    W = 800;
    H = 500;
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize (W, H); 
    glutInitWindowPosition ((glutGet(GLUT_SCREEN_WIDTH )-W)/2,
                            (glutGet(GLUT_SCREEN_HEIGHT)-H)/2);
    glutCreateWindow (argv[0]);
    init ();
    glutDisplayFunc( display ); 
    glutReshapeFunc( reshape ); 
    glutMouseFunc( mouseclick );
    glutMotionFunc( mousemove );
    glutPassiveMotionFunc( mousepassivemove );
    glutKeyboardFunc( keyboard );
    glutMainLoop();

    /* Free memory */
/*
    int l;
    for (l = 0; l < MAX_NLINES; l++)
        free( line_pts[l] );

    return 0;
*/
}
