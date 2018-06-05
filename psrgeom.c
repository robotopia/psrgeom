#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "psrgeom.h"

#define MAX_NLINES  1000
#define MAX_NPOINTS 1000
static int nlines;
static int npoints[MAX_NLINES];
static point foot_pts[MAX_NLINES];
static point foot_pts_mag[MAX_NLINES];
static pulsar psr;
static point *line_pts[MAX_NLINES];
static double step_size;
static gamma_distr gd;

// Mouse-related variables
static int mouse_old_x;
static int mouse_old_y;
static int button_down;

// Determine which feature is nearer
static int nearest_feature;
static int selected_feature;

static int highlight;

static double P4_scale; // Degrees (of circular P4 arrow) per second

static double t;

enum
{
    NO_FEATURE,
    ALPHA_LINE,
    ZETA_LINE,
    CSL_CIRCLE,
    SPARK_CIRCLE,
    P4_ARROW
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
} view;

typedef struct graph_t
{
    double xmin, xmax;
    double ymin, ymax;
    int logx, logy;
    double xtics, ytics;
    double xmtics, ymtics;
    char xlabel[256], ylabel[256];
} graph;

#define NVIEWS   8
#define NSCENES  7
static view  views[NVIEWS];
static scene scenes[NSCENES];
static graph gamma_graph;

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
    VIEW_LEFT_MAIN    = 5,
    VIEW_RIGHT_TOP    = 6,
    VIEW_RIGHT_BOTTOM = 7,
    VIEW_STATUS       = 8
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

void calculate_fieldlines()
{
    int l;
    for (l = 0; l < nlines; l++)
    {
        mag_to_obs_frame( &foot_pts_mag[l], &psr, NULL, &foot_pts[l] );
        copy_point( &foot_pts[l], &line_pts[l][0] );
        npoints[l] = 1;
        do
        {
            Bstep(  &(line_pts[l][npoints[l]-1]),
                    &psr, step_size, DIR_OUTWARD,
                    &(line_pts[l][npoints[l]]) );
            set_point_xyz( &(line_pts[l][npoints[l]]),
                    line_pts[l][npoints[l]].x[0],
                    line_pts[l][npoints[l]].x[1],
                    line_pts[l][npoints[l]].x[2],
                    POINT_SET_ALL );
            npoints[l]++;
        }
        while ( (npoints[l] < MAX_NPOINTS) &&
                (line_pts[l][npoints[l]-1].r     > psr.r) &&
                (line_pts[l][npoints[l]-1].rhosq < psr.rL2) );
    }
}

void regenerate_footpoints()
{
    int i;
    for (i = 0; i < nlines; i++)
        random_spark_footpoint( NULL, &foot_pts_mag[i], &psr, 0.0 );
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
    int status_bar_height = 30;

    // Thumbnails
    int ntn = 5;          // Number of thumnails
    int tn_s = W/(2*ntn); // Size of thumbnails

    views[VIEW_THUMBNAIL_1].x = 0;
    views[VIEW_THUMBNAIL_1].y = H-tn_s;
    views[VIEW_THUMBNAIL_1].h = tn_s;
    views[VIEW_THUMBNAIL_1].w = tn_s;

    views[VIEW_THUMBNAIL_2].x = tn_s;
    views[VIEW_THUMBNAIL_2].y = H-tn_s;
    views[VIEW_THUMBNAIL_2].h = tn_s;
    views[VIEW_THUMBNAIL_2].w = tn_s;

    views[VIEW_THUMBNAIL_3].x = 2*tn_s;
    views[VIEW_THUMBNAIL_3].y = H-tn_s;
    views[VIEW_THUMBNAIL_3].h = tn_s;
    views[VIEW_THUMBNAIL_3].w = tn_s;

    views[VIEW_THUMBNAIL_4].x = 3*tn_s;
    views[VIEW_THUMBNAIL_4].y = H-tn_s;
    views[VIEW_THUMBNAIL_4].h = tn_s;
    views[VIEW_THUMBNAIL_4].w = tn_s;

    views[VIEW_THUMBNAIL_5].x = 4*tn_s;
    views[VIEW_THUMBNAIL_5].y = H-tn_s;
    views[VIEW_THUMBNAIL_5].h = tn_s;
    views[VIEW_THUMBNAIL_5].w = tn_s;

    views[VIEW_LEFT_MAIN].x = 0;
    views[VIEW_LEFT_MAIN].y = status_bar_height;
    views[VIEW_LEFT_MAIN].h = H-tn_s - status_bar_height;
    views[VIEW_LEFT_MAIN].w = W/2;

    views[VIEW_RIGHT_TOP].x = W/2;
    views[VIEW_RIGHT_TOP].y = (H + status_bar_height)/2;
    views[VIEW_RIGHT_TOP].h = (H - status_bar_height)/2;
    views[VIEW_RIGHT_TOP].w = W/2;

    views[VIEW_RIGHT_BOTTOM].x = W/2;
    views[VIEW_RIGHT_BOTTOM].y = status_bar_height;
    views[VIEW_RIGHT_BOTTOM].h = (H - status_bar_height)/2;
    views[VIEW_RIGHT_BOTTOM].w = W/2;

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
    views[VIEW_THUMBNAIL_1 ].scene_num = SCENE_FOOTPTS;
    views[VIEW_THUMBNAIL_2 ].scene_num = SCENE_FIELDLINES;
    views[VIEW_THUMBNAIL_3 ].scene_num = SCENE_ANGLES;
    views[VIEW_THUMBNAIL_4 ].scene_num = SCENE_GAMMA;
    views[VIEW_THUMBNAIL_5 ].scene_num = SCENE_TELESCOPE;
    views[VIEW_LEFT_MAIN   ].scene_num = SCENE_FIELDLINES;
    views[VIEW_RIGHT_TOP   ].scene_num = SCENE_BEAM;
    views[VIEW_RIGHT_BOTTOM].scene_num = SCENE_PROFILE;
    views[VIEW_STATUS      ].scene_num = SCENE_STATUS;

    reshape_views();
}

void set_scene_properties()
{
    scene *scn;

    scn = &scenes[SCENE_FOOTPTS];
    scn->dim = 2;

    scn = &scenes[SCENE_FIELDLINES];
    scn->dim = 3;
    scn->near_clip = 0.0;
    scn->far_clip  = 2.0*psr.rL;

    scn = &scenes[SCENE_ANGLES];
    scn->dim = 2;
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
                vw->left = -0.1;
                vw->right = 1.5;
                vw->bottom = 0.0;
                vw->top = (vw->right - vw->left)*(double)vw->h/(double)vw->w;
                break;
            case SCENE_GAMMA:
                vw->left = 0.0;
                vw->right = 1.0;
                vw->bottom = 0.0;
                vw->top = 1.0;
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
    n_tmp = (int)floor((double)N * ph0.deg / 360.0 + 0.5);
    while (n_tmp <  0)  n_tmp += N;
    while (n_tmp >= N)  n_tmp -= N;

    if (S)
        copy_psr_angle( &th, S );
    if (n)
        *n = n_tmp;
    if (s)
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
    set_pulsar( &psr, NULL, NULL, 1.0, 1.0, &al, &ze );

    // Set up default carousel
    int nsparks = 7;
    psr_angle s, S;
    set_psr_angle_deg( &S, 1.0 );
    set_psr_angle_deg( &s, 0.1 );
    double P4_sec = -10.0;
    set_pulsar_carousel( &psr, nsparks, &s, &S, GAUSSIAN, P4_sec );
    P4_scale = 2.0;

    // Set up the default gamma distribution
    gd.type  = NORMAL;
    gd.mean  = 500.0;
    gd.std   = 10.0;
    gd.idx   = -6.2;
    gd.g_min = 450.0;
    gd.g_max = 550.0;

    // Set up the default gamma distribution graph
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

    // Set the rest of the scene properties based on the cameras
    init_views();
    init_scenes();
    set_scene_properties();
    set_view_properties();

    // Generate some footpoints
    nlines = 100;
    regenerate_footpoints();

    // Allocate memory for the line points arrays,
    // and calculate line points
    int l;
    for(l = 0; l < nlines; l++)
    {
        line_pts[l] = (point *)malloc( MAX_NPOINTS * sizeof(point) );
        npoints[l] = 0;
    }
    step_size = MAX_BSTEP*psr.rL;
    calculate_fieldlines();

    // Set time to zero
    t = 0;

    // Unset nearest_feature
    nearest_feature = NO_FEATURE;
    selected_feature = NO_FEATURE;
    highlight = VIEW_NONE;

    // Unset active_view
    active_view = SCENE_NONE;
}


void draw_graph_frame( graph *gr )
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

    // Calculate the tics
    if (gr->xtics > 0.0)
    {
        double xtic_min =  ceil( gr->xmin / gr->xtics ) * gr->xtics;
        double xtic_max = floor( gr->xmax / gr->xtics ) * gr->xtics;
        double xtic;
        for (xtic = xtic_min; xtic <= xtic_max; xtic += gr->xtics)
        {
            glVertex2d( xtic, gr->ymin );
            glVertex2d( xtic, 0.05*(gr->ymax - gr->ymin) + gr->ymin );
        }
    }
    if (gr->ytics > 0.0)
    {
        double ytic_min =  ceil( gr->ymin / gr->ytics ) * gr->ytics;
        double ytic_max = floor( gr->ymax / gr->ytics ) * gr->ytics;
        double ytic;
        for (ytic = ytic_min; ytic <= ytic_max; ytic += gr->ytics)
        {
            glVertex2d( 0.0,  ytic );
            glVertex2d( 0.05, ytic );
        }
    }

    glEnd();

    glPopMatrix();
}


void display_footpts( int view_num )
{
    apply_2D_camera( view_num );

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

    // Draw footpoints
    glPushMatrix();
    glColor3f( 1.0, 0.0, 0.0 );
    glRotatef( 360.0*t/psr.csl.P4 - 90.0, 0.0, 0.0, 1.0 );
    glPointSize( 3.0 );
    glBegin( GL_POINTS );
    for( l = 0; l < nlines; l++)
    {
        glVertex2d( foot_pts_mag[l].th.deg*foot_pts_mag[l].ph.cos,
                    foot_pts_mag[l].th.deg*foot_pts_mag[l].ph.sin );
    }
    glEnd();
    glPopMatrix();

    // Draw circle of selected feature, if any
    glColor3f( 0.0, 0.0, 1.0 );
    if (nearest_feature == CSL_CIRCLE)
    {
        glPushMatrix();
        draw_2D_circle( psr.csl.S.deg, 0.0, 0.0 );
        glPopMatrix();
    }
    else if (nearest_feature == SPARK_CIRCLE)
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
    double R = (psr.csl.S.deg + psr.csl.s.deg)*1.3;
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
}


void display_fieldlines( int view_num )
{
    apply_3D_camera( view_num );
    glPushMatrix();

    glColor3f( 0.65, 0.65, 0.65 );
    glutSolidSphere(psr.r, 100, 100);

    // Draw field lines
    glColor3f(1.0, 0.0, 0.0);
    //glRotatef( psr.al.deg, 0.0, 1.0, 0.0 );
    int l, p;
    for (l = 0; l < nlines; l++)
    {
        glBegin( GL_LINE_STRIP );
        for (p = 0; p < npoints[l]; p++)
        {
            glVertex3d( line_pts[l][p].x[0],
                    line_pts[l][p].x[1],
                    line_pts[l][p].x[2] );
        }
        glEnd();
    }

    // Draw the light cylinder
    glColor3f( 0.65, 0.65, 0.65 );
    double z;
    int zi;
    for (zi = -10; zi <= 10; zi++)
    {
        z = 0.5*zi*psr.rL;
        draw_3D_circle( psr.rL, 0.0, 0.0, z );
    }
    int nl = 12;
    glBegin( GL_LINES );
    double minz = -5.0*psr.rL;
    double maxz =  5.0*psr.rL;
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

    glPushMatrix();

    // Translate to "graph coords"
    glTranslated( 0.1, 0.5, 0.0 );
    glScaled( 0.85, 0.20, 1.0 );

    draw_graph_frame( &gamma_graph );

    // Now go to the graph's world coordinates
    glScaled( 1.0/(gamma_graph.xmax - gamma_graph.xmin), 1.0/(gamma_graph.ymax - gamma_graph.ymin), 1.0 );
    glTranslated( -gamma_graph.xmin, -gamma_graph.ymin, 0.0 );

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
    glPushMatrix();

    // Draw a unit circle at the origin
    glColor3f( 0.0, 0.0, 0.0 );
    draw_2D_circle( 1.0, 0.0, 0.0 );

    // Draw lines for the various angles
    glColor3f( 0.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 0.0, 1.5 );
    glEnd();

    if (nearest_feature == ALPHA_LINE)
        glColor3f( 1.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.5*sin(psr.al.rad), 1.5*cos(psr.al.rad) );
    glEnd();
    glColor3f( 0.0, 0.0, 0.0 );

    if (nearest_feature == ZETA_LINE)
        glColor3f( 1.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.5*sin(psr.ze.rad), 1.5*cos(psr.ze.rad) );
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

    mouse_old_x = x;
    mouse_old_y = y;

    double dr;

    if (active_view == VIEW_THUMBNAIL_1 ||
        active_view == VIEW_THUMBNAIL_2 ||
        active_view == VIEW_THUMBNAIL_3 ||
        active_view == VIEW_THUMBNAIL_4)
    {
        if (button == MOUSE_LEFT_BUTTON &&
            state  == GLUT_DOWN)
        {
            views[VIEW_LEFT_MAIN].scene_num = views[active_view].scene_num;
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
                    if (selected_feature == ALPHA_LINE)
                        calculate_fieldlines();
                    if (selected_feature == CSL_CIRCLE ||
                        selected_feature == SPARK_CIRCLE)
                    {
                        regenerate_footpoints();
                        calculate_fieldlines();
                        set_view_properties();
                    }
                    selected_feature = NO_FEATURE;
                }
                break;
            case MOUSE_SCROLL_UP:
                dr = -(cam->r - psr.r)/4.0;
                reposition_3D_camera( active_view, dr, 0.0, 0.0 );
                break;
            case MOUSE_SCROLL_DOWN:
                dr = (cam->r - psr.r)/3.0;
                reposition_3D_camera( active_view, dr, 0.0, 0.0 );
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
        active_view == VIEW_THUMBNAIL_4)
    {
        return;
    }

    double xw, yw;
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

            if (selected_feature == ALPHA_LINE)
                copy_psr_angle( &ph, &psr.al );
            else if (selected_feature == ZETA_LINE)
                copy_psr_angle( &ph, &psr.ze );
            break;
        case SCENE_FOOTPTS:
            screen2csl( x, y, &S, &s, NULL, active_view );
            if (selected_feature == CSL_CIRCLE)
                copy_psr_angle( &S, &psr.csl.S );
            else if (selected_feature == SPARK_CIRCLE)
                copy_psr_angle( &s, &psr.csl.s );
            else if (selected_feature == P4_ARROW)
            {
                if (ph.deg > 0.0)
                    psr.csl.P4 = ( 180.0 - ph.deg) / P4_scale;
                else
                    psr.csl.P4 = (-180.0 - ph.deg) / P4_scale;
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
    double xw, yw;
    double dal, dze, dS, ds, dP4;
    double P4x, P4y, P4r;
    psr_angle th, ph;
    psr_angle S, s;
    int n;

    // Behaviour is different depending on whether the mouse is over a
    // thumbnail, or over one of the larger panes

    if (view_num == VIEW_THUMBNAIL_1 ||
        view_num == VIEW_THUMBNAIL_2 ||
        view_num == VIEW_THUMBNAIL_3 ||
        view_num == VIEW_THUMBNAIL_4 )
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
                dal = fabs( ph.rad - psr.al.rad );
                dze = fabs( ph.rad - psr.ze.rad );
                nearest_feature = (dal <= dze ?  ALPHA_LINE : ZETA_LINE);
                if (nearest_feature == ALPHA_LINE && dal > 5.0*PI/180.0)
                    nearest_feature = NO_FEATURE;
                else if (nearest_feature == ZETA_LINE && dze > 5.0*PI/180.0)
                    nearest_feature = NO_FEATURE;
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
                    nearest_feature = CSL_CIRCLE;
                else if ((ds <= dP4) && (ds/psr.csl.S.deg < 0.05))
                    nearest_feature = SPARK_CIRCLE;
                else if (dP4/psr.csl.S.deg < 0.05)
                    nearest_feature = P4_ARROW;
                else
                    nearest_feature = NO_FEATURE;

                glutPostRedisplay();
        }
    }
}


void keyboard( unsigned char key, int x, int y )
{
    if (x || y) {} // Dummy line to avoid compiler warnings

    switch (key)
    {
        case 'a':
            regenerate_footpoints();
            calculate_fieldlines();
            glutPostRedisplay();
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
    int l;
    for (l = 0; l < MAX_NLINES; l++)
        free( line_pts[l] );

    return 0;
}
